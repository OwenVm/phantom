!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module inject
!
! Handles initial setup of stellar atmosphere with pulsating boundary layers.
!
! Variable-resolution atmosphere: each shell carries a number of particles
! proportional to the local density, so that all gas particles share a
! single common mass and all boundary particles share a (generally
! different) common mass.  This naturally accommodates steep density
! profiles (rho_power >> 2) without creating an excessive number of
! boundary particles.
!
! Gas and boundary shells are built on independent radial grids: the
! spacing of each grid is derived iteratively from the local per-shell
! particle count via get_fibonacci_spacing, so that the inter-particle
! distance is isotropic everywhere.  The ratio of boundary to gas
! particles per shell is controlled by boundary_fraction.
!
! :References: None
!
! :Owner: Owen Vermeulen
!
! :Runtime parameters:
!   - iboundary_spheres   : *number of boundary spheres (integer)*
!   - n_profile_points    : *number of points in stellar profile calculation (integer)*
!   - n_particles         : *total number of gas particles (integer, not used if n_shells > 0)*
!   - n_shells            : *total number of shells (integer, if <0 determined automatically from n_particles)*
!   - boundary_fraction   : *ratio N_boundary_per_shell / N_gas_per_shell (controls boundary resolution)*
!   - rho_power           : *density profile exponent: rho ~ r^(-rho_power)*
!   - r_min_on_rstar      : *inner radius as fraction of R_star*
!   - r_max_on_rstar      : *outer radius as fraction of R_star*
!   - pulsation_period    : *pulsation period (days)*
!   - pulsation_amplitude : *fractional pulsation amplitude*
!   - piston_velocity     : *piston velocity amplitude (km/s)*
!   - atmos_mass_fraction : *atmospheric mass as fraction of total stellar mass*
!   - surface_pressure    : *surface pressure (cgs)*
!   - iwind               : *wind type: 1=prescribed, 2=period from mass-radius relation*
!   - pulsation_timestep  : *pulsation timestep as fraction of pulsation period*
!   - phi0                : *initial phase offset (radians)*
!   - wss                 : *fraction of tangential and radial distance between particles*
!   - reinject_enabled    : *enable dynamic reinjection (logical)*
!   - reinject_period_days: *period between reinjections in days*
!   - mass_loss_start     : *start time for mass-loss calculation in years*
!   - mass_loss_end       : *end time for mass-loss calculation in years*
!   - check_radius_au     : *radius within which to count mass (AU)*
!   - meas_int_days       : *interval for mass measurements in days*
!
! :Dependencies: dim, eos, icosahedron, infile_utils, injectutils, io,
!   part, partinject, physcon, units, set_star
!
 use io, only:fatal
 implicit none
 character(len=*), parameter, public :: inject_type = 'atmosphere'

 public :: init_inject, inject_particles, write_options_inject, read_options_inject, &
           set_default_options_inject, update_injected_par
 private

 !------------------------------------------------------------------
 ! Runtime parameters (read from input file)
 !------------------------------------------------------------------
 integer :: iboundary_spheres    = 5
 integer :: n_profile_points     = 10000
 integer :: n_particles          = 500000
 integer :: n_shells             = 25
 ! Ratio N_boundary_per_shell / N_gas_per_shell at the interface (r_min).
 ! Each boundary shell gets particles proportional to local density scaled
 ! by this fraction relative to gas shells at the same density.
 real    :: boundary_fraction       = 0.1
 real    :: rho_power_in            = 2.0
 real    :: r_min_on_rstar          = 0.9
 real    :: r_max_on_rstar          = 1.4
 ! Inner edge of boundary region — also sets the inner edge of the
 ! hydrostatic profile integration, so no extrapolation occurs.
 real    :: r_b_min_on_rmax = 0.7
 real    :: dtpulsation          = huge(0.)
 real    :: pulsation_period_days= 300.0
 real    :: piston_velocity_km_s = 4.0
 real    :: time_puls            = -1.0
 real    :: atmos_mass_fraction  = 5e-5
 real    :: surface_pressure     = 0.001
 integer :: iwind                = 1
 real    :: pulsation_timestep   = 0.02
 real    :: phi0                 = -3.1415926536d0/2.0
 real    :: wss                  = 1.0
 logical :: var_boundary         = .false.

 ! Reinjection parameters
 logical :: reinject_enabled      = .true.
 real    :: reinject_period_days  = 10.0
 real    :: mass_loss_start       = 1.0
 real    :: mass_loss_end         = 3.0
 real    :: check_radius_au       = 3.0
 real    :: meas_int_days         = 10.0

 !------------------------------------------------------------------
 ! Global variables
 !------------------------------------------------------------------
 integer, parameter :: wind_emitting_sink = 1
 integer, parameter :: max_measurements   = 10000

 real :: omega_pulsation, deltaR_osc, pulsation_period, piston_velocity
 real :: Rstar, r_min, r_boundary_min
 real :: Mtotal, Matmos, Msink

 ! Per-shell particle counts (allocated to n_shells_total)
 integer, allocatable :: npart_per_shell(:)          ! gas shells (iboundary_spheres+1 : n_shells_total)
 integer, allocatable :: npart_per_boundary_shell(:) ! boundary shells (1 : iboundary_spheres)

 ! Independent radial grids for gas and boundary shells
 real, allocatable :: delta_r_gas(:)      ! spacing for gas shells
 real, allocatable :: delta_r_boundary(:) ! spacing for boundary shells (independent grid)
 real, allocatable :: shell_radii_gas(:)  ! central radii of gas shells
 real, allocatable :: shell_radii_bnd(:)  ! central radii of boundary shells

 ! Single particle masses (constant within each particle type)
 real :: mass_of_gas_particle      = 0.0
 real :: mass_of_boundary_particle = 0.0

 ! Combined delta_r array (boundary shells first, then gas) used by
 ! apply_pulsation and perform_reinjection which index by shell number.
 real, allocatable :: delta_r_radial(:)

 logical :: atmosphere_setup_complete = .false.
 integer :: n_shells_total   ! total gas shells
 integer :: n_shells_bnd     ! total boundary shells (== iboundary_spheres)

 ! Boundary particle tracking
 real, allocatable    :: r_boundary_equilibrium(:)
 integer, allocatable :: boundary_particle_ids(:)
 integer              :: n_boundary_particles
 integer              :: active_boundary_spheres

 logical :: reinjection_needed = .false.

 ! Continuous reinjection tracking
 real    :: time_last_reinject = 0.0
 real    :: reinject_period
 integer :: n_reinjections     = 0

 ! Mass-loss rate tracking
 real    :: mass_loss_start_time
 real    :: mass_loss_end_time
 real    :: mass_loss_check_radius
 real    :: measurement_interval
 real    :: time_next_measurement
 real    :: mass_previous_measurement
 real, allocatable :: mass_loss_rates(:)
 integer :: n_measurements              = 0
 real    :: mean_mass_loss_rate         = 0.0
 logical :: mass_loss_rate_calculated   = .false.
 logical :: measurement_active          = .false.
 integer :: particles_to_inject         = 0

 character(len=*), parameter :: label = 'inject_atmosphere'

contains

!-----------------------------------------------------------------------
!+
!  Set default options
!+
!-----------------------------------------------------------------------
subroutine set_default_options_inject(flag)
 integer, optional, intent(in) :: flag

 iboundary_spheres          = 5
 n_profile_points           = 10000
 n_particles                = 500000
 n_shells                   = 25
 boundary_fraction          = 0.1
 rho_power_in               = 2.0
 r_min_on_rstar             = 0.9
 r_max_on_rstar             = 1.4
 r_b_min_on_rmax    = 0.7
 dtpulsation          = huge(0.)
 atmos_mass_fraction  = 5e-5
 surface_pressure     = 0.001
 iwind                = 1
 pulsation_period_days= 300.0
 piston_velocity_km_s = 4.0
 time_puls            = -1.0
 pulsation_timestep   = 0.02
 phi0                 = -3.1415926536d0/2.0
 wss                  = 1.0
 var_boundary         = .false.
 reinject_enabled     = .true.
 reinject_period_days = 10.0
 mass_loss_start      = 1.0
 mass_loss_end        = 3.0
 check_radius_au      = 3.0
 meas_int_days        = 10.0

end subroutine set_default_options_inject

!-----------------------------------------------------------------------
!+
!  Compute per-shell particle counts so that every particle in a given
!  population (gas or boundary) has the same mass.
!
!  For a shell at radius r with local density rho(r), the shell mass is
!    M_shell = 4*pi * r^2 * dr * rho(r)      (thin-shell approximation)
!  We want M_shell / N_shell = m_particle = const, so
!    N_shell  proportional to  M_shell  proportional to  rho(r) * r^2 * dr
!
!  In practice we normalise so that sum(N_shell) = n_total_target, then
!  enforce N_shell >= 1.
!
!  shell_r(1:nsh)  - central radius of each shell  [code units]
!  shell_dr(1:nsh) - radial width of each shell     [code units]
!  nsh             - number of shells
!  n_total_target  - desired total particle count
!  npart_shell(:)  - output: integer particles per shell
!+
!-----------------------------------------------------------------------
subroutine calc_particles_per_shell(shell_r, shell_dr, nsh, n_total_target, npart_shell)
 use wind_pulsating, only:interp_stellar_profile
 use units,          only:unit_density

 real,    intent(in)  :: shell_r(:), shell_dr(:)
 integer, intent(in)  :: nsh, n_total_target
 integer, intent(out) :: npart_shell(nsh)

 integer :: i
 real    :: rho, P, u, T
 real    :: weight(nsh), weight_sum, scale

 ! Compute shell mass weights  ~ rho(r) * r^2 * dr
 weight_sum = 0.0
 do i = 1, nsh
    call interp_stellar_profile(shell_r(i), rho, P, u, T)
    ! rho returned in code units; we only need relative weights
    weight(i) = rho * shell_r(i)**2 * shell_dr(i)
    weight_sum = weight_sum + weight(i)
 enddo

 ! Scale so sum == n_total_target, then round
 scale = real(n_total_target) / weight_sum
 do i = 1, nsh
    npart_shell(i) = max(1, nint(weight(i) * scale))
 enddo

end subroutine calc_particles_per_shell

!-----------------------------------------------------------------------
!+
!  Initialize atmospheric setup and pulsation parameters
!+
!-----------------------------------------------------------------------
subroutine init_inject(ierr)
 use io,            only:fatal
 use physcon,       only:pi,days,au,solarm,km,years
 use icosahedron,   only:compute_matrices,compute_corners
 use eos,           only:gmw,gamma
 use units,         only:utime,umass,unit_velocity,unit_luminosity
 use part,          only:xyzmh_ptmass,massoftype,igas,iboundary,nptmass,iTeff,iReff,iLum
 use injectutils,   only:get_parts_per_sphere, get_fibonacci_spacing
 use wind_pulsating,only:setup_star,calc_stellar_profile,interp_stellar_profile
 use dust_formation,only:calc_kappa_max

 integer, intent(out) :: ierr
 real    :: Mstar_cgs, Rstar_cgs, Tstar, Lstar_cgs
 real    :: delta_r_tangential, current_radius
 integer :: shell_index, max_shells, temp_particles, particles_per_shell_ref
 integer :: expected_measurements, i
 integer, parameter  :: max_shells_tmp = 2000
 real    :: tmp_dr_gas(max_shells_tmp), tmp_r_gas(max_shells_tmp)
 real    :: tmp_dr_bnd(max_shells_tmp), tmp_r_bnd(max_shells_tmp)
 integer :: tmp_nbnd(max_shells_tmp)
 integer :: iter, n_bnd_iter, n_bnd_ref
 real    :: dr_bnd_iter, r_shell_bnd, rho_bnd, P_bnd, u_bnd, T_bnd
 real    :: gas_weight_ref, bnd_weight, gas_scale
 logical :: converged, file_exists

 ierr = 0

 if (nptmass < 1) call fatal(label,'need at least one sink particle for central star')

 Mtotal    = xyzmh_ptmass(4, wind_emitting_sink)
 Rstar     = xyzmh_ptmass(iReff, wind_emitting_sink)
 Rstar_cgs = Rstar * au
 Mstar_cgs = Mtotal * solarm
 Tstar     = xyzmh_ptmass(iTeff, wind_emitting_sink)
 Lstar_cgs = xyzmh_ptmass(iLum, wind_emitting_sink) * unit_luminosity

 call calc_kappa_max(Mstar_cgs, Lstar_cgs)

 Matmos = atmos_mass_fraction * Mtotal

 inquire(file='mass_loss_rate.dat', exist=file_exists)

 if (.not. file_exists) then
    Msink = Mtotal - Matmos
    xyzmh_ptmass(4, wind_emitting_sink) = Msink
 endif

 active_boundary_spheres = iboundary_spheres

 if (iwind == 2 .and. .not. file_exists) call calculate_period(Mtotal, Rstar, pulsation_period_days)

 pulsation_period = pulsation_period_days * (days / utime)
 omega_pulsation  = 2.0*pi / pulsation_period
 piston_velocity  = piston_velocity_km_s * (km / unit_velocity)
 deltaR_osc       = pulsation_period * piston_velocity / (2.0*pi)

 reinject_period        = reinject_period_days * (days / utime)
 mass_loss_start_time   = mass_loss_start  * (years / utime)
 mass_loss_end_time     = mass_loss_end    * (years / utime)
 mass_loss_check_radius = check_radius_au
 measurement_interval   = meas_int_days   * (days / utime)
 time_next_measurement  = mass_loss_start_time
 n_measurements         = 0

 expected_measurements = ceiling((mass_loss_end_time - mass_loss_start_time) / measurement_interval) + 1
 allocate(mass_loss_rates(expected_measurements))
 mass_loss_rates = 0.0

 if (n_shells > 0) then
    max_shells = n_shells
 else
    max_shells = 200
 endif

 call setup_star(Msink * umass, r_max_on_rstar * Rstar * au, r_b_min_on_rmax * Rstar * au, &
                 gmw, gamma, surface_pressure, Matmos * umass, rho_power_in)

 r_min          = r_min_on_rstar          * Rstar
 r_boundary_min = r_b_min_on_rmax * Rstar

 if (r_boundary_min >= r_min) &
    call fatal(label,'r_b_min_on_rmax must be less than r_min_on_rstar')

 ! Hydrostatic profile must exist before any interp calls below
 call calc_stellar_profile(n_profile_points)

 ! ================================================================
 ! STEP 1: Build gas shell grid outward from r_min
 !
 ! Converge on a reference per-shell count particles_per_shell_ref
 ! such that the shells fill [r_min, r_max] with the desired number
 ! of shells.  The spacing dr(i) = wss * r(i) * fibonacci(N_ref)
 ! uses a uniform N_ref at this stage — the actual density-weighted
 ! counts come in Step 2.  The grid only needs self-consistent
 ! isotropic spacing, so the converged N_ref is appropriate here.
 ! ================================================================
 if (n_shells < 0) then
    temp_particles          = n_particles
    particles_per_shell_ref = nint(real(temp_particles) / real(max_shells))
    converged = .false.
    do while (.not. converged)
       current_radius          = r_min
       particles_per_shell_ref = particles_per_shell_ref + 10
       shell_index             = 0
       do while (current_radius < r_max_on_rstar * Rstar)
          shell_index = shell_index + 1
          if (shell_index > max_shells_tmp) call fatal(label,'gas shell tmp array too small')
          delta_r_tangential      = current_radius * get_fibonacci_spacing(particles_per_shell_ref)
          tmp_dr_gas(shell_index) = wss * delta_r_tangential
          tmp_r_gas(shell_index)  = current_radius + 0.5*tmp_dr_gas(shell_index)
          current_radius          = current_radius + tmp_dr_gas(shell_index)
       enddo
       if (shell_index * particles_per_shell_ref >= temp_particles) converged = .true.
    enddo
 else
    temp_particles          = 100
    particles_per_shell_ref = nint(real(temp_particles) / real(max_shells))
    converged = .false.
    do while (.not. converged)
       current_radius          = r_min
       particles_per_shell_ref = particles_per_shell_ref + 10
       shell_index             = 0
       do while (current_radius < r_max_on_rstar * Rstar)
          shell_index = shell_index + 1
          if (shell_index > max_shells_tmp) call fatal(label,'gas shell tmp array too small')
          delta_r_tangential      = current_radius * get_fibonacci_spacing(particles_per_shell_ref)
          tmp_dr_gas(shell_index) = wss * delta_r_tangential
          tmp_r_gas(shell_index)  = current_radius + 0.5*tmp_dr_gas(shell_index)
          current_radius          = current_radius + tmp_dr_gas(shell_index)
       enddo
       if (shell_index >= max_shells) converged = .true.
    enddo
 endif

 n_shells_total = shell_index

 ! ================================================================
 ! STEP 2: Density-weighted per-shell particle counts for gas shells
 !
 ! N_gas(i) ∝ rho(r_i) * r_i^2 * dr_i  =>  equal mass per particle
 ! Normalised so sum == target gas count.
 ! ================================================================
 if (n_shells < 0) then
    temp_particles = n_particles
 else
    temp_particles = n_shells_total * particles_per_shell_ref
 endif

 allocate(npart_per_shell(n_shells_total))
 allocate(delta_r_gas(n_shells_total))
 allocate(shell_radii_gas(n_shells_total))
 delta_r_gas     = tmp_dr_gas(1:n_shells_total)
 shell_radii_gas = tmp_r_gas(1:n_shells_total)

 call calc_particles_per_shell(shell_radii_gas, delta_r_gas, n_shells_total, &
                                temp_particles, npart_per_shell)

 ! ================================================================
 ! STEP 3: Build boundary shell grid filling [r_boundary_min, r_min]
 !
 ! iboundary_spheres shells are placed inward from r_min down to
 ! r_boundary_min.  The spacing must be self-consistent with the
 ! local particle count via Fibonacci, AND the outermost boundary
 ! shell must honour the interface condition at r_min:
 !
 !   N_bnd_outer = boundary_fraction * N_gas_inner
 !
 ! where N_gas_inner is the innermost gas shell count.  This fixes
 ! the spacing of the outermost boundary shell, and the remaining
 ! shells scale by density inward from there.
 !
 ! We converge on a boundary reference count N_bnd_ref such that
 ! exactly iboundary_spheres shells of density-weighted spacing
 ! fill [r_boundary_min, r_min].  The loop mirrors the gas grid
 ! convergence loop exactly.
 ! ================================================================
 n_shells_bnd = iboundary_spheres

 if (n_shells_bnd < 1) then
    allocate(npart_per_boundary_shell(0))
    allocate(delta_r_boundary(0))
    allocate(shell_radii_bnd(0))
 else
    ! Interface condition: outermost boundary shell particle count.
    ! N_bnd_ref is the reference count for the boundary grid,
    ! analogous to particles_per_shell_ref for the gas grid.
    ! We initialise from the interface constraint and then converge.
    !
    ! N_bnd at interface = boundary_fraction * N_gas(innermost shell)
    ! but N_gas(innermost) is not yet available here since it is
    ! density-weighted; use npart_per_shell(1) computed in Step 2.
    !
    ! The boundary convergence loop: increase N_bnd_ref until
    ! iboundary_spheres shells of spacing
    !   dr(i) = wss * r(i) * fibonacci(N_bnd_at_shell_i)
    ! exactly fill [r_boundary_min, r_min].
    ! N_bnd at each shell scales with density relative to N_bnd_ref
    ! at r_min, consistent with boundary_fraction * gas weighting.

    ! gas_scale: gas particles per unit geometric weight at innermost gas shell
    gas_weight_ref = shell_radii_gas(1)**2 * delta_r_gas(1)
    gas_scale      = real(npart_per_shell(1)) / gas_weight_ref

    ! Starting N_bnd_ref from interface constraint
    n_bnd_ref = max(1, nint(boundary_fraction * real(npart_per_shell(1))))

    converged = .false.
    do while (.not. converged)
       current_radius = r_min
       shell_index    = 0

       do while (current_radius > r_boundary_min)
          shell_index = shell_index + 1
          if (shell_index > max_shells_tmp) call fatal(label,'boundary shell tmp array too small')

          ! N_bnd at this shell: density-weighted relative to N_bnd_ref at r_min.
          ! Use 3-iteration inner loop to self-consistently determine
          ! N and dr (they depend on each other via Fibonacci spacing).
          n_bnd_iter = n_bnd_ref
          do iter = 1, 3
             dr_bnd_iter = wss * current_radius * get_fibonacci_spacing(n_bnd_iter)
             r_shell_bnd = current_radius - 0.5 * dr_bnd_iter
             ! Clamp to profile domain to avoid extrapolation during convergence
             r_shell_bnd = max(r_shell_bnd, r_boundary_min + tiny(0.))
             call interp_stellar_profile(r_shell_bnd, rho_bnd, P_bnd, u_bnd, T_bnd)
             bnd_weight = r_shell_bnd**2 * dr_bnd_iter
             n_bnd_iter = max(1, nint(boundary_fraction * gas_scale * bnd_weight))
          enddo

          dr_bnd_iter             = wss * current_radius * get_fibonacci_spacing(n_bnd_iter)
          tmp_nbnd(shell_index)   = n_bnd_iter
          tmp_dr_bnd(shell_index) = dr_bnd_iter
          tmp_r_bnd(shell_index)  = current_radius - 0.5 * dr_bnd_iter
          current_radius          = current_radius - dr_bnd_iter
       enddo

       if (shell_index >= n_shells_bnd) then
          converged = .true.
       else
          ! Too few shells — reduce N_bnd_ref to increase spacing and
          ! pack more shells into [r_boundary_min, r_min]
          n_bnd_ref = max(1, n_bnd_ref - 1)
          if (n_bnd_ref == 1) then
             ! Cannot reduce further; warn and accept however many shells fit
             print *, 'Warning: could not fit iboundary_spheres=', n_shells_bnd, &
                      ' boundary shells in [r_boundary_min, r_min]; got ', shell_index
             converged = .true.
          endif
       endif
    enddo

    ! Use exactly iboundary_spheres outermost shells from the converged grid
    ! (extra shells beyond iboundary_spheres are discarded)
    n_shells_bnd = min(n_shells_bnd, shell_index)

    allocate(npart_per_boundary_shell(n_shells_bnd))
    allocate(delta_r_boundary(n_shells_bnd))
    allocate(shell_radii_bnd(n_shells_bnd))

    ! tmp arrays are outermost->innermost; reverse to innermost->outermost
    do i = 1, n_shells_bnd
       npart_per_boundary_shell(i) = tmp_nbnd(n_shells_bnd + 1 - i)
       delta_r_boundary(i)         = tmp_dr_bnd(n_shells_bnd + 1 - i)
       shell_radii_bnd(i)          = tmp_r_bnd(n_shells_bnd + 1 - i)
    enddo
 endif

 ! ================================================================
 ! STEP 4: Combined delta_r_radial (boundary first, then gas) so
 ! that routines indexing by shell number remain correct.
 ! ================================================================
 if (allocated(delta_r_radial)) deallocate(delta_r_radial)
 allocate(delta_r_radial(n_shells_bnd + n_shells_total))
 if (n_shells_bnd > 0) delta_r_radial(1:n_shells_bnd) = delta_r_boundary
 delta_r_radial(n_shells_bnd+1 : n_shells_bnd+n_shells_total) = delta_r_gas

 ! ================================================================
 ! STEP 5: Particle masses — constant within each type
 ! ================================================================
 mass_of_gas_particle      = Matmos / real(sum(npart_per_shell))
 if (n_shells_bnd > 0) then
    mass_of_boundary_particle = Matmos / real(sum(npart_per_boundary_shell))
 else
    mass_of_boundary_particle = mass_of_gas_particle
 endif

 massoftype(igas)      = mass_of_gas_particle
 massoftype(iboundary) = mass_of_boundary_particle

 if (file_exists) call read_mass_loss_data()

 print *, ''
 print *, ' rho_power                        :', rho_power_in
 print *, ' boundary_fraction                :', boundary_fraction
 print *, ' Boundary region  [r_bnd_min, r_min] / Rstar :', r_b_min_on_rmax, r_min_on_rstar
 print *, ' Gas region       [r_min,     r_max] / Rstar :', r_min_on_rstar, r_max_on_rstar
 print *, ' Gas shells                       :', n_shells_total
 print *, ' Boundary shells                  :', n_shells_bnd
 print *, ' Total gas particles              :', sum(npart_per_shell)
 if (n_shells_bnd > 0) then
    print *, ' Total boundary particles         :', sum(npart_per_boundary_shell)
    print *, ' Innermost boundary N_per_shell   :', npart_per_boundary_shell(1)
    print *, ' Outermost boundary N_per_shell   :', npart_per_boundary_shell(n_shells_bnd)
    print *, ' Innermost gas     N_per_shell    :', npart_per_shell(1)
    print *, ' Interface ratio (actual)         :', &
              real(npart_per_boundary_shell(n_shells_bnd)) / real(npart_per_shell(1))
 endif
 print *, ' Outermost gas N_per_shell        :', npart_per_shell(n_shells_total)
 print *, ' Gas particle mass (Msun)         :', mass_of_gas_particle
 print *, ' Boundary particle mass (Msun)    :', mass_of_boundary_particle
 print *, ' Boundary/gas mass ratio          :', mass_of_boundary_particle / mass_of_gas_particle
 print *, ''

end subroutine init_inject

!-----------------------------------------------------------------------
!+
!  Main routine: called at start to setup atmosphere, then each timestep
!+
!-----------------------------------------------------------------------
subroutine inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npart_old,npartoftype,dtinject)
 use part, only:igas,iboundary,iamtype

 real,    intent(in)    :: time,dtlast
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(inout) :: npart,npart_old
 integer, intent(inout) :: npartoftype(:)
 real,    intent(out)   :: dtinject

 dtinject = pulsation_timestep * pulsation_period

 if (npart > 0 .and. .not. atmosphere_setup_complete) then
    atmosphere_setup_complete = .true.
 endif

 if (.not. atmosphere_setup_complete) then
    print *, 'Setting up stellar atmosphere with ', n_shells_total, ' shells.'
    call setup_initial_atmosphere(xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npartoftype)
    atmosphere_setup_complete = .true.
    print *, 'Stellar atmosphere setup complete.'
    return
 endif

 if (atmosphere_setup_complete .and. .not. allocated(boundary_particle_ids)) then
    call reconstruct_boundary_info(time,xyzh,npart,xyzmh_ptmass)
    time_last_reinject = time - mod(time, reinject_period)
    call read_mass_loss_data()
 endif

 if (reinject_enabled .and. .not. mass_loss_rate_calculated) then
    call take_periodic_mass_measurements(time,xyzh,npart,xyzmh_ptmass,npartoftype)
 endif

 if (reinject_enabled .and. mass_loss_rate_calculated) then
    call check_continuous_reinject(time,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npartoftype)
    if (reinjection_needed) then
       call perform_reinjection(time,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npartoftype)
       time_last_reinject = time
       reinjection_needed = .false.
    endif
 endif

 call apply_pulsation(time,xyzh,vxyzu,npart,xyzmh_ptmass,vxyz_ptmass)

end subroutine inject_particles

!-----------------------------------------------------------------------
!+
!  Take periodic mass measurements and calculate mean mass-loss rate
!+
!-----------------------------------------------------------------------
subroutine take_periodic_mass_measurements(time,xyzh,npart,xyzmh_ptmass,npartoftype)
 use part,   only:igas,iboundary,iphase,iamtype
 use units,  only:utime,umass
 use physcon,only:solarm,years,days

 real,    intent(in) :: time
 real,    intent(in) :: xyzh(:,:),xyzmh_ptmass(:,:)
 integer, intent(in) :: npart
 integer, intent(in) :: npartoftype(:)

 real    :: current_mass_within_radius, mass_lost, rate_this_interval
 real    :: sink_mass, x0(3), dx, dy, dz, r
 real    :: sum_rates
 integer :: i

 if (.not. measurement_active .and. time >= mass_loss_start_time) then
    measurement_active = .true.
    x0        = xyzmh_ptmass(1:3, wind_emitting_sink)
    sink_mass = xyzmh_ptmass(4,   wind_emitting_sink)
    mass_previous_measurement = sink_mass
    do i = 1, npart
       dx = xyzh(1,i) - x0(1)
       dy = xyzh(2,i) - x0(2)
       dz = xyzh(3,i) - x0(3)
       r  = sqrt(dx**2 + dy**2 + dz**2)
       if (r <= mass_loss_check_radius) then
          ! Use appropriate particle mass per type
          if (iamtype(iphase(i)) == iboundary) then
             mass_previous_measurement = mass_previous_measurement + mass_of_boundary_particle
          else
             mass_previous_measurement = mass_previous_measurement + mass_of_gas_particle
          endif
       endif
    enddo
    time_next_measurement = time + measurement_interval
 endif

 if (measurement_active .and. time >= time_next_measurement .and. time < mass_loss_end_time) then
    x0        = xyzmh_ptmass(1:3, wind_emitting_sink)
    sink_mass = xyzmh_ptmass(4,   wind_emitting_sink)
    current_mass_within_radius = sink_mass
    do i = 1, npart
       dx = xyzh(1,i) - x0(1)
       dy = xyzh(2,i) - x0(2)
       dz = xyzh(3,i) - x0(3)
       r  = sqrt(dx**2 + dy**2 + dz**2)
       if (r <= mass_loss_check_radius) then
          if (iamtype(iphase(i)) == iboundary) then
             current_mass_within_radius = current_mass_within_radius + mass_of_boundary_particle
          else
             current_mass_within_radius = current_mass_within_radius + mass_of_gas_particle
          endif
       endif
    enddo
    mass_lost          = mass_previous_measurement - current_mass_within_radius
    rate_this_interval = mass_lost / measurement_interval
    n_measurements     = n_measurements + 1
    mass_loss_rates(n_measurements) = rate_this_interval
    mass_previous_measurement = current_mass_within_radius
    time_next_measurement     = time + measurement_interval
 endif

 if (measurement_active .and. .not. mass_loss_rate_calculated .and. time >= mass_loss_end_time) then
    if (n_measurements > 0) then
       sum_rates = 0.0
       do i = 1, n_measurements
          sum_rates = sum_rates + mass_loss_rates(i)
       enddo
       mean_mass_loss_rate = sum_rates / real(n_measurements)
       ! Reinjected particles are gas particles placed just outside the boundary
       particles_to_inject = nint((mean_mass_loss_rate * reinject_period) / mass_of_gas_particle)
       if (particles_to_inject < 1) particles_to_inject = 1
       mass_loss_rate_calculated = .true.
       call write_mass_loss_data()
    endif
 endif

end subroutine take_periodic_mass_measurements

!-----------------------------------------------------------------------
!+
!  Check if it is time for continuous reinjection
!+
!-----------------------------------------------------------------------
subroutine check_continuous_reinject(time,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npartoftype)
 use units,  only:utime
 use physcon,only:days

 real,    intent(in)    :: time
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)

 if ((time - time_last_reinject) < reinject_period .and. time >= mass_loss_start_time) return

 print *, ''
 print *, '-----------------------------------------'
 print *, 'Reinjection triggered at time: ', time
 print *, 'Time since last reinject: ', (time - time_last_reinject)*utime/days
 print *, '-----------------------------------------'
 print *, ''

 reinjection_needed = .true.

end subroutine check_continuous_reinject

!-----------------------------------------------------------------------
!+
!  Perform reinjection: inject new gas particles to replenish lost mass
!+
!-----------------------------------------------------------------------
subroutine perform_reinjection(time,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npartoftype)
 use part,           only:igas,iboundary,iamtype,set_particle_type
 use injectutils,    only:inject_fibonacci_sphere
 use wind_pulsating, only:interp_stellar_profile
 use physcon,        only:pi

 real,    intent(in)    :: time
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)

 integer :: old_npart, i
 real    :: r_inject, phase, r_dot, rho, u, T, P
 real    :: x0(3), v0(3)
 real    :: mass_injected

 x0 = xyzmh_ptmass(1:3, wind_emitting_sink)
 v0 = vxyz_ptmass(1:3,  wind_emitting_sink)

 phase = omega_pulsation * time + phi0
 r_dot = piston_velocity * cos(phase)

 if (allocated(r_boundary_equilibrium) .and. n_boundary_particles > 0) then
    r_inject = r_boundary_equilibrium(n_boundary_particles) + delta_r_radial(iboundary_spheres + 1)
    r_inject = r_inject + deltaR_osc * sin(phase)
 else
    r_inject = r_min
    do i = 1, iboundary_spheres + 1
       r_inject = r_inject + delta_r_radial(i)
    enddo
    r_inject = r_inject + deltaR_osc * sin(phase)
 endif

 call interp_stellar_profile(r_inject, rho, P, u, T)

 old_npart    = npart
 n_reinjections = n_reinjections + 1

 call inject_fibonacci_sphere(n_shells_total + n_reinjections, npart + 1, particles_to_inject, &
                               r_inject, r_dot, u, rho, &
                               npart, npartoftype, xyzh, vxyzu, igas, x0, v0)

 ! Reinjected particles are gas particles
 mass_injected = real(npart - old_npart) * mass_of_gas_particle
 xyzmh_ptmass(4, wind_emitting_sink) = xyzmh_ptmass(4, wind_emitting_sink) - mass_injected

 print *, ''
 print *, ' Particles injected         :', (npart - old_npart)
 print *, ' Injection radius           :', r_inject
 print *, ' New total particles        :', npart
 print *, 'Reinjection complete.'
 print *, ''

end subroutine perform_reinjection

!-----------------------------------------------------------------------
!+
!  Setup initial atmosphere with all shells at t=0
!+
!-----------------------------------------------------------------------
subroutine setup_initial_atmosphere(xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npartoftype)
 use part,           only:igas,iboundary,iphase,iamtype
 use injectutils,    only:inject_fibonacci_sphere
 use wind_pulsating, only:interp_stellar_profile
 use physcon,        only:pi,km,au

 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real,    intent(in)    :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)

 integer :: i, j, first_particle, ipart_type, nboundary
 real    :: r, r_cur, rho, u, T, P, x0(3), v0(3), v_radial
 logical :: is_boundary

 x0 = xyzmh_ptmass(1:3, wind_emitting_sink)
 v0 = vxyz_ptmass(1:3,  wind_emitting_sink)

 npart = 0

 ! ---- Boundary shells (independent grid, innermost first) ----
 r_cur = shell_radii_bnd(1) - 0.5*delta_r_boundary(1)  ! inner edge of innermost bnd shell
 do i = 1, n_shells_bnd
    r     = shell_radii_bnd(i)
    r_cur = r_cur + delta_r_boundary(i)   ! advance to shell centre
    call interp_stellar_profile(r, rho, P, u, T)
    v_radial      = 0.0
    first_particle = npart + 1
    call inject_fibonacci_sphere(i, first_particle, npart_per_boundary_shell(i), r, v_radial, u, rho, &
                                 npart, npartoftype, xyzh, vxyzu, iboundary, x0, v0)
 enddo

 ! ---- Gas shells (outward from r_min) ----
 do i = 1, n_shells_total
    r = shell_radii_gas(i)
    call interp_stellar_profile(r, rho, P, u, T)
    v_radial      = 0.0
    first_particle = npart + 1
    call inject_fibonacci_sphere(n_shells_bnd + i, first_particle, npart_per_shell(i), r, v_radial, u, rho, &
                                 npart, npartoftype, xyzh, vxyzu, igas, x0, v0)
 enddo

 nboundary           = npartoftype(iboundary)
 n_boundary_particles = nboundary

 print *, 'Boundary particles : ', nboundary
 print *, 'Gas particles      : ', npartoftype(igas)

 if (nboundary > 0) then
    allocate(r_boundary_equilibrium(nboundary))
    allocate(boundary_particle_ids(nboundary))

    j = 0
    do i = 1, npart
       if (j >= nboundary) exit
       if (iamtype(iphase(i)) == iboundary) then
          j = j + 1
          boundary_particle_ids(j) = i
          r_boundary_equilibrium(j) = sqrt( (xyzh(1,i)-x0(1))**2 + &
                                            (xyzh(2,i)-x0(2))**2 + &
                                            (xyzh(3,i)-x0(3))**2 ) &
                                      - deltaR_osc * sin(phi0)
       endif
    enddo
 endif

end subroutine setup_initial_atmosphere

!-----------------------------------------------------------------------
!+
!  Reconstruct boundary particle information when resuming from dump
!+
!-----------------------------------------------------------------------
subroutine reconstruct_boundary_info(time,xyzh,npart,xyzmh_ptmass)
 use part,   only:iboundary,iphase,iamtype
 use physcon,only:pi

 real,    intent(in) :: time
 real,    intent(in) :: xyzh(:,:),xyzmh_ptmass(:,:)
 integer, intent(in) :: npart
 integer :: i, j
 real    :: x0(3), r_current, phase

 x0    = xyzmh_ptmass(1:3, wind_emitting_sink)
 phase = omega_pulsation * time + phi0

 n_boundary_particles = 0
 do i = 1, npart
    if (iamtype(iphase(i)) == iboundary) n_boundary_particles = n_boundary_particles + 1
 enddo

 if (n_boundary_particles > 0) then
    allocate(r_boundary_equilibrium(n_boundary_particles))
    allocate(boundary_particle_ids(n_boundary_particles))

    j = 0
    do i = 1, npart
       if (iamtype(iphase(i)) == iboundary) then
          j = j + 1
          boundary_particle_ids(j) = i
          r_current = sqrt((xyzh(1,i)-x0(1))**2 + &
                           (xyzh(2,i)-x0(2))**2 + &
                           (xyzh(3,i)-x0(3))**2)
          r_boundary_equilibrium(j) = r_current - deltaR_osc * sin(phase)
       endif
    enddo

    print *, 'Reconstructed boundary particle info:'
    print *, 'Boundary particles: ', n_boundary_particles
 endif

end subroutine reconstruct_boundary_info

!-----------------------------------------------------------------------
!+
!  Apply pulsation to boundary particles
!+
!-----------------------------------------------------------------------
subroutine apply_pulsation(time,xyzh,vxyzu,npart,xyzmh_ptmass,vxyz_ptmass)
 use physcon,        only:pi
 use wind_pulsating, only:interp_stellar_profile

 real,    intent(in)    :: time
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(in)    :: npart

 integer :: i, ipart
 real    :: r_eq, r_new, r_current, phase, piston_velocity_n, deltaR_osc_n
 real    :: x_hat(3), r_dot, x0(3), v0(3)
 real    :: x, y, z, rho, u, T, P

 if (.not. allocated(boundary_particle_ids)) return
 if (n_boundary_particles == 0) return

 x0 = xyzmh_ptmass(1:3, wind_emitting_sink)
 v0 = vxyz_ptmass(1:3,  wind_emitting_sink)

 phase = omega_pulsation * time + phi0

 if (time < pulsation_period * time_puls .and. time_puls > 0) then
    piston_velocity_n = time * piston_velocity / (pulsation_period * time_puls)
 else
    piston_velocity_n = piston_velocity
 endif

 deltaR_osc_n = pulsation_period * piston_velocity_n / (2.0*acos(-1.0))
 r_dot        = piston_velocity_n * cos(phase)

 do i = 1, n_boundary_particles
    ipart = boundary_particle_ids(i)
    r_eq  = r_boundary_equilibrium(i)
    r_new = r_eq + deltaR_osc_n * sin(phase)

    x = xyzh(1,ipart) - x0(1)
    y = xyzh(2,ipart) - x0(2)
    z = xyzh(3,ipart) - x0(3)
    r_current = sqrt(x**2 + y**2 + z**2)

    x_hat(1) = x / r_current
    x_hat(2) = y / r_current
    x_hat(3) = z / r_current

    xyzh(1,ipart) = r_new * x_hat(1) + x0(1)
    xyzh(2,ipart) = r_new * x_hat(2) + x0(2)
    xyzh(3,ipart) = r_new * x_hat(3) + x0(3)

    vxyzu(1,ipart) = r_dot * x_hat(1) + v0(1)
    vxyzu(2,ipart) = r_dot * x_hat(2) + v0(2)
    vxyzu(3,ipart) = r_dot * x_hat(3) + v0(3)

    if (var_boundary) then
       call interp_stellar_profile(r_new, rho, P, u, T)
       vxyzu(4,ipart) = u
       xyzh(4,ipart)  = (mass_of_boundary_particle / rho)**(1./3.)
    endif
 enddo

end subroutine apply_pulsation

subroutine update_injected_par
 ! placeholder
end subroutine update_injected_par

!-----------------------------------------------------------------------
!+
!  Write mass-loss rate data to file for restart
!+
!-----------------------------------------------------------------------
subroutine write_mass_loss_data()
 use io, only:iprint
 integer :: iunit, ierr, i

 if (.not. mass_loss_rate_calculated) return

 open(newunit=iunit, file='mass_loss_rate.dat', status='replace', iostat=ierr)
 if (ierr /= 0) then
    write(iprint,*) 'Could not write mass_loss_rate.dat'
    return
 endif

 write(iunit,*) '# Mass-loss rate data for restart'
 write(iunit,*) mass_loss_rate_calculated
 write(iunit,*) mean_mass_loss_rate
 write(iunit,*) Mtotal
 write(iunit,*) particles_to_inject
 write(iunit,*) n_measurements
 write(iunit,*) mass_of_gas_particle
 write(iunit,*) mass_of_boundary_particle
 do i = 1, n_measurements
    write(iunit,*) mass_loss_rates(i)
 enddo
 close(iunit)

 write(iprint,*) 'Mass-loss rate data written to mass_loss_rate.dat'

end subroutine write_mass_loss_data

!-----------------------------------------------------------------------
!+
!  Read mass-loss rate data from file for restart
!+
!-----------------------------------------------------------------------
subroutine read_mass_loss_data()
 use io, only:iprint
 integer :: iunit, ierr, i
 logical :: file_exists

 inquire(file='mass_loss_rate.dat', exist=file_exists)
 if (.not. file_exists) return

 open(newunit=iunit, file='mass_loss_rate.dat', status='old', iostat=ierr)
 if (ierr /= 0) return

 read(iunit,*)   ! skip comment
 read(iunit,*, iostat=ierr) mass_loss_rate_calculated
 if (ierr /= 0) then; close(iunit); return; endif
 read(iunit,*, iostat=ierr) mean_mass_loss_rate
 read(iunit,*, iostat=ierr) Mtotal
 read(iunit,*, iostat=ierr) particles_to_inject
 read(iunit,*, iostat=ierr) n_measurements
 read(iunit,*, iostat=ierr) mass_of_gas_particle
 read(iunit,*, iostat=ierr) mass_of_boundary_particle

 if (n_measurements > 0) then
    if (.not. allocated(mass_loss_rates)) allocate(mass_loss_rates(n_measurements))
    do i = 1, n_measurements
       read(iunit,*, iostat=ierr) mass_loss_rates(i)
       if (ierr /= 0) exit
    enddo
 endif
 close(iunit)

 write(iprint,*) 'Mass-loss rate data read from mass_loss_rate.dat'
 write(iprint,*) ' Mean mass-loss rate          :', mean_mass_loss_rate
 write(iprint,*) ' Gas particle mass            :', mass_of_gas_particle
 write(iprint,*) ' Boundary particle mass       :', mass_of_boundary_particle
 write(iprint,*) ' Particles to inject          :', particles_to_inject

end subroutine read_mass_loss_data

!-----------------------------------------------------------------------
!+
!  Calculate pulsation period from stellar mass-radius relation
!+
!-----------------------------------------------------------------------
subroutine calculate_period(M, R, pulsation_period_days)
 real, intent(in)  :: M, R
 real, intent(out) :: pulsation_period_days
 real :: logP, logM, logR

 logM = log10(M)
 logR = log10(R * 215.032)
 logP = -1.92 - 0.73*logM + 1.86*logR
 pulsation_period_days = 10.0**logP

 print *, 'Calculated pulsation period (days): ', pulsation_period_days

end subroutine calculate_period

!-----------------------------------------------------------------------
!+
!  Write options to input file
!+
!-----------------------------------------------------------------------
subroutine write_options_inject(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(n_profile_points,        'n_profile_points',          'number of points in stellar profile',iunit)
 call write_inopt(iboundary_spheres,       'iboundary_spheres',         'number of boundary spheres (piston layers)',iunit)
 call write_inopt(n_particles,             'n_particles',               'target total gas particles',iunit)
 call write_inopt(n_shells,                'n_shells','number of gas shells (if <0 determined from n_particles)',iunit)
 call write_inopt(boundary_fraction,       'boundary_fraction', 'ratio N_boundary_per_shell/N_gas_per_shell at interface',iunit)
 call write_inopt(rho_power_in,            'rho_power','density profile exponent: rho ~ r^(-rho_power)',iunit)
 call write_inopt(r_b_min_on_rmax, 'r_b_min_on_rmax','inner edge of boundary region as fraction of R_star',iunit)
 call write_inopt(r_min_on_rstar,          'r_min_on_rstar',         'gas atmosphere inner radius as fraction of R_star',iunit)
 call write_inopt(r_max_on_rstar,          'r_max_on_rstar',         'gas atmosphere outer radius as fraction of R_star',iunit)
 call write_inopt(atmos_mass_fraction, 'atmos_mass_fraction', 'atmospheric mass as fraction of total stellar mass',iunit)
 call write_inopt(surface_pressure,    'surface_pressure',    'surface pressure (cgs)',iunit)
 call write_inopt(iwind,               'iwind',               'wind type: 1=prescribed, 2=period from mass-radius relation',iunit)
 call write_inopt(pulsation_period_days,'pulsation_period',   'pulsation period (days)',iunit)
 call write_inopt(piston_velocity_km_s,'piston_velocity',     'piston velocity amplitude (km/s)',iunit)
 call write_inopt(time_puls,           'time_puls',           'time for piston to ramp up (in periods, -1=instant)',iunit)
 call write_inopt(pulsation_timestep,  'pulsation_timestep',  'pulsation timestep as fraction of period',iunit)
 call write_inopt(phi0,                'phi0',                'initial phase offset (radians)',iunit)
 call write_inopt(wss,                 'wss',                 'radial/tangential spacing ratio',iunit)
 call write_inopt(var_boundary,        'var_boundary',        'update boundary thermo with pulsation (logical)',iunit)
 call write_inopt(reinject_enabled,    'reinject_enabled',    'enable dynamic reinjection (logical)',iunit)
 call write_inopt(reinject_period_days,'reinject_period_days','period between reinjections (days)',iunit)
 call write_inopt(mass_loss_start,     'mass_loss_start',     'start time for mass-loss calculation (years)',iunit)
 call write_inopt(mass_loss_end,       'mass_loss_end',       'end time for mass-loss calculation (years)',iunit)
 call write_inopt(check_radius_au,     'check_radius_au',     'mass-loss counting radius (AU)',iunit)
 call write_inopt(meas_int_days,       'meas_int_days',       'mass measurement interval (days)',iunit)

end subroutine write_options_inject

!-----------------------------------------------------------------------
!+
!  Read options from input file
!+
!-----------------------------------------------------------------------
subroutine read_options_inject(name,valstring,imatch,igotall,ierr)
 use io, only:fatal
 character(len=*), intent(in)  :: name, valstring
 logical,          intent(out) :: imatch, igotall
 integer,          intent(out) :: ierr

 integer, save      :: ngot = 0
 integer, parameter :: noptions = 25
 logical :: init_opt = .false.

 if (.not. init_opt) then
    init_opt = .true.
    call set_default_options_inject()
 endif

 imatch  = .true.
 igotall = .false.

 select case(trim(name))
 case('n_profile_points')
    read(valstring,*,iostat=ierr) n_profile_points
    ngot = ngot + 1
    if (n_profile_points <= 10) call fatal(label,'n_profile_points must be > 10')
 case('iboundary_spheres')
    read(valstring,*,iostat=ierr) iboundary_spheres
    ngot = ngot + 1
    if (iboundary_spheres < 0) call fatal(label,'iboundary_spheres must be >= 0')
 case('n_particles')
    read(valstring,*,iostat=ierr) n_particles
    ngot = ngot + 1
    if (n_particles < 1) call fatal(label,'n_particles must be >= 1')
 case('n_shells')
    read(valstring,*,iostat=ierr) n_shells
    ngot = ngot + 1
    if (n_shells < -10) call fatal(label,'n_shells must be >= -10')
 case('boundary_fraction')
    read(valstring,*,iostat=ierr) boundary_fraction
    ngot = ngot + 1
    if (boundary_fraction <= 0. .or. boundary_fraction > 1.0) &
       call fatal(label,'boundary_fraction must be in (0,1]')
 case('rho_power')
    read(valstring,*,iostat=ierr) rho_power_in
    ngot = ngot + 1
    if (rho_power_in <= 0.) call fatal(label,'rho_power must be > 0')
 case('r_min_on_rstar')
    read(valstring,*,iostat=ierr) r_min_on_rstar
    ngot = ngot + 1
    if (r_min_on_rstar <= 0. .or. r_min_on_rstar >= 1.0) &
       call fatal(label,'r_min_on_rstar must be in (0,1)')
 case('r_b_min_on_rmax')
    read(valstring,*,iostat=ierr) r_b_min_on_rmax
    ngot = ngot + 1
    if (r_b_min_on_rmax <= 0. .or. r_b_min_on_rmax >= 1.0) &
       call fatal(label,'r_b_min_on_rmax must be in (0,1)')
    if (r_b_min_on_rmax >= r_min_on_rstar) &
       call fatal(label,'r_b_min_on_rmax must be less than r_min_on_rstar')
 case('r_max_on_rstar')
    read(valstring,*,iostat=ierr) r_max_on_rstar
    ngot = ngot + 1
    if (r_max_on_rstar <= 0. .or. r_max_on_rstar > 10.0) &
       call fatal(label,'r_max_on_rstar must be in (0,10]')
 case('atmos_mass_fraction')
    read(valstring,*,iostat=ierr) atmos_mass_fraction
    ngot = ngot + 1
    if (atmos_mass_fraction <= 0. .or. atmos_mass_fraction >= 1.0) &
       call fatal(label,'atmos_mass_fraction must be in (0,1)')
 case('surface_pressure')
    read(valstring,*,iostat=ierr) surface_pressure
    ngot = ngot + 1
    if (surface_pressure < 0.) call fatal(label,'surface_pressure must be >= 0')
 case('iwind')
    read(valstring,*,iostat=ierr) iwind
    ngot = ngot + 1
    if (iwind /= 1 .and. iwind /= 2) call fatal(label,'iwind must be 1 or 2')
 case('pulsation_period')
    read(valstring,*,iostat=ierr) pulsation_period_days
    ngot = ngot + 1
    if (pulsation_period_days < 0.) call fatal(label,'pulsation_period must be >= 0')
 case('piston_velocity')
    read(valstring,*,iostat=ierr) piston_velocity_km_s
    ngot = ngot + 1
    if (piston_velocity_km_s < 0.) call fatal(label,'piston_velocity must be >= 0')
 case('time_puls')
    read(valstring,*,iostat=ierr) time_puls
    ngot = ngot + 1
    if (time_puls < -1) call fatal(label,'time_puls must be >= -1')
 case('pulsation_timestep')
    read(valstring,*,iostat=ierr) pulsation_timestep
    ngot = ngot + 1
    if (pulsation_timestep <= 0. .or. pulsation_timestep > 1.0) &
       call fatal(label,'pulsation_timestep must be in (0,1]')
 case('phi0')
    read(valstring,*,iostat=ierr) phi0
    ngot = ngot + 1
    if (phi0 < -3.1415926536d0 .or. phi0 > 3.1415926536d0) &
       call fatal(label,'phi0 must be in (-pi,pi)')
 case('wss')
    read(valstring,*,iostat=ierr) wss
    ngot = ngot + 1
    if (wss <= 0. .or. wss > 10.0) call fatal(label,'wss must be in (0,10]')
 case('var_boundary')
    read(valstring,*,iostat=ierr) var_boundary
    ngot = ngot + 1
 case('reinject_enabled')
    read(valstring,*,iostat=ierr) reinject_enabled
    ngot = ngot + 1
 case('reinject_period_days')
    read(valstring,*,iostat=ierr) reinject_period_days
    ngot = ngot + 1
    if (reinject_period_days <= 0.) call fatal(label,'reinject_period_days must be > 0')
 case('mass_loss_start')
    read(valstring,*,iostat=ierr) mass_loss_start
    ngot = ngot + 1
    if (mass_loss_start < 0.) call fatal(label,'mass_loss_start must be >= 0')
 case('mass_loss_end')
    read(valstring,*,iostat=ierr) mass_loss_end
    ngot = ngot + 1
    if (mass_loss_end <= 0.) call fatal(label,'mass_loss_end must be > 0')
 case('check_radius_au')
    read(valstring,*,iostat=ierr) check_radius_au
    ngot = ngot + 1
    if (check_radius_au <= 0.) call fatal(label,'check_radius_au must be > 0')
 case('meas_int_days')
    read(valstring,*,iostat=ierr) meas_int_days
    ngot = ngot + 1
    if (meas_int_days <= 0.) call fatal(label,'meas_int_days must be > 0')
 case default
    imatch = .false.
 end select

 igotall = (ngot >= noptions)

end subroutine read_options_inject

end module inject