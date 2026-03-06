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

 integer :: iboundary_spheres     = 5
 integer :: n_profile_points      = 10000
 integer :: n_particles           = 500000
 integer :: n_shells              = 25
 real    :: boundary_fraction     = 0.9
 real    :: rho_power_in          = 10.0
 real    :: r_min_on_rstar        = 0.9
 real    :: r_max_on_rstar        = 1.4
 real    :: r_b_min_on_rmax       = 0.8
 real    :: dtpulsation           = huge(0.)
 real    :: pulsation_period_days = 300.0
 real    :: piston_velocity_km_s  = 4.0
 real    :: time_puls             = -1.0
 real    :: atmos_mass_fraction   = 5e-5
 real    :: surface_pressure      = 0.001
 integer :: iwind                 = 1
 real    :: pulsation_timestep    = 0.02
 real    :: phi0                  = 3.1415926536d0/2.0
 real    :: wss                   = 1.0
 logical :: var_boundary          = .false.

 logical :: reinject_enabled      = .true.
 real    :: reinject_period_days  = 10.0
 real    :: mass_loss_start       = 1.0
 real    :: mass_loss_end         = 3.0
 real    :: check_radius_au       = 3.0
 real    :: meas_int_days         = 10.0

 integer, parameter :: wind_emitting_sink = 1
 integer, parameter :: max_measurements   = 10000

 real :: omega_pulsation, deltaR_osc, pulsation_period, piston_velocity
 real :: Rstar, r_min, r_boundary_min
 real :: Mtotal, Matmos, Msink

 integer, allocatable :: npart_per_shell(:)
 integer, allocatable :: npart_per_boundary_shell(:)

 real, allocatable :: delta_r_gas(:)
 real, allocatable :: delta_r_boundary(:)
 real, allocatable :: shell_radii_gas(:)
 real, allocatable :: shell_radii_bnd(:)

 real :: mass_of_gas_particle      = 0.0
 real :: mass_of_boundary_particle = 0.0

 real, allocatable :: delta_r_radial(:)

 logical :: atmosphere_setup_complete = .false.
 integer :: n_shells_total
 integer :: n_shells_bnd

 real, allocatable    :: r_boundary_equilibrium(:)
 integer, allocatable :: boundary_particle_ids(:)
 integer              :: n_boundary_particles
 integer              :: active_boundary_spheres

 logical :: reinjection_needed = .false.

 real    :: time_last_reinject = 0.0
 real    :: reinject_period
 integer :: n_reinjections     = 0

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

subroutine set_default_options_inject(flag)
 integer, optional, intent(in) :: flag

 iboundary_spheres     = 5
 n_profile_points      = 10000
 n_particles           = 500000
 n_shells              = 25
 boundary_fraction     = 0.9
 rho_power_in          = 10.0
 r_min_on_rstar        = 0.9
 r_max_on_rstar        = 1.4
 r_b_min_on_rmax       = 0.8
 dtpulsation           = huge(0.)
 atmos_mass_fraction   = 5e-5
 surface_pressure      = 0.001
 iwind                 = 1
 pulsation_period_days = 300.0
 piston_velocity_km_s  = 4.0
 time_puls             = -1.0
 pulsation_timestep    = 0.02
 phi0                  = 3.1415926536d0/2.0
 wss                   = 1.0
 var_boundary          = .false.
 reinject_enabled      = .true.
 reinject_period_days  = 10.0
 mass_loss_start       = 1.0
 mass_loss_end         = 3.0
 check_radius_au       = 3.0
 meas_int_days         = 10.0

end subroutine set_default_options_inject

subroutine init_inject(ierr)
 use io,            only:fatal
 use physcon,       only:pi,days,au,solarm,km,years
 use icosahedron,   only:compute_matrices,compute_corners
 use eos,           only:gmw,gamma
 use units,         only:utime,umass,unit_velocity,unit_luminosity
 use part,          only:xyzmh_ptmass,massoftype,igas,iboundary,nptmass,iTeff,iReff,iLum
 use injectutils,   only:get_parts_per_sphere, get_fibonacci_spacing
 use wind_pulsating,only:setup_star,calc_stellar_profile,region_mass,interp_stellar_profile
 use dust_formation,only:calc_kappa_max

 integer, intent(out) :: ierr
 real    :: Mstar_cgs, Rstar_cgs, Tstar, Lstar_cgs
 real    :: current_radius, dr, dr_new, r_c, rho, P, u, T, m_shell
 integer :: shell_index, max_shells, n_first, n_shell
 integer :: expected_measurements, i, iter, bisect_iter
 integer, parameter  :: max_shells_tmp = 2000
 integer, parameter  :: max_iter_dr    = 100
 real,    parameter  :: tol_dr         = 1.0e-6
 real    :: tmp_dr_gas(max_shells_tmp), tmp_r_gas(max_shells_tmp)
 integer :: tmp_n_gas(max_shells_tmp)
 real    :: tmp_dr_bnd(max_shells_tmp), tmp_r_bnd(max_shells_tmp)
 integer :: tmp_n_bnd(max_shells_tmp)
 real    :: m_bnd, m_bnd_lo, m_bnd_hi, m_bnd_mid
 integer :: n_bnd_count
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

 r_min          = r_min_on_rstar * Rstar
 r_boundary_min = r_b_min_on_rmax * Rstar

 if (r_boundary_min >= r_min) &
    call fatal(label,'r_b_min_on_rmax must be less than r_min_on_rstar')

 call calc_stellar_profile(n_profile_points)

 ! ================================================================
 ! GAS GRID: single forward-building loop from r_min.
 !
 ! Free parameter: n_first (particle count of the first shell).
 ! Outer loop increments n_first until the shell count or total
 ! particle budget target is reached.
 !
 ! Per shell, a fixed-point inner loop finds the self-consistent dr:
 !   N_shell = max(1, round(4*pi * rho(r_c) * r_c^2 * dr / m_gas))
 !   dr_new  = wss * r_inner * fibonacci(N_shell)
 ! where m_gas = M([r_min, r_max]) / n_particles (for n_shells<0)
 ! or is determined after convergence (for n_shells>0).
 !
 ! For n_shells>0: m_gas is not known a priori, so we use n_first
 ! directly as the particle count of the first shell and propagate
 ! forward using the density ratio between adjacent shell centres.
 ! ================================================================

 if (n_shells < 0) then
    ! m_gas is known: total mass / target particle count
    mass_of_gas_particle = region_mass(r_min, r_max_on_rstar * Rstar) / real(n_particles)

    n_first   = 1
    converged = .false.
    do while (.not. converged)
       n_first        = n_first + 1
       current_radius = r_min
       shell_index    = 0

       do while (current_radius < r_max_on_rstar * Rstar)
          shell_index = shell_index + 1
          if (shell_index > max_shells_tmp) call fatal(label,'gas shell tmp array too small')

          ! fixed-point loop for self-consistent dr
          dr = wss * current_radius * get_fibonacci_spacing(n_first)
          do iter = 1, max_iter_dr
             r_c    = current_radius + 0.5 * dr
             call interp_stellar_profile(r_c, rho, P, u, T)
             m_shell = 4.0 * pi * rho * r_c**2 * dr
             n_shell = max(1, nint(m_shell / mass_of_gas_particle))
             dr_new  = wss * current_radius * get_fibonacci_spacing(n_shell)
             if (abs(dr_new - dr) < tol_dr * dr) then
                dr = dr_new
                exit
             endif
             dr = dr_new
          enddo

          tmp_dr_gas(shell_index) = dr
          tmp_r_gas(shell_index)  = current_radius + 0.5 * dr
          tmp_n_gas(shell_index)  = n_shell
          current_radius          = current_radius + dr
       enddo

       if (sum(tmp_n_gas(1:shell_index)) >= n_particles) converged = .true.
    enddo

 else
    ! m_gas not known a priori; use n_first as first-shell count and
    ! propagate N forward via density ratio between shell centres.
    n_first   = 1
    converged = .false.
    do while (.not. converged)
       n_first        = n_first + 1
       current_radius = r_min
       shell_index    = 0

       do while (current_radius < r_max_on_rstar * Rstar)
          shell_index = shell_index + 1
          if (shell_index > max_shells_tmp) call fatal(label,'gas shell tmp array too small')

          if (shell_index == 1) then
             n_shell = n_first
          else
             ! propagate N using density ratio: N_i/N_{i-1} = rho_i*r_i^2*dr_i / (rho_{i-1}*r_{i-1}^2*dr_{i-1})
             ! approximated at shell centres with equal-mass condition
             call interp_stellar_profile(tmp_r_gas(shell_index-1), rho, P, u, T)
             m_shell = 4.0 * pi * rho * tmp_r_gas(shell_index-1)**2 * tmp_dr_gas(shell_index-1)
             ! m_gas estimate from previous shell
             mass_of_gas_particle = m_shell / real(tmp_n_gas(shell_index-1))
          endif

          ! fixed-point loop for self-consistent dr
          dr = wss * current_radius * get_fibonacci_spacing(n_shell)
          do iter = 1, max_iter_dr
             r_c    = current_radius + 0.5 * dr
             call interp_stellar_profile(r_c, rho, P, u, T)
             m_shell = 4.0 * pi * rho * r_c**2 * dr
             if (shell_index > 1) then
                n_shell = max(1, nint(m_shell / mass_of_gas_particle))
             endif
             dr_new  = wss * current_radius * get_fibonacci_spacing(n_shell)
             if (abs(dr_new - dr) < tol_dr * dr) then
                dr = dr_new
                exit
             endif
             dr = dr_new
          enddo

          tmp_dr_gas(shell_index) = dr
          tmp_r_gas(shell_index)  = current_radius + 0.5 * dr
          tmp_n_gas(shell_index)  = n_shell
          current_radius          = current_radius + dr
       enddo

       if (shell_index >= max_shells) converged = .true.
    enddo

    ! Now set m_gas from the converged grid
    mass_of_gas_particle = region_mass(r_min, r_max_on_rstar * Rstar) / real(sum(tmp_n_gas(1:shell_index)))
 endif

 n_shells_total = shell_index

 allocate(npart_per_shell(n_shells_total))
 allocate(delta_r_gas(n_shells_total))
 allocate(shell_radii_gas(n_shells_total))
 do i = 1, n_shells_total
    npart_per_shell(i)  = tmp_n_gas(i)
    delta_r_gas(i)      = tmp_dr_gas(i)
    shell_radii_gas(i)  = tmp_r_gas(i)
 enddo

 ! ================================================================
 ! BOUNDARY GRID: bisect on m_bnd until exactly iboundary_spheres
 ! shells fit stepping inward from r_min.
 !
 ! Higher m_bnd -> wider dr -> fewer shells; lower -> more shells.
 ! Bracket: find m_bnd_lo (too many shells) and m_bnd_hi (too few),
 ! then bisect to 60 iterations (~2^-60 relative precision).
 !
 ! Interface condition: the gap between the outermost boundary shell
 ! centre and the innermost gas shell centre equals
 !   0.5*dr_bnd(n_shells_bnd) + 0.5*dr_gas(1)
 ! which is automatically satisfied since both grids tile from r_min.
 ! ================================================================
 n_shells_bnd = iboundary_spheres

 if (n_shells_bnd < 1) then
    allocate(npart_per_boundary_shell(0))
    allocate(delta_r_boundary(0))
    allocate(shell_radii_bnd(0))
 else
    ! Initial bracket anchor from gas particle mass
    m_bnd_lo = mass_of_gas_particle * 1.0e-6
    m_bnd_hi = mass_of_gas_particle * 1.0e6

    ! Widen lo until it gives strictly more than iboundary_spheres shells
    do bisect_iter = 1, 200
       if (count_bnd_shells(m_bnd_lo, max_shells_tmp, tmp_dr_bnd, tmp_r_bnd, tmp_n_bnd, &
                             r_min, r_boundary_min, pi, tol_dr, max_iter_dr) > iboundary_spheres) exit
       m_bnd_lo = m_bnd_lo * 0.1
    enddo
    ! Widen hi until it gives strictly fewer than iboundary_spheres shells
    do bisect_iter = 1, 200
       if (count_bnd_shells(m_bnd_hi, max_shells_tmp, tmp_dr_bnd, tmp_r_bnd, tmp_n_bnd, &
                             r_min, r_boundary_min, pi, tol_dr, max_iter_dr) < iboundary_spheres) exit
       m_bnd_hi = m_bnd_hi * 10.0
    enddo

    ! Bisect: maintain invariant lo gives >=iboundary_spheres, hi gives <iboundary_spheres
    do bisect_iter = 1, 60
       m_bnd_mid   = 0.5 * (m_bnd_lo + m_bnd_hi)
       n_bnd_count = count_bnd_shells(m_bnd_mid, max_shells_tmp, tmp_dr_bnd, tmp_r_bnd, tmp_n_bnd, &
                                       r_min, r_boundary_min, pi, tol_dr, max_iter_dr)
       if (n_bnd_count >= iboundary_spheres) then
          m_bnd_lo = m_bnd_mid
       else
          m_bnd_hi = m_bnd_mid
       endif
    enddo
    m_bnd = m_bnd_lo  ! m_bnd_lo is the largest m_bnd that still gives >= iboundary_spheres shells
    mass_of_boundary_particle = m_bnd

    ! Final pass with converged m_bnd to get shell arrays
    n_bnd_count = count_bnd_shells(m_bnd, max_shells_tmp, tmp_dr_bnd, tmp_r_bnd, tmp_n_bnd, &
                                    r_min, r_boundary_min, pi, tol_dr, max_iter_dr)
    n_shells_bnd = min(n_shells_bnd, n_bnd_count)

    allocate(npart_per_boundary_shell(n_shells_bnd))
    allocate(delta_r_boundary(n_shells_bnd))
    allocate(shell_radii_bnd(n_shells_bnd))

    ! Reverse: tmp arrays are outermost->innermost; store innermost->outermost
    do i = 1, n_shells_bnd
       delta_r_boundary(i)         = tmp_dr_bnd(n_shells_bnd + 1 - i)
       shell_radii_bnd(i)          = tmp_r_bnd(n_shells_bnd + 1 - i)
       npart_per_boundary_shell(i) = tmp_n_bnd(n_shells_bnd + 1 - i)
    enddo
 endif

 ! Combined delta_r_radial (boundary first, then gas)
 if (allocated(delta_r_radial)) deallocate(delta_r_radial)
 allocate(delta_r_radial(n_shells_bnd + n_shells_total))
 if (n_shells_bnd > 0) delta_r_radial(1:n_shells_bnd) = delta_r_boundary
 delta_r_radial(n_shells_bnd+1 : n_shells_bnd+n_shells_total) = delta_r_gas

 ! Particle masses
 mass_of_gas_particle = region_mass(r_min, r_max_on_rstar * Rstar) / real(sum(npart_per_shell))
 if (n_shells_bnd > 0) then
    mass_of_boundary_particle = region_mass(r_boundary_min, r_min) / real(sum(npart_per_boundary_shell))
 else
    mass_of_boundary_particle = mass_of_gas_particle
 endif

 massoftype(igas)      = mass_of_gas_particle
 massoftype(iboundary) = mass_of_boundary_particle

 if (file_exists) call read_mass_loss_data()

 print *, ''
 print *, ' rho_power                        :', rho_power_in
 print *, ' Boundary region  [r_bnd_min, r_min] / Rstar :', r_b_min_on_rmax, r_min_on_rstar
 print *, ' Gas region       [r_min,     r_max] / Rstar :', r_min_on_rstar, r_max_on_rstar
 print *, ' M_gas   (Msun)                   :', region_mass(r_min,          r_max_on_rstar * Rstar)
 print *, ' M_bnd   (Msun)                   :', region_mass(r_boundary_min, r_min)
 print *, ' Gas shells                       :', n_shells_total
 print *, ' Boundary shells                  :', n_shells_bnd
 print *, ' Total gas particles              :', sum(npart_per_shell)
 if (n_shells_bnd > 0) then
    print *, ' Total boundary particles         :', sum(npart_per_boundary_shell)
    print *, ' N_bnd / N_gas (actual)           :', &
              real(sum(npart_per_boundary_shell)) / real(sum(npart_per_shell))
    print *, ' Innermost boundary N_per_shell   :', npart_per_boundary_shell(1)
    print *, ' Outermost boundary N_per_shell   :', npart_per_boundary_shell(n_shells_bnd)
    print *, ' Innermost gas     N_per_shell    :', npart_per_shell(1)
 endif
 print *, ' Outermost gas N_per_shell        :', npart_per_shell(n_shells_total)
 print *, ' Gas particle mass (Msun)         :', mass_of_gas_particle
 print *, ' Boundary particle mass (Msun)    :', mass_of_boundary_particle
 print *, ' Boundary/gas mass ratio          :', mass_of_boundary_particle / mass_of_gas_particle
 print *, ''

 do i = 1, n_shells_bnd
    print *, 'Boundary shell ', i, ': r=', shell_radii_bnd(i)/Rstar, ' Rstar, dr=', delta_r_boundary(i)/Rstar, &
             ' Rstar, N_particles=', npart_per_boundary_shell(i)
 enddo
 do i = 1, n_shells_total
    print *, 'Gas shell      ', i, ': r=', shell_radii_gas(i)/Rstar, ' Rstar, dr=', delta_r_gas(i)/Rstar, &
             ' Rstar, N_particles=', npart_per_shell(i)
 enddo

end subroutine init_inject

! Helper: count how many boundary shells fit stepping inward from r_min
! with a given m_bnd. Also fills tmp arrays for the final pass.
integer function count_bnd_shells(m_bnd, max_shells_tmp, tmp_dr, tmp_r, tmp_n, &
                                   r_min_in, r_bnd_min, pi_in, tol, max_iter)
 use injectutils,    only:get_fibonacci_spacing
 use wind_pulsating, only:interp_stellar_profile

 real,    intent(in)    :: m_bnd, r_min_in, r_bnd_min, pi_in, tol
 integer, intent(in)    :: max_shells_tmp, max_iter
 real,    intent(out)   :: tmp_dr(max_shells_tmp), tmp_r(max_shells_tmp)
 integer, intent(out)   :: tmp_n(max_shells_tmp)

 real    :: current_radius, dr, dr_new, r_c, rho, P, u, T, m_shell
 integer :: n_shell, idx, iter

 current_radius = r_min_in
 idx = 0

 do while (current_radius > r_bnd_min)
    idx = idx + 1
    if (idx > max_shells_tmp) then
       idx = idx - 1
       exit
    endif

    dr = min(wss * current_radius * get_fibonacci_spacing(1), &
             (r_min_in - r_bnd_min) * 0.5)

    do iter = 1, max_iter
       r_c    = current_radius - 0.5 * dr
       if (r_c <= r_bnd_min) then
          dr = (current_radius - r_bnd_min) * 0.9
          r_c = current_radius - 0.5 * dr
       endif
       call interp_stellar_profile(r_c, rho, P, u, T)
       m_shell = 4.0 * pi_in * rho * r_c**2 * dr
       n_shell = max(1, nint(m_shell / m_bnd))
       dr_new  = wss * current_radius * get_fibonacci_spacing(n_shell)
       if (abs(dr_new - dr) < tol * dr) then
          dr = dr_new
          exit
       endif
       dr = dr_new
    enddo

    ! If this shell would overshoot, stop
    if (current_radius - dr < r_bnd_min) then
       if (current_radius - 0.5*dr > r_bnd_min) then
          tmp_dr(idx) = dr
          tmp_r(idx)  = current_radius - 0.5 * dr
          tmp_n(idx)  = n_shell
       else
          idx = idx - 1
       endif
       exit
    endif

    tmp_dr(idx) = dr
    tmp_r(idx)  = current_radius - 0.5 * dr
    tmp_n(idx)  = n_shell
    current_radius = current_radius - dr
 enddo

 count_bnd_shells = idx

end function count_bnd_shells


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
       particles_to_inject = nint((mean_mass_loss_rate * reinject_period) / mass_of_gas_particle)
       if (particles_to_inject < 1) particles_to_inject = 1
       mass_loss_rate_calculated = .true.
       call write_mass_loss_data()
    endif
 endif

end subroutine take_periodic_mass_measurements

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

 old_npart      = npart
 n_reinjections = n_reinjections + 1

 call inject_fibonacci_sphere(n_shells_total + n_reinjections, npart + 1, particles_to_inject, &
                               r_inject, r_dot, u, rho, &
                               npart, npartoftype, xyzh, vxyzu, igas, x0, v0)

 mass_injected = real(npart - old_npart) * mass_of_gas_particle
 xyzmh_ptmass(4, wind_emitting_sink) = xyzmh_ptmass(4, wind_emitting_sink) - mass_injected

 print *, ''
 print *, ' Particles injected         :', (npart - old_npart)
 print *, ' Injection radius           :', r_inject
 print *, ' New total particles        :', npart
 print *, 'Reinjection complete.'
 print *, ''

end subroutine perform_reinjection

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

 r_cur = shell_radii_bnd(1) - 0.5*delta_r_boundary(1)
 do i = 1, n_shells_bnd
    r     = shell_radii_bnd(i)
    r_cur = r_cur + delta_r_boundary(i)
    call interp_stellar_profile(r, rho, P, u, T)
    v_radial       = 0.0
    first_particle = npart + 1
    call inject_fibonacci_sphere(i, first_particle, npart_per_boundary_shell(i), r, v_radial, u, rho, &
                                 npart, npartoftype, xyzh, vxyzu, iboundary, x0, v0)
 enddo

 do i = 1, n_shells_total
    r = shell_radii_gas(i)
    call interp_stellar_profile(r, rho, P, u, T)
    v_radial       = 0.0
    first_particle = npart + 1
    call inject_fibonacci_sphere(n_shells_bnd + i, first_particle, npart_per_shell(i), r, v_radial, u, rho, &
                                 npart, npartoftype, xyzh, vxyzu, igas, x0, v0)
 enddo

 nboundary            = npartoftype(iboundary)
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
end subroutine update_injected_par

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

subroutine read_mass_loss_data()
 use io, only:iprint
 integer :: iunit, ierr, i
 logical :: file_exists

 inquire(file='mass_loss_rate.dat', exist=file_exists)
 if (.not. file_exists) return

 open(newunit=iunit, file='mass_loss_rate.dat', status='old', iostat=ierr)
 if (ierr /= 0) return

 read(iunit,*)
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

subroutine write_options_inject(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(n_profile_points,     'n_profile_points',    'number of points in stellar profile',iunit)
 call write_inopt(iboundary_spheres,    'iboundary_spheres',   'number of boundary spheres (piston layers)',iunit)
 call write_inopt(n_particles,          'n_particles',         'target total gas particles',iunit)
 call write_inopt(n_shells,             'n_shells',            'number of gas shells (if <0 determined from n_particles)',iunit)
 call write_inopt(boundary_fraction,    'boundary_fraction',   'ratio N_boundary_per_shell/N_gas_per_shell at interface',iunit)
 call write_inopt(rho_power_in,         'rho_power',           'density profile exponent: rho ~ r^(-rho_power)',iunit)
 call write_inopt(r_b_min_on_rmax,      'r_b_min_on_rmax',     'inner edge of boundary region as fraction of R_star',iunit)
 call write_inopt(r_min_on_rstar,       'r_min_on_rstar',      'gas atmosphere inner radius as fraction of R_star',iunit)
 call write_inopt(r_max_on_rstar,       'r_max_on_rstar',      'gas atmosphere outer radius as fraction of R_star',iunit)
 call write_inopt(atmos_mass_fraction,  'atmos_mass_fraction', 'atmospheric mass as fraction of total stellar mass',iunit)
 call write_inopt(surface_pressure,     'surface_pressure',    'surface pressure (cgs)',iunit)
 call write_inopt(iwind,                'iwind',               'wind type: 1=prescribed, 2=period from mass-radius relation',iunit)
 call write_inopt(pulsation_period_days,'pulsation_period',    'pulsation period (days)',iunit)
 call write_inopt(piston_velocity_km_s, 'piston_velocity',     'piston velocity amplitude (km/s)',iunit)
 call write_inopt(time_puls,            'time_puls',           'time for piston to ramp up (in periods, -1=instant)',iunit)
 call write_inopt(pulsation_timestep,   'pulsation_timestep',  'pulsation timestep as fraction of period',iunit)
 call write_inopt(phi0,                 'phi0',                'initial phase offset (radians)',iunit)
 call write_inopt(wss,                  'wss',                 'radial/tangential spacing ratio',iunit)
 call write_inopt(var_boundary,         'var_boundary',        'update boundary thermo with pulsation (logical)',iunit)
 call write_inopt(reinject_enabled,     'reinject_enabled',    'enable dynamic reinjection (logical)',iunit)
 call write_inopt(reinject_period_days, 'reinject_period_days','period between reinjections (days)',iunit)
 call write_inopt(mass_loss_start,      'mass_loss_start',     'start time for mass-loss calculation (years)',iunit)
 call write_inopt(mass_loss_end,        'mass_loss_end',       'end time for mass-loss calculation (years)',iunit)
 call write_inopt(check_radius_au,      'check_radius_au',     'mass-loss counting radius (AU)',iunit)
 call write_inopt(meas_int_days,        'meas_int_days',       'mass measurement interval (days)',iunit)

end subroutine write_options_inject

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