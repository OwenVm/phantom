!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module inject
!
! Handles pulsating AGB stars
!
! :References: None
!
! :Owner: Owen Vermeulen
!
! :Runtime parameters:
!   - iboundary_spheres      : *number of boundary spheres (integer)*
!   - n_profile_points       : *number of points in stellar profile calculation (integer)*
!   - min_mass_fraction      : *criteria determining when to stop building shells, i.e., when M_shell / M_tot < mass_fraction*
!   - rho_power_in           : *density profile exponent: rho ~ r^(-rho_power)*
!   - r_min_on_rstar         : *inner radius as fraction of R_star*
!   - dtpulsation            : *pulsation timestep as fraction of pulsation period*
!   - rho_inner              : *density at inner boundary r_min (cgs)*
!   - iwind                  : *wind type: 1=prescribed, 2=period from mass-radius relation*
!   - pulsation_period_days  : *pulsation period (days)*
!   - piston_velocity_km_s   : *piston velocity amplitude (km/s)*
!   - phi0                   : *initial phase offset (radians)*
!   - wss                    : *fraction of tangential and radial distance between particles*
!   - save_period            : *wether to save dumps as an multiple of the pulsation period (0=off, 1=on)*
!   - dumps_p_period         : *how many dumps to save every period, if save_period is activated*
!   - reinject_enabled       : *enable reinjection (logical)*
!   - n_inject_period        : *number of reinjections per period*
!   - mass_loss_start        : *start time for mass-loss calculation in years*
!   - mass_loss_end          : *end time for mass-loss calculation in years*
!   - update_L               : *wether to update the luminosity of the sink particle with the pulsation period (0=off, 1=on)*
!   - use_file_mdot          : *skip measurement phase and use mass_loss_rate.dat directly (0=off, 1=on)*
!
! :Dependencies: dim, eos, icosahedron, infile_utils, injectutils, io,
!   part, partinject, physcon, units
!
 use io, only:fatal
 implicit none
 character(len=*), parameter, public :: inject_type = 'atmosphere'

 public :: init_inject, inject_particles, write_options_inject, read_options_inject, &
           set_default_options_inject, update_injected_par
 private

 integer :: iboundary_spheres        = 5
 integer :: n_profile_points         = 10000
 real    :: min_mass_fraction        = 0.01
 real    :: rho_power_in             = 6.0
 real    :: r_min_on_rstar           = 1.0
 real    :: dtpulsation              = huge(0.)
 real    :: pulsation_period_days    = 300.0
 real    :: piston_velocity_km_s     = 4.0
 real    :: rho_inner                = 1.0e-12
 integer :: iwind                    = 1
 real    :: phi0                     = -3.1415926536d0/2.0
 real    :: wss                      = 1.0
 integer :: save_period              = 0
 integer :: dumps_p_period           = 10

 integer :: reinject_enabled         = 1
 integer :: n_inject_period          = 40
 real    :: mass_loss_start          = 6.0
 real    :: mass_loss_end            = 8.0
 integer :: update_L                 = 0
 integer :: verbose                  = 1
 integer :: use_file_mdot            = 0   ! when 1: read mass_loss_rate.dat at init and skip measurement

 integer, parameter :: wind_emitting_sink = 1
 integer, parameter :: companion_sink     = 2
 integer, parameter :: max_measurements   = 10000

 real :: omega_pulsation, deltaR_osc, pulsation_period, piston_velocity
 real :: Rstar, r_min, r_max, Mstar, Mtotal

 integer, allocatable :: npart_per_shell(:)
 integer, allocatable :: npart_per_boundary_shell(:)

 real, allocatable :: delta_r_gas(:)
 real, allocatable :: delta_r_boundary(:)
 real, allocatable :: shell_radii_gas(:)
 real, allocatable :: shell_radii_bnd(:)

 real :: mass_of_gas_particle      = 0.0   ! in code units; set during init
 real :: pulsation_timestep        = 0.02

 real, allocatable :: delta_r_radial(:)

 logical :: atmosphere_setup_complete = .false.
 integer :: n_shells_total
 integer :: n_shells_bnd

 real, allocatable    :: r_boundary_equilibrium(:)
 integer, allocatable :: boundary_particle_ids(:)
 integer              :: n_boundary_particles

 logical :: reinjection_needed = .false.

 real    :: time_last_reinject = 0.0
 real    :: reinject_period
 integer :: n_reinjections     = 0
 integer :: n_escaping_prev = 0

 real    :: mass_loss_start_time
 real    :: mass_loss_end_time
 real    :: measurement_interval
 real    :: time_next_measurement
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

 iboundary_spheres        = 5
 n_profile_points         = 10000
 min_mass_fraction        = 0.01
 rho_power_in             = 6.0
 r_min_on_rstar           = 1.0
 dtpulsation              = huge(0.)
 rho_inner                = 1.0e-12
 iwind                    = 1
 pulsation_period_days    = 300.0
 piston_velocity_km_s     = 4.0
 phi0                     = -3.1415926536d0/2.0
 wss                      = 1.0
 save_period              = 0
 dumps_p_period           = 10
 reinject_enabled         = 1
 n_inject_period          = 40
 mass_loss_start          = 6.0
 mass_loss_end            = 8.0
 update_L                 = 0
 verbose                  = 1
 use_file_mdot            = 0

end subroutine set_default_options_inject

!----------------------------------------------------------------
!+
!  Derive n_particles_first from the target particle mass and the
!  power-law density profile at r_min.  The first-shell inter-particle
!  spacing dr = wss * r_min * get_fibonacci_spacing(N) implies a shell
!  mass M_shell(N), so we solve N = M_shell(N)/m_p iteratively.
!
!  All arithmetic is done in CGS so that no call to setup_star / region_mass
!  is needed yet (those require r_max which is only known after the loop).
!+
!----------------------------------------------------------------
subroutine derive_n_particles_first(r_min_cgs, rho_inner_cgs, rho_power, &
                                     wss_in, m_particle_cgs, n_first)
 use physcon,     only:pi
 use injectutils, only:get_fibonacci_spacing

 real,    intent(in)  :: r_min_cgs, rho_inner_cgs, rho_power, wss_in, m_particle_cgs
 integer, intent(out) :: n_first

 integer, parameter :: max_iter = 100
 real,    parameter :: tol      = 0.01   ! converged when |N_new - N_old| <= tol

 integer :: iter, n_old, n_new
 real    :: C_rho, dr_cgs, r_out_cgs, exponent, M_shell_cgs

 C_rho    = rho_inner_cgs * r_min_cgs**rho_power
 exponent = 3.0 - rho_power

 ! Start from a reasonable guess: treat the shell as infinitesimally thin
 ! so M ~ 4*pi*r^2 * rho(r) * dr with dr ~ r * spacing(1000)
 n_old = 1000
 do iter = 1, max_iter
    dr_cgs    = wss_in * r_min_cgs * get_fibonacci_spacing(n_old)
    r_out_cgs = r_min_cgs + dr_cgs

    ! Exact mass of the shell from the analytic power-law integral
    M_shell_cgs = 4.0*pi * C_rho * (r_out_cgs**exponent - r_min_cgs**exponent) / exponent

    n_new = max(1, nint(M_shell_cgs / m_particle_cgs))

    if (abs(real(n_new - n_old)) <= tol) exit

    n_old = (n_old + n_new) / 2
 enddo

 n_first = n_new

 if (verbose == 1) then
    print *, ''
    print *, ' derive_n_particles_first: converged in ', iter, ' iterations'
    print *, ' Target particle mass (cgs)    :', m_particle_cgs
    print *, ' First-shell particle count    :', n_first
    print *, ' First-shell dr / r_min        :', dr_cgs / r_min_cgs
    print *, ''
 endif

end subroutine derive_n_particles_first

!----------------------------------------------------------------
!+
!  Initialize everything
!+
!----------------------------------------------------------------
subroutine init_inject(ierr)
 use io,            only:fatal
 use physcon,       only:pi,days,au,solarm,km,years
 use eos,           only:gmw,gamma
 use units,         only:utime,umass,udist,unit_velocity,unit_luminosity
 use part,          only:xyzmh_ptmass,massoftype,igas,iboundary,nptmass,iTeff,iReff,iLum,npartoftype
 use injectutils,   only:get_parts_per_sphere, get_fibonacci_spacing, find_optimal_rotation
 use wind_pulsating,only:setup_star,calc_stellar_profile,region_mass,interp_stellar_profile
 use dust_formation,only:calc_kappa_max
 use timestep,      only:dtmax

 integer, intent(out) :: ierr
 real    :: Mstar_cgs, Rstar_cgs, Tstar, Lstar_cgs
 real    :: current_radius, dr
 integer :: shell_index, max_shells, n_shell
 integer :: n_particles_first              ! derived locally — no longer a module parameter
 integer :: expected_measurements, i
 integer, parameter  :: max_shells_tmp = 2000
 integer, parameter  :: max_iter_dr    = 100
 real,    parameter  :: tol_dr         = 1.0e-6
 real    :: tmp_dr(max_shells_tmp), tmp_r(max_shells_tmp)
 integer :: tmp_n(max_shells_tmp), n_tot
 logical :: file_exists
 integer :: iunit
 real    :: r_min_cgs, m_particle_cgs

 ierr = 0

 if (nptmass < 1) call fatal(label,'need at least one sink particle for central star')

 Mstar     = xyzmh_ptmass(4, wind_emitting_sink)
 Rstar     = xyzmh_ptmass(iReff, wind_emitting_sink)
 Rstar_cgs = Rstar * au
 Mstar_cgs = Mstar * solarm
 Tstar     = xyzmh_ptmass(iTeff, wind_emitting_sink)
 Lstar_cgs = xyzmh_ptmass(iLum, wind_emitting_sink) * unit_luminosity

 call calc_kappa_max(Mstar_cgs, Lstar_cgs)

 inquire(file='mass_loss_rate.dat', exist=file_exists)

 if ( npartoftype(igas) < 100 .and. file_exists .and. use_file_mdot == 0) then
       print *, 'Existing mass loss data file found, but this is a fresh start, so delete'
       open(newunit=iunit, file='mass_loss_rate.dat', status='old', iostat=ierr)
       close(iunit, status='delete')
       file_exists = .false.
 endif

 if (iwind == 2 .and. .not. file_exists) call calculate_period(Mstar, Rstar, pulsation_period_days)

 pulsation_period = pulsation_period_days * (days / utime)
 omega_pulsation  = 2.0*pi / pulsation_period
 piston_velocity  = piston_velocity_km_s * (km / unit_velocity)
 deltaR_osc       = pulsation_period * piston_velocity / (2.0*pi)

 r_min = r_min_on_rstar * Rstar + deltaR_osc * sin(phi0)

 if (save_period == 1) then
    dtmax = 1. / (dumps_p_period) * pulsation_period
    print *, 'dtmax: ', dtmax
 endif

 if (r_min <= 0.) call fatal(label,'r_min must be > 0')

 mass_of_gas_particle = massoftype(igas)          ! code units (Msun when G=1, dist=au)
 m_particle_cgs       = mass_of_gas_particle * umass
 r_min_cgs            = r_min * udist

 call derive_n_particles_first(r_min_cgs, rho_inner, rho_power_in, &
                                wss, m_particle_cgs, n_particles_first)

 current_radius = r_min
 shell_index    = 0
 n_tot          = 0

 do
    shell_index = shell_index + 1
    if (shell_index > max_shells_tmp) &
       call fatal(label,'max_shells_tmp exceeded; increase max_shells_tmp')
    if (shell_index == 1) then
       n_shell = n_particles_first
    else
       n_shell = max(1, nint(real(tmp_n(shell_index-1)) * (current_radius / tmp_r(shell_index-1))**( 2.*(3.-rho_power_in) / 3.)))
    endif

    dr = wss * current_radius * get_fibonacci_spacing(n_shell)

    if (shell_index > 2 .and. ( real(n_tot + n_shell) / real(n_tot) - 1) < real(min_mass_fraction)) then
       shell_index = shell_index - 1
       r_max = current_radius - dr
       exit
    endif

    tmp_dr(shell_index) = dr
    tmp_r(shell_index)  = current_radius
    tmp_n(shell_index)  = n_shell
    current_radius      = current_radius + dr
    n_tot               = n_tot + n_shell
 enddo

 call setup_star(Mstar * umass, Tstar, r_max * au, r_min * au, gmw, gamma, rho_inner, rho_power_in)

 call calc_stellar_profile(n_profile_points)

 reinject_period        = pulsation_period / real(n_inject_period)
 measurement_interval   = reinject_period
 mass_loss_start_time   = mass_loss_start * pulsation_period
 mass_loss_end_time     = mass_loss_end   * pulsation_period
 time_next_measurement  = mass_loss_start_time
 n_measurements         = 0

 expected_measurements = ceiling((mass_loss_end_time - mass_loss_start_time) / measurement_interval) + 1
 allocate(mass_loss_rates(expected_measurements))
 mass_loss_rates = 0.0

 if (verbose == 1) then
   print *, ''
   print *, 'Calculated reinject period:', reinject_period
   print *, 'Measurement period:', measurement_interval
   print *, 'Mass loss measurement start time:', mass_loss_start_time
   print *, 'Mass loss measurement end time  :', mass_loss_end_time
   print *, 'Rmax                            :', r_max
   print *, 'Expected number of measurements :', expected_measurements
   print *, ''
 endif

 n_shells_total = shell_index
 n_shells_bnd   = min(iboundary_spheres, n_shells_total)

 allocate(npart_per_boundary_shell(n_shells_bnd))
 allocate(delta_r_boundary(n_shells_bnd))
 allocate(shell_radii_bnd(n_shells_bnd))
 do i = 1, n_shells_bnd
    npart_per_boundary_shell(i) = tmp_n(i)
    delta_r_boundary(i)         = tmp_dr(i)
    shell_radii_bnd(i)          = tmp_r(i)
 enddo

 allocate(npart_per_shell(n_shells_total - n_shells_bnd))
 allocate(delta_r_gas(n_shells_total - n_shells_bnd))
 allocate(shell_radii_gas(n_shells_total - n_shells_bnd))
 do i = 1, n_shells_total - n_shells_bnd
    npart_per_shell(i)  = tmp_n(n_shells_bnd + i)
    delta_r_gas(i)      = tmp_dr(n_shells_bnd + i)
    shell_radii_gas(i)  = tmp_r(n_shells_bnd + i)
 enddo
 n_shells_total = n_shells_total - n_shells_bnd

 if (allocated(delta_r_radial)) deallocate(delta_r_radial)
 allocate(delta_r_radial(n_shells_bnd + n_shells_total))
 if (n_shells_bnd > 0) delta_r_radial(1:n_shells_bnd) = delta_r_boundary
 delta_r_radial(n_shells_bnd+1 : n_shells_bnd+n_shells_total) = delta_r_gas

 massoftype(igas)      = mass_of_gas_particle
 massoftype(iboundary) = mass_of_gas_particle

 Mtotal = Mstar + sum(tmp_n(1:n_shells_total)) * mass_of_gas_particle

 ! If use_file_mdot=1, read mass_loss_rate.dat now and mark measurement as done,
 ! so the simulation goes straight to reinjection without any measurement phase.
 ! A hard error is raised if the file is absent, since the user explicitly asked for it.
 if (use_file_mdot == 1) then
    if (.not. file_exists) &
       call fatal(label,'use_file_mdot=1 but mass_loss_rate.dat not found')
    call read_mass_loss_data()
    mass_loss_rate_calculated = .true.
    call find_optimal_rotation(particles_to_inject)
    if (verbose == 1) then
       print *, ''
       print *, 'use_file_mdot=1: skipping measurement phase.'
       print *, 'particles_to_inject read from file:', particles_to_inject
       print *, ''
    endif
 elseif (file_exists) then
    call read_mass_loss_data()
    if (mass_loss_rate_calculated) call find_optimal_rotation(particles_to_inject)
 endif

 if (verbose == 1) then
    print *, ''
    print *, ' rho_power                        :', rho_power_in
    print *, ' rho_inner (cgs)                  :', rho_inner
    print *, ' Particle mass from setup (Msun)  :', massoftype(igas)
    print *, ' n_particles_first (derived)       :', n_particles_first
    print *, ' Atmosphere [r_min, r_max] Rstar  :', r_min, r_max
    print *, ' M_atmos / M_total                :', region_mass(r_min, r_max) / Mtotal
    print *, ' M_atmos (Msun)                   :', region_mass(r_min, r_max)
    print *, ' Boundary shells                  :', n_shells_bnd
    print *, ' Gas shells                       :', n_shells_total
    print *, ' Total boundary particles         :', sum(npart_per_boundary_shell)
    print *, ' Total gas particles              :', sum(npart_per_shell)
    print *, ' Innermost boundary N_per_shell   :', npart_per_boundary_shell(1)
    print *, ' Outermost boundary N_per_shell   :', npart_per_boundary_shell(n_shells_bnd)
    print *, ' Innermost gas      N_per_shell   :', npart_per_shell(1)
    print *, ' Outermost gas      N_per_shell   :', npart_per_shell(n_shells_total)
    print *, ' Particle mass (Msun)             :', mass_of_gas_particle
    print *, ''
 endif

 if (verbose == 1) then
    do i = 1, n_shells_bnd
       print *, 'Boundary shell ', i, ': r=', shell_radii_bnd(i)/Rstar, ' Rstar, dr=', delta_r_boundary(i)/Rstar, &
                ' Rstar, N_particles=', npart_per_boundary_shell(i)
    enddo
    do i = 1, n_shells_total
       print *, 'Gas shell      ', i, ': r=', shell_radii_gas(i)/Rstar, ' Rstar, dr=', delta_r_gas(i)/Rstar, &
                ' Rstar, N_particles=', npart_per_shell(i)
    enddo
 endif

end subroutine init_inject

!----------------------------------------------------------------
!+
!  The actual function that is called by phantom
!+
!----------------------------------------------------------------
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
    print *, ''
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

 if (reinject_enabled == 1 .and. .not. mass_loss_rate_calculated) then
    call take_periodic_mass_measurements(time,xyzh,vxyzu,npart,xyzmh_ptmass,vxyz_ptmass,npartoftype)
 endif

 if (reinject_enabled == 1 .and. mass_loss_rate_calculated) then
    if ((time - time_last_reinject) >= reinject_period .and. time >= mass_loss_start_time) then
       call perform_reinjection(time,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npartoftype)
       time_last_reinject = time
       reinjection_needed = .false.
    endif
 endif

 call apply_pulsation(time,xyzh,vxyzu,npart,xyzmh_ptmass,vxyz_ptmass)

end subroutine inject_particles

!----------------------------------------------------------------
!+
!  Checks how much mass the star has lost, and calculates the mass loss rate
!+
!----------------------------------------------------------------
subroutine take_periodic_mass_measurements(time,xyzh,vxyzu,npart,xyzmh_ptmass,vxyz_ptmass,npartoftype)
 use part,        only:igas,iboundary,iamtype
 use injectutils, only:find_optimal_rotation

 real,    intent(in) :: time
 real,    intent(in) :: xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(in) :: npart
 integer, intent(in) :: npartoftype(:)

 real    :: rate_this_interval, sum_rates
 real    :: x_agb(3), v_agb(3)
 real    :: dx1, dy1, dz1, r1
 real    :: M1, e_therm, e_kin, e_pot
 integer :: i, n_escaping, newly_unbound
 logical :: unbound

 M1    = xyzmh_ptmass(4, wind_emitting_sink)
 x_agb = xyzmh_ptmass(1:3, wind_emitting_sink)
 v_agb = vxyz_ptmass(1:3,  wind_emitting_sink)

 if (.not. measurement_active .and. time >= mass_loss_start_time) then
    measurement_active    = .true.
    time_next_measurement = time + measurement_interval

    n_escaping_prev = 0
    do i = 1, npart
       dx1   = xyzh(1,i) - x_agb(1)
       dy1   = xyzh(2,i) - x_agb(2)
       dz1   = xyzh(3,i) - x_agb(3)
       r1    = sqrt(dx1**2 + dy1**2 + dz1**2)
       if (r1 <= 0.) cycle

       e_pot = - xyzmh_ptmass(4, 1) / r1
       e_kin = 0.5 * ( (vxyzu(1,i) - vxyz_ptmass(1, 1))**2 &
                     + (vxyzu(2,i) - vxyz_ptmass(2, 1))**2 &
                     + (vxyzu(3,i) - vxyz_ptmass(3, 1))**2 )
       e_therm = vxyzu(4,i)
       unbound   = .false.
       if (e_kin + e_therm + e_pot > 0.) unbound = .true.
       if (unbound) n_escaping_prev = n_escaping_prev + 1
    enddo

    if (verbose == 1) print *, ' Baseline unbound particles at measurement start:', n_escaping_prev
 endif

 if (measurement_active .and. time >= time_next_measurement .and. time < mass_loss_end_time) then

    n_escaping = 0
    do i = 1, npart
       dx1   = xyzh(1,i) - x_agb(1)
       dy1   = xyzh(2,i) - x_agb(2)
       dz1   = xyzh(3,i) - x_agb(3)
       r1    = sqrt(dx1**2 + dy1**2 + dz1**2)
       if (r1 <= 0.) cycle

       e_pot = - xyzmh_ptmass(4, 1) / r1
       e_kin = 0.5 * ( (vxyzu(1,i) - vxyz_ptmass(1, 1))**2 &
                     + (vxyzu(2,i) - vxyz_ptmass(2, 1))**2 &
                     + (vxyzu(3,i) - vxyz_ptmass(3, 1))**2 )
       e_therm = vxyzu(4,i)
       unbound   = .false.
       if (e_kin + e_therm + e_pot > 0.) unbound = .true.
       if (unbound) n_escaping = n_escaping + 1
    enddo

    newly_unbound           = max(0, n_escaping - n_escaping_prev)
    rate_this_interval      = real(newly_unbound) * mass_of_gas_particle / measurement_interval
    n_escaping_prev         = n_escaping
    n_measurements          = n_measurements + 1
    mass_loss_rates(n_measurements) = rate_this_interval
    time_next_measurement   = time + measurement_interval

    if (verbose == 1) then
       print *, ' Unbound+outflowing particles     :', n_escaping
       print *, ' Newly unbound since last snapshot:', newly_unbound
       print *, ' Rate this interval               :', rate_this_interval
    endif
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
       call find_optimal_rotation(particles_to_inject)
    endif
 endif

end subroutine take_periodic_mass_measurements

!----------------------------------------------------------------
!+
!  Inject particles throughout the simulation
!+
!----------------------------------------------------------------
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

 call inject_fibonacci_sphere(n_shells_total + n_reinjections, npart + 1, particles_to_inject, r_inject, r_dot, u, rho, &
                               npart, npartoftype, xyzh, vxyzu, igas, x0, v0)

 mass_injected = real(npart - old_npart) * mass_of_gas_particle
 xyzmh_ptmass(4, wind_emitting_sink) = xyzmh_ptmass(4, wind_emitting_sink) - mass_injected

 if (verbose == 1) then
    print *, ''
    print *, ' Particles injected         :', (npart - old_npart)
    print *, ' Injection radius           :', r_inject
    print *, ' New total particles        :', npart
    print *, 'Reinjection complete.'
    print *, ''
 endif

end subroutine perform_reinjection

!----------------------------------------------------------------
!+
!  Build the initial atmospheric setup (i.e. build the shells)
!+
!----------------------------------------------------------------
subroutine setup_initial_atmosphere(xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npartoftype)
 use part,           only:igas,iboundary,iamtype
 use injectutils,    only:inject_fibonacci_sphere
 use wind_pulsating, only:interp_stellar_profile
 use physcon,        only:pi,km,au

 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real,    intent(in)    :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)

 integer :: i, first_particle, nboundary
 real    :: r, rho, u, T, P, x0(3), v0(3), v_radial

 x0 = xyzmh_ptmass(1:3, wind_emitting_sink)
 v0 = vxyz_ptmass(1:3,  wind_emitting_sink)

 npart = 0

 do i = 1, n_shells_bnd
    r     = shell_radii_bnd(i)
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

    do i = 1, nboundary
       r_boundary_equilibrium(i) = sqrt( (xyzh(1,i)-x0(1))**2 + &
                                         (xyzh(2,i)-x0(2))**2 + &
                                         (xyzh(3,i)-x0(3))**2 ) &
                                      - deltaR_osc * sin(phi0)
    enddo
 endif

end subroutine setup_initial_atmosphere

!----------------------------------------------------------------
!+
!  Reconstructs boundary particle info after resuming from a dump
!+
!----------------------------------------------------------------
subroutine reconstruct_boundary_info(time,xyzh,npart,xyzmh_ptmass)
 use part,   only:iboundary,iamtype,npartoftype
 use physcon,only:pi

 real,    intent(in) :: time
 real,    intent(inout) :: xyzh(:,:),xyzmh_ptmass(:,:)
 integer, intent(in) :: npart
 integer :: i
 real    :: x0(3), r_current, phase

 x0    = xyzmh_ptmass(1:3, wind_emitting_sink)
 phase = omega_pulsation * time + phi0

 n_boundary_particles = npartoftype(3)

 if (n_boundary_particles > 0) then
    allocate(r_boundary_equilibrium(n_boundary_particles))
    allocate(boundary_particle_ids(n_boundary_particles))

    do i = 1, n_boundary_particles
       r_current = sqrt((xyzh(1,i)-x0(1))**2 + &
                        (xyzh(2,i)-x0(2))**2 + &
                        (xyzh(3,i)-x0(3))**2)
       r_boundary_equilibrium(i) = r_current - deltaR_osc * sin(phase)
    enddo

    print *, 'Reconstructed boundary particle info:'
    print *, 'Boundary particles: ', n_boundary_particles
 endif

end subroutine reconstruct_boundary_info

!----------------------------------------------------------------
!+
!  Applies the pulsation to the boundary layers
!+
!----------------------------------------------------------------
subroutine apply_pulsation(time,xyzh,vxyzu,npart,xyzmh_ptmass,vxyz_ptmass)
 use physcon,        only:pi,solarl
 use wind_pulsating, only:interp_stellar_profile
 use part,           only:iTeff,iLum,iReff

 real,    intent(in)    :: time
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(in)    :: npart

 integer :: i
 real    :: r_eq, r_new, r_current, phase, deltaR_osc
 real    :: x_hat(3), r_dot, x0(3), v0(3)
 real    :: x, y, z
 real    :: Reff, Teff, Lum

 if (.not. allocated(boundary_particle_ids)) return
 if (n_boundary_particles == 0) return

 x0 = xyzmh_ptmass(1:3, wind_emitting_sink)
 v0 = vxyz_ptmass(1:3,  wind_emitting_sink)

 phase      = omega_pulsation * time + phi0
 deltaR_osc = pulsation_period * piston_velocity / (2.0 * pi)
 r_dot      = piston_velocity * cos(phase)

 do i = 1, n_boundary_particles
    r_eq  = r_boundary_equilibrium(i)
    r_new = r_eq + deltaR_osc * sin(phase)

    x = xyzh(1,i) - x0(1)
    y = xyzh(2,i) - x0(2)
    z = xyzh(3,i) - x0(3)
    r_current = sqrt(x**2 + y**2 + z**2)

    ! Radial unit vector from current position
    x_hat(1) = x / r_current
    x_hat(2) = y / r_current
    x_hat(3) = z / r_current

    xyzh(1,i) = r_new + x0(1)
    xyzh(2,i) = r_new + x0(2)
    xyzh(3,i) = r_new + x0(3)

    vxyzu(1,i) = r_dot + v0(1)
    vxyzu(2,i) = r_dot + v0(2)
    vxyzu(3,i) = r_dot + v0(3)

    if (update_L == 1) then
       Reff = xyzmh_ptmass(iReff,1) + deltaR_osc * sin(phase)
       Teff = xyzmh_ptmass(iTeff,1)
       Lum  = xyzmh_ptmass(iLum,1)
       call get_lum(Lum, Teff, Reff)
       xyzmh_ptmass(iLum,1) = Lum
    endif

 enddo

end subroutine apply_pulsation

!----------------------------------------------------------------
!+
!  Placeholder function
!+
!----------------------------------------------------------------
subroutine update_injected_par

end subroutine update_injected_par

!----------------------------------------------------------------
!+
!  Get luminosity
!+
!----------------------------------------------------------------
subroutine get_lum(Lum,Teff,Reff)
 use physcon, only:au,steboltz,solarl,pi
 use units,   only:unit_luminosity
 real, intent(inout) :: Lum
 real, intent(in)    :: Reff, Teff
 real :: lum_lsun

 lum_lsun = 4.*pi*steboltz*Teff**4*(Reff*au)**2/solarl
 Lum  = lum_lsun*(solarl/unit_luminosity)

end subroutine get_lum

!----------------------------------------------------------------
!+
!  Write mass-loss information to file for resume from dump
!+
!----------------------------------------------------------------
subroutine write_mass_loss_data()
 use io, only:iprint
 use physcon, only:solarm,years
 use units,   only:umass,utime
 integer :: iunit, ierr, i

 if (.not. mass_loss_rate_calculated) return

 open(newunit=iunit, file='mass_loss_rate.dat', status='replace', iostat=ierr)
 if (ierr /= 0) then
    write(iprint,*) 'Could not write mass_loss_rate.dat'
    return
 endif

 write(iunit,*) '# Mass-loss rate data for restart'
 write(iunit,*) r_max
 write(iunit,*) r_max
 write(iunit,*) mass_loss_rate_calculated
 write(iunit,*) mean_mass_loss_rate / (solarm / umass) / (utime / years)
 write(iunit,*) Mtotal
 write(iunit,*) particles_to_inject
 write(iunit,*) n_measurements
 write(iunit,*) mass_of_gas_particle
 do i = 1, n_measurements
    write(iunit,*) mass_loss_rates(i)
 enddo
 close(iunit)

 print *, ''
 write(iprint,*) 'Mass-loss rate data written to mass_loss_rate.dat'
 print *, ' '

end subroutine write_mass_loss_data

!----------------------------------------------------------------
!+
!  Read mass-loss information from file after resuming from dump
!+
!----------------------------------------------------------------
subroutine read_mass_loss_data()
 use io, only:iprint
 integer :: iunit, ierr, i
 logical :: file_exists

 inquire(file='mass_loss_rate.dat', exist=file_exists)
 if (.not. file_exists) return

 open(newunit=iunit, file='mass_loss_rate.dat', status='old', iostat=ierr)
 if (ierr /= 0) return

 read(iunit,*)
 read(iunit,*, iostat=ierr) r_max
 read(iunit,*, iostat=ierr) r_max
 read(iunit,*, iostat=ierr) mass_loss_rate_calculated
 if (ierr /= 0) then; close(iunit); return; endif
 read(iunit,*, iostat=ierr) mean_mass_loss_rate
 read(iunit,*, iostat=ierr) Mtotal
 read(iunit,*, iostat=ierr) particles_to_inject
 read(iunit,*, iostat=ierr) n_measurements
 read(iunit,*, iostat=ierr) mass_of_gas_particle

 if (n_measurements > 0) then
    if (.not. allocated(mass_loss_rates)) allocate(mass_loss_rates(n_measurements))
    do i = 1, n_measurements
       read(iunit,*, iostat=ierr) mass_loss_rates(i)
       if (ierr /= 0) exit
    enddo
 endif
 close(iunit)

 if (verbose == 1) then
    write(iprint,*) 'Mass-loss rate data read from mass_loss_rate.dat'
    write(iprint,*) ' Mean mass-loss rate          :', mean_mass_loss_rate
    write(iprint,*) ' Gas particle mass            :', mass_of_gas_particle
    write(iprint,*) ' Boundary particle mass       :', mass_of_gas_particle
    write(iprint,*) ' Particles to inject          :', particles_to_inject
 endif

end subroutine read_mass_loss_data

!----------------------------------------------------------------
!+
!  Use mass-period relation to estimate the pulsation period
!+
!----------------------------------------------------------------
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

!----------------------------------------------------------------
!+
!  Write options to .in file
!+
!----------------------------------------------------------------
subroutine write_options_inject(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(n_profile_points,      'n_profile_points',   'number of points in stellar profile',iunit)
 call write_inopt(iboundary_spheres,     'iboundary_spheres',  'number of boundary spheres (piston layers)',iunit)
 call write_inopt(min_mass_fraction,     'min_mass_fraction',  'minimum mass fraction per shell',iunit)
 call write_inopt(rho_power_in,          'rho_power',          'density profile exponent: rho ~ r^(-rho_power)',iunit)
 call write_inopt(r_min_on_rstar,        'r_min_on_rstar',     'gas atmosphere inner radius as fraction of R_star',iunit)
 call write_inopt(rho_inner,             'rho_inner',          'inner boundary density at r_min (cgs)',iunit)
 call write_inopt(iwind,                 'iwind',              'wind type: 1=prescribed, 2=period from mass-radius relation',iunit)
 call write_inopt(pulsation_period_days, 'pulsation_period',   'pulsation period (days)',iunit)
 call write_inopt(piston_velocity_km_s,  'piston_velocity',    'piston velocity amplitude (km/s)',iunit)
 call write_inopt(phi0,                  'phi0',               'initial phase offset (radians)',iunit)
 call write_inopt(wss,                   'wss',                'radial/tangential spacing ratio',iunit)
 call write_inopt(save_period,           'save_period',        'wether to save dumps as fraction of period (0=off, 1=on)',iunit)
 call write_inopt(dumps_p_period,        'dumps_p_period',     'number of dumps per period (if save_period = 1)',iunit)
 call write_inopt(reinject_enabled,      'reinject_enabled',   'enable dynamic reinjection (0=off, 1=on)',iunit)
 call write_inopt(n_inject_period,       'n_inject_period',    'period between reinjections (periods)',iunit)
 call write_inopt(mass_loss_start,       'mass_loss_start',    'start time for mass-loss calculation (periods)',iunit)
 call write_inopt(mass_loss_end,         'mass_loss_end',      'end time for mass-loss calculation (periods)',iunit)
 call write_inopt(use_file_mdot,         'use_file_mdot',      'skip measurement phase (0=off, 1=on)',iunit)
 call write_inopt(update_L,              'update_L',           'update luminosity with pulsation (0=off, 1=on)',iunit)
 call write_inopt(verbose,               'verbose',            'enable verbose output (0=off, 1=on)',iunit)

end subroutine write_options_inject

!----------------------------------------------------------------
!+
!  Read options from .in file
!+
!----------------------------------------------------------------
subroutine read_options_inject(name,valstring,imatch,igotall,ierr)
 use io, only:fatal
 character(len=*), intent(in)  :: name, valstring
 logical,          intent(out) :: imatch, igotall
 integer,          intent(out) :: ierr

 integer, save      :: ngot = 0
 integer, parameter :: noptions = 20
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
 case('min_mass_fraction')
    read(valstring,*,iostat=ierr) min_mass_fraction
    ngot = ngot + 1
    if (min_mass_fraction < 0) call fatal(label,'min_mass_fraction must be >= 0')
 case('rho_power')
    read(valstring,*,iostat=ierr) rho_power_in
    ngot = ngot + 1
    if (rho_power_in <= 0.) call fatal(label,'rho_power must be > 0')
 case('r_min_on_rstar')
    read(valstring,*,iostat=ierr) r_min_on_rstar
    ngot = ngot + 1
    if (r_min_on_rstar <= 0. .or. r_min_on_rstar >= 2.0) &
       call fatal(label,'r_min_on_rstar must be in (0,2)')
 case('rho_inner')
    read(valstring,*,iostat=ierr) rho_inner
    ngot = ngot + 1
    if (rho_inner <= 0.) call fatal(label,'rho_inner must be > 0')
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
 case('phi0')
    read(valstring,*,iostat=ierr) phi0
    ngot = ngot + 1
    if (phi0 < -3.1415926536d0 .or. phi0 > 3.1415926536d0) &
       call fatal(label,'phi0 must be in (-pi,pi)')
 case('wss')
    read(valstring,*,iostat=ierr) wss
    ngot = ngot + 1
    if (wss <= 0. .or. wss > 10.0) call fatal(label,'wss must be in (0,10]')
 case('save_period')
    read(valstring,*,iostat=ierr) save_period
    ngot = ngot + 1
    if (save_period /= 0 .and. save_period /= 1) call fatal(label,'save_period must be 0 or 1')
 case('dumps_p_period')
    read(valstring,*,iostat=ierr) dumps_p_period
    ngot = ngot + 1
    if (dumps_p_period < 0) call fatal(label,'dumps_p_period must be > 0')
 case('reinject_enabled')
    read(valstring,*,iostat=ierr) reinject_enabled
    ngot = ngot + 1
    if (reinject_enabled /= 0 .and. reinject_enabled /= 1) call fatal(label,'reinject_enabled must be 0 or 1')
 case('n_inject_period')
    read(valstring,*,iostat=ierr) n_inject_period
    ngot = ngot + 1
    if (n_inject_period <= 0) call fatal(label,'n_inject_period must be > 0')
 case('mass_loss_start')
    read(valstring,*,iostat=ierr) mass_loss_start
    ngot = ngot + 1
    if (mass_loss_start < 0.) call fatal(label,'mass_loss_start must be >= 0')
 case('mass_loss_end')
    read(valstring,*,iostat=ierr) mass_loss_end
    ngot = ngot + 1
    if (mass_loss_end <= 0.) call fatal(label,'mass_loss_end must be > 0')
 case('use_file_mdot')
    read(valstring,*,iostat=ierr) use_file_mdot
    ngot = ngot + 1
    if (use_file_mdot /= 0 .and. use_file_mdot /= 1) call fatal(label,'use_file_mdot must be 0 or 1')
 case('update_L')
    read(valstring,*,iostat=ierr) update_L
    ngot = ngot + 1
    if (update_L /= 0 .and. update_L /= 1) call fatal(label,'update_L must be 0 or 1')
 case('verbose')
    read(valstring,*,iostat=ierr) verbose
    ngot = ngot + 1
    if (verbose /= 0 .and. verbose /= 1) call fatal(label,'verbose must be 0 or 1')
 case default
    imatch = .false.
 end select

 igotall = (ngot >= noptions)

end subroutine read_options_inject

end module inject