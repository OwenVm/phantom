!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module inject
!
! Handles initial setup of stellar atmosphere with pulsating boundary layers
! Modified to calculate mass-loss rate at regular intervals and use mean for reinjection
!
! :References: None
!
! :Owner: Owen Vermeulen
!
! :Runtime parameters:
!   - iboundary_spheres  : *number of boundary spheres (integer)*
!   - n_profile_points   : *number of points in stellar profile calculation (integer)*
!   - n_particles        : *total number of particles (integer, not used if n_shells > 0)*
!   - n_shells          : *total number of shells (integer, if <0 determined automatically from n_particles)*
!   - r_min_on_rstar     : *inner radius as fraction of R_star*
!   - r_max_on_rstar     : *outer radius as fraction of R_star*
!   - pulsation_period   : *pulsation period (days)*
!   - pulsation_amplitude: *fractional pulsation amplitude*
!   - piston_velocity    : *piston velocity amplitude (km/s)*
!   - atmos_mass_fraction: *atmospheric mass as fraction of total stellar mass*
!   - surface_pressure   : *surface pressure (cgs)*
!   - iwind              : *wind type: 1=prescribed, 2=period from mass-radius relation*
!   - pulsation_timestep : *pulsation timestep as fraction of pulsation period*
!   - phi0               : *initial phase offset (radians) (best taken to be -pi/2 to start at minimum radius)*
!   - wss                : *fraction of tangential and radial distance between particles in initial atmosphere setup*
!   - reinject_enabled   : *enable dynamic reinjection (logical)*
!   - reinject_period_days: *period between reinjections in days (for continuous mode)*
!   - mass_loss_start: *start time for mass-loss calculation in years*
!   - mass_loss_end: *end time for mass-loss calculation in years*
!   - check_radius_au: *radius within which to count mass (AU)*
!   - meas_int_days: *interval for mass measurements in days*
!
! :Dependencies: dim, eos, icosahedron, infile_utils, injectutils, io,
!   part, partinject, physcon, units, set_star
!
 use io,            only:fatal
 implicit none
 character(len=*), parameter, public :: inject_type = 'atmosphere'

 public :: init_inject,inject_particles,write_options_inject,read_options_inject,&
           set_default_options_inject,update_injected_par
 private

!--runtime settings for this module
!
! Read from input file
 integer :: iboundary_spheres = 5 ! Number of boundary spheres 
 integer :: n_profile_points = 10000 ! Number of points in stellar profile calculation
 integer :: n_particles = 500000 ! Total number of particles (not used if n_shells > 0)
 integer :: n_shells = 25 ! Total number of shells (if <0 determined automatically from n_particles)
 real    :: r_min_on_rstar = 0.9 ! Inner radius (R_eq, not R_min) as fraction of Rstar
 real    :: r_max_on_rstar = 1.4 ! Outer radius as fraction of Rstar
 real    :: dtpulsation = huge(0.)
 real    :: pulsation_period_days = 300.0  ! Pulsation period in days
 real    :: piston_velocity_km_s = 4.0     ! Piston velocity (in km/s)
 real    :: time_puls = -1.0 ! Time for the piston to accelerate from 0 to max velocity (in periods)
 real    :: atmos_mass_fraction = 5e-5  ! Atmosphere mass as fraction of total mass
 real    :: surface_pressure = 0.001  ! Surface pressure in cgs units
 integer :: iwind = 1  ! Wind type: 1=prescribed, 2=period from mass-radius relation
 real    :: pulsation_timestep = 0.02
 real    :: phi0 = -3.1415926536d0/2.0  ! Initial phase offset (-pi/2 for starting at minimal radius)
 real    :: wss = 1.0 ! Fraction of the tangential and radial distance between particles in the initial setup
 logical :: var_boundary = .false.

 ! Reinjection parameters
 logical :: reinject_enabled = .true.
 real    :: reinject_period_days = 10.0  ! Period between reinjections (days)
 real    :: mass_loss_start = 1.0  ! Start time for mass-loss calculation (years)
 real    :: mass_loss_end = 2.0    ! End time for mass-loss calculation (years)
 real    :: check_radius_au = 3.0  ! Radius within which to count mass (AU)
 real    :: meas_int_days = 10.0  ! Interval for mass measurements (days)


! global variables
 integer, parameter :: wind_emitting_sink = 1
 integer, parameter :: max_measurements = 10000  ! Maximum number of measurements to store
 real :: omega_pulsation, deltaR_osc, pulsation_period, piston_velocity
 real :: Rstar, r_min, mass_of_particles
 real :: Mtotal, Matmos, Msink  ! Total, atmosphere, and sink masses
 real, allocatable :: delta_r_radial(:)
 integer :: particles_per_sphere
 logical :: atmosphere_setup_complete = .false.
 real, allocatable :: shell_radii(:)  ! Store radii for each shell
 integer :: n_shells_total
 
 ! Store boundary particle information
 real, allocatable    :: r_boundary_equilibrium(:)
 integer, allocatable :: boundary_particle_ids(:)
 integer              :: n_boundary_particles
 integer              :: active_boundary_spheres  ! Current number of active boundary spheres
 logical              :: reinjection_needed = .false.  ! Flag to track if reinjection should happen
 
 ! Continuous reinjection tracking
 real    :: time_last_reinject = 0.0  ! Time of last reinjection
 real    :: reinject_period  ! Period in code units
 integer :: n_reinjections = 0  ! Number of reinjections performed so far
 
 ! Mass-loss rate tracking with periodic measurements
 real    :: mass_loss_start_time  ! Start time in code units
 real    :: mass_loss_end_time    ! End time in code units
 real    :: mass_loss_check_radius  ! Radius in code units
 real    :: measurement_interval  ! Measurement interval in code units
 real    :: time_next_measurement  ! Time for next mass measurement
 real    :: mass_previous_measurement  ! Mass at previous measurement
 real, allocatable :: mass_loss_rates(:)  ! Array to store individual mass-loss rates
 integer :: n_measurements  ! Number of measurements taken
 real    :: mean_mass_loss_rate = 0.0  ! Mean mass-loss rate in code units (Msun/time)
 logical :: mass_loss_rate_calculated = .false.
 logical :: measurement_active = .false.
 integer :: particles_to_inject = 0  ! Number of particles to inject based on mass-loss rate

 character(len=*), parameter :: label = 'inject_atmosphere'

contains

!-----------------------------------------------------------------------
!+
!  Set default options
!+
!-----------------------------------------------------------------------
subroutine set_default_options_inject(flag)
 integer, optional, intent(in) :: flag

 iboundary_spheres = 5
 n_profile_points = 10000
 n_particles = 500000
 n_shells = 25
 r_min_on_rstar = 0.9
 r_max_on_rstar = 1.4
 dtpulsation = huge(0.)
 atmos_mass_fraction = 5e-5
 surface_pressure = 0.001
 iwind = 1
 pulsation_period_days = 300.0
 piston_velocity_km_s = 4.0
 time_puls = -1.0
 pulsation_timestep = 0.02
 phi0 = -3.1415926536d0/2.0
 wss = 1.0
 var_boundary = .false.
 reinject_enabled = .true.
 reinject_period_days = 10.0
 mass_loss_start = 1.0
 mass_loss_end = 3.0
 check_radius_au = 3.0
 meas_int_days = 10.0

end subroutine set_default_options_inject

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
 use wind_pulsating,only:setup_star,calc_stellar_profile
 use dust_formation,only:calc_kappa_max

 integer, intent(out) :: ierr
 real :: Mstar_cgs, Rstar_cgs, Tstar, Lstar_cgs, delta_r_tangential, current_radius
 integer :: shell_index, particles_per_shell, distributed_particles, max_shells, temp_particles
 integer :: expected_measurements
 logical :: converged, file_exists

 ierr = 0

 if (nptmass < 1) then
    call fatal(label,'need at least one sink particle for central star')
 endif

 ! Get stellar properties from sink particle
 Mtotal    = xyzmh_ptmass(4,wind_emitting_sink)
 Rstar     = xyzmh_ptmass(iReff,wind_emitting_sink)
 Rstar_cgs = Rstar * au 
 Mstar_cgs = Mtotal * solarm 
 Tstar     = xyzmh_ptmass(iTeff,wind_emitting_sink)
 Lstar_cgs = xyzmh_ptmass(iLum,wind_emitting_sink) * unit_luminosity 

 call calc_kappa_max(Mstar_cgs, Lstar_cgs)

 ! Calculate mass distribution
 Matmos = atmos_mass_fraction * Mtotal
 
 inquire(file='mass_loss_rate.dat', exist=file_exists)
   
 if (.not. file_exists) then
    Msink = Mtotal - Matmos
    xyzmh_ptmass(4,wind_emitting_sink) = Msink
 endif

 ! Initialize active boundary spheres
 active_boundary_spheres = iboundary_spheres

 if (iwind == 2 .and. .not. file_exists) then
    call calculate_period(Mtotal, Rstar, pulsation_period_days)
 endif

 ! Setup pulsation parameters
 pulsation_period = pulsation_period_days * (days / utime)
 omega_pulsation = 2.0*pi / pulsation_period
 piston_velocity = piston_velocity_km_s * (km / unit_velocity)
 deltaR_osc = pulsation_period * piston_velocity / (2.0*pi)
 
 ! Setup continuous reinjection period
 reinject_period = reinject_period_days * (days / utime)
 
 ! Setup mass-loss calculation parameters
 mass_loss_start_time = mass_loss_start * (years / utime)
 mass_loss_end_time = mass_loss_end * (years / utime)
 mass_loss_check_radius = check_radius_au  ! Already in code units (AU)
 measurement_interval = meas_int_days * (days / utime)
 time_next_measurement = mass_loss_start_time  ! First measurement at start time
 n_measurements = 0
 
 ! Allocate array for mass-loss rate measurements
 ! Calculate expected number of measurements
 expected_measurements = ceiling((mass_loss_end_time - mass_loss_start_time) / measurement_interval) + 1
 allocate(mass_loss_rates(expected_measurements))
 mass_loss_rates = 0.0
 
 print *, ''
 print *, 'Initializing pulsating atmosphere injection:'
 print *, 'pulsation period: ', pulsation_period
 print *, 'piston velocity: ', piston_velocity
 print *, 'deltaR_osc: ', deltaR_osc
 if (reinject_enabled) then
     print *, '  Reinjection period (days): ', reinject_period_days
     print *, '  Mass-loss calculation window:'
     print *, '    Start time (years): ', mass_loss_start
     print *, '    End time (years): ', mass_loss_end
     print *, '    Check radius (AU): ', check_radius_au
     print *, '    Measurement interval (days): ', meas_int_days
     print *, '    Expected number of measurements: ', expected_measurements
     print *, '  Will calculate mass-loss rate at regular intervals'
 endif
 print *, ''

 if (n_shells > 0) then
    max_shells = n_shells
 else
    max_shells = 200
 endif

 ! Setup stellar structure calculation
 call setup_star(Msink * umass, r_max_on_rstar * Rstar * au, r_min_on_rstar * Rstar * au, gmw, gamma,&
                  surface_pressure, Matmos * umass)
 

 ! Allocate delta_r_radial array
 if (allocated(delta_r_radial)) deallocate(delta_r_radial)
 if (allocated(shell_radii)) deallocate(shell_radii)
 allocate(delta_r_radial(max_shells))
 allocate(shell_radii(max_shells))

 r_min = r_min_on_rstar * Rstar
 current_radius = r_min

 if (n_shells < 0) then
   
   temp_particles = n_particles

   particles_per_shell = nint(real(temp_particles) / real(max_shells))

   converged = .false.

   do while (.not. converged)
      current_radius = r_min
      particles_per_shell = particles_per_shell + 10
      shell_index = 0 

      do while (current_radius < r_max_on_rstar * Rstar)
         shell_index = shell_index + 1
         delta_r_tangential = current_radius * get_fibonacci_spacing(particles_per_shell)
         delta_r_radial(shell_index) = wss * delta_r_tangential
         current_radius = current_radius + delta_r_radial(shell_index)
      end do
      
      distributed_particles = shell_index * particles_per_shell

      if (distributed_particles >= temp_particles) then
         converged = .true.
      endif
   end do 
 else
   temp_particles = 100

   particles_per_shell = nint(real(temp_particles) / real(max_shells))

   converged = .false.

   do while (.not. converged)
      current_radius = r_min
      particles_per_shell = particles_per_shell + 10
      shell_index = 0 

      do while (current_radius < r_max_on_rstar * Rstar)
         shell_index = shell_index + 1
         delta_r_tangential = current_radius * get_fibonacci_spacing(particles_per_shell)
         delta_r_radial(shell_index) = wss * delta_r_tangential
         current_radius = current_radius + delta_r_radial(shell_index)
      end do
      
      distributed_particles = shell_index * particles_per_shell

      if (shell_index >= max_shells) then
         converged = .true.
      endif
   end do 

 endif

 n_shells_total = shell_index
 particles_per_sphere = particles_per_shell
 
 ! Calculate stellar profile for the atmosphere
 call calc_stellar_profile(n_profile_points)

 ! Calculate particle mass from atmospheric mass
 ! Total atmospheric mass distributed over all particles
 mass_of_particles = Matmos / real(n_shells_total * particles_per_sphere)

 if (file_exists) then
    call read_mass_loss_data()
 endif 

 print *, ''
 print *, 'Atmospheric particle mass (Msun): ', mass_of_particles
 print *, 'Total number of shells: ', n_shells_total
 print *, 'Particles per sphere: ', particles_per_sphere
 print *, 'Amount of particles: ', n_shells_total * particles_per_sphere
 print *, ''

 massoftype(igas) = mass_of_particles
 massoftype(iboundary) = mass_of_particles

end subroutine init_inject

!-----------------------------------------------------------------------
!+
!  Main routine: called at the start to setup atmosphere,
!  then called each timestep to handle pulsation and reinjection
!+
!-----------------------------------------------------------------------
subroutine inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npart_old,npartoftype,dtinject)
 use part,        only:igas,iboundary,iamtype

 real,    intent(in)    :: time,dtlast
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(inout) :: npart,npart_old
 integer, intent(inout) :: npartoftype(:)
 real,    intent(out)   :: dtinject

 ! Set timestep constraint for pulsation
 dtinject = pulsation_timestep * pulsation_period

 ! This is neccesary to not re-setup the atmosphere when resuming from a dump
 if (npart > 0 .and. .not. atmosphere_setup_complete) then
    atmosphere_setup_complete = .true.
 endif

 ! Initial setup: create all shells
 if (.not. atmosphere_setup_complete) then
    print *, 'Setting up stellar atmosphere with ', n_shells_total, ' shells.'
    call setup_initial_atmosphere(xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npartoftype)
    atmosphere_setup_complete = .true.
    
    print *, 'Stellar atmosphere setup complete.'
    return
 endif

 ! Reconstruct boundary particle info if resuming from dump and read injection parameters
 if (atmosphere_setup_complete .and. .not. allocated(boundary_particle_ids)) then
    call reconstruct_boundary_info(time, xyzh,npart,xyzmh_ptmass)
    time_last_reinject = time - mod(time, reinject_period)
    ! Try to read mass-loss rate data from file if available
    call read_mass_loss_data()
 endif

 ! Take periodic mass measurements if we're in the measurement period
 if (reinject_enabled .and. .not. mass_loss_rate_calculated) then
    call take_periodic_mass_measurements(time, xyzh, npart, xyzmh_ptmass, npartoftype)
 endif

 ! Check if reinjection is needed (only after mass-loss rate is calculated)
 if (reinject_enabled .and. mass_loss_rate_calculated) then
    
    ! Continuous reinjection mode: check time-based trigger
    call check_continuous_reinject(time,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npartoftype)
    
    ! If reinjection is needed, perform it immediately
    if (reinjection_needed) then
       call perform_reinjection(time,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npartoftype)
       
       ! Update last reinjection time
       time_last_reinject = time
       
       ! Reset the flag after reinjection to prevent continuous reinjection
       reinjection_needed = .false.
       
    endif
 endif

 ! Every subsequent call, move the boundary particles
 call apply_pulsation(time,xyzh,vxyzu,npart,xyzmh_ptmass,vxyz_ptmass)

end subroutine inject_particles

!-----------------------------------------------------------------------
!+
!  Take periodic mass measurements and calculate mean mass-loss rate
!+
!-----------------------------------------------------------------------
subroutine take_periodic_mass_measurements(time, xyzh, npart, xyzmh_ptmass, npartoftype)
 use part,   only:igas,iboundary,iphase,iamtype
 use units,  only:utime,umass
 use physcon,only:solarm,years,days
 
 real,    intent(in) :: time
 real,    intent(in) :: xyzh(:,:), xyzmh_ptmass(:,:)
 integer, intent(in) :: npart
 integer, intent(in) :: npartoftype(:)
 
 real :: current_mass_within_radius, mass_lost, time_elapsed, rate_this_interval
 real :: sink_mass, x0(3), dx, dy, dz, r
 real :: sum_rates, std_dev, mean_rate
 integer :: i
 
 ! Start measurements at the start time
 if (.not. measurement_active .and. time >= mass_loss_start_time) then
    measurement_active = .true.
    
    ! Get sink position
    x0 = xyzmh_ptmass(1:3, wind_emitting_sink)
    sink_mass = xyzmh_ptmass(4, wind_emitting_sink)
    
    ! Count mass of particles within radius for first measurement
    mass_previous_measurement = sink_mass  ! Start with sink mass
    
    do i = 1, npart
       ! Calculate distance from sink
       dx = xyzh(1,i) - x0(1)
       dy = xyzh(2,i) - x0(2)
       dz = xyzh(3,i) - x0(3)
       r = sqrt(dx**2 + dy**2 + dz**2)
       
       ! Only count particles within the check radius
       if (r <= mass_loss_check_radius) then
          mass_previous_measurement = mass_previous_measurement + mass_of_particles
       endif
    enddo
    
    print *, ''
    print *, '========================================='
    print *, 'MASS-LOSS MEASUREMENTS - START'
    print *, '========================================='
    print *, 'Time (years): ', time * utime / years
    print *, 'Mass within ', mass_loss_check_radius, ' AU: ', mass_previous_measurement, ' Msun'
    print *, 'Sink mass: ', sink_mass, ' Msun'
    print *, 'Will take measurements every ', meas_int_days, ' days'
    print *, '========================================='
    print *, ''
    
    ! Set time for next measurement
    time_next_measurement = time + measurement_interval
 endif
 
 ! Take measurement if it's time
 if (measurement_active .and. time >= time_next_measurement .and. time < mass_loss_end_time) then
    
    ! Get sink position
    x0 = xyzmh_ptmass(1:3, wind_emitting_sink)
    sink_mass = xyzmh_ptmass(4, wind_emitting_sink)
    
    ! Count mass of particles within radius
    current_mass_within_radius = sink_mass  ! Start with sink mass
    
    do i = 1, npart
       ! Calculate distance from sink
       dx = xyzh(1,i) - x0(1)
       dy = xyzh(2,i) - x0(2)
       dz = xyzh(3,i) - x0(3)
       r = sqrt(dx**2 + dy**2 + dz**2)
       
       ! Only count particles within the check radius
       if (r <= mass_loss_check_radius) then
          current_mass_within_radius = current_mass_within_radius + mass_of_particles
       endif
    enddo
    
    ! Calculate mass lost since last measurement
    mass_lost = mass_previous_measurement - current_mass_within_radius
    
    ! Calculate instantaneous mass-loss rate for this interval
    rate_this_interval = mass_lost / measurement_interval
    
    ! Store this rate
    n_measurements = n_measurements + 1
    mass_loss_rates(n_measurements) = rate_this_interval
    
    print *, 'Measurement #', n_measurements, ' at t = ', time * utime / years, ' years:'
    print *, '  Mass: ', current_mass_within_radius, ' Msun'
    print *, '  Mass lost since last: ', mass_lost, ' Msun'
    print *, '  Instantaneous rate: ', rate_this_interval * utime / years, ' Msun/yr'
    
    ! Update for next interval
    mass_previous_measurement = current_mass_within_radius
    time_next_measurement = time + measurement_interval
 endif
 
 ! Calculate mean mass-loss rate at end time
 if (measurement_active .and. .not. mass_loss_rate_calculated .and. time >= mass_loss_end_time) then
    
    if (n_measurements > 0) then
       ! Calculate mean of all measurements
       sum_rates = 0.0
       do i = 1, n_measurements
          sum_rates = sum_rates + mass_loss_rates(i)
       enddo
       mean_mass_loss_rate = sum_rates / real(n_measurements)
        
       ! Calculate number of particles to inject per reinjection event
       particles_to_inject = nint((mean_mass_loss_rate * reinject_period) / mass_of_particles)
       
       ! Ensure at least 1 particle is injected
       if (particles_to_inject < 1) particles_to_inject = 1
       
       ! Mark as calculated
       mass_loss_rate_calculated = .true.
       
       print *, ''
       print *, '========================================='
       print *, 'MASS-LOSS RATE CALCULATION COMPLETE'
       print *, '========================================='
       print *, 'Measurement period:'
       print *, '  Start time (years): ', mass_loss_start_time * utime / years
       print *, '  End time (years): ', mass_loss_end_time * utime / years
       print *, '  Duration (years): ', (mass_loss_end_time - mass_loss_start_time) * utime / years
       print *, '  Check radius (AU): ', mass_loss_check_radius
       print *, '  Number of measurements: ', n_measurements
       print *, '  Mean mass-loss rate (Msun/yr): ', mean_mass_loss_rate * utime / years
       print *, ''
       print *, 'Reinjection parameters:'
       print *, '  Reinjection period (days): ', reinject_period_days
       print *, '  Particles to inject per event: ', particles_to_inject
       print *, '  Mass injected per event (Msun): ', particles_to_inject * mass_of_particles
       print *, '========================================='
       print *, ''
       
       ! Write data to file for restart capability
       call write_mass_loss_data()
       
    endif
    
 endif
 
end subroutine take_periodic_mass_measurements

!-----------------------------------------------------------------------
!+
!  Check if it's time for continuous reinjection (time-based)
!+
!-----------------------------------------------------------------------
subroutine check_continuous_reinject(time,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npartoftype)
 use part, only:igas,iboundary,iamtype
 use units, only:utime
 use physcon, only:days 
 
 real,    intent(in)    :: time
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 
 ! Check if enough time has passed since last reinjection
 if ( (time - time_last_reinject) < reinject_period  .and. time >= mass_loss_start_time) then
    return
 endif

 print *, ''
 print *, '========================================='
 print *, 'REINJECTION TRIGGERED'
 print *, 'Time: ', time
 print *, 'Time since last reinject: ', (time - time_last_reinject)*utime/days
 print *, '========================================='
 print *, ''
 
 reinjection_needed = .true.
 
end subroutine check_continuous_reinject
 

!-----------------------------------------------------------------------
!+
!  Perform reinjection: inject new gas particles to replenish lost mass
!+
!-----------------------------------------------------------------------
subroutine perform_reinjection(time,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npartoftype)
 use part,        only:igas,iboundary,iamtype,set_particle_type
 use injectutils, only:inject_fibonacci_sphere
 use wind_pulsating, only:interp_stellar_profile
 use physcon,     only:pi
 
 real,    intent(in)    :: time
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 
 integer :: old_npart, i
 real    :: r_inject, phase, r_dot, rho, u, T, P
 real    :: x0(3), v0(3)
 real    :: total_mass_before, particle_mass_before, sink_mass_before
 real    :: mass_injected
 
 sink_mass_before = xyzmh_ptmass(4, wind_emitting_sink)
 particle_mass_before = npartoftype(igas) * mass_of_particles + npartoftype(iboundary) * mass_of_particles
 total_mass_before = sink_mass_before + particle_mass_before
 
 ! Get sink particle position
 x0 = xyzmh_ptmass(1:3,wind_emitting_sink)
 v0 = vxyz_ptmass(1:3,wind_emitting_sink)
 
 ! Get current phase for pulsation
 phase = omega_pulsation * time + phi0
 r_dot = piston_velocity * cos(phase)
 
 ! Inject new gas particles just outside the outermost boundary sphere
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
 
 ! Get stellar properties at injection radius
 call interp_stellar_profile(r_inject, rho, P, u, T)
 
 ! Store npart before injection
 old_npart = npart

 ! Increase shell number for rotation purposes
 n_reinjections = n_reinjections + 1

 ! Use calculated particles_to_inject based on mean mass-loss rate
 call inject_fibonacci_sphere(n_shells_total + n_reinjections, npart + 1, particles_to_inject, &
                                r_inject, r_dot, u, rho, &
                                npart, npartoftype, xyzh, vxyzu, igas, x0, v0)
 
 ! Calculate mass injected
 mass_injected = (npart - old_npart) * mass_of_particles
 
 ! Subtract this mass from the central sink particle to conserve mass
 xyzmh_ptmass(4, wind_emitting_sink) = xyzmh_ptmass(4, wind_emitting_sink) - mass_injected
 
 print *, ''
 print *, 'Reinjection performed:'
 print *, '  Number of particles injected: ', (npart - old_npart)
 print *, '  Injection radius: ', r_inject
 print *, '  Mass injected (Msun): ', mass_injected
 print *, '  New total particles: ', npart
 print *, '  Relative to original particle count: ', real(npart)/real(particles_per_sphere * n_shells_total)
 print *, 'Reinjection complete.'
 print *, ''
 
end subroutine perform_reinjection


!-----------------------------------------------------------------------
!+
!  Setup initial atmosphere with all shells at t=0
!+
!-----------------------------------------------------------------------
subroutine setup_initial_atmosphere(xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npartoftype)
 use part,        only:igas,iboundary,iphase,iamtype
 use injectutils, only:inject_fibonacci_sphere
 use wind_pulsating, only:interp_stellar_profile
 use physcon,     only:pi,km, au

 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real,    intent(in)    :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)

 integer :: i,j,first_particle,ipart_type,nboundary
 real    :: r,dr(n_shells_total),rho,u,T,P,x0(3),v0(3),GM,v_radial, r_previous
 logical :: is_boundary

 ! Get sink particle position
 x0 = xyzmh_ptmass(1:3,wind_emitting_sink)
 v0 = vxyz_ptmass(1:3,wind_emitting_sink)
 GM = xyzmh_ptmass(4,wind_emitting_sink)

 ! Shell spacing
 dr = delta_r_radial

 r_previous = r_min 

 ! Create shells from inner to outer
 npart = 0
 do i = 1, n_shells_total

    r = (r_previous + delta_r_radial(i))
    r_previous = r

    ! Determine if this is a boundary or free shell
    is_boundary = (i <= iboundary_spheres)
    
    ! Get stellar properties at this radius from 1D stellar profile
    call interp_stellar_profile(r, rho, P, u, T)

    v_radial = 0  ! Initial radial velocity at t=0
    
    ! Set particle type - this tagging ensures forces are handled correctly
    if (is_boundary) then
       ipart_type = iboundary
    else
       ipart_type = igas
    endif
    
    first_particle = npart + 1
    
    call inject_fibonacci_sphere(i, first_particle, particles_per_sphere, r, v_radial, u, rho, &
                                    npart, npartoftype, xyzh, vxyzu, ipart_type, x0, v0)
 enddo

 ! Store information about boundary particles for pulsation
 nboundary = npartoftype(iboundary)
 n_boundary_particles = nboundary

 print *, 'Number of boundary particles: ', nboundary
 print *, 'Number of gas particles: ', npartoftype(igas)
 
 if (nboundary > 0) then
    allocate(r_boundary_equilibrium(nboundary))
    allocate(boundary_particle_ids(nboundary))
    
    ! Store equilibrium radii and IDs of boundary particles
    j = 0
    do i = 1, npart
       if (j >= nboundary) then
          return
       elseif (iamtype(iphase(i)) == iboundary) then
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
 use part, only:iboundary,iphase,iamtype
 use physcon, only:pi
 real,    intent(in) :: time
 real,    intent(in) :: xyzh(:,:),xyzmh_ptmass(:,:)
 integer, intent(in) :: npart
 integer :: i,j
 real :: x0(3), r_current, phase
 
 x0 = xyzmh_ptmass(1:3,wind_emitting_sink)
 
 ! Calculate current phase to remove pulsation displacement
 phase = omega_pulsation * time + phi0
 
 ! Count boundary particles
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
          ! Calculate current radius
          r_current = sqrt((xyzh(1,i)-x0(1))**2 + &
                          (xyzh(2,i)-x0(2))**2 + &
                          (xyzh(3,i)-x0(3))**2)
          ! Remove current pulsation displacement to get equilibrium radius
          ! r_current = r_eq + deltaR_osc * sin(phase)
          r_boundary_equilibrium(j) = r_current - deltaR_osc * sin(phase)
       endif
    enddo
    
    print *, 'Reconstructed boundary particle info from dump:'
    print *, '  Number of boundary particles: ', n_boundary_particles
 endif
 
end subroutine reconstruct_boundary_info

!-----------------------------------------------------------------------
!+
!  Apply pulsation to boundary particles
!+
!-----------------------------------------------------------------------

subroutine apply_pulsation(time,xyzh,vxyzu,npart,xyzmh_ptmass,vxyz_ptmass)
 use physcon, only:pi
 use wind_pulsating, only:interp_stellar_profile

 real,    intent(in)    :: time
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(in)    :: npart

 integer :: i,ipart
 real    :: r_eq,r_new,r_current,phase,piston_velocity_n,deltaR_osc_n
 real    :: x_hat(3),r_dot
 real    :: x0(3),v0(3),GM
 real    :: x, y, z
 real    :: rho,u,T,P  
 
 if (.not. allocated(boundary_particle_ids)) return
 if (n_boundary_particles == 0) return

 ! Get sink particle position
 x0 = xyzmh_ptmass(1:3,wind_emitting_sink)
 v0 = vxyz_ptmass(1:3,wind_emitting_sink)
 GM = xyzmh_ptmass(4,wind_emitting_sink)

 phase = omega_pulsation * time + phi0
 
 if (time < pulsation_period * time_puls .and. time_puls > 0) then
     piston_velocity_n = time * (piston_velocity) / (pulsation_period * time_puls)
 else
     piston_velocity_n = piston_velocity
 endif

 deltaR_osc_n = pulsation_period * piston_velocity_n / (2.0*pi)
 
 ! Pulsation amplitude and velocity
 r_dot = piston_velocity_n * cos(phase)

 ! Update each boundary particle
 do i = 1, n_boundary_particles
    ipart = boundary_particle_ids(i)
    
    ! Equilibrium radius for this particle
    r_eq = r_boundary_equilibrium(i)
    
    ! New radius with pulsation
    r_new = r_eq + deltaR_osc_n * sin(phase)

    x = xyzh(1,ipart) - x0(1)
    y = xyzh(2,ipart) - x0(2)
    z = xyzh(3,ipart) - x0(3)

    r_current = sqrt(x**2 + y**2 + z**2)
    
    ! Radial unit vector
    x_hat(1) = x / r_current
    x_hat(2) = y / r_current
    x_hat(3) = z / r_current

    xyzh(1,ipart) = r_new * x_hat(1) + x0(1)
    xyzh(2,ipart) = r_new * x_hat(2) + x0(2)
    xyzh(3,ipart) = r_new * x_hat(3) + x0(3)

    ! Update velocity (radial pulsation velocity)
    vxyzu(1,ipart) = r_dot * x_hat(1) + v0(1)
    vxyzu(2,ipart) = r_dot * x_hat(2) + v0(2)
    vxyzu(3,ipart) = r_dot * x_hat(3) + v0(3)
    
    ! Update thermodynamic variables based on new radius
    if (var_boundary) then
       call interp_stellar_profile(r_new, rho, P, u, T)
       vxyzu(4,ipart) = u
       xyzh(4,ipart) = (mass_of_particles / rho)**(1./3.)
    endif

 enddo

end subroutine apply_pulsation

subroutine update_injected_par
 ! -- placeholder function
end subroutine update_injected_par

!-----------------------------------------------------------------------
!+
!  Write mass-loss rate data to file for restart
!+
!-----------------------------------------------------------------------
subroutine write_mass_loss_data()
 use io, only:iprint
 integer :: iunit, ierr, i
 
 ! Only write if we have calculated the rate
 if (.not. mass_loss_rate_calculated) return
 
 open(newunit=iunit, file='mass_loss_rate.dat', status='replace', iostat=ierr)
 if (ierr /= 0) then
    write(iprint,*) 'WARNING: Could not write mass_loss_rate.dat'
    return
 endif
 
 write(iunit,*) '# Mass-loss rate data for restart'
 write(iunit,*) mass_loss_rate_calculated
 write(iunit,*) mean_mass_loss_rate
 write(iunit,*) Mtotal                
 write(iunit,*) particles_to_inject
 write(iunit,*) n_measurements
 write(iunit,*) mass_of_particles
 
 ! Write all individual measurements
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
 
 ! Skip comment line
 read(iunit,*)
 
 read(iunit,*, iostat=ierr) mass_loss_rate_calculated
 if (ierr /= 0) then
    close(iunit)
    return
 endif
 
 read(iunit,*, iostat=ierr) mean_mass_loss_rate
 read(iunit,*, iostat=ierr) Mtotal
 read(iunit,*, iostat=ierr) particles_to_inject
 read(iunit,*, iostat=ierr) n_measurements
 read(iunit,*, iostat=ierr) mass_of_particles
 
 ! Allocate and read individual measurements
 if (n_measurements > 0) then
    if (.not. allocated(mass_loss_rates)) allocate(mass_loss_rates(n_measurements))
    do i = 1, n_measurements
       read(iunit,*, iostat=ierr) mass_loss_rates(i)
       if (ierr /= 0) exit
    enddo
 endif
 
 close(iunit)
 
 write(iprint,*) 'Mass-loss rate data read from mass_loss_rate.dat'
 write(iprint,*) '  Mean mass-loss rate: ', mean_mass_loss_rate
 write(iprint,*) '  Total mass: ', Mtotal
 write(iprint,*) '  Particles to inject: ', particles_to_inject
 write(iprint,*) '  Number of measurements: ', n_measurements
 write(iprint,*) '  Mass of particles: ', mass_of_particles
end subroutine read_mass_loss_data

!-----------------------------------------------------------------------
!+
!  Calculate pulsation period based on stellar mass and radius
!+
!-----------------------------------------------------------------------
subroutine calculate_period(M, R, pulsation_period_days)
 real, intent(in)  :: M, R
 real              :: logP, logM, logR
 real, intent(out) :: pulsation_period_days

 print *, 'Calculating pulsation period from mass-radius relation:'
 print *, 'Stellar mass (Msun): ', M
 print *, 'Stellar radius (Rsun): ', R

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

 call write_inopt(n_profile_points,'n_profile_points', 'number of points in stellar profile',iunit)
 call write_inopt(iboundary_spheres,'iboundary_spheres', 'number of boundary spheres (inner layers)',iunit)
 call write_inopt(n_particles,'n_particles', 'number of particles per sphere (if using Fibonacci lattice)',iunit)
 call write_inopt(n_shells,'n_shells', 'number of shells (if <0 determined from n_particles)',iunit)
 call write_inopt(r_min_on_rstar,'r_min_on_rstar', 'inner radius as fraction of R_star',iunit)
 call write_inopt(r_max_on_rstar,'r_max_on_rstar', 'outer radius as fraction of R_star',iunit)
 call write_inopt(atmos_mass_fraction,'atmos_mass_fraction', 'atmospheric mass as fraction of total stellar mass',iunit)
 call write_inopt(surface_pressure,'surface_pressure', 'surface pressure (cgs)',iunit)
 call write_inopt(iwind,'iwind','wind type: 1=prescribed, 2=period from mass-radius relation',iunit)
 call write_inopt(pulsation_period_days,'pulsation_period','pulsation period (days) (if iwind == 2 this is overwritten)',iunit)
 call write_inopt(piston_velocity_km_s,'piston_velocity','piston velocity amplitude (km/s)',iunit)
 call write_inopt(time_puls,'time_puls','time for piston to accelerate from 0 to max velocity (in periods)',iunit)
 call write_inopt(pulsation_timestep,'pulsation_timestep','pulsation timestep as fraction of pulsation period',iunit)
 call write_inopt(phi0,'phi0','initial phase offset (radians) (set to 0. if time_puls > 0)',iunit)
 call write_inopt(wss,'wss','fraction of radial to tangential distance between particles in initial setup',iunit)
 call write_inopt(var_boundary,'var_boundary','allow boundary particles to vary thermodynamic properties (logical)',iunit)
 call write_inopt(reinject_enabled,'reinject_enabled','enable dynamic reinjection of boundary spheres (logical)',iunit)
 call write_inopt(reinject_period_days,'reinject_period_days','period between reinjections in days (for continuous mode)',iunit)
 call write_inopt(mass_loss_start,'mass_loss_start','start time for mass-loss calculation in years',iunit)
 call write_inopt(mass_loss_end,'mass_loss_end','end time for mass-loss calculation in years',iunit)
 call write_inopt(check_radius_au,'check_radius_au','radius in which to count mass for loss rate (AU)',iunit) 
 call write_inopt(meas_int_days,'meas_int_days',&
                  'interval for mass measurements in days',iunit)

end subroutine write_options_inject

!-----------------------------------------------------------------------
!+
!  Read options from input file
!+
!-----------------------------------------------------------------------
subroutine read_options_inject(name,valstring,imatch,igotall,ierr)
 use io, only:fatal
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr

 integer, save :: ngot = 0
 integer, parameter :: noptions = 22
 logical :: init_opt = .false.

 if (.not. init_opt) then
    init_opt = .true.
    call set_default_options_inject()
 endif
 imatch = .true.
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
 case('r_min_on_rstar')
    read(valstring,*,iostat=ierr) r_min_on_rstar
    ngot = ngot + 1
    if (r_min_on_rstar <= 0. .or. r_min_on_rstar >= 1.0) &
       call fatal(label,'r_min_on_rstar must be in range (0,1)')
 case('r_max_on_rstar')
    read(valstring,*,iostat=ierr) r_max_on_rstar
    ngot = ngot + 1
    if (r_max_on_rstar <= 0. .or. r_max_on_rstar > 10.0) &
       call fatal(label,'r_max_on_rstar must be in range (0,10]')
 case('atmos_mass_fraction')
    read(valstring,*,iostat=ierr) atmos_mass_fraction
    ngot = ngot + 1
    if (atmos_mass_fraction <= 0. .or. atmos_mass_fraction >= 1.0) &
       call fatal(label,'atmos_mass_fraction must be in range (0,1)')
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
    if (pulsation_timestep <= 0. .or. pulsation_timestep > 1.0) call fatal(label,'pulsation_timestep must be in range (0,1]')
 case('phi0')
    read(valstring,*,iostat=ierr) phi0
    ngot = ngot + 1
    if (phi0 < -3.1415926536d0 .or. phi0 > 3.1415926536d0) call fatal(label,'phi0 must be in range (-pi,pi)')
 case('wss')
    read(valstring,*,iostat=ierr) wss
    ngot = ngot + 1
    if (wss <= 0. .or. wss > 10.0) call fatal(label,'wss must be in range (0,10]')
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