!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module inject
!
! Handles initial setup of stellar atmosphere with pulsating boundary layers
! and dynamic reinjection when mass drops below thresholds
!
! :References: None
!
! :Owner: Owen Vermeulen
!
! :Runtime parameters:
!   - iboundary_spheres  : *number of boundary spheres (integer)*
!   - n_shells_total     : *total number of atmospheric shells*
!   - iwind_resolution   : *geodesic sphere resolution*
!   - r_min_on_rstar     : *inner radius as fraction of R_star*
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
!   - mass_loss_fraction : *fraction of initial atmospheric mass loss before reinjection*
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
 integer :: iboundary_spheres = 5
 integer :: n_shells_total = 30
 integer :: n_profile_points = 10000
 integer :: iwind_resolution = 30
 integer :: N_particles = 10000
 logical :: use_fibonacci = .true.
 real    :: r_min_on_rstar = 0.9
 real    :: dtpulsation = huge(0.)
 real    :: pulsation_period_days = 300.0  ! Pulsation period in days
 real    :: piston_velocity_km_s = 4.0     ! Piston velocity (in km/s)
 real    :: atmos_mass_fraction = 0.005  ! Atmosphere mass as fraction of total mass
 real    :: surface_pressure = 0.001  ! Surface pressure in cgs units
 integer :: iwind = 1  ! Wind type: 1=prescribed, 2=period from mass-radius relation
 real    :: pulsation_timestep = 0.02
 real    :: phi0 = -3.1415926536d0/2.0  ! Initial phase offset (-pi/2 for starting at minimal radius)
 real    :: wss = 1.0 ! Fraction of the tangential and radial distance between particles in the initial setup
 
 ! Reinjection parameters
 logical :: reinject_enabled = .true.
 real    :: mass_loss_fraction = 0.1  ! Reinject when 10% of initial mass is lost
 real    :: mass_check_radius = -1.0  ! Radius to check for mass loss (negative = use r_max)

! global variables
 integer, parameter :: wind_emitting_sink = 1
 real :: geodesic_R(0:19,3,3), geodesic_v(0:11,3)
 real :: omega_pulsation, deltaR_osc, pulsation_period, piston_velocity
 real :: Rstar, r_min, r_max, mass_of_particles
 real :: Mtotal, Matmos, Msink  ! Total, atmosphere, and sink masses
 real :: Matmos_initial  ! Initial atmospheric mass for tracking
 real, allocatable :: delta_r_radial(:)
 integer :: particles_per_sphere, iresolution
 logical :: atmosphere_setup_complete = .false.
 
 ! Store boundary particle information
 real, allocatable    :: r_boundary_equilibrium(:)
 integer, allocatable :: boundary_particle_ids(:)
 integer              :: n_boundary_particles
 integer              :: active_boundary_spheres  ! Current number of active boundary spheres
 logical              :: reinjection_needed = .false.  ! Flag to track if reinjection should happen
 real                 :: last_phase = 0.0  ! Track last phase to detect maxima

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
 n_shells_total = 30
 n_profile_points = 10000
 iwind_resolution = 30
 N_particles = 10000
 use_fibonacci = .true.
 r_min_on_rstar = 0.9
 dtpulsation = huge(0.)
 atmos_mass_fraction = 0.005
 surface_pressure = 0.001
 iwind = 1
 pulsation_period_days = 300.0
 piston_velocity_km_s = 4.0
 pulsation_timestep = 0.02
 phi0 = -3.1415926536d0/2.0
 wss = 1.0
 reinject_enabled = .true.
 mass_loss_fraction = 0.1
 mass_check_radius = -1.0

end subroutine set_default_options_inject

!-----------------------------------------------------------------------
!+
!  Initialize atmospheric setup and pulsation parameters
!+
!-----------------------------------------------------------------------
subroutine init_inject(ierr)
 use io,            only:fatal
 use physcon,       only:pi,days,au,solarm,km
 use icosahedron,   only:compute_matrices,compute_corners
 use eos,           only:gmw,gamma
 use units,         only:utime,umass,unit_velocity
 use part,          only:xyzmh_ptmass,massoftype,igas,iboundary,nptmass,iTeff,iReff
 use injectutils,   only:get_parts_per_sphere, get_neighb_distance, get_fibonacci_spacing
 use wind_pulsating,only:setup_star,calc_stellar_profile

 integer, intent(out) :: ierr
 real :: Mstar_cgs, Rstar_cgs, Tstar, delta_r_tangential, current_radius
 integer :: i 

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

 ! Calculate mass distribution
 Matmos = atmos_mass_fraction * Mtotal
 Msink  = Mtotal - Matmos
 
 ! Note: Matmos_initial will be set after atmosphere setup
 ! to account for the actual mass within the check radius

 ! Initialize active boundary spheres
 active_boundary_spheres = iboundary_spheres

 if (iwind == 2) then
    call calculate_period(Mtotal, Rstar, pulsation_period_days)
 endif

 ! Setup pulsation parameters
 pulsation_period = pulsation_period_days * (days / utime)
 omega_pulsation = 2.0*pi / pulsation_period
 piston_velocity = piston_velocity_km_s * (km / unit_velocity)
 deltaR_osc = pulsation_period * piston_velocity / (2.0*pi)

 print *, ''
 print *, 'Initializing pulsating atmosphere injection:'
 print *, 'pulsation period: ', pulsation_period
 print *, 'piston velocity: ', piston_velocity
 print *, 'deltaR_osc: ', deltaR_osc
 if (reinject_enabled) then
    print *, 'Reinjection enabled: will reinject when mass drops by', mass_loss_fraction*100., '%'
 endif
 print *, ''
 
 ! Setup sphere parameters
 if (use_fibonacci) then
    print *, 'Using Fibonacci lattice for sphere injection'
    particles_per_sphere = N_particles
 else
    print *, 'Using Geodesic sphere for sphere injection'
    iresolution = iwind_resolution
    call compute_matrices(geodesic_R)
    call compute_corners(geodesic_v)
    particles_per_sphere = get_parts_per_sphere(iresolution)
 endif 

 ! Allocate delta_r_radial array
 if (allocated(delta_r_radial)) deallocate(delta_r_radial)
 allocate(delta_r_radial(n_shells_total))

 r_min = r_min_on_rstar * Rstar
 current_radius = r_min
 
 do i = 1, n_shells_total
    if (use_fibonacci) then
       delta_r_tangential = current_radius * get_fibonacci_spacing(N_particles)
    else
       delta_r_tangential = current_radius * get_neighb_distance(iresolution)
    endif
    delta_r_radial(i) = wss * delta_r_tangential
    current_radius = current_radius + delta_r_radial(i)
 enddo

 r_max = current_radius

 ! Setup stellar structure calculation
 call setup_star(Msink * umass, r_max * au, r_min * au, gmw, gamma, &
                 n_shells_total,surface_pressure, Matmos * umass)
 
 ! Calculate stellar profile for the atmosphere
 call calc_stellar_profile(n_profile_points)

 ! Calculate particle mass from atmospheric mass
 ! Total atmospheric mass distributed over all particles
 mass_of_particles = Matmos / real(n_shells_total * particles_per_sphere)

 print *, ''
 print *, 'Atmospheric particle mass (Msun): ', mass_of_particles
 print *, 'Particles per sphere: ', particles_per_sphere
 print *, 'Amount of particles: ', n_shells_total * particles_per_sphere
 print *, ''

 massoftype(igas) = mass_of_particles
 massoftype(iboundary) = mass_of_particles

 xyzmh_ptmass(4,wind_emitting_sink) = Msink

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
    
    ! Initialize Matmos_initial based on actual mass within check radius
    if (reinject_enabled) then
       call initialize_atmospheric_mass(xyzh,npart,xyzmh_ptmass)
    endif
    
    print *, 'Stellar atmosphere setup complete.'
    return
 endif

 ! Reconstruct boundary particle info if resuming from dump
 if (atmosphere_setup_complete .and. .not. allocated(boundary_particle_ids)) then
    call reconstruct_boundary_info(xyzh,npart,xyzmh_ptmass)
 endif

 ! Check if reinjection is needed
 if (reinject_enabled) then
    call check_and_reinject(time,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npartoftype)
    
    ! If reinjection is needed, perform it immediately
    if (reinjection_needed) then
       call perform_reinjection(time,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npartoftype)
       
       ! CRITICAL: Reset the flag after reinjection to prevent continuous reinjection
       reinjection_needed = .false.
       
       ! Recalculate Matmos_initial after reinjection to include new particles
       call initialize_atmospheric_mass(xyzh,npart,xyzmh_ptmass)
    endif
 endif

 ! Every subsequent call, move the boundary particles
 call apply_pulsation(time,xyzh,vxyzu,npart,xyzmh_ptmass,vxyz_ptmass)

end subroutine inject_particles

!-----------------------------------------------------------------------
!+
!  Check atmospheric mass and trigger reinjection if needed
!+
!-----------------------------------------------------------------------
subroutine check_and_reinject(time,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npartoftype)
 use part, only:igas,iboundary,iphase,iamtype,maxvxyzu,massoftype
 
 real,    intent(in)    :: time
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 
 real :: current_atmos_mass, mass_loss, r_check, r2_check
 real :: x0(3), dx, dy, dz, r2
 integer :: i, itype, n_counted
 integer, save :: ncalls = 0
 logical, save :: first_call = .true.
 
 ncalls = ncalls + 1
 
 ! Radius within which to check for atmospheric mass (in code units)
 ! If mass_check_radius is negative, use r_max (outer edge of initial atmosphere)
 if (mass_check_radius < 0.0) then
    r_check = r_max
 else
    r_check = mass_check_radius
 endif
 r2_check = r_check**2
 
 ! Get sink particle position
 x0 = xyzmh_ptmass(1:3,wind_emitting_sink)
 
 ! Calculate current atmospheric mass within r_check
 current_atmos_mass = 0.0
 n_counted = 0
 do i = 1, npart
    ! Get particle type
    itype = iamtype(iphase(i))
    
    ! Only count gas and boundary particles (not accreted particles, etc.)
    if (itype == igas .or. itype == iboundary) then
       ! Calculate distance from sink
       dx = xyzh(1,i) - x0(1)
       dy = xyzh(2,i) - x0(2)
       dz = xyzh(3,i) - x0(3)
       r2 = dx*dx + dy*dy + dz*dz
       
       ! If within r_check, add to atmospheric mass
       if (r2 < r2_check) then
          current_atmos_mass = current_atmos_mass + massoftype(itype)
          n_counted = n_counted + 1
       endif
    endif
 enddo
 
 ! Calculate mass loss
 mass_loss = Matmos_initial - current_atmos_mass
 
 ! Print diagnostic info on first call or every 100 calls
 if (first_call .or. mod(ncalls, 100) == 0) then
    print *, ''
    print *, '========================================='
    if (first_call) then
       print *, 'REINJECTION MASS CHECK (first call)'
    else
       print *, 'REINJECTION MASS CHECK (call #', ncalls, ')'
    endif
    print *, 'Time: ', time
    print *, 'Check radius: ', r_check
    print *, 'Initial atmospheric mass: ', Matmos_initial
    print *, 'Current atmospheric mass: ', current_atmos_mass
    print *, 'Particles counted: ', n_counted
    print *, 'Mass lost: ', mass_loss
    print *, 'Mass loss fraction: ', mass_loss / Matmos_initial
    print *, 'Trigger threshold: ', mass_loss_fraction
    print *, 'Reinjection enabled: ', reinject_enabled
    if (mass_loss >= mass_loss_fraction * Matmos_initial) then
       print *, 'THRESHOLD MET - will trigger reinjection'
    else
       print *, 'THRESHOLD NOT MET - need to lose ', &
                mass_loss_fraction * Matmos_initial - mass_loss, ' more mass'
    endif
    print *, '========================================='
    print *, ''
    first_call = .false.
 endif
 
 ! Check if we've lost enough mass to trigger reinjection
 if (mass_loss >= mass_loss_fraction * Matmos_initial) then
    if (.not. reinjection_needed) then
       print *, ''
       print *, '========================================='
       print *, 'REINJECTION THRESHOLD MET'
       print *, 'Time: ', time
       print *, 'Check radius: ', r_check
       print *, 'Initial atmospheric mass: ', Matmos_initial
       print *, 'Current atmospheric mass: ', current_atmos_mass
       print *, 'Particles counted: ', n_counted
       print *, 'Mass lost: ', mass_loss
       print *, 'Mass loss fraction: ', mass_loss / Matmos_initial
       print *, '========================================='
       print *, ''
    endif
    
    ! Set flag - actual reinjection will happen at next pulsation maximum
    reinjection_needed = .true.
 endif
 
end subroutine check_and_reinject

!-----------------------------------------------------------------------
!+
!  Initialize Matmos_initial based on actual mass within check radius
!+
!-----------------------------------------------------------------------
subroutine initialize_atmospheric_mass(xyzh,npart,xyzmh_ptmass)
 use part, only:igas,iboundary,iphase,iamtype,massoftype
 
 real,    intent(in) :: xyzh(:,:), xyzmh_ptmass(:,:)
 integer, intent(in) :: npart
 
 real :: r_check, r2_check, r2
 real :: x0(3), dx, dy, dz
 integer :: i, itype, n_counted
 
 ! Determine check radius
 if (mass_check_radius < 0.0) then
    r_check = r_max
 else
    r_check = mass_check_radius
 endif
 r2_check = r_check**2
 
 ! Get sink particle position
 x0 = xyzmh_ptmass(1:3,wind_emitting_sink)
 
 ! Calculate atmospheric mass within check radius
 Matmos_initial = 0.0
 n_counted = 0
 do i = 1, npart
    itype = iamtype(iphase(i))
    
    if (itype == igas .or. itype == iboundary) then
       dx = xyzh(1,i) - x0(1)
       dy = xyzh(2,i) - x0(2)
       dz = xyzh(3,i) - x0(3)
       r2 = dx*dx + dy*dy + dz*dz
       
       if (r2 < r2_check) then
          Matmos_initial = Matmos_initial + massoftype(itype)
          n_counted = n_counted + 1
       endif
    endif
 enddo
 
 print *, ''
 print *, '========================================='
 print *, 'REINJECTION: Initialized atmospheric mass'
 print *, 'Check radius: ', r_check
 print *, 'Particles within radius: ', n_counted, ' / ', npart
 print *, 'Matmos_initial: ', Matmos_initial
 print *, '========================================='
 print *, ''
 
end subroutine initialize_atmospheric_mass


!-----------------------------------------------------------------------
!+
!  Perform reinjection: inject new gas particles to replenish lost mass
!+
!-----------------------------------------------------------------------
subroutine perform_reinjection(time,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npartoftype)
 use part,        only:igas,iboundary,iphase,iamtype,set_particle_type
 use injectutils, only:inject_geodesic_sphere, inject_fibonacci_sphere
 use wind_pulsating, only:interp_stellar_profile
 use physcon,     only:pi
 
 real,    intent(in)    :: time
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 
 integer :: old_npart, i
 real    :: r_inject, phase, r_dot, rho, u, T, P
 real    :: x0(3), v0(3)
 real    :: total_mass_before, total_mass_after, particle_mass_before, sink_mass_before
 real    :: mass_injected
 
 ! Store total mass before reinjection for verification
 sink_mass_before = xyzmh_ptmass(4, wind_emitting_sink)
 particle_mass_before = npartoftype(igas) * mass_of_particles + &
                        npartoftype(iboundary) * mass_of_particles
 total_mass_before = sink_mass_before + particle_mass_before
 
 ! Get sink particle position
 x0 = xyzmh_ptmass(1:3,wind_emitting_sink)
 v0 = vxyz_ptmass(1:3,wind_emitting_sink)
 
 ! Get current phase for pulsation
 phase = omega_pulsation * time + phi0
 r_dot = piston_velocity * cos(phase)
 
 ! Inject new gas particles at the position of the original first gas layer
 ! This is just outside the outermost boundary sphere
 ! Calculate radius: r_min + sum of boundary shell spacings
 r_inject = r_min
 do i = 1, iboundary_spheres
    r_inject = r_inject + delta_r_radial(i)
 enddo
 
 ! Add the pulsation offset to match current phase
 r_inject = r_inject + deltaR_osc * sin(phase)
 
 ! Get stellar properties at injection radius
 call interp_stellar_profile(r_inject, rho, P, u, T)
 
 ! Store npart before injection
 old_npart = npart

 print *, ''
 print *, '========================================='
 print *, 'REINJECTION: Adding new gas particles'
 print *, 'Time: ', time
 print *, 'Boundary spheres (unchanged): ', active_boundary_spheres
 print *, 'Injecting new gas shell at radius: ', r_inject
 print *, 'Phase: ', phase
 print *, 'Pulsation offset: ', deltaR_osc * sin(phase)
 print *, 'Outermost boundary equilibrium radius: ', r_boundary_equilibrium(n_boundary_particles)
 print *, '========================================='
 print *, ''
 
 ! Inject new gas shell
 if (use_fibonacci) then
    call inject_fibonacci_sphere(n_shells_total + 1, npart + 1, particles_per_sphere, &
                                r_inject, r_dot, u, rho, &
                                npart, npartoftype, xyzh, vxyzu, igas, x0, v0)
 else
    call inject_geodesic_sphere(n_shells_total + 1, npart + 1, iresolution, &
                                r_inject, r_dot, u, rho, &
                                geodesic_R, geodesic_V, npart, npartoftype, &
                                xyzh, vxyzu, igas, x0, v0)
 endif
 
 ! Calculate mass injected
 mass_injected = (npart - old_npart) * mass_of_particles
 
 ! Subtract this mass from the central sink particle to conserve mass
 xyzmh_ptmass(4, wind_emitting_sink) = xyzmh_ptmass(4, wind_emitting_sink) - mass_injected
 
 print *, 'Injected ', npart - old_npart, ' new gas particles'
 print *, 'Mass injected: ', mass_injected
 print *, 'Sink mass before: ', sink_mass_before
 print *, 'Sink mass after: ', xyzmh_ptmass(4, wind_emitting_sink)
 print *, 'Reinjection complete.'
 
 ! Verify mass conservation
 total_mass_after = xyzmh_ptmass(4, wind_emitting_sink) + &
                    (npartoftype(igas) * mass_of_particles + &
                     npartoftype(iboundary) * mass_of_particles)
 
 print *, ''
 print *, 'Mass conservation check:'
 print *, 'Total mass before: ', total_mass_before
 print *, 'Total mass after:  ', total_mass_after
 print *, 'Mass difference:   ', abs(total_mass_after - total_mass_before)
 if (abs(total_mass_after - total_mass_before) > 1.0e-10) then
    print *, 'WARNING: Mass not conserved in reinjection!'
 else
    print *, 'Mass conserved successfully.'
 endif
 print *, ''
 
end subroutine perform_reinjection


!-----------------------------------------------------------------------
!+
!  Setup initial atmosphere with all shells at t=0
!+
!-----------------------------------------------------------------------
subroutine setup_initial_atmosphere(xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npartoftype)
 use part,        only:igas,iboundary,iphase,iamtype
 use injectutils, only:inject_geodesic_sphere, inject_fibonacci_sphere
 use wind_pulsating, only:interp_stellar_profile
 use physcon,     only:pi,km, au
 use units,       only:udist, unit_density, unit_ergg, unit_pressure

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
 print *, 'Shell spacing dr:', dr

 r_previous = r_min 

 ! Create shells from inner to outer
 npart = 0
 do i = 1, n_shells_total

    ! Calculate radius for this shell
    r = (r_previous + delta_r_radial(i))
    r_previous = r

    ! Determine if this is a boundary or free shell
    is_boundary = (i <= iboundary_spheres)
    
    ! Get stellar properties at this radius from 1D stellar profile
    call interp_stellar_profile(r, rho, P, u, T)

    v_radial = 0.0  ! Boundary shells have zero initial velocity
    
    ! Set particle type
    if (is_boundary) then
       ipart_type = iboundary
    else
       ipart_type = igas
    endif
    
    ! Inject this shell
    first_particle = npart + 1
    if (use_fibonacci) then
        call inject_fibonacci_sphere(i, first_particle, particles_per_sphere, r, v_radial, u, rho, &
                                    npart, npartoftype, xyzh, vxyzu, ipart_type, x0, v0)
    else
        call inject_geodesic_sphere(i, first_particle, iresolution, r, v_radial, u, rho, &
                                geodesic_R, geodesic_V, npart, npartoftype, &
                                xyzh, vxyzu, ipart_type, x0, v0)
    endif
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
subroutine reconstruct_boundary_info(xyzh,npart,xyzmh_ptmass)
 use part, only:iboundary,iphase,iamtype
 real,    intent(in) :: xyzh(:,:),xyzmh_ptmass(:,:)
 integer, intent(in) :: npart
 integer :: i,j
 real :: x0(3)
 
 x0 = xyzmh_ptmass(1:3,wind_emitting_sink)
 
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
          ! Store current radius as equilibrium (assumes restart at a similar phase)
          r_boundary_equilibrium(j) = sqrt((xyzh(1,i)-x0(1))**2 + &
                                           (xyzh(2,i)-x0(2))**2 + &
                                           (xyzh(3,i)-x0(3))**2)
       endif
    enddo
    
    print *, 'Reconstructed boundary particle info from dump:'
    print *, '  Number of boundary particles: ', n_boundary_particles
 endif
 
end subroutine reconstruct_boundary_info

!-----------------------------------------------------------------------
!+
!  Apply radial pulsation to boundary particles
!+
!-----------------------------------------------------------------------
subroutine apply_pulsation(time,xyzh,vxyzu,npart,xyzmh_ptmass,vxyz_ptmass)
 use physcon, only:pi

 real,    intent(in)    :: time
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(in)    :: npart

 integer :: i,ipart
 real    :: r_eq,r_new,r_current,phase
 real    :: x_hat(3),r_dot
 real    :: x0(3),v0(3),GM
 real    :: x, y, z
 
 if (.not. allocated(boundary_particle_ids)) return
 if (n_boundary_particles == 0) return

 ! Get sink particle position
 x0 = xyzmh_ptmass(1:3,wind_emitting_sink)
 v0 = vxyz_ptmass(1:3,wind_emitting_sink)
 GM = xyzmh_ptmass(4,wind_emitting_sink)

 phase = omega_pulsation * time + phi0
 
 ! Pulsation velocity
 r_dot = piston_velocity * cos(phase)

 ! Update each boundary particle
 do i = 1, n_boundary_particles
    ipart = boundary_particle_ids(i)
    
    ! Equilibrium radius for this particle
    r_eq = r_boundary_equilibrium(i)
    
    ! New radius with pulsation
    r_new = r_eq + deltaR_osc * sin(phase)

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

    ! Update velocity (radial pulsation velocity + orbital velocity of sink)
    vxyzu(1,ipart) = r_dot * x_hat(1) + v0(1)
    vxyzu(2,ipart) = r_dot * x_hat(2) + v0(2)
    vxyzu(3,ipart) = r_dot * x_hat(3) + v0(3)
 enddo

end subroutine apply_pulsation

subroutine update_injected_par
 ! -- placeholder function
end subroutine update_injected_par

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
 call write_inopt(n_shells_total,'n_shells_total', 'total number of atmospheric shells',iunit)
 call write_inopt(iboundary_spheres,'iboundary_spheres', 'number of boundary spheres (inner layers)',iunit)
 call write_inopt(iwind_resolution,'iwind_resolution', 'geodesic sphere resolution (integer)',iunit)
 call write_inopt(N_particles,'N_particles', 'number of particles per sphere (if using Fibonacci lattice)',iunit)
 call write_inopt(use_fibonacci,'use_fibonacci', 'use Fibonacci lattice for sphere injection (logical)',iunit)
 call write_inopt(r_min_on_rstar,'r_min_on_rstar', 'inner radius as fraction of R_star',iunit)
 call write_inopt(atmos_mass_fraction,'atmos_mass_fraction', 'atmospheric mass as fraction of total stellar mass',iunit)
 call write_inopt(surface_pressure,'surface_pressure', 'surface pressure (cgs)',iunit)
 call write_inopt(iwind,'iwind','wind type: 1=prescribed, 2=period from mass-radius relation',iunit)
 call write_inopt(pulsation_period_days,'pulsation_period','pulsation period (days) (if iwind == 2 this is overwritten)',iunit)
 call write_inopt(piston_velocity_km_s,'piston_velocity','piston velocity amplitude (km/s)',iunit)
 call write_inopt(pulsation_timestep,'pulsation_timestep','pulsation timestep as fraction of pulsation period',iunit)
 call write_inopt(phi0,'phi0','initial phase offset (radians)',iunit)
 call write_inopt(wss,'wss','fraction of radial to tangential distance between particles in initial setup',iunit)
 call write_inopt(reinject_enabled,'reinject_enabled','enable dynamic reinjection of boundary spheres (logical)',iunit)
 call write_inopt(mass_loss_fraction,'mass_loss_fraction','fraction of atmospheric mass loss before reinjection',iunit)
 call write_inopt(mass_check_radius,'mass_check_radius','radius to check for mass loss (-1 = r_max)',iunit)
 
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
 integer, parameter :: noptions = 17
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
 case('n_shells_total')
    read(valstring,*,iostat=ierr) n_shells_total
    ngot = ngot + 1
    if (n_shells_total <= 0) call fatal(label,'n_shells_total must be > 0')
 case('iboundary_spheres')
    read(valstring,*,iostat=ierr) iboundary_spheres
    ngot = ngot + 1
    if (iboundary_spheres < 0) call fatal(label,'iboundary_spheres must be >= 0')
    if (iboundary_spheres > n_shells_total) &
       call fatal(label,'iboundary_spheres must be <= n_shells_total')
 case('iwind_resolution')
    read(valstring,*,iostat=ierr) iwind_resolution
    ngot = ngot + 1
    if (iwind_resolution < 1) call fatal(label,'iwind_resolution must be >= 1')
 case('N_particles')
    read(valstring,*,iostat=ierr) N_particles
    ngot = ngot + 1
    if (N_particles < 1) call fatal(label,'N_particles must be >= 1')
 case('use_fibonacci')
    read(valstring,*,iostat=ierr) use_fibonacci
    ngot = ngot + 1
 case('r_min_on_rstar')
    read(valstring,*,iostat=ierr) r_min_on_rstar
    ngot = ngot + 1
    if (r_min_on_rstar <= 0. .or. r_min_on_rstar >= 1.0) &
       call fatal(label,'r_min_on_rstar must be in range (0,1)')
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
 case('reinject_enabled')
    read(valstring,*,iostat=ierr) reinject_enabled
    ngot = ngot + 1
 case('mass_loss_fraction')
    read(valstring,*,iostat=ierr) mass_loss_fraction
    ngot = ngot + 1
    if (mass_loss_fraction <= 0. .or. mass_loss_fraction >= 1.0) &
       call fatal(label,'mass_loss_fraction must be in range (0,1)')
 case('mass_check_radius')
    read(valstring,*,iostat=ierr) mass_check_radius
    ngot = ngot + 1
 case default
    imatch = .false.
 end select

 igotall = (ngot >= noptions)

end subroutine read_options_inject

end module inject