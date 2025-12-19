!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module inject
!
! Handles initial setup of stellar atmosphere with pulsating boundary layers
! Extended to support non-radial pulsations via spherical harmonics
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
!   - enable_nonradial   : *enable non-radial pulsations (0=off, 1=on)*
!   - l_mode             : *spherical harmonic degree l*
!   - m_mode             : *spherical harmonic order m*
!   - nonradial_amplitude: *amplitude of non-radial mode as fraction of radial amplitude*
!   - nonradial_phase    : *phase offset for non-radial mode (radians)*
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
 integer :: iboundary_spheres = 10
 integer :: n_shells_total = 50
 integer :: n_profile_points = 10000
 integer :: iwind_resolution = 30
 real    :: r_min_on_rstar = 0.9
 real    :: dtpulsation = huge(0.)
 real    :: pulsation_period_days = 300.0  ! Pulsation period in days
 real    :: piston_velocity_km_s = 4.0     ! Piston velocity (in km/s)
 real    :: atmos_mass_fraction = 0.005  ! Atmosphere mass as fraction of total mass
 real    :: surface_pressure = 0.001  ! Surface pressure in cgs units
 integer :: iwind = 1  ! Wind type: 1=prescribed, 2=period from mass-radius relation
 real    :: pulsation_timestep = 0.02
 real    :: phi0 = -3.1415926536d0/2.0  ! Initial phase offset (-pi/2 for starting at minimal radius)
 real    :: wss = 2.0 ! Fraction of the tangential and radial distance between particles in the initial setup

 ! Non-radial pulsation parameters
 integer :: enable_nonradial = 0  ! 0=disabled, 1=enabled
 integer :: l_mode = 2  ! Spherical harmonic degree
 integer :: m_mode = 0  ! Spherical harmonic order
 real    :: nonradial_amplitude = 0.1  ! Amplitude relative to radial pulsation
 real    :: nonradial_phase = 0.0  ! Phase offset for non-radial mode

! global variables
 integer, parameter :: wind_emitting_sink = 1
 real :: geodesic_R(0:19,3,3), geodesic_v(0:11,3)
 real :: omega_pulsation, deltaR_osc, pulsation_period, piston_velocity
 real :: Rstar, r_min, r_max, mass_of_particles
 real :: Mtotal, Matmos, Msink  ! Total, atmosphere, and sink masses
 real, allocatable :: delta_r_radial(:)
 integer :: particles_per_sphere, iresolution
 logical :: atmosphere_setup_complete = .false.
 
 ! Store boundary particle information
 real, allocatable    :: r_boundary_equilibrium(:)
 real, allocatable    :: theta_boundary(:), phi_boundary(:)  ! Spherical coords
 integer, allocatable :: boundary_particle_ids(:)
 integer              :: n_boundary_particles

 character(len=*), parameter :: label = 'inject_atmosphere'

contains

!-----------------------------------------------------------------------
!+
!  Set default options
!+
!-----------------------------------------------------------------------
subroutine set_default_options_inject(flag)
 integer, optional, intent(in) :: flag

 iboundary_spheres = 10
 n_shells_total = 50
 n_profile_points = 10000
 iwind_resolution = 30
 r_min_on_rstar = 0.9
 dtpulsation = huge(0.)
 atmos_mass_fraction = 0.005
 surface_pressure = 0.001
 iwind = 1
 pulsation_period_days = 300.0
 piston_velocity_km_s = 4.0
 pulsation_timestep = 0.02
 phi0 = -3.1415926536d0/2.0
 wss = 2.0
 
 ! Non-radial defaults
 enable_nonradial = 0
 l_mode = 2
 m_mode = 0
 nonradial_amplitude = 0.1
 nonradial_phase = 0.0

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
 use injectutils,   only:get_parts_per_sphere, get_neighb_distance
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
 
 if (enable_nonradial == 1) then
    print *, ''
    print *, 'Non-radial pulsations ENABLED:'
    print *, '  l mode: ', l_mode
    print *, '  m mode: ', m_mode
    print *, '  amplitude (fraction of radial): ', nonradial_amplitude
    print *, '  phase offset: ', nonradial_phase
 else
    print *, 'Non-radial pulsations DISABLED'
 endif
 print *, ''
 
 ! Setup geodesic sphere parameters
 iresolution = iwind_resolution
 call compute_matrices(geodesic_R)
 call compute_corners(geodesic_v)
 particles_per_sphere = get_parts_per_sphere(iresolution)

 ! Allocate delta_r_radial array
 if (allocated(delta_r_radial)) deallocate(delta_r_radial)
 allocate(delta_r_radial(n_shells_total))

 r_min = r_min_on_rstar * Rstar
 current_radius = r_min
 
 do i = 1, n_shells_total
    delta_r_tangential = current_radius * get_neighb_distance(iresolution)
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
!  then called each timestep to handle pulsation
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

 ! Every subsequent call, move the boundary particles
 if (enable_nonradial == 1) then
    call apply_pulsation_nonradial(time,xyzh,vxyzu,npart,xyzmh_ptmass,vxyz_ptmass)
 else
    call apply_pulsation(time,xyzh,vxyzu,npart,xyzmh_ptmass,vxyz_ptmass)
 endif

end subroutine inject_particles

!-----------------------------------------------------------------------
!+
!  Setup initial atmosphere with all shells at t=0
!+
!-----------------------------------------------------------------------
subroutine setup_initial_atmosphere(xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npartoftype)
 use part,        only:igas,iboundary,iphase,iamtype
 use injectutils, only:inject_geodesic_sphere
 use wind_pulsating, only:interp_stellar_profile
 use physcon,     only:pi,km, au
 use units,       only:udist, unit_density, unit_ergg, unit_pressure

 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real,    intent(in)    :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)

 integer :: i,j,first_particle,ipart_type,nboundary
 real    :: r,dr(n_shells_total),rho,u,T,P,x0(3),v0(3),GM,v_radial, r_previous
 real    :: theta, phi_angle, x_rel, y_rel, z_rel
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

    v_radial = 0.0
    
    ! Set particle type
    if (is_boundary) then
       ipart_type = iboundary
    else
       ipart_type = igas
    endif
    
    ! Inject this shell using geodesic sphere
    first_particle = npart + 1
    call inject_geodesic_sphere(i, first_particle, iresolution, r, v_radial, u, rho, &
                                geodesic_R, geodesic_V, npart, npartoftype, &
                                xyzh, vxyzu, ipart_type, x0, v0)
 enddo

 ! Store information about boundary particles for pulsation
 nboundary = npartoftype(iboundary)
 n_boundary_particles = nboundary

 print *, 'Number of boundary particles: ', nboundary
 print *, 'Number of gas particles: ', npartoftype(igas)
 
 if (nboundary > 0) then
    allocate(r_boundary_equilibrium(nboundary))
    allocate(theta_boundary(nboundary))
    allocate(phi_boundary(nboundary))
    allocate(boundary_particle_ids(nboundary))
    
    ! Store equilibrium radii, IDs, and spherical coordinates of boundary particles
    j = 0
    do i = 1, npart
       if (j >= nboundary) then
          return
       elseif (iamtype(iphase(i)) == iboundary) then
          j = j + 1
          boundary_particle_ids(j) = i
          
          ! Calculate relative position
          x_rel = xyzh(1,i) - x0(1)
          y_rel = xyzh(2,i) - x0(2)
          z_rel = xyzh(3,i) - x0(3)
          
          r_boundary_equilibrium(j) = sqrt(x_rel**2 + y_rel**2 + z_rel**2) &
                                            - deltaR_osc * sin(phi0)
          
          ! Calculate spherical coordinates (theta, phi)
          theta = acos(z_rel / sqrt(x_rel**2 + y_rel**2 + z_rel**2))
          phi_angle = atan2(y_rel, x_rel)
          
          theta_boundary(j) = theta
          phi_boundary(j) = phi_angle
       endif
    enddo
 endif

end subroutine setup_initial_atmosphere

!-----------------------------------------------------------------------
!+
!  Apply radial pulsation to boundary particles (original version)
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
 
 ! Pulsation amplitude and velocity
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

    ! Update velocity (radial pulsation velocity)
    vxyzu(1,ipart) = r_dot * x_hat(1) + v0(1)
    vxyzu(2,ipart) = r_dot * x_hat(2) + v0(2)
    vxyzu(3,ipart) = r_dot * x_hat(3) + v0(3)
 enddo

end subroutine apply_pulsation

!-----------------------------------------------------------------------
!+
!  Apply radial + non-radial pulsation to boundary particles
!+
!-----------------------------------------------------------------------
subroutine apply_pulsation_nonradial(time,xyzh,vxyzu,npart,xyzmh_ptmass,vxyz_ptmass)
 use physcon, only:pi

 real,    intent(in)    :: time
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(in)    :: npart

 integer :: i,ipart
 real    :: r_eq,r_new,r_current,phase
 real    :: x_hat(3),r_dot
 real    :: x0(3),v0(3),GM
 real    :: x, y, z
 real    :: theta, phi_angle
 real    :: Y_lm, dY_lm_dt
 real    :: delta_r_nonradial, v_nonradial
 
 if (.not. allocated(boundary_particle_ids)) return
 if (n_boundary_particles == 0) return

 ! Get sink particle position
 x0 = xyzmh_ptmass(1:3,wind_emitting_sink)
 v0 = vxyz_ptmass(1:3,wind_emitting_sink)
 GM = xyzmh_ptmass(4,wind_emitting_sink)

 phase = omega_pulsation * time + phi0
 
 ! Radial pulsation amplitude and velocity
 r_dot = piston_velocity * cos(phase)

 ! Update each boundary particle
 do i = 1, n_boundary_particles
    ipart = boundary_particle_ids(i)
    
    ! Equilibrium radius for this particle
    r_eq = r_boundary_equilibrium(i)
    
    ! Get spherical coordinates for this particle
    theta = theta_boundary(i)
    phi_angle = phi_boundary(i)
    
    ! Calculate spherical harmonic for this position
    call spherical_harmonic(l_mode, m_mode, theta, phi_angle, Y_lm)
    
    ! Non-radial perturbation
    ! delta_r = A * Y_lm * sin(omega*t + phi_nr)
    delta_r_nonradial = nonradial_amplitude * deltaR_osc * Y_lm * &
                        sin(omega_pulsation * time + nonradial_phase)
    
    ! Time derivative for velocity
    v_nonradial = nonradial_amplitude * piston_velocity * Y_lm * &
                  cos(omega_pulsation * time + nonradial_phase)
    
    ! Combined radius: radial + non-radial
    r_new = r_eq + deltaR_osc * sin(phase) + delta_r_nonradial

    x = xyzh(1,ipart) - x0(1)
    y = xyzh(2,ipart) - x0(2)
    z = xyzh(3,ipart) - x0(3)

    r_current = sqrt(x**2 + y**2 + z**2)
    
    ! Radial unit vector
    x_hat(1) = x / r_current
    x_hat(2) = y / r_current
    x_hat(3) = z / r_current

    ! Update position with combined pulsation
    xyzh(1,ipart) = r_new * x_hat(1) + x0(1)
    xyzh(2,ipart) = r_new * x_hat(2) + x0(2)
    xyzh(3,ipart) = r_new * x_hat(3) + x0(3)

    ! Update velocity with combined pulsation velocity
    vxyzu(1,ipart) = (r_dot + v_nonradial) * x_hat(1) + v0(1)
    vxyzu(2,ipart) = (r_dot + v_nonradial) * x_hat(2) + v0(2)
    vxyzu(3,ipart) = (r_dot + v_nonradial) * x_hat(3) + v0(3)
 enddo

end subroutine apply_pulsation_nonradial

!-----------------------------------------------------------------------
!+
!  Calculate real spherical harmonic Y_l^m(theta, phi)
!  Uses simplified formulas for common modes
!+
!-----------------------------------------------------------------------
subroutine spherical_harmonic(l, m, theta, phi, Y_lm)
 use physcon, only:pi
 integer, intent(in) :: l, m
 real, intent(in)    :: theta, phi
 real, intent(out)   :: Y_lm
 
 real :: cos_theta, sin_theta, P_lm
 real :: norm
 
 cos_theta = cos(theta)
 sin_theta = sin(theta)
 
 ! Calculate associated Legendre polynomial
 call associated_legendre(l, abs(m), cos_theta, P_lm)
 
 ! Normalization factor
 norm = legendre_normalization(l, abs(m))
 
 ! Combine with azimuthal part
 if (m > 0) then
    Y_lm = norm * P_lm * cos(m * phi)
 elseif (m < 0) then
    Y_lm = norm * P_lm * sin(abs(m) * phi)
 else
    Y_lm = norm * P_lm
 endif

end subroutine spherical_harmonic

!-----------------------------------------------------------------------
!+
!  Calculate associated Legendre polynomial P_l^m(x)
!  Using recurrence relations for efficiency
!+
!-----------------------------------------------------------------------
subroutine associated_legendre(l, m, x, P_lm)
 integer, intent(in) :: l, m
 real, intent(in)    :: x
 real, intent(out)   :: P_lm
 
 real :: P_mm, P_mm1, P_ll
 integer :: i, ll
 real :: somx2
 
 if (m < 0 .or. m > l .or. abs(x) > 1.0) then
    P_lm = 0.0
    return
 endif
 
 ! Calculate P_m^m
 P_mm = 1.0
 if (m > 0) then
    somx2 = sqrt((1.0 - x) * (1.0 + x))
    do i = 1, m
       P_mm = P_mm * (-1.0) * (2.0*i - 1.0) * somx2
    enddo
 endif
 
 if (l == m) then
    P_lm = P_mm
    return
 endif
 
 ! Calculate P_m^(m+1)
 P_mm1 = x * (2.0*m + 1.0) * P_mm
 
 if (l == m + 1) then
    P_lm = P_mm1
    return
 endif
 
 ! Use recurrence for P_l^m with l > m+1
 do ll = m + 2, l
    P_ll = (x * (2.0*ll - 1.0) * P_mm1 - (ll + m - 1.0) * P_mm) / real(ll - m)
    P_mm = P_mm1
    P_mm1 = P_ll
 enddo
 
 P_lm = P_mm1

end subroutine associated_legendre

!-----------------------------------------------------------------------
!+
!  Calculate normalization factor for spherical harmonics
!+
!-----------------------------------------------------------------------
function legendre_normalization(l, m) result(norm)
 use physcon, only:pi
 integer, intent(in) :: l, m
 real :: norm
 integer :: i
 real :: factor
 
 ! Calculate (l-m)! / (l+m)!
 factor = 1.0
 do i = l - m + 1, l + m
    factor = factor / real(i)
 enddo
 
 ! Full normalization: sqrt((2l+1)/(4*pi) * (l-m)!/(l+m)!)
 norm = sqrt((2.0*l + 1.0) / (4.0 * pi) * factor)
 
 ! For real spherical harmonics, add factor of sqrt(2) for m != 0
 if (m /= 0) then
    norm = norm * sqrt(2.0)
 endif

end function legendre_normalization

subroutine update_injected_par
 ! -- placeholder function
end subroutine update_injected_par

!-----------------------------------------------------------------------
!+
!  Calculate pulsation period based on stellar mass and radius
!  Using empirical relation from Ostlie & Cox (1986)
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
 call write_inopt(r_min_on_rstar,'r_min_on_rstar', 'inner radius as fraction of R_star',iunit)
 call write_inopt(atmos_mass_fraction,'atmos_mass_fraction', 'atmospheric mass as fraction of total stellar mass',iunit)
 call write_inopt(surface_pressure,'surface_pressure', 'surface pressure (cgs)',iunit)
 call write_inopt(iwind,'iwind','wind type: 1=prescribed, 2=period from mass-radius relation',iunit)
 call write_inopt(pulsation_period_days,'pulsation_period','pulsation period (days) (if iwind == 2 this is overwritten)',iunit)
 call write_inopt(piston_velocity_km_s,'piston_velocity','piston velocity amplitude (km/s)',iunit)
 call write_inopt(pulsation_timestep,'pulsation_timestep','pulsation timestep as fraction of pulsation period',iunit)
 call write_inopt(phi0,'phi0','initial phase offset (radians)',iunit)
 call write_inopt(wss,'wss','fraction of radial to tangential distance between particles in initial setup',iunit)
 
 ! Non-radial options
 call write_inopt(enable_nonradial,'enable_nonradial','enable non-radial pulsations (0=off, 1=on)',iunit)
 call write_inopt(l_mode,'l_mode','spherical harmonic degree l',iunit)
 call write_inopt(m_mode,'m_mode','spherical harmonic order m',iunit)
 call write_inopt(nonradial_amplitude,'nonradial_amplitude','amplitude of non-radial mode as fraction of radial amplitude',iunit)
 call write_inopt(nonradial_phase,'nonradial_phase','phase offset for non-radial mode (radians)',iunit)

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
 integer, parameter :: noptions = 17  ! Updated from 12
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
 ! Non-radial options
 case('enable_nonradial')
    read(valstring,*,iostat=ierr) enable_nonradial
    ngot = ngot + 1
    if (enable_nonradial /= 0 .and. enable_nonradial /= 1) &
       call fatal(label,'enable_nonradial must be 0 or 1')
 case('l_mode')
    read(valstring,*,iostat=ierr) l_mode
    ngot = ngot + 1
    if (l_mode < 0) call fatal(label,'l_mode must be >= 0')
 case('m_mode')
    read(valstring,*,iostat=ierr) m_mode
    ngot = ngot + 1
    if (abs(m_mode) > l_mode) call fatal(label,'|m_mode| must be <= l_mode')
 case('nonradial_amplitude')
    read(valstring,*,iostat=ierr) nonradial_amplitude
    ngot = ngot + 1
    if (nonradial_amplitude < 0. .or. nonradial_amplitude > 5.0) &
       call fatal(label,'nonradial_amplitude must be in range [0,5]')
 case('nonradial_phase')
    read(valstring,*,iostat=ierr) nonradial_phase
    ngot = ngot + 1
    if (nonradial_phase < -3.1415926536d0 .or. nonradial_phase > 3.1415926536d0) &
       call fatal(label,'nonradial_phase must be in range (-pi,pi)')
 case default
    imatch = .false.
 end select

 igotall = (ngot >= noptions)

end subroutine read_options_inject

end module inject