!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module inject
!
! Handles pulsating AGB stars
!
! The boundary layers follow a radial pulsation plus an optional sum of
! non-radial modes (real spherical harmonics), all displacing particles
! radially along their unit vector from the wind-emitting sink:
!
!   r_i(t) = r_eq,i + dR sin(w t + phi0) + sum_k dR_k Y_k(theta_i,phi_i) sin(w_k t + phi_k)
!   v_i(t) = v_sink + [ v_p cos(w t + phi0) + sum_k w_k dR_k Y_k cos(w_k t + phi_k) ] rhat_i
!
! with dR = v_p/w and dR_k = A_k dR (A_k = mode_amplitudes). theta, phi are measured about the
! z-axis of the simulation, centred on the sink. The gas shells built at setup
! and every reinjected shell are shifted by the same non-radial displacement
! as the boundary, so they keep a constant spacing above the deformed surface.
!
! :References: None
!
! :Owner: Owen Vermeulen
!
! :Runtime parameters:
!   - dumps_p_period    : *number of dumps per period (if save_period = 1)*
!   - iboundary_spheres : *number of boundary spheres (piston layers)*
!   - iwind             : *wind type: 1=prescribed, 2=period from mass-radius relation*
!   - l_modes           : *spherical harmonic degree l per mode (comma-separated)*
!   - m_modes           : *spherical harmonic order m per mode (comma-separated)*
!   - mode_amplitudes   : *amplitude per mode, relative to the radial amplitude dR*
!   - mass_loss_end     : *end time for mass-loss calculation (periods)*
!   - mass_loss_start   : *start time for mass-loss calculation (periods)*
!   - min_mass_fraction : *minimum mass fraction per shell*
!   - mode_periods      : *period per mode (days), <=0 uses the radial period*
!   - mode_phases       : *phase offset per mode (radians)*
!   - n_inject_period   : *period between reinjections (periods)*
!   - n_modes           : *number of non-radial modes (0=radial pulsation only)*
!   - n_profile_points  : *number of points in stellar profile*
!   - phi0              : *initial phase offset (radians)*
!   - piston_velocity   : *piston velocity amplitude (km/s)*
!   - pulsation_period  : *pulsation period (days)*
!   - r_min_on_rstar    : *gas atmosphere inner radius as fraction of R_star*
!   - reinject_enabled  : *enable dynamic reinjection (0=off, 1=on)*
!   - rho_inner         : *inner boundary density at r_min (cgs)*
!   - rho_power         : *density profile exponent: rho ~ r^(-rho_power)*
!   - save_period       : *wether to save dumps as fraction of period (0=off, 1=on)*
!   - update_L          : *update luminosity with pulsation (0=off, 1=on)*
!   - use_file_mdot     : *skip measurement phase (0=off, 1=on)*
!   - verbose           : *enable verbose output (0=off, 1=on)*
!   - wss               : *radial/tangential spacing ratio*
!
! :Dependencies: dust_formation, eos, infile_utils, injectutils, io, part,
!   mathfunc, physcon, timestep, units, wind_pulsating
!
 use io,      only:fatal
 use physcon, only:pi
 implicit none
 character(len=*), parameter, public :: inject_type = 'atmosphere'

 public :: init_inject, inject_particles, write_options_inject, read_options_inject, &
           set_default_options_inject, update_injected_par
 private

 integer :: iboundary_spheres     = 5
 integer :: n_profile_points      = 10000
 real    :: min_mass_fraction     = 0.01
 real    :: rho_power_in          = 6.0
 real    :: r_min_on_rstar        = 1.0
 real    :: rho_inner             = 1.0e-12
 integer :: iwind                 = 1
 real    :: pulsation_period_days = 300.0
 real    :: piston_velocity_km_s  = 4.0
 real    :: phi0                  = -0.5*pi
 real    :: wss                   = 1.0
 integer :: save_period           = 0
 integer :: dumps_p_period        = 10
 integer :: reinject_enabled      = 1
 integer :: n_inject_period       = 40
 real    :: mass_loss_start       = 6.0
 real    :: mass_loss_end         = 8.0
 integer :: update_L              = 0
 integer :: verbose               = 1
 integer :: use_file_mdot         = 0

 ! non-radial modes; max_modes is limited by the length of a value in the .in file
 integer, parameter :: max_modes = 8
 integer :: n_modes                         = 0
 integer :: l_modes(max_modes)              = 0
 integer :: m_modes(max_modes)              = 0
 real    :: mode_amplitudes(max_modes)      = 0.1   ! amplitudes relative to the radial pulsation
 real    :: mode_period_days(max_modes)     = 0.
 real    :: mode_phases(max_modes)          = 0.

 integer,          parameter :: wind_emitting_sink = 1
 real,             parameter :: pulsation_timestep = 0.02
 integer,          parameter :: maxlen_modelist    = 95
 character(len=*), parameter :: label = 'inject_atmosphere'

 real    :: omega_pulsation, deltaR_osc, pulsation_period, piston_velocity
 real    :: shortest_period   ! shortest of the radial and (non-zero) mode periods, sets dtinject
 real    :: r_min, r_max, Mtotal, mass_of_gas_particle
 integer :: n_shells_total, n_shells_bnd
 integer, allocatable :: npart_per_shell(:), npart_per_boundary_shell(:)
 real,    allocatable :: shell_radii_gas(:), shell_radii_bnd(:), delta_r_radial(:)

 logical :: atmosphere_setup_complete = .false.
 logical :: boundary_info_ready       = .false.
 integer :: n_boundary_particles      = 0
 real, allocatable :: r_boundary_equilibrium(:)
 real    :: mode_omega(max_modes), mode_deltaR(max_modes), mode_velocity(max_modes)
 real, allocatable :: ylm_boundary(:,:)   ! Y_k(theta_i,phi_i) for each mode k and boundary particle i

 real    :: reinject_period
 real    :: time_last_reinject  = 0.
 integer :: n_reinjections      = 0
 integer :: particles_to_inject = 0

 real    :: mass_loss_start_time, mass_loss_end_time, time_next_measurement
 real    :: mean_mass_loss_rate       = 0.
 integer :: n_measurements            = 0
 integer :: n_escaping_prev           = 0
 logical :: mass_loss_rate_calculated = .false.
 logical :: measurement_active        = .false.
 real, allocatable :: mass_loss_rates(:)

contains

subroutine set_default_options_inject(flag)
 integer, optional, intent(in) :: flag

 iboundary_spheres     = 5
 n_profile_points      = 10000
 min_mass_fraction     = 0.01
 rho_power_in          = 6.0
 r_min_on_rstar        = 1.0
 rho_inner             = 1.0e-12
 iwind                 = 1
 pulsation_period_days = 300.0
 piston_velocity_km_s  = 4.0
 phi0                  = -0.5*pi
 wss                   = 1.0
 save_period           = 0
 dumps_p_period        = 10
 reinject_enabled      = 1
 n_inject_period       = 40
 mass_loss_start       = 6.0
 mass_loss_end         = 8.0
 update_L              = 0
 verbose               = 1
 use_file_mdot         = 0
 n_modes               = 0
 l_modes               = 0
 m_modes               = 0
 mode_amplitudes       = 0.1
 mode_period_days      = 0.
 mode_phases           = 0.

end subroutine set_default_options_inject

!----------------------------------------------------------------
!+
!  Initialize everything
!+
!----------------------------------------------------------------
subroutine init_inject(ierr)
 use physcon,        only:days,au,solarm,km
 use eos,            only:gmw,gamma
 use units,          only:utime,umass,udist,unit_velocity,unit_luminosity
 use part,           only:xyzmh_ptmass,massoftype,igas,iboundary,nptmass,iTeff,iReff,iLum,npartoftype
 use injectutils,    only:get_neighb_distance,seed_random
 use wind_pulsating, only:setup_star,calc_stellar_profile,region_mass
 use dust_formation, only:calc_kappa_max
 use timestep,       only:dtmax
 integer, intent(out) :: ierr
 integer, parameter :: max_shells_tmp = 2000
 real    :: tmp_dr(max_shells_tmp), tmp_r(max_shells_tmp)
 integer :: tmp_n(max_shells_tmp)
 real    :: Mstar, Rstar, Tstar, current_radius, dr, max_nr
 integer :: nshells, n_shell, n_tot, i, k, iunit
 integer :: n_particles_first              ! derived locally — no longer a module parameter
 logical :: file_exists

 ierr = 0
 seed_random = -1  ! reset seed_random so the shell rotations are reproducible
 if (nptmass < 1) call fatal(label,'need at least one sink particle for central star')

 Mstar = xyzmh_ptmass(4,wind_emitting_sink)
 Rstar = xyzmh_ptmass(iReff,wind_emitting_sink)
 Tstar = xyzmh_ptmass(iTeff,wind_emitting_sink)
 call calc_kappa_max(Mstar*solarm, xyzmh_ptmass(iLum,wind_emitting_sink)*unit_luminosity)

 inquire(file='mass_loss_rate.dat', exist=file_exists)
 if (npartoftype(igas) < 100 .and. file_exists .and. use_file_mdot == 0) then
    print*,'Existing mass loss data file found, but this is a fresh start, so delete'
    open(newunit=iunit,file='mass_loss_rate.dat',status='old',iostat=ierr)
    close(iunit,status='delete')
    file_exists = .false.
 endif

 if (iwind == 2 .and. .not. file_exists) call calculate_period(Mstar, Rstar, pulsation_period_days)

 pulsation_period = pulsation_period_days*days/utime
 omega_pulsation  = 2.*pi/pulsation_period
 piston_velocity  = piston_velocity_km_s*km/unit_velocity
 deltaR_osc       = piston_velocity/omega_pulsation

 do k = 1,n_modes
    if (mode_period_days(k) > 0.) then
       mode_omega(k) = 2.*pi/(mode_period_days(k)*days/utime)
    else
       mode_omega(k) = omega_pulsation
    endif
    mode_deltaR(k)   = mode_amplitudes(k)*deltaR_osc
    mode_velocity(k) = mode_omega(k)*mode_deltaR(k)
 enddo
 shortest_period = pulsation_period
 do k = 1,n_modes
    if (mode_velocity(k) > 0.) shortest_period = min(shortest_period,2.*pi/mode_omega(k))
 enddo

 r_min = r_min_on_rstar*Rstar + deltaR_osc*sin(phi0)
 if (r_min <= 0.) call fatal(label,'r_min must be > 0')

 if (save_period == 1) then
    dtmax = pulsation_period/dumps_p_period
    print*,'dtmax: ',dtmax
 endif

 mass_of_gas_particle  = massoftype(igas)          ! code units (Msun when G=1, dist=au)
 massoftype(iboundary) = mass_of_gas_particle

 call derive_n_particles_first(r_min*udist, rho_inner, rho_power_in, wss, &
                               mass_of_gas_particle*umass, n_particles_first)

 current_radius = r_min
 nshells = 0
 n_tot   = 0
 do
    nshells = nshells + 1
    if (nshells > max_shells_tmp) call fatal(label,'max_shells_tmp exceeded; increase max_shells_tmp')
    if (nshells == 1) then
       n_shell = n_particles_first
    else
       n_shell = max(1, nint(tmp_n(nshells-1)*(current_radius/tmp_r(nshells-1))**(2.*(3.-rho_power_in)/3.)))
    endif
    dr = wss*current_radius*get_neighb_distance(n_shell)
    if (nshells > 2 .and. real(n_shell)/real(n_tot) < min_mass_fraction) then
       nshells = nshells - 1
       r_max   = current_radius - dr
       exit
    endif
    tmp_dr(nshells) = dr
    tmp_r(nshells)  = current_radius
    tmp_n(nshells)  = n_shell
    current_radius  = current_radius + dr
    n_tot           = n_tot + n_shell
 enddo

 call setup_star(Mstar*umass, Tstar, r_max*au, r_min*au, gmw, gamma, rho_inner, rho_power_in)
 call calc_stellar_profile(n_profile_points)

 n_shells_bnd   = min(iboundary_spheres, nshells)
 n_shells_total = nshells - n_shells_bnd
 allocate(npart_per_boundary_shell(n_shells_bnd), shell_radii_bnd(n_shells_bnd))
 allocate(npart_per_shell(n_shells_total), shell_radii_gas(n_shells_total))
 allocate(delta_r_radial(nshells))
 npart_per_boundary_shell = tmp_n(1:n_shells_bnd)
 shell_radii_bnd          = tmp_r(1:n_shells_bnd)
 npart_per_shell          = tmp_n(n_shells_bnd+1:nshells)
 shell_radii_gas          = tmp_r(n_shells_bnd+1:nshells)
 delta_r_radial           = tmp_dr(1:nshells)

 Mtotal = Mstar + sum(tmp_n(1:n_shells_total))*mass_of_gas_particle

 reinject_period      = pulsation_period/n_inject_period
 mass_loss_start_time = mass_loss_start*pulsation_period
 mass_loss_end_time   = mass_loss_end*pulsation_period
 allocate(mass_loss_rates(ceiling((mass_loss_end_time - mass_loss_start_time)/reinject_period) + 1))
 mass_loss_rates = 0.

 if (verbose == 1) then
    print *, ''
    print *, 'Calculated reinject period:', reinject_period
    print *, 'Measurement period:', reinject_period
    print *, 'Mass loss measurement start time:', mass_loss_start_time
    print *, 'Mass loss measurement end time  :', mass_loss_end_time
    print *, 'Rmax                            :', r_max
    print *, 'Expected number of measurements :', size(mass_loss_rates)
    print *, ''
 endif

 ! If use_file_mdot=1, read mass_loss_rate.dat now and mark measurement as done,
 ! so the simulation goes straight to reinjection without any measurement phase.
 ! A hard error is raised if the file is absent, since the user explicitly asked for it.
 if (use_file_mdot == 1 .and. .not. file_exists) &
    call fatal(label,'use_file_mdot=1 but mass_loss_rate.dat not found')
 if (file_exists) then
    call read_mass_loss_data()
    if (use_file_mdot == 1) mass_loss_rate_calculated = .true.
    if (use_file_mdot == 1 .and. verbose == 1) then
       print *, ''
       print *, 'use_file_mdot=1: skipping measurement phase.'
       print *, 'particles_to_inject read from file:', particles_to_inject
       print *, ''
    endif
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
       print *, 'Boundary shell ', i, ': r=', shell_radii_bnd(i)/Rstar, ' Rstar, dr=', delta_r_radial(i)/Rstar, &
                ' Rstar, N_particles=', npart_per_boundary_shell(i)
    enddo
    do i = 1, n_shells_total
       print *, 'Gas shell      ', i, ': r=', shell_radii_gas(i)/Rstar, ' Rstar, dr=', delta_r_radial(n_shells_bnd+i)/Rstar, &
                ' Rstar, N_particles=', npart_per_shell(i)
    enddo
 endif

 if (n_modes > 0) then
    call check_real_ylm()
    print *, ''
    print *, 'Non-radial pulsation modes: ', n_modes
    do k = 1,n_modes
       print "(a,i2,a,i3,a,i3,a,es11.4,a,es11.4,a,es11.4,a,es11.4,a,es11.4,a,es11.4,a)", '  mode ',k,': l=',l_modes(k), &
             ' m=',m_modes(k),' period=',2.*pi/mode_omega(k)*utime/days,' d, A_k=',mode_amplitudes(k), &
             ', dR_k=',mode_deltaR(k),' (',mode_deltaR(k)/Rstar,' Rstar), v_k=',mode_velocity(k)*unit_velocity/km, &
             ' km/s, phase=',mode_phases(k),' rad'
    enddo
    if (piston_velocity <= 0.) print *, ' WARNING: piston_velocity = 0, so all non-radial modes have zero amplitude'
    ! new shells follow the deformation, but already injected gas does not
    max_nr = max_nonradial_amplitude()
    print *, ' max over the sphere of sum_k |dR_k Y_k| :', max_nr, ' (', max_nr/Rstar, ' Rstar)'
    if (n_shells_bnd < nshells) then
       if (max_nr > delta_r_radial(n_shells_bnd+1)) then
          print *, ' WARNING: non-radial displacement exceeds the spacing above the outermost boundary shell'
          print *, '          (', delta_r_radial(n_shells_bnd+1), '), so the boundary can push into gas injected earlier'
       endif
    endif
    print *, ''
 endif

end subroutine init_inject

!----------------------------------------------------------------
!+
!  Derive n_particles_first from the target particle mass and the
!  power-law density profile at r_min.
!+
!----------------------------------------------------------------
subroutine derive_n_particles_first(r_min_cgs,rho_inner_cgs,rho_power,wss_in,m_particle_cgs,n_first)
 use injectutils, only:get_neighb_distance
 real,    intent(in)  :: r_min_cgs,rho_inner_cgs,rho_power,wss_in,m_particle_cgs
 integer, intent(out) :: n_first
 integer, parameter :: max_iter = 100
 integer :: iter,n_old
 real    :: C_rho,exponent,dr_cgs

 C_rho    = rho_inner_cgs*r_min_cgs**rho_power
 exponent = 3. - rho_power
 n_old    = 1000
 do iter = 1,max_iter
    dr_cgs  = wss_in*r_min_cgs*get_neighb_distance(n_old)
    n_first = max(1, nint(4.*pi*C_rho*((r_min_cgs + dr_cgs)**exponent - r_min_cgs**exponent) &
                          /exponent/m_particle_cgs))
    if (n_first == n_old) exit
    n_old = (n_old + n_first)/2
 enddo

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
!  The actual function that is called by phantom
!+
!----------------------------------------------------------------
subroutine inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npart_old,npartoftype,dtinject)
 real,    intent(in)    :: time,dtlast
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(inout) :: npart,npart_old
 integer, intent(inout) :: npartoftype(:)
 real,    intent(out)   :: dtinject

 dtinject = pulsation_timestep*shortest_period

 if (.not. atmosphere_setup_complete) then
    atmosphere_setup_complete = .true.
    if (npart == 0) then
       print *, ''
       print *, 'Setting up stellar atmosphere with ', n_shells_total, ' shells.'
       call setup_initial_atmosphere(time,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npartoftype)
       print *, 'Stellar atmosphere setup complete.'
       return
    endif
 endif

 if (.not. boundary_info_ready) then
    call reconstruct_boundary_info(time,xyzh,xyzmh_ptmass)
    time_last_reinject = time - mod(time,reinject_period)
 endif

 if (reinject_enabled == 1) then
    if (.not. mass_loss_rate_calculated) &
       call take_periodic_mass_measurements(time,xyzh,vxyzu,npart,xyzmh_ptmass,vxyz_ptmass)
    if (mass_loss_rate_calculated .and. time - time_last_reinject >= reinject_period &
        .and. time >= mass_loss_start_time) then
       call perform_reinjection(time,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npartoftype)
       time_last_reinject = time
    endif
 endif

 call apply_pulsation(time,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass)

end subroutine inject_particles

!----------------------------------------------------------------
!+
!  Counts particles that are unbound from the AGB star
!+
!----------------------------------------------------------------
function count_unbound(xyzh,vxyzu,npart,xyzmh_ptmass,vxyz_ptmass) result(n)
 real,    intent(in) :: xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(in) :: npart
 integer :: n,i
 real    :: x0(3),v0(3),m,r

 x0 = xyzmh_ptmass(1:3,wind_emitting_sink)
 v0 = vxyz_ptmass(1:3,wind_emitting_sink)
 m  = xyzmh_ptmass(4,wind_emitting_sink)
 n  = 0
 do i = 1,npart
    r = norm2(xyzh(1:3,i) - x0)
    if (r <= 0.) cycle
    if (0.5*sum((vxyzu(1:3,i) - v0)**2) + vxyzu(4,i) - m/r > 0.) n = n + 1
 enddo

end function count_unbound

!----------------------------------------------------------------
!+
!  Checks how much mass the star has lost, and calculates the mass loss rate
!+
!----------------------------------------------------------------
subroutine take_periodic_mass_measurements(time,xyzh,vxyzu,npart,xyzmh_ptmass,vxyz_ptmass)
 real,    intent(in) :: time,xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(in) :: npart
 integer :: n_escaping,newly_unbound

 if (.not. measurement_active .and. time >= mass_loss_start_time) then
    measurement_active    = .true.
    time_next_measurement = time + reinject_period
    n_escaping_prev       = count_unbound(xyzh,vxyzu,npart,xyzmh_ptmass,vxyz_ptmass)
    if (verbose == 1) print*,' Baseline unbound particles at measurement start:',n_escaping_prev
 endif

 if (measurement_active .and. time >= time_next_measurement .and. time < mass_loss_end_time) then
    n_escaping      = count_unbound(xyzh,vxyzu,npart,xyzmh_ptmass,vxyz_ptmass)
    newly_unbound   = max(0, n_escaping - n_escaping_prev)
    n_escaping_prev = n_escaping
    n_measurements  = n_measurements + 1
    mass_loss_rates(n_measurements) = newly_unbound*mass_of_gas_particle/reinject_period
    time_next_measurement = time + reinject_period
    if (verbose == 1) then
       print *, ' Unbound+outflowing particles     :',n_escaping
       print*,' Newly unbound since last snapshot:',newly_unbound
       print*,' Rate this interval               :',mass_loss_rates(n_measurements)
    endif
 endif

 if (measurement_active .and. time >= mass_loss_end_time .and. n_measurements > 0) then
    mean_mass_loss_rate       = sum(mass_loss_rates(1:n_measurements))/n_measurements
    particles_to_inject       = max(1, nint(mean_mass_loss_rate*reinject_period/mass_of_gas_particle))
    mass_loss_rate_calculated = .true.
    call write_mass_loss_data()
 endif

end subroutine take_periodic_mass_measurements

!----------------------------------------------------------------
!+
!  Inject particles throughout the simulation
!+
!----------------------------------------------------------------
subroutine perform_reinjection(time,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npartoftype)
 use part,           only:igas
 use injectutils,    only:inject_geodesic_sphere
 use wind_pulsating, only:interp_stellar_profile
 real,    intent(in)    :: time
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(inout) :: npart,npartoftype(:)
 integer :: old_npart
 real    :: r_inject,phase,rho,u,T,P,x0(3),v0(3)

 x0    = xyzmh_ptmass(1:3,wind_emitting_sink)
 v0    = vxyz_ptmass(1:3,wind_emitting_sink)
 phase = omega_pulsation*time + phi0

 if (n_boundary_particles > 0) then
    r_inject = r_boundary_equilibrium(n_boundary_particles) + delta_r_radial(iboundary_spheres+1)
 else
    r_inject = r_min + sum(delta_r_radial(1:iboundary_spheres+1))
 endif
 r_inject = r_inject + deltaR_osc*sin(phase)
 call interp_stellar_profile(r_inject,rho,P,u,T)

 old_npart      = npart
 n_reinjections = n_reinjections + 1
 call inject_geodesic_sphere(n_shells_total + n_reinjections, npart+1, particles_to_inject, r_inject, &
                             piston_velocity*cos(phase), u, rho, npart, npartoftype, xyzh, vxyzu, igas, x0, v0, &
                             wind_emitting_sink)
 ! follow the boundary deformation below the new shell, in position and velocity
 if (n_modes > 0) call displace_nonradial(time,x0,old_npart+1,npart,xyzh,vxyzu,.true.)

 xyzmh_ptmass(4,wind_emitting_sink) = xyzmh_ptmass(4,wind_emitting_sink) - (npart - old_npart)*mass_of_gas_particle

 if (verbose == 1) then
    print *, ''
    print *, ' Injection time             :', time
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
subroutine setup_initial_atmosphere(time,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npartoftype)
 use part,           only:igas,iboundary
 use injectutils,    only:inject_geodesic_sphere
 use wind_pulsating, only:interp_stellar_profile
 real,    intent(in)    :: time,vxyz_ptmass(:,:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:)
 integer, intent(inout) :: npart,npartoftype(:)
 integer :: i
 real    :: rho,u,T,P,x0(3),v0(3)

 x0 = xyzmh_ptmass(1:3,wind_emitting_sink)
 v0 = vxyz_ptmass(1:3,wind_emitting_sink)

 npart = 0
 do i = 1,n_shells_bnd
    call interp_stellar_profile(shell_radii_bnd(i),rho,P,u,T)
    call inject_geodesic_sphere(i, npart+1, npart_per_boundary_shell(i), shell_radii_bnd(i), 0., u, rho, &
                                npart, npartoftype, xyzh, vxyzu, iboundary, x0, v0, wind_emitting_sink)
 enddo
 do i = 1,n_shells_total
    call interp_stellar_profile(shell_radii_gas(i),rho,P,u,T)
    call inject_geodesic_sphere(n_shells_bnd+i, npart+1, npart_per_shell(i), shell_radii_gas(i), 0., u, rho, &
                                npart, npartoftype, xyzh, vxyzu, igas, x0, v0, wind_emitting_sink)
 enddo

 n_boundary_particles = npartoftype(iboundary)
 allocate(r_boundary_equilibrium(n_boundary_particles))
 do i = 1,n_boundary_particles
    r_boundary_equilibrium(i) = norm2(xyzh(1:3,i) - x0) - deltaR_osc*sin(phi0)
 enddo
 boundary_info_ready = .true.

 ! shells are built at the radial-only position; add the non-radial displacement
 ! and set consistent velocities so the first pulsation step does not jump
 if (n_modes > 0) then
    call compute_boundary_ylm(xyzh,x0)
    call apply_pulsation(time,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass)
    call print_boundary_displacement(time)

    ! the gas shells start on the deformed surface; gas velocities are left at zero,
    ! as for the radial term
    call displace_nonradial(time,x0,n_boundary_particles+1,npart,xyzh,vxyzu,.false.)
    print*,'Gas shells displaced to follow the initial boundary deformation'
 endif

 print*,'Boundary particles : ',n_boundary_particles
 print*,'Gas particles      : ',npartoftype(igas)

end subroutine setup_initial_atmosphere

!----------------------------------------------------------------
!+
!  Reconstructs boundary particle info after resuming from a dump
!+
!----------------------------------------------------------------
subroutine reconstruct_boundary_info(time,xyzh,xyzmh_ptmass)
 use part, only:iboundary,npartoftype
 real, intent(in) :: time,xyzh(:,:),xyzmh_ptmass(:,:)
 integer :: i
 real    :: x0(3),osc,dr_nr,sterm(max_modes),cterm(max_modes)

 x0  = xyzmh_ptmass(1:3,wind_emitting_sink)
 osc = deltaR_osc*sin(omega_pulsation*time + phi0)

 n_boundary_particles = npartoftype(iboundary)
 allocate(r_boundary_equilibrium(n_boundary_particles))
 ! boundary particles only move radially, so Y_k can be recovered from their current direction
 if (n_modes > 0) call compute_boundary_ylm(xyzh,x0)
 call get_mode_terms(time,sterm,cterm)
 ! subtract the full displacement at the restart time (radial + all modes)
 do i = 1,n_boundary_particles
    dr_nr = 0.
    if (n_modes > 0) dr_nr = dot_product(ylm_boundary(1:n_modes,i),sterm(1:n_modes))
    r_boundary_equilibrium(i) = norm2(xyzh(1:3,i) - x0) - osc - dr_nr
 enddo
 boundary_info_ready = .true.

 if (n_boundary_particles > 0) then
    print *, 'Reconstructed boundary particle info:'
    print *, 'Boundary particles: ', n_boundary_particles
 endif
 if (n_modes > 0) call print_boundary_displacement(time)

end subroutine reconstruct_boundary_info

!----------------------------------------------------------------
!+
!  Applies the pulsation to the boundary layers
!+
!----------------------------------------------------------------
subroutine apply_pulsation(time,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass)
 use part, only:iTeff,iLum,iReff
 real, intent(in)    :: time,vxyz_ptmass(:,:)
 real, intent(inout) :: xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:)
 integer :: i
 real    :: phase,osc,r_dot,x0(3),v0(3),x_hat(3)
 real    :: dr_nr,dv_nr,sterm(max_modes),cterm(max_modes)

 if (n_boundary_particles == 0) return

 x0    = xyzmh_ptmass(1:3,wind_emitting_sink)
 v0    = vxyz_ptmass(1:3,wind_emitting_sink)
 phase = omega_pulsation*time + phi0
 osc   = deltaR_osc*sin(phase)
 r_dot = piston_velocity*cos(phase)
 call get_mode_terms(time,sterm,cterm)

 do i = 1,n_boundary_particles
    ! Radial unit vector from current position
    x_hat        = xyzh(1:3,i) - x0
    x_hat        = x_hat/norm2(x_hat)
    dr_nr        = 0.
    dv_nr        = 0.
    if (n_modes > 0) then
       dr_nr = dot_product(ylm_boundary(1:n_modes,i),sterm(1:n_modes))
       dv_nr = dot_product(ylm_boundary(1:n_modes,i),cterm(1:n_modes))
    endif
    xyzh(1:3,i)  = x0 + (r_boundary_equilibrium(i) + osc + dr_nr)*x_hat
    vxyzu(1:3,i) = v0 + (r_dot + dv_nr)*x_hat
 enddo

 ! the luminosity follows the radial term only: l >= 1 modes integrate to zero
 ! over the sphere, so they do not change the mean radius to first order
 ! (an l = 0 mode would, but it is not included here either)
 if (update_L == 1) call get_lum(xyzmh_ptmass(iLum,wind_emitting_sink), xyzmh_ptmass(iTeff,wind_emitting_sink), &
                                 xyzmh_ptmass(iReff,wind_emitting_sink) + osc)

end subroutine apply_pulsation

!----------------------------------------------------------------
!+
!  Per-mode time factors: sterm = dR_k sin(w_k t + phi_k) and
!  cterm = w_k dR_k cos(w_k t + phi_k), to be multiplied by Y_k
!+
!----------------------------------------------------------------
subroutine get_mode_terms(time,sterm,cterm)
 real, intent(in)  :: time
 real, intent(out) :: sterm(max_modes),cterm(max_modes)
 integer :: k
 real    :: phase_k

 sterm = 0.
 cterm = 0.
 do k = 1,n_modes
    phase_k  = mode_omega(k)*time + mode_phases(k)
    sterm(k) = mode_deltaR(k)*sin(phase_k)
    cterm(k) = mode_velocity(k)*cos(phase_k)
 enddo

end subroutine get_mode_terms

!----------------------------------------------------------------
!+
!  Store Y_k(theta_i,phi_i) for every boundary particle and mode.
!  Boundary particles only move radially, so this is done once.
!+
!----------------------------------------------------------------
subroutine compute_boundary_ylm(xyzh,x0)
 real, intent(in) :: xyzh(:,:),x0(3)
 integer :: i,k
 real    :: x_hat(3)

 if (allocated(ylm_boundary)) deallocate(ylm_boundary)
 allocate(ylm_boundary(n_modes,n_boundary_particles))
 do i = 1,n_boundary_particles
    x_hat = xyzh(1:3,i) - x0
    x_hat = x_hat/norm2(x_hat)
    do k = 1,n_modes
       ylm_boundary(k,i) = real_ylm(l_modes(k),m_modes(k),x_hat)
    enddo
 enddo

end subroutine compute_boundary_ylm

!----------------------------------------------------------------
!+
!  Shift gas particles i1..i2 radially by the non-radial displacement
!  of the boundary in their direction, so a shell keeps a constant
!  spacing above the deformed surface (u and h are unchanged).
!  If add_velocity, also add the matching non-radial radial velocity.
!+
!----------------------------------------------------------------
subroutine displace_nonradial(time,x0,i1,i2,xyzh,vxyzu,add_velocity)
 real,    intent(in)    :: time,x0(3)
 integer, intent(in)    :: i1,i2
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 logical, intent(in)    :: add_velocity
 integer :: i,k
 real    :: x_hat(3),r,ylm,dr_nr,dv_nr,sterm(max_modes),cterm(max_modes)

 call get_mode_terms(time,sterm,cterm)
 do i = i1,i2
    x_hat = xyzh(1:3,i) - x0
    r     = norm2(x_hat)
    x_hat = x_hat/r
    dr_nr = 0.
    dv_nr = 0.
    do k = 1,n_modes
       ylm   = real_ylm(l_modes(k),m_modes(k),x_hat)
       dr_nr = dr_nr + ylm*sterm(k)
       dv_nr = dv_nr + ylm*cterm(k)
    enddo
    xyzh(1:3,i) = x0 + (r + dr_nr)*x_hat
    if (add_velocity) vxyzu(1:3,i) = vxyzu(1:3,i) + dv_nr*x_hat
 enddo

end subroutine displace_nonradial

!----------------------------------------------------------------
!+
!  Print the actual non-radial displacement of the boundary particles
!  (the Y_lm coefficient alone is not the peak displacement)
!+
!----------------------------------------------------------------
subroutine print_boundary_displacement(time)
 real, intent(in) :: time
 integer :: i,k
 real    :: sterm(max_modes),cterm(max_modes),dr_nr,dmin,dmax,dabs_max

 if (n_boundary_particles == 0) return
 call get_mode_terms(time,sterm,cterm)
 print *, 'Non-radial displacement of the ', n_boundary_particles, ' boundary particles:'
 do k = 1,n_modes
    print "(a,i2,a,2(1x,es11.4))", '  mode ',k,': min/max of dR_k*Y_k              :', &
          minval(mode_deltaR(k)*ylm_boundary(k,:)),maxval(mode_deltaR(k)*ylm_boundary(k,:))
 enddo
 dmin     = huge(dmin)
 dmax     = -huge(dmax)
 dabs_max = 0.
 do i = 1,n_boundary_particles
    dr_nr    = dot_product(ylm_boundary(1:n_modes,i),sterm(1:n_modes))
    dmin     = min(dmin,dr_nr)
    dmax     = max(dmax,dr_nr)
    dabs_max = max(dabs_max,sum(abs(mode_deltaR(1:n_modes)*ylm_boundary(1:n_modes,i))))
 enddo
 print "(a,2(1x,es11.4))", '  total displacement now, min/max          :',dmin,dmax
 print "(a,1x,es11.4)",    '  total upper bound, max_i sum_k|dR_k*Y_k| :',dabs_max

end subroutine print_boundary_displacement

!----------------------------------------------------------------
!+
!  Orthonormal real spherical harmonic Y_lm for a unit vector x_hat,
!  with theta, phi measured about the z-axis:
!
!    Y_l0 =         N_l0  P_l(cos theta)
!    Y_lm = sqrt(2) N_lm  P_l^m(cos theta) cos(m phi)      (m > 0)
!    Y_lm = sqrt(2) N_l|m| P_l^|m|(cos theta) sin(|m| phi) (m < 0)
!
!  with N_lm = sqrt((2l+1)/(4 pi) (l-m)!/(l+m)!) and P_l^m WITHOUT the
!  Condon-Shortley phase. legendre_associated includes (-1)^m, which is
!  removed here, so that Y_11 ~ +x, Y_1-1 ~ +y and Y_10 ~ +z.
!+
!----------------------------------------------------------------
real function real_ylm(l,m,x_hat)
 use mathfunc, only:legendre_associated
 integer, intent(in) :: l,m
 real,    intent(in) :: x_hat(3)
 integer :: i,ma
 real    :: plm(0:l),norm,phi

 ma = abs(m)
 call legendre_associated(l,ma,max(-1.,min(1.,x_hat(3))),plm)
 norm = (2.*l + 1.)/(4.*pi)
 do i = l-ma+1,l+ma
    norm = norm/i
 enddo
 norm = sqrt(norm)*(-1)**ma
 if (m == 0) then
    real_ylm = norm*plm(l)
 else
    phi = atan2(x_hat(2),x_hat(1))
    if (m > 0) then
       real_ylm = sqrt(2.)*norm*plm(l)*cos(ma*phi)
    else
       real_ylm = sqrt(2.)*norm*plm(l)*sin(ma*phi)
    endif
 endif

end function real_ylm

!----------------------------------------------------------------
!+
!  Check real_ylm against closed forms for l <= 2
!+
!----------------------------------------------------------------
subroutine check_real_ylm()
 integer, parameter :: ndir = 4
 real    :: dirs(3,ndir),x,y,z,c1,c2,expected(9),got(9)
 integer :: idir,l,m,j

 dirs(:,1) = (/0.3,-0.5,0.8/)
 dirs(:,2) = (/-0.6,0.2,-0.1/)
 dirs(:,3) = (/0.,0.,1./)
 dirs(:,4) = (/-0.2,-0.7,0./)
 c1 = sqrt(3./(4.*pi))
 c2 = sqrt(15./pi)
 do idir = 1,ndir
    dirs(:,idir) = dirs(:,idir)/norm2(dirs(:,idir))
    x = dirs(1,idir)
    y = dirs(2,idir)
    z = dirs(3,idir)
    ! ordered as (l,m) = (0,0),(1,-1),(1,0),(1,1),(2,-2),(2,-1),(2,0),(2,1),(2,2)
    expected = (/0.5/sqrt(pi), c1*y, c1*z, c1*x, 0.5*c2*x*y, 0.5*c2*y*z, &
                 0.25*sqrt(5./pi)*(3.*z**2 - 1.), 0.5*c2*x*z, 0.25*c2*(x**2 - y**2)/)
    j = 0
    do l = 0,2
       do m = -l,l
          j = j + 1
          got(j) = real_ylm(l,m,dirs(:,idir))
       enddo
    enddo
    if (any(abs(got - expected) > 1.e-12)) then
       print *, 'Y_lm check failed: got     ', got
       print *, '                   expected', expected
       call fatal(label,'real spherical harmonics do not match closed forms')
    endif
 enddo

end subroutine check_real_ylm

!----------------------------------------------------------------
!+
!  Worst-case non-radial displacement, max over the sphere of
!  sum_k |dR_k Y_k|, sampled on a theta-phi grid
!+
!----------------------------------------------------------------
real function max_nonradial_amplitude()
 integer, parameter :: ntheta = 181, nphi = 360
 integer :: i,j,k
 real    :: theta,phi,x_hat(3),amp

 max_nonradial_amplitude = 0.
 do i = 1,ntheta
    theta = pi*(i-1)/(ntheta-1)
    do j = 1,nphi
       phi   = 2.*pi*(j-1)/nphi
       x_hat = (/sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)/)
       amp   = 0.
       do k = 1,n_modes
          amp = amp + abs(mode_deltaR(k)*real_ylm(l_modes(k),m_modes(k),x_hat))
       enddo
       max_nonradial_amplitude = max(max_nonradial_amplitude,amp)
    enddo
 enddo

end function max_nonradial_amplitude

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
 use physcon, only:au,steboltz
 use units,   only:unit_luminosity
 real, intent(out) :: Lum
 real, intent(in)  :: Teff,Reff

 Lum = 4.*pi*steboltz*Teff**4*(Reff*au)**2/unit_luminosity

end subroutine get_lum

!----------------------------------------------------------------
!+
!  Write mass-loss information to file for resume from dump
!+
!----------------------------------------------------------------
subroutine write_mass_loss_data()
 use io,      only:iprint
 use physcon, only:solarm,years
 use units,   only:umass,utime
 integer :: iunit,ierr,i

 open(newunit=iunit,file='mass_loss_rate.dat',status='replace',iostat=ierr)
 if (ierr /= 0) then
    write(iprint,*) 'Could not write mass_loss_rate.dat'
    return
 endif

 write(iunit,*) '# Mass-loss rate data for restart'
 write(iunit,*) r_max
 write(iunit,*) r_max
 write(iunit,*) mass_loss_rate_calculated
 write(iunit,*) mean_mass_loss_rate/(solarm/umass)/(utime/years)
 write(iunit,*) Mtotal
 write(iunit,*) particles_to_inject
 write(iunit,*) n_measurements
 write(iunit,*) mass_of_gas_particle
 do i = 1,n_measurements
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
 integer :: iunit,ierr,i

 open(newunit=iunit,file='mass_loss_rate.dat',status='old',iostat=ierr)
 if (ierr /= 0) return

 read(iunit,*)
 read(iunit,*,iostat=ierr) r_max
 read(iunit,*,iostat=ierr) r_max
 read(iunit,*,iostat=ierr) mass_loss_rate_calculated
 if (ierr /= 0) then
    close(iunit)
    return
 endif
 read(iunit,*,iostat=ierr) mean_mass_loss_rate
 read(iunit,*,iostat=ierr) Mtotal
 read(iunit,*,iostat=ierr) particles_to_inject
 read(iunit,*,iostat=ierr) n_measurements
 read(iunit,*,iostat=ierr) mass_of_gas_particle

 if (n_measurements > 0) then
    if (.not. allocated(mass_loss_rates)) allocate(mass_loss_rates(n_measurements))
    do i = 1,n_measurements
       read(iunit,*,iostat=ierr) mass_loss_rates(i)
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
subroutine calculate_period(M,R,pulsation_period_days)
 real, intent(in)  :: M,R
 real, intent(out) :: pulsation_period_days

 pulsation_period_days = 10.**(-1.92 - 0.73*log10(M) + 1.86*log10(R*215.032))
 print*,'Calculated pulsation period (days): ',pulsation_period_days

end subroutine calculate_period

!----------------------------------------------------------------
!+
!  Write options to .in file
!+
!----------------------------------------------------------------
subroutine write_options_inject(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(n_profile_points,      'n_profile_points',  'number of points in stellar profile',iunit)
 call write_inopt(iboundary_spheres,     'iboundary_spheres', 'number of boundary spheres (piston layers)',iunit)
 call write_inopt(min_mass_fraction,     'min_mass_fraction', 'minimum mass fraction per shell',iunit)
 call write_inopt(rho_power_in,          'rho_power',         'density profile exponent: rho ~ r^(-rho_power)',iunit)
 call write_inopt(r_min_on_rstar,        'r_min_on_rstar',    'gas atmosphere inner radius as fraction of R_star',iunit)
 call write_inopt(rho_inner,             'rho_inner',         'inner boundary density at r_min (cgs)',iunit)
 call write_inopt(iwind,                 'iwind',             'wind type: 1=prescribed, 2=period from mass-radius relation',iunit)
 call write_inopt(pulsation_period_days, 'pulsation_period',  'pulsation period (days)',iunit)
 call write_inopt(piston_velocity_km_s,  'piston_velocity',   'piston velocity amplitude (km/s)',iunit)
 call write_inopt(phi0,                  'phi0',              'initial phase offset (radians)',iunit)
 call write_inopt(wss,                   'wss',               'radial/tangential spacing ratio',iunit)
 call write_inopt(save_period,           'save_period',       'wether to save dumps as fraction of period (0=off, 1=on)',iunit)
 call write_inopt(dumps_p_period,        'dumps_p_period',    'number of dumps per period (if save_period = 1)',iunit)
 call write_inopt(reinject_enabled,      'reinject_enabled',  'enable dynamic reinjection (0=off, 1=on)',iunit)
 call write_inopt(n_inject_period,       'n_inject_period',   'period between reinjections (periods)',iunit)
 call write_inopt(mass_loss_start,       'mass_loss_start',   'start time for mass-loss calculation (periods)',iunit)
 call write_inopt(mass_loss_end,         'mass_loss_end',     'end time for mass-loss calculation (periods)',iunit)
 call write_inopt(use_file_mdot,         'use_file_mdot',     'skip measurement phase (0=off, 1=on)',iunit)
 call write_inopt(update_L,              'update_L',          'update luminosity with pulsation (0=off, 1=on)',iunit)
 call write_inopt(verbose,               'verbose',           'enable verbose output (0=off, 1=on)',iunit)
 call write_inopt(n_modes,               'n_modes',           'number of non-radial modes (0=radial pulsation only)',iunit)
 if (n_modes > 0) then
    call write_mode_list(int_list(l_modes(1:n_modes)),             'l_modes',        &
                         'spherical harmonic degree l per mode',iunit)
    call write_mode_list(int_list(m_modes(1:n_modes)),             'm_modes',        &
                         'spherical harmonic order m per mode',iunit)
    call write_mode_list(real_list(mode_amplitudes(1:n_modes)),    'mode_amplitudes', &
                         'amplitude per mode, relative to the radial amplitude',iunit)
    call write_mode_list(real_list(mode_period_days(1:n_modes)),   'mode_periods',   &
                         'period per mode (days), <=0 uses the radial period',iunit)
    call write_mode_list(real_list(mode_phases(1:n_modes)),        'mode_phases',    &
                         'phase offset per mode (radians)',iunit)
 endif

end subroutine write_options_inject

!----------------------------------------------------------------
!+
!  Write a comma-separated list of mode values to the .in file
!+
!----------------------------------------------------------------
subroutine write_mode_list(str,name,descript,iunit)
 use infile_utils, only:write_inopt
 use io,           only:error
 character(len=*), intent(in) :: str,name,descript
 integer,          intent(in) :: iunit

 ! the .in reader keeps 120 characters per line and the value starts at column 25
 if (len(str) > maxlen_modelist) call error(label,'list of '//name//' too long for the .in file: '//str)
 call write_inopt(str,name,descript,iunit)

end subroutine write_mode_list

!----------------------------------------------------------------
!+
!  Comma-separated list of integers
!+
!----------------------------------------------------------------
function int_list(arr) result(str)
 integer, intent(in) :: arr(:)
 character(len=:), allocatable :: str
 character(len=12) :: tmp
 integer :: i

 str = ''
 do i = 1,size(arr)
    write(tmp,'(i0)') arr(i)
    if (i > 1) str = str//','
    str = str//trim(tmp)
 enddo

end function int_list

!----------------------------------------------------------------
!+
!  Comma-separated list of reals with 6 significant digits and
!  trailing zeros stripped (stable under write -> read -> write);
!  exponent notation only outside 1e-3 <= |x| < 1e6
!+
!----------------------------------------------------------------
function real_list(arr) result(str)
 real, intent(in) :: arr(:)
 character(len=:), allocatable :: str
 character(len=24) :: tmp,expo
 integer :: i,j,iexp,iexpo

 str = ''
 do i = 1,size(arr)
    write(tmp,'(es12.5e2)') arr(i)
    tmp   = adjustl(tmp)
    iexpo = index(tmp,'E')
    read(tmp(iexpo+1:),*) iexp
    if (arr(i) == 0. .or. (iexp >= -3 .and. iexp < 6)) then
       ! fixed-point notation with the same 6 significant digits
       write(tmp,'(f24.'//achar(iachar('0')+max(1,5-iexp))//')') arr(i)
       tmp  = adjustl(tmp)
       expo = ''
    else
       expo = tmp(iexpo:)
       tmp  = tmp(1:iexpo-1)
    endif
    j = len_trim(tmp)
    do while (tmp(j:j) == '0' .and. tmp(j-1:j-1) /= '.')
       tmp(j:j) = ' '
       j = j - 1
    enddo
    if (i > 1) str = str//','
    str = str//trim(tmp)//trim(expo)
 enddo

end function real_list

!----------------------------------------------------------------
!+
!  Read options from .in file
!+
!----------------------------------------------------------------
subroutine read_options_inject(db,nerr)
 use infile_utils, only:inopts,read_inopt
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr
 integer :: k

 call read_inopt(n_profile_points,      'n_profile_points',  db,errcount=nerr,min=11)
 call read_inopt(iboundary_spheres,     'iboundary_spheres', db,errcount=nerr,min=0)
 call read_inopt(min_mass_fraction,     'min_mass_fraction', db,errcount=nerr,min=0.)
 call read_inopt(rho_power_in,          'rho_power',         db,errcount=nerr,min=epsilon(rho_power_in))
 call read_inopt(r_min_on_rstar,        'r_min_on_rstar',    db,errcount=nerr,min=epsilon(r_min_on_rstar),max=2.)
 if (r_min_on_rstar >= 2.) call fatal(label,'r_min_on_rstar must be in (0,2)')
 call read_inopt(rho_inner,             'rho_inner',         db,errcount=nerr,min=epsilon(rho_inner))
 call read_inopt(iwind,                 'iwind',             db,errcount=nerr,min=1,max=2)
 call read_inopt(pulsation_period_days, 'pulsation_period',  db,errcount=nerr,min=0.)
 call read_inopt(piston_velocity_km_s,  'piston_velocity',   db,errcount=nerr,min=0.)
 call read_inopt(phi0,                  'phi0',              db,errcount=nerr,min=-pi,max=pi)
 call read_inopt(wss,                   'wss',               db,errcount=nerr,min=epsilon(wss),max=10.)
 call read_inopt(save_period,           'save_period',       db,errcount=nerr,min=0,max=1)
 call read_inopt(dumps_p_period,        'dumps_p_period',    db,errcount=nerr,min=0)
 call read_inopt(reinject_enabled,      'reinject_enabled',  db,errcount=nerr,min=0,max=1)
 call read_inopt(n_inject_period,       'n_inject_period',   db,errcount=nerr,min=1)
 call read_inopt(mass_loss_start,       'mass_loss_start',   db,errcount=nerr,min=0.)
 call read_inopt(mass_loss_end,         'mass_loss_end',     db,errcount=nerr,min=epsilon(mass_loss_end))
 call read_inopt(use_file_mdot,         'use_file_mdot',     db,errcount=nerr,min=0,max=1)
 call read_inopt(update_L,              'update_L',          db,errcount=nerr,min=0,max=1)
 call read_inopt(verbose,               'verbose',           db,errcount=nerr,min=0,max=1)
 call read_inopt(n_modes,               'n_modes',           db,errcount=nerr,min=0,max=max_modes,default=0)
 if (n_modes > 0) then
    call read_mode_list('l_modes',        db,nerr,ival=l_modes)
    call read_mode_list('m_modes',        db,nerr,ival=m_modes)
    call read_mode_list('mode_amplitudes',db,nerr,rval=mode_amplitudes)
    call read_mode_list('mode_periods',   db,nerr,rval=mode_period_days)
    call read_mode_list('mode_phases',    db,nerr,rval=mode_phases)
    do k = 1,n_modes
       if (l_modes(k) < 0 .or. abs(m_modes(k)) > l_modes(k)) then
          print "(a,i2,a,i3,a,i3,a)", ' ERROR: mode ',k,' has l=',l_modes(k),', m=',m_modes(k),'; need l >= 0 and |m| <= l'
          nerr = nerr + 1
       endif
       if (mode_amplitudes(k) < 0.) then
          print "(a,i2,a)", ' ERROR: mode ',k,' has a negative amplitude'
          nerr = nerr + 1
       endif
    enddo
 endif

end subroutine read_options_inject

!----------------------------------------------------------------
!+
!  Read a comma-separated list of n_modes integers (ival) or reals (rval)
!+
!----------------------------------------------------------------
subroutine read_mode_list(tag,db,nerr,ival,rval)
 use infile_utils, only:inopts,read_inopt
 character(len=*),  intent(in)    :: tag
 type(inopts),      intent(inout) :: db(:)
 integer,           intent(inout) :: nerr
 integer, optional, intent(inout) :: ival(:)
 real,    optional, intent(inout) :: rval(:)
 character(len=120) :: valstring
 integer :: i,nval,ierr

 call read_inopt(valstring,tag,db,err=ierr)
 if (ierr /= 0) then
    print "(a)", ' ERROR: '//tag//' not found (required when n_modes > 0)'
    nerr = nerr + 1
    return
 endif
 nval = 0
 if (len_trim(valstring) > 0) nval = 1
 do i = 1,len_trim(valstring)
    if (valstring(i:i) == ',') nval = nval + 1
 enddo
 if (nval /= n_modes) then
    print "(a,i2,a,i2)", ' ERROR: '//tag//' = '//trim(valstring)//' has ',nval,' entries, but n_modes = ',n_modes
    nerr = nerr + 1
    return
 endif
 if (present(ival)) read(valstring,*,iostat=ierr) ival(1:n_modes)
 if (present(rval)) read(valstring,*,iostat=ierr) rval(1:n_modes)
 if (ierr /= 0) then
    print "(a)", ' ERROR: could not read '//tag//' = '//trim(valstring)
    nerr = nerr + 1
 endif

end subroutine read_mode_list

end module inject
