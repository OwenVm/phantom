!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module wind_pulsating
!
! driver to integrate the hydrostatic equilibrium equations to set the outer layers of an AGB star
!
! :References: None
!
! :Owner: Owen Vermeulen
!
! :Runtime parameters: None
!
! :Dependencies: eos, io, physcon, table_utils, units
!
 implicit none
 public :: setup_star
 public :: stellar_state,save_stellarprofile,interp_stellar_profile,calc_stellar_profile
 public :: region_mass

 private
 ! Density profile exponent: rho ~ r^(-rho_power)
 ! This is now a runtime parameter set via setup_star.
 ! Typical values: 2.0 (default, analytic), 10.0 (grid-code steep profile)
 real :: rho_power = 2.0

 ! input parameters
 real :: Mstar_cgs, Rstar_cgs, r_inner, Star_gamma, Star_mu, P0, rho0, Menv_cgs
 real, dimension(:,:), allocatable, public :: stellar_1D

 ! stellar properties
 type stellar_state
    real :: r, r0, Rstar, rho, P, u, T
    integer :: nsteps
    logical :: error
 end type stellar_state

contains

subroutine setup_star(Mstar_in, Rstar_in, r_min, mu_in, gamma_in, surface_pressure, M_env_in, rho_power_in)
 use physcon, only:au, solarm

 real, intent(in)           :: Mstar_in, Rstar_in, r_min, mu_in, gamma_in, surface_pressure, M_env_in
 real, intent(in), optional :: rho_power_in

 Mstar_cgs  = Mstar_in
 Rstar_cgs  = Rstar_in
 r_inner    = r_min
 Star_gamma = gamma_in
 Star_mu    = mu_in
 P0         = surface_pressure
 Menv_cgs   = M_env_in

 ! Set density profile exponent; default to 2 if not supplied
 if (present(rho_power_in)) then
    rho_power = rho_power_in
 else
    rho_power = 2.0
 endif

 print *, "Setting up star with parameters:"
 print *, "Rstar_cgs  :", Rstar_cgs
 print *, "Mstar_cgs  :", Mstar_cgs
 print *, "Menv_cgs   :", Menv_cgs
 print *, "r_inner    :", r_inner
 print *, "rho_power  :", rho_power

end subroutine setup_star

!-----------------------------------------------------------------------
!
!  Normalization constant for the power-law density profile
!
!  rho(r) = C_rho / r^rho_power
!
!-----------------------------------------------------------------------
real function calc_C_rho()
 use physcon, only:pi

 real :: exponent, integral

 exponent = 3.0 - rho_power

 if (abs(exponent) > 1.e-6) then
    integral = (Rstar_cgs**exponent - r_inner**exponent) / exponent
 else
    ! rho_power == 3: integral of 1/r from r_inner to Rstar
    integral = log(Rstar_cgs / r_inner)
 endif

 calc_C_rho = Menv_cgs / (4.0 * pi * integral)

end function calc_C_rho

!-----------------------------------------------------------------------
!
!  Enclosed envelope mass from r_inner to r (used in hydrostatic step)
!
!  M_env(r) = 4*pi * C_rho * int_{r_inner}^{r} r'^(2-rho_power) dr'
!
!-----------------------------------------------------------------------
real function enclosed_env_mass(r, C_rho)
 use physcon, only:pi

 real, intent(in) :: r, C_rho
 real :: exponent, integral

 exponent = 3.0 - rho_power

 if (abs(exponent) > 1.e-6) then
    integral = (r**exponent - r_inner**exponent) / exponent
 else
    integral = log(r / r_inner)
 endif

 enclosed_env_mass = 4.0 * pi * C_rho * integral

end function enclosed_env_mass

!-----------------------------------------------------------------------
!
!  Mass of the envelope between two radii r_a and r_b (r_a < r_b),
!
!-----------------------------------------------------------------------
real function region_mass(r_a_code, r_b_code)
 use units,   only:udist, umass
 use physcon, only:au

 real, intent(in) :: r_a_code, r_b_code  
 real :: r_a_cgs, r_b_cgs, C_rho

 r_a_cgs = r_a_code * udist
 r_b_cgs = r_b_code * udist
 C_rho   = calc_C_rho()

 region_mass = (enclosed_env_mass(r_b_cgs, C_rho) &
              - enclosed_env_mass(r_a_cgs, C_rho)) / umass

end function region_mass

!-----------------------------------------------------------------------
!
!  Initialize variables for stellar profile integration
!
!-----------------------------------------------------------------------
subroutine init_atmosphere(state)
 use physcon, only:pi, Rg, kboltz, mass_proton_cgs
 type(stellar_state), intent(out) :: state
 real :: C_rho

 ! Initialize at stellar surface
 state%r0    = Rstar_cgs
 state%r     = Rstar_cgs
 state%Rstar = Rstar_cgs
 state%P     = P0

 C_rho     = calc_C_rho()
 rho0      = C_rho / Rstar_cgs**rho_power
 state%rho = rho0
 state%u   = state%P / (state%rho * (Star_gamma - 1.))
 state%T   = Star_mu * mass_proton_cgs / kboltz * (Star_gamma - 1.) * state%u

 print *, ""
 print *, "Initial stellar surface conditions:"
 print *, " mu          :", Star_mu
 print *, " gamma       :", Star_gamma
 print *, " rho_power   :", rho_power
 print *, " rho (surf)  :", rho0
 print *, ""

 state%nsteps = 1
 state%error  = .false.

end subroutine init_atmosphere

!-----------------------------------------------------------------------
!
!  Integrate hydrostatic equilibrium over one radial step
!
!-----------------------------------------------------------------------
subroutine stellar_step(state, r_new)
 use physcon, only:Gg, pi, Rg, kboltz, mass_proton_cgs

 type(stellar_state), intent(inout) :: state
 real, intent(in) :: r_new
 real :: dr, r_mid, rho_mid, mr_mid, dP, C_rho

 dr    = r_new - state%r
 r_mid = 0.5 * (state%r + r_new)

 C_rho   = calc_C_rho()
 rho_mid = C_rho / r_mid**rho_power
 ! Total enclosed mass at midpoint = stellar core mass + envelope mass from r_inner to r_mid
 mr_mid  = Mstar_cgs + enclosed_env_mass(r_mid, C_rho)

 ! dP/dr = -rho * G * M(r) / r^2
 dP = -(Gg * mr_mid * rho_mid / r_mid**2) * dr

 state%r   = r_new
 state%P   = state%P + dP
 state%rho = C_rho / state%r**rho_power

 state%u = state%P / (state%rho * (Star_gamma - 1.))
 state%T = Star_mu * mass_proton_cgs / kboltz * (Star_gamma - 1.) * state%u

 state%nsteps = state%nsteps + 1

end subroutine stellar_step

!-----------------------------------------------------------------------
!
!  Integrate the hydrostatic equilibrium equation
!
!-----------------------------------------------------------------------
subroutine calc_stellar_profile(n)
 integer, intent(in) :: n
 type(stellar_state) :: state
 integer :: i
 real :: r_new, dr

 call init_atmosphere(state)

 if (allocated(stellar_1D)) deallocate(stellar_1D)
 allocate(stellar_1D(5, n))

 stellar_1D(1, n) = state%r
 stellar_1D(2, n) = state%rho
 stellar_1D(3, n) = state%P
 stellar_1D(4, n) = state%u
 stellar_1D(5, n) = state%T

 dr = (Rstar_cgs - r_inner) / real(n-1)

 do i = n-1, 1, -1
    r_new = Rstar_cgs - real(n-i)*dr
    call stellar_step(state, r_new)

    stellar_1D(1, i) = state%r
    stellar_1D(2, i) = state%rho
    stellar_1D(3, i) = state%P
    stellar_1D(4, i) = state%u
    stellar_1D(5, i) = state%T
 enddo

 call save_stellarprofile(n, 'stellar_profile1D.dat')

end subroutine calc_stellar_profile

!-----------------------------------------------------------------------
!
!  Interpolate stellar profile at given radius
!
!-----------------------------------------------------------------------
subroutine interp_stellar_profile(r, rho, P, u, T)
 use units,       only:udist, unit_density, unit_ergg, unit_pressure
 use table_utils, only:find_nearest_index, interp_1d
 use io,          only:fatal

 real, intent(in)  :: r
 real, intent(out) :: rho, P, u, T
 real :: r_cgs
 integer :: indx, n
 character(len=*), parameter :: label = 'interp_stellar_profile'

 if (.not. allocated(stellar_1D)) then
    call fatal(label, 'stellar_1D not allocated. Call setup_star first.')
 endif

 n     = size(stellar_1D, 2)
 r_cgs = r * udist

 if (r_cgs <= stellar_1D(1,1)) then
    print *, 'Warning: r_cgs < stellar_1D(1,1). Extrapolating...'
    rho = stellar_1D(2, 1) / unit_density
    P   = stellar_1D(3, 1) / unit_pressure
    u   = stellar_1D(4, 1) / unit_ergg
    T   = stellar_1D(5, 1)
    return
 elseif (r_cgs >= stellar_1D(1, n)) then
    rho = stellar_1D(2, n) / unit_density
    P   = stellar_1D(3, n) / unit_pressure
    u   = stellar_1D(4, n) / unit_ergg
    T   = stellar_1D(5, n)
    return
 endif

 call find_nearest_index(stellar_1D(1,:), r_cgs, indx)

 rho = interp_1d(r_cgs, stellar_1D(1,indx), stellar_1D(1,indx+1), &
                 stellar_1D(2,indx), stellar_1D(2,indx+1)) / unit_density
 P   = interp_1d(r_cgs, stellar_1D(1,indx), stellar_1D(1,indx+1), &
                 stellar_1D(3,indx), stellar_1D(3,indx+1)) / unit_pressure
 u   = interp_1d(r_cgs, stellar_1D(1,indx), stellar_1D(1,indx+1), &
                 stellar_1D(4,indx), stellar_1D(4,indx+1)) / unit_ergg
 T   = interp_1d(r_cgs, stellar_1D(1,indx), stellar_1D(1,indx+1), &
                 stellar_1D(5,indx), stellar_1D(5,indx+1))

end subroutine interp_stellar_profile

!-----------------------------------------------------------------------
!
!  Save stellar profile to file
!
!-----------------------------------------------------------------------
subroutine save_stellarprofile(n, filename)
 use physcon, only:au
 use io,      only:iverbose
 integer, intent(in) :: n
 character(*), intent(in) :: filename
 integer :: i, nwrite
 integer, parameter :: iunit = 1338

 if (iverbose >= 1) then
    write(*,'("Saving 1D stellar model to ",A)') trim(filename)
 endif

 open(unit=iunit, file=filename, status='replace')
 call filewrite_stellar_header(iunit, nwrite)

 do i = 1, n
    call filewrite_stellar_state(iunit, nwrite, i)
 enddo
 close(iunit)

end subroutine save_stellarprofile

subroutine filewrite_stellar_header(iunit, nwrite)
 integer, intent(in)  :: iunit
 integer, intent(out) :: nwrite
 character(len=20) :: fmt

 nwrite = 5
 write(fmt,*) nwrite
 write(iunit,'('// adjustl(fmt) //'(a15))') 'r','rho','P','u','T'

end subroutine filewrite_stellar_header

subroutine state_to_array(i, array)
 integer, intent(in)  :: i
 real,    intent(out) :: array(:)

 array(1) = stellar_1D(1, i)
 array(2) = stellar_1D(2, i)
 array(3) = stellar_1D(3, i)
 array(4) = stellar_1D(4, i)
 array(5) = stellar_1D(5, i)

end subroutine state_to_array

subroutine filewrite_stellar_state(iunit, nwrite, i)
 integer, intent(in) :: iunit, nwrite, i
 real :: array(nwrite)
 character(len=20) :: fmt

 call state_to_array(i, array)
 write(fmt,*) nwrite
 write(iunit,'('// adjustl(fmt) //'(1x,es14.6E3:))') array(1:nwrite)

end subroutine filewrite_stellar_state

end module wind_pulsating