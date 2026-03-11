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
! Pressure is evaluated analytically at every radius by integrating
! dP/dr = -rho(r)*g(r) from the outer boundary (P=0 at Rstar) inward
! to r, assuming a power-law density profile rho = C_rho / r^rho_power
! and M(r) = Mstar + Menv(r).  This avoids the numerical precision
! problems that arise when forward-stepping with a very steep profile.
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
 real :: rho_power = 2.0

 ! input parameters
 real :: Mstar_cgs, Rstar_cgs, r_inner, Star_gamma, Star_mu, rho_inner_cgs
 real, dimension(:,:), allocatable, public :: stellar_1D

 type stellar_state
    real :: r, r0, Rstar, rho, P, u, T
    integer :: nsteps
    logical :: error
 end type stellar_state

contains

subroutine setup_star(Mstar_in, Rstar_in, r_min, mu_in, gamma_in, rho_inner_in, rho_power_in)
 use physcon, only:au, solarm

 real, intent(in)           :: Mstar_in, Rstar_in, r_min, mu_in, gamma_in, rho_inner_in
 real, intent(in), optional :: rho_power_in

 Mstar_cgs     = Mstar_in
 Rstar_cgs     = Rstar_in
 r_inner       = r_min
 Star_gamma    = gamma_in
 Star_mu       = mu_in
 rho_inner_cgs = rho_inner_in

 if (present(rho_power_in)) then
    rho_power = rho_power_in
 else
    rho_power = 2.0
 endif

 print *, "Setting up star with parameters:"
 print *, "Rstar_cgs  :", Rstar_cgs
 print *, "Mstar_cgs  :", Mstar_cgs
 print *, "rho_inner  :", rho_inner_cgs
 print *, "r_inner    :", r_inner
 print *, "rho_power  :", rho_power

end subroutine setup_star

!-----------------------------------------------------------------------
!
!  Normalization constant for the power-law density profile
!
!  rho(r) = C_rho / r^rho_power
!  Anchored at r_inner: C_rho = rho_inner_cgs * r_inner^rho_power
!
!-----------------------------------------------------------------------
real function calc_C_rho()
 calc_C_rho = rho_inner_cgs * r_inner**rho_power
end function calc_C_rho

!-----------------------------------------------------------------------
!
!  Enclosed envelope mass from r_inner to r
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
!  Mass of the envelope between two radii r_a and r_b (r_a < r_b)
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
!  Analytic pressure at radius r, anchored at P(Rstar) = 0.
!
!  Integrates dP/dr = -rho(r)*[G*M(r)/r^2] analytically from Rstar to r:
!
!    P(r) = G*Mstar*C_rho/(p+1) * (r^{-(p+1)} - Rstar^{-(p+1)})
!         + 4*pi*G*C_rho^2/(3-p) * (I1(r,Rstar) - r^{3-p} * I2(r,Rstar))
!
!  where:
!    I1 = (Rstar^{2-2p} - r^{2-2p}) / (2-2p)     [p=1: log(Rstar/r)]
!    I2 = (r^{-(p+1)}   - Rstar^{-(p+1)}) / (p+1) [same as stellar term coeff]
!
!  Degenerate cases (p=1 for I1, p=3 for envelope self-gravity) are
!  handled explicitly.
!
!-----------------------------------------------------------------------
real function analytic_pressure(r, C_rho)
 use physcon, only:pi, Gg

 real, intent(in) :: r, C_rho
 real :: Pval, I1val, I2val, pp

 pp = rho_power

 ! stellar point-mass term
 Pval = Gg * Mstar_cgs * C_rho / (pp + 1.0) * &
        (r**(-(pp+1.0)) - Rstar_cgs**(-(pp+1.0)))

 ! envelope self-gravity term
 if (abs(2.0 - 2.0*pp) > 1.e-6) then
    I1val = (Rstar_cgs**(2.0-2.0*pp) - r**(2.0-2.0*pp)) / (2.0 - 2.0*pp)
 else
    I1val = log(Rstar_cgs / r)
 endif
 I2val = (r**(-(pp+1.0)) - Rstar_cgs**(-(pp+1.0))) / (pp + 1.0)

 if (abs(3.0 - pp) > 1.e-6) then
    Pval = Pval + 4.0*pi * Gg * C_rho**2 / (3.0 - pp) * &
           (I1val - r**(3.0-pp) * I2val)
 else
    ! pp=3: enclosed mass diverges logarithmically; skip self-gravity term
    print *, "Warning: rho_power=3, envelope self-gravity term skipped"
 endif

 analytic_pressure = max(Pval, 0.0)

end function analytic_pressure

!-----------------------------------------------------------------------
!
!  Evaluate the full stellar state (rho, P, u, T) at radius r_cgs (cgs).
!
!-----------------------------------------------------------------------
subroutine eval_stellar_state(r_cgs, C_rho, rho, P, u, T)
 use physcon, only:kboltz, mass_proton_cgs

 real, intent(in)  :: r_cgs, C_rho
 real, intent(out) :: rho, P, u, T

 rho = C_rho / r_cgs**rho_power
 P   = analytic_pressure(r_cgs, C_rho)
 u   = P / (rho * (Star_gamma - 1.))
 T   = Star_mu * mass_proton_cgs / kboltz * (Star_gamma - 1.) * u

end subroutine eval_stellar_state

!-----------------------------------------------------------------------
!
!  Initialize the stellar state at r_inner (for compatibility/printing).
!
!-----------------------------------------------------------------------
subroutine init_atmosphere(state)
 use physcon, only:kboltz, mass_proton_cgs
 type(stellar_state), intent(out) :: state
 real :: C_rho, rho, P, u, T

 C_rho = calc_C_rho()
 call eval_stellar_state(r_inner, C_rho, rho, P, u, T)

 state%r0    = r_inner
 state%r     = r_inner
 state%Rstar = Rstar_cgs
 state%rho   = rho
 state%P     = P
 state%u     = u
 state%T     = T
 state%nsteps = 1
 state%error  = .false.

 print *, ""
 print *, "Initial inner boundary conditions (analytic):"
 print *, " mu          :", Star_mu
 print *, " gamma       :", Star_gamma
 print *, " rho_power   :", rho_power
 print *, " rho (inner) :", rho
 print *, " P   (inner) :", P
 print *, " T   (inner) :", T
 print *, ""

end subroutine init_atmosphere

!-----------------------------------------------------------------------
!
!  Compute the full stellar profile by direct analytic evaluation at
!  each grid point.  No numerical stepping involved.
!
!-----------------------------------------------------------------------
subroutine calc_stellar_profile(n)
 integer, intent(in) :: n
 real :: C_rho, r_cgs, dr, rho, P, u, T
 integer :: i

 call init_atmosphere_print()

 C_rho = calc_C_rho()
 dr    = (Rstar_cgs - r_inner) / real(n-1)

 if (allocated(stellar_1D)) deallocate(stellar_1D)
 allocate(stellar_1D(5, n))

 do i = 1, n
    r_cgs = r_inner + real(i-1) * dr
    call eval_stellar_state(r_cgs, C_rho, rho, P, u, T)
    stellar_1D(1, i) = r_cgs
    stellar_1D(2, i) = rho
    stellar_1D(3, i) = P
    stellar_1D(4, i) = u
    stellar_1D(5, i) = T
 enddo

 call save_stellarprofile(n, 'stellar_profile1D.dat')

end subroutine calc_stellar_profile

! helper: print inner BC without duplicating init_atmosphere logic
subroutine init_atmosphere_print()
 type(stellar_state) :: state
 call init_atmosphere(state)
end subroutine init_atmosphere_print

!-----------------------------------------------------------------------
!
!  Interpolate stellar profile at given radius (code units)
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
 use io, only:iverbose
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