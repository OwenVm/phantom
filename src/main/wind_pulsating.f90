module wind_pulsating

 implicit none
 public :: setup_star
 public :: stellar_state, save_stellarprofile, interp_stellar_profile, calc_stellar_profile
 public :: region_mass
 private

 real :: T_power
 real :: Mstar_cgs
 real :: r_outer
 real :: Rstar_cgs
 real :: r_inner
 real :: Star_gamma
 real :: Star_mu
 real :: rho_inner_cgs
 real :: T_inner_cgs

 real, public :: rho_power_fit = 2.0

 real, dimension(:,:), allocatable, public :: stellar_1D

 type stellar_state
    real    :: r, r0, Rstar, rho, P, u, T
    integer :: nsteps
    logical :: error
 end type stellar_state

contains

subroutine setup_star(Mstar_in, r_min, r_max, mu_in, gamma_in, rho_inner_in, rho_power_in, T_inner_in, T_power_in)

 real, intent(in) :: Mstar_in, r_min, r_max
 real, intent(in) :: mu_in, gamma_in
 real, intent(in) :: rho_inner_in, rho_power_in  ! rho_power_in kept for interface compatibility but not used
 real, intent(in) :: T_inner_in, T_power_in

 Mstar_cgs     = Mstar_in
 Rstar_cgs     = r_min
 r_outer       = r_max
 r_inner       = r_min
 Star_gamma    = gamma_in
 Star_mu       = mu_in
 rho_inner_cgs = rho_inner_in
 T_inner_cgs   = T_inner_in
 T_power       = T_power_in

 print *, "Setting up star with parameters:"
 print *, "  r_inner        :", r_inner
 print *, "  r_outer        :", r_outer
 print *, "  Mstar_cgs      :", Mstar_cgs
 print *, "  rho_inner      :", rho_inner_cgs
 print *, "  T_inner        :", T_inner_cgs
 print *, "  T_power        :", T_power
 print *, "  (rho_power_in ignored; will be fitted from hydrostatic solution)"

end subroutine setup_star

!-----------------------------------------------------------------------
! Temperature profile T(r) = T_inner * (r_inner/r)^T_power
!-----------------------------------------------------------------------
real function T_profile(r_cgs)
 real, intent(in) :: r_cgs
 T_profile = T_inner_cgs * (r_inner / r_cgs)**T_power
end function T_profile

!-----------------------------------------------------------------------
! RHS of hydrostatic ODE: dP/dr = -rho(P,r) * G * M / r^2
! With ideal gas rho = P * mu * m_H / (k_B * T(r)), this becomes
! dP/dr = -(mu * m_H / k_B) * G * M * P / (T(r) * r^2)
!-----------------------------------------------------------------------
real function dPdr_rhs(r_cgs, P_cgs)
 use physcon, only:kboltz, mass_proton_cgs, Gg
 real, intent(in) :: r_cgs, P_cgs
 real :: T, rho
 T      = T_profile(r_cgs)
 rho    = P_cgs * Star_mu * mass_proton_cgs / (kboltz * T)
 dPdr_rhs = -rho * Gg * Mstar_cgs / r_cgs**2
end function dPdr_rhs

!-----------------------------------------------------------------------
! Single RK4 step for the pressure ODE
!-----------------------------------------------------------------------
subroutine rk4_step(r, P, dr, P_new)
 real, intent(in)  :: r, P, dr
 real, intent(out) :: P_new
 real :: k1, k2, k3, k4
 k1 = dPdr_rhs(r,              P)
 k2 = dPdr_rhs(r + 0.5*dr,     P + 0.5*dr*k1)
 k3 = dPdr_rhs(r + 0.5*dr,     P + 0.5*dr*k2)
 k4 = dPdr_rhs(r + dr,         P +     dr*k3)
 P_new = P + (dr/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4)
end subroutine rk4_step

!-----------------------------------------------------------------------
! Fit a power law rho ~ r^(-alpha) to the hydrostatic profile in log-log
! space.  Sets the module-level rho_power_fit.
!-----------------------------------------------------------------------
subroutine fit_rho_powerlaw(n)
 integer, intent(in) :: n
 integer :: i, nfit
 real    :: sum_x, sum_y, sum_xx, sum_xy, x, y, alpha, rho_min_threshold
 real    :: rho_ref

 ! Only fit over the region where density is physically meaningful
 ! (not underflowed to zero or negative)
 rho_ref           = stellar_1D(2, 1)
 rho_min_threshold = rho_ref * 1.0e-20

 sum_x  = 0.0;  sum_y  = 0.0
 sum_xx = 0.0;  sum_xy = 0.0
 nfit   = 0

 do i = 1, n
    if (stellar_1D(2, i) <= rho_min_threshold) cycle
    x = log(stellar_1D(1, i) / r_inner)   ! log(r/r_inner)
    y = log(stellar_1D(2, i))              ! log(rho)
    sum_x  = sum_x  + x
    sum_y  = sum_y  + y
    sum_xx = sum_xx + x*x
    sum_xy = sum_xy + x*y
    nfit   = nfit + 1
 enddo

 ! Linear regression in log-log: y = alpha*x + const  =>  alpha = slope
 alpha = (real(nfit)*sum_xy - sum_x*sum_y) / (real(nfit)*sum_xx - sum_x**2)

 ! rho_power_fit is defined as the positive exponent in rho ~ r^(-rho_power_fit)
 rho_power_fit = -alpha

 print *, ''
 print *, '  Hydrostatic rho power-law fit:'
 print *, '    alpha (rho ~ r^-alpha) :', rho_power_fit
 print *, '    fit points             :', nfit, ' of ', n
 print *, ''

end subroutine fit_rho_powerlaw

!-----------------------------------------------------------------------
! Build the 1D stellar profile by integrating the hydrostatic ODE.
! Replaces the old power-law eval_stellar_state.
!-----------------------------------------------------------------------
subroutine calc_stellar_profile(n)
 use physcon, only:kboltz, mass_proton_cgs

 integer, intent(in) :: n
 real    :: r_cgs, dr, P_cgs, P_next, rho, T, u
 integer :: i

 dr = (r_outer - r_inner) / real(n - 1)

 if (allocated(stellar_1D)) deallocate(stellar_1D)
 allocate(stellar_1D(5, n))

 ! Inner boundary condition from ideal gas + prescribed rho_inner
 T     = T_profile(r_inner)
 P_cgs = rho_inner_cgs * kboltz * T / (Star_mu * mass_proton_cgs)

 r_cgs = r_inner
 do i = 1, n
    T   = T_profile(r_cgs)
    rho = P_cgs * Star_mu * mass_proton_cgs / (kboltz * T)
    u   = P_cgs / (rho * (Star_gamma - 1.0))

    stellar_1D(1, i) = r_cgs
    stellar_1D(2, i) = rho
    stellar_1D(3, i) = P_cgs
    stellar_1D(4, i) = u
    stellar_1D(5, i) = T

    if (i < n) then
       call rk4_step(r_cgs, P_cgs, dr, P_next)
       P_cgs = P_next
    endif
    r_cgs = r_cgs + dr
 enddo

 print *, ""
 print *, "Stellar profile (inner boundary):"
 print *, "  rho :", stellar_1D(2,1), " g/cm^3"
 print *, "  T   :", stellar_1D(5,1), " K"
 print *, "  P   :", stellar_1D(3,1), " dyn/cm^2"
 print *, "Stellar profile (outer boundary):"
 print *, "  rho :", stellar_1D(2,n), " g/cm^3"
 print *, "  T   :", stellar_1D(5,n), " K"
 print *, "  P   :", stellar_1D(3,n), " dyn/cm^2"
 print *, ""

 call fit_rho_powerlaw(n)

 call save_stellarprofile(n, 'stellar_profile1D.dat')

end subroutine calc_stellar_profile

subroutine interp_stellar_profile(r, rho, P, u, T)
 use units,       only:udist, unit_density, unit_ergg, unit_pressure
 use table_utils, only:find_nearest_index, interp_1d
 use io,          only:fatal

 real, intent(in)  :: r
 real, intent(out) :: rho, P, u, T
 real    :: r_cgs
 integer :: indx, n
 character(len=*), parameter :: label = 'interp_stellar_profile'

 if (.not. allocated(stellar_1D)) call fatal(label, 'stellar_1D not allocated. Call calc_stellar_profile first.')

 n     = size(stellar_1D, 2)
 r_cgs = r * udist

 if (r_cgs <= stellar_1D(1,1)) then
    rho = stellar_1D(2,1) / unit_density
    P   = stellar_1D(3,1) / unit_pressure
    u   = stellar_1D(4,1) / unit_ergg
    T   = stellar_1D(5,1)
    return
 elseif (r_cgs >= stellar_1D(1,n)) then
    rho = stellar_1D(2,n) / unit_density
    P   = stellar_1D(3,n) / unit_pressure
    u   = stellar_1D(4,n) / unit_ergg
    T   = stellar_1D(5,n)
    return
 endif

 call find_nearest_index(stellar_1D(1,:), r_cgs, indx)

 rho = interp_1d(r_cgs, stellar_1D(1,indx), stellar_1D(1,indx+1), stellar_1D(2,indx), stellar_1D(2,indx+1)) / unit_density
 P   = interp_1d(r_cgs, stellar_1D(1,indx), stellar_1D(1,indx+1), stellar_1D(3,indx), stellar_1D(3,indx+1)) / unit_pressure
 u   = interp_1d(r_cgs, stellar_1D(1,indx), stellar_1D(1,indx+1), stellar_1D(4,indx), stellar_1D(4,indx+1)) / unit_ergg
 T   = interp_1d(r_cgs, stellar_1D(1,indx), stellar_1D(1,indx+1), stellar_1D(5,indx), stellar_1D(5,indx+1))

end subroutine interp_stellar_profile

!-----------------------------------------------------------------------
! Enclosed envelope mass — still useful for region_mass calls elsewhere.
! Uses the fitted power law as an approximation.
!-----------------------------------------------------------------------
real function region_mass(r_a_code, r_b_code)
 use physcon, only:pi
 use units,   only:udist, umass

 real, intent(in) :: r_a_code, r_b_code
 real :: C_rho, r_a, r_b, exponent

 C_rho    = rho_inner_cgs * r_inner**rho_power_fit
 r_a      = r_a_code * udist
 r_b      = r_b_code * udist
 exponent = 3.0 - rho_power_fit

 if (abs(exponent) > 1.e-6) then
    region_mass = 4.0*pi * C_rho * (r_b**exponent - r_a**exponent) / exponent / umass
 else
    region_mass = 4.0*pi * C_rho * log(r_b / r_a) / umass
 endif

end function region_mass

subroutine save_stellarprofile(n, filename)
 use io, only:iverbose

 integer,      intent(in) :: n
 character(*), intent(in) :: filename
 integer :: i
 integer, parameter :: iunit = 1338

 if (iverbose >= 1) write(*,'("Saving 1D stellar model to ",A)') trim(filename)

 open(unit=iunit, file=filename, status='replace')
 write(iunit,'(5(a15))') 'r','rho','P','u','T'
 do i = 1, n
    write(iunit,'(5(1x,es14.6E3))') stellar_1D(1,i), stellar_1D(2,i), stellar_1D(3,i), stellar_1D(4,i), stellar_1D(5,i)
 enddo
 close(iunit)

end subroutine save_stellarprofile

end module wind_pulsating