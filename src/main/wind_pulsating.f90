module wind_pulsating

 implicit none
 public :: setup_star
 public :: stellar_state, save_stellarprofile, interp_stellar_profile, calc_stellar_profile
 public :: region_mass, eval_boundary_state
 private

 real :: rho_power
 real :: T_power
 real :: Mstar_cgs
 real :: r_outer
 real :: Rstar_cgs
 real :: r_inner
 real :: Star_gamma
 real :: Star_mu
 real :: rho_inner_cgs
 real :: T_inner_cgs

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
 real, intent(in) :: rho_inner_in, rho_power_in
 real, intent(in) :: T_inner_in, T_power_in

 Mstar_cgs      = Mstar_in
 Rstar_cgs = r_min
 r_outer        = r_max
 r_inner        = r_min
 Star_gamma     = gamma_in
 Star_mu        = mu_in
 rho_inner_cgs  = rho_inner_in
 rho_power      = rho_power_in
 T_inner_cgs    = T_inner_in
 T_power        = T_power_in

 print *, "Setting up star with parameters:"
 print *, "  r_inner        :", r_inner
 print *, "  r_outer        :", r_outer
 print *, "  Mstar_cgs      :", Mstar_cgs
 print *, "  rho_inner      :", rho_inner_cgs
 print *, "  rho_power      :", rho_power
 print *, "  T_inner        :", T_inner_cgs
 print *, "  T_power        :", T_power

end subroutine setup_star

real function calc_C_rho()

 calc_C_rho = rho_inner_cgs * r_inner**rho_power

end function calc_C_rho

real function enclosed_env_mass(r, C_rho)
 use physcon, only:pi

 real, intent(in) :: r, C_rho
 real :: exponent

 exponent = 3.0 - rho_power

 if (abs(exponent) > 1.e-6) then
    enclosed_env_mass = 4.0*pi * C_rho * (r**exponent - r_inner**exponent) / exponent
 else
    enclosed_env_mass = 4.0*pi * C_rho * log(r / r_inner)
 endif

end function enclosed_env_mass

real function region_mass(r_a_code, r_b_code)
 use units, only:udist, umass

 real, intent(in) :: r_a_code, r_b_code
 real :: C_rho

 C_rho = calc_C_rho()
 region_mass = (enclosed_env_mass(r_b_code * udist, C_rho) - enclosed_env_mass(r_a_code * udist, C_rho)) / umass

end function region_mass

subroutine eval_stellar_state(r_cgs, rho, P, u, T)
 use physcon, only:kboltz, mass_proton_cgs

 real, intent(in)  :: r_cgs
 real, intent(out) :: rho, P, u, T
 real :: C_rho

 C_rho = calc_C_rho()

 rho = C_rho / r_cgs**rho_power
 T   = T_inner_cgs * (r_inner / r_cgs)**T_power
 P   = kboltz * rho * T / (Star_mu * mass_proton_cgs)
 u   = P / (rho * (Star_gamma - 1.0))

end subroutine eval_stellar_state

subroutine eval_boundary_state(r_cgs, rho, P, u, T)
 use physcon, only:kboltz, mass_proton_cgs, Gg
 use units,   only:unit_density, unit_ergg, unit_pressure

 real, intent(in)  :: r_cgs
 real, intent(out) :: rho, P, u, T
 real :: P_star, P_cgs, u_cgs

 T      = T_inner_cgs
 P_star = kboltz * rho_inner_cgs * T_inner_cgs / (Star_mu * mass_proton_cgs)
 P_cgs  = P_star + Gg * Mstar_cgs / (2.0 * Rstar_cgs**3) * rho_inner_cgs * (Rstar_cgs**2 - r_cgs**2)
 P_cgs  = max(P_cgs, P_star)

 ! u set from T_inner directly, not from P, so Phantom sees correct temperature
 u_cgs  = P_cgs / (rho_inner_cgs * (Star_gamma - 1.0))

 rho = rho_inner_cgs / unit_density
 P   = P_cgs        / unit_pressure
 u   = u_cgs        / unit_ergg

end subroutine eval_boundary_state

subroutine calc_stellar_profile(n)

 integer, intent(in) :: n
 real    :: r_cgs, dr, rho, P, u, T
 integer :: i

 dr = (r_outer - r_inner) / real(n-1)

 if (allocated(stellar_1D)) deallocate(stellar_1D)
 allocate(stellar_1D(5, n))

 do i = 1, n
    r_cgs = r_inner + real(i-1) * dr
    call eval_stellar_state(r_cgs, rho, P, u, T)
    stellar_1D(1, i) = r_cgs
    stellar_1D(2, i) = rho
    stellar_1D(3, i) = P
    stellar_1D(4, i) = u
    stellar_1D(5, i) = T
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