!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
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
! :Dependencies: dim, dust_formation, io, physcon, table_utils, units
!
 implicit none
 public :: setup_star
 public :: stellar_state,save_stellarprofile,interp_stellar_profile,calc_stellar_profile
 public :: region_mass

 private
 real :: rho_power = 2.0

 ! input parameters
 real :: Mstar_cgs, Rstar_cgs, Tstar_cgs, r_inner, Star_gamma, Star_mu, rho_inner_cgs
 real, dimension(:,:), allocatable, public :: stellar_1D

 ! columns of stellar_1D: r, rho, P, u, T, mu, gamma
 integer, parameter :: ncols = 7

 ! number of profile points where T, mu and gamma did not converge
 integer :: n_not_converged = 0

 ! default relative change in any quantity before a new profile point is stored
 real, parameter :: eps_default = 0.01

 type stellar_state
    real :: r, rho, P, u, T, mu, gamma
 end type stellar_state

contains

subroutine setup_star(Mstar_in, Tstar_in, Rstar_in, r_min, mu_in, gamma_in, rho_inner_in, rho_power_in)
 real, intent(in)           :: Mstar_in, Tstar_in, Rstar_in, r_min, mu_in, gamma_in, rho_inner_in
 real, intent(in), optional :: rho_power_in

 Mstar_cgs     = Mstar_in
 Rstar_cgs     = Rstar_in
 Tstar_cgs     = Tstar_in
 r_inner       = r_min
 Star_gamma    = gamma_in
 Star_mu       = mu_in
 rho_inner_cgs = rho_inner_in

 rho_power = 2.0
 if (present(rho_power_in)) rho_power = rho_power_in

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
 real :: exponent

 exponent = 3.0 - rho_power
 enclosed_env_mass = 4.0*pi * C_rho * (r**exponent - r_inner**exponent) / exponent

end function enclosed_env_mass

!-----------------------------------------------------------------------
!
!  Mass of the envelope between two radii r_a and r_b (r_a < r_b)
!
!-----------------------------------------------------------------------
real function region_mass(r_a_code, r_b_code)
 use units, only:udist, umass
 real, intent(in) :: r_a_code, r_b_code
 real :: C_rho

 C_rho       = calc_C_rho()
 region_mass = (enclosed_env_mass(r_b_code*udist, C_rho) - enclosed_env_mass(r_a_code*udist, C_rho)) / umass

end function region_mass

!-----------------------------------------------------------------------
!
!  Set T, mu, gamma and u from P and rho using the ideal gas law.
!  When mu and gamma depend on the H2 formation (ieos=5), T is iterated
!  with mu(rho,T) and gamma(rho,T) from calc_muGamma until consistent,
!  so that u = P/((gamma-1) rho) matches the equation of state of the run.
!
!-----------------------------------------------------------------------
subroutine set_u_and_T(state)
 use physcon,        only:kboltz, mass_proton_cgs
 use dim,            only:update_muGamma
 use dust_formation, only:calc_muGamma
 use units,          only:unit_density
 type(stellar_state), intent(inout) :: state
 integer, parameter :: itermax = 100
 real,    parameter :: tol = 1.e-4   ! calc_muGamma itself converges to 1e-3
 real    :: T_old, T_in, pH, pH_tot
 integer :: iter
 logical :: converged

 state%mu    = Star_mu
 state%gamma = Star_gamma
 state%T     = state%mu * mass_proton_cgs * state%P / (kboltz * state%rho)

 if (update_muGamma) then
    converged = .false.
    do iter = 1, itermax
       T_old = state%T
       T_in  = state%T   ! calc_muGamma may modify its temperature argument
       call calc_muGamma(state%rho/unit_density, T_in, state%mu, state%gamma, pH, pH_tot)
       state%T = state%mu * mass_proton_cgs * state%P / (kboltz * state%rho)
       if (abs(state%T - T_old) < tol*T_old) then
          converged = .true.
          exit
       endif
    enddo
    if (.not. converged) n_not_converged = n_not_converged + 1
 endif

 state%u = state%P / (state%rho * (state%gamma - 1.))

end subroutine set_u_and_T

!-----------------------------------------------------------------------
!
!  Initialize variables for stellar profile integration at r_outer (Rstar).
!  Anchors pressure at the outer boundary using the ideal gas law at Tstar
!  (with the mu of the equation of state, so that T = Tstar there).
!  Integration then proceeds inward.
!
!-----------------------------------------------------------------------
subroutine init_atmosphere(state)
 use physcon, only:kboltz, mass_proton_cgs
 type(stellar_state), intent(out) :: state
 integer :: iter

 ! Anchor at outer boundary with ideal gas law
 state%r   = Rstar_cgs
 state%rho = calc_C_rho() / Rstar_cgs**rho_power
 state%P   = state%rho * kboltz * Tstar_cgs / (Star_mu * mass_proton_cgs)
 call set_u_and_T(state)
 ! if mu depends on T (ieos=5), rescale P until the outer temperature is Tstar
 do iter = 1, 50
    if (abs(state%T - Tstar_cgs) < 1.e-6*Tstar_cgs) exit
    state%P = state%P * Tstar_cgs / state%T
    call set_u_and_T(state)
 enddo

 print *, ""
 print *, "Outer boundary conditions (integration start):"
 print *, " r    (outer) :", state%r
 print *, " rho  (outer) :", state%rho
 print *, " P    (outer) :", state%P
 print *, " T    (outer) :", state%T
 print *, ""

end subroutine init_atmosphere

!-----------------------------------------------------------------------
!
!  Integrate hydrostatic equilibrium over one radial step.
!  Works for both inward (dr < 0) and outward (dr > 0) steps.
!
!-----------------------------------------------------------------------
subroutine stellar_step(state, r_new)
 use physcon, only:Gg
 type(stellar_state), intent(inout) :: state
 real, intent(in) :: r_new
 real :: r_mid, rho_mid, C_rho, M_enc

 r_mid   = 0.5 * (state%r + r_new)
 C_rho   = calc_C_rho()
 rho_mid = C_rho / r_mid**rho_power
 M_enc   = Mstar_cgs + enclosed_env_mass(r_mid, C_rho)

 state%P   = state%P - Gg * (M_enc * rho_mid / r_mid**2) * (r_new - state%r)
 state%r   = r_new
 state%rho = C_rho / state%r**rho_power
 call set_u_and_T(state)

end subroutine stellar_step

!-----------------------------------------------------------------------
!
!  Pack a stellar state into an array with the column order of stellar_1D
!
!-----------------------------------------------------------------------
function state_to_array(state) result(array)
 type(stellar_state), intent(in) :: state
 real :: array(ncols)

 array = [state%r, state%rho, state%P, state%u, state%T, state%mu, state%gamma]

end function state_to_array

!-----------------------------------------------------------------------
!
!  Integrate the hydrostatic equilibrium equation inward from Rstar to r_inner,
!  then reverse the array so it runs from r_inner to Rstar.
!
!  n is the number of integration steps. A point is only stored once any
!  quantity (r, rho, P, u, T, mu, gamma) has changed by more than a relative amount
!  eps (default eps_default) since the last stored point, so the stored
!  profile has at most n points. The innermost point is always stored.
!
!-----------------------------------------------------------------------
subroutine calc_stellar_profile(n, eps_in)
 use io, only:warning
 integer, intent(in)        :: n
 real,    intent(in), optional :: eps_in
 type(stellar_state) :: state
 real, allocatable   :: tmp(:,:)
 real    :: dr, eps, new(ncols)
 integer :: i, nstored

 eps = eps_default
 if (present(eps_in)) eps = eps_in

 n_not_converged = 0
 call init_atmosphere(state)
 allocate(tmp(ncols, n))

 ! dr is negative — stepping inward
 dr = (r_inner - Rstar_cgs) / real(n-1)

 nstored   = 1
 tmp(:, 1) = state_to_array(state)
 do i = 2, n
    call stellar_step(state, Rstar_cgs + real(i-1) * dr)
    new = state_to_array(state)
    if (i == n .or. any(abs(new - tmp(:, nstored)) > eps * abs(tmp(:, nstored)))) then
       nstored = nstored + 1
       tmp(:, nstored) = new
    endif
 enddo

 ! Reverse so stellar_1D runs from r_inner (index 1) to Rstar (index nstored)
 if (allocated(stellar_1D)) deallocate(stellar_1D)
 allocate(stellar_1D(ncols, nstored))
 do i = 1, nstored
    stellar_1D(:, i) = tmp(:, nstored+1-i)
 enddo
 deallocate(tmp)

 print *, ""
 print *, "Inner boundary conditions (after inward integration):"
 print *, " r    (inner) :", stellar_1D(1, 1)
 print *, " rho  (inner) :", stellar_1D(2, 1)
 print *, " P    (inner) :", stellar_1D(3, 1)
 print *, " T    (inner) :", stellar_1D(5, 1)
 print *, " points stored:", nstored, " of ", n, " integration steps (eps =", eps, ")"
 if (n_not_converged > 0) then
    call warning('calc_stellar_profile','no consistent T, mu, gamma at some profile points '// &
                 '(mu(T) is discontinuous there)',var='n_points',ival=n_not_converged)
 endif
 print *, ""

 call save_stellarprofile(nstored, 'stellar_profile1D.dat')

end subroutine calc_stellar_profile

!-----------------------------------------------------------------------
!
!  Interpolate stellar profile at given radius
!
!-----------------------------------------------------------------------
subroutine interp_stellar_profile(r, rho, P, u, T, mu, gamma)
 use units,       only:udist, unit_density, unit_ergg, unit_pressure
 use table_utils, only:find_nearest_index, interp_1d
 use io,          only:fatal
 real, intent(in)  :: r
 real, intent(out) :: rho, P, u, T
 real, intent(out), optional :: mu, gamma
 real    :: r_cgs, vals(ncols)
 integer :: indx, n, j
 character(len=*), parameter :: label = 'interp_stellar_profile'

 if (.not. allocated(stellar_1D)) call fatal(label, 'stellar_1D not allocated. Call setup_star first.')

 n     = size(stellar_1D, 2)
 r_cgs = r * udist

 if (r_cgs <= stellar_1D(1, 1)) then
    vals = stellar_1D(:, 1)
 elseif (r_cgs >= stellar_1D(1, n)) then
    vals = stellar_1D(:, n)
 else
    call find_nearest_index(stellar_1D(1,:), r_cgs, indx)
    do j = 2, ncols
       vals(j) = interp_1d(r_cgs, stellar_1D(1,indx), stellar_1D(1,indx+1), stellar_1D(j,indx), stellar_1D(j,indx+1))
    enddo
 endif

 rho = vals(2) / unit_density
 P   = vals(3) / unit_pressure
 u   = vals(4) / unit_ergg
 T   = vals(5)
 if (present(mu))    mu    = vals(6)
 if (present(gamma)) gamma = vals(7)

end subroutine interp_stellar_profile

!-----------------------------------------------------------------------
!
!  Save stellar profile to file
!
!-----------------------------------------------------------------------
subroutine save_stellarprofile(n, filename)
 use io, only:iverbose
 integer,      intent(in) :: n
 character(*), intent(in) :: filename
 integer :: i, iunit

 if (iverbose >= 1) write(*,'("Saving 1D stellar model to ",A)') trim(filename)

 open(newunit=iunit,file=filename,status='replace')
 write(iunit,'(7(a15))') 'r','rho','P','u','T','mu','gamma'
 do i = 1, n
    write(iunit,'(7(1x,es14.6E3:))') stellar_1D(:, i)
 enddo
 close(iunit)

end subroutine save_stellarprofile

end module wind_pulsating
