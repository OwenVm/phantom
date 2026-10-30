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
! With itgas_profile=1 the gas is given a prescribed temperature instead of the
! one needed to hold up the density profile. The difference between the
! hydrostatic and the gas pressure gradients is tabulated as the acceleration
! a_sup, which stands in for the support by pulsation-driven shocks.
!
! :References: None
!
! :Owner: Owen Vermeulen
!
! :Runtime parameters: None
!
! :Dependencies: dust_formation, eos, io, physcon, table_utils, units
!
 implicit none
 public :: setup_star
 public :: stellar_state,save_stellarprofile,interp_stellar_profile,calc_stellar_profile
 public :: region_mass,get_support_table,print_gas_profile

 private
 real :: rho_power = 2.0

 ! input parameters
 real :: Mstar_cgs, Rstar_cgs, Tstar_cgs, r_inner, Star_gamma, Star_mu, rho_inner_cgs
 ! Rstar_cgs is the outer radius of the profile, Reff_cgs the radius of the star
 real    :: Reff_cgs = 0., tgas_exp = 0.5
 integer :: itgas_profile = 0
 real, dimension(:,:), allocatable, public :: stellar_1D

 ! columns of stellar_1D, in cgs. P is P_HSE, u, T, mu and gamma are those of the gas
 integer, parameter :: ir = 1, irho = 2, iP = 3, iu = 4, iT = 5, iPgas = 6, imu = 7, igam = 8, iasup = 9
 integer, parameter :: ncols = 9

 ! default relative change in any quantity before a new profile point is stored
 real, parameter :: eps_default = 0.01

 type stellar_state
    real :: r, rho, P, u, T, P_gas, mu, gamma
 end type stellar_state

contains

subroutine setup_star(Mstar_in, Tstar_in, Rstar_in, r_min, mu_in, gamma_in, rho_inner_in, rho_power_in, &
                      Reff_in, itgas_profile_in, tgas_exp_in)
 real,    intent(in)           :: Mstar_in, Tstar_in, Rstar_in, r_min, mu_in, gamma_in, rho_inner_in
 real,    intent(in), optional :: rho_power_in, Reff_in, tgas_exp_in
 integer, intent(in), optional :: itgas_profile_in

 Mstar_cgs     = Mstar_in
 Rstar_cgs     = Rstar_in
 Tstar_cgs     = Tstar_in
 r_inner       = r_min
 Star_gamma    = gamma_in
 Star_mu       = mu_in
 rho_inner_cgs = rho_inner_in

 rho_power = 2.0
 if (present(rho_power_in)) rho_power = rho_power_in

 itgas_profile = 0
 if (present(itgas_profile_in)) itgas_profile = itgas_profile_in
 if (present(tgas_exp_in))      tgas_exp      = tgas_exp_in
 if (present(Reff_in))          Reff_cgs      = Reff_in

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
!  Set the gas properties. With itgas_profile=0, u and T follow from P and
!  rho using the ideal gas law. With itgas_profile=1 the temperature is
!  prescribed, and u uses the conversion in eos.f90 for ieos=5:
!
!    T = mu (gamma-1) u / Rg
!
!-----------------------------------------------------------------------
subroutine set_u_and_T(state)
 use physcon,        only:kboltz, mass_proton_cgs, Rg
 use dust_formation, only:init_muGamma
 type(stellar_state), intent(inout) :: state

 if (itgas_profile == 1) then
    state%T = Tstar_cgs * (Reff_cgs / state%r)**tgas_exp
    call init_muGamma(state%rho, state%T, state%mu, state%gamma)
    state%u     = Rg * state%T / (state%mu * (state%gamma - 1.))
    state%P_gas = (state%gamma - 1.) * state%rho * state%u
 else
    state%mu    = Star_mu
    state%gamma = Star_gamma
    state%u     = state%P / (state%rho * (Star_gamma - 1.))
    state%T     = Star_mu * mass_proton_cgs / kboltz * (Star_gamma - 1.) * state%u
    state%P_gas = state%P
 endif

end subroutine set_u_and_T

!-----------------------------------------------------------------------
!
!  Initialize variables for stellar profile integration at r_outer (Rstar).
!  Anchors pressure at the outer boundary using the ideal gas law at Tstar.
!  Integration then proceeds inward.
!
!-----------------------------------------------------------------------
subroutine init_atmosphere(state)
 use physcon, only:kboltz, mass_proton_cgs
 type(stellar_state), intent(out) :: state

 ! Anchor at outer boundary with ideal gas law
 state%r   = Rstar_cgs
 state%rho = calc_C_rho() / Rstar_cgs**rho_power
 state%P   = state%rho * kboltz * Tstar_cgs / (Star_mu * mass_proton_cgs)
 call set_u_and_T(state)

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
!  The run has no self-gravity, so with itgas_profile=1 only the sink mass
!  counts, as in the simulation.
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
 M_enc   = Mstar_cgs
 if (itgas_profile == 0) M_enc = M_enc + enclosed_env_mass(r_mid, C_rho)

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

 array = [state%r, state%rho, state%P, state%u, state%T, state%P_gas, state%mu, state%gamma, 0.]

end function state_to_array

!-----------------------------------------------------------------------
!
!  Integrate the hydrostatic equilibrium equation inward from Rstar to r_inner,
!  then reverse the array so it runs from r_inner to Rstar.
!
!  n is the number of integration steps. A point is only stored once any
!  quantity (r, rho, P, u, T, ...) has changed by more than a relative amount
!  eps (default eps_default) since the last stored point, so the stored
!  profile has at most n points. The innermost point is always stored.
!
!-----------------------------------------------------------------------
subroutine calc_stellar_profile(n, eps_in)
 integer, intent(in)        :: n
 real,    intent(in), optional :: eps_in
 type(stellar_state) :: state
 real, allocatable   :: tmp(:,:)
 real    :: dr, eps, new(ncols)
 integer :: i, nstored

 eps = eps_default
 if (present(eps_in)) eps = eps_in

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

 call calc_support_accel(nstored)

 print *, ""
 print *, "Inner boundary conditions (after inward integration):"
 print *, " r    (inner) :", stellar_1D(ir, 1)
 print *, " rho  (inner) :", stellar_1D(irho, 1)
 print *, " P    (inner) :", stellar_1D(iP, 1)
 print *, " T    (inner) :", stellar_1D(iT, 1)
 print *, " points stored:", nstored, " of ", n, " integration steps (eps =", eps, ")"
 print *, ""

 call save_stellarprofile(nstored, 'stellar_profile1D.dat')

end subroutine calc_stellar_profile

!-----------------------------------------------------------------------
!
!  First derivative of f on the (non-uniform) grid r, second order
!  accurate: centred in the interior and one-sided at both ends. The grid
!  is non-uniform since the profile is thinned, for which a plain centred
!  difference would only be first order.
!
!-----------------------------------------------------------------------
subroutine radial_derivative(r, f, dfdr)
 real, intent(in)  :: r(:), f(:)
 real, intent(out) :: dfdr(:)
 real    :: h1, h2
 integer :: i, n

 n = size(r)
 do i = 2, n-1
    h1 = r(i) - r(i-1)
    h2 = r(i+1) - r(i)
    dfdr(i) = (h1**2*f(i+1) - h2**2*f(i-1) + (h2**2 - h1**2)*f(i)) / (h1*h2*(h1 + h2))
 enddo
 h1 = r(2) - r(1)
 h2 = r(3) - r(2)
 dfdr(1) = -(2.*h1 + h2)/(h1*(h1 + h2))*f(1) + (h1 + h2)/(h1*h2)*f(2) - h1/(h2*(h1 + h2))*f(3)
 h1 = r(n-1) - r(n-2)
 h2 = r(n) - r(n-1)
 dfdr(n) = h2/(h1*(h1 + h2))*f(n-2) - (h1 + h2)/(h1*h2)*f(n-1) + (h1 + 2.*h2)/(h2*(h1 + h2))*f(n)

end subroutine radial_derivative

!-----------------------------------------------------------------------
!
!  Acceleration that holds up the profile against gravity with the gas
!  pressure, positive outward:
!
!    a_sup = (dP_gas/dr - dP_HSE/dr)/rho
!
!  from -(1/rho) dP_gas/dr - g + a_sup = 0 and dP_HSE/dr = -rho g.
!
!-----------------------------------------------------------------------
subroutine calc_support_accel(n)
 use io,      only:fatal, warning
 use physcon, only:Gg
 integer, intent(in) :: n
 real    :: dPgas(n), dPhse(n), g(n), residual
 character(len=*), parameter :: label = 'calc_support_accel'

 if (n < 3) call fatal(label, 'need at least 3 points in the stellar profile', var='n', ival=n)

 call radial_derivative(stellar_1D(ir,1:n), stellar_1D(iP,1:n), dPhse)
 call radial_derivative(stellar_1D(ir,1:n), stellar_1D(iPgas,1:n), dPgas)
 stellar_1D(iasup,1:n) = (dPgas - dPhse) / stellar_1D(irho,1:n)

 if (itgas_profile == 1) then
    g = Gg * Mstar_cgs / stellar_1D(ir,1:n)**2
    residual = maxval(abs(-dPgas/stellar_1D(irho,1:n) - g + stellar_1D(iasup,1:n)) / g)
    print *, " max force residual / g :", residual
    print *, " a_sup / g, min and max :", minval(stellar_1D(iasup,1:n)/g), maxval(stellar_1D(iasup,1:n)/g)
    if (residual > 1.e-3) call warning(label, 'large force residual, increase n_profile_points', val=residual)
 endif
 if (any(stellar_1D(iasup,1:n) < 0.)) &
    call warning(label, 'a_sup < 0: gas pressure gradient exceeds gravity', var='min(a_sup)', &
                 val=minval(stellar_1D(iasup,1:n)))

end subroutine calc_support_accel

!-----------------------------------------------------------------------
!
!  The support acceleration table in code units, for the force routine
!
!-----------------------------------------------------------------------
subroutine get_support_table(r, a_sup)
 use units, only:udist, utime
 real, allocatable, intent(out) :: r(:), a_sup(:)

 r     = stellar_1D(ir,:) / udist
 a_sup = stellar_1D(iasup,:) * utime**2 / udist

end subroutine get_support_table

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

 rho = vals(irho) / unit_density
 P   = vals(iPgas) / unit_pressure
 u   = vals(iu) / unit_ergg
 T   = vals(iT)
 if (present(mu))    mu    = vals(imu)
 if (present(gamma)) gamma = vals(igam)

end subroutine interp_stellar_profile

!-----------------------------------------------------------------------
!
!  Check that u_gas gives back T_gas through what the code does for
!  ieos=5: equationofstate in cons2prim, then calc_muGamma in the substeps.
!  Prints the worst relative error over the profile, then the gas
!  properties at the radii r (code units) for comparison.
!
!-----------------------------------------------------------------------
subroutine print_gas_profile(r)
 use dust_formation, only:calc_muGamma
 use eos,            only:ieos, equationofstate
 use io,             only:warning
 use units,          only:unit_density, unit_ergg
 real, intent(in) :: r(:)
 real    :: rho, P, u, T, mu, gamma, gamma_i, mu_i, ponrho, cs, T_eos, T_sub, pH, pH_tot, err
 integer :: i
 character(len=*), parameter :: label = 'print_gas_profile'

 if (itgas_profile == 1) then
    err = 0.
    do i = 1, size(stellar_1D, 2)
       rho     = stellar_1D(irho, i)
       mu_i    = stellar_1D(imu, i)
       gamma_i = stellar_1D(igam, i)
       call equationofstate(ieos, ponrho, cs, rho/unit_density, 0., 0., 0., T_eos, &
                            eni=stellar_1D(iu, i)/unit_ergg, gamma_local=gamma_i, mu_local=mu_i)
       T_sub = T_eos
       call calc_muGamma(rho, T_sub, mu_i, gamma_i, pH, pH_tot)
       err = max(err, abs(T_eos/stellar_1D(iT, i) - 1.), abs(T_sub/stellar_1D(iT, i) - 1.))
    enddo
    print *, ""
    print *, " max |T_eos/T_gas - 1| over the profile :", err
    if (err > 1.e-6) call warning(label, 'u_gas does not return T_gas', val=err)
 endif

 print "(a)", "        r       rho(cgs)     T_gas        mu       gamma       T_eos    "
 do i = 1, size(r)
    call interp_stellar_profile(r(i), rho, P, u, T, mu, gamma)
    mu_i    = mu
    gamma_i = gamma
    call equationofstate(ieos, ponrho, cs, rho, 0., 0., 0., T_eos, eni=u, gamma_local=gamma_i, mu_local=mu_i)
    print "(6(1x,es11.4))", r(i), rho*unit_density, T, mu, gamma, T_eos
 enddo
 print *, ""

end subroutine print_gas_profile

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
 write(iunit,'(9(a15))') 'r','rho','P_HSE','u_gas','T_gas','P_gas','mu','gamma','a_sup'
 do i = 1, n
    write(iunit,'(9(1x,es14.6E3:))') stellar_1D(:, i)
 enddo
 close(iunit)

end subroutine save_stellarprofile

end module wind_pulsating
