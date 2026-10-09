!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module cooling_functions
!
! A library of cooling functions that can be handled by cooling_solver
!  Contributed by Lionel Siess and Ward Homan
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: physcon
!
 implicit none

 real, public :: bowen_Cprime     = 3.000d-5
 real, public :: lambda_shock_cgs = 1.d0
 real, public :: T1_factor = 20., T0_value = 0.
 real, public :: kappa_dust_min = 1e-3  ! dust opacity value below which dust cooling is not calculated
 real, public :: CO_abun = 3.e-4        ! n_CO/n_H2 used for CO rotational line cooling

 ! parameters of the Decin et al. (2006) heating/cooling terms
 real, public :: H2O_abun    = 1.25e-4    ! n_H2O/n_H2 used for H2O rotational line cooling
 real, public :: dust_to_gas = 0.01     ! dust-to-gas mass ratio psi
 real, public :: grain_amin  = 5.e-7    ! minimum grain size (cm)
 real, public :: grain_amax  = 2.5e-5   ! maximum grain size (cm)
 real, public :: grain_rhos  = 3.3      ! specific density of the dust grains (g/cm^3)
 real, public :: G0_UV       = 1.       ! interstellar far-UV field in Habing units
 real, public :: r_half_CO   = 0.       ! CO photodissociation radius r_1/2 (cm), 0 = no dissociation
 real, public :: alpha_CO    = 2.5      ! exponent of the CO photodissociation profile (Mamon et al. 1988)

 real, parameter, public :: He_abun = 1.04e-1 ! n_He/n_H, as in dust_formation

 public :: cool_dust_discrete_contact, cool_coulomb, &
           cool_HI, cool_H_ionisation, cool_He_ionisation, &
           cool_H2_rovib, cool_H2_dissociation, cool_CO_rovib, &
           cool_H2O_rovib, cool_OH_rot, heat_dust_friction, &
           heat_dust_photovoltaic_soft, heat_CosmicRays, &
           heat_H2_recombination, &
           cool_dust_full_contact, cool_dust_radiation, &
           cooling_neutral_hydrogen, cool_metal_ions, &
           cool_thermal_bremsstrahlung, &
           heat_Compton, heat_recombination, &
           heat_dust_photovoltaic_hard, &
           piecewise_law, cooling_high_temp, &
           cooling_Bowen_relaxation, &
           cooling_dust_collision, &
           cooling_radiative_relaxation, &
           cooling_H2, cooling_CO_rot, &
           testing_cooling_functions, &
           decin_densities, heating_drift_decin, cooling_H2O_rot_decin, &
           cooling_CO_rot_decin, cooling_H2_vib_decin, &
           heating_dust_gas_decin, heating_cosmic_rays_decin, &
           heating_photoelectric_decin

 private
 real, parameter  :: xH = 0.7, xHe = 0.28 !assumed H and He mass fractions
 integer, parameter :: ngrain    = 32     ! grain-size bins for the Decin et al. (2006) integrals
 real,    parameter :: v_sputter = 2.e6   ! drift velocity above which grains are sputtered (cm/s)

contains
!-----------------------------------------------------------------------
!+
!  Piecewise cooling law for simple shock problem (Creasey et al. 2011)
!+
!-----------------------------------------------------------------------
subroutine piecewise_law(T, T0, rho_cgs, ndens, Q, dlnQ)

 real, intent(in)  :: T, T0, rho_cgs, ndens
 real, intent(out) :: Q, dlnQ
 real :: T1,Tmid !,dlnT,fac

 T1 = T1_factor*T0
 Tmid = 0.5*(T0+T1)
 if (T < T0) then
    Q    = 0.
    dlnQ = 0.
 elseif (T >= T0 .and. T <= Tmid) then
    !dlnT = (T-T0)/(T0/100.)
    Q = -lambda_shock_cgs*ndens**2/rho_cgs*(T-T0)/T0
    !fac = 2./(1.d0 + exp(dlnT))
    dlnQ = 1./(T-T0+epsilon(0.))
 elseif (T >= Tmid .and. T <= T1) then
    Q = -lambda_shock_cgs*ndens**2/rho_cgs*(T1-T)/T0
    dlnQ = -1./(T1-T+epsilon(0.))
 else
    Q    = 0.
    dlnQ = 0.
 endif
 !derivatives are discontinuous!

end subroutine piecewise_law

!-----------------------------------------------------------------------
!+
!  Bowen 1988 cooling prescription
!+
!-----------------------------------------------------------------------
subroutine cooling_Bowen_relaxation(T, Tdust, rho_cgs, mu, gamma, Q_cgs, dlnQ_dlnT)

 use physcon, only:Rg

 real, intent(in)  :: T, Tdust, rho_cgs, mu, gamma
 real, intent(out) :: Q_cgs, dlnQ_dlnT

 Q_cgs     = Rg/((gamma-1.)*mu)*rho_cgs*(Tdust-T)/bowen_Cprime
 dlnQ_dlnT = -T/(Tdust-T+1.d-10)

end subroutine cooling_Bowen_relaxation

!-----------------------------------------------------------------------
!+
!  collisionnal cooling
!+
!-----------------------------------------------------------------------
subroutine cooling_dust_collision(T, Tdust, rho, K2, mu, Q_cgs, dlnQ_dlnT)

 use physcon, only:kboltz, mass_proton_cgs, pi

 real, intent(in)  :: T, Tdust, rho, K2, mu
 real, intent(out) :: Q_cgs, dlnQ_dlnT

 real, parameter   :: f = 0.15, a0 = 1.28e-8
 real              :: A

 A = 2. * f * kboltz * a0**2/(mass_proton_cgs**2*mu) &
         * (1.05/1.54) * sqrt(2.*pi*kboltz/mass_proton_cgs) * 2.*K2 * rho
 Q_cgs = A * sqrt(T) * (Tdust-T)
 if (Q_cgs  >  1.d6) then
    print *, f, kboltz, a0, mass_proton_cgs, mu
    print *, mu, K2, rho, T, Tdust, A, Q_cgs
    stop 'cooling'
 else
    dlnQ_dlnT = 0.5+T/(Tdust-T+1.d-10)
 endif

end subroutine cooling_dust_collision

!-----------------------------------------------------------------------
!+
!  Woitke (2006 A&A) cooling term
!+
!-----------------------------------------------------------------------
subroutine cooling_radiative_relaxation(T, Tdust, kappa, Q_cgs, dlnQ_dlnT)

 use physcon, only:steboltz

 real, intent(in)  :: T, Tdust, kappa
 real, intent(out) :: Q_cgs, dlnQ_dlnT

 Q_cgs     = 4.*steboltz*(Tdust**4-T**4)*kappa
 dlnQ_dlnT = -4.*T**4/(Tdust**4-T**4+1.d-10)

end subroutine cooling_radiative_relaxation

!-----------------------------------------------------------------------
!+
!  Cooling due to electron excitation of neutral H (Spitzer 1978)
!+
!-----------------------------------------------------------------------
subroutine cooling_neutral_hydrogen(T, rho_cgs, Q_cgs, dlnQ_dlnT)

 use physcon, only:mass_proton_cgs

 real, intent(in)  :: T, rho_cgs
 real, intent(out) :: Q_cgs,dlnQ_dlnT

 real, parameter   :: f = 1.0d0
 real              :: ne,nH

 if (T > 3000. .and. T < 1.2e+04) then
    nH = rho_cgs/(1.4*mass_proton_cgs)
    ne = min(calc_eps_e(T),1.)*nH
    !the term 1/(1+sqrt(T)) comes from Cen (1992, ApjS, 78, 341)
    Q_cgs  = -f*7.3d-19*ne*nH*exp(-118400./T)/rho_cgs/(1.+sqrt(T/1.d5))
    dlnQ_dlnT = -118400./T+log(nH*calc_eps_e(1.001*T)/ne)/log(1.001) &
         - 0.5*sqrt(T/1.d5)/(1.+sqrt(T/1.d5))
 else
    Q_cgs = 0.
    dlnQ_dlnT = 0.
 endif

end subroutine cooling_neutral_hydrogen

!-----------------------------------------------------------------------
!+
!  high temperatures cooling & heating due to recombination, metal ions,
!  thermal bremsstrahlung & Compton contribution (Mathews & Doane 1990)
!+
!-----------------------------------------------------------------------
subroutine cooling_high_temp(T, rho_cgs, Q_cgs, dlnQ_dlnT)
 use physcon,        only: mass_proton_cgs
 real, intent(in)  :: T,rho_cgs
 real, intent(out) :: Q_cgs,dlnQ_dlnT

 real, parameter :: b1 = 1.57e-27, b2 = 1.25e-09, p = 1.2, q = 1.85
 real, parameter :: A = 1.28e-19, B = 1.71e-06, C = 3.02e-49
 real, parameter :: lambda = 2.4e-27, T_Compton = 2e+07
 real            :: mu, ne ! mean molecular weight & LTE electron density

 if (T > 1.1d4) then
    ! solve Saha equations for H, H2 and He to get the mean molecular weight
    call nelectron_mu(T, rho_cgs, 0., 0., ne, mu)
    ! Mathews & Doane 1990 (eq. 1)
    Q_cgs = (A * T**(-0.8) * max(0., 1. - B * T) - b1 * T**p / (1. + b2 * T**q) &
         - lambda * T**(1./2.) + C * (T_Compton - T) / rho_cgs) &
         * rho_cgs / (mu * mass_proton_cgs)**2
    ! could compute analytical formula but not needed for implicit cooling
    dlnQ_dlnT = 0.
 else
    Q_cgs = 0.
    dlnQ_dlnT = 0.
 endif

end subroutine cooling_high_temp

!-----------------------------------------------------------------------
!+
!  Cooling by H2 molecules
!
! :References:
!   Groenenwegen (1994), A&A 290, 531
!+
!-----------------------------------------------------------------------
subroutine cooling_H2(T, rho_cgs, Q_cgs, dlnQ_dlnT)

 use physcon, only:mass_proton_cgs

 real, intent(in)  :: T, rho_cgs
 real, intent(out) :: Q_cgs, dlnQ_dlnT

 real, parameter   :: f = 1.0d0
 real              :: nH2

 if (T < 900.) then
    nH2 = 0.5 * rho_cgs/(1.4*mass_proton_cgs)
    Q_cgs = -f*2.61111e-21 * nH2 * (T/1000.)**(4.74) / rho_cgs
    dlnQ_dlnT = 4.74
 else
    Q_cgs = 0.
    dlnQ_dlnT = 0.
 endif

end subroutine cooling_H2

!-----------------------------------------------------------------------
!+
!  Cooling by CO rotational lines in AGB outflows
!
!  log10(Lambda) = a + sum_i b_i log10(T)^i + sum_i c_i log10(n_H2)^i
!                    + sum_i d_i log10(n_CO/n_H2)^i + sum_i e_i log10(n_CO/|div v|)^i
!  with Lambda in W/m^3 and all inputs in SI. Inputs are clamped to the
!  validity range of the fit, and the cooling is switched off below Tmin
!
! :References:
!   Ceulemans, De Ceuster, Vermeulen & Decin (2026), RASTI 5, rzag009
!+
!-----------------------------------------------------------------------
subroutine cooling_CO_rot(T, rho_cgs, divv_cgs, Q_cgs, dlnQ_dlnT)

 use physcon, only:mass_proton_cgs

 real, intent(in)  :: T, rho_cgs, divv_cgs   ! divv in s^-1
 real, intent(out) :: Q_cgs, dlnQ_dlnT

 real, parameter :: a    = -17.424938947711073
 real, parameter :: b(3) = [1.08678135, -0.06837165, 0.0308211]
 real, parameter :: c(4) = [-6.13253761, 1.01649680, -5.32517202e-2, 9.04689294e-4]
 real, parameter :: d(1) = [0.9611932140011241]
 real, parameter :: e(2) = [0.22329612, -0.00630373]
 ! validity range of the fit (SI)
 real, parameter :: Tmin = 14., Tmax = 3000.
 real, parameter :: nH2_min = 1.4e8, nH2_max = 1.9e15
 real, parameter :: xCO_min = 1.e-4, xCO_max = 1.e-3
 real, parameter :: ndv_min = 1.5e13, ndv_max = 4.2e23
 real :: nH2, logT, lognH2, logxCO, logndv, logL
 integer :: i

 Q_cgs     = 0.
 dlnQ_dlnT = 0.
 if (T < Tmin) return

 nH2    = 0.5*rho_cgs/(1.4*mass_proton_cgs)*1.e6   ! m^-3
 logT   = log10(min(T, Tmax))
 lognH2 = log10(min(max(nH2, nH2_min), nH2_max))
 logxCO = log10(min(max(CO_abun, xCO_min), xCO_max))
 if (abs(divv_cgs) > tiny(0.)) then
    logndv = log10(min(max(CO_abun*nH2/abs(divv_cgs), ndv_min), ndv_max))
 else
    logndv = log10(ndv_max)
 endif

 logL = a
 do i = 1,size(b)
    logL = logL + b(i)*logT**i
 enddo
 do i = 1,size(c)
    logL = logL + c(i)*lognH2**i
 enddo
 do i = 1,size(d)
    logL = logL + d(i)*logxCO**i
 enddo
 do i = 1,size(e)
    logL = logL + e(i)*logndv**i
 enddo

 ! W/m^3 -> erg/s/cm^3 (x10), then per unit mass
 Q_cgs = -10.*10.**logL/rho_cgs
 if (T <= Tmax) then
    do i = 1,size(b)
       dlnQ_dlnT = dlnQ_dlnT + i*b(i)*logT**(i-1)
    enddo
 endif

end subroutine cooling_CO_rot

!-----------------------------------------------------------------------
!+
!  Number densities of the collision partners used in the Decin et al.
!  (2006) terms. The atomic fraction of H nuclei is recovered from mu
!  assuming a neutral gas, and rho = m_H n_H (1+4 f_He) as in Decin et al.
!  The fractions f_H = n(H)/n(H2) of the paper are avoided by working
!  with n(H) and n(H2) directly, e.g. n(H2)(f_H+2) = n_H
!+
!-----------------------------------------------------------------------
subroutine decin_densities(rho_cgs, mu, fHe, nH, nHI, nH2, nHe)
 use physcon, only:mass_proton_cgs
 real, intent(in)  :: rho_cgs, mu, fHe
 real, intent(out) :: nH, nHI, nH2, nHe
 real :: y

 nH  = rho_cgs/(mass_proton_cgs*(1.+4.*fHe))
 nHe = fHe*nH
 ! (1+4 f_He)/mu = y + f_He + (1-y)/2, with y = n(H)/n_H
 y   = min(max(2.*((1.+4.*fHe)/mu - fHe) - 1., 0.), 1.)
 nHI = y*nH
 nH2 = 0.5*(1.-y)*nH

end subroutine decin_densities

!-----------------------------------------------------------------------
!+
!  Drift velocity of a grain of size a (Decin et al. 2006, Eq. 4)
!
!  v_K^2 = v Q(a) L/(Mdot c) with Mdot = 4 pi r^2 rho v. The flux-mean
!  efficiency is taken in the Rayleigh regime, Q(a) proportional to a,
!  and normalised so that the MRN distribution exerts the radiative
!  acceleration arad that phantom applies to the gas, which gives
!  v_K^2 = 4 rho_s a arad / (3 psi rho)
!+
!-----------------------------------------------------------------------
real function drift_velocity_decin(a, T, mu, rho_cgs, arad)
 use physcon, only:kboltz,mass_proton_cgs
 real, intent(in) :: a, T, mu, rho_cgs, arad
 real :: vK2, vT, x

 vK2 = 4.*grain_rhos*a*arad/(3.*dust_to_gas*rho_cgs)
 if (vK2 <= 0.) then
    drift_velocity_decin = 0.
    return
 endif
 vT = 0.75*sqrt(3.*kboltz*T/(mu*mass_proton_cgs))
 x  = 0.5*vT**2/vK2
 ! sqrt(1+x^2)-x written in a form that does not cancel for large x
 drift_velocity_decin = sqrt(vK2/(sqrt(1.+x**2)+x))

end function drift_velocity_decin

!-----------------------------------------------------------------------
!+
!  A(r) n_H for the MRN size distribution n_d(a) = A a^-3.5 n_H, set
!  from the local dust-to-gas mass ratio psi
!+
!-----------------------------------------------------------------------
real function grain_norm_decin(rho_cgs)
 use physcon, only:pi
 real, intent(in) :: rho_cgs

 grain_norm_decin = 3.*dust_to_gas*rho_cgs/(8.*pi*grain_rhos*(sqrt(grain_amax)-sqrt(grain_amin)))

end function grain_norm_decin

!-----------------------------------------------------------------------
!+
!  Gas-grain collisional (drift) heating
!  H_gg = pi/2 A m_H n_H^2 (1+4f_He) int a^-1.5 v_drift^3 da
!  Grains drifting faster than 20 km/s are sputtered and do not contribute
!
! :References:
!   Decin et al. (2006), A&A 456, 549, Eqs. 4, 7, 8
!   Goldreich & Scoville (1976), ApJ 205, 144
!+
!-----------------------------------------------------------------------
subroutine heating_drift_decin(T, rho_cgs, mu, arad, Q_cgs, dlnQ_dlnT)
 use physcon, only:pi,kb_on_mh
 real, intent(in)  :: T, rho_cgs, mu, arad
 real, intent(out) :: Q_cgs, dlnQ_dlnT
 real    :: dlna, a, vd, w, x, vT2, vK2, s, ds
 integer :: k

 Q_cgs     = 0.
 dlnQ_dlnT = 0.
 if (arad <= 0. .or. dust_to_gas <= 0.) return

 vT2  = 0.75**2*3.*kb_on_mh*T/mu
 dlna = log(grain_amax/grain_amin)/(ngrain-1)
 s    = 0.
 ds   = 0.
 do k = 1,ngrain
    a  = grain_amin*exp((k-1)*dlna)
    vd = drift_velocity_decin(a, T, mu, rho_cgs, arad)
    if (vd > v_sputter .or. vd <= 0.) cycle
    w  = 1.
    if (k == 1 .or. k == ngrain) w = 0.5
    ! integrate in ln a: a^-1.5 v^3 da = a^-0.5 v^3 dln a
    s  = s + w*vd**3/sqrt(a)
    vK2 = 4.*grain_rhos*a*arad/(3.*dust_to_gas*rho_cgs)
    x  = 0.5*vT2/vK2
    ! dln v_drift/dln T = -x/(2 sqrt(1+x^2))
    ds = ds - 1.5*w*vd**3/sqrt(a)*x/sqrt(1.+x**2)
 enddo
 if (s <= 0.) return

 ! H_gg/rho, with rho = m_H n_H (1+4f_He)
 Q_cgs     = 0.5*pi*grain_norm_decin(rho_cgs)*s*dlna
 dlnQ_dlnT = ds/s

end subroutine heating_drift_decin

!-----------------------------------------------------------------------
!+
!  Heat exchange between dust and gas
!  H = 4 pi k sqrt(8k/(pi m_H)) A n_H^2 alpha_T sqrt(T) (T_d-T) int a^-1.5 da
!  Sputtered grains (v_drift > 20 km/s) do not contribute
!
! :References:
!   Decin et al. (2006), A&A 456, 549, Eqs. 12, 13
!   Burke & Hollenbach (1983), ApJ 265, 223
!+
!-----------------------------------------------------------------------
subroutine heating_dust_gas_decin(T, Tdust, rho_cgs, mu, nH, arad, Q_cgs, dlnQ_dlnT)
 use physcon, only:pi,kboltz,mass_proton_cgs
 real, intent(in)  :: T, Tdust, rho_cgs, mu, nH, arad
 real, intent(out) :: Q_cgs, dlnQ_dlnT
 real    :: dlna, a, w, s, alphaT
 integer :: k

 Q_cgs     = 0.
 dlnQ_dlnT = 0.
 if (dust_to_gas <= 0. .or. abs(Tdust-T) < tiny(T)) return

 dlna = log(grain_amax/grain_amin)/(ngrain-1)
 s    = 0.
 do k = 1,ngrain
    a = grain_amin*exp((k-1)*dlna)
    if (arad > 0.) then
       if (drift_velocity_decin(a, T, mu, rho_cgs, arad) > v_sputter) cycle
    endif
    w = 1.
    if (k == 1 .or. k == ngrain) w = 0.5
    s = s + w/sqrt(a)    ! a^-1.5 da = a^-0.5 dln a
 enddo
 s = s*dlna

 alphaT = 0.35*exp(-sqrt((Tdust+T)/500.)) + 0.1
 Q_cgs  = pi*grain_norm_decin(rho_cgs)*nH*2.*kboltz*sqrt(8.*kboltz*T/(pi*mass_proton_cgs)) &
          *alphaT*(Tdust-T)*s/rho_cgs
 dlnQ_dlnT = 0.5 - T/(Tdust-T)

end subroutine heating_dust_gas_decin

!-----------------------------------------------------------------------
!+
!  Heating by cosmic rays
!  H_cr = 6.4e-28 n(H2) (1+f_H/2) (1+4f_He)  erg/s/cm^3
!
! :References:
!   Decin et al. (2006), A&A 456, 549, Eq. 15
!   Goldsmith & Langer (1978), ApJ 222, 881
!+
!-----------------------------------------------------------------------
subroutine heating_cosmic_rays_decin(rho_cgs, nH2, nHI, fHe, Q_cgs, dlnQ_dlnT)
 real, intent(in)  :: rho_cgs, nH2, nHI, fHe
 real, intent(out) :: Q_cgs, dlnQ_dlnT

 Q_cgs     = 6.4e-28*(nH2+0.5*nHI)*(1.+4.*fHe)/rho_cgs
 dlnQ_dlnT = 0.

end subroutine heating_cosmic_rays_decin

!-----------------------------------------------------------------------
!+
!  Photoelectric heating from dust grains (Bakes & Tielens 1994 scaled
!  by 0.2 for a_min = 50 A), attenuated by the circumstellar UV extinction
!  tau_uv = 1.8 A_v, A_v = 1.6e-22/0.01 psi N_H
!
!  Electrons come from CO -> C + O -> C+ + e, so n_e is the dissociated
!  part of the CO abundance (Mamon et al. 1988 profile). The radial column
!  from the outer edge is approximated by N_H = n_H r (rho ~ r^-2)
!
! :References:
!   Decin et al. (2006), A&A 456, 549, Eqs. 9, 10, 22
!+
!-----------------------------------------------------------------------
subroutine heating_photoelectric_decin(T, rho_cgs, nH, nH2, r, Q_cgs, dlnQ_dlnT)
 real, intent(in)  :: T, rho_cgs, nH, nH2, r
 real, intent(out) :: Q_cgs, dlnQ_dlnT
 real :: ne, Av, xx

 Q_cgs     = 0.
 dlnQ_dlnT = 0.
 if (r_half_CO <= 0. .or. r <= 0.) return

 ne = CO_abun*nH2*(1.-exp(-log(2.)*(r/r_half_CO)**alpha_CO))
 if (ne <= 0.) return
 Av = 1.6e-22/0.01*dust_to_gas*nH*r
 xx = 2.e-4*G0_UV*sqrt(T)/ne
 Q_cgs     = 0.2e-24*3.e-2/(1.+xx)*nH*G0_UV*exp(-1.8*Av)/rho_cgs
 dlnQ_dlnT = -0.5*xx/(1.+xx)

end subroutine heating_photoelectric_decin

!-----------------------------------------------------------------------
!+
!  Rotational line cooling in the three-level approximation of
!  Goldreich & Scoville (1976), as used by Justtanont et al. (1994)
!
!  cooling = n_mol C h nu21 [exp(-h nu21/kT) - exp(-h nu21/kTx)]
!  with nu21 = nu0 Tx^0.5, and the excitation temperature Tx solving
!  C h nu21/k (1/Tx-1/T) = beta21 A21
!                          + eps W A31/2 h nu21/k exp(-h nu31/kT*) (1/T*-1/Tx)
!  beta21 A21 = betacoef/n_mol (v/r)(1+eps/2) Tx^pexp
!  For a spherical outflow (v/r)(1+eps/2) = div v/2
!+
!-----------------------------------------------------------------------
subroutine cooling_rot_J94(T, rho_cgs, Ccoll, nmol, nu0, nu31, A31, betacoef, pexp, &
                           divv_cgs, eps, W, Tstar, Q_cgs, dlnQ_dlnT)
 use physcon, only:planckh,kboltz
 real, intent(in)  :: T, rho_cgs, Ccoll, nmol, nu0, nu31, A31, betacoef, pexp
 real, intent(in)  :: divv_cgs, eps, W, Tstar
 real, intent(out) :: Q_cgs, dlnQ_dlnT
 integer, parameter :: itermax = 60
 real,    parameter :: tol = 1.e-7
 real    :: theta0, B, P, Tx, Tlo, Thi, g, dg, Tnew, E, eT, eTx
 integer :: iter

 Q_cgs     = 0.
 dlnQ_dlnT = 0.
 if (nmol <= 0. .or. Ccoll <= 0.) return

 theta0 = planckh*nu0/kboltz
 B = betacoef/nmol*0.5*abs(divv_cgs)
 P = 0.
 if (Tstar > 0. .and. eps > 0.) P = 0.5*eps*W*A31*exp(-planckh*nu31/(kboltz*Tstar))
 if (B <= 0. .and. P <= 0.) return   ! lines fully trapped: Tx = T

 ! g(Tx) decreases monotonically, g(0+) > 0 and g(max(T,T*)) <= 0
 Tlo = 1.e-3
 Thi = max(T, Tstar)
 Tx  = T
 do iter = 1,itermax
    g  = Ccoll*theta0*(1./sqrt(Tx) - sqrt(Tx)/T) - B*Tx**pexp
    dg = -0.5*Ccoll*theta0*(Tx**(-1.5) + 1./(sqrt(Tx)*T)) - pexp*B*Tx**(pexp-1.)
    if (P > 0.) then
       g  = g  - P*theta0*(sqrt(Tx)/Tstar - 1./sqrt(Tx))
       dg = dg - 0.5*P*theta0*(1./(sqrt(Tx)*Tstar) + Tx**(-1.5))
    endif
    if (g > 0.) then
       Tlo = Tx
    else
       Thi = Tx
    endif
    Tnew = Tx - g/dg
    if (Tnew <= Tlo .or. Tnew >= Thi) Tnew = sqrt(Tlo*Thi)
    if (abs(Tnew-Tx) < tol*Tx) exit
    Tx = Tnew
 enddo
 Tx = Tnew

 E   = theta0*sqrt(Tx)          ! h nu21/k
 eT  = exp(-E/T)
 eTx = exp(-E/Tx)
 Q_cgs = -nmol*Ccoll*kboltz*E*(eT-eTx)/rho_cgs
 ! at fixed Tx: C ~ T^0.5
 if (abs(eT-eTx) > tiny(eT)) dlnQ_dlnT = 0.5 + E/T*eT/(eT-eTx)

end subroutine cooling_rot_J94

!-----------------------------------------------------------------------
!+
!  Cooling by rotational excitation of H2O
!  <sigma v>(He-H2O) = 0.21e-11 T^0.5, H a factor 1.16 larger, H2 from
!  Phillips et al. (1996) ratios, ortho:para = 3:1
!
! :References:
!   Decin et al. (2006), A&A 456, 549, Eq. 17
!   Justtanont et al. (1994), ApJ 435, 852, Eqs. 11, 12
!+
!-----------------------------------------------------------------------
subroutine cooling_H2O_rot_decin(T, rho_cgs, nHI, nH2, nHe, divv_cgs, eps, W, Tstar, &
                                 Q_cgs, dlnQ_dlnT)
 real, intent(in)  :: T, rho_cgs, nHI, nH2, nHe, divv_cgs, eps, W, Tstar
 real, intent(out) :: Q_cgs, dlnQ_dlnT
 real :: Ccoll

 Ccoll = 0.21e-11*sqrt(T)*(0.83*nHI + 0.715*nHe + 4.5*nH2)
 call cooling_rot_J94(T, rho_cgs, Ccoll, H2O_abun*nH2, 2.6e11, 1.13e14, 34., 6.57, 3., &
                      divv_cgs, eps, W, Tstar, Q_cgs, dlnQ_dlnT)

end subroutine cooling_H2O_rot_decin

!-----------------------------------------------------------------------
!+
!  Cooling by rotational excitation of CO
!  <sigma v>(CO-H2) = 4e-12 T^0.5, <sigma v>(CO-He) = 4.5e-13 T^0.5,
!  H a factor 1.16 larger than He. CO is photodissociated following
!  Mamon et al. (1988) when r_half_CO > 0
!
! :References:
!   Decin et al. (2006), A&A 456, 549, Eqs. 21, 22
!   Justtanont et al. (1994), ApJ 435, 852, Eqs. 16, 17
!+
!-----------------------------------------------------------------------
subroutine cooling_CO_rot_decin(T, rho_cgs, nHI, nH2, nHe, r, divv_cgs, eps, W, Tstar, &
                                Q_cgs, dlnQ_dlnT)
 real, intent(in)  :: T, rho_cgs, nHI, nH2, nHe, r, divv_cgs, eps, W, Tstar
 real, intent(out) :: Q_cgs, dlnQ_dlnT
 real :: Ccoll, nCO

 nCO = CO_abun*nH2
 if (r_half_CO > 0.) nCO = nCO*exp(-log(2.)*(r/r_half_CO)**alpha_CO)
 Ccoll = 4.5e-13*sqrt(T)*(1.16*nHI + nHe + 10.*nH2)
 call cooling_rot_J94(T, rho_cgs, Ccoll, nCO, 4.89e10, 6.43e13, 34., 6.60, 2.5, &
                      divv_cgs, eps, W, Tstar, Q_cgs, dlnQ_dlnT)

end subroutine cooling_CO_rot_decin

!-----------------------------------------------------------------------
!+
!  Cooling by vibrational excitation of H2 (v=1 -> 0)
!
! :References:
!   Decin et al. (2006), A&A 456, 549, Eqs. 23-27
!   Hollenbach & McKee (1979, 1989)
!+
!-----------------------------------------------------------------------
subroutine cooling_H2_vib_decin(T, rho_cgs, nHI, nH2, Q_cgs, dlnQ_dlnT)
 use physcon, only:eV,kboltz
 real, intent(in)  :: T, rho_cgs, nHI, nH2
 real, intent(out) :: Q_cgs, dlnQ_dlnT
 real, parameter :: A10 = 3.e-7, E10 = 0.6*eV
 real :: e, alpha_n, n1

 Q_cgs     = 0.
 dlnQ_dlnT = 0.
 if (nH2 <= 0.) return

 e       = exp(-E10/(kboltz*T))
 alpha_n = nHI*1.0e-12*sqrt(T)*exp(-1000./T) + nH2*1.4e-12*sqrt(T)*exp(-18100./(T+1200.))
 n1      = nH2*alpha_n*e/(alpha_n*(1.+e) + A10)
 Q_cgs   = -A10*E10*n1/rho_cgs
 ! dominant (Boltzmann factor) part of the temperature dependence
 dlnQ_dlnT = E10/(kboltz*T)*(alpha_n+A10)/(alpha_n*(1.+e)+A10)

end subroutine cooling_H2_vib_decin

!-----------------------------------------------------------------------
!+
!  compute electron equilibrium abundance per nH atom (Palla et al 1983)
!+
!-----------------------------------------------------------------------
real function calc_eps_e(T)

 real, intent(in) :: T

 real             :: k1, k2, k3, k8, k9, p, q

 k1 = 1.88d-10 / T**6.44e-1
 k2 = 1.83d-18 * T
 k3 = 1.35d-9
 k8 = 5.80d-11 * sqrt(T) * exp(-1.58d5/T)
 k9 = 1.7d-4 * k8
 p  = .5*k8/k9
 q  = k1*(k2+k3)/(k3*k9)
 calc_eps_e = min(1.,(p + sqrt(q+p**2))/q) !must be <= 1

end function calc_eps_e

!-----------------------------------------------------------------------
!+
!  cooling functions with analytical solutions (for analysis)
!+
!-----------------------------------------------------------------------
subroutine testing_cooling_functions(ifunct, T, Q, dlnQ_dlnT)

 integer, intent(in)  :: ifunct
 real,    intent(in)  :: T
 real,    intent(out) :: Q,dlnQ_dlnT

 select case(ifunct)
 case (0)
    !test1 : du/dt = cst --> linear decrease in time
    Q = -1e13
    dlnQ_dlnT = 0.
 case(1)
    !test2 : du/dt = -a*u --> exponential decrease in time
    Q = -1e7*T
    dlnQ_dlnT = 1.
 case(3)
    !test3 : du/dt = -a*u**3 --> powerlaw decrease in time**(-1/2)
    Q = -1e-5*T**3
    dlnQ_dlnT = 3.
 case(-3)
    !test4 : du/dt = -a*u**-3 --> powerlaw decrease in time**(1/4)
    Q = -1e20/T**3
    dlnQ_dlnT = -3.
 case default
    Q = 0.
    dlnQ_dlnT = 0.
 end select

end subroutine testing_cooling_functions

!-----------------------------------------------------------------------
!+
!  ADDITIONAL PHYSICS: compute LTE electron density from SAHA equations
!                      (following D'Angelo & Bodenheimer 2013)
!+
!-----------------------------------------------------------------------
subroutine nelectron_mu(T_gas, rho_gas, nH, nHe, n_e, mu)

 use physcon, only:kboltz, mass_proton_cgs, mass_electron_cgs, planckhbar, pi

 real, intent(in)  :: T_gas, rho_gas, nH, nHe
 real, intent(out) :: n_e
 real, intent(out), optional :: mu

 real, parameter  :: H2_diss = 7.178d-12    !  4.48 eV in erg
 real, parameter  :: H_ion   = 2.179d-11    ! 13.60 eV in erg
 real, parameter  :: He_ion  = 3.940d-11    ! 24.59 eV in erg
 real, parameter  :: He2_ion = 8.720d-11    ! 54.42 eV in erg
 real             :: KH, KH2, xx, yy, KHe, KHe2, z1, z2, cst

 cst = mass_proton_cgs/rho_gas*sqrt(mass_electron_cgs*kboltz*T_gas/(2.*pi*planckhbar**2))**3
 if (T_gas > 1.d5) then
    ! all hydrogen is ionized
    xx = 1.
 else
    KH   = cst/xH * exp(-H_ion /(kboltz*T_gas))
    ! solution to quadratic SAHA equations (Eq. 16 in D'Angelo et al 2013)
    xx   = 0.5 * (-KH + sqrt(KH**2+4.*KH))
 endif

 if (T_gas > 1.d4) then
    ! all H2 has been dissociated
    yy = 1.
 else
    KH2  = 0.5*cst/xH * sqrt(0.5*mass_proton_cgs/mass_electron_cgs)**3 * exp(-H2_diss/(kboltz*T_gas))
    ! solution to quadratic SAHA equations (Eq. 15 in D'Angelo et al 2013)
    yy   = 0.5 * (-KH2 + sqrt(KH2**2+4.*KH2))
 endif

 if (T_gas > 4.d5) then
    ! all helium has been ionized twice
    z1 = 1.
    z2 = 1.
 else
    KHe    = 4.*cst * exp(-He_ion/(kboltz*T_gas))
    KHe2   =    cst * exp(-He2_ion/(kboltz*T_gas))
    ! solution to quadratic SAHA equations (Eq. 17 in D'Angelo et al 2013)
    z1     = (2./XHe) * (-KHe-xH + sqrt((KHe+xH)**2+KHe*xHe))
    ! solution to quadratic SAHA equations (Eq. 18 in D'Angelo et al 2013)
    z2     = (2./xHe) * (-KHe2-xH-xHe/4. + sqrt((KHe2+xH+xHe/4.)**2+KHe2*xHe))
 endif

 n_e       = xx * nH + z1*(1.+z2) * nHe
 mu        = 4./(2.*xH*(1.+yy+2.*xx*yy)+xHe*(1+z1+z1*z2))

end subroutine nelectron_mu

!-----------------------------------------------------------------------
!+
!  ADDITIONAL PHYSICS: compute mean thermal speed of molecules
!+
!-----------------------------------------------------------------------
real function v_th(T_gas,mu)

 use physcon, only:kboltz, mass_proton_cgs

 real, intent(in) :: T_gas, mu

 v_th = sqrt((3.*kboltz*T_gas)/(mu*mass_proton_cgs))

end function v_th

!-----------------------------------------------------------------------
!+
!  ADDITIONAL PHYSICS: compute fraction of gas that has speeds lower than v_crit
!                      from the cumulative distribution function of the
!                      Maxwell-Boltzmann distribution
! doi : 10.4236/ijaa.2020.103010
!-----------------------------------------------------------------------
real function MaxBol_cumul(T_gas, mu,  v_crit)

 use physcon, only:kboltz, mass_proton_cgs, pi

 real, intent(in) :: T_gas, mu, v_crit

 real             :: a

 a            = sqrt(2.*kboltz*T_gas/(mu*mass_proton_cgs))
 MaxBol_cumul = erf(v_crit/a) - 2./sqrt(pi) * v_crit/a *exp(-(v_crit/a)**2)

end function MaxBol_cumul

!-----------------------------------------------------------------------
!+
!  ADDITIONAL PHYSICS: compute dust number density from dust-to-gas mass ratio,
!                      mean grain size a, and specific density of the grain
!+
!-----------------------------------------------------------------------
real function n_dust(rho_gas, d2g, a, rho_grain)

 use physcon, only:pi

 real, intent(in) :: rho_gas,d2g,a,rho_grain

 n_dust = ( rho_gas*d2g ) / ( (4./3.)*pi*a**3.*rho_grain )

end function n_dust

!=======================================================================
!=======================================================================
!=======================================================================
!
!  Cooling functions    **** ALL IN cgs  ****
!
!=======================================================================

!-----------------------------------------------------------------------
!+
!  DUST:  Full contact cooling (Bowen 1988)
!+
!-----------------------------------------------------------------------
real function cool_dust_full_contact(T_gas, rho_gas, mu, T_dust, kappa_dust)

 use physcon, only:Rg

 real, intent(in) :: T_gas, rho_gas, mu
 real, intent(in) :: T_dust, kappa_dust

 if (kappa_dust > kappa_dust_min) then
    cool_dust_full_contact = (3.*Rg)/(2.*mu*bowen_Cprime)*rho_gas*(T_gas-T_dust)
 else
    cool_dust_full_contact = 0.0
 endif
end function cool_dust_full_contact

!-----------------------------------------------------------------------
!+
!  DUST: Discrete contact cooling (Hollenbach & McKee 1979)
!+
!-----------------------------------------------------------------------
real function cool_dust_discrete_contact(T_gas, rho_gas, mu, T_dust, d2g, a, rho_grain, kappa_dust)

 use physcon, only:kboltz, mass_proton_cgs, pi

 real, intent(in) :: T_gas, rho_gas, mu
 real, intent(in) :: T_dust, d2g, a, rho_grain, kappa_dust

 real, parameter   :: alpha = 0.33  ! See Burke & Hollenbach 1983
 real              :: n_gas, sigma_dust

 if (kappa_dust > kappa_dust_min) then
    sigma_dust                 = 2.*pi*a**2
    n_gas                      = rho_gas/(mu*mass_proton_cgs)
    cool_dust_discrete_contact = alpha*n_gas*n_dust(rho_gas,d2g,a,rho_grain)*sigma_dust*v_th(T_gas,mu)*kboltz*(T_gas-T_dust)
 else
    cool_dust_discrete_contact = 0.0
 endif
end function cool_dust_discrete_contact

!-----------------------------------------------------------------------
!+
!  DUST: Radiative cooling (Woitke 2006) - DO NOT USE, PHYSICALLY INCORRECT
!+
!-----------------------------------------------------------------------
real function cool_dust_radiation(T_gas, kappa_gas, T_dust, kappa_dust)

 use physcon, only:steboltz

 real, intent(in) :: T_gas, kappa_gas
 real, intent(in) :: T_dust, kappa_dust

 if (kappa_dust > kappa_dust_min) then
    cool_dust_radiation = 4.*steboltz*(kappa_gas*T_gas**4-kappa_dust*T_dust**4)
 else
    cool_dust_radiation = 0.0
 endif
end function cool_dust_radiation

!-----------------------------------------------------------------------
!+
!  DUST: Friction heating caused by dust-drift through gas (Golreich & Scoville 1976)
!+
!-----------------------------------------------------------------------
real function heat_dust_friction(rho_gas, v_drift, d2g, a, rho_grain, kappa_dust)
 use physcon, only:pi
 real, intent(in) :: rho_gas
 real, intent(in) :: v_drift, d2g, a, rho_grain, kappa_dust

 real              :: sigma_dust
 real, parameter   :: alpha = 0.33                            ! see Burke & Hollenbach 1983

 ! Warning, alpha depends on the type of dust
 if (kappa_dust > kappa_dust_min) then
    sigma_dust         = 2.*pi*a**2
    heat_dust_friction = n_dust(rho_gas,d2g,a,rho_grain)*sigma_dust*v_drift*alpha*0.5*rho_gas*v_drift**2
 else
    heat_dust_friction = 0.0
 endif

end function heat_dust_friction

!-----------------------------------------------------------------------
!+
!  DUST: photovoltaic heating by soft UV field (Weingartner & Draine 2001)
!+
!-----------------------------------------------------------------------
real function heat_dust_photovoltaic_soft(T_gas, rho_gas, mu, nH, nHe, kappa_dust)

 real, intent(in) :: T_gas, rho_gas, mu, nH, nHe
 real, intent(in) :: kappa_dust

 real              :: x,n_e
 real, parameter   :: G=1.68 ! ratio of true background UV field to Habing field
 real, parameter   :: C0=5.45, C1=2.50, C2=0.00945, C3=0.01453, C4=0.147, C5=0.623, C6=0.511 ! see Table 2 in Weingartner & Draine 2001, last line

 if (kappa_dust > kappa_dust_min) then
    call nelectron_mu(T_gas, rho_gas, nH, nHe, n_e)
    x = G*sqrt(T_gas)/n_e
    heat_dust_photovoltaic_soft = 1.d-26*G*nH*(C0+C1*T_gas**C4)/(1.+C2*x**C5*(1.+C3*x**C6))
 else
    heat_dust_photovoltaic_soft = 0.0
 endif

end function heat_dust_photovoltaic_soft

!-----------------------------------------------------------------------
!+
!  DUST: photovoltaic heating by hard UV field (Inoue & Kamaya 2010)
!+
!-----------------------------------------------------------------------
real function heat_dust_photovoltaic_hard(T_gas, nH, d2g, kappa_dust, JL)

 real, intent(in) :: T_gas, nH
 real, intent(in) :: d2g, kappa_dust
 real, intent(in) :: JL       ! mean intensity of background UV radiation at hydrogen Lyman limit (91.2 nm)

 if (kappa_dust > kappa_dust_min) then
    heat_dust_photovoltaic_hard = 1.2d-34*(d2g   /1.d-4) &
                                         *(nH   /1.d-5 )**(4. /3.) &
                                         *(T_gas/1.d4  )**(-1./6.) &
                                         *(JL   /1.d-21)**(2. /3.)
 else
    heat_dust_photovoltaic_hard = 0.0
 endif

end function heat_dust_photovoltaic_hard

!-----------------------------------------------------------------------
!+
!  PARTICLE: Coulomb cooling via electron scattering (Weingartner & Draine 2001)
!+
!-----------------------------------------------------------------------
real function cool_coulomb(T_gas, rho_gas, mu, nH, nHe)

 real, intent(in) :: T_gas, rho_gas, mu, nH, nHe

 real              :: x, n_e
 real, parameter   :: G=1.68 ! ratio of true background UV field to Habing field
 real, parameter   :: D0=0.4255, D1=2.457, D2=-6.404, D3=1.513, D4=0.05343 ! see Table 3 in Weingartner & Draine 2001, last line

 if (T_gas > 1000.) then !. .and. T_gas < 1.e4) then
    call nelectron_mu(T_gas, rho_gas, nH, nHe, n_e)
    x  = log(G*sqrt(T_gas)/n_e)
    cool_coulomb = 1.d-28*n_e*nH*T_gas**(D0+D1/x)*exp(D2+D3*x-D4*x**2)
 else
    cool_coulomb = 0.0
 endif

end function cool_coulomb

!-----------------------------------------------------------------------
!+
!  PARTICLE: Cosmic ray heating (Jonkheid et al. 2004)
!+
!-----------------------------------------------------------------------
real function heat_CosmicRays(nH, nH2)

 real, intent(in) :: nH, nH2
 real, parameter  :: Rcr = 5.0d-17  !cosmic ray ionisation rate [s^-1]

 heat_CosmicRays = Rcr*(5.5d-12*nH+2.5d-11*nH2)

end function heat_CosmicRays

!-----------------------------------------------------------------------
!+
!  ATOMIC: Cooling due to electron excitation of neutral H (Spitzer 1978, Black 1982, Cen 1992)
!+
!-----------------------------------------------------------------------
real function cool_HI(T_gas, rho_gas, mu, nH, nHe)

 use physcon, only:mass_proton_cgs

 real, intent(in) :: T_gas, rho_gas, mu, nH, nHe
 real              :: n_gas,n_e

 ! all hydrogen atomic, so nH = n_gas
 ! Dalgarno & McCray (1972) provide data starting at 3000K
 ! (1+sqrt(T_gas/1.d5))**(-1) correction factor added by Cen 1992
 if (T_gas > 200000.) then
    n_gas   = rho_gas/(mu*mass_proton_cgs)
    !nH      = XH*n_gas
    call nelectron_mu(T_gas, rho_gas, nH, nHe, n_e)
    cool_HI = 7.3d-19*n_e*n_gas/(1.+sqrt(T_gas/1.d5))*exp(-118400./T_gas)
 else
    cool_HI = 0.0
 endif

end function cool_HI

!-----------------------------------------------------------------------
!+
!  ATOMIC: Cooling due to collisional ionisation of neutral H (Black 1982, Cen 1992)
!+
!-----------------------------------------------------------------------
real function cool_H_ionisation(T_gas, rho_gas, mu, nH, nHe)

 use physcon, only:mass_proton_cgs

 real, intent(in) :: T_gas, rho_gas, mu, nH, nHe
 real              :: n_gas, n_e

 ! all hydrogen atomic, so nH = n_gas
 ! (1+sqrt(T_gas/1.d5))**(-1) correction factor added by Cen 1992
 if (T_gas > 4000.) then
    n_gas   = rho_gas/(mu*mass_proton_cgs)
    !nH      = XH*n_gas
    call nelectron_mu(T_gas, rho_gas, nH, nHe, n_e)
    cool_H_ionisation = 1.27d-21*n_e*n_gas*sqrt(T_gas)/(1.+sqrt(T_gas/1.d5))*exp(-157809./T_gas)
 else
    cool_H_ionisation = 0.0
 endif

end function cool_H_ionisation

!-----------------------------------------------------------------------
!+
!  ATOMIC: Cooling due to collisional ionisation of neutral He (Black 1982, Cen 1992)
!+
!-----------------------------------------------------------------------
real function cool_He_ionisation(T_gas, rho_gas, mu, nH, nHe)

 use physcon, only:mass_proton_cgs

 real, intent(in) :: T_gas, rho_gas, mu, nH, nHe
 real              :: n_gas, n_e

 ! all hydrogen atomic, so nH = n_gas
 ! (1+sqrt(T_gas/1.d5))**(-1) correction factor added by Cen 1992
 if (T_gas > 4000.) then
    n_gas   = rho_gas/(mu*mass_proton_cgs)
    !nH      = XH*n_gas
    call nelectron_mu(T_gas, rho_gas, nH, nHe, n_e)
    cool_He_ionisation = 9.38d-22*n_e*nHe*sqrt(T_gas)*(1+sqrt(T_gas/1.d5))**(-1)*exp(-285335./T_gas)
 else
    cool_He_ionisation = 0.0
 endif

end function cool_He_ionisation

!-----------------------------------------------------------------------
!+
!  CHEMICAL: Cooling due to ro-vibrational excitation of H2 (Lepp & Shull 1983)
!            (Smith & Rosen, 2003, MNRAS, 339)
!+
!-----------------------------------------------------------------------
real function cool_H2_rovib(T_gas, nH, nH2)

 real, intent(in) :: T_gas, nH, nH2
 real              :: kH_01, kH2_01
 real              :: Lvh, Lvl, Lrh, Lrl
 real              :: x, Qn

 if (T_gas < 1635.) then
    kH_01 = 1.4d-13*exp((T_gas/125.)-(T_gas/577.)**2)
 else
    kH_01 = 1.0d-12*sqrt(T_gas)*exp(-1000./T_gas)
 endif
 kH2_01 = 1.45d-12*sqrt(T_gas)*exp(-28728./(T_gas+1190.))
 Lvh    = 1.1d-18*exp(-6744./T_gas)
 Lvl    = 8.18d-13*(nH*kH_01+nH2*kH2_01)*exp(-6840./T_gas)

 x   = log10(T_gas/1.0d4)
 if (T_gas < 1087.) then
    Lrh = 10.**(-19.24+0.474*x-1.247*x**2)
 else
    Lrh = 3.9d-19*exp(-6118./T_gas)
 endif

 Qn = nH2**0.77+1.2*nH**0.77
 if (T_gas > 4031.) then
    Lrl = 10.**(-22.9-0.553*x-1.148*x**2)*Qn
 else
    Lrl = 1.38d-22*exp(-9243./T_gas)*Qn
 endif

 cool_H2_rovib = nH2*( Lvh/(1.+(Lvh/Lvl)) + Lrh/(1.+(Lrh/Lrl)) )

end function cool_H2_rovib

!-----------------------------------------------------------------------
!+
!  CHEMICAL: H2 dissociation cooling (Shapiro & Kang 1987, Smith & Rosen 2003)
!+
!-----------------------------------------------------------------------
real function cool_H2_dissociation(T_gas, rho_gas, mu, nH, nH2)

 use physcon, only:mass_proton_cgs

 real, intent(in) :: T_gas, rho_gas, mu, nH, nH2

 real              :: n_gas
 real              :: x, n1, n2, beta
 real              :: kD_H, kD_H2

 n_gas = rho_gas/(mu*mass_proton_cgs)
 x     = log10(T_gas/1.0d4)
 n1    = 10.**(4.0   -0.416*x -0.327*x**2)
 n2    = 10.**(4.845 -1.3*x   +1.62*x**2)
 beta  = 1./(1.+n_gas*(2.*nH2/n_gas*((1./n2)-(1./n1))+1./n1))
 kD_H  = 1.2d-9*exp(-52400/T_gas)*(0.0933*exp(-17950./T_gas))**beta
 kD_H2 = 1.3d-9*exp(-53300/T_gas)*(0.0908*exp(-16200./T_gas))**beta

 cool_H2_dissociation = 7.18d-12*(nH2**2*kD_H2+nH*nH2*kD_H)

end function cool_H2_dissociation

!-----------------------------------------------------------------------
!+
!  CHEMICAL: H2 recombination heating (Hollenbach & Mckee 1979)
!            for an overview, see Wakelam et al. 2017, Smith & Rosen 2003
!+
!-----------------------------------------------------------------------
real function heat_H2_recombination(T_gas, rho_gas, mu, nH, nH2, T_dust)

 use physcon, only:mass_proton_cgs

 real, intent(in) :: T_gas, rho_gas, mu, nH, nH2, T_dust

 real              :: n_gas
 real              :: x, n1, n2, beta
 real              :: xi, fa, k_rec

 n_gas  = rho_gas/(mu*mass_proton_cgs)
 x      = log10(T_gas/1.0d4)
 n1     = 10.**(4.0   -0.416*x -0.327*x**2)
 n2     = 10.**(4.845 -1.3*x   +1.62*x**2)
 beta   = 1./(1.+n_gas*(2.*nH2/n_gas*((1./n2)-(1./n1))+1./n1))
 xi     = 7.18d-12*n_gas*nH*(1.-beta)

 fa     = 1./(1.+1.d4*exp(-600./T_dust))    ! eq 3.4
 k_rec  = 3.d-18*(sqrt(T_gas)*fa)/(1.+0.04*sqrt(T_gas+T_dust)+2.d-3*T_gas+8.d-6*T_gas**2) ! eq 3.8

 heat_H2_recombination = k_rec*xi

end function heat_H2_recombination

!-----------------------------------------------------------------------
!+
!  RADIATIVE: optically thin CO ro-vibrational cooling (Hollenbach & McKee 1979, McKee et al. 1982)
!+
!-----------------------------------------------------------------------
real function cool_CO_rovib(T_gas, rho_gas, mu, nH, nH2, nCO)

 use physcon, only:kboltz, mass_proton_cgs

 real, intent(in) :: T_gas, rho_gas, mu, nH, nH2, nCO

 real              :: Qrot, QvibH2, QvibH
 real              :: n_gas, n_crit, sigma
 real              :: v_crit, nfCO

! CO bond dissociation energy = 11.11 eV = 1.78e-11 erg
! use cumulative distribution of Maxwell-Boltzmann
! to account for collisions that destroy CO

 if (T_gas > 3000. .or. T_gas < 250.) then
    cool_CO_rovib = 0.
    return
 endif
 v_crit = sqrt( 2.*1.78d-11/(mu*mass_proton_cgs) )  ! kinetic energy
 nfCO   = MaxBol_cumul(T_gas, mu,  v_crit) * nCO

 n_gas  = rho_gas/(mu*mass_proton_cgs)
 n_crit = 3.3d6*(T_gas/1000.)**0.75      !McKee et al. 1982 eq. 5.3
 sigma  = 3.d-16*(T_gas/1000.)**(-0.25)  !McKee et al. 1982 eq. 5.4
 !v_th = sqrt((8.*kboltz*T_gas)/(pi*mH2_cgs)) !3.1
 Qrot   = 0.5*n_gas*nfCO*kboltz*T_gas*sigma*v_th(T_gas, mu) / (1. + (n_gas/n_crit) + 1.5*sqrt(n_gas/n_crit))
!McKee et al. 1982 eq. 5.2

 QvibH2 = 1.83d-26*nH2*nfCO*T_gas*exp(-3080./T_gas)*exp(-68./(T_gas**(1./3.)))  !Smith & Rosen
 QvibH  = 1.28d-24*nH *nfCO*sqrt(T_gas)*exp(-3080./T_gas)*exp(-(2000./T_gas)**3.43) !Smith & Rosen

 cool_CO_rovib = Qrot+QvibH+QvibH2

end function cool_CO_rovib

!-----------------------------------------------------------------------
!+
!  RADIATIVE: H20 ro-vibrational cooling (Hollenbach & McKee 1989, Neufeld & Kaufman 1993)
!+
!-----------------------------------------------------------------------
real function cool_H2O_rovib(T_gas, rho_gas, mu, nH, nH2, nH2O)

 use physcon, only:mass_proton_cgs

 real, intent(in) :: T_gas, rho_gas, mu, nH, nH2, nH2O

 real              :: Qrot, QvibH2, QvibH
 real              :: alpha, lambdaH2O
 real              :: v_crit, nfH2O

! Binding energy of singular O-H bond = 5.151 eV = 8.25e-12 erg
! use cumulative distribution of Maxwell-Boltzmann
! to account for collisions that destroy H2O

 v_crit = sqrt( 2.*8.25d-12/(mu*mass_proton_cgs) )  ! kinetic energy
 nfH2O  = MaxBol_cumul(T_gas, mu,  v_crit) * nH2O

 alpha     = 1.35 - 0.3*log10(T_gas/1000.)          ! Neufeld & Kaufmann 1993
 lambdaH2O = 1.32d-23*(T_gas/1000.)**alpha
 Qrot      = (nH2+1.39*nH)*nfH2O*LambdaH2O

 QvibH2    = 1.03d-26*nH2*nfH2O*T_gas*exp(-2352./T_gas)*exp(-47.5/(T_gas**(1./3.)))            !Hollenbach & McKee 1989 eq 2.14b
 QvibH     = 7.40d-27*nH *nfH2O*lambdaH2O*T_gas*exp(-2352./T_gas)*exp(-34.5/(T_gas**(1./3.)))  !Hollenbach & McKee 1989 eq 2.14a

 cool_H2O_rovib = Qrot+QvibH+QvibH2

end function cool_H2O_rovib

!-----------------------------------------------------------------------
!+
!  RADIATIVE: OH rotational cooling (Hollenbach & McKee 1979, McKee et al. 1982)
!+
!-----------------------------------------------------------------------
real function cool_OH_rot(T_gas, rho_gas, mu, nOH)
 use physcon, only:kboltz, mass_proton_cgs

 real, intent(in) :: T_gas, rho_gas, mu, nOH

 real              :: n_gas
 real              :: sigma, n_crit
 real              :: v_crit, nfOH

! Binding energy of singular O-H bond = 5.151 eV = 8.25e-12 erg
! use cumulative distribution of Maxwell-Boltzmann
! to account for collisions that destroy OH

 v_crit = sqrt( 2.*8.25d-12/(mu*mass_proton_cgs) )  ! kinetic energy
 nfOH   = MaxBol_cumul(T_gas, mu,  v_crit) * nOH

 n_gas     = rho_gas/(mu*mass_proton_cgs)
 sigma     = 2.0d-16
 !n_crit    = 1.33d7*sqrt(T_gas)
 n_crit    = 1.5d10*sqrt(T_gas/1000.) !table 3 Hollenbach & McKee 1989

 cool_OH_rot = n_gas*nfOH*(kboltz*T_gas*sigma*v_th(T_gas, mu)) / (1 + n_gas/n_crit + 1.5*sqrt(n_gas/n_crit))  !McKee et al. 1982 eq. 5.2

end function cool_OH_rot

!-----------------------------------------------------------------------
!+
!  RADIATIVE: heating due to recombination (Mathews & Doane 1990)
!+
!-----------------------------------------------------------------------
real function heat_recombination(T_gas)

 real, intent(in) :: T_gas
 real, parameter :: A = 1.28e-19, B = 1.71e-06

 if (T_gas > 1e+04) then
    heat_recombination = A * T_gas**(-0.8) * max(0., 1. - B * T_gas) ! Mathews & Doane 1990 eq. (1)
 else
    heat_recombination = 0.0
 endif

end function heat_recombination

!-----------------------------------------------------------------------
!+
!  RADIATIVE: radiative cooling & heating (Raymond et al. 1976, Mathews & Doane 1990)
!+
!-----------------------------------------------------------------------
real function cool_metal_ions(T_gas)

 real, intent(in) :: T_gas
 real, parameter :: b1 = 1.53e-27, b2 = 1.25e-9, p = 1.2, q = 1.85

 if (T_gas > 1e+04) then
    cool_metal_ions = b1 * T_gas**p / (1. + b2 * T_gas**q) ! Mathews & Doane 1990 eq. (1)
 else
    cool_metal_ions = 0.0
 endif

end function cool_metal_ions

!-----------------------------------------------------------------------
!+
!  RADIATIVE: cooling due to thermal bremsstrahlung (Mathews & Doane 1990)
!+
!-----------------------------------------------------------------------
real function cool_thermal_bremsstrahlung(T_gas)

 real, intent(in) :: T_gas
 real, parameter :: lambda = 2.4e-27

 if (T_gas > 1e+04) then
    cool_thermal_bremsstrahlung = lambda * sqrt(T_gas) ! Mathews & Doane 1990 eq. (1)
 else
    cool_thermal_bremsstrahlung = 0.0
 endif

end function cool_thermal_bremsstrahlung

!-----------------------------------------------------------------------
!+
!  RADIATIVE: Compton heating (Mathews & Doane 1990)
!+
!-----------------------------------------------------------------------
real function heat_Compton(T_gas, rho_gas)

 real, intent(in) :: T_gas, rho_gas
 real, parameter :: C = 3.02e-49, T_Compton = 2e+07

 if (T_gas > 1e+04) then
    heat_Compton = C * (T_Compton - T_gas) / rho_gas ! Mathews & Doane 1990 eq. (1)
 else
    heat_Compton = 0.0
 endif

end function heat_Compton

end module cooling_functions
