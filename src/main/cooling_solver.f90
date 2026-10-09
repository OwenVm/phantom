!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module cooling_solver
!
! Generic module handling analytic cooling functions
!  with known derivatives. These can be solved either
!  explicitly, implicitly or with the Townsend (2009)
!  exact method. Implementation by Lionel Siess.
!
! :References:
!   Townsend (2009), ApJS 181, 391-397
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - CO_abun        : *CO abundance n_CO/n_H2 for CO cooling*
!   - CO_cooling     : *CO rotational line cooling (0=off, 1=Ceulemans+26 fit, 2=Decin+06)*
!   - CR_heating     : *cosmic ray heating, Decin+06 (1=on/0=off)*
!   - G0_UV          : *interstellar far-UV field (Habing units)*
!   - H2O_abun       : *H2O abundance n_H2O/n_H2 for H2O cooling*
!   - H2O_cooling    : *H2O rotational line cooling, Decin+06 (1=on/0=off)*
!   - H2_cooling     : *H2 cooling (1=on/0=off)*
!   - H2vib_cooling  : *H2 vibrational cooling, Decin+06 (1=on/0=off)*
!   - PE_heating     : *photoelectric heating, Decin+06 (1=on/0=off)*
!   - T1_factor      : *factor by which T0 is increased (T1= T1_factor*T0)*
!   - bowen_Cprime   : *radiative cooling rate (g.s/cm³)*
!   - alpha_CO       : *exponent of the CO photodissociation profile*
!   - drift_heating  : *gas-grain drift heating, Decin+06 (1=on/0=off)*
!   - dust_collision : *dust collision (0=off, 1=on, 2=Decin+06 heat exchange)*
!   - dust_to_gas    : *dust-to-gas mass ratio of the MRN grains*
!   - grain_amax     : *maximum grain size (cm)*
!   - grain_amin     : *minimum grain size (cm)*
!   - grain_rhos     : *specific density of the dust grains (g/cm^3)*
!   - excitation_HI  : *cooling via electron excitation of HI (1=on/0=off)*
!   - high_temp      : *radiative cooling for high temperatures (1=on/0=off)*
!   - icool_method   : *integration method (0=implicit, 1=explicit, 2=exact solution)*
!   - lambda_shock   : *Cooling rate parameter for analytic shock solution*
!   - r_half_CO      : *CO photodissociation radius (cm), 0=no dissociation*
!   - relax_bowen    : *Bowen (diffusive) relaxation (1=on/0=off)*
!   - relax_stefan   : *radiative relaxation (1=on/0=off)*
!   - shock_problem  : *piecewise formulation for analytic shock solution (1=on/0=off)*
!
! :Dependencies: cooling_functions, infile_utils, io, physcon, timestep,
!   units
!

 use cooling_functions, only:bowen_Cprime,lambda_shock_cgs,T0_value,T1_factor,CO_abun,&
                              H2O_abun,dust_to_gas,grain_amin,grain_amax,grain_rhos,&
                              G0_UV,r_half_CO,alpha_CO,He_abun
 implicit none
 character(len=*), parameter :: label = 'cooling_library'
 integer, public :: excitation_HI = 0, relax_Bowen = 0, dust_collision = 0, &
                    relax_Stefan = 0, shock_problem = 0, H2_cooling = 0, CO_cooling = 0
 integer, public :: H2O_cooling = 0, H2vib_cooling = 0, drift_heating = 0, &
                    CR_heating = 0, PE_heating = 0
 integer, public :: icool_method  = 0, high_temp = 0
 integer, parameter :: nTg  = 64
 real :: Tref = 1.d7 !higher value of the temperature grid (for exact cooling)
 real :: Tgrid(nTg)

 public :: init_cooling_solver,read_options_cooling_solver,write_options_cooling_solver
 public :: energ_cooling_solver,calc_cooling_rate, calc_Q
 public :: testfunc,print_cooling_rates
 public :: T0_value,lambda_shock_cgs ! expose to cooling module
 logical, public :: Townsend_test = .false. !for analysis_cooling

 private
 real,    parameter :: Tcap = 1.d3 !Townsend cap temperature

contains
!-----------------------------------------------------------------------
!+
!   Initialise cooling functions and check mix of options are sensible
!+
!-----------------------------------------------------------------------
subroutine init_cooling_solver(ierr)
 use io, only:error
 integer, intent(out) :: ierr

 ierr = 0
 !you can't have cool_relaxation_Stefan and cool_relaxation_Bowen at the same time
 if (relax_Bowen == 1 .and. relax_Stefan == 1) then
    call error(label,'you can"t have bowen and stefan cooling at the same time')
    ierr = 1
 endif
 !if no cooling flag activated, disable cooling
 if ( (excitation_HI+relax_Bowen+dust_collision+relax_Stefan+shock_problem+high_temp+H2_cooling+CO_cooling &
       +H2O_cooling+H2vib_cooling+drift_heating+CR_heating+PE_heating) == 0) then
    print *,'ERROR: no cooling prescription activated'
    ierr = 2
 endif

 call set_Tgrid()

end subroutine init_cooling_solver

!-----------------------------------------------------------------------
!+
!   Get right hand side of energy equation for the desired
!   cooling prescription and choice of solver
!+
!-----------------------------------------------------------------------
subroutine energ_cooling_solver(ui,dudt,rho,dt,mu,gamma,Tdust,K2,kappa, i, divv, xyzi, vxyzi, arad_in)
 real, intent(in)  :: ui,rho,dt
 real, intent(in)  :: Tdust,mu,gamma,K2,kappa
 real, intent(out) :: dudt
 integer, optional, intent(in) :: i !index for saving cooling rate in dump (i)
 real,    optional, intent(in) :: divv !velocity divergence (code units)
 real,    optional, intent(in) :: xyzi(3),vxyzi(3) !position and velocity (code units)
 real,    optional, intent(in) :: arad_in !radial radiative acceleration on the gas (code units)

 if (icool_method == 2) then
    call exact_cooling(ui,dudt,rho,dt,mu,gamma,Tdust,K2,kappa, i, divv, xyzi, vxyzi, arad_in)
 elseif (icool_method == 0) then
    call implicit_cooling(ui,dudt,rho,dt,mu,gamma,Tdust,K2,kappa, i, divv, xyzi, vxyzi, arad_in)
 else
    call explicit_cooling(ui,dudt,rho,dt,mu,gamma,Tdust,K2,kappa, i, divv, xyzi, vxyzi, arad_in)
 endif

end subroutine energ_cooling_solver

!-----------------------------------------------------------------------
!+
!   explicit cooling
!+
!-----------------------------------------------------------------------
subroutine explicit_cooling (ui, dudt, rho, dt, mu, gamma, Tdust, K2, kappa, i, divv, xyzi, vxyzi, arad_in)

 use physcon, only:Rg
 use units,   only:unit_ergg

 real, intent(in)  :: ui, rho, dt, Tdust, mu, gamma !code units
 real, intent(in)  :: K2, kappa                     !cgs
 real, intent(out) :: dudt                          !code units
 integer, optional, intent(in) :: i !index for saving cooling rate in dump (i)
 real,    optional, intent(in) :: divv !velocity divergence (code units)
 real,    optional, intent(in) :: xyzi(3),vxyzi(3) !position and velocity (code units)
 real,    optional, intent(in) :: arad_in !radial radiative acceleration on the gas (code units)

 real              :: u,Q,dlnQ_dlnT,T,T_on_u

 T_on_u = (gamma-1.)*mu*unit_ergg/Rg
 T      = T_on_u*ui
 call calc_cooling_rate(Q, dlnQ_dlnT, rho, T, Tdust, mu, gamma, K2, kappa, i, divv, xyzi, vxyzi, arad_in)
 if (ui + Q*dt < 0.) then   ! assume thermal equilibrium
    if (Townsend_test) then
       !special fix for Townsend benchmark
       u = Tcap/T_on_u
    else
       u = Tdust/T_on_u     ! set T=Tdust
    endif
    dudt = (u-ui)/dt
 else
    dudt = Q
 endif

end subroutine explicit_cooling

!-----------------------------------------------------------------------
!+
!   implicit cooling
!+
!-----------------------------------------------------------------------
subroutine implicit_cooling (ui, dudt, rho, dt, mu, gamma, Tdust, K2, kappa, i, divv, xyzi, vxyzi, arad_in)
 use physcon, only:Rg
 use units,   only:unit_ergg

 real, intent(in)  :: ui, rho, dt, mu, gamma
 real, intent(in)  :: Tdust, K2, kappa
 real, intent(out) :: dudt
 integer, optional, intent(in) :: i !index for saving cooling rate in dump (i)
 real,    optional, intent(in) :: divv !velocity divergence (code units)
 real,    optional, intent(in) :: xyzi(3),vxyzi(3) !position and velocity (code units)
 real,    optional, intent(in) :: arad_in !radial radiative acceleration on the gas (code units)

 real, parameter    :: tol = 1.d-6, Tmin = 1.
 integer, parameter :: iter_max = 40
 real               :: u,Q,dlnQ_dlnT,T_on_u,Qi,f0,fi,fmid,T,T0,dx,Tmid
 integer            :: iter

 u       = ui
 T_on_u  = (gamma-1.)*mu*unit_ergg/Rg
 T       = ui*T_on_u
 call calc_cooling_rate(Q,dlnQ_dlnT, rho, T, Tdust, mu, gamma, K2, kappa, i, divv, xyzi, vxyzi, arad_in)
 !cooling negligible, return
 if (abs(Q) < tiny(0.)) then
    dudt = 0.
    return
 endif
 T0   = T
 f0   = -Q*dt*T_on_u
 fi   = f0
 iter = 0
 !define bisection interval for function f(T) = T^(n+1)-T^n-Q*dt*T_on_u
 do while (((f0 > 0. .and. fi > 0.) .or. (f0 < 0. .and. fi < 0.)) .and. iter < iter_max)
    Tmid = max(T+Q*dt*T_on_u,Tmin)
    call calc_cooling_rate(Qi,dlnQ_dlnT, rho, Tmid, Tdust, mu, gamma, K2, kappa, i, divv, xyzi, vxyzi, arad_in)
    fi = Tmid-T0-Qi*dt*T_on_u
    T  = Tmid
    iter = iter+1
 enddo
 !Temperature is between T0 and Tmid
 if (iter > iter_max) stop '[implicit_cooling] cannot bracket cooling function'
 iter = 0
 if (Tmid > T0) then
    T = T0
 else
    if (Townsend_test) then
       !special fix for Townsend benchmark
       T = max(Tcap,Tmid)
    else
       T = Tmid
    endif
 endif
 dx = abs(Tmid-T0)
 do while (dx/T0 > tol .and. iter < iter_max)
    dx = dx*.5
    Tmid = T+dx
    call calc_cooling_rate(Qi,dlnQ_dlnT, rho, Tmid, Tdust, mu, gamma, K2, kappa, i, divv, xyzi, vxyzi, arad_in)
    fmid = Tmid-T0-Qi*dt*T_on_u
    if (Townsend_test) then
       !special fix for Townsend benchmark
       if (fmid <= 0.) Tmid = max(Tcap,Tmid)
    else
       if (fmid <= 0.) T = Tmid
    endif
    iter = iter + 1
    !print *,iter,fmid,T,Tmid
 enddo
 u = Tmid/T_on_u
 dudt =(u-ui)/dt
 if (u < 0. .or. isnan(u)) then
    print *,u
    stop '[implicit_cooling] u<0'
 endif

end subroutine implicit_cooling

!-----------------------------------------------------------------------
!+
!   cooling using Townsend (2009), ApJS 181, 391-397 method with
!   analytical cooling rate prescriptions
!+
!-----------------------------------------------------------------------
subroutine exact_cooling(ui, dudt, rho, dt, mu, gamma, Tdust, K2, kappa, i, divv, xyzi, vxyzi, arad_in)

 use physcon, only:Rg
 use units,   only:unit_ergg

 real, intent(in)  :: ui, rho, dt, Tdust, mu, gamma
 real, intent(in)  :: K2, kappa
 real, intent(out) :: dudt
 integer, optional, intent(in) :: i !index for saving cooling rate in dump (i)
 real,    optional, intent(in) :: divv !velocity divergence (code units)
 real,    optional, intent(in) :: xyzi(3),vxyzi(3) !position and velocity (code units)
 real,    optional, intent(in) :: arad_in !radial radiative acceleration on the gas (code units)

 real, parameter :: tol = 1.d-12
 real            :: Qref,dlnQref_dlnT,Q,dlnQ_dlnT,Y,Yk,Yinv,Temp,dy,T,T_on_u,T_floor,Qi
 integer         :: k

 if (Townsend_test) then
    T_floor = Tcap
 else
    T_floor = 10.
 endif
 T_on_u = (gamma-1.)*mu*unit_ergg/Rg
 T      = T_on_u*ui

 if (T < T_floor) then
    Temp = T_floor
 elseif (T > Tref) then
    call calc_cooling_rate(Q, dlnQ_dlnT, rho, T, Tdust, mu, gamma, K2, kappa, i, divv, xyzi, vxyzi, arad_in)
    Temp = T+T_on_u*Q*dt
 else
    call calc_cooling_rate(Qref,dlnQref_dlnT, rho, Tref, Tdust, mu, gamma, K2, kappa, i, divv, xyzi, vxyzi, arad_in)
    Qi = Qref
    Y         = 0.
    k         = nTg
    Q         = Qref          ! default value if Tgrid < T for all k
    dlnQ_dlnT = dlnQref_dlnT  ! default value if Tgrid < T for all k
    do while (Tgrid(k) > T)
       k = k-1
       call calc_cooling_rate(Q, dlnQ_dlnT, rho, Tgrid(k), Tdust, mu, gamma, K2, kappa, i, divv, xyzi, vxyzi, arad_in)

       if ((Qi /= 0.) .and. (Q /= 0.)) then
          dlnQ_dlnT = log(Qi/Q)/log(Tgrid(k+1)/Tgrid(k))
          dlnQ_dlnT = sign(min(50.,abs(dlnQ_dlnT)),dlnQ_dlnT)
       else
          dlnQ_dlnT = 0.
       endif
       Q = Q-1.d-80 !enforce Q /=0

       Qi = Q
       ! eqs A6 to get Yk
       if (abs(dlnQ_dlnT-1.) < tol) then
          y = y - Qref*Tgrid(k)/(Q*Tref)*log(Tgrid(k)/Tgrid(k+1))
       else
          y = y - Qref*Tgrid(k)/(Q*Tref*(1.-dlnQ_dlnT))*(1.-(Tgrid(k)/Tgrid(k+1))**(dlnQ_dlnT-1.))
       endif
    enddo
    !eqs A5 for Y(T)
    yk = y
    if (abs(dlnQ_dlnT-1.) < tol) then
       y = yk + Qref*Tgrid(k)/(Q*Tref)*log(Tgrid(k)/T)
    else
       y = yk + Qref*Tgrid(k)/((Q*Tref)*(1.-dlnQ_dlnT))*(1.-(Tgrid(k)/T)**(dlnQ_dlnT-1.))
    endif
    !argument of Y^(-1) in eq 26
    dy = -Qref*dt*T_on_u/Tref
    y  = y + dy
    !find new k for eq A7 (not necessarily the same as k for eq A5)
    do while(y>yk .AND. k>1)
       k = k-1
       call calc_cooling_rate(Q, dlnQ_dlnT, rho, Tgrid(k), Tdust, mu, gamma, K2, kappa, i, divv, xyzi, vxyzi, arad_in)

       if ((Qi /= 0.) .and. (Q /= 0.)) then
          dlnQ_dlnT = log(Qi/Q)/log(Tgrid(k+1)/Tgrid(k))
          dlnQ_dlnT = sign(min(50.,abs(dlnQ_dlnT)),dlnQ_dlnT)
       else
          dlnQ_dlnT = 0.
       endif
       Q = Q-1.d-80 !enforce Q /=0

       Qi = Q
       ! eqs A6 to get Yk
       if (abs(dlnQ_dlnT-1.) < tol) then
          yk = yk - Qref*Tgrid(k)/(Q*Tref)*log(Tgrid(k)/Tgrid(k+1))
       else
          yk = yk - Qref*Tgrid(k)/(Q*Tref*(1.-dlnQ_dlnT))*(1.-(Tgrid(k)/Tgrid(k+1))**(dlnQ_dlnT-1.))
       endif
    enddo
    !compute Yinv (eqs A7)
    if (abs(dlnQ_dlnT-1.) < tol) then
       Temp = max(Tgrid(k)*exp(-Q*Tref*(y-yk)/(Qref*Tgrid(k))),T_floor)
    else
       Yinv = 1.-(1.-dlnQ_dlnT)*Q*Tref/(Qref*Tgrid(k))*(y-yk)
       if (Yinv > 0.) then
          Temp = max(Tgrid(k)*(Yinv**(1./(1.-dlnQ_dlnT))),T_floor)
       else
          Temp = T_floor
       endif
    endif
 endif

 dudt = (Temp-T)/T_on_u/dt
 !note that u = Temp/T_on_u

end subroutine exact_cooling

!-----------------------------------------------------------------------
!+
!  calculate cooling rates
!+
!-----------------------------------------------------------------------
subroutine calc_cooling_rate(Q, dlnQ_dlnT, rho, T, Teq, mu, gamma, K2, kappa, i, divv, xyzi, vxyzi, arad_in)
 use units,   only:unit_ergg,unit_density,utime,udist,unit_velocity
 use physcon, only:mass_proton_cgs
 use part,    only:nptmass,xyzmh_ptmass,vxyz_ptmass,iReff,iTeff
 use cooling_functions, only:cooling_neutral_hydrogen,&
     cooling_Bowen_relaxation,cooling_dust_collision,&
     cooling_radiative_relaxation,piecewise_law, &
     cooling_high_temp,testing_cooling_functions,cooling_H2,cooling_CO_rot,&
     decin_densities,heating_drift_decin,cooling_H2O_rot_decin,cooling_CO_rot_decin,&
     cooling_H2_vib_decin,heating_dust_gas_decin,heating_cosmic_rays_decin,&
     heating_photoelectric_decin

 real, intent(in)  :: rho, T, Teq     !rho in code units
 real, intent(in)  :: mu, gamma
 real, intent(in)  :: K2, kappa
 real, intent(out) :: Q, dlnQ_dlnT    !code units
 integer, optional, intent(in) :: i !index for saving cooling rate in dump (i)
 real,    optional, intent(in) :: divv !velocity divergence (code units)
 real,    optional, intent(in) :: xyzi(3),vxyzi(3) !position and velocity (code units)
 real,    optional, intent(in) :: arad_in !radial radiative acceleration on the gas (code units)

 real :: Q_cgs
 real :: Q_H0, Q_relax_Bowen,  Q_col_dust, Q_relax_Stefan, Q_molec, Q_shock, Q_H2, Q_CO
 real :: dlnQ_H0, dlnQ_relax_Bowen, dlnQ_col_dust, dlnQ_relax_Stefan
 real :: dlnQ_molec, dlnQ_shock, dlnQ_H2, dlnQ_CO
 real :: Q_hightemp, dlnQ_hightemp
 real :: Q_decin, dlnQ_decin, Qx, dlnQx
 real :: rho_cgs, ndens, divv_cgs, r, vr, arad, W, Tstar, eps, dx(3), Rstar
 real :: nH, nHI, nH2, nHe

 rho_cgs = rho*unit_density
 ndens   = rho_cgs/mass_proton_cgs

 Q_H0              = 0.
 Q_relax_Bowen     = 0.
 Q_col_dust        = 0.
 Q_relax_Stefan    = 0.
 Q_shock           = 0.
 Q_molec           = 0.
 Q_hightemp        = 0.
 Q_H2              = 0.
 Q_CO              = 0.

 dlnQ_H0           = 0.
 dlnQ_relax_Bowen  = 0.
 dlnQ_col_dust     = 0.
 dlnQ_relax_Stefan = 0.
 dlnQ_shock        = 0.
 dlnQ_molec        = 0.
 dlnQ_hightemp     = 0.
 dlnQ_H2           = 0.
 dlnQ_CO           = 0.

 if (excitation_HI      == 1) call cooling_neutral_hydrogen(T, rho_cgs, Q_H0, dlnQ_H0)

 if (relax_Bowen        == 1) call cooling_Bowen_relaxation(T, Teq, rho_cgs, mu, gamma,Q_relax_Bowen, dlnQ_relax_Bowen)

 if (H2_cooling         == 1) call cooling_H2(T, rho_cgs, Q_H2, dlnQ_H2)

 if (CO_cooling == 1 .and. present(divv)) call cooling_CO_rot(T, rho_cgs, divv/utime, Q_CO, dlnQ_CO)

 ! Decin et al. (2006) terms; Q_decin sums them, dlnQ_decin holds sum Q_k dlnQ_k
 Q_decin    = 0.
 dlnQ_decin = 0.
 if (present(xyzi) .and. present(vxyzi) .and. nptmass > 0) then
    arad  = 0.
    if (present(arad_in)) arad = arad_in*udist/utime**2
    ! distance, radial velocity, dilution factor and T_eff w.r.t. the star (sink 1)
    dx    = xyzi - xyzmh_ptmass(1:3,1)
    r     = sqrt(dot_product(dx,dx))
    vr    = 0.
    if (r > 0.) vr = dot_product(vxyzi - vxyz_ptmass(1:3,1),dx)/r*unit_velocity
    Rstar = xyzmh_ptmass(iReff,1)
    W     = 0.5
    if (r > Rstar) W = 0.5*(1.-sqrt(1.-(Rstar/r)**2))
    Tstar = xyzmh_ptmass(iTeff,1)
    r     = r*udist
    divv_cgs = 0.
    if (present(divv)) divv_cgs = divv/utime
    ! eps = dln v/dln r from div v = dv/dr + 2v/r (spherical outflow)
    eps = 0.
    if (vr > 0. .and. r > 0.) eps = max(r*divv_cgs/vr - 2., 0.)
    call decin_densities(rho_cgs, mu, He_abun, nH, nHI, nH2, nHe)

    if (drift_heating == 1) then
       call heating_drift_decin(T, rho_cgs, mu, arad, Qx, dlnQx)
       Q_decin = Q_decin + Qx; dlnQ_decin = dlnQ_decin + Qx*dlnQx
    endif
    if (dust_collision == 2) then
       call heating_dust_gas_decin(T, Teq, rho_cgs, mu, nH, arad, Qx, dlnQx)
       Q_decin = Q_decin + Qx; dlnQ_decin = dlnQ_decin + Qx*dlnQx
    endif
    if (CR_heating == 1) then
       call heating_cosmic_rays_decin(rho_cgs, nH2, nHI, He_abun, Qx, dlnQx)
       Q_decin = Q_decin + Qx; dlnQ_decin = dlnQ_decin + Qx*dlnQx
    endif
    if (PE_heating == 1) then
       call heating_photoelectric_decin(T, rho_cgs, nH, nH2, r, Qx, dlnQx)
       Q_decin = Q_decin + Qx; dlnQ_decin = dlnQ_decin + Qx*dlnQx
    endif
    if (H2O_cooling == 1) then
       call cooling_H2O_rot_decin(T, rho_cgs, nHI, nH2, nHe, divv_cgs, eps, W, Tstar, Qx, dlnQx)
       Q_decin = Q_decin + Qx; dlnQ_decin = dlnQ_decin + Qx*dlnQx
    endif
    if (CO_cooling == 2) then
       call cooling_CO_rot_decin(T, rho_cgs, nHI, nH2, nHe, r, divv_cgs, eps, W, Tstar, Qx, dlnQx)
       Q_decin = Q_decin + Qx; dlnQ_decin = dlnQ_decin + Qx*dlnQx
    endif
    if (H2vib_cooling == 1) then
       call cooling_H2_vib_decin(T, rho_cgs, nHI, nH2, Qx, dlnQx)
       Q_decin = Q_decin + Qx; dlnQ_decin = dlnQ_decin + Qx*dlnQx
    endif
 endif

 if (dust_collision == 1 .and. K2 > 0.) call cooling_dust_collision(T, Teq, rho_cgs, K2, &
                                                        mu, Q_col_dust, dlnQ_col_dust)
 if (relax_Stefan   == 1) call cooling_radiative_relaxation(T, Teq, kappa, Q_relax_Stefan, &
                                                        dlnQ_relax_Stefan)
 if (shock_problem  == 1) call piecewise_law(T, T0_value, rho_cgs, ndens, Q_H0, dlnQ_H0)
 if (high_temp      == 1) call cooling_high_temp(T, rho_cgs, Q_hightemp, dlnQ_hightemp)

 if (excitation_HI  == 99) call testing_cooling_functions(int(K2), T, Q_H0, dlnQ_H0)
 !if (do_molecular_cooling) call calc_cool_molecular(T, r, rho_cgs, Q_molec, dlnQ_molec)

 Q_cgs = Q_H0 + Q_relax_Bowen + Q_col_dust + Q_relax_Stefan + Q_molec + Q_shock + Q_hightemp + Q_H2 + Q_CO &
       + Q_decin

 if (Q_cgs == 0.) then
    dlnQ_dlnT = 0.
 else
    dlnQ_dlnT = (Q_H0*dlnQ_H0 + Q_relax_Bowen*dlnQ_relax_Bowen + Q_col_dust*dlnQ_col_dust &
               + Q_relax_Stefan*dlnQ_relax_Stefan + Q_molec*dlnQ_molec + Q_shock*dlnQ_shock&
   + Q_hightemp*dlnQ_hightemp &
               + Q_H2*dlnQ_H2 + Q_CO*dlnQ_CO + dlnQ_decin) / Q_cgs
 endif

 !limit exponent to prevent overflow
 dlnQ_dlnT = sign(min(50.,abs(dlnQ_dlnT)),dlnQ_dlnT)
 Q         = Q_cgs/(unit_ergg/utime)

 !call testfunc()
 !call exit

end subroutine calc_cooling_rate

!-----------------------------------------------------------------------
!+
!  UTILITY: Total cooling function
!+
!-----------------------------------------------------------------------
real function calc_Q(T_gas, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas, &
                     T_dust, v_drift, d2g, a, rho_grain, kappa_dust, JL)
 use cooling_functions
 real, intent(in) :: T_gas, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas
 real, intent(in) :: T_dust, v_drift, d2g, a, rho_grain, kappa_dust, JL

 calc_Q =  cool_dust_discrete_contact(T_gas, rho_gas, mu, T_dust, d2g, a, rho_grain, kappa_dust) &
!     + cool_dust_full_contact(T_gas, rho_gas, mu, T_dust, kappa_dust) &
!     + cool_dust_radiation(T_gas, kappa_gas, T_dust, kappa_dust) &
    + cool_coulomb(T_gas, rho_gas, mu, nH, nHe) &
    + cool_HI(T_gas, rho_gas, mu, nH, nHe) &
    + cool_H_ionisation(T_gas, rho_gas, mu, nH, nHe) &
    + cool_He_ionisation(T_gas, rho_gas, mu, nH, nHe) &
    + cool_H2_rovib(T_gas, nH, nH2) &
    + cool_H2_dissociation(T_gas, rho_gas, mu, nH, nH2) &
    + cool_CO_rovib(T_gas, rho_gas, mu, nH, nH2, nCO) &
    + cool_H2O_rovib(T_gas, rho_gas, mu, nH, nH2, nH2O) &
    + cool_OH_rot(T_gas, rho_gas, mu, nOH) &
    + cool_metal_ions(T_gas) &
    + cool_thermal_bremsstrahlung(T_gas) &
    - heat_recombination(T_gas) &
    - heat_Compton(T_gas, rho_gas) &
    - heat_dust_friction(rho_gas, v_drift, d2g, a, rho_grain, kappa_dust) &
    - heat_dust_photovoltaic_soft(T_gas, rho_gas, mu, nH, nHe, kappa_dust) &
!     - heat_dust_photovoltaic_hard(T_gas, nH, d2g, kappa_dust, JL) &
    - heat_CosmicRays(nH, nH2)
!     - heat_H2_recombination(T_gas, rho_gas, mu, nH, nH2, T_dust)

!  call print_cooling_rates(T_gas, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas, &
!                     T_dust, v_drift, d2g, a, rho_grain, kappa_dust, JL)

end function calc_Q

!-----------------------------------------------------------------------
!+
!  UTILITY: numerical estimate of dlnQ/dlnT
!+
!-----------------------------------------------------------------------
real function calc_dlnQdlnT(T_gas, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas, &
                            T_dust, v_drift, d2g, a, rho_grain, kappa_dust, JL)

 use timestep, only:bignumber

 real, intent(in) :: T_gas, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas
 real, intent(in) :: T_dust, v_drift, d2g, a, rho_grain, kappa_dust
 real, intent(in) :: JL

 real, parameter    :: tolQ    = 1.d-4
 real               :: Qtot, dlnQ_dlnT, dT, Q1, Q2, dQdT
 integer, parameter :: itermax = 20
 integer            :: iter

 Qtot      = calc_Q(T_gas, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas, &
                    T_dust, v_drift, d2g, a, rho_grain, kappa_dust, JL)
 dlnQ_dlnT = 0.

! centered finite order approximation for the numerical derivative
 dT = T_gas/100.
 Q1 = calc_Q(T_gas+dT, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas, &
              T_dust, v_drift, d2g, a, rho_grain, kappa_dust, JL)
 Q2 = calc_Q(T_gas-dT, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas, &
              T_dust, v_drift, d2g, a, rho_grain, kappa_dust, JL)

 iter    = 0
 do while ( abs(Q1/Qtot-1.) > tolQ .and. abs(Q2/Qtot-1.) > tolQ .and. iter < itermax )
    dT   = dT/2.
    Q1   = calc_Q(T_gas+dT, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas, &
                  T_dust, v_drift, d2g, a, rho_grain, kappa_dust, JL)
    Q2   = calc_Q(T_gas-dT, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas, &
                  T_dust, v_drift, d2g, a, rho_grain, kappa_dust, JL)
    iter = iter + 1
 enddo

 dQdT      = (Q1-Q2)/(2.*dT)
 dlnQ_dlnT = (T_gas/Qtot)*dQdT

 if (dlnQ_dlnT > bignumber) then
    dlnQ_dlnT = 0.
 endif

 calc_dlnQdlnT = dlnQ_dlnT

end function calc_dlnQdlnT

!-----------------------------------------------------------------------
!+
!  Set Temperature grid for exact cooling: between 10K and Tref
!+
!-----------------------------------------------------------------------
subroutine set_Tgrid
 integer :: i
 real    :: dlnT,T1

 logical, parameter :: logscale = .true.
 real :: Tmin = 10.

 if (shock_problem  == 1) then
    T1 =   T1_factor * T0_value
    Tref = T1 - (T1 - T0_value)/10000. !slightly below T1 so Qref /= 0
    !Tref = (T1_factor - (T1_factor - 1.))*T0_value/10000.
 endif

 if (logscale) then
    dlnT = log(Tref/Tmin)/(nTg-1)
    do i = 1,nTg
       Tgrid(i) = Tmin*exp((i-1)*dlnT)
       !print *,i,Tgrid(i)
    enddo
 else
    dlnT = (Tref-Tmin)/(nTg-1)
    do i = 1,nTg
       Tgrid(i) = Tmin+(i-1)*dlnT
       !print *,i,Tgrid(i)
    enddo
 endif
end subroutine set_Tgrid

!-----------------------------------------------------------------------
!+
!  writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_cooling_solver(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 !call write_options_molecularcooling(iunit)
 call write_inopt(icool_method,'icool_method','integration method (0=implicit, 1=explicit, 2=exact solution)',iunit)
 call write_inopt(excitation_HI,'excitation_HI','cooling via electron excitation of HI (1=on/0=off)',iunit)
 call write_inopt(relax_bowen,'relax_bowen','Bowen (diffusive) relaxation (1=on/0=off)',iunit)
 call write_inopt(relax_stefan,'relax_stefan','radiative relaxation (1=on/0=off)',iunit)
 call write_inopt(dust_collision,'dust_collision','dust collision (0=off, 1=on, 2=Decin+06 heat exchange)',iunit)
 call write_inopt(shock_problem,'shock_problem','piecewise formulation for analytic shock solution (1=on/0=off)',iunit)
 call write_inopt(high_temp,'high_temp','radiative cooling for high temperatures (1=on/0=off)',iunit)
 call write_inopt(H2_cooling,'H2_cooling','H2 cooling (1=on/0=off)',iunit)
 call write_inopt(CO_cooling,'CO_cooling','CO rotational line cooling (0=off, 1=Ceulemans+26 fit, 2=Decin+06)',iunit)
 if (CO_cooling > 0 .or. PE_heating == 1) &
    call write_inopt(CO_abun,'CO_abun','CO abundance n_CO/n_H2 for CO cooling',iunit)
 call write_inopt(H2O_cooling,'H2O_cooling','H2O rotational line cooling, Decin+06 (1=on/0=off)',iunit)
 if (H2O_cooling == 1) call write_inopt(H2O_abun,'H2O_abun','H2O abundance n_H2O/n_H2 for H2O cooling',iunit)
 call write_inopt(H2vib_cooling,'H2vib_cooling','H2 vibrational cooling, Decin+06 (1=on/0=off)',iunit)
 call write_inopt(drift_heating,'drift_heating','gas-grain drift heating, Decin+06 (1=on/0=off)',iunit)
 call write_inopt(CR_heating,'CR_heating','cosmic ray heating, Decin+06 (1=on/0=off)',iunit)
 call write_inopt(PE_heating,'PE_heating','photoelectric heating, Decin+06 (1=on/0=off)',iunit)
 if (drift_heating == 1 .or. dust_collision == 2 .or. PE_heating == 1) then
    call write_inopt(dust_to_gas,'dust_to_gas','dust-to-gas mass ratio of the MRN grains',iunit)
    call write_inopt(grain_amin,'grain_amin','minimum grain size (cm)',iunit)
    call write_inopt(grain_amax,'grain_amax','maximum grain size (cm)',iunit)
    call write_inopt(grain_rhos,'grain_rhos','specific density of the dust grains (g/cm^3)',iunit)
 endif
 if (PE_heating == 1) call write_inopt(G0_UV,'G0_UV','interstellar far-UV field (Habing units)',iunit)
 if (PE_heating == 1 .or. CO_cooling == 2) then
    call write_inopt(r_half_CO,'r_half_CO','CO photodissociation radius (cm), 0=no dissociation',iunit)
    call write_inopt(alpha_CO,'alpha_CO','exponent of the CO photodissociation profile',iunit)
 endif
 if (shock_problem == 1) then
    call write_inopt(lambda_shock_cgs,'lambda_shock','Cooling rate parameter for analytic shock solution',iunit)
    call write_inopt(T1_factor,'T1_factor','factor by which T0 is increased (T1= T1_factor*T0)',iunit)
    call write_inopt(T0_value,'T0','temperature to cool towards (do not modify! set by setup)',iunit)
 endif
 call write_inopt(bowen_Cprime,'bowen_Cprime','radiative cooling rate (g.s/cm³)',iunit)

end subroutine write_options_cooling_solver

!-----------------------------------------------------------------------
!+
!  reads input options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_cooling_solver(db,nerr)
 use infile_utils, only:inopts,read_inopt
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr

 call read_inopt(icool_method,'icool_method',db,errcount=nerr,min=0,max=2)
 call read_inopt(excitation_HI,'excitation_HI',db,errcount=nerr,min=0,max=1)
 call read_inopt(relax_bowen,'relax_bowen',db,errcount=nerr,min=0,max=1)
 call read_inopt(relax_stefan,'relax_stefan',db,errcount=nerr,min=0,max=1)
 call read_inopt(dust_collision,'dust_collision',db,errcount=nerr,min=0,max=2)
 call read_inopt(shock_problem,'shock_problem',db,errcount=nerr,min=0,max=1)
 call read_inopt(high_temp,'high_temp',db,errcount=nerr,min=0,max=1,default=high_temp)
 call read_inopt(H2_cooling,'H2_cooling',db,errcount=nerr,min=0,max=1,default=H2_cooling)
 call read_inopt(CO_cooling,'CO_cooling',db,errcount=nerr,min=0,max=2,default=CO_cooling)
 call read_inopt(H2O_cooling,'H2O_cooling',db,errcount=nerr,min=0,max=1,default=H2O_cooling)
 call read_inopt(H2vib_cooling,'H2vib_cooling',db,errcount=nerr,min=0,max=1,default=H2vib_cooling)
 call read_inopt(drift_heating,'drift_heating',db,errcount=nerr,min=0,max=1,default=drift_heating)
 call read_inopt(CR_heating,'CR_heating',db,errcount=nerr,min=0,max=1,default=CR_heating)
 call read_inopt(PE_heating,'PE_heating',db,errcount=nerr,min=0,max=1,default=PE_heating)
 if (CO_cooling > 0 .or. PE_heating == 1) &
    call read_inopt(CO_abun,'CO_abun',db,errcount=nerr,min=0.,default=CO_abun)
 if (H2O_cooling == 1) call read_inopt(H2O_abun,'H2O_abun',db,errcount=nerr,min=0.,default=H2O_abun)
 if (drift_heating == 1 .or. dust_collision == 2 .or. PE_heating == 1) then
    call read_inopt(dust_to_gas,'dust_to_gas',db,errcount=nerr,min=0.,default=dust_to_gas)
    call read_inopt(grain_amin,'grain_amin',db,errcount=nerr,min=0.,default=grain_amin)
    call read_inopt(grain_amax,'grain_amax',db,errcount=nerr,min=0.,default=grain_amax)
    call read_inopt(grain_rhos,'grain_rhos',db,errcount=nerr,min=0.,default=grain_rhos)
 endif
 if (PE_heating == 1) call read_inopt(G0_UV,'G0_UV',db,errcount=nerr,min=0.,default=G0_UV)
 if (PE_heating == 1 .or. CO_cooling == 2) then
    call read_inopt(r_half_CO,'r_half_CO',db,errcount=nerr,min=0.,default=r_half_CO)
    call read_inopt(alpha_CO,'alpha_CO',db,errcount=nerr,min=0.,default=alpha_CO)
 endif
 if (shock_problem == 1) then
    call read_inopt(lambda_shock_cgs,'lambda_shock',db,errcount=nerr,min=0.)
    call read_inopt(T1_factor,'T1_factor',db,errcount=nerr,min=0.)
    call read_inopt(T0_value,'T0',db,errcount=nerr,min=0.)
 endif
 call read_inopt(bowen_Cprime,'bowen_Cprime',db,errcount=nerr,min=0.,default=bowen_Cprime)

end subroutine read_options_cooling_solver

!=======================================================================
!
!  Test routine for cooling functions
!
!=======================================================================

subroutine testfunc()

 use physcon, only:mass_proton_cgs

 real :: T_gas, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas
 real :: T_dust, v_drift, d2g, a, rho_grain, kappa_dust, JL
 real :: n_gas

 ! evaluate parameters using plausible values to test the cooling functions

 T_gas      = 4000.
 rho_gas    = 2.d-16
 mu         = 1.78

 n_gas      = rho_gas/(mu*mass_proton_cgs)

 nH         = 0.5 *n_gas
 nH2        = 0.5 *n_gas
 nHe        = 0.1 *n_gas
 nCO        = 1.d-4*n_gas
 nH2O       = 5.d-5*n_gas
 nOH        = 1.d-7*n_gas
 kappa_gas  = 2.d-4

 T_dust     = 1500.
 v_drift    = 1.d6
 d2g        = 1./200.
 a          = 1.d-5
 rho_grain  = 2.
 kappa_dust = 1.d-4

 JL = 2.5d-12     ! Value taken from Barstow et al. 1997
 call print_cooling_rates(T_gas, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas, &
                     T_dust, v_drift, d2g, a, rho_grain, kappa_dust, JL)

end subroutine testfunc

subroutine print_cooling_rates(T_gas, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas, &
                     T_dust, v_drift, d2g, a, rho_grain, kappa_dust, JL)
 use cooling_functions
 real, intent(in) :: T_gas, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas
 real, intent(in) :: T_dust, v_drift, d2g, a, rho_grain, kappa_dust
 real, intent(in) :: JL
 real :: Q1, Q2, Q3, Q4, Q5, Q6, Q7, Q8, Q9, Q10, Q11, Q12, Q13, Q14, Q15, Q16, Q17
 real :: Qtot, dlnQ_dlnT, nH_tot

 !nH_tot = nH+2.*nH2
 nH_tot = 1.
 print*, ' '
 print*, ' '
 print*, '-----------------------------------------------------------------------------'
 print*, ' '
 print*, ' '
 Q1  = cool_dust_discrete_contact(T_gas, rho_gas, mu, T_dust, d2g, a, rho_grain, kappa_dust)
 print*, 'Q1   = ', Q1/nH_tot
 Q2  = cool_dust_full_contact(T_gas, rho_gas, mu, T_dust, kappa_dust)
 print*, 'Q2   = ', Q2/nH_tot
 Q3  = cool_dust_radiation(T_gas, kappa_gas, T_dust, kappa_dust)
 print*, 'Q3   = ', Q3/nH_tot
 Q4  = cool_coulomb(T_gas, rho_gas, mu, nH, nHe)
 print*, 'Q4   = ', Q4/nH_tot
 Q5  = cool_HI(T_gas, rho_gas, mu, nH, nHe)
 print*, 'Q5   = ', Q5/nH_tot
 Q6  = cool_H_ionisation(T_gas, rho_gas, mu, nH, nHe)
 print*, 'Q6   = ', Q6/nH_tot
 Q7  = cool_He_ionisation(T_gas, rho_gas, mu, nH, nHe)
 print*, 'Q7   = ', Q7/nH_tot
 Q8  = cool_H2_rovib(T_gas, nH, nH2)
 print*, 'Q8   = ', Q8/nH_tot
 Q9  = cool_H2_dissociation(T_gas, rho_gas, mu, nH, nH2)
 print*, 'Q9   = ', Q9/nH_tot
 Q10 = cool_CO_rovib(T_gas, rho_gas, mu, nH, nH2, nCO)
 print*, 'Q10  = ', Q10/nH_tot
 Q11 = cool_H2O_rovib(T_gas, rho_gas, mu, nH, nH2, nH2O)
 print*, 'Q11  = ', Q11/nH_tot
 Q12 = cool_OH_rot(T_gas, rho_gas, mu, nOH)
 print*, 'Q12  = ', Q12/nH_tot
 Q13 = heat_dust_friction(rho_gas, v_drift, d2g, a, rho_grain, kappa_dust)
 print*, 'Q13  = ', Q13/nH_tot
 Q14 = heat_dust_photovoltaic_soft(T_gas, rho_gas, mu, nH, nHe, kappa_dust)
 print*, 'Q14  = ', Q14/nH_tot
 Q15 = heat_dust_photovoltaic_hard(T_gas, nH, d2g, kappa_dust, JL)
 print*, 'Q15  = ', Q15/nH_tot
 Q16 = heat_CosmicRays(nH, nH2)
 print*, 'Q16  = ', Q16/nH_tot
 Q17 = heat_H2_recombination(T_gas, rho_gas, mu, nH, nH2, T_dust)
 print*, 'Q17  = ', Q17/nH_tot

 Qtot = calc_Q(T_gas, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas, &
                     T_dust, v_drift, d2g, a, rho_grain, kappa_dust, JL)
 print*, 'Qtot = ', Qtot/nH_tot

 dlnQ_dlnT = calc_dlnQdlnT(T_gas, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas, &
                            T_dust, v_drift, d2g, a, rho_grain, kappa_dust, JL)
 print*, 'dlnQdlnT = ', dlnQ_dlnT

 print*, ' '
 print*, ' '
 print*, '------------------- exit in calc_cooling_rate --------------------------------'
 print*, ' '
 print*, ' '

end subroutine print_cooling_rates

end module cooling_solver
