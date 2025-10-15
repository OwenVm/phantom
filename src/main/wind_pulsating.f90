!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module wind_pulsating
!
! :Owner: Owen Vermeulen
!
! :Runtime parameters: None
!
 use dim, only : isothermal
 implicit none
 public :: setup_wind_pulsating

contains

subroutine setup_wind_pulsating
! Setup the wind pulsating module
end subroutine setup_wind_pulsating

!-----------------------------------------------------------------------
!
!  Interpolate 1D stellar profile
!
!-----------------------------------------------------------------------
subroutine interp_wind_profile(r, v)
 use units,      only:udist, utime, unit_velocity, unit_density, unit_ergg
 use part,       only:idgamma 
 use eos,        only:gamma
 use table_utils only:interp_1d

 real, intent(in)  

end subroutine interp_wind_profile

end module wind_pulsating