!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module set_star
!
! driver to integrate the hydrostatic equilibrium equations
!
! :References: None
!
! :Owner: Owen Vermeulen
!
! :Runtime parameters: None
!
! :Dependencies: eos, io, physcon, table_utils, units
!

 use dim, only:isothermal
 implicit none
 public :: setup_star
 public :: stellar_state,save_stellarprofile,interp_stellar_profile

 private
 ! Shared variables
 real, parameter :: rho_power = 2.0  ! Density profile exponent
 character(len=*), parameter :: label = 'set_star'

 ! input parameters
 real :: Mstar_cgs, Rstar_cgs, star_gamma, star_mu
 real, dimension(:,:), allocatable, public :: stellar_1D

 ! stellar properties
 type stellar_state
    real :: r, r0, Rstar, rho, P, u, T, mr
    integer :: nsteps
    logical :: error
 end type stellar_state

contains

subroutine setup_star(Mstar_in, Rstar_in, mu_in, gamma_in, n)
 use physcon, only:solarm, au
 use eos,     only:gamma, gmw
 use io,      only:iverbose

 real, intent(in)    :: Mstar_in, Rstar_in, mu_in, gamma_in
 integer, intent(in) :: n

 Mstar_cgs  = Mstar_in * solarm
 Rstar_cgs  = Rstar_in * au
 star_gamma = gamma_in
 star_mu    = mu_in

 if (iverbose >= 1) then
    print *, '========================================='
    print *, 'Setting up stellar hydrostatic profile'
    print *, '========================================='
    print *, 'M* = ', Mstar_in, ' Msun'
    print *, 'R* = ', Rstar_in, ' AU'
    print *, 'mu = ', mu_in
    print *, 'gamma = ', gamma_in
 endif

 ! Calculate stellar structure
 call calc_stellar_profile(n)

end subroutine setup_star

!-----------------------------------------------------------------------
!
!  Initialize variables for stellar profile integration
!
!-----------------------------------------------------------------------
subroutine init_stellar(n, state)
! all quantities in cgs
 use physcon, only:pi, Rg, kboltz, mass_proton_cgs
 integer, intent(in) :: n
 type(stellar_state), intent(out) :: state

 ! Initialize at stellar surface
 state%r0     = Rstar_cgs
 state%r      = Rstar_cgs
 state%Rstar  = Rstar_cgs
 state%P      = 0.  ! Zero pressure at surface
 state%rho    = Mstar_cgs / (4.*pi * Rstar_cgs**rho_power * Rstar_cgs)
 state%u      = state%P / (state%rho * (star_gamma - 1.))
 state%T      = star_mu * mass_proton_cgs / kboltz * (star_gamma - 1.) * state%u
 state%mr     = Mstar_cgs
 state%nsteps = 1
 state%error  = .false.

end subroutine init_stellar

!-----------------------------------------------------------------------
!
!  Integrate hydrostatic equilibrium over one radial step
!
!-----------------------------------------------------------------------
subroutine stellar_step(state, r_new)
! all quantities in cgs
 use physcon, only:Gg, pi, Rg, kboltz, mass_proton_cgs

 type(stellar_state), intent(inout) :: state
 real, intent(in) :: r_new
 real :: dr, r_mid, rho_mid, mr_mid, dP, C_rho

 dr = r_new - state%r
 r_mid = 0.5 * (state%r + r_new)

 ! Get density and enclosed mass at midpoint
 C_rho = Mstar_cgs / (4.*pi * Rstar_cgs)
 rho_mid = C_rho / r_mid**rho_power
 mr_mid = Mstar_cgs * r_mid / Rstar_cgs

 ! Integrate hydrostatic equilibrium: dP/dr = -rho * G * M(r) / r^2
 dP = -(Gg * mr_mid * rho_mid / r_mid**2) * dr

 ! Update state
 state%r   = r_new
 state%P   = state%P + dP
 state%rho = C_rho / state%r**rho_power
 state%mr  = Mstar_cgs * state%r / Rstar_cgs

 ! Calculate thermodynamic quantities
 state%u = state%P / (state%rho * (star_gamma - 1.))
 state%T = star_mu * mass_proton_cgs / kboltz * (star_gamma - 1.) * state%u

 state%nsteps = state%nsteps + 1

end subroutine stellar_step

!-----------------------------------------------------------------------
!
!  Integrate the hydrostatic equilibrium equation
!
!-----------------------------------------------------------------------
subroutine calc_stellar_profile(n)
! all quantities in cgs
 integer, intent(in) :: n
 type(stellar_state) :: state
 integer :: i
 real :: r_new

 ! Initialize stellar structure
 call init_stellar(n, state)

 ! Allocate storage for profile
 if (allocated(stellar_1D)) deallocate(stellar_1D)
 allocate(stellar_1D(5, n))

 ! Store surface values
 stellar_1D(1, n) = state%r
 stellar_1D(2, n) = state%rho
 stellar_1D(3, n) = state%P
 stellar_1D(4, n) = state%u
 stellar_1D(5, n) = state%T

 ! Integrate inward from surface to center
 do i = n-1, 1, -1
    r_new = Rstar_cgs * real(i) / real(n)
    call stellar_step(state, r_new)

    ! Store in profile
    stellar_1D(1, i) = state%r
    stellar_1D(2, i) = state%rho
    stellar_1D(3, i) = state%P
    stellar_1D(4, i) = state%u
    stellar_1D(5, i) = state%T
 enddo

 ! Save profile to file
!  call save_stellarprofile(n, 'stellar_profile1D.dat')

end subroutine calc_stellar_profile

!-----------------------------------------------------------------------
!
!  Interpolate stellar profile at given radius
!
!-----------------------------------------------------------------------
subroutine interp_stellar_profile(r, rho, P, u, T)
 !in/out variables in code units
 use units,       only:udist, unit_density, unit_ergg
 use table_utils, only:find_nearest_index, interp_1d
 use io,          only:fatal

 real, intent(in)  :: r  ! in code units
 real, intent(out) :: rho, P, u, T
 real :: r_cgs
 integer :: indx, n

 ! Check if profile exists
 if (.not. allocated(stellar_1D)) then
    call fatal(label, 'stellar_1D not allocated. Call setup_star first.')
 endif

 n = size(stellar_1D, 2)
 r_cgs = r * udist

 ! Handle boundary cases
 if (r_cgs <= stellar_1D(1,1)) then
    ! Below minimum radius - use center values
    rho = stellar_1D(2, 1) / unit_density
    P   = stellar_1D(3, 1)
    u   = stellar_1D(4, 1) / unit_ergg
    T   = stellar_1D(5, 1)
    return
 elseif (r_cgs >= stellar_1D(1, n)) then
    ! Above maximum radius - use surface values
    rho = stellar_1D(2, n) / unit_density
    P   = stellar_1D(3, n)
    u   = stellar_1D(4, n) / unit_ergg
    T   = stellar_1D(5, n)
    return
 endif

 ! Find nearest index and interpolate (in cgs)
 call find_nearest_index(stellar_1D(1,:), r_cgs, indx)

 rho = interp_1d(r_cgs, stellar_1D(1,indx), stellar_1D(1,indx+1), &
                        stellar_1D(2,indx), stellar_1D(2,indx+1))
 P   = interp_1d(r_cgs, stellar_1D(1,indx), stellar_1D(1,indx+1), &
                        stellar_1D(3,indx), stellar_1D(3,indx+1))
 u   = interp_1d(r_cgs, stellar_1D(1,indx), stellar_1D(1,indx+1), &
                        stellar_1D(4,indx), stellar_1D(4,indx+1))
 T   = interp_1d(r_cgs, stellar_1D(1,indx), stellar_1D(1,indx+1), &
                        stellar_1D(5,indx), stellar_1D(5,indx+1))

 ! Convert to code units
 rho = rho / unit_density
 u   = u / unit_ergg

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

 ! Write profile data
 do i = 1, n
    call filewrite_stellar_state(iunit, nwrite, i)
 enddo

 close(iunit)

 if (iverbose >= 1) then
    print *, 'Stellar profile saved successfully'
    print *, 'Points written: ', n
    print *, 'Surface: T = ', stellar_1D(5, n), ' K'
    print *, 'Center:  T = ', stellar_1D(5, 1), ' K'
    print *, 'Radius range: ', stellar_1D(1,1)/au, ' to ', stellar_1D(1,n)/au, ' AU'
 endif

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
 use physcon, only:pi
 integer, intent(in) :: i
 real, intent(out) :: array(:)
 real :: C_rho, mr

 ! Calculate auxiliary quantities
 C_rho = Mstar_cgs / (4.*pi * Rstar_cgs)
 mr = Mstar_cgs * stellar_1D(1,i) / Rstar_cgs

 array(1) = stellar_1D(1, i)  ! r
 array(2) = stellar_1D(2, i)  ! rho
 array(3) = stellar_1D(3, i)  ! P
 array(4) = stellar_1D(4, i)  ! u
 array(5) = stellar_1D(5, i)  ! T
 array(6) = mr                 ! enclosed mass
 array(7) = C_rho / stellar_1D(1,i)**rho_power  ! analytic density

end subroutine state_to_array

subroutine filewrite_stellar_state(iunit, nwrite, i)
 integer, intent(in) :: iunit, nwrite, i
 real :: array(nwrite)
 character(len=20) :: fmt

 call state_to_array(i, array)
 write(fmt,*) nwrite
 write(iunit,'('// adjustl(fmt) //'(1x,es14.6E3:))') array(1:nwrite)

end subroutine filewrite_stellar_state

end module set_star