!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Setup routine for pulsating AGB star with atmospheric wind
!
! :References: None
!
! :Owner: Owen Vermeulen
!
! :Runtime parameters:
!   - Mstar_Msun           : *stellar mass (solar masses)*
!   - Rstar_au             : *stellar radius (AU)*
!   - Tstar                : *effective temperature (K)*
!   - Psurf                : *pressure at the surface (dyn/cm^2)*
!   - atmos_mass_fraction  : *atmospheric mass as fraction of total*
!   - n_shells_total       : *total number of atmospheric shells*
!   - iboundary_spheres    : *number of inner boundary shells*
!   - r_min_on_rstar       : *inner atmosphere radius / R_star*
!   - pulsation_amplitude  : *amplitude of the pulsations in km/s*
!   - iwind_resolution     : *geodesic sphere resolution*
!
! :Dependencies: boundary, dim, infile_utils, io, kernel, options, part,
!   physcon, prompting, ptmass, setup_params, table_utils, timestep, units
!
 use part, only: nptmass
 implicit none
 
 public :: setpart
 
 ! Default parameters
 real    :: Mstar_Msun = 1.0
 real    :: Rstar_au = 1.0
 real    :: Tstar = 3000.0
 real    :: Psurf  = 300.0
 real    :: atmos_mass_fraction = 0.01
 integer :: n_shells_total = 50
 integer :: iboundary_spheres = 5
 real    :: r_min_on_rstar = 0.8
 real    :: pulsation_amplitude = 10
 integer :: iwind_resolution = 5
 
 private

contains

!-----------------------------------------------------------------------
!+
!  Calculate pulsation period based on stellar mass and radius
!  Using empirical relation from Ostlie & Cox (1986)
!+
!+-----------------------------------------------------------------------

subroutine calculate_period(M, R)
 use physcon, only:solarm,au,days
 real, intent(in)  :: M, R
 real              :: logP, logM, logR
 real, intent(out) :: pulsation_period_days

 logM = log10(M)
 logR = log10(R)
 logP = -1.92 - 0.73*logM + 1.86*logR
 pulsation_period_days = 10.0**logP

end subroutine calculate_period

!----------------------------------------------------------------
!+
!  Setup routine for pulsating AGB star with wind
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use physcon,      only:solarm,au,pi
 use units,        only:set_units,umass,udist,unit_luminosity
 use part,         only:nptmass,xyzmh_ptmass,vxyz_ptmass,ihacc,ihsoft
 use part,         only:igas,iTeff,iReff,ilum
 use io,           only:master,fatal
 use timestep,     only:dtmax,tmax
 use options,      only:iexternalforce
 use eos,          only:gmw
 use kernel,       only:hfact_default
 use ptmass,       only:h_acc,r_crit,icreate_sinks
 
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 
 real :: Lstar_Lsun, Lstar_cgs, M_atmos
 integer :: ierr
 
 ! Set units (mass in Msun, distance in AU, G=1)
 call set_units(mass=solarm,dist=au,G=1.d0)
 
 ! Equation of state
 gamma = 5./3.
 gmw = 2.381  ! Mean molecular weight
 polyk = 0.
 hfact = hfact_default
 time = 0.
 
 ! No particles initially - atmosphere will be created by inject module
 npart = 0
 npartoftype(:) = 0
 massoftype(:) = 0.
 
 ! Setup sink particle for central star
 nptmass = 1
 
 ! Position and velocity at origin
 xyzmh_ptmass(:,1) = 0.
 vxyz_ptmass(:,1) = 0.
 
 ! The inject module will split this into sink + atmosphere
 xyzmh_ptmass(4,1) = Mstar_Msun
 
 ! Accretion radius (small, inside inner atmosphere boundary)
 h_acc = 0.5 * r_min_on_rstar * Rstar_au
 xyzmh_ptmass(ihacc,1) = h_acc
 xyzmh_ptmass(ihsoft,1) = 0.0  ! No softening initially
 
 ! Stellar properties
 xyzmh_ptmass(iTeff,1) = Tstar
 xyzmh_ptmass(iReff,1) = Rstar_au
 
 ! Luminosity (Stefan-Boltzmann law)
 ! L = 4π R² σ T⁴
 Lstar_cgs = 4.0*pi * (Rstar_au*au)**2 * 5.67e-5 * Tstar**4
 Lstar_Lsun = Lstar_cgs / 3.839e33
 xyzmh_ptmass(ilum,1) = Lstar_Lsun * 3.839e33 / unit_luminosity
 
 ! Timestep and simulation parameters
 tmax = 1000.  ! Simulation time (code units ~ years)
 dtmax = 10.   ! Time between dumps
 
 ! No external forces initially
 iexternalforce = 0
 
 ! Sink creation off (we already have our sink)
 icreate_sinks = 0
 r_crit = 0.
 
 if (id==master) then
    print "(a)",'================================================================'
    print "(a,f8.3,a)",'  Stellar mass (total)      = ',Mstar_Msun,' M_sun'
    print "(a,f8.3,a)",'  Stellar radius             = ',Rstar_au,' AU'
    print "(a,f8.1,a)",'  Effective temperature      = ',Tstar,' K'
    print "(a,f8.3,a)",'  Luminosity                 = ',Lstar_Lsun,' L_sun'
    print "(a,f8.4,a)",'  Accretion radius           = ',h_acc,' AU'
    print "(a)",' '
    print "(a,f6.2,a)",'  Atmosphere mass fraction   = ',atmos_mass_fraction*100.,' %'
    print "(a,i5)",    '  Total shells               = ',n_shells_total
    print "(a,i5)",    '  Boundary shells            = ',iboundary_spheres
    print "(a,i5)",    '  Free shells                = ',n_shells_total-iboundary_spheres
    print "(a,f6.3)",  '  Inner radius / R_star      = ',r_min_on_rstar
    print "(a,i5)",    '  Wind resolution            = ',iwind_resolution
    print "(a)",' '
    print "(a,f8.2,a)",'  Pulsation period           = ',pulsation_period_days,' days'
    print "(a,f6.3)",  '  Pulsation amplitude        = ',pulsation_amplitude
    print "(a)",'================================================================'
    print "(a)",' '
    print "(a)",'NOTE: Atmosphere will be created on first timestep by inject module'
    print "(a,f6.2,a)",'      Atmosphere mass = ',atmos_mass_fraction*Mstar_Msun,' M_sun'
    print "(a,f6.2,a)",'      Sink mass after = ',(1.-atmos_mass_fraction)*Mstar_Msun,' M_sun'
    print "(a)",' '
 endif

end subroutine setpart

!----------------------------------------------------------------
!+
!  Write setup parameters to input file
!+
!----------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer :: iunit
 
 print "(a)",' Writing setup options file '//trim(filename)
 open(newunit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# Setup file for pulsating AGB star with wind'
 write(iunit,"(a)") '# Runtime parameters for this setup'
 
 write(iunit,"(/,a)") '# Stellar parameters'
 call write_inopt(Mstar_Msun,'Mstar_Msun','stellar mass (M_sun)',iunit)
 call write_inopt(Rstar_au,'Rstar_au','stellar radius (AU)',iunit)
 call write_inopt(Tstar,'Tstar','effective temperature (K)',iunit)
 
 write(iunit,"(/,a)") '# Atmosphere structure'
 call write_inopt(atmos_mass_fraction,'atmos_mass_fraction',&
      'atmospheric mass / total mass',iunit)
 call write_inopt(n_shells_total,'n_shells_total',&
      'total number of atmospheric shells',iunit)
 call write_inopt(iboundary_spheres,'iboundary_spheres',&
      'number of inner boundary shells',iunit)
 call write_inopt(r_min_on_rstar,'r_min_on_rstar',&
      'inner atmosphere radius / R_star',iunit)
 call write_inopt(iwind_resolution,'iwind_resolution',&
      'geodesic sphere resolution',iunit)
 
 write(iunit,"(/,a)") '# Pulsation parameters'
 call write_inopt(pulsation_period_days,'pulsation_period_days',&
      'pulsation period (days)',iunit)
 call write_inopt(pulsation_amplitude,'pulsation_amplitude',&
      'fractional pulsation amplitude (DeltaR/R)',iunit)
 
 close(iunit)
 
end subroutine write_setupfile

!----------------------------------------------------------------
!+
!  Read setup parameters from input file
!+
!----------------------------------------------------------------
subroutine read_setupfile(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 use io,           only:error,fatal
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 21
 type(inopts), allocatable :: db(:)
 
 print "(a)",' Reading setup options from '//trim(filename)
 call open_db_from_file(db,filename,iunit,ierr)
 
 ! Stellar parameters
 call read_inopt(Mstar_Msun,'Mstar_Msun',db,ierr)
 call read_inopt(Rstar_au,'Rstar_au',db,ierr)
 call read_inopt(Tstar,'Tstar',db,ierr)
 
 ! Atmosphere structure
 call read_inopt(atmos_mass_fraction,'atmos_mass_fraction',db,ierr)
 call read_inopt(n_shells_total,'n_shells_total',db,ierr)
 call read_inopt(iboundary_spheres,'iboundary_spheres',db,ierr)
 call read_inopt(r_min_on_rstar,'r_min_on_rstar',db,ierr)
 call read_inopt(iwind_resolution,'iwind_resolution',db,ierr)
 
 ! Pulsation parameters
 call read_inopt(pulsation_period_days,'pulsation_period_days',db,ierr)
 call read_inopt(pulsation_amplitude,'pulsation_amplitude',db,ierr)
 
 call close_db(db)
 
 ! Sanity checks
 if (Mstar_Msun <= 0.) call fatal('setup','Mstar_Msun must be > 0')
 if (Rstar_au <= 0.) call fatal('setup','Rstar_au must be > 0')
 if (Tstar <= 0.) call fatal('setup','Tstar must be > 0')
 if (atmos_mass_fraction <= 0. .or. atmos_mass_fraction >= 1.) &
    call fatal('setup','atmos_mass_fraction must be in (0,1)')
 if (n_shells_total <= 0) call fatal('setup','n_shells_total must be > 0')
 if (iboundary_spheres < 0) call fatal('setup','iboundary_spheres must be >= 0')
 if (iboundary_spheres > n_shells_total) &
    call fatal('setup','iboundary_spheres must be <= n_shells_total')
 if (r_min_on_rstar <= 0. .or. r_min_on_rstar >= 1.) &
    call fatal('setup','r_min_on_rstar must be in (0,1)')
 if (pulsation_period_days < 0.) &
    call fatal('setup','pulsation_period_days must be >= 0')
 if (pulsation_amplitude < 0.) &
    call fatal('setup','pulsation_amplitude must be >= 0')
 
end subroutine read_setupfile

end module setup