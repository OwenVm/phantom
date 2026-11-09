!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Setup of the outer run of a periodic track replay (see inject_replay):
! no gas, only the sink(s) of the inner run, at the time of the first dump
! of the replay period, so that the injection starts at the beginning of
! a replay cycle and fills the domain from r_inj outwards. On the first
! setup the physics options are copied from the .in file of the inner run.
!
! :References: None
!
! :Owner: Owen Vermeulen
!
! :Runtime parameters:
!   - buffer_width  : *width of the replay buffer r_inj - r_b (au, 0=auto from h)*
!   - merge_sinks   : *merge the sinks into one at their centre of mass (0=no, 1=yes)*
!   - r_inj         : *interface radius: particles are released here (au)*
!   - r_out         : *outer radius rkill of this run (au)*
!   - replay_dumps  : *first full dump of the replay period (path as seen from this run)*
!   - replay_period : *replay period = pulsation or orbital period (days)*
!
! :Dependencies: infile_utils, inject, injection, io, kernel, part,
!   physcon, prompting, readwrite_infile, replayutils, units
!
 implicit none
 public :: setpart

 private

 character(len=100) :: replay_dumps = '../init_run/wind_00000'
 real    :: r_inj_au    = 30.
 real    :: period_days = 300.
 real    :: buffer_au   = 0.
 real    :: r_out_au    = 1000.
 integer :: merge_sinks = 1

contains

subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use part,             only:xyzmh_ptmass,vxyz_ptmass,nptmass,igas,iboundary,ihacc,ihsoft,ilum,iTeff,iReff
 use units,            only:set_units,udist,utime
 use physcon,          only:au
 use kernel,           only:hfact_default
 use inject,           only:set_replay_options,read_inner_sinks
 use injection,        only:rkill
 use readwrite_infile, only:read_infile
 use io,               only:fatal,warning,master
 use replayutils,      only:isk_m,isk_lum,isk_teff,isk_reff,isk_hacc
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=*),  intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 character(len=120) :: setupfile,inner_in,logf,evf,dumpf
 real, allocatable  :: sinks(:,:),hsoft(:)
 real    :: t0,ud,um,ut,pm,gam,pk,mtot,xc(3),vc(3),ltot
 integer :: ierr,ns,k,iloc
 logical :: iexist_setup,iexist_in,iexist

 setupfile = trim(fileprefix)//'.setup'
 inquire(file=trim(setupfile),exist=iexist_setup)
 inquire(file=trim(fileprefix)//'.in',exist=iexist_in)
 if (iexist_setup) then
    call read_setupfile(setupfile,ierr)
    if (ierr /= 0) call fatal('setup','errors in '//trim(setupfile))
 elseif (id == master) then
    call interactive_setup()
    call write_setupfile(setupfile)
 endif

 ! sinks, time, units and particle mass of the first dump of the replay period
 call read_inner_sinks(trim(replay_dumps),t0,ud,um,ut,pm,gam,pk,ns,sinks,hsoft,ierr)
 if (ierr /= 0) call fatal('setup','cannot read '//trim(replay_dumps)//' (set replay_dumps in '// &
                           trim(setupfile)//')')
 if (ns < 1) call fatal('setup','no sink particles in '//trim(replay_dumps))

 ! physics options of the inner run, on the first setup only (afterwards the .in of this run is kept)
 if (.not.iexist_in) then
    iloc = index(replay_dumps,'_',back=.true.)
    inner_in = replay_dumps(1:max(iloc-1,1))//'.in'
    inquire(file=trim(inner_in),exist=iexist)
    if (iexist) then
       print "(a)",' copying the physics options from '//trim(inner_in)
       call read_infile(inner_in,logf,evf,dumpf)
    else
       call warning('setup',trim(inner_in)//' not found: check the physics options in '//trim(fileprefix)//'.in')
    endif
 endif

 call set_units(dist=ud,mass=um,G=1.)
 if (abs(utime/ut - 1.) > 1.e-6) call fatal('setup','the inner run does not use units with G = 1')

 ! the sink(s), optionally merged into one at their centre of mass
 xyzmh_ptmass(:,:) = 0.
 vxyz_ptmass(:,:)  = 0.
 nptmass = ns
 do k = 1,ns
    xyzmh_ptmass(1:3,k)    = sinks(1:3,k)
    vxyz_ptmass(1:3,k)     = sinks(4:6,k)
    xyzmh_ptmass(4,k)      = sinks(isk_m,k)
    xyzmh_ptmass(ihacc,k)  = sinks(isk_hacc,k)
    xyzmh_ptmass(ihsoft,k) = hsoft(k)
    xyzmh_ptmass(ilum,k)   = sinks(isk_lum,k)
    xyzmh_ptmass(iTeff,k)  = sinks(isk_teff,k)
    xyzmh_ptmass(iReff,k)  = sinks(isk_reff,k)
 enddo
 if (ns > 1 .and. merge_sinks == 1) then
    mtot = sum(xyzmh_ptmass(4,1:ns))
    ltot = sum(xyzmh_ptmass(ilum,1:ns))
    do k = 1,3
       xc(k) = sum(xyzmh_ptmass(4,1:ns)*xyzmh_ptmass(k,1:ns))/mtot
       vc(k) = sum(xyzmh_ptmass(4,1:ns)*vxyz_ptmass(k,1:ns))/mtot
    enddo
    xyzmh_ptmass(1:3,1)  = xc
    vxyz_ptmass(1:3,1)   = vc
    xyzmh_ptmass(4,1)    = mtot
    xyzmh_ptmass(ilum,1) = ltot
    xyzmh_ptmass(:,2:ns) = 0.
    vxyz_ptmass(:,2:ns)  = 0.
    nptmass = 1
    print "(a,i2,a,es12.4,a)",' merged ',ns,' sinks into one of mass ',mtot,' (code units)'
 endif

 ! no gas: the injection fills the domain from r_inj outwards
 npart          = 0
 npartoftype(:) = 0
 xyzh(:,:)      = 0.
 vxyzu(:,:)     = 0.
 massoftype(igas)      = pm
 massoftype(iboundary) = pm
 gamma = 5./3.
 if (gam > 1.) gamma = gam
 polyk = pk
 hfact = hfact_default
 time  = t0

 call set_replay_options(trim(replay_dumps),period_days,r_inj_au,buffer_au)
 rkill = r_out_au*au/udist

 print "(a)",' outer run of a periodic track replay'
 print "(a,es12.4,a)",' starts at t = ',time,' (time of '//trim(replay_dumps)//')'
 print "(a,i3,a,es12.4,a,es12.4,a)",' sinks: ',nptmass,', r_inj = ',r_inj_au,' au, rkill = ',r_out_au,' au'

end subroutine setpart

!-----------------------------------------------------------------------
!+
!  Ask for the options (first setup)
!+
!-----------------------------------------------------------------------
subroutine interactive_setup()
 use prompting, only:prompt

 call prompt('first full dump of the replay period (path as seen from this run)',replay_dumps)
 call prompt('interface radius r_inj (au)',r_inj_au,0.)
 call prompt('replay period = pulsation or orbital period (days)',period_days,0.)
 call prompt('width of the replay buffer (au, 0=auto from h)',buffer_au,0.)
 call prompt('outer radius rkill of this run (au)',r_out_au,r_inj_au)
 call prompt('merge the sinks into one at their centre of mass (0=no, 1=yes)',merge_sinks,0,1)

end subroutine interactive_setup

!-----------------------------------------------------------------------
!+
!  Write / read the setup file
!+
!-----------------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer :: iunit

 print "(a)",' writing setup options file '//trim(filename)
 open(newunit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for the outer run of a periodic track replay'
 call write_inopt(trim(replay_dumps),'replay_dumps','first full dump of the replay period (path as seen from this run)',iunit)
 call write_inopt(r_inj_au,'r_inj','interface radius: particles are released here (au)',iunit)
 call write_inopt(period_days,'replay_period','replay period = pulsation or orbital period (days)',iunit)
 call write_inopt(buffer_au,'buffer_width','width of the replay buffer r_inj - r_b (au, 0=auto from h)',iunit)
 call write_inopt(r_out_au,'r_out','outer radius rkill of this run (au)',iunit)
 call write_inopt(merge_sinks,'merge_sinks','merge the sinks into one at their centre of mass (0=no, 1=yes)',iunit)
 close(iunit)

end subroutine write_setupfile

subroutine read_setupfile(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 type(inopts), allocatable :: db(:)
 integer :: nerr

 nerr = 0
 print "(a)",' reading setup options from '//trim(filename)
 call open_db_from_file(db,filename,21,ierr)
 if (ierr /= 0) return
 call read_inopt(replay_dumps,'replay_dumps',db,errcount=nerr)
 call read_inopt(r_inj_au,'r_inj',db,errcount=nerr,min=0.)
 call read_inopt(period_days,'replay_period',db,errcount=nerr,min=0.)
 call read_inopt(buffer_au,'buffer_width',db,errcount=nerr,min=0.)
 call read_inopt(r_out_au,'r_out',db,errcount=nerr,min=0.)
 call read_inopt(merge_sinks,'merge_sinks',db,errcount=nerr,min=0,max=1)
 call close_db(db)
 ierr = nerr

end subroutine read_setupfile

end module setup
