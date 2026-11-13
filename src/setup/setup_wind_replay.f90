!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Setup of the outer run of a periodic track replay (see inject_replay):
! no gas, only the sink(s) of the inner run at the time of the first dump
! of the replay period, so that the injection starts at the beginning of
! a replay cycle and fills the domain from r_inj outwards. On the first
! setup the physics options are copied from the .in file of the inner run.
!
! :References: None
!
! :Owner: Owen Vermeulen
!
! :Runtime parameters:
!   - merge_sinks  : *merge the sinks into one at their centre of mass (0=no, 1=yes)*
!   - r_inj        : *injection radius of the outer run (au)*
!   - r_out        : *outer radius rkill of this run (au)*
!   - replay_dumps : *first full dump of the replay period (path as seen from this run)*
!
! :Dependencies: infile_utils, inject, injection, io, kernel, part,
!   physcon, prompting, readwrite_infile, replayutils, units
!
 implicit none
 public :: setpart

 private

 character(len=200) :: replay_dumps    = '../init_run/wind_00000'
 real               :: r_inj_au        = 30.
 real               :: r_out_au        = 1000.
 integer            :: merge_sinks_opt = 1

contains

!----------------------------------------------------------------
!+
!  Setup of the outer run
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use part,             only:xyzmh_ptmass,vxyz_ptmass,nptmass,igas,iboundary
 use units,            only:set_units,udist,utime
 use physcon,          only:au
 use kernel,           only:hfact_default
 use inject,           only:set_replay_options
 use replayutils,      only:inner_dump,read_inner_dump,inner_infile,merge_sinks
 use injection,        only:rkill
 use readwrite_infile, only:read_infile
 use io,               only:fatal,warning,master
 integer,          intent(in)    :: id
 integer,          intent(inout) :: npart
 integer,          intent(out)   :: npartoftype(:)
 real,             intent(out)   :: xyzh(:,:)
 real,             intent(out)   :: massoftype(:)
 real,             intent(out)   :: polyk,gamma,hfact
 real,             intent(inout) :: time
 character(len=*), intent(in)    :: fileprefix
 real,             intent(out)   :: vxyzu(:,:)
 type(inner_dump)   :: inner
 character(len=200) :: setupfile,inner_in,logf,evf,dumpf
 integer :: ierr,nsinks
 logical :: iexist

 setupfile = trim(fileprefix)//'.setup'
 inquire(file=trim(setupfile),exist=iexist)
 if (iexist) then
    call read_setupfile(setupfile,ierr)
    if (ierr /= 0) call fatal('setup','errors in '//trim(setupfile))
 elseif (id == master) then
    call interactive_setup()
    call write_setupfile(setupfile)
 endif

 ! the first dump of the replay period
 call read_inner_dump(trim(replay_dumps),inner,ierr)
 if (ierr /= 0) call fatal('setup','cannot read '//trim(replay_dumps)//' (set replay_dumps in '//trim(setupfile)//')')
 if (inner%nptmass < 1) call fatal('setup','no sink particles in '//trim(replay_dumps))

 ! physics options of the inner run, on the first setup only (afterwards the .in of this run is kept)
 inquire(file=trim(fileprefix)//'.in',exist=iexist)
 if (.not. iexist) then
    inner_in = inner_infile(replay_dumps)
    inquire(file=trim(inner_in),exist=iexist)
    if (iexist) then
       print *, 'Copying the physics options from '//trim(inner_in)
       call read_infile(inner_in,logf,evf,dumpf)
    else
       call warning('setup',trim(inner_in)//' not found: check the physics options in '//trim(fileprefix)//'.in')
    endif
 endif

 call set_units(dist=inner%udist,mass=inner%umass,G=1.)
 if (abs(utime/inner%utime - 1.) > 1.e-6) call fatal('setup','the inner run does not use units with G = 1')

 ! the sink(s) of the inner run
 nsinks  = inner%nptmass
 nptmass = nsinks
 xyzmh_ptmass(:,:)        = 0.
 vxyz_ptmass(:,:)         = 0.
 xyzmh_ptmass(:,1:nsinks) = inner%xyzmh_ptmass(:,1:nsinks)
 vxyz_ptmass(:,1:nsinks)  = inner%vxyz_ptmass(:,1:nsinks)
 if (nsinks > 1 .and. merge_sinks_opt == 1) then
    call merge_sinks(nptmass,xyzmh_ptmass,vxyz_ptmass)
    print *, 'Merged ',nsinks,' sinks into one of mass ',xyzmh_ptmass(4,1)
 endif

 ! no gas: the injection fills the domain from r_inj outwards
 npart          = 0
 npartoftype(:) = 0
 xyzh(:,:)      = 0.
 vxyzu(:,:)     = 0.
 massoftype(igas)      = inner%pmass
 massoftype(iboundary) = inner%pmass
 gamma = 5./3.
 if (inner%gamma > 1.) gamma = inner%gamma
 polyk = inner%polyk
 hfact = hfact_default
 time  = inner%time

 call set_replay_options(trim(replay_dumps),r_inj_au)
 rkill = r_out_au*au/udist

 print *, ''
 print *, ' Outer run of a periodic track replay'
 print *, ' Start time (time of the dump) :', time
 print *, ' Sinks                         :', nptmass
 print *, ' r_inj (au)                    :', r_inj_au
 print *, ' rkill (au)                    :', r_out_au
 print *, ''

end subroutine setpart

!----------------------------------------------------------------
!+
!  Ask for the options (first setup)
!+
!----------------------------------------------------------------
subroutine interactive_setup()
 use prompting, only:prompt

 call prompt('first full dump of the replay period (path as seen from this run)',replay_dumps)
 call prompt('injection radius r_inj (au)',r_inj_au,0.)
 call prompt('outer radius rkill of this run (au)',r_out_au,r_inj_au)
 call prompt('merge the sinks into one at their centre of mass (0=no, 1=yes)',merge_sinks_opt,0,1)

end subroutine interactive_setup

!----------------------------------------------------------------
!+
!  Write the setup file
!+
!----------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer :: iunit

 print *, 'Writing setup options file '//trim(filename)
 open(newunit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for the outer run of a periodic track replay'
 call write_inopt(trim(replay_dumps), 'replay_dumps', 'first full dump of the replay period (path as seen from this run)',iunit)
 call write_inopt(r_inj_au,           'r_inj',        'injection radius of the outer run (au)',iunit)
 call write_inopt(r_out_au,           'r_out',        'outer radius rkill of this run (au)',iunit)
 call write_inopt(merge_sinks_opt,    'merge_sinks',  'merge the sinks into one at their centre of mass (0=no, 1=yes)',iunit)
 close(iunit)

end subroutine write_setupfile

!----------------------------------------------------------------
!+
!  Read the setup file
!+
!----------------------------------------------------------------
subroutine read_setupfile(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 use io,           only:ireadin
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 type(inopts), allocatable :: db(:)
 integer :: nerr

 nerr = 0
 print *, 'Reading setup options from '//trim(filename)
 call open_db_from_file(db,filename,ireadin,ierr)
 if (ierr /= 0) return
 call read_inopt(replay_dumps,    'replay_dumps', db,errcount=nerr)
 call read_inopt(r_inj_au,        'r_inj',        db,errcount=nerr,min=0.)
 call read_inopt(r_out_au,        'r_out',        db,errcount=nerr,min=0.)
 call read_inopt(merge_sinks_opt, 'merge_sinks',  db,errcount=nerr,min=0,max=1)
 call close_db(db)
 ierr = nerr

end subroutine read_setupfile

end module setup
