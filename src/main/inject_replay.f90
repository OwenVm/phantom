!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module inject
!
! Injection by periodic replay of the particle tracks of an inner wind
! simulation: an outer run is stitched onto the inner one at r_inj.
! The track table and the replay itself are in utils_replay.
!
! :References: None
!
! :Owner: Owen Vermeulen
!
! :Runtime parameters:
!   - irotate         : *per-cycle transformation (0=none, 1=random rotation, 2=z-mirror)*
!   - r_inj           : *injection radius of the outer run (au)*
!   - replay_dumps    : *first full dump of the inner run of the replay period*
!   - replay_period   : *replay period (code units), set by the run*
!   - use_inner_dtmax : *dtmax = time between the dumps of the inner run (0=no, 1=yes)*
!
! :Dependencies: dim, dynamic_dtmax, infile_utils, io, part, partinject,
!   physcon, replayutils, timestep, units
!
 use io, only:fatal,warning
 implicit none
 character(len=*), parameter, public :: inject_type = 'replay'

 public :: init_inject, inject_particles, write_options_inject, read_options_inject, &
           set_default_options_inject, update_injected_par
 public :: set_replay_options   ! used by the setup of the outer run

 private

 character(len=200) :: replay_dumps    = 'wind_00000'
 real               :: r_inj_au        = 50.
 integer            :: irotate         = 0
 integer            :: use_inner_dtmax = 0

 character(len=*), parameter :: label = 'inject_replay'

contains

subroutine set_default_options_inject(flag)
 integer, optional, intent(in) :: flag

 replay_dumps    = 'wind_00000'
 r_inj_au        = 50.
 irotate         = 0
 use_inner_dtmax = 0

end subroutine set_default_options_inject

!----------------------------------------------------------------
!+
!  Set the options from the setup
!+
!----------------------------------------------------------------
subroutine set_replay_options(dumps,r_inj_in_au)
 character(len=*), intent(in) :: dumps
 real,             intent(in) :: r_inj_in_au

 replay_dumps = dumps
 r_inj_au     = r_inj_in_au

end subroutine set_replay_options

!----------------------------------------------------------------
!+
!  Initialize everything
!+
!----------------------------------------------------------------
subroutine init_inject(ierr)
 use physcon,       only:au
 use units,         only:udist
 use part,          only:massoftype,igas,iboundary,maxphase,maxp
 use dim,           only:mpi
 use timestep,      only:dtmax
 use dynamic_dtmax, only:get_dtmax_initial,dtmax_user
 use replayutils,   only:init_replay,particle_mass,n_dumps,n_sinks,dump_time
 integer, intent(out) :: ierr

 ierr = 0
 if (mpi) call fatal(label,'replay injection does not work with MPI')
 if (maxphase /= maxp) call fatal(label,'replay injection needs particle types (iphase)')

 ! the track table, read from file or built from the dumps of the inner run
 call init_replay(replay_dumps,r_inj_au*au/udist)

 ! the injected particles have the mass of those of the inner run
 if (massoftype(igas) > 0. .and. abs(massoftype(igas)/particle_mass - 1.) > 1.e-6) &
    call fatal(label,'particle mass differs from that of the inner run',var='massoftype',val=massoftype(igas))
 massoftype(igas)      = particle_mass
 massoftype(iboundary) = particle_mass

 ! dumps at the same times as those of the inner run: dtmax is the mean time between them
 if (use_inner_dtmax == 1 .and. n_dumps > 1) then
    dtmax = (dump_time(n_dumps) - dump_time(1))/(n_dumps - 1)
    if (any(abs(dump_time(2:n_dumps) - dump_time(1:n_dumps-1) - dtmax) > 0.01*dtmax)) &
       call warning(label,'the dumps of the inner run are not evenly spaced: dtmax is their mean spacing')
    call get_dtmax_initial(dtmax)
    print *, 'dtmax = time between the dumps of the inner run:',dtmax_user
 endif

 if (irotate == 1 .and. n_sinks > 1) call warning(label,'irotate=1 breaks the spiral pattern of a binary')

end subroutine init_inject

!----------------------------------------------------------------
!+
!  The actual function that is called by phantom
!+
!----------------------------------------------------------------
subroutine inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npart_old,npartoftype,dtinject)
 use part,        only:nptmass
 use partinject,  only:updated_particle
 use replayutils, only:dt_replay,window_start,replay_period,replay_started,get_centre_of_mass,start_replay, &
                       buffer_is_intact,rebuild_buffer,replay_tracks,kill_inner_gas,replay_luminosity
 real,    intent(in)    :: time,dtlast
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(inout) :: npart,npart_old
 integer, intent(inout) :: npartoftype(:)
 real,    intent(out)   :: dtinject
 real    :: xc(3),vc(3)
 integer :: ncycle

 dtinject         = dt_replay
 updated_particle = .true.   ! so that update_injected_particles applies dtinject on every step

 ! the tracks are relative to the centre of mass of the sinks, and repeat with the replay period
 call get_centre_of_mass(nptmass,xyzmh_ptmass,vxyz_ptmass,xc,vc)
 ncycle = floor((time - window_start)/replay_period)

 if (.not. replay_started) call start_replay(time,ncycle,npart,npartoftype,xyzh)
 if (.not. buffer_is_intact(npart,xyzh)) call rebuild_buffer(time,irotate,xc,vc,npart,npartoftype,xyzh,vxyzu)
 call replay_tracks(time,ncycle,irotate,xc,vc,npart,npartoftype,xyzh,vxyzu,xyzmh_ptmass)
 call kill_inner_gas(xc,npart,npartoftype,xyzh)
 call replay_luminosity(time,xyzmh_ptmass)

end subroutine inject_particles

!----------------------------------------------------------------
!+
!  Nothing to save between runs: the buffer is rebuilt from the
!  track table
!+
!----------------------------------------------------------------
subroutine update_injected_par
end subroutine update_injected_par

!----------------------------------------------------------------
!+
!  Write options to .in file
!+
!----------------------------------------------------------------
subroutine write_options_inject(iunit)
 use infile_utils, only:write_inopt
 use replayutils,  only:replay_period
 integer, intent(in) :: iunit

 call write_inopt(trim(replay_dumps), 'replay_dumps',    'first full dump of the inner run of the replay period',iunit)
 call write_inopt(r_inj_au,           'r_inj',           'injection radius of the outer run (au)',iunit)
 call write_inopt(irotate,            'irotate',         'per-cycle transformation (0=none, 1=random rotation, 2=z-mirror)',iunit)
 call write_inopt(use_inner_dtmax,    'use_inner_dtmax', 'dtmax = time between the dumps of the inner run (0=no, 1=yes)',iunit)
 ! not an input: the period found at the start of the run, read by a replay nested in this run
 call write_inopt(replay_period,      'replay_period',   'replay period (code units), set by the run',iunit)

end subroutine write_options_inject

!----------------------------------------------------------------
!+
!  Read options from .in file (all with defaults, so that the
!  setup can read the .in file of the inner run)
!+
!----------------------------------------------------------------
subroutine read_options_inject(db,nerr)
 use infile_utils, only:inopts,read_inopt
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr
 real    :: period_in_file
 integer :: ierr

 call read_inopt(replay_dumps,    'replay_dumps',    db,err=ierr)
 if (ierr /= 0 .or. len_trim(replay_dumps) == 0) replay_dumps = 'wind_00000'
 call read_inopt(r_inj_au,        'r_inj',           db,errcount=nerr,min=0.,default=50.)
 call read_inopt(irotate,         'irotate',         db,errcount=nerr,min=0,max=2,default=0)
 call read_inopt(use_inner_dtmax, 'use_inner_dtmax', db,errcount=nerr,min=0,max=1,default=0)
 ! written by write_options_inject, not an input: read only so that it is not an unknown variable
 call read_inopt(period_in_file,  'replay_period',   db,err=ierr)

end subroutine read_options_inject

end module inject
