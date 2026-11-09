!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module inject
!
! Injection by periodic replay of the particle tracks at the interface
! of an inner wind simulation (see replayutils).
!
! At the start of the run the track table is built from one period of
! full dumps of the inner run (replay_dumps, then the following dumps up
! to replay_period later) and saved to replay_file; a restart reads that
! file back. The table holds the path of every particle through the
! buffer shell r_b < r < r_inj, clipped to the period.
!
! On every call, i.e. every step:
!  - instances whose entry time has passed are spawned as boundary
!    particles (at r_b, or anywhere in the buffer at the start of a cycle)
!  - all buffer particles are moved to their state on the recorded track
!    at the current time (position, velocity, u, h, dust, mu, gamma, alpha)
!  - instances whose track has crossed r_inj become gas; the others are
!    removed when their track ends
!  - gas that has fallen back inside r_kill_inner is removed
!
! Boundary particles are moved by the drift between two calls but are
! never kicked or accreted, so they follow the tangent of their track
! until the next call.
!
! :References: None
!
! :Owner: Owen Vermeulen
!
! :Runtime parameters:
!   - buffer_width   : *width of the replay buffer r_inj - r_b (au, 0=auto from h)*
!   - dt_replay      : *maximum timestep while replaying (periods)*
!   - ideplete_sink  : *remove released mass from sink 1 (0=no, 1=yes)*
!   - include_drop   : *replay tracks that leave the buffer inwards (0=no, 1=yes)*
!   - ireplay_lum    : *replay the recorded sink luminosity (0=no, 1=yes)*
!   - irotate        : *per-cycle transformation (0=none, 1=random rotation, 2=z-mirror)*
!   - r_inj          : *interface radius: particles are released here (au)*
!   - r_kill_inner   : *gas inside this radius is removed (au, <=0: r_b)*
!   - replay_dumps   : *first full dump of the inner run to build the table from*
!   - replay_file    : *track table (built from the dumps if absent)*
!   - replay_period  : *replay period = pulsation or orbital period (days)*
!   - t_epoch        : *time at which replay cycle 0 starts (code units, <0: time of the first dump)*
!   - verbose_replay : *print a summary every cycle (0=no, 1=yes)*
!
! :Dependencies: dim, dump_utils, fileutils, infile_utils, io, kernel,
!   part, partinject, physcon, readwrite_dumps_common, replayutils, units
!
 use replayutils, only:dp,track_table,replay_schedule
 implicit none
 character(len=*), parameter, public :: inject_type = 'replay'

 public :: init_inject,inject_particles,write_options_inject,read_options_inject, &
           set_default_options_inject,update_injected_par
 public :: set_replay_options,read_inner_sinks

 private

 character(len=200) :: replay_dumps    = 'wind_00000'
 real               :: replay_period   = 300.
 real               :: r_inj_au        = 50.
 real               :: buffer_width_au = 0.
 character(len=120) :: replay_file     = 'replay_tracks.dat'
 real               :: t_epoch         = -1.
 real               :: r_kill_inner_au = 0.
 real               :: dt_replay       = 0.02
 integer            :: irotate         = 0
 integer            :: ideplete_sink   = 1
 integer            :: ireplay_lum     = 0
 integer            :: include_drop    = 1
 integer            :: iverbose_replay = 1

 integer, parameter :: iseed_rotate = 20261108
 character(len=*), parameter :: label = 'inject_replay'

 type(track_table)     :: tab
 type(replay_schedule) :: sch
 logical  :: started = .false.
 real(dp) :: t_ep = 0._dp
 real     :: dt_rep = 0., r_kill2 = 0., pmass = 0.
 ! active instances: track, cycle, particle index and particle id
 integer  :: nact = 0
 integer,         allocatable :: act_k(:),act_n(:),act_i(:)
 integer(kind=8), allocatable :: act_id(:)
 ! position of the recorded quantities in the table
 integer  :: is_nuc = 0, is_mu = 0, is_gam = 0, is_alpha = 0, is_tdust = 0, is_abund = 0
 ! rotation matrices of recent cycles
 integer, parameter :: ncache = 8
 integer  :: cache_n(ncache) = -huge(1)
 real(dp) :: cache_rot(3,3,ncache)
 ! diagnostics
 integer  :: mcycle_last = -huge(1), nrel_c = 0, nend_c = 0, nkill_c = 0, nspawn_c = 0, nlost = 0

 ! one full dump of the inner run (only what the table needs)
 type :: inner_dump
    real(dp) :: time = 0._dp, udist = 1._dp, umass = 1._dp, utime = 1._dp, pmass = 0._dp
    real(dp) :: gamma = 0._dp, polyk = 0._dp
    integer  :: n = 0, nsink = 0
    real(dp),        allocatable :: x(:,:),v(:,:),u(:),h(:),scal(:,:),sink(:,:),hsoft(:)
    integer(kind=8), allocatable :: iorig(:)
    integer(kind=1), allocatable :: itype(:)
    logical,         allocatable :: found(:)
    logical  :: got_iorig = .false., got_itype = .false.
 end type inner_dump
 ! per-particle quantities looked for in the dumps, with their names in the table
 integer :: ncand = 0
 character(len=16), allocatable :: cand_label(:)

contains

subroutine set_default_options_inject(flag)
 integer, intent(in), optional :: flag

 replay_dumps    = 'wind_00000'
 replay_period   = 300.
 r_inj_au        = 50.
 buffer_width_au = 0.
 replay_file     = 'replay_tracks.dat'
 t_epoch         = -1.
 r_kill_inner_au = 0.
 dt_replay       = 0.02
 irotate         = 0
 ideplete_sink   = 1
 ireplay_lum     = 0
 include_drop    = 1
 iverbose_replay = 1

end subroutine set_default_options_inject

!-----------------------------------------------------------------------
!+
!  Set the replay options (used by moddump_replay)
!+
!-----------------------------------------------------------------------
subroutine set_replay_options(dumps,period_days,r_au,buffer_au)
 character(len=*), intent(in) :: dumps
 real,             intent(in) :: period_days,r_au,buffer_au

 replay_dumps    = dumps
 replay_period   = period_days
 r_inj_au        = r_au
 buffer_width_au = buffer_au
 t_epoch         = -1.

end subroutine set_replay_options

!-----------------------------------------------------------------------
!+
!  Read the table, or build it from the dumps of the inner run
!+
!-----------------------------------------------------------------------
subroutine init_inject(ierr)
 use io,          only:fatal,warning,iprint
 use units,       only:udist,umass,utime
 use physcon,     only:au,days
 use part,        only:massoftype,igas,iboundary,maxphase,maxp,nucleation_label,n_nucleation, &
                       abundance_label
 use dim,         only:mpi
 use replayutils, only:read_track_table,write_track_table,schedule_init,find_label,percentile
 integer, intent(out) :: ierr
 real(dp), allocatable :: ratio(:)
 real(dp) :: period,vmag,t_rec,r_inj,t0
 integer  :: j,m,istride,ierr2
 logical  :: iexist,rebuild

 ierr = 0
 if (mpi) call fatal(label,'replay injection does not work with MPI')
 if (maxphase /= maxp) call fatal(label,'replay injection needs particle types (iphase)')

 t_rec = replay_period*days/utime
 r_inj = r_inj_au*au/udist
 inquire(file=trim(replay_file),exist=iexist)
 rebuild = .not.iexist
 if (iexist) then
    call read_track_table(trim(replay_file),tab,ierr)
    if (ierr /= 0) call fatal(label,'cannot read '//trim(replay_file))
    if (abs(tab%r_inj/r_inj - 1._dp) > 1.e-6_dp .or. abs(tab%t_rec/t_rec - 1._dp) > 1.e-6_dp) then
       write(iprint,"(a)") ' replay: '//trim(replay_file)//' was built for another r_inj or period: rebuilding'
       rebuild = .true.
    else
       call read_dump_time(trim(replay_dumps),t0,ierr2)
       if (ierr2 == 0 .and. abs(t0 - tab%t0) > 1.e-6_dp*t_rec) then
          write(iprint,"(a)") ' replay: '//trim(replay_file)//' was built from other dumps: rebuilding'
          rebuild = .true.
       endif
    endif
 endif
 if (rebuild) then
    call build_table_from_dumps(r_inj,t_rec)
    call write_track_table(trim(replay_file),tab,ierr)
    if (ierr /= 0) call warning(label,'could not write '//trim(replay_file))
    ierr = 0
 else
    write(iprint,"(/,a)") ' replay: read the track table from '//trim(replay_file)
 endif
 if (tab%ntracks < 1) call fatal(label,'no tracks in the table')
 if (abs(tab%udist/udist - 1._dp) > 1.e-6_dp .or. abs(tab%umass/umass - 1._dp) > 1.e-6_dp .or. &
     abs(tab%utime/utime - 1._dp) > 1.e-6_dp) call fatal(label,'units differ from those of the inner run')

 pmass = real(tab%pmass)
 if (massoftype(igas) > 0.) then
    if (abs(massoftype(igas)/pmass - 1.) > 1.e-6) &
       call fatal(label,'particle mass differs from that of the inner run',var='massoftype',val=massoftype(igas))
 endif
 massoftype(igas)      = pmass
 massoftype(iboundary) = pmass

 if (t_epoch < 0.) t_epoch = real(tab%t0)
 t_ep = t_epoch
 call schedule_init(sch,tab,t_ep)

 if (r_kill_inner_au > 0.) then
    r_kill2 = (r_kill_inner_au*au/udist)**2
 else
    r_kill2 = real(tab%r_b)**2
 endif

 ! timestep cap: dt_replay periods, and the time to cross 0.1 h in the buffer
 period = tab%t_rec
 m = 0
 istride = max(1,tab%nsamp/20000)
 allocate(ratio(tab%nsamp/istride + 1))
 do j = 1,tab%nsamp,istride
    vmag = norm2(real(tab%vs(:,j),dp))
    if (vmag > 0._dp .and. tab%hs(j) > 0.) then
       m = m + 1
       ratio(m) = tab%hs(j)/vmag
    endif
 enddo
 dt_rep = real(dt_replay*period)
 if (m > 0) dt_rep = min(dt_rep,real(0.1_dp*percentile(m,ratio,0.05_dp)))

 ! where the recorded quantities are in the table
 is_nuc   = find_label(tab%label,'nuc:'//trim(adjustl(nucleation_label(1))))
 if (is_nuc > 0) then
    if (is_nuc + n_nucleation - 1 > tab%nscal) is_nuc = 0
 endif
 is_mu    = find_label(tab%label,'eos:mu')
 is_gam   = find_label(tab%label,'eos:gamma')
 is_alpha = find_label(tab%label,'alpha')
 is_tdust = find_label(tab%label,'Tdust')
 is_abund = find_label(tab%label,'abund:'//trim(adjustl(abundance_label(1))))

 started = .false.
 nact    = 0
 write(iprint,"(a,i10,a,i12,a)") ' replay: ',tab%ntracks,' tracks, ',tab%nsamp,' samples'
 write(iprint,"(a,2(es12.4,a))") ' replay: r_b = ',tab%r_b*udist/au,' au, r_inj = ',tab%r_inj*udist/au,' au'
 write(iprint,"(a,3(es12.4,a))") ' replay: t0 = ',tab%t0,', t_rec = ',tab%t_rec,', t_epoch = ',t_ep,' (code units)'
 write(iprint,"(a,es12.4,a)") ' replay: dtinject = ',dt_rep,' (code units)'
 if (is_nuc == 0) call warning(label,'no dust moments in the table')
 if (irotate == 1 .and. tab%nsink > 1) call warning(label,'irotate=1 breaks the spiral pattern of a binary')

end subroutine init_inject

!-----------------------------------------------------------------------
!+
!  Build the track table from one period of dumps of the inner run,
!  starting at replay_dumps
!+
!-----------------------------------------------------------------------
subroutine build_table_from_dumps(r_inj,t_rec)
 use io,          only:fatal,warning,iprint
 use units,       only:udist,umass,utime
 use physcon,     only:au,solarm,years,km
 use kernel,      only:radkern
 use part,        only:igas
 use fileutils,   only:getnextfilename
 use replayutils, only:replay_builder,replay_record,builder_init,builder_add_record,builder_finalise, &
                       alloc_record,record_stats,percentile,sp,i8,nsinkcol,isk_m
 real(dp), intent(in) :: r_inj,t_rec
 type(inner_dump)     :: d
 type(replay_builder) :: b
 type(replay_record)  :: rec
 type(track_table)    :: head
 character(len=200)   :: fname
 real(dp), allocatable :: hl(:)
 integer,  allocatable :: col(:)
 real(dp) :: t0,tprev,xc(3),vc(3),mtot,r,delta,r_b,band_in,band_out,h90
 real(dp) :: mout,mdin,vrf,h90r,hmean,s_mout,s_mdin,s_vrm,s_h,umdot,uvel
 integer  :: i,k,m,nh,ierr,ndumps,nstat,nsc
 logical  :: iexist

 call set_candidates()
 fname = replay_dumps
 write(iprint,"(/,a)") ' replay: building the track table from '//trim(fname)//' and the following dumps'
 call read_inner_dump(trim(fname),d,ierr)
 if (ierr /= 0) call fatal(label,'cannot read '//trim(fname))
 if (abs(d%udist/udist - 1._dp) > 1.e-6_dp .or. abs(d%umass/umass - 1._dp) > 1.e-6_dp .or. &
     abs(d%utime/utime - 1._dp) > 1.e-6_dp) call fatal(label,'units of '//trim(fname)//' differ from those of this run')
 if (.not.d%got_iorig) call fatal(label,trim(fname)//' has no particle ids (iorig)')
 t0 = d%time

 ! buffer width from h at r_inj in the first dump
 call centre(d,xc,vc)
 allocate(hl(max(d%n,1)))
 nh = 0
 do i = 1,d%n
    if (.not.is_gas(i)) cycle
    r = norm2(d%x(:,i) - xc)
    if (abs(r - r_inj) < 0.05_dp*r_inj) then
       nh = nh + 1
       hl(nh) = d%h(i)
    endif
 enddo
 if (buffer_width_au > 0.) then
    delta = buffer_width_au*au/udist
 else
    if (nh < 20) call fatal(label,'too few particles at r_inj in '//trim(fname)//' to set the buffer width')
    h90   = percentile(nh,hl,0.9_dp)
    delta = 1.2_dp*radkern*h90
 endif
 r_b      = r_inj - delta
 band_in  = r_b - 0.5_dp*delta
 band_out = r_inj + 0.5_dp*delta
 if (r_b <= 0._dp) call fatal(label,'buffer wider than r_inj',var='buffer width (au)',val=real(delta*udist/au))

 ! per-particle quantities found in the dump
 nsc = count(d%found)
 allocate(col(max(nsc,1)))
 head%nscal = nsc
 allocate(head%label(max(nsc,1)))
 k = 0
 do i = 1,ncand
    if (d%found(i)) then
       k = k + 1
       col(k) = i
       head%label(k) = cand_label(i)
    endif
 enddo
 head%udist  = udist
 head%umass  = umass
 head%utime  = utime
 head%pmass  = d%pmass
 head%r_inj  = r_inj
 head%r_b    = r_b
 head%t0     = t0
 head%t_rec  = t_rec
 head%period = t_rec
 call builder_init(b,head,0._dp,include_drop == 1,.true.)

 ndumps = 0
 nstat  = 0
 s_mout = 0._dp; s_mdin = 0._dp; s_vrm = 0._dp; s_h = 0._dp
 tprev  = t0
 do
    ! the gas in the band, relative to the centre of mass of the sinks
    call centre(d,xc,vc)
    m = 0
    do i = 1,d%n
       if (.not.is_gas(i)) cycle
       r = norm2(d%x(:,i) - xc)
       if (r >= band_in .and. r <= band_out) m = m + 1
    enddo
    call alloc_record(rec,d%nsink,m,nsc)
    rec%time  = d%time
    rec%xc    = xc
    rec%vc    = vc
    rec%nsink = d%nsink
    if (d%nsink > 0) rec%sink(1:nsinkcol,1:d%nsink) = d%sink(1:nsinkcol,1:d%nsink)
    m = 0
    do i = 1,d%n
       if (.not.is_gas(i)) cycle
       r = norm2(d%x(:,i) - xc)
       if (r < band_in .or. r > band_out) cycle
       m = m + 1
       rec%iorig(m) = int(d%iorig(i),i8)
       rec%x(:,m)   = real(d%x(:,i) - xc,sp)
       rec%v(:,m)   = real(d%v(:,i) - vc,sp)
       rec%u(m)     = real(d%u(i),sp)
       rec%h(m)     = real(d%h(i),sp)
       if (nsc > 0) rec%scal(1:nsc,m) = real(d%scal(col(1:nsc),i),sp)
    enddo
    rec%nsel = m
    call builder_add_record(b,rec)
    ndumps = ndumps + 1
    if (d%time < t0 + t_rec*(1._dp - 1.e-6_dp)) then
       call record_stats(rec,d%pmass,r_inj,0.5_dp*delta,mout,mdin,vrf,h90r,hmean)
       nstat  = nstat + 1
       s_mout = s_mout + mout
       s_mdin = s_mdin + mdin
       s_vrm  = s_vrm + vrf*mout
       s_h    = s_h + hmean
    endif
    if (d%time >= t0 + t_rec*(1._dp - 1.e-6_dp)) exit

    fname = getnextfilename(fname)
    inquire(file=trim(fname),exist=iexist)
    if (.not.iexist) call fatal(label,'the dumps end before one replay period: '//trim(fname)//' is missing')
    call read_inner_dump(trim(fname),d,ierr)
    if (ierr /= 0) call fatal(label,'cannot read '//trim(fname))
    if (d%time <= tprev) call fatal(label,trim(fname)//' is not later than the previous dump')
    tprev = d%time
 enddo
 call builder_finalise(b,tab)

 umdot = (umass/solarm)/(utime/years)
 uvel  = udist/utime/km
 write(iprint,"(a,i6,a,a)") ' replay: read ',ndumps,' dumps, up to ',trim(fname)
 write(iprint,"(a,2(es12.4,a))") ' replay: window t0 = ',t0,', t_rec = ',t_rec,' (code units)'
 write(iprint,"(a,2(es12.4,a))") ' replay: buffer from r_b = ',r_b*udist/au,' to r_inj = ',r_inj*udist/au,' au'
 if (nstat > 0 .and. s_mout > 0._dp) then
    write(iprint,"(a,es12.4,a)") ' replay: Mdot out at r_inj : ',s_mout/nstat*umdot,' Msun/yr (shell estimate)'
    write(iprint,"(a,es12.4)")   ' replay: Mdot in / out     : ',s_mdin/s_mout
    write(iprint,"(a,es12.4,a)") ' replay: v_r at r_inj      : ',s_vrm/s_mout*uvel,' km/s'
    if (s_mdin/s_mout > 0.05_dp) call warning(label,'more than 5% of the mass flux at r_inj falls back: '// &
                                                    'r_inj is too close to the star')
    if (ndumps > 1) then
       if ((s_vrm/s_mout)*t_rec/(ndumps - 1) > 0.5_dp*s_h/nstat) &
          call warning(label,'particles move more than h/2 between dumps: write dumps more often')
    endif
 endif
 write(iprint,"(a,i10)") ' replay: release tracks    : ',b%nrel
 write(iprint,"(a,i10)") ' replay: drop tracks       : ',b%ndrop
 write(iprint,"(a,i10,a)") ' replay: in buffer at t0   : ',b%nappear,' (appear at the start of every cycle)'
 write(iprint,"(a,i10,a)") ' replay: in buffer at end  : ',b%nvanish,' (removed at the end of every cycle)'
 write(iprint,"(a,i10)") ' replay: samples           : ',tab%nsamp
 write(iprint,"(a,es12.4,a)") ' replay: replayed Mdot     : ',b%nrel*d%pmass/t_rec*umdot,' Msun/yr'
 if (b%nlost > 0) write(iprint,"(a,i8,a)") ' replay: WARNING: ',b%nlost, &
       ' particles left the band between two dumps without crossing r_inj or r_b'

contains

 logical function is_gas(ip)
  integer, intent(in) :: ip

  is_gas = .true.
  if (d%got_itype) is_gas = (d%itype(ip) == igas)

 end function is_gas

 subroutine centre(dd,xcen,vcen)
  type(inner_dump), intent(in)  :: dd
  real(dp),         intent(out) :: xcen(3),vcen(3)

  xcen = 0._dp
  vcen = 0._dp
  mtot = 0._dp
  do k = 1,dd%nsink
     if (dd%sink(isk_m,k) <= 0._dp) cycle
     xcen = xcen + dd%sink(isk_m,k)*dd%sink(1:3,k)
     vcen = vcen + dd%sink(isk_m,k)*dd%sink(4:6,k)
     mtot = mtot + dd%sink(isk_m,k)
  enddo
  if (mtot > 0._dp) then
     xcen = xcen/mtot
     vcen = vcen/mtot
  endif

 end subroutine centre

end subroutine build_table_from_dumps

!-----------------------------------------------------------------------
!+
!  Names of the per-particle quantities looked for in the dumps, as they
!  are called in the table
!+
!-----------------------------------------------------------------------
subroutine set_candidates()
 use part, only:n_nucleation,nucleation_label,nabundances,abundance_label
 integer :: j

 ncand = n_nucleation + 4 + nabundances
 if (allocated(cand_label)) deallocate(cand_label)
 allocate(cand_label(ncand))
 do j = 1,n_nucleation
    cand_label(j) = 'nuc:'//trim(adjustl(nucleation_label(j)))
 enddo
 cand_label(n_nucleation+1) = 'eos:mu'
 cand_label(n_nucleation+2) = 'eos:gamma'
 cand_label(n_nucleation+3) = 'alpha'
 cand_label(n_nucleation+4) = 'Tdust'
 do j = 1,nabundances
    cand_label(n_nucleation+4+j) = 'abund:'//trim(adjustl(abundance_label(j)))
 enddo

end subroutine set_candidates

!-----------------------------------------------------------------------
!+
!  Index in cand_label of a particle array of the dump (0 if not wanted).
!  The nucleation arrays (Jstar, K0..K3, mu, gamma, S, kappa, alphw) are
!  written one after the other, so their 'mu' and 'gamma' are told apart
!  from those of eos_vars by their position in that sequence.
!+
!-----------------------------------------------------------------------
integer function cand_index(tag,nuc_pos)
 use part, only:n_nucleation,nucleation_label,nabundances,abundance_label
 character(len=*), intent(in)    :: tag
 integer,          intent(inout) :: nuc_pos
 integer :: j

 cand_index = 0
 if (trim(tag) == trim(adjustl(nucleation_label(1)))) then
    nuc_pos = 1
    cand_index = 1
    return
 elseif (nuc_pos > 0 .and. nuc_pos < n_nucleation) then
    if (trim(tag) == trim(adjustl(nucleation_label(nuc_pos+1)))) then
       nuc_pos = nuc_pos + 1
       cand_index = nuc_pos
       return
    endif
 endif
 nuc_pos = 0
 select case(trim(tag))
 case('mu')
    cand_index = n_nucleation + 1
 case('gamma')
    cand_index = n_nucleation + 2
 case('alpha')
    cand_index = n_nucleation + 3
 case('Tdust')
    cand_index = n_nucleation + 4
 case default
    do j = 1,nabundances
       if (trim(tag) == trim(adjustl(abundance_label(j)))) cand_index = n_nucleation + 4 + j
    enddo
 end select

end function cand_index

!-----------------------------------------------------------------------
!+
!  Read the arrays needed for the table from one full dump, in one pass
!+
!-----------------------------------------------------------------------
subroutine read_inner_dump(fname,d,ierr)
 use dump_utils,  only:open_dumpfile_r,read_header,free_header,extract,read_global_block_header, &
                       read_block_header,dump_h,lentag,lenid,ndatatypes,i_int1,i_int8,i_real,i_real4,i_real8
 use readwrite_dumps_common, only:get_options_from_fileid
 use part,        only:maxtypes,igas
 use replayutils, only:nsinkcol,isk_m,isk_lum,isk_teff,isk_reff,isk_hacc
 character(len=*), intent(in)    :: fname
 type(inner_dump), intent(inout) :: d
 integer,          intent(out)   :: ierr
 integer, parameter :: maxarraylengths = 12
 type(dump_h) :: hdr
 character(len=lenid)  :: fileid
 character(len=lentag) :: tag
 integer(kind=8) :: number8(maxarraylengths)
 integer :: nums(ndatatypes,maxarraylengths)
 integer :: iunit,nblocks,narraylengths,iblock,i,j,k,np,ns,ioff,isoff,ic,nuc_pos,npt,nsk,ierr1
 real(dp) :: mot(maxtypes),tm
 real(dp),        allocatable :: rbuf(:)
 real(kind=4),    allocatable :: r4buf(:)
 logical  :: tagged,phantomdump,smalldump,use_dust

 ierr  = 0
 iunit = free_unit()
 call open_dumpfile_r(iunit,fname,fileid,ierr,requiretags=.true.)
 if (ierr /= 0) return
 call get_options_from_fileid(fileid,tagged,phantomdump,smalldump,use_dust,ierr)
 if (smalldump .or. .not.tagged) then
    ierr = 1
    close(iunit)
    return
 endif
 call read_header(iunit,hdr,ierr,tagged=.true.)
 if (ierr /= 0) then
    close(iunit)
    return
 endif
 call extract('time',tm,hdr,ierr1)
 d%time = tm
 call extract('udist',d%udist,hdr,ierr1)
 call extract('umass',d%umass,hdr,ierr1)
 call extract('utime',d%utime,hdr,ierr1)
 mot = 0._dp
 call extract('massoftype',mot,hdr,ierr1)
 d%pmass = mot(igas)
 call extract('nparttot',npt,hdr,ierr1)
 call extract('nptmass',nsk,hdr,ierr1,default=0)
 d%gamma = 0._dp
 d%polyk = 0._dp
 call extract('gamma',tm,hdr,ierr1)
 if (ierr1 == 0) d%gamma = tm
 call extract('RK2',tm,hdr,ierr1)
 if (ierr1 == 0) d%polyk = tm/1.5_dp
 call read_global_block_header(nblocks,narraylengths,hdr,iunit,ierr)
 call free_header(hdr,ierr1)
 if (ierr /= 0) then
    close(iunit)
    return
 endif

 call alloc_dump(d,npt,nsk)
 ioff    = 0
 isoff   = 0
 nuc_pos = 0
 do iblock = 1,nblocks
    call read_block_header(narraylengths,number8,nums,iunit,ierr)
    if (ierr /= 0) exit
    np = int(number8(1))
    ns = 0
    if (narraylengths >= 2) ns = int(number8(2))
    do j = 1,narraylengths
       do i = 1,ndatatypes
          do k = 1,nums(i,j)
             read(iunit,iostat=ierr) tag
             if (ierr /= 0) exit
             if (j == 1 .and. np > 0) then
                call read_particle_array(i,np)
             elseif (j == 2 .and. ns > 0) then
                call read_sink_array(i,ns)
             else
                read(iunit,iostat=ierr)
             endif
             if (ierr /= 0) exit
          enddo
       enddo
    enddo
    ioff  = ioff + np
    isoff = isoff + ns
 enddo
 close(iunit)
 d%n     = ioff
 d%nsink = isoff
 if (ierr /= 0) return
 ! u is needed; h, x and v always exist in a full dump
 if (.not.allocated(d%u)) ierr = 2

contains

 subroutine read_particle_array(itype,n)
  integer, intent(in) :: itype,n

  select case(itype)
  case(i_int1)
     if (trim(tag) == 'itype') then
        read(iunit,iostat=ierr) d%itype(ioff+1:ioff+n)
        d%got_itype = .true.
     else
        read(iunit,iostat=ierr)
     endif
  case(i_int8)
     if (trim(tag) == 'iorig') then
        read(iunit,iostat=ierr) d%iorig(ioff+1:ioff+n)
        d%got_iorig = .true.
     else
        read(iunit,iostat=ierr)
     endif
  case(i_real,i_real8)
     ic = cand_index(tag,nuc_pos)
     select case(trim(tag))
     case('x','y','z','vx','vy','vz','u')
        call read_real8(n)
        select case(trim(tag))
        case('x');  d%x(1,ioff+1:ioff+n) = rbuf(1:n)
        case('y');  d%x(2,ioff+1:ioff+n) = rbuf(1:n)
        case('z');  d%x(3,ioff+1:ioff+n) = rbuf(1:n)
        case('vx'); d%v(1,ioff+1:ioff+n) = rbuf(1:n)
        case('vy'); d%v(2,ioff+1:ioff+n) = rbuf(1:n)
        case('vz'); d%v(3,ioff+1:ioff+n) = rbuf(1:n)
        case('u');  d%u(ioff+1:ioff+n)   = rbuf(1:n)
        end select
     case default
        if (ic > 0) then
           call read_real8(n)
           d%scal(ic,ioff+1:ioff+n) = rbuf(1:n)
           d%found(ic) = .true.
        else
           read(iunit,iostat=ierr)
        endif
     end select
  case(i_real4)
     ic = cand_index(tag,nuc_pos)
     if (trim(tag) == 'h') then
        call read_real4(n)
        d%h(ioff+1:ioff+n) = r4buf(1:n)
     elseif (ic > 0) then
        call read_real4(n)
        d%scal(ic,ioff+1:ioff+n) = r4buf(1:n)
        d%found(ic) = .true.
     else
        read(iunit,iostat=ierr)
     endif
  case default
     read(iunit,iostat=ierr)
  end select

 end subroutine read_particle_array

 subroutine read_sink_array(itype,n)
  integer, intent(in) :: itype,n
  integer :: icol

  icol = 0
  if (itype == i_real .or. itype == i_real8) then
     select case(trim(tag))
     case('x');    icol = 1
     case('y');    icol = 2
     case('z');    icol = 3
     case('vx');   icol = 4
     case('vy');   icol = 5
     case('vz');   icol = 6
     case('m');    icol = isk_m
     case('lum');  icol = isk_lum
     case('Teff'); icol = isk_teff
     case('Reff'); icol = isk_reff
     case('h');    icol = isk_hacc
     end select
  endif
  if (icol > 0) then
     call read_real8(n)
     d%sink(icol,isoff+1:isoff+n) = rbuf(1:n)
  elseif ((itype == i_real .or. itype == i_real8) .and. trim(tag) == 'hsoft') then
     call read_real8(n)
     d%hsoft(isoff+1:isoff+n) = rbuf(1:n)
  else
     read(iunit,iostat=ierr)
  endif

 end subroutine read_sink_array

 subroutine read_real8(n)
  integer, intent(in) :: n

  if (allocated(rbuf)) then
     if (size(rbuf) < n) deallocate(rbuf)
  endif
  if (.not.allocated(rbuf)) allocate(rbuf(n))
  read(iunit,iostat=ierr) rbuf(1:n)

 end subroutine read_real8

 subroutine read_real4(n)
  integer, intent(in) :: n

  if (allocated(r4buf)) then
     if (size(r4buf) < n) deallocate(r4buf)
  endif
  if (.not.allocated(r4buf)) allocate(r4buf(n))
  read(iunit,iostat=ierr) r4buf(1:n)

 end subroutine read_real4

end subroutine read_inner_dump

subroutine alloc_dump(d,n,ns)
 use replayutils, only:nsinkcol
 type(inner_dump), intent(inout) :: d
 integer,          intent(in)    :: n,ns

 if (allocated(d%x)) then
    if (size(d%x,2) < n .or. size(d%sink,2) < max(ns,1)) then
       deallocate(d%x,d%v,d%u,d%h,d%scal,d%iorig,d%itype,d%sink,d%hsoft)
    endif
 endif
 if (.not.allocated(d%x)) then
    allocate(d%x(3,max(n,1)),d%v(3,max(n,1)),d%u(max(n,1)),d%h(max(n,1)),d%scal(ncand,max(n,1)), &
             d%iorig(max(n,1)),d%itype(max(n,1)),d%sink(nsinkcol,max(ns,1)),d%hsoft(max(ns,1)))
 endif
 if (.not.allocated(d%found)) allocate(d%found(ncand))
 d%found     = .false.
 d%got_iorig = .false.
 d%got_itype = .false.
 d%sink      = 0._dp
 d%hsoft     = 0._dp
 d%u         = 0._dp

end subroutine alloc_dump

!-----------------------------------------------------------------------
!+
!  Sinks, time, units, particle mass, gamma and polyk of a dump of the
!  inner run (used by the setup of the outer run). Sink columns as in
!  replayutils (x,y,z,vx,vy,vz,m,lum,Teff,Reff,hacc).
!+
!-----------------------------------------------------------------------
subroutine read_inner_sinks(fname,time,udist_d,umass_d,utime_d,pmass_d,gamma_d,polyk_d,nsink,sinks,hsoft,ierr)
 use replayutils, only:nsinkcol
 character(len=*),  intent(in)  :: fname
 real,              intent(out) :: time,udist_d,umass_d,utime_d,pmass_d,gamma_d,polyk_d
 integer,           intent(out) :: nsink,ierr
 real, allocatable, intent(out) :: sinks(:,:),hsoft(:)
 type(inner_dump) :: d

 call set_candidates()
 call read_inner_dump(fname,d,ierr)
 if (ierr /= 0) return
 time    = real(d%time)
 udist_d = real(d%udist)
 umass_d = real(d%umass)
 utime_d = real(d%utime)
 pmass_d = real(d%pmass)
 gamma_d = real(d%gamma)
 polyk_d = real(d%polyk)
 nsink   = d%nsink
 allocate(sinks(nsinkcol,max(nsink,1)),hsoft(max(nsink,1)))
 sinks = 0.
 hsoft = 0.
 if (nsink > 0) then
    sinks(:,1:nsink) = real(d%sink(:,1:nsink))
    hsoft(1:nsink)   = real(d%hsoft(1:nsink))
 endif

end subroutine read_inner_sinks

!-----------------------------------------------------------------------
!+
!  Time of a dump (from its header)
!+
!-----------------------------------------------------------------------
subroutine read_dump_time(fname,time,ierr)
 use dump_utils, only:open_dumpfile_r,read_header,free_header,extract,dump_h,lenid
 character(len=*), intent(in)  :: fname
 real(dp),         intent(out) :: time
 integer,          intent(out) :: ierr
 type(dump_h) :: hdr
 character(len=lenid) :: fileid
 integer :: iunit,ierr1

 time  = 0._dp
 iunit = free_unit()
 call open_dumpfile_r(iunit,fname,fileid,ierr,requiretags=.true.)
 if (ierr /= 0) return
 call read_header(iunit,hdr,ierr,tagged=.true.)
 if (ierr == 0) call extract('time',time,hdr,ierr)
 call free_header(hdr,ierr1)
 close(iunit)

end subroutine read_dump_time

integer function free_unit()
 logical :: isopen

 do free_unit = 200,999
    inquire(unit=free_unit,opened=isopen)
    if (.not.isopen) return
 enddo

end function free_unit

!-----------------------------------------------------------------------
!+
!  Called every step: spawn, move, release and remove particles
!+
!-----------------------------------------------------------------------
subroutine inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,&
                            npart,npart_old,npartoftype,dtinject)
 use part,        only:kill_particle,nptmass
 use partinject,  only:updated_particle
 use replayutils, only:schedule_next_spawn,interp_track,iend_release,iend_drop
 real,    intent(in)    :: time,dtlast
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(inout) :: npart,npart_old
 integer, intent(inout) :: npartoftype(:)
 real,    intent(out)   :: dtinject
 real(dp) :: xc(3),vc(3),tau,xr(3),vr(3),u,h,s(max(tab%nscal,1)),rot(3,3)
 integer  :: k,n,ia,i
 logical  :: found

 dtinject = dt_rep
 ! keeps update_injected_particles enforcing dtinject on every step
 updated_particle = .true.
 call get_centre(xyzmh_ptmass,vxyz_ptmass,xc,vc)
 call report_cycle(time)

 if (.not.started) then
    call start_replay(time,xc,vc,npart,npartoftype,xyzh,vxyzu)
    started = .true.
 endif
 call check_indices(npart,xyzh)

 ! spawn
 do
    call schedule_next_spawn(sch,real(time,dp),k,n,found)
    if (.not.found) exit
    if (include_drop == 0 .and. tab%iend(k) == iend_drop) cycle
    call spawn_instance(k,n,time,xc,vc,npart,npartoftype,xyzh,vxyzu)
 enddo

 ! move; release or remove at the end of the track
 ia = 1
 do while (ia <= nact)
    k   = act_k(ia)
    n   = act_n(ia)
    i   = act_i(ia)
    tau = time - t_ep - n*tab%t_rec + tab%t0
    call interp_track(tab,k,tau,xr,vr,u,h,s)
    call get_rotation(n,rot)
    xr = xc + matmul(rot,xr)
    vr = vc + matmul(rot,vr)
    if (tau >= tab%t_exit(k)) then
       if (tab%iend(k) == iend_release) then
          call release_particle(i,xr,vr,u,h,s,npart,npartoftype,xyzh,vxyzu)
          if (ideplete_sink == 1 .and. nptmass > 0) xyzmh_ptmass(4,1) = xyzmh_ptmass(4,1) - pmass
          nrel_c = nrel_c + 1
       else
          call kill_particle(i,npartoftype)
          nend_c = nend_c + 1
       endif
       call remove_instance(ia)
       cycle
    endif
    xyzh(1:3,i)  = real(xr)
    xyzh(4,i)    = real(h)
    vxyzu(1:3,i) = real(vr)
    if (size(vxyzu,1) >= 4) vxyzu(4,i) = real(u)
    call apply_scalars(i,s)
    ia = ia + 1
 enddo

 call kill_inner_gas(xc,npart,npartoftype,xyzh)
 if (ireplay_lum == 1) call replay_luminosity(time,xyzmh_ptmass)

end subroutine inject_particles

!-----------------------------------------------------------------------
!+
!  First call (fresh start or restart): the buffer is rebuilt from the
!  table, so nothing but t_epoch needs to be saved
!+
!-----------------------------------------------------------------------
subroutine start_replay(time,xc,vc,npart,npartoftype,xyzh,vxyzu)
 use part,        only:iphase,iamtype,iboundary,isdead_or_accreted,kill_particle
 use io,          only:iprint
 use replayutils, only:active_instances_at,schedule_reset,iend_drop
 real,     intent(in)    :: time
 real(dp), intent(in)    :: xc(3),vc(3)
 integer,  intent(inout) :: npart,npartoftype(:)
 real,     intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer, allocatable :: klist(:),nlist(:)
 integer :: i,j,nkill,na

 nkill = 0
 do i = 1,npart
    if (isdead_or_accreted(xyzh(4,i))) cycle
    if (iamtype(iphase(i)) == iboundary) then
       call kill_particle(i,npartoftype)
       nkill = nkill + 1
    endif
 enddo
 call active_instances_at(sch,tab,real(time,dp),na,klist,nlist)
 do j = 1,na
    if (include_drop == 0 .and. tab%iend(klist(j)) == iend_drop) cycle
    call spawn_instance(klist(j),nlist(j),time,xc,vc,npart,npartoftype,xyzh,vxyzu)
 enddo
 call schedule_reset(sch,real(time,dp))
 write(iprint,"(a,es12.4,a,i8,a,i8,a)") ' replay: start at t = ',time,': removed ',nkill, &
       ' boundary particles, spawned ',nact,' buffer particles'
 if (time < t_ep) write(iprint,"(a)") ' replay: note: the run starts before t_epoch'

end subroutine start_replay

!-----------------------------------------------------------------------
!+
!  Create the boundary particle of instance (k,n)
!+
!-----------------------------------------------------------------------
subroutine spawn_instance(k,n,time,xc,vc,npart,npartoftype,xyzh,vxyzu)
 use part,        only:iboundary,iorig,n_nucleation
 use dim,         only:do_nucleation
 use partinject,  only:add_or_update_particle
 use replayutils, only:interp_track
 integer,  intent(in)    :: k,n
 real,     intent(in)    :: time
 real(dp), intent(in)    :: xc(3),vc(3)
 integer,  intent(inout) :: npart,npartoftype(:)
 real,     intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real(dp) :: tau,xr(3),vr(3),u,h,s(max(tab%nscal,1)),rot(3,3)
 real     :: x(3),v(3)
 integer  :: ipart

 tau = time - t_ep - n*tab%t_rec + tab%t0
 call interp_track(tab,k,tau,xr,vr,u,h,s)
 call get_rotation(n,rot)
 x = real(xc + matmul(rot,xr))
 v = real(vc + matmul(rot,vr))
 ipart = npart + 1
 if (do_nucleation .and. is_nuc > 0) then
    call add_or_update_particle(iboundary,x,v,real(h),real(u),ipart,npart,npartoftype,xyzh,vxyzu, &
                                JKmuS=real(s(is_nuc:is_nuc+n_nucleation-1)))
 else
    call add_or_update_particle(iboundary,x,v,real(h),real(u),ipart,npart,npartoftype,xyzh,vxyzu)
 endif
 call apply_scalars(ipart,s)
 call add_instance(k,n,ipart,iorig(ipart))
 nspawn_c = nspawn_c + 1

end subroutine spawn_instance

!-----------------------------------------------------------------------
!+
!  The track of a buffer particle has crossed r_inj: it becomes gas
!+
!-----------------------------------------------------------------------
subroutine release_particle(i,x,v,u,h,s,npart,npartoftype,xyzh,vxyzu)
 use part,       only:igas,iphase,iunknown,n_nucleation
 use dim,        only:do_nucleation
 use partinject, only:add_or_update_particle
 integer,  intent(in)    :: i
 real(dp), intent(in)    :: x(3),v(3),u,h,s(:)
 integer,  intent(inout) :: npart,npartoftype(:)
 real,     intent(inout) :: xyzh(:,:),vxyzu(:,:)

 if (do_nucleation .and. is_nuc > 0) then
    call add_or_update_particle(igas,real(x),real(v),real(h),real(u),i,npart,npartoftype,xyzh,vxyzu, &
                                JKmuS=real(s(is_nuc:is_nuc+n_nucleation-1)))
 else
    call add_or_update_particle(igas,real(x),real(v),real(h),real(u),i,npart,npartoftype,xyzh,vxyzu)
 endif
 ! add_or_update_particle resets mu, gamma and alpha
 call apply_scalars(i,s)
 ! flag the particle so that update_injected_particles sets its timestep
 iphase(i) = int(iunknown,kind=kind(iphase))

end subroutine release_particle

!-----------------------------------------------------------------------
!+
!  Copy the recorded per-particle quantities to particle i
!+
!-----------------------------------------------------------------------
subroutine apply_scalars(i,s)
 use dim,  only:do_nucleation,update_muGamma,maxalpha,maxp,store_dust_temperature,h2chemistry
 use part, only:nucleation,n_nucleation,eos_vars,imu,igamma,alphaind,dust_temp,abundance,nabundances, &
                idmu,idgamma
 integer,  intent(in) :: i
 real(dp), intent(in) :: s(:)

 if (do_nucleation .and. is_nuc > 0) nucleation(1:n_nucleation,i) = real(s(is_nuc:is_nuc+n_nucleation-1))
 if (update_muGamma) then
    if (is_mu > 0) then
       eos_vars(imu,i) = real(s(is_mu))
    elseif (is_nuc > 0) then
       eos_vars(imu,i) = real(s(is_nuc+idmu-1))
    endif
    if (is_gam > 0) then
       eos_vars(igamma,i) = real(s(is_gam))
    elseif (is_nuc > 0) then
       eos_vars(igamma,i) = real(s(is_nuc+idgamma-1))
    endif
 endif
 if (maxalpha == maxp .and. is_alpha > 0) alphaind(1,i) = real(s(is_alpha),kind=kind(alphaind))
 if (store_dust_temperature .and. is_tdust > 0) dust_temp(i) = real(s(is_tdust))
 if (h2chemistry .and. is_abund > 0) then
    if (is_abund + nabundances - 1 <= tab%nscal) abundance(1:nabundances,i) = real(s(is_abund:is_abund+nabundances-1))
 endif

end subroutine apply_scalars

!-----------------------------------------------------------------------
!+
!  Particles are renumbered (shuffle_part) when dead particles are
!  removed before a dump: find the buffer particles again by iorig
!+
!-----------------------------------------------------------------------
subroutine check_indices(npart,xyzh)
 use part,        only:iorig,iphase,iamtype,iboundary,isdead_or_accreted
 use io,          only:warning
 use replayutils, only:sort_index_i8
 integer, intent(in) :: npart
 real,    intent(in) :: xyzh(:,:)
 integer(kind=8), allocatable :: ids(:)
 integer,         allocatable :: idx(:),ipos(:)
 integer :: ia,i,nb,lo,hi,mid,nmiss
 logical :: ok

 ok = .true.
 do ia = 1,nact
    i = act_i(ia)
    if (i < 1 .or. i > npart) then
       ok = .false.
       exit
    endif
    if (iorig(i) /= act_id(ia) .or. iamtype(iphase(i)) /= iboundary .or. isdead_or_accreted(xyzh(4,i))) then
       ok = .false.
       exit
    endif
 enddo
 if (ok) return

 nb = 0
 do i = 1,npart
    if (isdead_or_accreted(xyzh(4,i))) cycle
    if (iamtype(iphase(i)) == iboundary) nb = nb + 1
 enddo
 allocate(ids(max(nb,1)),ipos(max(nb,1)),idx(max(nb,1)))
 nb = 0
 do i = 1,npart
    if (isdead_or_accreted(xyzh(4,i))) cycle
    if (iamtype(iphase(i)) /= iboundary) cycle
    nb = nb + 1
    ids(nb)  = iorig(i)
    ipos(nb) = i
 enddo
 call sort_index_i8(nb,ids,idx)
 nmiss = 0
 ia = 1
 do while (ia <= nact)
    lo = 0
    hi = nb + 1
    do while (hi - lo > 1)
       mid = (lo + hi)/2
       if (ids(idx(mid)) <= act_id(ia)) then
          lo = mid
       else
          hi = mid
       endif
    enddo
    if (lo >= 1) then
       if (ids(idx(lo)) == act_id(ia)) then
          act_i(ia) = ipos(idx(lo))
          ia = ia + 1
          cycle
       endif
    endif
    nmiss = nmiss + 1
    call remove_instance(ia)
 enddo
 if (nmiss > 0) then
    nlost = nlost + nmiss
    call warning(label,'buffer particles lost',var='n',ival=nmiss)
 endif

end subroutine check_indices

!-----------------------------------------------------------------------
!+
!  Remove gas that has fallen back inside r_kill_inner
!+
!-----------------------------------------------------------------------
subroutine kill_inner_gas(xc,npart,npartoftype,xyzh)
 use part, only:iphase,iamtype,igas,isdead_or_accreted,kill_particle
 real(dp), intent(in)    :: xc(3)
 integer,  intent(in)    :: npart
 integer,  intent(inout) :: npartoftype(:)
 real,     intent(inout) :: xyzh(:,:)
 integer  :: i
 real(dp) :: dx(3)

 do i = 1,npart
    if (isdead_or_accreted(xyzh(4,i))) cycle
    if (iamtype(iphase(i)) /= igas) cycle
    dx = xyzh(1:3,i) - xc
    if (dot_product(dx,dx) < r_kill2) then
       call kill_particle(i,npartoftype)
       nkill_c = nkill_c + 1
    endif
 enddo

end subroutine kill_inner_gas

!-----------------------------------------------------------------------
!+
!  Recorded luminosity at the same phase of the window. If the sinks
!  were merged the luminosities of the inner sinks are summed.
!+
!-----------------------------------------------------------------------
subroutine replay_luminosity(time,xyzmh_ptmass)
 use part,        only:nptmass,ilum
 use replayutils, only:isk_lum
 real, intent(in)    :: time
 real, intent(inout) :: xyzmh_ptmass(:,:)
 real(dp) :: tph,f
 integer  :: lo,hi,mid,k

 if (tab%nsinkrec < 2 .or. nptmass < 1) return
 tph = tab%t0 + modulo(time - t_ep,tab%t_rec)
 if (tph <= tab%tsink(1)) then
    lo = 1
    hi = 1
 elseif (tph >= tab%tsink(tab%nsinkrec)) then
    lo = tab%nsinkrec
    hi = lo
 else
    lo = 1
    hi = tab%nsinkrec
    do while (hi - lo > 1)
       mid = (lo + hi)/2
       if (tab%tsink(mid) <= tph) then
          lo = mid
       else
          hi = mid
       endif
    enddo
 endif
 f = 0._dp
 if (hi > lo) f = (tph - tab%tsink(lo))/(tab%tsink(hi) - tab%tsink(lo))
 if (nptmass == tab%nsink) then
    do k = 1,nptmass
       xyzmh_ptmass(ilum,k) = real((1._dp - f)*tab%sink(isk_lum,k,lo) + f*tab%sink(isk_lum,k,hi))
    enddo
 else
    xyzmh_ptmass(ilum,1) = real(sum((1._dp - f)*tab%sink(isk_lum,:,lo) + f*tab%sink(isk_lum,:,hi)))
 endif

end subroutine replay_luminosity

!-----------------------------------------------------------------------
!+
!  Centre of mass of the sinks (the recorded positions are relative to it)
!+
!-----------------------------------------------------------------------
subroutine get_centre(xyzmh_ptmass,vxyz_ptmass,xc,vc)
 use part, only:nptmass
 real,     intent(in)  :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 real(dp), intent(out) :: xc(3),vc(3)
 real(dp) :: mtot
 integer  :: k

 xc   = 0._dp
 vc   = 0._dp
 mtot = 0._dp
 do k = 1,nptmass
    if (xyzmh_ptmass(4,k) <= 0.) cycle
    xc   = xc + xyzmh_ptmass(4,k)*xyzmh_ptmass(1:3,k)
    vc   = vc + xyzmh_ptmass(4,k)*vxyz_ptmass(1:3,k)
    mtot = mtot + xyzmh_ptmass(4,k)
 enddo
 if (mtot > 0._dp) then
    xc = xc/mtot
    vc = vc/mtot
 endif

end subroutine get_centre

!-----------------------------------------------------------------------
!+
!  Rotation of cycle n (cached)
!+
!-----------------------------------------------------------------------
subroutine get_rotation(n,rot)
 use replayutils, only:cycle_rotation
 integer,  intent(in)  :: n
 real(dp), intent(out) :: rot(3,3)
 integer :: j

 do j = 1,ncache
    if (cache_n(j) == n) then
       rot = cache_rot(:,:,j)
       return
    endif
 enddo
 call cycle_rotation(n,irotate,iseed_rotate,rot)
 j = modulo(n,ncache) + 1
 cache_n(j) = n
 cache_rot(:,:,j) = rot

end subroutine get_rotation

!-----------------------------------------------------------------------
!+
!  List of active instances
!+
!-----------------------------------------------------------------------
subroutine add_instance(k,n,i,id)
 integer,         intent(in) :: k,n,i
 integer(kind=8), intent(in) :: id
 integer,         allocatable :: itmp(:)
 integer(kind=8), allocatable :: i8tmp(:)
 integer :: ncap

 if (.not.allocated(act_k)) allocate(act_k(4096),act_n(4096),act_i(4096),act_id(4096))
 if (nact >= size(act_k)) then
    ncap = 2*size(act_k)
    allocate(itmp(ncap)); itmp(1:nact) = act_k(1:nact); call move_alloc(itmp,act_k)
    allocate(itmp(ncap)); itmp(1:nact) = act_n(1:nact); call move_alloc(itmp,act_n)
    allocate(itmp(ncap)); itmp(1:nact) = act_i(1:nact); call move_alloc(itmp,act_i)
    allocate(i8tmp(ncap)); i8tmp(1:nact) = act_id(1:nact); call move_alloc(i8tmp,act_id)
 endif
 nact = nact + 1
 act_k(nact)  = k
 act_n(nact)  = n
 act_i(nact)  = i
 act_id(nact) = id

end subroutine add_instance

subroutine remove_instance(ia)
 integer, intent(in) :: ia

 act_k(ia)  = act_k(nact)
 act_n(ia)  = act_n(nact)
 act_i(ia)  = act_i(nact)
 act_id(ia) = act_id(nact)
 nact = nact - 1

end subroutine remove_instance

!-----------------------------------------------------------------------
!+
!  One line per replay cycle
!+
!-----------------------------------------------------------------------
subroutine report_cycle(time)
 use io,      only:iprint
 use units,   only:umass,utime
 use physcon, only:solarm,years
 real, intent(in) :: time
 integer :: m

 m = floor((time - t_ep)/tab%t_rec)
 if (m == mcycle_last) return
 if (iverbose_replay > 0 .and. mcycle_last > -huge(1)) then
    write(iprint,"(a,i6,a,i8,a,i8,a,i8,a,i8,a,es10.3,a)") ' replay: cycle ',mcycle_last,': spawned ',nspawn_c, &
          ', released ',nrel_c,', ended in buffer ',nend_c,', removed inside ',nkill_c,', Mdot = ', &
          nrel_c*pmass/tab%t_rec*(umass/solarm)/(utime/years),' Msun/yr'
 endif
 mcycle_last = m
 nrel_c   = 0
 nend_c   = 0
 nkill_c  = 0
 nspawn_c = 0

end subroutine report_cycle

!-----------------------------------------------------------------------
!+
!  Nothing to save between runs: the buffer is rebuilt from the table
!+
!-----------------------------------------------------------------------
subroutine update_injected_par
end subroutine update_injected_par

!-----------------------------------------------------------------------
!+
!  Write / read options
!+
!-----------------------------------------------------------------------
subroutine write_options_inject(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# options controlling replay injection'
 call write_inopt(trim(replay_dumps),'replay_dumps','first full dump of the inner run to build the table from',iunit)
 call write_inopt(replay_period,'replay_period','replay period = pulsation or orbital period (days)',iunit)
 call write_inopt(r_inj_au,'r_inj','interface radius: particles are released here (au)',iunit)
 call write_inopt(buffer_width_au,'buffer_width','width of the replay buffer r_inj - r_b (au, 0=auto from h)',iunit)
 call write_inopt(trim(replay_file),'replay_file','track table (built from the dumps if absent)',iunit)
 call write_inopt(t_epoch,'t_epoch','time at which replay cycle 0 starts (code units, <0: time of the first dump)',iunit)
 call write_inopt(r_kill_inner_au,'r_kill_inner','gas inside this radius is removed (au, <=0: r_b)',iunit)
 call write_inopt(dt_replay,'dt_replay','maximum timestep while replaying (periods)',iunit)
 call write_inopt(irotate,'irotate','per-cycle transformation (0=none, 1=random rotation, 2=z-mirror)',iunit)
 call write_inopt(ideplete_sink,'ideplete_sink','remove released mass from sink 1 (0=no, 1=yes)',iunit)
 call write_inopt(ireplay_lum,'ireplay_lum','replay the recorded sink luminosity (0=no, 1=yes)',iunit)
 call write_inopt(include_drop,'include_drop','replay tracks that leave the buffer inwards (0=no, 1=yes)',iunit)
 call write_inopt(iverbose_replay,'verbose_replay','print a summary every cycle (0=no, 1=yes)',iunit)

end subroutine write_options_inject

subroutine read_options_inject(db,nerr)
 use infile_utils, only:inopts,read_inopt
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr
 integer :: ierr

 ! all options have defaults, so that moddump_replay can read the .in file of the inner run
 call read_inopt(replay_dumps,'replay_dumps',db,err=ierr)
 if (ierr /= 0 .or. len_trim(replay_dumps) == 0) replay_dumps = 'wind_00000'
 call read_inopt(replay_period,'replay_period',db,errcount=nerr,min=0.,default=300.)
 call read_inopt(r_inj_au,'r_inj',db,errcount=nerr,min=0.,default=50.)
 call read_inopt(buffer_width_au,'buffer_width',db,errcount=nerr,min=0.,default=0.)
 call read_inopt(replay_file,'replay_file',db,err=ierr)
 if (ierr /= 0 .or. len_trim(replay_file) == 0) replay_file = 'replay_tracks.dat'
 call read_inopt(t_epoch,'t_epoch',db,errcount=nerr,default=-1.)
 call read_inopt(r_kill_inner_au,'r_kill_inner',db,errcount=nerr,default=0.)
 call read_inopt(dt_replay,'dt_replay',db,errcount=nerr,min=0.,default=0.02)
 call read_inopt(irotate,'irotate',db,errcount=nerr,min=0,max=2,default=0)
 call read_inopt(ideplete_sink,'ideplete_sink',db,errcount=nerr,min=0,max=1,default=1)
 call read_inopt(ireplay_lum,'ireplay_lum',db,errcount=nerr,min=0,max=1,default=0)
 call read_inopt(include_drop,'include_drop',db,errcount=nerr,min=0,max=1,default=1)
 call read_inopt(iverbose_replay,'verbose_replay',db,errcount=nerr,min=0,max=1,default=1)

end subroutine read_options_inject

end module inject
