!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module replayutils
!
! Track table of the periodic track replay (see inject_replay), and the
! routines that replay it in the outer run.
!
! :References: None
!
! :Owner: Owen Vermeulen
!
! :Runtime parameters: None
!
! :Dependencies: dim, dump_utils, fileutils, infile_utils, io, kernel,
!   part, partinject, physcon, readwrite_dumps_common, sortutils, timing,
!   units
!
 use io, only:fatal,warning
 implicit none

 ! used by inject_replay
 public :: init_replay, start_replay, buffer_is_intact, rebuild_buffer, &
           replay_tracks, kill_inner_gas, replay_luminosity, get_centre_of_mass
 public :: dt_replay, window_start, replay_period, replay_started, particle_mass, n_dumps, n_sinks, dump_time
 ! used by the setup of the outer run
 public :: inner_dump, read_inner_dump, inner_infile, merge_sinks
 private

 integer,          parameter :: nbase         = 8   ! particle state: x(3), v(3), u, h, then the scalars
 integer,          parameter :: lenlabel      = 16
 integer,          parameter :: iend_release  = 1   ! track leaves through r_inj: becomes gas
 integer,          parameter :: iend_drop     = 2   ! track leaves through r_b: removed
 integer,          parameter :: iend_vanish   = 3   ! track still in the buffer at the end of the window: removed
 integer,          parameter :: iend_open     = 0
 integer,          parameter :: iend_discard  = -1
 integer,          parameter :: table_version = 2
 character(len=*), parameter :: table_magic   = 'PHANTOMREPLAYTAB'
 character(len=*), parameter :: table_file    = 'replay_tracks.dat'
 character(len=*), parameter :: label         = 'utils_replay'

 ! one full dump of the inner run
 type inner_dump
    real    :: time      = 0.
    real    :: udist     = 1.
    real    :: umass     = 1.
    real    :: utime     = 1.
    real    :: pmass     = 0.
    real    :: gamma     = 0.
    real    :: polyk     = 0.
    integer :: npart     = 0
    integer :: nptmass   = 0
    logical :: got_iorig = .false.
    real,            allocatable :: xyz(:,:), vxyz(:,:), u(:), h(:), scalars(:,:)
    integer(kind=8), allocatable :: iorig(:)
    integer(kind=1), allocatable :: itype(:)
    logical,         allocatable :: got_scalar(:)
    real,            allocatable :: xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 end type inner_dump

 ! the gas of one dump in the band around the buffer, in order of id
 type band
    real    :: time = 0.
    integer :: n    = 0
    integer(kind=8), allocatable :: id(:)
    real(kind=4),    allocatable :: state(:,:)
    integer,         allocatable :: track(:)   ! track of each particle (0 if none)
 end type band

 ! track table
 character(len=200) :: table_dumpfile = ''
 real    :: table_udist   = 1.
 real    :: table_umass   = 1.
 real    :: table_utime   = 1.
 real    :: particle_mass = 0.
 real    :: r_inj         = 0.
 real    :: r_b           = 0.
 real    :: window_start  = 0.
 real    :: replay_period = 0.
 integer :: n_scalars     = 0
 integer :: n_tracks      = 0
 integer :: n_samples     = 0
 integer :: n_dumps       = 0
 integer :: n_sinks       = 0
 character(len=lenlabel), allocatable :: scalar_names(:)
 real,         allocatable :: track_t_entry(:), track_t_exit(:), sample_time(:)
 integer,      allocatable :: track_end(:), track_first(:), track_last(:)
 real(kind=4), allocatable :: sample_state(:,:)
 real,         allocatable :: dump_time(:), dump_lum(:,:)
 integer,      allocatable :: sample_track(:)   ! track of each sample, until the table is finished

 ! checks of the inner dumps, printed when the table is built
 integer :: n_lost = 0, n_out = 0, n_in = 0, n_moved = 0
 real    :: sum_moved = 0.

 ! replaying the track table
 real    :: dt_replay      = 0.
 integer :: cycle_last     = 0
 logical :: replay_started = .false.
 ! the particle on each track (0 if none), the last cycle in which the track started, and the iorig of the particle
 integer,         allocatable :: track_particle(:), track_cycle(:)
 integer(kind=8), allocatable :: track_particle_id(:)
 ! rows of the recorded quantities in the particle state (0 if not recorded)
 integer :: irow_mu = 0, irow_gamma = 0, irow_alpha = 0, irow_tdust = 0, irow_abundance = 0
 ! counts of the current cycle
 integer :: n_spawned = 0, n_released = 0, n_removed = 0, n_killed = 0

contains

!----------------------------------------------------------------
!+
!  Read the track table, or build it if it was made from other
!  dumps or for another r_inj, and prepare the replay
!+
!----------------------------------------------------------------
subroutine init_replay(dumpfile,r_inj_wanted)
 use physcon, only:au
 use units,   only:udist,umass,utime
 use part,    only:abundance_label,nabundances
 character(len=*), intent(in) :: dumpfile
 real,             intent(in) :: r_inj_wanted
 real, allocatable :: h_over_v(:)
 real    :: vmag
 integer :: i,n,ierr

 call read_track_table(ierr)
 if (ierr == 0 .and. (table_dumpfile /= dumpfile .or. abs(r_inj/r_inj_wanted - 1.) > 1.e-6)) then
    print *, ''
    print *, table_file//' was built from other dumps or for another r_inj'
    ierr = 1
 endif
 if (ierr /= 0) then
    call build_track_table(dumpfile,r_inj_wanted)
    call write_track_table(ierr)
    if (ierr /= 0) call warning(label,'could not write '//table_file)
    if (ierr == 0) print *, 'Track table written to '//table_file
 else
    print *, ''
    print *, 'Track table read from '//table_file
 endif
 if (n_tracks < 1) call fatal(label,'no tracks in the table')
 if (abs(table_udist/udist - 1.) > 1.e-6 .or. abs(table_umass/umass - 1.) > 1.e-6 .or. &
     abs(table_utime/utime - 1.) > 1.e-6) call fatal(label,'units differ from those of the inner run')

 ! time between two calls: at most P/50, and the time to move 0.1 h in the buffer
 allocate(h_over_v(n_samples))
 n = 0
 do i = 1,n_samples
    vmag = norm2(real(sample_state(4:6,i)))
    if (vmag > 0. .and. sample_state(8,i) > 0.) then
       n = n + 1
       h_over_v(n) = sample_state(8,i)/vmag
    endif
 enddo
 dt_replay = 0.02*replay_period
 if (n > 0) dt_replay = min(dt_replay,0.1*percentile(h_over_v(1:n),0.05))

 ! rows of the recorded quantities in the particle state
 irow_mu        = find_scalar_row('eos:mu',1)
 irow_gamma     = find_scalar_row('eos:gamma',1)
 irow_alpha     = find_scalar_row('alpha',1)
 irow_tdust     = find_scalar_row('Tdust',1)
 irow_abundance = find_scalar_row('abund:'//trim(abundance_label(1)),nabundances)

 if (allocated(track_particle)) deallocate(track_particle,track_cycle,track_particle_id)
 allocate(track_particle(n_tracks),track_cycle(n_tracks),track_particle_id(n_tracks))
 track_particle = 0
 replay_started = .false.

 print *, ''
 print *, ' Tracks                     :', n_tracks
 print *, ' Samples                    :', n_samples
 print *, ' r_b (au)                   :', r_b*udist/au
 print *, ' r_inj (au)                 :', r_inj*udist/au
 print *, ' Start of the window        :', window_start
 print *, ' Replay period              :', replay_period
 print *, ' dtinject                   :', dt_replay
 print *, ''

end subroutine init_replay

!----------------------------------------------------------------
!+
!  Row in the particle state of the n quantities starting with
!  name (0 if they are not in the table)
!+
!----------------------------------------------------------------
integer function find_scalar_row(name,n)
 character(len=*), intent(in) :: name
 integer,          intent(in) :: n
 integer :: i

 find_scalar_row = 0
 do i = 1,n_scalars
    if (trim(scalar_names(i)) == trim(name)) then
       if (i + n - 1 <= n_scalars) find_scalar_row = nbase + i
       return
    endif
 enddo

end function find_scalar_row

!----------------------------------------------------------------
!+
!  Build the track table from one period of dumps of the inner run
!+
!----------------------------------------------------------------
subroutine build_track_table(dumpfile,r_inj_new)
 use physcon,   only:au,solarm,years
 use units,     only:udist,umass,utime
 use kernel,    only:radkern
 use part,      only:igas,iboundary,ilum
 use fileutils, only:getnextfilename,basename
 use timing,    only:getused
 use io,        only:stdout
 character(len=*), intent(in) :: dumpfile
 real,             intent(in) :: r_inj_new
 type(inner_dump) :: dump
 type(band)       :: prev,now
 character(len=200) :: filename
 character(len=lenlabel), allocatable :: all_names(:)
 integer, allocatable :: icolumn(:)
 real,    allocatable :: h_near(:)
 real    :: xc(3),vc(3),r,r_boundary,buffer_width,umdot
 real(kind=4) :: t_start,t_now   ! wall time (s), from getused
 integer :: i,n,n_dumps_read,n_release
 logical :: iexist

 call getused(t_start)
 filename = dumpfile
 print *, ''
 print *, 'Building the track table from '//trim(filename)//' and the following dumps'
 call read_dump_or_die(filename,dump)
 if (.not. dump%got_iorig) call fatal(label,trim(filename)//' has no particle ids (iorig)')
 if (abs(dump%udist/udist - 1.) > 1.e-6 .or. abs(dump%umass/umass - 1.) > 1.e-6 .or. &
     abs(dump%utime/utime - 1.) > 1.e-6) call fatal(label,'units of '//trim(filename)//' differ from those of this run')

 call clear_table()
 n_lost    = 0
 n_out     = 0
 n_in      = 0
 n_moved   = 0
 sum_moved = 0.
 replay_period = calculate_replay_period(dumpfile,dump)

 ! buffer width: 1.2 kernel radii of the 90th percentile of h at r_inj, and the outermost
 ! boundary particle of the inner run (for a nested replay: the r_inj of the inner run)
 call get_centre_of_mass(dump%nptmass,dump%xyzmh_ptmass,dump%vxyz_ptmass,xc,vc)
 allocate(h_near(dump%npart))
 n = 0
 r_boundary = 0.
 do i = 1,dump%npart
    r = norm2(dump%xyz(:,i) - xc)
    if (dump%itype(i) == igas .and. abs(r - r_inj_new) < 0.05*r_inj_new) then
       n = n + 1
       h_near(n) = dump%h(i)
    elseif (dump%itype(i) == iboundary) then
       r_boundary = max(r_boundary,r)
    endif
 enddo
 if (n < 20) call fatal(label,'too few particles at r_inj in '//trim(filename)//' to set the buffer width')
 buffer_width = 1.2*radkern*percentile(h_near(1:n),0.9)
 if (buffer_width >= r_inj_new) call fatal(label,'buffer wider than r_inj',var='buffer width (au)',val=buffer_width*udist/au)
 ! the band around the buffer must be gas in the inner run, outside its boundary particles
 if (r_inj_new - 1.5*buffer_width <= r_boundary) call fatal(label,'the band around the buffer reaches the '// &
    'boundary particles of the inner run: increase r_inj',var='r_inj (au) at least',val=(r_boundary + 1.5*buffer_width)*udist/au)

 ! what the table is built from, and the per-particle quantities found in the dump
 table_dumpfile = dumpfile
 table_udist    = udist
 table_umass    = umass
 table_utime    = utime
 particle_mass  = dump%pmass
 r_inj          = r_inj_new
 r_b            = r_inj_new - buffer_width
 window_start   = dump%time
 call get_scalar_names(all_names)
 n_scalars = count(dump%got_scalar)
 allocate(icolumn(n_scalars),scalar_names(n_scalars))
 n = 0
 do i = 1,size(all_names)
    if (dump%got_scalar(i)) then
       n = n + 1
       icolumn(n)      = i
       scalar_names(n) = all_names(i)
    endif
 enddo

 print *, ''
 print *, ' Window (code time)         :', window_start, ' to', window_start + replay_period
 print *, ' Buffer, r_b to r_inj (au)  :', r_b*udist/au, ' to', r_inj*udist/au
 print "(1x,a,99(1x,a))", ' Recorded per particle      : x v u h', (trim(scalar_names(i)), i=1,n_scalars)
 print *, ''
 print "(1x,a)", '        time % of P   in band  in buffer  released    samples  wall (s)   dump'

 ! add the dumps until one period has been covered, with one line of progress per dump
 n_dumps_read = 0
 do
    call get_band(dump,icolumn,r_inj_new - 1.5*buffer_width,r_inj_new + 0.5*buffer_width,now)
    call record_luminosity(dump%time,dump%xyzmh_ptmass(ilum,1:dump%nptmass))
    call add_dump_to_table(prev,now)
    prev = now
    n_dumps_read = n_dumps_read + 1
    call getused(t_now)
    print "(1x,f12.4,i6,'%',i10,i11,i10,i11,f10.2,3x,a)", dump%time, &
          nint(100.*(dump%time - window_start)/replay_period), now%n, count(now%track > 0), &
          count(track_end(1:n_tracks) == iend_release), n_samples, t_now - t_start, trim(basename(filename))
    flush(stdout)
    if (dump%time >= window_start + replay_period*(1. - 1.e-6)) exit
    filename = getnextfilename(filename)
    inquire(file=trim(filename),exist=iexist)
    if (.not. iexist) call fatal(label,'the dumps end before one replay period: '//trim(filename)//' is missing')
    call read_dump_or_die(filename,dump)
    if (dump%time <= prev%time) call fatal(label,trim(filename)//' is not later than the previous dump')
 enddo
 call finish_track_table(prev)

 umdot     = (umass/solarm)/(utime/years)
 n_release = count(track_end(1:n_tracks) == iend_release)
 print *, ''
 print *, ' Dumps read                 :', n_dumps_read
 print *, ' Release tracks             :', n_release
 print *, ' Drop tracks                :', count(track_end(1:n_tracks) == iend_drop)
 print *, ' In the buffer at the start :', count(track_t_entry(1:n_tracks) <= window_start)
 print *, ' In the buffer at the end   :', count(track_end(1:n_tracks) == iend_vanish)
 print *, ' Replayed Mdot (Msun/yr)    :', n_release*particle_mass/replay_period*umdot
 if (n_out > 0) print *, ' Crossings of r_inj in/out  :', real(n_in)/n_out
 call getused(t_now)
 print "(1x,a,f10.2)", ' Time to build (s)          :', t_now - t_start
 print *, ''
 if (n_in > 0.05*n_out) call warning(label,'more than 5% of the gas crossing r_inj falls back: '// &
                                           'r_inj is too close to the star')
 if (n_moved > 0) then
    if (sum_moved/n_moved > 0.5) call warning(label,'particles move more than h/2 between dumps: '// &
                                                    'write dumps more often')
 endif
 if (n_lost > 0) call warning(label,'particles left the band between two dumps without crossing r_inj or r_b', &
                              var='n',ival=n_lost)

end subroutine build_track_table

!----------------------------------------------------------------
!+
!  Replay period: that of the inner run if it is itself a replay
!  run (nested replay), otherwise the orbital period of a binary,
!  or the pulsation period in the .in file of the inner run
!+
!----------------------------------------------------------------
real function calculate_replay_period(dumpfile,dump) result(period)
 use physcon,      only:au,days,pi
 use units,        only:udist,utime,get_G_code
 use infile_utils, only:open_db_from_file,close_db,inopts,read_inopt
 use io,           only:ireadin
 character(len=*), intent(in) :: dumpfile
 type(inner_dump), intent(in) :: dump
 type(inopts), allocatable :: db(:)
 character(len=200) :: infile
 real    :: pulsation_period_days,pulsation_period,inner_replay_period,dx(3),dv(3),gm,one_over_a
 integer :: ierr

 ! in the .in file of the inner run (0 if not found): the pulsation period (with iwind=2 the
 ! inner run computes it, and writes it to its .in file at every dump), and the replay period
 ! (code units) that a replay run writes to its .in file
 period                = 0.
 pulsation_period_days = 0.
 inner_replay_period   = 0.
 infile                = inner_infile(dumpfile)
 call open_db_from_file(db,infile,ireadin,ierr)
 if (ierr == 0) then
    call read_inopt(pulsation_period_days,'pulsation_period',db,err=ierr)
    if (ierr /= 0) pulsation_period_days = 0.
    call read_inopt(inner_replay_period,'replay_period',db,err=ierr)
    if (ierr /= 0) inner_replay_period = 0.
 endif
 call close_db(db)
 pulsation_period = pulsation_period_days*days/utime

 ! nested replay: the flow of an inner replay run repeats with its replay period
 if (inner_replay_period > 0.) then
    period = inner_replay_period
    print *, 'Nested replay: replay period = that of the inner run (days):',period*utime/days
    return
 endif

 select case(dump%nptmass)
 case(1)
    if (pulsation_period <= 0.) call fatal(label,'single star, but no pulsation_period in '//trim(infile))
    period = pulsation_period
    print *, 'Single star: replay period = pulsation period (days):',pulsation_period_days
 case(2)
    ! semi-major axis from the energy of the relative orbit (vis-viva)
    dx         = dump%xyzmh_ptmass(1:3,2) - dump%xyzmh_ptmass(1:3,1)
    dv         = dump%vxyz_ptmass(1:3,2) - dump%vxyz_ptmass(1:3,1)
    gm         = get_G_code()*(dump%xyzmh_ptmass(4,1) + dump%xyzmh_ptmass(4,2))
    one_over_a = 2./norm2(dx) - dot_product(dv,dv)/gm
    if (one_over_a <= 0.) call fatal(label,'the two sinks of '//trim(dumpfile)//' are not bound')
    period = 2.*pi*sqrt(1./(one_over_a**3*gm))
    print *, 'Binary: replay period = orbital period (days):',period*utime/days
    print *, 'Semi-major axis (au)                         :',udist/(one_over_a*au)
    if (pulsation_period > 0. .and. abs(period/pulsation_period - nint(period/pulsation_period)) > 0.01) &
       call warning(label,'the pulsation is not periodic with the orbit',var='P_orb/P_puls',val=period/pulsation_period)
 case default
    call fatal(label,'the replay period is only known for one or two sinks',var='nptmass',ival=dump%nptmass)
 end select

end function calculate_replay_period

!----------------------------------------------------------------
!+
!  The .in file of the run that wrote a dump (wind_00440 -> wind.in)
!+
!----------------------------------------------------------------
function inner_infile(dumpfile) result(infile)
 character(len=*), intent(in) :: dumpfile
 character(len=200) :: infile

 infile = dumpfile(1:max(index(dumpfile,'_',back=.true.)-1,1))//'.in'

end function inner_infile

!----------------------------------------------------------------
!+
!  The gas of a dump in the band r_min <= r <= r_max, relative to
!  the centre of mass of the sinks, in order of id
!+
!----------------------------------------------------------------
subroutine get_band(dump,icolumn,r_min,r_max,now)
 use part,      only:igas
 use sortutils, only:indexx
 type(inner_dump), intent(in)  :: dump
 integer,          intent(in)  :: icolumn(:)
 real,             intent(in)  :: r_min,r_max
 type(band),       intent(out) :: now
 integer, allocatable :: iband(:),isort(:)
 real    :: xc(3),vc(3),r
 integer :: i,j,n

 call get_centre_of_mass(dump%nptmass,dump%xyzmh_ptmass,dump%vxyz_ptmass,xc,vc)
 allocate(iband(dump%npart))
 n = 0
 do i = 1,dump%npart
    r = norm2(dump%xyz(:,i) - xc)
    if (dump%itype(i) == igas .and. r >= r_min .and. r <= r_max) then
       n = n + 1
       iband(n) = i
    endif
 enddo
 allocate(isort(n))
 call indexx(n,dump%iorig(iband(1:n)),isort)

 now%time = dump%time
 now%n    = n
 allocate(now%id(n),now%state(nbase+size(icolumn),n),now%track(n))
 do j = 1,n
    i = iband(isort(j))
    now%id(j)             = dump%iorig(i)
    now%state(1:3,j)      = real(dump%xyz(:,i) - xc,kind=4)
    now%state(4:6,j)      = real(dump%vxyz(:,i) - vc,kind=4)
    now%state(7,j)        = real(dump%u(i),kind=4)
    now%state(8,j)        = real(dump%h(i),kind=4)
    now%state(nbase+1:,j) = real(dump%scalars(icolumn,i),kind=4)
 enddo
 now%track = 0

end subroutine get_band

!----------------------------------------------------------------
!+
!  Names of the per-particle quantities copied from the inner run,
!  in the order of the rows of dump%scalars
!+
!----------------------------------------------------------------
subroutine get_scalar_names(names)
 use part, only:nabundances,abundance_label
 character(len=lenlabel), allocatable, intent(out) :: names(:)
 integer :: i

 allocate(names(4+nabundances))
 names(1) = 'eos:mu'
 names(2) = 'eos:gamma'
 names(3) = 'alpha'
 names(4) = 'Tdust'
 do i = 1,nabundances
    names(4+i) = 'abund:'//trim(abundance_label(i))
 enddo

end subroutine get_scalar_names

!----------------------------------------------------------------
!+
!  Read what is needed from a full dump of the inner run
!+
!----------------------------------------------------------------
subroutine read_inner_dump(filename,dump,ierr)
 use dump_utils, only:open_dumpfile_r,read_header,free_header,extract,read_global_block_header, &
                      read_block_header,read_array,dump_h,lentag,lenid,ndatatypes
 use readwrite_dumps_common, only:get_options_from_fileid
 use io,         only:idisk1
 use part,       only:maxtypes,igas,nsinkproperties,xyzmh_ptmass_label,vxyz_ptmass_label, &
                      nabundances,abundance_label
 character(len=*), intent(in)  :: filename
 type(inner_dump), intent(out) :: dump   ! intent(out): emptied on entry
 integer,          intent(out) :: ierr
 integer, parameter :: maxarraylengths = 12
 type(dump_h) :: hdr
 character(len=lenid)  :: fileid
 character(len=lentag) :: tag
 integer(kind=8) :: number8(maxarraylengths)
 integer :: nums(ndatatypes,maxarraylengths)
 integer :: nblocks,narraylengths,iblock,j,ik,k,np,ns,i1,i2,is1,is2,ierr1
 real    :: massoftype(maxtypes),rk2
 logical :: tagged,phantomdump,smalldump,use_dust,matched
 logical :: gotx(3),gotv(3),gotu,goth,gott,gots(nsinkproperties),gotvs(3)

 call open_dumpfile_r(idisk1,filename,fileid,ierr,requiretags=.true.)
 if (ierr /= 0) return
 call get_options_from_fileid(fileid,tagged,phantomdump,smalldump,use_dust,ierr1)
 if (smalldump .or. .not. tagged) ierr = 1
 if (ierr == 0) call read_header(idisk1,hdr,ierr,tagged=.true.)
 if (ierr /= 0) then
    close(idisk1)
    return
 endif
 massoftype = 0.
 call extract('time',dump%time,hdr,ierr1)
 call extract('udist',dump%udist,hdr,ierr1)
 call extract('umass',dump%umass,hdr,ierr1)
 call extract('utime',dump%utime,hdr,ierr1)
 call extract('massoftype',massoftype,hdr,ierr1)
 call extract('gamma',dump%gamma,hdr,ierr1)
 call extract('RK2',rk2,hdr,ierr1)
 call extract('nparttot',np,hdr,ierr1)
 call extract('nptmass',ns,hdr,ierr1,default=0)
 dump%pmass = massoftype(igas)
 dump%polyk = rk2/1.5
 call read_global_block_header(nblocks,narraylengths,hdr,idisk1,ierr)
 call free_header(hdr,ierr1)
 if (ierr /= 0) then
    close(idisk1)
    return
 endif

 allocate(dump%xyz(3,np),dump%vxyz(3,np),dump%u(np),dump%h(np),dump%iorig(np),dump%itype(np))
 allocate(dump%scalars(4+nabundances,np),dump%got_scalar(4+nabundances))
 allocate(dump%xyzmh_ptmass(nsinkproperties,ns),dump%vxyz_ptmass(3,ns))
 dump%u            = 0.
 dump%itype        = int(igas,kind=1)
 dump%got_scalar   = .false.
 dump%xyzmh_ptmass = 0.
 dump%vxyz_ptmass  = 0.

 i2  = 0
 is2 = 0
 do iblock = 1,nblocks
    call read_block_header(narraylengths,number8,nums,idisk1,ierr)
    i1  = i2 + 1
    i2  = i2 + int(number8(1))
    is1 = is2 + 1
    if (narraylengths >= 2) is2 = is2 + int(number8(2))
    do j = 1,narraylengths
       do ik = 1,ndatatypes
          do k = 1,nums(ik,j)
             if (ierr /= 0) exit
             read(idisk1,iostat=ierr) tag
             matched = .false.
             if (j == 1) then
                ! particle arrays
                call read_array(dump%xyz,['x','y','z'],gotx,ik,i1,i2,0,idisk1,tag,matched,ierr)
                call read_array(dump%vxyz,['vx','vy','vz'],gotv,ik,i1,i2,0,idisk1,tag,matched,ierr)
                call read_array(dump%u,'u',gotu,ik,i1,i2,0,idisk1,tag,matched,ierr)
                call read_array(dump%h,'h',goth,ik,i1,i2,0,idisk1,tag,matched,ierr)
                call read_array(dump%iorig,'iorig',dump%got_iorig,ik,i1,i2,0,idisk1,tag,matched,ierr)
                call read_array(dump%itype,'itype',gott,ik,i1,i2,0,idisk1,tag,matched,ierr)
                call read_array(dump%scalars(1:4,:),['mu   ','gamma','alpha','Tdust'], &
                                dump%got_scalar(1:4),ik,i1,i2,0,idisk1,tag,matched,ierr)
                call read_array(dump%scalars(5:,:),abundance_label,dump%got_scalar(5:), &
                                ik,i1,i2,0,idisk1,tag,matched,ierr)
             elseif (j == 2) then
                ! sink arrays
                call read_array(dump%xyzmh_ptmass,xyzmh_ptmass_label,gots,ik,is1,is2,0,idisk1,tag,matched,ierr)
                call read_array(dump%vxyz_ptmass,vxyz_ptmass_label,gotvs,ik,is1,is2,0,idisk1,tag,matched,ierr)
             endif
             if (.not. matched) read(idisk1,iostat=ierr)
          enddo
       enddo
    enddo
 enddo
 close(idisk1)
 dump%npart   = i2
 dump%nptmass = is2

end subroutine read_inner_dump

!----------------------------------------------------------------
!+
!  Read a dump of the inner run, and stop if that fails
!+
!----------------------------------------------------------------
subroutine read_dump_or_die(filename,dump)
 character(len=*), intent(in)  :: filename
 type(inner_dump), intent(out) :: dump
 integer :: ierr

 call read_inner_dump(filename,dump,ierr)
 if (ierr /= 0) call fatal(label,'cannot read '//trim(filename))

end subroutine read_dump_or_die

!----------------------------------------------------------------
!+
!  Add the band particles of a dump to the table: they are matched
!  with those of the previous dump by their id
!+
!----------------------------------------------------------------
subroutine add_dump_to_table(prev,now)
 type(band), intent(in)    :: prev
 type(band), intent(inout) :: now
 integer(kind=8) :: id_prev,id_now
 integer :: i,j,k
 real    :: r

 ! walk through the particles of both dumps in order of id
 i = 1
 j = 1
 do while (i <= prev%n .or. j <= now%n)
    id_prev = huge(id_prev)
    id_now  = huge(id_now)
    if (i <= prev%n) id_prev = prev%id(i)
    if (j <= now%n)  id_now  = now%id(j)
    if (id_prev < id_now) then
       ! left the band without crossing r_b or r_inj: moved too far between dumps
       if (prev%track(i) > 0) then
          track_end(prev%track(i)) = iend_discard
          n_lost = n_lost + 1
       endif
       i = i + 1
    elseif (id_now < id_prev) then
       ! new in the band: a track starts if it is already in the buffer (first dump)
       r = norm2(real(now%state(1:3,j)))
       if (r > r_b .and. r < r_inj) then
          now%track(j) = new_track(now%time)
          call add_sample(now%track(j),now%time,now%state(:,j))
       endif
       j = j + 1
    else
       k = prev%track(i)
       call follow_particle(k,prev%time,prev%state(:,i),now%time,now%state(:,j))
       now%track(j) = k
       i = i + 1
       j = j + 1
    endif
 enddo

end subroutine add_dump_to_table

!----------------------------------------------------------------
!+
!  Path of one particle between two dumps (1: previous, 2: this
!  one): start, extend or end its track k (0 if it has none)
!+
!----------------------------------------------------------------
subroutine follow_particle(k,t1,state1,t2,state2)
 integer,      intent(inout) :: k
 real,         intent(in)    :: t1,t2
 real(kind=4), intent(in)    :: state1(:),state2(:)
 real    :: r1,r2,t_entry,t_exit
 logical :: entered

 r1 = norm2(real(state1(1:3)))
 r2 = norm2(real(state2(1:3)))

 ! diagnostics
 if (t1 < window_start + replay_period) then
    if (r1 <  r_inj .and. r2 >= r_inj) n_out = n_out + 1
    if (r1 >= r_inj .and. r2 <  r_inj) n_in  = n_in + 1
 endif
 if (r1 > r_b .and. r1 < r_inj .and. state1(8) > 0.) then
    sum_moved = sum_moved + norm2(real(state2(1:3) - state1(1:3)))/state1(8)
    n_moved   = n_moved + 1
 endif

 entered = .false.
 t_entry = t1
 if (k == 0) then
    if (r1 >= r_b .or. r2 < r_b) return    ! does not enter the buffer
    t_entry = crossing_time(r_b,t1,state1,t2,state2)
    k = new_track(t_entry)
    call add_sample(k,t1,state1)
    call add_interpolated_sample(k,t_entry,t1,state1,t2,state2)
    entered = .true.
 endif

 if (r2 >= r_inj) then
    ! leaves through r_inj: released
    t_exit = max(crossing_time(r_inj,t1,state1,t2,state2),t_entry)
    call add_interpolated_sample(k,t_exit,t1,state1,t2,state2)
    call add_sample(k,t2,state2)
    call end_track(k,t_exit,iend_release)
 elseif (r2 < r_b) then
    if (entered) then
       ! in and out of the buffer between two dumps: ignored
       track_end(k) = iend_discard
       k = 0
    else
       ! leaves through r_b: dropped
       t_exit = crossing_time(r_b,t1,state1,t2,state2)
       call add_interpolated_sample(k,t_exit,t1,state1,t2,state2)
       call add_sample(k,t2,state2)
       call end_track(k,t_exit,iend_drop)
    endif
 else
    call add_sample(k,t2,state2)
 endif

end subroutine follow_particle

!----------------------------------------------------------------
!+
!  Time at which the path between two dumps crosses radius r
!+
!----------------------------------------------------------------
real function crossing_time(r,t1,state1,t2,state2)
 real,         intent(in) :: r,t1,t2
 real(kind=4), intent(in) :: state1(:),state2(:)
 real    :: r1,r2,f
 logical :: found

 call find_sphere_crossing(real(state1(1:3)),real(state1(4:6)),real(state2(1:3)),real(state2(4:6)), &
                           t1,t2,r,crossing_time,found)
 if (.not. found) then
    ! no sign change found: linear in radius
    r1 = norm2(real(state1(1:3)))
    r2 = norm2(real(state2(1:3)))
    f  = 0.
    if (r2 /= r1) f = min(max((r - r1)/(r2 - r1),0.),1.)
    crossing_time = t1 + f*(t2 - t1)
 endif

end function crossing_time

!----------------------------------------------------------------
!+
!  Add to track k the interpolated state at time t
!+
!----------------------------------------------------------------
subroutine add_interpolated_sample(k,t,t1,state1,t2,state2)
 integer,      intent(in) :: k
 real,         intent(in) :: t,t1,t2
 real(kind=4), intent(in) :: state1(:),state2(:)
 real :: state(size(state1))

 call interpolate_state(t1,state1,t2,state2,t,state)
 call add_sample(k,t,real(state,kind=4))

end subroutine add_interpolated_sample

!----------------------------------------------------------------
!+
!  Track k ends at time t_exit: keep the part of it inside the window
!+
!----------------------------------------------------------------
subroutine end_track(k,t_exit,iend)
 integer, intent(inout) :: k
 real,    intent(in)    :: t_exit
 integer, intent(in)    :: iend
 real :: window_end

 window_end = window_start + replay_period
 if (t_exit < window_end) then
    track_t_exit(k) = t_exit
    track_end(k)    = iend
 elseif (track_t_entry(k) < window_end) then
    track_t_exit(k) = window_end
    track_end(k)    = iend_vanish
 else
    track_end(k)    = iend_discard
 endif
 k = 0

end subroutine end_track

!----------------------------------------------------------------
!+
!  Start a new track at time t_entry
!+
!----------------------------------------------------------------
integer function new_track(t_entry) result(k)
 real, intent(in) :: t_entry

 k = n_tracks + 1
 n_tracks = k
 call grow_real(track_t_entry,k)
 call grow_real(track_t_exit,k)
 call grow_int(track_end,k)
 track_t_entry(k) = t_entry
 track_t_exit(k)  = t_entry
 track_end(k)     = iend_open

end function new_track

!----------------------------------------------------------------
!+
!  Add a sample (time t, particle state) to track k
!+
!----------------------------------------------------------------
subroutine add_sample(k,t,state)
 integer,      intent(in) :: k
 real,         intent(in) :: t
 real(kind=4), intent(in) :: state(:)
 integer :: j

 j = n_samples + 1
 n_samples = j
 call grow_int(sample_track,j)
 call grow_real(sample_time,j)
 call grow_real4(sample_state,j,size(state))
 sample_track(j)   = k
 sample_time(j)    = t
 sample_state(:,j) = state

end subroutine add_sample

!----------------------------------------------------------------
!+
!  Record the luminosity of the sinks at a dump
!+
!----------------------------------------------------------------
subroutine record_luminosity(time,lum)
 real, intent(in) :: time,lum(:)
 real, allocatable :: time_tmp(:),lum_tmp(:,:)

 if (n_dumps == 0) then
    n_sinks = size(lum)
 elseif (size(lum) /= n_sinks) then
    return   ! the number of sinks changed
 endif
 allocate(time_tmp(n_dumps+1),lum_tmp(n_sinks,n_dumps+1))
 if (n_dumps > 0) then
    time_tmp(1:n_dumps)  = dump_time(1:n_dumps)
    lum_tmp(:,1:n_dumps) = dump_lum(:,1:n_dumps)
 endif
 n_dumps = n_dumps + 1
 time_tmp(n_dumps)  = time
 lum_tmp(:,n_dumps) = lum
 call move_alloc(time_tmp,dump_time)
 call move_alloc(lum_tmp,dump_lum)

end subroutine record_luminosity

!----------------------------------------------------------------
!+
!  Finish the table: tracks still open in the last dump end at the
!  end of the window, discarded tracks are dropped, and the samples
!  of each track are put next to each other
!+
!----------------------------------------------------------------
subroutine finish_track_table(last_dump)
 type(band), intent(in) :: last_dump
 integer,      allocatable :: new_number(:),n_samples_track(:),next(:),track_end_new(:),first(:),last(:)
 real,         allocatable :: t_entry(:),t_exit(:),time(:)
 real(kind=4), allocatable :: state(:,:)
 integer :: i,j,k,m

 do i = 1,last_dump%n
    k = last_dump%track(i)
    if (k > 0) call end_track(k,window_start + replay_period,iend_vanish)
 enddo

 ! new numbers of the tracks kept (0 if dropped), and their number of samples
 allocate(new_number(n_tracks))
 m = 0
 do k = 1,n_tracks
    new_number(k) = 0
    if (track_end(k) > 0) then
       m = m + 1
       new_number(k) = m
    endif
 enddo
 allocate(n_samples_track(0:m))
 n_samples_track = 0
 do j = 1,n_samples
    n_samples_track(new_number(sample_track(j))) = n_samples_track(new_number(sample_track(j))) + 1
 enddo

 ! the tracks kept, each with its range of samples
 allocate(t_entry(m),t_exit(m),track_end_new(m),first(m),last(m),next(m))
 j = 0
 do k = 1,n_tracks
    if (new_number(k) == 0) cycle
    m = new_number(k)
    t_entry(m)       = track_t_entry(k)
    t_exit(m)        = track_t_exit(k)
    track_end_new(m) = track_end(k)
    first(m)         = j + 1
    j                = j + n_samples_track(m)
    last(m)          = j
 enddo

 ! their samples, which were recorded in time order and stay in time order
 allocate(time(j),state(size(sample_state,1),j))
 next = first
 do j = 1,n_samples
    m = new_number(sample_track(j))
    if (m == 0) cycle
    time(next(m))    = sample_time(j)
    state(:,next(m)) = sample_state(:,j)
    next(m) = next(m) + 1
 enddo

 n_tracks  = size(t_entry)
 n_samples = size(time)
 call move_alloc(t_entry,track_t_entry)
 call move_alloc(t_exit,track_t_exit)
 call move_alloc(track_end_new,track_end)
 call move_alloc(first,track_first)
 call move_alloc(last,track_last)
 call move_alloc(time,sample_time)
 call move_alloc(state,sample_state)
 if (allocated(sample_track)) deallocate(sample_track)

end subroutine finish_track_table

!----------------------------------------------------------------
!+
!  State of track k at time t (recorded time). Outside its samples
!  the position is extrapolated linearly
!+
!----------------------------------------------------------------
subroutine get_track_state(k,t,state)
 integer, intent(in)  :: k
 real,    intent(in)  :: t
 real,    intent(out) :: state(:)
 integer :: lo,hi,mid

 lo = track_first(k)
 hi = track_last(k)
 if (t <= sample_time(lo) .or. lo == hi) then
    call extrapolate_sample(lo,t,state)
 elseif (t >= sample_time(hi)) then
    call extrapolate_sample(hi,t,state)
 else
    ! the two samples around t
    do while (hi - lo > 1)
       mid = (lo + hi)/2
       if (sample_time(mid) <= t) then
          lo = mid
       else
          hi = mid
       endif
    enddo
    call interpolate_state(sample_time(lo),sample_state(:,lo),sample_time(hi),sample_state(:,hi),t,state)
 endif

end subroutine get_track_state

!----------------------------------------------------------------
!+
!  State of sample i, moved in a straight line to time t
!+
!----------------------------------------------------------------
subroutine extrapolate_sample(i,t,state)
 integer, intent(in)  :: i
 real,    intent(in)  :: t
 real,    intent(out) :: state(:)

 state      = sample_state(:,i)
 state(1:3) = state(1:3) + state(4:6)*(t - sample_time(i))

end subroutine extrapolate_sample

!----------------------------------------------------------------
!+
!  State at time t between two samples: Hermite for x and v,
!  linear for u and the scalars, log-linear for h
!+
!----------------------------------------------------------------
subroutine interpolate_state(t1,state1,t2,state2,t,state)
 real,         intent(in)  :: t1,t2,t
 real(kind=4), intent(in)  :: state1(:),state2(:)
 real,         intent(out) :: state(:)
 real :: f

 call hermite_interpolation(real(state1(1:3)),real(state1(4:6)),real(state2(1:3)),real(state2(4:6)), &
                            t1,t2,t,state(1:3),state(4:6))
 f = 0.
 if (t2 > t1) f = min(max((t - t1)/(t2 - t1),0.),1.)
 state(7) = (1. - f)*state1(7) + f*state2(7)
 if (state1(8) > 0. .and. state2(8) > 0.) then
    state(8) = exp((1. - f)*log(real(state1(8))) + f*log(real(state2(8))))
 else
    state(8) = (1. - f)*state1(8) + f*state2(8)
 endif
 state(nbase+1:) = (1. - f)*state1(nbase+1:) + f*state2(nbase+1:)

end subroutine interpolate_state

!----------------------------------------------------------------
!+
!  Cubic Hermite interpolation of a path between (x1,v1) at t1
!  and (x2,v2) at t2. The velocity is the time derivative of x
!+
!----------------------------------------------------------------
pure subroutine hermite_interpolation(x1,v1,x2,v2,t1,t2,t,x,v)
 real, intent(in)  :: x1(3),v1(3),x2(3),v2(3),t1,t2,t
 real, intent(out) :: x(3),v(3)
 real :: dt,s,s2,s3

 dt = t2 - t1
 if (dt <= 0.) then
    x = x1 + v1*(t - t1)
    v = v1
    return
 endif
 s  = (t - t1)/dt
 s2 = s*s
 s3 = s2*s
 x = (2.*s3 - 3.*s2 + 1.)*x1 + (s3 - 2.*s2 + s)*dt*v1 &
   + (-2.*s3 + 3.*s2)*x2 + (s3 - s2)*dt*v2
 v = (6.*s2 - 6.*s)*(x1 - x2)/dt + (3.*s2 - 4.*s + 1.)*v1 + (3.*s2 - 2.*s)*v2

end subroutine hermite_interpolation

!----------------------------------------------------------------
!+
!  Time tc at which the Hermite path crosses the sphere of radius
!  r, by bisection on |x|^2 - r^2 (found is false if it has the
!  same sign at both ends)
!+
!----------------------------------------------------------------
subroutine find_sphere_crossing(x1,v1,x2,v2,t1,t2,r,tc,found)
 real,    intent(in)  :: x1(3),v1(3),x2(3),v2(3),t1,t2,r
 real,    intent(out) :: tc
 logical, intent(out) :: found
 real    :: g_lo,g_hi,g_mid,t_lo,t_hi,x(3),v(3)
 integer :: iter

 g_lo  = dot_product(x1,x1) - r*r
 g_hi  = dot_product(x2,x2) - r*r
 tc    = t1
 found = (g_lo <= 0. .or. g_hi <= 0.) .and. (g_lo >= 0. .or. g_hi >= 0.)
 if (.not. found .or. abs(g_lo) <= tiny(g_lo)) return
 if (abs(g_hi) <= tiny(g_hi)) then
    tc = t2
    return
 endif
 t_lo = t1
 t_hi = t2
 do iter = 1,200
    tc = 0.5*(t_lo + t_hi)
    call hermite_interpolation(x1,v1,x2,v2,t1,t2,tc,x,v)
    g_mid = dot_product(x,x) - r*r
    if (g_mid*g_lo < 0.) then
       t_hi = tc
    else
       t_lo = tc
       g_lo = g_mid
    endif
    if (t_hi - t_lo <= 1.e-14*(t2 - t1)) exit
 enddo
 tc = 0.5*(t_lo + t_hi)

end subroutine find_sphere_crossing

!----------------------------------------------------------------
!+
!  Transformation applied to replay cycle n: none (mode 0), a
!  random rotation (mode 1, reproducible, none for n = 0), or the
!  mirror z -> -z for odd n (mode 2)
!+
!----------------------------------------------------------------
subroutine get_cycle_rotation(n,mode,rot)
 integer, intent(in)  :: n,mode
 real,    intent(out) :: rot(3,3)
 integer(kind=8), parameter :: iseed = 20261108_8
 real,            parameter :: pi = 3.141592653589793
 integer(kind=8) :: state
 real    :: u1,u2,u3,qw,qx,qy,qz
 integer :: i

 rot = 0.
 do i = 1,3
    rot(i,i) = 1.
 enddo
 select case(mode)
 case(1)
    if (n == 0) return
    state = mod(iseed*7919_8 + abs(int(n,kind=8))*104729_8 + 12345_8,2147483646_8) + 1_8
    do i = 1,9
       u1 = park_miller(state)
    enddo
    u2 = park_miller(state)
    u3 = park_miller(state)
    ! uniform random unit quaternion (Shoemake 1992)
    qw = sqrt(1. - u1)*sin(2.*pi*u2)
    qx = sqrt(1. - u1)*cos(2.*pi*u2)
    qy = sqrt(u1)*sin(2.*pi*u3)
    qz = sqrt(u1)*cos(2.*pi*u3)
    rot(1,:) = [1. - 2.*(qy*qy + qz*qz), 2.*(qx*qy - qw*qz), 2.*(qx*qz + qw*qy)]
    rot(2,:) = [2.*(qx*qy + qw*qz), 1. - 2.*(qx*qx + qz*qz), 2.*(qy*qz - qw*qx)]
    rot(3,:) = [2.*(qx*qz - qw*qy), 2.*(qy*qz + qw*qx), 1. - 2.*(qx*qx + qy*qy)]
 case(2)
    if (mod(abs(n),2) == 1) rot(3,3) = -1.
 end select

end subroutine get_cycle_rotation

!----------------------------------------------------------------
!+
!  Park-Miller random number generator: next number in (0,1)
!+
!----------------------------------------------------------------
real function park_miller(state)
 integer(kind=8), intent(inout) :: state

 state = mod(16807_8*state,2147483647_8)
 park_miller = real(state)/2147483647.

end function park_miller

!----------------------------------------------------------------
!+
!  Empty the track table
!+
!----------------------------------------------------------------
subroutine clear_table()

 table_dumpfile = ''
 table_udist    = 1.
 table_umass    = 1.
 table_utime    = 1.
 particle_mass  = 0.
 r_inj          = 0.
 r_b            = 0.
 window_start   = 0.
 replay_period  = 0.
 n_scalars      = 0
 n_tracks       = 0
 n_samples      = 0
 n_dumps        = 0
 n_sinks        = 0
 if (allocated(scalar_names))  deallocate(scalar_names)
 if (allocated(track_t_entry)) deallocate(track_t_entry)
 if (allocated(track_t_exit))  deallocate(track_t_exit)
 if (allocated(track_end))     deallocate(track_end)
 if (allocated(track_first))   deallocate(track_first)
 if (allocated(track_last))    deallocate(track_last)
 if (allocated(sample_time))   deallocate(sample_time)
 if (allocated(sample_state))  deallocate(sample_state)
 if (allocated(dump_time))     deallocate(dump_time)
 if (allocated(dump_lum))      deallocate(dump_lum)
 if (allocated(sample_track))  deallocate(sample_track)

end subroutine clear_table

!----------------------------------------------------------------
!+
!  Write the track table to file
!+
!----------------------------------------------------------------
subroutine write_track_table(ierr)
 integer, intent(out) :: ierr
 integer :: iunit

 open(newunit=iunit,file=table_file,access='stream',form='unformatted',status='replace',action='write',iostat=ierr)
 if (ierr /= 0) return
 write(iunit) table_magic,table_version
 write(iunit) table_dumpfile,table_udist,table_umass,table_utime,particle_mass,r_inj,r_b,window_start,replay_period
 write(iunit) n_scalars,n_tracks,n_samples,n_dumps,n_sinks
 write(iunit) scalar_names(1:n_scalars)
 write(iunit) track_t_entry(1:n_tracks),track_t_exit(1:n_tracks),track_end(1:n_tracks), &
              track_first(1:n_tracks),track_last(1:n_tracks)
 write(iunit) sample_time(1:n_samples),sample_state(:,1:n_samples)
 write(iunit) dump_time(1:n_dumps),dump_lum(:,1:n_dumps)
 close(iunit)

end subroutine write_track_table

!----------------------------------------------------------------
!+
!  Read the track table from file
!+
!----------------------------------------------------------------
subroutine read_track_table(ierr)
 integer, intent(out) :: ierr
 character(len=len(table_magic)) :: magic
 integer :: iunit,version

 call clear_table()
 open(newunit=iunit,file=table_file,access='stream',form='unformatted',status='old',action='read',iostat=ierr)
 if (ierr /= 0) return
 read(iunit,iostat=ierr) magic,version
 if (ierr == 0 .and. (magic /= table_magic .or. version /= table_version)) ierr = 1
 if (ierr == 0) read(iunit,iostat=ierr) table_dumpfile,table_udist,table_umass,table_utime,particle_mass, &
                                        r_inj,r_b,window_start,replay_period
 if (ierr == 0) read(iunit,iostat=ierr) n_scalars,n_tracks,n_samples,n_dumps,n_sinks
 if (ierr == 0) then
    allocate(scalar_names(n_scalars),track_t_entry(n_tracks),track_t_exit(n_tracks),track_end(n_tracks), &
             track_first(n_tracks),track_last(n_tracks),sample_time(n_samples), &
             sample_state(nbase+n_scalars,n_samples),dump_time(n_dumps),dump_lum(n_sinks,n_dumps))
    read(iunit,iostat=ierr) scalar_names
 endif
 if (ierr == 0) read(iunit,iostat=ierr) track_t_entry,track_t_exit,track_end,track_first,track_last
 if (ierr == 0) read(iunit,iostat=ierr) sample_time,sample_state
 if (ierr == 0) read(iunit,iostat=ierr) dump_time,dump_lum
 close(iunit)

end subroutine read_track_table

!----------------------------------------------------------------
!+
!  Make room for n elements (the size is doubled when needed)
!+
!----------------------------------------------------------------
subroutine grow_real(a,n)
 real, allocatable, intent(inout) :: a(:)
 integer,           intent(in)    :: n
 real, allocatable :: tmp(:)

 if (.not. allocated(a)) allocate(a(0))
 if (n <= size(a)) return
 allocate(tmp(max(n,2*size(a),1024)))
 tmp(1:size(a)) = a
 call move_alloc(tmp,a)

end subroutine grow_real

subroutine grow_int(a,n)
 integer, allocatable, intent(inout) :: a(:)
 integer,              intent(in)    :: n
 integer, allocatable :: tmp(:)

 if (.not. allocated(a)) allocate(a(0))
 if (n <= size(a)) return
 allocate(tmp(max(n,2*size(a),1024)))
 tmp(1:size(a)) = a
 call move_alloc(tmp,a)

end subroutine grow_int

subroutine grow_real4(a,n,nrow)
 real(kind=4), allocatable, intent(inout) :: a(:,:)
 integer,                   intent(in)    :: n,nrow
 real(kind=4), allocatable :: tmp(:,:)

 if (.not. allocated(a)) allocate(a(nrow,0))
 if (n <= size(a,2)) return
 allocate(tmp(nrow,max(n,2*size(a,2),1024)))
 tmp(:,1:size(a,2)) = a
 call move_alloc(tmp,a)

end subroutine grow_real4

!----------------------------------------------------------------
!+
!  Replay the tracks up to cycle ncycle: normally only that cycle,
!  two cycles when one has ended since the last step
!+
!----------------------------------------------------------------
subroutine replay_tracks(time,ncycle,irotate,xc,vc,npart,npartoftype,xyzh,vxyzu,xyzmh_ptmass)
 real,    intent(in)    :: time,xc(3),vc(3)
 integer, intent(in)    :: ncycle,irotate
 integer, intent(inout) :: npart,npartoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:)
 integer :: icycle

 do icycle = cycle_last,ncycle
    call replay_cycle(icycle,time,irotate,xc,vc,npart,npartoftype,xyzh,vxyzu,xyzmh_ptmass)
    if (icycle < ncycle) call report_cycle(icycle)
 enddo
 cycle_last = ncycle

end subroutine replay_tracks

!----------------------------------------------------------------
!+
!  Create, move, release and remove the particles of the tracks
!  in cycle icycle
!+
!----------------------------------------------------------------
subroutine replay_cycle(icycle,time,irotate,xc,vc,npart,npartoftype,xyzh,vxyzu,xyzmh_ptmass)
 use part, only:iboundary,kill_particle
 integer, intent(in)    :: icycle,irotate
 real,    intent(in)    :: time,xc(3),vc(3)
 integer, intent(inout) :: npart,npartoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:)
 real    :: tau,rot(3,3)
 integer :: k

 tau = time - icycle*replay_period
 call get_cycle_rotation(icycle,irotate,rot)
 do k = 1,n_tracks
    if (track_cycle(k) < icycle .and. track_t_entry(k) <= tau) then
       ! the recorded particle enters the buffer
       track_cycle(k) = icycle
       if (tau < track_t_exit(k)) then
          call spawn_particle(k,tau,rot,xc,vc,npart,npartoftype,xyzh,vxyzu)
          n_spawned = n_spawned + 1
       elseif (track_end(k) == iend_release) then
          ! the whole track was passed since the last step
          call release_particle(npart+1,k,tau,rot,xc,vc,npart,npartoftype,xyzh,vxyzu,xyzmh_ptmass)
       endif
    elseif (track_particle(k) > 0) then
       if (tau < track_t_exit(k)) then
          call put_particle_on_track(track_particle(k),iboundary,k,tau,rot,xc,vc,npart,npartoftype,xyzh,vxyzu)
       elseif (track_end(k) == iend_release) then
          call release_particle(track_particle(k),k,tau,rot,xc,vc,npart,npartoftype,xyzh,vxyzu,xyzmh_ptmass)
          track_particle(k) = 0
       else
          call kill_particle(track_particle(k),npartoftype)
          track_particle(k) = 0
          n_removed = n_removed + 1
       endif
    endif
 enddo

end subroutine replay_cycle

!----------------------------------------------------------------
!+
!  First call of a run (start or restart): the buffer is rebuilt
!  from the table, so nothing needs to be saved in the dumps
!+
!----------------------------------------------------------------
subroutine start_replay(time,ncycle,npart,npartoftype,xyzh)
 real,    intent(in)    :: time
 integer, intent(in)    :: ncycle,npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: xyzh(:,:)
 real    :: tau
 integer :: k,n_boundary_removed

 call remove_boundary_particles(npart,npartoftype,xyzh,n_boundary_removed)
 tau = time - ncycle*replay_period
 do k = 1,n_tracks
    track_particle(k) = 0
    if (track_t_exit(k) <= tau) then
       track_cycle(k) = ncycle       ! already over in this cycle
    else
       track_cycle(k) = ncycle - 1   ! started when tau reaches its entry time
    endif
 enddo
 cycle_last     = ncycle
 replay_started = .true.

 print *, ''
 print *, 'Replay starts at t =',time,' in cycle',ncycle
 print *, 'Boundary particles removed:',n_boundary_removed

end subroutine start_replay

!----------------------------------------------------------------
!+
!  The particle arrays are reordered when dead particles are
!  removed (before a dump is written): check that the buffer
!  particles are still in place
!+
!----------------------------------------------------------------
logical function buffer_is_intact(npart,xyzh)
 use part, only:iorig,iphase,iamtype,iboundary,isdead_or_accreted
 integer, intent(in) :: npart
 real,    intent(in) :: xyzh(:,:)
 integer :: i,k

 buffer_is_intact = .true.
 do k = 1,n_tracks
    i = track_particle(k)
    if (i == 0) cycle
    if (i > npart) then
       buffer_is_intact = .false.
    elseif (iorig(i) /= track_particle_id(k) .or. iamtype(iphase(i)) /= iboundary .or. &
            isdead_or_accreted(xyzh(4,i))) then
       buffer_is_intact = .false.
    endif
    if (.not. buffer_is_intact) return
 enddo

end function buffer_is_intact

!----------------------------------------------------------------
!+
!  Create the buffer particles again
!+
!----------------------------------------------------------------
subroutine rebuild_buffer(time,irotate,xc,vc,npart,npartoftype,xyzh,vxyzu)
 real,    intent(in)    :: time,xc(3),vc(3)
 integer, intent(in)    :: irotate
 integer, intent(inout) :: npart,npartoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real    :: tau,rot(3,3)
 integer :: k,n_boundary_removed

 call remove_boundary_particles(npart,npartoftype,xyzh,n_boundary_removed)
 tau = time - cycle_last*replay_period
 call get_cycle_rotation(cycle_last,irotate,rot)
 do k = 1,n_tracks
    if (track_particle(k) > 0) call spawn_particle(k,tau,rot,xc,vc,npart,npartoftype,xyzh,vxyzu)
 enddo

end subroutine rebuild_buffer

!----------------------------------------------------------------
!+
!  New boundary particle on track k
!+
!----------------------------------------------------------------
subroutine spawn_particle(k,tau,rot,xc,vc,npart,npartoftype,xyzh,vxyzu)
 use part, only:iorig,iboundary
 integer, intent(in)    :: k
 real,    intent(in)    :: tau,rot(3,3),xc(3),vc(3)
 integer, intent(inout) :: npart,npartoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)

 call put_particle_on_track(npart+1,iboundary,k,tau,rot,xc,vc,npart,npartoftype,xyzh,vxyzu)
 track_particle(k)    = npart
 track_particle_id(k) = iorig(npart)

end subroutine spawn_particle

!----------------------------------------------------------------
!+
!  Particle i (npart+1 for a new one) on track k becomes gas,
!  taking its mass from the sink
!+
!----------------------------------------------------------------
subroutine release_particle(i,k,tau,rot,xc,vc,npart,npartoftype,xyzh,vxyzu,xyzmh_ptmass)
 use part, only:igas,nptmass
 integer, intent(in)    :: i,k
 real,    intent(in)    :: tau,rot(3,3),xc(3),vc(3)
 integer, intent(inout) :: npart,npartoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:)

 call put_particle_on_track(i,igas,k,tau,rot,xc,vc,npart,npartoftype,xyzh,vxyzu)
 if (nptmass > 0) xyzmh_ptmass(4,1) = xyzmh_ptmass(4,1) - particle_mass
 n_released = n_released + 1

end subroutine release_particle

!----------------------------------------------------------------
!+
!  Put particle i (npart+1 for a new one) on track k at time tau,
!  as a particle of type itype
!+
!----------------------------------------------------------------
subroutine put_particle_on_track(i,itype,k,tau,rot,xc,vc,npart,npartoftype,xyzh,vxyzu)
 use part,       only:iboundary,igas,iphase,iunknown
 use partinject, only:add_or_update_particle
 integer, intent(in)    :: i,itype,k
 real,    intent(in)    :: tau,rot(3,3),xc(3),vc(3)
 integer, intent(inout) :: npart,npartoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real :: state(nbase+n_scalars),x(3),v(3)

 call get_track_state(k,tau,state)
 x = xc + matmul(rot,state(1:3))
 v = vc + matmul(rot,state(4:6))
 if (itype == iboundary .and. i <= npart) then
    ! a boundary particle moved along its track (add_or_update_particle would not update h)
    xyzh(1:3,i)  = x
    xyzh(4,i)    = state(8)
    vxyzu(1:3,i) = v
    if (size(vxyzu,1) >= 4) vxyzu(4,i) = state(7)
 else
    call add_or_update_particle(itype,x,v,state(8),state(7),i,npart,npartoftype,xyzh,vxyzu)
 endif
 ! after add_or_update_particle, which resets mu, gamma and alpha
 call apply_scalars(i,state)
 ! so that update_injected_particles sets the timestep of the new gas particle
 if (itype == igas) iphase(i) = int(iunknown,kind=kind(iphase))

end subroutine put_particle_on_track

!----------------------------------------------------------------
!+
!  Copy the recorded per-particle quantities to particle i
!+
!----------------------------------------------------------------
subroutine apply_scalars(i,state)
 use dim,  only:update_muGamma,maxalpha,maxp,store_dust_temperature,h2chemistry
 use part, only:eos_vars,imu,igamma,alphaind,dust_temp,abundance,nabundances
 integer, intent(in) :: i
 real,    intent(in) :: state(:)

 if (update_muGamma .and. irow_mu > 0)    eos_vars(imu,i)    = state(irow_mu)
 if (update_muGamma .and. irow_gamma > 0) eos_vars(igamma,i) = state(irow_gamma)
 if (maxalpha == maxp .and. irow_alpha > 0) alphaind(1,i) = real(state(irow_alpha),kind=kind(alphaind))
 if (store_dust_temperature .and. irow_tdust > 0) dust_temp(i) = state(irow_tdust)
 if (h2chemistry .and. irow_abundance > 0) &
    abundance(1:nabundances,i) = state(irow_abundance:irow_abundance+nabundances-1)

end subroutine apply_scalars

!----------------------------------------------------------------
!+
!  Remove all boundary particles
!+
!----------------------------------------------------------------
subroutine remove_boundary_particles(npart,npartoftype,xyzh,n_boundary_removed)
 use part, only:iphase,iamtype,iboundary,isdead_or_accreted,kill_particle
 integer, intent(in)    :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: xyzh(:,:)
 integer, intent(out)   :: n_boundary_removed
 integer :: i

 n_boundary_removed = 0
 do i = 1,npart
    if (isdead_or_accreted(xyzh(4,i))) cycle
    if (iamtype(iphase(i)) == iboundary) then
       call kill_particle(i,npartoftype)
       n_boundary_removed = n_boundary_removed + 1
    endif
 enddo

end subroutine remove_boundary_particles

!----------------------------------------------------------------
!+
!  Remove the gas that has fallen back inside r_b
!+
!----------------------------------------------------------------
subroutine kill_inner_gas(xc,npart,npartoftype,xyzh)
 use part, only:iphase,iamtype,igas,isdead_or_accreted,kill_particle
 real,    intent(in)    :: xc(3)
 integer, intent(in)    :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: xyzh(:,:)
 real    :: dx(3)
 integer :: i

 do i = 1,npart
    if (isdead_or_accreted(xyzh(4,i))) cycle
    if (iamtype(iphase(i)) /= igas) cycle
    dx = xyzh(1:3,i) - xc
    if (dot_product(dx,dx) < r_b**2) then
       call kill_particle(i,npartoftype)
       n_killed = n_killed + 1
    endif
 enddo

end subroutine kill_inner_gas

!----------------------------------------------------------------
!+
!  Sink luminosity at the same phase of the recorded window (if
!  the sinks were merged, the luminosities of the inner sinks are
!  summed)
!+
!----------------------------------------------------------------
subroutine replay_luminosity(time,xyzmh_ptmass)
 use part, only:nptmass,ilum
 real, intent(in)    :: time
 real, intent(inout) :: xyzmh_ptmass(:,:)
 real    :: t,f
 integer :: j

 if (n_dumps < 2 .or. n_sinks < 1 .or. nptmass < 1) return
 t = window_start + modulo(time - window_start,replay_period)
 j = 1
 do while (j < n_dumps - 1 .and. dump_time(j+1) <= t)
    j = j + 1
 enddo
 f = min(max((t - dump_time(j))/(dump_time(j+1) - dump_time(j)),0.),1.)
 if (nptmass == n_sinks) then
    xyzmh_ptmass(ilum,1:nptmass) = (1. - f)*dump_lum(:,j) + f*dump_lum(:,j+1)
 else
    xyzmh_ptmass(ilum,1) = sum((1. - f)*dump_lum(:,j) + f*dump_lum(:,j+1))
 endif

end subroutine replay_luminosity

!----------------------------------------------------------------
!+
!  Centre of mass of the sinks (the tracks are relative to it)
!+
!----------------------------------------------------------------
subroutine get_centre_of_mass(nsink,xyzmh,vxyz,xc,vc)
 integer, intent(in)  :: nsink
 real,    intent(in)  :: xyzmh(:,:),vxyz(:,:)
 real,    intent(out) :: xc(3),vc(3)
 real    :: mtot
 integer :: k

 xc   = 0.
 vc   = 0.
 mtot = 0.
 do k = 1,nsink
    if (xyzmh(4,k) <= 0.) cycle
    xc   = xc + xyzmh(4,k)*xyzmh(1:3,k)
    vc   = vc + xyzmh(4,k)*vxyz(1:3,k)
    mtot = mtot + xyzmh(4,k)
 enddo
 if (mtot > 0.) then
    xc = xc/mtot
    vc = vc/mtot
 endif

end subroutine get_centre_of_mass

!----------------------------------------------------------------
!+
!  Merge the sinks into one at their centre of mass, with their
!  total mass and luminosity (the other properties are those of
!  sink 1)
!+
!----------------------------------------------------------------
subroutine merge_sinks(nptmass,xyzmh_ptmass,vxyz_ptmass)
 use part, only:ilum
 integer, intent(inout) :: nptmass
 real,    intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 real :: xc(3),vc(3)

 call get_centre_of_mass(nptmass,xyzmh_ptmass,vxyz_ptmass,xc,vc)
 xyzmh_ptmass(1:3,1)  = xc
 vxyz_ptmass(1:3,1)   = vc
 xyzmh_ptmass(4,1)    = sum(xyzmh_ptmass(4,1:nptmass))
 xyzmh_ptmass(ilum,1) = sum(xyzmh_ptmass(ilum,1:nptmass))
 xyzmh_ptmass(:,2:nptmass) = 0.
 vxyz_ptmass(:,2:nptmass)  = 0.
 nptmass = 1

end subroutine merge_sinks

!----------------------------------------------------------------
!+
!  Summary of a replay cycle
!+
!----------------------------------------------------------------
subroutine report_cycle(icycle)
 use units,   only:umass,utime
 use physcon, only:solarm,years
 integer, intent(in) :: icycle

 print *, ''
 print *, ' Replay cycle               :', icycle
 print *, ' Particles spawned          :', n_spawned
 print *, ' Particles released         :', n_released
 print *, ' Removed at end of track    :', n_removed
 print *, ' Fallen back inside r_b     :', n_killed
 print *, ' Replayed Mdot (Msun/yr)    :', n_released*particle_mass/replay_period*(umass/solarm)/(utime/years)
 print *, ''
 n_spawned  = 0
 n_released = 0
 n_removed  = 0
 n_killed   = 0

end subroutine report_cycle

!----------------------------------------------------------------
!+
!  The q-th quantile (0 <= q <= 1) of a
!+
!----------------------------------------------------------------
real function percentile(a,q)
 use sortutils, only:indexx
 real, intent(in) :: a(:),q
 integer, allocatable :: isort(:)
 integer :: n

 n = size(a)
 percentile = 0.
 if (n < 1) return
 allocate(isort(n))
 call indexx(n,a,isort)
 percentile = a(isort(min(max(1,nint(q*(n-1)) + 1),n)))

end function percentile

end module replayutils
