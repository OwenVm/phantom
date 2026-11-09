!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module testreplay
!
! Unit tests of the periodic track replay routines (replayutils)
!
! :References: None
!
! :Owner: Owen Vermeulen
!
! :Runtime parameters: None
!
! :Dependencies: io, replayutils, testutils
!
 use replayutils
 implicit none

 public :: test_replay

 private

contains

subroutine test_replay(ntests,npass)
 use io, only:id,master
 integer, intent(inout) :: ntests,npass

 if (id==master) write(*,"(/,a,/)") '--> TESTING REPLAY ROUTINES'

 call test_hermite(ntests,npass)
 call test_crossing(ntests,npass)
 call test_flux(ntests,npass)
 call test_rotation(ntests,npass)
 call test_schedule(ntests,npass)
 call test_table_io(ntests,npass)
 call test_pipeline(ntests,npass)

 if (id==master) write(*,"(/,a)") '<-- REPLAY TESTS COMPLETE'

end subroutine test_replay

!-----------------------------------------------------------------------
!+
!  Hermite interpolation is exact for cubic paths
!+
!-----------------------------------------------------------------------
subroutine test_hermite(ntests,npass)
 use testutils, only:checkval,update_test_scores
 integer, intent(inout) :: ntests,npass
 real(dp) :: c(3,0:3),ta,tb,t,xa(3),va(3),xb(3),vb(3),x(3),v(3),xe(3),ve(3)
 integer  :: nfailed(2)

 c(:,0) = [0.3_dp,-1.2_dp,0.7_dp]
 c(:,1) = [1.1_dp,0.4_dp,-0.5_dp]
 c(:,2) = [-0.6_dp,0.9_dp,0.2_dp]
 c(:,3) = [0.25_dp,-0.15_dp,0.4_dp]
 ta = 0.3_dp
 tb = 1.1_dp
 t  = 0.77_dp
 call cubic(ta,xa,va)
 call cubic(tb,xb,vb)
 call cubic(t,xe,ve)
 call hermite_eval(xa,va,xb,vb,ta,tb,t,x,v)
 nfailed = 0
 call checkval(3,x,xe,1.e-12,nfailed(1),'Hermite position for a cubic path')
 call checkval(3,v,ve,1.e-12,nfailed(2),'Hermite velocity for a cubic path')
 call update_test_scores(ntests,nfailed,npass)

contains

 subroutine cubic(tt,xx,vv)
  real(dp), intent(in)  :: tt
  real(dp), intent(out) :: xx(3),vv(3)

  xx = c(:,0) + c(:,1)*tt + c(:,2)*tt**2 + c(:,3)*tt**3
  vv = c(:,1) + 2._dp*c(:,2)*tt + 3._dp*c(:,3)*tt**2

 end subroutine cubic

end subroutine test_hermite

!-----------------------------------------------------------------------
!+
!  Crossing time of a straight line with a sphere
!+
!-----------------------------------------------------------------------
subroutine test_crossing(ntests,npass)
 use testutils, only:checkval,update_test_scores
 integer, intent(inout) :: ntests,npass
 real(dp) :: x0(3),v0(3),r,tc,te,a,b,c,ta,tb
 logical  :: found
 integer  :: nfailed(3)

 x0 = [0.2_dp,0.1_dp,0._dp]
 v0 = [1._dp,0.5_dp,0.2_dp]
 r  = 1._dp
 ta = 0._dp
 tb = 1.5_dp
 ! |x0 + v0 t|^2 = r^2
 a  = dot_product(v0,v0)
 b  = 2._dp*dot_product(x0,v0)
 c  = dot_product(x0,x0) - r*r
 te = (-b + sqrt(b*b - 4._dp*a*c))/(2._dp*a)
 call find_sphere_crossing(x0,v0,x0+v0*tb,v0,ta,tb,r,tc,found)
 nfailed = 0
 call checkval(found,.true.,nfailed(1),'crossing found')
 call checkval(tc,te,1.e-10_dp,nfailed(2),'crossing time of a straight line')
 call find_sphere_crossing(x0,v0,x0+v0*0.1_dp,v0,0._dp,0.1_dp,r,tc,found)
 call checkval(found,.false.,nfailed(3),'no crossing inside the sphere')
 call update_test_scores(ntests,nfailed,npass)

end subroutine test_crossing

!-----------------------------------------------------------------------
!+
!  Shell estimator of the mass flux
!+
!-----------------------------------------------------------------------
subroutine test_flux(ntests,npass)
 use testutils, only:checkval,update_test_scores
 integer, intent(inout) :: ntests,npass
 integer, parameter :: n = 1000
 real(dp) :: r(n),vr(n),mout,min_in,vrf
 integer  :: i,nin,nfailed(3)

 nin = 0
 do i = 1,n
    r(i)  = 0.5_dp + real(i,dp)/n
    vr(i) = 2._dp
    if (mod(i,10) == 0) vr(i) = -1._dp
 enddo
 nin = count(abs(r - 1._dp) < 0.05_dp .and. vr > 0._dp)
 call shell_flux(n,r,vr,0.01_dp,1._dp,0.1_dp,mout,min_in,vrf)
 nfailed = 0
 call checkval(mout,0.01_dp*2._dp*nin/0.1_dp,1.e-12_dp,nfailed(1),'outward mass flux')
 call checkval(vrf,2._dp,1.e-12_dp,nfailed(2),'flux-weighted outflow speed')
 call checkval(min_in,0.01_dp*count(abs(r - 1._dp) < 0.05_dp .and. vr < 0._dp)/0.1_dp, &
               1.e-12_dp,nfailed(3),'inward mass flux')
 call update_test_scores(ntests,nfailed,npass)

end subroutine test_flux

!-----------------------------------------------------------------------
!+
!  Per-cycle transformations
!+
!-----------------------------------------------------------------------
subroutine test_rotation(ntests,npass)
 use testutils, only:checkval,update_test_scores
 integer, intent(inout) :: ntests,npass
 real(dp) :: rot(3,3),rot2(3,3),eye(3,3),rr(3,3),det
 integer  :: i,nfailed(5)

 eye = 0._dp
 do i = 1,3
    eye(i,i) = 1._dp
 enddo
 nfailed = 0
 call cycle_rotation(0,irot_random,42,rot)
 call checkval(9,reshape(rot,[9]),reshape(eye,[9]),1.e-14,nfailed(1),'cycle 0 is not rotated')
 call cycle_rotation(3,irot_random,42,rot)
 rr  = matmul(transpose(rot),rot)
 det = rot(1,1)*(rot(2,2)*rot(3,3) - rot(2,3)*rot(3,2)) &
     - rot(1,2)*(rot(2,1)*rot(3,3) - rot(2,3)*rot(3,1)) &
     + rot(1,3)*(rot(2,1)*rot(3,2) - rot(2,2)*rot(3,1))
 call checkval(9,reshape(rr,[9]),reshape(eye,[9]),1.e-12,nfailed(2),'rotation is orthogonal')
 call checkval(det,1._dp,1.e-12_dp,nfailed(3),'rotation is proper')
 call cycle_rotation(3,irot_random,42,rot2)
 call checkval(9,reshape(rot2,[9]),reshape(rot,[9]),0.,nfailed(4),'rotation is reproducible')
 call cycle_rotation(1,irot_mirror,42,rot)
 eye(3,3) = -1._dp
 call checkval(9,reshape(rot,[9]),reshape(eye,[9]),0.,nfailed(5),'mirror for odd cycles')
 call update_test_scores(ntests,nfailed,npass)

end subroutine test_rotation

!-----------------------------------------------------------------------
!+
!  The schedule spawns every instance once and on time, and the active
!  set at any time equals the one maintained incrementally
!+
!-----------------------------------------------------------------------
subroutine test_schedule(ntests,npass)
 use testutils, only:checkval,update_test_scores
 integer, intent(inout) :: ntests,npass
 integer, parameter :: ntr = 300, nmaxcyc = 12
 type(track_table)     :: tab
 type(replay_schedule) :: sch
 logical  :: spawned(ntr,0:nmaxcyc),alive(ntr,0:nmaxcyc),found
 integer, allocatable :: klist(:),nlist(:)
 integer  :: k,n,nact,istep,ndup,nlate,nmismatch,nmissed,nfailed(4),i
 integer(i8) :: iseed
 real(dp) :: t,dt,t_ep,a,b,tstart,tend

 tab%t0    = 5._dp
 tab%t_rec = 1._dp
 tab%ntracks = ntr
 allocate(tab%t_entry(ntr),tab%t_exit(ntr))
 iseed = 12345_i8
 do k = 1,ntr
    a = -2.3_dp + 3.3_dp*rnd()
    b = max(a,0._dp) + (1._dp - max(a,0._dp))*rnd()
    if (b >= 1._dp) b = 0.999_dp
    tab%t_entry(k) = tab%t0 + a
    tab%t_exit(k)  = tab%t0 + b
 enddo
 t_ep   = 10._dp
 dt     = 0.003_dp
 tstart = t_ep - 0.4_dp
 tend   = t_ep + 6._dp
 call schedule_init(sch,tab,t_ep)

 spawned = .false.
 alive   = .false.
 call active_instances_at(sch,tab,tstart,nact,klist,nlist)
 do i = 1,nact
    spawned(klist(i),nlist(i)) = .true.
    alive(klist(i),nlist(i))   = .true.
 enddo
 call schedule_reset(sch,tstart)

 ndup = 0; nlate = 0; nmismatch = 0
 t = tstart
 istep = 0
 do while (t < tend)
    istep = istep + 1
    t = tstart + istep*dt
    do
       call schedule_next_spawn(sch,t,k,n,found)
       if (.not.found) exit
       if (n > nmaxcyc) cycle
       if (spawned(k,n)) ndup = ndup + 1
       spawned(k,n) = .true.
       alive(k,n)   = .true.
       ! spawn time must lie in (t-dt,t]
       a = t_ep + n*tab%t_rec + tab%t_entry(k) - tab%t0
       if (a > t .or. a <= t - dt - 1.e-12_dp) nlate = nlate + 1
    enddo
    do k = 1,ntr
       do n = 0,nmaxcyc
          if (alive(k,n)) then
             if (t >= t_ep + n*tab%t_rec + tab%t_exit(k) - tab%t0) alive(k,n) = .false.
          endif
       enddo
    enddo
    if (mod(istep,97) == 0) then
       call active_instances_at(sch,tab,t,nact,klist,nlist)
       if (nact /= count(alive)) nmismatch = nmismatch + 1
       do i = 1,nact
          if (nlist(i) <= nmaxcyc) then
             if (.not.alive(klist(i),nlist(i))) nmismatch = nmismatch + 1
          endif
       enddo
    endif
 enddo
 ! every instance that started before tend must have been spawned
 nmissed = 0
 do k = 1,ntr
    do n = 0,nmaxcyc
       a = t_ep + n*tab%t_rec + tab%t_entry(k) - tab%t0
       b = t_ep + n*tab%t_rec + tab%t_exit(k)  - tab%t0
       if (a <= tend - dt .and. b > tstart .and. .not.spawned(k,n)) nmissed = nmissed + 1
    enddo
 enddo
 nfailed = 0
 call checkval(ndup,0,0,nfailed(1),'no instance spawned twice')
 call checkval(nlate,0,0,nfailed(2),'instances spawned on time')
 call checkval(nmissed,0,0,nfailed(3),'no instance missed')
 call checkval(nmismatch,0,0,nfailed(4),'active set from schedule = incremental set')
 call update_test_scores(ntests,nfailed,npass)

contains

 real(dp) function rnd()
  iseed = mod(16807_i8*iseed,2147483647_i8)
  rnd = real(iseed,dp)/2147483647._dp
 end function rnd

end subroutine test_schedule

!-----------------------------------------------------------------------
!+
!  Write and read back a track table
!+
!-----------------------------------------------------------------------
subroutine test_table_io(ntests,npass)
 use testutils, only:checkval,update_test_scores
 integer, intent(inout) :: ntests,npass
 type(track_table) :: tab,tab2
 character(len=*), parameter :: fname = 'test_replay_table.dat'
 integer :: ierr,iunit,nfailed(8),j

 tab%udist = 1.5e13_dp; tab%umass = 2.e33_dp; tab%utime = 5.e6_dp
 tab%pmass = 1.e-10_dp; tab%r_inj = 50._dp; tab%r_b = 45._dp
 tab%t0 = 100._dp; tab%t_rec = 5._dp; tab%period = 5._dp
 tab%nscal = 2; tab%ntracks = 2; tab%nsamp = 5; tab%nsink = 1; tab%nsinkrec = 3
 allocate(tab%label(2))
 tab%label = ['nuc:K3   ','alpha    ']
 allocate(tab%t_entry(2),tab%t_exit(2),tab%iend(2),tab%isamp1(2),tab%nsampk(2))
 tab%t_entry = [99._dp,101._dp]
 tab%t_exit  = [100.5_dp,103._dp]
 tab%iend    = [iend_release,iend_drop]
 tab%isamp1  = [1,3]
 tab%nsampk  = [2,3]
 allocate(tab%ts(5),tab%xs(3,5),tab%vs(3,5),tab%us(5),tab%hs(5),tab%scal(2,5))
 do j = 1,5
    tab%ts(j)     = 98._dp + j
    tab%xs(:,j)   = real([j,2*j,3*j],sp)
    tab%vs(:,j)   = real([-j,j,0],sp)
    tab%us(j)     = real(j,sp)
    tab%hs(j)     = real(0.1_dp*j,sp)
    tab%scal(:,j) = real([j,j*j],sp)
 enddo
 allocate(tab%tsink(3),tab%sink(nsinkcol,1,3))
 tab%tsink = [99._dp,101._dp,104._dp]
 tab%sink  = 1._dp

 call write_track_table(fname,tab,ierr)
 call read_track_table(fname,tab2,ierr)
 nfailed = 0
 call checkval(ierr,0,0,nfailed(1),'read track table')
 call checkval(tab2%ntracks,2,0,nfailed(2),'number of tracks')
 call checkval(tab2%t0,tab%t0,0._dp,nfailed(3),'window start')
 call checkval(2,tab2%t_exit,tab%t_exit,0.,nfailed(4),'exit times')
 call checkval(5,tab2%ts,tab%ts,0.,nfailed(5),'sample times')
 call checkval(15,reshape(tab2%xs,[15]),reshape(tab%xs,[15]),0.,nfailed(6),'positions')
 call checkval(tab2%label(2),tab%label(2),nfailed(7),'labels')
 call checkval(tab2%nsinkrec,3,0,nfailed(8),'sink history')
 call update_test_scores(ntests,nfailed,npass)

 open(newunit=iunit,file=fname,status='old',iostat=ierr)
 if (ierr == 0) close(iunit,status='delete')

end subroutine test_table_io

!-----------------------------------------------------------------------
!+
!  Table from exactly one period of records of a synthetic radial
!  outflow, with the tracks clipped to the window. Checked: the release
!  tracks and their exit times, and that at any phase the replayed buffer
!  holds as many particles as the flow it was recorded from
!+
!-----------------------------------------------------------------------
subroutine test_pipeline(ntests,npass)
 use testutils, only:checkval,update_test_scores
 integer, intent(inout) :: ntests,npass
 integer,  parameter :: np = 600
 real(dp), parameter :: vel = 1._dp, r_inj = 1._dp, r_b = 0.8_dp, dtr = 0.01_dp
 real(dp), parameter :: t0 = 1._dp, trec = 1._dp, band_in = 0.6_dp, band_out = 1.1_dp
 type(track_table)     :: head,tab
 type(replay_builder)  :: b
 type(replay_record)   :: rec
 type(replay_schedule) :: sch
 real(dp) :: r0(np),dir(3,np),t,r,expected(np),texit(np),x(3),v(3),u,h,s(1),errmax,tau
 integer  :: i,j,k,nexp,nrec,nfailed(7),nmis,nact,nin,nappear,nvanish
 integer, allocatable :: idx(:),klist(:),nlist(:)
 integer(i8) :: iseed

 iseed = 4321_i8
 do i = 1,np
    r0(i) = -2._dp + 2.45_dp*rnd()
    dir(:,i) = [rnd() - 0.5_dp,rnd() - 0.5_dp,rnd() - 0.5_dp]
    dir(:,i) = dir(:,i)/norm2(dir(:,i))
 enddo

 head%nscal = 1
 allocate(head%label(1))
 head%label(1) = 'Tdust'
 head%pmass  = 1.e-3_dp
 head%r_inj  = r_inj
 head%r_b    = r_b
 head%t0     = t0
 head%t_rec  = trec
 head%period = trec
 call builder_init(b,head,0._dp,.true.,.true.)
 nrec = nint(trec/dtr)
 do j = 0,nrec
    call make_record(t0 + j*dtr)
    call builder_add_record(b,rec)
 enddo
 call builder_finalise(b,tab)

 ! exact: r = r0 + vel t crosses r_inj at t = (r_inj - r0)/vel
 nexp = 0
 nappear = 0
 nvanish = 0
 do i = 1,np
    t = (r_inj - r0(i))/vel
    if (t >= t0 .and. t < t0 + trec) then
       nexp = nexp + 1
       expected(nexp) = t
    endif
    r = r0(i) + vel*t0
    if (r > r_b .and. r < r_inj) nappear = nappear + 1
    r = r0(i) + vel*(t0 + trec)
    if (r > r_b .and. r < r_inj) nvanish = nvanish + 1
 enddo
 nfailed = 0
 call checkval(count(tab%iend(1:tab%ntracks) == iend_release),nexp,0,nfailed(1),'number of release tracks')
 call checkval(b%nappear,nappear,0,nfailed(2),'tracks already in the buffer at t0')
 call checkval(b%nvanish,nvanish,0,nfailed(3),'tracks still in the buffer at t0+t_rec')

 errmax = 0._dp
 if (count(tab%iend(1:tab%ntracks) == iend_release) == nexp .and. nexp > 0) then
    allocate(idx(nexp))
    call sort_index_dp(nexp,expected,idx)
    texit(1:nexp) = expected(idx)
    k = 0
    do j = 1,tab%ntracks
       if (tab%iend(j) /= iend_release) cycle
       k = k + 1
       expected(k) = tab%t_exit(j)
    enddo
    call sort_index_dp(nexp,expected,idx)
    do j = 1,nexp
       errmax = max(errmax,abs(expected(idx(j)) - texit(j)))
    enddo
 endif
 call checkval(errmax,0._dp,1.e-6_dp,nfailed(4),'exit times')

 ! the interpolated paths of the release tracks stay on r = r0 + vel t
 errmax = 0._dp
 do k = 1,tab%ntracks
    if (tab%iend(k) /= iend_release) cycle
    do j = 1,9
       tau = tab%t_entry(k) + 0.1_dp*j*(tab%t_exit(k) - tab%t_entry(k))
       call interp_track(tab,k,tau,x,v,u,h,s)
       errmax = max(errmax,abs(norm2(v) - vel))
       errmax = max(errmax,abs(norm2(x) - (r_inj - vel*(tab%t_exit(k) - tau))))
    enddo
 enddo
 call checkval(errmax,0._dp,5.e-5_dp,nfailed(5),'interpolated paths')

 ! at any phase, the replayed buffer holds the particles of the flow
 call schedule_init(sch,tab,t0)
 nmis = 0
 do j = 1,37
    t = t0 + (j - 0.5_dp)*trec/37._dp
    call active_instances_at(sch,tab,t,nact,klist,nlist)
    nin = 0
    do i = 1,np
       r = r0(i) + vel*t
       if (r > r_b .and. r < r_inj) nin = nin + 1
    enddo
    if (nact /= nin) nmis = nmis + 1
 enddo
 call checkval(nmis,0,0,nfailed(6),'replayed buffer = buffer of the flow at all phases')
 call checkval(all(tab%t_entry(1:tab%ntracks) >= t0 .and. tab%t_exit(1:tab%ntracks) <= t0 + trec), &
               .true.,nfailed(7),'all tracks lie inside the window')
 call update_test_scores(ntests,nfailed,npass)

contains

 subroutine make_record(tt)
  real(dp), intent(in) :: tt
  integer :: ii,m,jj

  rec%time = tt
  rec%xc = 0._dp
  rec%vc = 0._dp
  call alloc_record(rec,1,np,1)
  rec%nsink = 1
  rec%sink(:,1) = 0._dp
  rec%sink(isk_m,1) = 1._dp
  m = 0
  ! store in reverse order to exercise the sorting by iorig
  do jj = np,1,-1
     ii = jj
     r = r0(ii) + vel*tt
     if (r < band_in .or. r > band_out) cycle
     m = m + 1
     rec%iorig(m)  = int(ii,i8)
     rec%x(:,m)    = real(r*dir(:,ii),sp)
     rec%v(:,m)    = real(vel*dir(:,ii),sp)
     rec%u(m)      = 1._sp
     rec%h(m)      = 0.05_sp
     rec%scal(1,m) = real(r,sp)
  enddo
  rec%nsel = m

 end subroutine make_record

 real(dp) function rnd()
  iseed = mod(16807_i8*iseed,2147483647_i8)
  rnd = real(iseed,dp)/2147483647._dp
 end function rnd

end subroutine test_pipeline

end module testreplay
