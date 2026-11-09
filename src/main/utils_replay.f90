!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module replayutils
!
! Core routines for stitching an inner wind simulation to an outer one by
! periodic replay of particle tracks.
!
! The outer run reads one period of dumps of the inner run (records) and
! turns them into a track table: the path of every particle through the
! buffer shell r_b < r < r_inj, Hermite-interpolated between dumps, with
! exact entry and exit times. The table is replayed with period t_rec.
!
! A track starts when a particle crosses r_b outwards and ends when it
! crosses r_inj outwards (release: the particle becomes gas in the outer
! run) or r_b inwards (drop). With only one period of data the tracks are
! clipped to the window [t0,t0+t_rec]: a particle already in the buffer at
! t0 starts there (it appears at the start of every cycle) and a particle
! still in the buffer at t0+t_rec ends there (vanish). The replayed buffer
! then matches the inner one at every phase, and every crossing of r_inj
! in the window is replayed once per cycle.
!
! This module depends on no other Phantom module, so that it can also be
! used by the test suite.
!
! :References: None
!
! :Owner: Owen Vermeulen
!
! :Runtime parameters: None
!
! :Dependencies: None
!
 implicit none

 integer, parameter, public :: sp = selected_real_kind(6,37)
 integer, parameter, public :: dp = selected_real_kind(15,307)
 integer, parameter, public :: i4 = selected_int_kind(9)
 integer, parameter, public :: i8 = selected_int_kind(18)

 integer, parameter, public :: lenlabel = 16
 ! columns of the sink data stored in the log and the table
 integer, parameter, public :: nsinkcol = 11
 integer, parameter, public :: isk_m = 7, isk_lum = 8, isk_teff = 9, isk_reff = 10, isk_hacc = 11
 ! how a track ends
 integer, parameter, public :: iend_release = 1, iend_drop = 2, iend_vanish = 3
 ! per-cycle transformations of the replayed tracks
 integer, parameter, public :: irot_none = 0, irot_random = 1, irot_mirror = 2

 real(dp), parameter :: pi = 3.141592653589793_dp
 real(dp), parameter :: solarm_cgs = 1.9891e33_dp, year_cgs = 3.1556926e7_dp
 real(dp), parameter :: km_cgs = 1.e5_dp, au_cgs = 1.496e13_dp

 character(len=16), parameter :: tabmagic = 'PHANTOMREPLAYTAB'
 integer(i4),       parameter :: tabversion = 1

 !
 ! one record (dump): the sinks and the gas particles in the band
 ! (positions and velocities relative to the centre of mass of the sinks)
 !
 type, public :: replay_record
    real(dp) :: time = 0._dp
    real(dp) :: xc(3) = 0._dp, vc(3) = 0._dp
    integer  :: nsink = 0, nsel = 0
    real(dp),    allocatable :: sink(:,:)
    integer(i8), allocatable :: iorig(:)
    real(sp),    allocatable :: x(:,:),v(:,:),u(:),h(:),scal(:,:)
 end type replay_record

 !
 ! the track table
 !
 type, public :: track_table
    integer  :: nscal = 0, nsink = 0, ntracks = 0, nsamp = 0, nsinkrec = 0
    real(dp) :: udist = 1._dp, umass = 1._dp, utime = 1._dp, pmass = 0._dp
    real(dp) :: r_inj = 0._dp, r_b = 0._dp, t0 = 0._dp, t_rec = 0._dp, period = 0._dp
    character(len=lenlabel), allocatable :: label(:)
    ! per track
    real(dp), allocatable :: t_entry(:),t_exit(:)
    integer,  allocatable :: iend(:),isamp1(:),nsampk(:)
    ! per sample
    real(dp), allocatable :: ts(:)
    real(sp), allocatable :: xs(:,:),vs(:,:),us(:),hs(:),scal(:,:)
    ! sink history over the window
    real(dp), allocatable :: tsink(:),sink(:,:,:)
 end type track_table

 !
 ! spawn order of the replay instances
 !
 type, public :: replay_schedule
    real(dp) :: t_epoch = 0._dp, t0 = 0._dp, t_rec = 1._dp
    integer  :: ntr = 0, jnext = 1, mcycle = 0
    integer,  allocatable :: order(:),coff(:)
    real(dp), allocatable :: phase(:)
 end type replay_schedule

 !
 ! track being followed by the builder
 !
 type :: open_track
    integer  :: n = 0
    logical  :: used = .false., unknown = .false.
    real(dp) :: t_entry = 0._dp
    real(dp), allocatable :: t(:)
    real(sp), allocatable :: x(:,:),v(:,:),u(:),h(:),s(:,:)
 end type open_track

 type, public :: replay_builder
    real(dp) :: r_inj = 0._dp, r_b = 0._dp, t0 = 0._dp, t_rec = 0._dp, dt_table = 0._dp
    real(dp) :: tsink_lo = 0._dp, tsink_hi = 0._dp, tprev = 0._dp
    integer  :: nscal = 0, nprev = 0, nfree = 0
    logical  :: include_drop = .true., clip = .false.
    integer(i8), allocatable :: piorig(:)
    real(sp),    allocatable :: px(:,:),pv(:,:),pu(:),ph(:),ps(:,:)
    integer,     allocatable :: ptr(:),freelist(:)
    type(open_track), allocatable :: otr(:)
    type(track_table) :: tab
    ! statistics
    integer  :: nrel = 0, ndrop = 0, nunknown = 0, nlost = 0, nskipdrop = 0, ndip = 0, nopen = 0
    integer  :: nappear = 0, nvanish = 0
    real(dp) :: tres_min = huge(1._dp), tres_max = 0._dp, tres_sum = 0._dp, amin = huge(1._dp)
 end type replay_builder

 public :: hermite_eval,find_sphere_crossing,shell_flux
 public :: sort_index_i8,sort_index_dp,percentile,find_label
 public :: alloc_record,record_stats
 public :: builder_init,builder_add_record,builder_finalise
 public :: write_track_table,read_track_table,read_track_table_header
 public :: interp_track
 public :: schedule_init,schedule_reset,schedule_next_spawn,active_instances_at
 public :: cycle_rotation

 private

contains

!-----------------------------------------------------------------------
!+
!  Cubic Hermite interpolation of a path between two samples (x,v at ta
!  and tb). The velocity is the exact time derivative of the position.
!+
!-----------------------------------------------------------------------
pure subroutine hermite_eval(xa,va,xb,vb,ta,tb,t,x,v)
 real(dp), intent(in)  :: xa(3),va(3),xb(3),vb(3),ta,tb,t
 real(dp), intent(out) :: x(3),v(3)
 real(dp) :: dt,s,s2,s3,h00,h10,h01,h11,d00,d10,d11

 dt = tb - ta
 if (dt <= 0._dp) then
    x = xa + va*(t - ta)
    v = va
    return
 endif
 s   = (t - ta)/dt
 s2  = s*s
 s3  = s2*s
 h00 = 2._dp*s3 - 3._dp*s2 + 1._dp
 h10 = s3 - 2._dp*s2 + s
 h01 = -2._dp*s3 + 3._dp*s2
 h11 = s3 - s2
 d00 = 6._dp*s2 - 6._dp*s
 d10 = 3._dp*s2 - 4._dp*s + 1._dp
 d11 = 3._dp*s2 - 2._dp*s
 x = h00*xa + h10*dt*va + h01*xb + h11*dt*vb
 v = d00*(xa - xb)/dt + d10*va + d11*vb

end subroutine hermite_eval

!-----------------------------------------------------------------------
!+
!  Time at which the Hermite path between two samples crosses the sphere
!  of radius r, by bisection on g = |x|^2 - r^2. found is false if g has
!  the same sign at both ends.
!+
!-----------------------------------------------------------------------
subroutine find_sphere_crossing(xa,va,xb,vb,ta,tb,r,tc,found)
 real(dp), intent(in)  :: xa(3),va(3),xb(3),vb(3),ta,tb,r
 real(dp), intent(out) :: tc
 logical,  intent(out) :: found
 integer,  parameter   :: maxiter = 200
 real(dp) :: glo,ghi,gm,tlo,thi,tm,x(3),v(3)
 integer  :: iter

 glo = dot_product(xa,xa) - r*r
 ghi = dot_product(xb,xb) - r*r
 tc  = ta
 found = .false.
 if (glo > 0._dp .and. ghi > 0._dp) return
 if (glo < 0._dp .and. ghi < 0._dp) return
 found = .true.
 if (abs(glo) <= tiny(glo)) return
 if (abs(ghi) <= tiny(ghi)) then
    tc = tb
    return
 endif
 tlo = ta
 thi = tb
 do iter = 1,maxiter
    tm = 0.5_dp*(tlo + thi)
    call hermite_eval(xa,va,xb,vb,ta,tb,tm,x,v)
    gm = dot_product(x,x) - r*r
    if (gm*glo < 0._dp) then
       thi = tm
    else
       tlo = tm
       glo = gm
    endif
    if (thi - tlo <= 1.e-14_dp*(tb - ta)) exit
 enddo
 tc = 0.5_dp*(tlo + thi)

end subroutine find_sphere_crossing

!-----------------------------------------------------------------------
!+
!  State at time t between two samples: Hermite for x and v, linear for
!  u and the scalars, log-linear for h (scalars are not extrapolated)
!+
!-----------------------------------------------------------------------
subroutine interp_state(ta,xa,va,ua,ha,sa,tb,xb,vb,ub,hb,sb,t,x,v,u,h,s)
 real(dp), intent(in)  :: ta,tb,t
 real(sp), intent(in)  :: xa(3),va(3),ua,ha,sa(:),xb(3),vb(3),ub,hb,sb(:)
 real(dp), intent(out) :: x(3),v(3),u,h,s(:)
 real(dp) :: f

 call hermite_eval(real(xa,dp),real(va,dp),real(xb,dp),real(vb,dp),ta,tb,t,x,v)
 f = 0._dp
 if (tb > ta) f = min(max((t - ta)/(tb - ta),0._dp),1._dp)
 u = (1._dp - f)*ua + f*ub
 if (ha > 0. .and. hb > 0.) then
    h = exp((1._dp - f)*log(real(ha,dp)) + f*log(real(hb,dp)))
 else
    h = (1._dp - f)*ha + f*hb
 endif
 if (size(s) > 0) s = (1._dp - f)*sa + f*sb

end subroutine interp_state

!-----------------------------------------------------------------------
!+
!  Mass flux through a sphere of radius rshell from the particles within
!  dr/2 of it, Mdot = (1/dr) sum m v_r, separately for outflow and inflow
!  (both positive). vr_flux is the flux-weighted outflow speed and wflux
!  the outward flux of the quantity w (per unit mass).
!+
!-----------------------------------------------------------------------
pure subroutine shell_flux(n,r,vr,pmass,rshell,dr,mdot_out,mdot_in,vr_flux,w,wflux)
 integer,  intent(in)  :: n
 real(dp), intent(in)  :: r(:),vr(:),pmass,rshell,dr
 real(dp), intent(out) :: mdot_out,mdot_in,vr_flux
 real(dp), intent(in),  optional :: w(:)
 real(dp), intent(out), optional :: wflux
 real(dp) :: sum_vr,sum_vr2,sum_in,sum_w
 integer  :: i

 sum_vr = 0._dp; sum_vr2 = 0._dp; sum_in = 0._dp; sum_w = 0._dp
 do i = 1,n
    if (abs(r(i) - rshell) >= 0.5_dp*dr) cycle
    if (vr(i) > 0._dp) then
       sum_vr  = sum_vr  + vr(i)
       sum_vr2 = sum_vr2 + vr(i)**2
       if (present(w)) sum_w = sum_w + vr(i)*w(i)
    else
       sum_in  = sum_in  - vr(i)
    endif
 enddo
 mdot_out = 0._dp; mdot_in = 0._dp; vr_flux = 0._dp
 if (dr > 0._dp) then
    mdot_out = pmass*sum_vr/dr
    mdot_in  = pmass*sum_in/dr
 endif
 if (sum_vr > 0._dp) vr_flux = sum_vr2/sum_vr
 if (present(wflux)) then
    wflux = 0._dp
    if (dr > 0._dp) wflux = pmass*sum_w/dr
 endif

end subroutine shell_flux

!-----------------------------------------------------------------------
!+
!  Heapsort index of 8-byte integer keys (ascending)
!+
!-----------------------------------------------------------------------
subroutine sort_index_i8(n,key,idx)
 integer,     intent(in)  :: n
 integer(i8), intent(in)  :: key(:)
 integer,     intent(out) :: idx(:)
 integer     :: i,j,l,ir,it
 integer(i8) :: kt

 do i = 1,n
    idx(i) = i
 enddo
 if (n < 2) return
 l  = n/2 + 1
 ir = n
 do
    if (l > 1) then
       l  = l - 1
       it = idx(l)
       kt = key(it)
    else
       it = idx(ir)
       kt = key(it)
       idx(ir) = idx(1)
       ir = ir - 1
       if (ir == 1) then
          idx(1) = it
          exit
       endif
    endif
    i = l
    j = l + l
    do while (j <= ir)
       if (j < ir) then
          if (key(idx(j)) < key(idx(j+1))) j = j + 1
       endif
       if (kt < key(idx(j))) then
          idx(i) = idx(j)
          i = j
          j = j + j
       else
          j = ir + 1
       endif
    enddo
    idx(i) = it
 enddo

end subroutine sort_index_i8

!-----------------------------------------------------------------------
!+
!  Heapsort index of real keys (ascending)
!+
!-----------------------------------------------------------------------
subroutine sort_index_dp(n,key,idx)
 integer,  intent(in)  :: n
 real(dp), intent(in)  :: key(:)
 integer,  intent(out) :: idx(:)
 integer  :: i,j,l,ir,it
 real(dp) :: kt

 do i = 1,n
    idx(i) = i
 enddo
 if (n < 2) return
 l  = n/2 + 1
 ir = n
 do
    if (l > 1) then
       l  = l - 1
       it = idx(l)
       kt = key(it)
    else
       it = idx(ir)
       kt = key(it)
       idx(ir) = idx(1)
       ir = ir - 1
       if (ir == 1) then
          idx(1) = it
          exit
       endif
    endif
    i = l
    j = l + l
    do while (j <= ir)
       if (j < ir) then
          if (key(idx(j)) < key(idx(j+1))) j = j + 1
       endif
       if (kt < key(idx(j))) then
          idx(i) = idx(j)
          i = j
          j = j + j
       else
          j = ir + 1
       endif
    enddo
    idx(i) = it
 enddo

end subroutine sort_index_dp

!-----------------------------------------------------------------------
!+
!  The q-th quantile (0 <= q <= 1) of the first n values of a
!+
!-----------------------------------------------------------------------
real(dp) function percentile(n,a,q)
 integer,  intent(in) :: n
 real(dp), intent(in) :: a(:),q
 integer, allocatable :: idx(:)
 integer :: j

 percentile = 0._dp
 if (n < 1) return
 allocate(idx(n))
 call sort_index_dp(n,a,idx)
 j = min(max(1,nint(q*(n-1)) + 1),n)
 percentile = a(idx(j))

end function percentile

!-----------------------------------------------------------------------
!+
!  Index of a label in a list (0 if absent)
!+
!-----------------------------------------------------------------------
integer function find_label(label,name)
 character(len=*), intent(in) :: label(:),name
 integer :: i

 find_label = 0
 do i = 1,size(label)
    if (trim(adjustl(label(i))) == trim(adjustl(name))) then
       find_label = i
       return
    endif
 enddo

end function find_label

!-----------------------------------------------------------------------
!+
!  Make sure a record can hold nsink sinks and nsel particles
!+
!-----------------------------------------------------------------------
subroutine alloc_record(rec,nsink,nsel,nscal)
 type(replay_record), intent(inout) :: rec
 integer,             intent(in)    :: nsink,nsel,nscal
 integer :: ncap

 if (allocated(rec%sink)) then
    if (size(rec%sink,2) < nsink) deallocate(rec%sink)
 endif
 if (.not.allocated(rec%sink)) allocate(rec%sink(nsinkcol,max(nsink,1)))

 if (allocated(rec%iorig)) then
    if (size(rec%iorig) >= nsel .and. size(rec%scal,1) == max(nscal,1)) return
    deallocate(rec%iorig,rec%x,rec%v,rec%u,rec%h,rec%scal)
 endif
 ncap = max(nsel + nsel/5,16)
 allocate(rec%iorig(ncap),rec%x(3,ncap),rec%v(3,ncap),rec%u(ncap),rec%h(ncap),rec%scal(max(nscal,1),ncap))

end subroutine alloc_record

!-----------------------------------------------------------------------
!+
!  Set up the builder. head carries the header of the table to build:
!  units, particle mass, r_inj, r_b, window [t0,t0+t_rec], period and the
!  labels of the scalars. With clip, the tracks are clipped to the window
!  (see the module description).
!+
!-----------------------------------------------------------------------
subroutine builder_init(b,head,dt_table,include_drop,clip)
 type(replay_builder), intent(inout) :: b
 type(track_table),    intent(in)    :: head
 real(dp),             intent(in)    :: dt_table
 logical,              intent(in)    :: include_drop,clip

 b%r_inj        = head%r_inj
 b%r_b          = head%r_b
 b%t0           = head%t0
 b%t_rec        = head%t_rec
 b%dt_table     = dt_table
 b%include_drop = include_drop
 b%clip         = clip
 b%nscal        = head%nscal
 b%tsink_lo     = head%t0 - 0.05_dp*head%t_rec
 b%tsink_hi     = head%t0 + 1.05_dp*head%t_rec
 b%nprev        = 0

 b%tab%udist  = head%udist
 b%tab%umass  = head%umass
 b%tab%utime  = head%utime
 b%tab%pmass  = head%pmass
 b%tab%r_inj  = head%r_inj
 b%tab%r_b    = head%r_b
 b%tab%t0     = head%t0
 b%tab%t_rec  = head%t_rec
 b%tab%period = head%period
 b%tab%nscal  = head%nscal
 if (allocated(b%tab%label)) deallocate(b%tab%label)
 allocate(b%tab%label(max(head%nscal,0)))
 if (head%nscal > 0) b%tab%label = head%label(1:head%nscal)
 b%tab%ntracks  = 0
 b%tab%nsamp    = 0
 b%tab%nsinkrec = 0
 b%tab%nsink    = 0

end subroutine builder_init

!-----------------------------------------------------------------------
!+
!  Add one record to the builder. Particles are matched with the previous
!  record by iorig; for each particle the path between the two records is
!  checked for crossings of r_b and r_inj.
!+
!-----------------------------------------------------------------------
subroutine builder_add_record(b,rec)
 type(replay_builder), intent(inout) :: b
 type(replay_record),  intent(in)    :: rec
 integer,     allocatable :: idx(:),ctr(:)
 integer(i8), allocatable :: ciorig(:)
 real(sp),    allocatable :: cx(:,:),cv(:,:),cu(:),ch(:),cs(:,:)
 integer :: n,i,j,ns,itr

 n  = rec%nsel
 ns = max(b%nscal,1)
 allocate(idx(n),ctr(n),ciorig(n),cx(3,n),cv(3,n),cu(n),ch(n),cs(ns,n))
 call sort_index_i8(n,rec%iorig,idx)
 do j = 1,n
    ciorig(j) = rec%iorig(idx(j))
    cx(:,j)   = rec%x(:,idx(j))
    cv(:,j)   = rec%v(:,idx(j))
    cu(j)     = rec%u(idx(j))
    ch(j)     = rec%h(idx(j))
    cs(:,j)   = rec%scal(1:ns,idx(j))
 enddo
 ctr = 0

 i = 1
 j = 1
 do while (i <= b%nprev .or. j <= n)
    if (i <= b%nprev .and. j <= n) then
       if (b%piorig(i) == ciorig(j)) then
          itr = b%ptr(i)
          call process_interval(b,itr,b%tprev,b%px(:,i),b%pv(:,i),b%pu(i),b%ph(i),b%ps(:,i), &
                                rec%time,cx(:,j),cv(:,j),cu(j),ch(j),cs(:,j))
          ctr(j) = itr
          i = i + 1
          j = j + 1
       elseif (b%piorig(i) < ciorig(j)) then
          call left_band(i)
          i = i + 1
       else
          call entered_band(j)
          j = j + 1
       endif
    elseif (i <= b%nprev) then
       call left_band(i)
       i = i + 1
    else
       call entered_band(j)
       j = j + 1
    endif
 enddo

 call store_sinks(b,rec)

 call move_alloc(ciorig,b%piorig)
 call move_alloc(cx,b%px)
 call move_alloc(cv,b%pv)
 call move_alloc(cu,b%pu)
 call move_alloc(ch,b%ph)
 call move_alloc(cs,b%ps)
 call move_alloc(ctr,b%ptr)
 b%nprev = n
 b%tprev = rec%time

contains

 ! a particle that was followed has left the band: should not happen for a track
 subroutine left_band(ip)
  integer, intent(in) :: ip

  if (b%ptr(ip) > 0) then
     b%nlost = b%nlost + 1
     call free_track(b,b%ptr(ip))
  endif

 end subroutine left_band

 ! a particle seen for the first time; if it is already inside the buffer its
 ! entry was not recorded (first record, or a gap in the log)
 subroutine entered_band(jp)
  integer, intent(in) :: jp
  real(dp) :: r
  integer  :: k

  r = norm2(real(cx(:,jp),dp))
  if (r > b%r_b .and. r < b%r_inj) then
     k = new_track(b)
     b%otr(k)%unknown = .true.
     b%otr(k)%t_entry = rec%time
     call add_sample(b%otr(k),rec%time,cx(:,jp),cv(:,jp),cu(jp),ch(jp),cs(:,jp))
     ctr(jp) = k
  endif

 end subroutine entered_band

end subroutine builder_add_record

!-----------------------------------------------------------------------
!+
!  Path of one particle between two records (a,b): open, extend or close
!  its track
!+
!-----------------------------------------------------------------------
subroutine process_interval(b,itr,ta,xa,va,ua,ha,sa,tb,xb,vb,ub,hb,sb)
 type(replay_builder), intent(inout) :: b
 integer,              intent(inout) :: itr
 real(dp),             intent(in)    :: ta,tb
 real(sp),             intent(in)    :: xa(3),va(3),ua,ha,sa(:),xb(3),vb(3),ub,hb,sb(:)
 real(dp) :: ra,rb,te,tx
 logical  :: entered

 ra = norm2(real(xa,dp))
 rb = norm2(real(xb,dp))
 entered = .false.
 te = ta
 if (itr == 0) then
    if (ra < b%r_b .and. rb >= b%r_b) then
       te  = crossing_time(b%r_b)
       itr = new_track(b)
       b%otr(itr)%t_entry = te
       call add_sample(b%otr(itr),ta,xa,va,ua,ha,sa)
       call add_interp(te)
       entered = .true.
    else
       return
    endif
 endif

 if (rb >= b%r_inj) then
    tx = max(crossing_time(b%r_inj),te)
    call add_interp(tx)
    call add_sample(b%otr(itr),tb,xb,vb,ub,hb,sb)
    call close_track(b,itr,tx,iend_release)
 elseif (rb < b%r_b) then
    if (entered) then
       ! entered and left within one interval
       b%ndip = b%ndip + 1
       call free_track(b,itr)
       itr = 0
    else
       tx = crossing_time(b%r_b)
       call add_interp(tx)
       call add_sample(b%otr(itr),tb,xb,vb,ub,hb,sb)
       call close_track(b,itr,tx,iend_drop)
    endif
 else
    call add_sample(b%otr(itr),tb,xb,vb,ub,hb,sb)
 endif

contains

 real(dp) function crossing_time(r)
  real(dp), intent(in) :: r
  logical  :: found
  real(dp) :: f

  call find_sphere_crossing(real(xa,dp),real(va,dp),real(xb,dp),real(vb,dp),ta,tb,r,crossing_time,found)
  if (.not.found) then
     f = 0._dp
     if (abs(rb - ra) > 0._dp) f = min(max((r - ra)/(rb - ra),0._dp),1._dp)
     crossing_time = ta + f*(tb - ta)
  endif

 end function crossing_time

 subroutine add_interp(t)
  real(dp), intent(in) :: t
  real(dp) :: x(3),v(3),u,h,s(size(sa))

  call interp_state(ta,xa,va,ua,ha,sa,tb,xb,vb,ub,hb,sb,t,x,v,u,h,s)
  call add_sample(b%otr(itr),t,real(x,sp),real(v,sp),real(u,sp),real(h,sp),real(s,sp))

 end subroutine add_interp

end subroutine process_interval

!-----------------------------------------------------------------------
!+
!  Track bookkeeping
!+
!-----------------------------------------------------------------------
integer function new_track(b) result(itr)
 type(replay_builder), intent(inout) :: b
 type(open_track), allocatable :: tmp(:)
 integer, allocatable :: itmp(:)
 integer :: ncap,k

 if (.not.allocated(b%otr)) then
    allocate(b%otr(1024),b%freelist(1024))
    b%nfree = 1024
    do k = 1,1024
       b%freelist(k) = 1025 - k
    enddo
 endif
 if (b%nfree == 0) then
    ncap = size(b%otr)
    allocate(tmp(2*ncap))
    tmp(1:ncap) = b%otr(1:ncap)
    call move_alloc(tmp,b%otr)
    allocate(itmp(2*ncap))
    call move_alloc(itmp,b%freelist)
    do k = 1,ncap
       b%freelist(k) = 2*ncap + 1 - k
    enddo
    b%nfree = ncap
 endif
 itr = b%freelist(b%nfree)
 b%nfree = b%nfree - 1
 b%otr(itr)%used    = .true.
 b%otr(itr)%unknown = .false.
 b%otr(itr)%n       = 0
 b%otr(itr)%t_entry = 0._dp

end function new_track

subroutine free_track(b,itr)
 type(replay_builder), intent(inout) :: b
 integer,              intent(in)    :: itr

 if (.not.b%otr(itr)%used) return
 b%otr(itr)%used = .false.
 b%otr(itr)%n    = 0
 b%nfree = b%nfree + 1
 b%freelist(b%nfree) = itr

end subroutine free_track

subroutine add_sample(ot,t,x,v,u,h,s)
 type(open_track), intent(inout) :: ot
 real(dp),         intent(in)    :: t
 real(sp),         intent(in)    :: x(3),v(3),u,h,s(:)
 real(dp), allocatable :: rtmp(:)
 real(sp), allocatable :: tmp1(:),tmp2(:,:)
 integer :: ncap,ns

 ns = size(s)
 if (.not.allocated(ot%t)) then
    allocate(ot%t(32),ot%x(3,32),ot%v(3,32),ot%u(32),ot%h(32),ot%s(max(ns,1),32))
 endif
 if (ot%n >= size(ot%t)) then
    ncap = 2*size(ot%t)
    allocate(rtmp(ncap)); rtmp(1:ot%n) = ot%t(1:ot%n); call move_alloc(rtmp,ot%t)
    allocate(tmp2(3,ncap)); tmp2(:,1:ot%n) = ot%x(:,1:ot%n); call move_alloc(tmp2,ot%x)
    allocate(tmp2(3,ncap)); tmp2(:,1:ot%n) = ot%v(:,1:ot%n); call move_alloc(tmp2,ot%v)
    allocate(tmp1(ncap)); tmp1(1:ot%n) = ot%u(1:ot%n); call move_alloc(tmp1,ot%u)
    allocate(tmp1(ncap)); tmp1(1:ot%n) = ot%h(1:ot%n); call move_alloc(tmp1,ot%h)
    allocate(tmp2(size(ot%s,1),ncap)); tmp2(:,1:ot%n) = ot%s(:,1:ot%n); call move_alloc(tmp2,ot%s)
 endif
 ot%n = ot%n + 1
 ot%t(ot%n)   = t
 ot%x(:,ot%n) = x
 ot%v(:,ot%n) = v
 ot%u(ot%n)   = u
 ot%h(ot%n)   = h
 if (ns > 0) ot%s(1:ns,ot%n) = s

end subroutine add_sample

!-----------------------------------------------------------------------
!+
!  A track has ended at time tx: keep it if it ends inside the window
!+
!-----------------------------------------------------------------------
subroutine close_track(b,itr,tx,iend)
 type(replay_builder), intent(inout) :: b
 integer,              intent(inout) :: itr
 real(dp),             intent(in)    :: tx
 integer,              intent(in)    :: iend
 real(dp) :: tres

 if (b%clip .and. tx >= b%t0 + b%t_rec .and. b%otr(itr)%t_entry < b%t0 + b%t_rec) then
    ! left the buffer just after the window: the piece ends at the window end
    call keep_track(b%t0 + b%t_rec,iend_vanish)
 elseif (tx >= b%t0 .and. tx < b%t0 + b%t_rec) then
    if (b%otr(itr)%unknown .and. .not.b%clip) then
       if (iend == iend_release) b%nunknown = b%nunknown + 1
    elseif (iend == iend_drop .and. .not.b%include_drop) then
       b%nskipdrop = b%nskipdrop + 1
    else
       call keep_track(tx,iend)
    endif
 endif
 call free_track(b,itr)
 itr = 0

contains

 subroutine keep_track(tend,iendk)
  real(dp), intent(in) :: tend
  integer,  intent(in) :: iendk

  call append_track(b,b%otr(itr),tend,iendk)
  select case(iendk)
  case(iend_release)
     b%nrel = b%nrel + 1
  case(iend_drop)
     b%ndrop = b%ndrop + 1
  case default
     b%nvanish = b%nvanish + 1
  end select
  if (b%otr(itr)%unknown) b%nappear = b%nappear + 1
  tres = tend - b%otr(itr)%t_entry
  b%tres_min = min(b%tres_min,tres)
  b%tres_max = max(b%tres_max,tres)
  b%tres_sum = b%tres_sum + tres
  b%amin     = min(b%amin,b%otr(itr)%t_entry - b%t0)

 end subroutine keep_track

end subroutine close_track

!-----------------------------------------------------------------------
!+
!  Append a finished track to the table, keeping the sample before entry,
!  the entry and exit states, the sample after exit and interior samples
!  spaced by at least dt_table
!+
!-----------------------------------------------------------------------
subroutine append_track(b,ot,tx,iend)
 type(replay_builder), intent(inout) :: b
 type(open_track),     intent(in)    :: ot
 real(dp),             intent(in)    :: tx
 integer,              intent(in)    :: iend
 logical, allocatable :: keep(:)
 integer  :: n,j,nkeep,ns,i0
 real(dp) :: tl

 n  = ot%n
 ns = b%nscal
 allocate(keep(n))
 if (n <= 4 .or. b%dt_table <= 0._dp) then
    keep = .true.
 else
    keep = .false.
    keep(1:2) = .true.
    keep(n-1:n) = .true.
    tl = ot%t(2)
    do j = 3,n-2
       if (ot%t(j) - tl >= (1._dp - 1.e-6_dp)*b%dt_table) then
          keep(j) = .true.
          tl = ot%t(j)
       endif
    enddo
 endif
 nkeep = count(keep)
 call ensure_table_capacity(b%tab,b%tab%ntracks+1,b%tab%nsamp+nkeep)

 b%tab%ntracks = b%tab%ntracks + 1
 b%tab%t_entry(b%tab%ntracks) = ot%t_entry
 b%tab%t_exit(b%tab%ntracks)  = tx
 b%tab%iend(b%tab%ntracks)    = iend
 b%tab%isamp1(b%tab%ntracks)  = b%tab%nsamp + 1
 b%tab%nsampk(b%tab%ntracks)  = nkeep
 i0 = b%tab%nsamp
 do j = 1,n
    if (.not.keep(j)) cycle
    i0 = i0 + 1
    b%tab%ts(i0)   = ot%t(j)
    b%tab%xs(:,i0) = ot%x(:,j)
    b%tab%vs(:,i0) = ot%v(:,j)
    b%tab%us(i0)   = ot%u(j)
    b%tab%hs(i0)   = ot%h(j)
    if (ns > 0) b%tab%scal(1:ns,i0) = ot%s(1:ns,j)
 enddo
 b%tab%nsamp = i0

end subroutine append_track

subroutine ensure_table_capacity(tab,ntr,nsamp)
 type(track_table), intent(inout) :: tab
 integer,           intent(in)    :: ntr,nsamp
 real(dp), allocatable :: rtmp(:)
 integer,  allocatable :: itmp(:)
 real(sp), allocatable :: tmp1(:),tmp2(:,:)
 integer :: ncap,nold,ns

 ns = max(tab%nscal,1)
 if (.not.allocated(tab%t_entry)) then
    ncap = max(ntr,4096)
    allocate(tab%t_entry(ncap),tab%t_exit(ncap),tab%iend(ncap),tab%isamp1(ncap),tab%nsampk(ncap))
 elseif (ntr > size(tab%t_entry)) then
    ncap = max(ntr,2*size(tab%t_entry))
    nold = tab%ntracks
    allocate(rtmp(ncap)); rtmp(1:nold) = tab%t_entry(1:nold); call move_alloc(rtmp,tab%t_entry)
    allocate(rtmp(ncap)); rtmp(1:nold) = tab%t_exit(1:nold);  call move_alloc(rtmp,tab%t_exit)
    allocate(itmp(ncap)); itmp(1:nold) = tab%iend(1:nold);    call move_alloc(itmp,tab%iend)
    allocate(itmp(ncap)); itmp(1:nold) = tab%isamp1(1:nold);  call move_alloc(itmp,tab%isamp1)
    allocate(itmp(ncap)); itmp(1:nold) = tab%nsampk(1:nold);  call move_alloc(itmp,tab%nsampk)
 endif

 if (.not.allocated(tab%ts)) then
    ncap = max(nsamp,65536)
    allocate(tab%ts(ncap),tab%xs(3,ncap),tab%vs(3,ncap),tab%us(ncap),tab%hs(ncap),tab%scal(ns,ncap))
 elseif (nsamp > size(tab%ts)) then
    ncap = max(nsamp,2*size(tab%ts))
    nold = tab%nsamp
    allocate(rtmp(ncap));    rtmp(1:nold)   = tab%ts(1:nold);     call move_alloc(rtmp,tab%ts)
    allocate(tmp2(3,ncap));  tmp2(:,1:nold) = tab%xs(:,1:nold);   call move_alloc(tmp2,tab%xs)
    allocate(tmp2(3,ncap));  tmp2(:,1:nold) = tab%vs(:,1:nold);   call move_alloc(tmp2,tab%vs)
    allocate(tmp1(ncap));    tmp1(1:nold)   = tab%us(1:nold);     call move_alloc(tmp1,tab%us)
    allocate(tmp1(ncap));    tmp1(1:nold)   = tab%hs(1:nold);     call move_alloc(tmp1,tab%hs)
    allocate(tmp2(ns,ncap)); tmp2(:,1:nold) = tab%scal(:,1:nold); call move_alloc(tmp2,tab%scal)
 endif

end subroutine ensure_table_capacity

!-----------------------------------------------------------------------
!+
!  Keep the sink history over the window (for the luminosity replay)
!+
!-----------------------------------------------------------------------
subroutine store_sinks(b,rec)
 type(replay_builder), intent(inout) :: b
 type(replay_record),  intent(in)    :: rec
 real(dp), allocatable :: rtmp(:),stmp(:,:,:)
 integer :: ncap,n

 if (rec%nsink < 1) return
 if (rec%time < b%tsink_lo .or. rec%time > b%tsink_hi) return
 if (b%tab%nsinkrec == 0) then
    b%tab%nsink = rec%nsink
    if (allocated(b%tab%tsink)) deallocate(b%tab%tsink,b%tab%sink)
    allocate(b%tab%tsink(256),b%tab%sink(nsinkcol,rec%nsink,256))
 elseif (rec%nsink /= b%tab%nsink) then
    return
 endif
 n = b%tab%nsinkrec
 if (n >= size(b%tab%tsink)) then
    ncap = 2*size(b%tab%tsink)
    allocate(rtmp(ncap)); rtmp(1:n) = b%tab%tsink(1:n); call move_alloc(rtmp,b%tab%tsink)
    allocate(stmp(nsinkcol,b%tab%nsink,ncap)); stmp(:,:,1:n) = b%tab%sink(:,:,1:n)
    call move_alloc(stmp,b%tab%sink)
 endif
 n = n + 1
 b%tab%tsink(n)    = rec%time
 b%tab%sink(:,:,n) = rec%sink(:,1:b%tab%nsink)
 b%tab%nsinkrec    = n

end subroutine store_sinks

!-----------------------------------------------------------------------
!+
!  Finish: tracks still open could end inside the window if the recording
!  stopped too early, so they are counted
!+
!-----------------------------------------------------------------------
subroutine builder_finalise(b,tab)
 type(replay_builder), intent(inout) :: b
 type(track_table),    intent(out)   :: tab
 integer :: k

 b%nopen = 0
 if (allocated(b%otr)) then
    do k = 1,size(b%otr)
       if (b%otr(k)%used) then
          if (b%clip .and. b%otr(k)%t_entry < b%t0 + b%t_rec) then
             call append_track(b,b%otr(k),b%t0 + b%t_rec,iend_vanish)
             b%nvanish = b%nvanish + 1
             if (b%otr(k)%unknown) b%nappear = b%nappear + 1
          elseif (.not.b%otr(k)%unknown) then
             b%nopen = b%nopen + 1
          endif
          call free_track(b,k)
       endif
    enddo
 endif
 call ensure_table_capacity(b%tab,max(b%tab%ntracks,1),max(b%tab%nsamp,1))
 tab = b%tab

end subroutine builder_finalise

!-----------------------------------------------------------------------
!+
!  Mass flux through r_inj and smoothing length there, from one record:
!  h is averaged over the particles within hwin of r_inj, and the shell
!  estimator uses dr = mean h
!+
!-----------------------------------------------------------------------
subroutine record_stats(rec,pmass,r_inj,hwin,mout,mdin,vrf,h90,hmean)
 type(replay_record), intent(in)  :: rec
 real(dp),            intent(in)  :: pmass,r_inj,hwin
 real(dp),            intent(out) :: mout,mdin,vrf,h90,hmean
 real(dp), allocatable :: r(:),vr(:),hl(:)
 integer  :: i,n,nh
 real(dp) :: dr

 n = rec%nsel
 allocate(r(max(n,1)),vr(max(n,1)),hl(max(n,1)))
 nh = 0
 do i = 1,n
    r(i)  = norm2(real(rec%x(:,i),dp))
    vr(i) = 0._dp
    if (r(i) > 0._dp) vr(i) = dot_product(real(rec%x(:,i),dp),real(rec%v(:,i),dp))/r(i)
    if (abs(r(i) - r_inj) < hwin) then
       nh = nh + 1
       hl(nh) = rec%h(i)
    endif
 enddo
 if (nh >= 3) then
    hmean = sum(hl(1:nh))/nh
    h90   = percentile(nh,hl,0.9_dp)
    dr    = hmean
 else
    hmean = 0._dp
    h90   = 0._dp
    dr    = hwin
 endif
 call shell_flux(n,r,vr,pmass,r_inj,dr,mout,mdin,vrf)

end subroutine record_stats

!-----------------------------------------------------------------------
!+
!  Write / read the track table
!+
!-----------------------------------------------------------------------
subroutine write_track_table(fname,tab,ierr)
 character(len=*),  intent(in)  :: fname
 type(track_table), intent(in)  :: tab
 integer,           intent(out) :: ierr
 integer :: iunit,nt,ns,nsc

 open(newunit=iunit,file=fname,access='stream',form='unformatted',status='replace',action='write',iostat=ierr)
 if (ierr /= 0) return
 nt  = tab%ntracks
 ns  = tab%nsamp
 nsc = tab%nscal
 write(iunit) tabmagic,tabversion
 write(iunit) tab%udist,tab%umass,tab%utime,tab%pmass,tab%r_inj,tab%r_b,tab%t0,tab%t_rec,tab%period
 write(iunit) int(nsc,i4),int(tab%nsink,i4),int(nt,i4),int(ns,i4),int(tab%nsinkrec,i4)
 if (nsc > 0) write(iunit) tab%label(1:nsc)
 if (nt > 0) write(iunit) tab%t_entry(1:nt),tab%t_exit(1:nt),int(tab%iend(1:nt),i4), &
                          int(tab%isamp1(1:nt),i4),int(tab%nsampk(1:nt),i4)
 if (ns > 0) then
    write(iunit) tab%ts(1:ns),tab%xs(:,1:ns),tab%vs(:,1:ns),tab%us(1:ns),tab%hs(1:ns)
    if (nsc > 0) write(iunit) tab%scal(1:nsc,1:ns)
 endif
 if (tab%nsinkrec > 0 .and. tab%nsink > 0) write(iunit) tab%tsink(1:tab%nsinkrec), &
     tab%sink(:,1:tab%nsink,1:tab%nsinkrec)
 close(iunit)

end subroutine write_track_table

subroutine read_table_header_unit(iunit,tab,ierr)
 integer,           intent(in)    :: iunit
 type(track_table), intent(inout) :: tab
 integer,           intent(out)   :: ierr
 character(len=16) :: magic
 integer(i4) :: iversion,nsc,nsk,nt,ns,nsr

 read(iunit,iostat=ierr) magic,iversion
 if (ierr /= 0) return
 if (magic /= tabmagic .or. iversion /= tabversion) then
    ierr = 1
    return
 endif
 read(iunit,iostat=ierr) tab%udist,tab%umass,tab%utime,tab%pmass,tab%r_inj,tab%r_b,tab%t0,tab%t_rec,tab%period
 if (ierr /= 0) return
 read(iunit,iostat=ierr) nsc,nsk,nt,ns,nsr
 if (ierr /= 0) return
 tab%nscal    = nsc
 tab%nsink    = nsk
 tab%ntracks  = nt
 tab%nsamp    = ns
 tab%nsinkrec = nsr
 if (allocated(tab%label)) deallocate(tab%label)
 allocate(tab%label(max(nsc,0)))
 if (nsc > 0) read(iunit,iostat=ierr) tab%label(1:nsc)

end subroutine read_table_header_unit

subroutine read_track_table_header(fname,tab,ierr)
 character(len=*),  intent(in)  :: fname
 type(track_table), intent(out) :: tab
 integer,           intent(out) :: ierr
 integer :: iunit

 open(newunit=iunit,file=fname,access='stream',form='unformatted',status='old',action='read',iostat=ierr)
 if (ierr /= 0) return
 call read_table_header_unit(iunit,tab,ierr)
 close(iunit)

end subroutine read_track_table_header

subroutine read_track_table(fname,tab,ierr)
 character(len=*),  intent(in)  :: fname
 type(track_table), intent(out) :: tab
 integer,           intent(out) :: ierr
 integer(i4), allocatable :: itmp(:,:)
 integer :: iunit,nt,ns,nsc

 open(newunit=iunit,file=fname,access='stream',form='unformatted',status='old',action='read',iostat=ierr)
 if (ierr /= 0) return
 call read_table_header_unit(iunit,tab,ierr)
 if (ierr /= 0) then
    close(iunit)
    return
 endif
 nt  = tab%ntracks
 ns  = tab%nsamp
 nsc = tab%nscal
 allocate(tab%t_entry(max(nt,1)),tab%t_exit(max(nt,1)),tab%iend(max(nt,1)), &
          tab%isamp1(max(nt,1)),tab%nsampk(max(nt,1)),itmp(max(nt,1),3))
 allocate(tab%ts(max(ns,1)),tab%xs(3,max(ns,1)),tab%vs(3,max(ns,1)),tab%us(max(ns,1)), &
          tab%hs(max(ns,1)),tab%scal(max(nsc,1),max(ns,1)))
 if (nt > 0) then
    read(iunit,iostat=ierr) tab%t_entry(1:nt),tab%t_exit(1:nt),itmp(1:nt,1),itmp(1:nt,2),itmp(1:nt,3)
    if (ierr /= 0) goto 10
    tab%iend(1:nt)   = itmp(1:nt,1)
    tab%isamp1(1:nt) = itmp(1:nt,2)
    tab%nsampk(1:nt) = itmp(1:nt,3)
 endif
 if (ns > 0) then
    read(iunit,iostat=ierr) tab%ts(1:ns),tab%xs(:,1:ns),tab%vs(:,1:ns),tab%us(1:ns),tab%hs(1:ns)
    if (ierr /= 0) goto 10
    if (nsc > 0) read(iunit,iostat=ierr) tab%scal(1:nsc,1:ns)
    if (ierr /= 0) goto 10
 endif
 if (tab%nsinkrec > 0 .and. tab%nsink > 0) then
    allocate(tab%tsink(tab%nsinkrec),tab%sink(nsinkcol,tab%nsink,tab%nsinkrec))
    read(iunit,iostat=ierr) tab%tsink,tab%sink
 else
    tab%nsinkrec = 0
 endif
10 close(iunit)

end subroutine read_track_table

!-----------------------------------------------------------------------
!+
!  State of track k at recorded time tau. Outside the samples the
!  position is extrapolated linearly.
!+
!-----------------------------------------------------------------------
subroutine interp_track(tab,k,tau,x,v,u,h,s)
 type(track_table), intent(in)  :: tab
 integer,           intent(in)  :: k
 real(dp),          intent(in)  :: tau
 real(dp),          intent(out) :: x(3),v(3),u,h,s(:)
 integer :: i1,i2,lo,hi,mid,ns

 ns = tab%nscal
 i1 = tab%isamp1(k)
 i2 = i1 + tab%nsampk(k) - 1
 if (tau <= tab%ts(i1) .or. i2 == i1) then
    call sample_state(i1,tau - tab%ts(i1))
 elseif (tau >= tab%ts(i2)) then
    call sample_state(i2,tau - tab%ts(i2))
 else
    lo = i1
    hi = i2
    do while (hi - lo > 1)
       mid = (lo + hi)/2
       if (tab%ts(mid) <= tau) then
          lo = mid
       else
          hi = mid
       endif
    enddo
    call interp_state(tab%ts(lo),tab%xs(:,lo),tab%vs(:,lo),tab%us(lo),tab%hs(lo),tab%scal(1:ns,lo), &
                      tab%ts(hi),tab%xs(:,hi),tab%vs(:,hi),tab%us(hi),tab%hs(hi),tab%scal(1:ns,hi), &
                      tau,x,v,u,h,s(1:ns))
 endif

contains

 subroutine sample_state(i,dt)
  integer,  intent(in) :: i
  real(dp), intent(in) :: dt

  x = tab%xs(:,i) + tab%vs(:,i)*dt
  v = tab%vs(:,i)
  u = tab%us(i)
  h = tab%hs(i)
  if (ns > 0) s(1:ns) = tab%scal(1:ns,i)

 end subroutine sample_state

end subroutine interp_track

!-----------------------------------------------------------------------
!+
!  Spawn order of the replay instances. Instance (k,n) of track k in
!  cycle n is active for t_epoch + n t_rec + a_k <= t < t_epoch + n t_rec + b_k,
!  with a_k = t_entry - t0 and b_k = t_exit - t0. With the entry phase
!  phi_k = a_k mod t_rec and offset c_k = floor(a_k/t_rec), track k is
!  spawned in global cycle M at t_epoch + M t_rec + phi_k as instance
!  n = M - c_k, so that sorting by phi_k gives increasing spawn times.
!+
!-----------------------------------------------------------------------
subroutine schedule_init(sch,tab,t_epoch)
 type(replay_schedule), intent(out) :: sch
 type(track_table),     intent(in)  :: tab
 real(dp),              intent(in)  :: t_epoch
 integer  :: k
 real(dp) :: a

 sch%t_epoch = t_epoch
 sch%t0      = tab%t0
 sch%t_rec   = tab%t_rec
 sch%ntr     = tab%ntracks
 allocate(sch%order(max(sch%ntr,1)),sch%coff(max(sch%ntr,1)),sch%phase(max(sch%ntr,1)))
 do k = 1,sch%ntr
    a = tab%t_entry(k) - tab%t0
    sch%coff(k)  = floor(a/sch%t_rec)
    sch%phase(k) = a - sch%coff(k)*sch%t_rec
    if (sch%phase(k) >= sch%t_rec) then
       sch%phase(k) = sch%phase(k) - sch%t_rec
       sch%coff(k)  = sch%coff(k) + 1
    endif
    if (sch%phase(k) < 0._dp) sch%phase(k) = 0._dp
 enddo
 call sort_index_dp(sch%ntr,sch%phase,sch%order)
 sch%mcycle = 0
 sch%jnext  = 1

end subroutine schedule_init

!-----------------------------------------------------------------------
!+
!  Point the schedule at the first spawn strictly after time t
!+
!-----------------------------------------------------------------------
subroutine schedule_reset(sch,t)
 type(replay_schedule), intent(inout) :: sch
 real(dp),              intent(in)    :: t
 integer  :: m,lo,hi,mid
 real(dp) :: ph

 if (sch%ntr < 1) return
 m  = floor((t - sch%t_epoch)/sch%t_rec)
 ph = t - sch%t_epoch - m*sch%t_rec
 ! first j with phase(order(j)) > ph
 lo = 0
 hi = sch%ntr + 1
 do while (hi - lo > 1)
    mid = (lo + hi)/2
    if (sch%phase(sch%order(mid)) > ph) then
       hi = mid
    else
       lo = mid
    endif
 enddo
 if (hi > sch%ntr) then
    sch%mcycle = m + 1
    sch%jnext  = 1
 else
    sch%mcycle = m
    sch%jnext  = hi
 endif

end subroutine schedule_reset

!-----------------------------------------------------------------------
!+
!  Next instance (k,n) to spawn, if its spawn time is <= t
!+
!-----------------------------------------------------------------------
subroutine schedule_next_spawn(sch,t,k,n,found)
 type(replay_schedule), intent(inout) :: sch
 real(dp),              intent(in)    :: t
 integer,               intent(out)   :: k,n
 logical,               intent(out)   :: found

 found = .false.
 k = 0
 n = 0
 if (sch%ntr < 1) return
 do
    k = sch%order(sch%jnext)
    if (sch%t_epoch + sch%mcycle*sch%t_rec + sch%phase(k) > t) return
    n = sch%mcycle - sch%coff(k)
    sch%jnext = sch%jnext + 1
    if (sch%jnext > sch%ntr) then
       sch%jnext  = 1
       sch%mcycle = sch%mcycle + 1
    endif
    if (n >= 0) then
       found = .true.
       return
    endif
 enddo

end subroutine schedule_next_spawn

!-----------------------------------------------------------------------
!+
!  All instances (k,n) active at time t (start or restart)
!+
!-----------------------------------------------------------------------
subroutine active_instances_at(sch,tab,t,nact,klist,nlist)
 type(replay_schedule), intent(in)  :: sch
 type(track_table),     intent(in)  :: tab
 real(dp),              intent(in)  :: t
 integer,               intent(out) :: nact
 integer, allocatable,  intent(out) :: klist(:),nlist(:)
 integer, allocatable :: ktmp(:),ntmp(:)
 integer  :: k,n,nmin,nmax,ncap
 real(dp) :: a,b

 ncap = max(1024,sch%ntr)
 allocate(klist(ncap),nlist(ncap))
 nact = 0
 do k = 1,sch%ntr
    a = tab%t_entry(k) - tab%t0
    b = tab%t_exit(k)  - tab%t0
    nmax = floor((t - sch%t_epoch - a)/sch%t_rec)
    nmin = max(0,floor((t - sch%t_epoch - b)/sch%t_rec) + 1)
    do n = nmin,nmax
       if (nact >= ncap) then
          ncap = 2*ncap
          allocate(ktmp(ncap),ntmp(ncap))
          ktmp(1:nact) = klist(1:nact)
          ntmp(1:nact) = nlist(1:nact)
          call move_alloc(ktmp,klist)
          call move_alloc(ntmp,nlist)
       endif
       nact = nact + 1
       klist(nact) = k
       nlist(nact) = n
    enddo
 enddo

end subroutine active_instances_at

!-----------------------------------------------------------------------
!+
!  Transformation applied to the tracks of replay cycle n: none, a
!  random rotation (deterministic in n and iseed, identity for n = 0)
!  or the mirror z -> -z for odd n
!+
!-----------------------------------------------------------------------
subroutine cycle_rotation(n,mode,iseed,rot)
 integer,  intent(in)  :: n,mode,iseed
 real(dp), intent(out) :: rot(3,3)
 integer(i8) :: state
 real(dp) :: u1,u2,u3,qw,qx,qy,qz
 integer  :: i

 rot = 0._dp
 do i = 1,3
    rot(i,i) = 1._dp
 enddo
 select case(mode)
 case(irot_random)
    if (n == 0) return
    state = mod(abs(int(iseed,i8))*7919_i8 + abs(int(n,i8))*104729_i8 + 12345_i8,2147483646_i8) + 1_i8
    do i = 1,8
       u1 = park_miller(state)
    enddo
    u1 = park_miller(state)
    u2 = park_miller(state)
    u3 = park_miller(state)
    ! uniform random unit quaternion (Shoemake 1992)
    qw = sqrt(1._dp - u1)*sin(2._dp*pi*u2)
    qx = sqrt(1._dp - u1)*cos(2._dp*pi*u2)
    qy = sqrt(u1)*sin(2._dp*pi*u3)
    qz = sqrt(u1)*cos(2._dp*pi*u3)
    rot(1,:) = [1._dp - 2._dp*(qy*qy + qz*qz), 2._dp*(qx*qy - qw*qz), 2._dp*(qx*qz + qw*qy)]
    rot(2,:) = [2._dp*(qx*qy + qw*qz), 1._dp - 2._dp*(qx*qx + qz*qz), 2._dp*(qy*qz - qw*qx)]
    rot(3,:) = [2._dp*(qx*qz - qw*qy), 2._dp*(qy*qz + qw*qx), 1._dp - 2._dp*(qx*qx + qy*qy)]
 case(irot_mirror)
    if (mod(abs(n),2) == 1) rot(3,3) = -1._dp
 end select

contains

 real(dp) function park_miller(st)
  integer(i8), intent(inout) :: st

  st = mod(16807_i8*st,2147483647_i8)
  park_miller = real(st,dp)/2147483647._dp

 end function park_miller

end subroutine cycle_rotation

end module replayutils
