!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! Initial conditions of an outer run for periodic track replay, from a
! dump of the inner run (normally the first dump of the replay period):
! the boundary (piston) particles and all gas inside r_inj are removed,
! the sinks can be merged into one at their centre of mass, and rkill and
! the replay options are set. The outer run builds the track table from
! this dump and the following ones at its start.
!
! :References: None
!
! :Owner: Owen Vermeulen
!
! :Runtime parameters: None
!
! :Dependencies: inject, injection, io, part, physcon, prompting, units
!
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,      only:xyzmh_ptmass,vxyz_ptmass,nptmass,iphase,iamtype,igas,iboundary, &
                     kill_particle,shuffle_part,isdead_or_accreted,ilum,maxphase,maxp
 use inject,    only:set_replay_options
 use injection, only:rkill
 use prompting, only:prompt
 use units,     only:udist
 use physcon,   only:au
 use io,        only:fatal
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 character(len=200) :: dumps
 real    :: xc(3),vc(3),mtot,r_out_au,ltot,dx(3),r_inj_au,period_days,buffer_au,r_inj
 integer :: i,k,nbound,ngas
 logical :: merge

 if (maxphase /= maxp) call fatal('moddump_replay','particle types are needed')

 ! the dump given to phantommoddump is the natural start of the replay period
 call get_command_argument(1,dumps)
 call prompt('first dump of the replay period (path as seen from the outer run)',dumps)
 r_inj_au    = 50.
 period_days = 300.
 buffer_au   = 0.
 call prompt('interface radius r_inj (au)',r_inj_au,0.)
 call prompt('replay period = pulsation or orbital period (days)',period_days,0.)
 call prompt('width of the replay buffer (au, 0=auto from h)',buffer_au,0.)
 r_inj = r_inj_au*au/udist

 ! centre of mass of the sinks
 xc   = 0.
 vc   = 0.
 mtot = 0.
 do k = 1,nptmass
    xc   = xc + xyzmh_ptmass(4,k)*xyzmh_ptmass(1:3,k)
    vc   = vc + xyzmh_ptmass(4,k)*vxyz_ptmass(1:3,k)
    mtot = mtot + xyzmh_ptmass(4,k)
 enddo
 if (mtot > 0.) then
    xc = xc/mtot
    vc = vc/mtot
 endif

 ! remove the piston and the gas inside r_inj
 nbound = 0
 ngas   = 0
 do i = 1,npart
    if (isdead_or_accreted(xyzh(4,i))) cycle
    if (iamtype(iphase(i)) == iboundary) then
       call kill_particle(i,npartoftype)
       nbound = nbound + 1
    elseif (iamtype(iphase(i)) == igas) then
       dx = xyzh(1:3,i) - xc
       if (dot_product(dx,dx) <= r_inj**2) then
          call kill_particle(i,npartoftype)
          ngas = ngas + 1
       endif
    endif
 enddo
 call shuffle_part(npart)
 print "(a,i10,a,i10,a)",' removed ',nbound,' boundary particles and ',ngas,' gas particles inside r_inj'
 print "(a,i10,a)",' kept ',npart,' particles'

 ! merge the sinks
 if (nptmass > 1) then
    merge = .true.
    call prompt('merge the sinks into one at their centre of mass?',merge)
    if (merge) then
       ltot = sum(xyzmh_ptmass(ilum,1:nptmass))
       xyzmh_ptmass(1:3,1)  = xc
       vxyz_ptmass(1:3,1)   = vc
       xyzmh_ptmass(4,1)    = mtot
       xyzmh_ptmass(ilum,1) = ltot
       xyzmh_ptmass(:,2:nptmass) = 0.
       vxyz_ptmass(:,2:nptmass)  = 0.
       nptmass = 1
       print "(a,es12.4,a)",' merged the sinks, total mass ',mtot,' (code units)'
    endif
 endif

 r_out_au = 20.*r_inj_au
 call prompt('outer radius rkill of the outer run (au)',r_out_au,r_inj_au)
 rkill = r_out_au*au/udist

 massoftype(iboundary) = massoftype(igas)
 call set_replay_options(trim(dumps),period_days,r_inj_au,buffer_au)
 print "(a)",' set rkill and the replay options in the input file'

end subroutine modify_dump

end module moddump
