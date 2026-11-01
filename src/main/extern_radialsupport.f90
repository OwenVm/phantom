!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module extern_radialsupport
!
! This module contains routines relating to the computation of a
! tabulated radial force on the gas, which supports the extended
! atmosphere in the pulsating wind setup. The table is passed in by the
! inject module, which is compiled after this one. The force is called
! directly from get_force, so it does not need an iexternalforce.
!
! :References: None
!
! :Owner: Owen Vermeulen
!
! :Runtime parameters: None
!
! :Dependencies: physcon
!
 use physcon, only:pi
 implicit none

 public :: set_radial_support,radial_support_is_set,get_radial_support
 private

 integer, parameter :: nsupport = 4096
 logical :: support_is_set = .false.
 real    :: a_support(nsupport),r_support_in,r_support_out,dr_support
 real    :: r_support_start,r_support_width

contains

!-----------------------------------------------------------------------
!+
!  Set a radial force on the gas from a table of acceleration (positive
!  outward) against radius. The table is resampled on a uniform grid so
!  that no search is needed in the force loop. The force is on up to
!  rstart, then goes to zero as cos^2 over rwidth (never if rwidth <= 0).
!+
!-----------------------------------------------------------------------
subroutine set_radial_support(rtab,atab,rstart,rwidth)
 real, intent(in) :: rtab(:),atab(:),rstart,rwidth
 integer :: i,j,n
 real    :: ri

 n = size(rtab)
 r_support_in  = rtab(1)
 r_support_out = rtab(n)
 dr_support    = (rtab(n) - rtab(1))/real(nsupport - 1)
 j = 1
 do i=1,nsupport
    ri = rtab(1) + (i-1)*dr_support
    do while (j < n-1 .and. rtab(j+1) < ri)
       j = j + 1
    enddo
    a_support(i) = atab(j) + (ri - rtab(j))*(atab(j+1) - atab(j))/(rtab(j+1) - rtab(j))
 enddo
 r_support_start = rstart
 r_support_width = rwidth
 support_is_set  = .true.

end subroutine set_radial_support

!-----------------------------------------------------------------------
!+
!  Query function for whether the radial force has been set
!+
!-----------------------------------------------------------------------
logical function radial_support_is_set()

 radial_support_is_set = support_is_set

end function radial_support_is_set

!-----------------------------------------------------------------------
!+
!  Add the radial force at (xi,yi,zi), directed away from the point x0.
!  Below the table the innermost value is used. Above it the force falls
!  as 1/r^2, so that a particle leaving the atmosphere sees the same
!  fraction of gravity as at the edge, until the cos^2 decline removes it.
!+
!-----------------------------------------------------------------------
subroutine get_radial_support(xi,yi,zi,x0,fextxi,fextyi,fextzi)
 real, intent(in)    :: xi,yi,zi,x0(3)
 real, intent(inout) :: fextxi,fextyi,fextzi
 real    :: dx,dy,dz,r,rgrid,rtaper,acc
 integer :: i

 dx = xi - x0(1)
 dy = yi - x0(2)
 dz = zi - x0(3)
 r  = sqrt(dx*dx + dy*dy + dz*dz)
 if (r <= 0.) return

 rtaper = r - r_support_start
 if (r_support_width > 0. .and. rtaper >= r_support_width) return

 rgrid = (r - r_support_in)/dr_support
 if (rgrid <= 0.) then
    acc = a_support(1)
 elseif (rgrid >= nsupport - 1) then
    acc = a_support(nsupport)*(r_support_out/r)**2
 else
    i = int(rgrid) + 1
    acc = a_support(i) + (rgrid - (i - 1))*(a_support(i+1) - a_support(i))
 endif

 if (r_support_width > 0. .and. rtaper > 0.) acc = acc*cos(0.5*pi*rtaper/r_support_width)**2

 acc = acc/r
 fextxi = fextxi + acc*dx
 fextyi = fextyi + acc*dy
 fextzi = fextzi + acc*dz

end subroutine get_radial_support

end module extern_radialsupport
