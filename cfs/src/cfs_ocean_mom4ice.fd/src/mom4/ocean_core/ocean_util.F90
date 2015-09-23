!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify       !!
!! it and are expected to follow the terms of the GNU General Public !!
!! License as published by the Free Software Foundation.             !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS; if not, write to:                                 !!
!!          Free Software Foundation, Inc.                           !!
!!          59 Temple Place, Suite 330                               !!
!!          Boston, MA  02111-1307  USA                              !!
!! or see:                                                           !!
!!          http://www.gnu.org/licenses/gpl.txt                      !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module ocean_util_mod
!
!<CONTACT EMAIL="Ronald.Pacanowski@noaa.gov"> Ron Pacanowski
!</CONTACT>
!
!<CONTACT EMAIL="Matthew.Harrison@noaa.gov"> Matt Harrison 
!</CONTACT>
!
! <REVIEWER EMAIL="Stephen.Griffies@noaa.gov"> S. M. Griffies 
! </REVIEWER>
!
!<OVERVIEW>
! This module contains many routines of use for mom4. 
!</OVERVIEW>
!
!<DESCRIPTION>
! A utility module for mom4. 
!</DESCRIPTION>
!
use mpp_domains_mod,     only: mpp_update_domains, mpp_global_field
use mpp_mod,             only: stdout, stdlog, FATAL
use mpp_mod,             only: mpp_error, mpp_pe, mpp_root_pe, mpp_min, mpp_max
use time_manager_mod,    only: time_type, get_date
  
use ocean_domains_mod,   only: get_local_indices
use ocean_types_mod,     only: ocean_domain_type

implicit none

private

character(len=256) :: version='CVS $Id'
character(len=256) :: tagname='Tag $Name'

! for output
integer :: unit=6

#include <ocean_memory.h>

type(ocean_domain_type), pointer :: Dom =>NULL()

public iplot
public matrix
public ocean_util_init
public invtri
public min_max
public write_timestamp

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_util_init">
!
! <DESCRIPTION>
! Initialize mom4 utilities.
! </DESCRIPTION>
!
subroutine ocean_util_init (Domain)

  type(ocean_domain_type), intent(in), target :: Domain

  write( stdlog(),'(/a/)') trim(version)

#ifndef STATIC_MEMORY
  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
#endif

  Dom => Domain

end subroutine ocean_util_init
! </SUBROUTINE> NAME="ocean_util_init">


!#######################################################################
! <SUBROUTINE NAME="invtri">
!
! <DESCRIPTION>
! solve the vertical diffusion equation implicitly using the
! method of inverting a tridiagonal matrix as described in
! Numerical Recipes in Fortran, The art of Scientific Computing,
! Second Edition, Press, Teukolsky, Vetterling, Flannery, 1992
! pages 42,43.
! this routine assumes that the variables are defined at grid points
! and the top and bottom b.c. are flux conditions.
!
! inputs:
!
! z         = right hand side terms
!
! nk        = number of vertical levels
!
! topbc     = top boundary condition
!
! botbc     = bottom boundary condition
!
! dcb       = vertical mixing coeff at base of cell
!
! tdt       = 2 * timestep
!
! kmz       = level indicator
!
! mask      = land/sea mask
!
! outputs:
!
! z         = returned solution
!
! </DESCRIPTION>
!
subroutine invtri (z, topbc, botbc, dcb, tdt, kmz, mask, dh, dhw, aidif, nk)

  integer, intent(in)                                 :: nk
  real, intent(inout), dimension(isd:ied,jsd:jed, nk) :: z
  real, intent(in), dimension(isd:ied,jsd:jed,nk)     :: dcb
  real, intent(in), dimension(isd:ied,jsd:jed,nk)     :: dh
  real, intent(in), dimension(isd:ied,jsd:jed,nk)     :: mask
  real, intent(in), dimension(isd:ied,jsd:jed,0:nk)   :: dhw
  real, intent(in), dimension(isd:ied,jsd:jed)        :: topbc
  real, intent(in), dimension(isd:ied,jsd:jed)        :: botbc
  real, intent(in)                                    :: tdt
  integer, intent(in), dimension(isd:ied,jsd:jed)     :: kmz
  real, intent(in)                                    :: aidif

  real, dimension(isd:ied,0:nk) :: a, b, c, e, f
  real, dimension(isd:ied) :: bet

  integer :: i, j, k, km1, kp1
  real :: eps, factu, factl, tdt_aidif

  eps = 1.e-30
  tdt_aidif = tdt*aidif
  do j=jsc,jec
    do k=1,nk
      km1   = max(1,k-1)
      kp1   = min(k+1,nk)
      do i=isc,iec
        factu  = tdt_aidif/(dhw(i,j,k-1)*dh(i,j,k))
        factl  = tdt_aidif/(dhw(i,j,k)*dh(i,j,k))
        a(i,k) = -dcb(i,j,km1)*factu*mask(i,j,k)
        c(i,k) = -dcb(i,j,k)*factl*mask(i,j,kp1)
        f(i,k) = z(i,j,k)*mask(i,j,k) 
        b(i,k) = 1.0 - a(i,k) - c(i,k)
      enddo
    enddo

    do i=isc,iec
      a(i,1)  = 0.0
      c(i,nk) = 0.0
      b(i,1)  = 1.0 - a(i,1) - c(i,1)
      b(i,nk) = 1.0 - a(i,nk) - c(i,nk)

      ! top and bottom b.c.
      f(i,1)  = z(i,j,1) + topbc(i,j)*tdt_aidif*mask(i,j,1)/dh(i,j,1)
      k = max(2,kmz(i,j))
      f(i,k)   = z(i,j,k) - botbc(i,j)*tdt_aidif*mask(i,j,k)/dh(i,j,k)
    enddo

    ! decomposition and forward substitution
    do i=isc,iec
      bet(i) = mask(i,j,1)/(b(i,1) + eps)
      z(i,j,1) = f(i,1)*bet(i)
    enddo
    do k=2,nk
      do i=isc,iec
        e(i,k) = c(i,k-1)*bet(i)
        bet(i) = mask(i,j,k)/(b(i,k) - a(i,k)*e(i,k) + eps)
        z(i,j,k) = (f(i,k) - a(i,k)*z(i,j,k-1))*bet(i)
      enddo
    enddo

    ! back substitution
    do k=nk-1,1,-1
      do i=isc,iec
        z(i,j,k) = z(i,j,k) - e(i,k+1)*z(i,j,k+1)
      enddo
    enddo
  enddo

end subroutine invtri
! </SUBROUTINE> NAME="invtri">


!#######################################################################
! <SUBROUTINE NAME="iplot">
!
! <DESCRIPTION>
!
!  map integer array "iarray" into characters for printing with
!  format (a1) to provide a contour map of the integer field.
!  note: max number of unique characters = 80
!
! inputs:
!
! iarray = integer array to be plotted
!
! is     = starting index along inner dimension of "iarray"
!
! ie     = ending index along inner dimension of "iarray"
!
! js     = starting index along outer dimension of "iarray"
!
! je     = ending index along outer dimension of "iarray"
!
! output: prints contour map of "iarray"
!
! </DESCRIPTION>
!
subroutine iplot (iarray, is, ie, js, je, ni, nj)
  
  integer, intent(in) :: is, ie, js, je, ni, nj
  integer, intent(in), dimension(ni,nj) ::  iarray
  character*80 levels
  character*80 lev1
  integer :: i, j, il, jl, l, incr, ii, jj, inc, last, jinc, maxint, minint
  save levels
  write (stdout(),*) ' '

  ! set character markers
  lev1(1:51) = '.+*ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuv'
  levels = lev1(1:51)//'wxyz0123456789-=!@#$%<>[]{}()'

  ! find range of integers
  maxint = iarray(is,js)
  minint = iarray(is,js)
  do j=js,je
    do i=is,ie
      maxint = max(maxint,iarray(i,j))
      minint = min(minint,iarray(i,j))
    enddo
  enddo

  ! show mapping of integers into characters
  write (stdout(),*) ' '
  write (stdout(),*) ' "iplot" mapping of integers to characters is as follows:'
  inc  = 3
  last = min(minint+80-1,maxint) 
  do i=minint,last,inc
    ii = i-minint+1
    if (i+inc <= last) then
      jinc = inc
    else
      jinc = last-i+1
    endif
    write (stdout(),'(6(1x,i6,a,a,3x))')  (j+minint-1,' is printed as ',levels(j:j),j=ii,ii+jinc-1)
  enddo
  write (stdout(),*) ' '

  if (maxint - minint + 1 > 80) then
    write (stdout(),*) ' Note: there are ',maxint-minint+1,' integers in the field'
    write (stdout(),*) '       "iplot" cannot uniquely assign more than 80 characters for plotting symbols.'
    write (stdout(),*) '       therefore integers are represented by cyclicly reusing the list of plotting symbols'
    write (stdout(),*) ' '
  endif

  ! print character representation of integers
  inc=124
  il = ie-is+1
  jl = je-js+1
  do l=0,il,inc
    incr = min(inc,il-l)
    write (stdout(),8800) (l+i,i=1,incr,4)
    do jj=1,jl
      j = jl+1-jj
      write (stdout(),8900) j, (levels(mod(iarray(l+i+is-1,j+js-1)-minint+1-1,80)+1:&
             mod(iarray(l+i+is-1,j+js-1)-minint+1-1,80)+1),i=1,incr) 
    enddo
  enddo
  8800  format (/, 2x, 31i4)
  8900  format (1x,i3,1x, 124a1)

end subroutine iplot
! </SUBROUTINE> NAME="iplot">


!#######################################################################
! <SUBROUTINE NAME="matrix">
!
! <DESCRIPTION>
! matrix is a general two-dimensional array printing routine,
! input:
!
! array = the array to be printed
!
! istrt = the 1st element of the 1st dimension to be printed
!
! im    = the last element of the 1st dimension to be printed
!
! jstrt = the 1st element of the 2nd dimension to be printed
!
! jm    = the last element of the 2nd dimension to be printed
!         the 2nd dimension is printed in reverse order if both
!         jstrt and jm are negative
!
! scale = a scaling factor by which array is divided before
!         printing.  (if this is zero, no scaling is done.)
!         if scale=0, 10 columns are printed across in e format
!         if scale>0, 20 columns are printed across in f format
!
! output: print "array" as a matrix
!
! </DESCRIPTION>
!
subroutine matrix (array, istrt, im, jstrt, jm, scale)

  integer :: i, l, istrt, im, jstrt, jm, is, ie, js, je, jinc, unit
  real, dimension(istrt:im,abs(jstrt):abs(jm)) ::  array
  real :: scale, scaler

  unit = 6
  if (jstrt*jm < 0) then
    write (unit,999)  jstrt, jm
    call mpp_error(FATAL,'==>Error in ocean_util_mod (matrix): jstrt*jm < 0 found in matrix')
  endif

  ! allow for inversion of 2nd dimension
  if (jm < 0) then
    js   = -jm
    je   = -jstrt
    jinc = -1
  else
    js   = jstrt
    je   = jm
    jinc = 1
  endif

  if (scale == 0.0) then

    do is=istrt,im,10
      ie = min(is + 9,im)
      write (unit,9001) (i, i=is,ie)
      do l=js,je,jinc
        write (unit,9002) l, (array(i,l),i=is,ie)
      enddo
      write (unit,'(/)')
    enddo
  else
    scaler = 1.0/scale
    do is=istrt,im,20
      ie = min(is + 19,im)
      write (unit,9003) (i, i=is,ie)
      do l=js,je,jinc
        write (unit,9004) l, (array(i,l)*scaler,i=is,ie)
      enddo
      write (unit,'(/)')
    enddo
  endif
  999   format (1x,'jstrt=',i5,' jm=',i5,' in matrix')
  9001  format(10i13)
  9002  format(i3,10(1pe13.5))
  9003  format(3x,20i6)
  9004  format(1x,i3,1x,20f6.2)

end subroutine matrix
! </SUBROUTINE> NAME="matrix">


!#######################################################################
! <SUBROUTINE NAME="min_max">
!
! <DESCRIPTION>
! Determine min and max
! </DESCRIPTION>
!
subroutine min_max (a, nk0, mask, text, by_depth)

  integer, intent(in) :: nk0
  real, intent(in), dimension(isd:ied,jsd:jed,nk0) :: a, mask
  character*(*), intent(in) :: text
  character*(*), intent(in), optional :: by_depth
  real :: mina, maxa, mina0, maxa0, fudge
  integer :: i, j, k, unit, imin, jmin, kmin, imax, jmax, kmax
  
  unit  = 6
  fudge = 1.0 + 1.e-12*mpp_pe()  ! to distinguish processors when "a" is independent of processor
  if (present(by_depth)) then
    do k=1,nk0
      mina = 1.e20; maxa = -1.e20
      do j=jsc,jec
        do i=isc,iec
          if (mask(i,j,k) == 1.0) then
            if (mina > a(i,j,k)) then
              mina = a(i,j,k); imin=i; jmin=j
            endif
            if (maxa < a(i,j,k)) then
              maxa = a(i,j,k); imax=i; jmax=j
            endif
          endif       
        enddo
      enddo
      mina  = mina*fudge
      mina0 = mina
      call mpp_min(mina)
      if (mina0 == mina) then
        write (unit,*) trim(text),', k=',k,': min =',mina/fudge,' located at (i,j) = (',imin,',',jmin,')'
      endif
      maxa  = maxa*fudge
      maxa0 = maxa
      call mpp_max(maxa)
      if (maxa0 == maxa) then
        write (unit,*) trim(text),', k=',k,': max =',maxa/fudge,' located at (i,j) = (',imax,',',jmax,')'
      endif
    enddo
  else
    mina = 1.e20; maxa = -1.e20
    do k=1,nk0
      do j=jsc,jec
        do i=isc,iec
          if (mask(i,j,k) == 1.0) then
            if (mina > a(i,j,k)) then
              mina = a(i,j,k); imin=i; jmin=j; kmin=k
            endif
            if (maxa < a(i,j,k)) then
              maxa = a(i,j,k); imax=i; jmax=j; kmax=k
            endif
          endif       
        enddo
      enddo
    enddo
    mina  = mina*fudge
    mina0 = mina
    call mpp_min(mina)
    if (mina0 == mina) then
      if (nk0 == 1) then
        write (unit,*) trim(text),': min =',mina/fudge,' located at (i,j) = (',imin,',',jmin,')'
      else
        write (unit,*) trim(text),': min =',mina/fudge,' located at (i,j,k) = (',imin,',',jmin,',',kmin,')'
      endif
    endif
    maxa  = maxa*fudge
    maxa0 = maxa
    call mpp_max(maxa)
    if (maxa0 == maxa) then
      if (nk0 == 1) then
        write (unit,*) trim(text),': max =',maxa/fudge,' located at (i,j) = (',imax,',',jmax,')'
      else
        write (unit,*) trim(text),': max =',maxa/fudge,' located at (i,j,k) = (',imax,',',jmax,',',kmax,')'
      endif
    endif
  endif
  write (stdout(),*) ' '

end subroutine min_max
! </SUBROUTINE> NAME="min_max">


!#######################################################################
! <SUBROUTINE NAME="write_timestamp">
!
! <DESCRIPTION>
! Write the time stamp.
! </DESCRIPTION>
!
subroutine write_timestamp(Time)

  type(time_type), intent(in) :: Time

  integer :: yr, mon, day, hr, min, sec

  call get_date(Time, yr, mon, day, hr, min, sec)

  write(stdout(),'(a,i4,a,i2,a,i2,1x,i2,a,i2,a,i2)' ) &
       'yyyy/mm/dd hh:mm:ss = ', yr, '/',mon,'/',day, hr, ':',min,':', sec

end subroutine write_timestamp
! </SUBROUTINE> NAME="min_max">




end module ocean_util_mod
