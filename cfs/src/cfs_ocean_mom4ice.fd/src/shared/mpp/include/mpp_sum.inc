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
!#######################################################################

    subroutine MPP_SUM_SCALAR_( a, pelist )
!sums array a when only first element is passed: this routine just converts to a call to MPP_SUM_
      MPP_TYPE_, intent(inout) :: a
      integer, intent(in), optional :: pelist(:)
      MPP_TYPE_ :: b(1)

      b(1) = a
      if( debug )call mpp_error( NOTE, 'MPP_SUM_SCALAR_: calling MPP_SUM_ ...' )
      call MPP_SUM_( b, 1, pelist )
      a = b(1)
      return
    end subroutine MPP_SUM_SCALAR_

!#######################################################################
    subroutine MPP_SUM_2D_( a, length, pelist )
      MPP_TYPE_, intent(inout) :: a(:,:)
      integer, intent(in) :: length
      integer, intent(in), optional :: pelist(:)
      MPP_TYPE_ :: a1D(length)

#ifdef use_CRI_pointers
      pointer( ptr, a1D )
      ptr = LOC(a)
      call mpp_sum( a1D, length, pelist )
#else
      a1D = TRANSFER( a, a1D, length ) !faster than RESHAPE? length is probably redundant
!      a1D = RESHAPE( a, SHAPE(a1D) )
      call mpp_sum( a1D, length, pelist )
      a = RESHAPE( a1D, SHAPE(a) )
#endif
      return
    end subroutine MPP_SUM_2D_

!#######################################################################
    subroutine MPP_SUM_3D_( a, length, pelist )
      MPP_TYPE_, intent(inout) :: a(:,:,:)
      integer, intent(in) :: length
      integer, intent(in), optional :: pelist(:)
      MPP_TYPE_ :: a1D(length)

#ifdef use_CRI_pointers
      pointer( ptr, a1D )
      ptr = LOC(a)
      call mpp_sum( a1D, length, pelist )
#else
      a1D = TRANSFER( a, a1D, length ) !faster than RESHAPE? length is probably redundant
!      a1D = RESHAPE( a, SHAPE(a1D) )
      call mpp_sum( a1D, length, pelist )
      a = RESHAPE( a1D, SHAPE(a) )
#endif
      return
    end subroutine MPP_SUM_3D_

!#######################################################################
    subroutine MPP_SUM_4D_( a, length, pelist )
      MPP_TYPE_, intent(inout) :: a(:,:,:,:)
      integer, intent(in) :: length
      integer, intent(in), optional :: pelist(:)
      MPP_TYPE_ :: a1D(length)

#ifdef use_CRI_pointers
      pointer( ptr, a1D )
      ptr = LOC(a)
      call mpp_sum( a1D, length, pelist )
#else
      a1D = TRANSFER( a, a1D, length ) !faster than RESHAPE? length is probably redundant
!      a1D = RESHAPE( a, SHAPE(a1D) )
      call mpp_sum( a1D, length, pelist )
      a = RESHAPE( a1D, SHAPE(a) )
#endif
      return
    end subroutine MPP_SUM_4D_

!#######################################################################
    subroutine MPP_SUM_5D_( a, length, pelist )
      MPP_TYPE_, intent(inout) :: a(:,:,:,:,:)
      integer, intent(in) :: length
      integer, intent(in), optional :: pelist(:)
      MPP_TYPE_ :: a1D(length)

#ifdef use_CRI_pointers
      pointer( ptr, a1D )
      ptr = LOC(a)
      call mpp_sum( a1D, length, pelist )
#else
      a1D = TRANSFER( a, a1D, length ) !faster than RESHAPE? length is probably redundant
!      a1D = RESHAPE( a, SHAPE(a1D) )
      call mpp_sum( a1D, length, pelist )
      a = RESHAPE( a1D, SHAPE(a) )
#endif
      return
    end subroutine MPP_SUM_5D_
