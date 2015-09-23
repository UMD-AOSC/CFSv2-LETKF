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
    subroutine MPP_SUM_( a, length, pelist )
      use  mpp_datatype_mod, only: CAFPNTR_TYPE_1D_
!sums array a over the PEs in pelist (all PEs if this argument is omitted)
!result is also automatically broadcast: all PEs have the sum in a at the end
  !we are using f77-style call: array passed by address and not descriptor; further, 
  !the f90 conformance check is avoided.
      integer, intent(in) :: length
      integer, intent(in), optional :: pelist(:)
      MPP_TYPE_, intent(inout),target :: a(*)

      type(CAFPNTR_TYPE_1D_), allocatable, save :: cafptr[:]
      integer :: list,k,m,n
      MPP_TYPE_ :: l_work(length)
      logical, save :: first_time=.true.


      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_SUM: You must first call mpp_init.' )
      if(first_time)then
        first_time = .false.
        allocate(cafptr[0:*])
      endif

      n = get_peset(pelist); if( peset(n)%count.EQ.1 )return

      if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)

      l_work = 0.d0
      if(n==world_peset_num)then
         if(length == 1)then
           cafptr%pfield =>a(1:1)
           call mpp_sync(pelist,do_self=.false.)
           do k=0,npes-1
             l_work(1) = l_work(1) + cafptr[k]%pfield(1)
           end do
         else
           cafptr%pfield =>a(1:length)
           call mpp_sync(pelist,do_self=.false.)
           do k=0,npes-1
             do m=1,length
               l_work(m) = l_work(m) +  cafptr[k]%pfield(m)
             end do
           end do
         endif

      else

        list = size(peset(n)%list(:))
        if(length == 1)then
           cafptr%pfield =>a(1:1)
           call mpp_sync(pelist,do_self=.false.)
           do k=0,list-1
             l_work(1) = l_work(1) + cafptr[peset(n)%list(k+1)]%pfield(1)
           end do
         else
           cafptr%pfield =>a(1:length)
           call mpp_sync(pelist,do_self=.false.)
           do k=0,list-1
             do m=1,length
               l_work(m) = l_work(m) + cafptr[peset(n)%list(k+1)]%pfield(m)
             end do
           end do
         endif
      endif
      a(1:length) = l_work(1:length)

      call mpp_sync(pelist,do_self=.false.)
      if( current_clock.NE.0 )call increment_current_clock( EVENT_ALLREDUCE, length*MPP_TYPE_BYTELEN_ )

    end subroutine MPP_SUM_

!#######################################################################
#include <mpp_sum.inc>
