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
    subroutine MPP_REDUCE_( a, pelist )
      use  mpp_datatype_mod, only: CAFPNTR_TYPE_1D_
!find the max of scalar a the PEs in pelist (all PEs if this argument is omitted)
!result is also automatically broadcast to all PEs
      MPP_TYPE_, intent(inout) :: a
      integer, intent(in), optional :: pelist(:)

      type(CAFPNTR_TYPE_1D_), allocatable, save :: cafptr[:]
      integer :: list,m,n
      MPP_TYPE_,save :: l_work
      MPP_TYPE_, save,target :: work(1)
      logical, save :: first_time=.true.


      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_REDUCE: You must first call mpp_init.' )
      if(first_time)then
        first_time = .false.
        allocate(cafptr[0:*])
      endif

      n = get_peset(pelist); if( peset(n)%count.EQ.1 )return

      if( current_clock.NE.0 )call SYSTEM_CLOCK(start_tick)

      l_work = a; work(1) = a
      cafptr%pfield =>work(1:1)
      call mpp_sync(pelist,do_self=.false.)
      if(n==world_peset_num)then
        do m=0,npes-1
          if(cafptr[m]%pfield(1) CAF_REDUCE_ l_work)l_work = cafptr[m]%pfield(1)
        end do
      else
        list = size(peset(n)%list(:))
        do m=0,list-1
          if(cafptr[peset(n)%list(m+1)]%pfield(1) CAF_REDUCE_ l_work)l_work = cafptr[peset(n)%list(m+1)]%pfield(1)
        end do
      endif
      a = l_work

      call mpp_sync(pelist,do_self=.false.)
      if( current_clock.NE.0 )call increment_current_clock( EVENT_ALLREDUCE, MPP_TYPE_BYTELEN_ )

    end subroutine MPP_REDUCE_
