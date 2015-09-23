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
    subroutine MPP_DO_GLOBAL_FIELD_3Dnew_( local, gfield, d_comm )
!get a global field from a local field
!local field may be on compute OR data domain
      type(DomainCommunicator2D), intent(in) :: d_comm
      MPP_TYPE_, intent(in)  ::  local(:,:,:)   ! not needed for GSM; kept for interface consistency
      MPP_TYPE_, intent(inout) :: gfield(d_comm%domain%x%global%begin:,d_comm%domain%y%global%begin:,:)
      MPP_TYPE_ :: r_field(d_comm%isize_max*d_comm%jsize_max*d_comm%ke)
      pointer( r_ptr,r_field )
      integer :: i, j, k, m, n
      integer :: is, ie, js, je, ke, isize, jsize,istat
      integer :: ii,jj,isizeR,jsizeR
      integer :: list

      integer,save :: done_cnt
      logical, dimension(0:d_comm%Rlist_size-1) :: done

      integer(KIND(gfield)),save :: lsum
      integer(KIND(gfield)),allocatable :: lfield(:)


!     call mpp_error( NOTE, 'Calling gsm MPP_GLOBAL_FIELD' )

      isize = size(gfield,1)
      jsize = size(gfield,2)
      ke = size(gfield,3)

      if(d_comm%isize_out /= isize)then
        call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS: Incorrect I extent.' )
      elseif(d_comm%jsize_out /= jsize)then
        call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS: Incorrect J extent.' )
      elseif(d_comm%ke /= ke)then
        call mpp_error( FATAL, 'MPP_UPDATE_DOMAINS: Incorrect K extent.' )
      endif

      n = d_comm%Rlist_size
!     done = .false.; done_cnt = n
!     gsm_sync_start = Set_Sync_Logical(.true.)
      call mpp_sync(do_self=.false.)
!     do while(done_cnt /= 0); do list=0,n-1
      do list=0,n-1
!        if(done(list))cycle
!        if( Is_Remote_False(d_comm%sync_start_list(list)) )cycle
!        done(list) = .true.; done_cnt = done_cnt - 1
         r_ptr = d_comm%rem_addr(list)
         isizeR = d_comm%isizeR(list)
         jsizeR = d_comm%jsizeR(list)
         is=d_comm%recvis(1,list); ie=d_comm%recvie(1,list)
         js=d_comm%recvjs(1,list); je=d_comm%recvje(1,list)
         do k = 1,ke
           jj = d_comm%sendjsR(1,list)
           do j = js,je
             ii = d_comm%sendisR(1,list)
             do i = is,ie
               gfield(i,j,k) = r_field(ii+(jj-1)*isizeR+(k-1)*isizeR*jsizeR)
               ii = ii + 1
             end do
             jj = jj + 1
           end do
         end do
      end do
!     end do; end do

      if(debug_gsm)then
        allocate(lfield(d_comm%domain%x%global%size * d_comm%domain%y%global%size * d_comm%ke))
        lfield = transfer(gfield,lfield,size(gfield(:,:,:)))
        lsum = mpp_chksum(lfield(1:size(gfield(:,:,:))))
        if(mpp_pe() == mpp_root_pe())write(stdout(),*) 'Global Field checksum=',lsum
        deallocate(lfield)
      endif

      call mpp_sync(do_self=.false.)
!     gsm_sync_start = Set_Sync_Logical(.false.)
    end subroutine MPP_DO_GLOBAL_FIELD_3Dnew_
