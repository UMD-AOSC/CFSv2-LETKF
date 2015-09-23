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
    subroutine MPP_DO_REDISTRIBUTE_3Dnew_( f_in, f_out, d_comm, d_type )
      integer(LONG_KIND), intent(in)                 :: f_in(:), f_out(:)
      type(DomainCommunicator2D), intent(in), target :: d_comm
      MPP_TYPE_, intent(in)                          :: d_type

      MPP_TYPE_ :: Rfield_in(d_comm%isize_max*d_comm%jsize_max*d_comm%ke)
      pointer( ptr_field_in,Rfield_in )  ! pointer to remote field
      MPP_TYPE_ :: field_out(d_comm%domain_out%x%data%begin:d_comm%domain_out%x%data%end, &
                             d_comm%domain_out%y%data%begin:d_comm%domain_out%y%data%end,d_comm%ke)
      pointer( ptr_field_out, field_out)
      real(KIND(d_type)) :: r_field(d_comm%domain_out%x%data%size,d_comm%domain_out%y%data%size,d_comm%ke)
      pointer( r_ptr,r_field )
      real(KIND(d_type)),save :: d_field(1,1,1)=0.d0  ! dummy for use with mpp_chksum
      integer :: i, j, k, l, m, n, l_size
      integer :: is, ie, js, je, istat
      integer :: isize_in, jsize_in, ke_in, isize_out, jsize_out, ke_out
      integer :: ii,jj,isizeR,jsizeR
      integer :: list
      integer,save :: done_cnt
      logical, dimension(0:d_comm%Rlist_size-1) :: done


      l_size = size(f_out(:))
      isize_out = d_comm%isize_out
      jsize_out = d_comm%jsize_out
      ke_out = d_comm%ke

!recv    
      n = d_comm%Rlist_size
!     done = .false. .OR. .NOT.d_comm%R_do_buf(:)
!     done_cnt = n - count(done)
!     gsm_sync_start = Set_Sync_Logical(.true.)
      call mpp_sync(do_self=.false.)
!     do while(done_cnt /= 0); do list=0,n-1
      do list=0,n-1
         if( .NOT. d_comm%R_do_buf(list) )cycle
!        if(done(list))cycle
!        if( Is_Remote_False(d_comm%sync_start_list(list)) )cycle
!        done(list) = .true.; done_cnt = done_cnt - 1
         isizeR = d_comm%isizeR(list)
         jsizeR = d_comm%jsizeR(list)
         is=d_comm%recvis(1,list); ie=d_comm%recvie(1,list)
         js=d_comm%recvjs(1,list); je=d_comm%recvje(1,list)
         do l=1,l_size  ! loop over number of in/out fields
           ptr_field_in = d_comm%rem_addrl(l,list)
           ptr_field_out = f_out(l)
           do k = 1,ke_out
             jj = d_comm%sendjsR(1,list)
             do j = js,je
               ii = d_comm%sendisR(1,list)
               do i = is,ie
                  field_out(i,j,k) = Rfield_in(ii+(jj-1)*isizeR+(k-1)*isizeR*jsizeR)
                  ii = ii + 1
               end do
               jj = jj + 1
             end do
           end do
         end do
      end do
!     end do; end do
      call mpp_sync(do_self=.false.)
!     gsm_sync_start = Set_Sync_Logical(.false.)

      if(debug_gsm)then
        do l=1,l_size  ! loop over number of in/out fields
          if(ANY(d_comm%R_do_buf(:)))then
            r_ptr = f_out(l)
            write(stdout(),*) 'Domain checksum=',mpp_chksum(r_field)
          else
            write(stdout(),*) 'Domain checksum=',mpp_chksum(d_field)
          endif
        end do 
      endif
    end subroutine MPP_DO_REDISTRIBUTE_3Dnew_
