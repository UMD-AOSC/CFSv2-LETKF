      subroutine mpi_tracers_a_to_b(rgt_a,lats_nodes_a,global_lats_a,
     &                    for_gr_r_2,lats_nodes_r,global_lats_r,
     &                    k_pointer,i_print)

!     initialize for_gr_r_2 tracers to rgt_a
!     from gloopa by copy or send-receive

      use machine   , only : kind_evod
      use resol_def , only : latg,latr,levs,lonf,lonr,lots,ntrac
      use layout1   , only : me,nodes,
     &                       ipt_lats_node_a,lats_dim_a,lats_node_a,
     &                       ipt_lats_node_r,lats_dim_r,lats_node_r
      use mpi_def   , only : mc_comp,mpi_r_def,stat

      implicit none

      real(kind=kind_evod) rgt_a(lonf,levs,lats_dim_a,ntrac)  !  input

!mjr  real(kind=kind_evod) for_gr_r_2(lonrx,lots,lats_dim_r)  !  output
      real(kind=kind_evod) for_gr_r_2(lonr ,lots,lats_dim_r)  !  output

      integer  lats_nodes_a(nodes), global_lats_a(latg)
      integer  lats_nodes_r(nodes), global_lats_r(latr)

      integer  k_pointer
      integer  i_print  !  =1 for diagnostic printing

      integer  i
      integer  i1
      integer  i2
      integer  ierr
      integer  iprint
      integer  j
      integer  jj
      integer  k, ktem
      integer  lan_a
      integer  lan_b
      integer  lat
      integer  n, ntr
      integer  num_copy
      integer  num_recv
      integer  num_send

      integer  global_lats_a_sn(latg)
ccmr  integer  global_lats_r_sn(latr)

      integer req_recv(latg)
      integer req_send(latg)

!mjr  integer lan_b_recv(latg)

!mjr  real(kind=kind_evod) grid_recv(lonf,levs,3,lats_node_r)
      real(kind=kind_evod) grid_send(lonf,levs,ntrac,lats_node_a)

      integer, parameter :: itaga=2**5, itagb=2**4
      integer i_count
      save    i_count
      data    i_count / 0 /
              i_count = i_count + 1

!-----------------------------------------------------------------------

      i1=1
      i2=0
      do n=1,nodes
         j=0
         i2=i2+lats_nodes_a(n)
         do i=i1,i2
            j=j+1
            global_lats_a_sn(i)=global_lats_a(i2+1-j)
         enddo
         i1=i1+lats_nodes_a(n)
      enddo

!-----------------------------------------------------------------------

ccmr  i1=1
ccmr  i2=0
ccmr  do n=1,nodes
ccmr     j=0
ccmr     i2=i2+lats_nodes_r(n)
ccmr     do i=i1,i2
ccmr        j=j+1
ccmr        global_lats_r_sn(i)=global_lats_r(i2+1-j)
ccmr     enddo
ccmr     i1=i1+lats_nodes_r(n)
ccmr  enddo

!-----------------------------------------------------------------------

           iprint =  0
      if ( iprint.eq.1 .and. me.eq.0 .and. i_count.eq.1 ) then

         j=0
         do n=1,nodes
            if (lats_nodes_a(n).gt.0) then
               write(me+7500,'(" ")')
               do jj=1,lats_nodes_a(n)
                  j=j+1
                  write(me+7500,'(2i4,
     &                      "     global_lats_a("   ,i4,")=",i4,
     &                      "     global_lats_a_sn(",i4,")=",i4)')
     &                  n, jj, j, global_lats_a(j),
     &                         j, global_lats_a_sn(j)
               enddo
            endif
         enddo

ccmr           write(me+7500,'(" ")')
ccmr           write(me+7500,'("------------------------------")')
ccmr           write(me+7500,'(" ")')
ccmr     j=0
ccmr     do n=1,nodes
ccmr        if (lats_nodes_r(n).gt.0) then
ccmr           write(me+7500,'(" ")')
ccmr           do jj=1,lats_nodes_r(n)
ccmr              j=j+1
ccmr              write(me+7500,'(2i4,
ccmr &                      "     global_lats_r("   ,i4,")=",i4,
ccmr &                      "     global_lats_r_sn(",i4,")=",i4)')
ccmr &                  n, jj, j, global_lats_r(j),
ccmr &                         j, global_lats_r_sn(j)
ccmr           enddo
ccmr        endif
ccmr     enddo

         close(me+7500)

      endif

!-----------------------------------------------------------------------

cc    processor copies or receives values for its own latitudes
      num_copy=0  !  set number copy    count to zero
      num_recv=0  !  set number receive count to zero
      do lan_b=1,lats_node_r
cc       lat is latitude number from 1 to latr
         lat = global_lats_r(ipt_lats_node_r-1+lan_b)
         i=0  !  set pointer for all latitudes to zero
         do n=1,nodes  !  increment through all computational tasks
            do lan_a=1,lats_nodes_a(n)
               i=i+1  !  increment pointer for all latitudes
               if (lat .eq. global_lats_a_sn(i)) then
                  if (n-1.eq.me) then
cc                   values are in same processor
                     num_copy=num_copy+1  !  increment number copy    count
                     do ntr=1,ntrac
                       do k=1,levs
!mjr                    do i=1,min(lonf,lons_lat)
                         ktem = k_pointer-1+k+(ntr-1)*levs
                         do i=1,lonf
                           for_gr_r_2(i,ktem,lan_b) =
     &                          rgt_a(i,k,lan_a,ntr)
                         enddo
                       enddo
                     enddo
                  else
cc                   values are in different processors
                     num_recv=num_recv+1  !  increment number receive count
!mjr                 lan_b_recv(num_recv)=lan_b
                     call mpi_irecv(  !  nonblocking receive operation
!mjr &                    grid_recv(1,1,1,num_recv),      !  receive buffer
     &                    for_gr_r_2(1,k_pointer,lan_b),  !  receive buffer
     &                    lonf*levs*ntrac,  !  number of elements
     &                    mpi_r_def,  !  data type
     &                    n-1,  !  rank of source process
!Moor&                    1*10**8+(n-1)*10**4+lat,  !  message tag
     &                    1*itaga+(n-1)*itagb+lat,  !  message tag
     &                    mc_comp,  !  communicator handle
     &                    req_recv(num_recv),  !  communication request
     &                    ierr)  !  out: return code
                  endif
                  go to 12345
               endif
            enddo ! do lan_a
         enddo ! do n
12345    continue
      enddo ! do lan_b

!-----------------------------------------------------------------------

cc    processor sends values of its own latitudes
cc    to other processors' latitudes
      num_send=0  !  set number send count to zero
      i=0  !  set pointer for all latitudes to zero
      do n=1,nodes  !  increment through all computational tasks
        do lan_b=1,lats_nodes_r(n)
          i=i+1  !  increment pointer for all latitudes
          do lan_a=1,lats_node_a
cc          lat is latitude number from 1 to latg
            lat = global_lats_a_sn(ipt_lats_node_a-1+lan_a)
            if (lat .eq. global_lats_r(i)) then
              if (n-1.ne.me) then
                 num_send=num_send+1              !  increment number send count
                 do ntr=1,ntrac
                   grid_send(:,:,ntr,num_send) = rgt_a(:,:,lan_a,ntr)
                 enddo
                 call mpi_isend(                  !  nonblocking send operation
     &                grid_send(1,1,1,num_send),  !  send buffer
     &                lonf*levs*ntrac,            !  number of elements
     &                mpi_r_def,                  !  data type
     &                n-1,                        !  rank of destination process
!Moor&                1*10**8+me*10**4+lat,       !  message tag
     &                1*itaga+me*itagb+lat,       !  message tag
     &                mc_comp,                    !  communicator handle
     &                req_send(num_send),         !  communication request
     &                ierr)                       !  out: return code
                  endif
                  go to 23456
               endif
            enddo ! do lan_a
23456       continue
         enddo ! do lan_b
      enddo ! do n

!-----------------------------------------------------------------------

ccmr  call mpi_barrier (mc_comp,ierr)

      do j=1,num_recv  !  increment through number receive count
cc       waits for nonblocking receive operation to complete
         call mpi_wait(
     &        req_recv(j),  !  communication request to wait for
     &        stat,           !  out: status object
     &        ierr)           !  out: return code
      enddo

      do j=1,num_send  !  increment through number send count
cc       waits for nonblocking send operation to complete
         call mpi_wait(
     &        req_send(j),  !  communication request to wait for
     &        stat,           !  out: status object
     &        ierr)           !  out: return code
      enddo

!-----------------------------------------------------------------------

!mjr  do lan_b=1,num_recv
!mjr  do n=1,3
!mjr  do k=1,levs
!mjr  do i=1,min(lonf,lons_lat)
!mjr  do i=1,lonf
!mjr     for_gr_r_2(i,k_pointer-1+k+(n-1)*levs,lan_b_recv(lan_b))=
!mjr &    grid_recv(i,k,n,lan_b)
!mjr  enddo
!mjr  enddo
!mjr  enddo
!mjr  enddo

!-----------------------------------------------------------------------

      if ( i_print .eq. 1 ) then

         print*,' '
         print*,' mpi_tracers_a_to_b ',
     &          ' k_pointer=', k_pointer,
     &          ' i_count=',   i_count,
     &          ' num_copy=',  num_copy,
     &          ' num_recv=',  num_recv,
     &          ' num_send=',  num_send
         print*,' '

         do lan_a=1,lats_node_a
            lat = global_lats_a(ipt_lats_node_a-1+lan_a)
            write(*,'("lan_a=",i5," lat=",i5,"      rgt_a=",4e22.12)')
     &           lan_a, lat,
     &           rgt_a(1,     1,lan_a,:)
         enddo
         print*,' '

         do lan_b=1,lats_node_r
            lat = global_lats_r(ipt_lats_node_r-1+lan_b)
            write(*,'("lan_b=",i5," lat=",i5," for_gr_r_2=",4e22.12)')
     &           lan_b, lat,
     &           for_gr_r_2(1,k_pointer-1+     1,lan_b),
     &           for_gr_r_2(1,k_pointer-1+     2,lan_b),
     &           for_gr_r_2(1,k_pointer-1+levs-1,lan_b),
     &           for_gr_r_2(1,k_pointer-1+levs  ,lan_b)
         enddo
         print*,' '

      endif

!-----------------------------------------------------------------------

      return
      end
      subroutine mpi_tracers_b_to_a(
     &          bak_gr_r_2,lats_nodes_r,global_lats_r,
     &          rgt_h,lats_nodes_a,global_lats_a,
     &          k_pointer,i_print)

!     initialize rg1_h rg2_h rg3_h for gloopa
!     from bak_gr_r_2 tracers by copy or send-receive

      use machine             , only : kind_evod
      use resol_def           , only : latg,latr,levs,lonf,lonr,lota
     &,                                ntrac
      use layout1             , only : me,nodes,
     &                                 ipt_lats_node_a,lats_node_a,
     &                                 ipt_lats_node_r,lats_node_r,
     &                                 lats_dim_r
      use layout_lag          , only : lats_dim_h
      use layout_grid_tracers , only : xhalo,yhalo
      use mpi_def             , only : mc_comp,mpi_r_def,stat

      implicit none

!mjr  real(kind=kind_evod) bak_gr_r_2(lonrx,lota,lats_dim_r)  !  input
      real(kind=kind_evod) bak_gr_r_2(lonr ,lota,lats_dim_r)  !  input

      real(kind=kind_evod) rgt_h(lonf+1+2*xhalo,levs,lats_dim_h,ntrac) ! output

      integer  lats_nodes_r(nodes), global_lats_r(latr)
      integer  lats_nodes_a(nodes), global_lats_a(latg)

      integer  k_pointer
      integer  i_print  !  =1 for diagnostic printing

      integer  i
      integer  i1
      integer  i2
      integer  ierr
      integer  iprint
      integer  j
      integer  jj
      integer  k, kk, ktem
      integer  lan_a
      integer  lan_b
      integer  lat, ll
      integer  n, ntr
      integer  num_copy
      integer  num_recv
      integer  num_send

      integer  global_lats_a_sn(latg)
ccmr  integer  global_lats_r_sn(latr)

      integer req_recv(latg)
      integer req_send(latg)

      integer lan_a_recv(latg)

!mjr  real(kind=kind_evod) grid_recv(lonrx,levs,3,lats_node_a)
      real(kind=kind_evod) grid_recv(lonr ,levs,3,lats_node_a)

      integer, parameter :: itaga=2**5, itagb=2**4
      integer i_count
      save    i_count
      data    i_count / 0 /
              i_count = i_count + 1

!-----------------------------------------------------------------------

      i1=1
      i2=0
      do n=1,nodes
         j=0
         i2=i2+lats_nodes_a(n)
         do i=i1,i2
            j=j+1
            global_lats_a_sn(i)=global_lats_a(i2+1-j)
         enddo
         i1=i1+lats_nodes_a(n)
      enddo

!-----------------------------------------------------------------------

ccmr  i1=1
ccmr  i2=0
ccmr  do n=1,nodes
ccmr     j=0
ccmr     i2=i2+lats_nodes_r(n)
ccmr     do i=i1,i2
ccmr        j=j+1
ccmr        global_lats_r_sn(i)=global_lats_r(i2+1-j)
ccmr     enddo
ccmr     i1=i1+lats_nodes_r(n)
ccmr  enddo

!-----------------------------------------------------------------------

           iprint =  0
      if ( iprint.eq.1 .and. me.eq.0 .and. i_count.eq.1 ) then

         j=0
         do n=1,nodes
            if (lats_nodes_a(n).gt.0) then
               write(me+8500,'(" ")')
               do jj=1,lats_nodes_a(n)
                  j=j+1
                  write(me+8500,'(2i4,
     &                      "     global_lats_a("   ,i4,")=",i4,
     &                      "     global_lats_a_sn(",i4,")=",i4)')
     &                  n, jj, j, global_lats_a(j),
     &                         j, global_lats_a_sn(j)
               enddo
            endif
         enddo

ccmr           write(me+8500,'(" ")')
ccmr           write(me+8500,'("------------------------------")')
ccmr           write(me+8500,'(" ")')
ccmr     j=0
ccmr     do n=1,nodes
ccmr        if (lats_nodes_r(n).gt.0) then
ccmr           write(me+8500,'(" ")')
ccmr           do jj=1,lats_nodes_r(n)
ccmr              j=j+1
ccmr              write(me+8500,'(2i4,
ccmr &                      "     global_lats_r("   ,i4,")=",i4,
ccmr &                      "     global_lats_r_sn(",i4,")=",i4)')
ccmr &                  n, jj, j, global_lats_r(j),
ccmr &                         j, global_lats_r_sn(j)
ccmr           enddo
ccmr        endif
ccmr     enddo

      close(me+8500)

      endif

!-----------------------------------------------------------------------

cc    processor copies or receives values for its own latitudes
      num_copy=0  !  set number copy    count to zero
      num_recv=0  !  set number receive count to zero
      do lan_a=1,lats_node_a
cc       lat is latitude number from 1 to latg
         lat = global_lats_a_sn(ipt_lats_node_a-1+lan_a)
         i=0  !  set pointer for all latitudes to zero
         do n=1,nodes  !  increment through all computational tasks
            do lan_b=1,lats_nodes_r(n)
               i=i+1  !  increment pointer for all latitudes
               if (lat .eq. global_lats_r(i)) then
                  if (n-1.eq.me) then
cc                   values are in same processor
                     num_copy=num_copy+1  !  increment number copy    count
                     do ntr=1,ntrac
                       do k=1,levs
!mjr                    do i=1,min(lonf,lons_lat)
                         ktem = k_pointer-1+k+(ntr-1)*levs
                         do i=1,lonf
                            rgt_h(xhalo+i,levs+1-k,lan_a+yhalo,ntr) =
     &                          bak_gr_r_2(i,ktem,lan_b)
                         enddo
                       enddo
                     enddo
                  else
cc                   values are in different processors
                     num_recv=num_recv+1  !  increment number receive count
                     lan_a_recv(num_recv)=lan_a
                     call mpi_irecv(  !  nonblocking receive operation
     &                    grid_recv(1,1,1,num_recv),  !  receive buffer
!mjr &                    lonrx*levs*3,  !  number of elements
     &                    lonr *levs*ntrac,  !  number of elements
     &                    mpi_r_def,  !  data type
     &                    n-1,  !  rank of source process
!Moor&                    1*10**8+(n-1)*10**4+lat,  !  message tag
     &                    1*itaga+(n-1)*itagb+lat,  !  message tag
     &                    mc_comp,  !  communicator handle
     &                    req_recv(num_recv),  !  communication request
     &                    ierr)  !  out: return code
                  endif
                  go to 12345
               endif
            enddo ! do lan_b
         enddo ! do n
12345    continue
      enddo ! do lan_a

!-----------------------------------------------------------------------

cc    processor sends values of its own latitudes
cc    to other processors' latitudes
      num_send=0  !  set number send count to zero
      i=0  !  set pointer for all latitudes to zero
      do n=1,nodes  !  increment through all computational tasks
         do lan_a=1,lats_nodes_a(n)
            i=i+1  !  increment pointer for all latitudes
            do lan_b=1,lats_node_r
cc             lat is latitude number from 1 to latr
               lat = global_lats_r(ipt_lats_node_r-1+lan_b)
               if (lat .eq. global_lats_a_sn(i)) then
                  if (n-1.ne.me) then
                     num_send=num_send+1  !  increment number send count
                     call mpi_isend(  !  nonblocking send operation
     &                    bak_gr_r_2(1,k_pointer,lan_b),  !  send buffer
!mjr &                    lonrx*levs*3,  !  number of elements
     &                    lonr *levs*ntrac,  !  number of elements
     &                    mpi_r_def,  !  data type
     &                    n-1,  !  rank of destination process
!Moor&                    1*10**8+me*10**4+lat,  !  message tag
     &                    1*itaga+me*itagb+lat,  !  message tag
     &                    mc_comp,  !  communicator handle
     &                    req_send(num_send),  !  communication request
     &                    ierr)  !  out: return code
                  endif
                  go to 23456
               endif
            enddo ! do lan_b
23456       continue
         enddo ! do lan_a
      enddo ! do n

!-----------------------------------------------------------------------

ccmr  call mpi_barrier (mc_comp,ierr)

      do j=1,num_recv  !  increment through number receive count
cc       waits for nonblocking receive operation to complete
         call mpi_wait(
     &        req_recv(j),  !  communication request to wait for
     &        stat,           !  out: status object
     &        ierr)           !  out: return code
      enddo

      do j=1,num_send  !  increment through number send count
cc       waits for nonblocking send operation to complete
         call mpi_wait(
     &        req_send(j),  !  communication request to wait for
     &        stat,           !  out: status object
     &        ierr)           !  out: return code
      enddo

!-----------------------------------------------------------------------

      do lan_a=1,num_recv
         ll = lan_a_recv(lan_a) + yhalo
         do ntr=1,ntrac
           do k=1,levs
             kk = levs+1-k
!mjr         do i=1,min(lonf,lons_lat)
             do i=1,lonf
               rgt_h(xhalo+i,kk,ll,ntr) = grid_recv(i,k,ntr,lan_a)
             enddo
           enddo
         enddo
      enddo

!-----------------------------------------------------------------------

      if ( i_print .eq. 1 ) then

         print*,' '
         print*,' mpi_tracers_b_to_a ',
     &          ' k_pointer=', k_pointer,
     &          ' i_count=',   i_count,
     &          ' num_copy=',  num_copy,
     &          ' num_recv=',  num_recv,
     &          ' num_send=',  num_send
         print*,' '

         do lan_a=1,lats_node_a
            lat = global_lats_a(ipt_lats_node_a-1+lan_a)
            write(*,'("lan_a=",i5," lat=",i5,"      rg1_h=",4e22.12)')
     &           lan_a, lat,
     &           rgt_h(xhalo+1,     1,lan_a+yhalo,1),
     &           rgt_h(xhalo+1,     2,lan_a+yhalo,1),
     &           rgt_h(xhalo+1,levs-1,lan_a+yhalo,1),
     &           rgt_h(xhalo+1,levs  ,lan_a+yhalo,1)
         enddo
         print*,' '

         do lan_b=1,lats_node_r
            lat = global_lats_r(ipt_lats_node_r-1+lan_b)
            write(*,'("lan_b=",i5," lat=",i5," bak_gr_r_2=",4e22.12)')
     &           lan_b, lat,
     &           bak_gr_r_2(1,k_pointer-1+     1,lan_b),
     &           bak_gr_r_2(1,k_pointer-1+     2,lan_b),
     &           bak_gr_r_2(1,k_pointer-1+levs-1,lan_b),
     &           bak_gr_r_2(1,k_pointer-1+levs  ,lan_b)
         enddo
         print*,' '

      endif

!-----------------------------------------------------------------------

      return
      end
