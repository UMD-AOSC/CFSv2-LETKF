      subroutine set_halos(grid_h_2,
     &                     lats_nodes_h,global_lats_h_ns,
     &                     lonsperlat,
     &                     lon_dim_h,xhalo,yhalo,
     &                     lot,lot_dim)
cc
cc    In this version, latitudes of grid_h_2 are south to north.
cc    Subroutine set_halos sets x and y halos
cc    for real array grid_h_2
cc    dimensioned:   grid_h_2(lon_dim_h,lot_dim,lats_dim_h)
cc
cc    integer lats_nodes_h(nodes) - 
cc    number of latitudes for each processor.
cc
cc    integer global_lats_h_ns(latg+2*yhalo*nodes) -
cc    latitude numbers for all processors, north to south.
cc
cc    integer lonsperlat(latg) -
cc    number of longitudes for each latitude number.
cc
cc    integer xhalo -
cc    length of east-west halo.
cc
cc    integer yhalo -
cc    length of north-south halo.
cc
cc    integer lot -
cc    number of sets of data; lot .le. lot_dim
cc
      use machine    , only : kind_evod
      use resol_def  , only : latg,lonf
      use layout1    , only : me,nodes
      use mpi_def    , only : mc_comp,mpi_r_def,stat
      use layout_lag , only : ipt_lats_node_h,lats_dim_h,lats_node_h
cc
      implicit none
cc
      integer xhalo, yhalo, lot, lot_dim
      real(kind=kind_evod) grid_h_2(lon_dim_h,lot_dim,lats_dim_h)
cc
      integer lats_nodes_h(nodes)
      integer global_lats_h_ns(latg+2*yhalo*nodes)
      integer global_lats_h_sn(latg+2*yhalo*nodes)
      integer lonsperlat(latg)
      integer lon_dim_h
cc
      integer i, i1, i2, ierr, iprint, j, jj, k, lan, lan_n, lat, n 
     &,       lons_lat, num_recv, num_send
!mjr  integer lon_dim
cc
      integer req_recv(2*yhalo*nodes)
      integer req_send(2*yhalo*nodes)
cc
      integer, parameter :: itaga=2**5, itagb=2**4
      integer i_count
      save    i_count
      data    i_count / 0 /
              i_count = i_count + 1
cc
      if ( xhalo .eq. 0 .and. yhalo .eq. 0 ) return
cc
ccmr---------------------------------------------------
cc
      i1=1
      i2=0
      do n=1,nodes
         j=0
         i2=i2+lats_nodes_h(n)
         do i=i1,i2
            j=j+1
            global_lats_h_sn(i)=global_lats_h_ns(i2+1-j)
         enddo
         i1=i1+lats_nodes_h(n)
      enddo
cc
ccmr---------------------------------------------------
cc
           iprint =  1
      if ( iprint.eq.1 .and. me.eq.0 .and. i_count.eq.1 ) then
cc
         write(me+6500,'("set_halos south north   yhalo=",i3,
     &                 "   nodes=",i3/)')
     &                 yhalo,nodes
cc
         do n=1,nodes
            write(me+6500,'("lats_nodes_h(",i4,")=",i4)')
     &                    n, lats_nodes_h(n)
         enddo
cc
         j=0
         do n=1,nodes
            if (lats_nodes_h(n).gt.0) then
               write(me+6500,'(" ")')
               do jj=1,lats_nodes_h(n)
                  j=j+1
                  write(me+6500,'(2i4,
     &                      "     global_lats_h_ns(",i4,")=",i4,
     &                      "     global_lats_h_sn(",i4,")=",i4)')
     &                  n, jj, j, global_lats_h_ns(j),
     &                         j, global_lats_h_sn(j)
               enddo
            endif
         enddo
cc
      endif
cc
      close(me+6500)
cc
ccmr---------------------------------------------------
cc
      do lan=1+yhalo,lats_node_h-yhalo

         lat = global_lats_h_sn(ipt_lats_node_h-1+lan)
!mjr     lon_dim = lon_dims_h_grid(lan)
         lons_lat = lonsperlat(lat)

         do k=1,lot
!           initialize xhalo east endpoints
            i2=0
            do i=xhalo,1,-1
               grid_h_2(i                ,k,lan)=
     &         grid_h_2(xhalo+lons_lat-i2,k,lan)
               i2=i2+1
            enddo
!           initialize xhalo west endpoints
!mjr        do i=1,xhalo+(lon_dim  -lonf-2*xhalo)
            do i=1,xhalo+(lon_dim_h-lonf-2*xhalo)
               grid_h_2(xhalo+i+lons_lat,k,lan)=
     &         grid_h_2(xhalo+i         ,k,lan) 
            enddo
         enddo

!     fin initialization of xhalo endpoints
      enddo
      if(nodes.eq.1)return
cc
ccmr---------------------------------------------------
cc
cc    each processor has lats_node_h latitudes:
cc    yhalo south halo latitudes;
cc    lats_node_h-2*yhalo non-halo latitudes;
cc    yhalo north halo latitudes;
cc
cc    copy or receive for south halo and north halo;
cc    processor sets values for its own halo latitudes
      num_recv=0  !  set number receive count to zero
      do k=1,2  !  k=1 south halo ;  k=2 north halo 
         do j=1,yhalo  !  increment south to north in halo
            if (k.eq.1) lan=j  !  lan is pointer in processor's lats
            if (k.eq.2) lan=lats_node_h-yhalo+j
cc          lat is latitude number from 1 to latg
            lat = global_lats_h_sn(ipt_lats_node_h-1+lan)
cc          lon_dim is dimension of longitude circle
!mjr        lon_dim = lon_dims_h_grid(lan)
            i=0  !  set pointer for all latitudes (halo and non-halo) to zero
            do n=1,nodes  !  increment through all computational tasks
               i=i+yhalo  !  skip south halo
               do lan_n=1+yhalo,lats_nodes_h(n)-yhalo  !  non-halo lats
                  i=i+1  !  increment pointer for all latitudes
                  if (lat .eq. global_lats_h_sn(i)) then
                     if (n-1.eq.me) then
cc                      values are in same processor
                        grid_h_2(:,:,lan  )=
     &                  grid_h_2(:,:,lan_n)
                     else
cc                      values are in different processors
                        num_recv=num_recv+1  !  increment number receive count
                        call mpi_irecv(  !  nonblocking receive operation
     &                       grid_h_2(1,1,lan),  !  receive buffer
!mjr &                       lot*lon_dim,    !  number of elements
     &                       lot*lon_dim_h,  !  number of elements
     &                       mpi_r_def,  !  data type
     &                       n-1,  !  rank of source process
!Moo &                       1*10**8+(n-1)*10**4+lat,  !  message tag
     &                       1*itaga+(n-1)*itagb+lat,  !  message tag
     &                       mc_comp,  !  communicator handle
     &                       req_recv(num_recv),  !  communication request
     &                       ierr)  !  out: return code
                     endif
                     go to 12345
                  endif
               enddo
               i=i+yhalo  !  skip north halo
            enddo
12345       continue
         enddo
      enddo
cc
cc    processor sends values of its own non-halo latitudes
cc    to other processors' halo latitudes
      num_send=0  !  set number send count to zero
      i=0  !  set pointer for all latitudes (halo and non-halo) to zero
      do n=1,nodes  !  increment through all computational tasks
         do k=1,2  !  k=1 south halo ;  k=2 north halo 
            do j=1,yhalo  !  increment south to north in halo
               i=i+1  !  increment pointer for all latitudes
               do lan=1+yhalo,lats_node_h-yhalo  !  !  non-halo lats
cc                lat is latitude number from 1 to latg
                  lat = global_lats_h_sn(ipt_lats_node_h-1+lan)
                  if (lat .eq. global_lats_h_sn(i)) then
                     if (n-1.ne.me) then
cc                      lon_dim is dimension of longitude circle
!mjr                    lon_dim = lon_dims_h_grid(lan)
                        num_send=num_send+1  !  increment number send count
                        call mpi_isend(  !  nonblocking send operation
     &                       grid_h_2(1,1,lan),  !  send buffer
!mjr &                       lot*lon_dim,    !  number of elements
     &                       lot*lon_dim_h,  !  number of elements
     &                       mpi_r_def,  !  data type
     &                       n-1,  !  rank of destination process
!Moo &                       1*10**8+me*10**4+lat,  !  message tag
     &                       1*itaga+me*itagb+lat,  !  message tag
     &                       mc_comp,  !  communicator handle
     &                       req_send(num_send),  !  communication request
     &                       ierr)  !  out: return code
                     endif
                     go to 23456
                  endif
               enddo
23456          continue
            enddo
            if (k.eq.1) i=i+lats_nodes_h(n)-2*yhalo  !  skip non-halo
         enddo
      enddo
cc
ccmr  print*,' num_recv=',num_recv,' num_send=',num_send
cc
ccmr  call mpi_barrier (mc_comp,ierr)
cc
      do j=1,num_recv  !  increment through number receive count
cc       waits for nonblocking grid_h_2 receive operation to complete
         call mpi_wait(
     &        req_recv(j),  !  communication request to wait for
     &        stat,           !  out: status object
     &        ierr)           !  out: return code
      enddo
cc
      do j=1,num_send  !  increment through number send count
cc       waits for nonblocking grid_h_2 send operation to complete
         call mpi_wait(
     &        req_send(j),  !  communication request to wait for
     &        stat,           !  out: status object
     &        ierr)           !  out: return code
      enddo
cc
      return
      end
