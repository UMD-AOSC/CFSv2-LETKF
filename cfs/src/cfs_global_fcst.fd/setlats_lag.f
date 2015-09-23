      subroutine setlats_lag(lats_nodes_a,global_lats_a,
     &                     lats_nodes_h,global_lats_h,yhalo)
cc
      use resol_def , only : latg
      use layout1   , only : me,nodes
!mjr  use mpi_def
      implicit none
cc
      integer              yhalo
cc
      integer              lats_nodes_a(nodes)
      integer              lats_nodes_h(nodes)
cc
      integer              global_lats_a(latg)
      integer              global_lats_h(latg+2*yhalo*nodes)
cc
      integer              iprint
cc
      integer              jj,jpt_a,jpt_h,lat_val,nn,nodes_lats
cc
      lats_nodes_h=0
cc
      nodes_lats=0
      do nn=1,nodes
         if (lats_nodes_a(nn).gt.0) then
             lats_nodes_h(nn)=lats_nodes_a(nn)+2*yhalo
             nodes_lats=nodes_lats+1
         endif
      enddo
cc
      global_lats_h=0
cc
cc    set non-yhalo latitudes
      jpt_a=0
      jpt_h=yhalo
      do nn=1,nodes
         if        (lats_nodes_a(nn).gt.0) then
            do jj=1,lats_nodes_a(nn)
               jpt_a=jpt_a+1
               jpt_h=jpt_h+1
               global_lats_h(jpt_h)=global_lats_a(jpt_a)
            enddo
            jpt_h=jpt_h+2*yhalo
         endif
      enddo
cc
cc    set north pole yhalo
      do jj=yhalo,1,-1
         global_lats_h(jj)=global_lats_a(1)+(yhalo-jj)
      enddo
cc
cc    set south pole yhalo
      do jj=1,yhalo
         global_lats_h(latg+2*yhalo*nodes_lats-yhalo+jj)=
     &   global_lats_a(latg)-(jj-1)
      enddo
cc
      if (lats_nodes_a(1).ne.latg) then
cc
cc       set non-polar south yhalos
         jpt_h=0
         do nn=1,nodes-1
            jpt_h=jpt_h+lats_nodes_h(nn)
            lat_val=global_lats_h(jpt_h-yhalo)
            do jj=1,yhalo
               global_lats_h(jpt_h-yhalo+jj)=min(lat_val+jj,latg)
            enddo
         enddo
cc
cc       set non-polar north yhalos
         jpt_h=0
         do nn=1,nodes-1
            jpt_h=jpt_h+lats_nodes_h(nn)
            lat_val=global_lats_h(jpt_h+yhalo+1)
            do jj=1,yhalo
               global_lats_h(jpt_h+yhalo-(jj-1))=max(lat_val-jj,1)
            enddo
         enddo
cc
      endif
cc

      iprint = 0
!     iprint = 1
      if (iprint.eq.1 .and.me.eq.0) then
cc
         write(me+6000,'("setlats_h   yhalo=",i3,"   nodes=",i3/)')
     &                 yhalo,nodes
cc
         do nn=1,nodes
            write(me+6000,'("lats_nodes_a(",i4,")=",i4,"   ",
     &                   "   lats_nodes_h(",i4,")=",i4)')
     &                   nn, lats_nodes_a(nn),
     &                   nn, lats_nodes_h(nn)
         enddo
cc
         jpt_a=0
         do nn=1,nodes
            if (lats_nodes_a(nn).gt.0) then
               write(me+6000,'(" ")')
               do jj=1,lats_nodes_a(nn)
                  jpt_a=jpt_a+1
                  write(me+6000,'(2i4,"   global_lats_a(",i4,")=",i4)')
     &                     nn, jj, jpt_a, global_lats_a(jpt_a)
               enddo
            endif
         enddo
cc
         jpt_h=0
         do nn=1,nodes
            if (lats_nodes_h(nn).gt.0) then
               write(me+6000,'(" ")')
               do jj=1,lats_nodes_h(nn)
                  jpt_h=jpt_h+1
                  write(me+6000,'(2i4,"   global_lats_h(",i4,")=",i4)')
     &                     nn, jj, jpt_h, global_lats_h(jpt_h)
               enddo
            endif
         enddo
cc
        close(6000+me)
      endif
!     close(6000+me)
cc
      return
      end
