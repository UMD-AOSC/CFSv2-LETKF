      subroutine four2flnio(lslag,workdim,nvarsdim,nvars,four_gr,
     x                    ls_nodes,max_ls_nodes,
     x                    lats_nodes,global_lats,lon_dims,
     x                    lats_node,ipt_lats_node,dimg,
     x                    lat1s,londi,latl,latl2,
     x                    flnev,flnod,
     x                    plnev,plnod,ls_node)
cc
!#include "f_hpm.h"
cc
      use resol_def
      use layout1
      use mpi_def
      implicit none
cc
      integer              nvarsdim,latl2,latl
      integer              nvars
      integer              workdim
      integer              id,dimg,ilon,londi
      integer              lat1s(0:jcap)
      logical              lslag
      integer              lats_node,ipt_lats_node
      integer              lon_dims(latgd)
cc
      real(kind=kind_evod) four_gr( londi*nvarsdim, workdim )
cc
cc
      integer              ls_nodes(ls_dim,nodes)
      integer                 max_ls_nodes(nodes)
      integer                   lats_nodes(nodes)
      integer              global_lats(latl+dimg)
cc
      real(kind=kind_mpi) works(2,nvars,ls_dim*workdim,nodes)
      real(kind=kind_mpi) workr(2,nvars,ls_dim*workdim,nodes)
cc
      integer                    kpts(1+jcap)
      integer                    kptr(1+jcap)
      integer              sendcounts(1+jcap)
      integer              recvcounts(1+jcap)
      integer                 sdispls(1+jcap)
cc
      integer              ierr,ilat,ipt_ls,i3
      integer              lval,ndisp,node,nvar,ifin
cc
ccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
cc
cc
cc
cc
      real(kind=kind_evod) fp(latl2,2,nvars)
      real(kind=kind_evod) fm(latl2,2,nvars)
cc
      real(kind=kind_evod) flnev(len_trie_ls,2,nvars)
      real(kind=kind_evod) flnod(len_trio_ls,2,nvars)
cc
      real(kind=kind_evod) plnev(len_trie_ls,latl2)
      real(kind=kind_evod) plnod(len_trio_ls,latl2)
cc
      integer              ls_node(ls_dim,3)
cc
cc
cc    local scalars
cc    -------------
cc
      integer              j,k,L,n
      integer              ind,lat,lat1
      integer              indev1,indev2
      integer              indod1,indod2
      integer              num_threads
      integer              nvar_thread_max
      integer              nvar_1,nvar_2
      integer              thread
      integer              ipt_wr(4,latl2,ls_max_node)
cc    statement functions
cc    -------------------
cc
      integer              indlsev,jbasev
      integer              indlsod,jbasod
cc
      include 'function2'
cc
      real(kind=kind_evod) cons0     !constant
      real(kind=kind_evod) cons1     !constant
cc
ccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
cc
      cons0 = 0.d0     !constant
      cons1 = 1.d0     !constant
cc
      call f_hpmstart(22,"1st_four2fln")
cc
!!
      if (.not.lslag) then
        ifin=lats_node
        id=0
      else
        id=1
        if (me.ne.nodes-1) then
           ifin=lats_node
        else
           ifin=lats_node-jintmx
        endif
      endif
!!
      kpts   = 0
!$omp parallel do private(node,L,lval,j,ilon,lat,nvar,ndisp)
      do node=1,nodes
         do L=1,max_ls_nodes(node)
            lval=ls_nodes(L,node)+1
            do j=1,ifin
               if (.not.lslag) then
                 ilon=lon_dims(j)
               else
                 ilon=londi
               endif
               lat=global_lats(ipt_lats_node-1+j)
               if ( min(lat,latl-lat+1) .ge. lat1s(lval-1) ) then
                  kpts(node)=kpts(node)+1
                  do nvar=1,nvars
cc
cjfe                 ndisp = (nvar-1)*lon_dims_ext(j)
                     ndisp = (nvar-1)*ilon
cc
                       works(1,nvar,kpts(node),node) =
     x               four_gr(ndisp+2*lval-1+id  ,j)
cjfe x               four_gr(ndisp+2*lval-1,j)
cc
                       works(2,nvar,kpts(node),node) =
     x               four_gr(ndisp+2*lval+id,j)
cjfe x               four_gr(ndisp+2*lval  ,j)
cc
                  enddo
               endif
            enddo
         enddo
      enddo
cc
cc
      kptr   = 0
      do L=1,ls_max_node
         ilat   = 1
         do node=1,nodes
            if (node.ne.nodes.or..not.lslag) then
               ifin=lats_nodes(node)
            else
               ifin=lats_nodes(node)-jintmx
            endif
cjfe        do j=1,lats_nodes_ext(node)
            do j=1,ifin
               lat=global_lats(ilat)
               ipt_ls=min(lat,latl-lat+1)
                  if( ipt_ls .ge. lat1s(ls_nodes(L,me+1)) ) then
                     kptr(node)=kptr(node)+1
                     if ( lat .le. latl2 ) then
                        ipt_wr(1,ipt_ls,L)=kptr(node)
                        ipt_wr(2,ipt_ls,L)=     node
                     else
                        ipt_wr(3,ipt_ls,L)=kptr(node)
                        ipt_wr(4,ipt_ls,L)=     node
                     endif
                  endif
               ilat = ilat + 1
            enddo
         enddo
      enddo
cc
cc
      do node=1,nodes
         sendcounts(node) = kpts(node) * 2 * nvars
         recvcounts(node) = kptr(node) * 2 * nvars
            sdispls(node) = (node-1) * 2*ls_dim*workdim*nvars
      end do
cc
      call f_hpmstop(22)
cc
cc
      call f_hpmstart(85,"bar_four2fln_alltoal")
cc
      call mpi_barrier (MC_COMP,ierr)
cc
      call f_hpmstop(85)
cc
cc
      call f_hpmstart(23,"2nd_four2fln_alltoal")
cc
      call mpi_alltoallv(works,sendcounts,sdispls,MPI_R_MPI,
     x                   workr,recvcounts,sdispls,MPI_R_MPI,
     x                   MC_COMP,ierr)
cc
      call f_hpmstop(23)
cc
cjfe  call mpi_barrier (MC_COMP,ierr)
cc
cc
cc
cc
ccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
cc
      CALL countperf(0,7,0.)
      num_threads=NUM_PARTHDS()
      nvar_thread_max=(nvars+num_threads-1)/num_threads
cc
cc    -------------------------------------------------
cc    compute the coefficients of the expansion
cc    in spherical harmonics of the field at each level
cc    -------------------------------------------------
cc
      call f_hpmstart(24,"4th_four2fln_fl2eov")
cc
      do j = 1, ls_max_node   ! start of j loop ########################
cc
cc
cc
              L=ls_node(j,1)
         jbasev=ls_node(j,2)
         jbasod=ls_node(j,3)
cc
         lat1 = lat1s(L)
cc
         indev1 = indlsev(L,L)
cc
         indod1 = indlsod(L+1,L)
         if (mod(L,2).eq.mod(jcap+1,2)) then
            indev2 = indlsev(jcap+1,L)
            indod2 = indlsod(jcap  ,L)
         else
            indev2 = indlsev(jcap  ,L)
            indod2 = indlsod(jcap+1,L)
         endif
!$omp parallel do shared(fm,fp,workr,ipt_wr)
!$omp+shared(plnev,plnod,flnev,flnod)
!$omp+shared(indev1,indev2,indod1,indod2,jbasev,jbasod)
!$omp+shared(j,L,lat1,nvar_thread_max)
!$omp+private(thread,k,lat,nvar_1,nvar_2)
cc
         do thread=1,num_threads   ! start of thread loop ..............
            nvar_1=(thread-1)*nvar_thread_max+1
            nvar_2=min(nvar_1+nvar_thread_max-1,nvars)
ccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
cc
            DO k = nvar_1,nvar_2
              do lat = lat1, latl2
cc
cc
cc
                fp(lat,1,k) = workr(1,k,ipt_wr(1,lat,j),ipt_wr(2,lat,j))
     x                       +workr(1,k,ipt_wr(3,lat,j),ipt_wr(4,lat,j))
cc
                fp(lat,2,k) = workr(2,k,ipt_wr(1,lat,j),ipt_wr(2,lat,j))
     x                       +workr(2,k,ipt_wr(3,lat,j),ipt_wr(4,lat,j))
cc
                fm(lat,1,k) = workr(1,k,ipt_wr(1,lat,j),ipt_wr(2,lat,j))
     x                       -workr(1,k,ipt_wr(3,lat,j),ipt_wr(4,lat,j))
cc
                fm(lat,2,k) = workr(2,k,ipt_wr(1,lat,j),ipt_wr(2,lat,j))
     x                       -workr(2,k,ipt_wr(3,lat,j),ipt_wr(4,lat,j))
cc
               enddo
            enddo
cc
ccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
            if ( kind_evod .eq. 8 ) then !------------------------------
cc
cc          compute even real      expansion coefficients
cc          compute even imaginary expansion coefficients
cc
            call dgemm ('n', 'n', indev2-indev1+1, 2*(nvar_2-nvar_1+1),
     &                  latl2-lat1+1, cons1,                 !constant
     &                  plnev(indev1,lat1), len_trie_ls,
     &                  fp(lat1,1,nvar_1), latl2, cons0,     !constant
     &                  flnev(indev1,1,nvar_1), len_trie_ls)
cc
cc          compute odd real      expansion coefficients
cc          compute odd imaginary expansion coefficients
cc
            call dgemm ('n', 'n', indod2-indod1+1, 2*(nvar_2-nvar_1+1),
     &                  latl2-lat1+1, cons1,                 !constant
     &                  plnod(indod1,lat1), len_trio_ls,
     &                  fm(lat1,1,nvar_1), latl2, cons0,     !constant
     &                  flnod(indod1,1,nvar_1), len_trio_ls)
            else !------------------------------------------------------
cc
cc          compute even real      expansion coefficients
cc          compute even imaginary expansion coefficients
cc
            call sgemm ('n', 'n', indev2-indev1+1, 2*(nvar_2-nvar_1+1),
     &                  latl2-lat1+1, cons1,                 !constant
     &                  plnev(indev1,lat1), len_trie_ls,
     &                  fp(lat1,1,nvar_1), latl2, cons0,     !constant
     &                  flnev(indev1,1,nvar_1), len_trie_ls)
cc
cc          compute odd real      expansion coefficients
cc          compute odd imaginary expansion coefficients
cc
            call sgemm ('n', 'n', indod2-indod1+1, 2*(nvar_2-nvar_1+1),
     &                  latl2-lat1+1, cons1,                 !constant
     &                  plnod(indod1,lat1), len_trio_ls,
     &                  fm(lat1,1,nvar_1), latl2, cons0,     !constant
     &                  flnod(indod1,1,nvar_1), len_trio_ls)
            endif !-----------------------------------------------------
cc
            if (mod(L,2).eq.mod(jcap+1,2)) then
cc             set the even (n-L) terms of the top row to zero
               do k = nvar_1, nvar_2
                  flnev(indev2,1,k) = cons0     !constant
                  flnev(indev2,2,k) = cons0     !constant
               end do
            else
cc             set the  odd (n-L) terms of the top row to zero
               do k = nvar_1, nvar_2
                  flnod(indod2,1,k) = cons0     !constant
                  flnod(indod2,2,k) = cons0     !constant
               end do
            endif
cc
         end do   ! end of thread loop .................................
cc
      end do   ! end of do j loop ######################################
cc
      call f_hpmstop(24)
cc
      CALL countperf(1,7,0.)
      return
      end
