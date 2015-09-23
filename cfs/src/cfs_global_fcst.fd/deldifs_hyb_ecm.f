      subroutine deldifs_hyb_ecm(rte,we,qme,xe,ye,
     x                   rto,wo,qmo,xo,yo,deltim,
     x                   ls_node,coef00,k_level)
CC
      use machine        , only : kind_evod
      use resol_def      , only : jcap,levp1,levs,ntrac
      use layout1        , only : len_trie_ls,len_trio_ls,
     &                            ls_dim,ls_max_node
      use coordinate_def , only : ak5,bk5
      IMPLICIT NONE
CC

      real(kind=kind_evod) rte(len_trie_ls,2,levs,ntrac)
      real(kind=kind_evod)  we(len_trie_ls,2)
      real(kind=kind_evod) qme(len_trie_ls,2)
      real(kind=kind_evod)  xe(len_trie_ls,2)
      real(kind=kind_evod)  ye(len_trie_ls,2)
      real(kind=kind_evod) rto(len_trio_ls,2,levs,ntrac)
      real(kind=kind_evod)  wo(len_trio_ls,2)
      real(kind=kind_evod) qmo(len_trio_ls,2)
      real(kind=kind_evod)  xo(len_trio_ls,2)
      real(kind=kind_evod)  yo(len_trio_ls,2)
      real(kind=kind_evod) deltim
      integer              ls_node(ls_dim,3)
      real(kind=kind_evod) coef00(levs,ntrac)
      integer              k_level
      real(kind=kind_evod),allocatable :: rk_b2t(:)
      real(kind=kind_evod),allocatable :: dne(:)
      real(kind=kind_evod),allocatable :: dno(:)
      real(kind=kind_evod) sf(levs),si(levp1),sl(levs)
      integer              i,it,k
      integer              l,locl,n,n0,nd
      integer              indev
      integer              indod
      integer              indev1,indev2
      integer              indod1,indod2
      real(kind=kind_evod) pr_t2b_5(levp1)
      real(kind=kind_evod) pr_t2b(levs),tau,rkzero,rk_t2b(levs)         
      real(kind=kind_evod) zdilev,rj_jp1_sq,rfact,fact,realval                  
      real(kind=kind_evod) rjunk1,rjunk2                                        
      integer ist,lev
       
      real(kind=kind_evod) cons0,cons1,cons2      !constant
      integer              indlsev,jbasev
      integer              indlsod,jbasod
      include 'function2'
      save      dne , dno ,rk_b2t
      cons0     = 0.d0       !constant
      cons1     = 1.d0       !constant
      cons2     = 2.d0       !constant


      if(k_level.eq.0) then
         call countperf(0,15,0.)
       tau=-999999.
       if(jcap.eq.62)  tau=2700.
       if(jcap.eq.170) tau=5400.
       if(jcap.eq.254) tau=2500.
!sela  if(jcap.eq.382) tau=2500.  ! this is provisional !!!!!!!!!!!!!!$$$$$$$$$$$$
!sela  if(jcap.eq.382) tau=5000.  ! this is provisional !!!!!!!!!!!!!!$$$$$$$$$$$$
       if(jcap.eq.382) tau=4000.  ! this is provisional !!!!!!!!!!!!!!$$$$$$$$$$$$

       rj_jp1_sq=jcap*(jcap+1)
       rj_jp1_sq=rj_jp1_sq**2             

        print*,' deltim in deldif_hyb_ecm=',deltim
       rkzero=cons2*deltim/(tau*rj_jp1_sq)
       
100    format('rkzero=',e13.4,' jcap=',i4,
     . ' tau=',e12.5,' rj_jp1_sq=',e13.4)

      do k=1,levp1
       pr_t2b_5(k)= ak5(k)+bk5(k)*101.325  !ak(k) bk(k) go top to bottom
!sela  if(me.eq.0)print*,' pr_t2b_5=',k, pr_t2b_5(k)
      enddo
       
      do  k=1,levs
       pr_t2b(k)=0.5*(pr_t2b_5(k)+pr_t2b_5(k+1))
!sela  if(me.eq.0)print*,' pr_t2b=', pr_t2b(k),k
      enddo

         allocate ( rk_b2t(levs)     )
         allocate ( dne(len_trie_ls) )
         allocate ( dno(len_trio_ls) )


        if(pr_t2b(1).gt.0.0105)then
         do lev=1,levs
          ist=min(9,levs)
          rjunk1=1.3**(max(0,ist-lev))
          rjunk2=min(16.,rjunk1)        
          rk_t2b(lev)=rkzero*rjunk2
          rk_b2t(levp1-lev)=rk_t2b(lev)

!sela     if (me.eq.0)print*,'pr_t2b(first)=',pr_t2b(lev),' lev=',lev,
!sela.   ' ist=',ist,' rjunk1=',rjunk1,' rjunk2=',rjunk2,
!sela.   ' first branch  rk_t2b(lev)=',rk_t2b(lev)
         enddo
        endif
        do lev=1,levs
!sela   if (me.eq.0)print*,'rk_b2t(first)=',rk_b2t(lev),' lev=',lev
        enddo
       

        if(pr_t2b(1).lt.0.0100)then
          do lev=1,levs
           rjunk1=log10(pr_t2b(lev))
           zdilev=1.+7.5*(3.-log10(pr_t2b(lev)))
           rk_t2b(lev)=rkzero*(  min( 16.,max(1.,zdilev) )  )
           rk_b2t(levp1-lev)=rk_t2b(lev)


!sela     if (me.eq.0)print*,'pr_t2b(second)=',pr_t2b(lev),' lev=',lev,
!sela.   ' rjunk1=',rjunk1,' zdilev=',zdilev,
!sela.   ' second branch  rk_t2b(lev)=',rk_t2b(lev)
          enddo
        do lev=1,levs
!sela   if (me.eq.0)print*,'rk_b2t(second)=',rk_b2t(lev),' lev=',lev
        enddo

        endif



      do locl=1,ls_max_node
              l=ls_node(locl,1)
         jbasev=ls_node(locl,2)
         indev=indlsev(l,l)
         do n=l,jcap,2
            nd=n
            realval=nd*(nd+1)
            dne(indev)=realval**2
            indev=indev+1
         enddo
      enddo
      do locl=1,ls_max_node
              l=ls_node(locl,1)
         jbasev=ls_node(locl,2)
         if (mod(l,2).eq.mod(jcap+1,2)) then
            dne(indlsev(jcap+1,l))=cons0     !constant
         endif
      enddo


      do locl=1,ls_max_node
              l=ls_node(locl,1)
         jbasod=ls_node(locl,3)
         indod=indlsod(l+1,l)
         do n=l+1,jcap,2
            nd=n
            realval=nd*(nd+1)
            dno(indod)=realval**2
            indod=indod+1
         enddo
      enddo
      do locl=1,ls_max_node
              l=ls_node(locl,1)
         jbasod=ls_node(locl,3)
         if (mod(l,2).ne.mod(jcap+1,2)) then
            dno(indlsod(jcap+1,l))=cons0     !constant
         endif
      enddo
      call countperf(1,15,0.)
      return
      endif   !  if(k_level.eq.0) 


      call countperf(0,13,0.)
      k=k_level
         do locl=1,ls_max_node
                 l=ls_node(locl,1)
            jbasev=ls_node(locl,2)
            if (l.eq.0) then
                              n0=2
                        else
                              n0=l
            endif
            indev1 = indlsev(n0,l)
            if (mod(l,2).eq.mod(jcap+1,2)) then
               indev2 = indlsev(jcap+1,l)
            else
               indev2 = indlsev(jcap  ,l)
            endif
            do indev = indev1 , indev2
       
               fact=dne(indev)*rk_b2t(k)
               rfact  =cons1/(cons1+fact)                 !constant
       
               we(indev,1)=we(indev,1)*rfact
               we(indev,2)=we(indev,2)*rfact
       
               xe(indev,1)=xe(indev,1)*rfact
               xe(indev,2)=xe(indev,2)*rfact
       
               ye(indev,1)=ye(indev,1)*rfact
               ye(indev,2)=ye(indev,2)*rfact
       
               rte(indev,1,1,1)=rte(indev,1,1,1)*rfact
               rte(indev,2,1,1)=rte(indev,2,1,1)*rfact
       
            enddo
         enddo
         do locl=1,ls_max_node
                 l=ls_node(locl,1)
            jbasod=ls_node(locl,3)
            indod1 = indlsod(l+1,l)
            if (mod(l,2).eq.mod(jcap+1,2)) then
               indod2 = indlsod(jcap  ,l)
            else
               indod2 = indlsod(jcap+1,l)
            endif
            do indod = indod1 , indod2
       
               fact=dno(indod)*rk_b2t(k)
               rfact  =cons1/(cons1+fact)                 !constant
       
               wo(indod,1)=wo(indod,1)*rfact
               wo(indod,2)=wo(indod,2)*rfact
       
               xo(indod,1)=xo(indod,1)*rfact
               xo(indod,2)=xo(indod,2)*rfact
       
               yo(indod,1)=yo(indod,1)*rfact
               yo(indod,2)=yo(indod,2)*rfact
       
               rto(indod,1,1,1)=rto(indod,1,1,1)*rfact
               rto(indod,2,1,1)=rto(indod,2,1,1)*rfact
       
            enddo
         enddo
      do 90 it=2,ntrac
            do locl=1,ls_max_node
                    l=ls_node(locl,1)
               jbasev=ls_node(locl,2)
               if (l.eq.0) then
                                 n0=2
                           else
                                 n0=l
               endif
               indev1 = indlsev(n0,l)
               if (mod(l,2).eq.mod(jcap+1,2)) then
                  indev2 = indlsev(jcap+1,l)
               else
                  indev2 = indlsev(jcap  ,l)
               endif
               do indev = indev1 , indev2
       
               fact=dne(indev)*rk_b2t(k)
               rfact  =cons1/(cons1+fact)                 !constant
       
               rte(indev,1,1,it)= rte(indev,1,1,it)*rfact
               rte(indev,2,1,it)= rte(indev,2,1,it)*rfact
       
               enddo
            enddo
            do locl=1,ls_max_node
                    l=ls_node(locl,1)
               jbasod=ls_node(locl,3)
               indod1 = indlsod(l+1,l)
               if (mod(l,2).eq.mod(jcap+1,2)) then
                  indod2 = indlsod(jcap  ,l)
               else
                  indod2 = indlsod(jcap+1,l)
               endif
               do indod = indod1 , indod2
       
               fact=dno(indod)*rk_b2t(k)
               rfact  =cons1/(cons1+fact)                 !constant
       
               rto(indod,1,1,it)= rto(indod,1,1,it)*rfact
               rto(indod,2,1,it)= rto(indod,2,1,it)*rfact
       
               enddo
            enddo
   90 continue
      call countperf(1,13,0.)
      return
      end
