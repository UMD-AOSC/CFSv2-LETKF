
      subroutine get_tc_ecm(factor_b2t_ref,levs,me)
CC
      use machine        , only : kind_evod
!!!!  use resol_def
      use coordinate_def , only : ak5,bk5
      IMPLICIT NONE
CC

      real(kind=kind_evod) factor_b2t_ref(levs)
      real(kind=kind_evod) factor_t2b_ref(levs)
      real(kind=kind_evod) pr_t2b_5(levs+1),pr_t2b(levs)
      real(kind=kind_evod) t_rc(levs)
      real(kind=kind_evod) rtemp1,rtemp2
                                                                                
      real(kind=kind_evod) p_rs,t_rs,alfa,t_rt
      integer k,levs,me
                                                                                
                                                                                
      real(kind=kind_evod) cons1,consalfa,cons288,cons1013p2  !constant
      real(kind=kind_evod) cons216p5, cons0p5,cons0  !constant
                                                                                
                                                                                
      cons0        = 0.d0       !constant
      cons1        = 1.d0       !constant
      consalfa     = 5.256d0    !constant
      cons288      = 288.d0     !constant
      cons1013p2   =1013.2d0    !constant
      cons0p5      =0.5d0       !constant
      cons216p5    =216.5d0     !constant
                                                                                
      p_rs=cons1013p2
      t_rs=cons288
      alfa=consalfa
      t_rt=cons216p5
!sela  print*,'me in get_tc_ecm=',me                              
!sela  print*,'begin get_tc_ecm',p_rs,t_rs,alfa,t_rt,cons0,cons1,
!sela. consalfa,cons288,cons1013p2,cons0p5,cons216p5
                                                                                
      do k=1,levs+1
       pr_t2b_5(k)= ak5(k)+bk5(k)*101.325  !ak(k) bk(k) go top to bottom
!sela       if(me.eq.0)print*,'tc_ecm  pr_t2b_5=',k, pr_t2b_5(k)
       pr_t2b_5(k)= pr_t2b_5(k)*10.          ! to make h.pascal
!sela     if(me.eq.0)print*,'tc_ecm   times 10 pr_t2b_5=',k, pr_t2b_5(k)
      enddo
                                                                                
      do  k=1,levs
       pr_t2b(k)=cons0p5*(pr_t2b_5(k)+pr_t2b_5(k+1))
!sela  if(me.eq.0)print*,'tc_ecm  pr_t2b in h.pascal =', pr_t2b(k),k
      enddo
      do  k=1,levs
       t_rc(k)=t_rs*(pr_t2b(k)/p_rs)**alfa
!sela  if(me.eq.0)print*,'tc_ecm t_rc(k)=',t_rc(k),' k=',k
      enddo
      factor_b2t_ref=-99999999.
      factor_t2b_ref=-99999999.
                                                                                
      do  k=1,levs
         rtemp1=t_rc(k)
         if( rtemp1 .gt. t_rt )then

          rtemp2=cons0p5*(ak5(k)+bk5(k+1))*alfa*t_rc(k)*p_rs
          factor_t2b_ref(k)=rtemp2/pr_t2b(k)

!sela  if(me.eq.0)print*,'tc_ecm factor_t2b_ref(k)=',factor_t2b_ref(k),
!sela. 'if( rtemp1 .gt. t_rt ) k=',k
         else

          factor_t2b_ref(k)=cons0
!sela  if(me.eq.0)print*,'tc_ecm factor_t2b_ref(k)=',factor_t2b_ref(k),
!sela. 'if( rtemp1 .le. t_rt ) k=',k
         endif    
      enddo

      do k=1,levs
       factor_b2t_ref(k)=factor_t2b_ref(levs+1-k)
!sela  if(me.eq.0)print*,'tc_ecm factor_b2t_ref(k)=',factor_b2t_ref(k),
!sela. ' k=',k
      enddo
      return
      end
