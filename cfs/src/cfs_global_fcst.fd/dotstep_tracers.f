      SUBROUTINE do_tstep(deltim,kdt,PHOUR,
     &                 TRIE_LS,TRIO_LS,
     &                 LS_NODE,LS_NODES,MAX_LS_NODES,
     &                 LATS_NODES_A,GLOBAL_LATS_A,
     &                 LONSPERLAT,
     &                 LATS_NODES_R,GLOBAL_LATS_R,
     &                 LONSPERLAR,
!    &                 LATS_NODES_EXT,GLOBAL_LATS_EXT,
     &                 EPSE,EPSO,EPSEDN,EPSODN,
     &                 SNNP1EV,SNNP1OD,NDEXEV,NDEXOD,
     &                 PLNEV_A,PLNOD_A,PDDEV_A,PDDOD_A,
     &                 PLNEW_A,PLNOW_A,
     &                 PLNEV_R,PLNOD_R,PDDEV_R,PDDOD_R,
     &                 PLNEW_R,PLNOW_R,
     &                 SYN_LS_A,DYN_LS_A,
!    &                 SYN_GR_A_1,DYN_GR_A_1,ANL_GR_A_1,
!    &                 SYN_GR_A_2,DYN_GR_A_2,ANL_GR_A_2,
     &                 XLON,XLAT,COSZDG, sfc_fld, flx_fld, nst_fld,
     &                 HPRIME,SWH,HLW,FLUXR,SFALB,SLAG,SDEC,CDEC,
     &                 OZPLIN,JINDX1,JINDX2,DDY,PDRYINI,
     &                 phy_f3d,  phy_f2d,
     &                 ZHOUR,N1,N4,LSOUT,COLAT1,CFHOUR1,SPS,fscav)
!!
!#include "f_hpm.h"
      use machine             , only : kind_evod,kind_phys,kind_rad
      use resol_def           , only : latg,latg2,latr,latr2,levh,levs,
     &                                 lonr,lotd,lots,lsoil,nfxr,nmtvr,
     &                                 ntoz,ntrac,ncld,num_p2d,num_p3d,
     &                                 p_di,p_dim,p_q,p_qm,p_rm,p_rq,
     &                                 p_rt,p_te,p_tem,p_uln,p_vln,
     &                                 p_w,p_x,p_y,p_ze,p_zem,p_zq,lonf
      use layout1             , only : ipt_lats_node_r,lats_node_r,
     &                                 len_trie_ls,len_trio_ls,
     &                                 ls_dim,ls_max_node,
     &                                 me,me_l_0,nodes,lats_dim_a,
     .                                 ipt_lats_node_a,lats_node_a
      use vert_def            , only : am,bm,si,sl,sv,tov
      use date_def            , only : fhour,idate,shour,spdmax
      use namelist_def        , only : adiab,ens_nam,fhcyc,filta,
     &                                 gen_coord_hybrid,gg_tracers,
     &                                 hybrid, igen,explicit,mom4ice,
     &                                 ldiag3d,lsfwd,lslwr,lsswr,
     &                                 lggfs3d,fhgoc3d,ialb,nst_fcst,
     &                                 ngptc,nscyc,nsres,nszer,semilag,
     &                                 sl_epsln,nsout
      use mpi_def             , only : icolor,kind_mpi,liope,
     &                                 mc_comp,mpi_r_mpi,comp_task
      use ozne_def            , only : latsozp,levozp,pl_coeff,timeoz

!     use layout_grid_tracers , only : rgt_a


      use Sfc_Flx_ESMFMod
      use Nst_Var_ESMFMod
      use d3d_def

      use atm_cc              , only : Coupler_id

      IMPLICIT NONE
!!     
      TYPE(Sfc_Var_Data)        :: sfc_fld
      TYPE(Flx_Var_Data)        :: flx_fld
      TYPE(Nst_Var_Data)        :: nst_fld
      integer lat
      CHARACTER(16)             :: CFHOUR1
      INTEGER,INTENT(IN)        :: LONSPERLAT(LATG),N1,N4
!!     
      REAL(KIND=KIND_EVOD),INTENT(IN)    :: deltim,PHOUR
      REAL(KIND=KIND_EVOD),INTENT(INOUT) :: ZHOUR

      integer ifirst
      data ifirst /1/
      save ifirst
!
      real, allocatable   :: gzie_ln(:,:),gzio_ln(:,:),factor_b2t_ref(:)
      save gzie_ln,gzio_ln,factor_b2t_ref

      REAL(KIND=KIND_EVOD) TRIE_LS(LEN_TRIE_LS,2,11*LEVS+3*LEVH+6)
      REAL(KIND=KIND_EVOD) TRIO_LS(LEN_TRIO_LS,2,11*LEVS+3*LEVH+6)
!!
      integer              ls_node(ls_dim,3)
!!
!     ls_node(1,1) ... ls_node(ls_max_node,1) : values of L
!     ls_node(1,2) ... ls_node(ls_max_node,2) : values of jbasev
!     ls_node(1,3) ... ls_node(ls_max_node,3) : values of jbasod
!!
      INTEGER              LS_NODES(LS_DIM,NODES)
      INTEGER          MAX_LS_NODES(NODES)
      INTEGER               LATS_NODES_A(NODES)
!     INTEGER               LATS_NODES_EXT(NODES)
      INTEGER              GLOBAL_LATS_A(LATG)
!     INTEGER        GLOBAL_LATS_EXT(LATG+2*JINTMX+2*NYPT*(NODES-1))
      INTEGER               LATS_NODES_R(NODES)
      INTEGER              GLOBAL_LATS_R(LATR)
      INTEGER                 LONSPERLAR(LATR)
!
      integer               lats_nodes_r_old(nodes)
      integer              global_lats_r_old(latr)
      logical ifshuff
!
      real(kind=kind_evod) colat1
      REAL(KIND=KIND_EVOD)    EPSE(LEN_TRIE_LS)
      REAL(KIND=KIND_EVOD)    EPSO(LEN_TRIO_LS)
      REAL(KIND=KIND_EVOD)  EPSEDN(LEN_TRIE_LS)
      REAL(KIND=KIND_EVOD)  EPSODN(LEN_TRIO_LS)
      REAL(KIND=KIND_EVOD) SNNP1EV(LEN_TRIE_LS)
      REAL(KIND=KIND_EVOD) SNNP1OD(LEN_TRIO_LS)
      INTEGER               NDEXEV(LEN_TRIE_LS)
      INTEGER               NDEXOD(LEN_TRIO_LS)
      REAL(KIND=KIND_EVOD)   PLNEV_A(LEN_TRIE_LS,LATG2)
      REAL(KIND=KIND_EVOD)   PLNOD_A(LEN_TRIO_LS,LATG2)
      REAL(KIND=KIND_EVOD)   PDDEV_A(LEN_TRIE_LS,LATG2)
      REAL(KIND=KIND_EVOD)   PDDOD_A(LEN_TRIO_LS,LATG2)
      REAL(KIND=KIND_EVOD)   PLNEW_A(LEN_TRIE_LS,LATG2)
      REAL(KIND=KIND_EVOD)   PLNOW_A(LEN_TRIO_LS,LATG2)
      REAL(KIND=KIND_EVOD)   PLNEV_R(LEN_TRIE_LS,LATR2)
      REAL(KIND=KIND_EVOD)   PLNOD_R(LEN_TRIO_LS,LATR2)
      REAL(KIND=KIND_EVOD)   PDDEV_R(LEN_TRIE_LS,LATR2)
      REAL(KIND=KIND_EVOD)   PDDOD_R(LEN_TRIO_LS,LATR2)
      REAL(KIND=KIND_EVOD)   PLNEW_R(LEN_TRIE_LS,LATR2)
      REAL(KIND=KIND_EVOD)   PLNOW_R(LEN_TRIO_LS,LATR2)
      REAL(KIND=KIND_EVOD) SYN_LS_A(4*LS_DIM,LOTS,LATG2)
      REAL(KIND=KIND_EVOD) DYN_LS_A(4*LS_DIM,LOTD,LATG2)

!     REAL(KIND=KIND_EVOD) SYN_GR_A_1(LONFX*LOTS,LATS_DIM_EXT)
!     REAL(KIND=KIND_EVOD) DYN_GR_A_1(LONFX*LOTD,LATS_DIM_EXT)
!     REAL(KIND=KIND_EVOD) ANL_GR_A_1(LONFX*LOTA,LATS_DIM_EXT)
!     REAL(KIND=KIND_EVOD) SYN_GR_A_2(LONFX*LOTS,LATS_DIM_EXT)
!     REAL(KIND=KIND_EVOD) DYN_GR_A_2(LONFX*LOTD,LATS_DIM_EXT)
!     REAL(KIND=KIND_EVOD) ANL_GR_A_2(LONFX*LOTA,LATS_DIM_EXT)
!!     
      REAL (KIND=KIND_RAD) XLON(LONR,LATS_NODE_R),
     &                     XLAT(LONR,LATS_NODE_R),
     &                     COSZDG(LONR,LATS_NODE_R),
     &                     HPRIME(LONR,NMTVR,LATS_NODE_R),
     &                     FLUXR(LONR,nfxr,LATS_NODE_R),
     &                     SFALB(LONR,LATS_NODE_R),
     &                     SWH(LONR,LEVS,LATS_NODE_R),
     &                     HLW(LONR,LEVS,LATS_NODE_R)

      REAL (kind=kind_phys)
     &     phy_f3d(LONR,LEVS,num_p3d,lats_node_r),
     &     phy_f2d(lonr,num_p2d,lats_node_r),
!
     &     DDY(LATS_NODE_R), fscav(ntrac-ncld-1)

      INTEGER JINDX1(LATS_NODE_R),JINDX2(LATS_NODE_R)
!!     
      INTEGER LEV,LEVMAX
      REAL OZPLIN(LATSOZP,LEVOZP,pl_coeff,timeoz) !OZONE PL Coeff
      REAL (KIND=KIND_PHYS) PDRYINI
      REAL(KIND=KIND_EVOD) SLAG,SDEC,CDEC
!
!****************************************************************************
!$$$      INTEGER   P_GZ,P_ZEM,P_DIM,P_TEM,P_RM,P_QM
!$$$      INTEGER   P_ZE,P_DI,P_TE,P_RQ,P_Q,P_DLAM,P_DPHI,P_ULN,P_VLN
!$$$      INTEGER   P_W,P_X,P_Y,P_RT,P_ZQ
!$$$      PARAMETER(P_GZ  = 0*LEVS+0*LEVH+1,  !      GZE/O(LNTE/OD,2),
!$$$     X          P_ZEM = 0*LEVS+0*LEVH+2,  !     ZEME/O(LNTE/OD,2,LEVS),
!$$$     X          P_DIM = 1*LEVS+0*LEVH+2,  !     DIME/O(LNTE/OD,2,LEVS),
!$$$     X          P_TEM = 2*LEVS+0*LEVH+2,  !     TEME/O(LNTE/OD,2,LEVS),
!$$$     X          P_RM  = 3*LEVS+0*LEVH+2,  !      RME/O(LNTE/OD,2,LEVH),
!$$$     X          P_QM  = 3*LEVS+1*LEVH+2,  !      QME/O(LNTE/OD,2),
!$$$     X          P_ZE  = 3*LEVS+1*LEVH+3,  !      ZEE/O(LNTE/OD,2,LEVS),
!$$$     X          P_DI  = 4*LEVS+1*LEVH+3,  !      DIE/O(LNTE/OD,2,LEVS),
!$$$     X          P_TE  = 5*LEVS+1*LEVH+3,  !      TEE/O(LNTE/OD,2,LEVS),
!$$$     X          P_RQ  = 6*LEVS+1*LEVH+3,  !      RQE/O(LNTE/OD,2,LEVH),
!$$$     X          P_Q   = 6*LEVS+2*LEVH+3,  !       QE/O(LNTE/OD,2),
!$$$     X          P_DLAM= 6*LEVS+2*LEVH+4,  !  DPDLAME/O(LNTE/OD,2),
!$$$     X          P_DPHI= 6*LEVS+2*LEVH+5,  !  DPDPHIE/O(LNTE/OD,2),
!$$$     X          P_ULN = 6*LEVS+2*LEVH+6,  !     ULNE/O(LNTE/OD,2,LEVS),
!$$$     X          P_VLN = 7*LEVS+2*LEVH+6,  !     VLNE/O(LNTE/OD,2,LEVS),
!$$$     X          P_W   = 8*LEVS+2*LEVH+6,  !       WE/O(LNTE/OD,2,LEVS),
!$$$     X          P_X   = 9*LEVS+2*LEVH+6,  !       XE/O(LNTE/OD,2,LEVS),
!$$$     X          P_Y   =10*LEVS+2*LEVH+6,  !       YE/O(LNTE/OD,2,LEVS),
!$$$     X          P_RT  =11*LEVS+2*LEVH+6,  !      RTE/O(LNTE/OD,2,LEVH),
!$$$     X          P_ZQ  =11*LEVS+3*LEVH+6)  !      ZQE/O(LNTE/OD,2)
!****************************************************************************
      INTEGER   kdt,       IERR,J,K,L,LOCL,N
      real(kind=kind_evod)  batah
      REAL(KIND=kind_mpi)  coef00m(LEVS,ntrac)! temp. ozone clwater  
      REAL(KIND=kind_evod) coef00(LEVS,ntrac) ! temp. ozone clwater  
      INTEGER              INDLSEV,JBASEV
      INTEGER              INDLSOD,JBASOD
      integer iprint
      include 'function2'

      LOGICAL LSOUT, SPS, wrt_g3d
!
      real, PARAMETER:: RLAPSE=0.65E-2, omz1=10.0
!
! timings
      real(kind=kind_evod) global_times_a(latg,nodes)
     &,                    global_times_b(latr,nodes)
     &,                    global_times_r(latr,nodes)
      integer              tag,ireq1,ireq2,i
      real*8               rtc ,timer1,timer2,dt_warm, tem1, tem2
!
!     if(ifirst == 1)then
!      allocate ( factor_b2t_ref(levs), gzie_ln(len_trie_ls,2),
!    &            gzio_ln(len_trio_ls,2) )
!      ifirst=0
!     endif
!
      SHOUR = SHOUR + deltim

!-> Coupling insertion
      if (comp_task) then
        call ATM_DBG2(kdt,PHOUR,ZHOUR,SHOUR,3)
        CALL ATM_TSTEP_INIT(kdt)
      endif
!<- Coupling insertion

      if (comp_task) then
!
!     print *,' in do tstep SEMILAG=',semilag,' kdt=',kdt

        if (semilag) then    ! Joe Sela's Semi-Lagrangian Code

!         batah = 0.
!         batah = 1.                   ! Commented by Moorthi 11/23/2010
          batah = 1.0 + sl_epsln       ! Moorthi

          if(ifirst == 1) then
            allocate ( factor_b2t_ref(levs), gzie_ln(len_trie_ls,2),
     &                 gzio_ln(len_trio_ls,2) )
            call get_cd_hyb_slg(deltim,batah)

            k = 0
            CALL deldifs(
     .                TRIE_LS(1,1,P_RT+k-1), TRIE_LS(1,1,P_W+k-1),
     X                TRIE_LS(1,1,P_QM    ), TRIE_LS(1,1,P_X+k-1),
     X                TRIE_LS(1,1,P_Y +k-1), TRIE_LS(1,1,P_TEM+k-1),    ! hmhj
     X                TRIO_LS(1,1,P_RT+k-1), TRIO_LS(1,1,P_W+k-1),
     X                TRIO_LS(1,1,P_QM    ), TRIO_LS(1,1,P_X+k-1),
     X                TRIO_LS(1,1,P_Y +k-1), TRIO_LS(1,1,P_TEM+k-1),    ! hmhj
     X                deltim,SL,LS_NODE,coef00,k,hybrid,                ! hmhj
     &                gen_coord_hybrid)

            ifirst=0
          endif
!         if(kdt < 24) print*,'entering dotstep deltim=', deltim,
!    &                        ' kdt=',kdt
          global_times_a = 0.
          timer1         = rtc()
          call gloopa_hyb_slg
     &      (deltim,trie_ls,trio_ls,gzie_ln,gzio_ln,
     &       ls_node,ls_nodes,max_ls_nodes,
     &       lats_nodes_a,global_lats_a,
     &       lonsperlat,
     &       epse,epso,epsedn,epsodn,
     &       snnp1ev,snnp1od,ndexev,ndexod,
     &       plnev_a,plnod_a,pddev_a,pddod_a,plnew_a,plnow_a,
     &       global_times_a,kdt,batah,lsout)
          timer2 = rtc()

!         if (kdt.lt.4)then
!           print*,' gloopa timer = ',timer2-timer1,' kdt=',kdt
!         endif

          if(.not. adiab) then ! first if.not.adiab
            if (nscyc > 0 .and. mod(kdt,nscyc) == 1) then
!             if (me == 0) print*,' calling gcycle at kdt=',kdt
                CALL gcycle(me,LATS_NODE_R,LONSPERLAR,global_lats_r,
     &                      ipt_lats_node_r,idate,fhour,fhcyc,
     &                      XLON ,XLAT  , sfc_fld, ialb)
            endif
!
            if (num_p3d  == 3) then  ! Ferrier Microphysics initialization
              call init_micro(deltim,lonr,levs,num_p3d,lats_node_r,
     &                        phy_f3d(1,1,1,1),   fhour, me)
            endif
          endif              ! first if.not.adiab

!
!-> Coupling insertion

  ! lgetSSTICE_cc must be defined by this moment. It used to be an argument
  ! to ATM_GETSST, accessible here via USE SURFACE_cc. Now it is defined in
  ! ATM_TSTEP_INIT called above, and the USE is removed. (Even in the earlier
  ! version lgetSSTICE_cc did not have to be an actual argumnent, since
  ! it is in the module SURFACE_cc USEd by ATM_GETSST.)

          if (coupler_id >= 0) then
            do j = 1, lats_node_r
              do i = 1, lonr
                if (sfc_fld%slmsk(i,j) == 0 ) then
                  sfc_fld%tsea(i,j) = sfc_fld%tsea(i,j)
     &                 + (sfc_fld%oro(i,j)-sfc_fld%oro_uf(i,j))*rlapse
                endif
              enddo
            enddo
          endif

          call ATM_GETSSTICE(sfc_fld%TSEA,sfc_fld%TISFC,sfc_fld%FICE,
     &                       sfc_fld%HICE,sfc_fld%SHELEG,sfc_fld%SLMSK,
     &                       kdt)

          if (coupler_id >= 0) then
            do j = 1, lats_node_r
              do i = 1, lonr
                if (sfc_fld%slmsk(i,j) == 0 ) then
                  sfc_fld%tsea(i,j) = sfc_fld%tsea(i,j)
     &                 - (sfc_fld%oro(i,j)-sfc_fld%oro_uf(i,j))*rlapse
                endif
              enddo
            enddo
          endif

!<- Coupling insertion

!
          if (nst_fcst > 1) then                         ! update TSEA
            if (Coupler_id < 0 .or. .not. mom4ice) then  ! Standalone mode
              do j = 1, lats_node_r
                do i = 1, lonr
                  if (sfc_fld%slmsk(i,j) == 0 ) then
                    dt_warm = (nst_fld%xt(i,j)+nst_fld%xt(i,j))
     &                      /  nst_fld%xz(i,j)
                    sfc_fld%TSEA(i,j) = nst_fld%tref(i,j)
     &                   + dt_warm - nst_fld%dt_cool(i,j)
     &                   - (sfc_fld%oro(i,j)-sfc_fld%oro_uf(i,j))*rlapse
                  endif
                enddo
              enddo
            else                                         ! Coupled to MOM4 OM
              tem1 = 0.5 / omz1
              do j = 1, lats_node_r
                do i = 1, lonr
                  if (sfc_fld%slmsk(i,j) == 0 ) then
                    tem2 = 1.0 / nst_fld%xz(i,j)
                    sfc_fld%tsea(i,j) = sfc_fld%tsea(i,j)
     &                   + (sfc_fld%oro(i,j)-sfc_fld%oro_uf(i,j))*rlapse
                    dt_warm = (nst_fld%xt(i,j)+nst_fld%xt(i,j)) * tem2

                    if ( nst_fld%xz(i,j) > omz1) then
                      nst_fld%tref(i,j) = sfc_fld%tsea(i,j)
     &                 - (1.0-0.5*omz1*tem2) * dt_warm
     &                 + nst_fld%z_c(i,j)*nst_fld%dt_cool(i,j)*tem1
                    else
                     nst_fld%tref(i,j) = sfc_fld%tsea(i,j)
     &                 - (nst_fld%xz(i,j)*dt_warm
     &                 -  nst_fld%z_c(i,j)*nst_fld%dt_cool(i,j))*tem1
                    endif
                    sfc_fld%TSEA(i,j) = nst_fld%tref(i,j)
     &                  + dt_warm - nst_fld%dt_cool(i,j)
     &                  - (sfc_fld%oro(i,j)-sfc_fld%oro_uf(i,j))*rlapse
                  endif
                enddo
              enddo
            endif
          endif

          global_times_r = 0.0

          if (lsswr .or. lslwr) then         ! radiation call!
            if(.not. adiab) then             ! second  if.not.adiab
              call gloopr
     &          (trie_ls,trio_ls,
     &           ls_node,ls_nodes,max_ls_nodes,
     &           lats_nodes_a,global_lats_a,
     &           lats_nodes_r,global_lats_r,
     &           lonsperlar,
     &           epse,epso,epsedn,epsodn,
     &           snnp1ev,snnp1od,plnev_r,plnod_r,
     &           pddev_r,pddod_r,
     &           phour,
     &           xlon,xlat,coszdg,flx_fld%coszen,
     &           sfc_fld%slmsk,sfc_fld%sheleg,sfc_fld%SNCOVR,
     &           sfc_fld%SNOALB,sfc_fld%ZORL,sfc_fld%TSEA,
     &           HPRIME,SFALB,sfc_fld%ALVSF,sfc_fld%ALNSF,
     &           sfc_fld%ALVWF,sfc_fld%ALNWF,sfc_fld%FACSF,
     &           sfc_fld%FACWF,sfc_fld%CV,sfc_fld%CVT ,
     &           sfc_fld%CVB,SWH,HLW,flx_fld%SFCNSW,flx_fld%SFCDLW,
     &           sfc_fld%FICE,sfc_fld%TISFC,flx_fld%SFCDSW,
     &           flx_fld%sfcemis,                                    ! yth 4/09
     &           flx_fld%TSFLW,FLUXR,phy_f3d,SLAG,SDEC,CDEC,KDT,
     &           global_times_r)
!                if (iprint == 1) print*,' me = fin gloopr ',me
            endif                            ! second  if.not.adiab
          endif  !sswr .or. lslwr

!         if (iprint .eq. 1) print*,' me = beg gloopb ',me
!         if(kdt < 4)then
!           print*,' deltim in if(kdt.lt.4)=',deltim
!         endif

!$omp parallel do private(locl)
          do locl=1,ls_max_node
            call sicdife_hyb_slg(trie_ls(1,1,p_x  ), trie_ls(1,1,p_y ),
     x                         trie_ls(1,1,p_zq ), deltim/2.,
     x                         trie_ls(1,1,p_uln), trie_ls(1,1,p_vln),
     x                         ls_node,snnp1ev,ndexev,locl,batah)
            call sicdifo_hyb_slg(trio_ls(1,1,p_x  ), trio_ls(1,1,p_y ),
     x                         trio_ls(1,1,p_zq ), deltim/2.,
     x                         trio_ls(1,1,p_uln), trio_ls(1,1,p_vln),
     x                         ls_node,snnp1od,ndexod,locl,batah)
          enddo
          do j=1,len_trie_ls
            trie_ls(j,1,p_zq ) = trie_ls(j,1,p_zq)-gzie_ln(j,1)
            trie_ls(j,2,p_zq ) = trie_ls(j,2,p_zq)-gzie_ln(j,2)
          enddo
          do j=1,len_trio_ls
            trio_ls(j,1,p_zq ) = trio_ls(j,1,p_zq)-gzio_ln(j,1)
            trio_ls(j,2,p_zq ) = trio_ls(j,2,p_zq)-gzio_ln(j,2)
          enddo
!save n-1 values for diffusion, not really part of samilag scheme
          do j=1,len_trie_ls
            trie_ls(j,1,p_qm ) = trie_ls(j,1,p_zq)
            trie_ls(j,2,p_qm ) = trie_ls(j,2,p_zq)
          enddo
          do j=1,len_trio_ls
            trio_ls(j,1,p_qm ) = trio_ls(j,1,p_zq)
            trio_ls(j,2,p_qm ) = trio_ls(j,2,p_zq)
          enddo

          do k=1,levs
            do j=1,len_trie_ls
              trie_ls(j,1,p_tem+k-1) = trie_ls(j,1,p_y+k-1)
              trie_ls(j,2,p_tem+k-1) = trie_ls(j,2,p_y+k-1)
            enddo
          enddo
          do k=1,levs
            do j=1,len_trio_ls
              trio_ls(j,1,p_tem+k-1) = trio_ls(j,1,p_y+k-1)
              trio_ls(j,2,p_tem+k-1) = trio_ls(j,2,p_y+k-1)
            enddo
          enddo
!--------------------------------------------------------
          coef00(:,:) = 0.0
          IF ( ME .EQ. ME_L_0 ) THEN
            DO LOCL=1,LS_MAX_NODE
              l      = ls_node(locl,1)
              jbasev = ls_node(locl,2)
              IF ( L  ==  0 ) THEN
                N = 0
! 1 Corresponds to temperature,  2 corresponds to ozon, 3 to clwater
                DO K=1,LEVS
                  coef00(K,1) = TRIE_LS(INDLSEV(N,L),1,P_Y +K-1)
                  if (ntoz .gt. 1 .and.                                  ! hmhj
     &               .not. (hybrid.or.gen_coord_hybrid)) then            ! hmhj
                     coef00(K,ntoz) = TRIE_LS(INDLSEV(N,L),1,
     &                                  (ntoz-1)*levs+P_rt+K-1)
                  endif
                ENDDO
              ENDIF
            END DO
          END IF

          coef00m = coef00
          CALL MPI_BCAST(coef00m,levs*ntrac,MPI_R_MPI,ME_L_0,MC_COMP,
     &                                                       IERR)
          coef00=coef00m
          if( gen_coord_hybrid ) then                                    ! hmhj
            call updown_gc(sl,coef00(1,1))                               ! hmhj
          else                                                           ! hmhj
            call updown(sl,coef00(1,1))
          endif                                                          ! hmhj
          if (ntoz .gt. 1 .and. .not. (hybrid.or.gen_coord_hybrid)) then ! hmhj
            call updown(sl,coef00(1,ntoz))
          endif
          if (gg_tracers) then
!$omp parallel do shared(TRIE_LS,TRIO_LS)
!$omp+shared(deltim,SL,LS_NODE,coef00,hybrid)
            do k=1,levs
              CALL deldifs_tracers(
     .                TRIE_LS(1,1,P_RT+k-1), TRIE_LS(1,1,P_W+k-1),
     X                TRIE_LS(1,1,P_QM    ), TRIE_LS(1,1,P_X+k-1),
     X                TRIE_LS(1,1,P_Y +k-1), TRIE_LS(1,1,P_TEM+k-1),
     X                TRIO_LS(1,1,P_RT+k-1), TRIO_LS(1,1,P_W+k-1),
     X                TRIO_LS(1,1,P_QM    ), TRIO_LS(1,1,P_X+k-1),
     X                TRIO_LS(1,1,P_Y +k-1), TRIO_LS(1,1,P_TEM+k-1),
     X                deltim,SL,LS_NODE,coef00,k,hybrid,
     &                gen_coord_hybrid)
            enddo
          else
!
!$omp parallel do shared(TRIE_LS,TRIO_LS)
!$omp+shared(deltim,SL,LS_NODE,coef00,hybrid)
            do k=1,levs
              CALL deldifs(
     &                TRIE_LS(1,1,P_RT+k-1), TRIE_LS(1,1,P_W+k-1),
     &                TRIE_LS(1,1,P_QM    ), TRIE_LS(1,1,P_X+k-1),
     &                TRIE_LS(1,1,P_Y +k-1), TRIE_LS(1,1,P_TEM+k-1),
     &                TRIO_LS(1,1,P_RT+k-1), TRIO_LS(1,1,P_W+k-1),
     &                TRIO_LS(1,1,P_QM    ), TRIO_LS(1,1,P_X+k-1),
     &                TRIO_LS(1,1,P_Y +k-1), TRIO_LS(1,1,P_TEM+k-1),
     &                deltim,SL,LS_NODE,coef00,k,hybrid,
     &                gen_coord_hybrid)
            enddo
          endif
!--------------------------------------------------------
          do j=1,len_trie_ls
            trie_ls(j,1,p_q ) = trie_ls(j,1,p_zq)
            trie_ls(j,2,p_q ) = trie_ls(j,2,p_zq)
          enddo
          do j=1,len_trio_ls
            trio_ls(j,1,p_q ) = trio_ls(j,1,p_zq)
            trio_ls(j,2,p_q ) = trio_ls(j,2,p_zq)
          enddo
!         if (iprint .eq. 1) print*,' me = beg gloopb ',me
          timer1 = rtc()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! add impadj_slg to gloopb with batah, and set timetsteps to deltim
! add impadj_slg to gloopb with batah, and set timetsteps to deltim
! add impadj_slg to gloopb with batah, and set timetsteps to deltim
! add impadj_slg to gloopb with batah, and set timetsteps to deltim
! add impadj_slg to gloopb with batah, and set timetsteps to deltim
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          global_times_b = 0.0
          if(.not. adiab) then  ! third if.not.adiab
            call gloopb
     &        (trie_ls,trio_ls,
     &         ls_node,ls_nodes,max_ls_nodes,
     &         lats_nodes_a,global_lats_a,
     &         lats_nodes_r,global_lats_r,
     &         lonsperlar,
     &         epse,epso,epsedn,epsodn,
     &         snnp1ev,snnp1od,ndexev,ndexod,
     &         plnev_r,plnod_r,pddev_r,pddod_r,plnew_r,plnow_r,
     &         deltim,phour,sfc_fld, flx_fld, nst_fld, SFALB,
     &         xlon,
     &         swh,hlw,hprime,slag,sdec,cdec,
     &         ozplin,jindx1,jindx2,ddy,pdryini,
     &         phy_f3d,  phy_f2d, xlat,kdt,
     &         global_times_b,batah,lsout,fscav)
          endif            !  third if.not.adiab

!!$omp parallel do shared(trie_ls,ndexev,trio_ls,ndexod)
!!$omp+shared(sl,spdmax,deltim,ls_node)
!         do k=1,levs
!sela       call damp_speed(trie_ls(1,1,p_x+k-1), trie_ls(1,1,p_w +k-1),
!selax                  trie_ls(1,1,p_y+k-1), trie_ls(1,1,p_rt+k-1),
!selax                  ndexev,
!selax                  trio_ls(1,1,p_x+k-1), trio_ls(1,1,p_w +k-1),
!selax                  trio_ls(1,1,p_y+k-1), trio_ls(1,1,p_rt+k-1),
!selax                  ndexod,
!selax                  sl,spdmax(k),deltim,ls_node)
!         enddo

          do k=1,levs
            do j=1,len_trie_ls
              trie_ls(j,1,p_di+k-1) = trie_ls(j,1,p_x+k-1)
              trie_ls(j,2,p_di+k-1) = trie_ls(j,2,p_x+k-1)
              trie_ls(j,1,p_ze+k-1) = trie_ls(j,1,p_w+k-1)
              trie_ls(j,2,p_ze+k-1) = trie_ls(j,2,p_w+k-1)
              trie_ls(j,1,p_te+k-1) = trie_ls(j,1,p_y+k-1)
              trie_ls(j,2,p_te+k-1) = trie_ls(j,2,p_y+k-1)
            enddo
          enddo
          do k=1,levs
            do j=1,len_trio_ls
              trio_ls(j,1,p_di+k-1) = trio_ls(j,1,p_x+k-1)
              trio_ls(j,2,p_di+k-1) = trio_ls(j,2,p_x+k-1)
              trio_ls(j,1,p_ze+k-1) = trio_ls(j,1,p_w+k-1)
              trio_ls(j,2,p_ze+k-1) = trio_ls(j,2,p_w+k-1)
              trio_ls(j,1,p_te+k-1) = trio_ls(j,1,p_y+k-1)
              trio_ls(j,2,p_te+k-1) = trio_ls(j,2,p_y+k-1)
            enddo
          enddo
          if(.not. gg_tracers)then
            do k=1,levh
              do j=1,len_trie_ls
                trie_ls(j,1,p_rq+k-1) = trie_ls(j,1,p_rt+k-1)
                trie_ls(j,2,p_rq+k-1) = trie_ls(j,2,p_rt+k-1)
              enddo
            enddo
            do k=1,levh
              do j=1,len_trio_ls
                trio_ls(j,1,p_rq+k-1) = trio_ls(j,1,p_rt+k-1)
                trio_ls(j,2,p_rq+k-1) = trio_ls(j,2,p_rt+k-1)
              enddo
            enddo
          endif !  if(.not.gg_tracers)
!
!----------------------------------------------------------
        else                          ! Eulerian Dynamics
!----------------------------------------------------------
!!
          if(ifirst == 1) then
            k = 0
            CALL deldifs(
     &                TRIE_LS(1,1,P_RT+k-1), TRIE_LS(1,1,P_W+k-1),
     &                TRIE_LS(1,1,P_QM    ), TRIE_LS(1,1,P_X+k-1),
     &                TRIE_LS(1,1,P_Y +k-1), TRIE_LS(1,1,P_TEM+k-1),    ! hmhj
     &                TRIO_LS(1,1,P_RT+k-1), TRIO_LS(1,1,P_W+k-1),
     &                TRIO_LS(1,1,P_QM    ), TRIO_LS(1,1,P_X+k-1),
     &                TRIO_LS(1,1,P_Y +k-1), TRIO_LS(1,1,P_TEM+k-1),    ! hmhj
     &                deltim,SL,LS_NODE,coef00,k,hybrid,                ! hmhj
     &                gen_coord_hybrid)

            ifirst=0
          endif
          global_times_a=0.

!     print *,' Eulerian dynamics callling GLOOPA for kdt=',kdt

          CALL GLOOPA
     &      (deltim,TRIE_LS,TRIO_LS,
     &       LS_NODE,LS_NODES,MAX_LS_NODES,
     &       LATS_NODES_A,GLOBAL_LATS_A,
     &       LONSPERLAT,
     &       EPSE,EPSO,EPSEDN,EPSODN,
     &       SNNP1EV,SNNP1OD,NDEXEV,NDEXOD,
     &       PLNEV_A,PLNOD_A,PDDEV_A,PDDOD_A,PLNEW_A,PLNOW_A,
     &       global_times_a,kdt)
       
!
!
          iprint = 0
!         if (iprint .eq. 1) print*,' fin gloopa kdt = ',kdt
!
!my gather lat timings for load balancing
!sela     if (reshuff_lats_a .and. kdt .eq. 5) then
!sela       call redist_lats_a(kdt,global_times_a,
!selax                lats_nodes_a,global_lats_a,
!selax                lonsperlat,
!selax                lats_nodes_ext,global_lats_ext,iprint)
!sela     endif
!----------------------------------------------------------

!
          if(.not. adiab) then
            if (nscyc > 0 .and. mod(kdt,nscyc) == 1) then
             CALL gcycle(me,LATS_NODE_R,LONSPERLAR,global_lats_r,
     &                  ipt_lats_node_r,idate,fhour,fhcyc,
     &                  XLON ,XLAT  , sfc_fld, ialb)
            endif
!
            if (num_p3d  == 3) then  ! Ferrier Microphysics initialization
              call init_micro(deltim,lonr,levs,num_p3d,lats_node_r,
     &                        phy_f3d(1,1,1,1),   fhour, me)
            endif
          endif
!
!-> Coupling insertion

  ! lgetSSTICE_cc must be defined by this moment. It used to be an argument
  ! to ATM_GETSST, accessible here via USE SURFACE_cc. Now it is defined in
  ! ATM_TSTEP_INIT called above, and the USE is removed. (Even in the earlier
  ! version lgetSSTICE_cc did not have to be an actual argumnent, since
  ! it is in the module SURFACE_cc USEd by ATM_GETSST.)

          if (coupler_id >= 0) then
            do j = 1, lats_node_r
              do i = 1, lonr
                if (sfc_fld%slmsk(i,j) == 0 ) then
                  sfc_fld%tsea(i,j) = sfc_fld%tsea(i,j)
     &                 + (sfc_fld%oro(i,j)-sfc_fld%oro_uf(i,j))*rlapse
                endif
              enddo
            enddo
          endif

!     print *,' AM: before ATM_GETSSTICE coupler_id=',coupler_id,
!    &' TSEA=',sfc_fld%TSEA(1,1),' me=',me
          call ATM_GETSSTICE(sfc_fld%TSEA,sfc_fld%TISFC,sfc_fld%FICE,
     &                       sfc_fld%HICE,sfc_fld%SHELEG,sfc_fld%SLMSK,
     &                       kdt)
!     print *,' AM: after ATM_GETSSTICE coupler_id=',coupler_id,
!    & ' me=',me

          if (coupler_id >= 0) then
            do j = 1, lats_node_r
              do i = 1, lonr
                if (sfc_fld%slmsk(i,j) == 0 ) then
                  sfc_fld%tsea(i,j) = sfc_fld%tsea(i,j)
     &                 - (sfc_fld%oro(i,j)-sfc_fld%oro_uf(i,j))*rlapse
                endif
              enddo
            enddo
          endif

!<- Coupling insertion

!
          if (nst_fcst > 1) then                         ! update TSEA
            if (Coupler_id < 0 .or. .not. mom4ice) then  ! Standalone mode
              do j = 1, lats_node_r
                do i = 1, lonr
                  if (sfc_fld%slmsk(i,j) == 0 ) then
                    dt_warm = (nst_fld%xt(i,j)+nst_fld%xt(i,j))
     &                      /  nst_fld%xz(i,j)
                    sfc_fld%TSEA(i,j) = nst_fld%tref(i,j)
     &                  + dt_warm - nst_fld%dt_cool(i,j)
     &                  - (sfc_fld%oro(i,j)-sfc_fld%oro_uf(i,j))*rlapse
                  endif
                enddo
              enddo
            else                                         ! Coupled to MOM4 OM
              tem1 = 0.5 / omz1
              do j = 1, lats_node_r
                do i = 1, lonr
                  if (sfc_fld%slmsk(i,j) == 0 ) then
                    tem2 = 1.0 / nst_fld%xz(i,j)
                    sfc_fld%tsea(i,j) = sfc_fld%tsea(i,j)
     &                   + (sfc_fld%oro(i,j)-sfc_fld%oro_uf(i,j))*rlapse
                    dt_warm = (nst_fld%xt(i,j)+nst_fld%xt(i,j)) * tem2

                    if ( nst_fld%xz(i,j) > omz1) then
                      nst_fld%tref(i,j) = sfc_fld%tsea(i,j)
     &                 - (1.0-0.5*omz1*tem2) * dt_warm
     &                 + nst_fld%z_c(i,j)*nst_fld%dt_cool(i,j)*tem1
                    else
                     nst_fld%tref(i,j) = sfc_fld%tsea(i,j)
     &                 - (nst_fld%xz(i,j)*dt_warm
     &                 -  nst_fld%z_c(i,j)*nst_fld%dt_cool(i,j))*tem1
                    endif
                    sfc_fld%TSEA(i,j) = nst_fld%tref(i,j)
     &                  + dt_warm - nst_fld%dt_cool(i,j)
     &                  - (sfc_fld%oro(i,j)-sfc_fld%oro_uf(i,j))*rlapse
                  endif
                enddo
              enddo
            endif
          endif

!           do j = 1, lats_node_r
!             do i = 1, lonr
!               if (sfc_fld%slmsk(i,j) == 0 ) then
!                 dt_warm = (nst_fld%xt(i,j)+nst_fld%xt(i,j))
!    &                    /  nst_fld%xz(i,j)
!                 sfc_fld%TSEA(i,j) = sfc_fld%TSEA(i,j)
!    &                    + dt_warm - nst_fld%dt_cool(i,j)
!               endif
!             enddo
!           enddo
!         endif


!sela     if (me.eq.0) PRINT*,'COMPLETED GLOOPA IN do_tstep'

          global_times_r = 0.0               !my set to zero for every timestep

          if (lsswr .or. lslwr) then         ! Radiation Call!
            if(.not. adiab) then

!           if(.not.adiab .and. kdt > 1) then
!     print *,' before calling GLOOPR kdt=',kdt

              call gloopr
     &          (trie_ls,trio_ls,
     &           ls_node,ls_nodes,max_ls_nodes,
     &           lats_nodes_a,global_lats_a,
     &           lats_nodes_r,global_lats_r,
     &           lonsperlar,
     &           epse,epso,epsedn,epsodn,
     &           snnp1ev,snnp1od,plnev_r,plnod_r,
     &           pddev_r,pddod_r,
     &           phour,
     &           xlon,xlat,coszdg,flx_fld%coszen,
     &           sfc_fld%slmsk,sfc_fld%sheleg,sfc_fld%SNCOVR,
     &           sfc_fld%SNOALB,sfc_fld%ZORL,sfc_fld%TSEA,
     &           HPRIME,SFALB,sfc_fld%ALVSF,sfc_fld%ALNSF,
     &           sfc_fld%ALVWF,sfc_fld%ALNWF,sfc_fld%FACSF,
     &           sfc_fld%FACWF,sfc_fld%CV,sfc_fld%CVT ,
     &           sfc_fld%CVB,SWH,HLW,flx_fld%SFCNSW,flx_fld%SFCDLW,
     &           sfc_fld%FICE,sfc_fld%TISFC,flx_fld%SFCDSW,
     &           flx_fld%sfcemis,                                    ! yth 4/09
     &           flx_fld%TSFLW,FLUXR,phy_f3d,SLAG,SDEC,CDEC,KDT,
     &           global_times_r)
            endif               ! second  if.not.adiab
          endif  !sswr .or. lslwr

!     if (me == 0) then
!     print *,' aft gloopr HLW45=',hlw(1,:,45)
!     print *,' aft gloopr SWH45=',swh(1,:,45)
!     endif
!!
!         print *,' finished GLOOPR at kdt=',kdt
!         call mpi_quit(1111)

          if( .not. explicit ) then					! hmhj
!
            if( gen_coord_hybrid ) then                                 ! hmhj

!$omp parallel do private(locl)
              do locl=1,ls_max_node                                     ! hmhj
                call sicdife_hyb_gc(
     &                       trie_ls(1,1,P_dim), trie_ls(1,1,P_tem),    ! hmhj
     &                       trie_ls(1,1,P_qm ), trie_ls(1,1,P_x  ),    ! hmhj
     &                       trie_ls(1,1,P_y  ), trie_ls(1,1,P_zq ),    ! hmhj
     &                       trie_ls(1,1,P_di ), trie_ls(1,1,P_te ),    ! hmhj
     &                       trie_ls(1,1,P_q  ),deltim,                 ! hmhj
     &                       trie_ls(1,1,P_uln), trie_ls(1,1,P_vln),    ! hmhj
     &                       ls_node,snnp1ev,ndexev,locl)               ! hmhj

                call sicdifo_hyb_gc(
     &                       trio_ls(1,1,P_dim), trio_ls(1,1,P_tem),    ! hmhj
     &                       trio_ls(1,1,P_qm ), trio_ls(1,1,P_x  ),    ! hmhj
     &                       trio_ls(1,1,P_y  ), trio_ls(1,1,P_zq ),    ! hmhj
     &                       trio_ls(1,1,P_di ), trio_ls(1,1,P_te ),    ! hmhj
     &                       trio_ls(1,1,P_q  ),deltim,                 ! hmhj
     &                       trio_ls(1,1,P_uln), trio_ls(1,1,P_vln),    ! hmhj
     &                       ls_node,snnp1od,ndexod,locl)               ! hmhj
              enddo                                                     ! hmhj

            else if(hybrid)then                                         ! hmhj

!         print *,' calling sicdife_hyb at kdt=',kdt
!$omp parallel do private(locl)
              do locl=1,ls_max_node
                call sicdife_hyb(
     &                    trie_ls(1,1,P_dim), trie_ls(1,1,P_tem),
     &                    trie_ls(1,1,P_qm ), trie_ls(1,1,P_x  ),
     &                    trie_ls(1,1,P_y  ), trie_ls(1,1,P_zq ),
     &                    trie_ls(1,1,P_di ), trie_ls(1,1,P_te ),
     &                    trie_ls(1,1,P_q  ),deltim,
     &                    trie_ls(1,1,P_uln), trie_ls(1,1,P_vln),
     &                    ls_node,snnp1ev,ndexev,locl)

                 call sicdifo_hyb(
     &                    trio_ls(1,1,P_dim), trio_ls(1,1,P_tem),
     &                    trio_ls(1,1,P_qm ), trio_ls(1,1,P_x  ),
     &                    trio_ls(1,1,P_y  ), trio_ls(1,1,P_zq ),
     &                    trio_ls(1,1,P_di ), trio_ls(1,1,P_te ),
     &                    trio_ls(1,1,P_q  ),deltim,
     &                    trio_ls(1,1,P_uln), trio_ls(1,1,P_vln),
     &                    ls_node,snnp1od,ndexod,locl)
              enddo

!         print *,' after calling sicdife_hyb at kdt=',kdt
            else ! hybrid

!$omp parallel do private(locl)
              do locl=1,ls_max_node
                CALL SICDIFE_sig(
     &                    TRIE_LS(1,1,P_DIM), TRIE_LS(1,1,P_TEM),
     &                    TRIE_LS(1,1,P_QM ), TRIE_LS(1,1,P_X  ),
     &                    TRIE_LS(1,1,P_Y  ), TRIE_LS(1,1,P_ZQ ),
     &                    AM,BM,TOV,SV,deltim,
     &                    TRIE_LS(1,1,P_ULN), TRIE_LS(1,1,P_VLN),
     &                    LS_NODE,SNNP1EV,NDEXEV,locl,TRIE_LS(1,1,P_DI))

                CALL SICDIFO_sig(
     &                    TRIO_LS(1,1,P_DIM), TRIO_LS(1,1,P_TEM),
     &                    TRIO_LS(1,1,P_QM ), TRIO_LS(1,1,P_X  ),  
     &                    TRIO_LS(1,1,P_Y  ), TRIO_LS(1,1,P_ZQ ),
     &                    AM,BM,TOV,SV,deltim,
     &                    TRIO_LS(1,1,P_ULN), TRIO_LS(1,1,P_VLN),
     &                    LS_NODE,SNNP1OD,NDEXOD,locl,TRIO_LS(1,1,P_DI))
              enddo
            endif ! hybrid

          endif 		! not explicit	 			! hmhj

!
!----------------------------------------------------------
!sela     if (.NOT.LIOPE.or.icolor.ne.2) then
!sela       print*,'liope=',liope,' icolor=',icolor,' after loopa'
!sela       CALL RMS_spect(TRIE_LS(1,1,P_zq ), TRIE_LS(1,1,P_x  ),
!selaX             TRIE_LS(1,1,P_y  ), TRIE_LS(1,1,P_w  ),
!selaX             TRIE_LS(1,1,P_Rt ),
!selaX             TRIO_LS(1,1,P_zq ), TRIO_LS(1,1,P_x  ),
!selaX             TRIO_LS(1,1,P_y  ), TRIO_LS(1,1,P_w  ),
!selaX             TRIO_LS(1,1,P_Rt ),
!selaX             LS_NODES,MAX_LS_NODES)
!sela     endif
!----------------------------------------------------------

! hmhj compute coef00 for all, even for hybrid mode

          coef00(:,:) = 0.0
          IF ( ME .EQ. ME_L_0 ) THEN
            DO LOCL=1,LS_MAX_NODE
              l      = ls_node(locl,1)
              jbasev = ls_node(locl,2)
              IF ( L  ==  0 ) THEN
                N = 0
! 1 Corresponds to temperature,  2 corresponds to ozone, 3 to cloud condensate
                DO K=1,LEVS
                  coef00(K,1) = TRIE_LS(INDLSEV(N,L),1,P_Y +K-1)
!                 if (ntoz .gt. 1) then
                  if (ntoz .gt. 1 .and.                                   ! hmhj
     &               .not. (hybrid.or.gen_coord_hybrid)) then             ! hmhj
                     coef00(K,ntoz) = TRIE_LS(INDLSEV(N,L),1,
     &                                   (ntoz-1)*levs+P_rt+K-1)
                  endif
                ENDDO
              ENDIF
            END DO
          END IF
          coef00m = coef00
          CALL MPI_BCAST(coef00m,levs*ntrac,MPI_R_MPI,ME_L_0,MC_COMP,
     &                                                       IERR)
          coef00=coef00m
          if( gen_coord_hybrid ) then                                     ! hmhj
            call updown_gc(sl,coef00(1,1))                                ! hmhj
          else                                                            ! hmhj
            call updown(sl,coef00(1,1))
          endif                                                           ! hmhj
!         if (ntoz > 1) call updown(sl,coef00(1,ntoz))
          if (ntoz > 1 .and. .not. (hybrid.or.gen_coord_hybrid)) then     ! hmhj
               call updown(sl,coef00(1,ntoz))
          endif

!         print *,' calling deldifs at kdt=',kdt
!
!$omp parallel do shared(TRIE_LS,TRIO_LS)
!$omp+shared(deltim,SL,LS_NODE,coef00,hybrid,gen_coord_hybrid)
          do k=1,levs
            CALL deldifs(TRIE_LS(1,1,P_RT+k-1), TRIE_LS(1,1,P_W+k-1),
     X                   TRIE_LS(1,1,P_QM    ), TRIE_LS(1,1,P_X+k-1),
     X                   TRIE_LS(1,1,P_Y +k-1), TRIE_LS(1,1,P_TEM+k-1),   ! hmhj
     X                   TRIO_LS(1,1,P_RT+k-1), TRIO_LS(1,1,P_W+k-1),
     X                   TRIO_LS(1,1,P_QM    ), TRIO_LS(1,1,P_X+k-1),
     X                   TRIO_LS(1,1,P_Y +k-1), TRIO_LS(1,1,P_TEM+k-1),   ! hmhj
     X                   deltim,SL,LS_NODE,coef00,k,hybrid,               ! hmhj
     &                   gen_coord_hybrid)                                ! hmhj
          enddo
!         print *,' after calling deldifs at kdt=',kdt
!
!
!-------------------------------------------
          if(.not.lsfwd)then
!-------------------------------------------
            CALL FILTR1EO(TRIE_LS(1,1,P_TEM), TRIE_LS(1,1,P_TE ),
     &                    TRIE_LS(1,1,P_Y  ), TRIE_LS(1,1,P_DIM),
     &                    TRIE_LS(1,1,P_DI ), TRIE_LS(1,1,P_X  ),
     &                    TRIE_LS(1,1,P_ZEM), TRIE_LS(1,1,P_ZE ),
     &                    TRIE_LS(1,1,P_W  ), TRIE_LS(1,1,P_RM ),
     &                    TRIE_LS(1,1,P_RQ ), TRIE_LS(1,1,P_RT ),
     &                    TRIO_LS(1,1,P_TEM), TRIO_LS(1,1,P_TE ),
     &                    TRIO_LS(1,1,P_Y  ), TRIO_LS(1,1,P_DIM),
     &                    TRIO_LS(1,1,P_DI ), TRIO_LS(1,1,P_X  ),
     &                    TRIO_LS(1,1,P_ZEM), TRIO_LS(1,1,P_ZE ),
     &                    TRIO_LS(1,1,P_W  ), TRIO_LS(1,1,P_RM ),
     &                    TRIO_LS(1,1,P_RQ ), TRIO_LS(1,1,P_RT ),
     &                    FILTA,LS_NODE)

            CALL countperf(0,13,0.)
            DO J=1,LEN_TRIE_LS
              TRIE_LS(J,1,P_QM) = TRIE_LS(J,1,P_Q )
              TRIE_LS(J,2,P_QM) = TRIE_LS(J,2,P_Q )
              TRIE_LS(J,1,P_Q ) = TRIE_LS(J,1,P_ZQ)
              TRIE_LS(J,2,P_Q ) = TRIE_LS(J,2,P_ZQ)
            ENDDO
            DO J=1,LEN_TRIO_LS
              TRIO_LS(J,1,P_QM) = TRIO_LS(J,1,P_Q )
              TRIO_LS(J,2,P_QM) = TRIO_LS(J,2,P_Q )
              TRIO_LS(J,1,P_Q ) = TRIO_LS(J,1,P_ZQ)
              TRIO_LS(J,2,P_Q ) = TRIO_LS(J,2,P_ZQ)
            ENDDO
            CALL countperf(1,13,0.)

!-------------------------------------------
          else
!-------------------------------------------
            CALL countperf(0,13,0.)
            DO J=1,LEN_TRIE_LS
              TRIE_LS(J,1,P_Q) = TRIE_LS(J,1,P_ZQ)
              TRIE_LS(J,2,P_Q) = TRIE_LS(J,2,P_ZQ)
            ENDDO
            DO J=1,LEN_TRIO_LS
              TRIO_LS(J,1,P_Q) = TRIO_LS(J,1,P_ZQ)
              TRIO_LS(J,2,P_Q) = TRIO_LS(J,2,P_ZQ)
            ENDDO
            CALL countperf(1,13,0.)
!-------------------------------------------
          endif
!
!-------------------------------------------
!         if (iprint .eq. 1) print*,' me = beg gloopb ',me
!my set to zero for every timestep
          global_times_b = 0.0

          if(.not. adiab) then

!        print *,' before calling GLOOPB kdt=',kdt
            call gloopb
     &        (trie_ls,trio_ls,
     &         ls_node,ls_nodes,max_ls_nodes,
     &         lats_nodes_a,global_lats_a,
     &         lats_nodes_r,global_lats_r,
     &         lonsperlar,
     &         epse,epso,epsedn,epsodn,
     &         snnp1ev,snnp1od,ndexev,ndexod,
     &         plnev_r,plnod_r,pddev_r,pddod_r,plnew_r,plnow_r,
     &         deltim,phour,sfc_fld, flx_fld, nst_fld, SFALB,
     &         xlon,
     &         swh,hlw,hprime,slag,sdec,cdec,
     &         ozplin,jindx1,jindx2,ddy,pdryini,
     &         phy_f3d,  phy_f2d, xlat,kdt,
     &         global_times_b,batah,lsout,fscav)
!
!           if (kdt .eq. 1) call mpi_quit(222)
          endif ! not.adiab

!        print *,' after calling GLOOPB kdt=',kdt
!
!!$omp parallel do shared(TRIE_LS,NDEXEV,TRIO_LS,NDEXOD)
!!$omp+shared(SL,SPDMAX,deltim,LS_NODE)
!$omp parallel do private(k)
          do k=1,levs
            CALL damp_speed(TRIE_LS(1,1,P_X+k-1), TRIE_LS(1,1,P_W +k-1),
     &                      TRIE_LS(1,1,P_Y+k-1), TRIE_LS(1,1,P_RT+k-1),
     &                      NDEXEV,
     &                      TRIO_LS(1,1,P_X+k-1), TRIO_LS(1,1,P_W +k-1),
     &                      TRIO_LS(1,1,P_Y+k-1), TRIO_LS(1,1,P_RT+k-1),
     &                      NDEXOD,
     &                      SL,SPDMAX(k),deltim,LS_NODE)
          enddo
!
!--------------------------------------------
          if(.not. lsfwd)then
!--------------------------------------------
            CALL FILTR2EO(TRIE_LS(1,1,P_TEM), TRIE_LS(1,1,P_TE ),
     &                    TRIE_LS(1,1,P_Y  ), TRIE_LS(1,1,P_DIM),
     &                    TRIE_LS(1,1,P_DI ), TRIE_LS(1,1,P_X  ),
     &                    TRIE_LS(1,1,P_ZEM), TRIE_LS(1,1,P_ZE ),
     &                    TRIE_LS(1,1,P_W  ), TRIE_LS(1,1,P_RM ),
     &                    TRIE_LS(1,1,P_RQ ), TRIE_LS(1,1,P_RT ),
     &                    TRIO_LS(1,1,P_TEM), TRIO_LS(1,1,P_TE ),
     &                    TRIO_LS(1,1,P_Y  ), TRIO_LS(1,1,P_DIM),
     &                    TRIO_LS(1,1,P_DI ), TRIO_LS(1,1,P_X  ),
     &                    TRIO_LS(1,1,P_ZEM), TRIO_LS(1,1,P_ZE ),
     &                    TRIO_LS(1,1,P_W  ), TRIO_LS(1,1,P_RM ),
     &                    TRIO_LS(1,1,P_RQ ), TRIO_LS(1,1,P_RT ),
     &                    FILTA,LS_NODE)
!--------------------------------------------
          else
!--------------------------------------------
            CALL countperf(0,13,0.)
            DO K=1,LEVS
              DO J=1,LEN_TRIE_LS
                TRIE_LS(J,1,P_DI+K-1) = TRIE_LS(J,1,P_X+K-1)
                TRIE_LS(J,2,P_DI+K-1) = TRIE_LS(J,2,P_X+K-1)
                TRIE_LS(J,1,P_ZE+K-1) = TRIE_LS(J,1,P_W+K-1)
                TRIE_LS(J,2,P_ZE+K-1) = TRIE_LS(J,2,P_W+K-1)
                TRIE_LS(J,1,P_TE+K-1) = TRIE_LS(J,1,P_Y+K-1)
                TRIE_LS(J,2,P_TE+K-1) = TRIE_LS(J,2,P_Y+K-1)
              ENDDO
            ENDDO
            DO K=1,LEVS
              DO J=1,LEN_TRIO_LS
                TRIO_LS(J,1,P_DI+K-1) = TRIO_LS(J,1,P_X+K-1)
                TRIO_LS(J,2,P_DI+K-1) = TRIO_LS(J,2,P_X+K-1)
                TRIO_LS(J,1,P_ZE+K-1) = TRIO_LS(J,1,P_W+K-1)
                TRIO_LS(J,2,P_ZE+K-1) = TRIO_LS(J,2,P_W+K-1)
                TRIO_LS(J,1,P_TE+K-1) = TRIO_LS(J,1,P_Y+K-1)
                TRIO_LS(J,2,P_TE+K-1) = TRIO_LS(J,2,P_Y+K-1)
              ENDDO
            ENDDO
            DO K=1,LEVH
              DO J=1,LEN_TRIE_LS
                TRIE_LS(J,1,P_RQ+K-1) = TRIE_LS(J,1,P_RT+K-1)
                TRIE_LS(J,2,P_RQ+K-1) = TRIE_LS(J,2,P_RT+K-1)
              ENDDO
            ENDDO
            DO K=1,LEVH
              DO J=1,LEN_TRIO_LS
                TRIO_LS(J,1,P_RQ+K-1) = TRIO_LS(J,1,P_RT+K-1)
                TRIO_LS(J,2,P_RQ+K-1) = TRIO_LS(J,2,P_RT+K-1)
              ENDDO
            ENDDO
            CALL countperf(1,13,0.)
!--------------------------------------------
          endif
!         if (kdt .eq. 2) call mpi_quit(444)
!!
        endif                   ! if (semilag) then loop

      endif                     ! if(comp_task) then loop
!     endif                     !.NOT.LIOPE.or.icolor.ne.2
!
!--------------------------------------------
!--------------------------------------------
      IF (lsout) THEN
!!
        CALL f_hpmstart(32,"TWRITEEO")
!!
        CALL countperf(0,18,0.)
!
        wrt_g3d = MOD(kdt ,nsout) == 0 .or. phour == 0.0
        CALL WRTOUT(PHOUR,FHOUR,ZHOUR,IDATE,
     &              TRIE_LS,TRIO_LS,
     &              SL,SI,
     &              ls_node,LS_NODES,MAX_LS_NODES,
     &              sfc_fld, flx_fld, nst_fld,
     &              fluxr,pdryini,
     &              lats_nodes_r,global_lats_r,lonsperlar,
     &              COLAT1,CFHOUR1,pl_coeff,
     &              epsedn,epsodn,snnp1ev,snnp1od,plnev_r,plnod_r,
     &              plnew_r,plnow_r,'SIG.F','SFC.F','FLX.F',wrt_g3d)


!     endif
!
      CALL f_hpmstop(32)
!!
      CALL countperf(1,18,0.)
!!
!!
      IF (mod(kdt,nsres) == 0 .and. (.not. SPS)) THEN
!!
        CALL wrt_restart(TRIE_LS,TRIO_LS,
     &       sfc_fld, nst_fld,
     &       SI,SL,fhour,idate,
     &       igen,pdryini,
     x       ls_node,ls_nodes,max_ls_nodes,
     &       global_lats_r,lonsperlar,SNNP1EV,SNNP1OD,
     &       phy_f3d, phy_f2d, ngptc, adiab, ens_nam,
     &       nst_fcst,'SIGR1','SIGR2','SFCR','NSTR')
!
        ENDIF
      ENDIF ! if ls_out
!
      IF (mod(kdt,nszer) == 0 .and. lsout) THEN
        call flx_init(flx_fld,ierr)
        zhour = fhour
        FLUXR = 0.
!
        if (ldiag3d .or. lggfs3d) then
          call d3d_zero(ldiag3d,lggfs3d)
          if (fhour >= fhgoc3d) lggfs3d = .false.
        endif
      ENDIF
!
! Coupling insertion->
      if (comp_task) then
        CALL ATM_SENDFLUXES(sfc_fld%SLMSK)
      endif
!<- Coupling insertion

      RETURN
      END
