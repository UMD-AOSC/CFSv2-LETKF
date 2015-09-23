      SUBROUTINE tldfi(deltim,kdt,PHOUR,
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
!    &                 LSLAG,
     &                 XLON,XLAT,COSZDG,sfc_fld,flx_fld,nst_fld,
     &                 HPRIME,SWH,HLW,FLUXR,
     &                 SFALB,SLAG,SDEC,CDEC,
     &                 OZPLIN,JINDX1,JINDX2,
     &                 DDY,PDRYINI,
     &                 phy_f3d,  phy_f2d,
     &                 ZHOUR,N1,N4,LSOUT,COLAT1,
     &                 CFHOUR1,fscav)
       
      use machine             , only : kind_evod,kind_phys,kind_rad
      use resol_def           , only : latg,latg2,latr,latr2,levh,levs,
     &                                 lnte,lnto,lonf,lonr,lotd,lots,
     &                                 lsoil,nfxr,nmtvr,ntrac,ncld,
     &                                 num_p2d,num_p3d,
     &                                 p_di,p_dim,p_q,p_qm,p_rm,p_rq,
     &                                 p_te,p_tem,p_ze,p_zem
      use layout1             , only : lats_dim_a,lats_node_r,
     &                                 len_trie_ls,len_trio_ls,
     &                                 ls_dim,me,nodes
      use layout_grid_tracers , only : rgt_a
      use vert_def            , only : am,bm,sv,tov
      use date_def            , only : fhour
      use namelist_def        , only : fhdfi,gen_coord_hybrid,
     &                                 gg_tracers,hybrid,
     &                                 lscca,lsfwd,lslwr,lssav,lsswr,
     &                                 nsdfi,nslwr,nsout,nsswr,
     &                                 semilag,nsout_hf,fhmax_hf
      use mpi_def             , only : comp_task
      use ozne_def            , only : latsozp,levozp,pl_coeff,timeoz

!     use resol_def
!     use layout1
      use gg_def
!     use vert_def
!     use date_def
!     use namelist_def
!     use ozne_def
      use Sfc_Flx_ESMFMod
      use Nst_Var_ESMFMod
!
      IMPLICIT NONE
!     include 'mpif.h'      

!!     
      TYPE(Sfc_Var_Data)        :: sfc_fld
      TYPE(Flx_Var_Data)        :: flx_fld
      TYPE(Nst_Var_Data)        :: nst_fld

      CHARACTER(16)                     :: CFHOUR1
      INTEGER,INTENT(IN):: LONSPERLAT(LATG),N1,N4
!!     
      REAL(KIND=KIND_EVOD),INTENT(INOUT):: deltim,PHOUR,ZHOUR
!!     
      INTEGER I
      REAL(KIND=KIND_EVOD) TRIE_LS(LEN_TRIE_LS,2,11*LEVS+3*LEVH+6)
      REAL(KIND=KIND_EVOD) TRIO_LS(LEN_TRIO_LS,2,11*LEVS+3*LEVH+6)
      INTEGER              LS_NODE (LS_DIM)
      INTEGER              LS_NODES(LS_DIM,NODES)
      INTEGER          MAX_LS_NODES(NODES)
      INTEGER               LATS_NODES_A(NODES)
!     INTEGER               LATS_NODES_EXT(NODES)
      INTEGER              GLOBAL_LATS_A(LATG)
!     INTEGER        GLOBAL_LATS_EXT(LATG+2*JINTMX+2*NYPT*(NODES-1))
      INTEGER               LATS_NODES_R(NODES)
      INTEGER              GLOBAL_LATS_R(LATR)
      INTEGER                 LONSPERLAR(LATR)
!!
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

c$$$      INTEGER                LOTS,LOTD,LOTA
c$$$      PARAMETER            ( LOTS = 5*LEVS+1*LEVH+3 )
c$$$      PARAMETER            ( LOTD = 6*LEVS+2*LEVH+0 )
c$$$      PARAMETER            ( LOTA = 3*LEVS+1*LEVH+1 )

      REAL(KIND=KIND_EVOD) SYN_LS_A(4*LS_DIM,LOTS,LATG2)
      REAL(KIND=KIND_EVOD) DYN_LS_A(4*LS_DIM,LOTD,LATG2)

!     REAL(KIND=KIND_EVOD) SYN_GR_A_1(LONFX*LOTS,LATS_DIM_EXT)
!     REAL(KIND=KIND_EVOD) DYN_GR_A_1(LONFX*LOTD,LATS_DIM_EXT)
!     REAL(KIND=KIND_EVOD) ANL_GR_A_1(LONFX*LOTA,LATS_DIM_EXT)
!     REAL(KIND=KIND_EVOD) SYN_GR_A_2(LONFX*LOTS,LATS_DIM_EXT)
!     REAL(KIND=KIND_EVOD) DYN_GR_A_2(LONFX*LOTD,LATS_DIM_EXT)
!     REAL(KIND=KIND_EVOD) ANL_GR_A_2(LONFX*LOTA,LATS_DIM_EXT)
!!     
      REAL (KIND=KIND_RAD) XLON(LONR,LATS_NODE_R)
      REAL (KIND=KIND_RAD) XLAT(LONR,LATS_NODE_R)
      REAL (KIND=KIND_RAD) COSZDG(LONR,LATS_NODE_R),
     &                     HPRIME(LONR,NMTVR,LATS_NODE_R),
     &                     FLUXR(LONR,nfxr,LATS_NODE_R),
     &                     SFALB(LONR,LATS_NODE_R),
     &                     SWH(LONR,LEVS,LATS_NODE_R),
     &                     HLW(LONR,LEVS,LATS_NODE_R)

      REAL (kind=kind_phys)
     &     phy_f3d(lonr,LEVS,num_p3d,lats_node_r),
     &     phy_f2d(lonr,num_p2d,lats_node_r),
!
     &     DDY(LATS_NODE_R), fscav(ntrac-ncld-1)

      INTEGER JINDX1(LATS_NODE_R),JINDX2(LATS_NODE_R)
!!     
      INTEGER LEV,LEVMAX
      REAL OZPLIN(LATSOZP,LEVOZP,pl_coeff,timeoz) !OZONE Coeff
      REAL (KIND=KIND_PHYS) PDRYINI
      REAL(KIND=KIND_EVOD) SLAG,SDEC,CDEC

c$$$      INTEGER   P_GZ,P_ZEM,P_DIM,P_TEM,P_RM,P_QM
c$$$      INTEGER   P_ZE,P_DI,P_TE,P_RQ,P_Q,P_DLAM,P_DPHI,P_ULN,P_VLN
c$$$      INTEGER   P_W,P_X,P_Y,P_RT,P_ZQ
c$$$      PARAMETER(P_GZ  = 0*LEVS+0*LEVH+1,  !      GZE/O(LNTE/OD,2),
c$$$     X          P_ZEM = 0*LEVS+0*LEVH+2,  !     ZEME/O(LNTE/OD,2,LEVS),
c$$$     X          P_DIM = 1*LEVS+0*LEVH+2,  !     DIME/O(LNTE/OD,2,LEVS),
c$$$     X          P_TEM = 2*LEVS+0*LEVH+2,  !     TEME/O(LNTE/OD,2,LEVS),
c$$$     X          P_RM  = 3*LEVS+0*LEVH+2,  !      RME/O(LNTE/OD,2,LEVH),
c$$$     X          P_QM  = 3*LEVS+1*LEVH+2,  !      QME/O(LNTE/OD,2),
c$$$     X          P_ZE  = 3*LEVS+1*LEVH+3,  !      ZEE/O(LNTE/OD,2,LEVS),
c$$$     X          P_DI  = 4*LEVS+1*LEVH+3,  !      DIE/O(LNTE/OD,2,LEVS),
c$$$     X          P_TE  = 5*LEVS+1*LEVH+3,  !      TEE/O(LNTE/OD,2,LEVS),
c$$$     X          P_RQ  = 6*LEVS+1*LEVH+3,  !      RQE/O(LNTE/OD,2,LEVH),
c$$$     X          P_Q   = 6*LEVS+2*LEVH+3,  !       QE/O(LNTE/OD,2),
c$$$     X          P_DLAM= 6*LEVS+2*LEVH+4,  !  DPDLAME/O(LNTE/OD,2),
c$$$     X          P_DPHI= 6*LEVS+2*LEVH+5,  !  DPDPHIE/O(LNTE/OD,2),
c$$$     X          P_ULN = 6*LEVS+2*LEVH+6,  !     ULNE/O(LNTE/OD,2,LEVS),
c$$$     X          P_VLN = 7*LEVS+2*LEVH+6,  !     VLNE/O(LNTE/OD,2,LEVS),
c$$$     X          P_W   = 8*LEVS+2*LEVH+6,  !       WE/O(LNTE/OD,2,LEVS),
c$$$     X          P_X   = 9*LEVS+2*LEVH+6,  !       XE/O(LNTE/OD,2,LEVS),
c$$$     X          P_Y   =10*LEVS+2*LEVH+6,  !       YE/O(LNTE/OD,2,LEVS),
c$$$     X          P_RT  =11*LEVS+2*LEVH+6,  !      RTE/O(LNTE/OD,2,LEVH),
c$$$     X          P_ZQ  =11*LEVS+3*LEVH+6)  !      ZQE/O(LNTE/OD,2)

      INTEGER   kdt,IERR,J,K,L,LOCL,N
      integer  idt,mdt,jdt,kdtdfi
      logical lsout
      REAL(KIND=KIND_EVOD) TEE1(LEVS)
cjfe
      real(kind=kind_evod)  ye1(levs)
      REAL(KIND=KIND_EVOD)  coef00(LEVS,ntrac) ! temp. ozone clwater  
      INTEGER              INDLSEV,JBASEV
      INTEGER              INDLSOD,JBASOD
      include 'function2'

      real ,allocatable ::  qse(:,:)
      real ,allocatable :: dise(:,:,:)
      real ,allocatable ::  zes(:,:,:)
      real ,allocatable ::  tes(:,:,:)
      real ,allocatable :: rqse(:,:,:)
!
      real ,allocatable ::  qso(:,:)
      real ,allocatable :: diso(:,:,:)
      real ,allocatable ::  zos(:,:,:)
      real ,allocatable ::  tos(:,:,:)
      real ,allocatable :: rqso(:,:,:)

!     REAL QSE(lnte,2),DISE(lnte,2,LEVS),
!    &     ZES(lnte,2,LEVS),
!    &     TES(lnte,2,LEVS),RQSE(lnte,2,levh)
!     REAL QSO(lnto,2),DISO(lnto,2,LEVS),
!    &     ZOS(lnto,2,LEVS),
!    &     TOS(lnto,2,LEVS),RQSO(lnto,2,levh)

      real totsum
      real deltim_loc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   following allocs. are for gridded tracers
      real ,allocatable :: rgt_s(:,:,:,:)
cc
       allocate( rgt_s(lonf,levs,lats_dim_a,ntrac))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (.not.allocated(rgt_a))
     &         allocate (rgt_a(lonf,levs,lats_dim_a,ntrac))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      allocate (  qse(len_trie_ls,2)      )
      allocate ( dise(len_trie_ls,2,levs) )
      allocate (  zes(len_trie_ls,2,levs) )
      allocate (  tes(len_trie_ls,2,levs) )
      allocate ( rqse(len_trie_ls,2,levh) )
!
      allocate (  qso(len_trio_ls,2)      )
      allocate ( diso(len_trio_ls,2,levs) )
      allocate (  zos(len_trio_ls,2,levs) )
      allocate (  tos(len_trio_ls,2,levs) )
      allocate ( rqso(len_trio_ls,2,levh) )

      rgt_a(:,:,:,:) = 0.0

      if (semilag) then
        deltim_loc = deltim
      else
        deltim_loc = deltim*0.5
      endif

!     print*,' beg tldfi- semilag,deltim_loc = ',semilag,deltim_loc,
!    . 'gg_tracers=',gg_tracers


!     print *,' enter tldfi ' 					! hmhj
!
!  INCLUDE FIRST TWO TIME LEVELS
!!
      KDTDFI = KDT + NSDFI
      if (comp_task) then
        CALL DFINI(-NSDFI-1  ,NSDFI,
     &   trie_ls(1,1,P_q),trie_ls(1,1,P_di),
     &   trie_ls(1,1,P_ze),trie_ls(1,1,P_te),trie_ls(1,1,P_rq),
     &   trio_ls(1,1,P_q),trio_ls(1,1,P_di),
     &   trio_ls(1,1,P_ze),trio_ls(1,1,P_te),trio_ls(1,1,P_rq),
     &   TOTSUM,QSE,DISE,ZES,TES,RQSE,QSO,DISO,ZOS,TOS,RQSO,
     &   rgt_s,gg_tracers)
!!
      CALL DFINI(KDT-KDTDFI,NSDFI,
     &   trie_ls(1,1,P_q),trie_ls(1,1,P_di),
     &   trie_ls(1,1,P_ze),trie_ls(1,1,P_te),trie_ls(1,1,P_rq),
     &   trio_ls(1,1,P_q),trio_ls(1,1,P_di),
     &   trio_ls(1,1,P_ze),trio_ls(1,1,P_te),trio_ls(1,1,P_rq),
     &   TOTSUM,QSE,DISE,ZES,TES,RQSE,QSO,DISO,ZOS,TOS,RQSO,
     &   rgt_s,gg_tracers)
      endif

      KDT   = KDT + 1
      FHOUR = KDT*DELTIM/3600
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      lssav = .true. !always true, except in digital filter
      lsswr = .true. !ex short wave radaition, used in gloopr(astronomy)
      lslwr = .true. !ex long  wave radaition, used in gloopr(astronomy)
      lsfwd = .true. !true only during forward step
      lscca = .false.!get clouds from precp.(first step use fixio_R clds)
      lsout = MOD(KDT,NSOUT).EQ.0 .OR. PHOUR.EQ.0.
      if (nsout_hf > 0 .and. phour <= fhmax_hf)                         &
     &   lsout = MOD(kdt ,NSOUT_hf) == 0 .OR. lsout
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
      if(hybrid)then
       call get_cd_hyb(deltim/2.)
      else if( gen_coord_hybrid ) then				! hmhj
       call get_cd_hyb_gc(deltim/2.)				! hmhj
      else
       call get_cd_sig(am,bm,deltim/2.,tov,sv)
      endif
CC
      if (me == 0) print*,'kdt before forward step in tldfi=',kdt
      CALL do_tstep(deltim_loc,kdt,PHOUR,
     &           TRIE_LS,TRIO_LS,
     &           LS_NODE,LS_NODES,MAX_LS_NODES,
     &           LATS_NODES_A,GLOBAL_LATS_A,
     &           LONSPERLAT,
     &           LATS_NODES_R,GLOBAL_LATS_R,
     &           LONSPERLAR,
!    &           LATS_NODES_EXT,GLOBAL_LATS_EXT,
     &           EPSE,EPSO,EPSEDN,EPSODN,
     &           SNNP1EV,SNNP1OD,NDEXEV,NDEXOD,
     X           PLNEV_A,PLNOD_A,PDDEV_A,PDDOD_A,
     X           PLNEW_A,PLNOW_A,
     X           PLNEV_R,PLNOD_R,PDDEV_R,PDDOD_R,
     X           PLNEW_R,PLNOW_R,
     X           SYN_LS_A,DYN_LS_A,
!    X           SYN_GR_A_1,DYN_GR_A_1,ANL_GR_A_1,
!    X           SYN_GR_A_2,DYN_GR_A_2,ANL_GR_A_2,
!    &           LSLAG,
     &           XLON,XLAT,COSZDG,sfc_fld,flx_fld,nst_fld,
     &           HPRIME,SWH,HLW,FLUXR,
     &           SFALB,SLAG,SDEC,CDEC,
     &           OZPLIN,JINDX1,JINDX2,
     &           DDY,PDRYINI,
     &           phy_f3d,  phy_f2d,
     &           ZHOUR,N1,N4,LSOUT,COLAT1,
     &           CFHOUR1,.false.,fscav)
CC
      if (me == 0) print*,'kdt after forward step in digifilter=',kdt
!!
      PHOUR = FHOUR
      if (comp_task) then
        CALL DFINI(KDT-KDTDFI,NSDFI,
     &   trie_ls(1,1,P_q),trie_ls(1,1,P_di),
     &   trie_ls(1,1,P_ze),trie_ls(1,1,P_te),trie_ls(1,1,P_rq),
     &   trio_ls(1,1,P_q),trio_ls(1,1,P_di),
     &   trio_ls(1,1,P_ze),trio_ls(1,1,P_te),trio_ls(1,1,P_rq),
     &   TOTSUM,QSE,DISE,ZES,TES,RQSE,QSO,DISO,ZOS,TOS,RQSO,
     &   rgt_s,gg_tracers)
      endif
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  INCLUDE TIME LEVELS UP TO HOUR FHOUR+FHDFI
c.....
      if(hybrid)then
       call get_cd_hyb(deltim)
      else if( gen_coord_hybrid ) then				! hmhj
       call get_cd_hyb_gc(deltim)				! hmhj
      else
       call get_cd_sig(am,bm,deltim,tov,sv)
      endif
c.....

      LSFWD = .FALSE.
      LSSAV = .TRUE.
      IDT   = KDT + 1
      MDT   = KDTDFI
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO JDT=IDT,MDT
        KDT=KDT+1
        FHOUR = KDT*DELTIM/3600
        LSOUT = MOD(KDT,NSOUT) == 0
        if (nsout_hf > 0 .and. fhour <= fhmax_hf)                         &
     &  lsout = MOD(kdt ,NSOUT_hf) == 0
        LSCCA = MOD(KDT,NSSWR) == 0
        LSSWR = MOD(KDT,NSSWR) == 1
        LSLWR = MOD(KDT,NSLWR) == 1
cccccccccccccccccccccccccccc
        CALL do_tstep(deltim,kdt,PHOUR,
     &           TRIE_LS,TRIO_LS,
     &           LS_NODE,LS_NODES,MAX_LS_NODES,
     &           LATS_NODES_A,GLOBAL_LATS_A,
     &           LONSPERLAT,
     &           LATS_NODES_R,GLOBAL_LATS_R,
     &           LONSPERLAR,
!    &           LATS_NODES_EXT,GLOBAL_LATS_EXT,
     &           EPSE,EPSO,EPSEDN,EPSODN,
     &           SNNP1EV,SNNP1OD,NDEXEV,NDEXOD,
     X           PLNEV_A,PLNOD_A,PDDEV_A,PDDOD_A,
     X           PLNEW_A,PLNOW_A,
     X           PLNEV_R,PLNOD_R,PDDEV_R,PDDOD_R,
     X           PLNEW_R,PLNOW_R,
     X           SYN_LS_A,DYN_LS_A,
!    X           SYN_GR_A_1,DYN_GR_A_1,ANL_GR_A_1,
!    X           SYN_GR_A_2,DYN_GR_A_2,ANL_GR_A_2,
!    &           LSLAG,
     &           XLON,XLAT,COSZDG,sfc_fld,flx_fld,nst_fld,
     &           HPRIME,SWH,HLW,FLUXR,
     &           SFALB,SLAG,SDEC,CDEC,
     &           OZPLIN,JINDX1,JINDX2,
     &           DDY,PDRYINI,
     &           phy_f3d,  phy_f2d,
     &           ZHOUR,N1,N4,LSOUT,COLAT1,
     &           CFHOUR1,.false.,fscav)
!!
      if (comp_task) then
        CALL DFINI(KDT-KDTDFI,NSDFI,
     &   trie_ls(1,1,P_q),trie_ls(1,1,P_di),
     &   trie_ls(1,1,P_ze),trie_ls(1,1,P_te),trie_ls(1,1,P_rq),
     &   trio_ls(1,1,P_q),trio_ls(1,1,P_di),
     &   trio_ls(1,1,P_ze),trio_ls(1,1,P_te),trio_ls(1,1,P_rq),
     &   TOTSUM,QSE,DISE,ZES,TES,RQSE,QSO,DISO,ZOS,TOS,RQSO,
     &   rgt_s,gg_tracers)
      endif
        PHOUR = FHOUR
      ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  SAVE SURFACE CONDITIONS
!!
      CALL synchro()
!!
!  SAVE SURFACE CONDITIONS
      if (comp_task) then
        call fixwr(1,

     &   nst_fld%xt,     nst_fld%xs,     nst_fld%xu,
     &   nst_fld%xv,     nst_fld%xz,     nst_fld%zm,
     &   nst_fld%xtts,   nst_fld%xzts,   nst_fld%dt_cool,
     &   nst_fld%z_c,    nst_fld%c_0,    nst_fld%c_d,
     &   nst_fld%w_0,    nst_fld%w_d,    nst_fld%d_conv,
     &   nst_fld%ifd,    nst_fld%tref,   nst_fld%Qrain,

     &   sfc_fld%hice,   sfc_fld%fice,   sfc_fld%tisfc, sfc_fld%tsea,  
     &   sfc_fld%smc,    sfc_fld%sheleg, sfc_fld%stc,   sfc_fld%tg3,
     &   sfc_fld%zorl,   sfc_fld%cv,     sfc_fld%cvb,   sfc_fld%cvt,
     &   sfc_fld%alvsf,  sfc_fld%alvwf,  sfc_fld%alnsf, sfc_fld%alnwf,
     &   sfc_fld%vfrac,  sfc_fld%canopy, sfc_fld%f10m,  sfc_fld%vtype,
     &   sfc_fld%stype,  sfc_fld%facsf,  sfc_fld%facwf, sfc_fld%uustar,
     &   sfc_fld%ffmm,   sfc_fld%ffhh,   sfc_fld%tprcp, sfc_fld%srflag,
     +   sfc_fld%slc,    sfc_fld%snwdph, sfc_fld%slope, sfc_fld%shdmin,
     &   sfc_fld%shdmax, sfc_fld%snoalb, sfc_fld%sncovr)
      endif
   
c......................................................................
C  INCLUDE TIME LEVELS UP TO HOUR FHOUR+2*FHDFI
C  BUT DO NOT SAVE DIAGNOSTICS FOR THIS TIME
      LSSAV = .FALSE.
      LSOUT = .FALSE.
      IDT   = KDT + 1
      MDT   = KDT + NSDFI
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DO JDT=IDT,MDT
        KDT=KDT+1
        FHOUR = KDT*DELTIM/3600
        LSCCA = MOD(KDT,NSSWR) == 0
        LSSWR = MOD(KDT,NSSWR) == 1
        LSLWR = MOD(KDT,NSLWR) == 1
cccccccccccccccccccccccccccc
        if (me == 0)
     &  print *, ' calling do_tstep in second loop kdt=',kdt,' me=',me
        CALL do_tstep(deltim,kdt,PHOUR,
     &           TRIE_LS,TRIO_LS,
     &           LS_NODE,LS_NODES,MAX_LS_NODES,
     &           LATS_NODES_A,GLOBAL_LATS_A,
     &           LONSPERLAT,
     &           LATS_NODES_R,GLOBAL_LATS_R,
     &           LONSPERLAR,
!    &           LATS_NODES_EXT,GLOBAL_LATS_EXT,
     &           EPSE,EPSO,EPSEDN,EPSODN,
     &           SNNP1EV,SNNP1OD,NDEXEV,NDEXOD,
     X           PLNEV_A,PLNOD_A,PDDEV_A,PDDOD_A,
     X           PLNEW_A,PLNOW_A,
     X           PLNEV_R,PLNOD_R,PDDEV_R,PDDOD_R,
     X           PLNEW_R,PLNOW_R,
     X           SYN_LS_A,DYN_LS_A,
!    X           SYN_GR_A_1,DYN_GR_A_1,ANL_GR_A_1,
!    X           SYN_GR_A_2,DYN_GR_A_2,ANL_GR_A_2,
!    &           LSLAG,
     &           XLON,XLAT,COSZDG,sfc_fld,flx_fld,nst_fld,
     &           HPRIME,SWH,HLW,FLUXR,
     &           SFALB,SLAG,SDEC,CDEC,
     &           OZPLIN,JINDX1,JINDX2,
     &           DDY,PDRYINI,
     &           phy_f3d,  phy_f2d,
     &           ZHOUR,N1,N4,LSOUT,COLAT1,
     &           CFHOUR1,.false.,fscav)
!!
!!
      if (comp_task) then
        CALL DFINI(KDT-KDTDFI,NSDFI,
     &   trie_ls(1,1,P_q),trie_ls(1,1,P_di),
     &   trie_ls(1,1,P_ze),trie_ls(1,1,P_te),trie_ls(1,1,P_rq),
     &   trio_ls(1,1,P_q),trio_ls(1,1,P_di),
     &   trio_ls(1,1,P_ze),trio_ls(1,1,P_te),trio_ls(1,1,P_rq),
     &   TOTSUM,QSE,DISE,ZES,TES,RQSE,QSO,DISO,ZOS,TOS,RQSO,
     &   rgt_s,gg_tracers)
       endif
       PHOUR = FHOUR
       if (me == 0)  print*,' fhour in second loop of tldfi kdt=',
     &                        fhour,kdt
      ENDDO
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  DO FINAL DIGITAL FILTER, set (n-1)=(n), and run forward step in main
      if (comp_task) then
        CALL DFINI(NSDFI+1,NSDFI,
     &   trie_ls(1,1,P_q),trie_ls(1,1,P_di),
     &   trie_ls(1,1,P_ze),trie_ls(1,1,P_te),trie_ls(1,1,P_rq),
     &   trio_ls(1,1,P_q),trio_ls(1,1,P_di),
     &   trio_ls(1,1,P_ze),trio_ls(1,1,P_te),trio_ls(1,1,P_rq),
     &   TOTSUM,QSE,DISE,ZES,TES,RQSE,QSO,DISO,ZOS,TOS,RQSO,
     &   rgt_s,gg_tracers)

        if (me == 0) print*,'kdt after last dfini in digifilter=',kdt
!!
        do i=1,len_trie_ls
          trie_ls(i,1,P_qm) = trie_ls(i,1,P_q)
          trie_ls(i,2,P_qm) = trie_ls(i,2,P_q)
        enddo
        do i=1,len_trio_ls
          trio_ls(i,1,p_qm) = trio_ls(i,1,p_q)
          trio_ls(i,2,p_qm) = trio_ls(i,2,p_q)
        enddo
!
        do k=1,levs
          do i=1,len_trie_ls
            trie_ls(i,1,P_tem+k-1) = trie_ls(i,1,P_te+k-1)
            trie_ls(i,1,P_dim+k-1) = trie_ls(i,1,P_di+k-1)
            trie_ls(i,1,P_zem+k-1) = trie_ls(i,1,P_ze+k-1)
            trie_ls(i,2,P_tem+k-1) = trie_ls(i,2,P_te+k-1)
            trie_ls(i,2,P_dim+k-1) = trie_ls(i,2,P_di+k-1)
            trie_ls(i,2,P_zem+k-1) = trie_ls(i,2,P_ze+k-1)
          enddo
          do i=1,len_trio_ls
            trio_ls(i,1,p_tem+k-1) = trio_ls(i,1,p_te+k-1)
            trio_ls(i,1,p_dim+k-1) = trio_ls(i,1,p_di+k-1)
            trio_ls(i,1,p_zem+k-1) = trio_ls(i,1,p_ze+k-1)
            trio_ls(i,2,p_tem+k-1) = trio_ls(i,2,p_te+k-1)
            trio_ls(i,2,p_dim+k-1) = trio_ls(i,2,p_di+k-1)
            trio_ls(i,2,p_zem+k-1) = trio_ls(i,2,p_ze+k-1)
          enddo
        enddo
!
        if(gg_tracers)then
!         in case tracers need to be moved to other holding area
        else
          do k=1,levh
            do i=1,len_trie_ls
              trie_ls(i,1,p_rm+k-1) = trie_ls(i,1,p_rq+k-1)
              trie_ls(i,2,p_rm+k-1) = trie_ls(i,2,p_rq+k-1)
            enddo
            do i=1,len_trio_ls
              trio_ls(i,1,p_rm+k-1) = trio_ls(i,1,p_rq+k-1)
              trio_ls(i,2,p_rm+k-1) = trio_ls(i,2,p_rq+k-1)
            enddo
          enddo
        endif
!!
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  REPLACE SURFACE CONDITIONS with conditions written at mid ini segment.
c  forward step in main begins with values in the middle of the filter span
c  forward step in main begins with values at time fhdfi
        call fixwr(2,

     &   nst_fld%xt,     nst_fld%xs,     nst_fld%xu,
     &   nst_fld%xv,     nst_fld%xz,     nst_fld%zm,
     &   nst_fld%xtts,   nst_fld%xzts,   nst_fld%dt_cool,
     &   nst_fld%z_c,    nst_fld%c_0,    nst_fld%c_d,
     &   nst_fld%w_0,    nst_fld%w_d,    nst_fld%d_conv,
     &   nst_fld%ifd,    nst_fld%tref,   nst_fld%Qrain,

     &   sfc_fld%hice,   sfc_fld%fice,   sfc_fld%tisfc, sfc_fld%tsea,
     &   sfc_fld%smc,    sfc_fld%sheleg, sfc_fld%stc,   sfc_fld%tg3,
     &   sfc_fld%zorl,   sfc_fld%cv,     sfc_fld%cvb,   sfc_fld%cvt,
     &   sfc_fld%alvsf,  sfc_fld%alvwf,  sfc_fld%alnsf, sfc_fld%alnwf,
     &   sfc_fld%vfrac,  sfc_fld%canopy, sfc_fld%f10m,  sfc_fld%vtype,
     &   sfc_fld%stype,  sfc_fld%facsf,  sfc_fld%facwf, sfc_fld%uustar,
     &   sfc_fld%ffmm,   sfc_fld%ffhh,   sfc_fld%tprcp, sfc_fld%srflag,
     +   sfc_fld%slc,    sfc_fld%snwdph, sfc_fld%slope, sfc_fld%shdmin,
     &   sfc_fld%shdmax, sfc_fld%snoalb, sfc_fld%sncovr)
       endif
!!
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  RESET CLOCK AND OUTPUT initialized fields
      KDT   = KDTDFI
      FHOUR = KDT*DELTIM/3600 ! note that fhour also comes from last fixio
      if(me == 0) print*,'fhour after reset clock digifilter=',fhour,kdt
      LSOUT = MOD(KDT,NSOUT) == 0
      if (nsout_hf > 0 .and. fhour <= fhmax_hf)                         &
     &  lsout = MOD(kdt ,NSOUT_hf) == 0
!!
      IF (me == 0) THEN
        write(*,*)'Initialized values in digifilter'
        write(*,*)'************'
c$$$        CALL bar3(trie_ls(1,1,P_ze),trio_ls(1,1,P_ze),'ze ',levs)
c$$$        CALL bar3(trie_ls(1,1,P_di),trio_ls(1,1,P_di),'di ',levs)
c$$$        CALL bar3(trie_ls(1,1,P_te),trio_ls(1,1,P_te),'te ',levs)
c$$$        CALL bar3(trie_ls(1,1,P_rq),trio_ls(1,1,P_rq),'rq ',levs)
c$$$        CALL bar3(trie_ls(1,1,P_rq+levs),trio_ls(1,1,P_rq+levs),
c$$$     &            'oz1 ',levs)
c$$$        CALL bar3(trie_ls(1,1,P_rq+2*levs),trio_ls(1,1,P_rq+2*levs),
c$$$     &            'oz2 ',levs)
c$$$        CALL bar3(trie_ls(1,1,P_q),trio_ls(1,1,P_q),'q ',1)
c$$$        CALL bar3(trie_ls(1,1,P_gz),trio_ls(1,1,P_gz),'gz ',1)
!       print*,'P_qm =',P_qm ,' P_rm =',P_rm 
!sela if (.NOT.LIOPE.or.icolor.ne.2) then
c$$$        CALL RMS_spect(TRIE_LS(1,1,P_QM ), TRIE_LS(1,1,P_DIM),
c$$$     X             TRIE_LS(1,1,P_TEM), TRIE_LS(1,1,P_ZEM),
c$$$     X             TRIE_LS(1,1,P_RM ),
c$$$     X             TRIO_LS(1,1,P_QM ), TRIO_LS(1,1,P_DIM),
c$$$     X             TRIO_LS(1,1,P_TEM), TRIO_LS(1,1,P_ZEM),
c$$$     X             TRIO_LS(1,1,P_RM ),
c$$$     X             LS_NODES,MAX_LS_NODES)
!sela endif
!---------------------------------------------------------------

      ENDIF
!!
      PHOUR = FHOUR
!--------------------------------------------
cmy reset digifilter switch to zero for activiation of reshuffling lats loopa
      fhdfi = 0
!
       deallocate( rgt_s)

      deallocate ( qse  )
      deallocate ( dise )
      deallocate ( zes  )
      deallocate ( tes  )
      deallocate ( rqse )
!
      deallocate ( qso  )
      deallocate ( diso )
      deallocate ( zos  )
      deallocate ( tos  )
      deallocate ( rqso )
!
      if (me == 0)       print *,' exit tldif kdt=',kdt

      RETURN
      END
!!
      subroutine dfini(kstep,nstep,qe,die,ze,te,rqe,
     &                qo,dio,zo,to,rqo,
     & totsum,qse,dise,zes,tes,rqse,qso,diso,zos,tos,rqso,
     & rgt_s,gg_tracers)
cc
      use resol_def           , only : levh,levs,lnte,lnto,lonf,ntrac
      use layout1             , only : lats_dim_a,lats_node_a,
     &                                 len_trie_ls,len_trio_ls,me
      use layout_grid_tracers , only : rgt_a
      implicit none
cc
       real rgt_s(lonf,levs,lats_dim_a,ntrac)
cc

cc
      real qe(len_trie_ls,2),die(len_trie_ls,2,levs),
     &     ze(len_trie_ls,2,levs),
     &     te(len_trie_ls,2,levs),rqe(len_trie_ls,2,levh)
      real qo(len_trio_ls,2),dio(len_trio_ls,2,levs),
     &     zo(len_trio_ls,2,levs),
     &     to(len_trio_ls,2,levs),rqo(len_trio_ls,2,levh)
!!
      integer len_trie,len_trio
      real digfil,sx,wx,totsumi
      real qse(len_trie_ls,2),dise(len_trie_ls,2,levs),
     &     zes(len_trie_ls,2,levs),
     &     tes(len_trie_ls,2,levs),rqse(len_trie_ls,2,levh)
      real qso(len_trio_ls,2),diso(len_trio_ls,2,levs),
     &     zos(len_trio_ls,2,levs),
     &     tos(len_trio_ls,2,levs),rqso(len_trio_ls,2,levh)
      integer levl,kl
ccc      save totsum,qse,dise,zes,tes,rqse,qso,diso,zos,tos,rqso
c
      real totsum
      integer i,k,nstep,kstep,lan,nt
      logical gg_tracers

      if (me == 0) print *,' enter dfini ','kstep=',kstep,
     &                     ' nstep=',nstep                     ! hmhj

      if(kstep.lt.-nstep) then !++++++++++++++++++++++++++++++++++
        totsum = 0
        qse    = 0.
        qso    = 0.
        dise   = 0.
        zes    = 0.
        tes    = 0.
        rqse   = 0.
        diso   = 0.
        zos    = 0.
        tos    = 0.
        rqso   = 0.
        rgt_s  = 0.
        return
      endif
!
      if(kstep.le.nstep) then !++++++++++++++++++++++++++++++
!sela  print*,'arrived at elseif(kstep.le.nstep)',ktstep
        if(kstep.ne.0) then  !--------------------------------
          sx     = acos(-1.)*kstep/nstep
          wx     = acos(-1.)*kstep/(nstep+1)
          digfil = sin(wx)/wx*sin(sx)/sx
          if(me.eq.0)then
            print*,'in dfini sx=',sx,'wx=',wx,'digfil=',digfil,
     &      'at kstep=',kstep
          endif
        else                 !--------------------------------
!sela     print*,'arrived at if(kstep.ne.0) in elseif(kstep.le.nstep),
!sela                 ktstep= ntstep=',ktstep,ntstep
          sx     = acos(-1.)*kstep/nstep
          wx     = acos(-1.)*kstep/(nstep+1)
          digfil=1
          if(me.eq.0)then
            print*,'in dfini sx=',sx,'wx=',wx,'digfil=',digfil,
     &      'at kstep=',kstep
          endif
        endif                !--------------------------------

        totsum = totsum + digfil
        do k=1,2
          do i=1,len_trie_ls
           qse(i,k)=qse(i,k)+digfil*qe(i,k)
          enddo
        enddo
        do k=1,2
          do i=1,len_trio_ls
           qso(i,k)=qso(i,k)+digfil*qo(i,k)
          enddo
        enddo
        do kl=1,levs
          do k=1,2
            do i=1,len_trie_ls
              dise(i,k,kl)=dise(i,k,kl)+digfil*die(i,k,kl)
              zes(i,k,kl)=zes(i,k,kl)+digfil*ze(i,k,kl)
              tes(i,k,kl)=tes(i,k,kl)+digfil*te(i,k,kl)
            enddo
            do i=1,len_trio_ls
              diso(i,k,kl)=diso(i,k,kl)+digfil*dio(i,k,kl)
              zos(i,k,kl)=zos(i,k,kl)+digfil*zo(i,k,kl)
              tos(i,k,kl)=tos(i,k,kl)+digfil*to(i,k,kl)
            enddo
          enddo
        enddo

        if(gg_tracers)then

      if (me == 0) print *,' kstep=',kstep,' digfil=',digfil,
     & ' totsum=',totsum,' in dfini for gg_tracers'
        do nt=1,ntrac
          do lan=1,lats_node_a   !sela begin lan loop
            do k=1,levs
              do i=1,lonf
!               rgt_s(i,k,lan,nt) = rgt_a(i,k,lan,nt)
!    &                            + rgt_s(i,k,lan,nt)
                rgt_s(i,k,lan,nt) = rgt_s(i,k,lan,nt)
     &                            + digfil*rgt_a(i,k,lan,nt)
              enddo
            enddo
          enddo  ! lan loop
        enddo  ! nt loop

        else

          do kl=1,levh
            do k=1,2
              do i=1,len_trie_ls
                rqse(i,k,kl)=rqse(i,k,kl)+digfil*rqe(i,k,kl)
              enddo
              do i=1,len_trio_ls
                rqso(i,k,kl)=rqso(i,k,kl)+digfil*rqo(i,k,kl)
              enddo
            enddo
          enddo
c
        endif   ! if(gg_tracers)then
        return
      endif !++++++++++++++++++++++++++++++++++++++++++++++++++++

!sela  print*,'arrived at (kstep.lt.-nstep) in dfini
!sela& ktstep= ntstep=',ktstep,ntstep

      totsumi = 1.0 / totsum

      if (me == 0) print *,' totsum=',totsum,' totsumi=',totsumi

      do k=1,2
        do i=1,len_trie_ls
          qe(i,k)  = qse(i,k)  * totsumi
        enddo
        do i=1,len_trio_ls
          qo(i,k)  = qso(i,k)  * totsumi
        enddo
      enddo
      do kl=1,levs
        do k=1,2
          do i=1,len_trie_ls
            die(i,k,kl) = dise(i,k,kl) * totsumi
            ze(i,k,kl)  = zes(i,k,kl)  * totsumi
            te(i,k,kl)  = tes(i,k,kl)  * totsumi
          enddo
          do i=1,len_trio_ls
            dio(i,k,kl) = diso(i,k,kl) * totsumi
            zo(i,k,kl)  = zos(i,k,kl)  * totsumi
            to(i,k,kl)  = tos(i,k,kl)  * totsumi
          enddo
        enddo
      enddo

      if(gg_tracers)then
        do nt=1,ntrac
          do lan=1,lats_node_a   !sela begin lan loop
            do k=1,levs
              do i=1,lonf
                rgt_a(i,k,lan,nt) = rgt_s(i,k,lan,nt)*totsumi
              enddo
            enddo
          enddo  ! lan loop
        enddo  ! nt loop
        if (me == 0) then
          print *,' rgt_a=',rgt_a(1,:,1,1)
        endif
      else
        do kl=1,levh
          do k=1,2
            do i=1,len_trie_ls
              rqe(i,k,kl) = rqse(i,k,kl) * totsumi
            enddo
            do i=1,len_trio_ls
              rqo(i,k,kl) = rqso(i,k,kl) * totsumi
            enddo
          enddo
        enddo
      endif ! if(gg_tracers)then

      if (me == 0) print *,' leave dfini '                            ! hmhj
      end
      SUBROUTINE fixwr(iflag,

     & xt,xs,xu,xv,xz,zm,xtts,xzts,dt_cool,z_c,
     & c_0,c_d,w_0,w_d,d_conv,ifd,Tref,Qrain,

     & hice,fice,tisfc,                                ! FOR SEA-ICE - XW Nov04
     & tsea,smc,sheleg,stc,tg3,zorl,cv,cvb,cvt,
     & alvsf,alvwf,alnsf,alnwf,vfrac,canopy,f10m,vtype,stype,
     & facsf,facwf,uustar,ffmm,ffhh,tprcp,srflag,
     & slc,snwdph,slope,shdmin,shdmax,snoalb,sncovr)

c
c***********************************************************************
c     PURPOSE:
c      save or retrieve fixed fields in digifilt
c
c***********************************************************************
c
      use resol_def
      use layout1
      implicit none
      integer iflag
      real SMC(lonr,lsoil,lats_node_r),STC(lonr,lsoil,lats_node_r),
     &     HICE(lonr,lats_node_r),FICE(lonr,lats_node_r),  ! FOR SEA-ICE - NOV04
     &     TISFC(lonr,lats_node_r),

! li added for NST components
     &     xt     (lonr,lats_node_r), xs   (lonr,lats_node_r),
     &     xu     (lonr,lats_node_r), xv   (lonr,lats_node_r),
     &     xz     (lonr,lats_node_r), zm   (lonr,lats_node_r),
     &     xtts   (lonr,lats_node_r), xzts (lonr,lats_node_r),
     &     dt_cool(lonr,lats_node_r), z_c  (lonr,lats_node_r),
     &     c_0    (lonr,lats_node_r), c_d  (lonr,lats_node_r),
     &     w_0    (lonr,lats_node_r), w_d  (lonr,lats_node_r),
     &     d_conv (lonr,lats_node_r), ifd  (lonr,lats_node_r),
     &     Tref   (lonr,lats_node_r), Qrain(lonr,lats_node_r),
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     &     TSEA  (lonr,lats_node_r),SHELEG(lonr,lats_node_r),
     &     TG3   (lonr,lats_node_r),
     &     ZORL  (lonr,lats_node_r),CV    (lonr,lats_node_r),
     &     CVB   (lonr,lats_node_r),
     &     CVT   (lonr,lats_node_r),ALVSF (lonr,lats_node_r),
     &     ALVWF (lonr,lats_node_r),
     &     ALNSF (lonr,lats_node_r),ALNWF (lonr,lats_node_r),
     &     SLMSK (lonr,lats_node_r),
     &     VFRAC (lonr,lats_node_r),CANOPY(lonr,lats_node_r),
     &     F10M  (lonr,lats_node_r),
     &     VTYPE (lonr,lats_node_r),STYPE (lonr,lats_node_r),
     &     FACSF (lonr,lats_node_r),
     &     FACWF (lonr,lats_node_r),UUSTAR(lonr,lats_node_r),
     &     FFMM  (lonr,lats_node_r),
     &     FFHH  (lonr,lats_node_r)
Clu [+5L]: add (tprcp,srflag),(slc,snwdph,snoalb,slope,shdmin,shdmax)
     +,    TPRCP (lonr,lats_node_r),SRFLAG(lonr,lats_node_r)
     +,    SLC    (lonr,lsoil,lats_node_r)
     +,    SNWDPH (lonr,lats_node_r)
     +,    SNOALB (lonr,lats_node_r),SLOPE (lonr,lats_node_r)
     +,    SHDMIN (lonr,lats_node_r),SHDMAX(lonr,lats_node_r)
     +,    SNCOVR (lonr,lats_node_r)

      real , allocatable :: SMC1(:,:,:),STC1(:,:,:),
     &  HICE1(:,:),FICE1(:,:),TISFC1(:,:),                   ! FOR SEA-ICE - XW Nov04

! li added for NST components
     &  xt1(:,:),xs1(:,:),xu1(:,:),xv1(:,:),xz1(:,:),zm1(:,:),
     &  xtts1(:,:),xzts1(:,:),dt_cool1(:,:),z_c1(:,:),
     &  c_01(:,:),c_d1(:,:),w_01(:,:),w_d1(:,:),
     &  d_conv1(:,:),ifd1(:,:),Tref1(:,:),Qrain1(:,:),
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     &  TSEA1(:,:),SHELEG1(:,:),TG31(:,:),
     &  ZORL1(:,:),CV1(:,:),CVB1(:,:),
     &  CVT1(:,:),ALVSF1(:,:),ALVWF1(:,:),
     &  ALNSF1(:,:),ALNWF1(:,:),SLMSK1(:,:),
     &  VFRAC1(:,:),CANOPY1(:,:),F10M1(:,:),
     &  VTYPE1(:,:),STYPE1(:,:),FACSF1(:,:),
     &  FACWF1(:,:),UUSTAR1(:,:),FFMM1(:,:),
     &  FFHH1(:,:)
Clu [+3L]: add (tprcp1,srflag1),(slc1,snwdph1,slope1,shdmin1,shdmax1,snoalb1)
     +, TPRCP1(:,:),SRFLAG1(:,:)
     +, SLC1(:,:,:),SNWDPH1(:,:),SLOPE1(:,:)
     +, SHDMIN1(:,:),SHDMAX1(:,:),SNOALB1(:,:), SNCOVR1(:,:)

      logical first
      data first/.true./
      save   first,SMC1,STC1,TSEA1,SHELEG1,TG31,ZORL1,CV1,CVB1,CVT1
      save   HICE1,FICE1,TISFC1                        ! FOR SEA-ICE - XW Nov04
!                                                      ! FOR NST   - XL Dec0r97
      save   xt1,xs1,xu1,xv1,xz1,zm1,xtts1,xzts1,
     &       dt_cool1,z_c1,c_01,c_d1,w_01,w_d1,
     &       d_conv1,ifd1,Tref1,Qrain1
      save   ALVSF1,ALVWF1,ALNSF1,ALNWF1,SLMSK1,VFRAC1,CANOPY1,F10M1
      save   VTYPE1,STYPE1,FACSF1,FACWF1,UUSTAR1,FFMM1,FFHH1
Clu [+2L]: save (tprcp1,srflag1),(slc1,snwdph1,slope1,shdmin1,shdmax1,snoalb1)
      save   TPRCP1,SRFLAG1
      save   SLC1,SNWDPH1,SLOPE1,SHDMIN1,SHDMAX1,SNOALB1,SNCOVR1
  
      integer i,j,k

c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
!     print *,' enter fixwr ' 					! hmhj
      if (first) then
        allocate (SMC1(lonr,lsoil,lats_node_r))
        allocate (STC1(lonr,lsoil,lats_node_r))
        allocate (HICE1(lonr,lats_node_r))              ! FOR SEA-ICE - XW Nov04
        allocate (FICE1(lonr,lats_node_r))              ! FOR SEA-ICE - XW Nov04
        allocate (TISFC1(lonr,lats_node_r))             ! FOR SEA-ICE - XW Nov04
                                                        ! FOR NST     - XL Dec09
        allocate (xt1(lonr,lats_node_r))
        allocate (xs1(lonr,lats_node_r))
        allocate (xu1(lonr,lats_node_r))
        allocate (xv1(lonr,lats_node_r))
        allocate (xz1(lonr,lats_node_r))
        allocate (zm1(lonr,lats_node_r))
        allocate (xtts1(lonr,lats_node_r))
        allocate (xzts1(lonr,lats_node_r))

        allocate (dt_cool1(lonr,lats_node_r))
        allocate (z_c1(lonr,lats_node_r))
        allocate (c_01(lonr,lats_node_r))
        allocate (c_d1(lonr,lats_node_r))
        allocate (w_01(lonr,lats_node_r))
        allocate (w_d1(lonr,lats_node_r))
        allocate (d_conv1(lonr,lats_node_r))
        allocate (ifd1(lonr,lats_node_r))
        allocate (Tref1(lonr,lats_node_r))
        allocate (Qrain1(lonr,lats_node_r))

        allocate (TSEA1(lonr,lats_node_r))
        allocate (SHELEG1(lonr,lats_node_r))
        allocate (TG31(lonr,lats_node_r))
        allocate (ZORL1(lonr,lats_node_r))
        allocate (CV1(lonr,lats_node_r))
        allocate (CVB1(lonr,lats_node_r))
        allocate (CVT1(lonr,lats_node_r))
        allocate (ALVSF1(lonr,lats_node_r))
        allocate (ALVWF1(lonr,lats_node_r))
        allocate (ALNSF1(lonr,lats_node_r))
        allocate (ALNWF1(lonr,lats_node_r))
        allocate (SLMSK1(lonr,lats_node_r))
        allocate (VFRAC1(lonr,lats_node_r))
        allocate (CANOPY1(lonr,lats_node_r))
        allocate (F10M1(lonr,lats_node_r))
        allocate (VTYPE1(lonr,lats_node_r))
        allocate (STYPE1(lonr,lats_node_r))
        allocate (FACSF1(lonr,lats_node_r))
        allocate (FACWF1(lonr,lats_node_r))
        allocate (UUSTAR1(lonr,lats_node_r))
        allocate (FFMM1(lonr,lats_node_r))
        allocate (FFHH1(lonr,lats_node_r))
Clu [+8L]: allocate (tprcp,srflag),(slc,snwdph,slope,shdmin,shdmax,snoalb)
        allocate (TPRCP1(lonr,lats_node_r))
        allocate (SRFLAG1(lonr,lats_node_r))
        allocate (SLC1(lonr,lsoil,lats_node_r))
        allocate (SNWDPH1(lonr,lats_node_r))
        allocate (SLOPE1(lonr,lats_node_r))
        allocate (SHDMIN1(lonr,lats_node_r))
        allocate (SHDMAX1(lonr,lats_node_r))
        allocate (SNOALB1(lonr,lats_node_r))
        allocate (SNCOVR1(lonr,lats_node_r))
        first = .false.
      endif
!
      if(iflag == 1) then
        do k=1,lsoil
          do j=1,lats_node_r
            do i=1,lonr
              smc1(i,k,j) = smc(i,k,j)
              stc1(i,k,j) = stc(i,k,j)
              slc1(i,k,j) = slc(i,k,j)        !! Clu [+1L]: slc -> slc1
            enddo
          enddo
        enddo
        do j=1,lats_node_r
          do i=1,lonr
            hice1(i,j)    = hice(i,j)                   ! FOR SEA-ICE - XW Nov04
            fice1(i,j)    = fice(i,j)                   ! FOR SEA-ICE - XW Nov04
            tisfc1(i,j)   = tisfc(i,j)                  ! FOR SEA-ICE - XW Nov04
                                                        ! For NST
            xs1(i,j)      = xs(i,j)
            xu1(i,j)      = xu(i,j)
            xv1(i,j)      = xv(i,j)
            xz1(i,j)      = xz(i,j)
            zm1(i,j)      = zm(i,j)
            xtts1(i,j)    = xtts(i,j)
            xzts1(i,j)    = xzts(i,j)

            dt_cool1(i,j) = dt_cool(i,j)
            z_c1(i,j)     = z_c(i,j)
            c_01(i,j)     = c_0(i,j)
            c_d1(i,j)     = c_d(i,j)
            w_01(i,j)     = w_0(i,j)
            w_d1(i,j)     = w_d(i,j)
            d_conv1(i,j)  = d_conv(i,j)
            Tref1(i,j)    = Tref(i,j)
            Qrain1(i,j)   = Qrain(i,j)

            tsea1(i,j)    = tsea(i,j)
            sheleg1(i,j)  = sheleg(i,j)
            tg31(i,j)     = tg3(i,j)
            zorl1(i,j)    = zorl(i,j)
            cv1(i,j)      = cv(i,j)
            cvb1(i,j)     = cvb(i,j)
            cvt1(i,j)     = cvt(i,j)
            alvsf1(i,j)   = alvsf(i,j)
            alvwf1(i,j)   = alvwf(i,j)
            alnsf1(i,j)   = alnsf(i,j)
            alnwf1(i,j)   = alnwf(i,j)
            slmsk1(i,j)   = slmsk(i,j)
            vfrac1(i,j)   = vfrac(i,j)
            canopy1(i,j)  = canopy(i,j)
            f10m1(i,j)    = f10m(i,j)
            vtype1(i,j)   = vtype(i,j)
            stype1(i,j)   = stype(i,j)
            facsf1(i,j)   = facsf(i,j)
            facwf1(i,j)   = facwf(i,j)
            uustar1(i,j)  = uustar(i,j)
            ffmm1(i,j)    = ffmm(i,j)
            ffhh1(i,j)    = ffhh(i,j)
Clu [+7L]: add (tprcp,srflag),(snwdph,slope,shdmin,shdmax,snoalb)
            tprcp1(i,j)   = tprcp(i,j)
            srflag1(i,j)  = srflag(i,j)
            snwdph1(i,j)  = snwdph(i,j)
            slope1(i,j)   = slope(i,j)
            shdmin1(i,j)  = shdmin(i,j)
            shdmax1(i,j)  = shdmax(i,j)
            snoalb1(i,j)  = snoalb(i,j)
            sncovr1(i,j)  = sncovr(i,j)
          enddo
        enddo
      elseif(iflag == 2) then
        do k=1,lsoil
          do j=1,lats_node_r
            do i=1,lonr
              smc(i,k,j) = smc1(i,k,j)
              stc(i,k,j) = stc1(i,k,j)
              slc(i,k,j) = slc1(i,k,j)          !! Clu [+1L]: slc1 -> slc
            enddo
          enddo
        enddo
        do j=1,lats_node_r
          do i=1,lonr
            hice(i,j)    = hice1(i,j)                 ! FOR SEA-ICE - XW Nov04
            fice(i,j)    = fice1(i,j)                 ! FOR SEA-ICE - XW Nov04
            tisfc(i,j)   = tisfc1(i,j)                ! FOR SEA-ICE - XW Nov04

            xs(i,j)      = xs1(i,j)
            xu(i,j)      = xu1(i,j)
            xv(i,j)      = xv1(i,j)
            xz(i,j)      = xz1(i,j)
            zm(i,j)      = zm1(i,j)
            xtts(i,j)    = xtts1(i,j)
            xzts(i,j)    = xzts1(i,j)

            dt_cool(i,j) = dt_cool1(i,j)
            z_c(i,j)     = z_c1(i,j)
            c_0(i,j)     = c_01(i,j)
            c_d(i,j)     = c_d1(i,j)
            w_0(i,j)     = w_01(i,j)
            w_d(i,j)     = w_d1(i,j)
            d_conv(i,j)  = d_conv1(i,j)
            Tref(i,j)    = Tref1(i,j)
            Qrain(i,j)   = Qrain1(i,j)

            tsea(i,j)    = tsea1(i,j)
            sheleg(i,j)  = sheleg1(i,j)
            tg3(i,j)     = tg31(i,j)
            zorl(i,j)    = zorl1(i,j)
            cv(i,j)      = cv1(i,j)
            cvb(i,j)     = cvb1(i,j)
            cvt(i,j)     = cvt1(i,j)
            alvsf(i,j)   = alvsf1(i,j)
            alvwf(i,j)   = alvwf1(i,j)
            alnsf(i,j)   = alnsf1(i,j)
            alnwf(i,j)   = alnwf1(i,j)
            slmsk(i,j)   = slmsk1(i,j)
            vfrac(i,j)   = vfrac1(i,j)
            canopy(i,j)  = canopy1(i,j)
            f10m(i,j)    = f10m1(i,j)
            vtype(i,j)   = vtype1(i,j)
            stype(i,j)   = stype1(i,j)
            facsf(i,j)   = facsf1(i,j)
            facwf(i,j)   = facwf1(i,j)
            uustar(i,j)  = uustar1(i,j)
            ffmm(i,j)    = ffmm1(i,j)
            ffhh(i,j)    = ffhh1(i,j)
Clu [+7L]: add (tprcp,srflag),(snwdph,slope,shdmin,shdmax,snoalb)
            tprcp(i,j)   = tprcp1(i,j)
            srflag(i,j)  = srflag1(i,j)
            snwdph(i,j)  = snwdph1(i,j)
            slope(i,j)   = slope1(i,j)
            shdmin(i,j)  = shdmin1(i,j)
            shdmax(i,j)  = shdmax1(i,j)
            snoalb(i,j)  = snoalb1(i,j)
            sncovr(i,j)  = sncovr1(i,j)
          enddo
        enddo
      endif
!     print *,' leave fixwr '
      return
      end
