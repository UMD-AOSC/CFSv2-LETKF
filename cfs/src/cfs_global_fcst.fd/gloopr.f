       subroutine gloopr
     x    (trie_ls,trio_ls,
     x     ls_node,ls_nodes,max_ls_nodes,
     x     lats_nodes_a,global_lats_a,
     x     lats_nodes_r,global_lats_r,
     x     lonsperlar,
     x     epse,epso,epsedn,epsodn,
     x     snnp1ev,snnp1od, plnev_r,plnod_r,
     x     pddev_r,pddod_r,
!    x     snnp1ev,snnp1od,ndexev,ndexod,
!    x     plnev_r,plnod_r,pddev_r,pddod_r,plnew_r,plnow_r,
     x     phour,
     &     xlon,xlat,coszdg,COSZEN,
     &     SLMSK,SHELEG,SNCOVR,SNOALB,ZORL,TSEA,HPRIME,
!lu [+1L]: extract snow-free albedo (SFALB)
     +     SFALB,
     &     ALVSF,ALNSF ,ALVWF ,ALNWF,FACSF ,FACWF,CV,CVT ,
     &     CVB  ,SWH,HLW,SFCNSW,SFCDLW,
     &     FICE ,TISFC, SFCDSW, sfcemis,                ! FOR SEA-ICE - XW Nov04
     &     TSFLW,FLUXR ,       phy_f3d,slag,sdec,cdec,KDT,
     &     global_times_r)
!!
!#include "f_hpm.h"
!
      USE MACHINE              ,     ONLY : kind_phys
      USE FUNCPHYS             ,     ONLY : fpkap
      USE PHYSCONS, FV => con_fvirt, rerth => con_rerth 	! hmhj

      use module_radiation_driver,   only : radinit, grrad
      use module_radiation_astronomy,only : astronomy
!
!! ---  for optional spectral band heating outputs
!!    use module_radsw_parameters,   only : NBDSW
!!    use module_radlw_parameters,   only : NBDLW
!
      use resol_def
      use layout1
      use layout_grid_tracers
      use gg_def
      use vert_def
      use date_def
      use namelist_def
      use coordinate_def					! hmhj
      use tracer_const						! hmhj
      use d3d_def , only : cldcov
      use mersenne_twister, only : random_setseed, random_index,        &
     &                             random_stat
!
      implicit none
      include 'mpif.h'
!
      real (kind=kind_phys), parameter :: QMIN =1.0e-10
      real (kind=kind_evod), parameter :: Typical_pgr = 95.0
      real (kind=kind_evod), parameter :: cons0 = 0.0,  cons2 = 2.0
!
!  --- ...  inputs:
      integer, intent(in) :: ls_node(LS_DIM,3), ls_nodes(LS_DIM,NODES), &
     &                       max_ls_nodes(NODES), lats_nodes_r(NODES),  &
     &                       global_lats_r(LATR), lonsperlar(LATR)

      integer, intent(in) :: lats_nodes_a(nodes), global_lats_a(latg)


      real(kind=kind_evod), dimension(LEN_TRIE_LS), intent(in) ::       &
     &                       epse, epsedn, snnp1ev

      real(kind=kind_evod), dimension(LEN_TRIO_LS), intent(in) ::       &
     &                       epso, epsodn, snnp1od

      real(kind=kind_evod), intent(in) :: plnev_r(LEN_TRIE_LS, LATR2)
      real(kind=kind_evod), intent(in) :: plnod_r(LEN_TRIO_LS, LATR2)

      real (kind=kind_phys), dimension(LONR,LATS_NODE_R), intent(in) :: &
     &                       xlon, xlat, slmsk, sheleg, zorl, tsea,     &
     &                       alvsf, alnsf, alvwf, alnwf, facsf, facwf,  &
     &                       cv, cvt, cvb, FICE, tisfc, sncovr, snoalb

      real (kind=kind_phys), intent(in) ::                              &
     &                    hprime(LONR,NMTVR,LATS_NODE_R), phour,        &
     &                    phy_f3d(LONR,LEVS,NUM_P3D,LATS_NODE_R)
!
!  --- ...  input and output:
      real(kind=kind_evod), intent(inout) ::                            &
     &                    trie_ls(LEN_TRIE_LS,2,11*LEVS+3*LEVH+6),      &
     &                    trio_ls(LEN_TRIO_LS,2,11*LEVS+3*LEVH+6)
      integer              ipt_ls                                       ! hmhj
      real(kind=kind_evod) reall                                        ! hmhj
      real(kind=kind_evod) rlcs2(jcap1)                                 ! hmhj


      real (kind=kind_phys), intent(inout) ::                           &
     &                    fluxr (LONR,NFXR,LATS_NODE_R)

!  --- ...  inputs but not used anymore:
      real(kind=kind_evod), intent(in) :: pddev_r(LEN_TRIE_LS,LATR2),   &
     &                                    pddod_r(LEN_TRIO_LS,LATR2)    &
!    &                                    plnew_r(LEN_TRIE_LS,LATR2),   &
!    &                                    plnow_r(LEN_TRIO_LS,LATR2)
!    &                                    syn_ls_r(4*LS_DIM,LOTS,LATR2)

!     integer, intent(in) :: ndexev(LEN_TRIE_LS), ndexod(LEN_TRIO_LS)
      integer, intent(in) :: KDT
!  --- ...  outputs:
      real(kind=kind_evod), intent(inout) ::                            &
     &                    global_times_r(LATR,NODES)
      real(kind=kind_evod) ::                                           &
     &                    for_gr_r_1(LONRX,LOTS,LATS_DIM_R),            &
     &                    dyn_gr_r_1(lonrx,lotd,lats_dim_r),            ! hmhj
!mjr &                    for_gr_r_2(LONRX,LOTS,LATS_DIM_R),
     &                    for_gr_r_2(LONR ,LOTS,LATS_DIM_R),
!mjr &                    dyn_gr_r_2(lonrx,lotd,lats_dim_r)             ! hmhj
     &                    dyn_gr_r_2(lonr ,lotd,lats_dim_r)             ! hmhj

      real (kind=kind_phys), intent(out) ::                             &
     &                    swh(LONR,LEVS,LATS_NODE_R),                   &
     &                    hlw(LONR,LEVS,LATS_NODE_R)

      real (kind=kind_phys),dimension(LONR,LATS_NODE_R), intent(out) :: &
     &                    coszdg, coszen, sfcnsw, sfcdlw, tsflw,        &
     &                    sfcdsw, SFALB, sfcemis

      real (kind=kind_phys), intent(out) :: slag, sdec, cdec

!! --- ...  optional spectral band heating rates
!!    real (kind=kind_phys), optional, intent(out) ::                   &
!!   &                 htrswb(NGPTC,LEVS,NBDSW,NBLCK,LATS_NODE_R),      &
!!   &                 htrlwb(NGPTC,LEVS,NBDLW,NBLCK,LATS_NODE_R)

!  --- ...  locals:
!     real(kind=kind_phys) :: prsl(NGPTC,LEVS),  prdel(NGPTC,LEVS),     &
      real(kind=kind_phys) :: prsl(NGPTC,LEVS),  prslk(NGPTC,LEVS),     &
     &                        prsi(NGPTC,LEVP1), prsik(NGPTC,LEVP1),    &
     &                        hlw_v(NGPTC,LEVS), swh_v(NGPTC,LEVS)

      real (kind=kind_phys) :: si_loc(LEVR+1)

      real (kind=kind_phys) ::                                           &
!    &                      gu(NGPTC,LEVS),   gv1(NGPTC,LEVS),           &
!    &                      gt(NGPTC,LEVR),   gd (NGPTC,LEVS),           &
     &                      gt(NGPTC,LEVR),   gq(NGPTC),                 &
     &                      gr(NGPTC,LEVR),   gr1(NGPTC,LEVR,NTRAC-1),   &
!    &                      gphi(NGPTC),      glam(NGPTC), gq(NGPTC),    &
     &                      gtv(NGPTC,LEVR)
!    &                      sumq(NGPTC,LEVR), xcp(NGPTC,LEVR),           &! hmhj
!    &                      gtv(NGPTC,LEVR),  gtvx(NGPTC,LEVR),          &! hmhj
!    &                      gtvy(NGPTC,LEVR)                              ! hmhj
!    &,                     vvel(ngptc,levs)

      real (kind=kind_phys), allocatable ::  sumq(:,:), xcp(:,:),        &
     &                                       gtvx(:,:), gtvy(:,:),       &
!    &                                                  gd(:,:),         &
     &                                       vvel(:,:), gd(:,:),         &
     &                                       gu(:,:),   gv1(:,:),        &
     &                                       gphi(:),   glam(:)

      real (kind=kind_phys) :: f_ice(NGPTC,LEVS), f_rain(NGPTC,LEVS),   &
     &                         r_rime(NGPTC,LEVS)

      real (kind=kind_phys) :: cldcov_v(NGPTC,LEVS), fluxr_v(NGPTC,NFXR)
      real (kind=kind_phys) :: flgmin_v(ngptc), work1, work2

      real (kind=kind_phys), dimension(LONR,LATS_NODE_R) ::             &
     &                         sinlat_v, coslat_v

      real (kind=kind_phys) :: rinc(5), dtsw, dtlw, solcon, raddt

      real (kind=kind_phys), save :: facoz

!     integer :: njeff, lon, lan, lat, iblk, lon_dim, lons_lat, istrt
      integer :: njeff, lon, lan, lat,                lons_lat
      integer :: idat(8), jdat(8), DAYS(13), iday, imon, midmon, id
      integer :: lmax, nlnsp(LATS_NODE_R)

!  ---  variables used for random number generator (thread safe mode)
      type (random_stat) :: stat
      integer :: numrdm(LONR*LATR*2), ixseed(LONR,LATS_NODE_R,2)
      integer :: ipseed, icsdlw(NGPTC), icsdsw(NGPTC)
      integer, parameter :: ipsdlim = 1.0e8      ! upper limit for random seeds

      integer, save :: icwp, k1oz, k2oz, midm, midp, ipsd0

!  ---  number of days in a month
      data DAYS / 31,28,31,30,31,30,31,31,30,31,30,31,30 /

!  --- ...  control parameters: 
!           (some of the them may be moved into model namelist)

!  ---  ICTM=yyyy#, controls time sensitive external data (e.g. CO2, solcon, aerosols, etc)
!     integer, parameter :: ICTM =   -2 ! same as ICTM=0, but add seasonal cycle from
!                                       ! climatology. no extrapolation.
!     integer, parameter :: ICTM =   -1 ! use user provided external data set for the
!                                       ! forecast time. no extrapolation.
!     integer, parameter :: ICTM =    0 ! use data at initial cond time, if not
!                                       ! available, use latest, no extrapolation.
!!    integer, parameter :: ICTM =    1 ! use data at the forecast time, if not
!                                       ! available, use latest and extrapolation.
!     integer, parameter :: ICTM =yyyy0 ! use yyyy data for the forecast time,
!                                       ! no further data extrapolation.
!     integer, parameter :: ICTM =yyyy1 ! use yyyy data for the fcst. if needed, do
!                                       ! extrapolation to match the fcst time.

!  ---  ISOL controls solar constant data source
!!    integer, parameter :: ISOL = 0   ! use prescribed solar constant
!     integer, parameter :: ISOL = 1   ! use varying solar const with 11-yr cycle

!  ---  ICO2 controls co2 data source for radiation
!     integer, parameter :: ICO2 = 0   ! prescribed global mean value (old opernl)
!!    integer, parameter :: ICO2 = 1   ! use obs co2 annual mean value only
!     integer, parameter :: ICO2 = 2   ! use obs co2 monthly data with 2-d variation

!  ---  IALB controls surface albedo for sw radiation
!!    integer, parameter :: IALB = 0   ! use climatology alb, based on sfc type
!     integer, parameter :: IALB = 1   ! use modis derived alb (to be developed)

!  ---  IEMS controls surface emissivity and sfc air/ground temp for lw radiation
!        ab: 2-digit control flags. a-for sfc temperature;  b-for emissivity
!!    integer, parameter :: IEMS = 00  ! same air/ground temp; fixed emis = 1.0
!!    integer, parameter :: IEMS = 01  ! same air/ground temp; varying veg typ based emis
!!    integer, parameter :: IEMS = 10  ! diff air/ground temp; fixed emis = 1.0
!!    integer, parameter :: IEMS = 11  ! diff air/ground temp; varying veg typ based emis

!  ---  IAER  controls aerosols scheme selections
!     Old definition
!     integer, parameter :: IAER  = 1  ! opac climatology, without volc forcing
!     integer, parameter :: IAER  =11  ! opac climatology, with volcanic forcing
!     integer, parameter :: IAER  = 2  ! gocart prognostic, without volc forcing
!     integer, parameter :: IAER  =12  ! gocart prognostic, with volcanic forcing
!     New definition in this code IAER = abc (a:volcanic; b:lw; c:sw)
!                             b, c values: (0:none; 1:opac; 2:gocart)
!  IAER =   0 --> no aerosol effect at all (volc, sw, lw)
!       =   1 --> only tropospheric sw aerosols, no trop-lw and volc
!       =  10 --> only tropospheric lw aerosols, no trop-sw and volc
!       =  11 --> both trop-sw and trop-lw aerosols, no volc
!       = 100 --> only strato-volc aeros, no trop-sw and trop-lw
!       = 101 --> only sw aeros (trop + volc), no lw aeros
!       = 110 --> only lw aeros (trop + volc), no sw aeros
!       = 111 --> both sw and lw aeros (trop + volc) 
!

!  ---  IOVR controls cloud overlapping method in radiation:
!     integer, parameter :: IOVR_SW = 0  ! sw: random overlap clouds
!!    integer, parameter :: IOVR_SW = 1  ! sw: max-random overlap clouds

!     integer, parameter :: IOVR_LW = 0  ! lw: random overlap clouds
!!    integer, parameter :: IOVR_LW = 1  ! lw: max-random overlap clouds

!  ---  ISUBC controls sub-column cloud approximation in radiation:
!     integer, parameter :: ISUBC_SW = 0 ! sw: without sub-col clds approx
!     integer, parameter :: ISUBC_SW = 1 ! sw: sub-col clds with prescribed seeds
!     integer, parameter :: ISUBC_SW = 2 ! sw: sub-col clds with random seeds

!     integer, parameter :: ISUBC_LW = 0 ! lw: without sub-col clds approx
!     integer, parameter :: ISUBC_LW = 1 ! lw: sub-col clds with prescribed seeds
!     integer, parameter :: ISUBC_LW = 2 ! lw: sub-col clds with random seeds

!  ---  iflip indicates model vertical index direction:
!     integer, parameter :: IFLIP = 0    ! virtical profile index from top to bottom
      integer, parameter :: IFLIP = 1    ! virtical profile index from bottom to top
!
!    The following parameters are from gbphys
!
      real (kind=kind_phys), parameter :: dxmax=-16.118095651,          &
     &                dxmin=-9.800790154, dxinv=1.0/(dxmax-dxmin)

      integer :: kr, kt, kd, kq, ku, kv, ierr, dimg, kx, ky
      integer :: i, j, k, n, nt
      integer :: kdtphi,kdtlam,ks                                ! hmhj

      logical :: change
      logical, save :: first, sas_shal
      data  first / .true. /
!
!  ---  for debug test use
      real (kind=kind_phys) :: temlon, temlat, alon, alat
      integer :: ipt
      logical :: lprnt

!  ---  timers:
      real*8 :: rtc, timer1, timer2
!
!===> *** ...  begin here
!
!!
!cmr  ls_node(1,1) ... ls_node(ls_max_node,1) : values of L
!cmr  ls_node(1,2) ... ls_node(ls_max_node,2) : values of jbasev
!cmr  ls_node(1,3) ... ls_node(ls_max_node,3) : values of jbasod
!!
!$$$      integer                lots,lotd,lota
!$$$cc
!$$$      parameter            ( lots = 5*levs+1*levh+3 )
!$$$      parameter            ( lotd = 6*levs+2*levh+0 )
!$$$      parameter            ( lota = 3*levs+1*levh+1 )
!!
!!
      integer              kap,kar,kat,kau,kav,kdrlam
      integer              ksd,ksplam,kspphi,ksq,ksr,kst
      integer              ksu,ksv,ksz,node,item,jtem
!!
!     real(kind=kind_evod) spdlat(levs,lats_dim_r)
!Moor real(kind=kind_phys) slk(levs)
!     real(kind=kind_evod) spdmax_node (levs)
!     real(kind=kind_evod) spdmax_nodes(levs,nodes)
!!
!!
!!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!!
!!
!!................................................................
!!  syn(1, 0*levs+0*levh+1, lan)  ze
!!  syn(1, 1*levs+0*levh+1, lan)  di
!!  syn(1, 2*levs+0*levh+1, lan)  te
!!  syn(1, 3*levs+0*levh+1, lan)  rq
!!  syn(1, 3*levs+1*levh+1, lan)  q
!!  syn(1, 3*levs+1*levh+2, lan)  dpdlam
!!  syn(1, 3*levs+1*levh+3, lan)  dpdphi
!!  syn(1, 3*levs+1*levh+4, lan)  uln
!!  syn(1, 4*levs+1*levh+4, lan)  vln
!!................................................................
!!  dyn(1, 0*levs+0*levh+1, lan)  d(t)/d(phi)
!!  dyn(1, 1*levs+0*levh+1, lan)  d(rq)/d(phi)
!!  dyn(1, 1*levs+1*levh+1, lan)  d(t)/d(lam)
!!  dyn(1, 2*levs+1*levh+1, lan)  d(rq)/d(lam)
!!  dyn(1, 2*levs+2*levh+1, lan)  d(u)/d(lam)
!!  dyn(1, 3*levs+2*levh+1, lan)  d(v)/d(lam)
!!  dyn(1, 4*levs+2*levh+1, lan)  d(u)/d(phi)
!!  dyn(1, 5*levs+2*levh+1, lan)  d(v)/d(phi)
!!................................................................
!!  anl(1, 0*levs+0*levh+1, lan)  w     dudt
!!  anl(1, 1*levs+0*levh+1, lan)  x     dvdt
!!  anl(1, 2*levs+0*levh+1, lan)  y     dtdt
!!  anl(1, 3*levs+0*levh+1, lan)  rt    drdt
!!  anl(1, 3*levs+1*levh+1, lan)  z     dqdt
!!................................................................
!!
!!
!$$$      parameter(ksz     =0*levs+0*levh+1,
!$$$     x          ksd     =1*levs+0*levh+1,
!$$$     x          kst     =2*levs+0*levh+1,
!$$$     x          ksr     =3*levs+0*levh+1,
!$$$     x          ksq     =3*levs+1*levh+1,
!$$$     x          ksplam  =3*levs+1*levh+2,
!$$$     x          kspphi  =3*levs+1*levh+3,
!$$$     x          ksu     =3*levs+1*levh+4,
!$$$     x          ksv     =4*levs+1*levh+4)
!!
!$$$      parameter(kdtphi  =0*levs+0*levh+1,
!$$$     x          kdrphi  =1*levs+0*levh+1,
!$$$     x          kdtlam  =1*levs+1*levh+1,
!$$$     x          kdrlam  =2*levs+1*levh+1,
!$$$     x          kdulam  =2*levs+2*levh+1,
!$$$     x          kdvlam  =3*levs+2*levh+1,
!$$$     x          kduphi  =4*levs+2*levh+1,
!$$$     x          kdvphi  =5*levs+2*levh+1)
!!
!$$$      parameter(kau     =0*levs+0*levh+1,
!$$$     x          kav     =1*levs+0*levh+1,
!$$$     x          kat     =2*levs+0*levh+1,
!$$$     x          kar     =3*levs+0*levh+1,
!$$$     x          kap     =3*levs+1*levh+1)
!!
!!
!$$$      integer   P_gz,P_zem,P_dim,P_tem,P_rm,P_qm
!$$$      integer   P_ze,P_di,P_te,P_rq,P_q,P_dlam,P_dphi,P_uln,P_vln
!$$$      integer   P_w,P_x,P_y,P_rt,P_zq
!$$$cc
!$$$cc                                               old common /comfspec/
!$$$      parameter(P_gz  = 0*levs+0*levh+1,  !      gze/o(lnte/od,2),
!$$$     x          P_zem = 0*levs+0*levh+2,  !     zeme/o(lnte/od,2,levs),
!$$$     x          P_dim = 1*levs+0*levh+2,  !     dime/o(lnte/od,2,levs),
!$$$     x          P_tem = 2*levs+0*levh+2,  !     teme/o(lnte/od,2,levs),
!$$$     x          P_rm  = 3*levs+0*levh+2,  !      rme/o(lnte/od,2,levh),
!$$$     x          P_qm  = 3*levs+1*levh+2,  !      qme/o(lnte/od,2),
!$$$     x          P_ze  = 3*levs+1*levh+3,  !      zee/o(lnte/od,2,levs),
!$$$     x          P_di  = 4*levs+1*levh+3,  !      die/o(lnte/od,2,levs),
!$$$     x          P_te  = 5*levs+1*levh+3,  !      tee/o(lnte/od,2,levs),
!$$$     x          P_rq  = 6*levs+1*levh+3,  !      rqe/o(lnte/od,2,levh),
!$$$     x          P_q   = 6*levs+2*levh+3,  !       qe/o(lnte/od,2),
!$$$     x          P_dlam= 6*levs+2*levh+4,  !  dpdlame/o(lnte/od,2),
!$$$     x          P_dphi= 6*levs+2*levh+5,  !  dpdphie/o(lnte/od,2),
!$$$     x          P_uln = 6*levs+2*levh+6,  !     ulne/o(lnte/od,2,levs),
!$$$     x          P_vln = 7*levs+2*levh+6,  !     vlne/o(lnte/od,2,levs),
!$$$     x          P_w   = 8*levs+2*levh+6,  !       we/o(lnte/od,2,levs),
!$$$     x          P_x   = 9*levs+2*levh+6,  !       xe/o(lnte/od,2,levs),
!$$$     x          P_y   =10*levs+2*levh+6,  !       ye/o(lnte/od,2,levs),
!$$$     x          P_rt  =11*levs+2*levh+6,  !      rte/o(lnte/od,2,levh),
!$$$     x          P_zq  =11*levs+3*levh+6)  !      zqe/o(lnte/od,2)
!!
!!
!     print *,' in gloopr vertcoord_id =',vertcoord_id

!
      ksz     =0*levs+0*levh+1
      ksd     =1*levs+0*levh+1
      kst     =2*levs+0*levh+1
      ksq     =3*levs+0*levh+1
      ksplam  =3*levs+0*levh+2
      kspphi  =3*levs+0*levh+3
      ksu     =3*levs+0*levh+4
      ksv     =4*levs+0*levh+4
      ksr     =5*levs+0*levh+4

      kdtphi  =0*levs+0*levh+1                          ! hmhj
      kdtlam  =1*levs+1*levh+1                          ! hmhj
!!
      idat = 0
      idat(1) = idate(4)
      idat(2) = idate(2)
      idat(3) = idate(3)
      idat(5) = idate(1)
      rinc = 0.
      rinc(2) = fhour
      call w3movdat(rinc, idat, jdat)
!
      if (ntoz .le. 0) then                ! Climatological Ozone!
!
!     if(me .eq. 0) WRITE (6,989) jdat(1),jdat(2),jdat(3),jdat(5)
! 989 FORMAT(' UPDATING OZONE FOR ', I4,I3,I3,I3)
!
        IDAY   = jdat(3)
        IMON   = jdat(2)
        MIDMON = DAYS(IMON)/2 + 1
        CHANGE = FIRST .OR.
     &          ( (IDAY .EQ. MIDMON) .AND. (jdat(5).EQ.0) )
!
        IF (CHANGE) THEN
            IF (IDAY .LT. MIDMON) THEN
               K1OZ = MOD(IMON+10,12) + 1
               MIDM = DAYS(K1OZ)/2 + 1
               K2OZ = IMON
               MIDP = DAYS(K1OZ) + MIDMON
            ELSE
               K1OZ = IMON
               MIDM = MIDMON
               K2OZ = MOD(IMON,12) + 1
               MIDP = DAYS(K2OZ)/2 + 1 + DAYS(K1OZ)
            ENDIF
        ENDIF
!
        IF (IDAY .LT. MIDMON) THEN
           ID = IDAY + DAYS(K1OZ)
        ELSE
           ID = IDAY
        ENDIF
        FACOZ = real (ID-MIDM) / real (MIDP-MIDM)
      endif
!
      if (first) then
        sas_shal = sashal .and. (.not. ras)
!
        if( hybrid.or.gen_coord_hybrid ) then                             ! hmhj

          if( gen_coord_hybrid ) then                                     ! hmhj
            si_loc(levr+1) = si(levp1)                                    ! hmhj
            do k=1,levr                                                   ! hmhj
              si_loc(k) = si(k)                                           ! hmhj
            enddo                                                         ! hmhj
          else                                                            ! hmhj
!  ---  get some sigma distribution for radiation-cloud initialization
!sela   si(k)=(ak5(k)+bk5(k)*Typical_pgr)/Typical_pgr   !ak(k) bk(k) go top to botto
            si_loc(levr+1)= ak5(1)/typical_pgr+bk5(1)
            do k=1,levr
              si_loc(levr+1-k)= ak5(levp1-levr+k)/typical_pgr
     &                        + bk5(levp1-levr+k)
            enddo
          endif
        else
          do k = 1, levr
            si_loc(k) = si(k)
          enddo
          si_loc(levr+1) = si(levp1)
        endif       ! end_if_hybrid

!  --- determin prognostic/diagnostic cloud scheme

        icwp   = 0
        if (NTCW > 0) icwp = 1

        if( thermodyn_id.eq.3 ) then
          if (.not. allocated(xcp)) allocate (xcp(ngptc,levr))
          if (.not. allocated(sumq)) allocate (sumq(ngptc,levr))
        endif
        if( ntcw <= 0 ) then
          if( gen_coord_hybrid .and. vertcoord_id == 3.) then
            if (.not. allocated(gtvx)) allocate (gtvx(ngptc,levs))
            if (.not. allocated(gtvy)) allocate (gtvy(ngptc,levs))
          endif
          if (.not. allocated(gu))   allocate (gu(ngptc,levs))
          if (.not. allocated(gv1))  allocate (gv1(ngptc,levs))
          if (.not. allocated(gd))   allocate (gd(ngptc,levs))
          if (.not. allocated(vvel)) allocate (vvel(ngptc,levs))
          if (.not. allocated(gphi)) allocate (gphi(ngptc))
          if (.not. allocated(glam)) allocate (glam(ngptc))
        endif

!  ---  generate initial permutation seed for random number generator

        if ( ISUBC_LW==2 .or. ISUBC_SW==2 ) then
          ipsd0 = 17*idate(1) + 43*idate(2) + 37*idate(3) + 23*idate(4)
          if ( me == 0 ) then
            print *,'  Radiation sub-cloud initial seed =',ipsd0,       &
     &              ' idate =',idate
          endif
        endif
           
        first = .false.
           
      endif         ! end_if_first
!
!===> *** ...  radiation initialization
!
      dtsw  = 3600.0 * fhswr
      dtlw  = 3600.0 * fhlwr
      raddt = min(dtsw, dtlw)
                                                                                                            
      call radinit                                                      &
!  ---  input:
     &     ( si_loc, LEVR, IFLIP, idat, jdat, ICTM, ISOL, ICO2,         &
     &       IAER, IALB, IEMS, ICWP, NUM_P3D, ISUBC_SW, ISUBC_LW,       &
     &       IOVR_SW, IOVR_LW, me )
!  ---  output: ( none )

      do j = 1, LATS_NODE_R
        k   = global_lats_r(IPT_LATS_NODE_R-1+j)
        nlnsp(j) = lonsperlar(k)

        do i = 1, nlnsp(j)
          sinlat_v(i,j) = sinlat_r(k)
          coslat_v(i,j) = coslat_r(k)
        enddo

!  ---  padding spaces left
        n = nlnsp(j)
        do while (n < LONR)
          n = n + 1
          sinlat_v(n,j) = 0.0
          coslat_v(n,j) = 0.0
        enddo
      enddo
!
!===> *** ...  astronomy for sw radiation calculation.
!
!     print *,' calling astronomy'
      call astronomy                                                    &
!  ---  inputs:
     &     ( sinlat_v, coslat_v, xlon, fhswr, jdat,                     &
     &       LONR, LATS_NODE_R, nlnsp, lsswr, me,                       &
!  ---  outputs:
     &       solcon, slag, sdec, cdec, coszen, coszdg                   &
     &      )
!     print *,' returned from astro'

!
!===> *** ...  generate 2-d random seeds array for sub-grid cloud-radiation
!
      if ( ISUBC_LW==2 .or. ISUBC_SW==2 ) then
        ipseed = mod(nint(100.0*sqrt(fhour*3600)), ipsdlim) + 1 + ipsd0

        call random_setseed                                             &
!  ---  inputs:
     &     ( ipseed,                                                    &
!  ---  outputs:
     &       stat                                                       &
     &      )
        call random_index                                               &
!  ---  inputs:
     &     ( ipsdlim,                                                   &
!  ---  outputs:
     &       numrdm, stat                                               &
     &     )

        do k = 1, 2
          do j = 1, lats_node_r
            lat = global_lats_r(ipt_lats_node_r-1+j)

            do i = 1, LONR
              ixseed(i,j,k) = numrdm(i+(lat-1)*LONR+(k-1)*LATR)
            enddo
          enddo
        enddo
      endif

!
!===> *** ...  spectrum to grid transformation for radiation calculation.
!     -----------------------------------
!!
      call f_hpmstart(61,"gr delnpe")
      call delnpe(trie_ls(1,1,P_q   ),
     x            trio_ls(1,1,P_dphi),
     x            trie_ls(1,1,P_dlam),
     x            epse,epso,ls_node)
      call f_hpmstop(61)
!!
      call f_hpmstart(62,"gr delnpo")
      call delnpo(trio_ls(1,1,P_q   ),
     x            trie_ls(1,1,P_dphi),
     x            trio_ls(1,1,P_dlam),
     x            epse,epso,ls_node)
      call f_hpmstop(62)
!!
!     print *,' after delnpeo'
!!
      call f_hpmstart(63,"gr dezouv dozeuv")
!
!$omp parallel do shared(trie_ls,trio_ls)
!$omp+shared(epsedn,epsodn,snnp1ev,snnp1od,ls_node)
!$omp+private(k)
      do k=1,levs
         call dezouv(trie_ls(1,1,P_di +k-1), trio_ls(1,1,P_ze +k-1),
     x               trie_ls(1,1,P_uln+k-1), trio_ls(1,1,P_vln+k-1),
     x               epsedn,epsodn,snnp1ev,snnp1od,ls_node)
!!
         call dozeuv(trio_ls(1,1,P_di +k-1), trie_ls(1,1,P_ze +k-1),
     x               trio_ls(1,1,P_uln+k-1), trie_ls(1,1,P_vln+k-1),
     x               epsedn,epsodn,snnp1ev,snnp1od,ls_node)
      enddo
      call f_hpmstop(63)
!!
!sela print*,'completed call to dztouv'
!!
!!mr  call mpi_barrier (mpi_comm_world,ierr)
!!
      CALL countperf(0,5,0.)
      CALL synctime()
      CALL countperf(1,5,0.)
!!
      dimg=0
      CALL countperf(0,1,0.)
!!
      call f_hpmstart(67,"gr sumfln")
!!
!sela print*,'begining  call to sumfln'

      call sumfln_slg_gg(trie_ls(1,1,P_ze),
     x                   trio_ls(1,1,P_ze),
     x            lat1s_r,
     x            plnev_r,plnod_r,
     x            5*levs+3,ls_node,latr2,
     x            lats_dim_r,lots,for_gr_r_1,
     x            ls_nodes,max_ls_nodes,
     x            lats_nodes_r,global_lats_r,
!mjr x            lats_node_r,ipt_lats_node_r,lon_dims_r,
     x            lats_node_r,ipt_lats_node_r,lon_dim_r,
     x            lonsperlar,lonrx,latr,0)
!
      if(.not. gg_tracers) then
!       tracers grid values will be set from layout_grid_traces
        call sumfln_slg_gg(trie_ls(1,1,P_rq),
     x                     trio_ls(1,1,P_rq),
     x                     lat1s_r,
     x                     plnev_r,plnod_r,
     x                     levh,ls_node,latr2,
     x                     lats_dim_r,lots,for_gr_r_1,
     x                     ls_nodes,max_ls_nodes,
     x                     lats_nodes_r,global_lats_r,
!mjr x                     lats_node_r,ipt_lats_node_r,lon_dims_r,
     x                     lats_node_r,ipt_lats_node_r,lon_dim_r,
     x                     lonsperlar,lonrx,latr,5*levs+3)
      endif

!     print*,'completed call to sumfln'
!sela print*,'completed call to sumfln'
      call f_hpmstop(67)
!!
      CALL countperf(1,1,0.)
!!
! -----------------------------------
      if( vertcoord_id == 3. ) then
! -----------------------------------
        CALL countperf(0,1,0.)                                            ! hmhj
!
        call f_hpmstart(68,"gr sumder2")                                  ! hmhj
!
        do lan=1,lats_node_r
          lat = global_lats_r(ipt_lats_node_r-1+lan)
          lmax = min(jcap,lonsperlar(lat)/2)
          if ( (lmax+1)*2+1 .le. lonsperlar(lat)+2 ) then
            do k=levs+1,4*levs+2*levh
               do i = (lmax+1)*2+1, lonsperlar(lat)+2
                  dyn_gr_r_1(i,k,lan) = cons0    !constant
               enddo
            enddo
          endif
        enddo
!
        call sumder2_slg(trie_ls(1,1,P_te),                             ! hmhj
     x                   trio_ls(1,1,P_te),                             ! hmhj
     x                   lat1s_r,                                       ! hmhj
     x                   pddev_r,pddod_r,                               ! hmhj
     x                   levs,ls_node,latr2,                            ! hmhj
     x                   lats_dim_r,lotd,                               ! hmhj
     x                   dyn_gr_r_1,                                    ! hmhj
     x                   ls_nodes,max_ls_nodes,                         ! hmhj
     x                   lats_nodes_r,global_lats_r,                    ! hmhj
!mjr x                   lats_node_r,ipt_lats_node_r,lon_dims_r,        ! hmhj
     x                   lats_node_r,ipt_lats_node_r,lon_dim_r,         ! hmhj
     x                   lonsperlar,lonrx,latr,0)                       ! hmhj
!
        call f_hpmstop(68)                                                ! hmhj
!
        CALL countperf(1,1,0.)                                            ! hmhj
! --------------------------------
      endif     ! vertcoord_id=3
! --------------------------------
!
 
!!mr  call mpi_barrier (mpi_comm_world,ierr)

      if(gg_tracers .and. shuff_lats_r) then
        print*,' gloopr mpi_tracers_a_to_b shuff_lats_r',shuff_lats_r
        call mpi_tracers_a_to_b(
     x       rgt_a,lats_nodes_a,global_lats_a,
     x       for_gr_r_2(1,1,1),
     x       lats_nodes_r,global_lats_r,ksr,0)
      endif ! gg_tracers .and.  shuff_lats_r
 
      do lan=1,lats_node_r
         timer1 = rtc()
!!
         lat = global_lats_r(ipt_lats_node_r-1+lan)
!!
!        lon_dim = lon_dims_r(lan)
!!
         lons_lat = lonsperlar(lat)

! -------------------------------------------------------
         if( gen_coord_hybrid .and. vertcoord_id.eq.3. ) then
! -------------------------------------------------------
!
           lmax   = min(jcap,lons_lat/2)                                ! hmhj
           ipt_ls = min(lat,latr-lat+1)                                 ! hmhj

           do i=1,lmax+1                                                ! hmhj
             if ( ipt_ls .ge. lat1s_r(i-1) ) then                       ! hmhj
               reall    = i-1                                           ! hmhj
               rlcs2(i) = reall*rcs2_r(ipt_ls)/rerth                    ! hmhj
             else                                                       ! hmhj
               rlcs2(i) = cons0     !constant                           ! hmhj
             endif                                                      ! hmhj
           enddo                                                        ! hmhj
!
!$omp parallel do private(k,i,item,jtem)
          do k=1,levs                                                   ! hmhj
            item = kdtlam-1+k
            jtem = kst   -1+K
            do i=1,lmax+1                                               ! hmhj
!
!             d(t)/d(lam)                                               ! hmhj
               dyn_gr_r_1(i+i-1,item,lan) = -for_gr_r_1(i+i,jtem,lan)
     &                                    *  rlcs2(i)                   ! hmhj
               dyn_gr_r_1(i+i,item,lan)   =  for_gr_r_1(i+i-1,jtem,lan)
     &                                    *  rlcs2(i)                   ! hmhj
            enddo                                                       ! hmhj
          enddo                                                         ! hmhj
! --------------------
        endif       ! gc and vertcoord_id=3
! ---------------------
!
!!
         CALL countperf(0,6,0.)
!sela    print*,' beginning call four2grid',lan
!        print*,' beginning call four2grid',lan
         CALL FOUR_TO_GRID(for_gr_r_1(1,1,lan),for_gr_r_2(1,1,lan),
!mjr &                     lon_dim  ,lon_dim    ,lons_lat,5*levs+3)
     &                     lon_dim_r,lonr,lons_lat,5*levs+3)

!        print*,' after first call four2grid',lan
      if(gg_tracers)then
!
!   set tracers grid values from layout_grid_tracers
!
         if (.not.shuff_lats_r) then
!          set for_gr_r_2 to rgt_a  from gloopa
           do nt=1,ntrac
             do k=1,levs
               item = KSR - 1 + k +(nt-1)*levs
               jtem = lats_node_a+1-lan
               do i=1,min(lonf,lons_lat)
                 for_gr_r_2(i,item,lan) = rgt_a(i,k,jtem,nt)
               enddo
             enddo
           enddo
         endif ! not shuff_lats_r

      else
!     print *,' begin second call to FOUR_TO_GRID in gloopr'
         CALL FOUR_TO_GRID(for_gr_r_1(1,KSR,lan),
     &                     for_gr_r_2(1,KSR,lan),
!mjr &                     lon_dim  ,lon_dim    ,lons_lat,levh)
     &                     lon_dim_r,lonr,lons_lat,levh)
      endif
!     print *,' after second call to FOUR_TO_GRID in gloopr'

! -------------------------------------------------------
        if( gen_coord_hybrid.and.vertcoord_id.eq.3. ) then              ! hmhj
! -------------------------------------------------------
          CALL FOUR_TO_GRID(dyn_gr_r_1(1,1,lan),dyn_gr_r_2(1,1,lan),    ! hmhj
!mjr &                      lon_dim  ,lon_dim    ,lons_lat,levs)        ! hmhj
     &                      lon_dim_r,lonr,lons_lat,levs)               ! hmhj
          CALL FOUR_TO_GRID(dyn_gr_r_1(1,KDTLAM,lan),                   ! hmhj
     &                      dyn_gr_r_2(1,KDTLAM,lan),                   ! hmhj
!mjr &                      lon_dim  ,lon_dim    ,lons_lat,levs)        ! hmhj
     &                      lon_dim_r,lonr,lons_lat,levs)               ! hmhj
! -------------------------
        endif       ! gc and vertcoord_id=3
! -------------------------

!        print*,' completed call four2grid lan=',lan
!sela    print*,' completed call four2grid lan=',lan
         CALL countperf(1,6,0.)
!!
        if( .not. gen_coord_hybrid ) then                               ! hmhj

          do k = 1, LEVS
            item = KSR-1+k
            jtem = KST-1+k
            do j = 1, lons_lat
              if (for_gr_r_2(j,item,lan) <= qmin) then
                 for_gr_r_2(j,item,lan) = qmin
              endif
              for_gr_r_2(j,jtem,lan) = for_gr_r_2(j,jtem,lan)
     &                               / (1.0 + FV*for_gr_r_2(j,item,lan))
            enddo
          enddo
!     print *,' now do sfc pressure for lan=',lan
          do j = 1, lons_lat
            for_gr_r_2(j,KSQ,lan) = exp( for_gr_r_2(j,KSQ,lan) )
          enddo
!     print *,' after sfc pressure for lan=',lan

        endif                                                           ! hmhj
!
        timer2 = rtc()
        global_times_r(lat,me+1) = timer2 - timer1

!$$$    print*,'timeloopr',me,timer1,timer2,global_times_r(lat,me+1)
 
!!
      enddo   !lan
!
      call f_hpmstart(69,"gr lat_loop2")
!
!===> *** ...  starting latitude loop
!
      do lan=1,lats_node_r
! 
         lat = global_lats_r(ipt_lats_node_r-1+lan)
!
         lons_lat = lonsperlar(lat)

!!
!$omp parallel do schedule(dynamic,1) private(lon,i,j,k)
!$omp+private(vvel,gu,gv1,gd,gt,gr,gr1,gq,gphi,glam)
!$omp+private(gtv,gtvx,gtvy,sumq,xcp)
!$omp+private(cldcov_v,fluxr_v,f_ice,f_rain,r_rime)
!$omp+private(prslk,prsl,prsik,prsi,flgmin_v,hlw_v,swh_v)
!$omp+private(njeff,n,item,jtem,ks,work1,work2)
!$omp+private(icsdsw,icsdlw)
!$omp+private(lprnt,ipt)
!!!$omp+private(temlon,temlat,lprnt,ipt)

        DO lon=1,lons_lat,NGPTC
!!
          NJEFF   = MIN(NGPTC,lons_lat-lon+1)
!!
          lprnt = .false.
!
!  --- ...  for debug test
!     alon = 236.25
!     alat = 56.189
!     alon = 97.5
!     alat = -6.66
!     ipt = 0
!     do i = 1, njeff
!       item = lon + i - 1
!       temlon = xlon(item,lan) * 57.29578
!       if (temlon < 0.0) temlon = temlon + 360.0
!       temlat = xlat(item,lan) * 57.29578
!       lprnt = abs(temlon-alon) < 1.1 .and. abs(temlat-alat) < 1.1
!    &        .and. kdt > 0
!       if ( lprnt ) then
!         ipt = i
!         print *,' ipt=',ipt,' lon=',lon,' lan=',lan
!         exit
!       endif
!     enddo
!     lprnt = .false.
!!
          if (ntcw <= 0) then
            do k = 1, LEVS
              do j = 1, njeff
                jtem = lon-1+j
                gu (j,k)  = for_gr_r_2(jtem,KSU-1+k,lan)
                gv1(j,k)  = for_gr_r_2(jtem,KSV-1+k,lan)
                gd (j,k)  = for_gr_r_2(jtem,KSD-1+k,lan)
              enddo
            enddo
          endif

          if( gen_coord_hybrid ) then                                    ! hmhj
            do k=1,levr                                                  ! hmhj
              do j=1,NJEFF                                               ! hmhj
                jtem = lon-1+j
                gtv (j,k) = for_gr_r_2(jtem,KST-1+k,lan)
                gr  (j,k) = max(qmin, for_gr_r_2(jtem,KSR-1+k,lan))
!               gt  (j,k) = gtv(j,k) / (1.+fv*gr(j,k))
              enddo                                                      ! hmhj
            enddo
! --------------------------------------
            if( vertcoord_id == 3. .and. ntcw <= 0 ) then
! --------------------------------------
              do k=1,levs                                                ! hmhj
                do j=1,NJEFF                                             ! hmhj
                  jtem = lon-1+j
                  gtvx(j,k) = dyn_gr_r_2(jtem,KDTLAM-1+k,lan)
                  gtvy(j,k) = dyn_gr_r_2(jtem,KDTPHI-1+k,lan)
                enddo                                                    ! hmhj
              enddo
! -----------------------------
            endif
! -----------------------------
            if( thermodyn_id.eq.3 ) then ! get dry temperature from enthalpy
              do k=1,levr						! hmhj
                do j=1,njeff			    			! hmhj
                  sumq(j,k) = 0.0					! hmhj
                  xcp(j,k)  = 0.0					! hmhj
                enddo
              enddo
              do i=1,ntrac						! hmhj
                if( cpi(i).ne.0.0 ) then				! hmhj
                  ks = ksr+(i-1)*levs					! hmhj
                  do k=1,levr						! hmhj
                    item = ks-1+k
                    do j=1,njeff					! hmhj
                      jtem = lon-1+j
                      sumq(j,k) = sumq(j,k) + for_gr_r_2(jtem,item,lan)	! hmhj
                      xcp(j,k)  = xcp(j,k)
     &                          + cpi(i)*for_gr_r_2(jtem,item,lan)	! hmhj
                    enddo						! hmhj
                  enddo							! hmhj
                endif							! hmhj
              enddo							! hmhj
              do k=1,levr						! hmhj
                do j=1,njeff						! hmhj
                  xcp(j,k) = (1.-sumq(j,k))*cpi(0) + xcp(j,k)  		! hmhj
                  gt(j,k)  = gtv(j,k) / xcp(j,k) 			! hmhj
                enddo							! hmhj
              enddo							! hmhj
            else if( thermodyn_id.le.1 ) then				! hmhj
! get dry temperature from virtual temperature				! hmhj
              do k=1,levr                                               ! hmhj
                do j=1,njeff                                            ! hmhj
                  gt(j,k) = gtv(j,k) / (1.+fv*gr(j,k))  	        ! hmhj
                enddo                                                   ! hmhj
              enddo							! hmhj
            else
! get dry temperature from dry temperature           		        ! hmhj
              do k=1,levr                                               ! hmhj
                do j=1,njeff                                            ! hmhj
                  gt(j,k) = gtv(j,k)                                    ! hmhj
                enddo                                                   ! hmhj
              enddo 
            endif

          else                                                          ! hmhj
!
            do k = 1, levr
              do j = 1, njeff
                 jtem = lon-1+j
                 gt(j,k) = for_gr_r_2(jtem,KST-1+k,lan)
                 gr(j,k) = for_gr_r_2(jtem,KSR-1+k,lan)
              enddo
            enddo

          endif
!
!       Remaining tracers
!
          do n = 1, NTRAC-1
            do k = 1, LEVR
              item = KSR-1+k+n*levs
              do j = 1, njeff
                gr1(j,k,n) = for_gr_r_2(lon-1+j,item,lan)
              enddo
            enddo
          enddo
          if (ntcw > 0) then
            do j = 1, njeff
              gq  (j) = for_gr_r_2(lon-1+j,KSQ,lan)
            enddo
          else
            do j = 1, njeff
              jtem = lon-1+j
              gq  (j) = for_gr_r_2(jtem,KSQ   ,lan)
              gphi(j) = for_gr_r_2(jtem,KSPPHI,lan)
              glam(j) = for_gr_r_2(jtem,KSPLAM,lan)
            enddo
          endif
!!
!  ---  vertical structure variables:   del,si,sl,prslk,prdel
!
          if( gen_coord_hybrid ) then                                    ! hmhj
            call  hyb2press_gc(njeff,ngptc,gq,gtv,prsi,prsl,prsik,prslk) ! hmhj
            if (ntcw <= 0)
     &      call omegtes_gc(njeff,ngptc,rcs2_r(min(lat,latr-lat+1)),     ! hmhj
     &                     gq,gphi,glam,gtv,gtvx,gtvy,gd,gu,gv1,vvel)    ! hmhj
          elseif (hybrid) then
            call  hyb2press(njeff,ngptc,gq, prsi, prsl,prsik, prslk)
            if (ntcw <= 0)
     &      call omegtes(njeff,ngptc,rcs2_r(min(lat,latr-lat+1)),
     &                   gq,gphi,glam,gd,gu,gv1,vvel)
!    &                   gq,gphi,glam,gd,gu,gv1,vvel,ngptc,lprnt,ipt)
          else
            call  sig2press(njeff,ngptc,gq,sl,si,slk,sik,
     &                                        prsi,prsl,prsik,prslk)
            CALL countperf(0,12,0.)
            if (ntcw <= 0)
     &      call omegast3(njeff,ngptc,levs,
     &                    gphi,glam,gu,gv1,gd,del,
     &                    rcs2_r(min(lat,latr-lat+1)),vvel,gq,sl)
          endif
!.....
          if (levr .lt. levs) then
            do j=1,njeff
              prsi(j,levr+1)  = prsi(j,levp1)
              prsl(j,levr)    = (prsi(j,levp1) + prsi(j,levr)) * 0.5
              prsik(j,levr+1) = prslk(j,levp1)
              prslk(j,levr)   = fpkap(prsl(j,levr)*1000.0)
            enddo
          endif
!
          do k=1,nfxr
            do j=1,njeff
              fluxr_v(j,k) = fluxr(lon+j-1,k,lan)
            enddo
          enddo
!
          if (num_p3d == 3) then
            do k = 1, LEVR
              do j = 1, njeff
                jtem = lon+j-1
                f_ice (j,k) = phy_f3d(jtem,k,1,lan)
                f_rain(j,k) = phy_f3d(jtem,k,2,lan)
                r_rime(j,k) = phy_f3d(jtem,k,3,lan)
              enddo
            enddo

            work1 = (log(coslat_r(lat) / (lons_lat*latg)) - dxmin)*dxinv
            work1 = max(0.0, min(1.0,work1))
            work2 = flgmin(1)*work1 + flgmin(2)*(1.0-work1)
            do j=1,njeff
              flgmin_v(j) = work2
            enddo
          else
            do j=1,njeff
              flgmin_v(j) = 0.0
            enddo
          endif

!  *** ...  assign random seeds for sw and lw radiations

          if ( ISUBC_LW==2 .or. ISUBC_SW==2 ) then
            do j = 1, njeff
              icsdsw(j) = ixseed(lon+j-1,lan,1)
              icsdlw(j) = ixseed(lon+j-1,lan,2)
            enddo
          endif
!
!  *** ...  calling radiation driver

!
!     lprnt = me .eq. 0 .and. kdt .ge. 120
!     if (lprnt) then
!     if (kdt .gt. 85) then
!     print *,' calling grrad for me=',me,' lan=',lan,' lat=',lat
!    &,' num_p3d=',num_p3d
!      if (lan == 47) print *,' gt=',gt(1,:)
!      if (kdt > 3) call mpi_quit(5555)

!

          call grrad                                                    &
!  ---  inputs:
     &     ( prsi,prsl,prslk,gt,gr,gr1,vvel,slmsk(lon,lan),             &
     &       xlon(lon,lan),xlat(lon,lan),tsea(lon,lan),                 &
     &       sheleg(lon,lan),sncovr(lon,lan),snoalb(lon,lan),           &
     &       zorl(lon,lan),hprime(lon,1,lan),                           &
     &       alvsf(lon,lan),alnsf(lon,lan),alvwf(lon,lan),              &
     &       alnwf(lon,lan),facsf(lon,lan),facwf(lon,lan),              &
                                          ! fice FOR SEA-ICE XW Nov04
     &       fice(lon,lan),tisfc(lon,lan),                              &
     &       solcon,coszen(lon,lan),coszdg(lon,lan),k1oz,k2oz,facoz,    &
     &       cv(lon,lan),cvt(lon,lan),cvb(lon,lan),                     &
     &       IOVR_SW,IOVR_LW,f_ice,f_rain,r_rime,flgmin_v,              &
     &       icsdsw,icsdlw,NUM_P3D,NTCW-1,NCLD,NTOZ-1,NTRAC-1,NFXR,     &
     &       dtlw,dtsw,lsswr,lslwr,lssav,sas_shal,norad_precip,         &
     &       crick_proof, ccnorm,                                       &
     &       ngptc,njeff,LEVR,IFLIP, me, lprnt,ipt,kdt,                 &
!    &       ngptc,njeff,LEVR,IFLIP, me, lprnt,                         &
!  ---  outputs:
     &       swh_v,sfcnsw(lon,lan),sfcdsw(lon,lan),                     &
     &       sfalb(lon,lan),                                            &
     &       hlw_v,sfcdlw(lon,lan),tsflw(lon,lan),                      &
     &       sfcemis(lon,lan),cldcov_v,                                 &
!  ---  input/output:
     &       fluxr_v                                                    &
     &       )
!
!
!     if (lprnt) print *,' returned from grrad for me=',me,' lan=',
!    &lan,' lat=',lat,' kdt=',kdt
!     print *,' end gloopr HLW=',hlw(lon,:,lan),' lan=',lan
!
!        if (lprnt) print *,' swh_vg=',swh_v(ipt,:)
!
!
          if (lssav) then
            if (ldiag3d .or. lggfs3d) then
              do k=1,levr
                do j=1,njeff
                  cldcov(lon+j-1,k,lan) = cldcov(lon+j-1,k,lan)         &
     &                                  + cldcov_v(j,k) * raddt
                enddo
              enddo
            endif
          endif
          do k=1,nfxr
            do j=1,njeff
              fluxr(lon+j-1,k,lan) = fluxr_v(j,k)
            enddo
          enddo
          if (lslwr) then
            do k=1,levr
              do j=1,njeff
                jtem = lon + j - 1
                hlw(jtem,k,lan) = hlw_v(j,k)
                swh(jtem,k,lan) = swh_v(j,k)
              enddo
            enddo
            if (levr .lt. levs) then
              do k=levr+1,levs
                do j=1,njeff
                  jtem = lon + j - 1
                  hlw(jtem,k,lan) = hlw_v(j,levr)
                  swh(jtem,k,lan) = swh_v(j,levr)
                enddo
              enddo
            endif
          endif
!
!         if (lat == 45 .and. me == 0 .and. lon == 1) then
!           print *,' after grrad hlw_v=',hlw_v(1,:)
!           print *,' after grrad swh_v=',swh_v(1,:)
!         endif
!        if (lprnt) print *,' hlwg=',hlw(lon+ipt-1,:,lan)
!        if (lprnt) print *,' swhg=',swh(lon+ipt-1,:,lan)
!        if (lprnt) print *,' swh_vg=',swh_v(ipt,:)
 
!$$$          write(2900+lat,*) ' ilon = ',istrt
c$$$          write(2900+lat,'("swh",T16,"hlw")')
!$$$      do k=1,levs
!$$$         write(2900+lat,
!$$$     .         '(e10.3,T16,e10.3,T31,e10.3)')
!$$$     . swh(1,k,iblk,lan),hlw(1,k,iblk,lan)
!$$$       enddo
 
!!
!     print *,' completed grrad for lan=',lan,' istrt=',istrt
          CALL countperf(1,12,0.)
          ENDDO
!
      enddo
!!
      call f_hpmstop(69)
!!
      CALL countperf(0,5,0.)
      CALL synctime()
      CALL countperf(1,5,0.)
!sela print*,'completed gloopr_v kdt=',kdt
!!
      return
      end subroutine gloopr
