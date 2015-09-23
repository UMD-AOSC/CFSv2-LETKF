      subroutine gloopb
     x    (trie_ls,trio_ls,
     x     ls_node,ls_nodes,max_ls_nodes,
     x     lats_nodes_a,global_lats_a,
     x     lats_nodes_r,global_lats_r,
     x     lonsperlar,
     x     epse,epso,epsedn,epsodn,
     x     snnp1ev,snnp1od,ndexev,ndexod,
     x     plnev_r,plnod_r,pddev_r,pddod_r,plnew_r,plnow_r,
     &     tstep,phour,sfc_fld, flx_fld, nst_fld, SFALB,
     &     xlon,
     &     swh,hlw,hprime,slag,sdec,cdec,
     &     ozplin,jindx1,jindx2,ddy,pdryini,
     &     phy_f3d, phy_f2d,xlat,kdt,
     &     global_times_b,batah,lsout,fscav)
!!
!#include "f_hpm.h"
!!
      use machine             , only : kind_evod,kind_phys,kind_rad
      use resol_def           , only : jcap,jcap1,latg,latr,latr2,
     &                                 levh,levp1,levs,lnt2,
     &                                 lonf,lonr,lonrx,lota,lotd,lots,
     &                                 lsoil,ncld,nmtvr,nrcm,ntcw,ntoz,
     &                                 ntrac,num_p2d,num_p3d,
     &                                 p_di,p_dlam,p_dphi,p_q,
     &                                 p_rq,p_rt,p_te,p_uln,p_vln,
     &                                 p_w,p_x,p_y,p_ze,p_zq,
     &                                 thermodyn_id,sfcpress_id,nfxr

      use layout1             , only : ipt_lats_node_r,
     &                                 lat1s_r,lats_dim_r,
     &                                 lats_node_a,lats_node_r,
     &                                 len_trie_ls,len_trio_ls,
     &                                 lon_dim_r,ls_dim,ls_max_node,
     &                                 me,me_l_0,nodes
      use layout_grid_tracers , only : rgt_a,xhalo,
     &                                 rgt_h,yhalo
      use gg_def              , only : coslat_r,rcs2_r,sinlat_r,wgt_r
      use vert_def            , only : am,bm,del,si,sik,sl,slk,sv
      use date_def            , only : fhour,idate
      use namelist_def        , only : crtrh,fhswr,flgmin,
     &                                 gen_coord_hybrid,gg_tracers,
     &                                 hybrid,ldiag3d,lscca,lsfwd,
     &                                 lsm,lssav,lsswr,ncw,ngptc,
     &                                 old_monin,pre_rad,random_clds,
     &                                 ras,semilag,shuff_lats_r,
     &                                 sashal,ctei_rm,mom4ice,newsas,
     &                                 ccwf,cnvgwd,lggfs3d,trans_trac,
     &                                 mstrat,cal_pre,nst_fcst,dtphys,
     &                                 dlqf,moist_adj,cdmbgwd,
     &                                 bkgd_vdif_m, bkgd_vdif_h,
     &                                 bkgd_vdif_s,shal_cnv,
     &                                 psautco, prautco, evpco, wminco
      use coordinate_def      , only : ak5,bk5,vertcoord_id               ! hmhj
      use bfilt_def           , only : bfilte,bfilto
      use module_ras          , only : ras_init
      use physcons            , only :  grav => con_g,
     &                                 rerth => con_rerth,   ! hmhj
     &                                    fv => con_fvirt,   ! mjr
     &                                 rvrdm1 => con_FVirt,
     &                                    rd => con_rd
      use ozne_def            , only : latsozp,levozp,
     &                                 pl_coeff,pl_pres,timeoz
!-> Coupling insertion
      USE SURFACE_cc
!<- Coupling insertion

      use Sfc_Flx_ESMFMod
      use Nst_Var_ESMFMod
      use mersenne_twister
      use d3d_def
      use tracer_const
!
      implicit none
      include 'mpif.h'
!
      TYPE(Sfc_Var_Data)        :: sfc_fld
      TYPE(Flx_Var_Data)        :: flx_fld
      TYPE(Nst_Var_Data)        :: nst_fld
!
      real(kind=kind_phys), PARAMETER :: RLAPSE=0.65E-2
      real(kind=kind_evod), parameter :: cons_0=0.0,   cons_24=24.0
     &,                                  cons_99=99.0, cons_1p0d9=1.0E9


!$$$      integer n1rac, n2rac,nlons_v(ngptc)
!$$$      parameter (n1rac=ntrac-ntshft-1, n2rac=n1rac+1)
!
!     integer id,njeff,istrt,lon,kdt
      integer id,njeff,      lon,kdt
!!
      real(kind=kind_evod) save_qe_ls(len_trie_ls,2)
      real(kind=kind_evod) save_qo_ls(len_trio_ls,2)

      real(kind=kind_evod) sum_k_rqchange_ls(len_trie_ls,2)
      real(kind=kind_evod) sum_k_rqchango_ls(len_trio_ls,2)
!!
      real(kind=kind_evod) trie_ls(len_trie_ls,2,11*levs+3*levh+6)
      real(kind=kind_evod) trio_ls(len_trio_ls,2,11*levs+3*levh+6)
      real(kind=kind_evod) trie_ls_rqt(len_trie_ls,2,levs)
      real(kind=kind_evod) trio_ls_rqt(len_trio_ls,2,levs)
      real(kind=kind_evod) trie_ls_sfc(len_trie_ls,2)                   ! hmhj
      real(kind=kind_evod) trio_ls_sfc(len_trio_ls,2)                   ! hmhj
!!
      real(kind=kind_phys)    typdel(levs), batah
      real(kind=kind_phys)    prsl(ngptc,levs)
      real(kind=kind_phys)   prslk(ngptc,levs),dpshc(ngptc)
      real(kind=kind_phys)    prsi(ngptc,levs+1),phii(ngptc,levs+1)
      real(kind=kind_phys)   prsik(ngptc,levs+1),phil(ngptc,levs)
!!
      real (kind=kind_phys) gu(ngptc,levs),  gv1(ngptc,levs)
      real (kind=kind_phys) ugrd(ngptc,levs),vgrd(ngptc,levs)
      real (kind=kind_phys) gphi(ngptc),     glam(ngptc)
      real (kind=kind_phys) gq(ngptc),       gt(ngptc,levs), pgr(ngptc)
      real (kind=kind_phys) gr(ngptc,levs,ntrac)
      real (kind=kind_phys) gd(ngptc,levs)
      real (kind=kind_phys) adt(ngptc,levs), adr(ngptc,levs,ntrac)
      real (kind=kind_phys) adu(ngptc,levs), adv(ngptc,levs)
      real (kind=kind_phys) gtv(ngptc,levs)                              ! hmhj
      real (kind=kind_phys) gtvx(ngptc,levs), gtvy(ngptc,levs)           ! hmhj
      real (kind=kind_phys) sumq(ngptc,levs), xcp(ngptc,levs)
!
      real (kind=kind_phys) dt3dt_v(ngptc,levs,6),
     &                      dq3dt_v(ngptc,levs,5+pl_coeff),
     &                      du3dt_v(ngptc,levs,4),
     &                      dv3dt_v(ngptc,levs,4)
     &,                     upd_mf_v(ngptc,levs), dwn_mf_v(ngptc,levs)
     &,                     det_mf_v(ngptc,levs)
     &,                     dkh_v(ngptc,LEVS),    rnp_v(ngptc,levs)
!!
      real(kind=kind_evod) gq_save(lonr,lats_dim_r)
!!
      real (kind=kind_rad) slag,sdec,cdec,phour
      real (kind=kind_rad) xlon(lonr,lats_node_r)
      real (kind=kind_rad) xlat(lonr,lats_node_r)
      real (kind=kind_rad) coszdg(lonr,lats_node_r),
     &                     hprime(lonr,nmtvr,lats_node_r),
!    &                     fluxr(lonr,nfxr,lats_node_r),
     &                     sfalb(lonr,lats_node_r)
      real (kind=kind_rad) swh(lonr,levs,lats_node_r)
      real (kind=kind_rad) hlw(lonr,levs,lats_node_r)
!!
      real  (kind=kind_phys)
     &     phy_f3d(lonr,levs,num_p3d,lats_node_r),
     &     phy_f2d(lonr,num_p2d,lats_node_r), fscav(ntrac-ncld-1)
!
!
!     real (kind=kind_phys) exp,dtphys,dtp,dtf,sumed(2)
      real (kind=kind_phys) exp,       dtp,dtf,sumed(2)
      real (kind=kind_evod) tstep
      real (kind=kind_phys) pdryini,sigshc,rk
!!
      integer              ls_node(ls_dim,3)
cc
!     ls_node(1,1) ... ls_node(ls_max_node,1) : values of l
!     ls_node(1,2) ... ls_node(ls_max_node,2) : values of jbasev
!     ls_node(1,3) ... ls_node(ls_max_node,3) : values of jbasod
cc
      integer              ls_nodes(ls_dim,nodes)
cc
      integer              max_ls_nodes(nodes)
      integer              lats_nodes_a(nodes)
      integer              lats_nodes_r(nodes)
cc
      integer              global_lats_a(latg)
      integer              global_lats_r(latr)
      integer                 lonsperlar(latr)
      integer dimg
cc
      real(kind=kind_evod)    epse(len_trie_ls)
      real(kind=kind_evod)    epso(len_trio_ls)
      real(kind=kind_evod)  epsedn(len_trie_ls)
      real(kind=kind_evod)  epsodn(len_trio_ls)
cc
      real(kind=kind_evod) snnp1ev(len_trie_ls)
      real(kind=kind_evod) snnp1od(len_trio_ls)
cc
      integer               ndexev(len_trie_ls)
      integer               ndexod(len_trio_ls)
cc
      real(kind=kind_evod)   plnev_r(len_trie_ls,latr2)
      real(kind=kind_evod)   plnod_r(len_trio_ls,latr2)
      real(kind=kind_evod)   pddev_r(len_trie_ls,latr2)
      real(kind=kind_evod)   pddod_r(len_trio_ls,latr2)
      real(kind=kind_evod)   plnew_r(len_trie_ls,latr2)
      real(kind=kind_evod)   plnow_r(len_trio_ls,latr2)
cc
c$$$      integer                lots,lotd,lota
      integer lotn
c$$$cc
c$$$      parameter            ( lots = 5*levs+1*levh+3 )
c$$$      parameter            ( lotd = 6*levs+2*levh+0 )
c$$$      parameter            ( lota = 3*levs+1*levh+1 )
cc
      real(kind=kind_evod) for_gr_r_1(lonrx,lots,lats_dim_r)
      real(kind=kind_evod) dyn_gr_r_1(lonrx,lotd,lats_dim_r)            ! hmhj
      real(kind=kind_evod) bak_gr_r_1(lonrx,lota,lats_dim_r)
cc
!     real(kind=kind_evod) for_gr_r_2(lonrx*lots,lats_dim_r)
!     real(kind=kind_evod) dyn_gr_r_2(lonrx*lotd,lats_dim_r)            ! hmhj
!     real(kind=kind_evod) bak_gr_r_2(lonrx*lota,lats_dim_r)
!
      real(kind=kind_evod) for_gr_r_2(lonr,lots,lats_dim_r)
      real(kind=kind_evod) dyn_gr_r_2(lonr,lotd,lats_dim_r)             ! hmhj
      real(kind=kind_evod) bak_gr_r_2(lonr,lota,lats_dim_r)
cc
      integer              i,ierr,iter,j,k,kap,kar,kat,kau,kav,ksq,jj,kk
      integer              kst,kdtphi,kdtlam                            ! hmhj
      integer              l,lan,lan0,lat,lmax,locl,ii,lonrbm
!     integer              lon_dim,lons_lat,n,node
      integer                      lons_lat,n,node
      integer nsphys
!
      real(kind=kind_evod) pwatg(latr),pwatj(lats_node_r),
     &                     pwatp,ptotg(latr),sumwa,sumto,
     &                     ptotj(lats_node_r),pcorr,pdryg,
     &                     solhr,clstp
cc
      integer              ipt_ls                                       ! hmhj
      real(kind=kind_evod) reall                                        ! hmhj
      real(kind=kind_evod) rlcs2(jcap1)                                 ! hmhj
 
       real(kind=kind_evod) typical_pgr
c
!timers______________________________________________________---
 
      real*8 rtc ,timer1,timer2
      real(kind=kind_evod) global_times_b(latr,nodes)
 
!timers______________________________________________________---
cc
cc
c$$$      parameter(ksq     =0*levs+0*levh+1,
c$$$     x          ksplam  =0*levs+0*levh+2,
c$$$     x          kspphi  =0*levs+0*levh+3,
c$$$     x          ksu     =0*levs+0*levh+4,
c$$$     x          ksv     =1*levs+0*levh+4,
c$$$     x          ksz     =2*levs+0*levh+4,
c$$$     x          ksd     =3*levs+0*levh+4,
c$$$     x          kst     =4*levs+0*levh+4,
c$$$     x          ksr     =5*levs+0*levh+4)
cc
c$$$      parameter(kdtphi  =0*levs+0*levh+1,
c$$$     x          kdrphi  =1*levs+0*levh+1,
c$$$     x          kdtlam  =1*levs+1*levh+1,
c$$$     x          kdrlam  =2*levs+1*levh+1,
c$$$     x          kdulam  =2*levs+2*levh+1,
c$$$     x          kdvlam  =3*levs+2*levh+1,
c$$$     x          kduphi  =4*levs+2*levh+1,
c$$$     x          kdvphi  =5*levs+2*levh+1)
cc
c$$$      parameter(kau     =0*levs+0*levh+1,
c$$$     x          kav     =1*levs+0*levh+1,
c$$$     x          kat     =2*levs+0*levh+1,
c$$$     x          kar     =3*levs+0*levh+1,
c$$$     x          kap     =3*levs+1*levh+1)
cc
cc
c$$$      integer   p_gz,p_zem,p_dim,p_tem,p_rm,p_qm
c$$$      integer   p_ze,p_di,p_te,p_rq,p_q,p_dlam,p_dphi,p_uln,p_vln
c$$$      integer   p_w,p_x,p_y,p_rt,p_zq
c$$$cc
c$$$cc                                               old common /comfspec/
c$$$      parameter(p_gz  = 0*levs+0*levh+1,  !      gze/o(lnte/od,2),
c$$$     x          p_zem = 0*levs+0*levh+2,  !     zeme/o(lnte/od,2,levs),
c$$$     x          p_dim = 1*levs+0*levh+2,  !     dime/o(lnte/od,2,levs),
c$$$     x          p_tem = 2*levs+0*levh+2,  !     teme/o(lnte/od,2,levs),
c$$$     x          p_rm  = 3*levs+0*levh+2,  !      rme/o(lnte/od,2,levh),
c$$$     x          p_qm  = 3*levs+1*levh+2,  !      qme/o(lnte/od,2),
c$$$     x          p_ze  = 3*levs+1*levh+3,  !      zee/o(lnte/od,2,levs),
c$$$     x          p_di  = 4*levs+1*levh+3,  !      die/o(lnte/od,2,levs),
c$$$     x          p_te  = 5*levs+1*levh+3,  !      tee/o(lnte/od,2,levs),
c$$$     x          p_rq  = 6*levs+1*levh+3,  !      rqe/o(lnte/od,2,levh),
c$$$     x          p_q   = 6*levs+2*levh+3,  !       qe/o(lnte/od,2),
c$$$     x          p_dlam= 6*levs+2*levh+4,  !  dpdlame/o(lnte/od,2),
c$$$     x          p_dphi= 6*levs+2*levh+5,  !  dpdphie/o(lnte/od,2),
c$$$     x          p_uln = 6*levs+2*levh+6,  !     ulne/o(lnte/od,2,levs),
c$$$     x          p_vln = 7*levs+2*levh+6,  !     vlne/o(lnte/od,2,levs),
c$$$     x          p_w   = 8*levs+2*levh+6,  !       we/o(lnte/od,2,levs),
c$$$     x          p_x   = 9*levs+2*levh+6,  !       xe/o(lnte/od,2,levs),
c$$$     x          p_y   =10*levs+2*levh+6,  !       ye/o(lnte/od,2,levs),
c$$$     x          p_rt  =11*levs+2*levh+6,  !      rte/o(lnte/od,2,levh),
c$$$     x          p_zq  =11*levs+3*levh+6)  !      zqe/o(lnte/od,2)
cc
cc
      integer              indlsev,jbasev,n0
      integer              indlsod,jbasod
cc
      include 'function2'
cc
      real(kind=kind_evod) cons0,cons2     !constant
cc
      logical lsout
      logical, parameter :: flipv = .true.
cc
! for nasa/nrl ozone production and distruction rates:(input through fixio)
      real ozplin(latsozp,levozp,pl_coeff,timeoz)
      integer jindx1(lats_node_r),jindx2(lats_node_r)!for ozone interpolaton
      real ddy(lats_node_r)                          !for ozone interpolaton
      real ozplout(levozp,lats_node_r,pl_coeff)
!Moor real ozplout(lonr,levozp,pl_coeff,lats_node_r)
!!
      real(kind=kind_phys), allocatable :: acv(:,:),acvb(:,:),acvt(:,:)
      save acv,acvb,acvt
!!
!     integer, parameter :: maxran=5000
!     integer, parameter :: maxran=3000
!     integer, parameter :: maxran=6000, maxsub=6, maxrs=maxran/maxsub
      integer, parameter :: maxran=3000, maxsub=6, maxrs=maxran/maxsub
      type (random_stat) :: stat
      real (kind=kind_phys), allocatable, save :: rannum_tank(:,:,:)
      real (kind=kind_phys), allocatable       :: rannum(:)
!     integer iseed, nrc, seed0, kss, ksr, indxr(nrcm), iseedl, latseed
      integer iseed, nrc, seed0, kss, ksr, indxr(nrcm), iseedl, latseed
      integer nf0,nf1,ind,nt,indod,indev
      real(kind=kind_evod) fd2, wrk(1), wrk2(nrcm)

      logical first,ladj
      parameter (ladj=.true.)
      data first/.true./
!     save    krsize, first, nrnd,seed0
      save    first, seed0
!!
      integer nlons_v(ngptc)
      real(kind=kind_phys) smc_v(ngptc,lsoil),stc_v(ngptc,lsoil)
     &,                    slc_v(ngptc,lsoil)
     &,                    swh_v(ngptc,levs), hlw_v(ngptc,levs)
     &,                    vvel(ngptc,levs)
     &,                    hprime_v(ngptc,nmtvr)
      real(kind=kind_phys) phy_f3dv(ngptc,LEVS,num_p3d),
     &                     phy_f2dv(ngptc,num_p2d)
     &,                    rannum_v(ngptc,nrcm)
      real(kind=kind_phys) sinlat_v(ngptc),coslat_v(ngptc)
     &,                    ozplout_v(ngptc,levozp,pl_coeff)
      real (kind=kind_rad) rqtk(ngptc), rcs2_lan, rcs_lan
!     real (kind=kind_rad) rcs_v(ngptc), rqtk(ngptc), rcs2_lan
!
!--------------------------------------------------------------------
!     print *,' in gloopb vertcoord_id =',vertcoord_id

!     real(kind=kind_evod) sinlat_v(lonr),coslat_v(lonr),rcs2_v(lonr)
!     real(kind=kind_phys) dpshc(lonr)
      real (kind=kind_rad) work1, qmin, tem
      parameter (qmin=1.0e-10)
      integer              ksd,ksplam,kspphi
      integer              ksu,ksv,ksz,item,jtem,ktem,ltem,mtem
      integer              nn, nnr, nnrcm

!
!  ---  for debug test use
      real (kind=kind_phys) :: temlon, temlat, alon, alat
      integer :: ipt
      logical :: lprnt
!
!
      ksq     =0*levs+0*levh+1
      ksplam  =0*levs+0*levh+2
      kspphi  =0*levs+0*levh+3
      ksu     =0*levs+0*levh+4
      ksv     =1*levs+0*levh+4
      ksz     =2*levs+0*levh+4
      ksd     =3*levs+0*levh+4
      kst     =4*levs+0*levh+4
      ksr     =5*levs+0*levh+4
!
      kau     =0*levs+0*levh+1
      kav     =1*levs+0*levh+1
      kat     =2*levs+0*levh+1
      kap     =3*levs+0*levh+1
      kar     =3*levs+0*levh+2
!
!     ksq     = 0*levs + 0*levh + 1
!     kst     = 4*levs + 0*levh + 4                          ! hmhj
      kdtphi  = 0*levs + 0*levh + 1                          ! hmhj
      kdtlam  = 1*levs+1*levh+1                              ! hmhj
!
c$$$      kau     =0*levs+0*levh+1
c$$$      kav     =1*levs+0*levh+1
c$$$      kat     =2*levs+0*levh+1
c$$$      kar     =3*levs+0*levh+1
c$$$      kap     =3*levs+1*levh+1
cc
cc--------------------------------------------------------------------
cc
      save_qe_ls(:,:) = trie_ls(:,:,p_q)
      save_qo_ls(:,:) = trio_ls(:,:,p_q)


!
      if (first) then
        allocate (bfilte(lnt2),bfilto(lnt2))
!
!     initializations for the gloopb filter
!     *************************************
        nf0 = (jcap+1)*2/3  ! highest wavenumber gloopb filter keeps fully
        nf1 = (jcap+1)      ! lowest wavenumber gloopb filter removes fully
        fd2 = 1./(nf1-nf0)**2
        do locl=1,ls_max_node
             l   = ls_node(locl,1)
          jbasev = ls_node(locl,2)
!sela     if (l.eq.0) then
!sela       n0=2
!sela     else
!sela       n0=l
!sela     endif
!
!sela     indev = indlsev(n0,l)
          indev = indlsev(l,l)
!sela     do n=n0,jcap1,2
!mjr      do n=l,jcap1,2
          do n=l,jcap,2
            bfilte(indev) = max(1.-fd2*max(n-nf0,0)**2,cons_0)     !constant
            indev         = indev + 1
          enddo
          if (mod(L,2).eq.mod(jcap+1,2)) bfilte(indev) = 1.
        enddo
!!
        do locl=1,ls_max_node
              l  = ls_node(locl,1)
          jbasod = ls_node(locl,3)
          indod  = indlsod(l+1,l)
!mjr      do n=l+1,jcap1,2
          do n=l+1,jcap,2
            bfilto(indod) = max(1.-fd2*max(n-nf0,0)**2,cons_0)     !constant
            indod = indod+1
          enddo
          if (mod(L,2).ne.mod(jcap+1,2)) bfilto(indod) = 1.
        enddo
!!
!       call random_seed(size=krsize)
!       if (me.eq.0) print *,' krsize=',krsize
!       allocate (nrnd(krsize))
        allocate (acv(lonr,lats_node_r))
        allocate (acvb(lonr,lats_node_r))
        allocate (acvt(lonr,lats_node_r))
!
!
        if (.not. newsas) then  ! random number needed for RAS and old SAS
          if (random_clds) then ! create random number tank
!                                 -------------------------
            seed0 = idate(1) + idate(2) + idate(3) + idate(4)
            call random_setseed(seed0)
            call random_number(wrk)
            seed0 = seed0 + nint(wrk(1)*1000.0)
            if (me == 0) print *,' seed0=',seed0,' idate=',idate,
     &                           ' wrk=',wrk
!
            if (.not. allocated(rannum_tank))
     &                allocate (rannum_tank(lonr,maxran,lats_node_r))
            if (.not. allocated(rannum)) allocate (rannum(lonr*maxrs))
            lonrbm = lonr / maxsub
            if (me == 0) write(0,*)' maxran=',maxran,' maxrs=',maxrs,
     &          'maxsub=',maxsub,' lonrbm=',lonrbm,
     &          ' lats_node_r=',lats_node_r
            do j=1,lats_node_r
              iseedl = global_lats_r(ipt_lats_node_r-1+j) + seed0
              call random_setseed(iseedl)
              call random_number(rannum)
!$omp parallel do  shared(j,lonr,lonrbm,rannum,rannum_tank)
!$omp+private(nrc,nn,i,ii,k,kk)
              do nrc=1,maxrs
                nn = (nrc-1)*lonr
                do k=1,maxsub
                  kk = k - 1
                  do i=1,lonr
                    ii = kk*lonrbm + i
                    if (ii > lonr) ii = ii - lonr
                    rannum_tank(i,nrc+kk*maxrs,j) = rannum(ii+nn)
                  enddo
                enddo
              enddo
            enddo
            if (allocated(rannum)) deallocate (rannum)
          endif
        endif
!
        if (me .eq. 0) then
          if (num_p3d .eq. 3) print *,' USING Ferrier-MICROPHYSICS'
          if (num_p3d .eq. 4) print *,' USING ZHAO-MICROPHYSICS'
        endif
        if (fhour .eq. 0.0) then
          do j=1,lats_node_r
            do i=1,lonr
!             phy_f2d(i,j,num_p2d) = 0.0
              phy_f2d(i,num_p2d,j) = 0.0
            enddo
          enddo
        endif
       
        if (ras) call ras_init(levs, me)
       
        first = .false.
      endif                  ! if (first) done
!
!     print *,' after if(first) before if semilag'

      if (semilag) then
!       dtphys = 300.0
        nsphys = max(int(tstep/dtphys),1)
        dtp    = tstep / nsphys
        dtf    = dtp
      else
!       dtphys = 3600.
        nsphys = max(int((tstep+tstep)/dtphys+0.9999),1)
        dtp    = (tstep+tstep)/nsphys
        dtf    = 0.5*dtp
      endif
      nnrcm  = max(1, nrcm/nsphys)
!
      if(lsfwd) dtf = dtp
!
      solhr = mod(phour+idate(1),cons_24)

! **************  Ken Campana Stuff  ********************************
!...  set switch for saving convective clouds
      if(lscca.and.lsswr) then
        clstp = 1100+min(fhswr,fhour,cons_99)  !initialize,accumulate,convert
      elseif(lscca) then
        clstp = 0100+min(fhswr,fhour,cons_99)  !accumulate,convert
      elseif(lsswr) then
        clstp = 1100                           !initialize,accumulate
      else
        clstp = 0100                           !accumulate
      endif
! **************  Ken Campana Stuff  ********************************
!
!


      if (.not. newsas) then  ! random number needed for RAS and old SAS
        if (random_clds) then
          iseed = mod(100.0*sqrt(fhour*3600),cons_1p0d9) + 1 + seed0
          call random_setseed(iseed)
          call random_number(wrk2)
          do nrc=1,nrcm
            indxr(nrc) = max(1, min(nint(wrk2(nrc)*maxran)+1,maxran))
          enddo
        endif
      endif
!
! doing ozone i/o and latitudinal interpolation to local gauss lats
!     ifozphys=.true.
 
      if (ntoz .gt. 0) then
       call ozinterpol(me,lats_node_r,lats_node_r,idate,fhour,
     &                 jindx1,jindx2,ozplin,ozplout,ddy)

!Moor   call ozinterpol(lats_node_r,lats_node_r,idate,fhour,
!    &                  jindx1,jindx2,ozplin,ozplout,ddy,
!    &                  global_lats_r,lonsperlar)
      endif

!     if (me == 0) write(0,*)' after ozinterpol'
!!
c ----------------------------------------------------
cc................................................................
cc
cc
      call f_hpmstart(41,"gb delnpe")
      call delnpe(trie_ls(1,1,p_zq   ),
     x            trio_ls(1,1,p_dphi),
     x            trie_ls(1,1,p_dlam),
     x            epse,epso,ls_node)
      call f_hpmstop(41)
cc
      call f_hpmstart(42,"gb delnpo")
      call delnpo(trio_ls(1,1,p_zq   ),
     x            trie_ls(1,1,p_dphi),
     x            trio_ls(1,1,p_dlam),
     x            epse,epso,ls_node)
      call f_hpmstop(42)
cc
cc
      call f_hpmstart(43,"gb dezouv dozeuv")
!$OMP parallel do shared(trie_ls,trio_ls)
!$OMP+shared(epsedn,epsodn,snnp1ev,snnp1od,ls_node)
!$OMP+private(k)
      do k=1,levs
         call dezouv(trie_ls(1,1,p_x  +k-1), trio_ls(1,1,p_w  +k-1),
     x               trie_ls(1,1,p_uln+k-1), trio_ls(1,1,p_vln+k-1),
     x               epsedn,epsodn,snnp1ev,snnp1od,ls_node)
cc
         call dozeuv(trio_ls(1,1,p_x  +k-1), trie_ls(1,1,p_w  +k-1),
     x               trio_ls(1,1,p_uln+k-1), trie_ls(1,1,p_vln+k-1),
     x               epsedn,epsodn,snnp1ev,snnp1od,ls_node)
      enddo
      call f_hpmstop(43)
cc
!     call mpi_barrier (mpi_comm_world,ierr)
cc
      call countperf(0,4,0.)
      call synctime()
      call countperf(1,4,0.)
!!
      dimg=0
      call countperf(0,1,0.)
cc
!     call f_hpmstart(48,"gb syn_ls2lats")
cc
!     call f_hpmstop(48)
cc
      call f_hpmstart(49,"gb sumfln")
cc
      call sumfln_slg_gg(trie_ls(1,1,p_q),
     x                   trio_ls(1,1,p_q),
     x                   lat1s_r,
     x                   plnev_r,plnod_r,
     x                   5*levs+3,ls_node,latr2,
     x                   lats_dim_r,lots,for_gr_r_1,
     x                   ls_nodes,max_ls_nodes,
     x                   lats_nodes_r,global_lats_r,
!mjr x                   lats_node_r,ipt_lats_node_r,lon_dims_r,
     x                   lats_node_r,ipt_lats_node_r,lon_dim_r,
     x                   lonsperlar,lonrx,latr,0)
!!
      if(.not.gg_tracers)then
        call sumfln_slg_gg(trie_ls(1,1,p_rt),
     x                     trio_ls(1,1,p_rt),
     x                     lat1s_r,
     x                     plnev_r,plnod_r,
     x                     levh,ls_node,latr2,
     x                     lats_dim_r,lots,for_gr_r_1,
     x                     ls_nodes,max_ls_nodes,
     x                     lats_nodes_r,global_lats_r,
!mjr x                     lats_node_r,ipt_lats_node_r,lon_dims_r,
     x                     lats_node_r,ipt_lats_node_r,lon_dim_r,
     x                     lonsperlar,lonrx,latr,5*levs+3)
      endif ! if(.not.gg_tracers)then
cc
      call f_hpmstop(49)
cc
      call countperf(1,1,0.)
cc
!     print *,' in GLOOPB after second sumfln'
cc
      pwatg = 0.
      ptotg = 0.

!--------------------
!     if( vertcoord_id == 3. ) then         ! For sigms/p/theta
!--------------------
        call countperf(0,11,0.)
        CALL countperf(0,1,0.)                                          ! hmhj
        call f_hpmstart(50,"gb sumder2")                                ! hmhj
!
        do lan=1,lats_node_r
          timer1=rtc()
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
     x               lat1s_r,                                           ! hmhj
     x               pddev_r,pddod_r,                                   ! hmhj
     x               levs,ls_node,latr2,                                ! hmhj
     x               lats_dim_r,lotd,                                   ! hmhj
     x               dyn_gr_r_1,                                        ! hmhj
     x               ls_nodes,max_ls_nodes,                             ! hmhj
     x               lats_nodes_r,global_lats_r,                        ! hmhj
!mjr x               lats_node_r,ipt_lats_node_r,lon_dims_r,            ! hmhj
     x               lats_node_r,ipt_lats_node_r,lon_dim_r,             ! hmhj
     x               lonsperlar,lonrx,latr,0)                           ! hmhj
!
        call f_hpmstop(50)                                              ! hmhj
        CALL countperf(1,1,0.)
! -----------------
!     endif
! -----------------
cc
      do lan=1,lats_node_r
        timer1=rtc()
!       if (me == 0) print *,' In lan loop  lan=', lan
cc

        lat = global_lats_r(ipt_lats_node_r-1+lan)
!       lon_dim = lon_dims_r(lan)
cc
        lons_lat = lonsperlar(lat)
!-----------------------------------------
!       if( vertcoord_id == 3. ) then
!-----------------------------------------

!!        calculate t rq u v zonal derivs. by multiplication with i*l
!!        note rlcs2=rcs2*L/rerth

          lmax = min(jcap,lons_lat/2)                                   ! hmhj
!
          ipt_ls=min(lat,latr-lat+1)                                    ! hmhj

          do i=1,lmax+1                                                 ! hmhj
            if ( ipt_ls .ge. lat1s_r(i-1) ) then                        ! hmhj
               reall=i-1                                                ! hmhj
               rlcs2(i)=reall*rcs2_r(ipt_ls)/rerth                      ! hmhj
            else                                                        ! hmhj
               rlcs2(i)=cons0     !constant                             ! hmhj
            endif                                                       ! hmhj
          enddo                                                         ! hmhj
!
!$omp parallel do private(k,i)
          do k=1,levs                                                   ! hmhj
            do i=1,lmax+1                                               ! hmhj
!
!           d(t)/d(lam)                                                 ! hmhj
               dyn_gr_r_1(2*i-1,(kdtlam-1+k),lan)=                      ! hmhj
     x        -for_gr_r_1(2*i  ,(kst   -1+k),lan)*rlcs2(i)              ! hmhj
               dyn_gr_r_1(2*i  ,(kdtlam-1+k),lan)=                      ! hmhj
     x         for_gr_r_1(2*i-1,(kst   -1+k),lan)*rlcs2(i)              ! hmhj
            enddo                                                       ! hmhj
          enddo                                                         ! hmhj
! -----------------------
!       endif
! -----------------------

!     print *,' in GLOOPB before four_to_grid'
cc
        call countperf(0,6,0.)
        call four_to_grid(for_gr_r_1(1,1,lan),for_gr_r_2(1,1,lan),
!mjr &                    lon_dim  ,lon_dim    ,lons_lat,5*levs+3)
     &                    lon_dim_r,lon_dim_r-2,lons_lat,5*levs+3)

        if(.not.gg_tracers)then
          CALL FOUR_TO_GRID(for_gr_r_1(1,ksr,lan),
     &                      for_gr_r_2(1,ksr,lan),
!mjr &                      lon_dim  ,lon_dim    ,lons_lat,levh)
     &                      lon_dim_r,lon_dim_r-2,lons_lat,levh)
        else  ! gg_tracers
          if (.not.shuff_lats_r) then
!                  set for_gr_r_2 to rg1_a rg2_a rg3_a from gloopa
            do n=1,ntrac
              do k=1,levs
                do i=1,min(lonf,lons_lat)
                  for_gr_r_2(i,ksr-1+k+(n-1)*levs,lan) =
     &               rgt_a(i,k,lats_node_a+1-lan,n)
                enddo
              enddo
            enddo
          endif              ! not shuff_lats_r
        endif                ! gg_tracers

!     print *,' in GLOOPB after four_to_grid '
! ----------------------------------
!       if( vertcoord_id == 3. ) then
! ----------------------------------
          CALL FOUR_TO_GRID(dyn_gr_r_1(1,kdtphi,lan),                ! hmhj
     &                      dyn_gr_r_2(1,kdtphi,lan),                ! hmhj
!mjr &                      lon_dim  ,lon_dim    ,lons_lat,levs)     ! hmhj
     &                      lon_dim_r,lon_dim_r-2,lons_lat,levs)     ! hmhj
          CALL FOUR_TO_GRID(dyn_gr_r_1(1,kdtlam,lan),                ! hmhj
     &                      dyn_gr_r_2(1,kdtlam,lan),                ! hmhj
!mjr &                      lon_dim  ,lon_dim    ,lons_lat,levs)     ! hmhj
     &                      lon_dim_r,lon_dim_r-2,lons_lat,levs)     ! hmhj
! ----------------------------
!       endif
! ---------------------------
        call countperf(1,6,0.)
!
        timer2 = rtc()
        global_times_b(lat,me+1) = timer2-timer1
c$$$    print*,'timeloopb',me,timer1,timer2,global_times_b(lat,me+1)
!!
      enddo   !lan
cc

      if(gg_tracers .and. shuff_lats_r) then
!        print*,' gloopb mpi_tracers_a_to_b shuff_lats_r',shuff_lats_r
         call mpi_tracers_a_to_b(
     x        rgt_a,lats_nodes_a,global_lats_a,
     x        for_gr_r_2(1,1,1),
     x        lats_nodes_r,global_lats_r,ksr,0)
      endif                       ! gg_tracers .and.  shuff_lats_r

      call f_hpmstart(51,"gb lat_loop2")

      do lan=1,lats_node_r

!
        lat = global_lats_r(ipt_lats_node_r-1+lan)
!!
!       lon_dim  = lon_dims_r(lan)
        lons_lat = lonsperlar(lat)
        pwatp    = 0.
        rcs2_lan = rcs2_r(min(lat,latr-lat+1))
        rcs_lan  = sqrt(rcs2_lan)

!$omp parallel do  schedule(dynamic,1) private(lon)
!$omp+private(sumq,xcp,hprime_v,swh_v,hlw_v,stc_v,smc_v,slc_v)
!$omp+private(nlons_v,sinlat_v,coslat_v,ozplout_v,rannum_v)
!$omp+private(prslk,prsl,prsik,prsi,phil,phii,dpshc,work1,tem)
!$omp+private(gu,gv1,gd,gq,gphi,glam,gt,gtv,gr,vvel,gtvx,gtvy)
!$omp+private(adt,adr,adu,adv,pgr,ugrd,vgrd,rqtk)
!!$omp+private(adt,adr,adu,adv,pgr,rcs_v,ugrd,vgrd,rqtk)
!$omp+private(phy_f3dv,phy_f2dv)
!$omp+private(dt3dt_v,dq3dt_v,du3dt_v,dv3dt_v,upd_mf_v,dwn_mf_v)
!$omp+private(det_mf_v,dkh_v,rnp_v)
!$omp+private(njeff,item,jtem,ktem,i,j,k,kss,n,nn,nnr)
!!!$omp+private(temlon,temlat,lprnt,ipt)
        do lon=1,lons_lat,ngptc
!!
          njeff = min(ngptc,lons_lat-lon+1)
!!
!     lprnt = .false.
!
!  --- ...  for debug test
!     alon = 236.25
!     alat = 56.189
!     alon = 26.25
!     alat = 6.66
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
!         exit
!       endif
!     enddo
!     lprnt = .false.
!!
          do k = 1, LEVS
            do j = 1, njeff
             jtem = lon-1+j
             gu (j,k)  = for_gr_r_2(jtem,KSU-1+k,lan)
             gv1(j,k)  = for_gr_r_2(jtem,KSV-1+k,lan)
             gd (j,k)  = for_gr_r_2(jtem,KSD-1+k,lan)
            enddo
          enddo
!
! p in cb by finite difference from henry juang not ln(p)               ! hmhj
          if(.not.gen_coord_hybrid) then                                ! hmhj
            do j=1,njeff
              item = lon+j-1
              work1 = exp(for_gr_r_2(item,ksq,lan))
              for_gr_r_2(item,ksq,lan) = work1
!             for_gr_r_2(item,ksq,lan) = exp(for_gr_r_2(item,ksq,lan))
            enddo
          endif     ! .not.gen_coord_hybrid                             ! hmhj
          do i=1,njeff
            item = lon+i-1
            gq(i)   = for_gr_r_2(item,ksq,lan)
            gphi(i) = for_gr_r_2(item,kspphi,lan)
            glam(i) = for_gr_r_2(item,ksplam,lan)
          enddo
!  Tracers
          do n=1,ntrac
            do k=1,levs
              item = KSR-1+k+(n-1)*levs
              do j=1,njeff
                gr(j,k,n) = for_gr_r_2(lon-1+j,item,lan)
              enddo
            enddo
          enddo

!
!   For omega in gen_coord_hybrid                                         ! hmhj
!   the same variables for thermodyn_id=3 for enthalpy                    ! hmhj
          if( gen_coord_hybrid ) then
            do k=1,levs                                                   ! hmhj
              do j=1,njeff                                                ! hmhj
                gtv(j,k) = for_gr_r_2(lon-1+j,kst-1+k,lan) 
              enddo
            enddo
! --------------------------------------
            if( vertcoord_id.eq.3. ) then
! --------------------------------------
              do k=1,levs                                                ! hmhj
                do j=1,NJEFF                                             ! hmhj
                  jtem = lon-1+j
                  gtvx(j,k) = dyn_gr_r_2(jtem,kdtlam-1+k,lan)
                  gtvy(j,k) = dyn_gr_r_2(jtem,kdtphi-1+k,lan)
                enddo                                                    ! hmhj
              enddo
! -----------------------------
            endif
! -----------------------------
            if( thermodyn_id.eq.3 ) then                                ! hmhj
!               get dry temperature from enthalpy                       ! hmhj
              do k=1,levs
                do j=1,njeff
                  sumq(j,k) = 0.0
                  xcp(j,k)  = 0.0
                enddo
              enddo
              do n=1,ntrac                                                ! hmhj
                if( cpi(n).ne.0.0 ) then                                  ! hmhj
!                 kss = ksr+(n-1)*levs                                    ! hmhj
                  do k=1,levs                                             ! hmhj
!                   ktem = kss+k-1
                    do j=1,njeff                                          ! hmhj
                      sumq(j,k) = sumq(j,k) + gr(j,k,n)                   ! hmhj
                      xcp(j,k)  = xcp(j,k)  + cpi(n)*gr(j,k,n)            ! hmhj
                    enddo                                                 ! hmhj
                  enddo                                                   ! hmhj
                endif                                                     ! hmhj
              enddo                                                       ! hmhj
              do k=1,levs                                                 ! hmhj
                do j=1,njeff                                              ! hmhj
                  work1 = (1.-sumq(j,k))*cpi(0) + xcp(j,k)                ! hmhj
                  gt(j,k) = gtv(j,k) / work1                              ! hmhj
                enddo                                                     ! hmhj
              enddo                                                       ! hmhj
! get dry temperature from virtual temperature                            ! hmhj
            else if( thermodyn_id.le.1 ) then                             ! hmhj
              do k=1,levs
                do j=1,njeff
                  gt(j,k) = gtv(j,k) / (1.0 + fv*max(gr(j,k,1),qmin))
                enddo
              enddo
            else
! get dry temperature from dry temperature                              ! hmhj
              do k=1,levs                                               ! hmhj
                do j=1,njeff                                            ! hmhj
                  gt(j,k) = gtv(j,k)                                    ! hmhj
                enddo                                                   ! hmhj
              enddo
            endif               ! if(thermodyn_id.eq.3)
          else
            do k=1,levs
              do j=1,njeff
                gt(j,k) = for_gr_r_2(lon+j-1,kst+k-1,lan)
     &                  / (1.0 + fv*max(gr(j,k,1),qmin))
              enddo
            enddo

          endif               ! if(gen_coord_hybrid)
!
          do j=1,njeff
            item = lon+j-1
            gq_save(item,lan) = for_gr_r_2(item,ksq,lan)
          enddo
!
!     if (lprnt) then
!     print *,' gq=',gq(ipt),' gphi=',gphi(ipt),glam(ipt)
!     print *,' gd=',gd(ipt,:)
!     print *,' gu=',gu(ipt,:)
!     print *,' gv1=',gv1(ipt,:)
!     endif
!                                   hmhj for gen_coord_hybrid
          if( gen_coord_hybrid ) then                                   ! hmhj

            call hyb2press_gc(njeff,ngptc,gq, gtv, prsi,prsl            ! hmhj
     &,                                            prsik, prslk)
!           call omegtes_gc(njeff,ngptc,rcs2_r(min(lat,latr-lat+1)),    ! hmhj
            call omegtes_gc(njeff,ngptc,rcs2_lan,                       ! hmhj
     &                   gq,gphi,glam,gtv,gtvx,gtvy,gd,gu,gv1,vvel)    ! hmhj
          elseif( hybrid )then                                           ! hmhj
!               vertical structure variables:   del,si,sl

!       if (lprnt) print *,' ipt=',ipt,' ugrb=',gu(ipt,levs),
!    &' vgrb=',gv1(ipt,levs),' lon=',lon
!    &,' xlon=',xlon(lon+ipt-1,lan),' xlat=',xlat(lon+ipt-1,lan)

            call  hyb2press(njeff,ngptc,gq, prsi, prsl,prsik,prslk)
!           call omegtes(njeff,ngptc,rcs2_r(min(lat,latr-lat+1)),
            call omegtes(njeff,ngptc,rcs2_lan,
     &                   gq,gphi,glam,gd,gu,gv1,vvel)
!    &                   gq,gphi,glam,gd,gu,gv1,vvel,lprnt,ipt)

!     if (lprnt) then
!     print *,' vvel=',vvel(ipt,:)
!     call mpi_quit(9999)
!     endif
          else                            ! for sigma coordinate
            call  sig2press(njeff,ngptc,gq,sl,si,slk,sik,
     &                                    prsi,prsl,prsik,prslk)
            call omegast3(njeff,ngptc,levs,
     &              gphi,glam,gu,gv1,gd,del,rcs2_lan,vvel,gq,sl)
!    &              gphi,glam,gu,gv1,gd,del,
!    &              rcs2_r(min(lat,latr-lat+1)),vvel,gq,sl)

          endif
!
          do i=1,njeff
            phil(i,levs) = 0.0 ! forces calculation of geopotential in gbphys
            pgr(i)       = gq(i) * 1000.0  ! Convert from kPa to Pa for physics
            prsi(i,1)    = pgr(i)
            dpshc(i)     = 0.3  * prsi(i,1)
!
            nlons_v(i)   = lons_lat
            sinlat_v(i)  = sinlat_r(lat)
            coslat_v(i)  = coslat_r(lat)
!           rcs_v(i)     = sqrt(rcs2_lan)
!           rcs_v(i)     = sqrt(rcs2_r(min(lat,latr-lat+1)))
          enddo
          do k=1,levs
            do i=1,njeff
              ugrd(i,k)   = gu(i,k)     * rcs_lan
              vgrd(i,k)   = gv1(i,k)    * rcs_lan
!             ugrd(i,k)   = gu(i,k)     * rcs_v(i)
!             vgrd(i,k)   = gv1(i,k)    * rcs_v(i)
              work1       = prsl(i,k)   * 1000.0
              prsl(i,k)   = work1
              work1       = prsi(i,k+1) * 1000.0
              prsi(i,k+1) = work1
              work1       = vvel(i,k)   * 1000.0  ! Convert from Cb/s to Pa/s
              vvel(i,k)   = work1
!             prsl(i,k)   = prsl(i,k)   * 1000.0
!             prsi(i,k+1) = prsi(i,k+1) * 1000.0
!             vvel(i,k)   = vvel(i,k)   * 1000.0  ! Convert from Cb/s to Pa/s
            enddo
          enddo
          if (gen_coord_hybrid .and. thermodyn_id == 3) then
            do i=1,ngptc
              prslk(i,1) = 0.0 ! forces calculation of geopotential in gbphys
              prsik(i,1) = 0.0 ! forces calculation of geopotential in gbphys
            enddo
          endif
!
          if (ntoz .gt. 0) then
            do j=1,pl_coeff
              do k=1,levozp
                do i=1,ngptc
!Moorthi          ozplout_v(i,k,j) = ozplout(lon+i-1,k,j,lan)
                  ozplout_v(i,k,j) = ozplout(k,lan,j)
                enddo
              enddo
            enddo
          endif
          do k=1,lsoil
            do i=1,njeff
              item = lon+i-1
              smc_v(i,k) = sfc_fld%smc(item,k,lan)
              stc_v(i,k) = sfc_fld%stc(item,k,lan)
              slc_v(i,k) = sfc_fld%slc(item,k,lan)
            enddo
          enddo
          do k=1,nmtvr
            do i=1,njeff
              hprime_v(i,k) = hprime(lon+i-1,k,lan)
            enddo
          enddo
!!
          do j=1,num_p3d
            do k=1,levs
              do i=1,njeff
                phy_f3dv(i,k,j) = phy_f3d(lon+i-1,k,j,lan)
              enddo
            enddo
          enddo
          do j=1,num_p2d
            do i=1,njeff
              phy_f2dv(i,j) = phy_f2d(lon+i-1,j,lan)
            enddo
          enddo
          if (.not. newsas) then
            if (random_clds) then
              do j=1,nrcm
                do i=1,njeff
                  rannum_v(i,j) = rannum_tank(lon+i-1,indxr(j),lan)
                enddo
!     if (me == 0) print *,' rannum=',rannum_v(1:6,j),' j=',j
              enddo
            else
              do j=1,nrcm
                do i=1,njeff
                  rannum_v(i,j) = 0.6    ! This is useful for debugging
                enddo
              enddo
            endif
          endif
          do k=1,levs
            do i=1,njeff
              item = lon+i-1
              swh_v(i,k) = swh(item,k,lan)
              hlw_v(i,k) = hlw(item,k,lan)
            enddo
          enddo
          if (ldiag3d) then
            do n=1,6
             do k=1,levs
              do i=1,njeff
               dt3dt_v(i,k,n) = dt3dt(lon+i-1,k,n,lan)
              enddo
             enddo
            enddo
            do n=1,4
             do k=1,levs
              do i=1,njeff
               du3dt_v(i,k,n) = du3dt(lon+i-1,k,n,lan)
               dv3dt_v(i,k,n) = dv3dt(lon+i-1,k,n,lan)
              enddo
             enddo
            enddo
          endif
          if (ldiag3d .or. lggfs3d) then
            do n=1,5+pl_coeff
             do k=1,levs
              do i=1,njeff
               dq3dt_v(i,k,n) = dq3dt(lon+i-1,k,n,lan)
              enddo
             enddo
            enddo
            do k=1,levs
             do i=1,njeff
              upd_mf_v(i,k) = upd_mf(lon+i-1,k,lan)
              dwn_mf_v(i,k) = dwn_mf(lon+i-1,k,lan)
              det_mf_v(i,k) = det_mf(lon+i-1,k,lan)
             enddo
            enddo
          endif
          if (lggfs3d) then
            do k=1,levs
             do i=1,njeff
              dkh_v(i,k) = dkh(lon+i-1,k,lan)
              rnp_v(i,k) = rnp(lon+i-1,k,lan)
             enddo
            enddo
          endif
!
!     write(0,*)' calling gbphys kdt=',kdt,' lon=',lon,' lan=',lan
!    &,' nlons_v=',ntoz,ntcw,nmtvr,lonr,latr,jcap,ras
!    &,' tisfc=',sfc_fld%tisfc(lon,lan)
!     print *,' temp=',for_gr_r_2(lon,kst,lan)
!
          do nn=1,nsphys
            nnr = (nn-1)*nnrcm + 1
            call gbphys                                                 &
!  ---  inputs:
     &    ( njeff,ngptc,levs,lsoil,lsm,ntrac,ncld,ntoz,ntcw,            &
     &      nmtvr,nnrcm,levozp,lonr,latr,jcap,num_p3d,num_p2d,          &
     &      kdt,lat,me,pl_coeff,nlons_v,ncw,flgmin,crtrh,cdmbgwd,       &
     &      ccwf,dlqf,ctei_rm,clstp,dtp,dtf,fhour,solhr,                &
     &      slag,sdec,cdec,sinlat_v,coslat_v,pgr,ugrd,vgrd,             &
     &      gt,gr,vvel,prsi,prsl,prslk,prsik,phii,phil,                 &
     &      rannum_v(1,nnr),ozplout_v,pl_pres,dpshc,                    &
     &      hprime_v, xlon(lon,lan),xlat(lon,lan),                      &
     &      sfc_fld%slope (lon,lan),    sfc_fld%shdmin(lon,lan),        &
     &      sfc_fld%shdmax(lon,lan),    sfc_fld%snoalb(lon,lan),        &
     &      sfc_fld%tg3   (lon,lan),    sfc_fld%slmsk (lon,lan),        &
     &      sfc_fld%vfrac (lon,lan),    sfc_fld%vtype (lon,lan),        &
     &      sfc_fld%stype (lon,lan),    sfc_fld%uustar(lon,lan),        &
     &      sfc_fld%oro   (lon,lan),    sfc_fld%oro_uf(lon,lan),        &   
     &      flx_fld%coszen(lon,lan),                                    &
     &      flx_fld%sfcdsw(lon,lan),    flx_fld%sfcnsw(lon,lan),        &
     &      flx_fld%sfcdlw(lon,lan),    flx_fld%tsflw (lon,lan),        &
     &      flx_fld%sfcemis(lon,lan),   sfalb(lon,lan),                 &
     &      swh_v,                      hlw_v,                          &
     &      ras,pre_rad,ldiag3d,lggfs3d,lssav,lssav_cc,                 &
     &      bkgd_vdif_m,bkgd_vdif_h,bkgd_vdif_s,psautco,prautco, evpco, &
     &      wminco,                                                     &
     &      flipv,old_monin,cnvgwd,shal_cnv,sashal,newsas,cal_pre,      &
     &      mom4ice,mstrat,trans_trac,nst_fcst,moist_adj,fscav,         &
     &      thermodyn_id, sfcpress_id, gen_coord_hybrid,                &
!  ---  input/outputs:
     &      sfc_fld%hice  (lon,lan),    sfc_fld%fice  (lon,lan),        &
     &      sfc_fld%tisfc (lon,lan),    sfc_fld%tsea  (lon,lan),        &
     &      sfc_fld%tprcp (lon,lan),    sfc_fld%cv    (lon,lan),        &
     &      sfc_fld%cvb   (lon,lan),    sfc_fld%cvt   (lon,lan),        &
     &      sfc_fld%srflag(lon,lan),    sfc_fld%snwdph(lon,lan),        &
     &      sfc_fld%sheleg(lon,lan),    sfc_fld%sncovr(lon,lan),        &
     &      sfc_fld%zorl  (lon,lan),    sfc_fld%canopy(lon,lan),        &
     &      sfc_fld%ffmm  (lon,lan),    sfc_fld%ffhh  (lon,lan),        &
     &      sfc_fld%f10m  (lon,lan),    flx_fld%srunoff(lon,lan),       &
     &      flx_fld%evbsa (lon,lan),    flx_fld%evcwa (lon,lan),        &
     &      flx_fld%snohfa(lon,lan),    flx_fld%transa(lon,lan),        &
     &      flx_fld%sbsnoa(lon,lan),    flx_fld%snowca(lon,lan),        &
     &      flx_fld%soilm (lon,lan),    flx_fld%tmpmin(lon,lan),        &
     &      flx_fld%tmpmax(lon,lan),    flx_fld%dusfc (lon,lan),        &
     &      flx_fld%dvsfc (lon,lan),    flx_fld%dtsfc (lon,lan),        &
     &      flx_fld%dqsfc (lon,lan),    flx_fld%geshem(lon,lan),        &
     &      flx_fld%gflux (lon,lan),    flx_fld%dlwsfc(lon,lan),        & 
     &      flx_fld%ulwsfc(lon,lan),    flx_fld%suntim(lon,lan),        &
     &      flx_fld%runoff(lon,lan),    flx_fld%ep    (lon,lan),        &
     &      flx_fld%cldwrk(lon,lan),    flx_fld%dugwd (lon,lan),        &
     &      flx_fld%dvgwd (lon,lan),    flx_fld%psmean(lon,lan),        &
     &      flx_fld%bengsh(lon,lan),    flx_fld%spfhmin(lon,lan),       &
     &      flx_fld%spfhmax(lon,lan),                                   &
     &      dt3dt_v, dq3dt_v, du3dt_v, dv3dt_v,                         &
     &      acv(lon,lan), acvb(lon,lan), acvt(lon,lan),                 &
     &      slc_v, smc_v, stc_v,                                        &
     &      upd_mf_v, dwn_mf_v, det_mf_v, dkh_v, rnp_v,                 &
     &      phy_f3dv, phy_f2dv,                                         &
     &      DLWSFC_cc(lon,lan),  ULWSFC_cc(lon,lan),                    &
     &      DTSFC_cc(lon,lan),   SWSFC_cc(lon,lan),                     &
     &      DUSFC_cc(lon,lan),   DVSFC_cc(lon,lan),                     &
     &      DQSFC_cc(lon,lan),   PRECR_cc(lon,lan),                     &

     &      nst_fld%xt(lon,lan),        nst_fld%xs(lon,lan),            &
     &      nst_fld%xu(lon,lan),        nst_fld%xv(lon,lan),            &
     &      nst_fld%xz(lon,lan),        nst_fld%zm(lon,lan),            &
     &      nst_fld%xtts(lon,lan),      nst_fld%xzts(lon,lan),          &
     &      nst_fld%d_conv(lon,lan),    nst_fld%ifd(lon,lan),           &
     &      nst_fld%dt_cool(lon,lan),   nst_fld%Qrain(lon,lan),         &
!  ---  outputs:
     &      adt, adr, adu, adv,                                         &
     &      sfc_fld%t2m   (lon,lan),    sfc_fld%q2m   (lon,lan),        &
     &      flx_fld%u10m  (lon,lan),    flx_fld%v10m  (lon,lan),        &
     &      flx_fld%zlvl  (lon,lan),    flx_fld%psurf (lon,lan),        &
     &      flx_fld%hpbl  (lon,lan),    flx_fld%pwat  (lon,lan),        &
     &      flx_fld%t1    (lon,lan),    flx_fld%q1    (lon,lan),        &
     &      flx_fld%u1    (lon,lan),    flx_fld%v1    (lon,lan),        &
     &      flx_fld%chh   (lon,lan),    flx_fld%cmm   (lon,lan),        &
     &      flx_fld%dlwsfci(lon,lan),   flx_fld%ulwsfci(lon,lan),       &
     &      flx_fld%dswsfci(lon,lan),   flx_fld%uswsfci(lon,lan),       &
     &      flx_fld%dtsfci(lon,lan),    flx_fld%dqsfci(lon,lan),        &
     &      flx_fld%gfluxi(lon,lan),    flx_fld%epi   (lon,lan),        &
     &      flx_fld%smcwlt2(lon,lan),   flx_fld%smcref2(lon,lan),       &
!hchuang code change [+3L] 11/12/2007 : add 2D
     &     flx_fld%gsoil(lon,lan),      flx_fld%gtmp2m(lon,lan),        &
     &     flx_fld%gustar(lon,lan),     flx_fld%gpblh(lon,lan),         &
     &     flx_fld%gu10m(lon,lan),      flx_fld%gv10m(lon,lan),         &
     &     flx_fld%gzorl(lon,lan),      flx_fld%goro(lon,lan),          &

     &      XMU_cc(lon,lan), DLW_cc(lon,lan), DSW_cc(lon,lan),          &
     &      SNW_cc(lon,lan), LPREC_cc(lon,lan),                         &

     &      nst_fld%Tref(lon,lan),       nst_fld%z_c(lon,lan),          &
     &      nst_fld%c_0 (lon,lan),       nst_fld%c_d(lon,lan),          &
     &      nst_fld%w_0 (lon,lan),       nst_fld%w_d(lon,lan),          &
     &      rqtk                                                        &! rqtkD
!    &      bak_gr_r_2(lon,kap,lan),                                    &! rqtkD
     &      )
!!
            if (nn < nsphys) then
              gt   = adt
              gr   = adr
              ugrd = adu
              vgrd = adv
            endif
          enddo                                ! end of nsphys loop
          do k=1,levs
            do i=1,njeff
              item = lon + i - 1
              bak_gr_r_2(item,kau+k-1,lan) = adu(i,k) * rcs_lan
              bak_gr_r_2(item,kav+k-1,lan) = adv(i,k) * rcs_lan
!             bak_gr_r_2(item,kau+k-1,lan) = adu(i,k) * rcs_v(i)
!             bak_gr_r_2(item,kav+k-1,lan) = adv(i,k) * rcs_v(i)
              bak_gr_r_2(item,kat+k-1,lan) = adt(i,k)
            enddo
          enddo
          do n=1,ntrac
            do k=1,levs
              ktem = kar+k-1+(n-1)*levs
              do i=1,njeff
                item = lon + i - 1
                bak_gr_r_2(item,ktem,lan) = adr(i,k,n)
              enddo
            enddo
          enddo
          if (gg_tracers) then
            do i=1,njeff
              bak_gr_r_2(lon+i-1,kap,lan) = rqtk(i)
            enddo
          else
            do i=1,njeff
              bak_gr_r_2(lon+i-1,kap,lan) = 0.0
            enddo
          endif
!!
!<-- cpl insertion: instantanious variables
          do i=1,njeff
            item = lon+i-1
            U_BOT_cc(item,lan)  = adu(i,1)
            V_BOT_cc(item,lan)  = adv(i,1)
            Q_BOT_cc(item,lan)  = adr(i,1,1)
            P_BOT_cc(item,lan)  = prsl(i,1)
            P_SURF_cc(item,lan) = prsi(i,1)
          enddo

          do i=1,njeff
            item = lon+i-1
            T_BOT_cc(item,lan) = adt(i,1)
            tem                = adt(i,1)*(1+RVRDM1*adr(i,1,1))
            Z_BOT_cc(item,lan) = -(RD/grav)*tem
     &                         * LOG(prsl(i,1)/prsi(i,1))
!
            ffmm_cc(item,lan)  = sfc_fld%ffmm(item,lan)
            ffhh_cc(item,lan)  = sfc_fld%ffhh(item,lan)
            if (sfc_fld%SLMSK(item,lan) .lt. 0.01) then
              T_SFC_cc(item,lan) = sfc_fld%tsea(item,lan)
     &                           + (sfc_fld%oro(item,lan)
     &                           -  sfc_fld%oro_uf(item,lan))*RLAPSE
            else
              T_SFC_cc(item,lan) = sfc_fld%tsea(item,lan)
            end if
            FICE_SFC_cc(item,lan) = sfc_fld%fice(item,lan)
            HICE_SFC_cc(item,lan) = sfc_fld%hice(item,lan)
     &                            * sfc_fld%fice(item,lan)
          enddo
!         do i=istrt,istrt+njeff-1
!          if (ffmm_cc(i,lan).LT.1.0) print *,'ffmm_cc<1',ffmm_cc(i,lan)
!          if (ffhh_cc(i,lan).LT.1.0) print *,'ffhh_cc<1',ffmm_cc(i,lan)
!         enddo
!         if (me .eq. 0) then
!           call atm_maxmin(njeff,1,LPREC_cc(lon,lan),
!     >     'in gbphys_call, LPREC_cc')
!           print *,'after cpl,istrt=',istrt,'istrt+njeff-1=',
!     >       istrt+njeff-1,'lan=',lan
!         endif
!--> cpl insertion

          if( gen_coord_hybrid .and. thermodyn_id.eq.3 ) then           ! hmhj

!            convert dry temperature to enthalpy                        ! hmhj
            do k=1,levs
              do j=1,njeff
                item = lon+j-1
                sumq(j,k) = 0.0
                xcp(j,k)  = 0.0
              enddo
            enddo
            do i=1,ntrac                                                ! hmhj
              kss = kar+(i-1)*levs
              if( cpi(i).ne.0.0 ) then                                  ! hmhj
                do k=1,levs                                             ! hmhj
                  ktem = kss+k-1
                  do j=1,njeff                                          ! hmhj
                    item      = lon+j-1
                    work1     = bak_gr_r_2(item,ktem,lan)               ! hmhj
                    sumq(j,k) = sumq(j,k) + work1                	! hmhj
                    xcp(j,k)  = xcp(j,k)  + cpi(i)*work1   	        ! hmhj
                enddo                                                   ! hmhj
               enddo                                                    ! hmhj
             endif                                                      ! hmhj
            enddo                                                       ! hmhj
            do k=1,levs                                                 ! hmhj
              ktem = kat+k-1
              do j=1,njeff                                              ! hmhj
                item  = lon+j-1
                work1 = (1.-sumq(j,k))*cpi(0) + xcp(j,k)                ! hmhj
                bak_gr_r_2(item,ktem,lan) = bak_gr_r_2(item,ktem,lan)
     &                                    * work1                       ! hmhj
                adt(j,k) = adt(j,k)*work1
              enddo                                                     ! hmhj
            enddo                                                       ! hmhj

          else                                                          ! hmhj

!           convert dry temperture to virtual temperature               ! hmhj
            do k=1,levs                                                 ! hmhj
              ktem = kar+k-1
              jtem = kat+k-1
              do j=1,njeff                                              ! hmhj
                item  = lon+j-1
                work1 = 1.0 + fv * max(bak_gr_r_2(item,ktem,lan),qmin)  ! hmhj
                bak_gr_r_2(item,jtem,lan) = bak_gr_r_2(item,jtem,lan)
     &                                    * work1                       ! hmhj
                adt(j,k) = adt(j,k)*work1
              enddo                                                     ! hmhj
            enddo                                                       ! hmhj

          endif
          if( gen_coord_hybrid .and. vertcoord_id == 3. ) then          ! hmhj
            prsi  = prsi * 0.001                 ! Convert from Pa to kPa
            if( thermodyn_id == 3. ) then                               ! hmhj
              call gbphys_adv_h(njeff,ngptc,dtf,gtv,gu,gv1,gr , gq,    ! hmhj
     &                              adt,adu,adv,adr,prsi )
!             call gbphys_adv_h(njeff,ngptc,dtf,
!    &                  for_gr_r_2(lon,kst,lan),
!    &                  for_gr_r_2(lon,ksu,lan),
!    &                  for_gr_r_2(lon,ksv,lan),
!    &                  for_gr_r_2(lon,ksr,lan),
!    &                  for_gr_r_2(lon,ksq,lan),
!    &                  bak_gr_r_2(lon,kat,lan),
!    &                  bak_gr_r_2(lon,kau,lan),
!    &                  bak_gr_r_2(lon,kav,lan),
!    &                  bak_gr_r_2(lon,kar,lan),
!    &                  prsi )
            else
              call gbphys_adv(njeff,ngptc,dtf,gtv,gu,gv1,gr,gq,       ! hmhj
     &                              adt,adu,adv,adr,prsi )
!             call gbphys_adv(njeff,ngptc,dtf,
!    &                  for_gr_r_2(lon,kst,lan),
!    &                  for_gr_r_2(lon,ksu,lan),
!    &                  for_gr_r_2(lon,ksv,lan),
!    &                  for_gr_r_2(lon,ksr,lan),
!    &                  for_gr_r_2(lon,ksq,lan),
!    &                  bak_gr_r_2(lon,kat,lan),
!    &                  bak_gr_r_2(lon,kau,lan),
!    &                  bak_gr_r_2(lon,kav,lan),
!    &                  bak_gr_r_2(lon,kar,lan),
!    &                  prsi )
              endif                                                       ! hmhj
            endif                                                         ! hmhj
!!
            do k=1,lsoil
              do i=1,njeff
                item = lon + i - 1
                sfc_fld%smc(item,k,lan) = smc_v(i,k)
                sfc_fld%stc(item,k,lan) = stc_v(i,k)
                sfc_fld%slc(item,k,lan) = slc_v(i,k)
              enddo
            enddo
!!
            do j=1,num_p3d
              do k=1,levs
                do i=1,njeff
                  phy_f3d(lon+i-1,k,j,lan) = phy_f3dv(i,k,j)
                enddo
              enddo
            enddo
            do j=1,num_p2d
              do i=1,njeff
                phy_f2d(lon+i-1,j,lan) = phy_f2dv(i,j)
              enddo
            enddo
!
          if (ldiag3d) then
            do n=1,6
             do k=1,levs
              do i=1,njeff
               dt3dt(lon+i-1,k,n,lan) = dt3dt_v(i,k,n)
              enddo
             enddo
            enddo
            do n=1,4
             do k=1,levs
              do i=1,njeff
               du3dt(lon+i-1,k,n,lan) = du3dt_v(i,k,n)
               dv3dt(lon+i-1,k,n,lan) = dv3dt_v(i,k,n)
              enddo
             enddo
            enddo
          endif
          if (ldiag3d .or. lggfs3d) then
            do n=1,5+pl_coeff
             do k=1,levs
              do i=1,njeff
               dq3dt(lon+i-1,k,n,lan) = dq3dt_v(i,k,n)
              enddo
             enddo
            enddo
            do k=1,levs
             do i=1,njeff
              upd_mf(lon+i-1,k,lan) = upd_mf_v(i,k)
              dwn_mf(lon+i-1,k,lan) = dwn_mf_v(i,k)
              det_mf(lon+i-1,k,lan) = det_mf_v(i,k)
             enddo
            enddo
          endif
          if (lggfs3d) then
            do k=1,levs
             do i=1,njeff
              dkh(lon+i-1,k,lan) = dkh_v(i,k)
              rnp(lon+i-1,k,lan) = rnp_v(i,k)
             enddo
            enddo
          endif
!
        enddo   ! lon loop
!
!
!       CALL dscal(LEVS*lonr,rcs2_v,bak_gr_r_2(1,kau,lan),1)
!       CALL dscal(LEVS*lonr,rcs2_v,bak_gr_r_2(1,kav,lan),1)
!
!
        ptotj(lan) = 0.
        do j=1,lons_lat
          ptotj(lan) = ptotj(lan) + gq_save(j,lan)
          pwatp      = pwatp + flx_fld%pwat(j,lan)
!     print *,' kdt=',kdt,' pwatp=',pwatp,' pwat=',flx_fld%pwat(j,lan)
!    &,' j=',j
        enddo
        pwatj(lan) = pwatp*grav/(2.*lonsperlar(lat)*1.e3)
        ptotj(lan) = ptotj(lan)/(2.*lonsperlar(lat))
!!
!!
c$$$        if (kdt.eq.1) then
c$$$        do j=1,lons_lat
c$$$        do i=1,levs
c$$$        write(8700+lat,*)
c$$$     &    bak_gr_r_2(j,kat-1+i,lan),i,j
c$$$        write(8800+lat,*)
c$$$     &    bak_gr_r_2(j,kar-1+i,lan),i,j
c$$$        write(8900+lat,*)
c$$$     &    bak_gr_r_2(j,kau-1+i,lan),i,j
c$$$        write(8100+lat,*)
c$$$     &    bak_gr_r_2(j,kav-1+i,lan),i,j
c$$$        write(8200+lat,*)
c$$$     &    bak_gr_r_2(j,kar-1+i+levs,lan),i,j
c$$$        write(8300+lat,*)
c$$$     &    bak_gr_r_2(j,kar-1+i+2*levs,lan),i,j
c$$$        enddo
c$$$        enddo
c$$$        endif
!!
      enddo   ! lan loop
!
      call f_hpmstop(51)
!
!     lotn=3*levs+1*levh
!
      do lan=1,lats_node_r    !    four_to_grid lan loop
!
         lat = global_lats_r(ipt_lats_node_r-1+lan)
!        lon_dim = lon_dims_r(lan)
         lons_lat = lonsperlar(lat)
!
         call countperf(0,6,0.)
!
         call grid_to_four(bak_gr_r_2(1,1,lan),bak_gr_r_1(1,1,lan),
     &                     lon_dim_r-2,lon_dim_r,lons_lat,3*levs+1)
!
         if (.not. gg_tracers .or. lsout) then
            call grid_to_four(bak_gr_r_2(1,kar,lan),
     &                        bak_gr_r_1(1,kar,lan),
     &                        lon_dim_r-2,lon_dim_r,lons_lat,levh)
         endif
         call countperf(1,6,0.)

        if (gg_tracers) then
          if (.not.shuff_lats_r) then
            item = lats_node_a + 1 - lan + yhalo
            do n=1,ntrac
              do k=1,levs
                jtem = levs + 1 - k
                ktem = kar  - 1 + k + (n-1)*levs
                do i=1,min(lonf,lons_lat)
                  rgt_h(xhalo+i,jtem,item,n) = bak_gr_r_2(i,ktem,lan)

c$$$            if (kdt .eq. 1) write(888,*) 'rg1_h, = ',
c$$$     &      i,k,lan, rg1_h(xhalo+i,levs+1-k,lats_node_a+1-lan+yhalo)
                enddo
              enddo
            enddo
          endif ! .not.shuff_lats_r
        endif ! gg_tracers
!
      enddo  !    fin four_to_grid lan loop
!
      if (gg_tracers .and. shuff_lats_r) then
!        print*,' gloopb mpi_tracers_b_to_a shuff_lats_r',shuff_lats_r
ccmr     call mpi_barrier (mc_comp,ierr)
         call mpi_tracers_b_to_a(
     &          bak_gr_r_2(1,1,1),
     &          lats_nodes_r,global_lats_r,
     &          rgt_h,lats_nodes_a,global_lats_a,kar,0)
       endif ! gg_tracers .and. shuff_lats_r

      call countperf(1,11,0.)
!!
      call countperf(0,4,0.)
      call synctime()
      call countperf(1,4,0.)
!!
      call excha(lats_nodes_r,global_lats_r,ptotj,pwatj,ptotg,pwatg)
      sumwa = 0.
      sumto = 0.
      do lat=1,latr
         sumto = sumto + wgt_r(min(lat,latr-lat+1))*ptotg(lat)
         sumwa = sumwa + wgt_r(min(lat,latr-lat+1))*pwatg(lat)
!     print *,' kdt=',kdt,' lat=',lat,' sumwa=',sumwa,' sumto=',sumto,
!    &' ptotg=',ptotg(lat),' pwatg=',pwatg(lat)
      enddo
cjfe
cjfe  write(70+me,*) sumto,sumwa,kdt
      pdryg = sumto - sumwa
!!
      if(pdryini == 0.) pdryini = pdryg

      if( gen_coord_hybrid ) then                               ! hmhj
        pcorr = (pdryini-pdryg)         * sqrt(2.)              ! hmhj
      else                                                      ! hmhj
        pcorr = (pdryini-pdryg) / sumto * sqrt(2.)
      endif                                                     ! hmhj
!!
!     call f_hpmstart(53,"gb lats2ls")
cc
cc
      call countperf(0,1,0.)
cc
!     call f_hpmstop(53)
!!
!     call f_hpmstart(54,"gb fl2eov")
!     call f_hpmstop(54)
!
      call f_hpmstart(52,"gb four2fln")
!
      call four2fln_gg(lats_dim_r,lota,3*levs+1,bak_gr_r_1,
     x              ls_nodes,max_ls_nodes,
!mjr x              lats_nodes_r,global_lats_r,lon_dims_r,
     x              lats_nodes_r,global_lats_r,lon_dim_r,
     x              lats_node_r,ipt_lats_node_r,
     x              lat1s_r,lonrx,latr,latr2,
     x              trie_ls(1,1,p_ze), trio_ls(1,1,p_ze),
     x              plnew_r, plnow_r,
     x              ls_node,0,
     x              2*levs+1,3*levs+1)

      sum_k_rqchange_ls(:,:) = trie_ls(:,:,p_q)
      sum_k_rqchango_ls(:,:) = trio_ls(:,:,p_q)

      trie_ls(:,:,p_q)       = save_qe_ls(:,:)
      trio_ls(:,:,p_q)       = save_qo_ls(:,:)
cc
      if (.not. gg_tracers .or.lsout ) then
         call four2fln_gg(lats_dim_r,lota,levh,bak_gr_r_1,
     x              ls_nodes,max_ls_nodes,
!mjr x              lats_nodes_r,global_lats_r,lon_dims_r,
     x              lats_nodes_r,global_lats_r,lon_dim_r,
     x              lats_node_r,ipt_lats_node_r,
     x              lat1s_r,lonrx,latr,latr2,
     x              trie_ls(1,1,p_rq), trio_ls(1,1,p_rq),
     x              plnew_r, plnow_r,
     x              ls_node,3*levs+1,
     x              1,levh)
      endif
!
      call f_hpmstop(52)
!
      call f_hpmstart(55,"gb uveodz uvoedz")
!
!$OMP parallel do shared(trie_ls,trio_ls)
!$OMP+shared(epse,epso,ls_node)
!$OMP+private(k)
      do k=1,levs
         call uveodz(trie_ls(1,1,p_ze +k-1), trio_ls(1,1,p_di +k-1),
     x               trie_ls(1,1,p_uln+k-1), trio_ls(1,1,p_vln+k-1),
     x               epse,epso,ls_node)
cc
         call uvoedz(trio_ls(1,1,p_ze +k-1), trie_ls(1,1,p_di +k-1),
     x               trio_ls(1,1,p_uln+k-1), trie_ls(1,1,p_vln+k-1),
     x               epse,epso,ls_node)
      enddo
      call f_hpmstop(55)
!
!.............................................................
      do k=1,levs
        ktem = p_w   + k - 1
        jtem = p_vln + k - 1
        do i=1,len_trie_ls
          trie_ls(i,1,ktem) = trie_ls(i,1,ktem) + bfilte(i)*
     &                       (trie_ls(i,1,jtem)-trie_ls(i,1,ktem))

          trie_ls(i,2,ktem) = trie_ls(i,2,ktem) + bfilte(i)*
     &                       (trie_ls(i,2,jtem)-trie_ls(i,2,ktem))
        enddo
        do i=1,len_trio_ls
          trio_ls(i,1,ktem) = trio_ls(i,1,ktem) + bfilto(i)*
     &                       (trio_ls(i,1,jtem)-trio_ls(i,1,ktem))

          trio_ls(i,2,ktem) = trio_ls(i,2,ktem) + bfilto(i)*
     &                       (trio_ls(i,2,jtem)-trio_ls(i,2,ktem))
        enddo
      enddo
cc.............................................................
      if(.not.gg_tracers)then
        do k=1,levs
          ktem = p_rt + k - 1
          jtem = p_rq + k - 1
          do i=1,len_trie_ls
            tem = bfilte(i)*(trie_ls(i,1,jtem)-trie_ls(i,1,ktem))
            trie_ls_rqt(i,1,k) = tem
            trie_ls(i,1,ktem)  = trie_ls(i,1,ktem) + tem
!
            tem = bfilte(i)*(trie_ls(i,2,jtem)-trie_ls(i,2,ktem))
            trie_ls_rqt(i,2,k) = tem
            trie_ls(i,2,ktem)  = trie_ls(i,2,ktem) + tem
          enddo
!!
          do i=1,len_trio_ls
            tem = bfilto(i)*(trio_ls(i,1,jtem)-trio_ls(i,1,ktem))
            trio_ls_rqt(i,1,k)   = tem
            trio_ls(i,1,ktem)     = trio_ls(i,1,ktem) + tem
!
            tem = bfilto(i)*(trio_ls(i,2,jtem)-trio_ls(i,2,ktem))
            trio_ls_rqt(i,2,k)    = tem
            trio_ls(i,2,p_rt+k-1) = trio_ls(i,2,ktem) + tem
          enddo
        enddo
!
!!.............................................................
!
        do nt=2,ntrac
          do k=levs*(nt-2)+1,levs*(nt-1)
            ktem = p_rt + levs + k - 1
            jtem = p_rq + levs + k - 1
            do i=1,len_trie_ls
              trie_ls(i,1,ktem) = trie_ls(i,1,ktem) + bfilte(i)*
     &                           (trie_ls(i,1,jtem)-trie_ls(i,1,ktem))
              trie_ls(i,2,ktem) = trie_ls(i,2,ktem) + bfilte(i)*
     &                           (trie_ls(i,2,jtem)-trie_ls(i,2,ktem))
            enddo
            do i=1,len_trio_ls
              trio_ls(i,1,ktem) = trio_ls(i,1,ktem) + bfilto(i)*
     &                           (trio_ls(i,1,jtem)-trio_ls(i,1,ktem))
              trio_ls(i,2,ktem) = trio_ls(i,2,ktem) + bfilto(i)*
     &                           (trio_ls(i,2,jtem)-trio_ls(i,2,ktem))
            enddo
          enddo
        enddo
      endif  ! if(.not.gg_tracers)
!!
!----------------------------------------------------------------------
!!
      if(hybrid)then
 
!     get some sigma distribution and compute   typdel from it.
 
      typical_pgr=85.
!sela  si(k)=(ak5(k)+bk5(k)*typical_pgr)/typical_pgr   !ak(k) bk(k) go top to botto
        do k=1,levp1
          si(levs+2-k) =  ak5(k)/typical_pgr + bk5(k)
        enddo
      endif

      DO k=1,LEVS
        typDEL(k)= SI(k)-SI(k+1)
      ENDDO
 
!----------------------------------------------------------------------
 
      if (ladj) then
        trie_ls(:,:,p_zq) = 0.
        trio_ls(:,:,p_zq) = 0.
        if (me == me_l_0) then
          trie_ls(1,1,p_zq) = pcorr
        endif
!!
        if( gen_coord_hybrid ) then                                       ! hmhj
          trie_ls_sfc = 0.0                                               ! hmhj
          trio_ls_sfc = 0.0                                               ! hmhj
          do k=1,levs                                                     ! hmhj
            do i=1,len_trie_ls                                            ! hmhj
              trie_ls_sfc(i,1) = trie_ls_sfc(i,1)
     &                         + typdel(k)*trie_ls_rqt(i,1,k)             ! hmhj
             trie_ls_sfc(i,2)  = trie_ls_sfc(i,2)
     &                         + typdel(k)*trie_ls_rqt(i,2,k)             ! hmhj
            enddo                                                         ! hmhj
            do i=1,len_trio_ls                                            ! hmhj
              trio_ls_sfc(i,1) = trio_ls_sfc(i,1)
     &                         + typdel(k)*trio_ls_rqt(i,1,k)             ! hmhj
             trio_ls_sfc(i,2)  = trio_ls_sfc(i,2)
     &                         + typdel(k)*trio_ls_rqt(i,2,k)             ! hmhj
            enddo                                                         ! hmhj
          enddo                                                           ! hmhj

          do i=1,len_trie_ls                                              ! hmhj
            trie_ls(i,1,p_zq) = trie_ls(i,1,p_zq)                         ! hmhj
     &                        + trie_ls(i,1,p_q )*trie_ls_sfc(i,1)        ! hmhj
            trie_ls(i,2,p_zq) = trie_ls(i,2,p_zq)                         ! hmhj
     &                        + trie_ls(i,2,p_q )*trie_ls_sfc(i,2)        ! hmhj
          enddo
          do i=1,len_trio_ls                                              ! hmhj
            trio_ls(i,1,p_zq) = trio_ls(i,1,p_zq)                         ! hmhj
     &                        + trio_ls(i,1,p_q )*trio_ls_sfc(i,1)        ! hmhj
            trio_ls(i,2,p_zq) = trio_ls(i,2,p_zq)                         ! hmhj
     &                        + trio_ls(i,2,p_q )*trio_ls_sfc(i,2)        ! hmhj
          enddo

        else                  ! For hybrid or sigma coordinate

          if(gg_tracers)then
            do i=1,len_trie_ls
              trie_ls(i,1,p_zq) = trie_ls(i,1,p_zq)
     &                          + sum_k_rqchange_ls(i,1)
              trie_ls(i,2,p_zq) = trie_ls(i,2,p_zq)
     &                          + sum_k_rqchange_ls(i,2)
            enddo
            do i=1,len_trio_ls
              trio_ls(i,1,p_zq) = trio_ls(i,1,p_zq)
     &                          + sum_k_rqchango_ls(i,1)
              trio_ls(i,2,p_zq) = trio_ls(i,2,p_zq)
     &                          + sum_k_rqchango_ls(i,2)
            enddo
          else
            do k=1,levs
              do i=1,len_trie_ls
                trie_ls(i,1,p_zq) = trie_ls(i,1,p_zq)
     &                            + typdel(k)*trie_ls_rqt(i,1,k)
                trie_ls(i,2,p_zq) = trie_ls(i,2,p_zq)
     &                            + typdel(k)*trie_ls_rqt(i,2,k)
              enddo
              do i=1,len_trio_ls
                trio_ls(i,1,p_zq) = trio_ls(i,1,p_zq)
     &                            + typdel(k)*trio_ls_rqt(i,1,k)
                trio_ls(i,2,p_zq) = trio_ls(i,2,p_zq)
     &                            + typdel(k)*trio_ls_rqt(i,2,k)
              enddo
            enddo
          endif               !fin if(gg_tracers)

        endif                 !fin  if (gen_coord_hybrid)                 ! hmhj
!!
        do k=1,levs
          item = p_di+k-1
          jtem = p_uln+k-1
          ktem = p_x+k-1
          ltem = p_te+k-1
          mtem = p_y+k-1

          do i=1,len_trie_ls
            trie_ls(i,1,item) = bfilte(i)
     &                        * (trie_ls(i,1,jtem)-trie_ls(i,1,ktem))
            trie_ls(i,1,ltem) = bfilte(i)
     &                        * (trie_ls(i,1,ltem)-trie_ls(i,1,mtem))
            trie_ls(i,2,item) = bfilte(i)
     &                        * (trie_ls(i,2,jtem)-trie_ls(i,2,ktem))
            trie_ls(i,2,ltem) = bfilte(i)
     &                        * (trie_ls(i,2,ltem)-trie_ls(i,2,mtem))
          enddo
          do i=1,len_trio_ls
            trio_ls(i,1,item) = bfilto(i)
     &                        * (trio_ls(i,1,jtem)-trio_ls(i,1,ktem))
            trio_ls(i,1,ltem) = bfilto(i)
     &                        * (trio_ls(i,1,ltem)-trio_ls(i,1,mtem))
            trio_ls(i,2,item) = bfilto(i)
     &                        * (trio_ls(i,2,jtem)-trio_ls(i,2,ktem))
            trio_ls(i,2,ltem) = bfilto(i)
     &                        * (trio_ls(i,2,ltem)-trio_ls(i,2,mtem))
          enddo
        enddo
 
!---------------------------------------------------------
        if( gen_coord_hybrid ) then                                       ! hmhj

!$OMP parallel do private(locl)
          do locl=1,ls_max_node                                           ! hmhj

            call impadje_hyb_gc(trie_ls(1,1,p_x),trie_ls(1,1,p_y),        ! hmhj
     &                    trie_ls(1,1,p_q),trie_ls(1,1,p_di),             ! hmhj
     &                    trie_ls(1,1,p_te),trie_ls(1,1,p_zq),            ! hmhj
     &                      tstep,                                        ! hmhj
     &                    trie_ls(1,1,p_uln),trie_ls(1,1,p_vln),          ! hmhj
     &                    snnp1ev,ndexev,ls_node,locl)                    ! hmhj
!!
            call impadjo_hyb_gc(trio_ls(1,1,p_x),trio_ls(1,1,p_y),        ! hmhj
     &                    trio_ls(1,1,p_q),trio_ls(1,1,p_di),             ! hmhj
     &                    trio_ls(1,1,p_te),trio_ls(1,1,p_zq),            ! hmhj
     &                      tstep,                                        ! hmhj
     &                    trio_ls(1,1,p_uln),trio_ls(1,1,p_vln),          ! hmhj
     &                    snnp1od,ndexod,ls_node,locl)                    ! hmhj
          enddo                                                           ! hmhj
        elseif(hybrid) then           ! for sigma/p hybrid coordinate    ! hmhj
          if (.not. semilag) then     ! for Eulerian hybrid case

 
!$OMP parallel do private(locl)
            do locl=1,ls_max_node
              call impadje_hyb(trie_ls(1,1,p_x),trie_ls(1,1,p_y),
     &             trie_ls(1,1,p_q),trie_ls(1,1,p_di),
     &             trie_ls(1,1,p_te),trie_ls(1,1,p_zq),
     &                      tstep,
     &             trie_ls(1,1,p_uln),trie_ls(1,1,p_vln),
     &             snnp1ev,ndexev,ls_node,locl)
!!
              call impadjo_hyb(trio_ls(1,1,p_x),trio_ls(1,1,p_y),
     &             trio_ls(1,1,p_q),trio_ls(1,1,p_di),
     &             trio_ls(1,1,p_te),trio_ls(1,1,p_zq),
     &                      tstep,
     &             trio_ls(1,1,p_uln),trio_ls(1,1,p_vln),
     &             snnp1od,ndexod,ls_node,locl)
            enddo
          else        ! for semi-Lagrangian hybrid case
!$OMP parallel do private(locl)
            do locl=1,ls_max_node


              call impadje_slg(trie_ls(1,1,p_x),trie_ls(1,1,p_y),
     &             trie_ls(1,1,p_q),trie_ls(1,1,p_di),
     &             trie_ls(1,1,p_te),trie_ls(1,1,p_zq),
     &                      tstep,
     &             trie_ls(1,1,p_uln),trie_ls(1,1,p_vln),
     &             snnp1ev,ndexev,ls_node,locl,batah)
!!
              call impadjo_slg(trio_ls(1,1,p_x),trio_ls(1,1,p_y),
     &             trio_ls(1,1,p_q),trio_ls(1,1,p_di),
     &             trio_ls(1,1,p_te),trio_ls(1,1,p_zq),
     &                      tstep,
     &             trio_ls(1,1,p_uln),trio_ls(1,1,p_vln),
     &             snnp1od,ndexod,ls_node,locl,batah)
            enddo
          endif
 
        else  !  massadj in sigma coordinate
 
          call countperf(0,9,0.)
!$OMP parallel do private(locl)
          do locl=1,ls_max_node
            call impadje(trie_ls(1,1,p_x),trie_ls(1,1,p_y),
     &             trie_ls(1,1,p_q),trie_ls(1,1,p_di),
     &             trie_ls(1,1,p_te),trie_ls(1,1,p_zq),
     &             am,bm,sv,tstep,
     &             trie_ls(1,1,p_uln),trie_ls(1,1,p_vln),
     &             snnp1ev,ndexev,ls_node,locl)
!!
            call impadjo(trio_ls(1,1,p_x),trio_ls(1,1,p_y),
     &             trio_ls(1,1,p_q),trio_ls(1,1,p_di),
     &             trio_ls(1,1,p_te),trio_ls(1,1,p_zq),
     &             am,bm,sv,tstep,
     &             trio_ls(1,1,p_uln),trio_ls(1,1,p_vln),
     &             snnp1od,ndexod,ls_node,locl)
          enddo
 
          call countperf(1,9,0.)
 
        endif  ! fin massadj in sigma
!---------------------------------------------------------
 
      else  ! fin massadj,    following is with no masadj
        DO k=1,LEVS
          del(k) = typDEL(k)                 ! sela 4.5.07
        ENDDO
        if (me == me_l_0) then
          trie_ls(1,1,p_q) = trie_ls(1,1,p_q) + pcorr
        endif
!
! testing mass correction on sep 25
!!
        if(gg_tracers)then
          do i=1,len_trie_ls
            trie_ls(i,1,p_q) = trie_ls(i,1,p_q) + sum_k_rqchange_ls(i,1)
            trie_ls(i,2,p_q) = trie_ls(i,2,p_q) + sum_k_rqchange_ls(i,2)
          enddo
          do i=1,len_trio_ls
            trio_ls(i,1,p_q) = trio_ls(i,1,p_q) + sum_k_rqchango_ls(i,1)
            trio_ls(i,2,p_q) = trio_ls(i,2,p_q) + sum_k_rqchango_ls(i,2)
          enddo
        else
          do k=1,levs
            do i=1,len_trie_ls
             trie_ls(i,1,p_q)=trie_ls(i,1,p_q)+del(k)*trie_ls_rqt(i,1,k)
             trie_ls(i,2,p_q)=trie_ls(i,2,p_q)+del(k)*trie_ls_rqt(i,2,k)
            enddo
            do i=1,len_trio_ls
             trio_ls(i,1,p_q)=trio_ls(i,1,p_q)+del(k)*trio_ls_rqt(i,1,k)
             trio_ls(i,2,p_q)=trio_ls(i,2,p_q)+del(k)*trio_ls_rqt(i,2,k)
            enddo
          enddo
        endif
!
! testing mass correction on sep 25
!
        do k=1,levs
          item = p_di+k-1
          jtem = p_uln+k-1
          ktem = p_x+k-1
          ltem = p_te+k-1
          mtem = p_y+k-1
          do i=1,len_trie_ls
            trie_ls(i,1,ktem) = trie_ls(i,1,ktem) + bfilte(i)
     &                        *(trie_ls(i,1,jtem)-trie_ls(i,1,ktem))
            trie_ls(i,2,ktem) = trie_ls(i,2,ktem) + bfilte(i)
     &                        *(trie_ls(i,2,jtem)-trie_ls(i,2,ktem))
            trie_ls(i,1,mtem) = trie_ls(i,1,mtem) + bfilte(i)
     &                        *(trie_ls(i,1,ltem)-trie_ls(i,1,mtem))
            trie_ls(i,2,mtem) = trie_ls(i,2,mtem) + bfilte(i)
     &                        *(trie_ls(i,2,ltem)-trie_ls(i,2,mtem))
          enddo

          do i=1,len_trio_ls
            trio_ls(i,1,ktem) = trio_ls(i,1,ktem)+bfilto(i)
     &                        *(trio_ls(i,1,jtem)-trio_ls(i,1,ktem))
            trio_ls(i,2,ktem) = trio_ls(i,2,ktem) + bfilto(i)
     &                        *(trio_ls(i,2,jtem)-trio_ls(i,2,ktem))
            trio_ls(i,1,mtem) = trio_ls(i,1,mtem) + bfilto(i)
     &                        *(trio_ls(i,1,ltem)-trio_ls(i,1,mtem))
            trio_ls(i,2,mtem) = trio_ls(i,2,mtem) + bfilto(i)
     &                        *(trio_ls(i,2,ltem)-trio_ls(i,2,mtem))
          enddo
        enddo
      endif   ! fin no ladj (i.e. no massadj)
!!
      return
      end
