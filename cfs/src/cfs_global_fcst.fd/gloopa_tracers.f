       subroutine gloopa_hyb_slg
     &    (deltim,trie_ls,trio_ls,gzie_ln,gzio_ln,
     &     ls_node,ls_nodes,max_ls_nodes,
     &     lats_nodes_a,global_lats_a,
     &     lonsperlat,
     &     epse,epso,epsedn,epsodn,
     &     snnp1ev,snnp1od,ndexev,ndexod,
     &     plnev_a,plnod_a,pddev_a,pddod_a,plnew_a,plnow_a,
     &     global_times_a,kdt,batah,lsout)
      use machine             , only : kind_evod
      use resol_def           , only : jcap,jcap1,latg,latg2,ntrac,
     &                                 levh,levp1,levs,lnt2,lonf,lota,
     &                                 p_di,p_dlam,p_dphi,p_gz,p_q,
     &                                 p_rq,p_rt,p_te,p_uln,p_vln,
     &                                 p_w,p_x,p_y,p_ze
      use layout1             , only : ipt_lats_node_a,lat1s_a,
     &                                 lats_dim_a,lats_node_a,
     &                                 len_trie_ls,len_trio_ls,
     &                                 lon_dim_a,ls_dim,
     &                                 ls_max_node,me,nodes
      use gg_def              , only : rcs2_a
!mjr  use vert_def
      use date_def            , only : spdmax
      use namelist_def        , only : gg_tracers,hybrid,ref_temp,ngptc
      use mpi_def             , only : kind_mpi,mc_comp,
     &                                 mpi_r_mpi,mpi_real8
      use physcons            , only : rerth => con_rerth,
     &                                    rd => con_rd,
     &                                  grav => con_g,
     &                                 omega => con_omega
      use layout_lag          , only : lat1s_h,lats_dim_h,lon_dim_h
      use layout_grid_tracers , only : rgt_a,rgt_h,xhalo,yhalo
      use coordinate_def      , only : ak5,bk5
      use      pmgrid         , only : platd,plond
      use      slgshr         , only : dlam,lam,nlonex,ra,rdlam,rdlam6
!
      implicit none
!     character*20 hpmlab
      integer iter_slg,latco,maxkdt_slg,kdt_slg,maps,num_threads
      integer ifirst  , i_1,i_2,lon85,lat85
      real    ,allocatable ::            tb_k(:)
      integer ,allocatable ::    lats_nodes_h(:)
      integer ,allocatable ::   global_lats_h(:)
      real    ,allocatable :: lammp_a(:,:,:),phimp_a(:,:,:),
     &                        sigmp_a(:,:,:)
      real batah,r_rt ! used in ln(ps) 'linarization'

      real(kind=kind_evod) gzie_ln(len_trie_ls,2)
      real(kind=kind_evod) gzio_ln(len_trio_ls,2)

      real(kind=kind_evod) trie_ls(len_trie_ls,2,11*levs+3*levh+6)
      real(kind=kind_evod) trio_ls(len_trio_ls,2,11*levs+3*levh+6)
      integer              ls_node(ls_dim,3)
      integer              ls_nodes(ls_dim,nodes)
      integer              max_ls_nodes(nodes)
      integer              lats_nodes_a(nodes)
      real*8 rtc ,timer1,timer2
      real(kind=kind_evod) global_times_a(latg,nodes)
      integer              global_lats_a(latg)
      integer                 lonsperlat(latg)

      real(kind=kind_evod)    epse(len_trie_ls)
      real(kind=kind_evod)    epso(len_trio_ls)
      real(kind=kind_evod)  epsedn(len_trie_ls)
      real(kind=kind_evod)  epsodn(len_trio_ls)
      real(kind=kind_evod) snnp1ev(len_trie_ls)
      real(kind=kind_evod) snnp1od(len_trio_ls)
      integer               ndexev(len_trie_ls)
      integer               ndexod(len_trio_ls)
      real(kind=kind_evod)   plnev_a(len_trie_ls,latg2)
      real(kind=kind_evod)   plnod_a(len_trio_ls,latg2)
      real(kind=kind_evod)   pddev_a(len_trie_ls,latg2)
      real(kind=kind_evod)   pddod_a(len_trio_ls,latg2)
      real(kind=kind_evod)   plnew_a(len_trie_ls,latg2)
      real(kind=kind_evod)   plnow_a(len_trio_ls,latg2)
!
      real(kind=kind_evod) ,allocatable ::
     &            for_coef_a(:,:,:), dyn_coef_a(:,:,:),
     &            bak_coef_a(:,:,:), for_grid_a(:,:,:),
     &            dyn_grid_a(:,:,:), bak_grid_a(:,:),
     &            etaint(:),         etamid(:),
     &            unlm(:,:,:),       unla(:,:,:),
     &            vnlm(:,:,:),       vnla(:,:,:),
     &            tnlm(:,:,:),       tnla(:,:,:),
     &            z_m(:,:,:),        sd_m(:,:),
     &            ud3_a(:,:,:),      vd3_a(:,:,:),
     &            td3_a(:,:,:),      dpdt_a(:,:),
     &            dpdt_d3_a(:,:,:),  ug_h(:,:,:),
     &            ug_m(:,:,:),       vg_h(:,:,:),
     &            vg_m(:,:,:),       ww_h(:,:,:),
     &            ud3_h1(:,:,:),     ud3_h2(:,:,:),
     &            vd3_h1(:,:,:),     vd3_h2(:,:,:),
     &            td3_h1(:,:,:),     td3_h2(:,:,:),
     &            dpdt_d3_h(:,:,:),
     &            phige(:,:),        phigo(:,:),
     &            gz_grid(:,:),      grad_gzlam(:,:),
     &            grad_gzphi(:,:),   rlcs2(:,:),
     &            omega_v(:),        cosf1(:),
     &            rcos1(:),
!    &            lam_all_dp(:,:,:), phi_all_dp(:,:,:),
!    &            lam_all_ar(:,:,:), phi_all_ar(:,:,:),
     &            rgt_hl(:,:,:)

      real(kind=kind_evod) trisca(lnt2)

      integer     i,ierr,iter,j,k,kap,kar,kat,kau,kav,kdrlam,
     &            kdrphi,kdtlam,kdtphi,kdulam,kduphi,kdvlam,
     &            kdvphi,ksd,ksplam,kspphi,ksq,ksrc,ksrg,kst,
     &            ksu,ksv,ksz,l,lan,lat,lmax,locl,
     &            lons_lat,n,node,nvcn,ii,
     &            init,s_lat_in_pe,n_lat_in_pe,
     &            ibmsign,ipt_ls,n_ggp,s_ggp,iter_max,n_rgt_a
     &,           ntr,s1,s2,s3
!
      real(kind=kind_evod) deltim, scale_ibm, xvcn, reall
      real(kind=kind_evod), allocatable :: spdlat(:,:)
      real(kind=kind_evod ) spdmax_node(levs),spdmax_nodes(levs,nodes)
      real(kind=kind_mpi)  spdmax_nodem(levs),spdmax_nodesm(levs,nodes)
!
      real ak(levs),bk(levs),pkref(levs),tkref(levs)
!Moo  real(kind=kind_evod), parameter ::  gama=6.5/1000.,
!Moo &                                    expon=rd*gama/grav,
!Moo &                                    psref=101.325,  tzero=288.
      real(kind=kind_evod), parameter ::  gama=0.0/1000.,
     &                                    expon=rd*gama/grav,
     &                                    psref=101.325,  zero=0.0d0
!    &                                    psref=101.325,  tzero=300.
      real tref, tzero,rcs2_loc
!sela real(kind=kind_evod) rlcs2(jcap1,lats_node_a)
!
      integer kdt,kdt_save,indlsev,jbasev , indlsod,jbasod
      integer ini_slg,ini_dp
      integer njeff,lon,iblk,ngptcd,ngptcs,xlon,ylan
!
      include 'function2'
      real(kind=kind_evod) consv

!     real(kind=kind_evod) xxsi(levp1), xxsl(levs), delsig
!     real(kind=kind_evod) gg_ncar(lonf,latg)
!
      logical lsout,go_forward,lprnt
!
      save  gz_grid,phige,phigo,num_threads,r_rt,tb_k,omega_v
      save  grad_gzlam,grad_gzphi,lats_nodes_h,global_lats_h
      save  lammp_a,phimp_a,sigmp_a,z_m,sd_m,unlm,vnlm,tnlm
      save  etamid,etaint,cosf1,rcos1,ug_m,vg_m
      save  s_lat_in_pe,n_lat_in_pe,iter_max,rlcs2
      save  ifirst,kdt_save,ini_slg,ini_dp

!     character*10 label_1
      data ifirst /1/,kdt_save/0/
      real zdp(lonf,levs,latg)
      real sc,lat_phi,degrad,hfdt
      integer jtem, ktem, kk, nn, lanx
      integer num_parthds
!
      logical, save :: first 
      data first/.true./

!       print*,' begin gloopa_tracers kdt=',kdt,' ifirst=',ifirst

      xhalo = 1
      yhalo = 6
!mjr  lon_dim_h  = lonf + 1 + 2*xhalo  !  even or odd

      allocate ( for_coef_a(lon_dim_a,(4*levs+levh+3),lats_dim_a) )
      allocate ( for_grid_a(lonf     ,(2*levs+levh+3),lats_dim_a) )
      allocate ( dyn_coef_a(lon_dim_a,2*levs,lats_dim_a) )
      allocate ( dyn_grid_a(lonf     ,2*levs,lats_dim_a) )
      allocate ( bak_coef_a(lon_dim_a,lota,lats_dim_a) )
      allocate ( bak_grid_a(lonf ,lota) )
      allocate ( unla(lonf,levs,lats_dim_a) )
      allocate ( vnla(lonf,levs,lats_dim_a) )
      allocate ( tnla(lonf,levs,lats_dim_a) )
      allocate ( dpdt_a(lon_dim_a,lats_dim_a) )
      allocate ( ud3_a(lonf,levs,lats_dim_a) )
      allocate ( vd3_a(lonf,levs,lats_dim_a) )
      allocate ( td3_a(lonf,levs,lats_dim_a) )
!!
      if (.not. allocated(rgt_a))
     &          allocate (rgt_a(lonf,levs,lats_dim_a,ntrac))
!!
      allocate (  dpdt_d3_a(lonf,levs,lats_dim_a) )

!     allocate ( lam_all_dp(lonf,levs,lats_node_a))
!     allocate ( phi_all_dp(lonf,levs,lats_node_a))
!     allocate ( lam_all_ar(lonf,levs,lats_node_a))
!     allocate ( phi_all_ar(lonf,levs,lats_node_a))
!
!     print *,' IN GLOOPA lats_node_a=',lats_node_a,' lats_dim_a=',
!    &lats_dim_a,' me=',me

!     ngptcs = 25
!     ngptcd = 30

      lprnt = .false.
      ngptcs = ngptc
      ngptcd = ngptc

      if (me == 0 .and. kdt <10 ) print*,'ngptcd = ',ngptcd,
     &                                  ' ngptcs = ',ngptcs

      ksz     =                   1
      ksd     =                   1
      kst     =   levs          + 1
      ksq     = 2*levs          + 1
      ksplam  = 2*levs          + 2
      kspphi  = 2*levs          + 3
      ksu     = 2*levs          + 4
      ksv     = 3*levs          + 4
      ksrc    = 4*levs          + 4
      ksrg    = 2*levs          + 4
      kdtphi  =                   1
      kdrphi  =   levs          + 1
      kdtlam  =   levs          + 1
      kdrlam  = 2*levs +   levh + 1
      kdulam  = 2*levs + 2*levh + 1
      kdvlam  = 3*levs + 2*levh + 1
      kduphi  = 4*levs + 2*levh + 1
      kdvphi  = 5*levs + 2*levh + 1
      kau     =                   1
      kav     =   levs          + 1
      kat     = 2*levs          + 1
      kap     = 3*levs          + 1
      kar     = 3*levs          + 2

      if (ifirst == 1) then                      ! first ifirst
        rgt_a(:,:,:,1) = 0.
        allocate (tb_k(levs))
        tref   = ref_temp
        tzero  = ref_temp
        r_rt   = 1./(rd*tref)
!!!!    r_rt   = 0.0   ! for no ln(ps) linearization
!--------------------------------------------------
        degrad = 45./atan(1.)

        if (me == 0 .and. kdt < 10) print*,'in gloopa  expon=',expon,
     &' grav=',grav,' rd=',rd,' psref=',psref,' tzero=',tzero,
     &' r_rt=',r_rt

        do k=1,levs
          ak(k)    = 0.5*(ak5(k)+ak5(k+1))
          bk(k)    = 0.5*(bk5(k)+bk5(k+1))
          pkref(k) = ak(k) + bk(k)*psref
          tkref(k) = tzero*(pkref(k)/psref)**expon

          tb_k(k)  = bk(k)*psref*expon*tkref(k)/pkref(k) !order is top2bot
          tb_k(k)  = - tb_k(k)/(rd*tzero)                !check tzero or tref
!!!!      tb_k(k)  = 0.0 !  "for no temp. perturbation"
        enddo
        if(me == 0)then
          do k=1,levs
            print 111, k, pkref(k),tkref(k),tb_k(k),ak(k),bk(k)
          enddo
        endif

111     format('in gloopa k=',i3,2x,' pkref=',e10.3,2x,'tkref=',e10.3,
     . ' tb_k=',e12.4,2x,'ak=',e10.3,2x,'bk=',e10.3)
                                                                                               
        num_threads = num_parthds()
        allocate (omega_v(lats_dim_a) )
        allocate (  cosf1(lats_dim_a) )
        allocate (  rcos1(lats_dim_a) )
        allocate ( lammp_a(lonf,levs,lats_dim_a) )
        allocate ( phimp_a(lonf,levs,lats_dim_a) )
        allocate ( sigmp_a(lonf,levs,lats_dim_a) )
        allocate (    unlm(lonf,levs,lats_dim_a) )
        allocate (    vnlm(lonf,levs,lats_dim_a) )
        allocate (    tnlm(lonf,levs,lats_dim_a) )
        allocate (     z_m(lonf,levs,lats_dim_a) )
        allocate (    sd_m(lon_dim_a,lats_dim_a) )
109     format(1h ,'------------------------------------')
        print 109
        allocate (  etamid(levs)  )
        allocate (  etaint(levp1) )

        call set_pmgrid_h(lonf,latg,levs)

        allocate(dlam(platd),rdlam(platd),rdlam6(platd))
        allocate(lam(plond,platd),nlonex(platd))        

        ra = 1./rerth

        print*,' deltim in gloopa_slg=',deltim,' batah=',batah
        print 109

        allocate ( lats_nodes_h(nodes) )
        allocate ( lat1s_h(0:jcap) )
        allocate ( global_lats_h(latg+2*yhalo*nodes) )

        call getcon_lag(lats_nodes_a,global_lats_a,
     x                  lats_nodes_h,global_lats_h,
     x                  lonsperlat,xhalo,yhalo)

        print*,' xhalo=',xhalo,' yhalo=',yhalo,' lon_dim_h=',lon_dim_h

        allocate (    phige(len_trie_ls,2) )
        allocate (    phigo(len_trio_ls,2) )

        allocate ( rlcs2(jcap1,lats_node_a) )
        allocate (    gz_grid(lon_dim_a,lats_dim_a) )
        allocate ( grad_gzlam(lon_dim_a,lats_dim_a) )
        allocate ( grad_gzphi(lon_dim_a,lats_dim_a) )

        allocate (  ug_m(lon_dim_h,levs  ,lats_dim_h) )
        allocate (  vg_m(lon_dim_h,levs  ,lats_dim_h) )


        s_lat_in_pe = latg + 1
     &              - global_lats_a(ipt_lats_node_a-1+lats_node_a)
        s_lat_in_pe = s_lat_in_pe - yhalo + 2

!       print*,'me1z=',me,' s_lat_in_pe=',s_lat_in_pe
!    &,' lats_dim_h=',lats_dim_h

!       n_lat_in_pe = latg + 1 - global_lats_a(ipt_lats_node_a)
!       n_lat_in_pe = n_lat_in_pe + yhalo + 2

        n_lat_in_pe = s_lat_in_pe + lats_dim_h - 1

!       print*,'me1z=',me,' n_lat_in_pe=',n_lat_in_pe
!    &,' lats_dim_h=',lats_dim_h

!---------------------------------------------------
        do k=1,levp1
          etaint(k) = ak5(k)/psref + bk5(k)
          if (me == 0) print*,'etaint_def  etaint(k) k=',etaint(k),k
        enddo
        do k=1,levs
          etamid(k) = (etaint(k) + etaint(k+1) ) * 0.5
          if (me == 0) print*,'etamid_def  etamid(k) k=',etamid(k),k
        enddo
        if (me == 0) then
          do k=2,levs             ! note k=1 not used
            print*,'detamd_def k=',k, etamid(k)-etamid(k-1)
          enddo
        endif
!---------------------------------------------------
        iter_max = 0
        do lan=1,lats_node_a   ! begin  lan loop 0
          lat      = global_lats_a(ipt_lats_node_a-1+lan)
          lons_lat = lonsperlat(lat)
          iter_max = max ( iter_max , (lons_lat+ngptcd-1)/ngptcd )
          ipt_ls   = min(lat,latg-lat+1)
!--------------------------------------------
          lmax     = min(jcap,lons_lat/2)
          do i=1,lmax+1
            if ( ipt_ls >= lat1s_a(i-1) ) then
               reall        = i - 1
               rlcs2(i,lan) = reall*rcs2_a(ipt_ls)*ra
            else
               rlcs2(i,lan) = zero
            endif
!sela    print*,'rlcs2(i,lan)=',rlcs2(i,lan),' i=',i,' lan=',lan
          enddo
!--------------------------------------------
          cosf1(lan) = 1. / sqrt( rcs2_a(ipt_ls) )
          lat_phi    = degrad*acos(cosf1(lan))

!sela     print 297,lan,lat,ipt_ls,cosf1(lan),lat_phi
297   format(' lan=',i3,2x,' lat=',i3,' ipt_ls=',i4,2x,
     &  ' cosf1(lan)=',e12.4,' lat_phi=',f8.4)

          rcos1(lan)   = sqrt( rcs2_a(min(lat,latg-lat+1)) )
          omega_v(lan) = 2.*omega*rerth*cosf1(lan)
          consv        = 2.*omega*rerth
          if (me == 0) print*,' consv=2. omega rerth ',consv
!sela     omega_v(lan) = 0.
          lat_phi      = degrad*acos(1./rcos1(lan))

!sela     print*,'lan=',lan,'  omega_v=',omega_v(lan),
!sela.                     '  lat_phi=',lat_phi,
!sela.                     '  cosf1=',cosf1(lan),
!sela.                     '  rcos1=',rcos1(lan)

        enddo  !fin lan loop 0
c***************************************************
        call get_topo_grid_grad
     x          (cosf1,
     x          grad_gzlam,grad_gzphi,gzie_ln,gzio_ln,
     x          gz_grid,
     x          phige,phigo,
     x          for_coef_a,for_grid_a,
     x          trie_ls,trio_ls,
     x          ls_node,ls_nodes,max_ls_nodes,
     x          lats_nodes_a,global_lats_a,
     x          lonsperlat,
     x          epse,epso,
     x          plnev_a,plnod_a,
     x          0,r_rt) ! 1 for write ; 0 for no write
c***************************************************
        print*,' lats_dim_h = ',lats_dim_h
        ifirst = 0
      endif !  first ifirst

      allocate (   ug_h(lon_dim_h,levs  ,lats_dim_h) )
      allocate (   vg_h(lon_dim_h,levs  ,lats_dim_h) )
      allocate (   ww_h(lon_dim_h,levs  ,lats_dim_h) )
      allocate ( ud3_h1(lon_dim_h,levs  ,lats_dim_h) )
      allocate ( ud3_h2(lon_dim_h,levs  ,lats_dim_h) )
      allocate ( vd3_h1(lon_dim_h,levs  ,lats_dim_h) )
      allocate ( vd3_h2(lon_dim_h,levs  ,lats_dim_h) )
      allocate ( td3_h1(lon_dim_h,levs  ,lats_dim_h) )
      allocate ( td3_h2(lon_dim_h,levs  ,lats_dim_h) )
      allocate ( dpdt_d3_h(lon_dim_h,levs,lats_dim_h) )
!!

      if (.not.allocated(rgt_h))
     &         allocate (rgt_h(lon_dim_h,levs,lats_dim_h,ntrac))
!!
!for halo testing             ud3_h =0.     
!for halo testing             vd3_h =0.     
!for halo testing             td3_h =0.     
!for halo testing             dpdt_d3_h =0.     

!        if(kdt.lt.3)print*,' me =',me,'  iter_max=',iter_max

      allocate ( spdlat(levs,iter_max ) )
      call delnpe(trie_ls(1,1,p_q   ),
     x            trio_ls(1,1,p_dphi),
     x            trie_ls(1,1,p_dlam),
     x            epse,epso,ls_node)
      call delnpo(trio_ls(1,1,p_q   ),
     x            trie_ls(1,1,p_dphi),
     x            trio_ls(1,1,p_dlam),
     x            epse,epso,ls_node)
!$omp parallel do shared(trie_ls,trio_ls)
!$omp+shared(epsedn,epsodn,snnp1ev,snnp1od,ls_node)
!$omp+private(k)
      do k=1,levs
         call dezouv(trie_ls(1,1,p_di +k-1), trio_ls(1,1,p_ze +k-1),
     x               trie_ls(1,1,p_uln+k-1), trio_ls(1,1,p_vln+k-1),
     x               epsedn,epsodn,snnp1ev,snnp1od,ls_node)
         call dozeuv(trio_ls(1,1,p_di +k-1), trie_ls(1,1,p_ze +k-1),
     x               trio_ls(1,1,p_uln+k-1), trie_ls(1,1,p_vln+k-1),
     x               epsedn,epsodn,snnp1ev,snnp1od,ls_node)
      enddo
      call synctime()

      call sumfln_slg_gg(trie_ls(1,1,p_di),
     x                   trio_ls(1,1,p_di),
     x            lat1s_a,
     x            plnev_a,plnod_a,
     x            4*levs+3,ls_node,latg2,
     x            lats_dim_a,4*levs+levh+3,
     x            for_coef_a,
     x            ls_nodes,max_ls_nodes,
     x            lats_nodes_a,global_lats_a,
     x            lats_node_a,ipt_lats_node_a,lon_dim_a,
     x            lonsperlat,lon_dim_a,latg,0)

      if (.not. gg_tracers .or. kdt == 1 ) then
           call sumfln_slg_gg(trie_ls(1,1,p_rq),
     x                        trio_ls(1,1,p_rq),
     x            lat1s_a,
     x            plnev_a,plnod_a,
     x            levh,ls_node,latg2,
     x            lats_dim_a,4*levs+levh+3,
     x            for_coef_a,
     x            ls_nodes,max_ls_nodes,
     x            lats_nodes_a,global_lats_a,
     x            lats_node_a,ipt_lats_node_a,lon_dim_a,
     x            lonsperlat,lon_dim_a,latg,4*levs+3)
      endif

      call sumder2_slg(trie_ls(1,1,p_te),
     x                 trio_ls(1,1,p_te),
     x             lat1s_a,
     x             pddev_a,pddod_a,
     x             levs,ls_node,latg2,
     x             lats_dim_a,2*levs,
     x             dyn_coef_a,
     x             ls_nodes,max_ls_nodes,
     x             lats_nodes_a,global_lats_a,
     x             lats_node_a,ipt_lats_node_a,lon_dim_a,
     x             lonsperlat,lon_dim_a,latg,0)

      do k=1,levs
        spdmax_node(k) = zero
      enddo

      do lan=1,lats_node_a   !begin lan loop 1 (four and grid)
       timer1   = rtc()
       lat      = global_lats_a(ipt_lats_node_a-1+lan)
       lons_lat = lonsperlat(lat)
       lmax     = min(jcap,lons_lat/2)
!$omp parallel do private(k,i)
         do k=1,levs
            do i=1,lmax+1
               dyn_coef_a(2*i-1,kdtlam-1+k,lan)=
     x        -for_coef_a(2*i  ,kst   -1+k,lan)*rlcs2(i,lan)
               dyn_coef_a(2*i  ,kdtlam-1+k,lan)=
     x         for_coef_a(2*i-1,kst   -1+k,lan)*rlcs2(i,lan)
            enddo
            if (   (lmax+1)*2+1 .le. lons_lat+2 ) then
            do i = (lmax+1)*2+1,     lons_lat+2
             dyn_coef_a(i,kdtlam-1+k,lan)=zero
!            dyn_coef_a(i,kdtphi-1+k,lan)=zero
            enddo
            endif
         end do
!!!!  synthesize  grids in South to north order
         call four_to_grid(for_coef_a(1,ksd,lan),
     &                     for_grid_a(1,ksd,lats_node_a+1-lan),
     &                     lon_dim_a,lon_dim_a-2,lons_lat,2*levs+3)
         call four_to_grid(for_coef_a(1,ksu,lan),
     &                     ug_h(1+xhalo,1,lats_node_a+1-lan+yhalo),
     &                     lon_dim_a,
     &                     lon_dim_h,
     &                     lons_lat,
     &                     levs)
         call four_to_grid(for_coef_a(1,ksv,lan),
     &                     vg_h(1+xhalo,1,lats_node_a+1-lan+yhalo),
     &                     lon_dim_a,
     &                     lon_dim_h,
     &                     lons_lat,
     &                     levs)
         if ( .not. gg_tracers .or. kdt == 1) then
            call four_to_grid(for_coef_a(1,ksrc,lan),
     &                        for_grid_a(1,ksrg,lats_node_a+1-lan),
     &                        lon_dim_a,lon_dim_a-2,lons_lat,levh)
         endif
         call four_to_grid(dyn_coef_a(1,1,lan),
     &                     dyn_grid_a(1,1,lats_node_a+1-lan),
     &                     lon_dim_a,lon_dim_a-2,lons_lat,2*levs)
          timer2=rtc()
          global_times_a(lat,me+1)=timer2-timer1
      enddo !sela fin lan loop 1

      if (kdt == 1) then ! calculate max meridional velocity
        do lan=1,lats_node_a   !sela begin lan loop 2
          timer1   = rtc()
          lat      = global_lats_a(ipt_lats_node_a-1+lan)
          lons_lat = lonsperlat(lat)

!sela  ngptcd=(lons_lat+num_threads-1)/num_threads ! opt ngptc
!$omp parallel do schedule(dynamic,1) private(lon)
!$omp+private(njeff,iblk)
!$omp+private(nvcn,xvcn)
          do lon=1,lons_lat,ngptcd
            njeff = min(ngptcd,lons_lat-lon+1)
            iblk  = (lon-1)/ngptcd+ 1
            call gfidi_speed
     &        (lon_dim_a, lon_dim_h, njeff, lat,
     &        vg_h(xhalo+lon,1,lats_node_a+1-lan+yhalo),
     &        rcs2_a(min(lat,latg-lat+1)),spdlat(1,iblk))
          enddo   ! omp lon
          iblk=1
          do lon=1,lons_lat,ngptcd
            do k=1,levs
               spdmax_node(k)=max(spdmax_node(k),spdlat(k,iblk))
            enddo
            iblk=iblk+1
          enddo

        enddo ! fin  lan loop 2

        spdmax_nodem=spdmax_node
        call mpi_gather(spdmax_nodem,levs,mpi_r_mpi,
     x                  spdmax_nodesm,levs,mpi_r_mpi,
     x                  0,mc_comp,ierr)
        spdmax_nodes=spdmax_nodesm
        if ( me .eq. 0 ) then
          do k=1,levs
            spdmax(k) = zero
            do node=1,nodes
               spdmax(k)=max(spdmax(k),spdmax_nodes(k,node))
            enddo
            spdmax(k)=sqrt(spdmax(k))
          enddo
          print 101,(spdmax(k),k=1,levs)
101       format(' v_max(01:10)=',10f5.0,:/' v_max(11:20)=',10f5.0,
     x         :/' v_max(21:30)=',10f5.0,:/' v_max(31:40)=',10f5.0,
     x         :/' v_max(41:50)=',10f5.0,:/' v_max(51:60)=',10f5.0,
     x         :/' v_max(61:70)=',10f5.0,:/' v_max(71:80)=',10f5.0,
     x         :/' v_max(81:90)=',10f5.0,:/' v_max(91:00)=',10f5.0)
        endif
      endif ! fin meridional velocity

      do k=1,levs
       spdmax_node(k) = zero
      enddo

      go_forward = .false.
      iter_slg   = 2
!
       if (.not. allocated(rgt_hl))
     &           allocate ( rgt_hl(lon_dim_h,levs,ntrac))
!

      if (kdt == 1 .or. first) then       
          iter_slg   = 4         !if(kdt .eq. 1 )
          go_forward = .true.
!       if (kdt == 1) then
!         iter_slg   = 4         !if(kdt .eq. 1 )
!         go_forward = .true.
!       endif
        ini_slg      = 1
        ini_dp       = 1
        print*,' first time step kdt,kdt_save = ',kdt,kdt_save,
     &  'iter_slg=',iter_slg
        if (me == 0) print*,'ngptcd = ',ngptcd,' ngptcs = ',ngptcs
        first = .false.
      endif

      if (kdt > 1 .and. kdt < kdt_save) then
        iter_slg   = 4   !if(kdt .gt. 1 .and. kdt .lt. kdt_save)
        go_forward = .true.
        ini_slg    = 0
        ini_dp     = 1
        if (me == 0)
     &  print*,' middle of filter kdt,kdt_save = ',kdt,kdt_save,
     &  'iter_slg=',iter_slg, 'get here only once at mid digifilter'      
      endif 

      if ( ini_slg == 1 .or. ini_dp == 1) then
            i_1 = 1  ! this call of slgscan must be made from 
            i_2 = 2  ! an unthreaded region  for initializing.

            call slgscan_h(i_1, i_2, s_lat_in_pe, n_lat_in_pe,
     &        j, lat, lon_dim_h,
     &        deltim, iter_slg, etamid, lats_node_a ,
     &        ug_h, vg_h, ww_h, lammp_a, phimp_a, sigmp_a,me,
     &        rgt_h(1,1,1,1),rgt_a(1,1,1,1) ,
     &        ud3_h1(1,1,1),ud3_h2(1,1,1),ud3_a(1,1,1) ,
     &        vd3_h1(1,1,1),vd3_h2(1,1,1),vd3_a(1,1,1) ,
     &        td3_h1(1,1,1),td3_h2(1,1,1),td3_a(1,1,1) ,
     &        dpdt_d3_h(1,1,1), dpdt_d3_a(1,1,1) ,
     &        zdp,
!    &        lam_all_dp, phi_all_dp, lam_all_ar, phi_all_ar,
     &        global_lats_a, lonsperlat, ini_slg, ini_dp,.false.)
      endif

       ini_slg  = 0
       ini_dp   = 0

       kdt_save = kdt

!      if (kdt < 3)
!    & print*,' before slgscan kdt,kdt_save = ',kdt,kdt_save,
!    & 'iter_slg=',iter_slg,'s_lat_in_pe=',s_lat_in_pe,
!    & 'n_lat_in_pe=',n_lat_in_pe
!    &,' lon_dim_h=',lon_dim_h,' lats_dim_h=',lats_dim_h,' me=',me
!    &,' lats_node_a=',lats_node_a
!     if (kdt == 1) call mpi_quit(444)

      do lan=1,lats_node_a   !sela begin lan loop 3
        timer1   = rtc()
        lat      = global_lats_a(ipt_lats_node_a+lats_node_a-lan)
        lons_lat = lonsperlat(lat)
        ylan     = yhalo + lan
        rcs2_loc =  rcs2_a(min(lat,latg-lat+1))

!       if (kdt == 1) print *,' lan=',lan,' rcs2i=',sqrt(1.0/rcs2_loc)
!    &,' lat=',lat,' ipt_lats_node_a=',ipt_lats_node_a,' lats_node_a='
!    &,lats_node_a,' lons_lat=',lons_lat,' me=',me

        if(hybrid) then !--------------  hybrid  ----------------------

!$omp parallel do schedule(dynamic,1) private(lon)
!$omp+private(njeff,iblk,xlon)
!$omp+private(nvcn,xvcn)

          do lon=1,lons_lat,ngptcd
            njeff = min(ngptcd,lons_lat-lon+1)
            iblk  = (lon-1)/ngptcd+ 1
            xlon  = xhalo + lon
            call gfidi_hyb_resonan
     &             (lon_dim_a-2, lon_dim_h, njeff, lat,
     &        for_grid_a( lon,ksd   , lan ),
     &        for_grid_a( lon,kst   , lan ),
     &              ug_h(xlon,1     ,ylan),
     &              vg_h(xlon,1     ,ylan),
     &              ug_m(xlon,1     ,ylan),
     &              vg_m(xlon,1     ,ylan),
     &              ww_h(xlon,1     ,ylan),
     &        for_grid_a( lon,ksrg  , lan),
     &        for_grid_a( lon,kspphi, lan),
     &        for_grid_a( lon,ksplam, lan),
     &        for_grid_a( lon,ksq   , lan),
     &            rcs2_loc,
     &            etamid,etaint,
     &            spdlat(1,iblk),
     &            deltim,nvcn,xvcn,
     &        dyn_grid_a( lon,kdtphi, lan),
     &        dyn_grid_a( lon,kdtlam, lan),
     &             td3_h1(xlon,1     ,ylan),
     &             td3_h2(xlon,1     ,ylan),
     &             rgt_hl(xlon,1     ,1   ),
     &             ud3_h1(xlon,1     ,ylan),
     &             ud3_h2(xlon,1     ,ylan),
     &             vd3_h1(xlon,1     ,ylan),
     &             vd3_h2(xlon,1     ,ylan),
     &               sd_m( lon       , lan),
     &                z_m( lon,1     , lan),
     &               unlm( lon,1     , lan),
     &               vnlm( lon,1     , lan),
     &               tnlm( lon,1     , lan),
     &               unla( lon,1     , lan),
     &               vnla( lon,1     , lan),
     &               tnla( lon,1     , lan),
     &             dpdt_a( lon       , lan),
     &          dpdt_d3_h(xlon,1     ,ylan),
     &         grad_gzlam( lon       , lan),
     &         grad_gzphi( lon       , lan),
     &            gz_grid( lon       , lan),
     &             go_forward,
     &             batah,gg_tracers,tb_k,kdt,r_rt,omega_v(lan))
          enddo   ! omp lon
!
!     if (me == 1 .and. lat == 48) then
 !    if (lat == 47) then
!     if ( lan == 21 .and.  me == 2) then
!          print *,'OZGFIDI me=',me,' rgt_h=',rgt_hl(xhalo+1,1:5,2)
!    &,       ' lons_lat=',lons_lat,xhalo,yhalo,lan
!     endif

          do ntr=1,ntrac
            do k=1,levs

!!$omp parallel do schedule(dynamic,1) private(lon,iblk)
              do lon=1,lons_lat
                xlon  = xhalo + lon
                rgt_h(xlon,k,ylan,ntr) = rgt_hl(xlon,k,ntr)
              enddo
            enddo
          enddo

!     if (me == 1 .and. lat == 48) then
!     if (lat == 47) then
!     if ( lan == 21 .and.  me == 2) then
!     lanx = lan
!     print *,' OZONE aft gfidi for me=',me,' rgt_h=',
!    & rgt_h(xhalo+1,1:5,yhalo+lan,2),' lan=',lan,' lat=',lat
!    &,' rcs2i=',sqrt(1.0/rcs2_loc)
!     endif

        else  !-------------  sigma branch ---------------------------

!      It appears that the sigma option is not coded - Moorthi
!      -------------------------------------------------------
!$omp parallel do schedule(dynamic,1) private(lon)
!$omp+private(njeff,iblk)
!$omp+private(nvcn,xvcn)
          do lon=1,lons_lat,ngptcd
            njeff = min(ngptcd,lons_lat-lon+1)
            iblk  = (lon-1)/ngptcd+ 1
          enddo   !lon
        endif ! -----------------------  hybrid  ------------------

        iblk = 1
        do lon=1,lons_lat,ngptcd
           do k=1,levs
             spdmax_node(k) = max(spdmax_node(k),spdlat(k,iblk))
           enddo
          iblk = iblk + 1
        enddo
        timer2 = rtc()
        global_times_a(lat,me+1)=global_times_a(lat,me+1)+timer2-timer1

      enddo !sela fin lan loop 3

!pl   if (kdt .eq. 1 .and. nodes .eq. 1) then
!pl    call plot_85_a_sn(tnlm,'tnla',lonf,latg,levs,kdt,800)
!pl    call plot_85_a_sn(unlm,'unla',lonf,latg,levs,kdt,900)
!pl    call plot_85_a_sn(vnlm,'vnla',lonf,latg,levs,kdt,1000)
!pl    print*,' completed plotting after gfidi kdt=',kdt
!pl    stop
!pl   endif


!$$$  to use these prints halo quantities must be zeroed before gfidi
!$$$      do j=1,lats_dim_h   !  lan loop scan
!$$$        k = j
!$$$        print*, 'lat,k ud3_h = ',lat,k
!$$$        print*,'b ud3_h = ',k,minval(ud3_h(:,:,k)),maxval(ud3_h(:,:,k))
!$$$        print*,'b vd3_h = ',k,minval(vd3_h(:,:,k)),maxval(vd3_h(:,:,k))
!$$$        print*,'b td3_h = ',k,minval(td3_h(:,:,k)),maxval(td3_h(:,:,k))
!$$$        print*,'b dpdt_d3_h = ',k,minval(dpdt_d3_h(:,:,k)),
!$$$     .          maxval(dpdt_d3_h(:,:,k))
!$$$        print*,'b ug_h = ',k,minval(ug_h(:,:,k)),maxval(ug_h(:,:,k))
!$$$        print*,'b vg_h = ',k,minval(vg_h(:,:,k)),maxval(vg_h(:,:,k))
!$$$      enddo
!$$$      stop

      call  set_halos(ug_h,
     &               lats_nodes_h,global_lats_h,
     &               lonsperlat,
     &               lon_dim_h,xhalo,yhalo,
     &               levs,
     &               levs)
      call  set_halos(vg_h,
     &               lats_nodes_h,global_lats_h,
     &               lonsperlat,
     &               lon_dim_h,xhalo,yhalo,
     &               levs,
     &               levs)
      call set_halos(ww_h,
     &               lats_nodes_h,global_lats_h,
     &               lonsperlat,
     &               lon_dim_h,xhalo,yhalo,
     &               levs,
     &               levs)
      call set_halos(ud3_h1,
     &               lats_nodes_h,global_lats_h,
     &               lonsperlat,
     &               lon_dim_h,xhalo,yhalo,
     &               levs,
     &               levs)
      call set_halos(ud3_h2,
     &               lats_nodes_h,global_lats_h,
     &               lonsperlat,
     &               lon_dim_h,xhalo,yhalo,
     &               levs,
     &               levs)
      call set_halos(vd3_h1,
     &               lats_nodes_h,global_lats_h,
     &               lonsperlat,
     &               lon_dim_h,xhalo,yhalo,
     &               levs,
     &               levs)
      call set_halos(vd3_h2,
     &               lats_nodes_h,global_lats_h,
     &               lonsperlat,
     &               lon_dim_h,xhalo,yhalo,
     &               levs,
     &               levs)
      call set_halos(td3_h1,
     &               lats_nodes_h,global_lats_h,
     &               lonsperlat,
     &               lon_dim_h,xhalo,yhalo,
     &               levs,
     &               levs)
      call set_halos(td3_h2,
     &               lats_nodes_h,global_lats_h,
     &               lonsperlat,
     &               lon_dim_h,xhalo,yhalo,
     &               levs,
     &               levs)

!     if (me == 2) then
!     if (me == 1) then
!     print *,' OZONE BEF set_halos for me=',me,' rgt_h=',
!    & rgt_h(xhalo+1,1:5,yhalo+lanx,2)
!    &,' lon_dim_h=',lon_dim_h,' lats_dim=',lats_dim_a
!     endif

      do ntr=1,ntrac
        call set_halos(rgt_h(1,1,1,ntr),
     &                 lats_nodes_h,global_lats_h,
     &                 lonsperlat,
     &                 lon_dim_h,xhalo,yhalo,
     &                 levs,
     &                 levs)
      enddo

!     if (me == 1 .and. lat == 48) then
!     if (me == 2) then
!     if (me == 1) then
!     print *,' OZONE aft set_halos for me=',me,' rgt_h=',
!    & rgt_h(xhalo+1,1:5,yhalo+lanx,2)
!     endif

      call set_halos(dpdt_d3_h,
     &               lats_nodes_h,global_lats_h,
     &               lonsperlat,
     &               lon_dim_h,xhalo,yhalo,
     &               levs,
     &               levs)
200   format(' n.hemi me=',i3,2x,' n_ggp=',i3)
201   format(' s.hemi me=',i3,2x,' s_ggp=',i3)
      if (me == 0) then
        n_ggp = lats_nodes_h(1) - yhalo
        if (kdt == 1)print 200,me,n_ggp
        call sltini_n_uvw(n_ggp,lon_dim_h,ug_h,vg_h,ww_h,yhalo,
     &                     lonsperlat(latg))                     
        call sltini_n_uvw(n_ggp,lon_dim_h,ud3_h1,vd3_h1,ww_h,yhalo,
     &                     lonsperlat(latg))                     
cmy redundant for ww_h - for testing fix later
        call sltini_n_uvw(n_ggp,lon_dim_h,ud3_h2,vd3_h2,ww_h,yhalo,
     &                     lonsperlat(latg))                     
        call sltini_nh_scalar(n_ggp,lon_dim_h,td3_h1(1,1,1),yhalo,
     &                     lonsperlat(latg))                     
        call sltini_nh_scalar(n_ggp,lon_dim_h,td3_h2(1,1,1),yhalo,
     &                     lonsperlat(latg))                     
        call sltini_nh_scalar(n_ggp,lon_dim_h,dpdt_d3_h(1,1,1),yhalo,
     &                     lonsperlat(latg))                     
        do ntr=1,ntrac
          call sltini_nh_scalar(n_ggp,lon_dim_h,rgt_h(1,1,1,ntr),yhalo,
     &                     lonsperlat(latg))                     
        enddo
      endif
      if (me == nodes-1) then
        s_ggp=yhalo+1
        if (kdt == 1) print 201,me,s_ggp
        call sltini_s_uvw(s_ggp,lon_dim_h,ug_h,vg_h,ww_h,yhalo,
     &                     lonsperlat(1))                     
        call sltini_s_uvw(s_ggp,lon_dim_h,ud3_h1,vd3_h1,ww_h,yhalo,
     &                     lonsperlat(1))                     
cmy redundant for ww_h - for testing fix later
        call sltini_s_uvw(s_ggp,lon_dim_h,ud3_h2,vd3_h2,ww_h,yhalo,
     &                     lonsperlat(1))                     
        call sltini_sh_scalar(s_ggp,lon_dim_h,td3_h1(1,1,1),yhalo,
     &                     lonsperlat(1))                     
        call sltini_sh_scalar(s_ggp,lon_dim_h,td3_h2(1,1,1),yhalo,
     &                     lonsperlat(1))                     
        call sltini_sh_scalar(s_ggp,lon_dim_h,dpdt_d3_h(1,1,1),yhalo,
     &                     lonsperlat(1))                     
        do ntr=1,ntrac
          call sltini_sh_scalar(s_ggp,lon_dim_h,rgt_h(1,1,1,ntr),yhalo,
     &                     lonsperlat(1))                     
        enddo
      endif


!sela if(kdt.eq.1) then
!sela do j=1,lats_dim_h   !  lan loop scan
!sela
!sela   k = j
!sela   print*, 'lat,k ud3_h = ',lat,k
!sela   print*,'m ud3_h = ',k,minval(ud3_h(:,:,k)),maxval(ud3_h(:,:,k))
!sela   print*,'m vd3_h = ',k,minval(vd3_h(:,:,k)),maxval(vd3_h(:,:,k))
!sela   print*,'m td3_h = ',k,minval(td3_h(:,:,k)),maxval(td3_h(:,:,k))
!sela   print*,'m dpdt_d3_h = ',k,minval(dpdt_d3_h(:,:,k)),
!sela.          maxval(dpdt_d3_h(:,:,k))
!sela   print*,'m ug_h = ',k,minval(ug_h(:,:,k)),maxval(ug_h(:,:,k))
!sela   print*,'m vg_h = ',k,minval(vg_h(:,:,k)),maxval(vg_h(:,:,k))
!sela enddo
!sela stop
!sela endif
      
      if (me == 0 .and. kdt < 3) print*,'iter_slg = ',iter_slg 
!     if (me == 0) print*,'before call to slgscan iter_slg=',iter_slg 

!     if (kdt == 1) call mpi_quit(444)

      rgt_a = -9999.0
      do j=1,lats_node_a   !  lan loop scan
         lat      = latg+1-global_lats_a(ipt_lats_node_a+lats_node_a-j)
         lons_lat = lonsperlat(lat)

!$omp parallel do schedule(dynamic,1) private(lon)
!$omp+private(njeff,i_1,i_2)
         do lon=1,lons_lat,ngptcs
            njeff = min(ngptcs,lons_lat-lon+1)
            i_1   = lon
            i_2   = lon+njeff-1

!           if (kdt == 1 ) then
!           lprnt = latg+1-lat == 47 .and. lon == 1
!           lprnt = j == 21 .and.  me == 2 .and. lon == 1
!           lprnt = j == 19 .and.  me == 0 .and. lon == 1
!           lprnt = me == 0 .and. lon == 1
!           if (lprnt) then
!             print*,' before slgscan latj,i_1,i_2=',j,i_1,i_2, '  me=',
!    &                 me,' njeff=',njeff,' s_lat_in_pe=',s_lat_in_pe,
!    &' n_lat_in_pe=',n_lat_in_pe
!             print *,' lat=',lat,' rgt_h=',rgt_h(2,1:5,yhalo+j,2)
!    &,' me=',me,' lon=',lon,' lat=',lat,' lats_node_a=',lats_node_a
 
!           endif
!        print *,' BEFORE SLGSCAN rgt_htop=',rgt_h(1,levs,j,2),' j=',j

            call slgscan_h(i_1,i_2,s_lat_in_pe,n_lat_in_pe,
     &        j,lat,lon_dim_h,
     &        deltim , iter_slg , etamid  , lats_node_a ,
     &        ug_h,vg_h,ww_h,lammp_a,phimp_a,sigmp_a,me,
     &        rgt_h(1,1,1,1),rgt_a(1,1,1,1) ,
     &        ud3_h1(1,1,1),ud3_h2(1,1,1),ud3_a(1,1,1) ,
     &        vd3_h1(1,1,1),vd3_h2(1,1,1),vd3_a(1,1,1) ,
     &        td3_h1(1,1,1),td3_h2(1,1,1),td3_a(1,1,1) ,
     &        dpdt_d3_h(1,1,1),dpdt_d3_a(1,1,1) ,
     &        zdp,
!    &        lam_all_dp,phi_all_dp,lam_all_ar,phi_all_ar,
     &        global_lats_a,lonsperlat,ini_slg,ini_dp,lprnt)
         enddo   ! lon
           if (kdt.eq.1) print*,' completed slgscan j= ',j,' kdt=',kdt
     &,' me=',me

!        if (me == 0) then
!        print *,' AFTER SLGSCAN rgt_atop=',rgt_a(1,levs,j,2),' j=',j
!    &,' me=',me
!        print *,' lats_dim_h=',lats_dim_h,s_lat_in_pe,n_lat_in_pe
!        print *,' OZONE rgt_a=',rgt_a(1,levs-4:levs,j,2),' me=',me,
!    &' latg+1-lat=',latg+1-lat,' j=',j
!        endif

      enddo   !  fin lan loop scan

!     print *,' OZONE AT THE TOP : me=',me,'  rgt_a=',rgt_a(1,levs,:,2)

      hfdt = 0.5*deltim

      do lan=1,lats_node_a                      !sela begin lan loop 5
         jtem = lats_node_a+1-lan
         do k=1,levs
           ktem = levs+1-k
           do i=1,lonf
!sela       bak_grid_a(i,kau-1+k)= -omega_v(lan)+
             bak_grid_a(i,kau-1+k) = 
     &            (ud3_a(i,k,jtem) + hfdt*unla(i,ktem,jtem))*rcos1(lan)

             bak_grid_a(i,kav-1+k) = 
     &            (vd3_a(i,k,jtem) + hfdt*vnla(i,ktem,jtem))*rcos1(lan)

             bak_grid_a(i,kat-1+k) = td3_a(i,k,jtem)
     &                             + hfdt*tnla(i,ktem,jtem)

           enddo
         enddo
         if(.not.gg_tracers .or. lsout)then
           do ntr=1,ntrac
             nn = (ntr-1)*levs
             do k=1,levs
               kk = kar-1+nn+k
               do i=1,lonf
                 bak_grid_a(i,kk) = rgt_a(i,k,jtem,ntr)
               enddo
             enddo
           enddo
         endif
         do i=1,lonf
           bak_grid_a(i,kap) = dpdt_a(i,jtem)
         enddo
         do k=1,levs
           do i=1,lonf
             bak_grid_a(i,kap) = bak_grid_a(i,kap) +
     &                           dpdt_d3_a(i,levp1-k,jtem)
           enddo
         enddo
 
         lat      = latg+1 - global_lats_a(ipt_lats_node_a-1+lan)
         lons_lat = lonsperlat(lat)
!sela    if(kdt.eq.1)then
!sela    print*,' before grid_to_four lan,lat,lon_dim_a,lons_lat = ',
!sela.   lan,lat,lon_dim_a,lons_lat
!sela    endif

         call grid_to_four(bak_grid_a(1,kau),
     &                     bak_coef_a(1,kau,lan),
     &                     lon_dim_a-2,lon_dim_a,lons_lat,3*levs+1)

         if (.not. gg_tracers .or. lsout) then
            call grid_to_four(bak_grid_a(1,kar),
     &                        bak_coef_a(1,kar,lan),
     &                        lon_dim_a-2,lon_dim_a,lons_lat,levh)
         endif

      enddo                            ! fin lan loop 5
      call synctime()

      call four2fln_gg(lats_dim_a,lota,3*levs+1,bak_coef_a,
     &              ls_nodes,max_ls_nodes,
     &              lats_nodes_a,global_lats_a,lon_dim_a,
     &              lats_node_a,ipt_lats_node_a,
     &              lat1s_a,lon_dim_a,latg,latg2,
     &              trie_ls(1,1,p_w), trio_ls(1,1,p_w),
     &              plnew_a, plnow_a,
     &              ls_node,0,2*levs+1,3*levs+1)

      if(.not.gg_tracers .or. lsout)then
         call four2fln_gg(lats_dim_a,lota,levh,bak_coef_a,
     &              ls_nodes,max_ls_nodes,
     &              lats_nodes_a,global_lats_a,lon_dim_a,
     &              lats_node_a,ipt_lats_node_a,
     &              lat1s_a,lon_dim_a,latg,latg2,
     &              trie_ls(1,1,p_rt), trio_ls(1,1,p_rt),
     &              plnew_a, plnow_a,
     &              ls_node,3*levs+1,1,levh)
      endif

!$omp parallel do shared(trie_ls,trio_ls)
!$omp+shared(epse,epso,ls_node)
!$omp+private(k)
      do k=1,levs
         call uveodz(trie_ls(1,1,p_w  +k-1), trio_ls(1,1,p_x  +k-1),
     x               trie_ls(1,1,p_uln+k-1), trio_ls(1,1,p_vln+k-1),
     x               epse,epso,ls_node)
         call uvoedz(trio_ls(1,1,p_w  +k-1), trie_ls(1,1,p_x  +k-1),
     x               trio_ls(1,1,p_uln+k-1), trie_ls(1,1,p_vln+k-1),
     x               epse,epso,ls_node)
      enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!sela add tb to Y to complete temp. tendency. WATCH SIGN!!!!!!!!

      do k=1,levs
        ktem = p_y+k-1
        kk   = levp1-k
        do i=1,len_trie_ls
          trie_ls(i,1,ktem) = trie_ls(i,1,ktem) + tb_k(kk)*phige(i,1)
          trie_ls(i,2,ktem) = trie_ls(i,2,ktem) + tb_k(kk)*phige(i,2)
        enddo
        do i=1,len_trio_ls
          trio_ls(i,1,ktem) = trio_ls(i,1,ktem) + tb_k(kk)*phigo(i,1)
          trio_ls(i,2,ktem) = trio_ls(i,2,ktem) + tb_k(kk)*phigo(i,2)
        enddo
      enddo
!sela add tb to Y to complete temp. tendency. WATCH SIGN!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (r_rt == 0.) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do k=1,levs
          ktem = p_x  +k-1
          jtem = p_uln+k-1
          do i=1,len_trie_ls
            trie_ls(i,1,ktem) = trie_ls(i,1,jtem)
            trie_ls(i,2,ktem) = trie_ls(i,2,jtem)
          enddo
          do i=1,len_trio_ls
            trio_ls(i,1,ktem) = trio_ls(i,1,jtem)
            trio_ls(i,2,ktem) = trio_ls(i,2,jtem)
          enddo
        enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      else
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do k=1,levs
          ktem = p_x  +k-1
          jtem = p_uln+k-1
          do i=1,len_trie_ls
            trie_ls(i,1,ktem) = trie_ls(i,1,jtem)
     &                        - batah*hfdt *trie_ls(i,1,P_gz)
            trie_ls(i,2,ktem) = trie_ls(i,2,jtem)
     &                        - batah*hfdt *trie_ls(i,2,P_gz)
          enddo
          do i=1,len_trio_ls
            trio_ls(i,1,ktem) = trio_ls(i,1,jtem)
     &                        - batah*hfdt *trio_ls(i,1,P_gz)
            trio_ls(i,2,ktem) = trio_ls(i,2,jtem)
     &                        - batah*hfdt *trio_ls(i,2,P_gz)
          enddo
        enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do k=1,levs
        ktem = p_w  +k-1
        jtem = p_vln+k-1
        do i=1,len_trie_ls
          trie_ls(i,1,ktem) = trie_ls(i,1,jtem)
          trie_ls(i,2,ktem) = trie_ls(i,2,jtem)
        enddo
        do i=1,len_trio_ls
          trio_ls(i,1,ktem) = trio_ls(i,1,jtem)
          trio_ls(i,2,ktem) = trio_ls(i,2,jtem)
        enddo
      enddo
      do locl=1,ls_max_node
              l=ls_node(locl,1)
         jbasev=ls_node(locl,2)
         if ( l == 0 ) then
            n = 0
            do k=1,levs
               trie_ls(indlsev(n,l),1,p_w+k-1) = zero
               trie_ls(indlsev(n,l),2,p_w+k-1) = zero
            enddo
         endif
      end do
      spdmax_nodem = spdmax_node
      call mpi_gather(spdmax_nodem,levs,mpi_r_mpi,
     x                spdmax_nodesm,levs,mpi_r_mpi,
     x                0,mc_comp,ierr)
      spdmax_nodes = spdmax_nodesm
      if ( me .eq. 0 ) then
         do k=1,levs
            spdmax(k) = zero
            do node=1,nodes
               spdmax(k) = max(spdmax(k),spdmax_nodes(k,node))
            enddo
            spdmax(k) = sqrt(spdmax(k))
         enddo
!       print*,'in gloopa kdt=',kdt
!        print 100,(spdmax(k),k=1,levs)
!100     format(' spdmx(01:10)=',10f5.0,:/' spdmx(11:20)=',10f5.0,
!    x        :/' spdmx(21:30)=',10f5.0,:/' spdmx(31:40)=',10f5.0,
!    x        :/' spdmx(41:50)=',10f5.0,:/' spdmx(51:60)=',10f5.0,
!    x        :/' spdmx(61:70)=',10f5.0,:/' spdmx(71:80)=',10f5.0,
!    x        :/' spdmx(81:90)=',10f5.0,:/' spdmx(91:00)=',10f5.0)
      endif
      if (lsout) then
        do k=1,levs
          n_rgt_a = 0

          do lan=1,lats_node_a   !sela begin neg. moist loop
            lat      = global_lats_a(ipt_lats_node_a-1+lan)
            lons_lat = lonsperlat(lat)
            do lon=1,lons_lat
              if(rgt_a(lon,k,lan,1).lt.0.)then
                n_rgt_a = n_rgt_a + 1
              endif
            enddo   ! omp lon
          enddo     ! fin  lat loop

!         if(n_rgt_a.gt.0)then
!sela print*,'kdt=',kdt,' k=',k,' number of negative point =',n_rgt_a
!         endif 

        enddo       ! fin  levs loop
      endif         !lsout

      call mpi_bcast(spdmax,levs,mpi_real8,
     x               0,mc_comp,ierr)
      deallocate ( spdlat )
      deallocate ( for_coef_a )
      deallocate ( dyn_coef_a )
      deallocate ( bak_coef_a )
      deallocate ( for_grid_a )
      deallocate ( dyn_grid_a )
      deallocate ( bak_grid_a )
      deallocate ( ud3_a )
      deallocate ( vd3_a )
      deallocate ( td3_a )
      deallocate ( dpdt_a )
      deallocate ( dpdt_d3_a )
      deallocate ( ug_h )
      deallocate ( vg_h )
      deallocate ( ww_h )
      deallocate ( ud3_h1 )
      deallocate ( ud3_h2 )
      deallocate ( vd3_h1)
      deallocate ( vd3_h2)
      deallocate ( td3_h1 )
      deallocate ( td3_h2 )
      deallocate ( dpdt_d3_h )
      deallocate ( unla,vnla,tnla)
!     deallocate ( lam_all_dp,phi_all_dp)
!     deallocate ( lam_all_ar,phi_all_ar)
      deallocate ( rgt_hl)
!     print*,' completed gloopa_tracers kdt=',kdt
      return
      end
