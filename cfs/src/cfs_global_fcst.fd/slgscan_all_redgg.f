       subroutine slgscan_h(i_1, i_2, s_lat_in_pe, n_lat_in_pe,
     &  lan, lat, lon_dim_h,
     &  ztodt  , iter ,  etamid  , lats_node_a,
     &  uu_gr_h_2, vv_gr_h_2, ww_gr_h_2, lammp_a, phimp_a, sigmp_a, me,
     &  rgt_h, rgt_a, ud3_h1, ud3_h2, ud3_a, vd3_h1, vd3_h2, vd3_a, 
     &  td3_h1, td3_h2, td3_a, dpdt_d3_h, dpdt_d3_a, zdp_grid,
!    &  lam_all_dp, phi_all_dp, lam_all_ar, phi_all_ar,
     &  global_lats_a, lonsperlat, ini_slg, ini_dp,lprint)
!
      use machine             , only : kind_phys
      use resol_def           , only : latg,levs,lonf,ntrac
      use namelist_def        , only : redgg_a, phigs
     &,                                herm_x, herm_y, herm_z
      use pmgrid              , only :      plon,platd,
!sela use pmgrid              , only : pmap,plon,platd,
     &                                 plev,plat,plevp
      use slgshr              , only : phi,dphi,lbasiy,lam,rdlam,
     &                                 dlam,ra,lbasdy,nlonex,rdlam6
      use gg_def              , only : rcs2_a
!Moor use layout_grid_tracers , only : yhalo
      use layout1             , only : ipt_lats_node_a,lats_dim_a

      implicit none
      integer i_1,i_2,s_lat_in_pe,n_lat_in_pe
      integer lan,lat,lon_dim_h,lats_node_a
      real    ztodt, etamid(levs)           ! eta at levels
      real    lam0,pi,cos_f_a,sin_f_a ,mindetam,twopi,degrad
      real    cos_l_a(i_1:i_2), sin_l_a(i_1:i_2)

      integer             global_lats_a(latg)
      integer             lonsperlat(latg)
      integer             ini_slg,ini_dp
      integer             lan_loc,lat_loc
      real                dlam_con
      logical lprint

      integer ,allocatable :: kdpmpf(:)      ! mapping from artificial array
                                             ! to model level
      real    ,allocatable :: lbasdz(:,:,:)  ! basis funcs for vert
                                             ! deriv est. at level
      real    ,allocatable :: lbasiz(:,:,:)  ! lagrange cubic interp
                                             ! wghts (vert)
      real    ,allocatable :: detam(:)       ! delta eta at levels

      save kdpmpf, lbasdz, lbasiz, detam, rdel, pmap

!     ------------------------------------------------------------------
!     for  grdini_h


      real udmin,udmax
      integer jj, jj1, k1

!     ------------------------------------------------------------------
!     for  lcbas_h

      real x0mx1,               ! |
     $     x0mx2,               ! |
     $     x0mx3,               ! |- grid value differences used in
     $     x1mx2,               ! |  weights
     $     x1mx3,               ! |
     $     x2mx3                ! |

!     ------------------------------------------------------------------
!     for  lcdbas_h

      real x_1,                  ! |
     $     x_2,                  ! |- grid values
     $     x_3,                  ! |
     $     x_4,                  ! |
!    $     x1mx2,               !  |
!    $     x1mx3,               !  |
     $     x1mx4,               !  |- differences of grid values
!    $     x2mx3,               !  |
     $     x2mx4,               !  |
     $     x3mx4                !  |

      real x1mx2i, x1mx3i, x1mx4i, x2mx3i, x2mx4i, x3mx4i

!     ------------------------------------------------------------------
!     for  eta_stuff

      real    del,leps,dp

      integer imin
!     integer lats_dim_h

!     ------------------------------------------------------------------

      real   uu_gr_h_2(lon_dim_h*levs ,s_lat_in_pe:n_lat_in_pe)
      real   vv_gr_h_2(lon_dim_h*levs ,s_lat_in_pe:n_lat_in_pe)
      real   ww_gr_h_2(lon_dim_h*levs ,s_lat_in_pe:n_lat_in_pe)

      real         rgt_h(lon_dim_h,levs ,s_lat_in_pe:n_lat_in_pe,ntrac)
      real         ud3_h1(lon_dim_h,levs ,s_lat_in_pe:n_lat_in_pe)
      real         ud3_h2(lon_dim_h,levs ,s_lat_in_pe:n_lat_in_pe)
      real         vd3_h1(lon_dim_h,levs ,s_lat_in_pe:n_lat_in_pe)
      real         vd3_h2(lon_dim_h,levs ,s_lat_in_pe:n_lat_in_pe)
      real         td3_h1(lon_dim_h,levs ,s_lat_in_pe:n_lat_in_pe)
      real         td3_h2(lon_dim_h,levs ,s_lat_in_pe:n_lat_in_pe)
      real     dpdt_d3_h(lon_dim_h,levs ,s_lat_in_pe:n_lat_in_pe)

      real  lammp_a(lonf,levs,lats_node_a)
      real  phimp_a(lonf,levs,lats_node_a)
      real  sigmp_a(lonf,levs,lats_node_a)

!     real lam_all_dp(lonf,levs,lats_node_a)
!     real phi_all_dp(lonf,levs,lats_node_a)
!     real lam_all_ar(lonf,levs,lats_node_a)
!     real phi_all_ar(lonf,levs,lats_node_a)

      real     rgt_a(lonf,levs,lats_dim_a,ntrac)
      real     ud3_a(lonf,levs,lats_node_a)
      real     vd3_a(lonf,levs,lats_node_a)
      real     td3_a(lonf,levs,lats_node_a)
      real dpdt_d3_a(lonf,levs,lats_node_a)

      real lamdp(i_1:i_2,levs),       ! x-coord of dep pt
     &     phidp(i_1:i_2,levs),       ! y-coord of dep pt
     &     sigdp(i_1:i_2,levs)        ! z-coord of dep pt

      real cos_l_d(i_1:i_2,levs),
     &     sin_l_d(i_1:i_2,levs),
     &     cos_f_d(i_1:i_2,levs),
     &     sin_f_d(i_1:i_2,levs)

      integer i,j,k,n,jcen,iter,me,i_branch
      integer idp(i_1:i_2,levs,4),jdp(i_1:i_2,levs),kdp(i_1:i_2,levs)
      integer kkdp(i_1:i_2,plev),kdim,kdimm1,kdimm2

      real    rdel,dphibr,phibs
      real   fintx(i_1:i_2,plev,4,4),finty(i_1:i_2,plev,4)
      real     rbi(i_1:i_2,plev),     rbj(i_1:i_2,plev),
     &         rmm(i_1:i_2,plev)
      real     rww(i_1:i_2,plev),     rnn(i_1:i_2,plev)
      real  term1x(i_1:i_2,plev,4),term2x(i_1:i_2,plev,4),
     $      term3x(i_1:i_2,plev,4),term4x(i_1:i_2,plev,4)
      real      x2(i_1:i_2,plev,4),    x3(i_1:i_2,plev,4)
      real  term1y(i_1:i_2,plev),  term2y(i_1:i_2,plev),
     &      term3y(i_1:i_2,plev),  term4y(i_1:i_2,plev),
     &          yb(i_1:i_2,plev),      yt(i_1:i_2,plev)
      real      ht(i_1:i_2,plev),      hb(i_1:i_2,plev),
     &         dht(i_1:i_2,plev),   dhb(i_1:i_2,plev),
     $          zt(i_1:i_2,plev),    zb(i_1:i_2,plev),
     &      term1z(i_1:i_2,plev),term2z(i_1:i_2,plev),
     &      term3z(i_1:i_2,plev),term4z(i_1:i_2,plev),
     &       ud3_d(lonf,plev),vd3_d(lonf,plev)         


      integer kk, ierr, pmap, bisection, nt, jp

      real zdp_grid(plon,levs,plat)
      real ys(i_1:i_2,plev), yn(i_1:i_2,plev),
     &     hs(i_1:i_2,plev), hn(i_1:i_2,plev),
     &    dhs(i_1:i_2,plev),dhn(i_1:i_2,plev),
     &  rdphi(i_1:i_2,plev)
      real xl(i_1:i_2,plev,4), xr(i_1:i_2,plev,4),
     &     hl(i_1:i_2,plev,4), hr(i_1:i_2,plev,4),
     &    dhl(i_1:i_2,plev,4),dhr(i_1:i_2,plev,4)

      logical her_x,her_y,her_z,lin_xyz,wgt_cub_lin_xyz,lprnt
      real    tem
      parameter (lin_xyz=.true.,wgt_cub_lin_xyz=.true.)


      integer     i_count
      save        i_count
      data        i_count   / 0 /

      i_count = i_count + 1
!
      pi     = 4.*atan(1.)
      twopi  = pi + pi
      degrad = 180./pi

!    -------------------------------------------------------------------

      if ( ini_slg == 1 ) then

        if (.not.allocated(kdpmpf)) allocate (kdpmpf(plev-1))
        if (.not.allocated(lbasdz)) allocate (lbasdz(4,2,plev))
        if (.not.allocated(lbasiz)) allocate (lbasiz(4,2,plev))
        if (.not.allocated(lbasdy)) allocate (lbasdy(4,2,platd))
        if (.not.allocated(lbasiy)) allocate (lbasiy(4,2,platd))
        if (.not.allocated(phi))    allocate (phi(platd))
        if (.not.allocated(dphi))   allocate (dphi(platd))
        if (.not.allocated(detam))  allocate (detam(plev))

!       if (me == 0) print*,' calling grdini_h from slgscan_h '
!       call grdini_h(rcs2_a,lonsperlat)

!     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

!       inline  grdini_h

        print*,' reduced grid inside grdini_h  '


        do i=1,plat/2
          tem = acos(sqrt(1./rcs2_a(i)))
          phi(i+2)       = -tem
          phi(platd-1-i) =  tem
        enddo

        phi(1)      = -pi - phi(3)
        phi(2)      = -pi*0.5

        phi(plat+3) =  pi*0.5
        phi(platd)  = pi - phi(plat+2)

!       do j=1,platd
!          write(6,100)j,phi(j)*degrad
!       enddo
100     format(' in grdini    lat=',i4,2x,' phi=',f9.3)

        do jj = 2,plat+2
!         call lcdbas_h( phi(jj-1), lbasdy(1,1,jj), lbasdy(1,2,jj) )

!     . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 

!      inline  lcdbas_h(grd,dbas2,dbas3)

          jj1 = jj - 1
          x_1 = phi(jj1+0)
          x_2 = phi(jj1+1)
          x_3 = phi(jj1+2)
          x_4 = phi(jj1+3)

          x1mx2 = x_1 - x_2
          x1mx3 = x_1 - x_3
          x1mx4 = x_1 - x_4
          x2mx3 = x_2 - x_3
          x2mx4 = x_2 - x_4
          x3mx4 = x_3 - x_4

          x1mx2i = 1.0 / x1mx2
          x1mx3i = 1.0 / x1mx3
          x1mx4i = 1.0 / x1mx4
          x2mx3i = 1.0 / x2mx3
          x2mx4i = 1.0 / x2mx4
          x3mx4i = 1.0 / x3mx4

          lbasdy(1,1,jj) =   x2mx3 * x2mx4 * x1mx2i * x1mx3i * x1mx4i
          lbasdy(2,1,jj) = - x1mx2i + x2mx3i + x2mx4i
          lbasdy(3,1,jj) = - x1mx2 * x2mx4 * x1mx3i * x2mx3i * x3mx4i
          lbasdy(4,1,jj) =   x1mx2 * x2mx3 * x1mx4i * x2mx4i * x3mx4i

          lbasdy(1,2,jj) = - x2mx3 * x3mx4 * x1mx2i * x1mx3i * x1mx4i
          lbasdy(2,2,jj) =   x1mx3 * x3mx4 * x1mx2i * x2mx3i * x2mx4i
          lbasdy(3,2,jj) =   -x1mx3i - x2mx3i + x3mx4i
          lbasdy(4,2,jj) = - x1mx3 * x2mx3 * x1mx4i * x2mx4i * x3mx4i


!     . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 

!         call lcbas_h( phi(jj-1), lbasiy(1,1,jj), lbasiy(1,2,jj) )

!     . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 

!     inline  lcbas_h( grd, bas1, bas2 )

!
          lbasiy(1,1,jj) = x_1
          lbasiy(2,1,jj) = x_2
          lbasiy(3,1,jj) = x_3
          lbasiy(4,1,jj) = x_4

          lbasiy(1,2,jj) =   x1mx2i * x1mx3i * x1mx4i
          lbasiy(2,2,jj) = - x1mx2i * x2mx3i * x2mx4i
          lbasiy(3,2,jj) =   x1mx3i * x2mx3i * x3mx4i
          lbasiy(4,2,jj) = - x1mx4i * x2mx4i * x3mx4i

!     . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 

        end do  !  do jj

        do j = 1,platd-1
          dphi(j) = phi(j+1) - phi(j)
        end do

!my nlonex is extended lonsperlat and reversed so south to north sn
!my real lats

        do j=1,plat
           jj = platd -2 -j+1 
           nlonex(jj) = lonsperlat(j)
        enddo
!my sp
        nlonex(1) = lonsperlat(plat)
        nlonex(2) = lonsperlat(plat)
!my np
        nlonex(plat+3) = lonsperlat(1)
        nlonex(plat+4) = lonsperlat(1)

        lam0 = 0.0
        do j=1,platd
           dlam(j)   = twopi/float(nlonex(j))
           rdlam(j)  = 1.0/dlam(j)
           rdlam6(j) = (1.0/6.0)*rdlam(j)
           do i = 1,nlonex(j)+3
             lam(i,j) = float(i-2)*dlam(j) + lam0
           end do
        enddo

         if (me == 0) print*,' fini    grdini_h from slgscan_h '

!        call eta_stuff(lats_dim_a,lats_node_a,etamid,etaint,
!    .     kdpmpf,kdpmph,detam,detai,lbasdz,lbassd,lbasiz,lbassi)

!     inline  eta_stuff

        do k = 2,plev-2

!          call lcbas_h(  etamid(k-1), lbasiz(1,1,k),
!    .                               lbasiz(1,2,k) )
!          call lcdbas_h( etamid(k-1), lbasdz(1,1,k),
!    .                               lbasdz(1,2,k) )

!     . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 

!     inline  lcbas_h( grd, bas1, bas2 )

           k1 = k - 1

           x_1 = etamid(k1+0)
           x_2 = etamid(k1+1)
           x_3 = etamid(k1+2)
           x_4 = etamid(k1+3)

           x1mx2 = x_1 - x_2
           x1mx3 = x_1 - x_3
           x1mx4 = x_1 - x_4
           x2mx3 = x_2 - x_3
           x2mx4 = x_2 - x_4
           x3mx4 = x_3 - x_4

           x1mx2i = 1.0 / x1mx2
           x1mx3i = 1.0 / x1mx3
           x1mx4i = 1.0 / x1mx4
           x2mx3i = 1.0 / x2mx3
           x2mx4i = 1.0 / x2mx4
           x3mx4i = 1.0 / x3mx4

           lbasdz(1,1,k) =   x2mx3 * x2mx4 * x1mx2i * x1mx3i * x1mx4i
           lbasdz(2,1,k) = - x1mx2i + x2mx3i + x2mx4i
           lbasdz(3,1,k) = - x1mx2 * x2mx4 * x1mx3i * x2mx3i * x3mx4i
           lbasdz(4,1,k) =   x1mx2 * x2mx3 * x1mx4i * x2mx4i * x3mx4i

           lbasdz(1,2,k) = - x2mx3 * x3mx4 * x1mx2i * x1mx3i * x1mx4i
           lbasdz(2,2,k) =   x1mx3 * x3mx4 * x1mx2i * x2mx3i * x2mx4i
           lbasdz(3,2,k) = - x1mx3i - x2mx3i + x3mx4i
           lbasdz(4,2,k) = - x1mx3 * x2mx3 * x1mx4i * x2mx4i * x3mx4i

           lbasiz(1,1,k) = x_1
           lbasiz(2,1,k) = x_2
           lbasiz(3,1,k) = x_3
           lbasiz(4,1,k) = x_4

           lbasiz(1,2,k) =   x1mx2i * x1mx3i * x1mx4i
           lbasiz(2,2,k) = - x1mx2i * x2mx3i * x2mx4i
           lbasiz(3,2,k) =   x1mx3i * x2mx3i * x3mx4i
           lbasiz(4,2,k) = - x1mx4i * x2mx4i * x3mx4i

!     . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 

!     inline  lcdbas_h(grd,dbas2,dbas3)

!     . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 

        end do  !  do k
!
        mindetam = 999999.9999
        do k = 1,plev-1
          detam  (k) = etamid(k+1) - etamid(k)       !min(detam) < pmap
          if( detam(k) < mindetam ) mindetam = detam(k)
        end do
!
        k = plev + 1
        leps = 1.e-05
        pmap = (etamid(plev) - etamid(1)) / mindetam + 1
        del  = (etamid(plev) - etamid(1)) / float(pmap)
        rdel = float(pmap)/(etamid(levs) - etamid(1))

        if (me == 0) then
          print*,'calculated pmap value is = ',pmap
          print*,'mindetam,del=',mindetam,del,' plev=',plev
          print *,' etamid=',etamid(1:plev)
        endif
!
        kdpmpf(1) = 1
        k = 2
         do kk = 2,pmap
            dp = etamid  (1) + float(kk-1)*del
            if (dp > etamid(k)+leps) then
               kdpmpf(k) = kk
               k = k + 1
            endif
         enddo
      endif  ! fin ini_slg

!     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      if ( ini_dp .eq. 1 ) then

        do lan_loc=1,lats_node_a   !sela begin lan_loc over grids without halos
          lat_loc  = latg + 1
     &            - global_lats_a(ipt_lats_node_a+lats_node_a-lan_loc)
          dlam_con = dlam(2+lat_loc)
          do k=1,levs
            do i=1,lonsperlat(lat_loc)
               lammp_a(i,k,lan_loc) = float(i-1)*dlam_con
               phimp_a(i,k,lan_loc) = phi(2+lat_loc)
               sigmp_a(i,k,lan_loc) = etamid(k)
            enddo
          enddo
          do i=1,lonsperlat(lat_loc)
            sigmp_a(i,plev,lan_loc) = sigmp_a(i,plev,lan_loc) - 1.e-12
          enddo
        enddo

      endif ! fin ini_dp

      if ( ini_slg == 1 .or. ini_dp == 1) return

!sela  completed initialization of slg array and base functions

!    -------------------------------------------------------------------

      jcen = 2 + lat
      her_x = herm_x
      her_y = herm_y
      her_z = herm_z


!     if(lat.eq.1 .and. me.eq.0)
!     if(i_count.eq.2 .and. me.eq.0)
!    &  print*,'her_x=',her_x,' her_y=',her_y,' her_z=',her_z

      lam0 = 0.0

      do i=i_1,i_2
         cos_l_a(i) = cos(lam(i+1,jcen))
         sin_l_a(i) = sin(lam(i+1,jcen))
      enddo
!
!$$$      do k=1,levs
!$$$         do i=1,lonf
!$$$            lam_all_ar(i,k,lat) = lam(i+1)
!$$$            phi_all_ar(i,k,lat) = phi(2+lat)
!$$$         enddo
!$$$      enddo

       cos_f_a = cos(phi(jcen))
       sin_f_a = sin(phi(jcen))

      if(abs(phi(jcen)) >= phigs) then
        i_branch = 0

!$$$       print*,' phi(jcen)=',phi(jcen)*57.295,' jcen=',jcen,
!$$$     $ ' dep3dg i_branch=',i_branch

!$$$      if (me .eq. 0 .and. lan. eq. 1) then
!$$$         do k=1,plev
!$$$            write(133,*) ' k = ',k
!$$$            write(133,1340) phimp_a(:,k,1)
!$$$ 1340       format(4g14.8)
!$$$         enddo
!$$$         close(133)
!$$$      endif
!
         call dep3dg_h(i_1, i_2, s_lat_in_pe, n_lat_in_pe,
     &                 jcen, ztodt, iter, ra,
     &                 uu_gr_h_2(1,s_lat_in_pe),
     &                 vv_gr_h_2(1,s_lat_in_pe),
     &                 ww_gr_h_2(1,s_lat_in_pe),
     &                 lam, phi, dphi, etamid, detam, kdpmpf,
     &                 lammp_a(1,1,lan), phimp_a(1,1,lan),
     &                 sigmp_a(1,1,lan),
     &                 lamdp, phidp, sigdp, me,
     &                 sin_l_d, cos_l_d, sin_f_d, cos_f_d, rdel, pi)

!$$$         print*,' returned from dep3dg_h'
!$$$       do k=1,levs
!$$$         print*,'lat,i_branch,k,min,max of phidp(i_1,i_2,k) = ',
!$$$     . lat,i_branch,k,minval(phidp(i_1:i_2,k)),maxval(phidp(i_1:i_2,k))
!$$$       enddo

      else
        i_branch = 1

!$$$       print*,' phi(jcen)=',phi(jcen)*57.295,' jcen=',jcen,
!$$$     $ ' dep3ds i_branch=',i_branch

         call dep3ds_h(i_1, i_2, s_lat_in_pe, n_lat_in_pe,
     &                 jcen, ztodt, iter, ra,
     &                 uu_gr_h_2(1,s_lat_in_pe),
     &                 vv_gr_h_2(1,s_lat_in_pe),
     &                 ww_gr_h_2(1,s_lat_in_pe),
     &                 lam, phi, dphi, etamid, detam, kdpmpf,
     &                 lammp_a(1,1,lan), phimp_a(1,1,lan),
     &                 sigmp_a(1,1,lan),
     &                 lamdp, phidp, sigdp, me,
     &                 sin_l_d, cos_l_d, sin_f_d, cos_f_d, rdel, pi)

!$$$         print*,' returned from dep3ds_h'

      endif

!$$$      lam_all_dp(:,:,lat) = lamdp ! for plotting only
!$$$      phi_all_dp(:,:,lat) = phidp ! for plotting only

      dphibr = 1./( phi(platd/2+1) - phi(platd/2) )
      phibs  = phi(1)

!my compute jdp first before idp since rdlam is now lat dependent 

      do k = 1,levs
         do i=i_1,i_2
            jdp(i,k) = int ( (phidp(i,k) - phibs)*dphibr + 1. )
            if( phidp(i,k) >= phi(jdp(i,k)+1) )
     &          jdp(i,k) = jdp(i,k) + 1
         enddo
!$$$         print*,'lat,i_branch,k,min,max of jdp(i_1,i_2,k) = ',
!$$$     .    lat,i_branch,k,minval(jdp(i_1:i_2,k)),maxval(jdp(i_1:i_2,k))

      enddo
!
!sela do k = 1,levs
!sela    do i=i_1,i_2
!sela    if (me .eq. 46) write(4646,*) 'i,k,lamdp,jdp,rdlam(jdp-1) = ',
!sela.    i,k,lamdp(i,k),jdp(i,k)
!sela    enddo
!sela enddo
!sela if (me .eq. 46) close(4646)           

      do k = 1,levs
        do i=i_1,i_2
          idp(i,k,1) = 2 + int( lamdp(i,k) * rdlam(jdp(i,k)-1) )
          if (redgg_a) then
            idp(i,k,2) = 2 + int( lamdp(i,k) * rdlam(jdp(i,k)) )
            idp(i,k,3) = 2 + int( lamdp(i,k) * rdlam(jdp(i,k)+1) )
            idp(i,k,4) = 2 + int( lamdp(i,k) * rdlam(jdp(i,k)+2) )
          else
             idp(i,k,2) = idp(i,k,1)
             idp(i,k,3) = idp(i,k,1)
             idp(i,k,4) = idp(i,k,1)
          endif
!
          kdp(i,k) = bisection(kdpmpf,plev-1,
     &                 int((sigdp(i,k) - etamid(1))*rdel + 1. ))
          if(sigdp(i,k) >= etamid(kdp(i,k)+1)) then
            kdp(i,k) = kdp(i,k) + 1
          endif
!
          zdp_grid(i,k,lat) = (jdp(i,k) - jcen)
        enddo
      enddo
       
      call herxinit_h(i_1,i_2,idp,jdp,lamdp,xl,xr,hl,hr,dhl,dhr)
!
      call lagxinit_h(i_1,i_2,lam,lamdp,idp,jdp,x2,x3,
     &                    term1x,term2x,term3x,term4x)

       call heryinit_h(i_1,i_2,phi,dphi,phidp,jdp,ys,yn,
     &                 hs,hn,dhs,dhn,rdphi)

       call lagyinit_h(i_1,i_2,phi,dphi,lbasiy,phidp,jdp,
     &                 yb,yt,term1y,term2y,term3y,term4y)

      kdim   = plev
      kdimm1 = kdim - 1
      kdimm2 = kdim - 2
      do k = 1,plev
        do i = i_1,i_2
          kkdp(i,k) = min0( kdimm2, max0(2, kdp(i,k)) )
!sela     write(me+2000,123)i,k,jcen,kkdp(i,k)
        enddo
      enddo

!sela123   format(' i=',i6,' k=',i6,' jcen=',i6,' kkdp(i,k)=',i5)
!
      call herzinit_h(i_1,i_2,kdim,etamid,detam,sigdp,kdp,
     &                hb,ht,dhb,dht,zb,zt)
      call lagzinit_h_m(i_1,i_2,
     &                kdim,lbasiz,etamid,detam,sigdp,kdp,kkdp,
     &                hb,ht,dhb,dht,term1z,term2z,term3z,term4z)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!-----------------------------------------------------------------------
       do nt=1,ntrac
         if (her_x .and. her_y .and. her_z) then
           call her_xyz_in_h(i_1,i_2,s_lat_in_pe,n_lat_in_pe,
     &         idp,jdp,kdp,kkdp,xl,xr,hl,hr,dhl,dhr,
     &         rgt_h(1,1,s_lat_in_pe,nt),
     &         fintx,ys,yn,hs,hn,dhs,dhn,rdphi,finty,
     &         kdim,lbasdz,hb,ht,dhb,dht,detam,rgt_a(1,1,lan,nt))
         else
           if (her_x) then
!!!        call mpi_barrier (mc_comp,ierr)
!sela      print*,' entering herxin_h'
             call herxin_h(i_1,i_2,s_lat_in_pe,n_lat_in_pe,idp,jdp,kdp,
     &            kkdp,xl,xr,hl,hr,dhl,dhr,rgt_h(1,1,s_lat_in_pe,nt),
     &            fintx)
           else
             call lagxin_h(i_1,i_2,s_lat_in_pe,n_lat_in_pe,
     &                 kdim,rgt_h(1,1,s_lat_in_pe,nt),idp,jdp,kdp,kkdp,
     &                 x2,x3,term1x,term2x,term3x,term4x,fintx)
           endif
           if (her_y) then
             call heryin_h(i_1,i_2,jdp,ys,yn,hs,hn,dhs,dhn,rdphi,
     &                   fintx,finty)
           else
             call lagyin_h(i_1,i_2,fintx,finty,yb,yt,term1y,term2y,
     &                      term3y,term4y)
           endif

           if (her_z) then
             call herzin_h(i_1,i_2,
     &                   kdim,finty,lbasdz,kdp,hb,ht,dhb,dht,detam,
     &                   term1z,term2z,term3z,term4z,rgt_a(1,1,lan,nt))
           else
             call lagzin_h_m(i_1,i_2,
     &                   kdim,finty,lbasdz,kdp,hb,ht,dhb,dht,detam,
     &                   term1z,term2z,term3z,term4z,rgt_a(1,1,lan,nt))
           endif
         endif  !  if (her_x .and. her_y .and. her_z) then

!my weighted cubic linear interpolation for generalized number of tracers
        if (wgt_cub_lin_xyz) then
          call wgt_cub_lin_xyz_in_h(i_1,i_2,s_lat_in_pe,n_lat_ in_pe,
     &       rgt_h(1,1,s_lat_in_pe,nt),idp,jdp,kdp,xl,xr,ys,yn ,zb,zt,
     &       rgt_a(1,1,lan,nt))
        endif

       enddo
!-----------------------------------------------------------------------

c$$$       her_x = .true.                 !!!!!! switch interpolation
c$$$       her_y = .true.                 !!!!!! switch interpolation
!-----------------------------------------------------------------------
       if (her_x .and. her_y .and. her_z) then
         call her_xyz_in_h(i_1,i_2,s_lat_in_pe,n_lat_in_pe,
     &                     idp,jdp,kdp,kkdp,xl,xr,hl,hr,dhl,dhr,ud3_h1,
     &                     fintx,ys,yn,hs,hn,dhs,dhn,rdphi,finty,
     &                     kdim,lbasdz,hb,ht,dhb,dht,detam,ud3_d)

         call her_xyz_in_h(i_1,i_2,s_lat_in_pe,n_lat_in_pe,
     &                     idp,jdp,kdp,kkdp,xl,xr,hl,hr,dhl,dhr,vd3_h1,
     &                     fintx,ys,yn,hs,hn,dhs,dhn,rdphi,finty,
     &                     kdim,lbasdz,hb,ht,dhb,dht,detam,vd3_d)
       else
!                    For ud3
!                    -------
         if (her_x) then
           call herxin_h(i_1,i_2,s_lat_in_pe,n_lat_in_pe,idp,jdp,kdp,
     &                   kkdp,xl,xr,hl,hr,dhl,dhr,ud3_h1,fintx)
         else
           call lagxin_h(i_1,i_2,s_lat_in_pe,n_lat_in_pe,
     &                   kdim,ud3_h1,idp,jdp,kdp,kkdp,
     &                   x2,x3,term1x,term2x,term3x,term4x,fintx)
         endif

         if (her_y) then
           call heryin_h(i_1,i_2,jdp,ys,yn,hs,hn,dhs,dhn,rdphi,
     &                   fintx,finty)
         else
           call lagyin_h(i_1,i_2,fintx,finty,yb,yt,term1y,term2y,
     &                   term3y,term4y)
         endif
         if (her_z) then
             call herzin_h(i_1,i_2,
     &                    kdim,finty,lbasdz,kdp,hb,ht,dhb,dht,detam,
     &                    term1z,term2z,term3z,term4z,ud3_d)
         else
             call lagzin_h_m(i_1,i_2,
     &                    kdim,finty,lbasdz,kdp,hb,ht,dhb,dht,detam,
     &                    term1z,term2z,term3z,term4z,ud3_d)
          endif

!                    For vd3
!                    -------
          if (her_x) then
            call herxin_h(i_1,i_2,s_lat_in_pe,n_lat_in_pe,idp,jdp,kdp,
     &                    kkdp,xl,xr,hl,hr,dhl,dhr,vd3_h1,fintx)
          else
            call lagxin_h(i_1,i_2,s_lat_in_pe,n_lat_in_pe,
     &                    kdim,vd3_h1,idp,jdp,kdp,kkdp,
     &                    x2,x3,term1x,term2x,term3x,term4x,fintx)
          endif
          if (her_y) then
            call heryin_h(i_1,i_2,jdp,ys,yn,hs,hn,dhs,dhn,rdphi,
     &                    fintx,finty)
          else

            call lagyin_h(i_1,i_2,fintx,finty,yb,yt,term1y,term2y,
     &                    term3y,term4y)
          endif
          if (her_z) then
             call herzin_h(i_1,i_2,
     &                    kdim,finty,lbasdz,kdp,hb,ht,dhb,dht,detam,
     &                    term1z,term2z,term3z,term4z,vd3_d)
         else
             call lagzin_h_m(i_1,i_2,
     &                       kdim,finty,lbasdz,kdp,hb,ht,dhb,dht,detam,
     &                       term1z,term2z,term3z,term4z,vd3_d)
         endif
       endif  !  if (her_x.and.her_y.and.her_z) then
c
       if (her_x.and.her_y.and.her_z) then
          call her_xyz_in_h(i_1,i_2,s_lat_in_pe,n_lat_in_pe,
     &         idp,jdp,kdp,kkdp,xl,xr,hl,hr,dhl,dhr,td3_h1,
     &         fintx,ys,yn,hs,hn,dhs,dhn,rdphi,finty,
     &         kdim,lbasdz,hb,ht,dhb,dht,detam,td3_a(1,1,lan))
       else

         if (her_x) then
           call herxin_h(i_1,i_2,s_lat_in_pe,n_lat_in_pe,idp,jdp,kdp,
     &                   kkdp,xl,xr,hl,hr,dhl,dhr,td3_h1,fintx)
         else
           call lagxin_h(i_1,i_2,s_lat_in_pe,n_lat_in_pe,
     &                   kdim,td3_h1,idp,jdp,kdp,kkdp,
     &                   x2,x3,term1x,term2x,term3x,term4x,fintx)

         endif
         if (her_y) then
           call heryin_h(i_1,i_2,jdp,ys,yn,hs,hn,dhs,dhn,rdphi,
     &                   fintx,finty)
         else

           call lagyin_h(i_1,i_2,fintx,finty,yb,yt,term1y,term2y,
     &                   term3y,term4y)
         endif
         if (her_z) then
          call herzin_h(i_1,i_2,
     &                  kdim,finty,lbasdz,kdp,hb,ht,dhb,dht,detam,
     &                  term1z,term2z,term3z,term4z,td3_a(1,1,lan))
         else
           call lagzin_h_m(i_1,i_2,
     &                     kdim,finty,lbasdz,kdp,hb,ht,dhb,dht,detam,
     &                     term1z,term2z,term3z,term4z,td3_a(1,1,lan))
         endif
       endif  !  if (her_x.and.her_y.and.her_z) then


cmy weighted cubic linear interpolation
       if (wgt_cub_lin_xyz) then
          call wgt_cub_lin_xyz_in_h(i_1,i_2,s_lat_in_pe,n_lat_in_pe,
     &       ud3_h1,idp,jdp,kdp,xl,xr,ys,yn,zb,zt,
     &       ud3_d)

          call wgt_cub_lin_xyz_in_h(i_1,i_2,s_lat_in_pe,n_lat_in_pe,
     &       vd3_h1,idp,jdp,kdp,xl,xr,ys,yn,zb,zt,
     &       vd3_d)

          call wgt_cub_lin_xyz_in_h(i_1,i_2,s_lat_in_pe,n_lat_in_pe,
     &       td3_h1,idp,jdp,kdp,xl,xr,ys,yn,zb,zt,
     &       td3_a(1,1,lan))

       endif

!my use linear interpolation for u,v,t for nonlinear terms and add to advected terms ug,vg,tg

       if (lin_xyz) then
          call lin_xyz_in_h(i_1,i_2,s_lat_in_pe,n_lat_in_pe,
     &       ud3_h2,vd3_h2,td3_h2,idp,jdp,kdp,xl,xr,ys,yn,zb,zt,
     &       ud3_d,vd3_d,td3_a(1,1,lan))
       endif


!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!sela udmax=-9999.
!sela udmin= 9999.
!sela do k=1,levs
!sela  do i=1,lonf
!sela   if(ud3_d(i,k) > udmax)udmax=ud3_d(i,k)
!sela   if(ud3_d(i,k) < udmin)udmin=ud3_d(i,k)
!sela  enddo
!sela enddo
!sela print*,'lan=',lan,' lonf=',lonf,' udmax udmin=',udmax,udmin
!sela udmax=-9999.
!sela udmin= 9999.
!sela do k=1,levs
!sela  do i=1,lonf
!sela   if(vd3_d(i,k) > udmax)udmax=vd3_d(i,k)
!sela   if(vd3_d(i,k) < udmin)udmin=vd3_d(i,k)
!sela  enddo
!sela enddo
!sela print*,'lan=',lan,' lonf=',lonf,' vdmax vdmin=',udmax,udmin

      call rotate_uv_h(i_1,i_2,lon_dim_h,
     &                 ud3_d,vd3_d, ud3_a(1,1,lan),vd3_a(1,1,lan),
     &                 cos_l_a,sin_l_a,cos_f_a,sin_f_a,
     &                 cos_l_d,sin_l_d,cos_f_d,sin_f_d)

!sela udmax=-9999.
!sela udmin= 9999.
!sela do k=1,levs
!sela  do i=1,lonf
!sela   if(ud3_a(i,k,lan) > udmax)udmax=ud3_a(i,k,lan)
!sela   if(ud3_a(i,k,lan) < udmin)udmin=ud3_a(i,k,lan)
!sela  enddo
!sela enddo
!sela print*,'lan=',lan,' lonf=',lonf,' uamax uamin=',udmax,udmin
!sela udmax=-9999.
!sela udmin= 9999.
!sela do k=1,levs
!sela  do i=1,lonf
!sela   if(vd3_a(i,k,lan) > udmax)udmax=vd3_a(i,k,lan)
!sela   if(vd3_a(i,k,lan) < udmin)udmin=vd3_a(i,k,lan)
!sela  enddo
!sela enddo
!sela print*,'lan=',lan,' lonf=',lonf,' vamax vamin=',udmax,udmin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         call int2d_h_levs(i_1,i_2,s_lat_in_pe,n_lat_in_pe,kdim,
     &              dpdt_d3_h,idp,jdp,
     &         x2,x3,term1x,term2x,term3x,term4x,
     &               term1y,term2y,term3y,term4y,
     &              dpdt_d3_a(1,1,lan))

      return
      end
      subroutine dep3dg_h(i_1,i_2,s_lat_in_pe,n_lat_in_pe,
     &                  jcen,dt,iterdp,ra,ub,vb,wb,
     &                  lam,phib,dphib,etamid,detam,
     &                  kdpmpf,lammp_a,phimp_a,sigmp_a,
     &                  lamdp,phidp,sigdp,me,
     &                  sin_l_d,cos_l_d,sin_f_d,cos_f_d,rdel,pi)
!
      use      pmgrid        , only : mprec,plat,platd,
     &                                plev,plon,plond
      use slgshr             , only : dlam,rdlam
      use namelist_def       , only : redgg_a
!
      implicit none
      integer i_1,i_2,s_lat_in_pe,n_lat_in_pe
      integer jcen,iterdp,idp(i_1:i_2,plev,4), jdp(i_1:i_2,plev),
     $                    kdp(i_1:i_2,plev),me
!sela$                    kdp(i_1:i_2,plev),kdph(i_1:i_2,plev),me
      integer kdpmpf(plev-1)   ! artificial vert grid indic
      real ub (plond,plev ,s_lat_in_pe:n_lat_in_pe),
     &     vb (plond,plev ,s_lat_in_pe:n_lat_in_pe),
     &     wb (plond,plev ,s_lat_in_pe:n_lat_in_pe),
     &     uxl(i_1:i_2,plev ,2,platd),
     &     uxr(i_1:i_2,plev ,2,platd)
      real detam(plev)
      real etamid(plev)           ! eta at levels
      real dt,ra,lam(plond,platd),phib(platd),dphib(platd),
     &     lammp_a(plon,plev),phimp_a(plon,plev),sigmp_a(plon,plev),
     &     lamdp(i_1:i_2,plev),phidp(i_1:i_2,plev),sigdp(i_1:i_2,plev)

      real    sin_l_d(i_1:i_2,plev)   ,cos_l_d(i_1:i_2,plev),
     &        sin_f_d(i_1:i_2,plev)   ,cos_f_d(i_1:i_2,plev)

      real xl(i_1:i_2,plev,4),xr(i_1:i_2,plev,4),
     &     ys(i_1:i_2,plev),yn(i_1:i_2,plev)
      integer iter,i,k,twofld
      real junk1_address,junk2_address
      real ump(i_1:i_2,plev),vmp(i_1:i_2,plev),wmp(i_1:i_2,plev),
     &     upr(i_1:i_2,plev),vpr(i_1:i_2,plev),
     &     phicen,cphic,sphic,cphic2i,
     &     lampr(i_1:i_2,plev),phipr(i_1:i_2,plev),
     &     sigpr(i_1:i_2,plev),hfdt
      real sphisc,cphisc,clamsc,phibs,dphibr
      real fac
      parameter ( fac = 1. - 1.e-12 )
      real pi,twopi,pi2,sgnphi,sphipr,cphipr,clampr,slam2,phipi2,
     &     slampr(i_1:i_2,plev),dlamx(i_1:i_2,plev),coeff,distmx,
     &       dist(i_1:i_2,plev)
      real cdlam,clamp,cphimp,cphip,sdlam,slamp,sphimp,sphip
      real rdel
      integer kk
      integer bisection
!sela print*,'i_1,i_2 in dep3dg=',i_1,i_2
!sela print*,'lan,lats_node_a in dep3dg=',lan,lats_node_a
!sela print*,'jcen,dt,iterdp,ra in dep3dg=',jcen,dt,iterdp,ra

      twopi = pi  + pi
      pi2   = 0.5 * pi
      hfdt  = 0.5 * dt
      phicen = phib(jcen)
      cphic  = cos( phicen )
      sphic  = sin( phicen )
!
      do k = 1,plev
         do i = i_1,i_2
            sphisc = sin( phimp_a(i,k) )
            cphisc = cos( phimp_a(i,k) )
            clamsc = cos( lam(i+1,jcen) - lammp_a(i,k) )
            phipr(i,k) = asin( sphisc*cphic - cphisc*sphic*clamsc )
         end do
      end do
!

!$$$      if (me .eq. 0)  then
!$$$         do k=1,plev
!$$$            write(134,*) ' k = ',k
!$$$            write(134,1340) phimp_a(:,k)
!$$$ 1340       format(4g14.8)
!$$$         enddo
!$$$         close(134)
!$$$      endif
				   
      dphibr = 1./( phib(platd/2+1) - phib(platd/2) )
      phibs  = phib(1)

      do iter=1,iterdp


!my compute jdp first since rdlam is now lat dependent 


        do k = 1,plev

!$$$      print*,
!$$$     . 'bef97 dep3dg jcen,k,mi/ma phimp_a(i,k),phibs,dphibr,mi/ma jd=',
!$$$     .   jcen,k,minval(phimp_a(i_1:i_2,k)),maxval(phimp_a(i_1:i_2,k)),
!$$$     .   phibs,dphibr,minval(jdp(:,k)),maxval(jdp(:,k))

          do i=i_1,i_2
            jdp(i,k) = int ( (phimp_a(i,k) - phibs)*dphibr + 1. )

! Commented by Moorthi on 10/18/2011
!        if (phimp_a(i,k) < -900.0) then
!        write(999,*) 
!    .   'inside dep3dg i,k,phimp_a(i,k),phibs,dphibr,jdp = ',
!    .      i,k,phimp_a(i,k),phibs,dphibr,jdp(i,k)
!         close(999)
!         stop 999
!        endif

            if( phimp_a(i,k) >= phib(jdp(i,k)+1) )
     &          jdp(i,k) = jdp(i,k) + 1


!$$$      print*,
!$$$     . 'in dep3dg jcen,k,mi/ma phimp_a(i,k),phibs,dphibr,mi/ma jdp = ',
!$$$     .   jcen,k,minval(phimp_a(i_1:i_2,k)),maxval(phimp_a(i_1:i_2,k)),
!$$$     .   phibs,dphibr,minval(jdp(:,k)),maxval(jdp(:,k))


!111   format('jcen=',i5,' i=',i4,2x,' k=',i3,2x,' jdp-jcen=',i6)
!
            idp(i,k,1) = 2 + int( lammp_a(i,k) * rdlam(jdp(i,k)-1) )
            if (redgg_a) then
               idp(i,k,2) = 2 + int( lammp_a(i,k) * rdlam(jdp(i,k)) )
               idp(i,k,3) = 2 + int( lammp_a(i,k) * rdlam(jdp(i,k)+1) )
               idp(i,k,4) = 2 + int( lammp_a(i,k) * rdlam(jdp(i,k)+2) )
            else
               idp(i,k,2) = idp(i,k,1)
               idp(i,k,3) = idp(i,k,1)
               idp(i,k,4) = idp(i,k,1)
            endif
!
            kdp(i,k) = bisection(kdpmpf,plev-1,
     &                 int((sigmp_a(i,k) - etamid(1))*rdel + 1. ))
            if(sigmp_a(i,k) >= etamid(kdp(i,k)+1)) then
              kdp(i,k) = kdp(i,k) + 1
            end if
          enddo
        enddo
!
        call xywgts_h(i_1,i_2,
     &              lam,phib,dphib,idp,jdp,lammp_a,phimp_a,xl,xr,ys,yn)
!
        call int3dv_h(i_1,i_2,s_lat_in_pe,n_lat_in_pe,
     &                ub,vb,etamid,detam,sigmp_a,idp,jdp,kdp,
     &                xl,xr,ys,yn,ump,vmp,
     &                wb,wmp)
!sela&              wb,wmp,etaint,detai,kdph)
!
        do k = 1,plev
          do i = i_1,i_2
            ump(i,k) = ump(i,k)*ra
            vmp(i,k) = vmp(i,k)*ra
!
            sdlam  = sin( lam(i+1,jcen) - lammp_a(i,k) )
            cdlam  = cos( lam(i+1,jcen) - lammp_a(i,k) )
            sphimp = sin( phimp_a(i,k) )
            cphimp = cos( phimp_a(i,k) )
            sphip  = sphimp*cphic - cphimp*sphic*cdlam
            cphip  = cos( asin( sphip ) )
            slamp  = -sdlam*cphimp/cphip
            clamp  = cos( asin( slamp ) )
            vpr(i,k) = (vmp(i,k)*(cphimp*cphic + sphimp*sphic*cdlam) -
     &                  ump(i,k)*sphic*sdlam)/cphip
            upr(i,k) = (ump(i,k)*cdlam + vmp(i,k)*sphimp*sdlam +
     &                  vpr(i,k)*slamp*sphip)/clamp
!
            lampr(i,k) = -hfdt* upr(i,k) / cos( phipr(i,k) )
            phipr(i,k) = -hfdt* vpr(i,k)
          enddo
        enddo

        coeff  = (1.1*hfdt)**2
        distmx = (sign(pi2,phicen) - phicen)**2/coeff
        sgnphi = sign( 1., phicen )
!
        do k=1,plev
          do i=i_1,i_2
            sphipr      = sin( phipr(i,k) )
            cphipr      = cos( phipr(i,k) )
            slampr(i,k) = sin( lampr(i,k) )
            clampr      = cos( lampr(i,k) )
            phimp_a(i,k)=asin((sphipr*cphic + cphipr*sphic*clampr)*fac)
            if ( abs(phimp_a(i,k)) >= phib(3+plat)*fac )
     $      phimp_a(i,k) = sign( phib(3+plat),phimp_a(i,k) )*fac
            dlamx(i,k)=asin((slampr(i,k)*cphipr/cos(phimp_a(i,k)))*fac)
            dist(i,k) = upr(i,k)*upr(i,k) + vpr(i,k)*vpr(i,k)
!
            if (dist(i,k) > distmx) then
              slam2 = slampr(i,k)*slampr(i,k)
              phipi2 = asin((sqrt((slam2-1.)/(slam2-cphic2i)))*fac)
              if (sgnphi*phipr(i,k) > phipi2) then
                dlamx(i,k) = sign(pi,lampr(i,k)) - dlamx(i,k)
              end if
            end if
!
            lammp_a(i,k) = lam(i+1,jcen) + dlamx(i,k)
            if( lammp_a(i,k) >= twopi ) lammp_a(i,k)=lammp_a(i,k)-twopi
            if( lammp_a(i,k) <  0.0 )   lammp_a(i,k)=lammp_a(i,k)+twopi
!
               sigpr(i,k) =  -hfdt*wmp(i,k)
               sigmp_a(i,k) = etamid(k) + sigpr(i,k)
!
            if (sigmp_a(i,k) < etamid(1)) then
                sigmp_a(i,k) = etamid(1)
            end if
            if (sigmp_a(i,k) >= etamid(plev)) then
                sigmp_a(i,k) = etamid(plev) - mprec
            end if
          enddo
        enddo
      enddo          ! do iter=1,iterdp ends here
!
      do 110        k  = 1,plev
         do 109     i  = i_1,i_2
            lampr(i,k) = 2.*lampr(i,k)
            phipr(i,k) = 2.*phipr(i,k)
            sigdp(i,k) = etamid(k) + 2.*sigpr(i,k)
  109    continue
  110 continue
      coeff  = (1.1*dt)**2
      distmx = (sign(pi2,phicen) - phicen)**2/coeff
      sgnphi = sign( 1., phicen )
!
      do k=1,plev
         kk = plev+1-k
         do i=i_1,i_2
            sphipr      = sin( phipr(i,k) )
            cphipr      = cos( phipr(i,k) )
            slampr(i,k) = sin( lampr(i,k) )
            clampr      = cos( lampr(i,k) )
            phidp(i,k)=asin((sphipr*cphic + cphipr*sphic*clampr)*fac)
            if ( abs(phidp(i,k)) >= phib(3+plat)*fac )
     $           phidp(i,k) = sign( phib(3+plat),phidp(i,k) )*fac

            cos_f_d(i,kk)=cos(phidp(i,k)) ! for vector allignment
            sin_f_d(i,kk)=sin(phidp(i,k)) ! for vector allignment

            dlamx(i,k)=asin((slampr(i,k)*cphipr/cos(phidp(i,k)))*fac)
            dist(i,k) = upr(i,k)*upr(i,k) + vpr(i,k)*vpr(i,k)
!
            if (dist(i,k) > distmx) then
              slam2 = slampr(i,k)**2
              phipi2 = asin((sqrt((slam2-1.)/(slam2-cphic2i)))*fac)
              if (sgnphi*phipr(i,k) > phipi2) then
                dlamx(i,k) = sign(pi,lampr(i,k)) - dlamx(i,k)
              end if
            end if
!
            lamdp(i,k) = lam(i+1,jcen) + dlamx(i,k)
            if( lamdp(i,k) >= twopi ) lamdp(i,k) = lamdp(i,k) - twopi
            if( lamdp(i,k) < 0.0  )   lamdp(i,k) = lamdp(i,k) + twopi

            cos_l_d(i,kk)=cos(lamdp(i,k)) ! for vector allignment
            sin_l_d(i,kk)=sin(lamdp(i,k)) ! for vector allignment
!

           if (sigdp(i,k) < etamid(1)) sigdp(i,k) = etamid(1)
           if (sigdp(i,k) >= etamid(plev))
     &                                 sigdp(i,k) = etamid(plev) - mprec
        enddo
      enddo
      return
      end
      subroutine dep3ds_h(i_1,i_2,s_lat_in_pe,n_lat_in_pe,
     &                  jcen,dt,iterdp,ra,ub,vb,wb,
     &                  lam,phib,dphib,etamid,detam,
     &                  kdpmpf,lammp_a,phimp_a,sigmp_a,
     &                  lamdp,phidp,sigdp,me,
     &                  sin_l_d,cos_l_d,sin_f_d,cos_f_d,rdel,pi)
      use      pmgrid        , only : mprec,platd,plev,
     &                                plon,plond
      use slgshr             , only : dlam,rdlam
      use namelist_def       , only : redgg_a

      implicit none
      integer i_1,i_2,s_lat_in_pe,n_lat_in_pe
      integer jcen,iterdp,idp(i_1:i_2,plev,4), jdp(i_1:i_2,plev),
     &                    kdp(i_1:i_2,plev)                   ,me
!sela$                    kdp(i_1:i_2,plev),kdph(i_1:i_2,plev),me
      integer kdpmpf(plev-1)   ! artificial vert grid indic
      real ub (plond,plev ,s_lat_in_pe:n_lat_in_pe),
     &     vb (plond,plev ,s_lat_in_pe:n_lat_in_pe),
     &     wb (plond,plev ,s_lat_in_pe:n_lat_in_pe),
     &     uxl(i_1:i_2,plev ,2,platd),
     &     uxr(i_1:i_2,plev ,2,platd)
      real detam(plev)
      real etamid(plev)          ! eta at levels
      real dt,ra,lam(plond,platd),phib(platd),dphib(platd),
     &     lammp_a(plon,plev),phimp_a(plon,plev),sigmp_a(plon,plev),
     &     lamdp(i_1:i_2,plev),phidp(i_1:i_2,plev),sigdp(i_1:i_2,plev)

      real    sin_l_d(i_1:i_2,plev)   ,cos_l_d(i_1:i_2,plev),
     &        sin_f_d(i_1:i_2,plev)   ,cos_f_d(i_1:i_2,plev)
!
      real    xl(i_1:i_2,plev,4), xr(i_1:i_2,plev,4),
     &        ys(i_1:i_2,plev)  , yn(i_1:i_2,plev)
      integer iter,i,k,twofld
      real junk1_address,junk2_address
      real ump(i_1:i_2,plev),vmp(i_1:i_2,plev),  wmp(i_1:i_2,plev),
     &     phicen,lampr(i_1:i_2,plev),phipr(i_1:i_2,plev),
     &                     sigpr(i_1:i_2,plev),pi,twopi,hfdt
      real dphibr,phibs, rdel
      integer kk,bisection

!     pi     = 4.*atan(1.)
!     twopi  = 2.*pi
      twopi  = pi  + pi
      hfdt   = 0.5 * dt
      phicen = phib(jcen)
!
      dphibr = 1./( phib(platd/2+1) - phib(platd/2) )
      phibs  = phib(1)
!
      do iter=1,iterdp


!my compute jdp first since rdlam is now lat dependent 

        do k = 1,plev
          do i=i_1,i_2
            jdp(i,k) = int ( (phimp_a(i,k) - phibs)*dphibr + 1. )

! Commented by Moorthi on 10/18/2011
!        if (phimp_a(i,k) < -900.0) then
!        write(999,*) 
!    .   'inside dep3ds i,k,phimp_a(i,k),phibs,dphibr,jdp = ',
!    .      i,k,phimp_a(i,k),phibs,dphibr,jdp(i,k)
!         close(999)
!         stop 999
!        endif

!$$$      print*,'inside dep3ds jcen,i,k,phimp_a(i,k),phibs,dphibr,jdp = ',
!$$$     .      jcen,i,k,phimp_a(i,k),phibs,dphibr,jdp(i,k)

            if( phimp_a(i,k) >= phib(jdp(i,k)+1) )
     &          jdp(i,k) = jdp(i,k) + 1

!111   format('jcen=',i5,' i=',i4,2x,' k=',i3,2x,' jdp-jcen=',i6)

            idp(i,k,1) = 2 + int( lammp_a(i,k) * rdlam(jdp(i,k)-1) )
            if (redgg_a) then
               idp(i,k,2) = 2 + int( lammp_a(i,k) * rdlam(jdp(i,k)) )
               idp(i,k,3) = 2 + int( lammp_a(i,k) * rdlam(jdp(i,k)+1) )
               idp(i,k,4) = 2 + int( lammp_a(i,k) * rdlam(jdp(i,k)+2) )
            else
               idp(i,k,2) = idp(i,k,1)
               idp(i,k,3) = idp(i,k,1)
               idp(i,k,4) = idp(i,k,1)
            endif
!

            kdp(i,k) = bisection(kdpmpf,plev-1,
     &                 int((sigmp_a(i,k) - etamid(1))*rdel + 1. ))
            if(sigmp_a(i,k) >= etamid(kdp(i,k)+1)) then
              kdp(i,k) = kdp(i,k) + 1
            end if
          enddo
        enddo

        call xywgts_h(i_1,i_2,
     &              lam,phib,dphib,idp,jdp,lammp_a,phimp_a,xl,xr,ys,yn)

        call int3dv_h(i_1,i_2,s_lat_in_pe,n_lat_in_pe,
     &                ub,vb,etamid,detam,sigmp_a,idp,jdp,kdp,
     &                xl,xr,ys,yn,ump,vmp,
     &                wb,wmp)
!sela&                wb,wmp,etaint,detai,kdph)

        do k = 1,plev
          do i = i_1,i_2
            ump(i,k) = ump(i,k)*ra
            vmp(i,k) = vmp(i,k)*ra
!
            lampr(i,k) =  -hfdt*ump(i,k) / cos(phimp_a(i,k))
            phipr(i,k) =  -hfdt*vmp(i,k)
            sigpr(i,k) =  -hfdt*wmp(i,k)
            lammp_a(i,k) = lam(i+1,jcen) + lampr(i,k)
            phimp_a(i,k) = phicen        + phipr(i,k)
            sigmp_a(i,k) = etamid(k)     + sigpr(i,k)
            if(lammp_a(i,k) >= twopi) lammp_a(i,k)=lammp_a(i,k)-twopi
            if(lammp_a(i,k) <    0.0) lammp_a(i,k)=lammp_a(i,k)+twopi
!
            if (sigmp_a(i,k) < etamid(1)) sigmp_a(i,k) = etamid(1)
            if (sigmp_a(i,k) >= etamid(plev))
     &                                sigmp_a(i,k) = etamid(plev)-mprec
          enddo
        enddo
      enddo                                           !do iter=1,iterdp loop

      do k  = 1,plev
         kk = plev+1-k
         do i=i_1,i_2
            lamdp(i,k) = lam(i+1,jcen) + 2.*lampr(i,k)
            phidp(i,k) = phicen        + 2.*phipr(i,k)
            sigdp(i,k) = etamid(k)     + 2.*sigpr(i,k)

            if(lamdp(i,k) >= twopi) lamdp(i,k) = lamdp(i,k) - twopi
            if(lamdp(i,k) <   0.0)  lamdp(i,k) = lamdp(i,k) + twopi

            cos_l_d(i,kk) = cos(lamdp(i,k))           ! for vector allignment
            sin_l_d(i,kk) = sin(lamdp(i,k))           ! for vector allignment

            cos_f_d(i,kk) = cos(phidp(i,k))           ! for vector allignment
            sin_f_d(i,kk) = sin(phidp(i,k))           ! for vector allignment

            if (sigdp(i,k) < etamid(1)) sigdp(i,k) = etamid(1)
            if (sigdp(i,k) >= etamid(plev))
     &                                  sigdp(i,k) = etamid(plev)-mprec
        enddo
      enddo
      return
      end
      subroutine grdini_h(rcs2,lonsperlat)
      use      pmgrid        , only : plat,platd
      use slgshr             , only : dlam,dphi,lam,lbasdy,lbasiy,
     &                                nlonex,phi,rdlam,rdlam6
      implicit none
      real     rcs2  (plat/2)       !  1/cos**2 at gaussian colatitudes
      integer lonsperlat(plat)
      integer       imin,kk,newmap
      real del,dp,pi,degrad
      integer jj,i,j,k,ig
      real lam0
      external lcdbas_h,lcbas_h
!     print*,' reduced grid inside grdini_h  '
      pi = 4.*atan(1.)
      degrad = 180./pi

      do i=1,plat/2
        phi(i+2) =         -acos(sqrt(1./rcs2(i)))
        phi(platd -1-i) =  acos(sqrt(1./rcs2(i)))
      enddo

      phi(1) = -pi - phi(3)
      phi(2) = -pi/2.0

      phi(plat+3) =  pi/2.0
      phi(platd) = pi - phi(plat+2)
!     do j=1,platd
!       write(6,100)j,phi(j)*degrad
!     enddo
100   format(' in grdini    lat=',i3,2x,' phi=',f9.3)
      do jj = 2,plat+2
         call lcdbas_h( phi(jj-1), lbasdy(1,1,jj), lbasdy(1,2,jj) )
      end do
      do jj = 2,plat+2
         call lcbas_h( phi(jj-1), lbasiy(1,1,jj), lbasiy(1,2,jj) )
      end do
      do j = 1,platd-1
         dphi(j) = phi(j+1) - phi(j)
      end do
 
!
!my nlonex is extended lonsperlat and reversed so south to north sn
!
!my real lats
      do j=1,plat
         jj = platd -2 -j+1 
         nlonex(jj) = lonsperlat(j)
      enddo
!my sp
      nlonex(1) = lonsperlat(plat)
      nlonex(2) = lonsperlat(plat)
!my np
      nlonex(plat+3) = lonsperlat(1)
      nlonex(plat+4) = lonsperlat(1)

      lam0 = 0.0
      do j=1,platd
         dlam(j) = 2.*pi/float(nlonex(j))
         rdlam(j) = 1.0/dlam(j)
         rdlam6(j) = 1.0/(6.0*dlam(j))
         do i = 1,nlonex(j)+3
            lam(i,j)=float(i-2)*dlam(j) + lam0
         end do
      enddo

      return
      end
      subroutine herxinit_h(i_1,i_2,idp,jdp,lamdp,xl,xr,hl,hr,
     &                      dhl,dhr)
!
      use      pmgrid        , only : plev
      use slgshr             , only : lam,dlam,rdlam
      implicit none
!
      integer i_1,i_2
      integer idp(i_1:i_2,plev,4),jdp(i_1:i_2,plev)
      real    lamdp(i_1:i_2,plev     )       ! x-coord of dep pt

      real    xl(i_1:i_2,plev,4), xr(i_1:i_2,plev,4),
     &        hl(i_1:i_2,plev,4), hr(i_1:i_2,plev,4),
     &        dhl(i_1:i_2,plev,4),dhr(i_1:i_2,plev,4)

      integer jj(i_1:i_2,plev,4)

      integer i,k,n
      real    tem1, tem2
!
      do n=1,4
        do k=1,plev
          do  i=i_1,i_2
!$$$        jj(i,k,n) = jdp(i,k) - 2 + n
!$$$        xl(i,k,n) = (lam(idp(i,k,n)+1,jj(i,k,n)) - lamdp(i,k))*rdx(jj(i,k,n))
            jj(i,k,n) = jdp(i,k) - 2 + n
            xl(i,k,n) = (lam(idp(i,k,n)+1,jj(i,k,n)) - lamdp(i,k))*
     &                   rdlam(jj(i,k,n))
            xr(i,k,n) = 1. - xl(i,k,n)
          end do
        end do
      end do
      do n=2,3
        do k=1,plev
          do i=i_1,i_2
            tem1 = xl(i,k,n) * xl(i,k,n)
            tem2 = xr(i,k,n) * xr(i,k,n)
            hl (i,k,n) = ( 3.0 - 2.0*xl(i,k,n) ) * tem1
            hr (i,k,n) = ( 3.0 - 2.0*xr(i,k,n) ) * tem2

!$$$        dhl(i,k,n) = -dx(jj(i,k,n))*( xl(i,k,n) - 1. )*xl(i,k,n)**2
!$$$        dhr(i,k,n) =  dx(jj(i,k,n))*( xr(i,k,n) - 1. )*xr(i,k,n)**2

            dhl(i,k,n) =-dlam(jj(i,k,n)) * (xl(i,k,n)-1. ) * tem1
            dhr(i,k,n) = dlam(jj(i,k,n)) * (xr(i,k,n)-1. ) * tem2
          end do
        end do
      enddo

      return
      end
      subroutine int2d_h_levs(i_1,i_2,s_lat_in_pe,n_lat_in_pe,
     &                        kdim,fb1,idp,jdp,
     &                        x2,x3,term1x,term2x,term3x,term4x,
     &                        term1y,term2y,term3y,term4y, fdp1)
!
      use      pmgrid        , only : plev,plon,plond
      implicit none
!
      integer i_1,i_2,s_lat_in_pe,n_lat_in_pe
      integer kdim,idp(i_1:i_2,plev,4),jdp(i_1:i_2,plev)
      real fb1 (plond,kdim,s_lat_in_pe:n_lat_in_pe),
     *  fdp1(plon,plev)
      integer i,k, kk
      real f1,f2,f3,f4
      real term1x(i_1:i_2,plev,4),term2x(i_1:i_2,plev,4),
     $     term3x(i_1:i_2,plev,4),term4x(i_1:i_2,plev,4)
      real term1y(i_1:i_2,plev),term2y(i_1:i_2,plev),
     $     term3y(i_1:i_2,plev),term4y(i_1:i_2,plev)
      real x2(i_1:i_2,plev,4),x3(i_1:i_2,plev,4)
      integer ii1,ii2,ii3,ii4,jj
!
      do k=1,plev
        kk = plev+1-k
        do i=i_1,i_2
          ii1 = idp(i,k,1)
          ii2 = idp(i,k,2)
          ii3 = idp(i,k,3)
          ii4 = idp(i,k,4)
          jj = jdp(i,k)

          f1 = fb1 (ii1+1,k,jj-1) * x2    (i,k,1)
     &       - fb1 (ii1  ,k,jj-1) * x3    (i,k,1)
          f2 = fb1 (ii2-1,k,jj  ) * term1x(i,k,2)
     &       + fb1 (ii2  ,k,jj  ) * term2x(i,k,2)
     &       + fb1 (ii2+1,k,jj  ) * term3x(i,k,2)
     &       + fb1 (ii2+2,k,jj  ) * term4x(i,k,2)
          f3 = fb1 (ii3-1,k,jj+1) * term1x(i,k,3)
     &       + fb1 (ii3  ,k,jj+1) * term2x(i,k,3)
     &       + fb1 (ii3+1,k,jj+1) * term3x(i,k,3)
     &       + fb1 (ii3+2,k,jj+1) * term4x(i,k,3)
          f4 = fb1 (ii4+1,k,jj+2) * x2    (i,k,4)
     &       - fb1 (ii4  ,k,jj+2) * x3    (i,k,4)
!
          fdp1(i,kk) = f1*term1y(i,k) + f2*term2y(i,k)
     &               + f3*term3y(i,k) + f4*term4y(i,k)
        end do
      end do
      return
      end
      subroutine int3dv_h(i_1,i_2,s_lat_in_pe,n_lat_in_pe,
     &                    fb1,fb2,etamid,detam,zdp,idp,jdp,kdp,
     &                    xl,xr,ys,yn,fdp1,fdp2,
     &                    fb3,fdp3)
!sela&                    fb3,fdp3,etaint,detai,kdph)
!
      use      pmgrid        , only : plev,plon,plond
      implicit none
!
      real detam(plev)
      real etamid(plev)           ! eta at levels
      integer i_1,i_2,s_lat_in_pe,n_lat_in_pe
      integer idp(i_1:i_2,plev,4), jdp(i_1:i_2,plev),
     &        kdp(i_1:i_2,plev)
!sela&        kdp(i_1:i_2,plev),kdph(i_1:i_2,plev)
      real fb1(plond,plev ,s_lat_in_pe:n_lat_in_pe),
     &     fb2(plond,plev ,s_lat_in_pe:n_lat_in_pe),
     &     fb3(plond,plev ,s_lat_in_pe:n_lat_in_pe),
     &     zdp(plon,plev),
     &     fdp1(i_1:i_2,plev),
     &     fdp2(i_1:i_2,plev),
     &     fdp3(i_1:i_2,plev)
      integer i,jj,kk,jinc,k,l,ii2,ii3,jjp1,kkp1
      real xl(i_1:i_2,plev,4),xr(i_1:i_2,plev,4)
      real ys(i_1:i_2,plev),yn(i_1:i_2,plev)
      real zt(i_1:i_2,plev),zb(i_1:i_2,plev)
!
      do k=1,plev
        do i=i_1,i_2

          ii2  = idp(i,k,2)
          ii3  = idp(i,k,3)
          jj   = jdp(i,k)
          jjp1 = jj + 1
          kk   = kdp(i,k)
          kkp1 = kk + 1
!
          zt(i,k)   = (etamid(kkp1) - zdp(i,k)) / detam(kk)
          zb(i,k)   = 1. - zt(i,k)
!
          fdp1(i,k) =
     &    (( fb1 (ii2  ,kk  ,jj  ) * xl (i,k,2)
     &  +    fb1 (ii2+1,kk  ,jj  ) * xr (i,k,2) ) * ys(i,k)
     &  +  ( fb1 (ii3  ,kk  ,jjp1) * xl (i,k,3)
     &  +    fb1 (ii3+1,kk  ,jjp1) * xr (i,k,3) ) * yn(i,k) ) * zt(i,k)
     &  + (( fb1 (ii2  ,kkp1,jj  ) * xl (i,k,2)
     &  +    fb1 (ii2+1,kkp1,jj  ) * xr (i,k,2) ) * ys(i,k)
     &  +  ( fb1 (ii3  ,kkp1,jjp1) * xl (i,k,3)
     &  +    fb1 (ii3+1,kkp1,jjp1) * xr (i,k,3) ) * yn(i,k) ) * zb(i,k)
!
!$$$       print*,' inside int3dv i,k,i2,i3,kdp,jdp,xl2,xl3,xr,ys,yn,zt = ',
!$$$     . i,k,idp(i,k,2),idp(i,k,3),kdp(i,k),jdp(i,k),xl(i,k,2),xl(i,k,3),
!$$$     . ys(i,k),yn(i,k),zt(i,k)

          fdp2(i,k) =
     &    (( fb2 (ii2  ,kk  ,jj  ) * xl (i,k,2)
     &  +    fb2 (ii2+1,kk  ,jj  ) * xr (i,k,2) ) * ys(i,k)
     &  +  ( fb2 (ii3  ,kk  ,jjp1) * xl (i,k,3)
     &  +    fb2 (ii3+1,kk  ,jjp1) * xr (i,k,3) ) * yn(i,k) ) * zt(i,k)
     &  + (( fb2 (ii2  ,kkp1,jj  ) * xl (i,k,2)
     &  +    fb2 (ii2+1,kkp1,jj  ) * xr (i,k,2) ) * ys(i,k)
     &  +  ( fb2 (ii3  ,kkp1,jjp1) * xl (i,k,3)
     &  +    fb2 (ii3+1,kkp1,jjp1) * xr (i,k,3) ) * yn(i,k) ) * zb(i,k)
!
          fdp3(i,k) =
     &   (( fb3 (ii2  ,kk  ,jj  ) * xl (i,k,2)
     &  +    fb3 (ii2+1,kk  ,jj  ) * xr (i,k,2) ) * ys(i,k)
     &  +  ( fb3 (ii3  ,kk  ,jjp1) * xl (i,k,3)
     &  +    fb3 (ii3+1,kk  ,jjp1) * xr (i,k,3) ) * yn(i,k) ) * zt(i,k)
     &  + (( fb3 (ii2  ,kkp1,jj  ) * xl (i,k,2)
     &  +    fb3 (ii2+1,kkp1,jj  ) * xr (i,k,2) ) * ys(i,k)
     &  +  ( fb3 (ii3  ,kkp1,jjp1) * xl (i,k,3)
     &  +    fb3 (ii3+1,kkp1,jjp1) * xr (i,k,3) ) * yn(i,k) ) * zb(i,k)
        enddo
      enddo
      return
      end
      subroutine lagxin_h(i_1,i_2,s_lat_in_pe,n_lat_in_pe,
     &                    kdim,fb,idp,jdp,kdp,kkdp,
     &                    x2,x3,term1,term2,term3,term4,fint)
!sk 2/25/2012
!sk modified original subroutine for quasi-monotone 
!sk quasi-cubic Lagrange interpolation
!sk 
      use      pmgrid        , only : plev,plond,quamon
      implicit none
      integer i_1,i_2,s_lat_in_pe,n_lat_in_pe
      integer, parameter ::
     &                      ppdy=4,      ! length of interp. grid stencil in y
     &                      ppdz=4       ! length of interp. grid stencil in z
!
      integer kdim,idp(i_1:i_2,plev,4),jdp (i_1:i_2,plev),
     &             kdp(i_1:i_2,plev),kkdp(i_1:i_2,plev)
      real    fb(plond,plev,s_lat_in_pe:n_lat_in_pe),
     &        fint(i_1:i_2,plev,ppdy,ppdz)
      integer i,k,kdimm1,kdimm2
!     integer i,k,kdimm1,kdimm2,nclamp
      real    term1(i_1:i_2,plev,4),term2(i_1:i_2,plev,4),
     &        term3(i_1:i_2,plev,4),term4(i_1:i_2,plev,4)
      real    x2(i_1:i_2,plev,4),   x3(i_1:i_2,plev,4)
      integer ii1,ii2,ii3,ii4,jj,kk
      real    fmax(i_1:i_2), fmin(i_1:i_2)
!
!     nclamp = 93
      kdimm1 = kdim - 1
      kdimm2 = kdim - 2

      do k=1,plev
        do i=i_1,i_2
          ii1 = idp(i,k,1)
          ii2 = idp(i,k,2)

          ii3 = idp(i,k,3)
          ii4 = idp(i,k,4)

          jj = jdp(i,k)
          kk = kkdp(i,k)
!
          fint(i,k,2,1) = fb (ii2+1,kk-1,jj  ) * x2   (i,k,2)
     &                  - fb (ii2  ,kk-1,jj  ) * x3   (i,k,2)

          fint(i,k,3,1) = fb (ii3+1,kk-1,jj+1) * x2   (i,k,3)
     &                  - fb (ii3  ,kk-1,jj+1) * x3   (i,k,3)

          fint(i,k,1,2) = fb (ii1+1,kk  ,jj-1) * x2   (i,k,1)
     &                  - fb (ii1  ,kk  ,jj-1) * x3   (i,k,1)

          fint(i,k,2,2) = fb (ii2-1,kk  ,jj  ) * term1(i,k,2)
     &                  + fb (ii2  ,kk  ,jj  ) * term2(i,k,2)
     &                  + fb (ii2+1,kk  ,jj  ) * term3(i,k,2)
     &                  + fb (ii2+2,kk  ,jj  ) * term4(i,k,2)

          fint(i,k,3,2) = fb (ii3-1,kk  ,jj+1) * term1(i,k,3)
     &                  + fb (ii3  ,kk  ,jj+1) * term2(i,k,3)
     &                  + fb (ii3+1,kk  ,jj+1) * term3(i,k,3)
     &                  + fb (ii3+2,kk  ,jj+1) * term4(i,k,3)

          if (quamon) then
            fmax(i) = max(fb (ii2-1,kk  ,jj  ), fb (ii2  ,kk  ,jj  ),
     &                    fb (ii2+1,kk  ,jj  ), fb (ii2+2,kk  ,jj  ))
            fmin(i) = min(fb (ii2-1,kk  ,jj  ), fb (ii2  ,kk  ,jj  ),
     &                    fb (ii2+1,kk  ,jj  ), fb (ii2+2,kk  ,jj  ))
            fint(i,k,2,2) = max(fmin(i),min(fmax(i),fint(i,k,2,2)))
            fmax(i) = max(fb (ii3-1,kk  ,jj+1), fb (ii3  ,kk  ,jj+1),
     &                    fb (ii3+1,kk  ,jj+1), fb (ii3+2,kk  ,jj+1))
            fmin(i) = min(fb (ii3-1,kk  ,jj+1), fb (ii3  ,kk  ,jj+1),
     &                    fb (ii3+1,kk  ,jj+1), fb (ii3+2,kk  ,jj+1))
            fint(i,k,3,2) = max(fmin(i),min(fmax(i),fint(i,k,3,2)))

!           call quasim(fb(ii2,kk,jj),fb(ii2+1,kk,jj),fint(i,k,2,2))
!           call quasim(fb(ii3,kk,jj+1),fb(ii3+1,kk,jj+1),fint(i,k,3,2))
          endif

          fint(i,k,4,2) = fb (ii4+1,kk  ,jj+2) * x2   (i,k,4)
     &                  - fb (ii4  ,kk  ,jj+2) * x3   (i,k,4)

          fint(i,k,1,3) = fb (ii1+1,kk+1,jj-1) * x2   (i,k,1)
     &                  - fb (ii1 , kk+1,jj-1) * x3   (i,k,1)

          fint(i,k,2,3) = fb (ii2-1,kk+1,jj  ) * term1(i,k,2)
     &                  + fb (ii2  ,kk+1,jj  ) * term2(i,k,2)
     &                  + fb (ii2+1,kk+1,jj  ) * term3(i,k,2)
     &                  + fb (ii2+2,kk+1,jj  ) * term4(i,k,2)

          fint(i,k,3,3) = fb (ii3-1,kk+1,jj+1) * term1(i,k,3)
     &                  + fb (ii3  , kk+1,jj+1) * term2(i,k,3)
     &                  + fb (ii3+1,kk+1,jj+1) * term3(i,k,3)
     &                  + fb (ii3+2,kk+1,jj+1) * term4(i,k,3)

         if (quamon) then
           fmax(i) = max(fb (ii2-1,kk+1,jj  ), fb (ii2  ,kk+1,jj  ),
     &                   fb (ii2+1,kk+1,jj  ), fb (ii2+2,kk+1,jj  ))
           fmin(i) = min(fb (ii2-1,kk+1,jj  ), fb (ii2  ,kk+1,jj  ),
     &                   fb (ii2+1,kk+1,jj  ), fb (ii2+2,kk+1,jj  ))
           fint(i,k,2,3) = max(fmin(i), min(fmax(i),fint(i,k,2,3)))
           fmax(i) = max(fb (ii3-1,kk+1,jj+1), fb (ii3  , kk+1,jj+1),
     &                   fb (ii3+1,kk+1,jj+1), fb (ii3+2,kk+1,jj+1))
           fmin(i) = min(fb (ii3-1,kk+1,jj+1), fb (ii3  , kk+1,jj+1),
     &                   fb (ii3+1,kk+1,jj+1), fb (ii3+2,kk+1,jj+1))
           fint(i,k,3,3) = max(fmin(i), min(fmax(i),fint(i,k,3,3)))

!          call quasim(fb(ii2,kk+1,jj),   fb(ii2+1,kk+1,jj),
!    &                                    fint(i,k,2,3))
!          call quasim(fb(ii3,kk+1,jj+1), fb(ii3+1,kk+1,jj+1),
!    &                                    fint(i,k,3,3))
         endif

          fint(i,k,4,3) = fb (ii4+1,kk+1,jj+2) * x2   (i,k,4)
     &                  - fb (ii4  ,kk+1,jj+2) * x3   (i,k,4)

          fint(i,k,2,4) = fb (ii2+1,kk+2,jj  ) * x2   (i,k,2)
     &                  - fb (ii2  ,kk+2,jj  ) * x3   (i,k,2)

          fint(i,k,3,4) = fb (ii3+1,kk+2,jj+1) * x2   (i,k,3)
     &                  - fb (ii3  ,kk+2,jj+1) * x3   (i,k,3)

        enddo
      enddo
      do k=1,plev
        do i=i_1,i_2
          ii1 = idp(i,k,1)
          ii2 = idp(i,k,2)
          ii3 = idp(i,k,3)
          ii4 = idp(i,k,4)
          jj  = jdp(i,k)
!
          if(kdp (i,k) ==  1) then
            fint(i,k,2,1) = fint(i,k,2,4)
            fint(i,k,3,1) = fint(i,k,3,4)
            fint(i,k,1,3) = fint(i,k,1,2)
            fint(i,k,2,3) = fint(i,k,2,2)
            fint(i,k,3,3) = fint(i,k,3,2)
            fint(i,k,4,3) = fint(i,k,4,2)

            fint(i,k,1,2) = fb (ii1+1,1,jj-1) * x2   (i,k,1)
     &                    - fb (ii1  ,1,jj-1) * x3   (i,k,1)

            fint(i,k,2,2) = fb (ii2-1,1,jj  ) * term1(i,k,2)
     &                    + fb (ii2  ,1,jj  ) * term2(i,k,2)
     &                    + fb (ii2+1,1,jj  ) * term3(i,k,2)
     &                    + fb (ii2+2,1,jj  ) * term4(i,k,2)

            fint(i,k,3,2) = fb (ii3-1,1,jj+1) * term1(i,k,3)
     &                    + fb (ii3  ,1,jj+1) * term2(i,k,3)
     &                    + fb (ii3+1,1,jj+1) * term3(i,k,3)
     &                    + fb (ii3+2,1,jj+1) * term4(i,k,3)


            if (quamon) then
              fmax(i) = max(fb (ii2-1,1,jj  ), fb (ii2  ,1,jj  ),
     &                      fb (ii2+1,1,jj  ), fb (ii2+2,1,jj  ))
              fmin(i) = min(fb (ii2-1,1,jj  ), fb (ii2  ,1,jj  ),
     &                      fb (ii2+1,1,jj  ), fb (ii2+2,1,jj  ))
              fint(i,k,2,2) = max(fmin(i),
     &                        min(fmax(i),fint(i,k,2,2)))
              fmax(i) = max(fb (ii3-1,1,jj+1), fb (ii3  ,1,jj+1),
     &                      fb (ii3+1,1,jj+1), fb (ii3+2,1,jj+1))
              fmin(i) = min(fb (ii3-1,1,jj+1), fb (ii3  ,1,jj+1),
     &                      fb (ii3+1,1,jj+1), fb (ii3+2,1,jj+1))
              fint(i,k,3,2) = max(fmin(i),
     &                        min(fmax(i), fint(i,k,3,2)))
!
!             call quasim(fb(ii2,1,jj),   fb(ii2+1,1,jj),
!    &                                    fint(i,k,2,2))
!             call quasim(fb(ii3,1,jj+1), fb(ii3+1,1,jj+1),
!    &                                    fint(i,k,3,2))
            endif

            fint(i,k,4,2) = fb (ii4+1,1,jj+2) * x2   (i,k,4)
     &                    - fb (ii4  ,1,jj+2) * x3   (i,k,4)

            fint(i,k,2,4) = fb (ii2+1,3,jj  ) * x2   (i,k,2)
     &                    - fb (ii2  ,3,jj  ) * x3   (i,k,2)

            fint(i,k,3,4) = fb (ii3+1,3,jj+1) * x2   (i,k,3)
     &                    - fb (ii3  ,3,jj+1) * x3   (i,k,3)

          elseif(kdp (i,k) == kdimm1) then
            fint(i,k,2,4) = fint(i,k,2,1)
            fint(i,k,3,4) = fint(i,k,3,1)
            fint(i,k,1,2) = fint(i,k,1,3)
            fint(i,k,2,2) = fint(i,k,2,3)
            fint(i,k,3,2) = fint(i,k,3,3)
            fint(i,k,4,2) = fint(i,k,4,3)

            fint(i,k,2,1) = fb (ii2+1,kdimm2,jj  ) * x2   (i,k,2)
     &                    - fb (ii2  ,kdimm2,jj  ) * x3   (i,k,2)

            fint(i,k,3,1) = fb (ii3+1,kdimm2,jj+1) * x2   (i,k,3)
     &                    - fb (ii3  ,kdimm2,jj+1) * x3   (i,k,3)

            fint(i,k,1,3) = fb (ii1+1,kdim  ,jj-1) * x2   (i,k,1)
     &                    - fb (ii1  ,kdim  ,jj-1) * x3   (i,k,1)

            fint(i,k,2,3) = fb (ii2-1,kdim  ,jj  ) * term1(i,k,2)
     &                    + fb (ii2  ,kdim  ,jj  ) * term2(i,k,2)
     &                    + fb (ii2+1,kdim  ,jj  ) * term3(i,k,2)
     &                    + fb (ii2+2,kdim  ,jj  ) * term4(i,k,2)
            fint(i,k,3,3) = fb (ii3-1,kdim  ,jj+1) * term1(i,k,3)
     &                    + fb (ii3  ,kdim  ,jj+1) * term2(i,k,3)
     &                    + fb (ii3+1,kdim  ,jj+1) * term3(i,k,3)
     &                    + fb (ii3+2,kdim  ,jj+1) * term4(i,k,3)
            if (quamon) then
              fmax(i) = max(fb (ii2-1,kdim  ,jj  ),
     &                      fb (ii2  ,kdim  ,jj  ),
     &                      fb (ii2+1,kdim  ,jj  ),
     &                      fb (ii2+2,kdim  ,jj ))
              fmin(i) = min(fb (ii2-1,kdim  ,jj  ),
     &                      fb (ii2  ,kdim  ,jj ),
     &                      fb (ii2+1,kdim  ,jj  ),
     &                      fb (ii2+2,kdim  ,jj ))
              fint(i,k,2,3) = max(fmin(i),
     &                        min(fmax(i), fint(i,k,2,3)))
              fmax(i) = max(fb (ii3-1,kdim  ,jj+1),
     &                      fb (ii3  ,kdim ,jj+1),
     &                      fb (ii3+1,kdim  ,jj+1),
     &                      fb (ii3+2,kdim ,jj+1))
              fmin(i) = min(fb (ii3-1,kdim  ,jj+1),
     &                      fb (ii3  ,kdim ,jj+1),
     &                      fb (ii3+1,kdim  ,jj+1),
     &                      fb (ii3+2,kdim ,jj+1))
              fint(i,k,3,3) = max(fmin(i),
     &                        min(fmax(i), fint(i,k,3,3)))
!
!             call quasim(fb(ii2,kdim,jj),   fb(ii2+1,kdim,jj),
!    &                                       fint(i,k,2,3))
!             call quasim(fb(ii3,kdim,jj+1), fb(ii3+1,kdim,jj+1),
!    &                                       fint(i,k,3,3))
            endif

            fint(i,k,4,3) = fb (ii4+1,kdim  ,jj+2) * x2   (i,k,4)
     &                    - fb (ii4  ,kdim  ,jj+2) * x3   (i,k,4)
          endif
        enddo
      enddo
      return
      end



      subroutine lagxinit_h(i_1,i_2,
     &                      x,xdp,idp,jdp,x2,x3,term1,term2,term3,term4)
!
      use      pmgrid        , only : platd,plev,plond
      use slgshr             , only : rdlam
      implicit none
!
      integer i_1,i_2
      integer       idp(i_1:i_2,plev,4)
      real x(plond,platd),xdp(i_1:i_2,plev)
      integer jdp(i_1:i_2,plev)
!     integer jdp(i_1:i_2,plev),jj(i_1:i_2,plev,4)
      integer i,kk,jinc,k,l,kinc,kdimm1,kdimm2,n,jj
      real denom1,denom2,denom3,denom4,coef12,coef34
      real  term1(i_1:i_2,plev,4),term2(i_1:i_2,plev,4),
     $      term3(i_1:i_2,plev,4),term4(i_1:i_2,plev,4)
      real x2(i_1:i_2,plev,4),   x3(i_1:i_2,plev,4)

      denom1 = -1./6.
      denom2 =  0.5
      denom3 = -0.5
      denom4 =  1./6.
      
      do n=1,4
        do k=1,plev
          do i=i_1,i_2
            jj        = jdp(i,k) - 2 + n
            x2(i,k,n) = (xdp(i,k) - x(idp(i,k,n),jj)) * rdlam(jj)
            x3(i,k,n) = x2(i,k,n) - 1.
          end do
        end do
      end do
!
      do n=2,3
        do k=1,plev
          do i=i_1,i_2
            coef12       = x3(i,k,n)*(x2(i,k,n) - 2.)
            coef34       = (x2(i,k,n) + 1.)*x2(i,k,n)
            term1(i,k,n) = denom1*coef12*x2(i,k,n)
            term2(i,k,n) = denom2*coef12*(x2(i,k,n) + 1.)
            term3(i,k,n) = denom3*coef34*(x2(i,k,n) - 2.)
            term4(i,k,n) = denom4*coef34*x3(i,k,n)
          enddo
        enddo
      enddo
      return
      end
      subroutine xywgts_h(i_1,i_2,x,y,dy,idp,jdp,xdp,ydp,xl,xr,ys,yn)
!
      use      pmgrid        , only : platd,plev,plon,plond
      use layout1            , only : me
      implicit none
!
      integer i_1,i_2
      integer idp(i_1:i_2,plev,4),jdp(i_1:i_2,plev)
      real xdp(plon,plev),ydp(plon,plev),x(plond,platd)
      real y(platd),dy(platd)
      real xl(i_1:i_2,plev,4),xr(i_1:i_2,plev,4),
     $     ys(i_1:i_2,plev),yn(i_1:i_2,plev)
      integer i,j,k,jj2,jj3
      real rdx(platd)

      do j=1,platd
        rdx(j) = 1. / (x(2,j) - x(1,j))
      enddo
      do k=1,plev
        do i=i_1,i_2

          jj2 = jdp(i,k)
          jj3 = jdp(i,k) + 1
          xl(i,k,2) = (x(idp(i,k,2)+1,jj2) - xdp(i,k))*rdx(jj2)
          xl(i,k,3) = (x(idp(i,k,3)+1,jj3) - xdp(i,k))*rdx(jj3)

!$$$      if (me .eq. 43) print*,'i,k,jdp,i2,i3,x(i2j2),x(i3j3),xl2,xl3',
!$$$   .  i,k,jdp(i,k),idp(i,k,2),idp(i,k,3),x(idp(i,k,2)+1,jj2),
!$$$   .  x(idp(i,k,3)+1,jj3),xdp(idp),xl(i,k,2),xl(i,k,3)

          xr(i,k,2) = 1. - xl(i,k,2)
          xr(i,k,3) = 1. - xl(i,k,3)
          ys(i,k)   = (y(jdp(i,k)+1) - ydp(i,k)) / dy(jdp(i,k))
          yn(i,k)   = 1. - ys(i,k)
        enddo
      enddo

      return
      end
      subroutine rotate_uv_h(i_1,i_2, lon_dim_h,
     . ud3_d,vd3_d,ud3_a,vd3_a,
     . cos_l_a,sin_l_a,cos_f_a,sin_f_a,
     . cos_l_d,sin_l_d,cos_f_d,sin_f_d)

      use      pmgrid        , only : plev,plon
      implicit none
      integer i_1,i_2,i,k,lon_dim_h

      real cos_l_a(i_1:i_2),sin_l_a(i_1:i_2),
     &cos_f_a,sin_f_a
!watch lon_dim_h dimension


      real cos_l_d(i_1:i_2,plev),
     $     sin_l_d(i_1:i_2,plev),
     $     cos_f_d(i_1:i_2,plev),
     $     sin_f_d(i_1:i_2,plev),
     &       ud3_d(plon   ,plev),
     &       vd3_d(plon   ,plev),
     &     alfa_ua(i_1:i_2,plev),
     &     alfa_va(i_1:i_2,plev),
     &     beta_ua(i_1:i_2,plev),
     &     beta_va(i_1:i_2,plev)

      real ud3_a(plon,plev),vd3_a(plon,plev)
      real cos_l_diff,sin_l_diff

! the trigonometric functions are bottom to top, they are inverted
! in dep3dg and dep3ds, therefore alfas and betas are bottom to top.
      do k = 1,plev
       do i = i_1,i_2
          cos_l_diff = cos_l_a(i)*cos_l_d(i,k)+sin_l_a(i)*sin_l_d(i,k)
          sin_l_diff = sin_l_a(i)*cos_l_d(i,k)-sin_l_d(i,k)*cos_l_a(i)

          alfa_ua(i,k) = cos_l_diff
          alfa_va(i,k) = sin_f_d(i,k)*sin_l_diff

          beta_ua(i,k) = -sin_f_a*sin_l_diff
          beta_va(i,k) = cos_f_d(i,k)*cos_f_a + 
     .                   sin_f_d(i,k)*sin_f_a*cos_l_diff
       enddo
      enddo
!
!sela     do k = 1,plev
!sela      do i = i_1,i_2
!sela      alfa_ua(i,k)=1.
!sela      alfa_va(i,k)=0.
!sela
!sela      beta_ua(i,k)=0.
!sela      beta_va(i,k)=1.
!sela      enddo
!sela     enddo
! ud3_d, vd3_d are bottom to top, coming from lagzin.
! alfas,betas  are bottom to top, coming from previous loops.
!therefore:
! ud3_a, vd3_a are bottom to top, as expected in gloopa.
      do k = 1,plev
       do i = i_1,i_2
        ud3_a(i,      k) = alfa_ua(i,k) * ud3_d(i,k)
     &                   + alfa_va(i,k) * vd3_d(i,k)

        vd3_a(i,      k) = beta_ua(i,k) * ud3_d(i,k)
     &                   + beta_va(i,k) * vd3_d(i,k)
       end do
      end do
      return
      end
!
      function bisection(x,plev,u)
      implicit none
      integer u,i,j,k,plev,bisection
      integer x(plev)
       i = 1
       j = plev
       do
         k = (i+j)/2
         if (u < x(k)) then
           j = k
         else
           i = k
         end if
         if (i+1 >= j)  then
           bisection = i
           return
         endif
       end do
      end
      subroutine endrun_h
      implicit none
      print*,' **** endrun_h was called - execution stopped '
      stop13
      end
      subroutine lagyin_h(i_1,i_2,
     &                    fintx,finty,yb,yt,term1,term2,term3,term4)
!sk 2/25/2012
!sk modified original subroutine for quasi-monotone
!sk quasi-cubic Lagrange interpolation
!sk
      use      pmgrid        , only : plev,quamon
      implicit none
      integer i_1,i_2
      real fintx(i_1:i_2,plev,4   ,4   ) ! x-interpolants
      real finty(i_1:i_2,plev,4   )  ! interpolants at the horiz. depart
      integer
     &     i,                  ! index
     &     k                   ! index
      real
     &     term1(i_1:i_2,plev),         ! |
     &     term2(i_1:i_2,plev),         ! |
     &     term3(i_1:i_2,plev),         ! |
     &     term4(i_1:i_2,plev),         ! |
     &     yb   (i_1:i_2,plev),
     &     yt   (i_1:i_2,plev)
      real fmax(i_1:i_2), fmin(i_1:i_2)
!
      do k=1,plev
        do i=i_1,i_2
          finty(i,k,1) =   fintx(i,k,2,1)*yb   (i,k)
     &                   + fintx(i,k,3,1)*yt   (i,k)
          finty(i,k,2) =   fintx(i,k,1,2)*term1(i,k)
     &                   + fintx(i,k,2,2)*term2(i,k)
     &                   + fintx(i,k,3,2)*term3(i,k)
     &                   + fintx(i,k,4,2)*term4(i,k)
          finty(i,k,3) =   fintx(i,k,1,3)*term1(i,k)
     &                   + fintx(i,k,2,3)*term2(i,k)
     &                   + fintx(i,k,3,3)*term3(i,k)
     &                   + fintx(i,k,4,3)*term4(i,k)
          if (quamon) then
            fmax(i) = max(fintx(i,k,1,2), fintx(i,k,2,2),
     &                    fintx(i,k,3,2), fintx(i,k,4,2))
            fmin(i) = min(fintx(i,k,1,2), fintx(i,k,2,2),
     &                    fintx(i,k,3,2), fintx(i,k,4,2))
            finty(i,k,2) = max(fmin(i), min(fmax(i), finty(i,k,2)))
            fmax(i) = max(fintx(i,k,1,3), fintx(i,k,2,3),
     &                    fintx(i,k,3,3), fintx(i,k,4,3))
            fmin(i) = min(fintx(i,k,1,3), fintx(i,k,2,3),
     &                    fintx(i,k,3,3), fintx(i,k,4,3))
            finty(i,k,3) = max(fmin(i), min(fmax(i),finty(i,k,3)))

!           call quasim(fintx(i,k,2,2),fintx(i,k,3,2),finty(i,k,2))
!           call quasim(fintx(i,k,2,3),fintx(i,k,3,3),finty(i,k,3))
          endif
          finty(i,k,4) =   fintx(i,k,2,4)*yb   (i,k)
     &                   + fintx(i,k,3,4)*yt   (i,k)
        enddo
      enddo
      return
      end


      subroutine lagyinit_h(i_1,i_2,y,dy,lbasiy,ydp,jdp,
     &                      yb,yt,term1,term2,term3,term4)
!
      use      pmgrid        , only : platd,plev
      implicit none
      integer i_1,i_2
      real    y (platd),     dy(platd),
     &        lbasiy(4,2,platd),      ! y-interpolation weights
     &        ydp(i_1:i_2,plev)    ! y-coordinates of departure pts.
      integer jdp(i_1:i_2,plev)    ! j-index of departure point coord.
      integer i, k, m, jj
      real ymy1,                ! |
     &     ymy2,                ! |
     &     ymy3,                ! |
     &     ymy4,                ! |
     &     coef12,              ! |
     &     coef34,              ! | -- interpolation weights/coeffs.
     &     term1(i_1:i_2,plev), ! |
     &     term2(i_1:i_2,plev), ! |
     &     term3(i_1:i_2,plev), ! |
     &     term4(i_1:i_2,plev), ! |
     &     yb   (i_1:i_2,plev),
     &     yt   (i_1:i_2,plev)
      real dyj
!
113    format(5(i6,1x))

       do k=1,plev
         do i=i_1,i_2
           jj = jdp(i,k)
           dyj        =  dy(jj)
           yb(i,k)    = ( y(jj+1) - ydp(i,k) )/dyj
           yt(i,k)    = 1. - yb(i,k)
           ymy1       = ydp(i,k) - lbasiy(1,1,jj)
           ymy2       = ydp(i,k) - lbasiy(2,1,jj)
           ymy3       = ydp(i,k) - lbasiy(3,1,jj)
           ymy4       = ydp(i,k) - lbasiy(4,1,jj)
           coef12     = ymy3*ymy4
           coef34     = ymy1*ymy2
           term1(i,k) = coef12*ymy2*lbasiy(1,2,jj)
           term2(i,k) = coef12*ymy1*lbasiy(2,2,jj)
           term3(i,k) = coef34*ymy4*lbasiy(3,2,jj)
           term4(i,k) = coef34*ymy3*lbasiy(4,2,jj)
         end do
       end do

      return
      end
      subroutine lagzin_h(i_1,i_2,
     &                  kdim,finty,lbasdz,kdp,hb,ht,dhb,dht,detam,
     &                  term1,term2,term3,term4,fdp)
!
      use      pmgrid        , only : plev,plevp,plon
      implicit none
!
      integer i_1,i_2
      integer kdim,kdp(i_1:i_2,plev)
      real       finty(i_1:i_2,plev,4),
     &     lbasdz(4,2,kdim),fdp(plon,plev)
      integer i,k,m,kdimm1,kdimm2,kk,kr
      real  ht(i_1:i_2,plev), hb(i_1:i_2,plev),
     &     dht(i_1:i_2,plev),dhb(i_1:i_2,plev),
     &     ftop,fbot
      real term1(i_1:i_2,plev),    ! |
     &     term2(i_1:i_2,plev),    ! |
     &     term3(i_1:i_2,plev),    ! |
     &     term4(i_1:i_2,plev)     ! |
      real detam(plev)  ! delta eta at levels
!
      kdimm1 = kdim - 1
      kdimm2 = kdim - 2

      do k = 1,plev
        kr = plevp-k
        do i = i_1,i_2
          kk = kdp(i,k)
          if(kk == 1) then
            ftop       = (finty(i,k,3) - finty(i,k,2)) / detam(1)
            fbot       = lbasdz(1,1,2) * finty(i,k,2)
     &                 + lbasdz(2,1,2) * finty(i,k,3)
     &                 + lbasdz(3,1,2) * finty(i,k,4)
     &                 + lbasdz(4,1,2) * finty(i,k,1)
            fdp(i,kr) = finty(i,k,2)   * ht (i,k)
     &                + ftop           * dht(i,k)
     &                + finty(i,k,3)   * hb (i,k)
     &                + fbot           * dhb(i,k)
          elseif(kk == kdimm1) then
            ftop       = lbasdz(1,2,kdimm2) * finty(i,k,4)+
     &                   lbasdz(2,2,kdimm2) * finty(i,k,1)+
     &                   lbasdz(3,2,kdimm2) * finty(i,k,2)+
     &                   lbasdz(4,2,kdimm2) * finty(i,k,3)
!$$$        fbot       = 0.0
            fdp(i,kr) = finty(i,k,2) * ht (i,k)
     &                + ftop         * dht(i,k)
     &                + finty(i,k,3) * hb (i,k)
!$$$ &                + fbot         * dhb(i,k)
          else
            fdp(i,kr) = finty(i,k,1) * term1(i,k)
     &                + finty(i,k,2) * term2(i,k)
     &                + finty(i,k,3) * term3(i,k)
     &                + finty(i,k,4) * term4(i,k)
          endif
        enddo
      enddo
      return
      end
      subroutine lagzinit_h(i_1,i_2,kdim,lbasiz,
     &                      etamid,detam,sigdp,kdp,kkdp,
     &                      hb,ht,dhb,dht,term1z,term2z,term3z,term4z)
!
      use      pmgrid        , only : plev
      implicit none
!
      integer i_1,i_2
      integer kdim,kdp(i_1:i_2,plev),kkdp(i_1:i_2,plev)
      real lbasiz(4,2,kdim),
     &  etamid(kdim),detam(kdim),
     &  sigdp(i_1:i_2,plev),fdp(i_1:i_2,plev)
      integer i,k,m,kdimm1,kdimm2,kk
      real dzk,zt,zb,zt2,zb2,
     &     ht(i_1:i_2,plev), hb(i_1:i_2,plev),
     &    dht(i_1:i_2,plev),dhb(i_1:i_2,plev)
      real zmz1,                 ! |
     &     zmz2,                 ! |
     &     zmz3,                 ! |
     &     zmz4,                 ! |
     &     coef12,               ! |
     &     coef34,               ! | -- interpolation weights/coeffs.
     &     term1z(i_1:i_2,plev), ! |
     &     term2z(i_1:i_2,plev), ! |
     &     term3z(i_1:i_2,plev), ! |
     &     term4z(i_1:i_2,plev)  ! |
!
      do k=1,plev
        do i=i_1,i_2
          dzk  =  detam(kdp(i,k))
          zt   = ( etamid(kdp(i,k)+1) - sigdp(i,k) )/dzk
          zb   = 1. - zt
          zt2  = zt * zt
          zb2  = zb * zb
          kk   = kkdp(i,k)
          ht(i,k)     = ( 3.0 - 2.0*zt )  * zt2
          hb(i,k)     = ( 3.0 - 2.0*zb )  * zb2
          dht(i,k)    = -dzk*( zt - 1. ) * zt2
          dhb(i,k)    =  dzk*( zb - 1. ) * zb2
          zmz1        = sigdp(i,k) -  lbasiz(1,1,kk)
          zmz2        = sigdp(i,k) -  lbasiz(2,1,kk)
          zmz3        = sigdp(i,k) -  lbasiz(3,1,kk)
          zmz4        = sigdp(i,k) -  lbasiz(4,1,kk)
          coef12      = zmz3 * zmz4
          coef34      = zmz1 * zmz2
          term1z(i,k) = coef12 * zmz2 * lbasiz(1,2,kk)
          term2z(i,k) = coef12 * zmz1 * lbasiz(2,2,kk)
          term3z(i,k) = coef34 * zmz4 * lbasiz(3,2,kk)
          term4z(i,k) = coef34 * zmz3 * lbasiz(4,2,kk)
        end do
      end do
      return
      end
      subroutine lcbas_h( grd, bas1, bas2 )
      implicit none
      real grd(4)               ! grid stencil
      real bas1(4),             ! grid values on stencil
     &     bas2(4)              ! lagrangian basis functions
      real x0mx1,               ! |
     &     x0mx2,               ! |
     &     x0mx3,               ! |- grid value differences used in
     &     x1mx2,               ! |  weights
     &     x1mx3,               ! |
     &     x2mx3                ! |
      x0mx1   = grd(1) - grd(2)
      x0mx2   = grd(1) - grd(3)
      x0mx3   = grd(1) - grd(4)
      x1mx2   = grd(2) - grd(3)
      x1mx3   = grd(2) - grd(4)
      x2mx3   = grd(3) - grd(4)
      bas1(1) = grd(1)
      bas1(2) = grd(2)
      bas1(3) = grd(3)
      bas1(4) = grd(4)
      bas2(1) =  1./ ( x0mx1 * x0mx2 * x0mx3 )
      bas2(2) = -1./ ( x0mx1 * x1mx2 * x1mx3 )
      bas2(3) =  1./ ( x0mx2 * x1mx2 * x2mx3 )
      bas2(4) = -1./ ( x0mx3 * x1mx3 * x2mx3 )
      return
      end
      subroutine lcdbas_h(grd,dbas2,dbas3)
      implicit none
      real grd(4)               ! grid stencil
      real dbas2(4),            ! derivatives at grid point 2.
     $     dbas3(4)             ! derivatives at grid point 3.
      real x1,                  ! |
     $     x2,                  ! |- grid values
     $     x3,                  ! |
     $     x4,                  ! |
     $     x1mx2,               !  |
     $     x1mx3,               !  |
     $     x1mx4,               !  |- differences of grid values
     $     x2mx3,               !  |
     $     x2mx4,               !  |
     $     x3mx4                !  |
      x1 = grd(1)
      x2 = grd(2)
      x3 = grd(3)
      x4 = grd(4)
      x1mx2 = x1 - x2
      x1mx3 = x1 - x3
      x1mx4 = x1 - x4
      x2mx3 = x2 - x3
      x2mx4 = x2 - x4
      x3mx4 = x3 - x4
      dbas2(1) =   x2mx3 * x2mx4 / ( x1mx2 * x1mx3 * x1mx4 )
      dbas2(2) =   -1./x1mx2 + 1./x2mx3 + 1./x2mx4
      dbas2(3) = - x1mx2 * x2mx4 / ( x1mx3 * x2mx3 * x3mx4 )
      dbas2(4) =   x1mx2 * x2mx3 / ( x1mx4 * x2mx4 * x3mx4 )
      dbas3(1) = - x2mx3 * x3mx4 / ( x1mx2 * x1mx3 * x1mx4 )
      dbas3(2) =   x1mx3 * x3mx4 / ( x1mx2 * x2mx3 * x2mx4 )
      dbas3(3) =   -1./x1mx3 - 1./x2mx3 + 1./x3mx4
      dbas3(4) = - x1mx3 * x2mx3 / ( x1mx4 * x2mx4 * x3mx4 )
      return
      end
      subroutine set_pmgrid_h(lonf,latg,levs)
      use pmgrid , only : mprec,pgls,plat,platd,
     &                    plev,plevd,plevp,plon,plond
      implicit none
      integer lonf,latg,levs,nxpt,jintmx
!
      plon   = lonf
      plev   = levs
      plat   = latg
      plevp  = plev + 1
      nxpt   = 1
      jintmx = 1
      plond  = plon + 1 + 2*nxpt
      platd  = plat + 2*nxpt + 2*jintmx
      plevd  = plev
      mprec  = 1.e-12
      pgls   = plon*plev

      return
      end
      subroutine herzinit_h(i_1,i_2,kdim,
     &                      etamid,detam,sigdp,kdp,
     &                      hb,ht,dhb,dht,zb,zt)
!
      use      pmgrid        , only : plev
      implicit none
!
      integer i_1,i_2
      integer kdim,kdp(i_1:i_2,plev)
      real etamid(kdim),detam(kdim),
     &  sigdp(i_1:i_2,plev),fdp(i_1:i_2,plev)
      integer i,k,m,kdimm1,kdimm2
      real dzk,zb2,zt2,
     &     ht(i_1:i_2,plev), hb(i_1:i_2,plev),
     &    dht(i_1:i_2,plev),dhb(i_1:i_2,plev),
     &    zt(i_1:i_2,plev),zb(i_1:i_2,plev)
       do k=1,plev
         do i=i_1,i_2
            dzk  =  detam(kdp(i,k))
            zt(i,k)   = ( etamid(kdp(i,k)+1) - sigdp(i,k) )/dzk
            zb(i,k)   = 1. - zt(i,k)
            zt2       = zt(i,k) * zt(i,k)
            zb2       = zb(i,k) * zb(i,k)
             ht(i,k)  = ( 3.0 - 2.0*zt(i,k) ) * zt2
             hb(i,k)  = ( 3.0 - 2.0*zb(i,k) ) * zb2
            dht(i,k)  = -dzk*( zt(i,k) - 1. ) * zt2
            dhb(i,k)  =  dzk*( zb(i,k) - 1. ) * zb2
         enddo
       end do
      return
      end
      subroutine heryinit_h(i_1,i_2,phi,dphi,phidp,jdp,ys,yn,
     &                      hs,hn,dhs,dhn,rdphi)
     

      use      pmgrid        , only : platd,plev
      implicit none

      real phi (platd),
     &     dphi(platd),
     &     phidp(i_1:i_2,plev)    ! y-coordinates of departure pts.
      integer i_1,i_2,
     &     jdp(i_1:i_2,plev)    ! j-index of departure point coord.
      real ys(i_1:i_2,plev),yn(i_1:i_2,plev),
     &     hs(i_1:i_2,plev),hn(i_1:i_2,plev),
     &     dhs(i_1:i_2,plev),dhn(i_1:i_2,plev),
     &     rdphi(i_1:i_2,plev)
! locals
      real dyj,ys2,yn2
      integer i,k
!
      do k=1,plev
         do i=i_1,i_2
            dyj = dphi(jdp(i,k))
            rdphi(i,k) = 1.0/dyj
!
            ys(i,k)  = ( phi(jdp(i,k)+1) - phidp(i,k) )*rdphi(i,k)
            yn(i,k)  = 1. - ys(i,k)
            ys2      = ys(i,k) * ys(i,k)
            yn2      = yn(i,k) * yn(i,k)
!
            hs(i,k)  = (3.0 - 2.0*ys(i,k)) * ys2
            hn(i,k)  = (3.0 - 2.0*yn(i,k)) * yn2
!
            dhs(i,k) = -dyj * (ys(i,k) - 1.) * ys2
            dhn(i,k) =  dyj * (yn(i,k) - 1.) * yn2
        enddo
      enddo
!

      return
      end
      subroutine her_xyz_in_h(i_1,i_2,s_lat_in_pe,n_lat_in_pe,
     &                        idp,jdp,kdp,kkdp,xl,xr,hl,hr,dhl,dhr,fb,
     &                        fintx,ys,yn,hs,hn,dhs,dhn,rdphi,finty,
     &                        kdim,lbasdz,hb,ht,dhb,dht,detam,fdp)
!    &,        lprnt,jp,lan)
!-----------------------------------------------------------------------
      use      pmgrid        , only : plev,plevp,plon,plond
!     include 'slgshr.com.h'
      use      slgshr        , only : lbasdy,rdlam,rdlam6
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      integer i_1,i_2,lon_dim_h,s_lat_in_pe,n_lat_in_pe,jp,lan
      integer idp(i_1:i_2,plev,4),jdp(i_1:i_2,plev),kdp(i_1:i_2,plev)
      integer kkdp(i_1:i_2,plev),kdim,kdimm1,kdimm2
      real rdlamjj,rdlamjjp1,rdlam6jj,rdlam6jjp1
      real xl(i_1:i_2,plev,4), xr(i_1:i_2,plev,4),
     &     hl(i_1:i_2,plev,4), hr(i_1:i_2,plev,4),
     &    dhl(i_1:i_2,plev,4),dhr(i_1:i_2,plev,4)
      real fb(plond,plev,s_lat_in_pe:n_lat_in_pe)
      real fintx(i_1:i_2,plev,4,4)
      logical lprnt
! locals
      integer ii1,ii2,ii3,ii4,jj,kk,i,k
      real fxl,fxr,deli,tmp1,tmp2
      real fac
!-----------------------------------------------------------------------
      real ys(i_1:i_2,plev),
     &     yn(i_1:i_2,plev),
     &     hs(i_1:i_2,plev),
     &     hn(i_1:i_2,plev),
     &     dhs(i_1:i_2,plev),
     &     dhn(i_1:i_2,plev),
     &     rdphi(i_1:i_2,plev)
      real finty(i_1:i_2,plev,4)
! locals
      real ftop(i_1:i_2,plev,4),fbot(i_1:i_2,plev,4)
!-----------------------------------------------------------------------
      real lbasdz(4,2,kdim),fdp(plon,plev)
      integer m
      real  ht(i_1:i_2,plev), hb(i_1:i_2,plev),
     $     dht(i_1:i_2,plev),dhb(i_1:i_2,plev)
      real detam(plev)  ! delta eta at levels
!-----------------------------------------------------------------------
!     if (lprnt) then
!     print*,'hello  from her_xyz_in_h '
!     print*,' i_1=',i_1,' plon=',plon,' plev=',plev
!     print*,' i_2=',i_2
!     print*,' s_lat_in_pe=',s_lat_in_pe
!     print*,' n_lat_in_pe=',n_lat_in_pe
!     print*,' kdim=',kdim
!     print*,' fb=', fb(2,1:5,jp)
!!!   print*,' fdp=', fdp(1,1)
!!!   print*,' =',
!     print*,'goodby from her_xyz_in_h '
!     endif
!-----------------------------------------------------------------------
!
      fac  = 3.*(1. - 10.*epsilon(fac))
 !    print*,' her_xyz_in_h fac = ',fac

      kdimm1 = kdim - 1
      kdimm2 = kdim - 2
!
! PART 1:  x-interpolation
!
! Loop over fields.
! ..x interpolation at each height needed for z interpolation.
! ...x interpolation at each latitude needed for y interpolation.
!
      do k=1,plev
        do i=i_1,i_2
          ii1 = idp(i,k,1)
          ii2 = idp(i,k,2)
          ii3 = idp(i,k,3)
          ii4 = idp(i,k,4)
          jj  = jdp(i,k)
          kk  = kkdp(i,k)

          rdlamjj   = rdlam(jj)
          rdlamjjp1 = rdlam(jj+1)
!
          rdlam6jj   = rdlam6(jj)
          rdlam6jjp1 = rdlam6(jj+1)

!sela     if (ii2-1 .le. 0) then
!sela        print*,' zero index in herxin ii',ii2
!sela        stop221
!sela     endif
!sela     if (jj-1 .le. 0) then
!sela        print*,' zero index in herxin jj',jj
!sela        stop222
!sela     endif
!sela     if (kk-1 .le. 0) then
!sela        print*,' zero index in herxin kk',kk
!sela        stop223
!sela     endif
!

!
! Height level 1:  Linear interpolation on inner two latitudes only
!
!!!       fintx(i,k,1,1) = not used

          fintx(i,k,2,1) = fb(ii2  ,kk-1,jj  ) * xl(i,k,2) +
     &                     fb(ii2+1,kk-1,jj  ) * xr(i,k,2)
          fintx(i,k,3,1) = fb(ii3  ,kk-1,jj+1) * xl(i,k,3) +
     &                     fb(ii3+1,kk-1,jj+1) * xr(i,k,3)

!!!       fintx(i,k,4,1) = not used
!
! Height level 2
!
!   Latitude 1:  Linear interpolation
!
          fintx(i,k,1,2) = fb(ii1  ,kk,jj-1) * xl(i,k,1) +
     &                     fb(ii1+1,kk,jj-1) * xr(i,k,1)
!
!   Latitude 2:  Cubic interpolation
!
          fxl = (   - 2.*fb(ii2-1,kk,jj) 
     &              - 3.*fb(ii2  ,kk,jj) 
     &              + 6.*fb(ii2+1,kk,jj) 
     &              -    fb(ii2+2,kk,jj) )*rdlam6jj
          fxr = (        fb(ii2-1,kk,jj) 
     &              - 6.*fb(ii2  ,kk,jj) 
     &              + 3.*fb(ii2+1,kk,jj) 
     &              + 2.*fb(ii2+2,kk,jj) )*rdlam6jj
!
          deli = (       fb(ii2+1,kk,jj) - 
     &                   fb(ii2  ,kk,jj) )*rdlamjj
          tmp1 = fac*deli
          tmp2 = abs( tmp1 )
          if( deli*fxl   <= 0.0 ) fxl = 0.
          if( deli*fxr   <= 0.0 ) fxr = 0.
          if( abs( fxl ) > tmp2 ) fxl = tmp1
          if( abs( fxr ) > tmp2 ) fxr = tmp1
!
          fintx(i,k,2,2) = fb(ii2  ,kk,jj) * hl(i,k,2) 
     &                   + fb(ii2+1,kk,jj) * hr(i,k,2) 
     &                   + fxl*dhl(i,k,2) + fxr*dhr(i,k,2)
!
!   Latitude 3:  Cubic interpolation
!
          fxl = (   - 2.*fb(ii3-1,kk  ,jj+1) 
     &              - 3.*fb(ii3  ,kk  ,jj+1) 
     &              + 6.*fb(ii3+1,kk  ,jj+1) 
     &              -    fb(ii3+2,kk  ,jj+1) )*rdlam6jjp1
          fxr = (        fb(ii3-1,kk  ,jj+1) 
     &              - 6.*fb(ii3  ,kk  ,jj+1) 
     &              + 3.*fb(ii3+1,kk  ,jj+1) 
     &              + 2.*fb(ii3+2,kk  ,jj+1) )*rdlam6jjp1
!
          deli = (       fb(ii3+1,kk  ,jj+1) - 
     &                   fb(ii3  ,kk  ,jj+1) )*rdlamjjp1
          tmp1 = fac*deli
          tmp2 = abs( tmp1 )
          if( deli*fxl   <= 0.0 ) fxl = 0.
          if( deli*fxr   <= 0.0 ) fxr = 0.
          if( abs( fxl ) > tmp2 ) fxl = tmp1
          if( abs( fxr ) > tmp2 ) fxr = tmp1
!
          fintx(i,k,3,2) = fb(ii3  ,kk  ,jj+1) * hl(i,k,3) 
     &                   + fb(ii3+1,kk  ,jj+1) * hr(i,k,3) 
     &                   + fxl*dhl(i,k,3) + fxr*dhr(i,k,3)
!
!   Latitude 4:  Linear interpolation
!
          fintx(i,k,4,2) = fb(ii4  ,kk,jj+2) * xl(i,k,4) 
     &                   + fb(ii4+1,kk,jj+2) * xr(i,k,4)
!
! Height level 3
!
!   Latitude 1:  Linear interpolation
!
          fintx(i,k,1,3) = fb(ii1  ,kk+1,jj-1) * xl (i,k,1) 
     &                   + fb(ii1+1,kk+1,jj-1) * xr (i,k,1)
!
!   Latitude 2:  Cubic interpolation
!
          fxl = (   - 2.*fb(ii2-1,kk+1,jj  ) 
     &              - 3.*fb(ii2  ,kk+1,jj  ) 
     &              + 6.*fb(ii2+1,kk+1,jj  ) 
     &              -    fb(ii2+2,kk+1,jj  ) )*rdlam6jj
          fxr = (        fb(ii2-1,kk+1,jj  ) 
     &              - 6.*fb(ii2  ,kk+1,jj  ) 
     &              + 3.*fb(ii2+1,kk+1,jj  ) 
     &              + 2.*fb(ii2+2,kk+1,jj  ) )*rdlam6jj
!
          deli = (       fb(ii2+1,kk+1,jj  ) - 
     &                   fb(ii2  ,kk+1,jj  ) )*rdlamjj
          tmp1 = fac*deli
          tmp2 = abs( tmp1 )
          if( deli*fxl   <= 0.0 ) fxl = 0.
          if( deli*fxr   <= 0.0 ) fxr = 0.
          if( abs( fxl ) > tmp2 ) fxl = tmp1
          if( abs( fxr ) > tmp2 ) fxr = tmp1
!
          fintx(i,k,2,3) = fb(ii2  ,kk+1,jj ) * hl(i,k,2) 
     &                   + fb(ii2+1,kk+1,jj ) * hr(i,k,2) 
     &                   + fxl*dhl(i,k,2) + fxr*dhr(i,k,2)
!
!   Latitude 3:  Cubic interpolation
!
          fxl = (- 2.*fb(ii3-1,kk+1,jj+1) 
     &           - 3.*fb(ii3  ,kk+1,jj+1) 
     &           + 6.*fb(ii3+1,kk+1,jj+1) 
     &           -    fb(ii3+2,kk+1,jj+1) )*rdlam6jjp1
          fxr = (     fb(ii3-1,kk+1,jj+1) 
     &           - 6.*fb(ii3  ,kk+1,jj+1) 
     &           + 3.*fb(ii3+1,kk+1,jj+1) 
     &           + 2.*fb(ii3+2,kk+1,jj+1) )*rdlam6jjp1
!
          deli = (    fb(ii3+1,kk+1,jj+1) - 
     &                fb(ii3  ,kk+1,jj+1) )*rdlamjjp1
          tmp1 = fac*deli
          tmp2 = abs( tmp1 )
          if( deli*fxl  <= 0.0  ) fxl = 0.
          if( deli*fxr  <= 0.0  ) fxr = 0.
          if( abs( fxl ) > tmp2 ) fxl = tmp1
          if( abs( fxr ) > tmp2 ) fxr = tmp1
!
          fintx(i,k,3,3) = fb(ii3  ,kk+1,jj+1) * hl(i,k,3) 
     &                   + fb(ii3+1,kk+1,jj+1) * hr(i,k,3) 
     &                   + fxl*dhl(i,k,3) + fxr*dhr(i,k,3)
!
!   Latitude 4:  Linear interpolation
!
          fintx(i,k,4,3) = fb(ii4  ,kk+1,jj+2) * xl(i,k,4) 
     &                   + fb(ii4+1,kk+1,jj+2) * xr(i,k,4)
!
! Height level 4:  Linear interpolation on inner two latitudes only
!
!!!       fintx(i,k,1,4) = not used
          fintx(i,k,2,4) = fb(ii2  ,kk+2,jj  ) * xl(i,k,2) 
     &                   + fb(ii2+1,kk+2,jj  ) * xr(i,k,2)
          fintx(i,k,3,4) = fb(ii3  ,kk+2,jj+1) * xl(i,k,3) 
     &                   + fb(ii3+1,kk+2,jj+1) * xr(i,k,3)
!!!       fintx(i,k,4,4) = not used

! TOP interval
!
!The following loop computes x-derivatives for those cases when the
! departure point lies in either the top or bottom interval of the
! model grid.  In this special case, data are shifted up or down to
! keep the departure point in the middle interval of the 4-point
! stencil. Therefore, some derivatives that were computed above will
! be over-written.
          if(kdp (i,k) ==  1) then
!
! shift levels 4 and 2 data to levels 1 and 3, respectively
!
             fintx(i,k,2,1) = fintx(i,k,2,4)
             fintx(i,k,3,1) = fintx(i,k,3,4)
!
             fintx(i,k,1,3) = fintx(i,k,1,2)
             fintx(i,k,2,3) = fintx(i,k,2,2)
             fintx(i,k,3,3) = fintx(i,k,3,2)
             fintx(i,k,4,3) = fintx(i,k,4,2)
!
! Height level 1 (placed in level 2 of stencil):
!
!   Latitude 1:  Linear interpolation
!
!     print *,' i=',i,' k=',k,' ii1=',ii1,' jj=',jj,' in her_xyz'

             fintx(i,k,1,2) = fb(ii1  ,1,jj-1) * xl(i,k,1)
     &                      + fb(ii1+1,1,jj-1) * xr(i,k,1)
!
!   Latitude 2:  Cubic interpolation
!
             fxl = (- 2.*fb(ii2-1,1,jj) - 3.*fb(ii2  ,1,jj)
     &              + 6.*fb(ii2+1,1,jj) -    fb(ii2+2,1,jj) )*rdlam6jj

             fxr = (     fb(ii2-1,1,jj) - 6.*fb(ii2  ,1,jj)
     &              + 3.*fb(ii2+1,1,jj) + 2.*fb(ii2+2,1,jj) )*rdlam6jj
!
             deli = (    fb(ii2+1,1,jj) -    fb(ii2  ,1,jj) )*rdlamjj

             tmp1 = fac*deli
             tmp2 = abs( tmp1 )
             if( deli*fxl   <= 0.0  ) fxl = 0.
             if( deli*fxr   <= 0.0  ) fxr = 0.
             if( abs( fxl ) > tmp2 ) fxl = tmp1
             if( abs( fxr ) > tmp2 ) fxr = tmp1
!
             fintx(i,k,2,2) = fb(ii2  ,1,jj) * hl(i,k,2)
     &                      + fb(ii2+1,1,jj) * hr(i,k,2)
     &                      + fxl*dhl(i,k,2) + fxr*dhr(i,k,2)
!
!   Latitude 3:  Cubic interpolation
!
             fxl = (- 2.*fb(ii3-1,1,jj+1) - 3.*fb(ii3  ,1,jj+1)
     &              + 6.*fb(ii3+1,1,jj+1) -    fb(ii3+2,1,jj+1) )
     &                                    * rdlam6jjp1
             fxr = (  fb(ii3-1,1,jj+1) - 6.*fb(ii3  ,1,jj+1)
     &           + 3.*fb(ii3+1,1,jj+1) + 2.*fb(ii3+2,1,jj+1))*rdlam6jjp1
!
             deli = (fb (ii3+1,1,jj+1) - fb(ii3  ,1,jj+1))*rdlamjjp1
             tmp1 = fac*deli
             tmp2 = abs( tmp1 )
             if( deli*fxl   .le. 0.0  ) fxl = 0.
             if( deli*fxr   .le. 0.0  ) fxr = 0.
             if( abs( fxl ) > tmp2 ) fxl = tmp1
             if( abs( fxr ) > tmp2 ) fxr = tmp1
!
             fintx(i,k,3,2) = fb(ii3  ,1,jj+1) * hl(i,k,3)
     &                      + fb(ii3+1,1,jj+1) * hr(i,k,3)
     &                      + fxl*dhl(i,k,3) + fxr*dhr(i,k,3)
!
!   Latitude 4:  Linear interpolation
!
             fintx(i,k,4,2) = fb(ii4  ,1,jj+2) * xl(i,k,4)
     &                      + fb(ii4+1,1,jj+2) * xr(i,k,4)
!
! Height level 3 (placed in level 4 of stencil):
!  Linear interpolation on inner two latitudes only
!
!!!          fintx(i,k,1,4) = not used
             fintx(i,k,2,4) = fb(ii2  ,3,jj  ) * xl(i,k,2)
     &                      + fb(ii2+1,3,jj  ) * xr(i,k,2)
             fintx(i,k,3,4) = fb(ii3  ,3,jj+1) * xl(i,k,3)
     &                      + fb(ii3+1,3,jj+1) * xr(i,k,3)
!!!          fintx(i,k,4,4) = not used
!
! BOT interval
!
          elseif(kdp (i,k) == kdimm1) then
!
! shift levels 1 and 3 data to levels 4 and 2, respectively
!
             fintx(i,k,2,4) = fintx(i,k,2,1)
             fintx(i,k,3,4) = fintx(i,k,3,1)
!
             fintx(i,k,1,2) = fintx(i,k,1,3)
             fintx(i,k,2,2) = fintx(i,k,2,3)
             fintx(i,k,3,2) = fintx(i,k,3,3)
             fintx(i,k,4,2) = fintx(i,k,4,3)
!
! Height level 2 (placed in level 1 of stencil):
!  Linear interpolation on inner two latitudes only
!
!!!          fintx(i,k,1,1) =  not used
             fintx(i,k,2,1) = fb(ii2  ,kdimm2,jj  ) * xl (i,k,2)
     &                      + fb(ii2+1,kdimm2,jj  ) * xr (i,k,2)
             fintx(i,k,3,1) = fb(ii3  ,kdimm2,jj+1) * xl (i,k,3)
     &                      + fb(ii3+1,kdimm2,jj+1) * xr (i,k,3)
!!!          fintx(i,k,4,1) =  not used
!
! Height level 4 (placed in level 3 of stencil):
!
!   Latitude 1:  Linear interpolation
!
             fintx(i,k,1,3) = fb(ii1  ,kdim,jj-1) * xl (i,k,1)
     &                      + fb(ii1+1,kdim,jj-1) * xr (i,k,1)
!
!   Latitude 2:  Cubic interpolation
!
             fxl = (   - 2.*fb(ii2-1,kdim,jj  )
     &                 - 3.*fb(ii2  ,kdim,jj  )
     &                 + 6.*fb(ii2+1,kdim,jj  )
     &                 -    fb(ii2+2,kdim,jj  ) )*rdlam6jj
             fxr = (        fb(ii2-1,kdim,jj  )
     &                 - 6.*fb(ii2  ,kdim,jj  )
     &                 + 3.*fb(ii2+1,kdim,jj  )
     &                 + 2.*fb(ii2+2,kdim,jj  ) )*rdlam6jj
!
             deli = (       fb(ii2+1,kdim,jj  ) -
     &                      fb(ii2  ,kdim,jj  ) )*rdlamjj
             tmp1 = fac*deli
             tmp2 = abs( tmp1 )
             if( deli*fxl   <= 0.0 ) fxl = 0.
             if( deli*fxr   <= 0.0 ) fxr = 0.
             if( abs( fxl ) > tmp2 ) fxl = tmp1
             if( abs( fxr ) > tmp2 ) fxr = tmp1
!
                 fintx(i,k,2,3) = fb(ii2  ,kdim,jj  ) * hl (i,k,2)
     &                          + fb(ii2+1,kdim,jj  ) * hr (i,k,2)
     &                          + fxl*dhl(i,k,2) + fxr*dhr(i,k,2)
!
!   Latitude 3:  Cubic interpolation
!
             fxl = (   - 2.*fb(ii3-1,kdim,jj+1)
     &                 - 3.*fb(ii3  ,kdim,jj+1)
     &                 + 6.*fb(ii3+1,kdim,jj+1)
     &                 -    fb(ii3+2,kdim,jj+1) )*rdlam6jjp1
             fxr = (        fb(ii3-1,kdim,jj+1)
     &                 - 6.*fb(ii3  ,kdim,jj+1)
     &                 + 3.*fb(ii3+1,kdim,jj+1)
     &                 + 2.*fb(ii3+2,kdim,jj+1) )*rdlam6jjp1
!
             deli = (       fb(ii3+1,kdim,jj+1) -
     &                      fb(ii3  ,kdim,jj+1) )*rdlamjjp1
             tmp1 = fac*deli
             tmp2 = abs( tmp1 )
             if( deli*fxl   <= 0.0 ) fxl = 0.
             if( deli*fxr   <= 0.0 ) fxr = 0.
             if( abs( fxl ) > tmp2 ) fxl = tmp1
             if( abs( fxr ) > tmp2 ) fxr = tmp1
!
             fintx(i,k,3,3) = fb(ii3  ,kdim,jj+1) * hl(i,k,3)
     &                      + fb(ii3+1,kdim,jj+1) * hr(i,k,3)
     &                      + fxl*dhl(i,k,3) + fxr*dhr(i,k,3)
!
!   Latitude 4:  Linear interpolation
!
             fintx(i,k,4,3) = fb(ii4  ,kdim,jj+2) * xl(i,k,4)
     &                      + fb(ii4+1,kdim,jj+2) * xr(i,k,4)
          endif
        enddo
      enddo
!
!-----------------------------------------------------------------------
!
! y derivatives at the inner height levels (kk = 2,3) needed for
! z-interpolation
!
       do kk  = 2,3
         do k=1,plev
           do i=i_1,i_2
             jj = jdp(i,k)
             fbot(i,k,kk) = lbasdy(1,1,jj) * fintx(i,k,1,kk) 
     &                    + lbasdy(2,1,jj) * fintx(i,k,2,kk) 
     &                    + lbasdy(3,1,jj) * fintx(i,k,3,kk) 
     &                    + lbasdy(4,1,jj) * fintx(i,k,4,kk)
             ftop(i,k,kk) = lbasdy(1,2,jj) * fintx(i,k,1,kk) 
     &                    + lbasdy(2,2,jj) * fintx(i,k,2,kk) 
     &                    + lbasdy(3,2,jj) * fintx(i,k,3,kk) 
     &                    + lbasdy(4,2,jj) * fintx(i,k,4,kk)
!
! Apply SCM0 limiter to derivative estimates.
!
             deli = ( fintx(i,k,3,kk) - fintx(i,k,2,kk) )*rdphi(i,k)
             tmp1 = fac*deli
             tmp2 = abs( tmp1 )
             if( deli*fbot(i,k,kk)   <= 0.0 ) fbot(i,k,kk) = 0.
             if( deli*ftop(i,k,kk)   <= 0.0 ) ftop(i,k,kk) = 0.
             if( abs( fbot(i,k,kk) ) > tmp2 ) fbot(i,k,kk) = tmp1
             if( abs( ftop(i,k,kk) ) > tmp2 ) ftop(i,k,kk) = tmp1
           enddo
         enddo
       enddo
!
! PART 3:  y-interpolants
!
        do k=1,plev
           do i=i_1,i_2
              finty(i,k,1) = fintx(i,k,2,1)*ys (i,k) 
     &                     + fintx(i,k,3,1)*yn (i,k)

              finty(i,k,2) = fintx(i,k,2,2)*hs (i,k)
     &                     + fbot (i,k  ,2)*dhs(i,k)
     &                     + fintx(i,k,3,2)*hn (i,k)
     &                     + ftop (i,k  ,2)*dhn(i,k)

              finty(i,k,3) = fintx(i,k,2,3)*hs (i,k)
     &                     + fbot (i,k  ,3)*dhs(i,k)
     &                     + fintx(i,k,3,3)*hn (i,k)
     &                     + ftop (i,k  ,3)*dhn(i,k)

              finty(i,k,4) = fintx(i,k,2,4)*ys (i,k) 
     &                     + fintx(i,k,3,4)*yn (i,k)
           end do
        end do
!
!-----------------------------------------------------------------------
!
       do k = 1,plev
         do i = i_1,i_2
           kk = kdp (i,k)
           if(kk == 1) then

              ftop(i,k,1) = (finty(i,k,3) - finty(i,k,2))
     &                    / detam(1)
              fbot(i,k,1) = lbasdz(1,1,2) * finty(i,k,2)
     &                    + lbasdz(2,1,2) * finty(i,k,3)
     &                    + lbasdz(3,1,2) * finty(i,k,4)
     &                    + lbasdz(4,1,2) * finty(i,k,1)

           elseif(kk == kdimm1) then

              ftop(i,k,1) = lbasdz(1,2,kdimm2) * finty(i,k,4)
     &                    + lbasdz(2,2,kdimm2) * finty(i,k,1)
     &                    + lbasdz(3,2,kdimm2) * finty(i,k,2)
     &                    + lbasdz(4,2,kdimm2) * finty(i,k,3)
              fbot(i,k,1) = 0.0

           else

              ftop(i,k,1) = lbasdz(1,1,kk) * finty(i,k,1)
     &                    + lbasdz(2,1,kk) * finty(i,k,2)
     &                    + lbasdz(3,1,kk) * finty(i,k,3)
     &                    + lbasdz(4,1,kk) * finty(i,k,4)
              fbot(i,k,1) = lbasdz(1,2,kk) * finty(i,k,1)
     &                    + lbasdz(2,2,kk) * finty(i,k,2)
     &                    + lbasdz(3,2,kk) * finty(i,k,3)
     &                    + lbasdz(4,2,kk) * finty(i,k,4)

           endif
!
           deli = (finty(i,k,3) - finty(i,k,2)) / detam(kk)
           tmp1 = fac*deli
           tmp2 = abs(tmp1)

           if( deli*fbot(i,k,1)   <= 0.0 ) fbot(i,k,1) = 0.
           if( deli*ftop(i,k,1)   <= 0.0 ) ftop(i,k,1) = 0.
           if( abs( fbot(i,k,1) ) > tmp2 ) fbot(i,k,1) = tmp1
           if( abs( ftop(i,k,1) ) > tmp2 ) ftop(i,k,1) = tmp1

           fdp(i,plevp-k) = 
     &               finty(i,k,2)*ht(i,k) + ftop(i,k,1)*dht(i,k)
     &             + finty(i,k,3)*hb(i,k) + fbot(i,k,1)*dhb(i,k)


         end do
       enddo

!-----------------------------------------------------------------------
!     if (lprnt) then
!     print*,' fdp=', fdp(1,plevp-5:plev)
!     print*,'goodby from her_xyz_in_h '
!     endif
!-----------------------------------------------------------------------
      return
      end
      subroutine herxin_h(i_1,i_2,s_lat_in_pe,n_lat_in_pe,
     &                    idp,jdp,kdp,kkdp,xl,xr,hl,hr,dhl,dhr,fb,fintx)
!
      use      pmgrid        , only : plev,plond
      use slgshr             , only : rdlam,rdlam6
      implicit none
!
      integer i_1,i_2,lon_dim_h,s_lat_in_pe,n_lat_in_pe
      integer idp(i_1:i_2,plev,4),jdp(i_1:i_2,plev),kdp(i_1:i_2,plev)
      integer kkdp(i_1:i_2,plev),kdim,kdimm1,kdimm2
      real rdlamjj,rdlamjjp1,rdlam6jj,rdlam6jjp1
      real xl(i_1:i_2,plev,4), xr(i_1:i_2,plev,4),
     &     hl(i_1:i_2,plev,4), hr(i_1:i_2,plev,4),
     &    dhl(i_1:i_2,plev,4),dhr(i_1:i_2,plev,4)
      real fb(plond,plev,s_lat_in_pe:n_lat_in_pe)
      real fintx(i_1:i_2,plev,4,4)
! locals
      integer ii1,ii2,ii3,ii4,jj,kk,i,k
      real fxl,fxr,deli,tmp1,tmp2
      real fac
!
      fac  = 3.*(1. - 10.*epsilon(fac))
!
      kdim   = plev
      kdimm1 = kdim - 1
      kdimm2 = kdim - 2
!
! PART 1:  x-interpolation
!
! Loop over fields.
! ..x interpolation at each height needed for z interpolation.
! ...x interpolation at each latitude needed for y interpolation.
!
      do k=1,plev
        do i=i_1,i_2
          ii1 = idp(i,k,1)
          ii2 = idp(i,k,2)
          ii3 = idp(i,k,3)
          ii4 = idp(i,k,4)
          jj  = jdp(i,k)
          kk  = kkdp(i,k)
          rdlamjj   = rdlam(jj)
          rdlamjjp1 = rdlam(jj+1)
!
          rdlam6jj   = rdlam6(jj)
          rdlam6jjp1 = rdlam6(jj+1)

!sela     if (ii2-1 .le. 0) then
!sela       print*,' zero index in herxin ii',ii2
!sela       stop221
!sela     endif
!sela     if (jj-1 .le. 0) then
!sela        print*,' zero index in herxin jj',jj
!sela        stop222
!sela     endif
!sela     if (kk-1 .le. 0) then
!sela        print*,' zero index in herxin kk',kk
!sela        stop223
!sela     endif
!
!
! Height level 1:  Linear interpolation on inner two latitudes only
!
!!!       fintx(i,k,1,1) = not used
          fintx(i,k,2,1) = fb(ii2  ,kk-1,jj  ) * xl(i,k,2) +
     &                     fb(ii2+1,kk-1,jj  ) * xr(i,k,2)
          fintx(i,k,3,1) = fb(ii3  ,kk-1,jj+1) * xl(i,k,3) +
     &                   fb(ii3+1,kk-1,jj+1) * xr(i,k,3)
!!!       fintx(i,k,4,1) = not used
!
! Height level 2
!
!   Latitude 1:  Linear interpolation
!
          fintx(i,k,1,2) = fb(ii1  ,kk,jj-1) * xl(i,k,1) +
     &                     fb(ii1+1,kk,jj-1) * xr(i,k,1)
!
!   Latitude 2:  Cubic interpolation
!
          fxl = (- 2.*fb(ii2-1,kk,jj) - 3.*fb(ii2  ,kk,jj) 
     &           + 6.*fb(ii2+1,kk,jj) -    fb(ii2+2,kk,jj) )*rdlam6jj
          fxr = (     fb(ii2-1,kk,jj) - 6.*fb(ii2  ,kk,jj) 
     &           + 3.*fb(ii2+1,kk,jj) + 2.*fb(ii2+2,kk,jj) )*rdlam6jj
!
          deli = (    fb (ii2+1,kk,jj) -   fb (ii2  ,kk,jj) )*rdlamjj
          tmp1 = fac*deli
          tmp2 = abs( tmp1 )
          if( deli*fxl   <= 0.0 ) fxl = 0.
          if( deli*fxr   <= 0.0 ) fxr = 0.
          if( abs( fxl ) > tmp2 ) fxl = tmp1
          if( abs( fxr ) > tmp2 ) fxr = tmp1
!
          fintx(i,k,2,2) = fb(ii2  ,kk,jj) * hl (i,k,2) 
     &                   + fb(ii2+1,kk,jj) * hr (i,k,2) 
     &                   + fxl*dhl(i,k,2) + fxr*dhr(i,k,2)
!
!   Latitude 3:  Cubic interpolation
!
          fxl = (- 2.*fb(ii3-1,kk  ,jj+1)
     &           - 3.*fb(ii3  ,kk  ,jj+1) 
     &           + 6.*fb(ii3+1,kk  ,jj+1)
     &           -    fb(ii3+2,kk  ,jj+1) )*rdlam6jjp1
          fxr = (     fb(ii3-1,kk  ,jj+1) 
     &           - 6.*fb(ii3  ,kk  ,jj+1) 
     &           + 3.*fb(ii3+1,kk  ,jj+1) 
     &           + 2.*fb(ii3+2,kk  ,jj+1) )*rdlam6jjp1
!
          deli = (    fb(ii3+1,kk  ,jj+1) - 
     &                fb(ii3  ,kk  ,jj+1) )*rdlamjjp1
          tmp1 = fac*deli
          tmp2 = abs( tmp1 )
          if( deli*fxl   <= 0.0 ) fxl = 0.
          if( deli*fxr   <= 0.0 ) fxr = 0.
          if( abs( fxl ) > tmp2 ) fxl = tmp1
          if( abs( fxr ) > tmp2 ) fxr = tmp1
!
          fintx(i,k,3,2) = fb(ii3  ,kk  ,jj+1) * hl(i,k,3) 
     &                   + fb(ii3+1,kk  ,jj+1) * hr(i,k,3) 
     &                  + fxl*dhl(i,k,3) + fxr*dhr(i,k,3)
!
!   Latitude 4:  Linear interpolation
!
          fintx(i,k,4,2) = fb(ii4  ,kk,jj+2) * xl(i,k,4) 
     &                   + fb(ii4+1,kk,jj+2) * xr(i,k,4)
!
! Height level 3
!
!   Latitude 1:  Linear interpolation
!
          fintx(i,k,1,3) = fb(ii1  ,kk+1,jj-1) * xl(i,k,1) 
     &                   + fb(ii1+1,kk+1,jj-1) * xr(i,k,1)
!
!   Latitude 2:  Cubic interpolation
!
          fxl = (- 2.*fb(ii2-1,kk+1,jj  ) 
     &           - 3.*fb(ii2  ,kk+1,jj  ) 
     &           + 6.*fb(ii2+1,kk+1,jj  ) 
     &           -    fb(ii2+2,kk+1,jj  ) )*rdlam6jj
          fxr = (     fb(ii2-1,kk+1,jj  ) 
     &           - 6.*fb(ii2  ,kk+1,jj  ) 
     &           + 3.*fb(ii2+1,kk+1,jj  ) 
     &           + 2.*fb(ii2+2,kk+1,jj  ) )*rdlam6jj
!
          deli = (    fb(ii2+1,kk+1,jj  ) - 
     &                       fb (ii2  ,kk+1,jj  ) )*rdlamjj
          tmp1 = fac*deli
          tmp2 = abs( tmp1 )
          if( deli*fxl   <= 0.0 ) fxl = 0.
          if( deli*fxr   <= 0.0 ) fxr = 0.
          if( abs( fxl ) > tmp2 ) fxl = tmp1
          if( abs( fxr ) > tmp2 ) fxr = tmp1
!
          fintx(i,k,2,3) = fb(ii2  ,kk+1,jj  ) * hl (i,k,2) 
     &                   + fb(ii2+1,kk+1,jj  ) * hr (i,k,2) 
     &                   + fxl*dhl(i,k,2) + fxr*dhr(i,k,2)
!
!   Latitude 3:  Cubic interpolation
!
          fxl = (   - 2.*fb(ii3-1,kk+1,jj+1) 
     &              - 3.*fb(ii3  ,kk+1,jj+1) 
     &              + 6.*fb(ii3+1,kk+1,jj+1) 
     &              -    fb(ii3+2,kk+1,jj+1) )*rdlam6jjp1
          fxr = (        fb(ii3-1,kk+1,jj+1) 
     &              - 6.*fb(ii3  ,kk+1,jj+1) 
     &              + 3.*fb(ii3+1,kk+1,jj+1) 
     &              + 2.*fb(ii3+2,kk+1,jj+1) )*rdlam6jjp1
!
          deli = (       fb (ii3+1,kk+1,jj+1) - 
     &                   fb (ii3  ,kk+1,jj+1) )*rdlamjjp1
          tmp1 = fac*deli
          tmp2 = abs( tmp1 )
          if( deli*fxl   <= 0.0 ) fxl = 0.
          if( deli*fxr   <= 0.0 ) fxr = 0.
          if( abs( fxl ) > tmp2 ) fxl = tmp1
          if( abs( fxr ) > tmp2 ) fxr = tmp1

          fintx(i,k,3,3) = fb(ii3  ,kk+1,jj+1) * hl (i,k,3) 
     &                   + fb(ii3+1,kk+1,jj+1) * hr (i,k,3) 
     &                   + fxl*dhl(i,k,3) + fxr*dhr(i,k,3)

!   Latitude 4:  Linear interpolation
!
          fintx(i,k,4,3) = fb(ii4  ,kk+1,jj+2) * xl (i,k,4) 
     &                   + fb(ii4+1,kk+1,jj+2) * xr (i,k,4)
!
! Height level 4:  Linear interpolation on inner two latitudes only
!
!!!       fintx(i,k,1,4) = not used
          fintx(i,k,2,4) = fb(ii2  ,kk+2,jj  ) * xl (i,k,2) 
     &                   + fb(ii2+1,kk+2,jj  ) * xr (i,k,2)
          fintx(i,k,3,4) = fb(ii3  ,kk+2,jj+1) * xl (i,k,3) 
     &                   + fb(ii3+1,kk+2,jj+1) * xr (i,k,3)
!!!       fintx(i,k,4,4) = not used 
!
! The following loop computes x-derivatives for those cases when the
! departure point lies in either the top or bottom interval of the
! model grid.  In this special case, data are shifted up or down to
! keep the departure point in the middle interval of the 4-point
! stencil. Therefore, some derivatives that were computed above will
! be over-written.

! TOP interval
!
          if(kdp (i,k) .eq. 1) then
!
! shift levels 4 and 2 data to levels 1 and 3, respectively
!
            fintx(i,k,2,1) = fintx(i,k,2,4)
            fintx(i,k,3,1) = fintx(i,k,3,4)
!
            fintx(i,k,1,3) = fintx(i,k,1,2)
            fintx(i,k,2,3) = fintx(i,k,2,2)
            fintx(i,k,3,3) = fintx(i,k,3,2)
            fintx(i,k,4,3) = fintx(i,k,4,2)
!
! Height level 1 (placed in level 2 of stencil):
!
!   Latitude 1:  Linear interpolation
!
            fintx(i,k,1,2) = fb (ii1  ,1,jj-1)*xl (i,k,1)
     &                     + fb (ii1+1,1,jj-1)*xr (i,k,1)
!
!   Latitude 2:  Cubic interpolation
!
            fxl = (   - 2.*fb (ii2-1,1,jj  )
     &                - 3.*fb (ii2  ,1,jj  )
     &                + 6.*fb (ii2+1,1,jj  )
     &                -    fb (ii2+2,1,jj  ) )*rdlam6jj
            fxr = (        fb (ii2-1,1,jj  )
     &                - 6.*fb (ii2  ,1,jj  )
     &                + 3.*fb (ii2+1,1,jj  )
     &                + 2.*fb (ii2+2,1,jj  ) )*rdlam6jj
!
            deli = (       fb (ii2+1,1,jj  ) -
     &                      fb (ii2  ,1,jj  ) )*rdlamjj
            tmp1 = fac*deli
            tmp2 = abs( tmp1 )
            if( deli*fxl   <= 0.0 ) fxl = 0.
            if( deli*fxr   <= 0.0 ) fxr = 0.
            if( abs( fxl ) > tmp2 ) fxl = tmp1
            if( abs( fxr ) > tmp2 ) fxr = tmp1
!
            fintx(i,k,2,2) = fb (ii2  ,1,jj  )*hl (i,k,2)
     &                     + fb (ii2+1,1,jj  )*hr (i,k,2)
     &                     + fxl*dhl(i,k,2) + fxr*dhr(i,k,2)

!   Latitude 3:  Cubic interpolation
!
            fxl = (   - 2.*fb (ii3-1,1,jj+1)
     &                - 3.*fb (ii3  ,1,jj+1)
     &                + 6.*fb (ii3+1,1,jj+1)
     &                -    fb (ii3+2,1,jj+1) )*rdlam6jjp1
            fxr = (        fb (ii3-1,1,jj+1)
     &                - 6.*fb (ii3  ,1,jj+1)
     &                + 3.*fb (ii3+1,1,jj+1)
     &                + 2.*fb (ii3+2,1,jj+1) )*rdlam6jjp1

            deli = (       fb (ii3+1,1,jj+1) -
     &                     fb (ii3  ,1,jj+1) )*rdlamjjp1
            tmp1 = fac*deli
            tmp2 = abs( tmp1 )
            if( deli*fxl   <= 0.0 ) fxl = 0.
            if( deli*fxr   <= 0.0 ) fxr = 0.
            if( abs( fxl ) > tmp2 ) fxl = tmp1
            if( abs( fxr ) > tmp2 ) fxr = tmp1
!
            fintx(i,k,3,2) = fb (ii3  ,1,jj+1)*hl (i,k,3)
     &                     + fb (ii3+1,1,jj+1)*hr (i,k,3)
     &                     + fxl*dhl(i,k,3) + fxr*dhr(i,k,3)
!
!   Latitude 4:  Linear interpolation
!
            fintx(i,k,4,2) = fb (ii4  ,1,jj+2)*xl (i,k,4)
     &                     + fb (ii4+1,1,jj+2)*xr (i,k,4)

! Height level 3 (placed in level 4 of stencil):
!  Linear interpolation on inner two latitudes only
!
!!!         fintx(i,k,1,4) = not used
            fintx(i,k,2,4) = fb (ii2  ,3,jj  )*xl (i,k,2)
     &                     + fb (ii2+1,3,jj  )*xr (i,k,2)
            fintx(i,k,3,4) = fb (ii3  ,3,jj+1)*xl (i,k,3)
     &                     + fb (ii3+1,3,jj+1)*xr (i,k,3)
!!!         fintx(i,k,4,4) = not used

! BOT interval
!
          else if(kdp (i,k) .eq. kdimm1) then
!
! shift levels 1 and 3 data to levels 4 and 2, respectively
!
            fintx(i,k,2,4) = fintx(i,k,2,1)
            fintx(i,k,3,4) = fintx(i,k,3,1)
!
            fintx(i,k,1,2) = fintx(i,k,1,3)
            fintx(i,k,2,2) = fintx(i,k,2,3)
            fintx(i,k,3,2) = fintx(i,k,3,3)
            fintx(i,k,4,2) = fintx(i,k,4,3)

! Height level 2 (placed in level 1 of stencil):
!  Linear interpolation on inner two latitudes only
!
!!!         fintx(i,k,1,1) =  not used
            fintx(i,k,2,1) = fb (ii2  ,kdimm2,jj  )*xl (i,k,2)
     &                     + fb (ii2+1,kdimm2,jj  )*xr (i,k,2)
            fintx(i,k,3,1) = fb (ii3  ,kdimm2,jj+1)*xl (i,k,3)
     &                     + fb (ii3+1,kdimm2,jj+1)*xr (i,k,3)
!!!         fintx(i,k,4,1) =  not used
!
! Height level 4 (placed in level 3 of stencil):
!
!   Latitude 1:  Linear interpolation
!
            fintx(i,k,1,3) = fb (ii1  ,kdim,jj-1)*xl (i,k,1)
     &                     + fb (ii1+1,kdim,jj-1)*xr (i,k,1)
!
!   Latitude 2:  Cubic interpolation
!
            fxl = (   - 2.*fb (ii2-1,kdim,jj  )
     &                - 3.*fb (ii2  ,kdim,jj  )
     &                + 6.*fb (ii2+1,kdim,jj  )
     &                -    fb (ii2+2,kdim,jj  ) )*rdlam6jj
            fxr = (        fb (ii2-1,kdim,jj  )
     &                - 6.*fb (ii2  ,kdim,jj  )
     &                + 3.*fb (ii2+1,kdim,jj  )
     &                + 2.*fb (ii2+2,kdim,jj  ) )*rdlam6jj
!
            deli = (       fb (ii2+1,kdim,jj  ) -
     &                     fb (ii2  ,kdim,jj  ) )*rdlamjj
            tmp1 = fac*deli
            tmp2 = abs( tmp1 )
            if( deli*fxl   <= 0.0 ) fxl = 0.
            if( deli*fxr   <= 0.0 ) fxr = 0.
            if( abs( fxl ) > tmp2 ) fxl = tmp1
            if( abs( fxr ) > tmp2 ) fxr = tmp1

            fintx(i,k,2,3) = fb (ii2  ,kdim,jj  )*hl (i,k,2)
     &                     + fb (ii2+1,kdim,jj  )*hr (i,k,2)
     &                     + fxl*dhl(i,k,2) + fxr*dhr(i,k,2)

!   Latitude 3:  Cubic interpolation
!
            fxl = (   - 2.*fb (ii3-1,kdim,jj+1)
     &                - 3.*fb (ii3  ,kdim,jj+1)
     &                + 6.*fb (ii3+1,kdim,jj+1)
     &                -    fb (ii3+2,kdim,jj+1) )*rdlam6jjp1
            fxr = (        fb (ii3-1,kdim,jj+1)
     &                - 6.*fb (ii3  ,kdim,jj+1)
     &                + 3.*fb (ii3+1,kdim,jj+1)
     &                + 2.*fb (ii3+2,kdim,jj+1) )*rdlam6jjp1
!
            deli = (       fb (ii3+1,kdim,jj+1) -
     &                     fb (ii3  ,kdim,jj+1) )*rdlamjjp1
            tmp1 = fac*deli
            tmp2 = abs( tmp1 )
            if( deli*fxl   <= 0.0 ) fxl = 0.
            if( deli*fxr   <= 0.0 ) fxr = 0.
            if( abs( fxl ) > tmp2 ) fxl = tmp1
            if( abs( fxr ) > tmp2 ) fxr = tmp1
!
            fintx(i,k,3,3) = fb (ii3  ,kdim,jj+1)*hl (i,k,3)
     &                     + fb (ii3+1,kdim,jj+1)*hr (i,k,3)
     &                     + fxl*dhl(i,k,3) + fxr*dhr(i,k,3)
!
!   Latitude 4:  Linear interpolation
!
            fintx(i,k,4,3) = fb (ii4  ,kdim,jj+2)*xl (i,k,4)
     &                     + fb (ii4+1,kdim,jj+2)*xr (i,k,4)
          endif
        end do
      end do
!
      return
      end
      subroutine heryin_h(i_1,i_2,jdp,ys,yn,hs,hn,dhs,dhn,rdphi,
     &                    fintx,finty)
!
      use      pmgrid        , only : plev
      use      slgshr        , only : lbasdy
      implicit none

      integer i_1,i_2,
     &        jdp(i_1:i_2,plev)    ! j-index of departure point coord.
      real    ys(i_1:i_2,plev),
     &        yn(i_1:i_2,plev),
     &        hs(i_1:i_2,plev),
     &        hn(i_1:i_2,plev),
     &        dhs(i_1:i_2,plev),
     &        dhn(i_1:i_2,plev),
     &        rdphi(i_1:i_2,plev)
      real   fintx(i_1:i_2,plev,4,4),finty(i_1:i_2,plev,4)
! locals
      integer i,k,kk,jj
      real    fac,deli,tmp1,tmp2
      real    ftop(i_1:i_2,plev,4),fbot(i_1:i_2,plev,4)
!
      fac  = 3.*(1. - 10.*epsilon(fac))           

!
! y derivatives at the inner height levels (kk = 2,3) needed for
! z-interpolation
!
      do k=1,plev
        do kk  = 2,3
          do i=i_1,i_2
            jj = jdp(i,k)
            fbot(i,k,kk) = lbasdy(1,1,jj) * fintx(i,k,1,kk) 
     &                   + lbasdy(2,1,jj) * fintx(i,k,2,kk) 
     &                   + lbasdy(3,1,jj) * fintx(i,k,3,kk) 
     &                   + lbasdy(4,1,jj) * fintx(i,k,4,kk)
            ftop(i,k,kk) = lbasdy(1,2,jj) * fintx(i,k,1,kk) 
     &                   + lbasdy(2,2,jj) * fintx(i,k,2,kk) 
     &                   + lbasdy(3,2,jj) * fintx(i,k,3,kk) 
     &                   + lbasdy(4,2,jj) * fintx(i,k,4,kk)
          enddo
        enddo
      enddo
!
! Apply SCM0 limiter to derivative estimates.
!
      do kk  = 2,3
        do k=1,plev
          do i=i_1,i_2
            deli = ( fintx(i,k,3,kk) - fintx(i,k,2,kk) )*rdphi(i,k)
            tmp1 = fac*deli
            tmp2 = abs( tmp1 )
            if( deli*fbot(i,k,kk)   <= 0.0 ) fbot(i,k,kk) = 0.
            if( deli*ftop(i,k,kk)   <= 0.0 ) ftop(i,k,kk) = 0.
            if( abs( fbot(i,k,kk) ) > tmp2 ) fbot(i,k,kk) = tmp1
            if( abs( ftop(i,k,kk) ) > tmp2 ) ftop(i,k,kk) = tmp1
          enddo
        enddo
      enddo
!
! PART 3:  y-interpolants
!
      do k=1,plev
        do i=i_1,i_2
          finty(i,k,1) = fintx(i,k,2,1) * ys (i,k) 
     &                 + fintx(i,k,3,1) * yn (i,k)

          finty(i,k,2) = fintx(i,k,2,2) * hs (i,k)
     &                 + fbot (i,k  ,2) * dhs(i,k)
     &                 + fintx(i,k,3,2) * hn (i,k)
     &                 + ftop (i,k  ,2) * dhn(i,k)

          finty(i,k,3) = fintx(i,k,2,3) * hs (i,k)
     &                 + fbot (i,k  ,3) * dhs(i,k)
     &                 + fintx(i,k,3,3) * hn (i,k)
     &                 + ftop (i,k  ,3) * dhn(i,k)

          finty(i,k,4) = fintx(i,k,2,4) * ys (i,k) 
     &                 + fintx(i,k,3,4) * yn (i,k)
        enddo
      enddo

      return
      end
      subroutine herzin_h(i_1,i_2,
     &                    kdim,finty,lbasdz,kdp,hb,ht,dhb,dht,detam,
     &                    term1,term2,term3,term4,fdp)
!
      use      pmgrid        , only : plev,plevp,plon
      implicit none
!
      integer i_1,i_2
      integer kdim,kdp(i_1:i_2,plev)
      real    finty(i_1:i_2,plev,4),
     &        lbasdz(4,2,kdim),fdp(plon,plev),
     &        ftop(i_1:i_2,plev,4),fbot(i_1:i_2,plev,4)
      real    ht(i_1:i_2,plev), hb(i_1:i_2,plev),
     &        dht(i_1:i_2,plev),dhb(i_1:i_2,plev)
      real    term1(i_1:i_2,plev),    ! |
     &        term2(i_1:i_2,plev),    ! |
     &        term3(i_1:i_2,plev),    ! |
     &        term4(i_1:i_2,plev)     ! |
      real    detam(plev)  ! delta eta at levels
      real    deli,tmp1,tmp2,fac
      integer i,k,m,kdimm1,kdimm2,kk
!
      kdimm1 = kdim - 1
      kdimm2 = kdim - 2
!
      do k = 1,plev
        do i = i_1,i_2
          kk = kdp (i,k)
          if(kk  == 1) then

             ftop(i,k,1) = (finty(i,k,3) - finty(i,k,2)) / detam(kk)
             fbot(i,k,1) = lbasdz(1,1,2) * finty(i,k,2)
     &                   + lbasdz(2,1,2) * finty(i,k,3)
     &                   + lbasdz(3,1,2) * finty(i,k,4)
     &                   + lbasdz(4,1,2) * finty(i,k,1)
          elseif(kk == kdimm1) then
             ftop(i,k,1) = lbasdz(1,2,kdimm2) * finty(i,k,4)
     &                   + lbasdz(2,2,kdimm2) * finty(i,k,1)
     &                   + lbasdz(3,2,kdimm2) * finty(i,k,2)
     &                   + lbasdz(4,2,kdimm2) * finty(i,k,3)
             fbot(i,k,1) = 0.0
          else
             ftop(i,k,1) = lbasdz(1,1,kk)*finty(i,k,1)
     &                   + lbasdz(2,1,kk)*finty(i,k,2)
     &                   + lbasdz(3,1,kk)*finty(i,k,3)
     &                   + lbasdz(4,1,kk)*finty(i,k,4)
             fbot(i,k,1) = lbasdz(1,2,kk)*finty(i,k,1)
     &                   + lbasdz(2,2,kk)*finty(i,k,2)
     &                   + lbasdz(3,2,kk)*finty(i,k,3)
     &                   + lbasdz(4,2,kk)*finty(i,k,4)

          endif
        enddo
      enddo

      fac  = 3.*(1. - 10.*epsilon(fac))
!     print*,' herzin_h fac = ',fac
      do k=1,plev
        do i = i_1,i_2
          deli = (finty(i,k,3) - finty(i,k,2)) / detam(kdp(i,k))
          tmp1 = fac*deli
          tmp2 = abs(tmp1)
          if( deli*fbot(i,k,1)   <= 0.0 ) fbot(i,k,1) = 0.
          if( deli*ftop(i,k,1)   <= 0.0 ) ftop(i,k,1) = 0.
          if( abs( fbot(i,k,1) ) > tmp2 ) fbot(i,k,1) = tmp1
          if( abs( ftop(i,k,1) ) > tmp2 ) ftop(i,k,1) = tmp1

          fdp(i,plevp-k) = finty(i,k,2)*ht(i,k) + ftop(i,k,1)*dht(i,k)
     &                   + finty(i,k,3)*hb(i,k) + fbot(i,k,1)*dhb(i,k)
        end do
      enddo

      return
      end
      subroutine lin_xyz_in_h(i_1,i_2,s_lat_in_pe,n_lat_in_pe,
     &                        fb1,fb2,fb3,idp,jdp,kdp,xl,xr,ys,yn,zb,zt,
     &                        fdp1,fdp2,fdp3)
!my
      use      pmgrid        , only : plev,plevp,plon,plond
      implicit none
      integer i_1,i_2,s_lat_in_pe,n_lat_in_pe
      integer idp(i_1:i_2,plev,4), jdp(i_1:i_2,plev),
     &        kdp(i_1:i_2,plev)
      real fb1(plond,plev ,s_lat_in_pe:n_lat_in_pe),
     &     fb2(plond,plev ,s_lat_in_pe:n_lat_in_pe),
     &     fb3(plond,plev ,s_lat_in_pe:n_lat_in_pe),
     &     fdp1(plon,plev),
     &     fdp2(plon,plev),
     &     fdp3(plon,plev)
      integer i,k,ii2,ii3,jj,kk,kr
      real xl(i_1:i_2,plev,4),xr(i_1:i_2,plev,4)
      real ys(i_1:i_2,plev),yn(i_1:i_2,plev)
      real zt(i_1:i_2,plev),zb(i_1:i_2,plev)
!
!my  assumes lagrange or hermite interpolation performed on ug,vg,tg and 
!my  now add in the nonliner terms contribution using trilinear interpolation

      do k=1,plev
        kr = plevp - k
        do i=i_1,i_2
          ii2 = idp(i,k,2)
          ii3 = idp(i,k,3)
          kk  = kdp(i,k)
          jj  = jdp(i,k)

          fdp1(i,kr) = fdp1(i,kr)
     &               + (( fb1 (ii2  ,kk  ,jj  ) * xl (i,k,2)
     &               +    fb1 (ii2+1,kk  ,jj  ) * xr (i,k,2) )*ys(i,k)
     &               +  ( fb1 (ii3  ,kk  ,jj+1) * xl (i,k,3)
     &               +    fb1 (ii3+1,kk  ,jj+1) * xr (i,k,3) )*yn(i,k) )
     &                                                        *zt(i,k)
     &               + (( fb1 (ii2  ,kk+1,jj  ) * xl (i,k,2)
     &               +    fb1 (ii2+1,kk+1,jj  ) * xr (i,k,2) )*ys(i,k)
     &               +  ( fb1 (ii3  ,kk+1,jj+1) * xl (i,k,3)
     &               +    fb1 (ii3+1,kk+1,jj+1) * xr (i,k,3) )*yn(i,k) )
     &                                                        *zb(i,k)
!$$$       print*,
!$$$     . ' in lin_xyz_in_h i,k,i2,i3,kdp,jdp,xl2,xl3,xr,ys,yn,zt,zb = ',
!$$$     . i,k,idp(i,k,2),idp(i,k,3),kdp(i,k),jdp(i,k),xl(i,k,2),xl(i,k,3),
!$$$     . ys(i,k),yn(i,k),zt(i,k),zb(i,k)

          fdp2(i,kr) = fdp2(i,kr)
     &               + (( fb2 (ii2  ,kk  ,jj  )*xl (i,k,2)
     &               +    fb2 (ii2+1,kk  ,jj  )*xr (i,k,2) )*ys(i,k)
     &               +  ( fb2 (ii3  ,kk  ,jj+1)*xl (i,k,3)
     &               +    fb2 (ii3+1,kk  ,jj+1)*xr (i,k,3) )*yn(i,k) )
     &                                                        *zt(i,k)
     &               + (( fb2 (ii2  ,kk+1,jj  ) * xl (i,k,2)
     &               +    fb2 (ii2+1,kk+1,jj  ) * xr (i,k,2) )*ys(i,k)
     &               +  ( fb2 (ii3  ,kk+1,jj+1) * xl (i,k,3)
     &               +    fb2 (ii3+1,kk+1,jj+1) * xr (i,k,3) )*yn(i,k) )
     &                                                        *zb(i,k)

          fdp3(i,kr) = fdp3(i,kr)
     &               + (( fb3 (ii2  ,kk  ,jj  ) * xl (i,k,2)
     &               +    fb3 (ii2+1,kk  ,jj  ) * xr (i,k,2) )*ys(i,k)
     &               +  ( fb3 (ii3  ,kk  ,jj+1) * xl (i,k,3)
     &               +    fb3 (ii3+1,kk  ,jj+1) * xr (i,k,3) )*yn(i,k) )
     &                                                        *zt(i,k)
     &               + (( fb3 (ii2  ,kk+1,jj  ) * xl (i,k,2)
     &               +    fb3 (ii2+1,kk+1,jj  ) * xr (i,k,2) )*ys(i,k)
     &               +  ( fb3 (ii3  ,kk+1,jj+1) * xl (i,k,3)
     &               +    fb3 (ii3+1,kk+1,jj+1) * xr (i,k,3) )*yn(i,k) )
     &                                                        *zb(i,k)
        enddo
      enddo

      return
      end
      subroutine quasim(f1,f2,fint)
!sk   quasi-monotone correction of fint 
      implicit none
      real f1,f2,fint,fmin,fmax, f(2)
      f(1) = f1
      f(2) = f2
      fmax = maxval(f)
      fmin = minval(f)
      if (fint.lt.fmin) fint = fmin
      if (fint.gt.fmax) fint = fmax
      return
      end
      subroutine wgt_cub_lin_xyz_in_h(i_1,i_2,s_lat_in_pe,n_lat_in_pe,
     &                                fb1,idp,jdp,kdp,xl,xr,ys,yn,zb,zt,
     &                                fdp1)
      use      pmgrid        , only : plev,plevp,plon,plond,wgt
      implicit none
      integer i_1,i_2,s_lat_in_pe,n_lat_in_pe
      integer idp(i_1:i_2,plev,4), jdp(i_1:i_2,plev),
     &        kdp(i_1:i_2,plev)
      real fb1(plond,plev ,s_lat_in_pe:n_lat_in_pe),
     &     fdp1(plon,plev)
      integer i,k,ii2,ii3,jj,kk,kr
      real xl(i_1:i_2,plev,4),xr(i_1:i_2,plev,4)
      real ys(i_1:i_2,plev),yn(i_1:i_2,plev)
      real zt(i_1:i_2,plev),zb(i_1:i_2,plev)
!     real wgt(plev)
!     real alphau,alphav,alphat
!     parameter (alphau=0.3,alphav=0.3,alphat=0.3)
!
!my  copy of lin_xyz_in_h for one variable with weights added
!my  
!     do k = 1,34
!       wgt(k) = 0.1
!     enddo
!     do k = 35,43
!       wgt(k) = 0.5
!     enddo
!     do k = 44,plev
!       wgt(k) = 1.0
!     enddo

      do k=1,plev
        kr = plevp - k
        do i=i_1,i_2
          ii2 = idp(i,k,2)
          ii3 = idp(i,k,3)
          kk  = kdp(i,k)
          jj  = jdp(i,k)

          fdp1(i,kr) = wgt(k)       * fdp1(i,kr)
     &               + (1.0-wgt(k)) * (
     &    (( fb1 (ii2  ,kk  ,jj  ) * xl (i,k,2)
     &  +    fb1 (ii2+1,kk  ,jj  ) * xr (i,k,2) )*ys(i,k)
     &  +  ( fb1 (ii3  ,kk  ,jj+1) * xl (i,k,3)
     &  +    fb1 (ii3+1,kk  ,jj+1) * xr (i,k,3) )*yn(i,k) )
     &                                                        *zt(i,k)
     &  + (( fb1 (ii2  ,kk+1,jj  ) * xl (i,k,2)
     &  +    fb1 (ii2+1,kk+1,jj  ) * xr (i,k,2) )*ys(i,k)
     &  +  ( fb1 (ii3  ,kk+1,jj+1) * xl (i,k,3)
     &  +    fb1 (ii3+1,kk+1,jj+1) * xr (i,k,3) )*yn(i,k) )
     &                                                        *zb(i,k)
     &                                            )
        enddo
      enddo
      return
      end
      subroutine lagzin_h_m(i_1,i_2,
     &                  kdim,finty,lbasdz,kdp,hb,ht,dhb,dht,detam,
     &                  term1,term2,term3,term4,fdp)

!cmy renamed with _m
!SK linear interpolation at the top and bottom
!sk 2/25/2012
!sk modified original subroutine for quasi-monotone
!sk quasi-cubic Lagrange interpolation
!sk
      use      pmgrid        , only : plev,plevp,plon,quamon
      implicit none

      integer i_1,i_2
      integer kdim,kdp(i_1:i_2,plev)
      real       finty(i_1:i_2,plev,4),
     &     lbasdz(4,2,kdim),fdp(plon,plev)
      integer i,k,m,kdimm1,kdimm2
      real  ht(i_1:i_2,plev), hb(i_1:i_2,plev),
     &     dht(i_1:i_2,plev),dhb(i_1:i_2,plev),
     &     ftop,fbot
      real term1(i_1:i_2,plev),    ! |
     &     term2(i_1:i_2,plev),    ! |
     &     term3(i_1:i_2,plev),    ! |
     &     term4(i_1:i_2,plev)     ! |
      real detam(plev)  ! delta eta at levels
      real fmax(i_1:i_2), fmin(i_1:i_2)
!
      kdimm1 = kdim - 1
         do k = 1,plev
            do i = i_1,i_2
               if(kdp (i,k) .eq. 1) then
                  fdp(i,plevp-k) = finty(i,k,2)*ht (i,k)
     &                           + finty(i,k,3)*hb (i,k)
               elseif(kdp (i,k) .eq. kdimm1) then
                  fdp(i,plevp-k) = finty(i,k,2)*ht (i,k)
     &                           + finty(i,k,3)*hb (i,k)
               else
               fdp(i,plevp-k) =   finty(i,k,1)*term1(i,k)
     &                          + finty(i,k,2)*term2(i,k)
     &                          + finty(i,k,3)*term3(i,k)
     &                          + finty(i,k,4)*term4(i,k)
      if (quamon) then
        fmax(i) = max(finty(i,k,1), finty(i,k,2),
     &                finty(i,k,3), finty(i,k,4))
        fmin(i) = min(finty(i,k,1), finty(i,k,2),
     &                finty(i,k,3), finty(i,k,4))
        fdp(i,plev-k) = max(fmin(i), min(fmax(i),fdp(i,plev-k)))
!       call quasim(finty(i,k,2),finty(i,k,3),fdp(i,plev-k))
      endif
               endif
            end do
         end do
      return
      end
      subroutine lagzinit_h_m(i_1,i_2,kdim,lbasiz,
     & etamid,detam,sigdp,kdp,kkdp,
     & hb,ht,dhb,dht,term1z,term2z,term3z,term4z)
!SK linear-interpolation weights if dp falls between layers 1 and 2
!   or layers K-1 and K.
      use      pmgrid        , only : plev
      implicit none
      integer i_1,i_2
      integer kdim,kdp(i_1:i_2,plev),kkdp(i_1:i_2,plev)
      real lbasiz(4,2,kdim),
     &  etamid(kdim),detam(kdim),
     &  sigdp(i_1:i_2,plev),fdp(i_1:i_2,plev)
      integer i,k,m,kdimm1,kdimm2
      real dzk,zt,zb,
     &     ht(i_1:i_2,plev), hb(i_1:i_2,plev),
     &    dht(i_1:i_2,plev),dhb(i_1:i_2,plev)
      real zmz1,                 ! |
     &     zmz2,                 ! |
     &     zmz3,                 ! |
     &     zmz4,                 ! |
     &     coef12,               ! |
     &     coef34,               ! | -- interpolation weights/coeffs.
     &     term1z(i_1:i_2,plev), ! |
     &     term2z(i_1:i_2,plev), ! |
     &     term3z(i_1:i_2,plev), ! |
     &     term4z(i_1:i_2,plev)  ! |
         do k=1,plev
         do i=i_1,i_2
            dzk  =  detam(kdp(i,k))
            zt   = ( etamid(kdp(i,k)+1) - sigdp(i,k) )/dzk
            zb   = 1. - zt
            ht(i,k)  = ( 3.0 - 2.0*zt )*zt**2
            hb(i,k)  = ( 3.0 - 2.0*zb )*zb**2
            dht(i,k)  = -dzk*( zt - 1. )*zt**2
            dhb(i,k)  =  dzk*( zb - 1. )*zb**2
!SKar 12/23/2011
	    ht(i,k)  = zt
            hb(i,k)  = zb
!SKar
            zmz1       = sigdp(i,k) -  lbasiz(1,1,kkdp(i,k))
            zmz2       = sigdp(i,k) -  lbasiz(2,1,kkdp(i,k))
            zmz3       = sigdp(i,k) -  lbasiz(3,1,kkdp(i,k))
            zmz4       = sigdp(i,k) -  lbasiz(4,1,kkdp(i,k))
            coef12     = zmz3*zmz4
            coef34     = zmz1*zmz2
            term1z(i,k) = coef12*zmz2*lbasiz(1,2,kkdp(i,k))
            term2z(i,k) = coef12*zmz1*lbasiz(2,2,kkdp(i,k))
            term3z(i,k) = coef34*zmz4*lbasiz(3,2,kkdp(i,k))
            term4z(i,k) = coef34*zmz3*lbasiz(4,2,kkdp(i,k))
         end do
         end do
      return
      end
