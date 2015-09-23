      subroutine gfidi_hyb_resonan
     & (lon_dim_a,lon_dim_h,lons_lat,lat,
     &  dg,tg,ug,vg,ug_m,vg_m,ww,rgs,dqdphi,dqdlam,qg,
     &  rcl,etamid,etaint,spdmax,deltim,nvcn,xvcn,
     &  dtdf,dtdl,
     &  td3_h1,td3_h2,rgt_h,ud3_h1,ud3_h2,vd3_h1,vd3_h2,
     &  sd_m,z_m,
     &  unlm,vnlm,tnlm,unla,vnla,tnla,
     &  dpdt_a,dpdt_d3,
     &  grad_gzlam,grad_gzphi,gz_grid,go_forward,batah,gg_tracers,tb_k,
     &  kdt,r_rt,omega_v)

      use machine        , only : kind_evod
!mjr  use machine        , only : kind_phys
      use resol_def      , only : latg2,levm1,levp1,levs,ntrac
      use coordinate_def , only : ak5,bk5,ck,dbk,sv_ecm,t_ecm,y_ecm
      use namelist_def   , only : ref_temp
      use layout1        , only : me
      use physcons       , only : rerth => con_rerth,
     &                               rd => con_rd,
     &                               cp => con_cp,
     &                            omega => con_omega,
     &                             cvap => con_cvap,
     &                             grav => con_g   
      implicit none
      integer lon_dim_a,lon_dim_h,lons_lat
      integer nvcn,kdt
      real coriol,rcl,rk,sinra,deltim,xvcn,sinlat,speed2,omega_v
      real  ugsave,vgsave, tref , rdtref, r_rt 
      real etamid(levs),etaint(levp1)
!     real etamid(levs),etaint(levp1),wgtw(7)
      real grad_gzlam(lon_dim_a),grad_gzphi(lon_dim_a)
      real    gz_grid(lon_dim_a)
      real spdmax(levs),
     1    dg(lon_dim_a,levs),   tg(lon_dim_a,levs),
     1  dtdf(lon_dim_a,levs), dtdl(lon_dim_a,levs),
     1  sd_m(lon_dim_a     ),  z_m(lon_dim_a,levs),
     1dpdt_a(lon_dim_a     ),
     1  unlm(lon_dim_a,levs),
     1  vnlm(lon_dim_a,levs), tnlm(lon_dim_a,levs),
     1  unla(lon_dim_a,levs),
     1  vnla(lon_dim_a,levs), tnla(lon_dim_a,levs),
     1   rgs(lon_dim_a,levs,ntrac),
     1  dqdphi(lon_dim_a), dqdlam(lon_dim_a), qg(lon_dim_a)
      real  ww(lon_dim_h,levs),   wfilt(lon_dim_h,levp1),
     1      zz(lon_dim_h,levp1),  
     1      ug(lon_dim_h,levs),        vg(lon_dim_h,levs),
     1    ug_m(lon_dim_h,levs),      vg_m(lon_dim_h,levs),
     1   ud3_h1(lon_dim_h,levs),     vd3_h1(lon_dim_h,levs),
     1   ud3_h2(lon_dim_h,levs),     vd3_h2(lon_dim_h,levs),
     1     u_n(lon_dim_h,levs),       u_l(lon_dim_h,levs),
     1     v_n(lon_dim_h,levs),       v_l(lon_dim_h,levs),
     1     t_n(lon_dim_h,levs),       t_l(lon_dim_h,levs),
     1                            dpdt_d3(lon_dim_h,levs),
     1   td3_h1(lon_dim_h,levs),td3_h2(lon_dim_h,levs),
     1   tb(lon_dim_a,levs),
     1   rgt_h(lon_dim_h,levs,ntrac)
      real pk5(lons_lat,levp1), dpk(lons_lat,levs),sv(levs)
      real
     &     dot(lons_lat,levp1),   
     &    advz(lons_lat,levs),      cg(lons_lat,levs),
     &      cb(lons_lat,levs),      db(lons_lat,levs),
     &    zlam(lons_lat,levs),    zphi(lons_lat,levs),
     &   worka(lons_lat,levs),   workb(lons_lat,levs),
     &   workc(lons_lat,levs),
     &    phiu(lons_lat,levs),    phiv(lons_lat,levs),
     &    uprs(lons_lat,levs),    vprs(lons_lat,levs),
     &    cofa(lons_lat,levs),    cofb(lons_lat,levs),
     &    alfa(lons_lat,levs),    rlnp(lons_lat,levs),
     &    px1u(lons_lat,levs),    px1v(lons_lat,levs),
     &    px2u(lons_lat,levs),    px2v(lons_lat,levs),
     &     px2(lons_lat,levs),
     &    px3u(lons_lat,levs),    px3v(lons_lat,levs),
     &    px4u(lons_lat,levs),    px4v(lons_lat,levs),
     &    px5u(lons_lat,levs),    px5v(lons_lat,levs),
     &    uphi(lons_lat,levs),    vphi(lons_lat,levs),
     &    sadx(lons_lat,levs),    sadk(lons_lat,levs),
     &    expq(lons_lat),         sadt(lons_lat,levs),
     &    dqdt(lons_lat),
     1    sd_n(lons_lat)     ,   z_n(lons_lat  ,levs),
     & z_n_sum(lons_lat),
     &    rdel(lons_lat,levs),   rdel2(lons_lat,levs),
     &  sumdel(lons_lat),          del(lons_lat,levs),
     &    rk1,rkr,                 dif(lons_lat)
      real clog2
      real rmin,rmax,delta,delta1,batah,batah2
      logical gg_tracers,test_mode,go_forward
      real rcoslat,coslat,tbar
      real tb_k(levs)           
      integer i,j,k,n,ifirst,lat,k_row,k_col,kk,kk1
      real, parameter :: cons0=0.d0, cons0p5=0.5d0, cons1=1.d0,
     &                   cons2=2.d0
!     save clog2,ifirst,delta,delta1
!     data ifirst /1/
!!      data wgtw /0.015625,-0.09375,0.234375,0.6875,
!!     .           0.234375,-0.09375,0.015625/

!!    data wgtw /-.03125,0.,.28125,.5,.28125,0.,-.03125/

!Moor   data wgtw /0.015625, 0.09375,0.234375,0.3125,
!    .             0.234375, 0.09375,0.015625/

      rk= rd /cp
!--------------------------------------------------
      call top2bot_array(lon_dim_a,lons_lat,dg)
      call top2bot_array(lon_dim_a,lons_lat,tg)
      call top2bot_array(lon_dim_h,lons_lat,ug)
      call top2bot_array(lon_dim_h,lons_lat,vg)
      call top2bot_array(lon_dim_a,lons_lat,dtdf)
      call top2bot_array(lon_dim_a,lons_lat,dtdl)

      call top2bot_ntrac(lon_dim_a,lons_lat,rgs)
!
      if(gg_tracers .and. kdt == 1)then          ! do only once at start
        do k=1,levs
          do j=1,lons_lat
            rgt_h(j,k,1) = max(0.0, rgs(j,k,1))
          enddo
        enddo
        if (ntrac > 1) then
          do n=2,ntrac
            do k=1,levs
              do j=1,lons_lat
                rgt_h(j,k,n) = rgs(j,k,n)
              enddo
            enddo
          enddo
        endif
      endif                                      ! tracers at kdt.eq.1

      if (.not. gg_tracers) then
        do n=1,ntrac
          do k=1,levs
            do j=1,lons_lat
              rgt_h(j,k,n) = rgs(j,k,n)
            enddo
          enddo
        enddo

!     if (me == 1 .and. lat == 48) then
!     if (lat == 48) then
!     print *,' OZONE in GFIDI=',rgt_h(1,1:5,2),' at lat=',lat,' me=',me
!     print *,' rgs in GFIDI=',rgs(1,1:5,2),' at lat=',lat,' me=',me
!     endif

      endif

      sinra  = sqrt(cons1-cons1/rcl)
      coriol = cons2*omega*sinra
      if (lat > latg2) coriol = -coriol
!     coriol=0.

      rcoslat = sqrt(rcl)      ! =1/cos(lat)
      coslat  = 1.0 / rcoslat  ! cos(lat)
      clog2   = log(cons2)     ! constant
      delta   = cvap/cp  ! check these cpv cpd (at const p for vapor and dry
      delta1  = delta-cons1
      rk1     = rk + 1.e0
      rkr     = 1.0/rk


!     ifirst=0

      do j=1,lons_lat
        expq(j)=exp(qg(j))
      enddo
      do k=1,levp1
        do j=1,lons_lat
          pk5(j,k) = ak5(k) + bk5(k)*expq(j)
        enddo
      enddo
      do k=1,levs
        do j=1,lons_lat
          dpk(j,k)   = pk5(j,k+1) - pk5(j,k)
          rdel(j,k)  = cons1/dpk(j,k)
          rdel2(j,k) = 0.5*rdel(j,k)*gz_grid(j)
        enddo
      enddo
      k = 1
      do j=1,lons_lat
        alfa(j,k)=clog2
        rlnp(j,1)=cons0
      enddo
      do  k=2,levs
        do j=1,lons_lat
          rlnp(j,k) = log( pk5(j,k+1)/pk5(j,k) )
          alfa(j,k) = cons1-( pk5(j,k)/dpk(j,k) )*rlnp(j,k)
        enddo
      enddo
      spdmax = 0.

      do  k=1,levs
        do j=1,lons_lat
          ug(j,k)   = ug(j,k)*rcoslat ! ug at input is scaled U=u*cos(theta)
          vg(j,k)   = vg(j,k)*rcoslat ! vg at input is scaled V=v*cos(theta)

          dtdl(j,k) = dtdl(j,k)*coslat !dtdl in=(dt/dlamda  )/cos**2(theta)
          dtdf(j,k) = dtdf(j,k)*coslat !dtdf in=(dt/d(theta))/cos**2(theta)

          speed2    = ug(j,k)*ug(j,k)+vg(j,k)*vg(j,k)

          if (speed2  > spdmax(k))  spdmax(k) = speed2 
        enddo
      enddo

      do j=1,lons_lat
        dqdlam(j) = dqdlam(j)*rcoslat                 
        dqdphi(j) = dqdphi(j)*rcoslat                 
      enddo
         
      do k=1,levs
        do j=1,lons_lat
          cg(j,k) = ug(j,k)*dqdlam(j) + vg(j,k)*dqdphi(j)
        enddo
      enddo

      k = 1
      do j=1,lons_lat
        db(j,k) = dg(j,k)*dpk(j,k)
        cb(j,k) = cg(j,k)*dbk(k)
      enddo

      do k=1,levm1
        do j=1,lons_lat
          db(j,k+1) = db(j,k) + dg(j,k+1)*dpk(j,k+1)
          cb(j,k+1) = cb(j,k) + cg(j,k+1)*dbk(k+1)
        enddo
      enddo

!sela dqdt=partial(ln(p_s)/partial(t)
      do j=1,lons_lat
        dqdt(j)      = -db(j,levs)/expq(j)-cb(j,levs) !old eulerian dqdt
        dot(j,    1) = cons0
        dot(j,levp1) = cons0
      enddo
!sela dot=( eta_dot*d(p)/d(eta) )(k+1/2) as in Eulerian 
      do k=1,levs-1
        do j=1,lons_lat
          dot(j,k+1) = -expq(j)*(bk5(k+1)*dqdt(j)+cb(j,k)) -db(j,k)
        enddo
      enddo
!sela ww=eta_dot
      do k=2,levs
        do j=1,lons_lat
          zz(j,k) = dot(j,k)
     &     *( etaint(k+1)-etaint(k-1) )/( pk5(j,k+1)-pk5(j,k-1) )
        enddo
      enddo
      do j=1,lons_lat
        zz(j,    1) = 0.0
        zz(j,levp1) = 0.0
      enddo

!     do j=1,lons_lat
!$$$       wfilt(j,    1)=0.0
!$$$       wfilt(j,levp1)=0.0
!     enddo
!$$$!    filter ertical velocities
!$$$!     k=2
!$$$        do i=1,lons_lat
!$$$         wfilt(i,2)=-zz(i,3)*wgtw(1)
!$$$     .              -zz(i,2)*wgtw(2)
!$$$! zz(i,1)=0         +zz(i,1)*wgtw(3)
!$$$     .              +zz(i,2)*wgtw(4)
!$$$     .              +zz(i,3)*wgtw(5)
!$$$     .              +zz(i,4)*wgtw(6)
!$$$     .              +zz(i,5)*wgtw(7)
!$$$        enddo
!$$$!     k=3
!$$$        do i=1,lons_lat
!$$$         wfilt(i,3)=-zz(i,2)*wgtw(1)
!$$$! zz(i,1)=0         -zz(i,1)*wgtw(2)
!$$$     .              +zz(i,2)*wgtw(3)
!$$$     .              +zz(i,3)*wgtw(4)
!$$$     .              +zz(i,4)*wgtw(5)
!$$$     .              +zz(i,5)*wgtw(6)
!$$$     .              +zz(i,6)*wgtw(7)
!$$$        enddo
!$$$
!$$$       do k=4,levs-2
!$$$        do i=1,lons_lat
!$$$         wfilt(i,k)=zz(i,k-3)*wgtw(1)
!$$$     .             +zz(i,k-2)*wgtw(2)
!$$$     .             +zz(i,k-1)*wgtw(3)
!$$$     .             +zz(i,k  )*wgtw(4)
!$$$     .             +zz(i,k+1)*wgtw(5)
!$$$     .             +zz(i,k+2)*wgtw(6)
!$$$     .             +zz(i,k+3)*wgtw(7)
!$$$        enddo
!$$$       enddo
!$$$!     k=levs-1
!$$$        do i=1,lons_lat
!$$$         wfilt(i,levs-1)= zz(i,levs-4)*wgtw(1)
!$$$     .                   +zz(i,levs-3)*wgtw(2)
!$$$     .                   +zz(i,levs-2)*wgtw(3)
!$$$     .                   +zz(i,levs-1)*wgtw(4)
!$$$     .                   +zz(i,levs  )*wgtw(5)
!$$$! zz(levs+1)=0.          +zz(i,levs+1)*wgtw(6)
!$$$     .                   -zz(i,levs  )*wgtw(7)
!$$$        enddo
!$$$!     k=levs   
!$$$        do i=1,lons_lat
!$$$         wfilt(i,levs  )= zz(i,levs-3)*wgtw(1)
!$$$     .                   +zz(i,levs-2)*wgtw(2)
!$$$     .                   +zz(i,levs-1)*wgtw(3)
!$$$     .                   +zz(i,levs  )*wgtw(4)
!$$$! zz(levs+1)=0           +zz(i,levs+1)*wgtw(5)
!$$$     .                   -zz(i,levs  )*wgtw(6)
!$$$     .                   -zz(i,levs-1)*wgtw(7)
!$$$        enddo
!sela do k=1,levs
!sela    do j=1,lons_lat
!sela       ww(j,k)=0.5*(wfilt(j,k)+wfilt(j,k+1) )
!sela    enddo
!$$$! fin   filter ertical velocities

      tref   = ref_temp
      rdtref = rd*tref
!
      do k=1,levs
        do j=1,lons_lat
          ww(j,k) = 0.5*(zz(j,k)+zz(j,k+1) )

          tb(j,k) = tb_k(k)*gz_grid(j)                          
!
          sadx(j,k) = ug(j,k)*grad_gzlam(j) + vg(j,k)*grad_gzphi(j)     
          advz(j,k) = sadx(j,k)*r_rt  
        enddo
      enddo
!
!sela check sign of rdel2
!
      do j=1,lons_lat
        sadk(j,1) = rdel2(j,1) * dot(j,2) * (tb_k(2)-tb_k(1))
      enddo
      k = levs
      do j=1,lons_lat
        sadk(j,k) = rdel2(j,k) * dot(j,k) * (tb_k(k)-tb_k(k-1))
      enddo

      do k=2,levm1
        do j=1,lons_lat
          sadk(j,k) = rdel2(j,k) * (dot(j,k+1) * (tb_k(k+1)-tb_k(k))
     &                           +  dot(j,k)   * (tb_k(k)-tb_k(k-1)))
        enddo
      enddo
      do k=1,levs
        do j=1,lons_lat
          sadt(j,k) = tb_k(k)*sadx(j,k) + sadk(j,k)
        enddo
      enddo
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!!!for testing input fields=output or any locally computed fields=output
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      test_mode = .false.
      if (test_mode) then !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        print*,' enter test mode in gfidi_tracers with lin temp.'
        do k=1,levs
          do j=1,lons_lat
           ud3_h1(j,k)=ug(j,k)*rcl !watch scaling for test
           vd3_h1(j,k)=vg(j,k)*rcl !watch scaling for test
!
            dpdt_d3(j,k) = qg(j)
            tg(j,k)      =  sadt(j,k) !for testing with ncar graphics thru td3_h
          td3_h1(j,k)=tg(j,k)
          enddo
        enddo
        do j=1,lons_lat
          dpdt_a(j) = qg(j)
        enddo
        return
      endif
!      if (test_mode) return !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      batah2 = batah*0.5
      do k=1,levs
        sv(k) = batah2*sv_ecm(k)
      enddo

      z_n_sum = 0.
      do k=1,levs
        do j=1,lons_lat
          z_n(j,k)   = dbk(k)*( dqdt(j) + cg(j,k) + advz(j,k) )
          z_n_sum(j) = z_n_sum(j) + z_n(j,k)
        enddo
      enddo
      do j=1,lons_lat
        sd_n(j) = sv(1)*dg(j,1)
      enddo
      do k=2,levs
        do j=1,lons_lat
          sd_n(j) = sd_n(j) + sv(k)*dg(j,k)
        enddo
      enddo

      if (go_forward) then
        do j=1,lons_lat
          dpdt_a(j) = deltim*( 0.5*z_n_sum(j)+sd_n(j) ) ! at arrival points
          sd_m(j)   = sd_n(j)
        enddo
        do k=1,levs
          do j=1,lons_lat
            dpdt_d3(j,k) = deltim*( 0.5*z_n(j,k) )
     &                   + ( qg(j) + gz_grid(j)*r_rt )*dbk(k)
            z_m(j,k)     = z_n(j,k)
          enddo
        enddo

      else  ! go_forward

        do k=1,levs
          do j=1,lons_lat
            dpdt_d3(j,k) = deltim*( z_n(j,k)-0.5*z_m(j,k) )
     &                   + ( qg(j) + gz_grid(j)*r_rt )*dbk(k)
     &                   + dbk(k)*deltim*( sd_n(j) - sd_m(j) )
            z_m(j,k)     = z_n(j,k)
          enddo
        enddo
        do j=1,lons_lat
          dpdt_a(j) = deltim*( 0.5*z_n_sum(j)+sd_n(j) ) ! at arrival points
          sd_m(j)   = sd_n(j)
        enddo
      endif !go_forward
!
      do j=1,lons_lat
        cofb(j,1)    = -rdel(j,1) * (alfa(j,1)*dbk(1))
        px2(j,levs)  = cons0
        px3u(j,levs) = cons0
        px3v(j,levs) = cons0
      enddo
!
      do k=2,levs
        kk  = levp1 - k + 1
        kk1 = kk + 1
        do j=1,lons_lat
          cofb(j,k)    = -rdel(j,k)*(bk5(k)*rlnp(j,k)+alfa(j,k)*dbk(k))

          px2(j,kk-1)  = px2(j,kk)
     &       - rd*(bk5(kk1)/pk5(j,kk1)-bk5(kk)/pk5(j,kk)) * tg(j,kk)
!
          px3u(j,kk-1) = px3u(j,kk) - rd*rlnp(j,kk)*dtdl(j,kk)
          px3v(j,kk-1) = px3v(j,kk) - rd*rlnp(j,kk)*dtdf(j,kk)
        enddo
      enddo
      do k=1,levs
        do j=1,lons_lat
          u_n(j,k) =  vg(j,k)*coriol ! ug,vg are not scaled in this version
          v_n(j,k) = -ug(j,k)*coriol ! ug,vg are not scaled in this version
!sela     u_n(j,k)=0.
!sela     v_n(j,k)=0.
!
          uprs(j,k) = cofb(j,k)*rd*tg(j,k)*expq(j)*dqdlam(j)
          vprs(j,k) = cofb(j,k)*rd*tg(j,k)*expq(j)*dqdphi(j)
!
          cofa(j,k) = -rdel(j,k)*(
     &                 bk5(k+1)*pk5(j,k)/pk5(j,k+1) - bk5(k)
     &              +  rlnp(j,k) * (bk5(k)-pk5(j,k)*dbk(k)*rdel(j,k)) )
!
          px1u(j,k) = -grad_gzlam(j)  ! grad_gz is expected to be true grad
          px1v(j,k) = -grad_gzphi(j)  ! grad_gz is expected to be true grad
!
          px2u(j,k) = px2(j,k)*expq(j)*dqdlam(j)
          px2v(j,k) = px2(j,k)*expq(j)*dqdphi(j)
        enddo
      enddo

!     do k=1,levs
!     do j=1,lons_lat
!      px3u(j,k)=px3u(j,k)/rcl
!      px3v(j,k)=px3v(j,k)/rcl
!     enddo
!     enddo

      do k=1,levs
        do j=1,lons_lat
          px4u(j,k) = -rd*alfa(j,k)*dtdl(j,k) ! in scaled version /rcl
          px4v(j,k) = -rd*alfa(j,k)*dtdf(j,k) ! in scaled version /rcl
!
          px5u(j,k) = -cofa(j,k)*rd*tg(j,k)*expq(j)*dqdlam(j)
          px5v(j,k) = -cofa(j,k)*rd*tg(j,k)*expq(j)*dqdphi(j)
!
          uphi(j,k) = px1u(j,k)+px2u(j,k)+px3u(j,k)+px4u(j,k)+px5u(j,k)
          vphi(j,k) = px1v(j,k)+px2v(j,k)+px3v(j,k)+px4v(j,k)+px5v(j,k)
!
          u_n(j,k) = u_n(j,k) + uphi(j,k) + uprs(j,k)
          v_n(j,k) = v_n(j,k) + vphi(j,k) + vprs(j,k)
!
          worka(j,k) = rk*tg(j,k)/( cons1+delta1*rgt_h(j,k,1) )
     &                           * rdel(j,k)
        enddo
      enddo

      do j=1,lons_lat
        workb(j,1)=
     &  alfa(j,1)*( dg(j,1)*dpk(j,1)+expq(j)*cb(j,1)*dbk(1) )
!
        workc(j,1) = expq(j)*cg(j,1)*dbk(1)
      enddo

      do k=2,levs
        do j=1,lons_lat
          workb(j,k) = rlnp(j,k) * (db(j,k-1)+expq(j)*cb(j,k-1))
     &    + alfa(j,k) * (dg(j,k)*dpk(j,k)+expq(j)*cg(j,k)*dbk(k))
!
          workc(j,k) = expq(j) * cg(j,k)
     &               * ( dbk(k) + ck(k)*rlnp(j,k)*rdel(j,k) )
        enddo
      enddo

      do k=1,levs
        do j=1,lons_lat
          t_n(j,k) = worka(j,k) * ( -workb(j,k) + workc(j,k))
        enddo
      enddo

!$$$  do k=1,levs
!$$$  do j=1,lons_lat
!$$$   u_l(j,k)= 0.
!$$$   v_l(j,k)= 0.
!$$$   t_l(j,k)= 0.
!$$$  enddo
!$$$  enddo

!$$$  do k_row=1,levs
!$$$   do k_col=1,levs
!$$$    do j=1,lons_lat
!$$$     u_l(j,k_row)=u_l(j,k_row)+y_ecm(k_row,k_col)*dtdl(j,k_col)
!$$$     v_l(j,k_row)=v_l(j,k_row)+y_ecm(k_row,k_col)*dtdf(j,k_col)
!$$$     t_l(j,k_row)=t_l(j,k_row)+t_ecm(k_row,k_col) * dg(j,k_col)
!$$$    enddo
!$$$   enddo
!$$$  enddo

      if ( kind_evod .eq. 8 ) then !------------------------------------
         call dgemm ('n', 't',
     &               lons_lat, levs, levs, cons1,
     &               dtdl(1,1), lon_dim_a,
     &               y_ecm(1,1), levs, cons0,
     &               u_l(1,1), lon_dim_h)
         call dgemm ('n', 't',
     &               lons_lat, levs, levs, cons1,
     &               dtdf(1,1), lon_dim_a,
     &               y_ecm(1,1), levs, cons0,
     &               v_l(1,1), lon_dim_h)
         call dgemm ('n', 't',
     &               lons_lat, levs, levs, cons1,
     &               dg(1,1), lon_dim_a,
     &               t_ecm(1,1), levs, cons0,
     &               t_l(1,1), lon_dim_h)
      else !------------------------------------------------------------
         call sgemm ('n', 't',
     &               lons_lat, levs, levs, cons1,
     &               dtdl(1,1), lon_dim_a,
     &               y_ecm(1,1), levs, cons0,
     &               u_l(1,1), lon_dim_h)
         call sgemm ('n', 't',
     &               lons_lat, levs, levs, cons1,
     &               dtdf(1,1), lon_dim_a,
     &               y_ecm(1,1), levs, cons0,
     &               v_l(1,1), lon_dim_h)
         call sgemm ('n', 't',
     &               lons_lat, levs, levs, cons1,
     &               dg(1,1), lon_dim_a,
     &               t_ecm(1,1), levs, cons0,
     &               t_l(1,1), lon_dim_h)
      endif !-----------------------------------------------------------

      do k=1,levs
        do j=1,lons_lat
          u_l(j,k) = rdtref*dqdlam(j) + u_l(j,k) ! in scaled version /rcl
          v_l(j,k) = rdtref*dqdphi(j) + v_l(j,k) ! in scaled version /rcl
        enddo
      enddo


      if (go_forward) then
        do k=1,levs
          do j=1,lons_lat
            unla(j,k) = u_n(j,k) + batah*u_l(j,k)
            vnla(j,k) = v_n(j,k) + batah*v_l(j,k)
            tnla(j,k) = t_n(j,k) + batah*t_l(j,k) - sadt(j,k)

!sela   ud3_h(j,k)=ug(j,k)+omega_v+
        ud3_h1(j,k) = ug(j,k)
        ud3_h2(j,k) = deltim*(0.5*unla(j,k)-batah2*u_l(j,k))       

        vd3_h1(j,k) = vg(j,k)
        vd3_h2(j,k) = deltim*(0.5*vnla(j,k)-batah2*v_l(j,k))

        td3_h1(j,k) = tg(j,k)
        td3_h2(j,k) = deltim*(0.5*tnla(j,k)-batah2*t_l(j,k)) -tb(j,k)

            unlm(j,k)  = unla(j,k)
            vnlm(j,k)  = vnla(j,k)
            tnlm(j,k)  = tnla(j,k)

            ug_m(j,k)  = ug(j,k) !for midpoint value in dep. point calc
            vg_m(j,k)  = vg(j,k) !for midpoint value in dep. point calc
          enddo
        enddo

      else  ! go_forward

        do k=1,levs
          do j=1,lons_lat
            unla(j,k) = u_n(j,k) + batah*u_l(j,k)
            vnla(j,k) = v_n(j,k) + batah*v_l(j,k)
            tnla(j,k) = t_n(j,k) + batah*t_l(j,k)-sadt(j,k)
!sela   ud3_h(j,k)=ug(j,k)+omega_v+
        ud3_h1(j,k) = ug(j,k)
        ud3_h2(j,k) = deltim*(unla(j,k)-0.5*unlm(j,k)-batah2*u_l(j,k))
        vd3_h1(j,k) = vg(j,k)
        vd3_h2(j,k) = deltim*(vnla(j,k)-0.5*vnlm(j,k)-batah2*v_l(j,k))
        td3_h1(j,k) = tg(j,k) 
        td3_h2(j,k) = deltim*
     .               (tnla(j,k)-0.5*tnlm(j,k)-batah2*t_l(j,k)) - tb(j,k)
            unlm(j,k)  = unla(j,k)
            vnlm(j,k)  = vnla(j,k)
            tnlm(j,k)  = tnla(j,k)
          enddo
        enddo

      endif ! go_forward

      do k=1,levs
        do j=1,lons_lat

!!!!    ud3_h(j,k)=ud3_h(j,k)
!!!!    vd3_h(j,k)=vd3_h(j,k)

          ugsave    = ug(j,k)
          ug(j,k)   = 1.5*ug(j,k) - 0.5*ug_m(j,k)
          ug_m(j,k) = ugsave

          vgsave    = vg(j,k)
          vg(j,k)   = 1.5*vg(j,k) - 0.5*vg_m(j,k)
          vg_m(j,k) = vgsave
        enddo
      enddo

      call top2bot_array(1,1,spdmax)


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!!!for testing input fields=output or any locally computed fields=output
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      test_mode =.false.
      if (test_mode) then !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
       print*,' enter test mode in gfidi_tracers with lin temp.'

c$$$         do j=1,lons_lat
c$$$           dpdt_a(j)=qg(j)
c$$$         enddo
        do k=1,levs
          do j=1,lons_lat

            unlm(j,k) = unla(j,k) * coslat 
            vnlm(j,k) = vnla(j,k) * coslat
            tnlm(j,k) = tnla(j,k) ! for testing with ncar graphics thru td3_h
           enddo
        enddo
      endif
       if(test_mode) return !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

      return
      end

      subroutine top2bot_array(lon_dim,lons_lat,array)
!mjr  use machine   , only : kind_phys
      use resol_def , only : levp1,levs
!
      implicit none
      integer lon_dim,lons_lat,i,k,kk
      real array(lon_dim,levs),tempo(lons_lat)
!
      do k=1,levs/2
        kk = levp1-k
        do i=1,lons_lat
          tempo(i)    = array(i,kk)
          array(i,kk) = array(i,k)
          array(i,k)  = tempo(i)
        enddo
      enddo
      return
      end
      subroutine top2bot_ntrac(lon_dim,lons_lat,array)
!mjr  use machine   , only : kind_phys
      use resol_def , only : levp1,levs,ntrac
!
      implicit none
      integer lon_dim,lons_lat,i,k,n,kk
      real array(lon_dim,levs,ntrac),tempo(lons_lat)
      do n=1,ntrac
        do k=1,levs/2
          kk = levp1-k
          do i=1,lons_lat
            tempo(i)      = array(i,kk,n)
            array(i,kk,n) = array(i,k,n)
            array(i,k,n)  = tempo(i)
          enddo
        enddo
      enddo
      return
      end
