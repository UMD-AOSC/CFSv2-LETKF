      subroutine sltini_n_uvw(n_ggp,lon_dim_h,uu,vv,ww,yhalo,lons_lat)
      use resol_def          , only : levs
!mjr  use layout1
      use layout_lag         , only : lats_dim_h
!mjr  include 'pmgrid.com.h'
      use      pmgrid        , only : plev
      implicit none
      integer lons_lat
      real
     $        gg1_mal(lons_lat+3),     ! sin(lamda_at_last_gg_lat)
     $     singg1_mal(lons_lat),      ! sin(lamda_at_last_gg_lat)
     $     cosgg1_mal(lons_lat)       ! cos(lamda_at_last_gg_lat)
      real   ww(lon_dim_h , levs  ,lats_dim_h)
      real   uu(lon_dim_h , levs  ,lats_dim_h)
      real   vv(lon_dim_h , levs  ,lats_dim_h)
      integer i,j,jj,k             ! index
      real zave,zave1,zave2,zave3  ! accumulator for zonal averaging
      real zavec1,zavec2,          ! accumulator for wavenumber 1
     $     zaves1,zaves2           ! accumulator for wavenumber 1
      real lam_0,pi,degrad,dgg1_mal
!mjr  integer ig,lonh,n_ggp,lon_dim_h,yhalo
      integer ig,     n_ggp,lon_dim_h,yhalo
      integer ifirst
      data    ifirst /0/
      save    ifirst
      ifirst =ifirst + 1
!sela print*,'from sltini_nh_uvw ifirst=',ifirst,' n_ggp=',n_ggp,
!sela&       ' lon_dim_h=',lon_dim_h,' lota=',lota
!mjr  lonh=lon_dim_h
      lam_0 = 0.0
      pi = 4.*atan(1.)
      degrad = 180./pi
      dgg1_mal = 2.*pi/float(lons_lat)
      do i = 1,lons_lat+3
       gg1_mal(i) = float(i-2)*dgg1_mal + lam_0
      end do
      ig = 0
      do i = 2,lons_lat+1
       ig = ig + 1
       singg1_mal(ig) = sin( gg1_mal(i) )
       cosgg1_mal(ig) = cos( gg1_mal(i) )
      end do
      do k = 1,plev
       zavec1 = 0.0
       zavec2 = 0.0
       zaves1 = 0.0
       zaves2 = 0.0
       ig = 0
        do i = 2,lons_lat+1
         ig    = ig + 1
          zavec1 = zavec1 + vv(i,k,n_ggp)*cosgg1_mal(ig)
          zaves1 = zaves1 + vv(i,k,n_ggp)*singg1_mal(ig)
          zavec2 = zavec2 + uu(i,k,n_ggp)*cosgg1_mal(ig)
          zaves2 = zaves2 + uu(i,k,n_ggp)*singg1_mal(ig)
         end do
        zavec1 = zavec1/float(lons_lat/2)
        zaves1 = zaves1/float(lons_lat/2)
        zavec2 = zavec2/float(lons_lat/2)
        zaves2 = zaves2/float(lons_lat/2)
        ig    = 0
         do i = 2,lons_lat+1
          ig           = ig + 1
c$$$          vv(i,k,n_ggp+1)=0.5*(zavec1+zaves2)*cosgg1_mal(ig)
c$$$     &                            +0.5*(zaves1-zavec2)*singg1_mal(ig)
c$$$          uu(i,k,n_ggp+1)=0.5*(zavec2-zaves1)*cosgg1_mal(ig)
c$$$     &                            +0.5*(zaves2+zavec1)*singg1_mal(ig)

          vv(i,k,n_ggp+1)=zavec1*cosgg1_mal(ig)
     &                            +zaves1*singg1_mal(ig)
          uu(i,k,n_ggp+1)=zavec2*cosgg1_mal(ig)
     &                            +zaves2*singg1_mal(ig)

         end do
      end do
            do k = 1,plev
               do i = 2,1+(lons_lat/2)
                  vv(         i,k,n_ggp+2  ) =
     &           -vv((lons_lat/2)+i,k,n_ggp)

                  vv((lons_lat/2)+i,k,n_ggp+2  ) =
     &           -vv(         i,k,n_ggp)

                  uu(         i,k,n_ggp+2  ) =
     &           -uu((lons_lat/2)+i,k,n_ggp)

                  uu((lons_lat/2)+i,k,n_ggp+2  ) =
     &           -uu(         i,k,n_ggp)

               end do
            end do
      do k=1,plev
         zave = 0.0
         do i=2,lons_lat+1
            zave = zave + ww(i,k,n_ggp)
         end do

           zave = zave/float(lons_lat)

         do i=2,lons_lat+1
           ww(i,k,n_ggp+1) = zave
         end do

      end do !plev loop

            do k=1,plev
               do i=2,1+(lons_lat/2)
                  ww(         i,k,n_ggp+2  ) =
     &            ww((lons_lat/2)+i,k,n_ggp)

                  ww((lons_lat/2)+i,k,n_ggp+2  ) =
     &            ww(         i,k,n_ggp)
               end do
            end do

         do jj=1,2
                     j=n_ggp+jj
            do k=1,plev
               uu(     1,k,j) = uu(1+lons_lat,k,j)
               uu(lons_lat+2,k,j) = uu(2,k,j)
               uu(lons_lat+3,k,j) = uu(3,k,j)

               vv(     1,k,j) = vv(1+lons_lat,k,j)
               vv(lons_lat+2,k,j) = vv(2,k,j)
               vv(lons_lat+3,k,j) = vv(3,k,j)

               ww(     1,k,j) = ww(1+lons_lat,k,j)
               ww(lons_lat+2,k,j) = ww(2,k,j)
               ww(lons_lat+3,k,j) = ww(3,k,j)
            end do
          ww(     1,k,j) = ww(1+lons_lat,k,j)
          ww(lons_lat+2,k,j) = ww(2     ,k,j)
          ww(lons_lat+3,k,j) = ww(3     ,k,j)
         end do
      if ( ifirst .eq. 1 ) then
         do jj=1,2
                       j=n_ggp+jj
            do k=1,plev
               do i=1,lons_lat+3
               enddo
            enddo
               k=1
            do k=1,plev
               do i=1,lons_lat+3
               enddo
            enddo
         enddo
      endif
      return
      end
      subroutine sltini_nh_scalar(n_ggp,lon_dim_h,q3,yhalo,lons_lat)
      use resol_def          , only : levs
!mjr  use layout1
      use layout_lag         , only : lats_dim_h
!mjr  include 'pmgrid.com.h'
      use      pmgrid        , only : plev
      implicit none
      integer lons_lat
      real   q3(lon_dim_h , levs  ,lats_dim_h)
      integer i,j,jj,k             ! index
      real zave,zave1,zave2,zave3  ! accumulator for zonal averaging
      real zavec1,zavec2,          ! accumulator for wavenumber 1
     $     zaves1,zaves2           ! accumulator for wavenumber 1
!mjr  integer ig,lonh,n_ggp,lon_dim_h,yhalo
      integer ig,     n_ggp,lon_dim_h,yhalo
      integer ifirst
      data    ifirst /0/
      save    ifirst
      ifirst =ifirst + 1
!     print*,'hello from sltini_nh ifirst=',ifirst,' n_ggp=',n_ggp,
!    &       ' lon_dim_h=',lon_dim_h,' lota=',lota
!mjr  lonh=lon_dim_h
      do k=1,plev
         zave1 = 0.0
         do i=2,lons_lat+1
            zave1 = zave1 + q3(i,k ,n_ggp)
         end do
            zave1 = zave1/float(lons_lat)
         do i=2,lons_lat+1
            q3(i,k ,n_ggp+1) = zave1
         end do
      end do

            do k=1,plev
               do i=2,1+(lons_lat/2)
                  q3(i,k,n_ggp+2) = q3((lons_lat/2)+i,k ,n_ggp)
         q3((lons_lat/2)+i,k,n_ggp+2) = q3(i,k ,n_ggp)
               end do
            end do

         do jj=1,2
          j=n_ggp+jj
            do k=1,plev
               q3(     1,k ,j) = q3(1+lons_lat,k  ,j)
               q3(lons_lat+2,k ,j) = q3(2,     k  ,j)
               q3(lons_lat+3,k ,j) = q3(3,     k  ,j)
            end do
         end do

      if ( ifirst .eq. 1 ) then
         do jj=1,2
          j=n_ggp+jj
            do k=1,plev
               do i=1,lons_lat+3
               enddo
            enddo
         enddo
      endif
      return
      end
      subroutine sltini_nh_scalar1(n_ggp,lon_dim_h,q3,yhalo,lons_lat)
!mjr  use resol_def
!mjr  use layout1
      use layout_lag , only : lats_dim_h
!mjr  include 'pmgrid.com.h'
!mjr  use      pmgrid
      implicit none
      integer lons_lat
      real   q3(lon_dim_h ,        lats_dim_h)
      integer i,j,jj               ! index
      real zave,zave1,zave2,zave3  ! accumulator for zonal averaging
      real zavec1,zavec2,          ! accumulator for wavenumber 1
     $     zaves1,zaves2           ! accumulator for wavenumber 1
!mjr  integer ig,lonh,n_ggp,lon_dim_h,yhalo
      integer ig,     n_ggp,lon_dim_h,yhalo
      integer ifirst
      data    ifirst /0/
      save    ifirst
      ifirst =ifirst + 1
!     print*,'hello from sltini_nh ifirst=',ifirst,' n_ggp=',n_ggp,
!    &       ' lon_dim_h=',lon_dim_h,' lota=',lota
!mjr  lonh=lon_dim_h
         zave1 = 0.0
         do i=2,lons_lat+1
            zave1 = zave1 + q3(i,n_ggp)
         end do
            zave1 = zave1/float(lons_lat)
         do i=2,lons_lat+1
            q3(i,n_ggp+1) = zave1
         end do

               do i=2,1+(lons_lat/2)
         q3(i,n_ggp+2) = q3((lons_lat/2)+i,n_ggp )
         q3((lons_lat/2)+i,n_ggp+2) = q3(i,n_ggp )
               end do

         do jj=1,2
          j=n_ggp+jj
               q3(     1,j) = q3(1+lons_lat,j)
               q3(lons_lat+2,j) = q3(2,     j)
               q3(lons_lat+3,j) = q3(3,     j)
         end do

      if ( ifirst .eq. 1 ) then
         do jj=1,2
          j=n_ggp+jj
               do i=1,lons_lat+3
               enddo
         enddo
      endif
      return
      end
      subroutine sltini_s_uvw(s_ggp,lon_dim_h,uu,vv,ww,yhalo,lons_lat)
      use resol_def          , only : levs
!mjr  use layout1
      use layout_lag         , only : lats_dim_h
!mjr  include 'pmgrid.com.h'
      use      pmgrid        , only : plev
      implicit none
      integer lons_lat
      real
     $        gg1_mal(lons_lat+3),     ! sin(lamda_at_last_gg_lat)
     $     singg1_mal(lons_lat),      ! sin(lamda_at_last_gg_lat)
     $     cosgg1_mal(lons_lat)       ! cos(lamda_at_last_gg_lat)
      real   ww(lon_dim_h , levs  ,lats_dim_h)
      real   uu(lon_dim_h , levs  ,lats_dim_h)
      real   vv(lon_dim_h , levs  ,lats_dim_h)
      integer i,j,jj,k             ! index
      real zave,zave1,zave2,zave3  ! accumulator for zonal averaging
      real zavec1,zavec2,          ! accumulator for wavenumber 1
     $     zaves1,zaves2           ! accumulator for wavenumber 1
      real lam_0,pi,degrad,dgg1_mal
!mjr  integer ig,lonh,s_ggp,lon_dim_h,yhalo
      integer ig,     s_ggp,lon_dim_h,yhalo
      integer ifirst
      data    ifirst /0/
      save    ifirst
      ifirst =ifirst + 1
!     print*,'hello from sltini_sh ifirst=',ifirst,' s_ggp=',s_ggp,
!    &       ' lon_dim_h=',lon_dim_h,' lota=',lota
!mjr  lonh=lon_dim_h
      lam_0 = 0.0
      pi = 4.*atan(1.)
      degrad = 180./pi
      dgg1_mal = 2.*pi/float(lons_lat)
      do i = 1,lons_lat+3
       gg1_mal(i) = float(i-2)*dgg1_mal + lam_0
      end do
      ig = 0
      do i = 2,lons_lat+1
       ig = ig + 1
       singg1_mal(ig) = sin( gg1_mal(i) )
       cosgg1_mal(ig) = cos( gg1_mal(i) )
      end do
      do k = 1,plev
         zavec1 = 0.0
         zaves1 = 0.0
         zavec2 = 0.0
         zaves2 = 0.0
         ig = 0
         do i = 2,lons_lat+1
            ig    = ig + 1
            zavec1 = zavec1 + vv(i,k,s_ggp  )*cosgg1_mal(ig)
            zaves1 = zaves1 + vv(i,k,s_ggp  )*singg1_mal(ig)
            zavec2 = zavec2 + uu(i,k,s_ggp  )*cosgg1_mal(ig)
            zaves2 = zaves2 + uu(i,k,s_ggp  )*singg1_mal(ig)
         end do
         zavec1 = zavec1/float(lons_lat/2)
         zaves1 = zaves1/float(lons_lat/2)
         zavec2 = zavec2/float(lons_lat/2)
         zaves2 = zaves2/float(lons_lat/2)
         ig    = 0
       do i = 2,lons_lat+1
        ig           = ig + 1
c$$$        vv(i,k,s_ggp-1)=0.5*(zavec1-zaves2)*cosgg1_mal(ig)
c$$$     &                          +0.5*(zaves1+zavec2)*singg1_mal(ig)
c$$$        uu(i,k,s_ggp-1)=0.5*(zavec2+zaves1)*cosgg1_mal(ig)
c$$$     &                          +0.5*(zaves2-zavec1)*singg1_mal(ig)

        vv(i,k,s_ggp-1)=zavec1*cosgg1_mal(ig)
     &                          +zaves1*singg1_mal(ig)
        uu(i,k,s_ggp-1)=zavec2*cosgg1_mal(ig)
     &                          +zaves2*singg1_mal(ig)

       end do
      end do
            do k = 1,plev
               do i = 2,1+(lons_lat/2)
                  vv(         i,k,s_ggp-2  ) =
     &           -vv((lons_lat/2)+i,k,s_ggp)     

                  vv((lons_lat/2)+i,k,s_ggp-2  ) =
     &           -vv(         i,k,s_ggp)     

                  uu(         i,k,s_ggp-2  ) =
     &           -uu((lons_lat/2)+i,k,s_ggp)     

                  uu((lons_lat/2)+i,k,s_ggp-2  ) =
     &           -uu(         i,k,s_ggp)     
               end do
            end do
      do k=1,plev
         zave = 0.0
         do i = 2,lons_lat+1
            zave = zave + ww(i,k,s_ggp  )
         end do
            zave = zave/float(lons_lat)
         do i=2,lons_lat+1
            ww(i,k,s_ggp-1) = zave
         end do
      end do !plev loop

            do k=1,plev
               do i=2,1+(lons_lat/2)
                  ww(         i,k,s_ggp-2  ) =
     &            ww((lons_lat/2)+i,k,s_ggp)

                  ww((lons_lat/2)+i,k,s_ggp-2  ) =
     &            ww(         i,k,s_ggp)       
               end do
            end do

         do jj=1,2
                     j=s_ggp-jj
            do k=1,plev
               uu(     1,k,j) = uu(1+lons_lat,k,j)
               uu(lons_lat+2,k,j) = uu(2,k,j)
               uu(lons_lat+3,k,j) = uu(3,k,j)
               vv(     1,k,j) = vv(1+lons_lat,k,j)
               vv(lons_lat+2,k,j) = vv(2,k,j)
               vv(lons_lat+3,k,j) = vv(3,k,j)
               ww(     1,k,j) = ww(1+lons_lat,k,j)
               ww(lons_lat+2,k,j) = ww(2,k,j)
               ww(lons_lat+3,k,j) = ww(3,k,j)
            end do
          ww(     1,k,j) = ww(1+lons_lat,k,j)
          ww(lons_lat+2,k,j) = ww(2,k,j)
          ww(lons_lat+3,k,j) = ww(3,k,j)
         end do
      if ( ifirst .eq. 1 ) then
         do jj=1,2
                       j=s_ggp-jj
            do k=1,plev
               do i=1,lons_lat+3
               enddo
            enddo
               k=1
            do k=1,plev
               do i=1,lons_lat+3
               enddo
            enddo
         enddo
      endif
      return
      end
      subroutine sltini_sh_scalar(s_ggp,lon_dim_h,q3,yhalo,lons_lat)
      use resol_def          , only : levs
!mjr  use layout1
      use layout_lag         , only : lats_dim_h
!mjr  include 'pmgrid.com.h'
      use      pmgrid        , only : plev
      implicit none
      integer lons_lat
      real   q3(lon_dim_h , levs  ,lats_dim_h)
      integer i,j,jj,k             ! index
      real zave,zave1,zave2,zave3  ! accumulator for zonal averaging
      real zavec1,zavec2,          ! accumulator for wavenumber 1
     $     zaves1,zaves2           ! accumulator for wavenumber 1
!mjr  integer ig,lonh,s_ggp,lon_dim_h,yhalo
      integer ig,     s_ggp,lon_dim_h,yhalo
      integer ifirst
      data    ifirst /0/
      save    ifirst
      ifirst =ifirst + 1
!     print*,'hello from sltini_sh ifirst=',ifirst,' s_ggp=',s_ggp,
!    &       ' lon_dim_h=',lon_dim_h,' lota=',lota
!mjr  lonh=lon_dim_h
      do k=1,plev
         zave1 = 0.0
         do i = 2,lons_lat+1
            zave1 = zave1 + q3(i,k         ,s_ggp  )
         end do
            zave1 = zave1/float(lons_lat)
         do i=2,lons_lat+1
            q3(i,k         ,s_ggp-1) = zave1
         end do
      end do
            do k=1,plev
               do i=2,1+(lons_lat/2)
                  q3(         i,k         ,s_ggp-2  ) =
     &            q3((lons_lat/2)+i,k         ,s_ggp)       

                  q3((lons_lat/2)+i,k         ,s_ggp-2  ) =
     &            q3(         i,k         ,s_ggp)
               end do
            end do

         do jj=1,2
            j=s_ggp-jj
            do k=1,plev
               q3(     1,k         ,j) = q3(1+lons_lat,k         ,j)
               q3(lons_lat+2,k         ,j) = q3(2,     k         ,j)
               q3(lons_lat+3,k         ,j) = q3(3,     k         ,j)
            end do
         end do
      if ( ifirst .eq. 1 ) then
         do jj=1,2
            j=s_ggp-jj
            do k=1,plev
               do i=1,lons_lat+3
               enddo
            enddo
         enddo
      endif
      return
      end
      subroutine sltini_sh_scalar1(s_ggp,lon_dim_h,q3,yhalo,lons_lat)
!mjr  use resol_def
!mjr  use layout1
      use layout_lag , only : lats_dim_h
!mjr  include 'pmgrid.com.h'
!mjr  use      pmgrid
      implicit none
      integer lons_lat
      real   q3(lon_dim_h ,        lats_dim_h)
      integer i,j,jj               ! index
      real zave,zave1,zave2,zave3  ! accumulator for zonal averaging
      real zavec1,zavec2,          ! accumulator for wavenumber 1
     $     zaves1,zaves2           ! accumulator for wavenumber 1
!mjr  integer ig,lonh,s_ggp,lon_dim_h,yhalo
      integer ig,     s_ggp,lon_dim_h,yhalo
      integer ifirst
      data    ifirst /0/
      save    ifirst
      ifirst =ifirst + 1
!     print*,'hello from sltini_sh ifirst=',ifirst,' s_ggp=',s_ggp,
!    &       ' lon_dim_h=',lon_dim_h,' lota=',lota
!mjr  lonh=lon_dim_h
         zave1 = 0.0
         do i = 2,lons_lat+1
            zave1 = zave1 + q3(i,s_ggp  )
         end do
            zave1 = zave1/float(lons_lat)
         do i=2,lons_lat+1
            q3(i,s_ggp-1) = zave1
         end do

               do i=2,1+(lons_lat/2)
                  q3(         i,s_ggp-2  ) =
     &            q3((lons_lat/2)+i,s_ggp)
                  q3((lons_lat/2)+i,s_ggp-2  ) =
     &            q3(         i,s_ggp)
               end do

         do jj=1,2
               j=s_ggp-jj
               q3(     1,j) = q3(1+lons_lat,j)
               q3(lons_lat+2,j) = q3(2,     j)
               q3(lons_lat+3,j) = q3(3,     j)
         end do
      if ( ifirst .eq. 1 ) then
         do jj=1,2
               j=s_ggp-jj
               do i=1,lons_lat+3
               enddo
         enddo
      endif
      return
      end
