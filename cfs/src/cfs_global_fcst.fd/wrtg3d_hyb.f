      subroutine wrtg3d_hyb(ioproc,nog3d,zhour,fhour,idate,colat1,
     &                 global_lats_r,lonsperlar,pl_coeff,SECSWR)
      use resol_def
      use layout1
      use sig_io
      use namelist_def
      use d3d_def
      implicit none
!
!2008/10/08 HO-Chun Huang   Output only 3D fileds Specific for GFS-GOCART
!                           Including DQ3DT, upd_mf, dwn_mf, det_mf, dkh, rnp
!
!!
      integer ioproc, pl_coeff
      integer idq3(5+pl_coeff), idq3a(7), idq3b(9), icc
      integer icmf(3), idkh, irnp
      data    idq3a/249, 243, 245, 217, 141, 142, 143/
      data    idq3b/249, 243, 245, 217, 141, 142, 143, 139, 140/
! rnp use large scale precipate = 62,  dkh use Velocity potential = 36
! dwn_mf=230, upd_mf=231, det_mf=236
!      data    icmf/230,231,236/, idkh /36/, irnp /62/, icc /213/
! table 133 dwn_mf=209, upd_mf=202, det_mf=219, rnp=134
! table 182 dkh   =182
      data    icmf/202,209,219/, idkh /182/, irnp /134/, icc /213/
!
      save idq3a, idq3b, icmf, idkh, irnp, icc
!
      integer   isglev, iavg, ifhour
      parameter(isglev=109, IAVG=3, IFHOUR=1)
      real(kind=kind_io4) wrkga(lonr*latr),wrkgb(lonr*latr)
!
      logical(1) lbm(lonr*latr)
      character g(200+lonr*latr*(16+1)/8)
      integer     idate(4), ids(255)
      integer nv,il1,il2,j,i,k,itop,ibot,k4,l,nog3d
      integer ilpds,iyr,imo,ida,ihr,ifhr,ithr,lg,ierr
      real (kind=kind_io8) colat1,cl1,secswr,zhour,fhour,rtime,rtimsw

      real (kind=kind_io8)   glolal(lonr,lats_node_r)
      real (kind=kind_io8)   buffo(lonr,lats_node_r)
      real (kind=kind_io4)   buff1(lonr,latr)
!!
      integer kmsk0(lonr*latr)
      real (kind=kind_io8)   glolal1(lonr*latr)
      integer               lats_nodes_r(nodes)
      integer              global_lats_r(latr)
      integer                 lonsperlar(latr)
      integer lan,lat,lons_lat,lon
!
      if (pl_coeff .eq. 2) then
        idq3(1:7) = idq3a(1:7)
      else
        idq3(1:9) = idq3b(1:9)
      endif
!
      ids=0
      g=' '
!
      kmsk0=0
!
      call idsdef(1,ids)
!                     For the Ozone diagnostics
      ids(139) = 13
      ids(140) = 15
      ids(141) = 13
      ids(142) = 13
      ids(143) = 13
!
      IDS(217) = 10   ! This is for moisture change due to large-scale
!
      ids(202) = 5    ! Cumulus updraft massflux
      ids(209) = 5    ! Cumulus downdraft massflux
      ids(219) = 5    ! Cumulus detrainment massflux
!
      IDS(idkh) = 7   ! vertical diffusivity
      IDS(irnp) = 7   ! rain production rate
!
      ilpds=28
      if(icen2.eq.2) ilpds=45
      iens(1)=1
      iens(2)=ienst
      iens(3)=iensi
      iens(4)=1
      iens(5)=255
      iyr=idate(4)
      imo=idate(2)
      ida=idate(3)
      ihr=idate(1)
      ifhr=nint(zhour)
      ithr=nint(fhour)
      if(fhour.gt.zhour) then
        rtime=1./(3600.*(fhour-zhour))
      else
        rtime=0.
      endif
!hchuang WRITE(*,'(''wrtg3d_hyb: CHECK RTIME = '', E12.2, TR2, ''FHOUR = '',
!hchuang & F10.2, TR2, ''ZHOUR = '', F10.2)') RTIME, FHOUR, ZHOUR
      IF(SECSWR.GT.0.) THEN
        RTIMSW=1./SECSWR
      ELSE
        RTIMSW=1.
      ENDIF
      cl1=colat1
!..........................................................
!
!     Moisture tendencies -- use id 249
!     GFS-GOCART NEED THE SUM OF FIRST 4 - 249,243,245,217
!
!     do nv=1,5+pl_coeff
       do k=1,levs
         il1 = k
         il2 = il1
        glolal=0.
        do lan=1,lats_node_r
         lat = global_lats_r(ipt_lats_node_r-1+lan)
         lons_lat = lonsperlar(lat)
         do i=1,lons_lat
          do nv=1,4
            glolal(i,lan) = glolal(i,lan) + dq3dt(i,k,nv,lan)*rtime
          enddo
         enddo
        enddo
        call uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
        call unsplit2d(ioproc,wrkga,buffo,global_lats_r)
        if(me.eq.ioproc) then
         call gribit(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &            0,idq3(1),isglev,il1,il2,iyr,imo,ida,ihr,
     &            ifhour,ifhr,ithr,iavg,0,0,icen2,ids(idq3(1)),iens,
     &            0.,0.,0.,0.,0.,0.,g,lg,ierr)
         if(ierr.ne.0)print*,'wrtg3d_hyb gribit ierr=',ierr,'  ',
     &   ' k=',k
!    x   '01)moisture tendency (g/kg/s) '
        endif
        if(ierr.eq.0 .and. me.eq.ioproc) call wryte(nog3d,lg,g)
       enddo
!     enddo
!
!     cloud cover
!
      do k=1,levs
        il1 = k
        il2 = il1
         do j=1,lats_node_r
           do i=1,lonr
!            glolal(i,j) = cldcov(k,i,j) * (rtimsw * 100.0)
            glolal(i,j) = cldcov(i,k,j) * rtimsw
          enddo
        enddo
        call uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
        call unsplit2d(ioproc,wrkga,buffo,global_lats_r)
        if(me.eq.ioproc) then
         call gribit(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &            0,icc,isglev,il1,il2,iyr,imo,ida,ihr,
     &            ifhour,ifhr,ithr,iavg,0,0,icen2,ids(icc),iens,
     &            0.,0.,0.,0.,0.,0.,g,lg,ierr)
         if(ierr.ne.0)print*,'wrtg3d_hyb gribit ierr=',ierr,'  ',' k=',k
!    x   '02) 3D cloud cover (fraction)                                 '
        endif

        if(ierr.eq.0 .and. me.eq.ioproc) call wryte(nog3d,lg,g)
      enddo
!
! hchuang 03/07/08 Correct the icmf sequence icmf(1) is updraft
!                                            icmf(2) is downdraft
!
!     Convective updraft Massflux      Table 133 parameter ID 202
!
       do k=1,levs
        il1 = k
        il2 = il1
        glolal=0.
        do lan=1,LATS_NODE_R
         lat = global_lats_r(ipt_lats_node_r-1+lan)
         lons_lat = lonsperlar(lat)
         DO i=1,lons_lat
           glolal(i,lan) = upd_mf(i,k,lan)* RTIME
         enddo
        enddo
        CALL uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
        call unsplit2d(ioproc,wrkga,buffo,global_lats_r)
        if(me.eq.ioproc) then
         call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,133,ICEN,IGEN,
     &            0,ICMF(1),ISGLEV,il1,il2,IYR,IMO,IDA,IHR,
     &            IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(icmf(1)),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
         if(ierr.ne.0)print*,'wrtg3d_hyb gribit ierr=',ierr,'  ',' K=',K
!    x   '04)Convective Updraft Massflux (Pa/s) '
        endif

        IF(IERR.EQ.0 .and. me.eq.ioproc) CALL WRYTE(nog3d,LG,G)
      enddo
!
!     Convective downdraft Massflux      Table 133 parameter ID 209
!
       do k=1,levs
        il1 = k
        il2 = il1
        glolal=0.
        do lan=1,LATS_NODE_R
         lat = global_lats_r(ipt_lats_node_r-1+lan)
         lons_lat = lonsperlar(lat)
         DO i=1,lons_lat
           glolal(i,lan) = dwn_mf(i,k,lan)* RTIME
         enddo
        enddo
        CALL uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
        call unsplit2d(ioproc,wrkga,buffo,global_lats_r)
        if(me.eq.ioproc) then
         call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,133,ICEN,IGEN,
     &            0,ICMF(2),ISGLEV,il1,il2,IYR,IMO,IDA,IHR,
     &            IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(icmf(2)),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
         if(ierr.ne.0)print*,'wrtg3d_hyb gribit ierr=',ierr,'  ',' K=',K
!    x   '03)Convective Downdraft Massflux (Pa/s) '
        endif

        IF(IERR.EQ.0 .and. me.eq.ioproc) CALL WRYTE(nog3d,LG,G)
      enddo
!
!     Convective Detrainment Massflux      Table 133 parameter ID 219
!
       do k=1,levs
        il1 = k
        il2 = il1
        glolal=0.
        do lan=1,LATS_NODE_R
         lat = global_lats_r(ipt_lats_node_r-1+lan)
         lons_lat = lonsperlar(lat)
         DO i=1,lons_lat
           glolal(i,lan) = det_mf(i,k,lan)* RTIME
         enddo
        enddo
        CALL uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
        call unsplit2d(ioproc,wrkga,buffo,global_lats_r)
        if(me.eq.ioproc) then
         call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,133,ICEN,IGEN,
     &            0,ICMF(3),ISGLEV,il1,il2,IYR,IMO,IDA,IHR,
     &            IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(icmf(3)),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
         if(ierr.ne.0)print*,'wrtg3d_hyb gribit ierr=',ierr,'  ',' K=',K
!    x   '05)Convective Detrainment Massflux (Pa/s) '
        endif

        IF(IERR.EQ.0 .and. me.eq.ioproc) CALL WRYTE(nog3d,LG,G)
      enddo
!
!     Vertical Diffusion Diffusivity      Table 129 parameter ID 182
!
       do k=1,levs
        il1 = k
        il2 = il1
        glolal=0.
        do lan=1,LATS_NODE_R
         lat = global_lats_r(ipt_lats_node_r-1+lan)
         lons_lat = lonsperlar(lat)
         DO i=1,lons_lat
           glolal(i,lan) = dkh(i,k,lan)* RTIME
         enddo
        enddo
        CALL uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
        call unsplit2d(ioproc,wrkga,buffo,global_lats_r)
        if(me.eq.ioproc) then
         call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,129,ICEN,IGEN,
     &            0,idkh,ISGLEV,il1,il2,IYR,IMO,IDA,IHR,
     &            IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(idkh),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
         if(ierr.ne.0)print*,'wrtg3d_hyb gribit ierr=',ierr,'  ',' K=',K
!    x   '06)Vertical Diffusion Diffusivity (m**2/s) '
        endif

        IF(IERR.EQ.0 .and. me.eq.ioproc) CALL WRYTE(nog3d,LG,G)
      enddo
!
!     large-scale rain production rate in kg/kg/s
!     Table 133 parameter ID 134
!
       do k=1,levs
        il1 = k
        il2 = il1
        glolal=0.
        do lan=1,LATS_NODE_R
         lat = global_lats_r(ipt_lats_node_r-1+lan)
         lons_lat = lonsperlar(lat)
         DO i=1,lons_lat
           glolal(i,lan) = rnp(i,k,lan)* RTIME ! unit in kg/kg/s
           if ( glolal(i,lan) < 0. ) glolal(i,lan) = 0.
         enddo
        enddo
        CALL uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
        call unsplit2d(ioproc,wrkga,buffo,global_lats_r)
        if(me.eq.ioproc) then
         call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,133,ICEN,IGEN,
     &            0,irnp,ISGLEV,il1,il2,IYR,IMO,IDA,IHR,
     &            IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(irnp),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
         if(ierr.ne.0)print*,'wrtg3d_hyb gribit ierr=',ierr,'  ',' K=',K
!    x   '07)3D LargeG-Scale Precip rate (kg/kg/s) '
        endif

        IF(IERR.EQ.0 .and. me.eq.ioproc) CALL WRYTE(nog3d,LG,G)
      enddo

!
      return
      end
