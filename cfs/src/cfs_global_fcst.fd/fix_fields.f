
      SUBROUTINE fix_fields(
     &                  LONSPERLAR,GLOBAL_LATS_R,XLON,XLAT,sfc_fld,
     &                  nst_fld,HPRIME,JINDX1,JINDX2,DDY,OZPLIN,
     &                  CREAD,CREAD_NST)
!!     
      use machine , only : kind_rad
      use funcphys                         
      use module_progtm             
      use resol_def
      use namelist_def
      use layout1
      use ozne_def
      use Sfc_Flx_ESMFMod
      use Nst_Var_ESMFMod
      IMPLICIT NONE
!!     
      INTEGER NREAD, NREAD_NST
      TYPE(Sfc_Var_Data)  :: sfc_fld
      TYPE(Nst_Var_Data)  :: nst_fld
      CHARACTER (len=*)   :: CREAD
      CHARACTER (len=*)   :: CREAD_NST
      INTEGER JINDX1(LATS_NODE_R),JINDX2(LATS_NODE_R)
      REAL (KIND=KIND_RAD) DDY(LATS_NODE_R)
      REAL (KIND=KIND_RAD) HPRIME(LONR,NMTVR,LATS_NODE_R)

      INTEGER IOZONDP
      REAL (kind=kind_rad) OZPLIN(LATSOZP,LEVOZP,pl_coeff,timeoz)
     &,                    XLON(LONR,LATS_NODE_R)
     &,                    XLAT(LONR,LATS_NODE_R)
       
      INTEGER              GLOBAL_LATS_R(LATR)
      INTEGER                 LONSPERLAR(LATR)
      real, PARAMETER:: RLAPSE=0.65E-2
      real dt_warm
      integer needoro, i, j
!!     
      call gfuncphys
      if (lsm == 0) then ! For OSU LSM
         CALL GRDDF
         CALL GRDKT
      endif
!!     
      IOZONDP = 0
      if (ntoz .gt. 0) IOZONDP = 1
      NREAD   = 14
!     CREAD   = 'fort.14'
!     sfc_fld%ORO     = 0.
!     NEEDORO = 0
!
      NEEDORO = 1
      CALL read_mtn_hprim_oz(sfc_fld%SLMSK,HPRIME,NEEDORO,sfc_fld%ORO,
     &                       sfc_fld%oro_uf,
     &                       IOZONDP,OZPLIN, GLOBAL_LATS_R,LONSPERLAR)
!
      NEEDORO = 0
      if(.not.adiab)then
        if (fhini == fhrot) then
          if (me == 0) print *,' call read_sfc CREAD=',cread
          CALL read_sfc(sfc_fld,NEEDORO,NREAD,
     &                  CREAD,GLOBAL_LATS_R,LONSPERLAR)

          if (nst_fcst > 0) then
            if (me == 0) print *,' call read_nst nst_spinup : ',
     &                             nst_spinup
            nst_fld%slmsk = sfc_fld%slmsk
            if ( nst_spinup ) then
              CALL set_nst(sfc_fld%tsea,nst_fld)
            else
              NREAD_NST   = 15
              CALL read_nst(nst_fld,NREAD_NST,CREAD_NST,
     &                      GLOBAL_LATS_R,LONSPERLAR)
              if (  nst_fcst > 1 ) then
                do j = 1, lats_node_r
                  do i = 1, lonr
                    if ( sfc_fld%SLMSK(i,j) == 0.0 ) then
                      dt_warm = (nst_fld%xt(i,j)+nst_fld%xt(i,j))
     &                        /  nst_fld%xz(i,j)
                      sfc_fld%TSEA(i,j) = nst_fld%Tref(i,j)
     &                 + dt_warm - nst_fld%dt_cool(i,j)
     &                 - (sfc_fld%oro(i,j)-sfc_fld%oro_uf(i,j)) * rlapse
                    endif
                  enddo
                enddo
!
!    When AM and NST is not coupled, tsea (in surface file) ==> tref
!
              elseif (nst_fcst == 1) then
                nst_fld%tref = sfc_fld%tsea
              endif
!
!   Reset the non-water points, since no mask (sea ice) update for nstanl file
!   done at present
!
              call nst_reset_nonwater(sfc_fld%tsea,nst_fld)
            endif                        ! if ( nst_spinup ) then
          endif                          ! if ( nst_fcst > 0 ) then
        else
          if (me .eq. 0) print *,' call read_sfc_r CREAD=',cread
          CALL read_sfc_r(sfc_fld,NEEDORO,NREAD,
     &                    CREAD,GLOBAL_LATS_R,LONSPERLAR)

          if ( nst_fcst > 0 ) then
            nst_fld%slmsk = sfc_fld%slmsk
            NREAD_NST   = 15
            if (me .eq. 0) print *,' call read_nst_r CREAD=',cread_nst
            CALL read_nst_r(nst_fld,NREAD_NST,CREAD_NST,
     &                      GLOBAL_LATS_R,LONSPERLAR)
          endif
        endif
      endif
!     NEEDORO=1
!     CALL read_mtn_hprim_oz(sfc_fld%SLMSK,HPRIME,NEEDORO,sfc_fld%ORO,
!    &     IOZONDP,OZPLIN, GLOBAL_LATS_R,LONSPERLAR)
!      
      CALL SETINDXOZ(LATS_NODE_R,LATS_NODE_R,GLOBAL_LATS_R,
     &               JINDX1,JINDX2,DDY)
!      
      CALL LONLAT_PARA(GLOBAL_LATS_R,XLON,XLAT,LONSPERLAR)
!!     
      RETURN
      END
