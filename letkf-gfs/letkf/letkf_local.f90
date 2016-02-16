MODULE letkf_local
!=======================================================================
!
! [PURPOSE:] Module for LETKF with GFS
!
! [HISTORY:]
!   01/26/2009 Takemasa Miyoshi  created
!   10/04/2012 Guo-Yuan Lien     modified for GFS model
!   07/01/2013 Daisuke Hotta     ported EFSO code from Y.Ohta's code
!   01/01/2014 Guo-Yuan Lien     merged to GFS-LETKF main development
!
!=======================================================================
  use letkf_params
  USE common
  use common_obs
  USE common_mpi
  USE common_gfs
  USE common_mpi_gfs
  USE common_letkf
  USE letkf_obs
  use kdtree
  ! USE efso_nml
  ! USE efso_tools
  ! USE sigio_module

  IMPLICIT NONE

  public :: obs_local, nobstotal

  PRIVATE  

  integer,save :: nobstotal
  integer :: initialized = 0
  type(KD_ROOT) :: kdtree_root
  
!   INTEGER, PARAMETER :: lev_update_q = 30 !q and qc are only updated below and equal to this model level
!   REAL(r_size), PARAMETER :: q_sprd_max = 0.5 !GYL, maximum q (ensemble spread)/(ensemble mean)

!   REAL(r_size),PARAMETER :: var_local(nv3d+nv2d,6) = RESHAPE( (/ &
! !       U    V    T    Q   QC   PS   
!    & 1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,& ! U,V
!    & 1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,& ! T,Tv
!    & 1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,& ! Q,RH
!    & 1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,& ! PS
!    & 1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,& ! RAIN
!    & 1.d0,1.d0,1.d0,1.d0,1.d0,1.d0 & ! TC
!    & /),(/nv3d+nv2d,6/))
!   INTEGER,SAVE :: var_local_n2n(nv3d+nv2d)

CONTAINS

  
!-----------------------------------------------------------------------
! Project global observations to local
!     (hdxf_g,dep_g,rdiag_g) -> (hdxf,dep,rdiag)
! -- modified, using (rlon,rlat,rlev) instead of (ij,ilev), Guo-Yuan Lien
! -- optional oindex output, followed by D.Hotta
!-----------------------------------------------------------------------
SUBROUTINE obs_local(rlon,rlat,rlev,nvar,hdxf,rdiag,rloc,dep,nobsl,oindex)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: rlon,rlat,rlev
  INTEGER,INTENT(IN) :: nvar
  REAL(r_size),INTENT(OUT) :: hdxf(nobstotal,nbv)
  REAL(r_size),INTENT(OUT) :: rdiag(nobstotal)
  REAL(r_size),INTENT(OUT) :: rloc(nobstotal)
  REAL(r_size),INTENT(OUT) :: dep(nobstotal)
  INTEGER,INTENT(OUT) :: nobsl
  INTEGER,INTENT(OUT),OPTIONAL :: oindex(nobstotal)      ! DH
  REAL(r_size) :: dlon_zero,dlat_zero                    ! GYL
  REAL(r_size) :: minlon,maxlon,minlat,maxlat
  REAL(r_size) :: tmplon,tmplat,tmperr,tmpwgt(nlev)
  REAL(r_size) :: logrlev
  INTEGER,ALLOCATABLE:: nobs_use(:)
  INTEGER :: imin,imax,jmin,jmax,im,ichan
  INTEGER :: n,nn,iobs

  real(r_size) :: sigma_max_h_d0, dist_zero_h, sigma_h
  real(r_size) :: dist_zero_v, sigma_v
  integer :: idx(nobstotal), j, id
  real(r_size) :: dist(nobstotal), dlev
  real(r_size) :: loc_h, loc_t, loc_v
  
  !!------------------------------------------------------------
  !!initialize the KD search tree
  if (initialized == 0) then
     WRITE(6,*) "initializing the obs_local()"
     initialized = 1
     call kd_init(kdtree_root, obslon, obslat)
     write(6,*) "done constructing KD search tree"
     write(6,*) "nobstotal",nobstotal
  end if


  !! determine the max search radius for the initial obs search
  sigma_max_h_d0 = max(sigma_atm_h(1),sigma_ocn_h(1)) * sqrt(10.0/3.0) * 2

  
  !!------------------------------------------------------------
  !! query the KD tree
  call kd_search(kdtree_root, obslon, obslat, (/rlon, rlat/), &
       sigma_max_h_d0, idx, dist, nn)

  
  !! for each observation found in the radius, do the localization
  logrlev = log(rlev)  
  nobsl = 0
  do n = 1, nn
     loc_h = 0.0
     loc_v = 0.0
     loc_t = 0.0
     
     if(nobsl >= nobstotal) return

     !! determine domain specific parameters
     !! --------------------------------------------------
     id = obselm(idx(n))
     
     !! horizontal localization
     !! -----------------------------         
     if (id >= obsid_atm_min .and. id <= obsid_atm_max) then
        sigma_h = (1.0d0-abs(rlat)/90.0d0)*&
             (sigma_atm_h(1)-sigma_atm_h(2))+sigma_atm_h(2)
     else if (id >= obsid_ocn_min .and. id <= obsid_ocn_max) then
        sigma_h = (1.0d0-abs(rlat)/90.0d0)*&
             (sigma_ocn_h(1)-sigma_ocn_h(2))+sigma_ocn_h(2)
     else
        cycle
     end if
     dist_zero_h = sigma_h * sqrt(10.0/3.0) * 2
     if(dist(n) > dist_zero_h) cycle
     loc_h = (dist(n)/sigma_h)**2

     
     !! vertical localization into atmosphere
     !! ------------------------------
     if (id == obsid_atm_ps) then
        dlev = abs(log(obsdat(idx(n))) - logrlev)
        sigma_v = sigma_atm_v
     else if (id >= obsid_atm_min .and. id <= obsid_atm_max) then
        dlev = abs(log(obslev(idx(n))) - logrlev)
        sigma_v = sigma_atm_v        
     else if (id >= obsid_ocn_min .and. id <= obsid_ocn_max) then
        !!TODO, do this correctly
        dlev = abs(log(1013.0d2)-logrlev)
        sigma_v = sigma_ocn_v
     end if
     dist_zero_v = sigma_v * sqrt(10.0/3.0) * 2
     if (dlev > dist_zero_v) cycle
     loc_v = (dlev/sigma_v)**2        
               
     
     !! use the observation
     !! ------------------------------     
     nobsl = nobsl+1
     rloc(nobsl) = exp(-0.5d0 * (loc_h + loc_t + loc_v))
     hdxf(nobsl,:) = obshdxf(idx(n), :)
     dep(nobsl) = obsdep(idx(n))
     rdiag(nobsl) = obserr(idx(n))**2
  end do
  return
! !
! ! INITIALIZE
! !
!   IF( nobs > 0 ) THEN
!     ALLOCATE(nobs_use(nobs))
!   END IF
! !
! ! data search
! !
!   dlat_zero = MAX(dist_zero,dist_zero_rain) / pi / re * 180.0d0 ! GYL
!   dlon_zero = dlat_zero / COS(rlat*pi/180.0d0)                  ! GYL
!   minlon = rlon - dlon_zero
!   maxlon = rlon + dlon_zero
!   minlat = rlat - dlat_zero
!   maxlat = rlat + dlat_zero
!   IF(maxlon - minlon >= 360.0d0) THEN
!     minlon = 0.0d0
!     maxlon = 360.0d0
!   END IF

!   DO jmin=1,nlat-2
!     IF(minlat < lat(jmin+1)) EXIT
!   END DO
!   DO jmax=1,nlat-2
!     IF(maxlat < lat(jmax+1)) EXIT
!   END DO
!   nn = 1
!   IF(minlon >= 0 .AND. maxlon <= 360.0) THEN
!     DO imin=1,nlon-1
!       IF(minlon < lon(imin+1)) EXIT
!     END DO
!     DO imax=1,nlon-1
!       IF(maxlon < lon(imax+1)) EXIT
!     END DO
!     IF( nobs > 0 ) &
!     & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
!   ELSE IF(minlon >= 0 .AND. maxlon > 360.0) THEN
!     DO imin=1,nlon-1
!       IF(minlon < lon(imin+1)) EXIT
!     END DO
!     maxlon = maxlon - 360.0d0
!     IF(maxlon > 360.0d0) THEN
!       imin = 1
!       imax = nlon
!       IF( nobs > 0 ) &
!       & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
!     ELSE
!       DO imax=1,nlon-1
!         IF(maxlon < lon(imax+1)) EXIT
!       END DO
!       IF(imax > imin) THEN
!         imin = 1
!         imax = nlon
!         IF( nobs > 0 ) &
!         & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
!       ELSE
!         imin = 1
!         IF( nobs > 0 ) &
!         & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
!         DO imin=1,nlon-1
!           IF(minlon < lon(imin+1)) EXIT
!         END DO
!         imax = nlon
!         IF( nobs > 0 ) &
!         & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
!       END IF
!     END IF
!   ELSE IF(minlon < 0 .AND. maxlon <= 360.0d0) THEN
!     DO imax=1,nlon-1
!       IF(maxlon < lon(imax+1)) EXIT
!     END DO
!     minlon = minlon + 360.0d0
!     IF(minlon < 0) THEN
!       imin = 1
!       imax = nlon
!       IF( nobs > 0 ) &
!       & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
!     ELSE
!       DO imin=1,nlon-1
!         IF(minlon < lon(imin+1)) EXIT
!       END DO
!       IF(imin < imax) THEN
!         imin = 1
!         imax = nlon
!         IF( nobs > 0 ) &
!         & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
!       ELSE
!         imin = 1
!         IF( nobs > 0 ) &
!         & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
!         DO imin=1,nlon-1
!           IF(minlon < lon(imin+1)) EXIT
!         END DO
!         imax = nlon
!         IF( nobs > 0 ) &
!         & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
!       END IF
!     END IF
!   ELSE
!     maxlon = maxlon - 360.0d0
!     minlon = minlon + 360.0d0
!     IF(maxlon > 360.0 .OR. minlon < 0) THEN
!       imin = 1
!       imax = nlon
!       IF( nobs > 0 ) &
!       & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
!     ELSE
!       DO imin=1,nlon-1
!         IF(minlon < lon(imin+1)) EXIT
!       END DO
!       DO imax=1,nlon-1
!         IF(maxlon < lon(imax+1)) EXIT
!       END DO
!       IF(imin > imax) THEN
!         imin = 1
!         imax = nlon
!         IF( nobs > 0 ) &
!         & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
!       ELSE
!         IF( nobs > 0 ) &
!         & CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
!       END IF
!     END IF
!   END IF
!   nn = nn-1
!   IF(nn < 1) THEN
!     nobsl = 0
!     RETURN
!   END IF
! !
! ! CONVENTIONAL
! !
!   logrlev = LOG(rlev)
!   nobsl = 0
!   IF(nn > 0) THEN
!     DO n=1,nn
!       !
!       ! vertical localization
!       !
!       IF(NINT(obselm(nobs_use(n))) == obsid_atm_ps) THEN
!         dlev = ABS(LOG(obsdat(nobs_use(n))) - logrlev)
!         IF(dlev > dist_zerov) CYCLE
!       ELSE IF(NINT(obselm(nobs_use(n))) == obsid_atm_rain) THEN
!         dlev = ABS(LOG(base_obsv_rain) - logrlev)
!         IF(dlev > dist_zerov_rain) CYCLE
!       ! ELSE IF(NINT(obselm(nobs_use(n))) >= id_tclon_obs) THEN !TC track obs
!       !   dlev = 0.0d0
!       ELSE !! other (3D) variables
!         dlev = ABS(LOG(obslev(nobs_use(n))) - logrlev)
!         IF(dlev > dist_zerov) CYCLE
!       END IF
!       !
!       ! horizontal localization
!       !
!       CALL com_distll_1(obslon(nobs_use(n)),obslat(nobs_use(n)),rlon,rlat,dist)
!       IF(NINT(obselm(nobs_use(n))) == obsid_atm_rain) THEN
!         IF(dist > dist_zero_rain) CYCLE
!       ELSE
!         IF(dist > dist_zero) CYCLE
!       END IF
!       !
!       ! variable localization
!       !
!       IF(nvar > 0) THEN ! use variable localization only when nvar > 0
!         SELECT CASE(NINT(obselm(nobs_use(n))))
!         CASE(obsid_atm_u, obsid_atm_v)
!           iobs=1
!         CASE(obsid_atm_t, obsid_atm_tv)
!           iobs=2
!         CASE(obsid_atm_q, obsid_atm_rh)
!           iobs=3
!         CASE(obsid_atm_ps)
!           iobs=4
!         CASE(obsid_atm_rain)
!           iobs=5
!         ! CASE(id_tclon_obs)
!         !   iobs=6
!         ! CASE(id_tclat_obs)
!         !   iobs=6
!         ! CASE(id_tcmip_obs)
!         !   iobs=6
!         END SELECT
!         IF(var_local(nvar,iobs) < TINY(var_local)) CYCLE
!       END IF

!       nobsl = nobsl + 1
!       hdxf(nobsl,:) = obshdxf(nobs_use(n),:)
!       dep(nobsl)    = obsdep(nobs_use(n))
!       !
!       ! Observational localization
!       !
!       tmperr=obserr(nobs_use(n))
!       rdiag(nobsl) = tmperr * tmperr
!       IF(NINT(obselm(nobs_use(n))) == obsid_atm_rain) THEN                              ! GYL
!         rloc(nobsl) =EXP(-0.5d0 * ((dist/sigma_obs_rain)**2 + (dlev/sigma_obsv)**2)) ! GYL
!       ELSE                                                                           ! GYL
!         rloc(nobsl) =EXP(-0.5d0 * ((dist/sigma_obs)**2 + (dlev/sigma_obsv)**2))      ! GYL
!       END IF                                                                         ! GYL
!       IF(nvar > 0) THEN ! use variable localization only when nvar > 0
!         rloc(nobsl) = rloc(nobsl) * var_local(nvar,iobs)
!       END IF
!       IF(PRESENT(oindex)) oindex(nobsl) = nobs_use(n)      ! DH
!     END DO
!   END IF
! !
!   IF( nobsl > nobstotal ) THEN
!     WRITE(6,'(A,I5,A,I5)') 'FATAL ERROR, NOBSL=',nobsl,' > NOBSTOTAL=',nobstotal
!     WRITE(6,*) 'LON,LAT,LEV,NN=', rlon,rlat,rlev,nn
!     STOP 99
!   END IF
! !
!   IF( nobs > 0 ) THEN
!     DEALLOCATE(nobs_use)
!   END IF
! !
!   RETURN
END SUBROUTINE obs_local

END MODULE letkf_local
