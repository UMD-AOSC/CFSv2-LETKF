PROGRAM obsope
!=======================================================================
!
! [PURPOSE:] Main program of observation operator
!
! [HISTORY:]
!   04/03/2013 Takemasa Miyoshi  created
!   08/30/2013 Guo-Yuan Lien     modified for GFS model
!   09/28/2013 Guo-Yuan Lien     add accumulated precipitation obs
!
!=======================================================================
  USE common
  USE common_gfs
  USE common_obs_gfs
  use letkf_params

  IMPLICIT NONE
  INTEGER,PARAMETER :: nslots=7 ! number of time slots for 4D-LETKF
  INTEGER,PARAMETER :: nbslot=4 ! basetime slot
  REAL(r_size),PARAMETER :: hourslot=1.0d0 ! time interval between slots in hour
  CHARACTER(11) :: obsinfile='obsinTT.dat' !IN
  CHARACTER(11) :: obsinfile_mean='obsinmn.dat' !IN
  CHARACTER(10) :: guesfile='guesTT.grd'   !IN
  CHARACTER(10) :: obsoutfile='obsout.dat' !OUT
  REAL(r_size),ALLOCATABLE :: elem(:)
  REAL(r_size),ALLOCATABLE :: rlon(:)
  REAL(r_size),ALLOCATABLE :: rlat(:)
  REAL(r_size),ALLOCATABLE :: rlev(:)
  REAL(r_size),ALLOCATABLE :: odat(:)
  REAL(r_size),ALLOCATABLE :: oerr(:)
  REAL(r_size),ALLOCATABLE :: otyp(:)
  REAL(r_size),ALLOCATABLE :: tdif(:)
  REAL(r_size),ALLOCATABLE :: ohx(:)
  INTEGER,ALLOCATABLE :: oqc(:)
  REAL(r_size),allocatable :: v3d(:,:,:,:)
  REAL(r_size),allocatable :: v2d(:,:,:)
  REAL(r_size),allocatable :: pp(:,:)
  REAL(r_size),PARAMETER :: threshold_dz=500.0d0
  REAL(r_size) :: dz,tg,qg
  INTEGER :: nobslots(nslots+1)
  INTEGER :: nobs,maxnobs
  REAL(r_size) :: ri,rj,rk
  INTEGER :: n,islot
  LOGICAL :: firstwrite

  call params_init()
  
  CALL set_common_gfs

  DO islot=1,nslots
    WRITE(obsinfile(6:7),'(I2.2)') islot
    CALL get_nobs(obsinfile,7,nobslots(islot))
    WRITE(6,'(2A,I9,A)') obsinfile, ':', nobslots(islot), ' OBSERVATIONS'
  END DO
  CALL get_nobs(obsinfile_mean,7,nobslots(nslots+1))
  WRITE(6,'(2A,I9,A)') obsinfile_mean, ':', nobslots(nslots+1), ' OBSERVATIONS'
  nobs = SUM(nobslots)
  maxnobs = MAXVAL(nobslots)
  WRITE(6,'(A,I9,A)') 'TOTAL:      ', nobs, ' OBSERVATIONS'
  ALLOCATE( elem(maxnobs) )
  ALLOCATE( rlon(maxnobs) )
  ALLOCATE( rlat(maxnobs) )
  ALLOCATE( rlev(maxnobs) )
  ALLOCATE( odat(maxnobs) )
  ALLOCATE( oerr(maxnobs) )
  ALLOCATE( otyp(maxnobs) )
  ALLOCATE( tdif(maxnobs) )
  ALLOCATE( ohx(maxnobs) )
  ALLOCATE( oqc(maxnobs) )
  allocate( v3d(nlon,nlat,nlev,nv3dx) )
  allocate( v2d(nlon,nlat,nv2dx) )
  allocate( pp(nlon,nlat) )
  
  firstwrite = .TRUE.
  pp = 0.0d0

  DO islot=1,nslots
    IF(nobslots(islot) == 0) CYCLE
    WRITE(obsinfile(6:7),'(I2.2)') islot
    CALL read_obs(obsinfile,nobslots(islot),&
      & elem(1:nobslots(islot)),rlon(1:nobslots(islot)),&
      & rlat(1:nobslots(islot)),rlev(1:nobslots(islot)),&
      & odat(1:nobslots(islot)),oerr(1:nobslots(islot)),&
      & otyp(1:nobslots(islot)))
    tdif(1:nobslots(islot)) = REAL(islot-nbslot,r_size)*hourslot
    WRITE(guesfile(5:6),'(I2.2)') islot
    CALL read_grdx(guesfile,v3d,v2d)
    if (islot == 1 .or. islot == nslots) then
      pp = pp + 0.5d0 * hourslot * v2d(:,:,iv2d_tprcp)
    else
      pp = pp + hourslot * v2d(:,:,iv2d_tprcp)
    end if
    ohx=0.0d0
    oqc=0
    DO n=1,nobslots(islot)
      CALL phys2ijk(v3d(:,:,:,iv3d_p),elem(n),rlon(n),rlat(n),rlev(n),ri,rj,rk)

      if (nint(otyp(n)) == obsid_atm_ps) rk = 0.0d0

      
      IF(CEILING(ri) < 2 .OR. nlon+1 < CEILING(ri))  CYCLE
      IF(CEILING(rj) < 2 .OR. nlat < CEILING(rj))    CYCLE      
      IF(CEILING(rk) > nlev) CYCLE
      
      !! set the level to 0.0 for surface obs (Ps or special platforms)
      !! otherwise make sure the level is in range
      select case(nint(otyp(n)))
      case (8,9,11,15,19) !i'm sure theres more surface platforms...
         rk = 0.0d0
      case default
         if (ceiling(rk) < 2) cycle         
      end select

      if (nint(elem(n)) == obsid_atm_ps) then
         !! make sure the pressure's elevation is in range         
         if ( odat(n) < -100.0d0) cycle

         !! calculate pressure adjustment
         CALL itpl_2d(v2d(:,:,iv2d_orog),ri,rj,dz)
         rk = rlev(n) - dz
         IF(ABS(rk) > threshold_dz) cycle
      END IF

      
      !
      ! observational operator
      !
      CALL Trans_XtoY(elem(n),ri,rj,rk,v3d,v2d,ohx(n))
      oqc(n) = 1
   END DO

   
    IF(firstwrite) THEN
      CALL write_obs2(obsoutfile,nobslots(islot),&
        & elem(1:nobslots(islot)),rlon(1:nobslots(islot)),&
        & rlat(1:nobslots(islot)),rlev(1:nobslots(islot)),&
        & odat(1:nobslots(islot)),oerr(1:nobslots(islot)),&
        & otyp(1:nobslots(islot)),tdif(1:nobslots(islot)),&
        & ohx(1:nobslots(islot)),oqc(1:nobslots(islot)),0)
      firstwrite = .FALSE.
    ELSE
      CALL write_obs2(obsoutfile,nobslots(islot),&
        & elem(1:nobslots(islot)),rlon(1:nobslots(islot)),&
        & rlat(1:nobslots(islot)),rlev(1:nobslots(islot)),&
        & odat(1:nobslots(islot)),oerr(1:nobslots(islot)),&
        & otyp(1:nobslots(islot)),tdif(1:nobslots(islot)),&
        & ohx(1:nobslots(islot)),oqc(1:nobslots(islot)),1)
    END IF
  END DO

  !!TODO, need to fix rk for surface obs
  IF(nobslots(nslots+1) > 0) THEN
    v2d(:,:,iv2d_tprcp) = pp
    CALL read_obs(obsinfile_mean,nobslots(nslots+1),&
      & elem(1:nobslots(nslots+1)),rlon(1:nobslots(nslots+1)),&
      & rlat(1:nobslots(nslots+1)),rlev(1:nobslots(nslots+1)),&
      & odat(1:nobslots(nslots+1)),oerr(1:nobslots(nslots+1)),&
      & otyp(1:nobslots(nslots+1)))
    tdif(1:nobslots(nslots+1)) = 0.0d0
    ohx=0.0d0
    oqc=0
    DO n=1,nobslots(nslots+1)
      IF(elem(n) /= obsid_atm_rain) CYCLE
      CALL phys2ij(rlon(n),rlat(n),ri,rj)
      IF(CEILING(ri) < 2 .OR. nlon+1 < CEILING(ri)) THEN
!        WRITE(6,'(A)') '* X-coordinate out of range'
!        WRITE(6,'(A,F6.2,A,F6.2)') '*   ri=',ri,', rlon=',rlon(n)
        CYCLE
      END IF
      IF(CEILING(rj) < 2 .OR. nlat < CEILING(rj)) THEN
!        WRITE(6,'(A)') '* Y-coordinate out of range'
!       WRITE(6,'(A,F6.2,A,F6.2)') '*   rj=',rj,', rlat=',rlat(n)
        CYCLE
      END IF
      !
      ! observational operator
      !
      rk = 0.0d0
      CALL Trans_XtoY(elem(n),ri,rj,rk,v3d,v2d,ohx(n))
      oqc(n) = 1
    END DO
    IF(firstwrite) THEN
      CALL write_obs2(obsoutfile,nobslots(nslots+1),&
        & elem(1:nobslots(nslots+1)),rlon(1:nobslots(nslots+1)),&
        & rlat(1:nobslots(nslots+1)),rlev(1:nobslots(nslots+1)),&
        & odat(1:nobslots(nslots+1)),oerr(1:nobslots(nslots+1)),&
        & otyp(1:nobslots(nslots+1)),tdif(1:nobslots(nslots+1)),&
        & ohx(1:nobslots(nslots+1)),oqc(1:nobslots(nslots+1)),0)
      firstwrite = .FALSE.
    ELSE
      CALL write_obs2(obsoutfile,nobslots(nslots+1),&
        & elem(1:nobslots(nslots+1)),rlon(1:nobslots(nslots+1)),&
        & rlat(1:nobslots(nslots+1)),rlev(1:nobslots(nslots+1)),&
        & odat(1:nobslots(nslots+1)),oerr(1:nobslots(nslots+1)),&
        & otyp(1:nobslots(nslots+1)),tdif(1:nobslots(nslots+1)),&
        & ohx(1:nobslots(nslots+1)),oqc(1:nobslots(nslots+1)),1)
    END IF
  END IF

  DEALLOCATE( elem,rlon,rlat,rlev,odat,oerr,otyp,tdif,ohx,oqc )

END PROGRAM obsope
