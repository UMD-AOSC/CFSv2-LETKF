MODULE letkf_tools
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
  use letkf_local
  USE efso_nml
  USE efso_tools
  USE sigio_module

  IMPLICIT NONE

  PRIVATE
  PUBLIC ::  das_letkf, das_efso

!  INTEGER,SAVE :: nobstotal

!  real(r_size),parameter :: infl_rtps = 0.85
!  REAL(r_size),PARAMETER :: cov_infl_mul = 1.00d0
! > 0: globally constant covariance inflation
! < 0: 3D inflation values input from a GPV file "infl_mul.grd"
  REAL(r_size),PARAMETER :: sp_infl_add = 0.0d0 !additive inflation

  INTEGER, PARAMETER :: lev_update_q = 30 !q and qc are only updated below and equal to this model level
  REAL(r_size), PARAMETER :: q_sprd_max = 0.5 !GYL, maximum q (ensemble spread)/(ensemble mean)

  REAL(r_size),PARAMETER :: var_local(nv3d+nv2d,6) = RESHAPE( (/ &
!       U    V    T    Q   QC   PS   
   & 1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,& ! U,V
   & 1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,& ! T,Tv
   & 1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,& ! Q,RH
   & 1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,& ! PS
   & 1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,& ! RAIN
   & 1.d0,1.d0,1.d0,1.d0,1.d0,1.d0 & ! TC
   & /),(/nv3d+nv2d,6/))
  INTEGER,SAVE :: var_local_n2n(nv3d+nv2d)

CONTAINS
!-----------------------------------------------------------------------
! Data Assimilation
!-----------------------------------------------------------------------
SUBROUTINE das_letkf(gues3d,gues2d,anal3d,anal2d)
  IMPLICIT NONE
  CHARACTER(12) :: inflfile='infl_mul.grd'
  REAL(r_size),INTENT(INOUT) :: gues3d(nij1,nlev,nbv,nv3d) ! background ensemble
  REAL(r_size),INTENT(INOUT) :: gues2d(nij1,nbv,nv2d)      !  output: destroyed
  REAL(r_size),INTENT(OUT) :: anal3d(nij1,nlev,nbv,nv3d) ! analysis ensemble
  REAL(r_size),INTENT(OUT) :: anal2d(nij1,nbv,nv2d)
  REAL(r_size),ALLOCATABLE :: mean3d(:,:,:)
  REAL(r_size),ALLOCATABLE :: mean2d(:,:)
  REAL(r_size),ALLOCATABLE :: hdxf(:,:)
  REAL(r_size),ALLOCATABLE :: rdiag(:)
  REAL(r_size),ALLOCATABLE :: rloc(:)
  REAL(r_size),ALLOCATABLE :: dep(:)
  REAL(r_size),ALLOCATABLE :: work3d(:,:,:)
  REAL(r_size),ALLOCATABLE :: work2d(:,:)
  REAL(r_sngl),ALLOCATABLE :: work3dg(:,:,:,:)
  REAL(r_sngl),ALLOCATABLE :: work2dg(:,:,:)
  REAL(r_size),ALLOCATABLE :: tmptv(:,:)
  REAL(r_size),ALLOCATABLE :: pfull(:,:)
  REAL(r_size) :: parm
  REAL(r_size) :: trans(nbv,nbv,nv3d+nv2d)
  REAL(r_size) :: q_mean,q_sprd  ! GYL
  REAL(r_size) :: q_anal(nbv)    ! GYL
  LOGICAL :: ex
  INTEGER :: ij,ilev,n,m,i,j,k,nobsl,ierr,iret

  real(r_size), allocatable :: sprdg3d(:,:,:)
  real(r_size), allocatable :: sprdg2d(:,:)
  real(r_size), allocatable :: sprda3d(:,:,:)
  real(r_size), allocatable :: sprda2d(:,:)
  real(r_size), allocatable :: meana3d(:,:,:)
  real(r_size), allocatable :: meana2d(:,:)  

  WRITE(6,'(A)') 'Hello from das_letkf'
  WRITE(6,'(A,F15.2)') '  cov_infl_mul = ',infl_mult
  nobstotal = nobs
  WRITE(6,'(A,I8)') 'Target observation numbers : NOBS=',nobs
  !
  ! In case of no obs
  !
  IF(nobstotal == 0) THEN
    WRITE(6,'(A)') 'No observation assimilated'
    anal3d = gues3d
    anal2d = gues2d
    RETURN
  END IF
  !
  ! Variable localization
  !
  var_local_n2n(1) = 1
  DO n=2,nv3d+nv2d
    DO i=1,n
      var_local_n2n(n) = i
      IF(MAXVAL(ABS(var_local(i,:)-var_local(n,:))) < TINY(var_local)) EXIT
    END DO
  END DO
  !
  ! FCST PERTURBATIONS
  !
  ALLOCATE(mean3d(nij1,nlev,nv3d))
  ALLOCATE(mean2d(nij1,nv2d))
  CALL ensmean_grd(nbv,nij1,gues3d,gues2d,mean3d,mean2d)
  DO n=1,nv3d
    DO m=1,nbv
      DO k=1,nlev
        DO i=1,nij1
          gues3d(i,k,m,n) = gues3d(i,k,m,n) - mean3d(i,k,n)
        END DO
      END DO
    END DO
  END DO
  DO n=1,nv2d
    DO m=1,nbv
      DO i=1,nij1
        gues2d(i,m,n) = gues2d(i,m,n) - mean2d(i,n)
      END DO
    END DO
  END DO
  !
  ! multiplicative inflation
  !
  IF(infl_mult > 0.0d0) THEN ! fixed multiplicative inflation parameter
    ALLOCATE( work3d(nij1,nlev,nv3d) )
    ALLOCATE( work2d(nij1,nv2d) )
    work3d = infl_mult
    work2d = infl_mult
  END IF
  IF(infl_mult <= 0.0d0) THEN ! 3D parameter values are read-in
    ALLOCATE( work3dg(nlon,nlat,nlev,nv3d) )
    ALLOCATE( work2dg(nlon,nlat,nv2d) )
    ALLOCATE( work3d(nij1,nlev,nv3d) )
    ALLOCATE( work2d(nij1,nv2d) )
    IF(myrank == 0) THEN
      INQUIRE(FILE=inflfile,EXIST=ex)
      IF(ex) THEN
        WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is reading.. ',inflfile
        CALL read_grd4(inflfile,work3dg,work2dg,0)
      ELSE
        WRITE(6,'(2A)') '!!WARNING: no such file exist: ',inflfile
        work3dg = -1.0d0 * infl_mult
        work2dg = -1.0d0 * infl_mult
      END IF
    END IF
    CALL scatter_grd_mpi(0,work3dg,work2dg,work3d,work2d)
  END IF
  !
  ! p_full for background ensemble mean
  !
  ALLOCATE( tmptv(nij1,nlev) )
  ALLOCATE( pfull(nij1,nlev) )
  tmptv = mean3d(:,:,iv3d_t) * (1.0d0 + fvirt * mean3d(:,:,iv3d_q))
  call sigio_modprd(nij1,nij1,nlev,gfs_nvcoord,gfs_idvc,gfs_idsl, &
                    gfs_vcoord,iret,mean2d(:,iv2d_ps),tmptv,pm=pfull)
  DEALLOCATE(tmptv)
  !
  ! MAIN ASSIMILATION LOOP
  !
  ALLOCATE( hdxf(1:nobstotal,1:nbv),rdiag(1:nobstotal),rloc(1:nobstotal),dep(1:nobstotal) )
  DO ij=1,nij1
     !! doing everythin one column at a time should be faster   
    WRITE(6,*) 'ij',ij,"of",nij1
    DO ilev=1,nlev       
      DO n=1,nv3d
        IF(var_local_n2n(n) < n) THEN
          trans(:,:,n) = trans(:,:,var_local_n2n(n))
          work3d(ij,ilev,n) = work3d(ij,ilev,var_local_n2n(n))
        ELSE
          CALL obs_local(lon1(ij),lat1(ij),pfull(ij,ilev),n,hdxf,rdiag,rloc,dep,nobsl)
          parm = work3d(ij,ilev,n)
          CALL letkf_core(nobstotal,nobsl,hdxf,rdiag,rloc,dep,parm,trans(:,:,n))
          work3d(ij,ilev,n) = parm
        END IF
        IF((n == iv3d_q .OR. n == iv3d_qc) .AND. ilev > lev_update_q) THEN   ! GYL, do not update upper-level q,qc
          anal3d(ij,ilev,:,n) = mean3d(ij,ilev,n) + gues3d(ij,ilev,:,n)      ! GYL
        ELSE                                                                 ! GYL
          DO m=1,nbv                                                         ! GYL
            anal3d(ij,ilev,m,n) = mean3d(ij,ilev,n)
            DO k=1,nbv
              anal3d(ij,ilev,m,n) = anal3d(ij,ilev,m,n) &
                & + gues3d(ij,ilev,k,n) * trans(k,m,n)
            END DO
          END DO                                                             ! GYL
        END IF                                                               ! GYL
        IF(n == iv3d_q .AND. ilev <= lev_update_q) THEN                      ! GYL, limit the lower-level q spread
          q_mean = SUM(anal3d(ij,ilev,:,n)) / REAL(nbv,r_size)               ! GYL
          q_sprd = 0.0d0                                                     ! GYL
          DO m=1,nbv                                                         ! GYL
            q_anal(m) = anal3d(ij,ilev,m,n) - q_mean                         ! GYL
            q_sprd = q_sprd + q_anal(m)**2                                   ! GYL
          END DO                                                             ! GYL
          q_sprd = SQRT(q_sprd / REAL(nbv-1,r_size)) / q_mean                ! GYL
          IF(q_sprd > q_sprd_max) THEN                                       ! GYL
            DO m=1,nbv                                                       ! GYL
              anal3d(ij,ilev,m,n) = q_mean + q_anal(m) * q_sprd_max / q_sprd ! GYL
            END DO                                                           ! GYL
          END IF                                                             ! GYL
        END IF
      END DO
      IF(ilev == 1) THEN !update 2d variable at ilev=1
        DO n=1,nv2d
          IF(var_local_n2n(nv3d+n) < nv3d+n) THEN                  ! GYL
            trans(:,:,nv3d+n) = trans(:,:,var_local_n2n(nv3d+n))
            IF(var_local_n2n(nv3d+n) <= nv3d) THEN                 ! GYL
              work2d(ij,n) = work3d(ij,ilev,var_local_n2n(nv3d+n)) ! GYL
            ELSE                                                   ! GYL
              work2d(ij,n) = work2d(ij,var_local_n2n(nv3d+n)-nv3d) ! GYL
            END IF                                                 ! GYL
          ELSE
            CALL obs_local(lon1(ij),lat1(ij),pfull(ij,ilev),nv3d+n,hdxf,rdiag,rloc,dep,nobsl)
            parm = work2d(ij,n)
            CALL letkf_core(nobstotal,nobsl,hdxf,rdiag,rloc,dep,parm,trans(:,:,nv3d+n))
            work2d(ij,n) = parm
          END IF
          DO m=1,nbv
            anal2d(ij,m,n)  = mean2d(ij,n)
            DO k=1,nbv
              anal2d(ij,m,n) = anal2d(ij,m,n) + gues2d(ij,k,n) * trans(k,m,nv3d+n)
            END DO
          END DO
        END DO
      END IF
    END DO
  END DO
  DEALLOCATE(hdxf,rdiag,rloc,dep)
  !
  ! Compute analyses of observations (Y^a)
  !
  IF(obsanal_output) THEN
    call das_letkf_obs(work3dg,work2dg)
  END IF
  !
  ! Write updated inflation parameters
  !
  IF(infl_mult < 0.0d0) THEN
    CALL gather_grd_mpi(0,work3d,work2d,work3dg,work2dg)
    IF(myrank == 0) THEN
      WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing.. ',inflfile
      
      CALL write_grd4(inflfile,work3dg,work2dg,0)
    END IF
    DEALLOCATE(work3dg,work2dg,work3d,work2d)
  END IF
  !
  ! Additive inflation
  !
  IF(sp_infl_add > 0.0d0) THEN
    CALL read_ens_mpi('addi',nbv,gues3d,gues2d)
    ALLOCATE( work3d(nij1,nlev,nv3d) )
    ALLOCATE( work2d(nij1,nv2d) )
    CALL ensmean_grd(nbv,nij1,gues3d,gues2d,work3d,work2d)
    DO n=1,nv3d
      DO m=1,nbv
        DO k=1,nlev
          DO i=1,nij1
            gues3d(i,k,m,n) = gues3d(i,k,m,n) - work3d(i,k,n)
          END DO
        END DO
      END DO
    END DO
    DO n=1,nv2d
      DO m=1,nbv
        DO i=1,nij1
          gues2d(i,m,n) = gues2d(i,m,n) - work2d(i,n)
        END DO
      END DO
    END DO

    !! additive inflation
    DEALLOCATE(work3d,work2d)
    WRITE(6,'(A)') '===== Additive covariance inflation ====='
    WRITE(6,'(A,F10.4)') '  parameter:',sp_infl_add
    WRITE(6,'(A)') '========================================='
!    parm = 0.7d0
!    DO ilev=1,nlev
!      parm_infl_damp(ilev) = 1.0d0 + parm &
!        & + parm * REAL(1-ilev,r_size)/REAL(nlev_dampinfl,r_size)
!      parm_infl_damp(ilev) = MAX(parm_infl_damp(ilev),1.0d0)
!    END DO
    DO n=1,nv3d
      DO m=1,nbv
        DO ilev=1,nlev
          DO ij=1,nij1
            anal3d(ij,ilev,m,n) = anal3d(ij,ilev,m,n) &
              & + gues3d(ij,ilev,m,n) * sp_infl_add
          END DO
        END DO
      END DO
    END DO
    DO n=1,nv2d
      DO m=1,nbv
        DO ij=1,nij1
          anal2d(ij,m,n) = anal2d(ij,m,n) + gues2d(ij,m,n) * sp_infl_add
        END DO
      END DO
    END DO
  END IF

  !! RTPS
  !! ------------------------------------------------------------
  if (infl_rtps > 0.0d0) then
     write (6,*) "RTPS", infl_rtps
     ALLOCATE( meana3d(nij1,nlev,nv3d) )
     ALLOCATE( meana2d(nij1,nv2d) )
     ALLOCATE( sprdg3d(nij1,nlev,nv3d) )
     ALLOCATE( sprdg2d(nij1,nv2d) )
     ALLOCATE( sprda3d(nij1,nlev,nv3d) )
     ALLOCATE( sprda2d(nij1,nv2d) )

     !calculate the analysis mean
     CALL ensmean_grd(nbv,nij1,anal3d,anal2d,meana3d,meana2d)

     sprdg3d = 0
     sprdg2d = 0
     sprda3d = 0
     sprda2d = 0
     do m=1,nbv
        ! gues3d and gues2d have already had the mean removed earlier
        sprdg3d = sprdg3d + (gues3d(:,:,m,:))**2
        sprdg2d = sprdg2d + (gues2d(:,m,:))**2
        sprda3d = sprda3d + (anal3d(:,:,m,:)-meana3d)**2
        sprda2d = sprda2d + (anal2d(:,m,:)-meana2d)**2
     end do
     sprdg3d  = sqrt(sprdg3d / (nbv-1))
     sprdg2d  = sqrt(sprdg2d / (nbv-1))
     sprda3d = sqrt(sprda3d / (nbv-1))
     sprda2d = sqrt(sprda2d / (nbv-1))
     !expand the perturbations back toward the guess spread
     do m=1,nbv
        where (sprda3d > 0.0)
           anal3d(:,:,m,:) = (anal3d(:,:,m,:) - meana3d) * &
                (infl_rtps*(sprdg3d-sprda3d)/sprda3d   + 1)&
                + meana3d
        end where
        where (sprda2d > 0.0)
           anal2d(:,m,:)   = (anal2d(:,m,:)   - meana2d) * &
                (infl_rtps*(sprdg2d-sprda2d)/sprda2d   + 1)&
                + meana2d
        end where
     end do

     ! cleanup
     deallocate(sprdg3d, sprdg2d, sprda3d, sprda2d, meana3d, meana2d) 
  end if
  
  DEALLOCATE(mean3d,mean2d)
  DEALLOCATE(pfull)
  RETURN
END SUBROUTINE das_letkf
!-----------------------------------------------------------------------
! Data assimilation for observations: Compute analyses of observations (Y^a)
! * currently only support multiplicative and adaptive inflation
!  -- 01/01/2014, Guo-Yuan Lien, 
!-----------------------------------------------------------------------
SUBROUTINE das_letkf_obs(v3dinfl,v2dinfl)
  IMPLICIT NONE
  REAL(r_sngl),INTENT(IN) :: v3dinfl(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(IN) :: v2dinfl(nlon,nlat,nv2d)
  REAL(r_size),ALLOCATABLE :: v3dinflx(:,:,:,:)
  REAL(r_size),ALLOCATABLE :: v2dinflx(:,:,:)
  REAL(r_size),ALLOCATABLE :: v3dtmp(:,:,:,:)
  REAL(r_size),ALLOCATABLE :: v2dtmp(:,:,:)
  REAL(r_size),ALLOCATABLE :: tmpps(:)
  REAL(r_size),ALLOCATABLE :: tmptv(:,:)
  REAL(r_size),ALLOCATABLE :: tmpp(:,:)
  REAL(r_size),ALLOCATABLE :: obsanal(:,:)
  REAL(r_size),ALLOCATABLE :: obsanalmean(:)
  REAL(r_size) :: hdxf(nobstotal,nbv)
  REAL(r_size) :: rdiag(nobstotal)
  REAL(r_size) :: rloc(nobstotal)
  REAL(r_size) :: dep(nobstotal)
  REAL(r_size) :: ohx(nobs)
  REAL(r_size) :: parm
  REAL(r_size) :: trans(nbv,nbv)
  REAL(r_size) :: ri,rj,rk
  REAL(r_size) :: rlev,p_update_q
  REAL(r_size) :: q_sprd
  REAL(r_size) :: q_anal(nbv)
  INTEGER :: n,nn,m,k,nobsl,ierr,iret
  INTEGER :: inflelem,irank,nobsp,nobspmax
  CHARACTER(14) :: obsanalfile='obsanalNNN.dat'

  WRITE(6,'(A)') 'Hello from das_letkf_obs: Compute [Y^a]'
  !
  ! If adaptive inflation is used, prepare a global array of inflation parameter
  !
  IF(infl_mult <= 0.0d0) THEN
    ALLOCATE(v3dinflx(nlon,nlat,nlev,nv3dx))
    ALLOCATE(v2dinflx(nlon,nlat,nv2dx))
    IF(myrank == 0) THEN
      ALLOCATE(v3dtmp(nlon,nlat,nlev,nv3d))
      ALLOCATE(v2dtmp(nlon,nlat,nv2d))
      ALLOCATE(tmpps(nlon*nlat))
      ALLOCATE(tmptv(nlon*nlat,nlev))
      ALLOCATE(tmpp(nlon*nlat,nlev))
      CALL read_grd('gues_me.grd',v3dtmp,v2dtmp,0)  ! read ensemble mean into a temporary array
      CALL read_grdx('gues001.grd',v3dinflx,v2dinflx) ! only the orography is used, P will be recalulated
      v3dinflx(:,:,:,iv3d_u) = v3dinfl(:,:,:,iv3d_u)
      v3dinflx(:,:,:,iv3d_v) = v3dinfl(:,:,:,iv3d_v)
      v3dinflx(:,:,:,iv3d_t) = v3dinfl(:,:,:,iv3d_t)
      v3dinflx(:,:,:,iv3d_q) = v3dinfl(:,:,:,iv3d_q)
      v3dinflx(:,:,:,iv3d_qc) = v3dinfl(:,:,:,iv3d_qc)
!      v2dinflx(:,:,iv2d_ps) = v2dinfl(:,:,iv2d_ps)
      v2dinflx(:,:,iv2d_ps) = v3dinfl(:,:,1,iv3d_u)
      tmpps = reshape(v2dtmp(:,:,iv2d_ps),(/nlon*nlat/))
      tmptv = reshape(v3dtmp(:,:,:,iv3d_t) * (1.0d0 + fvirt * v3dtmp(:,:,:,iv3d_q)),(/nlon*nlat,nlev/))
      call sigio_modprd(nlon*nlat,nlon*nlat,nlev,gfs_nvcoord,gfs_idvc,gfs_idsl, &
                        gfs_vcoord,iret,tmpps,tmptv,pm=tmpp)
      v3dinflx(:,:,:,iv3d_p) = reshape(tmpp,(/nlon,nlat,nlev/))
      DEALLOCATE(v3dtmp,v2dtmp,tmpps,tmptv,tmpp)
    END IF
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_BCAST(v3dinflx,nlon*nlat*nlev*nv3dx,MPI_r_size,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(v2dinflx,nlon*nlat*nv2dx,MPI_r_size,0,MPI_COMM_WORLD,ierr)
  END IF
  !
  ! Define the partition of observations for parallel computation
  !
  nn = MOD(nobs,nprocs)
  nobspmax = (nobs - nn)/nprocs + 1
  IF(myrank < nn) THEN
    nobsp = nobspmax
  ELSE
    nobsp = nobspmax-1
  END IF
  WRITE(6,'(A,I3.3,A,I8)') 'MYRANK ',myrank,' process obs number=', nobsp
  !
  ! Main LETKF loop
  !
  ALLOCATE(obsanal(nobs,nbv))
  ALLOCATE(obsanalmean(nobs))
  obsanal = 0.0d0
  obsanalmean = 0.0d0
  nn = myrank+1
  DO
    IF(nn > nobs) EXIT
!    WRITE(6,'(A,I8)') 'nn = ',nn
    !
    ! The observation variable type is different from the grid variable type.
    ! To compute the analyses of observations as regular grids,
    ! what grid variable will the observation variable be regarded as?
    !
    ! Also determine the pressure level will the observation variable be regarded as?
    !
    SELECT CASE(NINT(obselm(nn)))
    CASE(obsid_atm_u)
      n = iv3d_u          ! for variable localization, what grid variable to be regarded as? 
      inflelem = obsid_atm_u ! for inflation parameter,   what grid variable to be regarded as?
      rlev = obslev(nn)
    CASE(obsid_atm_v)
      n = iv3d_v
      inflelem = obsid_atm_v
      rlev = obslev(nn)
    CASE(obsid_atm_t,obsid_atm_tv)
      n = iv3d_t
      inflelem = obsid_atm_t
      rlev = obslev(nn)
    CASE(obsid_atm_q,obsid_atm_rh)
      n = iv3d_q
      inflelem = obsid_atm_q
      rlev = obslev(nn)
    CASE(obsid_atm_ps)
      n = nv3d+iv2d_ps
      inflelem = obsid_atm_ps
      rlev = obsdat(nn)   ! for ps variable, use the observed pressure value
    ! CASE(obsid_atm_rain)
    !   n = 0
    !   inflelem = obsid_atm_q
    !   rlev = base_obsv_rain ! for precipitation, assigh the level 'base_obsv_rain'
    CASE DEFAULT
      n = 0
      IF(NINT(obselm(nn)) > 9999) THEN
        inflelem = obsid_atm_ps
        CALL itpl_2d(v3dinflx(:,:,1,iv3d_p),ri,rj,rlev)
      ELSE
        inflelem = obsid_atm_u
        rlev = obslev(nn)
      END IF
    END SELECT
    !
    ! Determine the inflation parameter
    !
    IF(infl_mult > 0.0d0) THEN
      parm = infl_mult
    ELSE
      CALL phys2ijk(v3dinflx(:,:,:,iv3d_p),real(inflelem,r_size),obslon(nn),obslat(nn),rlev,ri,rj,rk)
      IF(CEILING(rk) > nlev) THEN
        rk = REAL(nlev,r_size)
      END IF
      IF(CEILING(rk) < 2 .AND. inflelem /= obsid_atm_ps) THEN
        IF(inflelem > 9999) THEN
          rk = 0.0d0
        ELSE
          rk = 1.00001d0
        END IF
      END IF
      IF(inflelem == obsid_atm_ps) THEN
        CALL itpl_2d(v2dinflx(:,:,iv2d_orog),ri,rj,rk)
        rk = obslev(nn) - rk
      END IF
      CALL Trans_XtoY(real(inflelem,r_size),ri,rj,rk,v3dinflx,v2dinflx,parm)
    END IF
    !
    ! LETKF computation
    !
    CALL obs_local(obslon(nn),obslat(nn),rlev,n,hdxf,rdiag,rloc,dep,nobsl)
    CALL letkf_core(nobstotal,nobsl,hdxf,rdiag,rloc,dep,parm,trans)

    IF(n == iv3d_q .OR. n == iv3d_qc) THEN
      CALL itpl_2d(v3dinflx(:,:,lev_update_q,iv3d_p),ri,rj,p_update_q)
    END IF
    IF((n == iv3d_q .OR. n == iv3d_qc) .AND. obslev(nn) < p_update_q) THEN
      obsanal(nn,:) = obsdat(nn) - obsdep(nn) + obshdxf(nn,:)
      obsanalmean(nn) = obsdat(nn) - obsdep(nn)
    ELSE
      DO m=1,nbv
        obsanal(nn,m) = obsdat(nn) - obsdep(nn)
        DO k=1,nbv
          obsanal(nn,m) = obsanal(nn,m) + obshdxf(nn,k) * trans(k,m)
        END DO
        obsanalmean(nn) = obsanalmean(nn) + obsanal(nn,m)
      END DO
      obsanalmean(nn) = obsanalmean(nn) / real(nbv,r_size)
    END IF
    IF(n == iv3d_q .AND. obslev(nn) >= p_update_q) THEN
      q_sprd = 0.0d0
      DO m=1,nbv
        q_anal(m) = obsanal(nn,m) - obsanalmean(nn)
        q_sprd = q_sprd + q_anal(m)**2
      END DO
      q_sprd = SQRT(q_sprd / REAL(nbv-1,r_size)) / obsanalmean(nn)
      IF(q_sprd > q_sprd_max) THEN
        DO m=1,nbv
          obsanal(nn,m) = obsanalmean(nn) + q_anal(m) * q_sprd_max / q_sprd
        END DO
      END IF
    END IF

    nn = nn + nprocs
  END DO
  !
  ! MPI_REDUCE and output obsanalfiles
  !
  ! mean
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_REDUCE(obsanalmean,ohx,nobs,MPI_r_size,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  IF(myrank == 0) THEN
    WRITE(obsanalfile(8:10),'(A3)') '_me'
    WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing a file ',obsanalfile
    CALL write_obs2(obsanalfile,nobs,obselm,obslon,obslat,obslev, &
                    obsdat,obserr,obstyp,obsdif,ohx,obsqc,0)
  END IF
  ! members
  irank = 0
  DO m=1,nbv
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    CALL MPI_REDUCE(obsanal(:,m),ohx,nobs,MPI_r_size,MPI_SUM,irank,MPI_COMM_WORLD,ierr)
    IF(myrank == irank) THEN
      WRITE(obsanalfile(8:10),'(I3.3)') m
      WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing a file ',obsanalfile
      CALL write_obs2(obsanalfile,nobs,obselm,obslon,obslat,obslev, &
                      obsdat,obserr,obstyp,obsdif,ohx,obsqc,0)
    END IF
    irank = irank + 1
    IF(irank >= nprocs) irank = 0
  END DO

  DEALLOCATE(obsanal)
  IF(infl_mult <= 0.0d0) THEN
    DEALLOCATE(v3dinflx,v2dinflx)
  END IF
  RETURN
END SUBROUTINE das_letkf_obs
!-----------------------------------------------------------------------
! Subroutine for observation sensitivity computation
! Ported from Y.Ohta's SPEEDY-LETKF system by D.Hotta, 07/01/2013
! [ref: Eq.(6,7,9), Ota et al. 2013]
!-----------------------------------------------------------------------
! [INPUT]
!  gues3d,gues2d: xmean^g_0
!  fcst3d,fcst2d: C^(1/2)*X^f_t                    [(J/kg)^(1/2)]
!  fcer3d,fcer2d: C^(1/2)*[1/2(K-1)](e^f_t+e^g_t)  [(J/kg)^(1/2)]
! (save variables)
!  obshdxf:
! [OUTPUT]
!-----------------------------------------------------------------------
SUBROUTINE das_efso(gues3d,gues2d,fcst3d,fcst2d,fcer3d,fcer2d)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: gues3d(nij1,nlev,nv3d)     !
  REAL(r_size),INTENT(IN) :: gues2d(nij1,nv2d)          !
  REAL(r_size),INTENT(IN) :: fcst3d(nij1,nlev,nbv,nv3d) ! forecast ensemble
  REAL(r_size),INTENT(IN) :: fcst2d(nij1,nbv,nv2d)      !
  REAL(r_size),INTENT(IN) :: fcer3d(nij1,nlev,nv3d) ! forecast error
  REAL(r_size),INTENT(IN) :: fcer2d(nij1,nv2d)      !
  REAL(r_size),ALLOCATABLE :: hdxf(:,:)
  REAL(r_size),ALLOCATABLE :: hdxa_rinv(:,:)
  REAL(r_size),ALLOCATABLE :: rdiag(:)
  REAL(r_size),ALLOCATABLE :: rloc(:)
  REAL(r_size),ALLOCATABLE :: dep(:)
  REAL(r_size),ALLOCATABLE :: tmptv(:,:)
  REAL(r_size),ALLOCATABLE :: pfull(:,:)
  REAL(r_size),ALLOCATABLE :: djdy(:,:)
  REAL(r_size),ALLOCATABLE :: recbuf(:,:)
  REAL(r_size) :: work1(nterm,nbv)
  INTEGER,ALLOCATABLE :: oindex(:)
  INTEGER :: ij,k,ilev,m,nob,nobsl,ierr,iret,iterm

  WRITE(6,'(A)') 'Hello from das_obsense'
  nobstotal = nobs !+ ntvs
  WRITE(6,'(A,I8)') 'Target observation numbers : NOBS=',nobs!,', NTVS=',ntvs
  !
  ! In case of no obs
  !
  IF(nobstotal == 0) THEN
    WRITE(6,'(A)') 'No observation assimilated'
    RETURN
  END IF
  ALLOCATE(djdy(nterm,nobstotal))
  djdy = 0.0_r_size
  !
  ! p_full for background ensemble mean
  !
  ALLOCATE( tmptv(nij1,nlev) )
  ALLOCATE( pfull(nij1,nlev) )
  tmptv = gues3d(:,:,iv3d_t) * (1.0d0 + fvirt * gues3d(:,:,iv3d_q))
  call sigio_modprd(nij1,nij1,nlev,gfs_nvcoord,gfs_idvc,gfs_idsl, &
                    gfs_vcoord,iret,gues2d(:,iv2d_ps),tmptv,pm=pfull)
  DEALLOCATE(tmptv)
  !
  ! MAIN ASSIMILATION LOOP
  !
!$OMP PARALLEL PRIVATE(ij,ilev,k,hdxf,rdiag,rloc,dep,nobsl,oindex, &
!$                     work1,m,nob)
  ALLOCATE( hdxf(1:nobstotal,1:nbv),rdiag(1:nobstotal),rloc(1:nobstotal), &
       & dep(1:nobstotal) )
  ALLOCATE(oindex(1:nobstotal))
!--- For ILEV = 1 - NLEV
!$OMP DO SCHEDULE(DYNAMIC)
  DO ilev=1,nlev
    WRITE(6,'(A,I3)') 'ilev = ',ilev
    DO ij=1,nij1
      IF(ABS(locadv_rate) > TINY(locadv_rate)) THEN
        CALL obs_local(lon2(ij,ilev),lat2(ij,ilev),pfull(ij,ilev),0,hdxf,rdiag,rloc,dep,nobsl,oindex)
      ELSE
        CALL obs_local(lon1(ij),lat1(ij),pfull(ij,ilev),0,hdxf,rdiag,rloc,dep,nobsl,oindex)
      END IF
      IF( nobsl /= 0 ) THEN
        ! Forecast error
        work1 = 0.0_r_size
        DO k=1,nv3d
          SELECT CASE(k)
          CASE(iv3d_u,iv3d_v)
            iterm = 1
          CASE(iv3d_t)
            iterm = 2
          CASE(iv3d_q)
            iterm = 3
          CASE DEFAULT
            iterm = 0
          END SELECT
          IF(iterm > 0) THEN
            DO m=1,nbv
              work1(iterm,m) = work1(iterm,m) + fcst3d(ij,ilev,m,k) * fcer3d(ij,ilev,k)
            END DO
          END IF
        END DO
        IF(ilev == 1) THEN
          DO k=1,nv2d
            IF(k == iv2d_ps) THEN
              DO m=1,nbv
                work1(2,m) = work1(2,m) + fcst2d(ij,m,k) * fcer2d(ij,k)
              END DO
            END IF
          END DO
        END IF
        !!! work1: [1/2(K-1)](X^f_t)^T*C*(e^f_t+e^g_t)  [J/kg]
        ! Hdxa Rinv
        ALLOCATE(hdxa_rinv(nobsl,nbv))
        DO m=1,nbv
          DO nob=1,nobsl
            hdxa_rinv(nob,m) = hdxf(nob,m) / rdiag(nob) * rloc(nob)
          END DO
        END DO
        !!! hdxa_rinv: rho*R^(-1)*Y^a_0 = rho*R^(-1)*(H X^a_0)
        ! dJ/dy
        DO nob=1,nobsl
          DO m=1,nbv
            djdy(:,oindex(nob)) = djdy(:,oindex(nob)) + work1(:,m) * hdxa_rinv(nob,m)
          END DO
        END DO
        !!! djdy: [1/2(K-1)]rho*R^(-1)*Y^a_0*(X^f_t)^T*C*(e^f_t+e^g_t)
        DEALLOCATE(hdxa_rinv)
      END IF
    END DO
  END DO
!$OMP END DO
  DEALLOCATE(hdxf,rdiag,rloc,dep,oindex)
!$OMP END PARALLEL
  !
  ! Calculate observation sensitivity
  !
!$OMP PARALLEL PRIVATE(nob)
!$OMP DO
  DO nob=1,nobstotal
    obsense(:,nob) = djdy(:,nob) * obsdep(nob)
  END DO
  !!! obsense: delta e^{f-g}_t = [1/2(K-1)][y_o-H(xmean^b_0)]^T*rho*R^(-1)*Y^a_0*(X^f_t)^T*C*(e^f_t+e^g_t)
!$OMP END DO
!$OMP END PARALLEL
  ! Gather observation sensitivity informations to the root
  ALLOCATE(recbuf(nterm,nobstotal))
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_REDUCE(obsense(:,1:nobstotal),recbuf,nterm*nobstotal,MPI_r_size,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  IF(myrank == 0) obsense(:,1:nobstotal) = recbuf(:,:)
  DEALLOCATE(recbuf)
  DEALLOCATE(djdy)
  DEALLOCATE(pfull)
  RETURN
END SUBROUTINE das_efso
END MODULE letkf_tools
