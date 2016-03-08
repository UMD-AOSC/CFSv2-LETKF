MODULE common_mpi_mom4
!=======================================================================
!
! [PURPOSE:] MPI procedures
!
! [ATTENTION:]
!   DO NOT COMPILE WITH BOTH INLINE EXPANSION AND OMP OPTIONS TOGETHER
!    (Use ONE if you want, but DON'T USE BOTH AT THE SAME TIME)
!
! [HISTORY:]
!   01/23/2009 Takemasa Miyoshi  created
!   04/26/2011 Steve Penny converted to OCEAN for use with MOM4
!   01/18/2015 Steve Penny converted for use with MOM6
!   09/11/2015 Steve Penny converted back for use with MOM4 due to MPI upgrades from MOM6 version
!
!=======================================================================
  USE common
  USE common_mpi
  USE common_mom4
  use letkf_mom_params, only : nbv
  IMPLICIT NONE
  PUBLIC

  INTEGER,PARAMETER :: mpibufsize=1000 !600 !(this worked as 'safe' with 480 procs on Gaea) !200 !1000  !STEVE: this fixes the problem of bad output when using over 6 nodes default=1000,mom2(mpich2)=200
  INTEGER,SAVE :: nij1                  !STEVE: this is the number of gridpoints to run on this (myrank) processor
  INTEGER,SAVE :: nij1max               !STEVE: the largest number of gridpoints on any 1 processor
  INTEGER,ALLOCATABLE,SAVE :: nij1node(:)
  REAL(r_size),ALLOCATABLE,SAVE :: phi1(:)
  REAL(r_size),ALLOCATABLE,SAVE :: kmt1(:)         !(OCEAN)
  REAL(r_size),ALLOCATABLE,SAVE :: dx1(:),dy1(:)
  REAL(r_size),ALLOCATABLE,SAVE :: lon1(:),lat1(:)
  REAL(r_size),ALLOCATABLE,SAVE :: i1(:),j1(:)     !(OCEAN) splits grid coordinates out into list like ijs

CONTAINS


SUBROUTINE set_common_mpi_mom4
  REAL(r_sngl), ALLOCATABLE :: v3dg(:,:,:,:) !(nlon,nlat,nlev,nv3d) != 0 !STEVE: initializing
  REAL(r_sngl), ALLOCATABLE :: v2dg(:,:,:) !(nlon,nlat,nv2d) != 0      !STEVE: initializing
  REAL(r_size),ALLOCATABLE :: v3d(:,:,:)
  REAL(r_size),ALLOCATABLE :: v2d(:,:)
  INTEGER :: i,j,n                                 !(OCEAN)
  INTEGER :: l, ll, mstart, mend, nv0 !(OCEAN)
  LOGICAL :: dodebug = .true.

  WRITE(6,'(A)') 'Hello from set_common_mpi_mom4'
  i = MOD(nlon*nlat,nprocs)
  nij1max = (nlon*nlat - i)/nprocs + 1
  WRITE(6,*) "nij1max = ", nij1max
  WRITE(6,*) "mpibufsize = ", mpibufsize

  IF(mpibufsize > nij1max) THEN
    WRITE(6,*) "mpibufsize > nij1max :: ", mpibufsize, nij1max
    WRITE(6,*) "using scatter/gather grd_mpi_fast."
  ELSE
    WRITE(6,*) "mpibufsize > nij1max :: ", mpibufsize, nij1max
    WRITE(6,*) "using scatter/gather grd_mpi_safe."
  END IF

  ! First, identify the number of gridpoints to use on this processor
  IF(myrank < i) THEN
    nij1 = nij1max
  ELSE
    nij1 = nij1max - 1
  END IF
  WRITE(6,'(A,I3.3,A,I6)') 'MYRANK ',myrank,' number of grid points: nij1= ',nij1
  ALLOCATE(nij1node(nprocs))
  DO n=1,nprocs
    IF(n-1 < i) THEN
      nij1node(n) = nij1max
    ELSE
      nij1node(n) = nij1max - 1
    END IF
  END DO

  ! Next, allocate arrays to hold grid, land/sea mask, etc. information
  if (dodebug) WRITE(6,*) "ALLOCATING fields to convert to vectorized form..."
  ALLOCATE(phi1(nij1))
  ALLOCATE(kmt1(nij1))               !(OCEAN)
  ALLOCATE(dx1(nij1))
  ALLOCATE(dy1(nij1))
  ALLOCATE(lon1(nij1))
  ALLOCATE(lat1(nij1))
  ALLOCATE(i1(nij1))                 !(OCEAN)
  ALLOCATE(j1(nij1))                 !(OCEAN)

  nv0=8
  ALLOCATE(v2d(nij1,nv0))
  ALLOCATE(v2dg(nlon,nlat,nv0))

  ! Distribute that data to the appropriate processors
  if (dodebug) WRITE(6,*) "Converting dx, dy, lon, lat, i, j, phi0, and kmt0 to vectorized form..."
  v2dg=0.0
  DO j=1,nlat
    ! 2D and 1D Data stored in first layer:
    v2dg(:,j,1) = SNGL(dx(:,j))
    v2dg(:,j,2) = SNGL(dy(:,j))
    v2dg(:,j,3) = SNGL(lon(:))
    v2dg(:,j,4) = SNGL(lat(j)) !(Single value promoted to array)
    ! 2D Data stored in second layer:
    v2dg(:,j,5) = SNGL(kmt0(:,j))
    v2dg(:,j,6) = SNGL(phi0(:,j))  !STEVE: WARNING - this needs to be generalized to the specified dimensions
    !STEVE: For custom localization: (need to know how the grid points are distributed per node)
    DO i=1,nlon                             !(OCEAN)
      v2dg(i,j,7) = REAL(i,r_sngl)        !(OCEAN)
    ENDDO                                   !(OCEAN)
    v2dg(:,j,8) = REAL(j,r_sngl)          !(OCEAN)
  END DO
  if (dodebug) WRITE(6,*) "Calling scatter_grd_mpi_small..."
  CALL scatter_grd_mpi_small(0,v2dg,v2d,nlon,nlat,nv0)
! if (dodebug) WRITE(6,*) "Calling scatter_grd_mpi_smalltoall..."
! ll = CEILING(REAL(nbv)/REAL(nprocs))
! DO l=1,ll
!   mstart = 1 + (l-1)*nprocs
!   mend = MIN(l*nprocs, nbv)
!   if (dodebug) WRITE(6,*) "l,ll,mstart,mend,nprocs,nbv = ", l,ll,mstart,mend,nprocs,nbv
!   CALL scatter_grd_mpi_smalltoall(mstart,mend,nbv,v2dg,v2d,nlon,nlat,nv0)
! END DO
! if (dodebug) WRITE(6,*) "Finished scatter_grd_mpi_smalltoall."
  dx1(:)  = v2d(:,1)
  dy1(:)  = v2d(:,2)
  lon1(:) = v2d(:,3)
  lat1(:) = v2d(:,4)
  kmt1(:) = v2d(:,5)               !(OCEAN)
  phi1(:) = v2d(:,6)               !(OCEAN)
  i1(:)   = v2d(:,7)                 !(OCEAN)
  j1(:)   = v2d(:,8)                 !(OCEAN)

! DEALLOCATE(v3d,v2d,v3dg,v2dg)
  DEALLOCATE(v2d,v2dg)
  if (dodebug) WRITE(6,*) "Finished set_common_mpi_mom4..."

  RETURN
END SUBROUTINE set_common_mpi_mom4

!-----------------------------------------------------------------------
! Scatter gridded data to processes (nrank -> all)
!-----------------------------------------------------------------------
SUBROUTINE scatter_grd_mpi(nrank,v3dg,v2dg,v3d,v2d)
  INTEGER,INTENT(IN) :: nrank
  REAL(r_sngl),INTENT(IN) :: v3dg(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(IN) :: v2dg(nlon,nlat,nv2d)
  REAL(r_size),INTENT(OUT) :: v3d(nij1,nlev,nv3d)
  REAL(r_size),INTENT(OUT) :: v2d(nij1,nv2d)
  LOGICAL :: dodebug = .true.

  IF(mpibufsize > nij1max) THEN
    if (dodebug) then
      WRITE(6,*) "scatter_grd_mpi: calling scatter_grd_mpi_fast. mpibufsize, nij1max = ", mpibufsize, nij1max
    endif
    CALL scatter_grd_mpi_fast(nrank,v3dg,v2dg,v3d,v2d)
  ELSE
    if (dodebug) then
      WRITE(6,*) "scatter_grd_mpi: calling scatter_grd_mpi_safe. mpibufsize, nij1max = ", mpibufsize, nij1max
    endif
    CALL scatter_grd_mpi_safe(nrank,v3dg,v2dg,v3d,v2d)
  END IF
  if (dodebug) WRITE(6,*) "Finished scatter_grd_mpi..."

  RETURN
END SUBROUTINE scatter_grd_mpi

SUBROUTINE scatter_grd_mpi_safe(nrank,v3dg,v2dg,v3d,v2d)
  INTEGER,INTENT(IN) :: nrank
  REAL(r_sngl),INTENT(IN) :: v3dg(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(IN) :: v2dg(nlon,nlat,nv2d)
  REAL(r_size),INTENT(OUT) :: v3d(nij1,nlev,nv3d)
  REAL(r_size),INTENT(OUT) :: v2d(nij1,nv2d)
  REAL(r_sngl), ALLOCATABLE :: tmp(:,:) !(nij1max,nprocs)
  REAL(r_sngl), ALLOCATABLE :: bufs(:,:) !(mpibufsize,nprocs)
  REAL(r_sngl), ALLOCATABLE :: bufr(:) !(mpibufsize)
  INTEGER :: i,j,k,n,ierr,ns,nr
  INTEGER :: iter,niter

  ALLOCATE(tmp(nij1max,nprocs),bufs(mpibufsize,nprocs),bufr(mpibufsize))

  ns = mpibufsize
  nr = ns
  niter = CEILING(REAL(nij1max)/REAL(mpibufsize))

  DO n=1,nv3d
    DO k=1,nlev
      IF(myrank == nrank) CALL grd_to_buf(v3dg(:,:,k,n),tmp)
      DO iter=1,niter
        IF(myrank == nrank) THEN
          i = mpibufsize * (iter-1)
          DO j=1,mpibufsize
            i=i+1
            IF(i > nij1) EXIT
            bufs(j,:) = tmp(i,:)
          END DO
        END IF
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        CALL MPI_SCATTER(bufs,ns,MPI_REAL,&
                       & bufr,nr,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)
        i = mpibufsize * (iter-1)
        DO j=1,mpibufsize
          i=i+1
          IF(i > nij1) EXIT
          v3d(i,k,n) = REAL(bufr(j),r_size)
        END DO
      END DO
    END DO
  END DO

  DO n=1,nv2d
    IF(myrank == nrank) CALL grd_to_buf(v2dg(:,:,n),tmp)
    DO iter=1,niter
      IF(myrank == nrank) THEN
        i = mpibufsize * (iter-1)
        DO j=1,mpibufsize
          i=i+1
          IF(i > nij1) EXIT
          bufs(j,:) = tmp(i,:)
        END DO
      END IF
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      CALL MPI_SCATTER(bufs,ns,MPI_REAL,&
                     & bufr,nr,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)
      i = mpibufsize * (iter-1)
      DO j=1,mpibufsize
        i=i+1
        IF(i > nij1) EXIT
        v2d(i,n) = REAL(bufr(j),r_size)
      END DO
    END DO
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  DEALLOCATE(tmp,bufs,bufr)

  RETURN
END SUBROUTINE scatter_grd_mpi_safe

SUBROUTINE scatter_grd_mpi_fast(nrank,v3dg,v2dg,v3d,v2d)
  INTEGER,INTENT(IN) :: nrank
  REAL(r_sngl),INTENT(IN) :: v3dg(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(IN) :: v2dg(nlon,nlat,nv2d)
  REAL(r_size),INTENT(OUT) :: v3d(nij1,nlev,nv3d)
  REAL(r_size),INTENT(OUT) :: v2d(nij1,nv2d)
  REAL(r_sngl), ALLOCATABLE :: bufs(:,:,:) !(nij1max,nlevall,nprocs)
  REAL(r_sngl), ALLOCATABLE :: bufr(:,:) !(nij1max,nlevall)
  INTEGER :: j,k,n,ierr,ns,nr

  ALLOCATE(bufs(nij1max,nlevall,nprocs), bufr(nij1max,nlevall))

  ns = nij1max * nlevall
  nr = ns
  IF(myrank == nrank) THEN
    j=0
    DO n=1,nv3d
      DO k=1,nlev
        j = j+1
        CALL grd_to_buf(v3dg(:,:,k,n),bufs(:,j,:))
      END DO
    END DO

    DO n=1,nv2d
      j = j+1
      CALL grd_to_buf(v2dg(:,:,n),bufs(:,j,:))
    END DO
  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_SCATTER(bufs,ns,MPI_REAL,&
                 & bufr,nr,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)

  j=0
  DO n=1,nv3d
    DO k=1,nlev
      j = j+1
      v3d(:,k,n) = REAL(bufr(1:nij1,j),r_size)
    END DO
  END DO

  DO n=1,nv2d
    j = j+1
    v2d(:,n) = REAL(bufr(1:nij1,j),r_size)
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  DEALLOCATE(bufs,bufr)

  RETURN
END SUBROUTINE scatter_grd_mpi_fast
!-----------------------------------------------------------------------
! Gather gridded data (all -> nrank)
!-----------------------------------------------------------------------
SUBROUTINE gather_grd_mpi(nrank,v3d,v2d,v3dg,v2dg)
  INTEGER,INTENT(IN) :: nrank
  REAL(r_size),INTENT(IN) :: v3d(nij1,nlev,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nij1,nv2d)
  REAL(r_sngl),INTENT(OUT) :: v3dg(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(OUT) :: v2dg(nlon,nlat,nv2d)
  LOGICAL :: dodebug = .true.

  IF(mpibufsize > nij1max) THEN
    if (dodebug) then
      WRITE(6,*) "gather_grd_mpi: calling gather_grd_mpi_fast. mpibufsize,nij1max = ", mpibufsize, nij1max
    endif
    CALL gather_grd_mpi_fast(nrank,v3d,v2d,v3dg,v2dg)
  ELSE
    if (dodebug) then
      WRITE(6,*) "gather_grd_mpi: calling gather_grd_mpi_safe. mpibufsize,nij1max = ", mpibufsize, nij1max
    endif
    CALL gather_grd_mpi_safe(nrank,v3d,v2d,v3dg,v2dg)
  END IF

  RETURN
END SUBROUTINE gather_grd_mpi

SUBROUTINE gather_grd_mpi_safe(nrank,v3d,v2d,v3dg,v2dg)
  INTEGER,INTENT(IN) :: nrank
  REAL(r_size),INTENT(IN) :: v3d(nij1,nlev,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nij1,nv2d)
  REAL(r_sngl),INTENT(OUT) :: v3dg(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(OUT) :: v2dg(nlon,nlat,nv2d)
  REAL(r_sngl), ALLOCATABLE :: tmp(:,:) !(nij1max,nprocs)
  REAL(r_sngl), ALLOCATABLE :: bufs(:) !(mpibufsize)
  REAL(r_sngl), ALLOCATABLE :: bufr(:,:) !(mpibufsize,nprocs)
  INTEGER :: i,j,k,n,ierr,ns,nr
  INTEGER :: iter,niter

  ALLOCATE(tmp(nij1max,nprocs),bufs(mpibufsize),bufr(mpibufsize,nprocs))

  ns = mpibufsize
  nr = ns
  niter = CEILING(REAL(nij1max)/REAL(mpibufsize))

  DO n=1,nv3d
    DO k=1,nlev
      DO iter=1,niter
        i = mpibufsize * (iter-1)
        DO j=1,mpibufsize
          i=i+1
          IF(i > nij1) EXIT
          bufs(j) = REAL(v3d(i,k,n),r_sngl)
        END DO
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        CALL MPI_GATHER(bufs,ns,MPI_REAL,&
                      & bufr,nr,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)
        IF(myrank == nrank) THEN
          i = mpibufsize * (iter-1)
          DO j=1,mpibufsize
            i=i+1
            IF(i > nij1) EXIT
            tmp(i,:) = bufr(j,:)
          END DO
        END IF
      END DO
      IF(myrank == nrank) CALL buf_to_grd(tmp,v3dg(:,:,k,n))
    END DO
  END DO

  DO n=1,nv2d
    DO iter=1,niter
      i = mpibufsize * (iter-1)
      DO j=1,mpibufsize
        i=i+1
        IF(i > nij1) EXIT
        bufs(j) = REAL(v2d(i,n),r_sngl)
      END DO
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      CALL MPI_GATHER(bufs,ns,MPI_REAL,&
                    & bufr,nr,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)
      IF(myrank == nrank) THEN
        i = mpibufsize * (iter-1)
        DO j=1,mpibufsize
          i=i+1
          IF(i > nij1) EXIT
          tmp(i,:) = bufr(j,:)
        END DO
      END IF
    END DO
    IF(myrank == nrank) CALL buf_to_grd(tmp,v2dg(:,:,n))
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  DEALLOCATE(tmp,bufs,bufr)

  RETURN
END SUBROUTINE gather_grd_mpi_safe

SUBROUTINE gather_grd_mpi_fast(nrank,v3d,v2d,v3dg,v2dg)
  INTEGER,INTENT(IN) :: nrank
  REAL(r_size),INTENT(IN) :: v3d(nij1,nlev,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nij1,nv2d)
  REAL(r_sngl),INTENT(OUT) :: v3dg(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(OUT) :: v2dg(nlon,nlat,nv2d)
  REAL(r_sngl), ALLOCATABLE :: bufs(:,:) !(nij1max,nlevall)
  REAL(r_sngl), ALLOCATABLE :: bufr(:,:,:) !(nij1max,nlevall,nprocs)
  INTEGER :: j,k,n,ierr,ns,nr

  ALLOCATE(bufs(nij1max,nlevall),bufr(nij1max,nlevall,nprocs))

  ns = nij1max * nlevall
  nr = ns
  j=0
  DO n=1,nv3d
    DO k=1,nlev
      j = j+1
      bufs(1:nij1,j) = REAL(v3d(:,k,n),r_sngl)
    END DO
  END DO

  DO n=1,nv2d
    j = j+1
    bufs(1:nij1,j) = REAL(v2d(:,n),r_sngl)
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_GATHER(bufs,ns,MPI_REAL,&
                & bufr,nr,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)

  IF(myrank == nrank) THEN
    j=0
    DO n=1,nv3d
      DO k=1,nlev
        j = j+1
        CALL buf_to_grd(bufr(:,j,:),v3dg(:,:,k,n))
      END DO
    END DO

    DO n=1,nv2d
      j = j+1
      CALL buf_to_grd(bufr(:,j,:),v2dg(:,:,n))
    END DO
  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  DEALLOCATE(bufs,bufr)

  RETURN
END SUBROUTINE gather_grd_mpi_fast

!-----------------------------------------------------------------------
! Read ensemble data and distribute to processes
!-----------------------------------------------------------------------
SUBROUTINE read_ens_mpi(file,member,v3d,v2d)
  CHARACTER(4),INTENT(IN) :: file
  INTEGER,INTENT(IN) :: member
  REAL(r_size),INTENT(OUT) :: v3d(nij1,nlev,member,nv3d)
  REAL(r_size),INTENT(OUT) :: v2d(nij1,member,nv2d)
  REAL(r_sngl), ALLOCATABLE :: v3dg(:,:,:,:) !(nlon,nlat,nlev,nv3d)
  REAL(r_sngl), ALLOCATABLE :: v2dg(:,:,:) !(nlon,nlat,nv2d)
  INTEGER :: l,n,ll,im,mstart,mend
  CHARACTER(11) :: filename='file000'
  !STEVE: debug
  LOGICAL :: dodebug = .true.

  ALLOCATE(v3dg(nlon,nlat,nlev,nv3d),v2dg(nlon,nlat,nv2d))

  ll = CEILING(REAL(member)/REAL(nprocs))
  DO l=1,ll
    im = myrank+1 + (l-1)*nprocs
    IF(im <= member) THEN
      WRITE(filename(1:7),'(A4,I3.3)') file,im
      WRITE(6,'(A,I3.3,2A)') 'In common_mpi_mom4.f90::read_ens_mpi, MYRANK ',myrank,' is reading a file ',filename
!     CALL read_restart(filename,v3dg,v2dg,2) !STEVE: 20150317, trying this out...
      CALL read_grd4(filename,v3dg,v2dg) !STEVE: 20130709, trying this out...
!     CALL read_grd(filename,v3dg,v2dg)  !STEVE: causes type problem in scatter_grd_mpi
!     if (.false.) CALL write_grd4('test.'//filename,v3dg,v2dg) 
    END IF

    mstart = 1 + (l-1)*nprocs
    mend = MIN(l*nprocs, member)
    if (dodebug) WRITE(6,*) "In common_mpi_mom4.f90::read_ens_mpi, calling scatter_grd_mpi_alltoall..."
    CALL scatter_grd_mpi_alltoall(mstart,mend,member,v3dg,v2dg,v3d,v2d)

  END DO

  DEALLOCATE(v3dg,v2dg)

END SUBROUTINE read_ens_mpi

!-----------------------------------------------------------------------
! Write ensemble data after collecting data from processes
!-----------------------------------------------------------------------
SUBROUTINE write_ens_mpi(file,member,v3d,v2d)
! INCLUDE 'netcdf.inc' !STEVE: for NaN correction (OCEAN)
  CHARACTER(4),INTENT(IN) :: file
  INTEGER,INTENT(IN) :: member
  REAL(r_size),INTENT(IN) :: v3d(nij1,nlev,member,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nij1,member,nv2d)
  REAL(r_sngl), ALLOCATABLE :: v3dg(:,:,:,:) !(nlon,nlat,nlev,nv3d)
  REAL(r_sngl), ALLOCATABLE :: v2dg(:,:,:) !(nlon,nlat,nv2d)
  INTEGER :: l,n,ll,im
  CHARACTER(11) :: filename='file000'
  INTEGER :: i,j,k,m !STEVE: for debugging
  LOGICAL :: verbose = .true.
  INTEGER :: convcnt = 0
  INTEGER :: mstart,mend

  ALLOCATE(v3dg(nlon,nlat,nlev,nv3d),v2dg(nlon,nlat,nv2d))

  ll = CEILING(REAL(member)/REAL(nprocs))
  DO l=1,ll
    mstart = 1 + (l-1)*nprocs
    mend = MIN(l*nprocs, member)
    CALL gather_grd_mpi_alltoall(mstart,mend,member,v3d,v2d,v3dg,v2dg)

    im = myrank+1 + (l-1)*nprocs
    IF(im <= member) THEN
      WRITE(filename(1:7),'(A4,I3.3)') file,im
      WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing file: ',filename

      !STEVE: debug
      WRITE(6,*) "common_mpi_mom4.f90::write_ens_mpi:: MAXVAL(ABS(v3dg(:,:,:,iv3d_t))) = ", MAXVAL(ABS(v3dg(:,:,:,iv3d_t)))
      WRITE(6,*) "common_mpi_mom4.f90::write_ens_mpi:: MAXVAL(ABS(v3dg(:,:,:,iv3d_s))) = ", MAXVAL(ABS(v3dg(:,:,:,iv3d_s)))

!     CALL write_restart(filename,v3dg,v2dg)
      CALL write_grd4(filename,v3dg,v2dg)
    END IF
  END DO

  DEALLOCATE(v3dg,v2dg)

  RETURN
END SUBROUTINE write_ens_mpi

!!!!!!!!!
!STEVE: for debugging grid:
!!!!!!!!!
SUBROUTINE write_ens_mpi_grd(file,member,v3d,v2d)
  CHARACTER(4),INTENT(IN) :: file
  INTEGER,INTENT(IN) :: member
  REAL(r_size),INTENT(IN) :: v3d(nij1,nlev,member,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nij1,member,nv2d)
  REAL(r_sngl), ALLOCATABLE :: v3dg(:,:,:,:) !(nlon,nlat,nlev,nv3d)
  REAL(r_sngl), ALLOCATABLE :: v2dg(:,:,:) !(nlon,nlat,nv2d)
  INTEGER :: l,n,ll,im
  CHARACTER(11) :: filename='file000.grd'
  INTEGER :: i,j,k,m !STEVE: for debugging
  LOGICAL :: verbose = .false.
  INTEGER :: convcnt = 0

  ALLOCATE(v3dg(nlon,nlat,nlev,nv3d),v2dg(nlon,nlat,nv2d))

  ll = CEILING(REAL(member)/REAL(nprocs))
  DO l=1,ll
    DO n=0,nprocs-1
      im = n+1 + (l-1)*nprocs
      IF(im <= member) THEN
        CALL gather_grd_mpi(n,v3d(:,:,im,:),v2d(:,im,:),v3dg,v2dg)
      END IF
    END DO

    im = myrank+1 + (l-1)*nprocs
    IF(im <= member) THEN
      WRITE(filename(1:7),'(A4,I3.3)') file,im
      WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing a file ',filename

      !STEVE: debug
      print *, "common_mpi_mom4.f90::write_ens_mpi:: MAXVAL(ABS(v3dg(:,:,:,iv3d_t))) = ", MAXVAL(ABS(v3dg(:,:,:,iv3d_t)))

      CALL write_bingrd4(filename,v3dg,v2dg)
    END IF

  END DO

  DEALLOCATE(v3dg,v2dg)

  RETURN
END SUBROUTINE write_ens_mpi_grd

!-----------------------------------------------------------------------
! gridded data -> buffer
!-----------------------------------------------------------------------
SUBROUTINE grd_to_buf(grd,buf)
  REAL(r_sngl),INTENT(IN) :: grd(nlon,nlat)
  REAL(r_sngl),INTENT(OUT) :: buf(nij1max,nprocs)
  INTEGER :: i,j,m,ilon,ilat

  DO m=1,nprocs
    DO i=1,nij1node(m)
      j = m-1 + nprocs * (i-1)
      ilon = MOD(j,nlon) + 1
      ilat = (j-ilon+1) / nlon + 1
      buf(i,m) = grd(ilon,ilat)
    END DO
  END DO

  RETURN
END SUBROUTINE grd_to_buf
!-----------------------------------------------------------------------
! buffer -> gridded data
!-----------------------------------------------------------------------
SUBROUTINE buf_to_grd(buf,grd)
  REAL(r_sngl),INTENT(IN) :: buf(nij1max,nprocs)
  REAL(r_sngl),INTENT(OUT) :: grd(nlon,nlat)
  INTEGER :: i,j,m,ilon,ilat

  DO m=1,nprocs
    DO i=1,nij1node(m)
      j = m-1 + nprocs * (i-1)
      ilon = MOD(j,nlon) + 1
      ilat = (j-ilon+1) / nlon + 1
      grd(ilon,ilat) = buf(i,m)
    END DO
  END DO

  RETURN
END SUBROUTINE buf_to_grd
!-----------------------------------------------------------------------
! STORING DATA (ensemble mean and spread)
!-----------------------------------------------------------------------
SUBROUTINE write_ensmspr_mpi(file,member,v3d,v2d)
  CHARACTER(4),INTENT(IN) :: file
  INTEGER,INTENT(IN) :: member
  REAL(r_size),INTENT(IN) :: v3d(nij1,nlev,member,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nij1,member,nv2d)
  REAL(r_size), ALLOCATABLE :: v3dm(:,:,:) !(nij1,nlev,nv3d)
  REAL(r_size), ALLOCATABLE :: v2dm(:,:) !(nij1,nv2d)
  REAL(r_size), ALLOCATABLE :: v3ds(:,:,:) !(nij1,nlev,nv3d)
  REAL(r_size), ALLOCATABLE :: v2ds(:,:) !(nij1,nv2d)
  REAL(r_sngl), ALLOCATABLE :: v3dg(:,:,:,:) !(nlon,nlat,nlev,nv3d)
  REAL(r_sngl), ALLOCATABLE :: v2dg(:,:,:) !(nlon,nlat,nv2d)
  INTEGER :: i,k,m,n,j,l,ll,im,mstart,mend
  CHARACTER(11) :: filename='file000.grd'

  ALLOCATE(v3dm(nij1,nlev,nv3d),v2dm(nij1,nv2d))
  ALLOCATE(v3ds(nij1,nlev,nv3d),v2ds(nij1,nv2d))
  ALLOCATE(v3dg(nlon,nlat,nlev,nv3d),v2dg(nlon,nlat,nv2d))

  CALL ensmean_grd(member,nij1,v3d,v2d,v3dm,v2dm)

! ll = CEILING(REAL(member)/REAL(nprocs))
! DO l=1,ll
!   mstart = 1 + (l-1)*nprocs
!   mend = MIN(l*nprocs, member)
!   CALL gather_grd_mpi_alltoall(mstart,mend,member,v3d,v2d,v3dg,v2dg)
! ENDDO
  CALL gather_grd_mpi(0,v3dm,v2dm,v3dg,v2dg)

  IF(myrank == 0) THEN
    WRITE(filename(1:7),'(A4,A3)') file,'_me'
    WRITE(6,'(A,I3.3,2A)') 'write_ensmspr_mpi::MYRANK ',myrank,' is writing a file ',filename
    CALL write_bingrd4(filename,v3dg,v2dg)
  END IF

  DO n=1,nv3d
    DO k=1,nlev
      DO i=1,nij1
        v3ds(i,k,n) = (v3d(i,k,1,n)-v3dm(i,k,n))**2
        DO m=2,member
          v3ds(i,k,n) = v3ds(i,k,n) + (v3d(i,k,m,n)-v3dm(i,k,n))**2
        END DO
        v3ds(i,k,n) = SQRT(v3ds(i,k,n) / REAL(member-1,r_size))
      END DO
    END DO
  END DO

  DO n=1,nv2d
    DO i=1,nij1
      v2ds(i,n) = (v2d(i,1,n)-v2dm(i,n))**2
      DO m=2,member
        v2ds(i,n) = v2ds(i,n) + (v2d(i,m,n)-v2dm(i,n))**2
      END DO
      v2ds(i,n) = SQRT(v2ds(i,n) / REAL(member-1,r_size))
    END DO
  END DO

! DO l=1,ll
!   mstart = 1 + (l-1)*nprocs
!   mend = MIN(l*nprocs, member)
!   CALL gather_grd_mpi_alltoall(mstart,mend,member,v3d,v2d,v3dg,v2dg)
! ENDDO
  CALL gather_grd_mpi(0,v3ds,v2ds,v3dg,v2dg)

  IF(myrank == 0) THEN
    WRITE(filename(1:7),'(A4,A3)') file,'_sp'
    WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing a file ',filename
    CALL write_bingrd4(filename,v3dg,v2dg)
  END IF

  DEALLOCATE(v3dm,v2dm,v3ds,v2ds,v3dg,v2dg)

  RETURN
END SUBROUTINE write_ensmspr_mpi

!-----------------------------------------------------------------------
! Scatter gridded data using MPI_ALLTOALL(V) (mstart~mend -> all)
!-----------------------------------------------------------------------

SUBROUTINE scatter_grd_mpi_alltoall(mstart,mend,member,v3dg,v2dg,v3d,v2d)
  INTEGER,INTENT(IN) :: mstart,mend,member
  REAL(r_sngl),INTENT(IN) :: v3dg(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(IN) :: v2dg(nlon,nlat,nv2d)
  REAL(r_size),INTENT(OUT) :: v3d(nij1,nlev,member,nv3d)
  REAL(r_size),INTENT(OUT) :: v2d(nij1,member,nv2d)
  REAL(r_sngl),ALLOCATABLE :: bufs3(:,:,:) , bufr3(:,:,:)
  REAL(r_sngl),ALLOCATABLE :: bufs2(:,:)   , bufr2(:,:)
  INTEGER :: k,n,m,mcount,ierr
  INTEGER :: ns(nprocs),nst(nprocs),nr(nprocs),nrt(nprocs)
  LOGICAL :: dodebug = .false.

  mcount = mend - mstart + 1
  IF(mcount > nprocs .OR. mcount <= 0) STOP

  if (dodebug) WRITE(6,*) "ALLOCATE bufs3,bufr3,bufs2,bufr2..."
  ALLOCATE( bufs3(nij1max,nlev,nprocs) , bufr3(nij1max,nlev,mcount) )
  ALLOCATE( bufs2(nij1max,nprocs)    , bufr2(nij1max,mcount) )

  if (dodebug) WRITE(6,*) "Cycling n=1,nv3d..."
  DO n=1,nv3d
    if (dodebug) WRITE(6,*) "n = ",n
    IF(myrank < mcount) THEN
      DO k=1,nlev
        if (dodebug) WRITE(6,*) "grd_to_buf :: k = ", k
        CALL grd_to_buf(v3dg(:,:,k,n),bufs3(:,k,:)) !,nlon,nlat,nij1max,nij1node)
      END DO
    END IF

    if (dodebug) WRITE(6,*) "==MPI_BARRIER=="
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    IF(mcount == nprocs) THEN
      if (dodebug) WRITE(6,*) "==MPI_ALLTOALL=="
      CALL MPI_ALLTOALL(bufs3, nij1max*nlev, MPI_REAL, &
                        bufr3, nij1max*nlev, MPI_REAL, MPI_COMM_WORLD, ierr)
    ELSE
      if (dodebug) WRITE(6,*) "set_alltoallv_counts..."
      CALL set_alltoallv_counts(mcount,nij1max*nlev,nr,nrt,ns,nst)
      if (dodebug) WRITE(6,*) "==MPI_ALLTOALLV=="
      CALL MPI_ALLTOALLV(bufs3, ns, nst, MPI_REAL, &
                         bufr3, nr, nrt, MPI_REAL, MPI_COMM_WORLD, ierr)
    END IF

    DO m = mstart,mend
      DO k=1,nlev
         if (dodebug) WRITE(6,*) "m,k,n = ",m,k,n
         if (dodebug) WRITE(6,*) "Assigning v3d(:,k,m,n) = REAL(bufr3(1:nij1,k,m-mstart+1),r_size)"
         if (dodebug) WRITE(6,*) "nij1,m-mstart+1 = ", nij1,m-mstart+1
         v3d(:,k,m,n) = REAL(bufr3(1:nij1,k,m-mstart+1),r_size)
      END  DO
    END DO
  END DO

  if (dodebug) WRITE(6,*) "Cycling n=1,nv2d..."
  DO n=1,nv2d
    IF(myrank < mcount) THEN
      if (dodebug) WRITE(6,*) "grd_to_buf :: 2D"
      CALL grd_to_buf(v2dg(:,:,n),bufs2(:,:)) !,nlon,nlat,nij1max,nij1node)
    END IF

    if (dodebug) WRITE(6,*) "==MPI_BARRIER=="
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    IF(mcount == nprocs) THEN
      if (dodebug) WRITE(6,*) "==MPI_ALLTOALL=="
      CALL MPI_ALLTOALL(bufs2, nij1max, MPI_REAL, &
                        bufr2, nij1max, MPI_REAL, MPI_COMM_WORLD, ierr)
    ELSE
      if (dodebug) WRITE(6,*) "set_alltoallv_counts..."
      CALL set_alltoallv_counts(mcount,nij1max,nr,nrt,ns,nst)
      if (dodebug) WRITE(6,*) "==MPI_ALLTOALLV=="
      CALL MPI_ALLTOALLV(bufs2, ns, nst, MPI_REAL, &
                         bufr2, nr, nrt, MPI_REAL, MPI_COMM_WORLD, ierr)
    END IF

    DO m = mstart,mend
      if (dodebug) WRITE(6,*) "m,n = ", m,n
      if (dodebug) WRITE(6,*) "Assigning v2d(:,m,n) = REAL(bufr2(1:nij1,m-mstart+1),r_size)"
      if (dodebug) WRITE(6,*) "nij1,m-mstart+1 = ", nij1,m-mstart+1
      v2d(:,m,n) = REAL(bufr2(1:nij1,m-mstart+1),r_size)
    END DO
  END DO

  if (dodebug) WRITE(6,*) "==MPI_BARRIER=="
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  DEALLOCATE(bufr3,bufs3,bufr2,bufs2)
  RETURN
END SUBROUTINE scatter_grd_mpi_alltoall

!-----------------------------------------------------------------------
! Gather gridded data using MPI_ALLTOALL(V) (all -> mstart~mend)
!-----------------------------------------------------------------------

SUBROUTINE gather_grd_mpi_alltoall(mstart,mend,member,v3d,v2d,v3dg,v2dg)
  INTEGER,INTENT(IN) :: mstart,mend,member
  REAL(r_size),INTENT(IN) :: v3d(nij1,nlev,member,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nij1,member,nv2d)
  REAL(r_sngl),INTENT(OUT) :: v3dg(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(OUT) :: v2dg(nlon,nlat,nv2d)
  REAL(r_sngl),ALLOCATABLE :: bufs3(:,:,:) , bufr3(:,:,:)
  REAL(r_sngl),ALLOCATABLE :: bufs2(:,:)   , bufr2(:,:)
  INTEGER :: k,n,m,mcount,ierr
  INTEGER :: ns(nprocs),nst(nprocs),nr(nprocs),nrt(nprocs)

  mcount = mend - mstart + 1
  IF(mcount > nprocs .OR. mcount <= 0) STOP

  ALLOCATE( bufs3(nij1max,nlev,mcount) , bufr3(nij1max,nlev,nprocs) )
  ALLOCATE( bufs2(nij1max,mcount)    , bufr2(nij1max,nprocs) )

  DO n=1,nv3d
    DO m = mstart,mend
      DO k=1,nlev
        bufs3(1:nij1,k,m-mstart+1) = REAL(v3d(:,k,m,n),r_sngl)
      END DO
    END DO

    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    IF(mcount == nprocs) THEN
      CALL MPI_ALLTOALL(bufs3, nij1max*nlev, MPI_REAL, &
                        bufr3, nij1max*nlev, MPI_REAL, MPI_COMM_WORLD, ierr)
    ELSE
      CALL set_alltoallv_counts(mcount,nij1max*nlev,ns,nst,nr,nrt)
      CALL MPI_ALLTOALLV(bufs3, ns, nst, MPI_REAL, &
                         bufr3, nr, nrt, MPI_REAL, MPI_COMM_WORLD, ierr)
    END IF

    IF(myrank < mcount) THEN
      DO k=1,nlev
        CALL buf_to_grd(bufr3(:,k,:),v3dg(:,:,k,n))
      ENDDO
    END IF
  END DO

  DO n=1,nv2d
    DO m = mstart,mend
      DO k=1,nlev
        bufs2(1:nij1,m-mstart+1) = REAL(v2d(:,m,n),r_sngl)
      END DO
    END DO

    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    IF(mcount == nprocs) THEN
      CALL MPI_ALLTOALL(bufs2, nij1max, MPI_REAL, &
                        bufr2, nij1max, MPI_REAL, MPI_COMM_WORLD, ierr)
    ELSE
      CALL set_alltoallv_counts(mcount,nij1max,ns,nst,nr,nrt)
      CALL MPI_ALLTOALLV(bufs2, ns, nst, MPI_REAL, &
                         bufr2, nr, nrt, MPI_REAL, MPI_COMM_WORLD, ierr)
    END IF

    IF(myrank < mcount) THEN
      CALL buf_to_grd(bufr2(:,:),v2dg(:,:,n))
    END IF
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  DEALLOCATE(bufr3,bufs3,bufr2,bufs2)
  RETURN
END SUBROUTINE gather_grd_mpi_alltoall

!-----------------------------------------------------------------------
! Set the send/recieve counts of MPI_ALLTOALLV
!-----------------------------------------------------------------------
SUBROUTINE set_alltoallv_counts(mcount,ngpblock,n_ens,nt_ens,n_mem,nt_mem)
  INTEGER,INTENT(IN) :: mcount,ngpblock
  INTEGER,INTENT(OUT) :: n_ens(nprocs),nt_ens(nprocs),n_mem(nprocs),nt_mem(nprocs)
  INTEGER :: p

  n_ens = 0
  nt_ens = 0
  n_mem = 0
  nt_mem = 0
  DO p=1,mcount
    n_ens(p) = ngpblock
    IF(myrank+1 == p) THEN
      n_mem(:) = ngpblock
    END IF
  END DO
  DO p=2,nprocs
    nt_ens(p) = nt_ens(p-1) + n_ens(p-1)
    nt_mem(p) = nt_mem(p-1) + n_mem(p-1)
  END DO

  RETURN
END SUBROUTINE set_alltoallv_counts


!-----------------------------------------------------------------------
! Scatter a smaller sized gridded data using MPI_ALLTOALL(V) (mstart~mend -> all)
!-----------------------------------------------------------------------

SUBROUTINE scatter_grd_mpi_smalltoall(mstart,mend,member,v2dg,v2d,nx,ny,nv)
  INTEGER,INTENT(IN) :: mstart,mend,member
  INTEGER,INTENT(IN) :: nx,ny,nv
  REAL(r_sngl),INTENT(IN) :: v2dg(nx,ny,nv)
  REAL(r_size),INTENT(OUT) :: v2d(nij1,member,nv)
  REAL(r_sngl),ALLOCATABLE :: bufs2(:,:)   , bufr2(:,:)
  INTEGER :: k,n,m,mcount,ierr
  INTEGER :: ns(nprocs),nst(nprocs),nr(nprocs),nrt(nprocs)
  LOGICAL :: dodebug = .true.

  mcount = mend - mstart + 1
  IF(mcount > nprocs .OR. mcount <= 0) STOP

  if (dodebug) WRITE(6,*) "ALLOCATE bufs2,bufr2..."
  ALLOCATE( bufs2(nij1max,nprocs)    , bufr2(nij1max,mcount) )

  if (dodebug) WRITE(6,*) "Cycling n=1,nv..."
  DO n=1,nv
    IF(myrank < mcount) THEN
      if (dodebug) WRITE(6,*) "grd_to_buf :: 2D"
      CALL grd_to_buf(v2dg(:,:,n),bufs2(:,:)) !,nlon,nlat,nij1max,nij1node)
    END IF

    if (dodebug) WRITE(6,*) "==MPI_BARRIER=="
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    IF(mcount == nprocs) THEN
      if (dodebug) WRITE(6,*) "==MPI_ALLTOALL=="
      CALL MPI_ALLTOALL(bufs2, nij1max, MPI_REAL, &
                        bufr2, nij1max, MPI_REAL, MPI_COMM_WORLD, ierr)
    ELSE
      if (dodebug) WRITE(6,*) "set_alltoallv_counts..."
      CALL set_alltoallv_counts(mcount,nij1max,nr,nrt,ns,nst)
      if (dodebug) WRITE(6,*) "==MPI_ALLTOALLV=="
      CALL MPI_ALLTOALLV(bufs2, ns, nst, MPI_REAL, &
                         bufr2, nr, nrt, MPI_REAL, MPI_COMM_WORLD, ierr)
    END IF

    DO m = mstart,mend
      if (dodebug) WRITE(6,*) "m,n = ", m,n
      if (dodebug) WRITE(6,*) "Assigning v2d(:,m,n) = REAL(bufr2(1:nij1,m-mstart+1),r_size)"
      if (dodebug) WRITE(6,*) "nij1,m-mstart+1 = ", nij1,m-mstart+1
      v2d(:,m,n) = REAL(bufr2(1:nij1,m-mstart+1),r_size)
    END DO
  END DO

  if (dodebug) WRITE(6,*) "==MPI_BARRIER=="
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  DEALLOCATE(bufr2,bufs2)
  RETURN
END SUBROUTINE scatter_grd_mpi_smalltoall

SUBROUTINE scatter_grd_mpi_small(nrank,v2dg,v2d,nx,ny,nv)
  INTEGER,INTENT(IN) :: nrank
  INTEGER,INTENT(IN) :: nx,ny,nv
  REAL(r_sngl),INTENT(IN) :: v2dg(nx,ny,nv)
  REAL(r_size),INTENT(OUT) :: v2d(nij1,nv)
  REAL(r_sngl), ALLOCATABLE :: bufs(:,:,:) !(nij1max,nlevall,nprocs)
  REAL(r_sngl), ALLOCATABLE :: bufr(:,:) !(nij1max,nlevall)
  INTEGER :: j,k,n,ierr,ns,nr

  ALLOCATE(bufs(nij1max,nlevall,nprocs), bufr(nij1max,nlevall))

  ns = nij1max * nlevall
  nr = ns
  IF(myrank == nrank) THEN
    j=0
    DO n=1,nv
      j = j+1
      CALL grd_to_buf(v2dg(:,:,n),bufs(:,j,:))
    END DO
  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_SCATTER(bufs,ns,MPI_REAL,&
                 & bufr,nr,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)
  j=0
  DO n=1,nv
    j = j+1
    v2d(:,n) = REAL(bufr(1:nij1,j),r_size)
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  DEALLOCATE(bufs,bufr)

  RETURN
END SUBROUTINE scatter_grd_mpi_small

END MODULE common_mpi_mom4
