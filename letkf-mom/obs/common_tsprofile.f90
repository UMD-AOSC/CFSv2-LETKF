module common_tsprofile
  use netcdf
  use common, only : r_size, r_sngl
  use common_obs, only : obsid_ocn_t, obsid_ocn_s
  use vars_model, only : lon0, lonf, lat0, latf, wrapgap
  implicit none

  private

  public :: tsprofile_getnum
  public :: read_tsprofile
  public :: write_tsprofile_diag

  logical :: verbose = .TRUE.
  logical :: debug = .FALSE.
  !INTEGER,PARAMETER :: r_size=kind(0.0d0)
  !INTEGER,PARAMETER :: r_dble=kind(0.0d0)
  !INTEGER,PARAMETER :: r_sngl=kind(0.0e0)

  !integer, parameter :: obsid_ocn_t    = 2210  ! hardcoded (should be removed later)
  !integer, parameter :: obsid_ocn_s    = 2220  ! hardcoded

contains 

subroutine check(status)
  integer, intent(in) :: status
  if (status/=nf90_noerr) then
     print *, status
     print *, trim(nf90_strerror(status))
     stop 1
  end if
end subroutine check

function tsprofile_getnum(infile) result(nobs)
  implicit none

  character(*),intent(in) :: infile
  integer :: nobs

  integer :: ncid, did_p, did_lev
  integer :: np, nlev

  nobs = 0
  call check(nf90_open(trim(infile), nf90_nowrite, ncid))
  call check(nf90_inq_dimid(ncid, "count", did_p))
  call check(nf90_inquire_dimension(ncid, did_p, len=np))
  call check(nf90_inq_dimid(ncid, "grid_z", did_lev))
  call check(nf90_inquire_dimension(ncid, did_lev, len=nlev))
  call check(nf90_close(ncid))

  nobs = np * nlev
  if (verbose) write(6,*) "[msg] tsprofile_getnum: np, nz, nobs=", np, nlev, nobs

endfunction tsprofile_getnum


subroutine read_tsprofile(infile,nobs,elem,rlon,rlat,rlev,odat)
  implicit none

  character(*),intent(in) :: infile
  integer,     intent(in) :: nobs
  real(r_size),intent(out) :: elem(nobs)
  real(r_size),intent(out) :: rlon(nobs)
  real(r_size),intent(out) :: rlat(nobs)
  real(r_size),intent(out) :: rlev(nobs)
  real(r_size),intent(out) :: odat(nobs)

  integer :: ncid, did_lev, did_p, vid_lon, vid_lat, vid_lev, vid_var
  integer :: np, nlev
  integer :: ierr_sal,ierr_tmp
  integer :: i, k, n, n1, n2, ierr

  real(r_sngl),allocatable :: lon_p(:),lat_p(:) ! np
  real(r_sngl),allocatable :: lev_common(:) ! nlev
  real(r_sngl),allocatable :: var(:,:) ! nlev*np
  integer :: elem_common 


  call check(nf90_open(trim(infile), nf90_nowrite, ncid))
  call check(nf90_inq_dimid(ncid, "count", did_p))
  call check(nf90_inquire_dimension(ncid, did_p, len=np))
  call check(nf90_inq_dimid(ncid, "grid_z", did_lev))
  call check(nf90_inquire_dimension(ncid, did_lev, len=nlev))

  allocate(lon_p(np),lat_p(np))
  allocate(lev_common(nlev))
  allocate(var(nlev,np))

  call check(nf90_inq_varid(ncid,"grid_z",vid_lev))
  call check(nf90_get_var(ncid,vid_lev,lev_common))

  call check(nf90_inq_varid(ncid,"xlon",vid_lon))
  call check(nf90_get_var(ncid,vid_lon,lon_p))

  call check(nf90_inq_varid(ncid,"ylat",vid_lat))
  call check(nf90_get_var(ncid,vid_lat,lat_p))

  ierr = nf90_inq_varid(ncid, "salt", vid_var)
  if (ierr /= nf90_noerr) then
    if(verbose) write(6,*) "[msg] read_tsprofile: temp profile"
    ierr = nf90_inq_varid(ncid, "temp", vid_var)
    call check(nf90_get_var(ncid,vid_var,var))
    elem_common = obsid_ocn_t
  else
    if(verbose) write(6,*) "[msg] read_tsprofile: salt profile"
    call check(nf90_get_var(ncid,vid_var,var))
    elem_common = obsid_ocn_s
  endif
  call check(nf90_close(ncid))

  write(6,*) "*********lon0, lonf, wrapgap=", lon0, lonf, wrapgap

! convert 2d data to 1d
  do i = 1, np
     do k = 1, nlev
        n = nlev*(i-1) + k
        elem(n) = elem_common
        rlon(n) = lon_p(i)
        rlat(n) = lat_p(i)
        rlev(n) = lev_common(k)
        odat(n) = var(k,i) ! nlev*np

    !   change longitude to the range required by the model
        if (.true.) then
          IF(rlon(n) >= lonf) THEN
            ! Update the coordinate
            if (abs(rlon(n) - lonf) < wrapgap ) then
              ! First, handle observations that are just outside of the model grid
              !STEVE: shift it if it's just outside grid
              if (abs(rlon(n) - lonf) < wrapgap/2) then
                rlon(n) = lonf
              else
                rlon(n) = lon0
              endif
            else
              ! Otherwise, wrap the observation coordinate to be inside of the defined model grid coordinates
              !Wrap the coordinate
              rlon(n) = REAL(lon0 + abs(rlon(n) - lonf) - wrapgap,r_size)
            endif

          ELSE IF(rlon(n) < lon0) THEN

            ! Update the coordinate
            if (abs(lon0 - rlon(n)) < wrapgap ) then
              !STEVE: shift it if it's just outside grid
              if (abs(lon0 - rlon(n)) < wrapgap/2) then
                rlon(n) = lon0
              else
                rlon(n) = lonf 
              endif
            else
              !Wrap the coordinate
              rlon(n) = REAL(lonf - abs(lon0 - rlon(n)) + wrapgap,r_size)
            endif

          ENDIF

        endif 


        if (verbose.and.debug) then
           write(100,"(999(F20.13,1X))") elem(n), rlon(n), rlat(n), rlev(n), odat(n)
        endif

     enddo


     if (verbose.and.debug) then
        n1 = nlev*(i-1) + 1
        n2 = nlev*i
        write(6,*) "P elem: min, max=", minval(elem(n1:n2)), maxval(elem(n1:n2))
        write(6,*) "P rlon: min, max=", minval(rlon(n1:n2)), maxval(rlon(n1:n2))
        write(6,*) "P rlat: min, max=", minval(rlat(n1:n2)), maxval(rlat(n1:n2))
        write(6,*) "P rlev: min, max=", minval(rlev(n1:n2)), maxval(rlev(n1:n2))
        write(6,*) "P odat: min, max=", minval(odat(n1:n2)), maxval(odat(n1:n2))
     endif
  enddo

  if(verbose) then
    write(6,*) "[msg] read_tsprofile: TOTAL elem: min, max=", minval(elem), maxval(elem)
    write(6,*) "[msg] read_tsprofile: TOTAL rlon: min, max=", minval(rlon), maxval(rlon)
    write(6,*) "[msg] read_tsprofile: TOTAL rlat: min, max=", minval(rlat), maxval(rlat)
    write(6,*) "[msg] read_tsprofile: TOTAL rlev: min, max=", minval(rlev), maxval(rlev)
    write(6,*) "[msg] read_tsprofile: TOTAL odat: min, max=", minval(odat), maxval(odat)
  endif

  deallocate(lon_p, lat_p, lev_common, var)

endsubroutine read_tsprofile


subroutine write_tsprofile_diag(outfile,nobs,elem,rlon,rlat,rlev,odat,ohx,oqc)
  implicit none

  character(*),intent(in) :: outfile
  integer, intent(in) :: nobs
  real(r_size), intent(in) :: elem(nobs), rlon(nobs), rlat(nobs), rlev(nobs), &
                              odat(nobs), ohx(nobs)
  integer, intent(in) :: oqc(nobs)
 
  integer :: n, lout

  lout = 20
  open(lout,file=trim(outfile),action="write")
  do n = 1, nobs
     write(lout,"(100(F20.13,1x))") elem(n), rlon(n), rlat(n), rlev(n), odat(n), ohx(n), real(oqc(n),r_size)
  enddo
  close(lout)  

endsubroutine write_tsprofile_diag


subroutine examples()
  implicit none

  character(80) :: infile, outfile
  integer :: nobs

  REAL(r_size), ALLOCATABLE :: elem(:)
  REAL(r_size), ALLOCATABLE :: rlon(:)
  REAL(r_size), ALLOCATABLE :: rlat(:)
  REAL(r_size), ALLOCATABLE :: rlev(:)
  REAL(r_size), ALLOCATABLE :: odat(:)
  REAL(r_size), ALLOCATABLE :: ohx(:)
  INTEGER     , ALLOCATABLE :: oqc(:)


  call get_command_argument(1,infile)
  outfile="ascii."//trim(infile)
  print*, "infile=", trim(infile)

! get obs num
  nobs = tsprofile_getnum(trim(infile))
  print*, "nobs=", nobs 

  ALLOCATE( elem(nobs) )
  ALLOCATE( rlon(nobs) )
  ALLOCATE( rlat(nobs) )
  ALLOCATE( rlev(nobs) )
  ALLOCATE( odat(nobs) )
  ALLOCATE( ohx(nobs) )
  ALLOCATE( oqc(nobs) )

! read in obs
  CALL read_tsprofile(trim(infile),nobs,elem,rlon,rlat,rlev,odat)

! write out
  call write_tsprofile_diag(trim(outfile),nobs,elem,rlon,rlat,rlev,odat,ohx,oqc)

endsubroutine examples



endmodule common_tsprofile


!program main
!  use common_tsprofile
!  implicit none
!  call examples()
!endprogram
