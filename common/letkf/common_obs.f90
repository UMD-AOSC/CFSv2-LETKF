!------------------------------------------------------------
! CFS common observations
! This contains the definitions for observations used by the
! CFSv2-LETKF. For strongly coupled data assimilation, observation
! ID's are required to be unique between the ocean and atmosphere
!------------------------------------------------------------
module common_obs

  use common

  implicit none  
  public

  !!------------------------------------------------------------
  !!------------------------------------------------------------
  !! Observations
  !!------------------------------------------------------------
  !!------------------------------------------------------------

  integer, parameter :: dom_atm = 1000
  integer, parameter :: dom_ocn = 2000
  
  !! unique ID's for observations
  !!------------------------------
  integer, parameter :: obsid_num = 16

  !! atmosphere obs
  integer, parameter :: obsid_atm_min  = 1000
  integer, parameter :: obsid_atm_max  = 1999
  integer, parameter :: obsid_atm_num  = 8
  integer, parameter :: obsid_atm_offset = 0 
  
  integer, parameter :: obsid_atm_ps   = 1100
  integer, parameter :: obsid_atm_rain = 1110
  integer, parameter :: obsid_atm_t    = 1210
  integer, parameter :: obsid_atm_tv   = 1211
  integer, parameter :: obsid_atm_q    = 1220
  integer, parameter :: obsid_atm_rh   = 1221
  integer, parameter :: obsid_atm_u    = 1250
  integer, parameter :: obsid_atm_v    = 1251

  !! ocean obs
  integer, parameter :: obsid_ocn_min  = 2000
  integer, parameter :: obsid_ocn_max  = 2999
  integer, parameter :: obsid_ocn_num  = 8
  integer, parameter :: obsid_ocn_offset = 8

  integer, parameter :: obsid_ocn_ssh  = 2100  ! (ADT)
  integer, parameter :: obsid_ocn_eta  = 2101  ! (SLA)
  integer, parameter :: obsid_ocn_sst  = 2110
  integer, parameter :: obsid_ocn_sss  = 2120
  integer, parameter :: obsid_ocn_t    = 2210
  integer, parameter :: obsid_ocn_s    = 2220
  integer, parameter :: obsid_ocn_u    = 2250
  integer, parameter :: obsid_ocn_v    = 2251

  !! arrays holding all observation id's and names, for easy iteration
  !! in loops that want to print stats for obs
  integer :: obsid_ids(obsid_num) = (/&
       obsid_atm_ps, obsid_atm_rain, obsid_atm_t, obsid_atm_tv, &
       obsid_atm_q, obsid_atm_rh, obsid_atm_u, obsid_atm_v, &
       obsid_ocn_ssh, obsid_ocn_eta, obsid_ocn_sst, obsid_ocn_sss, &
       obsid_ocn_t, obsid_ocn_s, obsid_ocn_u, obsid_ocn_v/)
  character (len=10) :: obsid_names(obsid_num) = (/&
       "ATM_PS  ", "ATM_RAIN", "ATM_T   ", "ATM_TV  ", &
       "ATM_Q   ", "ATM_RH  ", "ATM_U   ", "ATM_V   ", &
       "OCN_SSH ", "OCN_ETA ", "OCN_SST ", "OCN_SSS ",&
       "OCN_T   ", "OCN_S   ", "OCN_U   ", "OCN_V   "/)


    
  !!
  !! ------------------------------------------------------------
  !! Structure to hold LETKF formated observations,
  !! see the wiki for further documentation. This is the
  !! common format produced by the observation operators
  !! ------------------------------------------------------------
  type :: Observation
     integer      :: id     !! one of the obsid_* values from above
     real(r_size) :: lon    !! longitude (degrees)
     real(r_size) :: lat    !! latitude (degrees)
     real(r_size) :: lev    !! height/depth (mb, meters)
     real(r_size) :: odat   !! observation value
     real(r_size) :: oerr   !! estimate observation error
     integer      :: platform !! platform ID (currently not used)
     
     real(r_size) :: time   !! time realtive to analysis
     real(r_size) :: ohx    !! model value in observation space
     integer      :: qc     !! quality control (1 = valid, 0 = invalid)
  end type Observation


  !!------------------------------------------------------------
  !!------------------------------------------------------------
  !! Platform types, used mainly for QC, atmosphere only now
  !!------------------------------------------------------------
  !!------------------------------------------------------------
  integer, parameter :: platform_num = 21
  character (6), parameter :: platform_name(platform_num)= (/&
       'ADPUPA', 'AIRCAR', 'AIRCFT', 'SATWND', 'PROFLR', &
       'VADWND', 'SATEMP', 'ADPSFC', 'SFCSHP', 'SFCBOG', &
       'SPSSMI', 'SYNDAT', 'ERS1DA', 'GOESND', 'QKSWND', &
       'MSONET', 'GPSIPW', 'RASSDA', 'WDSATR', 'ASCATW', &
       'TMPAPR'/)

  
contains

  
  function getDomain(obsid)
    integer, intent(in) :: obsid
    integer ::  getDomain

    if (obsid >= obsid_atm_min .and. obsid <= obsid_atm_max) then
       getDomain = dom_atm
    else if(obsid >= obsid_ocn_min .and. obsid <= obsid_ocn_max) then
       getDomain = dom_ocn
    else
       getDomain = -1
    end if
  end function getDomain

  
  !! ------------------------------------------------------------
  !! given an integer observation ID, returns its location in the
  !! above obsid_ids and obdsid_name arrays (and presumably
  !! from any user specified arrays involving obsids as well)
  !! returns -1 if not found
  !! ------------------------------------------------------------
  function obsid_id2idx(id)
    integer, intent(in) ::id
    integer ::  obsid_id2idx
    integer n

    obsid_id2idx = -1
    do n=1,obsid_num
       if (obsid_ids(n) .eq. id) then
          obsid_id2idx = n
          return
       end if
    end do
  end function obsid_id2idx



  !! ------------------------------------------------------------
  !! ------------------------------------------------------------
  function obs_getnum(file, extended) result(nobs)
    character(*), intent(in) :: file
    logical, intent(in), optional :: extended
    logical :: ext
    integer :: nobs, unit, ios, recl, n
    real(r_sngl),allocatable :: wk(:)
    integer :: obs_count(obsid_num),inv_count
    
    !! determine if short, or extended file format
    ext = .false.
    if (present(extended)) ext = extended
    if (ext) then
       recl = 10
    else
       recl = 7
    end if
    allocate(wk(recl))

    !! open the file and count the number of observations
    nobs = 0
    obs_count(:) = 0
    inv_count = 0
    open(newunit=unit, file=file, form='unformatted', access='sequential')
    do
       read(unit, iostat=ios) wk
       if (ios /= 0) exit  !! end of the file found, exit loop
       !! another record found, increment count
       if ( nint(wk(1)) .gt. 0) then          
          nobs = nobs + 1
          n = obsid_id2idx( nint(wk(1)) )
          if (n < 0 ) then
             inv_count = inv_count + 1
          else
             obs_count(n) = obs_count(n) +1
          end if
          
       else
          inv_count = inv_count + 1
       end if
       
    end do

    !!done, cleanup
    deallocate(wk)
    close(unit)

    !! print statistics
    write (6,*) ""    
    write (6,*) "--------------------------"
    write (6,*) "Observation count check"
    write (6,*) "File: ",file
    write (6,*) "--------------------------"
    do n=1,obsid_num
       if (obs_count(n) > 0) then
          write (6,*) obsid_names(n), obs_count(n)
       end if
    end do
    write (6,*) "--------------------------"    
    write (6,*) "Total Valid: ",nobs
    write (6,*) "Invalid:     ",inv_count
    write (6,*) ""
    
  end function obs_getnum


  
  !! ------------------------------------------------------------
  !! reads in an LETKF formatted observation file.
  !!  This file contains a series of observations each with
  !!  7 values, or 10 values (if an extended format)
  !!  The short format is what is read into the obsop programs.
  !!  The extended format is read into/out of the LETKF program
  !!  Memory for "obs" will be allocated with this subroutine.
  !!  When the user is finished with obs, they should deallocate it.
  !! ------------------------------------------------------------
  subroutine obs_read(file, obs, extended)
    character(*), intent(in) :: file
    type(Observation), allocatable, intent(out) :: obs(:)
    logical, intent(in), optional :: extended
    real(r_sngl),allocatable :: wk(:)
    logical :: ex, ext
    integer :: count, unit, ios, n, recl

    !! make sure the file exists
    inquire(file=file, exist=ex)
    if ( .not. ex) then
       write (6,*) file, ' does not exist -- skipped'
    end if

    !! determine if short, or extended file format
    ext = .false.
    if (present(extended)) ext = extended
    if (ext) then
       recl = 10
    else
       recl = 7
    end if
    allocate(wk(recl))

    !! count the number of records, and allocate space
    count = obs_getnum(file, ext)
    write (6,*) 'Reading in ', count,' observations'
    allocate( obs(count) )

    !! read in the observations
    open(newunit=unit, file=file, form='unformatted', access='sequential')
    do n=1,count
       read(unit) wk
       obs(n)%id        = nint(wk(1))
       obs(n)%lon       = wk(2)
       obs(n)%lat       = wk(3)
       obs(n)%lev       = wk(4)
       obs(n)%odat      = wk(5)
       obs(n)%oerr      = wk(6)
       obs(n)%platform  = wk(7)
       if (ext) then
          obs(n)%time   = wk(8)
          obs(n)%ohx    = wk(9)
          obs(n)%qc     = wk(10)
       else
          obs(n)%time   = 0
          obs(n)%ohx    = 0
          obs(n)%qc     = 0
       end if
       if (obsid_id2idx(obs(n)%id) < 0) then
          write (6,*) "ERROR: Invalid observation ID", obs(n)%id
          stop 1
       end if
    end do

    write (*,*) size(obs)

    !! all done, close up
    deallocate(wk)
    close(unit)
  end subroutine obs_read


  !! ------------------------------------------------------------
  !! Writes an array of observations into a file with the LETKF format
  !!  See "read_obs" for further comments.
  !! ------------------------------------------------------------
  subroutine obs_write(file, obs, extended)
    character(*), intent(in) :: file
    type(Observation), intent(in) :: obs(:)
    logical, intent(in), optional :: extended
    real(r_sngl),allocatable :: wk(:)
    logical :: ex, ext
    integer :: recl, n, unit

    !! see if file already exists
    inquire(file=file, exist=ex)
    if (ex) then
       write(6,*) "WARNING: ",file," is being overwritten"
    end if

    !! determine if short, or extended file format
    ext = .false.
    if (present(extended)) ext = extended
    if (ext) then
       recl = 10
    else
       recl = 7
    end if
    allocate(wk(recl))

    !! open file and start writing records to it
    open(newunit=unit, file=file, form='unformatted', access='sequential')
    do n=1,size(obs)
       wk(1)  = obs(n)%id
       wk(2)  = obs(n)%lon
       wk(3)  = obs(n)%lat
       wk(4)  = obs(n)%lev
       wk(5)  = obs(n)%odat
       wk(6)  = obs(n)%oerr
       wk(7)  = obs(n)%platform
       if (ext) then
          wk(8)  = obs(n)%time
          wk(9)  = obs(n)%ohx
          wk(10) = obs(n)%qc
       end if
       write(unit) wk
    end do
    close(unit)
    
  end subroutine obs_write

  
  !!------------------------------------------------------------
  !! converts from an array of Observations types to several
  !! arrays of each variables. mainly for dealing with legacy code.
  !!------------------------------------------------------------
  !! TODO, try to get rid of references to this in the code
  subroutine obs_new2old(obs, elem,rlon,rlat,rlev,odat,oerr,otyp)
    type(Observation),intent(in) :: obs(:)
    real(r_size), intent(inout) :: elem(:)
    real(r_size), intent(inout) :: rlon(:)
    real(r_size), intent(inout) :: rlat(:)
    real(r_size), intent(inout) :: rlev(:)
    real(r_size), intent(inout) :: odat(:)
    real(r_size), intent(inout) :: oerr(:)
    real(r_size), intent(inout) :: otyp(:)

    integer n

    do n=1,size(obs)
       elem(n) = obs(n)%id
       rlon(n) = obs(n)%lon
       rlat(n) = obs(n)%lat
       rlev(n) = obs(n)%lev
       odat(n) = obs(n)%odat
       oerr(n) = obs(n)%oerr
       otyp(n) = obs(n)%platform

       !! convert pressure for atmosphere to Pa
       select case(obs(n)%id)
       case(obsid_atm_u, obsid_atm_v, obsid_atm_t,&
            obsid_atm_tv,obsid_atm_q)
          rlev(n) = rlev(n) * 100.0
       case ( obsid_atm_rh)
          rlev(n) = rlev(n) * 100.0
          odat(n) = odat(n) * 0.01
          oerr(n) = oerr(n) * 0.01
       case(obsid_atm_ps)
          odat(n) = odat(n) * 100.0
          oerr(n) = oerr(n) * 100.0
       end select       
    end do   
  end subroutine obs_new2old
  
end module common_obs
