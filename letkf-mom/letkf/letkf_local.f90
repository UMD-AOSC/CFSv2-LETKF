MODULE letkf_local
  USE common
  use common_obs
  use kdtree
  USE letkf_obs !contains debug_hdxf_0, and nobsgrd
  USE letkf_mom_params
  !  USE vars_letkf,   ONLY: var_local_n2n
  !  USE common_mpi
  !  USE common_mom4
  !  USE common_mpi_mom4
  !  USE common_letkf

  implicit none

  PUBLIC :: obs_local
  
  PRIVATE

  integer :: initialized = 0
  type(KD_ROOT) :: kdtree_root

  !!TODO, document
  integer, allocatable :: idx(:)
  real(r_size), allocatable :: dist(:)
  real(r_size) :: prev_lon, prev_lat
  integer :: nn

CONTAINS

  
SUBROUTINE obs_local(ij,ilev,var_local,hdxf,rdiag,rloc,dep,nobsl,nobstotal, useatm)
!===============================================================================
! Project global observations to local
!     (hdxf_global,dep_global,rdiag_global) -> (hdxf,dep,rdiag)
!===============================================================================
    integer,      intent(in)  :: ij,ilev,nobstotal
    real(r_size), intent(in)  :: var_local(obsid_ocn_num)
    real(r_size), intent(out) :: hdxf(nobstotal,nbv)
    real(r_size), INTENT(out) :: rdiag(nobstotal)
    real(r_size), INTENT(out) :: rloc(nobstotal)
    real(r_size), INTENT(out) :: dep(nobstotal)
    integer,      INTENT(out) :: nobsl
    logical, intent(in) :: useatm

    real(r_size),parameter  :: loc_cutoff = 4*10/3

    real(r_size) :: sigma_max_h_d0, sigma_h, sigma_atm_h_ij, sigma_ocn_h_ij
    real(r_size) :: loc, loc_h, loc_atm_v, loc_ocn_v
    real(r_size) :: dlev    
    integer :: n, j, k
    logical :: b

    !! determine the maximum horizontal search radius for the initial obs search 
    sigma_ocn_h_ij = (1.0d0-abs(lat1(ij))/90.0d0)*(sigma_ocn_h(1)-sigma_ocn_h(2))+sigma_ocn_h(2)
    sigma_atm_h_ij = (1.0d0-abs(lat1(ij))/90.0d0)*(sigma_atm_h(1)-sigma_atm_h(2))+sigma_atm_h(2)    
    sigma_max_h_d0 = max(sigma_ocn_h_ij, sigma_atm_h_ij)  * SQRT(loc_cutoff)

    !!  ------------------------------------------------------------
    !! Initialize the KD search tree
    if (initialized == 0) then
       WRITE(6,*) "Initializing the obs_local()"
       initialized = 1
       call kd_init(kdtree_root, obslon, obslat)
       WRITE(6,*) "Done constructing KD search tree"
       write (6,*) "nobstotal", nobstotal
       allocate(dist(nobstotal))
       allocate(idx(nobstotal))
       prev_lon = -1e10
       prev_lat = -1e10
    end if

    
    !! ------------------------------------------------------------
    !! query the KD tree
    if (lon1(ij) /= prev_lon .or. lat1(ij) /= prev_lat) then
       call kd_search(kdtree_root, obslon, obslat, (/lon1(ij),lat1(ij)/),&
            sigma_max_h_d0, idx, dist, nn)
       prev_lon = lon1(ij)
       prev_lat = lat1(ij)
    end if


    !! ------------------------------------------------------------
    !! for each observation found in the radius, do the localization and remove obs that aren't
    !! going to be used
    nobsl=0        
    do n = 1, nn
       loc_atm_v = 0.0
       loc_ocn_v = 0.0
       loc_h = 0.0
       
       if(nobsl >= nobstotal) return

       j = obselm(idx(n))

       !! calculate domain specific parameters
       !! ------------------------------       
       if (useatm .and. atm_obs .and. getDomain(j) == dom_atm) then
          !!------------------------------
          !! atmospheric observation
          !!------------------------------
          !! Has the user enabled assimilation of atmospheric observations?


          !! Has the user explicitly defined a subset of atmospheric platforms to use?
          !! TODO: this should be moved to early in the LETKF program
          b = .false.
          if (atm_obs_plat(1) > 0) then
             do k=1,size(atm_obs_plat)
                if (atm_obs_plat(k) <=0 ) exit
                if (atm_obs_plat(k) == obspla(idx(n))) then
                   b = .true.
                   exit
                end if
             end do
             if (.not. b) cycle
          end if
          
          !! has user set a subset of obs types to use?
          !! TODO: this as well should be moved to earlier in the LETKF
          !! to prevent the obs from even being placed in the kdtree
          b = .false.
          if (atm_obs_type(1) > 0) then
             do k=1,size(atm_obs_type)
                if (atm_obs_type(k) <= 0) exit
                if (atm_obs_type(k) == obselm(idx(n))) then
                   b = .true.
                   exit
                end if
             end do
             if (.not. b) cycle
          end if
          

          !! vertical localization into the atmosphere
          if (j /= obsid_atm_ps) then
             !!TODO, this assumes an atmospheric pressure of 1013mb over the ocean,
             !! for better accuracy the actual ensemble mean atmospheric pressure
             !! should be read in and used, though I don't know how much of an actual
             !! difference this would make, I'm lazy right now.
             loc_atm_v =  ((abs(log(obslev(idx(n)))-log(1013.0d0))) / sigma_atm_v) ** 2

             !! domain localization
             loc_atm_v = loc_atm_v + atm_loc
             
             if (loc_atm_v > loc_cutoff) cycle             
          end if

          !! vertical localization into the ocean
          dlev = lev(ilev)

          !! horizontal localization
          sigma_h = sigma_atm_h_ij          

       !! ------------------------------          
       else if (ocn_obs .and. getDomain(j) == dom_ocn) then
          !! ------------------------------          
          !! ocean observation
          !!------------------------------
          !! vertical localization within ocean
          dlev = abs(lev(ilev)-obslev(idx(n)))

          !! horizontal localization
          sigma_h = sigma_ocn_h_ij

       !! ------------------------------
       else
          !! ------------------------------
          !! unknown Domain
          !! ------------------------------
          cycle
       end if
       

       !! vertical localization into the ocean
       if (sigma_ocn_v > 0) then
          loc_ocn_v = (dlev / sigma_ocn_v) ** 2
       end if
             

       !! horizontal localization cutoff
       loc_h = (dist(n)/sigma_h)**2       

       !! total localization
       loc = loc_h + loc_atm_v + loc_ocn_v
       if (loc > loc_cutoff) cycle
       
       !! use this observation!
       nobsl = nobsl+1              
       rloc(nobsl) = exp(-0.5d0 * loc)
       hdxf(nobsl,:) = obshdxf(idx(n),:)
       dep(nobsl)  = obsdep(idx(n))
       rdiag(nobsl) = obserr(idx(n))**2
    end do

END SUBROUTINE obs_local




!(OCEAN) STEVE: add checks for atlantic/pacific basin boundary
PURE SUBROUTINE atlpac (xlat, xlon_in, lxap)
  REAL(r_size), INTENT(IN) :: xlat, xlon_in
  REAL(r_size), INTENT(OUT) :: lxap
  REAL(r_size) :: xlon

  ! Ensure the longitude is specified on a 0-360 grid:
  xlon = modulo(xlon_in,360.0)

  ! STEVE: Stolen from SODA: 
  ! ISSUE: use until we have a general method for managing land-blocked ocean basins...
  !=================================================================
  ! X. Cao 12/9/99
  !
  !   to make a mark to the location of a point in Caribbean area
  ! (xlat.gt.-2..and.xlat.lt.32..and.xlon.gt.245..and.xlon.lt.295.)
  ! to indicate if it is in Atlantic ocean (1) or in Pacific ocean (2)
  ! or out of Caribbean area (0)
  !=================================================================
  !
  lxap=0
  !
  ! -- Atlantic? Pacific?
  !
  if(xlat.gt.-2.0.and.xlat.le.8.5) then
     if(xlon.lt.285.) then
        lxap=2
     else
        lxap=1
     endif
  endif

  if(xlat.gt.8.5.and.xlat.le.15.5) then
     if(xlon.lt.276.) then
        lxap=2
     else
        lxap=1
     endif
  endif

  if(xlat.gt.15.5.and.xlat.le.19.5) then
     if(xlon.lt.270.) then
        lxap=2
     else
        lxap=1
     endif
  endif

  if(xlat.gt.19.5.and.xlat.le.32.0) then
     if(xlon.lt.258.) then
        lxap=2
     else
        lxap=1
     endif
  endif
END SUBROUTINE atlpac

END MODULE letkf_local
