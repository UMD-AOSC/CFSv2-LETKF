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

  !! the list of observations found and their distances
  !! are cached here, outside of obs_local,
  !! this greatly speeds things up by not having to search the kd
  !! tree again if we do all the vertical levels first for each grid point.
  integer,allocatable :: idx(:)
  real(r_size),allocatable :: dist(:)
  real(r_size) :: prev_lon, prev_lat
  integer :: nn
  
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
  REAL(r_size),INTENT(IN)  :: rlon,rlat,rlev
  INTEGER,INTENT(IN)       :: nvar
  REAL(r_size),INTENT(OUT) :: hdxf(nobstotal,nbv)
  REAL(r_size),INTENT(OUT) :: rdiag(nobstotal)
  REAL(r_size),INTENT(OUT) :: rloc(nobstotal)
  REAL(r_size),INTENT(OUT) :: dep(nobstotal)
  INTEGER,INTENT(OUT)      :: nobsl
  INTEGER,INTENT(OUT),OPTIONAL :: oindex(nobstotal)      ! DH

  REAL(r_size) :: logrlev
  INTEGER :: id, n, k
  logical :: b

  real(r_size), parameter :: loc_cutoff = 4.0*10.0/3.0
  real(r_size) :: sigma_max_h_ij, sigma_ocn_h_ij, sigma_atm_h_ij, sigma_h
  real(r_size) :: dlev
  real(r_size) :: loc, loc_h, loc_t, loc_atm_v, loc_ocn_v
  
  !!------------------------------------------------------------
  !!initialize the KD search tree if it hasn't been yet
  if (initialized == 0) then
     WRITE(6,*) "initializing the obs_local()"
     initialized = 1
     call kd_init(kdtree_root, obslon, obslat)
     write(6,*) "done constructing KD search tree"
     write(6,*) "nobstotal",nobstotal
     allocate(dist(nobstotal))
     allocate(idx(nobstotal))
     prev_lon = -1e10
     prev_lat = -1e10
  end if


  !! determine the max search radius for the initial obs search
  sigma_ocn_h_ij = (1.0d0-abs(rlat)/90.0d0) * (sigma_ocn_h(1)-sigma_ocn_h(2))+sigma_ocn_h(2)
  sigma_atm_h_ij = (1.0d0-abs(rlat)/90.0d0) * (sigma_atm_h(1)-sigma_atm_h(2))+sigma_atm_h(2)
  sigma_max_h_ij = max(sigma_ocn_h_ij, sigma_atm_h_ij) * sqrt(loc_cutoff)

  
  !!------------------------------------------------------------
  !! query the KD tree
  !!  if this was the same lat/lon that we saw last time this function
  !!  was called, just return the cached observations that were found
  !!  last time. This way, it should go much faster if we do DA one
  !!  column at a time.
  if (rlon /= prev_lon .or. rlat /= prev_lat) then
     call kd_search(kdtree_root, obslon, obslat, (/rlon, rlat/), &
          sigma_max_h_ij, idx, dist, nn)
     prev_lon = rlon
     prev_lat = rlat
  end if

  
  !! for each observation found in the radius, do the localization
  !! ------------------------------------------------------------
  logrlev = log(rlev)  
  nobsl = 0
  do n = 1, nn
     loc_h = 0.0
     loc_atm_v = 0.0
     loc_ocn_v = 0.0
     loc_t = 0.0
     
     if(nobsl >= nobstotal) return

     !! determine domain specific parameters
     !! --------------------------------------------------
     id = obselm(idx(n))

     if (atm_obs .and. getDomain(id) == dom_atm) then
        !! Atmospheric observations
        !! ------------------------------

        !! has the user explicitly defined a subset of atmospheric platforms to use
        !! TODO: this should be moved to earlier in the LETKF program
        b = .false.
        if (atm_obs_plat(1) > 0) then
           do k=1, size(atm_obs_plat)
              if (atm_obs_plat(k) <=0 ) exit
              if (atm_obs_plat(k) == obstyp(idx(n))) then
                 b = .true.
                 exit
              end if
           end do
           if (.not. b) cycle
        end if
        
        sigma_h = sigma_atm_h_ij
        if (id == obsid_atm_ps) then
           dlev = abs(log(obsdat(idx(n))) - logrlev)
        else
           dlev = abs(log(obslev(idx(n))) - logrlev)
        end if
        loc_atm_v = (dlev/sigma_atm_v) **2
        
     else if (ocn_obs .and. getDomain(id) == dom_ocn) then
        !! Ocean observations
        !! ------------------------------
        sigma_h = sigma_ocn_h_ij
        loc_ocn_v = (abs(obslev(idx(n)))/sigma_ocn_v) ** 2
        !! TODO, calculate the reference ps correctly
        loc_atm_v = (abs(log(1013.0d2)-logrlev)/sigma_atm_v) ** 2
     end if

     !! horizontal localization
     loc_h = (dist(n)/sigma_h) ** 2
    
     !! check combined localization
     loc = loc_h + loc_t + loc_atm_v + loc_ocn_v
     if (loc > loc_cutoff) cycle
     
     !! use the observation
     !! ------------------------------     
     nobsl = nobsl+1
     rloc(nobsl) = exp(-0.5d0 * loc)
     hdxf(nobsl,:) = obshdxf(idx(n), :)
     dep(nobsl) = obsdep(idx(n))
     rdiag(nobsl) = obserr(idx(n))**2
  end do
  return

END SUBROUTINE obs_local

END MODULE letkf_local
