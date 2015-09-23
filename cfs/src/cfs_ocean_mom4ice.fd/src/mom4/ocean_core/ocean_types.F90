!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify       !!
!! it and are expected to follow the terms of the GNU General Public !!
!! License as published by the Free Software Foundation.             !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS; if not, write to:                                 !!
!!          Free Software Foundation, Inc.                           !!
!!          59 Temple Place, Suite 330                               !!
!!          Boston, MA  02111-1307  USA                              !!
!! or see:                                                           !!
!!          http://www.gnu.org/licenses/gpl.txt                      !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! The include file fms_platform.h will handle the conversion of POINTER to ALLOCATABLE arrays
! for derived type members. The conversion affects performance only and should not change 
! any numeric result. It is limited to member arrays that are used within MOM4 only
! and to arrays that are never associated (=>) with another array.

! Fortran 90 requires arrays members of derived type to have the POINTER attribute. 
! However, most Fortran 95 compilers now also support ALLOCATABLE array components 
! (a Fortran 2003 feature). This avoids the aliasing problem afflicting pointers. 
! Some compilers may require an additional switch (e.g. -fno-alias) to fully exploit
! the performance benefit of the conversion.
!
! Macros used from fms_platform.h:
! _ALLOCATABLE maps to either POINTER  or ALLOCATABLE
! _NULL        maps to either =>NULL() or "nothing"

#include <fms_platform.h>

module ocean_types_mod
!
!<CONTACT EMAIL="Matthew.Harrison@noaa.gov"> Matt Harrison 
!</CONTACT>
!
!<CONTACT EMAIL="Stephen.Griffies@noaa.gov">
! S. M. Griffies 
!</CONTACT>
!
!<OVERVIEW>
! This module contains type declarations and default values for ocean model.
!</OVERVIEW>
!
!<DESCRIPTION>
! This module contains type declarations and default values for ocean model.
! Multiple model realizations need to be distinguished by
! an ensemble_id for use by the diag_manager.
!</DESCRIPTION>
!
  use fms_mod,          only: write_version_number
  use mpp_domains_mod,  only: domain2d
  use mpp_mod,          only: FATAL, mpp_error
  use time_manager_mod, only: time_type

  implicit none

  private

  logical :: module_is_initialized=.false.
  character(len=128) :: version = &
     '$Id$'
  character (len=128) :: tagname = &
     '$Name$'

  ! parameters for choosing the time tendency
  integer, parameter, public :: TWO_LEVEL   = 2
  integer, parameter, public :: THREE_LEVEL = 3

  ! parameters for tracer advection 
  integer, parameter, public :: ADVECT_UPWIND     = 1
  integer, parameter, public :: ADVECT_2ND_ORDER  = 2
  integer, parameter, public :: ADVECT_4TH_ORDER  = 3
  integer, parameter, public :: ADVECT_6TH_ORDER  = 4 
  integer, parameter, public :: ADVECT_QUICKER    = 5
  integer, parameter, public :: ADVECT_QUICKMOM3  = 6
  integer, parameter, public :: ADVECT_MDFL_SUP_B = 7
  integer, parameter, public :: ADVECT_MDFL_SWEBY = 8
  
  real,    parameter, public :: missing_value=-1.e10

  integer, parameter, public :: TEMP_ID = 1, SALT_ID = 2, SFC_HEIGHT_ID = 3


  type, public :: obc_flux
     real, _ALLOCATABLE, dimension(:,:) :: flux _NULL     ! flux through a boundary 
  end type obc_flux

  type, public :: ocean_time_steps_type
     character(len=32)  :: time_tendency  ! either "twolevel" or "threelevel" 
     integer :: tendency ! either "TWO_LEVEL" or "THREE_LEVEL"
     real    :: aidif    ! aidif=1.0 for fully implicit vertical mixing and 0.0 for fully explicit
     real    :: acor     ! acor=0.0 for explicit Coriolis, 0.5 <= acor <= 1.0 for semi-implicit
     real    :: dtts     ! tracer timestep (seconds)
     real    :: dtuv     ! internal mode timestep (seconds)
     real    :: dtfs     ! external mode timestep (seconds)
     real    :: dteta    ! ocean volume time step (eta_t) (seconds)
     real    :: dtime_t  ! 2*dtts  if threelevel, 1*dtts  if twolevel
     real    :: dtime_u  ! 2*dtuv  if threelevel, 1*dtuv  if twolevel
     real    :: dtime_e  ! 2*dteta if threelevel, 1*dteta if twolevel
  end type ocean_time_steps_type


#include <ocean_memory.h>

#ifdef STATIC_MEMORY
!################################################################################################

  type, public :: ocean_thickness_type
     real, dimension(isd:ied,jsd:jed,nk,3) :: dht   ! thickness of T cell (m) at three time levels
     real, dimension(isd:ied,jsd:jed,nk)   :: dhtr  ! 1.0/dht(taup1) (m^-1)
     real, dimension(isd:ied,jsd:jed,0:nk) :: dhwt  ! vertical distance between T points (m) at time tau
     real, dimension(isd:ied,jsd:jed,nk,3) :: dhu   ! thickness of U cell (m) at three time levels
     real, dimension(isd:ied,jsd:jed,nk)   :: dhur  ! 1.0/dhu(taup1) (m^-1)
     real, dimension(isd:ied,jsd:jed,0:nk) :: dhwu  ! vertical distance between U points (m) at time tau
     real, dimension(isd:ied,jsd:jed,nk)   :: ztp   ! depth to T cell grid point (m)
  end type ocean_thickness_type


  type, public :: ocean_grid_type
     character(len=32) :: name

     ! geometry and topology and rotation
     logical                           :: cyclic           ! true if domain is cyclic in the i direction
     logical                           :: solid_walls      ! for solid walls 
     logical                           :: tripolar         ! has folded connectivity at the "top" row within bipolar Arctic 
     logical                           :: beta_plane       ! beta plane Cartesian 
     logical                           :: f_plane          ! f-plane Cartesian
     real                              :: f_plane_latitude ! latitude where f_plane is centered  
     real, dimension(isd:ied,jsd:jed)  :: f                ! coriolis parameter at u-cell points (sec^-1)
     real, dimension(isd:ied,jsd:jed)  :: beta             ! planetary vorticity gradient df/dy at u-cell points (sec^-1*m^-1)     
     real, dimension(isd:ied,jsd:jed)  :: beta_eff         ! planetary beta plus topographic beta at u-cell points (sec^-1*m^-1)     


     ! vertical grid information (time independent) 
   
     integer, dimension(isd:ied,jsd:jed) :: kmt    ! number of t-levels
     integer, dimension(isd:ied,jsd:jed) :: kmu    ! number of u-levels
     real, dimension(nk)                 :: zt     ! distance from surface down to grid point in level k (m) 
     real, dimension(nk)                 :: zw     ! distance from surface down to bottom of level k (m) 
     real, dimension(nk)                 :: dzt    ! vertical resolution of T or U grid cells (m) 
     real, dimension(0:nk)               :: dzw    ! vertical resolution of W grid cells (m) 
     real, dimension(0:nk)               :: dzwr   ! reciprocal of dzw (W cell vertical resolution)
     real, dimension(nk,0:1)             :: fracdz ! fractional distance between grid point and cell top/bot 
     real, dimension(isd:ied,jsd:jed)    :: ht     ! depth to bottom of ocean (m) on t-cells at time tau
     real, dimension(isd:ied,jsd:jed)    :: hu     ! depth to bottom of ocean (m) on u-cells at time tau 
     real, dimension(isd:ied,jsd:jed)    :: dht_dx ! d(ht)/dx on u-cells (m/m) at time tau
     real, dimension(isd:ied,jsd:jed)    :: dht_dy ! d(ht)/dy on u-cells (m/m) at time tau

     ! horizontal grid information (time independent)
   
     real, dimension(isd:ied,jsd:jed) :: xt         ! longitude of the T grid points in degrees.
     real, dimension(isd:ied,jsd:jed) :: xu         ! longitude of the U grid points in degrees.
     real, dimension(isd:ied,jsd:jed) :: yt         ! latitude of the T grid points in degrees.
     real, dimension(isd:ied,jsd:jed) :: yu         ! latitude of the U grid points in degrees.
     real, dimension(isd:ied,jsd:jed) :: phiu       ! latitude of U grid point in radians
     real, dimension(isd:ied,jsd:jed) :: phit       ! latitude of T grid point in radians
     real, dimension(isd:ied,jsd:jed) :: h1t        ! metric factors in i grid direction
     real, dimension(isd:ied,jsd:jed) :: h1u        ! metric factors in j grid jirection
     real, dimension(isd:ied,jsd:jed) :: h2t        ! metric factors in i grid direction
     real, dimension(isd:ied,jsd:jed) :: h2u        ! metric factors in j grid jirection
     real, dimension(isd:ied,jsd:jed) :: dh2dx      ! (1/delta_y)*d(delta_y)/dx (1/m)
     real, dimension(isd:ied,jsd:jed) :: dh1dy      ! (1/delta_x)*d(delta_x)/dy (1/m)
     real, dimension(isd:ied,jsd:jed) :: dxt        ! longitudinal width of T-cells at grid point (m)
     real, dimension(isd:ied,jsd:jed) :: dxu        ! longitudinal width of U-cells at grid point (m)
     real, dimension(isd:ied,jsd:jed) :: dyt        ! latitudinal width of T-cells at grid point (m) 
     real, dimension(isd:ied,jsd:jed) :: dyu        ! latitudinal width of U-cells at grid point (m)
     real, dimension(isd:ied,jsd:jed) :: dau        ! area of U-cells (m^2)
     real, dimension(isd:ied,jsd:jed) :: dat        ! area of T-cells (m^2)
     real, dimension(isd:ied,jsd:jed) :: dxte       ! longitudinal width between grid points at i+1 and i in T-cells (m)
     real, dimension(isd:ied,jsd:jed) :: dxtn       ! longitudinal width of north face of T-cells (m)
     real, dimension(isd:ied,jsd:jed) :: dyte       ! latitudinal width of east  face of T-cells (m)
     real, dimension(isd:ied,jsd:jed) :: dytn       ! latitudinal width between grid points at j+1 and j in T-cells (m)
     real, dimension(isd:ied,jsd:jed) :: dxue       ! longitudinal width between grid points at i+1 and i in U-cells (m)
     real, dimension(isd:ied,jsd:jed) :: dxun       ! longitudinal width of north face of U-cells (m)
     real, dimension(isd:ied,jsd:jed) :: dyue       ! latitudinal width of east  face of U-cells (m)
     real, dimension(isd:ied,jsd:jed) :: dyun       ! latitudinal width between grid points at j+1 and j in U-cells (m)
     real, dimension(isd:ied,jsd:jed) :: datnr      ! reciprocal area at north face of T-cell
     real, dimension(isd:ied,jsd:jed) :: dater      ! reciprocal area at east face of T-cell
     real, dimension(isd:ied,jsd:jed) :: dun        ! width from grid point to north face of U-cells (m)
     real, dimension(isd:ied,jsd:jed) :: dus        ! width from grid point to south face of U-cells (m)
     real, dimension(isd:ied,jsd:jed) :: duw        ! width from grid point to west  face of U-cells (m)
     real, dimension(isd:ied,jsd:jed) :: due        ! width from grid point to east  face of U-cells (m)
     real, dimension(isd:ied,jsd:jed) :: dtn        ! width from grid point to north face of T-cells (m)
     real, dimension(isd:ied,jsd:jed) :: dts        ! width from grid point to south face of T-cells (m)
     real, dimension(isd:ied,jsd:jed) :: dtw        ! width from grid point to west  face of T-cells (m)
     real, dimension(isd:ied,jsd:jed) :: dte        ! width from grid point to east  face of T-cells (m)
     real, dimension(isd:ied,jsd:jed) :: dxtr       ! 1/dxt
     real, dimension(isd:ied,jsd:jed) :: dxur       ! 1/dxu
     real, dimension(isd:ied,jsd:jed) :: dytr       ! 1/dyt
     real, dimension(isd:ied,jsd:jed) :: dyur       ! 1/dyu
     real, dimension(isd:ied,jsd:jed) :: daur       ! 1/[area of U-cells (m^2)]
     real, dimension(isd:ied,jsd:jed) :: datr       ! 1/[area of T-cells (m^2)]
     real, dimension(isd:ied,jsd:jed) :: dxter      ! 1/dxte
     real, dimension(isd:ied,jsd:jed) :: dytnr      ! 1/dytn
     real, dimension(isd:ied,jsd:jed) :: dxuer      ! 1/dxue
     real, dimension(isd:ied,jsd:jed) :: dxunr      ! 1/dxun
     real, dimension(isd:ied,jsd:jed) :: dyuer      ! 1/dyue
     real, dimension(isd:ied,jsd:jed) :: dyunr      ! 1/dyun
     real, dimension(isd:ied,jsd:jed) :: dyue_dxuer ! dyue/dxue 
     real, dimension(isd:ied,jsd:jed) :: dxun_dyunr ! dxun/dyun

     ! land/sea masks 
     real, dimension(isd:ied,jsd:jed)     :: obc_tmask ! land/sea mask for T cell diagnostics with OBC 
     real, dimension(isd:ied,jsd:jed)     :: obc_umask ! land/sea mask for U cell diagnostics with OBC
     real, dimension(isd:ied,jsd:jed,nk)  :: tmask     ! land/sea mask for T cells
     real, dimension(isd:ied,jsd:jed,nk)  :: umask     ! land/sea mask for U cells

     ! grid areas and volumes 
     real                 :: tcellv    ! T cell volume m^3 (entire ocean with eta_t=0)
     real                 :: ucellv    ! U cell volume m^3 (entire ocean with eta_u=0)
     real                 :: tcellsurf ! T cell surface area (k=1) in bitwise reproducible form
     real                 :: ucellsurf ! U cell surface area (k=1) in bitwise reproducible form
     real, dimension(nk)  :: tcella    ! T cell surface area m^2 (entire ocean)
     real, dimension(nk)  :: ucella    ! U cell surface area m^2 (entire ocean)

     ! sine and cosine of rotation angles (clockwise) of velocity for tripolar 
     real, dimension(isd:ied,jsd:jed) :: sin_rot  
     real, dimension(isd:ied,jsd:jed) :: cos_rot  

     ! 1-d grid coordinates for COARDS NetCDF files
     real, dimension(ni) :: grid_x_t  
     real, dimension(nj) :: grid_y_t 
     real, dimension(ni) :: grid_x_u  
     real, dimension(nj) :: grid_y_u  

     ! axes id for diagnostic manager 
     integer, dimension(3)  :: tracer_axes        
     integer, dimension(3)  :: vel_axes_uv         
     integer, dimension(3)  :: vel_axes_wu    
     integer, dimension(3)  :: vel_axes_wt    
     integer, dimension(3)  :: tracer_axes_wt 
     integer, dimension(3)  :: tracer_axes_flux_x  
     integer, dimension(3)  :: tracer_axes_flux_y  
     integer, dimension(3)  :: vel_axes_flux_x    
     integer, dimension(3)  :: vel_axes_flux_y    
     integer                :: nk     ! number of vertical grid points 
     integer                :: ni, nj ! number of global points in the two horizontal directions
  end type ocean_grid_type


  type, public :: ocean_domain_type
     type(domain2d) :: domain2d           ! fms variable, used by mpp routines
     integer :: isc, iec, jsc, jec        ! computational domain indices 
     integer :: isd, ied, jsd, jed        ! local indices including halo, consistent with domain2d
     integer :: isg, ieg, jsg, jeg        ! global indices
     integer :: isa, iea, jsa, jea        ! active indices (for minimizing comm2d calls)
     integer :: xhalo, yhalo              ! halo sizes 
     integer :: xflags, yflags, layout(2) ! options to mpp_define_domains
     integer :: ioff , joff               ! index offset to get absolute indices if using static allocations (0 otherwise)
  end type ocean_domain_type

  type, public :: ocean_time_type
     type(time_type) :: model_time     ! fms variable
     type(time_type) :: Time_step      ! time step for tracers (and for ocean model)
     type(time_type) :: Time_init      ! initial time 
     integer :: calendar               ! calendar type defined by time_manager_mod
     logical :: init                   ! true at beginning of run
     integer :: itt                    ! timestep counter
     integer :: taum1, tau, taup1      ! time level indices
     integer :: tau_m2, tau_m1, tau_m0 ! time level indices for Adams-Bashforth velocity advection 
  end type ocean_time_type

  type, public :: ocean_adv_vel_type
     real, dimension(isd:ied,jsd:jed,nk)   :: uh_et ! thickness-weighted advective velocity component (m^2/sec) on i-face of T-cell
     real, dimension(isd:ied,jsd:jed,nk)   :: vh_nt ! thickness-weighted advective velocity component (m^2/sec) on j-face of T-cell
     real, dimension(isd:ied,jsd:jed,nk)   :: uh_eu ! thickness-weighted advective velocity component (m^2/sec) on i-face of U-cell
     real, dimension(isd:ied,jsd:jed,nk)   :: vh_nu ! thickness-weighted advective velocity component (m^2/sec) on j-face of U-cell
     real, dimension(isd:ied,jsd:jed,0:nk) :: w_bt  ! vertical advective velocity component on bottom-face (m/sec) of T-cell
     real, dimension(isd:ied,jsd:jed,0:nk) :: w_bu  ! vertical advective velocity component on bottom-face (m/sec) of U-cell
  end type ocean_adv_vel_type

  type , public :: ocean_density_type
     real, dimension(isd:ied,jsd:jed,nk)   :: rho                ! in situ density (kg/m^3) at time tau   
     real, dimension(isd:ied,jsd:jed,nk)   :: rho_taum1          ! in situ density (kg/m^3) at time taum1 
     real, dimension(isd:ied,jsd:jed,nk)   :: rho_taup1          ! in situ density (kg/m^3) at time taup1 (used in twolevel scheme)
     real, dimension(isd:ied,jsd:jed,nk)   :: pressure_at_depth  ! hydrostatic pressure due to overlying fluid (including patm)
     real, dimension(isd:ied,jsd:jed,nk)   :: potrho             ! potential density (kg/m^3)
     integer, dimension(3)                 :: potrho_axes        ! axis ids for diagnosing potential density 
     integer, dimension(3)                 :: theta_axes         ! axis ids for potential temperature 
     real, _ALLOCATABLE, dimension(:)      :: potrho_ref _NULL   ! partition vertical into potrho-classes 
     real, _ALLOCATABLE, dimension(:)      :: theta_ref  _NULL   ! partition vertical into theta-classes 
  end type ocean_density_type

  
  type, public :: ocean_prog_tracer_type
     character(len=32)  :: name
     character(len=32)  :: units
     character(len=128) :: longname

     logical :: use_only_advection        ! for testing purposes, evolve using ONLY advection
     logical :: complete                  ! to determine if ready to do mpp updates
     logical :: neutral_physics_limit     ! revert neutral physics to horz diffusion where tracer out of bounds

     integer :: sfc_flux_id=-1            ! index for time_interp_external
     integer :: horz_advect_scheme=-1     ! id for horizontal advection scheme
     integer :: vert_advect_scheme=-1     ! id for vertical advection scheme

     type(obc_flux), dimension(:), _ALLOCATABLE :: otf   _NULL   ! flux through open boundaries, allocate nobc
     real, dimension(isd:ied,jsd:jed,nk,3) :: field            ! tracer concentration at 3 time levels
     real, dimension(isd:ied,jsd:jed,nk)   :: th_tendency      ! thickness weighted tendency
     real, dimension(isd:ied,jsd:jed,nk)   :: wrk1             ! work array
     real, dimension(isd:ied,jsd:jed,nk)   :: tmask_limit      ! for limiting fluxes from advection and/or neutral physics 
     real, dimension(isd:ied,jsd:jed,nk)   :: K33_implicit     ! m^2/sec vertical diffusivity from neutral diffusion 
     real, dimension(isd:ied,jsd:jed)      :: stf              ! surface tracer flux in units m/sec*tracer concentration
     real, dimension(isd:ied,jsd:jed)      :: btf              ! bottom tracer flux in units m/sec*tracer concentration
     real, dimension(isd:ied,jsd:jed)      :: tpme, triver     ! tracer concentration in precip-evap and river water 
     real, dimension(isd:ied,jsd:jed)      :: riverdiffuse     ! sets where to enhance diff_cbt according to rivers
     real, dimension(isd:ied,jsd:jed)      :: surface_smooth   ! tendency source from local surface height smoothing 

     real                                  :: conversion       ! conversion from model specific dimensions to other useful dimensions 
     real                                  :: offset           ! offset in dimensions (e.g., Celsius to Kelvin)
     real                                  :: min_tracer       ! min acceptable value--model brought down if less than min_tracer
     real                                  :: max_tracer       ! max acceptable value--model brought down if greater than max_tracer
     real                                  :: min_range        ! min value used for calls to diagnostic manager 
     real                                  :: max_range        ! max value used for calls to diagnostic manager 
     real                                  :: min_tracer_limit ! min value used to limit quicker and neutral fluxes
     real                                  :: max_tracer_limit ! max value used to limit quicker and neutral fluxes 
     real                                  :: min_flux_range   ! min and max values used for flux diagnostics
     real                                  :: max_flux_range   ! min and max values used for flux diagnostics
     real                                  :: scale_in         ! value by which to scale input restart data    
     real                                  :: additive_in      ! uniform value added to input restart data 
     real                                  :: const_init_value ! value used to initialize tracer when using constant_init_tracer
     logical                               :: const_init_tracer! to initialize tracer to constant value 
     logical                               :: init             ! true if the input restart file is an initial condition file
     character(len=32)                     :: flux_units       ! units for the tracer flux 
     character(len=128)                    :: file_in          ! name for input restart file
     character(len=128)                    :: file_out         ! name for output restart file
     character(len=32)                     :: name_in          ! name of variable to use from the input restart file
  end type ocean_prog_tracer_type

  
  type, public :: ocean_diag_tracer_type
     character(len=32)  :: name, units
     character(len=128) :: longname
     real, dimension(isd:ied,jsd:jed,nk) :: field              ! tracer concentration at single time level 
     real :: factor                                            ! allows for a possibly different time stepping schemes in ocn and sea ice
                                                               ! (useful in particular for frazil)
     real :: conversion                                        ! conversion from model specific dimensions to other useful dimensions 
     real :: offset                                            ! offset in dimensions (e.g., Celsius to Kelvin)
     real :: min_tracer, max_tracer                            ! min and max acceptable values used for error checking 
     real :: min_range, max_range                              ! min and max values used for diagnostics
     logical            :: init                                ! true if the input restart file is an initial condition file
     character(len=128) :: file_in                             ! name for input restart file
     character(len=128) :: file_out                            ! name for output restart file
     character(len=32)  :: name_in                             ! name of variable to use from the input restart file
     real               :: scale_in                            ! value by which to scale input restart data
     real               :: additive_in                         ! uniform value added to input restart data
     real               :: const_init_value                    ! value used to initialize tracer when using constant_init_tracer
     logical            :: const_init_tracer                   ! to initialize tracer to constant value 
  end type ocean_diag_tracer_type


  type, public :: ocean_velocity_type
     real, dimension(isd:ied,jsd:jed,nk,2,3)  :: u         ! horz velocity (m/sec) in i,j directions at 3 time levels
     real, dimension(isd:ied,jsd:jed,2)       :: smf       ! momentum flux per mass across ocean surface (m^2/sec^2)
     real, dimension(isd:ied,jsd:jed,2)       :: bmf       ! momentum flux per mass across ocean bottom  (m^2/sec^2)
     real, dimension(isd:ied,jsd:jed,nk,2,3)  :: advection ! thickness weighted tendency from velocity advection (m^2/sec^2) 
  end type ocean_velocity_type


  type, public :: ocean_external_mode_type  
     real, dimension(isd:ied,jsd:jed,3)   :: eta_t          ! surface height on tracer cell center (m) 
     real, dimension(isd:ied,jsd:jed,3)   :: eta_u          ! surface height on velocity cell center  (m) 
     real, dimension(isd:ied,jsd:jed,3)   :: eta_t_bar      ! surface height on tracer cell time avg over ext-mode time steps (m) 
     real, dimension(isd:ied,jsd:jed)     :: deta_dt        ! surface height time tendency on t-cell (m/s)
     real, dimension(isd:ied,jsd:jed,3)   :: convud_t       ! convergence of ud on T-cell center (m/s)
     real, dimension(isd:ied,jsd:jed,2,3) :: ud             ! vertically integrated horizontal velocity (m^2/s)
     real, dimension(isd:ied,jsd:jed,2)   :: forcing_fs     ! depth integrated time change of velocity (without coriolis)
     real, dimension(isd:ied,jsd:jed)     :: ps             ! surface pressure (pressure at z=0 due to eta_t) at time tau
     real, dimension(isd:ied,jsd:jed,2)   :: grad_ps        ! rho0r * surface pressure gradient at time tau
     real, dimension(isd:ied,jsd:jed)     :: eta_source     ! eta sources (not including freshwater flux)
     real, dimension(isd:ied,jsd:jed)     :: surface_smooth ! tendency (m/s) from local smoothing of surface height 
  end type ocean_external_mode_type


#else
!############################################################################################
! not STATIC_MEMORY

  type, public :: ocean_thickness_type
     real, dimension(:,:,:,:), _ALLOCATABLE :: dht   _NULL ! thickness of T cell (m) at three time levels
     real, dimension(:,:,:), _ALLOCATABLE   :: dhtr  _NULL ! 1.0/dht(taup1) (m^-1)
     real, dimension(:,:,:), _ALLOCATABLE   :: dhwt  _NULL ! vertical distance between T points (m) at time tau 
     real, dimension(:,:,:,:), _ALLOCATABLE :: dhu   _NULL ! thickness of U cell (m) at three time levels
     real, dimension(:,:,:), _ALLOCATABLE   :: dhur  _NULL ! 1.0/dhu(taup1) (m^-1)
     real, dimension(:,:,:), _ALLOCATABLE   :: dhwu  _NULL ! vertical distance between U points (m) at time tau 
     real, dimension(:,:,:), _ALLOCATABLE   :: ztp   _NULL ! depth to T cell grid point (m) at time tau 
  end type ocean_thickness_type


  type, public :: ocean_grid_type
     character(len=32) :: name

     ! geometry and topology and rotation
     logical                       :: cyclic            ! true if domain is cyclic in the i direction
     logical                       :: solid_walls       ! for solid walls 
     logical                       :: tripolar          ! has folded connectivity at the "top" row within bipolar Arctic 
     logical                       :: beta_plane        ! beta plane Cartesian 
     logical                       :: f_plane           ! f-plane Cartesian
     real                          :: f_plane_latitude  ! latitude where f_plane is centered  
     real, dimension(:,:), _ALLOCATABLE :: f        _NULL ! coriolis parameter at u-cell points (sec^-1)
     real, dimension(:,:), _ALLOCATABLE :: beta     _NULL ! planetary vorticity gradient df/dy at u-cell points (sec^-1*m^-1)     
     real, dimension(:,:), _ALLOCATABLE :: beta_eff _NULL ! planetary beta plus topographic beta at u-cell points (sec^-1*m^-1)     

     ! vertical grid information (time independent) 
    
     integer, dimension(:,:), _ALLOCATABLE :: kmt    _NULL ! number of t-levels
     integer, dimension(:,:), _ALLOCATABLE :: kmu    _NULL ! number of u-levels
     real,    dimension(:),   _ALLOCATABLE :: zt     _NULL ! distance from surface down to grid point in level k (m) 
     real,    dimension(:),   _ALLOCATABLE :: zw     _NULL ! distance from surface down to bottom of level k (m) 
     real,    dimension(:),   _ALLOCATABLE :: dzt    _NULL ! vertical resolution of T or U grid cells (m) 
     real,    dimension(:),   _ALLOCATABLE :: dzw    _NULL ! vertical resolution of W grid cells (m) 
     real,    dimension(:),   _ALLOCATABLE :: dzwr   _NULL ! reciprocal of dzw (W cell vertical resolution)
     real,    dimension(:,:), _ALLOCATABLE :: fracdz _NULL ! fractional distance between grid point and cell top/bot 
     real,    dimension(:,:), _ALLOCATABLE :: ht     _NULL ! depth to bottom of ocean (m) on t-cells at time tau 
     real,    dimension(:,:), _ALLOCATABLE :: hu     _NULL ! depth to bottom of ocean (m) on u-cells at time tau  
     real,    dimension(:,:), _ALLOCATABLE :: dht_dx _NULL ! d(ht)/dx on u-cells (m/m) at time tau  
     real,    dimension(:,:), _ALLOCATABLE :: dht_dy _NULL ! d(ht)/dy on u-cells (m/m) at time tau  

     ! horizontal grid information (time independent)
    
     real, dimension(:,:), _ALLOCATABLE :: xt         _NULL ! longitude of the T grid points in degrees.
     real, dimension(:,:), _ALLOCATABLE :: xu         _NULL ! longitude of the U grid points in degrees.
     real, dimension(:,:), _ALLOCATABLE :: yt         _NULL ! latitude of the T grid points in degrees.
     real, dimension(:,:), _ALLOCATABLE :: yu         _NULL ! latitude of the U grid points in degrees.
     real, dimension(:,:), _ALLOCATABLE :: phiu       _NULL ! latitude of U grid point in radians
     real, dimension(:,:), _ALLOCATABLE :: phit       _NULL ! latitude of T grid point in radians
     real, dimension(:,:), _ALLOCATABLE :: h1t        _NULL ! metric factors in i grid direction
     real, dimension(:,:), _ALLOCATABLE :: h1u        _NULL ! metric factors in j grid jirection
     real, dimension(:,:), _ALLOCATABLE :: h2t        _NULL ! metric factors in i grid direction
     real, dimension(:,:), _ALLOCATABLE :: h2u        _NULL ! metric factors in j grid jirection
     real, dimension(:,:), _ALLOCATABLE :: dh2dx      _NULL ! (1/delta_y)*d(delta_y)/dx (1/m)
     real, dimension(:,:), _ALLOCATABLE :: dh1dy      _NULL ! (1/delta_x)*d(delta_x)/dy (1/m)
     real, dimension(:,:), _ALLOCATABLE :: dxt        _NULL ! longitudinal width of T-cells at grid point (m)
     real, dimension(:,:), _ALLOCATABLE :: dxu        _NULL ! longitudinal width of U-cells at grid point (m)
     real, dimension(:,:), _ALLOCATABLE :: dyt        _NULL ! latitudinal width of T-cells at grid point (m) 
     real, dimension(:,:), _ALLOCATABLE :: dyu        _NULL ! latitudinal width of U-cells at grid point (m)
     real, dimension(:,:), _ALLOCATABLE :: dau        _NULL ! area of U-cells (m^2)
     real, dimension(:,:), _ALLOCATABLE :: dat        _NULL ! area of T-cells (m^2)
     real, dimension(:,:), _ALLOCATABLE :: dxte       _NULL ! longitudinal width between grid points at i+1 and i in T-cells (m)
     real, dimension(:,:), _ALLOCATABLE :: dxtn       _NULL ! longitudinal width of north face of T-cells (m)
     real, dimension(:,:), _ALLOCATABLE :: dyte       _NULL ! latitudinal width of east  face of T-cells (m)
     real, dimension(:,:), _ALLOCATABLE :: dytn       _NULL ! latitudinal width between grid points at j+1 and j in T-cells (m)
     real, dimension(:,:), _ALLOCATABLE :: dxue       _NULL ! longitudinal width between grid points at i+1 and i in U-cells (m)
     real, dimension(:,:), _ALLOCATABLE :: dxun       _NULL ! longitudinal width of north face of U-cells (m)
     real, dimension(:,:), _ALLOCATABLE :: dyue       _NULL ! latitudinal width of east  face of U-cells (m)
     real, dimension(:,:), _ALLOCATABLE :: dyun       _NULL ! latitudinal width between grid points at j+1 and j in U-cells (m)
     real, dimension(:,:), _ALLOCATABLE :: datnr      _NULL ! reciprocal area at north face of T-cell
     real, dimension(:,:), _ALLOCATABLE :: dater      _NULL ! reciprocal area at east face of T-cell
     real, dimension(:,:), _ALLOCATABLE :: dun        _NULL ! width from grid point to north face of U-cells (m)
     real, dimension(:,:), _ALLOCATABLE :: dus        _NULL ! width from grid point to south face of U-cells (m)
     real, dimension(:,:), _ALLOCATABLE :: duw        _NULL ! width from grid point to west  face of U-cells (m)
     real, dimension(:,:), _ALLOCATABLE :: due        _NULL ! width from grid point to east  face of U-cells (m)
     real, dimension(:,:), _ALLOCATABLE :: dtn        _NULL ! width from grid point to north face of T-cells (m)
     real, dimension(:,:), _ALLOCATABLE :: dts        _NULL ! width from grid point to south face of T-cells (m)
     real, dimension(:,:), _ALLOCATABLE :: dtw        _NULL ! width from grid point to west  face of T-cells (m)
     real, dimension(:,:), _ALLOCATABLE :: dte        _NULL ! width from grid point to east  face of T-cells (m)
     real, dimension(:,:), _ALLOCATABLE :: dxtr       _NULL ! 1/dxt
     real, dimension(:,:), _ALLOCATABLE :: dxur       _NULL ! 1/dxu
     real, dimension(:,:), _ALLOCATABLE :: dytr       _NULL ! 1/dyt
     real, dimension(:,:), _ALLOCATABLE :: dyur       _NULL ! 1/dyu
     real, dimension(:,:), _ALLOCATABLE :: daur       _NULL ! 1/[area of U-cells (m^2)]
     real, dimension(:,:), _ALLOCATABLE :: datr       _NULL ! 1/[area of T-cells (m^2)]
     real, dimension(:,:), _ALLOCATABLE :: dxter      _NULL ! 1/dxte
     real, dimension(:,:), _ALLOCATABLE :: dytnr      _NULL ! 1/dytn
     real, dimension(:,:), _ALLOCATABLE :: dxuer      _NULL ! 1/dxue
     real, dimension(:,:), _ALLOCATABLE :: dxunr      _NULL ! 1/dxun
     real, dimension(:,:), _ALLOCATABLE :: dyuer      _NULL ! 1/dyue
     real, dimension(:,:), _ALLOCATABLE :: dyunr      _NULL ! 1/dyun
     real, dimension(:,:), _ALLOCATABLE :: dyue_dxuer _NULL  ! dyue/dxue 
     real, dimension(:,:), _ALLOCATABLE :: dxun_dyunr _NULL  ! dxun/dyun

     ! land/sea masks 
     real, dimension(:,:),   _ALLOCATABLE :: obc_tmask _NULL ! land/sea mask for T cell diagnostics with OBC 
     real, dimension(:,:),   _ALLOCATABLE :: obc_umask _NULL ! land/sea mask for U cell diagnostics with OBC
     real, dimension(:,:,:), _ALLOCATABLE :: tmask     _NULL ! land/sea mask for T cells
     real, dimension(:,:,:), _ALLOCATABLE :: umask     _NULL ! land/sea mask for U cells

     ! grid areas and volumes 
     real                         :: tcellv             ! T cell volume m^3 (entire ocean) at time tau 
     real                         :: ucellv             ! U cell volume m^3 (entire ocean) at time tau 
     real                         :: tcellsurf          ! T cell surface area (k=1) in bitwise reproducible form
     real                         :: ucellsurf          ! U cell surface area (k=1) in bitwise reproducible form
     real, dimension(:), _ALLOCATABLE  :: tcella   _NULL  ! T cell surface area m^2 (entire ocean)
     real, dimension(:), _ALLOCATABLE  :: ucella   _NULL  ! U cell surface area m^2 (entire ocean)

     ! sine and cosine of rotation angles (clockwise) of velocity for tripolar 
     real, dimension(:,:), _ALLOCATABLE :: sin_rot  _NULL  
     real, dimension(:,:), _ALLOCATABLE :: cos_rot  _NULL  

     ! 1-d grid coordinates for COARDS NetCDF files
     real, dimension(:), _ALLOCATABLE  :: grid_x_t _NULL 
     real, dimension(:), _ALLOCATABLE  :: grid_y_t _NULL 
     real, dimension(:), _ALLOCATABLE  :: grid_x_u _NULL 
     real, dimension(:), _ALLOCATABLE  :: grid_y_u _NULL 

     ! axes id for diagnostic manager 
     integer, dimension(3)  :: tracer_axes        
     integer, dimension(3)  :: vel_axes_uv         
     integer, dimension(3)  :: vel_axes_wu    
     integer, dimension(3)  :: vel_axes_wt    
     integer, dimension(3)  :: tracer_axes_wt 
     integer, dimension(3)  :: tracer_axes_flux_x  
     integer, dimension(3)  :: tracer_axes_flux_y  
     integer, dimension(3)  :: vel_axes_flux_x    
     integer, dimension(3)  :: vel_axes_flux_y    
     integer                :: nk      ! number of vertical grid points 
     integer                :: ni, nj  ! number of global points in the two horizontal directions
  end type ocean_grid_type
  

  type, public :: ocean_domain_type
     type(domain2d) :: domain2d           ! fms variable, used by mpp routines
     integer :: isc, iec, jsc, jec        ! computational domain indices 
     integer :: isd, ied, jsd, jed        ! local indices, consistent with domain2d
     integer :: isg, ieg, jsg, jeg        ! global indices
     integer :: isa, iea, jsa, jea        ! active indices (for minimizing comm2d calls)
     integer :: xhalo, yhalo              ! halo sizes 
     integer :: xflags, yflags, layout(2) ! options to mpp_define_domains
     integer :: ioff, joff                ! offset to get absolute indices when
                                          ! using static allocations (0 otherwise)
  end type ocean_domain_type

  type, public :: ocean_time_type
     type(time_type) :: model_time     ! fms variable
     type(time_type) :: time_step      ! ocean tracer timestep
     type(time_type) :: Time_init      ! initial ocean time 
     integer :: calendar               ! calendar type defined by time_manager_mod
     logical :: init                   ! true at beginning of run
     integer :: itt                    ! timestep counter
     integer :: taum1, tau, taup1      ! time level indices
     integer :: tau_m2, tau_m1, tau_m0 ! time level indices for Adams-Bashforth velocity advection 
  end type ocean_time_type

  type, public :: ocean_adv_vel_type
     real, _ALLOCATABLE, dimension(:,:,:)  :: uh_et _NULL ! thickness-weight advective velocity (m^2/sec) on i-face of T-cell
     real, _ALLOCATABLE, dimension(:,:,:)  :: vh_nt _NULL ! thickness-weight advective velocity (m^2/sec) on j-face of T-cell
     real, _ALLOCATABLE, dimension(:,:,:)  :: uh_eu _NULL ! thickness-weight advective velocity (m^2/sec) on i-face of U-cell
     real, _ALLOCATABLE, dimension(:,:,:)  :: vh_nu _NULL ! thickness-weight advective velocity (m^2/sec) on j-face of U-cell
     real, _ALLOCATABLE, dimension(:,:,:)  :: w_bt  _NULL ! vertical advective velocity on bottom-face (m/sec) of T-cell
     real, _ALLOCATABLE, dimension(:,:,:)  :: w_bu  _NULL ! vertical advective velocity on bottom-face (m/sec) of U-cell
  end type ocean_adv_vel_type

  type, public ::  ocean_density_type
     real, _ALLOCATABLE, dimension(:,:,:) :: rho               _NULL ! in situ density (kg/m^3) at time tau      
     real, _ALLOCATABLE, dimension(:,:,:) :: rho_taum1         _NULL ! in situ density (kg/m^3) at time taum1 
     real, _ALLOCATABLE, dimension(:,:,:) :: rho_taup1         _NULL ! in situ density (kg/m^3) at time taup1 (used in twolevel scheme) 
     real, _ALLOCATABLE, dimension(:,:,:) :: pressure_at_depth _NULL ! hydrostatic pressure due to overlying fluid (including patm)
     real, _ALLOCATABLE, dimension(:,:,:) :: potrho            _NULL ! potential density (kg/m^3)
     real, _ALLOCATABLE, dimension(:)     :: potrho_ref        _NULL ! partition vertical into potrho-classes 
     real, _ALLOCATABLE, dimension(:)     :: theta_ref         _NULL ! partition vertical into theta-classes 
     integer, dimension(3)           :: potrho_axes                ! axis ids for diagnosing potential density 
     integer, dimension(3)           :: theta_axes                 ! axis ids for potential temperature 
  end type ocean_density_type
  
  type, public :: ocean_prog_tracer_type
     character(len=32)  :: name
     character(len=32)  :: units
     character(len=128) :: longname

     logical :: use_only_advection      ! for testing purposes, evolve using ONLY advection
     logical :: neutral_physics_limit   ! revert neutral physics to horz diffusion where tracer out of bounds
     logical :: complete                ! to determine if ready to do mpp updates

     integer :: sfc_flux_id=-1          ! index for time_interp_external
     integer :: horz_advect_scheme=-1   ! id for horizontal advection scheme
     integer :: vert_advect_scheme=-1   ! id for vertical advection scheme
     integer :: id_obc                  ! id to identify tracer in OBC-subroutines

     real, _ALLOCATABLE, dimension(:,:,:,:) :: field _NULL         ! tracer concentration at 3 time levels
     real, _ALLOCATABLE, dimension(:,:,:)   :: th_tendency _NULL   ! thickness weighted tracer tendency
     real, _ALLOCATABLE, dimension(:,:)     :: surface_smooth _NULL ! tendency source from local surface height smoothing 
     real, _ALLOCATABLE, dimension(:,:,:)   :: wrk1  _NULL         ! work array
     real, _ALLOCATABLE, dimension(:,:,:)   :: tmask_limit _NULL   ! for limiting fluxes from advection and/or neutral physics
     real, _ALLOCATABLE, dimension(:,:,:)   :: K33_implicit _NULL  ! m^2/sec vertical diffusivity from neutral diffusion 
     real, dimension(:,:), _ALLOCATABLE     :: stf   _NULL         ! surface tracer flux in units m/sec*tracer concentration
     real, dimension(:,:), _ALLOCATABLE     :: btf   _NULL         ! bottom tracer flux in units m/sec*tracer concentration
     real, dimension(:,:), _ALLOCATABLE     :: tpme  _NULL         ! tracer concentration in precip-evap and river water
     real, dimension(:,:), _ALLOCATABLE     :: triver _NULL        ! tracer concentration in precip-evap and river water
     real, dimension(:,:), _ALLOCATABLE     :: riverdiffuse _NULL  ! where to enhance diff_cbt according to rivers
     real, dimension(:,:), _ALLOCATABLE     :: flux_int _NULL      ! integrated sfc tracer flux for diagnostics
     type(obc_flux), dimension(:), _ALLOCATABLE   :: otf   _NULL   ! flux through open boundaries, allocate nobc

     real                              :: conversion            ! conversion from model specific dimensions to other useful dimensions 
     real                              :: offset                ! offset in dimensions (e.g., Celsius to Kelvin)
     real                              :: min_tracer            ! min acceptable value used for error checking 
     real                              :: max_tracer            ! max acceptable value used for error checking 
     real                              :: min_range             ! min value used for calls to diagnostic manager
     real                              :: max_range             ! max value used for calls to diagnostic manager
     real                              :: min_tracer_limit      ! min value used to limit quicker and neutral fluxes
     real                              :: max_tracer_limit      ! max value used to limit quicker and neutral fluxes 
     real                              :: min_flux_range        ! min and max values used for flux diagnostics
     real                              :: max_flux_range        ! min and max values used for flux diagnostics
     real                              :: const_init_value      ! value used to initialize tracer when using constant_init_tracer
     real                              :: scale_in              ! value by which to scale input restart data    
     real                              :: additive_in           ! uniform value added to input restart data 
     logical                           :: const_init_tracer     ! to initialize tracer to constant value 
     logical                           :: init                  ! true if the input restart file is an initial condition file
     character(len=32)                 :: flux_units            ! units for the tracer flux 
     character(len=128)                :: file_in               ! name for input restart file
     character(len=128)                :: file_out              ! name for output restart file
     character(len=32)                 :: name_in               ! name of variable to use from the input restart file
  end type ocean_prog_tracer_type

  type, public :: ocean_diag_tracer_type
     character(len=32)  :: name, units
     character(len=128) :: longname
     real, dimension(:,:,:), _ALLOCATABLE  :: field _NULL       ! tracer concentration at single time level 
     real :: factor                                  ! allows for a possibly different time stepping schemes in ocn and sea ice
                                                     ! (useful in particular for frazil)
     real :: conversion                              ! conversion from model specific dimensions to other useful dimensions 
     real :: offset                                  ! offset in dimensions (e.g., Celsius to Kelvin)
     real :: min_tracer, max_tracer                  ! min and max acceptable values used for error checking 
     real :: min_range, max_range                    ! min and max values used for diagnostics
     logical            :: init                      ! true if the input restart file is an initial condition file
     character(len=128) :: file_in                   ! name for input restart file
     character(len=128) :: file_out                  ! name for output restart file
     character(len=32)  :: name_in                   ! name of variable to use from the input restart file
     real               :: scale_in                  ! value by which to scale input restart data
     real               :: additive_in               ! uniform value added to input restart data
     real               :: const_init_value          ! value used to initialize tracer when using constant_init_tracer
     logical            :: const_init_tracer         ! to initialize tracer to constant value 
  end type ocean_diag_tracer_type



  type, public :: ocean_velocity_type
     real, _ALLOCATABLE, dimension(:,:,:,:,:)  :: u        _NULL ! horz velocity (m/sec) in i,j directions at 3 time levels
     real, _ALLOCATABLE, dimension(:,:,:)      :: smf      _NULL ! momentum flux per mass across ocean surface (m^2/sec^2)
     real, _ALLOCATABLE, dimension(:,:,:)      :: bmf      _NULL ! momentum flux per mass across ocean bottom  (m^2/sec^2)
     real, _ALLOCATABLE, dimension(:,:,:,:,:) :: advection _NULL ! thickness weighted tendency from velocity advection (m^2/sec^2) 
  end type ocean_velocity_type


  type, public :: ocean_external_mode_type
     real, _ALLOCATABLE, dimension(:,:,:)   :: eta_t          _NULL ! surface height on tracer cell center (m) 
     real, _ALLOCATABLE, dimension(:,:,:)   :: eta_u          _NULL ! surface height on tracer cell center (m)
     real, _ALLOCATABLE, dimension(:,:,:)   :: eta_t_bar      _NULL ! surface height on tracer cell time avg over ext-mode time steps (m)
     real, _ALLOCATABLE, dimension(:,:)     :: deta_dt        _NULL ! surface height time tendency on t-cell (m/s)
     real, _ALLOCATABLE, dimension(:,:,:)   :: convud_t       _NULL ! convergence of ud on T-cell center (m/s) 
     real, _ALLOCATABLE, dimension(:,:,:)   :: forcing_fs     _NULL ! depth integrated time change of velocity (without coriolis)
     real, _ALLOCATABLE, dimension(:,:,:,:) :: ud             _NULL ! vertically integrated horizontal velocity (m^2/s)
     real, _ALLOCATABLE, dimension(:,:,:)   :: grad_ps        _NULL ! rho0r * horizontal surface pressure gradient at time tau
     real, _ALLOCATABLE, dimension(:,:)     :: ps             _NULL ! surface pressure (pressure at z=0 due to eta_t) at time tau
     real, _ALLOCATABLE, dimension(:,:)     :: eta_source     _NULL ! eta sources (not including freshwater flux)
     real, _ALLOCATABLE, dimension(:,:)     :: surface_smooth _NULL ! tendency (m/s) from local smoothing of surface height 
  end type ocean_external_mode_type


#endif
!############################################################################################
! end of STATIC_MEMORY

  type, public :: ice_ocean_boundary_type
!
     real, pointer, dimension(:,:) :: au_flux      =>NULL() ! i-direction wind stress (Pa)
     real, pointer, dimension(:,:) :: av_flux      =>NULL() ! j-direction wind stress (Pa)
     real, pointer, dimension(:,:) :: at_flux      =>NULL() ! sensible heat flux (W/m2)
     real, pointer, dimension(:,:) :: aq_flux      =>NULL() ! specific humidity flux (kg/m2/s)
     real, pointer, dimension(:,:) :: alw_flux     =>NULL() ! long wave radiation (W/m2)
     real, pointer, dimension(:,:) :: asw_flux     =>NULL() ! short wave radiation (W/m2)
     real, pointer, dimension(:,:) :: ta1_ocn      =>NULL() ! air temparature (t_bot) (K)
!
     real, pointer, dimension(:,:) :: u_flux          =>NULL() ! i-direction wind stress (Pa) 
     real, pointer, dimension(:,:) :: v_flux          =>NULL() ! j-direction wind stress (Pa) 
     real, pointer, dimension(:,:) :: t_flux          =>NULL() ! sensible heat flux (W/m2) 
     real, pointer, dimension(:,:) :: q_flux          =>NULL() ! specific humidity flux (kg/m2/s)
     real, pointer, dimension(:,:) :: salt_flux       =>NULL() ! salt flux (kg/m2/s)
     real, pointer, dimension(:,:) :: lw_flux         =>NULL() ! long wave radiation (W/m2)
     real, pointer, dimension(:,:) :: sw_flux         =>NULL() ! short wave radiation (W/m2)
     real, pointer, dimension(:,:) :: sw_flux_vis     =>NULL() ! visible sw radiation (W/m2)
     real, pointer, dimension(:,:) :: sw_flux_dir     =>NULL() ! direct sw radiation (W/m2)
     real, pointer, dimension(:,:) :: sw_flux_dif     =>NULL() ! diffuse sw radiation (W/m2)
     real, pointer, dimension(:,:) :: sw_flux_vis_dir =>NULL() ! direct visible sw radiation (W/m2)
     real, pointer, dimension(:,:) :: sw_flux_vis_dif =>NULL() ! diffuse visible sw radiation (W/m2) 
     real, pointer, dimension(:,:) :: lprec           =>NULL() ! mass flux of liquid precip (kg/m2/s)
     real, pointer, dimension(:,:) :: fprec           =>NULL() ! mass flux of frozen precip (kg/m2/s)
     real, pointer, dimension(:,:) :: runoff          =>NULL() ! mass flux of liquid runoff (kg/m2/s) 
     real, pointer, dimension(:,:) :: calving         =>NULL() ! mass flux of frozen runoff (kg/m2/s) 
     real, pointer, dimension(:,:) :: p               =>NULL() ! pressure of overlying ice and atmosphere (Pa)
     integer :: xtype                                          ! REGRID, REDIST or DIRECT
  end type ice_ocean_boundary_type

  ! for communication with FMS coupler
  type, public ::  ocean_data_type 
     type(domain2d) :: Domain
     ! ESMF requires these to be pointers
     real, pointer, dimension(:,:)    :: t_surf  =>NULL()   ! SST on t-cell        (degrees Kelvin)
     real, pointer, dimension(:,:)    :: s_surf  =>NULL()   ! SSS on t-cell        (psu)
     real, pointer, dimension(:,:)    :: u_surf  =>NULL()   ! i-velocity on u-cell (m/s)
     real, pointer, dimension(:,:)    :: v_surf  =>NULL()   ! j-velocity on u-cell (m/s)
     real, pointer, dimension(:,:)    :: sea_lev =>NULL()   ! dzt(1) + eta_t + patm/rho0/grav (m) 
     real, pointer, dimension(:,:)    :: frazil  =>NULL()   ! accumulated heating (Joules/m^2) from frazil formation in ocean 

     integer, _ALLOCATABLE, dimension(:)   :: pelist  _NULL
     integer                          :: avg_kount
     logical                          :: pe
  end type ocean_data_type

public ocean_types_init
 
contains 


!#######################################################################
! <SUBROUTINE NAME="ocean_types_init">
!
! <DESCRIPTION>
! Initialize the ocean types. 
! </DESCRIPTION>
!
  subroutine ocean_types_init()

    if (module_is_initialized) then 
       call mpp_error( FATAL, '==>Error: ocean_types_init: module already initialized')
    endif
    module_is_initialized = .true.

    call write_version_number(version, tagname)

    return

  end subroutine ocean_types_init
! </SUBROUTINE> NAME="ocean_types_init"


end module ocean_types_mod


