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
! Modified: Dave Behringer
!           David.Behringer@noaa.gov
module ocean_neutral_physics_mod
! 
!<CONTACT EMAIL="stephen.griffies@noaa.gov"> Stephen M. Griffies
!</CONTACT>
!
!<CONTACT EMAIL="russell.fiedler@csiro.au"> Russell Fiedler 
!</CONTACT>
!
!<OVERVIEW>
! Thickness weighted time tendency for tracer 
! from Laplacian neutral diffusion + Laplacian GM skew-diffusion.
!</OVERVIEW>
!
!<DESCRIPTION>
! This module computes the cell thickness weighted tracer tendency from
! small angle Laplacian neutral diffusion + Laplacian GM skew-diffusion.
! It uses the "triad" algorithm also coded in MOM2.2, MOM3, and MOM3.1.
! The MOM4 algorithm accounts for partial bottom cells and generalized
! orthogonal horizontal coordinates.
!</DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! D.B. Chelton,  R.A. deSzoeke, M.G. Schlax, K.E. Naggar, N. Siwertz
! Geographical Variability of the First Baroclinic Rossby Radius of Deformation
! Journal of Physical Oceanography (1998) vol 28 pages 433-460 
! </REFERENCE>
!
! <REFERENCE>
! G. Danabasoglu and J. C. McWilliams
! Sensitivity of the global ocean circulation to 
! parameterizations of mesoscale tracer transports
! Journal of Climate (1995) vol 8 pages 2967--2987 
! </REFERENCE>
!
! <REFERENCE>
! S.M. Griffies, A. Gnanadesikan, R.C. Pacanowski, V. Larichev, 
! J.K. Dukowicz,  and R.D. Smith
! Isoneutral diffusion in a z-coordinate ocean model
! Journal of Physical Oceanography (1998) vol 28 pages 805-830
! </REFERENCE>
!
! <REFERENCE>
! S.M. Griffies 
! The Gent-McWilliams Skew-flux 
! Journal of Physical Oceanography (1998) vol 28 pages 831-841
! </REFERENCE>
!
! <REFERENCE>
! S.M. Griffies 
! Fundamentals of Ocean Climate Models (2003)
! NOAA/Geophysical Fluid Dynamics Laboratory
! </REFERENCE>
!
! <REFERENCE>
! I.M. Held and V.D. Larichev
! A scaling theory for horizontally homogeneous baroclinically 
! unstable flow on a beta plane
! Journal of Atmospheric Sciences (1996) vol 53 pages 946-952
! </REFERENCE>
!
! <REFERENCE>
! M. Visbeck, J.C. Marshall, T. Haine and M. Spall
! Specification of eddy transfer coefficients in coarse resolution ocean
! circulation models
! Journal of Physical Oceanography (1997) vol 27 pages 381--402
! </REFERENCE>
!
! <NOTE>
! This module represents perhaps the most complicated piece of code in 
! mom4. It is strongly suggested that some time be taken prior to 
! making major changes. If you find problems, either in documentation, 
! algorithm, or implementation, please feel free to contact Griffies.  
! </NOTE> 
!
! <NOTE>
! Numerical implementation follows the triad approach documented in the
! references and implemented in MOM2 and MOM3.  
! </NOTE> 
!
! <NOTE>
! Diffusivities can be flow dependent, where the length scale is related
! to the first baroclinic Rossby radius and the time scale to the Eady
! growth rate.  Within 5deg of equator, use is made of the equatorial
! Rossby radius. Length scale can also be determined by the width of 
! the baroclinic zone, as done in the Hadley Centre model. 
! If the namelist settings agm and aredi are not equal,
! then flow dependent diffusivities will be computed ONLY for the GM 
! skew-fluxes.
! </NOTE> 
!
! <NOTE>
! neutral_physics_simple=.true. requires the neutral diffusive
! diffusivity to be set equal to the GM skew-diffusive diffusivity
! (aredi=agm).  neutral_physics_simple=.true. results in down-gradient  
! horizontal flux components. This setting reduces the overall cost 
! of the neutral physics scheme significantly.  
! </NOTE> 
!
! <NOTE> 
! In steep slope regions, neutral fluxes are generally tapered to zero
! with the tanh taper of Danabasoglu and McWilliams (1995) or the 
! quadratic scheme of Gerdes, Koberle, and Willebrand.  However, if
! neutral_physics_simple=.false., the GM skew-diffusive fluxes 
! can remain nonzero if have neutral_linear_gm_taper=.true.
! </NOTE> 
!
! </INFO>
!
!<NAMELIST NAME="ocean_neutral_physics_nml">
!
!  <DATA NAME="neutral_physics_on" TYPE="logical">
!  Must be true to use this module. 
!  </DATA> 
!  <DATA NAME="neutral_physics_debug" TYPE="logical">
!  For printing starting and ending checksums for restarts
!  </DATA> 
!  <DATA NAME="diffusion_all_explicit" UNITS="dimensionless" TYPE="logical">
!  To compute all contributions from neutral diffusion explicitly in time, including
!  the K33 diagonal piece.  This approach is available only when have small time 
!  steps and/or running with just a single tracer.  It is for testing purposes. 
!  </DATA> 
!
!  <DATA NAME="neutral_physics_limit" TYPE="logical">
!  When tracer falls outside a specified range, revert to horizontal 
!  diffusive fluxes at this cell. This is an ad hoc and incomplete attempt
!  to maintain monotonicity with the neutral physics scheme.  
!  </DATA> 
!  <DATA NAME="tmask_neutral_on" TYPE="logical">
!  If .true. then this logical reduces the neutral fluxes to 
!  horizontal/vertical diffusion next to boundaries.  
!  This approach has been found to reduce spurious 
!  extrema resulting from truncation of triads used to compute 
!  a neutral flux component.   
!  </DATA> 
!
!  <DATA NAME="dm_taper" TYPE="logical">
!  Set to true to use the tanh tapering scheme of Danabasoglu and McWilliams.
!  Default is true. 
!  </DATA> 
!  <DATA NAME="gkw_taper" TYPE="logical">
!  Set to true to use the quadradic tapering scheme of Gerdes, Koberle, and Willebrand.
!  Default is false. 
!  </DATA> 
!  <DATA NAME="smax" TYPE="real">
!  Value of the maximum neutral direction slope above which the neutral fluxes are 
!  either tapered to zero or saturated.  Typical value is smax=0.01 or smaller. 
!  </DATA> 
!  <DATA NAME="swidth" TYPE="real">
!  Width in slope over which use tanh to taper fluxes in steep sloped regions. 
!  Typical value swidth=0.1*smax
!  </DATA> 
!
!  <DATA NAME="aredi" TYPE="real">
!  Neutral diffusivity used for experiments using a constant diffusivity. 
!  </DATA> 
!  <DATA NAME="agm" TYPE="real">
!  GM-skew diffusivity used for experiments using a constant diffusivity. 
!  </DATA> 
!
!  <DATA NAME="agm_lat_zones" TYPE="logical">
!  If true, will set agm_array as constant within two latitudinal zones.  
!  The idea is that one may wish to use a larger agm in the ACC than 
!  elsewhere. 
!  </DATA> 
!  <DATA NAME="agm_lat_zones_boundary" TYPE="real">
!  Boundary between agm in the south and north zones. 
!  </DATA> 
!  <DATA NAME="agm_lat_zones_ratio" TYPE="real">
!  Ratio between the large agm used in the southern latitudinal zone
!  to that used in the north.  
!  agm_array(north) = agm
!  agm_array(south) = agm*agm_lat_zones_ratio
!  </DATA> 
!
!  <DATA NAME="tracer_mix_micom" TYPE="logical">
!  If .true., then the GM-skew diffusivity is set according to a velocity scale 
!  times the grid spacing. 
!  </DATA> 
!  <DATA NAME="vel_micom" UNITS="m/sec" TYPE="real">
!  Velocity scale that is used for computing the MICOM diffusivity. 
!  </DATA> 
!
!  <DATA NAME="neutral_horz_mix_bdy" TYPE="logical">
!  If .true., then use a horizontal diffusivity in the neutral boundary layer. 
!  </DATA> 
!  <DATA NAME="vel_micom_bdy" UNITS="m/sec" TYPE="real">
!  Velocity scale that is used for computing the MICOM horizontal diffusivity 
!  within the neutral boundary layer. 
!  </DATA> 
!  <DATA NAME="ah_bdy" UNITS="m^2/sec" TYPE="real">
!  Constant horizontal diffusivity for the boundary layer.  
!  </DATA> 
!
!  <DATA NAME="bryan_lewis_aredi" TYPE="logical">
!  Set bryan_lewis_aredi=.true. when want to have aredi a function of depth
!  according to the Bryan and Lewis (1979) profile.                                
!  </DATA>
!  <DATA NAME="ahs" TYPE="real">
!  ahs = adjustable parameter at the surface for bryan_lewis_aredi   
!  </DATA>
!  <DATA NAME="ahb" TYPE="real">
!  ahb = adjustable parameter at the bottom for bryan_lewis_aredi   
!  </DATA>
!
!  <DATA NAME="agm_closure" TYPE="logical">
!  If .true. then will compute the GM-skew diffusivity as a function of the flow.
!  The length scale is determined by the Rossby radius and the time scale is 
!  determined by the Eady growth rate.  Diffusivities are depth independent.  
!  </DATA> 
!  <DATA NAME="agm_closure_max" UNITS="m^2/sec" TYPE="real">
!  Maximum GM diffusivity allowed when using agm_closure=.true. 
!  </DATA> 
!  <DATA NAME="agm_closure_min" UNITS="m^2/sec" TYPE="real">
!  Minimum GM diffusivity allowed when using agm_closure=.true. 
!  </DATA> 
!  <DATA NAME="agm_closure_scaling" UNITS="dimensionless" TYPE="logical">
!  Dimensionless tuning parameter for computing flow dependent diffusivities. 
!  </DATA> 
!  <DATA NAME="agm_closure_upper_depth" UNITS="m" TYPE="real">
!  Upper depth where start the depth integration to compute the Eady 
!  growth rate and/or baroclinicity.
!  </DATA> 
!  <DATA NAME="agm_closure_lower_depth" UNITS="m" TYPE="real">
!  Deeper depth where finish the depth integration to compute the Eady 
!  growth rate and/or baroclinicity.
!  </DATA> 
!
!  <DATA NAME="agm_closure_baroclinic" TYPE="logical">
!  For computing the agm coefficient using only the vertically
!  averaged magnitude of the horizontal density gradient 
!  (i.e., the "baroclinicity").
!  </DATA> 
!  <DATA NAME="agm_closure_buoy_freq" UNITS="sec^-1" TYPE="real">
!  For computing the agm coefficient using only the vertically
!  averaged horizontal density gradient, use a buoyancy frequency,
!  which is fixed for all space-time.
!  </DATA> 
!
!  <DATA NAME="agm_closure_length_fixed" TYPE="logical">
!  Used fixed length scale for computing agm_closure diffusivity 
!  </DATA> 
!  <DATA NAME="agm_closure_length" UNITS="meter" TYPE="real">
!  Fixed length scale for use with agm_closure_fixed_length
!  </DATA> 
!  <DATA NAME="agm_closure_length_rossby" TYPE="logical">
!  For computing the agm_closure length scale according to Rossby radius. 
!  </DATA> 
!  <DATA NAME="agm_closure_length_bczone" TYPE="logical">
!  For computing the agm_closure length scale according to radius of baroclinic zone. 
!  </DATA> 
!  <DATA NAME="bczone_max_pts" TYPE="integer">
!  Max number of horizontal grid points for use in computing the baroclinic zone radius.  
!  </DATA> 
!  <DATA NAME="agm_closure_bczone_crit_rate" UNITS="sec^-1" TYPE="real">
!  Critical growth rate for determining width of the baroclinic zone. 
!  </DATA> 
!  <DATA NAME="agm_closure_growth_scale" UNITS="dimensionless" TYPE="real">
!  Dimensionless scaling used to set a maximum for agm_growth. 
!  </DATA> 
!  <DATA NAME="rossby_radius_max" UNITS="meter" TYPE="real">
!  Maximum Rossby Radius used for agm_closure_length_rossby and 
!  the sine_taper scheme. Default = 100e3 m.
!  </DATA> 
!  <DATA NAME="rossby_radius_min" UNITS="meter" TYPE="real">
!  Minimum Rossby Radius used for agm_closure_length_rossby and 
!  the sine_taper scheme. Default = 15e3 m.
!  </DATA> 
!
!
!  <DATA NAME="neutral_physics_simple" TYPE="logical">
!  If .true. then must have aredi=agm.  The horizontal flux components are then 
!  computed as horizontal downgradient diffusive fluxes regardless the neutral slope.
!  This approach improves algorithm efficiency.  However, a .true. setting precludes
!  being able to have the GM-skew fluxes remain active in the steep sloped regions,
!  thus shutting off their effects to reduce the slopes of isopycnals in convective and 
!  mixed layer regimes.   
!  </DATA> 
!
!  <DATA NAME="neutral_linear_gm_taper" TYPE="logical">
!  If .true. then with neutral_physics_simple=.false., will linearly taper GM
!  skew fluxes towards the surface within regions of steep neutral slopes.  
!  This approach leads to a constant horizontal eddy-induced velocity in 
!  the steeply sloping regions. 
!  </DATA> 
!  <DATA NAME="neutral_taper_diagonal" UNITS="dimensionless" TYPE="logical">
!  For cases with neutral_physics_simple=.false., then neutral_taper_diagonal=.true.
!  will taper the diagonal pieces of the horizontal flux components when neutral slopes
!  are steep. With neutral_taper_diagonal=.false., then the horizontal flux components will 
!  remain enabled for all slopes, thus producing horizontal downgradient diffusion in 
!  regions of vertical neutral directions.
!  </DATA> 
!
!  <DATA NAME="neutral_sine_taper" TYPE="logical">
!  If .true. then with neutral_physics_simple=.false., will apply a sine-taper 
!  to GM and neutral diffusive fluxes in regions where the penetration depth 
!  of eddies is deeper than the grid point. This method is essential for 
!  fine vertical resolution grids. 
!  </DATA> 
!
!  <DATA NAME="gm_velocity_save" TYPE="logical">
!  If .true. then will compute the GM eddy-induce advection velocity components
!  for purposes of diagnostics.  This computation is nontrivial as it requires 
!  some extra calls to mpp_update_domains.  The algorithm used to compute the 
!  velocities may require modifications for efficiency, and suggestions are welcome.  
!  </DATA> 
!
!  <DATA NAME="neutral_blayer_diagnose" TYPE="logical">
!  Diagnose properties of the neutral physics boundary layer, whether have 
!  neutral_linear_gm_taper or neutral_sine_taper true or not.  
!  </DATA> 
!
!</NAMELIST>

use constants_mod,           only: epsln, radius, pi, radian, grav, rho0r, rho_cp, rho0
use diag_manager_mod,        only: register_diag_field, register_static_field, send_data, need_data
use fms_mod,                 only: FATAL, NOTE, file_exist, read_data, write_data
use fms_mod,                 only: open_namelist_file, check_nml_error, close_file, write_version_number
use mpp_domains_mod,         only: mpp_update_domains, mpp_define_domains, CGRID_NE, WUPDATE, SUPDATE, EUPDATE, NUPDATE
use mpp_domains_mod,         only: cyclic_global_domain, global_data_domain
use mpp_mod,                 only: mpp_pe, mpp_error, mpp_chksum, mpp_min, stdout, stdlog
use mpp_mod,                 only: mpp_clock_id, mpp_clock_begin, mpp_clock_end, CLOCK_ROUTINE
use time_manager_mod,        only: set_time, time_type, increment_time, get_date, operator ( + )

use ocean_density_mod,       only: density_derivs 
use ocean_domains_mod,       only: get_local_indices, set_ocean_domain
use ocean_operators_mod,     only: FAX, FAY, FMX, FMY, FDX_PT, FDY_PT, BDX_ET, BDY_NT
use ocean_sigma_diffuse_mod, only: tmask_sigma_on, tmask_sigma  
use ocean_types_mod,         only: ocean_grid_type, ocean_domain_type
use ocean_types_mod,         only: ocean_prog_tracer_type, ocean_thickness_type
use ocean_types_mod,         only: ocean_time_type, ocean_time_steps_type
use ocean_workspace_mod,     only: wrk1 

implicit none

public ocean_neutral_physics_init
public neutral_physics
public ocean_neutral_physics_rstrt
public ocean_neutral_physics_end

private slope_function_gm
private gm_velocity
private compute_eady_rate
private compute_baroclinicity
private compute_rossby_radius
private compute_bczone_radius
private compute_diffusivity
private fx_flux
private fy_flux
private fz_flux
private fz_terms
private neutral_slopes
private neutral_blayer

private 

type(ocean_grid_type), pointer   :: Grd => NULL()
type(ocean_domain_type), pointer :: Dom => NULL()
type(ocean_domain_type), save    :: Dom_flux
type(ocean_domain_type), save    :: BCzone_domain   

real    :: aredi  = 1.0e3  ! global constant neutral diffusion tracer diffusivity (m^2/sec)
real    :: agm    = 1.0e3  ! global constant gent-mcwilliams skew-diffusion diffusivity (m^2/sec)

real    :: agm_closure_scaling = 0.02          ! dimensionless tuning parameter for use in flow dependent diffusivity
real    :: agm_closure_max     = 2.e3          ! maximum GM diffusivity allowed when agm_closure=.true.
real    :: agm_closure_min     = 2.e2          ! minimum GM diffusivity allowed when agm_closure=.true.
real    :: agm_closure_length  = 50.e3         ! length scale (m) for use in agm_closure with agm_closure_length_fixed=.true.

logical :: agm_closure                  =.false. ! if true, then compute flow dependent diffusivity 
logical :: agm_closure_length_fixed     =.false. ! if true, then fix length scale to agm_closure_length for use with agm_closure
logical :: agm_closure_length_rossby    =.false. ! set length scale for agm_closure according to Rossby length 
logical :: agm_closure_length_bczone    =.false. ! set length scale for agm_closure according to width of baroclinic zone 
logical :: agm_closure_baroclinic       =.false. ! for using the vert ave horiz density to determine the agm_array
integer :: bczone_max_pts               =10      ! maximum number of horz grid points in a given direction for baroclinic zone search 
real    :: agm_closure_bczone_crit_rate =1.4e-6  ! critical growth rate above which we are within the "baroclinic zone" (sec^-1)
real    :: agm_closure_growth_scale     =0.5     ! dimensionless scaling used to set a maximum for agm_growth. 0.5 yields max=coriolis.
real    :: agm_closure_buoy_freq        =0.004   ! sec^-1 buoyancy frequency for use with agm_closure_baroclinic
real    :: agm_closure_upper_depth      =100.0   ! Upper depth (m) of the Eady/baroclinicity calculation
real    :: agm_closure_lower_depth      =2000.0  ! Lower depth (m) of the Eady/baroclinicity calculation

logical :: bryan_lewis_aredi = .false.      ! for Bryan-Lewis Redi diffusivity 
real    :: ahs               = 0.0          ! for Bryan-Lewis Redi diffusivity 
real    :: ahb               = 0.0          ! for Bryan-Lewis Redi diffusivity 

logical :: agm_lat_bands          = .false. ! for setting agm according to latitudinal zones 
real    :: agm_lat_bands_boundary = -999.   ! boundary between agm in the south and north zones
real    :: agm_lat_bands_ratio    = 1.0     ! ratio agm(south)/agm(north)

logical :: tracer_mix_micom     =.false. ! if true, agm and aredi set according to grid spacing
real    :: vel_micom            = 0.0    ! constant velocity scale (m/s) for setting micom diffusivity for agm and aredi

logical :: neutral_horz_mix_bdy =.false. ! if true, apply horizontal diffusivity in neutral bdy 
real    :: vel_micom_bdy        = 0.0    ! constant velocity scale (m/s) for setting horizontal micom diffusivity w/i bdy
real    :: ah_bdy               = 0.0    ! constant horiz diffusivity in neutral bdy

logical :: dm_taper         = .true.        ! set true to use tanh tapering scheme of Danabasoglu and McWilliams
logical :: gkw_taper        = .false.       ! set true to use quadratic tapering of Gerdes, Koberle, and Willebrand
real    :: smax             = 0.01          ! maximum neutral direction slope allowed before tapering fluxes
real    :: swidth           = 0.05*.01      ! slope width over which fluxes tapered using tanh function 
real    :: dm_taper_const   = 1.0           ! internally set to unity when dm_taper=.true.
real    :: gkw_taper_const  = 0.0           ! internally set to unity when gkw_taper=.true.

logical :: neutral_physics_simple=.false.   ! .true. is available only when aredi=agm 
logical :: neutral_linear_gm_taper=.true.   ! for maintaining horizontal GM velocity in bdy layer that is constant with depth
logical :: neutral_sine_taper=.true.        ! for tapering neutral physics over penetration depth of eddies
logical :: neutral_taper_diagonal=.false.   ! to remove diagonal elements to neutral diffusion with in the tapering region 
logical :: tmask_neutral_on=.false.         ! if true, reduces fluxes to horz/vert diffusion next to boundaries
logical :: diffusion_all_explicit=.false.   ! if true, compute K33 explicitly in time. Must be false for realistic cases.  
logical :: gm_velocity_save =.false.        ! for saving diagnostic gm velocity field.  Nontrivial cost.
logical :: neutral_physics_limit=.false.    ! Revert to horizontal diffusion when tracer falls outside specified range 
logical :: neutral_blayer_diagnose=.false.  ! For diagnosing properties at the base of the "neutral physics boundary layer"

integer :: num_prog_tracers = 0

! to write some output to all processors 
integer :: unit=6

! clock ids
integer :: id_clock_density_derivs
integer :: id_clock_neutral_slopes
integer :: id_clock_neutral_blayer
integer :: id_clock_fz_terms 
integer :: id_clock_fx_flux
integer :: id_clock_fy_flux
integer :: id_clock_fz_flux
integer :: id_clock_compute_eady_rate
integer :: id_clock_compute_baroclinicity
integer :: id_clock_compute_rossby_radius
integer :: id_clock_compute_bczone_radius
integer :: id_clock_compute_diffusivity 

! for diagnostic manager 
logical :: used
integer :: id_k33_explicit=-1
integer :: id_aredi=-1
integer :: id_agm=-1
integer :: id_aredi_3d=-1
integer :: id_agm_3d=-1
integer :: id_ah_bdy=-1
integer :: id_ustar=-1
integer :: id_vstar=-1
integer :: id_wstar=-1
integer :: id_rossby=-1
integer :: id_rossby_taper=-1
integer :: id_bczone=-1
integer :: id_gw_speed=-1
integer :: id_eady=-1
integer :: id_growth_rate_baroclinic=-1
integer :: id_baroclinicity=-1
integer :: id_agm_growth=-1
integer :: id_agm_length=-1
integer :: id_agm_qg=-1
integer :: id_buoy=-1
integer :: id_tx_trans_gm=-1
integer :: id_ty_trans_gm=-1
integer :: id_depth_blayer_base=-1
integer :: id_eddy_depth=-1
integer :: id_steep_depth=-1
integer :: id_slope_blayer_base=-1
integer :: id_slope31=-1
integer :: id_slope32=-1
integer, dimension(:), allocatable  :: id_neutral_physics
integer, dimension(:), allocatable  :: id_k33_implicit
integer, dimension(:), allocatable  :: id_flux_x       ! for i-directed heat flux from neutral physics 
integer, dimension(:), allocatable  :: id_flux_y       ! for j-directed heat flux from neutral physics 
integer, dimension(:), allocatable  :: id_flux_x_int_z ! for vertically integrated i-directed tracer flux 
integer, dimension(:), allocatable  :: id_flux_y_int_z ! for vertically integrated j-directed tracer flux 

real :: taper_diagonal=0.0      ! taper_diagonal=1 if neutral_taper_diagonal=.true. ; =0 otherwise 
real :: rossby_radius_max=100e3 ! maximum rossby radius (m)
real :: rossby_radius_min=15e3  ! minimum rossby radius (m)

real    :: gravrho0r
real    :: gravrho0r_buoyr
real    :: fivedeg
real    :: swidthr
real    :: smax_swidthr

real    :: dtime
real    :: two_dtime_inv

logical :: neutral_physics_on     = .false.
logical :: neutral_physics_debug  = .false.

#include <ocean_memory.h>

#ifdef STATIC_MEMORY

real, dimension(isd:ied,jsd:jed,0:1)    :: dtew, dtns         !grid distances from T-point to cell faces (m)
real, dimension(isd:ied,jsd:jed,0:1)    :: dtwedyt, dxtdtsn   !horizontal areas (m^2) of quarter cell
real, dimension(isd:ied,jsd:jed,nk,0:1) :: delqc              !quarter cell vertical thickness(m)
real, dimension(isd:ied,jsd:jed,0:nk)   :: dhwtr              !(1/dhwt)(m^-1)

real, dimension(isd:ied,jsd:jed,nk) :: aredi_array         !3D array of redi diffusivities (m^2/sec)     
real, dimension(isd:ied,jsd:jed,nk) :: agm_array           !3D array of gm diffusivities (m^2/sec)        
real, dimension(isd:ied,jsd:jed)    :: agm_micom           !2D array of micom gm diffusivities (m^2/sec)
real, dimension(isd:ied,jsd:jed)    :: ah_array            !2D array of micom horizontal diffusivities (m^2/sec)
real, dimension(isd:ied,jsd:jed)    :: rossby_radius       !first baroclinic Rossby radius (m) 
real, dimension(isd:ied,jsd:jed)    :: bczone_radius       !radius of baroclinic zone (m) 
real, dimension(isd:ied,jsd:jed)    :: eady_rate           !Eady growth rate (sec^-1)
real, dimension(isd:ied,jsd:jed)    :: eady_termx          !intermediate term for computing Eady growth rate 
real, dimension(isd:ied,jsd:jed)    :: eady_termy          !intermediate term for computing Eady growth rate 
real, dimension(isd:ied,jsd:jed)    :: baroclinicity       !vertically averaged magnitude of horiz density gradient (kg/m^4)
real, dimension(isd:ied,jsd:jed)    :: baroclinic_termx    !intermediate term for computing vert ave baroclinicity
real, dimension(isd:ied,jsd:jed)    :: baroclinic_termy    !intermediate term for computing vert ave baroclinicity 
real, dimension(isd:ied,jsd:jed)    :: coriolis_param      !absolute value of the Coriolis parameter (sec^-1)
real, dimension(isd:ied,jsd:jed)    :: beta_param          !beta = d(Coriolis)/dy (m^-1 sec^-1)
real, dimension(isd:ied,jsd:jed)    :: grid_length         !grid length scale (m)
real, dimension(isd:ied,jsd:jed)    :: count_x, count_y    !normalization factors for Eady and baroclinicity calculations

real, dimension(isd:ied,jsd:jed,nk)     :: drhodT          !drho/dtheta     (kg/m^3/C)
real, dimension(isd:ied,jsd:jed,nk)     :: drhodS          !drho/dsalinity  (kg/m^3/psu)
real, dimension(isd:ied,jsd:jed,nk,0:1) :: drhodzb         !vertical neutral density gradient (kg/m^3/m)

real, dimension(isd:ied,jsd:jed,nk,0:1,0:1)  :: tensor_31    !tracer independent portion of mixing tensor
real, dimension(isd:ied,jsd:jed,nk,0:1,0:1)  :: tensor_32    !tracer independent portion of mixing tensor
real, dimension(isd:ied,jsd:jed,nk)          :: K33_implicit !implicit in time diagonal term in redi tensor (m^2/sec) 
real, dimension(isd:ied,jsd:jed,nk)          :: K33_explicit !explicit in time diagonal term in redi tensor (m^2/sec) 

type  :: tracer_2d_type
  real, dimension(isd:ied,jsd:jed)      :: field
end type tracer_2d_type
type  :: tracer_3d_0_nk_type
  real, dimension(isd:ied,jsd:jed,0:nk) :: field
end type tracer_3d_0_nk_type
type  :: tracer_3d_1_nk_type
  real, dimension(isd:ied,jsd:jed,nk)   :: field
end type tracer_3d_1_nk_type

real, dimension(isd:ied,jsd:jed,nk)    :: ustar             !zonal Gent-McWilliams velocity at U-cell point (m/s)
real, dimension(isd:ied,jsd:jed,nk)    :: vstar             !meridional Gent-McWilliams velocity at U-cell point (m/s)
real, dimension(isd:ied,jsd:jed,0:nk)  :: wstar             !vertical Gent-McWilliams velocity at T-cell bottom (m/s)
real, dimension(isd:ied,jsd:jed,nk)    :: tx_trans_gm       !for diagnosing i-transport due to GM (Sv)
real, dimension(isd:ied,jsd:jed,nk)    :: ty_trans_gm       !for diagnosing j-transport due to GM (Sv)

real, dimension(isd:ied,jsd:jed)       :: depth_blayer_base ! depth (m) of boundary layer base for neutral physics
real, dimension(isd:ied,jsd:jed)       :: slope_blayer_base ! abs(slope) at base of neutral boundary layer 
real, dimension(isd:ied,jsd:jed)       :: eddy_depth        ! max of depth (m) mesoscale eddies penetrate & kpp bdy layer

#else

real, dimension(:,:,:),     allocatable :: dtew, dtns         !grid distances from T-point to cell faces (m)
real, dimension(:,:,:),     allocatable :: dtwedyt, dxtdtsn   !horizontal areas (m^2) of quarter cell
real, dimension(:,:,:,:),   allocatable :: delqc              !quarter cell vertical thickness(m)
real, dimension(:,:,:),     allocatable :: dhwtr              !(1/dhwt)(m^-1)

real, dimension(:,:,:),     allocatable :: aredi_array        !3D array of redi diffusivities (m^2/sec) 
real, dimension(:,:,:),     allocatable :: agm_array          !3D array of gm diffusivities (m^2/sec)     
real, dimension(:,:),       allocatable :: agm_micom          !2D array of micom gm diffusivities (m^2/sec)
real, dimension(:,:),       allocatable :: ah_array           !2D array of micom horizontal diffusivities (m^2/sec)
real, dimension(:,:),       allocatable :: rossby_radius      !first baroclinic Rossby radius (m) 
real, dimension(:,:),       allocatable :: bczone_radius      !radius of baroclinic zone (m) 
real, dimension(:,:),       allocatable :: eady_rate          !Eady growth rate (sec^-1)
real, dimension(:,:),       allocatable :: eady_termx         !intermediate term for computing Eady growth rate 
real, dimension(:,:),       allocatable :: eady_termy         !intermediate term for computing Eady growth rate 
real, dimension(:,:),       allocatable :: baroclinicity      !vertically averaged magnitude of horiz density gradient (kg/m^4)
real, dimension(:,:),       allocatable :: baroclinic_termx   !intermediate term for computing vert ave baroclinicity 
real, dimension(:,:),       allocatable :: baroclinic_termy   !intermediate term for computing vert ave baroclinicity 
real, dimension(:,:),       allocatable :: coriolis_param     !absolute value of the Coriolis parameter (sec^-1)
real, dimension(:,:),       allocatable :: beta_param         !beta = d(Coriolis)/dy (m^-1 sec^-1)
real, dimension(:,:),       allocatable :: grid_length        !grid length scale (m)
real, dimension(:,:),       allocatable :: count_x, count_y   !normalization factors for Eady and baroclinicity calculations

real, dimension(:,:,:),     allocatable :: drhodT             !drho/dtheta     (kg/m^3/C)
real, dimension(:,:,:),     allocatable :: drhodS             !drho/dsalinity  (kg/m^3/psu)
real, dimension(:,:,:,:),   allocatable :: drhodzb            !vertical neutral density gradient (kg/m^3/m)

real, dimension(:,:,:,:,:), allocatable :: tensor_31      !tracer independent portion of mixing tensor
real, dimension(:,:,:,:,:), allocatable :: tensor_32      !tracer independent portion of mixing tensor
real, dimension(:,:,:),     allocatable :: K33_implicit   !implicit in time diagonal term in redi tensor (m^2/sec) 
real, dimension(:,:,:),     allocatable :: K33_explicit   !explicit in time diagonal term in redi tensor (m^2/sec) 

type  :: tracer_2d_type
  real, dimension(:,:), pointer     :: field => NULL()
end type tracer_2d_type
type  :: tracer_3d_0_nk_type
  real, dimension(:,:,:), pointer   :: field => NULL()
end type tracer_3d_0_nk_type
type  :: tracer_3d_1_nk_type
  real, dimension(:,:,:), pointer   :: field => NULL()
end type tracer_3d_1_nk_type

real, dimension(:,:,:),  allocatable :: ustar             !zonal Gent-McWilliams velocity at U-cell point (m/s)
real, dimension(:,:,:),  allocatable :: vstar             !meridional Gent-McWilliams velocity at U-cell point (m/s)
real, dimension(:,:,:),  allocatable :: wstar             !vertical Gent-McWilliams velocity at T-cell bottom (m/s)
real, dimension(:,:,:),  allocatable :: tx_trans_gm       !for diagnosing i-transport due to GM (Sv)
real, dimension(:,:,:),  allocatable :: ty_trans_gm       !for diagnosing j-transport due to GM (Sv)
real, dimension(:,:),    allocatable :: depth_blayer_base ! depth (m) of boundary layer base for neutral physics
real, dimension(:,:),    allocatable :: slope_blayer_base ! abs(slope) at base of neutral boundary layer 
real, dimension(:,:),    allocatable :: eddy_depth        ! max of depth (m) mesoscale eddies penetrate & kpp bdy layer


#endif

! for computing radius of baroclinic zone for nonconstant diffusivity 
real, dimension(:,:), allocatable :: bczone_dxt    ! for determining baroclinic zone radius (m) 
real, dimension(:,:), allocatable :: bczone_dyt    ! for determining baroclinic zone radius (m) 
real, dimension(:,:), allocatable :: bczone_rate   ! for determining baroclinic zone radius (sec^-1)


! introduce following derived types so that do not need to know num_prog_tracers at compile time 

type(tracer_3d_1_nk_type), dimension(:), allocatable  :: dTdx    ! tracer partial derivative (tracer/m)
type(tracer_3d_1_nk_type), dimension(:), allocatable  :: dTdy    ! tracer partial derivative (tracer/m)
type(tracer_3d_0_nk_type), dimension(:), allocatable  :: dTdz    ! tracer partial derivative (tracer/m)

type(tracer_2d_type), dimension(:), allocatable       :: fz1     ! z-flux component for all the tracers at particular k 
type(tracer_2d_type), dimension(:), allocatable       :: fz2     ! z-flux component for all the tracers at particular k 

type(tracer_3d_1_nk_type), dimension(:), allocatable  :: flux_x  ! i-flux component for all the tracers for all k
type(tracer_3d_1_nk_type), dimension(:), allocatable  :: flux_y  ! j-flux component for all the tracers for all k


integer :: index_temp
integer :: index_salt

character(len=128) :: version=&
     '$Id$'
character (len=128) :: tagname = &
     '$Name$'

logical :: module_is_initialized = .FALSE.

namelist /ocean_neutral_physics_nml/ neutral_physics_on, aredi, agm, &
                                     tracer_mix_micom, vel_micom, &
                                     neutral_horz_mix_bdy, vel_micom_bdy, ah_bdy, &
                                     agm_closure, agm_closure_scaling, agm_closure_max, agm_closure_min, &
                                     agm_closure_growth_scale, agm_closure_length_fixed, agm_closure_length, &
                                     agm_closure_length_rossby, agm_closure_length_bczone, &
                                     bczone_max_pts, agm_closure_bczone_crit_rate, &
                                     agm_closure_baroclinic, agm_closure_buoy_freq, & 
                                     agm_closure_upper_depth, agm_closure_lower_depth, &  
                                     dm_taper, gkw_taper, smax, swidth, tmask_neutral_on, &
                                     neutral_physics_simple, neutral_taper_diagonal, &
                                     neutral_linear_gm_taper, neutral_sine_taper,& 
                                     gm_velocity_save, diffusion_all_explicit, &   
                                     bryan_lewis_aredi, ahs, ahb, &
                                     agm_lat_bands, agm_lat_bands_boundary, agm_lat_bands_ratio, & 
                                     neutral_physics_limit, neutral_physics_debug, &
                                     rossby_radius_max, rossby_radius_min, neutral_blayer_diagnose
    

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_neutral_physics_init">
!
! <DESCRIPTION>
! Initialize the neutral physics module by registering fields for 
! diagnostic output and performing some numerical checks to see 
! that namelist settings are appropriate. 
! </DESCRIPTION>
!
subroutine ocean_neutral_physics_init(Grid, Domain, Time, Time_steps, Thickness, T_prog, ens_ocean)

  type(ocean_grid_type), intent(in), target   :: Grid
  type(ocean_domain_type), intent(in), target :: Domain
  type(ocean_time_type), intent(in)           :: Time
  type(ocean_time_steps_type), intent(in)     :: Time_steps
  type(ocean_thickness_type), intent(in)      :: Thickness
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)
  logical, intent(in), optional               :: ens_ocean

  integer :: ioun, io_status, i, j, k, n, num, num_schemes 
  integer :: i_delta, j_delta, k_delta, ierr
  real    :: ft, deltaX, deltaY, delta, delta_iso, delta_iso0, delta_min, A_max, A_max0
  real    :: deg2m
  
  if ( module_is_initialized ) then 
    call mpp_error(FATAL, '==>Error from ocean_neutral_physics_mod (ocean_neutral_physics_init): module already initialized')
  endif 

  module_is_initialized = .TRUE.

  call write_version_number( version, tagname )

  num_prog_tracers = size(T_prog(:))
  dtime            = Time_steps%dtime_t

  write(stdout(),'(/a,f10.2)')'==> Note from ocean_neutral_physics_mod: using forward time step of (secs)', dtime 

  ! provide for namelist over-ride of default values
  ioun =  open_namelist_file()
  read (ioun,ocean_neutral_physics_nml,IOSTAT=io_status)
  write (stdout(),'(/)')
  write (stdout(),ocean_neutral_physics_nml)  
  write (stdlog(),ocean_neutral_physics_nml)
  ierr = check_nml_error(io_status,'ocean_neutral_physics_nml')

  call close_file (ioun)

#ifndef STATIC_MEMORY
  call get_local_indices(Domain,isd,ied,jsd,jed,isc,iec,jsc,jec)
  nk = Grid%nk
#endif

  if(neutral_physics_on) then 
    call mpp_error(NOTE, '==> Note from ocean_neutral_physics_mod: USING ocean_neutral_physics_mod.')
    write(stdout(),*)'==> Note from ocean_neutral_physics_mod: USING ocean_neutral_physics.'
  else 
    call mpp_error(NOTE, '==> Note from ocean_neutral_physics_mod: NOT using ocean_neutral_physics_mod.')
    write(stdout(),*)'==> Note from ocean_neutral_physics_mod: NOT using ocean_neutral_physics.'
    return
  endif 

  allocate( dTdx(num_prog_tracers) )
  allocate( dTdy(num_prog_tracers) )
  allocate( dTdz(num_prog_tracers) )
  allocate( fz1(num_prog_tracers) )
  allocate( fz2(num_prog_tracers) )
  allocate( flux_x(num_prog_tracers) )
  allocate( flux_y(num_prog_tracers) )
     
  if(Time_steps%aidif /= 1.0 .and. aredi /= 0) then 
      call mpp_error(FATAL, &
           '==>Error from ocean_neutral_physics_mod: stability requires aidif=1.0 since K33 must be handled implicitly.')
  endif
  if(neutral_physics_simple .and. aredi/=agm) then
      write(stdout(),'(a,e10.4,a,e10.4)')'In namelist, aredi = ', aredi,' and agm = ',agm
      call mpp_error(FATAL, &
           'Error from ocean_neutral_physics_mod: with neutral_physics_simple=.true., must set aredi=agm.')
  endif
  if(neutral_physics_simple .and. aredi==agm) then
      write(stdout(),'(1x,a)')' ==> Note from ocean_neutral_physics_mod: Running with neutral_physics_simple=.true.  ' 
      write(stdout(),'(7x,a)')'Horizontal SGS tracer fluxes are downgradient for all neutral slopes. ' 
      write(stdout(),'(7x,a)')'Method is analogous to that used in MOM3.' 
  endif

  if(neutral_physics_limit) then
      write(stdout(),'(1x,a)')' ==>Note from ocean_neutral_physics_mod: neutral_physics_limit=.true.' 
      write(stdout(),'(7x,a)')'Will revert to horizontal diffusion for points where tracer is outside specified range.' 
  endif

  if(dm_taper .and. .not. gkw_taper) then
      write(stdout(),'(1x,a)')' ==>Note from ocean_neutral_physics_mod: dm_taper=.true. Will use the tanh scheme ' 
      write(stdout(),'(7x,a)')'of Danabasoglu and McWilliams to taper neutral physics in steep sloped regions'
      dm_taper_const =1.0
      gkw_taper_const=0.0
  endif
  if(gkw_taper .and. .not. dm_taper) then
      write(stdout(),'(1x,a)')' ==>Note from ocean_neutral_physics_mod: gkw_taper=.true. Will use the quadratic scheme' 
      write(stdout(),'(7x,a)')'of Gerdes, Koberle, and Willebrand to taper neutral physics in steep sloped regions'
      dm_taper_const =0.0
      gkw_taper_const=1.0
  endif
  if(gkw_taper .and. dm_taper) then
      dm_taper_const =0.0
      gkw_taper_const=0.0
      call mpp_error(FATAL, &
           '==>Error from ocean_neutral_physics_mod: gkw_taper and dm_taper cannot both be set true--choose only one.')
  endif

  if(.not. neutral_physics_simple) then 
      if(neutral_linear_gm_taper) then
          write(stdout(),'(1x,a)')' ==> Note from ocean_neutral_physics_mod: neutral_linear_gm_taper=.true., so will linearly'
          write(stdout(),'(7x,a)')'taper GM towards the surface when reaching steep neutral slopes in surface bdy.'  
      else
          call mpp_error(NOTE, &
               ' ==> Note from ocean_neutral_physics_mod: GM exponentially tapered to zero in neutral bdy layer region.')
      endif
      if(neutral_sine_taper) then
          write(stdout(),'(1x,a)')' ==> Note: Running with neutral_sine_taper=.true., and so will use sine-taper'
          write(stdout(),'(7x,a)')'on fluxes where eddy penetration depth and/or KPP hblt exceeds grid depth.'
      endif
  endif

  if(neutral_blayer_diagnose) then 
      write(stdout(),'(1x,a)')' ==> Note: Running with neutral_blayer_diagnose=.true., so will diagnose properties'
      write(stdout(),'(7x,a)')'of the neutral physics boundary layer.'
  endif

  if(neutral_linear_gm_taper) then
      write(stdout(),'(1x,a)')' ==> WARNING: Running with a nontrivial GM transport in steep neutral slope regions.'
      write(stdout(),'(7x,a)')'Doing so has been found to require the use of convective_adjustment=.true.'
      write(stdout(),'(7x,a)')'Otherwise, GM may be exposed to gravitationally unstable parcels, and will'
      write(stdout(),'(7x,a)')'then act to further destabilize these parcels as it reduces baroclinicity.'
  endif

  if(diffusion_all_explicit .and. aredi > 0.0) then
      write(stdout(),'(1x,a)')' ==> Warning: Running w/ diffusion_all_explicit=.true., which means compute K33 contribution'
      write(stdout(),'(7x,a)')'to neutral diffusion explicitly in time.  This method is stable only if taking'
      write(stdout(),'(7x,a)')'very small time steps and/or running with just a single tracer.'
  endif

  if(neutral_taper_diagonal) then
      taper_diagonal = 1.0 
      write(stdout(),'(1x,a,f5.2)')' ==>Note: Running w/ neutral_taper_diagonal=.true. and so taper_diagonal = ',taper_diagonal 
      write(stdout(),'(7x,a)')'Will taper all pieces of horiz neutral flux components in steep neutral slope regions.'
      write(stdout(),'(7x,a)')'neutral_taper_diagonal=.false. will alternatively keep the diagonal pieces untapered.'
  endif
  if(.not. neutral_taper_diagonal) then
      taper_diagonal = 0.0 
      write(stdout(),'(1x,a,f5.2)')' ==> Note: Running w/ neutral_taper_diagonal=.false. so taper_diagonal = ',taper_diagonal 
      write(stdout(),'(7x,a)')'Diagonal pieces of horizontal flux components are untapered regardless the neutral slope.'
  endif

  if(agm /= 0.0 .and. gm_velocity_save) then 
      write(stdout(),'(1x,a)')' ==> Note: Computing Gent-McWilliams velocity for diagnostic purposes.'
      write(stdout(),'(7x,a)')'This computation requires extra time and memory.'
      write(stdout(),'(7x,a)')'Algorithm used for this diagnostic calculation is rudimentary.'
      write(stdout(),'(7x,a)')'Set gm_velocity_save=.false. if do not wish to save GM velocities.'
  endif
  if(agm == 0.0 .and. gm_velocity_save) then 
      call mpp_error(FATAL, '==>Error: Computing Gent-McWilliams velocity for diagnostic purposes yet agm=0')
  endif
  if(aredi == 0.0) then 
      write(stdout(),'(1x,a)')' ==> Note: aredi=0, which means will not include neutral diffusion.'
  endif
  if(aredi==0.0 .and. agm==0) then 
      write(stdout(),'(a)')'=>Warning: using ocean_neutral_physics.F90 with agm=0 and aredi=0.'
      write(stdout(),'(a)')'           Set one or both to nonzero to get nontrivial effects.'
  endif

  if(agm_lat_bands) then
      write(stdout(),'(a)')'=>Note: Setting agm_array according to agm_lat_bands.'
      write(stdout(),'(a,e10.5)')'      The ratio agm(south)/agm(north) is given by ',agm_lat_bands_boundary
      write(stdout(),'(a,e10.5)')'      with the latitude separating the bands given by ',agm_lat_bands_ratio
      if(agm_lat_bands_boundary <= -90.0) then 
          write(stdout(),'(a)')    '      Since agm_lat_bands_boundary <= -90, will default to globally constant agm value' 
      endif
  endif

  if(agm_closure) then

      num_schemes = 0 
      if(agm_closure_length_fixed) then 
          num_schemes = num_schemes + 1
          write(stdout(),'(a)')'=>Note: Computing local flow-dependent tracer diffusivity with agm_closure_length_fixed.'
      endif
      if(agm_closure_length_rossby) then 
          num_schemes = num_schemes + 1
          write(stdout(),'(a)')'=>Note: Computing local flow-dependent tracer diffusivity with agm_closure_length_rossby.'
      endif
      if(agm_closure_length_bczone) then 
          num_schemes = num_schemes + 1
          write(stdout(),'(a)')'=>Note: Computing local flow-dependent tracer diffusivity with agm_closure_length_bczone.'
      endif
      if(agm_closure_baroclinic) then 
          num_schemes = num_schemes + 1
          write(stdout(),'(a)')'=>Note: Computing local flow-dependent tracer diffusivity with agm_closure_baroclinic.'
      endif

      if(num_schemes == 0) then 
          call mpp_error(FATAL, '==>Error: with agm_closure=.true., must choose one of the "agm_closure_length" methods')
      endif
      if(num_schemes > 1) then 
          call mpp_error(FATAL, '==>Error: with agm_closure=.true., can choose only one of the available agm_closure methods')
      endif

      write(stdout(),'(a,e10.5)')'        The maximum allowable diffusivity (m^2/s) is given by ',agm_closure_max
      write(stdout(),'(a,e10.5)')'        The minimum allowable diffusivity (m^2/s) is given by ',agm_closure_min
      if(aredi == agm) then 
          write(stdout(),'(a)')'        Since aredi==agm, use same flow dependent diffusivity for diffusion & skew-diffusion.'
      endif
      if(aredi /= agm) then 
          write(stdout(),'(a)')'        Since aredi/=agm, flow dependent diffusivity used ONLY for GM skew-diffusion.'
          write(stdout(),'(a,e10.5)')'       Neutral diffusivity will be constant and given by namelist aredi(m^2/sec) ', aredi
      endif

      if(agm_closure_baroclinic) then 
          write(stdout(),'(a)')'        Length and time scales set by vertically averaged baroclinicity |grad(rho)|, '
          write(stdout(),'(a,e10.5)')'        as well as the constant buoyancy freq(sec^-1) = ',agm_closure_buoy_freq
          write(stdout(),'(a,e10.5)')'        and the constant length scale (m) = ',agm_closure_length
      else 
          write(stdout(),'(a)')'        Eady growth rate gives inverse time scale.'
          if(agm_closure_length_fixed) then 
              write(stdout(),'(a,e10.5,a)')'        Length scale set by nml parameter agm_closure_length =',agm_closure_length,' m'
          elseif(agm_closure_length_rossby) then 
              write(stdout(),'(a)')'        First baroclinic Rossby radius gives the length scale.'
          elseif(agm_closure_length_bczone) then 
              write(stdout(),'(a)')'        Radius of baroclinic zone gives the length scale.'
          endif
      endif

  endif   ! endif for agm_closure 


  ! Useful constants 
  gravrho0r       = grav*rho0r                       !for buoyancy frequency calculation  
  gravrho0r_buoyr = gravrho0r/agm_closure_buoy_freq  !for agm_closure_baroclinic
  two_dtime_inv   = 0.5/dtime                        !for explicit piece of K33 
  fivedeg         = 5.0*pi/180.0                     !for equatorial Rossby radius of deformation
  swidthr         = 1.0/(swidth + epsln)             !for slope taper function when dm_taper used 
  smax_swidthr    = smax*swidthr                     !for slope taper function when dm_taper used 

  Dom => Domain
  Grd => Grid

  do n=1,num_prog_tracers
    T_prog(n)%neutral_physics_limit = neutral_physics_limit
  enddo 

  call set_ocean_domain(Dom_flux, Grid, xhalo=Dom%xhalo, yhalo = Dom%yhalo,name='flux dom neutral')

#ifndef STATIC_MEMORY

  allocate (dtew(isd:ied,jsd:jed,0:1))
  allocate (dtns(isd:ied,jsd:jed,0:1))
  allocate (dtwedyt(isd:ied,jsd:jed,0:1))
  allocate (dxtdtsn(isd:ied,jsd:jed,0:1))
  allocate (grid_length(isd:ied,jsd:jed))
  allocate (delqc(isd:ied,jsd:jed,nk,0:1))
  allocate (dhwtr(isd:ied,jsd:jed,0:nk))
  allocate (aredi_array(isd:ied,jsd:jed,nk)) 
  allocate (agm_array(isd:ied,jsd:jed,nk))   
  allocate (agm_micom(isd:ied,jsd:jed))
  allocate (ah_array(isd:ied,jsd:jed))
  allocate (rossby_radius(isd:ied,jsd:jed))
  allocate (bczone_radius(isd:ied,jsd:jed))
  allocate (eady_rate(isd:ied,jsd:jed))
  allocate (eady_termx(isd:ied,jsd:jed))
  allocate (eady_termy(isd:ied,jsd:jed))
  allocate (baroclinicity(isd:ied,jsd:jed))
  allocate (baroclinic_termx(isd:ied,jsd:jed))
  allocate (baroclinic_termy(isd:ied,jsd:jed))
  allocate (count_x(isd:ied,jsd:jed))
  allocate (count_y(isd:ied,jsd:jed))
  allocate(tx_trans_gm(isd:ied,jsd:jed,nk))
  allocate(ty_trans_gm(isd:ied,jsd:jed,nk))

  allocate (coriolis_param(isd:ied,jsd:jed))
  allocate (beta_param(isd:ied,jsd:jed))
  allocate (drhodT(isd:ied,jsd:jed,nk))
  allocate (drhodS(isd:ied,jsd:jed,nk))
  allocate (drhodzb(isd:ied,jsd:jed,nk,0:1))
  allocate (tensor_31(isd:ied,jsd:jed,nk,0:1,0:1))
  allocate (tensor_32(isd:ied,jsd:jed,nk,0:1,0:1))
  allocate (K33_implicit(isd:ied,jsd:jed,nk)) 
  allocate (K33_explicit(isd:ied,jsd:jed,nk)) 

  if(gm_velocity_save) then 
    allocate (ustar(isd:ied,jsd:jed,nk))
    allocate (vstar(isd:ied,jsd:jed,nk))
    allocate (wstar(isd:ied,jsd:jed,0:nk))
  else
    allocate (ustar(1,1,1))
    allocate (vstar(1,1,1))
    allocate (wstar(1,1,1))
  endif 

  do n = 1, num_prog_tracers
    allocate ( dTdx(n)%field(isd:ied,jsd:jed,nk) )
    allocate ( dTdy(n)%field(isd:ied,jsd:jed,nk) )
    allocate ( dTdz(n)%field(isd:ied,jsd:jed,0:nk) )
    allocate ( fz1(n)%field(isd:ied,jsd:jed) )
    allocate ( fz2(n)%field(isd:ied,jsd:jed) )
    allocate ( flux_x(n)%field(isd:ied,jsd:jed,nk) )
    allocate ( flux_y(n)%field(isd:ied,jsd:jed,nk) )
  enddo 

  allocate(depth_blayer_base(isd:ied,jsd:jed))
  allocate(slope_blayer_base(isd:ied,jsd:jed))
  allocate(eddy_depth(isd:ied,jsd:jed))

#endif

  depth_blayer_base = 0.0
  slope_blayer_base = 0.0
  eddy_depth        = 0.0
  tx_trans_gm       = 0.0
  ty_trans_gm       = 0.0

  ! for computing the baroclinic zone radius using Hadley Centre search algorithm 
  call mpp_define_domains((/1,Grd%ni,1,Grd%nj/),Domain%layout,BCzone_domain%domain2d, &
         xflags = CYCLIC_GLOBAL_DOMAIN, xhalo=bczone_max_pts, yhalo=bczone_max_pts, name='bczone')
  allocate (bczone_rate(isc-bczone_max_pts:iec+bczone_max_pts,jsc-bczone_max_pts:jec+bczone_max_pts))
  allocate (bczone_dxt (isc-bczone_max_pts:iec+bczone_max_pts,jsc-bczone_max_pts:jec+bczone_max_pts))
  allocate (bczone_dyt (isc-bczone_max_pts:iec+bczone_max_pts,jsc-bczone_max_pts:jec+bczone_max_pts))
  bczone_dxt=0.0 ; bczone_dyt=0.0 ; bczone_rate=0.0
  bczone_dxt(isc:iec,jsc:jec) = Grd%dxt(isc:iec,jsc:jec)
  bczone_dyt(isc:iec,jsc:jec) = Grd%dyt(isc:iec,jsc:jec)
  call mpp_update_domains (bczone_dxt(:,:), BCzone_domain%domain2d)
  call mpp_update_domains (bczone_dyt(:,:), BCzone_domain%domain2d)


  dtew(:,:,0) = Grd%dtw(:,:)
  dtew(:,:,1) = Grd%dte(:,:)
  dtns(:,:,0) = Grd%dts(:,:)
  dtns(:,:,1) = Grd%dtn(:,:)

  dtwedyt(:,:,:) = 0.0
  dtwedyt(:,:,0) = Grd%dte(:,:)*Grd%dyt(:,:)
  do i=isc-1,iec
    dtwedyt(i,:,1) = Grd%dtw(i+1,:)*Grd%dyt(i+1,:)
  enddo

  dxtdtsn(:,:,:) = 0.0
  dxtdtsn(:,:,0) = Grd%dxt(:,:)*Grd%dtn(:,:)
  do j=jsc-1,jec
    dxtdtsn(:,j,1) = Grd%dxt(:,j+1)*Grd%dts(:,j+1)
  enddo

  grid_length(:,:) = 2.0*Grd%dxt(:,:)*Grd%dyt(:,:)/(Grd%dxt(:,:)+Grd%dyt(:,:))

  aredi_array(:,:,:)       = aredi
  agm_array(:,:,:)         = agm
  agm_micom(:,:)           = vel_micom*grid_length(:,:)
  ah_array(:,:)            = 0.0
  rossby_radius(:,:)       = 0.0
  bczone_radius(:,:)       = 0.0
  eady_rate(:,:)           = 0.0
  eady_termx(:,:)          = 0.0
  eady_termy(:,:)          = 0.0
  baroclinicity(:,:)       = 0.0
  baroclinic_termx(:,:)    = 0.0
  baroclinic_termy(:,:)    = 0.0
  count_x(:,:)             = 0.0
  count_y(:,:)             = 0.0

  do i=isc,iec
    do j=jsc,jec
      count_x(i,j) = min(Grd%tmask(i+1,j,1) + Grd%tmask(i-1,j,1), 1.0/(2.0*(Grd%tmask(i+1,j,1) + Grd%tmask(i-1,j,1) + epsln)))
      count_y(i,j) = min(Grd%tmask(i,j+1,1) + Grd%tmask(i,j-1,1), 1.0/(2.0*(Grd%tmask(i,j+1,1) + Grd%tmask(i,j-1,1) + epsln)))
    enddo
  enddo

  coriolis_param(:,:)  = 2.0*7.292e-5*abs(sin(Grd%phit(:,:)))
  beta_param(:,:)      = 2.28e-11*abs(cos(Grd%phit(:,:)))
  drhodT(:,:,:)        = 0.0
  drhodS(:,:,:)        = 0.0
  drhodzb(:,:,:,:)     = 0.0
  tensor_31(:,:,:,:,:) = 0.0
  tensor_32(:,:,:,:,:) = 0.0
  K33_implicit(:,:,:)  = 0.0           
  K33_explicit(:,:,:)  = 0.0           

  ustar(:,:,:) = 0.0
  vstar(:,:,:) = 0.0
  wstar(:,:,:) = 0.0

  do n = 1, num_prog_tracers 
    dTdx(n)%field(:,:,:)   = 0.0
    dTdy(n)%field(:,:,:)   = 0.0
    dTdz(n)%field(:,:,:)   = 0.0
    fz1(n)%field(:,:)      = 0.0
    fz2(n)%field(:,:)      = 0.0
    flux_x(n)%field(:,:,:) = 0.0
    flux_y(n)%field(:,:,:) = 0.0
  enddo  



  ! put Bryan-Lewis profile into aredi_array
  ! this is ad hoc and implemented only for legacy 
  if(bryan_lewis_aredi .and. neutral_physics_on) then
    write(stdout(),'(/a)')' Using Bryan-Lewis depth profile for the Redi diffusivity.'
    write(stdout(),'(a)')' The Bryan-Lewis profile is ad hoc and it has been used ONLY with agm==0.'
    write(stdout(),'(a)')' Vertical dependence to the neutral diffusivity has not been implemented'
    write(stdout(),'(a)')' according to a physical theory giving a vertical structure depending'
    write(stdout(),'(a)')' on the flow field. More theoretical work is required.'
    do k=1,nk
      aredi_array(:,:,k) = (ahb + (ahs - ahb)*exp(-Grd%zt(k)*100./50000.0))
    enddo
    do k=1,nk
      write (stdout(),'(a,i3,e16.8)') '  k, diffus = ', k, aredi_array(isc,jsc,k)
    enddo
  endif

  ! Set diffusivity according to agm_lat_bands
  if(agm_lat_bands .and. neutral_physics_on) then
      agm_array(:,:,:) = agm
      if(agm_lat_bands_boundary > -90.) then 
          do j=jsc-Dom%yhalo,jec+Dom%yhalo
             do i=isc-Dom%xhalo,iec+Dom%xhalo
                if(Grd%yt(i,j) <= agm_lat_bands_boundary) agm_array(i,j,:) = agm*agm_lat_bands_ratio
             enddo
          enddo
      endif
      if(agm==aredi) aredi_array(:,:,:) = agm_array(:,:,:)
  endif

  ! grid-scale dependent diffusivity suggested by that commonly used in MICOM
  ! vel_micom (m/s) sets velocity scale.
  ! space scale is set by grid size. 
  if(tracer_mix_micom .and. neutral_physics_on) then
      do k=1,nk
        agm_array(:,:,k) = agm_micom(:,:)
      enddo
      do j=jsc,jec
         write (stdout(),'(a,i4,a,e14.7,a)') ' Micom agm_array at (isc,',j,',1) = ',agm_array(isc,j,1),' m^2/s'
      enddo
      if(aredi == agm) then
          aredi_array(:,:,:) = agm_array(:,:,:) 
          write (stdout(),'(a)') ' aredi_array = agm_array as given by Micom grid scale dependent diffusivity'
      else
          write (stdout(),'(a)') ' Since (agm/=aredi), aredi_array=aredi'
      endif
  endif

  ! grid-scale dependent diffusivity suggested by that commonly used in MICOM
  ! vel_micom (m/s) sets velocity scale.
  ! space scale is set by grid size. 
  if(neutral_horz_mix_bdy .and. neutral_physics_on) then
      write (stdout(),'(a)') ' Adding horizontal diffusivity in neutral boundary.'
      write (stdout(),'(a)') ' This method is implemented only for the case with neutral_physics_simple=.false.'
      if(vel_micom_bdy > 0.0) then
         ah_array(:,:) = vel_micom_bdy*grid_length(:,:)
      else 
         ah_array(:,:) = ah_bdy
      endif 
      do j=jsc,jec
         write (stdout(),'(a,i4,a,e14.7,a)') ' ah_array in neutral bdy layer at (isc,',j,',1)= ',ah_array(isc,j),' m^2/s'
      enddo
  else 
      ah_array=0.0 
  endif


  if(neutral_physics_on) then 

      if(agm_closure) then 
          if(.NOT.file_exist('INPUT/ocean_neutral.res.nc')) then
              if (.NOT. Time%init) call mpp_error(FATAL,'Expecting file INPUT/ocean_neutral.res.nc to exist.&
          &This file was not found and Time%init=.false.')
              if(tracer_mix_micom) then 
                write (stdout(),'(a)') 'non-constant agm_array being initialized to MICOM values.'
                do k=1,nk   
                   agm_array(:,:,k) = agm_micom(:,:)
                enddo
              else 
                write (stdout(),'(a)') 'non-constant agm_array being initialized to constant value.'
                agm_array(:,:,:) = agm 
              endif 
          else
              call read_data('INPUT/ocean_neutral.res','agm_array',agm_array,domain=&
                              Dom%domain2d, timelevel=1, append_pelist_name=ens_ocean)
              call mpp_update_domains(agm_array,Dom%domain2d)
              write (stdout(),'(a)') 'non-constant agm_array being read from restart.'
              if(neutral_physics_debug) then 
                  write (stdout(), *) &
                  'checksum start agm_array', mpp_chksum(agm_array(isc:iec,jsc:jec,:)*Grd%tmask(isc:iec,jsc:jec,:))
              endif
          endif
      endif
      if(aredi==agm) aredi_array(:,:,:) = agm_array(:,:,:)

      if((agm_closure .and. agm_closure_length_rossby) .or. neutral_sine_taper) then 
          if(file_exist('INPUT/ocean_neutral.res.nc')) then 
              call read_data('INPUT/ocean_neutral.res','rossby_radius',rossby_radius,domain=&
                              Dom%domain2d, timelevel=1, append_pelist_name=ens_ocean)
              call mpp_update_domains(rossby_radius,Dom%domain2d)
              write (stdout(),'(a)') 'rossby_radius being read from restart.'
              if(neutral_physics_debug) then 
                write (stdout(), *) &
                'checksum start rossby_radius', mpp_chksum(rossby_radius(isc:iec,jsc:jec)*Grd%tmask(isc:iec,jsc:jec,1))
              endif
          endif
      endif

      if(agm_closure .and. agm_closure_length_bczone) then 
          if(file_exist('INPUT/ocean_neutral.res.nc')) then 
              call read_data('INPUT/ocean_neutral.res','bczone_radius',bczone_radius,domain=&
                              Dom%domain2d, timelevel=1, append_pelist_name=ens_ocean)
              call mpp_update_domains(bczone_radius,Dom%domain2d)
              write (stdout(),'(a)') 'bczone_radius being read from restart.'
              if(neutral_physics_debug) then 
                  write (stdout(), *) &
                  'checksum start bczone_radius', mpp_chksum(bczone_radius(isc:iec,jsc:jec)*Grd%tmask(isc:iec,jsc:jec,1))
              endif
          endif
      endif

  endif


  ! Compute maximum stable neutral slope available for neutral diffusion.

  ! Assume the namelist aredi is max realized even if using non-constant 
  ! diffusivity scheme. 
  i_delta = isc; j_delta = jsc; k_delta = 1
  ft = 0.5/(aredi*dtime + epsln)
  delta_iso = Thickness%dht(i_delta,j_delta,k_delta,1)*Grd%dxt(i_delta,j_delta)*ft
  do k=1,nk
    do j=jsc,jec
      do i=isc,iec
        if(Grd%tmask(i,j,k) /= 0.0) then 
          deltaX = Grd%dxt(i,j)*Thickness%dht(i,j,k,1)*ft
          deltaY = Grd%dyt(i,j)*Thickness%dht(i,j,k,1)*ft
          if (delta_iso >= deltaX .or. delta_iso >= deltaY) then
            i_delta = i; j_delta = j; k_delta = k
            delta_iso = min(deltaX,deltaY)
          endif
        endif  
      enddo
    enddo
  enddo
  delta_iso  = delta_iso + 1.e-6*mpp_pe() ! to separate redundancies
  delta_iso0 = delta_iso
  call mpp_min (delta_iso)

  ! show most unstable location
  if (delta_iso == delta_iso0 .and. neutral_physics_on) then
     write(unit,'(a)')
     write(unit,'(a)')'---Neutral direction slope check I for linear stability of neutral diffusion---'
     write(unit,'(a,e14.7)')'Assuming max neutral diffusivity aredi (m^2/sec) given by ', aredi
     write(unit,'(a,e14.7)')'and neutral physics time step (secs) of ', dtime
     write(unit,'(a)')'the most stringent linear stability constraint was found at the following ocean cell:'
     write(unit,'(a,i4,a,i4,a,e14.7)')'long(',i_delta,',',j_delta,')                   = ',Grd%xt(i_delta,j_delta)
     write(unit,'(a,i4,a,i4,a,e14.7)')'lat (',i_delta,',',j_delta,')                   = ',Grd%yt(i_delta,j_delta)
     write(unit,'(a,i4,a,i4,a,i4,a,e14.7)')'thick(',i_delta,',',j_delta,',',k_delta,') = ',Thickness%dht(i_delta,j_delta,k_delta,1)
     write(unit,'(a,e14.7,a)')'delta_iso           = ',delta_iso,' is the maximum neutral direction slope' 
     write(unit,'(a)')'available for linear stability of the neutral diffusion scheme.'
     write(unit,'(a)')'The namelist parameter smax should conservatively be <= delta_iso.'
     if(smax >= delta_iso) then 
        write(unit,'(a,f10.5,a)')'==> Warning: The namelist parameter smax= ',smax, ' is >= to delta_iso.'
        write(unit,'(a)')'Linear stability of the neutral diffusion scheme may be compromised.'
     endif
     write(unit,'(a)')
  endif

  ! Compute maximum diffusivity available given a maximum slope of smax
  i_delta = isc; j_delta = jsc; k_delta = 1
  ft = 0.5/(smax*dtime + epsln)
  delta_min = Thickness%dht(i_delta,j_delta,k_delta,1)*Grd%dxt(i_delta,j_delta)
  do k=1,nk
    do j=jsc,jec
      do i=isc,iec
        if(Grd%tmask(i,j,k) /= 0.0) then        
          delta = min(Grd%dxt(i,j),Grd%dyt(i,j))*Thickness%dht(i,j,k,1)
          if (delta_min > delta) then
            i_delta = i; j_delta = j; k_delta = k
            delta_min = delta
          endif
        endif   
      enddo
    enddo
  enddo
  A_max  = ft*delta_min + 1.e-6*mpp_pe() ! to separate redundancies
  A_max0 = A_max
  call mpp_min (A_max)

  ! show most unstable location
  if (A_max == A_max0 .and. neutral_physics_on) then
     write(unit,'(a)')
     write(unit,'(a)')'---Neutral direction slope check II for linear stability of neutral diffusion---'
     write(unit,'(a,e14.7)')'Assuming maximum Redi neutral diffusion slope of ', smax
     write(unit,'(a,e14.7)')'and neutral physics time step (secs) of ', dtime
     write(unit,'(a)')'the most stringent linear stability constraint was found at the following ocean cell:'
     write(unit,'(a,i4,a,i4,a,e14.7)')'long(',i_delta,',',j_delta,')                  = ',Grd%xt(i_delta,j_delta)
     write(unit,'(a,i4,a,i4,a,e14.7)')'lat (',i_delta,',',j_delta,')                  = ',Grd%yt(i_delta,j_delta)
     write(unit,'(a,i4,a,i4,a,i4,a,e14.7)')'thick(',i_delta,',',j_delta,',',k_delta,')= ',Thickness%dht(i_delta,j_delta,k_delta,1)
     write(unit,'(a,e14.7,a)')'A_max      = ',A_max,' (m^2/sec) is the maximum neutral diffusivity'
     write(unit,'(a)')'available for linear stability of the neutral diffusion scheme.'
     write(unit,'(a)')'Conservatively, neutral diffusivities used in the model should be less than A_max.'
     write(unit,'(a)')'--------------------------------------------------------------------------------'
     write(unit,'(a)')
  endif

  ! check to be sure can run the scheme
  index_temp=-1;index_salt=-1
  do n=1, num_prog_tracers
     if (T_prog(n)%name == 'temp') index_temp = n
     if (T_prog(n)%name == 'salt') index_salt = n
  enddo
  if (index_temp == -1 .or. index_salt == -1) then 
     call mpp_error(FATAL,'==>Error: temp and/or salt not identified in call to ocean_neutral_physics_init')
  endif 

  ! initialize clock ids 
  id_clock_density_derivs        = mpp_clock_id('(Ocean neutral: density derivs)' ,grain=CLOCK_ROUTINE)
  id_clock_neutral_slopes        = mpp_clock_id('(Ocean neutral: slopes)'         ,grain=CLOCK_ROUTINE)
  id_clock_neutral_blayer        = mpp_clock_id('(Ocean neutral: blayer)'         ,grain=CLOCK_ROUTINE)
  id_clock_fz_terms              = mpp_clock_id('(Ocean neutral: fz-terms)'       ,grain=CLOCK_ROUTINE)
  id_clock_fx_flux               = mpp_clock_id('(Ocean neutral: fx-flux)'        ,grain=CLOCK_ROUTINE)
  id_clock_fy_flux               = mpp_clock_id('(Ocean neutral: fy-flux)'        ,grain=CLOCK_ROUTINE)
  id_clock_fz_flux               = mpp_clock_id('(Ocean neutral: fz-flux)'        ,grain=CLOCK_ROUTINE)
  id_clock_compute_eady_rate     = mpp_clock_id('(Ocean neutral: eady rate)'      ,grain=CLOCK_ROUTINE)
  id_clock_compute_baroclinicity = mpp_clock_id('(Ocean neutral: baroclinic)'     ,grain=CLOCK_ROUTINE)
  id_clock_compute_rossby_radius = mpp_clock_id('(Ocean neutral: rossby radius)'  ,grain=CLOCK_ROUTINE)
  id_clock_compute_bczone_radius = mpp_clock_id('(Ocean neutral: bc zone radius)' ,grain=CLOCK_ROUTINE)
  id_clock_compute_diffusivity   = mpp_clock_id('(Ocean neutral: diffusivity)'    ,grain=CLOCK_ROUTINE)


  ! register fields for diagnostic output 

  id_k33_explicit = -1
  id_k33_explicit   = register_diag_field ('ocean_model', 'k33_explicit', Grd%tracer_axes_wt(1:3), &
                      Time%model_time, 'K33_explicit tensor element', 'm^2/sec',&
                      missing_value=-10.0, range=(/-10.0,1.e20/))

  id_agm = -1
  id_agm   = register_diag_field ('ocean_model', 'agm', Grd%tracer_axes(1:2), Time%model_time, &
             'GM diffusivity at surface', 'm^2/sec', missing_value=-10.0, range=(/-10.0,1.e10/))

  id_ah_bdy = -1
  id_ah_bdy = register_static_field ('ocean_model', 'ah_bdy', Grd%tracer_axes(1:2), &
             'Horz diffusivity in surface neutral bdy', 'm^2/sec', missing_value=-10.0, range=(/-10.0,1.e10/))
  if (id_ah_bdy > 0)  used = send_data(id_ah_bdy, ah_array(:,:), &
                             Time%model_time, rmask=Grd%tmask(:,:,1), &
                             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

  id_agm_3d = -1
  id_agm_3d   = register_diag_field ('ocean_model', 'agm_3d', Grd%tracer_axes(1:3), Time%model_time, &
             '3d GM diffusivity', 'm^2/sec', missing_value=-10.0, range=(/-10.0,1.e10/))

  id_buoy = -1
  id_buoy  = register_diag_field ('ocean_model', 'buoy', Grd%tracer_axes_wt(1:3), Time%model_time, &
             'squared buoyancy frequency', '1/s^2', missing_value=-10.0, range=(/-10.0,1.e10/))

  id_tx_trans_gm = -1
  id_tx_trans_gm = register_diag_field ('ocean_model', 'tx_trans_gm', Grd%tracer_axes_flux_x(1:3), Time%model_time, &
             'T-cell volume i-transport from GM', 'Sv', missing_value=-1.e8, range=(/-1e8,1.e8/))

  id_ty_trans_gm = -1
  id_ty_trans_gm = register_diag_field ('ocean_model', 'ty_trans_gm', Grd%tracer_axes_flux_y(1:3), Time%model_time, &
             'T-cell volume j-transport from GM', 'Sv', missing_value=-1.e8, range=(/-1e8,1.e8/))

  id_depth_blayer_base = -1
  id_depth_blayer_base = register_diag_field ('ocean_model','depth_blayer_base', Grd%tracer_axes(1:2), Time%model_time, &
             'Depth of neutral physics surface bdy-layer', 'm', missing_value=-1.0, range=(/-1.0,1.e8/))

  id_slope_blayer_base = -1
  id_slope_blayer_base = register_diag_field ('ocean_model','slope_blayer_base', Grd%tracer_axes(1:2), Time%model_time, &
             'abs(slope) at neutral physics bdy-layer base', 'dimensionless', missing_value=-1.0, range=(/-1.0,1.e8/))

  id_eddy_depth = -1
  id_eddy_depth = register_diag_field ('ocean_model','eddy_depth', Grd%tracer_axes(1:2), Time%model_time, &
             'eddy penetration depth', 'm', missing_value=-1.0, range=(/-1.0,1.e8/))

  id_steep_depth = -1
  id_steep_depth = register_diag_field ('ocean_model','steep_depth', Grd%tracer_axes(1:2), Time%model_time, &
             'Depth at base of steep slope surface bdy layer', 'm', missing_value=-1.0, range=(/-1.0,1.e8/))

  id_slope31 = -1
  id_slope31 = register_diag_field ('ocean_model', 'slope31', Grd%tracer_axes(1:3), Time%model_time, &
             'neutral slope -(rho_x/rho_z)', 'dimensionless', missing_value=-1.e10, range=(/-1.e10,1.e10/))

  id_slope32 = -1
  id_slope32 = register_diag_field ('ocean_model', 'slope32', Grd%tracer_axes(1:3), Time%model_time, &
             'neutral slope -(rho_y/rho_z)', 'dimensionless', missing_value=-1.e10, range=(/-1.e10,1.e10/))

  if(gm_velocity_save) then

    id_ustar = -1
    id_ustar = register_diag_field ('ocean_model', 'ugm', Grd%vel_axes_uv(1:3), Time%model_time, &
               'GM zonal velocity', 'm/sec', missing_value=-10.0, range=(/-10.0,10.0/))

    id_vstar = -1    
    id_vstar = register_diag_field ('ocean_model', 'vgm', Grd%vel_axes_uv(1:3), Time%model_time, &
               'GM merid velocity', 'm/sec', missing_value=-10.0, range=(/-10.0,10.0/))

    id_wstar = -1
    id_wstar = register_diag_field ('ocean_model', 'wgm', Grd%vel_axes_wt(1:3), Time%model_time, &
               'GM vert velocity (T-cell bottom)', 'm/sec', missing_value=-10.0, range=(/-10.0,10.0/))
  endif

  if(agm_closure) then

    id_rossby = -1      
    id_rossby = register_diag_field ('ocean_model', 'rossby', Grd%tracer_axes(1:2), Time%model_time, &
                 'Rossby radius', 'm', missing_value=-10.0, range=(/-10.0,1.e10/))

    id_bczone = -1      
    id_bczone = register_diag_field ('ocean_model', 'bczone', Grd%tracer_axes(1:2), Time%model_time, &
                 'Radius of baroclinic zone', 'm', missing_value=-10.0, range=(/-10.0,1.e10/))

    id_eady = -1          
    id_eady = register_diag_field ('ocean_model', 'eady', Grd%tracer_axes(1:2), Time%model_time, &
                 'Eady growth rate', 's^-1', missing_value=-10.0, range=(/-10.0,1.e10/))

    id_baroclinicity = -1          
    id_baroclinicity = register_diag_field ('ocean_model', 'baroclinicity', Grd%tracer_axes(1:2), Time%model_time, &
                 'Baroclinicity (vert ave horz dens gradient)', 'kg/m^4', missing_value=-1e10, range=(/-1e10,1.e10/))

    id_growth_rate_baroclinic = -1          
    id_growth_rate_baroclinic = register_diag_field ('ocean_model', 'growth_rate_baroclinic', &
                  Grd%tracer_axes(1:2), Time%model_time, &
                 'growth rate using baroclinicity', 's^-1', missing_value=-10.0, range=(/-10.0,1.e10/))

    id_agm_growth = -1          
    id_agm_growth = register_diag_field ('ocean_model', 'agm_growth', Grd%tracer_axes(1:2), Time%model_time, &
                 'effective growth rate for agm', 's^-1', missing_value=-10.0, range=(/-10.0,1.e10/))

    id_agm_length = -1          
    id_agm_length = register_diag_field ('ocean_model', 'agm_length', Grd%tracer_axes(1:2), Time%model_time, &
                 'effective length scale for agm', 'm', missing_value=-10.0, range=(/-10.0,1.e10/))

    id_agm_qg = -1          
    id_agm_qg = register_diag_field ('ocean_model', 'agm_qg', Grd%tracer_axes(1:2), Time%model_time, &
                 'agm from QG theory', 'm^2/s', missing_value=-10.0, range=(/-10.0,1.e10/))
  endif

  id_aredi = -1            
  id_aredi = register_diag_field ('ocean_model', 'aredi', Grd%tracer_axes(1:2), Time%model_time, &
               'neutral diffusivity at surface', 'm^2/sec', missing_value=0.0, range=(/-10.0,1.e20/))
  id_aredi_3d = -1            
  id_aredi_3d = register_diag_field ('ocean_model', 'aredi', Grd%tracer_axes(1:3), Time%model_time, &
               '3d neutral diffusivity', 'm^2/sec', missing_value=0.0, range=(/-10.0,1.e20/))

  allocate (id_neutral_physics(num_prog_tracers))
  allocate (id_k33_implicit(num_prog_tracers))
  id_k33_implicit    = -1
  id_neutral_physics = -1
  do n=1,num_prog_tracers
     if (T_prog(n)%name == 'temp') then
       id_neutral_physics(n) = register_diag_field ('ocean_model', 'neutral_'//trim(T_prog(n)%name), &
                               Grd%tracer_axes(1:3), Time%model_time, &
                               'thk wghtd explicit neutral-physics heating ', 'Watts/m^2', &
                               missing_value=-1.e10, range=(/-1.e10,1.e10/))
     else 
       id_neutral_physics(n) = register_diag_field ('ocean_model', 'neutral_'//trim(T_prog(n)%name), &
                               Grd%tracer_axes(1:3), Time%model_time, &
                               'thk wghtd explicit neutral-physics for '//trim(T_prog(n)%name), 'm*kg/sec', &
                               missing_value=-1.e10, range=(/-1.e10,1.e10/))
     endif 
     id_k33_implicit(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_k33_implicit', &
                          Grd%tracer_axes_wt(1:3), &
                          Time%model_time, 'K33_implicit tensor element for '//trim(T_prog(n)%name), &
                          'm^2/sec', missing_value=-10.0, range=(/-10.0,1.e10/))
  enddo 

  allocate (id_flux_x(num_prog_tracers))
  allocate (id_flux_y(num_prog_tracers))
  allocate (id_flux_x_int_z(num_prog_tracers))
  allocate (id_flux_y_int_z(num_prog_tracers))
  id_flux_x       = -1
  id_flux_y       = -1
  id_flux_x_int_z = -1
  id_flux_y_int_z = -1

  do n=1,num_prog_tracers
     if(n == index_temp) then 
         id_flux_x(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_xflux_neutral', &
              Grd%tracer_axes_flux_x(1:3), Time%model_time, 'rho_cp*neutral_xflux*dyt*dht*temp', &
              'Watt', missing_value=-1.e20, range=(/-1.e20,1.e20/))
         id_flux_y(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_yflux_neutral', &
              Grd%tracer_axes_flux_y(1:3), Time%model_time, 'rho_cp*neutral_yflux*dxt*dht*temp', &
              'Watt', missing_value=-1.e20, range=(/-1.e20,1.e20/))
         id_flux_x_int_z(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_xflux_neutral_int_z', &
              Grd%tracer_axes_flux_x(1:2), Time%model_time, 'z-integral rho_cp*neutral_xflux*dyt*temp', &
              'Watt', missing_value=-1.e20, range=(/-1.e20,1.e20/))
         id_flux_y_int_z(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_yflux_neutral_int_z', &
              Grd%tracer_axes_flux_y(1:2), Time%model_time, 'z-integral rho_cp*neutral_yflux*dyt*temp', &
              'Watt', missing_value=-1.e20, range=(/-1.e20,1.e20/))
     else
         id_flux_x(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_xflux_neutral', &
              Grd%tracer_axes_flux_x(1:3), Time%model_time, &
              'rho0*neutral_xflux*dyt*dht*tracer for'//trim(T_prog(n)%name),&
              trim(T_prog(n)%units)//' kg/sec', missing_value=-1.e20, range=(/-1.e20,1.e20/))
         id_flux_y(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_yflux_neutral', &
              Grd%tracer_axes_flux_y(1:3), Time%model_time, &
              'rho0*neutral_yflux*dxt*dht*tracer for'//trim(T_prog(n)%name),&
              trim(T_prog(n)%units)//' kg/sec', missing_value=-1.e20, range=(/-1.e20,1.e20/))
         id_flux_x_int_z(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_xflux_neutral_int_z', &
              Grd%tracer_axes_flux_x(1:2), Time%model_time, &
              'z-integral rho0*neutral_xflux*dyt*dht*tracer for'//trim(T_prog(n)%name),&
              trim(T_prog(n)%units)//' kg/sec', missing_value=-1.e20, range=(/-1.e20,1.e20/))
         id_flux_y_int_z(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_yflux_neutral_int_z', &
              Grd%tracer_axes_flux_y(1:2), Time%model_time, &
              'z-integral rho0*neutral_yflux*dxt*dht*tracer for'//trim(T_prog(n)%name),&
              trim(T_prog(n)%units)//' kg/sec', missing_value=-1.e20, range=(/-1.e20,1.e20/))
     endif

  enddo

end subroutine ocean_neutral_physics_init
! </SUBROUTINE>  NAME="ocean_neutral_physics_init"

!#######################################################################
! <FUNCTION NAME="neutral_physics">
!
! <DESCRIPTION>
! This function computes the thickness weighted time tendency for tracer 
! from neutral physics.  Full discussion and details are provided by 
! Griffies (2004).  Here is a brief summary.  
!
!---How the neutral diffusive flux components are computed:
!
! The vertical flux component is split into diagonal (3,3) and 
! off-diagonal (3,1) and (3,2) terms. The off-diagonal (3,1) and (3,2) 
! terms are included explicitly in time. The main contribution from the 
! (3,3) term to the time tendency is included implicitly in time 
! along with the usual contribution from diapycnal processes 
! (vertical mixing schemes).  This is the K33_implicit term.
! This approach is necessary with high vertical resolution, as 
! noted by Cox (1987).  However, splitting the vertical flux into 
! an implicit and explicit piece compromises the 
! integrity of the vertical flux component (see Griffies et al. 1998).
! So to minimize the disparity engendered by this split, the portion of 
! K33 that can be stably included explicitly in time is computed along 
! with the (3,1) and (3,2) terms. 
! 
! All other terms in the mixing tensor are included explicitly in time. 
!
! The off-diagonal terms in the horizontal flux components, and all terms
! in the vertical flux component, are tapered in regions of steep neutral
! slope according to the requirements of linear stability.  MOM4 allows for 
! choice of two tapering schemes:
! (a) the tanh taper of Danabasoglu and McWilliams (1995)
! (b) the quadratic scheme of Gerdes, Koberle, and Willebrand (1991)
! Linear stability is far less stringent on the diagonal (1,1) and (2,2)
! part of the horizontal flux.  Indeed, these terms in practice need
! not be tapered in steep sloped regions. The namelist 
! neutral_taper_diagonal=.false. keeps the diagnonal terms maintained 
! for all neutral slopes. This approach assists in reducing numerical
! noise in regions where the physical system experiences a lot of
! diapycnal mixing anyhow. 
!
!---How the skew diffusive flux components are computed:
!
! The GM skew flux components are purely off-diagonal.  
! They are generally tapered when neutral slope 
! is large (neutral_physics_simple=.false).
! Doing so maintains a nontrivial GM slumping effect even when the 
! neutral slopes are vertical.  The alternative neutral_physics_simple=.true.
! is the approach used in MOM3, whereby GM effects are removed 
! in steep sloped regions.  neutral_physics_simple=.false. is 
! less efficient, but has been seen to yield superior simulations.
!
! Comments: Much research as of 2003-2004 aims to clarify the 
! schemes used to meld neutral physics with physics occuring 
! in steep sloped regions (i.e., boundary layers). Methods used
! in mom4 are not expected to be the final word on these matters.
!
! </DESCRIPTION>
!
subroutine neutral_physics (Time, Thickness, rho, pressure_at_depth, T_prog, hblt_depth)

  type(ocean_time_type), intent(in)                    :: Time
  type(ocean_thickness_type), intent(in)               :: Thickness
  real, intent(in), dimension(isd:ied,jsd:jed,nk)      :: rho
  real, intent(in), dimension(isd:ied,jsd:jed,nk)      :: pressure_at_depth
  real, intent(in), dimension(isd:ied,jsd:jed)         :: hblt_depth
  type(ocean_prog_tracer_type), intent(inout)          :: T_prog(:)
  real, dimension(isd:ied,jsd:jed)                     :: tchg
  real, dimension(isd:ied,jsd:jed)                     :: tmp_flux

  integer :: i, j, k, n, nn
  integer :: tau, taum1, dtint
  integer :: kr, kpkr, kp1, km1
  real    :: dTdz_ijk1, dTdz_ijk2

  type(time_type) :: next_time
  type(time_type) :: time_inc
  type(time_type) :: cur_time

  if ( .not. module_is_initialized ) then 
     call mpp_error(FATAL, '==>Error from ocean_neutral_physics (neutral_physics): module must be initialized')
  endif 

  if (.not. neutral_physics_on) return

  if (size(T_prog(:)) /= num_prog_tracers) then 
    call mpp_error(FATAL,'==>Error from ocean_neutral_physics (neutral_physics): inconsistent size of T_prog')
  endif 

  taum1 = Time%taum1
  tau   = Time%tau

  do k=1,nk
    delqc(:,:,k,0) = Grd%fracdz(k,0)*Thickness%dht(:,:,k,tau)
    delqc(:,:,k,1) = Grd%fracdz(k,1)*Thickness%dht(:,:,k,tau)
  enddo

  do k=0,nk
    dhwtr(:,:,k) = 1.0/Thickness%dhwt(:,:,k) 
  enddo

  call mpp_clock_begin(id_clock_density_derivs)
  call density_derivs(rho, T_prog(index_salt)%field(:,:,:,taum1), &
                      T_prog(index_temp)%field(:,:,:,taum1), pressure_at_depth, &
                      Time, drhodT, drhodS) 
  call mpp_clock_end(id_clock_density_derivs)

  ! compute derivative operators 
  do n=1,num_prog_tracers
     wrk1(:,:,1) = T_prog(n)%field(:,:,1,taum1)
     do k=1,nk
        kp1 = min(k+1,nk)
        wrk1(:,:,2)          = T_prog(n)%field(:,:,k,taum1)
        dTdx(n)%field(:,:,k) = FDX_PT(wrk1(:,:,1:2),k)*FMX(Grd%tmask(:,:,k))
        dTdy(n)%field(:,:,k) = FDY_PT(wrk1(:,:,1:2),k)*FMY(Grd%tmask(:,:,k))
        dTdz(n)%field(:,:,k) = Grd%tmask(:,:,kp1)*(T_prog(n)%field(:,:,k,taum1)-T_prog(n)%field(:,:,kp1,taum1)) &
                               *dhwtr(:,:,k)
        wrk1(:,:,1)          = wrk1(:,:,2)
     enddo
  enddo

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           dTdz_ijk1 = dTdz(index_temp)%field(i,j,k) ; dTdz_ijk2 = dTdz(index_salt)%field(i,j,k)

           kr=0
           kpkr = min(k+kr,nk)
           drhodzb(i,j,k,kr)  = drhodT(i,j,kpkr)*dTdz_ijk1 + drhodS(i,j,kpkr)*dTdz_ijk2 - epsln

           kr=1
           kpkr = min(k+kr,nk)
           drhodzb(i,j,k,kr)  = drhodT(i,j,kpkr)*dTdz_ijk1 + drhodS(i,j,kpkr)*dTdz_ijk2 - epsln
           
        enddo
     enddo
  enddo

  call neutral_slopes(Time, T_prog)
  if(neutral_sine_taper .or. neutral_linear_gm_taper .or. neutral_blayer_diagnose) then 
    call neutral_blayer(Time, T_prog, hblt_depth)
  endif 
  call fz_terms(Time, Thickness, T_prog)

  if(gm_velocity_save) then
     dtint = nint(dtime)
     next_time =  increment_time(Time%model_time, dtint, 0)
     if(need_data(id_ustar, next_time) .or. &
        need_data(id_vstar, next_time) .or. &
        need_data(id_wstar, next_time)) then 
        call gm_velocity(Thickness, tau)
     endif
  endif

  ! 1.e-6 converts m^3/sec to Sverdrups 
  if(id_tx_trans_gm > 0) then 
      call mpp_update_domains(tx_trans_gm, Dom%domain2d,flags=EUPDATE)
      do k=1,nk-1
         do j=jsc,jec
            do i=isc,iec
               tx_trans_gm(i,j,k) = 1.e-6*Grd%dxuer(i,j)* &
                                   (tx_trans_gm(i+1,j,k)*Grd%due(i,j) + &
                                    tx_trans_gm(i,j,k)  *Grd%duw(i+1,j))
            enddo
         enddo
      enddo
      used = send_data (id_tx_trans_gm, tx_trans_gm(:,:,:), &
             Time%model_time, rmask=Grd%tmask(:,:,:), &
             is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif

  ! 1.e-6 converts m^3/sec to Sverdrups 
  if(id_ty_trans_gm > 0) then 
      call mpp_update_domains(ty_trans_gm, Dom%domain2d,flags=NUPDATE)
      do k=1,nk-1
         do j=jsc,jec
            do i=isc,iec
               ty_trans_gm(i,j,k) = 1.e-6*Grd%dyunr(i,j)* &
                                          (ty_trans_gm(i,j+1,k)*Grd%dun(i,j) +  &
                                           ty_trans_gm(i,j,k)  *Grd%dus(i,j+1)) 
                                        
            enddo
         enddo
      enddo
      used = send_data (id_ty_trans_gm, ty_trans_gm(:,:,:), &
             Time%model_time, rmask=Grd%tmask(:,:,:), &
             is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif
  
  do n=1,num_prog_tracers  
    fz1(n)%field(:,:)      = 0.0
    fz2(n)%field(:,:)      = 0.0
    flux_x(n)%field(:,:,:) = 0.0
    flux_y(n)%field(:,:,:) = 0.0
    T_prog(n)%wrk1(:,:,:)  = 0.0
  enddo 


  if(neutral_physics_simple) then 

      ! horizontal flux components are downgradient diffusive 
      do k=1,nk
         do nn=1,num_prog_tracers
            flux_x(nn)%field(:,:,k) = agm_array(:,:,k)*dTdx(nn)%field(:,:,k)*FMX(Thickness%dht(:,:,k,tau)*Grd%tmask(:,:,k))
            flux_y(nn)%field(:,:,k) = agm_array(:,:,k)*dTdy(nn)%field(:,:,k)*FMY(Thickness%dht(:,:,k,tau)*Grd%tmask(:,:,k))
         enddo
      enddo
      if (Grd%tripolar) then 
          do nn=1,num_prog_tracers
             call mpp_update_domains(flux_x(nn)%field(:,:,:), flux_y(nn)%field(:,:,:), Dom_flux%domain2d, &
                  gridtype=CGRID_NE, complete=T_prog(nn)%complete) 
          enddo
      endif
      do k=1,nk
         call mpp_clock_begin(id_clock_fz_flux)
         call fz_flux(T_prog,k)
         call mpp_clock_end(id_clock_fz_flux)
         do nn=1,num_prog_tracers
            tchg(:,:)  = (BDX_ET(flux_x(nn)%field(:,:,k)) + BDY_NT(flux_y(nn)%field(:,:,k))) 
            do j=jsc,jec
               do i=isc,iec
                  T_prog(nn)%wrk1(i,j,k)        = Grd%tmask(i,j,k)*(tchg(i,j) + (fz1(nn)%field(i,j)-fz2(nn)%field(i,j)))
                  T_prog(nn)%th_tendency(i,j,k) = T_prog(nn)%th_tendency(i,j,k) + T_prog(nn)%wrk1(i,j,k)
                  fz1(nn)%field(i,j)       = fz2(nn)%field(i,j)
                  flux_x(nn)%field(i,j,k)  = Grd%dyte(i,j)*flux_x(nn)%field(i,j,k)
                  flux_y(nn)%field(i,j,k)  = Grd%dxtn(i,j)*flux_y(nn)%field(i,j,k)
               enddo
            enddo
         enddo
      enddo

  else 

      do k=1,nk
         call fx_flux(Time, Thickness, T_prog, k)
         call fy_flux(Time, Thickness, T_prog, k)
      enddo
      if (Grd%tripolar) then 
          do nn=1,num_prog_tracers
             call mpp_update_domains(flux_x(nn)%field(:,:,:), flux_y(nn)%field(:,:,:), Dom_flux%domain2d, &
                  gridtype=CGRID_NE, complete=T_prog(nn)%complete) 
          enddo
      endif

      do k=1,nk
        call mpp_clock_begin(id_clock_fz_flux)
        call fz_flux(T_prog,k)
        call mpp_clock_end(id_clock_fz_flux)
        do nn=1,num_prog_tracers
           do j=jsc,jec
              do i=isc,iec
                 T_prog(nn)%wrk1(i,j,k) = & 
                   Grd%tmask(i,j,k)  &
                  *(fz1(nn)%field(i,j)-fz2(nn)%field(i,j) &
                    +( flux_x(nn)%field(i,j,k)-flux_x(nn)%field(i-1,j,k) &
                      +flux_y(nn)%field(i,j,k)-flux_y(nn)%field(i,j-1,k) )*Grd%datr(i,j) &
                   )
                 T_prog(nn)%th_tendency(i,j,k) = T_prog(nn)%th_tendency(i,j,k) + T_prog(nn)%wrk1(i,j,k)
                 fz1(nn)%field(i,j) = fz2(nn)%field(i,j)
              enddo
           enddo
        enddo
     enddo

  endif  ! endif for neutral_physics_simple


  ! compute Eady growth rate and baroclinicity for use in next time step 
  if(agm_closure) then
     call compute_eady_rate(Time%model_time)
     call compute_baroclinicity(Time%model_time)
  endif 

  ! compute rossby radius for use in next time step 
  if(agm_closure .or. neutral_sine_taper) then
     call compute_rossby_radius(Thickness, Time%model_time)
  endif 

  ! compute width of baroclinic zone for use in next time step 
  if(agm_closure .and. agm_closure_length_bczone) then
     call compute_bczone_radius(Time%model_time)
  endif 
  
  ! update closure-based diffusivity for next time step 
  if(agm_closure) then
      call compute_diffusivity(Time%model_time)
   endif

  do nn=1,num_prog_tracers

     if(id_neutral_physics(nn) > 0) then 
        used = send_data (id_neutral_physics(nn), T_prog(nn)%wrk1(:,:,:)*T_prog(nn)%conversion, &
               Time%model_time, rmask=Grd%tmask(:,:,:), &
               is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
     endif

     ! send fluxes to diag_manager 
     ! minus sign is due to a MOM-convention for physics fluxes  
     if(id_flux_x(nn) > 0) then 
        used = send_data (id_flux_x(nn), -1.0*T_prog(nn)%conversion*flux_x(nn)%field(:,:,:), &
               Time%model_time, rmask=Grd%tmask(:,:,:), &
               is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
     endif
     if(id_flux_y(nn) > 0) then 
        used = send_data (id_flux_y(nn), -1.0*T_prog(nn)%conversion*flux_y(nn)%field(:,:,:), &
               Time%model_time, rmask=Grd%tmask(:,:,:), &
               is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
     endif
     if(id_flux_x_int_z(nn) > 0) then 
        tmp_flux = 0.0
        do k=1,nk
            tmp_flux(isc:iec,jsc:jec) = tmp_flux(isc:iec,jsc:jec) + flux_x(nn)%field(isc:iec,jsc:jec,k)
        enddo        
        used = send_data (id_flux_x_int_z(nn), -1.0*T_prog(nn)%conversion*tmp_flux(:,:), &
               Time%model_time, rmask=Grd%tmask(:,:,1), &
               is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     endif
     if(id_flux_y_int_z(nn) > 0) then 
        tmp_flux = 0.0
        do k=1,nk
            tmp_flux(isc:iec,jsc:jec) = tmp_flux(isc:iec,jsc:jec) + flux_y(nn)%field(isc:iec,jsc:jec,k)
        enddo        
        used = send_data (id_flux_y_int_z(nn), -1.0*T_prog(nn)%conversion*tmp_flux(:,:), &
               Time%model_time, rmask=Grd%tmask(:,:,1), &
               is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     endif
     
  enddo

  ! more diagnostic output

  if (id_k33_explicit > 0) used = send_data (id_k33_explicit, K33_explicit(:,:,:), &
                                  Time%model_time, rmask=Grd%tmask(:,:,:), &
                                  is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

  if (id_agm > 0) used = send_data (id_agm, agm_array(:,:,1), &
                         Time%model_time, rmask=Grd%tmask(:,:,1), &
                         is_in=isc, js_in=jsc, ie_in=iec, je_in=jec) 

  if (id_agm_3d > 0) used = send_data (id_agm_3d, agm_array(:,:,:), &
                            Time%model_time, rmask=Grd%tmask(:,:,:), &
                            is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk) 

  if (id_aredi > 0) used = send_data (id_aredi, aredi_array(:,:,1), &
                           Time%model_time, rmask=Grd%tmask(:,:,1), &
                           is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

  if (id_aredi_3d > 0) used = send_data (id_aredi_3d, aredi_array(:,:,:), &
                              Time%model_time, rmask=Grd%tmask(:,:,:), &
                              is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

  if(id_buoy > 0) used = send_data (id_buoy, -0.5*gravrho0r*(drhodzb(:,:,:,0)+drhodzb(:,:,:,1)), &
                         Time%model_time, rmask=Grd%tmask(:,:,:), &
                         is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk) 

  if(gm_velocity_save) then 
     if(need_data(id_ustar, next_time)) then
        if (id_ustar > 0) used = send_data (id_ustar, ustar(:,:,:), &
                                 Time%model_time, rmask=Grd%umask(:,:,:), &
                                 is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
     endif
     if(need_data(id_vstar, next_time)) then
        if (id_vstar > 0) used = send_data (id_vstar, vstar(:,:,:), &
                                 Time%model_time, rmask=Grd%umask(:,:,:), &
                                 is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
     endif
     if(need_data(id_wstar, next_time)) then
        if (id_wstar > 0) used = send_data (id_wstar, wstar(:,:,:), &
                                 Time%model_time, rmask=Grd%tmask(:,:,:), &
                                 is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
     endif
  endif


end subroutine neutral_physics
! </FUNCTION> NAME="neutral_physics"


!#######################################################################
! <SUBROUTINE NAME="neutral_slopes">
!
! <DESCRIPTION>
! Subroutine computes the neutral slopes for the triads associated 
! with the vertical flux component.  No tapering applied here. 
! Array tensor_31 initially holds the x-slope used for flux component fz.
! Array tensor_32 initially holds the y-slope used for flux component fz.
! </DESCRIPTION>
!
subroutine neutral_slopes(Time, T_prog)

  type(ocean_time_type), intent(in)        :: Time
  type(ocean_prog_tracer_type), intent(in) :: T_prog(num_prog_tracers)

  integer :: i, j, k, kp1, ip, jq, kr, kpkr
  real :: drhodT_ijk, drhodS_ijk, drhodT_ijkp1, drhodS_ijkp1
  real :: tmask_ijkp1

  real :: dTdx_ijk1, dTdx_ijk2, dTdx_im1jk1, dTdx_im1jk2, dTdx_ijkp11, dTdx_ijkp12, dTdx_im1jkp11, dTdx_im1jkp12
  real :: dTdy_ijk1, dTdy_ijk2, dTdy_ijm1k1, dTdy_ijm1k2, dTdy_ijkp11, dTdy_ijkp12, dTdy_ijm1kp11, dTdy_ijm1kp12
  real :: drhodzb_ijk0, drhodzb_ijk1
  real :: drhodzbr_ijk0, drhodzbr_ijk1

  call mpp_clock_begin(id_clock_neutral_slopes)

  do k=1,nk
     kp1 = min(k+1,nk)
     do j=jsc,jec
        do i=isc,iec

           tmask_ijkp1    = Grd%tmask(i,j,kp1)

           drhodT_ijk     = drhodT(i,j,k)
           drhodS_ijk     = drhodS(i,j,k)
           drhodT_ijkp1   = drhodT(i,j,kp1)
           drhodS_ijkp1   = drhodS(i,j,kp1)

           dTdx_ijk1     = dTdx(index_temp)%field(i,j,k)
           dTdx_ijk2     = dTdx(index_salt)%field(i,j,k)
           dTdx_im1jk1   = dTdx(index_temp)%field(i-1,j,k)
           dTdx_im1jk2   = dTdx(index_salt)%field(i-1,j,k)
           dTdx_ijkp11   = dTdx(index_temp)%field(i,j,kp1)
           dTdx_ijkp12   = dTdx(index_salt)%field(i,j,kp1)
           dTdx_im1jkp11 = dTdx(index_temp)%field(i-1,j,kp1)
           dTdx_im1jkp12 = dTdx(index_salt)%field(i-1,j,kp1)

           dTdy_ijk1     = dTdy(index_temp)%field(i,j,k)
           dTdy_ijk2     = dTdy(index_salt)%field(i,j,k)
           dTdy_ijm1k1   = dTdy(index_temp)%field(i,j-1,k)
           dTdy_ijm1k2   = dTdy(index_salt)%field(i,j-1,k)
           dTdy_ijkp11   = dTdy(index_temp)%field(i,j,kp1)
           dTdy_ijkp12   = dTdy(index_salt)%field(i,j,kp1)
           dTdy_ijm1kp11 = dTdy(index_temp)%field(i,j-1,kp1)
           dTdy_ijm1kp12 = dTdy(index_salt)%field(i,j-1,kp1)

           drhodzbr_ijk0 = 1.0/drhodzb(i,j,k,0)
           drhodzbr_ijk1 = 1.0/drhodzb(i,j,k,1)


           ! ip=jq=0

           !   kr=0
           tensor_31(i,j,k,0,0) = -tmask_ijkp1 &
                *(drhodT_ijk*dTdx_im1jk1 + drhodS_ijk*dTdx_im1jk2) &
                *drhodzbr_ijk0
           tensor_32(i,j,k,0,0) = -tmask_ijkp1 &
                *(drhodT_ijk*dTdy_ijm1k1 + drhodS_ijk*dTdy_ijm1k2) &
                *drhodzbr_ijk0
           !   kr=1 
           tensor_31(i,j,k,0,1) = -tmask_ijkp1 &
                *(drhodT_ijkp1*dTdx_im1jkp11 + drhodS_ijkp1*dTdx_im1jkp12) &
                *drhodzbr_ijk1

           tensor_32(i,j,k,0,1) = -tmask_ijkp1 &
                *(drhodT_ijkp1*dTdy_ijm1kp11 + drhodS_ijkp1*dTdy_ijm1kp12) &
                *drhodzbr_ijk1


           ! ip=jq=1

           !   kr=0 
           tensor_31(i,j,k,1,0) = -tmask_ijkp1 &
                *(drhodT_ijk*dTdx_ijk1 + drhodS_ijk*dTdx_ijk2) &
                *drhodzbr_ijk0

           tensor_32(i,j,k,1,0) = -tmask_ijkp1 &
                *(drhodT_ijk*dTdy_ijk1 + drhodS_ijk*dTdy_ijk2) &
                *drhodzbr_ijk0
           !   kr=1 
           tensor_31(i,j,k,1,1) = -tmask_ijkp1 &
                *(drhodT_ijkp1*dTdx_ijkp11 + drhodS_ijkp1*dTdx_ijkp12) &
                *drhodzbr_ijk1

           tensor_32(i,j,k,1,1) = -tmask_ijkp1 &
                *(drhodT_ijkp1*dTdy_ijkp11 + drhodS_ijkp1*dTdy_ijkp12) &
                *drhodzbr_ijk1

        enddo
     enddo
  enddo

  ! send to diagnostic manager 
  if (id_slope31 > 0) then 
       wrk1(isc:iec,jsc:jec,:) = 0.25* &
         (tensor_31(isc:iec,jsc:jec,:,0,0) + tensor_31(isc:iec,jsc:jec,:,0,1) + &
          tensor_31(isc:iec,jsc:jec,:,1,0) + tensor_31(isc:iec,jsc:jec,:,1,1))
       used = send_data (id_slope31, wrk1(:,:,:), &
              Time%model_time, rmask=Grd%tmask(:,:,:), &
              is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif 
  if (id_slope32 > 0) then 
       wrk1(isc:iec,jsc:jec,:) = 0.25* &
         (tensor_32(isc:iec,jsc:jec,:,0,0) + tensor_32(isc:iec,jsc:jec,:,0,1) + &
          tensor_32(isc:iec,jsc:jec,:,1,0) + tensor_32(isc:iec,jsc:jec,:,1,1))
       used = send_data (id_slope32, wrk1(:,:,:), &
              Time%model_time, rmask=Grd%tmask(:,:,:), &
              is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif  

  call mpp_clock_end(id_clock_neutral_slopes)

end subroutine neutral_slopes
! </SUBROUTINE> NAME="neutral_slopes"


!#######################################################################
! <SUBROUTINE NAME="neutral_blayer">
!
! <DESCRIPTION>
! Subroutine computes the boundary layer as determined by 
! 1. steep neutral slopes
! 2. depth within which typical mesoscale eddies are partially outcropped
! 3. depth within which vertical mixing scheme (e.g., kpp) computes a boundary layer
!
! Note: Only consider surface boundary layers here.  
! 
! </DESCRIPTION>
!
subroutine neutral_blayer(Time, T_prog, hblt_depth)

  type(ocean_time_type), intent(in)            :: Time
  type(ocean_prog_tracer_type), intent(in)     :: T_prog(num_prog_tracers)
  real, intent(in), dimension(isd:ied,jsd:jed) :: hblt_depth

  logical :: slopeok(isc:iec)
  integer :: i, j, k, ip, jq, kr, kk, kkpkr
  real :: depth, slope, absslope 
  real :: steep_depth(isd:ied,jsd:jed)
  real :: slope_blayer_baseA(isd:ied,jsd:jed)
  real :: slope_blayer_baseB(isd:ied,jsd:jed)
  real :: thick_31(isd:ied,jsd:jed,0:1,0:1)
  real :: thick_32(isd:ied,jsd:jed,0:1,0:1)
  real :: thick_13(isd:ied,jsd:jed,0:1,0:1)
  real :: thick_23(isd:ied,jsd:jed,0:1,0:1)
  real :: slope_31(isd:ied,jsd:jed,0:1,0:1)
  real :: slope_32(isd:ied,jsd:jed,0:1,0:1)
  real :: slope_13(isd:ied,jsd:jed,0:1,0:1)
  real :: slope_23(isd:ied,jsd:jed,0:1,0:1)

  call mpp_clock_begin(id_clock_neutral_blayer)

  thick_31 = Grd%zw(1)
  thick_32 = Grd%zw(1)
  thick_13 = Grd%zw(1)
  thick_23 = Grd%zw(1)
  slope_31 = 0.0
  slope_32 = 0.0
  slope_13 = 0.0
  slope_23 = 0.0

  slope_blayer_baseA = 0.0
  slope_blayer_baseB = 0.0
  steep_depth        = 0.0

  eddy_depth         = 0.0
  slope_blayer_base  = 0.0 
  depth_blayer_base  = 0.0


  if(neutral_sine_taper .or. neutral_blayer_diagnose) then 

  ! Determine depth over which mesoscale eddies feel the ocean 
  ! surface.  This depth is a function of the neutral slope 
  ! and the Rossby radius.  This depth is called "eddy_depth".
  ! The algorithm for computing this depth is taken from 
  ! the appendix to Large etal, 1997 JPO vol 27, 2418-2447. 
  ! 
  ! In addition to considering mesoscale eddy lengths,
  ! include the possibility that the diabatic vertical
  ! mixing (e.g., kpp) produces a mixed layer depth that is 
  ! deeper than the depth that mesoscale eddies feel the ocean 
  ! surface.  Include this hblt_depth in the considerations so 
  ! to determine the depth of this generalized "boundary layer" 
  ! and the neutral slope at the base of the boundary layer. 

      ! 31-triads 
      do ip=0,1
         do kr=0,1
            do j=jsc,jec
               do i=isc,iec 

                  do kk=1,nk-1
                     absslope = abs(tensor_31(i,j,kk,ip,kr))
                     depth    = max(rossby_radius(i,j)*absslope, hblt_depth(i,j))
                     if(depth > Grd%zw(kk)) then
                         thick_31(i,j,ip,kr) = Grd%zw(kk)
                         slope_31(i,j,ip,kr) = absslope
                     else 
                         exit
                     endif
                  enddo

               enddo
            enddo
         enddo
      enddo

      ! 13-triads 
      do kr=0,1
         do ip=0,1 
            do j=jsc,jec
               slopeok(:)=.false.

kkloop13a:     do kk=1,nk
                  kkpkr = min(kk+kr,nk)
                  do i=isc,iec

                     slope = -Grd%tmask(i+ip,j,kkpkr)&
                          *(drhodT(i+ip,j,kk)*dTdx(index_temp)%field(i,j,kk)   +&
                          drhodS(i+ip,j,kk)*dTdx(index_salt)%field(i,j,kk)) &
                          /(drhodT(i+ip,j,kk)*dTdz(index_temp)%field(i+ip,j,kk-1+kr) +&
                          drhodS(i+ip,j,kk)*dTdz(index_salt)%field(i+ip,j,kk-1+kr) - epsln)                           
                     absslope = abs(slope)  

                     depth = max(rossby_radius(i,j)*absslope, hblt_depth(i,j))
                     if(depth > Grd%zw(kk) .and. .not.slopeok(i)) then 
                         thick_13(i,j,ip,kr) = Grd%zw(kk)
                         slope_13(i,j,ip,kr) = absslope
                     else 
                         slopeok(i)=.true.
                     endif

                  enddo
                  if(all(slopeok(:))) exit kkloop13a
               enddo kkloop13a

            enddo
         enddo
      enddo

      ! 32-triads 
      do jq=0,1
         do kr=0,1
            do j=jsc,jec
               do i=isc,iec 

                  do kk=1,nk-1
                     absslope = abs(tensor_32(i,j,kk,jq,kr))
                     depth    = max(rossby_radius(i,j)*absslope, hblt_depth(i,j))
                     if(depth > Grd%zw(kk)) then
                         thick_32(i,j,jq,kr) = Grd%zw(kk)
                         slope_32(i,j,jq,kr) = absslope
                     else 
                         exit
                     endif
                  enddo

               enddo
            enddo
         enddo
      enddo


      ! 23-triads 
      do kr=0,1
         do jq=0,1  
            do j=jsc,jec
               slopeok(:)=.false.

kkloop23a:     do kk=1,nk
                  do i=isc,iec
                     kkpkr = min(kk+kr,nk)

                     slope = -Grd%tmask(i,j+jq,kkpkr)&
                          *(drhodT(i,j+jq,kk)*dTdy(index_temp)%field(i,j,kk)+ &
                          drhodS(i,j+jq,kk)*dTdy(index_salt)%field(i,j,kk)) & 
                          /(drhodT(i,j+jq,kk)*dTdz(index_temp)%field(i,j+jq,kk-1+kr) +&
                          drhodS(i,j+jq,kk)*dTdz(index_salt)%field(i,j+jq,kk-1+kr) - epsln)  
                     absslope = abs(slope)  

                     depth = max(rossby_radius(i,j)*absslope, hblt_depth(i,j))
                     if(depth > Grd%zw(kk) .and. .not.slopeok(i)) then 
                         thick_23(i,j,jq,kr) = Grd%zw(kk)
                         slope_23(i,j,jq,kr) = absslope
                     else
                         slopeok(i)=.true.
                     endif

                  enddo
                  if(all(slopeok(:))) exit kkloop23a
               enddo kkloop23a
            enddo
         enddo
      enddo

      ! max of triad depth defines eddy_depth. average slope
      ! defines slope at boundary layer base.  This method 
      ! maximizes eddy_depth, and regulates slope at base.
      do j=jsc,jec
         do i=isc,iec

            eddy_depth(i,j) = max(Grd%zw(1), &
                                  thick_31(i,j,0,0), thick_31(i,j,0,1), &
                                  thick_31(i,j,1,0), thick_31(i,j,1,1), &
                                  thick_32(i,j,0,0), thick_32(i,j,0,1), &
                                  thick_32(i,j,1,0), thick_32(i,j,1,1), &
                                  thick_13(i,j,0,0), thick_13(i,j,0,1), &
                                  thick_13(i,j,1,0), thick_13(i,j,1,1), &
                                  thick_23(i,j,0,0), thick_23(i,j,0,1), &
                                  thick_23(i,j,1,0), thick_23(i,j,1,1))

            ! make sure that none of the slopes are larger than smax 
            do kr=0,1
              do ip=0,1      
                 slope_31(i,j,ip,kr) = min(smax, slope_31(i,j,ip,kr)) 
                 slope_13(i,j,ip,kr) = min(smax, slope_13(i,j,ip,kr)) 
                 slope_32(i,j,ip,kr) = min(smax, slope_32(i,j,ip,kr)) 
                 slope_23(i,j,ip,kr) = min(smax, slope_23(i,j,ip,kr)) 
               enddo
            enddo

            ! average of the 16 slopes defines slope at boundary layer base  
            slope_blayer_baseA(i,j) = 0.0
            do kr=0,1
              do ip=0,1      
                slope_blayer_baseA(i,j) = slope_blayer_baseA(i,j) + &
                                          slope_31(i,j,ip,kr) + slope_13(i,j,ip,kr) + &
                                          slope_32(i,j,ip,kr) + slope_23(i,j,ip,kr)
               enddo
            enddo
            slope_blayer_baseA(i,j) = slope_blayer_baseA(i,j)/16.0

         enddo
      enddo

  endif  ! endif for neutral_sine_taper .or. neutral_blayer_diagnose


  
  if(neutral_linear_gm_taper .or. neutral_blayer_diagnose) then 

  ! Determine depth of surface boundary layer as defined 
  ! by depth where neutral slopes are larger than smax.
      
      thick_31 = Grd%zw(1)
      thick_32 = Grd%zw(1)
      thick_13 = Grd%zw(1)
      thick_23 = Grd%zw(1)

      ! 31-triads 
      do ip=0,1
         do kr=0,1
            do j=jsc,jec
               do i=isc,iec 

                  do kk=1,nk-1
                     absslope = abs(tensor_31(i,j,kk,ip,kr))
                     if(absslope > smax) then
                         thick_31(i,j,ip,kr) = Grd%zw(kk)
                         slope_31(i,j,ip,kr) = absslope
                     else 
                         exit
                     endif
                  enddo

               enddo
            enddo
         enddo
      enddo

      ! 13-triads 
      do kr=0,1
         do ip=0,1 
            do j=jsc,jec
               slopeok(:)=.false.

kkloop13b:     do kk=1,nk
                  kkpkr = min(kk+kr,nk)
                  do i=isc,iec

                     slope = -Grd%tmask(i+ip,j,kkpkr)&
                          *(drhodT(i+ip,j,kk)*dTdx(index_temp)%field(i,j,kk)   +&
                          drhodS(i+ip,j,kk)*dTdx(index_salt)%field(i,j,kk)) &
                          /(drhodT(i+ip,j,kk)*dTdz(index_temp)%field(i+ip,j,kk-1+kr) +&
                          drhodS(i+ip,j,kk)*dTdz(index_salt)%field(i+ip,j,kk-1+kr) - epsln)                           
                     absslope = abs(slope) 
 
                     if(absslope > smax .and. .not.slopeok(i)) then 
                         thick_13(i,j,ip,kr) = Grd%zw(kk)
                         slope_13(i,j,ip,kr) = absslope
                     else 
                         slopeok(i)=.true.
                     endif

                  enddo
                  if(all(slopeok(:))) exit kkloop13b
               enddo kkloop13b

            enddo

         enddo
      enddo

      ! 32-triads 
      do jq=0,1
         do kr=0,1
            do j=jsc,jec
               do i=isc,iec 

                  do kk=1,nk-1
                     absslope = abs(tensor_32(i,j,kk,jq,kr))
                     if(absslope > smax) then
                         thick_32(i,j,jq,kr) = Grd%zw(kk)
                         slope_32(i,j,jq,kr) = absslope
                     else 
                         exit
                     endif
                  enddo

               enddo
            enddo
         enddo
      enddo

      ! 23-triads 
      do kr=0,1
         do jq=0,1  
            do j=jsc,jec
               slopeok(:)=.false.

kkloop23b:     do kk=1,nk
                  do i=isc,iec
                     kkpkr = min(kk+kr,nk)

                     slope = -Grd%tmask(i,j+jq,kkpkr)&
                          *(drhodT(i,j+jq,kk)*dTdy(index_temp)%field(i,j,kk)+ &
                          drhodS(i,j+jq,kk)*dTdy(index_salt)%field(i,j,kk)) & 
                          /(drhodT(i,j+jq,kk)*dTdz(index_temp)%field(i,j+jq,kk-1+kr) +&
                          drhodS(i,j+jq,kk)*dTdz(index_salt)%field(i,j+jq,kk-1+kr) - epsln)  
                     absslope = abs(slope)  

                     if(absslope > smax .and. .not.slopeok(i)) then 
                         thick_23(i,j,jq,kr) = Grd%zw(kk)
                         slope_23(i,j,jq,kr) = absslope
                     else
                         slopeok(i)=.true.
                     endif

                  enddo
                  if(all(slopeok(:))) exit kkloop23b
               enddo kkloop23b

            enddo
         enddo
      enddo

      ! max of triad depth defines steep_depth. average slope
      ! defines slope at boundary layer base.  This method 
      ! maximizes steep_depth, and regulates slope at base.
      do j=jsc,jec
         do i=isc,iec

            steep_depth(i,j) =  max(Grd%zw(1), &
                                    thick_31(i,j,0,0), thick_31(i,j,0,1), &
                                    thick_31(i,j,1,0), thick_31(i,j,1,1), &
                                    thick_32(i,j,0,0), thick_32(i,j,0,1), &
                                    thick_32(i,j,1,0), thick_32(i,j,1,1), &
                                    thick_13(i,j,0,0), thick_13(i,j,0,1), &
                                    thick_13(i,j,1,0), thick_13(i,j,1,1), &
                                    thick_23(i,j,0,0), thick_23(i,j,0,1), &
                                    thick_23(i,j,1,0), thick_23(i,j,1,1))

            ! make sure that none of the slopes are larger than smax 
            do kr=0,1
              do ip=0,1      
                 slope_31(i,j,ip,kr) = min(smax, slope_31(i,j,ip,kr)) 
                 slope_13(i,j,ip,kr) = min(smax, slope_13(i,j,ip,kr)) 
                 slope_32(i,j,ip,kr) = min(smax, slope_32(i,j,ip,kr)) 
                 slope_23(i,j,ip,kr) = min(smax, slope_23(i,j,ip,kr)) 
               enddo
            enddo

            ! average of the 16 slopes defines slope at boundary layer base  
            slope_blayer_baseB(i,j) = 0.0
            do kr=0,1
              do ip=0,1      
                slope_blayer_baseB(i,j) = slope_blayer_baseB(i,j) + &
                                          slope_31(i,j,ip,kr) + slope_13(i,j,ip,kr) + &
                                          slope_32(i,j,ip,kr) + slope_23(i,j,ip,kr)
               enddo
            enddo
            slope_blayer_baseB(i,j) = slope_blayer_baseB(i,j)/16.0

         enddo
      enddo

  ! The maximum of steep_depth and eddy_depth define 
  ! neutral_blayer_depth.  For depths shallower than 
  ! neutral_blayer_depth, horizontal gm velocity is depth  
  ! independent and vertical gm velocity linearly decreases
  ! to zero as the surface is reached.   

      do j=jsc,jec
         do i=isc,iec

            if(steep_depth(i,j) > eddy_depth(i,j)) then 
                slope_blayer_base(i,j) = slope_blayer_baseB(i,j)
                depth_blayer_base(i,j) = steep_depth(i,j)
            else 
                slope_blayer_base(i,j) = slope_blayer_baseA(i,j)
                depth_blayer_base(i,j) = eddy_depth(i,j)
            endif

         enddo
      enddo

  endif  ! endif for neutral_linear_gm_taper .or. neutral_blayer_diagnose


  if (id_depth_blayer_base >  0) then 
       used = send_data (id_depth_blayer_base, depth_blayer_base(:,:), &
              Time%model_time, rmask=Grd%tmask(:,:,1), &
              is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif 

  if (id_slope_blayer_base >  0) then 
       used = send_data (id_slope_blayer_base, slope_blayer_base(:,:), &
              Time%model_time, rmask=Grd%tmask(:,:,1), &
              is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif 

  if (id_eddy_depth >  0) then 
       used = send_data (id_eddy_depth, eddy_depth(:,:), &
              Time%model_time, rmask=Grd%tmask(:,:,1), &
              is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif 

  if (id_steep_depth >  0) then 
       used = send_data (id_steep_depth, steep_depth(:,:), &
              Time%model_time, rmask=Grd%tmask(:,:,1), &
              is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif 


  ! update domains or zero the depths, depending on the chosen method 
  if(neutral_linear_gm_taper) then 
      call mpp_update_domains(depth_blayer_base, Dom%domain2d) 
      call mpp_update_domains(slope_blayer_base, Dom%domain2d) 
  else 
      depth_blayer_base = 0.0
  endif
  if(neutral_sine_taper) then 
      call mpp_update_domains(eddy_depth, Dom%domain2d) 
  else 
      eddy_depth = 0.0
  endif


  call mpp_clock_end(id_clock_neutral_blayer)

end subroutine neutral_blayer
! </SUBROUTINE> NAME="neutral_blayer"




!#######################################################################
! <SUBROUTINE NAME="fz_terms">
!
! <DESCRIPTION>
! Subroutine computes the tracer independent pieces of the vertical 
! flux component. As a result of this routine, 
! Array tensor_31 = x-diffusivity*slope (m^2/sec) for fz
! Array tensor_32 = y-diffusivity*slope (m^2/sec) for fz 
!
! K33 is the (3,3) term in small angle Redi diffusion tensor.
! It is broken into an explicit in time piece and implicit 
! in time piece.  
!
! K33 has units m^2/sec.
!
! Also will compute the squared Eady growth rate, with the maximum slope
! contributing to this growth rate set by smax.
! </DESCRIPTION>
!
subroutine fz_terms(Time, Thickness, T_prog)

  type(ocean_time_type), intent(in)           :: Time
  type(ocean_thickness_type), intent(in)      :: Thickness
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(num_prog_tracers)

  integer :: i, j, k, n, kk, ip, jq, kr, tau

  real :: eady_triad, baroclinic_triad, K33, K33_crit
  real :: sumx(0:1), sumy(0:1)

  real :: slope, absslope, depth_ratio
  real :: taperA, taperB
  real :: taper_slope, taper_slope2, gm_taper_slope 

  real :: aredi_scalar, agm_scalar, aredi_plus_agm 
  real :: absdrhodzb(0:1)
  real :: mindelqc(0:1,0:1), delqc_ijk_1, delqc_ijkp1_0

  call mpp_clock_begin(id_clock_fz_terms)

  tau = Time%tau
  
  eady_termx(:,:)       = 0.0
  eady_termy(:,:)       = 0.0
  baroclinic_termx(:,:) = 0.0
  baroclinic_termy(:,:) = 0.0


! Main loop. Short ip,jq,kr loops are explicitly unrolled
! in order to expose independence and allow vectorisation

  do k=1,nk-1
     do j=jsc,jec
        do i=isc,iec

           aredi_scalar    = aredi_array(i,j,k)
           agm_scalar      = agm_array(i,j,k)
           aredi_plus_agm  = aredi_scalar + agm_scalar

           absdrhodzb(0) = abs(drhodzb(i,j,k,0))
           absdrhodzb(1) = abs(drhodzb(i,j,k,1))

           delqc_ijk_1      = delqc(i,j,k,1)
           delqc_ijkp1_0    = delqc(i,j,k+1,0)
           
           !mindelqc(ip,kr) = min(delqc(i-1+ip,j,k+kr,1-kr),delqc(i+ip,j,k+kr,1-kr)) 
           ip=0 ; kr=0
           mindelqc(ip,kr)  = min(delqc(i-1,j,k,1),delqc_ijk_1) 
           ip=0 ; kr=1
           mindelqc(ip,kr)  = min(delqc(i-1,j,k+1,0),delqc_ijkp1_0) 
           ip=1 ; kr=0
           mindelqc(ip,kr)  = min(delqc_ijk_1,delqc(i+1,j,k,1)) 
           ip=1 ; kr=1
           mindelqc(ip,kr)  = min(delqc_ijkp1_0,delqc(i+1,j,k+1,0)) 

           tx_trans_gm(i,j,k) = 0.0
           eady_triad         = 0.0 
           baroclinic_triad   = 0.0

           ip=0   
              sumx(ip) = 0.0

              kr=0
                 slope    = tensor_31(i,j,k,ip,kr)
                 absslope = abs(slope)
                 
                 ! taper for steep slope regions 
                 if(absslope > smax) then 
                     taperA = 0.5*(1.0 + tanh(smax_swidthr-absslope*swidthr))             
                     taperA = dm_taper_const*taperA + gkw_taper_const*(smax/absslope)**2
                 else 
                     taperA = 1.0
                 endif

                 ! taper for grid depths shallower than eddy_depth
                 if(eddy_depth(i,j) >= Grd%zt(k)) then                                           
                     depth_ratio     = Grd%zt(k)/(epsln+eddy_depth(i,j))
                     taperB          = 0.5*(1.0 + sin(pi*(depth_ratio-0.5)))
                 else 
                     taperB          = 1.0
                 endif

                 ! taper times slope for use with neutral diffusion
                 taper_slope  = taperA*taperB*slope
                 taper_slope2 = slope*taper_slope

                 ! taper times slope for use with GM skewsion
                 if(depth_blayer_base(i,j) >= Grd%zt(k)) then      
                     gm_taper_slope = (Grd%zt(k)/depth_blayer_base(i,j))*slope_blayer_base(i,j)*sign(1.0,slope) 
                 else 
                     gm_taper_slope = taper_slope
                 endif


                 tx_trans_gm(i,j,k)     = tx_trans_gm(i,j,k)       + agm_scalar*gm_taper_slope       
                 tensor_31(i,j,k,ip,kr) = aredi_scalar*taper_slope + agm_scalar*gm_taper_slope
                 sumx(ip)               = sumx(ip)                 + mindelqc(ip,kr)*taper_slope2
                 eady_triad             = eady_triad               + absdrhodzb(kr)*taper_slope2
                 baroclinic_triad       = baroclinic_triad         + absdrhodzb(kr)*abs(taper_slope)

              kr=1
                 slope    = tensor_31(i,j,k,ip,kr)
                 absslope = abs(slope)
                           
                 ! taper for steep slope regions 
                 if(absslope > smax) then 
                     taperA = 0.5*(1.0 + tanh(smax_swidthr-absslope*swidthr))             
                     taperA = dm_taper_const*taperA + gkw_taper_const*(smax/absslope)**2
                 else 
                     taperA = 1.0
                 endif

                 ! taper for grid depths shallower than eddy_depth
                 if(eddy_depth(i,j) >= Grd%zt(k)) then                                           
                     depth_ratio     = Grd%zt(k)/(epsln+eddy_depth(i,j))
                     taperB          = 0.5*(1.0 + sin(pi*(depth_ratio-0.5)))
                 else 
                     taperB          = 1.0
                 endif

                 ! taper times slope for use with neutral diffusion
                 taper_slope  = taperA*taperB*slope
                 taper_slope2 = slope*taper_slope

                 ! taper times slope for use with GM skewsion
                 if(depth_blayer_base(i,j) >= Grd%zt(k)) then      
                     gm_taper_slope = (Grd%zt(k)/depth_blayer_base(i,j))*slope_blayer_base(i,j)*sign(1.0,slope) 
                 else 
                     gm_taper_slope = taper_slope
                 endif

                 taper_slope2           = slope*taper_slope
                 tx_trans_gm(i,j,k)     = tx_trans_gm(i,j,k)       + agm_scalar*gm_taper_slope      
                 tensor_31(i,j,k,ip,kr) = aredi_scalar*taper_slope + agm_scalar*gm_taper_slope
                 sumx(ip)               = sumx(ip)                 + mindelqc(ip,kr)*taper_slope2
                 eady_triad             = eady_triad               + absdrhodzb(kr)*taper_slope2
                 baroclinic_triad       = baroclinic_triad         + absdrhodzb(kr)*abs(taper_slope)
                 sumx(ip)               = sumx(ip)*dtew(i,j,ip)

           ip=1   
              sumx(ip) = 0.0

              kr=0
                 slope        = tensor_31(i,j,k,ip,kr)
                 absslope     = abs(slope)
                           
                 ! taper for steep slope regions 
                 if(absslope > smax) then 
                     taperA = 0.5*(1.0 + tanh(smax_swidthr-absslope*swidthr))             
                     taperA = dm_taper_const*taperA + gkw_taper_const*(smax/absslope)**2
                 else 
                     taperA = 1.0
                 endif

                 ! taper for grid depths shallower than eddy_depth
                 if(eddy_depth(i,j) >= Grd%zt(k)) then                                           
                     depth_ratio     = Grd%zt(k)/(epsln+eddy_depth(i,j))
                     taperB          = 0.5*(1.0 + sin(pi*(depth_ratio-0.5)))
                 else 
                     taperB          = 1.0
                 endif

                 ! taper times slope for use with neutral diffusion
                 taper_slope  = taperA*taperB*slope
                 taper_slope2 = slope*taper_slope

                 ! taper times slope for use with GM skewsion
                 if(depth_blayer_base(i,j) >= Grd%zt(k)) then      
                     gm_taper_slope = (Grd%zt(k)/depth_blayer_base(i,j))*slope_blayer_base(i,j)*sign(1.0,slope) 
                 else 
                     gm_taper_slope = taper_slope
                 endif

                 taper_slope2           = slope*taper_slope
                 tx_trans_gm(i,j,k)     = tx_trans_gm(i,j,k)       + agm_scalar*gm_taper_slope    
                 tensor_31(i,j,k,ip,kr) = aredi_scalar*taper_slope + agm_scalar*gm_taper_slope
                 sumx(ip)               = sumx(ip)                 + mindelqc(ip,kr)*taper_slope2
                 eady_triad             = eady_triad               + absdrhodzb(kr)*taper_slope2
                 baroclinic_triad       = baroclinic_triad         + absdrhodzb(kr)*abs(taper_slope)

              kr=1
                 slope        = tensor_31(i,j,k,ip,kr)
                 absslope     = abs(slope)
                                                     
                 ! taper for steep slope regions 
                 if(absslope > smax) then 
                     taperA = 0.5*(1.0 + tanh(smax_swidthr-absslope*swidthr))             
                     taperA = dm_taper_const*taperA + gkw_taper_const*(smax/absslope)**2
                 else 
                     taperA = 1.0
                 endif

                 ! taper for grid depths shallower than eddy_depth
                 if(eddy_depth(i,j) >= Grd%zt(k)) then                                           
                     depth_ratio     = Grd%zt(k)/(epsln+eddy_depth(i,j))
                     taperB          = 0.5*(1.0 + sin(pi*(depth_ratio-0.5)))
                 else 
                     taperB          = 1.0
                 endif

                 ! taper times slope for use with neutral diffusion
                 taper_slope  = taperA*taperB*slope
                 taper_slope2 = slope*taper_slope

                 ! taper times slope for use with GM skewsion
                 if(depth_blayer_base(i,j) >= Grd%zt(k)) then      
                     gm_taper_slope = (Grd%zt(k)/depth_blayer_base(i,j))*slope_blayer_base(i,j)*sign(1.0,slope) 
                 else 
                     gm_taper_slope = taper_slope
                 endif

                 taper_slope2           = slope*taper_slope
                 tx_trans_gm(i,j,k)     = tx_trans_gm(i,j,k)       + agm_scalar*gm_taper_slope   
                 tensor_31(i,j,k,ip,kr) = aredi_scalar*taper_slope + agm_scalar*gm_taper_slope
                 sumx(ip)               = sumx(ip)                 + mindelqc(ip,kr)*taper_slope2
                 eady_triad             = eady_triad               + absdrhodzb(kr)*taper_slope2
                 baroclinic_triad       = baroclinic_triad         + absdrhodzb(kr)*abs(taper_slope)
                 sumx(ip)               = sumx(ip)*dtew(i,j,ip)

           if(Grd%zt(k) >= agm_closure_upper_depth .and. Grd%zt(k) <= agm_closure_lower_depth) then 
             eady_termx(i,j)       = eady_termx(i,j)       + eady_triad*Thickness%dhwt(i,j,k)
             baroclinic_termx(i,j) = baroclinic_termx(i,j) + baroclinic_triad*Thickness%dhwt(i,j,k)
           endif 

           tx_trans_gm(i,j,k)    = 0.25*tx_trans_gm(i,j,k)*Grd%dyu(i,j)

           !mindelqc(jq,kr) = min(delqc(i,j-1+jq,k+kr,1-kr),delqc(i,j+jq,k+kr,1-kr)) 
           jq=0 ; kr=0
           mindelqc(jq,kr)  = min(delqc(i,j-1,k,1),delqc_ijk_1) 
           jq=0 ; kr=1
           mindelqc(jq,kr)  = min(delqc(i,j-1,k+1,0),delqc_ijkp1_0) 
           jq=1 ; kr=0
           mindelqc(jq,kr)  = min(delqc_ijk_1,delqc(i,j+1,k,1)) 
           jq=1 ; kr=1
           mindelqc(jq,kr)  = min(delqc_ijkp1_0,delqc(i,j+1,k+1,0)) 

           ty_trans_gm(i,j,k) = 0.0
           eady_triad         = 0.0 
           baroclinic_triad   = 0.0 

           jq=0   
              sumy(jq) = 0.0

               kr=0
                 slope        = tensor_32(i,j,k,jq,kr)
                 absslope     = abs(slope)
                           
                 ! taper for steep slope regions 
                 if(absslope > smax) then 
                     taperA = 0.5*(1.0 + tanh(smax_swidthr-absslope*swidthr))             
                     taperA = dm_taper_const*taperA + gkw_taper_const*(smax/absslope)**2
                 else 
                     taperA = 1.0
                 endif

                 ! taper for grid depths shallower than eddy_depth
                 if(eddy_depth(i,j) >= Grd%zt(k)) then                                           
                     depth_ratio     = Grd%zt(k)/(epsln+eddy_depth(i,j))
                     taperB          = 0.5*(1.0 + sin(pi*(depth_ratio-0.5)))
                 else 
                     taperB          = 1.0
                 endif

                 ! taper times slope for use with neutral diffusion
                 taper_slope  = taperA*taperB*slope
                 taper_slope2 = slope*taper_slope

                 ! taper times slope for use with GM skewsion
                 if(depth_blayer_base(i,j) >= Grd%zt(k)) then      
                     gm_taper_slope = (Grd%zt(k)/depth_blayer_base(i,j))*slope_blayer_base(i,j)*sign(1.0,slope) 
                 else 
                     gm_taper_slope = taper_slope
                 endif

                 taper_slope2           = slope*taper_slope
                 ty_trans_gm(i,j,k)     = ty_trans_gm(i,j,k)       + agm_scalar*gm_taper_slope       
                 tensor_32(i,j,k,jq,kr) = aredi_scalar*taper_slope + agm_scalar*gm_taper_slope
                 sumy(jq)               = sumy(jq)                 + mindelqc(jq,kr)*taper_slope2
                 eady_triad             = eady_triad               + absdrhodzb(kr)*taper_slope2
                 baroclinic_triad       = baroclinic_triad         + absdrhodzb(kr)*abs(taper_slope)

              kr=1
                 slope        = tensor_32(i,j,k,jq,kr)
                 absslope     = abs(slope)
                           
                 ! taper for steep slope regions 
                 if(absslope > smax) then 
                     taperA = 0.5*(1.0 + tanh(smax_swidthr-absslope*swidthr))             
                     taperA = dm_taper_const*taperA + gkw_taper_const*(smax/absslope)**2
                 else 
                     taperA = 1.0
                 endif

                 ! taper for grid depths shallower than eddy_depth
                 if(eddy_depth(i,j) >= Grd%zt(k)) then                                           
                     depth_ratio     = Grd%zt(k)/(epsln+eddy_depth(i,j))
                     taperB          = 0.5*(1.0 + sin(pi*(depth_ratio-0.5)))
                 else 
                     taperB          = 1.0
                 endif

                 ! taper times slope for use with neutral diffusion
                 taper_slope  = taperA*taperB*slope
                 taper_slope2 = slope*taper_slope

                 ! taper times slope for use with GM skewsion
                 if(depth_blayer_base(i,j) >= Grd%zt(k)) then      
                     gm_taper_slope = (Grd%zt(k)/depth_blayer_base(i,j))*slope_blayer_base(i,j)*sign(1.0,slope) 
                 else 
                     gm_taper_slope = taper_slope
                 endif

                 taper_slope2           = slope*taper_slope
                 ty_trans_gm(i,j,k)     = ty_trans_gm(i,j,k)       + agm_scalar*gm_taper_slope      
                 tensor_32(i,j,k,jq,kr) = aredi_scalar*taper_slope + agm_scalar*gm_taper_slope
                 sumy(jq)               = sumy(jq)                 + mindelqc(jq,kr)*taper_slope2
                 eady_triad             = eady_triad               + absdrhodzb(kr)*taper_slope2
                 baroclinic_triad       = baroclinic_triad         + absdrhodzb(kr)*abs(taper_slope)
                 sumy(jq)               = sumy(jq)*dtns(i,j,jq)

           jq=1   
              sumy(jq) = 0.0

               kr=0
                 slope    = tensor_32(i,j,k,jq,kr)
                 absslope = abs(slope)

                 ! taper for steep slope regions 
                 if(absslope > smax) then 
                     taperA = 0.5*(1.0 + tanh(smax_swidthr-absslope*swidthr))             
                     taperA = dm_taper_const*taperA + gkw_taper_const*(smax/absslope)**2
                 else 
                     taperA = 1.0
                 endif

                 ! taper for grid depths shallower than eddy_depth
                 if(eddy_depth(i,j) >= Grd%zt(k)) then                                           
                     depth_ratio     = Grd%zt(k)/(epsln+eddy_depth(i,j))
                     taperB          = 0.5*(1.0 + sin(pi*(depth_ratio-0.5)))
                 else 
                     taperB          = 1.0
                 endif

                 ! taper times slope for use with neutral diffusion
                 taper_slope  = taperA*taperB*slope
                 taper_slope2 = slope*taper_slope

                 ! taper times slope for use with GM skewsion
                 if(depth_blayer_base(i,j) >= Grd%zt(k)) then      
                     gm_taper_slope = (Grd%zt(k)/depth_blayer_base(i,j))*slope_blayer_base(i,j)*sign(1.0,slope) 
                 else 
                     gm_taper_slope = taper_slope
                 endif

                 taper_slope2           = slope*taper_slope
                 ty_trans_gm(i,j,k)     = ty_trans_gm(i,j,k)       + agm_scalar*gm_taper_slope     
                 tensor_32(i,j,k,jq,kr) = aredi_scalar*taper_slope + agm_scalar*gm_taper_slope
                 sumy(jq)               = sumy(jq)                 + mindelqc(jq,kr)*taper_slope2
                 eady_triad             = eady_triad               + absdrhodzb(kr)*taper_slope2
                 baroclinic_triad       = baroclinic_triad         + absdrhodzb(kr)*abs(taper_slope)

              kr=1
                 slope    = tensor_32(i,j,k,jq,kr)
                 absslope = abs(slope)

                 ! taper for steep slope regions 
                 if(absslope > smax) then 
                     taperA = 0.5*(1.0 + tanh(smax_swidthr-absslope*swidthr))             
                     taperA = dm_taper_const*taperA + gkw_taper_const*(smax/absslope)**2
                 else 
                     taperA = 1.0
                 endif

                 ! taper for grid depths shallower than eddy_depth
                 if(eddy_depth(i,j) >= Grd%zt(k)) then                                           
                     depth_ratio     = Grd%zt(k)/(epsln+eddy_depth(i,j))
                     taperB          = 0.5*(1.0 + sin(pi*(depth_ratio-0.5)))
                 else 
                     taperB          = 1.0
                 endif

                 ! taper times slope for use with neutral diffusion
                 taper_slope  = taperA*taperB*slope
                 taper_slope2 = slope*taper_slope

                 ! taper times slope for use with GM skewsion
                 if(depth_blayer_base(i,j) >= Grd%zt(k)) then      
                     gm_taper_slope = (Grd%zt(k)/depth_blayer_base(i,j))*slope_blayer_base(i,j)*sign(1.0,slope) 
                 else 
                     gm_taper_slope = taper_slope
                 endif

                 taper_slope2           = slope*taper_slope
                 ty_trans_gm(i,j,k)     = ty_trans_gm(i,j,k)       + agm_scalar*gm_taper_slope      
                 tensor_32(i,j,k,jq,kr) = aredi_scalar*taper_slope + agm_scalar*gm_taper_slope
                 sumy(jq)               = sumy(jq)                 + mindelqc(jq,kr)*taper_slope2
                 eady_triad             = eady_triad               + absdrhodzb(kr)*taper_slope2
                 baroclinic_triad       = baroclinic_triad         + absdrhodzb(kr)*abs(taper_slope)

                 sumy(jq) = sumy(jq)*dtns(i,j,jq)

           if(Grd%zt(k) >= agm_closure_upper_depth .and. Grd%zt(k) <= agm_closure_lower_depth) then
             eady_termy(i,j)       = eady_termy(i,j)       + eady_triad*Thickness%dhwt(i,j,k)
             baroclinic_termy(i,j) = baroclinic_termy(i,j) + baroclinic_triad*Thickness%dhwt(i,j,k)
           endif 

           ty_trans_gm(i,j,k)    = 0.25*ty_trans_gm(i,j,k)*Grd%dxu(i,j)
           
           K33 = aredi_scalar*Grd%tmask(i,j,k+1)*(Grd%dxtr(i,j)*(sumx(0)+sumx(1)) &
                            + Grd%dytr(i,j)*(sumy(0)+sumy(1)))*dhwtr(i,j,k)

           ! determine part of K33 included explicitly in time and that part implicitly. 
           ! explicit calculation is more accurate than implicit, so aim to compute as much 
           ! as stably possible via the explicit method. 
           K33_explicit(i,j,k) = K33
           K33_implicit(i,j,k) = 0.0 
           if(.not. diffusion_all_explicit) then 
               K33_crit = two_dtime_inv*Thickness%dht(i,j,k,tau)**2
               if(K33 >= K33_crit) then 
                   K33_explicit(i,j,k) = K33_crit
                   K33_implicit(i,j,k) = K33-K33_crit
               endif
           endif

        enddo
     enddo
  enddo

  if(tmask_sigma_on) then
     do j=jsc,jec
        do i=isc,iec
           if(tmask_sigma(i,j) > 0.0) then 
              k = Grd%kmt(i,j)-1
              K33_implicit(i,j,k) = K33_implicit(i,j,k)*(1.0-tmask_sigma(i,j))
              K33_explicit(i,j,k) = K33_explicit(i,j,k)*(1.0-tmask_sigma(i,j))
           endif
        enddo
     enddo
 endif
 
 if(tmask_neutral_on) then
     do j=jsc,jec
        do i=isc,iec
           K33_implicit(i,j,1) = 0.0  
           K33_explicit(i,j,1) = 0.0  
           tx_trans_gm(i,j,1)  = 0.0  
           ty_trans_gm(i,j,1)  = 0.0  
           if(Grd%kmt(i,j) > 1) then 
               k = Grd%kmt(i,j)-1
               K33_implicit(i,j,k) = 0.0
               K33_explicit(i,j,k) = 0.0
               tx_trans_gm(i,j,k)  = 0.0  
               ty_trans_gm(i,j,k)  = 0.0  
           endif
        enddo
     enddo
 endif

 do n=1,num_prog_tracers 
   T_prog(n)%K33_implicit(:,:,:) = K33_implicit(:,:,:)
 enddo
 if(neutral_physics_limit) then 
     do n=1,num_prog_tracers 
        do k=1,nk
           do j=jsc,jec
              do i=isc,iec
                if(T_prog(n)%tmask_limit(i,j,k)==1.0) T_prog(n)%K33_implicit(i,j,k) = 0.0
              enddo
           enddo
        enddo
     enddo
 endif


  do n=1,num_prog_tracers
     if (id_k33_implicit(n) > 0) then
        used = send_data (id_k33_implicit(n), T_prog(n)%K33_implicit(:,:,:), &
               Time%model_time, rmask=Grd%tmask(:,:,:), &
               is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
     endif 
  enddo 


  call mpp_clock_end(id_clock_fz_terms)

end subroutine fz_terms
! </SUBROUTINE>  NAME="fz_terms"


!#######################################################################
! <SUBROUTINE NAME="fx_flux">
!
! <DESCRIPTION>
! Subroutine computes the zonal neutral physics tracer flux component.
! Compute this component for all tracers at level k.
! fx has physical dimensions (area*diffusivity*tracer gradient)
! </DESCRIPTION>
!
subroutine fx_flux(Time, Thickness, T_prog, k)

  type(ocean_time_type), intent(in)        :: Time
  type(ocean_thickness_type), intent(in)   :: Thickness
  type(ocean_prog_tracer_type), intent(in) :: T_prog(num_prog_tracers)
  integer, intent(in) :: k

  integer :: nn, i, j, ip, tau
  integer :: kk, kr, kpkr
  real    :: tensor_11(isd:ied,jsd:jed,0:1,0:1)
  real    :: tensor_13(isd:ied,jsd:jed,0:1,0:1)
  real    :: sumz(isd:ied,jsd:jed,0:1,num_prog_tracers)
  real    :: taperA, taperB
  real    :: slope, absslope, taper_slope, gm_taper_slope
  real    :: depth_ratio

  call mpp_clock_begin(id_clock_fx_flux)

  tau = Time%tau

  ! initialize arrays to zero 
  tensor_11(:,:,:,:) = 0.0
  tensor_13(:,:,:,:) = 0.0

  ! tracer-independent part of the calculation 
  do kr=0,1
     kpkr = min(k+kr,nk)
     do ip=0,1 
        do j=jsc,jec
           do i=isc-1,iec

              slope = -Grd%tmask(i+ip,j,kpkr)*&
                   (drhodT(i+ip,j,k)*dTdx(index_temp)%field(i,j,k) + &
                   drhodS(i+ip,j,k)*dTdx(index_salt)%field(i,j,k)) &
                   /(drhodT(i+ip,j,k)*dTdz(index_temp)%field(i+ip,j,k-1+kr) +&
                   drhodS(i+ip,j,k)*dTdz(index_salt)%field(i+ip,j,k-1+kr) -epsln)
              absslope = abs(slope)

              ! taper for steep slope regions 
              if(absslope > smax) then 
                  taperA = 0.5*(1.0 + tanh(smax_swidthr-absslope*swidthr))             
                  taperA = dm_taper_const*taperA + gkw_taper_const*(smax/absslope)**2
              else 
                  taperA = 1.0
              endif

              ! taper for grid depths shallower than eddy_depth
              if(eddy_depth(i,j) >= Grd%zt(k)) then                                           
                  depth_ratio     = Grd%zt(k)/(epsln+eddy_depth(i,j))
                  taperB          = 0.5*(1.0 + sin(pi*(depth_ratio-0.5)))
              else 
                  taperB          = 1.0
              endif

              ! taper times slope for use with neutral diffusion
              taper_slope  = taperA*taperB*slope

              ! taper times slope for use with GM skewsion
              if(depth_blayer_base(i,j) >= Grd%zt(k)) then      
                  gm_taper_slope = (Grd%zt(k)/depth_blayer_base(i,j))*slope_blayer_base(i,j)*sign(1.0,slope) 
              else 
                  gm_taper_slope = taper_slope
              endif

              tensor_11(i,j,ip,kr) = aredi_array(i,j,k)*( (1.0-taper_diagonal) + taperA*taperB*taper_diagonal)  &
                                   + ah_array(i,j)*(1.0-taperB)
              tensor_13(i,j,ip,kr) = aredi_array(i,j,k)*taper_slope-agm_array(i,j,k)*gm_taper_slope

           enddo
        enddo
     enddo
  enddo

  ! tracer-dependent part of the calculation
  do nn=1,num_prog_tracers
     sumz(:,:,:,:) = 0.0
     do kr=0,1
        do ip=0,1
           do j=jsc,jec
              do i=isc-1,iec
                sumz(i,j,kr,nn) = sumz(i,j,kr,nn) + dtwedyt(i,j,ip)*Grd%tmask(i,j,k)*Grd%tmask(i+1,j,k) &
                *(tensor_11(i,j,ip,kr)*dTdx(nn)%field(i,j,k) + tensor_13(i,j,ip,kr)*dTdz(nn)%field(i+ip,j,k-1+kr)) &
                * min(delqc(i,j,k,kr),delqc(i+1,j,k,kr))
              enddo
           enddo
        enddo
     enddo
     do j=jsc,jec
        do i=isc-1,iec
           flux_x(nn)%field(i,j,k) = Grd%dxter(i,j)*(sumz(i,j,0,nn)+sumz(i,j,1,nn))
        enddo
     enddo
  enddo

  ! apply some masks 
  if(tmask_neutral_on) then
     do nn=1,num_prog_tracers
        do j=jsc,jec
           do i=isc-1,iec
              if(k==Grd%kmt(i,j)) then
                flux_x(nn)%field(i,j,k) = aredi_array(i,j,k)*dTdx(nn)%field(i,j,k)*Grd%dyte(i,j)* &
                                          min(Thickness%dht(i,j,k,tau),Thickness%dht(i+1,j,k,tau))
              endif
           enddo
        enddo
     enddo
     if(k==1) then
         do nn=1,num_prog_tracers
            flux_x(nn)%field(:,:,k) = aredi_array(:,:,k)*dTdx(nn)%field(:,:,k)*Grd%dyte(:,:)*FMX(Thickness%dht(:,:,k,tau))
         enddo
     endif
  endif
  if(tmask_sigma_on) then
     do nn=1,num_prog_tracers
        do j=jsc,jec
           do i=isc-1,iec
              if(k==Grd%kmt(i,j)) flux_x(nn)%field(i,j,k) = flux_x(nn)%field(i,j,k)*(1.0-tmask_sigma(i,j))
           enddo
        enddo
     enddo
 endif

 if(neutral_physics_limit) then 
     do nn=1,num_prog_tracers
        do j=jsc,jec
           do i=isc-1,iec
              if(T_prog(nn)%tmask_limit(i,j,k)==1.0) then 
                  flux_x(nn)%field(i,j,k) = aredi_array(i,j,k)*dTdx(nn)%field(i,j,k)*Grd%dyte(i,j) &
                                            *min(Thickness%dht(i,j,k,tau),Thickness%dht(i+1,j,k,tau))
              endif
           enddo
        enddo
     enddo
 endif

 call mpp_clock_end(id_clock_fx_flux) 

end subroutine fx_flux
! </SUBROUTINE> NAME="fx_flux"



!#######################################################################
! <SUBROUTINE NAME="fy_flux">
!
! <DESCRIPTION>
! Subroutine computes the meridional neutral physics tracer flux component.
! Compute this component for all tracers at level k.
! fy has physical dimensions (area*diffusivity*tracer gradient)
! </DESCRIPTION>
!
subroutine fy_flux(Time, Thickness, T_prog, k)

  type(ocean_time_type), intent(in)        :: Time
  type(ocean_thickness_type), intent(in)   :: Thickness
  type(ocean_prog_tracer_type), intent(in) :: T_prog(num_prog_tracers)
  integer, intent(in)                      :: k

  integer :: nn, i, j, jq, tau
  integer :: kk, kr, kpkr
  real :: tensor_22(isd:ied,jsd:jed,0:1,0:1)
  real :: tensor_23(isd:ied,jsd:jed,0:1,0:1)
  real :: sumz(isd:ied,jsd:jed,0:1,num_prog_tracers)
  real :: taperA, taperB
  real :: slope, absslope, taper_slope, gm_taper_slope
  real :: depth_ratio

  call mpp_clock_begin(id_clock_fy_flux)

  tau = Time%tau 

  ! initialize arrays to zero 
  tensor_22(:,:,:,:) = 0.0
  tensor_23(:,:,:,:) = 0.0

  ! tracer-independent part of the calculation 
  do kr=0,1
     kpkr = min(k+kr,nk)
     do jq=0,1  
       do j=jsc-1,jec
          do i=isc,iec
 
              slope = -Grd%tmask(i,j+jq,kpkr)&
                     *(drhodT(i,j+jq,k)*dTdy(index_temp)%field(i,j,k) + &
                     drhodS(i,j+jq,k)*dTdy(index_salt)%field(i,j,k)) & 
                     /(drhodT(i,j+jq,k)*dTdz(index_temp)%field(i,j+jq,k-1+kr) + &
                     drhodS(i,j+jq,k)*dTdz(index_salt)%field(i,j+jq,k-1+kr) - epsln)  
              absslope = abs(slope)  

              ! taper for steep slope regions 
              if(absslope > smax) then 
                  taperA = 0.5*(1.0 + tanh(smax_swidthr-absslope*swidthr))             
                  taperA = dm_taper_const*taperA + gkw_taper_const*(smax/absslope)**2
              else 
                  taperA = 1.0
              endif

              ! taper for grid depths shallower than eddy_depth
              if(eddy_depth(i,j) >= Grd%zt(k)) then                                           
                  depth_ratio     = Grd%zt(k)/(epsln+eddy_depth(i,j))
                  taperB          = 0.5*(1.0 + sin(pi*(depth_ratio-0.5)))
              else 
                  taperB          = 1.0
              endif

              ! taper times slope for use with neutral diffusion
              taper_slope  = taperA*taperB*slope

              ! taper times slope for use with GM skewsion
              if(depth_blayer_base(i,j) >= Grd%zt(k)) then      
                  gm_taper_slope = (Grd%zt(k)/depth_blayer_base(i,j))*slope_blayer_base(i,j)*sign(1.0,slope) 
              else 
                  gm_taper_slope = taper_slope
              endif

              tensor_22(i,j,jq,kr) = aredi_array(i,j,k)*( (1.0-taper_diagonal) + taperA*taperB*taper_diagonal) &
                                   + ah_array(i,j)*(1.0-taperB)
              tensor_23(i,j,jq,kr) = aredi_array(i,j,k)*taper_slope-agm_array(i,j,k)*gm_taper_slope

           enddo
        enddo
     enddo
  enddo

  ! tracer-dependent part of the calculation
  do nn=1,num_prog_tracers
     sumz(:,:,:,:) = 0.0
     do kr=0,1
        do jq=0,1
           do j=jsc-1,jec
              do i=isc,iec
                 sumz(i,j,kr,nn) = sumz(i,j,kr,nn) + dxtdtsn(i,j,jq)*Grd%tmask(i,j,k)*Grd%tmask(i,j+1,k) &
                      *(tensor_22(i,j,jq,kr)*dTdy(nn)%field(i,j,k) + tensor_23(i,j,jq,kr)*dTdz(nn)%field(i,j+jq,k-1+kr)) &
                      * min(delqc(i,j,k,kr),delqc(i,j+1,k,kr))
              enddo
           enddo
        enddo
     enddo
     do j=jsc-1,jec
        do i=isc,iec
           flux_y(nn)%field(i,j,k) = Grd%dytnr(i,j)*(sumz(i,j,0,nn)+sumz(i,j,1,nn))
        enddo
     enddo
  enddo


  ! apply some masks 
  if(tmask_neutral_on) then
     do nn=1,num_prog_tracers
        do j=jsc-1,jec
           do i=isc,iec
              if(k==Grd%kmt(i,j)) then
                flux_y(nn)%field(i,j,k) = aredi_array(i,j,k)*dTdy(nn)%field(i,j,k)*Grd%dxtn(i,j) &
                                          *min(Thickness%dht(i,j,k,tau),Thickness%dht(i,j+1,k,tau))
              endif
           enddo
        enddo
     enddo
     if(k==1) then
         do nn=1,num_prog_tracers
            flux_y(nn)%field(:,:,k) = aredi_array(:,:,k)*dTdy(nn)%field(:,:,k)*Grd%dxtn(:,:)*FMY(Thickness%dht(:,:,k,tau))
         enddo
     endif
  endif
  if(tmask_sigma_on) then
      do nn=1,num_prog_tracers
         do j=jsc-1,jec
            do i=isc,iec
               if(k==Grd%kmt(i,j)) flux_y(nn)%field(i,j,k) = flux_y(nn)%field(i,j,k)*(1.0-tmask_sigma(i,j))
            enddo
         enddo
      enddo
  endif

  if(neutral_physics_limit) then 
      do nn=1,num_prog_tracers
         do j=jsc-1,jec
            do i=isc,iec
               if(T_prog(nn)%tmask_limit(i,j,k)==1.0) then 
                   flux_y(nn)%field(i,j,k) = aredi_array(i,j,k)*dTdy(nn)%field(i,j,k)*Grd%dxtn(i,j) &
                                             *min(Thickness%dht(i,j,k,tau),Thickness%dht(i,j+1,k,tau))
               endif
            enddo
         enddo
      enddo
  endif

  call mpp_clock_end(id_clock_fy_flux)

end subroutine fy_flux
! </SUBROUTINE> NAME="fy_flux"


!#######################################################################
! <SUBROUTINE NAME="fz_flux">
!
! <DESCRIPTION>
! Subroutine computes the vertical neutral physics tracer flux component.
! Compute this component for all tracers at level k.
! Surface and bottom boundary condition fz(k=0)=fz(k=kmt(i,j))=0
! fz has physical dimensions (diffusivity*tracer gradient) (usual flux units)
! </DESCRIPTION>
!
subroutine fz_flux(T_prog, k)

  type(ocean_prog_tracer_type), intent(in) :: T_prog(num_prog_tracers)
  integer, intent(in) :: k

  integer :: nn, i, j, ip, jq, kr
  real :: sumx_0, sumx_1, sumy_0, sumy_1
  real :: temparray31(isc:iec,jsc:jec,0:1,0:1)
  real :: temparray32(isc:iec,jsc:jec,0:1,0:1)
  real :: dxtr_ij, dytr_ij, dhwtr_ijk
  real :: tmaskkp1
 
  if(tmask_neutral_on .and. k==1) then
     do nn = 1, num_prog_tracers
       fz2(nn)%field(:,:) = 0.0 
      enddo
     return 
  endif

  if(k==nk) then 

     do nn = 1, num_prog_tracers
        do j=jsc,jec
           do i=isc,iec
              fz2(nn)%field(i,j) = 0.0
           enddo
        enddo
     enddo
     return  

  elseif(k < nk) then  

      ! tracer-independent part of the calculation 
      do kr=0,1
         do ip=0,1
            do j=jsc,jec
               do i=isc,iec
                  temparray31(i,j,ip,kr) = tensor_31(i,j,k,ip,kr)*dtew(i,j,ip) &
                       *min(delqc(i-1+ip,j,k+kr,1-kr),delqc(i+ip,j,k+kr,1-kr))
               enddo
            enddo
         enddo
         do jq=0,1
            do j=jsc,jec
               do i=isc,iec
                  temparray32(i,j,jq,kr) = tensor_32(i,j,k,jq,kr)*dtns(i,j,jq) &
                       *min(delqc(i,j-1+jq,k+kr,1-kr),delqc(i,j+jq,k+kr,1-kr))  
               enddo
            enddo
         enddo
      enddo

      ! tracer-dependent part of the calculation  
      do nn=1,num_prog_tracers
         do j=jsc,jec
            do i=isc,iec
               sumx_0 =  temparray31(i,j,0,0)*dTdx(nn)%field(i-1,j,k) &
                      +  temparray31(i,j,0,1)*dTdx(nn)%field(i-1,j,k+1)
               sumx_1 =  temparray31(i,j,1,0)*dTdx(nn)%field(i,j,k) &
                      +  temparray31(i,j,1,1)*dTdx(nn)%field(i,j,k+1)
               sumy_0 =  temparray32(i,j,0,0)*dTdy(nn)%field(i,j-1,k) &
                      +  temparray32(i,j,0,1)*dTdy(nn)%field(i,j-1,k+1)
               sumy_1 =  temparray32(i,j,1,0)*dTdy(nn)%field(i,j,k) &
                      +  temparray32(i,j,1,1)*dTdy(nn)%field(i,j,k+1)

               ! compute time explicit portion of the vertical flux
               fz2(nn)%field(i,j) =   Grd%tmask(i,j,k+1) &
                    *( Grd%dxtr(i,j)*(sumx_0+sumx_1) &
                      +Grd%dytr(i,j)*(sumy_0+sumy_1)) &
                    *dhwtr(i,j,k) &
                    + K33_explicit(i,j,k)*dTdz(nn)%field(i,j,k)

            enddo
         enddo
      enddo
      
      if(tmask_neutral_on) then
          do nn=1,num_prog_tracers
             do j=jsc,jec
                do i=isc,iec
                   if(k==(Grd%kmt(i,j)-1))  fz2(nn)%field(i,j) = 0.0
                enddo
             enddo
          enddo
      endif

      if(tmask_sigma_on) then
          do nn=1,num_prog_tracers
             do j=jsc,jec
                do i=isc,iec
                   if(tmask_sigma(i,j) > 0.0 .and. k==(Grd%kmt(i,j)-1) ) then
                       fz2(nn)%field(i,j) = fz2(nn)%field(i,j)*(1.0-tmask_sigma(i,j))
                    endif
                enddo
             enddo
          enddo
      endif
 
      if(neutral_physics_limit) then 
          do nn=1,num_prog_tracers
             do j=jsc,jec
                do i=isc,iec
                   if(T_prog(nn)%tmask_limit(i,j,k)==1.0) fz2(nn)%field(i,j) = 0.0  
                enddo
             enddo
          enddo
      endif

  endif  !if-test for k-level 


end subroutine fz_flux
! </SUBROUTINE> NAME="fz_flux"


!#######################################################################
! <SUBROUTINE NAME="compute_eady_rate">
!
! <DESCRIPTION>
! Finish computing eady growth rate.
! </DESCRIPTION>
!
subroutine compute_eady_rate(model_time)

  type(time_type), intent(in) :: model_time
  integer :: i, j
  real    :: thickness 

  call mpp_clock_begin(id_clock_compute_eady_rate)
  
  thickness = epsln + agm_closure_lower_depth - agm_closure_upper_depth

  do j=jsc,jec
     do i=isc,iec
        eady_rate(i,j) = Grd%tmask(i,j,1)*sqrt((eady_termx(i,j)*count_x(i,j)+eady_termy(i,j)* &
             count_y(i,j))*gravrho0r/thickness)
     enddo
  enddo

  if (id_eady > 0) used = send_data (id_eady, eady_rate(:,:), &
                          model_time, rmask=Grd%tmask(:,:,1), &
                          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

  call mpp_clock_end(id_clock_compute_eady_rate)

end subroutine compute_eady_rate
! </SUBROUTINE> NAME="compute_eady_rate"



!#######################################################################
! <SUBROUTINE NAME="compute_baroclinicity">
!
! <DESCRIPTION>
! Finish computing baroclinicity, which is defined to be the vertically
! averaged magnitude of the horizontal density gradient.
! </DESCRIPTION>
!
subroutine compute_baroclinicity(model_time)

  type(time_type), intent(in) :: model_time
  integer :: i, j
  real    :: thickness 

  call mpp_clock_begin(id_clock_compute_baroclinicity)
  
  thickness = epsln + agm_closure_lower_depth - agm_closure_upper_depth

  do j=jsc,jec
     do i=isc,iec
        baroclinicity(i,j) = Grd%tmask(i,j,1)* &
             (baroclinic_termx(i,j)*count_x(i,j)+baroclinic_termy(i,j)*count_y(i,j)) &
             /thickness
     enddo
  enddo

  if (id_baroclinicity > 0) used = send_data (id_baroclinicity, baroclinicity(:,:), &
                                   model_time, rmask=Grd%tmask(:,:,1), &
                                   is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

  call mpp_clock_end(id_clock_compute_baroclinicity)

end subroutine compute_baroclinicity
! </SUBROUTINE> NAME="compute_baroclinicity"


!#######################################################################
! <SUBROUTINE NAME="compute_rossby_radius">
!
! <DESCRIPTION>
! Subroutine computes the first baroclinic Rossby radius of deformation. 
! Employ WKB approach described by Chelton et al.  In particular, 
! use formulae (2.2), (2.3a) and (2.3b) from their paper. 
!
! Place a max and min value on the Rossby radius, as such is useful
! for computing the diffusivity as well as for the sine_taper scheme.  
!
! Compute buoyancy frequency in terms of vertical gradient of 
! locally referenced potential density.  Place the reference point
! at the interface between the tracer cells, which is also where 
! the vertical derivative of neutral density is located.  This amounts 
! to a centered difference computation similar to that used by 
! Chelton et al. equation (B.4). 
! </DESCRIPTION>
!
subroutine compute_rossby_radius(Thickness, model_time)

  type(ocean_thickness_type), intent(in) :: Thickness
  type(time_type), intent(in)            :: model_time
  real, dimension(isc:iec,jsc:jec)       :: gravity_wave_speed 
  real                                   :: drhodzb_speed
  real                                   :: rossby_less5
  integer                                :: i,j,k

  call mpp_clock_begin(id_clock_compute_rossby_radius)
  
  gravity_wave_speed(:,:) = 0.0
  do k=1,nk-1
     do j=jsc,jec
        do i=isc,iec
           drhodzb_speed  = 0.5*( (drhodT(i,j,k) + drhodT(i,j,k+1))*dTdz(index_temp)%field(i,j,k) &
                +(drhodS(i,j,k) + drhodS(i,j,k+1))*dTdz(index_salt)%field(i,j,k) )
           gravity_wave_speed(i,j) = gravity_wave_speed(i,j) + Thickness%dhwt(i,j,k)*sqrt(abs(gravrho0r*drhodzb_speed))
        enddo
     enddo
  enddo

  do j=jsc,jec
     do i=isc,iec
        rossby_radius(i,j) = gravity_wave_speed(i,j)/(pi*coriolis_param(i,j) + epsln)
        if(abs(Grd%phit(i,j)) < fivedeg) then 
            rossby_less5 = sqrt(gravity_wave_speed(i,j)/(2.0*beta_param(i,j)))
            rossby_radius(i,j) = (2.0*rossby_less5*rossby_radius(i,j))/(rossby_less5+rossby_radius(i,j)+epsln)
        endif

     enddo
  enddo

  do j=jsc,jec
     do i=isc,iec
        rossby_radius(i,j) = min(rossby_radius_max,max(rossby_radius_min,rossby_radius(i,j)))
     enddo
  enddo

  if (id_rossby > 0) used = send_data (id_rossby, rossby_radius(:,:), &
                            model_time, rmask=Grd%tmask(:,:,1), &
                            is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

  call mpp_clock_end(id_clock_compute_rossby_radius)

end subroutine compute_rossby_radius
! </SUBROUTINE> NAME="compute_rossby_radius"


!#######################################################################
! <SUBROUTINE NAME="compute_bczone_radius">
!
! <DESCRIPTION>
! Subroutine computes the radius of the baroclinic zone in a manner 
! suggested by the Hadley Centre approach (Malcolm Roberts, personal 
! communication).  Algorithm was used in MOM3.
! </DESCRIPTION>
!
subroutine compute_bczone_radius(model_time)

  type(time_type), intent(in) :: model_time
  integer :: i,j,k
  integer :: ii, jj
  real    :: n_zone, e_zone, s_zone, w_zone, fract, nstot, ewtot
  integer :: ip, jq

  call mpp_clock_begin(id_clock_compute_bczone_radius)

  bczone_rate(isc:iec,jsc:jec) = eady_rate(isc:iec,jsc:jec)
  call mpp_update_domains (bczone_rate(:,:), BCzone_domain%domain2d)

  do j=jsc,jec
    do i=isc,iec
      if (bczone_rate(i,j) > agm_closure_bczone_crit_rate) then

        ! search northward 
        n_zone = bczone_dyt(i,j)
        do jq=j+1,j+bczone_max_pts
          if (bczone_rate(i,jq) > agm_closure_bczone_crit_rate) then
            n_zone = n_zone + bczone_dyt(i,jq)
          else
            exit
          endif
        enddo

        ! search southward 
        s_zone = bczone_dyt(i,j)
        do jq=j-1,j-bczone_max_pts,-1
          if (bczone_rate(i,jq) > agm_closure_bczone_crit_rate) then
            s_zone = s_zone + bczone_dyt(i,jq)
          else
            exit
          endif
        enddo

        ! search eastward 
        e_zone = bczone_dxt(i,j)
        do ip=i+1,i+bczone_max_pts
          if (bczone_rate(ip,j) > agm_closure_bczone_crit_rate) then
            e_zone = e_zone + bczone_dxt(ip,j)
          else
            exit
          endif
        enddo

        ! search westward 
        w_zone = bczone_dxt(i,j)
        do ip=i-1,i-bczone_max_pts,-1
          if (bczone_rate(ip,j) > agm_closure_bczone_crit_rate) then
            w_zone = w_zone + bczone_dxt(ip,j)
          else
            exit
          endif
        enddo

        ! total radius (subtraction accounts for double-counting central point)
        nstot=n_zone+s_zone-bczone_dyt(i,j)
        ewtot=e_zone+w_zone-bczone_dxt(i,j)

        if (nstot < ewtot) then 
          fract=min(n_zone,s_zone)/max(n_zone,s_zone) 
          bczone_radius(i,j)=fract*nstot 
        else   
          fract=min(e_zone,w_zone)/max(e_zone,w_zone)
          bczone_radius(i,j)=fract*ewtot
        endif
      endif
          
    enddo
  enddo

  if (id_bczone > 0) used = send_data (id_bczone, bczone_radius(:,:), &
                            model_time, rmask=Grd%tmask(:,:,1), &
                            is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

  call mpp_clock_end(id_clock_compute_bczone_radius)

end subroutine compute_bczone_radius
! </SUBROUTINE> NAME="compute_bczone_radius"


!#######################################################################
! <SUBROUTINE NAME="compute_diffusivity">
!
! <DESCRIPTION>
! Subroutine computes flow dependent diffusivity.
! Allow for an added dimensionless tuning factor as well as a 
! minimum and maximum diffusivity. 
! </DESCRIPTION>
!
subroutine compute_diffusivity(model_time)

  type(time_type), intent(in) :: model_time
  integer :: i, j, k
  real    :: agm_growth_rate(isd:ied,jsd:jed)
  real    :: agm_growth_rate_max(isd:ied,jsd:jed) 
  real    :: agm_length(isd:ied,jsd:jed)
  real    :: agm_qg(isd:ied,jsd:jed)
  real    :: growth_rate(isd:ied,jsd:jed)

  call mpp_clock_begin(id_clock_compute_diffusivity)

  agm_growth_rate     = 0.0
  agm_growth_rate_max = 0.0
  agm_length          = 0.0
  agm_qg              = 0.0
  growth_rate         = 0.0

  if(agm_closure_length_fixed .or. agm_closure_baroclinic) then 
      do j=jsc,jec
         do i=isc,iec
            agm_length(i,j) = agm_closure_length
         enddo
      enddo
  elseif(agm_closure_length_rossby) then 
      do j=jsc,jec
         do i=isc,iec
            agm_length(i,j) = 2.0*grid_length(i,j)*rossby_radius(i,j)/(grid_length(i,j)+rossby_radius(i,j)+epsln)
         enddo
      enddo
  elseif(agm_closure_length_bczone) then 
      do j=jsc,jec
         do i=isc,iec
            agm_length(i,j) = 2.0*grid_length(i,j)*bczone_radius(i,j)/(grid_length(i,j)+bczone_radius(i,j)+epsln)
         enddo
      enddo
  endif

  if(agm_closure_baroclinic) then 
      do j=jsc,jec
         do i=isc,iec
            growth_rate(i,j) = gravrho0r_buoyr*baroclinicity(i,j)
         enddo
      enddo
  else 
      do j=jsc,jec
         do i=isc,iec
            growth_rate(i,j) = eady_rate(i,j)
         enddo
      enddo
  endif   

  !max agm_growth_rate is 2.0*agm_closure_growth_scale*coriolis_param
  !if agm_closure_growth_scale=0.5, then agm_growth_rate is <= coriolis_param
  agm_growth_rate_max(isc:iec,jsc:jec) = agm_closure_growth_scale*coriolis_param(isc:iec,jsc:jec)
  agm_array(:,:,:) = 0.0   
  k=1
  do j=jsc,jec
     do i=isc,iec
        agm_growth_rate(i,j) = 2.0*growth_rate(i,j)*agm_growth_rate_max(i,j) &
                               /(growth_rate(i,j)+agm_growth_rate_max(i,j)+epsln)  
        agm_qg(i,j)          = agm_closure_scaling*agm_growth_rate(i,j)*agm_length(i,j)**2 
        agm_array(i,j,k)     = agm_qg(i,j) + agm_micom(i,j) 
        agm_array(i,j,k)     = max(agm_closure_min, min(agm_closure_max, agm_array(i,j,k)))      
     enddo
  enddo

  do k=2,nk   
     do j=jsc,jec
        do i=isc,iec
           agm_array(i,j,k) = agm_array(i,j,1) 
        enddo
     enddo
  enddo
  call mpp_update_domains (agm_array(:,:,:), Dom%domain2d) 

  if(aredi==agm) then 
      aredi_array(:,:,:) = agm_array(:,:,:) 
  endif

  if (id_agm_growth > 0) used = send_data (id_agm_growth, agm_growth_rate(:,:), &
                                model_time, rmask=Grd%tmask(:,:,1), &
                                is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  if (id_agm_length > 0) used = send_data (id_agm_length, agm_length(:,:), &
                                model_time, rmask=Grd%tmask(:,:,1), &
                                is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  if (id_agm_qg > 0)     used = send_data (id_agm_qg, agm_qg(:,:), &
                                model_time, rmask=Grd%tmask(:,:,1), &
                                is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  if (id_growth_rate_baroclinic > 0) used = send_data (id_growth_rate_baroclinic, growth_rate(:,:), &
                                            model_time, rmask=Grd%tmask(:,:,1), &
                                            is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

  call mpp_clock_end(id_clock_compute_diffusivity)

end subroutine compute_diffusivity
! </SUBROUTINE> NAME="compute_diffusivity"



!#######################################################################
! <SUBROUTINE NAME="gm_velocity">
!
! <DESCRIPTION>
! Subroutine computes GM eddy-induced velocity field for diagnostics.
! Compute ustar and vstar at U-cell point, and wstar at T-cell bottom.   
!
! Do a two-point average rather than more democratic four-point avg
! in order to avoid having to call mpp_update domains on tensor_31 and 
! tensor_32.  The 0.5 factor is due to the two-point average.
!   
! Note that this algorithm is ad hoc.  Researchers interested in this 
! field may wish to test alternatives.  
! </DESCRIPTION>
!
subroutine gm_velocity(Thickness, tau)

  type(ocean_thickness_type), intent(in) :: Thickness
  integer, intent(in)                    :: tau 
  integer                                :: i, j, k, km1, kp1
  real, dimension(0:1)                   :: agmSx, agmSy
  real, dimension(isd:ied,jsd:jed)       :: ugm, vgm
  real                                   :: normalize=0.5

  ugm(:,:) = 0.0 ; vgm(:,:) = 0.0

  do k=1,nk         
    km1 = max(1,k-1) 
    kp1 = min(nk,k+1) 

    do j=jsc,jec
      do i=isc,iec

        agmSx(0) =  agm_array(i,j,k)*(                            &   
                    slope_function_gm(tensor_31(i,j,km1,1,0))   + &
                    slope_function_gm(tensor_31(i,j,km1,1,1)) )   
        agmSx(1) =  agm_array(i,j,k)*(                            &   
                    slope_function_gm(tensor_31(i,j,k,1,0))     + &
                    slope_function_gm(tensor_31(i,j,k,1,1)) )     
        ugm(i,j) = normalize*(agmSx(0)-agmSx(1))/Thickness%dht(i,j,k,tau)*Grd%tmask(i,j,k)*Grd%tmask(i+1,j,k) 

        agmSy(0) =  agm_array(i,j,k)*(                            &   
                    slope_function_gm(tensor_32(i,j,km1,1,0))   + &
                    slope_function_gm(tensor_32(i,j,km1,1,1)) )   
        agmSy(1) =  agm_array(i,j,k)*(                            & 
                    slope_function_gm(tensor_32(i,j,k,1,0))     + &
                    slope_function_gm(tensor_32(i,j,k,1,1)) )     
        vgm(i,j) = normalize*(agmSy(0)-agmSy(1))/Thickness%dht(i,j,k,tau)*Grd%tmask(i,j,k)*Grd%tmask(i,j+1,k) 

      enddo
    enddo

    call mpp_update_domains (ugm(:,:), Dom%domain2d, flags=EUPDATE) 
    call mpp_update_domains (vgm(:,:), Dom%domain2d, flags=NUPDATE) 
    ustar(:,:,k) = FAY(ugm(:,:))
    vstar(:,:,k) = FAX(vgm(:,:))
  enddo

  call mpp_update_domains (ustar(:,:,:), vstar(:,:,:), Dom%domain2d, gridtype=CGRID_NE)
  do k = 1, nk
     wstar(:,:,k) = wstar(:,:,k-1) + BDX_ET(ustar(:,:,k)*Thickness%dhu(:,:,k,tau)) + &
                                     BDY_NT(vstar(:,:,k)*Thickness%dhu(:,:,k,tau))
  enddo

end subroutine gm_velocity
! </SUBROUTINE> NAME="gm_velocity"



!#######################################################################
! <FUNCTION NAME="slope_function_gm">
!
! <DESCRIPTION>
! Function for defining effective slope in diagnostic GM velocity
! calculation. Used only for diagnostic purposes.  
! </DESCRIPTION>
!
function slope_function_gm (slope)

  real, intent(in) :: slope
  real             :: absslope, slope_function_gm

  absslope = abs(slope) 

  if(absslope >= smax) then
      slope_function_gm = smax*sign(1.0,slope)
  else
      slope_function_gm = slope
  endif

end function slope_function_gm
! </FUNCTION> NAME="slope_function_gm"


!#######################################################################
! <SUBROUTINE NAME="ocean_neutral_physics_rstrt">
!
! <DESCRIPTION>
! Write to restart for case when running with agm_closure
! </DESCRIPTION>
!
subroutine ocean_neutral_physics_rstrt(Time, ens_ocean)

  type(ocean_time_type), intent(in) :: Time
  logical, intent(in), optional :: ens_ocean

  integer            :: yr, mon, day, hr, min, sec
  character(len=10)  :: rdte

  if(.not. neutral_physics_on) return

  if(.not. agm_closure .and. .not. neutral_sine_taper) return

  if (.not. module_is_initialized ) then
     call mpp_error(FATAL, '==>Error from ocean_neutral_physics (ocean_neutral_physics_rstrt): module must be initialized')
  endif

  call get_date(Time%model_time,yr,mon,day,hr,min,sec)
  write(rdte,'(i4,3i2.2)') yr,mon,day,hr

  if(agm_closure ) then
    call write_data('IRESTART/'//rdte//'ocean_neutral.res','agm_array',agm_array,&
                     domain=Dom%domain2d, append_pelist_name=ens_ocean)
  endif
  if((agm_closure .and. agm_closure_length_rossby) .or. neutral_sine_taper) then
    call write_data('IRESTART/'//rdte//'ocean_neutral.res','rossby_radius',rossby_radius,&
                     domain=Dom%domain2d, append_pelist_name=ens_ocean)
  endif
  if(agm_closure .and. agm_closure_length_bczone) then
    call write_data('IRESTART/'//rdte//'ocean_neutral.res','bczone_radius',bczone_radius,&
                     domain=Dom%domain2d, append_pelist_name=ens_ocean)
  endif

  if(neutral_physics_debug) then
    write (stdout(), *) 'checksum ending agm_array', mpp_chksum(agm_array(isc:iec,jsc:jec,:)*Grd%tmask(isc:iec,jsc:jec,:))
    write (stdout(), *) 'checksum ending rossby_radius', mpp_chksum(rossby_radius(isc:iec,jsc:jec)*Grd%tmask(isc:iec,jsc:jec,1))
    write (stdout(), *) 'checksum ending bczone_radius', mpp_chksum(bczone_radius(isc:iec,jsc:jec)*Grd%tmask(isc:iec,jsc:jec,1))
  endif

end subroutine ocean_neutral_physics_rstrt
! </SUBROUTINE> NAME="ocean_neutral_physics_rstrt"


!#######################################################################
! <SUBROUTINE NAME="ocean_neutral_physics_end">
!
! <DESCRIPTION>
! Write to restart for case when running with agm_closure 
! </DESCRIPTION>
!
subroutine ocean_neutral_physics_end(ens_ocean)

  logical, intent(in), optional :: ens_ocean
  
  if(.not. neutral_physics_on) return

  if(.not. agm_closure .and. .not. neutral_sine_taper) return

  if (.not. module_is_initialized ) then 
     call mpp_error(FATAL, '==>Error from ocean_neutral_physics (ocean_neutral_physics_end): module must be initialized')
  endif 

  if(agm_closure ) then 
    call write_data('RESTART/ocean_neutral.res','agm_array',agm_array,&
                     domain=Dom%domain2d, append_pelist_name=ens_ocean)
  endif 
  if((agm_closure .and. agm_closure_length_rossby) .or. neutral_sine_taper) then 
    call write_data('RESTART/ocean_neutral.res','rossby_radius',rossby_radius,&
                     domain=Dom%domain2d, append_pelist_name=ens_ocean)
  endif 
  if(agm_closure .and. agm_closure_length_bczone) then 
    call write_data('RESTART/ocean_neutral.res','bczone_radius',bczone_radius,&
                     domain=Dom%domain2d, append_pelist_name=ens_ocean)
  endif 

  if(neutral_physics_debug) then 
    write (stdout(), *) 'checksum ending agm_array', mpp_chksum(agm_array(isc:iec,jsc:jec,:)*Grd%tmask(isc:iec,jsc:jec,:))
    write (stdout(), *) 'checksum ending rossby_radius', mpp_chksum(rossby_radius(isc:iec,jsc:jec)*Grd%tmask(isc:iec,jsc:jec,1))
    write (stdout(), *) 'checksum ending bczone_radius', mpp_chksum(bczone_radius(isc:iec,jsc:jec)*Grd%tmask(isc:iec,jsc:jec,1))
  endif 

end subroutine ocean_neutral_physics_end
! </SUBROUTINE> NAME="ocean_neutral_physics_end"

end module ocean_neutral_physics_mod
