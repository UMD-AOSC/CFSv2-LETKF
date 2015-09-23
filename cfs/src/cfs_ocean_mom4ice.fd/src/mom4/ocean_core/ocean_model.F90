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
! Modified: Xingren Wu
!           Xingren.Wu@noaa.gov
! Modified: Dave Behringer
!           David.Behringer@noaa.gov
module ocean_model_mod
!
! <CONTACT EMAIL="Stephen.Griffies@noaa.gov"> Stephen M. Griffies 
! </CONTACT>
!
!<CONTACT EMAIL="Matthew.Harrison@noaa.gov"> Matt Harrison 
!</CONTACT>
!
!<OVERVIEW>
! Time step the ocean model using either a twolevel staggered scheme
! (the default) or threelevel leap-frog scheme (the older approach).
!</OVERVIEW>
!
!<DESCRIPTION>
! Top level module for ocean model.  Contains routines for 
! initialization, termination and update of ocean model state.
!
! Design consideration: declarations of top level ocean variables
! are private to this module and hence are only available to other routines
! through argument lists.  For instance, timestep information is passed to
! the various modules on the initialization call and stored internally
! in the respective modules.  This is a crucial design consideration sinces
! it maintains modularity and hence maintainability of the code.  (mjh)
!</DESCRIPTION>
!
!<NAMELIST NAME="ocean_model_nml">
!  <DATA NAME="layout" TYPE="integer">
!  Processor domain layout for ocean model. 
!  </DATA> 
!
!  <DATA NAME="time_tendency" TYPE="character">
!
!  Possible time stepping schemes are the following. 
!
!  1. "threelevel" has the following characteristics
!
!     leap-frog for the time tendency which means the 
!     inviscid/nondissipative processes are at time tau.  
!
!     forward for lateral mixing processes (dissipation at taum1)
!
!     implicit for vertical dissipative (with aidif = 1.0)
!
!     semi-implicit for Coriolis (with acor>0) 
!
!     Because of the need to apply time filters to suppress 
!     leap-frog splitting, the threelevel time stepping scheme
!     does not conserve total tracer content in the model.  
!
!  2. "twolevel" has the following characteristics: 
!
!     staggered 2nd order forward time tendency, which means 
!     that tracer advection, lateral tracer and velocity mixing, 
!     are at time tau, and pressure gradients at taup1.  
!
!     Adams-Bashforth (either 2nd or 3rd order) for velocity advection  
!
!     implicit vertical mixing (with aidif = 1.0)
!
!     semi-implicit for Coriolis (with acor > 0) 
!
!     This scheme conserves total volume and tracer in the model.  
!
!  </DATA> 
!
!  <DATA NAME="threelevel_time_filter" TYPE="logical">
!  For running with time filtering applied to the threelevel time
!  tendency scheme.  Without time filtering, the model soon goes unstable.
!  </DATA> 
!  <DATA NAME="robert" TYPE="real">
!  Parameter setting the strength of the time filtering applied via the 
!  Robert-Asselin filter when using threelevel_time_filter=.true.
!  Typically taken as robert ~ 0.05 for global models. 
!  </DATA> 
!
!  <DATA NAME="baroclinic_split" TYPE="integer">
!  Ratio baroclinic_split = dtts/dtuv.  
!  Transients corrupted if baroclinic_split > 1.
!  </DATA> 
!  <DATA NAME="barotropic_split" TYPE="integer">
!  Ratio barotropic_split = dtuv/dtfs.  Must be 
!  large enough to resolve the barotropic gravity waves 
!  captured by the barotropic part of the model. 
!  </DATA> 
!  <DATA NAME="surface_height_split" TYPE="integer">
!  Ratio surface_height_split = dtts/dteta.  
!  Typically unity for models where baroclinic_split=1,
!  but something larger when baroclinic_split order 10.  dteta is the time 
!  step used for update of eta_t. If surface_height_split is not equal to 
!  unity, then tracer conservation properties are compromised. 
!  </DATA> 
!
!  <DATA NAME="aidif" TYPE="real">
!  aidif=1 for implicit in time solution of the vertical mixing equation.
!  aidif=0 for explicit in time solution of the vertical mixing equation.
!  semi-implicit method with 0 < aidif < 1 is not fully supported in mom4. 
!  </DATA> 
!  <DATA NAME="acor" TYPE="real">
!  acor=0.0 means explicit Coriolis force.  0.5 le acor le 1.0 means semi-implicit.
!  semi-implicit removes dtuv time step constraint associated with intertial oscillations,
!  but it leads to Coriolis force affecting energy balances.  
!  </DATA> 
!
!  <DATA NAME="velocity_change_diag_freq" TYPE="integer">
!  perform velocity_change_check every n timesteps (1 == every_tstep). 
!  </DATA> 
!  <DATA NAME="velocity_change_check" TYPE="logical">
!  For checking to see whether the abs(vel(tau)-vel(tam1))
!  is greater than velocity_delta_max.  If so, then we likely have 
!  problems with leap-frog noise.  
!  </DATA> 
!  <DATA NAME="velocity_change_max" TYPE="real" UNITS="meter/sec" >
!  For checking to see whether the abs(vel(tau)-vel(taum1)) >velocity_change_max.
!  If so, then will bring the model down if have problems in too many places.  
!  </DATA> 
!  <DATA NAME="velocity_change_max_num" TYPE="integer">
!  Maximum number of points allowed where abs(vel(tau)-vel(taum1)) >velocity_change_max.
!  </DATA> 
!
!  <DATA NAME="energy_diag_freq" TYPE="integer">
!  perform energy analysis every n timesteps (1 == every_tstep). 
!  </DATA> 
!
!  <DATA NAME="debug" TYPE="logical">
!  For overall model debugging. Set true to print cksums at 
!  each timestep for debugging purposes.
!  </DATA> 
!</NAMELIST>

use constants_mod,            only: grav, rho0, rho0r, rho_cp, epsln, hlv
use diag_manager_mod,         only: register_diag_field, send_data
use fms_mod,                  only: FATAL, NOTE, WARNING, file_exist, read_data, write_data
use fms_mod,                  only: write_version_number, open_namelist_file, close_file, check_nml_error
use fms_mod,                  only: clock_flag_default
use mpp_domains_mod,          only: mpp_global_sum, domain2d, BITWISE_EXACT_SUM
use mpp_domains_mod,          only: mpp_update_domains, BGRID_NE
use mpp_mod,                  only: mpp_error, mpp_chksum, mpp_pe, mpp_npes, stdlog, stdout
use mpp_mod,                  only: mpp_clock_id, mpp_clock_begin, mpp_clock_end
use mpp_mod,                  only: CLOCK_COMPONENT, CLOCK_SUBCOMPONENT, CLOCK_MODULE
use time_interp_external_mod, only: time_interp_external_init
use time_manager_mod,         only: JULIAN, get_date, get_time
use time_manager_mod,         only: time_type, operator( /= ), operator( < )
use time_manager_mod,         only: set_time, operator( + ), operator( .eq. )

use ocean_advection_velocity_mod, only: ocean_advection_velocity_init, ocean_advection_velocity
use ocean_bbc_mod,                only: ocean_bbc_init, get_ocean_bbc
use ocean_bih_friction_mod,       only: ocean_bih_friction_init    
use ocean_convect_mod,            only: ocean_convect_init
use ocean_coriolis_mod,           only: ocean_coriolis_init
use ocean_density_mod,            only: ocean_density_init, ocean_density_rstrt, ocean_density_end
use ocean_density_mod,            only: update_ocean_density, update_ocean_density_taup1
use ocean_diagnostics_mod,        only: ocean_diag_init, ocean_diagnostics
use ocean_domains_mod,            only: set_ocean_domain, get_local_indices
use ocean_freesurf_mod,           only: ocean_freesurf_init, ocean_freesurf_rstrt, ocean_freesurf_end
use ocean_freesurf_mod,           only: update_ocean_surface_height, update_ocean_barotropic
use ocean_freesurf_mod,           only: ocean_barotropic_forcing, surface_height_tendency, ocean_surface_smooth
use ocean_grids_mod,              only: ocean_grids_init, set_ocean_grid_size 
use ocean_grids_mod,              only: set_ocean_hgrid_arrays, set_ocean_vgrid_arrays
use ocean_horz_diffuse_mod,       only: ocean_horz_diffuse_init
use ocean_lap_friction_mod,       only: ocean_lap_friction_init
use ocean_neutral_physics_mod,    only: ocean_neutral_physics_init, neutral_physics
use ocean_obc_mod,                only: ocean_obc_init, ocean_obc_end
use ocean_obc_mod,                only: ocean_obc_update_boundary, ocean_obc_prepare
use ocean_operators_mod,          only: ocean_operators_init, REMAP_BT_TO_BU
use ocean_overflow_mod,           only: ocean_overflow_init, overflow
use ocean_polar_filter_mod,       only: ocean_polar_filter_init, polar_filter_tracers
use ocean_pressure_mod,           only: ocean_pressure_init, pressure_in_dbars
use ocean_rivermix_mod,           only: ocean_rivermix_init, rivermix
use ocean_riverspread_mod,        only: ocean_riverspread_init
use ocean_sbc_mod,                only: ocean_sbc_init, initialize_ocean_sfc, ocean_sfc_rstrt, ocean_sfc_end
use ocean_sbc_mod,                only: sum_ocean_sfc, avg_ocean_sfc, get_ocean_sbc, flux_adjust
use ocean_shortwave_pen_mod,      only: ocean_shortwave_pen_init, sw_source
use ocean_sigma_diffuse_mod,      only: ocean_sigma_diffuse_init, sigma_diffusion
use ocean_sponge_mod,             only: ocean_sponge_init, sponge_tracer_source
use ocean_thickness_mod,          only: ocean_thickness_init, update_tcell_thickness, update_ucell_thickness 
use ocean_topog_mod,              only: ocean_topog_init
use ocean_tracer_advect_mod,      only: tracer_advection_init
use ocean_tracer_mod,             only: ocean_prog_tracer_init, ocean_diag_tracer_init, ocean_tracer_rstrt
use ocean_tracer_mod,             only: update_ocean_tracer, ocean_tracer_end, compute_tmask_limit
use ocean_tracer_util_mod,        only: ocean_tracer_util_init
use ocean_tpm_mod,                only: ocean_tpm_sbc, ocean_tpm_source, ocean_tpm_bbc, ocean_tpm_tracer
use ocean_tpm_mod,                only: ocean_tpm_end, ocean_tpm_start
use ocean_types_mod,              only: ocean_time_type, ocean_time_steps_type
use ocean_types_mod,              only: ocean_grid_type, ocean_thickness_type
use ocean_types_mod,              only: ocean_domain_type, ice_ocean_boundary_type, ocean_data_type
use ocean_types_mod,              only: ocean_external_mode_type, ocean_adv_vel_type
use ocean_types_mod,              only: ocean_velocity_type, ocean_density_type 
use ocean_types_mod,              only: ocean_prog_tracer_type, ocean_diag_tracer_type
use ocean_types_mod,              only: ocean_types_init, missing_value 
use ocean_types_mod,              only: TWO_LEVEL, THREE_LEVEL
use ocean_util_mod,               only: ocean_util_init, write_timestamp
use ocean_velocity_advect_mod,    only: ocean_velocity_advect_init
use ocean_velocity_diag_mod,      only: velocity_change
use ocean_velocity_mod,           only: ocean_velocity_init, update_ocean_velocity
use ocean_velocity_mod,           only: ocean_velocity_rstrt, ocean_velocity_end
use ocean_velocity_mod,           only: ocean_explicit_accel_a, ocean_explicit_accel_b
use ocean_velocity_mod,           only: energy_analysis, ocean_implicit_coriolis, ocean_implicit_friction
use ocean_vert_mix_coeff_mod,     only: ocean_vert_mix_coeff_init, vertical_mix_coeff
use ocean_vert_mix_mod,           only: ocean_vert_mix_init
use ocean_workspace_mod,          only: ocean_workspace_init, wrk1, wrk1_v
use ocean_xlandinsert_mod,        only: ocean_xlandinsert_init, xlandinsert
use ocean_xlandmix_mod,           only: ocean_xlandmix_init, xlandmix


#ifdef ENABLE_ODA    
  use oda_driver_mod, only : init_oda
#endif

implicit none

private


#include <ocean_memory.h>

#ifdef STATIC_MEMORY

  real, dimension(isd:ied,jsd:jed,nk,2) :: accel      ! thickness weighted acceleration (m^2/s^2)
  real, dimension(isd:ied,jsd:jed,nk,2) :: diff_cbt   ! diffusion coefficient at base of tracer cells (m^2/sec)  
                                                      ! Note: separate coefficient for temp and salinity 
  real, dimension(isd:ied,jsd:jed,nk)   :: visc_cbu   ! viscosity at base of momentum cells (m^2/sec)
  real, dimension(isd:ied,jsd:jed)      :: pme        ! volume flux per horz area from from precip minus evap (m/s)
  real, dimension(isd:ied,jsd:jed,2)    :: upme       ! horizontal velocity of precip minus evap (m/s)
  real, dimension(isd:ied,jsd:jed)      :: river      ! volume flux of fresh water per horz area from river (m/s)
  real, dimension(isd:ied,jsd:jed,2)    :: uriver     ! horizontal velocity from river 
  real, dimension(isd:ied,jsd:jed)      :: patm       ! pressure at ocean top from overlying atmosphere and/or ice (Pa) 
  real, dimension(isd:ied,jsd:jed)      :: swflx      ! short wave radiation flux 
  real, dimension(isd:ied,jsd:jed,0:nk) :: sw_frac_zt ! short wave radiation flux fraction
  real, dimension(isd:ied,jsd:jed)      :: hblt_depth ! boundary layer depth from vertical mixing scheme (m)

#else

  real, pointer, dimension(:,:,:,:) :: accel      =>NULL() ! thickness weighted acceleration (m^2/s^2)
  real, pointer, dimension(:,:,:,:) :: diff_cbt   =>NULL() ! diffusion coefficient at base of tracer cells (m^2/sec) 
                                                           ! Note: separate coefficient for temp and salinity 
  real, pointer, dimension(:,:,:)   :: visc_cbu   =>NULL() ! viscosity at base of momentum cells (m^2/sec)
  real, pointer, dimension(:,:)     :: pme        =>NULL() ! volume flux per horz area from from precip minus evap (m/s)
  real, pointer, dimension(:,:,:)   :: upme       =>NULL() ! horizontal velocity of precip minus evap (m/s)
  real, pointer, dimension(:,:)     :: river      =>NULL() ! volume flux of fresh water per horz area from river (m/s) 
  real, pointer, dimension(:,:,:)   :: uriver     =>NULL() ! horizontal velocity from river (m/s)
  real, pointer, dimension(:,:)     :: patm       =>NULL() ! pressure at ocean top from sea ice and/or atmosphere (Pa) 
  real, pointer, dimension(:,:)     :: swflx      =>NULL() ! short wave radiation flux   
  real, pointer, dimension(:,:,:)   :: sw_frac_zt =>NULL() ! short wave radiation flux fraction
  real, pointer, dimension(:,:)     :: hblt_depth =>NULL() ! boundary layer depth from vertical mixing scheme (m)

#endif

  ! time step related variables 
  real :: dtts=0   ! tracer timestep (seconds)
  real :: dtuv=0   ! internal mode timestep (seconds)
  real :: dtfs=0   ! external mode timestep (seconds)
  real :: dteta=0  ! ocean volume time step (eta_t) (seconds)
  real :: p5dttsr  ! 1/(2dtts) 
  real :: p5dtuvr  ! 1/(2dtuv)
  real :: dtime_t  ! 2*dtts  for threelevel and dtts  for twolevel
  real :: dtime_u  ! 2*dtuv  for threelevel and dtvu  for twolevel
  real :: dtime_e  ! 2*dteta for threelevel and dteta for twolevel

  ! setting the number of prognostic and diagnostic tracers  
  ! both determined by reading field_table 
  integer :: num_prog_tracers=-1 ! (e.g., temp, salt, age)
  integer :: num_diag_tracers=-1 ! (e.g., frazil, pH) 

  ! for setting model time steps 
  integer :: baroclinic_split    =1  ! ratio of tracer timestep dtts to baroclinic timestep dtuv
  integer :: surface_height_split=1  ! ratio of tracer timestep dtts to the eta_t timestep dteta 
  integer :: barotropic_split    =30 ! ratio of baroclinic timestep dtuv to barotropic timestep dtfs

  ! for setting how terms in the equations are time stepped
  character(len=32) :: time_tendency='twolevel' ! selecting time tendency discretization ("threelevel" or "twolevel")
  integer :: tendency=0                         ! integer corresponding to the choice of time tendency 
  logical :: threelevel_time_filter =.true.     ! to employ time filtering for threelevel tendency. Required for stability.  
  real    :: robert=0.0                         ! dimensionless filtering coefficient used with threelevel tendency
  real    :: aidif=1.0                          ! aidif=1.0 for fully implicit vertical mixing and 0.0 for fully explicit
  real    :: acor=0.0                           ! acor=0.0 for explicit Coriolis, 0.5 <= acor <= 1.0 for semi-implicit

  ! for diagnostics sent to ascii output 
  integer :: energy_diag_freq=-1          ! perform energy analysis every n timesteps (1==every_tstep)
  integer :: velocity_change_diag_freq=-1 ! velocity_change_check every n timesteps (1==every_tstep)

  ! for diagnosing leap-frog noise in velocity field (look at the equator for amazingly large noise levels)
  logical :: velocity_change_check   =.false. ! for checking abs of the single time step change in velocity
  real    :: velocity_change_max     = 10.0   ! (m/s)
  integer :: velocity_change_max_num = 10     ! max times permit abs(vel(tau)-vel(taum1)) > velocity_change_max

  ! for diagnostics sent to diagnostic manager 
  integer, allocatable, dimension(:) :: id_robert_tracer
  integer :: id_robert_velocity(2)=-1
  integer :: id_tfilter_etat=-1
  integer :: id_tfilter_etau=-1
  integer :: id_forward_u
  integer :: id_forward_v
  logical :: used

  character(len=128) :: version = '$Id$'
  character(len=128) :: tagname = '$Name$'

  type(ocean_external_mode_type),      save                 :: Ext_mode
  type(ocean_adv_vel_type),            save                 :: Adv_vel
  type(ocean_density_type),            target, save         :: Dens
  type(ocean_domain_type),             target, save         :: Domain

  type(ocean_grid_type),               target, save         :: Grid
  type(ocean_thickness_type),          target, save         :: Thickness

  type(ocean_time_type),               target, save         :: Time
  type(ocean_time_steps_type),         target, save         :: Time_steps
  type(ocean_velocity_type),           target, save         :: Velocity

  type(ocean_prog_tracer_type), dimension(:), pointer, save :: T_prog =>NULL()
  type(ocean_diag_tracer_type), dimension(:), pointer, save :: T_diag =>NULL() 

  integer :: index_temp=-1     ! index for potential temperature
  integer :: index_salt=-1     ! index for salinity (model units psu)
  integer :: layout(2)=(/1,1/) ! domain layout for parallel processors. for npe=1, layout(2)=(/1,1/)
  logical :: debug = .false.   ! to print cksums for debugging purposes

  ! identification numbers for mpp clocks
  integer :: id_init
  integer :: id_advect
  integer :: id_vmix
  integer :: id_neutral
  integer :: id_sw
  integer :: id_sponge
  integer :: id_sigma
  integer :: id_tracer
  integer :: id_bbc
  integer :: id_sbc
  integer :: id_flux_adjust
  integer :: id_explicit_accel_a
  integer :: id_explicit_accel_b
  integer :: id_implicit_friction
  integer :: id_implicit_coriolis
  integer :: id_surface_smooth
  integer :: id_eta_tendency
  integer :: id_eta_update
  integer :: id_barotropic_forcing
  integer :: id_barotropic_update
  integer :: id_velocity
  integer :: id_polar_filter
  integer :: id_xlandinsert
  integer :: id_xlandmix
  integer :: id_overflow
  integer :: id_rivermix
  integer :: id_density
  integer :: id_density_taup1
  integer :: id_otpm_sbc
  integer :: id_otpm_source
  integer :: id_otpm_bbc
  integer :: id_otpm_tracer
  integer :: id_energy_analysis
  integer :: id_diagnostics
  integer :: id_time_filter
  integer :: id_tcell_thickness
  integer :: id_ucell_thickness
  integer :: id_update_halo_tracer
  integer :: id_update_halo_velocity
  integer :: id_oda
  integer :: id_ocean_sfc
  integer :: id_ocean_seg_end
  integer :: id_tmask_limit
  integer :: id_ocean
           
  logical :: module_is_initialized =.false.
  logical :: ens_ocean             =.false.
  logical :: have_obc              =.false.   
  
  public ocean_model_init
  public ocean_model_end
  public ocean_model_rstrt
  public update_ocean_model
  public get_ocean_domain
  public get_ocean_grid_size
  public ocean_data_type
  public ice_ocean_boundary_type
  public read_ice_ocean_boundary
  public write_ice_ocean_boundary
  public init_default_ice_ocean_boundary

  private time_filter

  namelist /ocean_model_nml/ time_tendency, threelevel_time_filter, robert,            &
                             baroclinic_split, barotropic_split, surface_height_split, &
                             aidif, acor, layout, debug, energy_diag_freq,             &
                             velocity_change_max, velocity_change_max_num,             &
                             velocity_change_check, velocity_change_diag_freq

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_model_init">
!
! <DESCRIPTION>
! Initialize the ocean model. 
! </DESCRIPTION>
!
subroutine ocean_model_init(Ocean, Time_init, Time_in, Time_step_ocean, ensemble_ocean)
    
    type(ocean_data_type), intent(inout) :: Ocean
    type(time_type), intent(in)          :: Time_init
    type(time_type), intent(in)          :: Time_in 
    type(time_type), intent(in)          :: Time_step_ocean
    logical, intent(in), optional        :: ensemble_ocean
    
    integer :: tau, secs, days
    integer :: i, n, m
    integer :: ioun, io_status, ierr

    if (module_is_initialized) then 
      call mpp_error(FATAL,'==>Error in ocean_model_mod (ocean_model_init): module already initialized')
    endif 

    module_is_initialized = .true.

    if (PRESENT(ensemble_ocean)) ens_ocean = ensemble_ocean

    if (energy_diag_freq == 0)          energy_diag_freq          = 1
    if (velocity_change_diag_freq == 0) velocity_change_diag_freq = 1

    ! set clock ids
    id_ocean               = mpp_clock_id( 'Ocean', flags=clock_flag_default, grain=CLOCK_COMPONENT )
    id_init                = mpp_clock_id('(Ocean initialization) '         ,grain=CLOCK_SUBCOMPONENT)
    id_oda                 = mpp_clock_id('(Ocean ODA)'                     ,grain=CLOCK_SUBCOMPONENT)
    id_advect              = mpp_clock_id('(Ocean advection velocity) '     ,grain=CLOCK_MODULE)
    id_density             = mpp_clock_id('(Ocean density) '                ,grain=CLOCK_MODULE)    
    id_density_taup1       = mpp_clock_id('(Ocean density at taup1) '       ,grain=CLOCK_MODULE)    
    id_vmix                = mpp_clock_id('(Ocean vertical mixing coeff) '  ,grain=CLOCK_MODULE)
    id_neutral             = mpp_clock_id('(Ocean neutral physics) '        ,grain=CLOCK_MODULE)
    id_sw                  = mpp_clock_id('(Ocean shortwave) '              ,grain=CLOCK_MODULE)
    id_sponge              = mpp_clock_id('(Ocean sponge) '                 ,grain=CLOCK_MODULE)
    id_xlandinsert         = mpp_clock_id('(Ocean xlandinsert) '            ,grain=CLOCK_MODULE)
    id_xlandmix            = mpp_clock_id('(Ocean xlandmix) '               ,grain=CLOCK_MODULE)
    id_rivermix            = mpp_clock_id('(Ocean rivermix) '               ,grain=CLOCK_MODULE)
    id_overflow            = mpp_clock_id('(Ocean overflow) '               ,grain=CLOCK_MODULE)
    id_sigma               = mpp_clock_id('(Ocean sigma diffusion) '        ,grain=CLOCK_MODULE)
    id_tracer              = mpp_clock_id('(Ocean tracer update) '          ,grain=CLOCK_MODULE)
    id_sbc                 = mpp_clock_id('(Ocean surface flux) '           ,grain=CLOCK_MODULE)
    id_bbc                 = mpp_clock_id('(Ocean bottom flux) '            ,grain=CLOCK_MODULE)
    id_flux_adjust         = mpp_clock_id('(Ocean restoring flux) '         ,grain=CLOCK_MODULE)
    id_otpm_sbc            = mpp_clock_id('(Ocean TPM sbc) '                ,grain=CLOCK_MODULE)
    id_otpm_source         = mpp_clock_id('(Ocean TPM source) '             ,grain=CLOCK_MODULE)
    id_otpm_bbc            = mpp_clock_id('(Ocean TPM bbc) '                ,grain=CLOCK_MODULE)
    id_otpm_tracer         = mpp_clock_id('(Ocean TPM tracer) '             ,grain=CLOCK_MODULE)
    id_explicit_accel_a    = mpp_clock_id('(Ocean explicit accel_a) '       ,grain=CLOCK_MODULE)
    id_explicit_accel_b    = mpp_clock_id('(Ocean explicit accel_b) '       ,grain=CLOCK_MODULE)
    id_implicit_friction   = mpp_clock_id('(Ocean implicit friction) '      ,grain=CLOCK_MODULE)
    id_implicit_coriolis   = mpp_clock_id('(Ocean implicit Coriolis) '      ,grain=CLOCK_MODULE)
    id_eta_tendency        = mpp_clock_id('(Ocean surface height tendency)' ,grain=CLOCK_MODULE)
    id_surface_smooth      = mpp_clock_id('(Ocean surface height smoothing)',grain=CLOCK_MODULE)
    id_barotropic_forcing  = mpp_clock_id('(Ocean freesurface forcing) '    ,grain=CLOCK_MODULE)
    id_eta_update          = mpp_clock_id('(Ocean surface height update)'   ,grain=CLOCK_MODULE)
    id_barotropic_update   = mpp_clock_id('(Ocean barotropic dynamics)'     ,grain=CLOCK_MODULE)
    id_velocity            = mpp_clock_id('(Ocean velocity update) '        ,grain=CLOCK_MODULE)
    id_polar_filter        = mpp_clock_id('(Ocean polar filter) '           ,grain=CLOCK_MODULE)
    id_energy_analysis     = mpp_clock_id('(Ocean energy analysis) '        ,grain=CLOCK_MODULE)
    id_diagnostics         = mpp_clock_id('(Ocean diagnostics)'             ,grain=CLOCK_MODULE)
    id_time_filter         = mpp_clock_id('(Ocean time filter for 3-level)' ,grain=CLOCK_MODULE)
    id_tcell_thickness     = mpp_clock_id('(Ocean update T-cell thickness)' ,grain=CLOCK_MODULE)
    id_ucell_thickness     = mpp_clock_id('(Ocean update U-cell thickness)' ,grain=CLOCK_MODULE)
    id_update_halo_tracer  = mpp_clock_id('(Ocean tracer halo updates)'     ,grain=CLOCK_MODULE)
    id_update_halo_velocity= mpp_clock_id('(Ocean velocity halo update)'    ,grain=CLOCK_MODULE)
    id_ocean_sfc           = mpp_clock_id('(Ocean sum ocean surface)'       ,grain=CLOCK_MODULE)
    id_ocean_seg_end       = mpp_clock_id('(Ocean average state)'           ,grain=CLOCK_MODULE)
    id_tmask_limit         = mpp_clock_id('(Ocean tracer tmask limit)'      ,grain=CLOCK_MODULE)

    call mpp_clock_begin(id_init)
    call write_version_number( version, tagname )

    write(stdout(),'(/54x,a/)') '======== STARTING MOM4 INITIALIZATION ========'

#ifdef STATIC_MEMORY
    write(stdout(),*) ' '
    write(stdout(),*)'==>NOTE: Using mom4 STATIC_MEMORY ifdef option. This option improves efficiency on SGI machines.'
    write(stdout(),*)'         However, model grid size and computer domain decomposition must be chosen at compile time.'
#else
    write(stdout(),*) ' '
    write(stdout(),*)'==>NOTE: Using dynamically allocated array option in mom4, which has proven to be slower on SGI machines.'
#endif 

    call time_interp_external_init()

    ! provide for namelist over-ride of defaults 
    ioun = open_namelist_file()
    read  (ioun, ocean_model_nml,iostat=io_status)
    write (stdout(),'(/)')
    write (stdout(), ocean_model_nml)
    write (stdlog(), ocean_model_nml)
    ierr = check_nml_error(io_status,'ocean_model_nml')
    call close_file (ioun)

    write (stdout(),*) ' ' 
    write (stdout(),*) ' ==> Note: When run as a z-coordinate ocean model, mom4p0c is Boussinesq.' 
    write (stdout(),*) '     Model kinematics conserves volume.  Mass is not conserved.'
    write (stdout(),*) '     Hence, the surface height is not affected by steric effects.'
    write (stdout(),*) '     The mom4p0 non_boussinesq option has been eliminated.'
    write (stdout(),*) ' '

    write (stdout(),'(/a,i6,a/)') ' ==>Note: Running mom4 using',mpp_npes(),' computer processors.'  

    ! initialize ocean time type information

    Time%calendar  = JULIAN

    if(time_tendency=='threelevel') then 
      tendency        = THREE_LEVEL 
      Time%taum1      = 1
      Time%tau        = 2
      Time%taup1      = 3
      Time%tau_m2     = 1
      Time%tau_m1     = 2
      Time%tau_m0     = 3
      write (stdout(),*) ' ' 
      write (stdout(),*) '==>Note: Running mom4p0c with a leap frog to discretize the time tendency.'
      write (stdout(),*) '         Unfortunately, this method does not conserve volume and tracer because'
      write (stdout(),*) '         it is necessary to use time filtering.  Use the "twolevel" scheme to conserve.'
      write (stdout(),*) ' '
      if(robert==0.0) then 
         call mpp_error(WARNING,'==>Note from ocean_model_mod: robert=0.0 is unstable with time_tendency==threelevel.')
      endif 
    elseif(time_tendency=='twolevel') then  
      tendency        = TWO_LEVEL 
      Time%taum1      = 2
      Time%tau        = 2
      Time%taup1      = 3
      Time%tau_m2     = 1
      Time%tau_m1     = 2
      Time%tau_m0     = 3
      robert          = 0.0
      write (stdout(),*) ' ' 
      write (stdout(),*) '==>Note: Running mom4p0c with staggered twotime level scheme to compute time tendencies.'
      write (stdout(),*) '         This is the default time stepping in mom4p0c since volume and tracer are conserved.'
      write (stdout(),*) ' '
    else
      call mpp_error(FATAL,'==>Error from ocean_model_mod: time_tendency must be either "twolevel" or "threelevel".')
    endif 

    Time%init       = Time_in .eq. Time_init
    Time%Time_init  = Time_init
    Time%Time_step  = Time_step_ocean
    Time%model_time = Time_in
    Time%itt        = 0

    call get_time(Time_step_ocean, secs, days)

    write(stdout(),'(/a)')' ==>Note: Time%Time_init = time stamp at very start of the mom4 experiment is given by'
    call write_timestamp(Time%Time_init)

    write(stdout(),'(/a)')' ==>Note: Time%model_time = time stamp at start of this leg of the mom4 experiment is'
    call write_timestamp(Time%model_time)

    if(Time%init) then 
      write(stdout(),'(/a)')' ==>Note: Time%init=.true. =>mom4 will start from user specified initial conditions.' 
    endif 
    if(.not. Time%init) then 
      write(stdout(),'(/a)')' ==>Note: Time%init=.false. =>mom4 will start from restart conditions from previous leg of experiment.' 
    endif 

    dtts    = secs + days*86400
    dtuv    = dtts/baroclinic_split
    dteta   = dtts/surface_height_split
    dtfs    = dtuv/barotropic_split
    p5dttsr = 0.5/(epsln+dtts)
    p5dtuvr = 0.5/(epsln+dtuv)

    if(tendency==THREE_LEVEL) then 
      dtime_t = 2.0*dtts 
      dtime_u = 2.0*dtuv 
      dtime_e = 2.0*dteta
    elseif(tendency==TWO_LEVEL) then  
      dtime_t = dtts 
      dtime_u = dtuv
      dtime_e = dteta
    endif 

    Time_steps%time_tendency = time_tendency 
    Time_steps%tendency      = tendency 
    Time_steps%aidif         = aidif 
    Time_steps%acor          = acor
    Time_steps%dtts          = dtts
    Time_steps%dtuv          = dtuv
    Time_steps%dtfs          = dtfs
    Time_steps%dteta         = dteta
    Time_steps%dtime_t       = dtime_t
    Time_steps%dtime_u       = dtime_u
    Time_steps%dtime_e       = dtime_e

    if(tendency==TWO_LEVEL .and. acor==0.0) then  
       call mpp_error(WARNING,'==>ocean_model_mod: acor=0.0 and twolevel tendency is weakly unstable. Recommend 0.5 <= acor <=1.0')
    endif 

    write(stdout(),'(/a)')' ==> Note: time steps (seconds) used for mom4' 
    write(stdout(),'(a,f10.2)')'  dtts  (tracer)         = ',dtts
    write(stdout(),'(a,f10.2)')'  dtuv  (baroclinic)     = ',dtuv
    write(stdout(),'(a,f10.2)')'  dteta (surface height) = ',dteta
    write(stdout(),'(a,f10.2)')'  dtfs  (barotropic)     = ',dtfs

    if(baroclinic_split < 1) then 
       call mpp_error(FATAL,'==>Error from ocean_model_mod(ocean_model_init): baroclinic_split must be an integer >= 1')
    endif 
    if(baroclinic_split > 1) then 
       call mpp_error(NOTE,&
       '==>Note from ocean_model_mod(ocean_model_init): baroclinic_split > 1 will corrupt transients')
    endif 
    if(surface_height_split < 1) then 
       call mpp_error(FATAL,'==>Error from ocean_model_mod(ocean_model_init): surface_height_split must be an integer >= 1')
    endif 
    if(surface_height_split > 1) then 
       call mpp_error(NOTE, &
       '==>Note from ocean_model_mod(ocean_model_init): surface_height_split > 1 corrupts tracer conservation and transients')
    endif  

    if(barotropic_split < 1) then 
       call mpp_error(FATAL, &
       '==>Error from ocean_model_mod(ocean_model_init): barotropic_split must be an integer >= 1')
    endif 

    if(.not. threelevel_time_filter .and. tendency==THREE_LEVEL) then 
       call mpp_error(WARNING, &
       '==>Note from ocean_model_mod: threelevel_time_filter=.false. results in unstable leap-frog algorithm.')
    endif 

    if (nint(dtuv) /= nint(dtfs)) then
        write (stdout(),'(/1x,a/)') ' ==> Note: The velocity equations will be split into baroclinic and barotropic pieces.'
    elseif (nint(dtuv) /= 0) then
        write (stdout(),'(/1x,a/)') ' ==> Warning: The velocity equations will not be split into baroclinic and barotropic pieces.'
    endif

    if(0.0 < Time_steps%aidif .and. Time_steps%aidif < 1.0) then 
      call mpp_error(FATAL,'==>Error: ocean_model_init: set aidif==0.0 OR aidif==1.0. Intermediate values not supported.')  
    endif

   
    ! initialize grid and domain information
    call ocean_grids_init(debug)
    call set_ocean_grid_size(Grid, 'INPUT/grid_spec') 
    call set_ocean_domain(Domain, Grid, layout=layout)
    print *,'in ocean_model_init, after set_ocean_domain'
    call ocean_workspace_init(Domain, Grid)
    print *,'in ocean_model_init, after ocean_workspace_init'
    call set_ocean_hgrid_arrays(Time, Domain, Grid)
    print *,'in ocean_model_init, after set_ocean_hgrid_arrays'
    call ocean_topog_init(Domain, Grid, 'INPUT/grid_spec')
    print *,'in ocean_model_init, ocean_topog_init'
    call ocean_obc_init(have_obc, Time, Time_steps, Domain, Grid, debug)
    print *,'in ocean_model_init, after ocean_obc_init'
    call set_ocean_vgrid_arrays(Time, Domain, Grid, have_obc)
    print *,'in ocean_model_init, after set_ocean_vgrid_arrays'
    call ocean_util_init(Domain)
    print *,'in ocean_model_init, after ocean_util-init'


#ifndef STATIC_MEMORY
    call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
    nk = Grid%nk

    allocate(accel(isd:ied,jsd:jed,nk,2))
    allocate(visc_cbu(isd:ied,jsd:jed,nk))
    allocate(diff_cbt(isd:ied,jsd:jed,nk,2))
    allocate(pme(isd:ied,jsd:jed))
    allocate(upme(isd:ied,jsd:jed,2))
    allocate(river(isd:ied,jsd:jed))
    allocate(uriver(isd:ied,jsd:jed,2))
    allocate(patm(isd:ied,jsd:jed))
    allocate(swflx(isd:ied,jsd:jed))    
    allocate(sw_frac_zt(isd:ied,jsd:jed,0:nk))    
    allocate(hblt_depth(isd:ied,jsd:jed))    
#endif
    accel      = 0.0
    visc_cbu   = 0.0
    diff_cbt   = 0.0
    pme        = 0.0
    upme       = 0.0
    river      = 0.0
    uriver     = 0.0
    patm       = 0.0
    swflx      = 0.0    
    sw_frac_zt = 0.0    
    hblt_depth = 0.0

    ! initialize remaining modules     
    call ocean_types_init()
    call ocean_tracer_util_init(Grid, Domain)
    call ocean_coriolis_init(Grid, Domain, Time, Time_steps)
    call ocean_velocity_init(Grid, Domain, Time, Time_steps, Velocity, have_obc, debug, ens_ocean)
    call ocean_freesurf_init(Grid, Domain, Time, Time_steps, Ext_mode, have_obc, debug, ens_ocean)
    call ocean_thickness_init(Time, Domain, Grid, Ext_mode, Thickness, debug)
    call ocean_operators_init(Grid, Domain, Thickness)

    T_prog => ocean_prog_tracer_init(Grid, Thickness, Domain, Time, Time_steps, index_temp, index_salt, &
                                     num_prog_tracers, have_obc, debug, ens_ocean)
    T_diag => ocean_diag_tracer_init(Time, Thickness, num_diag_tracers, debug, ens_ocean)

    call ocean_advection_velocity_init(Grid, Domain, Time, Time_steps, Adv_vel, have_obc, debug)
    call ocean_density_init(Grid, Domain, Time, Time_steps, Thickness, T_prog(:), Dens, debug, ens_ocean)
    call ocean_pressure_init(Grid, Domain, Time, have_obc)

    print *,'in ocean_model_init, after ocean_pressure_init'

    call ocean_horz_diffuse_init(Grid, Domain, Time, T_prog(:), dtime_t, have_obc)
    call ocean_sigma_diffuse_init(Grid, Domain, Time, Thickness, T_prog(:), dtime_t)
    call ocean_neutral_physics_init(Grid, Domain, Time, Time_steps, Thickness, T_prog(:), ens_ocean)
    call ocean_lap_friction_init(Grid, Domain, Time, dtime_u)
    call ocean_bih_friction_init(Grid, Domain, Time, dtime_u)    
    call ocean_vert_mix_coeff_init(Grid, Domain, Time, Time_steps, T_prog(:), T_diag(:))
    call ocean_vert_mix_init(Grid, Domain, Time, Time_steps, T_prog(:))
    call tracer_advection_init(Grid, Domain, Time, T_prog(:), have_obc)
    call ocean_velocity_advect_init(Grid, Domain, Time, have_obc, debug)
    call ocean_convect_init(Grid, Domain, Time, T_prog(:), dtime_t, index_temp, index_salt)
    call ocean_sbc_init(Grid, Domain, Time, T_prog(:), T_diag(:), Velocity, Ocean, time_tendency)
    call ocean_bbc_init(Grid, Domain, T_prog(:))
    call ocean_shortwave_pen_init(Grid, Domain, Time)
    call ocean_sponge_init(Grid, Domain, T_prog(:), dtime_t)
    call ocean_xlandmix_init(Grid, Domain, Time, Thickness, T_prog(:), dtime_t)    
    call ocean_xlandinsert_init(Grid, Domain, Time, T_prog(:), dtts)    
    call ocean_riverspread_init(Grid, Domain)
    call ocean_rivermix_init(Grid, Domain, Time, Time_steps, T_prog(:), debug)    
    call ocean_overflow_init(Grid, Domain, Time, T_prog(:), debug)    
    call ocean_polar_filter_init(Grid, Domain, Time, T_prog(:), dtime_t, index_temp)

    call initialize_ocean_sfc(Time, Thickness, T_prog(:), Velocity, patm, Ocean)
    call ocean_tpm_start(T_prog(:), T_diag(:))
    call ocean_diag_init(Grid, Domain, Time, Time_steps, Adv_vel, T_prog(:), T_diag(:), &
                         Velocity, Ext_mode, Dens, have_obc)

    print *,'in ocean_model_init, after ocean_tpm_start'
#ifdef ENABLE_ODA    
    call init_oda(Domain, Grid, T_prog(:), Velocity, Ext_mode)
#endif

   ! register diagnostics 

   id_robert_velocity(1) = register_diag_field ('ocean_model', 'robert_u', &
          Grid%vel_axes_uv(1:3), Time%model_time, &
          'u-accel from robert time filter', 'm/s^2',&
           missing_value=missing_value, range=(/-1e6,1e6/))
   id_robert_velocity(2) = register_diag_field ('ocean_model', 'robert_v', &
          Grid%vel_axes_uv(1:3), Time%model_time, &
          'v-accel from robert time filter', 'm/s^2',&
           missing_value=missing_value, range=(/-1e6,1e6/))

   id_tfilter_etat = register_diag_field ('ocean_model', 'tfilter_etat', &
          Grid%tracer_axes(1:2), Time%model_time, &
          'eta_t change from time filter', 'm',&
           missing_value=missing_value, range=(/-1e6,1e6/))

   id_tfilter_etau = register_diag_field ('ocean_model', 'tfilter_etau', &
          Grid%vel_axes_uv(1:2), Time%model_time, &
          'eta_u change from time filter', 'm',&
           missing_value=missing_value, range=(/-1e6,1e6/))

   id_forward_u = register_diag_field ('ocean_model', 'forward_u', &
          Grid%vel_axes_uv(1:3), Time%model_time, &
          'u(tau)-u(taum1)', 'm/s',&
           missing_value=missing_value, range=(/-1e6,1e6/))
   id_forward_v = register_diag_field ('ocean_model', 'forward_v', &
          Grid%vel_axes_uv(1:3), Time%model_time, &
          'v(tau)-v(taum1)', 'm/s',&
           missing_value=missing_value, range=(/-1e6,1e6/))


    allocate( id_robert_tracer(num_prog_tracers) )
    id_robert_tracer(:) = -1
    do n=1,num_prog_tracers
       if (n == index_temp) then
           id_robert_tracer(n) = register_diag_field ('ocean_model', 'robert_'//trim(T_prog(n)%name), &
                Grid%tracer_axes(1:3), Time%model_time, &
                'thk wghtd heating from robert filter', 'Watts/m^2',&
                missing_value=missing_value, range=(/-1e6,1e6/))
       else 
           id_robert_tracer(n) = register_diag_field ('ocean_model', 'robert_'//trim(T_prog(n)%name), &
                Grid%tracer_axes(1:3), Time%model_time, &
                'thk wghtd tendency from robert filter for '//trim(T_prog(n)%longname), 'm*kg/sec',&
                missing_value=missing_value, range=(/-1e6,1e6/))
       endif
    enddo

    call mpp_clock_end(id_init) 

    write(stdout(),'(/52x,a/)') '======== COMPLETED MOM4 INITIALIZATION ========'

  end subroutine ocean_model_init
! </SUBROUTINE> NAME="ocean_model_init"


 !#######################################################################
! <SUBROUTINE NAME="update_ocean_model">
!
! <DESCRIPTION>
! Update in time the ocean model fields. 
! </DESCRIPTION>
!
 subroutine update_ocean_model(Ice_ocean_boundary, Ocean, ocean_seg_start, ocean_seg_end, num_ocean_calls)

    type(ice_ocean_boundary_type), intent(inout) :: Ice_ocean_boundary
    type(ocean_data_type), intent(inout)         :: Ocean
    logical, intent(in)                          :: ocean_seg_start
    logical, intent(in)                          :: ocean_seg_end
    integer, intent(in)                          :: num_ocean_calls

    integer :: taum1, tau, taup1
    integer :: i, j, k, n

    call mpp_clock_begin(id_ocean)

    ! increment ocean time and time labels 
    Time%model_time = Time%model_time + Time%Time_step
    Time%itt        = Time%itt+1

    Time%taum1      = mod(Time%itt+0,3)+1
    Time%tau        = mod(Time%itt+1,3)+1
    Time%taup1      = mod(Time%itt+2,3)+1

    Time%tau_m2     = mod(Time%itt+0,3)+1
    Time%tau_m1     = mod(Time%itt+1,3)+1
    Time%tau_m0     = mod(Time%itt+2,3)+1

    if(tendency==TWO_LEVEL) then  
      Time%taum1 = Time%tau
    endif 

    taum1 = Time%taum1
    tau   = Time%tau
    taup1 = Time%taup1

    ! initialize eta source to zero
    do j=jsd,jed
       do i=isd,ied
          Ext_mode%eta_source(i,j) = 0.0 
       enddo
    enddo
    ! initialize thickness weighted tracer tendency to zero
    do n=1,num_prog_tracers
       do k=1,nk
          do j=jsd,jed
             do i=isd,ied
                T_prog(n)%th_tendency(i,j,k) = 0.0
             enddo
          enddo
       enddo
    enddo

    ! initialize vertical diffusivity
    diff_cbt(:,:,:,:) = 0.0

    ! calculate tracer tmask_limit based on tracer values at time tau
    call mpp_clock_begin(id_tmask_limit)
    call compute_tmask_limit(Time, T_prog(1:num_prog_tracers))
    call mpp_clock_end(id_tmask_limit)

    ! calculate surface boundary fluxes using boundary field from coupler
    call mpp_clock_begin(id_sbc)
    call get_ocean_sbc(Time, Ice_ocean_boundary, Ext_mode, T_prog(1:num_prog_tracers), Velocity, &
                       pme, upme, river, uriver, swflx, patm)
    call mpp_clock_end(id_sbc)

    ! set ocean surface boundary conditions for the tracer packages
    call mpp_clock_begin(id_otpm_sbc)
    call ocean_tpm_sbc(robert)
    call mpp_clock_end(id_otpm_sbc)

    ! compute "flux adjustments" (e.g., surface restoring)
    call mpp_clock_begin(id_flux_adjust)
    call flux_adjust(Time, T_prog(1:num_prog_tracers), Velocity, pme)
    call mpp_clock_end(id_flux_adjust)

    ! calculate bottom momentum and bottom tracer fluxes
    call mpp_clock_begin(id_bbc)
    call get_ocean_bbc(Time, Velocity, T_prog(1:num_prog_tracers))
    call mpp_clock_end(id_bbc)
    
    ! get prescribed OBC data from files
    if (have_obc) call ocean_obc_prepare(Time, Ext_mode, T_prog(1:num_prog_tracers))

    ! compute quantities related to ocean density at time level tau
    call mpp_clock_begin(id_density)
    call update_ocean_density(Time, Thickness, T_prog(index_temp), T_prog(index_salt), Ext_mode, patm, Dens)
    call mpp_clock_end(id_density)

    ! compute advective velocity components on faces of T-cells and U-cells 
    call mpp_clock_begin(id_advect)
    call ocean_advection_velocity(Ext_mode, Velocity, Time, Thickness, Adv_vel)
    call mpp_clock_end(id_advect)

    ! add thickness weighted tendency to tracer source from shortwave heating
    call mpp_clock_begin(id_sw)
    call sw_source(Time, swflx, T_prog(index_temp), sw_frac_zt)   
    call mpp_clock_end(id_sw)

    ! compute vertical mixing coefficients
    ! if using kpp, then add thickness weighted nonlocal tendency to tracer source
    call mpp_clock_begin(id_vmix)    
    call vertical_mix_coeff(aidif, Time, Thickness, Velocity, T_prog(1:num_prog_tracers), &
                            T_diag(1:num_diag_tracers), Dens, swflx, sw_frac_zt, pme, &
                            river, visc_cbu, diff_cbt, hblt_depth)
    call mpp_clock_end(id_vmix)

    ! add thickness weighted tendency from neutral physics to tracer source.
    ! also add K33 contribution to vertical diffusivity for implicit update.
    ! hblt needed from kpp to help define "neutral physics boundary layer". 
    call mpp_clock_begin(id_neutral)
    call neutral_physics(Time, Thickness, Dens%rho_taum1, Dens%pressure_at_depth, T_prog(1:num_prog_tracers), hblt_depth)
    call mpp_clock_end(id_neutral)

    ! set the ocean source terms for the tracer packages
    call mpp_clock_begin(id_otpm_source)
    call ocean_tpm_source(Thickness)
    call mpp_clock_end(id_otpm_source)

    ! set the ocean surface boundary conditions for the tracer packages
    call mpp_clock_begin(id_otpm_bbc)
    call ocean_tpm_bbc
    call mpp_clock_end(id_otpm_bbc)
  
    ! add thickness weighted tendency from sponges to tracer source
    call mpp_clock_begin(id_sponge)
    call sponge_tracer_source(Time, Thickness, T_prog(1:num_prog_tracers)) 
    call mpp_clock_end(id_sponge)

    ! add thickness weighted tendency from cross land mixing to tracer source.
    ! add tendency to surface height source  
    call mpp_clock_begin(id_xlandmix)
    call xlandmix (Time, Thickness, T_prog(1:num_prog_tracers), Ext_mode)
    call mpp_clock_end(id_xlandmix)

    ! add thickness weighted tendency from cross land insertion to tracer source.
    ! add tendency to surface height source  
    call mpp_clock_begin(id_xlandinsert)
    call xlandinsert (Time, Thickness, T_prog(1:num_prog_tracers), Ext_mode)
    call mpp_clock_end(id_xlandinsert)

    ! add thickness weighted tendency from river discharge to tracer source 
    ! and/or enhance diff_cbt next to river mouths 
    call mpp_clock_begin(id_rivermix)
    call rivermix (Time, Thickness, T_prog(1:num_prog_tracers), Ext_mode, river, diff_cbt, index_temp, index_salt)
    call mpp_clock_end(id_rivermix)
    
    ! add thickness weighted tendency to tracer source due to discharge of dense 
    ! shelf water into abyss, with associated upstream tracer advection
    call mpp_clock_begin(id_overflow)
    call overflow (Time, Thickness, T_prog(1:num_prog_tracers), Dens, index_temp, index_salt, dtime_t)
    call mpp_clock_end(id_overflow)
    
    ! add thickness weighted tendency from sigma diffusion to tracer source
    call mpp_clock_begin(id_sigma)
    call sigma_diffusion(Time, Thickness, T_prog(1:num_prog_tracers), Adv_vel) 
    call mpp_clock_end(id_sigma)

    ! compute surface smoother source for eta and tracer 
    call mpp_clock_begin(id_surface_smooth)
    call ocean_surface_smooth(Time, Thickness, T_prog(1:num_prog_tracers), Ext_mode) 
    call mpp_clock_end(id_surface_smooth)

    ! compute surface height tendency (including contribution from eta_source)
    call mpp_clock_begin(id_eta_tendency)
    call surface_height_tendency(Time, Adv_vel, pme, river, Ext_mode)
    call mpp_clock_end(id_eta_tendency)

    ! update taup1 value of the ocean free surface height eta_t using "big time step"
    call mpp_clock_begin(id_eta_update)
    call update_ocean_surface_height(Time, Ext_mode, Dens, patm, pme, river)
    call mpp_clock_end(id_eta_update)

    ! update taup1 value of the tracer grid cell thickness
    call mpp_clock_begin(id_tcell_thickness)
    call update_tcell_thickness(Time, Ext_mode, Grid, Thickness) 
    call mpp_clock_end(id_tcell_thickness)

    ! update taup1 value of tracer concentrations (compatible with eta_t update)   
    call mpp_clock_begin(id_tracer)
    call update_ocean_tracer(Time, Dens, Adv_vel, Thickness, Ext_mode, pme, diff_cbt, &
                             T_prog(1:num_prog_tracers), T_diag(1:num_diag_tracers))
    call mpp_clock_end(id_tracer)

    ! apply polar filter to taup1 tracer concentrations
    call mpp_clock_begin(id_polar_filter) 
    call polar_filter_tracers(Time, Thickness, T_prog(1:num_prog_tracers), index_temp)   
    call mpp_clock_end(id_polar_filter) 

    ! perform extra calculations for the ocean tracer packages
    call mpp_clock_begin(id_otpm_tracer)
    call ocean_tpm_tracer
    call mpp_clock_end(id_otpm_tracer)

    ! fill halos for tracers(taup1) as this is needed prior to rho_taup1 computation 
    call mpp_clock_begin(id_update_halo_tracer)
    do n=1,num_prog_tracers
       call mpp_update_domains(T_prog(n)%field(:,:,:,taup1), Domain%domain2d, complete=T_prog(n)%complete)
    enddo
    do n=1,num_prog_tracers 
       if(have_obc)call ocean_obc_update_boundary(T_prog(n)%field(:,:,:,taup1), 'T')
    enddo
    call mpp_clock_end(id_update_halo_tracer)

    ! taup1 density to get pressure gradient for staggered time stepping scheme
    if(tendency == TWO_LEVEL) then 
      call mpp_clock_begin(id_density_taup1) 
      call update_ocean_density_taup1(Time, Thickness, T_prog(index_temp), T_prog(index_salt), Ext_mode, Dens, patm)   
      call mpp_clock_end(id_density_taup1) 
    endif 

    ! time explicit contributions to thickness weighted acceleration  
    ! compute just those pieces needed to force barotropic dynamics 
    call mpp_clock_begin(id_explicit_accel_a)
    if(tendency == TWO_LEVEL) then 
      call ocean_explicit_accel_a(Velocity, Time, Adv_vel, Thickness, Dens%rho_taup1, pme, river, upme, uriver, accel)
    elseif(tendency == THREE_LEVEL) then 
      call ocean_explicit_accel_a(Velocity, Time, Adv_vel, Thickness, Dens%rho, pme, river, upme, uriver, accel)
    endif 
    call mpp_clock_end(id_explicit_accel_a)

    ! vertical integral of forcing used for barotropic dynamics 
    call mpp_clock_begin(id_barotropic_forcing)
    call ocean_barotropic_forcing(Time, Velocity, accel, Ext_mode) 
    call mpp_clock_end(id_barotropic_forcing)

    ! update (ud,vd) and eta_t_bar using barotropic timesteps 
    call mpp_clock_begin(id_barotropic_update)
    call update_ocean_barotropic (Time, Dens, Ext_mode, patm, pme, river)
    call mpp_clock_end(id_barotropic_update)

    ! update taup1 value of the velocity grid cell thickness
    call mpp_clock_begin(id_ucell_thickness)
    call update_ucell_thickness(Time, Ext_mode, Grid, Thickness) 
    call mpp_clock_end(id_ucell_thickness)

    ! add remaining time explicit contributions to thickness weighted acceleration
    call mpp_clock_begin(id_explicit_accel_b)
    call ocean_explicit_accel_b(visc_cbu, Time, Velocity, Thickness, accel)
    call mpp_clock_end(id_explicit_accel_b)

    ! thickness weighted acceleration due to implicit vertical friction
    call mpp_clock_begin(id_implicit_friction)
    call ocean_implicit_friction(Time, Thickness, visc_cbu, Velocity, accel) 
    call mpp_clock_end(id_implicit_friction)

    ! thickness weighted acceleration due to implicit coriolis force 
    call mpp_clock_begin(id_implicit_coriolis)
    call ocean_implicit_coriolis(Time, Velocity, acor, accel) 
    call mpp_clock_end(id_implicit_coriolis)

    ! update taup1 value of ocean velocity
    call mpp_clock_begin(id_velocity)
    call update_ocean_velocity(Time, Thickness, accel, barotropic_split, Ext_mode, Velocity) 
    call mpp_clock_end(id_velocity)

    ! perform energy analysis 
    call mpp_clock_begin(id_energy_analysis)
    if (energy_diag_freq > 0 .and. mod(Time%itt, energy_diag_freq) == 0) then
        if(tendency == TWO_LEVEL) then 
           call energy_analysis (Time, Thickness, Dens%rho_taup1, pme, river, upme, uriver, &
                                 visc_cbu, Dens, Ext_mode, Adv_vel, Velocity, acor)
        elseif(tendency == THREE_LEVEL) then 
           call energy_analysis (Time, Thickness, Dens%rho, pme, river, upme, uriver, &
                                 visc_cbu, Dens, Ext_mode, Adv_vel, Velocity, acor)
        endif 
    endif
    call mpp_clock_end(id_energy_analysis)

    ! perform other numerical diagnostics 
    call mpp_clock_begin(id_diagnostics)
    call ocean_diagnostics(Time, Thickness, T_prog(1:num_prog_tracers), T_diag(1:num_diag_tracers), &
                           Velocity, Adv_vel, Ext_mode, Dens, pme, river, visc_cbu)
    call mpp_clock_end(id_diagnostics)


    ! fill velocity halos
    call mpp_clock_begin(id_update_halo_velocity)
    call mpp_update_domains(Velocity%u(:,:,:,1,taup1), Velocity%u(:,:,:,2,taup1), Domain%domain2d,gridtype=BGRID_NE)
    if(have_obc) then
       call ocean_obc_update_boundary(Velocity%u(:,:,:,1,taup1), 'M','s')
       call ocean_obc_update_boundary(Velocity%u(:,:,:,2,taup1), 'M','i')
       call ocean_obc_update_boundary(Velocity%u(:,:,:,1,taup1), 'Z','i')
       call ocean_obc_update_boundary(Velocity%u(:,:,:,2,taup1), 'Z','s')
    endif
    call mpp_clock_end(id_update_halo_velocity)
    
    ! apply time filter when using leap-frog time tendency 
    call mpp_clock_begin(id_time_filter)
    if(threelevel_time_filter .and. tendency==THREE_LEVEL) then 
      call time_filter(Time, Grid, Thickness, T_prog(1:num_prog_tracers), Velocity, Ext_mode)
    endif 
    call mpp_clock_end(id_time_filter)

    ! save previous density
    Dens%rho_taum1=Dens%rho      


    ! modifications to prognostic variables using ocean data assimilation 
#ifdef ENABLE_ODA
    call mpp_clock_begin(id_oda)
    call oda(T_prog(1:num_prog_tracers), Velocity, Ext_mode)
    call mpp_clock_end(id_oda)
#endif

    ! sum ocean sfc state over coupling interval
    call mpp_clock_begin(id_ocean_sfc)
    call sum_ocean_sfc(Time, Thickness, T_prog(1:num_prog_tracers), T_diag(1:num_diag_tracers), Velocity, Ocean)
    call mpp_clock_end(id_ocean_sfc)

    ! at end of coupling interval, pass averaged ocean state to other component models 
    if (ocean_seg_end) then 
       call mpp_clock_begin(id_ocean_seg_end)
       call avg_ocean_sfc(Time, Thickness, T_prog(1:num_prog_tracers), Velocity, patm, Ocean)
       call mpp_clock_end(id_ocean_seg_end)
    endif 


    call mpp_clock_end(id_ocean)

    return

  end subroutine update_ocean_model
! </SUBROUTINE> NAME="update_ocean_model"



!#######################################################################
! <SUBROUTINE NAME="get_ocean_grid_size">
!
! <DESCRIPTION>
! Obtain the ocean grid size. 
! </DESCRIPTION>
!
  subroutine get_ocean_grid_size(num_lon, num_lat, num_z)

    integer, intent(out) :: num_lon, num_lat
    integer, intent(out), optional :: num_z

    num_lon = Grid%ni
    num_lat = Grid%nj
    if (PRESENT(num_z)) num_z   = Grid%nk

    return
    
  end subroutine get_ocean_grid_size
! </SUBROUTINE> NAME="get_ocean_grid_size"



!#######################################################################
! <SUBROUTINE NAME="get_ocean_domain">
!
! <DESCRIPTION>
! Obtain the ocean domain size. 
! </DESCRIPTION>
!
  subroutine get_ocean_domain(Ocean_domain)

    type(domain2d), intent(out) :: Ocean_domain

    Ocean_domain = Domain%domain2d

    return

  end subroutine get_ocean_domain
! </SUBROUTINE> NAME="get_ocean_domain"

!#######################################################################
! <SUBROUTINE NAME="ocean_model_rstrt">
!
! <DESCRIPTION>
! Make an ocean restart.
! Do not close down the ocean model
! </DESCRIPTION>
!
  subroutine ocean_model_rstrt(Ocean)

    type(ocean_data_type), intent(in) :: Ocean

    call ocean_tracer_rstrt(Time, T_prog(:), T_diag(:), ens_ocean)
    call ocean_velocity_rstrt(Time, Velocity, ens_ocean)
    call ocean_freesurf_rstrt(Time, Ext_mode, ens_ocean)
    call ocean_density_rstrt(Time, Dens, ens_ocean)
    call ocean_sfc_rstrt(Time, Ocean, ens_ocean)

    return

  end subroutine ocean_model_rstrt
! </SUBROUTINE> NAME="ocean_model_rstrt"

!#######################################################################
! <SUBROUTINE NAME="ocean_model_end">
!
! <DESCRIPTION>
! Close down the ocean model 
! </DESCRIPTION>
!
  subroutine ocean_model_end(Ocean)

    type(ocean_data_type), intent(in) :: Ocean
    integer                           :: ocean_t_points 
    integer                           :: ocean_u_points 
    integer                           :: total_points

    call ocean_tracer_end(Time, T_prog(:), T_diag(:), ens_ocean)
    call ocean_tpm_end(Thickness)
    call ocean_velocity_end(Time, Velocity, ens_ocean)
    call ocean_freesurf_end(Time, Ext_mode, ens_ocean)
    call ocean_density_end(Time, Dens, ens_ocean)
    if(have_obc) call ocean_obc_end(have_obc)
    call ocean_sfc_end(Ocean, ens_ocean)

    total_points   = Grid%ni*Grid%nj*Grid%nk
    ocean_t_points = nint(mpp_global_sum(Domain%domain2d,Grid%tmask(:,:,:), BITWISE_EXACT_SUM)) 
    ocean_u_points = nint(mpp_global_sum(Domain%domain2d,Grid%umask(:,:,:), BITWISE_EXACT_SUM)) 

    write (stdout(),'(//,1x,a)') '==================Summary of completed mom4 integration===================='
    write (stdout(),'(1x,a,1x,i12)')  ' number of time steps                             = ',Time%itt
    write (stdout(),'(1x,a,1x,i12)')  ' number of prog-tracers                           = ',num_prog_tracers
    write (stdout(),'(1x,a,1x,i12)')  ' number of diag-tracers                           = ',num_diag_tracers
    write (stdout(),'(1x,a,1x,i12)')  ' number of i-points(ni)                           = ',Grid%ni
    write (stdout(),'(1x,a,1x,i12)')  ' number of j-points(nj)                           = ',Grid%nj
    write (stdout(),'(1x,a,1x,i12)')  ' number of k-points(nk)                           = ',Grid%nk
    write (stdout(),'(1x,a,1x,i12)')  ' number of computed ocean tracer points(ni*nj*nk) = ',total_points
    write (stdout(),'(1x,a,1x,i12)')  ' number of wet ocean tracer points                = ',ocean_t_points
    write (stdout(),'(1x,a,1x,i12)')  ' number of wet ocean velocity points              = ',ocean_u_points
    write (stdout(),'(1x,a/)')   '==========================================================='
    
    return

  end subroutine ocean_model_end
! </SUBROUTINE> NAME="ocean_model_end"



!#######################################################################
! <SUBROUTINE NAME="write_ice_ocean_boundary">
!
! <DESCRIPTION>
! Write the surface boundary conditions for use in coupled modeling.
! </DESCRIPTION>
!
  subroutine write_ice_ocean_boundary(file_name,iob,Ocean)

    character(LEN=*),             intent(IN) :: file_name
    type(ice_ocean_boundary_type),intent(IN) :: iob
    type(ocean_data_type),        intent(IN) :: Ocean

    call write_data(file_name,'u_flux',   iob%u_flux,   Ocean%Domain )
    call write_data(file_name,'v_flux',   iob%v_flux,   Ocean%Domain )
    call write_data(file_name,'t_flux',   iob%t_flux,   Ocean%Domain )
    call write_data(file_name,'q_flux',   iob%q_flux,   Ocean%Domain )
    call write_data(file_name,'salt_flux',iob%salt_flux,Ocean%Domain )
    call write_data(file_name,'lw_flux',  iob%lw_flux,  Ocean%Domain )
    call write_data(file_name,'sw_flux',  iob%sw_flux,  Ocean%Domain )
    call write_data(file_name,'lprec',    iob%lprec,    Ocean%Domain )
    call write_data(file_name,'fprec',    iob%fprec,    Ocean%Domain )
    call write_data(file_name,'runoff',   iob%runoff,   Ocean%Domain )
    call write_data(file_name,'calving',  iob%calving,  Ocean%Domain )
    call write_data(file_name,'p',        iob%p,        Ocean%Domain )

! iob%xtype set by flux_exchange_init

  end subroutine write_ice_ocean_boundary
! </SUBROUTINE> NAME="write_ice_ocean_boundary"


!#######################################################################
! <SUBROUTINE NAME="read_ice_ocean_boundary">
!
! <DESCRIPTION>
! Read the surface boundary conditions for use in coupled modeling.
! </DESCRIPTION>
!
  subroutine read_ice_ocean_boundary(file_name,iob,Ocean)

    character(LEN=*),             intent(IN)    :: file_name  
    type(ice_ocean_boundary_type),intent(INOUT) :: iob
    type(ocean_data_type),        intent(IN)    :: Ocean

    if (file_exist(file_name)) then
        call read_data(file_name,'u_flux',   iob%u_flux,   Ocean%Domain)
        call read_data(file_name,'v_flux',   iob%v_flux,   Ocean%Domain)
        call read_data(file_name,'t_flux',   iob%t_flux,   Ocean%Domain)
        call read_data(file_name,'q_flux',   iob%q_flux,   Ocean%Domain)
        call read_data(file_name,'salt_flux',iob%salt_flux,Ocean%Domain)
        call read_data(file_name,'lw_flux',  iob%lw_flux,  Ocean%Domain)
        call read_data(file_name,'sw_flux',  iob%sw_flux,  Ocean%Domain)
        call read_data(file_name,'lprec',    iob%lprec,    Ocean%Domain)
        call read_data(file_name,'fprec',    iob%fprec,    Ocean%Domain)
        call read_data(file_name,'runoff',   iob%runoff,   Ocean%Domain)
        call read_data(file_name,'calving',  iob%calving,  Ocean%Domain)
        call read_data(file_name,'p',        iob%p,        Ocean%Domain)
    else
        iob%u_flux=0.0
        iob%v_flux=0.0
        iob%t_flux=0.0
        iob%q_flux=0.0
        iob%salt_flux=0.0
        iob%lw_flux=0.0
        iob%sw_flux=0.0
        iob%lprec=0.0
        iob%fprec=0.0
        iob%runoff=0.0
        iob%calving=0.0
        iob%p=0.0
    endif
    

! iob%xtype set by flux_exchange_init             

  end subroutine read_ice_ocean_boundary
! </SUBROUTINE> NAME="read_ice_ocean_boundary"


!#######################################################################
! <SUBROUTINE NAME="init_default_ice_ocean_boundary">
!
! <DESCRIPTION>
! Default surface boundary conditions for use in coupled modeling.
! </DESCRIPTION>
!
  subroutine init_default_ice_ocean_boundary(iob)

    type(ice_ocean_boundary_type),intent(INOUT) :: iob

    iob%u_flux = 0.0
    iob%v_flux = 0.0
    iob%t_flux = 0.0
    iob%q_flux = 0.0
    iob%salt_flux = 0.0
    iob%lw_flux = 0.0
    iob%sw_flux = 0.0
    iob%lprec = 0.0
    iob%fprec = 0.0
    iob%runoff = 0.0
    iob%calving = 0.0
    iob%p = 0.0

! iob%xtype set by flux_exchange_init

  end subroutine init_default_ice_ocean_boundary
! </SUBROUTINE> NAME="read_ice_ocean_boundary"


!#######################################################################
! <SUBROUTINE NAME="time_filter">
!
! <DESCRIPTION>
! When threelevel time tendency is used, perform Robert-Asselin time 
! filtering on velocity and tracer.
!
! Also perform time filter on surface height by replacing eta_t(tau)
! with eta_t_bar(tau).
!
! </DESCRIPTION>
!
  subroutine time_filter(Time, Grid, Thickness, T_prog, Velocity, Ext_mode)

    type(ocean_time_type), intent(in)             :: Time
    type(ocean_grid_type), intent(in)             :: Grid
    type(ocean_thickness_type), intent(inout)     :: Thickness
    type(ocean_prog_tracer_type), intent(inout)   :: T_prog(:)
    type(ocean_velocity_type), intent(inout)      :: Velocity
    type(ocean_external_mode_type), intent(inout) :: Ext_mode

    real, dimension(isd:ied,jsd:jed) :: tmp              
    integer                          :: i, j, k, n
    integer                          :: tau, taum1, taup1

    taum1 = Time%taum1
    tau   = Time%tau
    taup1 = Time%taup1    


!---surface height and thickness of surface cell 

    tmp=0.0
    if(id_tfilter_etat > 0) then 
        do j=jsc,jec
           do i=isc,iec
              tmp(i,j) = Ext_mode%eta_t_bar(i,j,tau) - Ext_mode%eta_t(i,j,tau)
           enddo
        enddo
        used = send_data(id_tfilter_etat, tmp(:,:), &
               Time%model_time, rmask=Grid%tmask(:,:,1), &
               is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    endif
    if(id_tfilter_etau > 0) then 
        do j=jsc,jec
           do i=isc,iec
              tmp(i,j) = Ext_mode%eta_u(i,j,tau)
           enddo
        enddo
    endif

    Ext_mode%eta_t(:,:,tau) = Ext_mode%eta_t_bar(:,:,tau)
    Ext_mode%eta_u(:,:,tau) = Grid%umask(:,:,1)*REMAP_BT_TO_BU(Ext_mode%eta_t(:,:,tau))
    call mpp_update_domains (Ext_mode%eta_u(:,:,tau), Domain%domain2d)
    Thickness%dht(:,:,1,tau) = Grid%dzt(1) + Ext_mode%eta_t(:,:,tau)
    Thickness%dhu(:,:,1,tau) = Grid%dzt(1) + Ext_mode%eta_u(:,:,tau)

    if(id_tfilter_etau > 0) then 
        do j=jsc,jec
           do i=isc,iec
              tmp(i,j) = Ext_mode%eta_u(i,j,tau) - tmp(i,j)
           enddo
        enddo
        used = send_data(id_tfilter_etau, tmp(:,:), &
               Time%model_time, rmask=Grid%umask(:,:,1), &
               is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    endif


!---tracers 

    wrk1 = 0.0
    do n=1,num_prog_tracers

       if(id_robert_tracer(n) > 0) then 
           do k=1,nk
              do j=jsc,jec
                 do i=isc,iec
                    wrk1(i,j,k) = T_prog(n)%field(i,j,k,tau)
                 enddo
              enddo
           enddo
       endif

       do k=1,nk
          do j=jsd,jed
             do i=isd,ied
                T_prog(n)%field(i,j,k,tau) = T_prog(n)%field(i,j,k,tau) + &
                     robert*(0.5*(T_prog(n)%field(i,j,k,taup1) + T_prog(n)%field(i,j,k,taum1)) - T_prog(n)%field(i,j,k,tau))
             enddo
          enddo
       enddo

       if(id_robert_tracer(n) > 0) then 
           do k=1,nk
              do j=jsc,jec
                 do i=isc,iec
                    wrk1(i,j,k) = (T_prog(n)%field(i,j,k,tau)-wrk1(i,j,k)) &
                                   *p5dttsr*Thickness%dht(i,j,k,tau)*T_prog(n)%conversion 
                 enddo
              enddo
           enddo
           used = send_data(id_robert_tracer(n), wrk1(:,:,:), &
                  Time%model_time, rmask=Grid%tmask(:,:,:), &
                  is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
       endif

    enddo

!---velocity 

    wrk1_v=0.0 

    if(id_robert_velocity(1) > 0 .or. id_robert_velocity(2) > 0) then 
        do n=1,2
           do k=1,nk
              do j=jsc,jec
                 do i=isc,iec
                    wrk1_v(i,j,k,n) = Velocity%u(i,j,k,n,tau)
                 enddo
              enddo
           enddo
        enddo
    endif

    do n=1,2
       do k=1,nk
          do j=jsd,jed
             do i=isd,ied
                Velocity%u(i,j,k,n,tau) = Velocity%u(i,j,k,n,tau) + &
                     robert*(0.5*(Velocity%u(i,j,k,n,taup1) + Velocity%u(i,j,k,n,taum1)) - Velocity%u(i,j,k,n,tau))
             enddo
          enddo
       enddo
    enddo

    if(id_robert_velocity(1) > 0 .or. id_robert_velocity(2) > 0) then 
        do n=1,2
           do k=1,nk
              do j=jsc,jec
                 do i=isc,iec
                    wrk1_v(i,j,k,n) = p5dtuvr*(Velocity%u(i,j,k,n,tau)-wrk1_v(i,j,k,n))
                 enddo
              enddo
           enddo
        enddo
        if(id_robert_velocity(1) > 0) then 
           used = send_data(id_robert_velocity(1), wrk1_v(:,:,:,1), &
                  Time%model_time, rmask=Grid%umask(:,:,:), &
                  is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
        endif 
        if(id_robert_velocity(2) > 0) then 
           used = send_data(id_robert_velocity(2), wrk1_v(:,:,:,2), &
                  Time%model_time, rmask=Grid%umask(:,:,:), &
                  is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
        endif 
    endif

    ! check leap-frog noise remaining after Robert filter applied 
    if(id_forward_u > 0) then 
        do k=1,nk
           do j=jsc,jec
              do i=isc,iec
                 wrk1_v(i,j,k,1) = Velocity%u(i,j,k,1,tau)-Velocity%u(i,j,k,1,taum1)
              enddo
           enddo
        enddo
        used = send_data(id_forward_u, wrk1_v(:,:,:,1), &
               Time%model_time, rmask=Grid%umask(:,:,:), &
               is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
    endif
    
    if(id_forward_v > 0 ) then 
        do k=1,nk
           do j=jsc,jec
              do i=isc,iec
                 wrk1_v(i,j,k,2) = Velocity%u(i,j,k,2,tau)-Velocity%u(i,j,k,2,taum1)
              enddo
           enddo
        enddo
        used = send_data(id_forward_v, wrk1_v(:,:,:,2), &
               Time%model_time, rmask=Grid%umask(:,:,:), &
               is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
    endif

    if (velocity_change_check .and. velocity_change_diag_freq > 0 .and. &
       mod(Time%itt, velocity_change_diag_freq) == 0) then 
      call velocity_change(Time, Velocity, velocity_change_max, velocity_change_max_num)
    endif 


  end subroutine time_filter
! </SUBROUTINE> NAME="time_filter"


  
end module ocean_model_mod
  
  
 
