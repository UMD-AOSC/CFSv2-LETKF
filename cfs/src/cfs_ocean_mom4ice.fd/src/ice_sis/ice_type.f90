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
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_type_mod - maintains the sea ice data, reads/writes restarts, reads the  !
!                namelist and initializes diagnostics. - Mike Winton           !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Modified: Xingren Wu
!           Xingren.Wu@noaa.gov
! Modified: Shrinivas Moorthi
!           Shrinivas.Moorthi@noaa.gov
module ice_type_mod

  use mpp_mod,          only: mpp_pe, mpp_root_pe, mpp_sum, mpp_clock_id, CLOCK_COMPONENT
  use mpp_domains_mod,  only: domain2D, mpp_update_domains
  use fms_mod,          only: file_exist, open_namelist_file, check_nml_error
  use fms_mod,          only: read_data, write_data, close_file
  use fms_mod,          only: stdlog, error_mesg, FATAL, WARNING, clock_flag_default
  use diag_manager_mod, only: diag_axis_init, register_diag_field
  use diag_manager_mod, only: register_static_field, send_data
  use time_manager_mod, only: time_type, get_time
  use constants_mod,    only: Tfreeze
  use ice_grid_mod,     only: set_ice_grid, t_to_uv, ice_grid_end
  use ice_grid_mod,     only: Domain, isc, iec, jsc, jec, isd, ied, jsd, jed, im, jm, km
  use ice_grid_mod,     only: geo_lon, geo_lat, cell_area, sin_rot, cos_rot, wett, xb1d, yb1d
  use ice_thm_mod,      only: ice_thm_param
  use ice_dyn_mod,      only: ice_dyn_param
  use time_manager_mod, only: get_date

  implicit none
  private

public :: ice_data_type, ice_model_init, ice_model_end, kmelt, mom_rough_ice, &
          heat_rough_ice, atmos_winds, hlim, slab_ice, spec_ice, verbose,     &
          ice_bulk_salin, do_ice_restore, do_ice_limit, max_ice_limit,        &
          ice_restore_timescale, do_init, h2o, heat,                          &
          cm2_bugs, conservation_check, ice_model_rstrt

public  :: id_cn, id_hi, id_hs, id_t1, id_t2, id_ts
public  :: id_mi, id_sh, id_lh, id_sw, id_lw, id_snofl, id_rain, id_runoff,    &
           id_calving, id_evap, id_saltf, id_tmelt, id_bmelt, id_bheat, id_e2m,&
           id_frazil, id_alb, id_xprt, id_lsrc, id_lsnk, id_bsnk, id_strna,    &
           id_sigi, id_sigii, id_stren, id_ui, id_vi, id_fax, id_fay, id_fix,  &
           id_fiy, id_fcx, id_fcy, id_fwx, id_fwy, id_swdn, id_lwdn, id_sn2ic, &
           id_slp, id_ext, id_sst, id_sss, id_ssh, id_uo, id_vo, id_ta, id_obi,&
           id_qfres, id_qflim, id_ix_trans, id_iy_trans

public  :: iceClock

  !---- id for diagnositics -------------------
  integer :: id_xb, id_xt, id_yb, id_yt, id_ct, id_xv, id_yv
  integer :: id_cn, id_hi, id_hs, id_t1, id_t2, id_ts
  integer :: id_mi, id_sh, id_lh, id_sw, id_lw, id_snofl, id_rain, id_runoff
  integer :: id_calving, id_evap, id_saltf, id_tmelt, id_bmelt, id_bheat, id_e2m
  integer :: id_frazil, id_alb, id_xprt, id_lsrc, id_lsnk, id_bsnk, id_strna
  integer :: id_sigi, id_sigii, id_stren, id_ui, id_vi, id_fax, id_fay, id_fix
  integer :: id_fiy, id_fcx, id_fcy, id_fwx, id_fwy, id_swdn, id_lwdn, id_sn2ic
  integer :: id_slp, id_ext, id_sst, id_sss, id_ssh, id_uo, id_vo, id_ta, id_obi
  integer :: id_qfres, id_qflim, id_ix_trans, id_iy_trans
  integer :: id_sw_vis, id_sw_dir, id_sw_dif, id_sw_vis_dir, id_sw_vis_dif

  !--- namelist interface --------------
  real    :: mom_rough_ice  = 1.0e-4     ! momentum same, cd10=(von_k/ln(10/z0))^2
  real    :: heat_rough_ice = 1.0e-4     ! heat roughness length
  real    :: kmelt          = 6e-5*4e6   ! ocean/ice heat flux constant
  real    :: alb_sno        = 0.85       ! snow albedo (less if melting)
  real    :: alb_ice        = 0.5826     ! ice albedo (less if melting)
  real    :: pen_ice        = 0.3        ! part unreflected solar penetrates ice
  real    :: opt_dep_ice    = 0.67       ! ice optical depth
  real    :: t_range_melt   = 1.0        ! melt albedos scaled in over T range
  real    :: ice_bulk_salin = 0.0        ! ice bulk salinity (for ocean salt flux)
  real    :: p0             = 2.75e4     ! ice strength parameter
  real    :: c0             = 20.0       ! another ice strength parameter
  real    :: cdw            = 3.24e-3    ! water/ice drag coefficient
  real    :: wd_turn        = 25.0       ! water/ice drag turning angle
  integer :: nsteps_dyn     = 432        ! dynamics steps per slow timestep
  integer :: nsteps_adv     = 8          ! advection steps per slow timestep
  integer :: num_part       = 6          ! number of ice grid partitions
                                         ! partition 1 is open water
                                         ! partitions 2 to num_part-1 are
                                         !   thickness limited ice categories
                                         ! partition num_part is unlimited ice
  logical :: atmos_winds = .true.        ! wind stress from atmosphere model over t points and has wrong sign
  logical :: slab_ice    = .false.       ! do old-style GFDL slab ice?
  logical :: spec_ice    = .false.       ! old-style GFDL slab ice with SST, ice thickness and conc. from data
  logical :: do_ice_restore  = .false.   ! restore sea-ice toward climatology
  logical :: do_ice_limit    = .false.   ! limit sea ice to max_ice_limit
  real    :: max_ice_limit   = 4.0       ! maximum sea ice height(m),
                                         ! if do_ice_limit is true
                                         ! TK: default chosen based on observed
                                         !     ice thickness data used by climate
                                         !     group, which range up to 7.31 m
  real    :: ice_restore_timescale = 5.0 ! time scale for restoring ice (days)
  logical :: conservation_check = .true. ! check for heat and h2o conservation
  logical :: cm2_bugs           = .false.! keep cm2 bugs for reproducibility        
  logical :: verbose            = .false.! control printing message, will slow model down when turn true
  integer :: layout(2)          = (/0, 0/)

  namelist /ice_model_nml/ mom_rough_ice, heat_rough_ice, p0, c0, cdw, wd_turn,  &
                           kmelt, alb_sno, alb_ice, pen_ice, opt_dep_ice,        &
                           nsteps_dyn, nsteps_adv, num_part, atmos_winds,        &
                           slab_ice, spec_ice, ice_bulk_salin, layout,           &
                           do_ice_restore, do_ice_limit, max_ice_limit,          &
                           ice_restore_timescale, conservation_check,            &
                           t_range_melt, cm2_bugs, verbose

  logical :: do_init = .false.
  real    :: hlim(8) = (/ 0.0, 0.1, 0.3, 0.7, 1.1, 1.5, 2.0, 2.5 /) ! thickness limits 1...num_part-1
  real    :: h2o(4), heat(4) ! for conservation analysis
                             ! 1 - initial ice h2o/heat content
                             ! 2 - h2o/heat flux down at top of ice
                             ! 3 - h2o/heat flux down at bottom of ice
                             ! 4 - final ice h2o/heat content

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! This structure contains the ice model data (some used by calling routines);  !
  ! the third index is partition (1 is open water; 2 is ice cover)               !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  type ice_data_type
     type(domain2D)                     :: Domain
     type (time_type)                   :: Time_Init, Time
     type (time_type)                   :: Time_step_fast, Time_step_slow
     integer                            :: avg_kount
     logical, pointer, dimension(:,:)   :: mask                =>NULL() ! where ice can be
     logical, pointer, dimension(:,:,:) :: ice_mask            =>NULL() ! where ice actually is
     real,    pointer, dimension(:,:,:) :: fice                =>NULL()
     real,    pointer, dimension(:,:,:) :: part_size           =>NULL()
     real,    pointer, dimension(:,:,:) :: part_size_uv        =>NULL()
     real,    pointer, dimension(:,:,:) :: albedo              =>NULL()
     real,    pointer, dimension(:,:,:) :: albedo_vis_dir      =>NULL()
     real,    pointer, dimension(:,:,:) :: albedo_nir_dir      =>NULL()
     real,    pointer, dimension(:,:,:) :: albedo_vis_dif      =>NULL()
     real,    pointer, dimension(:,:,:) :: albedo_nir_dif      =>NULL()
     real,    pointer, dimension(:,:,:) :: rough_mom           =>NULL()
     real,    pointer, dimension(:,:,:) :: rough_heat          =>NULL()
     real,    pointer, dimension(:,:,:) :: rough_moist         =>NULL()
     real,    pointer, dimension(:,:,:) :: t_surf              =>NULL()
     real,    pointer, dimension(:,:,:) :: u_surf              =>NULL()
     real,    pointer, dimension(:,:,:) :: v_surf              =>NULL()
     real,    pointer, dimension(:,:)   :: sea_lev             =>NULL()
     real,    pointer, dimension(:,:)   :: s_surf              =>NULL()
     real,    pointer, dimension(:,:)   :: u_ocn               =>NULL()
     real,    pointer, dimension(:,:)   :: v_ocn               =>NULL()
     real,    pointer, dimension(:,:,:) :: flux_u_top          =>NULL()
     real,    pointer, dimension(:,:,:) :: flux_v_top          =>NULL()
     real,    pointer, dimension(:,:,:) :: flux_t_top          =>NULL()
     real,    pointer, dimension(:,:,:) :: flux_q_top          =>NULL()
     real,    pointer, dimension(:,:,:) :: flux_sw_top         =>NULL()
     real,    pointer, dimension(:,:,:) :: flux_lw_top         =>NULL()
     real,    pointer, dimension(:,:,:) :: flux_au_top         =>NULL()
     real,    pointer, dimension(:,:,:) :: flux_av_top         =>NULL()
     real,    pointer, dimension(:,:,:) :: hice_top            =>NULL()
     real,    pointer, dimension(:,:,:) :: hsno_top            =>NULL()
     real,    pointer, dimension(:,:,:) :: flux_at_top         =>NULL()
     real,    pointer, dimension(:,:,:) :: flux_aq_top         =>NULL()
     real,    pointer, dimension(:,:,:) :: flux_asw_top        =>NULL()
     real,    pointer, dimension(:,:,:) :: flux_alw_top        =>NULL()
     real,    pointer, dimension(:,:,:) :: ta1_ocn_top         =>NULL()
     real,    pointer, dimension(:,:,:) :: flux_sw_vis_top     =>NULL()
     real,    pointer, dimension(:,:,:) :: flux_sw_dir_top     =>NULL()
     real,    pointer, dimension(:,:,:) :: flux_sw_dif_top     =>NULL()
     real,    pointer, dimension(:,:,:) :: flux_sw_vis_dir_top =>NULL()
     real,    pointer, dimension(:,:,:) :: flux_sw_vis_dif_top =>NULL()
     real,    pointer, dimension(:,:,:) :: flux_lh_top         =>NULL()
     real,    pointer, dimension(:,:,:) :: lprec_top           =>NULL()
     real,    pointer, dimension(:,:,:) :: fprec_top           =>NULL()
     real,    pointer, dimension(:,:  ) :: flux_u              =>NULL()
     real,    pointer, dimension(:,:  ) :: flux_v              =>NULL()
     real,    pointer, dimension(:,:  ) :: flux_t              =>NULL()
     real,    pointer, dimension(:,:  ) :: flux_q              =>NULL()
     real,    pointer, dimension(:,:  ) :: flux_sw             =>NULL()
     real,    pointer, dimension(:,:  ) :: flux_lw             =>NULL()
     real,    pointer, dimension(:,:  ) :: flux_au             =>NULL()
     real,    pointer, dimension(:,:  ) :: flux_av             =>NULL()
     real,    pointer, dimension(:,:  ) :: flux_at             =>NULL()
     real,    pointer, dimension(:,:  ) :: flux_aq             =>NULL()
     real,    pointer, dimension(:,:  ) :: flux_asw            =>NULL()
     real,    pointer, dimension(:,:  ) :: flux_alw            =>NULL()
     real,    pointer, dimension(:,:  ) :: ta1_ocn             =>NULL()
     real,    pointer, dimension(:,:  ) :: flux_sw_vis         =>NULL()
     real,    pointer, dimension(:,:  ) :: flux_sw_dir         =>NULL()
     real,    pointer, dimension(:,:  ) :: flux_sw_dif         =>NULL()
     real,    pointer, dimension(:,:  ) :: flux_sw_vis_dir     =>NULL()
     real,    pointer, dimension(:,:  ) :: flux_sw_vis_dif     =>NULL()
     real,    pointer, dimension(:,:  ) :: flux_lh             =>NULL()
     real,    pointer, dimension(:,:  ) :: lprec               =>NULL()
     real,    pointer, dimension(:,:  ) :: fprec               =>NULL()
     real,    pointer, dimension(:,:  ) :: p_surf              =>NULL()
     real,    pointer, dimension(:,:  ) :: runoff              =>NULL()
     real,    pointer, dimension(:,:  ) :: calving             =>NULL()
     real,    pointer, dimension(:,:  ) :: flux_salt           =>NULL()
     real,    pointer, dimension(:,:)   :: lwdn                =>NULL()
     real,    pointer, dimension(:,:  ) :: swdn                =>NULL() ! downward long/shortwave
     real,    pointer, dimension(:,:,:) :: pen                 =>NULL()
     real,    pointer, dimension(:,:,:) :: trn                 =>NULL() ! ice optical parameters
     real,    pointer, dimension(:,:,:) :: tmelt               =>NULL()
     real,    pointer, dimension(:,:,:) :: bmelt               =>NULL()
     real,    pointer, dimension(:,:,:) :: h_snow              =>NULL()
     real,    pointer, dimension(:,:,:) :: h_ice               =>NULL()
     real,    pointer, dimension(:,:,:) :: t_ice1              =>NULL()
     real,    pointer, dimension(:,:,:) :: t_ice2              =>NULL()
     real,    pointer, dimension(:,:)   :: u_ice               =>NULL()
     real,    pointer, dimension(:,:)   :: v_ice               =>NULL()
     real,    pointer, dimension(:,:)   :: sig11               =>NULL()
     real,    pointer, dimension(:,:)   :: sig22               =>NULL()
     real,    pointer, dimension(:,:)   :: sig12               =>NULL()
     real,    pointer, dimension(:,:)   :: frazil              =>NULL()
     real,    pointer, dimension(:,:)   :: bheat               =>NULL()
     real,    pointer, dimension(:,:)   :: qflx_lim_ice        =>NULL()
     real,    pointer, dimension(:,:)   :: qflx_res_ice        =>NULL()
  end type ice_data_type

  integer :: iceClock

  contains

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! ice_model_init - initializes ice model data, parameters and diagnostics      !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine ice_model_init (Ice, Time_Init, Time, Time_step_fast, Time_step_slow )

    type (ice_data_type), intent(inout) :: Ice
    type (time_type)    , intent(in)    :: Time_Init      ! starting time of model integration
    type (time_type)    , intent(in)    :: Time           ! current time
    type (time_type)    , intent(in)    :: Time_step_fast ! time step for the ice_model_fast
    type (time_type)    , intent(in)    :: Time_step_slow ! time step for the ice_model_slow

    integer           :: io, ierr, nlon, nlat, npart, unit, log_unit, k
    integer           :: sc, dy, i, j 
    real              :: dt_slow
    character(len=22) :: restart = 'INPUT/ice_model.res.nc'
    !
    ! read namelist and write to logfile
    !
    unit = open_namelist_file()
    read  (unit, ice_model_nml,iostat=io)
    write (stdlog(), ice_model_nml)
    ierr = check_nml_error(io, 'ice_model_nml')
    call close_file(unit)

    if (spec_ice) then
       slab_ice = .true.
       nsteps_dyn = 0
       nsteps_adv = 0
    end if
    if (slab_ice) num_part = 2 ! open water and ice ... but never in same place
    if (num_part>size(hlim(:))+1) &
         call error_mesg ('ice_model_init', 'not enough thickness limits', FATAL)

    call get_time(Time_step_slow, sc, dy); dt_slow=864e2*dy+sc

    call set_ice_grid(Ice%domain, dt_slow, nsteps_dyn, nsteps_adv, num_part, layout )

    call ice_dyn_param(p0, c0, cdw, wd_turn, slab_ice)
    call ice_thm_param(alb_sno, alb_ice, pen_ice, opt_dep_ice, slab_ice, t_range_melt)

    allocate ( Ice % mask     (isc:iec, jsc:jec)       , &
         Ice % ice_mask       (isc:iec, jsc:jec, km)   , &
         Ice % t_surf         (isc:iec, jsc:jec, km)   , &
         Ice % s_surf         (isc:iec, jsc:jec)       , &
         Ice % sea_lev        (isd:ied, jsd:jed)       , &
         Ice % fice           (isd:ied, jsd:jed, km)   , &
         Ice % part_size      (isd:ied, jsd:jed, km)   , &
         Ice % part_size_uv   (isc:iec, jsc:jec, km)   , &
         Ice % u_surf         (isc:iec, jsc:jec, km)   , &
         Ice % v_surf         (isc:iec, jsc:jec, km)   , &
         Ice % u_ocn          (isd:ied, jsd:jed)       , &
         Ice % v_ocn          (isd:ied, jsd:jed)       , &
         Ice % rough_mom      (isc:iec, jsc:jec, km)   , &
         Ice % rough_heat     (isc:iec, jsc:jec, km)   , &
         Ice % rough_moist    (isc:iec, jsc:jec, km)   , &
         Ice % albedo         (isc:iec, jsc:jec, km)   , &                
         Ice % albedo_vis_dir (isc:iec, jsc:jec, km)   , &
         Ice % albedo_nir_dir (isc:iec, jsc:jec, km)   , &
         Ice % albedo_vis_dif (isc:iec, jsc:jec, km)   , &
         Ice % albedo_nir_dif (isc:iec, jsc:jec, km)   )

    allocate ( Ice % flux_u_top   (isd:ied, jsd:jed, km) ,       &
         Ice % flux_v_top         (isd:ied, jsd:jed, km) ,       &
         Ice % flux_t_top         (isc:iec, jsc:jec, km) ,       &
         Ice % flux_q_top         (isc:iec, jsc:jec, km) ,       &
         Ice % flux_sw_top        (isc:iec, jsc:jec, km) ,       &
         Ice % hice_top           (isc:iec, jsc:jec, km) ,       &
         Ice % hsno_top           (isc:iec, jsc:jec, km) ,       &
         Ice % flux_au_top        (isd:ied, jsd:jed, km) ,       &
         Ice % flux_av_top        (isd:ied, jsd:jed, km) ,       &
         Ice % flux_at_top        (isc:iec, jsc:jec, km) ,       &
         Ice % flux_aq_top        (isc:iec, jsc:jec, km) ,       &
         Ice % flux_asw_top       (isc:iec, jsc:jec, km) ,       &
         Ice % flux_alw_top       (isc:iec, jsc:jec, km) ,       &
         Ice % ta1_ocn_top        (isc:iec, jsc:jec, km) ,       &
         Ice % flux_sw_vis_top    (isc:iec, jsc:jec, km) ,       &
         Ice % flux_sw_dir_top    (isc:iec, jsc:jec, km) ,       &
         Ice % flux_sw_dif_top    (isc:iec, jsc:jec, km) ,       &
         Ice % flux_sw_vis_dir_top(isc:iec, jsc:jec, km) ,       &
         Ice % flux_sw_vis_dif_top(isc:iec, jsc:jec, km) ,       &
         Ice % flux_lw_top        (isc:iec, jsc:jec, km) ,       &
         Ice % flux_lh_top        (isc:iec, jsc:jec, km) ,       &
         Ice % lprec_top          (isc:iec, jsc:jec, km) ,       &
         Ice % fprec_top          (isc:iec, jsc:jec, km)   )

    allocate ( Ice % flux_u    (isc:iec, jsc:jec ) ,       &
         Ice % flux_v          (isc:iec, jsc:jec ) ,       &
         Ice % flux_t          (isc:iec, jsc:jec ) ,       &
         Ice % flux_q          (isc:iec, jsc:jec ) ,       &
         Ice % flux_sw         (isc:iec, jsc:jec ) ,       &
         Ice % flux_au         (isc:iec, jsc:jec ) ,       &
         Ice % flux_av         (isc:iec, jsc:jec ) ,       &
         Ice % flux_at         (isc:iec, jsc:jec ) ,       &
         Ice % flux_aq         (isc:iec, jsc:jec ) ,       &
         Ice % flux_asw        (isc:iec, jsc:jec ) ,       &
         Ice % flux_alw        (isc:iec, jsc:jec ) ,       &
         Ice % ta1_ocn         (isc:iec, jsc:jec ) ,       &
         Ice % flux_sw_vis     (isc:iec, jsc:jec) ,        &
         Ice % flux_sw_dir     (isc:iec, jsc:jec) ,        &
         Ice % flux_sw_dif     (isc:iec, jsc:jec) ,        &
         Ice % flux_sw_vis_dir (isc:iec, jsc:jec) ,        &
         Ice % flux_sw_vis_dif (isc:iec, jsc:jec) ,        &
         Ice % flux_lw         (isc:iec, jsc:jec ) ,       &
         Ice % flux_lh         (isc:iec, jsc:jec ) ,       &
         Ice % lprec           (isc:iec, jsc:jec ) ,       &
         Ice % fprec           (isc:iec, jsc:jec ) ,       &
         Ice % p_surf          (isc:iec, jsc:jec ) ,       &
         Ice % runoff          (isc:iec, jsc:jec ) ,       &
         Ice % calving         (isc:iec, jsc:jec ) ,       &
         Ice % flux_salt       (isc:iec, jsc:jec ) ,       &
         Ice % lwdn            (isc:iec, jsc:jec ) ,       &
         Ice % swdn            (isc:iec, jsc:jec )         )
    allocate ( Ice % frazil (isc:iec, jsc:jec), Ice % bheat  (isc:iec, jsc:jec), &
               Ice % u_ice  (isd:ied, jsd:jed), Ice % v_ice  (isd:ied, jsd:jed), &
               Ice % sig11  (isd:ied, jsd:jed), Ice % sig22  (isd:ied, jsd:jed), &
               Ice % sig12  (isd:ied, jsd:jed)                               )
    allocate ( Ice % tmelt  (isc:iec, jsc:jec, 2:km), Ice % bmelt  (isc:iec, jsc:jec, 2:km) , &
               Ice % pen    (isc:iec, jsc:jec, 2:km), Ice % trn    (isc:iec, jsc:jec, 2:km) , &
               Ice % h_snow (isd:ied, jsd:jed, 2:km), Ice % h_ice  (isd:ied, jsd:jed, 2:km) , &
               Ice % t_ice1 (isd:ied, jsd:jed, 2:km), Ice % t_ice2 (isd:ied, jsd:jed, 2:km)   )
    allocate ( Ice % qflx_lim_ice  (isc:iec, jsc:jec) , Ice % qflx_res_ice  (isc:iec, jsc:jec)   )

    Ice % flux_sw_vis     =0.
    Ice % flux_sw_dir     =0.
    Ice % flux_sw_dif     =0.
    Ice % flux_sw_vis_dir =0.
    Ice % flux_sw_vis_dif =0.
    Ice % flux_lh         =0. 
    Ice % lwdn            =0.
    Ice % swdn            =0.
    Ice % flux_u_top      =0. 
    Ice % flux_v_top      =0.
    Ice % sea_lev         =0.
    Ice % part_size       =0.
    Ice % u_ocn           =0.
    Ice % v_ocn           =0.
    Ice % u_ice           =0.
    Ice % v_ice           =0.
    Ice % sig11           =0.
    Ice % sig12           =0.
    Ice % sig22           =0.
    Ice % h_snow          =0.
    Ice % h_ice           =0.
    Ice % t_ice1          =0.
    Ice % t_ice2          =0.

    do j = jsc, jec
       do i = isc, iec
          Ice % mask(i,j) = wett(i,j)
       enddo
    enddo

    Ice % Time           = Time
    Ice % Time_Init      = Time_Init
    Ice % Time_step_fast = Time_step_fast
    Ice % Time_step_slow = Time_step_slow

    Ice % avg_kount      = 0
    !
    ! read restart
    !
    !! Need to create new restart version including the needed albedos
    !  that have been added ??

    if (file_exist(restart)) then
       !Balaji: netCDF restart via fms_io
       call read_data( restart, 'part_size',   Ice%part_size,        domain )
       call read_data( restart, 'albedo',      Ice%albedo,           domain )
       call read_data( restart, 'rough_mom',   Ice%rough_mom,        domain )
       call read_data( restart, 'rough_heat',  Ice%rough_heat,       domain )
       call read_data( restart, 'rough_moist', Ice%rough_moist,      domain )
       call read_data( restart, 't_surf',      Ice%t_surf,           domain )
       call read_data( restart, 'h_snow',      Ice%h_snow(:,:,2:km), domain )
       call read_data( restart, 'h_ice',       Ice%h_ice (:,:,2:km), domain )
       call read_data( restart, 't_ice1',      Ice%t_ice1(:,:,2:km), domain )
       call read_data( restart, 't_ice2',      Ice%t_ice2(:,:,2:km), domain )
       call read_data( restart, 'u_ice',       Ice%u_ice,            domain )
       call read_data( restart, 'v_ice',       Ice%v_ice,            domain )
       call read_data( restart, 'sig11',       Ice%sig11,            domain )
       call read_data( restart, 'sig22',       Ice%sig22,            domain )
       call read_data( restart, 'sig12',       Ice%sig12,            domain )
       call read_data( restart, 'flux_u',      Ice%flux_u,           domain )
       call read_data( restart, 'flux_v',      Ice%flux_v,           domain )
       call read_data( restart, 'flux_t',      Ice%flux_t,           domain )
       call read_data( restart, 'flux_q',      Ice%flux_q,           domain )
       call read_data( restart, 'flux_salt',   Ice%flux_salt,        domain )
       call read_data( restart, 'flux_sw',     Ice%flux_sw,          domain )
       call read_data( restart, 'flux_lw',     Ice%flux_lw,          domain )
       call read_data( restart, 'lprec',       Ice%lprec,            domain )
       call read_data( restart, 'fprec',       Ice%fprec,            domain )
       call read_data( restart, 'runoff',      Ice%runoff,           domain )
       call read_data( restart, 'calving',     Ice%calving,          domain )
       call read_data( restart, 'p_surf',      Ice%p_surf,           domain )

       !--- update to data domain
       call mpp_update_domains(Ice%part_size, Domain)
       call mpp_update_domains(Ice%h_snow(:,:,2:km), Domain )
       call mpp_update_domains(Ice%h_ice (:,:,2:km), Domain )
       call mpp_update_domains(Ice%t_ice1(:,:,2:km), Domain )
       call mpp_update_domains(Ice%t_ice2(:,:,2:km), Domain )
       call mpp_update_domains(Ice%u_ice, Domain )
       call mpp_update_domains(Ice%v_ice, Domain )
       call mpp_update_domains(Ice%sig11, Domain )
       call mpp_update_domains(Ice%sig22, Domain )
       call mpp_update_domains(Ice%sig12, Domain )
    else ! no restart => no ice
       Ice % part_size    = 0.0
       !   where (Ice%mask) Ice % part_size (:,:,1) = 1.0  - flux_exchange won't allow
       Ice % part_size (:,:,1) = 1.0
       Ice % albedo            = 0.0
       Ice % rough_mom         = mom_rough_ice
       Ice % rough_heat        = heat_rough_ice
       Ice % rough_moist       = heat_rough_ice
       Ice % t_surf            = Tfreeze-5.0
       Ice % h_snow            = 0.0
       Ice % h_ice             = 0.0
       Ice % t_ice1            = -5.0
       Ice % t_ice2            = -5.0
       Ice % u_ice             = 0.0
       Ice % v_ice             = 0.0
       Ice % sig11             = 0.0
       Ice % sig22             = 0.0
       Ice % sig12             = 0.0
       Ice % flux_u            = 0.0 
       Ice % flux_v            = 0.0
       Ice % flux_t            = 0.0 
       Ice % flux_q            = 0.0
       Ice % flux_sw           = 0.0 
       Ice % flux_lw           = 0.0
       Ice % flux_salt         = 0.0 
       Ice % lprec             = 0.0
       Ice % fprec             = 0.0
       Ice % p_surf            = 0.0
       Ice % runoff            = 0.0
       Ice % calving           = 0.0
       Ice % frazil            = 0.0
       do_init = .true. ! done in ice_model
    end if

    Ice % hice_top(:,:,1) = 0.0
    Ice % hsno_top(:,:,1) = 0.0
    Ice % flux_au   =  Ice % flux_u
    Ice % flux_av   =  Ice % flux_v
    Ice % flux_at   =  Ice % flux_t
    Ice % flux_aq   =  Ice % flux_q
    Ice % flux_asw  =  Ice % flux_sw
    Ice % flux_alw  =  Ice % flux_lw
    Ice % ta1_ocn   =  273.0

  !! INITIALIZE the dir, dif and vis, nir albedoes appropriately.
!!! THIS NEEDS TO BE DEFINED PROPERLY -- ARE THEY NEEDED ON RESTART ??
!!!  FOR NOW, simply set all values to that of Ice%albedo, In general,
!! the ability to read from restart file and initialize when none
!! present should be incorporated above.
    Ice%albedo_vis_dir = 0.0         
    Ice%albedo_nir_dir = 0.0       
    Ice%albedo_vis_dif = 0.0       
    Ice%albedo_nir_dif = 0.0       

    Ice % tmelt       = 0.0
    Ice % bmelt       = 0.0

    Ice % qflx_lim_ice = 0.0
    Ice % qflx_res_ice = 0.0

    Ice%part_size_uv(:,:,1) = 1.0
    do k=2,km
       call t_to_uv(Ice%part_size(:,:,k),Ice%part_size_uv(:,:,k))
       Ice%part_size_uv (:,:,1) = Ice%part_size_uv(:,:,1)-Ice%part_size_uv (:,:,k)
    end do

    call ice_diagnostics_init(Ice)
    !Balaji
    iceClock = mpp_clock_id( 'Ice', flags=clock_flag_default, grain=CLOCK_COMPONENT )
 
  end subroutine ice_model_init

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! ice_model_rstrt - writes the restart file                                      !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine ice_model_rstrt (Ice)
    type (ice_data_type), intent(inout) :: Ice
    integer           :: k

    integer           :: unit

    integer :: yr, mon, day, hr, min, sec
    character(len=128) :: restart
    character(len=10) :: rdte
!   character(len=22) :: restart='RESTART/ice_model.res'
!
    call get_date(Ice%Time,yr,mon,day,hr,min,sec)
    write(rdte,'(i4,3i2.2)') yr,mon,day,hr

    restart = 'IRESTART/'//rdte//'ice_model.res'
!

    if (conservation_check) then
       do k=1,4
          call mpp_sum(h2o(k))
          call mpp_sum(heat(k))
       end do
       if (mpp_pe()==mpp_root_pe()) then
          print *
          print '(a10,5a13)',   'ICE MODEL ','   AT START  ', &
               ' TOP FLUX DN.', &
               ' BOT FLUX DN.', &
               '   AT END    ', &
               '   ERROR     '
          print '(a10,5es13.5)','WATER     ', h2o , h2o (4)-(h2o (1)+h2o (2)-h2o (3))
          print '(a10,5es13.5)','HEAT      ', heat, heat(4)-(heat(1)+heat(2)-heat(3))
          print *
       end if
    end if
    !Balaji: netCDF restart via fms_io
    call write_data( restart, 'part_size',   Ice%part_size,        domain )
    call write_data( restart, 'albedo',      Ice%albedo,           domain )
    call write_data( restart, 'rough_mom',   Ice%rough_mom,        domain )
    call write_data( restart, 'rough_heat',  Ice%rough_heat,       domain )
    call write_data( restart, 'rough_moist', Ice%rough_moist,      domain )
    call write_data( restart, 't_surf',      Ice%t_surf,           domain )
    call write_data( restart, 'h_snow',      Ice%h_snow(:,:,2:km), domain )
    call write_data( restart, 'h_ice',       Ice%h_ice (:,:,2:km), domain )
    call write_data( restart, 't_ice1',      Ice%t_ice1(:,:,2:km), domain )
    call write_data( restart, 't_ice2',      Ice%t_ice2(:,:,2:km), domain )
    call write_data( restart, 'u_ice',       Ice%u_ice,            domain )
    call write_data( restart, 'v_ice',       Ice%v_ice,            domain )
    call write_data( restart, 'sig11',       Ice%sig11,            domain )
    call write_data( restart, 'sig22',       Ice%sig22,            domain )
    call write_data( restart, 'sig12',       Ice%sig12,            domain )
    call write_data( restart, 'flux_u',      Ice%flux_u,           domain )
    call write_data( restart, 'flux_v',      Ice%flux_v,           domain )
    call write_data( restart, 'flux_t',      Ice%flux_t,           domain )
    call write_data( restart, 'flux_q',      Ice%flux_q,           domain )
    call write_data( restart, 'flux_salt',   Ice%flux_salt,        domain )
    call write_data( restart, 'flux_sw',     Ice%flux_sw,          domain )
    call write_data( restart, 'flux_lw',     Ice%flux_lw,          domain )
    call write_data( restart, 'lprec',       Ice%lprec,            domain )
    call write_data( restart, 'fprec',       Ice%fprec,            domain )
    call write_data( restart, 'runoff',      Ice%runoff,           domain )
    call write_data( restart, 'calving',     Ice%calving,          domain )
    call write_data( restart, 'p_surf',      Ice%p_surf,           domain )

  end subroutine ice_model_rstrt
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! ice_model_end - writes the restart file                                      !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine ice_model_end (Ice)
    type (ice_data_type), intent(inout) :: Ice
    integer           :: k

    integer           :: unit
    character(len=22) :: restart='RESTART/ice_model.res'

    if (conservation_check) then
       do k=1,4
          call mpp_sum(h2o(k))
          call mpp_sum(heat(k))
       end do
       if (mpp_pe()==mpp_root_pe()) then
          print *
          print '(a10,5a13)',   'ICE MODEL ','   AT START  ', &
               ' TOP FLUX DN.', &
               ' BOT FLUX DN.', &
               '   AT END    ', &
               '   ERROR     '
          print '(a10,5es13.5)','WATER     ', h2o , h2o (4)-(h2o (1)+h2o (2)-h2o (3))
          print '(a10,5es13.5)','HEAT      ', heat, heat(4)-(heat(1)+heat(2)-heat(3))
          print *
       end if
    end if
    !Balaji: netCDF restart via fms_io
    call write_data( restart, 'part_size',   Ice%part_size,        domain )
    call write_data( restart, 'albedo',      Ice%albedo,           domain )
    call write_data( restart, 'rough_mom',   Ice%rough_mom,        domain )
    call write_data( restart, 'rough_heat',  Ice%rough_heat,       domain )
    call write_data( restart, 'rough_moist', Ice%rough_moist,      domain )
    call write_data( restart, 't_surf',      Ice%t_surf,           domain )
    call write_data( restart, 'h_snow',      Ice%h_snow(:,:,2:km), domain )
    call write_data( restart, 'h_ice',       Ice%h_ice (:,:,2:km), domain )
    call write_data( restart, 't_ice1',      Ice%t_ice1(:,:,2:km), domain )
    call write_data( restart, 't_ice2',      Ice%t_ice2(:,:,2:km), domain )
    call write_data( restart, 'u_ice',       Ice%u_ice,            domain )
    call write_data( restart, 'v_ice',       Ice%v_ice,            domain )
    call write_data( restart, 'sig11',       Ice%sig11,            domain )
    call write_data( restart, 'sig22',       Ice%sig22,            domain )
    call write_data( restart, 'sig12',       Ice%sig12,            domain )
    call write_data( restart, 'flux_u',      Ice%flux_u,           domain )
    call write_data( restart, 'flux_v',      Ice%flux_v,           domain )
    call write_data( restart, 'flux_t',      Ice%flux_t,           domain )
    call write_data( restart, 'flux_q',      Ice%flux_q,           domain )
    call write_data( restart, 'flux_salt',   Ice%flux_salt,        domain )
    call write_data( restart, 'flux_sw',     Ice%flux_sw,          domain )
    call write_data( restart, 'flux_lw',     Ice%flux_lw,          domain )
    call write_data( restart, 'lprec',       Ice%lprec,            domain )
    call write_data( restart, 'fprec',       Ice%fprec,            domain )
    call write_data( restart, 'runoff',      Ice%runoff,           domain )
    call write_data( restart, 'calving',     Ice%calving,          domain )
    call write_data( restart, 'p_surf',      Ice%p_surf,           domain )

    !--- release memory ------------------------------------------------
    call ice_grid_end()

    deallocate(Ice % mask, Ice % ice_mask, Ice % t_surf, Ice % s_surf, Ice % sea_lev )
    deallocate(Ice % part_size, Ice % part_size_uv, Ice % u_surf, Ice % v_surf )
    deallocate(Ice % u_ocn, Ice % v_ocn ,  Ice % rough_mom, Ice % rough_heat )
    deallocate(Ice % rough_moist, Ice % albedo, Ice % flux_u_top, Ice % flux_v_top )
    deallocate(Ice % flux_t_top, Ice % flux_q_top, Ice % flux_sw_top, Ice % flux_lw_top )
    deallocate(Ice % flux_lh_top, Ice % lprec_top, Ice % fprec_top, Ice % flux_u )
    deallocate(Ice % flux_v, Ice % flux_t, Ice % flux_q, Ice % flux_sw, Ice % flux_lw )
    deallocate(Ice % flux_lh, Ice % lprec, Ice % fprec, Ice % p_surf, Ice % runoff ) 
    deallocate(Ice % hice_top, Ice % hsno_top, Ice % ta1_ocn_top )
    deallocate(Ice % flux_au_top, Ice % flux_av_top, Ice % flux_at_top )
    deallocate(Ice % flux_aq_top, Ice % flux_asw_top, Ice % flux_alw_top )
    deallocate(Ice % flux_au, Ice % flux_av, Ice % flux_at, Ice % flux_aq )
    deallocate(Ice % flux_asw, Ice % flux_alw, Ice % ta1_ocn )
    deallocate(Ice % calving)
    deallocate(Ice % flux_salt)
    deallocate( Ice % lwdn)
    deallocate( Ice % swdn)
    deallocate( Ice % frazil )
    deallocate( Ice % fice )
    deallocate(Ice % bheat, Ice % u_ice, Ice % v_ice, Ice % sig11, Ice % sig22 )
    deallocate(Ice % sig12, Ice % tmelt, Ice % bmelt, Ice % pen, Ice % trn )
    deallocate(Ice % h_snow, Ice % h_ice, Ice % t_ice1, Ice % t_ice2  )
    deallocate(Ice % qflx_lim_ice, Ice % qflx_res_ice )

  end subroutine ice_model_end

  !#######################################################################

  subroutine ice_diagnostics_init(Ice)
    type (ice_data_type), intent(in) :: Ice

    real, parameter       :: missing = -1e34
    integer, dimension(2) :: axt, axv, axtv, axvt
    integer, dimension(3) :: axt2
    integer               :: id_geo_lon, id_geo_lat, id_sin_rot, id_cos_rot, id_cell_area
    logical               :: sent

    !
    ! diagnostics MUST use a domain without halos otherwise same as the
    ! regular domain:  Domain (see ice_grid.f90)
    !
    id_xv = diag_axis_init('xv', xb1d(2:im+1), 'degrees_E', 'X','longitude', set_name='ice', Domain2=Domain )
    id_yv = diag_axis_init('yv', yb1d(2:jm+1), 'degrees_N', 'Y','latitude',  set_name='ice', Domain2=Domain )
    id_xb = diag_axis_init('xb', xb1d, 'degrees_E', 'X', 'longitude', set_name='ice', Domain2=Domain )
    id_yb = diag_axis_init('yb', yb1d, 'degrees_N', 'Y', 'latitude', set_name='ice', Domain2=Domain )
    id_xt = diag_axis_init('xt', (xb1d(1:im)+xb1d(2:im+1))/2, 'degrees_E', 'X', &
            'longitude',set_name='ice',edges=id_xb,Domain2=Domain)
    id_yt = diag_axis_init('yt', (yb1d(1:jm)+yb1d(2:jm+1))/2, 'degrees_N', 'Y', &
            'latitude',set_name='ice', edges=id_yb,Domain2=Domain)
    id_ct = diag_axis_init('ct', hlim(1:num_part-1), 'meters','Z', 'thickness')
    axv  = (/ id_xv, id_yv       /)
    axt  = (/ id_xt, id_yt       /)
    axt2 = (/ id_xt, id_yt, id_ct/)
    axtv = (/ id_xt, id_yv /); ! for north faces of t-cells
    axvt = (/ id_xv, id_yt /); ! for east  faces of t-cells

    ! note: set all require=.true. when diag_manager writes these
    !       fields only to "ice_model"
    id_sin_rot   = register_static_field('ice_model', 'SINROT', axt,              &
                   '-SINROT,COSROT points north', 'none', require=.false.)
    id_cos_rot   = register_static_field('ice_model', 'COSROT', axt,              &
                   'COSROT,SINROT points east','none', require=.false.)
    id_geo_lon   = register_static_field('ice_model', 'GEOLON', axt, 'longitude', &
                   'degrees', require=.false.)
    id_geo_lat   = register_static_field('ice_model', 'GEOLAT', axt, 'latitude',  &
                   'degrees', require=.false.)
    id_cell_area = register_static_field('ice_model', 'CELL_AREA', axt,           &
                   'cell area', 'sphere',require=.false.)
    id_ext       = register_diag_field('ice_model', 'MOI', axt, Ice%Time,   &
                   'ice modeled', '0 or 1', missing_value=missing)
    if (id_ext > 0 ) then
       call error_mesg ('ice_model_init', &
            'Diagnostic MOI has been renamed EXT.  Change your diag_table.', WARNING)
    else
       id_ext = register_diag_field('ice_model', 'EXT', axt, Ice%Time, &
                'ice modeled', '0 or 1', missing_value=missing)
    end if
    id_mi       = register_diag_field('ice_model', 'MI', axt, Ice%Time,                  &
                 'ice mass', 'kg/m^2', missing_value=missing)
    id_cn       = register_diag_field('ice_model', 'CN', axt2, Ice%Time,                 &
                 'ice concentration', '0-1', missing_value=missing)
    id_hs       = register_diag_field('ice_model', 'HS', axt, Ice%Time,                  &
                 'snow thickness', 'm-snow', missing_value=missing)
    id_hi       = register_diag_field('ice_model', 'HI', axt, Ice%Time,                  &
                 'ice thickness', 'm-ice', missing_value=missing)
    id_t1       = register_diag_field('ice_model', 'T1', axt, Ice%Time,                  &
                 'upper ice layer temperature', 'C',  missing_value=missing)
    id_t2       = register_diag_field('ice_model', 'T2', axt, Ice%Time,                  &
                 'lower ice layer temperature', 'C',  missing_value=missing)
    id_ts       = register_diag_field('ice_model', 'TS', axt, Ice%Time,                  &
                 'surface temperature', 'C', missing_value=missing)
    id_sh       = register_diag_field('ice_model','SH' ,axt, Ice%Time,                   &
                 'sensible heat flux', 'W/m^2',  missing_value=missing)
    id_lh       = register_diag_field('ice_model','LH' ,axt, Ice%Time,                   &
                 'latent heat flux', 'W/m^2', missing_value=missing)
    id_sw       = register_diag_field('ice_model','SW' ,axt, Ice%Time,                   &
                 'short wave heat flux', 'W/m^2', missing_value=missing)
    id_lw       = register_diag_field('ice_model','LW' ,axt, Ice%Time,                   &
                 'long wave heat flux over ice', 'W/m^2', missing_value=missing)
    id_snofl    = register_diag_field('ice_model','SNOWFL' ,axt, Ice%Time,               &
                 'rate of snow fall', 'kg/(m^2*s)', missing_value=missing)
    id_rain     = register_diag_field('ice_model','RAIN' ,axt, Ice%Time,                 &
                 'rate of rain fall', 'kg/(m^2*s)', missing_value=missing)
    id_runoff   = register_diag_field('ice_model','RUNOFF' ,axt, Ice%Time,               &
                 'liquid runoff', 'kg/(m^2*s)', missing_value=missing)
    id_calving  = register_diag_field('ice_model','CALVING',axt, Ice%Time,               &
                 'frozen runoff', 'kg/(m^2*s)', missing_value=missing)
    id_evap     = register_diag_field('ice_model','EVAP',axt, Ice%Time,                  &
                 'evaporation', 'kg/(m^2*s)', missing_value=missing)
    id_saltf    = register_diag_field('ice_model','SALTF' ,axt, Ice%Time,                &
                 'ice to ocean salt flux', 'kg/(m^2*s)', missing_value=missing)
    id_sn2ic    = register_diag_field('ice_model','SN2IC'  ,axt,Ice%Time,                &
                 'rate of snow to ice conversion', 'kg/(m^2*s)', missing_value=missing)
    id_tmelt    = register_diag_field('ice_model','TMELT'  ,axt, Ice%Time,               &
                 'upper surface melting energy flux', 'W/m^2', missing_value=missing)
    id_bmelt    = register_diag_field('ice_model','BMELT'  ,axt, Ice%Time,               &
                 'bottom surface melting energy flux', 'W/m^2', missing_value=missing)
    id_bheat    = register_diag_field('ice_model','BHEAT'  ,axt, Ice%Time,               &
                 'ocean to ice heat flux', 'W/m^2', missing_value=missing)
    id_e2m      = register_diag_field('ice_model','E2MELT' ,axt, Ice%Time,               &
                 'heat needed to melt ice', 'J/m^2', missing_value=missing)
    id_frazil   = register_diag_field('ice_model','FRAZIL' ,axt, Ice%Time,               &
                 'energy flux of frazil formation', 'W/m^2', missing_value=missing)
    id_alb      = register_diag_field('ice_model','ALB',axt, Ice%Time,                   &
                 'surface albedo','0-1', missing_value=missing )
    id_xprt     = register_diag_field('ice_model','XPRT',axt, Ice%Time,                  &
                 'frozen water transport convergence', 'kg/(m^2*yr)', missing_value=missing)
    id_lsrc     = register_diag_field('ice_model','LSRC', axt, Ice%Time,                 &
                 'frozen water local source', 'kg/(m^2*yr)', missing_value=missing)
    id_lsnk     = register_diag_field('ice_model','LSNK',axt, Ice%Time,                  &
                 'frozen water local sink', 'kg/(m^2*yr)', missing_value=missing)
    id_bsnk     = register_diag_field('ice_model','BSNK',axt, Ice%Time,                  &
                 'frozen water local bottom sink', 'kg/(m^2*yr)', missing_value=missing)
    id_qfres    = register_diag_field('ice_model', 'QFLX_RESTORE_ICE', axt, Ice%Time,    &
                 'Ice Restoring heat flux', 'W/m^2', missing_value=missing)
    id_qflim    = register_diag_field('ice_model', 'QFLX_LIMIT_ICE', axt, Ice%Time,      &
                 'Ice Limit heat flux', 'W/m^2', missing_value=missing)
    id_strna    = register_diag_field('ice_model','STRAIN_ANGLE', axt,Ice%Time,          &
                 'strain angle', 'none', missing_value=missing)
    id_sigi     = register_diag_field('ice_model','SIGI' ,axt, Ice%Time,                 &
                 'first stress invariant', 'none', missing_value=missing)
    id_sigii    = register_diag_field('ice_model','SIGII' ,axt, Ice%Time,                &
                 'second stress invariant', 'none', missing_value=missing)
    id_stren    = register_diag_field('ice_model','STRENGTH' ,axt, Ice%Time,             &
                 'ice strength', 'Pa*m', missing_value=missing)
    id_ui       = register_diag_field('ice_model', 'UI', axv, Ice%Time,                  &
                 'ice velocity - x component', 'm/s', missing_value=missing)
    id_vi       = register_diag_field('ice_model', 'VI', axv, Ice%Time,                  &
                 'ice velocity - y component', 'm/s', missing_value=missing)
    id_ix_trans = register_diag_field('ice_model', 'IX_TRANS', axvt, Ice%Time,           &
                 'x-direction ice transport', 'kg/s', missing_value=missing)
    id_iy_trans = register_diag_field('ice_model', 'IY_TRANS', axtv, Ice%Time,           &
                 'y-direction ice transport', 'kg/s', missing_value=missing)
    id_fax      = register_diag_field('ice_model', 'FA_X', axv, Ice%Time,                &
                 'air stress on ice - x component', 'Pa', missing_value=missing)
    id_fay      = register_diag_field('ice_model', 'FA_Y', axv, Ice%Time,                &
                 'air stress on ice - y component', 'Pa', missing_value=missing)
    id_fix      = register_diag_field('ice_model', 'FI_X', axv, Ice%Time,                &
                 'ice internal stress - x component', 'Pa', missing_value=missing)
    id_fiy      = register_diag_field('ice_model', 'FI_Y', axv, Ice%Time,                &
                 'ice internal stress - y component', 'Pa', missing_value=missing)
    id_fcx      = register_diag_field('ice_model', 'FC_X', axv, Ice%Time,                &
                 'coriolis force - x component', 'Pa', missing_value=missing)
    id_fcy      = register_diag_field('ice_model', 'FC_Y', axv, Ice%Time,                &
                 'coriolis force - y component', 'Pa', missing_value=missing)
    id_fwx      = register_diag_field('ice_model', 'FW_X', axv, Ice%Time,                &
                 'water stress on ice - x component', 'Pa', missing_value=missing)
    id_fwy      = register_diag_field('ice_model', 'FW_Y', axv, Ice%Time,                &
                 'water stress on ice - y component', 'Pa', missing_value=missing)
    id_uo       = register_diag_field('ice_model', 'UO', axv, Ice%Time,                  &
                 'surface current - x component', 'm/s', missing_value=missing)
    id_vo       = register_diag_field('ice_model', 'VO', axv, Ice%Time,                  &
                 'surface current - y component', 'm/s', missing_value=missing)
    !
    ! diagnostics for quantities produced outside the ice model
    !
    id_swdn  = register_diag_field('ice_model','SWDN' ,axt, Ice%Time,       &
               'downward shortwave flux', 'W/m^2', missing_value=missing)
    id_lwdn  = register_diag_field('ice_model','LWDN' ,axt, Ice%Time,       &
               'downward longwave flux', 'W/m^2', missing_value=missing)
    id_ta    = register_diag_field('ice_model', 'TA', axt, Ice%Time,        &
               'surface air temperature', 'C', missing_value=missing)
    id_slp   = register_diag_field('ice_model', 'SLP', axt, Ice%Time,       &
               'sea level pressure', 'Pa', missing_value=missing)
    id_sst   = register_diag_field('ice_model', 'SST', axt, Ice%Time,       &
               'sea surface temperature', 'deg-C', missing_value=missing)
    id_sss   = register_diag_field('ice_model', 'SSS', axt, Ice%Time,       &
               'sea surface salinity', 'psu', missing_value=missing)
    id_ssh   = register_diag_field('ice_model', 'SSH', axt, Ice%Time,       &
               'sea surface height', 'm', missing_value=missing)
    id_obi   = register_diag_field('ice_model', 'OBI', axt, Ice%Time,       &
         'ice observed', '0 or 1', missing_value=missing)

    if (id_sin_rot>0)   sent=send_data(id_sin_rot, sin_rot(isc:iec,jsc:jec), Ice%Time);
    if (id_cos_rot>0)   sent=send_data(id_cos_rot, cos_rot(isc:iec,jsc:jec), Ice%Time);
    if (id_geo_lon>0)   sent=send_data(id_geo_lon, geo_lon, Ice%Time);
    if (id_geo_lat>0)   sent=send_data(id_geo_lat, geo_lat, Ice%Time);
    if (id_cell_area>0) sent=send_data(id_cell_area, cell_area, Ice%Time);

  end subroutine ice_diagnostics_init


end module ice_type_mod
