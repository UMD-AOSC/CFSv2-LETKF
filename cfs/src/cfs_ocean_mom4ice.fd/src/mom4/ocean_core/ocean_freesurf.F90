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
module ocean_freesurf_mod
!
! <CONTACT EMAIL="Stephen.Griffies@noaa.gov"> S.M. Griffies
! </CONTACT>
!
! <CONTACT EMAIL="Ronald.Pacanowski@noaa.gov"> R.C. Pacanowski
! </CONTACT>
!
! <CONTACT EMAIL="Zhi.Liang@noaa.gov"> Zhi Liang (OBC)
! </CONTACT>
!
! <CONTACT EMAIL="martin.schmidt@io-warnemuende.de"> Martin Schmidt (OBC)
! </CONTACT>
!
! <CONTACT EMAIL="Harper.Simmons@noaa.gov"> Harper Simmons (tides)
! </CONTACT>
!
! <REVIEWER EMAIL="Matthew.Harrison@noaa.gov"> M.J. Harrison 
! </REVIEWER>
!
!<OVERVIEW>
! Solve the vertically integrated dynamics using one of two 
! explicit free surface algorithms.
!</OVERVIEW>
!
!<DESCRIPTION>
!
! This module time steps the vertically integrated dynamics, 
! including the surface height and the vertically integrated 
! horizontal velocity. Two explicit time stepping schemes are 
! available:
!
! A. Leap-frog 
! B. Predictor-Corrector with adjustable dissipation
!
! Both use a split-explicit method with an explicit free 
! surface. There is no rigid lid available in mom4. 
!
! </DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! S.M. Griffies, R.C. Pacanowski, R.M. Schmidt, and V. Balaji
! Tracer Conservation with an Explicit Free Surface Method for 
! Z-coordinate Ocean Models
! Monthly Weather Review (2001) vol 129 pages 1081--1098
! </REFERENCE>
!
! <REFERENCE>
! S.M. Griffies
! Fundamentals of Ocean Climate Models
! Princeton University Press (2004)
! </REFERENCE>
!
! <REFERENCE>
! S.M. Griffies, M.J. Harrison, R.C. Pacanowski, and A. Rosati 
! A Technical Guide to MOM4, including supplement for MOM4p1 (2004)
! </REFERENCE>
!
! </INFO>
!<NAMELIST NAME="ocean_freesurf_nml">
!
!  <DATA NAME="zero_tendency" TYPE="logical">
!  If true, will not integrate the free surface.
!  </DATA> 
!  <DATA NAME="zero_eta_ic" TYPE="logical">
!  To initialize eta_t to zero. 
!  </DATA> 
!
!  <DATA NAME="barotropic_leap_frog" TYPE="logical">
!  Use leap-frog scheme for barotropic time stepping. 
!  Not the recommended method, since it requires smaller time 
!  steps.  It is maintained for legacy purposes.  
!  </DATA> 
!  <DATA NAME="robert_asselin" TYPE="logical">
!  Robert time filter for use with leap-frog scheme for barotropic.
!  </DATA> 
!
!  <DATA NAME="barotropic_pred_corr" TYPE="logical">
!  Use preditor-corrector for barotropic time stepping. 
!  This is the recommended method.
!  </DATA> 
!  <DATA NAME="pred_corr_gamma" UNITS="dimensionless" TYPE="real">
!  Dimensionless dissipation parameter for the preditor-corrector
!  scheme.  Setting pred_corr_gamma=0.0 reduces the scheme to a 
!  forward-backward, but it has been found to be unstable.  
!  So pred_corr_gamma > 0.0 is recommended.  Note that 
!  pred_corr_gamma > 0.25 may be over-dissipated and so may 
!  go unstable. 
!  </DATA> 
!
!  <DATA NAME="smooth_eta_t_fs_laplacian" TYPE="logical">
!  For spatially smoothing the eta_t field at each barotropic time step 
!  with a Laplacian operator.  May not be necessary when running with 
!  barotropic_pred_corr=.true. and pred_corr_gamma > 0.0, since 
!  predictor-corrector has dissipation from pred_corr_gamma > 0.0.
!  </DATA> 
!  <DATA NAME="smooth_eta_t_fs_biharmonic" TYPE="logical">
!  For spatially smoothing the eta_t field at each barotropic time step 
!  with a biharmonic operator. May not be necessary when running with 
!  barotropic_pred_corr=.true. and pred_corr_gamma > 0.0, since 
!  predictor-corrector has dissipation from pred_corr_gamma > 0.0. 
!  </DATA> 
!  <DATA NAME="smooth_eta_t_laplacian" TYPE="logical">
!  For spatially smoothing the eta_t field on the big time step
!  by using a laplacian operator. For compatibility 
!  and global conservation, must also introduce a mixing 
!  to the thickness weighted tracer concentration in the k=1 cell. 
!  </DATA> 
!  <DATA NAME="smooth_eta_t_biharmonic" TYPE="logical">
!  For spatially smoothing the eta_t field on the big time step
!  by using a biharmonic operator. For compatibility 
!  and global conservation, must also introduce a mixing 
!  to the thickness weighted tracer concentration in the k=1 cell. 
!  </DATA> 
!  <DATA NAME="vel_micom_lap" UNITS="m/sec" TYPE="real">
!  Velocity scale that is used for computing the MICOM Laplacian mixing 
!  coefficient used in the Laplacian smoothing of surface height.
!  </DATA>
!  <DATA NAME="vel_micom_bih" UNITS="m/sec" TYPE="real">
!  Velocity scale that is used for computing the MICOM biharmonic mixing 
!  coefficient used in the Laplacian smoothing of surface height.
!  </DATA>
!
!  <DATA NAME="tidal_forcing_m2" TYPE="logical">
!  Forces from lunar M2 tidal constituent. 
!  </DATA> 
!  <DATA NAME="tidal_forcing_8" TYPE="logical">
!  Forces from 8 lunar and solar tidal constituents. 
!  </DATA> 
!
!  <DATA NAME="truncate_eta" TYPE="logical">
!  To truncate the surface height deviation so to ensure positive thickness 
!  within the top cell. This method will not conserve volume or tracer. 
!  It is coded for cases when conservation is not critical but wish to 
!  run with large free surface height deviations, such as when running with 
!  tides and very fine vertical resolution. 
!  </DATA> 
!  <DATA NAME="verbose_truncate" TYPE="logical">
!  For verbose printout on truncate_eta
!  </DATA> 
!  <DATA NAME="frac_crit_cell_height" UNITS="dimensionless" TYPE="real">
!  The top model tracer grid cell has thickness dht(i,j,1) = dzt(1) + eta_t(i,j).
!  0 < frac_crit_cell_height <= 1 sets the fraction of dzt(1) that is allowed
!  prior to bringing the model down due to overly small dht(i,j,1).
!  </DATA> 
!  <DATA NAME="eta_max" UNITS="meter" TYPE="real">
!  The maximum positive eta_t allowed when truncate_eta is true. 
!  </DATA> 
!
!  <DATA NAME="debug_fs" TYPE="logical">
!  For debugging.
!  </DATA> 
!  <DATA NAME="verbose_init" TYPE="logical">
!  For brief or full printout on initialization
!  </DATA> 
!  <DATA NAME="diag_freq" TYPE="integer">
!  Frequency for output of freesurf diagnostics. 
!  </DATA> 
!
!</NAMELIST>

use constants_mod,    only: epsln, grav, radian, rho0
use diag_manager_mod, only: register_diag_field, register_static_field, send_data, need_data
use fms_mod,          only: write_version_number, write_data, read_data, FATAL, WARNING, NOTE
use fms_mod,          only: open_namelist_file, check_nml_error, close_file, file_exist
use fms_io_mod,       only: field_size
use mpp_domains_mod,  only: BGRID_NE, CGRID_NE, BITWISE_EXACT_SUM
use mpp_domains_mod,  only: mpp_update_domains, mpp_global_field, mpp_global_sum, NUPDATE, EUPDATE
use mpp_domains_mod,  only: mpp_get_domain_components,mpp_get_layout, domain1d,mpp_get_pelist
use mpp_mod,          only: mpp_chksum, mpp_max, mpp_min, mpp_sum, mpp_root_pe, mpp_pe
use mpp_mod,          only: mpp_error, mpp_broadcast, stdout, stdlog
use mpp_mod,          only: mpp_send, mpp_recv, mpp_sync_self 
use time_manager_mod, only: time_type, increment_time, get_time, get_date
 
use ocean_coriolis_mod,        only: coriolis
use ocean_domains_mod,         only: get_local_indices, get_global_indices, set_ocean_domain
use ocean_obc_mod,             only: ocean_obc_update_boundary, ocean_obc_freesurf
use ocean_obc_mod,             only: ocean_obc_adjust_forcing_fs
use ocean_operators_mod,       only: FDX_PT_flat, FDY_PT_flat, BDX_ET, BDY_NT, FMX, FMY
use ocean_operators_mod,       only: DIV_UD, REMAP_BT_TO_BU, GRAD_SURF_P, LAP_T
use ocean_types_mod,           only: ocean_domain_type, ocean_grid_type, ocean_thickness_type
use ocean_types_mod,           only: ocean_time_type, ocean_time_steps_type
use ocean_types_mod,           only: ocean_external_mode_type
use ocean_types_mod,           only: ocean_velocity_type, ocean_density_type
use ocean_types_mod,           only: ocean_adv_vel_type, ocean_prog_tracer_type
use ocean_types_mod,           only: TWO_LEVEL, THREE_LEVEL
use ocean_util_mod,            only: write_timestamp

implicit none

public ocean_freesurf_init
public update_ocean_surface_height
public update_ocean_barotropic
public ocean_barotropic_forcing
public surface_height_tendency
public ocean_surface_smooth
public ocean_freesurf_rstrt
public ocean_freesurf_end

private leap_frog_barotropic
private pred_corr_barotropic
private tidal_forcing_init
private get_tidal_forcing
private read_freesurf

private eta_check 
private freesurf_energy
private freesurf_integrals
private ocean_eta_chksum
private ocean_ud_chksum
private psi_compute
private eta_truncate
private maximum_convud

private

logical  :: zero_tendency=.false.     ! If true, will freeze the external mode fields 
logical  :: zero_eta_ic=.false.       ! If true, will initialize eta to zero.

logical  :: truncate_eta = .false.    ! if true, will keep eta_t small enough to ensure dzt(1) + eta_t > 0.
logical  :: verbose_truncate = .true. ! for full printout of the truncate_eta points. 
real     :: frac_crit_cell_height=.20 ! fraction of dzt(1) that is deemed too thin for upper level dht(i,j,1)
real     :: eta_max = 5.0             ! max amplitude surface height fluctuation (meter)

logical  :: debug_fs=.false.          ! for debugging--prints out lots of checksums 
logical  :: verbose_init=.true.       ! for printouts during initialization  
integer  :: diag_freq=-1              ! for printing ascii diagnostics 

! for tidal forcing 
logical  :: tidal_forcing_m2=.false.  ! for tidal forcing by just the M2 constituent 
logical  :: tidal_forcing_8=.false.   ! for tidal forcing by 8 lunar/solar constituents 

! for smoothing surface height either on each small barotropic time steps or the big-time step
logical  :: smooth_eta_t_fs_laplacian=.false.  ! to smooth eta_t_fs with Laplacian operator 
logical  :: smooth_eta_t_fs_biharmonic=.false. ! to smooth eta_t_fs with biharmonic operator 
logical  :: smooth_eta_t_laplacian=.false.     ! to smooth eta_t with Laplacian operator 
logical  :: smooth_eta_t_biharmonic=.true.     ! to smooth eta_t with biharmonic operator 
real     :: vel_micom_lap=.05                  ! velocity (m/s) to set Laplacian mixing to smooth eta.
real     :: vel_micom_bih=.01                  ! velocity (m/s) to set biharmonic mixing to smooth eta.

! barotropic leap-frog specific 
logical  :: barotropic_leap_frog=.false. ! if wish to use leap-frog time scheme
logical  :: splitting=.true.             ! internally set false if dtuv=dtfs
real     :: robert_asselin=0.0           ! dimensionless parameter to damp leap-frog splitting mode if dtfs=dtuv  
real     :: twodt=0.0                    ! internally set to 2*dtfs
 
! barotropic predictor-corrector specific 
logical  :: barotropic_pred_corr=.false. ! if wish to use predictor-corrector scheme 
real     :: pred_corr_gamma=0.2          ! dissipation parameter for predictor-corrector 

! internally set 
logical  :: module_is_initialized=.false. ! to indicate that the module has been properly initialized 
logical  :: have_obc=.false.              ! for running with OBC

integer  :: tendency      ! for discretization of time tendency ("threelevel" or "twolevel")
integer  :: nts=1         ! internally set to the number of free surface timesteps per dtuv timestep
real     :: dtts          ! tracer time step (secs)
real     :: dtuv          ! baroclinic velocity time step (secs)
real     :: dtfs          ! barotropic time step (secs)
real     :: dteta         ! timestep (secs) determining update of surface height
real     :: dtime         ! timestep (secs) (dtime=2*dteta for threelevel and dtime=dteta for twolevel)
real     :: dtimer        ! 1/dtime
real     :: area_total_r  ! inverse area (1/m) of k=1 tracer cells 

! for psi_compute
integer, allocatable :: pelist_x(:)
integer, allocatable :: pelist_y(:)
integer  :: layout(2)
integer  :: myid_x
integer  :: myid_y
integer  :: size_x
integer  :: size_y
integer  :: this_pe
type(domain1d),save :: Domx
type(domain1d),save :: Domy


character(len=128) :: &
     version='$Id$'
character (len=128) :: tagname = &
     '$Name$'

#include <ocean_memory.h>

#ifdef STATIC_MEMORY

real, dimension(isd:ied,jsd:jed,3)    :: eta_t_fs   ! sea surface height(m) on T-cells on free surface timesteps
real, dimension(isd:ied,jsd:jed,2)    :: grad_ps    ! gradient of surface pressure 
real, dimension(isd:ied,jsd:jed,2,3)  :: ud_fs      ! vertical average of depth*velocity on U-cells

real, dimension(isd:ied,jsd:jed)      :: eta_mix_lap ! mixing coefficient (m^2/s) for laplacian smoothing of eta_t
real, dimension(isd:ied,jsd:jed)      :: eta_mix_bih ! mixing coefficient (m^4/s) for biharmonic smoothing of eta_t
real, dimension(isd:ied,jsd:jed)      :: rho_g       ! rho(:,:,k=1)*grav
real, dimension(isd:ied,jsd:jed)      :: area_t      ! t-cell areas for k=1 including tmask 
real, dimension(isd:ied,jsd:jed)      :: tmp         ! temporary array for various calculations 

#else

real, dimension(:,:,:), allocatable   :: eta_t_fs    ! sea surface height(m) on T-cells on free surface timesteps
real, dimension(:,:,:), allocatable   :: grad_ps     ! gradient of surface pressure 
real, dimension(:,:,:,:), allocatable :: ud_fs       ! depth averaged velocity * ocean depth on U-cells

real, dimension(:,:), allocatable     :: eta_mix_lap ! mixing coefficient (m^2/s) for laplacian smoothing of eta_t
real, dimension(:,:), allocatable     :: eta_mix_bih ! mixing coefficient (m^4/s) for biharmonic smoothing of eta_t
real, dimension(:,:), allocatable     :: rho_g       ! rho(:,:,k=1)*grav 
real, dimension(:,:), allocatable     :: area_t      ! t-cell areas for k=1 including tmask 
real, dimension(:,:), allocatable     :: tmp         ! temporary array for various calculations 

#endif

! for diagnostics manager 
logical :: used
integer :: id_eta_t=-1
integer :: id_eta_t_fs=-1
integer :: id_eta_mix_lap=-1
integer :: id_eta_mix_bih=-1
integer :: id_deta_dt=-1
integer :: id_eta_u=-1
integer :: id_ud=-1
integer :: id_vd=-1
integer :: id_psiu=-1
integer :: id_psiv=-1
integer :: id_ps=-1
integer :: id_psx=-1
integer :: id_psy=-1
integer :: id_convud_t=-1
integer :: id_eta_source=-1
integer :: id_surface_smooth=-1
integer :: id_eta_eq_tidal =-1
integer :: id_ke_fs=-1
integer :: id_pe_fs=-1
integer :: id_pmei=-1
integer :: id_riveri=-1
integer :: id_etati=-1

! for ascii output
integer :: unit=6

type(ocean_domain_type), pointer :: Dom =>NULL()
type(ocean_grid_type), pointer   :: Grd =>NULL()
type(ocean_domain_type), save    :: Dom_flux


! for tidal forcing 
real    :: alphat                         ! ocean loading and self-attraction from tidal forcing (=0.948)
real    :: rradian                        ! inverse radians
real    :: tidal_omega_M2,Love_M2,amp_M2
real    :: tidal_omega_S2,Love_S2,amp_S2
real    :: tidal_omega_N2,Love_N2,amp_N2
real    :: tidal_omega_K2,Love_K2,amp_K2
real    :: tidal_omega_K1,Love_K1,amp_K1
real    :: tidal_omega_O1,Love_O1,amp_O1
real    :: tidal_omega_P1,Love_P1,amp_P1
real    :: tidal_omega_Q1,Love_Q1,amp_Q1
 
#ifdef STATIC_MEMORY
real, private, dimension(isd:ied,jsd:jed)     :: eta_eq_tidal! equilibrium tidal forcing (m)
real, private, dimension(isd:ied,jsd:jed)     :: cos2lon     ! cosine   of 2*longitude on T-cells
real, private, dimension(isd:ied,jsd:jed)     :: coslon      ! cosine   of   longitude on T-cells
real, private, dimension(isd:ied,jsd:jed)     :: coslat2     ! cosine^2 of   longitude on T-cells
real, private, dimension(isd:ied,jsd:jed)     :: sin2lon     !   sine   of 2*longitude on T-cells
real, private, dimension(isd:ied,jsd:jed)     :: sinlon      !   sine   of   longitude on T-cells
real, private, dimension(isd:ied,jsd:jed)     :: sin2lat     !   sine   of 2*latitude  on T-cells
#else
real, private, dimension(:,:), allocatable    :: eta_eq_tidal! equilibrium tidal forcing (m)
real, private, dimension(:,:), allocatable    :: cos2lon     ! cosine   of 2*longitude on T-cells
real, private, dimension(:,:), allocatable    :: coslon      ! cosine   of   longitude on T-cells
real, private, dimension(:,:), allocatable    :: coslat2     ! cosine^2 of   longitude on T-cells
real, private, dimension(:,:), allocatable    :: sin2lon     !   sine   of 2*longitude on T-cells
real, private, dimension(:,:), allocatable    :: sinlon      !   sine   of   longitude on T-cells
real, private, dimension(:,:), allocatable    :: sin2lat     !   sine   of 2*latitude  on T-cells
#endif

namelist /ocean_freesurf_nml/ zero_tendency, zero_eta_ic,                                     &
                              tidal_forcing_m2, tidal_forcing_8,                              &
                              barotropic_pred_corr, pred_corr_gamma,                          &
                              barotropic_leap_frog, robert_asselin,                           &
                              smooth_eta_t_fs_laplacian, smooth_eta_t_fs_biharmonic,          &
                              smooth_eta_t_laplacian, smooth_eta_t_biharmonic,                &
                              vel_micom_lap, vel_micom_bih,                                   & 
                              truncate_eta, verbose_truncate, eta_max, frac_crit_cell_height, &
                              verbose_init, debug_fs, diag_freq
                                
contains


!#######################################################################
! <SUBROUTINE NAME="ocean_freesurf_init">
!
! <DESCRIPTION>
! Initialize the free surface module
! </DESCRIPTION>
!
subroutine ocean_freesurf_init(Grid, Domain, Time, Time_steps, Ext_mode, obc, debug, ens_ocean)

  type(ocean_domain_type), intent(in),   target :: Domain
  type(ocean_grid_type), intent(in),     target :: Grid
  type(ocean_time_type), intent(in)             :: Time
  type(ocean_time_steps_type), intent(in)       :: Time_steps 
  type(ocean_external_mode_type), intent(inout) :: Ext_mode
  logical,                           intent(in) :: obc
  logical,                 intent(in), optional :: debug
  logical,                 intent(in), optional :: ens_ocean

  integer :: tau, taup1
  integer :: ioun, io_status
  real    :: crit, fudge
  integer :: n
  real    :: cfl_grid_factor
  real    :: cgmax, gridmin, cgext, gridsp, dtcg
  real    :: max_dt_for_cgext, max_dt_for_cgext0
  integer :: i, j, k, icg, jcg, ierr
  logical :: cfl_error=.false.

  if ( module_is_initialized ) then 
    call mpp_error(FATAL, '==>Error in ocean_freesurf_mod (ocean_freesurf_init): module already initialized.')
  endif 

  module_is_initialized = .TRUE.

  if (diag_freq == 0) diag_freq = 1

  call write_version_number( version, tagname )

  ioun = open_namelist_file()
  read  (ioun, ocean_freesurf_nml,iostat=io_status)
  write (stdout(),'(/)')
  write (stdout(), ocean_freesurf_nml)  
  write (stdlog(), ocean_freesurf_nml)
  ierr = check_nml_error(io_status, 'ocean_freesurf_nml')
  call close_file(ioun)

  if (PRESENT(debug) .and. .not. debug_fs) then
    debug_fs = debug
  endif 

  have_obc = obc
  
  Dom => Domain
  Grd => Grid
  call set_ocean_domain(Dom_flux, Grid, name='horz diff flux')


#ifndef STATIC_MEMORY

  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  call get_global_indices(Domain, isg, ieg, jsg, jeg)
  nk = Grid%nk

  allocate (Ext_mode%eta_t(isd:ied,jsd:jed,3))
  allocate (Ext_mode%eta_u(isd:ied,jsd:jed,3))
  allocate (Ext_mode%eta_t_bar(isd:ied,jsd:jed,3))
  allocate (Ext_mode%deta_dt(isd:ied,jsd:jed)) 
  allocate (Ext_mode%convud_t(isd:ied,jsd:jed,3))
  allocate (Ext_mode%ud(isd:ied,jsd:jed,2,3))
  allocate (Ext_mode%forcing_fs(isd:ied,jsd:jed,2))
  allocate (Ext_mode%ps(isd:ied,jsd:jed))
  allocate (Ext_mode%grad_ps(isd:ied,jsd:jed,2))
  allocate (Ext_mode%eta_source(isd:ied,jsd:jed))
  allocate (Ext_mode%surface_smooth(isd:ied,jsd:jed))
   
  allocate (tmp(isd:ied, jsd:jed))
  allocate (area_t(isd:ied, jsd:jed))
  allocate (eta_t_fs(isd:ied,jsd:jed,3))
  allocate (ud_fs(isd:ied,jsd:jed,2,3))
  allocate (rho_g(isd:ied,jsd:jed))
  
#endif
  
  tendency = Time_steps%tendency
  dtts     = Time_steps%dtts
  dtuv     = Time_steps%dtuv
  dtfs     = Time_steps%dtfs
  dteta    = Time_steps%dteta
  dtime    = Time_steps%dtime_e
  dtimer   = 1.0/(dtime+epsln)
  
  Ext_mode%eta_t(:,:,:)        = 0.0
  Ext_mode%eta_u(:,:,:)        = 0.0
  Ext_mode%convud_t(:,:,:)     = 0.0
  Ext_mode%ud(:,:,:,:)         = 0.0
  Ext_mode%forcing_fs(:,:,:)   = 0.0
  Ext_mode%ps(:,:)             = 0.0
  Ext_mode%grad_ps(:,:,:)      = 0.0
  Ext_mode%deta_dt(:,:)        = 0.0
  Ext_mode%eta_source(:,:)     = 0.0
  Ext_mode%surface_smooth(:,:) = 0.0
  Ext_mode%eta_t_bar(:,:,:)    = 0.0

  eta_t_fs(:,:,:)            = 0.0
  ud_fs(:,:,:,:)             = 0.0
  rho_g(:,:)                 = 0.0
  tmp(:,:)                   = 0.0
  area_t(:,:)                = Grd%dat(:,:)*Grd%tmask(:,:,1)
  area_total_r               = 1.0/(epsln+Grd%tcella(1)) 

#ifndef STATIC_MEMORY    
  allocate (eta_mix_lap(isd:ied,jsd:jed))
  allocate (eta_mix_bih(isd:ied,jsd:jed))
#endif    
  eta_mix_lap(:,:) = Grd%tmask(:,:,1)*vel_micom_lap*(2.0*Grd%dxt(:,:)*Grd%dyt(:,:)/(Grd%dxt(:,:)+Grd%dyt(:,:)))
  eta_mix_bih(:,:) = Grd%tmask(:,:,1)*vel_micom_bih*(2.0*Grd%dxt(:,:)*Grd%dyt(:,:)/(Grd%dxt(:,:)+Grd%dyt(:,:)))**3

  ! static fields for diagnostic manager 
  id_eta_mix_lap = register_static_field ('ocean_model', 'eta_mix_lap', Grd%tracer_axes(1:2), &
                                'laplacian mixing coefficient for eta smoothing', 'm^2/s',&
                                missing_value=0.0, range=(/0.0,1e12/))
  if (id_eta_mix_lap > 0) used = send_data(id_eta_mix_lap, eta_mix_lap(isc:iec,jsc:jec), Time%model_time)

  id_eta_mix_bih = register_static_field ('ocean_model', 'eta_mix_bih', Grd%tracer_axes(1:2), &
                                'biharmonic mixing coefficient for eta smoothing', 'm^4/s',&
                                missing_value=0.0, range=(/0.0,1e12/))
  if (id_eta_mix_bih > 0) used = send_data(id_eta_mix_bih, eta_mix_bih(isc:iec,jsc:jec), Time%model_time)


  ! convert eta_mix_bih to its sqrt form so to apply it to both parts of the iterated Laplacian 
  eta_mix_bih(:,:) = sqrt(eta_mix_bih(:,:))


  ! register time dependent fields for diagnostic manager 

  id_eta_t  = register_diag_field ('ocean_model', 'eta_t', Grd%tracer_axes(1:2), Time%model_time,&
                                'surface height on T cells', 'meter',&
                                missing_value=-10.0, range=(/-10.0,10.0/))

  id_eta_t_fs  = register_diag_field ('ocean_model', 'eta_t_fs', Grd%tracer_axes(1:2), Time%model_time,&
                                'surface height on T cells w/i barotropic loop', 'meter',&
                                missing_value=-10.0, range=(/-10.0,10.0/))

  id_deta_dt  = register_diag_field ('ocean_model', 'deta_dt', Grd%tracer_axes(1:2), Time%model_time,&
                                'tendency for eta_t', 'm/s', missing_value=-1e6, range=(/-1e6,1e6/))

  id_eta_u  = register_diag_field ('ocean_model', 'eta_u', Grd%vel_axes_uv(1:2), Time%model_time, &
                                'surface height on U cells', 'meter', missing_value=-10.0, range=(/-10.0,10.0/))

  id_ps   = register_diag_field ('ocean_model', 'ps', Grd%tracer_axes(1:2), Time%model_time, &
                                'surf+atmos pressure', 'Pascal',missing_value=-1e36, range=(/-1e5,1e5/))

  id_psx  = register_diag_field ('ocean_model', 'psx', Grd%vel_axes_uv(1:2), Time%model_time, &
                                'zonal surf pressure force', 'm/s^2', missing_value=-1e36, range=(/-1e5,1e5/))

  id_psy  = register_diag_field ('ocean_model', 'psy', Grd%vel_axes_uv(1:2), Time%model_time, &
                                'meridional surf pressure force', 'm/s^2', missing_value=-1e36, range=(/-1e5,1e5/))

  id_ud     = register_diag_field ('ocean_model', 'ud', Grd%vel_axes_uv(1:2), Time%model_time, &
                                'depth weighted u', 'meters*meters/sec', missing_value=-10.0, range=(/-10.0,10.0/))

  id_vd     = register_diag_field ('ocean_model', 'vd', Grd%vel_axes_uv(1:2), Time%model_time, &
                                'depth weighted v', 'meters*meters/sec', missing_value=-10.0, range=(/-10.0,10.0/))

  id_psiu   = register_diag_field ('ocean_model', 'psiu', Grd%vel_axes_uv(1:2), Time%model_time, &
                                'quasi-barotropic strmfcn psiu', 'Sv', missing_value=-1e5, range=(/-1e5,1e5/))

  id_psiv   = register_diag_field ('ocean_model', 'psiv', Grd%vel_axes_uv(1:2), Time%model_time, &
                                'quasi-barotropic strmfcn psiv', 'Sv', missing_value=-1e5, range=(/-1e5,1e5/))

  id_convud_t = register_diag_field ('ocean_model', 'convud_t', Grd%tracer_axes(1:2), Time%model_time, &
                                'convergence ud on T cells', 'm/sec', missing_value=-10.0, range=(/-10.0,10.0/))

  id_eta_source  = register_diag_field ('ocean_model', 'eta_source', Grd%tracer_axes(1:2), Time%model_time, &
                 'total eta source', 'm/sec', missing_value=-1e6, range=(/-1e6,1e6/))  

  id_surface_smooth  = register_diag_field ('ocean_model', 'eta_surf_smooth', Grd%tracer_axes(1:2), Time%model_time, &
                 'surface eta_t smoother', 'm/s', missing_value=-1e6, range=(/-1e6,1e6/))  

  id_ke_fs = register_diag_field ('ocean_model', 'ke_fs', Time%model_time, &
             'Free surface ke', '10^-15 joules', missing_value=0.0, range=(/0.0,1000.0/))

  id_pe_fs = register_diag_field ('ocean_model', 'pe_fs', Time%model_time, &
             'Free surface pe', '10^-15 joules', missing_value=0.0, range=(/0.0,1000.0/))

  id_pmei  = register_diag_field ('ocean_model', 'pmei', Time%model_time, &
             'pme integral', 'm/sec', missing_value=0.0, range=(/0.0,1000.0/))

  id_riveri  = register_diag_field ('ocean_model', 'riveri', Time%model_time, &
             'river integral', 'm/sec', missing_value=0.0, range=(/0.0,1000.0/))

  id_etati = register_diag_field ('ocean_model', 'etati', Time%model_time, &
             'eta_t integral','m/sec', missing_value=0.0, range=(/0.0,1000.0/))

  if (.NOT. file_exist('INPUT/ocean_freesurf.res.nc')) then
      if (.NOT.Time%init) call mpp_error(FATAL,'Expecting file INPUT/ocean_freesurf.res.nc to exist.&
           &This file was not found and Time%init=.false.')

      if( file_exist('INPUT/eta_t_ic.nc') ) then
          call read_data('INPUT/eta_t_ic.nc', 'eta_t',Ext_mode%eta_t(:,:,1), &
               Dom%domain2d,timelevel=1,  append_pelist_name=ens_ocean)
          call mpp_update_domains (Ext_mode%eta_t(:,:,1), Dom%domain2d)
          Ext_mode%eta_t(:,:,2)  = Ext_mode%eta_t(:,:,1)
          Ext_mode%eta_t(:,:,3)  = Ext_mode%eta_t(:,:,1)
          eta_t_fs(:,:,1)        = Ext_mode%eta_t(:,:,1)
          eta_t_fs(:,:,2)        = Ext_mode%eta_t(:,:,1)
          eta_t_fs(:,:,3)        = Ext_mode%eta_t(:,:,1)
          Ext_mode%eta_t_bar(:,:,1) = Ext_mode%eta_t(:,:,1)
          Ext_mode%eta_t_bar(:,:,2) = Ext_mode%eta_t(:,:,1)
          Ext_mode%eta_t_bar(:,:,3) = Ext_mode%eta_t(:,:,1)
      endif
  else
      call read_freesurf(Time, Ext_mode, ens_ocean)
  endif

  if(Time%init .and. zero_eta_ic) then 
      call mpp_error(NOTE,'==>Setting initial freesurface fields to zero. This overrides any initial or restart files.')
      Ext_mode%eta_t(:,:,:)     = 0.0
      Ext_mode%convud_t(:,:,:)  = 0.0
      Ext_mode%eta_t_bar(:,:,:) = 0.0
      Ext_mode%eta_u(:,:,:)     = 0.0
      Ext_mode%ud(:,:,:,:)      = 0.0
  endif

  if (zero_tendency) then
     call mpp_error(NOTE,'Using zero_tendency in ocean_freesurf_mod.  Will freeze external mode fields')
  endif 

  if (barotropic_leap_frog) then
     call mpp_error(NOTE,'Using barotropic_leap_frog=.true. for integrating barotropic dynamics.')
     write (stdout(),'(/a,f6.2)')' Robert-Asselin filter on barotropic dynamics has value= ',robert_asselin
  endif 
  if (barotropic_pred_corr) then
     call mpp_error(NOTE,'Using barotropic_pred_corr=.true. for integrating barotropic dynamics.')
     write (stdout(),'(/a,f6.2)')' Predictor-Corrector time filter on barotropic dynamics has value= ',pred_corr_gamma
     if(have_obc) then 
        call mpp_error(WARNING,'OBC has not been generalized to barotropic_pred_corr=.true.')
     endif 
  endif 
  if (barotropic_pred_corr .and. barotropic_leap_frog) then
     call mpp_error(FATAL,'==>Error from ocean_freesurf_mod: barotropic_pred_corr and barotropic_leap_frog cannot both be true.')
  endif 
  if (.not. barotropic_pred_corr .and. .not. barotropic_leap_frog) then
     call mpp_error(FATAL,'==>Error from ocean_freesurf_mod: barotropic_pred_corr and barotropic_leap_frog cannot both be false.')
  endif 

  ! timestep consistancy checks
  if (nint(dtuv) > 0) then

    if (nint(dtfs) > nint(dtuv)) then
      call mpp_error(FATAL,'==>Error from ocean_freesurf_mod: dtfs > dtuv is not allowed.')

    elseif (nint(dtfs) == 0) then
      call mpp_error(FATAL,'==>Error from ocean_freesurf_mod: dtfs=0 when dtuv > 0 is not implemented.')

    elseif (nint(dtuv) == nint(dtfs)) then

      if (barotropic_pred_corr) then
        call mpp_error(FATAL,'==>Error from ocean_freesurf_mod: dtfs=dtuv not allowed with barotropic_pred_corr')
      else
        call mpp_error(NOTE,'==>From ocean_freesurf_mod: dtfs=dtuv not fully supported. User beware.')
      endif 

      splitting = .false.

      nts = 1
      twodt = (2.0*dtuv)
      write (stdout(),'(/a)')' Warning from ocean_freesurf_mod: un-split momentum equations &
                               &imply one free surface timestep per velocity timestep.'
    else

      splitting = .true.

      nts  = nint((2.0*dtuv)/dtfs)
      crit = 1.e-6
      if (abs((2.0*dtuv)/nts - dtfs) > crit) then
        write (unit,'(/1x,a/)') ' ==>Error: Splitting the momentum equations requires &
             &that mod(2*dtuv,dtfs)=0. It is not!'
        write (unit,*)      '          Current setting: '
        write (unit,*)      '             dtuv =  ',dtuv,' sec.'
        write (unit,*)      '             dtfs =  ',dtfs,' sec.'
        write (unit,*)      '          New setting should be:  '
        write (unit,*)      '             dtfs =  ',(2.0*dtuv)/nts,' sec.'
        call mpp_error(FATAL,'==>Error from ocean_freesurf_mod: badly chosen time steps')
      endif

      if(barotropic_leap_frog) twodt = 2.0*(2.0*dtuv)/nts

      write (stdout(),'(/1x,a,i4,a)')' ==> Note: The barotropic dynamics integrate ',&
           nts, ' timesteps for every one baroclinic timestep.'
    endif
  endif

  ! Determine maximum timestep allowable to ensure barotropic 
  ! time stepping is CFL stable with external gravity waves

  cgmax            = 1.0
  icg              = isc
  jcg              = jsc
  gridmin          = 1.0e20
  max_dt_for_cgext = 1.e6

  if(barotropic_leap_frog) cfl_grid_factor=0.5
  if(barotropic_pred_corr) cfl_grid_factor=1.0
  do j=jsc,jec
     do i=isc,iec
        if (Grd%kmu(i,j) > 0) then
            cgext  = sqrt(grav*(Grd%zw(Grd%kmu(i,j))))
            gridsp = 1.0/(Grd%dxu(i,j)*Grd%dxu(i,j)) + 1.0/(Grd%dyu(i,j)*Grd%dyu(i,j))
            gridsp = sqrt(1.0/gridsp) 
            dtcg   = cfl_grid_factor*gridsp/cgext
            if (dtcg < max_dt_for_cgext) then
                gridmin = gridsp; cgmax = cgext
                max_dt_for_cgext = dtcg; icg = i; jcg = j
            endif
        endif
     enddo
  enddo

  if (nint(dtfs) > 0) then

    cgmax = nint(cgmax)
    max_dt_for_cgext = int(max_dt_for_cgext)

    fudge = 1 + 1.e-12*mpp_pe()  ! to distinguish redundancies between processors
    max_dt_for_cgext  = max_dt_for_cgext*fudge
    max_dt_for_cgext0 = max_dt_for_cgext
    call mpp_min (max_dt_for_cgext)
    call mpp_max (cgmax)

    ! show the most unstable location for barotropic gravity waves
    if (max_dt_for_cgext == max_dt_for_cgext0) then

      if (dtfs <= max_dt_for_cgext) then
        write (unit,'(/a,i4,a,i4,a,f6.2,a,f6.2,a)')' Note: Barotropic stability most nearly violated at U-cell (i,j) = (',&
        icg+Dom%ioff,',',jcg+Dom%joff,'), (lon,lat) = (',Grd%xu(icg,jcg),',',Grd%yu(icg,jcg),').'
        write(unit,'(a,i6)')    '         The number of kmu-levels at this point is ',Grd%kmu(icg,jcg) 
        write(unit,'(a,e12.6)') '         The dxu grid spacing at this point is ',Grd%dxu(icg,jcg) 
        write(unit,'(a,e12.6)') '         The dyu grid spacing at this point is ',Grd%dyu(icg,jcg)
        cfl_error=.false.
      else
        write (unit,'(/a,i4,a,i4,a,f6.2,a,f6.2,a)')'=>Error: Barotropic stability violated at U-cell (i,j) = (',&
        icg+Dom%ioff,',',jcg+Dom%joff,'), (lon,lat) = (',Grd%xu(icg,jcg),',',Grd%yu(icg,jcg),').'
        write(unit,'(a,i6)')    '         The number of kmu-levels at this point is ',Grd%kmu(icg,jcg) 
        write(unit,'(a,e12.6)') '         The dxu grid spacing at this point is ',Grd%dxu(icg,jcg) 
        write(unit,'(a,e12.6)') '         The dyu grid spacing at this point is ',Grd%dyu(icg,jcg) 
        cfl_error=.true.  
      endif

      write (unit,'(a,f5.1,a/a,f8.3,a,f8.3,a)')&
      '         where the barotropic gravity wave speed is ~',cgmax,' m/s.',&
      '         "dtfs" must be less than ',max_dt_for_cgext/fudge,' sec.   dtfs = ',dtfs,' sec.'

    endif

    if(cfl_error) then 
        call mpp_error(FATAL,'==>Error from ocean_freesurf_mod: time step instability &
                              &detected for external gravity waves. Time step too large.' )
    endif 

    max_dt_for_cgext = nint(max_dt_for_cgext)

    if (dtts /= dtuv .and. .not. splitting) then
      call mpp_error(FATAL,'==>Error from ocean_freesurf_mod: Solving the asynchronous &
                            & unsplit equations has not been implemented.')
    endif

  endif

  if(smooth_eta_t_fs_laplacian) then 
    write(stdout(),'(/1x,a)') ' ==> Using smooth_eta_t_fs_laplacian to smooth eta_t_fs on each barotropic time step (dtfs).'
  endif 
  if(smooth_eta_t_fs_biharmonic) then 
    write(stdout(),'(/1x,a)') ' ==> Using smooth_eta_t_fs_biharmonic to smooth eta_t_fs on each barotropic time step (dtfs).'
  endif 
  if(smooth_eta_t_fs_laplacian .and. smooth_eta_t_fs_biharmonic) then 
    call mpp_error(FATAL, '==>ocean_freesurf_mod: smooth_eta_t_fs_laplacian & smooth_eta_t_fs_biharmonic cannot both be true.')
  endif 

  if(smooth_eta_t_laplacian) then 
    write(stdout(),'(/1x,a)') ' ==> Using smooth_eta_t_laplacian to smooth eta_t on surface height time step (dteta).'
  endif 
  if(smooth_eta_t_biharmonic) then 
    write(stdout(),'(/1x,a)') ' ==> Using smooth_eta_t_biharmonic to smooth eta_t on surface height time step (dteta).'
  endif 
  if(smooth_eta_t_laplacian .and. smooth_eta_t_biharmonic) then 
    call mpp_error(FATAL, '==>ocean_freesurf_mod: smooth_eta_t_laplacian & smooth_eta_t_biharmonic cannot both be true.')
  endif 

  if(truncate_eta) then 
    write(stdout(),'(/a)')'==>Warning: truncate_eta=.true.' 
    write(stdout(),'(a,f10.4,a,f10.4)') &
     '   Will enforce dzt(1)+eta_t (meters) > ',Grd%dzt(1)*frac_crit_cell_height, '  and eta_t (meters) < ', eta_max
  endif 

  call tidal_forcing_init(Time)

  if(debug_fs) then 
      tau   = Time%tau
      taup1 = Time%taup1
      write(stdout(),*) ' '
      write(stdout(),*) 'From ocean_freesurf_mod: initial external mode checksums'
      call write_timestamp(Time%model_time)
      write(stdout(),*) 'chksum for eta_t(tau)       = ',  mpp_chksum(Ext_mode%eta_t(isc:iec,jsc:jec,tau))
      write(stdout(),*) 'chksum for eta_t(taup1)     = ',  mpp_chksum(Ext_mode%eta_t(isc:iec,jsc:jec,taup1))
      write(stdout(),*) 'chksum for convud_t(tau)    = ',  mpp_chksum(Ext_mode%convud_t(isc:iec,jsc:jec,tau))
      write(stdout(),*) 'chksum for convud_t(taup1)  = ',  mpp_chksum(Ext_mode%convud_t(isc:iec,jsc:jec,taup1))
      write(stdout(),*) 'chksum for eta_t_bar(tau)   = ',  mpp_chksum(Ext_mode%eta_t_bar(isc:iec,jsc:jec,tau))
      write(stdout(),*) 'chksum for eta_t_bar(taup1) = ',  mpp_chksum(Ext_mode%eta_t_bar(isc:iec,jsc:jec,taup1))
      write(stdout(),*) 'chksum for eta_u(tau)       = ',  mpp_chksum(Ext_mode%eta_u(isc:iec,jsc:jec,tau))
      write(stdout(),*) 'chksum for eta_u(taup1)     = ',  mpp_chksum(Ext_mode%eta_u(isc:iec,jsc:jec,taup1))
      write(stdout(),*) 'chksum for ud(tau)          = ',  mpp_chksum(Ext_mode%ud(isc:iec,jsc:jec,1,tau))
      write(stdout(),*) 'chksum for ud(taup1)        = ',  mpp_chksum(Ext_mode%ud(isc:iec,jsc:jec,1,taup1))
      write(stdout(),*) 'chksum for vd(tau)          = ',  mpp_chksum(Ext_mode%ud(isc:iec,jsc:jec,2,tau))
      write(stdout(),*) 'chksum for vd(taup1)        = ',  mpp_chksum(Ext_mode%ud(isc:iec,jsc:jec,2,taup1))
  endif

  ! get pelist for psi_compute
  call mpp_get_layout (Dom%domain2d, layout)
  allocate ( pelist_x(layout(1)), pelist_y(layout(2)))
  call mpp_get_domain_components (Dom%domain2d, Domx, Domy)
  call mpp_get_pelist            ( Domx, pelist_x )
  call mpp_get_pelist            ( Domy, pelist_y )
  size_x = size(pelist_x(:))
  size_y = size(pelist_y(:))
  this_pe = mpp_pe()
  myid_y = 0    
  do i = 1,size_y
     if (this_pe == pelist_y(i)) then
        myid_y = i
        exit
     endif
  enddo
  if(myid_y == 0) call mpp_error(FATAL,'==>Error in ocean_freesurf psi_compute myid_y')
  myid_x = 0
  do i = 1,size_x
     if (this_pe == pelist_x(i)) then
        myid_x = i
        exit
     endif
  enddo
  if(myid_x == 0) call mpp_error(FATAL,'==>Error in ocean_freesurf psi_compute myid_x')


end subroutine ocean_freesurf_init
! </SUBROUTINE> NAME="ocean_freesurf_init"


!#######################################################################
! <SUBROUTINE NAME="update_ocean_surface_height">
!
! <DESCRIPTION>
! Time step the surface height using a "big" time step. 
! </DESCRIPTION>
!
subroutine update_ocean_surface_height (Time, Ext_mode, Dens, patm, pme, river)

  type(ocean_time_type), intent(in)              :: Time
  type(ocean_external_mode_type), intent(inout)  :: Ext_mode
  type(ocean_density_type), intent(in)           :: Dens
  real, intent(in), dimension(isd:ied,jsd:jed)   :: patm 
  real, dimension(isd:ied,jsd:jed), intent(in)   :: pme
  real, dimension(isd:ied,jsd:jed), intent(in)   :: river
  integer                                        :: taum1, tau, taup1, i, j

  taum1 = Time%taum1
  tau   = Time%tau
  taup1 = Time%taup1

  ! update fields to time taup1
  if (splitting) then

      ! update eta_t with a "big leap-frog" or "big forward"
      do j=jsc,jec
         do i=isc,iec
            Ext_mode%eta_t(i,j,taup1) = Ext_mode%eta_t(i,j,taum1) + dtime*Ext_mode%deta_dt(i,j)
         enddo
      enddo

      ! open boundary condition at big leap-frog time step
      ! Ext_mode%deta_dt should be correct, since the velocity gradient  
      ! across the boundary is zero
      ! the additional contribution from the radiation condition to
      ! Ext_mode%deta_dt is diagnosed also
      if(have_obc) call ocean_obc_freesurf(Time, Ext_mode)   

      if(truncate_eta) call eta_truncate(Time, Ext_mode)

      call mpp_update_domains (Ext_mode%eta_t(:,:,taup1), Dom%domain2d)

  endif

  ! send to diag_manager     
  if (id_eta_t > 0)    used = send_data (id_eta_t, Ext_mode%eta_t(:,:,tau), &
                              Time%model_time, rmask=Grd%tmask(:,:,1), &
                              is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

  ! find max and min of some quantities 
  if (diag_freq > 0 .and. mod(Time%itt, diag_freq) == 0) then 
     call eta_check(Time, Ext_mode)
     call maximum_convud(Time, Ext_mode)
  endif 

  ! compute some global integrals for diagnostics 
  call freesurf_integrals (Time, Ext_mode, pme, river)

  ! compute some checksums for debugging 
  if (debug_fs) then
     call ocean_eta_chksum(Time, Ext_mode, taup1)
  endif


end subroutine update_ocean_surface_height
! </SUBROUTINE> NAME="update_ocean_surface_height"


!#######################################################################
! <SUBROUTINE NAME="surface_height_tendency">
!
! <DESCRIPTION>
!  Compute tendency for surface height. 
! </DESCRIPTION>
!
subroutine surface_height_tendency(Time, Adv_vel, pme, river, Ext_mode)
  
  type(ocean_time_type), intent(in)             :: Time
  type(ocean_adv_vel_type), intent(in)          :: Adv_vel
  real, dimension(isd:ied,jsd:jed), intent(in)  :: pme
  real, dimension(isd:ied,jsd:jed), intent(in)  :: river
  type(ocean_external_mode_type), intent(inout) :: Ext_mode

  integer :: i, j, tau
  tau = Time%tau

  if (zero_tendency) then 
      Ext_mode%deta_dt(:,:) = 0.0
  else 
      do j=jsc,jec
         do i=isc,iec       
            Ext_mode%deta_dt(i,j) = Grd%tmask(i,j,1)                                    &
                                 *( Ext_mode%convud_t(i,j,tau) + pme(i,j) + river(i,j)  &
                                   +Ext_mode%eta_source(i,j)   + Ext_mode%surface_smooth(i,j) )
         enddo
      enddo
  endif

  if (id_eta_source > 0) used = send_data (id_eta_source, Ext_mode%eta_source(:,:), &
                                Time%model_time, rmask=Grd%tmask(:,:,1), &
                                is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)  

  if (id_deta_dt > 0)    used = send_data (id_deta_dt, Ext_mode%deta_dt(:,:), &
                                Time%model_time, rmask=Grd%tmask(:,:,1), &
                                is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

end subroutine surface_height_tendency
! </SUBROUTINE> NAME="surface_height_tendency"


!#######################################################################
! <SUBROUTINE NAME="ocean_surface_smooth">
!
! <DESCRIPTION>
!  Compute tendency for surface diffusion in both eta and tracer. 
!  Use either a laplacian or a biharmonic smoother. 
! </DESCRIPTION>
!
subroutine ocean_surface_smooth(Time, Thickness, T_prog, Ext_mode)

  type(ocean_time_type), intent(in)             :: Time
  type(ocean_thickness_type), intent(in)        :: Thickness
  type(ocean_prog_tracer_type), intent(inout)   :: T_prog(:)
  type(ocean_external_mode_type), intent(inout) :: Ext_mode

  integer :: i, j, n, nprog, taum1

  if(.not. smooth_eta_t_laplacian .and. .not. smooth_eta_t_biharmonic) return 

  taum1 = Time%taum1
  nprog = size(T_prog(:))  

  do j = jsd, jed
     do i = isd, ied
        tmp(i,j)   =  Ext_mode%eta_t(i,j,taum1)
     enddo
  enddo

  if(smooth_eta_t_laplacian) then 
    tmp(:,:) =  LAP_T(tmp(:,:),eta_mix_lap(:,:), update=.true.)
  else 
    tmp(:,:) = -LAP_T(tmp(:,:),eta_mix_bih(:,:), update=.true.)
    call mpp_update_domains (tmp, Dom%domain2d)
    tmp(:,:) =  LAP_T(tmp(:,:),eta_mix_bih(:,:), update=.true.)
  endif 
  do j=jsc,jec
     do i=isc,iec
        Ext_mode%surface_smooth(i,j) = Grd%tmask(i,j,1)*tmp(i,j)
     enddo
  enddo

  if(smooth_eta_t_laplacian) then 
      do n=1,nprog
         do j=jsd,jed
            do i=isd,ied
               tmp(i,j) =  Thickness%dht(i,j,1,taum1)*T_prog(n)%field(i,j,1,taum1)
            enddo
         enddo
         tmp(:,:) =  LAP_T(tmp(:,:),eta_mix_lap(:,:), update=.true.)
         do j=jsc,jec
            do i=isc,iec
               T_prog(n)%surface_smooth(i,j) = tmp(i,j)
            enddo
         enddo
      enddo
  else 
      do n=1,nprog
         do j=jsd,jed
            do i=isd,ied
               tmp(i,j) =  Thickness%dht(i,j,1,taum1)*T_prog(n)%field(i,j,1,taum1)
            enddo
         enddo
         tmp(:,:) = -LAP_T(tmp(:,:),eta_mix_bih(:,:), update=.true.)
         call mpp_update_domains (tmp, Dom%domain2d)
         tmp(:,:) =  LAP_T(tmp(:,:),eta_mix_bih(:,:), update=.true.)
         do j=jsc,jec
            do i=isc,iec
               T_prog(n)%surface_smooth(i,j) = tmp(i,j)
            enddo
         enddo
      enddo
  endif

  if (id_surface_smooth > 0) used = send_data (id_surface_smooth, Ext_mode%surface_smooth(:,:), &
                                    Time%model_time, rmask=Grd%tmask(:,:,1), &
                                    is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)  


end subroutine ocean_surface_smooth
! </SUBROUTINE> NAME="ocean_surface_smooth"



!#######################################################################
! <SUBROUTINE NAME="update_ocean_barotropic">
!
! <DESCRIPTION>
! Time step the external mode fields using either a leap-frog scheme
! or a predictor-corrector scheme.  Time average these fields to 
! update the vertically integrated velocity (ud,vd) and the time 
! averaged surface height eta_t_bar. 
! </DESCRIPTION>
!
subroutine update_ocean_barotropic (Time, Dens, Ext_mode, patm, pme, river)

  type(ocean_time_type), intent(in)              :: Time
  type(ocean_density_type), intent(in)           :: Dens
  type(ocean_external_mode_type), intent(inout)  :: Ext_mode
  real, intent(in), dimension(isd:ied,jsd:jed)   :: patm 
  real, intent(in), dimension(isd:ied,jsd:jed)   :: pme
  real, intent(in), dimension(isd:ied,jsd:jed)   :: river

  type(time_type)                                :: next_time 
  real, dimension(isd:ied,jsd:jed)               :: psiu
  real, dimension(isd:ied,jsd:jed)               :: psiv
  real                                           :: rnts
  integer                                        :: taum1, tau, taup1
  integer                                        :: i, j, n

  taum1 = Time%taum1
  tau   = Time%tau
  taup1 = Time%taup1

  if(tendency==THREE_LEVEL) then 
      do j=jsd,jed
         do i=isd,ied
            rho_g(i,j) = grav*Dens%rho(i,j,1)
         enddo
      enddo
  elseif(tendency==TWO_LEVEL) then 
      do j=jsd,jed
         do i=isd,ied
            rho_g(i,j) = grav*Dens%rho_taup1(i,j,1)
         enddo
      enddo
  endif

  if (nint(dtfs) == 0 .or. zero_tendency) then 
      do n = 1,2
         do j = jsd,jed
            do i = isd,ied
               Ext_mode%ud(i,j,n,taup1)     = Ext_mode%ud(i,j,n,tau)
            enddo
         enddo
      enddo
      do j = jsd,jed
         do i = isd,ied
            Ext_mode%eta_t(i,j,taup1)    = Ext_mode%eta_t(i,j,tau)
            Ext_mode%eta_u(i,j,taup1)    = Ext_mode%eta_u(i,j,tau)
            Ext_mode%convud_t(i,j,taup1) = Ext_mode%convud_t(i,j,tau)
            Ext_mode%ps(i,j)             = patm(i,j) + rho_g(i,j)*Ext_mode%eta_t(i,j,tau)
         enddo
      enddo
      Ext_mode%grad_ps(:,:,:)      = GRAD_SURF_P(Ext_mode%ps(:,:))
  else
      if(barotropic_leap_frog) call leap_frog_barotropic (Time, Ext_mode, patm, pme, river)
      if(barotropic_pred_corr) call pred_corr_barotropic (Time, Ext_mode, patm, pme, river)
  endif

  ! update fields to time taup1 if splitting 
  if (splitting) then

      rnts = 1.0/(float(nts)+1.0)      
      do j = jsd, jed
         do i = isd, ied
            Ext_mode%eta_t_bar(i,j,taup1) = rnts*Ext_mode%eta_t_bar(i,j,taup1) 
         enddo
      enddo
      call mpp_update_domains(Ext_mode%eta_t_bar(:,:,taup1), Dom%domain2d)

      Ext_mode%eta_u(:,:,taup1) = Grd%umask(:,:,1)*REMAP_BT_TO_BU(Ext_mode%eta_t_bar(:,:,taup1))
      call mpp_update_domains (Ext_mode%eta_u(:,:,taup1), Dom%domain2d)

      do n = 1,2
         do j = jsd, jed
            do i = isd, ied
               Ext_mode%ud(i,j,n,taup1)      = rnts*Ext_mode%ud(i,j,n,taup1)
            enddo
         enddo
      enddo
      call mpp_update_domains(Ext_mode%ud(:,:,1,taup1),Ext_mode%ud(:,:,2,taup1),Dom%domain2d,gridtype=BGRID_NE) 

      Ext_mode%convud_t(:,:,taup1)  = -DIV_UD(Ext_mode%ud(:,:,:,taup1))
      call mpp_update_domains (Ext_mode%convud_t(:,:,taup1), Dom%domain2d)

  endif

  ! Construct surface pressure gradient.  This gradient is needed  
  ! for velocity eqn when dtuv=dtfs, or when check_numerics=true.
  ! in which case we use it for energy_analysis. 
  do j=jsd,jed
     do i=isd,ied
        Ext_mode%ps(i,j) = patm(i,j) + rho_g(i,j)*Ext_mode%eta_t(i,j,tau)
     enddo
  enddo
  Ext_mode%grad_ps(:,:,:) = GRAD_SURF_P(Ext_mode%ps(:,:))

  ! send diagnostics to diag_manager     
  if (id_eta_u > 0)    used = send_data (id_eta_u, Ext_mode%eta_u(:,:,taup1), &
                              Time%model_time, rmask=Grd%umask(:,:,1), &
                              is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  if (id_ud > 0)       used = send_data (id_ud, Ext_mode%ud(:,:,1,taup1), &
                              Time%model_time, rmask=Grd%umask(:,:,1), &
                              is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  if (id_vd > 0)       used = send_data (id_vd, Ext_mode%ud(:,:,2,taup1), &
                              Time%model_time, rmask=Grd%umask(:,:,1), &
                              is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  if (id_convud_t > 0) used = send_data (id_convud_t, Ext_mode%convud_t(:,:,tau), &
                              Time%model_time, rmask=Grd%tmask(:,:,1), &
                              is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  if (id_ps > 0)       used = send_data (id_ps, Ext_mode%ps(:,:), &
                              Time%model_time, rmask=Grd%tmask(:,:,1), &
                              is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  if (id_psx > 0)      used = send_data (id_psx, Ext_mode%grad_ps(:,:,1), &
                              Time%model_time, rmask=Grd%umask(:,:,1), &
                              is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  if (id_psy > 0)      used = send_data (id_psy, Ext_mode%grad_ps(:,:,2), &
                              Time%model_time, rmask=Grd%umask(:,:,1), &
                              is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

  ! compute barotropic energetics  
  next_time = increment_time(Time%model_time, int(dtts), 0)
  if (need_data(id_ke_fs,next_time) .or. need_data(id_pe_fs,next_time)) then
    call freesurf_energy (Time, Ext_mode)
  endif 

  ! diagnose quasi-streamfunctions for vertically integrated transport
  psiu=0.0 ; psiv=0.0
  if (need_data(id_psiu, next_time) .or. need_data(id_psiv, next_time))  then 
      call psi_compute(Time, Ext_mode, psiu, psiv)
      if(id_psiu>0) used = send_data (id_psiu, psiu(:,:), &
                           Time%model_time, rmask=Grd%umask(:,:,1), &
                           is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
      if(id_psiv>0) used = send_data (id_psiv, psiv(:,:), &
                           Time%model_time, rmask=Grd%umask(:,:,1), &
                           is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif 

  ! compute checksums for debugging 
  if (debug_fs) then
     call ocean_ud_chksum(Time, Ext_mode, taup1)
  endif

end subroutine update_ocean_barotropic
! </SUBROUTINE> NAME="update_ocean_barotropic"


!#######################################################################
! <SUBROUTINE NAME="ocean_barotropic_forcing">
!
! <DESCRIPTION>
! Construct the vertically integrated forcing. This forcing is to be  
! held constant over the barotropic timesteps. At the time of calling
! this routine, accel only has the contributions from explicit-time 
! forcing, except for the following:
!
! 1. Coriolis force is updated on the barotropic time steps when 
!    integrating the barotropic dynamics.  So it should not 
!    be included in forcing_fs.  
!
! 2. Contributions from smf and bmf are added to forcing_fs to allow 
!    for simpler handling of vertical friction implicitly. 
!
! 3. The accel array is already thickness weighted, so a vertical 
!    integral is here just a vertical sum.
!
! </DESCRIPTION>
!
subroutine ocean_barotropic_forcing(Time, Velocity, accel, Ext_mode)

  type(ocean_time_type), intent(in)                 :: Time 
  type(ocean_velocity_type), intent(in)             :: Velocity
  real, intent(in), dimension(isd:ied,jsd:jed,nk,2) :: accel
  type(ocean_external_mode_type), intent(inout)     :: Ext_mode

  integer :: i, j, k, n, tau
  
  tau = Time%tau

  do n=1,2
     do j = jsd,jed
        do i = isd,ied
           Ext_mode%forcing_fs(i,j,n) = 0.0
        enddo
     enddo
     do k=1,nk 
        do j=jsc,jec
           do i=isc,iec
              Ext_mode%forcing_fs(i,j,n) = Ext_mode%forcing_fs(i,j,n) + Grd%umask(i,j,k)*accel(i,j,k,n)
           enddo
        enddo
     enddo

     do j=jsc,jec
        do i=isc,iec
           Ext_mode%forcing_fs(i,j,n) =  Ext_mode%forcing_fs(i,j,n)  &
                                       + Grd%umask(i,j,1)*(Velocity%smf(i,j,n)-Velocity%bmf(i,j,n))
        enddo
     enddo

  enddo
  
  if (have_obc) call ocean_obc_adjust_forcing_fs(Ext_mode)

  call mpp_update_domains(Ext_mode%forcing_fs(:,:,1),Ext_mode%forcing_fs(:,:,2), Dom%domain2d, gridtype=BGRID_NE)

  if (debug_fs) then
    write(stdout(),*) 'from ocean_freesurf_mod: forcing_fs chksum ==>', &
                       mpp_chksum(Ext_mode%forcing_fs(isc:iec,jsc:jec,:))
  endif

end subroutine ocean_barotropic_forcing
! </SUBROUTINE> NAME="ocean_barotropic_forcing"


!#######################################################################
! <SUBROUTINE NAME="leap_frog_barotropic">
!
! <DESCRIPTION>
! Integrate barotropic dynamics for "nts" timesteps using leapfrog.
! </DESCRIPTION>
!
subroutine leap_frog_barotropic (Time, Ext_mode, patm, pme, river)

  type(ocean_time_type), intent(in)              :: Time 
  type(ocean_external_mode_type), intent(inout)  :: Ext_mode
  real, intent(in), dimension(isd:ied,jsd:jed)   :: patm 
  real, intent(in), dimension(isd:ied,jsd:jed)   :: pme
  real, intent(in), dimension(isd:ied,jsd:jed)   :: river

  type(time_type)                  :: time_fs 
  real, dimension(isd:ied,jsd:jed) :: depthu
  real, dimension(isd:ied,jsd:jed) :: steady_forcing

  integer :: fstaum1, fstau, fstaup1 
  integer :: itime, i, j, n
  integer :: tau, taum1, taup1
  integer :: day, sec 
  real    :: dayr    
  real    :: dtime

  if ( .not. module_is_initialized ) then 
      call mpp_error(FATAL, 'Error from ocean_freesurf_mod (leap_frog_barotropic): module must be initialized')
  endif

  depthu=0.0 ; eta_eq_tidal=0.0

  taum1 = Time%taum1
  tau   = Time%tau 
  taup1 = Time%taup1

  itime   = 1
  fstaum1 = mod(itime-1,3) + 1  
  fstau   = mod(itime+0,3) + 1

  if (splitting) then
      do j=jsd,jed
         do i=isd,ied
            eta_t_fs(i,j,fstaum1)         = Ext_mode%eta_t_bar(i,j,tau)
            eta_t_fs(i,j,fstau)           = Ext_mode%eta_t_bar(i,j,tau)
            Ext_mode%eta_t_bar(i,j,taup1) = Ext_mode%eta_t(i,j,tau)     ! for accumulation of eta_t_fs
         enddo
      enddo
      do n=1,2
         do j=jsd,jed
            do i=isd,ied
               ud_fs(i,j,n,fstaum1)     = Ext_mode%ud(i,j,n,tau)
               ud_fs(i,j,n,fstau)       = Ext_mode%ud(i,j,n,tau)
               Ext_mode%ud(i,j,n,taup1) = Ext_mode%ud(i,j,n,tau)      ! for accumulation of ud_fs  
            enddo
         enddo
      enddo
  else
      do j=jsd,jed
         do i=isd,ied
            eta_t_fs(i,j,fstaum1) = Ext_mode%eta_t(i,j,taum1)
            eta_t_fs(i,j,fstau)   = Ext_mode%eta_t(i,j,tau)
         enddo
      enddo
      do n =1,2
         do j=jsd,jed
            do i=isd,ied
               ud_fs(i,j,n,fstaum1)  = Ext_mode%ud(i,j,n,taum1)
               ud_fs(i,j,n,fstau)    = Ext_mode%ud(i,j,n,tau)
            enddo
         enddo
      enddo
  endif

  ! set total depth on U cells
  do j=jsd,jed
     do i=isd,ied
        depthu(i,j) = (Grd%hu(i,j) + Ext_mode%eta_u(i,j,tau))*Grd%umask(i,j,1)
     enddo
  enddo

  if(tidal_forcing_m2 .or. tidal_forcing_8) then 
      time_fs = increment_time(Time%model_time, int(dteta),0) 
      call get_time(time_fs, sec, day)                                 
      dayr = day + sec/86400.0     
  endif

  ! compute eta_t_fs forcing which is constant over barotropic integration 
  do j=jsd,jed
     do i=isd,ied 
        steady_forcing(i,j) = pme(i,j) + river(i,j) + Ext_mode%eta_source(i,j)
     enddo
  enddo


  ! step over the nts barotropic time steps 
  do itime=1,nts

     if (tidal_forcing_m2 .or. tidal_forcing_8) then
         dayr = dayr + dtfs/86400.0 
         call get_tidal_forcing(Time, Ext_mode, dayr)
     endif

     ! set leap-frog time level indices
     fstaum1 = mod(itime-1,3) + 1 
     fstau   = mod(itime+0,3) + 1
     fstaup1 = mod(itime+1,3) + 1

     if (itime == 1) then
         if(splitting) then 
             dtime = 0.5*twodt ! first timestep is euler forward
         else 
             dtime = twodt 
         endif
     else
         dtime = twodt     ! remaining timesteps are leapfrog
     endif

     ! update free surface height on barotropic time step
     tmp(:,:) = -DIV_UD(ud_fs(:,:,:,fstau))
     do j=jsd,jed
        do i=isd,ied      
           eta_t_fs(i,j,fstaup1) = eta_t_fs(i,j,fstaum1) + dtime*(tmp(i,j) + steady_forcing(i,j))
        enddo
     enddo

     if(have_obc) then 
     ! open boundary condition at barotropic time step
     ! This message passing is needed only, if the boundary scheme needs 
     ! values from two domains. A check should be introduced, which 
     ! switches off this call, if not needed, since message passing is  
     ! done here at high frequency.
         call ocean_obc_freesurf(eta_t_fs, fstaum1, fstau, fstaup1, dtime)
     endif

     if(smooth_eta_t_fs_laplacian) then
         eta_t_fs(:,:,fstaup1) = eta_t_fs(:,:,fstaup1) &
           + dtime*LAP_T(eta_t_fs(:,:,fstaum1),eta_mix_lap(:,:))
     endif
     if(smooth_eta_t_fs_biharmonic) then
         tmp(:,:) = -LAP_T(eta_t_fs(:,:,fstaum1), eta_mix_bih(:,:))
         call mpp_update_domains (tmp(:,:), Dom%domain2d)
         eta_t_fs(:,:,fstaup1) = eta_t_fs(:,:,fstaup1)    + dtime*LAP_T(tmp(:,:), eta_mix_bih(:,:))
     endif

     call mpp_update_domains (eta_t_fs(:,:,fstaup1), Dom%domain2d)   

     if (tidal_forcing_m2 .or. tidal_forcing_8) then
         do j = jsd,jed
            do i = isd, ied
               Ext_mode%ps(i,j) =  rho_g(i,j)*(alphat*eta_t_fs(i,j,fstau)-eta_eq_tidal(i,j)) + patm(i,j)
            enddo
         enddo
     else 
         do j = jsd,jed
            do i = isd, ied
               Ext_mode%ps(i,j) =  rho_g(i,j)*eta_t_fs(i,j,fstau) + patm(i,j)
            enddo
         enddo
     endif

     Ext_mode%grad_ps(:,:,:) = GRAD_SURF_P(Ext_mode%ps(:,:))
     ! handle redundancies along bipolar fold 
     if (Grd%tripolar) then 
         call mpp_update_domains (Ext_mode%grad_ps(:,:,1),Ext_mode%grad_ps(:,:,2), Dom%domain2d, gridtype=BGRID_NE)
     endif

     ! update vertically integrated horizontal velocity using leap-frog time step 
     do j=jsd,jed
        do i=isd,ied
           ud_fs(i,j,1,fstaup1) = ud_fs(i,j,1,fstaum1) + dtime*( Grd%f(i,j)*ud_fs(i,j,2,fstau) - &
                   depthu(i,j)*Ext_mode%grad_ps(i,j,1) + Ext_mode%forcing_fs(i,j,1))
           ud_fs(i,j,2,fstaup1) = ud_fs(i,j,2,fstaum1) + dtime*(-Grd%f(i,j)*ud_fs(i,j,1,fstau) - &
                   depthu(i,j)*Ext_mode%grad_ps(i,j,2) + Ext_mode%forcing_fs(i,j,2))
        enddo
     enddo

     ! update barotropic velocity on the global halo 
     ! gradient accross boundary = 0 of cross boundary flux 
     ! minimize +/- structure for along boundary flux
     if(have_obc) then
!!$      call mpp_update_domains(ud_fs(:,:,:,fstaup1), Dom%domain2d)
         call ocean_obc_update_boundary(ud_fs(:,:,1,fstaup1),'M','s')
         call ocean_obc_update_boundary(ud_fs(:,:,2,fstaup1),'M','i')
         call ocean_obc_update_boundary(ud_fs(:,:,1,fstaup1),'Z','i')
         call ocean_obc_update_boundary(ud_fs(:,:,2,fstaup1),'Z','s')
     endif

     if(robert_asselin > 0.0) then 
         do j = jsd,jed
            do i = isd,ied
               eta_t_fs(i,j,fstau)  = eta_t_fs(i,j,fstau) + &
                     robert_asselin*(0.5*(eta_t_fs(i,j,fstaup1)+eta_t_fs(i,j,fstaum1))-eta_t_fs(i,j,fstau))
               ud_fs(i,j,:,fstau)   = ud_fs(i,j,:,fstau)  + &
                     robert_asselin*(0.5*(ud_fs(i,j,:,fstaup1)+ud_fs(i,j,:,fstaum1))-ud_fs(i,j,:,fstau))
            enddo
         enddo
     endif 
  
     if (splitting) then 
         do j=jsd,jed
            do i=isd,ied
               Ext_mode%eta_t_bar(i,j,taup1) = Ext_mode%eta_t_bar(i,j,taup1) + eta_t_fs(i,j,fstaup1)
               Ext_mode%ud(i,j,1,taup1)      = Ext_mode%ud(i,j,1,taup1)      + ud_fs(i,j,1,fstaup1)
               Ext_mode%ud(i,j,2,taup1)      = Ext_mode%ud(i,j,2,taup1)      + ud_fs(i,j,2,fstaup1)
               Ext_mode%eta_t(i,j,tau)       = eta_t_fs(i,j,fstau)
            enddo
         enddo
     else
         do j=jsd,jed
            do i=isd,ied
               Ext_mode%eta_t(i,j,taup1) = eta_t_fs(i,j,fstaup1)
               Ext_mode%eta_t(i,j,tau)   = eta_t_fs(i,j,fstau)
               Ext_mode%ud(i,j,1,taup1)  = ud_fs(i,j,1,fstaup1)
               Ext_mode%ud(i,j,2,taup1)  = ud_fs(i,j,2,fstaup1)
               Ext_mode%ud(i,j,1,tau)    = ud_fs(i,j,1,fstau)
               Ext_mode%ud(i,j,2,tau)    = ud_fs(i,j,2,fstau)
            enddo
         enddo
     endif

     if (id_eta_t_fs > 0 .and. itime==int(nts/2) ) then
       used = send_data (id_eta_t_fs, eta_t_fs(:,:,tau), &
              Time%model_time, rmask=Grd%tmask(:,:,1), &
              is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     endif 

  enddo ! end of barotropic integration itime=1,nts 


  ! closure for eta_u if not splitting barotropic from baroclinic 
  if (.not. splitting) then 
      Ext_mode%eta_u(:,:,taup1) = Grd%umask(:,:,1)*REMAP_BT_TO_BU(Ext_mode%eta_t(:,:,taup1))
      call mpp_update_domains (Ext_mode%eta_u(:,:,taup1), Dom%domain2d)
  endif


end subroutine leap_frog_barotropic
! </SUBROUTINE> NAME="leap_frog_barotropic"



!#######################################################################
! <SUBROUTINE NAME="pred_corr_barotropic">
!
! <DESCRIPTION>
! Integrate barotropic dynamics for "nts" timesteps using predictor-corrector.
!
! This scheme is more stable than the leap_frog since it can run with longer
! time steps to resolve external mode gravity waves.  It also provides 
! some extra smoothing when pred_corr_gamma > 0 and so the options 
! smooth_eta_t_fs_laplacian and smooth_eta_t_fs_biharmonic may not be
! needed.  
!
! Note that OBC has not been tested for use with predictor-corrector.
! 
! </DESCRIPTION>
!
subroutine pred_corr_barotropic (Time, Ext_mode, patm, pme, river)

  type(ocean_time_type), intent(in)              :: Time 
  type(ocean_external_mode_type), intent(inout)  :: Ext_mode
  real, intent(in), dimension(isd:ied,jsd:jed)   :: patm 
  real, intent(in), dimension(isd:ied,jsd:jed)   :: pme
  real, intent(in), dimension(isd:ied,jsd:jed)   :: river

  type(time_type)                    :: time_fs 
  real, dimension(isd:ied,jsd:jed)   :: depthu 
  real, dimension(isd:ied,jsd:jed)   :: tmp

  integer :: itime, i, j
  integer :: tau, taup1
  integer :: fstau, fstaup1 
  integer :: day, sec 
  real    :: dayr, tmp1, tmp2, ud_tmp1, ud_tmp2

  if ( .not. module_is_initialized ) then 
      call mpp_error(FATAL, 'Error from ocean_freesurf_mod (pred_corr_barotropic): module must be initialized')
  endif
 
  eta_eq_tidal=0.0

  tau   = Time%tau 
  taup1 = Time%taup1

  itime   = 1
  fstau   = mod(itime+0,2) + 1
  fstaup1 = mod(itime+1,2) + 1

  ! initialize the barotropic time stepping 
  do  j=jsd,jed
     do i=isd,ied

        eta_t_fs(i,j,fstau)   = Ext_mode%eta_t_bar(i,j,tau)
        eta_t_fs(i,j,fstaup1) = Ext_mode%eta_t_bar(i,j,tau)
        ud_fs(i,j,1,fstau)    = Ext_mode%ud(i,j,1,tau)
        ud_fs(i,j,2,fstau)    = Ext_mode%ud(i,j,2,tau)

        ! for accumulating to take time average 
        Ext_mode%eta_t_bar(i,j,taup1) = Ext_mode%eta_t(i,j,tau) 
        Ext_mode%ud(i,j,1,taup1)      = Ext_mode%ud(i,j,1,tau)  
        Ext_mode%ud(i,j,2,taup1)      = Ext_mode%ud(i,j,2,tau)  
 
        ! set total depth on U cells
        depthu(i,j) = (Grd%hu(i,j) + Ext_mode%eta_u(i,j,tau))*Grd%umask(i,j,1)

     enddo
  enddo

  if(tidal_forcing_m2 .or. tidal_forcing_8) then 
      time_fs = increment_time(Time%model_time, int(dteta),0) 
      call get_time(time_fs, sec, day)                                 
      dayr = day + sec/86400.0     
  endif

  ! step over the nts barotropic time steps 
  do itime=1,nts

     ! set predictor-corrector time level indices
     fstau   = mod(itime+0,2) + 1
     fstaup1 = mod(itime+1,2) + 1

     ! predictor step for eta
     if(pred_corr_gamma > 0.0) then 

         tmp(:,:) = -DIV_UD(ud_fs(:,:,:,fstau))          
         do j=jsd,jed
            do i=isd,ied
               eta_t_fs(i,j,fstaup1) = eta_t_fs(i,j,fstau) &
                + pred_corr_gamma*dtfs*(tmp(i,j) + pme(i,j) + river(i,j) + Ext_mode%eta_source(i,j))
            enddo
         enddo

         if(have_obc) then 
         ! open boundary condition at barotropic time step
         ! This message passing is needed only, if the boundary scheme needs 
         ! values from two domains. A check should be introduced, which 
         ! switches off this call, if not needed, since message passing is  
         ! done here at high frequency.
           call ocean_obc_freesurf(eta_t_fs, fstau, fstau, fstaup1, pred_corr_gamma*dtfs)
         endif

         call mpp_update_domains (eta_t_fs(:,:,fstaup1), Dom%domain2d)      

     endif

     if (tidal_forcing_m2 .or. tidal_forcing_8) then
         dayr = dayr + dtfs/86400.0 
         call get_tidal_forcing(Time, Ext_mode, dayr)
         do j = jsd,jed
            do i = isd,ied
               Ext_mode%ps(i,j) = rho_g(i,j)*(alphat*eta_t_fs(i,j,fstaup1)-eta_eq_tidal(i,j)) + patm(i,j)
            enddo
         enddo
     else 
         do j = jsd,jed
            do i = isd,ied
               Ext_mode%ps(i,j) = rho_g(i,j)*eta_t_fs(i,j,fstaup1) + patm(i,j)
            enddo
         enddo
     endif

     Ext_mode%grad_ps(:,:,:) = GRAD_SURF_P(Ext_mode%ps(:,:))
     ! handle redundancies along bipolar fold 
     if (Grd%tripolar) then 
         call mpp_update_domains (Ext_mode%grad_ps(:,:,1),Ext_mode%grad_ps(:,:,2), Dom%domain2d, gridtype=BGRID_NE)
     endif

     ! update with explicit time pieces and implicit Coriolis 
     do j = jsd,jed
        do i = isd, ied
           ud_tmp1 = ud_fs(i,j,1,fstau) + dtfs*( 0.5*Grd%f(i,j)*ud_fs(i,j,2,fstau) - &
                     depthu(i,j)*Ext_mode%grad_ps(i,j,1) + Ext_mode%forcing_fs(i,j,1))
           ud_tmp2 = ud_fs(i,j,2,fstau) + dtfs*(-0.5*Grd%f(i,j)*ud_fs(i,j,1,fstau) - &
                     depthu(i,j)*Ext_mode%grad_ps(i,j,2) + Ext_mode%forcing_fs(i,j,2))
           tmp1    = 0.5*dtfs*Grd%f(i,j) 
           tmp2    = 1.0/(1.0 + tmp1**2)
           ud_fs(i,j,1,fstaup1) = tmp2*(ud_tmp1 + tmp1*ud_tmp2)
           ud_fs(i,j,2,fstaup1) = tmp2*(ud_tmp2 - tmp1*ud_tmp1)
        enddo
     enddo
     ! update barotropic velocity on the global halo 
     ! gradient accross boundary = 0 of cross boundary flux 
     ! minimize +/- structure for along boundary flux
     if(have_obc) then
         call ocean_obc_update_boundary(ud_fs(:,:,1,fstaup1),'M','s')
         call ocean_obc_update_boundary(ud_fs(:,:,2,fstaup1),'M','i')
         call ocean_obc_update_boundary(ud_fs(:,:,1,fstaup1),'Z','i')
         call ocean_obc_update_boundary(ud_fs(:,:,2,fstaup1),'Z','s')
     endif

     ! corrector step for eta 
     tmp(:,:) = -DIV_UD(ud_fs(:,:,:,fstaup1))
     do j=jsd,jed
        do i=isd,ied          
           eta_t_fs(i,j,fstaup1) = eta_t_fs(i,j,fstau) + dtfs*(tmp(i,j) + pme(i,j) + river(i,j) + Ext_mode%eta_source(i,j))
        enddo
     enddo
     if(have_obc) then 
     ! open boundary condition at barotropic time step
     ! This message passing is needed only, if the boundary scheme needs 
     ! values from two domains. A check should be introduced, which 
     ! switches off this call, if not needed, since message passing is  
     ! done here at high frequency.
       call ocean_obc_freesurf(eta_t_fs, fstau, fstau, fstaup1, dtfs)
     endif

     if(smooth_eta_t_fs_laplacian) then
         eta_t_fs(:,:,fstaup1) = eta_t_fs(:,:,fstaup1)  &
                                 + dtfs*LAP_T(eta_t_fs(:,:,fstau), eta_mix_lap(:,:))
     endif
     if(smooth_eta_t_fs_biharmonic) then
         tmp(:,:) = -LAP_T(eta_t_fs(:,:,fstau), eta_mix_bih(:,:))
         call mpp_update_domains (tmp(:,:), Dom%domain2d)
         eta_t_fs(:,:,fstaup1) = eta_t_fs(:,:,fstaup1) + dtfs*LAP_T(tmp(:,:),eta_mix_bih(:,:))
     endif
 
     call mpp_update_domains (eta_t_fs(:,:,fstaup1), Dom%domain2d)      

     ! accumulate for time average 
     do j=jsd,jed
        do i=isd,ied   
           Ext_mode%eta_t_bar(i,j,taup1) = Ext_mode%eta_t_bar(i,j,taup1) + eta_t_fs(i,j,fstaup1)
           Ext_mode%ud(i,j,1,taup1)      = Ext_mode%ud(i,j,1,taup1)      + ud_fs(i,j,1,fstaup1)
           Ext_mode%ud(i,j,2,taup1)      = Ext_mode%ud(i,j,2,taup1)      + ud_fs(i,j,2,fstaup1)
        enddo
     enddo

     if (id_eta_t_fs > 0 .and. itime==int(nts/2) ) then
       used = send_data (id_eta_t_fs, eta_t_fs(:,:,tau), &
              Time%model_time, rmask=Grd%tmask(:,:,1), &
              is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     endif 

  enddo ! end of barotropic integration itime=1,nts 


end subroutine pred_corr_barotropic
! </SUBROUTINE> NAME="pred_corr_barotropic"



!#######################################################################
! <SUBROUTINE NAME="freesurf_integrals">
!
! <DESCRIPTION>
! Compute area averaged fresh water and surface height 
! </DESCRIPTION>
!
subroutine freesurf_integrals (Time, Ext_mode, pme, river)
  
  type(ocean_time_type), intent(in)            :: Time
  type(ocean_external_mode_type), intent(in)   :: Ext_mode
  real, intent(in), dimension(isd:ied,jsd:jed) :: pme
  real, intent(in), dimension(isd:ied,jsd:jed) :: river 

  real    :: pmei, riveri, etati
  integer :: i, j, tau
  type(time_type) :: next_time

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL,'==>Error from ocean_freesurf_mod (freesurf_integrals): module must be initialized')
  endif 

  next_time = increment_time(Time%model_time, int(dtts), 0)  
  
  if (diag_freq > 0 .and. mod(Time%itt, diag_freq) == 0&
       .or. need_data(id_pmei,next_time).or. need_data(id_riveri,next_time)&
       .or. need_data(id_etati,next_time)) then  
      pmei   = 0.0       
      riveri = 0.0 
      etati  = 0.0 

      tau = Time%tau

      if(have_obc) then
        do j=jsc,jec
          do i=isc,iec
            pmei   = pmei   + Grd%dat(i,j)*Grd%tmask(i,j,1)*pme(i,j)*Grd%obc_tmask(i,j)
            riveri = riveri + Grd%dat(i,j)*Grd%tmask(i,j,1)*river(i,j)*Grd%obc_tmask(i,j)
            etati  = etati  + Grd%dat(i,j)*Grd%tmask(i,j,1)*Ext_mode%eta_t(i,j,tau)*Grd%obc_tmask(i,j)
          enddo
        enddo
      else
        do j=jsc,jec
          do i=isc,iec
            pmei   = pmei   + Grd%dat(i,j)*Grd%tmask(i,j,1)*pme(i,j)
            riveri = riveri + Grd%dat(i,j)*Grd%tmask(i,j,1)*river(i,j)
            etati  = etati  + Grd%dat(i,j)*Grd%tmask(i,j,1)*Ext_mode%eta_t(i,j,tau)
          enddo
        enddo
      endif
      call mpp_sum (pmei)
      call mpp_sum (riveri)
      call mpp_sum (etati)

      pmei   = pmei*dtuv/Grd%tcella(1)
      riveri = riveri*dtuv/Grd%tcella(1)
      etati  = etati/Grd%tcella(1)

      if (id_riveri > 0) used = send_data (id_riveri, riveri, Time%model_time)
      if (id_etati > 0)  used = send_data (id_etati,  etati,  Time%model_time)
      if (id_pmei > 0)   used = send_data (id_pmei,   pmei,   Time%model_time)

  endif
  
end subroutine freesurf_integrals
! </SUBROUTINE> NAME="freesurf_integrals"



!#######################################################################
! <SUBROUTINE NAME="freesurf_energy">
!
! <DESCRIPTION>
!  Compute free surface energetics 
! </DESCRIPTION>
!
subroutine freesurf_energy (Time, Ext_mode)

  type(ocean_time_type), intent(in)          :: Time
  type(ocean_external_mode_type), intent(in) :: Ext_mode
  real :: ke_fs, pe_fs, depth
  integer :: i, j, tau

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL,'==>Error from ocean_freesurf_mod (freesurf_energy): module must be initialized')
  endif 

  tau = Time%tau

  ke_fs = 0.0 
  pe_fs = 0.0
  if(have_obc) then
    do j=jsc,jec
      do i=isc,iec
        depth = (Grd%hu(i,j) + Ext_mode%eta_u(i,j,tau))*Grd%obc_umask(i,j) + epsln
        ke_fs = ke_fs + Grd%dau(i,j)*((Ext_mode%ud(i,j,1,tau))**2 + (Ext_mode%ud(i,j,2,tau))**2)*Grd%obc_umask(i,j)/depth
        pe_fs = pe_fs + Grd%dat(i,j)*Grd%tmask(i,j,1)*Ext_mode%eta_t(i,j,tau)**2 * Grd%obc_tmask(i,j)
      enddo
    enddo
  else
    do j=jsc,jec
      do i=isc,iec
        depth = Grd%hu(i,j) + Ext_mode%eta_u(i,j,tau) + epsln
        ke_fs = ke_fs + Grd%dau(i,j)*((Ext_mode%ud(i,j,1,tau))**2 + (Ext_mode%ud(i,j,2,tau))**2)/depth
        pe_fs = pe_fs + Grd%dat(i,j)*Grd%tmask(i,j,1)*Ext_mode%eta_t(i,j,tau)**2
      enddo
    enddo
  endif

  call mpp_sum (ke_fs)
  call mpp_sum (pe_fs)

  ke_fs = ke_fs*0.5*rho0
  pe_fs = pe_fs*grav*rho0

  if (id_ke_fs > 0) used = send_data (id_ke_fs, ke_fs/1.e15, Time%model_time)
  if (id_pe_fs > 0) used = send_data (id_pe_fs, pe_fs/1.e15, Time%model_time)

end subroutine freesurf_energy
! </SUBROUTINE> NAME="freesurf_energy"


!#######################################################################
! <SUBROUTINE NAME="read_freesurf">
!
! <DESCRIPTION>
!  Read in external mode fields from restart file.
! </DESCRIPTION>
subroutine read_freesurf(Time, Ext_mode, ens_ocean)

  type(ocean_time_type), intent(in)             :: Time
  type(ocean_external_mode_type), intent(inout) :: Ext_mode
  logical, intent(in),                 optional :: ens_ocean  

  character*128 file_name
  integer               :: tau, taum1, taup1
  integer, dimension(4) :: siz 

  tau   = Time%tau
  taum1 = Time%taum1
  taup1 = Time%taup1

  file_name = 'INPUT/ocean_freesurf.res.nc'

  if(tendency==THREE_LEVEL) then 

      write (stdout(),'(a)') '  Reading THREE_LEVEL restart for velocity from INPUT/ocean_freesurf.res.nc'
      write (stdout(),'(a)') '  If not an initial condition, then expect two time records for each restart field.'

      Ext_mode%eta_t(:,:,taum1) = 0.0
      call read_data(file_name, 'eta_t',Ext_mode%eta_t(:,:,tau),&
           domain=Dom%domain2d,timelevel=1, append_pelist_name=ens_ocean)  
      call read_data(file_name, 'eta_t',Ext_mode%eta_t(:,:,taup1),&
           domain=Dom%domain2d,timelevel=2, append_pelist_name=ens_ocean)
      call mpp_update_domains(Ext_mode%eta_t(:,:,:), Dom%domain2d)

      Ext_mode%convud_t(:,:,taum1) = 0.0
      call read_data(file_name, 'convud_t',Ext_mode%convud_t(:,:,tau),&
           domain=Dom%domain2d,timelevel=1, append_pelist_name=ens_ocean)
      call read_data(file_name, 'convud_t',Ext_mode%convud_t(:,:,taup1),&
           domain=Dom%domain2d,timelevel=2, append_pelist_name=ens_ocean)
      call mpp_update_domains(Ext_mode%convud_t(:,:,:), Dom%domain2d)

      Ext_mode%eta_t_bar(:,:,taum1) = 0.0
      call read_data(file_name, 'eta_t_bar',Ext_mode%eta_t_bar(:,:,tau),&
           domain=Dom%domain2d,timelevel=1, append_pelist_name=ens_ocean)
      call read_data(file_name, 'eta_t_bar',Ext_mode%eta_t_bar(:,:,taup1),&
           domain=Dom%domain2d,timelevel=2, append_pelist_name=ens_ocean)
      call mpp_update_domains(Ext_mode%eta_t_bar(:,:,:), Dom%domain2d)

      Ext_mode%eta_u(:,:,taum1) = 0.0
      call read_data(file_name, 'eta_u',Ext_mode%eta_u(:,:,tau),&
           domain=Dom%domain2d,timelevel=1, append_pelist_name=ens_ocean)
      call read_data(file_name, 'eta_u',Ext_mode%eta_u(:,:,taup1),&
           domain=Dom%domain2d,timelevel=2, append_pelist_name=ens_ocean)
      call mpp_update_domains(Ext_mode%eta_u(:,:,:), Dom%domain2d)

      Ext_mode%ud(:,:,:,taum1) = 0.0
      call read_data(file_name, 'ud',Ext_mode%ud(:,:,1,tau),&
           domain=Dom%domain2d,timelevel=1, append_pelist_name=ens_ocean)
      call read_data(file_name, 'ud',Ext_mode%ud(:,:,1,taup1),&
           domain=Dom%domain2d,timelevel=2, append_pelist_name=ens_ocean)
      call read_data(file_name, 'vd',Ext_mode%ud(:,:,2,tau),&
           domain=Dom%domain2d,timelevel=1, append_pelist_name=ens_ocean)
      call read_data(file_name, 'vd',Ext_mode%ud(:,:,2,taup1),&
           domain=Dom%domain2d,timelevel=2, append_pelist_name=ens_ocean)
      call mpp_update_domains(Ext_mode%ud(:,:,1,tau),   Ext_mode%ud(:,:,2,tau),   Dom%domain2d, gridtype=BGRID_NE)         
      call mpp_update_domains(Ext_mode%ud(:,:,1,taup1), Ext_mode%ud(:,:,2,taup1), Dom%domain2d, gridtype=BGRID_NE)

  elseif(tendency==TWO_LEVEL) then 

      write (stdout(),'(/a)') '  Reading TWO_LEVEL restart for velocity from INPUT/ocean_freesurf.res.nc'
      write (stdout(),'(a)')  '  Expecting only one time record for each restart field.'

      call field_size(file_name,'eta_t', siz)
      if (siz(4) > 1) then
        write(stdout(),'(/a)') '==>ERROR: Attempt to read ocean_freesurf.res.nc from a 3-level time scheme (2 time records)'
        write(stdout(),'(a)')  '          when running mom4 with 2-level timestepping (only need 1 time record in restart).'
        write(stdout(),'(a)')  '          Reduce restart file to only a single time record in order to avoid confusion.'
        call mpp_error(FATAL, &
                       'Reading 3-time level ocean_freesurf.res.nc (w/ 2 time records) while using 2-level (needs only 1 record)')
      endif

      Ext_mode%eta_t(:,:,:) = 0.0
      call read_data(file_name, 'eta_t',Ext_mode%eta_t(:,:,taup1),&
           domain=Dom%domain2d,timelevel=1, append_pelist_name=ens_ocean)
      call mpp_update_domains(Ext_mode%eta_t(:,:,:), Dom%domain2d)

      Ext_mode%convud_t(:,:,:) = 0.0
      call read_data(file_name, 'convud_t',Ext_mode%convud_t(:,:,taup1),&
           domain=Dom%domain2d,timelevel=1, append_pelist_name=ens_ocean)
      call mpp_update_domains(Ext_mode%convud_t(:,:,:), Dom%domain2d)

      Ext_mode%eta_t_bar(:,:,:) = 0.0
      call read_data(file_name, 'eta_t_bar',Ext_mode%eta_t_bar(:,:,taup1),&
           domain=Dom%domain2d,timelevel=1, append_pelist_name=ens_ocean)
      call mpp_update_domains(Ext_mode%eta_t_bar(:,:,:), Dom%domain2d)

      Ext_mode%eta_u(:,:,:) = 0.0
      call read_data(file_name, 'eta_u',Ext_mode%eta_u(:,:,taup1),&
           domain=Dom%domain2d,timelevel=1, append_pelist_name=ens_ocean)
      call mpp_update_domains(Ext_mode%eta_u(:,:,:), Dom%domain2d)

      Ext_mode%ud(:,:,:,:) = 0.0
      call read_data(file_name, 'ud',Ext_mode%ud(:,:,1,taup1),&
           domain=Dom%domain2d,timelevel=1, append_pelist_name=ens_ocean)
      call read_data(file_name, 'vd',Ext_mode%ud(:,:,2,taup1),&
           domain=Dom%domain2d,timelevel=1, append_pelist_name=ens_ocean)
      call mpp_update_domains(Ext_mode%ud(:,:,1,taup1), Ext_mode%ud(:,:,2,taup1), Dom%domain2d, gridtype=BGRID_NE)

  endif


  if(have_obc) then
     call ocean_obc_update_boundary(Ext_mode%ud(:,:,1,tau),'M','s')
     call ocean_obc_update_boundary(Ext_mode%ud(:,:,2,tau),'M','i')
     call ocean_obc_update_boundary(Ext_mode%ud(:,:,1,tau),'Z','i')
     call ocean_obc_update_boundary(Ext_mode%ud(:,:,2,tau),'Z','s')
     call ocean_obc_update_boundary(Ext_mode%ud(:,:,1,taup1),'M','s')
     call ocean_obc_update_boundary(Ext_mode%ud(:,:,2,taup1),'M','i')
     call ocean_obc_update_boundary(Ext_mode%ud(:,:,1,taup1),'Z','i')
     call ocean_obc_update_boundary(Ext_mode%ud(:,:,2,taup1),'Z','s')
  endif

end subroutine read_freesurf
! </SUBROUTINE> NAME="read_freesurf"



!#######################################################################
! <SUBROUTINE NAME="ocean_freesurf_rstrt">
!
! <DESCRIPTION>
!  Write out external mode fields to restart file.
! </DESCRIPTION>
subroutine ocean_freesurf_rstrt(Time, Ext_mode, ens_ocean)

  type(ocean_time_type), intent(in)             :: Time
  type(ocean_external_mode_type), intent(inout) :: Ext_mode
  logical, intent(in), optional                 :: ens_ocean

  integer :: tau, taup1
  integer :: yr, mon, day, hr, min, sec
  character(len=128) :: file_name
  character(len=10)  :: rdte

  if (.not.module_is_initialized) then
    call mpp_error(FATAL,'==>Error from ocean_freesurf_mod (ocean_freesurf_rstrt): module must be initialized')
  endif

  tau   = Time%tau
  taup1 = Time%taup1

  call get_date(Time%model_time,yr,mon,day,hr,min,sec)
  write(rdte,'(i4,3i2.2)') yr,mon,day,hr
  file_name = 'IRESTART/'//rdte//'ocean_freesurf.res'

  write(stdout(),*) ' '
  write(stdout(),*) 'Ending external mode checksums'
  call write_timestamp(Time%model_time)

  if(tendency==THREE_LEVEL) then

      call write_data(file_name,'eta_t',    Ext_mode%eta_t(:,:,tau),      domain=Dom%domain2d,append_pelist_name=ens_ocean)
      call write_data(file_name,'eta_t',    Ext_mode%eta_t(:,:,taup1),    domain=Dom%domain2d,append_pelist_name=ens_ocean)
      call write_data(file_name,'convud_t', Ext_mode%convud_t(:,:,tau),   domain=Dom%domain2d,append_pelist_name=ens_ocean)
      call write_data(file_name,'convud_t', Ext_mode%convud_t(:,:,taup1), domain=Dom%domain2d,append_pelist_name=ens_ocean)
      call write_data(file_name,'eta_t_bar',Ext_mode%eta_t_bar(:,:,tau),  domain=Dom%domain2d,append_pelist_name=ens_ocean)
      call write_data(file_name,'eta_t_bar',Ext_mode%eta_t_bar(:,:,taup1),domain=Dom%domain2d,append_pelist_name=ens_ocean)
      call write_data(file_name,'eta_u',    Ext_mode%eta_u(:,:,tau),      domain=Dom%domain2d,append_pelist_name=ens_ocean)
      call write_data(file_name,'eta_u',    Ext_mode%eta_u(:,:,taup1),    domain=Dom%domain2d,append_pelist_name=ens_ocean)
      call write_data(file_name,'ud',       Ext_mode%ud(:,:,1,tau),       domain=Dom%domain2d,append_pelist_name=ens_ocean)
      call write_data(file_name,'ud',       Ext_mode%ud(:,:,1,taup1),     domain=Dom%domain2d,append_pelist_name=ens_ocean)
      call write_data(file_name,'vd',       Ext_mode%ud(:,:,2,tau),       domain=Dom%domain2d,append_pelist_name=ens_ocean)
      call write_data(file_name,'vd',       Ext_mode%ud(:,:,2,taup1),     domain=Dom%domain2d,append_pelist_name=ens_ocean)

      write(stdout(),*) 'chksum for eta_t(tau)       = ',  mpp_chksum(Ext_mode%eta_t(isc:iec,jsc:jec,tau))
      write(stdout(),*) 'chksum for eta_t(taup1)     = ',  mpp_chksum(Ext_mode%eta_t(isc:iec,jsc:jec,taup1))
      write(stdout(),*) 'chksum for convud_t(tau)    = ',  mpp_chksum(Ext_mode%convud_t(isc:iec,jsc:jec,tau))
      write(stdout(),*) 'chksum for convud_t(taup1)  = ',  mpp_chksum(Ext_mode%convud_t(isc:iec,jsc:jec,taup1))
      write(stdout(),*) 'chksum for eta_t_bar(tau)   = ',  mpp_chksum(Ext_mode%eta_t_bar(isc:iec,jsc:jec,tau))
      write(stdout(),*) 'chksum for eta_t_bar(taup1) = ',  mpp_chksum(Ext_mode%eta_t_bar(isc:iec,jsc:jec,taup1))
      write(stdout(),*) 'chksum for eta_u(tau)       = ',  mpp_chksum(Ext_mode%eta_u(isc:iec,jsc:jec,tau))
      write(stdout(),*) 'chksum for eta_u(taup1)     = ',  mpp_chksum(Ext_mode%eta_u(isc:iec,jsc:jec,taup1))
      write(stdout(),*) 'chksum for ud(tau)          = ',  mpp_chksum(Ext_mode%ud(isc:iec,jsc:jec,1,tau))
      write(stdout(),*) 'chksum for ud(taup1)        = ',  mpp_chksum(Ext_mode%ud(isc:iec,jsc:jec,1,taup1))
      write(stdout(),*) 'chksum for vd(tau)          = ',  mpp_chksum(Ext_mode%ud(isc:iec,jsc:jec,2,tau))
      write(stdout(),*) 'chksum for vd(taup1)        = ',  mpp_chksum(Ext_mode%ud(isc:iec,jsc:jec,2,taup1))

  elseif(tendency==TWO_LEVEL) then

      call write_data(file_name,'eta_t',    Ext_mode%eta_t(:,:,taup1),    domain=Dom%domain2d,append_pelist_name=ens_ocean)
      call write_data(file_name,'convud_t', Ext_mode%convud_t(:,:,taup1), domain=Dom%domain2d,append_pelist_name=ens_ocean)
      call write_data(file_name,'eta_t_bar',Ext_mode%eta_t_bar(:,:,taup1),domain=Dom%domain2d,append_pelist_name=ens_ocean)
      call write_data(file_name,'eta_u',    Ext_mode%eta_u(:,:,taup1),    domain=Dom%domain2d,append_pelist_name=ens_ocean)
      call write_data(file_name,'ud',       Ext_mode%ud(:,:,1,taup1),     domain=Dom%domain2d,append_pelist_name=ens_ocean)
      call write_data(file_name,'vd',       Ext_mode%ud(:,:,2,taup1),     domain=Dom%domain2d,append_pelist_name=ens_ocean)

      write(stdout(),*) 'chksum for eta_t(taup1)     = ',  mpp_chksum(Ext_mode%eta_t(isc:iec,jsc:jec,taup1))
      write(stdout(),*) 'chksum for convud_t(taup1)  = ',  mpp_chksum(Ext_mode%convud_t(isc:iec,jsc:jec,taup1))
      write(stdout(),*) 'chksum for eta_t_bar(taup1) = ',  mpp_chksum(Ext_mode%eta_t_bar(isc:iec,jsc:jec,taup1))
      write(stdout(),*) 'chksum for eta_u(taup1)     = ',  mpp_chksum(Ext_mode%eta_u(isc:iec,jsc:jec,taup1))
      write(stdout(),*) 'chksum for ud(taup1)        = ',  mpp_chksum(Ext_mode%ud(isc:iec,jsc:jec,1,taup1))
      write(stdout(),*) 'chksum for vd(taup1)        = ',  mpp_chksum(Ext_mode%ud(isc:iec,jsc:jec,2,taup1))

  endif

end subroutine ocean_freesurf_rstrt
! </SUBROUTINE> NAME="ocean_freesurf_rstrt"


!#######################################################################
! <SUBROUTINE NAME="ocean_freesurf_end">
!
! <DESCRIPTION>
!  Write out external mode fields to restart file. 
! </DESCRIPTION>
subroutine ocean_freesurf_end(Time, Ext_mode, ens_ocean)
  
  type(ocean_time_type), intent(in)             :: Time
  type(ocean_external_mode_type), intent(inout) :: Ext_mode
  logical, intent(in), optional                 :: ens_ocean

  character*128 file_name
  integer :: tau, taup1

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL,'==>Error from ocean_freesurf_mod (ocean_freesurf_end): module must be initialized')
  endif 

  tau   = Time%tau
  taup1 = Time%taup1

  file_name = 'RESTART/ocean_freesurf.res'

  write(stdout(),*) ' '
  write(stdout(),*) 'Ending external mode checksums'
  call write_timestamp(Time%model_time)

  if(tendency==THREE_LEVEL) then 

      call write_data(file_name,'eta_t',    Ext_mode%eta_t(:,:,tau),      domain=Dom%domain2d,append_pelist_name=ens_ocean)
      call write_data(file_name,'eta_t',    Ext_mode%eta_t(:,:,taup1),    domain=Dom%domain2d,append_pelist_name=ens_ocean)
      call write_data(file_name,'convud_t', Ext_mode%convud_t(:,:,tau),   domain=Dom%domain2d,append_pelist_name=ens_ocean)
      call write_data(file_name,'convud_t', Ext_mode%convud_t(:,:,taup1), domain=Dom%domain2d,append_pelist_name=ens_ocean)
      call write_data(file_name,'eta_t_bar',Ext_mode%eta_t_bar(:,:,tau),  domain=Dom%domain2d,append_pelist_name=ens_ocean)
      call write_data(file_name,'eta_t_bar',Ext_mode%eta_t_bar(:,:,taup1),domain=Dom%domain2d,append_pelist_name=ens_ocean)
      call write_data(file_name,'eta_u',    Ext_mode%eta_u(:,:,tau),      domain=Dom%domain2d,append_pelist_name=ens_ocean)
      call write_data(file_name,'eta_u',    Ext_mode%eta_u(:,:,taup1),    domain=Dom%domain2d,append_pelist_name=ens_ocean)
      call write_data(file_name,'ud',       Ext_mode%ud(:,:,1,tau),       domain=Dom%domain2d,append_pelist_name=ens_ocean)
      call write_data(file_name,'ud',       Ext_mode%ud(:,:,1,taup1),     domain=Dom%domain2d,append_pelist_name=ens_ocean)
      call write_data(file_name,'vd',       Ext_mode%ud(:,:,2,tau),       domain=Dom%domain2d,append_pelist_name=ens_ocean)
      call write_data(file_name,'vd',       Ext_mode%ud(:,:,2,taup1),     domain=Dom%domain2d,append_pelist_name=ens_ocean)

      write(stdout(),*) 'chksum for eta_t(tau)       = ',  mpp_chksum(Ext_mode%eta_t(isc:iec,jsc:jec,tau))
      write(stdout(),*) 'chksum for eta_t(taup1)     = ',  mpp_chksum(Ext_mode%eta_t(isc:iec,jsc:jec,taup1))
      write(stdout(),*) 'chksum for convud_t(tau)    = ',  mpp_chksum(Ext_mode%convud_t(isc:iec,jsc:jec,tau))
      write(stdout(),*) 'chksum for convud_t(taup1)  = ',  mpp_chksum(Ext_mode%convud_t(isc:iec,jsc:jec,taup1))
      write(stdout(),*) 'chksum for eta_t_bar(tau)   = ',  mpp_chksum(Ext_mode%eta_t_bar(isc:iec,jsc:jec,tau))
      write(stdout(),*) 'chksum for eta_t_bar(taup1) = ',  mpp_chksum(Ext_mode%eta_t_bar(isc:iec,jsc:jec,taup1))
      write(stdout(),*) 'chksum for eta_u(tau)       = ',  mpp_chksum(Ext_mode%eta_u(isc:iec,jsc:jec,tau))
      write(stdout(),*) 'chksum for eta_u(taup1)     = ',  mpp_chksum(Ext_mode%eta_u(isc:iec,jsc:jec,taup1))
      write(stdout(),*) 'chksum for ud(tau)          = ',  mpp_chksum(Ext_mode%ud(isc:iec,jsc:jec,1,tau))
      write(stdout(),*) 'chksum for ud(taup1)        = ',  mpp_chksum(Ext_mode%ud(isc:iec,jsc:jec,1,taup1))
      write(stdout(),*) 'chksum for vd(tau)          = ',  mpp_chksum(Ext_mode%ud(isc:iec,jsc:jec,2,tau))
      write(stdout(),*) 'chksum for vd(taup1)        = ',  mpp_chksum(Ext_mode%ud(isc:iec,jsc:jec,2,taup1))

  elseif(tendency==TWO_LEVEL) then 

      call write_data(file_name,'eta_t',    Ext_mode%eta_t(:,:,taup1),    domain=Dom%domain2d,append_pelist_name=ens_ocean)
      call write_data(file_name,'convud_t', Ext_mode%convud_t(:,:,taup1), domain=Dom%domain2d,append_pelist_name=ens_ocean)
      call write_data(file_name,'eta_t_bar',Ext_mode%eta_t_bar(:,:,taup1),domain=Dom%domain2d,append_pelist_name=ens_ocean)
      call write_data(file_name,'eta_u',    Ext_mode%eta_u(:,:,taup1),    domain=Dom%domain2d,append_pelist_name=ens_ocean)
      call write_data(file_name,'ud',       Ext_mode%ud(:,:,1,taup1),     domain=Dom%domain2d,append_pelist_name=ens_ocean)
      call write_data(file_name,'vd',       Ext_mode%ud(:,:,2,taup1),     domain=Dom%domain2d,append_pelist_name=ens_ocean)

      write(stdout(),*) 'chksum for eta_t(taup1)     = ',  mpp_chksum(Ext_mode%eta_t(isc:iec,jsc:jec,taup1))
      write(stdout(),*) 'chksum for convud_t(taup1)  = ',  mpp_chksum(Ext_mode%convud_t(isc:iec,jsc:jec,taup1))
      write(stdout(),*) 'chksum for eta_t_bar(taup1) = ',  mpp_chksum(Ext_mode%eta_t_bar(isc:iec,jsc:jec,taup1))
      write(stdout(),*) 'chksum for eta_u(taup1)     = ',  mpp_chksum(Ext_mode%eta_u(isc:iec,jsc:jec,taup1))
      write(stdout(),*) 'chksum for ud(taup1)        = ',  mpp_chksum(Ext_mode%ud(isc:iec,jsc:jec,1,taup1))
      write(stdout(),*) 'chksum for vd(taup1)        = ',  mpp_chksum(Ext_mode%ud(isc:iec,jsc:jec,2,taup1))

  endif 

  module_is_initialized = .FALSE.

  nullify(Grd)
  nullify(Dom)

end subroutine ocean_freesurf_end
! </SUBROUTINE> NAME="ocean_freesurf_end"



!#######################################################################
! <SUBROUTINE NAME="eta_truncate">
!
! <DESCRIPTION>
!
!  Truncate eta_t to keep 
!
!  dzt(1) + eta_t >= frac_crit_cell_height*dzt(1)
!
!  and 
!
!  eta_t <= eta_max
!
! </DESCRIPTION>
!
subroutine eta_truncate(Time, Ext_mode)

  type(ocean_time_type), intent(in)             :: Time
  type(ocean_external_mode_type), intent(inout) :: Ext_mode

  integer :: i, j
  integer :: taup1
  real    :: cell_thickness

  if (.not. module_is_initialized) then 
    call mpp_error(FATAL,'==>Error from ocean_freesurf_mod (eta_truncate): module must be initialized')
  endif 

  taup1 = Time%taup1

  do j=jsc,jec
     do i=isc,iec
        cell_thickness = Grd%tmask(i,j,1)*Ext_mode%eta_t(i,j,taup1)+Grd%dzt(1)
        if (cell_thickness < Grd%dzt(1)*frac_crit_cell_height) then
            if (verbose_truncate) then 
                write(*,*) 'WARNING from ocean_freesurf_mod: eta_t(',i+Dom%ioff,',',j+Dom%joff,') = ', &
                Ext_mode%eta_t(i,j,taup1), '(m) has been truncated to', -Grd%dzt(1)*(1.0-frac_crit_cell_height)
            endif
            Ext_mode%eta_t(i,j,taup1) = -Grd%dzt(1)*(1.0-frac_crit_cell_height)
        elseif (Grd%tmask(i,j,1)*Ext_mode%eta_t(i,j,taup1) > eta_max) then
            if (verbose_truncate) then
                write(*,*) 'WARNING from ocean_freesurf_mod: eta_t(',i+Dom%ioff,',',j+Dom%joff,') = ', &
                Ext_mode%eta_t(i,j,taup1), '(m) has been truncated to', eta_max
            endif
            Ext_mode%eta_t(i,j,taup1) = eta_max
        endif
     enddo
  enddo


end subroutine eta_truncate
! </SUBROUTINE> NAME="eta_truncate"


!#######################################################################
! <SUBROUTINE NAME="maximum_convud">
!
! <DESCRIPTION>
! Compute maximum convergence(ud,vd).
! </DESCRIPTION>
!
subroutine maximum_convud(Time, Ext_mode)

  type(ocean_time_type), intent(in)             :: Time
  type(ocean_external_mode_type), intent(inout) :: Ext_mode

  real :: convudt0, convudu0, convudt, convudu, fudge
  integer :: i, j, k
  integer :: iconvudt, jconvudt, iconvudu, jconvudu
  integer :: tau

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL,'==>Error ocean_freesurf_mod (maximum_convud): module needs initialization ')
  endif 

  tau = Time%tau

  ! find max convergence(ud,vd) on T-cells

  fudge = 1 + 1.e-12*mpp_pe() ! to distinguish processors when convud is independent of processor
  convudt=0.0; iconvudt=isc; jconvudt=jsc

  do j=jsc,jec
    do i=isc,iec
      
      if (abs(Ext_mode%convud_t(i,j,tau)) > abs(convudt) .and. Grd%tmask(i,j,1) == 1.0) then
        convudt  = Ext_mode%convud_t(i,j,tau)
        iconvudt = i
        jconvudt = j
      endif
    enddo
  enddo

  ! find max convergence(ud,vd) on U-cells
  tmp = Grd%umask(:,:,1)*REMAP_BT_TO_BU(Ext_mode%convud_t(:,:,tau))

  convudu=0.0; iconvudu=isc; jconvudu=jsc
  do j=jsc,jec
    do i=isc,iec
      if (abs(tmp(i,j)) > abs(convudu) .and. Grd%umask(i,j,1) == 1.0) then 
        convudu  = tmp(i,j)
        iconvudu = i
        jconvudu = j
      endif
    enddo
  enddo
  
  write (stdout(),'(//60x,a/)') ' Convergence of depth integrated horz velocity summary:'

  convudt  = convudt*fudge
  convudu  = convudu*fudge
  convudt0 = convudt
  convudu0 = convudu
  convudt  = abs(convudt)
  convudu  = abs(convudu)
  call mpp_max(convudt)
  call mpp_max(convudu)

  if (abs(convudt0) == convudt .and. abs(convudt0) /= 0.0) then
    convudt = convudt0
    write (unit,9112) convudt/fudge, iconvudt, jconvudt, Grd%xt(iconvudt,jconvudt), Grd%yt(iconvudt,jconvudt)
  endif
  
  if (abs(convudu0) == convudu .and. abs(convudu0) /= 0.0) then
    convudu = convudu0
    write (unit,9113) convudu/fudge, iconvudu, jconvudu, Grd%xu(iconvudu,jconvudu), Grd%yu(iconvudu,jconvudu)
  endif


9112  format(/' Maximum at T-cell convud_t (',es10.3,' m/s) at (i,j) = ','(',i4,',',i4,'),',&
          ' (lon,lat) = (',f7.2,',',f7.2,')')
9113  format(' Maximum at U-cell convud_u (',es10.3,' m/s) at (i,j) = ','(',i4,',',i4,'),',&
          ' (lon,lat) = (',f7.2,',',f7.2,')'/)

end subroutine maximum_convud
! </SUBROUTINE> NAME="maximum_convud"


!#######################################################################
! <SUBROUTINE NAME="ocean_ud_chksum">
!
! <DESCRIPTION>
!  Compute checksum for ud.
! </DESCRIPTION>
!
subroutine ocean_ud_chksum(Time, Ext_mode, index)

  type(ocean_time_type), intent(in)          :: Time
  type(ocean_external_mode_type), intent(in) :: Ext_mode
  integer, intent(in) :: index

  write(stdout(),*) '===Ext_mode%ud chksum follows==='
  call write_timestamp(Time%model_time)
  tmp(isc:iec,jsc:jec) = Ext_mode%ud(isc:iec,jsc:jec,1,index)*Grd%umask(isc:iec,jsc:jec,1)
  write(stdout(),*) 'ud chksum    = ',  mpp_chksum(tmp(isc:iec,jsc:jec))
  tmp(isc:iec,jsc:jec) = Ext_mode%ud(isc:iec,jsc:jec,2,index)*Grd%umask(isc:iec,jsc:jec,1)
  write(stdout(),*) 'vd chksum    = ',  mpp_chksum(tmp(isc:iec,jsc:jec))

  return

end subroutine ocean_ud_chksum
! </SUBROUTINE> NAME="ocean_ud_chksum"


!#######################################################################
! <SUBROUTINE NAME="ocean_eta_chksum">
!
! <DESCRIPTION>
!  Compute checksum for surface height. 
! </DESCRIPTION>
!
subroutine ocean_eta_chksum(Time, Ext_mode, index)

  type(ocean_time_type), intent(in)          :: Time
  type(ocean_external_mode_type), intent(in) :: Ext_mode
  integer, intent(in) :: index

  write(stdout(),*) '===Ext_mode%eta chksum follows==='
  call write_timestamp(Time%model_time)
  tmp(isc:iec,jsc:jec) = Ext_mode%eta_t(isc:iec,jsc:jec,index)*Grd%tmask(isc:iec,jsc:jec,1)
  write(stdout(),*) 'eta_t chksum = ',  mpp_chksum(tmp(isc:iec,jsc:jec))
  tmp(isc:iec,jsc:jec) = Ext_mode%eta_u(isc:iec,jsc:jec,index)*Grd%umask(isc:iec,jsc:jec,1)
  write(stdout(),*) 'eta_u chksum = ',  mpp_chksum(tmp(isc:iec,jsc:jec))

  return

end subroutine ocean_eta_chksum
! </SUBROUTINE> NAME="ocean_eta_chksum"


!#######################################################################
! <SUBROUTINE NAME="psi_compute">
!
! <DESCRIPTION>
!  Compute quasi-barotropic streamfunctions for diagnostic purposes.
!  When no fresh water and steady state, these two streamfunctions 
!  will be equal, and they will be equal to the rigid lid barotropic
!  streamfunction.  Note that for plotting purposes, it is necessary 
!  to remove a global constant, usually taken as the value over 
!  the American continent. 
!
!  Original algorithm: Stephen.Griffies@noaa.gov
!  Modifications for parallel efficiency: Giang.Nong@noaa.gov
!
! </DESCRIPTION>
!
subroutine psi_compute(Time, Ext_mode, psiu, psiv)

  type(ocean_time_type), intent(in)          :: Time
  type(ocean_external_mode_type), intent(in) :: Ext_mode
  real, intent(inout)                        :: psiu(isd:ied,jsd:jed)
  real, intent(inout)                        :: psiv(isd:ied,jsd:jed)

  real                                       :: psiu2(isd:ied), psiv2(jsd:jed) 
  real                                       :: mkstoSv
  integer                                    :: lenx,leny
  integer                                    :: i,j,ii
  integer                                    :: tau

  tau     = Time%tau
  mkstoSv = 1.e-6
  psiu    = 0.0
  psiv    = 0.0
  psiu2   = 0.0
  psiv2   = 0.0 
  lenx    = iec-isc+1
  leny    = jec-jsc+1
  
  do j=jsc,jec
     do i=isc,iec
        psiu(i,j) = -Grd%dyu(i,j)*Ext_mode%ud(i,j,1,tau)*Grd%umask(i,j,1)*mkstoSv
        psiv(i,j) =  Grd%dxu(i,j)*Ext_mode%ud(i,j,2,tau)*Grd%umask(i,j,1)*mkstoSv
     enddo
  enddo

  if(Dom%joff == 0 .and. jsc==1) psiu(:,jsc)=0.0  
  do i=isc,iec
     do j=jsc+1,jec
        psiu(i,j) = psiu(i,j) + psiu(i,j-1)
     enddo
  enddo
  psiu2(:) =  psiu(:,jec)

  ! send psiu2 to other PEs in the same column with higher rank 
  if(myid_y<size_y) then
     do ii = myid_y+1,size_y
        call mpp_send(psiu2(isc:iec),lenx,pelist_y(ii))
     enddo
  endif
  call mpp_sync_self

  ! receive psiu2 from other PEs in the same column with lower rank
  if(myid_y>1) then
     do ii = 1,myid_y-1
        call mpp_recv(psiu2(isc:iec),lenx,pelist_y(ii))
         do i=isc,iec
            do j=jsc,jec
               psiu(i,j) = psiu(i,j) + psiu2(i)
            enddo
         enddo  
     enddo
  endif

  ! compute psiv
  if(Dom%ioff == 0 .and. isc==1 ) psiv(1,:) = psiu(1,:) 
  do j=jsc,jec
     do i=isc+1,iec
        psiv(i,j) = psiv(i,j)+psiv(i-1,j)
     enddo
  enddo
  psiv2(:) = psiv(iec,:)

  ! send psiv2 to other PEs in the same row with higher rank
  if(myid_x<size_x) then
     do ii = myid_x+1,size_x
        call mpp_send(psiv2(jsc:jec),leny,pelist_x(ii))
     enddo
  endif
  call mpp_sync_self

  ! receive psiv2 from other PEs in the same row with lower rank 
  if(myid_x>1) then
     do ii = 1,myid_x-1
        call mpp_recv(psiv2(jsc:jec),leny,pelist_x(ii))         
            do j=jsc,jec
               do i=isc,iec 
               psiv(i,j) = psiv(i,j) + psiv2(j)
            enddo
         enddo  
     enddo
  endif


end subroutine psi_compute
! </SUBROUTINE> NAME="psi_compute"



!#######################################################################
! <SUBROUTINE NAME="eta_check">
!
! <DESCRIPTION>
!  Perform diagnostic check on top cell thickness.
! </DESCRIPTION>
subroutine eta_check(Time, Ext_mode)

  type(ocean_time_type), intent(in)          :: Time
  type(ocean_external_mode_type), intent(in) :: Ext_mode

  integer :: i_thin_t_cell, j_thin_t_cell, i_thick_t_cell, j_thick_t_cell
  integer :: i_thin_u_cell, j_thin_u_cell, i_thick_u_cell, j_thick_u_cell
  integer :: i, j
  real :: thin_t_cell, thin_u_cell, thick_t_cell, thick_u_cell, cell_height, crit_cell_height
  real :: thin_t0, thin_u0, thick_t0, thick_u0, fudge
  integer :: tau, taup1

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL,'==>Error from ocean_freesurf_mod (eta_check): module must be initialized')
  endif 

  tau   = Time%tau
  taup1 = Time%taup1

  thin_t_cell=Grd%dzt(1); thick_t_cell=Grd%dzt(1); thin_u_cell=Grd%dzt(1); thick_u_cell=Grd%dzt(1)
  
  i_thin_t_cell=isc;  j_thin_t_cell=jsc; i_thick_t_cell=isc; j_thick_t_cell=jsc
  i_thin_u_cell=isc;  j_thin_u_cell=jsc; i_thick_u_cell=isc; j_thick_u_cell=jsc
  crit_cell_height = frac_crit_cell_height*Grd%dzt(1)

  do j=jsc,jec
    do i=isc,iec
      cell_height = Ext_mode%eta_t(i,j,taup1)+Grd%dzt(1)    
      if (cell_height < thin_t_cell .and. Grd%tmask(i,j,1) == 1) then
        i_thin_t_cell = i
        j_thin_t_cell = j
        thin_t_cell   = cell_height
      endif 
    
      if (cell_height > thick_t_cell) then
        i_thick_t_cell = i
        j_thick_t_cell = j
        thick_t_cell   = cell_height
      endif 

      cell_height = Ext_mode%eta_u(i,j,taup1)+Grd%dzt(1)    
      if (cell_height < thin_u_cell .and. Grd%umask(i,j,1) == 1) then
        i_thin_u_cell = i
        j_thin_u_cell = j
        thin_u_cell   = cell_height
      endif 
    
      if (cell_height > thick_u_cell) then
        i_thick_u_cell = i
        j_thick_u_cell = j
        thick_u_cell   = cell_height
      endif         
    enddo
  enddo

  write (stdout(),'(//60x,a/)') ' Surface cell thickness summary:  '
  fudge = 1 + 1.e-12*mpp_pe()  ! to distinguish processors when tracer is independent of processor
  thin_t_cell = thin_t_cell*fudge
  thin_u_cell = thin_u_cell*fudge
  thick_t_cell = thick_t_cell*fudge
  thick_u_cell = thick_u_cell*fudge
  
  thin_t0 = thin_t_cell
  thin_u0 = thin_u_cell
  thick_t0 = thick_t_cell
  thick_u0 = thick_u_cell
  call mpp_min(thin_t_cell)
  call mpp_min(thin_u_cell)
  call mpp_max(thick_t_cell)
  call mpp_max(thick_u_cell)

  if (thin_t0 == thin_t_cell) then
    write (unit,7000) ' Minimum thickness surface T cell = ',thin_t_cell/fudge, &
                        i_thin_t_cell+Dom%ioff, j_thin_t_cell+Dom%joff, &
                        Grd%xt(i_thin_t_cell,j_thin_t_cell), Grd%yt(i_thin_t_cell,j_thin_t_cell)
  endif
  
  if (thick_t0 == thick_t_cell) then
    write (unit,7000) ' Maximum thickness surface T cell = ',thick_t_cell/fudge, &
                        i_thick_t_cell+Dom%ioff, j_thick_t_cell+Dom%joff, &
                        Grd%xt(i_thick_t_cell,j_thick_t_cell), Grd%yt(i_thick_t_cell,j_thick_t_cell)
  endif

  if (thin_u0 == thin_u_cell) then
    write (unit,7000) ' Minimum thickness surface U cell = ',thin_u_cell/fudge, &
                        i_thin_u_cell+Dom%ioff, j_thin_u_cell+Dom%joff, &
                        Grd%xu(i_thin_u_cell,j_thin_u_cell), Grd%yu(i_thin_u_cell,j_thin_u_cell)
  endif

  if (thick_u0 == thick_u_cell) then
    write (unit,7000) ' Maximum thickness surface U cell = ',thick_u_cell/fudge, &
                        i_thick_u_cell+Dom%ioff, j_thick_u_cell+Dom%joff, &
                        Grd%xu(i_thick_u_cell,j_thick_u_cell), Grd%yu(i_thick_u_cell,j_thick_u_cell)
  endif

  if (thin_t_cell < crit_cell_height) then
    write (stdout(),*) '==>Error: minimum surface cell thickness is less than ',crit_cell_height,' meters.'
    call mpp_error(FATAL, '==>Error: minimum surface cell thickness is below critical value.')
  endif 

  7000 format (1x,a,e14.7,' m at (i,j) = (',i4,',',i4,'), (lon,lat) = (',f7.2,',',f7.2,')')

end subroutine eta_check
! </SUBROUTINE> NAME="eta_check"



!#######################################################################
! <SUBROUTINE NAME="tidal_forcing_init">
!
! <DESCRIPTION>
! Initialize fields needed for lunar and solar tidal forcing.  
! </DESCRIPTION>
!
subroutine tidal_forcing_init(Time)

  type(ocean_time_type), intent(in) :: Time

#ifndef STATIC_MEMORY
     allocate (eta_eq_tidal(isd:ied,jsd:jed))
     allocate (cos2lon(isd:ied,jsd:jed))
     allocate (coslon(isd:ied,jsd:jed))
     allocate (coslat2(isd:ied,jsd:jed))
     allocate (sin2lon(isd:ied,jsd:jed))
     allocate (sinlon(isd:ied,jsd:jed))
     allocate (sin2lat(isd:ied,jsd:jed))
#endif
     eta_eq_tidal=0.0 ; cos2lon=0.0 ; coslon=0.0 ; coslat2=0.0 ; sin2lon=0.0 ; sinlon=0.0 ; sin2lat=0.0

     if(tidal_forcing_m2) then
       call mpp_error(NOTE, &
        '==>Note from tidal_forcing_init: adding M2 lunar tidal forcing to ext-mode')
     endif 
     if(tidal_forcing_8) then
       call mpp_error(NOTE, &
         '==>Note from tidal_forcing_init: adding tidal forcing to ext-mode from 8-lunar/solar constituents')
     endif 
     if(tidal_forcing_m2 .and. tidal_forcing_8) then 
       call mpp_error(NOTE, &
        '==>Note from tidal_forcing_init: both tidal_forcing_m2 and tidal_forcing_8 = .true. Will use tidal_forcing_8.')
     endif 

     alphat = 0.948
     tidal_omega_K1 = 0.72921e-4*3600.*24.     ! K1-tide
     tidal_omega_O1 = 0.67598e-4*3600.*24.     ! O1-tide
     tidal_omega_P1 = 0.72523e-4*3600.*24.     ! P1-tide
     tidal_omega_Q1 = 0.64959e-4*3600.*24.     ! Q1-tide

     tidal_omega_M2 = 1.40519e-4*3600.*24.     ! M2-tide
     tidal_omega_S2 = 1.45444e-4*3600.*24.     ! S2-tide
     tidal_omega_N2 = 1.37880e-4*3600.*24.     ! N2-tide
     tidal_omega_K2 = 1.45842e-4*3600.*24.     ! K2-tide


     Love_M2 =   0.693                  ! (1+k-h), Love numbers k~=0.3 and h~=0.6 - 
     Love_S2 =   0.693                  ! (1+k-h), Love numbers k~=0.3 and h~=0.6 - 
     Love_N2 =   0.693                  ! (1+k-h), Love numbers k~=0.3 and h~=0.6 - 
     Love_K2 =   0.693                  ! (1+k-h), Love numbers k~=0.3 and h~=0.6 - 
     Love_K1 =   1.0 + 0.256 - 0.520    ! (1+k-h), Love numbers k~=0.3 and h~=0.6 - 
     Love_O1 =   1.0 + 0.298 - 0.603    ! (1+k-h), Love numbers k~=0.3 and h~=0.6 - 
     Love_P1 =   1.0 + 0.287 - 0.581    ! (1+k-h), Love numbers k~=0.3 and h~=0.6 - 
     Love_Q1 =   1.0 + 0.298 - 0.603    ! (1+k-h), Love numbers k~=0.3 and h~=0.6 - 

     amp_M2 =   0.242334 ! amplitude of M2 equilibrium tide (m)
     amp_S2 =   0.112743 ! amplitude of S2 equilibrium tide (m)
     amp_N2 =   0.046397 ! amplitude of N2 equilibrium tide (m)
     amp_K2 =   0.030684 ! amplitude of K2 equilibrium tide (m)
     amp_K1 =   0.141565 ! amplitude of K1 equilibrium tide (m)
     amp_O1 =   0.100661 ! amplitude of O1 equilibrium tide (m)
     amp_P1 =   0.046848 ! amplitude of P1 equilibrium tide (m)
     amp_Q1 =   0.019273 ! amplitude of Q1 equilibrium tide (m)

     rradian = 1./radian
     cos2lon(:,:) = cos(2.0*Grd%xt(:,:)*rradian)
     coslon (:,:) = cos(    Grd%xt(:,:)*rradian)
     coslat2(:,:) = cos(    Grd%xt(:,:)*rradian)**2
     sin2lon(:,:) = sin(2.0*Grd%xt(:,:)*rradian)
     sinlon (:,:) = sin(    Grd%xt(:,:)*rradian)
     sin2lat(:,:) = sin(2.0*Grd%xt(:,:)*rradian)

     id_eta_eq_tidal  = register_diag_field ('ocean_model', 'eta_eq_tidal', Grd%tracer_axes(1:2), &
                        Time%model_time, 'equilibrium tidal potential', 'meter',&
                        missing_value=-10.0, range=(/-1e3,1e3 /)) 
 
end subroutine tidal_forcing_init
! </SUBROUTINE> NAME="tidal_forcing_init"


!#######################################################################
! <SUBROUTINE NAME="get_tidal_forcing">
!
! <DESCRIPTION>
! Compute equilibrium tidal forcing.
! </DESCRIPTION>
!
subroutine get_tidal_forcing(Time, Ext_mode, dayr)

  type(ocean_time_type), intent(in)          :: Time 
  type(ocean_external_mode_type), intent(in) :: Ext_mode
  real, intent(in)                           :: dayr

  real    :: cosomegat_M2,sinomegat_M2
  real    :: cosomegat_S2,sinomegat_S2
  real    :: cosomegat_N2,sinomegat_N2
  real    :: cosomegat_K2,sinomegat_K2
  real    :: cosomegat_K1,sinomegat_K1
  real    :: cosomegat_O1,sinomegat_O1
  real    :: cosomegat_P1,sinomegat_P1
  real    :: cosomegat_Q1,sinomegat_Q1

  sinomegat_K1=sin(tidal_omega_K1*dayr)
  sinomegat_O1=sin(tidal_omega_O1*dayr)
  sinomegat_P1=sin(tidal_omega_P1*dayr)
  sinomegat_Q1=sin(tidal_omega_Q1*dayr)

  cosomegat_K1=cos(tidal_omega_K1*dayr)
  cosomegat_O1=cos(tidal_omega_O1*dayr)
  cosomegat_P1=cos(tidal_omega_P1*dayr)
  cosomegat_Q1=cos(tidal_omega_Q1*dayr)

  sinomegat_M2=sin(tidal_omega_M2*dayr)
  sinomegat_s2=sin(tidal_omega_S2*dayr)
  sinomegat_n2=sin(tidal_omega_N2*dayr)
  sinomegat_K2=sin(tidal_omega_K2*dayr)

  cosomegat_M2=cos(tidal_omega_M2*dayr)
  cosomegat_s2=cos(tidal_omega_S2*dayr)
  cosomegat_n2=cos(tidal_omega_N2*dayr)
  cosomegat_K2=cos(tidal_omega_K2*dayr)


! M2 constituent only
  if(tidal_forcing_m2) then 
    eta_eq_tidal(:,:)  =  Love_M2*amp_M2*(coslat2(:,:))*(cosomegat_M2*cos2lon(:,:) - sinomegat_M2*sin2lon(:,:))
  endif 

! 8 principle constituents
  if(tidal_forcing_8) then 
     eta_eq_tidal(:,:) =  Love_M2*amp_M2*(coslat2(:,:))*(cosomegat_M2*cos2lon(:,:) - sinomegat_M2*sin2lon(:,:))+&
                          Love_S2*amp_S2*(coslat2(:,:))*(cosomegat_S2*cos2lon(:,:) - sinomegat_S2*sin2lon(:,:))+&
                          Love_N2*amp_N2*(coslat2(:,:))*(cosomegat_N2*cos2lon(:,:) - sinomegat_N2*sin2lon(:,:))+&
                          Love_K2*amp_K2*(coslat2(:,:))*(cosomegat_K2*cos2lon(:,:) - sinomegat_K2*sin2lon(:,:))+&
                          Love_K1*amp_K1*(sin2lat(:,:))*(cosomegat_K1*coslon (:,:) - sinomegat_K1*sinlon (:,:))+&
                          Love_O1*amp_O1*(sin2lat(:,:))*(cosomegat_O1*coslon (:,:) - sinomegat_O1*sinlon (:,:))+&
                          Love_P1*amp_P1*(sin2lat(:,:))*(cosomegat_P1*coslon (:,:) - sinomegat_P1*sinlon (:,:))+&
                          Love_Q1*amp_Q1*(sin2lat(:,:))*(cosomegat_Q1*coslon (:,:) - sinomegat_Q1*sinlon (:,:))
  endif 

  if (id_eta_eq_tidal > 0) then 
    used = send_data (id_eta_eq_tidal, eta_eq_tidal(:,:), &
           Time%model_time, rmask=Grd%tmask(:,:,1), &
           is_in=isc, js_in=jsc, ie_in=iec, je_in=jec) 
  endif 

end subroutine get_tidal_forcing
! </SUBROUTINE> NAME="get_tidal_forcing"



end module ocean_freesurf_mod


