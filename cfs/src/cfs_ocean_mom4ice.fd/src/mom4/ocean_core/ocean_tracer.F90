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
module ocean_tracer_mod
!
!<CONTACT EMAIL="Stephen.Griffies@noaa.gov"> Stephen Griffies 
!</CONTACT>
!
!<CONTACT EMAIL="Matthew.Harrison@noaa.gov"> Matt Harrison 
!</CONTACT>
!
!<CONTACT EMAIL="Richard.Slater@noaa.gov"> Richard D. Slater (initialization)
!</CONTACT>
!
!<CONTACT EMAIL="Ronald.Pacanowski@noaa.gov"> Ron Pacanowski
!</CONTACT>
!
!<OVERVIEW>
! This module time steps the tracer fields.
!</OVERVIEW>
!
!<DESCRIPTION>
! This module time steps the tracer fields.
! Initialization for the tracer packages is done as well. 
!</DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! R.C. Pacanowski and S.M. Griffies
! The MOM3 Manual (1999)
! </REFERENCE>
!
! <REFERENCE>
! S.M. Griffies, M.J. Harrison, R.C. Pacanowski, and A. Rosati 
! A Technical Guide to MOM4 (2004)
! </REFERENCE>
!
! <REFERENCE>
! S.M. Griffies, Fundamentals of ocean climate models (2004)
! </REFERENCE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_tracer_nml">
!
!  <DATA NAME="zero_tendency" TYPE="logical">
!  If true, then will freeze the tracer fields.
!  </DATA> 
!
!  <DATA NAME="convective_adjust_on" TYPE="logical">
!  To use the Rahmstorf full convection scheme. 
!  </DATA> 
!
!  <DATA NAME="use_frazil" TYPE="logical">
!  If true, then will compute ocean heating in the top model grid cell
!  due to the production of frazil ice. 
!  </DATA> 
!
!  <DATA NAME="t_min" UNITS="deg C" TYPE="real">
!  Minimum potential temperature below which we gracefully bring down the model. 
!  </DATA> 
!  <DATA NAME="t_max" UNITS="deg C" TYPE="real">
!  Maximum potential temperature above which we gracefully bring down the model.  
!  </DATA> 
!  <DATA NAME="s_min" UNITS="psu" TYPE="real">
!  Minimum salinity below which we gracefully bring down the model. 
!  </DATA> 
!  <DATA NAME="s_max" UNITS="psu" TYPE="real">
!  Maximum salinity below which we gracefully bring down the model. 
!  </DATA> 
!
!  <DATA NAME="t_min_limit" UNITS="deg C" TYPE="real">
!  Minimum potential temperature below which will employ upwind advection 
!  instead of quicker, and horizontal diffusion instead of neutral physics. 
!  </DATA> 
!  <DATA NAME="t_max_limit" UNITS="deg C" TYPE="real">
!  Maximum potential temperature above which will employ upwind advection 
!  instead of quicker, and horizontal diffusion instead of neutral physics. 
!  </DATA> 
!  <DATA NAME="s_min_limit" UNITS="psu" TYPE="real">
!  Minimum salinity below which will employ upwind advection instead
!  of quicker, and horizontal diffusion instead of neutral physics.  
!  </DATA> 
!  <DATA NAME="s_max_limit" UNITS="psu" TYPE="real">
!  Maximum salinity below which will employ upwind advection instead
!  of quicker, and horizontal diffusion instead of neutral physics.  
!  </DATA> 
!
!  <DATA NAME="debug_tracer" TYPE="logical">
!  For debugging the tracer module
!  </DATA> 
!  <DATA NAME="ocean_tpm_debug" TYPE="logical">
!  For debugging ocean tracer package manager.  
!  </DATA> 
!</NAMELIST>
!
use constants_mod,     only: cp_ocean, epsln, rho_cp, kelvin, rho0, rho0r
use diag_manager_mod,  only: register_diag_field, send_data
use field_manager_mod, only: fm_string_len, fm_type_name_len
use field_manager_mod, only: fm_get_length, fm_get_type, fm_get_value
use field_manager_mod, only: fm_dump_list, fm_new_list, fm_change_list, fm_loop_over_list
use fms_mod,           only: read_data, write_data, open_namelist_file, check_nml_error, close_file
use fms_mod,           only: FATAL, WARNING, NOTE, stdout, stdlog
use fms_io_mod,        only: field_size
use mpp_domains_mod,   only: mpp_update_domains
use mpp_io_mod,        only: mpp_open, fieldtype
use mpp_io_mod,        only: MPP_NETCDF, MPP_OVERWR, MPP_ASCII, MPP_RDONLY, MPP_SINGLE, MPP_MULTI 
use mpp_mod,           only: mpp_error, mpp_pe, mpp_root_pe, mpp_broadcast, ALL_PES, mpp_max, mpp_min, mpp_chksum 
use mpp_mod,           only: mpp_clock_id, mpp_clock_begin, mpp_clock_end
use mpp_mod,           only: CLOCK_COMPONENT, CLOCK_SUBCOMPONENT, CLOCK_MODULE
use platform_mod,      only: i8_kind
use time_manager_mod,  only: time_type, set_time, increment_time, get_date, operator( + )

use ocean_convect_mod,          only: convection
use ocean_domains_mod,          only: get_local_indices
use ocean_horz_diffuse_mod,     only: horz_diffuse 
use ocean_neutral_physics_mod,  only: ocean_neutral_physics_rstrt, ocean_neutral_physics_end
use ocean_obc_mod,              only: ocean_obc_tracer, ocean_obc_update_boundary, ocean_obc_tracer_init
use ocean_tpm_mod,              only: ocean_tpm_init
use ocean_tpm_util_mod,         only: otpm_set_tracer_package, otpm_set_prog_tracer
use ocean_tpm_util_mod,         only: otpm_set_diag_tracer, otpm_get_string_array
use ocean_tpm_util_mod,         only: otpm_check_for_bad_fields, otpm_set_value
use ocean_tracer_advect_mod,    only: horz_advect_tracer, vert_advect_tracer
use ocean_tracer_util_mod,      only: tracer_prog_chksum, tracer_diag_chksum, tracer_min_max
use ocean_types_mod,            only: ocean_grid_type, ocean_domain_type, ocean_thickness_type
use ocean_types_mod,            only: ocean_time_type, ocean_time_steps_type
use ocean_types_mod,            only: ocean_prog_tracer_type, ocean_diag_tracer_type
use ocean_types_mod,            only: ocean_adv_vel_type, ocean_external_mode_type, ocean_density_type
use ocean_types_mod,            only: ADVECT_UPWIND, ADVECT_2ND_ORDER, ADVECT_4TH_ORDER, ADVECT_6TH_ORDER
use ocean_types_mod,            only: ADVECT_QUICKER, ADVECT_QUICKMOM3, ADVECT_MDFL_SUP_B, ADVECT_MDFL_SWEBY
use ocean_types_mod,            only: missing_value
use ocean_types_mod,            only: TWO_LEVEL, THREE_LEVEL
use ocean_util_mod,             only: invtri
use ocean_vert_mix_mod,         only: vert_diffuse 
use ocean_workspace_mod,        only: wrk1, wrk2


implicit none

private 

logical :: prog_module_initialized = .false.
logical :: diag_module_initialized = .false.

character(len=256) :: version='CVS $Id$'
character(len=256) :: tagname='Tag $Name$'
character(len=48), parameter          :: mod_name = 'ocean_tracer_mod'

integer :: num_tracers       =0
integer :: num_prog_tracers  =0
integer :: num_diag_tracers  =0
integer :: num_family_tracers=0

#include <ocean_memory.h>

type(ocean_grid_type), pointer   :: Grd =>NULL()
type(ocean_domain_type), pointer :: Dom =>NULL()

! for time steps and implicit vertical mixing 
integer :: tendency = 0
real    :: dtime    = 0.0
real    :: dtimer   = 0.0
real    :: dtts     = 0.0
real    :: dtuv     = 0.0
real    :: aidif    = 1.0

! for obc 
logical :: have_obc=.false.
 
integer :: index_temp  =-1
integer :: index_salt  =-1  
integer :: index_frazil=-1

! identification numbers for mpp clocks
integer :: id_clock_horz_diffuse
integer :: id_clock_vert_diffuse
integer :: id_clock_tracer_advect

! for diagnostics 
logical :: used
integer :: id_frazil_2d=-1
integer, allocatable, dimension(:) :: id_surface_smooth
integer, allocatable, dimension(:) :: id_prog
integer, allocatable, dimension(:) :: id_vdiff_impl
integer, allocatable, dimension(:) :: id_tendency
integer, allocatable, dimension(:) :: id_tendency_array
integer, allocatable, dimension(:) :: id_surf_tracer
integer, allocatable, dimension(:) :: id_diag
integer, allocatable, dimension(:) :: id_tmask_limit

! for ascii output
integer :: unit=6

public  update_ocean_tracer
public  ocean_prog_tracer_init
public  ocean_diag_tracer_init
public  ocean_tracer_rstrt
public  ocean_tracer_end
public  compute_tmask_limit

private update_advection_only

!---------------nml settings---------------

logical :: zero_tendency        = .false.
logical :: debug_tracer         = .false.
logical :: convective_adjust_on = .true.
logical :: use_frazil           = .false.

! for tracer package manager 
logical :: ocean_tpm_debug      = .false.

! set min/max temp and salinity range valid
! for equation of state. 
! if solution falls outside this range, model will be 
! brought down.
real    :: t_min=-5.0
real    :: t_max=40.0
real    :: s_min=-1.0
real    :: s_max=45.0

! min/max temp and salinity beyond which employ upwind 
! advection if set limit_tracer_range_quick=.true.
! (set in ocean_tracer_advect_nml) and/or 
! horizontal diffusion if limit_tracer_range_neutral=.true. 
! (set in ocean_neutral_physics_nml)
real    :: t_min_limit=-2.0
real    :: t_max_limit=32.0
real    :: s_min_limit=1.0
real    :: s_max_limit=42.0

namelist /ocean_tracer_nml/ zero_tendency, convective_adjust_on, &
                            t_min, t_max, s_min, s_max, use_frazil, debug_tracer, &
                            t_min_limit, t_max_limit, s_min_limit, s_max_limit,         &
                            ocean_tpm_debug

contains


!#######################################################################
! <FUNCTION NAME="ocean_prog_tracer_init">
!
! <DESCRIPTION>
! Initialization code for prognostic tracers, returning a pointer to 
! the T_prog array.
! </DESCRIPTION>
!
function ocean_prog_tracer_init (Grid, Thickness, Domain, Time, Time_steps, itemp, isalt, &
                                 num_prog, obc, debug, ens_ocean)                         &
                                 result (T_prog)  
  
  type(ocean_grid_type), intent(in), target   :: Grid
  type(ocean_thickness_type), intent(in)      :: Thickness
  type(ocean_domain_type), intent(in), target :: Domain
  type(ocean_time_type), intent(in)           :: Time
  type(ocean_time_steps_type), intent(in)     :: Time_steps
  integer, intent(out)                        :: itemp
  integer, intent(out)                        :: isalt
  integer, intent(out)                        :: num_prog
  logical, intent(in)                         :: obc
  logical, intent(in), optional               :: debug
  logical, intent(in), optional               :: ens_ocean

  ! return value 
  type(ocean_prog_tracer_type), dimension(:), pointer :: T_prog

  integer               :: i, j, k, m, n, kb
  integer               :: ierr, index, unit
  integer               :: tau, taum1, taup1
  integer               :: ndim, nvar, natt, ntime
  integer               :: ioun, io_status
  integer, dimension(4) :: siz 
  real                  :: fact
  character(len=32)     :: name, longname, units, control

  character(len=48),  parameter :: sub_name = 'ocean_prog_tracer_init'
  character(len=256), parameter :: error_header = '==>Error from ' // trim(mod_name) //   &
                                                  '(' // trim(sub_name) // '): '
  character(len=256), parameter :: warn_header = '==>Warning from ' // trim(mod_name) //  &
                                                 '(' // trim(sub_name) // '): '
  character(len=256), parameter :: note_header = '==>Note from ' // trim(mod_name) //     &
                                                 '(' // trim(sub_name) // '): '

  ! variables for tracer package
  real                                  :: conversion
  character(len=fm_string_len)          :: file_in
  character(len=fm_string_len)          :: file_out
  character(len=fm_string_len)          :: flux_units
  character(len=fm_type_name_len)       :: fm_type
  integer                               :: ind
  logical                               :: init_tracer
  real                                  :: max_range
  real                                  :: min_range
  real                                  :: max_flux_range
  real                                  :: min_flux_range
  real                                  :: max_tracer
  real                                  :: min_tracer
  real                                  :: min_tracer_limit
  real                                  :: max_tracer_limit
  character(len=fm_string_len)          :: name_in
  real                                  :: offset
  real, dimension(2)                    :: range_array
  real                                  :: scale_in
  real                                  :: additive_in
  logical                               :: use_only_advection
  logical                               :: const_init_tracer 
  real                                  :: const_init_value  
  character(len=fm_string_len)          :: string_fm
  character(len=fm_type_name_len)       :: typ
  character(len=fm_string_len), pointer, dimension(:) :: good_list


  if (prog_module_initialized) then
    call mpp_error(FATAL, trim(error_header) // ' Prognostic tracers already initialized')
  endif
  
  nullify(T_prog)

  write( stdlog(),'(/a/)') trim(version)

  have_obc = obc
  dtts     = Time_steps%dtts
  dtuv     = Time_steps%dtuv
  aidif    = Time_steps%aidif 

  ! provide for namelist over-ride
  ioun = open_namelist_file()
  read  (ioun, ocean_tracer_nml,iostat=io_status)
  write (stdout(),'(/)')
  write (stdout(), ocean_tracer_nml)  
  write (stdlog(), ocean_tracer_nml)
  ierr = check_nml_error(io_status,'ocean_tracer_nml')
  call close_file (ioun)

  if (PRESENT(debug) .and. .not. debug_tracer) then
    debug_tracer = debug
  endif 

  if(use_frazil) then 
    call mpp_error(NOTE, trim(note_header) // ' USING frazil heating in top model cell.')
  endif 
  if(zero_tendency) then 
    call mpp_error(NOTE, trim(note_header) // ' zero_tendency=true so will not time step tracer fields.')
  endif 
  if(convective_adjust_on) then 
    call mpp_error(NOTE, trim(note_header) // ' Using convective adjustment for vertically unstable water columns.')
  endif 

  ! some time step information 
  if (dtts /= dtuv) then
     write (stdout(),'(/a)') trim(warn_header) // ' Asynchronous timesteps (dtts > dtuv) imply inaccurate transients.'
     write (stdout(),'(a)')  '          and total tracer (i.e. heat content) is not conserved.'
     write (stdout(),'(a,f5.2,a/)') '            dtts =',dtts/dtuv,' times larger than dtuv.'         
  else
     call mpp_error(NOTE, trim(note_header) // ' Synchronous timesteps have been specified (dtts = dtuv).')
  endif

  tendency = Time_steps%tendency
  dtime    = Time_steps%dtime_t
  dtimer   = 1.0/(dtime+epsln)

  tau   = Time%tau
  taum1 = Time%taum1
  taup1 = Time%taup1

  ! initialize clock ids 
  id_clock_horz_diffuse  = mpp_clock_id('(Ocean tracer: horz diffuse) ',grain=CLOCK_MODULE)
  id_clock_vert_diffuse  = mpp_clock_id('(Ocean tracer: vert diffuse) ',grain=CLOCK_MODULE)
  id_clock_tracer_advect = mpp_clock_id('(Ocean tracer: advection)    ',grain=CLOCK_MODULE)

  ! perform requested debugging for ocean tracer package manager
  if (ocean_tpm_debug) then  
    write (stdout(),*)
    write (stdout(),*) 'Dumping field tree at start of ocean_prog_tracer_init'
    if (.not. fm_dump_list('/', recursive = .true.)) then  
      call mpp_error(FATAL, trim(error_header) // ' Problem dumping tracer tree')
    endif
  endif  


! initialize the tracer_packages, prog_tracers, diag_tracers, namelists, xland_mix and other GOOD lists
  
  ! make sure that /ocean_mod exists, just in case there were no inputs in the field table
  if (fm_new_list('/ocean_mod') .le. 0) then  
    call mpp_error(FATAL, trim(error_header) // ' Could not set the "ocean_mod" list')
  endif  

  if (fm_new_list('/ocean_mod/GOOD') .le. 0) then  
    call mpp_error(FATAL, trim(error_header) // ' Could not set the "GOOD" list')
  endif  
  call otpm_set_value('/ocean_mod/GOOD/good_ocean_mod_list', 'GOOD', append = .true.)

  if (fm_new_list('/ocean_mod/tracer_packages') .le. 0) then 
    call mpp_error(FATAL, trim(error_header) // ' Could not set "tracer packages" list')
  endif  
  call otpm_set_value('/ocean_mod/GOOD/good_ocean_mod_list', 'tracer_packages', append = .true.)

  if (fm_new_list('/ocean_mod/prog_tracers') .le. 0) then  
    call mpp_error(FATAL, trim(error_header) // ' Could not set the "prog_tracers" list')
  endif  
  call otpm_set_value('/ocean_mod/GOOD/good_ocean_mod_list', 'prog_tracers', append = .true.)

  if (fm_new_list('/ocean_mod/diag_tracers') .le. 0) then  
    call mpp_error(FATAL, trim(error_header) // ' Could not set the "diag_tracers" list')
  endif  
  call otpm_set_value('/ocean_mod/GOOD/good_ocean_mod_list', 'diag_tracers', append = .true.)

  if (fm_new_list('/ocean_mod/namelists') .le. 0) then  
    call mpp_error(FATAL, trim(error_header) // ' Could not set the "namelists" list')
  endif  
  call otpm_set_value('/ocean_mod/GOOD/good_ocean_mod_list', 'namelists', append = .true.)
    
  if (fm_new_list('/ocean_mod/xland_mix') .le. 0) then  
    call mpp_error(FATAL, trim(error_header) // ' Could not set the "xland_mix" list')
  endif  
  call otpm_set_value('/ocean_mod/GOOD/good_ocean_mod_list', 'xland_mix', append = .true.)

  if (fm_new_list('/ocean_mod/xland_insert') .le. 0) then  
    call mpp_error(FATAL, trim(error_header) // ' Could not set the "xland_insert" list')
  endif  
  call otpm_set_value('/ocean_mod/GOOD/good_ocean_mod_list', 'xland_insert', append = .true.)

  if (fm_new_list('/ocean_mod/diff_cbt_enhance') .le. 0) then  
    call mpp_error(FATAL, trim(error_header) // ' Could not set the "diff_cbt_enhance" list')
  endif  
  call otpm_set_value('/ocean_mod/GOOD/good_ocean_mod_list', 'diff_cbt_enhance', append = .true.)
  
  if (fm_new_list('/ocean_mod/spread_horz') .le. 0) then  
     call mpp_error(FATAL, trim(error_header) // ' Could not set the "spread_horz" list')
  endif
  call otpm_set_value('/ocean_mod/GOOD/good_ocean_mod_list', 'spread_horz', append = .true.)

  ! Initialize the good_namelists variable as it may not be otherwise set
  call otpm_set_value('/ocean_mod/GOOD/good_namelists', ' ', index = 0)


  ! Check for any errors in the number of fields in the ocean_mod list

  good_list => otpm_get_string_array('/ocean_mod/GOOD/good_ocean_mod_list',   &
       caller = trim(mod_name) // '(' // trim(sub_name) // ')')
  if (associated(good_list)) then  
    call otpm_check_for_bad_fields('/ocean_mod', good_list,       &
         caller = trim(mod_name) // '(' // trim(sub_name) // ')')
    deallocate(good_list)
  else  
    call mpp_error(FATAL,trim(error_header) // ' Empty "good_ocean_mod_list" list')
  endif 
 
  ! call routine to initialize the required tracer package
  if (otpm_set_tracer_package('required', caller=trim(mod_name)//'('//trim(sub_name)//')') .le. 0) then  
    call mpp_error(FATAL, trim(error_header) // ' Could not set required packages list')
  endif  

  ! initialize potential temperature and salinity
  index_temp = otpm_set_prog_tracer('temp', 'required',                                    &
       caller=trim(mod_name)//'('//trim(sub_name)//')',                                    &
       longname='Potential temperature', units='deg_C',                                    &
       conversion=rho_cp, offset=kelvin, min_tracer=t_min, max_tracer=t_max,               &
       min_range=-10.0, max_range=100.0, flux_units='Watts/m^2', min_flux_range=-1.0e+16,  &
       max_flux_range=1.0e+16, min_tracer_limit=t_min_limit, max_tracer_limit=t_max_limit, &
       file_in='INPUT/ocean_temp_salt.res.nc', file_out='RESTART/ocean_temp_salt.res.nc')

  index_salt = otpm_set_prog_tracer('salt', 'required',                                      &
       caller=trim(mod_name)//'('//trim(sub_name)//')',                                      &
       longname='Salinity', units='psu',                                                     &
       conversion=rho0*0.001, min_tracer=s_min, max_tracer=s_max,                            &
       min_range=-10.0, max_range=100.0, flux_units='kg/(sec*m^2)', min_flux_range=-1.0e+05, &
       max_flux_range=1.0e+05, min_tracer_limit=s_min_limit, max_tracer_limit=s_max_limit,   &
       file_in='INPUT/ocean_temp_salt.res.nc', file_out='RESTART/ocean_temp_salt.res.nc')


  ! initialize frazil as a diagnostic tracer (diagnosed from values of temp and salinity)
  ! allocate frazil even if use_frazil=.false--this alleviates problems with T_diag, as 
  ! otherwise it would have zero elements.  
!!$  if(use_frazil) then   
 index_frazil = otpm_set_diag_tracer('frazil',                                             &
     caller=trim(mod_name)//'('//trim(sub_name)//')',                                      &
     longname='frazil heating', units='J/m^2', name_in='frazil',                           &
     conversion=1.0, offset=0.0, min_tracer=0.0, max_tracer=1.e20,                         &
     min_range=-10.0, max_range=100.0, const_init_tracer=.true.,const_init_value=0.0,      &
     file_in='INPUT/ocean_frazil.res.nc',file_out='RESTART/ocean_frazil.res.nc')
!!$  endif 

  ! perform requested debugging for ocean tracer package manager
  if (ocean_tpm_debug) then  !{
    write (stdout(),*)
    write (stdout(),*) 'Dumping /ocean_mod field tree before ocean_tpm_init'
    if (.not. fm_dump_list('/ocean_mod', recursive = .true.)) then  !{
      call mpp_error(FATAL, trim(error_header) // ' Problem dumping tracer tree')
    endif  !}
  endif  !}


  ! call the initialization routine for the ocean tracer packages
  call ocean_tpm_init(Domain, Grid, Time, dtts) 

  ! perform requested debugging for ocean tracer package manager
  if (ocean_tpm_debug) then  
    write (stdout(),*)
    write (stdout(),*) 'Dumping /ocean_mod field tree after ocean_tpm_init'
    if (.not. fm_dump_list('/ocean_mod', recursive = .true.)) then  
      call mpp_error(FATAL, trim(error_header) // ' Problem dumping tracer tree')
    endif  
  endif  


  ! check for any errors in the number of fields in the tracer_packages list
  good_list => otpm_get_string_array('/ocean_mod/GOOD/good_tracer_packages',   &
       caller = trim(mod_name) // '(' // trim(sub_name) // ')')
  if (associated(good_list)) then  
    call otpm_check_for_bad_fields('/ocean_mod/tracer_packages', good_list,    &
         caller = trim(mod_name) // '(' // trim(sub_name) // ')')
    deallocate(good_list)
  else  
    call mpp_error(FATAL,trim(error_header) // ' Empty "good_tracer_packages" list')
  endif  

  ! check for any errors in the number of fields in the prog_tracers list
  good_list => otpm_get_string_array('/ocean_mod/GOOD/good_prog_tracers',         &
       caller = trim(mod_name) // '(' // trim(sub_name) // ')')
  if (associated(good_list)) then  
    call otpm_check_for_bad_fields('/ocean_mod/prog_tracers', good_list,          &
         caller = trim(mod_name) // '(' // trim(sub_name) // ')')
    deallocate(good_list)
  else  
    call mpp_error(FATAL,trim(error_header) // ' Empty "good_prog_tracers" list')
    endif  

  ! check for any errors in the number of fields in the namelists list
  good_list => otpm_get_string_array('/ocean_mod/GOOD/good_namelists',            &
       caller = trim(mod_name) // '(' // trim(sub_name) // ')')
  if (associated(good_list)) then  
    call otpm_check_for_bad_fields('/ocean_mod/namelists', good_list,             &
         caller = trim(mod_name) // '(' // trim(sub_name) // ')')
    deallocate(good_list)
  !else  
    !call mpp_error(FATAL,trim(error_header) // ' Empty "good_namelists" list')
  endif  

  ! get the number of tracers
  ! for now, we will not try to dynamically allocate the tracer arrays here
  num_prog_tracers = fm_get_length('/ocean_mod/prog_tracers')
  num_prog=num_prog_tracers

 
  ! allocate arrays based on the number of prognostic tracers
  allocate( T_prog            (num_prog_tracers) )
  allocate( id_surface_smooth (num_prog_tracers) )
  allocate( id_prog           (num_prog_tracers) )
  allocate( id_vdiff_impl     (num_prog_tracers) )
  allocate( id_tendency       (num_prog_tracers) )
  allocate( id_tendency_array (num_prog_tracers) )
  allocate( id_surf_tracer    (num_prog_tracers) )

  id_surface_smooth (:) = -1
  id_prog           (:) = -1
  id_vdiff_impl     (:) = -1
  id_tendency       (:) = -1
  id_tendency_array (:) = -1
  id_surf_tracer    (:) = -1

  ! set logical to determine when to mpp_update tracers
  do n=1,num_prog_tracers-1
    T_prog(n)%complete=.false.
  enddo
  T_prog(num_prog_tracers)%complete=.true.

  !  dump the lists for the tracer packages
  write (stdout(),*)
  write (stdout(),*) 'Dumping tracer_packages tracer tree'
  if (.not. fm_dump_list('/ocean_mod/tracer_packages', recursive = .true.)) then  
    call mpp_error(FATAL, trim(error_header) // ' Problem dumping tracer_packages tracer tree')
  endif  

  write (stdout(),*)
  write (stdout(),*) 'Dumping prog_tracers tracer tree'
  if (.not. fm_dump_list('/ocean_mod/prog_tracers', recursive = .true.)) then  
    call mpp_error(FATAL, trim(error_header) // ' Problem dumping prog_tracers tracer tree')
  endif  

  write (stdout(),*)
  write (stdout(),*) 'Dumping namelists tracer tree'
  if (.not. fm_dump_list('/ocean_mod/namelists', recursive = .true.)) then 
    call mpp_error(FATAL, trim(error_header) // ' Problem dumping namelists tracer tree')
  endif  

  !### finished with initializing t_prog arrays
  !### finished with call to tracer startup routine


  ! set local array indices

  Grd => Grid
  Dom => Domain

#ifndef STATIC_MEMORY
  call get_local_indices(Dom, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk=Grd%nk
#endif

  do n=1,num_prog_tracers
#ifndef STATIC_MEMORY
    allocate( T_prog(n)%field(isd:ied,jsd:jed,nk,3) ) 
    allocate( T_prog(n)%th_tendency(isd:ied,jsd:jed,nk)  )
    allocate( T_prog(n)%surface_smooth(isd:ied,jsd:jed)  )
    allocate( T_prog(n)%wrk1(isd:ied,jsd:jed,nk)    )
    allocate( T_prog(n)%tmask_limit(isd:ied,jsd:jed,nk) )
    allocate( T_prog(n)%K33_implicit(isd:ied,jsd:jed,nk))
#endif
    T_prog(n)%field(:,:,:,:)        = 0.0
    T_prog(n)%th_tendency(:,:,:)    = 0.0
    T_prog(n)%surface_smooth(:,:)   = 0.0
    T_prog(n)%wrk1(:,:,:)           = 0.0
    T_prog(n)%tmask_limit(:,:,:)    = 0.0
    T_prog(n)%K33_implicit(:,:,:)   = 0.0
    T_prog(n)%neutral_physics_limit = .false.
  enddo

  itemp = index_temp
  isalt = index_salt
  if (itemp == -1 .or. isalt == -1 ) then 
    call mpp_error(FATAL, trim(error_header) // ' temp or salt not included in table entries.  mom4 needs both.')
  endif 


  ! fill the field table entries
  n = 0
  do while (fm_loop_over_list('/ocean_mod/prog_tracers', name, typ, ind))  !{

     if (typ .ne. 'list') then  !{

         call mpp_error(FATAL, trim(error_header) // ' ' // trim(name) // ' is not a list')

     else  !}{

         n = n + 1  ! increment the array index

         if (n .ne. ind) then  !{
             write (stdout(),*) trim(warn_header), ' Tracer index, ', ind,   &
                  ' does not match array index, ', n, ' for ', trim(name)
         endif  !}

         ! save the name
         T_prog(n)%name = name

         if (.not. fm_change_list('/ocean_mod/prog_tracers/' // trim(name))) then  !{
             call mpp_error(FATAL, trim(error_header) // ' Problem changing to ' // trim(name))
         endif  !}

         ! save the units
         fm_type = fm_get_type('units')
         if (fm_type .eq. 'string') then  !{
             if (fm_get_value('units', units)) then  !{
                 T_prog(n)%units = units
             else  !}{
                 call mpp_error(FATAL, trim(error_header) // ' Problem getting "units" for ' // trim(name))
             endif  !}
         elseif (fm_type .eq. ' ') then  !}{
             T_prog(n)%units = ' '
             units = ' '
         else  !}{
             call mpp_error(FATAL, trim(error_header) // ' Wrong type "units" for ' // trim(name))
         endif  !}


         ! save the longname
         fm_type = fm_get_type('longname')
         if (fm_type .eq. 'string') then  !{
             if (fm_get_value('longname', longname)) then  !{
                 T_prog(n)%longname = longname
             else  !}{
                 call mpp_error(FATAL, trim(error_header) // ' Problem getting "longname" for ' // trim(name))
             endif  !}
         elseif (fm_type .eq. ' ') then  !}{
             T_prog(n)%longname = name
         else  !}{
             call mpp_error(FATAL, trim(error_header) // ' Wrong type "longname" for ' // trim(name))
         endif  !}

         ! save the conversion
         fm_type = fm_get_type('conversion')
         if (fm_type .eq. 'real') then  !{
             if (fm_get_value('conversion', conversion)) then  !{
                 T_prog(n)%conversion = conversion
             else  !}{
                 call mpp_error(FATAL, trim(error_header) // ' Problem getting "conversion" for ' // trim(name))
             endif  !}
         elseif (fm_type .eq. ' ') then  !}{
             T_prog(n)%conversion = 1.0
         else  !}{
             call mpp_error(FATAL, trim(error_header) // ' Wrong type "conversion" for ' // trim(name))
         endif  !}

         ! save the offset
         fm_type = fm_get_type('offset')
         if (fm_type .eq. 'real') then  !{
             if (fm_get_value('offset', offset)) then  !{
                 T_prog(n)%offset = offset
             else  !}{
                 call mpp_error(FATAL, trim(error_header) // ' Problem getting "offset" for ' // trim(name))
             endif  !}
         elseif (fm_type .eq. ' ') then  !}{
             T_prog(n)%offset = 0.0
         else  !}{
             call mpp_error(FATAL, trim(error_header) // ' Wrong type "offset" for ' // trim(name))
         endif  !}

         ! get the min and max of the tracer
         fm_type = fm_get_type('min_tracer')
         if (fm_type .eq. 'real') then  !{
             if (fm_get_value('min_tracer', min_tracer)) then  !{
                 fm_type = fm_get_type('max_tracer')
                 if (fm_type .eq. 'real') then  !{
                     if (.not. fm_get_value('max_tracer', max_tracer)) then  !{
                         call mpp_error(FATAL, trim(error_header) // ' Problem getting "max_tracer" for ' // trim(name))
                     endif  !}
                 elseif (fm_type .eq. ' ') then  !}{
                     call mpp_error(FATAL, trim(error_header) // ' "min_tracer" specified, but not "max_tracer" for ' // trim(name))
                 else  !}{
                     call mpp_error(FATAL, trim(error_header) // ' Wrong type "max_tracer" for ' // trim(name))
                 endif  !}
             else  !}{
                 call mpp_error(FATAL, trim(error_header) // ' Problem getting "min_tracer" for ' // trim(name))
             endif  !}
         elseif (fm_type .eq. ' ') then  !}{
             if (fm_get_type('max_tracer') .eq. ' ') then  !{
                 call mpp_error(NOTE, trim(note_header) // ' "min/max_tracer" not specified, using defaults for ' // trim(name))
                 min_tracer = -1.0e+20
                 max_tracer = 1.0e+20
             else  !}{
                 call mpp_error(FATAL, trim(error_header) // ' "max_tracer" specified, but not "min_tracer" for ' // trim(name))
             endif  !}
         else  !}{
             call mpp_error(FATAL, trim(error_header) // ' Wrong type "min_tracer" for ' // trim(name))
         endif  !}
         T_prog(n)%min_tracer = min_tracer
         T_prog(n)%max_tracer = max_tracer

         ! get the min and max of the range for analysis
         fm_type = fm_get_type('min_range')
         if (fm_type .eq. 'real') then  !{
             if (fm_get_value('min_range', min_range)) then  !{
                 fm_type = fm_get_type('max_range')
                 if (fm_type .eq. 'real') then  !{
                     if (.not. fm_get_value('max_range', max_range)) then  !{
                         call mpp_error(FATAL, trim(error_header) // ' Problem getting "max_range" for ' // trim(name))
                     endif  !}
                 elseif (fm_type .eq. ' ') then  !}{
                     call mpp_error(FATAL, trim(error_header) // ' "min_range" specified, but not "max_range" for ' // trim(name))
                 else  !}{
                     call mpp_error(FATAL, trim(error_header) // ' Wrong type "max_range" for ' // trim(name))
                 endif  !}
             else  !}{
                 call mpp_error(FATAL, trim(error_header) // ' Problem getting "min_range" for ' // trim(name))
             endif  !}
         elseif (fm_type .eq. ' ') then  !}{
             if (fm_get_type('max_range') .eq. ' ') then  !{
                 call mpp_error(NOTE, trim(note_header) // ' "min/max_range" not specified, using defaults for ' // trim(name))
                 min_range = 1.0
                 max_range = 0.0
             else  !}{
                 call mpp_error(FATAL, trim(error_header) // ' "max_range" specified, but not "min_range" for ' // trim(name))
             endif  !}
         else  !}{
             call mpp_error(FATAL, trim(error_header) // ' Wrong type "min_range" for ' // trim(name))
         endif  !}
         T_prog(n)%min_range = min_range
         T_prog(n)%max_range = max_range

         ! get the flux unit
         fm_type = fm_get_type('flux_units')
         if (fm_type .eq. 'string') then  !{
             if (fm_get_value('flux_units', flux_units)) then  !{
                 T_prog(n)%flux_units = flux_units
             else  !}{
                 call mpp_error(FATAL, trim(error_header) // ' Problem getting "flux_units" for ' // trim(name))
             endif  !}
         elseif (fm_type .eq. ' ') then  !}{
             T_prog(n)%flux_units = ' '
             flux_units = ' '
         else  !}{
             call mpp_error(FATAL, trim(error_header) // ' Wrong type "flux_units" for ' // trim(name))
         endif  !}

         ! get the min and max of the flux range for analysis
         fm_type = fm_get_type('min_flux_range')
         if (fm_type .eq. 'real') then  !{
             if (fm_get_value('min_flux_range', min_flux_range)) then  !{
                 fm_type = fm_get_type('max_flux_range')
                 if (fm_type .eq. 'real') then  !{
                     if (.not. fm_get_value('max_flux_range', max_flux_range)) then  !{
                         call mpp_error(FATAL, trim(error_header) //&
                                        ' Problem getting "max_flux_range" for ' // trim(name))
                     endif  !}
                 elseif (fm_type .eq. ' ') then  !}{
                     call mpp_error(FATAL, trim(error_header) // &
                                    ' "min_flux_range" specified, but not "max_flux_range" for ' // trim(name))
                 else  !}{
                     call mpp_error(FATAL, trim(error_header) // ' Wrong type "max_flux_range" for ' // trim(name))
                 endif  !}
             else  !}{
                 call mpp_error(FATAL, trim(error_header) // ' Problem getting "min_flux_range" for ' // trim(name))
             endif  !}
         elseif (fm_type .eq. ' ') then  !}{
             if (fm_get_type('max_flux_range') .eq. ' ') then  !{
                 call mpp_error(NOTE, trim(note_header) // &
                                ' "min/max_flux_range" not specified, using defaults for ' // trim(name))
                 min_flux_range = 1.0
                 max_flux_range = 0.0
             else  !}{
                 call mpp_error(FATAL, trim(error_header) //&
                                ' "max_flux_range" specified, but not "min_flux_range" for ' // trim(name))
             endif  !}
         else  !}{
             call mpp_error(FATAL, trim(error_header) // ' Wrong type "min_flux_range" for ' // trim(name))
         endif  !}
         T_prog(n)%min_flux_range = min_flux_range
         T_prog(n)%max_flux_range = max_flux_range

         ! save the input restart file
         fm_type = fm_get_type('file_in')
         if (fm_type .eq. 'string') then  !{
             if (fm_get_value('file_in', file_in)) then  !{
                 T_prog(n)%file_in = file_in
             else  !}{
                 call mpp_error(FATAL, trim(error_header) // ' Problem getting "file_in" for ' // trim(name))
             endif  !}
         else  !}{
             call mpp_error(FATAL, trim(error_header) // ' Wrong type or non-existent "file_in" for ' // trim(name))
         endif  !}

         ! save the input restart field name
         fm_type = fm_get_type('name_in')
         if (fm_type .eq. 'string') then  !{
             if (fm_get_value('name_in', name_in)) then  !{
                 T_prog(n)%name_in = name_in
             else  !}{
                 call mpp_error(FATAL, trim(error_header) // ' Problem getting "name_in" for ' // trim(name))
             endif  !}
         else  !}{
             call mpp_error(FATAL, trim(error_header) // ' Wrong type or non-existent "name_in" for ' // trim(name))
         endif  !}

         ! save the input scale factor
         fm_type = fm_get_type('scale_in')
         if (fm_type .eq. 'real') then  !{
             if (fm_get_value('scale_in', scale_in)) then  !{
                 T_prog(n)%scale_in = scale_in
             else  !}{
                 call mpp_error(FATAL, trim(error_header) // ' Problem getting "scale_in" for ' // trim(name))
             endif  !}
         else  !}{
             call mpp_error(FATAL, trim(error_header) // ' Wrong type or non-existent "scale_in" for ' // trim(name))
         endif  !}

         ! save the input additive factor
         fm_type = fm_get_type('additive_in')
         if (fm_type .eq. 'real') then  !{
             if (fm_get_value('additive_in', additive_in)) then  !{
                 T_prog(n)%additive_in = additive_in
             else  !}{
                 call mpp_error(FATAL, trim(error_header) // ' Problem getting "additive_in" for ' // trim(name))
             endif  !}
         elseif (fm_type .eq. ' ') then  !}{
             call mpp_error(NOTE, trim(note_header) // ' "additive_in" not specified, using default for ' // trim(name))
             T_prog(n)%additive_in = 0.0
         else  !}{
             call mpp_error(FATAL, trim(error_header) // ' Wrong type for "additive_in" for ' // trim(name))
         endif  !}

         ! save the output restart file
         fm_type = fm_get_type('file_out')
         if (fm_type .eq. 'string') then  !{
             if (fm_get_value('file_out', file_out)) then  !{
                 T_prog(n)%file_out = file_out
             else  !}{
                 call mpp_error(FATAL, trim(error_header) // ' Problem getting "file_out" for ' // trim(name))
             endif  !}
         else  !}{
             call mpp_error(FATAL, trim(error_header) // ' Wrong type or non-existent "file_out" for ' // trim(name))
         endif  !}

         ! save flag for whether to initialize this tracer
         fm_type = fm_get_type('init')
         if (fm_type .eq. 'logical') then  !{
             if (fm_get_value('init', init_tracer)) then  !{
                 T_prog(n)%init = init_tracer
             else  !}{
                 call mpp_error(FATAL, trim(error_header) // ' Problem getting "init" for ' // trim(name))
             endif  !}
         else  !}{
             call mpp_error(FATAL, trim(error_header) // ' Wrong type or non-existent "init" for ' // trim(name))
         endif  !}

         ! save flag for whether to globally initialize this tracer (optional)
         fm_type = fm_get_type('const_init_tracer')
         if (fm_type .eq. 'logical') then  !{
             if (fm_get_value('const_init_tracer', const_init_tracer)) then  !{
                 T_prog(n)%const_init_tracer = const_init_tracer
             else  !}{
                 call mpp_error(FATAL, trim(error_header) // ' Problem getting "const_init_tracer" for ' // trim(name))
             endif  !}
         elseif (fm_type .eq. ' ') then  !}{
             call mpp_error(NOTE, trim(note_header) // ' Using default (F) for "const_init_tracer" for ' // trim(name))
             T_prog(n)%const_init_tracer = .false.
         else  !}{
             call mpp_error(FATAL, trim(error_header) // ' Wrong type "const_init_tracer" for ' // trim(name))
         endif  !}

         ! save value to globally initialize this tracer (optional)
         fm_type = fm_get_type('const_init_value')
         if (fm_type .eq. 'real') then  !{
             if (fm_get_value('const_init_value', const_init_value)) then  !{
                 T_prog(n)%const_init_value = const_init_value
             else  !}{
                 call mpp_error(FATAL, trim(error_header) // ' Problem getting "const_init_value" for ' // trim(name))
             endif  !}
         elseif (fm_type .eq. ' ') then  !}{
             call mpp_error(NOTE, trim(note_header) // ' Using default (0) for "const_init_value" for ' // trim(name))
             T_prog(n)%const_init_value = 0.0
         else  !}{
             call mpp_error(FATAL, trim(error_header) // ' Wrong type "const_init_value" for ' // trim(name))
         endif  !}

         ! get the horizontal-advection-scheme
         fm_type = fm_get_type('horizontal-advection-scheme')
         if (fm_type .eq. 'string') then  !{
             if (.not. fm_get_value('horizontal-advection-scheme', string_fm)) then  !{
                 call mpp_error(FATAL, trim(error_header) // ' Bad horizontal-advection-scheme chosen for ' // trim(name))
             endif  !}
         elseif (fm_type .eq. ' ') then  !}{
             call mpp_error(FATAL, trim(error_header) // ' horz-advect-scheme unspecified for ' // trim(name))
         else  !}{
             call mpp_error(FATAL, trim(error_header) // ' Invalid horz-advect-scheme for ' // trim(name))
         endif  !}

         select case (trim(string_fm)) 
         case ('upwind')
             T_prog(n)%horz_advect_scheme = ADVECT_UPWIND
         case ('2nd_order')
             T_prog(n)%horz_advect_scheme = ADVECT_2ND_ORDER
         case ('4th_order')
             T_prog(n)%horz_advect_scheme = ADVECT_4TH_ORDER
         case ('quicker')
             T_prog(n)%horz_advect_scheme = ADVECT_QUICKER
         case ('quickMOM3')
             T_prog(n)%horz_advect_scheme = ADVECT_QUICKMOM3
         case ('6th_order')
             T_prog(n)%horz_advect_scheme = ADVECT_6TH_ORDER              
         case ('mdfl_sup_b')
             T_prog(n)%horz_advect_scheme = ADVECT_MDFL_SUP_B
         case ('mdfl_sweby')
             T_prog(n)%horz_advect_scheme = ADVECT_MDFL_SWEBY
         case default
             call mpp_error(FATAL, trim(error_header) // ' Invalid horz-advect-scheme for ' // trim(name))
         end select

         ! get the vertical-advection-scheme
         fm_type = fm_get_type('vertical-advection-scheme')
         if (fm_type .eq. 'string') then  !{
             if (.not. fm_get_value('vertical-advection-scheme', string_fm)) then  !{
                 call mpp_error(FATAL, trim(error_header) // ' Bad vertical-advection-scheme for ' // trim(name))
             endif  !}
         elseif (fm_type .eq. ' ') then  !}{
             call mpp_error(FATAL, trim(error_header) // ' Vertical advection scheme not specified for ' // trim(name))
         else  !}{
             call mpp_error(FATAL, trim(error_header) // ' Invalid vert-advect-scheme for ' // trim(name))
         endif  !}

         select case (trim(string_fm))    
         case ('upwind')
             T_prog(n)%vert_advect_scheme = ADVECT_UPWIND
         case ('2nd_order')
             T_prog(n)%vert_advect_scheme = ADVECT_2ND_ORDER
         case ('4th_order')
             T_prog(n)%vert_advect_scheme = ADVECT_4TH_ORDER
         case ('6th_order')
             T_prog(n)%vert_advect_scheme = ADVECT_6TH_ORDER              
         case ('quicker')
             T_prog(n)%vert_advect_scheme = ADVECT_QUICKER
         case ('quickMOM3')
             T_prog(n)%vert_advect_scheme = ADVECT_QUICKMOM3
         case ('mdfl_sup_b')
             T_prog(n)%vert_advect_scheme = ADVECT_MDFL_SUP_B
         case ('mdfl_sweby')
             T_prog(n)%vert_advect_scheme = ADVECT_MDFL_SWEBY
         case default
             call mpp_error(FATAL, trim(error_header) // ' Invalid vert-advect scheme for ' // trim(name))
         end select

         ! save the max_tracer_limit
         fm_type = fm_get_type('max_tracer_limit')
         if (fm_type .eq. 'real') then  !{
             if (fm_get_value('max_tracer_limit', max_tracer_limit)) then  !{
                 T_prog(n)%max_tracer_limit = max_tracer_limit
             else  !}{
                 call mpp_error(FATAL, trim(error_header) // ' Problem getting "max_tracer_limit" for ' // trim(name))
             endif  !}
         elseif (fm_type .eq. ' ') then  !}{
             T_prog(n)%max_tracer_limit = +1.0e+20
         else  !}{
             call mpp_error(FATAL, trim(error_header) // ' Wrong type "max_tracer_limit" for ' // trim(name))
         endif  !}

         ! save the min_tracer_limit
         fm_type = fm_get_type('min_tracer_limit')
         if (fm_type .eq. 'real') then  !{
             if (fm_get_value('min_tracer_limit', min_tracer_limit)) then  !{
                 T_prog(n)%min_tracer_limit = min_tracer_limit
             else  !}{
                 call mpp_error(FATAL, trim(error_header) // ' Problem getting "min_tracer_limit" for ' // trim(name))
             endif  !}
         elseif (fm_type .eq. ' ') then  !}{
             T_prog(n)%min_tracer_limit = -1.0e+20
         else  !}{
             call mpp_error(FATAL, trim(error_header) // ' Wrong type "min_tracer_limit" for ' // trim(name))
         endif  !}

         ! save flag for whether to transport tracer using only advection 
         fm_type = fm_get_type('use_only_advection')
         if (fm_type .eq. 'logical') then  !{
             if (fm_get_value('use_only_advection', use_only_advection)) then  !{
                 T_prog(n)%use_only_advection = use_only_advection
             else !}{
                 call mpp_error(FATAL, trim(error_header) // ' Problem getting "use_only_advection" for ' // trim(name))
             endif !} 
         elseif (fm_type .eq. ' ') then  !}{
             use_only_advection=.false.
             T_prog(n)%use_only_advection = use_only_advection
         else !}{
             call mpp_error(FATAL, trim(error_header) // ' Wrong type or non-existent "use_only_advection" for ' // trim(name))
         endif  !}
         if(T_prog(n)%use_only_advection) then
             call mpp_error(NOTE, trim(note_header) // ' Will evolve '// trim(T_prog(n)%name) // ' via advection alone.') 
         endif

         ! check consistency in selection of three-dimensional advection schemes 
         if(T_prog(n)%horz_advect_scheme==ADVECT_MDFL_SUP_B .and. T_prog(n)%vert_advect_scheme /= ADVECT_MDFL_SUP_B) then 
             call mpp_error(FATAL,&
                  trim(error_header) // ' mdfl_sup_b advection for ' // trim(name) // 'must be for **BOTH** horz & vert') 
         endif
         if(T_prog(n)%horz_advect_scheme==ADVECT_MDFL_SWEBY .and. T_prog(n)%vert_advect_scheme /= ADVECT_MDFL_SWEBY) then 
             call mpp_error(FATAL,&
                  trim(error_header) // ' mdfl_sweby advection for ' // trim(name) // 'must be for **BOTH** horz & vert') 
         endif
         if(T_prog(n)%vert_advect_scheme==ADVECT_MDFL_SUP_B .and. T_prog(n)%horz_advect_scheme /= ADVECT_MDFL_SUP_B) then 
             call mpp_error(FATAL,&
                  trim(error_header) // ' mdfl_sup_b advection for ' // trim(name) // 'must be for **BOTH** horz & vert') 
         endif
         if(T_prog(n)%vert_advect_scheme==ADVECT_MDFL_SWEBY .and. T_prog(n)%horz_advect_scheme /= ADVECT_MDFL_SWEBY) then 
             call mpp_error(FATAL,&
                  trim(error_header) // ' mdfl_sweby advection for ' // trim(name) // 'must be for **BOTH** horz & vert') 
         endif


     endif  !}
  enddo  !}


  ! register diagnostics 
  id_frazil_2d = register_diag_field ('ocean_model', 'frazil_2d', Grd%tracer_axes(1:2), &
         Time%model_time, 'ocn frazil heat flux over time step', 'W/m^2',&
         missing_value=missing_value, range=(/-1.e10,1.e10/))  

  allocate(id_tmask_limit(num_prog_tracers))

  do n = 1, num_prog_tracers  !{

   id_tmask_limit(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_tmask_limit', &
                       Grd%tracer_axes(1:3), Time%model_time, &
                       'use upwind (not quicker) and horz diff (not neutral)', &
                       'none', missing_value=-1.e2, range=(/-2.0,2.0/))

    if (T_prog(n)%max_range .gt. T_prog(n)%min_range) then  
      range_array(1) = T_prog(n)%min_range
      range_array(2) = T_prog(n)%max_range
      id_prog(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name), &
           Grd%tracer_axes(1:3),                                             &
           Time%model_time, trim(T_prog(n)%longname), trim(T_prog(n)%units), &
           missing_value=missing_value, range=range_array)
      id_surf_tracer(n) = register_diag_field ('ocean_model',                &
           'surface_'//trim(T_prog(n)%name),                                 &
           Grd%tracer_axes(1:2),                                             &
           Time%model_time, trim(T_prog(n)%longname), trim(T_prog(n)%units), &
           missing_value=missing_value, range=range_array)
    else  
      id_prog(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name), &
           Grd%tracer_axes(1:3),                                             &
           Time%model_time, trim(T_prog(n)%longname), trim(T_prog(n)%units), &
           missing_value=missing_value)
      id_surf_tracer(n) = register_diag_field ('ocean_model',                &
           'surface_'//trim(T_prog(n)%name),                                 &
           Grd%tracer_axes(1:2),                                             &
           Time%model_time, trim(T_prog(n)%longname), trim(T_prog(n)%units), &
           missing_value=missing_value)
    endif  
 

    ! register diagnostics for some tracer tendency terms

    if (T_prog(n)%max_flux_range .gt. T_prog(n)%min_flux_range) then  
      range_array(1) = T_prog(n)%min_flux_range
      range_array(2) = T_prog(n)%max_flux_range
      id_vdiff_impl(n) = register_diag_field ('ocean_model',                          &
           trim(T_prog(n)%name)//'_vdiff_impl', Grd%tracer_axes(1:3),               &
           Time%model_time, 'implicit vertical diffusion of '                         & 
           // trim(T_prog(n)%longname), trim(T_prog(n)%flux_units),                   &
           missing_value=missing_value, range=range_array)
      id_tendency(n) = register_diag_field ('ocean_model',                            & 
           trim(T_prog(n)%name)//'_totflux', Grd%tracer_axes(1:3),                  &
           Time%model_time, 'total tracer boundary flux', trim(T_prog(n)%flux_units), &
           missing_value=missing_value, range=range_array)
      id_tendency_array(n) = register_diag_field ('ocean_model',                      &
           trim(T_prog(n)%name)//'_th_tendency', Grd%tracer_axes(1:3),              &
           Time%model_time, 'th_tendency for ' // trim(T_prog(n)%longname),           &
           trim(T_prog(n)%flux_units),                                                &
           missing_value=missing_value, range=range_array)    
      id_surface_smooth(n) = register_diag_field ('ocean_model',                      &
           trim(T_prog(n)%name)//'_surface_smooth', Grd%tracer_axes(1:2),           &
           Time%model_time, 'surface smoother for ' // trim(T_prog(n)%name),          &
           trim(T_prog(n)%flux_units),                                                &
           missing_value=missing_value, range=range_array)    
    else  
      id_vdiff_impl(n) = register_diag_field ('ocean_model',                          &
           trim(T_prog(n)%name)//'_vdiff_impl', Grd%tracer_axes(1:3),               &
           Time%model_time, 'implicit vertical diffusion of ' //                      &
           trim(T_prog(n)%longname), trim(T_prog(n)%flux_units),                      &
           missing_value=missing_value)
      id_tendency(n) = register_diag_field ('ocean_model',                            &
           trim(T_prog(n)%name)//'_totflux', Grd%tracer_axes(1:3),                  &
           Time%model_time, 'total tracer boundary flux', trim(T_prog(n)%flux_units), &
           missing_value=missing_value)
      id_tendency_array(n) = register_diag_field ('ocean_model',                      &
           trim(T_prog(n)%name)//'_th_tendency', Grd%tracer_axes(1:3),              &
           Time%model_time, 'th_tendency for ' // trim(T_prog(n)%longname),           &
           trim(T_prog(n)%flux_units),                                                &
           missing_value=missing_value)    
      id_surface_smooth(n) = register_diag_field ('ocean_model',                      &
           trim(T_prog(n)%name)//'_surface_smooth', Grd%tracer_axes(1:2),           &
           Time%model_time, 'surface smoother for ' // trim(T_prog(n)%name),          &
           trim(T_prog(n)%flux_units),                                                &
           missing_value=missing_value)    
    endif  

  enddo  !} n-loop 

  ! read prognostic tracer intial conditions or restarts 
  write (stdout(),*) trim(note_header), ' Reading prognostic tracer initial conditions or restarts'

  do n=1,num_prog_tracers  !{

    write (stdout(),*)
    write (stdout(),*) 'Read occurring for tracer number', n,' at time level tau. This tracer is called ',trim(T_prog(n)%name)

    if(.not. (Time%init .or. T_prog(n)%init)) then 
        if(tendency==TWO_LEVEL) then 
            write (stdout(),'(/a)')'Expecting only one time record from the tracer restart.' 
        elseif(tendency==THREE_LEVEL) then  
            write (stdout(),'(/a)')'Expecting two time records from the tracer restart.' 
        endif
    endif

    T_prog(n)%field(:,:,:,taum1) = 0.0

    write (stdout(),*) trim(T_prog(n)%name)
    write (stdout(),*) trim(T_prog(n)%file_in)
    write (stdout(),*) trim(T_prog(n)%name_in)
    write (stdout(),*) T_prog(n)%scale_in
    write (stdout(),*) T_prog(n)%additive_in

    if( (Time%init .or. T_prog(n)%init) .and. T_prog(n)%const_init_tracer) then 
       write (stdout(),*) 'Initializing the tracer ',trim(T_prog(n)%name),' to the constant ',T_prog(n)%const_init_value 
       T_prog(n)%field(:,:,:,:) = T_prog(n)%const_init_value 
    else

       call field_size(T_prog(n)%file_in,T_prog(n)%name_in, siz)
       if (tendency==TWO_LEVEL .and. siz(4) > 1) then
          write(stdout(),'(/a)') '==>ERROR: Attempt to read restart for tracer from a 3-level time scheme (2 time records)'
          write(stdout(),'(a)')  '          when running mom4 with 2-level timestepping (only need 1 time record in restart).'
          write(stdout(),'(a)')  '          Reduce restart file to a single time record in order to avoid confusion.'
          call mpp_error(FATAL, &
                       'Reading 3-time level ocean tracer (w/ 2 time records) while using 2-level (needs only 1 record)')
       endif

       call read_data(T_prog(n)%file_in, T_prog(n)%name_in, T_prog(n)%field(:,:,:,tau), Domain%domain2d, timelevel=1)
       write (stdout(),*) 'Completed initialization of tracer ',trim(T_prog(n)%name),' at time level tau'
    endif 

    ! modify initial conditions 
    if (Time%init .or. T_prog(n)%init) then  

        ! table driven scaling and addition
        if (T_prog(n)%scale_in .ne. 1.0) then  
            write (stdout(),*) 'Scaling the initial condition for ', trim(T_prog(n)%name),' by the value ',T_prog(n)%scale_in
            T_prog(n)%field(:,:,:,tau) = T_prog(n)%field(:,:,:,tau) * T_prog(n)%scale_in
        endif
        if (T_prog(n)%additive_in .ne. 0.0) then  
            write (stdout(),*) 'Adding the constant ',T_prog(n)%additive_in,' to initial condition for ',trim(T_prog(n)%name)
            T_prog(n)%field(:,:,:,tau) = T_prog(n)%field(:,:,:,tau) + T_prog(n)%additive_in
        endif

        ! modifications due to partial bottom cell
        write (stdout(),*) 'After reading ic, linearly interpolate ',trim(T_prog(n)%name),' to partial cell bottom.'
        do j=jsc,jec
           do i=isc,iec
              kb = Grd%kmt(i,j)
              if (kb .gt. 1) then
                  fact = Thickness%dhwt(i,j,kb-1)/Grd%dzw(kb-1)
                  T_prog(n)%field(i,j,kb,tau) =                  &
                       T_prog(n)%field(i,j,kb-1,tau) -           &
                       fact*(T_prog(n)%field(i,j,kb-1,tau) -     &
                       T_prog(n)%field(i,j,kb,tau))
              endif

              ! zero out tracers in land points
              do k = kb+1,nk
                 T_prog(n)%field(i,j,k,tau) = 0.0
              enddo

           enddo
        enddo

    endif

    ! fill halos for tau tracer value 
    call mpp_update_domains(T_prog(n)%field(:,:,:,tau), Dom%domain2d)

    ! initialize tracer at taup1 to tracer at tau 
    T_prog(n)%field(:,:,:,taup1) = T_prog(n)%field(:,:,:,tau) 

    ! read in second time level when running with threelevel scheme 
    if(tendency==THREE_LEVEL .and. .not.Time%init .and. .not.T_prog(n)%init) then 
      write (stdout(),*) 'Since Time%init=.false. and using threelevel tendency, read taup1 for ',trim(T_prog(n)%name)
      call read_data(T_prog(n)%file_in, T_prog(n)%name_in, T_prog(n)%field(:,:,:,taup1), Domain%domain2d, timelevel=2)
      write (stdout(),*) 'Completed read of restart for ', trim(T_prog(n)%name),' at time taup1'
      call mpp_update_domains(T_prog(n)%field(:,:,:,taup1), Dom%domain2d)
    endif 

  enddo  !} n-loop 

  write (stdout(),*) trim(note_header), ' finished reading prognostic tracer restarts.'

  if(have_obc) call ocean_obc_tracer_init(Time, T_prog, num_prog_tracers, debug_tracer)

  prog_module_initialized = .true.

  if (debug_tracer) then
    do n= 1, num_prog_tracers
      write(stdout(),*) 'Checksum from ocean_tracer_mod: initial tracer checksum (tau)   ===>'
      call tracer_prog_chksum(Time, T_prog(n), tau)
      write(stdout(),*) 'Checksum from ocean_tracer_mod: initial tracer checksum (taup1) ===>'
      call tracer_prog_chksum(Time, T_prog(n), taup1)
      call tracer_min_max(Time, T_prog(n))
    enddo
  endif

end function ocean_prog_tracer_init  
! </FUNCTION> NAME="ocean_prog_tracer_init">


!#######################################################################
! <FUNCTION NAME="ocean_diag_tracer_init">
!
! <DESCRIPTION>
! Initialization code for diagnostic tracers, returning a pointer to the T_diag array
! </DESCRIPTION>
!
function ocean_diag_tracer_init (Time, Thickness, num_diag, debug, ens_ocean)   &
     result (T_diag)  !{
  
  type(ocean_time_type), intent(in), target   :: Time
  type(ocean_thickness_type), intent(in)      :: Thickness
  integer, intent(out)                        :: num_diag
  logical, intent(in), optional               :: debug
  logical, intent(in), optional               :: ens_ocean

!
!       Return type
!

  type(ocean_diag_tracer_type), dimension(:), pointer :: T_diag

  integer :: kb, i, j, k, ierr, index, unit
  integer :: ndim, nvar, natt, ntime
  integer :: ioun, io_status, n, m
  character(len=32) :: name, longname, units, control
  real :: fact
  character(len=48), parameter  :: sub_name = 'ocean_diag_tracer_init'
  character(len=256), parameter :: error_header = '==>Error from ' // trim(mod_name) //   &
                                                  '(' // trim(sub_name) // '): '
  character(len=256), parameter :: warn_header = '==>Warning from ' // trim(mod_name) //  &
                                                 '(' // trim(sub_name) // '): '
  character(len=256), parameter :: note_header = '==>Note from ' // trim(mod_name) //     &
                                                 '(' // trim(sub_name) // '): '

  ! variables for tracer package
  real                                  :: conversion
  character(len=fm_string_len)          :: file_in
  character(len=fm_string_len)          :: file_out
  character(len=fm_type_name_len)       :: fm_type
  integer                               :: ind
  integer                               :: index_diag_tracers
  logical                               :: init_tracer
  real                                  :: max_range
  real                                  :: min_range
  real                                  :: max_tracer
  real                                  :: min_tracer
  character(len=fm_string_len)          :: name_in
  real                                  :: offset
  real, dimension(2)                    :: range_array
  real                                  :: scale_in
  real                                  :: additive_in
  logical                               :: const_init_tracer 
  real                                  :: const_init_value  
  character(len=fm_string_len)          :: string_fm
  character(len=fm_type_name_len)       :: typ
  character(len=fm_string_len), pointer, dimension(:) :: good_list

  
  if (diag_module_initialized) then
    call mpp_error(FATAL, trim(error_header) // ' Diagnostic tracers already initialized')
  endif

  nullify(T_diag)
  
!
!       Check for any errors in the number of fields in the diag_tracers list
!

  good_list => otpm_get_string_array('/ocean_mod/GOOD/good_diag_tracers',         &
       caller = trim(mod_name) // '(' // trim(sub_name) // ')')
  if (associated(good_list)) then  !{
    call otpm_check_for_bad_fields('/ocean_mod/diag_tracers', good_list,             &
         caller = trim(mod_name) // '(' // trim(sub_name) // ')')
    deallocate(good_list)
  endif  !}


! Get the number of tracers
! For now, we will not try to dynamically allocate the tracer arrays here

  num_diag_tracers = fm_get_length('/ocean_mod/diag_tracers')

  num_diag = num_diag_tracers

!
!       Only do the following if there are diag tracers
!

  if (num_diag_tracers .gt. 0) then  !{

    write (stdout(),*)
    write (stdout(),*) trim(note_header), ' ', num_diag_tracers, ' diagnostic tracers requested.'

!
!       Allocate the T_diag array
!

    allocate( T_diag(num_diag_tracers) )
    allocate( id_diag(num_diag_tracers) )
    id_diag(:) = -1

!
!       Dump the diagnostic tracer tree
!

    write (stdout(),*)
    write (stdout(),*) 'Dumping ocean diag field tree after reading diag tracer tree'
    if (.not. fm_dump_list('/ocean_mod/diag_tracers', recursive = .true.)) then  
      call mpp_error(FATAL, trim(error_header) // ' Problem dumping diag_tracers tracer tree')
    endif  

  !### finished with initializing t_diag arrays
  !### finished with call to tracer startup routine


! allocate tracer type arrays and fill in tracer table information  

    if (num_diag_tracers .gt. 0 ) then 
      do n=1,num_diag_tracers 
#ifndef STATIC_MEMORY
        allocate(T_diag(n)%field(isd:ied,jsd:jed,nk))
#endif
        T_diag(n)%field(:,:,:) = 0.0
      enddo 
    endif 

!
!       process the T_diag array
!

    n = 0
    do while (fm_loop_over_list('/ocean_mod/diag_tracers', name, typ, ind))  !{
     
      if (typ .ne. 'list') then  !{
         
        call mpp_error(FATAL, trim(error_header) // ' ' // trim(name) // ' is not a list')
         
      else  !}{

        n = n + 1

        if (n .ne. ind) then  
          write (stdout(),*) trim(warn_header), ' Tracer index, ', ind,   &
               ' does not match array index, ', n, ' for ', trim(name)
        endif  

        ! save the name (required)
        T_diag(n)%name = name

        if (.not. fm_change_list('/ocean_mod/diag_tracers/' // trim(name))) then  
          call mpp_error(FATAL, trim(error_header) // ' Problem changing to ' // trim(name))
        endif  

        ! save the units (optional)
        fm_type = fm_get_type('units')
        if (fm_type .eq. 'string') then  !{
          if (fm_get_value('units', units)) then  !{
            T_diag(n)%units = units
          else  !}{
            call mpp_error(FATAL, trim(error_header) // ' Problem getting "units" for ' // trim(name))
          endif  !}
        elseif (fm_type .eq. ' ') then  !}{
          call mpp_error(NOTE, trim(note_header) // ' Using empty "units" for ' // trim(name))
          T_diag(n)%units = ' '
          units = ' '
        else  !}{
          call mpp_error(FATAL, trim(error_header) // ' Wrong type "units" for ' // trim(name))
        endif  !}


        ! save the longname (optional)
        fm_type = fm_get_type('longname')
        if (fm_type .eq. 'string') then  !{
          if (fm_get_value('longname', longname)) then  !{
            T_diag(n)%longname = longname
          else  !}{
            call mpp_error(FATAL, trim(error_header) // ' Problem getting "longname" for ' // trim(name))
          endif  !}
        elseif (fm_type .eq. ' ') then  !}{
          call mpp_error(NOTE, trim(note_header) // ' Using "name" as default for "longname" for ' // trim(name))
          T_diag(n)%longname = name
        else  !}{
          call mpp_error(FATAL, trim(error_header) // 'Wrong type "longname" for ' // trim(name))
        endif  !}

        ! save the conversion (optional)
        fm_type = fm_get_type('conversion')
        if (fm_type .eq. 'real') then  !{
          if (fm_get_value('conversion', conversion)) then  !{
            T_diag(n)%conversion = conversion
          else  !}{
            call mpp_error(FATAL, trim(error_header) // ' Problem getting "conversion" for ' // trim(name))
          endif  !}
        elseif (fm_type .eq. ' ') then  !}{
          call mpp_error(NOTE, trim(note_header) // ' Using 1.0 for "conversion" for ' // trim(name))
          T_diag(n)%conversion = 1.0
        else  !}{
         call mpp_error(FATAL, trim(error_header) // 'Wrong type "conversion" for ' // trim(name))
        endif  !}

        ! save the offset (optional)
        fm_type = fm_get_type('offset')
        if (fm_type .eq. 'real') then  !{
          if (fm_get_value('offset', offset)) then  !{
            T_diag(n)%offset = offset
          else  !}{
           call mpp_error(FATAL, trim(error_header) // ' Problem getting "offset" for ' // trim(name))
          endif  !}
        elseif (fm_type .eq. ' ') then  !}{
          call mpp_error(NOTE, trim(note_header) // ' Using 0.0 for "offset" for ' // trim(name))
          T_diag(n)%offset = 0.0
        else  !}{
          call mpp_error(FATAL, trim(error_header) // 'Wrong type "offset" for ' // trim(name))
        endif  !}

        ! get the min and max of the tracer for analysis (optional)
        fm_type = fm_get_type('min_tracer')
        if (fm_type .eq. 'real') then  !{
          if (fm_get_value('min_tracer', min_tracer)) then  !{
            fm_type = fm_get_type('max_tracer')
            if (fm_type .eq. 'real') then  !{
              if (.not. fm_get_value('max_tracer', max_tracer)) then  !{
                call mpp_error(FATAL, trim(error_header) // ' Problem getting "max_tracer" for ' // trim(name))
              endif  !}
            elseif (fm_type .eq. ' ') then  !}{
              call mpp_error(FATAL, trim(error_header) // '"min_tracer" specified, but not "max_tracer" for ' // trim(name))
            else  !}{
              call mpp_error(FATAL, trim(error_header) // 'Wrong type "max_tracer" for ' // trim(name))
            endif  !}
          else  !}{
            call mpp_error(FATAL, trim(error_header) // ' Problem getting "min_tracer" for ' // trim(name))
          endif  !}
        elseif (fm_type .eq. ' ') then  !}{
          if (fm_get_type('max_tracer') .eq. ' ') then  !{
            call mpp_error(NOTE, trim(note_header) // '"min/max_tracer" not specified, using defaults for ' // trim(name))
            min_tracer = -1.0e+20
            max_tracer = 1.0e+20
          else  !}{
            call mpp_error(FATAL, trim(error_header) // '"max_tracer" specified, but not "min_tracer" for ' // trim(name))
          endif  !}
        else  !}{
          call mpp_error(FATAL, trim(error_header) // 'Wrong type "min_tracer" for ' // trim(name))
        endif  !}
        T_diag(n)%min_tracer = min_tracer
        T_diag(n)%max_tracer = max_tracer

        ! get the min and max of the range for analysis (optional)
        fm_type = fm_get_type('min_range')
        if (fm_type .eq. 'real') then  !{
          if (fm_get_value('min_range', min_range)) then  !{
            fm_type = fm_get_type('max_range')
            if (fm_type .eq. 'real') then  !{
              if (.not. fm_get_value('max_range', max_range)) then  !{
                call mpp_error(FATAL, trim(error_header) // ' Problem getting "max_range" for ' // trim(name))
              endif  !}
            elseif (fm_type .eq. ' ') then  !}{
              call mpp_error(FATAL, trim(error_header) // '"min_range" specified, but not "max_range" for ' // trim(name))
            else  !}{
              call mpp_error(FATAL, trim(error_header) // 'Wrong type "max_range" for ' // trim(name))
            endif  !}
          else  !}{
            call mpp_error(FATAL, trim(error_header) // ' Problem getting "min_range" for ' // trim(name))
          endif  !}
        elseif (fm_type .eq. ' ') then  !}{
          if (fm_get_type('max_range') .eq. ' ') then  !{
            call mpp_error(NOTE, trim(note_header) // '"min/max_range" not specified, using defaults for ' // trim(name))
            min_range = 1.0
            max_range = 0.0
          else  !}{
            call mpp_error(FATAL, trim(error_header) // '"max_range" specified, but not "min_range" for ' // trim(name))
          endif  !}
        else  !}{
          call mpp_error(FATAL, trim(error_header) // 'Wrong type "min_range" for ' // trim(name))
        endif  !}
        T_diag(n)%min_range = min_range
        T_diag(n)%max_range = max_range

        ! save the input restart file (optional)
        fm_type = fm_get_type('file_in')
        if (fm_type .eq. 'string') then  !{
          if (fm_get_value('file_in', file_in)) then  !{
            T_diag(n)%file_in = file_in
          else  !}{
            call mpp_error(FATAL, trim(error_header) // ' Problem getting "file_in" for ' // trim(name))
          endif  !}
        elseif (fm_type .eq. ' ') then  !}{
          call mpp_error(NOTE, trim(note_header) // ' Using " " for "file_in" for ' // trim(name))
          T_diag(n)%file_in = ' '
        else  !}{
          call mpp_error(FATAL, trim(error_header) // 'Wrong type "file_in" for ' // trim(name))
        endif  !}

        ! save the input restart field name (optional)
        fm_type = fm_get_type('name_in')
        if (fm_type .eq. 'string') then  !{
          if (fm_get_value('name_in', name_in)) then  !{
            T_diag(n)%name_in = name_in
          else  !}{
            call mpp_error(FATAL, trim(error_header) // ' Problem getting "name_in" for ' // trim(name))
          endif  !}
        elseif (fm_type .eq. ' ') then  !}{
          call mpp_error(NOTE, trim(note_header) // ' Using name for "name_in" for ' // trim(name))
          T_diag(n)%name_in = name
        else  !}{
          call mpp_error(FATAL, trim(error_header) // 'Wrong type "name_in" for ' // trim(name))
        endif  !}

        ! save the input scale factor (optional)
        fm_type = fm_get_type('scale_in')
        if (fm_type .eq. 'real') then  !{
          if (fm_get_value('scale_in', scale_in)) then  !{
            T_diag(n)%scale_in = scale_in
          else  !}{
            call mpp_error(FATAL, trim(error_header) // ' Problem getting "scale_in" for ' // trim(name))
          endif  !}
        elseif (fm_type .eq. ' ') then  !}{
          call mpp_error(NOTE, trim(note_header) // ' Using 1.0 for "scale_in" for ' // trim(name))
          T_diag(n)%scale_in = 1.0
        else  !}{
          call mpp_error(FATAL, trim(error_header) // 'Wrong type "scale_in" for ' // trim(name))
        endif  !}

        ! save the input additive factor (optional)
        fm_type = fm_get_type('additive_in')
        if (fm_type .eq. 'real') then  !{
          if (fm_get_value('additive_in', additive_in)) then  !{
            T_diag(n)%additive_in = additive_in
          else  !}{
            call mpp_error(FATAL, trim(error_header) // ' Problem getting "additive_in" for ' // trim(name))
          endif  !}
        elseif (fm_type .eq. ' ') then  !}{
          call mpp_error(NOTE, trim(note_header) // ' Using 0.0 for "additive_in" for ' // trim(name))
          T_diag(n)%additive_in = 0.0
        else  !}{
          call mpp_error(FATAL, trim(error_header) // 'Wrong type "additive_in" for ' // trim(name))
        endif  !}

        ! save the output restart file (optional)
        fm_type = fm_get_type('file_out')
        if (fm_type .eq. 'string') then  !{
          if (fm_get_value('file_out', file_out)) then  !{
            T_diag(n)%file_out = file_out
          else  !}{
            call mpp_error(FATAL, trim(error_header) // ' Problem getting "file_out" for ' // trim(name))
          endif  !}
        elseif (fm_type .eq. ' ') then  !}{
          call mpp_error(NOTE, trim(note_header) // ' Using " " for "file_out" for ' // trim(name))
          T_diag(n)%file_out = ' '
        else  !}{
          call mpp_error(FATAL, trim(error_header) // ' Wrong type "file_out" for ' // trim(name))
        endif  !}

        ! save flag for whether to initialize this tracer (optional)
        fm_type = fm_get_type('init')
        if (fm_type .eq. 'logical') then  !{
          if (fm_get_value('init', init_tracer)) then  !{
            T_diag(n)%init = init_tracer
          else  !}{
            call mpp_error(FATAL, trim(error_header) // ' Problem getting "init" for ' // trim(name))
          endif  !}
        elseif (fm_type .eq. ' ') then  !}{
          call mpp_error(NOTE, trim(note_header) // ' Using false for "init" for ' // trim(name))
          T_diag(n)%init = .false.
        else  !}{
          call mpp_error(FATAL, trim(error_header) // ' Wrong type "init" for ' // trim(name))
        endif  !}

        ! save flag for whether to globally initialize this tracer (optional)
        fm_type = fm_get_type('const_init_tracer')
        if (fm_type .eq. 'logical') then  !{
          if (fm_get_value('const_init_tracer', const_init_tracer)) then  !{
            T_diag(n)%const_init_tracer = const_init_tracer
          else  !}{
            call mpp_error(FATAL, trim(error_header) // ' Problem getting "const_init_tracer" for ' // trim(name))
          endif  !}
        elseif (fm_type .eq. ' ') then  !}{
          call mpp_error(NOTE, trim(note_header) // ' Using false for "const_init_tracer" for ' // trim(name))
          T_diag(n)%const_init_tracer = .false.
        else  !}{
          call mpp_error(FATAL, trim(error_header) // ' Wrong type "const_init_tracer" for ' // trim(name))
        endif  !}

        ! save value to globally initialize this tracer (optional)
        fm_type = fm_get_type('const_init_value')
        if (fm_type .eq. 'real') then  !{
          if (fm_get_value('const_init_value', const_init_value)) then  !{
            T_diag(n)%const_init_value = const_init_value
          else  !}{
            call mpp_error(FATAL, trim(error_header) // ' Problem getting "const_init_value" for ' // trim(name))
          endif  !}
        elseif (fm_type .eq. ' ') then  !}{
          call mpp_error(NOTE, trim(note_header) // ' Using 0.0 for "const_init_value" for ' // trim(name))
          T_diag(n)%const_init_value = 0.0
        else  !}{
          call mpp_error(FATAL, trim(error_header) // ' Wrong type "const_init_value" for ' // trim(name))
        endif  !}

      endif  !}
    enddo  !}


    ! register diagnostics for the tracer
    do n = 1, num_diag_tracers  !{

      if (T_diag(n)%max_range .gt. T_diag(n)%min_range) then  !{
        range_array(1) = T_diag(n)%min_range
        range_array(2) = T_diag(n)%max_range
        id_diag(n) = register_diag_field ('ocean_model', trim(T_diag(n)%name), &
             Grd%tracer_axes(1:3),                                             & 
             Time%model_time, trim(T_diag(n)%longname), trim(T_diag(n)%units), &
             missing_value=missing_value, range=range_array)
      else  !}{
        id_diag(n) = register_diag_field ('ocean_model', trim(T_diag(n)%name), &
             Grd%tracer_axes(1:3),                                             &
             Time%model_time, trim(T_diag(n)%longname), trim(T_diag(n)%units), &
             missing_value=missing_value)
      endif  !}
    enddo  !} n


    ! read diagnostic tracer restarts 
    write (stdout(),*) trim(note_header), ' Reading diagnostic tracer initial conditions and/or restarts'
    do n = 1, num_diag_tracers  !{

      write (stdout(),*)

      if (T_diag(n)%file_in .eq. ' ') then  !{
        write (stdout(),*) 'Skipping tracer ', trim(T_diag(n)%name)
      else  !}{

        write (stdout(),*) 'Reading for tracer ', trim(T_diag(n)%name)
        write (stdout(),*) trim(T_diag(n)%file_in)
        write (stdout(),*) trim(T_diag(n)%name_in)
        write (stdout(),*) T_diag(n)%scale_in
        write (stdout(),*) T_diag(n)%additive_in

        if ( (Time%init .or. T_diag(n)%init) .and. T_diag(n)%const_init_tracer) then   !{
          write (stdout(),*) 'Initializing diagnostic tracer ', trim(T_diag(n)%name),' to constant ', T_diag(n)%const_init_value 
          T_diag(n)%field(:,:,:) = T_diag(n)%const_init_value 
        else  !}{
          call read_data(T_diag(n)%file_in, T_diag(n)%name_in, T_diag(n)%field(:,:,:), Dom%domain2d, timelevel=1)
        endif  !}

        write (stdout(),*) 'Completed read of ic/restart for diagnostic tracer ', trim(T_diag(n)%name)

        if (Time%init .or. T_diag(n)%init) then  !{

            ! table driven scaling and addition
            if (T_diag(n)%scale_in .ne. 1.0) then  !{
                write (stdout(),*) 'Scaling the initial condition for ', trim(T_diag(n)%name),' by the value ',T_diag(n)%scale_in
                T_diag(n)%field(:,:,:) = T_diag(n)%field(:,:,:) * T_diag(n)%scale_in
            endif  !}
            if (T_diag(n)%additive_in .ne. 0.0) then  !{
                write (stdout(),*) 'Adding the constant ',T_diag(n)%additive_in,' to initial condition for ',trim(T_diag(n)%name)
                T_diag(n)%field(:,:,:) = T_diag(n)%field(:,:,:) + T_diag(n)%additive_in
            endif  !}
 
            ! linearly interpolate tracers to partial bottom cells
            write (stdout(),*) 'Linearly interpolate tracer ',trim(T_diag(n)%name),' to partial cell bottom.'
            do j = jsc, jec  !{
               do i = isc, iec  !{
                  kb = Grd%kmt(i,j)
                  if (kb .gt. 1) then  !{
                      fact = Thickness%dhwt(i,j,kb-1)/Grd%dzw(kb-1)
                      T_diag(n)%field(i,j,kb) =                      &
                            (1.-fact)*T_diag(n)%field(i,j,kb-1) + fact*T_diag(n)%field(i,j,kb)
!                      T_diag(n)%field(i,j,kb) =                      &
!                           T_diag(n)%field(i,j,kb-1) -               &
!                           fact*(T_diag(n)%field(i,j,kb-1) -         &
!                           T_diag(n)%field(i,j,kb))
                  endif  !}

                  ! zero out tracers in land points
                  do k = kb+1,nk  !{
                     T_diag(n)%field(i,j,k) = 0.0
                  enddo  !} k

               enddo  !} i
            enddo  !} j

        endif  !}

        call mpp_update_domains(T_diag(n)%field(:,:,:), Dom%domain2d)

      endif  !}

    enddo  !} n
    write (stdout(),*) trim(note_header), ' Finished reading diagnostic tracer restarts.'

    diag_module_initialized = .true.

    if (debug_tracer) then
      do n= 1, num_diag_tracers
        if (T_diag(n)%file_in .ne. ' ') then  !{
          write(stdout(),*) 'Checksum from ocean_tracer_mod: initial tracer checksum (diag)  ===>'
          call tracer_diag_chksum(Time, T_diag(n))
        endif  !}
      enddo
    endif

  else  !}{

    write (stdout(),*)
    write (stdout(),*) trim(note_header), 'No diagnostic tracers requested.'

  endif  !}

end function ocean_diag_tracer_init  !}
! </FUNCTION> NAME="ocean_diag_tracer_init">


!#######################################################################
! <SUBROUTINE NAME="update_ocean_tracer">
!
! <DESCRIPTION>
! Update value of tracer concentration to time taup1.
! </DESCRIPTION>
!
subroutine update_ocean_tracer (Time, Dens, Adv_vel, Thickness, Ext_mode, pme, diff_cbt, T_prog, T_diag)

  type(ocean_time_type), intent(in)                    :: Time 
  type(ocean_density_type), intent(in)                 :: Dens
  type(ocean_adv_vel_type), intent(in)                 :: Adv_vel
  type(ocean_thickness_type), intent(in)               :: Thickness
  type(ocean_external_mode_type), intent(in)           :: Ext_mode
  real, intent(in), dimension(isd:ied,jsd:jed)         :: pme
  real, intent(in), dimension(isd:ied,jsd:jed,nk,2)    :: diff_cbt
  type(ocean_prog_tracer_type), intent(inout)          :: T_prog(num_prog_tracers)
  type(ocean_diag_tracer_type), intent(inout)          :: T_diag(num_diag_tracers)
 
  type(time_type) :: next_time
  type(time_type) :: time_step
  integer         :: i, j, k, n, nt2
  integer         :: taum1, tau, taup1, now
  
  taum1 = Time%taum1
  tau   = Time%tau
  taup1 = Time%taup1

  if (size(T_prog(:)) < num_prog_tracers) then 
    call mpp_error(FATAL,'==>Error from ocean_tracer_mod (update_ocean_tracer): T_prog array size too small')
  endif 

  if (zero_tendency) then

      ! hold tracer fields constant in time
      do n=1,num_prog_tracers
         T_prog(n)%field(:,:,:,taup1) = T_prog(n)%field(:,:,:,tau)
      enddo

  else 

     time_step = set_time(int(dtts),0)
     next_time = Time%model_time + time_step
     if(have_obc) then
        do k=1,nk
           do j=jsc,jec
              do i=isc,iec
                 wrk1(i,j,k) = Thickness%dht(i,j,k,taup1) + epsln
                 wrk2(i,j,k) = 1.0/wrk1(i,j,k)
              enddo
           enddo
        enddo
     endif

     ! compute thickness weighted time tendency for tracer arising from explicit tendency contributions. 
     do n=1,num_prog_tracers
 
        call mpp_clock_begin(id_clock_tracer_advect)
          call horz_advect_tracer(Time, Adv_vel, Thickness, T_prog(1:num_prog_tracers), T_prog(n), n, dtime)
          call vert_advect_tracer(Time, Adv_vel, T_prog(n), n)
        call mpp_clock_end(id_clock_tracer_advect)
        
        call mpp_clock_begin(id_clock_horz_diffuse)
          call horz_diffuse(Time, Thickness, T_prog(n), n)
        call mpp_clock_end(id_clock_horz_diffuse)

        call mpp_clock_begin(id_clock_vert_diffuse)
          call vert_diffuse(Time, Thickness, aidif, n, T_prog(n), T_prog(index_temp), &
                            T_prog(index_salt), diff_cbt)
        call mpp_clock_end(id_clock_vert_diffuse)
        
        ! --add pme contributions to the k=1 cell
        ! --river contributions are in T_prog%th_tendency from ocean_rivermix_mod
        ! --add surface smoother contributions due to surface height smoothing
        do j=jsc,jec
           do i=isc,iec
              T_prog(n)%th_tendency(i,j,1) = T_prog(n)%th_tendency(i,j,1) + &
                                   Grd%tmask(i,j,1)*(pme(i,j)*T_prog(n)%tpme(i,j)  &
                                   + T_prog(n)%surface_smooth(i,j) ) 
           enddo
        enddo

        ! update tracer concentration to taup1 using explicit tendency update 
        do k=1,nk
           do j=jsc,jec
              do i=isc,iec 
                 T_prog(n)%field(i,j,k,taup1) = &
                    (Thickness%dht(i,j,k,taum1)*T_prog(n)%field(i,j,k,taum1) + dtime*T_prog(n)%th_tendency(i,j,k)) &
                    *Thickness%dhtr(i,j,k)
              enddo
           enddo
        enddo

        if(have_obc) then
           call ocean_obc_tracer(T_prog(n)%field, Adv_vel%uh_et, Adv_vel%vh_nt, &
           wrk2(isc:iec,jsc:jec,:), taum1, tau, taup1, Time%model_time, T_prog(n)%name, n)
        endif

        if (debug_tracer) then
           write(stdout(),*) 'tracers after explicit tendency update (taup1) ==>'
           call tracer_prog_chksum(Time, T_prog(n), taup1)
           call tracer_min_max(Time, T_prog(n))
        endif

        if (id_tendency_array(n)> 0) then
            used = send_data(id_tendency_array(n), T_prog(n)%th_tendency(:,:,:), &
                   Time%model_time, rmask=Grd%tmask(:,:,:), &
                   is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
        endif
        if (id_surface_smooth(n)> 0) then
            used = send_data(id_surface_smooth(n), T_prog(n)%surface_smooth(:,:)*T_prog(n)%conversion, &
                   Time%model_time, rmask=Grd%tmask(:,:,1), &
                   is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
        endif

     enddo  ! enddo for n=1,num_prog_tracers


     ! compute ocean heating due to formation of frazil-ice (Joules/m^2)
     ! do this step prior to implicit vertical mixing or convection. 
     ! note that "frazil_factor" accounts for possibly different time 
     ! stepping used in ocean model and the sea ice model.  With mom4 
     ! using a leap-frog, and the GFDL ocean model SIS using forward,
     ! then frazil_factor=0.5. 
     ! If use twolevel in mom4, then frazil_factor=1.0
     if (use_frazil) then
         do j=jsc,jec
            do i=isc,iec 
               if(T_prog(index_temp)%field(i,j,1,taup1) < -0.054*T_prog(index_salt)%field(i,j,1,taup1)) then  
                   T_diag(index_frazil)%field(i,j,1) = &
                        (-0.054*T_prog(index_salt)%field(i,j,1,taup1)-T_prog(index_temp)%field(i,j,1,taup1)) &
                        *Thickness%dht(i,j,1,taup1)*rho0*cp_ocean*T_diag(index_frazil)%factor 
                   T_prog(index_temp)%field(i,j,1,taup1) = -0.054*T_prog(index_salt)%field(i,j,1,taup1)
               else 
                   T_diag(index_frazil)%field(i,j,1) = 0.0
               endif
            enddo
         enddo
     endif

     if (id_frazil_2d > 0) then 
          used = send_data(id_frazil_2d, T_diag(index_frazil)%field(:,:,1)*dtimer, &
                 Time%model_time, rmask=Grd%tmask(:,:,1), &
                 is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     endif 
    
     ! add implicit portion of vertical diffusion into tracer(tau+1)
     if (aidif > 0.0) then 

         wrk1=0.0 ; wrk2=0.0

         do n=1,num_prog_tracers

            if (n==index_salt) then
                nt2=2
            else
                nt2=1
            endif

            do k=1,nk
               do j=jsc,jec
                  do i=isc,iec
                     wrk1(i,j,k) = T_prog(n)%field(i,j,k,taup1)
                     wrk2(i,j,k) = diff_cbt(i,j,k,nt2) + T_prog(n)%K33_implicit(i,j,k)
                  enddo
               enddo
            enddo

            call invtri (T_prog(n)%field(:,:,:,taup1), T_prog(n)%stf, &
                         T_prog(n)%btf, wrk2(:,:,:), dtime, Grd%kmt,&
                         Grd%tmask, Thickness%dht(:,:,:,taup1), Thickness%dhwt, aidif, nk)  

            if (id_vdiff_impl(n) > 0) then
                do k=1,nk
                   do j=jsc,jec
                      do i=isc,iec  
                         wrk1(i,j,k) = dtimer*Thickness%dht(i,j,k,taup1)*T_prog(n)%conversion &
                                       *(T_prog(n)%field(i,j,k,taup1)-wrk1(i,j,k))
                      enddo
                   enddo
                enddo
                used = send_data(id_vdiff_impl(n), wrk1(:,:,:), &
                       Time%model_time, rmask=Grd%tmask(:,:,:), &
                       is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
            endif

            if (debug_tracer) then
                write(stdout(),*) 'tracers in update_ocean_tracer after diffusion invtri (taup1) ===>'
                call tracer_prog_chksum(Time, T_prog(n), taup1)
                call tracer_min_max(Time, T_prog(n))
            endif

         enddo

     endif  ! endif for aidif > 0

     ! convectively adjust water columns (diagnostics calculated inside convection subroutine)
     if (convective_adjust_on) then
         call convection (Time, Thickness, T_prog, Dens) 
         if (debug_tracer) then
             write(stdout(),*) 'tracers in update_ocean_tracer after convection (taup1) ===>'
             call tracer_prog_chksum(Time, T_prog(index_temp), taup1)
             call tracer_min_max(Time, T_prog(index_temp))
             call tracer_prog_chksum(Time, T_prog(index_salt), taup1)
             call tracer_min_max(Time, T_prog(index_salt))
         endif
     endif

     if(have_obc) then
         do n=1,num_prog_tracers
            call mpp_update_domains(T_prog(n)%field(:,:,:,taup1), Dom%domain2d, complete=T_prog(n)%complete)
         enddo
         do n=1,num_prog_tracers
            call ocean_obc_update_boundary(T_prog(n)%field(:,:,:,taup1), 'T')
         enddo
     endif

     ! for testing tracer advection
     do n=1,num_prog_tracers
        if(T_prog(n)%use_only_advection) then 
          call update_advection_only(Time, T_prog, Adv_vel, Thickness, Ext_mode, n)
        endif 
     enddo    


  endif  ! endif for zero_tendency


  ! send diagnostics to diag_manager 
  do n=1,num_prog_tracers

     if (id_prog(n) > 0) used = send_data (id_prog(n), T_prog(n)%field(:,:,:,tau), &
                                Time%model_time,rmask=Grd%tmask(:,:,:), &
                                is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

     if (id_surf_tracer(n) > 0) used = send_data (id_surf_tracer(n), T_prog(n)%field(:,:,1,tau), &
                                Time%model_time, rmask=Grd%tmask(:,:,1), &
                                is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

     if (id_tendency(n) > 0) then
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  wrk1(i,j,k) = dtimer*T_prog(n)%conversion &
                     *(T_prog(n)%field(i,j,k,taup1)*Thickness%dht(i,j,k,taup1) &
                      -T_prog(n)%field(i,j,k,taum1)*Thickness%dht(i,j,k,taum1))
                     
               enddo
            enddo
         enddo
         used = send_data(id_tendency(n),wrk1(:,:,:), &
                Time%model_time, rmask=Grd%tmask(:,:,:), &
                is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
     endif

  enddo

  do n=1,num_diag_tracers
     if (id_diag(n) > 0) used = send_data (id_diag(n), T_diag(n)%field(:,:,:), &
                                Time%model_time, rmask=Grd%tmask(:,:,:), &
                                is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  enddo


end subroutine update_ocean_tracer
! </SUBROUTINE> NAME="update_ocean_tracer">


!#######################################################################
! <SUBROUTINE NAME="update_advection_only">
!
! <DESCRIPTION>
! Redo tracers that are evolved using only advection--nothing else.
! Useful for testing advection schemes in Boussinesq simulations. 
!
! T_prog(n)%use_only_advection==.true. ignores all boundary forcing
! and sources, so if T_prog(n)%stf or pme, rivers, sources
! are nonzero, tracer diagnostics will spuriously indicate 
! non-conservation.  
!
! Assume for these tests that 
!  (1) vertical advection is done fully explictly
!  (2) pme, rivers, stf, btf, and other sources are zero 
!
! </DESCRIPTION>
!
subroutine update_advection_only(Time, T_prog, Adv_vel, Thickness, Ext_mode, n)

  type(ocean_time_type), intent(in)            :: Time
  type(ocean_prog_tracer_type), intent(inout)  :: T_prog(num_prog_tracers)
  type(ocean_adv_vel_type), intent(in)         :: Adv_vel
  type(ocean_thickness_type), intent(in)       :: Thickness
  type(ocean_external_mode_type), intent(in)   :: Ext_mode
  integer, intent(in)                          :: n
  integer                                      :: i, j, k, nn
  integer                                      :: taup1, taum1

  taup1 = Time%taup1  
  taum1 = Time%taum1

  if(n==1) then 
      do nn=1,num_prog_tracers
         T_prog(n)%th_tendency = 0.0 
      enddo
  endif

  call horz_advect_tracer(Time, Adv_vel, Thickness, T_prog(1:num_prog_tracers), T_prog(n), n, dtime, store_flux=.FALSE.)
  call vert_advect_tracer(Time, Adv_vel, T_prog(n), n)  

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           T_prog(n)%field(i,j,k,taup1) =  &
                (Thickness%dht(i,j,k,taum1)*T_prog(n)%field(i,j,k,taum1) + dtime*T_prog(n)%th_tendency(i,j,k)) &
                *Thickness%dhtr(i,j,k)
        enddo
     enddo
  enddo

end subroutine update_advection_only
! </SUBROUTINE> NAME="update_advection_only">



!#######################################################################
! <SUBROUTINE NAME="ocean_tracer_rstrt">
!
! <DESCRIPTION>
! Write ocean tracer restarts
! </DESCRIPTION>
!
subroutine ocean_tracer_rstrt(Time, T_prog, T_diag, ens_ocean)

  type(ocean_time_type)         :: Time
  type(ocean_prog_tracer_type)  :: T_prog(:)
  type(ocean_diag_tracer_type)  :: T_diag(:)
  logical, intent(in), optional :: ens_ocean

  integer            :: n, tau, taup1, unit
  integer            :: yr, mon, day, hr, min, sec
  integer(i8_kind)   :: chksum
  character(len=128) :: filename
  character(len=10)  :: rdte

  tau   = Time%tau
  taup1 = Time%taup1
  call get_date(Time%model_time,yr,mon,day,hr,min,sec)
  write(rdte,'(i4,3i2.2)') yr,mon,day,hr

  ! open ascii file containing checksums used to check for reproducibility
  filename = 'IRESTART/'//rdte//'ocean_tracer.res'
  call mpp_open(unit, trim(filename),form=MPP_ASCII,&
     action=MPP_OVERWR,threading=MPP_SINGLE,fileset=MPP_SINGLE,nohdrs=.true.)

  ! write the restart fields to the file T_prog(n)%file_out
  ! write the chksums to the ascii file RESTART/ocean_tracer.res
  do n=1,num_prog_tracers

     filename = 'IRESTART/'//rdte//T_prog(n)%file_out(9:)

     if(tendency==THREE_LEVEL) then
         call write_data(filename, trim(T_prog(n)%name), T_prog(n)%field(:,:,:,tau), domain = Dom%domain2d,&
                         append_pelist_name = ens_ocean)

         write(stdout(),*) 'Ending ',trim(T_prog(n)%name), ' chksum (tau) ==>'
         call tracer_prog_chksum(Time, T_prog(n), tau, chksum)
         if (mpp_pe() == mpp_root_pe()) write(unit,*) trim(T_prog(n)%name)//'   tau chksum =', chksum
     endif

     call write_data(filename, trim(T_prog(n)%name), T_prog(n)%field(:,:,:,taup1), domain = Dom%domain2d,&
                     append_pelist_name = ens_ocean)

     write(stdout(),*) 'Ending ',trim(T_prog(n)%name), ' chksum (taup1) ==>'
     call tracer_prog_chksum(Time, T_prog(n), taup1, chksum)
     if (mpp_pe() == mpp_root_pe()) write(unit,*) trim(T_prog(n)%name)//' taup1 chksum =', chksum

  enddo

  ! write the restart fields to the file T_diag(n)%file_out
  ! write the chksums to the ascii file RESTART/ocean_tracer.res where writing restart
  do n=1,num_diag_tracers

    filename = 'IRESTART/'//rdte//T_diag(n)%file_out(9:)

    if (filename .ne. ' ') then
      call write_data(filename, trim(T_diag(n)%name), T_diag(n)%field(:,:,:), domain = Dom%domain2d,&
                      append_pelist_name = ens_ocean)

      write(stdout(),*) 'Ending ',trim(T_diag(n)%name), ' chksum       ==>'
      call tracer_diag_chksum(Time, T_diag(n), chksum)
      if (mpp_pe() == mpp_root_pe()) then
        write(unit,*) trim(T_diag(n)%name)//'       chksum =', chksum
      endif

    endif

  enddo

  ! may need restart for neutral physics
  call ocean_neutral_physics_rstrt(Time,ens_ocean)

  return

end subroutine ocean_tracer_rstrt
! </SUBROUTINE NAME="ocean_tracer_rstrt">


!#######################################################################
! <SUBROUTINE NAME="ocean_tracer_end">
!
! <DESCRIPTION>
! Write ocean tracer restarts
! </DESCRIPTION>
!
subroutine ocean_tracer_end(Time, T_prog, T_diag, ens_ocean)

  type(ocean_time_type)         :: Time
  type(ocean_prog_tracer_type)  :: T_prog(:)
  type(ocean_diag_tracer_type)  :: T_diag(:)
  logical, intent(in), optional :: ens_ocean

  integer            :: n, tau, taup1, unit
  integer(i8_kind)   :: chksum
  character(len=128) :: filename

  tau   = Time%tau
  taup1 = Time%taup1

  ! open ascii file containing checksums used to check for reproducibility
  filename = 'RESTART/ocean_tracer.res'
  call mpp_open(unit, trim(filename),form=MPP_ASCII,&
     action=MPP_OVERWR,threading=MPP_SINGLE,fileset=MPP_SINGLE,nohdrs=.true.)

  ! write the restart fields to the file T_prog(n)%file_out
  ! write the chksums to the ascii file RESTART/ocean_tracer.res
  do n=1,num_prog_tracers

     filename = T_prog(n)%file_out

     if(tendency==THREE_LEVEL) then 
         call write_data(filename, trim(T_prog(n)%name), T_prog(n)%field(:,:,:,tau), domain = Dom%domain2d,&
                         append_pelist_name = ens_ocean)

         write(stdout(),*) 'Ending ',trim(T_prog(n)%name), ' chksum (tau) ==>'
         call tracer_prog_chksum(Time, T_prog(n), tau, chksum)
         if (mpp_pe() == mpp_root_pe()) write(unit,*) trim(T_prog(n)%name)//'   tau chksum =', chksum     
     endif

     call write_data(filename, trim(T_prog(n)%name), T_prog(n)%field(:,:,:,taup1), domain = Dom%domain2d,&
                     append_pelist_name = ens_ocean)

     write(stdout(),*) 'Ending ',trim(T_prog(n)%name), ' chksum (taup1) ==>'
     call tracer_prog_chksum(Time, T_prog(n), taup1, chksum)
     if (mpp_pe() == mpp_root_pe()) write(unit,*) trim(T_prog(n)%name)//' taup1 chksum =', chksum

  enddo

  ! write the restart fields to the file T_diag(n)%file_out
  ! write the chksums to the ascii file RESTART/ocean_tracer.res where writing restart
  do n=1,num_diag_tracers

    filename = T_diag(n)%file_out

    if (filename .ne. ' ') then  
      call write_data(filename, trim(T_diag(n)%name), T_diag(n)%field(:,:,:), domain = Dom%domain2d,&
                      append_pelist_name = ens_ocean)

      write(stdout(),*) 'Ending ',trim(T_diag(n)%name), ' chksum       ==>'
      call tracer_diag_chksum(Time, T_diag(n), chksum)
      if (mpp_pe() == mpp_root_pe()) then  
        write(unit,*) trim(T_diag(n)%name)//'       chksum =', chksum     
      endif 

    endif  

  enddo

  ! may need restart for neutral physics
  call ocean_neutral_physics_end(ens_ocean)

  return

end subroutine ocean_tracer_end
! </SUBROUTINE> NAME="ocean_tracer_end">


!#######################################################################
! <SUBROUTINE NAME="compute_tmask_limit">
!
! <DESCRIPTION>
! Provide for possibility that quicker advection and/or neutral physics 
! reverts to first order upwind when tracer is outside a specified range.
! For this purpose, we define a mask which is set to unity where 
! fluxes revert to first order upwind advection 
! (if using quicker) and horizontal diffusion (if using neutral). 
!
! This method is very ad hoc.  What is preferred for advection is to use
! a monotonic scheme, such as mdfl_sweby.  For neutral physics, no 
! analogous monotonic scheme has been implemented in mom4.  Such could be 
! useful, especially for passive tracers. In the meantime, tmask_limit 
! provides a very rough limiter for neutral physics to help keep tracers 
! within specified bounds.
!
! </DESCRIPTION>
!
subroutine compute_tmask_limit(Time, T_prog)

  type(ocean_time_type), intent(in)           :: Time
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)

  integer                                     :: tau
  integer                                     :: i, j, k, n 

  tau = Time%tau
  
  do n = 1,num_prog_tracers 

     do k=1,nk 
        do j=jsc-1,jec
           do i=isc-1,iec

              T_prog(n)%tmask_limit(i,j,k) = 0.0   

              ! is tracer outside range? 
              if(Grd%tmask(i,j,k) == 1.0) then 
                  if(T_prog(n)%field(i,j,k,tau) < T_prog(n)%min_tracer_limit .or. &
                     T_prog(n)%field(i,j,k,tau) > T_prog(n)%max_tracer_limit) then 
                     T_prog(n)%tmask_limit(i,j,k) = Grd%tmask(i,j,k)
                  endif
              endif

              ! is east tracer outside range? 
              if(Grd%tmask(i+1,j,k) == 1.0) then 
                  if(T_prog(n)%field(i+1,j,k,tau) < T_prog(n)%min_tracer_limit .or.  &
                     T_prog(n)%field(i+1,j,k,tau) > T_prog(n)%max_tracer_limit) then 
                     T_prog(n)%tmask_limit(i,j,k) = Grd%tmask(i,j,k)
                  endif
              endif

              ! is north tracer outside range? 
              if(Grd%tmask(i,j+1,k) == 1.0) then 
                  if(T_prog(n)%field(i,j+1,k,tau) < T_prog(n)%min_tracer_limit .or.  &
                     T_prog(n)%field(i,j+1,k,tau) > T_prog(n)%max_tracer_limit) then 
                     T_prog(n)%tmask_limit(i,j,k) = Grd%tmask(i,j,k)
                  endif
              endif

           enddo
        enddo
     enddo

     ! is deeper tracer outside range? 
     do k=1,nk-1 
        do j=jsc-1,jec
           do i=isc-1,iec
              if(Grd%tmask(i,j,k+1) == 1.0) then 
                  if(T_prog(n)%field(i,j,k+1,tau) < T_prog(n)%min_tracer_limit .or.  &
                     T_prog(n)%field(i,j,k+1,tau) > T_prog(n)%max_tracer_limit) then 
                     T_prog(n)%tmask_limit(i,j,k) = Grd%tmask(i,j,k)
                  endif
              endif
           enddo
        enddo
     enddo

  enddo ! end of loop n=1,num_prog_tracers  


  ! insist that if temperature or salinity tmask_limit is 1.0, then both must be 1.0
  do k=1,nk 
     do j=jsc-1,jec
        do i=isc-1,iec
           if(T_prog(index_temp)%tmask_limit(i,j,k)==1.0 .or. T_prog(index_salt)%tmask_limit(i,j,k)==1.0) then 
              T_prog(index_temp)%tmask_limit(i,j,k)=1.0 
              T_prog(index_salt)%tmask_limit(i,j,k)=1.0   
           endif
        enddo
     enddo
  enddo

  ! debugging and diagnostics  
  do n=1,num_prog_tracers 

     if(debug_tracer) then 
         write(stdout(),*) ' ' 
         write(stdout(),*) 'tmask_limit chksum for ',trim(T_prog(n)%name),' = ',&
              mpp_chksum(T_prog(n)%tmask_limit(isc:iec,jsc:jec,:))
     endif

     if (id_tmask_limit(n) > 0) then 
         used = send_data(id_tmask_limit(n), T_prog(n)%tmask_limit(:,:,:), &
                Time%model_time, rmask=Grd%tmask(:,:,:), &
                is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
     endif

  enddo


end subroutine compute_tmask_limit
! </SUBROUTINE> NAME="compute_tmask_limit"


     
end module ocean_tracer_mod
