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
! Modified: Xingren Wu
!           Xingren.Wu@noaa.gov
module ocean_sbc_mod
!  
!<CONTACT EMAIL="Matthew.Harrison@noaa.gov"> M.J. Harrison
!</CONTACT>
!
! <REVIEWER EMAIL="Tony.Rosati@noaa.gov"> A. Rosati 
! </REVIEWER>
!
! <REVIEWER EMAIL="Stephen.Griffies@noaa.gov"> S. M. Griffies 
! </REVIEWER>
!
! <REVIEWER EMAIL="V.Balaji@noaa.gov">
! V. Balaji
! </REVIEWER>
!
!<OVERVIEW>
! Set up the surface boundary conditions for mom4. 
!</OVERVIEW>
!
!<DESCRIPTION>
! This module sets up the surface boundary conditions for the model. 
! Also fill Ocean_sfc derived-type used to pass information to other 
! component models.
!</DESCRIPTION>
!
!<NAMELIST NAME="ocean_sbc_nml">
!
!  <DATA NAME="use_waterflux" TYPE="logical">
!  Set to true when wish to use real fresh water flux as opposed to virtual 
!  salt fluxes.   
!  </DATA> 
!  <DATA NAME="waterflux_tavg" TYPE="logical">
!  Set to true when aiming to suppress the leap-frog computational mode
!  by setting pme and river equal to a time averaged value over the 
!  present and previous time step.  This method requires an extra
!  restart file.  
!  </DATA> 
!
!  <DATA NAME="temp_restore_tscale" UNITS="day" TYPE="real">
!  Time scale in days for restoring temperature within the top model 
!  grid cell. 
!  </DATA> 
!  <DATA NAME="salinity_ref" UNITS="psu" TYPE="real">
!  Reference salinity used for converting fresh water flux
!  to salt flux. 
!  </DATA> 
!  <DATA NAME="salt_restore_tscale" UNITS="day" TYPE="real">
!  Time scale in days for restoring salinity within the top model 
!  grid cell. 
!  </DATA> 
!  <DATA NAME="salt_restore_under_ice" TYPE="logical">
!  Logical indicating whether to restore salinity under sea ice or not.
!  When .false. then will not restore salinity  in regions where we 
!  use a "frazil" condition as a proxy for where sea-ice is present.
!  Do not use sea ice extent from a sea ice model since we generally do 
!  not pass information regarding ice extent between the sea ice model 
!  and the ocean model.     
!  </DATA> 
!  <DATA NAME="zero_net_salt_restore" TYPE="logical">
!  Logical indicating whether to remove the area mean of the salinity 
!  restore flux so there is a net zero input of salt to the ocean
!  associated with restoring.    
!  </DATA> 
!  <DATA NAME="zero_net_water_restore" TYPE="logical">
!  Logical indicating whether to remove the area mean of the water 
!  restore flux so there is a net zero input of water to the ocean
!  associated with restoring.    
!  </DATA> 
!  <DATA NAME="zero_net_water_coupler" TYPE="logical">
!  Logical indicating whether to remove the area mean of the water 
!  passed through the coupler so there is a net zero input of 
!  fresh water to the ocean associated with p-e+r. Do so by removing 
!  area mean from pme--keep river values unchanged. Note that a choice
!  must be made whether to remove the area mean from rivers or pme.  
!  We choose pme since it is more evenly distributed than rivers.   
!  Setting zero_net_water_coupler to true may be appropriate when 
!  running an ice-ocean model using a bulk formulae to compute
!  evaporation (e.g., omip) and when only providing a weak (or zero)
!  salinity restoring.  It is not appropriate when running a coupled
!  ocean-atmosphere model, where the moisture budget should be 
!  conserved without an artificial removal of the global mean.  
!  </DATA> 
!
!  <DATA NAME="rotate_winds" TYPE="logical">
!  Set to true when need to rotate the winds onto the ocean model grid.
!  This is needed for cases where the winds are on a spherical grid and 
!  the ocean model uses tripolar=.true.  If generate the wind data on 
!  the ocean model grid, then do not need to rotate, since the rotation 
!  has already been done.  
!  </DATA> 
!
!  <DATA NAME="max_ice_thickness" UNITS="m" TYPE="real">
!  When coupling mom4 to an ice model, the sea ice thickness may need
!  to be restricted to prevent vanishing top-level in mom4. Set 
!  max_ice_thickness (meters) < dzt(k=1) to restrict. This truncation 
!  avoids the numerical problem but we loose mass conservation in the coupled
!  sea ice and ocean system. We also alter the pressure felt on the ocean 
!  as applied by the sea ice. Different vertical coordinates are needed 
!  to do the problem more realistically.   
!  </DATA> 
!
!  <DATA NAME="riverspread" TYPE="logical">
!  Set to true if wish to use the spread_river_horz algorithm to spread 
!  out the river discharge over a wide area into the ocean.  
!  </DATA> 
!
!  <DATA NAME="frazil_factor" UNITS="dimensionless" TYPE="real">
!  This factor accounts for possibly different time stepping used 
!  in the sea ice model relative to the ocean model.  If sea-ice 
!  and ocean use same time stepping schemes, then frazil_factor=1.0.
!  If sea-ice uses a twolevel scheme and ocean a threelevel leap-frog,
!  then frazil_factor=0.5. Default is 1.0 since the GFDL sea ice model 
!  SIS uses  a two-level time stepping scheme and mom4 defaults to 
!  a staggered two-level scheme. 
!  </DATA> 
!
!  <DATA NAME="avg_sfc_velocity" TYPE="logical">
!  If set to true, the u and v fields passed up to the sea ice
!  are averaged over a coupling interval. TRUE by default.
!  </DATA> 
!  <DATA NAME="avg_sfc_temp_salt_eta" TYPE="logical">
!  If set to true, the t, s and sea_level fields passed up to the sea ice
!  are averaged over a coupling interval. TRUE by default.
!  </DATA> 
!
!</NAMELIST>
!
use constants_mod,            only: epsln, grav, rho0, rho0r, rho_cp, hlv, hlf, Tfreeze
use diag_manager_mod,         only: register_diag_field, send_data
use fms_mod,                  only: open_namelist_file, check_nml_error, file_exist
use fms_mod,                  only: close_file, read_data, write_data
use mpp_domains_mod,          only: mpp_update_domains, BGRID_NE, mpp_define_domains, mpp_get_compute_domain
use mpp_domains_mod,          only: mpp_global_sum, BITWISE_EXACT_SUM
use mpp_mod,                  only: mpp_error, mpp_chksum, FATAL, stdout, stdlog
use time_interp_external_mod, only: time_interp_external, init_external_field
use time_manager_mod,         only: time_type, get_date

use ocean_domains_mod,        only: get_local_indices
use ocean_riverspread_mod,    only: spread_river_horz
use ocean_types_mod,          only: ocean_grid_type, ocean_domain_type, ocean_data_type
use ocean_types_mod,          only: ocean_time_type, ocean_thickness_type 
use ocean_types_mod,          only: ocean_prog_tracer_type, ocean_diag_tracer_type
use ocean_types_mod,          only: ocean_external_mode_type, ocean_velocity_type 
use ocean_types_mod,          only: missing_value, ice_ocean_boundary_type
  
implicit none

private

integer, allocatable, dimension(:) :: id_restore
integer :: index_temp=-1
integer :: index_salt=-1
integer :: index_frazil=-1
integer :: memuse
integer :: num_prog_tracers
integer :: num_diag_tracers

! ids for diagnostic manager 
integer :: id_tau_x=-1
integer :: id_tau_y=-1
integer :: id_river=-1
integer :: id_swflx=-1
integer, allocatable, dimension(:) :: id_stf
integer, allocatable, dimension(:) :: id_stf_adj
integer, allocatable, dimension(:) :: id_stf_total
integer, allocatable, dimension(:) :: id_stf_river
integer :: id_pme            =-1
integer :: id_pme_adj        =-1
integer :: id_pme_total      =-1
integer :: id_ice_mask       =-1
integer :: id_open_ocean_mask=-1

#include <ocean_memory.h>

#ifdef  STATIC_MEMORY
real, dimension(isd:ied,jsd:jed) :: data
real, dimension(isd:ied,jsd:jed) :: pme_taum1   ! pme from coupler at the taum1 time step 
real, dimension(isd:ied,jsd:jed) :: river_taum1 ! river from coupler at the taum1 time step 
#else
real, allocatable, dimension(:,:) :: data
real, allocatable, dimension(:,:) :: pme_taum1   ! pme from coupler at the taum1 time step 
real, allocatable, dimension(:,:) :: river_taum1 ! river from coupler at the taum1 time step 
#endif

! ice-ocean-boundary fields are allocated using absolute
! indices (regardless of whether ocean allocations are static)
integer :: isc_bnd, iec_bnd, jsc_bnd, jec_bnd
integer :: i_shift, j_shift                      ! shift isc_bnd to isc and jsc_bnd to jsc

type(ocean_grid_type), pointer   :: Grd =>NULL()
type(ocean_domain_type), pointer :: Dom =>NULL()

public :: ocean_sbc_init
public :: sum_ocean_sfc
public :: avg_ocean_sfc
public :: zero_ocean_sfc
public :: flux_adjust
public :: get_ocean_sbc
public :: initialize_ocean_sfc
public :: ocean_sfc_rstrt
public :: ocean_sfc_end

logical, public :: use_waterflux=.true.
logical :: waterflux_tavg=.false.
logical :: rotate_winds=.false.
logical :: riverspread=.false.
logical :: salt_restore_under_ice=.true.
logical :: zero_net_salt_restore=.false.
logical :: zero_net_water_restore=.false.
logical :: zero_net_water_coupler=.false.
real    :: temp_restore_tscale=30.
real    :: salt_restore_tscale=30.
real    :: max_ice_thickness=5.0 
real    :: frazil_factor=1.0
real    :: salinity_ref=35.0
real    :: temp_damp_factor
real    :: salt_damp_factor
logical :: avg_sfc_velocity=.true.
logical :: avg_sfc_temp_salt_eta=.true.

namelist /ocean_sbc_nml/ temp_restore_tscale, salt_restore_tscale, salt_restore_under_ice, &
                         rotate_winds, use_waterflux, waterflux_tavg, max_ice_thickness, riverspread, frazil_factor, &
                         salinity_ref, zero_net_salt_restore, zero_net_water_restore, zero_net_water_coupler, &
                         avg_sfc_velocity, avg_sfc_temp_salt_eta

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_sbc_init">
!
! <DESCRIPTION>
! Initialize the ocean sbc module 
! </DESCRIPTION>
!
subroutine ocean_sbc_init(Grid, Domain, Time, T_prog, T_diag, Velocity, Ocean_sfc, time_tendency)

  type(ocean_grid_type), intent(in), target           :: Grid
  type(ocean_domain_type), intent(in), target         :: Domain
  type(ocean_time_type), intent(in)                   :: Time
  type(ocean_prog_tracer_type), intent(inout), target :: T_prog(:)
  type(ocean_diag_tracer_type), intent(inout), target :: T_diag(:)
  type(ocean_velocity_type), intent(in), target       :: Velocity
  type(ocean_data_type), intent(inout)                :: Ocean_sfc
  character(len=32), intent(in)                       :: time_tendency 

  integer            :: taup1, n, ioun, ierr, io_status
  real               :: secday
  character(len=128) :: name

  ioun = open_namelist_file()
  read  (ioun, ocean_sbc_nml,iostat=io_status)
  write (stdout(),'(/)')
  write (stdout(), ocean_sbc_nml)  
  write (stdlog(), ocean_sbc_nml)
  ierr = check_nml_error(io_status,'ocean_sbc_nml')
  call close_file (ioun)

#ifndef STATIC_MEMORY  
  call get_local_indices(Domain,isd, ied, jsd, jed, isc, iec, jsc, jec)
  allocate(data(isd:ied,jsd:jed))
  allocate(pme_taum1(isd:ied,jsd:jed))
  allocate(river_taum1(isd:ied,jsd:jed))
#endif
  data       =0.0
  pme_taum1  =0.0
  river_taum1=0.0

  Dom => Domain
  Grd => Grid
  
  call mpp_define_domains((/1,Grd%ni,1,Grd%nj/),Dom%layout,Ocean_sfc%Domain)  
  call mpp_get_compute_domain(Ocean_sfc%Domain, isc_bnd, iec_bnd, &
                              jsc_bnd, jec_bnd)

  i_shift = isc - isc_bnd
  j_shift = jsc - jsc_bnd

  allocate ( Ocean_sfc%t_surf (isc_bnd:iec_bnd,jsc_bnd:jec_bnd), &
             Ocean_sfc%s_surf (isc_bnd:iec_bnd,jsc_bnd:jec_bnd), &
             Ocean_sfc%u_surf (isc_bnd:iec_bnd,jsc_bnd:jec_bnd), &
             Ocean_sfc%v_surf (isc_bnd:iec_bnd,jsc_bnd:jec_bnd), &
             Ocean_sfc%sea_lev(isc_bnd:iec_bnd,jsc_bnd:jec_bnd), &
             Ocean_sfc%frazil (isc_bnd:iec_bnd,jsc_bnd:jec_bnd))

  Ocean_sfc%t_surf  = 0.0  ! time averaged sst (Kelvin) passed to atmosphere/ice model
  Ocean_sfc%s_surf  = 0.0  ! time averaged sss (psu) passed to atmosphere/ice models
  Ocean_sfc%u_surf  = 0.0  ! time averaged u-current (m/sec) passed to atmosphere/ice models
  Ocean_sfc%v_surf  = 0.0  ! time averaged v-current (m/sec)  passed to atmosphere/ice models 
  Ocean_sfc%sea_lev = 0.0  ! time averaged thickness of top model grid cell (m) plus patm/rho0/grav 
  Ocean_sfc%frazil  = 0.0  ! time accumulated frazil (J/m^2) passed to ice model

  num_prog_tracers = size(T_prog)
  num_diag_tracers = size(T_diag)

  allocate( id_restore  (num_prog_tracers) )
  allocate( id_stf      (num_prog_tracers) )
  allocate( id_stf_adj  (num_prog_tracers) )
  allocate( id_stf_total(num_prog_tracers) )
  allocate( id_stf_river(num_prog_tracers) )

  id_restore  (:) = -1
  id_stf      (:) = -1
  id_stf_adj  (:) = -1
  id_stf_total(:) = -1
  id_stf_river(:) = -1

  do n=1, num_prog_tracers
     if (T_prog(n)%name == 'temp') index_temp = n
     if (T_prog(n)%name == 'salt') index_salt = n
  enddo

  do n=1, num_diag_tracers
     T_diag(n)%factor = 1.0
     if (T_diag(n)%name == 'frazil') then 
       index_frazil = n
       T_diag(n)%factor = frazil_factor  
       write(stdout(),'(/a,f5.2)') ' ==>Note from ocean_sbc_mod: using frazil_factor = ',frazil_factor
       if(time_tendency =='twolevel' .and. frazil_factor==0.5) then 
          call mpp_error(FATAL,'==>ocean_sbc_mod: with 2-level tendencies for mom4 and SIS, set frazil_factor==1.0 ')
          write(stdout(),'(/a)') '==>Error in ocean_sbc_mod: SIS is two-time ice level model.  You are also running mom4 with'
          write(stdout(),'(a) ') '   two-time level scheme. In this case, frazil_factor==1.0 should be set.'
          write(stdout(),'(a/) ') '   If using another ice model with different time stepping, set frazil_factor accordingly.'
       endif 
       if(time_tendency =='threelevel' .and. frazil_factor==1.0) then 
          call mpp_error(FATAL,'==>ocean_sbc_mod: with 3-level tendency for mom4 and 2-level for SIS, set frazil_factor==.50 ')
          write(stdout(),'(/a)') '==>Error in ocean_sbc_mod: SIS is two-time level ice model.  You are running mom4 with '
          write(stdout(),'(a) ') '   three-time level leap-frog. In this case, frazil_factor==.50 should be set.'
          write(stdout(),'(a/) ') '   If using another ice model with different time stepping, set frazil_factor accordingly.'
       endif 
     endif 
  enddo


  if (index_temp == -1 .or. index_salt == -1) then 
    call mpp_error(FATAL,'==>Error in ocean_sbc_mod (ocean_sbc_init): temp and/or salt not identified')
  endif 
  
  taup1 = Time%taup1

  do n = 1, num_prog_tracers
     ! init_external_field(file_name,field_name,domain)
#ifndef STATIC_MEMORY         
        allocate(T_prog(n)%stf(isd:ied,jsd:jed))
        allocate(T_prog(n)%tpme(isd:ied,jsd:jed))
        allocate(T_prog(n)%triver(isd:ied,jsd:jed))
        allocate(T_prog(n)%riverdiffuse(isd:ied,jsd:jed))
#endif
        name = 'INPUT/'//trim(T_prog(n)%name)//'_sfc_restore.nc'
        if (file_exist(trim(name))) then
            id_restore(n) = init_external_field(name, &
                 T_prog(n)%name, domain=Dom%domain2d)
            write(stdout(),*) '==>Note from ocean_sbc_mod: applying surface restoring to '//trim(T_prog(n)%name)
            if (id_restore(n) == -1) call mpp_error(FATAL,'==>Error in ocean_sbc_mod: failure to find sfc_restore field') 
        endif
        T_prog(n)%stf          = 0.0
        T_prog(n)%tpme         = 0.0
        T_prog(n)%triver       = 0.0
        T_prog(n)%riverdiffuse = 0.0
        if (n == index_temp) then
            id_stf(n) = register_diag_field('ocean_model','sfc_hflux', &
                 Grd%tracer_axes(1:2),&
                 Time%model_time, 'surface heat flux', &
                 'Watts/m^2' ,&
                 missing_value=missing_value,range=(/-1.e4,1.e4/))
            id_stf_adj(n) = register_diag_field('ocean_model','sfc_hflux_adj', &
                 Grd%tracer_axes(1:2),&
                 Time%model_time, 'surface heat flux adjustment', &
                 'Watts/m^2' ,&
                 missing_value=missing_value,range=(/-1.e4,1.e4/))            
            id_stf_total(n) = register_diag_field('ocean_model','sfc_hflux_total', &
                 Grd%tracer_axes(1:2),&
                 Time%model_time, 'total surface heat flux', &
                 'Watts/m^2' ,&
                 missing_value=missing_value,range=(/-1.e4,1.e4/))            
            id_stf_river(n) = register_diag_field('ocean_model','sfc_hflux_river', &
                 Grd%tracer_axes(1:2),&
                 Time%model_time, 'heat flux from calving and runoff', &
                 'Watts/m^2' ,&
                 missing_value=missing_value,range=(/-1.e4,1.e4/))            
        else           
            name = 'sfc_'//trim(T_prog(n)%name)//'_flux'
            id_stf(n) = register_diag_field('ocean_model',trim(name), &
                 Grd%tracer_axes(1:2),&
                 Time%model_time, trim(name), &
                 'kg/m^2/sec' ,&
                 missing_value=missing_value,range=(/-1.e4,1.e4/))
            name = 'sfc_'//trim(T_prog(n)%name)//'_flux_adj'
            id_stf_adj(n) = register_diag_field('ocean_model',&
                 trim(name), &
                 Grd%tracer_axes(1:2),&
                 Time%model_time, trim(name), &
                 'kg/m^2/sec' ,&
                 missing_value=missing_value,range=(/-1.e4,1.e4/))            
            name = 'sfc_'//trim(T_prog(n)%name)//'_flux_total'
            id_stf_total(n) = register_diag_field('ocean_model',&
                 trim(name), &
                 Grd%tracer_axes(1:2),&
                 Time%model_time, trim(name), &
                 'kg/m^2/sec' ,&
                 missing_value=missing_value,range=(/-1.e4,1.e4/))            
            name = 'sfc_'//trim(T_prog(n)%name)//'_flux_river'
            id_stf_river(n) = register_diag_field('ocean_model',&
                 trim(name), &
                 Grd%tracer_axes(1:2),&
                 Time%model_time, trim(name), &
                 'kg/m^2/sec' ,&
                 missing_value=missing_value,range=(/-1.e4,1.e4/))            
        endif
  enddo

  secday = 1.0/(60.0*1440.0)
  if (temp_restore_tscale > 0.0) then
     temp_damp_factor = Grd%dzt(1)*secday/temp_restore_tscale
  else
     temp_damp_factor = -1.0
     write(stdout(),*) '==>Note from ocean_sbc_mod: temp_restore_tscale < 0. will not apply surface restoring to temp'
  endif

  if (salt_restore_tscale > 0.0) then
     salt_damp_factor = Grd%dzt(1)*secday/salt_restore_tscale

     if(salt_restore_under_ice) then 
       write(stdout(),*) '==>Note from ocean_sbc_mod: salt_restore_under_ice=.true. will restore salinity even under sea ice.'
     else 
       write(stdout(),*) '==>Note from ocean_sbc_mod: salt_restore_under_ice=.false. will not restore salinity under sea ice.'
     endif 

     if(use_waterflux) then 
       if(zero_net_water_restore) then 
          write(stdout(),*) '==>Note from ocean_sbc_mod: zero_net_water_restore=.true.=>zero net restoring water put in ocean.'
       else 
          write(stdout(),*) '==>Note from ocean_sbc_mod: zero_net_water_restore=.false.=>nonzero net restoring water put in ocean.'
       endif 
       if(zero_net_water_coupler) then 
          write(stdout(),*) '==>Note from ocean_sbc_mod: zero_net_water_coupler=.true.=>zero net water into ocean via coupler.'
       else 
          write(stdout(),*) '==>Note from ocean_sbc_mod: zero_net_water_coupler=.false.=>nonzero net water into ocean via coupler.'
       endif 
     else 
       if(zero_net_salt_restore) then 
          write(stdout(),*) '==>Note from ocean_sbc_mod: zero_net_salt_restore=.true.=>zero net restoring salt put in ocean.'
       else 
          write(stdout(),*) '==>Note from ocean_sbc_mod: zero_net_salt_restore=.false.=>nonzero net restoring salt put in ocean.'
       endif 
     endif 

  else

     salt_damp_factor = -1.0
     write(stdout(),*) '==>Note from ocean_sbc_mod: salt_restore_tscale < 0 so will not apply surface restoring to salt.'

  endif

  if(avg_sfc_velocity) then 
    write(stdout(),*) '==>If coupling, then avg_sfc_velocity=.true. means will pass averaged ocean velocity to ice model.'
  else 
    write(stdout(),*) '==>If coupling, then avg_sfc_velocity=.false. means will pass most recent ocean velocity to ice model.'
  endif 
  if(avg_sfc_temp_salt_eta) then 
    write(stdout(),*) '==>If coupling, then avg_sfc_temp_salt_eta=.true. means will pass averaged sst, sss, eta to ice model.'
  else 
    write(stdout(),*) '==>If coupling, then avg_sfc_velocity=.false. means will pass most recent sst, sss, eta to ice model.'
  endif 

  id_tau_x = register_diag_field('ocean_model','tau_x', Grd%vel_axes_uv(1:2),&
       Time%model_time, 'i-directed wind stress', 'N/m^2',&
       missing_value=missing_value,range=(/-10.,10./))

  id_tau_y = register_diag_field('ocean_model','tau_y', Grd%vel_axes_uv(1:2),&
       Time%model_time, 'j-directed wind stress', 'N/m^2' ,&
       missing_value=missing_value,range=(/-10.,10./)) 

  id_ice_mask = register_diag_field('ocean_model','ice_mask', Grd%tracer_axes(1:2),&
       Time%model_time, 'ice mask according to near-frazil condition', 'none' ,&
       missing_value=missing_value,range=(/-10.,10./))

  id_open_ocean_mask = register_diag_field('ocean_model','open_ocean_mask', Grd%tracer_axes(1:2),&
       Time%model_time, 'open-ocean mask according to near-frazil condition', 'none' ,&
       missing_value=missing_value,range=(/-10.,10./))

  id_river = register_diag_field('ocean_model','river', Grd%tracer_axes(1:2),&
       Time%model_time, 'river water flux', 'm/sec' ,&
       missing_value=missing_value,range=(/-10.,10./))

  id_pme = register_diag_field('ocean_model','pme', Grd%tracer_axes(1:2),&
       Time%model_time, 'precip-evap (liquid, frozen, evaporation)', 'm/sec' ,&
       missing_value=missing_value,range=(/-10.,10./))

  id_pme_adj = register_diag_field('ocean_model','pme_adj', Grd%tracer_axes(1:2),&
       Time%model_time, 'precip-evap (from restoring)', 'm/sec' ,&
       missing_value=missing_value,range=(/-10.,10./))   

  id_pme_total = register_diag_field('ocean_model','pme_total', Grd%tracer_axes(1:2),&
       Time%model_time, 'precip-evap (total including restore)', 'm/sec' ,&
       missing_value=missing_value,range=(/-10.,10./))   

  id_swflx = register_diag_field('ocean_model','swflx', Grd%tracer_axes(1:2),&
       Time%model_time, 'shortwave flux into ocean', 'W/m^2' ,&
       missing_value=missing_value,range=(/-1.e10,1.e10/))   

  if(waterflux_tavg) then 
    write(stdout(),*) '==>Note: waterflux_tavg sets pme+river = avg of ice_ocean_boundary values over'
    write(stdout(),*) '         tau and taum1 time steps. This may damp splitting between leap-frog modes.'

    if(time_tendency=='twolevel') then 
      write(stdout(),'(/a)') '==>Warning in ocean_sbc_mod: waterflux_tavg=.true. unnecessary with time_tendency==twolevel.'
    endif 

    if (file_exist('INPUT/ocean_waterflux.res.nc')) then
        call read_data('INPUT/ocean_waterflux.res','pme_taum1',   pme_taum1,  Domain%domain2d)
        call read_data('INPUT/ocean_waterflux.res','river_taum1', river_taum1,Domain%domain2d)
    endif

  endif 
  
  return

end subroutine ocean_sbc_init
! </SUBROUTINE> NAME="ocean_sbc_init"


 !#######################################################################
! <SUBROUTINE NAME="initialize_ocean_sfc">
!
! <DESCRIPTION>
! Initialize the ocean surface type, which passes information between ocean 
! and other component models. 
!
!  Ocean_sfc%t_surf  = time averaged sst (Kelvin) passed to atmosphere/ice model
!  Ocean_sfc%s_surf  = time averaged sss (psu) passed to atmosphere/ice models
!  Ocean_sfc%u_surf  = time averaged u-current (m/sec) passed to atmosphere/ice models
!  Ocean_sfc%v_surf  = time averaged v-current (m/sec)  passed to atmosphere/ice models 
!  Ocean_sfc%sea_lev = time averaged thickness of top model grid cell (m) plus patm/rho0/grav 
!  Ocean_sfc%frazil  = time accumulated frazil (J/m^2) passed to ice model.  time averaging 
!                      not performed, since ice model needs the frazil accumulated over the 
!                      ocean time steps.  Note that Ocean_sfc%frazil is accumulated, whereas 
!                      T_diag%frazil (saved in diagnostic tracer restart file) is instantaneous. 
!
! </DESCRIPTION>
!
subroutine initialize_ocean_sfc(Time, Thickness, T_prog, Velocity, patm, Ocean_sfc)

  type(ocean_time_type), intent(in)            :: Time
  type(ocean_thickness_type), intent(in)       :: Thickness
  type(ocean_prog_tracer_type), intent(in)     :: T_prog(num_prog_tracers)
  type(ocean_velocity_type), intent(in)        :: Velocity
  real, dimension(isd:ied,jsd:jed), intent(in) :: patm
  type(ocean_data_type), intent(inout), target :: Ocean_sfc
  
  integer :: taup1

  taup1 = Time%taup1
  
  Ocean_sfc%avg_kount = 0
  Ocean_sfc%t_surf    = Tfreeze 

  where (Grd%tmask(isc:iec,jsc:jec,1) == 1.0)
    Ocean_sfc%t_surf(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) = T_prog(index_temp)%field(isc:iec,jsc:jec,1,taup1) + Tfreeze
    Ocean_sfc%s_surf(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) = T_prog(index_salt)%field(isc:iec,jsc:jec,1,taup1)
    Ocean_sfc%u_surf(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) = Velocity%u(isc:iec,jsc:jec,1,1,taup1)
    Ocean_sfc%v_surf(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) = Velocity%u(isc:iec,jsc:jec,1,2,taup1)
    Ocean_sfc%sea_lev(isc_bnd:iec_bnd,jsc_bnd:jec_bnd)= Thickness%dht(isc:iec,jsc:jec,1,taup1) + patm(isc:iec,jsc:jec)/grav/rho0
    Ocean_sfc%frazil(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) = 0.0
  end where 

  if (file_exist('INPUT/ocean_sbc.res.nc')) then
      call read_data('INPUT/ocean_sbc.res','t_surf', Ocean_sfc%t_surf,Ocean_sfc%Domain)
      call read_data('INPUT/ocean_sbc.res','s_surf', Ocean_sfc%s_surf,Ocean_sfc%Domain)
      call read_data('INPUT/ocean_sbc.res','u_surf', Ocean_sfc%u_surf,Ocean_sfc%Domain)
      call read_data('INPUT/ocean_sbc.res','v_surf', Ocean_sfc%v_surf,Ocean_sfc%Domain)
      call read_data('INPUT/ocean_sbc.res','sea_lev',Ocean_sfc%sea_lev,Ocean_sfc%Domain)
      call read_data('INPUT/ocean_sbc.res','frazil', Ocean_sfc%frazil,Ocean_sfc%Domain)
  endif

  return

end subroutine initialize_ocean_sfc
! </SUBROUTINE> NAME="initialize_ocean_sfc"


!#######################################################################
! <SUBROUTINE NAME="sum_ocean_sfc">
!
! <DESCRIPTION>
! Accumulate the ocean_sfc derived type over the course of the 
! ocean component sub-cycling used when coupling to other models. 
! </DESCRIPTION>
!  
subroutine sum_ocean_sfc(Time, Thickness, T_prog, T_diag, Velocity, Ocean_sfc)

  type(ocean_time_type), intent(in)            :: Time
  type(ocean_thickness_type), intent(in)       :: Thickness
  type(ocean_prog_tracer_type), intent(in)     :: T_prog(num_prog_tracers)
  type(ocean_diag_tracer_type), intent(in)     :: T_diag(num_diag_tracers)
  type(ocean_velocity_type), intent(in)        :: Velocity
  type(ocean_data_type), intent(inout), target :: Ocean_sfc
  
  integer :: taup1, i, j, ii, jj

  if (Ocean_sfc%avg_kount == 0) call zero_ocean_sfc(Ocean_sfc)

  taup1 = Time%taup1

  Ocean_sfc%avg_kount = Ocean_sfc%avg_kount + 1

  do j = jsc_bnd,jec_bnd
     do i = isc_bnd,iec_bnd
        ii = i + i_shift
        jj = j + j_shift
        Ocean_sfc%t_surf(i,j)  = Ocean_sfc%t_surf(i,j) + T_prog(index_temp)%field(ii,jj,1,taup1) 
        Ocean_sfc%s_surf(i,j)  = Ocean_sfc%s_surf(i,j) + T_prog(index_salt)%field(ii,jj,1,taup1)
        Ocean_sfc%u_surf(i,j)  = Ocean_sfc%u_surf(i,j) + Velocity%u(ii,jj,1,1,taup1)
        Ocean_sfc%v_surf(i,j)  = Ocean_sfc%v_surf(i,j) + Velocity%u(ii,jj,1,2,taup1)
        Ocean_sfc%sea_lev(i,j) = Ocean_sfc%sea_lev(i,j)+ Thickness%dht(ii,jj,1,taup1)
     enddo
  enddo

  if(index_frazil > 0) then 
     do j = jsc_bnd,jec_bnd
        do i = isc_bnd,iec_bnd
           ii = i + i_shift
           jj = j + j_shift
           Ocean_sfc%frazil(i,j)  = Ocean_sfc%frazil(i,j) + T_diag(index_frazil)%field(ii,jj,1)
        enddo
     enddo
  endif 
      
end subroutine sum_ocean_sfc
! </SUBROUTINE> NAME="sum_ocean_sfc"


!#######################################################################
! <SUBROUTINE NAME="zero_ocean_sfc">
!
! <DESCRIPTION>
! Zero the elements of the Ocean_sfc derived type.  
! </DESCRIPTION>
! 
subroutine zero_ocean_sfc(Ocean_sfc)

  type(ocean_data_type), intent(inout), target :: Ocean_sfc

  integer :: i, j

  Ocean_sfc%avg_kount = 0

  do j = jsc_bnd, jec_bnd
     do i = isc_bnd, iec_bnd
        Ocean_sfc%t_surf(i,j) = 0.0
        Ocean_sfc%s_surf(i,j) = 0.0
        Ocean_sfc%u_surf(i,j) = 0.0
        Ocean_sfc%v_surf(i,j) = 0.0
        Ocean_sfc%sea_lev(i,j)= 0.0
        Ocean_sfc%frazil(i,j) = 0.0
     enddo
  enddo

end subroutine zero_ocean_sfc
! </SUBROUTINE> NAME="zero_ocean_sfc"


!#######################################################################
! <SUBROUTINE NAME="avg_ocean_sfc">
!
! <DESCRIPTION>
! Compute average of ocean surface quantities.  This is for coupling, 
! where pass time averaged information from ocean to other component
! models. Note that Ocean_sfc%frazil is NOT time averaged.  Rather, it 
! is accumulated from T_diag(index_frazil)%field in subroutine sum_ocean_sfc.
! Doing so is necessary for heat conservation between ocean and sea 
! ice systems.  Since it is not time averaged, frazil is not part of 
! this averaging subroutine.  
!
! Note that if one removes the averaging, then we take only the 
! latest values of the surface fields.  This approach has been 
! found useful to stabilize the "concurrent" coupling approach.  
!
! </DESCRIPTION>
!
subroutine avg_ocean_sfc(Time, Thickness, T_prog, Velocity, patm, Ocean_sfc)

  type(ocean_time_type), intent(in)            :: Time
  type(ocean_thickness_type), intent(in)       :: Thickness
  type(ocean_prog_tracer_type), intent(in)     :: T_prog(num_prog_tracers)
  type(ocean_velocity_type), intent(in)        :: Velocity
  real, dimension(isd:ied,jsd:jed), intent(in) :: patm
  type(ocean_data_type), intent(inout), target :: Ocean_sfc
  
  integer :: taup1, i, j, ii, jj
  real    :: divid

  taup1 = Time%taup1

  if ( Ocean_sfc%avg_kount == 0) then 
    call mpp_error (FATAL,'==>Error from ocean_sbc_mod (avg_ocean_sfc): no ocean surface quantities have been time averaged')
  endif 
                        
  divid = 1./float(Ocean_sfc%avg_kount)

  do j = jsc_bnd,jec_bnd
     do i = isc_bnd,iec_bnd
        ii = i + i_shift
        jj = j + j_shift
        if(Grd%tmask(ii,jj,1) == 1.0) then
           Ocean_sfc%t_surf(i,j)  = Ocean_sfc%t_surf(i,j)*divid + Tfreeze  !C --> K
           Ocean_sfc%s_surf(i,j)  = Ocean_sfc%s_surf(i,j)*divid
           Ocean_sfc%u_surf(i,j)  = Ocean_sfc%u_surf(i,j)*divid
           Ocean_sfc%v_surf(i,j)  = Ocean_sfc%v_surf(i,j)*divid
           Ocean_sfc%sea_lev(i,j) = Ocean_sfc%sea_lev(i,j)*divid + patm(ii,jj)/grav/rho0 !salt has patm=0.0
        endif
     enddo
  enddo

  if(.NOT. avg_sfc_temp_salt_eta) then   !replace time-averaged t,s,sealev with latest value
     do j = jsc_bnd,jec_bnd
        do i = isc_bnd,iec_bnd
           ii = i + i_shift
           jj = j + j_shift
           if(Grd%tmask(ii,jj,1) == 1.0) then
              Ocean_sfc%t_surf(i,j) = T_prog(index_temp)%field(ii,jj,1,taup1) + Tfreeze
              Ocean_sfc%s_surf(i,j) = T_prog(index_salt)%field(ii,jj,1,taup1)
              Ocean_sfc%sea_lev(i,j)= Thickness%dht(ii,jj,1,taup1) + patm(ii,jj)/grav/rho0
           endif
        enddo
     enddo
  end if

  if(.NOT. avg_sfc_velocity) then   !replace time-averaged u,v with latest value
     do j = jsc_bnd,jec_bnd
        do i = isc_bnd,iec_bnd
           ii = i + i_shift
           jj = j + j_shift
           if(Grd%tmask(ii,jj,1) == 1.0) then
              Ocean_sfc%u_surf(i,j) = Velocity%u(ii,jj,1,1,taup1)
              Ocean_sfc%v_surf(i,j) = Velocity%u(ii,jj,1,2,taup1)
           endif
        enddo
     enddo
  end if

  !set count to zero and surface quantities will be zeroed out before next sum
  Ocean_sfc%avg_kount = 0   
  
end subroutine avg_ocean_sfc
! </SUBROUTINE> NAME="avg_ocean_sfc"


!#######################################################################
! <SUBROUTINE NAME="ocean_sfc_rstrt">
!
! <DESCRIPTION>
! Save information from Ocean_sfc to restarts. Note that it is
! important in general to distinguish the time accumulated quantity
! Ocean_sfc%frazil, saved here, from the instantaneous quantity
! T_diag%frazil, which is saved in the diagnostic tracer restart file.
! </DESCRIPTION>
!
subroutine ocean_sfc_rstrt(Time, Ocean_sfc, ens_ocean)

  type(ocean_time_type) :: Time
  type(ocean_data_type), intent(in), target :: Ocean_sfc
  logical, intent(in), optional :: ens_ocean
  integer            :: yr, mon, day, hr, min, sec
  character(len=10)  :: rdte

    call get_date(Time%model_time,yr,mon,day,hr,min,sec)
    write(rdte,'(i4,3i2.2)') yr,mon,day,hr

    call write_data('IRESTART/'//rdte//'ocean_sbc.res','t_surf',Ocean_sfc%t_surf,Ocean_sfc%Domain,&
             append_pelist_name = ens_ocean)
    call write_data('IRESTART/'//rdte//'ocean_sbc.res','s_surf',Ocean_sfc%s_surf,Ocean_sfc%Domain,&
             append_pelist_name = ens_ocean)
    call write_data('IRESTART/'//rdte//'ocean_sbc.res','u_surf',Ocean_sfc%u_surf,Ocean_sfc%Domain,&
             append_pelist_name = ens_ocean)
    call write_data('IRESTART/'//rdte//'ocean_sbc.res','v_surf',Ocean_sfc%v_surf,Ocean_sfc%Domain,&
             append_pelist_name = ens_ocean)
    call write_data('IRESTART/'//rdte//'ocean_sbc.res','sea_lev',Ocean_sfc%sea_lev,Ocean_sfc%Domain,&
             append_pelist_name = ens_ocean)
    call write_data('IRESTART/'//rdte//'ocean_sbc.res','frazil',Ocean_sfc%frazil,Ocean_sfc%Domain,&
             append_pelist_name = ens_ocean)

    if(waterflux_tavg) call ocean_sbc_rstrt(Time)

end subroutine ocean_sfc_rstrt
! </SUBROUTINE> NAME="ocean_sfc_rstrt"


!#######################################################################
! <SUBROUTINE NAME="ocean_sfc_end">
!
! <DESCRIPTION>
! Save information from Ocean_sfc to restarts. Note that it is 
! important in general to distinguish the time accumulated quantity 
! Ocean_sfc%frazil, saved here, from the instantaneous quantity 
! T_diag%frazil, which is saved in the diagnostic tracer restart file.  
! </DESCRIPTION>
!
subroutine ocean_sfc_end(Ocean_sfc, ens_ocean)

  type(ocean_data_type), intent(in), target :: Ocean_sfc
  logical, intent(in), optional :: ens_ocean
  
    call write_data('RESTART/ocean_sbc.res','t_surf',Ocean_sfc%t_surf,Ocean_sfc%Domain,&
             append_pelist_name = ens_ocean)
    call write_data('RESTART/ocean_sbc.res','s_surf',Ocean_sfc%s_surf,Ocean_sfc%Domain,&
             append_pelist_name = ens_ocean)
    call write_data('RESTART/ocean_sbc.res','u_surf',Ocean_sfc%u_surf,Ocean_sfc%Domain,&
             append_pelist_name = ens_ocean)
    call write_data('RESTART/ocean_sbc.res','v_surf',Ocean_sfc%v_surf,Ocean_sfc%Domain,&
             append_pelist_name = ens_ocean)
    call write_data('RESTART/ocean_sbc.res','sea_lev',Ocean_sfc%sea_lev,Ocean_sfc%Domain,&
             append_pelist_name = ens_ocean)
    call write_data('RESTART/ocean_sbc.res','frazil',Ocean_sfc%frazil,Ocean_sfc%Domain,&
             append_pelist_name = ens_ocean) 

    if(waterflux_tavg) call ocean_sbc_end

end subroutine ocean_sfc_end
! </SUBROUTINE> NAME="ocean_sfc_end"


!#######################################################################
! <SUBROUTINE NAME="get_ocean_sbc">
!
! <DESCRIPTION>
! Subroutine to get the surface fluxes passed into the ocean from 
! other component models. 
! </DESCRIPTION>
!
subroutine get_ocean_sbc(Time, Ice_ocean_boundary, Ext_mode, T_prog, Velocity, pme, upme, river, uriver, swflx, patm )

  type(ocean_time_type), intent(in)                 :: Time 
  type(ice_ocean_boundary_type), intent(in)         :: Ice_ocean_boundary
  type(ocean_external_mode_type), intent(in)        :: Ext_mode
  type(ocean_prog_tracer_type), intent(inout)       :: T_prog(num_prog_tracers)
  type(ocean_velocity_type), intent(inout)          :: Velocity
  real, dimension(isd:ied,jsd:jed), intent(inout)   :: pme, river, swflx, patm
  real, dimension(isd:ied,jsd:jed,2), intent(inout) :: upme, uriver
  real, dimension(isd:ied,jsd:jed)                  :: pme_river

  real    :: tmp_x, tmp_y
  real    :: pme_river_total, var
  integer :: tau, n, i, j, ii, jj
  logical :: used

  tau   = Time%tau

! i and j momentum flux per mass crossing ocean surface (m^2/sec^2)

  do j = jsc_bnd,jec_bnd
     do i = isc_bnd,iec_bnd
        ii = i + i_shift
        jj = j + j_shift
        Velocity%smf(ii,jj,1) = Ice_ocean_boundary%u_flux(i,j)*rho0r*Grd%umask(ii,jj,1)
        Velocity%smf(ii,jj,2) = Ice_ocean_boundary%v_flux(i,j)*rho0r*Grd%umask(ii,jj,1)
     enddo
  enddo

  if (rotate_winds) then
     do j=jsc,jec
        do i=isc,iec
           tmp_x =  Grd%cos_rot(i,j)*Velocity%smf(i,j,1) + Grd%sin_rot(i,j)*Velocity%smf(i,j,2)
           tmp_y = -Grd%sin_rot(i,j)*Velocity%smf(i,j,1) + Grd%cos_rot(i,j)*Velocity%smf(i,j,2)
           Velocity%smf(i,j,1) = tmp_x
           Velocity%smf(i,j,2) = tmp_y
        enddo
     enddo
  endif

  call mpp_update_domains(Velocity%smf(:,:,1),Velocity%smf(:,:,2),Dom%domain2d,gridtype=BGRID_NE)

  ! Set temperature for evap and precip to ocean surface value. 
  ! For all tracers, by default keep them set to zero.
  do j=jsc,jec
     do i=isc,iec
       T_prog(index_temp)%tpme(i,j) = T_prog(index_temp)%field(i,j,1,tau)
     enddo
  enddo

  ! Set temperature for river water to ocean surface value, but no
  ! less than zero degrees C. For other tracers, by default keep them
  ! set to zero concentration in rivers.
  do j=jsc,jec
     do i=isc,iec
        T_prog(index_temp)%triver(i,j) = max(0.0,T_prog(index_temp)%field(i,j,1,tau))
     enddo
  enddo


  if (use_waterflux) then

      ! freshwater flux in kg/(m^2 sec) from liquid, frozen, and q_flux (evaporation)  
      ! rho0r converts mass flux to (m/s) 
      if(waterflux_tavg) then 
         do j = jsc_bnd,jec_bnd
            do i = isc_bnd,iec_bnd
               ii = i + i_shift
               jj = j + j_shift
               pme(ii,jj) = 0.5*pme_taum1(ii,jj)
               pme_taum1(ii,jj) = (Ice_ocean_boundary%lprec(i,j) + Ice_ocean_boundary%fprec(i,j)-&
                                   Ice_ocean_boundary%q_flux(i,j))*rho0r*Grd%tmask(ii,jj,1) 
               pme(ii,jj) = pme(ii,jj) + 0.5*pme_taum1(ii,jj)
            enddo
         enddo
      else 
         do j = jsc_bnd,jec_bnd
            do i = isc_bnd,iec_bnd
               ii = i + i_shift
               jj = j + j_shift
                  pme(ii,jj) = (Ice_ocean_boundary%lprec(i,j) + Ice_ocean_boundary%fprec(i,j)-&
                                Ice_ocean_boundary%q_flux(i,j))*rho0r*Grd%tmask(ii,jj,1) 
            enddo
         enddo
      endif 

      ! freshwater flux in kg/(m^2 sec) from rivers and calving land glaciers  
      ! rho0r converts mass flux to (m/s) 
      if(waterflux_tavg) then 
         do j = jsc_bnd,jec_bnd
            do i = isc_bnd,iec_bnd
               ii = i + i_shift
               jj = j + j_shift
                  river(ii,jj) = 0.5*river_taum1(ii,jj)
                  river_taum1(ii,jj) = (Ice_ocean_boundary%runoff(i,j)+Ice_ocean_boundary%calving(i,j)) &
                                       *rho0r*Grd%tmask(ii,jj,1)
                  river(ii,jj) = river(ii,jj) + 0.5*river_taum1(ii,jj)
            enddo
         enddo
      else 
         do j = jsc_bnd,jec_bnd
            do i = isc_bnd,iec_bnd
               ii = i + i_shift
               jj = j + j_shift
               river(ii,jj) = (Ice_ocean_boundary%runoff(i,j)+Ice_ocean_boundary%calving(i,j)) &
                              *rho0r*Grd%tmask(ii,jj,1)
            enddo
         enddo

      endif 

      if(riverspread) call spread_river_horz(river)

      ! output river information (m/s)
      if (id_river > 0) used =  send_data(id_river, river(:,:),&
                                Time%model_time, rmask=Grd%tmask(:,:,1), &
                                is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)  

      ! set riverdiffuse to determine where to enhance diff_cbt inside ocean_rivermix_mod
      do n=1,num_prog_tracers
         do j = jsc, jec
            do i = isc, iec
               T_prog(n)%riverdiffuse(i,j) = river(i,j) 
            enddo
         enddo
      enddo

      ! produce a zero area average of pme + river 
      pme_river = 0.0 
      if(zero_net_water_coupler) then 
         do j = jsc, jec
            do i = isc, iec
               pme_river(i,j) = pme(i,j) + river(i,j)
            enddo
         enddo
         pme_river_total            = mpp_global_sum(Dom%domain2d,pme_river(:,:)*Grd%dat(:,:)*Grd%tmask(:,:,1),&
                                                     BITWISE_EXACT_SUM)/Grd%tcellsurf
         do j = jsc, jec
            do i = isc, iec
               pme(i,j)       = pme(i,j) - pme_river_total*Grd%tmask(i,j,1)
            enddo
         enddo
      endif 

      ! allow for nonzero salt flux from, say, ice melt.  
      ! the flux from ice-melt is in units of (kg/m^2/sec)
      ! and so should be converted to model units of 
      ! [salinity(psu)*m/sec], which requires multiplication by 
      ! rho0r*1000.0
      do j = jsc_bnd,jec_bnd
         do i = isc_bnd,iec_bnd
            ii = i + i_shift
            jj = j + j_shift      
               T_prog(index_salt)%stf(ii,jj) = -Ice_ocean_boundary%salt_flux(i,j)*rho0r*1000.0
         enddo
      enddo
  else
     
      ! freshwater mass flux (kg/m^2/sec) from rivers and calving land glaciers
      do j = jsc_bnd,jec_bnd
         do i = isc_bnd,iec_bnd
            ii = i + i_shift
            jj = j + j_shift  
            river(ii,jj) = (Ice_ocean_boundary%runoff(i,j)+Ice_ocean_boundary%calving(i,j)) &
                           *Grd%tmask(ii,jj,1)
         enddo
      enddo
      if(riverspread) call spread_river_horz(river)

      ! convert freshwater mass fluxes (kg/m^2/sec) into virtual salt flux [salinity(psu)*m/sec]
      ! also convert salt flux from ice to [salinity(psu)*m/sec] using rho0r*1000.0 
         do j = jsc_bnd, jec_bnd
            do i = isc_bnd, iec_bnd
               ii = i + i_shift
               jj = j + j_shift
               T_prog(index_salt)%stf(ii,jj) = -Ice_ocean_boundary%salt_flux(i,j)*rho0r*1000.0 -&
                                             (Ice_ocean_boundary%lprec(i,j) +                 &
                                             Ice_ocean_boundary%fprec(i,j) + river(ii,jj)   - &
                                             Ice_ocean_boundary%q_flux(i,j)) *salinity_ref*rho0r*Grd%tmask(ii,jj,1) 
         enddo
      enddo
      ! set the riverdiffuse "mask" to determine where to enhance diff_cbt
      do n=1,num_prog_tracers
         do j = jsc,jec
            do i = isc, iec
               T_prog(n)%riverdiffuse(i,j) = river(i,j) 
            enddo
         enddo
      enddo

      ! output river information (m/s) before set river array to zero 
      if (id_river > 0) used =  send_data(id_river, rho0r*river(:,:),&
                                Time%model_time, rmask=Grd%tmask(:,:,1), &
                                is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)  
 
      ! set pme and river to zero for case with use_water false 
      pme   = 0.0
      river = 0.0

  endif

  call mpp_update_domains(pme(:,:), Dom%domain2d)
  call mpp_update_domains(river(:,:), Dom%domain2d)
  
  do j = jsc_bnd, jec_bnd
     do i = isc_bnd, iec_bnd
        ii = i + i_shift
        jj = j + j_shift  
           T_prog(index_temp)%stf(ii,jj) = (Ice_ocean_boundary%sw_flux(i,j) +     &
                                            Ice_ocean_boundary%lw_flux(i,j) -     &
                                           (Ice_ocean_boundary%fprec(i,j) +       &
                                            Ice_ocean_boundary%calving(i,j))*hlf -&
                                            Ice_ocean_boundary%t_flux(i,j)  -     &
                                            Ice_ocean_boundary%q_flux(i,j)*hlv    &
                                           )/rho_cp*Grd%tmask(ii,jj,1) 
     enddo
  enddo

  ! output total surface tracer flux for temp and salt, 
  ! and temp/salt flux associated with rivers        

  if (id_stf(index_temp) > 0) then
      used = send_data(id_stf(index_temp), T_prog(index_temp)%stf(:,:)*rho_cp,&
             Time%model_time, rmask=Grd%tmask(:,:,1), &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif
  if (id_stf_river(index_temp) > 0) then
      used = send_data(id_stf_river(index_temp), river(:,:)*T_prog(index_temp)%triver(:,:)*rho_cp,&
             Time%model_time, rmask=Grd%tmask(:,:,1), &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif

  if (id_stf(index_salt) > 0) then
      used = send_data(id_stf(index_salt), T_prog(index_salt)%stf(:,:)*rho0*0.001,&
             Time%model_time, rmask=Grd%tmask(:,:,1), &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif
  if (id_stf_river(index_salt) > 0) then
      used = send_data(id_stf_river(index_salt), river(:,:)*T_prog(index_salt)%triver(:,:)*rho0*0.001,&
             Time%model_time, rmask=Grd%tmask(:,:,1), &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif
  
  if (id_pme > 0) then
      used = send_data(id_pme, pme(:,:), &
             Time%model_time, rmask=Grd%tmask(:,:,1), &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif


  ! set velocity of pme and river water to that of upper ocean cell.   
  ! generalizations may be suitable with refined component models.   
  do n = 1, size(upme,3)
     do j = jsc, jec
        do i = isc, iec
           upme(i,j,n)   = Velocity%u(i,j,1,n,tau)   ! velocity of precip-evap water
           uriver(i,j,n) = Velocity%u(i,j,1,n,tau)   ! velocity of river water 
        enddo
     enddo
  enddo
  ! apply ice pressure to ocean only when there is mass 
  ! transfer allowed between ocean and ice via fresh water. 
  ! if do not use fresh water, then cannot allow ocean to feel weight of ice.  
  if(use_waterflux) then 
     var = grav*rho0*max_ice_thickness
     do j = jsc_bnd,jec_bnd
        do i = isc_bnd,iec_bnd
           ii = i + i_shift
           jj = j + j_shift
           patm(ii,jj) = min(Ice_ocean_boundary%p(i,j),var)
        enddo
     enddo
     call mpp_update_domains(patm(:,:), Dom%domain2d)
  endif 

  swflx(isc:iec,jsc:jec) = Ice_ocean_boundary%sw_flux(isc_bnd:iec_bnd,jsc_bnd:jec_bnd)/rho_cp*&
       Grd%tmask(isc:iec,jsc:jec,1)

  if (id_swflx > 0) used =  send_data(id_swflx, rho_cp*swflx(:,:),&
                            Time%model_time, rmask=Grd%tmask(:,:,1), &
                            is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

  if (id_tau_x > 0) used =  send_data(id_tau_x, Velocity%smf(:,:,1)*rho0,&
                            Time%model_time, rmask=Grd%umask(:,:,1), &
                            is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

  if (id_tau_y > 0) used =  send_data(id_tau_y, Velocity%smf(:,:,2)*rho0,&
                            Time%model_time, rmask=Grd%umask(:,:,1), &
                            is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)  

end subroutine get_ocean_sbc
! </SUBROUTINE> NAME="get_ocean_sbc"


!#######################################################################
! <SUBROUTINE NAME="flux_adjust">
!
! <DESCRIPTION>
! Subroutine to compute the surface fluxes derived from a 
! restoring condition. We use a convention whereby a positive 
! flux enters the ocean:  (+) down convention. 
!
! Note that the term "flux adjust" should not be confused with 
! a similar term used in coupled modeling, where the atmosphere/ocean
! fluxes realized in a coupled atmosphere/ocean model are supplemented
! by fluxes diagnosed from an uncoupled model.  
! </DESCRIPTION>
!
subroutine flux_adjust(Time, T_prog, Velocity, pme)

  type(ocean_time_type), intent(in)                                        :: Time
  type(ocean_prog_tracer_type), dimension(num_prog_tracers), intent(inout) :: T_prog
  type(ocean_velocity_type), intent(inout)                                 :: Velocity
  real, dimension(isd:ied,jsd:jed), intent(inout)                          :: pme 

  real, dimension(isd:ied,jsd:jed) :: pme_adj, flx_adj, open_ocean_mask
  real                             :: pme_adj_total, flx_adj_total 
  integer                          :: i, j, n, tau, taum1
  logical                          :: used
  
  tau      = Time%tau
  taum1    = Time%taum1  
  flx_adj  = 0.0 
  pme_adj  = 0.0
  open_ocean_mask(isd:ied,jsd:jed) = Grd%tmask(isd:ied,jsd:jed,1)

! add restoring to fluxes from coupled model or
! some other type of flux correction
! NOTE: (+) down convention

  if (id_restore(index_salt) > 0 .and. salt_damp_factor > 0.0) then
      call time_interp_external(id_restore(index_salt), Time%model_time, data)

      ! use near-frazil condition as a proxy for where sea-ice is present 
      if(.not. salt_restore_under_ice) then 
          do j=jsc,jec
             do i=isc,iec
                if(Grd%tmask(i,j,1) == 1.0) then 
                    if(T_prog(index_temp)%field(i,j,1,tau) <= -0.0539*T_prog(index_salt)%field(i,j,1,tau)) then 
                        open_ocean_mask(i,j) = 0.0
                    endif
                endif
             enddo
          enddo
      endif

      if (use_waterflux) then

          ! put all fresh water adjustment into pme_adj (no river_adj has been coded) 
          do j = jsc,jec
             do i = isc,iec
                pme_adj(i,j) = salt_damp_factor*open_ocean_mask(i,j)*Grd%tmask(i,j,1)* &
                               (T_prog(index_salt)%field(i,j,1,taum1) -                &
                               data(i,j))/(T_prog(index_salt)%field(i,j,1,taum1)+epsln)
             enddo
          enddo     

          ! produce a zero area average so there is no net input 
          ! of fresh water to the ocean associated with the restoring 
          if(zero_net_water_restore) then 
             pme_adj_total = mpp_global_sum(Dom%domain2d,pme_adj(:,:)*Grd%dat(:,:)*Grd%tmask(:,:,1), &
                             BITWISE_EXACT_SUM)/Grd%tcellsurf
             pme_adj(isc:iec,jsc:jec) = pme_adj(isc:iec,jsc:jec) - pme_adj_total*Grd%tmask(isc:iec,jsc:jec,1)
          endif 

          ! add pme_adj to pme
          pme(isc:iec,jsc:jec) = pme(isc:iec,jsc:jec) + pme_adj(isc:iec,jsc:jec)
          call mpp_update_domains(pme(:,:), Dom%domain2d)

          flx_adj=0.0
      else

          pme_adj = 0.0
          do j = jsc,jec
             do i = isc,iec
                flx_adj(i,j) = salt_damp_factor*open_ocean_mask(i,j)*Grd%tmask(i,j,1)* &
                               (data(i,j) - T_prog(index_salt)%field(i,j,1,taum1))
             enddo
          enddo
          ! produce a zero area average so there is no net input
          ! of salt to the ocean associated with the restoring 
          if(zero_net_salt_restore) then 
             flx_adj_total =  &
               mpp_global_sum(Dom%domain2d,flx_adj(:,:)*Grd%dat(:,:)*Grd%tmask(:,:,1), BITWISE_EXACT_SUM)/Grd%tcellsurf
             flx_adj(isc:iec,jsc:jec) = flx_adj(isc:iec,jsc:jec) - flx_adj_total*Grd%tmask(isc:iec,jsc:jec,1)
          endif 

          ! add flx_adj to stf
          do j = jsc,jec
             do i = isc,iec
                T_prog(index_salt)%stf(i,j) = T_prog(index_salt)%stf(i,j) + flx_adj(i,j)
             enddo
          enddo

      endif
  else
      flx_adj = 0.0
      pme_adj = 0.
  endif

  ! send diagnostics 

  if (id_stf_adj(index_salt) > 0) then
      used = send_data(id_stf_adj(index_salt), flx_adj(:,:)*rho0*0.001,&
             Time%model_time, rmask=Grd%tmask(:,:,1), &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif
  if (id_stf_total(index_salt) > 0) then
      used = send_data(id_stf_total(index_salt), T_prog(index_salt)%stf(:,:)*rho0*0.001,&
             Time%model_time, rmask=Grd%tmask(:,:,1), &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif
  if (id_pme_adj > 0) then
      used = send_data(id_pme_adj, pme_adj(:,:), &
             Time%model_time, rmask=Grd%tmask(:,:,1), &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif
  if (id_pme_total > 0) then
      used = send_data(id_pme_total, pme(:,:), &
             Time%model_time, rmask=Grd%tmask(:,:,1), &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif
  if (id_ice_mask > 0) then
      used = send_data(id_ice_mask, (1.0-open_ocean_mask(:,:))*Grd%tmask(:,:,1),&
             Time%model_time,rmask=Grd%tmask(:,:,1), &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif
  if (id_open_ocean_mask > 0) then
      used = send_data(id_open_ocean_mask, open_ocean_mask(:,:),&
             Time%model_time,rmask=Grd%tmask(:,:,1), &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif


  ! temp restoring   
  if (id_restore(index_temp) > 0 .and. temp_damp_factor > 0.0) then
     call time_interp_external(id_restore(index_temp), Time%model_time, data)
     do j = jsc,jec
        do i = isc,iec
           flx_adj(i,j) = temp_damp_factor*Grd%tmask(i,j,1)*(data(i,j)-T_prog(index_temp)%field(i,j,1,taum1))
           T_prog(index_temp)%stf(i,j) = T_prog(index_temp)%stf(i,j) + flx_adj(i,j)
        enddo
     enddo
  else
     flx_adj = 0.0
  endif

  if (id_stf_adj(index_temp) > 0) then
      used = send_data(id_stf_adj(index_temp), flx_adj(:,:)*rho_cp, &
             Time%model_time,rmask=Grd%tmask(:,:,1), &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif
  
  if (id_stf_total(index_temp) > 0) then
      used = send_data(id_stf_total(index_temp), T_prog(index_temp)%stf(:,:)*rho_cp, &
             Time%model_time, rmask=Grd%tmask(:,:,1),&
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif


  return


end subroutine flux_adjust
! </SUBROUTINE> NAME="flux_adjust"


!#######################################################################
! <SUBROUTINE NAME="ocean_sbc_rstrt">
!
! <DESCRIPTION>
! Save pme_taum1 and river_taum1 to restart.
! </DESCRIPTION>
!
subroutine ocean_sbc_rstrt(Time, ens_ocean)

    type(ocean_time_type) :: Time
    logical, intent(in), optional :: ens_ocean
    integer            :: yr, mon, day, hr, min, sec
    character(len=10)  :: rdte

    call get_date(Time%model_time,yr,mon,day,hr,min,sec)
    write(rdte,'(i4,3i2.2)') yr,mon,day,hr

    call write_data('IRESTART/'//rdte//'ocean_waterflux.res','pme_taum1',pme_taum1,Dom%domain2d,&
                     append_pelist_name = ens_ocean)
    call write_data('IRESTART/'//rdte//'ocean_waterflux.res','river_taum1',river_taum1,Dom%domain2d,&
                     append_pelist_name = ens_ocean)

end subroutine ocean_sbc_rstrt
! </SUBROUTINE> NAME="ocean_sbc_rstrt"


!#######################################################################
! <SUBROUTINE NAME="ocean_sbc_end">
!
! <DESCRIPTION>
! Save pme_taum1 and river_taum1 to restart.  
! </DESCRIPTION>
!
subroutine ocean_sbc_end(ens_ocean)
    logical, intent(in), optional :: ens_ocean

    call write_data('RESTART/ocean_waterflux.res','pme_taum1',pme_taum1,Dom%domain2d,&
                     append_pelist_name = ens_ocean) 
    call write_data('RESTART/ocean_waterflux.res','river_taum1',river_taum1,Dom%domain2d,&
                     append_pelist_name = ens_ocean) 

end subroutine ocean_sbc_end
! </SUBROUTINE> NAME="ocean_sbc_end"



end module ocean_sbc_mod
