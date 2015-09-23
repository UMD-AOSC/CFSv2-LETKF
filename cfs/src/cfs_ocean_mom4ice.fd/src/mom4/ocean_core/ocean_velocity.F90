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
module ocean_velocity_mod
!
!<CONTACT EMAIL="Ronald.Pacanowski@noaa.gov"> R.C. Pacanowski
!</CONTACT>
!
!<CONTACT EMAIL="Tony.Rosati@noaa.gov"> A. Rosati
!</CONTACT>
!
! <REVIEWER EMAIL="Matthew.Harrison@noaa.gov">
! M.J. Harrison 
! </REVIEWER>
!
! <REVIEWER EMAIL="Stephen.Griffies@noaa.gov">
! S.M. Griffies 
! </REVIEWER>
!
!<OVERVIEW>
! Time step velocity 
!</OVERVIEW>
!
!<DESCRIPTION>
! This module steps the velocity field forward in time using a 
! leap-frog time stepping scheme.
!</DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! Durran, Numerical Methods for Wave Equations in Geophysical
! Fluid Dynamics (1999).
! </REFERENCE>
!
! <REFERENCE>
! R.C. Pacanowski and S.M. Griffies, The MOM3 Manual (1999).
! NOAA/Geophysical Fluid Dynamics Laboratory
! </REFERENCE>
!
! <REFERENCE>
! S.M. Griffies, M.J. Harrison, R.C. Pacanowski, and 
! A. Rosati, A Technical Guide to MOM4 (2004).
! NOAA/Geophysical Fluid Dynamics Laboratory
! </REFERENCE>
!
! <REFERENCE>
! S.M. Griffies, Fundamentals of Ocean Climate Models (2004).
! Princeton University Press.
! </REFERENCE>
!
! <REFERENCE>
! S.M. Griffies, R.C. Pacanowski, R.M. Schmidt, and V. Balaji
! Tracer Conservation with an Explicit Free Surface Method for 
! Z-coordinate Ocean Models
! Monthly Weather Review (2001) vol 129 pages 1081--1098
! </REFERENCE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_velocity_nml">
!  <DATA NAME="zero_tendency" TYPE="logical">
!  For debugging. Will freeze the baroclinic velocity  fields. 
!  </DATA> 
!  <DATA NAME="zero_tendency_explicit" TYPE="logical">
!  For debugging. Will not use explicit part of the tendency. 
!  </DATA> 
!  <DATA NAME="zero_tendency_implicit" TYPE="logical">
!  For debugging. Will not use implicit part of the tendency. 
!  </DATA> 
!  <DATA NAME="truncate_velocity" TYPE="logical">
!  Truncate the velocity to a maximum value.  Useful for cases where
!  the initial spin-up initiates spuriously large model velocities that would
!  otherwise cause the model to blow-up. 
!  </DATA> 
!  <DATA NAME="vel_max" UNITS="meter/sec" TYPE="real">
!  Truncation velocity
!  </DATA> 
!  <DATA NAME="max_cgint" TYPE="real">
!  Maximum internal gravity wave speed--used for diagnosing conservative
!  estimate of stable time steps.  
!  </DATA>
!
!  <DATA NAME="adams_bashforth_epsilon" UNITS="dimensionless" TYPE="real">
!  Dimensionless parameter for 2nd order Adams-Bashforth implementation of 
!  velocity advection.  Values between 0.5 and 1.0 are recommended.  
!  Value of 0.5 leads to second order accurate, but it is formally 
!  weakly unstable (Durran, Section 2.3.4).
!  </DATA> 
!  <DATA NAME="adams_bashforth_third" TYPE="logicall">
!  For a third order treatment of the velocity advection.  
!  This is stable and so needs no temporal dissipation 
!  (Section 2.3.6 of Durran).  This is the model default. 
!  </DATA> 

!  <DATA NAME="truncate_verbose" TYPE="logical">
!  For verbose printout 
!  </DATA> 
!</NAMELIST>

use constants_mod,    only: epsln, rho0, grav, rho0r, c2dbars
use diag_manager_mod, only: register_diag_field, send_data
use fms_mod,          only: write_version_number, check_nml_error, close_file, open_namelist_file
use fms_mod,          only: read_data, write_data, file_exist, FATAL, WARNING, NOTE
use fms_io_mod,       only: field_size
use mpp_domains_mod,  only: mpp_update_domains, mpp_global_sum, BGRID_NE, BITWISE_EXACT_SUM
use mpp_mod,          only: mpp_chksum, mpp_error, mpp_pe, mpp_min, mpp_max, mpp_sum
use mpp_mod,          only: mpp_broadcast, stdout, stdlog
use time_manager_mod, only: get_date

use ocean_bih_friction_mod,    only: horz_bih_friction
use ocean_coriolis_mod,        only: coriolis
use ocean_domains_mod,         only: get_local_indices
use ocean_lap_friction_mod,    only: horz_lap_friction
use ocean_obc_mod,             only: ocean_obc_update_boundary
use ocean_operators_mod,       only: DIV_UD, REMAP_BT_TO_BU, FAX, FAY, FDX_NT, FDY_ET, BDX_ET, BDY_NT
use ocean_pressure_mod,        only: pressure_gradient
use ocean_types_mod,           only: ocean_domain_type, ocean_grid_type, ocean_thickness_type
use ocean_types_mod,           only: ocean_time_type, ocean_time_steps_type
use ocean_types_mod,           only: ocean_density_type, ocean_thickness_type, ocean_velocity_type
use ocean_types_mod,           only: ocean_adv_vel_type, ocean_external_mode_type
use ocean_types_mod,           only: TWO_LEVEL, THREE_LEVEL
use ocean_util_mod,            only: invtri, write_timestamp
use ocean_velocity_advect_mod, only: horz_advection_of_velocity, vertical_advection_of_velocity
use ocean_velocity_diag_mod,   only: kinetic_energy, potential_energy 
use ocean_vert_mix_mod,        only: vert_frict
use ocean_workspace_mod,       only: wrk1, wrk2, wrk1_v  

implicit none

private

! for diagnostics 
integer :: id_u=-1
integer :: id_v=-1
integer :: id_usurf=-1
integer :: id_vsurf=-1
integer :: id_accel(2)=-1
integer :: id_speed=-1
integer :: id_vfrict_impl(2)=-1
integer :: id_ucori_impl=-1
integer :: id_vcori_impl=-1
integer :: id_advectionu=-1
integer :: id_advectionv=-1
logical :: used

integer :: unit=6

character(len=128) :: &
     version='$Id$'
character (len=128) :: tagname = &
     '$Name$'

logical :: debug_velocity = .false.
logical :: have_obc       = .false.

#include <ocean_memory.h>

#ifdef STATIC_MEMORY
real, dimension(isd:ied,jsd:jed,nk,2) :: vfrict_impl
#else
real, allocatable, dimension(:,:,:,:) :: vfrict_impl
#endif

real    :: aidif
real    :: dtime 
real    :: dtimer 
integer :: time_index 
integer :: tendency

type(ocean_domain_type), pointer :: Dom =>NULL()
type(ocean_grid_type), pointer   :: Grd =>NULL()

public ocean_velocity_init
public ocean_velocity_rstrt
public ocean_velocity_end
public ocean_explicit_accel_a
public ocean_explicit_accel_b
public update_ocean_velocity
public ocean_implicit_friction
public ocean_implicit_coriolis
public energy_analysis 

private ocean_velocity_chksum
private check_gravity_wave_cfl

! Adams-Bashforth treatment of velocity advection 
logical :: adams_bashforth_third   = .true.
real    :: adams_bashforth_epsilon = 0.6
real    :: abtau_m0 =  23.0/12.0
real    :: abtau_m1 = -16.0/12.0
real    :: abtau_m2 =  5.0/12.0

logical :: module_is_initialized  = .false.
logical :: zero_tendency          = .false.
logical :: zero_tendency_explicit = .false.
logical :: zero_tendency_implicit = .false.
logical :: truncate_velocity      = .false.
logical :: truncate_verbose       = .false.
real    :: vel_max                = 1.0

! for CFL check on gravity waves 
real    :: max_cgint = 2.0   ! speed (m/s) used to compute CFL check for internal gravity waves  
real    :: max_dt_for_cgint  ! maximum time step (s) allowed by CFL to resolve internal gravity waves 


namelist /ocean_velocity_nml/ zero_tendency, zero_tendency_explicit, zero_tendency_implicit, &
                              truncate_velocity, truncate_verbose, vel_max, max_cgint,       &
                              adams_bashforth_third, adams_bashforth_epsilon

contains



!#######################################################################
! <SUBROUTINE NAME="ocean_velocity_init">
!
! <DESCRIPTION>
! Initialize terms for the velocity equation. 
! </DESCRIPTION>
!
subroutine ocean_velocity_init (Grid, Domain, Time, Time_steps, Velocity, obc, debug, ens_ocean)

  type(ocean_grid_type), target, intent(in)   :: Grid
  type(ocean_domain_type), target, intent(in) :: Domain
  type(ocean_time_type), intent(in)           :: Time
  type(ocean_time_steps_type), intent(in)     :: Time_steps 
  type(ocean_velocity_type), intent(inout)    :: Velocity
  logical, intent(in)            :: obc
  logical, intent(in), optional  :: debug
  logical, intent(in), optional  :: ens_ocean

  integer               :: ioun, io_status, ierr
  integer               :: m, n
  integer               :: taum1, tau, taup1
  integer               :: tau_m0, tau_m1
  integer, dimension(4) :: siz 
  character(len=128)    :: filename

  if ( module_is_initialized ) then 
    call mpp_error(FATAL, '==>Error from ocean_velocity_mod (ocean_velocity_init): module already initialized')
  endif 

  module_is_initialized = .TRUE.

  call write_version_number( version, tagname )

  have_obc = obc
  tendency = Time_steps%tendency
  aidif    = Time_steps%aidif 
  dtime    = Time_steps%dtime_u
  dtimer   = 1.0/dtime

  ! provide for namelist over-ride of defaults 
  ioun = open_namelist_file()
  read  (ioun, ocean_velocity_nml,iostat=io_status)
  write (stdout(),'(/)')
  write (stdout(), ocean_velocity_nml)  
  write (stdlog(), ocean_velocity_nml)
  ierr = check_nml_error(io_status, 'ocean_velocity_nml')
  call close_file(ioun)

  if (PRESENT(debug) .and. .not. debug_velocity) then 
    debug_velocity = debug
  endif 

  Dom => Domain
  Grd => Grid

#ifndef STATIC_MEMORY
  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk=Grid%nk
  allocate (vfrict_impl(isd:ied,jsd:jed,nk,2))
  allocate (Velocity%u(isd:ied,jsd:jed,nk,2,3))
  allocate (Velocity%smf(isd:ied,jsd:jed,2))
  allocate (Velocity%bmf(isd:ied,jsd:jed,2))
  allocate (Velocity%advection(isd:ied,jsd:jed,nk,2,3))
#endif

  Velocity%smf       = 0.0
  Velocity%bmf       = 0.0
  Velocity%advection = 0.0

  ! register fields for diagnostic output

  id_u    = register_diag_field ('ocean_model', 'u', Grd%vel_axes_uv(1:3), Time%model_time, &
     'zonal current', 'm/sec', missing_value=-10.0, range=(/-10.0,10.0/))
  id_v    = register_diag_field ('ocean_model', 'v', Grd%vel_axes_uv(1:3), Time%model_time, &
     'meridional current', 'm/sec', missing_value=-10.0, range=(/-10.0,10.0/))
  id_usurf  = register_diag_field ('ocean_model', 'usurf', Grd%vel_axes_uv(1:2), Time%model_time, &
     'zonal surface current', 'm/sec', missing_value=-10.0, range=(/-10.0,10.0/))
  id_vsurf  = register_diag_field ('ocean_model', 'vsurf', Grd%vel_axes_uv(1:2), Time%model_time, &
     'meridional surface current', 'm/sec', missing_value=-10.0, range=(/-10.0,10.0/))
  id_speed = register_diag_field ('ocean_model', 'speed', Grd%vel_axes_uv(1:3), &
     Time%model_time, 'speed of horizontal current', 'm/sec', missing_value=-10.0, range=(/-10.0,10.0/))
  id_accel(1) = register_diag_field ('ocean_model', 'accel_u', Grd%vel_axes_uv(1:3), &
     Time%model_time, 'baroclinic forcing of u', 'm/s^2', missing_value=-1e10, range=(/-1e10,1e10/))
  id_accel(2) = register_diag_field ('ocean_model', 'accel_v', Grd%vel_axes_uv(1:3), Time%model_time, &
     'baroclinic forcing of v', 'm/s^2', missing_value=-1e10, range=(/-1e10,1e10/))
  id_vfrict_impl(1) = register_diag_field ('ocean_model', 'vfrict_impl_u', Grd%vel_axes_uv(1:3), Time%model_time, &
     'implicit vertical u-mixing', 'm/s^2', missing_value=-1e10, range=(/-1e10,1e10/))
  id_vfrict_impl(2) = register_diag_field ('ocean_model', 'vfrict_impl_v', Grd%vel_axes_uv(1:3), Time%model_time, &
     'implicit vertical v-mixing', 'm/s^2', missing_value=-1e10, range=(/-1e10,1e10/))
  id_ucori_impl = register_diag_field ('ocean_model', 'ucori_impl', Grd%vel_axes_uv(1:3), Time%model_time, &
     'implicit Coriolis force in i-direction', 'm/s^2', missing_value=-1e10, range=(/-1e10,1e10/))
  id_vcori_impl = register_diag_field ('ocean_model', 'vcori_impl', Grd%vel_axes_uv(1:3), Time%model_time, &
     'implicit Coriolis force in j-direction', 'm/s^2', missing_value=-1e10, range=(/-1e10,1e10/))
  id_advectionu = register_diag_field ('ocean_model', 'advectionu', Grd%vel_axes_uv(1:3), Time%model_time, &
     'thickness weighted advection of u', 'm^2/s^2', missing_value=-1e10, range=(/-1e10,1e10/))
  id_advectionv = register_diag_field ('ocean_model', 'advectionv', Grd%vel_axes_uv(1:3), Time%model_time, &
     'thickness weighted advection of v', 'm^2/s^2', missing_value=-1e10, range=(/-1e10,1e10/))

  if (truncate_velocity) then
      write (stdout(),'(/a)')     '=>Warning: truncate_velocity=.true. Model will truncate baroclinc velocity'
      write (stdout(),'(a,f5.3)') '           so that max velocity magnitude in any direction is (m/s) ',vel_max 
      write (stdout(),'(a)')      '           This option may be useful during early spin-up phase of experiment.'
      write (stdout(),'(a/)')     '           Note that this option will not reproduce across restarts.  '
  endif

  if(tendency==TWO_LEVEL) then 
     write(stdout(),'(/a)')'==>Note from ocean_velocity_mod: use of twolevel time_tendency' 
     write(stdout(),'(a)')'   necessitates an Adams-Bashforth treatment of velocity advection.'

     if(adams_bashforth_third) then 
       write(stdout(),'(/a)')'   Using 3rd order Adams-Bashforth for velocity advection.'
       write(stdout(),'(a/)')'   This is the model default. '
       abtau_m0= 23.0/12.0
       abtau_m1=-16.0/12.0
       abtau_m2= 5.0/12.0
     else 
       write(stdout(),'(/a,f8.3/)')'   Using 2nd order Adams-Bashforth with dimensionless AB parameter = ',adams_bashforth_epsilon 
       abtau_m0= 1.0 + adams_bashforth_epsilon
       abtau_m1= 0.0 - adams_bashforth_epsilon
       abtau_m2= 0.0
       if(adams_bashforth_epsilon < 0.5 .or. adams_bashforth_epsilon > 1.0) then 
         call mpp_error(FATAL,'==>Error from ocean_velocity_mod: adams_bashforth parameter must be between 0.5 and 1.0')
       endif 
     endif  

  endif 

  taum1 = Time%taum1
  tau   = Time%tau
  taup1 = Time%taup1

  tau_m1 = Time%tau_m1
  tau_m0 = Time%tau_m0
  
  do m=1,3 ; do n=1,2  
     Velocity%u(:,:,:,n,m) = epsln*Grd%umask(:,:,:)
  enddo; enddo


  filename = 'INPUT/ocean_velocity.res.nc'
  if (file_exist(trim(filename))) then 

      if(tendency==THREE_LEVEL) then

          write (stdout(),'(a)') '  Reading THREE_LEVEL restart for velocity from INPUT/ocean_velocity.res.nc'
          write (stdout(),'(a)') '  Expecting two time records for each restart field.'

          call read_data(filename,'u',Velocity%u(:,:,:,1,tau),&
               Dom%domain2d,timelevel=1,append_pelist_name=ens_ocean)
          call read_data(filename,'u',Velocity%u(:,:,:,1,taup1),&
               Dom%domain2d,timelevel=2,append_pelist_name=ens_ocean)
          call read_data(filename,'v',Velocity%u(:,:,:,2,tau),&
               Dom%domain2d,timelevel=1,append_pelist_name=ens_ocean)
          call read_data(filename,'v',Velocity%u(:,:,:,2,taup1),&
               Dom%domain2d,timelevel=2,append_pelist_name=ens_ocean)

          call mpp_update_domains(Velocity%u(:,:,:,1,tau), Velocity%u(:,:,:,2,tau),&
               Dom%domain2d,gridtype=BGRID_NE)
          call mpp_update_domains(Velocity%u(:,:,:,1,taup1), Velocity%u(:,:,:,2,taup1),&
               Dom%domain2d,gridtype=BGRID_NE)

      elseif (tendency==TWO_LEVEL) then

          write (stdout(),'(/a)') '  Reading TWO_LEVEL restart for velocity from INPUT/ocean_velocity.res.nc'
          write (stdout(),'(a)')  '  Expecting only one time record for each restart field.'

          call field_size(filename,'u', siz)
          if (siz(4) > 1) then
            write(stdout(),'(/a)') '==>ERROR: Attempt to read ocean_velocity.res.nc from a 3-level time scheme (2 time records)'
            write(stdout(),'(a)')  '          when running mom4 with 2-level timestepping (only need 1 time record in restart).'
            write(stdout(),'(a)')  '          Reduce restart file to a single time record in order to avoid confusion.'
            call mpp_error(FATAL, &
                           'Reading 3-time lev ocean_velocity.res.nc (w/ 2 time records) while using 2-lev (needs only 1 record)')
          endif

          call read_data(filename,'u',Velocity%u(:,:,:,1,taup1),&
               Dom%domain2d,timelevel=1,append_pelist_name=ens_ocean)
          call read_data(filename,'v',Velocity%u(:,:,:,2,taup1),&
               Dom%domain2d,timelevel=1,append_pelist_name=ens_ocean)

          call mpp_update_domains(Velocity%u(:,:,:,1,taup1), Velocity%u(:,:,:,2,taup1),&
               Dom%domain2d,gridtype=BGRID_NE)


          write (stdout(),'(a)') '  Finished reading restart for velocity field'

          if(.not. Time%init) then  

              filename = 'INPUT/ocean_velocity_advection.res.nc'
              write (stdout(),'(/a)') '  Reading TWO_LEVEL restart from INPUT/ocean_velocity_advection.res.nc'
              write (stdout(),'(a)')  '  Expecting two time records for the advection velocity fields.'

              call read_data(filename,'advectionu',Velocity%advection(:,:,:,1,tau_m1),&
                   Dom%domain2d,timelevel=1,append_pelist_name=ens_ocean)
              call read_data(filename,'advectionu',Velocity%advection(:,:,:,1,tau_m0),&
                   Dom%domain2d,timelevel=2,append_pelist_name=ens_ocean)
              call read_data(filename,'advectionv',Velocity%advection(:,:,:,2,tau_m1),&
                   Dom%domain2d,timelevel=1,append_pelist_name=ens_ocean)
              call read_data(filename,'advectionv',Velocity%advection(:,:,:,2,tau_m0),&
                   Dom%domain2d,timelevel=2,append_pelist_name=ens_ocean)

              write (stdout(),'(a)') '  Finished reading restart for velocity advection operator'

              call mpp_update_domains(Velocity%advection(:,:,:,1,tau_m1),Dom%domain2d)
              call mpp_update_domains(Velocity%advection(:,:,:,1,tau_m0),Dom%domain2d)
              call mpp_update_domains(Velocity%advection(:,:,:,2,tau_m1),Dom%domain2d)
              call mpp_update_domains(Velocity%advection(:,:,:,2,tau_m0),Dom%domain2d)

          endif

      endif

      if(have_obc) then
          call ocean_obc_update_boundary(Velocity%u(:,:,:,1,taup1), 'M','s')
          call ocean_obc_update_boundary(Velocity%u(:,:,:,2,taup1), 'M','i')
          call ocean_obc_update_boundary(Velocity%u(:,:,:,1,taup1), 'Z','i')
          call ocean_obc_update_boundary(Velocity%u(:,:,:,2,taup1), 'Z','s')
      endif

 else

     if (Time%init) then
          call mpp_error(NOTE,&
         '==> From ocean_velocity_mod: Initializing velocity to zero since Time%init=.true. &
          &yet did not find INPUT/ocean_velocity.res.nc.')
     endif 

     if (.NOT. Time%init) then
          call mpp_error(FATAL,&
         '==> Error from ocean_velocity_mod: Expecting INPUT/ocean_velocity.res.nc to exist.&
          &This file was not found and Time%init=.false.')
     endif 

 endif
 
  if(debug_velocity) then
        write(stdout(),*) 'Initial Velocity chksum (tau) ==>'
        call ocean_velocity_chksum(Time, Velocity, tau, write_advection=.false.)
        write(stdout(),*) 'Initial Velocity chksum (taup1) ==>'
        call ocean_velocity_chksum(Time, Velocity, taup1, write_advection=.false.)
  endif

  if(zero_tendency) then
     call mpp_error(NOTE,'==>zero_tendency=.true. so will not time step the velocity field ')
  endif 
  if(zero_tendency_explicit) then
     call mpp_error(NOTE,'==>zero_tendency_explicit=.true. so will not include explicit part of the velocity tendency ')
  endif 
  if(zero_tendency_implicit) then
     call mpp_error(NOTE,'==>zero_tendency_implicit=.true. so will not include implicit part of the velocity tendency ')
  endif 

  call check_gravity_wave_cfl()
  
end subroutine ocean_velocity_init
! </SUBROUTINE> NAME="ocean_velocity_init"



!#######################################################################
! <SUBROUTINE NAME="check_gravity_wave_cfl">
!
! <DESCRIPTION>
! Check CFL for internal gravity waves. 
! </DESCRIPTION>
!
 subroutine check_gravity_wave_cfl()

    logical :: cfl_error=.false.
    real    :: gridsp, dtcg, max_dt_for_cgint0
    real    :: cfl_grid_factor, dtuv
    integer :: i, j, icg, jcg
       
   ! estimate the maximum timestep allowable for resolving internal gravity waves
    if(tendency==TWO_LEVEL)   then 
       cfl_grid_factor=1.0
       dtuv = dtime 
    elseif(tendency==THREE_LEVEL) then
       cfl_grid_factor=0.5
       dtuv = 0.5*dtime
    endif 
    icg = isc; jcg = jsc; gridsp = 1.0e20; max_dt_for_cgint = 1.0e6
    do j=jsc,jec
       do i=isc,iec
          if (Grd%kmu(i,j) > 0) then
              gridsp = 1.0/(Grd%dxu(i,j)*Grd%dxu(i,j)) + 1.0/(Grd%dyu(i,j)*Grd%dyu(i,j))
              gridsp = sqrt(1.0/gridsp) 
              dtcg = cfl_grid_factor*gridsp/(epsln+max_cgint)                       
              if (dtcg < max_dt_for_cgint) then
                  max_dt_for_cgint = dtcg; icg  = i; jcg  = j
              endif
          endif
       enddo
    enddo
    max_dt_for_cgint = nint(max_dt_for_cgint)
    max_dt_for_cgint = max_dt_for_cgint + 0.001*mpp_pe() ! to separate redundancies
    max_dt_for_cgint0 = max_dt_for_cgint
    call mpp_min (max_dt_for_cgint)

    ! show the most unstable location for baroclinic gravity waves
    if (max_dt_for_cgint == max_dt_for_cgint0 .and. nint(dtuv) > 0) then
        if (dtuv <= max_dt_for_cgint) then
            write (unit,'(/a,i4,a,i4,a,f6.2,a,f6.2,a)')' Baroclinic time step stability most nearly violated at U-cell (i,j) = (',&
                 icg+Dom%ioff,',',jcg+Dom%joff,'), (lon,lat) = (',Grd%xu(icg,jcg),',',Grd%yu(icg,jcg),').'
            write(unit,'(a,i6)')    '         The number of kmu-levels  at this point is ',Grd%kmu(icg,jcg) 
            write(unit,'(a,e12.6)') '         The dxu grid distance (m) at this point is ',Grd%dxu(icg,jcg) 
            write(unit,'(a,e12.6)') '         The dyu grid distance (m) at this point is ',Grd%dyu(icg,jcg) 
            cfl_error=.false.
        else
            write (unit,'(/a,i4,a,i4,a,f6.2,a,f6.2,a)')'=>Error: Baroclinic time step stability violated at U-cell (i,j) = (',&
                 icg+Dom%ioff,',',jcg+Dom%joff,'), (lon,lat) = (',Grd%xu(icg,jcg),',',Grd%yu(icg,jcg),').'
            write(unit,'(a,i6)')    '         The number of kmu-levels  at this point is ',Grd%kmu(icg,jcg) 
            write(unit,'(a,e12.6)') '         The dxu grid distance (m) at this point is ',Grd%dxu(icg,jcg) 
            write(unit,'(a,e12.6)') '         The dyu grid distance (m) at this point is ',Grd%dyu(icg,jcg) 
            cfl_error=.true.
        endif
        write (unit,'(a,f5.2,a/a,f6.0,a,f6.0,a)')&
             '         Due to a specified maximum baroclinic gravity wave speed of ',max_cgint,' m/s.',&
             '         "dtuv" must be less than ',max_dt_for_cgint,' sec. "dtuv" = ', dtuv,' sec.'
        if(cfl_error) then  
          call mpp_error(FATAL,'==>Error: time step instability detected for baroclinic gravity waves in ocean_model_mod')
        endif 
        write (stdlog(),'(/a/)') ' Note: A more appropriate maximum &
             &baroclinic gravity wave speed can be specified via namelist.'
    endif
    max_dt_for_cgint = nint(max_dt_for_cgint)
    
  end subroutine check_gravity_wave_cfl
! </SUBROUTINE> NAME="check_gravity_wave_cfl"


 
!#######################################################################
! <SUBROUTINE NAME="ocean_explicit_accel_a">
!
! <DESCRIPTION>
! Time explicit contributions to thickness weighted acceleration. 
! Omit here the Coriolis force and verrtical friction here. 
! </DESCRIPTION>
!
subroutine ocean_explicit_accel_a(Velocity, Time, Adv_vel, Thickness, rho, pme, river, upme, uriver, accel)

  type(ocean_velocity_type), intent(inout)             :: Velocity
  type(ocean_time_type), intent(in)                    :: Time
  type(ocean_adv_vel_type), intent(in)                 :: Adv_vel
  type(ocean_thickness_type), intent(in)               :: Thickness
  real, intent(in), dimension(isd:ied,jsd:jed,nk)      :: rho
  real, intent(in), dimension(isd:ied,jsd:jed)         :: pme
  real, intent(in), dimension(isd:ied,jsd:jed)         :: river
  real, intent(in), dimension(isd:ied,jsd:jed,2)       :: upme
  real, intent(in), dimension(isd:ied,jsd:jed,2)       :: uriver 
  real, intent(inout), dimension(isd:ied,jsd:jed,nk,2) :: accel

  integer :: i, j, k, n
  integer :: taum1, tau, taup1
  integer :: tau_m0, tau_m1, tau_m2

  taum1 = Time%taum1
  tau   = Time%tau
  taup1 = Time%taup1

  tau_m2 = Time%tau_m2
  tau_m1 = Time%tau_m1
  tau_m0 = Time%tau_m0

  accel = 0.0

  if (.not. zero_tendency .and. .not. zero_tendency_explicit) then 

      if(tendency==THREE_LEVEL) then 

          accel(isc:iec,jsc:jec,:,:)  =  &
               -pressure_gradient(Time, Velocity, Thickness, rho(:,:,:))                          &
               -horz_advection_of_velocity(Time, Thickness, Velocity, Adv_vel)                    &
               -vertical_advection_of_velocity(Time, Velocity, Adv_vel, pme, river, upme, uriver) &
               +horz_lap_friction(Time, Thickness, Velocity)                                      &
               +horz_bih_friction(Time, Thickness, Velocity)                               

      elseif(tendency==TWO_LEVEL) then   

          Velocity%advection(isc:iec,jsc:jec,:,:,tau_m0) =                      &
               horz_advection_of_velocity(Time, Thickness, Velocity, Adv_vel)   &
               +vertical_advection_of_velocity(Time, Velocity, Adv_vel, pme, river, upme, uriver)

          do n=1,2
             do k=1,nk
                do j=jsc,jec
                   do i=isc,iec
                      accel(i,j,k,n) = -abtau_m0*Velocity%advection(i,j,k,n,tau_m0) & 
                                       -abtau_m1*Velocity%advection(i,j,k,n,tau_m1) & 
                                       -abtau_m2*Velocity%advection(i,j,k,n,tau_m2)  
                   enddo
                enddo
             enddo
          enddo

          accel(isc:iec,jsc:jec,:,:)  =   accel(isc:iec,jsc:jec,:,:)  &
                                        + horz_lap_friction(Time, Thickness, Velocity) &
                                        + horz_bih_friction(Time, Thickness, Velocity) &
                                        - pressure_gradient(Time, Velocity, Thickness, rho(:,:,:)) 

      endif

  endif

  if (debug_velocity) then
      write(stdout(),*) 'explicit accel_a(1) ==> ',mpp_chksum(accel(isc:iec,jsc:jec,:,1))
      write(stdout(),*) 'explicit accel_a(2) ==> ',mpp_chksum(accel(isc:iec,jsc:jec,:,2))
  endif

end subroutine ocean_explicit_accel_a
! </SUBROUTINE> NAME="ocean_explicit_accel_a"

 
!#######################################################################
! <SUBROUTINE NAME="ocean_explicit_accel_b">
!
! <DESCRIPTION>
! Add Coriolis force and explicit vertical friction to explicit-time
! thickness weighted acceleration. 
! </DESCRIPTION>
!
subroutine ocean_explicit_accel_b(visc_cbu, Time, Velocity, Thickness, accel)

  real, intent(in), dimension(isd:ied,jsd:jed,nk)      :: visc_cbu
  type(ocean_time_type), intent(in)                    :: Time
  type(ocean_velocity_type), intent(in)                :: Velocity
  type(ocean_thickness_type), intent(in)               :: Thickness
  real, intent(inout), dimension(isd:ied,jsd:jed,nk,2) :: accel

  if (.not. zero_tendency .and. .not. zero_tendency_implicit) then 
      accel(isc:iec,jsc:jec,:,:)  =  accel(isc:iec,jsc:jec,:,:)              &
                                   + coriolis(Time, Velocity, Thickness)     &
                                   + vert_frict(Time, Thickness, aidif, Velocity, visc_cbu)
  endif

  if (debug_velocity) then
      write(stdout(),*) 'explicit accel_b(1) ==> ',mpp_chksum(accel(isc:iec,jsc:jec,:,1))
      write(stdout(),*) 'explicit accel_b(2) ==> ',mpp_chksum(accel(isc:iec,jsc:jec,:,2))
  endif

end subroutine ocean_explicit_accel_b
! </SUBROUTINE> NAME="ocean_explicit_accel_b"



!#######################################################################
! <SUBROUTINE NAME="update_ocean_velocity">
!
! <DESCRIPTION>
! Update velocity components
! </DESCRIPTION>
!
subroutine update_ocean_velocity(Time, Thickness, accel, barotropic_split, Ext_mode, Velocity)

  type(ocean_time_type), intent(in)                 :: Time
  type(ocean_thickness_type), intent(in)            :: Thickness
  real, dimension(isd:ied,jsd:jed,nk,2), intent(in) :: accel
  integer, intent(in)                               :: barotropic_split
  type(ocean_external_mode_type), intent(inout)     :: Ext_mode
  type(ocean_velocity_type), intent(inout)          :: Velocity

  real, dimension(isd:ied,jsd:jed) :: tmp
  logical :: velocity_truncated 
  integer :: i, j, k, n, taum1, tau, taup1

  taum1 = Time%taum1
  tau   = Time%tau
  taup1 = Time%taup1

  tmp   = 0.0

  if(zero_tendency) then 

      Velocity%u(isc:iec,jsc:jec,:,:,taup1) = Velocity%u(isc:iec,jsc:jec,:,:,tau)

  else 

      if (barotropic_split > 1) then  

          do n=1,2

             ! compute ustar 
             ! smg: keep /Thickness%dhu instead of *Thickness%dhur to bitwise reproduce old results
             do k=1,nk
                do j=jsc,jec
                   do i=isc,iec
                      Velocity%u(i,j,k,n,taup1) = &
                       (Thickness%dhu(i,j,k,taum1)*Velocity%u(i,j,k,n,taum1) + dtime*accel(i,j,k,n)*Grd%umask(i,j,k))  &
                                                   /Thickness%dhu(i,j,k,taup1)
                   enddo
                enddo
             enddo

             if(truncate_velocity) then
                 velocity_truncated = .false. 
                 wrk1 = 0.0
                 do k=1,nk
                    do j=jsc,jec
                       do i=isc,iec
                          if (abs(Velocity%u(i,j,k,n,taup1)) > vel_max) then
                              Velocity%u(i,j,k,n,taup1) = sign(vel_max,Velocity%u(i,j,k,n,taup1))
                              wrk1(i,j,k)               = 1.0
                              velocity_truncated        = .true.  
                          endif
                       enddo
                    enddo
                 enddo
                 if(truncate_verbose .and. velocity_truncated) then 
                    do k=1,nk
                       do j=jsc,jec
                          do i=isc,iec
                             if (wrk1(i,j,k) ==1.0) then 
                               write(stdout(),*) 'WARNING: truncated baroclinic velocity component(',n,') at (i,j,k) = '&
                                ,i+Dom%ioff,j+Dom%joff,k &
                                ,'which is at (x,y,z) = ',Grd%dxu(i,j), ',', Grd%dyu(i,j), ',', Grd%kmu(i,j)
                             endif
                          enddo
                       enddo
                    enddo
                 endif
             endif

             ! sum ustar*dhu
             tmp=0.0
             do k=1,nk
                do j=jsc,jec
                   do i=isc,iec  
                      tmp(i,j) = tmp(i,j) + Velocity%u(i,j,k,n,taup1)*Thickness%dhu(i,j,k,taup1)*Grd%umask(i,j,k)
                   enddo
                enddo
             enddo

             ! (vertical mean) - (value from barotropic integration) 
             do j=jsc,jec
                do i=isc,iec  
                   tmp(i,j) = Grd%umask(i,j,1)*(tmp(i,j)-Ext_mode%ud(i,j,n,taup1))/(Grd%hu(i,j)+Ext_mode%eta_u(i,j,taup1)+epsln) 
                enddo
             enddo

             ! replace vertical mean with value from barotropic integration 
             do k=1,nk 
                do j=jsc,jec
                   do i=isc,iec  
                      Velocity%u(i,j,k,n,taup1) = (Velocity%u(i,j,k,n,taup1)-tmp(i,j))*Grd%umask(i,j,k) 
                   enddo
                enddo
             enddo

          enddo  ! finish the n=1,2 do-loop 


      elseif(barotropic_split==1) then   ! no splitting case: estimate u(tau+1) using surface pressure gradient

          do n=1,2

             do k=1,nk 
                do j=jsc,jec
                   do i=isc,iec 
                      Velocity%u(i,j,k,n,taup1) = (Thickness%dhu(i,j,k,taum1)*Velocity%u(i,j,k,n,taum1) + &
                           dtime*(accel(i,j,k,n) - Ext_mode%grad_ps(i,j,n)*Grd%umask(i,j,k)))*Thickness%dhur(i,j,k)
                   enddo
                enddo
             enddo

             Ext_mode%ud(:,:,:,taup1) = 0.0
             do k=1,nk
                do j=jsc,jec
                   do i=isc,iec
                      Ext_mode%ud(i,j,n,taup1) = Ext_mode%ud(i,j,n,taup1) + Velocity%u(i,j,k,n,taup1)*Thickness%dhu(i,j,k,taup1)
                   enddo
                enddo
             enddo

          enddo

          call mpp_update_domains (Ext_mode%ud(:,:,1,taup1),Ext_mode%ud(:,:,2,taup1), Dom%domain2d, gridtype=BGRID_NE)

          ! construct convergence of ud at (tau+1)
          Ext_mode%convud_t(:,:,taup1)  = - DIV_UD(Ext_mode%ud(:,:,:,taup1))
          call mpp_update_domains (Ext_mode%convud_t(:,:,taup1), Dom%domain2d)

      endif   ! endif for splitting 

  endif   ! endif for zero_tendency

  !--- update velocity at the global halo points to make the gradient is 0 accross boundary
  if(have_obc) then
!     call mpp_update_domains (Velocity%u(:,:,:,:,taup1), Dom%domain2d)
     call ocean_obc_update_boundary(Velocity%u(:,:,:,1,taup1), 'M','s')
     call ocean_obc_update_boundary(Velocity%u(:,:,:,2,taup1), 'M','i')
     call ocean_obc_update_boundary(Velocity%u(:,:,:,1,taup1), 'Z','i')
     call ocean_obc_update_boundary(Velocity%u(:,:,:,2,taup1), 'Z','s')
  endif
  if (debug_velocity) then        
     write(stdout(),*) 'From update_ocean_velocity: velocity(taup1) chksum  ===>'
     call ocean_velocity_chksum(Time, Velocity, taup1, write_advection=.false.)
  endif

  if (id_u > 0) used = send_data (id_u, Velocity%u(:,:,:,1,taup1), &
                       Time%model_time, rmask=Grd%umask(:,:,:), &
                       is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  if (id_v > 0) used = send_data (id_v, Velocity%u(:,:,:,2,taup1), &
                       Time%model_time, rmask=Grd%umask(:,:,:), &
                       is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  if (id_usurf > 0) used = send_data (id_usurf, Velocity%u(:,:,1,1,taup1), &
                           Time%model_time, rmask=Grd%umask(:,:,1), &
                           is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  if (id_vsurf > 0) used = send_data (id_vsurf, Velocity%u(:,:,1,2,taup1), &
                           Time%model_time, rmask=Grd%umask(:,:,1), &
                           is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

end subroutine update_ocean_velocity
! </SUBROUTINE> NAME="update_ocean_velocity"


!#######################################################################
! <SUBROUTINE NAME="ocean_implicit_friction">
!
! <DESCRIPTION>
! Contributions to thickness weighted acceleration from implicit
! vertical friction. 
! </DESCRIPTION>
!
subroutine ocean_implicit_friction(Time, Thickness, visc_cbu, Velocity, accel)

  type(ocean_time_type), intent(in)                    :: Time  
  type(ocean_thickness_type), intent(in)               :: Thickness
  real, intent(in), dimension(isd:ied,jsd:jed,nk)      :: visc_cbu
  type(ocean_velocity_type), intent(in)                :: Velocity
  real, intent(inout), dimension(isd:ied,jsd:jed,nk,2) :: accel

  integer :: taum1, tau, taup1
  integer :: i, j, k, n
  integer :: now

  taum1 = Time%taum1
  tau   = Time%tau
  taup1 = Time%taup1

  vfrict_impl = 0.0
  wrk1_v      = 0.0
  wrk1        = 0.0

  ! acceleration due to time-implicit vertical friction 
  if (aidif /= 0.0) then  

      ! construct updated velocity due to time-explicit contributions, sans grad_ps  
      ! smg: keep /Thickness%dhu instead of *Thickness%dhur to bitwise reproduce old results
      do n=1,2
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  wrk1_v(i,j,k,n) = (Thickness%dhu(i,j,k,taum1)*Velocity%u(i,j,k,n,taum1) + dtime*accel(i,j,k,n)) &
                                    /Thickness%dhu(i,j,k,taup1)
               enddo
            enddo
         enddo
      enddo

      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               vfrict_impl(i,j,k,:) = wrk1_v(i,j,k,:)
               wrk1(i,j,k)          = Thickness%dhu(i,j,k,taup1)
            enddo
         enddo
      enddo

      ! call invtri to invert tridiagonal and update taup1 velocity. include bmf and smf       
      call invtri (wrk1_v(:,:,:,1), Velocity%smf(:,:,1), Velocity%bmf(:,:,1), visc_cbu, dtime, Grd%kmu, &
                   Grd%umask, wrk1, Thickness%dhwu, aidif, nk)
      call invtri (wrk1_v(:,:,:,2), Velocity%smf(:,:,2), Velocity%bmf(:,:,2), visc_cbu, dtime, Grd%kmu, &
                   Grd%umask, wrk1, Thickness%dhwu, aidif, nk)

      ! compute time-implicit tendency for diagnostics 
      ! update thickness weighted acceleration due to implicit vertical friction
      do n=1,2
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  vfrict_impl(i,j,k,n) = dtimer*(wrk1_v(i,j,k,n)-vfrict_impl(i,j,k,n))
                  accel(i,j,k,n)       = accel(i,j,k,n) + Thickness%dhu(i,j,k,taup1)*vfrict_impl(i,j,k,n)
               enddo
            enddo
         enddo
      enddo

  endif  !endif for aidif /= 0

  if (id_vfrict_impl(1) > 0) used = send_data(id_vfrict_impl(1), vfrict_impl(:,:,:,1), &
                                    Time%model_time, rmask=Grd%umask(:,:,:), &
                                    is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  if (id_vfrict_impl(2) > 0) used = send_data(id_vfrict_impl(2), vfrict_impl(:,:,:,2), &
                                    Time%model_time, rmask=Grd%umask(:,:,:), &
                                    is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)  

  if (debug_velocity) then
    write(stdout(),*) 'explicit+implicit ==> ',mpp_chksum(accel(isc:iec,jsc:jec,:,:))
    write(stdout(),*) 'velocity after calling invtri ==> ',mpp_chksum(wrk1_v(isc:iec,jsc:jec,:,:))
  endif

end subroutine ocean_implicit_friction
! </SUBROUTINE> NAME="ocean_implicit_friction"


!#######################################################################
! <SUBROUTINE NAME="ocean_implicit_coriolis">
!
! <DESCRIPTION>
! Contributions to acceleration from time-implicit Coriolis force.
! </DESCRIPTION>
!
subroutine ocean_implicit_coriolis(Time, Velocity, acor, accel)

  type(ocean_time_type), intent(in)                    :: Time
  type(ocean_velocity_type), intent(in)                :: Velocity
  real, intent(in)                                     :: acor
  real, intent(inout), dimension(isd:ied,jsd:jed,nk,2) :: accel

  integer :: i,j,k
  real    :: dtimeacor
  real    :: lambda,factor
  
  if(acor > 0.0) then 

      dtimeacor = dtime*acor 

      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk1_v(i,j,k,:) = accel(i,j,k,:)  
               lambda          = dtimeacor*Grd%f(i,j)
               factor          = 1.0/(1.0 + lambda*lambda) 
               wrk1(i,j,k)     = (accel(i,j,k,1) + lambda*accel(i,j,k,2))*factor
               wrk2(i,j,k)     = (accel(i,j,k,2) - lambda*accel(i,j,k,1))*factor
            enddo
         enddo
      enddo

      do k=1,nk
         do j=jsc,jec
            do i=isc,iec       
               accel(i,j,k,1)  = wrk1(i,j,k) 
               accel(i,j,k,2)  = wrk2(i,j,k) 
               wrk1_v(i,j,k,:) = accel(i,j,k,:) - wrk1_v(i,j,k,:) 
            enddo
         enddo
      enddo

  endif

  ! send to diagnostics manager  
  if (id_accel(1) > 0) used = send_data(id_accel(1), accel(:,:,:,1), &
                              Time%model_time, rmask=Grd%umask(:,:,:), &
                              is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

  if (id_accel(2) > 0) used = send_data(id_accel(2), accel(:,:,:,2), &
                              Time%model_time, rmask=Grd%umask(:,:,:), &
                              is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

  if (id_ucori_impl > 0) used = send_data(id_ucori_impl, wrk1_v(:,:,:,1), &
                                Time%model_time, rmask=Grd%umask(:,:,:), &
                                is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)  

  if (id_vcori_impl > 0) used = send_data(id_vcori_impl, wrk1_v(:,:,:,2), &
                                Time%model_time, rmask=Grd%umask(:,:,:), &
                                is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)  

  if (debug_velocity) then
      write(stdout(),*) 'explicit+implicit ==> ',mpp_chksum(accel(isc:iec,jsc:jec,:,:))
      write(stdout(),*) 'velocity after calling implicit-coriolis ==> ',mpp_chksum(wrk1_v(isc:iec,jsc:jec,:,:))
  endif

end subroutine ocean_implicit_coriolis
! </SUBROUTINE> NAME="ocean_implicit_coriolis"


!#######################################################################
! <SUBROUTINE NAME="ocean_velocity_rstrt">
!
! <DESCRIPTION>
! Write the velocity field to a restart
! </DESCRIPTION>
subroutine ocean_velocity_rstrt(Time, Velocity, ens_ocean)

  type(ocean_time_type), intent(in)     :: Time
  type(ocean_velocity_type), intent(in) :: Velocity
  logical, intent(in), optional         :: ens_ocean

  integer :: n
  integer :: tau, taup1
  integer :: tau_m0, tau_m1
  integer :: yr, mon, day, hr, min, sec
  character(len=128) :: filename
  character(len=10) :: rdte

  tau    = Time%tau
  taup1  = Time%taup1

  tau_m1 = Time%tau_m1
  tau_m0 = Time%tau_m0

  call get_date(Time%model_time,yr,mon,day,hr,min,sec)
  write(rdte,'(i4,3i2.2)') yr,mon,day,hr

  filename = 'IRESTART/'//rdte//'ocean_velocity.res'

  if(tendency==THREE_LEVEL) then

      call write_data(filename, 'u', Velocity%u(:,:,:,1,tau), &
           domain=Dom%domain2d,append_pelist_name=ens_ocean)
      call write_data(filename, 'u', Velocity%u(:,:,:,1,taup1), &
           domain=Dom%domain2d,append_pelist_name=ens_ocean)
      call write_data(filename, 'v', Velocity%u(:,:,:,2,tau), &
           domain=Dom%domain2d,append_pelist_name=ens_ocean)
      call write_data(filename, 'v', Velocity%u(:,:,:,2,taup1), &
           domain=Dom%domain2d,append_pelist_name=ens_ocean)

      write(stdout(),*) 'Ending Velocity chksum (tau) ==>'
      call ocean_velocity_chksum(Time, Velocity, tau, write_advection=.false.)
      write(stdout(),*) 'Ending Velocity chksum (taup1) ==>'
      call ocean_velocity_chksum(Time, Velocity, taup1, write_advection=.false.)

  elseif(tendency==TWO_LEVEL) then

      call write_data(filename, 'u', Velocity%u(:,:,:,1,taup1), &
           domain=Dom%domain2d,append_pelist_name=ens_ocean)
      call write_data(filename, 'v', Velocity%u(:,:,:,2,taup1), &
           domain=Dom%domain2d,append_pelist_name=ens_ocean)

      filename = 'IRESTART/'//rdte//'ocean_velocity_advection.res.nc'
      call write_data(filename,'advectionu',Velocity%advection(:,:,:,1,tau_m1),&
           Dom%domain2d,append_pelist_name=ens_ocean)
      call write_data(filename,'advectionu',Velocity%advection(:,:,:,1,tau_m0),&
           Dom%domain2d,append_pelist_name=ens_ocean)
      call write_data(filename,'advectionv',Velocity%advection(:,:,:,2,tau_m1),&
           Dom%domain2d,append_pelist_name=ens_ocean)
      call write_data(filename,'advectionv',Velocity%advection(:,:,:,2,tau_m0),&
           Dom%domain2d,append_pelist_name=ens_ocean)

      write(stdout(),*) 'Ending Velocity chksum (taup1) ==>'
      call ocean_velocity_chksum(Time, Velocity, taup1, write_advection=.true.)

  endif

end subroutine ocean_velocity_rstrt
! </SUBROUTINE> NAME="ocean_velocity_rstrt"



!#######################################################################
! <SUBROUTINE NAME="ocean_velocity_end">
!
! <DESCRIPTION>
! Write the velocity field to a restart 
! </DESCRIPTION>
subroutine ocean_velocity_end(Time, Velocity, ens_ocean)

  type(ocean_time_type), intent(in)     :: Time
  type(ocean_velocity_type), intent(in) :: Velocity
  logical, intent(in), optional         :: ens_ocean
  
  integer :: n
  integer :: tau, taup1
  integer :: tau_m0, tau_m1
  character(len=128) :: filename
  
  tau    = Time%tau
  taup1  = Time%taup1

  tau_m1 = Time%tau_m1
  tau_m0 = Time%tau_m0

  filename = 'RESTART/ocean_velocity.res'

  if(tendency==THREE_LEVEL) then 

      call write_data(filename, 'u', Velocity%u(:,:,:,1,tau), &
           domain=Dom%domain2d,append_pelist_name=ens_ocean)
      call write_data(filename, 'u', Velocity%u(:,:,:,1,taup1), &
           domain=Dom%domain2d,append_pelist_name=ens_ocean)
      call write_data(filename, 'v', Velocity%u(:,:,:,2,tau), &
           domain=Dom%domain2d,append_pelist_name=ens_ocean)
      call write_data(filename, 'v', Velocity%u(:,:,:,2,taup1), &
           domain=Dom%domain2d,append_pelist_name=ens_ocean)

      write(stdout(),*) 'Ending Velocity chksum (tau) ==>'
      call ocean_velocity_chksum(Time, Velocity, tau, write_advection=.false.)
      write(stdout(),*) 'Ending Velocity chksum (taup1) ==>'
      call ocean_velocity_chksum(Time, Velocity, taup1, write_advection=.false.)

  elseif(tendency==TWO_LEVEL) then 

      call write_data(filename, 'u', Velocity%u(:,:,:,1,taup1), &
           domain=Dom%domain2d,append_pelist_name=ens_ocean)
      call write_data(filename, 'v', Velocity%u(:,:,:,2,taup1), &
           domain=Dom%domain2d,append_pelist_name=ens_ocean)

      filename = 'RESTART/ocean_velocity_advection.res.nc'
      call write_data(filename,'advectionu',Velocity%advection(:,:,:,1,tau_m1),&
           Dom%domain2d,append_pelist_name=ens_ocean)
      call write_data(filename,'advectionu',Velocity%advection(:,:,:,1,tau_m0),&
           Dom%domain2d,append_pelist_name=ens_ocean)
      call write_data(filename,'advectionv',Velocity%advection(:,:,:,2,tau_m1),&
           Dom%domain2d,append_pelist_name=ens_ocean)
      call write_data(filename,'advectionv',Velocity%advection(:,:,:,2,tau_m0),&
           Dom%domain2d,append_pelist_name=ens_ocean)

      write(stdout(),*) 'Ending Velocity chksum (taup1) ==>'
      call ocean_velocity_chksum(Time, Velocity, taup1, write_advection=.true.)

  endif

  
  return

end subroutine ocean_velocity_end
! </SUBROUTINE> NAME="ocean_velocity_end"


!#######################################################################
! <SUBROUTINE NAME="energy_analysis">
!
! <DESCRIPTION>
! Perform energy analysis by taking scalar product of horizontal
! velocity with the velocity equations and integrating over the ocean volume.
! </DESCRIPTION>
!
subroutine energy_analysis (Time, Thickness, rho, pme, river, upme, uriver, &
                            visc_cbu, Dens, Ext_mode, Adv_vel, Velocity, acor)

  type(ocean_time_type), intent(in)          :: Time
  type(ocean_thickness_type), intent(in)     :: Thickness
  type(ocean_density_type), intent(in)       :: Dens
  type(ocean_adv_vel_type), intent(in)       :: Adv_vel
  type(ocean_external_mode_type), intent(in) :: Ext_mode
  type(ocean_velocity_type), intent(in)      :: Velocity

  real, intent(in) :: rho(isd:ied,jsd:jed,nk)
  real, intent(in) :: pme(isd:ied,jsd:jed)
  real, intent(in) :: river(isd:ied,jsd:jed)
  real, intent(in) :: upme(isd:ied,jsd:jed,2)
  real, intent(in) :: uriver(isd:ied,jsd:jed,2)
  real, intent(in) :: visc_cbu(isd:ied,jsd:jed,nk)
  real, intent(in) :: acor

  integer :: i, j, k, n, kz
  real    :: boxvol, boxarea, ucell_volume_k(nk), ucell_volume
  real    :: uint, uext, term, vel, fx
  real    :: ke_tot, pe_tot
  
  real, dimension(isd:ied,jsd:jed,2,3) :: ubar
  real, dimension(10) :: engint            ! volume averaged internal mode energy integral components
  real, dimension(10) :: engext            ! volume averaged external mode energy integral components
  real                :: buoy              ! volume averaged buoyancy
  real                :: kinetic, plicin, plicex, buoerr, enleak
  integer             :: tau, taum1, taup1

  real, dimension(isd:ied,jsd:jed) :: pme_u
  real, dimension(isd:ied,jsd:jed) :: river_u
  real, dimension(isd:ied,jsd:jed) :: convud_u
  
  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, '==>Error in ocean_velocity_mod (energy_analysis): module must be initialized')
  endif 

  tau   = Time%tau
  taum1 = Time%taum1
  taup1 = Time%taup1

  ! obtain volume of the ocean (time dependent due to eta_t)
  ucell_volume = Grd%ucellv
  
  buoy      = 0.0
  engint(:) = 0.0
  engext(:) = 0.0
  ubar(:,:,1,tau) = Grd%umask(:,:,1)*Ext_mode%ud(:,:,1,tau)/(Grd%hu(:,:)+Ext_mode%eta_u(:,:,tau)+epsln)
  ubar(:,:,2,tau) = Grd%umask(:,:,1)*Ext_mode%ud(:,:,2,tau)/(Grd%hu(:,:)+Ext_mode%eta_u(:,:,tau)+epsln)

  ! u dot the time rate of change 
  do n=1,2
    do k=1,nk
      do j=jsc,jec
        do i=isc,iec
          boxarea   = Grd%umask(i,j,k)*Grd%dau(i,j)
          uext      = ubar(i,j,n,tau)
          uint      = Velocity%u(i,j,k,n,tau) - uext
          term      = (Thickness%dhu(i,j,k,taup1)*Velocity%u(i,j,k,n,taup1) &
                      -Thickness%dhu(i,j,k,taum1)*Velocity%u(i,j,k,n,taum1)) &
                      *dtimer*boxarea
          engint(1) = engint(1) + uint*term
          engext(1) = engext(1) + uext*term
        enddo
      enddo
    enddo
  enddo

  ! work done by horizontal pressure gradients
  ! pressure_gradient returns thickness weighted pressure gradient (m^2/s^2)
  wrk1_v(isc:iec,jsc:jec,:,:) = pressure_gradient(Time, Velocity, Thickness, rho(:,:,:), .false.)
  do n=1,2
    do k=1,nk
      do j=jsc,jec
        do i=isc,iec
          boxarea   = Grd%dau(i,j)
          uext      = ubar(i,j,n,tau)
          uint      = Velocity%u(i,j,k,n,tau) - uext  
          term      = wrk1_v(i,j,k,n)*boxarea
          engint(6) = engint(6) - uint*term
          engext(6) = engext(6) - uext*term
        enddo
      enddo
    enddo
  enddo

  ! Add contribution from gradients in the surface pressure.  Note  
  ! that uint*grad_ps vertically integrates to zero, by definition.  
  ! So its contribution is not computed here. 
  do j=jsc,jec
    do i=isc,iec
      boxvol    = Grd%dau(i,j)*(Grd%hu(i,j)+Ext_mode%eta_u(i,j,tau))
      term      = Ext_mode%grad_ps(i,j,1)*ubar(i,j,1,tau) + Ext_mode%grad_ps(i,j,2)*ubar(i,j,2,tau)            
      engext(6) = engext(6) - term*boxvol
    enddo
  enddo

  ! work done by horizontal advection 
  ! horz_advection_of_velocity returns thickness weighted horizontal advection (m^2/s^2)
  wrk1_v(isc:iec,jsc:jec,:,:) = horz_advection_of_velocity(Time, Thickness, Velocity, Adv_vel, .false.)
  do n=1,2
    do k=1,nk
      do j=jsc,jec
        do i=isc,iec
          boxarea   = Grd%dau(i,j)
          uext      = ubar(i,j,n,tau)
          uint      = Velocity%u(i,j,k,n,tau) - uext 
          term      = wrk1_v(i,j,k,n)*boxarea  
          engint(2) = engint(2) - uint*term
          engext(2) = engext(2) - uext*term
        enddo
      enddo
    enddo
  enddo

  ! work done by vertical advection
  ! vertical_advection_of_velocity returns thickness weighted vertical advection (m^2/s^2) 
  wrk1_v(isc:iec,jsc:jec,:,:) = vertical_advection_of_velocity(Time, Velocity, Adv_vel, pme, river, upme, uriver, .false.) 
  do n=1,2
     do k=1,nk
        do j=jsc,jec
           do i=isc,iec
              boxarea   = Grd%dau(i,j)
              uext      = ubar(i,j,n,tau)
              uint      = Velocity%u(i,j,k,n,tau) - uext   
              term      = wrk1_v(i,j,k,n)*boxarea
              engint(3) = engint(3) - uint*term
              engext(3) = engext(3) - uext*term
           enddo
        enddo
     enddo
  enddo

  ! work done by fresh water forcing and convergence 
  pme_u(:,:)    = REMAP_BT_TO_BU(pme)
  river_u(:,:)  = REMAP_BT_TO_BU(river)
  convud_u(:,:) = Grd%umask(:,:,1)*REMAP_BT_TO_BU(Ext_mode%convud_t(:,:,tau)) 
  do n=1,2
    do j=jsc,jec
      do i=isc,iec
        boxarea    = Grd%dau(i,j)*Grd%umask(i,j,1)
        vel        = Velocity%u(i,j,1,n,tau)
        uext       = ubar(i,j,n,tau)
        uint       = vel - uext   
        term       =  (pme_u(i,j)*upme(i,j,n) + river_u(i,j)*uriver(i,j,n))  &
                    + (pme_u(i,j)             + river_u(i,j))*vel            &
                    + vel*convud_u(i,j)
        term       = 0.5*term*boxarea 
        engint(10) = engint(10) + uint*term
        engext(10) = engext(10) + uext*term
      enddo
    enddo
  enddo

  ! work done by horz friction
  ! horz_lap_friction and horz_bih_friction return thickness weighted horizontal friction (m^2/s^2)
  wrk1_v(isc:iec,jsc:jec,:,:) = horz_lap_friction(Time, Thickness, Velocity, diag_flag=.false.) &
                              + horz_bih_friction(Time, Thickness, Velocity, diag_flag=.false.)
  do n=1,2
    do k=1,nk
      do j=jsc,jec
        do i=isc,iec
          boxarea   = Grd%dau(i,j)
          uext      = ubar(i,j,n,tau)
          uint      = Velocity%u(i,j,k,n,tau) - uext
          term      = wrk1_v(i,j,k,n)*boxarea    
          engint(4) = engint(4) + uint*term
          engext(4) = engext(4) + uext*term
        enddo
      enddo
    enddo
  enddo

  ! work done by vertical friction (contributions from wind and bottom drag removed below)
  ! implicit acceleration has been saved from previous call to ocean_implicit_accel
  ! vert_frict returns thickness weighted vertical friction (m^2/s^2)
  wrk1_v(isc:iec,jsc:jec,:,:) =  vert_frict(Time, Thickness, aidif, Velocity,visc_cbu, .false.) &
                               + vfrict_impl(isc:iec,jsc:jec,:,:)
  do n=1,2
    do k=1,nk
      do j=jsc,jec
        do i=isc,iec
          boxarea   = Grd%dau(i,j)
          uext      = ubar(i,j,n,tau)
          uint      = Velocity%u(i,j,k,n,tau) - uext  
          term      = wrk1_v(i,j,k,n)*boxarea 
          engint(5) = engint(5) + uint*term
          engext(5) = engext(5) + uext*term
        enddo
      enddo
    enddo
  enddo

  ! work due to Coriolis force (zero on B_grid only when acor=0.0)
  ! coriolis returns thickness weighted Coriolis force (m^2/s^2) 
  wrk1_v(isc:iec,jsc:jec,:,:) = coriolis(Time, Velocity, Thickness, .false., .true.)
  do n=1,2
    do k=1,nk
      do j=jsc,jec
        do i=isc,iec
          boxarea   = Grd%dau(i,j)
          uext      = ubar(i,j,n,tau)
          uint      = Velocity%u(i,j,k,n,tau) - uext  
          term      = wrk1_v(i,j,k,n)*boxarea 
          engint(9) = engint(9) + uint*term   
          engext(9) = engext(9) + uext*term   
        enddo
      enddo
    enddo
  enddo

  ! work done by wind stress and bottom drag 
  do n=1,2

    k = 1

    ! wind stress 
    do j=jsc,jec
      do i=isc,iec
        uext = ubar(i,j,n,tau)
        uint = Velocity%u(i,j,k,n,tau) - uext
        term = Grd%umask(i,j,k)*Velocity%smf(i,j,n)*Grd%dau(i,j)  
        engint(7) = engint(7) + uint*term    
        engext(7) = engext(7) + uext*term
      enddo
    enddo

    ! bottom drag
    do j=jsc,jec
      do i=isc,iec
        k = Grd%kmu(i,j)
        if (k /= 0) then
          uext = ubar(i,j,n,tau)
          uint = Velocity%u(i,j,k,n,tau) - uext    
          term = Grd%umask(i,j,k)*Velocity%bmf(i,j,n)*Grd%dau(i,j)
          engint(8) = engint(8) - uint*term    
          engext(8) = engext(8) - uext*term
        endif
      enddo
    enddo

  enddo


  ! subtract wind and bottom drag contributions from vertical friction 
   engint(5) = engint(5) - engint(7) - engint(8) 
   engext(5) = engext(5) - engext(7) - engext(8) 


  ! work done by buoyancy
  buoy = 0.0
  do j=jsc,jec
    do i=isc,iec
      kz = Grd%kmt(i,j)
      if (kz /= 0) then
        term = -Grd%dat(i,j)*rho0r*Adv_vel%w_bt(i,j,0)*(grav*rho(i,j,1)*Grd%dzw(0)+Ext_mode%ps(i,j))
        buoy = buoy + term
        fx   = Grd%dat(i,j)*grav*rho0r*0.5    ! 0.5 arises from vertical average on density 
        do k=2,kz
          term =-fx*Adv_vel%w_bt(i,j,k-1)*Thickness%dhwt(i,j,k-1)*(rho(i,j,k-1) + rho(i,j,k))
          buoy = buoy + term
        enddo
      endif
    enddo
  enddo

  ! sigma-contribution to buoyancy arising from dhwt being a function of x and y
  do k=2,nk
    wrk1_v(:,:,k,1) = BDX_ET((Adv_vel%uh_et(:,:,k))*FAX(rho(:,:,k)))*Grd%dat(:,:)
    wrk1_v(:,:,k,2) = BDY_NT((Adv_vel%vh_nt(:,:,k))*FAY(rho(:,:,k)))*Grd%dat(:,:)
    fx = grav*rho0r
    do j=jsc,jec
      do i=isc,iec
        term = -fx*(wrk1_v(i,j,k,1) + wrk1_v(i,j,k,2))*Thickness%ztp(i,j,k)
        buoy = buoy + term
      enddo
    enddo
 enddo

  call mpp_sum( buoy)
  call mpp_sum( engint, size(engint(:)))
  call mpp_sum( engext, size(engext(:)))

  buoy      = buoy/ucell_volume
  engint(:) = engint(:)/ucell_volume
  engext(:) = engext(:)/ucell_volume

  plicin = engint(1) - engint(2) - engint(3) - engint(4) - engint(5) - engint(6) - engint(7) - engint(8) - engint(9)
  plicex = engext(1) - engext(2) - engext(3) - engext(4) - engext(5) - engext(6) - engext(7) - engext(8) - engext(9)
  buoerr = buoy - engint(6) - engext(6)
  enleak = engint(2) + engint(3) - engint(10) + engext(2) + engext(3) - engext(10)

  call kinetic_energy(Time, Thickness, Velocity, Dens, ke_tot, .false., .false.)
  call potential_energy(Time, Thickness, Dens, pe_tot, .false., .false.) 

  write (stdout(),'(//60x,a)') ' Globally averaged energy analysis  '

  if(acor > 0) then 
    write (stdout(),'(1x,/a)')  ' ==> Note: acor > 0 means Coriolis force will not conserve energy.'
  endif 
  if(Grd%tripolar) then 
    write (stdout(),'(/1x,/a)')  ' ==> NOTE: Energy conversion errors are known to be poor when using the tripolar grid.'
    write (stdout(),'(7x,a)')   'The problem is related ONLY to the energy diagnostic, not the model prognostic equations.'
    write (stdout(),'(7x,a)')   'It is caused by the diagnostic not having extra mpp calls to set redundancies at j=nj.'
  endif 
  write(stdout(),'(/,1x,a,e20.12)') 'Potential energy relative to first time (Joules) = ', pe_tot
  write(stdout(),'(1x,a,e20.12)') 'Kinetic energy  (Joules)                         = ', ke_tot
  write (stdout(),9100) 'Globally integrated U dot momentum eqns per unit volume (m^2/s^3 = W/kg)', ucell_volume, Grd%ucella(1)
  write (stdout(),9101) ' time rate of change ',  engint(1)+engext(1),engint(1), engext(1)
  write (stdout(),9101) ' horizontal advection',  engint(2)+engext(2),engint(2), engext(2)
  write (stdout(),9101) ' vertical advection  ',  engint(3)+engext(3),engint(3), engext(3)
  write (stdout(),9101) ' horizontal friction ',  engint(4)+engext(4),engint(4), engext(4)
  write (stdout(),9101) ' vertical friction   ',  engint(5)+engext(5),engint(5), engext(5)
  write (stdout(),9101) ' pressure forces     ',  engint(6)+engext(6),engint(6), engext(6)
  write (stdout(),9101) ' coriolis forces     ',  engint(9)+engext(9),engint(9), engext(9)
  write (stdout(),9101) ' work by wind        ',  engint(7)+engext(7),engint(7), engext(7)
  write (stdout(),9101) ' work by bottom drag ',  engint(8)+engext(8),engint(8), engext(8)
  write (stdout(),9101) ' imbalance           ',  plicin+plicex,      plicin,    plicex     
  write (stdout(),'(1x,/a)') ' Note the "imbalance" term is small only if run without baroclinc/barotropic split'  
  write (stdout(),9110) buoy, buoerr, enleak
  

9100 format(/1x,a,1x,/,' ocean volume       =',e16.9,' m^3',/, ' ocean surface area =',e16.9,&
          ' m^2',//,32x,'total                   internal mode              external mode')
9101 format(a21,3(3x,es24.17))
9102 format(a21,3x,es24.17)
9110 format(/ ' power per mass contributed by buoyancy: "w*rho*g"      =',es24.17,&
            /,' mismatch between "-u dot grad(p)" and  "w*rho*g"       =',es24.17,&
            /,' mismatch between "-u dot grad (v u)" and surface terms =',es24.17/)
  

end subroutine energy_analysis
! </SUBROUTINE> NAME="energy_analysis"



!#######################################################################
! <SUBROUTINE NAME="ocean_velocity_chksum">
!
! <DESCRIPTION>
! Compute checksum for velocity components 
! </DESCRIPTION>
subroutine ocean_velocity_chksum(Time, Velocity, index, write_advection)
 
  type(ocean_time_type), intent(in)     :: Time
  type(ocean_velocity_type), intent(in) :: Velocity
  integer, intent(in)                   :: index
  logical, intent(in), optional         :: write_advection

  real :: a, tmax, tmin, xtmax, xtmin, ytmax, ytmin, ztmax, ztmin, tmax0, tmin0
  integer :: i,j,k, itmax, jtmax, ktmax, itmin, jtmin, ktmin, n
  real :: fudge
  logical :: continue, writeadvection

  if (PRESENT(write_advection)) then 
    writeadvection = write_advection
  else 
    writeadvection = .true.
  endif 

  write(stdout(),*) '=== From ocean_velocity_mod (ocean_velocity_chksum): velocity checksum follows ==='
  call write_timestamp(Time%model_time)

  wrk1_v = 0.0

  wrk1_v(isc:iec,jsc:jec,:,1) = Velocity%u(isc:iec,jsc:jec,:,1,index)*Grd%umask(isc:iec,jsc:jec,:)
  wrk1_v(isc:iec,jsc:jec,:,2) = Velocity%u(isc:iec,jsc:jec,:,2,index)*Grd%umask(isc:iec,jsc:jec,:)

  write(stdout(),*) 'Zonal velocity chksum      = ',  mpp_chksum(wrk1_v(isc:iec,jsc:jec,:,1))
  write(stdout(),*) 'Meridional velocity chksum = ',  mpp_chksum(wrk1_v(isc:iec,jsc:jec,:,2))

  if(tendency==TWO_LEVEL .and. writeadvection) then 
    wrk1_v(isc:iec,jsc:jec,:,1) = Velocity%advection(isc:iec,jsc:jec,:,1,index)*Grd%umask(isc:iec,jsc:jec,:)
    wrk1_v(isc:iec,jsc:jec,:,2) = Velocity%advection(isc:iec,jsc:jec,:,2,index)*Grd%umask(isc:iec,jsc:jec,:)
    write(stdout(),*) 'Advection of u chksum = ',  mpp_chksum(wrk1_v(isc:iec,jsc:jec,:,1))
    write(stdout(),*) 'Advection of v chksum = ',  mpp_chksum(wrk1_v(isc:iec,jsc:jec,:,2))
  endif 

  return
  
end subroutine ocean_velocity_chksum
! </SUBROUTINE> NAME="ocean_velocity_chksum"


end module ocean_velocity_mod
