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
module ocean_advection_velocity_mod
!  
!<CONTACT EMAIL="Ronald.Pacanowski@noaa.gov"> R.C. Pacanowski
!</CONTACT>
!
!<CONTACT EMAIL="Stephen.Griffies@noaa.gov"> S.M. Griffies 
!</CONTACT>
!
!<OVERVIEW>
! Advection velocity components for tracer and momenta transport    
!</OVERVIEW>
!
!<DESCRIPTION>
! The module computes the horizontal and vertical components to the 
! advection velocity on the face of tracer and velocity cells.
! The three components are related by continuity. 
! The horizontal components to the advection velocity are
! thickness weighted.  Some diagnostics related to volume transport
! classified according to both depth and density classes are also computed.
! </DESCRIPTION>
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
! A Technical Guide to MOM4 (2003)
! </REFERENCE>
!
! <NOTE>
! The expressions for the horizontal components for tracer advection
! allow for a proper conversion between pressure work and buoyancy.
! </NOTE>
!
! <NOTE>
! The remapping operators are derived from considerations of linear
! interpolation and volume conservation.  A "remapping error" is
! computed to determine the consistency between the tracer and velocity
! grid advection velocities.  This error is roundoff only for cases
! where the horizontal tracer and velocity grids are linearly related,
! as is the case for the spherical coordinate version of mom4.  The
! tripolar version of mom4 does not have tracer and velocity grids
! related linearly, and so the "remapping error" is nontrivial.  The 
! significance of this error is unclear.  No adverse effects have been 
! identified. 
! </NOTE>
!
! <NOTE>
! The vertical velocity components for both the tracer and velocity cells
! are diagnosed via continuity (either volume or mass conservation depending
! on the use of the Boussinesq or non-Boussinesq versions of mom4).
! </NOTE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_advection_velocity_nml">
!
!  <DATA NAME="max_advection_velocity" UNITS="meter/sec" TYPE="real">
!  This is a check value used to determine if the time steps will result in 
!  linearly stable advection.  If set to a number < 0, then model will estimate the
!  value as a function of maximum grid size.
!  Note that this time step check is not rigorous, and it depends on the details of
!  the advection scheme.  Nonetheless, it provides some useful warning for setting the 
!  time steps in the model.  
!  </DATA> 
!
!</NAMELIST>

use axis_utils_mod,      only: nearest_index 
use constants_mod,       only: rho0r
use diag_manager_mod,    only: send_data, register_diag_field
use fms_mod,             only: write_version_number, FATAL, open_namelist_file, close_file, check_nml_error
use mpp_domains_mod,     only: mpp_update_domains
use mpp_mod,             only: mpp_error, mpp_chksum, stdout, stdlog
use mpp_mod,             only: mpp_min, mpp_max, mpp_error, mpp_pe

use ocean_domains_mod,   only: get_local_indices
use ocean_operators_mod, only: BAX, BAY, BDX_ET, BDY_NT, BDX_EU, BDY_NU
use ocean_operators_mod, only: REMAP_ET_TO_EU, REMAP_NT_TO_NU, REMAP_BT_TO_BU
use ocean_types_mod,     only: ocean_domain_type, ocean_grid_type
use ocean_types_mod,     only: ocean_grid_type, ocean_thickness_type
use ocean_types_mod,     only: ocean_adv_vel_type, ocean_velocity_type
use ocean_types_mod,     only: ocean_time_type, ocean_time_steps_type
use ocean_types_mod,     only: ocean_density_type, ocean_external_mode_type
use ocean_workspace_mod, only: wrk1 
use ocean_obc_mod,       only: ocean_obc_adjust_advel

implicit none

private

real    :: max_advection_velocity=-1.0 
integer :: unit=6

! for diagnostics 
integer :: id_w_bt=-1
integer :: id_uh_et=-1
integer :: id_vh_nt=-1
integer :: id_w_bu=-1
integer :: id_uh_eu=-1
integer :: id_vh_nu=-1
integer :: id_courant_uet=-1
integer :: id_courant_vnt=-1
integer :: id_courant_wbt=-1
integer :: id_courant_ueu=-1
integer :: id_courant_vnu=-1
integer :: id_courant_wbu=-1

real    :: dtts
real    :: dtuv

logical :: used

#include <ocean_memory.h>

type(ocean_domain_type), pointer :: Dom =>NULL()
type(ocean_grid_type), pointer   :: Grd =>NULL()

public ocean_advection_velocity_init
public ocean_advection_velocity

character(len=128) :: version=&
     '$Id$'
character (len=128) :: tagname = &
     '$Name$'

logical :: have_obc                 = .false.
logical :: debug_advection_velocity = .false.
logical :: module_is_initialized    = .FALSE.

namelist /ocean_advection_velocity_nml/ debug_advection_velocity, max_advection_velocity

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_advection_velocity_init">
!
! <DESCRIPTION>
! Initialize the advection velocity module
! </DESCRIPTION>
!
subroutine ocean_advection_velocity_init(Grid, Domain, Time, Time_steps, Adv_vel, obc, debug)

  type(ocean_grid_type), target, intent(in)       :: Grid
  type(ocean_domain_type), target, intent(in)     :: Domain  
  type(ocean_time_type), target, intent(in)       :: Time
  type(ocean_time_steps_type), intent(in)         :: Time_steps
  type(ocean_adv_vel_type), intent(inout)         :: Adv_vel
  logical, intent(in)                             :: obc
  logical, intent(in), optional                   :: debug

  integer :: ioun, io_status, ierr
  real    :: gridmin, gridmax, dtadv, a, b
  integer :: i, j, icg, jcg
  real    :: vertical_factor = 2.0
  real    :: max_dt_for_advection
  real    :: max_dt_for_advection0
  real, dimension(7) :: res = (/0.0625, 0.125, 0.25, 0.50, 1.00, 2.0, 4.0/) ! max resolution (deg)
  real, dimension(7) :: vel = (/2.0,    1.8,   1.6,  1.45, 1.25, 0.9, 0.4/) ! estimated max advection velocity (m/sec)

  if ( module_is_initialized ) call mpp_error(FATAL,'==>Error from ocean_advection_velocity_mod: module already initialized.')

  call write_version_number( version, tagname )

  module_is_initialized = .TRUE.

  if (PRESENT(debug)) debug_advection_velocity = debug
  have_obc       = obc

  ! provide for namelist over-ride of defaults 
  ioun = open_namelist_file ()
  read  (ioun, ocean_advection_velocity_nml,iostat=io_status)
  write (stdout(),'(/)')
  write (stdout(), ocean_advection_velocity_nml)
  write (stdlog(), ocean_advection_velocity_nml)
  ierr = check_nml_error(io_status,'ocean_advection_velocity_nml')
  call close_file(ioun)

  Dom => Domain
  Grd => Grid

  dtts = Time_steps%dtts
  dtuv = Time_steps%dtuv
  
#ifndef STATIC_MEMORY
  call get_local_indices(Dom, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk = Grd%nk
 
  ! allocate and initialze variables

  allocate ( Adv_vel%uh_et(isd:ied,jsd:jed,nk) )
  allocate ( Adv_vel%vh_nt(isd:ied,jsd:jed,nk) )
  allocate ( Adv_vel%uh_eu(isd:ied,jsd:jed,nk) )
  allocate ( Adv_vel%vh_nu(isd:ied,jsd:jed,nk) )
  allocate ( Adv_vel%w_bt(isd:ied,jsd:jed,0:nk) )
  allocate ( Adv_vel%w_bu(isd:ied,jsd:jed,0:nk) )

#endif

  Adv_vel%uh_et(:,:,:)   = 0.0
  Adv_vel%vh_nt(:,:,:)   = 0.0
  Adv_vel%w_bt(:,:,:)    = 0.0
  Adv_vel%uh_eu(:,:,:)   = 0.0
  Adv_vel%vh_nu(:,:,:)   = 0.0
  Adv_vel%w_bu(:,:,:)    = 0.0

  ! estimate the max advective speed for the given resolution

  if (max_advection_velocity < 0.0) then
    gridmax = 0.0
    do j=jsc,jec
      do i=isc,iec
        if (Grd%kmt(i,j) /= 0) then
          gridmax = max(gridmax,Grd%dyt(i,j),Grd%dxt(i,j))
        endif
      enddo
    enddo
    call mpp_max (gridmax)
    gridmax = gridmax/111.324e3 ! convert to degrees
    i = nearest_index(gridmax, res)
    if (gridmax <= res(1)) then
      max_advection_velocity =  vel(1)
    elseif (gridmax >= res(7)) then
      max_advection_velocity =  vel(7)
    elseif (gridmax < res(i)) then
      a = res(i)-gridmax
      b = gridmax-res(i-1)
      max_advection_velocity = (a*vel(i-1) + b*vel(i))/(a+b)      
    elseif (gridmax > res(i)) then
      a = gridmax-res(i)
      b = res(i+1)-gridmax
      max_advection_velocity = (a*vel(i+1) + b*vel(i))/(a+b)      
    endif
    write (stdout(),'(/a,f5.2,a/a/)')' Note: The max_advection_velocity is estimated to be ',max_advection_velocity,&
                   ' m/s based on grid resolution.','       A more appropriate value can be specified via namelist.'
  endif

  icg = isc; jcg = jsc; gridmin = 1.0e20; max_dt_for_advection = 1.e6; gridmax = 0.0
  do j=jsc,jec
    do i=isc,iec
      if (Grd%kmt(i,j) /= 0) then
        gridmin = min(gridmin,Grd%dyt(i,j),Grd%dxt(i,j))
        gridmax = max(gridmax,Grd%dyt(i,j),Grd%dxt(i,j))
        dtadv = 0.5*gridmin/max_advection_velocity
        if (dtadv < max_dt_for_advection) then
          max_dt_for_advection = dtadv; icg  = i; jcg  = j
        endif
      endif
    enddo
  enddo
  
  max_dt_for_advection = nint(max_dt_for_advection/vertical_factor) ! Questions? -> rcp@gfdl.noaa.gov
  max_dt_for_advection =  max_dt_for_advection + 0.001*mpp_pe()          ! to separate redundancies
  max_dt_for_advection0 = max_dt_for_advection
  call mpp_min (max_dt_for_advection)
  call mpp_max (gridmax)
  
  if (max_dt_for_advection == max_dt_for_advection0) then
    if (max(dtts,dtuv) <= max_dt_for_advection) then
      write (unit,'(/a,i4,a,i4,a,f6.2,a,f6.2,a)')' Note: Advection stability most nearly violated at T-cell (i,j) = (',&
      icg,',',jcg,'), (lon,lat) = (',Grd%xt(icg,jcg),',',Grd%yt(icg,jcg),').'
    else
      write (unit,'(/a,i4,a,i4,a,f6.2,a,f6.2,a)')'=>Error: Advection stability violated at T-cell (i,j) = (',&
      icg,',',jcg,'), (lon,lat) = (',Grd%xt(icg,jcg),',',Grd%yt(icg,jcg),').'
      call mpp_error(FATAL,'==>Error from ocean_advection_velocity_mod: tracer advection stability is violated.')
    endif
    write (unit,'(a,f7.2,a/a,f10.0,a,f10.0,a)')&
     '         due to a specified max_advection_velocity of ',max_advection_velocity,' m/s.',&
     '         "max(dtts,dtuv)" must be less than ',max_dt_for_advection,' sec. "max(dtts,dtuv)" = ', max(dtts,dtuv),' sec.'
  endif
  max_dt_for_advection = nint(max_dt_for_advection)

  ! register advective velocity components for diagnostic output
  id_w_bt = register_diag_field ('ocean_model', 'wt', Grd%vel_axes_wt(1:3), Time%model_time, &
    'vertical velocity T-points', 'm/sec', missing_value=-10.0, range=(/-10.0,10.0/))
  id_uh_et = register_diag_field ('ocean_model', 'uh_et', Grd%tracer_axes_flux_x(1:3), Time%model_time, &
    'uh_et on horz face of T-cells', 'm^2/sec', missing_value=-100.0, range=(/-100.0,100.0/))
  id_vh_nt = register_diag_field ('ocean_model', 'vh_nt', Grd%tracer_axes_flux_y(1:3), Time%model_time, &
    'vh_nt on horz face of T-cells', 'm^2/sec', missing_value=-100.0, range=(/-100.0,100.0/))
  id_w_bu = register_diag_field ('ocean_model', 'wu', Grd%vel_axes_wu(1:3), Time%model_time, &
    'vertical velocity U-points', 'm/sec', missing_value=-10.0, range=(/-10.0,10.0/))
  id_uh_eu = register_diag_field ('ocean_model', 'uh_eu', Grd%tracer_axes(1:3), Time%model_time, &
    'uh_eu on horz face of U-cells', 'm^2/sec', missing_value=-100.0, range=(/-100.0,100.0/))
  id_vh_nu = register_diag_field ('ocean_model', 'vh_nu', Grd%tracer_axes(1:3), Time%model_time, &
    'vh_nu on horz face of U-cells', 'm^2/sec', missing_value=-100.0, range=(/-100.0,100.0/))

  ! register Courant numbers for diagnostic output
  id_courant_uet = register_diag_field ('ocean_model', 'courant_uet', Grd%tracer_axes(1:3), Time%model_time, &
    'Courant number [uet*dtts/dxt]', 'none', missing_value=-1.e10, range=(/-1.e10,1.e10/))
  id_courant_vnt = register_diag_field ('ocean_model', 'courant_vnt', Grd%tracer_axes(1:3), Time%model_time, &
    'Courant number [vnt*dtts/dyt]', 'none', missing_value=-1.e10, range=(/-1.e10,1.e10/))
  id_courant_wbt = register_diag_field ('ocean_model', 'courant_wbt', Grd%tracer_axes(1:3), Time%model_time, &
    'Courant number [wbt*dtts/dht]', 'none', missing_value=-1.e10, range=(/-1.e10,1.e10/))
  id_courant_ueu = register_diag_field ('ocean_model', 'courant_ueu', Grd%vel_axes_uv(1:3), Time%model_time, &
    'Courant number [ueu*dtuv/dxu]', 'none', missing_value=-1.e10, range=(/-1.e10,1.e10/))
  id_courant_vnu = register_diag_field ('ocean_model', 'courant_vnu', Grd%vel_axes_uv(1:3), Time%model_time, &
    'Courant number [vnu*dtuv/dyu]', 'none', missing_value=-1.e10, range=(/-1.e10,1.e10/))
  id_courant_wbu = register_diag_field ('ocean_model', 'courant_wbu', Grd%vel_axes_uv(1:3), Time%model_time, &
    'Courant number [wnu*dtuv/dhu]', 'none', missing_value=-1.e10, range=(/-1.e10,1.e10/))


end subroutine ocean_advection_velocity_init
! </SUBROUTINE> NAME="ocean_advection_velocity_init"


!#######################################################################
! <SUBROUTINE NAME="ocean_advection_velocity">
!
! <DESCRIPTION>
! Compute thickness weighted advection velocities normal to
! side faces of U-cells and T-cells on a b-grid. Also compute
! vertical advection velocity normal to bottom faces of U-cells and T-cells. 
! </DESCRIPTION>
!
subroutine ocean_advection_velocity (Ext_mode, Velocity, Time, Thickness, Adv_vel)
  
  type(ocean_external_mode_type), intent(in) :: Ext_mode
  type(ocean_velocity_type), intent(in)      :: Velocity
  type(ocean_time_type), intent(in)          :: Time
  type(ocean_thickness_type), intent(in)     :: Thickness
  type(ocean_adv_vel_type), intent(inout)    :: Adv_vel

  integer ::  i, j, k, tau

  if ( .not. module_is_initialized ) &
   call mpp_error(FATAL, '==>Error from ocean_advection_velocity_mod (ocean_advection_velocity): module must be initialized')

  tau  = Time%tau
  wrk1 = 0.0

  ! compute thickness weighted C-grid advection velocity (uh_et, vh_nt) 
  ! and the corresponding advection velocity (uh_eu,vh_nu) to advect B-grid velocity 
  do k= 1, nk
     Adv_vel%uh_et(:,:,k) = BAY(Velocity%u(:,:,k,1,tau)*Grd%dyu(:,:)*Thickness%dhu(:,:,k,tau))/Grd%dyte(:,:) 
     Adv_vel%vh_nt(:,:,k) = BAX(Velocity%u(:,:,k,2,tau)*Grd%dxu(:,:)*Thickness%dhu(:,:,k,tau))/Grd%dxtn(:,:) 
     Adv_vel%uh_eu(:,:,k) = REMAP_ET_TO_EU(Adv_vel%uh_et(:,:,k))                                       
     Adv_vel%vh_nu(:,:,k) = REMAP_NT_TO_NU(Adv_vel%vh_nt(:,:,k))                                       
  enddo

  ! w_bt(k=0) is a placeholder for convud_t.  It should not be confused with the 
  ! vertical advection velocity at the ocean surface, which is -(pme + river).
  ! Defining w_bt(k=0) = convud_t allows for cleaner definitions of the vertical
  ! advection velocities at depth, as well as for some diagnostic purposes. 
  Adv_vel%w_bt(:,:,0) = Grd%tmask(:,:,1)*Ext_mode%convud_t(:,:,tau)     
  call mpp_update_domains (Adv_vel%w_bt(:,:,0), Dom%domain2d)
  Adv_vel%w_bu(:,:,0) = Grd%umask(:,:,1)*REMAP_BT_TO_BU(Adv_vel%w_bt(:,:,0))
  do k=1,nk

     Adv_vel%w_bt(:,:,k) = Adv_vel%w_bt(:,:,k-1)  + &
          BDX_ET(Adv_vel%uh_et(:,:,k)) + &
          BDY_NT(Adv_vel%vh_nt(:,:,k)) 

     Adv_vel%w_bu(:,:,k) = Adv_vel%w_bu(:,:,k-1)  + &
          BDX_EU(Adv_vel%uh_eu(:,:,k)) + &
          BDY_NU(Adv_vel%vh_nu(:,:,k)) 

  enddo


  if (have_obc) then
     call ocean_obc_adjust_advel(Adv_vel) 
   endif

  if(debug_advection_velocity) then 
      write(stdout(),*) 'Adv_vel%uh_et   ==> ', mpp_chksum(Adv_vel%uh_et(isc:iec,jsc:jec,:))
      write(stdout(),*) 'Adv_vel%vh_nt   ==> ', mpp_chksum(Adv_vel%vh_nt(isc:iec,jsc:jec,:))
      write(stdout(),*) 'Adv_vel%uh_eu   ==> ', mpp_chksum(Adv_vel%uh_eu(isc:iec,jsc:jec,:))
      write(stdout(),*) 'Adv_vel%vh_nu   ==> ', mpp_chksum(Adv_vel%vh_nu(isc:iec,jsc:jec,:))
      write(stdout(),*) 'Adv_vel%w_bt    ==> ', mpp_chksum(Adv_vel%w_bt(isc:iec,jsc:jec,:))
      write(stdout(),*) 'Adv_vel%w_bu    ==> ', mpp_chksum(Adv_vel%w_bu(isc:iec,jsc:jec,:))
      write(stdout(),*) 'Adv_vel%w_bt(0) ==> ', mpp_chksum(Adv_vel%w_bt(isc:iec,jsc:jec,0))
      write(stdout(),*) 'Adv_vel%w_bu(0) ==> ', mpp_chksum(Adv_vel%w_bu(isc:iec,jsc:jec,0))
  endif

  ! send advective velocity components to diagnostic manager 
  if (id_uh_et > 0) used = send_data (id_uh_et, Adv_vel%uh_et(:,:,:), &
                           Time%model_time, rmask=Grd%tmask(:,:,:), &
                           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  if (id_vh_nt > 0) used = send_data (id_vh_nt, Adv_vel%vh_nt(:,:,:), &
                           Time%model_time, rmask=Grd%tmask(:,:,:), &
                           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  if (id_w_bt > 0) used =  send_data (id_w_bt, Adv_vel%w_bt(:,:,1:nk), &
                           Time%model_time, rmask=Grd%tmask(:,:,:),&
                           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  if (id_uh_eu > 0) used = send_data (id_uh_eu, Adv_vel%uh_eu(:,:,:), &
                           Time%model_time, rmask=Grd%umask(:,:,:), &
                           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  if (id_vh_nu > 0) used = send_data (id_vh_nu, Adv_vel%vh_nu(:,:,:), &
                           Time%model_time, rmask=Grd%umask(:,:,:),&
                           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  if (id_w_bu > 0) used =  send_data (id_w_bu, Adv_vel%w_bu(:,:,1:nk), &
                           Time%model_time, rmask=Grd%umask(:,:,:), &
                           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

  ! send Courant numbers to diagnostic manager 
  if (id_courant_uet > 0) then 
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk1(i,j,k) = dtts*Adv_vel%uh_et(i,j,k)/(Thickness%dht(i,j,k,tau)*Grd%dxt(i,j))
            enddo
         enddo
      enddo
      used = send_data (id_courant_uet, wrk1(:,:,:), &
             Time%model_time, rmask=Grd%tmask(:,:,:), &
             is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif
  if (id_courant_vnt > 0) then 
      do k=1,nk  
         do j=jsc,jec
            do i=isc,iec
               wrk1(i,j,k) = dtts*Adv_vel%vh_nt(i,j,k)/(Thickness%dht(i,j,k,tau)*Grd%dyt(i,j))
            enddo
         enddo
      enddo
      used = send_data (id_courant_vnt, wrk1(:,:,:), &
             Time%model_time, rmask=Grd%tmask(:,:,:), &
             is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif
  if (id_courant_wbt > 0) then 
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk1(i,j,k) = dtts*Adv_vel%w_bt(i,j,k)/Thickness%dht(i,j,k,tau)
            enddo
         enddo
      enddo
      used = send_data (id_courant_wbt, wrk1(:,:,:), &
             Time%model_time, rmask=Grd%tmask(:,:,:), &
             is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif
  if (id_courant_ueu > 0) then 
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk1(i,j,k) = dtuv*Adv_vel%uh_eu(i,j,k)/(Thickness%dhu(i,j,k,tau)*Grd%dxu(i,j))
            enddo
         enddo
      enddo
      used = send_data (id_courant_ueu, wrk1(:,:,:), &
             Time%model_time, rmask=Grd%umask(:,:,:), &
             is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif 
  if (id_courant_vnu > 0) then 
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk1(i,j,k) = dtuv*Adv_vel%vh_nu(i,j,k)/(Thickness%dhu(i,j,k,tau)*Grd%dyu(i,j))
            enddo
         enddo
      enddo
      used = send_data (id_courant_ueu, wrk1(:,:,:), &
             Time%model_time, rmask=Grd%umask(:,:,:), &
             is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif
  if (id_courant_wbu > 0) then 
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               wrk1(i,j,k) = dtuv*Adv_vel%w_bu(i,j,k)/Thickness%dhu(i,j,k,tau)
            enddo
         enddo
      enddo
      used = send_data (id_courant_wbu, wrk1(:,:,:), &
             Time%model_time, rmask=Grd%umask(:,:,:), &
             is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif


end subroutine ocean_advection_velocity
! </SUBROUTINE> NAME="advection_velocity"


end module ocean_advection_velocity_mod
