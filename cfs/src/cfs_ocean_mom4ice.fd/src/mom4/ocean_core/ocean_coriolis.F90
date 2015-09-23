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
module ocean_coriolis_mod
!  
!<CONTACT EMAIL="Ronald.Pacanowski@noaa.gov"> R.C. Pacanowski
!</CONTACT>
!
!<CONTACT EMAIL="Tony.Rosati@noaa.gov"> A. Rosati 
!</CONTACT>
!
!<OVERVIEW>
! Coriolis acceleration 
!</OVERVIEW>

!<DESCRIPTION>
! This module computes Coriolis acceleration on a B-grid. 
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
! A Technical Guide to MOM4 (2003)
! </REFERENCE>
!
! </INFO>
!

use diag_manager_mod, only: register_static_field, register_diag_field, send_data
use fms_mod,          only: write_version_number, mpp_error, FATAL
use mpp_mod,          only: mpp_max, stdout, stdlog

use ocean_domains_mod, only: get_local_indices
use constants_mod,     only: radius, radian, pi, epsln, omega
use ocean_types_mod,   only: ocean_grid_type, ocean_domain_type
use ocean_types_mod,   only: ocean_time_type, ocean_time_steps_type
use ocean_types_mod,   only: ocean_velocity_type, ocean_thickness_type

implicit none

private


! parameter that determines how to time discretize Coriolis force 
real :: acor

! for diagnostics 
integer :: id_coriolis=-1
integer :: id_beta=-1
integer :: id_beta_eff=-1
integer :: id_cor_u=-1
integer :: id_cor_v=-1
integer :: id_h_cor_u=-1
integer :: id_h_cor_v=-1
logical :: used

#include <ocean_memory.h>

public coriolis
public ocean_coriolis_init

type(ocean_grid_type), pointer :: Grd =>NULL()

interface coriolis
 module procedure coriolis_field
 module procedure coriolis_slice
end interface

character(len=128) :: version = &
     '$Id$'
character (len=128) :: tagname = &
     '$Name$'

logical :: module_is_initialized = .FALSE.

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_coriolis_init">
!
! <DESCRIPTION>
! Initialize the Coriolis module 
! </DESCRIPTION>
!
subroutine ocean_coriolis_init(Grid, Domain, Time, Time_steps)

  type(ocean_grid_type), intent(inout), target :: Grid
  type(ocean_domain_type), intent(in), target  :: Domain  
  type(ocean_time_type), intent(in)            :: Time
  type(ocean_time_steps_type), intent(in)      :: Time_steps

  real    :: deg2m, sin1, y, fmax
  real    :: max_dt_for_inertial_oscillation
  integer :: i, j
  integer :: ioun, io_status, ierr

  if ( module_is_initialized ) then 
    call mpp_error(FATAL, '==>Error in ocean_coriolis_mod (ocean_coriolis_init): module already initialized.')
  endif 

  module_is_initialized = .TRUE.

  call write_version_number( version, tagname )

  Grd=> Grid

#ifndef STATIC_MEMORY
  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk = Grid%nk

  allocate (Grid%f(isd:ied,jsd:jed))
  allocate (Grid%beta(isd:ied,jsd:jed))
  allocate (Grid%beta_eff(isd:ied,jsd:jed))
#endif

  ! parameter determining how to time discretize Coriolis force 
  acor = Time_steps%acor

  if (Grid%beta_plane .or. Grid%f_plane) then

    ! beta plane with f = f0 + beta*y where f0 is at f_plane_latitude
    ! if f_plane then beta = 0

    if (Grid%f_plane) then
      Grid%beta(:,:) = 0.0
    else
      Grid%beta(:,:) = 2.0*omega*cos(Grid%f_plane_latitude*pi/180.0)/radius
    endif
    deg2m = radius/radian
    sin1  = sin(Grid%f_plane_latitude*pi/180.0)
    do j=jsd,jed
      do i=isd,ied
        y      = (Grid%yu(i,j)-Grid%f_plane_latitude)*deg2m
        Grid%f(i,j) = 2.0*omega*sin1  + y*Grid%beta(i,j)
      enddo
    enddo
    if (Grid%f_plane) then
      write (stdout(),'(//,a,f6.2,a//)') ' Note: "f plane" approximation set using f0 at latitude =', Grid%f_plane_latitude,'deg'
    else
      write (stdout(),'(//,a,f6.2,a,g14.7//)') ' Note: "beta plane" approximation set using f0 at latitude=',&
                                             Grid%f_plane_latitude,'deg and beta =',Grid%beta(isc,jsc)
    endif
  else
    Grid%f(:,:)    = 2.0*omega*sin(Grid%phiu(:,:))
    Grid%beta(:,:) = 2.0*omega*cos(Grid%phiu(:,:))/radius
  endif

  ! beta_eff centered on u-cell
  ! beta_eff = H*|grad(f/H)| --> |(f/hu)*dht_dx| + |beta-(f/hu)*dht_dy|
  Grid%beta_eff(:,:) = 0.0
  do j=jsc,jec
    do i=isc,iec
      Grid%beta_eff(i,j) = sqrt( (Grid%f(i,j)/(epsln+Grid%hu(i,j))*Grid%dht_dx(i,j))**2 &
                                +(Grid%beta(i,j)-Grid%f(i,j)/(epsln+Grid%hu(i,j))*Grid%dht_dy(i,j))**2 )
    enddo
  enddo 

  ! check for marginally resolved inertial oscillation

  fmax = 2.0*omega*0.0002 ! ~ 0.01 deg latitude
  do j=jsc,jec
    do i=isc,iec
      if (Grid%kmu(i,j) /= 0) then
        fmax = max(fmax,abs(Grid%f(i,j)))
      endif
    enddo
  enddo
  call mpp_max (fmax)
  max_dt_for_inertial_oscillation = int(1.0/fmax) ! assume 2*pi timesteps to resolve oscillation
  write (stdout(),'(/1x, a,f8.0,a)')&
  ' ==> Note: 2*pi timesteps/(min inertial period) implies a maximum dtuv =',max_dt_for_inertial_oscillation,' sec.'
  if(acor == 0.0) then 
    write (stdout(),'(6x,a)') 'Coriolis force is treated explicitly in time since the acor parameter is set to zero.'
  endif 
  if(acor > 0.0 .and. acor < 0.5) then 
    call mpp_error(FATAL,'==> Error in ocean_coriolis_mod: "acor" must be set either to acor=0.0 or 0.5 .le. acor .le. 1.0')
  endif 
  if(acor .ge. 0.5 .and. acor .le. 1.0) then 
    write (stdout(),'(6x,a,f6.3)') ' ==> Note: Coriolis force is treated semi-implicitly in time with acor set to ',acor
    write (stdout(),'(6x,a)') 'This approach removes the constraint on dtuv associated with inertial oscillations.'
  endif 
  if (Time_steps%dtuv > max_dt_for_inertial_oscillation .and. acor==0.0) then
    call mpp_error(FATAL,&
    '==> Error in ocean_coriolis_mod: inertial oscillation not resolved. Reduce "dtuv" or set 0.5 .le. acor .le. 1.0.')
  endif

  ! diagnostic manager registers and sends for static fields 

  id_coriolis  = register_static_field ('ocean_model', 'f_coriolis', Grd%vel_axes_uv(1:2), &
                                        'Coriolis frequency', '1/s',&
                                         missing_value=-10.0, range=(/-10.0,10.0/))
  if (id_coriolis > 0) used = send_data (id_coriolis, Grid%f(:,:), &
                              Time%model_time, rmask=Grd%umask(:,:,1), &
                              is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

  id_beta  = register_static_field ('ocean_model', 'beta', Grd%vel_axes_uv(1:2), &
                                    'planetary beta', '1/(m*s)',&
                                    missing_value=-10.0, range=(/-10.0,10.0/))
  if (id_beta > 0) used = send_data (id_beta, Grd%beta(:,:), &
                          Time%model_time, rmask=Grd%umask(:,:,1), &
                          is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

  id_beta_eff  = register_static_field ('ocean_model', 'beta_eff', Grd%vel_axes_uv(1:2), &
                                        'effective beta', '1/(m*s)',&
                                        missing_value=-10.0, range=(/-10.0,10.0/))
  if (id_beta_eff > 0) used = send_data (id_beta_eff, Grd%beta_eff(:,:), &
                              Time%model_time, rmask=Grd%umask(:,:,1), &
                              is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

  ! diagnostic manager registers for dynamic fields 
  id_cor_u =  register_diag_field ('ocean_model', 'cor_u', Grd%vel_axes_uv(1:3), Time%model_time, &
     'explicit coriolis accel in i-direct', 'm/s^2', missing_value=-1e10, range=(/-1e10,1e10/))
  id_cor_v =  register_diag_field ('ocean_model', 'cor_v', Grd%vel_axes_uv(1:3), Time%model_time, &
     'explicit coriolis accel in j-direct', 'm/s^2', missing_value=-1e10, range=(/-1e10,1e10/))
  id_h_cor_u =  register_diag_field ('ocean_model', 'h_cor_u', Grd%vel_axes_uv(1:3), Time%model_time, &
     'Thickness wghtd explicit coriolis accel in i-direct', 'm^2/s^2', missing_value=-1e10, range=(/-1e10,1e10/))
  id_h_cor_v =  register_diag_field ('ocean_model', 'h_cor_v', Grd%vel_axes_uv(1:3), Time%model_time, &
     'Thickness wghtd explicit coriolis accel in j-direct', 'm^2/s^2', missing_value=-1e10, range=(/-1e10,1e10/))


end subroutine ocean_coriolis_init
! </SUBROUTINE> NAME="ocean_coriolis_init"


!#######################################################################
! <FUNCTION NAME="coriolis_field">
!
! <DESCRIPTION>
! Compute thickness weighted acceleration due to Coriolis term on a B-grid
! </DESCRIPTION>
!
function coriolis_field(Time, Velocity, Thickness, diag_flag, full_flag)

  type(ocean_time_type), intent(in)      :: Time
  type(ocean_velocity_type), intent(in)  :: Velocity
  type(ocean_thickness_type), intent(in) :: Thickness
  logical, intent(in), optional          :: diag_flag 
  logical, intent(in), optional          :: full_flag 

  real, dimension(isc:iec,jsc:jec,nk,2) :: coriolis_field 
  integer                               :: k, itime, tau, taum1, taup1
  logical                               :: send_diagnostics 
  logical                               :: full_coriolis=.false.
  
  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, '==>Error in ocean_coriolis_mod (coriolis_field): module must be initialized')
  endif 

  tau   = Time%tau
  taum1 = Time%taum1
  taup1 = Time%taup1
  if(acor == 0.0) then 
    itime = tau
  else 
    itime = taum1
  endif 

  ! assume send_diagnostics=.true., unless diag_flag says it is false.
  if (PRESENT(diag_flag)) then 
    send_diagnostics = diag_flag
  else 
    send_diagnostics = .true.
  endif 

  ! assume full_flag=.false., unless diag_flag says it is true.
  if(acor > 0) then 
    if (PRESENT(full_flag)) then 
      full_coriolis = full_flag
    else 
      full_coriolis = .false.
    endif 
  endif 

  ! when acor > 0, Coriolis force dissipates energy.  To compute the 
  ! amount of dissipation in the energy analysis, we need to compute 
  ! the Coriolis force as given here.  
  if(full_coriolis) then 
    do k=1,nk
      coriolis_field(isc:iec,jsc:jec,k,1) =  Grd%f(isc:iec,jsc:jec)*( &
           Velocity%u(isc:iec,jsc:jec,k,2,taum1)*(1.0-acor) + Velocity%u(isc:iec,jsc:jec,k,2,taup1))
      coriolis_field(isc:iec,jsc:jec,k,2) = -Grd%f(isc:iec,jsc:jec)*( &
           Velocity%u(isc:iec,jsc:jec,k,1,itime)*(1.0-acor)  + Velocity%u(isc:iec,jsc:jec,k,1,taup1))
    enddo
  else 
    do k=1,nk
      coriolis_field(isc:iec,jsc:jec,k,1) =  Grd%f(isc:iec,jsc:jec)*Velocity%u(isc:iec,jsc:jec,k,2,itime)
      coriolis_field(isc:iec,jsc:jec,k,2) = -Grd%f(isc:iec,jsc:jec)*Velocity%u(isc:iec,jsc:jec,k,1,itime)
    enddo
  endif 
 
  if(send_diagnostics) then 
     if (id_cor_u > 0) used = send_data( id_cor_u, coriolis_field(isc:iec,jsc:jec,:,1), &
                              Time%model_time, rmask=Grd%umask(isc:iec,jsc:jec,:))
     if (id_cor_v > 0) used = send_data( id_cor_v, coriolis_field(isc:iec,jsc:jec,:,2), &
                              Time%model_time, rmask=Grd%umask(isc:iec,jsc:jec,:))
  endif 

  ! thickness weight the acceleration 
  do k=1,nk
     coriolis_field(isc:iec,jsc:jec,k,1) = coriolis_field(isc:iec,jsc:jec,k,1)*Thickness%dhu(isc:iec,jsc:jec,k,tau)
     coriolis_field(isc:iec,jsc:jec,k,2) = coriolis_field(isc:iec,jsc:jec,k,2)*Thickness%dhu(isc:iec,jsc:jec,k,tau)
  enddo  

  if(send_diagnostics) then 
     if (id_h_cor_u > 0) used = send_data( id_h_cor_u, coriolis_field(isc:iec,jsc:jec,:,1), &
                                Time%model_time, rmask=Grd%umask(isc:iec,jsc:jec,:))
     if (id_h_cor_v > 0) used = send_data( id_h_cor_v, coriolis_field(isc:iec,jsc:jec,:,2), &
                                Time%model_time, rmask=Grd%umask(isc:iec,jsc:jec,:))
  endif 


end function coriolis_field
! </FUNCTION> NAME="coriolis_field"


!#######################################################################
! <FUNCTION NAME="coriolis_slice">
!
! <DESCRIPTION>
! Compute acceleration due to Coriolis term on a B-grid for a particular
! depth level and particular direction.  
! </DESCRIPTION>
!
function coriolis_slice(Time, Velocity, Thickness, k, n)
  
  type(ocean_time_type), intent(in)      :: Time
  type(ocean_velocity_type), intent(in)  :: Velocity
  type(ocean_thickness_type), intent(in) :: Thickness
  integer, intent(in)                    :: k, n
  real, dimension(isc:iec,jsc:jec)       :: coriolis_slice 

  integer          :: itime, tau
  character(len=4) :: char_n

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, '==>Error in ocean_coriolis_mod (coriolis_slice): module must be initialized')
  endif 

  tau = Time%tau
 
  if(acor == 0.0) then 
    itime = Time%tau
  else 
    itime = Time%taum1
  endif 

  if (n == 1) then
    coriolis_slice(isc:iec,jsc:jec) =  Grd%f(isc:iec,jsc:jec)*Velocity%u(isc:iec,jsc:jec,k,2,itime)
  elseif (n == 2) then
    coriolis_slice(isc:iec,jsc:jec) = -Grd%f(isc:iec,jsc:jec)*Velocity%u(isc:iec,jsc:jec,k,1,itime)
  else
    write( char_n,'(i4)' ) n
    call mpp_error(FATAL, '==>Error in coriolis_slice with n='//trim(char_n) )
  endif

  ! thickness weight the acceleration 
  coriolis_slice(isc:iec,jsc:jec) = coriolis_slice(isc:iec,jsc:jec)*Thickness%dhu(isc:iec,jsc:jec,k,tau)


end function coriolis_slice
! </FUNCTION> NAME="coriolis_slice"


end module ocean_coriolis_mod
