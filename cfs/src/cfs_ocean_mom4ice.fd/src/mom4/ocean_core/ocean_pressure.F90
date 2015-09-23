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
module ocean_pressure_mod
!  
!<CONTACT EMAIL="Ronald.Pacanowski@noaa.gov"> R.C. Pacanowski
!</CONTACT>
!
! <REVIEWER EMAIL="Tony.Rosati@noaa.gov">
! A. Rosati 
! </REVIEWER>
!
! <REVIEWER EMAIL="Stephen.Griffies@noaa.gov">
! S.M. Griffies 
! </REVIEWER>
!
!<OVERVIEW>
! Hydrostatic pressure and its gradient
!</OVERVIEW>
!
!<DESCRIPTION>
! This module computes hydrostatic pressure and its horizontal gradient for use in 
! forcing the linear momentum.
!</DESCRIPTION>
!
!<NAMELIST NAME="ocean_pressure_nml">
!  <DATA NAME="debug_pressure" TYPE="logical">
!  For debugging. 
!  </DATA> 
!
!</NAMELIST>

use constants_mod,    only: rho0, grav, c2dbars, rho0r 
use diag_manager_mod, only: register_diag_field, send_data
use fms_mod,          only: open_namelist_file, check_nml_error, close_file, write_version_number
use mpp_io_mod,       only: mpp_open, mpp_close, MPP_RDONLY, MPP_ASCII
use mpp_mod,          only: mpp_error, FATAL, stdout, stdlog

use ocean_domains_mod,   only: get_local_indices
use ocean_operators_mod, only: FAY, FAX, FDX_NT, FDY_ET 
use ocean_types_mod,     only: ocean_domain_type, ocean_grid_type, ocean_time_type
use ocean_types_mod,     only: ocean_velocity_type, ocean_thickness_type 
use ocean_workspace_mod, only: wrk1, wrk2
use ocean_obc_mod,       only: store_ocean_obc_pressure_grad

implicit none

private

public ocean_pressure_init
public pressure_in_dbars
public pressure_gradient
private hydrostatic_pressure

#include <ocean_memory.h>

#ifdef STATIC_MEMORY
  real, dimension(isd:ied,jsd:jed,2)  :: tmp
#else
  real, allocatable, dimension(:,:,:) :: tmp
#endif


type(ocean_domain_type), pointer :: Dom =>NULL()
type(ocean_grid_type), pointer   :: Grd =>NULL()

character(len=128) :: version = &
     '$Id$'
character (len=128) :: tagname = &
     '$Name$'

! for diagnostics 
integer :: id_pgrad(2)  =-1
integer :: id_h_pgrad(2)=-1
logical :: used

logical :: debug_pressure=.false. 
logical :: module_is_initialized=.false.
logical :: have_obc=.false. 

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_pressure_init">
!
! <DESCRIPTION>
! Initialize the pressure module
! </DESCRIPTION>
!
subroutine ocean_pressure_init(Grid, Domain, Time, obc)

  type(ocean_grid_type), target, intent(in)   :: Grid
  type(ocean_domain_type), target, intent(in) :: Domain
  type(ocean_time_type), intent(in)           :: Time
  logical, intent(in)                         :: obc

  integer :: ioun, io_status, ierr
  namelist /ocean_pressure_nml/ debug_pressure 

  if ( module_is_initialized ) then 
    call mpp_error( FATAL, '==>Error: ocean_pressure_mod (ocean_pressure_init): module already initialized')
  endif 

  module_is_initialized = .TRUE.

  call write_version_number( version, tagname )

  have_obc = obc

  ioun = open_namelist_file()
  read  (ioun, ocean_pressure_nml,iostat=io_status)
  write (stdout(),'(/)')
  write (stdout(), ocean_pressure_nml)
  write (stdlog(), ocean_pressure_nml)
  ierr = check_nml_error(io_status, 'ocean_pressure_nml')
  call close_file(ioun)

#ifndef STATIC_MEMORY
  call get_local_indices(Domain,isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk = Grid%nk
  allocate(tmp(isd:ied,jsd:jed,2))
#endif

  Dom => Domain
  Grd => Grid

  ! register fields 
  id_pgrad(1) = register_diag_field ('ocean_model', 'pgrad_u', Grd%vel_axes_uv(1:3), Time%model_time, &
     'i-baroclinic pressure force', 'm/s^2', missing_value=-1e10, range=(/-1e10,1e10/))
  id_pgrad(2) = register_diag_field ('ocean_model', 'pgrad_v', Grd%vel_axes_uv(1:3), Time%model_time, &
     'j-baroclinic pressure force', 'm/s^2', missing_value=-1e10, range=(/-1e10,1e10/))
  id_h_pgrad(1) = register_diag_field ('ocean_model', 'h_pgrad_u', Grd%vel_axes_uv(1:3), Time%model_time, &
     'i-baroclinic pressure force', 'm^2/s^2', missing_value=-1e10, range=(/-1e10,1e10/))
  id_h_pgrad(2) = register_diag_field ('ocean_model', 'h_pgrad_v', Grd%vel_axes_uv(1:3), Time%model_time, &
     'j-baroclinic pressure force', 'm^2/s^2', missing_value=-1e10, range=(/-1e10,1e10/))
 

end subroutine ocean_pressure_init
! </SUBROUTINE> NAME="ocean_pressure_init"


!#######################################################################
! <FUNCTION NAME="pressure_in_dbars">
!
! <DESCRIPTION>
! Compute pressure (dbars) exerted at T cell grid point by weight of water
! column between z=0 and grid point 
!
! ro = density in kg/m^3
!
! psurf = surface pressure in kg/m/sec^2 = hydrostatic pressure 
! at z=0 associated with fluid between z=0 and z=eta_t 
! as well as atmospheric pressure patm. 
!
! </DESCRIPTION>
function pressure_in_dbars (Thickness, rho, psurf)

  type(ocean_thickness_type), intent(in)          :: Thickness
  real, intent(in), dimension(isd:ied,jsd:jed,nk) :: rho   ! density
  real, intent(in), dimension(isd:ied,jsd:jed)    :: psurf ! surface pressure
  real, dimension(isd:ied,jsd:jed,nk)             :: pressure_in_dbars
  integer :: k

  if ( .not. module_is_initialized ) then   
    call mpp_error(FATAL, '==>Error in ocean_pressure_mod (pressure_in_dbars): module must be initialized')
  endif 

  pressure_in_dbars(:,:,:) = hydrostatic_pressure (Thickness, rho(:,:,:))*c2dbars

  ! add pressure from free surface height and loading from atmosphere and/or sea ice
  do k=1,nk
     pressure_in_dbars(:,:,k) = pressure_in_dbars(:,:,k) + psurf(:,:)*c2dbars 
  enddo


end function pressure_in_dbars
! </FUNCTION> NAME="pressure_in_dbars"


!#######################################################################
! <FUNCTION NAME="hydrostatic_pressure">
!
! <DESCRIPTION>
! Hydrostatic pressure [newton/m^2] at T cell grid points.
! Integration here is from z=0 to depth of grid point.  This is 
! the so-called "baroclinic" pressure.  If input density "rho" is an anomoly,
! p will be a hydrostatic pressure anomoly. If "rho" is full density, 
! p will be a full hydrostatic pressure.
! </DESCRIPTION>
!
function hydrostatic_pressure (Thickness, rho)

  type(ocean_thickness_type), intent(in)          :: Thickness
  real, intent(in), dimension(isd:ied,jsd:jed,nk) :: rho
  real    :: p5grav
  integer :: i, j, k, kb
  real, dimension(isd:ied,jsd:jed,nk) :: hydrostatic_pressure

  hydrostatic_pressure(:,:,1) = rho(:,:,1)*grav*Grd%dzw(0)

  p5grav = 0.5*grav
  do k=2,nk
    hydrostatic_pressure(:,:,k) = hydrostatic_pressure(:,:,k-1) &
          + Grd%tmask(:,:,k)*(rho(:,:,k-1)+rho(:,:,k))*p5grav*Thickness%dhwt(:,:,k-1) 
  enddo

end function hydrostatic_pressure
! </FUNCTION> NAME="hydrostatic_pressure"


!#######################################################################
! <FUNCTION NAME="pressure_gradient">
!
! <DESCRIPTION>
! Gradient of hydrostatic pressure excluding the surface and atmospheric 
! pressures (i.e., we are computing here the gradient of the baroclinic 
! pressure).  Account is taken of variable partial cell thickness.
! 1 = dp/dx; 2 = dp/dy
! Thickness weight since this is what we wish to use in update of 
! the velocity. 
! </DESCRIPTION>
!
function pressure_gradient(Time, Velocity, Thickness, rho, diag_flag)

  type(ocean_time_type), intent(in)               :: Time
  type(ocean_velocity_type), intent(in)           :: Velocity
  type(ocean_thickness_type), intent(in)          :: Thickness
  real, intent(in), dimension(isd:ied,jsd:jed,nk) :: rho
  logical, intent(in), optional                   :: diag_flag 
  logical                                         :: send_diagnostics

  real, dimension(isc:iec,jsc:jec,nk,2) :: pressure_gradient
  real, dimension(isd:ied,jsd:jed)      :: diff_x, diff_y

  integer :: i, j, k, tau

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, '==>Error in ocean_pressure_mod (pressure_gradient): module must be initialized')
  endif 

  tau = Time%tau

  ! assume send_diagnostics=.true., unless diag_flag says it is false.
  if (PRESENT(diag_flag)) then 
    send_diagnostics = diag_flag
  else 
    send_diagnostics = .true.
  endif 

  wrk2(:,:,:) = rho(:,:,:)-rho0                               ! subtract rho0 to improve accuracy 
  wrk1(:,:,:) = hydrostatic_pressure(Thickness, wrk2(:,:,:))  ! hydrostatic pressure anomaly (Pascals)
  diff_x(:,:) = 0.0
  diff_y(:,:) = 0.0
  do k=1,nk
    do j=jsd,jed
      do i=isd,iec
        diff_x(i,j) = (Thickness%ztp(i+1,j,k)-Thickness%ztp(i,j,k))
      enddo
    enddo
    do j=jsd,jec
      do i=isd,ied
        diff_y(i,j) = (Thickness%ztp(i,j+1,k)-Thickness%ztp(i,j,k))
      enddo
    enddo
    tmp(:,:,1) = rho0r*Grd%umask(:,:,k)*( FDX_NT(FAY(wrk1(:,:,k))) - grav*FAY(FAX(wrk2(:,:,k))*diff_x(:,:))*Grd%dxur(:,:) )
    tmp(:,:,2) = rho0r*Grd%umask(:,:,k)*( FDY_ET(FAX(wrk1(:,:,k))) - grav*FAX(FAY(wrk2(:,:,k))*diff_y(:,:))*Grd%dyur(:,:) )

    pressure_gradient(isc:iec,jsc:jec,k,:) = tmp(isc:iec,jsc:jec,:)

  enddo
  
  if(send_diagnostics) then 
     if (id_pgrad(1) > 0) used = send_data( id_pgrad(1), pressure_gradient(isc:iec,jsc:jec,:,1), &
                          Time%model_time, rmask=Grd%umask(isc:iec,jsc:jec,:))
     if (id_pgrad(2) > 0) used = send_data( id_pgrad(2), pressure_gradient(isc:iec,jsc:jec,:,2), &
                          Time%model_time, rmask=Grd%umask(isc:iec,jsc:jec,:))
  endif 

  ! thickness weight
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           pressure_gradient(i,j,k,1) = pressure_gradient(i,j,k,1)*Thickness%dhu(i,j,k,tau)
           pressure_gradient(i,j,k,2) = pressure_gradient(i,j,k,2)*Thickness%dhu(i,j,k,tau)
        enddo
     enddo
  enddo
  if(send_diagnostics) then 
     if (id_h_pgrad(1) > 0) used = send_data( id_h_pgrad(1), pressure_gradient(isc:iec,jsc:jec,:,1), &
                                   Time%model_time, rmask=Grd%umask(isc:iec,jsc:jec,:))
     if (id_h_pgrad(2) > 0) used = send_data( id_h_pgrad(2), pressure_gradient(isc:iec,jsc:jec,:,2), &
                                   Time%model_time, rmask=Grd%umask(isc:iec,jsc:jec,:))
  endif 


  if (have_obc) call store_ocean_obc_pressure_grad(Thickness, pressure_gradient,tau) 


end function pressure_gradient
! </FUNCTION> NAME="pressure_gradient"

end module ocean_pressure_mod
