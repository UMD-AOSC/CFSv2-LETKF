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
module ocean_velocity_advect_mod
!
!<CONTACT EMAIL="Ronald.Pacanowski@noaa.gov"> R.C. Pacanowski 
!</CONTACT>
!
!<CONTACT EMAIL="Stephen.Griffies@noaa.gov"> S.M. Griffies 
!</CONTACT>
!
!<OVERVIEW>
! Velocity advective transport 
!</OVERVIEW>
!
!<DESCRIPTION>
! This module computes advection of velocity using a 
! second order centered scheme. 
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
! S.M. Griffies, M.J. Harrison, A. Rosati, and R.C. Pacanowski 
! A Guide to MOM4 for Users and Developers (2002)
! </REFERENCE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_velocity_advect_nml">
!  <DATA NAME="debug_velocity_advect" TYPE="logical">
!  For debugging
!  </DATA> 
!</NAMELIST>

use diag_manager_mod, only: register_diag_field, send_data
use fms_mod,          only: mpp_error, FATAL, NOTE, stdout, stdlog, write_version_number
use fms_mod,          only: read_data, open_namelist_file, check_nml_error, close_file
use mpp_domains_mod,  only: mpp_update_domains, BGRID_NE
use mpp_mod,          only: mpp_chksum, mpp_sync

use ocean_domains_mod,   only: get_local_indices
use ocean_obc_mod,       only: ocean_obc_update_boundary
use ocean_operators_mod, only: FAX, FAY, BDX_EU, BDY_NU, REMAP_BT_TO_BU
use ocean_types_mod,     only: ocean_grid_type, ocean_domain_type
use ocean_types_mod,     only: ocean_velocity_type, ocean_adv_vel_type
use ocean_types_mod,     only:  ocean_time_type, ocean_thickness_type

implicit none

private

public horz_advection_of_velocity
public vertical_advection_of_velocity
public ocean_velocity_advect_init

! for diagnostics 
logical :: used 
integer :: id_hadv_u=-1
integer :: id_hadv_v=-1
integer :: id_vadv_u=-1
integer :: id_vadv_v=-1
integer :: id_surf_accel(2)=-1

#include <ocean_memory.h>

#ifdef STATIC_MEMORY
real, dimension(isd:ied,jsd:jed)  :: tmp, tmp1, tmp2
#else
real, dimension(:,:), allocatable :: tmp, tmp1, tmp2
#endif

type(ocean_grid_type), pointer   :: Grd =>NULL()
type(ocean_domain_type), pointer :: Dom =>NULL()

character(len=128) :: version=&
  '$Id$'
character (len=128) :: tagname = &
     '$Name$'

logical :: module_is_initialized  = .FALSE.
logical :: debug_velocity_advect  = .false.
logical :: have_obc               = .false.

namelist /ocean_velocity_advect_nml/ debug_velocity_advect

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_velocity_advect_init">
!
! <DESCRIPTION>
! Initialize the velocity advection module 
! </DESCRIPTION>
!
subroutine ocean_velocity_advect_init(Grid, Domain, Time, obc, debug)

  type(ocean_grid_type), intent(in), target   :: Grid
  type(ocean_domain_type), intent(in), target :: Domain
  type(ocean_time_type), intent(in), target   :: Time
  logical, intent(in) :: obc
  logical, intent(in), optional :: debug

  integer :: ioun, io_status, ierr

  if ( module_is_initialized ) then 
    call mpp_error(FATAL,'==>Error from ocean_velocity_advect_mod (ocean_velocity_advect_init): module already initialized')
  endif 

  module_is_initialized = .TRUE.

  call write_version_number( version, tagname )

  Grd => Grid
  Dom => Domain

  have_obc       = obc

  ! provide for namelist over-ride
  ioun = open_namelist_file()
  read  (ioun, ocean_velocity_advect_nml,iostat=io_status)
  write (stdout(),'(/)')
  write (stdout(), ocean_velocity_advect_nml)  
  write (stdlog(), ocean_velocity_advect_nml)
  ierr = check_nml_error(io_status,'ocean_velocity_advect_nml')
  call close_file (ioun)

  if (PRESENT(debug)) debug_velocity_advect = debug
  
#ifndef STATIC_MEMORY
  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk=Grid%nk

  allocate(tmp(isd:ied,jsd:jed),tmp1(isd:ied,jsd:jed),tmp2(isd:ied,jsd:jed))
#endif


  id_hadv_u = register_diag_field ('ocean_model', 'hadv_u', Grd%vel_axes_uv(1:3), Time%model_time, &
     'Thickness wghtd horz advection of u', 'm^2/s^2', missing_value=-1e10, range=(/-1e10,1e10/))
  id_hadv_v = register_diag_field ('ocean_model', 'hadv_v', Grd%vel_axes_uv(1:3), Time%model_time, &
     'Thickness wghtd horz advection of v', 'm^2/s^2', missing_value=-1e10, range=(/-1e10,1e10/))
  id_vadv_u = register_diag_field ('ocean_model', 'vadv_u', Grd%vel_axes_uv(1:3), Time%model_time, &
     'Thickness wghtd vert advection of u', 'm^2/s^2', missing_value=-1e10, range=(/-1e10,1e10/))
  id_vadv_v = register_diag_field ('ocean_model', 'vadv_v', Grd%vel_axes_uv(1:3), Time%model_time, &
     'Thickness wghtd vert advection of v', 'm^2/s^2', missing_value=-1e10, range=(/-1e10,1e10/))
  id_surf_accel(1) = register_diag_field ('ocean_model', 'surf_accel_u', Grd%vel_axes_uv(1:2), &
     Time%model_time, 'uh-forcing from fresh water', 'm^2/s^2', missing_value=-1e10, range=(/-1e10,1e10/))
  id_surf_accel(2) = register_diag_field ('ocean_model', 'surf_accel_v', Grd%vel_axes_uv(1:2), &
     Time%model_time, 'vh-forcing from fresh water', 'm^2/s^2', missing_value=-1e10, range=(/-1e10,1e10/))


end subroutine ocean_velocity_advect_init
! </SUBROUTINE> NAME="ocean_velocity_advect_init"


!#######################################################################
! <FUNCTION NAME="horz_advection_of_velocity">
!
! <DESCRIPTION>
! Compute thickness weighted acceleration m^2/s^2 due to horzizontal 
! advection of velocity.
! </DESCRIPTION>
!
function horz_advection_of_velocity(Time, Thickness, Velocity, Adv_vel, diag_flag)

  type(ocean_time_type), intent(in)      :: Time
  type(ocean_thickness_type), intent(in) :: Thickness
  type(ocean_velocity_type), intent(in)  :: Velocity
  type(ocean_adv_vel_type), intent(in)   :: Adv_vel
  logical, intent(in), optional          :: diag_flag 
  logical                                :: send_diagnostics
 
  real, dimension(isc:iec,jsc:jec,nk,2)  :: horz_advection_of_velocity
  real                                   :: temp1, temp2
  integer                                :: i, j, k, n, tau

  if ( .not. module_is_initialized ) then 
     call mpp_error(FATAL,&
      '==>Error from ocean_velocity_advect_mod (horz_advection_of_velocity): module not yet initialized')
  endif 

  ! assume send_diagnostics=.true., unless diag_flag says it is false.
  if (PRESENT(diag_flag)) then 
    send_diagnostics = diag_flag
  else 
    send_diagnostics = .true.
  endif 

  tau = Time%tau

  do k=1,nk
     do n=1,2

        tmp1(:,:) = Adv_vel%uh_eu(:,:,k)*FAX(Velocity%u(:,:,k,n,tau))
        tmp2(:,:) = Adv_vel%vh_nu(:,:,k)*FAY(Velocity%u(:,:,k,n,tau))

        if(have_obc) then
          call ocean_obc_update_boundary(tmp1(:,:), 'M','s')
          call ocean_obc_update_boundary(tmp2(:,:), 'Z','s')
        endif

        tmp(:,:) =  BDX_EU(tmp1) + BDY_NU(tmp2) 
        do j=jsc,jec
          do i=isc,iec
            tmp(i,j) = tmp(i,j) + (3-2*n)*Thickness%dhu(i,j,k,tau) &
             *Velocity%u(i,j,k,3-n,tau)*(Grd%dh1dy(i,j)*Velocity%u(i,j,k,1,tau)-Grd%dh2dx(i,j)*Velocity%u(i,j,k,2, tau))
            horz_advection_of_velocity(i,j,k,n) = Grd%umask(i,j,k)*tmp(i,j)
          enddo
        enddo

     enddo
  enddo

  if(debug_velocity_advect) then 
      write(stdout(),*) 'horz_advection_of_velocity(1) ==> ',&
                         mpp_chksum(horz_advection_of_velocity(isc:iec,jsc:jec,:,1))
      write(stdout(),*) 'horz_advection_of_velocity(2) ==> ',&
                         mpp_chksum(horz_advection_of_velocity(isc:iec,jsc:jec,:,2))
  endif

  if(send_diagnostics) then 
     if (id_hadv_u > 0) used = send_data( id_hadv_u, horz_advection_of_velocity(isc:iec,jsc:jec,:,1), &
                               Time%model_time, rmask=Grd%umask(isc:iec,jsc:jec,:))
     if (id_hadv_v > 0) used = send_data( id_hadv_v, horz_advection_of_velocity(isc:iec,jsc:jec,:,2), &
                               Time%model_time, rmask=Grd%umask(isc:iec,jsc:jec,:))
  endif 

end function horz_advection_of_velocity
! </FUNCTION> NAME="horz_advection_of_velocity"


!#######################################################################
! <FUNCTION NAME="vertical_advection_of_velocity">
!
! <DESCRIPTION>
! Compute thickness weighted acceleration (m^2/s^2) due to vertical 
! advection of velocity.  Include vertical advection due to fresh
! water entering the surface cells.  
! </DESCRIPTION>
! 
function vertical_advection_of_velocity(Time, Velocity, Adv_vel, pme, river, upme, uriver, diag_flag)

  type(ocean_time_type), intent(in)               :: Time
  type(ocean_velocity_type), intent(in)           :: Velocity
  type(ocean_adv_vel_type), intent(in)            :: Adv_vel
  real, intent(in), dimension(isd:ied,jsd:jed)    :: pme
  real, intent(in), dimension(isd:ied,jsd:jed)    :: river
  real, intent(in), dimension(isd:ied,jsd:jed,2)  :: upme
  real, intent(in), dimension(isd:ied,jsd:jed,2)  :: uriver 
  logical, intent(in), optional                   :: diag_flag 
  logical                                         :: send_diagnostics

  real,dimension(isc:iec,jsc:jec,nk,2) :: vertical_advection_of_velocity 

  real, dimension(isd:ied,jsd:jed)   :: river_u
  real, dimension(isd:ied,jsd:jed)   :: pme_u 
  real,dimension(isd:ied,jsd:jed)    :: ft1, ft2
  real, dimension(isd:ied,jsd:jed,2) :: surf_accel 
  real,dimension(nk)                 :: kmask
  integer                            :: i, j, k, kp1, n, tau

  if ( .not. module_is_initialized ) then 
     call mpp_error(FATAL,&
      '==>Error from ocean_velocity_advect_mod (vertical_advection_of_velocity): module not yet initialized')
  endif 

  ! assume send_diagnostics=.true., unless diag_flag says it is false.
  if (PRESENT(diag_flag)) then 
    send_diagnostics = diag_flag
  else 
    send_diagnostics = .true.
  endif 

  tau           = Time%tau
  kmask(1:nk-1) = 1.0
  kmask(nk)     = 0.0

  pme_u     = 0.0
  river_u   = 0.0
  surf_accel = 0.0

  pme_u   =  REMAP_BT_TO_BU(pme(:,:))
  river_u =  REMAP_BT_TO_BU(river(:,:))
  if(id_surf_accel(1) > 0 .or. id_surf_accel(2) > 0) then   
      do n=1,2
         do j=jsc,jec
            do i=isc,iec
               surf_accel(i,j,n) = Grd%umask(i,j,1)*(pme_u(i,j)*upme(i,j,n) + river_u(i,j)*uriver(i,j,n))
            enddo
         enddo
      enddo
  endif

  do n=1,2
     do j=jsc,jec
        do i=isc,iec
           ft1(i,j) = -pme_u(i,j)*upme(i,j,n) -river_u(i,j)*uriver(i,j,n)
           ft2(i,j) = 0.0
        enddo
     enddo
     do k=1,nk
        kp1 = min(k+1,nk)
        do j=jsc,jec
           do i=isc,iec
              ft2(i,j) = Adv_vel%w_bu(i,j,k)*0.5*(Velocity%u(i,j,k,n,tau) + kmask(k)*Velocity%u(i,j,kp1,n,tau))
              vertical_advection_of_velocity(i,j,k,n) = Grd%umask(i,j,k)*(ft1(i,j)-ft2(i,j))
              ft1(i,j) = ft2(i,j)
           enddo
        enddo
     enddo
  enddo

  if(debug_velocity_advect) then 
      write(stdout(),*) 'vertical_advection_of_velocity(1) ==> ',&
                         mpp_chksum(vertical_advection_of_velocity(isc:iec,jsc:jec,:,1))
      write(stdout(),*) 'vertical_advection_of_velocity(2) ==> ',&
                         mpp_chksum(vertical_advection_of_velocity(isc:iec,jsc:jec,:,2))
  endif

  if(send_diagnostics) then 
     if (id_vadv_u > 0) used = send_data( id_vadv_u, vertical_advection_of_velocity(isc:iec,jsc:jec,:,1), &
                               Time%model_time, rmask=Grd%umask(isc:iec,jsc:jec,:))
     if (id_vadv_v > 0) used = send_data( id_vadv_v, vertical_advection_of_velocity(isc:iec,jsc:jec,:,2), &
                               Time%model_time, rmask=Grd%umask(isc:iec,jsc:jec,:))
  endif 

end function vertical_advection_of_velocity
! </FUNCTION> NAME="vertical_advection_of_velocity"


end module ocean_velocity_advect_mod
