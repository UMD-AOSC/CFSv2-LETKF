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
module ocean_thickness_mod
!
!<CONTACT EMAIL="Stephen.Griffies@noaa.gov"> S.M. Griffies 
!</CONTACT>
!
!<OVERVIEW>
! Determine the thickness of grid cells.  
!</OVERVIEW>
!
!<DESCRIPTION>
! This module determines the thickness of grid cells. 
! Thicknesses are generally time dependent.  
!</DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! S.M. Griffies, Developments toward mom4p1 (2004)
! </REFERENCE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_thickness_nml">
!  <DATA NAME="debug_thickness" TYPE="logical">
!  For debugging.
!  </DATA> 
! </NAMELIST>
!
use constants_mod,     only: epsln
use diag_manager_mod,  only: register_diag_field, send_data
use fms_mod,           only: write_version_number, error_mesg, FATAL
use fms_mod,           only: open_namelist_file, close_file, check_nml_error
use mpp_domains_mod,   only: mpp_update_domains, mpp_global_field, domain2d
use mpp_mod,           only: stdout, stdlog, mpp_error

use ocean_domains_mod, only: get_local_indices
use ocean_grids_mod,   only: update_boundaries 
use ocean_types_mod,   only: ocean_time_type, ocean_domain_type, ocean_external_mode_type
use ocean_types_mod,   only: ocean_grid_type, ocean_thickness_type 

implicit none

private

#include <ocean_memory.h>

! for diagnostics 
integer :: id_dht=-1
integer :: id_dhu=-1
logical :: used

! for debugging 
logical :: debug_thickness=.false.

character(len=128) :: version=&
     '$Id$'
character (len=128) :: tagname = &
     '$Name$'

public ocean_thickness_init
public update_tcell_thickness
public update_ucell_thickness

logical :: module_is_initialized = .FALSE.

namelist /ocean_thickness_nml/ debug_thickness

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_thickness_init">
!
! <DESCRIPTION>
! Initialize the thickness and the depth of partial cell point. 
! </DESCRIPTION>
!
subroutine ocean_thickness_init (Time, Domain, Grid, Ext_mode, Thickness, debug)

  type(ocean_time_type), intent(in)          :: Time
  type(ocean_domain_type), intent(inout)     :: Domain
  type(ocean_grid_type), intent(in)          :: Grid  
  type(ocean_external_mode_type), intent(in) :: Ext_mode
  type(ocean_thickness_type), intent(inout)  :: Thickness
  logical, intent(in), optional              :: debug

  real, allocatable, dimension(:,:) :: tmp
  integer                           :: ioun, io_status, ierr
  integer                           :: i, j, k, n, kb

  if ( module_is_initialized ) then
    call mpp_error(FATAL, '==>Error from ocean_grids_mod (ocean_grids_init): module has been initialized')
  endif 

  module_is_initialized = .TRUE.

  call write_version_number(version, tagname)

  if (PRESENT(debug)) debug_thickness = debug

  ! provide for namelist over-ride of defaults 
  ioun = open_namelist_file()
  read (ioun,ocean_thickness_nml,IOSTAT=io_status)
  write (stdout(),'(/)')
  write (stdout(),ocean_thickness_nml)  
  write (stdlog(),ocean_thickness_nml)
  ierr = check_nml_error(io_status, 'ocean_thickness_nml')
  call close_file(ioun)
 
#ifndef STATIC_MEMORY
  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk = Grid%nk

  allocate (Thickness%dht(isd:ied,jsd:jed,nk,3) )
  allocate (Thickness%dhtr(isd:ied,jsd:jed,nk) )
  allocate (Thickness%dhwt(isd:ied,jsd:jed,0:nk) )
  allocate (Thickness%dhu(isd:ied,jsd:jed,nk,3) )
  allocate (Thickness%dhur(isd:ied,jsd:jed,nk) )
  allocate (Thickness%dhwu(isd:ied,jsd:jed,0:nk) )
  allocate (Thickness%ztp(isd:ied,jsd:jed,nk) )
#endif
 
 Thickness%dhu =0.0
 Thickness%dhwu=0.0

 Thickness%dhwt(:,:,0)  = Grid%dzw(0)
 Thickness%dhwu(:,:,0)  = Grid%dzw(0)
 do k=1,nk
    Thickness%dht(:,:,k,:) = Grid%dzt(k)
    Thickness%ztp(:,:,k)   = Grid%zt(k)
    Thickness%dhwt(:,:,k)  = Grid%dzw(k)
 enddo

 ! modifications for partial cells
 do j=jsd,jed
    do i=isd,ied
       kb=Grid%kmt(i,j)
       if (kb .gt. 1) then
           Thickness%dht(i,j,kb,:)  = Grid%ht(i,j)  - Grid%zw(kb-1)
           Thickness%ztp(i,j,kb)    = Grid%zw(kb-1) + Grid%fracdz(kb,0)*Thickness%dht(i,j,kb,1)
           Thickness%dhwt(i,j,kb-1) = Thickness%ztp(i,j,kb) - Grid%zt(kb-1)
       endif
    enddo
 enddo

  ! construct dhu and dhwu as minimum of surrounding dht and dhwt
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           Thickness%dhu(i,j,k,:) = min(Thickness%dht(i,j,k,:),   Thickness%dht(i+1,j,k,:), &
                                        Thickness%dht(i,j+1,k,:), Thickness%dht(i+1,j+1,k,:))
           Thickness%dhwu(i,j,k)  = min(Thickness%dhwt(i,j,k),    Thickness%dhwt(i+1,j,k),  &
                                        Thickness%dhwt(i,j+1,k),  Thickness%dhwt(i+1,j+1,k))
        enddo
     enddo
  enddo

  allocate(tmp(Grid%ni,Grid%nj))
  do n=1,3
     do k=1,nk
        call mpp_global_field(Domain%domain2d,Thickness%dhu(:,:,k,n), tmp)
        call update_boundaries (Domain, Grid, Thickness%dhu(:,:,k,n), tmp, -1, -1) 
     enddo
  enddo
  do k=1,nk
     call mpp_global_field(Domain%domain2d,Thickness%dhwu(:,:,k), tmp)
     call update_boundaries (Domain, Grid, Thickness%dhwu(:,:,k), tmp, -1, -1)
  enddo
  deallocate(tmp)


  ! initialize top cell
  do n=1,3
    Thickness%dht(:,:,1,n) = Grid%dzt(1) + Ext_mode%eta_t(:,:,n)
    Thickness%dhu(:,:,1,n) = Grid%dzt(1) + Ext_mode%eta_u(:,:,n)
  enddo

  ! compute inverse thicknesses (m^-1)
  do k=1,nk
     do j=jsd,jed
        do i=isd,ied
           Thickness%dhtr(i,j,k)  = 1.0/(Thickness%dht(i,j,k,1)+epsln)
           Thickness%dhur(i,j,k)  = 1.0/(Thickness%dhu(i,j,k,1)+epsln)
        enddo
     enddo
  enddo

  id_dhu = register_diag_field ('ocean_model', 'dhu', Grid%tracer_axes(1:3), Time%model_time, &
                                'static ocean u-cell thickness', 'm',&
                                 missing_value=-1.e10, range=(/-1.e10,1.e10/))
  id_dht = register_diag_field ('ocean_model', 'dht', Grid%tracer_axes(1:3), Time%model_time, &
                                'static ocean t-cell thickness', 'm',&
                                 missing_value=-1.e10, range=(/-1.e10,1.e10/))

end subroutine ocean_thickness_init
! </SUBROUTINE> NAME="ocean_thickness_init"



!#######################################################################
! <SUBROUTINE NAME="update_tcell_thickness">
!
! <DESCRIPTION>
! Update time dependent thickness of T cells.
! </DESCRIPTION>
!
subroutine update_tcell_thickness (Time, Ext_mode, Grid, Thickness)

  type(ocean_time_type), intent(in)          :: Time
  type(ocean_external_mode_type), intent(in) :: Ext_mode
  type(ocean_grid_type), intent(in)          :: Grid  
  type(ocean_thickness_type), intent(inout)  :: Thickness

  integer :: i, j
  integer :: tau, taup1

  tau   = Time%tau
  taup1 = Time%taup1

  do j=jsd,jed
     do i=isd,ied
        Thickness%dht(i,j,1,taup1) = Grid%dzt(1) + Ext_mode%eta_t(i,j,taup1)
        Thickness%dhwt(i,j,0)      = Grid%dzw(0) + Ext_mode%eta_t(i,j,taup1)
        Thickness%dhtr(i,j,1)      = 1.0/(Thickness%dht(i,j,1,taup1)+epsln)
     enddo
  enddo

  if (id_dht > 0) used = send_data (id_dht, Thickness%dht(:,:,:,tau), &
                         Time%model_time, rmask=Grid%tmask(:,:,:), &
                         is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

end subroutine update_tcell_thickness
! </SUBROUTINE> NAME="update_tcell_thickness"


!#######################################################################
! <SUBROUTINE NAME="update_ucell_thickness">
!
! <DESCRIPTION>
! Update time dependent thickness of U cells.
! </DESCRIPTION>
!
subroutine update_ucell_thickness (Time, Ext_mode, Grid, Thickness)

  type(ocean_time_type), intent(in)          :: Time
  type(ocean_external_mode_type), intent(in) :: Ext_mode
  type(ocean_grid_type), intent(in)          :: Grid  
  type(ocean_thickness_type), intent(inout)  :: Thickness

  integer :: i, j
  integer :: tau, taup1

  tau   = Time%tau
  taup1 = Time%taup1

  do j=jsd,jed
     do i=isd,ied
       Thickness%dhu(i,j,1,taup1) = Grid%dzt(1) + Ext_mode%eta_u(i,j,taup1)
       Thickness%dhur(i,j,1)      = 1.0/Thickness%dhu(i,j,1,taup1)
       Thickness%dhwu(i,j,0)      = Grid%dzw(0) + Ext_mode%eta_u(i,j,taup1)
     enddo
  enddo

  if (id_dhu > 0) used = send_data (id_dhu, Thickness%dhu(:,:,:,tau), &
                         Time%model_time, rmask=Grid%umask(:,:,:),    &
                         is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

end subroutine update_ucell_thickness
! </SUBROUTINE> NAME="update_ucell_thickness"


end module ocean_thickness_mod
