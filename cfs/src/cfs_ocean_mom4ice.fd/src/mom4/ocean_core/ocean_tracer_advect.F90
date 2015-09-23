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
module ocean_tracer_advect_mod
!
!<CONTACT EMAIL="Matthew.Harrison@noaa.gov"> Matt Harrison 
!</CONTACT>
!
!<CONTACT EMAIL="Ronald.Pacanowski@noaa.gov"> Ron Pacanowski
!</CONTACT>
!
! <CONTACT EMAIL="Stephen.Griffies@noaa.gov"> S. M. Griffies 
! </CONTACT>
!
! <CONTACT EMAIL="John.Dunne@noaa.gov"> John Dunne 
! </CONTACT>
!
! <CONTACT EMAIL="Alistair.Adcroft@noaa.gov"> Alistair Adcroft 
! </CONTACT>
!
! <CONTACT EMAIL="Jennifer.Simeon@noaa.gov"> J. Simeon 
! </CONTACT>
!
!<OVERVIEW>
! This module computes thickness weighted tracer advection tendencies. 
!</OVERVIEW>
!
!<DESCRIPTION>
! This module computes thickness weighted tracer advection tendencies
! using a variety of advection schemes.  
!</DESCRIPTION>
!
!<NOTE>
! Contains a version of quicker with MOM3 masking.
!</NOTE>
!
!<NAMELIST NAME="ocean_tracer_advect_nml">
!  <DATA NAME="limit_with_upwind" TYPE="logical">
!  If true, will compute tracer fluxes entering a cell using upwind 
!  if the tracer value is outside a specified range. Implemented
!  only for quick at this time. This is an ad hoc and incomplete attempt
!  to maintain monotonicity with the quicker scheme.  
!  </DATA> 
!  <DATA NAME="debug_tracer_advect" TYPE="logical">
!  For debugging 
!  </DATA> 
!</NAMELIST>

use fms_mod,             only: write_version_number, check_nml_error, close_file, open_namelist_file
use mpp_domains_mod,     only: mpp_define_domains, mpp_update_domains
use mpp_domains_mod,     only: BGRID_NE, CGRID_NE, domain2d, XUPDATE, YUPDATE
use mpp_mod,             only: mpp_error, mpp_chksum, FATAL, NOTE, stdout, stdlog
use mpp_mod,             only: mpp_clock_id, mpp_clock_begin, mpp_clock_end, CLOCK_ROUTINE
use diag_manager_mod,    only: register_diag_field, send_data

use ocean_domains_mod,   only: get_local_indices, set_ocean_domain
use ocean_obc_mod,       only: store_ocean_obc_tracer_flux, store_ocean_obc_tracer_advect
use ocean_topog_mod,     only: ocean_topog_init
use ocean_types_mod,     only: ocean_domain_type, ocean_grid_type, ocean_time_type
use ocean_types_mod,     only: ocean_prog_tracer_type, ocean_thickness_type, ocean_adv_vel_type
use ocean_types_mod,     only: ADVECT_UPWIND, ADVECT_2ND_ORDER, ADVECT_4TH_ORDER, ADVECT_6TH_ORDER
use ocean_types_mod,     only: ADVECT_QUICKER, ADVECT_QUICKMOM3, ADVECT_MDFL_SUP_B, ADVECT_MDFL_SWEBY

implicit none

private
integer :: id0, id1, now
integer :: id_clock_up_horz
integer :: id_clock_up_vert
integer :: id_clock_2nd_horz
integer :: id_clock_2nd_vert
integer :: id_clock_4th_horz
integer :: id_clock_4th_vert
integer :: id_clock_6th_horz
integer :: id_clock_6th_vert
integer :: id_clock_quick_horz
integer :: id_clock_quick_vert
integer :: id_clock_quickmom3_horz
integer :: id_clock_quickmom3_vert
integer :: id_clock_mdfl_sup_b
integer :: id_clock_mdfl_sweby
integer :: id_clock_mdfl_sweby_all

#include <ocean_memory.h>

#ifdef STATIC_MEMORY

real, dimension(isd:ied,jsd:jed,nk) :: flux_x
real, dimension(isd:ied,jsd:jed,nk) :: flux_y
real, dimension(isd:ied,jsd:jed,nk) :: flux_z

! some fields for quicker 
real, dimension(isd:ied,jsd:jed,nk)         :: quick_x
real, dimension(isd:ied,jsd:jed,nk)         :: quick_y                    
real, dimension(nk,2)                       :: quick_z                             
real, dimension(isd:ied,jsd:jed,3)          :: curv_xp
real, dimension(isd:ied,jsd:jed,3)          :: curv_xn
real, dimension(isd:ied,jsd:jed,3)          :: curv_yp
real, dimension(isd:ied,jsd:jed,3)          :: curv_yn  
real, dimension(nk,3)                       :: curv_zp                   
real, dimension(nk,3)                       :: curv_zn                    
real, dimension(isc-2:iec+2,jsc-2:jec+2)    :: dxt_quick
real, dimension(isc-2:iec+2,jsc-2:jec+2)    :: dyt_quick
real, dimension(isc-2:iec+2,jsc-2:jec+2,nk) :: tmask_quick
real, dimension(isc-2:iec+2,jsc-2:jec+2,nk) :: tracer_quick

! masks and fields for other advection schemes 
real, dimension(isc-2:iec+2,jsc-2:jec+2,nk) :: tmask_fourth
real, dimension(isc-3:iec+3,jsc-3:jec+3,nk) :: tmask_sixth
real, dimension(isc-2:iec+2,jsc-2:jec+2,nk) :: tmask_mdfl
real, dimension(isc-2:iec+2,jsc-2:jec+2,nk) :: tracer_fourth
real, dimension(isc-3:iec+3,jsc-3:jec+3,nk) :: tracer_sixth
real, dimension(isc-2:iec+2,jsc-2:jec+2,nk) :: tracer_mdfl

! for doing all tracers with sweby 
type  :: tracer_mdfl_type
  real, dimension(isc-2:iec+2,jsc-2:jec+2,nk) :: field
end type tracer_mdfl_type

#else

real, dimension(:,:,:), allocatable :: flux_x
real, dimension(:,:,:), allocatable :: flux_y
real, dimension(:,:,:), allocatable :: flux_z

! some flow independent quantities for quicker 
real, dimension(:,:,:), allocatable :: quick_x, quick_y 
real, dimension(:,:), allocatable   :: quick_z 
real, dimension(:,:,:), allocatable :: curv_xp, curv_xn, curv_yp, curv_yn
real, dimension(:,:), allocatable   :: curv_zp, curv_zn
real, dimension(:,:), allocatable   :: dxt_quick
real, dimension(:,:), allocatable   :: dyt_quick
real, dimension(:,:,:), allocatable :: tmask_quick
real, dimension(:,:,:), allocatable :: tracer_quick

! masks for other advection schemes 
real, dimension(:,:,:), allocatable :: tmask_fourth
real, dimension(:,:,:), allocatable :: tmask_sixth
real, dimension(:,:,:), allocatable :: tmask_mdfl
real, dimension(:,:,:), allocatable :: tracer_fourth
real, dimension(:,:,:), allocatable :: tracer_sixth
real, dimension(:,:,:), allocatable :: tracer_mdfl

! for doing all tracers with sweby 
type  :: tracer_mdfl_type
  real, dimension(:,:,:), pointer   :: field => NULL()
end type tracer_mdfl_type

#endif

integer :: num_prog_tracers = 0

character(len=256) :: version='CVS $Id$'
character(len=256) :: tagname='Tag $Name$'


type(ocean_domain_type), pointer :: Dom =>NULL()
type(ocean_grid_type)  , pointer :: Grd =>NULL()
type(ocean_domain_type), save    :: Dom_quicker
type(ocean_domain_type), save    :: Dom_fourth
type(ocean_domain_type), save    :: Dom_sixth
type(ocean_domain_type), save    :: Dom_mdfl
type(ocean_domain_type), save    :: Dom_flux

type(tracer_mdfl_type), dimension(:), allocatable :: tracer_mdfl_all  

real, parameter :: a4=7.0/12.0, b4=-1.0/12.0
real, parameter :: a6=37.0/60.0, b6=-2.0/15.0, c6=1.0/60.0
real, parameter :: onesixth=1.0/6.0

logical :: module_is_initialized=.false.
logical :: limit_with_upwind=.false.
logical :: debug_tracer_advect=.false.
logical :: advect_sweby_all=.false.
logical :: have_obc=.false.

logical :: used
integer, dimension(:), allocatable :: id_sweby_advect
integer, dimension(:), allocatable :: id_vert_advect
integer, dimension(:), allocatable :: id_horz_advect
integer, dimension(:), allocatable :: id_adv_flux_x
integer, dimension(:), allocatable :: id_adv_flux_y
integer, dimension(:), allocatable :: id_adv_flux_z
integer, dimension(:), allocatable :: id_adv_flux_x_int_z
integer, dimension(:), allocatable :: id_adv_flux_y_int_z

public  horz_advect_tracer
private horz_advect_tracer_upwind
private horz_advect_tracer_2nd_order
private horz_advect_tracer_4th_order
private horz_advect_tracer_6th_order
private horz_advect_tracer_quicker
private horz_advect_tracer_quickmom3

public  vert_advect_tracer
private vert_advect_tracer_upwind
private vert_advect_tracer_2nd_order
private vert_advect_tracer_4th_order
private vert_advect_tracer_6th_order
private vert_advect_tracer_quicker
private vert_advect_tracer_quickmom3

private advect_tracer_mdfl_sub_b
private advect_tracer_mdfl_sweby
private advect_tracer_sweby_all

private quicker_init
private mdfl_init
private fourth_sixth_init

public tracer_advection_init

namelist /ocean_tracer_advect_nml/ debug_tracer_advect, limit_with_upwind, advect_sweby_all

contains


!#######################################################################
! <SUBROUTINE NAME="tracer_advection_init">
!
! <DESCRIPTION>
! Initialize the tracer advection module.
! </DESCRIPTION>
!
subroutine tracer_advection_init (Grid, Domain, Time, T_prog, obc)

  type(ocean_grid_type), intent(in), target   :: Grid
  type(ocean_domain_type), intent(in), target :: Domain
  type(ocean_time_type), intent(in)           :: Time
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)
  
  integer :: n, index_temp
  integer :: ioun, io_status, ierr
  logical, intent(in) :: obc
  
  write( stdlog(),'(/a/)') trim(version)

  module_is_initialized = .true.

  have_obc = obc

  call write_version_number( version, tagname )

  ! provide for namelist over-ride of defaults 
  ioun = open_namelist_file()
  read  (ioun, ocean_tracer_advect_nml,iostat=io_status)
  write (stdout(),'(/)')
  write (stdout(), ocean_tracer_advect_nml)  
  write (stdlog(), ocean_tracer_advect_nml)
  ierr = check_nml_error(io_status, 'ocean_tracer_advect_nml')
  call close_file(ioun)

  if(limit_with_upwind) then 
    call mpp_error(NOTE,'==>ocean_tracer_advect_mod: limit_with_upwind reverts quicker to upwind if tracer outside limits')
  endif 

  if(advect_sweby_all) then
    write(stdout(),'(/a)')'==>ocean_tracer_advect_mod: advect_sweby_all=.true. so all tracers advected'
    write(stdout(),'(a)') '   with mdfl_sweby, regardless of the settings in the field table.'
    write(stdout(),'(a/)')'   This method exploits mpp_update_domain capabilities for added efficiency.'
    if(have_obc) then 
       call mpp_error(FATAL,'==>ocean_tracer_advect_mod: advect_sweby_all=.true. has not been implemented with obc. ')
    endif 
  endif 

#ifndef STATIC_MEMORY
  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk  = Grid%nk
  allocate(flux_x(isd:ied,jsd:jed,nk))
  allocate(flux_y(isd:ied,jsd:jed,nk))
  allocate(flux_z(isd:ied,jsd:jed,nk))  
#endif
 
  flux_x = 0.0
  flux_y = 0.0
  flux_z = 0.0

  Dom => Domain
  Grd => Grid

  num_prog_Tracers = size(T_prog(:))
  do n=1,num_prog_tracers
     if (trim(T_prog(n)%name) == 'temp') index_temp = n
  enddo

  call mpp_define_domains( (/1,Grid%ni,1,Grid%nj/), Domain%layout, Dom_flux%domain2d&
             , xflags = Domain%xflags, yflags = Domain%yflags, xhalo=1, yhalo=1,name='flux')


  ! initialize schemes other than 2nd order centered and first order upwind 
  call fourth_sixth_init
  call quicker_init
  call mdfl_init

  ! register diagnostic output 
  allocate (id_sweby_advect(num_prog_tracers))
  allocate (id_horz_advect(num_prog_tracers))
  allocate (id_vert_advect(num_prog_tracers))
  allocate (id_adv_flux_x(num_prog_tracers))
  allocate (id_adv_flux_y(num_prog_tracers))
  allocate (id_adv_flux_z(num_prog_tracers))
  allocate (id_adv_flux_x_int_z(num_prog_tracers))
  allocate (id_adv_flux_y_int_z(num_prog_tracers))
  id_sweby_advect=-1
  id_horz_advect=-1
  id_vert_advect=-1
  id_adv_flux_x =-1
  id_adv_flux_y =-1
  id_adv_flux_z =-1
  id_adv_flux_x_int_z=-1
  id_adv_flux_y_int_z=-1

  do n=1,num_prog_tracers 

    if(n==index_temp) then 
      id_sweby_advect(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_sweby_advect', &
                   Grd%tracer_axes_flux_x(1:3), Time%model_time, 'rho_cp*dht*sweby advect tendency', &
                   'Watts/m^2', missing_value=-1.e20, range=(/-1.e20,1.e20/))
      id_horz_advect(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_horz_advect', &
                   Grd%tracer_axes_flux_x(1:3), Time%model_time, 'rho_cp*dht*horz advect tendency', &
                   'Watts/m^2', missing_value=-1.e20, range=(/-1.e20,1.e20/))
      id_vert_advect(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_vert_advect', &
                   Grd%tracer_axes_flux_x(1:3), Time%model_time, 'rho_cp*dht*vert advect tendency', &
                   'Watts/m^2', missing_value=-1.e20, range=(/-1.e20,1.e20/))
      id_adv_flux_x(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_xflux_adv', &
                   Grd%tracer_axes_flux_x(1:3), Time%model_time, 'rho_cp*u*dyt*dht*temp', &
                   'Watts', missing_value=-1.e20, range=(/-1.e20,1.e20/))
      id_adv_flux_y(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_yflux_adv', &
                   Grd%tracer_axes_flux_y(1:3), Time%model_time, 'rho_cp*v*dxt*dht*temp', &
                   'Watts', missing_value=-1.e20, range=(/-1.e20,1.e20/))
      id_adv_flux_z(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_zflux_adv', &
                   Grd%tracer_axes_wt(1:3), Time%model_time, 'rho_cp*wt*dxt*dyt*temp', &
                   'Watts', missing_value=-1.e20, range=(/-1.e20,1.e20/))      
      id_adv_flux_x_int_z(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_adv_flux_x_int_z', &
                   Grd%tracer_axes_flux_x(1:2), Time%model_time, 'z-integral of rho_cp*u*dyt*temp', &
                   'Watts', missing_value=-1.e20, range=(/-1.e20,1.e20/))
      id_adv_flux_y_int_z(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_adv_flux_y_int_z', &
                   Grd%tracer_axes_flux_y(1:2), Time%model_time, 'z-integral of rho_cp*v*dxt*temp', &
                   'Watts', missing_value=-1.e20, range=(/-1.e20,1.e20/))
    else
      id_sweby_advect(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_sweby_advect', &
                   Grd%tracer_axes_flux_x(1:3), Time%model_time, 'rho0*dht*sweby advect tendency', &
                   'm*kg/s)', missing_value=-1.e20, range=(/-1.e20,1.e20/))
      id_horz_advect(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_horz_advect', &
                   Grd%tracer_axes_flux_x(1:3), Time%model_time, 'rho0*dht*horz advect tendency', &
                   'm*kg/s', missing_value=-1.e20, range=(/-1.e20,1.e20/))
      id_vert_advect(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_vert_advect', &
                   Grd%tracer_axes_flux_x(1:3), Time%model_time, 'rho0*dht*vert advect tendency', &
                   'm*kg/s', missing_value=-1.e20, range=(/-1.e20,1.e20/))
      id_adv_flux_x(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_xflux_adv', &
                   Grd%tracer_axes_flux_x(1:3), Time%model_time, 'rho0*u*dyt*dht*tracer',&
                   &' kg/sec', missing_value=-1.e20, range=(/-1.e20,1.e20/))
      id_adv_flux_y(n)   = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_yflux_adv', &
                   Grd%tracer_axes_flux_y(1:3), Time%model_time, 'rho0*v*dyt*dht*tracer',&
                    &' kg/sec', missing_value=-1.e20, range=(/-1.e20,1.e20/))
      id_adv_flux_z(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_zflux_adv', &
                   Grd%tracer_axes_wt(1:3), Time%model_time, 'rho0*wt*dxt*dyt*temp', &
                   'kg/sec', missing_value=-1.e20, range=(/-1.e20,1.e20/))            
      id_adv_flux_x_int_z(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_adv_flux_x_int_z', &
                   Grd%tracer_axes_flux_x(1:2), Time%model_time, 'z-integral of rho0*u*dyt*dht*tracer',&
                   &' kg/sec', missing_value=-1.e20, range=(/-1.e20,1.e20/))
      id_adv_flux_y_int_z(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_adv_flux_y_int_z', &
                   Grd%tracer_axes_flux_y(1:2), Time%model_time, 'z-integral of rho0*v*dxt*dht*tracer',&
                   &' kg/sec', missing_value=-1.e20, range=(/-1.e20,1.e20/))
    endif 

  enddo 

  ! initialize clock ids 
  id_clock_up_horz        = mpp_clock_id('(Ocean advect: horz up)      ',grain=CLOCK_ROUTINE)
  id_clock_up_vert        = mpp_clock_id('(Ocean advect: vert up)      ',grain=CLOCK_ROUTINE)
  id_clock_2nd_horz       = mpp_clock_id('(Ocean advect: horz 2nd)     ',grain=CLOCK_ROUTINE)
  id_clock_2nd_vert       = mpp_clock_id('(Ocean advect: vert 2nd)     ',grain=CLOCK_ROUTINE)

end subroutine tracer_advection_init
! </SUBROUTINE> NAME="tracer_advection_init"



!#######################################################################
! <SUBROUTINE NAME="quicker_init">
!
! <DESCRIPTION>
! Initialize quicker specific fields. 
! </DESCRIPTION>
!
subroutine quicker_init

  integer :: i, j, k, kp1, km1, kp2

#ifndef STATIC_MEMORY
  allocate(quick_x(isd:ied,jsd:jed,2))
  allocate(quick_y(isd:ied,jsd:jed, 2))
  allocate(quick_z(nk,2))
  allocate(curv_xp(isd:ied,jsd:jed,3))
  allocate(curv_xn(isd:ied,jsd:jed, 3))
  allocate(curv_yp(isd:ied,jsd:jed,3))
  allocate(curv_yn(isd:ied,jsd:jed, 3))
  allocate(curv_zp(nk,3))
  allocate(curv_zn(nk,3))
  allocate(dxt_quick(isc-2:iec+2,jsc-2:jec+2))
  allocate(dyt_quick(isc-2:iec+2,jsc-2:jec+2))
  allocate(tmask_quick(isc-2:iec+2,jsc-2:jec+2,nk))    
  allocate(tracer_quick(isc-2:iec+2,jsc-2:jec+2,nk))
#endif
 
  quick_x       = 0.0
  quick_y       = 0.0
  quick_z       = 0.0
  curv_xp       = 0.0
  curv_xn       = 0.0
  curv_yp       = 0.0
  curv_yn       = 0.0
  curv_zp       = 0.0
  curv_zn       = 0.0
  dxt_quick     = 0.0
  dyt_quick     = 0.0
  tmask_quick   = 0.0
  tracer_quick  = 0.0

  call set_ocean_domain(Dom_quicker,Grd,xhalo=2,yhalo=2,name='quicker')

  do i=isc,iec
     do j=jsc,jec
        dxt_quick(i,j) = Grd%dxt(i,j)
        dyt_quick(i,j) = Grd%dyt(i,j)
        do k=1,nk
           tmask_quick(i,j,k)  = Grd%tmask(i,j,k)
        enddo
     enddo
  enddo
  call mpp_update_domains(tmask_quick,Dom_quicker%domain2d)  


  ! fill in boundaries (needed for non-cyclic case to prevent blow-ups)
  do i=isc-2,isc-1
     dxt_quick(i,jsc:jec) = dxt_quick(isc,jsc:jec)
     dyt_quick(i,:) = dyt_quick(isc,:)
  enddo

  do i=iec+1,iec+2
     dxt_quick(i,jsc:jec) = dxt_quick(iec,jsc:jec)
     dyt_quick(i,:) = dyt_quick(iec,:)
  enddo

  do j=jsc-2,jsc-1
     dxt_quick(:,j) = dxt_quick(:,jsc)
     dyt_quick(:,j) = dyt_quick(:,jsc)
  enddo

  do j=jec+1,jec+2
     dxt_quick(:,j) = dxt_quick(:,jec)
     dyt_quick(:,j) = dyt_quick(:,jec)
  enddo
  call mpp_update_domains(dxt_quick,Dom_quicker%domain2d)
  call mpp_update_domains(dyt_quick,Dom_quicker%domain2d)

  
! calculate quicker weights on computational domain (not filled in halos)
  
  do i=isc-1,iec
     do j=jsc-1,jec
        quick_x(i,j,1) = dxt_quick(i+1,j)/(dxt_quick(i+1,j)+dxt_quick(i,j))
        quick_x(i,j,2) = dxt_quick(i,j)/(dxt_quick(i+1,j)+dxt_quick(i,j))
        quick_y(i,j,1) = dyt_quick(i,j+1)/(dyt_quick(i,j+1)+dyt_quick(i,j))
        quick_y(i,j,2) = dyt_quick(i,j)/(dyt_quick(i,j+1)+dyt_quick(i,j))
        
        curv_xp(i,j,1) = dxt_quick(i,j)*dxt_quick(i+1,j)/&
             ((dxt_quick(i-1,j)+2.0*dxt_quick(i,j)+dxt_quick(i+1,j))&
             *(dxt_quick(i,j)+dxt_quick(i+1,j)))
        curv_xp(i,j,2) = -dxt_quick(i,j)*dxt_quick(i+1,j)/&
             ((dxt_quick(i,j)+dxt_quick(i+1,j))&
             *(dxt_quick(i-1,j)+dxt_quick(i,j)))        
        curv_xp(i,j,3) = dxt_quick(i,j)*dxt_quick(i+1,j)/&
             ((dxt_quick(i-1,j)+2.0*dxt_quick(i,j)+dxt_quick(i+1,j))&
             *(dxt_quick(i-1,j)+dxt_quick(i,j)))

        curv_xn(i,j,1) = dxt_quick(i,j)*dxt_quick(i+1,j)/&
             ((dxt_quick(i,j)+2.0*dxt_quick(i+1,j)+dxt_quick(i+2,j))&
             *(dxt_quick(i+1,j)+dxt_quick(i+2,j)))
        curv_xn(i,j,2) = -dxt_quick(i,j)*dxt_quick(i+1,j)/&
             ((dxt_quick(i+1,j)+dxt_quick(i+2,j))&
             *(dxt_quick(i,j)+dxt_quick(i+1,j)))        
        curv_xn(i,j,3) = dxt_quick(i,j)*dxt_quick(i+1,j)/&
             ((dxt_quick(i,j)+2.0*dxt_quick(i+1,j)+dxt_quick(i+2,j))&
             *(dxt_quick(i,j)+dxt_quick(i+1,j)))

        curv_yp(i,j,1) = dyt_quick(i,j)*dyt_quick(i,j+1)/&
             ((dyt_quick(i,j-1)+2.0*dyt_quick(i,j)+dyt_quick(i,j+1))&
             *(dyt_quick(i,j)+dyt_quick(i,j+1)))
        curv_yp(i,j,2) = -dyt_quick(i,j)*dyt_quick(i,j+1)/&
             ((dyt_quick(i,j)+dyt_quick(i,j+1))&
             *(dyt_quick(i,j-1)+dyt_quick(i,j)))        
        curv_yp(i,j,3) = dyt_quick(i,j)*dyt_quick(i,j+1)/&
             ((dyt_quick(i,j-1)+2.0*dyt_quick(i,j)+dyt_quick(i,j+1))&
             *(dyt_quick(i,j-1)+dyt_quick(i,j)))

        curv_yn(i,j,1) = dyt_quick(i,j)*dyt_quick(i,j+1)/&
             ((dyt_quick(i,j)+2.0*dyt_quick(i,j+1)+dyt_quick(i,j+2))&
             *(dyt_quick(i,j+1)+dyt_quick(i,j+2)))
        curv_yn(i,j,2) = -dyt_quick(i,j)*dyt_quick(i,j+1)/&
             ((dyt_quick(i,j+1)+dyt_quick(i,j+2))&
             *(dyt_quick(i,j)+dyt_quick(i,j+1)))        
        curv_yn(i,j,3) = dyt_quick(i,j)*dyt_quick(i,j+1)/&
             ((dyt_quick(i,j)+2.0*dyt_quick(i,j+1)+dyt_quick(i,j+2))&
             *(dyt_quick(i,j)+dyt_quick(i,j+1)))
        
     enddo
  enddo

  do k= 1,nk
     kp2 = min(k+2,nk)
     kp1 = min(k+1,nk)
     km1 = max(k-1,1)
     quick_z(k,1) = Grd%dzt(kp1)/(Grd%dzt(kp1)+Grd%dzt(k))
     quick_z(k,2) = Grd%dzt(k  )/(Grd%dzt(kp1)+Grd%dzt(k))
     curv_zp(k,1) = Grd%dzt(k)*Grd%dzt(kp1)/((Grd%dzt(km1)+2.0*Grd%dzt(k)+Grd%dzt(kp1))*(Grd%dzt(k)+Grd%dzt(kp1)))
     curv_zp(k,2) =-Grd%dzt(k)*Grd%dzt(kp1)/((Grd%dzt(k)+Grd%dzt(kp1))*(Grd%dzt(km1)+Grd%dzt(k)))
     curv_zp(k,3) = Grd%dzt(k)*Grd%dzt(kp1)/((Grd%dzt(km1)+2.0*Grd%dzt(k)+Grd%dzt(kp1))*(Grd%dzt(km1)+Grd%dzt(k)))
     curv_zn(k,1) = Grd%dzt(k)*Grd%dzt(kp1)/((Grd%dzt(k)+2.0*Grd%dzt(kp1)+Grd%dzt(kp2))*(Grd%dzt(kp1)+Grd%dzt(kp2)))
     curv_zn(k,2) =-Grd%dzt(k)*Grd%dzt(kp1)/((Grd%dzt(kp1)+Grd%dzt(kp2))*(Grd%dzt(k)+Grd%dzt(kp1)))
     curv_zn(k,3) = Grd%dzt(k)*Grd%dzt(kp1)/((Grd%dzt(k)+2.0*Grd%dzt(kp1)+Grd%dzt(kp2))*(Grd%dzt(k)+Grd%dzt(kp1)))
  enddo

  ! initialize clock ids 
  id_clock_quick_horz     = mpp_clock_id('(Ocean advect: horz quk)     ',grain=CLOCK_ROUTINE)
  id_clock_quick_vert     = mpp_clock_id('(Ocean advect: vert quk)     ',grain=CLOCK_ROUTINE)
  id_clock_quickmom3_horz = mpp_clock_id('(Ocean advect: horz qukmom3) ',grain=CLOCK_ROUTINE)
  id_clock_quickmom3_vert = mpp_clock_id('(Ocean advect: vert qukmom3) ',grain=CLOCK_ROUTINE)

end subroutine quicker_init
! </SUBROUTINE> NAME="quicker_init"


!#######################################################################
! <SUBROUTINE NAME="fourth_sixth_init">
!
! <DESCRIPTION>
! Initialize the fourth order and sixth order advection fields.
! </DESCRIPTION>
!
subroutine fourth_sixth_init

  integer :: i, j, k

#ifndef STATIC_MEMORY
  allocate(tmask_fourth (isc-2:iec+2,jsc-2:jec+2,nk))    
  allocate(tracer_fourth(isc-2:iec+2,jsc-2:jec+2,nk))    
  allocate(tmask_sixth  (isc-3:iec+3,jsc-3:jec+3,nk))    
  allocate(tracer_sixth (isc-3:iec+3,jsc-3:jec+3,nk))    
#endif
  tmask_fourth  = 0.0
  tracer_fourth = 0.0
  tmask_sixth   = 0.0
  tracer_sixth  = 0.0
  
  call mpp_define_domains( (/1,Grd%ni,1,Grd%nj/), Dom%layout, Dom_fourth%domain2d &
             , xflags = Dom%xflags, yflags = Dom%yflags, xhalo=2, yhalo=2,name='fourth')

  call mpp_define_domains( (/1,Grd%ni,1,Grd%nj/), Dom%layout, Dom_sixth%domain2d &
             , xflags = Dom%xflags, yflags = Dom%yflags, xhalo=3, yhalo=3,name='sixth')

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           tmask_fourth(i,j,k) = Grd%tmask(i,j,k)
           tmask_sixth(i,j,k)  = Grd%tmask(i,j,k)
        enddo
     enddo
  enddo
  call mpp_update_domains(tmask_fourth,Dom_fourth%domain2d)  
  call mpp_update_domains(tmask_sixth ,Dom_sixth%domain2d)  

  ! initialize clock ids 
  id_clock_4th_horz       = mpp_clock_id('(Ocean advect: horz 4th)     ',grain=CLOCK_ROUTINE)
  id_clock_4th_vert       = mpp_clock_id('(Ocean advect: vert 4th)     ',grain=CLOCK_ROUTINE)
  id_clock_6th_horz       = mpp_clock_id('(Ocean advect: horz 6th)     ',grain=CLOCK_ROUTINE)
  id_clock_6th_vert       = mpp_clock_id('(Ocean advect: vert 6th)     ',grain=CLOCK_ROUTINE)


end subroutine fourth_sixth_init
! </SUBROUTINE> NAME="fourth_sixth_init"


!#######################################################################
! <SUBROUTINE NAME="mdfl_init">
!
! <DESCRIPTION>
! Initialize mdfl specific fields.
! </DESCRIPTION>
!
subroutine mdfl_init 

  integer :: i, j, k, n 

  call set_ocean_domain(Dom_mdfl,Grd,xhalo=2,yhalo=2,name='mdfl')

  allocate(tracer_mdfl_all(num_prog_tracers))
#ifndef STATIC_MEMORY
  allocate(tmask_mdfl (isc-2:iec+2,jsc-2:jec+2,nk))    
  allocate(tracer_mdfl(isc-2:iec+2,jsc-2:jec+2,nk))    
  do n=1,num_prog_tracers
     allocate (tracer_mdfl_all(n)%field(isc-2:iec+2,jsc-2:jec+2,nk))
  enddo
#endif
  tmask_mdfl  = 0.0
  tracer_mdfl = 0.0 
  do n=1,num_prog_tracers
     tracer_mdfl_all(n)%field(:,:,:) = 0.0
  enddo

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           tmask_mdfl(i,j,k) = Grd%tmask(i,j,k)
        enddo
     enddo
  enddo
  call mpp_update_domains(tmask_mdfl,Dom_mdfl%domain2d)  


  ! initialize clock ids 
  id_clock_mdfl_sup_b     = mpp_clock_id('(Ocean advect: MDFL-sup-b)       ',grain=CLOCK_ROUTINE)
  id_clock_mdfl_sweby     = mpp_clock_id('(Ocean advect: MDFL-sweby)       ',grain=CLOCK_ROUTINE)
  id_clock_mdfl_sweby_all = mpp_clock_id('(Ocean advect: MDFL-sweby-all)   ',grain=CLOCK_ROUTINE)

end subroutine mdfl_init
! </SUBROUTINE> NAME="mdfl_init"


!#######################################################################
! <SUBROUTINE NAME="horz_advect_tracer">
!
! <DESCRIPTION>
! Compute horizontal advection of tracers 
! </DESCRIPTION>
!
subroutine horz_advect_tracer(Time, Adv_vel, Thickness, T_prog, Tracer, ntracer, dtime, store_flux)

  type(ocean_time_type), intent(in)           :: Time
  type(ocean_adv_vel_type), intent(in)        :: Adv_vel
  type(ocean_thickness_type), intent(in)      :: Thickness
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)
  type(ocean_prog_tracer_type), intent(inout) :: Tracer
  integer, intent(in)                         :: ntracer
  real, intent(in)                            :: dtime
  logical, intent(in), optional               :: store_flux

  real,dimension(isd:ied,jsd:jed)             :: tmp_flux
  integer                                     :: i, j, k, n
  logical                                     :: store

  store = .TRUE.
  if(present(store_flux)) store = store_flux

  if(.not. advect_sweby_all) then 

      do k=1,nk
         do j=jsd,jed
            do i=isd,ied
               Tracer%wrk1(i,j,k) = 0.0
            enddo
         enddo
      enddo

      select case (Tracer%horz_advect_scheme)

      case (ADVECT_UPWIND)
          Tracer%wrk1(isc:iec,jsc:jec,:) =  -horz_advect_tracer_upwind(Time, Adv_vel, Tracer)
      case (ADVECT_2ND_ORDER)
          Tracer%wrk1(isc:iec,jsc:jec,:) =  -horz_advect_tracer_2nd_order(Time, Adv_vel, Tracer)
      case (ADVECT_4TH_ORDER)
          Tracer%wrk1(isc:iec,jsc:jec,:) =  -horz_advect_tracer_4th_order(Time, Adv_vel, Tracer)
      case (ADVECT_6TH_ORDER)
          Tracer%wrk1(isc:iec,jsc:jec,:) =  -horz_advect_tracer_6th_order(Time, Adv_vel, Tracer)
      case (ADVECT_QUICKER)
          Tracer%wrk1(isc:iec,jsc:jec,:) =  -horz_advect_tracer_quicker(Time, Adv_vel, Tracer)
      case (ADVECT_QUICKMOM3)
          Tracer%wrk1(isc:iec,jsc:jec,:) =  -horz_advect_tracer_quickmom3(Time, Adv_vel, Tracer)
      case (ADVECT_MDFL_SUP_B)
          Tracer%wrk1(isc:iec,jsc:jec,:) =  -advect_tracer_mdfl_sub_b(Time, Adv_vel, Tracer, Thickness, dtime)
      case (ADVECT_MDFL_SWEBY)
          Tracer%wrk1(isc:iec,jsc:jec,:) =  -advect_tracer_mdfl_sweby(Time, Adv_vel, Tracer, Thickness, dtime)

      case default
          call mpp_error(FATAL,'==>Error from ocean_tracer_advect_mod (horz_advect_tracer): invalid advection scheme chosen')
      end select
 
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               Tracer%th_tendency(i,j,k) = Tracer%th_tendency(i,j,k) + Tracer%wrk1(i,j,k)
            enddo
         enddo
      enddo

      ! diagnostics 

      if (id_horz_advect(ntracer) > 0) used = send_data(id_horz_advect(ntracer), Tracer%conversion*Tracer%wrk1(:,:,:), &
           Time%model_time, rmask=Grd%tmask(:,:,:), &
           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
      if (id_adv_flux_x(ntracer) > 0) used = send_data(id_adv_flux_x(ntracer), Tracer%conversion*flux_x(:,:,:), &
           Time%model_time, rmask=Grd%tmask(:,:,:), &
           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
      if (id_adv_flux_y(ntracer) > 0) used = send_data(id_adv_flux_y(ntracer), Tracer%conversion*flux_y(:,:,:), &
           Time%model_time, rmask=Grd%tmask(:,:,:), &
           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

      if (id_adv_flux_x_int_z(ntracer) > 0) then 
          tmp_flux(:,:) = 0.0
          do k=1,nk
             tmp_flux(isc:iec,jsc:jec) = tmp_flux(isc:iec,jsc:jec) +  flux_x(isc:iec,jsc:jec,k) 
          enddo
          used = send_data(id_adv_flux_x_int_z(ntracer), Tracer%conversion*tmp_flux(:,:), &
               Time%model_time, rmask=Grd%tmask(:,:,1), &
               is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
      endif

      if (id_adv_flux_y_int_z(ntracer) > 0) then 
          tmp_flux(:,:) = 0.0
          do k=1,nk
             tmp_flux(isc:iec,jsc:jec) = tmp_flux(isc:iec,jsc:jec) +  flux_y(isc:iec,jsc:jec,k) 
          enddo
          used = send_data(id_adv_flux_y_int_z(ntracer), Tracer%conversion*tmp_flux(:,:), &
               Time%model_time, rmask=Grd%tmask(:,:,1), &
               is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
      endif

  endif   ! endif for (.not. advect_sweby_all)


  if(advect_sweby_all .and. ntracer==1) then 
    call advect_tracer_sweby_all(Time, Adv_vel, T_prog, Thickness, dtime)
  endif 

  if (have_obc .and. store) then
     call mpp_update_domains(flux_x, Dom%domain2d)
     call mpp_update_domains(flux_y, Dom%domain2d)
     call store_ocean_obc_tracer_flux(Tracer,flux_x,flux_y)
     call store_ocean_obc_tracer_advect(flux_x,flux_y)
  endif

end subroutine horz_advect_tracer
! </SUBROUTINE> NAME="horz_advect_tracer"


!#######################################################################
! <SUBROUTINE NAME="vert_advect_tracer">
!
! <DESCRIPTION>

! Compute vertical advection of tracers.  Note that 
! the mdfl-schemes are three-dimensional, with full 
! advection tendency computed in horz_advect_tracer. 
! So no need to do anything here if they are chosen.
!
! </DESCRIPTION>
!
subroutine vert_advect_tracer(Time, Adv_vel, Tracer, ntracer)

  type(ocean_time_type), intent(in)           :: Time
  type(ocean_adv_vel_type), intent(in)        :: Adv_vel
  type(ocean_prog_tracer_type), intent(inout) :: Tracer
  integer, intent(in)                         :: ntracer
  integer                                     :: i, j, k, n

  if(.not. advect_sweby_all) then 

      do k=1,nk
         do j=jsd,jed
            do i=isd,ied
               Tracer%wrk1(i,j,k) = 0.0
            enddo
         enddo
      enddo

      select case (Tracer%vert_advect_scheme)

      case (ADVECT_UPWIND)
          Tracer%wrk1(isc:iec,jsc:jec,:) = -vert_advect_tracer_upwind(Time, Adv_vel, Tracer)
      case (ADVECT_2ND_ORDER)
          Tracer%wrk1(isc:iec,jsc:jec,:) = -vert_advect_tracer_2nd_order(Time, Adv_vel, Tracer)
      case (ADVECT_4TH_ORDER)
          Tracer%wrk1(isc:iec,jsc:jec,:) = -vert_advect_tracer_4th_order(Time, Adv_vel, Tracer)
      case (ADVECT_6TH_ORDER)
          Tracer%wrk1(isc:iec,jsc:jec,:) = -vert_advect_tracer_6th_order(Time, Adv_vel, Tracer)
      case (ADVECT_QUICKER)
          Tracer%wrk1(isc:iec,jsc:jec,:) = -vert_advect_tracer_quicker(Time, Adv_vel, Tracer)
      case (ADVECT_QUICKMOM3)
          Tracer%wrk1(isc:iec,jsc:jec,:) = -vert_advect_tracer_quickmom3(Time, Adv_vel, Tracer)

      ! the mdfl-schemes are three-dimensional, with full 
      ! advection tendency computed in horz_advect_tracer
      case (ADVECT_MDFL_SUP_B)
      case (ADVECT_MDFL_SWEBY)

      case default
        call mpp_error(FATAL,'==>Error from ocean_tracer_advect_mod (vert_advect_tracer): invalid advection scheme chosen')
      end select

      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               Tracer%th_tendency(i,j,k) = Tracer%th_tendency(i,j,k) + Tracer%wrk1(i,j,k)
            enddo
         enddo
      enddo

      ! diagnostics 

      if (id_vert_advect(ntracer) > 0) used = send_data(id_vert_advect(ntracer), Tracer%conversion*Tracer%wrk1(:,:,:), &
           Time%model_time, rmask=Grd%tmask(:,:,:), &
           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

      if (id_adv_flux_z(ntracer) > 0) used = send_data(id_adv_flux_z(ntracer), Tracer%conversion*flux_z(:,:,:), &
                                         Time%model_time, rmask=Grd%tmask(:,:,:), &
                                         is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

  endif  ! endif for (.not. advect_sweby_all)

     
end subroutine vert_advect_tracer
! </SUBROUTINE> NAME="vert_advect_tracer"


!#######################################################################
! <FUNCTION NAME="horz_advect_tracer_upwind">
!
! <DESCRIPTION>
! Compute horizontal advection of tracers from first order upwind.  
! This scheme is positive definite but very diffusive.  
! </DESCRIPTION>
!
function horz_advect_tracer_upwind(Time, Adv_vel, Tracer)

  type(ocean_time_type), intent(in)        :: Time
  type(ocean_adv_vel_type), intent(in)     :: Adv_vel
  type(ocean_prog_tracer_type), intent(in) :: Tracer
  real,dimension(isc:iec,jsc:jec,nk)       :: horz_advect_tracer_upwind
  
  real, dimension(isd:ied,jsd:jed) :: fe, fn 
  real                             :: velocity, upos, uneg
  integer                          :: i, j, k
  integer                          :: taum1

  call mpp_clock_begin(id_clock_up_horz)

  if ( .not. module_is_initialized ) then 
    call  mpp_error(FATAL, &
     '==>Error from ocean_tracer_advect (horz_advect_tracer_upwind): ocean_tracer_advect_mod not yet initialized')
  endif 

  taum1 = Time%taum1
  
  do k=1,nk 

     ! i-flux
     do j=jsc,jec
        do i=isc-1,iec
           velocity = 0.5*Adv_vel%uh_et(i,j,k)
           upos     = velocity + abs(velocity)
           uneg     = velocity - abs(velocity)
           fe(i,j)  = Grd%dyte(i,j)*(upos*Tracer%field(i,j,k,taum1) + uneg*Tracer%field(i+1,j,k,taum1)) &
                      *Grd%tmask(i,j,k)*Grd%tmask(i+1,j,k)  
           flux_x(i,j,k) = fe(i,j)
        enddo
     enddo

     ! j-flux
     do j=jsc-1,jec
        do i=isc,iec
           velocity = 0.5*Adv_vel%vh_nt(i,j,k)
           upos     = velocity + abs(velocity)
           uneg     = velocity - abs(velocity)
           fn(i,j)  = Grd%dxtn(i,j)*(upos*Tracer%field(i,j,k,taum1) + uneg*Tracer%field(i,j+1,k,taum1)) &
                      *Grd%tmask(i,j,k)*Grd%tmask(i,j+1,k)
           flux_y(i,j,k) = fn(i,j)
        enddo
     enddo

     do j=jsc,jec
        do i=isc,iec
           horz_advect_tracer_upwind(i,j,k) = Grd%tmask(i,j,k)*(fe(i,j)-fe(i-1,j)+fn(i,j)-fn(i,j-1))*Grd%datr(i,j)
        enddo
     enddo

  enddo

  call mpp_clock_end(id_clock_up_horz)
  
  
end function horz_advect_tracer_upwind
! </FUNCTION> NAME="horz_advect_tracer_upwind"


!#######################################################################
! <FUNCTION NAME="horz_advect_tracer_2nd_order">
!
! <DESCRIPTION>
! Compute horizontal advection of tracers from
! second order centered differences.  
! </DESCRIPTION>
!
function horz_advect_tracer_2nd_order(Time, Adv_vel, Tracer)

  type(ocean_time_type), intent(in)        :: Time
  type(ocean_adv_vel_type), intent(in)     :: Adv_vel
  type(ocean_prog_tracer_type), intent(in) :: Tracer
  real,dimension(isc:iec,jsc:jec,nk)       :: horz_advect_tracer_2nd_order
  
  real, dimension(isd:ied) :: fs, fn
  integer                  :: i, j, k
  real                     :: fe, fw
  integer                  :: tau

  call mpp_clock_begin(id_clock_2nd_horz)

  if ( .not. module_is_initialized ) then 
    call  mpp_error(FATAL, &
     '==>Error from ocean_tracer_advect (horz_advect_tracer_2nd_order): ocean_tracer_advect_mod not yet initialized')
  endif 

  tau = Time%tau
  
  fs(:) = 0.0; fn(:) = 0.0
  do k=1,nk 
    j     = jsc-1
    fs(:) = Grd%dxtn(:,j)*Adv_vel%vh_nt(:,j,k)*0.5*(Tracer%field(:,j+1,k,tau)+Tracer%field(:,j,k,tau))
    do j=jsc,jec
      i = isc-1
      fw = Grd%dyte(i,j)*Adv_vel%uh_et(i,j,k)*0.5*(Tracer%field(i+1,j,k,tau)+Tracer%field(i,j,k,tau))
      do i=isc,iec
        fe    = Grd%dyte(i,j)*Adv_vel%uh_et(i,j,k)*0.5*(Tracer%field(i+1,j,k,tau)+Tracer%field(i,j,k,tau))
        fn(i) = Grd%dxtn(i,j)*Adv_vel%vh_nt(i,j,k)*0.5*(Tracer%field(i,j+1,k,tau)+Tracer%field(i,j,k,tau))
        flux_x(i,j,k) = fe
        flux_y(i,j,k) = fn(i)
        horz_advect_tracer_2nd_order(i,j,k) = Grd%tmask(i,j,k)*(fe - fw + fn(i) - fs(i))*Grd%datr(i,j)
        fw = fe
      enddo
      fs(:) = fn(:)
    enddo
  enddo

  call mpp_clock_end(id_clock_2nd_horz)

end function horz_advect_tracer_2nd_order
! </FUNCTION> NAME="horz_advect_tracer_2nd_order"



!#######################################################################
! <FUNCTION NAME="horz_advect_tracer_4th_order">
!
! <DESCRIPTION>
! Compute horizontal advection of tracers from fourth order centered
! differences. NOTE: this code does not account for non-uniform grids.   
! </DESCRIPTION>
!
function horz_advect_tracer_4th_order(Time, Adv_vel, Tracer)

  type(ocean_time_type), intent(in)        :: Time
  type(ocean_adv_vel_type), intent(in)     :: Adv_vel
  type(ocean_prog_tracer_type), intent(in) :: Tracer
  real,dimension(isc:iec,jsc:jec,nk)       :: horz_advect_tracer_4th_order

  integer :: i, j, k
  integer :: im1, ip2, jm1, jp2
  integer :: tau

  call mpp_clock_begin(id_clock_4th_horz)

  if ( .not. module_is_initialized ) then 
    call  mpp_error(FATAL, &
     '==>Error from ocean_tracer_advect (horz_advect_tracer_4th_order): ocean_tracer_advect_mod not yet initialized')
  endif 

  tau = Time%tau

  tracer_fourth = 0.0
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           tracer_fourth(i,j,k) = Tracer%field(i,j,k,tau)
        enddo
     enddo
  enddo
  call mpp_update_domains (tracer_fourth, Dom_fourth%domain2d)

  flux_x = 0.0
  flux_y = 0.0

  do k=1,nk

     do j=jsc,jec
        do i=isc-1,iec
           im1   = tmask_fourth(i-1,j,k)*(i-1) + (1.0-tmask_fourth(i-1,j,k))*i
           ip2   = tmask_fourth(i+2,j,k)*(i+2) + (1.0-tmask_fourth(i+2,j,k))*(i+1)    
           flux_x(i,j,k) = Grd%dyte(i,j)*Adv_vel%uh_et(i,j,k)*  &
                (a4*(tracer_fourth(i+1,j,k)+tracer_fourth(i,j,k)) + &
                 b4*(tracer_fourth(ip2,j,k)+tracer_fourth(im1,j,k)))
        enddo
     enddo

     do j=jsc-1,jec
        do i=isc,iec
           jm1   = tmask_fourth(i,j-1,k)*(j-1) + (1.0-tmask_fourth(i,j-1,k))*j
           jp2   = tmask_fourth(i,j+2,k)*(j+2) + (1.0-tmask_fourth(i,j+2,k))*(j+1)    
           flux_y(i,j,k) = Grd%dxtn(i,j)*Adv_vel%vh_nt(i,j,k)*&
                (a4*(tracer_fourth(i,j+1,k)+tracer_fourth(i,j,k)) + &
                 b4*(tracer_fourth(i,jp2,k)+tracer_fourth(i,jm1,k)))
        enddo
     enddo

  enddo

  ! to handle redundancies across the bipolar fold 
  if (Grd%tripolar) call mpp_update_domains(flux_x, flux_y, Dom_flux%domain2d, gridtype=CGRID_NE) 

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           horz_advect_tracer_4th_order(i,j,k)  &
                = Grd%tmask(i,j,k)*(flux_x(i,j,k) - flux_x(i-1,j,k) + flux_y(i,j,k) - flux_y(i,j-1,k))*Grd%datr(i,j)
        enddo
     enddo
  enddo

  call mpp_clock_end(id_clock_4th_horz)

end function horz_advect_tracer_4th_order
! </FUNCTION> NAME="horz_advect_tracer_4th_order"


!#######################################################################
! <FUNCTION NAME="horz_advect_tracer_6th_order">
!
! <DESCRIPTION>
! Compute horizontal advection of tracers from sixth order centered
! differences. NOTE: this code does not account for non-uniform grids.   
! </DESCRIPTION>
!
function horz_advect_tracer_6th_order(Time, Adv_vel, Tracer)

  type(ocean_time_type), intent(in)        :: Time
  type(ocean_adv_vel_type), intent(in)     :: Adv_vel
  type(ocean_prog_tracer_type), intent(in) :: Tracer
  real,dimension(isc:iec,jsc:jec,nk)       :: horz_advect_tracer_6th_order

  real, dimension(isc-3:iec+3,jsc-3:jec+3) :: tracr, mask
  real, dimension(isd:ied,jsd:jed)         :: qn, qe

  integer                                  :: i, j, k
  integer                                  :: im1, ip2, jm1, jp2
  integer                                  :: tau

  call mpp_clock_begin(id_clock_6th_horz)

  if ( .not. module_is_initialized ) then 
    call  mpp_error(FATAL, &
     '==>Error from ocean_tracer_advect (horz_advect_tracer_6th_order): ocean_tracer_advect_mod not yet initialized')
  endif 

  tau = Time%tau

  tracer_sixth=0.0
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           tracer_sixth(i,j,k) = Tracer%field(i,j,k,tau)
        enddo
     enddo
  enddo
  call mpp_update_domains (tracer_sixth, Dom_sixth%domain2d)

  flux_x = 0.0
  flux_y = 0.0
  do k=1,nk

     do j=jsc,jec
        do i=isc-1,iec
           if (tmask_sixth(i-2,j,k)*tmask_sixth(i-1,j,k) == 0 .or. tmask_sixth(i+3,j,k)*tmask_sixth(i+2,j,k) == 0) then
               im1   = tmask_sixth(i-1,j,k)*(i-1) + (1.0-tmask_sixth(i-1,j,k))*i
               ip2   = tmask_sixth(i+2,j,k)*(i+2) + (1.0-tmask_sixth(i+2,j,k))*(i+1)    
               flux_x(i,j,k) = Grd%dyte(i,j)*Adv_vel%uh_et(i,j,k) &
                    *(a4*(tracer_sixth(i+1,j,k)+tracer_sixth(i,j,k)) + b4*(tracer_sixth(ip2,j,k)+tracer_sixth(im1,j,k)))
           else
               flux_x(i,j,k) = Grd%dyte(i,j)*Adv_vel%uh_et(i,j,k) &
                    *(a6*(tracer_sixth(i+1,j,k)+tracer_sixth(i,j,k)) + b6*(tracer_sixth(i+2,j,k)+tracer_sixth(i-1,j,k))&
                    + c6*(tracer_sixth(i+3,j,k)+tracer_sixth(i-2,j,k)))
           endif
        enddo
     enddo

     do j=jsc-1,jec
        do i=isc,iec
           if (tmask_sixth(i,j-2,k)*tmask_sixth(i,j-1,k) == 0 .or. tmask_sixth(i,j+3,k)*tmask_sixth(i,j+2,k) == 0) then
               jm1   = tmask_sixth(i,j-1,k)*(j-1) + (1.0-tmask_sixth(i,j-1,k))*j
               jp2   = tmask_sixth(i,j+2,k)*(j+2) + (1.0-tmask_sixth(i,j+2,k))*(j+1)    
               flux_y(i,j,k) = Grd%dxtn(i,j)*Adv_vel%vh_nt(i,j,k) &
                    *(a4*(tracer_sixth(i,j+1,k)+tracer_sixth(i,j,k)) + b4*(tracer_sixth(i,jp2,k)+tracer_sixth(i,jm1,k)))
           else
               flux_y(i,j,k) = Grd%dxtn(i,j)*Adv_vel%vh_nt(i,j,k) &
                    *(a6*(tracer_sixth(i,j+1,k)+tracer_sixth(i,j,k)) + b6*(tracer_sixth(i,j+2,k)+tracer_sixth(i,j-1,k))&
                    + c6*(tracer_sixth(i,j+3,k)+tracer_sixth(i,j-2,k)))
           endif
        enddo
     enddo

  enddo

  ! to handle redundancies across the bipolar fold 
  if (Grd%tripolar) call mpp_update_domains(flux_x, flux_y, Dom_flux%domain2d, gridtype=CGRID_NE) 

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           flux_x(i,j,k) = qe(i,j)
           flux_y(i,j,k) = qn(i,j)
           horz_advect_tracer_6th_order(i,j,k) &
                = Grd%tmask(i,j,k)*(flux_x(i,j,k) - flux_x(i-1,j,k) + flux_y(i,j,k) - flux_y(i,j-1,k))*Grd%datr(i,j)
        enddo
     enddo
  enddo

  call mpp_clock_end(id_clock_6th_horz)

end function horz_advect_tracer_6th_order
! </FUNCTION> NAME="horz_advect_tracer_6th_order"


!#######################################################################
! <FUNCTION NAME="horz_advect_tracer_quicker">
!
! <DESCRIPTION>
! Compute horizontal advection of tracers from quicker. 
! </DESCRIPTION>
!
function horz_advect_tracer_quicker(Time, Adv_vel, Tracer)

  type(ocean_time_type), intent(in)        :: Time
  type(ocean_adv_vel_type), intent(in)     :: Adv_vel
  type(ocean_prog_tracer_type), intent(in) :: Tracer
  real,dimension(isc:iec,jsc:jec,nk)       :: horz_advect_tracer_quicker
  real, dimension(isd:ied,jsd:jed)         :: qn, qe
  integer                                  :: i, j, k
  integer                                  :: tau, taum1
  real                                     :: fe, fw, rnormsk,soutmsk, eastmsk, westmsk
  real                                     :: upos, uneg, tm1, tp2, vel

  call mpp_clock_begin(id_clock_quick_horz)

  if ( .not. module_is_initialized ) then 
    call  mpp_error(FATAL, &
     '==>Error from ocean_tracer_advect (horz_advect_tracer_quicker): ocean_tracer_advect_mod not yet initialized')
  endif 

  taum1 = Time%taum1
  tau   = Time%tau

  tracer_quick = 0.0
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           tracer_quick(i,j,k) = Tracer%field(i,j,k,taum1)
        enddo
     enddo
  enddo
  call mpp_update_domains (tracer_quick, Dom_quicker%domain2d)

  flux_x = 0.0
  flux_y = 0.0
  do k=1,nk

     do j=jsc-1,jec
        do i=isc,iec
           vel  = Grd%dxtn(i,j)*Adv_vel%vh_nt(i,j,k)
           upos = 0.5*(vel + abs(vel))
           uneg = 0.5*(vel - abs(vel))
           rnormsk = tmask_quick(i,j,k)*(1-tmask_quick(i,j-1,k))
           soutmsk = tmask_quick(i,j+1,k)*(1-tmask_quick(i,j+2,k))
           flux_y(i,j,k) = vel*(quick_y(i,j,1)*Tracer%field(i,j,k,tau)+&
                quick_y(i,j,2)*Tracer%field(i,j+1,k,tau))&
                - upos*(curv_yp(i,j,1)*tracer_quick(i,j+1,k) + &
                curv_yp(i,j,2)*tracer_quick(i,j,k)   + &
                curv_yp(i,j,3)*(tracer_quick(i,j-1,k)*(1-rnormsk) + tracer_quick(i,j,k)*rnormsk)) &
                - uneg*(curv_yn(i,j,1)*(tracer_quick(i,j+2,k)*(1-soutmsk)+tracer_quick(i,j+1,k)*soutmsk) + &
                curv_yn(i,j,2)*tracer_quick(i,j+1,k)   + &
                curv_yn(i,j,3)*tracer_quick(i,j,k))        
        enddo
     enddo

     do j=jsc,jec
        do i=isc-1,iec
           vel  = Grd%dyte(i,j)*Adv_vel%uh_et(i,j,k)
           upos = 0.5*(vel + abs(vel))
           uneg = 0.5*(vel - abs(vel))
           eastmsk = tmask_quick(i,j,k)*(1-tmask_quick(i-1,j,k))
           westmsk = tmask_quick(i+1,j,k)*(1-tmask_quick(i+2,j,k))
           flux_x(i,j,k) = vel*(quick_x(i,j,1)*Tracer%field(i,j,k,tau)+quick_x(i,j,2)*Tracer%field(i+1,j,k,tau))&
                - upos*(curv_xp(i,j,1)*tracer_quick(i+1,j,k) &
                + curv_xp(i,j,2)*tracer_quick(i,j,k) + &
                curv_xp(i,j,3)*(tracer_quick(i-1,j,k)*(1-eastmsk)+tracer_quick(i,j,k)*eastmsk)) &
                - uneg*(curv_xn(i,j,1)*(tracer_quick(i+2,j,k)*(1-westmsk)+tracer_quick(i+1,j,k)*westmsk) &
                +       curv_xn(i,j,2)*tracer_quick(i+1,j,k) &
                +       curv_xn(i,j,3)*tracer_quick(i,j,k))
        enddo
     enddo


     ! recompute fluxes as upwind if tmask_limit flag satisfied 
     if(limit_with_upwind) then 
         do j=jsc,jec
            do i=isc-1,iec
               if(Tracer%tmask_limit(i,j,k)==1.0) then   
                   vel  = Adv_vel%uh_et(i,j,k)
                   upos = 0.5*(vel + abs(vel))
                   uneg = 0.5*(vel - abs(vel))
                   flux_x(i,j,k)  = Grd%dyte(i,j)*(upos*Tracer%field(i,j,k,taum1) + uneg*Tracer%field(i+1,j,k,taum1)) &
                        *Grd%tmask(i,j,k)*Grd%tmask(i+1,j,k)  
               endif
            enddo
         enddo
         do j=jsc-1,jec
            do i=isc,iec
               if(Tracer%tmask_limit(i,j,k)==1.0) then   
                   vel  = Adv_vel%vh_nt(i,j,k)
                   upos = 0.5*(vel + abs(vel))
                   uneg = 0.5*(vel - abs(vel))
                   flux_y(i,j,k) = Grd%dxtn(i,j)*(upos*Tracer%field(i,j,k,taum1) + uneg*Tracer%field(i,j+1,k,taum1)) &
                        *Grd%tmask(i,j,k)*Grd%tmask(i,j+1,k) 
               endif
            enddo
         enddo
     endif

  enddo  ! end of k-loop 

  ! to handle redundancies across the bipolar fold 
  if (Grd%tripolar) call mpp_update_domains(flux_x, flux_y, Dom_flux%domain2d, gridtype=CGRID_NE) 

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           horz_advect_tracer_quicker(i,j,k) = &
              Grd%tmask(i,j,k)*(flux_x(i,j,k) - flux_x(i-1,j,k) + flux_y(i,j,k) - flux_y(i,j-1,k))*Grd%datr(i,j)
        enddo
     enddo
  enddo

  call mpp_clock_end(id_clock_quick_horz)

end function horz_advect_tracer_quicker
! </FUNCTION> NAME="horz_advect_tracer_quicker"


!#######################################################################
! <FUNCTION NAME="horz_advect_tracer_quickmom3">
!
! <DESCRIPTION>
! Compute horizontal advection of tracers from quicker using MOM3 masking. 
! This method has proven useful for reproducing MOM3 results with MOM4.  
! </DESCRIPTION>
!
function horz_advect_tracer_quickmom3(Time, Adv_vel, Tracer)

  type(ocean_time_type), intent(in)        :: Time
  type(ocean_adv_vel_type), intent(in)     :: Adv_vel
  type(ocean_prog_tracer_type), intent(in) :: Tracer
  real,dimension(isc:iec,jsc:jec,nk)       :: horz_advect_tracer_quickmom3
  real, dimension(isd:ied,jsd:jed)         :: qn, qe
  integer                                  :: i, j, k, tau, taum1, jp1, jp2, jp2_mod
  real                                     :: fe, fw, rnormsk,soutmsk, eastmsk, westmsk
  real                                     :: upos, uneg, tm1, tp2, vel

  call mpp_clock_begin(id_clock_quickmom3_horz)

  if ( .not. module_is_initialized ) then 
    call  mpp_error(FATAL, &
     '==>Error from ocean_tracer_advect (horz_advect_tracer_quickmom3): ocean_tracer_advect_mod not yet initialized')
  endif 

  taum1 = Time%taum1
  tau   = Time%tau
  
  tracer_quick = 0.0
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           tracer_quick(i,j,k) = Tracer%field(i,j,k,taum1)
        enddo
     enddo
  enddo
  call mpp_update_domains (tracer_quick, Dom_quicker%domain2d)

  flux_x = 0.0
  flux_y = 0.0

  do k=1,nk

     do j=jsc-1,jec
        do i=isc,iec
           jp1 = min(j+1, jec)
           jp2 = min(j+2, jec)
           jp2_mod = nint(tmask_quick(i,jp2,k)*jp2 + (1.0-tmask_quick(i,jp2,k))*jp1)

           vel  = Grd%dxtn(i,j)*Adv_vel%vh_nt(i,j,k)
           upos = 0.5*(vel + abs(vel))*tmask_quick(i,j-1,k)*tmask_quick(i,j,k)*tmask_quick(i,j+1,k)
           uneg = 0.5*(vel - abs(vel))*tmask_quick(i,jp2,k)*tmask_quick(i,j+1,k)*tmask_quick(i,j,k)

           flux_y(i,j,k) = vel*(quick_y(i,j,1)*Tracer%field(i,j,k,tau)+  &
                quick_y(i,j,2)*Tracer%field(i,j+1,k,tau))                &
                - upos*(curv_yp(i,j,1)*tracer_quick(i,j+1,k) +           &
                curv_yp(i,j,2)*tracer_quick(i,j,k)   +                   &
                curv_yp(i,j,3)*tracer_quick(i,j-1,k))                    &
                - uneg*(curv_yn(i,j,1)*tracer_quick(i,jp2,k) +           &
                curv_yn(i,j,2)*tracer_quick(i,j+1,k)   +                 &
                curv_yn(i,j,3)*tracer_quick(i,j,k))        
        enddo
     enddo

     do j=jsc,jec
        do i=isc-1,iec
           vel  = Grd%dyte(i,j)*Adv_vel%uh_et(i,j,k)
           upos = 0.5*(vel + abs(vel))*tmask_quick(i-1,j,k)*tmask_quick(i,j,k)*tmask_quick(i+1,j,k)
           uneg = 0.5*(vel - abs(vel))*tmask_quick(i+2,j,k)*tmask_quick(i+1,j,k)*tmask_quick(i,j,k)

           flux_x(i,j,k) = vel*(quick_x(i,j,1)*Tracer%field(i,j,k,tau)+ &
                quick_x(i,j,2)*Tracer%field(i+1,j,k,tau))               &
                - upos*(curv_xp(i,j,1)*tracer_quick(i+1,j,k)            &
                + curv_xp(i,j,2)*tracer_quick(i,j,k) +                  &
                curv_xp(i,j,3)*tracer_quick(i-1,j,k))                   &
                - uneg*(curv_xn(i,j,1)*tracer_quick(i+2,j,k)            &
                +       curv_xn(i,j,2)*tracer_quick(i+1,j,k)            &
                +       curv_xn(i,j,3)*tracer_quick(i,j,k))
        enddo
     enddo


     ! recompute fluxes as upwind if tmask_limit flag satisfied 
     if(limit_with_upwind) then 
         do j=jsc,jec
            do i=isc-1,iec
               if(Tracer%tmask_limit(i,j,k)==1.0) then   
                   vel  = Adv_vel%uh_et(i,j,k)
                   upos = 0.5*(vel + abs(vel))
                   uneg = 0.5*(vel - abs(vel))
                   flux_x(i,j,k) = Grd%dyte(i,j)*(upos*Tracer%field(i,j,k,taum1) + uneg*Tracer%field(i+1,j,k,taum1)) &
                        *Grd%tmask(i,j,k)*Grd%tmask(i+1,j,k)  
               endif
            enddo
         enddo
         do j=jsc-1,jec
            do i=isc,iec
               if(Tracer%tmask_limit(i,j,k)==1.0) then   
                   vel  = Adv_vel%vh_nt(i,j,k)
                   upos = 0.5*(vel + abs(vel))
                   uneg = 0.5*(vel - abs(vel))
                   flux_y(i,j,k) = Grd%dxtn(i,j)*(upos*Tracer%field(i,j,k,taum1) + uneg*Tracer%field(i,j+1,k,taum1)) &
                        *Grd%tmask(i,j,k)*Grd%tmask(i,j+1,k) 
               endif
            enddo
         enddo
     endif

  enddo ! end of k-loop

  ! to handle redundancies across the bipolar fold 
  if (Grd%tripolar) call mpp_update_domains(flux_x, flux_y, Dom_flux%domain2d, gridtype=CGRID_NE) 

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           horz_advect_tracer_quickmom3(i,j,k) = &
                Grd%tmask(i,j,k)*(flux_x(i,j,k) - flux_x(i-1,j,k) + flux_y(i,j,k) - flux_y(i,j-1,k))*Grd%datr(i,j)
        enddo
     enddo
  enddo

  call mpp_clock_end(id_clock_quickmom3_horz)

end function horz_advect_tracer_quickmom3
! </FUNCTION> NAME="horz_advect_tracer_quickmom3"


!#######################################################################
! <FUNCTION NAME="vert_advect_tracer_upwind">
!
! <DESCRIPTION>
! Compute vertical advection of tracers from first order upwind. 
! This scheme is positive definite, but very diffusive. 
! </DESCRIPTION>
!
function vert_advect_tracer_upwind(Time, Adv_vel, Tracer)

  type(ocean_time_type), intent(in)        :: Time
  type(ocean_adv_vel_type), intent(in)     :: Adv_Vel
  type(ocean_prog_tracer_type), intent(in) :: Tracer
  real,dimension(isc:iec, jsc:jec, nk)     :: vert_advect_tracer_upwind
  real,dimension(isc:iec, jsc:jec)         :: ft1, ft2 
  real                                     :: velocity, wpos, wneg
  integer                                  :: k, i, j, kp1
  integer                                  :: taum1

  call mpp_clock_begin(id_clock_up_vert)

  taum1 = Time%taum1
 
  ft1  = 0.0
  do k=1,nk
    kp1 = min(k+1,nk)
    do j=jsc,jec
      do i=isc,iec
        velocity = 0.5*Adv_vel%w_bt(i,j,k)
        wpos     = velocity + abs(velocity) 
        wneg     = velocity - abs(velocity) 
        ft2(i,j) = (wneg*Tracer%field(i,j,k,taum1) + wpos*Tracer%field(i,j,kp1,taum1)) &
                   *Grd%tmask(i,j,k)*Grd%tmask(i,j,kp1) 
        flux_z(i,j,k) = Grd%dat(i,j)*ft2(i,j)
        vert_advect_tracer_upwind(i,j,k) = Grd%tmask(i,j,k)*(ft1(i,j)-ft2(i,j))
        ft1(i,j) = ft2(i,j)
      enddo
    enddo
  enddo

  call mpp_clock_end(id_clock_up_vert)

end function vert_advect_tracer_upwind
! </FUNCTION> NAME="vert_advect_tracer_upwind"


!#######################################################################
! <FUNCTION NAME="vert_advect_tracer_2nd_order">
!
! <DESCRIPTION>
! Compute vertical advection of tracers from second order centered
! differences. 
! </DESCRIPTION>
!
function vert_advect_tracer_2nd_order(Time, Adv_vel, Tracer)

  type(ocean_time_type), intent(in)        :: Time
  type(ocean_adv_vel_type), intent(in)     :: Adv_Vel
  type(ocean_prog_tracer_type), intent(in) :: Tracer
  real,dimension(isc:iec, jsc:jec, nk)     :: vert_advect_tracer_2nd_order 
  integer                                  :: k, i, j, kp
  real,dimension(isc:iec,jsc:jec)          :: ft1, ft2
  integer                                  :: tau

  call mpp_clock_begin(id_clock_2nd_vert)

  tau = Time%tau
  ft1 = 0.0
  
  do k=1,nk
    kp = min(k+1,nk)
    ft2(isc:iec,jsc:jec) = Adv_vel%w_bt(isc:iec,jsc:jec,k)*0.5*(Tracer%field(isc:iec,jsc:jec,k,tau)+&
                           Tracer%field(isc:iec,jsc:jec,kp,tau))
    flux_z(isc:iec,jsc:jec,k) = Grd%dat(isc:iec,jsc:jec)*ft2(isc:iec,jsc:jec)
    vert_advect_tracer_2nd_order(isc:iec,jsc:jec,k) = Grd%tmask(isc:iec,jsc:jec,k)* &
                                                         (ft1(isc:iec,jsc:jec)-ft2(isc:iec,jsc:jec))
    ft1(:,:) = ft2(:,:)
  enddo

  call mpp_clock_end(id_clock_2nd_vert)

end function vert_advect_tracer_2nd_order
! </FUNCTION> NAME="vert_advect_tracer_2nd_order"


!#######################################################################
! <FUNCTION NAME="vert_advect_tracer_4th_order">
!
! <DESCRIPTION>
! Compute vertical advection of tracers from fourth order centered
! differences. NOTE: this code does not account for non-uniform grids.   
! </DESCRIPTION>
!
function vert_advect_tracer_4th_order(Time, Adv_vel, Tracer)

  type(ocean_time_type), intent(in)        :: Time
  type(ocean_adv_vel_type), intent(in)     :: Adv_Vel
  type(ocean_prog_tracer_type), intent(in) :: Tracer
  real,dimension(isc:iec, jsc:jec, nk)     :: vert_advect_tracer_4th_order 
  integer                                  :: k, i, j, kp2, km1, p1, p2
  real,dimension(isc:iec,jsc:jec)          :: ft1, ft2
  integer                                  :: tau

  call mpp_clock_begin(id_clock_4th_vert)

  tau = Time%tau

  ft1(isc:iec,jsc:jec) = 0.0
  do k=1,nk
    km1 = max(k-1,1)
    p1  = min(k+1,nk)
    p2  = min(k+2,nk)
    do j=jsc,jec
      do i=isc,iec
        kp2   = nint(Grd%tmask(i,j,p2)*p2 + (1.0-Grd%tmask(i,j,p2))*p1)    
        ft2(i,j) = Adv_vel%w_bt(i,j,k)*(a4*(Tracer%field(i,j,k,tau)+&
                   Tracer%field(i,j,p1,tau)) + b4*(Tracer%field(i,j,km1,tau)+&
                   Tracer%field(i,j,kp2,tau)))
        flux_z(i,j,k) = Grd%dat(i,j)*ft2(i,j)
        vert_advect_tracer_4th_order(i,j,k) = Grd%tmask(i,j,k)*(ft1(i,j)-ft2(i,j))
      enddo
    enddo
    ft1(isc:iec,jsc:jec) = ft2(isc:iec,jsc:jec)
  enddo

  call mpp_clock_end(id_clock_4th_vert)

end function vert_advect_tracer_4th_order
! </FUNCTION> NAME="vert_advect_tracer_4th_order"


!#######################################################################
! <FUNCTION NAME="vert_advect_tracer_6th_order">
!
! <DESCRIPTION>
! Compute vertical advection of tracers from sixth order centered
! differences. NOTE: this code does not account for non-uniform grids.   
! </DESCRIPTION>
!
function vert_advect_tracer_6th_order(Time, Adv_vel, Tracer)

  type(ocean_time_type), intent(in)        :: Time
  type(ocean_adv_vel_type), intent(in)     :: Adv_vel
  type(ocean_prog_tracer_type), intent(in) :: Tracer

  real,dimension(isc:iec,jsc:jec)    :: ft1, ft2
  real,dimension(isc:iec,jsc:jec,nk) :: vert_advect_tracer_6th_order 
  integer                            :: i, j, k, km1, kp2, p1, p2, tau

  call mpp_clock_begin(id_clock_6th_vert)

  tau = Time%tau

  ft1(isc:iec,jsc:jec) = 0.0
  do k=1,nk
    km1 = max(k-1,1)
    p1 = min(k+1,nk)
    p2 = min(k+2,nk)
    if (k <= 2 .or. k >= nk-1) then
      do j=jsc,jec
        do i=isc,iec
          kp2   = nint(Grd%tmask(i,j,p2)*p2 + (1.0-Grd%tmask(i,j,p2))*p1)    
          ft2(i,j) = Adv_vel%w_bt(i,j,k)*(a4*(Tracer%field(i,j,k,tau)+&
                     Tracer%field(i,j,p1,tau)) + &
                     b4*(Tracer%field(i,j,km1,tau)+Tracer%field(i,j,kp2,tau)))
          flux_z(i,j,k) = Grd%dat(i,j)*ft2(i,j)          
          vert_advect_tracer_6th_order(i,j,k) = Grd%tmask(i,j,k)*(ft1(i,j)-ft2(i,j))
        enddo
      enddo
    else
      do j=jsc,jec
        do i=isc,iec
          kp2   = nint(Grd%tmask(i,j,p2)*p2 + (1.0-Grd%tmask(i,j,p2))*p1)    
          if (k < Grd%kmt(i,j)-2) then
            ft2(i,j) = Adv_vel%w_bt(i,j,k)*(a6*(Tracer%field(i,j,k,tau)+&
                       Tracer%field(i,j,p1,tau)) + &
                       b6*(Tracer%field(i,j,km1,tau)+ &
                       Tracer%field(i,j,kp2,tau))&
                       + c6*(Tracer%field(i,j,k-2,tau)+ &
                       Tracer%field(i,j,k+3,tau)))
          else
            ft2(i,j) = Adv_vel%w_bt(i,j,k)*(a4*(Tracer%field(i,j,k,tau)+&
                       Tracer%field(i,j,p1,tau)) +&
                       b4*(Tracer%field(i,j,km1,tau)+&
                       Tracer%field(i,j,kp2,tau)))
        endif
          flux_z(i,j,k) = Grd%dat(i,j)*ft2(i,j)        
          vert_advect_tracer_6th_order(i,j,k) = Grd%tmask(i,j,k)*(ft1(i,j)-ft2(i,j))
        enddo
      enddo
    endif
    ft1(isc:iec,jsc:jec) = ft2(isc:iec,jsc:jec)
  enddo

  call mpp_clock_end(id_clock_6th_vert)

end function vert_advect_tracer_6th_order
! </FUNCTION> NAME="vert_advect_tracer_6th_order"


!#######################################################################
! <FUNCTION NAME="vert_advect_tracer_quicker">
!
! <DESCRIPTION>
! Compute vertical advection of tracers from quicker.
! </DESCRIPTION>
!
function vert_advect_tracer_quicker(Time, Adv_vel, Tracer)

  type(ocean_time_type), intent(in)        :: Time
  type(ocean_adv_vel_type), intent(in)     :: Adv_vel
  type(ocean_prog_tracer_type), intent(in) :: Tracer
  real,dimension(isc:iec,jsc:jec,nk)       :: vert_advect_tracer_quicker 
  integer                                  :: k, i, j, kp, tau, taum1
  real,dimension(isc:iec,jsc:jec)          :: ft1, ft2
  real                                     :: upos, uneg, tm1, tp2, vel, upmsk
  integer                                  :: km1, kp2, kp1, p2
  
  call mpp_clock_begin(id_clock_quick_vert)

  tau   = Time%tau
  taum1 = Time%taum1

  ft1(isc:iec,jsc:jec) = 0.0
  do k=1,nk
    km1 = max(k-1,1)
    kp1 = min(k+1,nk)
    p2  = min(k+2,nk)
    do j=jsc,jec
      do i=isc,iec
        vel  = Adv_vel%w_bt(i,j,k)
        upos = 0.5*(vel + abs(vel))
        uneg = 0.5*(vel - abs(vel))
        if(Tracer%tmask_limit(i,j,k)==1.0) then 
          ft2(i,j) = (uneg*Tracer%field(i,j,k,taum1) + upos*Tracer%field(i,j,kp1,taum1)) &
                          *Grd%tmask(i,j,k)*Grd%tmask(i,j,kp1)
        else 
          kp2   = nint(Grd%tmask(i,j,p2)*p2 + (1.0-Grd%tmask(i,j,p2))*kp1)
          upmsk = Grd%tmask(i,j,k)*(1-Grd%tmask(i,j,km1))
          ft2(i,j) = vel*(quick_z(k,1)*Tracer%field(i,j,k,tau) + &
                        quick_z(k,2)*Tracer%field(i,j,kp1,tau)) &
                    -uneg*(curv_zp(k,1)*Tracer%field(i,j,kp1,taum1) &
                          +curv_zp(k,2)*Tracer%field(i,j,k,taum1) &
                          +curv_zp(k,3)*Tracer%field(i,j,km1,taum1)) &
                    -upos*(curv_zn(k,1)*Tracer%field(i,j,kp2,taum1) &
                          +curv_zn(k,2)*Tracer%field(i,j,kp1,taum1) &
                          +curv_zn(k,3)*(Tracer%field(i,j,k,taum1)*(1-upmsk)+Tracer%field(i,j,kp1,taum1)*upmsk))

        endif 
        flux_z(i,j,k) = Grd%dat(i,j)*ft2(i,j)
        vert_advect_tracer_quicker(i,j,k) = Grd%tmask(i,j,k)*(ft1(i,j)-ft2(i,j))
      enddo
    enddo
    ft1(isc:iec,jsc:jec) = ft2(isc:iec,jsc:jec)
  enddo

  call mpp_clock_end(id_clock_quick_vert)
  
end function vert_advect_tracer_quicker
! </FUNCTION> NAME="vert_advect_tracer_quicker"


!#######################################################################
! <FUNCTION NAME="vert_advect_tracer_quickmom3">
!
! <DESCRIPTION>
! Compute vertical advection of tracers from quicker using MOM3 masking.
! This method has proven useful for reproducing MOM3 results with MOM4.  
! </DESCRIPTION>
!
function vert_advect_tracer_quickmom3(Time, Adv_vel, Tracer)

  type(ocean_time_type), intent(in)        :: Time
  type(ocean_adv_vel_type), intent(in)     :: Adv_vel
  type(ocean_prog_tracer_type), intent(in) :: Tracer
  real,dimension(isc:iec,jsc:jec,nk)       :: vert_advect_tracer_quickmom3 
  integer                                  :: k, i, j, kp, tau, taum1
  real,dimension(isc:iec,jsc:jec)          :: ft1, ft2
  real                                     :: upos, uneg, tm1, tp2, vel, upmsk
  integer                                  :: kp1
  
  call mpp_clock_begin(id_clock_quickmom3_vert)

  tau   = Time%tau
  taum1 = Time%taum1

  ft1(isc:iec,jsc:jec) = 0.0
  do k=1,nk

     if (k == 1) then
      do j=jsc,jec
        do i=isc,iec
          vel  = Adv_vel%w_bt(i,j,k)
          upos = 0.5*(vel + abs(vel))*Grd%tmask(i,j,k+2)*Grd%tmask(i,j,k+1)*Grd%tmask(i,j,k)
          if(Tracer%tmask_limit(i,j,k)==1.0) then 
            kp1 = min(k+1,nk)
          ft2(i,j) = (uneg*Tracer%field(i,j,k,taum1) +               &
                      upos*Tracer%field(i,j,kp1,taum1))              &
                     *Grd%tmask(i,j,k)*Grd%tmask(i,j,kp1)
          else 
            ft2(i,j) = vel*(quick_z(k,1)*Tracer%field(i,j,k,tau) +   &
                        quick_z(k,2)*Tracer%field(i,j,k+1,tau))      &
                    -upos*(curv_zn(k,1)*Tracer%field(i,j,k+2,taum1)  &
                          +curv_zn(k,2)*Tracer%field(i,j,k+1,taum1)  &
                          +curv_zn(k,3)*Tracer%field(i,j,k,taum1))
          endif 
         flux_z(i,j,k) = Grd%dat(i,j)*ft2(i,j)
         vert_advect_tracer_quickmom3(i,j,k) = Grd%tmask(i,j,k)*(ft1(i,j)-ft2(i,j))
       enddo
     enddo

    elseif (k == nk-1) then
      do j=jsc,jec
        do i=isc,iec
         vel  = Adv_vel%w_bt(i,j,k)
         uneg = 0.5*(vel - abs(vel))*Grd%tmask(i,j,k-1)*Grd%tmask(i,j,k)*Grd%tmask(i,j,k+1)
          if(Tracer%tmask_limit(i,j,k)==1.0) then 
            kp1 = min(k+1,nk)
            ft2(i,j) = (uneg*Tracer%field(i,j,k,taum1) +     &
                        upos*Tracer%field(i,j,kp1,taum1))    &
                        *Grd%tmask(i,j,k)*Grd%tmask(i,j,kp1)
          else
            ft2(i,j) = vel*(quick_z(k,1)*Tracer%field(i,j,k,tau) +   &
                        quick_z(k,2)*Tracer%field(i,j,k+1,tau))      &
                    -uneg*(curv_zp(k,1)*Tracer%field(i,j,k+1,taum1)  &
                          +curv_zp(k,2)*Tracer%field(i,j,k,taum1)    &
                          +curv_zp(k,3)*Tracer%field(i,j,k-1,taum1))
          endif 
         flux_z(i,j,k) = Grd%dat(i,j)*ft2(i,j)
         vert_advect_tracer_quickmom3(i,j,k) = Grd%tmask(i,j,k)*(ft1(i,j)-ft2(i,j))
        enddo
      enddo

    else
      do j=jsc,jec
        do i=isc,iec
         vel  = Adv_vel%w_bt(i,j,k)
         upos = 0.5*(vel + abs(vel))*Grd%tmask(i,j,k+2)*Grd%tmask(i,j,k+1)*Grd%tmask(i,j,k)
         uneg = 0.5*(vel - abs(vel))*Grd%tmask(i,j,k-1)*Grd%tmask(i,j,k)*Grd%tmask(i,j,k+1)
          if(Tracer%tmask_limit(i,j,k)==1.0) then 
            kp1 = min(k+1,nk)
            ft2(i,j) = (uneg*Tracer%field(i,j,k,taum1) +     &
                        upos*Tracer%field(i,j,kp1,taum1))    &
                        *Grd%tmask(i,j,k)*Grd%tmask(i,j,kp1)
          else


            ft2(i,j) = vel*(quick_z(k,1)*Tracer%field(i,j,k,tau) +   &
                        quick_z(k,2)*Tracer%field(i,j,k+1,tau))      &
                    -upos*(curv_zn(k,1)*Tracer%field(i,j,k+2,taum1)  &
                          +curv_zn(k,2)*Tracer%field(i,j,k+1,taum1)  &
                          +curv_zn(k,3)*Tracer%field(i,j,k,taum1))   &
                    -uneg*(curv_zp(k,1)*Tracer%field(i,j,k+1,taum1)  &
                          +curv_zp(k,2)*Tracer%field(i,j,k,taum1)    &
                          +curv_zp(k,3)*Tracer%field(i,j,k-1,taum1))
          endif 
         flux_z(i,j,k) = Grd%dat(i,j)*ft2(i,j)
         vert_advect_tracer_quickmom3(i,j,k) = Grd%tmask(i,j,k)*(ft1(i,j)-ft2(i,j))
        enddo
      enddo
     endif


    ft1(isc:iec,jsc:jec) = ft2(isc:iec,jsc:jec)
  enddo

  call mpp_clock_end(id_clock_quickmom3_vert)
  
end function vert_advect_tracer_quickmom3
! </FUNCTION> NAME="vert_advect_tracer_quickmom3"


!#######################################################################
! <FUNCTION NAME="advect_tracer_mdfl_sub_b">
!
! <DESCRIPTION>
! Compute tendency due to 3D advection of tracers using a multi-dimensional flux-limited
! method.  This method differs from other methods in the following ways:
!
! 1) Horizontal and vertical advection are combined
! 2) Calculations of the three coordinates (Z, X, Y) are performed sequentially as updates
!        to the tracer, so that the advection components for X and Y depend on Z and Z,Y
!        respectively... This helps limit the flux.
! 3) During the update for each direction, the 2nd order Super-B flux limiter is applied.
! 4) Flux divergence is included within the calculation to also help limit the flux:
!          - During the update for each direction, the divergence in each direction is added.
!          - During the overall tendency calculated, the divergence in all three directions 
!                 is removed.
! 5) All fluxes are functions of the tracer field at the taum1 time step.  This 
!    means that this method is ideally suited for the twolevel time stepping scheme,
!    in which "taum1=tau", thus enabling twice the tracer time step available for the 
!    threelevel scheme. 
!
!The calculation proceeds as follows:
!
! IMPORTANT NOTE: If this scheme is used at all, it must be used as the option for BOTH
! horizontal and vertical advection.  In the the tracer tendency, it is applied as the 
! horizontal term, but applies to vertical as well, for which case the vertical term in
! the tracer tendency equation is set to zero.
!
! This scheme was ported to mom4 from the MIT-GCM by John Dunne and Alistair Adcroft  
! during Summer 2003 
!
! smg: 3/3/2004
!      need to generalize to thickness weighted algorithm.  Presently is in terms
!      of pure advection tendency.  Not clean, though may be workable.  
!
! </DESCRIPTION>
!
function advect_tracer_mdfl_sub_b(Time, Adv_vel, Tracer, Thickness, dtime)

  type(ocean_time_type), intent(in)        :: Time
  type(ocean_adv_vel_type), intent(in)     :: Adv_vel
  type(ocean_prog_tracer_type), intent(in) :: Tracer
  type(ocean_thickness_type), intent(in)   :: Thickness
  real, intent(in)                         :: dtime

  real,dimension(isc  :iec,jsc-1:jec)      :: fn
  real,dimension(isc-1:iec,jsc  :jec)      :: fe
  real,dimension(isc:iec,jsc:jec,nk)       :: advect_tracer_mdfl_sub_b 
  integer                                  :: k, i, j, kp1, kp2, km1, tau, taum1
  real,dimension(isc:iec,jsc:jec)          :: ftp, fbt, wkm1
  real,dimension(isc-2:iec+2,jsc-2:jec+2)  :: tmp_updated_tracer
  real                                     :: Rjm, Rj, Rjp, Cr, cfl, volflux

  call mpp_clock_begin(id_clock_mdfl_sup_b)

  ! initialize to zero local arrays 
  tmp_updated_tracer = 0.0
  ftp                = 0.0
  fbt                = 0.0
  wkm1               = 0.0
  fn                 = 0.0
  fe                 = 0.0

  ! set time indices 
  tau   = Time%tau
  taum1 = Time%taum1

!
!Boundary condition at surface
!
    do j=jsc,jec
      do i=isc,iec
        ftp(i,j) = 0.0
        wkm1(i,j) = 0.0
      enddo !i
    enddo !j
!
! Main loop over all depths
!
  do k=1,nk

!
! Get tracer field on computational grid
!
    do j=jsc,jec
      do i=isc,iec
        tmp_updated_tracer(i,j) = Tracer%field(i,j,k,taum1)
      enddo !i
    enddo !j
!
! Account for boundaries
!
    kp1 = min(k+1,nk)
    kp2 = min(k+2,nk)
    km1 = max(k-1,1)
!
! Calculate flux at bottom of box and update the tracer
!
    do j=jsc,jec
      do i=isc,iec
        Rjp = ( Tracer%field(i,j,km1,taum1) - Tracer%field(i,j,k,taum1) )  &
                  * tmask_mdfl(i,j,km1) *       tmask_mdfl(i,j,k)
        Rj  = ( Tracer%field(i,j,k,taum1) - Tracer%field(i,j,kp1,taum1) )    &
                  * tmask_mdfl(i,j,k) *       tmask_mdfl(i,j,kp1)
        Rjm  = ( Tracer%field(i,j,kp1,taum1) - Tracer%field(i,j,kp2,taum1) )   &
                  * tmask_mdfl(i,j,kp1) *        tmask_mdfl(i,j,kp2)
        volflux = Grd%dat(i,j) * Adv_vel%w_bt(i,j,k)
        cfl = abs(Adv_vel%w_bt(i,j,k) * dtime  / Thickness%dht(i,j,k,tau))
        if ( Rj .NE. 0.0) then
         if ( volflux .GT. 0.0) then
           Cr = Rjm / Rj
         else
           Cr = Rjp / Rj
         endif
        else
         if ( volflux .GT. 0.0) then
           Cr = Rjm * 1.0e20
         else
           Cr = Rjp * 1.0e20
         endif
        endif
        Cr = max(0.0, max(min(1.0, 2.0 * Cr), min(2.0, Cr)))
        fbt(i,j) = 0.5 * ( volflux                                            &
                 * ( Tracer%field(i,j,k,taum1) + Tracer%field(i,j,kp1,taum1) )&
                 - abs(volflux) * ( 1.0 + (cfl - 1 ) * Cr ) * Rj              &
                  ) * tmask_mdfl(i,j,kp1) * tmask_mdfl(i,j,k)
        tmp_updated_tracer(i,j) = tmp_updated_tracer(i,j)                     &
               + dtime / Thickness%dht(i,j,k,tau) * (                         &
                 Grd%datr(i,j) * ( fbt(i,j) - ftp(i,j))                       &
                 + Tracer%field(i,j,k,taum1) *                                &
                 (wkm1(i,j) - Adv_vel%w_bt(i,j,k)) )
        flux_z(i,j,k) = fbt(i,j)
        ftp(i,j) = fbt(i,j)
      enddo !i
    enddo !j
!
! Update the two-point tracer halo
!
    call mpp_update_domains (tmp_updated_tracer, Dom_mdfl%domain2d)
!
! Calculate flux at the eastern wall of the boxes
!
    do j=jsc,jec
      do i=isc-1,iec
        Rjp = ( tmp_updated_tracer(i+2,j) - tmp_updated_tracer(i+1,j) )       &
                      * tmask_mdfl(i+2,j,k) *       tmask_mdfl(i+1,j,k)
        Rj =  ( tmp_updated_tracer(i+1,j) - tmp_updated_tracer( i ,j) )       &
                      * tmask_mdfl(i+1,j,k) *       tmask_mdfl( i ,j,k)
        Rjm = ( tmp_updated_tracer( i ,j) - tmp_updated_tracer(i-1,j) )       &
                      * tmask_mdfl( i ,j,k) *       tmask_mdfl(i-1,j,k)
        volflux = Grd%dyte(i,j) * Adv_vel%uh_et(i,j,k)
        cfl = abs( Adv_vel%uh_et(i,j,k) * dtime * 2.0                         &
        / ((Thickness%dht(i,j,k,tau) + Thickness%dht(i+1,j,k,tau)) * Grd%dxte(i,j)))
        if ( Rj .NE. 0.0) then
         if ( volflux .GT. 0.0) then
           Cr = Rjm / Rj
         else
           Cr = Rjp / Rj
         endif
        else
         if ( volflux .GT. 0.0) then
           Cr = Rjm * 1.0e20
         else
           Cr = Rjp * 1.0e20
         endif
        endif
        Cr = max(0.0, max(min(1.0, 2.0 * Cr), min(2.0, Cr)))
        fe(i,j) = 0.5 * ( volflux                                             &
                   * (tmp_updated_tracer(i+1,j) + tmp_updated_tracer(i,j) )   &
                 - abs(volflux) * ( 1.0 + (cfl - 1 ) * Cr ) * Rj              &
                  ) * tmask_mdfl(i,j,k) * tmask_mdfl(i+1,j,k)
      enddo !i
    enddo !j
!
! Update the tracer
!
    do j=jsc,jec
      do i=isc,iec
        tmp_updated_tracer(i,j) = tmp_updated_tracer(i,j)                           &
            + dtime * tmask_mdfl(i,j,k) * Grd%datr(i,j) / Thickness%dht(i,j,k,tau)  &
            * (fe(i-1,j) - fe(i,j)                                                  &
            + Tracer%field(i,j,k,taum1)* (                                          &
                      Grd%dyte(i,j) *   Adv_vel%uh_et( i ,j,k)                      &
                    - Grd%dyte(i-1,j) * Adv_vel%uh_et(i-1,j,k) ) )
        flux_x(i,j,k) = fe(i,j)
      enddo !i
    enddo !j
!
! Update the two-point tracer halo
!
    call mpp_update_domains (tmp_updated_tracer, Dom_mdfl%domain2d)
!
! Calculate flux at the northern wall of the boxes
!
    do j=jsc-1,jec
      do i=isc,iec
        Rjp = ( tmp_updated_tracer(i,j+2) - tmp_updated_tracer(i,j+1) )       &
                      * tmask_mdfl(i,j+2,k) *       tmask_mdfl(i,j+1,k)
        Rj =  ( tmp_updated_tracer(i,j+1) - tmp_updated_tracer(i, j ) )       &
                      * tmask_mdfl(i,j+1,k) *       tmask_mdfl(i, j ,k)
        Rjm = ( tmp_updated_tracer(i, j ) - tmp_updated_tracer(i,j-1) )       &
                      * tmask_mdfl(i, j ,k) *       tmask_mdfl(i,j-1,k)
        volflux = Grd%dxtn(i,j) * Adv_vel%vh_nt(i,j,k)
        cfl = abs(Adv_vel%vh_nt(i,j,k) * dtime * 2.0                          &
        / ((Thickness%dht(i,j,k,tau) + Thickness%dht(i,j+1,k,tau)) * Grd%dytn(i,j)))
        if ( Rj .NE. 0.0) then
         if ( volflux .GT. 0.0) then
           Cr = Rjm / Rj
         else
           Cr = Rjp / Rj
         endif
        else
         if ( volflux .GT. 0.0) then
           Cr = Rjm * 1.0e20
         else
           Cr = Rjp * 1.0e20
         endif
        endif
        Cr = max(0.0, max(min(1.0, 2.0 * Cr), min(2.0, Cr)))
        fn(i,j) = 0.5 * ( volflux                                             &
                   * (tmp_updated_tracer(i,j+1) + tmp_updated_tracer(i,j) )   &
                 - abs(volflux) * ( 1.0 + (cfl - 1.0 ) * Cr ) * Rj            &
                  ) * tmask_mdfl(i,j,k) * tmask_mdfl(i,j+1,k)
      enddo !i
    enddo !j
!
! Update the tracer
!
    do j=jsc,jec
      do i=isc,iec
        tmp_updated_tracer(i,j) = tmp_updated_tracer(i,j)                       &
         + dtime * Grd%tmask(i,j,k) * Grd%datr(i,j) / Thickness%dht(i,j,k,tau)  &
         * (fn(i,j-1) - fn(i,j))
        flux_y(i,j,k) = fn(i,j)
      enddo !i
    enddo !j
!
! Calculate the overall tendency
!
    do j=jsc,jec
      do i=isc,iec
        tmp_updated_tracer(i,j) = tmp_updated_tracer(i,j)                     &
           + dtime * Tracer%field(i,j,k,taum1) / Thickness%dht(i,j,k,tau) * ( &
                      Adv_vel%w_bt(i,j,k) - wkm1(i,j)                         &
             + Grd%datr(i,j)*(                                                &
                         ( Grd%dyte(i-1,j) * Adv_vel%uh_et(i-1,j,k)           &
                         - Grd%dyte( i ,j) * Adv_vel%uh_et( i ,j,k))))
!
! By convention, advection is applied as a negative tendency (i.e., convergence)
!
        advect_tracer_mdfl_sub_b(i,j,k) = - Thickness%dht(i,j,k,tau) *     &
          ( tmp_updated_tracer(i,j) - Tracer%field(i,j,k,taum1) )             &
            / dtime * tmask_mdfl(i,j,k)
      enddo !i
    enddo !j
!
! Update vertical velocity for next level
!
    do j=jsc,jec
      do i=isc,iec
        wkm1(i,j) = Adv_vel%w_bt(i,j,k)
      enddo !i
    enddo !j

  enddo !k
  
  call mpp_clock_end(id_clock_mdfl_sup_b)

end function advect_tracer_mdfl_sub_b
! </FUNCTION> NAME="advect_tracer_mdfl_sub_b"



!#######################################################################
! <FUNCTION NAME="advect_tracer_mdfl_sweby">
!
! <DESCRIPTION>
! Compute tendency due to 3D advection of tracers using a multi-dimensional flux-limited
! method.  This method differs from other methods in the following ways:
!
! 1) Horizontal and vertical advection are combined
! 2) Calculations of the three coordinates (Z, X, Y) are performed sequentially as updates
!        to the tracer, so that the advection components for X and Y depend on Z and Z,Y
!        respectively... This helps limit the flux.
! 3) During the update for each direction, the 3rd order Sweby flux limiter is applied.
! 4) Flux divergence is included within the calculation to also help limit the flux:
!          - During the update for each direction, the divergence in each direction is added.
!          - During the overall tendency calculated, the divergence in all three directions
!                is removed.
! 5) All fluxes are functions of the tracer field at the taum1 time step.  This 
!    means that this method is ideally suited for the twolevel time stepping scheme,
!    in which "taum1=tau", thus enabling twice the tracer time step available for the 
!    threelevel scheme. 
!
! The calculation proceeds as follows:
!
! IMPORTANT NOTE: If this scheme is used at all, it must be used as the option for BOTH
! horizontal and vertical advection.  In the the tracer tendency, it is applied as the 
! horizontal term, but applies to vertical as well, for which case the vertical term in
! the tracer tendency equation is set to zero.
!
! This scheme was ported to mom4 from the MIT-GCM by John Dunne and Alistair Adcroft  
! during Summer 2003 
!
!
! Griffies: 5/27/04
! Optimized by filling 3d arrays prior to sending for mpp_update_domains
! </DESCRIPTION>
!
function advect_tracer_mdfl_sweby(Time, Adv_vel, Tracer, Thickness, dtime)

  type(ocean_time_type), intent(in)        :: Time
  type(ocean_adv_vel_type), intent(in)     :: Adv_vel
  type(ocean_prog_tracer_type), intent(in) :: Tracer
  type(ocean_thickness_type), intent(in)   :: Thickness
  real, intent(in)                         :: dtime

  real,dimension(isc:iec,jsc:jec)          :: ftp
  real,dimension(isc:iec,jsc:jec)          :: fbt
  real,dimension(isc:iec,jsc:jec)          :: wkm1

  real,dimension(isc:iec,jsc:jec,nk)       :: advect_tracer_mdfl_sweby 

  integer                                  :: i, j, k
  integer                                  :: kp1, kp2, km1
  integer                                  :: tau, taum1
  real                                     :: Rjm, Rj, Rjp, cfl, volflux
  real                                     :: d0, d1, thetaP, psiP 
  real                                     :: thetaM, psiM

  call mpp_clock_begin(id_clock_mdfl_sweby)

  ftp         = 0.0
  fbt         = 0.0
  wkm1        = 0.0
  tracer_mdfl = 0.0
  flux_x      = 0.0
  flux_y      = 0.0

  tau   = Time%tau
  taum1 = Time%taum1

  do k=1,nk

     do j=jsc,jec
        do i=isc,iec
           tracer_mdfl(i,j,k) = Tracer%field(i,j,k,taum1)
        enddo
     enddo

     kp1 = min(k+1,nk)
     kp2 = min(k+2,nk)
     km1 = max(k-1,1)

     ! calculate flux at bottom of box and update the tracer
     do j=jsc,jec
        do i=isc,iec

           Rjp = (Tracer%field(i,j,km1,taum1) - Tracer%field(i,j,k,taum1))    &
                 *tmask_mdfl(i,j,km1)*tmask_mdfl(i,j,k)
           Rj  = (Tracer%field(i,j,k,taum1) - Tracer%field(i,j,kp1,taum1))    &
                 *tmask_mdfl(i,j,k)*tmask_mdfl(i,j,kp1)
           Rjm = (Tracer%field(i,j,kp1,taum1) - Tracer%field(i,j,kp2,taum1))  &
                 *tmask_mdfl(i,j,kp1)*tmask_mdfl(i,j,kp2)

           volflux = Grd%dat(i,j) * Adv_vel%w_bt(i,j,k)
           cfl = abs(Adv_vel%w_bt(i,j,k) * dtime  / Thickness%dht(i,j,k,tau))
           d0 = (2.0 - cfl) * (1.0 - cfl) * onesixth
           d1 = (1.0 - cfl * cfl) * onesixth

           thetaP = Rjm / (1.0e-30 + Rj)
           thetaM = Rjp / (1.0e-30 + Rj)

           psiP = max(0.0, min(min(1.0, d0 + d1 * thetaP),       &
                (1.0 - cfl) / (1.0e-30 + cfl) * thetaP ) )
           psiM = max(0.0, min(min(1.0, d0 + d1 * thetaM),       &
                (1.0 - cfl) / (1.0e-30 + cfl) * thetaM ) )

           fbt(i,j) =  0.5 * ( ( volflux + abs(volflux) )        &
                * ( Tracer%field(i,j,kp1,taum1) + psiP * Rj )    &
                + ( volflux - abs(volflux) )                     &
                * ( Tracer%field(i,j, k ,taum1) - psiM * Rj ) )  &
                * tmask_mdfl(i,j,kp1) * tmask_mdfl(i,j,k)

           tracer_mdfl(i,j,k) = tracer_mdfl(i,j,k)               &
                + dtime / Thickness%dht(i,j,k,tau) * (           &
                Grd%datr(i,j) * ( fbt(i,j) - ftp(i,j))           &
                + Tracer%field(i,j,k,taum1) *                    &
                (wkm1(i,j) - Adv_vel%w_bt(i,j,k)) )

           flux_z(i,j,k) = fbt(i,j)
           ftp(i,j)      = fbt(i,j)

        enddo
     enddo

     ! update vertical velocity for next k-level
     do j=jsc,jec
        do i=isc,iec
           wkm1(i,j) = Adv_vel%w_bt(i,j,k)
        enddo
     enddo

  enddo ! end of k-loop
  call mpp_update_domains (tracer_mdfl, Dom_mdfl%domain2d, flags=XUPDATE)


  ! calculate flux at the eastern wall of the boxes
  do k=1,nk

     do j=jsc,jec
        do i=isc-1,iec

           Rjp = (tracer_mdfl(i+2,j,k) - tracer_mdfl(i+1,j,k))    &
                 *tmask_mdfl(i+2,j,k)*tmask_mdfl(i+1,j,k)
           Rj  = (tracer_mdfl(i+1,j,k) - tracer_mdfl( i ,j,k) )   &
                 *tmask_mdfl(i+1,j,k)*tmask_mdfl(i,j,k)
           Rjm = (tracer_mdfl(i,j,k) - tracer_mdfl(i-1,j,k))      &
                 *tmask_mdfl(i,j,k)*tmask_mdfl(i-1,j,k)

           volflux = Grd%dyte(i,j) * Adv_vel%uh_et(i,j,k)
           cfl = abs(Adv_vel%uh_et(i,j,k) * dtime * 2.0           &
           / ((Thickness%dht(i,j,k,tau) + Thickness%dht(i+1,j,k,tau)) * Grd%dxte(i,j)))
           d0 = (2.0 - cfl) * (1.0 - cfl) * onesixth
           d1 = (1.0 - cfl * cfl) * onesixth

           thetaP = Rjm / (1.0e-30 + Rj)
           thetaM = Rjp / (1.0e-30 + Rj)

           psiP = max(0.0, min(min(1.0, d0 + d1 * thetaP),      &
                (1.0 - cfl) / (1.0e-30 + cfl) * thetaP ) )
           psiM = max(0.0, min(min(1.0, d0 + d1 * thetaM),      &
                (1.0 - cfl) / (1.0e-30 + cfl) * thetaM ) )

           flux_x(i,j,k) =  0.5 * ( ( volflux + abs(volflux) )  &
                * ( tracer_mdfl( i ,j,k) + psiP * Rj )          &
                + ( volflux - abs(volflux) )                    &
                * ( tracer_mdfl(i+1,j,k) - psiM * Rj ) )        &
                * tmask_mdfl(i,j,k) * tmask_mdfl(i+1,j,k)
        enddo
     enddo

     ! update the tracer
     do j=jsc,jec
        do i=isc,iec
           tracer_mdfl(i,j,k) = tracer_mdfl(i,j,k)                                &
          + dtime * tmask_mdfl(i,j,k) * Grd%datr(i,j) / Thickness%dht(i,j,k,tau)  &
                * (flux_x(i-1,j,k) - flux_x(i,j,k)                                &
                + Tracer%field(i,j,k,taum1)* (                                    &
                Grd%dyte(i,j) *   Adv_vel%uh_et( i ,j,k)                          &
                - Grd%dyte(i-1,j) * Adv_vel%uh_et(i-1,j,k) ) )
        enddo
     enddo

  enddo ! end of k-loop
  call mpp_update_domains (tracer_mdfl, Dom_mdfl%domain2d, flags=YUPDATE)


  ! calculate flux at the northern wall of the boxes
  do j=jsc,jec
     do i=isc,iec
        wkm1(i,j) = 0.0
     enddo
  enddo

  do k=1,nk

     do j=jsc-1,jec
        do i=isc,iec

           Rjp = (tracer_mdfl(i,j+2,k) - tracer_mdfl(i,j+1,k))       &
                 *tmask_mdfl(i,j+2,k)*tmask_mdfl(i,j+1,k)
           Rj  = (tracer_mdfl(i,j+1,k) - tracer_mdfl(i,j,k))         &
                 *tmask_mdfl(i,j+1,k)*tmask_mdfl(i,j,k)
           Rjm = (tracer_mdfl(i,j,k) - tracer_mdfl(i,j-1,k))         &
                *tmask_mdfl(i,j,k)*tmask_mdfl(i,j-1,k)

           volflux = Grd%dxtn(i,j) * Adv_vel%vh_nt(i,j,k)
           cfl = abs(Adv_vel%vh_nt(i,j,k) * dtime * 2.0              &
           / ((Thickness%dht(i,j,k,tau) + Thickness%dht(i,j+1,k,tau)) * Grd%dytn(i,j)))
           d0 = (2.0 - cfl) * (1.0 - cfl) * onesixth
           d1 = (1.0 - cfl * cfl) * onesixth

           thetaP = Rjm / (1.0e-30 + Rj)
           thetaM = Rjp / (1.0e-30 + Rj)

           psiP = max(0.0, min(min(1.0, d0 + d1 * thetaP),           &
                (1.0 - cfl) / (1.0e-30 + cfl) * thetaP ) )
           psiM = max(0.0, min(min(1.0, d0 + d1 * thetaM),           &
                (1.0 - cfl) / (1.0e-30 + cfl) * thetaM ) )

           flux_y(i,j,k) =  0.5 * ( ( volflux + abs(volflux) )       &
                * ( tracer_mdfl(i,j,k) + psiP * Rj )                 &
                + ( volflux - abs(volflux) )                         &
                * ( tracer_mdfl(i,j+1,k) - psiM * Rj ) )             &
                * tmask_mdfl(i,j,k) * tmask_mdfl(i,j+1,k)

        enddo
     enddo

     ! update tracer
     do j=jsc,jec
        do i=isc,iec
           tracer_mdfl(i,j,k) = tracer_mdfl(i,j,k)                                      &
                + dtime * tmask_mdfl(i,j,k) * Grd%datr(i,j) / Thickness%dht(i,j,k,tau)  &
                * (flux_y(i,j-1,k) - flux_y(i,j,k))
        enddo
     enddo

    ! calculate the overall tendency
    do j=jsc,jec
      do i=isc,iec

        tracer_mdfl(i,j,k) = tracer_mdfl(i,j,k)                               &
           + dtime * Tracer%field(i,j,k,taum1) / Thickness%dht(i,j,k,tau) * ( &
                      Adv_vel%w_bt(i,j,k) - wkm1(i,j)                         &
             + Grd%datr(i,j)*(                                                &
                         ( Grd%dyte(i-1,j) * Adv_vel%uh_et(i-1,j,k)           &
                         - Grd%dyte( i ,j) * Adv_vel%uh_et( i ,j,k))))

        ! advection is applied as a convergence
        advect_tracer_mdfl_sweby(i,j,k) = -Thickness%dht(i,j,k,tau) *         &
          ( tracer_mdfl(i,j,k) - Tracer%field(i,j,k,taum1) )                  &
            / dtime * tmask_mdfl(i,j,k)
      enddo 
    enddo 

   ! update vertical velocity for next level
    do j=jsc,jec
      do i=isc,iec
        wkm1(i,j) = Adv_vel%w_bt(i,j,k)
      enddo 
    enddo 

  enddo ! end of k-loop

  
  call mpp_clock_end(id_clock_mdfl_sweby)

end function advect_tracer_mdfl_sweby
! </FUNCTION> NAME="advect_tracer_mdfl_sweby"



!#######################################################################
! <SUBROUTINE NAME="advect_tracer_sweby_all">
!
! <DESCRIPTION>
!
! Griffies: June 2004
!
! Sweby scheme optimized by doing all tracers at once, so 
! can send larger packet to mpp_update_domains.
!
! This scheme is available ONLY when advecting all tracers with the 
! sweby scheme (which is the common case at GFDL).
!
! This scheme has not been ported to the OBC, though such port 
! is likely trivial to do.  
!
! </DESCRIPTION>
!
subroutine advect_tracer_sweby_all(Time, Adv_vel, T_prog, Thickness, dtime)

  type(ocean_time_type), intent(in)           :: Time
  type(ocean_adv_vel_type), intent(in)        :: Adv_vel
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)
  type(ocean_thickness_type), intent(in)      :: Thickness
  real, intent(in)                            :: dtime

  real,dimension(isc:iec,jsc:jec)             :: ftp
  real,dimension(isc:iec,jsc:jec)             :: fbt
  real,dimension(isc:iec,jsc:jec)             :: wkm1
  real,dimension(isd:ied,jsd:jed)             :: tmp_flux

  integer                                     :: i, j, k, n
  integer                                     :: kp1, kp2, km1
  integer                                     :: tau, taum1
  real                                        :: Rjm, Rj, Rjp, cfl, volflux
  real                                        :: d0, d1, thetaP, psiP 
  real                                        :: thetaM, psiM

  call mpp_clock_begin(id_clock_mdfl_sweby_all)

  tau     = Time%tau
  taum1   = Time%taum1
  flux_x  = 0.0
  flux_y  = 0.0
  flux_z  = 0.0

  ! calculate flux at bottom face of the T-cells
  do n=1,num_prog_tracers

     ftp  = 0.0
     fbt  = 0.0
     wkm1 = 0.0

     do k=1,nk
        do j=jsd,jed
           do i=isd,ied
              T_prog(n)%wrk1(i,j,k) = 0.0
           enddo
        enddo
     enddo

     do k=1,nk

        do j=jsc,jec
           do i=isc,iec
              tracer_mdfl_all(n)%field(i,j,k) = T_prog(n)%field(i,j,k,taum1)
           enddo
        enddo

        kp1 = min(k+1,nk)
        kp2 = min(k+2,nk)
        km1 = max(k-1,1)

        do j=jsc,jec
           do i=isc,iec

              Rjp = (T_prog(n)%field(i,j,km1,taum1) - T_prog(n)%field(i,j,k,taum1))    &
                   *tmask_mdfl(i,j,km1)*tmask_mdfl(i,j,k)
              Rj  = (T_prog(n)%field(i,j,k,taum1) - T_prog(n)%field(i,j,kp1,taum1))    &
                   *tmask_mdfl(i,j,k)*tmask_mdfl(i,j,kp1)
              Rjm = (T_prog(n)%field(i,j,kp1,taum1) - T_prog(n)%field(i,j,kp2,taum1))  &
                   *tmask_mdfl(i,j,kp1)*tmask_mdfl(i,j,kp2)

              volflux = Grd%dat(i,j) * Adv_vel%w_bt(i,j,k)
              cfl = abs(Adv_vel%w_bt(i,j,k) * dtime  / Thickness%dht(i,j,k,tau))
              d0 = (2.0 - cfl) * (1.0 - cfl) * onesixth
              d1 = (1.0 - cfl * cfl) * onesixth

              thetaP = Rjm / (1.0e-30 + Rj)
              thetaM = Rjp / (1.0e-30 + Rj)

              psiP = max(0.0, min(min(1.0, d0 + d1 * thetaP),       &
                   (1.0 - cfl) / (1.0e-30 + cfl) * thetaP ) )
              psiM = max(0.0, min(min(1.0, d0 + d1 * thetaM),       &
                   (1.0 - cfl) / (1.0e-30 + cfl) * thetaM ) )

              fbt(i,j) =  0.5 * ( ( volflux + abs(volflux) )           &
                   * ( T_prog(n)%field(i,j,kp1,taum1) + psiP * Rj )    &
                   + ( volflux - abs(volflux) )                        &
                   * ( T_prog(n)%field(i,j, k ,taum1) - psiM * Rj ) )  &
                   * tmask_mdfl(i,j,kp1) * tmask_mdfl(i,j,k)

              tracer_mdfl_all(n)%field(i,j,k) = tracer_mdfl_all(n)%field(i,j,k)  &
                   + dtime / Thickness%dht(i,j,k,tau) * (                        &
                   Grd%datr(i,j) * ( fbt(i,j) - ftp(i,j))                        &
                   + T_prog(n)%field(i,j,k,taum1) *                              &
                   (wkm1(i,j) - Adv_vel%w_bt(i,j,k)) )

              flux_z(i,j,k) = fbt(i,j)
              ftp(i,j)      = fbt(i,j)

           enddo
        enddo

        ! update vertical velocity for next k-level
        do j=jsc,jec
           do i=isc,iec
              wkm1(i,j) = Adv_vel%w_bt(i,j,k)
           enddo
        enddo

     enddo ! end of k-loop

     call mpp_update_domains (tracer_mdfl_all(n)%field(:,:,:), Dom_mdfl%domain2d, flags=XUPDATE, complete=T_prog(n)%complete)

     if (id_adv_flux_z(n) > 0) used = send_data(id_adv_flux_z(n), T_prog(n)%conversion*flux_z(:,:,:), &
         Time%model_time, rmask=Grd%tmask(:,:,:), &
         is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

  enddo ! end of n-loop for tracers  


  ! calculate flux at the eastern face of the T-cells
  do n=1,num_prog_tracers

     do k=1,nk

        do j=jsc,jec
           do i=isc-1,iec

              Rjp = (tracer_mdfl_all(n)%field(i+2,j,k) - tracer_mdfl_all(n)%field(i+1,j,k))    &
                   *tmask_mdfl(i+2,j,k)*tmask_mdfl(i+1,j,k)
              Rj  = (tracer_mdfl_all(n)%field(i+1,j,k) - tracer_mdfl_all(n)%field( i ,j,k) )   &
                   *tmask_mdfl(i+1,j,k)*tmask_mdfl(i,j,k)
              Rjm = (tracer_mdfl_all(n)%field(i,j,k) - tracer_mdfl_all(n)%field(i-1,j,k))      &
                   *tmask_mdfl(i,j,k)*tmask_mdfl(i-1,j,k)

              volflux = Grd%dyte(i,j) * Adv_vel%uh_et(i,j,k)
              cfl = abs(Adv_vel%uh_et(i,j,k) * dtime * 2.0    &
                   / ((Thickness%dht(i,j,k,tau) + Thickness%dht(i+1,j,k,tau)) * Grd%dxte(i,j)))
              d0 = (2.0 - cfl) * (1.0 - cfl) * onesixth
              d1 = (1.0 - cfl * cfl) * onesixth

              thetaP = Rjm / (1.0e-30 + Rj)
              thetaM = Rjp / (1.0e-30 + Rj)

              psiP = max(0.0, min(min(1.0, d0 + d1 * thetaP),  &
                   (1.0 - cfl) / (1.0e-30 + cfl) * thetaP ) )
              psiM = max(0.0, min(min(1.0, d0 + d1 * thetaM),  &
                   (1.0 - cfl) / (1.0e-30 + cfl) * thetaM ) )

              flux_x(i,j,k) =  0.5 * ( ( volflux + abs(volflux) )        &
                   * ( tracer_mdfl_all(n)%field( i ,j,k) + psiP * Rj )   &
                   + ( volflux - abs(volflux) )                          &
                   * ( tracer_mdfl_all(n)%field(i+1,j,k) - psiM * Rj ) ) &
                   * tmask_mdfl(i,j,k) * tmask_mdfl(i+1,j,k)
           enddo
        enddo

        ! update the tracer
        do j=jsc,jec
           do i=isc,iec
              tracer_mdfl_all(n)%field(i,j,k) = tracer_mdfl_all(n)%field(i,j,k)            &
                   + dtime * tmask_mdfl(i,j,k) * Grd%datr(i,j) / Thickness%dht(i,j,k,tau)  &
                   * (flux_x(i-1,j,k) - flux_x(i,j,k)                                      &
                   + T_prog(n)%field(i,j,k,taum1)* (                                       &
                     Grd%dyte(i,j) *   Adv_vel%uh_et( i ,j,k)                              &
                   - Grd%dyte(i-1,j) * Adv_vel%uh_et(i-1,j,k) ) )
           enddo
        enddo

     enddo ! end of k-loop

     call mpp_update_domains (tracer_mdfl_all(n)%field(:,:,:), Dom_mdfl%domain2d, flags=YUPDATE, complete=T_prog(n)%complete)

     if (id_adv_flux_x(n) > 0) used = send_data(id_adv_flux_x(n), T_prog(n)%conversion*flux_x(:,:,:), &
          Time%model_time, rmask=Grd%tmask(:,:,:), &
          is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

     if (id_adv_flux_x_int_z(n) > 0) then 
         tmp_flux(:,:) = 0.0
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  tmp_flux(i,j) = tmp_flux(i,j) +  flux_x(i,j,k) 
               enddo
            enddo
         enddo
         used = send_data(id_adv_flux_x_int_z(n), T_prog(n)%conversion*tmp_flux(:,:), &
              Time%model_time, rmask=Grd%tmask(:,:,1), &
              is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     endif

  enddo ! end of n-loop for tracers  


  ! calculate flux at the northern face of the T-cells 
  do n=1,num_prog_tracers 

     do j=jsc,jec
        do i=isc,iec
           wkm1(i,j) = 0.0
        enddo
     enddo

     do k=1,nk

        do j=jsc-1,jec
           do i=isc,iec

              Rjp = (tracer_mdfl_all(n)%field(i,j+2,k) - tracer_mdfl_all(n)%field(i,j+1,k))   &
                   *tmask_mdfl(i,j+2,k)*tmask_mdfl(i,j+1,k)
              Rj  = (tracer_mdfl_all(n)%field(i,j+1,k) - tracer_mdfl_all(n)%field(i,j,k))     &
                   *tmask_mdfl(i,j+1,k)*tmask_mdfl(i,j,k)
              Rjm = (tracer_mdfl_all(n)%field(i,j,k)   - tracer_mdfl_all(n)%field(i,j-1,k))   &
                   *tmask_mdfl(i,j,k)*tmask_mdfl(i,j-1,k)

              volflux = Grd%dxtn(i,j) * Adv_vel%vh_nt(i,j,k)
              cfl = abs(Adv_vel%vh_nt(i,j,k) * dtime * 2.0              &
                   / ((Thickness%dht(i,j,k,tau) + Thickness%dht(i,j+1,k,tau)) * Grd%dytn(i,j)))
              d0 = (2.0 - cfl) * (1.0 - cfl) * onesixth
              d1 = (1.0 - cfl * cfl) * onesixth

              thetaP = Rjm / (1.0e-30 + Rj)
              thetaM = Rjp / (1.0e-30 + Rj)

              psiP = max(0.0, min(min(1.0, d0 + d1 * thetaP),   &
                   (1.0 - cfl) / (1.0e-30 + cfl) * thetaP ) )
              psiM = max(0.0, min(min(1.0, d0 + d1 * thetaM),   &
                   (1.0 - cfl) / (1.0e-30 + cfl) * thetaM ) )

              flux_y(i,j,k) =  0.5 * ( ( volflux + abs(volflux) )         &
                   * ( tracer_mdfl_all(n)%field(i,j,k) + psiP * Rj )      &
                   + ( volflux - abs(volflux) )                           &
                   * ( tracer_mdfl_all(n)%field(i,j+1,k) - psiM * Rj ) )  &
                   * tmask_mdfl(i,j,k) * tmask_mdfl(i,j+1,k)

           enddo
        enddo

        ! update tracer
        do j=jsc,jec
           do i=isc,iec
              tracer_mdfl_all(n)%field(i,j,k) = tracer_mdfl_all(n)%field(i,j,k)            &
                   + dtime * tmask_mdfl(i,j,k) * Grd%datr(i,j) / Thickness%dht(i,j,k,tau)  &
                   * (flux_y(i,j-1,k) - flux_y(i,j,k))
           enddo
        enddo

        ! calculate the overall tendency
        do j=jsc,jec
           do i=isc,iec

              tracer_mdfl_all(n)%field(i,j,k) = tracer_mdfl_all(n)%field(i,j,k)          &
                   + dtime * T_prog(n)%field(i,j,k,taum1) / Thickness%dht(i,j,k,tau) * ( &
                   Adv_vel%w_bt(i,j,k) - wkm1(i,j)                                       &
                   + Grd%datr(i,j)*(                                                     &
                   ( Grd%dyte(i-1,j) * Adv_vel%uh_et(i-1,j,k)                            &
                   - Grd%dyte( i ,j) * Adv_vel%uh_et( i ,j,k))))

              T_prog(n)%wrk1(i,j,k) = &
                    Thickness%dht(i,j,k,tau)*(tracer_mdfl_all(n)%field(i,j,k)-T_prog(n)%field(i,j,k,taum1))/dtime  &
                   *tmask_mdfl(i,j,k) 
              T_prog(n)%th_tendency(i,j,k) = T_prog(n)%th_tendency(i,j,k) + T_prog(n)%wrk1(i,j,k)

           enddo
        enddo

        ! update vertical velocity for next level
        do j=jsc,jec
           do i=isc,iec
              wkm1(i,j) = Adv_vel%w_bt(i,j,k)
           enddo
        enddo

     enddo ! end of k-loop

     if (id_adv_flux_y(n) > 0) used = send_data(id_adv_flux_y(n), T_prog(n)%conversion*flux_y(:,:,:), &
          Time%model_time, rmask=Grd%tmask(:,:,:), &
          is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

     if (id_adv_flux_y_int_z(n) > 0) then 
         tmp_flux(:,:) = 0.0
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  tmp_flux(i,j) = tmp_flux(i,j) +  flux_y(i,j,k) 
               enddo
            enddo
         enddo
         used = send_data(id_adv_flux_y_int_z(n), T_prog(n)%conversion*tmp_flux(:,:), &
              Time%model_time, rmask=Grd%tmask(:,:,1), &
              is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
     endif

     if (id_sweby_advect(n) > 0) used = send_data(id_sweby_advect(n), T_prog(n)%conversion*T_prog(n)%wrk1(:,:,:), &
          Time%model_time, rmask=Grd%tmask(:,:,:), &
          is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)


  enddo ! end of n-loop
  
  call mpp_clock_end(id_clock_mdfl_sweby_all)


end subroutine advect_tracer_sweby_all
! </SUBROUTINE> NAME="advect_tracer_sweby_all"




end module ocean_tracer_advect_mod
