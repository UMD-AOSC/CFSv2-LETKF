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
! Modified: Xingren Wu
!           Xingren.Wu@noaa.gov
module ocean_grids_mod
!
!<CONTACT EMAIL="Matthew.Harrison@noaa.gov"> M.J. Harrison 
!</CONTACT>
!
!<CONTACT EMAIL="Ronald.Pacanowski@noaa.gov"> R. C. Pacanowski 
!</CONTACT>
!
!<CONTACT EMAIL="Zhi.Liang@noaa.gov"> Zhi Liang 
!</CONTACT>
!
! <REVIEWER EMAIL="Stephen.Griffies@noaa.gov">
! S.M. Griffies 
! </REVIEWER>
!
!<OVERVIEW>
! Set up the ocean model grid spacing 
!</OVERVIEW>
!
!<DESCRIPTION>
! This module sets up the ocean model grid based on information read in 
! from the grid_spec.nc file. It translates the generic names from the 
! grid_spec.nc file to the names used by mom4. 
!</DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! S.M. Griffies, M.J. Harrison, A. Rosati, and R.C. Pacanowski 
! A Technical Guide to MOM4 (2003)
! </REFERENCE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_grids_nml">
!  <DATA NAME="debug_grid" TYPE="logical">
!  For debugging. Note that most of the debugging stuff 
!  has been removed, but keep flag around in case need in future.
!  </DATA> 
!  <DATA NAME="verbose_init" TYPE="logical">
!  Prints out lots of initial checksums.  Useful to have on, so 
!  defaulted to true. 
!   </DATA> 
! </NAMELIST>
!
use constants_mod,     only: pi, radian, radius, epsln
use diag_manager_mod,  only: diag_manager_init, diag_axis_init
use diag_manager_mod,  only: register_static_field, register_diag_field, send_data
use fms_mod,           only: write_version_number, error_mesg, FATAL, NOTE
use fms_mod,           only: open_namelist_file, close_file, check_nml_error
use mpp_domains_mod,   only: mpp_update_domains, mpp_global_field, domain2d
use mpp_domains_mod,   only: mpp_global_sum, BITWISE_EXACT_SUM
use mpp_io_mod,        only: axistype, fieldtype, atttype
use mpp_io_mod,        only: MPP_RDONLY, MPP_NETCDF, MPP_MULTI, MPP_SINGLE 
use mpp_io_mod,        only: mpp_open, mpp_read, mpp_get_atts, mpp_get_axis_data 
use mpp_io_mod,        only: mpp_get_info, mpp_get_axes, mpp_close, mpp_get_fields
use mpp_mod,           only: stdout, stdlog, mpp_error, mpp_chksum, mpp_sum

use ocean_domains_mod, only: get_local_indices, get_global_indices, get_halo_sizes, get_domain_offsets
use ocean_types_mod,   only: ocean_grid_type, ocean_time_type, ocean_domain_type

implicit none

private

#include <ocean_memory.h>

#ifdef STATIC_MEMORY
integer, parameter :: xhalo=1
integer, parameter :: yhalo=1
#else
integer, private   :: xhalo
integer, private   :: yhalo
#endif

! for diagnostics 
integer :: id_ht=-1
integer :: id_hu=-1
integer :: id_dht_dx=-1
integer :: id_dht_dy=-1
logical :: used

logical :: debug_grid   =.false.
logical :: verbose_init =.true.
logical :: beta_plane   = .false.
logical :: f_plane      = .false.

character(len=128) :: version=&
     '$Id$'
character (len=128) :: tagname = &
     '$Name$'

public ocean_grids_init
public set_ocean_grid_size
public set_ocean_hgrid_arrays
public set_ocean_vgrid_arrays
public update_boundaries

private axes_info

logical :: module_is_initialized = .FALSE.
logical :: is_new_grid

  namelist /ocean_grids_nml/ debug_grid, verbose_init, beta_plane, f_plane

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_grids_init">
!
! <DESCRIPTION>
! Initialize the grids module. 
! </DESCRIPTION>
!
subroutine ocean_grids_init(debug)

  logical, intent(in), optional    :: debug
  integer :: ioun, io_status, ierr

  if ( module_is_initialized ) then
    call mpp_error(FATAL, '==>Error from ocean_grids_mod (ocean_grids_init): module has been initialized')
  endif 

  module_is_initialized = .TRUE.

  call write_version_number( version, tagname )

  if (PRESENT(debug)) debug_grid = debug

  ! provide for namelist over-ride of defaults 
  ioun = open_namelist_file()
  read (ioun,ocean_grids_nml,IOSTAT=io_status)
  write (stdout(),'(/)')
  write (stdout(),ocean_grids_nml)  
  write (stdlog(),ocean_grids_nml)
  ierr = check_nml_error(io_status, 'ocean_grids_nml')
  call close_file(ioun)

end subroutine ocean_grids_init
! </SUBROUTINE> NAME="ocean_grids_init"


!#######################################################################
! <SUBROUTINE NAME="set_ocean_grid_size">
!
! <DESCRIPTION>
! Set the ocean grid size.  Model expects the grid specification file
! to be called grid_spec.nc.  
! </DESCRIPTION>
!
subroutine set_ocean_grid_size(Grid, grid_file, grid_name)

  type(ocean_grid_type), intent(inout)      :: Grid
  character(len=*), intent(in), optional    :: grid_file
  character(len=*), intent(in), optional    :: grid_name

  character(len=128)                        :: name
  character(len=128)                        :: grd_file
  integer                                   :: i, unit, ndim, nvar, natt, ntime, len
  type(axistype), dimension(:), allocatable :: axes
  type(atttype), dimension(:), allocatable  :: global_atts
  
  Grid%ni=0 ; Grid%nj=0 ; Grid%nk=0 

  grd_file = 'INPUT/grid_spec'
  if (PRESENT(grid_file)) grd_file = trim(grid_file) 

  call mpp_open(unit,trim(grd_file),MPP_RDONLY,MPP_NETCDF,threading=MPP_MULTI,fileset=MPP_SINGLE)
  call mpp_get_info(unit, ndim, nvar, natt, ntime)
  allocate(axes(ndim))
  call mpp_get_axes(unit,axes)

  !--- determine if it is new_grid or old grid
  do i=1, ndim
     call mpp_get_atts(axes(i),name=name)  
        select case(trim(name))
        case ('grid_x_T')
          is_new_grid = .true.
        case ('gridlon_t')
          is_new_grid = .false.
        end select
  enddo        

  if(is_new_grid) then
    call mpp_error(NOTE, '==>Note from ocean_grids_mod (set_ocean_grid_size): running with is_new_grid=.true.')
  else
    call mpp_error(NOTE, '==>Note from ocean_grids_mod (set_ocean_grid_size): running with is_new_grid=.false.')
  endif  

  do i=1, ndim
     call mpp_get_atts(axes(i),name=name,len=len)       
     if(is_new_grid) then
        select case(trim(name))
        case ('grid_y_T')
           Grid%nj = len
        case ('grid_x_T')
           Grid%ni = len
        case ( 'zt' )
           Grid%nk = len
        end select
     else
        select case(trim(name))
        case ('gridlat_t')
           Grid%nj = len
        case ('gridlon_t')
           Grid%ni = len
        case ( 'zt' )
           Grid%nk = len
        end select
     endif
  enddo

  Grid%beta_plane = beta_plane
  Grid%f_plane = f_plane
  
  if (Grid%ni == 0 .or. Grid%nj == 0 .or. Grid%nk == 0) then
     write(stdout(),*) '==>Error reading grid information from ',trim(grid_file),'. Make sure file exists'
     call mpp_error(FATAL,'==>Error reading grid information from grid file.  Are you sure file exists?')
  endif

  Grid%cyclic=.false.;Grid%solid_walls=.false.;Grid%tripolar=.false.

  allocate(global_atts(natt))
  call mpp_get_atts(unit,global_atts)
  do i=1,natt
     select case (trim(global_atts(i)%name))
     case ('x_boundary_type')
         if (trim(global_atts(i)%catt) == 'cyclic') then
             Grid%cyclic = .true.
             call mpp_error(NOTE,'==>Note from ocean_grids_mod(set_ocean_grid_size): x_boundary_type is cyclic')
         else if (trim(global_atts(i)%catt) == 'solid_walls') then
             Grid%solid_walls = .true.
             call mpp_error(NOTE,'==>Note from ocean_grids_mod(set_ocean_grid_size): x_boundary_type is solid_walls')
         endif
     case ('y_boundary_type')
         if (trim(global_atts(i)%catt) == 'fold_north_edge') then
             Grid%tripolar = .true.
             call mpp_error(NOTE,'==>Note from ocean_grids_mod(set_ocean_grid_size): y_boundary_type is fold_north_edge')
         else
             call mpp_error(NOTE,'==>Note from ocean_grids_mod(set_ocean_grid_size): y_boundary_type is solid_walls')
         endif 
     end select
  end do

  if (PRESENT(grid_name)) then
     Grid%name = grid_name
  else
     Grid%name = 'ocean'
  endif

  deallocate(axes, global_atts)
  call mpp_close(unit)

  if(Grid%tripolar) then 
    write (stdout(),'(1x,/a)')  ' ==> Note: Energy conversion errors are nontrivial when using tripolar=.true.'
    write (stdout(),'(7x,a)')   'The cause is related to the need to update redundantly computed information' 
    write (stdout(),'(7x,a)')   'across the Arctic bipolar fold in a bit-wise exact manner for terms contributing'
    write (stdout(),'(7x,a)')   'to the energy conversion analysis.  The extra code and mpp calls have not been'
    write (stdout(),'(7x,a)')   'implemented. '
  endif 

end subroutine set_ocean_grid_size
! </SUBROUTINE> NAME="set_ocean_grid_size"


!#######################################################################
! <SUBROUTINE NAME="set_ocean_hgrid_arrays">
!
! <DESCRIPTION>
! Define horizontal (and some vertical) grid arrays.
! </DESCRIPTION>
subroutine set_ocean_hgrid_arrays(Time, Domain, Grid)

  type(ocean_time_type), intent(in)          :: Time
  type(ocean_domain_type), intent(inout)     :: Domain
  type(ocean_grid_type), intent(inout)       :: Grid

  real, dimension(:), allocatable            :: data
  real, dimension(:,:), allocatable          :: tmp
  real, dimension(:,:), allocatable          :: sin_tmp
  real, dimension(:,:), allocatable          :: cos_tmp
  real, dimension(:,:), allocatable          :: tmp_local
  real, dimension(:,:), allocatable          :: tmp1_local
  real, dimension(:,:), allocatable          :: tmp2_local
  character(len=128)                         :: name
  character(len=128)                         :: name2
  type(axistype), dimension(:), allocatable  :: axes
  type(fieldtype), dimension(:), allocatable :: fields
  
  real    :: sum, angle, lon_scale, del
  integer :: i, j, k, ii, kp1, kp2, km1, io, len
  integer :: unit, ndim, nvar, natt, ntime, ioff, joff
  logical :: found_zt, found_zw, found_grid_x_t, found_grid_y_t, found_grid_x_u, found_grid_y_u
  logical :: found_xt, found_yt, found_xu, found_yu, found_dxt, found_dyt, found_dxu, found_dyu
  logical :: found_dxtn, found_dytn, found_dxte, found_dyte, found_dxun, found_dyun
  logical :: found_dxue, found_dyue, found_duw, found_due, found_dus, found_dun, found_dtw
  logical :: found_dte, found_dts, found_dtn, found_angle, found_sin_rot, found_cos_rot

  ! set the grid points coordinates (degrees) and grid spacing (degrees)

#ifndef STATIC_MEMORY

  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  call get_global_indices(Domain, isg, ieg, jsg, jeg)
  call get_halo_sizes(Domain, xhalo, yhalo)

  ni = Grid%ni;nj = Grid%nj; nk = Grid%nk
 
  allocate (Grid%xt(isd:ied,jsd:jed))
  allocate (Grid%yt(isd:ied,jsd:jed))
  allocate (Grid%xu(isd:ied,jsd:jed))
  allocate (Grid%yu(isd:ied,jsd:jed))
  allocate (Grid%grid_x_t(ni))
  allocate (Grid%grid_y_t(nj))
  allocate (Grid%grid_x_u(ni))
  allocate (Grid%grid_y_u(nj))
  allocate (Grid%zt(nk))
  allocate (Grid%zw(nk))

  allocate(Grid%phit(isd:ied,jsd:jed))
  allocate(Grid%phiu(isd:ied,jsd:jed))
  allocate(Grid%h1t(isd:ied,jsd:jed))
  allocate(Grid%h2t(isd:ied,jsd:jed))
  allocate(Grid%h1u(isd:ied,jsd:jed))
  allocate(Grid%h2u(isd:ied,jsd:jed))

  allocate (Grid%dxt(isd:ied,jsd:jed))
  allocate (Grid%dxu(isd:ied,jsd:jed))
  allocate (Grid%dyt(isd:ied,jsd:jed))
  allocate (Grid%dyu(isd:ied,jsd:jed))
  allocate (Grid%dat(isd:ied,jsd:jed))
  allocate (Grid%dau(isd:ied,jsd:jed))

  allocate (Grid%dxtn(isd:ied,jsd:jed))
  allocate (Grid%dytn(isd:ied,jsd:jed))
  allocate (Grid%dxte(isd:ied,jsd:jed))
  allocate (Grid%dyte(isd:ied,jsd:jed))
  allocate (Grid%dxun(isd:ied,jsd:jed))
  allocate (Grid%dyun(isd:ied,jsd:jed))
  allocate (Grid%dxue(isd:ied,jsd:jed))
  allocate (Grid%dyue(isd:ied,jsd:jed))

  allocate (Grid%dxtr(isd:ied,jsd:jed))
  allocate (Grid%dxur(isd:ied,jsd:jed))
  allocate (Grid%dytr(isd:ied,jsd:jed))
  allocate (Grid%dyur(isd:ied,jsd:jed))
  allocate (Grid%dxuer(isd:ied,jsd:jed))
  allocate (Grid%dyuer(isd:ied,jsd:jed))
  allocate (Grid%dxunr(isd:ied,jsd:jed))
  allocate (Grid%dyunr(isd:ied,jsd:jed))
  allocate (Grid%dxter(isd:ied,jsd:jed))
  allocate (Grid%dytnr(isd:ied,jsd:jed))
  allocate (Grid%dyue_dxuer(isd:ied,jsd:jed))
  allocate (Grid%dxun_dyunr(isd:ied,jsd:jed))
  allocate (Grid%datr(isd:ied,jsd:jed))
  allocate (Grid%daur(isd:ied,jsd:jed))
  allocate (Grid%dater(isd:ied,jsd:jed))
  allocate (Grid%datnr(isd:ied,jsd:jed))

  allocate (Grid%dh1dy(isd:ied,jsd:jed))
  allocate (Grid%dh2dx(isd:ied,jsd:jed))

  allocate ( Grid%duw(isd:ied,jsd:jed) )
  allocate ( Grid%due(isd:ied,jsd:jed) )
  allocate ( Grid%dus(isd:ied,jsd:jed) )
  allocate ( Grid%dun(isd:ied,jsd:jed) )
  allocate ( Grid%dtw(isd:ied,jsd:jed) )
  allocate ( Grid%dte(isd:ied,jsd:jed) )
  allocate ( Grid%dts(isd:ied,jsd:jed) )
  allocate ( Grid%dtn(isd:ied,jsd:jed) )

  allocate(Grid%sin_rot(isd:ied,jsd:jed))
  allocate(Grid%cos_rot(isd:ied,jsd:jed))

  allocate (Grid%obc_tmask(isd:ied,jsd:jed))
  allocate (Grid%obc_umask(isd:ied,jsd:jed))

#endif

  call get_domain_offsets(Domain,ioff, joff)  
  
  call mpp_open(unit,'INPUT/grid_spec',MPP_RDONLY,MPP_NETCDF,threading=MPP_MULTI,fileset=MPP_SINGLE)
   print *,'in mom4,ocean_grids.F90, before open INPUT/grid_spec,unit=',unit
  call mpp_get_info(unit, ndim, nvar, natt, ntime)
  allocate(axes(ndim))
  call mpp_get_axes(unit,axes)
  allocate(fields(nvar))
  call mpp_get_fields(unit, fields)

   print *,'in mom4,ocean_grids.F90, before open get field'

  Grid%zt=0.0 ; Grid%zw=0.0 ; Grid%grid_x_t=0.0 ; Grid%grid_y_t=0.0 ; Grid%grid_x_u=0.0 ; Grid%grid_y_u=0.0
  found_zt       = .false.;  found_zw       = .false. 
  found_grid_x_t = .false.;  found_grid_y_t = .false.
  found_grid_x_u = .false.;  found_grid_y_u = .false.

  do i=1, ndim
     call mpp_get_atts(axes(i),name=name,len=len)
     if (len > 0) then
         allocate(data(len))
         call mpp_get_axis_data(axes(i),data=data)
         if(is_new_grid) then
            select case (trim(name))
            case ('zt')
                Grid%zt=data
                found_zt = .true.
            case ('zb')
                Grid%zw = data
                found_zw = .true.
            case ('grid_x_T')
                Grid%grid_x_t = data
                found_grid_x_t = .true.
            case ('grid_y_T')
                Grid%grid_y_t = data
                found_grid_y_t = .true.
            case ('grid_x_C')
                Grid%grid_x_u = data
                found_grid_x_u = .true.
            case ('grid_y_C')
                Grid%grid_y_u = data
                found_grid_y_u = .true.
            end select
         else
            select case (trim(name))
            case ('zt')
                Grid%zt=data
                found_zt = .true.
            case ('zw')
                Grid%zw = data
                found_zw = .true.
            case ('gridlon_t')
                Grid%grid_x_t = data
                found_grid_x_t = .true.
            case ('gridlat_t')
                Grid%grid_y_t = data
                found_grid_y_t = .true.
            case ('gridlon_vert_t')
                Grid%grid_x_u = data(2:ni+1)
                found_grid_x_u = .true.
            case ('gridlat_vert_t')
                Grid%grid_y_u = data(2:nj+1)
                found_grid_y_u = .true.
            end select
         endif
         deallocate(data)
     endif
  enddo
  
  if(is_new_grid) then
     if(.not. found_zt)       call mpp_error(FATAL,'ocean_grids_mod: axis zt is not in file INPUT/grid_spec.nc')
     if(.not. found_zw)       call mpp_error(FATAL,'ocean_grids_mod: axis zb is not in file INPUT/grid_spec.nc')
     if(.not. found_grid_x_t) call mpp_error(FATAL,'ocean_grids_mod: axis grid_x_T is not in file INPUT/grid_spec.nc')
     if(.not. found_grid_y_t) call mpp_error(FATAL,'ocean_grids_mod: axis grid_y_T is not in file INPUT/grid_spec.nc')
     if(.not. found_grid_x_u) call mpp_error(FATAL,'ocean_grids_mod: axis grid_x_C is not in file INPUT/grid_spec.nc')
     if(.not. found_grid_y_u) call mpp_error(FATAL,'ocean_grids_mod: axis grid_y_C is not in file INPUT/grid_spec.nc')
  else 
     if(.not. found_zt)       call mpp_error(FATAL,'ocean_grids_mod: axis zt is not in file INPUT/grid_spec.nc')
     if(.not. found_zw)       call mpp_error(FATAL,'ocean_grids_mod: axis zw is not in file INPUT/grid_spec.nc')
     if(.not. found_grid_x_t) call mpp_error(FATAL,'ocean_grids_mod: axis gridlon_t is not in file INPUT/grid_spec.nc')
     if(.not. found_grid_y_t) call mpp_error(FATAL,'ocean_grids_mod: axis gridlat_t is not in file INPUT/grid_spec.nc')
     if(.not. found_grid_x_u) call mpp_error(FATAL,'ocean_grids_mod: axis gridlon_vert_t is not in file INPUT/grid_spec.nc')
     if(.not. found_grid_y_u) call mpp_error(FATAL,'ocean_grids_mod: axis gridlat_vert_t is not in file INPUT/grid_spec.nc')
  endif

  Grid%xt=0.0;Grid%yt=0.0;Grid%xu=0.0;Grid%yu=0.0
  allocate(tmp(ni,nj))
  found_xt = .false. ; found_yt = .false.; found_xu = .false. ; found_yu = .false.

   print *,'after found_xt'

  do i=1, nvar
     call mpp_get_atts(fields(i),name=name)
     if(is_new_grid) then
        select case (trim(name))
            case ('x_T')
                call mpp_read(unit,fields(i), tmp)
                Grid%xt(isc:iec,jsc:jec) = tmp(isc+ioff:iec+ioff,jsc+joff:jec+joff)
                call update_boundaries(Domain,Grid,Grid%xt,tmp,0,0)
                found_xt = .true.
            case ('y_T')
                call mpp_read(unit,fields(i),tmp)
                Grid%yt(isc:iec,jsc:jec) = tmp(isc+ioff:iec+ioff,jsc+joff:jec+joff)
                call update_boundaries(Domain,Grid,Grid%yt,tmp,0,0)
                found_yt = .true.
            case ('x_C')
                call mpp_read(unit,fields(i), tmp)
                Grid%xu(isc:iec,jsc:jec) = tmp(isc+ioff:iec+ioff,jsc+joff:jec+joff)
                call update_boundaries(Domain,Grid,Grid%xu,tmp,-1,-1)
                found_xu = .true.
            case ('y_C')
                call mpp_read(unit,fields(i),tmp)
                Grid%yu(isc:iec,jsc:jec) = tmp(isc+ioff:iec+ioff,jsc+joff:jec+joff)
                call update_boundaries(Domain,Grid,Grid%yu,tmp,-1,-1)
                found_yu = .true.
        end select
     else
        select case (trim(name))
            case ('geolon_t')
                call mpp_read(unit,fields(i), tmp)
                Grid%xt(isc:iec,jsc:jec) = tmp(isc+ioff:iec+ioff,jsc+joff:jec+joff)
                call update_boundaries(Domain,Grid,Grid%xt,tmp,0,0)
                found_xt = .true.
            case ('geolat_t')
                call mpp_read(unit,fields(i),tmp)
                Grid%yt(isc:iec,jsc:jec) = tmp(isc+ioff:iec+ioff,jsc+joff:jec+joff)
                call update_boundaries(Domain,Grid,Grid%yt,tmp,0,0)
                found_yt = .true.
            case ('geolon_c')
                call mpp_read(unit,fields(i), tmp)
                Grid%xu(isc:iec,jsc:jec) = tmp(isc+ioff:iec+ioff,jsc+joff:jec+joff)
                call update_boundaries(Domain,Grid,Grid%xu,tmp,-1,-1)
                found_xu = .true.
            case ('geolat_c')
                call mpp_read(unit,fields(i),tmp)
                Grid%yu(isc:iec,jsc:jec) = tmp(isc+ioff:iec+ioff,jsc+joff:jec+joff)
                call update_boundaries(Domain,Grid,Grid%yu,tmp,-1,-1)
                found_yu = .true.
        end select
     endif
  enddo
   print *,'after do loop mpp_read'

  if(is_new_grid) then
     if(.not. found_xt) call mpp_error(FATAL,'ocean_grids_mod: field x_T is not in file INPUT/grid_spec.nc')  
     if(.not. found_yt) call mpp_error(FATAL,'ocean_grids_mod: field y_T is not in file INPUT/grid_spec.nc') 
     if(.not. found_xu) call mpp_error(FATAL,'ocean_grids_mod: field x_C is not in file INPUT/grid_spec.nc') 
     if(.not. found_yu) call mpp_error(FATAL,'ocean_grids_mod: field y_C is not in file INPUT/grid_spec.nc') 
  else
     if(.not. found_xt) call mpp_error(FATAL,'ocean_grids_mod: field geolon_t is not in file INPUT/grid_spec.nc')  
     if(.not. found_yt) call mpp_error(FATAL,'ocean_grids_mod: field geolat_t is not in file INPUT/grid_spec.nc') 
     if(.not. found_xu) call mpp_error(FATAL,'ocean_grids_mod: field geolon_c is not in file INPUT/grid_spec.nc') 
     if(.not. found_yu) call mpp_error(FATAL,'ocean_grids_mod: field geolat_c is not in file INPUT/grid_spec.nc') 
  endif
   print *,'after is_new_grid'

  if (debug_grid .or. verbose_init) then
      print *,'in debug_grid, isc=',isc,'iec=',iec,'jsc=',jsc,'jec=',jec,Grid%xt(isc,jsc)
!--> cpl deletion
!     write(stdout(),*) 'xt chksum = ',mpp_chksum(Grid%xt(isc:iec,jsc:jec))
!     write(stdout(),*) 'xu chksum = ',mpp_chksum(Grid%xu(isc:iec,jsc:jec))
!     write(stdout(),*) 'yt chksum = ',mpp_chksum(Grid%yt(isc:iec,jsc:jec))
!     write(stdout(),*) 'yu chksum = ',mpp_chksum(Grid%yu(isc:iec,jsc:jec))
!<-- cpl deletion
  endif

   print *,'after debug_grid'
 if (Grid%beta_plane .or. Grid%f_plane) then
    Grid%phiu(:,:) = Grid%f_plane_latitude/radian
    Grid%phit(:,:) = Grid%phiu(:,:)
 else
    Grid%phit(:,:) = Grid%yt(:,:)/radian
    Grid%phiu(:,:) = Grid%yu(:,:)/radian
 endif

   print *,'after grid'
 Grid%h1t(:,:) = radius*cos(Grid%phit(:,:))
 where (cos(Grid%phit) == 0.0) Grid%h1t = radius*abs(epsln)
 Grid%h1u(:,:) = radius*cos(Grid%phiu(:,:))
 where (cos(Grid%phiu) == 0.0) Grid%h1u = radius*abs(epsln)
 Grid%h2t(:,:) = radius
 Grid%h2u(:,:) = radius


  ! set cell widths (meters) and area (meters^2)

 Grid%dxt=0.0 ; Grid%dxu=0.0 ; Grid%dyt=0.0 ; Grid%dyu=0.0
 found_dxt = .false.; found_dyt = .false.;  found_dxu = .false.; found_dyu = .false.

    print *,'after set cell widths'
 do i=1, nvar
    call mpp_get_atts(fields(i),name=name)
    if(is_new_grid) then
       select case (trim(name))
       case ('ds_01_21_T')
           call mpp_read(unit, fields(i), tmp)
           Grid%dxt(isc:iec,jsc:jec) = tmp(isc+ioff:iec+ioff,jsc+joff:jec+joff)
           call update_boundaries(Domain,Grid,Grid%dxt,tmp,0,0)
           found_dxt = .true.
       case ('ds_01_21_C')
           call mpp_read(unit, fields(i), tmp)
           Grid%dxu(isc:iec,jsc:jec) = tmp(isc+ioff:iec+ioff,jsc+joff:jec+joff)
           call update_boundaries(Domain,Grid,Grid%dxu,tmp,-1,-1)
           found_dxu = .true.
       case ('ds_10_12_T')
           call mpp_read(unit, fields(i), tmp)
           Grid%dyt(isc:iec,jsc:jec) = tmp(isc+ioff:iec+ioff,jsc+joff:jec+joff)
           call update_boundaries(Domain,Grid,Grid%dyt,tmp,0,0)
           found_dyt = .true.
       case ('ds_10_12_C')
           call mpp_read(unit, fields(i), tmp)
           Grid%dyu(isc:iec,jsc:jec) = tmp(isc+ioff:iec+ioff,jsc+joff:jec+joff)
           call update_boundaries(Domain,Grid,Grid%dyu,tmp,-1,-1)
           found_dyu = .true.
       end select
    else
       select case (trim(name))
       case ('dxt')
           call mpp_read(unit, fields(i), tmp)
           Grid%dxt(isc:iec,jsc:jec) = tmp(isc+ioff:iec+ioff,jsc+joff:jec+joff)
           call update_boundaries(Domain,Grid,Grid%dxt,tmp,0,0)
           found_dxt = .true.
       case ('dxu')
           call mpp_read(unit, fields(i), tmp)
           Grid%dxu(isc:iec,jsc:jec) = tmp(isc+ioff:iec+ioff,jsc+joff:jec+joff)
           call update_boundaries(Domain,Grid,Grid%dxu,tmp,-1,-1)
           found_dxu = .true.
       case ('dyt')
           call mpp_read(unit, fields(i), tmp)
           Grid%dyt(isc:iec,jsc:jec) = tmp(isc+ioff:iec+ioff,jsc+joff:jec+joff)
           call update_boundaries(Domain,Grid,Grid%dyt,tmp,0,0)
           found_dyt = .true.
       case ('dyu')
           call mpp_read(unit, fields(i), tmp)
           Grid%dyu(isc:iec,jsc:jec) = tmp(isc+ioff:iec+ioff,jsc+joff:jec+joff)
           call update_boundaries(Domain,Grid,Grid%dyu,tmp,-1,-1)
           found_dyu = .true.
       end select
    endif
 enddo

  if(is_new_grid) then
     if(.not. found_dxt) call mpp_error(FATAL,'ocean_grids_mod: field ds_01_21_T is not in file INPUT/grid_spec.nc')  
     if(.not. found_dyt) call mpp_error(FATAL,'ocean_grids_mod: field ds_10_12_T is not in file INPUT/grid_spec.nc') 
     if(.not. found_dxu) call mpp_error(FATAL,'ocean_grids_mod: field ds_01_21_C is not in file INPUT/grid_spec.nc') 
     if(.not. found_dyu) call mpp_error(FATAL,'ocean_grids_mod: field ds_10_12_C is not in file INPUT/grid_spec.nc') 
  else
     if(.not. found_dxt) call mpp_error(FATAL,'ocean_grids_mod: field dxt is not in file INPUT/grid_spec.nc')  
     if(.not. found_dyt) call mpp_error(FATAL,'ocean_grids_mod: field dyt is not in file INPUT/grid_spec.nc') 
     if(.not. found_dxu) call mpp_error(FATAL,'ocean_grids_mod: field dxu is not in file INPUT/grid_spec.nc') 
     if(.not. found_dyu) call mpp_error(FATAL,'ocean_grids_mod: field dyu is not in file INPUT/grid_spec.nc') 
  endif

 Grid%dat(:,:)  = Grid%dxt(:,:)*Grid%dyt(:,:)
 Grid%dau(:,:)  = Grid%dxu(:,:)*Grid%dyu(:,:)

 if (debug_grid .or. verbose_init) then
   print *,'in debug_grid, dxt,isc=',isc,'iec=',iec,'jsc=',jsc,'jec=',jec,Grid%dxt(isc,jsc) 
!--> cpl deletion
!    write(stdout(),*) 'dxt chksum = ',mpp_chksum(Grid%dxt(isc:iec,jsc:jec))
!    write(stdout(),*) 'dxu chksum = ',mpp_chksum(Grid%dxu(isc:iec,jsc:jec))
!    write(stdout(),*) 'dyt chksum = ',mpp_chksum(Grid%dyt(isc:iec,jsc:jec))
!    write(stdout(),*) 'dyu chksum = ',mpp_chksum(Grid%dyu(isc:iec,jsc:jec))
!    write(stdout(),*) 'dat chksum = ',mpp_chksum(Grid%dat(isc:iec,jsc:jec))
!    write(stdout(),*) 'dau chksum = ',mpp_chksum(Grid%dau(isc:iec,jsc:jec))
!<-- cpl deletion
 endif

  ! set lengths at edges of grid cells (meters)

   print *,'after set lengths at edges of grid cells'
 allocate(tmp1_local(isd:ied,jsd:jed), tmp2_local(isd:ied,jsd:jed))

 Grid%dxtn=0.0; Grid%dytn=0.0; Grid%dxte=0.0; Grid%dyte=0.0
 Grid%dxun=0.0; Grid%dyun=0.0; Grid%dxue=0.0; Grid%dyue=0.0

 found_dxtn = .false. ; found_dytn = .false. ;  found_dxte = .false. ; found_dyte = .false. 
 found_dxun = .false. ; found_dyun = .false. ;  found_dxue = .false. ; found_dyue = .false. 
 do i=1, nvar
    call mpp_get_atts(fields(i),name=name)
    if(is_new_grid) then
       select case (trim(name))
       case ('ds_02_22_T')
           call mpp_read(unit, fields(i), tmp)
           Grid%dxtn(isc:iec,jsc:jec) = tmp(isc+ioff:iec+ioff,jsc+joff:jec+joff)
           call update_boundaries(Domain,Grid,Grid%dxtn,tmp,0,-1)
           found_dxtn = .true.
       case ('ds_00_02_C')
           call mpp_read(unit, fields(i), tmp)
           Grid%dytn(isc:iec,jsc:jec) = tmp(isc+ioff:iec+ioff,jsc+joff:jec+joff)
           call update_boundaries(Domain,Grid,Grid%dytn,tmp,0,-1)
           found_dytn = .true.
       case ('ds_00_20_C')
           call mpp_read(unit, fields(i), tmp)
           Grid%dxte(isc:iec,jsc:jec) = tmp(isc+ioff:iec+ioff,jsc+joff:jec+joff)
           call update_boundaries(Domain,Grid,Grid%dxte,tmp,-1,0) 
           found_dxte = .true.
       case ('ds_20_22_T')
           call mpp_read(unit, fields(i), tmp)
           Grid%dyte(isc:iec,jsc:jec) = tmp(isc+ioff:iec+ioff,jsc+joff:jec+joff)
           call update_boundaries(Domain,Grid,Grid%dyte,tmp,-1,0)
           found_dyte = .true.
       case ('ds_02_22_C')
           call mpp_read(unit, fields(i), tmp)
           Grid%dxun(isc:iec,jsc:jec) = tmp(isc+ioff:iec+ioff,jsc+joff:jec+joff)
           call update_boundaries(Domain,Grid,Grid%dxun,tmp,-1,-2)
           found_dxun = .true.
       case ('ds_11_12_C')
           call mpp_read(unit, fields(i), tmp) 
           tmp1_local(isc:iec,jsc:jec) = tmp(isc+ioff:iec+ioff,jsc+joff:jec+joff)
           do ii=1, nvar
              call mpp_get_atts(fields(ii),name=name2)
              if (trim(name2) == 'ds_10_11_C') call mpp_read(unit,fields(ii),Domain%domain2d,tmp2_local)
           enddo        
           call update_boundaries(Domain,Grid,tmp2_local,tmp,-1,-1)
           Grid%dyun(isc:iec,jsc:jec) = tmp1_local(isc:iec,jsc:jec)+tmp2_local(isc:iec,jsc+1:jec+1) 
           call mpp_global_field(Domain%domain2d, Grid%dyun, tmp)
           call update_boundaries(Domain,Grid,Grid%dyun,tmp,-1,-2)
       case ('ds_11_21_C')
           call mpp_read(unit, fields(i), tmp) 
           tmp1_local(isc:iec,jsc:jec) = tmp(isc+ioff:iec+ioff,jsc+joff:jec+joff)
           do ii=1, nvar
              call mpp_get_atts(fields(ii),name=name2)
              if (trim(name2) == 'ds_01_11_C') call mpp_read(unit,fields(ii),Domain%domain2d,tmp2_local)
           enddo        
           call update_boundaries(Domain,Grid,tmp2_local,tmp,-1,-1)
           Grid%dxue(isc:iec,jsc:jec) = tmp1_local(isc:iec,jsc:jec)+tmp2_local(isc+1:iec+1,jsc:jec) 
           call mpp_global_field(Domain%domain2d, Grid%dxue, tmp)
           call update_boundaries(Domain,Grid,Grid%dxue,tmp,-2,-1)
       case ('ds_20_22_C')
           call mpp_read(unit, fields(i), tmp)
           Grid%dyue(isc:iec,jsc:jec) = tmp(isc+ioff:iec+ioff,jsc+joff:jec+joff)
           call update_boundaries(Domain,Grid,Grid%dyue,tmp,-2,-1) 
           found_dyue = .true.
       end select
    else
       select case (trim(name))
       case ('dxtn')
           call mpp_read(unit, fields(i), tmp)
           Grid%dxtn(isc:iec,jsc:jec) = tmp(isc+ioff:iec+ioff,jsc+joff:jec+joff)
           call update_boundaries(Domain,Grid,Grid%dxtn,tmp,0,-1)
           found_dxtn = .true.
       case ('dytn')
           call mpp_read(unit, fields(i), tmp)
           Grid%dytn(isc:iec,jsc:jec) = tmp(isc+ioff:iec+ioff,jsc+joff:jec+joff)
           call update_boundaries(Domain,Grid,Grid%dytn,tmp,0,-1)
           found_dytn = .true.
       case ('dxte')
           call mpp_read(unit, fields(i), tmp)
           Grid%dxte(isc:iec,jsc:jec) = tmp(isc+ioff:iec+ioff,jsc+joff:jec+joff)
           call update_boundaries(Domain,Grid,Grid%dxte,tmp,-1,0) 
           found_dxte = .true.
       case ('dyte')
           call mpp_read(unit, fields(i), tmp)
           Grid%dyte(isc:iec,jsc:jec) = tmp(isc+ioff:iec+ioff,jsc+joff:jec+joff)
           call update_boundaries(Domain,Grid,Grid%dyte,tmp,-1,0)
           found_dyte = .true.
       case ('dxun')
           call mpp_read(unit, fields(i), tmp)
           Grid%dxun(isc:iec,jsc:jec) = tmp(isc+ioff:iec+ioff,jsc+joff:jec+joff)
           call update_boundaries(Domain,Grid,Grid%dxun,tmp,-1,-2)
           found_dxun = .true.
       case ('dyun')
           call mpp_read(unit, fields(i), tmp)
           Grid%dyun(isc:iec,jsc:jec) = tmp(isc+ioff:iec+ioff,jsc+joff:jec+joff)
           call update_boundaries(Domain,Grid,Grid%dyun,tmp,-1,-2)
           found_dyun = .true.
       case ('dxue')
           call mpp_read(unit, fields(i), tmp)
           Grid%dxue(isc:iec,jsc:jec) = tmp(isc+ioff:iec+ioff,jsc+joff:jec+joff)
           call update_boundaries(Domain,Grid,Grid%dxue,tmp,-2,-1) 
           found_dxue = .true.
       case ('dyue')
           call mpp_read(unit, fields(i), tmp)
           Grid%dyue(isc:iec,jsc:jec) = tmp(isc+ioff:iec+ioff,jsc+joff:jec+joff)
           call update_boundaries(Domain,Grid,Grid%dyue,tmp,-2,-1) 
           found_dyue = .true.
       end select
    endif
 enddo 

 if(is_new_grid) then
    if(.not. found_dxtn) call mpp_error(FATAL,'ocean_grids_mod: field ds_02_22_T is not in file INPUT/grid_spec.nc')  
    if(.not. found_dytn) call mpp_error(FATAL,'ocean_grids_mod: field ds_00_02_C is not in file INPUT/grid_spec.nc') 
    if(.not. found_dxte) call mpp_error(FATAL,'ocean_grids_mod: field ds_00_20_C is not in file INPUT/grid_spec.nc')  
    if(.not. found_dyte) call mpp_error(FATAL,'ocean_grids_mod: field ds_20_22_T is not in file INPUT/grid_spec.nc') 
    if(.not. found_dxun) call mpp_error(FATAL,'ocean_grids_mod: field ds_02_22_C is not in file INPUT/grid_spec.nc') 
    if(.not. found_dyue) call mpp_error(FATAL,'ocean_grids_mod: field ds_20_22_C is not in file INPUT/grid_spec.nc') 
 else
    if(.not. found_dxtn) call mpp_error(FATAL,'ocean_grids_mod: field dxtn is not in file INPUT/grid_spec.nc')  
    if(.not. found_dytn) call mpp_error(FATAL,'ocean_grids_mod: field dytn is not in file INPUT/grid_spec.nc') 
    if(.not. found_dxte) call mpp_error(FATAL,'ocean_grids_mod: field dxte is not in file INPUT/grid_spec.nc')  
    if(.not. found_dyte) call mpp_error(FATAL,'ocean_grids_mod: field dyte is not in file INPUT/grid_spec.nc') 
    if(.not. found_dxun) call mpp_error(FATAL,'ocean_grids_mod: field dxun is not in file INPUT/grid_spec.nc') 
    if(.not. found_dyun) call mpp_error(FATAL,'ocean_grids_mod: field dyun is not in file INPUT/grid_spec.nc') 
    if(.not. found_dxue) call mpp_error(FATAL,'ocean_grids_mod: field dxue is not in file INPUT/grid_spec.nc') 
    if(.not. found_dyue) call mpp_error(FATAL,'ocean_grids_mod: field dyue is not in file INPUT/grid_spec.nc') 
 endif

 if (debug_grid .or. verbose_init) then
!---cpl deletion
!    write(stdout(),*) 'dxtn chksum = ',mpp_chksum(Grid%dxtn(isc:iec,jsc:jec))
!    write(stdout(),*) 'dytn chksum = ',mpp_chksum(Grid%dytn(isc:iec,jsc:jec))
!    write(stdout(),*) 'dxte chksum = ',mpp_chksum(Grid%dxte(isc:iec,jsc:jec))
!    write(stdout(),*) 'dyte chksum = ',mpp_chksum(Grid%dyte(isc:iec,jsc:jec))
!    write(stdout(),*) 'dxun chksum = ',mpp_chksum(Grid%dxun(isc:iec,jsc:jec))
!    write(stdout(),*) 'dyun chksum = ',mpp_chksum(Grid%dyun(isc:iec,jsc:jec))
!    write(stdout(),*) 'dxue chksum = ',mpp_chksum(Grid%dxue(isc:iec,jsc:jec))
!    write(stdout(),*) 'dyue chksum = ',mpp_chksum(Grid%dyue(isc:iec,jsc:jec))
!<-- cpl deltion
 endif

  ! set reciprocals and related quantities
   print *,'set reciprocals and related quantities'
 do j=jsd,jed
    do i=isd,ied
       Grid%dxtr(i,j)       = 1.0/(Grid%dxt(i,j)+epsln)
       Grid%dxur(i,j)       = 1.0/(Grid%dxu(i,j)+epsln)
       Grid%dytr(i,j)       = 1.0/(Grid%dyt(i,j)+epsln)
       Grid%dyur(i,j)       = 1.0/(Grid%dyu(i,j)+epsln)
       Grid%dxuer(i,j)      = 1.0/(Grid%dxue(i,j)+epsln)
       Grid%dyuer(i,j)      = 1.0/(Grid%dyue(i,j)+epsln)
       Grid%dxunr(i,j)      = 1.0/(Grid%dxun(i,j)+epsln)
       Grid%dyunr(i,j)      = 1.0/(Grid%dyun(i,j)+epsln)
       Grid%dxter(i,j)      = 1.0/(Grid%dxte(i,j)+epsln)
       Grid%dytnr(i,j)      = 1.0/(Grid%dytn(i,j)+epsln)
       Grid%dyue_dxuer(i,j) = Grid%dyue(i,j)/(Grid%dxue(i,j)+epsln)
       Grid%dxun_dyunr(i,j) = Grid%dxun(i,j)/(Grid%dyun(i,j)+epsln)
       Grid%datr(i,j)       = 1.0/(Grid%dat(i,j)+epsln)
       Grid%daur(i,j)       = 1.0/(Grid%dau(i,j)+epsln)
       Grid%dater(i,j)      = 1.0/(Grid%dxte(i,j)*Grid%dyte(i,j)+epsln)
       Grid%datnr(i,j)      = 1.0/(Grid%dxtn(i,j)*Grid%dytn(i,j)+epsln)
    enddo
 enddo

  ! set quantities which account for rotation of basis vectors

 allocate(tmp_local(isd:ied,jsd:jed))
 tmp_local=1.0
 
 ! expanded version of dh2dx(:,:) = BDX_EU(tmp(:,:))
 do j=jsc,jec
    do i=isc,iec
       Grid%dh2dx(i,j) = (Grid%dyue(i,j)*tmp_local(i,j) - Grid%dyue(i-1,j)*tmp_local(i-1,j))*Grid%daur(i,j)
    enddo
 enddo
 if (.not. Grid%tripolar .and. jec+joff == nj) then
   do i=isc,iec
     Grid%dh2dx(i,jec) = 0.0
   enddo
 endif

 call mpp_global_field(Domain%domain2d,Grid%dh2dx,tmp)
 call update_boundaries(Domain,Grid,Grid%dh2dx,tmp,-1,-1)

 ! expanded version of dh1dy(:,:) = BDY_NU(tmp(:,:))
 do j=jsc,jec
    do i=isc,iec
       Grid%dh1dy(i,j) = (Grid%dxun(i,j)*tmp_local(i,j) - Grid%dxun(i,j-1)*tmp_local(i,j-1))*Grid%daur(i,j)
    enddo
 enddo
 if (.not. Grid%tripolar .and. jec+joff == nj) then
   do i=isc,iec
     Grid%dh1dy(i,jec) = 0.0
   enddo
 endif

 call mpp_global_field(Domain%domain2d,Grid%dh1dy,tmp)
 call update_boundaries(Domain,Grid,Grid%dh1dy,tmp,-1,-1)

    print *,'after update_boundaries'

 ! build distances from grid point to cell faces {meters}
 Grid%duw=0.0 ; Grid%due=0.0 ; Grid%dus=0.0 ; Grid%dun=0.0 ; Grid%dtw=0.0 ; Grid%dte=0.0 ; Grid%dts=0.0; Grid%dtn=0.0
 found_duw = .false.; found_due = .false.; found_dus = .false.; found_dun = .false. 
 found_dtw = .false.; found_dte = .false.; found_dts = .false.; found_dtn = .false.
 do i=1, nvar
    call mpp_get_atts(fields(i),name=name)
    if(is_new_grid) then
       select case (trim(name))
       case ('ds_01_11_C')
           call mpp_read(unit, fields(i), Domain%domain2d, Grid%duw)
           do ii=1, nvar
              call mpp_get_atts(fields(ii),name=name2)
              if (trim(name2) == 'ds_11_21_C') call mpp_read(unit,fields(ii),tmp)
           enddo
           call update_boundaries(Domain,Grid,Grid%duw,tmp,-1,-1)
           found_duw = .true.
       case ('ds_11_21_C')
           call mpp_read(unit, fields(i), Domain%domain2d, Grid%due)
           do ii=1, nvar
              call mpp_get_atts(fields(ii),name=name2)
              if (trim(name2) == 'ds_01_11_C') call mpp_read(unit,fields(ii),tmp)
           enddo
           call update_boundaries(Domain,Grid,Grid%due,tmp,-1,-1)
           found_due = .true.
       case ('ds_10_11_C')
           call mpp_read(unit, fields(i), Domain%domain2d, Grid%dus)
           do ii=1, nvar
              call mpp_get_atts(fields(ii),name=name2)
              if (trim(name2) == 'ds_11_12_C') call mpp_read(unit,fields(ii),tmp)
           enddo
           call update_boundaries(Domain,Grid,Grid%dus,tmp,-1,-1)
           found_dus = .true.
       case ('ds_11_12_C')
           call mpp_read(unit, fields(i), Domain%domain2d, Grid%dun)
           do ii=1, nvar
              call mpp_get_atts(fields(ii),name=name2)
              if (trim(name2) == 'ds_10_11_C') call mpp_read(unit,fields(ii),tmp)
           enddo
           call update_boundaries(Domain,Grid,Grid%dun,tmp,-1,-1)
           found_dun = .true.
       case ('ds_01_11_T')
           call mpp_read(unit, fields(i), Domain%domain2d, Grid%dtw)
           do ii=1, nvar
              call mpp_get_atts(fields(ii),name=name2)
              if (trim(name2) == 'ds_11_21_T') call mpp_read(unit,fields(ii),tmp)
           enddo
           call update_boundaries(Domain,Grid,Grid%dtw,tmp,0,0)
           found_dtw = .true.
       case ('ds_11_21_T')
           call mpp_read(unit, fields(i), Domain%domain2d, Grid%dte)
           do ii=1, nvar
              call mpp_get_atts(fields(ii),name=name2)
              if (trim(name2) == 'ds_01_11_T') call mpp_read(unit,fields(ii),tmp)
           enddo
           call update_boundaries(Domain,Grid,Grid%dte,tmp,0,0)
           found_dte = .true.
       case ('ds_10_11_T')
           call mpp_read(unit, fields(i), Domain%domain2d, Grid%dts)
           do ii=1, nvar
              call mpp_get_atts(fields(ii),name=name2)
              if (trim(name2) == 'ds_11_12_T') call mpp_read(unit,fields(ii),tmp)
           enddo
           call update_boundaries(Domain,Grid,Grid%dts,tmp,0,0)
           found_dts = .true.
       case ('ds_11_12_T')
           call mpp_read(unit, fields(i), Domain%domain2d, Grid%dtn)
           do ii=1, nvar
              call mpp_get_atts(fields(ii),name=name2)
              if (trim(name2) == 'ds_10_11_T') call mpp_read(unit,fields(ii),tmp)
           enddo
           call update_boundaries(Domain,Grid,Grid%dtn,tmp,0,0)
           found_dtn = .true.
       end select
    else
       select case (trim(name))
       case ('duw')
           call mpp_read(unit, fields(i), Domain%domain2d, Grid%duw)
           do ii=1, nvar
              call mpp_get_atts(fields(ii),name=name2)
              if (trim(name2) == 'due') call mpp_read(unit,fields(ii),tmp)
           enddo
           call update_boundaries(Domain,Grid,Grid%duw,tmp,-1,-1)
           found_duw = .true.
       case ('due')
           call mpp_read(unit, fields(i), Domain%domain2d, Grid%due)
           do ii=1, nvar
              call mpp_get_atts(fields(ii),name=name2)
              if (trim(name2) == 'duw') call mpp_read(unit,fields(ii),tmp)
           enddo
           call update_boundaries(Domain,Grid,Grid%due,tmp,-1,-1)
           found_due = .true.
       case ('dus')
           call mpp_read(unit, fields(i), Domain%domain2d, Grid%dus)
           do ii=1, nvar
              call mpp_get_atts(fields(ii),name=name2)
              if (trim(name2) == 'dun') call mpp_read(unit,fields(ii),tmp)
           enddo
           call update_boundaries(Domain,Grid,Grid%dus,tmp,-1,-1)
           found_dus = .true.
       case ('dun')
           call mpp_read(unit, fields(i), Domain%domain2d, Grid%dun)
           do ii=1, nvar
              call mpp_get_atts(fields(ii),name=name2)
              if (trim(name2) == 'dus') call mpp_read(unit,fields(ii),tmp)
           enddo
           call update_boundaries(Domain,Grid,Grid%dun,tmp,-1,-1)
           found_dun = .true.
       case ('dtw')
           call mpp_read(unit, fields(i), Domain%domain2d, Grid%dtw)
           do ii=1, nvar
              call mpp_get_atts(fields(ii),name=name2)
              if (trim(name2) == 'dte') call mpp_read(unit,fields(ii),tmp)
           enddo
           call update_boundaries(Domain,Grid,Grid%dtw,tmp,0,0)
           found_dtw = .true.
       case ('dte')
           call mpp_read(unit, fields(i), Domain%domain2d, Grid%dte)
           do ii=1, nvar
              call mpp_get_atts(fields(ii),name=name2)
              if (trim(name2) == 'dtw') call mpp_read(unit,fields(ii),tmp)
           enddo
           call update_boundaries(Domain,Grid,Grid%dte,tmp,0,0)
           found_dte = .true.
       case ('dts')
           call mpp_read(unit, fields(i), Domain%domain2d, Grid%dts)
           do ii=1, nvar
              call mpp_get_atts(fields(ii),name=name2)
              if (trim(name2) == 'dtn') call mpp_read(unit,fields(ii),tmp)
           enddo
           call update_boundaries(Domain,Grid,Grid%dts,tmp,0,0)
           found_dts = .true.
       case ('dtn')
           call mpp_read(unit, fields(i), Domain%domain2d, Grid%dtn)
           do ii=1, nvar
              call mpp_get_atts(fields(ii),name=name2)
              if (trim(name2) == 'dts') call mpp_read(unit,fields(ii),tmp)
           enddo
           call update_boundaries(Domain,Grid,Grid%dtn,tmp,0,0)
           found_dtn = .true.
       end select
    endif
 enddo  

 if(is_new_grid) then
    if(.not. found_duw) call mpp_error(FATAL,'ocean_grids_mod: field ds_01_11_C is not in file INPUT/grid_spec.nc')  
    if(.not. found_due) call mpp_error(FATAL,'ocean_grids_mod: field ds_11_21_C is not in file INPUT/grid_spec.nc') 
    if(.not. found_dus) call mpp_error(FATAL,'ocean_grids_mod: field ds_10_11_C is not in file INPUT/grid_spec.nc')  
    if(.not. found_dun) call mpp_error(FATAL,'ocean_grids_mod: field ds_11_12_C is not in file INPUT/grid_spec.nc') 
    if(.not. found_dtw) call mpp_error(FATAL,'ocean_grids_mod: field ds_01_11_T is not in file INPUT/grid_spec.nc')  
    if(.not. found_dte) call mpp_error(FATAL,'ocean_grids_mod: field ds_11_21_T is not in file INPUT/grid_spec.nc') 
    if(.not. found_dts) call mpp_error(FATAL,'ocean_grids_mod: field ds_10_11_T is not in file INPUT/grid_spec.nc')  
    if(.not. found_dtn) call mpp_error(FATAL,'ocean_grids_mod: field ds_11_12_T is not in file INPUT/grid_spec.nc') 
 else
    if(.not. found_duw) call mpp_error(FATAL,'ocean_grids_mod: field duw is not in file INPUT/grid_spec.nc')  
    if(.not. found_due) call mpp_error(FATAL,'ocean_grids_mod: field due is not in file INPUT/grid_spec.nc') 
    if(.not. found_dus) call mpp_error(FATAL,'ocean_grids_mod: field dus is not in file INPUT/grid_spec.nc')  
    if(.not. found_dun) call mpp_error(FATAL,'ocean_grids_mod: field dun is not in file INPUT/grid_spec.nc') 
    if(.not. found_dtw) call mpp_error(FATAL,'ocean_grids_mod: field dtw is not in file INPUT/grid_spec.nc')  
    if(.not. found_dte) call mpp_error(FATAL,'ocean_grids_mod: field dte is not in file INPUT/grid_spec.nc') 
    if(.not. found_dts) call mpp_error(FATAL,'ocean_grids_mod: field dts is not in file INPUT/grid_spec.nc')  
    if(.not. found_dtn) call mpp_error(FATAL,'ocean_grids_mod: field dtn is not in file INPUT/grid_spec.nc') 
 endif


 if (debug_grid .or. verbose_init) then
!--> cpl deletion
!    write(stdout(),*) 'dtn+dts chksum = ',mpp_chksum(Grid%dtn(isc:iec,jsc:jec)+Grid%dts(isc:iec,jsc:jec))
!    write(stdout(),*) 'dun+dus chksum = ',mpp_chksum(Grid%dun(isc:iec,jsc:jec)+Grid%dus(isc:iec,jsc:jec))
!    write(stdout(),*) 'dte+dtw chksum = ',mpp_chksum(Grid%dte(isc:iec,jsc:jec)+Grid%dtw(isc:iec,jsc:jec))
!    write(stdout(),*) 'due+duw chksum = ',mpp_chksum(Grid%due(isc:iec,jsc:jec)+Grid%duw(isc:iec,jsc:jec))
!    write(stdout(),*) 'dte chksum = ',mpp_chksum(Grid%dte(isc:iec,jsc:jec) )
!    write(stdout(),*) 'dtw chksum = ',mpp_chksum(Grid%dtw(isc:iec,jsc:jec) )
!    write(stdout(),*) 'due chksum = ',mpp_chksum(Grid%due(isc:iec,jsc:jec) )
!    write(stdout(),*) 'duw chksum = ',mpp_chksum(Grid%duw(isc:iec,jsc:jec) )
!<-- cpl deletion
 endif


! ensure the sphere is tiled properly by u-cells at the southern boundary
 if (jsc+joff == 1) then
   
     call mpp_error(NOTE, '==>Note from ocean_grids_mod (set_ocean_hgrid_arrays): altering U-grid arrays at j=0')

     Grid%dxun(:,0) = Grid%dxun(:,1)*Grid%h1t(:,1)/Grid%h1t(:,2)
     Grid%dun(:,0)  = Grid%dyte(:,1) - Grid%dus(:,1)

     if (Grid%beta_plane .or. Grid%f_plane) then
         Grid%phiu(:,0) = Grid%f_plane_latitude/radian
     else
         Grid%phiu(:,0) = Grid%yu(:,1)/radian - Grid%dyte(:,1)/radius
     endif

     Grid%phiu(:,0) = max(min(Grid%phiu(:,0),90.0/radian),-90.0/radian)

     where (cos(Grid%phiu(:,0)) == 0.0)
         Grid%h1u(:,0) = radius*abs(epsln)
     elsewhere
         Grid%h1u(:,0) = radius*cos(Grid%phiu(:,0))
     end where

     Grid%due(:,0) = Grid%due(:,1)*Grid%h1u(:,0)/Grid%h1u(:,1)
     Grid%duw(:,0) = Grid%duw(:,1)*Grid%h1u(:,0)/Grid%h1u(:,1)
     Grid%dus(:,0) = Grid%dau(:,0)/(Grid%duw(:,0)+Grid%due(:,0)+epsln) - Grid%dun(:,0)

 endif

  ! calculate rotation angles on velocity points
  Grid%sin_rot=0.0;Grid%cos_rot=0.0
  found_angle = .false.; found_sin_rot = .false.; found_cos_rot = .false.
  if(is_new_grid) then
     do i=1, nvar
        call mpp_get_atts(fields(i),name=name)
        select case (trim(name))
        case ('angle_C')
           call mpp_read(unit,fields(i),tmp)
        end select
     enddo
     allocate(sin_tmp(ni,nj), cos_tmp(ni,nj))

     sin_tmp = sin(tmp)
     cos_tmp = cos(tmp)
     Grid%sin_rot(isc:iec,jsc:jec) = sin(tmp(isc+ioff:iec+ioff,jsc+joff:jec+joff))
     Grid%cos_rot(isc:iec,jsc:jec) = cos(tmp(isc+ioff:iec+ioff,jsc+joff:jec+joff))
     call update_boundaries(Domain,Grid,Grid%sin_rot,sin_tmp,-1,-1)
     call update_boundaries(Domain,Grid,Grid%cos_rot,cos_tmp,-1,-1)
     deallocate(sin_tmp, cos_tmp)
     found_angle = .true.
  else
     do i=1, nvar
        call mpp_get_atts(fields(i),name=name)
        select case (trim(name))
        case ('sin_rot')
           call mpp_read(unit,fields(i),tmp)
           Grid%sin_rot(isc:iec,jsc:jec) = tmp(isc+ioff:iec+ioff,jsc+joff:jec+joff)
           call update_boundaries(Domain,Grid,Grid%sin_rot,tmp,-1,-1)
           found_sin_rot = .true.
        case ('cos_rot')
           call mpp_read(unit,fields(i),tmp)
           Grid%cos_rot(isc:iec,jsc:jec) = tmp(isc+ioff:iec+ioff,jsc+joff:jec+joff)
           call update_boundaries(Domain,Grid,Grid%cos_rot,tmp,-1,-1)
           found_cos_rot = .true.
        end select
     enddo
  endif

  if(is_new_grid) then
     if(.not. found_angle)   call mpp_error(FATAL,'ocean_grids_mod: field angle_C is not in file INPUT/grid_spec.nc')    
  else
     if(.not. found_sin_rot) call mpp_error(FATAL,'ocean_grids_mod: field sin_rot is not in file INPUT/grid_spec.nc')
     if(.not. found_cos_rot) call mpp_error(FATAL,'ocean_grids_mod: field cos_rot is not in file INPUT/grid_spec.nc')
  endif

   print *,'endof set_ocean_hgrid_arrays,unit=',unit
  
  deallocate(tmp, tmp_local, tmp1_local, tmp2_local, axes, fields)


end subroutine set_ocean_hgrid_arrays
! </SUBROUTINE> NAME="set_ocean_hgrid_arrays"


!#######################################################################
! <SUBROUTINE NAME="set_ocean_vgrid_arrays">
!
! <DESCRIPTION>
! Compute vertical (and some horizontal) grids for ocean model. 
! Also compute axes information for diagnostic manager.  
! </DESCRIPTION>
!
subroutine set_ocean_vgrid_arrays (Time, Domain, Grid, obc)

  type(ocean_time_type), intent(in)         :: Time
  type(ocean_domain_type), intent(in)       :: Domain
  type(ocean_grid_type), intent(inout)      :: Grid
  logical, intent(in), optional             :: obc

  real    :: ocnp
  real    :: ht_fax, ht_fay, ht_jp1_fax, ht_ip1_fay
  integer :: i, j, k, n, kzt, kzu, kb
  integer :: id_ztp
  integer :: id_geolon_t
  integer :: id_geolat_t
  integer :: id_geolon_uv
  integer :: id_geolat_uv
  integer :: id_area_t
  integer :: id_area_u

  logical :: used
  logical :: have_obc=.false.
  
  if (PRESENT(obc)) have_obc = obc
  
#ifndef STATIC_MEMORY
  allocate (Grid%tcella(nk))
  allocate (Grid%ucella(nk))
  allocate (Grid%tmask(isd:ied,jsd:jed,nk))
  allocate (Grid%umask(isd:ied,jsd:jed,nk))

  allocate (Grid%dht_dx(isd:ied,jsd:jed))
  allocate (Grid%dht_dy(isd:ied,jsd:jed))

  allocate (Grid%dzt(nk))
  allocate (Grid%dzw(0:nk))

  allocate (Grid%fracdz(nk,0:1)) 
  allocate (Grid%dzwr(0:nk))
#endif

  ! set depth arrays
  Grid%dzw(0)         = Grid%zt(1)
  Grid%dzw(1)         = Grid%zt(2)-Grid%zt(1)
  Grid%dzwr(0)        = 1.0/Grid%dzw(0)
  Grid%dzwr(1)        = 1.0/Grid%dzw(1)
  Grid%dzt(1)         = Grid%zw(1)
  Grid%fracdz(1,0)    = Grid%zt(1)/Grid%dzt(1)  
  Grid%fracdz(1,1)    = (Grid%zw(1) - Grid%zt(1))/Grid%dzt(1)  
  do k=2,nk-1
     Grid%dzt(k)      = Grid%zw(k)-Grid%zw(k-1)
     Grid%dzw(k)      = Grid%zt(k+1)-Grid%zt(k)
     Grid%dzwr(k)     = 1.0/Grid%dzw(k)
     Grid%fracdz(k,0) = (Grid%zt(k) - Grid%zw(k-1))/Grid%dzt(k)
     Grid%fracdz(k,1) = (Grid%zw(k) - Grid%zt(k))/Grid%dzt(k)
  enddo
  Grid%dzt(nk)        = Grid%zw(nk)-Grid%zw(nk-1)
  Grid%dzw(nk)        = Grid%zw(nk) - Grid%zt(nk)  
  Grid%dzwr(nk)       = 1.0/Grid%dzw(nk)
  Grid%fracdz(nk,0)   = (Grid%zt(nk) - Grid%zw(nk-1))/Grid%dzt(nk)
  Grid%fracdz(nk,1)   = (Grid%zw(nk) - Grid%zt(nk))/Grid%dzt(nk)


 ! construct T cell and U cell land/sea masks
 do k=1,nk
    do j=jsd,jed
       do i=isd,ied
          if (Grid%kmt(i,j) .ge. k) then
             Grid%tmask(i,j,k) = 1.0
          else
             Grid%tmask(i,j,k) = 0.0
          endif
          if (Grid%kmu(i,j) .ge. k) then
             Grid%umask(i,j,k) = 1.0
          else
             Grid%umask(i,j,k) = 0.0
          endif
       enddo
    enddo
 enddo

  ! compute surface area and volume of ocean (T cells and U cells)
  ocnp           = 0.0 !number of ocean T cells
  Grid%tcellv    = 0.0 !total ocean volume on T cells
  Grid%ucellv    = 0.0 !total ocean volume on U cells
  Grid%tcella(:) = 0.0 !total ocean surface area on T cells in level k
  Grid%ucella(:) = 0.0 !total ocean surface area on U cells in level k

  ! exclude points at open boundaries from diagnostics
  do j=jsc,jec
    do i=isc,iec
       kzt = Grid%kmt(i,j)
       if(have_obc) kzt = kzt * Grid%obc_tmask(i,j)
       if (kzt .gt. 0) then
           do k=1,kzt
              Grid%tcella(k) = Grid%tcella(k) + Grid%dat(i,j)
           enddo
           Grid%tcellv = Grid%tcellv + Grid%dat(i,j)*Grid%ht(i,j)
           ocnp   = ocnp + kzt
       endif
    enddo
 enddo

 do j=jsc,jec
    do i=isc,iec
       kzu = Grid%kmu(i,j)
       if(have_obc) kzu = kzu * Grid%obc_umask(i,j)
       if (kzu .gt. 0) then
           do k=1,kzu
              Grid%ucella(k) = Grid%ucella(k) + Grid%dau(i,j)
           enddo
           Grid%ucellv = Grid%ucellv + Grid%dau(i,j)*Grid%hu(i,j)
       endif
    enddo
 enddo

 call mpp_sum (Grid%tcella, nk)
 call mpp_sum (Grid%ucella, nk)
 call mpp_sum (Grid%tcellv)
 call mpp_sum (Grid%ucellv)
 call mpp_sum (ocnp)

 if (debug_grid .or. verbose_init) then
     write (stdout(),'(/a,es24.1,a)')  '  Number of wet ocean tracer points      =', ocnp
     write (stdout(),'(/a,es24.17,a)') '  Ocean volume with eta_t=0.0 (T-cells)  =', Grid%tcellv,    ' m^3 (not bit reproducible)'
     write (stdout(),'(a,es24.17,a)')  '  Ocean surface area          (T-cells)  =', Grid%tcella(1), ' m^2 (not bit reproducible)'
     write (stdout(),'(/a,es24.17,a)') '  Ocean volume with eta_u=0.0 (U-cells)  =', Grid%ucellv,    ' m^3 (not bit reproducible)'
     write (stdout(),'(a,es24.17,a/)') '  Ocean surface area          (U-cells)  =', Grid%ucella(1), ' m^2 (not bit reproducible)'
 endif
 
 ! compute bitwise reproducible surface areas
 if(have_obc) then 
   Grid%tcellsurf = mpp_global_sum(Domain%domain2d,Grid%dat(:,:)*Grid%tmask(:,:,1)*Grid%obc_tmask(:,:), BITWISE_EXACT_SUM)
   Grid%ucellsurf = mpp_global_sum(Domain%domain2d,Grid%dau(:,:)*Grid%umask(:,:,1)*Grid%obc_umask(:,:), BITWISE_EXACT_SUM)
 else 
   Grid%tcellsurf = mpp_global_sum(Domain%domain2d,Grid%dat(:,:)*Grid%tmask(:,:,1), BITWISE_EXACT_SUM)
   Grid%ucellsurf = mpp_global_sum(Domain%domain2d,Grid%dau(:,:)*Grid%umask(:,:,1), BITWISE_EXACT_SUM)
 endif 

  ! compute topographic slopes at U-cell. Note: cannot use 
  ! operators as ocean_operators_init has not yet been called.
  do j=jsc,jec
    do i=isc,iec
      ht_fax           = 0.5*( Grid%ht(i,j) + Grid%ht(i+1,j))
      ht_fay           = 0.5*( Grid%ht(i,j) + Grid%ht(i,j+1))
      ht_jp1_fax       = 0.5*( Grid%ht(i,j+1) + Grid%ht(i+1,j+1))
      ht_ip1_fay       = 0.5*( Grid%ht(i+1,j) + Grid%ht(i+1,j+1))
      Grid%dht_dx(i,j) = (ht_ip1_fay-ht_fay)*Grid%dxur(i,j)*Grid%umask(i,j,1)
      Grid%dht_dy(i,j) = (ht_jp1_fax-ht_fax)*Grid%dyur(i,j)*Grid%umask(i,j,1)
    enddo
  enddo


  ! set up the axis definition for variables
  call axes_info(Domain, Grid)


  ! output grid information to the diagnostic manager
  id_geolon_t = register_static_field('ocean_model','geolon_t',&
       Grid%tracer_axes(1:2), 'tracer longitude','degrees_E',&
       range=(/-281.,361./))
  if (id_geolon_t > 0)  used = send_data(id_geolon_t,Grid%xt(isc:iec,jsc:jec),  Time%model_time)

  id_geolat_t = register_static_field('ocean_model','geolat_t',&
       Grid%tracer_axes(1:2), 'tracer latitude','degrees_N',&
       range=(/-91.,91./))  
  if (id_geolat_t > 0)  used = send_data(id_geolat_t,Grid%yt(isc:iec,jsc:jec),  Time%model_time)

  id_geolon_uv = register_static_field('ocean_model','geolon_c',&
       Grid%vel_axes_uv(1:2), 'uv longitude','degrees_E',&
       range=(/-281.,361./))
  if (id_geolon_uv > 0) used = send_data(id_geolon_uv,Grid%xu(isc:iec,jsc:jec), Time%model_time)

  id_geolat_uv = register_static_field('ocean_model','geolat_c',&
       Grid%vel_axes_uv(1:2), 'uv latitude','degrees_N',&
       range=(/-91.,91./))
  if (id_geolat_uv > 0) used = send_data(id_geolat_uv,Grid%yu(isc:iec,jsc:jec), Time%model_time)

  id_area_t = register_static_field('ocean_model','area_t',&
       Grid%tracer_axes(1:2), 'tracer cell area','m^2',&
       range=(/0.,1.e10/))
  if (id_area_t > 0)    used = send_data(id_area_t,Grid%dat(isc:iec,jsc:jec),   Time%model_time)

  id_area_u = register_static_field('ocean_model','area_u',&
       Grid%vel_axes_uv(1:2), 'velocity cell area','m^2',&
       range=(/0.,1.e10/))  
  if (id_area_u > 0)    used = send_data(id_area_u,Grid%dau(isc:iec,jsc:jec),   Time%model_time)    

  id_ht  = register_static_field ('ocean_model', 'ht', Grid%tracer_axes(1:2), 'ocean depth on t-cells', 'm',&
                                  missing_value=-1.e10, range=(/-1.e10,1.e10/))
  if (id_ht > 0) used = send_data (id_ht, Grid%ht(isc:iec,jsc:jec), &
                        Time%model_time, rmask=Grid%tmask(isc:iec,jsc:jec,1))

  id_hu  = register_static_field ('ocean_model', 'hu', Grid%vel_axes_uv(1:2), 'ocean depth on u-cells', 'm',&
                                  missing_value=-1.e10, range=(/-1.e10,1.e10/))
  if (id_hu > 0) used = send_data (id_hu, Grid%hu(isc:iec,jsc:jec), &
                        Time%model_time, rmask=Grid%umask(isc:iec,jsc:jec,1))

  id_dht_dx  = register_static_field ('ocean_model', 'dht_dx', Grid%vel_axes_uv(1:2), 'd(ht)/dx on u-cells', 'm/m',&
                                      missing_value=-1.e10, range=(/-1.e10,1.e10/))
  if (id_dht_dx > 0) used = send_data (id_dht_dx, Grid%dht_dx(isc:iec,jsc:jec), &
                            Time%model_time, rmask=Grid%umask(isc:iec,jsc:jec,1))

  id_dht_dy  = register_static_field ('ocean_model', 'dht_dy', Grid%vel_axes_uv(1:2), 'd(ht)/dy on u-cells', 'm/m',&
                                      missing_value=-1.e10, range=(/-1.e10,1.e10/))
  if (id_dht_dy > 0) used = send_data (id_dht_dy, Grid%dht_dy(isc:iec,jsc:jec), &
                            Time%model_time, rmask=Grid%umask(isc:iec,jsc:jec,1))


end subroutine set_ocean_vgrid_arrays
! </SUBROUTINE> NAME="set_ocean_vgrid_arrays"



!#######################################################################
! <SUBROUTINE NAME="axes_info">
!
! <DESCRIPTION>
! Set up the definitions of the axes 
! </DESCRIPTION>
!
subroutine axes_info(Domain, Grid)

  type(ocean_domain_type), intent(in)       :: Domain
  type(ocean_grid_type), intent(inout)      :: Grid

  integer :: id_xt , id_yt, id_xu , id_yu
  integer :: id_zt_bounds, id_zt,  id_zw_bounds, id_zw
  integer :: i, j, k
  real, allocatable, dimension(:) :: zt_bounds
  real, allocatable, dimension(:) :: zw_bounds

  ! horizontal axes
  ! Note: 1d grid arrays are not in general representative of the grid.
  ! Need to generalize diagnostics to employ CF conventions.

  id_xt = diag_axis_init ('xt_'//trim(Grid%name),Grid%grid_x_t,'degrees_E','x','tcell longitude',set_name='ocean',&
                                  Domain2=Domain%domain2d)
  id_yt = diag_axis_init ('yt_'//trim(Grid%name),Grid%grid_y_t,'degrees_N','y','tcell latitude',set_name='ocean',&
                                  Domain2=Domain%domain2d)
  id_xu = diag_axis_init ('xu_'//trim(Grid%name),Grid%grid_x_u,'degrees_E','x','ucell longitude',set_name='ocean',&
                                  Domain2=Domain%domain2d)
  id_yu = diag_axis_init ('yu_'//trim(Grid%name),Grid%grid_y_u,'degrees_N','y','ucell latitude',set_name='ocean',&
                                  Domain2=Domain%domain2d)


  ! vertical axes (needs to be re-thought for partial cells)

  allocate(zt_bounds(Grid%nk+1))
  allocate(zw_bounds(Grid%nk+1))
  zt_bounds(1) = Grid%zt(1)-Grid%dzt(1)/2.0
  zw_bounds(1) = Grid%zw(1)-Grid%dzw(1)/2.0
  do k=2,Grid%nk+1
     zt_bounds(k)=zt_bounds(k-1)+Grid%dzt(k-1)
     zw_bounds(k)=zw_bounds(k-1)+Grid%dzw(k-1)
  enddo

  id_zt_bounds = diag_axis_init ('zt_edges_'//trim(Grid%name), zt_bounds, 'meters', 'z',&
                                 'tcell depth edges',direction=-1, set_name='ocean')

  id_zt = diag_axis_init ('zt_'//trim(Grid%name), Grid%zt, 'meters', 'z','tcell depth',&
                                 edges = id_zt_bounds,direction=-1, set_name='ocean')

  id_zw_bounds = diag_axis_init ( 'zw_edges_'//trim(Grid%name), zw_bounds, 'meters', 'z',&
                                 'ucell depth edges',direction=-1, set_name='ocean')

  id_zw = diag_axis_init ( 'zw_'//trim(Grid%name), Grid%zw, 'meters', 'z','ucell depth',&
                                 edges = id_zw_bounds,direction=-1, set_name='ocean')

  ! attributes for variables

  Grid%tracer_axes        = (/ id_xt, id_yt, id_zt /)
  Grid%tracer_axes_flux_x = (/ id_xu, id_yt, id_zt /)
  Grid%tracer_axes_flux_y = (/ id_xt, id_yu, id_zt /)

  Grid%vel_axes_uv        = (/ id_xu, id_yu, id_zt /)
  Grid%vel_axes_flux_x    = (/ id_xt, id_yu, id_zt /)
  Grid%vel_axes_flux_y    = (/ id_xu, id_yt, id_zt /)

  Grid%tracer_axes_wt     = (/ id_xt, id_yt, id_zw /)
  Grid%vel_axes_wu        = (/ id_xu, id_yu, id_zw /)
  Grid%vel_axes_wt        = (/ id_xt, id_yt, id_zw /)

  deallocate(zt_bounds, zw_bounds)

end subroutine axes_info
! </SUBROUTINE> NAME="axes_info"


!#######################################################################
! <SUBROUTINE NAME="update_boundaries">
!
! <DESCRIPTION>
! Set halo points at model boundaries equal to values at boundaries 
! if no grid connectivity exists. If model is connected along
! boundary then use mpp_update_domains.
! </DESCRIPTION>
!
subroutine update_boundaries(Domain, Grid, field, global_field, i_offset, j_offset)  

type(ocean_domain_type), intent(inout)          :: Domain
type(ocean_grid_type), intent(in)               :: Grid
real, intent(inout), dimension(isd:ied,jsd:jed) :: field
real, intent(in), dimension(isg:ieg,jsg:jeg)    :: global_field
integer, intent(in)                             :: i_offset, j_offset
integer                                         :: i,j, ii

call mpp_update_domains(field,Domain%domain2d)

! fill e/w halo points along computational rows
if (iec + Domain%ioff == ieg .and. .not. Grid%cyclic) then
   do i=iec+1,iec+xhalo
      field(i,jsc:jec) = field(iec,jsc:jec)
   enddo
endif
if (isc + Domain%ioff == isg .and. .not. Grid%cyclic) then
   do i=isc-xhalo,isc-1
      field(i,jsc:jec) = field(isc,jsc:jec)
   enddo
endif

! fill southernmost halo rows
if (jsc + Domain%joff == jsg) then
   do j=jsc-yhalo,jsc-1
      field(:,j) = field(:,jsc)
   enddo
endif

! fill northernmost halo rows
if (jec + Domain%joff == jeg) then
   if (.not. Grid%tripolar) then
      do j=jec+1,jec+yhalo
         field(:,j) = field(:,jec)
      enddo
   else
      do j=jec+1,jec+yhalo
         do i=isc-xhalo,iec+xhalo
            ii = ni-(i+Domain%ioff)+1+i_offset;if (ii.le.0) ii=ii+ni;if (ii.gt.ni) ii=ii-ni
            field(i,j) = global_field(ii,2*nj-(j+Domain%joff)+1+j_offset)
         enddo
      enddo
   endif
endif

! fill halo points at west boundary along south halo rows 
if (.not. Grid%cyclic .and. isc + Domain%ioff == isg .and. jsc + Domain%joff /= jsg ) then 
   do j=jsc-yhalo,jsc-1
      do i=isc-xhalo,isc-1
         field(i,j) = field(isc,j)
      enddo
   enddo
endif

! fill halo points at west boundary along north halo rows
if (.not. Grid%cyclic .and. isc + Domain%ioff == isg .and. jec + Domain%joff /= jeg ) then 
   do j=jec+1,jec+yhalo
      do i=isc-xhalo,isc-1
         field(i,j) = field(isc,j)
      enddo
   enddo
endif

! fill halo points at east boundary along south halo rows
if (.not. Grid%cyclic .and. iec + Domain%ioff == ieg .and. jsc + Domain%joff /= jsg) then 
   do j=jsc-yhalo,jsc-1
      do i=iec+1,iec+xhalo
         field(i,j) = field(iec,j)
      enddo
   enddo
endif

! fill halo points east boundary along north halo rows (I know, it is confusing !!)
if (.not. Grid%cyclic .and. iec + Domain%ioff == ieg .and. jec + Domain%joff /= jeg) then 
   do j=jec+1,jec+yhalo
      do i=iec+1,iec+xhalo
         field(i,j) = field(iec,j)
      enddo
   enddo
endif

return

end subroutine update_boundaries
! </SUBROUTINE> NAME="update_boundaries"





end module ocean_grids_mod
