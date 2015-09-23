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
program river_regrid
  !-----------------------------------------------------------------------
  !                   GNU General Public License                        
  !                                                                      
  ! This program is free software; you can redistribute it and/or modify it and  
  ! are expected to follow the terms of the GNU General Public License  
  ! as published by the Free Software Foundation; either version 2 of   
  ! the License, or (at your option) any later version.                 
  !                                                                      
  ! MOM is distributed in the hope that it will be useful, but WITHOUT    
  ! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  
  ! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    
  ! License for more details.                                           
  !                                                                      
  ! For the full text of the GNU General Public License,                
  ! write to: Free Software Foundation, Inc.,                           
  !           675 Mass Ave, Cambridge, MA 02139, USA.                   
  ! or see:   http://www.gnu.org/licenses/gpl.html                      
  !-----------------------------------------------------------------------
  !
  !<CONTACT EMAIL="Zhi.Liang@noaa.gov"> Zhi Liang  </CONTACT>
  !<CONTACT EMAIL="Kirsten.Findell@noaa.gov"> Kirsten Findell </CONTACT>

  !<OVERVIEW>
  !  This program can remap river network data from spherical grid onto another
  !  spherical grid.
  !</OVERVIEW>
  !<DESCRIPTION>
  !  The program expects to read river network data from a netcdf file, which 
  !  is specified by the namelist variable "river_input_file". This file should 
  !  contains field 'cellarea', 'tocell', 'fromcell', 'basin', 'basincells', 
  !  'order', 'travel', 'subL', 'subA','disttomouth', 'disttoocean' and 'dx'. 
  !  The program will remap the data in "river_input_file" onto the grid, which
  !  is specified by namelist land_grid_file. The grid file should contains the same
  !  grid as the land modle to be run. The output is stored in a netcdf file, which
  !  is specified by river_output_file.
  !</DESCRIPTION>

  use mpp_mod,          only : mpp_error, FATAL, NOTE
  use mpp_io_mod,       only : mpp_open, mpp_close, axistype, fieldtype, mpp_write, mpp_write_meta
  use mpp_io_mod,       only : MPP_OVERWR, MPP_NETCDF, MPP_MULTI, MPP_RDONLY, MPP_ASCII, MPP_SINGLE 
  use mpp_io_mod,       only : mpp_get_atts, mpp_get_axes, mpp_get_axis_data, mpp_get_fields
  use mpp_io_mod,       only : mpp_read, mpp_get_info
  use fms_mod,          only : open_namelist_file, check_nml_error, stdout, stdlog, close_file
  use fms_mod,          only : file_exist,  write_version_number, error_mesg, write_data
  use fms_mod,          only : fms_init, fms_end
  use axis_utils_mod,   only : get_axis_cart, get_axis_bounds
  use constants_mod,    only : PI, radius, constants_init
  use fms_io_mod,       only : fms_io_exit

  implicit none

  !--- namelist interface ----------------------------------------------
  !<NAMELIST NAME="river_regrid_nml">
  !  <DATA NAME="river_input_file" TYPE="character(len=128)" >
  !    river data source file.
  !  </DATA>
  !  <DATA NAME="land_grid_file" TYPE="character(len=128)" >
  !    the grid file that contains land grid information.
  !  </DATA>
  !  <DATA NAME="river_output_file" TYPE="character(len=128)" >
  !    The output river data file after coupled with land grid.
  !  </DATA>
  !</NAMELIST>
  character(len=128) :: river_input_file  = 'river_input.nc'
  character(len=128) :: land_grid_file    = 'grid_spec.nc'
  character(len=128) :: river_output_file = 'river_output.nc'
!  character(len=128) :: edit_table       = 'river_table'

  namelist /river_regrid_nml/ river_input_file, river_output_file, land_grid_file !, edit_table

  real :: missing_value = -9999.

  !--- derived data type ------------------------------------------------
  type river_regrid_type
     real, dimension(:),      pointer :: lon         => NULL()   ! longitude (in degree)
     real, dimension(:),      pointer :: lat         => NULL()   ! lattitude (in degree)
     real, dimension(:),      pointer :: lonb        => NULL()   ! longitude edges (in degree)
     real, dimension(:),      pointer :: latb        => NULL()   ! latitude edges (in degree)
     real, dimension(:,:),    pointer :: cellarea    => NULL()   ! cell area
     real, dimension(:,:),    pointer :: subL        => NULL()   ! subbasin length
     real, dimension(:,:),    pointer :: subA        => NULL()   ! subbasin area
     real, dimension(:,:),    pointer :: disttomouth => NULL()   ! distant to mouth
     real, dimension(:,:),    pointer :: disttoocean => NULL()   ! distant to ocean
     real, dimension(:,:),    pointer :: celllength  => NULL()   ! cell size
     integer, dimension(:,:), pointer :: order       => NULL()   ! river order
     integer, dimension(:,:), pointer :: basinid     => NULL()   ! basin id
     integer, dimension(:,:), pointer :: tocell      => NULL()   ! tocell information
     integer, dimension(:,:), pointer :: fromcell    => NULL()   ! fromcell information
     integer, dimension(:,:), pointer :: basincells  => NULL()   ! basincells
     integer, dimension(:,:), pointer :: travel      => NULL()   ! number of travel to ocean
     real, dimension(:,:),    pointer :: mask        => NULL()   ! land/sea mask(land=1)
     real, dimension(:,:,:),  pointer :: dist        => NULL()   ! distance to neighbor cell
     integer, dimension(:,:), pointer :: neighbor    => NULL()   ! neighbor information
     integer, dimension(:,:), pointer :: land_order  => NULL()   ! land_order
     integer                          :: nlon, nlat              ! grid size
  end type river_regrid_type

  !--- version information ---------------------------------------------
  character(len=128) :: version = '$ID$'
  character(len=128) :: tagname = '$Name$'

  !--- other variables
  integer                           :: nlon_lnd, nlat_lnd
  real, dimension(:),   allocatable :: lonb_lnd, latb_lnd
  real, dimension(:,:), allocatable :: lnd_mask
  type(river_regrid_type),save      :: Source
  type(river_regrid_type),save      :: River
  real, parameter                   :: lnd_thresh = 0.1
  real, parameter                   :: epsln      = 1.0e-10  ! some small number
  logical                           :: module_is_initialized = .FALSE.
  integer                           :: next_basinid = 0
  real                              :: D2R
  integer                           :: nlon, nlat  ! grid size of extended river grid

  !--- begin the program

  call fms_init
  call constants_init
  call river_regrid_init()
  call river_regrid_data(River)
  call write_river_regrid_data(River)
  call river_regrid_end(River)
  call fms_end

contains

  !--- initialization routine, reading namelist, read river source 
  !--- data and land grid information.
  subroutine river_regrid_init

    integer :: unit, ierr, io_status

    !--- read namelist and write out namelist --------------------------
    unit = open_namelist_file()
    read  (unit, river_regrid_nml,iostat=io_status)
    write (stdout(),'(/)')
    write (stdout(), river_regrid_nml)
    write (stdlog(), river_regrid_nml)
    ierr = check_nml_error(io_status,'river_regrid_nml')
    call close_file (unit)

    D2R = PI/180.0

    !--- write out version information ---------------------------------
    call write_version_number( version, tagname )

    !--- read the river data from river_input_file
    call read_river_src_data( )

    !--- read land grid from land_grid_file
    call get_land_grid_data( )

    return

  end subroutine river_regrid_init

  !#####################################################################

  !--- remap and define river routing network.
  subroutine river_regrid_data(River)
    type(river_regrid_type), intent(out) :: River

    !--- define River to make land and river has the same range and initialize
    !--- the river data.
    call define_new_grid(River)

    !--- remove previously land point
    call trim_river(River)

    !--- define land order for new river grid
    call define_new_land_order(River)  

    !--- define new land point data
    call define_new_land_data(River)

  end subroutine river_regrid_data

  !#####################################################################
  !--- write out river routing data
  subroutine write_river_regrid_data(River)
    type(river_regrid_type), intent(in) :: River

    integer                    :: unit, i, j
    type(axistype)             :: axis_x, axis_y, axis_xt, axis_yt, axis_xb, axis_yb
    type(fieldtype)            :: fld_basin, fld_disttomouth, fld_basincells, fld_cellarea
    type(fieldtype)            :: fld_disttoocean, fld_fromcell, fld_order, fld_subA
    type(fieldtype)            :: fld_subL, fld_tocell, fld_travel, fld_celllength, fld_mask
    type(fieldtype)            :: fld_neighbor, fld_land_order, fld_land_mask
    real, allocatable          :: tmp(:,:), tmp1d(:)

    call mpp_open(unit,trim(river_output_file), MPP_OVERWR, MPP_NETCDF, threading = MPP_MULTI)

    !--- write out axis meta data ---------------------------------------------
    call mpp_write_meta(unit, axis_x,'lon','degree_east','Nominal Longitude of T-cell center', &
         cartesian ='X', data = River%lon )

    call mpp_write_meta(unit, axis_y,'lat','degree_north','Nominal Latitude of T-cell center', &
         cartesian ='Y', data = River%lat )

    allocate(tmp1d(nlon_lnd))
    do i = 1, nlon_lnd
       tmp1d(i) = 0.5*(lonb_lnd(i) + lonb_lnd(i+1) ) 
    enddo

    call mpp_write_meta(unit, axis_xt, 'lont_land', 'degree_east',  &
      'Nominal Longitude of land grid center', cartesian ='X', data = tmp1d )

    deallocate(tmp1d)
    allocate(tmp1d(nlat_lnd))
    do j = 1, nlat_lnd
       tmp1d(j) = 0.5*(latb_lnd(j) + latb_lnd(j+1) ) 
    enddo

    call mpp_write_meta(unit, axis_yt, 'latt_land', 'degree_north',  &
      'Nominal Latitude of land grid center', cartesian ='Y', data = tmp1d )
    deallocate(tmp1d)

    call mpp_write_meta(unit, axis_xb, 'lonb_land', 'degree_east',  &
      'Nominal Longitude of land grid edge', cartesian ='X', data = lonb_lnd )

    call mpp_write_meta(unit, axis_yb, 'latb_land', 'degree_north',  &
      'Nominal Latitude of land grid edge', cartesian ='Y', data = latb_lnd )

    !--- write out field meta data ---------------------------------------------
    call mpp_write_meta(unit, fld_basin, (/axis_x, axis_y/), 'basin', &
         'none','river basin id', pack=1, missing=missing_value )

    call mpp_write_meta(unit, fld_disttomouth, (/axis_x, axis_y/), 'disttomouth', &
         'm','distance to mouth', pack=1, missing=missing_value)

    call mpp_write_meta(unit, fld_basincells, (/axis_x, axis_y/), 'basincells', &
         'none','number of cells in subbasin', pack=1, missing=missing_value)

    call mpp_write_meta(unit, fld_cellarea, (/axis_x, axis_y/), 'cellarea', &
         'm2','cell area', pack=1, missing=missing_value)

    call mpp_write_meta(unit, fld_disttoocean, (/axis_x, axis_y/), 'disttoocean', &
         'm','distance to ocean', pack=1, missing=missing_value)

    call mpp_write_meta(unit, fld_fromcell, (/axis_x, axis_y/), 'fromcell', &
         'none','sum of directions from upstream cells', pack=1, missing=missing_value)

    call mpp_write_meta(unit, fld_order, (/axis_x, axis_y/), 'order', &
         'none','Strahler stream order', pack=1, missing=missing_value)

    call mpp_write_meta(unit, fld_subA, (/axis_x, axis_y/), 'subA', &
         'm2','subbasin area', pack=1, missing=missing_value)

    call mpp_write_meta(unit, fld_subL, (/axis_x, axis_y/), 'subL', &
         'm','subbasin length', pack=1, missing=missing_value)

    call mpp_write_meta(unit, fld_tocell, (/axis_x, axis_y/), 'tocell', &
         'none','direction to downstream cell', pack=1, missing=missing_value)

    call mpp_write_meta(unit, fld_travel, (/axis_x, axis_y/), 'travel', &
         'none','cells left to travel before reaching ocean', pack=1, missing=missing_value)

    call mpp_write_meta(unit, fld_celllength, (/axis_x, axis_y/), 'celllength', &
         'm','cell length', pack=1, missing=missing_value)
    call mpp_write_meta(unit, fld_mask, (/axis_x, axis_y/), 'mask', &
         'none','land/sea mask(land = 1)', pack=1, missing=missing_value)
    call mpp_write_meta(unit, fld_neighbor, (/axis_x, axis_y/), 'neighbor', &
         'none','neighbor point', pack=1, missing=missing_value)
    call mpp_write_meta(unit, fld_land_order, (/axis_x, axis_y/), 'land_order', &
         'none','neighbor point', pack=1, missing=missing_value)

    call mpp_write_meta(unit, fld_land_mask, (/axis_xt, axis_yt/), 'mask_land', &
         'none','land/sea mask of land grid', pack=1)

    call mpp_write(unit, axis_x)
    call mpp_write(unit, axis_y)
    call mpp_write(unit, axis_xt)
    call mpp_write(unit, axis_yt)
    call mpp_write(unit, axis_xb)
    call mpp_write(unit, axis_yb)
    allocate(tmp(nlon, nlat) )
    tmp = River%basinid(1:nlon,1:nlat)
    call mpp_write(unit, fld_basin, tmp)
    call mpp_write(unit, fld_disttomouth, River%disttomouth(1:nlon,1:nlat))
    tmp = River%basincells(1:nlon,1:nlat)
    call mpp_write(unit, fld_basincells, tmp)
    call mpp_write(unit, fld_cellarea, River%cellarea(1:nlon,1:nlat))
    call mpp_write(unit, fld_disttoocean, River%disttoocean(1:nlon,1:nlat))
    tmp = River%fromcell(1:nlon,1:nlat)
    call mpp_write(unit, fld_fromcell, tmp)
    tmp = River%order(1:nlon,1:nlat)
    call mpp_write(unit, fld_order, tmp)
    call mpp_write(unit, fld_subA, River%subA(1:nlon,1:nlat))
    call mpp_write(unit, fld_subL, River%subL(1:nlon,1:nlat))
    tmp = River%tocell(1:nlon,1:nlat)
    call mpp_write(unit, fld_tocell, tmp)
    tmp = River%travel(1:nlon,1:nlat)
    call mpp_write(unit, fld_travel, tmp)
    call mpp_write(unit, fld_celllength, River%celllength(1:nlon,1:nlat))
    call mpp_write(unit, fld_mask, River%mask(1:nlon,1:nlat))
    tmp = River%neighbor(1:nlon,1:nlat)
    call mpp_write(unit, fld_neighbor, tmp)
    tmp = River%land_order(1:nlon,1:nlat)
    call mpp_write(unit, fld_land_order, tmp)
    call mpp_write(unit, fld_land_mask, lnd_mask(1:nlon_lnd,1:nlat_lnd) )
    call mpp_close(unit)
    deallocate(tmp)

  end subroutine write_river_regrid_data

  !#####################################################################
  ! release memory
  subroutine river_regrid_end(River)
    type(river_regrid_type), intent(inout) :: River

    module_is_initialized = .FALSE.
    deallocate(River%lonb, River%latb, River%cellarea, River%subL, River%subA )
    deallocate(River%disttomouth, River%disttoocean, River%celllength, River%order )
    deallocate(River%basinid, River%tocell, River%fromcell, River%basincells )
    deallocate(River%travel, River%mask, River%neighbor,  River%land_order )
    deallocate(River%dist, River%lon, River%lat )

  end subroutine river_regrid_end

  !#####################################################################
  !--- read river source data
  subroutine read_river_src_data

    type(axistype), dimension(:),  allocatable :: axes, axes_bounds
    type(fieldtype), dimension(:), allocatable :: fields
    real, dimension(:,:),          allocatable :: tmp
    real               :: old_missing_value
    integer            :: unit, i, j, n, ndim, nvar, natt, ntime, len_axes, nlon_src, nlat_src
    character(len=128) :: name
    character(len=1)   :: cart
    logical            :: found_cellarea, found_tocell, found_fromcell, found_basin
    logical            :: found_basincells, found_order, found_travel, found_subL
    logical            :: found_subA, found_disttomouth, found_disttoocean, found_dx

    write(stdout(),*) ' read river src data from file '//trim(river_input_file)

    if(.not. file_exist(trim(river_input_file))) &
         call mpp_error(FATAL, 'river_regrid: file '//trim(river_input_file)//' does not exist')

    !--- get the grid size and grids
    call mpp_open(unit,trim(river_input_file),action=MPP_RDONLY,form=MPP_NETCDF)
    call mpp_get_info(unit, ndim, nvar, natt, ntime)
    allocate(axes(ndim), axes_bounds(ndim), fields(nvar))
    call mpp_get_axes(unit, axes)     
    do n = 1, ndim
       cart = 'N'
       call get_axis_cart(axes(n), cart)
       call mpp_get_atts(axes(n),len=len_axes)
       if (cart == 'N') call mpp_error(FATAL,'river_routing_mod: couldnt recognize axis atts.')
       select case (cart)
       case ('X')
          nlon_src = len_axes
          allocate(Source%lonb(nlon_src+1), Source%lon(nlon_src))
          call mpp_get_axis_data(axes(n), Source%lon)
          call get_axis_bounds(axes(n), axes_bounds(n), axes)
          call mpp_get_axis_data(axes_bounds(n), Source%lonb)
       case ('Y')
          nlat_src = len_axes
          allocate(Source%latb(nlat_src+1), Source%lat(nlat_src))
          call mpp_get_axis_data(axes(n), Source%lat)
          call get_axis_bounds(axes(n), axes_bounds(n), axes)
          call mpp_get_axis_data(axes_bounds(n), Source%latb)
       end select
    enddo

    Source%nlon = nlon_src
    Source%nlat = nlat_src
    !--- read the field
    allocate(Source%cellarea   (nlon_src,nlat_src), Source%tocell     (nlon_src,nlat_src) )
    allocate(Source%fromcell   (nlon_src,nlat_src), Source%basinid    (nlon_src,nlat_src) )
    allocate(Source%basincells (nlon_src,nlat_src), Source%travel     (nlon_src,nlat_src) )
    allocate(Source%subL  (nlon_src,nlat_src), Source%subA  (nlon_src,nlat_src) )
    allocate(Source%disttomouth(nlon_src,nlat_src), Source%disttoocean(nlon_src,nlat_src) )
    allocate(Source%celllength (nlon_src,nlat_src), Source%order      (nlon_src,nlat_src) )
    allocate(tmp               (nlon_src,nlat_src) )

    call mpp_get_fields(unit,fields)
    found_cellarea = .FALSE.    ; found_tocell = .FALSE.
    found_fromcell = .FALSE.    ; found_basin = .FALSE.
    found_basincells = .FALSE.  ; found_order = .FALSE.
    found_travel = .FALSE.      ; found_subL = .FALSE.
    found_subA = .FALSE.        ; found_disttomouth = .FALSE.
    found_disttoocean = .FALSE. ; found_dx = .FALSE.
    do n = 1, nvar
       call mpp_get_atts(fields(n),name=name, missing = old_missing_value)
       select case (trim(name))
       case ('cellarea')
          found_cellarea = .TRUE.
          call mpp_read(unit,fields(n), Source%cellarea )
          where(Source%cellarea == old_missing_value) Source%cellarea = missing_value
       case ('tocell')
          found_tocell = .TRUE.
          call mpp_read(unit,fields(n), tmp )
          Source%tocell = tmp               
          where(Source%tocell == old_missing_value) Source%tocell = missing_value
       case ('fromcell')
          found_fromcell = .TRUE.
          call mpp_read(unit,fields(n), tmp ) 
          Source%fromcell = tmp
          where(Source%fromcell == old_missing_value) Source%fromcell = missing_value
       case ('basin')
          found_basin = .TRUE.
          call mpp_read(unit,fields(n), tmp ) 
          Source%basinid = tmp
          where(Source%basinid == old_missing_value) Source%basinid = missing_value
       case ('basincells')
          found_basincells = .TRUE.
          call mpp_read(unit,fields(n), tmp) 
          Source%basincells = tmp 
          where(Source%basincells == old_missing_value) Source%basincells = missing_value           
       case ('order')
          found_order = .TRUE.
          call mpp_read(unit,fields(n), tmp ) 
          Source%order = tmp
          where(Source%order == old_missing_value) Source%order = missing_value
       case ('travel')
          found_travel = .TRUE.
          call mpp_read(unit,fields(n), tmp )
          Source%travel = tmp
          where(Source%travel == old_missing_value) Source%travel = missing_value
       case ('subL')
          found_subL = .TRUE.
          call mpp_read(unit,fields(n), Source%subL ) 
          where(Source%subL == old_missing_value) Source%subL = missing_value
       case ('subA')
          found_subA = .TRUE.
          call mpp_read(unit,fields(n), Source%subA ) 
          where(Source%subA == old_missing_value) Source%subA = missing_value
       case ('disttomouth')
          found_disttomouth = .TRUE.
          call mpp_read(unit,fields(n), Source%disttomouth ) 
          where(Source%disttomouth == old_missing_value) Source%disttomouth = missing_value
       case ('disttoocean')
          found_disttoocean = .TRUE.
          call mpp_read(unit,fields(n), Source%disttoocean ) 
          where(Source%disttoocean == old_missing_value) Source%disttoocean = missing_value
       case ('dx')
          found_dx = .TRUE.
          call mpp_read(unit,fields(n), Source%celllength ) 
          where(Source%celllength == old_missing_value) Source%celllength = missing_value
       end select
    enddo
    deallocate(tmp, axes, axes_bounds, fields )

    if(.not. found_cellarea) call mpp_error(FATAL,'river_routing_mod: cellarea is not in file '//trim(river_input_file) )
    if(.not. found_tocell) call mpp_error(FATAL,'river_routing_mod: tocell is not in file '//trim(river_input_file) )
    if(.not. found_fromcell) call mpp_error(FATAL,'river_routing_mod: fromcell is not in file '//trim(river_input_file) )
    if(.not. found_basin) call mpp_error(FATAL,'river_routing_mod: basin is not in file '//trim(river_input_file) )
    if(.not. found_basincells) call mpp_error(FATAL,'river_routing_mod: basincells is not in file '//trim(river_input_file) )
    if(.not. found_order) call mpp_error(FATAL,'river_routing_mod: order is not in file '//trim(river_input_file) )
    if(.not. found_travel) call mpp_error(FATAL,'river_routing_mod: travel is not in file '//trim(river_input_file) )
    if(.not. found_subL) call mpp_error(FATAL,'river_routing_mod: subL is not in file '//trim(river_input_file) )
    if(.not. found_subA) call mpp_error(FATAL,'river_routing_mod: subA is not in file '//trim(river_input_file) )
    if(.not. found_disttomouth) call mpp_error(FATAL,'river_routing_mod: disttomouth is not in file '//trim(river_input_file) )
    if(.not. found_disttoocean) call mpp_error(FATAL,'river_routing_mod: disttoocean is not in file '//trim(river_input_file) )
    if(.not. found_dx) call mpp_error(FATAL,'river_routing_mod: dx is not in file '//trim(river_input_file) )

    call mpp_close(unit)

  end subroutine read_river_src_data

  !#####################################################################
  !--- read land grid information
  subroutine get_land_grid_data

    integer                      :: unit, ndim, nvar, natt, ntime, len, i, j, n
    type(axistype), allocatable  :: axes(:)
    type(fieldtype), allocatable :: fields(:)
    real,            allocatable :: area_lnd_cell(:,:), area_lnd(:,:)
    character(len=128)           :: name
    logical                      :: found_xbl, found_ybl, found_area_lnd_cell, found_area_lnd

    write(stdout(),*) 'read land grid information from file '//trim(land_grid_file)
    if(.not. file_exist(trim(land_grid_file))) &
         call mpp_error(FATAL, 'river_regrid: file '//trim(land_grid_file)//' does not exist')

    !--- get the grid size and grids
    call mpp_open(unit,trim(land_grid_file),action=MPP_RDONLY,form=MPP_NETCDF)   

    call mpp_get_info(unit, ndim, nvar, natt, ntime)
    allocate(axes(ndim), fields(nvar))
    call mpp_get_axes(unit,axes)

    found_xbl = .false. ; found_ybl = .false.
    do n=1, ndim
       call mpp_get_atts(axes(n),name=name,len=len)
       select case (trim(name))
       case ('xbl')
          found_xbl = .true.
          nlon_lnd = len - 1
          allocate(lonb_lnd(len))
          call mpp_get_axis_data(axes(n),data=lonb_lnd)
       case ('ybl')
          found_ybl = .true.
          nlat_lnd = len - 1
          allocate(latb_lnd(len))
          call mpp_get_axis_data(axes(n),data=latb_lnd)
       end select
    enddo

    if(.not. found_xbl) call mpp_error(FATAL,'xbl is not in the file '//trim(land_grid_file) )
    if(.not. found_ybl) call mpp_error(FATAL,'ybl is not in the file '//trim(land_grid_file) )

    !--- get the land/sea mask
    allocate( area_lnd_cell(nlon_lnd,nlat_lnd), area_lnd(nlon_lnd, nlat_lnd), lnd_mask(-1:nlon_lnd+2,-1:nlat_lnd+2))
    call mpp_get_fields(unit,fields)
    found_area_lnd_cell = .false.; found_area_lnd = .false.
    do n = 1, nvar
       call mpp_get_atts(fields(n),name=name)
       select case (trim(name))
       case ('AREA_LND_CELL')
          found_area_lnd_cell = .true.
          call mpp_read(unit,fields(n), area_lnd_cell )
       case ('AREA_LND')
          found_area_lnd = .true. 
          call mpp_read(unit,fields(n), area_lnd )
       end select
    enddo

    if(.not. found_area_lnd_cell) call mpp_error(FATAL,'AREA_LND_CELL is not in the file '//trim(land_grid_file) )
    if(.not. found_area_lnd)      call mpp_error(FATAL,'AREA_LND is not in the file '//trim(land_grid_file) )

    lnd_mask = 0.0
    lnd_mask(1:nlon_lnd,1:nlat_lnd) = area_lnd(:,:)/area_lnd_cell(:,:)
    call update_halo(lnd_mask(-1:nlon_lnd+2,1:nlat_lnd),2)
    call mpp_close(unit)
    deallocate(axes, fields, area_lnd_cell, area_lnd )

  end subroutine get_land_grid_data

  !#####################################################################
  !--- extend river grid to match land grid
  !--- define land/sea mask for each grid.
  subroutine define_new_grid(River)
    type(river_regrid_type), intent(out) :: River

    integer :: num_lon_begin, num_lon_end, num_lat_begin, num_lat_end, nlon_src, nlat_src
    integer :: i, j, m, n, i_trans, j_trans, npes, im1, ip1, jm1, jp1, layout(2) = (/1,0/)
    integer :: i_1, j_1, i_2, j_2, num_missing_nbr, nlon_lnd, nlat_lnd
    real    :: dx, dy

    write(stdout(),*) 'define new grid After coupled with land grid'
    !--- We assume river source grid always is uniform
    dx = Source%lonb(2) - Source%lonb(1)
    dy = Source%latb(2) - Source%latb(1)  
    nlon_src = Source%nlon; nlat_src = Source%nlat
    nlon_lnd = size(lonb_lnd); nlat_lnd = size(latb_lnd)

    !--- The river grid should be inside the land grid
    if(Source%lonb(1) .lt. lonb_lnd(1) .or. Source%lonb(nlon_src+1) .gt. lonb_lnd(nlon_lnd)  .or. &
         Source%latb(1) .lt. latb_lnd(1) .or. Source%latb(nlat_src+1) .gt. latb_lnd(nlat_lnd))    &
         call mpp_error(FATAL,'river_regrid: static river grid is not inside land grid')

    !--- find how many points need to be added at the begining and end
    call find_grid_addition(lonb_lnd(1), lonb_lnd(nlon_lnd), Source%lonb(1), &
         Source%lonb(size(Source%lonb)), dx, num_lon_begin, num_lon_end)

    call find_grid_addition(latb_lnd(1), latb_lnd(nlat_lnd), Source%latb(1), &
         Source%latb(size(Source%latb)), dy, num_lat_begin, num_lat_end)  

    River%nlon = nlon_src + num_lon_begin + num_lon_end
    River%nlat = nlat_src + num_lat_begin + num_lat_end
    nlon = River%nlon
    nlat = River%nlat

    !--- memory allocation 
    allocate(River%lonb(nlon+1), River%latb(nlat+1), River%cellarea(nlon,nlat) )
    allocate(River%subL(0:nlon+1,0:nlat+1), River%subA(0:nlon+1,0:nlat+1) )
    allocate(River%disttomouth(nlon,nlat), River%disttoocean(nlon,nlat) )
    allocate(River%celllength(nlon,nlat), River%order(0:nlon+1,0:nlat+1) )
    allocate(River%basinid(nlon,nlat), River%tocell(0:nlon+1,0:nlat+1) )
    allocate(River%fromcell(nlon,nlat), River%basincells(nlon,nlat) )
    allocate(River%travel(nlon,nlat), River%mask(-1:nlon+2,-1:nlat+2) )
    allocate(River%neighbor(nlon,nlat),  River%land_order(nlon,nlat) )
    allocate(River%dist(nlon, nlat, 8), River%lon(nlon), River%lat(nlat) )

    !--- get the expanded river grid
    call grid_expansion(River%lonb, Source%lonb, lonb_lnd(1), &
         lonb_lnd(nlon_lnd), dx, num_lon_begin, num_lon_end)
    call grid_expansion(River%latb, Source%latb, latb_lnd(1), &
         latb_lnd(nlat_lnd), dy, num_lat_begin, num_lat_end)

    !--- each river grid box should be inside a land grid box.
    do j = 1, nlat
       do n = 1, nlat_lnd
          if(River%latb(j) .lt. latb_lnd(n)- epsln .and. River%latb(j+1) .gt. latb_lnd(n) + epsln) then
             write(stdout(),*) 'river grid box at j= ',j,' with latitude bound = ',  River%latb(j:j+1), &
                  ' is accross land grid box at n = ', n, ' with latitude = ', latb_lnd(n)
             call mpp_error(FATAL,'river_regrid: Not all the river grid box lies inside a land grid boc')
          endif
       enddo
    enddo

    do i =1, nlon
       do m = 1, nlon_lnd
          if(River%lonb(i) .lt. lonb_lnd(m)- epsln  .and. River%lonb(i+1) .gt. lonb_lnd(m) + epsln) then
             write(stdout(),*) 'river grid box at i= ',i,' with lontitude bound = ',  River%lonb(i:i+1), &
                  ' is accross land grid box at m = ', m, ' with lontitude = ', lonb_lnd(m)
             call mpp_error(FATAL,'river_regrid: Not all the river grid box lies inside a land grid boc')
          endif
       enddo
    enddo

    write(stdout(),*) num_lon_begin, 'points is added at the west bound.  ',  num_lon_end, 'points is added at the east bound'
    write(stdout(),*) num_lat_begin, 'points is added at the south bound. ',  num_lat_end, 'points is added at the north bound'

    do i = 1, nlon
       River%lon(i) = ( River%lonb(i) + River%lonb(i+1) ) * 0.5
    enddo
    do j = 1, nlat
       River%lat(j) = ( River%latb(j) + River%latb(j+1) ) * 0.5
    enddo

    !--- fill it with nan value
    River%cellarea    = missing_value
    River%subL   = missing_value  
    River%subA   = missing_value
    River%disttomouth = missing_value  
    River%disttoocean = missing_value
    River%celllength  = missing_value
    River%order       = missing_value
    River%basinid     = missing_value
    River%tocell      = missing_value
    River%fromcell    = missing_value
    River%basincells  = missing_value
    River%travel      = missing_value

    !--- file it with orginal grid data
    do j = 1, nlat_src
       do i = 1, nlon_src
          i_trans = num_lon_begin + i
          j_trans = num_lat_begin + j
          River%cellarea(i_trans,j_trans)    = Source%cellarea(i,j)
          River%subL(i_trans,j_trans)   = Source%subL(i,j)  
          River%subA(i_trans,j_trans)   = Source%subA(i,j)
          River%disttomouth(i_trans,j_trans) = Source%disttomouth(i,j)  
          River%disttoocean(i_trans,j_trans) = Source%disttoocean(i,j)
          River%celllength(i_trans,j_trans)  = Source%celllength(i,j)
          River%order(i_trans,j_trans)       = Source%order(i,j)
          River%basinid(i_trans,j_trans)     = Source%basinid(i,j)
          River%tocell(i_trans,j_trans)      = Source%tocell(i,j)
          River%fromcell(i_trans,j_trans)    = Source%fromcell(i,j)
          River%basincells(i_trans,j_trans)  = Source%basincells(i,j)
          River%travel(i_trans,j_trans)      = Source%travel(i,j)
       enddo
    enddo

    !--- we assume the cyclic condition
    River%order(0,1:nlat) = River%order(nlon,1:nlat)
    River%order(nlon+1,1:nlat) = River%order(1,1:nlat)
    River%subL(0,1:nlat) = River%subL(nlon,1:nlat)
    River%subL(nlon+1,1:nlat) = River%subL(1,1:nlat)
    River%subA(0,1:nlat) = River%subA(nlon,1:nlat)
    River%subA(nlon+1,1:nlat) = River%subA(1,1:nlat)

    !--- define the mask of river grid according to land grid mask.

    call define_river_mask(River)

    !--- here we suppose the cyclic condition exists
    call update_halo(River%mask(-1:nlon+2,1:nlat), 2)

    River%tocell(0,1:nlat)      = River%tocell(nlon,1:nlat)
    River%tocell(nlon+1,1:nlat) = River%tocell(1,1:nlat)

    !--- define the distance to the neighbor at each point
    River%dist = 1.0e20      
    do j = 1, nlat
       do i = 1, nlon
          ip1 = i + 1; im1 = i - 1
          jp1 = j + 1; jm1 = j - 1
          if(im1 == 0) im1 = nlon
          if(ip1 == nlon+1) ip1 = 1
          River%dist(i,j,1) = dist(River%lon(i), River%lat(j), River%lon(ip1), River%lat(j) )
          if(jm1 .ne. 0 ) then
             River%dist(i,j,2) = dist(River%lon(i), River%lat(j), River%lon(ip1), River%lat(jm1) )
             River%dist(i,j,3) = dist(River%lon(i), River%lat(j), River%lon(i  ), River%lat(jm1) )
             River%dist(i,j,4) = dist(River%lon(i), River%lat(j), River%lon(im1), River%lat(jm1) )
          endif
          River%dist(i,j,5) = dist(River%lon(i), River%lat(j), River%lon(im1), River%lat(j) )
          if(jp1 .ne. nlat+1) then
             River%dist(i,j,6) = dist(River%lon(i), River%lat(j), River%lon(im1), River%lat(jp1) )
             River%dist(i,j,7) = dist(River%lon(i), River%lat(j), River%lon(i  ), River%lat(jp1) )
             River%dist(i,j,8) = dist(River%lon(i), River%lat(j), River%lon(ip1), River%lat(jp1) )
          endif
       enddo
    enddo

    !--- ??? next_basinid should be the maximum basinid + 1
    next_basinid = maxval(Source%basinid)

  end subroutine define_new_grid

  !#####################################################################
  !--- find number of points needed to be extended.
  subroutine find_grid_addition(lnd_start, lnd_end, river_start, river_end, dx, num_begin, num_end)
    real, intent(in)     :: lnd_start, lnd_end, river_start, river_end, dx
    integer, intent(out) :: num_begin, num_end
    real                 :: diff, min_size 

    min_size =  0.05
    if(dx .lt. min_size) call mpp_error(FATAL,'river_routing_mod: river grid size is too small')
    !---first check how many points we need to add at the begining.
    num_begin = 0
    diff = river_start - lnd_start
    if( diff .gt. min_size) then
       num_begin = int(diff / dx)
       if(diff-num_begin*dx .gt. min_size) num_begin = num_begin + 1
    endif

    num_end = 0
    diff = lnd_end - river_end
    if( diff .gt. min_size) then
       num_end = int(diff / dx)
       if(diff-num_end*dx .gt. min_size) num_end = num_end + 1
    endif

  end subroutine find_grid_addition

  !#####################################################################
  !--- define new river grid
  subroutine grid_expansion(river_grid, orig_grid, lnd_start,lnd_end, dx, num_begin, num_end)
    real, dimension(:),    intent(out) :: river_grid
    real, dimension(:),     intent(in) :: orig_grid
    real,                   intent(in) :: lnd_start, lnd_end, dx
    integer,                intent(in) :: num_begin, num_end
    integer                            :: num_orig, i

    num_orig = size(orig_grid)

    river_grid(num_begin+1:num_begin+num_orig) = orig_grid(:)

    do i=num_begin,2,-1
       river_grid(i) = orig_grid(1) - (num_begin-i+1)*dx
    enddo
    if(num_begin .gt. 0) river_grid(1) = lnd_start

    do i=1,num_end-1
       river_grid(num_begin+num_orig+i) = orig_grid(num_orig)+dx*i
    enddo

    if(num_end .gt. 0) river_grid(size(river_grid)) = lnd_end

    return
  end subroutine grid_expansion

  !#####################################################################
  !--- define land order for new river grid.
  subroutine define_new_land_order(River)
    type(river_regrid_type), intent(inout) :: River

    integer, dimension(:,:), allocatable :: land_order
    integer                              :: neighbor, i, j, num_order
    logical                              :: found_new_land_pt

    write(stdout(),*) 'define land order for new river grid.'

    allocate(land_order(0:nlon+1,0:nlat+1) )

    !---set the initial value to some dummy number
    land_order = -1
    where( River%mask(1:nlon,1:nlat) .gt. lnd_thresh .and. River%travel(1:nlon,1:nlat) .eq. missing_value ) 
       land_order(1:nlon,1:nlat) = 0  ! new land point
    end where

    !--- continue if there is any new land point to check

    !--- first find the coastal points
    River%neighbor = 0
    do j = 1, nlat
       do i = 1, nlon
          if(land_order(i,j) == 0) then
             River%neighbor(i,j) = get_neighbor( River%mask(i-1:i+1,j-1:j+1) .lt. lnd_thresh )
             if(River%neighbor(i,j) .gt. 0) land_order(i,j) = 1
          endif
       enddo
    enddo

    land_order(0,1:nlat) = land_order(nlon,1:nlat)
    land_order(nlon+1,1:nlat) = land_order(1,1:nlat)

    num_order = 2
    do while(any(land_order(1:nlon,1:nlat) == 0) )
       found_new_land_pt = .FALSE.
       do j = 1, nlat
          do i = 1, nlon
             if(land_order(i,j) == 0) then
                River%neighbor(i,j) = get_neighbor(land_order(i-1:i+1,j-1:j+1) .eq. num_order-1)
                if( River%neighbor(i,j) .gt. 0 ) then
                   found_new_land_pt = .true.
                   land_order(i,j) = num_order
                   if(i == 1) land_order(nlon+1,j) = num_order
                   if(i == nlon) land_order(0,j) = num_order
                endif
             endif
          enddo
       enddo
       if(.not. found_new_land_pt ) exit
       num_order = num_order + 1
    enddo

    River%land_order(1:nlon,1:nlat) = land_order(1:nlon,1:nlat)
    deallocate(land_order)

    return
  end subroutine define_new_land_order

  !#####################################################################
  !--- define river routing data for new river grid
  subroutine define_new_land_data(River)

    type(river_regrid_type), intent(inout) :: River

    integer :: max_adj_ocean, max_land_order, num_nbrs, num_max_order, num_max_subL
    integer :: i, j, k, n, i_lo, the_basinid, ii, jj, num_max_subA
    logical :: is_neighbor(8)
    integer :: i_nbr, j_nbr, i_array(8), j_array(8), i_in(8), j_in(8), i_out(8), j_out(8)

    write(stdout(),*) 'define the new river grid data '

    !--- first define the outflow direction for new land grid
    do j = 1, nlat
       do i = 1, nlon
          if(River%land_order(i,j) .ge. 1) then
             do k = 1,8
                is_neighbor(k) = btest(River%neighbor(i,j),k-1)
             enddo
             River%tocell(i,j) = outflow_direction(is_neighbor, River%dist(i,j,:))
          endif
       enddo
    enddo

    !--- then define the inflow direction for new land grid
    do j = 1, nlat
       do i = 1, nlon
          if(River%land_order(i,j) .ge. 1) then
             River%fromcell(i,j) = inflow_direction(River%tocell(i-1:i+1,j-1:j+1))
          endif
       enddo
    enddo
    ! ****************************************************************
    ! *  Update properties of new land cells                         *
    ! *  (and attributes of other cells impacted by the new land)    *
    ! ****************************************************************

    ! First: properties that are only about each individual cell
    !        only have to define these for the new points added to the system
    do j = 1, nlat
       do i = 1, nlon
          if(River%land_order(i,j) .ge. 1) then
             River%cellarea(i,j) = area(River%lonb(i:i+1), River%latb(j:j+1) )
             call get_nbr_index(River%tocell(i,j), i, j, i_array, j_array, num_nbrs)
             if(num_nbrs .gt. 1) call mpp_error(FATAL,'River can only flow to one direction')
             if(num_nbrs .eq. 0) call mpp_error(FATAL,'river_regrid: should flow to some place')
             i_nbr = i_array(1)
             j_nbr = j_array(1)
             River%celllength(i,j) = dist(River%lon(i),River%lat(j), River%lon(i_nbr), River%lat(j_nbr) )
          endif
       enddo
    enddo
    ! Now update properties that are dependent on upstream cells (fromcell)
    !     Need to do this from the inner-most new land points downstream toward the ocean
    ! Loop: for new points, from max(River%land_order(i,j)) to 1, counting down
    max_land_order = maxval(River%land_order)
    do i_lo = max_land_order, 1, -1
       do j = 1, nlat
          do i = 1, nlon
             if (River%land_order(i,j) .eq. i_lo) then
                !--- find neighbors --------------
                call get_nbr_index(River%fromcell(i,j), i,j,i_array, j_array, num_nbrs)
                ! No neighbors feed this cell (fromcell = 0)
                if (num_nbrs .eq. 0) then
                   River%basinid(i,j) = next_basinid
                   next_basinid = next_basinid + 1
                   River%order(i,j) = 1
                   River%subL(i,j) = River%celllength(i,j)
                   River%subA(i,j) = River%cellarea(i,j)
                   River%basincells(i,j) = 1
                   ! only one neighbor feeds this cell
                else if(num_nbrs == 1  ) then
                   i_nbr = i_array(1)
                   j_nbr = j_array(1)
                   River%order(i,j) = River%order(i_nbr,j_nbr)
                   River%basinid(i,j) = River%basinid(i_nbr,j_nbr)
                   River%subL(i,j) = River%celllength(i,j) + River%subL(i_nbr,j_nbr)
                   River%subA(i,j) = River%cellarea(i,j) + River%subA(i_nbr,j_nbr)
                   River%basincells(i,j) = 1 + River%basincells(i_nbr,j_nbr)

                else       ! multiple neighbors feed this cell
                   ! one neighbor is the clear mainstem (by stream order)
                   i_in(1:num_nbrs) = i_array(1:num_nbrs)
                   j_in(1:num_nbrs) = j_array(1:num_nbrs)
                   call find_max_int ( River%order, i_in, j_in, num_nbrs, i_out, j_out, num_max_order ) 
                   if( num_max_order == 1) then
                      i_nbr = i_out(1)
                      j_nbr = j_out(1)
                      River%order(i,j) = River%order(i_nbr,j_nbr)
                      ! no neighbor is the clear mainstem (by stream order)
                   else   ! multiple neighbors feeding in have the same, max(order)) 
                      ! determine mainstem by largest subL, then largest subA, then arbitrary
                      i_in(1:num_max_order) = i_out(1:num_max_order)
                      j_in(1:num_max_order) = j_out(1:num_max_order)
                      call find_max_real( River%subL, i_in, j_in, num_max_order, i_out, j_out, num_max_subL )
                      if(num_max_subL == 1) then
                         i_nbr  = i_out(1)
                         j_nbr  = j_out(1)
                         River%order(i,j) = 1 + River%order(i_nbr,j_nbr)
                      else 
                         i_in(1:num_max_subL) = i_out(1:num_max_subL)
                         j_in(1:num_max_subL) = j_out(1:num_max_subL)
                         call find_max_real( River%subA, i_in, j_in, num_max_subL, i_out, j_out, num_max_subA )
                         if(num_max_subA == 1) then
                            i_nbr  = i_out(1)
                            j_nbr  = j_out(1)
                            River%order(i,j) = 1 + River%order(i_nbr,j_nbr)
                         else
                            i_nbr  = i_out(1)
                            j_nbr  = j_out(1)
                            River%order(i,j) = 1 + River%order(i_nbr,j_nbr)
                            write(stdout(),*) 'i,j, num_max_subA, i_out, j_out', i,j, &
                                 num_max_subA, i_out(1:num_max_subA), j_out(1:num_max_subA)
                            call mpp_error(NOTE,'river_routing_mod: that is impossible')
                         endif
                      endif
                   endif
                   ! These properties depend on only one input cell
                   River%basinid(i,j) = River%basinid(i_nbr,j_nbr)
                   River%subL(i,j) = River%celllength(i,j) + River%subL(i_nbr,j_nbr)
                   River%subA(i,j) = River%cellarea(i,j)
                   River%basincells(i,j) = 1
                   do n = 1, num_nbrs
                      River%subA(i,j)  = River%subA(i,j)  + River%subA(i_array(n),j_array(n))
                      River%basincells(i,j) = River%basincells(i,j) + River%basincells(i_array(n),j_array(n))
                   enddo
                   ! These properties depend on all the input neighbor cells 
                   !                   River%subA(i,j) = River%cellarea(i,j) + SUMOVER_NEIGHBORS(River%subA(?,?))
                   !                   River%basincells(i,j) = 1 + SUMOVER_NEIGHBORS(River%basincells(?,?))
                   ! update basinid of other neighbors of point (i,j)
                   do n = 1, num_nbrs
                      where(River%basinid == River%basinid(i_array(n),j_array(n)) )  &
                           River%basinid = River%basinid(i,j)
                   enddo
                endif
                if(i == 1) then
                   River%order(nlon+1,j) = River%order(i,j)
                   River%subL(nlon+1,j) = River%subL(i,j)
                else if(i == nlon) then
                   River%order(0,j) = River%order(i,j)
                   River%subL(0,j) = River%subL(i,j)
                endif
             endif
          enddo
       enddo
    enddo
    ! Now update properties that are dependent on downstream cells (tocell)
    ! Loop: start at the coast and move upstream
    do i_lo = 1, max_land_order
       do j = 1, nlat
          do i = 1, nlon
             if (River%land_order(i,j) .eq. i_lo) then
                River%travel(i,j) = River%land_order(i,j) - 1
                if (River%travel(i,j) .eq. 0) then
                   River%disttoocean(i,j) = River%celllength(i,j)
                   River%disttomouth(i,j) = River%celllength(i,j)
                else
                   call get_nbr_index(River%tocell(i,j), i, j, i_array, j_array, num_nbrs)
                   i_nbr = i_array(1)
                   j_nbr = j_array(1)
                   River%disttoocean(i,j) = River%celllength(i,j) + River%disttoocean(i_nbr,j_nbr)

                   if (River%order(i,j) .eq. River%order(i_nbr,j_nbr)) then
                      River%disttomouth(i,j) = River%celllength(i,j) + River%disttomouth(i_nbr,j_nbr)
                   else
                      River%disttomouth(i,j) = River%celllength(i,j) 
                   endif
                endif
             endif
          enddo
       enddo
    enddo

    ! Now must update pre-existing network (anywhere that it used to go to the ocean and now it goes to land)
    ! Only the three cell attributes which depend on downstream cells need to be updated.
    do j = 1, nlat
       do i = 1, nlon
          if ((River%travel(i,j) .eq. 0) .and. (River%land_order(i,j) .lt. 0) ) then
             call get_nbr_index(River%tocell(i,j), i, j, i_array, j_array, num_nbrs)
             i_nbr = i_array(1)
             j_nbr = j_array(1)
             if(River%land_order(i_nbr, j_nbr) .gt. 0) then
                River%travel(i,j) = River%travel(i,j) + River%travel(i_nbr, j_nbr)
                River%disttoocean(i,j) = River%disttoocean(i,j) + River%disttoocean(i_nbr, j_nbr)
                if(River%order(i,j) .EQ. River%order(i_nbr,j_nbr) )   &
                     River%disttomouth(i,j) = River%disttomouth(i,j) + River%disttomouth(i_nbr,j_nbr)
                the_basinid = River%basinid(i,j)
                do jj = 1, nlat
                   do ii = 1, nlon
                      if(River%land_order(ii,jj) .lt. 0 .and. River%basinid(ii,jj) .EQ. the_basinid) then
                         if( is_connected(River%tocell(1:nlon,1:nlat), ii, jj, i, j) )then
                            River%travel(ii,jj) = River%travel(ii,jj) + River%travel(i, j)
                            River%disttoocean(ii,jj) = River%disttoocean(ii,jj) + River%disttoocean(i, j)
                            if( River%order(ii,jj) .EQ. River%order(i,j))  &
                                 River%disttomouth(ii,jj) = River%disttomouth(ii,jj) + River%disttomouth(i,j)
                         endif
                      endif
                   enddo
                enddo
             endif
          endif
       enddo
    enddo

    return

  end subroutine define_new_land_data

  !#####################################################################
  !--- this function will return an integer of power of 2.
  function outflow_direction(neighbor, dist)
    logical, dimension(:), intent(in) :: neighbor
    real, dimension(:),    intent(in) :: dist
    integer                           :: outflow_direction

    integer :: max_count, max_adj, num_set, start(4)
    integer :: i, j, ii, n, outflow, num_min_points, index_min(8)
    logical :: keep_going
    real    :: min_dist

    !--- first find the maximum adjacent ocean points
    start = 0
    max_adj = 0
    num_set = 1
    i = 1

    if(count(neighbor) == 8) then  ! surrouding points are all ocean
       max_adj = 8
       start = 1
    else
       do while (i .le. 8)
          max_count = 0
          if(neighbor(i)) then
             keep_going = .true.
             max_count = 1
             j = i + 1
             do while(keep_going) 
                if(j .gt. 8) j = j-8
                if(neighbor(j)) then
                   max_count = max_count + 1
                else
                   keep_going = .false.
                   exit
                endif
                j = j + 1
             enddo
             if(max_count .gt. max_adj) then
                num_set = 1
                max_adj = max_count
                start(num_set) = i
             else if(max_count .eq. max_adj) then
                num_set = num_set  + 1
                start(num_set) = i
             endif
          endif
          i = i + max_count + 1
       enddo
    endif

    !--- find the outflow direction --------------------------------------
    if(num_set .eq. 1 .and. mod(max_adj,2) == 1) then
       !--- There will be a center point if max_adj is odd  --------------
       outflow = start(1) + max_count/2
       if(outflow .gt. 8) outflow = outflow - 8
    else !--- find the shortest distance point
       min_dist = dist(start(1))
       num_min_points = 1
       index_min(num_min_points) = start(num_set)
       do n = 1, num_set 
          do i = 1, max_adj - 1
             ii = start(n)+i
             if(ii .gt. 8 ) ii = ii - 8
             if(dist(ii) .lt. min_dist) then
                min_dist = dist(ii)
                num_min_points = 1
                index_min(1)   = ii
             else if(dist(ii) .eq. min_dist) then
                num_min_points = num_min_points + 1
                index_min(num_min_points) = ii
             endif
          enddo
       enddo
       if(num_min_points .eq. 1) then   
          outflow = index_min(num_min_points)
       else
          outflow = random_select(index_min(1:num_min_points))
       endif
    endif

    outflow_direction = 2 ** (outflow - 1)

    return

  end function outflow_direction

  !#####################################################################
  !--- calculate inflow direction from outflow data.
  function inflow_direction(outflow)
    integer, dimension(:,:), intent(in) :: outflow
    integer :: inflow_direction

    integer :: i, j
    i = 2;    j = 2    

    inflow_direction = 0
    if( outflow(i+1,j) == 16 )   inflow_direction = inflow_direction + 1
    if( outflow(i+1,j-1) == 32)  inflow_direction = inflow_direction + 2
    if( outflow(i,j-1) == 64)    inflow_direction = inflow_direction + 4
    if( outflow(i-1,j-1) == 128) inflow_direction = inflow_direction + 8
    if( outflow(i-1,j) == 1)     inflow_direction = inflow_direction + 16
    if( outflow(i-1,j+1) == 2)   inflow_direction = inflow_direction + 32
    if( outflow(i,j+1) == 4)     inflow_direction = inflow_direction + 64
    if( outflow(i+1,j+1) == 8)   inflow_direction = inflow_direction + 128

    return
  end function inflow_direction

  !#####################################################################
  !--- find neighbor points.
  function get_neighbor(is_neighbor )
    logical, dimension(:,:), intent(in) :: is_neighbor
    integer                             :: get_neighbor
    integer                             :: i,j

    i=2; j=2

    get_neighbor = 0
    if(is_neighbor(i+1,j))   get_neighbor = get_neighbor + 1
    if(is_neighbor(i+1,j-1)) get_neighbor = get_neighbor + 2
    if(is_neighbor(i,j-1))   get_neighbor = get_neighbor + 4
    if(is_neighbor(i-1,j-1)) get_neighbor = get_neighbor + 8
    if(is_neighbor(i-1,j))   get_neighbor = get_neighbor + 16
    if(is_neighbor(i-1,j+1)) get_neighbor = get_neighbor + 32
    if(is_neighbor(i,j+1))   get_neighbor = get_neighbor + 64
    if(is_neighbor(i+1,j+1)) get_neighbor = get_neighbor + 128

    return
  end function get_neighbor


  !#####################################################################
  ! --- this is only for rectangular grid --------------------------------------
  function area(x, y)
    real, dimension(2), intent(in) :: x,y
    real :: area

    real :: dx
    dx = (x(2)-x(1)) * D2R

    if(dx > PI)  dx = dx - 2.0*PI
    if(dx < -PI) dx = dx + 2.0*PI

    !??? Do we need to multiple radius*radius  
    area = radius*radius*dx*(sin(y(2)*D2R) - sin(y(1)*D2R))

    return

  end function area

  !#####################################################################
  !--- find the distance of any two grid.
  function dist(lon1, lat1, lon2, lat2)
    real, intent(in)           :: lon1,lat1,lon2,lat2
    real                       :: dist

    real :: s1, s2, dx

    dx = lon2 - lon1
    if(dx .gt. 180.) dx = dx - 360.
    if(dx .lt. -180.) dx = dx + 360.

    if(lon1 == lon2) then
       dist = abs(lat2 - lat1) * D2R
    else if(lat1 == lat2) then
       dist = abs(dx) * D2R * cos(lat1*D2R)
    else    ! diagonal distance
       s1 =  abs(dx)  * D2R * cos(lat1*D2R)
       s2 =  abs(lat2 - lat1) * D2R
       dist = sqrt(s1*s1+s2*s2)
    endif
    dist = radius * dist

    return
  end function dist

  !#####################################################################
  !--- find the max value of real data array
  subroutine find_max_real(data, i_in, j_in, num_in, i_out, j_out, num_out)
    real, dimension(0:,0:),  intent(in)  :: data
    integer, dimension(:), intent(in)  :: i_in, j_in
    integer,               intent(in)  :: num_in
    integer, dimension(:), intent(out) :: i_out, j_out
    integer,               intent(out) :: num_out

    integer, dimension(num_in) :: index
    integer                    :: n, m

    i_out(1) = i_in(1)
    j_out(1) = j_in(1)
    index(1) = 1   
    num_out  = 1  

    do n = 2, num_in
       if(data(i_in(n),j_in(n)) .gt. data(i_out(num_out), j_out(num_out))) then
          num_out = 1
          i_out(num_out) = i_in(n)
          j_out(num_out) = j_in(n)        
          index(num_out) = n
       else if(data(i_in(n),j_in(n)) .eq. data(i_out(num_out), j_out(num_out))) then
          num_out = num_out + 1
          i_out(num_out) = i_in(n)
          j_out(num_out) = j_in(n) 
          index(num_out) = n         
       endif
    enddo

    return
  end subroutine find_max_real

  !#####################################################################
  !--- find the max value of integer data array
  subroutine find_max_int(data, i_in, j_in, num_in, i_out, j_out, num_out)
    integer, dimension(0:,0:), intent(in)  :: data
    integer, dimension(:),   intent(in)  :: i_in, j_in
    integer,                 intent(in)  :: num_in
    integer, dimension(:),   intent(out) :: i_out, j_out
    integer,                 intent(out) :: num_out

    real, dimension(size(data,1), size(data,2) ) :: real_data
    real_data = data

    call find_max_real(real_data, i_in, j_in, num_in, i_out, j_out, num_out)

    return
  end subroutine find_max_int

  !#####################################################################
  !--- remove land point in source data but is ocean recognized by land model.
  subroutine trim_river(River)
    type(river_regrid_type), intent(inout) :: River

    integer, dimension(:,:), allocatable :: land_order
    integer :: i, j
    logical :: need_to_modify

    write(stdout(),*) 'remove previously land point ( is ocean after coupled with land grid). '

    !---set the initial value to some dummy number
    allocate(land_order(nlon, nlat) )
    land_order = -1
    !---find points that have river network info but the land model thinks are ocean points
    do j = 1, nlat
       do i = 1, nlon
          if(River%mask(i,j) .lt. lnd_thresh .and. River%travel(i,j) .ne. missing_value ) land_order(i,j) = 1  
       enddo
    enddo

    !--- continue if there is any land point in the river network to change to ocean
    do while(any(land_order == 1) )
       write(stdout(), *)' number of points land_order == 1', count(land_order == 1)
       write(stdout(), *)' number of points land_order == 1 and travel == 0', &
            count(land_order == 1.and. River%travel==0)
       if(.not.( any((land_order == 1) .and. (River%travel .eq. 0) ) .or. &
            any (land_order == 1 .and. River%travel .gt. 0 .and. River%fromcell .le. 0) )  )then
          do j = 1, nlat
             do i = 1, nlon
                if(land_order(i,j) == 1) then
                   write(stdout(),*)i,j, River%lon(i), &
                        River%lat(j), land_order(i,j), River%travel(i,j)
                endif
             enddo
          enddo
          exit
       endif
       do j =1, nlat
          do i =1, nlon
             !--- first situation: The previously land point is at coast.
             if(land_order(i,j) == 1) then
                need_to_modify = .false.
                if(River%travel(i,j) .eq. 0) then
                   call update_upstream(River,i,j)
                   need_to_modify = .true.
                   !--- second situation: The previously land point is not at coast and is start point of a river
                else if(River%travel(i,j) .gt. 0 .and. River%fromcell(i,j) .le. 0) then
                   call update_downstream(River,i,j)
                   need_to_modify = .true.
                endif
                if(need_to_modify) then
                   land_order(i,j) = -1        !--- return it to dummy value
                   !--- fill it with nan value
                   River%cellarea(i,j)    = missing_value
                   River%subL(i,j)   = missing_value  
                   River%subA(i,j)   = missing_value
                   River%disttomouth(i,j) = missing_value  
                   River%disttoocean(i,j) = missing_value
                   River%celllength(i,j)  = missing_value
                   River%order(i,j)       = missing_value
                   River%basinid(i,j)     = missing_value
                   River%tocell(i,j)      = missing_value
                   River%fromcell(i,j)    = missing_value
                   River%basincells(i,j)  = missing_value
                   River%travel(i,j)      = missing_value
                endif
             endif
          enddo
       enddo
    enddo

    do while(any(land_order == 1) )
       write(stdout(), *)' number of points land_order == 1', count(land_order == 1)
       if(.not. any(land_order == 1 .and. River%travel .gt. 0 .and. River%fromcell .gt. 0) ) &
            call mpp_error(FATAL,'infinity loop')

       do j =1, nlat
          do i =1, nlon
             !--- first situation: The previously land point is at coast.
             if(land_order(i,j) == 1) then
                if(River%travel(i,j) .gt. 0 .and. River%fromcell(i,j) .gt. 0) then
                   call update_upstream(River,i,j)

                   call update_downstream(River,i,j)
                endif
                land_order(i,j) = -1        !--- return it to dummy value
                !--- fill it with nan value
                River%cellarea(i,j)    = missing_value
                River%subL(i,j)   = missing_value  
                River%subA(i,j)   = missing_value
                River%disttomouth(i,j) = missing_value  
                River%disttoocean(i,j) = missing_value
                River%celllength(i,j)  = missing_value
                River%order(i,j)       = missing_value
                River%basinid(i,j)     = missing_value
                River%tocell(i,j)      = missing_value
                River%fromcell(i,j)    = missing_value
                River%basincells(i,j)  = missing_value
                River%travel(i,j)      = missing_value
             endif
          enddo
       enddo
    enddo
  end subroutine trim_river

  !#####################################################################
  !--- update downstream points.
  subroutine update_downstream(River,i,j)
    type(river_regrid_type), intent(inout) :: River
    integer,                 intent(in)    :: i,j
    integer :: i_array(8), j_array(8), i_in(8), j_in(8), i_out(8), j_out(8)
    integer :: i_nbr,j_nbr, ii,jj, iii, jjj, num_nbrs, i_other, j_other
    integer :: num_max_subL, num_max_subA, num_max_order, n, nn
    logical :: keep_going

    call get_nbr_index(River%tocell(i,j), i, j, i_array, j_array, num_nbrs)
    if(num_nbrs .eq. 1) then
       i_nbr = i_array(1); j_nbr = j_array(1)
       call get_nbr_index(River%fromcell(i_nbr,j_nbr), i_nbr, j_nbr, i_array, j_array, num_nbrs)

       ii = i_nbr; jj = j_nbr
       iii = i;    jjj = j
       do while (River%tocell(ii,jj) .gt. 0) 
          call get_nbr_index(River%fromcell(ii,jj), ii, jj, i_array, j_array, num_nbrs)
          if(num_nbrs == 1) then
             if(i_array(1) == i .and. j_array(1) == j) then
                River%order(ii,jj) = 1
             else
                River%order(ii,jj) =  River%order(i_array(1), j_array(1))
             endif
             River%subL(ii,jj) = River%subL(ii,jj) - River%subL(i,j)
          else 
             !--- first find the largest order
             i_in(1:num_nbrs) = i_array(1:num_nbrs)
             j_in(1:num_nbrs) = j_array(1:num_nbrs)
             call find_max_int ( River%order, i_in, j_in, num_nbrs, i_out, j_out, num_max_order ) 
             ! if (iii,jjj) is not among the points with max order, we can stop
             keep_going = .false.
             do n = 1, num_max_order
                if(i_out(n) .eq. iii .and. j_out(n) .eq. jjj) keep_going = .true.
             enddo
             if(.not. keep_going) exit

             if( num_max_order == 1) then
                River%subL(ii,jj) = River%subL(ii,jj) - River%subL(i,j)
                nn = 0
                do n = 1, num_nbrs
                   if( .not. (i_in(n) == i .and. j_in(n) == j ) ) then
                      nn = nn + 1
                      i_in(nn) = i_in(n)
                      j_in(nn) = j_in(n)
                   endif
                enddo
                call find_max_int(River%order, i_in, j_in, nn, i_out, j_out, num_max_order) 
                if(num_max_order == 1 )then
                   River%order(ii,jj) = River%order(i_out(1), j_out(1) )
                else
                   River%order(ii,jj) = River%order(i_out(1), j_out(1) ) + 1
                endif
                ! no neighbor is the clear mainstem (by stream order)
             else   ! multiple neighbors feeding in have the same, max(order)) 
                ! determine mainstem by largest subL, then largest subA, then arbitrary

                !--- if num_max_order is greater than 2, do not need to modify order.
                if(num_max_order .gt. 2) exit

                if(iii == i .and. jjj == j) then   !--- (iii, jjj ) is the point we are going to remove
                   keep_going = .false.
                   do n = 1, num_max_order
                      if( i_out(n) .eq. i .and. j_out(n) .eq. j) then
                         River%order(ii,jj) = River%order(i_out(1), j_out(1) )
                         keep_going = .true.
                      endif
                   enddo
                else                               !--- (iii, jjj ) is not the point we are going to remove
                   if(River%order(ii,jj) == River%order(i_out(1), j_out(1)) +1 ) then
                      keep_going = .false.
                   else
                      River%order(ii,jj) = River%order(i_out(1), j_out(1) ) + 1
                   endif
                endif
                if(.not. keep_going) exit

                !--- change subL.
                i_in(1:num_max_order) = i_out(1:num_max_order)
                j_in(1:num_max_order) = j_out(1:num_max_order)
                call find_max_real( River%subL, i_in, j_in, num_max_order, i_out, j_out, num_max_subL )
                if(num_max_subL == 1) then
                   if(i_out(1) .eq. iii .and. j_out(1) .eq. jjj) then
                      River%subL(ii,jj) = River%subL(ii,jj) - River%subL(i,j)
                   else
                      exit
                   endif
                else 
                   i_in(1:num_max_subL) = i_out(1:num_max_subL)
                   j_in(1:num_max_subL) = j_out(1:num_max_subL)
                   call find_max_real( River%subA, i_in, j_in, num_max_subL, i_out, j_out, num_max_subA )
                   if(num_max_subA == 1) then
                      if(i_out(1) .eq. iii .and. j_out(1) .eq. jjj) then
                         River%subL(ii,jj) = River%subL(ii,jj) - River%subL(i,j)
                      else
                         exit
                      endif
                   else
                      exit
                      call mpp_error(FATAL,'river_routing_mod: that is impossible')
                   endif
                endif
             endif
          endif
          call get_nbr_index(River%tocell(ii,jj), ii, jj, i_array, j_array, num_nbrs)                   
          iii = ii;         jjj = jj
          ii  = i_array(1); jj = j_array(1)
       enddo

       ii = i_nbr; jj = j_nbr
       do while(River%tocell(ii,jj) > 0)
          River%subA(ii,jj) = River%subA(ii,jj) - River%subA(i,j)
          River%basincells(ii,jj) = River%basincells(ii,jj) - River%basincells(i,j)                    
          call get_nbr_index(River%tocell(ii,jj), ii, jj, i_array, j_array, num_nbrs)
          if(num_nbrs .ne. 1) then
             call mpp_error(FATAL,'num_nbrs should be 1 here' )
          endif
          ii = i_array(1); jj = j_array(1)
       enddo

       if(River%tocell(i,j) .lt. 16) then
          River%fromcell(i_nbr,j_nbr) = River%fromcell(i_nbr,j_nbr) - River%tocell(i,j)*16
       else if(River%tocell(i,j) .ge. 16) then
          River%fromcell(i_nbr,j_nbr) = River%fromcell(i_nbr,j_nbr) - River%tocell(i,j)/16
       endif
    endif

  end subroutine update_downstream

  !#######################################################################
  !--- update upstream points
  subroutine update_upstream(River, i, j)
    type(river_regrid_type), intent(inout) :: River
    integer,                 intent(in)    :: i,j

    integer :: i_array(8), j_array(8)
    integer :: ii,jj, num_nbrs, n
    logical :: connect

    call get_nbr_index(River%fromcell(i,j), i, j, i_array, j_array, num_nbrs)
    if(num_nbrs == 1) then
       do jj = 1, nlat
          do ii = 1, nlon
             connect = is_connected(River%tocell(1:nlon,1:nlat), ii, jj, i, j)
             if( connect) then
                if (River%basinid(ii,jj).EQ. River%basinid(i,j) .and. (.not.(ii == i .and. jj ==j) ) ) then
                   River%travel(ii,jj) = River%travel(ii,jj)-1
                   River%disttoocean(ii,jj) = River%disttoocean(ii,jj) - River%disttoocean(i,j)
                endif

                if (River%basinid(ii,jj).EQ. River%basinid(i,j) .and. River%order(ii,jj) .EQ. River%order(i,j) &
                     .and. (.not.(ii == i .and. jj ==j) ) )then
                   River%disttomouth(ii,jj) = River%disttomouth(ii,jj) - River%disttomouth(i,j)
                endif
             endif
          enddo
       enddo
    else  if(num_nbrs > 1) then
       do n = 1, num_nbrs
          River%basinid(i_array(n),j_array(n)) = next_basinid
          River%travel(i_array(n),j_array(n)) = River%travel(i_array(n),j_array(n)) - 1
          River%disttoocean(i_array(n),j_array(n)) = &
               River%disttoocean(i_array(n),j_array(n)) - River%disttoocean(i,j)
          if( River%order(i_array(n),j_array(n)) .EQ. River%order(i,j)) then
             River%disttomouth(i_array(n),j_array(n)) = &
                  River%disttomouth(i_array(n),j_array(n)) - River%disttomouth(i,j)
          endif

          !--- Search all the points to see which point is connected with i_array(n), j_array(n)
          do jj = 1, nlat
             do ii = 1, nlon
                if(River%basinid(ii,jj) .eq. River%basinid(i,j) ) then
                   connect = is_connected(River%tocell(1:nlon,1:nlat), ii, jj, i_array(n), j_array(n))
                   if( connect) then
                      River%basinid(ii,jj) = next_basinid
                      River%travel(ii,jj) = River%travel(ii,jj) - 1
                      River%disttoocean(ii,jj) = River%disttoocean(ii,jj) - River%disttoocean(i,j)
                      if( River%order(ii,jj) .EQ. River%order(i,j)) then
                         River%disttomouth(ii,jj) = River%disttomouth(ii,jj) - River%disttomouth(i,j)
                      endif
                   endif
                endif
             enddo
          enddo
          next_basinid = next_basinid + 1
       enddo
    endif

  end subroutine update_upstream

  !#####################################################################
  !--- check if two point are connected by a river
  function is_connected(tocell, i, j, i_dst, j_dst)
    integer, dimension(:,:), intent(in) :: tocell
    integer,                 intent(in) :: i, j, i_dst, j_dst            
    logical                             :: is_connected

    integer :: cell, ii, jj, m, n, num_nbr, i_array(8), j_array(8)

    is_connected = .false.
    cell = tocell(i,j)
    ii = i; jj = j
    do while(cell .gt. 0)
       call get_nbr_index(cell, ii, jj, i_array, j_array, num_nbr)
       if(num_nbr .ne. 1) call mpp_error(FATAL,'river_regrid: num_nbr should be 1')
       m = i_array(1); n = j_array(1)
       if(m == i_dst .and. n == j_dst) then
          is_connected = .true.
          return
       else
          cell = tocell(m,n)
          ii = m; jj = n
       endif

    enddo

  end function is_connected

  !#####################################################################
  !--- random select a direction.
  function random_select(array)
    integer, dimension(:), intent(in) :: array
    integer :: random_select

    integer, save :: prev = 0
    integer :: i, num_ele
    num_ele = size(array)

    do 
       prev = prev + 1
       if(prev .gt. 8) prev = prev - 8
       do i = 1, num_ele
          if(array(i) .eq. prev) then
             random_select = prev
             return
          endif
       enddo

    enddo

  end function random_select

  !#####################################################################
  !--- get the index of neighbor point
  subroutine get_nbr_index(number, i_in, j_in, i_out, j_out, num_nbr)
    integer,               intent(in)  :: number        !number should be between 0 to 255
    integer,               intent(in)  :: i_in, j_in    
    integer, dimension(:), intent(out) :: i_out, j_out
    integer,               intent(out) :: num_nbr
    integer                            :: n

    num_nbr = 0
    n = number

    if(n .ge. 128) then
       num_nbr = num_nbr + 1
       n = n - 128
       i_out(num_nbr) = i_in+1; j_out(num_nbr) = j_in + 1
    endif
    if(n .ge. 64) then
       num_nbr = num_nbr + 1
       n = n - 64
       i_out(num_nbr) = i_in; j_out(num_nbr) = j_in + 1
    endif
    if(n .ge. 32) then
       num_nbr = num_nbr + 1
       n = n - 32
       i_out(num_nbr) = i_in-1; j_out(num_nbr) = j_in + 1
    endif
    if(n .ge. 16) then
       num_nbr = num_nbr + 1
       n = n - 16
       i_out(num_nbr) = i_in-1; j_out(num_nbr) = j_in 
    endif
    if(n .ge. 8) then
       num_nbr = num_nbr + 1
       n = n - 8
       i_out(num_nbr) = i_in-1; j_out(num_nbr) = j_in - 1
    endif
    if(n .ge. 4) then
       num_nbr = num_nbr + 1
       n = n - 4
       i_out(num_nbr) = i_in; j_out(num_nbr) = j_in - 1
    endif
    if(n .ge. 2) then
       num_nbr = num_nbr + 1
       n = n - 2
       i_out(num_nbr) = i_in+1; j_out(num_nbr) = j_in - 1
    endif
    if(n .ge. 1) then
       num_nbr = num_nbr + 1
       n = n - 1
       i_out(num_nbr) = i_in+1; j_out(num_nbr) = j_in
    endif

    do n = 1, num_nbr
       if(i_out(n) .eq. nlon+1) i_out(n) = 1
       if(i_out(n) .eq. 0)    i_out(n) = nlon
       if(j_out(n) .eq. 0 .or. j_out(n) .eq. nlat+1) &
            call mpp_error(FATAL,'river_regrid: neighbor latitude index should be between 1 and nj')
    enddo

    return
  end subroutine get_nbr_index

  !#####################################################################
  ! define river grid mask. We assume each rive grid box is inside a land grid box.
  subroutine define_river_mask(River)
    type(river_regrid_type), intent(inout) :: River
    integer                                :: nlon_lnd, nlat_lnd
    integer, dimension(:), allocatable     :: i_lnd, j_lnd     ! index in land grid of each river grid
    integer, dimension(:,:), allocatable   :: num_river        ! number of river grid in each land grid box
    integer, dimension(:,:,:), allocatable :: i_river, j_river ! index of all river grid contained in each land grid
    integer, dimension(:), allocatable     :: i2, j2           ! index of the river points to be land grid
    real,    dimension(:,:), allocatable   :: land_first_nbr   ! sum of land percentage of first neighbor points
    real,    dimension(:,:), allocatable   :: land_second_nbr  ! sum of land percentage of second neighbor points
    real,    dimension(:), allocatable     :: river_first_nbr  ! sum of land percentage of first neighbor points
    real,    dimension(:), allocatable     :: river_second_nbr ! sum of land percentage of second neighbor points
    logical, dimension(:,:), allocatable   :: is_not_defined   ! indicate if the mask of a land grid is defined or not.
    integer                                :: i, j, m, n, i1,j1, ii, jj, is,ie,js,je
    integer                                :: max_num, num, num_found, unit
    real                                   :: tmpsum, land_left
    logical                                :: start, flag
    character(len=128)                     :: txt

    nlon_lnd = size(lonb_lnd) - 1
    nlat_lnd = size(latb_lnd) - 1

    allocate(i_lnd(-1:nlon+2),       j_lnd(-1:nlat+2))
    allocate(num_river(nlon_lnd, nlat_lnd) )


    River%mask = missing_value  ! set the initial value to some dummy value

    !--- define the points that preset to be ocean in edit_table.
!    if( file_exist(trim(edit_table))) then
!       call mpp_error(NOTE, 'read predefined points from file '//trim(edit_table))

!       call mpp_open(unit,trim(edit_table),MPP_RDONLY,MPP_ASCII,threading=MPP_SINGLE,fileset=MPP_SINGLE)
!       do 
!          read(unit,'(a)',end=99,err=99) txt
!          call parse_edits(txt,is,ie,js,je,flag)
!          if(flag) then
!             write(stdout(),*)' points (',is,':',ie,',',js,':',je,') is set to be ocean'
!             do j = js,je
!                do i = is,ie
!                   if(i < 0 .or. i > nlon .or. j < 0 .or. j > nlat) call mpp_error(FATAL,'indices exceed grid bounds')
!                   River%mask(i,j) = 0.0
!                enddo
!             enddo
!          endif
!       enddo
!99     call mpp_close(unit)
!    else
!       call mpp_error(NOTE, 'river_regrid: No edit_table exists')
!    endif

    ! ---- find the index in land grid for each river grid
    do i = 1, nlon
       do m = 1, nlon_lnd
          if( River%lon(i) .ge. lonb_lnd(m) .and. River%lon(i) .le. lonb_lnd(m+1) ) then
             i_lnd(i) = m
             exit
          endif
       enddo
    enddo

    do j = 1, nlat
       do n = 1, nlat_lnd
          if( River%lat(j) .ge. latb_lnd(n) .and. River%lat(j) .le. latb_lnd(n+1) ) then
             j_lnd(j) = n
             exit
          endif
       enddo
    enddo

    !--- we assume the cyclic condition
    i_lnd(-1) = i_lnd(nlon-1)
    i_lnd(0)  = i_lnd(nlon)
    i_lnd(nlon+1) = i_lnd(1)
    i_lnd(nlon+2) = i_lnd(2)
    j_lnd(-1:0)          = missing_value
    j_lnd(nlat+1:nlat+2) = missing_value

    ! --- calculate number of river grids on each land grid
    do j = 1, nlat_lnd
       do i = 1, nlon_lnd
          num_river(i,j) = count(i_lnd(:) == i) * count(j_lnd(:) == j)
       enddo
    enddo

    max_num = maxval(num_river)

    ! --- get the index of river grids on each land grid
    allocate( i_river(nlon_lnd, nlat_lnd, max_num), j_river(nlon_lnd, nlat_lnd, max_num) )
    allocate(river_first_nbr(max_num), river_second_nbr(max_num) )
    allocate(i2(max_num), j2(max_num) )

    num_river = 0
    do j = 1, nlat
       do i = 1, nlon
          num_river(i_lnd(i),j_lnd(j)) = num_river(i_lnd(i),j_lnd(j)) + 1
          i_river(i_lnd(i),j_lnd(j),num_river(i_lnd(i),j_lnd(j))) = i
          j_river(i_lnd(i),j_lnd(j),num_river(i_lnd(i),j_lnd(j))) = j
       enddo
    enddo

    !--- define  mask of river grid that lies inside of a full land grid or a full ocean grid
    ! --- calculate the percentage of first and second neighbor land points of each partial cell land grid

    allocate(land_first_nbr(nlon_lnd, nlat_lnd), land_second_nbr(nlon_lnd, nlat_lnd))
    allocate(is_not_defined(nlon_lnd, nlat_lnd) )
    is_not_defined = .false.
    land_first_nbr  = 0
    land_second_nbr = 0

    do j =1, nlat_lnd
       do i = 1, nlon_lnd
          if( lnd_mask(i,j) .ge. 1.0-epsln) then      ! full land grid cell
             do n = 1, num_river(i,j)
                River%mask(i_river(i,j,n),j_river(i,j,n)) = 1.0
             enddo
          else if( lnd_mask(i,j) .le. epsln) then      ! full land grid cell
             do n = 1, num_river(i,j)
                River%mask(i_river(i,j,n),j_river(i,j,n)) = 0.0
             enddo
          else   ! partial land cell
             is_not_defined(i,j)  = .true.
             tmpsum               =  sum(lnd_mask(i-1:i+1,j-1:j+1))
             land_first_nbr (i,j) =  tmpsum - lnd_mask(i,j)
             land_second_nbr(i,j) =  sum(lnd_mask(i-2:i+2,j-2:j+2)) - tmpsum
          endif
       enddo
    enddo

    !--- update halo points
    call update_halo(River%mask(-1:nlon+2,1:nlat), 2)

    !--- Then define mask of river grid that lies inside of a partial land grid.

    do while ( any(is_not_defined) )

       write(stdout(),*) 'number of land points not defined is ', count(is_not_defined)
       ! --- first find the land points with maximum first and second neighbor land percentage 
       ! --- we suppose only one such point exists.
       call search_max_land( land_first_nbr, land_second_nbr, i1, j1 )

       ! --- calculate how many points is land inside the land point
       land_left = lnd_mask(i1,j1) * num_river(i1,j1)

       num = num_river(i1,j1)

       start = .true.
       do while(land_left .gt. 0.0)
          ! --- define first and second neighbor land percentage of each point
          river_first_nbr  = 0.0
          river_second_nbr = 0.0
          do n = 1, num
             i = i_river(i1,j1,n)
             j = j_river(i1,j1,n)
             if(River%mask(i,j) == missing_value) then
                do jj = j-1, j+1
                   do ii = i-1 ,i+1
                      if(River%mask(ii,jj) .ne. missing_value) then
                         river_first_nbr(n)  = river_first_nbr(n) + River%mask(ii,jj)
                      else
                         if(i_lnd(ii) .ne. missing_value .and. j_lnd(jj) .ne. missing_value) then
                            river_first_nbr(n) = river_first_nbr(n) + lnd_mask(i_lnd(ii),j_lnd(jj))
                         endif
                      endif
                   enddo
                enddo
                do jj = j-2, j+2
                   do ii = i-2, i+2
                      if(River%mask(ii,jj) .ne. missing_value) then
                         river_second_nbr(n)  = river_second_nbr(n) + River%mask(ii,jj)
                      else
                         if(i_lnd(ii) .ne. missing_value .and. j_lnd(jj) .ne. missing_value) then
                            river_second_nbr(n) = river_second_nbr(n) + lnd_mask(i_lnd(ii),j_lnd(jj))
                         endif
                      endif
                   enddo
                enddo               
             endif 
          enddo

          ! --- find the river points with maximum first and second neighbor land points
          ! --- if the land point is an isolated land point, we started from the center
          if(land_first_nbr(i1,j1) .le. epsln .and. start) then
             start = .false.
             if( mod(num,2) == 1) then
               num_found = 1
               i2(1) = i_river(i1,j1,num/2+1)
               j2(1) = j_river(i1,j1,num/2+1)                          
             else
               num_found = 2
               i2(1) = i_river(i1,j1,num/2)
               j2(1) = j_river(i1,j1,num/2)   
               i2(2) = i_river(i1,j1,num/2+1)
               j2(2) = j_river(i1,j1,num/2+1)  
             endif    
          else

             call search_max_river(river_first_nbr(1:num), river_second_nbr(1:num), i_river(i1,j1,1:num), &
                 j_river(i1,j1,1:num), i2, j2, num_found)
          endif

          if( land_left .ge. real(num_found)) then
             do m = 1, num_found
                River%mask(i2(m),j2(m)) = 1.0
                land_left = land_left - 1.0
             enddo
          else
             do m = 1, num_found
                River%mask(i2(m),j2(m)) = land_left/real(num_found)
             enddo
             land_left = 0.0
          endif
          call update_halo(River%mask(-1:nlon+2,1:nlat), 2)
       enddo

       do n =1, num
          if(River%mask(i_river(i1,j1,n),j_river(i1,j1,n)) == missing_value) then
             River%mask(i_river(i1,j1,n),j_river(i1,j1,n)) = 0.0
          endif
       enddo

       call update_halo(River%mask(-1:nlon+2,1:nlat), 2)

       is_not_defined(i1,j1)  = .false.     ! all the points inside this land grid is defined.
       land_first_nbr(i1,j1)  = 0.0
       land_second_nbr(i1,j1) = 0.0
    enddo

    if(any(River%mask(:,1:nlat) == missing_value)) call mpp_error(FATAL, 'still some river mask is not defined')

    !--- set mask to large positive value for the purpose of defining neighbor points.
    River%mask(:,-1:0) = 999.0
    River%mask(:,nlat+1:nlat+2) = 999.0

  end subroutine define_river_mask

  !#####################################################################
  !--- find land points with maximum first and second nbr land portion.
  subroutine search_max_land(first_nbr, second_nbr, i_out, j_out)
    real, dimension(:,:),   intent(in) :: first_nbr, second_nbr
    integer,               intent(out) :: i_out, j_out  
    integer, dimension(:), allocatable :: i_array1, j_array1, i_array2, j_array2 
    integer                            :: i, j, n, num, nlon_lnd, nlat_lnd
    integer                            :: num_max_first, num_max_second, num_max_mask
    real                               :: max_first_nbr, max_second_nbr, max_mask

    nlon_lnd = size(lonb_lnd) - 1
    nlat_lnd = size(latb_lnd) - 1
    max_second_nbr = -1.0
    num_max_second = 1
    max_first_nbr  = maxval(first_nbr)
    num_max_first  = count( first_nbr ==  max_first_nbr ) 
    if(num_max_first .gt. 1) then
       allocate(i_array1(num_max_first),j_array1(num_max_first)) 
       allocate(i_array2(num_max_first),j_array2(num_max_first)) 
    endif

    num = 1
    do j = 1, nlat_lnd
       do i = 1, nlon_lnd
          if(first_nbr(i,j) == max_first_nbr) then
             if(num_max_first == 1) then
                i_out = i
                j_out = j 
                return
             else
                i_array1(num) = i
                j_array1(num) = j
                num = num + 1
             endif
          endif
       enddo
    enddo

    ! --- if there are more than one points with maximum first negibor,
    ! --- find the point with the maximum second neighbor.

    do n = 1, num_max_first
       if(second_nbr(i_array1(n), j_array1(n)) .gt. max_second_nbr) then
          max_second_nbr = second_nbr(i_array1(n), j_array1(n)) 
          num_max_second = 1
          i_array2(1) = i_array1(n)
          j_array2(1) = j_array1(n)
       else if(second_nbr(i_array1(n), j_array1(n)) == max_second_nbr) then
          num_max_second = num_max_second + 1
          i_array2(num_max_second) = i_array1(n)
          j_array2(num_max_second) = j_array1(n)
       endif
    enddo

    if(num_max_second == 1) then
       i_out = i_array2(1)
       j_out = j_array2(1)
       return
    endif

    ! --- then find the point with maximum land percentage (or smallest land/sea mask)
    max_mask = 0.0
    do n = 1, num_max_second
       if(lnd_mask(i_array2(n), j_array2(n)) .gt. max_mask) then
          max_mask = lnd_mask(i_array2(n), j_array2(n))
          num_max_mask = 1
          i_out = i_array1(n)
          j_out = j_array1(n)
       else if(lnd_mask(i_array2(n), j_array2(n)) == max_mask ) then
          max_mask = max_mask + 1
       endif
    enddo

    if(num_max_mask .gt. 1) call mpp_error(FATAL,'river_regrid: The points with maximum '// &
         'first and second neighbor land and percentage of land part is not unique ')

  end subroutine search_max_land

  !#####################################################################
  !--- find the river point to be defined.
  subroutine search_max_river(first_nbr, second_nbr, i_river, j_river, i_out, j_out, num_out)
    real,    dimension(:), intent(in)  :: first_nbr, second_nbr
    integer, dimension(:), intent(in)  :: i_river, j_river
    integer, dimension(:), intent(out) :: i_out, j_out 
    integer,               intent(out) :: num_out
    integer, dimension(size(i_out))    :: i_array, j_array  
    integer                            :: n, m
    integer                            :: num_max_first, num_max_second
    real                               :: max_first_nbr, max_second_nbr    

    !--- first fine the points with maximum first nbr
    max_first_nbr = -1.0
    max_second_nbr = -1
    num_max_first = 0

    do n = 1, size(first_nbr)
       if(first_nbr(n) > max_first_nbr + epsln) then
          num_max_first = 1
          i_out(1)   = i_river(n)
          j_out(1)   = j_river(n)
          max_first_nbr = first_nbr(n)
       else if(abs(first_nbr(n)-max_first_nbr) .le. epsln) then
          num_max_first = num_max_first + 1
          i_out(num_max_first) = i_river(n)
          j_out(num_max_first) = j_river(n)
       endif
    enddo

    if(num_max_first == 1)  then
       num_out = 1
       return
    endif

    i_array = i_out
    j_array = j_out

    ! --- if there are more than one points with maximum first negibor,
    ! --- find the point with the maximum second neighbor.

    do n = 1, num_max_first
       do m = 1, size(first_nbr)
          if(i_river(m) == i_array(n) .and. j_river(m) == j_array(n)) then
             if(second_nbr(m) .gt. max_second_nbr+ epsln) then
                max_second_nbr = second_nbr(m) 
                num_max_second = 1
                i_out(1) = i_array(n)
                j_out(1) = j_array(n)
             else if(abs(second_nbr(m)-max_second_nbr) .le. epsln) then
                num_max_second = num_max_second + 1
                i_out(num_max_second) = i_array(n)
                j_out(num_max_second) = j_array(n)
             endif
          endif
       enddo
    enddo

    num_out = num_max_second

  end subroutine search_max_river

  !#####################################################################
  !--- update point on the global halo.
  subroutine update_halo(data, halo)
    integer,                       intent(in) :: halo
    real, dimension(1-halo:,:), intent(inout) :: data

    integer                             :: n, ni

    ni = size(data,1) - 2*halo
    do n = 1, halo

       data(1-n,:)      = data(ni+1-n,:)
       data(ni+n,:) = data(halo,:)

    enddo

  end subroutine update_halo

  !#####################################################################
  !--- read each line of the file grid_edits
  subroutine parse_edits(txt,is,ie,js,je,flag)

    character(len=*), intent(in) :: txt
    integer, intent(inout)       :: is,ie,js,je
    logical, intent(inout)       :: flag

    integer :: i1,i2, i3
    character(len=128) :: txt2

    flag = .true.
    i1 = scan(txt,',')
    if (i1 <= 0) goto 90

    txt2 = txt(1:i1-1)

    i2 = scan(txt2,':')
    if (i2 <= 0) then
       read(txt2,*,err=90) is
       ie = is
    else
       read(txt2(1:i2-1),*,err=90) is
       read(txt2(i2+1:),*,err=90) ie
    endif

    txt2 = txt(i1+1:)
 
    i3 = scan(txt2,':')
     if (i3 <= 0) then
       read(txt2,*,err=90) js
       je = js
    else
       read(txt2(1:i3-1),*,err=90) js
       read(txt2(i3+1:),*,err=90) je
    endif

    return

90  continue

    flag = .false.

    return
  end subroutine parse_edits

  !#####################################################################


end program river_regrid
