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
module land_model_mod

  use time_manager_mod, only: time_type

  use mpp_domains_mod,  only: domain1D, domain2d, mpp_get_layout, mpp_define_layout, &
                              mpp_define_domains, mpp_get_compute_domain,            &
                              CYCLIC_GLOBAL_DOMAIN, mpp_get_compute_domains,         &
                              mpp_get_domain_components, mpp_get_pelist

  use fms_mod,          only: write_version_number, error_mesg, FATAL, NOTE, &
                              mpp_pe, mpp_npes, mpp_root_pe

  use diag_manager_mod, only : send_data
  
implicit none

private

type land_data_type
   type(domain2d) :: domain  ! our computation domain
   real, pointer, dimension(:,:,:)   :: &  ! (lon, lat, tile)
        tile_size =>NULL(),       & ! fractional coverage of cell by tile, dimensionless
        t_surf =>NULL(),          & ! ground surface temperature, degK
        t_ca =>NULL(),            & ! canopy air temperature, degK
        q_ca =>NULL(),            & ! canopy air specific humidity, kg/kg
        albedo =>NULL(),          & ! snow-adjusted land albedo
        albedo_vis_dir =>NULL(),  & ! albedo for direct visible-band radiation
        albedo_nir_dir =>NULL(),  & ! albedo for direct nir-band radiation 
        albedo_vis_dif =>NULL(),  & ! albedo for diffuse visible-band radiation 
        albedo_nir_dif =>NULL(),  & ! albedo for diffuse nir-band radiation
        rough_mom =>NULL(),       & ! momentum roughness length, m
        rough_heat =>NULL(),      & ! roughness length for tracers (heat and water), m
        rough_scale =>NULL()        ! roughness length for drag scaling

! --> for cpl
   real, pointer, dimension(:,:,:)   :: &  ! (lon, lat, tile)
        t_sst =>null(),           &
        cice  =>null(),           &
        fice  =>null(),           &
        hice  =>null(),           &
        hsno  =>null()
! --> end for cpl

   real, pointer, dimension(:,:) :: &  ! (lon, lat)
        discharge =>NULL(),       & ! flux from surface drainage network out of land model
        discharge_snow =>NULL()     ! snow analogue of discharge

   logical, pointer, dimension(:,:,:):: &
        mask =>NULL()                ! true if land

   integer :: axes(2)      ! axes IDs for diagnostics  

end type land_data_type

type atmos_land_boundary_type
   ! data passed from the atmosphere to the surface
   real, dimension(:,:,:), pointer :: &
        t_flux =>NULL(),  &
        q_flux =>NULL(),  &
        lw_flux =>NULL(), &
        lwdn_flux =>NULL(), &
        sw_flux =>NULL(), &
        sw_flux_down_vis_dir =>NULL(), &
        sw_flux_down_total_dir =>NULL(), &
        sw_flux_down_vis_dif =>NULL(), &
        sw_flux_down_total_dif =>NULL(), &
        lprec =>NULL(),   &
        fprec =>NULL()
   ! derivatives of the fluxes
   real, dimension(:,:,:), pointer :: &
        dhdt =>NULL(),    &
        dedt =>NULL(),    &
        dedq =>NULL(),    &
        drdt =>NULL()
   real, dimension(:,:,:), pointer :: &
        cd_m => NULL(),      &   ! drag coefficient for momentum, dimensionless
        cd_t => NULL(),      &   ! drag coefficient for tracers, dimensionless
        ustar => NULL(),     &   ! turbulent wind scale, m/s
        bstar => NULL(),     &   ! turbulent bouyancy scale, m/s
        wind => NULL(),      &   ! abs wind speed at the bottom of the atmos, m/s
        z_bot => NULL(),     &   ! height of the bottom atmos. layer above surface, m
        drag_q =>NULL(),  &
        p_surf =>NULL()

   real, dimension(:,:,:), pointer :: &
        data =>NULL() ! collective field for "named" fields above
   integer :: xtype             !REGRID, REDIST or DIRECT
end type atmos_land_boundary_type

! ==== public interfaces =====================================================
public land_model_init          ! initialize the land model
public land_model_end           ! finish land model calculations

public update_land_model_fast   ! fast time-scale update of the land
public update_land_model_slow   ! slow time-scale update of the land
public send_averaged_data       ! send tile-averaged data to diag manager
! ==== end of public interfaces ==============================================

! ==== public data type =====================================================
public land_data_type, atmos_land_boundary_type


! ==== some names, for information only ======================================
logical                       :: module_is_initialized = .FALSE.
character(len=*),   parameter :: module_name = 'land_mod'
character(len=128), parameter :: version     = ''
character(len=128), parameter :: tagname     = ''



! ==== local module variables ================================================
integer            :: n_tiles = 1  ! number of tiles

! ---- interface to tile-averaged diagnostic routines ------------------------
interface send_averaged_data
   module procedure send_averaged_data2d
   module procedure send_averaged_data3d
end interface

contains

! ============================================================================
subroutine land_model_init &
     ( atmos2land, Land_bnd, time_init, time, dt_fast, dt_slow, atmos_domain)

  use fms_mod, only : field_size, read_data
  use mpp_domains_mod, only : mpp_define_layout, mpp_define_domains,&
                              CYCLIC_GLOBAL_DOMAIN, mpp_get_data_domain
  use diag_manager_mod, only : diag_axis_init
! ============================================================================
! initialiaze land model using grid description file as an input. This routine
! reads land grid boundaries and area of land from a grid description file

! NOTES: theoretically, the grid description file can specify any regular
! rectangular grid for land, not jast lon/lat grid. Therefore the variables
! "xbl" and "ybl" in NetCDF grid spec file are not necessarily lon and lat
! boundaries of the grid.
!   However, at this time the module land_properties assumes that grid _is_
! lon/lat and therefore the entire module also have to assume that the land 
! grid is lon/lat.
!   lon/lat grid is also assumed for the diagnostics, but this is probably not
! so critical. 

  ! ---- arguments -----------------------------------------------------------
  type(atmos_land_boundary_type), intent(inout) :: atmos2land ! land boundary data
  type(land_data_type), intent(inout) :: Land_bnd ! land boundary data
  type(time_type), intent(in)    :: time_init ! initial time of simulation (?)
  type(time_type), intent(in)    :: time      ! current time
  type(time_type), intent(in)    :: dt_fast   ! fast time step
  type(time_type), intent(in)    :: dt_slow   ! slow time step
  type(domain2d),  intent(in), optional :: atmos_domain ! domain of computations

  real, dimension(:), allocatable :: glonb, glatb, glon, glat
  real, dimension(:,:), allocatable ::  gfrac, garea
  integer :: siz(4), layout(2), is, ie, js, je, i, j, k
  
  
  call field_size('INPUT/grid_spec','AREA_LND',siz)
  allocate(gfrac(siz(1),siz(2)))
  allocate(garea(siz(1),siz(2)))
  allocate(glonb(siz(1)+1))
  allocate(glatb(siz(2)+1))
  allocate(glon(siz(1)))
  allocate(glat(siz(2)))  
  call read_data('INPUT/grid_spec','xbl',glonb)
  call read_data('INPUT/grid_spec','ybl',glatb)
  call read_data('INPUT/grid_spec','AREA_LND',gfrac)
  call read_data('INPUT/grid_spec','AREA_LND_CELL',garea)

  glonb(:) = glonb(:)
  glatb(:) = glatb(:)

  glon(1:siz(1)) = (glonb(1:siz(1))+glonb(2:siz(1)+1))/2.0
  glat(1:siz(2)) = (glatb(1:siz(2))+glatb(2:siz(2)+1))/2.0

  gfrac = gfrac/garea
  
  call mpp_define_layout((/1,siz(1),1,siz(2)/), mpp_npes(), layout)
  call mpp_define_domains((/1,siz(1),1,siz(2)/), layout, Land_bnd%domain, &
     xflags = CYCLIC_GLOBAL_DOMAIN)

  call mpp_get_data_domain(Land_bnd%domain,is,ie,js,je)

  allocate ( &
       Land_bnd % tile_size      (is:ie,js:je,n_tiles), & 
       Land_bnd % t_surf         (is:ie,js:je,n_tiles), &
       Land_bnd % t_ca           (is:ie,js:je,n_tiles), &
       Land_bnd % q_ca           (is:ie,js:je,n_tiles), &
       Land_bnd % albedo         (is:ie,js:je,n_tiles), & 
       Land_bnd % albedo_vis_dir (is:ie,js:je,n_tiles), &
       Land_bnd % albedo_nir_dir (is:ie,js:je,n_tiles), &
       Land_bnd % albedo_vis_dif (is:ie,js:je,n_tiles), &
       Land_bnd % albedo_nir_dif (is:ie,js:je,n_tiles), &
       Land_bnd % rough_mom      (is:ie,js:je,n_tiles), & 
       Land_bnd % rough_heat     (is:ie,js:je,n_tiles), & 
       Land_bnd % rough_scale    (is:ie,js:je,n_tiles), & 
       Land_bnd % discharge      (is:ie,js:je),   & 
       Land_bnd % discharge_snow (is:ie,js:je),   &
       Land_bnd % mask       (is:ie,js:je,n_tiles)    )
!--> cpl insertion
  allocate(  Land_bnd % t_sst         (is:ie,js:je,n_tiles))
  allocate(  Land_bnd % cice (is:ie,js:je,n_tiles))
  allocate(  Land_bnd % fice (is:ie,js:je,n_tiles))
  allocate(  Land_bnd % hice (is:ie,js:je,n_tiles))
  allocate(  Land_bnd % hsno (is:ie,js:je,n_tiles))
!<-- cpl insertion
  
  do i=is,ie
     do j=js,je
        do k= 1, n_tiles
           if (gfrac(i,j) > 0.0) then
               Land_bnd%tile_size(i,j,k)  = 1.0/n_tiles
               Land_bnd%mask(i,j,k) = .true.
           else
               Land_bnd%tile_size(i,j,k) = 0.0
               Land_bnd%mask(i,j,k) = .false.
           endif
        enddo
     enddo
  enddo

!-->cpl insertion
  Land_bnd%t_sst = 273.0
  Land_bnd%cice  = 0.0
  Land_bnd%fice  = 0.0
  Land_bnd%hice  = 0.0
  Land_bnd%hsno  = 0.0
!<-- cpl insertion

  Land_bnd%t_surf = 273.0
  Land_bnd%t_ca = 273.0
  Land_bnd%q_ca = 0.0
  Land_bnd%albedo = 0.0
  Land_bnd % albedo_vis_dir = 0.0
  Land_bnd % albedo_nir_dir = 0.0
  Land_bnd % albedo_vis_dif = 0.0
  Land_bnd % albedo_nir_dif = 0.0
  Land_bnd%rough_mom = 0.01
  Land_bnd%rough_heat = 0.01
  Land_bnd%rough_scale = 1.0
  Land_bnd%discharge = 0.0
  Land_bnd%discharge_snow = 0.0
  Land_bnd%mask = .true.
  
  Land_bnd%axes(1) = diag_axis_init('lon',glon,'degrees_E','X','longitude',&
       set_name='land',domain2 = Land_bnd%domain)

  Land_bnd%axes(2) = diag_axis_init('lat',glon,'degrees_N','Y','latitude',&
       set_name='land',domain2 = Land_bnd%domain)  
  
  allocate( atmos2land % t_flux  (is:ie,js:je,n_tiles) )
  allocate( atmos2land % q_flux  (is:ie,js:je,n_tiles) )
  allocate( atmos2land % lw_flux (is:ie,js:je,n_tiles) )
  allocate( atmos2land % sw_flux (is:ie,js:je,n_tiles) )
  allocate( atmos2land % lprec   (is:ie,js:je,n_tiles) )
  allocate( atmos2land % fprec   (is:ie,js:je,n_tiles) )
  allocate( atmos2land % dhdt    (is:ie,js:je,n_tiles) )
  allocate( atmos2land % dedt    (is:ie,js:je,n_tiles) )
  allocate( atmos2land % dedq    (is:ie,js:je,n_tiles) )
  allocate( atmos2land % drdt    (is:ie,js:je,n_tiles) )
  allocate( atmos2land % drag_q  (is:ie,js:je,n_tiles) )
  allocate( atmos2land % p_surf  (is:ie,js:je,n_tiles) )
  allocate( atmos2land % sw_flux_down_vis_dir   (is:ie,js:je,n_tiles) )
  allocate( atmos2land % sw_flux_down_total_dir (is:ie,js:je,n_tiles) )
  allocate( atmos2land % sw_flux_down_vis_dif   (is:ie,js:je,n_tiles) )
  allocate( atmos2land % sw_flux_down_total_dif (is:ie,js:je,n_tiles) )

  return

end subroutine land_model_init






! ============================================================================
subroutine land_model_end ( atmos2land, bnd )
! ============================================================================
! destruct the land model data

  ! ---- arguments -----------------------------------------------------------
  type(atmos_land_boundary_type), intent(inout) :: atmos2land
  type(land_data_type), intent(inout) :: bnd

  module_is_initialized = .FALSE.

!  deallocate boundary exchange data
!  call deallocate_boundary_data ( bnd )
  
end subroutine land_model_end



! ============================================================================
subroutine update_land_model_fast ( atmos2land, bnd )
! ============================================================================
! updates state of the land model on the fast time scale

  ! ---- arguments -----------------------------------------------------------
  type(atmos_land_boundary_type), intent(inout)    :: atmos2land
  type(land_data_type),  intent(inout) :: bnd ! state to update

  return

end subroutine update_land_model_fast



! ============================================================================
subroutine update_land_model_slow ( atmos2land, bnd )
! ============================================================================
! updates land on slow time scale

  ! ---- arguments -----------------------------------------------------------
  type(atmos_land_boundary_type), intent(inout) :: atmos2land
  type(land_data_type), intent(inout) :: bnd


  return
  
end subroutine update_land_model_slow


! ============================================================================
subroutine update_land_bnd_fast ( bnd )
! ============================================================================
! updates land boundary data (the ones that atmosphere sees) on the fast time 
! scale this routine does not update tiling structure, because it is assumed
! that the tiling does not change on fast time scale

  ! ---- arguments -----------------------------------------------------------
  type(land_data_type), intent(inout) :: bnd

  return

end subroutine update_land_bnd_fast


! ============================================================================
subroutine update_land_bnd_slow ( bnd )
! ============================================================================
! updates land boundary data for the atmosphere on the slow time scale. This
! subroutine is responsible for the changing of tiling structure, if necessary,
! as well as for changing albedo, drag coefficients and such.

! NOTE: if tiling structure has been modified, then probably the distribution
! of other boundary values, such as temperature and surface humidity, should be
! modified too, not yet clear how.

  ! ---- arguments -----------------------------------------------------------
  type(land_data_type), intent(inout) :: bnd

  return

end subroutine update_land_bnd_slow

! ============================================================================
function send_averaged_data2d ( id, field, area, time, mask )
! ============================================================================
! average the data over tiles and send then to diagnostics

  ! --- return value ---------------------------------------------------------
  logical                      :: send_averaged_data2d
  ! --- arguments ------------------------------------------------------------
  integer, intent(in)          :: id             ! id od the diagnostic field 
  real,    intent(in)          :: field(:,:,:)   ! field to average and send
  real,    intent(in)          :: area (:,:,:)   ! area of tiles (== averaging 
                                                 ! weights), arbitrary units
  type(time_type), intent(in)  :: time           ! current time
  logical, intent(in),optional :: mask (:,:,:)   ! land mask

  ! --- local vars -----------------------------------------------------------
  real  :: out(size(field,1), size(field,2))

  call average_tiles( field, area, mask, out )
  send_averaged_data2d = send_data( id, out, time, mask=ANY(mask,DIM=3) )
end function send_averaged_data2d


! ============================================================================
function send_averaged_data3d( id, field, area, time, mask )
! ============================================================================
! average the data over tiles and send then to diagnostics

  ! --- return value ---------------------------------------------------------
  logical                      :: send_averaged_data3d
  ! --- arguments ------------------------------------------------------------
  integer, intent(in)          :: id              ! id of the diagnostic field
  real,    intent(in)          :: field(:,:,:,:)  ! (lon, lat, tile, lev) field 
                                                  ! to average and send
  real,    intent(in)          :: area (:,:,:)    ! (lon, lat, tile) tile areas 
                                                  ! ( == averaging weights), 
                                                  ! arbitrary units
  type(time_type), intent(in)  :: time            ! current time
  logical, intent(in),optional :: mask (:,:,:)    ! (lon, lat, tile) land mask

  ! --- local vars -----------------------------------------------------------
  real    :: out(size(field,1), size(field,2), size(field,4))
  logical :: mask3(size(field,1), size(field,2), size(field,4))
  integer :: it

  do it=1,size(field,4)
     call average_tiles( field(:,:,:,it), area, mask, out(:,:,it) )
  enddo

  mask3(:,:,1) = ANY(mask,DIM=3)
  do it = 2, size(field,4)
     mask3(:,:,it) = mask3(:,:,1)
  enddo

  send_averaged_data3d = send_data( id, out, time, mask=mask3 )
end function send_averaged_data3d

subroutine average_tiles ( x, area, mask, out )
! ============================================================================
! average 2-dimensional field over tiles
  ! --- arguments ------------------------------------------------------------
  real,    intent(in)  :: x   (:,:,:) ! (lon, lat, tile) field to average
  real,    intent(in)  :: area(:,:,:) ! (lon, lat, tile) fractional area
  logical, intent(in)  :: mask(:,:,:) ! (lon, lat, tile) land mask
  real,    intent(out) :: out (:,:)   ! (lon, lat)       result of averaging

  ! --- local vars -----------------------------------------------------------------
  integer  :: it                      ! iterator over tile number
  real     :: s(size(x,1),size(x,2))  ! area accumulator

  s(:,:)   = 0.0
  out(:,:) = 0.0

  do it = 1,size(area,3)
     where (mask(:,:,it)) 
        out(:,:) = out(:,:) + x(:,:,it)*area(:,:,it)
        s(:,:)   = s(:,:) + area(:,:,it)
     endwhere
  enddo

  where( s(:,:) > 0 ) &
       out(:,:) = out(:,:)/s(:,:)

end subroutine average_tiles

end module land_model_mod
   
