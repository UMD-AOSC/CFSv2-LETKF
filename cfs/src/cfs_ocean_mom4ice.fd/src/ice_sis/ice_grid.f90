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
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_grid_mod - sets up grid, processor domain, and does advection-Michael.Winton@noaa.gov!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
module ice_grid_mod

  use constants_mod,   only: radius, omega, pi
  use mpp_mod,         only: mpp_pe, mpp_npes, mpp_root_pe, mpp_error, NOTE
  use mpp_mod,         only: mpp_sync_self, mpp_send, mpp_recv
  use mpp_domains_mod, only: mpp_define_domains, CYCLIC_GLOBAL_DOMAIN, FOLD_NORTH_EDGE
  use mpp_domains_mod, only: mpp_update_domains, domain2D, mpp_global_field
  use mpp_domains_mod, only: mpp_get_compute_domain, mpp_get_data_domain
  use mpp_io_mod,      only: mpp_open, mpp_close, mpp_get_info, mpp_get_atts
  use mpp_io_mod,      only: atttype, MPP_RDONLY, MPP_NETCDF 
  use fms_mod,         only: error_mesg, FATAL, field_exist, field_size, read_data

  implicit none
  include 'netcdf.inc'
  private

  public :: set_ice_grid, dt_evp, evp_sub_steps, g_sum, ice_avg, all_avg
  public :: Domain, isc, iec, jsc, jec, isd, ied, jsd, jed, im, jm, km
  public :: dtw, dte, dts, dtn, dxt, dxv, dyt, dyv, cor, xb1d, yb1d
  public :: geo_lon, geo_lat, sin_rot, cos_rot, cell_area, wett, wetv
  public :: dTdx, dTdy, t_on_uv, t_to_uv, uv_to_t, ice_advect
  public :: ice_line, vel_t_to_uv, cut_check, latitude, slab_ice_advect
  public :: dxdy, dydx, ice_grid_end
  public :: tripolar_grid

  type(domain2D), save :: Domain

  integer                              :: isc, iec, jsc, jec ! compute domain
  integer                              :: isd, ied, jsd, jed ! data domain
  integer                              :: im, jm, km         ! global domain and vertical size
  logical, allocatable, dimension(:,:) :: wett, wetv         ! t and v cell masks
  !
  ! grid geometry
  !
  logical                           ::  x_cyclic           ! x boundary condition
  logical                           ::  tripolar_grid      ! y boundary condition
  real, allocatable, dimension(:,:) ::  dtw, dte, dts, dtn ! size of t cell sides
  real, allocatable, dimension(:,:) ::  dxt, dxv           ! x-extent of t and v cells
  real, allocatable, dimension(:,:) ::  dyt, dyv           ! y-extent of t and v cells
  real, allocatable, dimension(:,:) ::  dxdy, dydx         
  real, allocatable, dimension(:,:) ::  latitude           ! latitude of t cells
  real, allocatable, dimension(:,:) ::  cor                ! coriolis on v cells
  real, allocatable, dimension(:,:) ::  geo_lat            ! true latitude (rotated grid)
  real, allocatable, dimension(:,:) ::  geo_lon            ! true longitude              
  real, allocatable, dimension(:  ) ::  xb1d, yb1d         ! 1d global grid for diag_mgr
  real, allocatable, dimension(:,:) ::  sin_rot, cos_rot   ! sin/cos of vector rotation angle
  real, allocatable, dimension(:,:) ::  cell_area          ! grid cell area; sphere frac.
  !
  ! timestep parameters
  !
  integer            :: evp_sub_steps                ! evp subcycles / timestep
  real               :: dt_evp = 0.0                 ! evp timestep (sec)
  integer            :: adv_sub_steps                ! advection steps / timestep
  real               :: dt_adv = 0.0                 ! advection timestep (sec)
  integer            :: comm_pe                      ! pe to be communicated with

contains
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! ice_avg - take area weighted average over ice partiions                      !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  function ice_avg(x,part)
    real, dimension(:,:,:),    intent(in) :: x
    real, dimension(:,:,:),    intent(in) :: part
    real, dimension(size(x,1),size(x,2) ) :: ice_avg, ice_part

    integer :: i, j, k

    ice_part = 0.0
    ice_avg  = 0.0
    if (size(x,3)==km) then
       do k=2,km
          do j = 1, size(x,2)
             do i = 1, size(x,1)
                ice_avg(i,j) = ice_avg(i,j) + part(i,j,k)*x(i,j,k)
             enddo
          enddo
       enddo
    else if (size(x,3)==km-1) then
       do k=2,km
          do j = 1, size(x,2)
             do i = 1, size(x,1)
                ice_avg(i,j) = ice_avg(i,j) + part(i,j,k)*x(i,j,k-1)
             enddo
          enddo
       enddo
    end if

    do k=2,km
       ice_part = ice_part + part(:,:,k)
    enddo

    do j = 1, size(x,2)
       do i = 1, size(x,1)
          if(ice_part(i,j) > 0 ) then
             ice_avg(i,j) = ice_avg(i,j)/ice_part(i,j)
          else
             ice_avg(i,j) = 0
          endif
       enddo
    enddo
    return
  end function ice_avg

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! all_avg - take area weighted average over all partiions                      !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  function all_avg(x,part)
    real, dimension(:,:, :),             intent(in) :: x
    real, dimension(isc:,jsc:,:), intent(in) :: part
    real, dimension(isc:iec,jsc:jec)         :: all_avg

    integer :: k

    all_avg = 0.0
    if (size(x,3)==km) then
       do k=1,km
          all_avg = all_avg + part(:,:,k)*x(:,:,k)
       end do
    else
       do k=2,km
          all_avg = all_avg + part(:,:,k)*x(:,:,k-1)
       end do
    end if
    return
  end function all_avg

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! g_sum - returns the global sum of a real array                               !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  real function g_sum(x)
    real, dimension(:,:) :: x

    real, dimension(1 :im, 1 :jm) :: g_x

    call mpp_global_field(Domain, x, g_x)
    g_sum = sum(g_x)

    return
  end function g_sum

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! uv_to_t - average v-values to t-points and apply mask                        !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine uv_to_t(uv, t)
    real,    intent(in   ), dimension(isd:,jsd:) :: uv
    real,    intent(  out), dimension(isc:,jsc:) :: t

    integer :: i, j

    do j = jsc, jec
       do i = isc, iec
          if(wett(i,j) ) then
             t(i,j) = 0.25*(uv(i,j) + uv(i,j-1) + uv(i-1,j) + uv(i-1,j-1) )   
          else
             t(i,j) = 0.0
          endif
       enddo
    enddo

  end subroutine uv_to_t

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! dTdx - eastward difference of tracer, result on uv cells                     !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  function dTdx(T)
    real, intent(in   ), dimension(isd:,jsd:) :: T 
    real, dimension(isc:iec,jsc:jec)          :: dTdx

    integer :: i, j

    do j = jsc, jec
       do i = isc, iec
          DTdx(i,j) = 0.5*(T(i+1,j+1) - T(i,j+1) + T(i+1,j) - T(i,j) )
       enddo
    enddo

    return
  end function dTdx

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! dTdy - northward difference of tracer, result on uv cells                    !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  function dTdy(T)
    real, intent(in   ), dimension(isd:,jsd:) :: T 
    real, dimension(isc:iec,jsc:jec)          :: dTdy

    integer :: i, j

    do j = jsc, jec
       do i = isc, iec
          DTdy(i,j) = 0.5*(T(i+1,j+1) - T(i+1,j) + T(i,j+1) - T(i,j) )
       enddo
    enddo

    return
  end function dTdy

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! t_on_uv - average tracer to uv cells                                      !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  function t_on_uv(T)
    real, intent(in   ), dimension(isd:,jsd:) :: T 
    real, dimension(isc:iec,jsc:jec)          :: t_on_uv

    integer :: i, j

    do j = jsc, jec
       do i = isc, iec
          t_on_uv(i,j) = 0.25*(T(i+1,j+1)+T(i+1,j)+T(i,j+1)+T(i,j) )
       enddo
    enddo

    return
  end function t_on_uv

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! t_to_uv - average t-values to v-points and apply mask                        !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine t_to_uv(t, uv)
    real,    intent(in ), dimension(isd:,jsd:) :: t
    real,    intent(out), dimension(isc:,jsc:) :: uv

    integer :: i, j

    do j = jsc, jec
       do i = isc, iec
          if(wetv(i,j) ) then
             uv(i,j) = 0.25*( t(i+1, j+1)+t(i+1, j) + t(i,j+1)+t(i,j) )
          else
             uv(i,j) = 0.0
          endif
       enddo
    enddo

  end subroutine t_to_uv

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! vel_t_to_uv - average vel component on t-points to v-points and apply mask   !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine vel_t_to_uv(tx, ty, uvx, uvy)
    real,    intent(in   ), dimension(isd:,jsd:) :: tx, ty
    real,    intent(  out), dimension(isc:,jsc:) :: uvx, uvy

    integer :: i, j

    do j = jsc, jec
       do i = isc, iec
          if( wetv(i,j) ) then
             uvx(i,j) = 0.25*( tx(i+1, j+1)+tx(i+1, j) + tx(i,j+1)+tx(i,j) )
             uvy(i,j) = 0.25*( ty(i+1, j+1)+ty(i+1, j) + ty(i,j+1)+ty(i,j) )
          else
             uvx(i,j) = 0.0
             uvy(i,j) = 0.0
          endif
       enddo
    enddo

  end subroutine vel_t_to_uv

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! set_ice_grid - initialize sea ice grid for dynamics and transport            !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine set_ice_grid(ice_domain, dt_slow, dyn_sub_steps_in, &
                          adv_sub_steps_in, km_in, layout )
    type(domain2D),        intent(inout) :: ice_domain
    real,                  intent(in)    :: dt_slow
    integer,               intent(in)    :: dyn_sub_steps_in
    integer,               intent(in)    :: adv_sub_steps_in
    integer,               intent(in)    :: km_in
    integer, dimension(2), intent(inout) :: layout

    real                                :: angle, lon_scale
    integer                             :: i, j, unit
    integer                             :: dims(4), ndim, nvar, natt, ntime
    real, allocatable, dimension(:,:)   :: rmask
    real, allocatable, dimension(:,:,:) :: x_vert_t, y_vert_t
    real, allocatable,   dimension(:,:) :: geo_latv   
    real, allocatable,   dimension(:,:) :: geo_lonv   
    character(len=80)                   :: domainname
    character(len=128)                  :: grid_file
    type(atttype), allocatable, dimension(:) :: global_atts

    grid_file = 'INPUT/grid_spec.nc'
    call mpp_open(unit,trim(grid_file),MPP_RDONLY,MPP_NETCDF)
    call mpp_get_info(unit, ndim, nvar, natt, ntime)

    !--- get boundary condition
    x_cyclic = .true.; tripolar_grid = .true.
    allocate(global_atts(natt))
    call mpp_get_atts(unit,global_atts)
    do i=1,natt
       select case (trim(global_atts(i)%name))
       case ('x_boundary_type')
          if (trim(global_atts(i)%catt) == 'cyclic') then
             x_cyclic = .true.
             call mpp_error(NOTE,'==>Note from ice_grid_mod: x_boundary_type is cyclic')
          else if (trim(global_atts(i)%catt) == 'solid_walls') then
             x_cyclic = .false.
             call mpp_error(NOTE,'==>Note from ice_grid_mod: x_boundary_type is solid_walls')
          endif
       case ('y_boundary_type')
          if (trim(global_atts(i)%catt) == 'fold_north_edge') then
             tripolar_grid = .true.
             call mpp_error(NOTE,'==>Note from ice_grid_mod: y_boundary_type is fold_north_edge')
          else
             tripolar_grid = .false.
             call mpp_error(NOTE,'==>Note from ice_grid_mod: y_boundary_type is solid_walls')
          endif
       end select
    end do

    call mpp_close(unit)

    !--- get the grid size
    call field_size(grid_file, 'wet', dims)
    im = dims(1)
    jm = dims(2)

    ! default is merdional domain decomp. to load balance xgrid
    if( layout(1)==0 .and. layout(2)==0 ) layout=(/ mpp_npes(), 1 /)
    if( layout(1)/=0 .and. layout(2)==0 ) layout(2) = mpp_npes()/layout(1)
    if( layout(1)==0 .and. layout(2)/=0 ) layout(1) = mpp_npes()/layout(2)
    domainname = 'ice model'
    if(tripolar_grid) then    
       call mpp_define_domains( (/1,im,1,jm/), layout, Domain,        &
                                xflags=CYCLIC_GLOBAL_DOMAIN, xhalo=1, &
                                yflags=FOLD_NORTH_EDGE, yhalo=1, name=domainname )
    else if(x_cyclic) then
       call mpp_define_domains( (/1,im,1,jm/), layout, Domain,        &
                                xflags=CYCLIC_GLOBAL_DOMAIN, xhalo=1, yhalo=1, name=domainname )
    else
       call mpp_define_domains( (/1,im,1,jm/), layout, Domain,        &
                                xhalo=1, yhalo=1, name=domainname )
    endif

    call mpp_get_compute_domain( Domain, isc, iec, jsc, jec )
    call mpp_get_data_domain( Domain, isd, ied, jsd, jed )

    call mpp_define_domains( (/1,im,1,jm/), layout, ice_domain) ! domain without halo

    allocate ( geo_lonv(1:im+1, 1:jm+1), geo_latv(1:im+1, 1:jm+1) )

    allocate ( wett   (isd:ied,jsd:jed), wetv   ( isc:iec,jsc:jec),  &
               dxt    (isd:ied,jsd:jed), dxv     (isd:ied,jsd:jed),  &
               dyt    (isd:ied,jsd:jed), dyv     (isd:ied,jsd:jed),  &
               dtw    (isc:iec,jsc:jec), dte     (isd:ied,jsd:jed),  &
               dts    (isc:iec,jsc:jec), dtn     (isd:ied,jsd:jed),  &
               dxdy   (isc:iec,jsc:jec), dydx    (isc:iec,jsc:jec),  &
               geo_lat(isc:iec,jsc:jec), geo_lon (isc:iec,jsc:jec),  &
               cor    (isc:iec,jsc:jec), sin_rot (isd:ied,jsd:jed),  &
               cos_rot(isd:ied,jsd:jed), latitude(isc:iec,jsc:jec) )
    allocate ( rmask (isc:iec, jsc:jec) )
    !--- read data from grid_spec.nc
    call read_data(grid_file, 'wet', rmask, Domain)
    wett = .false.
    do j = jsc, jec
       do i = isc, iec
          if( rmask(i,j) > 0.5) wett(i,j) = .true.
       enddo
    enddo

    call mpp_update_domains(wett, Domain)

    do j = jsc, jec
       do i = isc, iec
          if( wett(i,j) .and. wett(i,j+1) .and. wett(i+1,j) .and. wett(i+1,j+1) ) then
             wetv(i,j) = .true.
          else
             wetv(i,j) = .false.
          endif
       enddo
    enddo

    deallocate ( rmask )

    if(tripolar_grid) then
       if (jsc==1.and.any(wett(:,jsc))) call error_mesg ('ice_model_mod', &
            'ice model requires southernmost row of land', FATAL);
    endif

    evp_sub_steps = dyn_sub_steps_in
    if (evp_sub_steps>0) then
       dt_evp = dt_slow/evp_sub_steps
    else                            ! evp_sub_steps==0 means no dynamics
       dt_evp = dt_slow              ! but set dt>0 to avoid divide by zero
    end if

    adv_sub_steps = adv_sub_steps_in
    if (adv_sub_steps>0) then
       dt_adv = dt_slow/adv_sub_steps
    else                            ! adv_sub_steps==0 means no advection
       dt_adv = dt_slow              ! but set dt>0 to avoid divide by zero
    end if

    km = km_in

    if(field_exist(grid_file, 'x_vert_T')) then
       allocate ( x_vert_t(im,jm,4), y_vert_t(im,jm,4) )
       call read_data(grid_file, 'x_vert_T', x_vert_t)
       call read_data(grid_file, 'y_vert_T', y_vert_t)
       geo_lonv(1:im,1:jm) = x_vert_t(1:im,1:jm,1)
       geo_lonv(im+1,1:jm) = x_vert_t(im,1:jm,2)
       geo_lonv(1:im,jm+1) = x_vert_t(1:im,jm,4)
       geo_lonv(im+1,jm+1) = x_vert_t(im,jm,3)
       geo_latv(1:im,1:jm) = y_vert_t(1:im,1:jm,1)
       geo_latv(im+1,1:jm) = y_vert_t(im,1:jm,2)
       geo_latv(1:im,jm+1) = y_vert_t(1:im,jm,4)
       geo_latv(im+1,jm+1) = y_vert_t(im,jm,3)
       deallocate(x_vert_t, y_vert_t)
    else if(field_exist(grid_file, 'geolon_vert_t')) then
       call read_data(grid_file, 'geolon_vert_t', geo_lonv)
       call read_data(grid_file, 'geolat_vert_t', geo_latv)      
    else
       call error_mesg('ice_grid_mod', 'both x_vert_T and geolon_vert_t is not in the grid file ' &
                            //trim(grid_file), FATAL )
    endif

    allocate ( cell_area(isc:iec,jsc:jec) )
    call read_data(grid_file, 'AREA_OCN', cell_area, Domain)

    allocate ( xb1d (im+1), yb1d (jm+1) )

    xb1d = sum(geo_lonv,2)/(jm+1);
    yb1d = sum(geo_latv,1)/(im+1);

    dte = 0.0
    dtn = 0.0
    cos_rot = 0.0
    sin_rot = 0.0

    do j=jsc,jec
       do i=isc,iec
          dts(i,j)     = edge_length(geo_lonv(i,j),  geo_latv(i,j), geo_lonv(i+1,j),geo_latv(i+1,j))
          dtn(i,j)     = edge_length(geo_lonv(i,j+1),geo_latv(i,j+1), geo_lonv(i+1,j+1),geo_latv(i+1,j+1))
          dtw(i,j)     = edge_length(geo_lonv(i,j),  geo_latv(i,j), geo_lonv(i,j+1),geo_latv(i,j+1))
          dte(i,j)     = edge_length(geo_lonv(i+1,j),geo_latv(i+1,j),geo_lonv(i+1,j+1),geo_latv(i+1,j+1))
          lon_scale    = cos((geo_latv(i,j  )+geo_latv(i+1,j  )+geo_latv(i,j+1)+geo_latv(i+1,j+1))*atan(1.0)/180)
          angle        = atan2((geo_lonv(i,j+1)+geo_lonv(i+1,j+1)-geo_lonv(i,j)-geo_lonv(i+1,j))&
                         *lon_scale, geo_latv(i,j+1)+geo_latv(i+1,j+1)-geo_latv(i,j)-geo_latv(i+1,j) )
          sin_rot(i,j) = sin(angle) ! angle is the clockwise angle from lat/lon to ocean
          cos_rot(i,j) = cos(angle) ! grid (e.g. angle of ocean "north" from true north)
       end do
    end do

    call mpp_update_domains(dte, Domain)
    call mpp_update_domains(dtn, Domain)
    call mpp_update_domains(cos_rot, Domain)
    call mpp_update_domains(sin_rot, Domain)

    dxt = 0.0
    dyt = 0.0

    do j = jsc, jec
       do i = isc, iec
          dxt(i,j) = (dts(i,j) + dtn(i,j) )/2
          if(cell_area(i,j) > 0.0) then
             dyt(i,j) = cell_area(i,j)*4*pi*radius*radius/dxt(i,j)
          else
             dyt(i,j) = (dtw(i,j) + dte(i,j) )/2
          endif
       enddo
    enddo

    call mpp_update_domains(dxt, Domain )
    call mpp_update_domains(dyt, Domain )

    dxv = 1.0
    dyv = 1.0

    dxv(isc:iec,jsc:jec) = t_on_uv(dxt)
    dyv(isc:iec,jsc:jec) = t_on_uv(dyt)

    call mpp_update_domains(dxv, Domain )
    call mpp_update_domains(dyv, Domain )

    !--- dxdy and dydx to be used by ice_dyn_mod.
    dydx = dTdx(dyt)
    dxdy = dTdy(dxt)

    do j=jsc,jec
       do i = isc,iec
          geo_lon(i,j) = lon_avg( (/ geo_lonv(i,j  ), geo_lonv(i+1,j  ), &
                         geo_lonv(i,j+1), geo_lonv(i+1,j+1) /) )
       end do
    end do

    geo_lat  = (geo_latv(isc:iec,jsc:jec) + geo_latv(isc+1:iec+1,jsc:jec) &
               +geo_latv(isc:iec,jsc+1:jec+1)+geo_latv(isc+1:iec+1,jsc+1:jec+1))/4
    cor      = 2*omega*sin(geo_latv(isc+1:iec+1,jsc+1:jec+1)*pi/180)
    latitude = geo_lat(isc:iec,jsc:jec)*pi/180

    deallocate (geo_lonv, geo_latv)

!    comm_pe = mpp_npes() - mpp_pe() + 2*mpp_root_pe() - 1
!   consider arbitray layout
    comm_pe = mpp_pe() + layout(1) - 2*mod(mpp_pe()-mpp_root_pe(),layout(1)) - 1

    return
  end subroutine set_ice_grid

  !#####################################################################
  !--- release memory
  subroutine ice_grid_end

     deallocate(wett, wetv, dtw, dte, dts, dtn, dxt, dxv, dyt, dyv, latitude )
     deallocate(cor, geo_lat, geo_lon, xb1d, yb1d, sin_rot, cos_rot, cell_area )


  end subroutine ice_grid_end

  !#####################################################################

  real function lon_avg(lons)
    real, dimension(:), intent(in) :: lons

    real, dimension(size(lons(:))) :: lons2 ! lons relative to lon(1)
    integer                        :: i

    lons2(1) = 0.0
    do i=2,size(lons(:))
       lons2(i) = lons(i)-lons(1)
       if (lons2(i) >  180) lons2(i) = lons2(i) - 360;
       if (lons2(i) < -180) lons2(i) = lons2(i) + 360;
    end do
    lon_avg = lons(1)+sum(lons2)/size(lons(:))
  end function lon_avg

  !#####################################################################
  function edge_length(x1, y1, x2, y2)
    real, intent(in) :: x1, x2, y1, y2 ! end-point coordinates in degrees
    real             :: edge_length
    real             :: dx, dy

    dx = (x2-x1)*cos((atan(1.0)/45)*(y2+y1)/2)
    dy = y2-y1
    edge_length = radius*(atan(1.0)/45)*(dx*dx+dy*dy)**0.5
  end function edge_length

  !#####################################################################
  subroutine ice_line(year, day, second, cn, sst)
    integer,                               intent(in) :: year, day, second
    real, dimension(isc:,jsc:,1:),   intent(in) :: cn
    real, dimension(isc:,jsc:),      intent(in) :: sst

    real, dimension(isc:iec,jsc:jec) :: x
    real                             :: gx(3)
    integer                          :: i, k

    do i=-1,1,2
       x = 0.0
       where (cn(:,:,1)<0.85 .and. i*geo_lat>0.0) x=cell_area
       gx((i+3)/2) = g_sum(x)*4*pi*radius*radius/1e12
    end do
    gx(3) = g_sum(sst*cell_area)/(g_sum(cell_area)+1e-10)
    !
    ! print info every 5 days
    !
    if ( mpp_pe()==0 .and. second==0 .and. mod(day,5)==0 ) &
       print '(a,2I4,3F10.5)','ICE y/d (SH_ext NH_ext SST):', year, day, gx
  end subroutine ice_line

  !#####################################################################
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! ice_advect - take adv_sub_steps upstream advection timesteps                 !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine ice_advect(ui, vi, trc, uf, vf)
    real, intent(in   ),           dimension(isd:,jsd:) :: ui, vi  ! advecting velocity
    real, intent(inout),           dimension(isd:,jsd:) :: trc     ! tracer to advect
    real, optional, intent(inout), dimension(isc:,jsc:) :: uf, vf

    integer                          :: l, i, j
    real, dimension(isd:ied,jsd:jed) :: ue, vn
    real, dimension(isd:ied,jsd:jed) ::  uflx, vflx

    if (adv_sub_steps==0) return;

    ue(:,jsd) = 0.0; ue(:,jed) = 0.0
    vn(:,jsd) = 0.0; vn(:,jed) = 0.0

    do j = jsc, jec
       do i = isc, iec
          ue(i,j) = 0.5 * ( ui(i,j-1) + ui(i,j) )
          vn(i,j) = 0.5 * ( vi(i-1,j) + vi(i,j) )
       enddo
    enddo

    call mpp_update_domains(ue, Domain)
    call mpp_update_domains(vn, Domain)

    if (present(uf)) uf = 0.0
    if (present(vf)) vf = 0.0

    uflx = 0.0
    vflx = 0.0

    do l=1,adv_sub_steps
       do j = jsd, jec
          do i = isd, iec
             if( ue(i,j) > 0.0 ) then
                uflx(i,j) = ue(i,j) * trc(i,j) * dte(i,j)
             else
                uflx(i,j) = ue(i,j) * trc(i+1,j) * dte(i,j)
             endif

             if( vn(i,j) > 0.0 ) then
                vflx(i,j) = vn(i,j) * trc(i,j) * dtn(i,j)
             else
                vflx(i,j) = vn(i,j) * trc(i,j+1) * dtn(i,j)
             endif
          enddo
       enddo

       do j = jsc, jec
          do i = isc, iec
             trc(i,j) = trc(i,j) + dt_adv * ( uflx(i-1,j) - uflx(i,j) + &
                        vflx(i,j-1) - vflx(i,j) )/ ( dxt(i,j) * dyt(i,j) ) 
          enddo
       enddo

       call mpp_update_domains(trc, Domain)

       if (present(uf)) then
          do j = jsc, jec
             do i = isc, iec       
                uf(i,j) = uf(i,j) + uflx(i,j)
             enddo
          enddo
       endif

       if (present(vf)) then
          do j = jsc, jec
             do i = isc, iec       
                vf(i,j) = vf(i,j) + vflx(i,j)
             enddo
          enddo
       endif

    end do

    if (present(uf)) uf = uf/adv_sub_steps;
    if (present(vf)) vf = vf/adv_sub_steps;

  end subroutine ice_advect

  !#####################################################################
  subroutine slab_ice_advect(ui, vi, trc, stop_lim)
    real, intent(in   ), dimension(isd:,jsd:) :: ui, vi       ! advecting velocity
    real, intent(inout), dimension(isd:,jsd:) :: trc          ! tracer to advect
    real, intent(in   )                             :: stop_lim

    integer                          :: l, i, j
    real, dimension(isd:ied,jsd:jed) :: ue, vn, uflx, vflx
    real                             :: avg, dif

    if (adv_sub_steps==0) return;

    ue(:,jsd) = 0.0; ue(:,jed) = 0.0
    vn(:,jsd) = 0.0; vn(:,jed) = 0.0

    do j = jsc, jec
       do i = isc, iec
          ue(i,j) = 0.5 * ( ui(i,j-1) + ui(i,j) )
          vn(i,j) = 0.5 * ( vi(i-1,j) + vi(i,j) )
       enddo
    enddo

    call mpp_update_domains(ue, Domain)
    call mpp_update_domains(vn, Domain)

    do l=1,adv_sub_steps
       do j = jsd, jec
          do i = isd, iec
             avg = ( trc(i,j) + trc(i+1,j) )/2
             dif = trc(i+1,j) - trc(i,j)
             if( avg > stop_lim .and. ue(i,j) * dif > 0.0) then
                uflx(i,j) = 0.0
             else if( ue(i,j) > 0.0 ) then
                uflx(i,j) = ue(i,j) * trc(i,j) * dte(i,j)
             else
                uflx(i,j) = ue(i,j) * trc(i+1,j) * dte(i,j)
             endif

             avg = ( trc(i,j) + trc(i,j+1) )/2
             dif = trc(i,j+1) - trc(i,j)
             if( avg > stop_lim .and. vn(i,j) * dif > 0.0) then
                vflx(i,j) = 0.0
             else if( vn(i,j) > 0.0 ) then
                vflx(i,j) = vn(i,j) * trc(i,j) * dtn(i,j)
             else
                vflx(i,j) = vn(i,j) * trc(i,j+1) * dtn(i,j)
             endif
          enddo
       enddo

       do j = jsc, jec
          do i = isc, iec
             trc(i,j) = trc(i,j) + dt_adv * ( uflx(i-1,j)-uflx(i,j) + vflx(i,j-1)-vflx(i,j) ) / (dxt(i,j)*dyt(i,j))
          enddo
       enddo

       call mpp_update_domains(trc, Domain)

    end do
  end subroutine slab_ice_advect

  !#####################################################################
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! cut_check - check for mirror symmetry of uv field along bipolar grid cut     !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine cut_check(mesg, uv)
    character(len=*), intent(in) :: mesg
    real, intent(inout), dimension(isd:,jsd:) :: uv

!    real, dimension(1:im, 1:jm)  :: global_uv

    real    :: fix
    integer :: i, j, l
    logical :: cut_error = .false.
    real, dimension(iec-isc+1) :: send_buffer, recv_buffer

!    call mpp_sync_self()

    l = 0
    do i = iec-1, isd, -1
       l = l+1
       send_buffer(l) = uv(i,jec)
    enddo

    call mpp_send(send_buffer(1), plen=(iec-isc+1), to_pe=comm_pe)
    call mpp_recv(recv_buffer(1), glen=(iec-isc+1), from_pe=comm_pe)

    ! cut_check only at the north pole.
    if( jec == jm) then
       l = 0
       do i = isc, iec
          l=l+1
          ! excludes the point at im/2, and im.
          if(i == im/2 .or. i == im) cycle
          uv(i,jec) = (uv(i,jec) - recv_buffer(l))/2
       enddo
    endif

    call mpp_sync_self()

    !    call mpp_global_field ( Domain, uv, global_uv )

    !    do i=1,im/2-1
    !       if (global_uv(i,jm)/=-global_uv(im-i,jm)) then
    !      cut_error =.true.
    !     if (mpp_pe()==0) &
         !       print *, mesg, i, im-i, global_uv(i,jm), global_uv(im-i,jm)
    !          fix = (global_uv(i,jm)-global_uv(im-i,jm))/2
    !          global_uv(i   ,jm) =  fix
    !          global_uv(im-i,jm) = -fix
    !       end if
    !    end do

    !    uv(isc:iec,jsc:jec) = global_uv(isc:iec,jsc:jec)

    ! if (cut_error) call error_mesg ('ice_model_mod', &
         !                     'cut check of mirror anti-symmetry failed', FATAL)
  end subroutine cut_check
  !#####################################################################

end module ice_grid_mod
