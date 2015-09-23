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
module ocean_obc_mod
  !
  ! <CONTACT EMAIL="mschmidt@io-warnemuende.de"> Martin Schmidt </CONTACT>
  ! <CONTACT EMAIL="Zhi.Liang@noaa.gov"> Zhi Liang </CONTACT>

  ! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>
  ! 
  ! <OVERVIEW>
  !   Open Boundary condition for mom4 ocean model.
  ! </OVERVIEW>
  !
  ! <DESCRIPTION>
  !   This module can extrapolate data on the open lateral boundaries for mom4 ocean model. Tracer and surface height 
  !   are extrapolated on the boundary by using implicit radiation boundary conditions, velocities are calculated
  !   on the boundary from a linear equation (omitted advection equation). The gradient of each field is supposed 
  !   to be zero accross boundary.
  ! </DESCRIPTION>

  use constants_mod,            only: grav, rho0
  use data_override_mod,        only: data_override
  use diag_manager_mod,         only: register_diag_field, send_data
  use field_manager_mod,        only: MODEL_OCEAN, fm_get_length
  use fms_mod,                  only: write_version_number, open_namelist_file, close_file
  use fms_mod,                  only: stdout, stdlog, check_nml_error, FATAL, NOTE 
  use mpp_domains_mod,          only: mpp_get_compute_domain, mpp_get_data_domain, mpp_get_global_domain
  use mpp_domains_mod,          only: domain2d, mpp_update_domains
  use mpp_mod,                  only: mpp_error, mpp_pe, mpp_chksum
  use mpp_mod,                  only: CLOCK_COMPONENT, CLOCK_SUBCOMPONENT, CLOCK_MODULE
  use mpp_mod,                  only: mpp_clock_id, mpp_clock_begin, mpp_clock_end
  use time_interp_external_mod, only: time_interp_external, init_external_field, get_external_field_size
  use time_manager_mod,         only: time_type
  use tracer_manager_mod,       only: get_tracer_names, get_tracer_indices, get_number_tracers

  use ocean_domains_mod,        only: get_local_indices, get_domain_offsets
  use ocean_types_mod,          only: ocean_grid_type, ocean_domain_type, ocean_thickness_type
  use ocean_types_mod,          only: ocean_prog_tracer_type, ocean_adv_vel_type
  use ocean_types_mod,          only: ocean_external_mode_type, obc_flux
  use ocean_types_mod,          only: ocean_time_type, ocean_time_steps_type

  implicit none
  private

  public :: ocean_obc_init, ocean_obc_end, ocean_obc_tracer, ocean_obc_freesurf
  public :: ocean_obc_update_boundary, ocean_obc_prepare
  public :: ocean_obc_tracer_flux, ocean_obc_mass_flux
  public :: store_ocean_obc_tracer_flux, store_ocean_obc_tracer_advect
  public :: ocean_obc_tracer_init, ocean_obc_adjust_advel
  public :: store_ocean_obc_pressure_grad, ocean_obc_adjust_forcing_fs

  !--- some module variables -------------------------------------------
  integer, parameter :: max_obc = 4             ! maximum number of open boundaries (increase if want more)
  integer, parameter :: max_prog_tracers=3      ! maximum number of prognostics tracers (increase if want more)

  integer, parameter :: WEST  = 1               ! western boundary
  integer, parameter :: EAST  = 2               ! eastern boundary
  integer, parameter :: SOUTH = 3               ! western boundary
  integer, parameter :: NORTH = 4               ! eastern boundary

                                                ! it can be increased if needed.
  real, parameter    :: small = 1.0e-8          ! some small number
  real               :: dtuv, dtts, dtfs, dteta ! time step.
  real               :: dtime_e                 ! timestep (secs) (dtime_e=2*dteta for threelevel and dtime_e=dteta for twolevel)
  real               :: dtime_t                 ! timestep (secs) (dtime_t=2*dtts for threelevel and dtime_t=dtts for twolevel)
  integer            :: isc, iec, jsc, jec      ! computational domain decompsition
  integer            :: isd, ied, jsd, jed      ! data domain decompsition
  integer            :: nk                      ! number of vertical levels

  type(ocean_domain_type), pointer :: Dom => NULL()  ! ocean domain
  type(ocean_grid_type), pointer   :: Grd => NULL()  ! ocean grid

  !--- derived types to specify the domain and grid --------------------
  type obc_bound_type
     integer           :: is, ie, js, je         ! index specifying bound location on current pe
     integer           :: isd,ied,jsd,jed        ! index specifying bound location on current pe of data domain
     integer           :: np                     ! specifies total number of boundary points
     integer           :: direction              ! its value can be "w"(west), "e"(east), "s"(south), "n"(north)     
     logical           :: on_bound               ! true means there are some points on the boundary on current pe.
     real,    pointer  :: ctrop(:)   => NULL()   ! barotropic phase speed
     real,    pointer  :: pgrad(:)   => NULL()   ! barotropic phase speed
     real,    pointer  :: hadv (:,:) => NULL()   ! cross boundary horizontal tracer advection, updated in the tracer loop, so no tracer number needed
     logical, pointer  :: mask(:)    => NULL()   ! logical varible to indicate the grid is open or not.
     character(len=16) :: name                   ! 
     real              :: ctrop_max              ! maximum barotropic phase speed (should be 1 .. 2) 
     real              :: ctrop_min              ! minimum barotropic phase speed (should be 0 .. 1) 
     logical           :: obc_consider_convu     ! true to consider convu for d eta/dt
     logical           :: obc_vert_advel_t       ! consider vertical advection for tracers 
     logical           :: obc_vert_advel_u       ! consider vertical advection for velocity 
     logical           :: obc_adjust_forcing_fs  ! remove baroclinic press gradients from forcing_fs
     logical           :: obc_relax_eta          ! true to relax sea level
     logical           :: obc_relax_eta_profile  ! true to relax sea level comparing mean sea level with reference value
     logical           :: obc_relax_tracer       ! true to relax at least one tracer
     logical           :: obc_tracer_orlanski    ! true to relax all tracers with Orlanski scheme
  end type obc_bound_type

  type obc_data_type_2d
     real, pointer     :: data(:,:) => NULL() ! boundary values of eta for relaxation
     character(len=128):: file, field         ! name of the inputfile and inputfield
     real              :: rel_coef            ! coefficient to data
     integer           :: id                  ! the data id for time interpolation
  end type obc_data_type_2d

  type obc_data_type_3d
     real, pointer     :: data(:,:,:) => NULL() ! boundary values of eta for relaxation
     character(len=128):: file, field         ! name of the inputfile and inputfield
     real              :: rel_coef            ! coefficient to data
     integer           :: id                  ! the data id for time interpolation
  end type obc_data_type_3d

  type ocean_obc_type
     integer                       :: nobc               ! number of boundaries
     type(obc_bound_type), pointer :: bound(:) => NULL() ! contains boundary information.
     type(obc_data_type_2d),  pointer :: eta(:)      => NULL() ! contains boundary data information for sea level.
     type(obc_data_type_3d),  pointer :: tracer(:,:) => NULL() ! contains boundary data information for tracers.
  end type ocean_obc_type

  type(ocean_obc_type), save :: Obc

  ! <INTERFACE NAME="ocean_obc_freesurf">
  !   <DESCRIPTION>
  !      Extrapolate surface height on the open boundaries for ocean model.
  !   </DESCRIPTION>
  !   <IN NAME="taum1, tau, taup1" >
  !      Time step index.
  !   </IN>
  !   <INOUT NAME="eta" >
  !      surface height
  !   </INOUT>
  !   <INOUT NAME="etataum1" >
  !      surface height at time step taum1
  !   </INOUT>
  interface ocean_obc_freesurf
     module procedure ocean_obc_freesurf_dtfs
     module procedure ocean_obc_freesurf_dteta
  end interface
  !</INTERFACE>

  !<INTERFACE NAME="obc_update_boundary">
  !   <DESCRIPTION>
  !   update field on the halo points at the global boundaries.
  !   </DESCRIPTION>
  !   <INOUT NAME="field">
  !      field to be update on the boundary
  !   </INOUT>
  interface ocean_obc_update_boundary
     module procedure ocean_obc_update_boundary_2d
     module procedure ocean_obc_update_boundary_3d
     module procedure ocean_obc_update_boundary_4d
  end interface
  !</INTERFACE>

  !--- namelist interface ---------------------------------------------
  !<NAMELIST NAME="ocean_obc_nml">
  !  <DATA NAME="nobc" TYPE="integer" >
  !    number of open boundary condition. Its value should be less than max_obc. Increase max_obc if needed.
  !  </DATA>
  !  <DATA NAME="direction" TYPE="character(len=10), dimension(max_obc)" >
  !    open boundary direction. Each element value should be west, east, south or north.
  !  </DATA>
  !  <DATA NAME="is, ie, js, je" TYPE="integer, dimension(max_obc)" >
  !    open boundary position. 
  !  </DATA>
  !  <DATA NAME="name" TYPE="character(len=32), dimension(max_obc)" >
  !    type of open bounday.
  !  </DATA>
  !  <DATA NAME="obc_relax_eta" TYPE="logical, dimension(max_obc)" >
  !    logical variable that decide whether relax eta or not. Default value is .false.
  !  </DATA>
  !  <DATA NAME="obc_consider_convu" TYPE="logical, dimension(max_obc)" >
  !    logical variable that decide whether to account for one 
  !    component of convu within the boundary. The appropriate behavior
  !    depends on the model resolution.
  !    Default value is .false.
  !  </DATA>
  !  <DATA NAME="obc_relax_eta_profile" TYPE="logical, dimension(max_obc)" >
  !    logical variable that decide whether relax eta to a prscribed profile or not. 
  !    Default value is .false.
  !  </DATA>
  !  <DATA NAME="obc_relax_tracer" TYPE="logical, dimension(max_obc)" >
  !    logical variable that decide whether relax tracer or not. Default value is .false.
  !  </DATA>
  !  <DATA NAME="obc_tracer_orlanski" TYPE="logical, dimension(max_obc)" >
  !    logical variable that decide whether apply orlanski obc on tracer or not. Default value is .false.
  !  </DATA>
  !  <DATA NAME="" TYPE="real, dimension(max_obc)" >
  !    Relaxation coefficient. 
  !  </DATA>
  !  <DATA NAME="adjust_topog" TYPE="logical">
  !    Control adjusting topography on the open boundary. Default value is true. If true, no topographic 
  !    gradient accross boundary.
  !  </DATA>
  !  <DATA NAME="nt" TYPE="integer">
  !    number of tracers to use open boundary condition
  !  </DATA>
  !  <DATA NAME="debug_phase_speed" TYPE="logical">
  !  Includes the phase speed into the model output.
  !  </DATA> 
  !  <DATA NAME="debug_obc" TYPE="logical">
  !  For debugging.
  !  </DATA> 
  !</NAMELIST>

  integer                               :: nobc= 0
  character(len=10), dimension(max_obc) :: direction='' ! open directions; to be consistent with is, ie, js, je
  integer, dimension(max_obc)           :: is=-999, ie=-999, js=-999, je=-999     ! boundary position
  character(len=32), dimension(max_obc) :: name = ''
  logical, dimension(max_obc)           :: obc_relax_eta = .FALSE.
  logical, dimension(max_obc)           :: obc_relax_eta_profile = .FALSE.
  logical, dimension(max_obc)           :: obc_consider_convu = .FALSE.
  logical, dimension(max_obc)           :: obc_vert_advel_t= .TRUE.
  logical, dimension(max_obc)           :: obc_vert_advel_u= .TRUE.
  logical, dimension(max_obc)           :: obc_adjust_forcing_fs = .FALSE.
  real,    dimension(max_obc)           :: rel_coef_eta = 0.0               ! relaxation coefficient, specify in namelist
  real,    dimension(max_obc)           :: ctrop_max = 1.5                  ! will be multiplied with +/- sqrt(g*H)
  real,    dimension(max_obc)           :: ctrop_min = 1.0                  ! will be multiplied with +/- sqrt(g*H)
  character(len=128), dimension(max_obc):: filename_eta = ''
  character(len=32),  dimension(max_obc):: fieldname_eta = ''
  logical, dimension(max_obc)           :: obc_relax_tracer = .FALSE.
  logical, dimension(max_obc)           :: obc_tracer_orlanski = .FALSE.
  real,               dimension(max_obc,max_prog_tracers):: rel_coef_tracer = 0.0     ! relaxation coefficient
  character(len=128), dimension(max_obc,max_prog_tracers):: filename_tracer
  character(len=32),  dimension(max_obc,max_prog_tracers):: fieldname_tracer
  logical                               :: adjust_topog = .TRUE.
  logical                               :: debug_obc = .FALSE.
  logical                               :: debug_phase_speed = .FALSE.
  integer, parameter                    :: mobcm1 = max_obc-1
  integer, parameter                    :: ntm2   = max_prog_tracers-2
  integer, parameter                    :: tobcm2 = max_obc*(max_prog_tracers-2)
  integer                               :: ne, te
  
  integer                               :: id_obc
  
  data (name(ne), ne=1,max_obc)          /'test_obc', mobcm1*'none'/
  data (fieldname_eta(ne), ne=1,max_obc) /'eta_t', mobcm1*'none'/
  data (filename_eta(ne), ne=1,max_obc)  /'obc_eta_t.nc', mobcm1*'none'/
  data ((fieldname_tracer(ne,te),  ne=1, max_obc), te=1, max_prog_tracers) &
           /max_obc*'temp_obc', max_obc*'salt_obc', tobcm2*'none'/
  data ((filename_tracer(ne,te),   ne=1, max_obc), te=1, max_prog_tracers) &
           /max_obc*'INPUT/obc_tr.nc', max_obc*'INPUT/obc_tr.nc' ,tobcm2*'none'/

  namelist /ocean_obc_nml/ nobc, direction, is, ie, js, je, name, &
                           obc_relax_eta, obc_relax_eta_profile,  rel_coef_eta,  &
                           ctrop_max, ctrop_min, &
                           filename_eta, fieldname_eta, &
                           obc_consider_convu, obc_adjust_forcing_fs, &
                           obc_vert_advel_t, obc_vert_advel_u, &
                           obc_relax_tracer, obc_tracer_orlanski, &
                           rel_coef_tracer, filename_tracer, fieldname_tracer,   &
                           adjust_topog, debug_phase_speed, debug_obc

  character(len=128)   :: version = '$ID$'
  character (len=128)  :: tagname = '$Name$'
  logical              :: module_is_initialized = .FALSE.
  integer              :: nt  = 0                         ! number of tracers
  integer              :: id_transport
  integer, allocatable :: id_tracer_flux(:)
  integer              :: id_ctrop_p, id_ctrop_a, id_cclin
  real, allocatable    :: wrk(:,:,:)      ! needed for output of phase speed

contains

  !#######################################################################
  !  <SUBROUTINE NAME="ocean_obc_init" >
  !   <DESCRIPTION>
  !      Allocates space and initializes a derived-type variable that 
  !      contains domain decompostion and grid information.
  !   </DESCRIPTION>
  !   <IN NAME="dtts, dtuv, dtfs, dteta" TYPE="real"></IN>
  !    time step.
  !   <IN NAME="Domain" TYPE="type(ocean_domain_type)" >
  !      A derived data type that contains domain information for mom4.
  !   </IN>
  !   <IN NAME="Grid" TYPE="type(ocean_grid_type)" >
  !      A derived data type that contains grid information for mom4.
  !   </IN>
  !   <INOUT NAME="have_obc" TYPE="logical" >
  !      logical variable to indicate if there is any open boundary condition. 
  !      if true, open boudanry exists. 
  !   </INOUT>
  !<PUBLICROUTINE>
  subroutine ocean_obc_init(have_obc, Time, Time_steps, Domain, Grid, debug)
  !</PUBLICROUTINE>

    logical, intent(inout)                       :: have_obc
    type(ocean_time_type), intent(in)            :: Time
    type(ocean_time_steps_type), intent(in)      :: Time_steps
    type (ocean_domain_type),intent(in),  target :: Domain
    type(ocean_grid_type), intent(inout), target :: Grid
    logical, intent(in), optional                :: debug

    !--- some local variables ------------------------------------------
    integer :: m, n, i, j, k, unit, io_status, ierr, ioff, joff

    if ( module_is_initialized ) then 
       call mpp_error(FATAL, '==>Error in ocean_obc_mod (ocean_obc_init): module already initialized')
    endif

    module_is_initialized = .TRUE.

    Dom => Domain
    Grd => Grid

    if (PRESENT(debug)) debug_obc = debug

    id_obc  = mpp_clock_id('(Ocean open boundaries) ',grain=CLOCK_MODULE)
    


    !--- read namelist and write out namelist --------------------------
    unit = open_namelist_file()
    read  (unit, ocean_obc_nml,iostat=io_status)
    write (stdout(),'(/)')
    write (stdout(), ocean_obc_nml)
    write (stdlog(), ocean_obc_nml)
    ierr = check_nml_error(io_status,'ocean_obc_nml')
    call close_file (unit)

    !--- write out version information ---------------------------------
    call write_version_number( version, tagname )

    !--- if there is no open boundary, just return
    Obc%nobc = nobc
    if(Obc%nobc == 0) then
       return
    endif
    have_obc = .true.

    call mpp_clock_begin(id_obc)
    !--- number of boundaries should not be over the max_obc
    if(nobc .gt. max_obc) call mpp_error(FATAL,'==>Error in ocean_obc_mod: nml nobc is greater than maximum'// &
         ' allowed bounds max_obc. Modify nobc or increase "max_obc" at top of ocean_obc_mod' )

    !--- cyclic condition and open boundary condition along zonal direction 
    !--- these two boundary conditions cannot exist at the same place 
    do m = 1, nobc
       if((trim(direction(m)) == 'west' .or. trim(direction(m)) == 'east')  .and. Grid%cyclic) &
            call mpp_error(FATAL, "==>Error in ocean_obc_mod: when west or east boundary is open, cyclic condition cannot exist")
    enddo

    !--- get the domain decomposition -----------------------------------
    call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
    call get_domain_offsets(Domain, ioff, joff)
    do m = 1, nobc
       is(m) = is(m) - ioff
       ie(m) = ie(m) - ioff
       js(m) = js(m) - joff
       je(m) = je(m) - joff
    enddo

    !--- get time step -------------------------------------------------
    dtts    = Time_steps%dtts
    dtuv    = Time_steps%dtuv
    dtfs    = Time_steps%dtfs
    dteta   = Time_steps%dteta
    dtime_t = Time_steps%dtime_t
    dtime_e = Time_steps%dtime_e

    ! number of vertical grid points 
    nk      = size(Grid%zw)

    !--- Initialize Obc data type --------------------------------------
    allocate(Obc%bound (nobc))

    Obc%bound%is = -999; Obc%bound%ie = -999    ! dummy number to indicate not on the boundary
    Obc%bound%js = -999; Obc%bound%je = -999
    Obc%bound%on_bound = .FALSE.
    do m = 1, nobc
       Obc%bound(m)%name                  = trim(name(m))
       select case(trim(direction(m)))
       case('west')
          Obc%bound(m)%direction = WEST
       case('east')
          Obc%bound(m)%direction = EAST
       case('south')
          Obc%bound(m)%direction = SOUTH
       case('north')
          Obc%bound(m)%direction = NORTH
       case default
          call mpp_error(FATAL,'each element of nml direction should be west, east, south or north')
       end select

       Obc%bound(m)%ctrop_max             = ctrop_max(m)
       Obc%bound(m)%ctrop_min             = ctrop_min(m)
       Obc%bound(m)%obc_consider_convu    = obc_consider_convu(m)
       Obc%bound(m)%obc_adjust_forcing_fs = obc_adjust_forcing_fs(m)
       Obc%bound(m)%obc_vert_advel_t      = obc_vert_advel_t(m)
       Obc%bound(m)%obc_vert_advel_u      = obc_vert_advel_u(m)
       Obc%bound(m)%obc_relax_eta         = obc_relax_eta(m)
       Obc%bound(m)%obc_relax_eta_profile = obc_relax_eta_profile(m)
       Obc%bound(m)%obc_relax_tracer      = obc_relax_tracer(m)
       Obc%bound(m)%obc_tracer_orlanski   = obc_tracer_orlanski(m)
       if (Obc%bound(m)%ctrop_max .lt. Obc%bound(m)%ctrop_min) then
         call mpp_error(FATAL,' ctrop_max <  ctrop_min for open boundary '//trim(Obc%bound(m)%name))
       endif
       
       if(Obc%bound(m)%direction == WEST .or. Obc%bound(m)%direction == EAST) then
          !--- define mask to indicate if it is open or not.
          allocate(Obc%bound(m)%mask(jsd:jed))
          Obc%bound(m)%mask = .FALSE.
          if(is(m) .le. iec .and. is(m) .ge. isc) then
             do j = jsd, jed
                if(j .ge. js(m) .and. j .le. je(m)) Obc%bound(m)%mask(j) = .TRUE.
             enddo
             if((jsc .ge. js(m) .and. jsc .le. je(m)) .or. (jec .ge. js(m) .and. jec .le. je(m) ) ) then
                Obc%bound(m)%on_bound = .true.
                Obc%bound(m)%is       = is(m)
                Obc%bound(m)%ie       = is(m)
                Obc%bound(m)%js       = max(js(m),jsc)
                Obc%bound(m)%je       = min(je(m),jec)
                Obc%bound(m)%np       = je(m) - js(m) + 1
                Obc%bound(m)%isd      = Obc%bound(m)%is
                Obc%bound(m)%ied      = Obc%bound(m)%ie
                Obc%bound(m)%jsd      = Obc%bound(m)%js
                Obc%bound(m)%jed      = Obc%bound(m)%je
                if(Obc%bound(m)%js > js(m) ) Obc%bound(m)%jsd = max(js(m),jsd)
                if(Obc%bound(m)%je < je(m) ) Obc%bound(m)%jed = min(je(m),jed)
                allocate(Obc%bound(m)%ctrop(Obc%bound(m)%js:Obc%bound(m)%je))
                allocate(Obc%bound(m)%pgrad(Obc%bound(m)%js:Obc%bound(m)%je))
             endif
          endif
       else    ! north or south direction
          allocate(Obc%bound(m)%mask(isd:ied))
          Obc%bound(m)%mask = .FALSE.
          if(js(m) .le. jec .and. js(m) .ge. jsc) then
             do i = isd, ied
                if(i .ge. is(m) .and. i .le. ie(m)) Obc%bound(m)%mask(i) = .TRUE.
             enddo
             if((isc .ge. is(m) .and. isc .le. ie(m)) .or. (iec .ge. is(m) .and. iec .le. ie(m) ) )then
                Obc%bound(m)%on_bound = .true.
                Obc%bound(m)%js       = js(m)
                Obc%bound(m)%je       = js(m)
                Obc%bound(m)%is       = max(is(m),isc)
                Obc%bound(m)%ie       = min(ie(m),iec)
                Obc%bound(m)%np       = ie(m) - is(m) + 1
                Obc%bound(m)%isd      = Obc%bound(m)%is
                Obc%bound(m)%ied      = Obc%bound(m)%ie
                Obc%bound(m)%jsd      = Obc%bound(m)%js
                Obc%bound(m)%jed      = Obc%bound(m)%je
                if(Obc%bound(m)%is > is(m) ) Obc%bound(m)%isd = max(is(m),isd)
                if(Obc%bound(m)%ie < ie(m) ) Obc%bound(m)%ied = min(ie(m),ied)
                allocate(Obc%bound(m)%ctrop(Obc%bound(m)%is:Obc%bound(m)%ie))
                allocate(Obc%bound(m)%pgrad(Obc%bound(m)%is:Obc%bound(m)%ie))
             endif
          endif
       endif
    enddo

    write(stdout(),*) '-----------------------------------------------------------------'
    write(stdout(),*) 'The following setup for OBC has been found:'
    write(stdout(),*) 'Total number of OBC: ', nobc
    do m = 1, nobc
       write(stdout(),*) ' Setup of OBC ',m,', ',Obc%bound(m)%name,':'
       write(stdout(),*) '  direction     : ', trim(direction(m))
       write(stdout(),*) '  indizees      :', Obc%bound(m)%is, Obc%bound(m)%ie, Obc%bound(m)%js, Obc%bound(m)%je
       write(stdout(),*) '  points        :', Obc%bound(m)%np
       if (Obc%bound(m)%obc_vert_advel_t) then 
          write(stdout(),*) ' with vertical tracer advection'
       else  
          write(stdout(),*) ' no vertical tracer advection'
       endif
       if (Obc%bound(m)%obc_vert_advel_u) then 
          write(stdout(),*) ' with vertical momentum advection'
       else  
          write(stdout(),*) ' no vertical momentum advection'
       endif
       if (Obc%bound(m)%obc_adjust_forcing_fs) then
          write(stdout(),*) ' remove baroclinic cross boundary pressure gradients from'
          write(stdout(),*) ' vertically integrated forcing'
       else  
          write(stdout(),*) ' do not remove baroclinic cross boundary pressure gradients from'
          write(stdout(),*) ' vertically integrated forcing'
       endif
       if (Obc%bound(m)%obc_consider_convu) then 
          write(stdout(),*) ' consider convu for d eta/dt'
       else  
          write(stdout(),*) ' do not consider convu for d eta/dt'
       endif
       write(stdout(),*) '  relax eta     :', Obc%bound(m)%obc_relax_eta
       if (Obc%bound(m)%obc_relax_eta_profile) then
          write(stdout(),*) '  method        :', ' profile'
       else  
          write(stdout(),*) '  method        :', ' average'
       endif
       write(stdout(),*) '  relax tracer  :', Obc%bound(m)%obc_relax_tracer
       write(stdout(),*) '  t - Orlanski  :', Obc%bound(m)%obc_tracer_orlanski
    enddo
    write(stdout(),*) '-----------------------------------------------------------------'

    !--- adjust topography at the open bounday
    if(adjust_topog) call ocean_obc_adjust_topog(Grid%ht, Grid%hu, Grid%kmt, Grid%kmu)
    
    call ocean_obc_set_mask

    
    
    call mpp_clock_end(id_obc)
    
    return

  end subroutine ocean_obc_init
  !  </SUBROUTINE>  

  !#######################################################################
  !  <SUBROUTINE NAME="ocean_obc_tracer_init" >
  !   <DESCRIPTION>
  !      Allocates space and initializes all stuff for tracers at OBC
  !   </DESCRIPTION>
  !   <IN NAME="debug" TYPE="logical"></IN>
  !<PUBLICROUTINE>
  subroutine ocean_obc_tracer_init(Time, T_prog, num_prog_tracers, debug)
    !</PUBLICROUTINE>
    type(ocean_time_type), intent(in)                   :: Time
    type(ocean_prog_tracer_type), intent(inout), target :: T_prog(:)
    integer, intent(in)                                 :: num_prog_tracers
    logical, intent(in),        optional                :: debug
    integer :: fld_size(4)

    !--- some local variables ------------------------------------------
    integer :: m, n, taum1, tau, taup1
    
    if (PRESENT(debug)) debug_obc = (debug.or.debug_obc)

    !   Allocate space to store tracer flux through open boundary       
    !--- get the number of prog tracers
    nt = num_prog_tracers
    
    if(debug_obc) write(stdout(),*) 'Using ', nt,' tracers in OBC'
    
    if(nt .gt. max_prog_tracers) call mpp_error(FATAL,'==>Error in ocean_obc_mod: num_prog_tracers is greater than '//&
         'maximum number of prog tracers allocated in obc. For more tracers, increase "max_prog_tracers" at top of ocean_obc_mod')

    if(nt == 0) then 
       call mpp_error(NOTE,'==>NOTE in ocean_obc_mod: number of prognostics tracers is 0')
       return
    endif

    call mpp_clock_begin(id_obc)    

    do n=1, nt
       allocate(T_prog(n)%otf(nobc))
    enddo
    do m = 1, nobc
       if(.not. Obc%bound(m)%on_bound) cycle

       select case(Obc%bound(m)%direction )
       case(WEST)
          do n=1, nt
             allocate(T_prog(n)%otf(m)%flux(Obc%bound(m)%js:Obc%bound(m)%je,nk))
          enddo
       case(EAST)
          do n=1, nt
             allocate(T_prog(n)%otf(m)%flux(Obc%bound(m)%js:Obc%bound(m)%je,nk))
          enddo
       case(SOUTH)
          do n=1, nt
             allocate(T_prog(n)%otf(m)%flux(Obc%bound(m)%is:Obc%bound(m)%ie,nk))
          enddo
       case(NORTH)
          do n=1, nt
             allocate(T_prog(n)%otf(m)%flux(Obc%bound(m)%is:Obc%bound(m)%ie,nk))
          enddo
       end select
    enddo
!   allocate space to store advection across the boundary
    do m = 1, nobc
       if(.not. Obc%bound(m)%on_bound) cycle
      
       select case( Obc%bound(m)%direction )
       case(WEST)
         allocate(Obc%bound(m)%hadv(Obc%bound(m)%js:Obc%bound(m)%je,nk))
       case(EAST)
         allocate(Obc%bound(m)%hadv(Obc%bound(m)%js:Obc%bound(m)%je,nk))
       case(SOUTH)
         allocate(Obc%bound(m)%hadv(Obc%bound(m)%is:Obc%bound(m)%ie,nk))
       case(NORTH)
         allocate(Obc%bound(m)%hadv(Obc%bound(m)%is:Obc%bound(m)%ie,nk))
       end select
    enddo
    
    write(stdout(),*) '-----------------------------------------------------------------'
    write(stdout(),*) 'The following setup for time interpolation of OBC data has been found:'

    allocate(Obc%tracer(nobc,nt))
    allocate(Obc%eta   (nobc))  
    do m = 1, nobc
       if (Obc%bound(m)%on_bound) then
          !         data needed for eta
          Obc%eta(m)%rel_coef                   = rel_coef_eta(m)
          Obc%eta(m)%file                       = trim(filename_eta(m))
          Obc%eta(m)%field                      = trim(fieldname_eta(m))
          !         data needed for tracer
          do n=1, nt
             Obc%tracer(m,n)%file                = trim(filename_tracer(m,n))
             Obc%tracer(m,n)%field               = trim(fieldname_tracer(m,n))
             Obc%tracer(m,n)%rel_coef            = rel_coef_tracer(m,n)
             write(stdout(),*) trim(Obc%tracer(m,n)%file)
             write(stdout(),*) trim(Obc%tracer(m,n)%field)
          enddo
          if(Obc%bound(m)%direction == WEST .or. Obc%bound(m)%direction == EAST) then
             if(Obc%bound(m)%obc_relax_eta) then 
                allocate(Obc%eta(m)%data(1,Obc%bound(m)%js:Obc%bound(m)%je))
                Obc%eta(m)%id   = init_external_field( &
                     Obc%eta(m)%file,Obc%eta(m)%field,desired_units='m', &
                     verbose=.true.)
                fld_size        = get_external_field_size(Obc%eta(m)%id)
                if(fld_size(2) .ne. size(Obc%eta(m)%data,2)) then
                   call mpp_error(FATAL,'invalid dimension 2 of input field in '//trim(Obc%eta(m)%file))
                endif
             endif
             !           allocate variables for time interpolation of external data
             !           if relaxation or upstream advection is switched on
             if(Obc%bound(m)%obc_relax_tracer.or.(.not.Obc%bound(m)%obc_tracer_orlanski)) then 
                do n=1, nt
                   allocate(Obc%tracer(m,n)%data(1,Obc%bound(m)%js:Obc%bound(m)%je,nk))
                   Obc%tracer(m,n)%id = init_external_field(Obc%tracer(m,n)%file,Obc%tracer(m,n)%field,verbose=.true.)
                   fld_size = get_external_field_size(Obc%tracer(m,n)%id)
                   if(fld_size(2) .ne. size(Obc%tracer(m,n)%data,2)) then
                      call mpp_error(FATAL,'invalid dimension 2 of input field in '//trim(Obc%tracer(m,n)%file))
                   endif
                   if(fld_size(3) .ne. size(Obc%tracer(m,n)%data,3)) then
                      call mpp_error(FATAL,'invalid dimension 3 of input field in '//trim(Obc%tracer(m,n)%file))
                   endif
                enddo
             endif
          else
             if(Obc%bound(m)%obc_relax_eta) then 
                allocate(Obc%eta(m)%data(Obc%bound(m)%is:Obc%bound(m)%ie,1))
                Obc%eta(m)%id   = init_external_field( &
                     Obc%eta(m)%file,Obc%eta(m)%field,desired_units='m', &
                     verbose=.true.)
                fld_size = get_external_field_size(Obc%eta(m)%id)
                if(fld_size(1) .ne. size(Obc%eta(m)%data,1)) then
                   call mpp_error(FATAL,'invalid dimension 1 of input field in '//trim(Obc%eta(m)%file))
                endif
             endif
             !           allocate variables for time interpolation of external data
             !           if relaxation or upstream advection is switched on
             if(Obc%bound(m)%obc_relax_tracer.or.(.not.Obc%bound(m)%obc_tracer_orlanski)) then 
                do n=1, nt
                   allocate(Obc%tracer(m,n)%data(Obc%bound(m)%is:Obc%bound(m)%ie,1,nk))
                   Obc%tracer(m,n)%id = init_external_field(Obc%tracer(m,n)%file,Obc%tracer(m,n)%field,verbose=.true.)
                   fld_size = get_external_field_size(Obc%tracer(m,n)%id)
                   if(fld_size(1) .ne. size(Obc%tracer(m,n)%data,1)) then
                      call mpp_error(FATAL,'invalid dimension 1 of input field in '//trim(Obc%tracer(m,n)%file))
                   endif
                   if(fld_size(3) .ne. size(Obc%tracer(m,n)%data,3)) then
                      call mpp_error(FATAL,'invalid dimension 3 of input field in '//trim(Obc%tracer(m,n)%file))
                   endif
                enddo
             endif
          endif
          write(stdout(),*) '-----------------------------------------------------------------'
          write(stdout(),*) 'The following scheme for OBC data input has been found:'
          if (Obc%bound(m)%obc_relax_eta) then
             write(stdout(),*) '  sea level relaxation'
             write(stdout(),*) '  rel_coef  :', Obc%eta(m)%rel_coef 
             write(stdout(),*) '  input-file:', trim(Obc%eta(m)%file) 
             write(stdout(),*) '  input-name:', trim(Obc%eta(m)%field)
          else
             write(stdout(),*) '  no sea level relaxation'
          endif
          if (Obc%bound(m)%obc_relax_tracer) then
             write(stdout(),*) '  tracer relaxation'
             write(stdout(),*) '  relax tracer: nb rel_coef input-file       input-name'
             do n=1, nt
                write(stdout(),'(17x,I2.2,1x,F8.5,2x,a,2x,a)') n, Obc%tracer(m,n)%rel_coef,& 
                     trim(obc%tracer(m,n)%file), trim(obc%tracer(m,n)%field)
             enddo
          else
             write(stdout(),*) '  no tracer relaxation'
          endif
       endif
    enddo

    ! register diagnostic output 
    allocate(id_tracer_flux(nt))
    do n = 1, nt
       if(trim(T_prog(n)%name) == 'temp') then
          id_tracer_flux(n) = register_diag_field ('ocean_model', 'obc_'//trim(T_prog(n)%name)//'_adv_flux_int_z',  &
               Grd%tracer_axes(1:2), Time%model_time, 'z-integral of horizontal tracer advection flux on obc',      &
               'Watt/m', missing_value=-1.e20, range=(/-1.e20,1.e20/))
       else
          id_tracer_flux(n) = register_diag_field ('ocean_model', 'obc_'//trim(T_prog(n)%name)//'_adv_flux_int_z', &
               Grd%tracer_axes(1:2), Time%model_time, 'z-integral of horizontal tracer advection flux on obc',     &
               &' kg/m/sec', missing_value=-1.e20, range=(/-1.e20,1.e20/))
       endif
    enddo

    id_transport = register_diag_field ('ocean_model', 'obc_transport_int_z',                 &
         Grd%vel_axes_uv(1:2), Time%model_time, 'z-integral of horizontal transport on obc',  &
         'm^2/s', missing_value=-1.e20, range=(/-1.e20,1.e20/))

!  Set tracer poutside the boundaries to defined values to avoid artificial diffusion   
    taum1 = Time%taum1
    tau   = Time%tau
    taup1 = Time%taup1
    
    do n = 1, nt
      call ocean_obc_update_boundary(T_prog(n)%field(:,:,:,taum1), 'T')
      call ocean_obc_update_boundary(T_prog(n)%field(:,:,:,tau), 'T')
      call ocean_obc_update_boundary(T_prog(n)%field(:,:,:,taup1), 'T')
    enddo
    
    call mpp_clock_end(id_obc)
    if (debug_phase_speed) then
       id_ctrop_p = register_diag_field ('ocean_model', 'ctrop_p',                     &
            Grd%tracer_axes(1:2), Time%model_time, 'barotr phase speed on open bounds',  &
            'm/s', missing_value=-1.e20, range=(/-1.e20,1.e20/))
       id_ctrop_a = register_diag_field ('ocean_model', 'ctrop_a',                     &
            Grd%tracer_axes(1:2), Time%model_time, 'barotr phase speed on open bounds',  &
            'm/s', missing_value=-1.e20, range=(/-1.e20,1.e20/))
       allocate(wrk(isd:ied,jsd:jed,nk))
    endif
    return

  end subroutine ocean_obc_tracer_init
  !  </SUBROUTINE>

  !#####################################################################
  !  <SUBROUTINE NAME="ocean_obc_prepare" >
  !   <DESCRIPTION>
  !      Prepares OBC  
  !      
  !   </DESCRIPTION>
  !<PUBLICROUTINE>

  subroutine ocean_obc_prepare(Time, Ext_mode, T_prog)
    type(ocean_time_type), intent(in)           :: Time 
    type(ocean_external_mode_type),  intent(in) :: Ext_mode
    type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)

    integer                                     :: m, n, i, j, taum1, tau
    logical                                     :: used
  !</PUBLICROUTINE>

    call mpp_clock_begin(id_obc)

!   prepare the data needed for time interpolation of external data.
    do m = 1, nobc

       !--- if on current pe there is no point on the bound, then just return
       if(.not. Obc%bound(m)%on_bound) cycle

!      Prepare sea level data for relaxation         
       if(Obc%bound(m)%obc_relax_eta) then
         
         call time_interp_external(Obc%eta(m)%id,Time%model_time,Obc%eta(m)%data,verbose=.true.)

       endif
       do n=1, nt
         if(Obc%bound(m)%obc_relax_tracer) then         
          
           call time_interp_external(Obc%tracer(m,n)%id,Time%model_time,Obc%tracer(m,n)%data,verbose=.true.)

         endif
       enddo

       !--- set the tracer_flux to 0
       do n =1, nt
         T_prog(n)%otf(m)%flux(:,:) = 0.0
      enddo

    enddo

!   Now initialize the barotropic phase speed smoothing    
    
    taum1=Time%taum1; tau=Time%tau
    call barotropic_phase_speed_init(Ext_mode%eta_t, taum1, tau)
    
    if(id_ctrop_p > 0) then

      wrk = 0.
      do m = 1, nobc       
       if(.not. Obc%bound(m)%on_bound) cycle

       if(Obc%bound(m)%direction == WEST .or. Obc%bound(m)%direction == EAST) then
         do j = Obc%bound(m)%js, Obc%bound(m)%je
           wrk(Obc%bound(m)%is,j,1) = Obc%bound(m)%ctrop(j)
         enddo
       endif
       if(Obc%bound(m)%direction == SOUTH .or. Obc%bound(m)%direction == NORTH) then
         do i = Obc%bound(m)%is, Obc%bound(m)%ie
           wrk(i,Obc%bound(m)%js,1) = Obc%bound(m)%ctrop(i)
         enddo
       endif
      enddo

      used = send_data(id_ctrop_p, wrk(:,:,1), &
             Time%model_time, rmask=Grd%tmask(:,:,1), &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    endif

    call mpp_clock_end(id_obc)

    return

  end subroutine ocean_obc_prepare
  !  </SUBROUTINE>

  !#####################################################################
  !update surface height and the vertically integrated horizontal velocity
  !<SUBROUTINE NAME="ocean_obc_freesurf_dtfs" INTERFACE="ocean_obc_freesurf">
  !   <IN NAME="taum1, tau, taup1" TYPE="integer"> </IN>
  !   <INOUT NAME="eta" TYPE="real, dimension(isd:,jsd:,:)"></INOUT>

  !<PUBLICROUTINE INTERFACE="ocean_obc_freesurf" >
  subroutine ocean_obc_freesurf_dtfs(eta, taum1, tau, taup1, tstep)
    !</PUBLICROUTINE>

    real, dimension(isd:,jsd:,:),   intent(inout) :: eta
    integer, intent(in)                           :: taum1, tau, taup1
    real, intent(in)                              :: tstep

    real    :: cdtdxr, cdtdyr, etasum
    integer :: i, j, m, taust

    call mpp_clock_begin(id_obc)

    !--- calculate the phase speed -------------------------------------
    call barotropic_phase_speed(eta, taum1, tau, taup1, tstep)  

    do m = 1, nobc

       !--- if on current pe there is no point on the bound, then just return
       if(.not. Obc%bound(m)%on_bound) cycle

       if(Obc%bound(m)%obc_consider_convu) then
         taust=taup1
       else
         taust=taum1
       endif
       
       select case( Obc%bound(m)%direction )
       case(WEST)
          i =  Obc%bound(m)%is
          do j = Obc%bound(m)%js, Obc%bound(m)%je
             cdtdxr = Obc%bound(m)%ctrop(j)*tstep/Grd%dxu(i,j)
             eta(i,j,taup1) = (eta(i,j,taust) - cdtdxr*eta(i+1,j,taup1))/(1. - cdtdxr)
          enddo
          if(Obc%bound(m)%obc_relax_eta) then
            if(Obc%bound(m)%obc_relax_eta_profile) then
 !            Relax towards a prescribed profile 
              do j = Obc%bound(m)%js, Obc%bound(m)%je
                 eta(i,j,taup1) = eta(i,j,taup1) + (Obc%eta(m)%data(1,j) - eta(i,j,taup1)) & 
                                * Obc%eta(m)%rel_coef * tstep
              enddo
            else
!             Relax mean sea level leaving the gradient profile in the boundary unchanged         
              etasum = sum(eta(i+1,Obc%bound(m)%js:Obc%bound(m)%je,taup1))/Obc%bound(m)%np
              do j = Obc%bound(m)%js, Obc%bound(m)%je
                eta(i,j,taup1) = eta(i,j,taup1) + (Obc%eta(m)%data(1,j) - etasum) & 
                               * Obc%eta(m)%rel_coef * tstep
              enddo
            endif
          endif
       case (EAST)
          i =  Obc%bound(m)%is
          do j = Obc%bound(m)%js, Obc%bound(m)%je
             cdtdxr = Obc%bound(m)%ctrop(j)*tstep/Grd%dxu(i,j)
             eta(i,j,taup1) = (eta(i,j,taust) + cdtdxr*eta(i-1,j,taup1))/(1. + cdtdxr)
          enddo
          if(Obc%bound(m)%obc_relax_eta) then
            if(Obc%bound(m)%obc_relax_eta_profile) then
!             Relax towards a prescribed profile         
              do j = Obc%bound(m)%js, Obc%bound(m)%je
                 eta(i,j,taup1) = eta(i,j,taup1) + (Obc%eta(m)%data(1,j) - eta(i,j,taup1)) & 
                                * Obc%eta(m)%rel_coef * tstep
              enddo
            else
!             Relax mean sea level leaving the gradient profile in the boundary unchanged         
              etasum = sum(eta(i-1,Obc%bound(m)%js:Obc%bound(m)%je,taup1))/Obc%bound(m)%np
              do j = Obc%bound(m)%js, Obc%bound(m)%je
                eta(i,j,taup1) = eta(i,j,taup1) + (Obc%eta(m)%data(1,j) - etasum) & 
                               * Obc%eta(m)%rel_coef * tstep
              enddo
            endif
          endif
       case(SOUTH)
          j = Obc%bound(m)%js
          do i = Obc%bound(m)%is, Obc%bound(m)%ie
             cdtdyr = Obc%bound(m)%ctrop(i)*tstep/Grd%dyu(i,j)
             eta(i,j,taup1) = (eta(i,j,taust) - cdtdyr*eta(i,j+1,taup1))/(1. - cdtdxr)
          enddo
          if(Obc%bound(m)%obc_relax_eta) then
            if(Obc%bound(m)%obc_relax_eta_profile) then
!             Relax towards a prescribed profile         
              do i = Obc%bound(m)%is, Obc%bound(m)%ie
                eta(i,j,taup1) = eta(i,j,taup1) + (Obc%eta(m)%data(i,1) - eta(i,j,taup1)) & 
                               * Obc%eta(m)%rel_coef * tstep
              enddo
            else
!             Relax mean sea level leaving the gradient profile in the boundary unchanged         
              etasum = sum(eta(Obc%bound(m)%is:Obc%bound(m)%ie,j+1,taup1))/Obc%bound(m)%np
              do i = Obc%bound(m)%is, Obc%bound(m)%ie
                eta(i,j,taup1) = eta(i,j,taup1) + (Obc%eta(m)%data(i,1) - etasum) & 
                               * Obc%eta(m)%rel_coef * tstep
              enddo
            endif
          endif
       case(NORTH)
          j = Obc%bound(m)%js
          do i = Obc%bound(m)%is, Obc%bound(m)%ie
             cdtdyr = Obc%bound(m)%ctrop(i)*tstep/Grd%dyu(i,j)
             eta(i,j,taup1) = (eta(i,j,tau) + cdtdyr*eta(i,j-1,taup1))/(1. + cdtdxr)
          enddo
          if(Obc%bound(m)%obc_relax_eta) then
            if(Obc%bound(m)%obc_relax_eta_profile) then
!             Relax towards a prescribed profile         
              do i = Obc%bound(m)%is, Obc%bound(m)%ie
                eta(i,j,taup1) = eta(i,j,taust) + (Obc%eta(m)%data(i,1) - eta(i,j,taup1)) & 
                               * Obc%eta(m)%rel_coef * tstep
              enddo
            else  
!             Relax mean sea level leaving the gradient profile in the boundary unchanged         
              etasum = sum(eta(Obc%bound(m)%is:Obc%bound(m)%ie,j-1,taup1))/Obc%bound(m)%np
              do i = Obc%bound(m)%is, Obc%bound(m)%ie
                eta(i,j,taup1) = eta(i,j,taup1) + (Obc%eta(m)%data(i,1) - etasum) & 
                               * Obc%eta(m)%rel_coef * tstep
              enddo
            endif
          endif
       end select
    enddo

    !--- update eta at global halo point to make the gradient accross boundary is 0
    call ocean_obc_update_boundary(eta(:,:,taup1), 'T')    

    if(debug_obc) then
       write(stdout(),*) 'After ocean_obc_freesurf_dtfs , eta chksum ===>', mpp_chksum(eta(isc:iec,jsc:jec,taup1))
    endif

    call mpp_clock_end(id_obc)

    return
  end subroutine ocean_obc_freesurf_dtfs
  !</SUBROUTINE>
 
  !#####################################################################
  !<SUBROUTINE NAME="ocean_obc_freesurf_dteta" INTERFACE="ocean_obc_freesurf">
  !   <INOUT NAME="Ext_mode" TYPE="ocean_external_mode_type"></INOUT>
  !<PUBLICROUTINE INTERFACE="ocean_obc_freesurf" >
  subroutine ocean_obc_freesurf_dteta(Time, Ext_mode)
    type(ocean_time_type), intent(in)             :: Time
    type(ocean_external_mode_type), intent(inout) :: Ext_mode
    !</PUBLICROUTINE>
    integer                                       :: taum1, tau, taup1, taust

  ! Find also an approximation for Ext_mode%deta_dt at the boundary
  ! taum1 may be equal to taup1, hence, prev is needed

    real :: cdtdxr, cdtdyr, prev
    integer :: i, ii, j, m
    logical :: used

    call mpp_clock_begin(id_obc)

    taum1=Time%taum1;tau=Time%tau;taup1=Time%taup1

    !--- calculate eta at the leap-frog time step implicitly
    !    updating from taum1 to taup1 implies wave contribution only
    !    updating from taup1 to taup1 implies, that div U along the boundary is
    !    also taken into account 
    do m = 1, nobc

       !--- if on current pe there is no point on the bound, then just return
      if(.not.Obc%bound(m)%on_bound) cycle
       
      if(Obc%bound(m)%obc_consider_convu) then
         taust=taup1
      else
         taust=taum1
      endif
         
      select case( Obc%bound(m)%direction )
      case(WEST)
          i =  Obc%bound(m)%is
          do j = Obc%bound(m)%js, Obc%bound(m)%je
             cdtdxr = Obc%bound(m)%ctrop(j)*dtime_e/Grd%dxu(i,j)
!             prev   = Ext_mode%eta_t(i,j,taust)
             Ext_mode%eta_t(i,j,taup1) = (Ext_mode%eta_t(i,j,taust) - cdtdxr*Ext_mode%eta_t(i+1,j,taup1))/(1. - cdtdxr)
!             Ext_mode%deta_dt(i,j) = Ext_mode%deta_dt(i,j) + (Ext_mode%eta_t(i,j,taup1) - prev)/dtime_e
          enddo
      case(EAST)
          i =  Obc%bound(m)%is
          do j = Obc%bound(m)%js, Obc%bound(m)%je
             cdtdxr = Obc%bound(m)%ctrop(j)*dtime_e/Grd%dxu(i,j)
!             prev   = Ext_mode%eta_t(i,j,taust)
             Ext_mode%eta_t(i,j,taup1) = (Ext_mode%eta_t(i,j,taust) + cdtdxr*Ext_mode%eta_t(i-1,j,taup1))/(1. + cdtdxr)
!             Ext_mode%deta_dt(i,j) = Ext_mode%deta_dt(i,j) + (Ext_mode%eta_t(i,j,taup1) - prev)/dtime_e
          enddo
       case(SOUTH)
          j = Obc%bound(m)%js
          do i = Obc%bound(m)%is, Obc%bound(m)%ie
             cdtdyr = Obc%bound(m)%ctrop(i)*dtime_e/Grd%dyu(i,j)
!             prev   = Ext_mode%eta_t(i,j,taust)
             Ext_mode%eta_t(i,j,taup1) = (Ext_mode%eta_t(i,j,taust) - cdtdyr*Ext_mode%eta_t(i,j+1,taup1))/(1. - cdtdxr)
!             Ext_mode%deta_dt(i,j) = Ext_mode%deta_dt(i,j) + (Ext_mode%eta_t(i,j,taup1) - prev)/dtime_e
          enddo
       case(NORTH)
          j = Obc%bound(m)%js
          do i = Obc%bound(m)%is, Obc%bound(m)%ie
             cdtdyr = Obc%bound(m)%ctrop(i)*dtime_e/Grd%dyu(i,j)
!             prev   = Ext_mode%eta_t(i,j,taust)
             Ext_mode%eta_t(i,j,taup1) = (Ext_mode%eta_t(i,j,taust) + cdtdyr*Ext_mode%eta_t(i,j-1,taup1))/(1. + cdtdxr)
!             Ext_mode%deta_dt(i,j) = Ext_mode%deta_dt(i,j) + (Ext_mode%eta_t(i,j,taup1) - prev)/dtime_e
          enddo
       end select
    enddo

    !--- update eta at global halo point to make the gradient accross boundary is 0
    call ocean_obc_update_boundary(Ext_mode%eta_t(:,:,taup1), 'T')

    if(id_ctrop_a > 0) then

      wrk = 0.
      do m = 1, nobc       
       if(.not. Obc%bound(m)%on_bound) cycle

       if(Obc%bound(m)%direction == WEST .or. Obc%bound(m)%direction == EAST) then
         do j = Obc%bound(m)%js, Obc%bound(m)%je
           wrk(Obc%bound(m)%is,j,1) = Obc%bound(m)%ctrop(j)
         enddo
       endif
       if(Obc%bound(m)%direction == SOUTH .or. Obc%bound(m)%direction == NORTH) then
         do i = Obc%bound(m)%is, Obc%bound(m)%ie
           wrk(i,Obc%bound(m)%js,1) = Obc%bound(m)%ctrop(i)
         enddo
       endif
      enddo

      used = send_data(id_ctrop_a, wrk(:,:,1), &
             Time%model_time, rmask=Grd%tmask(:,:,1), &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
    endif

    if(debug_obc) then
       write(stdout(),*) 'After ocean_obc_freesurf_dteta , eta chksum ===>', &
                          mpp_chksum(Ext_mode%eta_t(isc:iec,jsc:jec,taup1))
    endif

    call mpp_clock_end(id_obc)

    return
  end subroutine ocean_obc_freesurf_dteta

  !</SUBROUTINE>

  !#####################################################################
  !<SUBROUTINE NAME="ocean_obc_adjust_advel">
  !   <DESCRIPTION>
  !      Subtract wrong vertical bottom velocity 
  !   </DESCRIPTION>
  !   <INOUT NAME="Adv_vel" TYPE="ocean_adv_vel_type" >
  !      Advection velocities 
  !   </INOUT>
  !<PUBLICROUTINE>
  subroutine ocean_obc_adjust_advel(Adv_vel)
    !</PUBLICROUTINE>

    type(ocean_adv_vel_type), intent(inout) :: Adv_vel

    integer :: i, j, k, m, kmt

    call mpp_clock_begin(id_obc)

    do m = 1, nobc
       !      If obc_vert_advel_t is true (default) calculate a reasonale approximation for 
       !      Adv_vel%w_bt. This will be used for Ext_mode%deta_dt and tracer
       !      advection as well. Adv_vel%w_bu will be calculated for baroclinic velocity.

       !      If obc_vert_advel_t is false (from namelist) calculate a reasonale approximation for
       !      Adv_vel%w_bt for Adv_vel%w_bu. Then set Adv_vel%w_bt to zero.

       !--- if on current pe there is no point on the bound, then just return
       if(Obc%bound(m)%on_bound) then
          !--- calculate the phase speed at west or east boundary -------------------
         if( Obc%bound(m)%obc_vert_advel_t .or. Obc%bound(m)%obc_vert_advel_u) then
           do j = Obc%bound(m)%jsd, Obc%bound(m)%jed
             do i = Obc%bound(m)%isd, Obc%bound(m)%ied
                kmt = Grd%kmt(i,j) 
                do k=0, kmt
                   Adv_vel%w_bt(i,j,k) = Adv_vel%w_bt(i,j,k) - Adv_vel%w_bt(i,j,kmt)   !*float(k)/float(kmt) 
                enddo
             enddo
           enddo
          endif
       endif

 !      call mpp_update_domains(Adv_vel%w_bt(:,:,:), Dom%domain2d)
       call ocean_obc_update_boundary(Adv_vel%w_bt(:,:,:), 'T')

       if(.not. Obc%bound(m)%on_bound) cycle

       select case( Obc%bound(m)%direction )
       case(WEST)
          i = Obc%bound(m)%is
          if( Obc%bound(m)%obc_vert_advel_u ) then
             do j = Obc%bound(m)%js, Obc%bound(m)%je
                do k=0, Grd%kmu(i,j)
                   Adv_vel%w_bu(i,j,k) = (Adv_vel%w_bt(i,j,k)*Grd%dte(i,j)*Grd%dus(i,j)     + &
                        Adv_vel%w_bt(i+1,j,k)*Grd%dtw(i+1,j)*Grd%dus(i,j) + &
                        Adv_vel%w_bt(i,j+1,k)*Grd%dte(i,j+1)*Grd%dun(i,j) + &
                        Adv_vel%w_bt(i+1,j+1,k)*Grd%dtw(i+1,j+1)*Grd%dun(i,j))*Grd%daur(i,j)
                enddo
             enddo
          else
             do j = Obc%bound(m)%jsd, Obc%bound(m)%jed
                Adv_vel%w_bu(i,j,:) = 0. 
             enddo
          endif
       case(EAST)
          i = Obc%bound(m)%is-1
          if( Obc%bound(m)%obc_vert_advel_u ) then
             do j = Obc%bound(m)%js, Obc%bound(m)%je
                do k=0, Grd%kmu(i,j)
                   Adv_vel%w_bu(i,j,k) = (Adv_vel%w_bt(i,j,k)*Grd%dte(i,j)*Grd%dus(i,j)     + &
                        Adv_vel%w_bt(i+1,j,k)*Grd%dtw(i+1,j)*Grd%dus(i,j) + &
                        Adv_vel%w_bt(i,j+1,k)*Grd%dte(i,j+1)*Grd%dun(i,j) + &
                        Adv_vel%w_bt(i+1,j+1,k)*Grd%dtw(i+1,j+1)*Grd%dun(i,j))*Grd%daur(i,j)
                enddo
             enddo
          else
             do j = Obc%bound(m)%jsd, Obc%bound(m)%jed
                Adv_vel%w_bu(i,j,:) = 0. 
             enddo
          endif
       case(SOUTH)
          j = Obc%bound(m)%js
          if( Obc%bound(m)%obc_vert_advel_u ) then
             do i = Obc%bound(m)%is, Obc%bound(m)%ie
                do k=0, Grd%kmu(i,j)
                   Adv_vel%w_bu(i,j,k) = (Adv_vel%w_bt(i,j,k)*Grd%dte(i,j)*Grd%dus(i,j)     +  &
                        Adv_vel%w_bt(i+1,j,k)*Grd%dtw(i+1,j)*Grd%dus(i,j) +  &
                        Adv_vel%w_bt(i,j+1,k)*Grd%dte(i,j+1)*Grd%dun(i,j) +  &
                        Adv_vel%w_bt(i+1,j+1,k)*Grd%dtw(i+1,j+1)*Grd%dun(i,j))*Grd%daur(i,j)
                enddo
             enddo
          else
             do j = Obc%bound(m)%jsd, Obc%bound(m)%jed
                Adv_vel%w_bu(i,j,:) = 0. 
             enddo
          endif
       case(NORTH)
          j = Obc%bound(m)%js-1
          if( Obc%bound(m)%obc_vert_advel_u ) then
             do i = Obc%bound(m)%is, Obc%bound(m)%ie
                do k=0, Grd%kmu(i,j)
                   Adv_vel%w_bu(i,j,k) = (Adv_vel%w_bt(i,j,k)*Grd%dte(i,j)*Grd%dus(i,j)     + &
                        Adv_vel%w_bt(i+1,j,k)*Grd%dtw(i+1,j)*Grd%dus(i,j) + &
                        Adv_vel%w_bt(i,j+1,k)*Grd%dte(i,j+1)*Grd%dun(i,j) + &
                        Adv_vel%w_bt(i+1,j+1,k)*Grd%dtw(i+1,j+1)*Grd%dun(i,j))*Grd%daur(i,j)
                enddo
             enddo
          else
             do j = Obc%bound(m)%jsd, Obc%bound(m)%jed
                Adv_vel%w_bu(i,j,:) = 0. 
             enddo
          endif
       end select

       if( .not. Obc%bound(m)%obc_vert_advel_t ) then
          do j = Obc%bound(m)%jsd, Obc%bound(m)%jed
             do i = Obc%bound(m)%isd, Obc%bound(m)%ied
                Adv_vel%w_bt(i,j,:) = 0. 
             enddo
          enddo
       endif

       if(debug_obc) then
          write(stdout(),*) '==> ocean_obc_adjust_advel: '
          write(stdout(),*) '    The vertical bottom velocity has been removed for boundary',m,'!'
       endif
    enddo

    call mpp_clock_end(id_obc)

    return
  end subroutine ocean_obc_adjust_advel


  !</SUBROUTINE>

  !#####################################################################
  !<SUBROUTINE NAME="ocean_obc_adjust_forcing_fs">
  !   <DESCRIPTION>
  !      Add wrong pressure gradient
  !   </DESCRIPTION>
  !   <INOUT NAME="Ext_mode" TYPE="ocean_external_mode_type" > </INOUT>
  !<PUBLICROUTINE>  
  subroutine ocean_obc_adjust_forcing_fs(Ext_mode)
    !</PUBLICROUTINE>
    type(ocean_external_mode_type), intent(inout) :: Ext_mode

    integer :: i, j, m
 
    call mpp_clock_begin(id_obc)

    do m = 1, nobc

       if(.not. Obc%bound(m)%on_bound) cycle
       if(.not. Obc%bound(m)%obc_adjust_forcing_fs) cycle
       
       select case( Obc%bound(m)%direction )
       case(WEST)
         i = Obc%bound(m)%is
         do j = Obc%bound(m)%js, Obc%bound(m)%je
           Ext_mode%forcing_fs(i,j,1) = Ext_mode%forcing_fs(i,j,1) + Obc%bound(m)%pgrad(j) 
         enddo
       case(EAST)
         i = Obc%bound(m)%is-1
         do j = Obc%bound(m)%js, Obc%bound(m)%je
           Ext_mode%forcing_fs(i,j,1) = Ext_mode%forcing_fs(i,j,1) + Obc%bound(m)%pgrad(j) 
         enddo
       case(SOUTH)
         j = Obc%bound(m)%js
         do i = Obc%bound(m)%is, Obc%bound(m)%ie
           Ext_mode%forcing_fs(i,j,2) = Ext_mode%forcing_fs(i,j,2) + Obc%bound(m)%pgrad(i) 
         enddo
       case(NORTH)
         j = Obc%bound(m)%js-1
         do i = Obc%bound(m)%is, Obc%bound(m)%ie
           Ext_mode%forcing_fs(i,j,2) = Ext_mode%forcing_fs(i,j,2) + Obc%bound(m)%pgrad(i) 
         enddo
       end select

       
       if(debug_obc) then
          write(stdout(),*) '==> ocean_obc_adjust_forcing_fs: '
          write(stdout(),*) '    The baroclinic pressure gradient across boundary',m
          write(stdout(),*) '    has been removed from the free surface mode!'
       endif
    enddo      

    call mpp_clock_end(id_obc)

    return
  end subroutine ocean_obc_adjust_forcing_fs
  !</SUBROUTINE>

  !#####################################################################
  !<SUBROUTINE NAME="ocean_obc_tracer">
  !   <DESCRIPTION>
  !      Extrapolate tracer on the open boundaries for ocean model and regional atmosphere model.
  !   </DESCRIPTION>
  !   <IN NAME="wrk2" TYPE="real, dimension(isc:,jsc:,:)"> 
  !     contains rho0r/dht from update_tracer
  !   </IN>
  !   <IN NAME="taum1, tau, taup1" TYPE="integer" >
  !     time step index
  !   </IN>
  !   <IN NAME="time" TYPE="type(time_type)">
  !     model time
  !   </IN>
  !   <IN NAME="name" TYPE="character(len=*)">
  !     tracer name.
  !   </IN>
  !   <IN NAME="n" TYPE="integer">
  !     tracer number
  !   </IN>
  !   <INOUT NAME="tracer" TYPE="real, dimension(isd:,jsd:,:,:)"> 
  !      Tracer field
  !   </INOUT>
  !<PUBLICROUTINE>
  subroutine ocean_obc_tracer(tracer, adv_vet, adv_vnt, wrk2, taum1, tau, taup1, time, name, n)
    !</PUBLICROUTINE>
    real, dimension(isd:,jsd:,:,:), intent(inout) :: tracer            !tracer
    real, dimension(isd:,jsd:,:),      intent(in) :: adv_vet           !advection velocity on east face of t-cell
    real, dimension(isd:,jsd:,:),      intent(in) :: adv_vnt           !advection velocity on north face of t-cell
    real, dimension(isc:,jsc:,:),      intent(in) :: wrk2              !rho0r/dht
    integer,                           intent(in) :: taum1, tau, taup1 ! time step index
    type(time_type),                   intent(in) :: time              ! model time
    character(len=*),                  intent(in) :: name              ! name of the tracer
    integer, intent(in)                           :: n                 ! only when n=1, the max phase speed will be calculated

    !--- local variables -----------------------------------------------
    integer :: it, iu, jt, ju, m

    call mpp_clock_begin(id_obc)

    do m = 1, nobc

       !--- if on current pe there is no point on the bound, then just return
       if(.not.Obc%bound(m)%on_bound) cycle

       !--- compute the barotropic phase speed for the Orlanki condition -- ?? should be taum1 or taup1.

       if(Obc%bound(m)%direction .eq. WEST .or. Obc%bound(m)%direction .eq. EAST ) then
          it = Obc%bound(m)%is
          if(Obc%bound(m)%direction .eq. WEST ) then
            iu = it
          else
            iu = it - 1
          endif
          call ocean_obc_tracer_zonal(Obc%bound(m),Obc%tracer(m,n), tracer, adv_vet(iu,:,:), wrk2(it,:,:), &
               taum1, tau, taup1, n)
       else
          jt = Obc%bound(m)%js
          if(Obc%bound(m)%direction .eq. SOUTH ) then
            ju = jt
          else
            ju = jt - 1
          endif
          call ocean_obc_tracer_merid(Obc%bound(m),Obc%tracer(m,n), tracer, adv_vnt(:,ju,:), wrk2(:,jt,:), &
               taum1, tau, taup1, n)
       endif
    enddo
 
   !--- update the tracer at global halo point to make the gradient accross boundary is 0
!    call ocean_obc_update_boundary(tracer(:,:,:,taup1), 'T')

    if(debug_obc) then
       write(stdout(),*) 'After ocean_obc_tracer, tracer chksum ===>', mpp_chksum(tracer(isc:iec,jsc:jec,:,taup1))
    endif

    call mpp_clock_end(id_obc)

    return
  end subroutine ocean_obc_tracer
  !</SUBROUTINE>

  !#######################################################################
  subroutine ocean_obc_tracer_zonal(Bound, Data, tracer, adv_vet, wrk2, taum1, tau, taup1, n)

    type(obc_bound_type),              intent(in) :: Bound             ! lateral boundary data
    type(obc_data_type_3d) ,           intent(in) :: Data              ! lateral boundary data
    real, dimension(isd:,jsd:,:,:), intent(inout) :: tracer            ! tracer
    real, dimension(jsd:,:),           intent(in) :: adv_vet           ! advection velocity
    real, dimension(jsc:,:),           intent(in) :: wrk2              
    integer,                           intent(in) :: taum1, tau, taup1 !time step index
    integer,                           intent(in) :: n                 ! tracer index

    !--- local variables -----------------------------------------------
    integer :: k, i, j, sign, tlevel
    real    :: cgrid, var, var1, var2, cmax, uout, uin
    real, dimension(jsc:jec,nk) :: cclin, adv_tx_obc

    i = Bound%is
    tlevel = taum1
    !----------------------------------------------------------------------
    !--- compute the advective velocity "uad" ----------------------------- 
    !--- at the east or west face of the T cells and ----------------------
    !--- compute phase velocity at the western or eastern boundary: cclin -
    if(debug_obc) then
      write(stdout(),*) 'Doing update of tracer ',n,' at obc ', Bound%direction
      write(stdout(),*) ' i = ', i
    endif
    
    sign = 1
    if(Bound%direction == EAST) sign = -1
    do j = Bound%js, Bound%je
       cgrid = Grd%dxu(i+sign,j)/dtime_t
       if(Bound%direction == WEST)  then
         cmax  = max(Bound%ctrop(j),-cgrid)
       else
         cmax  = min(Bound%ctrop(j),cgrid)
       endif
       do k = 1, nk
          var2 = sign*( tracer(i+2*sign,j,k,tau)-tracer(i+sign,j,k,tau) )
          if(abs(var2) .lt. small) then
             cclin(j,k) = cmax
          else
             cclin(j,k) = -(tracer(i+sign,j,k,taup1)-tracer(i+sign,j,k,taum1))*cgrid/var2
             if(cclin(j,k)*sign .gt. 0.0) then
               if( Bound%obc_relax_tracer ) then
                  cclin(j,k) = 0.
               else
                  cclin(j,k) = 0.5*cmax
               endif
             endif
          endif
          if(Bound%direction == WEST)  then
             cclin(j,k) = max(cmax,cclin(j,k))
          else
             cclin(j,k) = min(cmax,cclin(j,k))
          endif
       enddo
    enddo
    if(debug_obc) then
      k = 1
      if(Bound%direction == WEST) then
        j = jec
      else
        j = jsc + 1
      endif
      write(stdout(),*) 'k= 1, j = ',j
      write(stdout(),*) 'cgrid   = ',cgrid
      write(stdout(),*) 'cmax    = ',cmax
      write(stdout(),*) 'ctrop   = ',Bound%ctrop(j)
      write(stdout(),*) 'clin    = ',cclin(j,k)
    endif

    !----------------------------------------------------------------
    !--- recalculate lateral advection of tracers
    !--- radiation condition at the western or eastern wall
    !----------------------------------------------------------------

    ! orlanski   : upwind advection for outflow, no advection for inflow
    ! with relax : upwind advection in any case, 
    !              for inflow with prescribed tracer
    do j= Bound%js, Bound%je
       do k= 1, nk
          uout = 0.5*(adv_vet(j,k) - sign * abs(adv_vet(j,k)))
          if(Bound%obc_tracer_orlanski) then
             adv_tx_obc(j,k) = - sign * Grd%tmask(i,j,k) * Grd%datr(i,j) &
                               * uout* ( tracer(i,j,k,tlevel)*Grd%dyte(i,j)-tracer(i+sign,j,k,tlevel)*Grd%dyte(i+sign,j) )
          else
             uin  = 0.5*(adv_vet(j,k) + sign * abs(adv_vet(j,k)))
             adv_tx_obc(j,k) = - sign * Grd%tmask(i,j,k) * Grd%datr(i,j)&
                               *(                                       &
                                  uin * ( Data%data(1,j,k)*Grd%dyte(i-sign,j)-tracer(i,j,k,tlevel)*Grd%dyte(i,j) ) &
                                + uout* ( tracer(i,j,k,tlevel)*Grd%dyte(i,j)-tracer(i+sign,j,k,tlevel)*Grd%dyte(i+sign,j) )&
                                 )
          endif
       enddo
    enddo

    if(debug_obc) then
      write(stdout(),*) 'tlevel, tm1, t, tp1 = ',tlevel, taum1, tau, taup1
      if(Bound%direction == WEST) then
        k = 1
        j = jec
        write(stdout(),*) 'k=1,  j   = ',j
        write(stdout(),*) 't(i,tm1)  = ',tracer(isd:i+2,j,k,taum1)
        write(stdout(),*) 't(i,t)    = ',tracer(isd:i+2,j,k,tau)
        write(stdout(),*) 't(i,tp1)  = ',tracer(isd:i+2,j,k,taup1)
        write(stdout(),*) 'var2      = ',sign * dtime_t / Grd%dxu(i,j) * cclin(j,k)
      else
        k = 1
        j = jsc + 1
        write(stdout(),*) 'k=1,  j   = ',j
        write(stdout(),*) 't(i,tm1)  = ',tracer(iec-3:i+1,j,k,taum1)
        write(stdout(),*) 't(i,t)    = ',tracer(iec-3:i+1,j,k,tau)
        write(stdout(),*) 't(i,tp1)  = ',tracer(iec-3:i+1,j,k,taup1)
        write(stdout(),*) 'var2      = ',sign * dtime_t / Grd%dxu(i,j) * cclin(j,k)
      endif
      write(stdout(),*) 'hadv_old  = ',Bound%hadv(j,k)
      write(stdout(),*) 'hadv_new  = ',adv_tx_obc(j,k)
      write(stdout(),*) 'wrk2      = ',wrk2(j,k)
    endif
    !----------------------------------------------------------------
    !--- calculate the new tracer quantities allowing for implicit
    !--- treatment of vertical diffusion
    !----------------------------------------------------------------
    var1 = dtime_t*Data%rel_coef
    do j = Bound%js, Bound%je
       var = dtime_t / Grd%dxu(i,j)
       do k = 1, nk
          tracer(i,j,k,taup1) = tracer(i,j,k,taup1) + dtime_t * ( &
                          Bound%hadv(j,k)  & ! remove original advection
                        - adv_tx_obc(j,k)) * Grd%tmask(i,j,k) &
                        * rho0 * wrk2(j,k)
          !---  Add the wave-like contribution implicitly ----------------
          var2 = sign * var * cclin(j,k)
          tracer(i,j,k,taup1) = (tracer(i,j,k,taup1) - tracer(i+sign,j,k,taup1)*var2)/(1. - var2)
          if( Bound%obc_relax_tracer ) then
             !---  Add the relaxation implicitly -------------------------
             tracer(i,j,k,taup1) = (tracer(i,j,k,taup1)+ var1*Data%data(1,j,k))/(1.0+var1)
          endif
       enddo
    enddo
    if(debug_obc) then
      write(stdout(),*) ' after updating open boundary'
      if(Bound%direction == WEST) then
        k = 1
        j = jec
        write(stdout(),*) 'k=1,  j   = ',j
        write(stdout(),*) 't(i,tp1)  = ',tracer(isd:i+2,j,k,taup1)
      else
        k = 1
        j = jsc + 1
        write(stdout(),*) 'k=1,  j   = ',j
        write(stdout(),*) 't(i,tp1)  = ',tracer(iec-3:i+1,j,k,taup1)
      endif
    endif
    return
  end subroutine ocean_obc_tracer_zonal

  !#######################################################################
  subroutine ocean_obc_tracer_merid(Bound, Data, tracer, adv_vnt, wrk2, taum1, tau, taup1, n)

    type(obc_bound_type),              intent(in) :: Bound             ! lateral boundary data
    type(obc_data_type_3d) ,           intent(in) :: Data              ! lateral boundary data
    real, dimension(isd:,jsd:,:,:), intent(inout) :: tracer            ! tracer
    real, dimension(isd:,:),           intent(in) :: adv_vnt           ! advection velocity
    real, dimension(isc:,:),           intent(in) :: wrk2              
    integer,                           intent(in) :: taum1, tau, taup1 ! time step index
    integer,                           intent(in) :: n                 ! tracer index

    !--- local variables -----------------------------------------------
    integer :: k, i, j, sign, tlevel
    real    :: cgrid, var, var1, var2, cmax, uout, uin
    real, dimension(isc:iec,nk) :: cclin, adv_ty_obc

    j = Bound%js
    !------------------------------------------------------------------------
    !--- compute the advective velocity "uad" -------------------------------
    !--- at the south or north face of the T cells and ----------------------
    !--- compute phase velocity at the southern or northern boundary: cclin -
    sign = 1
    if(Bound%direction == NORTH) sign = -1

    do i = Bound%is, Bound%ie
       cgrid = Grd%dyu(i,j+sign)/dtime_t
       if(Bound%direction == NORTH)  then
         cmax  = max(Bound%ctrop(i),-cgrid)
       else
         cmax  = min(Bound%ctrop(i),cgrid)
       endif
       do k = 1, nk
          var2 = sign*( tracer(i,j+2*sign,k,taum1)-tracer(i,j+sign,k,taum1) )
          if(abs(var2) .lt. small) then
             cclin(i,k) = cmax
          else
             cclin(i,k) = -(tracer(i,j+sign,k,taup1)-tracer(i,j+sign,k,taum1))*cgrid/var2
             if(cclin(i,k)*sign .gt. 0.0) then
               if( Bound%obc_relax_tracer ) then
                  cclin(i,k) = 0.
               else
                  cclin(i,k) = 0.5*cmax
               endif
             endif
          endif
         if(Bound%direction == SOUTH) then
             cclin(i,k) = max(cmax,cclin(i,k))
         else
             cclin(i,k) = min(cmax,cclin(i,k))
         endif
       enddo
    enddo

    !----------------------------------------------------------------
    !--- recalculate lateral advection of tracers
    !--- radiation condition at the western or eastern wall
    !----------------------------------------------------------------

    ! orlanski   : upwind advection for outflow, no advection for inflow
    ! with relax : upwind advection in any case, 
    !              for inflow with prescribed tracer
    tlevel = taum1
    do j= Bound%is, Bound%ie
       do k= 1, nk
          uout = 0.5*(adv_vnt(i,k) - sign * abs(adv_vnt(i,k)))
          if(Bound%obc_tracer_orlanski) then
             adv_ty_obc(i,k) = - sign * Grd%tmask(i,j,k) * Grd%datr(i,j) &
                               * uout* ( tracer(i,j,k,tlevel)*Grd%dxtn(i,j)-tracer(i,j+sign,k,tlevel)*Grd%dxtn(i,j+sign) )
          else
             uin  = 0.5*(adv_vnt(i,k) + sign * abs(adv_vnt(i,k)))
             adv_ty_obc(i,k) = - sign * Grd%tmask(i,j,k) * Grd%datr(i,j)&
                               *(                                       &
                                  uin * ( Data%data(1,j,k)*Grd%dytn(i,j-sign)-tracer(i,j,k,tlevel)*Grd%dxtn(i,j) ) &
                                + uout* ( tracer(i,j,k,tlevel)*Grd%dxtn(i,j)-tracer(i,j+sign,k,tlevel)*Grd%dxtn(i,j+sign) )&
                                 )
          endif
       enddo
    enddo

    !----------------------------------------------------------------
    !--- calculate the new tracer quantities allowing for implicit
    !--- treatment of vertical diffusion
    !----------------------------------------------------------------
    var1 = dtime_t*Data%rel_coef
    do i = Bound%is, Bound%ie
       var = dtime_t / Grd%dyu(i,j)
       do k = 1, nk
          tracer(i,j,k,taup1) = tracer(i,j,k,taup1) + dtime_t * ( &
                          Bound%hadv(i,k)  & ! remove original advection
                        - adv_ty_obc(i,k)) * Grd%tmask(i,j,k) &
                        * rho0 * wrk2(i,k)
          !---  Add the wave-like contribution implicitely ----------------
          var2 = sign * var * cclin(i,k)
          tracer(i,j,k,taup1) = (tracer(i,j,k,taup1) - tracer(i,j+sign,k,taup1)*var2)/(1. - var2)
          if(Bound%obc_relax_tracer ) then
             !---  Add the relaxation implicitely -------------------------
             tracer(i,j,k,taup1) = (tracer(i,j,k,taup1)+ var1* Data%data(i,1,k))/(1.+var1)
          endif
       enddo
    enddo

    return

  end subroutine ocean_obc_tracer_merid

  !######################################################################
  !--- update field on the halo points at the boundaries
  !<SUBROUTINE NAME="ocean_obc_update_boundary_2d" INTERFACE="ocean_obc_update_boundary">
  !   <INOUT NAME="field" TYPE="real, dimension(:,:)"></INOUT> 
  !   <IN NAME="grid_type" TYPE=" character(len=1)"></IN> 
  !   <IN NAME="update_type" TYPE=" character(len=1)"></IN> 
  !<DESCRIPTION>
  ! This subroutine sets values at the halo points beyond open boundaries
  ! grid_type can be:    'T' for variables at tracer point or 
  !                      'C' for variables at velocity points 
  !                      'Z' velocity points at zonal (north, south) boundaries
  !                      'M' velocity points at meridional (east, west) boundaries
  ! 
  ! update_type can be: 's' for setting the field beyond the 
  !                         boundary to the value at the boundary 
  !                         - for a velocity 'Z' at northern and southern 
  !                           boundary
  !                         - for a velocity 'M' at western and eastern 
  !                           boundary
  !                         - for gridtype 'T' and 'C' at all boundaries
  !                     'x' for linear extrapolation of the field beyond the 
  !                         boundary from the value at the boundary and the 
  !                         first internal point
  !                     'i' for setting the field beyond the boundary to   
  !                         the value at the first internal point:
  !                         - for a velocity 'Z' at northern and southern 
  !                           boundary
  !                         - for a velocity 'M' at western and eastern 
  !                           boundary
  !                         - for gridtype 'T' and 'C' at all boundaries
  !                        
  !</DESCRIPTION>
  !<PUBLICROUTINE INTERFACE="ocean_obc_update_boundary" >
  ! 

  subroutine ocean_obc_update_boundary_2d(field,grid_type,update_type)
    !</PUBLICROUTINE>

    real, dimension(isd:,jsd:), intent(inout)  :: field
    character(len=1),               intent(in) :: grid_type
    character(len=1), optional,     intent(in) :: update_type

    integer :: i, j, m, istr, jstr, iend, jend, iref, jref, i1, i2, j1, j2
    real    :: u1, u2, x1, x2, v1, v2, y1, y2
    character(len=1) :: uptype          ! set field at points outside the boundary to 
                                        ! the boundary value
    uptype = 's' 
    if(PRESENT(update_type)) uptype=update_type
    
    if(uptype .ne. 's' .and. uptype .ne. 'x'.and. uptype .ne. 'i') call mpp_error(FATAL, &
         'ocean_obc_mod(ocean_obc_update_boundary) : update_type= '// update_type//' should be either "s", "x" or "i" ' )
    
    if(grid_type .ne. 'T' .and. grid_type .ne. 'C'.and. grid_type .ne. 'Z'.and. grid_type .ne. 'M') call mpp_error(FATAL, &
         'ocean_obc_mod(ocean_obc_update_boundary) : grid_type= '// grid_type//' should be either "T", "C", "Z" or "M" ' )
      
    do m = 1, nobc
      !--- if on current pe there is no point on the bound, then just return
      if(.not.Obc%bound(m)%on_bound) cycle
            
      select case( Obc%bound(m)%direction )
      
      case (WEST)  
        if( grid_type == 'Z' ) cycle  
        
        if(debug_obc) then
          write(stdout(),*) 'ocean_obc_update_boundary_2d: updating ',trim(Obc%bound(m)%name), ' ',&
                       Obc%bound(m)%direction, ' grid_type: ',grid_type,' update_type: ',uptype
        endif
        
        istr = max(isd+1,Obc%bound(m)%is-1)
        iend   = Obc%bound(m)%is-1

        if (uptype == 'x') then
          i1 = Obc%bound(m)%is
          i2 = Obc%bound(m)%is+1
          do j = Obc%bound(m)%jsd, Obc%bound(m)%jed
            if(Obc%bound(m)%mask(j)) then
              u1 = field(i1,j)
              u2 = field(i2,j)
              x1 = Grd%xt(i1,j)
              x2 = Grd%xt(i2,j)
              do i = istr, iend
                 field(i,j) = u1 + (u2-u1)*(Grd%xt(i,j) - x1)/(x2-x1)
              enddo
            endif
          enddo
        else
          if (uptype == 's') then
            iref = iend+1
          else
            iref = iend+2
          endif
           do j = Obc%bound(m)%jsd, Obc%bound(m)%jed
            if(Obc%bound(m)%mask(j)) then
              do i = istr, iend
                 field(i,j) = field(iref,j)
              enddo
            endif
          enddo
        endif
        
      case (EAST)
        if( grid_type == 'Z' ) cycle  

        if(debug_obc) then
          write(stdout(),*) 'ocean_obc_update_boundary_2d: updating ',trim(Obc%bound(m)%name), ' ',&
                       Obc%bound(m)%direction, ' grid_type: ',grid_type,' update_type: ',uptype
        endif
                
        iend = max(ied-1,Obc%bound(m)%ie+1)
        
        if (uptype == 'x') then
          if(grid_type == 'T') then
            i1 = Obc%bound(m)%ie
            i2 = Obc%bound(m)%ie-1
            istr = Obc%bound(m)%ie+1    
          else
            i1 = Obc%bound(m)%ie-1
            i2 = Obc%bound(m)%ie-2
            istr = Obc%bound(m)%ie   
          endif
          do j = Obc%bound(m)%jsd, Obc%bound(m)%jed
            if(Obc%bound(m)%mask(j)) then
              u1 = field(i1,j)
              u2 = field(i2,j)
              x1 = Grd%xt(i1,j)
              x2 = Grd%xt(i2,j)
              do i = istr, iend
                field(i,j) = u1 + (u2-u1)*(Grd%xt(i,j) - x1)/(x2-x1) 
              enddo
             endif
          enddo
        else
          istr = Obc%bound(m)%ie    ! for velocity points
        
          if (uptype == 's') then
            iref = istr-1
          else
            iref = istr-2
          endif
          if(grid_type == 'T') then   ! adjust for tracer points
            istr = istr+1
            iref = iref+1
          endif
          do j = Obc%bound(m)%jsd, Obc%bound(m)%jed
            if(Obc%bound(m)%mask(j)) then
              do i = istr, iend
                field(i,j) = field(iref,j)
              enddo
            endif
          enddo
        endif

      case (SOUTH)
        if( grid_type == 'M' ) cycle  

        if(debug_obc) then
          write(stdout(),*) 'ocean_obc_update_boundary_2d: updating ',trim(Obc%bound(m)%name), ' ',&
                       Obc%bound(m)%direction, ' grid_type: ',grid_type,' update_type: ',uptype
        endif
        jstr = max(jsd+1,Obc%bound(m)%js-1)
        jend   = Obc%bound(m)%js-1
        
        if (uptype == 'x') then
          j1 = Obc%bound(m)%js
          j2 = Obc%bound(m)%js+1
          do i = Obc%bound(m)%isd, Obc%bound(m)%ied
             if(Obc%bound(m)%mask(i)) then
                v1 = field(i,j1)
                v2 = field(i,j2)
                y1 = Grd%xt(i,j1)
                y2 = Grd%xt(i,j2)
                do j = jstr, jend
                   field(i,j) = v1 + (v2-v1)*(Grd%yt(i,j) - y1)/(y2-y1)
                enddo
             endif
          enddo
        else
          if (uptype == 's') then
            jref = jend+1
          else
            jref = jend+2
          endif
          do i = Obc%bound(m)%isd, Obc%bound(m)%ied
            if(Obc%bound(m)%mask(i)) then
              do j = jstr, jend
                field(i,j) = field(i,jref) 
              enddo
            endif
          enddo
        endif
        
      case (NORTH)
        if( grid_type == 'M' ) cycle  

        if(debug_obc) then
          write(stdout(),*) 'ocean_obc_update_boundary_2d: updating ',trim(Obc%bound(m)%name), ' ',&
                       Obc%bound(m)%direction, ' grid_type: ',grid_type,' update_type: ',uptype
        endif
        jend = max(Obc%bound(m)%je+1,jed-1)
        if (uptype == 'x') then
          if(grid_type == 'T') then
            j1 = Obc%bound(m)%je
            j2 = Obc%bound(m)%je-1
            jstr = Obc%bound(m)%je+1
          else
            j1 = Obc%bound(m)%je-1
            j2 = Obc%bound(m)%je-2
            jstr = Obc%bound(m)%je
          endif
          do i = Obc%bound(m)%isd, Obc%bound(m)%ied
            if(Obc%bound(m)%mask(i)) then
              v1 = field(i,j1)
              v2 = field(i,j2)
              y1 = Grd%xt(i,j1)
              y2 = Grd%xt(i,j2)
              do j = jstr, jend
                field(i,j) = v1 + (v2-v1)*(Grd%yt(i,j) - y1)/(y2-y1)
              enddo
            endif
          enddo
        else
          jstr = Obc%bound(m)%je
          if (uptype == 's') then
            jref = jstr-1
          else
            jref = jstr-2
          endif
          if(grid_type == 'T') then   ! adjust for tracer points
            jstr = jstr+1
            jref = jref+1
          endif
          do i = Obc%bound(m)%isd, Obc%bound(m)%ied
            if(Obc%bound(m)%mask(i)) then
              do j = jstr, jend
                field(i,j) = field(i,jref) 
              enddo
            endif
          enddo
        endif
        
      end select

    enddo

    return

  end subroutine ocean_obc_update_boundary_2d
  !</SUBROUTINE>

  !#####################################################################
  !--- update field on the halo points at the boundaries 
  !<SUBROUTINE NAME="ocean_obc_update_boundary_3d" INTERFACE="ocean_obc_update_boundary">
  !   <INOUT NAME="field" TYPE="real, dimension(:,:,:)"></INOUT> 
  subroutine ocean_obc_update_boundary_3d(field,grid_type,update_type)
    real, dimension(isd:,jsd:,:), intent(inout) :: field
    character(len=1),                intent(in) :: grid_type
    character(len=1), optional,      intent(in) :: update_type
    integer :: k, n3

    n3 = size(field,3)
    if(PRESENT(update_type)) then
      do k = 1, n3
        call ocean_obc_update_boundary_2d(field(:,:,k),grid_type,update_type)
      enddo
    else
      do k = 1, n3
        call ocean_obc_update_boundary_2d(field(:,:,k),grid_type)
      enddo
    endif
    return
  end subroutine ocean_obc_update_boundary_3d
  !</SUBROUTINE>

  !#####################################################################
  !--- update field on the halo points at the boundaries 
  !<SUBROUTINE NAME="ocean_obc_update_boundary_4d" INTERFACE="obc_update_boundary">
  !   <INOUT NAME="field" TYPE="real, dimension(:,:,:,:)"></INOUT> 
  subroutine ocean_obc_update_boundary_4d(field,grid_type,update_type)
    real, dimension(isd:,jsd:,:,:), intent(inout) :: field
    character(len=1),                  intent(in) :: grid_type
    character(len=1), optional,        intent(in) :: update_type
    integer :: n3, n4, k, n

    call mpp_clock_begin(id_obc)
    
    n3 = size(field,3)
    n4 = size(field,4)
    
    if(PRESENT(update_type)) then
      do k = 1, n3
        do n = 1, n4
          call ocean_obc_update_boundary_2d(field(:,:,k,n),grid_type,update_type)
        enddo
      enddo
    else
      do k = 1, n3
        do n = 1, n4
          call ocean_obc_update_boundary_2d(field(:,:,k,n),grid_type)
        enddo
      enddo
    endif

    call mpp_clock_end(id_obc)

    return
  end subroutine ocean_obc_update_boundary_4d
  !</SUBROUTINE>

  !#####################################################################
  subroutine ocean_obc_adjust_topog(ht, hu, kmt, kmu)

    real, dimension(isd:,jsd:),    intent(inout) :: ht
    real, dimension(isd:,jsd:),    intent(inout) :: hu
    integer, dimension(isd:,jsd:), intent(inout) :: kmt
    integer, dimension(isd:,jsd:), intent(inout) :: kmu

    integer :: i, ib, istr, iend, j, jb, jstr, jend, m, kmt_temp
    real    :: ht_temp

    do m = 1, nobc
             
       if(.not. Obc%bound(m)%on_bound) cycle

       if(debug_obc) write(stdout(),*) 'adjust topog ', Obc%bound(m)%direction
       
       select case( Obc%bound(m)%direction )
       case (WEST) 
          call mpp_error(NOTE,'ocean_obc_mod: topography is modified on the west open boundary') 
          ib   = Obc%bound(m)%is
          istr = max(isd+1, ib-1)
          do j = jsd, jed 
             if(Obc%bound(m)%mask(j)) then
                kmt_temp = min(kmt(ib,j),kmt(ib+1,j),kmt(ib+2,j))
                if(kmt_temp.eq.0) then 
                   ht_temp = 0.0
                else
                   ht_temp = Grd%zw(kmt_temp)
                endif
                do i = istr, ib+2
                   kmt(i,j)=kmt_temp
                   ht(i,j)=ht_temp
                enddo
             endif
          enddo
          do j = jsc, jec
             if(Obc%bound(m)%mask(j)) then
                do i=istr, ib+2
                   kmu(i,j) = min(kmt(i,j), kmt(i+1,j), kmt(i,j+1), kmt(i+1,j+1))
                   hu(i,j)  = min(ht(i,j),  ht(i+1,j),  ht(i,j+1),  ht(i+1,j+1))
                enddo
             endif
          enddo
       case (EAST) 
          call mpp_error(NOTE,'ocean_obc_mod: topography is modified on the east open boundary') 
          ib   = Obc%bound(m)%ie
          iend = min(ied-1, ib+1)
          do j = jsd, jed
             if(Obc%bound(m)%mask(j)) then
                kmt_temp = min(kmt(ib,j),kmt(ib-1,j),kmt(ib-2,j))
                if(kmt_temp.eq.0) then 
                   ht_temp = 0.0
                else
                   ht_temp = Grd%zw(kmt_temp)
                endif
                do i = ib-2,iend
                   kmt(i,j)=kmt_temp
                   ht(i,j)=ht_temp
                enddo
             endif
          enddo
          do j = jsc, jec
             if(Obc%bound(m)%mask(j)) then
                do i=ib-3,iend-1
                   kmu(i,j) = min(kmt(i,j), kmt(i+1,j), kmt(i,j+1), kmt(i+1,j+1))
                   hu(i,j)  = min(ht(i,j),  ht(i+1,j),  ht(i,j+1),  ht(i+1,j+1))
                enddo
                kmu(iend,j) = kmu(iend-1,j)
                hu(iend,j)  = hu(iend-1,j)
             endif
          enddo
       case (SOUTH) 
          call mpp_error(NOTE,'ocean_obc_mod: topography is modified on the south open boundary') 
          jb   = Obc%bound(m)%js
          jstr = max(jsd, jb-2)
          do i = isd, ied
             if(Obc%bound(m)%mask(i)) then
                kmt_temp = min(kmt(i,jb),kmt(i,jb+1),kmt(i,jb+2))
                if(kmt_temp.eq.0) then 
                   ht_temp = 0.0
                else
                   ht_temp = Grd%zw(kmt_temp)
                endif
                do j = jb, jb+2
                   kmt(i,j)=kmt_temp
                   ht(i,j)=ht_temp
                enddo
             endif
          enddo
          do i = isd, ied
             if(Obc%bound(m)%mask(i)) then
                do j = jb, jb+2
                   kmu(i,j) = min(kmt(i,j), kmt(i+1,j), kmt(i,j+1), kmt(i+1,j+1))
                   hu(i,j)  = min(ht(i,j),  ht(i+1,j),  ht(i,j+1),  ht(i+1,j+1))
                enddo
             endif
          enddo
       case (NORTH) 
          jb   = Obc%bound(m)%js
          jend = min(jed, jb+2)
          call mpp_error(NOTE,'ocean_obc_mod: topography is modified on the north open boundary') 
          do i = isd, ied
             if(Obc%bound(m)%mask(i)) then
                kmt_temp = min(kmt(i,jb),kmt(i,jb-1),kmt(i,jb-2))
                if(kmt_temp.eq.0) then 
                   ht_temp = 0.0
                else
                   ht_temp = Grd%zw(kmt_temp)
                endif
                do j = jb-2, jend
                   kmt(i,j)=kmt_temp
                   ht(i,j)=ht_temp
                enddo
             endif
          enddo
          do i = isd, ied
             if(Obc%bound(m)%mask(i)) then
                do j = jb-3, jend-1
                   kmu(i,j) = min(kmt(i,j), kmt(i+1,j), kmt(i,j+1), kmt(i+1,j+1))
                   hu(i,j)  = min(ht(i,j),  ht(i+1,j),  ht(i,j+1),  ht(i+1,j+1))
                enddo
                kmu(i,jend) = kmu(i,jend-1)
                hu(i,jend)  = hu(i,jend-1)
             endif
          enddo
       end select
    enddo

    call mpp_update_domains(kmt(:,:),Dom%domain2d)
    call mpp_update_domains(ht(:,:),Dom%domain2d)
    call mpp_update_domains(kmu(:,:),Dom%domain2d)
    call mpp_update_domains(hu(:,:),Dom%domain2d)

  end subroutine ocean_obc_adjust_topog


  
  subroutine ocean_obc_set_mask
    integer :: i, j, m

    do j=jsd,jed
      do i=isd,ied
         Grd%obc_tmask(i,j) = min(1.0, float(Grd%kmt(i,j)))
         Grd%obc_umask(i,j) = min(1.0, float(Grd%kmu(i,j)))
      enddo
    enddo

!   Zero out OBC points in Grd%obc_xmask to exclude OBC points from diagnostics    

    do m = 1, nobc
      if(.not. Obc%bound(m)%on_bound) cycle
      
      select case( Obc%bound(m)%direction )
      case(WEST)
        do i = isd, Obc%bound(m)%is
          do j = Obc%bound(m)%js, Obc%bound(m)%je
            Grd%obc_tmask(i,j) = 0.
            Grd%obc_umask(i,j) = 0.
          enddo
        enddo
      case(EAST)
        do i = Obc%bound(m)%is, ied
          do j = Obc%bound(m)%js, Obc%bound(m)%je
            Grd%obc_tmask(i,j)   = 0.
            Grd%obc_umask(i-1,j) = 0.
          enddo
        enddo
      case(SOUTH)
        do j = jsd, Obc%bound(m)%js
          do i = Obc%bound(m)%is, Obc%bound(m)%ie
            Grd%obc_tmask(i,j) = 0.
            Grd%obc_umask(i,j) = 0.
          enddo
        enddo
      case(NORTH)
        do j = Obc%bound(m)%js, jed
          do i = Obc%bound(m)%is, Obc%bound(m)%ie
            Grd%obc_tmask(i,j)   = 0.
            Grd%obc_umask(i,j-1) = 0.
          enddo
        enddo
      end select
    enddo

    return
    
  end subroutine ocean_obc_set_mask

  !#####################################################################
  !--- release memory --------------------------------------------------
  !<SUBROUTINE NAME="ocean_obc_end" >
  ! <DESCRIPTION>
  !    Destructor routine. Release memory.
  !   </DESCRIPTION>
  !   <OUT NAME="have_obc" TYPE="logical">
  !      Contains open boundary information
  !   </OUT>

  subroutine ocean_obc_end(have_obc)

    logical, intent(inout)      :: have_obc
    deallocate( Obc%bound )

    module_is_initialized = .FALSE.
    have_obc = .FALSE.
    
    return

  end subroutine ocean_obc_end
  !</SUBROUTINE>  

  !#######################################################################
  subroutine barotropic_phase_speed(eta, taum1, tau, taup1, tstep )

    real, dimension(isd:,jsd:,:), intent(in) :: eta
    integer,                      intent(in) :: taum1, tau, taup1
    real, intent(in)                         :: tstep

    real    :: var, detai, cmax, cmin, c1tmp
    integer :: i, j, sign, m

    do m = 1, nobc
       !--- if on current pe there is no point on the bound, then just return
       if(.not. Obc%bound(m)%on_bound) cycle

       sign = 1   ! on south or west direction, sign = 1

       if(Obc%Bound(m)%direction == WEST .or. Obc%Bound(m)%direction == EAST) then
          !--- calculate the phase speed at west or east boundary -------------------
          if(Obc%Bound(m)%direction == EAST)  sign = -1
          i = Obc%Bound(m)%is

          do j= Obc%Bound(m)%js, Obc%Bound(m)%je
             var = Grd%dxu(i+sign,j)/tstep
             detai =( eta(i+2*sign,j,tau) - eta(i+sign,j,tau) ) * sign
             !         clip with shallow water phase speed and CFL-phase speed
             cmax  = - sign * min(Obc%bound(m)%ctrop_max*sqrt(grav*Grd%ht(i,j)), 0.5 * var)
             if(Obc%bound(m)%obc_relax_eta) then
               cmin = 0.
             else
               cmin  = - sign * Obc%bound(m)%ctrop_min * sqrt(grav*Grd%ht(i,j))
             endif

             if (ABS(detai).lt. small) then
                !            badly defined phases -> do (about) nothing
                c1tmp = 0.99*Obc%Bound(m)%ctrop(j)
             else
                !            implicit time scheme for internal points
                c1tmp = - var*(eta(i+sign,j,taup1)-eta(i+sign,j,taum1))/detai
             endif
             !         If incoming waves are detected:
             !         - if relaxation is done, Obc%Bound(m)%ctrop(j) is about zero    
             !           -> boundary follows relaxation and generates incoming waves in the internal area
             !         - if no relaxation is done, Obc%Bound(m)%ctrop(j) is about cmin 
             !           -> replace undesired incoming waves by outgoing waves
             if (c1tmp*sign .ge. 0.0) c1tmp = 0.5*Obc%Bound(m)%ctrop(j)   
             if(sign == 1) then
                c1tmp = max(cmax,c1tmp)                    ! cmax is the phase speed limit
                c1tmp = min(cmin,c1tmp)                    ! ensure minimum phase speed
             else
                c1tmp = min(cmax,c1tmp)                    ! cmax is the phase speed limit
                c1tmp = max(cmin,c1tmp)                    ! ensure minimum phase speed
             endif
             Obc%Bound(m)%ctrop(j) = 0.7*Obc%Bound(m)%ctrop(j) + 0.3*c1tmp    ! mix with previous phase values
          enddo

       else     ! south or north direction
          if(Obc%Bound(m)%direction == NORTH) sign = -1
          j = Obc%Bound(m)%js

          do i=Obc%Bound(m)%is, Obc%Bound(m)%ie
             var = Grd%dyu(i,j+sign)/tstep
             detai =( eta(i,j+2*sign,tau) - eta(i,j+sign,tau) ) * sign
             !         clip with shallow water phase speed and CFL-phase speed
             cmax  = - sign * min(Obc%bound(m)%ctrop_max*sqrt(grav*Grd%ht(i,j)), 0.5 * var)
             if(Obc%bound(m)%obc_relax_eta) then
               cmin = 0.
             else
               cmin  = - sign * Obc%bound(m)%ctrop_min * sqrt(grav*Grd%ht(i,j))
             endif
             if (ABS(detai).lt. small) then
                !            badly defined phases -> do (about) nothing
                c1tmp = 0.99*Obc%Bound(m)%ctrop(i)
             else
                !            implicit time scheme for internal points
                c1tmp = - var*(eta(i,j+sign,taup1)-eta(i,j+sign,taum1))/detai
             endif
             !         If incoming waves are detected:
             !         - if relaxation is done, Obc%Bound(m)%ctrop(j) is about zero    
             !           -> boundary follows relaxation and generates incoming waves in the internal area
             !         - if no relaxation is done, Obc%Bound(m)%ctrop(j) is about cmin 
             !           -> replace undesired incoming waves by outgoing waves
             if (c1tmp*sign .ge. 0.0) c1tmp = 0.5*Obc%Bound(m)%ctrop(i) 
             if(sign == 1) then 
                c1tmp = max(cmax,c1tmp)                    ! cmax is the phase speed limit
                c1tmp = min(cmin,c1tmp)                    ! ensure minimum phase speed
             else
                c1tmp = min(cmax,c1tmp)                    ! cmax is the phase speed limit
                c1tmp = max(cmin,c1tmp)                    ! ensure minimum phase speed
             endif
             Obc%Bound(m)%ctrop(i) = 0.7*Obc%Bound(m)%ctrop(i) + 0.3*c1tmp    ! mix with previous phase values
          enddo
       endif
    enddo

    return
  end subroutine barotropic_phase_speed

  !#####################################################################
  subroutine barotropic_phase_speed_init( eta, taum1, tau )

!   initialize the barotropic phase speed in the first call of the barotropic
!   OBC. Usually, eta is Ext_mode%eta_t or Ext_mode%eta_t_bar with splitting on  

    real, dimension(isd:,jsd:,:), intent(in) :: eta
    integer,                      intent(in) :: taum1, tau

    real    :: var, detai, cmax, cmin, c1tmp, cgrid, csur 
    integer :: i, j, m, sign


    do m = 1, nobc

       !--- if on current pe there is no point on the bound, then just return
       if(.not. Obc%bound(m)%on_bound) cycle

       sign = 1   ! on south or west direction, sign = 1
       if(Obc%bound(m)%direction == WEST .or. Obc%bound(m)%direction == EAST) then
          !--- calculate the phase speed at west or east boundary -------------------
          if(Obc%bound(m)%direction == EAST)  sign = -1
          
          i = Obc%bound(m)%is

          do j= Obc%bound(m)%js, Obc%bound(m)%je
!            Use the barotropic timestep for the CFL-clipping,          
             var = Grd%dxu(i+sign,j)/dtfs 
             cgrid = Grd%dxu(i+sign,j)/dteta
             detai =( eta(i+2*sign,j,tau) - eta(i+sign,j,tau) ) * sign
!            clip with shallow water phase speed and CFL-phase speed
             csur  = sqrt(grav*Grd%ht(i,j))
             cmax  = - sign * min(Obc%bound(m)%ctrop_max * csur, 0.5 * var)
             cmin  = - sign * Obc%bound(m)%ctrop_min * csur 
             csur  = - sign * csur

             if (ABS(detai).lt. small) then
!             badly defined phases -> move this feature out with barotropic phase speed
               c1tmp = csur
             else
!             explicit time scheme for internal points
               c1tmp = - cgrid*(eta(i+sign,j,tau)-eta(i+sign,j,taum1))/detai 
             endif
!            If incoming waves are detected:
!            - if relaxation is done, Obc%bound(m)%ctrop(j) must be about zero    
!              -> boundary follows relaxation and generates incoming waves in the internal area
!            - if no relaxation is done, Obc%bound(m)%ctrop(j) is about cmin 
!              -> replace undesired incoming waves by outgoing waves
             if (c1tmp*sign .ge. 0.0) then ! incoming waves
               if(Obc%bound(m)%obc_relax_eta) then
                 c1tmp = 0     !boundary point follows relaxation as fast as possible 
               else
                 c1tmp = cmin  !there must not be incoming waves 
               endif
             else
               if(sign == 1) then
                 c1tmp = max(cmax,c1tmp)                    ! cmax is the phase speed limit
                 c1tmp = min(cmin,c1tmp)                    ! ensure minimum phase speed
               else
                 c1tmp = min(cmax,c1tmp)                    ! cmax is the phase speed limit
                 c1tmp = max(cmin,c1tmp)                    ! ensure minimum phase speed
               endif
             endif
             Obc%bound(m)%ctrop(j) = c1tmp 
          enddo

       else     ! south or north direction
          if(Obc%bound(m)%direction == NORTH) sign = -1
          j = Obc%bound(m)%js

          do i=Obc%bound(m)%is, Obc%bound(m)%ie
!            Use the barotropic timestep for the CFL-clipping,          
             var = Grd%dyu(i,j+sign)/dtfs
             cgrid = Grd%dyu(i,j+sign)/dteta
             detai =( eta(i,j+2*sign,tau) - eta(i,j+sign,tau) ) * sign
!            clip with shallow water phase speed and CFL-phase speed
             cmax  = - sign * min(Obc%bound(m)%ctrop_max*sqrt(grav*Grd%ht(i,j)), 0.5 * var)
             cmin  = - sign * Obc%bound(m)%ctrop_min * sqrt(grav*Grd%ht(i,j))

             if (ABS(detai).lt. small) then
!               badly defined phases -> move this feature out with maximum phase speed
                c1tmp = cmax
             else
!               implicit time scheme for internal points
                c1tmp = - cgrid*(eta(i,j+sign,tau)-eta(i,j+sign,taum1))/detai 
             endif
!            If incoming waves are detected:
!            - if relaxation is done, Obc%bound(m)%ctrop(j) must be  about zero    
!              -> boundary follows relaxation and generates incoming waves in the internal area
!            - if no relaxation is done, Obc%bound(m)%ctrop(j) is about cmin 
!              -> replace undesired incoming waves by outgoing waves
             if (c1tmp*sign .ge. 0.0) then
               if(Obc%bound(m)%obc_relax_eta) then
                 c1tmp = 0     !boundary point follows relaxation as fast as possible 
               else
                 c1tmp = cmin  !there must not be incoming waves
               endif
             else
               if(sign == 1) then 
                 c1tmp = max(cmax,c1tmp)                    ! cmax is the phase speed limit
                 c1tmp = min(cmin,c1tmp)                    ! ensure minimum phase speed
               else
                 c1tmp = min(cmax,c1tmp)                    ! cmax is the phase speed limit
                 c1tmp = max(cmin,c1tmp)                    ! ensure minimum phase speed
               endif
             endif
             Obc%bound(m)%ctrop(i) = c1tmp   
          enddo
       endif    
    enddo    
             
    return   
  end subroutine barotropic_phase_speed_init


  !###############################################################################
  subroutine ocean_obc_mass_flux(Time, Ext_mode, mass_flux)
    type(ocean_time_type), intent(in)          :: Time
    type(ocean_external_mode_type), intent(in) :: Ext_mode
    real, dimension(isd:,jsd:), intent(inout)  :: mass_flux
    integer                                    :: i, j, m, tau
    real                                       :: uh, uhjm, vh, vhim
    logical                                    :: used

    tau = Time%tau 

    mass_flux = 0.
    do m = 1, nobc
      if(.not. Obc%bound(m)%on_bound) cycle

      select case( Obc%bound(m)%direction )
      case(WEST)        
        i = Obc%bound(m)%is
        j = Obc%bound(m)%js-1
        uhjm = Ext_mode%ud(i,j,1,tau)*Grd%dyu(i,j)
        do j = Obc%bound(m)%js, Obc%bound(m)%je
          uh = Ext_mode%ud(i,j,1,tau)*Grd%dyu(i,j)
          mass_flux(i,j) = 0.5*(uh+uhjm)
          uhjm = uh
        enddo
      case(EAST)
        i = Obc%bound(m)%is-1
        j = Obc%bound(m)%js-1
        uhjm = Ext_mode%ud(i,j,1,tau)*Grd%dyu(i,j)
        do j = Obc%bound(m)%js, Obc%bound(m)%je
          uh = Ext_mode%ud(i,j,1,tau)*Grd%dyu(i,j)
          mass_flux(i,j) = - 0.5*(uh+uhjm)
          uhjm = uh
        enddo
      case(SOUTH)
        j = Obc%bound(m)%js
        i = Obc%bound(m)%is-1
        vhim = Ext_mode%ud(i,j,2,tau)*Grd%dxu(i,j)
        do i = Obc%bound(m)%is, Obc%bound(m)%ie
          vh = Ext_mode%ud(i,j,2,tau)*Grd%dxu(i,j)
          mass_flux(i,j) = 0.5*(vh+vhim)
          vhim = vh
        enddo
      case(NORTH)
        j = Obc%bound(m)%js-1
        i = Obc%bound(m)%is-1
        do i = Obc%bound(m)%is, Obc%bound(m)%ie
          vh = Ext_mode%ud(i,j,2,tau)*Grd%dxu(i,j)
          mass_flux(i,j) = - 0.5*(vh+vhim)
          vhim = vh
        enddo
      end select
    enddo

    !--- send out diagnostics data
    if(id_transport > 0) then
        used = send_data(id_transport, mass_flux(isc:iec,jsc:jec), &
               Time%model_time, rmask = Grd%umask(isc:iec,jsc:jec,1))
    endif
    
    return   
  end subroutine ocean_obc_mass_flux

  !#####################################################################

  subroutine ocean_obc_tracer_flux(Time, Tracer, tracer_flux, n, send_out)

    type(ocean_time_type), intent(in)         :: Time
    type(ocean_prog_tracer_type),  intent(in) :: Tracer
    real, dimension(isd:,jsd:), intent(inout) :: tracer_flux
    logical,                       intent(in) :: send_out
    integer,                       intent(in) :: n
    integer                                   :: i, j, k, m
    logical                                   :: used

    tracer_flux = 0.
    do m = 1, nobc
      if(.not. Obc%bound(m)%on_bound) cycle
      select case( Obc%bound(m)%direction )
      case(WEST)        
        i = Obc%bound(m)%is+1
        do k=1,nk
          do j = Obc%bound(m)%js, Obc%bound(m)%je
            tracer_flux(i,j) = tracer_flux(i,j) + Tracer%otf(m)%flux(j,k)
          enddo
        enddo
      case(EAST)
        i = Obc%bound(m)%is-1
        do k=1,nk
          do j = Obc%bound(m)%js, Obc%bound(m)%je
            tracer_flux(i,j) = tracer_flux(i,j) + Tracer%otf(m)%flux(j,k)
          enddo
        enddo
      case(SOUTH)
        j = Obc%bound(m)%js+1
        do k=1,nk
          do i = Obc%bound(m)%is, Obc%bound(m)%ie
            tracer_flux(i,j) = tracer_flux(i,j) + Tracer%otf(m)%flux(i,k)
          enddo
        enddo
      case(NORTH)
        j = Obc%bound(m)%js-1
        do k=1,nk
          do i = Obc%bound(m)%is, Obc%bound(m)%ie
            tracer_flux(i,j) = tracer_flux(i,j) + Tracer%otf(m)%flux(i,k)
          enddo
        enddo
      end select
    enddo

    !--- send out diagnostics data
    if(send_out .and. id_tracer_flux(n) > 0) then
        used = send_data(id_tracer_flux(n), Tracer%conversion*tracer_flux(isc:iec,jsc:jec), &
               Time%model_time, rmask = Grd%tmask(isc:iec,jsc:jec,1))
    endif
    
    return   
  end subroutine ocean_obc_tracer_flux

  !#####################################################################

  subroutine store_ocean_obc_tracer_flux(Tracer, flux_e, flux_n)
    type(ocean_prog_tracer_type), intent(inout) :: Tracer
    real, dimension(isd:,jsd:,:), intent(in)    :: flux_e, flux_n
    integer   :: i, j, m

    do m = 1, nobc
      if(.not. Obc%bound(m)%on_bound) cycle
      select case( Obc%bound(m)%direction )
      case(WEST)        
        i = Obc%bound(m)%is
        do j = Obc%bound(m)%js, Obc%bound(m)%je
          Tracer%otf(m)%flux(j,:) = Tracer%otf(m)%flux(j,:) + flux_e(i,j,:)
        enddo
      case(EAST)
        i = Obc%bound(m)%is-1
        do j = Obc%bound(m)%js, Obc%bound(m)%je
          Tracer%otf(m)%flux(j,:) = Tracer%otf(m)%flux(j,:) - flux_e(i,j,:) ! negative for eastward flux
        enddo
      case(SOUTH)
        j = Obc%bound(m)%js
        do i = Obc%bound(m)%is, Obc%bound(m)%ie
          Tracer%otf(m)%flux(i,:) = Tracer%otf(m)%flux(i,:) + flux_n(i,j,:)
        enddo
      case(NORTH)
        j = Obc%bound(m)%js-1
        do i = Obc%bound(m)%is, Obc%bound(m)%ie
          Tracer%otf(m)%flux(i,:) = Tracer%otf(m)%flux(i,:) - flux_n(i,j,:) ! negative for northward flux
        enddo
      end select
    enddo
    
    return   
  end subroutine store_ocean_obc_tracer_flux

  !#####################################################################

  subroutine store_ocean_obc_tracer_advect(flux_e, flux_n)
!   If vertical tracer advection is disabled, also all horizontal
!   advection must be removed for tracer conservation

    real, dimension(isd:,jsd:,:), intent(in)    :: flux_e, flux_n
    integer   :: i, j, m

    do m = 1, nobc
      if(.not. Obc%bound(m)%on_bound) cycle
      select case( Obc%bound(m)%direction )
      case(WEST)        
        i = Obc%bound(m)%is
        do j = Obc%bound(m)%js, Obc%bound(m)%je
          Obc%bound(m)%hadv(j,:) = (flux_e(i,j,:)-flux_e(i-1,j,:)) * Grd%datr(i,j)
        enddo
        if( .not. Obc%bound(m)%obc_vert_advel_t) then 
          do j = Obc%bound(m)%js, Obc%bound(m)%je
            Obc%bound(m)%hadv(j,:) = Obc%bound(m)%hadv(j,:) &
                                 + (flux_n(i,j,:)-flux_n(i,j-1,:)) * Grd%datr(i,j)
          enddo
        endif
      case(EAST)
        i = Obc%bound(m)%is
        do j = Obc%bound(m)%js, Obc%bound(m)%je
          Obc%bound(m)%hadv(j,:) = (flux_e(i,j,:)-flux_e(i-1,j,:)) * Grd%datr(i,j)
        enddo
        if( .not. Obc%bound(m)%obc_vert_advel_t) then 
          do j = Obc%bound(m)%js, Obc%bound(m)%je
            Obc%bound(m)%hadv(j,:) = Obc%bound(m)%hadv(j,:) &
                                 + (flux_n(i,j,:)-flux_n(i,j-1,:)) * Grd%datr(i,j)
          enddo
        endif
      case(SOUTH)
        j = Obc%bound(m)%js
        do i = Obc%bound(m)%is, Obc%bound(m)%ie
          Obc%bound(m)%hadv(i,:) = (flux_n(i,j,:)-flux_n(i,j-1,:)) * Grd%datr(i,j)
        enddo
        if( .not. Obc%bound(m)%obc_vert_advel_t) then 
          do j = Obc%bound(m)%js, Obc%bound(m)%je
            Obc%bound(m)%hadv(j,:) = Obc%bound(m)%hadv(j,:) &
                                 + (flux_e(i,j,:)-flux_e(i-1,j,:)) * Grd%datr(i,j)
          enddo
        endif
      case(NORTH)
        j = Obc%bound(m)%js
        do i = Obc%bound(m)%is, Obc%bound(m)%ie
          Obc%bound(m)%hadv(i,:) = (flux_n(i,j,:)-flux_n(i,j-1,:)) * Grd%datr(i,j)
        enddo
        if( .not. Obc%bound(m)%obc_vert_advel_t) then 
          do j = Obc%bound(m)%js, Obc%bound(m)%je
            Obc%bound(m)%hadv(j,:) = Obc%bound(m)%hadv(j,:) &
                                 + (flux_e(i,j,:)-flux_e(i-1,j,:)) * Grd%datr(i,j)
          enddo
        endif
      end select 
    enddo
    
    return   
  end subroutine store_ocean_obc_tracer_advect


  !#####################################################################
  ! Store the pressure gradient across the boundary.
  subroutine store_ocean_obc_pressure_grad(Thickness, pressure_gradient, tau)

   type(ocean_thickness_type), intent(in)            :: Thickness
   real, dimension(isc:iec,jsc:jec,nk,2), intent(in) :: pressure_gradient
   integer, intent(in)                               :: tau
   integer                                           :: i, j, k, m

    do m = 1, nobc
      if(.not. Obc%bound(m)%on_bound) cycle

      Obc%bound(m)%pgrad(:) = 0

      select case( Obc%bound(m)%direction )
      case(WEST)        
        i = Obc%bound(m)%is
        do j = Obc%bound(m)%js, Obc%bound(m)%je
          do k=1,nk ! also remove the coriolis term from the integral                                             
            Obc%bound(m)%pgrad(j) = Obc%bound(m)%pgrad(j) + &
                                      pressure_gradient(i,j,k,1) * &
                                      Thickness%dhu(i,j,k,tau) * Grd%umask(i,j,k)
          enddo
        enddo
      case(EAST)
        i = Obc%bound(m)%is-1
        do j = Obc%bound(m)%js, Obc%bound(m)%je
          do k=1,nk ! also remove the coriolis term from the integral                                             
            Obc%bound(m)%pgrad(j) = Obc%bound(m)%pgrad(j) + &
                                      pressure_gradient(i,j,k,1) * &
                                      Thickness%dhu(i,j,k,tau) * Grd%umask(i,j,k)
          enddo
        enddo
      case(SOUTH)
        j = Obc%bound(m)%js
        do i = Obc%bound(m)%is, Obc%bound(m)%ie
          do k=1,nk ! also remove the coriolis term from the integral                                             
            Obc%bound(m)%pgrad(i) = Obc%bound(m)%pgrad(j) + &
                                      pressure_gradient(i,j,k,2) * &
                                      Thickness%dhu(i,j,k,tau) * Grd%umask(i,j,k)
          enddo
        enddo
      case(NORTH)
        j = Obc%bound(m)%js
        do i = Obc%bound(m)%is, Obc%bound(m)%ie
          do k=1,nk ! also remove the coriolis term from the integral                                             
            Obc%bound(m)%pgrad(i) = Obc%bound(m)%pgrad(j) + &
                                      pressure_gradient(i,j,k,2) * &
                                      Thickness%dhu(i,j,k,tau) * Grd%umask(i,j,k)
          enddo
        enddo
      end select
    enddo
    
    return   
  end subroutine store_ocean_obc_pressure_grad

  !#####################################################################

end module ocean_obc_mod
