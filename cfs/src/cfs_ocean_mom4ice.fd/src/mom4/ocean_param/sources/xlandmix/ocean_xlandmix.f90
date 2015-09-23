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
module ocean_xlandmix_mod
!  
!<CONTACT EMAIL="Stephen.Griffies@noaa.gov"> S.M. Griffies 
!</CONTACT>
!
!<REVIEWER EMAIL="Matthew.Harrison@noaa.gov"> M.J. Harrison 
!</REVIEWER>
!
!<REVIEWER EMAIL="Keith.Dixon@noaa.gov"> K. Dixon 
!</REVIEWER>
!
!<OVERVIEW>
! Tracer and eta source from cross-land mixing
!</OVERVIEW>
!
!<DESCRIPTION>
! Compute thickness weighted tendency [tracer*meter/sec]
! of tracer associated with cross-land mixing of points. 
! In particular, this module provides interaction between bodies of water
! separated by land in coarse models (e.g., Med-Atlantic).  
! Interaction consists of mixing tracer and volume between water 
! columns found on each side of the land barrier.
! To conserve total tracer when using xlandmix with k=1, 
! also need to connect volumes between remote bodies of water
! via xlandmix of eta.  
!</DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! S.M. Griffies, M.J. Harrison, R. C. Pacanowski, and A. Rosati
! A Guide to MOM4 (2002)
! NOAA/Geophysical Fluid Dynamics Laboratory
! </REFERENCE>
!
! <NOTE>
! Algorithm ensures both total tracer and total volume is conserved.
! </NOTE>
!
! <NOTE>
! 2D domain decomposition is implemented according to following
! notions.  If the two points in an xlandmix pair live within halo 
! of each other (i.e., within same local_domain), 
! then no added mpp communication required. However, nore generally 
! the two points live further away.  In this case, xland_domain
! is defined so that its halos incorporate the maximum separation 
! of xlandmix points.  New tracer, eta, and grid arrays
! are defined over this extended xland_domain.  This added domain
! size will come at some computational cost, so it is advantageous
! to limit the separation between points within an xlandmix pair. 
! </NOTE>
!
! <NOTE>
! The current implementation of xlandinsert has not been generalized 
! to allow for communication between points separated by the tripolar fold.  
! The problem is in the logic used to compute xland_domain.  
! There is nothing fundamental limiting a generalization of the
! logic used to generate xland_domain.
! </NOTE>
!
! <NOTE>
! Many of the user specified values given in USER INPUT are
! model dependent since unresolved straits can become resolved 
! in finer mesh models. 
! </NOTE>
!
! <NOTE>
! Algorithm originally developed by Keith Dixon (kd@gfdl.noaa.gov) 
! for rigid lid and full cells and one processor (i.e., MOM1).  
! Algorithm extended to free surface and partial cells and 2D 
! domain decomposition by S.M.Griffies (smg@gfdl.noaa.gov).
! </NOTE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_xlandmix_nml">
!  <DATA NAME="use_xlandmix" TYPE="logical">
!  Needs to be true in order to use this scheme. 
!  </DATA> 
!  <DATA NAME="verbose_init" TYPE="logical">
!  For verbose initialization information.
!  </DATA> 
!</NAMELIST>

use constants_mod,     only: epsln
use diag_manager_mod,  only: register_diag_field, send_data
use field_manager_mod, only: MODEL_OCEAN, parse, find_field_index, get_field_methods, method_type, get_field_info
use fms_mod,           only: stdout, stdlog, FATAL, NOTE, WARNING
use fms_mod,           only: write_version_number, open_namelist_file, check_nml_error, close_file
use mpp_domains_mod,   only: domain2D, mpp_define_domains, mpp_update_domains
use mpp_domains_mod,   only: cyclic_global_domain, global_data_domain 
use mpp_mod,           only: mpp_error

use ocean_domains_mod, only: get_local_indices
use ocean_types_mod,   only: ocean_domain_type, ocean_grid_type
use ocean_types_mod,   only: ocean_thickness_type, ocean_time_type
use ocean_types_mod,   only: ocean_prog_tracer_type, ocean_external_mode_type


implicit none

private 

public ocean_xlandmix_init
private xlandmix_eta
public xlandmix

type(ocean_domain_type), pointer :: Dom => NULL()
type(ocean_grid_type),   pointer :: Grd => NULL()


! Variables set inside subroutine ocean_xlandmix_init
! See comments in ocean_xlandmix_init for description
integer, dimension(:,:), allocatable :: ixland     
integer, dimension(:,:), allocatable :: jxland     
integer, dimension(:,:), allocatable :: kxland     
real, dimension(:), allocatable      :: vxland     

! number of xland mixing pairs 
integer :: nxland=0 

! Variables calculated from user input
real, dimension(:,:), allocatable    :: depth_static   ! time independent depth of columns (m)
real, dimension(:,:), allocatable    :: gamma_xland    ! inverse damping time (s^-1)

integer                              :: xland_halo     ! halo size needed to perform xlandmix
integer                              :: kxland_min     ! min of the k-levels for the nxland pairs
integer                              :: kxland_max     ! max of the k-levels for the nxland pairs
type(ocean_domain_type), save        :: Xland_domain   ! domain for xlandmix 

! variables defined with halos given by xland_halo
real, dimension(:,:), allocatable    :: xland_eta      ! eta source associated with xlandmix (m/sec)
real, dimension(:,:), allocatable    :: xland_eta_pair ! eta source on space of xlandmix pairs (m/sec)
real, dimension(:,:,:), allocatable  :: tracerx        ! tracer concentration at taum1
real, dimension(:,:), allocatable    :: etax           ! surface height (m)
real, dimension(:,:), allocatable    :: xtx            ! zonal position of tracer points (degrees)
real, dimension(:,:), allocatable    :: ytx            ! meridional position of tracer points (degrees)
real, dimension(:,:), allocatable    :: datx           ! area of tracer cell (m^2)
real, dimension(:,:,:), allocatable  :: dhtx           ! partial cell thickness of tracer cell (m)

integer                              :: unit=6         !processor zero writes to unit 6

logical, dimension(:), allocatable   :: error_xland    !for checking that all xland points are OK.  

! for diagnostics 
logical :: used
integer, dimension(:), allocatable  :: id_xland
integer :: id_xland_eta=-1
integer :: num_prog_tracers=0

character(len=128) :: version=&
       '$Id$'
character (len=128) :: tagname = &
     '$Name$'

integer :: isd, ied, jsd, jed, isc, iec, jsc, jec, nk

integer :: ixl, jxl, ixl2, jxl2

logical :: verbose_init = .true.
logical :: use_xlandmix = .false.
logical :: module_is_initialized = .FALSE.

namelist /ocean_xlandmix_nml/ verbose_init, use_xlandmix 

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_xlandmix_init">
!
! <DESCRIPTION>
!    Initial set up for crossland mixing of tracers 
!
!    (i,j,k) locations of crossland mixing points are set via 
!    field table entries.
!
!    Time invariant crossland mixing rates (in units of m**3/sec) are 
!    set via field table entries. 
!
!    Checks are performed to ensure that the crossland mixing
!    grid locations are valid according to model configuration.
!
!    A summary of the locations of crossland mixing points is written out.
!
!    User specified inputs in "USER INPUT" section:
!
!    ixland and jxland = user specified nxland pairs of i,j grid locations
!
!    kxland = user specified uppermost (ktop=kxland(n,1)) and deepest 
!             (kbot=kxland(n,2)) levels over which crossland mixing will 
!             be done for each pair of crossland mixing points.
!
!    vxland = user specified time invariant rates of crossland mixing (m3/sec)
!             Equivalent to "the flow to the east = the flow to the west"
!             Dynamic vxland is not available in MOM4. 
!
!    NOTE: for ixland, jxland, and kxland pairs, values of the
!    (n,1) element should be < the corresponding (n,2) value.
! </DESCRIPTION>
!
subroutine ocean_xlandmix_init(Grid, Domain, Time, Thickness, T_prog, dtime)

  type(ocean_grid_type), intent(in), target   :: Grid
  type(ocean_domain_type), intent(in), target :: Domain
  type(ocean_time_type), intent(in)           :: Time
  type(ocean_thickness_type), intent(in)      :: Thickness
  type(ocean_prog_tracer_type), intent(in)    :: T_prog(:)
  real, intent(in)                            :: dtime

  integer :: n, model, imax, imin
  integer :: nxl, lx, xland_halo_test,i
  type(method_type), allocatable, dimension(:) :: xland_methods
  character(len=32) :: fld_type, fld_name
  integer :: parse_ok
  integer :: io_status, ioun, ierr
  real    :: ztop 

  if (module_is_initialized ) then 
    call mpp_error(FATAL, '==>Error in ocean_xlandmix_mod (ocean_xlandmix_init): module already initialized')
  endif 

  module_is_initialized = .TRUE.

  call write_version_number( version, tagname )

  call get_local_indices(Domain,isd,ied,jsd,jed,isc,iec,jsc,jec)
  nk = Grid%nk

  Dom => Domain
  Grd => Grid
  
  n = find_field_index(MODEL_OCEAN,'xland_mix')

  num_prog_tracers = size(T_prog(:))

  if (n < 1) return

  ! provide for namelist over-ride of default values
  ioun =  open_namelist_file()
  read (ioun,ocean_xlandmix_nml,IOSTAT=io_status)
  write (stdlog(),ocean_xlandmix_nml)
  write (stdout(),'(/)') 
  write (stdout(),ocean_xlandmix_nml)
  ierr = check_nml_error(io_status,'ocean_xlandmix_nml')
  call close_file (ioun)

  call get_field_info(n,fld_type,fld_name,model,nxland)

  if(use_xlandmix) then 
    call mpp_error(NOTE, '==>Note from ocean_xlandmix_mod (ocean_xlandmix_init): USING ocean_xlandmix_mod')
    write(stdout(),*)'==>Note from ocean_xlandmix_mod (ocean_xlandmix_init): USING ocean_xlandmix_mod'
    if(nxland == 0) then 
       call mpp_error(FATAL,'==>NOTE from ocean_xlandmix_mod: No cross-land mixing points chosen')
    endif
  else
    call mpp_error(NOTE, '==> Note from ocean_xlandmix_mod (ocean_xlandmix_init):: NOT using ocean_xlandmix_mod')
    write(stdout(),*)'==> Note from ocean_xlandmix_mod (ocean_xlandmix_init):: NOT using ocean_xlandmix_mod'
    if(nxland > 0) then 
      call mpp_error(WARNING, &
       '==>Warning: nxland > 0, but use_xlandmix=.false. Remove xland table entries to reduce memory')
    endif 
    return
  endif 

  write(stdout(),'(/a,f10.2)')'==> Note from ocean_xlandmix_mod: using forward time step of (secs)', dtime 

  allocate(xland_methods(nxland))

  allocate (ixland(nxland,2))
  ixland(:,:) = 0
  allocate (jxland(nxland,2))
  jxland(:,:) = 0
  allocate (kxland(nxland,2))
  kxland(:,:) = 0
  allocate (vxland(nxland))
  vxland(:)   = 0.0
  allocate (depth_static(nxland,2))
  depth_static(:,:) = 0.0
  allocate (gamma_xland(nxland,2))
  gamma_xland(:,:) = 0.0
  allocate (error_xland(nxland))
  error_xland(:) = .false.

  call get_field_methods(n,xland_methods)
  do i=1, nxland
     parse_ok = parse(xland_methods(i)%method_control,'ixland_1',ixland(i,1))
     if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_xlandmix_mod: xland table entry "ixland_1" error')
     parse_ok = parse(xland_methods(i)%method_control,'ixland_2',ixland(i,2))
     if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_xlandmix_mod: xland table entry "ixland_2" error')
     parse_ok = parse(xland_methods(i)%method_control,'jxland_1',jxland(i,1))
     if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_xlandmix_mod: xland table entry "jxland_1" error')
     parse_ok = parse(xland_methods(i)%method_control,'jxland_2',jxland(i,2))
     if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_xlandmix_mod: xland table entry "jxland_2" error')
     parse_ok = parse(xland_methods(i)%method_control,'kxland_1',kxland(i,1))
     if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_xlandmix_mod: xland table entry "kxland_1" error')
     parse_ok = parse(xland_methods(i)%method_control,'kxland_2',kxland(i,2))
     if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_xlandmix_mod: xland table entry "kxland_2" error')
     parse_ok   = parse(xland_methods(i)%method_control,'vxland',vxland(i))
     if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_xlandmix_mod: xland table entry "vxland" error')
  enddo

  if(Grid%tripolar) then 
     write(stdout(),'(a)')'==>Warning: xlandmix has not been implemented to connect points across the tripolar fold.'
  endif

  write(stdout(), *)'Using xlandmix to connect tracers and surface height between non-local ocean cells.'


  xland_halo = min(isc-isd,ied-iec,jsc-jsd,jed-jec)
  kxland_min = nk
  kxland_max = 0

  do nxl=1,nxland

     ! check for invalid crossland mixing grid locations

     if (kxland(nxl,1) > kxland(nxl,2)) then
        write (unit,994) nxl, kxland(nxl,1), nxl, kxland(nxl,2)
        error_xland(nxl) = .true.
     endif
     if (kxland(nxl,1) < kxland_min ) kxland_min = kxland(nxl,1)  
     if (kxland(nxl,2) > kxland_max ) kxland_max = kxland(nxl,2)  

     do lx = 1,2
        if (ixland(nxl,lx) < 1 .or. ixland(nxl,lx) > Grd%ni) then
           write(unit,991) nxl, lx, ixland(nxl,lx)
           error_xland(nxl) = .true.
        endif
        if (jxland(nxl,lx) < 1 .or. jxland(nxl,lx) > Grd%nj ) then
           write(unit,992) nxl, lx, jxland(nxl,lx)
           error_xland(nxl) = .true.
        endif
        if (kxland(nxl,lx) < 1 .or. kxland(nxl,lx) > Grd%nk ) then
           write(unit,993) nxl, lx, kxland(nxl,lx)
           error_xland(nxl) = .true.
        endif
     enddo

     ! determine size for xland_halo
     if(Grd%cyclic) then
       imax = max(ixland(nxl,1),ixland(nxl,2))
       imin = min(ixland(nxl,1),ixland(nxl,2))
       xland_halo_test = min( abs(imax-imin), abs(imax-Grd%ni-imin))
     else
       xland_halo_test = abs(ixland(nxl,1)-ixland(nxl,2))
     endif
     if (xland_halo_test > xland_halo ) xland_halo = xland_halo_test

     xland_halo_test = abs(jxland(nxl,1)-jxland(nxl,2))
     if (xland_halo_test > xland_halo ) xland_halo = xland_halo_test

  enddo  ! end of do-loop nxl=1,nxland

  if(xland_halo == 0) then 
     call mpp_error(FATAL,&
     '==>Error in ocean_xlandmix_mod (ocean_xlandmix_init): xland_halo=0 => problems defining xlandmix points.')
  endif
 
  if(xland_halo > Dom%xhalo .or. xland_halo > Dom%yhalo) then 
    write(stdout(),*)'Defining extra tracer arrays for xland_domain over k-levels ',kxland_min,' through ',kxland_max 
  endif 

  allocate (xland_eta_pair(nxland,2))
  xland_eta_pair(:,:)  = 0.0


  ! define arrays and xland_domain
  if(xland_halo <= Dom%xhalo .and. xland_halo <= Dom%yhalo) then 

     write(stdout(),*)'The model local computational domain has a x,y halo = ',Dom%xhalo,Dom%yhalo
     write(stdout(),*)'This is large enough for chosen points in xlandmix,'
     write(stdout(),*)'where we need to have a halo of size = ', xland_halo

     allocate (dhtx(isd:ied,jsd:jed,nk))
     allocate (datx(isd:ied,jsd:jed))
     allocate (xtx(isd:ied,jsd:jed))
     allocate (ytx(isd:ied,jsd:jed))
     xtx(:,:)    = Grd%xt(:,:)
     ytx(:,:)    = Grd%yt(:,:)
     dhtx(:,:,:) = Thickness%dht(:,:,:,1)
     datx(:,:)   = Grd%dat(:,:)

     allocate (xland_eta(isd:ied,jsd:jed))
     allocate (tracerx(1,1,1))
     allocate (etax(1,1))
     xland_eta(:,:) = 0.0
     tracerx(:,:,:) = 0.0
     etax(:,:)      = 0.0
  endif 

  if(xland_halo > Dom%xhalo .or. xland_halo > Dom%yhalo) then 
     write(stdout(),*)'The model local computational domain has a x,y halo = ',Dom%xhalo, Dom%yhalo
     write(stdout(),*)'This is smaller than the halo required for xlandmix.'
     write(stdout(),*)'For xlandmix, will define a new domain type with halo = ', xland_halo

     if(xland_halo > (Dom%xhalo+6)) then 
       write(stdout(),*)'=>Warning: ocean_xlandmix_init determined xland_halo = ', xland_halo
       write(stdout(),*)'           Are you sure you wish to connect such distant locations?'
       write(stdout(),*)'           Are you instead trying to connect points across a tripolar fold?'
       write(stdout(),*)'           If so, then presently this has yet to be coded into xlandmix.'
       write(stdout(),*)'           Alternative xlandmix points must be taken, or new code implemented.'
     endif

     call mpp_define_domains((/1,Grid%ni,1,Grid%nj/),Domain%layout,Xland_domain%domain2d, &
           xflags = CYCLIC_GLOBAL_DOMAIN, xhalo = xland_halo, yhalo = xland_halo, name='xlandmix')

     ! time independent arrays defined over the xland_domain region 
     allocate (dhtx(isc-xland_halo:iec+xland_halo,jsc-xland_halo:jec+xland_halo,nk))
     allocate (datx(isc-xland_halo:iec+xland_halo,jsc-xland_halo:jec+xland_halo))
     allocate (xtx(isc-xland_halo:iec+xland_halo,jsc-xland_halo:jec+xland_halo))
     allocate (ytx(isc-xland_halo:iec+xland_halo,jsc-xland_halo:jec+xland_halo))
     dhtx(:,:,:)    = 0.0
     datx(:,:)      = 0.0
     xtx(:,:)       = 0.0
     ytx(:,:)       = 0.0

     dhtx(isc:iec,jsc:jec,:)    = Thickness%dht(isc:iec,jsc:jec,:,1)
     datx(isc:iec,jsc:jec)      = Grd%dat(isc:iec,jsc:jec)
     xtx(isc:iec,jsc:jec)       = Grd%xt(isc:iec,jsc:jec)
     ytx(isc:iec,jsc:jec)       = Grd%yt(isc:iec,jsc:jec)
     call mpp_update_domains (dhtx(:,:,:), Xland_domain%domain2d)
     call mpp_update_domains (datx(:,:), Xland_domain%domain2d)
     call mpp_update_domains (xtx(:,:), Xland_domain%domain2d)
     call mpp_update_domains (ytx(:,:), Xland_domain%domain2d)

     ! time dependent arrays defined over the xland_domain region 
     allocate (xland_eta(isc-xland_halo:iec+xland_halo,jsc-xland_halo:jec+xland_halo))
     allocate (tracerx(isc-xland_halo:iec+xland_halo,jsc-xland_halo:jec+xland_halo,kxland_min:kxland_max))
     allocate (etax(isc-xland_halo:iec+xland_halo,jsc-xland_halo:jec+xland_halo))
     xland_eta(:,:) = 0.0
     tracerx(:,:,:) = 0.0
     etax(:,:)      = 0.0

  endif


  do nxl=1,nxland

     ! check for attempts to mix land rather than sea
     do lx = 1,2
        if(on_comp_domain(nxl, lx)) then
            ixl = ixland(nxl,lx)-Dom%ioff
            jxl = jxland(nxl,lx)-Dom%joff
           if (kxland(nxl,2) > Grid%kmt(ixl,jxl)) then
              write (unit,192) nxl, kxland(nxl,2), ixland(nxl,lx),jxland(nxl,lx),Grid%kmt(ixl,jxl)
              error_xland(nxl)=.true.
           endif
        endif
     enddo

     ! write out summary information for this pair of crossland mixing points
     if (kxland(nxl,1) == 1) then
        ztop = 0.0
     else
        ztop = Grid%zw(kxland(nxl,1)-1)
     endif

     if(at_least_one_in_comp_domain(nxl)) then
         ixl = ixland(nxl,1)-Dom%ioff
         jxl = jxland(nxl,1)-Dom%joff
         ixl2 = ixland(nxl,2)-Dom%ioff
         jxl2 = jxland(nxl,2)-Dom%joff
         if(verbose_init) then 
             write(unit,191) nxl, ixland(nxl,1), jxland(nxl,1) &
              ,xtx(ixl,jxl), ytx(ixl,jxl) &
              ,ixland(nxl,2), jxland(nxl,2) &
              ,xtx(ixl2,jxl2), ytx(ixl2,jxl2) &
              ,kxland(nxl,1), kxland(nxl,2), ztop, Grid%zw(kxland(nxl,2))
         endif 
         call xland_check (nxl, dtime, error_xland(nxl))
     endif

  enddo

  ! Bring model down here if there are problems with xlandmix points.  
  do nxl=1,nxland
    if(error_xland(nxl)) then
      call mpp_error(FATAL,'==>Error detected in ocean_xlandmix_mod: Use "grep -i error" to find ALL errors.')
    endif 
  enddo 

 
  ! register for diag_manager 
  allocate (id_xland(num_prog_tracers))
  id_xland = -1

  do n=1,num_prog_tracers
     if(T_prog(n)%name == 'temp') then 
         id_xland(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_xland', &
              Grd%tracer_axes(1:3), Time%model_time, &
              'rho_cp*xlandmix*dht*temp', &
              'Watt/m^2', missing_value=-1.e10, range=(/-1.e10,1.e10/))
     else
         id_xland(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_xland', &
              Grd%tracer_axes(1:3), Time%model_time, &
              'rho0*xlandmix*dht*tracer for'//trim(T_prog(n)%name),&
              trim(T_prog(n)%units)//' m*kg/sec', missing_value=-1.e10, range=(/-1.e10,1.e10/))
     endif
  enddo

  id_xland_eta  = register_diag_field ('ocean_model', 'eta_xland', &
                 Grd%tracer_axes(1:2), Time%model_time, &
                'eta xland source', 'm/sec', missing_value=-1e5, range=(/-1.e5,1.e5/))


191 format(/' ===== from ocean_xland_init ====='/                                  &
       ' for crossland sea connection pair number',i3,/                            &
       ' mix  (i,j) gridpoint (',i4,',',i4,') [long=',f8.3,' lat=',f8.3,']',/      &
       ' with (i,j) gridpoint (',i4,',',i4,') [long=',f8.3,' lat=',f8.3,']',/      & 
       ' from level',i3,' to',i3,' [depths of ',f10.3,' to ',f10.3,'m]' ) 
192 format(/'=>Error: Problem with crossland mixing: improper k-levels requested.'/&
       'kxland(',i2,',2) =',i4, ' exceeds kmt(',i4,',',i4,') = ',i4)
991 format(/'=>Error: problem with crossland tracer mixing:',/,                    &
       ' out of bounds i grid location requested'/  &
       '      ixland(',i4,',',i1,') was set incorrectly set to = ',i8,/)
992 format(/'=>Error: problem with crossland tracer mixing:',/,                    &
       ' out of bounds j-row grid location requested'/ &
       '      jxland(',i4,',',i1,') was set incorrectly set to = ',i8,/)
993 format(/'=>Error: problem with crossland tracer mixing:',/,                    &
       ' out of bounds k-level grid location requested'/ &
       '      kxland(',i4,',',i1,') was set incorrectly set to = ',i8,/)
994 format(/'=>Error: problem with crossland tracer mixing:',/,                    &
       ' improper k-values requested'/' kxland(',i3,',1)=',i5,                     &
       ' and is not less than or equal to kxland(',i3,',2)=',i5,/)

end subroutine ocean_xlandmix_init
! </SUBROUTINE> NAME="ocean_xlandmix_init"


!#######################################################################
! <SUBROUTINE NAME="xlandmix">
!
! <DESCRIPTION>
! Compute the sources for tracer concentration and free surface height.
! The two sources are related via requirements of tracer and volume 
! conservation. Logic is a bit tricky in order to allow each (i,j,k) point
! to participate in an arbitrary number of xlandmix pairs.  
!
! Compute thickness weighted tracer source with units tracer*meter/sec.
!
! For generalized vertical coordinates, need to update dthx each time 
! step since it may now generally have time dependent thicknesses. 
! This work remains to be done. 
!
! </DESCRIPTION>
!
subroutine xlandmix (Time, Thickness, T_prog, Ext_mode)

  type(ocean_time_type), intent(in)             :: Time
  type(ocean_thickness_type), intent(in)        :: Thickness
  type(ocean_prog_tracer_type), intent(inout)   :: T_prog(:)
  type(ocean_external_mode_type), intent(inout) :: Ext_mode
 
  real, dimension(2)  :: top_thickness       ! surface cell thickness at time taum1
  integer             :: k, nxl, lx, taum1, tau, n

  if (.not. module_is_initialized ) then 
    call mpp_error(FATAL, '==>Error in ocean_xlandmix_mod (xlandmix): module must be initialized')
  endif 

  if (nxland < 1) return
  if (.not. use_xlandmix) return

  taum1 = Time%taum1
  tau   = Time%tau

  xland_eta(:,:)      = 0.0
  xland_eta_pair(:,:) = 0.0
     
  do n = 1, num_prog_tracers

  T_prog(n)%wrk1(:,:,:) = 0.0

  if(xland_halo > Dom%xhalo .or. xland_halo > Dom%yhalo) then 

    do k=kxland_min,kxland_max
      tracerx(isc:iec,jsc:jec,k) = T_prog(n)%field(isc:iec,jsc:jec,k,taum1) 
    enddo
    call mpp_update_domains (tracerx(:,:,:), Xland_domain%domain2d)
    if(kxland_min==1) then
      etax(isc:iec,jsc:jec)   = Ext_mode%eta_t(isc:iec,jsc:jec,taum1)
      call mpp_update_domains (etax(:,:), Xland_domain%domain2d)
    endif 
  
    do nxl=1,nxland 

       if(at_least_one_in_comp_domain(nxl)) then  

           do lx=1,2
              ixl = ixland(nxl,lx)-Dom%ioff; jxl = jxland(nxl,lx)-Dom%joff
              gamma_xland(nxl,lx) = (2.0*vxland(nxl))/  &
                   (datx(ixl,jxl)*(depth_static(nxl,1)+depth_static(nxl,2)))
           enddo

           do lx=1,2 
              if(on_comp_domain(nxl,lx)) then
                 ixl = ixland(nxl,lx)-Dom%ioff
                 jxl = jxland(nxl,lx)-Dom%joff
                 ixl2 = ixland(nxl,3-lx)-Dom%ioff
                 jxl2 = jxland(nxl,3-lx)-Dom%joff
                 do k=max(2,kxland(nxl,1)),kxland(nxl,2)  
                    T_prog(n)%wrk1(ixl,jxl,k) = T_prog(n)%wrk1(ixl,jxl,k) &
                   + gamma_xland(nxl,lx)*(dhtx(ixl2,jxl2,k)*tracerx(ixl2,jxl2,k)-dhtx(ixl,jxl,k)*tracerx(ixl,jxl,k))
                 enddo
              endif
           enddo

           if(kxland(nxl,1)==1) then 

              k=1
              do lx=1,2
                 ixl = ixland(nxl,lx)-Dom%ioff
                 jxl = jxland(nxl,lx)-Dom%joff                 
                 top_thickness(lx) = Grd%dzt(k) + etax(ixl,jxl)
              enddo

              if(n==1) then 
                do lx=1,2
                   ixl = ixland(nxl,lx)-Dom%ioff
                   jxl = jxland(nxl,lx)-Dom%joff                   
                   xland_eta_pair(nxl,lx)  = gamma_xland(nxl,lx)*( top_thickness(3-lx)-top_thickness(lx) )
                   xland_eta(ixl,jxl)      = xland_eta(ixl,jxl) + xland_eta_pair(nxl,lx) 
                enddo
              endif 

              do lx=1,2  
                 if(on_comp_domain(nxl,lx)) then
                     ixl = ixland(nxl,lx)-Dom%ioff
                     jxl = jxland(nxl,lx)-Dom%joff
                     ixl2 = ixland(nxl,3-lx)-Dom%ioff
                     jxl2 = jxland(nxl,3-lx)-Dom%joff
                     T_prog(n)%wrk1(ixl,jxl,k) = T_prog(n)%wrk1(ixl,jxl,k) &
                     + gamma_xland(nxl,lx)*(top_thickness(3-lx)*tracerx(ixl2,jxl2,k)-top_thickness(lx)*tracerx(ixl,jxl,k)) 
                 endif
              enddo

           endif

       endif ! end of at_least_one_in_comp_domain(nxl) if-test 

    enddo  ! end of the nxland do-loop 

  endif ! end of xland_halo > halo if-test 

  if(xland_halo <= Dom%xhalo .and. xland_halo<= Dom%yhalo) then 

    do nxl=1,nxland 

       if(at_least_one_in_comp_domain(nxl)) then  

           do lx=1,2
              ixl = ixland(nxl,lx)-Dom%ioff
              jxl = jxland(nxl,lx)-Dom%joff 
              gamma_xland(nxl,lx) = (2.0*vxland(nxl))/  &
                   (Grd%dat(ixl,jxl)*(depth_static(nxl,1)+depth_static(nxl,2)))
           enddo

           do lx=1,2 
              if(on_comp_domain(nxl,lx)) then
                  ixl = ixland(nxl,lx)-Dom%ioff
                  jxl = jxland(nxl,lx)-Dom%joff
                  ixl2 = ixland(nxl,3-lx)-Dom%ioff
                  jxl2 = jxland(nxl,3-lx)-Dom%joff
                  do k=max(2,kxland(nxl,1)),kxland(nxl,2)  
                    T_prog(n)%wrk1(ixl,jxl,k) = T_prog(n)%wrk1(ixl,jxl,k)&
                         + gamma_xland(nxl,lx) &
                           *(Thickness%dht(ixl2,jxl2,k,tau)*T_prog(n)%field(ixl2,jxl2,k,taum1) &  
                            -Thickness%dht(ixl,jxl,k,tau)  *T_prog(n)%field(ixl,jxl,k,taum1))
                 enddo
              endif
           enddo

           if(kxland(nxl,1)==1) then 

              k=1
              do lx=1,2
                 ixl = ixland(nxl,lx)-Dom%ioff
                 jxl = jxland(nxl,lx)-Dom%joff
                 top_thickness(lx) = Grd%dzt(k) + Ext_mode%eta_t(ixl,jxl,taum1)
              enddo

              if(n==1) then  
                do lx=1,2
                   ixl = ixland(nxl,lx)-Dom%ioff
                   jxl = jxland(nxl,lx)-Dom%joff                   
                   xland_eta_pair(nxl,lx) = gamma_xland(nxl,lx)*( top_thickness(3-lx)-top_thickness(lx) )
                   xland_eta(ixl,jxl)     = xland_eta(ixl,jxl) + xland_eta_pair(nxl,lx)
                 enddo
              endif 

              do lx=1,2  
                 if(on_comp_domain(nxl,lx)) then
                     ixl = ixland(nxl,lx)-Dom%ioff
                     jxl = jxland(nxl,lx)-Dom%joff
                     ixl2 = ixland(nxl,3-lx)-Dom%ioff
                     jxl2 = jxland(nxl,3-lx)-Dom%joff
                     T_prog(n)%wrk1(ixl,jxl,k) = T_prog(n)%wrk1(ixl,jxl,k) &
                          +gamma_xland(nxl,lx)*(top_thickness(3-lx)*T_prog(n)%field(ixl2,jxl2,k,taum1)  &
                          -top_thickness(lx)*T_prog(n)%field(ixl,jxl,k,taum1))  
                 endif
              enddo

           endif

       endif ! end of at_least_one_in_comp_domain(nxl) if-test 

    enddo  ! end of the nxland do-loop 

  endif ! end of xland_halo > halo if-test 


! fill source array and send diagnostics
  T_prog(n)%th_tendency(isc:iec,jsc:jec,:) = T_prog(n)%th_tendency(isc:iec,jsc:jec,:) &
                                             + T_prog(n)%wrk1(isc:iec,jsc:jec,:)
  if(id_xland(n) > 0) then 
    used = send_data (id_xland(n), T_prog(n)%wrk1(:,:,:), &
           Time%model_time, rmask=Grd%tmask(:,:,:), &
           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif 


  enddo ! end of num_prog_tracers do-loop 

  ! compute effects on eta 
  call xlandmix_eta(Time, Ext_mode)

end subroutine xlandmix
! </SUBROUTINE> NAME="xlandmix"


!#######################################################################
! <SUBROUTINE NAME="xlandmix_eta">
!
! <DESCRIPTION>
! Compute the source for surface height with units meter/sec. 
! Note that xlandmix has already been called, so xland_eta has been 
! filled.  
! </DESCRIPTION>
!
subroutine xlandmix_eta (Time, Ext_mode)

  type(ocean_time_type), intent(in)             :: Time
  type(ocean_external_mode_type), intent(inout) :: Ext_mode
 
  integer :: nxl, lx

  if (.not. module_is_initialized ) then 
    call mpp_error(FATAL, '==>Error in ocean_xlandmix_mod (xlandmix_eta): module must be initialized')
  endif 

  Ext_mode%eta_source(isc:iec,jsc:jec) = Ext_mode%eta_source(isc:iec,jsc:jec) + xland_eta(isc:iec,jsc:jec)

  call mpp_update_domains (Ext_mode%eta_source(:,:), Dom%domain2d)

  if(id_xland_eta > 0) then 
    used = send_data (id_xland_eta, xland_eta(isc:iec,jsc:jec), &
           Time%model_time, rmask=Grd%tmask(isc:iec,jsc:jec,1)) 
  endif  

end subroutine xlandmix_eta 
! </SUBROUTINE> NAME="xlandmix_eta"


!#######################################################################
! <SUBROUTINE NAME="xland_check">
!
! <DESCRIPTION>
! Check if prescribed too much mixing
!
! In this routine the crossland mixing rate vxland(nxl) for
! the pair of points associated with index number nxl is
! converted into the fraction of the model grid boxes to be mixed
! per second, and checked.  These checks ensure that the rate of
! crossland mixing requested is valid in that it can be realized
! given the timestep length and column volumes involved.
! </DESCRIPTION>
!
! <IN NAME="nxl" TYPE="integer">
! Integer labeling the particular xlandmix pair
! </IN>
! <INOUT NAME="error" TYPE="logical">
! Error flag indicates whether initialization was performed successfully.  
! </INOUT>
!
subroutine xland_check(nxl, dtime, error)

  logical, intent(inout) :: error
  integer, intent(in)    :: nxl
  real, intent(in)       :: dtime

  integer                   :: k, lx, ll
  real, dimension(2)        :: colvol
  real                      :: factor, tslim
  real, dimension(nxland,2) :: bxland 

  if (.not. module_is_initialized ) then 
    call mpp_error(FATAL, '==>Error in ocean_xlandmix_mod (xland_check): module must be initialized')
  endif 

  ! Calculate depth range over which crossland mixing will be applied, 
  ! ignoring free surface for this diagnostic calculation.    
  do lx=1,2
    depth_static(nxl,lx) = 0.0
    do k=max(2,kxland(nxl,1)),kxland(nxl,2)
       ixl = ixland(nxl,lx)-Dom%ioff
       jxl = jxland(nxl,lx)-Dom%joff                            
       depth_static(nxl,lx) = depth_static(nxl,lx) + dhtx(ixl,jxl,k)
    enddo
    if(kxland(nxl,1) == 1) then
      depth_static(nxl,lx) = depth_static(nxl,lx) + Grd%dzt(1)
    endif
  enddo

  do lx=1,2

! For each of the two (i,j) points to be mixed, convert from volume
! to be mixed per second vxland(nxl) to model box fraction to be
! mixed per second bxland(nxl,l)
     ixl = ixland(nxl,lx)-Dom%ioff
     jxl = jxland(nxl,lx)-Dom%joff                          
     colvol(lx)     = datx(ixl,jxl)*depth_static(nxl,lx)
     bxland(nxl,lx) = vxland(nxl) / (epsln + colvol(lx))

    ! check for attempts to mix at a rate so fast that it can not
    ! be stably achieved in a single forward timestep
    ! (for the two boxes to be thoroughly mixed, no more than
    ! one-half of the volume of a given grid box can be
    ! stably transported in effect into the other gridbox)
    factor = bxland(nxl,lx)*dtime
    if (factor > 0.5) then
      tslim = 0.5/bxland(nxl,lx)
      write (unit,997) nxl, factor, ixland(nxl,lx), jxland(nxl,lx), bxland(nxl,lx), tslim, &
             vxland(nxl), kxland(nxl,1), kxland(nxl,2), colvol(lx)
      error=.true.
      return
    endif

  enddo

  if(verbose_init) then 
    write (unit,197) nxl, ixland(nxl,1), jxland(nxl,1), ixland(nxl,2), jxland(nxl,2), kxland(nxl,1), &
         kxland(nxl,2), depth_static(nxl,1), (colvol(ll),ll=1,2), vxland(nxl), (bxland(nxl,ll),ll=1,2)
  endif 

197   format(/' ===== from xlandvchk ====='/                             &
     ' for crossland sea communication pair number',i3,/                 &
     ' mix I,J gridpoints (',i4,',',i4,') and (',i4,',',i4,')'/          &
     ' from level',i3,' to',i3,' (a depth range of ',e12.6,' m)'/        &
     ' column volumes =',e12.6,' and ',e12.6,' m^3'/                     & 
     ' simulated flow in = flow out = ',e12.6,' m^3/sec,',/,             &
     ' so mix ',e12.6,' fraction of 1st column with ',e12.6,             &
     ' of 2nd column per sec'/)
997   format(/' Error => problem with crossland tracer mixing',          &
     ' pair #',i3,/,' Error => attempting to mix at too fast a rate:'/   &
     ' asking for more than half (',g12.6,') the volume of a grid box',  &
     ' at i, j, k: ',2i5,/,' to be mixed into its neighbor on',          &
     ' leapfrog timesteps'/                                              &  
     ' As specified now, fraction of column mixed per sec = ',e12.6,/,   &
     ' Potential Fixes: shorten dtime to be less than',e13.6,/           &
     t18,' reduce mixing rate vxland (now ',e12.6,' m^3/sec)'/           &
     t18,' or mix over a larger depth range (now from k=',i3,' to ',i3   &
     ,' => water column volume = ',e12.6,' m^3)'/)

end subroutine xland_check
! </SUBROUTINE> NAME="xland_check"


!#######################################################################
! <FUNCTION NAME="at_least_one_in_comp_domain">
!
! <DESCRIPTION>
! Function to see if at least one of the two points in a crossland pair
! is within the computational domain for the processor. 
! </DESCRIPTION>
!
! <IN NAME="nxl" TYPE="integer">
! Integer labeling the particular xlandmix pair
! </IN>
!
function at_least_one_in_comp_domain(nxl)

  integer, intent(in) :: nxl
  integer             :: lx
  logical             :: in_domain(2), at_least_one_in_comp_domain

  do lx=1,2
    if(isc+Dom%ioff <= ixland(nxl,lx) .and. ixland(nxl,lx) <= iec+Dom%ioff .and. &
       jsc+Dom%joff <= jxland(nxl,lx) .and. jxland(nxl,lx) <= jec+Dom%joff) then
       in_domain(lx) = .true.
    else 
       in_domain(lx) = .false.  
    endif
  enddo

  at_least_one_in_comp_domain = .false.
  if(in_domain(1) .and. in_domain(2)     ) at_least_one_in_comp_domain = .true.
  if(in_domain(1) .and. .not.in_domain(2)) at_least_one_in_comp_domain = .true.
  if(.not.in_domain(1) .and. in_domain(2)) at_least_one_in_comp_domain = .true.


end function at_least_one_in_comp_domain 
! </FUNCTION> NAME="at_least_one_in_comp_domain"


!#######################################################################
! <FUNCTION NAME="both_in_local_domain">
!
! <DESCRIPTION>
! Determine if both points in a crossland pair are within the
! local domain for the processor. 
! </DESCRIPTION>
!
! <IN NAME="nxl" TYPE="integer">
! Integer labeling the particular xlandmix pair
! </IN>
!
function both_in_local_domain(nxl)

  integer, intent(in) :: nxl
  integer             :: lx, ilo, ihi, jlo, jhi
  logical             :: in_domain(2), both_in_local_domain

  ilo = isd+Dom%ioff
  ihi = ied+Dom%ioff
  jlo = jsd+Dom%joff
  jhi = jed+Dom%joff

  do lx=1,2
    if(ilo <= ixland(nxl,lx) .and. ixland(nxl,lx) <= ihi .and. &
       jlo <= jxland(nxl,lx) .and. jxland(nxl,lx) <= jhi) then
       in_domain(lx) = .true.
    else 
       in_domain(lx) = .false.  
    endif
  enddo
  if(in_domain(1) .and. in_domain(2)) then
    both_in_local_domain = .true.
  else 
    both_in_local_domain = .false.
  endif

end function both_in_local_domain 
! </FUNCTION> NAME="both_in_local_domain"


!#######################################################################
! <FUNCTION NAME="on_comp_domain">
!
! <DESCRIPTION>
! Determine if the point is in comp-domain for the processor
! </DESCRIPTION>
!
! <IN NAME="nxl" TYPE="integer">
! Integer labeling the particular xlandmix pair
! </IN>
! <IN NAME="lx" TYPE="integer">
! lx=1,2 labels the point within an xlandmix pair
! </IN>
!
function on_comp_domain(nxl, lx)

  integer, intent(in) :: nxl, lx
  logical             :: on_comp_domain

  if(isc+Dom%ioff <= ixland(nxl,lx) .and. ixland(nxl,lx) <= iec+Dom%ioff .and. &
     jsc+Dom%joff <= jxland(nxl,lx) .and. jxland(nxl,lx) <= jec+Dom%joff) then
     on_comp_domain = .true.
  else 
     on_comp_domain = .false.  
  endif

end function on_comp_domain
! </FUNCTION> NAME="on_comp_domain"


!#######################################################################
! <FUNCTION NAME="on_processor">
!
! <DESCRIPTION>
! Determine if the point is on processor 
! </DESCRIPTION>
!
function on_processor(nxl, lx)

  integer, intent(in) :: nxl, lx
  logical             :: on_processor

  if(isc+Dom%ioff-xland_halo <= ixland(nxl,lx) .and. ixland(nxl,lx) <= iec+Dom%ioff+xland_halo .and. &
     jsd+Dom%joff-xland_halo <= jxland(nxl,lx) .and. jxland(nxl,lx) <= jed+Dom%joff+xland_halo) then
     on_processor = .true.
  else 
     on_processor = .false.  
  endif

end function on_processor
! </FUNCTION> NAME="on_processor"


end module ocean_xlandmix_mod
      
      




