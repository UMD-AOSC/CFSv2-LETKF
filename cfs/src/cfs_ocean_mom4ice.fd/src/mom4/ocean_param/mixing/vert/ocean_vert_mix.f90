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
module ocean_vert_mix_mod
!
!<CONTACT EMAIL="Ronald.Pacanowski@noaa.gov"> R. C. Pacanowski 
!</CONTACT>
!
!<REVIEWER EMAIL="Tony.Rosati@noaa.gov"> A. Rosati 
!</REVIEWER>
!
!<REVIEWER EMAIL="Stephen.Griffies@noaa.gov"> S.M. Griffies 
!</REVIEWER>
!
!<OVERVIEW>
! Time tendency from vertical mixing 
!</OVERVIEW>
!
!<DESCRIPTION>
! This module computes thickness weighted time tendency for tracer
! due to vertical diffusion processes, and the thickness weighted
! acceleration for velocity due to vertical friction processes.
! It also adds any background diffusivity.
!</DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! Kirk Bryan and L. J. Lewis
! A water mass model of the world ocean
! Journal of Geophysical Research (1979) vol 84, pages 2503--2517
! </REFERENCE>
!
! <NOTE> 
! The Bryan-Lewis vertical diffusivity is small in the upper ocean and 
! increases with depth according to an inverse tangent profile.  The default
! values are from roughly 0.05e-5 m^2/sec to roughly 1.0e-4 m^2/sec.
! </NOTE>

! </INFO>
!
!<NAMELIST NAME="ocean_vert_mix_nml">
!  <DATA NAME="bryan_lewis_diffusivity" TYPE="logical">
!  If .true. then add a Bryan-Lewis background to the 
!  diffusivity.  This background is a time-independent function
!  of depth.  
!  </DATA> 
!  <DATA NAME="bryan_lewis_lat_depend" TYPE="logical">
!  If .true. then allow for Bryan-Lewis background to be different 
!  outside of a tropical band than inside the band. 
!  </DATA> 
!  <DATA NAME="bryan_lewis_lat_transition" TYPE="real">
!  North/South latitude where transition from Bryan-Lewis values
!  in the tropic to those in the higher latitudes. 
!  </DATA> 
!  <DATA NAME="afkph_90, dfkph_90, sfkph_90, zfkph_90" UNITS="dimensionless" TYPE="real">
!  Parameters setting the Bryan-Lewis vertical diffusivity profile. 
!  When use bryan_lewis_lat_depend, these are the values used in the pole.
!  </DATA> 
!  <DATA NAME="afkph_00, dfkph_00, sfkph_00, zfkph_00" UNITS="dimensionless" TYPE="real">
!  Parameters setting the Bryan-Lewis vertical diffusivity profile in the tropics. 
!  When use bryan_lewis_lat_depend=.true. , these are the values used in the tropics.  
!  When use bryan_lewis_lat_depend=.false., these are the values used globally. 
!  </DATA> 
!  <DATA NAME="verbose_init" TYPE="logical">
!  For verbose writes during initialization. 
!  </DATA> 
!  <DATA NAME="use_diff_cbt_table" TYPE="logical">
!  If .true., then read in a table that specifies (i,j,ktop-->kbottom) 
!  and the diffusivity. This method is useful when aiming to mix vertically
!  at points where do cross-land insertion or where may wish to enhance 
!  mixing at river mouths.  
!  </DATA> 
!  <DATA NAME="linear_taper_diff_cbt_table" TYPE="logical">
!  If .true., then linear taper the diff_cbt_table value from 
!  so that it gets smaller with depth. 
!  </DATA> 
!</NAMELIST>

use constants_mod,       only: epsln, pi
use diag_manager_mod,    only: register_diag_field, register_static_field, send_data
use field_manager_mod,   only: MODEL_OCEAN, parse, find_field_index, get_field_methods, method_type, get_field_info
use fms_mod,             only: open_namelist_file, check_nml_error, close_file, write_version_number
use fms_mod,             only: FATAL, NOTE, stdout, stdlog
use mpp_mod,             only: mpp_error

use ocean_domains_mod,   only: get_local_indices
use ocean_types_mod,     only: ocean_grid_type, ocean_domain_type, ocean_thickness_type
use ocean_types_mod,     only: ocean_time_type, ocean_time_steps_type
use ocean_types_mod,     only: ocean_velocity_type, ocean_prog_tracer_type
use ocean_types_mod,     only: missing_value
use ocean_workspace_mod, only: wrk1, wrk2

implicit none

private

public ocean_vert_mix_init
public vert_frict
public vert_diffuse

! for Bryan-Lewis background vertical diffusivity profile 
logical :: vert_diffuse_on=.true.                                        ! must be true to use explicit vertical diffusion 
logical :: bryan_lewis_diffusivity=.false.                               ! depth dependent Bryan-Lewis vertical diffusivity 
logical :: bryan_lewis_lat_depend=.false.                                ! for different Bryan-Lewis in the tropics 
real    :: bryan_lewis_lat_transition=35.0                               ! latitude where center transition Bryan-Lewis values
real    :: afkph_00=0.55,dfkph_00=1.05,sfkph_00=4.5e-5,zfkph_00=2500.0e2 ! Bryan-Lewis parameters for equator (or global)
real    :: afkph_90=0.55,dfkph_90=1.05,sfkph_90=4.5e-5,zfkph_90=2500.0e2 ! Bryan-Lewis parameters for pole
real, dimension(:), allocatable :: diff_bryan_lewis_00                   ! Bryan-Lewis diffusivity at equator (or global)
real, dimension(:), allocatable :: diff_bryan_lewis_90                   ! Bryan-Lewis diffusivity at pole

! for reading in points where wish to set a static background diffusivity 
logical :: use_diff_cbt_table=.false.                     
logical :: linear_taper_diff_cbt_table=.false.
integer, dimension(:), allocatable   :: itable
integer, dimension(:), allocatable   :: jtable     
integer, dimension(:,:), allocatable :: ktable     
real, dimension(:), allocatable      :: diff_cbt_table

! allocatable array for static background vertical diffusivity
real, dimension(:,:,:), allocatable, public :: diff_cbt_back  

! for diagnostics 
logical :: used
integer :: id_vfrict_expl_u       =-1
integer :: id_vfrict_expl_v       =-1
integer :: id_diff_cbt_t          =-1
integer :: id_diff_cbt_s          =-1
integer :: id_diff_bryan_lewis_00 =-1
integer :: id_diff_bryan_lewis_90 =-1
integer :: id_diff_cbt_table      =-1
integer :: id_diff_cbt_back       =-1
integer :: id_fric_coeff          =-1
integer, dimension(:), allocatable  :: id_zflux_diff
integer, dimension(:), allocatable  :: id_v_diffuse

type(ocean_grid_type), pointer   :: Grd => NULL()
type(ocean_domain_type), pointer :: Dom => NULL()

integer :: isd, ied, jsd, jed, isc, iec, jsc, jec, nk

real, dimension(:,:,:), allocatable :: flux_z

character(len=128) :: version = &
     '$Id$'
character (len=128) :: tagname = &
     '$Name$'

logical :: module_is_initialized = .FALSE.
logical :: verbose_init          = .true.

namelist /ocean_vert_mix_nml/ vert_diffuse_on, verbose_init, use_diff_cbt_table, linear_taper_diff_cbt_table, &
                              bryan_lewis_diffusivity, bryan_lewis_lat_depend, bryan_lewis_lat_transition, &
                              afkph_90, dfkph_90, sfkph_90, zfkph_90, afkph_00, dfkph_00, sfkph_00, zfkph_00

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_vert_mix_init">
!
! <DESCRIPTION>
! Initialization for the vertical mixing module
! </DESCRIPTION>
!
subroutine ocean_vert_mix_init (Grid, Domain, Time, Time_steps, T_prog)

  type(ocean_grid_type), intent(in), target   :: Grid
  type(ocean_domain_type), intent(in), target :: Domain
  type(ocean_time_type), intent(in)           :: Time
  type(ocean_time_steps_type), intent(in)     :: Time_steps
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)

  integer :: ntable, ntable_points, model, parse_ok
  integer :: itbl, jtbl, ntbl
  integer :: ioun, io_status, ierr
  integer :: i, j, k, n, num_prog_tracers=0
  real    :: taper_table, thickness_table, weight  

  character(len=32) :: fld_type, fld_name
  type(method_type), allocatable, dimension(:) :: diff_cbt_table_methods

  if ( module_is_initialized ) then 
    call mpp_error(FATAL,'==>Error from ocean_vert_mix_mod (ocean_vert_mix_init)" module already initialized')
  endif 

  module_is_initialized = .TRUE.

  num_prog_tracers = size(T_prog(:))

  call write_version_number( version, tagname )

  ioun = open_namelist_file()
  read  (ioun, ocean_vert_mix_nml,iostat=io_status)
  write (stdout(),'(/)')
  write (stdout(), ocean_vert_mix_nml)  
  write (stdlog(), ocean_vert_mix_nml)
  ierr = check_nml_error(io_status, 'ocean_vert_mix_nml')
  call close_file(ioun)


  ! for reading in diffusivities from a table 
  ntable = find_field_index(MODEL_OCEAN,'diff_cbt_enhance')
  if (ntable < 1) then 
    call mpp_error(NOTE,'==>Warning: ocean_vert_mix_init did NOT find a table for enhanced diff_cbt. Will not read diffusivities.')  
  endif
  if(ntable > 1 .and. .not. use_diff_cbt_table) then 
    call mpp_error(NOTE,'==>Warning: ocean_vert_mix_init found field table for enhanced diff_cbt, yet use_diff_cbt_table=.false.')  
  endif 
  if(ntable > 1 .and. use_diff_cbt_table) then 
    write(stdout(),*)'==>Note: ocean_vert_mix_init will read background diffusivities from diff_cbt_table.'  
    if(linear_taper_diff_cbt_table) then 
       write(stdout(),*)'==>Note: ocean_vert_mix_init will linearly taper diff_cbt_table values to zero with depth. '
    endif
  endif 

  Dom => Domain
  Grd => Grid

  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk = Grid%nk

  ! register acceleration due to vertical friction for diagnostic output 
  id_vfrict_expl_u = register_diag_field('ocean_model','vfrict_expl_u',Grd%vel_axes_uv(1:3), Time%model_time, &
       'explicit vert friction on u-velocity', 'm^2/s^2', missing_value=-1e10, range=(/-1e10,1e10/))
  id_vfrict_expl_v = register_diag_field('ocean_model','vfrict_expl_v',Grd%vel_axes_uv(1:3), Time%model_time, &
       'explicit vert friction on v-velocity', 'm^2/s^2', missing_value=-1e10, range=(/-1e10,1e10/))

  ! register vertical viscosity 
  id_fric_coeff  = register_diag_field ('ocean_model', 'vert_fric_coeff', Grid%vel_axes_wu(1:3), Time%model_time, &
       'total vertical viscosity', 'm^2/s',missing_value=-10.0, range=(/-10.0,1e6/))

  if(.not. vert_diffuse_on) then 
     write(stdout(),*)'=>Note: Not using explicit vertical diffusion of tracer.'
     return
  endif 

  allocate (diff_bryan_lewis_00(nk))
  diff_bryan_lewis_00(:) = 0.0 
  allocate (diff_bryan_lewis_90(nk))
  diff_bryan_lewis_90(:) = 0.0 

  allocate (diff_cbt_back(isd:ied,jsd:jed,nk))
  diff_cbt_back(:,:,:) = 0.0 

  allocate(flux_z(isd:ied,jsd:jed,nk))
  flux_z = 0.0

  if(bryan_lewis_diffusivity) then
      write(stdout(),*)'=>Note: USING Bryan Lewis vertical background tracer diffusivity.  Scheme adds a diffusivity'
      write(stdout(),*)'        on top of all other diffusivity schemes.  For example, if using constant vertical'
      write(stdout(),*)'        diffusivity, will have full diffusivity given by diff_cbt=kappa_h+diff_bryan_lewis'
      write(stdout(),*)'        If ONLY want the Bryan-Lewis background, set kappa_h=0.0 in namelist.'

      if(bryan_lewis_lat_depend) then 
          write(stdout(),*)'=>Note: USING bryan_lewis_lat_depend.  This will modify the background latitudinally, '
          write(stdout(),*)'        with a transition latitude bryan_lewis_lat_transition = ',bryan_lewis_lat_transition
      endif

      do k=1,nk

         diff_bryan_lewis_90(k) = 1.e-4*(afkph_90 + (dfkph_90/pi)*(atan(sfkph_90*(100.0*Grid%zw(k) - zfkph_90))))
         diff_bryan_lewis_00(k) = 1.e-4*(afkph_00 + (dfkph_00/pi)*(atan(sfkph_00*(100.0*Grid%zw(k) - zfkph_00))))

         if(diff_bryan_lewis_00(k) <= 0.0) then
             call mpp_error(FATAL,'==>Error in ocean_vert_mix_mod: problems w/ Bryan-Lewis: parameters give diffusivities < 0')
         endif
         if ((Time_steps%dtime_t*diff_bryan_lewis_00(k))/Grid%dzt(k)**2 >= 0.5  .and. Time_steps%aidif /= 1.0) then
             write (stdout(),'(a,a,i3)')'==> Warning: vertical diffusive criteria exceeded for "diff_bryan_lewis".',&
                  ' use a smaller "dtts" and/or "diff_bryan_lewis" at level k=',k
         endif
      enddo

      if(verbose_init) then 
          write(stdout(),*) ' '
          do k=1,nk
             write(stdout(),'(a,i4,a,e14.7)')'diff_bryan_lewis_00(m^2/sec)(',k,') = ',diff_bryan_lewis_00(k)
          enddo
          if(bryan_lewis_lat_depend) then 
             write(stdout(),*) ' '
             do k=1,nk
                write(stdout(),'(a,i4,a,e14.7)')'diff_bryan_lewis_90(m^2/sec)(',k,') = ',diff_bryan_lewis_90(k)
             enddo
          endif  
      endif

      ! add to background diffusivity 
      if(bryan_lewis_lat_depend) then 
          do k=1,nk
             do j=jsc,jec
                do i=isc,iec
                   weight = (1+tanh((abs(Grd%yt(i,j))-bryan_lewis_lat_transition)/5.0))/2.0
                   diff_cbt_back(i,j,k) = diff_cbt_back(i,j,k) +  &
                        Grd%tmask(i,j,k)*((1-weight)*diff_bryan_lewis_00(k) + weight*diff_bryan_lewis_90(k))
                enddo
             enddo
          enddo
      else 
          do k=1,nk
             do j=jsc,jec
                do i=isc,iec
                   diff_cbt_back(i,j,k) = diff_cbt_back(i,j,k) +  Grd%tmask(i,j,k)*diff_bryan_lewis_00(k)
                enddo
             enddo
          enddo
      endif

  endif


! read in static diffusivities specified from the field table "diff_cbt_table"
  wrk2(:,:,:) = 0.0 
  if(use_diff_cbt_table .and. ntable > 1) then 

      call get_field_info(ntable,fld_type,fld_name,model,ntable_points)
      allocate(diff_cbt_table_methods(ntable_points))

      allocate (itable(ntable_points))
      itable(:) = 0
      allocate (jtable(ntable_points))
      jtable(:) = 0
      allocate (ktable(ntable_points,2))
      ktable(:,:) = 0
      allocate (diff_cbt_table(ntable_points))
      diff_cbt_table(:) = 0.0

      call get_field_methods(ntable,diff_cbt_table_methods)
      do ntbl=1,ntable_points
         parse_ok = parse(diff_cbt_table_methods(ntbl)%method_control,'itable',itable(ntbl))
         if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_vert_mix_init: diff_cbt_table entry "itable" error')
         parse_ok = parse(diff_cbt_table_methods(ntbl)%method_control,'jtable',jtable(ntbl))
         if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_vert_mix_init: diff_cbt_table entry "jtable" error')
         parse_ok = parse(diff_cbt_table_methods(ntbl)%method_control,'ktable_1',ktable(ntbl,1))
         if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_vert_mix_init: diff_cbt_table entry "ktable_1" error')
         parse_ok = parse(diff_cbt_table_methods(ntbl)%method_control,'ktable_2',ktable(ntbl,2))
         if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_vert_mix_init: diff_cbt_table entry "ktable_2" error')
         parse_ok   = parse(diff_cbt_table_methods(ntbl)%method_control,'diff_cbt_table',diff_cbt_table(ntbl))
         if (parse_ok == 0) call mpp_error(FATAL,'==>Error ocean_vert_mix_init: diff_cbt_table entry "diff_cbt_table" error')
         if(ktable(ntbl,2) < ktable(ntbl,1) ) then 
             call mpp_error(FATAL,'==>Error ocean_vert_mix_init: diff_cbt_table entry must have ktable_2 >= ktable_1')
         endif
      enddo

      do ntbl=1,ntable_points

         if(on_comp_domain(ntbl)) then

             itbl=itable(ntbl)-Dom%ioff
             jtbl=jtable(ntbl)-Dom%joff

             if(linear_taper_diff_cbt_table) then  
                 thickness_table = Grd%zw(ktable(ntbl,2)) - Grd%zw(ktable(ntbl,1))
             endif

             do k=ktable(ntbl,1),ktable(ntbl,2)

                if(Grd%tmask(itbl,jtbl,k) == 0.0) then 
                    write(*,'(a,i4,a,i4,a,i4,a)')'==>Warning from ocean_vert_mix_init: ignored nonzero diff_cbt_table(', &
                         itable(ntbl),',',jtable(ntbl),',',k,') set over land. Is your table correct?'
                endif

                if(linear_taper_diff_cbt_table .and. thickness_table > 0.0) then              
                    taper_table = (Grd%zw(ktable(ntbl,2)) -Grd%zw(k))/thickness_table 
                else 
                    taper_table = 1.0
                endif

                wrk2(itbl,jtbl,k)          = taper_table*Grd%tmask(itbl,jtbl,k)*diff_cbt_table(ntbl)  
                diff_cbt_back(itbl,jtbl,k) = diff_cbt_back(itbl,jtbl,k) + wrk2(itbl,jtbl,k)

             enddo

         endif

      enddo

  else 

      ntable_points=1 
      allocate (itable(ntable_points))
      itable(:) = 0
      allocate (jtable(ntable_points))
      itable(:) = 0
      allocate (ktable(ntable_points,2))
      ktable(:,:) = 0
      allocate (diff_cbt_table(ntable_points))
      diff_cbt_table(:) = 0.0

  endif


  ! register thickness weighted vertical diffusion tendency for diagnostic output 
  allocate (id_zflux_diff(num_prog_tracers))
  allocate (id_v_diffuse(num_prog_tracers))
  id_zflux_diff=-1
  id_v_diffuse =-1

  do n=1,num_prog_tracers

     if (trim(T_prog(n)%name) == 'temp') then
         id_zflux_diff(n) = register_diag_field ('ocean_model', &
              trim(T_prog(n)%name)//'_zflux_diff', Grd%tracer_axes(1:3), &
              Time%model_time, 'z-diffusive heat flux', 'Watts/m^2',&
              missing_value=missing_value, range=(/-1.e16,1.e16/))
         id_v_diffuse(n) = register_diag_field ('ocean_model', &
              trim(T_prog(n)%name)//'_v_diffuse', Grd%tracer_axes(1:3), &
              Time%model_time, 'thk wghtd expl vert-diff of heat', 'Watts/m^2',&
              missing_value=missing_value, range=(/-1.e16,1.e16/)) 
     else 
         id_zflux_diff(n) = register_diag_field ('ocean_model',&
         trim(T_prog(n)%name)//'_zflux_diff', Grd%tracer_axes(1:3), &
              Time%model_time, 'k-diffusive flux of '//trim(T_prog(n)%name), 'm*kg/sec',&
              missing_value=missing_value, range=(/-1.e10,1.e10/))
         id_v_diffuse(n) = register_diag_field ('ocean_model',&
         trim(T_prog(n)%name)//'_v_diffuse', Grd%tracer_axes(1:3), &
              Time%model_time, 'thk wghtd expl vert-diff of '//trim(T_prog(n)%name), 'm*kg/sec',&
              missing_value=missing_value, range=(/-1.e10,1.e10/))         
    endif  

  enddo
  
  ! register vertical diffusivity
  id_diff_cbt_t = -1
  id_diff_cbt_t    = register_diag_field ('ocean_model', 'diff_cbt_t', Grid%tracer_axes(1:3), Time%model_time, &
       'total vert diff_cbt(temp) (w/o neutral included)', 'm^2/s',missing_value=-10.0, range=(/-10.0,1e6/))

  id_diff_cbt_s = -1
  id_diff_cbt_s    = register_diag_field ('ocean_model', 'diff_cbt_s', Grid%tracer_axes(1:3), Time%model_time, &
       'total vert diff_cbt(salt) (w/o neutral included)', 'm^2/s',missing_value=-10.0, range=(/-10.0,1e6/))

  ! register and send static diffusivities 

  id_diff_bryan_lewis_00 = -1
  id_diff_bryan_lewis_00 = register_static_field ('ocean_model', 'diff_bryan_lewis_00', Grid%tracer_axes(1:3), &
       'Bryan-Lewis vertical diffusivity 00', 'm^2/s',missing_value=-10.0, range=(/-10.0,1e6/))
  if (id_diff_bryan_lewis_00 > 0) then 
      do j=jsc,jec
         do i=isc,iec
            wrk1(i,j,:) = diff_bryan_lewis_00(:)*Grd%tmask(i,j,:)
         enddo
      enddo
      used = send_data (id_diff_bryan_lewis_00, wrk1(isc:iec,jsc:jec,:), &
             Time%model_time, rmask=Grd%tmask(isc:iec,jsc:jec,:))
  endif

  id_diff_bryan_lewis_90 = -1
  id_diff_bryan_lewis_90 = register_static_field ('ocean_model', 'diff_bryan_lewis_90', Grid%tracer_axes(1:3), &
       'Bryan-Lewis vertical diffusivity 90', 'm^2/s',missing_value=-10.0, range=(/-10.0,1e6/))
  if (id_diff_bryan_lewis_90 > 0) then 
      do j=jsc,jec
         do i=isc,iec
            wrk1(i,j,:) = diff_bryan_lewis_90(:)*Grd%tmask(i,j,:)
         enddo
      enddo
      used = send_data (id_diff_bryan_lewis_90, wrk1(isc:iec,jsc:jec,:), &
             Time%model_time, rmask=Grd%tmask(isc:iec,jsc:jec,:))
  endif

  id_diff_cbt_table = -1
  id_diff_cbt_table = register_static_field ('ocean_model', 'diff_cbt_table', Grid%tracer_axes(1:3), &
       'diff_cbt from table', 'm^2/s',missing_value=-10.0, range=(/-10.0,1e6/))
  if (id_diff_cbt_table > 0) then 
    used = send_data (id_diff_cbt_table, wrk2(isc:iec,jsc:jec,:), &
           Time%model_time, rmask=Grd%tmask(isc:iec,jsc:jec,:))
  endif

  id_diff_cbt_back = -1
  id_diff_cbt_back = register_static_field ('ocean_model', 'diff_cbt_back', Grid%tracer_axes(1:3), &
       'static background diff_cbt', 'm^2/s',missing_value=-10.0, range=(/-10.0,1e6/))
  if (id_diff_cbt_back > 0) then 
    used = send_data (id_diff_cbt_back, diff_cbt_back(isc:iec,jsc:jec,:), &
           Time%model_time, rmask=Grd%tmask(isc:iec,jsc:jec,:))
  endif



end subroutine ocean_vert_mix_init
! </SUBROUTINE> NAME="ocean_vert_mix_init"


!#######################################################################
! <SUBROUTINE NAME="vert_diffuse">
!
! <DESCRIPTION>
! This subroutine computes the thickness weighted time tendency 
! for tracer associated with vertical diffusion. 
! Only support aidif==0.0 or aidif==1.0.
! Do not support case with 0.0 < aidif < 1.0.
! </DESCRIPTION>
!
subroutine vert_diffuse (Time, Thickness, aidif, ntracer, Tracer, Temp, Salt, diff_cbt, diag_flag) 

  type(ocean_time_type), intent(in)                  :: Time 
  type(ocean_thickness_type), intent(in)             :: Thickness
  real, intent(in)                                   :: aidif              ! aidif=1.0=>implicit.  aidif=0.0=>explicit
  integer, intent(in)                                :: ntracer            ! integer for the tracer 
  type(ocean_prog_tracer_type), intent(inout)        :: Tracer             ! Tracer for vertical diffusion calculation
  type(ocean_prog_tracer_type), intent(inout)        :: Temp               ! temp type 
  type(ocean_prog_tracer_type), intent(inout)        :: Salt               ! salinity type 
  real, intent(in), dimension(isd:ied,jsd:jed,nk,2)  :: diff_cbt           !vertical diffusivity (m^2/sec)
  logical, intent(in), optional                      :: diag_flag          
  logical                                            :: send_diagnostics 

  real,dimension(isd:ied,jsd:jed)                    :: ft1, ft2
  integer                                            :: taum1, nmix
  integer                                            :: i, j, k, n, nt2, kp

  if(.not. vert_diffuse_on) then 
    return 
  endif 

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL,'==>Error from ocean_vert_mix_mod (vert_diffuse): module must be initialized')
  endif 

  ! assume send_diagnostics=.true., unless diag_flag says it is false.
  if (PRESENT(diag_flag)) then 
    send_diagnostics = diag_flag
  else 
    send_diagnostics = .true.
  endif 

  taum1 = Time%taum1

  if (Tracer%name == 'salt') then
     nmix = 2
  else
     nmix = 1
  endif

  ! fully implicit
  if (aidif==1.0) then  
    wrk1 = 0.0

  ! fully explicit
  elseif (aidif < 1.0) then

    ft1(isc:iec,jsc:jec) = Tracer%stf(isc:iec,jsc:jec)
    k = 1
    kp = min(k+1,nk)
    ft2(isc:iec,jsc:jec) = Grd%tmask(isc:iec,jsc:jec,kp)*diff_cbt(isc:iec,jsc:jec,k,nmix)*&
                           (Tracer%field(isc:iec,jsc:jec,k,taum1)-Tracer%field(isc:iec,jsc:jec,kp,taum1))&
                           /Thickness%dhwt(isc:iec,jsc:jec,k)
    wrk1(isc:iec,jsc:jec,k) = (ft1(isc:iec,jsc:jec)-ft2(isc:iec,jsc:jec))
    ft1(isc:iec,jsc:jec) = ft2(isc:iec,jsc:jec)
    do k=2,nk
      kp = min(k+1,nk)
      ft2(isc:iec,jsc:jec) = Grd%tmask(isc:iec,jsc:jec,kp)*diff_cbt(isc:iec,jsc:jec,k,nmix)*&
                           (Tracer%field(isc:iec,jsc:jec,k,taum1)-Tracer%field(isc:iec,jsc:jec,kp,taum1))&
                           /Thickness%dhwt(isc:iec,jsc:jec,k)
      if (id_zflux_diff(ntracer) > 0) then
          flux_z(isc:iec,jsc:jec,k) = ft2(isc:iec,jsc:jec)*Tracer%conversion
      endif
      wrk1(isc:iec,jsc:jec,k) = (ft1(isc:iec,jsc:jec)-ft2(isc:iec,jsc:jec))
      ft1(isc:iec,jsc:jec) = ft2(isc:iec,jsc:jec)
    enddo
  endif

  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           Tracer%th_tendency(i,j,k) = Tracer%th_tendency(i,j,k) + wrk1(i,j,k)
        enddo
     enddo
  enddo

  if(send_diagnostics) then 
    if (id_diff_cbt_t > 0 .and. Tracer%name=='temp') then 
      used = send_data (id_diff_cbt_t, diff_cbt(:,:,:,nmix), &
             Time%model_time, rmask=Grd%tmask(:,:,:), &
             is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
    endif 
    if (id_diff_cbt_s > 0 .and. Tracer%name=='salt') then 
      used = send_data (id_diff_cbt_s, diff_cbt(:,:,:,nmix), &
             Time%model_time, rmask=Grd%tmask(:,:,:), &
             is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
    endif 
    if (id_zflux_diff(ntracer) > 0) then 
      used = send_data (id_zflux_diff(ntracer), flux_z(:,:,:), &
             Time%model_time, rmask=Grd%tmask(:,:,:), &
             is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
    endif 
    if (id_v_diffuse(ntracer) > 0) then
        used = send_data( id_v_diffuse(ntracer), &
               wrk1(isc:iec,jsc:jec,:)*Tracer%conversion, Time%model_time, &
               rmask = Grd%tmask(isc:iec,jsc:jec,:))           
    endif
  endif 

end subroutine vert_diffuse
! </SUBROUTINE> NAME="vert_diffuse"


!#######################################################################
! <FUNCTION NAME="vert_frict">
!
! <DESCRIPTION>
! This subroutine computes the thickness weighted acceleration (m^2/s^2)
! associated with vertical friction.
! </DESCRIPTION>
!
function vert_frict (Time, Thickness, aidif, Velocity, visc_cbu, diag_flag)

  type(ocean_time_type), intent(in)               :: Time
  type(ocean_thickness_type), intent(in)          :: Thickness
  real, intent(in)                                :: aidif            !aidif=1.0=>implicit.  aidif=0.0=>explicit
  type(ocean_velocity_type), intent(in)           :: Velocity         !velocity type 
  real, intent(in), dimension(isd:ied,jsd:jed,nk) :: visc_cbu         !vertical viscosity (m^2/s)
  logical, intent(in), optional                   :: diag_flag 
  logical                                         :: send_diagnostics 

  real,dimension(isd:ied,jsd:jed)                 :: ft1, ft2
  real,dimension(isc:iec,jsc:jec,nk,2)            :: vert_frict 
  integer                                         :: k, i, j, n, taum1

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL,'==>Error from ocean_vert_mix_mod (vert_frict): module must be initialized')
  endif 

  ! assume send_diagnostics=.true., unless diag_flag says it is false.
  if (PRESENT(diag_flag)) then 
    send_diagnostics = diag_flag
  else 
    send_diagnostics = .true.
  endif

  taum1 = Time%taum1

  if (aidif == 1.0) then     ! fully implicit

      vert_frict(isc:iec,jsc:jec,:,:) = 0.0

  elseif (aidif < 1.0) then ! explicit

      do n=1,2

         do j=jsc,jec
            do i=isc,iec
               ft1(i,j) = Velocity%smf(i,j,n)
            enddo
         enddo

         do k=1,nk-1
            do j=jsc,jec
               do i=isc,iec
                  ft2(i,j) = (1.0-Grd%umask(i,j,k+1))*Velocity%bmf(i,j,n) &
                             + Grd%umask(i,j,k+1)*visc_cbu(i,j,k)         &
                             *(Velocity%u(i,j,k,n,taum1)-Velocity%u(i,j,k+1,n,taum1))/Thickness%dhwu(i,j,k)
                  vert_frict(i,j,k,n) = Grd%umask(i,j,k)*(ft1(i,j)-ft2(i,j))
                  ft1(i,j) = ft2(i,j)
               enddo
            enddo
         enddo
         do j=jsc,jec
            do i=isc,iec
               vert_frict(i,j,nk,n) = Grd%umask(i,j,nk)*(ft1(i,j)-Velocity%bmf(i,j,n))
            enddo
         enddo

      enddo

  endif

  if (aidif > 0.0 .and. aidif < 1.0) then ! mixed implicit/explicit
    vert_frict(isc:iec,jsc:jec,:,:) = vert_frict(isc:iec,jsc:jec,:,:)*(1.0-aidif)
  endif

  if(send_diagnostics) then 
     if (id_fric_coeff > 0) used = send_data (id_fric_coeff, visc_cbu(:,:,:), &
                                   Time%model_time, rmask=Grd%umask(:,:,:), &
                                   is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

     if (id_vfrict_expl_u > 0) used = send_data (id_vfrict_expl_u, vert_frict(isc:iec,jsc:jec,:,1), &  
                                      Time%model_time, rmask=Grd%umask(isc:iec,jsc:jec,:))

     if (id_vfrict_expl_v > 0) used = send_data (id_vfrict_expl_v, vert_frict(isc:iec,jsc:jec,:,2), &
                                      Time%model_time, rmask=Grd%umask(isc:iec,jsc:jec,:))
  endif 

end function vert_frict
! </FUNCTION> NAME="vert_frict"



!#######################################################################
! <FUNCTION NAME="on_comp_domain">
!
! <DESCRIPTION>
! Determine if the point is in comp-domain for the processor
! </DESCRIPTION>
!
! <IN NAME="nxl" TYPE="integer">
! Integer labeling the particular xlandinsert pair
! </IN>
function on_comp_domain(ntable)

  integer, intent(in) :: ntable
  logical             :: on_comp_domain

  if(isc+Dom%ioff <= itable(ntable) .and. itable(ntable) <= iec+Dom%ioff .and. &
     jsc+Dom%joff <= jtable(ntable) .and. jtable(ntable) <= jec+Dom%joff) then
     on_comp_domain = .true.
  else 
     on_comp_domain = .false.  
  endif

end function on_comp_domain
! </FUNCTION> NAME="on_comp_domain"




end module ocean_vert_mix_mod



