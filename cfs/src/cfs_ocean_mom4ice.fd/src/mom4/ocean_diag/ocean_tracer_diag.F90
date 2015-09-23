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
module ocean_tracer_diag_mod
!  
!<CONTACT EMAIL="Stephen.Griffies@noaa.gov"> S.M. Griffies 
!</CONTACT>
!
!<CONTACT EMAIL="Ronald.Pacanowski@noaa.gov"> R.C. Pacanowski
!</CONTACT>
!
!<OVERVIEW>
! Routines for tracer diagnostics 
!</OVERVIEW>
!
!<DESCRIPTION>
! Routines for tracer diagnostics.  Some are printed to ascii output, some are sent 
! to diagnostic manager. 
! </DESCRIPTION>
!
!<NAMELIST NAME="ocean_tracer_diag_nml">
!  <DATA NAME="tracer_conserve_days" UNITS="days" TYPE="real">
!  Number of days between which compute the tracer conservation diagnostics. 
!  </DATA> 
!  <DATA NAME="diag_freq" UNITS="dimensionless" TYPE="integer">
!  Number of time steps between which compute the diagnostics.
!  </DATA> 
!  <DATA NAME="diagnose_mixing_simple" TYPE="logical">
!  Set true for using diagnose_mixing_simple. 
!  </DATA> 
!  <DATA NAME="diagnose_mixing_sort" TYPE="logical">
!  Set true for using diagnose_mixing_sort. 
!  </DATA> 
!  <DATA NAME="debug_diagnose_mixingA" TYPE="logical">
!  Set true for help with debugging the diagnostic for mixing.
!  </DATA> 
!  <DATA NAME="debug_diagnose_mixingB" TYPE="logical">
!  Set true for more help with debugging the diagnostic for mixing.
!  </DATA> 
!  <DATA NAME="smooth_kappa_sort" TYPE="logical">
!  Set true to smooth the diagnosed mixing from the sorted approach.
!  </DATA> 
!  <DATA NAME="rho_grad_min" UNITS="kg/m^3/m" TYPE="real">
!  min vertical density gradient (kg/m^3/m) used in computing kappa sorted
!  in the diagnostic mixing sorted. 
!  </DATA> 
!  <DATA NAME="rho_grad_max" UNITS="kg/m^3/m" TYPE="real">
!  max vertical density gradient (kg/m^3/m) used in computing kappa sorted
!  </DATA>
!  <DATA NAME="buoyancy_crit" UNITS="m^2/sec" TYPE="real">
!  Critical buoyancy difference relative to surface for computing mixed layer depth. 
!  </DATA>
!  <DATA NAME="psu2ppt" TYPE="real">
!  The realistic EOS used in mom4 requires salinity to 
!  use the Practical Salinity Scale (pss).  This scale is 
!  also known as the Practical Salinity Unit (psu).
!  Additionally, ocean measurements use the psu scale
!  Hence, mom4 interprets its salinity as psu.  
!
!  However, salinity as an absolute concentration in 
!  parts per thousand is more convenient to use when 
!  performing budget analyses such as in this module.  
!  Conversion between pss and ppt depends on the precise
!  ratio of ions in the seawater. Hence, the conversion
!  is not constant. However, it is close to a constant,
!  as reported in Jackett etal (2004).  For purposes of 
!  budgets, we take this conversion as a constant.  
!  The conversion is 
!
!  s(ppt) = psu2ppt * s(psu) 
!
!  where again s(psu) is what mom4 carries as its 
!  prognostic salinity field. 
!
!  Jackett etal (2004), correcting a type in equation (53) 
!  of Feistel (2003), report that 
!
!  s(ppt) = 1.004867 * s(psu)
!
!  </DATA>
!  <DATA NAME="do_bitwise_exact_sum" TYPE="logical">
!  Set true to do bitwise exact global sum. When it is false, the global
!  sum will be non-bitwise_exact, but will significantly increase efficiency. 
!  The default value is false.
!  </DATA> 
!</NAMELIST>
!
use constants_mod,         only: grav, rho0, epsln
use diag_manager_mod,      only: register_diag_field, send_data, need_data
use fms_mod,               only: open_namelist_file, check_nml_error, close_file
use fms_mod,               only: write_version_number, FATAL, stdout, stdlog
use mpp_domains_mod,       only: mpp_global_field, mpp_global_sum, BITWISE_EXACT_SUM
use mpp_domains_mod,       only: NON_BITWISE_EXACT_SUM
use mpp_mod,               only: mpp_error, mpp_max
use mpp_mod,               only: mpp_clock_id, mpp_clock_begin, mpp_clock_end, CLOCK_ROUTINE
use time_manager_mod,      only: time_type, increment_time

use ocean_density_mod,     only: density_delta_sfc
use ocean_domains_mod,     only: get_local_indices
use ocean_obc_mod,         only: ocean_obc_tracer_flux, ocean_obc_mass_flux
use ocean_tracer_util_mod, only: tracer_min_max
use ocean_types_mod,       only: ocean_domain_type, ocean_grid_type
use ocean_types_mod,       only: ocean_prog_tracer_type, ocean_diag_tracer_type
use ocean_types_mod,       only: ocean_external_mode_type, ocean_density_type
use ocean_types_mod,       only: ocean_time_type, ocean_time_steps_type, ocean_thickness_type
use ocean_types_mod,       only: TWO_LEVEL, THREE_LEVEL
use ocean_util_mod,        only: write_timestamp
use ocean_workspace_mod,   only: wrk1, wrk2 

implicit none

private

#include <ocean_memory.h>

real    :: dtts
real    :: dteta
real    :: dtime
real    :: dtimer

integer :: num_prog_tracers
integer :: num_diag_tracers
logical :: have_obc 

integer :: index_temp
integer :: index_salt
integer :: index_frazil

! for diagnostics clocks 
integer :: id_mixed_layer_depth
integer :: id_potrho_mixed_layer
integer :: id_diagnose_depth_of_potrho
integer :: id_diagnose_depth_of_theta
integer :: id_tracer_numerical
integer :: id_mixing
integer :: id_conservation
integer :: id_total_tracer
integer :: id_tracer_integrals
integer :: id_tracer_change
integer :: id_tracer_land_cell_check

type(ocean_grid_type), pointer   :: Grd =>NULL()
type(ocean_domain_type), pointer :: Dom =>NULL()

logical :: module_is_initialized = .FALSE.
integer :: tendency=0

character(len=128) :: version=&
     '$Id$'
character (len=128) :: tagname = &
     '$Name$'

! for tracer conservation: need to know num_prog_tracers to allocate
real, dimension(:,:,:), allocatable :: tracer_flux
real, dimension(:,:,:), allocatable :: tracer_source
real, dimension(:,:,:), allocatable :: tracer_smooth
real, dimension(:,:)  , allocatable :: material_flux
real, dimension(:,:)  , allocatable :: eta_smooth
real, dimension(:,:)  , allocatable :: area_t
real, dimension(:), allocatable     :: tracer_start 
real, dimension(:), allocatable     :: tracer_final 
real                                :: mass_start=0.0
real                                :: mass_final=0.0
real                                :: volume_start=0.0
real                                :: volume_final=0.0
real                                :: tracer_conserve_days = 30.0 ! number of days over which compute tracer conservation
integer                             :: itts_tracer, itte_tracer    ! itt starting and ending tracer conservation diagnostic
integer                             :: itts_mass, itte_mass        ! itt starting and ending mass conservation diagnostic
integer                             :: itts_volume, itte_volume    ! itt starting and ending volume conservation diagnostic

! for saving out the total tracer content and mass/volume  
integer, allocatable, dimension(:) :: id_ttracer 
integer :: id_volume_seawater

! for diagnosing mixing between temperature classes. 
! assumes flat ocean bottom and single computer processor 
logical :: diagnose_mixing_simple=.false.
logical :: diagnose_mixing_sort  =.false.
logical :: debug_diagnose_mixingA=.false.
logical :: debug_diagnose_mixingB=.false.
real    :: rho_grad_max=1e28      !max vertical density gradient (kg/m^3/m) used in computing kappa
real    :: rho_grad_min=1e-5      !min vertical density gradient (kg/m^3/m) used in computing kappa 
real    :: alpha=0.255            !from linear equation of state 
integer :: smooth_kappa_sort=0    !for smoothing the diagnosed diffusivity 
integer :: id_kappa_simple=-1
integer :: id_kappa_sort=-1

integer :: global_sum_flag      ! flag for mpp_global_sum

! frequency for computing certain diagnostics 
integer :: diag_freq = -1

! for diagnostic manager 
logical :: used

! for output
integer :: unit=6

! for mixed layer depth 
integer :: id_mld=-1

! for depth of isopycnal and potential temperature surfaces 
integer :: id_depth_of_potrho=-1
integer :: id_depth_of_theta =-1

! for potential density mixed layer 
integer :: id_potrho_mix_depth=-1
integer :: id_potrho_mix_base=-1
real    :: buoyancy_crit=0.0003 ! (m^2/sec) critical buoyancy difference relative to surface 
logical :: do_bitwise_exact_sum = .false. !

! for salt budgets
real :: psu2ppt=1.004867   ! conversion from mom4's salinity in psu to ppt concentration 

public ocean_tracer_diag_init
public ocean_tracer_diagnostics 

private total_tracer
private tracer_change
private tracer_integrals 
private tracer_land_cell_check
private tracer_conservation
private volume_conservation
private mixed_layer_depth 
private diagnose_kappa_simple
private diagnose_kappa_sort
private diagnose_depth_of_potrho
private diagnose_depth_of_theta
private potrho_mixed_layer 
private total_volume
private send_total_tracer

namelist /ocean_tracer_diag_nml/ tracer_conserve_days , diag_freq, psu2ppt,      &
                                 diagnose_mixing_simple, diagnose_mixing_sort,   &
                                 debug_diagnose_mixingA, debug_diagnose_mixingB, &
                                 smooth_kappa_sort, rho_grad_min, rho_grad_max,  &
                                 buoyancy_crit, do_bitwise_exact_sum 

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_tracer_diag_init">
!
! <DESCRIPTION>
! Initialize the ocean_tracer_diag module containing subroutines
! diagnosing tracer related properties of the simulation.  These are 
! not terms in the equations, but rather they are diagnosed from 
! terms. 
! </DESCRIPTION>
!
  subroutine ocean_tracer_diag_init(Grid, Domain, Time, Time_steps, T_prog, T_diag, Dens, obc)

    type(ocean_grid_type), target, intent(in)   :: Grid
    type(ocean_domain_type), target, intent(in) :: Domain
    type(ocean_time_type), intent(in)           :: Time
    type(ocean_time_steps_type), intent(in)     :: Time_steps 
    type(ocean_prog_tracer_type), intent(in)    :: T_prog(:)
    type(ocean_diag_tracer_type), intent(in)    :: T_diag(:)
    type(ocean_density_type), intent(in)        :: Dens
    logical, intent(in)                         :: obc

    integer           :: n, ioun, io_status, ierr

    if (module_is_initialized) return

    module_is_initialized = .TRUE.

    num_prog_tracers = size(T_prog,1)
    num_diag_tracers = size(T_diag,1)

    call write_version_number(version, tagname)

    ioun = open_namelist_file()
    read(ioun, ocean_tracer_diag_nml, iostat=io_status)
    write (stdlog(), ocean_tracer_diag_nml)
    write (stdout(),'(/)')
    write (stdout(), ocean_tracer_diag_nml)
    ierr = check_nml_error(io_status,'ocean_tracer_diag_nml')
    call close_file(ioun)

    if (diag_freq == 0) diag_freq = 1

    if(do_bitwise_exact_sum) then
       global_sum_flag = BITWISE_EXACT_SUM
    else
       global_sum_flag = NON_BITWISE_EXACT_SUM
    endif

#ifndef STATIC_MEMORY
    call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
    ni = Grid%ni
    nj = Grid%nj
    nk = Grid%nk
#endif

    Dom => Domain
    Grd => Grid

    dtts     = Time_steps%dtts
    dteta    = Time_steps%dteta
    dtime    = Time_steps%dtime_t
    tendency = Time_steps%tendency
    dtimer   = 1.0/(dtime+epsln)
    have_obc = obc

    do n=1,num_prog_tracers
       if(T_prog(n)%name == 'temp' ) index_temp = n
       if(T_prog(n)%name == 'salt' ) index_salt = n
    enddo

    index_frazil=-1
    do n=1,num_diag_tracers
       if(T_diag(n)%name == 'frazil' ) index_frazil = n
    enddo

    allocate(area_t(isd:ied,jsd:jed) )    
    area_t(:,:) = Grd%dat(:,:)*Grd%tmask(:,:,1)
    if(have_obc) area_t(:,:) = area_t(:,:)*Grd%obc_tmask(:,:)

    if(tendency==THREE_LEVEL) then
        write (stdout(),'(/a)') 'Note: Will perform tracer and volume conservation tests based on time_tendency==threelevel.'
    elseif(tendency==TWO_LEVEL) then
        write (stdout(),'(/a)') 'Note: Will perform tracer and volume conservation tests based on time_tendency==twolevel.'
    endif

    if(diagnose_mixing_sort) then
        write (stdout(),'(/a)') 'Note: Will perform sorting to diagnose the mixing between temperature classes.'
        write (stdout(),'(a)') '      This algorithm is relevant only for cases with a linear equation of state.'
        write (stdout(),'(a)') '      It has been coded only for a single processor.' 
    endif
    if(diagnose_mixing_sort) then
        write(stdout(),'(/a)')'==> Warning: To use diagnose_mixing_sort, it is necessary to compile'
        write(stdout(),'(a)')'              with NAG routines in subroutine diagnose_kappa_sort.'
        write(stdout(),'(a)')'              Presently, those routines are commented out.'
        write(stdout(),'(a/)')'             Remove this error stop when have the correct routines linked.'
        call mpp_error(FATAL, '==>Error from ocean_tracer_diag: diagnose_mixing_sort needs tuning. See code for routines to link.')
    endif 

    if(diagnose_mixing_simple) then
        write (stdout(),'(/a)') 'Note: Will perform horz averaging to diagnose the mixing between temperature classes.'
        write (stdout(),'(a)') '      This algorithm is relevant only for cases with a linear equation of state.'
        write (stdout(),'(a)') '      It has been coded only for a single processor.' 
    endif

! registers to diagnostic manager 

    id_volume_seawater = register_diag_field ('ocean_model', 'volume_seawater', Time%model_time, &
         'total volume of seawater', 'm^3', missing_value=-1e1, range=(/-1e1,1e25/))

    allocate(id_ttracer(num_prog_tracers))
    id_ttracer=-1
    do n=1,num_prog_tracers
       if(T_prog(n)%name == 'temp') then 
           id_ttracer(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_total', Time%model_time, &
                'volume integrated '//trim(T_prog(n)%name), 'J/1e25', missing_value=-1e2, range=(/-1e2,1e2/))
       elseif(T_prog(n)%name(1:3) =='age') then 
           id_ttracer(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_total', Time%model_time, &
                'volume integrated '//trim(T_prog(n)%name), 'yr', missing_value=-1e2, range=(/-1e2,1e10/))
       else 
           id_ttracer(n) = register_diag_field ('ocean_model', trim(T_prog(n)%name)//'_total', Time%model_time, &
                'volume integrated '//trim(T_prog(n)%name), 'kg/1e18', missing_value=-1e6, range=(/-1e6,1e6/))
       endif
    enddo

    id_kappa_sort = register_diag_field ('ocean_model','kappa_sort', &
         Grd%tracer_axes(3:3), Time%model_time, &
         'kappa from sorting', 'm^2/sec', &
         missing_value=-1e6, range=(/-1e6,1e6/))
    id_kappa_simple = register_diag_field ('ocean_model','kappa_simple', &
         Grd%tracer_axes(3:3), Time%model_time, &
         'kappa from horz avg', 'm^2/sec', &
         missing_value=-1e6, range=(/-1e6,1e6/))

    id_mld = register_diag_field ('ocean_model', 'mld', &
         Grd%tracer_axes(1:2), Time%model_time, &
         'mixed layer depth ', 'm', &
         missing_value=0.0, range=(/0.0,1.e6/))   

    id_depth_of_potrho = register_diag_field ('ocean_model', 'depth_of_potrho', &
         Dens%potrho_axes(1:3), Time%model_time,    &
         'depth of potential density surface', 'm', &
         missing_value=-10.0, range=(/0.0,1.e10/))   

    id_depth_of_theta = register_diag_field ('ocean_model', 'depth_of_theta', &
         Dens%theta_axes(1:3), Time%model_time,  &
         'depth of potential temp surface', 'm', &
         missing_value=-10.0, range=(/0.0,1.e10/))   

    id_potrho_mix_depth = register_diag_field ('ocean_model','potrho_mix_depth',  &
         Grd%tracer_axes(1:2), Time%model_time, &
         'Depth of potential density mixed layer','m', &
         missing_value=-1.e6, range=(/-1e6,1e6/))

    id_potrho_mix_base = register_diag_field ('ocean_model','potrho_mix_base',  &
         Grd%tracer_axes(1:2), Time%model_time, &
         'Potential density at mixed layer base','kg/m^3', &
         missing_value=-1.e6, range=(/-1e6,1e6/))

    ! define clock integers 
    id_mixed_layer_depth        = mpp_clock_id('(Ocean tracer_diag: mld)'          ,grain=CLOCK_ROUTINE)
    id_potrho_mixed_layer       = mpp_clock_id('(Ocean tracer_diag: rho_mld)'      ,grain=CLOCK_ROUTINE)
    id_diagnose_depth_of_potrho = mpp_clock_id('(Ocean tracer_diag: potrho depth)' ,grain=CLOCK_ROUTINE)
    id_diagnose_depth_of_theta  = mpp_clock_id('(Ocean tracer_diag: theta depth)'  ,grain=CLOCK_ROUTINE)
    id_tracer_numerical         = mpp_clock_id('(Ocean tracer_diag: numerical)'    ,grain=CLOCK_ROUTINE)
    id_mixing                   = mpp_clock_id('(Ocean tracer_diag: mixing)'       ,grain=CLOCK_ROUTINE)
    id_conservation             = mpp_clock_id('(Ocean tracer_diag: conserve)'     ,grain=CLOCK_ROUTINE)
    id_total_tracer             = mpp_clock_id('(Ocean tracer_diag: total tracer)' ,grain=CLOCK_ROUTINE)
    id_tracer_integrals         = mpp_clock_id('(Ocean tracer_diag: integrals)'    ,grain=CLOCK_ROUTINE)
    id_tracer_change            = mpp_clock_id('(Ocean tracer_diag: change)'       ,grain=CLOCK_ROUTINE)
    id_tracer_land_cell_check   = mpp_clock_id('(Ocean tracer_diag: land check)'   ,grain=CLOCK_ROUTINE)

    end subroutine ocean_tracer_diag_init
! </SUBROUTINE>  NAME="ocean_tracer_diag_init"


!#######################################################################
! <SUBROUTINE NAME="ocean_tracer_diagnostics">
!
! <DESCRIPTION>
! Call diagnostics related to the tracer fields. 
! </DESCRIPTION>

subroutine ocean_tracer_diagnostics(Time, Thickness, T_prog, T_diag, Dens, Ext_mode, pme, river)

  type(ocean_time_type), intent(in)            :: Time
  type(ocean_thickness_type), intent(in)       :: Thickness
  type(ocean_prog_tracer_type), intent(in)     :: T_prog(:)
  type(ocean_diag_tracer_type), intent(in)     :: T_diag(:)
  type(ocean_density_type), intent(in)         :: Dens
  type(ocean_external_mode_type), intent(in)   :: Ext_mode
  real, intent(in), dimension(isd:ied,jsd:jed) :: pme
  real, intent(in), dimension(isd:ied,jsd:jed) :: river

  type(time_type)                              :: next_time 

  integer  :: tau, taum1, taup1
  integer  :: n

  if(size(T_prog,1) /= num_prog_tracers) then 
    call mpp_error(FATAL,'==>Error from ocean_tracer_diag_mod (ocean_tracer_diagnostics): T_prog has wrong dimensions')
  endif 

  tau   = Time%tau
  taum1 = Time%taum1
  taup1 = Time%taup1

  next_time = increment_time(Time%model_time, int(dtts), 0)
  
  ! compute mixed layer depth and send to diagnostics manager

  call mpp_clock_begin(id_mixed_layer_depth)
  if(need_data(id_mld, next_time)) then 
    call mixed_layer_depth(T_prog(index_salt)%field(isd:ied,jsd:jed,:,tau),T_prog(index_temp)%field(isd:ied,jsd:jed,:,tau),&
                           Dens%rho(isd:ied,jsd:jed,:), Dens%pressure_at_depth(isd:ied,jsd:jed,:), Time%model_time)
  endif
  call mpp_clock_end(id_mixed_layer_depth)


  call mpp_clock_begin(id_potrho_mixed_layer)
  call potrho_mixed_layer(Time, Thickness, Dens)
  call mpp_clock_end(id_potrho_mixed_layer)


 ! compute depth of isopycnal surfaces and send to diagnostics manager
  call mpp_clock_begin(id_diagnose_depth_of_potrho)
  if(need_data(id_depth_of_potrho, next_time)) then 
    call diagnose_depth_of_potrho(Time, Dens)
  endif
  call mpp_clock_end(id_diagnose_depth_of_potrho)

 ! compute depth of potential temp surfaces and send to diagnostics manager
  call mpp_clock_begin(id_diagnose_depth_of_theta)
  if(need_data(id_depth_of_theta, next_time)) then 
    call diagnose_depth_of_theta(Time, Dens, T_prog)
  endif
  call mpp_clock_end(id_diagnose_depth_of_theta)

  call mpp_clock_begin(id_tracer_numerical)
  if (diag_freq > 0) then
      if (mod(Time%itt, diag_freq) == 0) then
          do n = 1,num_prog_tracers
             call tracer_min_max(Time, T_prog(n))
             call tracer_integrals(Time, Thickness, T_prog(n), Dens, Ext_mode)         
          enddo
          call tracer_change(Time, Thickness, T_prog, T_diag, Ext_mode, pme, river)
          call tracer_land_cell_check (Time, T_prog)
      endif
  endif
  call mpp_clock_end(id_tracer_numerical)


  ! calls for diagnosing effective mixing coefficients 

  call mpp_clock_begin(id_mixing)
  if(diagnose_mixing_simple) call diagnose_kappa_simple(Time, T_prog(index_temp))
  if(diagnose_mixing_sort)   call diagnose_kappa_sort(Time, T_prog(index_temp))
  call mpp_clock_end(id_mixing)
  
  call mpp_clock_begin(id_conservation)
  if (nint(dtts) /= 0) then 
    call volume_conservation (Time, Thickness, Ext_mode, pme, river)
    call tracer_conservation (Time, Thickness, T_prog, T_diag, Ext_mode, Dens, pme, river)
  endif 
  call mpp_clock_end(id_conservation)

  ! send total tracer to diag_manager 
  call mpp_clock_begin(id_total_tracer)
  do n=1,num_prog_tracers 
    if(need_data(id_ttracer(n), next_time)) then 
      call send_total_tracer(Time, Thickness, T_prog(n), n)
    endif 
  enddo
  call mpp_clock_end(id_total_tracer)

end subroutine ocean_tracer_diagnostics
! </SUBROUTINE>  NAME="ocean_tracer_diagnostics"


!#######################################################################
! <SUBROUTINE NAME="mixed_layer_depth">
!
! <DESCRIPTION>
! Diagnose mixed layer depth (m), which is defined as the depth ( > 0 )
! where the buoyancy difference with respect to the surface level is
! equal to sfchmxl (=0.0003) (m/s2). Use dbsfc, which is defined on zt levels.
! Note that the diagnosed mixed layer depth is only used for diagnostics.            
! </DESCRIPTION>
!
subroutine mixed_layer_depth(salinity, theta, rho, pressure, model_time)

  real, intent(in), dimension(isd:ied,jsd:jed,nk) :: salinity
  real, intent(in), dimension(isd:ied,jsd:jed,nk) :: theta
  real, intent(in), dimension(isd:ied,jsd:jed,nk) :: rho
  real, intent(in), dimension(isd:ied,jsd:jed,nk) :: pressure 
  type(time_type), intent(in)                     :: model_time

  real, dimension(isd:ied,jsd:jed) :: hmxl
  real, parameter :: epsln   = 1.0e-20  ! for divisions 
  real, parameter :: sfchmxl = 0.0003   ! (m/s2) critical value that is used to compute mixed layer depth
  integer         :: i, j, k, km1

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL,'==>Error from ocean_tracer_diag_mod (mixed_layer_depth): module needs initialization ')
  endif 

  wrk2(:,:,1) = 0.0
  wrk1(:,:,:) = density_delta_sfc( rho(:,:,:), salinity(:,:,:), theta(:,:,:), pressure(:,:,:))
  do k=2,nk
     wrk2(:,:,k) = -grav*wrk1(:,:,k-1)/rho(:,:,k)
  enddo

  do j=jsc,jec
     do i=isc,iec
        if (Grd%kmt(i,j) == 0) then
            hmxl(i,j) = 0.0
        else
            hmxl(i,j) = Grd%zt(Grd%kmt(i,j))
        end if
     enddo
  enddo
  do k=2,nk
     km1 = k-1  
     do j=jsc,jec
        do i=isc,iec
           if (Grd%kmt(i,j) == 0) then
               hmxl(i,j) = 0.0
           else
               if ( wrk2(i,j,k) >= sfchmxl .and. hmxl(i,j) == Grd%zt(Grd%kmt(i,j)) ) then
                   hmxl(i,j) = Grd%zt(km1) - (Grd%zt(km1) - Grd%zt(k))  &
                        * (sfchmxl-wrk2(i,j,km1))        &
                        / (wrk2(i,j,k) - wrk2(i,j,km1) + epsln)
               endif
               hmxl(i,j) = hmxl(i,j) * Grd%tmask(i,j,1)
           endif
        enddo
     enddo
  enddo
  used = send_data (id_mld, hmxl(:,:), &
         model_time, rmask=Grd%tmask(:,:,1), &
         is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

end subroutine mixed_layer_depth
! </SUBROUTINE>  NAME="mixed_layer_depth"


!#######################################################################
! <SUBROUTINE NAME="tracer_change">
!
! <DESCRIPTION>
!
! Compute change in tracer over a time step.  Compate global change
! in tracer to global integral of boundary forcing.  Compute error
! due to any mis-match.  This routine is very useful for detecting
! bugs in tracer routines.  
! 
! </DESCRIPTION>
subroutine tracer_change (Time, Thickness, T_prog, T_diag, Ext_mode, pme, river)

  type(ocean_time_type), intent(in)            :: Time
  type(ocean_thickness_type), intent(in)       :: Thickness
  type(ocean_prog_tracer_type), intent(in)     :: T_prog(:) 
  type(ocean_diag_tracer_type), intent(in)     :: T_diag(:) 
  type(ocean_external_mode_type), intent(in)   :: Ext_mode
  real, intent(in), dimension(isd:ied,jsd:jed) :: pme
  real, intent(in), dimension(isd:ied,jsd:jed) :: river

  real, dimension(isd:ied,jsd:jed) :: tmp
  real, dimension(isd:ied,jsd:jed) :: area_k
  real, dimension(isd:ied,jsd:jed) :: volume_t
  real, dimension(isd:ied,jsd:jed) :: thick_t
  real :: volume_error, tracer_error, flux_error
  real :: eta_taum1, eta_tau, eta_taup1
  real :: tracer_input_stf, tracer_input_btf, tracer_input_pme
  real :: tracer_input_otf, tracer_input_smooth
  real :: tracer_input_river, tracer_input_total
  real :: pme_input, river_input, obc_input
  real :: dvoltracer(0:nk), eta_t_avg
  real :: frazil_heat 
  real :: tracer_sources(0:nk) 
  real :: ttracer 
  real :: volume_total, volume_source, volume_change, volume_surface 
  real :: volume_smooth
  real :: conversion, convert 

  integer :: k, n
  integer :: taum1, tau, taup1

  call mpp_clock_begin(id_tracer_change)

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL,'==>Error from ocean_tracer_diag_mod (tracer_change): module needs initialization ')
  endif 

  if(size(T_prog,1) /= num_prog_tracers) then 
    call mpp_error(FATAL,'==>Error from ocean_tracer_diag_mod (tracer_change): passing T_prog with wrong dimensions')
  endif 

  write (stdout(),'(//50x,a/)') ' Volume, mass, and tracer change summary (over one timestep):'

  ! compute total volume, total mass, and their changes over a leap-frog time step 
  pme_input       = 0.0 ! total mass of evap-precip input to ocean over leap-frog time step (kg)
  river_input     = 0.0 ! total mass of river water input to ocean over leap-frog time step (kg)
  obc_input       = 0.0 ! total mass of water input through open boundaries to ocean over leap-frog time step (kg)
  volume_surface  = 0.0 ! volume of ocean within the k=1 model grid cells (m^3)
  eta_t_avg       = 0.0 ! area average of eta_t (m)
  volume_total    = 0.0 ! volume of ocean on tracer cells (m^3)
  volume_smooth   = 0.0 ! from eta_t surface filtering (m^3)
  volume_change   = 0.0 ! change in volume from tau-1 to tau+1 (m^3)
   
  taum1 = Time%taum1
  tau   = Time%tau
  taup1 = Time%taup1
 
  ! surface level terms  
  k=1

  tmp(:,:)       = area_t(:,:)*Thickness%dht(:,:,k,tau)
  volume_surface = mpp_global_sum(Dom%domain2d,tmp(:,:), global_sum_flag)
 
  tmp(:,:)  = dtime*pme(:,:)*area_t(:,:)
  pme_input = mpp_global_sum(Dom%domain2d,tmp(:,:), global_sum_flag)

  tmp(:,:)    = dtime*river(:,:)*area_t(:,:)
  river_input = mpp_global_sum(Dom%domain2d,tmp(:,:), global_sum_flag)

  tmp(:,:)       = area_t(:,:)*(Ext_mode%eta_t(:,:,taup1)-Ext_mode%eta_t(:,:,taum1))
  volume_change  = mpp_global_sum(Dom%domain2d,tmp(:,:), global_sum_flag)

  tmp(:,:)   = area_t(:,:)*Ext_mode%eta_t(:,:,tau)
  eta_t_avg  = mpp_global_sum(Dom%domain2d,tmp(:,:), global_sum_flag)/Grd%tcellsurf

  tmp(:,:)      = dtime*area_t(:,:)*Ext_mode%eta_source(:,:)
  volume_source = mpp_global_sum(Dom%domain2d,tmp(:,:), global_sum_flag)

  tmp(:,:)      = dtime*area_t(:,:)*Ext_mode%surface_smooth(:,:)
  volume_smooth = mpp_global_sum(Dom%domain2d,tmp(:,:), global_sum_flag)


  ! compute volume passing across open boundaries 
  if (have_obc) then 
    call ocean_obc_mass_flux(Time, Ext_mode, tmp)
    tmp(:,:)  = tmp(:,:)*dtime
    obc_input = mpp_global_sum(Dom%domain2d,tmp(:,:), global_sum_flag) 
  endif 

  volume_error = volume_change-pme_input-river_input-obc_input-volume_source-volume_smooth

  do k=1,nk
    volume_t(:,:) = Thickness%dht(:,:,k,tau)*Grd%dat(:,:)*Grd%tmask(:,:,k)
    if(have_obc) volume_t(:,:) = volume_t(:,:)*Grd%obc_tmask(:,:)
    volume_total  = volume_total + mpp_global_sum(Dom%domain2d,volume_t(:,:), global_sum_flag)
  enddo
  
  write (stdout(),'(a,es24.17,a)') ' Area average surface height on tracer cells (tau)    = ',eta_t_avg,' m'
  write (stdout(),'(a,es24.17,a)') ' Surface area of k=1 tracer cells                     = ',Grd%tcellsurf,' m^2'
  write (stdout(),'(a,es24.17,a)') ' Area integral of surface height on tracer cells      = ',eta_t_avg*Grd%tcellsurf,' m^3'
  write (stdout(),'(a,es24.17,a)') ' Volume of k=1 tracer cells (tau)                     = ',volume_surface,' m^3'
  write (stdout(),'(a,es24.17,a)') ' Volume change from taum1 to taup1                    = ',volume_change,' m^3'
  write (stdout(),'(a,es24.17,a)') ' Volume of ocean tracer cells (tau)                   = ',volume_total,' m^3'
  
  write (stdout(),'(a,es24.17,a)') ' Volume from sources from taum1 to taup1              = ',volume_source,' m^3'
  write (stdout(),'(a,es24.17,a)') ' Volume from surface smoothing from taum1 to taup1    = ',volume_smooth,' m^3'
  write (stdout(),'(a,es24.17,a)') ' Volume of precip-evap input from taum1 to taup1      = ',pme_input,' m^3'
  write (stdout(),'(a,es24.17,a)') ' Volume of river water input from taum1 to taup1      = ',river_input,' m^3'
  write (stdout(),'(a,es24.17,a)') ' Volume flux through open bnds from taum1 to taup1    = ',obc_input,' m^3'
  write (stdout(),'(a,es24.17,a)') ' Mismatch between volume change and fresh water input = ',volume_error,' m^3'


  ! compute total tracer and changes over a time step 
  do n=1,num_prog_tracers

     dvoltracer         = 0.0 ! volume(taup1)*tracer(taup1)-volume(taum1)*tracer(taum1))
     tracer_input_stf   = 0.0 ! tracer quantity input to ocean via stf (positive is tracer entering ocean)
     tracer_input_btf   = 0.0 ! tracer quantity input to ocean via btf (positive is tracer leaving ocean)
     tracer_input_otf   = 0.0 ! tracer quantity input through open boundaries (positive is tracer entering ocean) 
     tracer_input_pme   = 0.0 ! tracer quantity input to ocean via pme flux 
     tracer_input_river = 0.0 ! tracer quantity input to ocean via rivers
     tracer_input_total = 0.0 ! total tracer quantity input to ocean 
     frazil_heat        = 0.0 ! heat (Joule) from frazil formation  
     tracer_input_smooth= 0.0 ! tracer input due to smoothing the surface height eta_t 
     tracer_sources     = 0.0 ! contribution from thickness weighted sources, including 
                              ! neutral diffusion, sigma diffusion, overflow, xlandmix, 
                              ! sponges, rivers, and other more traditional sources.  

     conversion = T_prog(n)%conversion     
     k=1

     tmp(:,:) = area_t(:,:)*T_prog(n)%conversion* &
          (Thickness%dht(:,:,k,taup1)*T_prog(n)%field(:,:,k,taup1) &
          -Thickness%dht(:,:,k,taum1)*T_prog(n)%field(:,:,k,taum1))  
     dvoltracer(k) = mpp_global_sum(Dom%domain2d,tmp(:,:), global_sum_flag) 

     tmp(:,:) = area_t(:,:)*conversion*dtime*T_prog(n)%stf(:,:)
     tracer_input_stf = mpp_global_sum(Dom%domain2d,tmp(:,:), global_sum_flag)

     tmp(:,:) = area_t(:,:)*conversion*dtime*T_prog(n)%btf(:,:)
     tracer_input_btf = -mpp_global_sum(Dom%domain2d,tmp(:,:), global_sum_flag)

     tmp(:,:) = area_t(:,:)*conversion*dtime*pme(:,:)*T_prog(n)%tpme(:,:)
     tracer_input_pme = mpp_global_sum(Dom%domain2d,tmp(:,:), global_sum_flag) 

     tmp(:,:) = area_t(:,:)*conversion*dtime*river(:,:)*T_prog(n)%triver(:,:)
     tracer_input_river = mpp_global_sum(Dom%domain2d,tmp(:,:), global_sum_flag)

     tmp(:,:) = area_t(:,:)*conversion*dtime*T_prog(n)%surface_smooth(:,:)
     tracer_input_smooth = mpp_global_sum(Dom%domain2d,tmp(:,:), global_sum_flag)
     
     tmp(:,:) = area_t(:,:)*conversion*dtime*T_prog(n)%th_tendency(:,:,k)
     if(have_obc) tmp(:,:) = tmp(:,:)*Grd%obc_tmask(:,:)
     tracer_sources(k) = mpp_global_sum(Dom%domain2d,tmp(:,:), global_sum_flag)

     ! for cleaner diagnostics, remove river, pme, and smoother contributions,
     ! each of which are added to T_prog%th_tendency inside of ocean_tracer.F90
     tracer_sources(k) = tracer_sources(k) - tracer_input_river - tracer_input_pme - tracer_input_smooth

     ! ocean_obc_tracer_flux returns vertically integrated horizontal tracer flux     
     ! at last internal points and zero otherwise     
     if(have_obc) then 
       call ocean_obc_tracer_flux(Time, T_prog(n), tmp, n, send_out =.true.) 
       tmp(:,:) = tmp(:,:)*conversion*dtime
       tracer_input_otf = mpp_global_sum(Dom%domain2d,tmp(:,:), global_sum_flag)
     endif 

     ! note that river contribution is added via T_prog%th_tendency inside ocean_rivermix_mod
     tracer_input_total = tracer_input_stf   + tracer_input_btf + tracer_input_otf + &
                          tracer_input_river + tracer_input_pme + tracer_input_smooth

     frazil_heat = 0.0
     if(T_prog(n)%name =='temp' .and. index_frazil > 0) then
       tmp(:,:)    = T_diag(index_frazil)%field(:,:,1)*area_t(:,:)/T_diag(index_frazil)%factor 
       frazil_heat = mpp_global_sum(Dom%domain2d,tmp(:,:), global_sum_flag) 
       tracer_input_total = tracer_input_total + frazil_heat
     endif 

     do k=2,nk

        area_k(:,:) = Grd%dat(:,:)*Grd%tmask(:,:,k)  
        if(have_obc) area_k(:,:) = area_k(:,:)*Grd%obc_tmask(:,:)

        tmp(:,:) = area_k(:,:)*T_prog(n)%conversion* &
             (Thickness%dht(:,:,k,taup1)*T_prog(n)%field(:,:,k,taup1) &
             -Thickness%dht(:,:,k,taum1)*T_prog(n)%field(:,:,k,taum1))  
        dvoltracer(k) = mpp_global_sum(Dom%domain2d,tmp(:,:), global_sum_flag) 

        tmp(:,:) = area_k(:,:)*conversion*dtime*T_prog(n)%th_tendency(:,:,k)
        if(have_obc) tmp(:,:) = tmp(:,:)*Grd%obc_tmask(:,:)
        tracer_sources(k) = mpp_global_sum(Dom%domain2d,tmp(:,:), global_sum_flag)
     enddo

     do k=1,nk ! sum integrals over all levels 
        tracer_sources(0) = tracer_sources(0) + tracer_sources(k)
        dvoltracer(0)     = dvoltracer(0)     + dvoltracer(k) 
     enddo

     ttracer      = total_tracer(T_prog(n),Thickness, tau)
     tracer_error = dvoltracer(0) -tracer_input_total-tracer_sources(0)
     flux_error   = tracer_error/(dtime*Grd%tcellsurf + epsln)

     if (T_prog(n)%name =='temp') then

         write (stdout(),'(/a)') ' ----Single time step diagnostics for tracer '//trim(T_prog(n)%name)//'----' 
         write (stdout(),'(a,es24.17,a)') ' Total heat in ocean at (tau) referenced to 0degC   = ',ttracer,' J'
         write (stdout(),'(a,es24.17,a)') ' Total heat input to ocean referenced to 0degC      = ',tracer_input_total,' J'
         write (stdout(),'(a,es24.17,a)') ' Heat input via surface tracer fluxes               = ',tracer_input_stf,' J'
         write (stdout(),'(a,es24.17,a)') ' Heat input via bottom tracer fluxes                = ',tracer_input_btf,' J'
         write (stdout(),'(a,es24.17,a)') ' Heat input via open boundaries                     = ',tracer_input_otf,' J'
         write (stdout(),'(a,es24.17,a)') ' Heat input via precip-evap                         = ',tracer_input_pme,' J'
         write (stdout(),'(a,es24.17,a)') ' Heat input via rivers and calving                  = ',tracer_input_river,' J'
         write (stdout(),'(a,es24.17,a)') ' Heat input via frazil formation                    = ',frazil_heat,' J'
         write (stdout(),'(a,es24.17,a)') ' Heat input via sources                             = ',tracer_sources(0),' J'
         write (stdout(),'(a,es24.17,a)') ' Heat input via eta_t smoother                      = ',tracer_input_smooth,' J'
         write (stdout(),'(a,es24.17,a)') ' cp*[T(taup1)*vol(taup1)-T(taum1)*vol(taum1))]      = ',dvoltracer(0),' J'
         write (stdout(),'(a,es24.17,a)') ' Tracer mismatch: cp*d(vol*T)-input                 = ',tracer_error,' J'
         write (stdout(),'(a,es24.17,a)') ' Mismatch converted to a surface flux               = ',flux_error,' W/m^2'

     elseif (T_prog(n)%name(1:3) =='age') then

         write (stdout(),'(/a)') ' ----Single time step diagnostics for tracer '//trim(T_prog(n)%name)//'----' 
         write (stdout(),'(a,es24.17,a)') ' Total ocean age at (tau)                            = ',&     
                                                                                         ttracer/volume_total,' yr'
         write (stdout(),'(a,es24.17,a)') ' Total age input to ocean                            = ',&
                                                                                         tracer_input_total/volume_total,' yr'
         write (stdout(),'(a,es24.17,a)') ' Age input via surface fluxes                        = ',&
                                                                                         tracer_input_stf/volume_total,' yr'
         write (stdout(),'(a,es24.17,a)') ' Age input via bottom fluxes                         = ',&
                                                                                         tracer_input_btf/volume_total,' yr'
         write (stdout(),'(a,es24.17,a)') ' Age input via open boundaries                       = ',&
                                                                                         tracer_input_otf/volume_total,' yr'
         write (stdout(),'(a,es24.17,a)') ' Age input via precip-evap                           = ',&
                                                                                         tracer_input_pme/volume_total,' yr'
         write (stdout(),'(a,es24.17,a)') ' Age input via rivers and calving                    = ',&
                                                                                         tracer_input_river/volume_total,' yr'
         write (stdout(),'(a,es24.17,a)') ' Age input via sources                               = ',&
                                                                                         tracer_sources(0)/volume_total,' yr'
         write (stdout(),'(a,es24.17,a)') ' Age input via eta_t smoother                        = ',&
                                                                                         tracer_input_smooth/volume_total,' yr'
         write (stdout(),'(a,es24.17,a)') ' [T(taup1)*vol(tau+1)-T(taum1)*vol(tau-1))]          = ',&
                                                                                         dvoltracer(0)/volume_total,' yr'
         write (stdout(),'(a,es24.17,a)') ' Tracer mismatch: d(vol*T)-input                     = ',&
                                                                                         tracer_error/volume_total,' yr'
         write (stdout(),'(a,es24.17,a)') ' Mismatch converted to a surface flux                = ',flux_error,' m'

     else 

         write (stdout(),'(/a)') '----Single time step diagnostics for tracer '//trim(T_prog(n)%name)//'----' 
         write (stdout(),'(a,es24.17,a)') ' Total tracer in ocean at time  (tau)                = ',ttracer,' kg'
         write (stdout(),'(a,es24.17,a)') ' Total tracer input to ocean                         = ',tracer_input_total,' kg'
         write (stdout(),'(a,es24.17,a)') ' Tracer input via surface fluxes                     = ',tracer_input_stf,' kg'
         write (stdout(),'(a,es24.17,a)') ' Tracer input via bottom fluxes                      = ',tracer_input_btf,' kg'
         write (stdout(),'(a,es24.17,a)') ' Tracer input via open boundaries                    = ',tracer_input_otf,' kg'
         write (stdout(),'(a,es24.17,a)') ' Tracer input via precip-evap                        = ',tracer_input_pme,' kg'
         write (stdout(),'(a,es24.17,a)') ' Tracer input via rivers and calving                 = ',tracer_input_river,' kg'
         write (stdout(),'(a,es24.17,a)') ' Tracer input via sources                            = ',tracer_sources(0),' kg'
         write (stdout(),'(a,es24.17,a)') ' Tracer input via eta_t smoother                     = ',tracer_input_smooth,' kg'

         if (T_prog(n)%name =='salt') then

         write (stdout(),'(a,es24.17,a)') ' .001*[S(taup1)*vol(taup1)-S(taum1)*vol(taum1))]     = ', &
                                                                                          psu2ppt*dvoltracer(0),' kg'
        write (stdout(),'(a,es24.17,a)')  ' Tracer mismatch: 0.001*d(vol*T)-input               = ', &
                                                                                          psu2ppt*tracer_error,' kg '
         write (stdout(),'(a,es24.17,a)') ' Mismatch converted to a surface flux                = ', &
                                                                                          psu2ppt*flux_error,' kg/(m^2 sec) '

         else  

         write (stdout(),'(a,es24.17,a)') ' [tracer(taup1)*vol(taup1)-tracer(taum1)*vol(taum1))]= ',dvoltracer(0),' kg'
         write (stdout(),'(a,es24.17,a)') ' Tracer mismatch: d(vol*T)-input                     = ',tracer_error,' kg'
         write (stdout(),'(a,es24.17,a)') ' Mismatch converted to a surface flux                = ',flux_error,' kg/(m^2 sec) '

         endif

     endif

  enddo  ! end n-loop

  call mpp_clock_end(id_tracer_change)

end subroutine tracer_change
! </SUBROUTINE>  NAME="tracer_change"

     
!#######################################################################
! <FUNCTION NAME="total_tracer">
!
! <DESCRIPTION>
! Compute integrated tracer in model. 
! </DESCRIPTION>
!
function total_tracer (Tracer, Thickness, index)

  type(ocean_prog_tracer_type), intent(in) :: Tracer
  type(ocean_thickness_type), intent(in)   :: Thickness
  integer, intent(in)                      :: index 
  real                                     :: total_tracer

  real, dimension(isd:ied,jsd:jed) :: tracer_k
  integer                          :: k

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL,'==>Error from ocean_tracer_diag_mod (total_tracer): module needs initialization ')
  endif 

  total_tracer=0.0
  do k=1,nk
    tracer_k(:,:) = Tracer%conversion*Grd%tmask(:,:,k)*Grd%dat(:,:)*Thickness%dht(:,:,k,index)*Tracer%field(:,:,k,index)
    if(have_obc) tracer_k(:,:) = tracer_k(:,:)*Grd%obc_tmask(:,:)
    total_tracer  = total_tracer + mpp_global_sum(Dom%domain2d,tracer_k(:,:), global_sum_flag)
  enddo

end function total_tracer
! </FUNCTION>  NAME="total_tracer"



!#######################################################################
! <FUNCTION NAME="total_volume">
!
! <DESCRIPTION>
! Compute total ocean tracer cell volume.
! </DESCRIPTION>
!
function total_volume (Thickness, index)

  type(ocean_thickness_type), intent(in) :: Thickness
  integer, intent(in)                    :: index

  real, dimension(isd:ied,jsd:jed) :: volume_t
  real                             :: total_volume
  integer                          :: k

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL,'==>Error from ocean_tracer_diag_mod (total_volume): module needs initialization ')
  endif 

  total_volume = 0.0
  do k=1,nk
    volume_t(:,:) = Thickness%dht(:,:,k,index)*Grd%dat(:,:)*Grd%tmask(:,:,k)
    if(have_obc) volume_t(:,:) = volume_t(:,:)*Grd%obc_tmask(:,:)
    total_volume  = total_volume + mpp_global_sum(Dom%domain2d,volume_t(:,:), global_sum_flag)
  enddo

end function total_volume
! </FUNCTION>  NAME="total_volume"



!#######################################################################
! <SUBROUTINE NAME="tracer_integrals">
!
! <DESCRIPTION>
! Compute some integrated tracer diagnostics. 
! </DESCRIPTION>
!
subroutine tracer_integrals (Time, Thickness, Tracer, Dens, Ext_mode)

  type(ocean_time_type), intent(in)          :: Time
  type(ocean_thickness_type), intent(in)     :: Thickness
  type(ocean_prog_tracer_type), intent(in)   :: Tracer
  type(ocean_density_type), intent(in)       :: Dens
  type(ocean_external_mode_type), intent(in) :: Ext_mode
  
  real, dimension(isd:ied,jsd:jed) :: tmp
  real, dimension(isd:ied,jsd:jed) :: volume_t
  real    :: tbar, tvar, tchg, ttot, t_cell_vol
  real    :: volume_total 
  integer :: i, j, k, n, tau, taum1, taup1

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL,'==>Error from ocean_tracer_diag_mod (tracer_integrals): module needs initialization ')
  endif 

  call mpp_clock_begin(id_tracer_integrals)

  taum1 = Time%taum1
  tau   = Time%tau
  taup1 = Time%taup1
  
  write(stdout(),*) ' '
  call write_timestamp(Time%model_time)
  write(stdout(),'(/a,a)') 'Tracer integrals for ' , trim(Tracer%name)
  
  volume_total = 0.0; tbar = 0.0; tvar = 0.0; tchg = 0.0; ttot = 0.0

  do k=1,nk

    volume_t(:,:) = Thickness%dht(:,:,k,tau)*Grd%dat(:,:)*Grd%tmask(:,:,k)
    if(have_obc) volume_t(:,:) = volume_t(:,:) * Grd%obc_tmask(:,:)
    volume_total  = volume_total + mpp_global_sum(Dom%domain2d,volume_t(:,:), global_sum_flag)

    tmp(:,:)      = Tracer%field(:,:,k,tau)*volume_t(:,:)
    tbar          = tbar + mpp_global_sum(Dom%domain2d,tmp(:,:), global_sum_flag)

    tmp(:,:)      = Tracer%field(:,:,k,tau)*Tracer%field(:,:,k,tau)*volume_t(:,:)
    tvar          = tvar + mpp_global_sum(Dom%domain2d,tmp(:,:), global_sum_flag)
    
    tmp(:,:)      = abs((Tracer%field(:,:,k,taup1)-Tracer%field(:,:,k,taum1))*dtimer)
    tchg          = tchg +  mpp_global_sum(Dom%domain2d,tmp(:,:), global_sum_flag)

  enddo 

  ttot = total_tracer(Tracer,Thickness,tau)
  tbar = tbar/volume_total
  tvar = tvar/volume_total - tbar**2

  write(stdout(), '(/a,a,a,es24.17)') 'Average  ',trim(Tracer%name),' = ', tbar
  write(stdout(), '(a,a,a,es24.17)')  'Variance ',trim(Tracer%name),' = ', tvar
  write(stdout(), '(a,a,a,es24.17)')  '|dT/dt|  ',trim(Tracer%name),' = ', tchg
  write(stdout(), '(a,a,a,es24.17/)') 'Total    ',trim(Tracer%name),' = ', ttot

  call mpp_clock_end(id_tracer_integrals)
  
end subroutine tracer_integrals
! </SUBROUTINE>  NAME="tracer_integrals"



!#######################################################################
! <SUBROUTINE NAME="tracer_land_cell_check">
!
! <DESCRIPTION>
! Check to be sure ocean tracer is zero over land 
! </DESCRIPTION>
!
subroutine tracer_land_cell_check (Time, T_prog)

  type(ocean_time_type), intent(in)        :: Time
  type(ocean_prog_tracer_type), intent(in) :: T_prog(:)
  integer :: i, j, k, n, num, taup1

  call mpp_clock_begin(id_tracer_land_cell_check)

  if(size(T_prog,1) /= num_prog_tracers) then 
    call mpp_error(FATAL,'==>Error from ocean_tracer_diag_mod (tracer_land_cell_check): passing T_prog with wrong dimensions')
  endif 

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL,'==>Error from ocean_tracer_diag_mod (tracer_land_cell_check): module needs initialization ')
  endif 

  taup1 = Time%taup1

  write (stdout(),'(1x,/a/)')'Locations (if any) where land cell tracer is non-zero...'
  num = 0
  do n=1,num_prog_tracers 
     do k=1,nk
        do j=jsc,jec
           do i=isc,iec
              if (Grd%tmask(i,j,k) == 0 .and. (T_prog(n)%field(i,j,k,taup1) /= 0.0) ) then
                  num = num + 1
                  if (num < 10) then
                    write(unit,9100) i+Dom%ioff,j+Dom%joff,k,Grd%xt(i,j),Grd%yt(i,j),Grd%zt(k),n,T_prog(n)%field(i,j,k,taup1)
                  endif 
              endif
           enddo
        enddo
     enddo
  enddo

  call mpp_max(num)
  if (num > 0) call mpp_error(FATAL)

9100  format(/' " =>Error: Land cell at (i,j,k) = ','(',i4,',',i4,',',i4,'),',&
          ' (lon,lat,dpt) = (',f7.2,',',f7.2,',',f7.0,'m) has tracer(n=',i2,') =',e10.3)

  call mpp_clock_end(id_tracer_land_cell_check)


end subroutine tracer_land_cell_check
! </SUBROUTINE>  NAME="tracer_land_cell_check"


!#######################################################################
! <SUBROUTINE NAME="volume_conservation">
!
! <DESCRIPTION>
! Compute change in volume over many time steps, and compare to the 
! input of volume through surface to check for volume conservation.
!
!============================================================
!
! threelevel scheme 
! 
! Here is the logic for the accumulation of the fluxes and 
! comparisons between volumes at the start and the end. 
!
! Consider accumulation over four leap-frog time steps. 
! Ignore time filtering.  
!
! eta(2) = eta(0) + 2deta*F(1)  taup1=2, taum1=0, tau=1 
!
! eta(3) = eta(1) + 2deta*F(2)  taup1=3, taum1=1, tau=2 
!
! eta(4) = eta(2) + 2deta*F(3)  taup1=4, taum1=2, tau=3 
!
! eta(5) = eta(3) + 2deta*F(4)  taup1=5, taum1=3, tau=4 
! 
! Hence,
!
! [eta(4) + eta(5)] = [eta(0) + eta(1)] + 2deta*[F(1)+F(2)+F(3)+F(4)]
! 
! For this example, we have 
!
! itts_volume=1 through itte_volume=4 for accumulating fluxes
!
! itt=itts_volume=1=tau we use taum1=0 and tau=1 to get starting volume
!
! itt=itte_volume=4=tau we use tau=4 and taup1=5 to get the final volume 
!
!============================================================
!
! twolevel scheme
! 
! Here is the logic for the accumulation of the fluxes and 
! comparisons between volumes at the start and the end. 
!
! Consider accumulation over four time steps. 
!
! eta(3/2) = eta(1/2) + deta*F(1)   taup1=3/2, taum1=1/2, tau=1 
!
! eta(5/2) = eta(3/2) + deta*F(2)   taup1=5/2, taum1=3/2, tau=2 
!
! eta(7/2) = eta(5/2) + deta*F(3)   taup1=7/2, taum1=5/2, tau=3 
!
! eta(9/2) = eta(7/2) + deta*F(4)   taup1=9/2, taum1=7/2, tau=4 
! 
! Hence,
!
! eta(9/2) = eta(1/2) + deta*[F(1)+F(2)+F(3)+F(4)]
! 
! For this example, we have 
!
! itts_volume=1 through itte_volume=4 for accumulating fluxes
!
! itt=itts_volume=1=tau we use taum1=1/2 to get starting volume
!
! itt=itte_volume=4=tau we use taup1=9/2 to get the final volume 
!
! </DESCRIPTION>
!
subroutine volume_conservation (Time, Thickness, Ext_mode, pme, river)

  type(ocean_time_type), intent(in)            :: Time
  type(ocean_thickness_type), intent(in)       :: Thickness
  type(ocean_external_mode_type), intent(in)   :: Ext_mode
  real, intent(in), dimension(isd:ied,jsd:jed) :: pme
  real, intent(in), dimension(isd:ied,jsd:jed) :: river

  real :: total_time, darea
  real :: volume_error, volume_error_rate
  real :: volume_input, volume_smooth, volume_chg

  integer :: i, j, k, n
  integer :: itt, taum1, tau, taup1

  logical, save :: first=.true.  

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL,'==>Error from ocean_tracer_diag (volume_conservation): module needs initialization ')
  endif 

  itt   = Time%itt
  taum1 = Time%taum1
  tau   = Time%tau
  taup1 = Time%taup1

  if (first) then
    first = .false.
    allocate (material_flux(isd:ied,jsd:jed))
    allocate (eta_smooth(isd:ied,jsd:jed))
    material_flux  = 0.0
    eta_smooth     = 0.0
    volume_start   = 0.0
    volume_final   = 0.0
    itts_volume    = itt+2
    itte_volume    = itts_volume-2 + nint(tracer_conserve_days*86400.0/dtts)
    if(itts_volume >= itte_volume) then 
      write(stdout(),*)'==>Warning: increase tracer_conserve_days to properly compute volume and tracer conservation.'
    endif 
  endif

  if(tendency==THREE_LEVEL) then 
      if (itt == itts_volume) then
          volume_start = 0.5*(total_volume(Thickness, taum1)+total_volume(Thickness, tau))    
      endif
      if (itt == itte_volume) then
          volume_final = 0.5*(total_volume(Thickness, tau)+total_volume(Thickness, taup1))
      endif
  elseif(tendency==TWO_LEVEL) then 
      if (itt == itts_volume) then
          volume_start = total_volume(Thickness, taum1)
      endif
      if (itt == itte_volume) then
          volume_final = total_volume(Thickness, taup1)
      endif
  endif

  if(itt >= itts_volume .and. itt <= itte_volume) then 
      do j=jsc,jec
         do i=isc,iec
            material_flux(i,j) = material_flux(i,j) + dteta*Grd%tmask(i,j,1)*Grd%dat(i,j) &
              *(pme(i,j) + river(i,j) + Ext_mode%eta_source(i,j) + Ext_mode%surface_smooth(i,j)  )
            eta_smooth(i,j) = eta_smooth(i,j) + dteta*Grd%tmask(i,j,1)*Grd%dat(i,j)*Ext_mode%surface_smooth(i,j)
         enddo
      enddo
  endif

  if (itt==itte_volume) then 

      if(have_obc) material_flux(:,:) = material_flux(:,:)*Grd%obc_tmask(:,:)
      if(have_obc) eta_smooth(:,:)    = eta_smooth(:,:)*Grd%obc_tmask(:,:)
      volume_input  = mpp_global_sum(Dom%domain2d,material_flux(:,:), global_sum_flag)
      volume_smooth = mpp_global_sum(Dom%domain2d,eta_smooth(:,:), global_sum_flag)

      total_time        = (itte_volume-itts_volume+1)*dtts+epsln
      volume_chg        = volume_final - volume_start
      volume_error      = volume_chg-volume_input
      volume_error_rate = volume_error/(total_time*Grd%tcellsurf) 

      write (stdout(),'(/50x,a/)') ' Measures of global integrated volume conservation over multiple time steps'
      write (stdout(),'(a,i10,a,es24.17,a)') ' Ocean volume at timestep ',itts_volume,' = ',volume_start,' m^3.'
      write (stdout(),'(a,i10,a,es24.17,a)') ' Ocean volume at timestep ',itte_volume,' = ',volume_final,' m^3.'

      write (stdout(),'(a,es24.17,a)')       ' Volume input via surface fluxes over time interval = ',volume_input,' m^3.'
      write (stdout(),'(a,es24.17,a)')       ' Volume input via eta_t smoother over time interval = ',volume_smooth,' m^3.'
      write (stdout(),'(a,es24.17,a)')       ' Change in ocean volume over time interval          = ',volume_chg,' m^3.'
      write (stdout(),'(a,es10.3,a)') ' Error in volume content change         = ',volume_error,' m^3.'
      write (stdout(),'(a,es10.3,a)') ' Error in rate of volume content change = ',volume_error_rate,' m/s.'
      write (stdout(),'(/)')

  endif

end subroutine volume_conservation
! </SUBROUTINE>  NAME="volume_conservation"



!#######################################################################
! <SUBROUTINE NAME="tracer_conservation">
!
! <DESCRIPTION>
! Compute change in global integrated tracer over many time steps,
! and compare to the input of tracer through the boundaries to 
! check for total tracer conservation.
!
! Accumulate fluxes as in the volume_conservation diagnostic. 
!
! </DESCRIPTION>
!
subroutine tracer_conservation (Time, Thickness, T_prog, T_diag, Ext_mode, Dens, pme, river)

  type(ocean_time_type), intent(in)            :: Time
  type(ocean_thickness_type), intent(in)       :: Thickness
  type(ocean_prog_tracer_type), intent(in)     :: T_prog(:)
  type(ocean_diag_tracer_type), intent(in)     :: T_diag(:)
  type(ocean_external_mode_type), intent(in)   :: Ext_mode
  type(ocean_density_type), intent(in)         :: Dens
  real, intent(in), dimension(isd:ied,jsd:jed) :: pme, river
 
  real, dimension(isd:ied,jsd:jed) :: tmp
  real :: tracer_chg, tracer_error, tracer_error_rate
  real :: tracer_source_input, tracer_smooth_input, tracer_flux_input
  real :: total_time, darea

  integer :: i, j, k, n
  integer :: itt
  integer :: tau, taum1, taup1

  logical, save :: first =.true.  

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL,'==>Error from ocean_tracer_diag (tracer_conservation): module needs initialization ')
  endif 

  if(size(T_prog,1) /= num_prog_tracers) then 
    call mpp_error(FATAL,'==>Error from ocean_tracer_diag (tracer_conservation): passing T_prog with wrong dimensions')
  endif 

  tau   = Time%tau
  taum1 = Time%taum1
  taup1 = Time%taup1
  itt   = Time%itt

  if (first) then
    first  = .false.
    allocate (tracer_flux(isd:ied,jsd:jed,num_prog_tracers))
    allocate (tracer_source(isd:ied,jsd:jed,num_prog_tracers))
    allocate (tracer_smooth(isd:ied,jsd:jed,num_prog_tracers))
    allocate (tracer_start(num_prog_tracers))
    allocate (tracer_final(num_prog_tracers))
    tracer_flux(:,:,:)   = 0.0
    tracer_source(:,:,:) = 0.0
    tracer_smooth(:,:,:) = 0.0
    tracer_start(:)      = 0.0
    tracer_final(:)      = 0.0
    itts_tracer          = itt+2
    itte_tracer          = itts_tracer-2 + nint(tracer_conserve_days*86400.0/dtts)
  endif

  if(tendency==THREE_LEVEL) then 
      if (itt == itts_tracer) then
          do n=1,num_prog_tracers
             tracer_start(n) = 0.5*(  total_tracer(T_prog(n),Thickness,taum1) &
                                    + total_tracer(T_prog(n),Thickness,tau))      
          enddo
      endif
      if (itt == itte_tracer) then
          do n=1,num_prog_tracers
             tracer_final(n) = 0.5*(  total_tracer(T_prog(n),Thickness,tau)  &
                                    + total_tracer(T_prog(n),Thickness,taup1))      
          enddo
      endif
  elseif(tendency==TWO_LEVEL) then 
      if (itt == itts_tracer) then
          do n=1,num_prog_tracers
             tracer_start(n) = total_tracer(T_prog(n),Thickness,taum1)    
          enddo
      endif
      if (itt == itte_tracer) then
          do n=1,num_prog_tracers
             tracer_final(n) = total_tracer(T_prog(n),Thickness,taup1)
          enddo
      endif
  endif

  if(itt >= itts_tracer .and. itt <= itte_tracer) then 

      do n=1,num_prog_tracers

         ! include river and pme contributions here 
         ! (must remove from tracer_source below) 
         do j=jsc,jec
            do i=isc,iec
               tracer_flux(i,j,n) = tracer_flux(i,j,n) + dtts*T_prog(n)%conversion*Grd%tmask(i,j,1)*Grd%dat(i,j) &
                    *(-T_prog(n)%btf(i,j)           + T_prog(n)%stf(i,j) + &
                       T_prog(n)%tpme(i,j)*pme(i,j) + T_prog(n)%triver(i,j)*river(i,j) )
            enddo
         enddo

         ! global integral of th_tendency leaves just the sources 
         ! (and other stuff to be removed below)
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  tracer_source(i,j,n) = tracer_source(i,j,n) + &
                                         dtts*T_prog(n)%conversion*Grd%tmask(i,j,k)*Grd%dat(i,j)*T_prog(n)%th_tendency(i,j,k)
               enddo
            enddo
         enddo

         do j=jsc,jec
            do i=isc,iec
               tracer_smooth(i,j,n) = tracer_smooth(i,j,n) + &
                    dtts*T_prog(n)%conversion*Grd%tmask(i,j,1)*Grd%dat(i,j)*T_prog(n)%surface_smooth(i,j)
            enddo
         enddo

         ! remove rivers, pme, and smooth from tracer source 
         ! (they are added to T_prog(n)%th_tendency inside ocean_rivermix_mod and ocean_tracers)
         k=1 
         do j=jsc,jec
            do i=isc,iec
               tracer_source(i,j,n) = tracer_source(i,j,n) &
                    -dtts*T_prog(n)%conversion*Grd%tmask(i,j,k)*Grd%dat(i,j)           &
                    *(T_prog(n)%triver(i,j)*river(i,j) + T_prog(n)%tpme(i,j)*pme(i,j))  &
                    -tracer_smooth(i,j,n)
            enddo
         enddo


         if(n==index_temp .and. index_frazil > 0) then 
             do j=jsc,jec
                do i=isc,iec
                   tracer_flux(i,j,n) = tracer_flux(i,j,n) + T_diag(index_frazil)%field(i,j,1)*Grd%dat(i,j)*Grd%tmask(i,j,1)
                enddo
             enddo
         endif

         if (have_obc) then
             call ocean_obc_tracer_flux(Time, T_prog(n), tmp, n, send_out = .false.)
             tracer_flux(:,:,n)   = (tracer_flux(:,:,n) + dtts*T_prog(n)%conversion*tmp(:,:))*Grd%obc_tmask(:,:)
             tracer_source(:,:,n) = tracer_source(:,:,n)*Grd%obc_tmask(:,:)
             tracer_smooth(:,:,n) = tracer_smooth(:,:,n)*Grd%obc_tmask(:,:)
         endif

      enddo

  endif

  if (itt==itte_tracer) then 

    write (stdout(),'(/50x,a/)') ' Measures of global integrated tracer conservation over multiple time steps'

    total_time = (itte_volume-itts_volume+1)*dtts+epsln

    do n=1,num_prog_tracers

      tracer_flux_input   = mpp_global_sum(Dom%domain2d,tracer_flux(:,:,n), global_sum_flag)
      tracer_source_input = mpp_global_sum(Dom%domain2d,tracer_source(:,:,n), global_sum_flag)
      tracer_smooth_input = mpp_global_sum(Dom%domain2d,tracer_smooth(:,:,n), global_sum_flag)
      tracer_chg          = tracer_final(n)-tracer_start(n)
      tracer_error        = tracer_chg-tracer_flux_input-tracer_source_input-tracer_smooth_input
      tracer_error_rate   = tracer_error/(total_time*Grd%tcellsurf) 

      if (n==index_temp) then 

        write (stdout(),'(a,i10,a,es24.17,a)') ' Ocean heat content at timestep ',itts_tracer,' = ',tracer_start(n),' J.'
        write (stdout(),'(a,i10,a,es24.17,a)') ' Ocean heat content at timestep ',itte_tracer,' = ',tracer_final(n),' J.'
        write (stdout(),'(a,es24.17,a)') ' Heat input by eta_t smoother over time interval   = ',tracer_smooth_input,' J.'
        write (stdout(),'(a,es24.17,a)') ' Heat input by sources over time interval          = ',tracer_source_input,' J.'
        write (stdout(),'(a,es24.17,a)') ' Heat input through boundaries over time interval  = ',tracer_flux_input,' J.'
        write (stdout(),'(a,es24.17,a)') ' Change in ocean heat content over time interval   = ',tracer_chg,' J.'
        write (stdout(),'(a,es10.3,a)') ' Error in         heat content change               =',tracer_error,' J.'
        write (stdout(),'(a,es10.3,a)') ' Error in rate of heat content change               =',tracer_error_rate,' W/m^2.'
        write (stdout(),'(/)')

      else

        write (stdout(),'(a,i10,a,es24.17,a)') ' Total '//trim(T_prog(n)%name)// &
                         ' mass at timestep ',itts_tracer,' = ',tracer_start(n),' kg'
        write (stdout(),'(a,i10,a,es24.17,a)') ' Total '//trim(T_prog(n)%name)// &
                         ' mass at timestep ',itte_tracer,' = ',tracer_final(n),' kg'
        write (stdout(),'(1x,a,es24.17,a)') trim(T_prog(n)%name)// &
                         ' input by eta_t smoother over time interval     = ',tracer_smooth_input,' kg.'
        write (stdout(),'(1x,a,es24.17,a)') trim(T_prog(n)%name)// &
                         ' input by sources over time interval            = ',tracer_source_input,' kg.'
        write (stdout(),'(1x,a,es24.17,a)') trim(T_prog(n)%name)// &
                         ' input through boundaries over time interval    = ',tracer_flux_input,' kg.'
        write (stdout(),'(a,es24.17,a)') ' Change in ocean '//trim(T_prog(n)%name)// &
                         ' mass over time interval  = ',tracer_chg,' kg.'
        write (stdout(),'(a,es10.3,a)') ' Error in         '//trim(T_prog(n)%name)//' mass change =', &
                                          tracer_error,' kg'
        write (stdout(),'(a,es10.3,a)') ' Error in rate of '//trim(T_prog(n)%name)//' mass change =', &
                                          tracer_error_rate,' kg/(m^2 sec)'
        write (stdout(),'(/)')
      endif 

    enddo 
  endif

end subroutine tracer_conservation
! </SUBROUTINE>  NAME="tracer_conservation"


  
!#######################################################################
! <SUBROUTINE NAME="diagnose_kappa_sort">
!
! <DESCRIPTION>
! Routine to diagnose the amount of mixing between classes of a 
! particular tracer.  Temperature is used as default.
! Method follows that used in the paper 
!
! Spurious diapycnal mixing associated with advection in a
! z-coordinate ocean model, 2000: S.M. Griffies, R.C.
! Pacanowski, and R.W. Hallberg. Monthly Weather Review, vol 128, 538--564.
!
! This diagnostic is most useful when computing the levels of 
! effective dia-tracer mixing occuring in a model.
!
! Algorithm notes:
!
! -assumes flat ocean bottom--non-flat bottoms loose the precise relation 
!  between sorted depth and true ocean depth.  This is a minor inconvenience.
!
! -defines some global arrays, so requires large memory.
!  this restriction can be removed if parallel sort is 
!  implemented.  so far, such has not been done.  
!  volunteers are welcome to contribute a parallel sort.  
!
! -Effective kappa is set to zero at bottom of bottom-most cell
! and top of top-most cell in order to ensure zero flux 
! conditions at the column boundaries.  This is not appropriate
! when running with surface and/or bottom no-flux conditions.  
!
! -Robust effective kappas require > 4 time steps--steps 1-4 corrupted.  
! -Uncomment the !!$ lines when link to the proper NAG routines. 
!
! </DESCRIPTION>
!
subroutine diagnose_kappa_sort(Time, Theta)

  type(ocean_time_type), intent(in)        :: Time
  type(ocean_prog_tracer_type), intent(in) :: Theta
  integer                                  :: tau, taum1, taup1

  real :: z_top(nk)      ! distance to top of resting T-cell from flat ocean bottom (m)
  real :: z_center(nk)   ! distance to center of resting T-cell from flat ocean bottom (m)
  real :: delta_z        ! vertical distance (m) between interpolated sorted parcels
  real :: kappa_plot(nk) ! diagnosed dianeutral diffusivity (m^2/sec) for plotting 
  real :: kappa_lay(nk)  ! diagnosed dianeutral diffusivity (m^2/sec) in layer 
  real :: thick_lay(nk)  ! thickness (m) of a density layer 
  real :: deriv_lay(nk)  ! vertical density gradient centered at top of sorted tracer cell
  real :: flux_lay(nk)   ! dianuetral tracer flux centered at top of sorted tracer cell
  real :: rho_lay(nk,3)  ! sorted density averaged onto layers 
  integer :: num_lay(nk) ! number of parcels within a density layer
  integer :: nsortpts    ! number points sorted = number grid points in global domain

  real, dimension(:,:), allocatable   ::  global_dat     ! for global area 
  real, dimension(:,:,:), allocatable ::  global_tmask   ! for global mask
  real, dimension(:,:,:), allocatable ::  global_tracer  ! for global tracer 

  real, dimension(:), allocatable ::  vol_sort  ! sorted volume of grid cell (m^3) associated with each fluid parcel 
  real, dimension(:), allocatable ::  rho_sort  ! sorted density (kg/m^3) of parcels
  real, dimension(:), allocatable ::  rho_int   ! interpolated density (kg/m^3)
  real, dimension(:), allocatable ::  zstar_int ! distance from bottom (m) of sorted density as interpolated
  real, dimension(:), allocatable ::  z_star    ! distance from bottom (m) of sorted density 
  real, dimension(:), allocatable ::  irank     ! array needed for NAG sorting routine 
  real, dimension(:), allocatable ::  deriv     ! vertical density gradient array 
  real, dimension(:), allocatable ::  deriv_dum ! dummy array

  integer :: i, j, k, m, itau, nsort, ifail, numzero, numlay
  real    :: thicklay, cellarea_r, kappa_prev, tmp 

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL,'==>Error from ocean_tracer_diag (diagnose_kappa_sort): module needs initialization ')
  endif 

  tau   = Time%tau
  taum1 = Time%taum1
  taup1 = Time%taup1

! NOTE: minimum vertical density gradient rho_grad_min is necessary to 
! avoid errors with truncation in the division by drho/dz when compute kappa.
! rho_grad_min corresponds roughly to the precision of the computation.
! Physically, with
!      
! N^2 = -(g/rho0)(drho/dz)
!      
! then rho_grad_min sets a minimum N^2 resolved. 
! This corresponds to a frequency f=N/2pi.  The typical 
! period of inertial oscillations in the deep ocean is 6hrs
! (Pickard and Emery, page 55-56).  In the upper ocean, it is
! 10-30 minutes, in pycnocline it is smaller still.  
! So to cover the majority of the ocean's stratification,
! we will want to set rho_grad_min to something smaller than 9e-6. 
!
! To bin the effective diffusivity, it is also useful to have a max
! vertical density gradient.        

! vertical position of center and top of cell relative to flat bottom   
  z_center(1) = Grd%dzw(nk)
  z_top(1)    = Grd%dzt(nk)
  do k=2,nk
     z_center(k) = z_center(k-1) + Grd%dzw(nk-k+1)
     z_top(k)    = z_top(k-1)    + Grd%dzt(nk-k+1)
  enddo
  if(debug_diagnose_mixingA) then
      write(stdout(),'(a)') ' ' 
      do k=nk,1,-1
         write(stdout(),'(a,i3,a,e16.9)')' z_center(',k,') = ',z_center(k) 
      enddo
      do k=nk,1,-1
         write(stdout(),'(a,i3,a,e16.9)')' z_top(',k,')    = ',z_top(k) 
      enddo
  endif

  allocate (global_tmask(ni,nj,nk)) ; global_tmask=0.0
  call mpp_global_field(Dom%domain2d, Grd%tmask, global_tmask)

! count the number of ocean points in domain 
  nsortpts=0 
  do k=1,nk      
     do j=1,nj
        do i=1,ni
           if(global_tmask(i,j,k) > 0.0) nsortpts=nsortpts+1
        enddo
     enddo
  enddo

! area factors 
  allocate (global_dat(ni,nj)) ; global_dat=0.0
  call mpp_global_field(Dom%domain2d, Grd%dat, global_dat)
  if(debug_diagnose_mixingA) then
      write(stdout(),'(a)') ' ' 
      write(stdout(),'(a,e16.9)')' in sort, total cross area of domain = ',Grd%tcellsurf 
      write(stdout(),'(a,e16.9)')' in sort, total volume of domain     = ',Grd%tcellsurf*Grd%zw(nk)
  endif
  cellarea_r = 1.0/(epsln + Grd%tcellsurf)


! do-loop over the three time levels.
! this needed for robust results.
! alternative of simply updating taup1 is not robust. 
  do itau=1,3


! initialize volume and density of sorted parcels to unsorted values.
! note: do not include rho0 in density calculation in order to preserve
! precision of temperature value.
     allocate (global_tracer(ni,nj,nk)) ; global_tracer=0.0
     call mpp_global_field(Dom%domain2d, Theta%field(:,:,:,itau), global_tracer)
     allocate (vol_sort(nsortpts))        
     allocate (rho_sort(nsortpts))        
     vol_sort=0.0 ; rho_sort=0.0
     nsort=0 
     do k=1,nk      
        do j=1,nj
           do i=1,ni
              if(global_tmask(i,j,k) == 1.0) then 
                  nsort=nsort+1
                  vol_sort(nsort) = global_dat(i,j)*Grd%dzt(k)
                  rho_sort(nsort) = -alpha*global_tracer(i,j,k)
              endif
           enddo
        enddo
     enddo
     deallocate(global_tracer)

     if(debug_diagnose_mixingB) then
         write(stdout(),'(a)') ' ' 
         do nsort=1,nsortpts 
            write(stdout(),'(a,i7,a,e16.9)')' before sorting rho_sort(',nsort,') = ',rho0+rho_sort(nsort) 
         enddo
         do nsort=1,nsortpts 
            write(stdout(),'(a,i7,a,e16.9)')' before sorting vol_sort(',nsort,') = ',vol_sort(nsort) 
         enddo
     endif

! Rank density in descending order
     allocate (irank(nsortpts))        
     irank=0 ; ifail=0
! Determine irank according to rho_sort
!!$     call M01DAF(rho_sort,1,nsortpts,'D',irank,ifail)
! Sort density in descending order
     ifail=0
!!$     call M01EAF(rho_sort,1,nsortpts,irank,ifail)
! Sort volume according to the density ranking
     ifail=0
!!$     call M01EAF(vol_sort,1,nsortpts,irank,ifail)
     deallocate (irank)        

     if(debug_diagnose_mixingB) then
         write(stdout(),'(a,i2)') 'diagnostic for itau= ',itau 
         do nsort=1,nsortpts 
            write(stdout(),'(a,i7,a,e16.9)')' after sorting rho_sort(',nsort,') = ',rho0+rho_sort(nsort) 
         enddo
         do nsort=1,nsortpts 
            write(stdout(),'(a,i7,a,e16.9)')' after sorting vol_sort(',nsort,') = ',vol_sort(nsort) 
         enddo
     endif

! Compute vertical height above the bottom for the sorted parcels.
! This is the value of the vertical coordinate for the center of
! mass of a parcel in the sorted state.
! zstar(1) = center of mass of the most dense parcel,
! whose density is rho_sort(1). 
! zstar(nsortpts) = center of mass for the least dense parcel,
! whose density is rho_sort(nsortpts).

     allocate (z_star(nsortpts))          
     z_star=0.0
     z_star(1) = 0.5*vol_sort(1)*cellarea_r
     do nsort=2,nsortpts
        z_star(nsort) = z_star(nsort-1) + 0.5*(vol_sort(nsort-1)+vol_sort(nsort))*cellarea_r 
     enddo
     if(debug_diagnose_mixingB) then
         write(stdout(),'(a,i2)') 'diagnostic for itau= ',itau 
         do nsort=1,nsortpts 
            write(stdout(),'(a,i7,a,e16.9)')' z_star(',nsort,') = ',z_star(nsort) 
         enddo
     endif

! vertical height to which the sorted density field is interpolated
     allocate (zstar_int(nsortpts))
     zstar_int=0.0
     delta_z=Grd%zw(nk)/nsortpts
     zstar_int(1)= delta_z
     do nsort=2,nsortpts  
        zstar_int(nsort)= zstar_int(nsort-1) + delta_z
     enddo
     if(debug_diagnose_mixingB) then
         write(stdout(),'(a,i2)') 'diagnostic for itau= ',itau 
         do nsort=1,nsortpts 
            write(stdout(),'(a,i7,a,e16.9)')' zstar_int(',nsort,') = ',zstar_int(nsort) 
         enddo
     endif

! Piecewise cubic Hermite interpolation of the sorted density to zstar_int
! This interpolation preserves monotonicity.

! The resulting "_int" variables represent equal-sized parcels
! of density rho_int stably stacked into a column separated by 
! distance delta_z.

! redefine endpoints for interpolation purposes 
     if(z_star(1) > zstar_int(1))               z_star(1) = 0.0
     if(z_star(nsortpts) < zstar_int(nsortpts)) z_star(nsortpts) = epsln+zstar_int(nsortpts)

     allocate (deriv(nsortpts))        
     allocate (deriv_dum(nsortpts))        
     allocate (rho_int(nsortpts))        
     deriv=0.0 ; deriv_dum=0.0 ; rho_int=0.0
     ifail=0
!!$     call E01BEF(nsortpts,z_star,rho_sort,deriv,ifail)
     ifail=0
!!$     call E01BGF(nsortpts,z_star,rho_sort,deriv,nsortpts,zstar_int,rho_int,deriv_dum,ifail) 
     deallocate (deriv)        
     deallocate (deriv_dum)        
     deallocate (rho_sort)        
     deallocate (z_star)        
     if(debug_diagnose_mixingB) then
         write(stdout(),'(a,i2)') 'diagnostic for itau= ',itau 
         do nsort=1,nsortpts 
            write(stdout(),'(a,i7,a,e16.9)')' rho_int(',nsort,') = ',rho_int(nsort) 
         enddo
     endif

! compute density averaged onto layers
     rho_lay(:,itau) = 0.0
     thick_lay(:)    = 0.0 
     num_lay(:)      = 0
     do nsort=1,nsortpts 
        do k=1,1
           if(zstar_int(nsort)  < z_top(k)) then
               num_lay(k)       = num_lay(k) + 1
               thick_lay(k)     = thick_lay(k) + delta_z
               rho_lay(k,itau)  = rho_lay(k,itau) + rho_int(nsort)
           endif
        enddo
        do k=2,nk-1
           if(zstar_int(nsort) >= z_top(k-1) .and. zstar_int(nsort) < z_top(k)) then
               num_lay(k)       = num_lay(k) + 1
               thick_lay(k)     = thick_lay(k) + delta_z
               rho_lay(k,itau)  = rho_lay(k,itau) + rho_int(nsort)
           endif
         enddo
        do k=nk,nk
           if(zstar_int(nsort) > z_top(k-1)) then
               num_lay(k)       = num_lay(k) + 1
               thick_lay(k)     = thick_lay(k) + delta_z
               rho_lay(k,itau)  = rho_lay(k,itau) + rho_int(nsort)
           endif
        enddo
     enddo
     deallocate(rho_int)
     deallocate(zstar_int)
     deallocate(vol_sort)

     numlay=0
     thicklay=0.0
     do k=1,nk
        numlay   = numlay   + num_lay(k)
        thicklay = thicklay + thick_lay(k)
     enddo

     do k=1,nk
        if(num_lay(k) > 0) then  
            rho_lay(k,itau) = rho_lay(k,itau)*delta_z/thick_lay(k)
        else
            rho_lay(k,itau) = 0.0
        endif
     enddo

  enddo ! end of itau do-loop 

  if(debug_diagnose_mixingA) then
      write(stdout(),'(/a,i7)') ' total number of parcels                         = ',nsortpts 
      write(stdout(),'(a,i7)') ' total number of parcels accounted for in layers = ',numlay
      write(stdout(),'(a,e16.9)') ' total layer thickness                           = ',thicklay
      write(stdout(),'(a)')' '
      do k=nk,1,-1
         write(stdout(),'(a,i3,a,i6,a,i3,a,e16.9)')' num_lay(',k,')= ',num_lay(k), ' thick_lay(',k,')= ',thick_lay(k)
      enddo
      write(stdout(),'(/a)') 'density at three time levels' 
      do k=nk,1,-1
         write(stdout(),'(a,i3,a,e16.9,a,i3,a,e16.9,a,i3,a,e16.9)')' rho_lay(taum1,',k,')= ',rho_lay(k,taum1)+rho0, &
              ' rho_lay(tau(',k,')= ',rho_lay(k,tau)+rho0,&
              ' rho_lay(taup1(',k,')= ',rho_lay(k,taup1)+rho0
      enddo
  endif

! compute derivatives across averaged layers.
! derivatives are defined at top of tracer cell,
! with deriv_lay(k=1) at top of bottom-most cell,
! and deriv_lay(k=nk)=0 at top of top-most cell. 
  deriv_lay(:) = 0.0
  do k=1,nk-1
     deriv_lay(k) = (rho_lay(k+1,taum1)-rho_lay(k,taum1))*Grd%dzwr(nk-k)
  enddo

  if(debug_diagnose_mixingA) then
      write(stdout(),'(/a)') 'vertical sorted density derivative across a layer' 
      do k=nk,1,-1
         write(stdout(),'(a,i3,a,e16.9)')' deriv_sort(',k,') = ',deriv_lay(k) 
      enddo
  endif

! compute effective diffusivity by diagnosing the  
! diapycnal flux of sorted density. The flux is 
! is defined on the top face of the t-cell.  
! start integration from the ocean bottom and work upwards. 
! Assume zero flux entering through the ocean bottom.
  flux_lay(:) = 0.0
  k=1
  flux_lay(k) =  -(rho_lay(k,taup1)-rho_lay(k,taum1))*Grd%dzt(nk-k+1)*dtimer
  do k=2,nk-1
     flux_lay(k) = flux_lay(k-1) -(rho_lay(k,taup1)-rho_lay(k,taum1))*Grd%dzt(nk-k+1)*dtimer           
  enddo
  if(debug_diagnose_mixingA) then
      write(stdout(),'(/a)')'vertical flux across a layer' 
      do k=nk,1,-1
         write(stdout(),'(a,i3,a,e16.9)')' flux_sort(',k,') = ',flux_lay(k) 
      enddo
  endif

  kappa_lay(:) = 0.0
  numzero = 0  
  do k=1,nk
     if(abs(deriv_lay(k)) >= rho_grad_min .and. abs(deriv_lay(k)) <= rho_grad_max) then
         kappa_lay(k) = -flux_lay(k)/deriv_lay(k)
     else 
         numzero = numzero + 1
         kappa_lay(k) = 0.0
     endif
  enddo

  if(smooth_kappa_sort>0) then 
    do m=1,smooth_kappa_sort
      kappa_prev = 0.25*kappa_lay(1)
      do k=2,nk-2
         tmp          =  kappa_lay(k)
         kappa_lay(k) =  kappa_prev + 0.5*kappa_lay(k) + 0.25*kappa_lay(k+1)
         kappa_prev   =  0.25*tmp
      enddo
    enddo
  endif

  if (id_kappa_sort > 0) then 
      do k=1,nk
         kappa_plot(k) = kappa_lay(nk-k+1)
      enddo
      used = send_data (id_kappa_sort, kappa_plot(:), Time%model_time)
  endif

if(debug_diagnose_mixingA) then
    write (stdout(),'(/a)') ' Effective diapycnal diffusivity (m^2/sec)'
    do k=nk,1,-1
       write(stdout(),'(a,i3,a,e16.9)')' kappa_sort(',k,') = ',kappa_lay(k) 
    enddo
endif

end subroutine diagnose_kappa_sort
! </SUBROUTINE>  NAME="diagnose_kappa_sort"

  
!#######################################################################
! <SUBROUTINE NAME="diagnose_kappa_simple">
!
! <DESCRIPTION>
! Routine to diagnose the amount of mixing between classes of a 
! particular tracer.  Temperature is used as default.
! Compute horizontal average of temp to define a stable profile.
! Evolution of this profile defines an effective diffusity. 
! This diffusivity is different than the one diagnosed
! from the adiabatic sorting approach.  The sorting approach is 
! more relevant.  The two approaches agree when there 
! is zero baroclinicity, and the present simple scheme is 
! useful ONLY for debugging the more complex sorting routine. 
! </DESCRIPTION>
!
subroutine diagnose_kappa_simple(Time, Theta)

  type(ocean_time_type), intent(in)        :: Time
  type(ocean_prog_tracer_type), intent(in) :: Theta
  integer                                  :: tau, taum1, taup1

  real :: kappa_simple(nk)  ! diagnosed dianeutral diffusivity (m^2/sec)
  real :: deriv_simple(nk)  ! vertical density gradient centered at top of sorted tracer cell
  real :: flux_simple(nk)   ! dianuetral tracer flux centered at top of sorted tracer cell
  real :: rho_simple(nk,3)  ! density 

  real, dimension(:,:), allocatable     ::  global_dat     ! for global area 
  real, dimension(:,:,:), allocatable   ::  global_tmask   ! for global mask
  real, dimension(:,:,:,:), allocatable ::  global_tracer  ! for global tracer 
  
  integer :: i, j, k
  real    :: cellarea_r 

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL,'==>Error from ocean_tracer_diag (diagnose_kappa_simple): module needs initialization ')
  endif 

  tau   = Time%tau
  taum1 = Time%taum1
  taup1 = Time%taup1

! NOTE: minimum vertical density gradient rho_grad_min is necessary to 
! avoid errors with truncation in the division by drho/dz when compute kappa.
! rho_grad_min corresponds roughly to the precision of the computation.
! Physically, with
!      
! N^2 = -(g/rho0)(drho/dz)
!      
! then rho_grad_min sets a minimum N^2 resolved. 
! This corresponds to a frequency f=N/2pi.  The typical 
! period of inertial oscillations in the deep ocean is 6hrs
! (Pickard and Emery, page 56).  In the upper ocean, it is
! 10-30 minutes.  So to cover the majority of the ocean's stratification,
! we will want to set rho_grad_min to something smaller than 9e-6. 
!
! To bin the effective diffusivity, it is also useful to have a max
! vertical density gradient.        

  allocate (global_dat(ni,nj)) ; global_dat=0.0
  call mpp_global_field(Dom%domain2d, Grd%dat, global_dat)
  allocate (global_tmask(ni,nj,nk)) ; global_tmask=0.0
  call mpp_global_field(Dom%domain2d, Grd%tmask, global_tmask)

! horizontal tracer area 
  if(debug_diagnose_mixingA) then
      write(stdout(),'(a)') ' ' 
      write(stdout(),'(a,e16.9)')' in simple, total cross area of domain = ',Grd%tcellsurf
      write(stdout(),'(a,e16.9)')' in simple, total volume of domain     = ',Grd%tcellsurf*Grd%zw(nk)
  endif
  cellarea_r = 1.0/(epsln + Grd%tcellsurf)

! horizontally average to get vertical density profile 
  allocate (global_tracer(ni,nj,nk,3)) ; global_tracer=0.0
  call mpp_global_field(Dom%domain2d, Theta%field(:,:,:,:), global_tracer)
  rho_simple(:,:) = 0.0 
  do k=1,nk      
     do j=1,nj
        do i=1,ni
           if(global_tmask(i,j,k) == 1.0) then 
             rho_simple(k,:) = -alpha*global_tracer(i,j,k,:)*global_dat(i,j) + rho_simple(k,:) 
           endif
        enddo
     enddo
     rho_simple(k,:) = rho_simple(k,:)*cellarea_r
  enddo
  deallocate(global_tracer)

! compute derivatives across averaged layers.
! derivatives are defined at top of tracer cell,
! with deriv_lay(k=1) at top of top-most cell,
! and deriv_lay(k=nk)=0 at top of bottom-most cell. 
  deriv_simple(:) = 0.0
  do k=2,nk
     deriv_simple(k) = (rho_simple(k-1,taum1)-rho_simple(k,taum1))*Grd%dzwr(k-1)
  enddo
  if(debug_diagnose_mixingA) then
    write(stdout(),'(/a)') 'density at three time levels' 
    do k=1,nk
      write(stdout(),'(a,i3,a,e16.9,a,i3,a,e16.9,a,i3,a,e16.9)')' rho_simple(taum1,',k,')= ',rho_simple(k,taum1)+rho0, &
                                                              ' rho_simple(tau(',k,')= ',rho_simple(k,tau)+rho0,&
                                                              ' rho_simple(taup1(',k,')= ',rho_simple(k,taup1)+rho0
    enddo 
    write(stdout(),'(/a)') 'vertical derivative of density' 
    do k=1,nk
      write(stdout(),'(a,i3,a,e16.9)')' deriv_simple(',k,') = ',deriv_simple(k) 
    enddo 
  endif 

! compute effective diffusivity by diagnosing the  
! diapycnal flux. The flux is defined on the top
! face of the t-cell. start integration from the 
! ocean bottom and work up. 
  flux_simple(:) = 0.0
  k=nk
  flux_simple(k) = -(rho_simple(k,taup1)-rho_simple(k,taum1))*Grd%dzt(k)*dtimer
  do k=nk-1,2,-1
     flux_simple(k) = flux_simple(k+1) - (rho_simple(k,taup1)-rho_simple(k,taum1))*Grd%dzt(k)*dtimer
  enddo
  if(debug_diagnose_mixingA) then
    write(stdout(),'(/a)')'vertical flux across a layer' 
    do k=1,nk
      write(stdout(),'(a,i3,a,e16.9)')' flux_simple(',k,') = ',flux_simple(k) 
    enddo 
  endif 

  do k=1,nk
     if(abs(deriv_simple(k)) >= rho_grad_min .and. abs(deriv_simple(k)) <= rho_grad_max) then
         kappa_simple(k) = -flux_simple(k)/deriv_simple(k)
     else 
         kappa_simple(k) = 0.0
     endif
  enddo

  if (id_kappa_simple > 0) then 
     used = send_data (id_kappa_simple, kappa_simple(:), Time%model_time)
  endif 

  if(debug_diagnose_mixingA) then
    write (stdout(),'(/a)') ' Effective diapycnal diffusivity (m^2/sec)'
    do k=1,nk
      write(stdout(),'(a,i3,a,e16.9)')' kappa_simple(',k,') = ',kappa_simple(k) 
    enddo 
  endif 

end subroutine diagnose_kappa_simple
! </SUBROUTINE>  NAME="diagnose_kappa_simple"



!#######################################################################
! <SUBROUTINE NAME="diagnose_depth_of_potrho">
! <DESCRIPTION>
! Diagnose depth (m) of a potential density surface.  Method uses linear 
! interpolation to find the depth of a potential rho surface.
! Scheme currently does not forward (backwards) interpolate if 
! rho surface lies within lowest (uppermost) grid cell. Nor does 
! it account for partial bottom cells or free surface.  Hence, depths
! for isopycnals near surface or bottom boundaries may be corrupted. 
! Diagnostic only makes sense when potrho is monotonically
! increasing with depth.
!
! Author: Harper.Simmons@noaa.gov
!         Zhi.Liang@noaa.gov
! </DESCRIPTION>
!
subroutine diagnose_depth_of_potrho(Time, Dens)

  type(ocean_time_type), intent(in)    :: Time
  type(ocean_density_type), intent(in) :: Dens

  real      :: w1,w2
  integer   :: potrho_nk
  integer   :: i, j, k, n
  real, dimension(isd:ied,jsd:jed,size(Dens%potrho_ref(:))) :: depth_of_potrho 

  potrho_nk = size(Dens%potrho_ref(:))

  if (.not.module_is_initialized) then 
      call mpp_error(FATAL,'==>Error from ocean_tracer_diag (diagnose_depth_of_potrho): module needs initialization ')
  endif

  depth_of_potrho(:,:,:)=-10.0
  do n=1,potrho_nk
     do j=jsc,jec
        do i=isc,iec
kloop:     do k=nk-1,1,-1
              if(    Dens%potrho(i,j,k) < Dens%potrho_ref(n)) then
                  if(Dens%potrho_ref(n) < Dens%potrho(i,j,k+1)) then
                      if(Grd%tmask(i,j,k+1) > 0) then 
                          W1= Dens%potrho_ref(n)   - Dens%potrho(i,j,k)
                          W2= Dens%potrho(i,j,k+1) - Dens%potrho_ref(n)
                          depth_of_potrho(i,j,n) = (Grd%zt(k+1)*W1 + Grd%zt(k)*W2)/(W1 + W2 + epsln)
                          exit kloop 
                      endif
                  endif
              endif
           enddo kloop
        enddo
     enddo
  enddo
  used = send_data (id_depth_of_potrho, depth_of_potrho(:,:,:), &
         Time%model_time, &
         is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=potrho_nk)

end subroutine diagnose_depth_of_potrho
! </SUBROUTINE>  NAME="diagnose_depth_of_potrho"



!#######################################################################
! <SUBROUTINE NAME="diagnose_depth_of_theta">
! <DESCRIPTION>
! Diagnose depth (m) of a potential temperature surface.  Method uses  
! linear interpolation to find the depth of a potential temp surface.
! Scheme currently does not forward (backwards) interpolate if 
! theta surface lies within lowest (uppermost) grid cell. Nor does 
! it account for partial bottom cells or free surface.  Hence, depths
! for theta surfaces near ocean surface or bottom boundaries may be 
! corrupted. Diagnostic only makes sense when theta is monotonically
! decreasing with depth.
!
! Author: Stephen.Griffies@noaa.gov
!         Zhi.Liang@noaa.gov 
! based on "diagnose_depth_of_potrho" by Harper.Simmons@noaa.gov
! </DESCRIPTION>
!
subroutine diagnose_depth_of_theta(Time, Dens, T_prog)

  type(ocean_time_type), intent(in)        :: Time
  type(ocean_density_type), intent(in)     :: Dens
  type(ocean_prog_tracer_type), intent(in) :: T_prog(:)

  real      :: w1,w2
  integer   :: theta_nk
  integer   :: i, j, k, n, tau
  real, dimension(isd:ied,jsd:jed,size(Dens%theta_ref(:))) :: depth_of_theta 

  theta_nk = size(Dens%theta_ref(:))
  tau      = Time%tau
 
  if (.not.module_is_initialized) then 
    call mpp_error(FATAL,'==>Error from ocean_tracer_diag (diagnose_depth_of_theta): module needs initialization ')
  endif 

  depth_of_theta(:,:,:)=-10.0
  do n=1,theta_nk
     do j=jsc,jec
        do i=isc,iec
kloop:     do k=nk-1,1,-1
              if(    Dens%theta_ref(n) < T_prog(index_temp)%field(i,j,k,tau)  ) then
                  if(T_prog(index_temp)%field(i,j,k+1,tau) < Dens%theta_ref(n)) then
                      if(Grd%tmask(i,j,k+1) > 0) then
                          W1= T_prog(index_temp)%field(i,j,k,tau) - Dens%theta_ref(n) 
                          W2= Dens%theta_ref(n) - T_prog(index_temp)%field(i,j,k+1,tau) 
                          depth_of_theta(i,j,n) = (Grd%zt(k+1)*W1 + Grd%zt(k)*W2)/(W1 + W2 + epsln)
                          exit kloop
                      endif
                  endif
              endif
           enddo kloop
        enddo
     enddo
  enddo
  used = send_data (id_depth_of_theta, depth_of_theta(:,:,:), &
         Time%model_time, &
         is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=theta_nk)

end subroutine diagnose_depth_of_theta
! </SUBROUTINE>  NAME="diagnose_depth_of_theta"


!#######################################################################
! <SUBROUTINE NAME="potrho_mixed_layer">
!
! <DESCRIPTION>
! Determine mixed layer depth and potential density at mixed layer base  
! according to depth at which buoyancy is greater than buoyancy_crit
! relative to the surface. Compute the buoyancy using potential 
! density, rather than the insitu density, since we aim for this 
! diagnostic to be comparable to diagnostics from isopcynal models. 
! </DESCRIPTION>
!
subroutine potrho_mixed_layer (Time, Thickness, Dens)

  type(ocean_time_type), intent(in)      :: Time 
  type(ocean_thickness_type), intent(in) :: Thickness
  type(ocean_density_type), intent(in)   :: Dens
  type(time_type)                        :: next_time 

  integer :: i, j, k
  real    :: potrho_mix_depth(isd:ied,jsd:jed)
  real    :: potrho_mix_base(isd:ied,jsd:jed) 
  real    :: depth, potrho_top, crit

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL,'==>Error from ocean_tracer_diag (potrho_mixed_layer): module needs initialization ')
  endif 

  crit = buoyancy_crit*rho0/grav 
  potrho_mix_depth = 0.0
  potrho_mix_base  = 0.0

  do j=jsc,jec
     do i=isc,iec
        depth=Thickness%dhwt(i,j,0)
        potrho_mix_base(i,j) = Dens%potrho(i,j,1) + crit
     enddo
  enddo
  if (id_potrho_mix_base > 0) then 
      used = send_data (id_potrho_mix_base, potrho_mix_base(:,:), &
             Time%model_time, rmask=Grd%tmask(:,:,1), &
             is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
  endif

  next_time = increment_time(Time%model_time, int(dtts), 0)
  if (need_data(id_potrho_mix_depth,next_time)) then
      potrho_mix_depth = 0.0
      do j=jsc,jec
         do i=isc,iec
            depth=max(0.0,Thickness%dhwt(i,j,0))
            do k=2,Grd%kmt(i,j)
               if( Dens%potrho(i,j,k) < potrho_mix_base(i,j)) then 
                   depth=depth+Thickness%dhwt(i,j,k)
               else
                   potrho_mix_depth(i,j) = depth + Thickness%dhwt(i,j,k) &
                                       *(potrho_mix_base(i,j)-Dens%potrho(i,j,k-1))&
                                       /(Dens%potrho(i,j,k)  -Dens%potrho(i,j,k-1)) 
                   exit  
               endif
            enddo
         enddo
      enddo
      if (id_potrho_mix_depth > 0) then 
          used = send_data (id_potrho_mix_depth, potrho_mix_depth(:,:), &
                 Time%model_time, rmask=Grd%tmask(:,:,1), &
                 is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
      endif
  endif

end subroutine potrho_mixed_layer
! </SUBROUTINE> NAME="potrho_mixed_layer"



!#######################################################################
! <SUBROUTINE NAME="send_total_tracer">
!
! <DESCRIPTION>
! Send total tracer to diagnostic manager.
! </DESCRIPTION>
!
subroutine send_total_tracer(Time, Thickness, Tracer, ntracer)

  type(ocean_time_type), intent(in)        :: Time
  type(ocean_thickness_type), intent(in)   :: Thickness
  type(ocean_prog_tracer_type), intent(in) :: Tracer
  integer, intent(in)                      :: ntracer 

  integer  :: tau
  real     :: ttracer

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL,'==>Error from ocean_tracer_diag (send_total_tracer): module needs initialization ')
  endif 

  tau   = Time%tau

  ttracer = total_tracer(Tracer,Thickness,tau)
  if(Tracer%name=='temp') then 
      ttracer = ttracer*1e-25
  elseif(Tracer%name(1:3)=='age') then 
      ttracer = ttracer
  else 
      ttracer = ttracer*1e-18
  endif
  used = send_data (id_ttracer(ntracer), ttracer, Time%model_time)

end subroutine send_total_tracer
! </SUBROUTINE>  NAME="send_total_tracer"


end module ocean_tracer_diag_mod


