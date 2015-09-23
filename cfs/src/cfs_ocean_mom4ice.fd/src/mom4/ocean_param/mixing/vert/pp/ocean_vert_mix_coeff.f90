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
module ocean_vert_mix_coeff_mod
!
!<CONTACT EMAIL="Ronald.Pacanowski@noaa.gov"> R. C. Pacanowski 
!</CONTACT>
!
!<OVERVIEW>
! Vertical viscosity and diffusivity according Pacanowski and Philander (1981)
!</OVERVIEW>
!
!<DESCRIPTION>
! This module computes vertical viscosity and diffusivity according to 
! Pacanowski and Philander (1981).  This scheme is most effective for 
! studies of the tropical circulation.  It computes the vertical mixing
! coefficient based on the Richardson number.  
!</DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! R.C. Pacanowski and G. Philander 
! Parametrization of vertical mixing in numerical models of the tropical ocean
! Journal of Physical Oceanography (1981) vol 11, pages 1442--1451
! </REFERENCE>
!
! <NOTE>
! This parameterization was designed for equatorial models
! and may not do a good job in mid or high latitudes. Simulations
! in these regions (where vertical shear is small) are improved with
! the addition of solar short wave penetration into the ocean which 
! reduces buoyancy and enhances vertical mixing.
! </NOTE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_vert_mix_coeff_pp_nml">
!  <DATA NAME="wndmix" UNITS="m^2/sec" TYPE="real">
!  Minimum viscosity at bottom of 1st level to simulate 
!  missing high frequency windstress components.
!  </DATA> 
!  <DATA NAME="fricmx" UNITS="m^2/sec" TYPE="real">
!  Maximum mixing
!  </DATA> 
!  <DATA NAME="diff_cbt_back_pp" UNITS="m^2/sec" TYPE="real">
!  Space-time independent background vertical diffusivity 
!  thought to be that arising from internal waves. Note that 
!  if using Bryan-Lewis background diffusivity, then should 
!  set diff_cbt_back_pp=0.0. 
!  </DATA> 
!  <DATA NAME="visc_cbu_back_pp" UNITS="m^2/sec" TYPE="real">
!  Background vertical viscosity
!  </DATA> 
!</NAMELIST>

use constants_mod,       only: grav, pi, epsln, rho0r
use diag_manager_mod,    only: register_diag_field, send_data
use fms_mod,             only: open_namelist_file, check_nml_error, close_file, write_version_number
use fms_mod,             only: FATAL, stdout, stdlog
use mpp_mod,             only: mpp_error, mpp_chksum

use ocean_density_mod,   only: density_delta_z
use ocean_domains_mod,   only: get_local_indices
use ocean_types_mod,     only: ocean_grid_type, ocean_domain_type
use ocean_types_mod,     only: ocean_prog_tracer_type, ocean_diag_tracer_type
use ocean_types_mod,     only: ocean_velocity_type, ocean_density_type
use ocean_types_mod,     only: ocean_time_steps_type, ocean_time_type, ocean_thickness_type
use ocean_types_mod,     only: missing_value
use ocean_vert_mix_mod,  only: diff_cbt_back
use ocean_workspace_mod, only: wrk1 

implicit none

public ocean_vert_mix_coeff_init
public vertical_mix_coeff
private ri_for_pp

private

real, private, dimension(:,:,:), allocatable :: riu  ! Richardson number at base of U-cells
real, private, dimension(:,:,:), allocatable :: rit  ! Richardson number at base of T-cells

type(ocean_grid_type), pointer   :: Grd => NULL()
type(ocean_domain_type), pointer :: Dom => NULL()

integer :: isd, ied, jsd, jed, isc, iec, jsc, jec, nk

integer :: index_temp = -1 , index_salt = -1, num_prog_tracers = -1, num_diag_tracers = -1

character(len=128) :: version = &
     '$Id$'

character (len=128) :: tagname = &
     '$Name$'

! for diagnostics 
integer :: id_diff_cbt_pp_t= -1
integer :: id_diff_cbt_pp_s= -1

logical :: module_is_initialized = .FALSE.

real :: fricmx           = 50.0e-4 ! (m^2/s) max vertical mixing coefficient
real :: wndmix           = 10.0e-4 ! (m^2/s) min vertical mixing in level 1 to simulate wind mixing
real :: diff_cbt_back_pp = 1.0e-5  ! (m^2/s) background "diff_cbt"
real :: visc_cbu_back_pp = 1.0e-4  ! (m^2/s) background "visc_cbu"
real :: visc_cbu_limit   = 1.0e2   ! (m^2/s) largest allowable "visc_cbu" (reset below)
real :: diff_cbt_limit   = 1.0e2   ! (m^2/s) largest allowable "diff_cbt" (reset below)

namelist /ocean_vert_mix_coeff_pp_nml/  wndmix, fricmx, diff_cbt_back_pp, visc_cbu_back_pp

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_vert_mix_coeff_init">
!
! <DESCRIPTION>
! Initialization for the Pacanowski/Philander vertical mixing scheme
!
! input:
!
!   dzt    = thickness of vertical levels (m)
!
!   nk     = number of vertical levels
!
!   yt     = latitude of grid points (deg)
!
!   nj     = number of latitudes
!
!   error  = logical to signal problems
!
! output:
!
!   wndmix = min value for mixing at surface to simulate high freq
!
!            wind mixing (if absent in forcing). (m^2/sec)
!
!   fricmx = maximum mixing (m^2/sec)
!
!   diff_cbt_back_pp = background "diff_cbt" (m^2/sec)
!
!   visc_cbu_back_pp = background "visc_cbu" (m^2/sec)
!
!   diff_cbt_limit = largest "diff_cbt" (m^2/sec)
!
!   visc_cbu_limit = largest "visc_cbu" (m^2/sec)
!
!   error  = true if some inconsistency was found
!
! </DESCRIPTION>
!
subroutine ocean_vert_mix_coeff_init (Grid, Domain, Time, Time_steps, T_prog, T_diag)

  type(ocean_grid_type), intent(in), target   :: Grid
  type(ocean_domain_type), intent(in), target :: Domain
  type(ocean_time_type), intent(in)           :: Time
  type(ocean_time_steps_type), intent(in)     :: Time_steps
  type(ocean_prog_tracer_type), intent(in)    :: T_prog(:)
  type(ocean_diag_tracer_type), intent(in)    :: T_diag(:)

  real :: dzmin, extlat
  integer :: k, j, ioun, io_status, ierr, n
  
  if ( module_is_initialized ) then 
     call mpp_error(FATAL, '==>Error from ocean_vert_mix_coeff_mod (ocean_vert_mix_coeff_init): module is already initialized')
  endif 

  module_is_initialized = .TRUE.

  call write_version_number( version, tagname )

  ! provide for namelist over-ride of defaults 
  ioun = open_namelist_file()
  read  (ioun, ocean_vert_mix_coeff_pp_nml,iostat=io_status)
  write (stdout(),'(/)')
  write (stdout(), ocean_vert_mix_coeff_pp_nml)  
  write (stdlog(), ocean_vert_mix_coeff_pp_nml)
  ierr = check_nml_error(io_status, 'ocean_vert_mix_coeff_pp_nml')
  call close_file(ioun)

  write(stdout(),'(/a/)')'==>USING Pacanowski-Philander (pp) vertical mixing scheme.'
  write(stdout(),'(/a,f10.2)')'==>Note from ocean_vert_mix_coeff_mod: using forward time step for vert-frict of (secs)', &
                               Time_steps%dtime_u 
  write(stdout(),'(/a,f10.2)')'==>Note from ocean_vert_mix_coeff_mod: using forward time step for vert-diff  of (secs)', &
                               Time_steps%dtime_t 

  num_prog_tracers = size(T_prog)
  num_diag_tracers = size(T_diag)

  do n = 1, num_prog_tracers
     if (T_prog(n)%name == 'temp') index_temp = n
     if (T_prog(n)%name == 'salt') index_salt = n
  enddo

  if (index_temp < 1 .or. index_salt < 1) then 
    call mpp_error(FATAL,'==>Error in ocean_vert_mix_coeff_mod: temp and/or salt not present in tracer array')
  endif 

  Dom => Domain
  Grd => Grid

  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk = Grid%nk

  allocate (riu(isd:ied,jsd:jed,nk))
  riu(:,:,:)    = 0.0

  allocate (rit(isd:ied,jsd:jed,nk))
  rit(:,:,:)    = 0.0

  if (Time_steps%aidif == 1.0) then  ! Implicit mixing
    visc_cbu_limit = fricmx
    diff_cbt_limit = 1.0e2  ! m^2/s
  else                      ! Explicit mixing: in regions of gravitational instability set mixing limits to the maximum
                            ! allowed by CFL criterion. convective adjustment will also act on the instability.
    visc_cbu_limit = fricmx
    diff_cbt_limit = fricmx
  endif

  dzmin  = 1.e10  ! meters
  do k=1,nk
    dzmin = min(dzmin,Grid%dzt(k))
  enddo
  if (dzmin >= 25.0) then
    write (stdout(),'(/,(1x,a))') '==> Warning: "ppmix" may not work well with coarse vertical resolution'
  endif

  extlat = 0.0
  do j=jsc,jec
    extlat = max(abs(Grid%yt(isc,j)),extlat)
  enddo
  if (extlat > 10.0) then
    write (stdout(),'(/,(1x,a))')'==> Warning: "ppmix" may not work well outside the tropics   '&
  ,'              where vertical shear is small unless solar shortwave penetration into the ocean is '&
  ,'              accounted for by enabeling  "shortwave"         '
  endif

  if (Time_steps%aidif < 1.0) then ! Explicit mixing 
    do k=1,nk
      if ((Time_steps%dtime_t*fricmx)/Grid%dzt(k)**2 >= 0.5) then
        write (stdout(),'(/,(1x,a))')'==> Error: vertical diffusive criteria exceeded for "fricmx".  use a smaller'&
       ,'            "dtts" and/or  "fricmx" .... or use "aidif=1.0"'
        write (stdout(),'(a48,i3)') ' at level =',k
        call mpp_error(FATAL, '==>Error in ocean_vert_mix_coeff_mod: vertical diffusive criteria exceeded for "fricmx" ')
      endif
      if ((Time_steps%dtime_t*diff_cbt_limit)/Grid%dzt(k)**2 >= 0.5) then
        write (stdout(),'(/,(1x,a))')'==> Error: vertical diffusive criteria exceeded for "diff_cbt_limit". use a smaller'&
        ,'            "dtts" and/or  "diff_cbt_limit" ...or use "aidif=1.0"'
        write (stdout(),'(a48,i3)') ' at level =',k
        call mpp_error(FATAL, '==> Error in ocean_vert_mix_coeff_mod: vertical diffusive criteria exceeded for "diff_cbt_limit" ')
      endif
    enddo

    if ((Time_steps%dtime_u*fricmx)/dzmin**2 >= 0.5) then
      write (stdout(),'(/,(1x,a))')'==> Error: vertical diffusive criteria exceeded for "fricmx". use a smaller'&
      ,'            "dtuv" and/or "fricmx"  ...or use "aidif=1.0"'
      call mpp_error(FATAL,'==> Error in ocean_vert_mix_coeff_mod: vertical diffusive criteria exceeded for "fricmx" ')
    endif

    if ((Time_steps%dtime_u*visc_cbu_limit)/dzmin**2 >= 0.5) then
      write (stdout(),'(/,(1x,a))')'==> Error: vertical diffusive criteria exceeded for "visc_cbu_limit". use a smaller'&
      ,'            "dtuv" or "visc_cbu_limit"  ...or use "aidif=1.0"'
      call mpp_error(FATAL, '==> Error in ocean_vert_mix_coeff_mod: vertical diffusive criteria exceeded for "visc_cbu_limit" ')
    endif
  else  ! Implicit mixing
    write (stdout(),'(/,(1x,a))')'==> Warning: using "aidif=1.0" with "ppmix" uses '&
    ,'             variables defined at "tau" rather than at "tau-1"'
  endif

  id_diff_cbt_pp_t = register_diag_field('ocean_model','diff_cbt_pp_t',Grd%tracer_axes_wt(1:3), &
       Time%model_time, 'vert diffusivity from pp for temp', 'm^2/sec',&
       missing_value = missing_value, range=(/-1.e5,1.e5/))

  id_diff_cbt_pp_s = register_diag_field('ocean_model','diff_cbt_pp_s',Grd%tracer_axes_wt(1:3), &
       Time%model_time, 'vert diffusivity from pp for salt', 'm^2/sec',&
       missing_value = missing_value, range=(/-1.e5,1.e5/))


end subroutine ocean_vert_mix_coeff_init
! </SUBROUTINE>  NAME="ocean_vert_mix_coeff_init"


!#######################################################################
! <SUBROUTINE NAME="vertical_mix_coeff">
!
! <DESCRIPTION>
! This subroutine computes the vertical diffusivity and viscosity
! according to the Pacanowski and Philander scheme. Mixing coefficients  
! are space and time dependent. 
!
! inputs:
!
!  nk              = number of vertical levels                               <BR/>
!  grav            = gravity (m/sec^2)                                       <BR/>  
!  fricmx          = max viscosity (m^2/sec)                                 <BR/> 
!  wndmix          = min viscosity at bottom of 1st level to simulate        <BR/>
!                    missing high frequency windstress components (m^2/sec)  <BR/>
!  visc_cbu_back_pp = background "visc_cbu" (m^2/sec)                        <BR/>  
!  diff_cbt_back_pp = background "diff_cbt" (m^2/sec)                        <BR/>  
!  visc_cbu_limit  = largest "visc_cbu" in regions of gravitational          <BR/> 
!                    instability (m^2/sec)                                   <BR/> 
!  diff_cbt_limit  = largest "diff_cbt" in regions of gravitational          <BR/>  
!                    instability (m^2/sec)                                   <BR/>
!  riu             = richardson number at bottom of U cells                  <BR/>   
!  rit             = richardson number at bottom of T cells                  <BR/>
!
! outputs:
!
!  visc_cbu = viscosity at bottom of U cells (m^2/s)  <BR/>
!  diff_cbt = diffusion at bottom of T cells (m^2/s) 
!
! </DESCRIPTION>
!
  subroutine vertical_mix_coeff(aidif, Time, Thickness, Velocity, T_prog, T_diag, Dens, swflx, sw_frac, pme, &
                                river, visc_cbu, diff_cbt, hblt_depth)

  real, intent(in)                                     :: aidif
  type(ocean_time_type), intent(in)                    :: Time
  type(ocean_thickness_type), intent(in)               :: Thickness
  type(ocean_velocity_type), intent(in)                :: Velocity
  type(ocean_prog_tracer_type), intent(in)             :: T_prog(:)
  type(ocean_diag_tracer_type), intent(in)             :: T_diag(:)
  type(ocean_density_type), intent(in)                 :: Dens
  real, intent(in), dimension(isd:ied,jsd:jed)         :: swflx
  real, intent(in), dimension(isd:ied,jsd:jed,0:nk)    :: sw_frac
  real, intent(in), dimension(isd:ied,jsd:jed)         :: pme
  real, intent(in), dimension(isd:ied,jsd:jed)         :: river
  real, intent(in), dimension(isd:ied,jsd:jed)         :: hblt_depth
  real, intent(inout), dimension(isd:ied,jsd:jed,nk)   :: visc_cbu
  real, intent(inout), dimension(isd:ied,jsd:jed,nk,2) :: diff_cbt

  integer :: i, j, k, n
  real    :: t1, t2
  logical :: used

  if ( .not. module_is_initialized ) then 
     call mpp_error(FATAL, '==>Error from ocean_vert_mix_coeff (vertical_mix_coeff): module must be initialized')
  endif 

  ! compute gradient richardson number at base of T cells and U cells

  call ri_for_pp (Time, Thickness, aidif, Velocity, T_prog(index_temp), T_prog(index_salt), &
                  Dens%pressure_at_depth, Dens%rho_taum1)

  ! viscosity and diffusivity are on bottom of T and U cells.

  do k=1,nk-1
    do j=jsc,jec
      do i=isc,iec
        t1                = 1.0/(1.0 + 5.0*riu(i,j,k))
        visc_cbu(i,j,k)   = fricmx*t1**2 + visc_cbu_back_pp
        t2                = 1.0/(1.0 + 5.0*rit(i,j,k))
        diff_cbt(i,j,k,:) = fricmx*t2**3 + diff_cbt_back_pp 
      enddo
    enddo
  enddo

  ! limit mixing coeffs on bottom of cells in unstable regions

  do k=1,nk-1
    do j=jsc,jec
      do i=isc,iec
        if (riu(i,j,k) < 0.0) visc_cbu(i,j,k)   = visc_cbu_limit
        if (rit(i,j,k) < 0.0) diff_cbt(i,j,k,:) = diff_cbt_limit
      enddo
    enddo
  enddo

  ! approximation for high freq wind mixing near the surface
  ! set no flux through bottom of bottom level "nk"

  k = 1
  where (diff_cbt(isc:iec,jsc:jec,k,:) < wndmix) diff_cbt(isc:iec,jsc:jec,k,:) = wndmix
  diff_cbt(isc:iec,jsc:jec,nk,:) = 0.0
  where (visc_cbu(isc:iec,jsc:jec,k) < wndmix) visc_cbu(isc:iec,jsc:jec,k) = wndmix
  visc_cbu(isc:iec,jsc:jec,nk) = 0.0

!-----------------------------------------------------------------------
!     incorporate background diffusivities that could be from 
!     bryan-lewis or from a table 
!-----------------------------------------------------------------------
  do k=1,nk-1
     do i=isc,iec
        do j=jsc,jec
           diff_cbt(i,j,k,1) = max(diff_cbt(i,j,k,1),diff_cbt_back(i,j,k))
        enddo
     enddo
  enddo
  diff_cbt(isc:iec,jsc:jec,:,2) = diff_cbt(isc:iec,jsc:jec,:,1)

  if (id_diff_cbt_pp_t > 0) used = send_data(id_diff_cbt_pp_t, diff_cbt(:,:,:,1), &
                                   Time%model_time, rmask=Grd%tmask(:,:,:), &
                                   is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  if (id_diff_cbt_pp_s > 0) used = send_data(id_diff_cbt_pp_s, diff_cbt(:,:,:,2), &
                                   Time%model_time, rmask=Grd%tmask(:,:,:), &
                                   is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)

end subroutine vertical_mix_coeff
! </SUBROUTINE> NAME="vertical_mix_coeff"


subroutine ri_for_pp (Time, Thickness, aidif, Velocity, Temp, Salt, pressure_at_depth, rho_taum1)

! compute richardson number for the pp scheme

  type(ocean_time_type), intent(in)               :: Time
  type(ocean_thickness_type), intent(in)          :: Thickness
  real, intent(in)                                :: aidif
  type(ocean_velocity_type), intent(in)           :: Velocity
  type(ocean_prog_tracer_type), intent(in)        :: Temp
  type(ocean_prog_tracer_type), intent(in)        :: Salt
  real, intent(in), dimension(isd:ied,jsd:jed,nk) :: pressure_at_depth
  real, intent(in), dimension(isd:ied,jsd:jed,nk) :: rho_taum1

  integer :: i, j, k, taum1
  real    :: active_cells, fx, t1

  fx   = -0.25*grav*rho0r

  taum1 = Time%taum1

  ! compute density difference across bottom of T cells at tau-1

  wrk1(:,:,:) = density_delta_z(rho_taum1(:,:,:), Salt%field(:,:,:,taum1), &
                                  Temp%field(:,:,:,taum1), pressure_at_depth(:,:,:))

  ! compute richardson numbers on bottom of U cells

  do k=1,nk-1
    do j=jsd,jed-1
      do i=isd,ied-1
        t1 = fx*Thickness%dhwt(i,j,k)
        riu(i,j,k) = t1*Grd%umask(i,j,k+1)*(wrk1(i,j+1,k) + wrk1(i+1,j+1,k) + wrk1(i,j,k)   + wrk1(i+1,j,k)) /&
                      ((Velocity%u(i,j,k,1,taum1) - Velocity%u(i,j,k+1,1,taum1))**2 + &
                       (Velocity%u(i,j,k,2,taum1) - Velocity%u(i,j,k+1,2,taum1))**2 + epsln)
      enddo
    enddo
  enddo

  ! compute richardson numbers on bottom of T cells as average
  ! of four nearest richardson numbers on bottom of U cells.
  ! (do not consider land cells in the average... only active ones)

  do k=1,nk-1
    do j=jsc,jec
      do i=isc,iec
        active_cells = Grd%umask(i,j,k+1) + Grd%umask(i-1,j,k+1) + Grd%umask(i,j-1,k+1) + Grd%umask(i-1,j-1,k+1) + epsln
        rit(i,j,k)   = (riu(i,j,k) + riu(i-1,j,k) + riu(i,j-1,k) + riu(i-1,j-1,k))/active_cells

        ! make sure no static instability exists (one that is not seen
        ! by the Richardson number).  This may happen due to
        ! horizontal averaging used in calculating the Richardson
        ! number.

        if (rit(i,j,k) > 0.0 .and. wrk1(i,j,k) > 0.0) then
          rit(i,j,k) = -10.0
        endif
      enddo
    enddo
  enddo

end subroutine ri_for_pp

end module ocean_vert_mix_coeff_mod
