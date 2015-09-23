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
module ocean_horz_diffuse_mod
! 
!<CONTACT EMAIL="Ronald.Pacanowski@noaa.gov"> R. C. Pacanowski 
!</CONTACT>
!
!<REVIEWER EMAIL="Stephen.Griffies@noaa.gov"> Stephen M. Griffies
!</REVIEWER>
!
!<OVERVIEW>
! Thickness weighted time tendency for tracer 
! from horizontal laplacian diffusion  
!</OVERVIEW>
!
!<DESCRIPTION>
! The diffusivity  used to determine the strength of the tendency can be 
! a general function of space yet it is constant in time.  A namelist 
! option exists that determines this diffusivity as a local function 
! of the grid spacing. 
!</DESCRIPTION>
!
! <INFO>
!
! <NOTE>
! The numerical implementation requires no calls to mpp_update_domains.  
! </NOTE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_horz_diffuse_lap_nml">
!  <DATA NAME="horz_diffuse_on" TYPE="logical">
!  Must be true to use this module
!  </DATA> 
!  <DATA NAME="alap" UNITS="m^2/sec" TYPE="real">
!  This is the value for the space-time constant Laplacian diffusivity. 
!  </DATA> 
!  <DATA NAME="tracer_mix_micom" TYPE="logical">
!  If .true., then the diffusivity is set according to a velocity scale times
!  the grid spacing. It is based on an approach recommended by 
!  Eric Chassignet that is used in the Miami Isopycnal Model (MICOM).  
!  </DATA> 
!  <DATA NAME="vel_micom" UNITS="m/sec" TYPE="real">
!  Velocity scale that is used for computing the MICOM diffusivity. 
!  </DATA> 
!  <DATA NAME="verbose_init" TYPE="logical">
!  For verbose writes during initialization 
!  </DATA> 
!</NAMELIST>
!
use constants_mod,       only: epsln
use diag_manager_mod,    only: register_diag_field, send_data
use fms_mod,             only: write_version_number, open_namelist_file, close_file, check_nml_error
use fms_mod,             only: FATAL, NOTE, stdout, stdlog
use mpp_domains_mod,     only: mpp_update_domains, CGRID_NE
use mpp_mod,             only: mpp_error

use ocean_domains_mod,   only: get_local_indices, set_ocean_domain
use ocean_obc_mod,       only: store_ocean_obc_tracer_flux  
use ocean_operators_mod, only: FMX, FMY, FDX_PT, FDY_PT, BDX_ET, BDY_NT
use ocean_types_mod,     only: ocean_grid_type, ocean_domain_type, ocean_time_type
use ocean_types_mod,     only: ocean_thickness_type, ocean_prog_tracer_type
use ocean_types_mod,     only: missing_value
use ocean_workspace_mod, only: wrk1 

implicit none

private

public ocean_horz_diffuse_init
public horz_diffuse

real    :: alap              = 0.0e3  ! horizontal tracer diffusivity (m^2/sec)
real    :: vel_micom         = 0.0    ! constant velocity scale (m/s) for setting micom diffusivity  
logical :: tracer_mix_micom  =.false. ! if true, diffusivity made a function of the grid spacing
logical :: verbose_init      =.true.  ! for verbose writes during initialization 

! for diagnostics 
logical :: used
integer, dimension(:), allocatable  :: id_xflux_diff
integer, dimension(:), allocatable  :: id_yflux_diff
integer, dimension(:), allocatable  :: id_h_diffuse

real, dimension(:,:,:), allocatable :: diff_cet ! horizontal diffusivity for eastern face of T cell (m^2/s)
real, dimension(:,:,:), allocatable :: diff_cnt ! horizontal diffusivity for northern face of T cell (m^2/s)

character(len=128) :: version=&
     '$Id$'
character (len=128) :: tagname = &
     '$Name$'

type(ocean_grid_type), pointer   :: Grd => NULL()
type(ocean_domain_type), pointer :: Dom => NULL()
type(ocean_domain_type), save    :: Dom_flux

real, dimension(:,:,:), allocatable :: fx, fy
integer :: isd, ied, jsd, jed, isc, iec, jsc, jec, nk
integer :: num_prog_tracers      = 0
logical :: horz_diffuse_on       = .false.
logical :: module_is_initialized = .FALSE.
logical :: have_obc=.false.

namelist /ocean_horz_diffuse_lap_nml/ horz_diffuse_on, alap, tracer_mix_micom, vel_micom, verbose_init

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_horz_diffuse_init">
!
! <DESCRIPTION>
! Initialize the horizontal laplacian diffusion module by 
! registering fields for diagnostic output and performing some 
! numerical checks to see that diffusivity is set appropriately.
! </DESCRIPTION>
!
subroutine ocean_horz_diffuse_init(Grid, Domain, Time, T_prog, dtime, obc, tmask_sigma)

  type(ocean_grid_type), intent(in), target    :: Grid
  type(ocean_domain_type),intent(in),target    :: Domain
  type(ocean_time_type), intent(in)            :: Time
  type(ocean_prog_tracer_type), intent(inout)  :: T_prog(:)
  real, intent(in)                             :: dtime
  logical, intent(in)                          :: obc
  real, intent(in), optional                   :: tmask_sigma(Domain%isc:Domain%iec,Domain%jsc:Domain%jec)

  integer :: ioun, io_status, i, j, k, n, num, ierr
  
  if ( module_is_initialized ) then 
     call mpp_error(FATAL, '==>Error from ocean_horz_diffuse_mod (ocean_horz_diffuse_init): module already initialized')
  endif 

  module_is_initialized = .TRUE.

  have_obc = obc

  num_prog_tracers = size(T_prog(:))

  call write_version_number( version, tagname )

  ! provide for namelist over-ride of defaults 
  ioun  = open_namelist_file()
  read (ioun,ocean_horz_diffuse_lap_nml,IOSTAT=io_status)
  write (stdout(),'(/)')
  write (stdout(),ocean_horz_diffuse_lap_nml)  
  write (stdlog(),ocean_horz_diffuse_lap_nml)
  ierr = check_nml_error(io_status, 'ocean_horz_diffuse_lap_nml')
  call close_file(ioun)

  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk = Grid%nk

  if(horz_diffuse_on) then 
    call mpp_error(NOTE, '==> NOTE: USING ocean_horz_diffuse_mod from lap.')
  else 
    call mpp_error(NOTE, '==> NOTE: NOT using ocean_horz_diffuse_mod.')
    return
  endif 

  if(alap > 0 .or. (tracer_mix_micom .and. vel_micom > 0)) then 
     write (stdout(),'(/1x,a)') ' ==> Note: USING horizontal Laplacian tracer diffusion'
     write (stdout(),'(/1x,a)') ' ==> Warning: horizontal Laplacian diffusion will generally contribute to unphysically'
     write (stdout(),'(7x,a)') 'large diapycnal mixing. Such may be damaging for longer simulations.'
     write (stdout(),'(7x,a)') 'It is recommended that the horizontal tracer diffusivity be set to zero, '
     write (stdout(),'(7x,a)') 'and one instead uses a nontrivial neutral mixing with ocean_neutral_physics.F90'
  endif

  write(stdout(),'(/a,f10.2)')'==> Note from ocean_horz_diffuse_mod: using forward time step of (secs)', dtime 

  Dom => Domain
  Grd => Grid

  call set_ocean_domain(Dom_flux, Grid, name='horz diff flux')

  allocate (diff_cet(isd:ied,jsd:jed,nk))
  allocate (diff_cnt(isd:ied,jsd:jed,nk))
  allocate (fx(isd:ied,jsd:jed,nk))
  allocate (fy(isd:ied,jsd:jed,nk))  

  diff_cet(:,:,:) = alap
  diff_cnt(:,:,:) = alap
  fx = 0.0
  fy = 0.0
  
  ! Micom background diffusivity
  ! Space scale set by grid size
  ! Velocity scale input via namelist 
  if(tracer_mix_micom) then
    do k=1,nk
      diff_cet(:,:,k) = vel_micom*(2.0*Grd%dxt(:,:)*Grd%dyte(:,:)/(Grd%dxt(:,:)+Grd%dyte(:,:)))
      diff_cnt(:,:,k) = vel_micom*(2.0*Grd%dxtn(:,:)*Grd%dyt(:,:)/(Grd%dxtn(:,:)+Grd%dyt(:,:)))
    enddo
    if(verbose_init) then 
      do j=jsc,jec
        write (stdout(),'(a,i4,a,e14.7,a)') ' MICOM Laplacian diffusivity at (isc,',j,',1) = ',diff_cet(isc,j,1),' m^2/s'
      enddo
    endif 
  endif 

  if(PRESENT(tmask_sigma)) then 
     do j=jsc,jec
        do i=isc,iec
           if(tmask_sigma(i,j) > 0.0) then   
              k = Grd%kmt(i,j)
              diff_cet(i,j,k) = diff_cet(i,j,k)*(1.0-tmask_sigma(i,j))
              diff_cnt(i,j,k) = diff_cnt(i,j,k)*(1.0-tmask_sigma(i,j))
           endif
        enddo
     enddo
  endif
  
  if(PRESENT(tmask_sigma) .or. tracer_mix_micom) then
    call mpp_update_domains (diff_cet(:,:,:), Dom%domain2d)
    call mpp_update_domains (diff_cnt(:,:,:), Dom%domain2d)
  endif

  ! register for diagnostics manager 
  allocate (id_xflux_diff(num_prog_tracers))
  allocate (id_yflux_diff(num_prog_tracers))
  allocate (id_h_diffuse(num_prog_tracers))
  id_xflux_diff = -1
  id_yflux_diff = -1
  id_h_diffuse  = -1

  do n=1,num_prog_tracers

     if (trim(T_prog(n)%name) == 'temp') then
         id_xflux_diff(n) = register_diag_field ('ocean_model', &
              trim(T_prog(n)%name)//'_xflux_diff', Grd%tracer_axes(1:3), &
              Time%model_time, 'thk wghtd i-diffusive heat flux', 'Watts/m',&
              missing_value=missing_value, range=(/-1.e16,1.e16/))
         id_yflux_diff(n) = register_diag_field ('ocean_model', &
              trim(T_prog(n)%name)//'_yflux_diff', Grd%tracer_axes(1:3), &
              Time%model_time, 'thk wghtd j-diffusive heat flux', 'Watts/m',&
              missing_value=missing_value, range=(/-1.e16,1.e16/))         
         id_h_diffuse(n) = register_diag_field ('ocean_model', &
              trim(T_prog(n)%name)//'_h_diffuse', Grd%tracer_axes(1:3), &
              Time%model_time, 'thk wghtd horz-diffusion of heat', 'Watts/m^2',&
              missing_value=missing_value, range=(/-1.e16,1.e16/))         
     else
         id_xflux_diff(n) = register_diag_field ('ocean_model',&
         trim(T_prog(n)%name)//'_xflux_diff', Grd%tracer_axes(1:3), &
              Time%model_time, 'thk wghtd i-diffusive flux of '//trim(T_prog(n)%name), 'm^2*kg/sec',&
              missing_value=missing_value, range=(/-1.e10,1.e10/))
         id_yflux_diff(n) = register_diag_field ('ocean_model',&
         trim(T_prog(n)%name)//'_yflux_diff', Grd%tracer_axes(1:3), &
              Time%model_time, 'thk wghtd j-diffusive flux of '//trim(T_prog(n)%name), 'm^2*kg/sec',&
              missing_value=missing_value, range=(/-1.e10,1.e10/))         
         id_h_diffuse(n) = register_diag_field ('ocean_model',&
         trim(T_prog(n)%name)//'_h_diffuse', Grd%tracer_axes(1:3), &
              Time%model_time, 'thk wghtd horz-diffusion of '//trim(T_prog(n)%name), 'm*kg/sec',&
              missing_value=missing_value, range=(/-1.e10,1.e10/))         
     endif

  enddo

  

end subroutine ocean_horz_diffuse_init
! </SUBROUTINE>  NAME="ocean_horz_diffuse_init"


!#######################################################################
! <SUBROUTINE NAME="horz_diffuse">
!
! <DESCRIPTION>
! This function computes the thickness weighted time tendency for tracer 
! from horizontal laplacian diffusion. 
! </DESCRIPTION>
!
subroutine horz_diffuse (Time, Thickness, Tracer, ntracer, diag_flag)

  type(ocean_time_type), intent(in)                  :: Time
  type(ocean_thickness_type), intent(in)             :: Thickness
  type(ocean_prog_tracer_type), intent(inout)        :: Tracer
  integer, intent(in)                                :: ntracer
  logical, intent(in), optional                      :: diag_flag 
  logical                                            :: send_diagnostics 

  real, dimension(isd:ied,jsd:jed)    :: tchg
  real, dimension(isd:ied,jsd:jed,2)  :: tracr
  integer                             :: i, j, k
  integer                             :: tau, taum1
 
  if(.not. horz_diffuse_on) then 
    return 
  endif 

  ! assume send_diagnostics=.true., unless diag_flag says it is false.
  if (PRESENT(diag_flag)) then 
    send_diagnostics = diag_flag
  else 
    send_diagnostics = .true.
  endif 

  if (.not. module_is_initialized ) then 
     call mpp_error(FATAL, '==>Error from ocean_horz_diffuse_mod (horz_diffuse): module must be initialized')
  endif 

  taum1 = Time%taum1
  tau   = Time%tau

  ! fx = flux component through "eastern"  face of T-cells at level k
  ! fy = flux component through "northern" face of T-cells at level k
  ! tracr(:,:,1) = Tracer%field at level k-1
  ! tracr(:,:,2) = Tracer%field at level k

  tracr(:,:,1) = Tracer%field(:,:,1,taum1)
  do k=1,nk
     tracr(:,:,2) = Tracer%field(:,:,k,taum1)
     fx(:,:,k)    = diff_cet(:,:,k)*FDX_PT(tracr(:,:,1:2),k)*FMX(Thickness%dht(:,:,k,tau)*Grd%tmask(:,:,k))
     fy(:,:,k)    = diff_cnt(:,:,k)*FDY_PT(tracr(:,:,1:2),k)*FMY(Thickness%dht(:,:,k,tau)*Grd%tmask(:,:,k))
     tracr(:,:,1) = tracr(:,:,2)
  enddo

  if(Grd%tripolar) call mpp_update_domains(fx, fy, Dom_flux%domain2d, gridtype=CGRID_NE)

  do k=1,nk
     wrk1(:,:,k) = (BDX_ET(fx(:,:,k)) + BDY_NT(fy(:,:,k)))
     do j=jsc,jec
        do i=isc,iec
           Tracer%th_tendency(i,j,k) = Tracer%th_tendency(i,j,k) + Grd%tmask(i,j,k)*wrk1(i,j,k)
        enddo
     enddo
  enddo

  ! send fluxes to diag_manager 
  ! minus sign is due to MOM-convention for physics fluxes. 
  if(send_diagnostics) then 

      if (id_xflux_diff(ntracer) > 0) then
          used = send_data(id_xflux_diff(ntracer), -1.0*fx(:,:,:)*Tracer%conversion, &
                 Time%model_time, rmask=Grd%tmask(:,:,:), &
                 is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)           
      endif
      if (id_yflux_diff(ntracer) > 0) then
          used = send_data(id_yflux_diff(ntracer), -1.0*fy(:,:,:)*Tracer%conversion, &
                 Time%model_time, rmask=Grd%tmask(:,:,:), &
                 is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)           
      endif
      if (id_h_diffuse(ntracer) > 0) then
          used = send_data(id_h_diffuse(ntracer), wrk1(isc:iec,jsc:jec,:)*Tracer%conversion, &
                 Time%model_time, rmask=Grd%tmask(isc:iec,jsc:jec,:))           
      endif

  endif

  if (have_obc) call store_ocean_obc_tracer_flux(Tracer,-1.0*fx, -1.0*fy)
  
end subroutine horz_diffuse
! </SUBROUTINE> NAME="horz_diffuse"


end module ocean_horz_diffuse_mod
      
      
