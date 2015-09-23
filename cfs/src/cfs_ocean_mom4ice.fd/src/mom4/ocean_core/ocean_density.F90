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
! Modified: Dave Behringer
!           David.Behringer@noaa.gov
module ocean_density_mod
!
!<CONTACT EMAIL="Stephen.Griffies@noaa.gov"> S.M. Griffies 
!</CONTACT>
!
!<OVERVIEW>
! Compute density and related quantities.
!</OVERVIEW>
!
!<DESCRIPTION>
! This module computes the in-situ density and its partial derivatives with 
! respect to potential temperature and with respect to salinity.  
!
! Based on McDougall, Wright, Jackett, and Feistel (2002).  This 
! equation of state is valid over the range: 
!
! 0psu <= salinity <= 40 psu
!
! -3C <= theta <= 40C
!
! 0dbar <= pressure <= 8000dbar 
!
!  Input variables are the following:
!
!  salinity in psu
!
!  potential temperature (theta) in deg C
!
!  pressure in dbars  (1bar = 10dbar = 10^5 Newton/m^2 = 10^5 Pascals). 
!  Note that in the ocean, pressure increases roughly by 1dbar for each meter depth.
!  Also note that pressure is the "gauge" pressure, which is the absolute pressure
!  minus the pressure of a standard atmosphere, which is 10.1325 dbars.
!
! check values (kindly provided by David Jackett)                     <BR/>  
!  rho(s=20psu,theta=20C,p=1000dbar)   = 1017.72674313979 (kg/m^3)    <BR/>
!  alpha(s=20psu,theta=20C,p=1000dbar) = 2.524181985549684e-4 (1/C)   <BR/>
!  beta(s=20psu,theta=20C,p=1000dbar)  = 7.382804621244401e-4 (1/psu) <BR/> 
!
! This equation of state should be suitable for all purposes of realistic 
! ocean climate modeling. 
!
! B. Linear equation for use in idealized Boussinesq studies
! 
! This equation renders density a linear function of potential 
! temperature.  All nonlinearities are ignored, as are salinity and 
! pressure effects.  Since there are no compressibility effects in 
! this equations of state, it is only appropriate for Boussinesq
! studies.
!
! The valid range for T and S is arbitrary for linearized density.
! However the range is restricted to the range for the standard EOS
! to keep density gradients within reasonable limits.
! So valid ranges are restricted to s=0 to 50 psu, t=-10 to 50 deg C
!
! </DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! Feistel (2003)
! A new extended Gibbs thermodynamic potential of seawater
! Progress in Oceanography. vol 58, pages 43-114.
! </REFERENCE>
!
! <REFERENCE>
! McDougall, Jackett, Wright, and Feistel (2002)
! Accurate and computationally efficient algorithms for 
! potential temperatue and density of seawater
! Journal of Atmospheric and Oceanic Technology, submitted 2002
! </REFERENCE>
!
! <REFERENCE>
! Jackett, McDougall, Feistel, Wright, and Griffies (2004)
! Updated algorithms for density, potential temperature, 
! conservative temperature, and freezing temperature of 
! seawater.  
! Journal of Atmospheric and Oceanic Technology, 2004 submitted
! </REFERENCE>
!
! <REFERENCE>
! S.M. Griffies, M.J. Harrison,  R.C. Pacanowski, and A. Rosati
! A Technical Guide to MOM4 (2003)
! </REFERENCE>
!
! <REFERENCE>
! S.M. Griffies, R.C. Pacanowski, R.M. Schmidt, and V. Balaji
! Tracer Conservation with an Explicit Free Surface Method for 
! Z-coordinate Ocean Models
! Monthly Weather Review (2001) vol 129 pages 1081--1098
! </REFERENCE>
!
! <NOTE>
! Density is computed as a function of potential temperature (C), salinity (psu),
! and in-situ pressure (dbar).  The pressure contribution includes that from 
! the free surface height and the atmospheric pressure.  Because the baroclinic
! component of the hydrostatic pressure is not known until the density is known, 
! the baroclinic pressure contribution to density is lagged by a single time step.  
! rho(tau) = rho[theta(tau),s(tau), p_atm(tau) + p_fs(tau) + p_baroclinic(tau-1)]  
! </NOTE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_density_nml">
!
!  <DATA NAME="s_test" UNITS="psu" TYPE="real">
!  Salinity for testing the EOS.
!  </DATA> 
!
!  <DATA NAME="press_standard" UNITS="dbar" TYPE="real">
!  Standard atmospheric pressure (dbar).  The realistic 
!  EOS used in mom4 requires gauge pressuer as an argument
!  rather than absolute pressure.  Gauge pressure is 
!  absolute pressure minus a standard atmospheric pressure 
!  of 10.1325dbar.  
!  For models that do have a realistic atmospheric loading, then it
!  is appropriate to remove 10.1325dbar prior to computing the EOS.
!  For those cases with zero atmospheric pressure, then it is not
!  necessary to remove the standard atmosphere.  As most model are
!  presently run with zero atmospheric pressure, the default for the 
!  press_standard is 0.0.   
!  </DATA> 
!
!  <DATA NAME="t_test" UNITS="C" TYPE="real">
!  Potential temperature for testing the EOS.
!  </DATA> 
!  <DATA NAME="p_test" UNITS="dbar" TYPE="real">
!  Gauge pressure for testing the EOS.
!  </DATA> 
!
!  <DATA NAME="linear_eos" TYPE="logical">
!  Set to true if wish to use the linear equation of state.  
!  </DATA>
!  <DATA NAME="alpha_linear_eos" TYPE="real">
!  Constant "thermal expansion coefficient" for EOS 
!  rho = rho0 - alpha_linear_eos*theta + beta_linear_eos*salinity
!  </DATA>
!  <DATA NAME="beta_linear_eos" TYPE="real">
!  Constant "saline contraction coefficient" for EOS 
!  rho = rho0 - alpha_linear_eos*theta + beta_linear_eos*salinity
!  </DATA>
!
!  <DATA NAME="potrho_press" UNITS="dbar" TYPE="real">
!  Gauge pressure for computing diagnostic potential density 
!  </DATA> 
!  <DATA NAME="potrho_min" UNITS="kg/m^3" TYPE="real">
!  Minimum potential density used to partition vertical according to potential density.  
!  </DATA> 
!  <DATA NAME="potrho_max" UNITS="kg/m^3" TYPE="real">
!  Maximum potential density used to partition vertical according to potential density.  
!  </DATA> 
!
!  <DATA NAME="theta_min" UNITS="C" TYPE="real">
!  Minimum potential temperature used to partition vertical according to theta.
!  </DATA> 
!  <DATA NAME="theta_max" UNITS="C" TYPE="real">
!  Maximum potential temperature used to partition vertical according to theta.
!  </DATA> 
!  <DATA NAME="layer_nk" TYPE="integer">
!  Number of classes used to partition vertical according to potential density
!  or potential temperature. Used for diagnostics. 
!  </DATA> 
!
!  <DATA NAME="debug_density" TYPE="logical">
!  For debugging nonlinear equation of state 
!  </DATA>
!
!</NAMELIST>

use constants_mod,       only: epsln, rho0, grav, c2dbars, rho0r
use diag_manager_mod,    only: register_diag_field, diag_axis_init
use diag_manager_mod,    only: need_data, send_data
use fms_mod,             only: write_data, read_data, write_version_number, mpp_error, FATAL
use fms_mod,             only: open_namelist_file, check_nml_error, close_file, file_exist
use mpp_domains_mod,     only: mpp_update_domains
use mpp_mod,             only: stdout, stdlog, mpp_chksum
use platform_mod,        only: i8_kind
use time_manager_mod,    only: time_type, increment_time, get_date

use ocean_domains_mod,   only: get_local_indices
use ocean_operators_mod, only: REMAP_T_TO_V_NOCONSERVE
use ocean_pressure_mod,  only: pressure_in_dbars
use ocean_types_mod,     only: ocean_domain_type, ocean_grid_type, ocean_thickness_type
use ocean_types_mod,     only: ocean_time_type, ocean_time_steps_type
use ocean_types_mod,     only: ocean_prog_tracer_type, ocean_density_type, ocean_external_mode_type
use ocean_util_mod,      only: write_timestamp
use ocean_workspace_mod, only: wrk1, wrk2

implicit none

private

real, dimension(:), allocatable     :: a, b                  ! polynomial coefficients in the realistic EOS
real, dimension(:,:,:), allocatable :: denr                  ! reciprocal of denominator for relastic EOS
logical                             :: linear_eos=.false.        

real :: a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11      ! polynomial coefficients in the realistic EOS
real :: b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12 ! polynomial coefficients in the realistic EOS
real :: two_a2, three_a3, two_a6, two_a8, two_a11             ! multiplied constants for density partial derivs
real :: two_b2, three_b3, four_b4, three_b7                   ! multiplied constants for density partial derivs
real :: onep5_b8, onep5_b9, two_b9, three_b11                 ! multiplied constants for density partial derivs

! some test settings 
real :: s_test=20.0
real :: t_test=20.0
real :: p_test=1000.0
real :: mjwf_rho  =0.101772674313979e4 
real :: mjwf_alpha=0.2524181985549684e-3 
real :: mjwf_beta =0.7382804621244401e-3

! for linear EOS 
real :: alpha_linear_eos=0.255
real :: beta_linear_eos =0.0

! time steps 
real :: dtts =0.0

! for diagnostics manager 
integer :: id_drhodtheta=-1
integer :: id_drhodsalt=-1
integer :: id_press=-1
integer :: id_rho=-1
integer :: id_pot_rho=-1
integer :: id_pot_rho_0=-1
integer :: id_pot_rho_et=-1
integer :: id_pot_rho_nt=-1
integer :: id_pot_rho_wt=-1
logical :: used

#include <ocean_memory.h>

type(ocean_domain_type), pointer :: Dom =>NULL()
type(ocean_grid_type), pointer   :: Grd =>NULL()

character(len=128) :: version=&
       '$Id$'
character (len=128) :: tagname = &
       '$Name$'

public ocean_density_init
public ocean_density_rstrt
public ocean_density_end
public update_ocean_density
public update_ocean_density_taup1
public density
public potential_density
public density_sfc
public density_derivs
public density_delta_z
public density_delta_sfc

private ocean_density_chksum

! interfaces for density 
interface density
   module procedure density_field
   module procedure density_level
   module procedure density_line
   module procedure density_point
end interface

! interfaces for density derivatives 
interface density_derivs
   module procedure density_derivs_field
   module procedure density_derivs_point
end interface

integer :: index_temp
integer :: index_salt
integer :: num_prog_tracers


! nml settings 

! for diagnostic partitioning of vertical according to potential density or theta classes 
integer :: layer_nk = 80       ! # of classes used to compute diagnostics associated with 
                               ! potential density or potential temperature classes. 

! for diagnostic partitioning of vertical according to potential density classes 
real    :: potrho_press    = 2000.0  ! gauge pressure (dbars) for potential density computed for diagnostics
real    :: potrho_min      = 1028.0  ! (kg/m^3)         
real    :: potrho_max      = 1038.0  ! (kg/m^3)           

! for diagnostic partitioning of vertical according to potential temperature classes 
real    :: theta_min = -2.0  ! (degrees C)         
real    :: theta_max = 30.0  ! (degrees C)           

! standard atmospheric pressure (dbar).
! Should be set to 0.0 if assume zero pressure for the 
! overlying atmosphere.  But if have a realistic atmospheric
! pressure loading, then set press_standard=10.1325.  
real :: press_standard = 0.0 

logical :: debug_density         = .false. 
logical :: module_is_initialized = .false.

namelist /ocean_density_nml/ s_test, t_test, p_test, press_standard, &
                             linear_eos, alpha_linear_eos, beta_linear_eos, & 
                             potrho_press, potrho_min, potrho_max, &
                             layer_nk, &
                             theta_min, theta_max, &
                             debug_density
contains

!#######################################################################
! <SUBROUTINE NAME="ocean_density_init">
!
! <DESCRIPTION>
! Initialize the density module
! </DESCRIPTION>
!
  subroutine ocean_density_init (Grid, Domain, Time, Time_steps, Thickness, T_prog, Dens, debug, ens_ocean)

    type(ocean_domain_type), intent(in), target :: Domain
    type(ocean_grid_type), intent(in), target   :: Grid
    type(ocean_time_type), intent(in)           :: Time
    type(ocean_time_steps_type), intent(in)     :: Time_steps
    type(ocean_thickness_type), intent(in)      :: Thickness
    type(ocean_prog_tracer_type), intent(in)    :: T_prog(:)
    type(ocean_density_type), intent(inout)     :: Dens

    logical, intent(in), optional :: debug
    logical, intent(in), optional :: ens_ocean
    
    integer :: ioun, io_status, krho, k, ierr, tau, taup1, n
    real    :: rho_test, alpha_test, beta_test, drho_dtheta_test, drho_dsal_test, diff_test

    integer                         :: id_potrho_bounds, id_potrho_axis
    integer                         :: id_potrho_xt, id_potrho_yt

    integer                         :: id_theta_bounds, id_theta_axis
    integer                         :: id_theta_xt, id_theta_yt

    real, allocatable, dimension(:) :: potrho_bounds
    real, allocatable, dimension(:) :: theta_bounds
    real, allocatable, dimension(:) :: zt_bounds
    real, allocatable, dimension(:) :: zw_bounds     

    real                            :: potrho_interval 
    real                            :: theta_interval 

    character*128 filename

    if ( module_is_initialized ) then 
      call mpp_error(FATAL, '==>Error in ocean_density_mod (ocean_density_init): module already initialized.')
    endif 

    module_is_initialized = .TRUE.

    if (PRESENT(debug)) debug_density = debug

    tau   = Time%tau
    taup1 = Time%taup1
    
    index_temp=-1;index_salt=-1
    num_prog_tracers = size(T_prog(:))
    do n= 1, num_prog_tracers
       if (T_prog(n)%name == 'temp') index_temp = n
       if (T_prog(n)%name == 'salt') index_salt = n
    enddo
       
    call write_version_number( version, tagname )

    ! provide for namelist over-ride of defaults
    ioun = open_namelist_file()
    read (ioun,ocean_density_nml,IOSTAT=io_status)
    write (stdout(),'(/)')
    write (stdout(),ocean_density_nml)
    write (stdlog(),ocean_density_nml)
    ierr = check_nml_error(io_status, 'ocean_density_nml')
    call close_file(ioun)

    if(linear_eos) then
          write (stdout(),'(/1x,a)') ' ==> Note: USING linear EOS designed for idealized Boussinesq simulations.'
          write (stdout(),'(7x,a)') 'It is a linear function of potential temperature and salinity.'
          write (stdout(),'(7x,a)') 'There is no pressure dependence.'
      else
          write (stdout(),'(/1x,a)') ' ==> Note: USING full EOS, as relevant for realistic ocean climate simulations.' 
          write (stdout(),'(1x,a,f12.6,a)') ' Subtracting standard atmosphere of ',press_standard,' dbar for EOS calculation.'  
      endif


    ! 25 coefficients in the realistic equation of state 
    a0  =  9.99843699e+02
    a1  =  7.35212840e+00
    a2  = -5.45928211e-02
    a3  =  3.98476704e-04
    a4  =  2.96938239e+00
    a5  = -7.23268813e-03
    a6  =  2.12382341e-03
    a7  =  1.04004591e-02
    a8  =  1.03970529e-07
    a9  =  5.18761880e-06
    a10 = -3.24041825e-08
    a11 = -1.23869360e-11

    two_a2 = 2.0*a2
    three_a3 = 3.0*a3
    two_a6   = 2.0*a6
    two_a8   = 2.0*a8
    two_a11  = 2.0*a11

    b0  =  1.00000000e+00 
    b1  =  7.28606739e-03
    b2  = -4.60835542e-05 
    b3  =  3.68390573e-07
    b4  =  1.80809186e-10
    b5  =  2.14691708e-03
    b6  = -9.27062484e-06
    b7  = -1.78343643e-10
    b8  =  4.76534122e-06
    b9  =  1.63410736e-09
    b10 =  5.30848875e-06
    b11 = -3.03175128e-16
    b12 = -1.27934137e-17

    two_b2    = 2.0*b2
    three_b3  = 3.0*b3
    four_b4   = 4.0*b4
    three_b7  = 3.0*b7
    onep5_b8  = 1.5*b8
    onep5_b9  = 1.5*b9
    two_b9    = 2.0*b9
    three_b11 = 3.0*b11

#ifndef STATIC_MEMORY
    call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
    nk = Grid%nk
    allocate(Dens%rho(isd:ied,jsd:jed,nk))
    allocate(Dens%rho_taum1(isd:ied,jsd:jed,nk))
    allocate(Dens%rho_taup1(isd:ied,jsd:jed,nk))
    allocate(Dens%potrho(isd:ied,jsd:jed,nk))
    allocate(Dens%pressure_at_depth(isd:ied,jsd:jed,nk))
#endif

    do k=1,nk 
       Dens%rho(:,:,k)       = rho0*Grid%tmask(:,:,k)
       Dens%rho_taum1(:,:,k) = rho0*Grid%tmask(:,:,k)
       Dens%rho_taup1(:,:,k) = rho0*Grid%tmask(:,:,k)
       Dens%potrho(:,:,k)    = rho0*Grid%tmask(:,:,k)
       Dens%pressure_at_depth(:,:,k) = rho0*grav*Thickness%ztp(:,:,k)*c2dbars
    enddo

    Dom  => Domain
    Grd  => Grid
    dtts = Time_steps%dtts
    
    allocate(denr(isd:ied,jsd:jed,nk))
    denr(:,:,:) = 0.0

    filename = 'INPUT/ocean_density.res.nc'    
    if (.NOT.file_exist(trim(filename))) then
       Dens%rho_taum1(:,:,:) = density(T_prog(index_salt)%field(:,:,:,tau),&
                                      T_prog(index_temp)%field(:,:,:,tau),&
                                      Dens%pressure_at_depth(:,:,:))
       Dens%rho(:,:,:)       = Dens%rho_taum1(:,:,:)
       Dens%rho_taup1(:,:,:) = Dens%rho_taum1(:,:,:)
    else
       call read_data(filename,'rho_taum1', Dens%rho_taum1, domain=Dom%domain2d, append_pelist_name=ens_ocean)
       call mpp_update_domains(Dens%rho_taum1,Dom%domain2d)
   endif

    
   ! Test values for EOS 
   if(debug_density) then 

      write (stdout(),'(/,a)') 'EQUATION OF STATE TEST VALUES'
      write (stdout(),'(a,f6.2,a,f6.2,a,f8.2)') 's_test(psu) = ',s_test,', t_test(C) = ',t_test,', p_test(dbar) = ',p_test  
      rho_test = density(s_test,t_test,p_test) 
      write (stdout(),' (a,f6.2,a,f6.2,a,f8.2,a,e22.16,a)') 'rho  (',s_test,',',t_test,',',p_test,') = ',rho_test,' kg/m^3'

      if(.not. linear_eos) then 
         diff_test = rho_test-mjwf_rho
         write (stdout(),' (a,e22.16,a)') 'diff from MJWF = ',diff_test,' kg/m^3'
      endif 

      call density_derivs(rho_test, s_test,t_test,p_test,drho_dtheta_test, drho_dsal_test)

      alpha_test = -drho_dtheta_test/(epsln+rho_test)
      write (stdout(),' (a,f6.2,a,f6.2,a,f8.2,a,e22.16,a)') 'alpha(',s_test,',',t_test,',',p_test,') = ',alpha_test,' 1/C'

      if(.not. linear_eos) then 
         diff_test = alpha_test-mjwf_alpha
         write (stdout(),' (a,e22.16,a)') 'diff from MJWF = ',diff_test,' 1/C'
      endif 

      beta_test  =  drho_dsal_test  /(epsln+rho_test)
      write (stdout(),' (a,f6.2,a,f6.2,a,f8.2,a,e22.16,a)') 'beta (',s_test,',',t_test,',',p_test,') = ',beta_test,' 1/psu'

      if(.not. linear_eos) then 
         diff_test = beta_test-mjwf_beta
         write (stdout(),' (a,e22.16,a)') 'diff from MJWF = ',diff_test,' 1/psu'
      endif

  endif 

  ! define vertical axes according to potential density classes  
  allocate ( Dens%potrho_ref(layer_nk))
  allocate ( potrho_bounds(layer_nk+1))
  potrho_bounds(1) = potrho_min
  potrho_interval  = (potrho_max-potrho_min)/(epsln+layer_nk)
  do k=2,layer_nk+1
    potrho_bounds(k)=potrho_bounds(k-1)+potrho_interval
  enddo 
  do k=1,layer_nk
    Dens%potrho_ref(k)=potrho_bounds(k)+0.5*potrho_interval
  enddo 

  ! define vertical axes according to potential temperature classes. 
  allocate ( Dens%theta_ref(layer_nk))
  allocate ( theta_bounds(layer_nk+1))
  theta_bounds(1)  = theta_min
  theta_interval   = (theta_max-theta_min)/(epsln+layer_nk)
  do k=2,layer_nk+1
    theta_bounds(k) =theta_bounds(k-1) +theta_interval
  enddo 
  do k=1,layer_nk
    Dens%theta_ref(k) =theta_bounds(k) +0.5*theta_interval
  enddo 

! assuming that if the Grd%name = 'ocean' than we will
! be using this diagnostic axis , otherwise it is not used
! This is done to prevent duplicate axis names for multiple
! grid objects
  
  if (trim(Grd%name) == 'ocean') then
    id_potrho_xt = diag_axis_init ('rho_xt_'//trim(Grd%name),Grd%grid_x_t,'degrees_E','x','tcell longitude',set_name='ocean',&
                                  Domain2=Dom%domain2d)
    id_potrho_yt = diag_axis_init ('rho_yt_'//trim(Grd%name),Grd%grid_y_t,'degrees_N','y','tcell latitude',set_name='ocean',&
                                  Domain2=Dom%domain2d)
    id_potrho_bounds = diag_axis_init &
     ( 'potrho_edges',potrho_bounds, 'kg/m^3', 'z','potential density edges',direction=-1, set_name='ocean')
    id_potrho_axis   = diag_axis_init &
     ( 'potrho',Dens%potrho_ref,&
       'kg/m^3', 'z','potential density',edges=id_potrho_bounds,direction=-1,set_name='ocean')

    Dens%potrho_axes = (/ id_potrho_xt, id_potrho_yt, id_potrho_axis /)

    id_theta_xt  = diag_axis_init ('theta_xt_'//trim(Grd%name),Grd%grid_x_t,'degrees_E','x','tcell longitude',set_name='ocean',&
                                  Domain2=Dom%domain2d)
    id_theta_yt  = diag_axis_init ('theta_yt_'//trim(Grd%name),Grd%grid_y_t,'degrees_N','y','tcell latitude',set_name='ocean',&
                                  Domain2=Dom%domain2d)
    id_theta_bounds = diag_axis_init &
     ( 'theta_edges',potrho_bounds, 'C', 'z','potential temperature edges',direction=1, set_name='ocean')
    id_theta_axis   = diag_axis_init &
     ( 'theta',Dens%theta_ref,&
       'C', 'z','potential temperature',edges=id_theta_bounds,direction=1,set_name='ocean')

    Dens%theta_axes = (/ id_theta_xt,  id_theta_yt,  id_theta_axis  /)

  endif


! register diagnostic fields for diag_manager 

  id_press = register_diag_field ('ocean_model', 'press', Grd%tracer_axes(1:3),&
        Time%model_time, 'absolute pressure', 'dbar', missing_value=-10.0, range=(/-10.,6000./))
  id_rho = register_diag_field ('ocean_model', 'rho', Grd%tracer_axes(1:3), &
       Time%model_time, 'in situ density', 'kg/m^3', missing_value=-10.0, range=(/-10.0,1e5/))
  id_pot_rho = register_diag_field ('ocean_model', 'pot_rho', Grd%tracer_axes(1:3), &
       Time%model_time, 'potential density', 'kg/m^3', &
       missing_value=-10.0, range=(/-10.0,1e5/))    
  id_pot_rho_0 = register_diag_field ('ocean_model', 'pot_rho_0', Grd%tracer_axes(1:3), &
       Time%model_time, 'potential density ref to 0dbar', 'kg/m^3', &
       missing_value=-10.0, range=(/-10.0,1e5/))    
  id_pot_rho_et = register_diag_field ('ocean_model', 'pot_rho_et', Grd%tracer_axes_flux_x(1:3), &
       Time%model_time, 'potential density at east face of tracer cell', 'kg/m^3', &
       missing_value=-10.0, range=(/-10.0,1e5/))    
  id_pot_rho_nt = register_diag_field ('ocean_model', 'pot_rho_nt', Grd%tracer_axes_flux_y(1:3), &
       Time%model_time, 'potential density at north face of tracer cell', 'kg/m^3', &
       missing_value=-10.0, range=(/-10.0,1e5/))    
  id_pot_rho_wt = register_diag_field ('ocean_model', 'pot_rho_wt', Grd%tracer_axes_wt(1:3), &
       Time%model_time, 'potential density at wt points', 'kg/m^3', &
       missing_value=-10.0, range=(/-10.0,1e5/))    
  id_drhodtheta = register_diag_field ('ocean_model', 'drhodtheta', Grd%tracer_axes(1:3), Time%model_time, &
         'd(rho)/d(theta)', 'kg/m^3/C', missing_value=-1e10, range=(/-1e10,1e10/))
  id_drhodsalt   = register_diag_field ('ocean_model', 'drhodsalinity', Grd%tracer_axes(1:3), Time%model_time, &
         'd(rho)/d(salinity)', 'kg/m^3/psu', missing_value=-1e10, range=(/-1e10,1e10/))


  end subroutine ocean_density_init
! </SUBROUTINE> NAME="ocean_density_init"


!#######################################################################
! <SUBROUTINE NAME="update_ocean_density">
!
! <DESCRIPTION>
! Compute ocean density and related fields.
! </DESCRIPTION>
!
  subroutine update_ocean_density(Time, Thickness, Temp, Salt, Ext_mode, patm,  Dens)

  type(ocean_time_type), intent(in)            :: Time
  type(ocean_thickness_type), intent(in)       :: Thickness
  type(ocean_prog_tracer_type), intent(in)     :: Temp
  type(ocean_prog_tracer_type), intent(in)     :: Salt
  type(ocean_external_mode_type), intent(in)   :: Ext_mode
  type(ocean_density_type), intent(inout)      :: Dens
  real, intent(in), dimension(isd:ied,jsd:jed) :: patm

  type(time_type)                              :: next_time 
  
  logical :: used
  integer :: i, j, k
  integer :: tau, taum1, taup1
  real    :: volume, volumeip1, volumejp1

  taum1 = Time%taum1
  tau   = Time%tau
  taup1 = Time%taup1  

  ! compute pressure at depth, including pressure from surface height and loading from atmosphere and/or sea ice
  Dens%pressure_at_depth(:,:,:) = &
  pressure_in_dbars(Thickness, Dens%rho_taum1, patm+grav*Dens%rho_taum1(:,:,1)*Ext_mode%eta_t(:,:,tau))

  ! compute in situ density 
  Dens%rho(:,:,:) = density(Salt%field(:,:,:,tau), Temp%field(:,:,:,tau), Dens%pressure_at_depth(:,:,:))

  ! compute potential density for diagnostic purposes 
  Dens%potrho(:,:,:) = potential_density(Salt%field(:,:,:,tau), Temp%field(:,:,:,tau), potrho_press)
  
  if (id_press > 0) &
   used = send_data (id_press, Dens%pressure_at_depth(:,:,:), &
          Time%model_time, rmask=Grd%tmask(:,:,:), &
          is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  if (id_rho > 0) &
   used = send_data (id_rho, Dens%rho(:,:,:), &
          Time%model_time, rmask=Grd%tmask(:,:,:), &
          is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  if (id_pot_rho > 0) & 
    used = send_data (id_pot_rho, Dens%potrho(:,:,:), &
           Time%model_time, rmask=Grd%tmask(:,:,:), &
           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  if (id_pot_rho_0 > 0) then 
    wrk1(:,:,:) = potential_density(Salt%field(:,:,:,tau), Temp%field(:,:,:,tau), 0.0)
    used = send_data (id_pot_rho_0, wrk1(:,:,:), &
           Time%model_time, rmask=Grd%tmask(:,:,:), &
           is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif 

  next_time = increment_time(Time%model_time, int(dtts), 0)
  if (need_data(id_pot_rho_wt,next_time)) then 
      wrk1(:,:,:) = 0.0
      do j=jsc,jec
         do i=isc,iec
            wrk1(i,j,nk) = Dens%potrho(i,j,nk)
            do k=1,nk-1
               if(Grd%tmask(i,j,k) > 0.0) then   
                  wrk1(i,j,k) = ( Grd%tmask(i,j,k+1)*Dens%potrho(i,j,k+1)*Thickness%dht(i,j,k+1,tau) &
                                                    +Dens%potrho(i,j,k)  *Thickness%dht(i,j,k,tau)) &
                       /( epsln + Thickness%dht(i,j,k,tau) + Grd%tmask(i,j,k+1)*Thickness%dht(i,j,k+1,tau) ) 
               endif   
            enddo
         enddo
      enddo
      used = send_data (id_pot_rho_wt, wrk1(:,:,:), &
             Time%model_time, rmask=Grd%tmask(:,:,:), &
             is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif
  if (need_data(id_pot_rho_et,next_time)) then 
      wrk1(:,:,:) = 0.0
      do j=jsc,jec
         do i=isc,iec
            do k=1,nk
              volume      = Thickness%dht(i,j,k,tau)  *Grd%dat(i,j)  *Grd%tmask(i,j,k) 
              volumeip1   = Thickness%dht(i+1,j,k,tau)*Grd%dat(i+1,j)*Grd%tmask(i+1,j,k) 
              wrk1(i,j,k) = ( Dens%potrho(i,j,k)*volume + Dens%potrho(i+1,j,k)*volumeip1 ) &
                            /( epsln + volume + volumeip1 ) 
            enddo
         enddo
      enddo
      used = send_data (id_pot_rho_et, wrk1(:,:,:), &
             Time%model_time, rmask=Grd%tmask(:,:,:), &
             is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif
  if (need_data(id_pot_rho_nt,next_time)) then 
      wrk1(:,:,:) = 0.0
      do j=jsc,jec
         do i=isc,iec
            do k=1,nk
              volume      = Thickness%dht(i,j,k,tau)  *Grd%dat(i,j)  *Grd%tmask(i,j,k) 
              volumejp1   = Thickness%dht(i,j+1,k,tau)*Grd%dat(i,j+1)*Grd%tmask(i,j+1,k) 
              wrk1(i,j,k) = ( Dens%potrho(i,j,k)*volume + Dens%potrho(i,j+1,k)*volumejp1 ) &
                            /( epsln + volume + volumejp1 ) 
            enddo
         enddo
      enddo
      used = send_data (id_pot_rho_nt, wrk1(:,:,:), &
             Time%model_time, rmask=Grd%tmask(:,:,:), &
             is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  endif

end subroutine update_ocean_density
! </SUBROUTINE> NAME="update_ocean_density"


!#######################################################################
! <SUBROUTINE NAME="update_ocean_density_taup1">
!
! <DESCRIPTION>
! Compute ocean density at taup1 
! </DESCRIPTION>
!
  subroutine update_ocean_density_taup1(Time, Thickness, Temp, Salt, Ext_mode, Dens, patm)

  type(ocean_time_type), intent(in)            :: Time
  type(ocean_thickness_type), intent(in)       :: Thickness
  type(ocean_prog_tracer_type), intent(in)     :: Temp
  type(ocean_prog_tracer_type), intent(in)     :: Salt
  type(ocean_external_mode_type), intent(in)   :: Ext_mode
  type(ocean_density_type), intent(inout)      :: Dens
  real, intent(in), dimension(isd:ied,jsd:jed) :: patm
  
  logical :: used
  integer :: i, j, k, tau, taup1

  tau   = Time%tau
  taup1 = Time%taup1  
 
  wrk1(:,:,:) = pressure_in_dbars(Thickness, Dens%rho(:,:,:), patm+grav*Dens%rho(:,:,1)*Ext_mode%eta_t(:,:,taup1))
  Dens%rho_taup1(:,:,:) = density(Salt%field(:,:,:,taup1), Temp%field(:,:,:,taup1), wrk1(:,:,:))

end subroutine update_ocean_density_taup1
! </SUBROUTINE> NAME="update_ocean_density_taup1"



!#######################################################################
! <FUNCTION NAME="density_field">
!
! <DESCRIPTION>
! Compute density for all grid points.  
! Note that pressure here is 
! gauge pressure = absolute pressure - press_standard 
! and salinity is in model units (psu).  
! </DESCRIPTION>
!
  function density_field (salinity, theta, press)

    integer :: i, j, k
    real    :: t1, t2, s1, sp5, p1, p1t1
    real    :: num, den

    real, intent(in), dimension(isd:ied,jsd:jed,nk) :: salinity, theta, press
    real, dimension(isd:ied,jsd:jed,nk)             :: density_field

    if ( .not. module_is_initialized ) then 
      call mpp_error(FATAL, '==>Error in ocean_density_mod (density_field): module must be initialized')
    endif 

    if(linear_eos) then

       do k=1,nk
          do j=jsd,jed
             do i=isd,ied
                density_field(i,j,k) = rho0 - alpha_linear_eos*theta(i,j,k) + beta_linear_eos*salinity(i,j,k)
             enddo
          enddo
       enddo

    else 

       do k=1,nk
          do j=jsd,jed
             do i=isd,ied

                t1  = theta(i,j,k)
                t2  = t1*t1

                s1  = salinity(i,j,k)
                sp5 = sqrt(s1) 

                p1   = press(i,j,k) - press_standard
                p1t1 = p1*t1

                num = a0 + t1*(a1 + t1*(a2+a3*t1) )  &
                     + s1*(a4 + a5*t1  + a6*s1)        &
                     + p1*(a7 + a8*t2 + a9*s1 + p1*(a10+a11*t2))

                den = b0 + t1*(b1 + t1*(b2 + t1*(b3 + t1*b4)))      &
                     + s1*(b5 + t1*(b6 + b7*t2) + sp5*(b8 + b9*t2)) &  
                     + p1*(b10 + p1t1*(b11*t2 + b12*p1))
 
                denr(i,j,k) = 1.0/(epsln+den) 

                density_field(i,j,k) = num*denr(i,j,k)

             enddo
          enddo
       enddo
    endif

  end function density_field
! </FUNCTION> NAME="density_field"


!#######################################################################
! <FUNCTION NAME="density_level">
!
! <DESCRIPTION>
! Compute density at a particular k-level. Note that pressure here is 
! the gauge pressure = absolute pressure - press_standard
! </DESCRIPTION>
!
  function density_level(salinity, theta, press)

    integer :: i, j
    real    :: t1, t2, s1, sp5, p1, p1t1
    real    :: num, den

    real, intent(in), dimension(isd:ied,jsd:jed) :: salinity, theta, press
    real, dimension(isd:ied,jsd:jed)             :: density_level 

    if ( .not. module_is_initialized ) then 
      call mpp_error(FATAL, '==>Error in ocean_density_mod (density_level): module must be initialized')
    endif 

    if(linear_eos) then

       do j=jsd,jed
          do i=isd,ied
             density_level(i,j) = rho0 - alpha_linear_eos*theta(i,j) + beta_linear_eos*salinity(i,j)
          enddo
       enddo

    else 

       do j=jsd,jed
          do i=isd,ied

             t1  = theta(i,j)
             t2  = t1*t1

             s1  = salinity(i,j)
             sp5 = sqrt(s1) 

             p1   = press(i,j) - press_standard
             p1t1 = p1*t1

             num = a0 + t1*(a1 + t1*(a2+a3*t1) )  &
                  + s1*(a4 + a5*t1  + a6*s1)        &
                  + p1*(a7 + a8*t2 + a9*s1 + p1*(a10+a11*t2))

             den = b0 + t1*(b1 + t1*(b2 + t1*(b3 + t1*b4)))      &
                  + s1*(b5 + t1*(b6 + b7*t2) + sp5*(b8 + b9*t2)) &  
                  + p1*(b10 + p1t1*(b11*t2 + b12*p1))

             density_level(i,j) = num/(epsln+den)

          enddo
       enddo
    endif

  end function density_level
! </FUNCTION> NAME="density_level"

!#######################################################################
! <FUNCTION NAME="density_line">
!
! <DESCRIPTION>
! Compute density at a particular k-level and j index.  This scheme
! is used in the vectorized version of the full convection scheme. 
! Note that pressure here is the 
! gauge pressure = absolute pressure - press_standard
! </DESCRIPTION>
!
  function density_line(salinity, theta, press)

    integer :: i, j
    real    :: t1, t2, s1, sp5, p1, p1t1
    real    :: num, den

    real, intent(in), dimension(isd:ied) :: salinity, theta, press
    real, dimension(isd:ied)             :: density_line 

    if ( .not. module_is_initialized ) then 
      call mpp_error(FATAL, '==>Error in ocean_density_mod (density_line): module must be initialized')
    endif 

    if(linear_eos) then

        do i=isd,ied
           density_line(i) = rho0 - alpha_linear_eos*theta(i) + beta_linear_eos*salinity(i)
        enddo

    else 

        do i=isd,ied

           t1  = theta(i)
           t2  = t1*t1

           s1  = salinity(i)
           sp5 = sqrt(s1) 

           p1   = press(i) - press_standard
           p1t1 = p1*t1

           num = a0 + t1*(a1 + t1*(a2+a3*t1) )  &
                + s1*(a4 + a5*t1  + a6*s1)        &
                + p1*(a7 + a8*t2 + a9*s1 + p1*(a10+a11*t2))

           den = b0 + t1*(b1 + t1*(b2 + t1*(b3 + t1*b4)))      &
                  + s1*(b5 + t1*(b6 + b7*t2) + sp5*(b8 + b9*t2)) &  
                  + p1*(b10 + p1t1*(b11*t2 + b12*p1))

           density_line(i) = num/(epsln+den)
       enddo
    endif

  end function density_line
! </FUNCTION> NAME="density_line"

!#######################################################################
! <FUNCTION NAME="potential_density">
!
! <DESCRIPTION>
! Compute potential density referenced to some given gauge pressure. 
! Note that potential density referenced to the surface (i.e., sigma_0)
! has a zero gauge pressure, so pressure=0.0 should be the argument
! to the function. 
! </DESCRIPTION>
!
  function potential_density (salinity, theta, pressure)

    real, intent(in), dimension(isd:ied,jsd:jed,nk) :: salinity, theta
    real, intent(in) :: pressure

    integer :: i, j, k
    real    :: t1, t2, s1, sp5, p1, p1t1
    real    :: num, den
    real, dimension(isd:ied,jsd:jed,nk) :: potential_density

    if ( .not. module_is_initialized ) then 
      call mpp_error(FATAL, '==>Error in ocean_density_mod (potential_density): module must be initialized')
    endif 

    if(linear_eos) then

       do k=1,nk
          do j=jsd,jed
             do i=isd,ied
                potential_density(i,j,k) = rho0 - alpha_linear_eos*theta(i,j,k) + beta_linear_eos*salinity(i,j,k)
             enddo
          enddo
       enddo

    else 

       p1 = pressure

       if(p1 > 0.0) then 

          do k=1,nk
             do j=jsd,jed
                do i=isd,ied
                   t1  = theta(i,j,k)
                   t2  = t1*t1

                   s1  = salinity(i,j,k)
                   sp5 = sqrt(s1) 

                   p1t1 = p1*t1

                   num = a0 + t1*(a1 + t1*(a2+a3*t1) )    &
                        + s1*(a4 + a5*t1  + a6*s1)        &
                        + p1*(a7 + a8*t2 + a9*s1 + p1*(a10+a11*t2))

                   den = b0 + t1*(b1 + t1*(b2 + t1*(b3 + t1*b4)))      &
                        + s1*(b5 + t1*(b6 + b7*t2) + sp5*(b8 + b9*t2)) &  
                        + p1*(b10 + p1t1*(b11*t2 + b12*p1))

                   potential_density(i,j,k) = num/(epsln+den)
                enddo
             enddo
          enddo

       elseif(p1==0.0) then ! for sigma_0

          do k=1,nk
             do j=jsd,jed
                do i=isd,ied
                   t1  = theta(i,j,k)
                   t2  = t1*t1

                   s1  = salinity(i,j,k)
                   sp5 = sqrt(s1) 

                   num = a0 + t1*(a1 + t1*(a2+a3*t1) )  &
                        + s1*(a4 + a5*t1  + a6*s1)        
                        
                   den = b0 + t1*(b1 + t1*(b2 + t1*(b3 + t1*b4)))      &
                        + s1*(b5 + t1*(b6 + b7*t2) + sp5*(b8 + b9*t2)) 

                   potential_density(i,j,k) = num/(epsln+den)
                enddo
             enddo
          enddo

       elseif(p1 < 0.0) then 
          call mpp_error(FATAL,'==>Error in ocean_density_mod: potential density at negative pressure is not defined')

       endif  ! endif for value of pressure 

    endif  ! endif for linear_eos 

  end function potential_density
! </FUNCTION> NAME="potential_density"


!#######################################################################
! <FUNCTION NAME="density_sfc">
!
! <DESCRIPTION>
! Compute density as a function of surface salinity, 
! surface theta, and insitu gauge pressure. 
! For use in KPP mixed layer scheme 
! </DESCRIPTION>
!
  function density_sfc (salinity, theta, press)

    integer :: i, j, k
    real    :: t1, t2, s1, sp5, p1, p1t1
    real    :: num, den

    real, intent(in), dimension(isd:ied,jsd:jed,nk) :: salinity, theta, press
    real, dimension(isd:ied,jsd:jed,nk) :: density_sfc

    if ( .not. module_is_initialized ) then 
       call mpp_error(FATAL, '==>Error in ocean_density_mod (density_sfc): module must be initialized')
    endif 

    if(linear_eos) then

       do k=1,nk
          do j=jsd,jed
             do i=isd,ied
                density_sfc(i,j,k) = rho0 - alpha_linear_eos*theta(i,j,1) + beta_linear_eos*salinity(i,j,1)
             enddo
          enddo
       enddo

    else 

       do k=1,nk
          do j=jsd,jed
             do i=isd,ied

                t1  = theta(i,j,1)
                t2  = t1*t1

                s1  = salinity(i,j,1)
                sp5 = sqrt(s1) 

                p1   = press(i,j,k) - press_standard 
                p1t1 = p1*t1

                num = a0 + t1*(a1 + t1*(a2+a3*t1) )  &
                     + s1*(a4 + a5*t1  + a6*s1)        &
                     + p1*(a7 + a8*t2 + a9*s1 + p1*(a10+a11*t2))

                den = b0 + t1*(b1 + t1*(b2 + t1*(b3 + t1*b4)))      &
                     + s1*(b5 + t1*(b6 + b7*t2) + sp5*(b8 + b9*t2)) &  
                     + p1*(b10 + p1t1*(b11*t2 + b12*p1))

                density_sfc(i,j,k) = num/(epsln + den)

             enddo
          enddo
       enddo

    endif

  end function density_sfc
! </FUNCTION> NAME="density_sfc"


!#######################################################################
! <FUNCTION NAME="density_point">
!
! <DESCRIPTION>
! Compute density at a single model grid point. Note that pressure here  
! is the gauge pressure = absolute pressure - press_standard 
! </DESCRIPTION>
!
  function density_point (s1, t1, p1_dbars)

    real, intent(in) :: s1, t1, p1_dbars
    real             :: t2, sp5, p1, p1t1
    real             :: num, den
    real             :: density_point

    if ( .not. module_is_initialized ) then 
      call mpp_error(FATAL, '==>Error in ocean_density_mod (density_point): module must be initialized')
    endif 

    if(linear_eos) then

       density_point = rho0 - alpha_linear_eos*t1 + beta_linear_eos*s1

    else 

       t2  = t1*t1
       sp5 = sqrt(s1) 

       p1   = p1_dbars - press_standard 
       p1t1 = p1*t1

       num = a0 + t1*(a1 + t1*(a2+a3*t1))                     &
            + s1*(a4 + a5*t1  + a6*s1)                        &
            + p1*(a7 + a8*t2 + a9*s1 + p1*(a10+a11*t2))

       den = b0 + t1*(b1 + t1*(b2 + t1*(b3 + t1*b4 )))      &
            + s1*(b5 + t1*(b6 + b7*t2) + sp5*(b8 + b9*t2))  &  
            + p1*(b10 + p1t1*(b11*t2 + b12*p1))

       density_point = num/(epsln+den)

    endif

  end function density_point
! </FUNCTION> NAME="density_point"


!#######################################################################
! <SUBROUTINE NAME="density_derivs_field">
!
! <DESCRIPTION>
! Compute partial derivative of density with respect to potential 
! temperature and with respect to salinity.  Hold pressure constant.  
! Pressure here is gauge pressure = absolute press - press_standard
! </DESCRIPTION>
!
  subroutine density_derivs_field (rho, salinity, theta, press, Time, density_theta, density_salinity)

    type(ocean_time_type), intent(in) :: Time
    real, intent(in),  dimension(isd:ied,jsd:jed,nk) :: rho
    real, intent(in),  dimension(isd:ied,jsd:jed,nk) :: salinity
    real, intent(in),  dimension(isd:ied,jsd:jed,nk) :: theta
    real, intent(in),  dimension(isd:ied,jsd:jed,nk) :: press
    real, intent(out), dimension(isd:ied,jsd:jed,nk) :: density_theta
    real, intent(out), dimension(isd:ied,jsd:jed,nk) :: density_salinity

    integer :: i, j, k
    real :: t1, t2, s1, sp5, p1, p1t1
    real :: dnum_dtheta,    dden_dtheta 
    real :: dnum_dsalinity, dden_dsalinity

    if ( .not. module_is_initialized ) then 
      call mpp_error(FATAL, '==>Error in ocean_density_mod (density_derivs_field): module must be initialized')
    endif 

    if(linear_eos) then

       do k=1,nk
          do j=jsd,jed
             do i=isd,ied
                density_theta(i,j,k)    = -alpha_linear_eos 
                density_salinity(i,j,k) =  beta_linear_eos
             enddo
          enddo
       enddo

    else 

       do k=1,nk
          do j=jsd,jed
             do i=isd,ied

                t1  = theta(i,j,k)
                t2  = t1*t1
                s1  = salinity(i,j,k)
                sp5 = sqrt(s1) 

                p1   = press(i,j,k) - press_standard 
                p1t1 = p1*t1

                dnum_dtheta = a1 + t1*(two_a2 + three_a3*t1) &
                     + a5*s1                                 &
                     + p1t1*(two_a8 + two_a11*p1)    
                dden_dtheta = b1 + t1*(two_b2 + t1*(three_b3 + four_b4*t1)) &
                     + s1*(b6 + t1*(three_b7*t1 + two_b9*sp5))              &
                     + p1*p1*(three_b11*t2 + b12*p1)

                dnum_dsalinity = a4 + a5*t1 + two_a6*s1 + a9*p1
                dden_dsalinity = b5 + t1*(b6 + b7*t2) + sp5*(onep5_b8 + onep5_b9*t2)

                density_theta(i,j,k)    = denr(i,j,k)*(dnum_dtheta    - rho(i,j,k)*dden_dtheta)
                density_salinity(i,j,k) = denr(i,j,k)*(dnum_dsalinity - rho(i,j,k)*dden_dsalinity)

             enddo
          enddo
       enddo

    endif

    if (id_drhodtheta > 0) used = send_data (id_drhodtheta, density_theta(:,:,:), &
                                  Time%model_time, rmask=Grd%tmask(:,:,:), &
                                  is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
    if (id_drhodsalt  > 0) used = send_data (id_drhodsalt, density_salinity(:,:,:), &
                                  Time%model_time, rmask=Grd%tmask(:,:,:), &
                                  is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)


  end subroutine density_derivs_field
! </SUBROUTINE> NAME="density_derivs_field"


!#######################################################################
! <SUBROUTINE NAME="density_derivs_point">
!
! <DESCRIPTION>
! Compute partial derivative of density with respect to potential 
! temperature and with respect to salinity.  Do so here for a point. 
! Pressure here is gauge pressure = absolute pressure - press_standard 
! </DESCRIPTION>
!
  subroutine density_derivs_point (rho, salinity, theta, press, density_theta, density_salinity)

    real, intent(in)  :: rho, salinity, theta, press
    real, intent(out) :: density_theta, density_salinity
    integer :: i, j, k
    real    :: t1, t2, s1, sp5, p1, p1t1
    real    :: dnum_dtheta,    dden_dtheta 
    real    :: dnum_dsalinity, dden_dsalinity
    real    :: den, denrecip

    if ( .not. module_is_initialized ) then 
       call mpp_error(FATAL, '==>Error in ocean_density_mod (density_derivs_point): module must be initialized')
    endif 

    if(linear_eos) then

       density_theta    = -alpha_linear_eos 
       density_salinity =  beta_linear_eos

    else 

       t1  = theta
       t2  = t1*t1
       s1  = salinity
       sp5 = sqrt(s1) 

       p1   = press - press_standard 
       p1t1 = p1*t1

       den = b0 + t1*(b1 + t1*(b2 + t1*(b3 + t1*b4 )))      &
            + s1*(b5 + t1*(b6 + b7*t2) + sp5*(b8 + b9*t2))  &  
            + p1*(b10 + p1t1*(b11*t2 + b12*p1))

       denrecip = 1.0/(den + epsln)

       dnum_dtheta = a1 + t1*(2.0*a2 + 3.0*a3*t1) &
            + a5*s1                                   &
            + 2.0*p1t1*(a8 + a11*p1)    
       dden_dtheta = b1 + t1*(2.0*b2 + t1*(3.0*b3 + 4.0*b4*t1)) &
            + s1*(b6 + t1*(3.0*b7*t1 + 2.0*b9*sp5))               &
            + p1*p1*(3.0*b11*t2 + b12*p1)

       dnum_dsalinity = a4 + a5*t1 + two_a6*s1 + a9*p1
       dden_dsalinity = b5 + t1*(b6 + b7*t2) + sp5*(onep5_b8 + onep5_b9*t2)
 
       density_theta    = denrecip*(dnum_dtheta    - rho*dden_dtheta)
       density_salinity = denrecip*(dnum_dsalinity - rho*dden_dsalinity)

    endif

  end subroutine density_derivs_point
! </SUBROUTINE> NAME="density_derivs_point"


!#######################################################################
! <FUNCTION NAME="density_delta_z">
!
! <DESCRIPTION>
! rho(k)-rho(k+1) for all i,j with both temperatures referenced to the 
! deeper pressure depth.  Of use for KPP scheme. 
! </DESCRIPTION>
!
function density_delta_z (rho_initial, salinity, theta, press)

  real, intent(in), dimension(isd:ied,jsd:jed,nk) :: rho_initial, salinity, theta, press
  real, dimension(isd:ied,jsd:jed,nk) :: density_delta_z
  integer :: k

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, '==>Error in ocean_density_mod (density_delta_z): module must be initialized')
  endif 

  do k=1,nk
    wrk1(:,:,k) = press(:,:,min(k+1,nk))
  enddo
  wrk2(:,:,:) = density(salinity(:,:,:), theta(:,:,:), wrk1(:,:,:))
  do k=1,nk
    density_delta_z(:,:,k) = (wrk2(:,:,k) - rho_initial(:,:,min(k+1,nk)))*Grd%tmask(:,:,min(k+1,nk))
  enddo

end function density_delta_z
! </FUNCTION> NAME="density_delta_z"



!#######################################################################
! <FUNCTION NAME="density_delta_sfc">
!
! <DESCRIPTION>
! rho(1)-rho(k+1) for all i,j. Of use for KPP scheme. 
! </DESCRIPTION>
!
function density_delta_sfc (rho_initial, salinity, theta, press)

  real, intent(in), dimension(isd:ied,jsd:jed,nk) :: rho_initial, salinity, theta, press
  real, dimension(isd:ied,jsd:jed,nk) :: density_delta_sfc
  integer :: k

  if ( .not. module_is_initialized ) then 
    call mpp_error(FATAL, '==>Error in ocean_density_mod (density_delta_sfc): module must be initialized')
  endif 

  do k=1,nk
    wrk1(:,:,k) = press(:,:,min(k+1,nk))
  enddo
  wrk2(:,:,:) = density_sfc(salinity, theta, wrk1)
  do k=1,nk
    density_delta_sfc(:,:,k) = (wrk2(:,:,k) - rho_initial(:,:,min(k+1,nk)))*Grd%tmask(:,:,min(k+1,nk))
  enddo

end function density_delta_sfc
! </FUNCTION> NAME="density_delta_sfc"


!#######################################################################
! <SUBROUTINE NAME="ocean_density_rstrt">
!
! <DESCRIPTION>
! Write density field to a restart.
! </DESCRIPTION>
subroutine ocean_density_rstrt(Time, Dens, ens_ocean)

  type(ocean_time_type), intent(in)    :: Time
  type(ocean_density_type), intent(in) :: Dens
  logical, intent(in), optional        :: ens_ocean

  integer            :: yr, mon, day, hr, min, sec
  character(len=128) :: filename
  character(len=10)  :: rdte

  if ( .not. module_is_initialized ) then
     call mpp_error(FATAL, '==>Error in ocean_density_mod (ocean_density_rstrt): module must be initialized')
  endif

  call get_date(Time%model_time,yr,mon,day,hr,min,sec)
  write(rdte,'(i4,3i2.2)') yr,mon,day,hr
  filename = 'IRESTART/'//rdte//'ocean_density.res'
  call write_data(filename,'rho_taum1', Dens%rho_taum1, domain=Dom%domain2d, append_pelist_name=ens_ocean)

  write(stdout(),*) ' '
  write(stdout(),*) '=== Ending density checksum follows ==='
  call ocean_density_chksum(Time, Dens)

end subroutine ocean_density_rstrt
! </SUBROUTINE> NAME="ocean_density_rstrt"


!#######################################################################
! <SUBROUTINE NAME="ocean_density_end">
!
! <DESCRIPTION>
! Write density field to a restart.
! </DESCRIPTION>
subroutine ocean_density_end(Time, Dens, ens_ocean)

  type(ocean_time_type), intent(in)    :: Time
  type(ocean_density_type), intent(in) :: Dens
  logical, intent(in), optional        :: ens_ocean
  
  character*128 filename

  if ( .not. module_is_initialized ) then 
     call mpp_error(FATAL, '==>Error in ocean_density_mod (ocean_density_end): module must be initialized')
  endif 

  filename = 'RESTART/ocean_density.res'
  call write_data(filename,'rho_taum1', Dens%rho_taum1, domain=Dom%domain2d, append_pelist_name=ens_ocean)

  write(stdout(),*) ' '
  write(stdout(),*) '=== Ending density checksum follows ==='
  call ocean_density_chksum(Time, Dens)

  module_is_initialized = .FALSE.

end subroutine ocean_density_end
! </SUBROUTINE> NAME="ocean_density_end"


!#######################################################################
! <SUBROUTINE NAME="ocean_density_chksum">
!
! <DESCRIPTION>
! Compute checksums for density. 
! </DESCRIPTION>
subroutine ocean_density_chksum(Time, Dens, chksum)

  type(ocean_time_type), intent(in)         :: Time
  type(ocean_density_type), intent(in)      :: Dens
  integer(i8_kind), optional, intent(inout) :: chksum
  integer(i8_kind)                          :: chk_sum
  
  if (.not. module_is_initialized) then 
    call mpp_error(FATAL,'==>Error from ocean_density_mod (ocean_density_chksum): module not yet initialized ')
  endif 

  call write_timestamp(Time%model_time)

  wrk1(isc:iec,jsc:jec,:) = Dens%rho_taum1(isc:iec,jsc:jec,:)*Grd%tmask(isc:iec,jsc:jec,:)
  chk_sum                 = mpp_chksum(wrk1(isc:iec,jsc:jec,:))

  write(stdout(),*) 'Density rho_taum1 chksum = ',  chk_sum

end subroutine ocean_density_chksum
! </SUBROUTINE>  NAME="ocean_density_chksum"


end module ocean_density_mod





