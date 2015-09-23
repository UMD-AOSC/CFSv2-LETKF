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
#include <fms_platform.h>

!
! 
!<CONTACT EMAIL="Richard.Slater@noaa.gov"> Richard D. Slater
!</CONTACT>
!
!<REVIEWER EMAIL="John.Dunne@noaa.gov"> John P. Dunne
!</REVIEWER>
!
!<OVERVIEW>
! Ocean Carbon Model Intercomparison Study II: Biotic module
!</OVERVIEW>
!
!<DESCRIPTION>
!       Implementation of routines to solve the OCMIP-2 Biotic
!       simulations as outlined in the Biotic-HOWTO documentation,
!       revision 1.7, 1999/10/05.
!</DESCRIPTION>
!
! <INFO>
! <REFERENCE>
! http://www.ipsl.jussieu.fr/OCMIP/phase2/simulations/Biotic/HOWTO-Biotic.html
! </REFERENCE>
!
! <REFERENCE>
! Press, W. H., S. A. Teukosky, W. T. Vetterling, B. P. Flannery, 1992. 
! Numerical Recipes in FORTRAN, Second Edition, Cambridge University Press. 
! </REFERENCE>
!
! <REFERENCE>
! Enting, I.G., T. M. L. Wigley, M. Heimann, 1994. Future Emissions 
! and concentrations of carbon dioxide: key ocean / atmosphere / 
! land analyses, CSIRO Aust. Div. Atmos. Res. Tech. Pap. No. 31, 
! 118 pp.
! </REFERENCE>
! </INFO>
!
! $Id$
!

!
!------------------------------------------------------------------
!
!       Module ocmip2_biotic_mod
!
!       Implementation of routines to solve the OCMIP-2 Biotic
!       simulations as outlined in the Biotic-HOWTO documentation,
!       revision 1.7, 1999/10/05.
!
!------------------------------------------------------------------
!

module  ocmip2_biotic_mod  !{

!
!------------------------------------------------------------------
!
!       Global definitions
!
!------------------------------------------------------------------
!

!
!----------------------------------------------------------------------
!
!       Modules
!
!----------------------------------------------------------------------
!

use diag_manager_mod,         only: send_data
use field_manager_mod,        only: fm_field_name_len, fm_path_name_len, fm_string_len
use field_manager_mod,        only: fm_get_length, fm_get_value, fm_new_value
use fms_mod,                  only: write_data
use mpp_mod,                  only: stdout, stdlog, mpp_error, mpp_sum, FATAL
use time_manager_mod,         only: get_date
use time_interp_external_mod, only: time_interp_external

use ocean_tpm_util_mod, only: otpm_set_tracer_package, otpm_set_prog_tracer
use ocean_tpm_util_mod, only: otpm_check_for_bad_fields, otpm_set_value
use ocean_tpm_util_mod, only: otpm_get_string, otpm_get_logical, otpm_get_integer, otpm_get_real
use ocean_tpm_util_mod, only: otpm_get_logical_array, otpm_get_real_array, otpm_get_string_array
use ocean_tpm_util_mod, only: otpm_start_namelist, otpm_end_namelist
use ocean_tpm_util_mod, only: domain, grid, time, dtts
use ocean_tpm_util_mod, only: isc, iec, jsc, jec, nk, isd, ied, jsd, jed 
use ocean_tpm_util_mod, only: taum1, tau, taup1 
use ocean_tpm_util_mod, only: t_prog, t_diag
use ocean_tpm_util_mod, only: indsal, indtemp
use ocean_tpm_util_mod, only: end_of_year, end_of_month

use ocean_sbc_mod,      only: use_waterflux
use ocean_types_mod,    only: ocean_thickness_type

!
!----------------------------------------------------------------------
!
!       force all variables to be "typed"
!
!----------------------------------------------------------------------
!

implicit none

!
!----------------------------------------------------------------------
!
!       Make all routines and variables private by default
!
!----------------------------------------------------------------------
!

private

!
!----------------------------------------------------------------------
!
!       Public routines
!
!----------------------------------------------------------------------
!

public  :: ocmip2_biotic_bbc
public  :: ocmip2_biotic_end
public  :: ocmip2_biotic_init
public  :: ocmip2_biotic_sbc
public  :: ocmip2_biotic_source
public  :: ocmip2_biotic_start
public  :: ocmip2_biotic_tracer

!
!----------------------------------------------------------------------
!
!       Private routines
!
!----------------------------------------------------------------------
!

private :: allocate_arrays
private :: bc_interp
private :: locate
private :: read_c14atm
private :: read_co2atm
private :: set_array

!public :: km_c
!public :: compensation_depth
!public :: ind_biotic_dop
!public :: pert_time
!public :: pert_starttime
!public :: pert_starttime_mod
!public :: bc_ptr_prev_mo
!public :: bc_ptr_next_mo
!public :: bc_frac_prev_mo
!public :: bc_frac_next_mo
!public :: biotic_flux_dop_a
!public :: biotic_flux_dop_d
!public :: impvd_dop
!public :: conv_dop

!
!----------------------------------------------------------------------
!
!       Private parameters
!
!----------------------------------------------------------------------
!

character(len=fm_field_name_len), parameter     :: package_name = 'ocmip2_biotic'
character(len=48), parameter                    :: mod_name = 'ocmip2_biotic_mod'
character(len=fm_string_len), parameter         :: default_file_in = 'INPUT/ocmip2_biotic.res.nc'
character(len=fm_string_len), parameter         :: default_file_out = 'RESTART/ocmip2_biotic.res.nc'

!
!----------------------------------------------------------------------
!
!       Private types
!
!----------------------------------------------------------------------
!
 
!
!  pco2atm_const        = constant atmospheric CO2 concentration (ppm)
!  sio2_const           = SiO2 concentration (mol/m^3)
!  dic_global           = global annual surface mean DIC concentration
!  dic_global_wrk       = work variable used in calculation of
!                         dic_global
!  alk_global           = global annual surface mean ALK concentration
!  alk_global_wrk       = work variable used in calculation of
!                         alk_global
!  sal_global           = surface global annual mean salinity
!                         concentration (PSU)
!  sal_global_wrk       = work variable used in calculation of
!                         sal_global
!  do_dic_virtual_flux  = true to compute virtual flux for DIC
!  do_alk_virtual_flux  = true to compute virtual flux for ALK
!  do_po4_virtual_flux  = true to compute virtual flux for PO4
!  do_dop_virtual_flux  = true to compute virtual flux for DOP
!  do_o2_virtual_flux   = true to compute virtual flux for O2
!  po4_global           = global annual surface mean PO4 concentration
!  po4_global_wrk       = work variable used in calculation of
!                         po4_global
!  dop_global           = global annual surface mean DOP concentration
!  dop_global_wrk       = work variable used in calculation of
!                         dop_global
!  o2_global            = global annual surface mean O2 concentration
!  o2_global_wrk        = work variable used in calculation of
!                         o2_global
!  atm_co2_flux         : global flux of CO_2 into the atmosphere
!                         (mol/cm^2/s) at each time-step
!  pco2atm_global       : atmospheric concentration of CO_2 (uatm)
!  total_atm_co2(3)     : total mass of CO_2 in the atmosphere (g)
!                         for the three time levels
!  add_phosphate        : if true, then add sufficient PO4 to keep
!                         the predicted PO4 the same as if no depletion
!                         or changed uptake rate were in effect
!  global_atmosphere    : true if use global atmosphere,
!                         false, otherwise
!

type mask_region_type
  real, dimension(:,:,:), pointer       :: mask => NULL()
  real, dimension(:), pointer           :: elon => NULL()
  real, dimension(:), pointer           :: nlat => NULL()
  real, dimension(:), pointer           :: slat => NULL()
  real, dimension(:), pointer           :: wlon => NULL()
  logical                               :: coastal_only
  real                                  :: factor
  logical, dimension(:), pointer        :: t_mask => NULL()
end type mask_region_type

type biotic_type  !{

  logical                               :: soft_tissue_pump
  logical                               :: add_phosphate
  real                                  :: alk_global = 0.0
  real                                  :: alk_global_wrk = 0.0
  real                                  :: atm_co2_flux
  real                                  :: bio_tau
  real                                  :: c_2_p
  real                                  :: ca_remin_depth
  real                                  :: caco3_2_c
  real                                  :: compensation_depth
  real, _ALLOCATABLE, dimension(:,:)         :: csat _NULL
  real, _ALLOCATABLE, dimension(:,:)         :: csat_csurf _NULL
  real, _ALLOCATABLE, dimension(:,:)         :: csurf _NULL
  real                                  :: dic_global = 0.0
  real                                  :: dic_global_wrk = 0.0
  logical                               :: do_alk_virtual_flux
  logical                               :: do_dic_virtual_flux
  logical                               :: do_dop_virtual_flux
  logical                               :: do_o2_virtual_flux
  logical                               :: do_po4_virtual_flux
  real                                  :: dop_global = 0.0
  real                                  :: dop_global_wrk = 0.0
  real, _ALLOCATABLE, dimension(:,:)         :: dpco2 _NULL
  real, _ALLOCATABLE, dimension(:,:,:)       :: fc _NULL
  real, _ALLOCATABLE, dimension(:,:,:)       :: fca _NULL
  character(len=fm_string_len)          :: file_in
  character(len=fm_string_len)          :: file_out
  real, _ALLOCATABLE, dimension(:,:)         :: flux_caco3 _NULL
  real, _ALLOCATABLE, dimension(:,:)         :: flux_poc _NULL
  logical                               :: global_atmosphere
  real                                  :: global_wrk_duration = 0.0
  real, _ALLOCATABLE, dimension(:,:,:)       :: htotal _NULL
  integer                               :: id_csat = -1
  integer                               :: id_csurf = -1
  integer                               :: id_dpco2 = -1
  integer                               :: id_fc = -1
  integer                               :: id_fca = -1
  integer                               :: id_flux_caco3 = -1
  integer                               :: id_flux_poc = -1
  integer                               :: id_htotal = -1
  integer                               :: id_jca = -1
  integer                               :: id_jdop = -1
  integer                               :: id_jo2 = -1
  integer                               :: id_jpo4 = -1
  integer                               :: id_jpo4_add = -1
  integer                               :: id_jprod = -1
  integer                               :: id_pco2atm = -1
  integer                               :: id_pco2surf = -1
  integer                               :: id_vstf_alk = -1
  integer                               :: id_vstf_dic = -1
  integer                               :: id_vstf_dop = -1
  integer                               :: id_vstf_o2 = -1
  integer                               :: id_vstf_po4 = -1
  integer                               :: ind_alk
  integer                               :: ind_dic
  integer                               :: ind_dop
  integer                               :: ind_o2
  integer                               :: ind_po4
  logical                               :: init
  logical                               :: init_global_atmosphere
  real, _ALLOCATABLE, dimension(:,:,:)       :: jca _NULL
  real, _ALLOCATABLE, dimension(:,:,:)       :: jdop _NULL
  real, _ALLOCATABLE, dimension(:,:,:)       :: jo2 _NULL
  real, _ALLOCATABLE, dimension(:,:,:)       :: jpo4 _NULL
  real, _ALLOCATABLE, dimension(:,:,:)       :: jpo4_add _NULL
  real, _ALLOCATABLE, dimension(:,:,:)       :: jprod _NULL
  real, _ALLOCATABLE, dimension(:,:,:)       :: jprod_norm _NULL
  real                                  :: kappa
  integer                               :: km_c
  real                                  :: martin_coeff
  real                                  :: martin_coeff_alt
  real                                  :: n_2_p
  character(len=fm_field_name_len)      :: name
  real                                  :: o2_global = 0.0
  real                                  :: o2_global_wrk = 0.0
  real                                  :: o2_min
  real                                  :: o_2_p
  real, _ALLOCATABLE, dimension(:,:)         :: pco2atm _NULL
  real                                  :: pco2atm_const
  real, _ALLOCATABLE, dimension(:,:)         :: pco2surf _NULL
  real                                  :: pert_start_year              ! this is real because read_data does not support integer
  real                                  :: pert_start_year_model        ! this is real because read_data does not support integer
  real                                  :: pert_time
  real                                  :: po4_global = 0.0
  real                                  :: po4_global_wrk = 0.0
  real                                  :: r_bio_tau
  type(mask_region_type)                :: r_bio_tau_a
  type(mask_region_type)                :: nut_depl
  type(mask_region_type)                :: norm_remin
  type(mask_region_type)                :: no_caco3
  real                                  :: sal_global = 0.0
  real                                  :: sal_global_wrk = 0.0
  real                                  :: sigma
  real, _ALLOCATABLE, dimension(:,:,:)       :: sio2 _NULL
  real                                  :: sio2_const
  real                                  :: stp_alkalinity
  real                                  :: stp_salinity
  real                                  :: stp_temperature
  real, dimension(3)                    :: total_atm_co2
  real, _ALLOCATABLE, dimension(:,:)         :: vstf_alk _NULL
  real, _ALLOCATABLE, dimension(:,:)         :: vstf_dic _NULL
  real, _ALLOCATABLE, dimension(:,:)         :: vstf_dop _NULL
  real, _ALLOCATABLE, dimension(:,:)         :: vstf_o2 _NULL
  real, _ALLOCATABLE, dimension(:,:)         :: vstf_po4 _NULL
  real, _ALLOCATABLE, dimension(:)           :: zfca _NULL
  real, _ALLOCATABLE, dimension(:)           :: zforg _NULL
  real, _ALLOCATABLE, dimension(:)           :: zforg_alt _NULL

end type biotic_type  !}

!
!----------------------------------------------------------------------
!
!       Public variables
!
!----------------------------------------------------------------------
!

logical, public :: do_ocmip2_biotic

!
!----------------------------------------------------------------------
!
!       Private variables
!
!----------------------------------------------------------------------
!

!logical                                :: no_sbc
!logical                                :: no_source
integer                                 :: package_index

!
!----------------------------------------------------------------------
!
!       Input parameters:
!
!  atmco2_t             = atmospheric CO2 concentration from data (ppm)
!  atmc14_t             = atmospheric 14C/12C ratio from data (permil)
!  htotal_in            = default value for htotal for an initial run
!  htotal_scale_lo      = scaling parameter to chose htotallo
!  htotal_scale_hi      = scaling parameter to chose htotalhi
!
!----------------------------------------------------------------------
!

real, dimension(3)                      :: atmc14_t
real                                    :: atmco2_t
real                                    :: htotal_in
real, allocatable, dimension(:,:)       :: htotal_scale_hi
real                                    :: htotal_scale_hi_in
real, allocatable, dimension(:,:)       :: htotal_scale_lo
real                                    :: htotal_scale_lo_in

!
!----------------------------------------------------------------------
!
!       Calculated parameters (with possible initial input values):
!
!  global_wrk_duration  = total time during calculation of global
!                         variables
!  equilibrium          = T if this is an equilibrium run
!  control              = T if this is a control run (one of the
!                         scenarios below, but with constant atmospheric
!                         values)
!  historical           = T if this is an historical run
!  scenario             = character string indicating the scenario
!                         to be performed in this run. Requires that
!                         a perturbation run has been started
!
!----------------------------------------------------------------------
!

logical                         :: control
logical                         :: equilibrium
logical                         :: historical
character(len=32)               :: scenario

!
!----------------------------------------------------------------------
!
!       Global atmosphere model
!
!  total_atm_mass       : total mass (g) of atmosphere
!  atm_mol_wgt          : atmospheric molecular weight (g)
!  co2_mol_wgt          : CO_2 molecular weight (g)
!
!----------------------------------------------------------------------
!

real, parameter                 :: atm_mol_wgt = 28.964 ! g
real, parameter                 :: co2_mol_wgt = 12.0 ! g
real, parameter                 :: total_atm_mass = 5.119e+21 ! g

!
!     VARIABLES:
!        fice        : sea ice fraction (0..1) 
!        xkw         : defined as Xconv * 0.337 * (u^2 + v), where Xconv
!                      is a conversion factor, u is the instantaneous
!                      SSMI winds, and v is the variance of the winds. 
!                      units are [cm/s]
!        patm        : atmospheric pressure. units are [atm]
!
!        biotic_cumul   : cumulative input of BIOTIC_DIC [mol/m2]
!
!        seaicefract_file  : bc file containing sea-ice fraction
!        pistonveloc_file  : bc file containing xkw 
!        atmpress_file     : bc file containing atm. pressure
!        po4_star_file     : bc file containing PO4 data
!  
!        sc_co2            : Schmidt number for the CO2 (1)  
!        kw_co2            : Piston velocity for the CO2 [cm/s]
!        sc_o2             : Schmidt number for the O2 (1)  
!        kw_o2             : Piston velocity for the O2 [cm/s]
!
!                            
!      
! =====================================================================
!

integer                                 :: atmpress_id
character*128                           :: atmpress_file    
character*32                            :: atmpress_name    
!real, allocatable, dimension(:,:,:)    :: biotic_cumul
real, allocatable, dimension(:,:)       :: fice_t
integer                                 :: id_kw_co2 = -1
integer                                 :: id_kw_o2 = -1
integer                                 :: id_o2_sat = -1
integer                                 :: id_sc_co2 = -1
integer                                 :: id_sc_o2 = -1
integer                                 :: km_c_max
real, allocatable, dimension(:,:)       :: kw_co2    
real, allocatable, dimension(:,:)       :: kw_o2    
real, allocatable, dimension(:,:)       :: patm_t
integer                                 :: pistonveloc_id
integer                                 :: po4_star_id
real, allocatable, dimension(:,:,:)     :: po4_star_t
integer                                 :: seaicefract_id
character*128                           :: pistonveloc_file
character*32                            :: pistonveloc_name
character*128                           :: po4_star_file    
character*32                            :: po4_star_name    
real, allocatable, dimension(:,:)       :: sc_co2
real, allocatable, dimension(:,:)       :: sc_o2
character*128                           :: seaicefract_file
character*32                            :: seaicefract_name
real, allocatable, dimension(:,:)       :: xkw_t
type(biotic_type), allocatable, dimension(:)    :: biotic
real, allocatable, dimension(:,:,:)             :: htotalhi
real, allocatable, dimension(:,:,:)             :: htotallo
integer                                         :: instances
real, allocatable, dimension(:,:)               :: o2_saturation
real, allocatable, dimension(:)                 :: tk
real, allocatable, dimension(:)                 :: ts
real, allocatable, dimension(:)                 :: ts2
real, allocatable, dimension(:)                 :: ts3
real, allocatable, dimension(:)                 :: ts4
real, allocatable, dimension(:)                 :: ts5
real, allocatable, dimension(:)                 :: tt



!
!-----------------------------------------------------------------------
!
!       Subroutine and function definitions
!
!-----------------------------------------------------------------------
!

contains

!#######################################################################
! <SUBROUTINE NAME="allocate_arrays">
!
! <DESCRIPTION>
!     Dynamically allocate arrays
! </DESCRIPTION>
!

subroutine allocate_arrays  !{

!
!       arguments
!

!
!       local variables
!

integer :: m
integer :: k
integer :: j
integer :: i
integer :: n

!
!-----------------------------------------------------------------------
!     start executable code
!-----------------------------------------------------------------------
!     

!
!       global variables
!

allocate( xkw_t(isd:ied,jsd:jed) )
allocate( patm_t(isd:ied,jsd:jed) )
allocate( fice_t(isd:ied,jsd:jed) )

allocate( sc_co2(isc:iec,jsc:jec) )
allocate( kw_co2(isc:iec,jsc:jec) )
allocate( sc_o2(isc:iec,jsc:jec) )
allocate( kw_o2(isc:iec,jsc:jec) )
allocate( htotal_scale_lo(isc:iec,jsc:jec) )
allocate( htotal_scale_hi(isc:iec,jsc:jec) )
allocate( htotallo(isc:iec,jsc:jec,1) )
allocate( htotalhi(isc:iec,jsc:jec,1) )
allocate( o2_saturation(isc:iec,jsc:jec) )
allocate( tt(isc:iec) )
allocate( tk(isc:iec) )
allocate( ts(isc:iec) )
allocate( ts2(isc:iec) )
allocate( ts3(isc:iec) )
allocate( ts4(isc:iec) )
allocate( ts5(isc:iec) )
!allocate( po4_star_t(isd:ied,jsd:jed,km_c_max) )
!       this should be dimensioned as above, but the time_interp routine
!       requires that the array dimensions match the datasets dimensions
allocate( po4_star_t(isd:ied,jsd:jed,nk) )

!
!       initialize some arrays
!

xkw_t(:,:) = 0.0
patm_t(:,:) = 0.0
fice_t(:,:) = 0.0

sc_co2(:,:) = 0.0
kw_co2(:,:) = 0.0
sc_o2(:,:) = 0.0
kw_o2(:,:) = 0.0
htotal_scale_lo(:,:) = 0.0
htotal_scale_hi(:,:) = 0.0
htotallo(:,:,:) = 0.0
htotalhi(:,:,:) = 0.0
o2_saturation(:,:) = 0.0
tt(:) = 0.0
tk(:) = 0.0
ts(:) = 0.0
ts2(:) = 0.0
ts3(:) = 0.0
ts4(:) = 0.0
ts5(:) = 0.0
po4_star_t(:,:,:) = 0.0

!
!       allocate biotic array elements
!

do n = 1, instances  !{

  allocate( biotic(n)%pco2atm(isc:iec,jsc:jec) )
  allocate( biotic(n)%htotal(isc:iec,jsc:jec,1) )
  allocate( biotic(n)%csurf(isc:iec,jsc:jec) )
  allocate( biotic(n)%csat(isc:iec,jsc:jec) )
  allocate( biotic(n)%csat_csurf(isc:iec,jsc:jec) )
  allocate( biotic(n)%pco2surf(isc:iec,jsc:jec) )
  allocate( biotic(n)%dpco2(isc:iec,jsc:jec) )
  allocate( biotic(n)%sio2(isc:iec,jsc:jec,1) )
  if (biotic(n)%do_dic_virtual_flux) then  !{
    allocate( biotic(n)%vstf_dic(isc:iec,jsc:jec) )
  endif  !}
  if (biotic(n)%do_alk_virtual_flux) then  !{
    allocate( biotic(n)%vstf_alk(isc:iec,jsc:jec) )
  endif  !}
  if (biotic(n)%do_po4_virtual_flux) then  !{
    allocate( biotic(n)%vstf_po4(isc:iec,jsc:jec) )
  endif  !}
  if (biotic(n)%do_dop_virtual_flux) then  !{
    allocate( biotic(n)%vstf_dop(isc:iec,jsc:jec) )
  endif  !}
  if (biotic(n)%do_o2_virtual_flux) then  !{
    allocate( biotic(n)%vstf_o2(isc:iec,jsc:jec) )
  endif  !}
  allocate( biotic(n)%zforg(nk) )
  allocate( biotic(n)%zforg_alt(nk) )
  allocate( biotic(n)%zfca(nk) )
  allocate( biotic(n)%jprod(isc:iec,jsc:jec,nk) )
  allocate( biotic(n)%jpo4(isc:iec,jsc:jec,nk) )
  allocate( biotic(n)%jdop(isc:iec,jsc:jec,nk) )
  allocate( biotic(n)%jo2(isc:iec,jsc:jec,nk) )
  allocate( biotic(n)%jca(isc:iec,jsc:jec,nk) )
  allocate( biotic(n)%fc(isc:iec,jsc:jec,nk) )
  allocate( biotic(n)%fca(isc:iec,jsc:jec,nk) )
  allocate( biotic(n)%flux_poc(isc:iec,jsc:jec) )
  allocate( biotic(n)%flux_caco3(isc:iec,jsc:jec) )
  allocate( biotic(n)%nut_depl%mask(isc:iec,jsc:jec,12) )
  allocate( biotic(n)%no_caco3%mask(isc:iec,jsc:jec,12) )
  allocate( biotic(n)%norm_remin%mask(isc:iec,jsc:jec,12) )
  allocate( biotic(n)%r_bio_tau_a%mask(isc:iec,jsc:jec,12) )

  if (biotic(n)%add_phosphate) then  !{
    allocate( biotic(n)%jprod_norm(isc:iec,jsc:jec,nk) )
    allocate( biotic(n)%jpo4_add(isc:iec,jsc:jec,nk) )
  endif  !}

enddo  !}

!if (.not. equilibrium) then  !{
!allocate( biotic_cumul(isc:iec,jsc:jec,3) )
!endif  !}

!
!       initialize some arrays
!

do n = 1, instances  !{

  biotic(n)%jprod(:,:,:) = 0.0
  biotic(n)%jpo4(:,:,:)  = 0.0
  biotic(n)%jdop(:,:,:)  = 0.0
  biotic(n)%jo2(:,:,:)   = 0.0
  biotic(n)%jca(:,:,:)   = 0.0
  biotic(n)%fc(:,:,:)    = 0.0
  biotic(n)%fca(:,:,:)   = 0.0
  if (.not. use_waterflux) then  !{
    if (biotic(n)%do_dic_virtual_flux) then  !{
      biotic(n)%vstf_dic(:,:) = 0.0
    endif  !}
    if (biotic(n)%do_alk_virtual_flux) then  !{
      biotic(n)%vstf_alk(:,:) = 0.0
    endif  !}
    if (biotic(n)%do_po4_virtual_flux) then  !{
      biotic(n)%vstf_po4(:,:) = 0.0
    endif  !}
    if (biotic(n)%do_dop_virtual_flux) then  !{
      biotic(n)%vstf_dop(:,:) = 0.0
    endif  !}
    if (biotic(n)%do_o2_virtual_flux) then  !{
      biotic(n)%vstf_o2(:,:) = 0.0
    endif  !}
  endif  !}

  if (biotic(n)%add_phosphate) then  !{
    biotic(n)%jprod_norm(:,:,:) = 0.0
    biotic(n)%jpo4_add(:,:,:)  = 0.0
  endif  !}

enddo  !} n



return
end subroutine  allocate_arrays  !}
! </SUBROUTINE> NAME="allocate_arrays"


!#######################################################################
! <SUBROUTINE NAME="bc_interp">
!
! <DESCRIPTION>
!     Interpolates atmospheric CO2 and C-14 to the timestep of the model
!
!     ARGUMENT LIST -
!
!     Note: Variable type is given in square brackets (below)
!     (r-real, i-integer, l-logical, c-character; s-scaler, a-array).
!
!     INPUT:
!
!     [rs] - year_model = decimal year of model (e.g., 1990.67), as
!                         computed from the timestep and the year of
!                         the initialization of the simulation. This
!                         information is necessary to interpolate
!                         atmospheric levels of CO2 (atmco2_t) and 
!                         C-14 (atmc14_t) from the historical records 
!                         chosen for OCMIP-2 (from Enting et al. (1994).
!
!     [cs] - futr_scen  = Either 'S450' or 'S650' (in upper CASE) 
!                         Choses from the IPCC future scenarios by the 
!                         same names.  These scenarios are standard for
!                         OCMIP-2 future runs.  The choice of "futr_scen"
!                         does NOT affect "atmc14_t". Likewise, it does
!                         NOT affect "atmco2_t" prior to 1990.5.
!
!     OUTPUT: 
!
!     [rs] - atmco2_t   =  Atmospheric CO2 (ppm) at "year_model"
!
!     [ra] - atmc14_t   =  3-member array of atmospheric C-14 at
!                         year_model".  Sequentially, the 3 values
!                         correspond to forcing in 3 latitudinal bands:
!                         (1) 90S - 20S,
!                         (2) 20S - 20N, and
!                         (3) 20N - 90N.
!
!     James Orr, LSCE/CEA-CNRS, Saclay, France, 20 April 1999
! </DESCRIPTION>
!

subroutine bc_interp(year_model, futr_scen, atmco2_t, atmc14_t)  !{

INTEGER maxrec
PARAMETER (maxrec=1200)

INTEGER maxrec14, nzon14
PARAMETER (maxrec14=300, nzon14=3)

real, intent(in)        :: year_model
CHARACTER(len=4), intent(in)    :: futr_scen
real, intent(out)       :: atmco2_t, atmc14_t(nzon14)

integer, save           :: nco2rec, nc14rec(nzon14)
real, save      :: yrco2rec(maxrec), atmco2rec(maxrec)
real, save      :: yrc14rec(maxrec14,nzon14)
real, save      :: atmc14rec(maxrec14,nzon14)

real    :: year

INTEGER iz
INTEGER n
real    :: x


integer, save   :: ientry = 0

real, parameter :: yrco2min = 1765.0
real, parameter :: yrco2max = 2300.5
real, parameter :: yrc14min = 1765.0
real, parameter :: yrc14max = 1995.5

!     Counter "ientry" for number of entries in this SUBROUTINE
!     (Read atmospheric records on 1st entry ONLY)
!     ------------------------------------------------------------------

ientry = ientry + 1

IF (ientry .EQ. 1) THEN
    CALL read_co2atm(futr_scen(1:4), nco2rec, yrco2rec, atmco2rec)
    CALL read_c14atm (nc14rec, yrc14rec, atmc14rec)
ENDIF

!     ------------------------------------------------------------------
!     Interpolate atmospheric CO2 and C-14 to timestep (decimal years)
!     of the model
!     ------------------------------------------------------------------

!     CO2:
!     ----

IF (year_model .lt. yrco2min) THEN 
    year = yrco2min
ELSE IF(year_model .ge. yrco2min  .AND.  year_model .le. yrco2max) THEN 
    year = year_model
ELSE IF (year_model .gt. yrco2max) THEN 
    year = yrco2max
END IF 

!     Find relative POSITION n for year_model in CO2 record

call locate(yrco2rec, nco2rec, year, n)

!     Determine linear interpolation factor "x"

x = (year - yrco2rec(n)) / (yrco2rec(n+1) - yrco2rec(n))

!     Perform temporal interpolation for atmospheric CO2

atmco2_t = atmco2rec(n) * (1. - x) + atmco2rec(n+1) * x

!     C-14:
!     -----

IF (year_model .lt. yrc14min) THEN 
    year = yrc14min
ELSE IF(year_model .ge. yrc14min .AND. year_model .le. yrc14max) THEN 
    year = year_model
ELSE IF (year_model .gt. yrc14max) THEN 
    year = yrc14max
END IF         

!     Find relative POSITION n for year_model in C-14 record

call locate(yrc14rec(1,1), nc14rec(1), year, n)
    
!     Determine linear interpolation factor "x"

x = (year - yrc14rec(n,1)) / (yrc14rec(n+1,1) - yrc14rec(n,1))
    
!     Perform temporal linear interpolation for atmospheric C-14 
!     (on all 3 latitudinal zones for C-14)

DO iz=1,3
  atmc14_t(iz) = atmc14rec(n,iz) * (1. - x) + atmc14rec(n+1,iz) * x
END DO

RETURN
end subroutine  bc_interp  !}
! </SUBROUTINE> NAME="bc_interp"


!#######################################################################
! <SUBROUTINE NAME="locate">
!
! <DESCRIPTION>
!     After Numerical recipes:
!
!     Given an array XX of length N, and a given value of X, returns a
!     value of J such that X is between XX(J) and XX(J+1).  XX must be
!     monotonic, either increasing or decreasing. J=0 or J=N is
!     returned to indicate that X is out of range.      

!       New features:
!
!       If "period" is specified, then the array, xx, is considered
!       to be periodic with a period of "period". If "x_in" is out
!       of range, then add or subtract "period" once to attempt to 
!       make "x_in" be in range.
!
!       If "nearest" is specified, and true, then return "j" such
!       that it is the element of "xx" which is nearest to the value
!       of "x_in" (where "x_in" may have been modified by the value
!       "period", above). With this option, "j" will be in
!       the range 1 <= j <= n.
! </DESCRIPTION>
!

subroutine locate(xx , n, x_in, j, period, nearest)  !{

integer, intent(in)                             :: n
real, intent(in)                        :: x_in
real, dimension(n), intent(in)  :: xx
integer, intent(out)                            :: j
real, optional, intent(in)              :: period
logical, optional, intent(in)                   :: nearest

!
!       local variables
!

integer :: jl, ju, jm
real    :: x, xt
logical :: increasing

increasing = xx(1) .lt. xx(n)

if (present(period)) then  !{
  if (increasing) then  !{

! increasing array

    if (x_in .lt. xx(1)) then  !{

! original value less than start, therefore add period

      xt = x_in + period
      if (xt .gt. xx(n)) then  !{

! new value greater than end

        if (abs(x_in - xx(1)) .gt. abs(xt - xx(n))) then  !{

! new value closer to end than original value to start
! use new value

          x = xt
        else  !}{

! original value closer to start than new value to end
! use original value

          x = x_in
        endif  !}
      else  !}{

! new value in range
! use new value

        x = xt
      endif  !}
    elseif (x_in .gt. xx(n)) then  !}{

! original value greater than end, therefore subtract period

      xt = x_in - period
      if (xt .lt. xx(1)) then  !{

! new value less than start

        if (abs(xt - xx(1)) .lt. abs(x_in - xx(n))) then  !{

! new value closer to start than original value to end
! use new value

          x = xt
        else  !}{

! original value closer to end than new value to start
! use original value

          x = x_in
        endif  !}
      else  !}{

! new value in range
! use new value

        x = xt
      endif  !}
    else  !}{

! original value in range
! use original value

      x = x_in
    endif  !}
  else  !}{

! decreasing array

    if (x_in .gt. xx(1)) then  !{

! original value greater than start, therefore subtract period

      xt = x_in - period
      if (xt .lt. xx(n)) then  !{

! new value less than end

        if (abs(x_in - xx(1)) .gt. abs(xt - xx(n))) then  !{

! new value closer to end than original value to start
! use new value

          x = xt
        else  !}{

! original value closer to start than new value to end
! use original value

          x = x_in
        endif  !}
      else  !}{

! new value in range
! use new value

        x = xt
      endif  !}
    elseif (x_in .lt. xx(n)) then  !}{

! original value less than end, therefore add period

      xt = x_in + period
      if (xt .gt. xx(1)) then  !{

! new value greater than start

        if (abs(xt - xx(1)) .lt. abs(x_in - xx(n))) then  !{

! new value closer to start than original value to end
! use new value

          x = xt
        else  !}{

! original value closer to end than new value to start
! use original value

          x = x_in
        endif  !}
      else  !}{

! new value in range
! use new value

        x = xt
      endif  !}
    else  !}{

! original value in range
! use original value

      x = x_in
    endif  !}
  endif  !}
else  !}{

! no period specified
! use original value

  x = x_in
endif  !}

jl = 0
ju = n+1
10 continue
if (ju - jl .gt. 1) then
  jm = (ju + jl) / 2
  if (increasing .eqv. (x .gt. xx(jm))) then
    jl = jm
  else
    ju = jm
  endif
  go to 10
endif
j = jl

if (present(nearest)) then  !{
  if (nearest) then  !{
    if (j .eq. 0) then  !{
      j = 1
    elseif (j .lt. n) then  !}{
      if (abs(x - xx(j)) .gt. abs(x - xx(j+1))) then  !{
        j = j + 1
      endif  !}
    endif  !}
  endif  !}
endif  !}

return
end subroutine  locate  !}
! </SUBROUTINE> NAME="locate"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_biotic_bbc">
!
! <DESCRIPTION>
!     calculate the surface boundary conditions
! </DESCRIPTION>
!

subroutine ocmip2_biotic_bbc  !{

!
!-----------------------------------------------------------------------
!     modules (have to come first)
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!     arguments
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'ocmip2_biotic_bbc'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header =                                &
     '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

!
!-----------------------------------------------------------------------
!     local variables
!-----------------------------------------------------------------------
!

integer  :: i, j, n, kz

!
! =====================================================================
!     begin executable code
! =====================================================================
!

!
!   set the bottom flux of the column for phosphate to reflect a
!   regenerative flux from the sediments where the compensation
!   depth is greater than the bottom depth
!

do n = 1, instances  !{
  do j = jsc, jec  !{
    do i = isc, iec  !{
      kz = grid%kmt(i,j)
      if (kz .lt. biotic(n)%km_c .and. kz .gt. 0) then  !{
        t_prog(biotic(n)%ind_po4)%btf(i,j) =                    &
             -biotic(n)%flux_poc(i,j)
        t_prog(biotic(n)%ind_o2)%btf(i,j)  =                    &
             biotic(n)%o_2_p * biotic(n)%flux_poc(i,j)
        t_prog(biotic(n)%ind_dic)%btf(i,j) =                    &
             -(biotic(n)%c_2_p * (1.0 + biotic(n)%caco3_2_c) *  &
             biotic(n)%flux_poc(i,j))
        t_prog(biotic(n)%ind_alk)%btf(i,j) =                    &
             (biotic(n)%n_2_p - 2.0 * biotic(n)%c_2_p *         &
             biotic(n)%caco3_2_c) * biotic(n)%flux_poc(i,j)
      endif  !}
    enddo  !} i
  enddo  !} j
enddo  !} n

return

end subroutine  ocmip2_biotic_bbc  !}
! </SUBROUTINE> NAME="ocmip2_biotic_bbc"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_biotic_end">
!
! <DESCRIPTION>
!     Clean up various BIOTIC quantities for this run.
! </DESCRIPTION>
!

subroutine ocmip2_biotic_end(Thickness)  !{ 

type(ocean_thickness_type), intent(in) :: Thickness

!
!-----------------------------------------------------------------------
!     modules (have to come first)
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'ocmip2_biotic_end'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header =                                &
     '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

!
!-----------------------------------------------------------------------
!     local variables
!-----------------------------------------------------------------------
!

integer :: i
integer :: j
integer :: k
integer :: lun
integer :: n
real    :: total_alkalinity
real    :: total_dic
real    :: total_dop
real    :: total_o2
real    :: total_phosphate
real    :: total_salinity
real    :: total_temp

!
!-----------------------------------------------------------------------
!     statement functions
!-----------------------------------------------------------------------
!
!
! =====================================================================
!     begin executable code
! =====================================================================
!

!
!       integrate the total concentrations of some tracers
!       for the end of the run
!

total_phosphate = 0.0
total_dop = 0.0
total_o2 = 0.0
total_dic = 0.0
total_alkalinity = 0.0
total_temp = 0.0
total_salinity = 0.0

!
!       Use tau time index for the start of a run, and taup1 time
!       index for the end of a run so that we are integrating the
!       same time level and should therefore get identical results
!

do k = 1, nk  !{
  do j = jsc, jec  !{
    do i = isc, iec  !{
      total_temp = total_temp +                                 &
           t_prog(indtemp)%field(i,j,k,taup1) *                 &
           grid%dat(i,j) * grid%tmask(i,j,k) * Thickness%dht(i,j,k,tau)
      total_salinity = total_salinity +                         &
           t_prog(indsal)%field(i,j,k,taup1) *                  &
           grid%dat(i,j) * grid%tmask(i,j,k) * Thickness%dht(i,j,k,tau)
    enddo  !} i
  enddo  !} j
enddo  !} k

call mpp_sum(total_temp)
call mpp_sum(total_salinity)

write (stdout(),*) trim(note_header),                           &
     'Global integrals at end of run'
write (stdout(),'(/'' Total temperature  = '',es19.12,          &
                  &'' Gdeg-C m^3'')')                           &
            total_temp * 1.0e-15
write (stdout(),'(/'' Total salinity  = '',es19.12,             &
                  &'' Gsal m^3'')')                             &
            total_salinity * 1.0e-15

do n = 1, instances  !{
  do k = 1,nk  !{
    do j = jsc, jec  !{
      do i = isc, iec  !{
        total_phosphate = total_phosphate +                     &
             t_prog(biotic(n)%ind_po4)%field(i,j,k,taup1) *     &
             grid%dat(i,j) * grid%tmask(i,j,k) * Thickness%dht(i,j,k,tau)
        total_dop = total_dop +                                 &
             t_prog(biotic(n)%ind_dop)%field(i,j,k,taup1) *     &
             grid%dat(i,j) * grid%tmask(i,j,k) * Thickness%dht(i,j,k,tau)
        total_o2 = total_o2 +                                   &
             t_prog(biotic(n)%ind_o2)%field(i,j,k,taup1) *      &
             grid%dat(i,j) * grid%tmask(i,j,k) * Thickness%dht(i,j,k,tau)
        total_dic = total_dic +                                 &
             t_prog(biotic(n)%ind_dic)%field(i,j,k,taup1) *     &
             grid%dat(i,j) * grid%tmask(i,j,k) * Thickness%dht(i,j,k,tau)
        total_alkalinity = total_alkalinity +                   &
             t_prog(biotic(n)%ind_alk)%field(i,j,k,taup1) *     &
             grid%dat(i,j) * grid%tmask(i,j,k) * Thickness%dht(i,j,k,tau)
      enddo  !} i
    enddo  !} j
  enddo  !} k

  call mpp_sum(total_phosphate)
  call mpp_sum(total_dop)
  call mpp_sum(total_o2)
  call mpp_sum(total_dic)
  call mpp_sum(total_alkalinity)

  write (stdout(),*) '  Instance ', trim(biotic(n)%name)
  write (stdout(),                                              &
       '(/'' Total phosphate  = '',es19.12,'' Gmol-P m^3'')')   &
       total_phosphate * 1.0e-09
  write (stdout(),                                              &
       '(/'' Total DOP  = '',es19.12,'' Gmol-P m^3'')')         &
       total_DOP * 1.0e-09
  write (stdout(),                                              &
       '(/'' Total O2  = '',es19.12,'' Gmol-O m^3'')')          &
       total_o2 * 1.0e-09
  write (stdout(),                                              &
       '(/'' Total DIC  = '',es19.12,'' Gmol-C m^3'')')         &
       total_dic * 1.0e-09
  write (stdout(),                                              &
       '(/'' Total alkalinity  = '',es19.12,'' Geq m^3'')')     &
       total_alkalinity * 1.0e-09
  write (stdout(),                                              &
       '(/'' Total phosphorous  = '',es19.12,'' Gmol-P m^3'')') &
       (total_phosphate + total_dop) * 1.0e-09
  write (stdout(),                                              &
       '(/'' Total real O2  = '',es19.12,'' Gmol-O m^3'')')     &
       (total_o2 + biotic(n)%o_2_p * total_phosphate) * 1.0e-09
  write (stdout(),                                              &
       '(/'' Total Carbon  = '',es19.12,'' Gmol-C m^3'')')      &
       (total_dic + biotic(n)%c_2_p * total_dop) * 1.0e-09
  write (stdout(),                                              &
       '(/'' Total real alkalinity  = '',es19.12,'' Geq m^3'')')&
       (total_alkalinity + biotic(n)%n_2_p * total_phosphate) * 1.0e-09
  if (biotic(n)%global_atmosphere) then  !{
    write (stdout(),                                            &
         '(/'' Total Carbon (atm) = '',es19.12,'' Pg-C m^3'')') &
         biotic(n)%total_atm_co2(taup1) * 1.0e-15
    write (stdout(),                                            &
         '(/'' Total Carbon (ocn) = '',es19.12,'' Pg-C m^3'')') &
         (total_dic + biotic(n)%c_2_p * total_dop) *            &
         co2_mol_wgt * 1.0e-15
    write (stdout(),                                            &
         '(/'' Total Carbon (atm+ocn) = '',es19.12,'' Pg-C m^3'')') &
        ((total_dic + biotic(n)%c_2_p * total_dop) *            &
         co2_mol_wgt + biotic(n)%total_atm_co2(taup1)) * 1.0e-15
  endif  !}
enddo  !} n

if (.not. equilibrium) then  !{

!
!-----------------------------------------------------------------------
!     save out the cumulative input and the perturbation
!      time
!-----------------------------------------------------------------------
!

  write (stdout(),*)
  write (stdout(),*) trim(note_header),                         &
       'Saving additional output to a restart file...'

  do n = 1, instances  !{

    call write_data(biotic(n)%file_out,                         &
         'pert_start_year', biotic(n)%pert_start_year)
    call write_data(biotic(n)%file_out,                         &
         'pert_start_year_model',                               &
         biotic(n)%pert_start_year_model)
    call write_data(biotic(n)%file_out,                         &
         'pert_time', biotic(n)%pert_time)
!write (lun) pert_starttime, pert_starttime_mod, pert_time,     &
!biotic_cumul
!write (lun) gmm_pco2atm,                                       &
!mgf_biotic_dic,                                        &
!mgvf_biotic_dic, mgvf_biotic_alk,              &
!mgv_biotic_dic, mgv_biotic_alk,                        &
!mgs_biotic_dic, mgs_biotic_alk,                        &
!mg_biotic_pco2surf, mg_biotic_dpco2,           &
!gcf_biotic_dic,                                        &
!gcvf_biotic_dic, gcvf_biotic_alk

  enddo  !} n

endif  !}

!
!-----------------------------------------------------------------------
!       save out additional information for a restart
!-----------------------------------------------------------------------
!

write(stdout(),*)

do n = 1, instances  !{

  write(stdout(),*) trim(note_header),                          &
       'Writing additional restart information for instance ',  &
       trim(biotic(n)%name)

  if (biotic(n)%global_atmosphere) then  !{
    call write_data(biotic(n)%file_out,                         &
         'total_atm_co2', biotic(n)%total_atm_co2(tau))
    call write_data(biotic(n)%file_out,                         &
         'total_atm_co2', biotic(n)%total_atm_co2(taup1))
    write (stdout(),'(/1x,a,es16.9,a,a)')                       &
          'Total atmospheric mass of CO2 (tau  ) = ',           &
          biotic(n)%total_atm_co2(tau  ) * 1.0e-15,             &
          ' (Pg) for instance ', trim(biotic(n)%name)
    write (stdout(),'(/1x,a,es16.9,a,a)')                       &
          'Total atmospheric mass of CO2 (taup1) = ',           &
          biotic(n)%total_atm_co2(taup1) * 1.0e-15,             &
          ' (Pg) for instance ', trim(biotic(n)%name)
  endif  !}

  call write_data(biotic(n)%file_out,                           &
       'htotal', biotic(n)%htotal,                              &
       domain=Domain%domain2d)

  if (.not. use_waterflux) then  !{

    if (biotic(n)%do_dic_virtual_flux .or.                      &
        biotic(n)%do_alk_virtual_flux .or.                      &
        biotic(n)%do_po4_virtual_flux .or.                      &
        biotic(n)%do_dop_virtual_flux .or.                      &
        biotic(n)%do_o2_virtual_flux) then  !{
      call write_data(biotic(n)%file_out,                       &
           'sal_global', biotic(n)%sal_global)
      call write_data(biotic(n)%file_out,                       &
           'sal_global_wrk', biotic(n)%sal_global_wrk)
      call write_data(biotic(n)%file_out,                       &
           'global_wrk_duration',                               &
           biotic(n)%global_wrk_duration)
    endif  !}
    if (biotic(n)%do_dic_virtual_flux) then  !{
      call write_data(biotic(n)%file_out,                       &
           'dic_global', biotic(n)%dic_global)
      call write_data(biotic(n)%file_out,                       &
           'dic_global_wrk', biotic(n)%dic_global_wrk)
    endif  !}
    if (biotic(n)%do_alk_virtual_flux) then  !{
      call write_data(biotic(n)%file_out,                       &
           'alk_global_wrk', biotic(n)%alk_global_wrk)
      call write_data(biotic(n)%file_out,                       &
           'alk_global', biotic(n)%alk_global)
    endif  !}
    if (biotic(n)%do_po4_virtual_flux) then  !{
      call write_data(biotic(n)%file_out,                       &
           'po4_global_wrk', biotic(n)%po4_global_wrk)
      call write_data(biotic(n)%file_out,                       &
           'po4_global', biotic(n)%po4_global)
    endif  !}
    if (biotic(n)%do_dop_virtual_flux) then  !{
      call write_data(biotic(n)%file_out,                       &
           'dop_global_wrk', biotic(n)%dop_global_wrk)
      call write_data(biotic(n)%file_out,                       &
           'dop_global', biotic(n)%dop_global)
    endif  !}
    if (biotic(n)%do_o2_virtual_flux) then  !{
      call write_data(biotic(n)%file_out,                       &
           'o2_global_wrk', biotic(n)%o2_global_wrk)
      call write_data(biotic(n)%file_out,                       &
           'o2_global', biotic(n)%o2_global)
    endif  !}

    if (biotic(n)%do_dic_virtual_flux .or.                      &
        biotic(n)%do_alk_virtual_flux .or.                      &
        biotic(n)%do_po4_virtual_flux .or.                      &
        biotic(n)%do_dop_virtual_flux .or.                      &
        biotic(n)%do_o2_virtual_flux) then  !{
      write (stdout(),'(/1x,a,es16.9,a,a)')                     &
            'Annual, global, surface mean salinity = ',         &
            biotic(n)%sal_global, ' (PSU) for instance ',       &
            trim(biotic(n)%name)
    endif  !}
    if (biotic(n)%do_dic_virtual_flux) then  !{
      write (stdout(),'(1x,a,es16.9,a,a)')                      &
            'Annual, global, surface mean DIC      = ',         &
            biotic(n)%dic_global, ' (mol/m^3) for instance ',   &
            trim(biotic(n)%name)
    endif  !}
    if (biotic(n)%do_alk_virtual_flux) then  !{
      write (stdout(),'(1x,a,es16.9,a,a)')                      &
            'Annual, global, surface mean ALK      = ',         &
            biotic(n)%alk_global, ' (eq/m^3) for instance ',    &
            trim(biotic(n)%name)
    endif  !}
    if (biotic(n)%do_po4_virtual_flux) then  !{
      write (stdout(),'(1x,a,es16.9,a,a)')                      &
            'Annual, global, surface mean PO4      = ',         &
            biotic(n)%po4_global, ' (mol/m^3) for instance ',   &
            trim(biotic(n)%name)
    endif  !}
    if (biotic(n)%do_dop_virtual_flux) then  !{
      write (stdout(),'(1x,a,es16.9,a,a)')                      &
            'Annual, global, surface mean DOP      = ',         &
            biotic(n)%dop_global, ' (mol/m^3) for instance ',   &
            trim(biotic(n)%name)
    endif  !}
    if (biotic(n)%do_o2_virtual_flux) then  !{
      write (stdout(),'(1x,a,es16.9,a,a)')                      &
            'Annual, global, surface mean O2       = ',         &
            biotic(n)%o2_global, ' (mol/m^3) for instance ',    &
            trim(biotic(n)%name)
    endif  !}

  endif  !}

  write (stdout(),*) trim(note_header),                         &
       'Done writing additional restart information for instance ',&
       trim(biotic(n)%name)

enddo  !} n

return
end subroutine  ocmip2_biotic_end  !}
! </SUBROUTINE> NAME="ocmip2_biotic_end"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_biotic_sbc">
!
! <DESCRIPTION>
!     Calculate the surface boundary conditions
! </DESCRIPTION>
!

subroutine ocmip2_biotic_sbc(robert)  !{

!
!-----------------------------------------------------------------------
!     modules (have to come first)
!-----------------------------------------------------------------------
!

use ocmip2_co2calc_mod
use mpp_mod, only : mpp_sum
use time_interp_external_mod, only: time_interp_external
use time_manager_mod, only: days_in_year, days_in_month,        &
     get_date, set_date

!
!-----------------------------------------------------------------------
!     arguments
!-----------------------------------------------------------------------
!

real, intent(in)        :: robert       ! robert time-filter coefficient

!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'ocmip2_biotic_sbc'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header =                                &
     '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

!
!       coefficients for O2 saturation
!

real, parameter :: a_0 = 2.00907
real, parameter :: a_1 = 3.22014
real, parameter :: a_2 = 4.05010
real, parameter :: a_3 = 4.94457
real, parameter :: a_4 = -2.56847e-01
real, parameter :: a_5 = 3.88767
real, parameter :: b_0 = -6.24523e-03
real, parameter :: b_1 = -7.37614e-03
real, parameter :: b_2 = -1.03410e-02
real, parameter :: b_3 = -8.17083e-03
real, parameter :: c_0 = -4.88682e-07

!
!-----------------------------------------------------------------------
!     local variables
!-----------------------------------------------------------------------
!

integer :: i
integer :: j
integer :: k
integer :: n
integer :: m
integer :: kz
real    :: pco2atm_global
integer :: hour
integer :: day
real    :: days_in_this_year
integer :: minute
integer :: month
integer :: num_days
integer :: second
integer :: year

!
! =====================================================================
!     begin executable code
! =====================================================================
!

!if (no_sbc) then  !{
  !write (stdout(),*) trim(warn_header), 'Skipping calculations'
  !return
!else  !}{
  !write (stdout(),*) trim(note_header), 'Doing calculations'
!endif  !}

if (use_waterflux) then  !{

  do n = 1, instances  !{

    if (biotic(n)%do_dic_virtual_flux) then  !{
      do j = jsc, jec  !{
        do i = isc, iec  !{
          T_prog(biotic(n)%ind_dic)%tpme(i,j) =                 &
               T_prog(biotic(n)%ind_dic)%field(i,j,1,tau)
          T_prog(biotic(n)%ind_dic)%triver(i,j) =               &
               T_prog(biotic(n)%ind_dic)%field(i,j,1,tau)
        enddo  !} i
      enddo  !} j
    endif  !}
    if (biotic(n)%do_alk_virtual_flux) then  !{
      do j = jsc, jec  !{
        do i = isc, iec  !{
          T_prog(biotic(n)%ind_alk)%tpme(i,j) =                 &
               T_prog(biotic(n)%ind_alk)%field(i,j,1,tau)
          T_prog(biotic(n)%ind_alk)%triver(i,j) =               &
               T_prog(biotic(n)%ind_alk)%field(i,j,1,tau)
        enddo  !} i
      enddo  !} j
    endif  !}
    if (biotic(n)%do_po4_virtual_flux) then  !{
      do j = jsc, jec  !{
        do i = isc, iec  !{
          T_prog(biotic(n)%ind_po4)%tpme(i,j) =                 &
               T_prog(biotic(n)%ind_po4)%field(i,j,1,tau)
          T_prog(biotic(n)%ind_po4)%triver(i,j) =               &
               T_prog(biotic(n)%ind_po4)%field(i,j,1,tau)
        enddo  !} i
      enddo  !} j
    endif  !}
    if (biotic(n)%do_dop_virtual_flux) then  !{
      do j = jsc, jec  !{
        do i = isc, iec  !{
          T_prog(biotic(n)%ind_dop)%tpme(i,j) =                 &
               T_prog(biotic(n)%ind_dop)%field(i,j,1,tau)
          T_prog(biotic(n)%ind_dop)%triver(i,j) =               &
               T_prog(biotic(n)%ind_dop)%field(i,j,1,tau)
        enddo  !} i
      enddo  !} j
    endif  !}
    if (biotic(n)%do_o2_virtual_flux) then  !{
      do j = jsc, jec  !{
        do i = isc, iec  !{
          T_prog(biotic(n)%ind_o2)%tpme(i,j) =                 &
               T_prog(biotic(n)%ind_o2)%field(i,j,1,tau)
          T_prog(biotic(n)%ind_o2)%triver(i,j) =               &
               T_prog(biotic(n)%ind_o2)%field(i,j,1,tau)
        enddo  !} i
      enddo  !} j
    endif  !}

  enddo  !} n

endif  !}

if (.not. equilibrium .and. .not. control) then  !{

!
!-----------------------------------------------------------------------
!     Update pert_time
!-----------------------------------------------------------------------
!

  call get_date(time%model_time,                                &
                year, month, day, hour, minute, second)
  num_days = days_in_year(time%model_time)
  days_in_this_year = 0.0
  do m = 1, month - 1
    days_in_this_year = days_in_this_year +                     &
                       days_in_month(set_date(year, m, 1))
  enddo
  days_in_this_year = days_in_this_year + day - 1 + hour/24.0 + &
                      minute/1440.0 + second/86400.0

!
!---------------------------------------------------------------------
!     calculate atmospheric pressures by calling routine bc_interp
!---------------------------------------------------------------------
!

  do n = 1, instances  !{

!
!-----------------------------------------------------------------------
!     assign perturbation time (real years e.g. 1935.67)
!-----------------------------------------------------------------------
!	

    biotic(n)%pert_time = year + days_in_this_year/num_days -   &
         biotic(n)%pert_start_year_model + biotic(n)%pert_start_year

    call bc_interp(biotic(n)%pert_time, scenario(1:4),          &
                   atmco2_t, atmc14_t)

!
!-----------------------------------------------------------------------
!     "interpolate" in space
!-----------------------------------------------------------------------
!

    biotic(n)%pco2atm(:,:) = atmco2_t

!
!-----------------------------------------------------------------------
!     info
!-----------------------------------------------------------------------
!

    if (end_of_month) then  !{
      write (stdout(), '(a,a,a)') trim(note_header),     &
           'Atmospheric concentrations : '
      write (stdout(),'(a,f9.4)') '    BIOTIC-DIC : ', atmco2_t
    endif  !}

  enddo  !} n

endif  !}

!
!---------------------------------------------------------------------
!  Compute the Schmidt number of CO2 in seawater using the 
!  formulation presented by Wanninkhof (1992, J. Geophys. Res., 97,
!  7373-7382).
!---------------------------------------------------------------------
!

do j = jsc, jec  !{
  do i = isc, iec  !{
    sc_co2(i,j) = 2073.1 + t_prog(indtemp)%field(i,j,1,taum1) * &
         (-125.62 + t_prog(indtemp)%field(i,j,1,taum1) *        &
          (3.6276 + t_prog(indtemp)%field(i,j,1,taum1) *        &
           (-0.043219))) * grid%tmask(i,j,1)
  enddo  !} i
enddo  !} j 

!
!---------------------------------------------------------------------
!  Compute the Schmidt number of O2 in seawater using the 
!  formulation proposed by Keeling et al. (1998, Global Biogeochem.
!  Cycles, 12, 141-163).
!---------------------------------------------------------------------
!

do j = jsc, jec  !{
  do i = isc, iec  !{
    sc_o2(i,j) = 1638.0 + t_prog(indtemp)%field(i,j,1,taum1) *  &
         (-81.83 + t_prog(indtemp)%field(i,j,1,taum1) *         &
          (1.483 + t_prog(indtemp)%field(i,j,1,taum1) *         &
           (-0.008004))) * grid%tmask(i,j,1)
  enddo  !} i
enddo  !} j 

!
!---------------------------------------------------------------------
!     calculate interpolated xkw, seaice fraction and atmospheric
!       pressure
!---------------------------------------------------------------------
!

call time_interp_external(pistonveloc_id, time%model_time, xkw_t)
call time_interp_external(seaicefract_id, time%model_time, fice_t)
call time_interp_external(atmpress_id, time%model_time, patm_t)

!
!---------------------------------------------------------------------
!     calculate csurf, csat and csat - csurf via the routine co2calc
!        input and output units are in mol/m^3
!---------------------------------------------------------------------
!
!
!       first calculate some derived quantities
!

do n = 1, instances  !{

  if (biotic(n)%global_atmosphere) then  !{

!
!       determine the atmospheric pCO2 concetration for this time-step
!

    pco2atm_global = (biotic(n)%total_atm_co2(taum1) /          &
                      total_atm_mass) *                         &
                     (atm_mol_wgt / co2_mol_wgt) * 1.0e+06
    biotic(n)%pco2atm(:,:) = pco2atm_global

  endif  !}

enddo  !} n

call init_ocmip2_co2calc(                                       &
     time%model_time, isc, iec, jsc, jec, 1,                    &
     t_prog(indtemp)%field(isc:iec,jsc:jec,1,taum1),            &
     t_prog(indsal)%field(isc:iec,jsc:jec,1,taum1))

do n = 1, instances  !{
  do k = 1, 1  !{
    do j = jsc, jec  !{
      do i = isc, iec  !{
        htotallo(i,j,k) = htotal_scale_lo(i,j) * biotic(n)%htotal(i,j,k)
        htotalhi(i,j,k) = htotal_scale_hi(i,j) * biotic(n)%htotal(i,j,k)
      enddo  !} i
    enddo  !} j 
  enddo  !} k 

  call ocmip2_co2calc(isc, iec, jsc, jec, 1,                    &
       grid%tmask(isc:iec,jsc:jec,1),                           &
       t_prog(biotic(n)%ind_dic)%field(isc:iec,jsc:jec,1,taum1),&
       t_prog(biotic(n)%ind_alk)%field(isc:iec,jsc:jec,1,taum1),&
       t_prog(biotic(n)%ind_po4)%field(isc:iec,jsc:jec,1,taum1),&
       biotic(n)%sio2,                                          &
       htotallo, htotalhi, biotic(n)%htotal,                    &
       biotic(n)%pco2atm, patm_t(isc:iec,jsc:jec),              &
       biotic(n)%csurf, biotic(n)%csat,                         &
       biotic(n)%csat_csurf,                                    &
       biotic(n)%pco2surf, biotic(n)%dpco2)
enddo  !} n 

!
!---------------------------------------------------------------------
!  Compute the oxygen saturation concentration at 1 atm total
!  pressure in mol/m^3 given the temperature (t, in deg C) and
!  the salinity (s, in permil)
!
!  From Garcia and Gordon (1992), Limnology and Oceonography.
!  The formula used is from page 1310, eq (8).
!
!  *** Note: the "a3*ts^2" term (in the paper) is incorrect. ***
!  *** It shouldn't be there.                                ***
!
!  o2_saturation is defiend between T(freezing) <= T <= 40 deg C and
!                                   0 permil <= S <= 42 permil
!
! check value: T = 10 deg C, S = 35 permil,
!              o2_saturation = 0.282015 mol/m^3
!---------------------------------------------------------------------
!

do j = jsc, jec  !{
  do i = isc, iec  !{
    tt(i) = 298.15 - t_prog(indtemp)%field(i,j,1,taum1)
    tk(i) = 273.15 + t_prog(indtemp)%field(i,j,1,taum1)
    ts(i) = log(tt(i) / tk(i))
    ts2(i) = ts(i) * ts(i)
    ts3(i) = ts2(i) * ts(i)
    ts4(i) = ts3(i) * ts(i)
    ts5(i) = ts4(i) * ts(i)
    o2_saturation(i,j) =                                        &
         exp(a_0 + a_1*ts(i) + a_2*ts2(i) +                     &
             a_3*ts3(i) + a_4*ts4(i) + a_5*ts5(i) +             &
             t_prog(indsal)%field(i,j,1,taum1) *                &
             (b_0 + b_1*ts(i) + b_2*ts2(i) + b_3*ts3(i) +       &
              c_0*t_prog(indsal)%field(i,j,1,taum1)))
  enddo  !} i
enddo  !} j 

!
!       convert from ml/l to mol/m^3
!

do j = jsc, jec  !{
  do i = isc, iec  !{
    o2_saturation(i,j) = o2_saturation(i,j) * (1000.0/22391.6)
  enddo  !} i
enddo  !} j 

!
!---------------------------------------------------------------------
!     calculate piston-velocities
!      includomg  effect of sea-ice
!      xkw is given in cm/s (converted in read_biotic_bc), therefore
!      kw is also in cm/s
!---------------------------------------------------------------------
!

do j = jsc, jec  !{
  do i = isc, iec  !{
    kw_co2(i,j) = (1.0 - fice_t(i,j)) * xkw_t(i,j) *            &
                   sqrt(660.0/sc_co2(i,j)) * grid%tmask(i,j,1)
    kw_o2(i,j) = (1.0 - fice_t(i,j)) * xkw_t(i,j) *             &
                  sqrt(660.0/sc_o2(i,j)) * grid%tmask(i,j,1)
  enddo  !} i
enddo  !} j 

!
!---------------------------------------------------------------------
!     calculate surface fluxes for BIOTICs
!       kw is in cm/s, csat and t are in mol/cm3, therefore
!       stf is in mol/cm2/s
!---------------------------------------------------------------------
!

do n = 1, instances  !{
  do j = jsc, jec  !{
    do i = isc, iec  !{
      t_prog(biotic(n)%ind_dic)%stf(i,j) = kw_co2(i,j) *        &
            biotic(n)%csat_csurf(i,j)
      t_prog(biotic(n)%ind_o2)%stf(i,j) = kw_o2(i,j) *          &
            (o2_saturation(i,j) * patm_t(i,j) -                 &
             t_prog(biotic(n)%ind_o2)%field(i,j,1,taum1))
    enddo  !} i
  enddo  !} j 
enddo  !} n 

!
!---------------------------------------------------------------------
!     add in the virtual fluxes as defined by equations (2) and (3)
!     in the biotic HOWTO.
!       Note: the factor of 1000 is to convert the delta salinity from
!             model units to PSU
!---------------------------------------------------------------------
!

if (.not. use_waterflux) then  !{

  do n = 1, instances  !{

    if (biotic(n)%do_dic_virtual_flux) then  !{
      do j = jsc, jec  !{
        do i = isc, iec  !{
          biotic(n)%vstf_dic(i,j) =                             &
               t_prog(indsal)%stf(i,j) *                        &
               biotic(n)%dic_global / biotic(n)%sal_global
          t_prog(biotic(n)%ind_dic)%stf(i,j) =                  &
               t_prog(biotic(n)%ind_dic)%stf(i,j) +             &
               biotic(n)%vstf_dic(i,j)
        enddo  !} i
      enddo  !} j
    endif  !}
    if (biotic(n)%do_alk_virtual_flux) then  !{
      do j = jsc, jec  !{
        do i = isc, iec  !{
          biotic(n)%vstf_alk(i,j) =                             &
               t_prog(indsal)%stf(i,j) *                        &
               biotic(n)%alk_global / biotic(n)%sal_global
          t_prog(biotic(n)%ind_alk)%stf(i,j) =                  &
               biotic(n)%vstf_alk(i,j)
        enddo  !} i
      enddo  !} j
    endif  !}
    if (biotic(n)%do_po4_virtual_flux) then  !{
      do j = jsc, jec  !{
        do i = isc, iec  !{
          biotic(n)%vstf_po4(i,j) =                             &
               t_prog(indsal)%stf(i,j) *                        &
               biotic(n)%po4_global / biotic(n)%sal_global
          t_prog(biotic(n)%ind_po4)%stf(i,j) =                  &
               biotic(n)%vstf_po4(i,j)
        enddo  !} i
      enddo  !} j
    endif  !}
    if (biotic(n)%do_dop_virtual_flux) then  !{
      do j = jsc, jec  !{
        do i = isc, iec  !{
          biotic(n)%vstf_dop(i,j) =                             &
               t_prog(indsal)%stf(i,j) *                        &
               biotic(n)%dop_global / biotic(n)%sal_global
          t_prog(biotic(n)%ind_dop)%stf(i,j) =                  &
               biotic(n)%vstf_dop(i,j)
        enddo  !} i
      enddo  !} j
    endif  !}
    if (biotic(n)%do_o2_virtual_flux) then  !{
      do j = jsc, jec  !{
        do i = isc, iec  !{
          biotic(n)%vstf_o2(i,j) =                              &
               t_prog(indsal)%stf(i,j) *                        &
               biotic(n)%o2_global / biotic(n)%sal_global
          t_prog(biotic(n)%ind_o2)%stf(i,j) =                   &
               biotic(n)%vstf_o2(i,j)
        enddo  !} i
      enddo  !} j
    endif  !}

  enddo  !} n 

endif  !}

!if (.not. equilibrium) then  !{

!
!---------------------------------------------------------------------
!     calculate cumulative input [mol/cm2]
!---------------------------------------------------------------------
!

!do j = jsc, jec  !{
!do i = isc, iec  !{
!biotic_cumul(i,j,1) =                          &
!biotic_cumul(i,j,1) +                  &
!t_prog(ind_dic(n))%stf(i,j) * c2dtts * 0.5
!biotic_cumul(i,j,2) =                          &
!biotic_cumul(i,j,2) +                  &
!virtual_stf(i,j,1) * c2dtts * 0.5
!biotic_cumul(i,j,3) =                          &
!biotic_cumul(i,j,3) +                  &
!virtual_stf(i,j,2) * c2dtts * 0.5
!enddo  !} i 
!enddo  !} j 

!endif  !}


do n = 1, instances  !{

  if (biotic(n)%global_atmosphere) then  !{

!
!       integrate the surface flux of CO_2
!

!
!       initialize the accumulator to zero 
!

    biotic(n)%atm_co2_flux = 0.0

!
!       integrate over the local domain
!

    do j = jsc, jec  !{
      do i = isc, iec  !{
        biotic(n)%atm_co2_flux = biotic(n)%atm_co2_flux +       &
             t_prog(biotic(n)%ind_dic)%stf(i,j) *               &
             grid%dat(i,j) * grid%tmask(i,j,1)
      enddo  !} i
    enddo  !} j

!
!       sum over all domains
!

    call mpp_sum(biotic(n)%atm_co2_flux)

!
!       calculate the new atmospheric mass
!       first change the units from mol/s to g/s
!

    biotic(n)%atm_co2_flux = biotic(n)%atm_co2_flux * co2_mol_wgt
    biotic(n)%total_atm_co2(taup1) =                            &
         biotic(n)%total_atm_co2(taum1) -                       &
         biotic(n)%atm_co2_flux * dtts

!
!       do the Robert filter
!

    biotic(n)%total_atm_co2(tau) =                              &
         biotic(n)%total_atm_co2(tau) +                         &
         robert*(0.5*(biotic(n)%total_atm_co2(taup1) +          &
         biotic(n)%total_atm_co2(taum1)) -                      &
         biotic(n)%total_atm_co2(tau))

  endif  !}

enddo  !} n

return

end subroutine  ocmip2_biotic_sbc  !}
! </SUBROUTINE> NAME="ocmip2_biotic_sbc"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_biotic_init">
!
! <DESCRIPTION>
!       Set up any extra fields needed by the tracer packages
!
!       Save pointers to various "types", such as Grid and Domains.
! </DESCRIPTION>

subroutine ocmip2_biotic_init  !{

!
!       local parameters
!

character(len=64), parameter    :: sub_name = 'ocmip2_biotic_init'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header =                                &
     '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

real, parameter :: rho_avg = 1024.5
real, parameter :: sperd = 24.0 * 3600.0
real, parameter :: spery = 365.25 * sperd

!
!-----------------------------------------------------------------------
!       arguments
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!       local variables
!-----------------------------------------------------------------------
!

integer                                                 :: n
character(len=fm_field_name_len)                        :: name
character(len=fm_path_name_len)                         :: path_to_names
character(len=fm_field_name_len+1)                      :: suffix
character(len=fm_string_len)                            :: string
character(len=fm_field_name_len+3)                      :: long_suffix
logical, dimension(12)                                  :: t_mask
character(len=256)                                      :: caller_str
character(len=fm_string_len), pointer, dimension(:)     :: good_list

!
!       Initialize the ocmip2 biotic package
!

package_index = otpm_set_tracer_package(package_name,            &
     file_in = default_file_in, file_out = default_file_out,     &
     caller = trim(mod_name) // '(' // trim(sub_name) // ')')

!
!       Check whether to use this package
!

path_to_names = '/ocean_mod/tracer_packages/' // trim(package_name) // '/names'
instances = fm_get_length(path_to_names)
if (instances .lt. 0) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not get number of instances')
endif  !}

!
!       Check some things
!

write (stdout(),*)
if (instances .eq. 0) then  !{
  write (stdout(),*) trim(note_header), ' No instances'
  do_ocmip2_biotic = .false.
else  !}{
  if (instances .eq. 1) then  !{
    write (stdout(),*) trim(note_header), ' ', instances, ' instance'
  else  !}{
    write (stdout(),*) trim(note_header), ' ', instances, ' instances'
  endif  !}
  do_ocmip2_biotic = .true.
endif  !}

!
!       Return if we don't want to use this package,
!       after changing the list back
!

if (.not. do_ocmip2_biotic) then  !{
  return
endif  !}

!
!       allocate storage for biotic array
!

allocate ( biotic(instances) )

!
!       loop over the names, saving them into the biotic array
!

do n = 1, instances  !{

  if (fm_get_value(path_to_names, name, index = n)) then  !{
    biotic(n)%name = name
  else  !}{
    write (name,*) n
    call mpp_error(FATAL, trim(error_header) //        &
         'Bad field name for index ' // trim(name))
  endif  !}

enddo  !}

!
!       Set up the field input
!

do n = 1, instances  !{

  name = biotic(n)%name
  if (name(1:1) .eq. '_') then  !{
    suffix = ' '
    long_suffix = ' '
  else  !}{
    suffix = '_' // name
    long_suffix = ' (' // trim(name) // ')'
  endif  !}

!
!       PO4
!

  biotic(n)%ind_po4 = otpm_set_prog_tracer('po4' // suffix,     &
       package_name,                                            &
       longname = 'Phosphate' // trim(long_suffix),             &
       units = 'mol/m^3', flux_units = 'mol/m^2/s',             &
       caller = trim(mod_name)//'('//trim(sub_name)//')')

!
!       DOP
!

  biotic(n)%ind_dop = otpm_set_prog_tracer('dop' // suffix,     &
       package_name,                                            &
       longname = 'DOP' // trim(long_suffix),                   &
       units = 'mol/m^3', flux_units = 'mol/m^2/s',             &
       caller = trim(mod_name)//'('//trim(sub_name)//')')

!
!       DIC
!

  biotic(n)%ind_dic = otpm_set_prog_tracer('dic' // suffix,     &
       package_name,                                            &
       longname = 'DIC' // trim(long_suffix),                   &
       units = 'mol/m^3', flux_units = 'mol/m^2/s',             &
       caller = trim(mod_name)//'('//trim(sub_name)//')')

!
!       O2
!

  biotic(n)%ind_o2 = otpm_set_prog_tracer('o2' // suffix,       &
       package_name,                                            &
       longname = 'Oxygen' // trim(long_suffix),                &
       units = 'mol/m^3', flux_units = 'mol/m^2/s',             &
       caller = trim(mod_name)//'('//trim(sub_name)//')')

!
!       ALK
!

  biotic(n)%ind_alk = otpm_set_prog_tracer('alk' // suffix,     &
       package_name,                                            &
       longname = 'Alkalinity' // trim(long_suffix),            &
       units = 'mol/m^3', flux_units = 'mol/m^2/s',             &
       caller = trim(mod_name)//'('//trim(sub_name)//')')

enddo  !} n

!
!-----------------------------------------------------------------------
!       Process the namelists
!-----------------------------------------------------------------------
!

!
!       Add the package name to the list of good namelists, to be used
!       later for a consistency check
!

if (fm_new_value('/ocean_mod/GOOD/good_namelists', package_name, append = .true.) .le. 0) then  !{
  call mpp_error(FATAL, trim(error_header) //                           &
       ' Could not add ' // trim(package_name) // ' to "good_namelists" list')
endif  !}

!
!-----------------------------------------------------------------------
!       Set up the *global* namelist
!-----------------------------------------------------------------------
!

caller_str = trim(mod_name) // '(' // trim(sub_name) // ')'

call otpm_start_namelist(package_name, '*global*', caller = caller_str, no_overwrite = .true., &
     check = .true.)

call otpm_set_value('atmpress_file', 'INPUT/atmpress_ocmip2.nc')
call otpm_set_value('atmpress_name', 'atmpress')
call otpm_set_value('pistonveloc_file', 'INPUT/pistonveloc_ocmip2.nc')
call otpm_set_value('pistonveloc_name', 'pistonveloc')
call otpm_set_value('seaicefract_file', 'INPUT/f_ice_ocmip2.nc')
call otpm_set_value('seaicefract_name', 'f_ice')
call otpm_set_value('po4_star_file', 'INPUT/po4_star_ocmip2.nc')
call otpm_set_value('po4_star_name', 'po4_star')
call otpm_set_value('pert_first', .false.)
call otpm_set_value('pert_start_year', 1765)
call otpm_set_value('pert_start_month', 1)
call otpm_set_value('pert_start_day', 1)
call otpm_set_value('pert_start_hour', 0)
call otpm_set_value('pert_start_minute', 0)
call otpm_set_value('pert_start_second', 0)
call otpm_set_value('control', .false.)
call otpm_set_value('historical ', .false.)
call otpm_set_value('scenario', 'Equilibrium')
call otpm_set_value('htotal_scale_lo_in', 0.01 )
call otpm_set_value('htotal_scale_hi_in', 100.0)
call otpm_set_value('htotal_in', 1.0e-08)
call otpm_set_value('sal_global', 35.0)             ! PSU

call otpm_end_namelist(package_name, '*global*', caller = caller_str, check = .true.)

!
!-----------------------------------------------------------------------
!       Set up the instance namelists
!-----------------------------------------------------------------------
!

t_mask(:) = .true.

do n = 1, instances  !{

!
!       create the instance namelist
!

  call otpm_start_namelist(package_name, biotic(n)%name, caller = caller_str, no_overwrite = .true., &
       check = .true.)

  call otpm_set_value('pco2atm_const', 278.0)                      ! ppm
  call otpm_set_value('compensation_depth', 75.0)                  ! m
  call otpm_set_value('add_phosphate', .false.)
  call otpm_set_value('do_dic_virtual_flux', .false.)
  call otpm_set_value('do_alk_virtual_flux', .false.)
  call otpm_set_value('do_po4_virtual_flux', .false.)
  call otpm_set_value('do_dop_virtual_flux', .false.)
  call otpm_set_value('do_o2_virtual_flux', .false.)
  call otpm_set_value('martin_coeff', 0.9)
  call otpm_set_value('martin_coeff_alt', 0.0)
  call otpm_set_value('ca_remin_depth', 3500.0)
  call otpm_set_value('soft_tissue_pump', .false.)
  call otpm_set_value('stp_temperature', 10.0)
  call otpm_set_value('stp_salinity', 34.7)
  call otpm_set_value('stp_alkalinity', 2370.0 * 1024.5 * 1.0e-06)    ! alkalinity is in ueq/kg, converted to eq/m^3
  call otpm_set_value('sio2_const', 7.7e-03)                          ! mol/m^3
  call otpm_set_value('dic_global', 2.0)                              ! mol/m^3
  call otpm_set_value('alk_global', 2370.0 * rho_avg * 1.0e-06)       ! eq/m^3
  call otpm_set_value('po4_global', 2.17 * rho_avg * 1.0e-06)         ! mol/m^3
  call otpm_set_value('dop_global', 0.02 * rho_avg * 1.0e-06)         ! mol/m^3
  call otpm_set_value('o2_global', 170.0 * rho_avg * 1.0e-06)         ! mol/m^3
  call otpm_set_value('file_in', default_file_in)
  call otpm_set_value('file_out', default_file_out)
  call otpm_set_value('global_atmosphere', .false.)
  call otpm_set_value('init_global_atmosphere', .false.)
  call otpm_set_value('init', .false.)
  call otpm_set_value('n_2_p', 16.0)
  call otpm_set_value('c_2_p', 117.0)
  call otpm_set_value('o_2_p', 170.0)
  call otpm_set_value('o2_min', 4.0 * rho_avg * 1.0e-06)
  call otpm_set_value('bio_tau', 30.0 * sperd)
  call otpm_set_value('sigma', 0.67)
  call otpm_set_value('kappa', 1.0 / 0.5 / spery)
  call otpm_set_value('caco3_2_c', 0.07)

  call otpm_end_namelist(package_name, biotic(n)%name, check = .true., caller = caller_str)

!
!       create some sub-namelists
!

  call otpm_start_namelist(trim(package_name), trim(biotic(n)%name) // '+norm_remin',     &
       caller = caller_str, no_overwrite = .true., &
       check = .true.)

  call otpm_set_value('factor', 0.0)
  call otpm_set_value('coastal_only', .false.)
  call otpm_set_value('t_mask', t_mask, size(t_mask))
  call otpm_set_value('wlon', 0.0, index = 0)
  call otpm_set_value('elon', 0.0, index = 0)
  call otpm_set_value('slat', 0.0, index = 0)
  call otpm_set_value('nlat', 0.0, index = 0)

  call otpm_end_namelist(trim(package_name), trim(biotic(n)%name) // '+norm_remin', caller = caller_str)

  call otpm_start_namelist(trim(package_name), trim(biotic(n)%name) // '+no_caco3',     &
       caller = caller_str, no_overwrite = .true., &
       check = .true.)

  call otpm_set_value('factor', 0.0)
  call otpm_set_value('coastal_only', .false.)
  call otpm_set_value('t_mask', t_mask, size(t_mask))
  call otpm_set_value('wlon', 0.0, index = 0)
  call otpm_set_value('elon', 0.0, index = 0)
  call otpm_set_value('slat', 0.0, index = 0)
  call otpm_set_value('nlat', 0.0, index = 0)

  call otpm_end_namelist(trim(package_name), trim(biotic(n)%name) // '+no_caco3', caller = caller_str)

  call otpm_start_namelist(trim(package_name), trim(biotic(n)%name) // '+nut_depl',     &
       caller = caller_str, no_overwrite = .true., &
       check = .true.)

  call otpm_set_value('factor', 0.0)
  call otpm_set_value('coastal_only', .false.)
  call otpm_set_value('t_mask', t_mask, size(t_mask))
  call otpm_set_value('wlon', 0.0, index = 0)
  call otpm_set_value('elon', 0.0, index = 0)
  call otpm_set_value('slat', 0.0, index = 0)
  call otpm_set_value('nlat', 0.0, index = 0)

  call otpm_end_namelist(trim(package_name), trim(biotic(n)%name) // '+nut_depl', caller = caller_str)

  call otpm_start_namelist(trim(package_name), trim(biotic(n)%name) // '+r_bio_tau_a',     &
       caller = caller_str, no_overwrite = .true., &
       check = .true.)

  call otpm_set_value('factor', 0.0)
  call otpm_set_value('coastal_only', .false.)
  call otpm_set_value('t_mask', t_mask, size(t_mask))
  call otpm_set_value('wlon', 0.0, index = 0)
  call otpm_set_value('elon', 0.0, index = 0)
  call otpm_set_value('slat', 0.0, index = 0)
  call otpm_set_value('nlat', 0.0, index = 0)

  call otpm_end_namelist(trim(package_name), trim(biotic(n)%name) // '+r_bio_tau_a', caller = caller_str)

enddo  !} n

!
!       Check for any errors in the number of fields in the namelists for this package
!

good_list => otpm_get_string_array('/ocean_mod/GOOD/namelists/' // trim(package_name) // '/good_values',   &
     caller = trim(mod_name) // '(' // trim(sub_name) // ')')
if (associated(good_list)) then  !{
  call otpm_check_for_bad_fields('/ocean_mod/namelists/' // trim(package_name), good_list,       &
       caller = trim(mod_name) // '(' // trim(sub_name) // ')')
  deallocate(good_list)
else  !}{
  call mpp_error(FATAL,trim(error_header) // ' Empty "' // trim(package_name) // '" list')
endif  !}

return

end subroutine ocmip2_biotic_init  !}
! </SUBROUTINE> NAME="ocmip2_biotic_init"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_biotic_source">
!
! <DESCRIPTION>
!     compute the source terms for the BIOTICs, including boundary
!     conditions (not done in setvbc, to minimize number
!     of hooks required in MOM base code)
! </DESCRIPTION>
!

subroutine ocmip2_biotic_source(Thickness)  !{

type(ocean_thickness_type), intent(in) :: Thickness

!
!-----------------------------------------------------------------------
!     modules (have to come first)
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'ocmip2_biotic_source'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header =                                &
     '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

!
!-----------------------------------------------------------------------
!     arguments
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!     local variables
!-----------------------------------------------------------------------
!

integer :: i
integer :: j
integer :: k
integer :: n
integer :: ind_po4
integer :: ind_dop
integer :: ind_dic
integer :: ind_alk
integer :: ind_o2
integer :: km_c
logical :: used
integer :: day
integer :: month
integer :: year
integer :: hour
integer :: minute
integer :: second

!
! =====================================================================
!     begin executable code
! =====================================================================
!

!if (no_source) then  !{
  !write (stdout(),*) trim(warn_header), 'Skipping calculations'
  !return
!else  !}{
  !write (stdout(),*) ', trim(note_header), 'Doing calculations'
!endif  !}

!
!       get the model month
!

call get_date(time%model_time, year, month, day,                &
              hour, minute, second)

!
!-----------------------------------------------------------------------
!     calculate the source terms for BIOTICs
!-----------------------------------------------------------------------
!

!
!---------------------------------------------------------------------
!     calculate interpolated PO4_star
!---------------------------------------------------------------------
!

call time_interp_external(po4_star_id, time%model_time, po4_star_t)

!
!       Loop over multiple instances
!

do n = 1, instances  !{

  ind_po4 = biotic(n)%ind_po4
  ind_dic = biotic(n)%ind_dic
  ind_dop = biotic(n)%ind_dop
  ind_o2 = biotic(n)%ind_o2
  ind_alk = biotic(n)%ind_alk
  km_c = biotic(n)%km_c

!
!-----------------------------------------------------------------------
!     Production
!-----------------------------------------------------------------------
!

!
!       compute PO4 restoring term and correct for partial
!       production in the bottom box
!

  do k = 1, km_c  !{
    do j = jsc, jec  !{
      do i = isc, iec  !{
        if (t_prog(ind_po4)%field(i,j,k,taum1) .gt.             &
            po4_star_t(i,j,k) *                                 &
            biotic(n)%nut_depl%mask(i,j,month)) then  !{
          biotic(n)%jprod(i,j,k) =                              &
               (t_prog(ind_po4)%field(i,j,k,taum1) -            &
                po4_star_t(i,j,k) *                             &
                biotic(n)%nut_depl%mask(i,j,month)) *           &
               biotic(n)%r_bio_tau_a%mask(i,j,month) * grid%tmask(i,j,k)
        else  !}{
          biotic(n)%jprod(i,j,k) = 0.0
        endif  !}
      enddo  !} i
    enddo  !} j
  enddo  !} k

  do j = jsc, jec  !{
    do i = isc, iec  !{
      biotic(n)%jprod(i,j,km_c) = biotic(n)%jprod(i,j,km_c) *   &
           (min(grid%ht(i,j),biotic(n)%compensation_depth) - grid%zw(km_c-1)) /   &
           Thickness%dht(i,j,km_c,tau)
    enddo  !} i
  enddo  !} j

!
!-----------------------------------------------------------------------
!     Normal production (used to maintain PO4 concentrations assuming no
!       other changes)
!-----------------------------------------------------------------------
!

  if (biotic(n)%add_phosphate) then  !{

!
!       compute PO4 restoring term and correct for partial
!       production in the bottom box
!

    do k = 1, km_c  !{
      do j = jsc, jec  !{
        do i = isc, iec  !{
          if (t_prog(ind_po4)%field(i,j,k,taum1) .gt.           &
              po4_star_t(i,j,k)) then  !{
            biotic(n)%jprod_norm(i,j,k) =                       &
                 (t_prog(ind_po4)%field(i,j,k,taum1) -          &
                  po4_star_t(i,j,k)) *                          &
                 biotic(n)%r_bio_tau * grid%tmask(i,j,k)
          else  !}{
            biotic(n)%jprod_norm(i,j,k) = 0.0
          endif  !}
        enddo  !} i
      enddo  !} j
    enddo  !} k

    do j = jsc, jec  !{
      do i = isc, iec  !{
        biotic(n)%jprod_norm(i,j,km_c) =                        &
             biotic(n)%jprod_norm(i,j,km_c) *                   &
             (min(grid%ht(i,j),biotic(n)%compensation_depth) - grid%zw(km_c-1)) / &
             Thickness%dht(i,j,km_c,tau)
      enddo  !} i
    enddo  !} j
  endif  !}

!
!-----------------------------------------------------------------------
!     Particle flux
!-----------------------------------------------------------------------
!

  do j = jsc, jec  !{
    do i = isc, iec  !{
      biotic(n)%flux_poc(i,j) = (1.0 - biotic(n)%sigma) *       &
           biotic(n)%jprod(i,j,1) * Thickness%dht(i,j,1,tau)
    enddo  !} i
  enddo  !} j

  do k = 2, km_c  !{
    do j = jsc, jec  !{
      do i = isc, iec  !{
        biotic(n)%flux_poc(i,j) = biotic(n)%flux_poc(i,j) +     &
             (1.0 - biotic(n)%sigma) * biotic(n)%jprod(i,j,k) * &
             Thickness%dht(i,j,k,tau)
      enddo  !} i
    enddo  !} j
  enddo  !} k

!
!       calculate the flux at the base of each layer below the
!       compensation depth
!

  do k = km_c, nk  !{
    do j = jsc, jec  !{
      do i = isc, iec  !{
        biotic(n)%fc(i,j,k) = biotic(n)%flux_poc(i,j) *         &
             grid%tmask(i,j,k) *                                &
             (biotic(n)%norm_remin%mask(i,j,month) *            &
              biotic(n)%zforg(k) +                              &
              (1.0 - biotic(n)%norm_remin%mask(i,j,month)) *    &
              biotic(n)%zforg_alt(k))
      enddo  !} i
    enddo  !} j
  enddo  !} k

!
!-----------------------------------------------------------------------
!     Calcium Carbonate
!-----------------------------------------------------------------------
!

!
!       calculate the formation of CaCO3 above the compensation
!       depth
!

  do k = 1, km_c  !{
    do j = jsc, jec  !{
      do i = isc, iec  !{
        biotic(n)%jca(i,j,k) =                                  &
             -biotic(n)%caco3_2_c * biotic(n)%c_2_p *           &
             (1.0 - biotic(n)%sigma) *                          &
             biotic(n)%jprod(i,j,k) * grid%tmask(i,j,k) *       &
             biotic(n)%no_caco3%mask(i,j,month)
      enddo  !} i
    enddo  !} j
  enddo  !} k

!
!       calculate the flux of CaCO3 at compensation depth and at
!       bottom of each level below
!

  do j = jsc, jec  !{
    do i = isc, iec  !{
      biotic(n)%flux_caco3(i,j) =                               &
           biotic(n)%caco3_2_c * biotic(n)%c_2_p *              &
           biotic(n)%flux_poc(i,j) * biotic(n)%no_caco3%mask(i,j,month)
    enddo  !} i
  enddo  !} j

  do k=km_c,nk  !{
    do j = jsc, jec  !{
      do i = isc, iec  !{
        biotic(n)%fca(i,j,k) = biotic(n)%flux_caco3(i,j) *      &
             biotic(n)%zfca(k) * grid%tmask(i,j,k)
      enddo  !} i
    enddo  !} j
  enddo  !} k

!
!       calculate the dissolution of CaCO3 below the compensation
!       depth
!

  do j = jsc, jec  !{
    do i = isc, iec  !{
      biotic(n)%jca(i,j,km_c) = biotic(n)%jca(i,j,km_c) +       &
           (biotic(n)%flux_caco3(i,j) * grid%tmask(i,j,km_c) -  &
            biotic(n)%fca(i,j,km_c) * grid%tmask(i,j,km_c+1)) / &
           Thickness%dht(i,j,km_c,tau)
    enddo  !} i
  enddo  !} j

  do k = km_c + 1, nk - 1  !{
    do j = jsc, jec  !{
      do i = isc, iec  !{
        biotic(n)%jca(i,j,k) =                                  &
             (biotic(n)%fca(i,j,k-1) * grid%tmask(i,j,k) -      &
              biotic(n)%fca(i,j,k) * grid%tmask(i,j,k+1)) / Thickness%dht(i,j,k,tau)
      enddo  !} i
    enddo  !} j
  enddo  !} k

  do j = jsc, jec  !{
    do i = isc, iec  !{
      biotic(n)%jca(i,j,nk) = biotic(n)%fca(i,j,nk-1) *         &
           grid%tmask(i,j,nk) / Thickness%dht(i,j,nk,tau)
    enddo  !} i
  enddo  !} j

!
!-----------------------------------------------------------------------
!     PO4
!-----------------------------------------------------------------------
!

  if (biotic(n)%add_phosphate) then  !{
    do k = 1, km_c  !{
      do j = jsc, jec  !{
        do i = isc, iec  !{
          biotic(n)%jpo4_add(i,j,k) =                           &
               biotic(n)%jprod(i,j,k) - biotic(n)%jprod_norm(i,j,k)
        enddo  !} i
      enddo  !} j
    enddo  !} k
  endif  !}

  do k = 1, km_c  !{
    do j = jsc, jec  !{
      do i = isc, iec  !{
        biotic(n)%jpo4(i,j,k) = -biotic(n)%jprod(i,j,k) +       &
             biotic(n)%kappa * t_prog(ind_dop)%field(i,j,k,taum1)
      enddo  !} i
    enddo  !} j
  enddo  !} k

  do j = jsc, jec  !{
    do i = isc, iec  !{
      biotic(n)%jpo4(i,j,km_c) = biotic(n)%jpo4(i,j,km_c) +     &
           (biotic(n)%flux_poc(i,j) * grid%tmask(i,j,km_c) -    &
            biotic(n)%fc(i,j,km_c) * grid%tmask(i,j,km_c+1)) /  &
           Thickness%dht(i,j,km_c,tau)
    enddo  !} i
  enddo  !} j

  do k = km_c + 1, nk - 1  !{
    do j = jsc, jec  !{
      do i = isc, iec  !{
        biotic(n)%jpo4(i,j,k) =                                 &
             (biotic(n)%fc(i,j,k-1) * grid%tmask(i,j,k) -       &
              biotic(n)%fc(i,j,k) * grid%tmask(i,j,k+1)) /      &
             Thickness%dht(i,j,k,tau) +                              &
             biotic(n)%kappa * t_prog(ind_dop)%field(i,j,k,taum1)
      enddo  !} i
    enddo  !} j
  enddo  !} k

  do j = jsc, jec  !{
    do i = isc, iec  !{
      biotic(n)%jpo4(i,j,nk) =                                  &
           biotic(n)%fc(i,j,nk-1) * grid%tmask(i,j,nk) /        &
           Thickness%dht(i,j,nk,tau) +                               &
           biotic(n)%kappa * t_prog(ind_dop)%field(i,j,nk,taum1)
    enddo  !} i
  enddo  !} j

  do k = 1, nk  !{
    do j = jsc, jec  !{
      do i = isc, iec  !{
        t_prog(ind_po4)%th_tendency(i,j,k) =                         &
             t_prog(ind_po4)%th_tendency(i,j,k) +                    &
             biotic(n)%jpo4(i,j,k) * Thickness%dht(i,j,k,tau)
      enddo  !} i
    enddo  !} j
  enddo  !} k

  if (biotic(n)%add_phosphate) then  !{
    do k = 1, km_c  !{
      do j = jsc, jec  !{
        do i = isc, iec  !{
          t_prog(ind_po4)%th_tendency(i,j,k) =                       &
               t_prog(ind_po4)%th_tendency(i,j,k) +                  &
               biotic(n)%jpo4_add(i,j,k) * Thickness%dht(i,j,k,tau)
        enddo  !} i
      enddo  !} j
    enddo  !} k
  endif  !}

!
!-----------------------------------------------------------------------
!     DOP
!-----------------------------------------------------------------------
!

  do k = 1, km_c  !{
    do j = jsc, jec  !{
      do i = isc, iec  !{
        biotic(n)%jdop(i,j,k) =                                 &
             biotic(n)%sigma * biotic(n)%jprod(i,j,k) -         &
             biotic(n)%kappa * t_prog(ind_dop)%field(i,j,k,taum1)
      enddo  !} i
    enddo  !} j
  enddo  !} k

  do k = km_c + 1, nk  !{
    do j = jsc, jec  !{
      do i = isc, iec  !{
        biotic(n)%jdop(i,j,k) = -biotic(n)%kappa *              &
                                t_prog(ind_dop)%field(i,j,k,taum1)
      enddo  !} i
    enddo  !} j
  enddo  !} k

  do k = 1, nk  !{
    do j = jsc, jec  !{
      do i = isc, iec  !{
        t_prog(ind_dop)%th_tendency(i,j,k) =                         &
             t_prog(ind_dop)%th_tendency(i,j,k) +                    &
             biotic(n)%jdop(i,j,k) * Thickness%dht(i,j,k,tau)
      enddo  !} i
    enddo  !} j
  enddo  !} k

!
!-----------------------------------------------------------------------
!     O2
!-----------------------------------------------------------------------
!

  do k = 1, nk  !{
    do j = jsc, jec  !{
      do i = isc, iec  !{
        if (t_prog(ind_o2)%field(i,j,k,taum1) .gt.              &
            biotic(n)%o2_min) then  !{
          biotic(n)%jo2(i,j,k) = -biotic(n)%o_2_p *             &
                                 biotic(n)%jpo4(i,j,k)
        else  !}{
          biotic(n)%jo2(i,j,k) = 0.0
        endif  !}
      enddo  !} i
    enddo  !} j
  enddo  !} k

  do k = 1, nk  !{
    do j = jsc, jec  !{
      do i = isc, iec  !{
        t_prog(ind_o2)%th_tendency(i,j,k) =                     &
             t_prog(ind_o2)%th_tendency(i,j,k) +                &
             biotic(n)%jo2(i,j,k) * Thickness%dht(i,j,k,tau)
      enddo  !} i
    enddo  !} j
  enddo  !} k

!
!-----------------------------------------------------------------------
!     DIC
!-----------------------------------------------------------------------
!

  do k = 1, nk  !{
    do j = jsc, jec  !{
      do i = isc, iec  !{
        t_prog(ind_dic)%th_tendency(i,j,k) =                    &
             t_prog(ind_dic)%th_tendency(i,j,k) +               &
             (biotic(n)%c_2_p * biotic(n)%jpo4(i,j,k) +         &
              biotic(n)%jca(i,j,k)) * Thickness%dht(i,j,k,tau)
      enddo  !} i
    enddo  !} j
  enddo  !} k

!
!-----------------------------------------------------------------------
!     ALK
!-----------------------------------------------------------------------
!

  do k = 1, nk  !{
    do j = jsc, jec  !{
      do i = isc, iec  !{
        t_prog(ind_alk)%th_tendency(i,j,k) =                    &
             t_prog(ind_alk)%th_tendency(i,j,k) +               &
             (-biotic(n)%n_2_p * biotic(n)%jpo4(i,j,k) +        &
              2.0 * biotic(n)%jca(i,j,k)) * Thickness%dht(i,j,k,tau)
      enddo  !} i
    enddo  !} j
  enddo  !} k

enddo  !} n

!
!-----------------------------------------------------------------------
!       Save variables for diagnostics
!-----------------------------------------------------------------------
!

if (id_sc_co2 .gt. 0) then
  used = send_data(id_sc_co2, sc_co2(isc:iec,jsc:jec),          &
       time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
endif
if (id_sc_o2 .gt. 0) then
  used = send_data(id_sc_o2, sc_o2(isc:iec,jsc:jec),            &
       time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
endif
if (id_o2_sat .gt. 0) then
  used = send_data(id_o2_sat, o2_saturation(isc:iec,jsc:jec),   &
       time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
endif
if (id_kw_co2 .gt. 0) then
  used = send_data(id_kw_co2, kw_co2(isc:iec,jsc:jec),          &
       time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
endif
if (id_kw_o2 .gt. 0) then
  used = send_data(id_kw_o2, kw_o2(isc:iec,jsc:jec),            &
       time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
endif

do n = 1, instances  !{

  if (biotic(n)%id_csat .gt. 0) then
    used = send_data(biotic(n)%id_csat,                         &
         biotic(n)%csat(isc:iec,jsc:jec),                       &
         time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
  endif
  if (biotic(n)%id_csurf .gt. 0) then
    used = send_data(biotic(n)%id_csurf,                        &
         biotic(n)%csurf(isc:iec,jsc:jec),                      &
         time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
  endif
  if (biotic(n)%id_pco2surf .gt. 0) then
    used = send_data(biotic(n)%id_pco2surf,                     &
         biotic(n)%pco2surf(isc:iec,jsc:jec),                   &
         time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
  endif
  if (biotic(n)%id_dpco2 .gt. 0) then
    used = send_data(biotic(n)%id_dpco2,                        &
         biotic(n)%dpco2(isc:iec,jsc:jec),                      &
         time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
  endif
  if (biotic(n)%id_flux_poc .gt. 0) then
    used = send_data(biotic(n)%id_flux_poc,                     &
         biotic(n)%flux_poc(isc:iec,jsc:jec),                   &
         time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
  endif
  if (biotic(n)%id_flux_caco3 .gt. 0) then
    used = send_data(biotic(n)%id_flux_caco3,                   &
         biotic(n)%flux_caco3(isc:iec,jsc:jec),                 &
         time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
  endif
  if (biotic(n)%id_pco2atm .gt. 0) then
    used = send_data(biotic(n)%id_pco2atm,                      &
         biotic(n)%pco2atm(isc:iec,jsc:jec),                    &
         time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
  endif
  if (biotic(n)%id_htotal .gt. 0) then
    used = send_data(biotic(n)%id_htotal,                       &
         biotic(n)%htotal(isc:iec,jsc:jec,1),                   &
         time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
  endif

  if (.not. use_waterflux) then  !{

    if (biotic(n)%do_dic_virtual_flux) then  !{
      if (biotic(n)%id_vstf_dic .gt. 0) then
        used = send_data(biotic(n)%id_vstf_dic,                 &
             biotic(n)%vstf_dic(isc:iec,jsc:jec),               &
             time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
      endif
    endif  !}
    if (biotic(n)%do_alk_virtual_flux) then  !{
      if (biotic(n)%id_vstf_alk .gt. 0) then
        used = send_data(biotic(n)%id_vstf_alk,                 &
             biotic(n)%vstf_alk(isc:iec,jsc:jec),               &
             time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
      endif
    endif  !}
    if (biotic(n)%do_po4_virtual_flux) then  !{
      if (biotic(n)%id_vstf_po4 .gt. 0) then
        used = send_data(biotic(n)%id_vstf_po4,                 &
             biotic(n)%vstf_po4(isc:iec,jsc:jec),               &
             time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
      endif
    endif  !}
    if (biotic(n)%do_dop_virtual_flux) then  !{
      if (biotic(n)%id_vstf_dop .gt. 0) then
        used = send_data(biotic(n)%id_vstf_dop,                 &
             biotic(n)%vstf_dop(isc:iec,jsc:jec),               &
             time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
      endif
    endif  !}
    if (biotic(n)%do_o2_virtual_flux) then  !{
      if (biotic(n)%id_vstf_o2 .gt. 0) then
        used = send_data(biotic(n)%id_vstf_o2,                 &
             biotic(n)%vstf_o2(isc:iec,jsc:jec),               &
             time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
      endif
    endif  !}

  endif  !}

  if (biotic(n)%id_jprod .gt. 0) then
    used = send_data(biotic(n)%id_jprod,                        &
         biotic(n)%jprod(isc:iec,jsc:jec,:),                    &
         time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,:))
  endif
  if (biotic(n)%id_jca .gt. 0) then
    used = send_data(biotic(n)%id_jca,                          &
         biotic(n)%jca(isc:iec,jsc:jec,:),                      &
         time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,:))
  endif
  if (biotic(n)%id_jpo4 .gt. 0) then
    used = send_data(biotic(n)%id_jpo4,                         &
         biotic(n)%jpo4(isc:iec,jsc:jec,:),                     &
         time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,:))
  endif
  if (biotic(n)%id_jdop .gt. 0) then
    used = send_data(biotic(n)%id_jdop,                         &
         biotic(n)%jdop(isc:iec,jsc:jec,:),                     &
         time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,:))
  endif
  if (biotic(n)%id_jo2 .gt. 0) then
    used = send_data(biotic(n)%id_jo2,                          &
         biotic(n)%jo2(isc:iec,jsc:jec,:),                      &
         time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,:))
  endif
  if (biotic(n)%id_fc .gt. 0) then
    used = send_data(biotic(n)%id_fc,                           &
         biotic(n)%fc(isc:iec,jsc:jec,:),                       &
         time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,:))
  endif
  if (biotic(n)%id_fca .gt. 0) then
    used = send_data(biotic(n)%id_fca,                          &
         biotic(n)%fca(isc:iec,jsc:jec,:),                      &
         time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,:))
  endif

  if (biotic(n)%add_phosphate) then  !{

    if (biotic(n)%id_jpo4_add .gt. 0) then
      used = send_data(biotic(n)%id_jpo4_add,                   &
           biotic(n)%jpo4_add(isc:iec,jsc:jec,1:km_c),          &
           time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1:km_c))
    endif

  endif  !}

enddo  !} n

return

end subroutine  ocmip2_biotic_source  !}
! </SUBROUTINE> NAME="ocmip2_biotic_source"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_biotic_start">
!
! <DESCRIPTION>
! Initialize variables, read in namelists, calculate constants
! for a given run and allocate diagnostic arrays
! </DESCRIPTION>
!

subroutine ocmip2_biotic_start  !{

!
!-----------------------------------------------------------------------
!       modules (have to come first)
!-----------------------------------------------------------------------
!

use time_manager_mod, only: days_in_year, days_in_month,        &
     get_date, set_date
use time_interp_external_mod, only: init_external_field
use diag_manager_mod, only: register_diag_field, diag_axis_init
use fms_mod, only : read_data

!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'ocmip2_biotic_start'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header =                                &
     '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

!
!-----------------------------------------------------------------------
!       local variables
!-----------------------------------------------------------------------
!

!
!----------------------------------------------------------------------
!
!       Global values to apply the following inhibitions
!       and depletions
!
!  coastal_only : if true, then only apply the changes in
!                 coastal boxes
!  t_mask_len   : parameter giving the number of elements in
!                 the time mask per year (eg., 12 would
!                 imply monthly)
!  t_mask_array : logical array controlling whether to apply
!                 the following inhibitions and depletions to
!                 each time-period (true means set the masks,
!                 false means use the defaults everywhere)
!  num_reg      : number of regions
!  factor       : factor by which to scale the field
!               : in the selected regions
!  wlon : western longitude of region
!  elon : eastern longitude of region
!  slat : southern latitude of region
!  nlat : northern latitude of region
!  mask(imt,jmt)  : mask array (0.0 - alternate, 1.0 - normal)
!
!       Set up a mask array using wlon,elon,nlat,slat
!       (any box with its lon,lat inside the box bounded by
!       wlon,elon,nlat,slat value in mask set to factor).
!  
!----------------------------------------------------------------------
!

real                                    :: n_2_p
real                                    :: c_2_p
real                                    :: o_2_p
real                                    :: o2_min
real                                    :: bio_tau
real                                    :: r_bio_tau
real                                    :: sigma
real                                    :: kappa
real                                    :: compensation_depth
real                                    :: ca_remin_depth
real                                    :: caco3_2_c
character(len=fm_string_len)            :: file_in
character(len=fm_string_len)            :: file_out
integer                                 :: index
logical                                 :: init_global_atmosphere
logical                                 :: init
logical                                 :: global_atmosphere
logical                                 :: do_alk_virtual_flux
logical                                 :: do_dic_virtual_flux
logical                                 :: do_dop_virtual_flux
logical                                 :: do_o2_virtual_flux
logical                                 :: do_po4_virtual_flux
logical                                 :: add_phosphate
logical                                 :: soft_tissue_pump
real                                    :: martin_coeff
real                                    :: martin_coeff_alt
real                                    :: sio2_const
real                                    :: stp_temperature
real                                    :: stp_salinity
real                                    :: stp_alkalinity
real                                    :: sal_global
real                                    :: dic_global
real                                    :: alk_global
real                                    :: po4_global
real                                    :: dop_global
real                                    :: o2_global
integer                                 :: check
integer                                 :: day
real                                    :: days_in_this_year
integer                                 :: done
logical                                 :: found_error
integer                                 :: hour
integer                                 :: i
integer                                 :: j
integer                                 :: k
integer                                 :: l
integer                                 :: lun
integer                                 :: m
integer                                 :: n
integer                                 :: minute
integer                                 :: month
character(len=fm_field_name_len)        :: name
integer                                 :: num_days
integer                                 :: second
character(len=fm_field_name_len+1)      :: str
integer                                 :: year
logical                                 :: pert_first
integer                                 :: pert_start_day
integer                                 :: pert_start_hour
integer                                 :: pert_start_minute
integer                                 :: pert_start_month
integer                                 :: pert_start_second
integer                                 :: pert_start_year
real, allocatable, dimension(:)         :: zt_bounds
integer, dimension(:)                   :: tracer_axes_1(3)
integer                                 :: id_zt_1
integer                                 :: id_zt_1_bounds
character(len=256)                      :: caller_str
integer                                 :: len_w
integer                                 :: len_e
integer                                 :: len_s
integer                                 :: len_n


!
! =====================================================================
!       begin of executable code
! =====================================================================
!
!
!-----------------------------------------------------------------------
!       give info
!-----------------------------------------------------------------------
!

write(stdout(),*) 
write(stdout(),*) trim(note_header),                     &
                  'Starting ', trim(package_name), ' module'

!
!-----------------------------------------------------------------------
!     dynamically allocate the global BIOTIC arrays
!-----------------------------------------------------------------------
!

call allocate_arrays

!
!-----------------------------------------------------------------------
!       save the *global* namelist values
!-----------------------------------------------------------------------
!

caller_str = trim(mod_name) // '(' // trim(sub_name) // ')'

call otpm_start_namelist(package_name, '*global*', caller = caller_str)

atmpress_file      =  otpm_get_string ('atmpress_file', scalar = .true.)
atmpress_name      =  otpm_get_string ('atmpress_name', scalar = .true.)
pistonveloc_file   =  otpm_get_string ('pistonveloc_file', scalar = .true.)
pistonveloc_name   =  otpm_get_string ('pistonveloc_name', scalar = .true.)
seaicefract_file   =  otpm_get_string ('seaicefract_file', scalar = .true.)
seaicefract_name   =  otpm_get_string ('seaicefract_name', scalar = .true.)
po4_star_file      =  otpm_get_string ('po4_star_file', scalar = .true.)
po4_star_name      =  otpm_get_string ('po4_star_name', scalar = .true.)
pert_first         =  otpm_get_logical('pert_first', scalar = .true.)
pert_start_year    =  otpm_get_integer('pert_start_year', scalar = .true.)
pert_start_month   =  otpm_get_integer('pert_start_month', scalar = .true.)
pert_start_day     =  otpm_get_integer('pert_start_day', scalar = .true.)
pert_start_hour    =  otpm_get_integer('pert_start_hour', scalar = .true.)
pert_start_minute  =  otpm_get_integer('pert_start_minute', scalar = .true.)
pert_start_second  =  otpm_get_integer('pert_start_second', scalar = .true.)
control            =  otpm_get_logical('control', scalar = .true.)
historical         =  otpm_get_logical('historical', scalar = .true.)
scenario           =  otpm_get_string ('scenario', scalar = .true.)
htotal_scale_lo_in =  otpm_get_real   ('htotal_scale_lo_in', scalar = .true.)
htotal_scale_hi_in =  otpm_get_real   ('htotal_scale_hi_in', scalar = .true.)
htotal_in          =  otpm_get_real   ('htotal_in', scalar = .true.)
sal_global         =  otpm_get_real   ('sal_global', scalar = .true.)

call otpm_end_namelist(package_name, '*global*', caller = caller_str)
      
!
!       save pert_start_year into the biotic array
!

do n = 1, instances  !{
  biotic(n)%pert_start_year = pert_start_year
enddo  !} n

!
!-----------------------------------------------------------------------
!       Open up the PO4 file for restoring
!-----------------------------------------------------------------------
!

po4_star_id = init_external_field(po4_star_file,                &
                                  po4_star_name,                &
                                  domain = Domain%domain2d)
if (po4_star_id .eq. 0) then  !{
  call mpp_error(FATAL, trim(error_header) //                   &
       'Could not open po4_star file: ' //                      &
       trim(po4_star_file))
endif  !}

!
!-----------------------------------------------------------------------
!       Open up the files for boundary conditions
!-----------------------------------------------------------------------
!

atmpress_id = init_external_field(atmpress_file,                &
                                  atmpress_name,                &
                                  domain = Domain%domain2d)
if (atmpress_id .eq. 0) then  !{
  call mpp_error(FATAL, trim(error_header) //                   &
       'Could not open atmpress file: ' //                      &
       trim(atmpress_file))
endif  !}

pistonveloc_id = init_external_field(pistonveloc_file,          &
                                     pistonveloc_name,          &
                                     domain = Domain%domain2d)
if (pistonveloc_id .eq. 0) then  !{
  call mpp_error(FATAL, trim(error_header) //                   &
       'Could not open pistonveloc file: ' //                   &
       trim(pistonveloc_file))
endif  !}

seaicefract_id = init_external_field(seaicefract_file,          &
                                     seaicefract_name,          &
                                     domain = Domain%domain2d)
if (seaicefract_id .eq. 0) then  !{
  call mpp_error(FATAL, trim(error_header) //                   &
       'Could not open seaicefract file: ' //                   &
       trim(seaicefract_file))
endif  !}

!
!-----------------------------------------------------------------------
!       determine whether an equilibrium run
!-----------------------------------------------------------------------
!

equilibrium = scenario .eq. 'Equilibrium'

!
!-----------------------------------------------------------------------
!       check some things
!-----------------------------------------------------------------------
!

if (equilibrium) then  !{
  if (historical) then  !{
    call mpp_error(FATAL, trim(error_header) //                 &
         'Equilibrium scenario specified for historical run')
  endif  !}
else  !}{
  if (scenario .ne. "S650" .and.                                &
      scenario .ne. "CIS9") then  !{
    call mpp_error(FATAL, trim(error_header) //                 &
         'Bad scenario (' // trim(scenario) //                  &
         ') for historical or future run')
  endif  !}
endif  !}

!
!-----------------------------------------------------------------------
!       determine number of records for historical or future runs
!-----------------------------------------------------------------------
!

!if (historical) then  !{
  !num_records = (1999 - 1765 + 1) * 12
!elseif (equilibrium) then  !}{
  !num_records = 0
!elseif (scenario .eq. "S650") then  !}{
  !if (control) then  !{
    !num_records = (2299 - 2000 + 1) * 12
  !else  !}{
    !num_records = (2299 - 1990 + 1) * 12
  !endif  !}
!else  !}{
  !if (control) then  !{
    !num_records = (2299 - 2000 + 1) * 12
  !else  !}{
    !num_records = (2099 - 1990 + 1) * 12
  !endif  !}
!endif  !}

!
! set default values for htotal_scale bounds
!

htotal_scale_lo(:,:) = htotal_scale_lo_in
htotal_scale_hi(:,:) = htotal_scale_hi_in

!
!-----------------------------------------------------------------------
!       read in the namelists for each instance
!-----------------------------------------------------------------------
!

do n = 1, instances  !{

  call otpm_start_namelist(package_name, biotic(n)%name, caller = caller_str)

  biotic(n)%pco2atm_const          = otpm_get_real   ('pco2atm_const', scalar = .true.)
  biotic(n)%compensation_depth     = otpm_get_real   ('compensation_depth', scalar = .true.)
  biotic(n)%add_phosphate          = otpm_get_logical('add_phosphate', scalar = .true.)
  biotic(n)%do_dic_virtual_flux    = otpm_get_logical('do_dic_virtual_flux', scalar = .true.)
  biotic(n)%do_alk_virtual_flux    = otpm_get_logical('do_alk_virtual_flux', scalar = .true.)
  biotic(n)%do_po4_virtual_flux    = otpm_get_logical('do_po4_virtual_flux', scalar = .true.)
  biotic(n)%do_dop_virtual_flux    = otpm_get_logical('do_dop_virtual_flux', scalar = .true.)
  biotic(n)%do_o2_virtual_flux     = otpm_get_logical('do_o2_virtual_flux', scalar = .true.)
  biotic(n)%martin_coeff           = otpm_get_real   ('martin_coeff', scalar = .true.)
  biotic(n)%martin_coeff_alt       = otpm_get_real   ('martin_coeff_alt', scalar = .true.)
  biotic(n)%ca_remin_depth         = otpm_get_real   ('ca_remin_depth', scalar = .true.)
  biotic(n)%soft_tissue_pump       = otpm_get_logical('soft_tissue_pump', scalar = .true.)
  biotic(n)%stp_temperature        = otpm_get_real   ('stp_temperature', scalar = .true.)
  biotic(n)%stp_salinity           = otpm_get_real   ('stp_salinity', scalar = .true.)
  biotic(n)%stp_alkalinity         = otpm_get_real   ('stp_alkalinity', scalar = .true.)
  biotic(n)%sio2_const             = otpm_get_real   ('sio2_const', scalar = .true.)
  biotic(n)%file_in                = otpm_get_string ('file_in', scalar = .true.)
  biotic(n)%file_out               = otpm_get_string ('file_out', scalar = .true.)
  biotic(n)%global_atmosphere      = otpm_get_logical('global_atmosphere', scalar = .true.)
  biotic(n)%init_global_atmosphere = otpm_get_logical('init_global_atmosphere', scalar = .true.)
  biotic(n)%init                   = otpm_get_logical('init', scalar = .true.)
  biotic(n)%dic_global             = otpm_get_real   ('dic_global', scalar = .true.)
  biotic(n)%alk_global             = otpm_get_real   ('alk_global', scalar = .true.)
  biotic(n)%po4_global             = otpm_get_real   ('po4_global', scalar = .true.)
  biotic(n)%dop_global             = otpm_get_real   ('dop_global', scalar = .true.)
  biotic(n)%o2_global              = otpm_get_real   ('o2_global', scalar = .true.)
  biotic(n)%n_2_p                  = otpm_get_real   ('n_2_p', scalar = .true.)
  biotic(n)%c_2_p                  = otpm_get_real   ('c_2_p', scalar = .true.)
  biotic(n)%o_2_p                  = otpm_get_real   ('o_2_p', scalar = .true.)
  biotic(n)%o2_min                 = otpm_get_real   ('o2_min', scalar = .true.)
  biotic(n)%bio_tau                = otpm_get_real   ('bio_tau', scalar = .true.)
  biotic(n)%sigma                  = otpm_get_real   ('sigma', scalar = .true.)
  biotic(n)%kappa                  = otpm_get_real   ('kappa', scalar = .true.)
  biotic(n)%caco3_2_c              = otpm_get_real   ('caco3_2_c', scalar = .true.)

  call otpm_end_namelist(package_name, biotic(n)%name, caller = caller_str)

  biotic(n)%r_bio_tau = 1.0 / biotic(n)%bio_tau

  call otpm_start_namelist(trim(package_name), trim(biotic(n)%name) // '+norm_remin', caller = caller_str)

  biotic(n)%norm_remin%factor        =  otpm_get_real          ('factor', scalar = .true.)
  biotic(n)%norm_remin%coastal_only  =  otpm_get_logical       ('coastal_only', scalar = .true.)
  biotic(n)%norm_remin%wlon          => otpm_get_real_array    ('wlon')
  biotic(n)%norm_remin%elon          => otpm_get_real_array    ('elon')
  biotic(n)%norm_remin%slat          => otpm_get_real_array    ('slat')
  biotic(n)%norm_remin%nlat          => otpm_get_real_array    ('nlat')
  biotic(n)%norm_remin%t_mask        => otpm_get_logical_array ('t_mask')

  call otpm_end_namelist(trim(package_name), trim(biotic(n)%name) // '+norm_remin', caller = caller_str)

  call otpm_start_namelist(trim(package_name), trim(biotic(n)%name) // '+no_caco3', caller = caller_str)

  biotic(n)%no_caco3%factor          =  otpm_get_real          ('factor', scalar = .true.)
  biotic(n)%no_caco3%coastal_only    =  otpm_get_logical       ('coastal_only', scalar = .true.)
  biotic(n)%no_caco3%wlon            => otpm_get_real_array    ('wlon')
  biotic(n)%no_caco3%elon            => otpm_get_real_array    ('elon')
  biotic(n)%no_caco3%slat            => otpm_get_real_array    ('slat')
  biotic(n)%no_caco3%nlat            => otpm_get_real_array    ('nlat')
  biotic(n)%no_caco3%t_mask          => otpm_get_logical_array ('t_mask')

  call otpm_end_namelist(trim(package_name), trim(biotic(n)%name) // '+no_caco3', caller = caller_str)

  call otpm_start_namelist(trim(package_name), trim(biotic(n)%name) // '+nut_depl', caller = caller_str)

  biotic(n)%nut_depl%factor          =  otpm_get_real          ('factor', scalar = .true.)
  biotic(n)%nut_depl%coastal_only    =  otpm_get_logical       ('coastal_only', scalar = .true.)
  biotic(n)%nut_depl%wlon            => otpm_get_real_array    ('wlon')
  biotic(n)%nut_depl%elon            => otpm_get_real_array    ('elon')
  biotic(n)%nut_depl%slat            => otpm_get_real_array    ('slat')
  biotic(n)%nut_depl%nlat            => otpm_get_real_array    ('nlat')
  biotic(n)%nut_depl%t_mask          => otpm_get_logical_array ('t_mask')

  call otpm_end_namelist(trim(package_name), trim(biotic(n)%name) // '+nut_depl', caller = caller_str)

  call otpm_start_namelist(trim(package_name), trim(biotic(n)%name) // '+r_bio_tau_a', caller = caller_str)

  biotic(n)%r_bio_tau_a%factor       =  otpm_get_real          ('factor', scalar = .true.)
  biotic(n)%r_bio_tau_a%coastal_only =  otpm_get_logical       ('coastal_only', scalar = .true.)
  biotic(n)%r_bio_tau_a%wlon         => otpm_get_real_array    ('wlon')
  biotic(n)%r_bio_tau_a%elon         => otpm_get_real_array    ('elon')
  biotic(n)%r_bio_tau_a%slat         => otpm_get_real_array    ('slat')
  biotic(n)%r_bio_tau_a%nlat         => otpm_get_real_array    ('nlat')
  biotic(n)%r_bio_tau_a%t_mask       => otpm_get_logical_array ('t_mask')

  call otpm_end_namelist(trim(package_name), trim(biotic(n)%name) // '+r_bio_tau_a', caller = caller_str)


enddo  !} n

!
!-----------------------------------------------------------------------
!     calculate the index for the box containing the compensation depth
!-----------------------------------------------------------------------
!

km_c_max = 0
do n = 1, instances  !{
  call locate(grid%zw, nk, biotic(n)%compensation_depth,        &
              biotic(n)%km_c, nearest = .true.)
  if (grid%zw(biotic(n)%km_c) .lt.                              &
      biotic(n)%compensation_depth) then  !{
    biotic(n)%km_c = biotic(n)%km_c + 1
  endif  !}

  write (stdout(),*) trim(note_header),                         &
                     'The compensation depth for instance ',    &
                     n, ', ', biotic(n)%compensation_depth,     &
                     ' m, occurs in box ', biotic(n)%km_c ,     &
                     ' between depths ',                        &
                     grid%zw(biotic(n)%km_c-1),                 &
                     ' m and ', grid%zw(biotic(n)%km_c), ' m'
  km_c_max = max(km_c_max, biotic(n)%km_c)
enddo  !} n

!
!-----------------------------------------------------------------------
!       read in the norm_remin namelist data
!-----------------------------------------------------------------------
!

do n = 1, instances  !{

!
!       Check some things
!

  if (associated(biotic(n)%norm_remin%wlon)) then  !{
    len_w = size(biotic(n)%norm_remin%wlon)
  else  !}{
    len_w = 0
  endif  !}
  if (associated(biotic(n)%norm_remin%elon)) then  !{
    len_e = size(biotic(n)%norm_remin%elon)
  else  !}{
    len_e = 0
  endif  !}
  if (associated(biotic(n)%norm_remin%slat)) then  !{
    len_s = size(biotic(n)%norm_remin%slat)
  else  !}{
    len_s = 0
  endif  !}
  if (associated(biotic(n)%norm_remin%nlat)) then  !{
    len_n = size(biotic(n)%norm_remin%nlat)
  else  !}{
    len_n = 0
  endif  !}

  if (len_e .ne. len_w .or. len_w .ne. len_s .or. len_s .ne. len_n) then  !{
    call mpp_error(FATAL, trim(error_header) // ' Region sizes are not equal for ' // trim(biotic(n)%name))
  endif  !}

  if (size(biotic(n)%norm_remin%t_mask) .ne. 12) then  !{
    call mpp_error(FATAL, trim(error_header) // ' t_mask size is not 12 for ' // trim(biotic(n)%name))
  endif  !}

!
!       set all of the values to the default
!

  biotic(n)%norm_remin%mask(:,:,:) = 1.0

  if (len_w .gt. 0) then  !{

    write (stdout(),*)
    write (stdout(),*) trim(note_header), 'Process norm_remin array for ', trim(biotic(n)%name)
    write (stdout(),*)

!
!       set values for this time-level
!

    done = 0
    do l = 1, 12  !{
      if (biotic(n)%norm_remin%t_mask(l)) then  !{
        if (done .eq. 0) then  !{

!
!       set the values via the input values, saving this time index
!       afterwards
!

          write (stdout(),*) 'Assigning month ', l
          call set_array(biotic(n)%norm_remin%mask(:,:,l),        &
                  isd, ied, jsd, jed, grid%xt, grid%yt, grid%kmt, &
                  len_w, biotic(n)%norm_remin%wlon,               &
                  biotic(n)%norm_remin%elon,                      &
                  biotic(n)%norm_remin%slat,                      &
                  biotic(n)%norm_remin%nlat,                      &
                  biotic(n)%norm_remin%factor, 1.0,               &
                  'Normal remineralization', biotic(n)%norm_remin%coastal_only)
          done = l
        else  !}{

!
!       Duplicate the values for a previous time-level
!
          write (stdout(),*) 'Duplicating month ', done, ' as ', l
          biotic(n)%norm_remin%mask(:,:,l) = biotic(n)%norm_remin%mask(:,:,done)
        endif  !}
      endif  !}
    enddo  !} l
  endif  !}

enddo  !} n

!
!-----------------------------------------------------------------------
!       read in the no_caco3 namelist data
!-----------------------------------------------------------------------
!

do n = 1, instances  !{

!
!       Check some things
!

  if (associated(biotic(n)%no_caco3%wlon)) then  !{
    len_w = size(biotic(n)%no_caco3%wlon)
  else  !}{
    len_w = 0
  endif  !}
  if (associated(biotic(n)%no_caco3%elon)) then  !{
    len_e = size(biotic(n)%no_caco3%elon)
  else  !}{
    len_e = 0
  endif  !}
  if (associated(biotic(n)%no_caco3%slat)) then  !{
    len_s = size(biotic(n)%no_caco3%slat)
  else  !}{
    len_s = 0
  endif  !}
  if (associated(biotic(n)%no_caco3%nlat)) then  !{
    len_n = size(biotic(n)%no_caco3%nlat)
  else  !}{
    len_n = 0
  endif  !}

  if (len_e .ne. len_w .or. len_w .ne. len_s .or. len_s .ne. len_n) then  !{
    call mpp_error(FATAL, trim(error_header) // ' Region sizes are not equal for ' // trim(biotic(n)%name))
  endif  !}

  if (size(biotic(n)%no_caco3%t_mask) .ne. 12) then  !{
    call mpp_error(FATAL, trim(error_header) // ' t_mask size is not 12 for ' // trim(biotic(n)%name))
  endif  !}

!
!       set all of the values to the default
!

  biotic(n)%no_caco3%mask(:,:,:) = 1.0

  if (len_w .gt. 0) then  !{

    write (stdout(),*)
    write (stdout(),*) trim(note_header), 'Process no_caco3 array for ', trim(biotic(n)%name)
    write (stdout(),*)

!
!       set values for this time-level
!

    done = 0
    do l = 1, 12  !{
      if (biotic(n)%no_caco3%t_mask(l)) then  !{
        if (done .eq. 0) then  !{

!
!       set the values via the input values, saving this time index
!       afterwards
!

          write (stdout(),*) 'Assigning month ', l
          call set_array(biotic(n)%no_caco3%mask(:,:,l),                &
                  isd, ied, jsd, jed, grid%xt, grid%yt, grid%kmt,       &
                  len_w, biotic(n)%no_caco3%wlon,                       &
                  biotic(n)%no_caco3%elon,                              &
                  biotic(n)%no_caco3%slat,                              &
                  biotic(n)%no_caco3%nlat,                              &
                  biotic(n)%no_caco3%factor, 1.0,                       &
                  'Carbonate inhibition', biotic(n)%no_caco3%coastal_only)
          done = l
        else  !}{

!
!       Duplicate the values for a previous time-level
!
          write (stdout(),*) 'Duplicating month ', done, ' as ', l
          biotic(n)%no_caco3%mask(:,:,l) = biotic(n)%no_caco3%mask(:,:,done)
        endif  !}
      endif  !}
    enddo  !} l
  endif  !}

enddo  !} n

!
!-----------------------------------------------------------------------
!       read in the nut_depl namelist data
!-----------------------------------------------------------------------
!

do n = 1, instances  !{

!
!       Check some things
!

  if (associated(biotic(n)%nut_depl%wlon)) then  !{
    len_w = size(biotic(n)%nut_depl%wlon)
  else  !}{
    len_w = 0
  endif  !}
  if (associated(biotic(n)%nut_depl%elon)) then  !{
    len_e = size(biotic(n)%nut_depl%elon)
  else  !}{
    len_e = 0
  endif  !}
  if (associated(biotic(n)%nut_depl%slat)) then  !{
    len_s = size(biotic(n)%nut_depl%slat)
  else  !}{
    len_s = 0
  endif  !}
  if (associated(biotic(n)%nut_depl%nlat)) then  !{
    len_n = size(biotic(n)%nut_depl%nlat)
  else  !}{
    len_n = 0
  endif  !}

  if (len_e .ne. len_w .or. len_w .ne. len_s .or. len_s .ne. len_n) then  !{
    call mpp_error(FATAL, trim(error_header) // ' Region sizes are not equal for ' // trim(biotic(n)%name))
  endif  !}

  if (size(biotic(n)%nut_depl%t_mask) .ne. 12) then  !{
    call mpp_error(FATAL, trim(error_header) // ' t_mask size is not 12 for ' // trim(biotic(n)%name))
  endif  !}

!
!       set all of the values to the default
!

  biotic(n)%nut_depl%mask(:,:,:) = 1.0

  if (len_w .gt. 0) then  !{

    write (stdout(),*)
    write (stdout(),*) trim(note_header), 'Process nut_depl array for ', trim(biotic(n)%name)
    write (stdout(),*)

!
!       set values for this time-level
!

    done = 0
    do l = 1, 12  !{
      if (biotic(n)%nut_depl%t_mask(l)) then  !{
        if (done .eq. 0) then  !{

!
!       set the values via the input values, saving this time index
!       afterwards
!

          write (stdout(),*) 'Assigning month ', l
          call set_array(biotic(n)%nut_depl%mask(:,:,l),                &
                  isd, ied, jsd, jed, grid%xt, grid%yt, grid%kmt,       &
                  len_w, biotic(n)%nut_depl%wlon,                       &
                  biotic(n)%nut_depl%elon,                              &
                  biotic(n)%nut_depl%slat,                              &
                  biotic(n)%nut_depl%nlat,                              &
                  biotic(n)%nut_depl%factor, 1.0,                       &
                  'Nutrient depletion', biotic(n)%nut_depl%coastal_only)
          done = l
        else  !}{

!
!       Duplicate the values for a previous time-level
!
          write (stdout(),*) 'Duplicating month ', done, ' as ', l
          biotic(n)%nut_depl%mask(:,:,l) = biotic(n)%nut_depl%mask(:,:,done)
        endif  !}
      endif  !}
    enddo  !} l
  endif  !}

enddo  !} n

!
!-----------------------------------------------------------------------
!       read in the r_bio_tau_a namelist data
!-----------------------------------------------------------------------
!

do n = 1, instances  !{

!
!       Check some things
!

  if (associated(biotic(n)%r_bio_tau_a%wlon)) then  !{
    len_w = size(biotic(n)%r_bio_tau_a%wlon)
  else  !}{
    len_w = 0
  endif  !}
  if (associated(biotic(n)%r_bio_tau_a%elon)) then  !{
    len_e = size(biotic(n)%r_bio_tau_a%elon)
  else  !}{
    len_e = 0
  endif  !}
  if (associated(biotic(n)%r_bio_tau_a%slat)) then  !{
    len_s = size(biotic(n)%r_bio_tau_a%slat)
  else  !}{
    len_s = 0
  endif  !}
  if (associated(biotic(n)%r_bio_tau_a%nlat)) then  !{
    len_n = size(biotic(n)%r_bio_tau_a%nlat)
  else  !}{
    len_n = 0
  endif  !}

  if (len_e .ne. len_w .or. len_w .ne. len_s .or. len_s .ne. len_n) then  !{
    call mpp_error(FATAL, trim(error_header) // ' Region sizes are not equal for ' // trim(biotic(n)%name))
  endif  !}

  if (size(biotic(n)%r_bio_tau_a%t_mask) .ne. 12) then  !{
    call mpp_error(FATAL, trim(error_header) // ' t_mask size is not 12 for ' // trim(biotic(n)%name))
  endif  !}

!
!       set all of the values to the default
!

  biotic(n)%r_bio_tau_a%mask(:,:,:) = 1.0

  if (len_w .gt. 0) then  !{

    write (stdout(),*)
    write (stdout(),*) trim(note_header), 'Process r_bio_tau_a array for ', trim(biotic(n)%name)
    write (stdout(),*)

!
!       set values for this time-level
!

    done = 0
    do l = 1, 12  !{
      if (biotic(n)%r_bio_tau_a%t_mask(l)) then  !{
        if (done .eq. 0) then  !{

!
!       set the values via the input values, saving this time index
!       afterwards
!

          write (stdout(),*) 'Assigning month ', l
          call set_array(biotic(n)%r_bio_tau_a%mask(:,:,l),             &
                  isd, ied, jsd, jed, grid%xt, grid%yt, grid%kmt,       &
                  len_w, biotic(n)%r_bio_tau_a%wlon,                    &
                  biotic(n)%r_bio_tau_a%elon,                           &
                  biotic(n)%r_bio_tau_a%slat,                           &
                  biotic(n)%r_bio_tau_a%nlat,                           &
                  biotic(n)%r_bio_tau_a%factor, 1.0,                    &
                  'Primary production limitation', biotic(n)%r_bio_tau_a%coastal_only)
          done = l
        else  !}{

!
!       Duplicate the values for a previous time-level
!
          write (stdout(),*) 'Duplicating month ', done, ' as ', l
          biotic(n)%r_bio_tau_a%mask(:,:,l) = biotic(n)%r_bio_tau_a%mask(:,:,done)
        endif  !}
      endif  !}
    enddo  !} l
  endif  !}

enddo  !} n

!
!       multiply by the restoring factor
!

do n = 1, instances  !{
  biotic(n)%r_bio_tau_a%mask(:,:,:) =                                &
       biotic(n)%r_bio_tau * biotic(n)%r_bio_tau_a%mask(:,:,:)
enddo  !} n

!
!-----------------------------------------------------------------------
!     initialize special arrays for remineralization
!-----------------------------------------------------------------------
!

do n = 1, instances  !{
  do k = 1, nk  !{
    biotic(n)%zforg(k) =                                        &
         (grid%zw(k) /                                          &
          biotic(n)%compensation_depth) ** (-biotic(n)%martin_coeff)
    biotic(n)%zfca(k) =                                         &
         exp(-(grid%zw(k) -                                     &
               biotic(n)%compensation_depth) / biotic(n)%ca_remin_depth)
  enddo  !} k
enddo  !} n

!
!-----------------------------------------------------------------------
!     initialize special arrays for alternate remineralization
!-----------------------------------------------------------------------
!

do n = 1, instances  !{
  do k = 1, nk  !{
    biotic(n)%zforg_alt(k) = (grid%zw(k) /                      &
         biotic(n)%compensation_depth) ** (-biotic(n)%martin_coeff_alt)
  enddo  !} k
enddo  !} n

!
!-----------------------------------------------------------------------
!     do some calendar calculation for the next section,
!      calculate the number of days in this year as well as
!      those that have already passed
!-----------------------------------------------------------------------
!

call get_date(time%model_time,                                  &
              year, month, day, hour, minute, second)
num_days = days_in_year(time%model_time)
days_in_this_year = 0.0
do m = 1, month - 1
  days_in_this_year = days_in_this_year +                       &
                     days_in_month(set_date(year, m, 1))
enddo
days_in_this_year = days_in_this_year + day - 1 + hour/24.0 +   &
                   minute/1440.0 + second/86400.0

if (.not. equilibrium) then  !{

  do n = 1, instances  !{

!     
!-----------------------------------------------------------------------
!     if this is the beginning of a perturbation run then calculate 
!         pert_starttime_mod (is in real years)
!     set cumulative input to zero for historical run, otherwise
!         initialize from the input file
!-----------------------------------------------------------------------
!

    if (pert_first) then  !{

      if (month .ne. pert_start_month .or.                      &
          day .ne. pert_start_day .or.                          &
          hour .ne. pert_start_hour .or.                        &
          minute .ne. pert_start_minute .or.                    &
          second .ne. pert_start_second) then  !{
        call mpp_error(FATAL, trim(error_header) //             &
             'Perturbation start time is at different' //       &
             ' point in year for instance ' // trim(biotic(n)%name))
      endif  !}

      biotic(n)%pert_start_year_model = year

      write (stdout(),*)
      write (stdout(),*) trim(note_header),                     &
                         'Instance ', trim(biotic(n)%name)
      write (stdout(),*) 'This is the beginning of a perturbation run'
      write (stdout(),*) '   perturbation start year: ',        &
                         biotic(n)%pert_start_year
      write (stdout(),*) '   perturbation start year (model): ',&
                         biotic(n)%pert_start_year_model
      write (stdout(),*) '   itt : ', time%itt

!if (historical) then  !{
!biotic_cumul(:,:,:) = 0.0
!else  !}{
!write(stdout(),*)
!write(stdout(),*) trim(note_header), 'Getting cumulative',&
!' input from restart file...'
!call getunit (lun, biotic_cumul_infile,                &
!'unformatted sequential rewind')

!read (lun) dummy1, dummy2, dummy3, biotic_cumul
!call relunit (lun)
!write (stdout(),*) '  read cumulative input from:'
!write (stdout(),*) '     pert_time         =', dummy3
!write (stdout(),*) '     pert_starttime    =', dummy1
!write (stdout(),*) '     pert_starttime_mod=', dummy2
!write (stdout(),*)
!endif  !}

!
!-----------------------------------------------------------------------
!     otherwise read info and cumulative input from a file
!       check also consistency of input 
!-----------------------------------------------------------------------
!

    else  !}{

      write(stdout(),*)
      write(stdout(),*) trim(note_header),                      &
           'Getting pert values from restart file',             &
           ' for instance ', trim(biotic(n)%name)

      call read_data(biotic(n)%file_in,                         &
           'pert_start_year', biotic(n)%pert_start_year,        &
           domain=Domain%domain2d, timelevel=1)
      call read_data(biotic(n)%file_in,                         &
           'pert_start_year_model',                             &
           biotic(n)%pert_start_year_model,                     &
           domain=Domain%domain2d, timelevel=1)
      call read_data(biotic(n)%file_in,                         &
           'pert_time', biotic(n)%pert_time,                    &
           domain=Domain%domain2d, timelevel=1)
!biotic_cumul
!read (lun) gmm_pco2atm,                                        &
!mgf_biotic_dic,                                &
!mgvf_biotic_dic, mgvf_biotic_alk,              &
!mgv_biotic_dic, mgv_biotic_alk,                &
!mgs_biotic_dic, mgs_biotic_alk,                &
!mg_biotic_pco2surf, mg_biotic_dpco2,           &
!gcf_biotic_dic,                                &
!gcvf_biotic_dic, gcvf_biotic_alk
        
      check = year - biotic(n)%pert_start_year_model +          &
              biotic(n)%pert_start_year +                       &
              days_in_this_year / num_days

      if (check .ne. biotic(n)%pert_time) then  !{

        write (stdout(),*) 
        write (stdout(),*) trim(error_header)
        write (stdout(),*) 'Model time does not match ',        &
                           ' perturbation time!'
        write (stdout(),*) '     pert_start_year            = ',&
                           biotic(n)%pert_start_year
        write (stdout(),*) '     pert_start_year_model      = ',&
                           biotic(n)%pert_start_year_model
        write (stdout(),*) '     pert_time from restart    = ', &
                           biotic(n)%pert_time
        write (stdout(),*) '     pert_time calc from model = ', &
                           check
        call mpp_error(FATAL, trim(error_header) // 'Time mismatch')

      endif  !}

      write (stdout(),*) 
      write (stdout(),*) trim(note_header)
      write (stdout(),*) '     pert_start_year       = ',       &
                         biotic(n)%pert_start_year
      write (stdout(),*) '     pert_start_year_model = ',       &
                         biotic(n)%pert_start_year_model
      write (stdout(),*) '     pert_time             = ',       &
                         biotic(n)%pert_time
      write (stdout(),*) 

!
!-----------------------------------------------------------------------
!     end of check for pert_first
!-----------------------------------------------------------------------
!

    endif  !}

  enddo  !} n

endif  !}

!
!-----------------------------------------------------------------------
!       read in additional information for a restart
!-----------------------------------------------------------------------
!

write(stdout(),*)

do n = 1, instances  !{

  if (biotic(n)%global_atmosphere) then  !{

    if (biotic(n)%init_global_atmosphere) then  !{

      write (stdout(),*) trim(note_header),                     &
           'Setting initial atmosphere for instance ',          &
           trim(biotic(n)%name)

      biotic(n)%total_atm_co2(tau) = total_atm_mass *           &
           (co2_mol_wgt / atm_mol_wgt) *                        &
           (biotic(n)%pco2atm_const * 1.0e+03)
      biotic(n)%total_atm_co2(taup1) = biotic(n)%total_atm_co2(tau)

    else  !}{

      write (stdout(),*) trim(note_header),                     &
           'Reading atmospheric information for instance ',     &
           trim(biotic(n)%name)

      call read_data(biotic(n)%file_in,                         &
           'total_atm_co2', biotic(n)%total_atm_co2(tau),       &
           domain=Domain%domain2d, timelevel=1)
      call read_data(biotic(n)%file_in,                         &
           'total_atm_co2', biotic(n)%total_atm_co2(taup1),     &
           domain=Domain%domain2d, timelevel=2)

    endif  !}

    write (stdout(),'(/1x,a,es16.9,a,a)')                       &
          'Total atmospheric mass of CO2 (tau  ) = ',           &
          biotic(n)%total_atm_co2(tau  ) * 1.0e+15,             &
          ' (Pg) for instance ', trim(biotic(n)%name)
    write (stdout(),'(/1x,a,es16.9,a,a)')                       &
          'Total atmospheric mass of CO2 (taup1) = ',           &
          biotic(n)%total_atm_co2(taup1) * 1.0e+15,             &
          ' (Pg) for instance ', trim(biotic(n)%name)

  endif  !}

  if (biotic(n)%init) then  !{

    write (stdout(),*) trim(note_header),                       &
         'Initializing instance ', trim(biotic(n)%name)

    biotic(n)%htotal(:,:,:) = htotal_in

    if (.not. use_waterflux) then  !{
      if (biotic(n)%do_dic_virtual_flux .or.                    &
          biotic(n)%do_alk_virtual_flux .or.                    &
          biotic(n)%do_po4_virtual_flux .or.                    &
          biotic(n)%do_dop_virtual_flux .or.                    &
          biotic(n)%do_o2_virtual_flux) then  !{
        biotic(n)%sal_global_wrk = 0.0
        biotic(n)%global_wrk_duration = 0.0
      endif  !}
      if (biotic(n)%do_dic_virtual_flux) then  !{
        biotic(n)%dic_global_wrk = 0.0
      endif  !}
      if (biotic(n)%do_alk_virtual_flux) then  !{
        biotic(n)%alk_global_wrk = 0.0
      endif  !}
      if (biotic(n)%do_po4_virtual_flux) then  !{
        biotic(n)%po4_global_wrk = 0.0
      endif  !}
      if (biotic(n)%do_dop_virtual_flux) then  !{
        biotic(n)%dop_global_wrk = 0.0
      endif  !}
      if (biotic(n)%do_o2_virtual_flux) then  !{
        biotic(n)%o2_global_wrk = 0.0
      endif  !}
    endif  !}

  else  !}{

    write (stdout(),*) trim(note_header),                       &
         'Reading additional information for instance ',        &
         ': Initializing instance ', trim(biotic(n)%name)

    call read_data(biotic(n)%file_in,                           &
         'htotal', biotic(n)%htotal,                            &
         domain=Domain%domain2d, timelevel=1)

    if (.not. use_waterflux) then  !{

      if (biotic(n)%do_dic_virtual_flux .or.                    &
          biotic(n)%do_alk_virtual_flux .or.                    &
          biotic(n)%do_po4_virtual_flux .or.                    &
          biotic(n)%do_dop_virtual_flux .or.                    &
          biotic(n)%do_o2_virtual_flux) then  !{
        call read_data(biotic(n)%file_in,                       &
             'sal_global', biotic(n)%sal_global,                &
             domain=Domain%domain2d, timelevel=1)
        call read_data(biotic(n)%file_in,                       &
             'sal_global_wrk', biotic(n)%sal_global_wrk,        &
             domain=Domain%domain2d, timelevel=1)
        call read_data(biotic(n)%file_in,                       &
             'global_wrk_duration',                             &
             biotic(n)%global_wrk_duration,                     &
             domain=Domain%domain2d, timelevel=1)
      endif  !}
      if (biotic(n)%do_dic_virtual_flux) then  !{
        call read_data(biotic(n)%file_in,                       &
             'dic_global', biotic(n)%dic_global,                &
             domain=Domain%domain2d, timelevel=1)
        call read_data(biotic(n)%file_in,                       &
             'dic_global_wrk', biotic(n)%dic_global_wrk,        &
             domain=Domain%domain2d, timelevel=1)
      endif  !}
      if (biotic(n)%do_alk_virtual_flux) then  !{
        call read_data(biotic(n)%file_in,                       &
             'alk_global_wrk', biotic(n)%alk_global_wrk,        &
             domain=Domain%domain2d, timelevel=1)
        call read_data(biotic(n)%file_in,                       &
             'alk_global', biotic(n)%alk_global,                &
             domain=Domain%domain2d, timelevel=1)
      endif  !}
      if (biotic(n)%do_po4_virtual_flux) then  !{
        call read_data(biotic(n)%file_in,                       &
             'po4_global_wrk', biotic(n)%po4_global_wrk,        &
             domain=Domain%domain2d, timelevel=1)
        call read_data(biotic(n)%file_in,                       &
             'po4_global', biotic(n)%po4_global,                &
             domain=Domain%domain2d, timelevel=1)
      endif  !}
      if (biotic(n)%do_dop_virtual_flux) then  !{
        call read_data(biotic(n)%file_in,                       &
             'dop_global_wrk', biotic(n)%dop_global_wrk,        &
             domain=Domain%domain2d, timelevel=1)
        call read_data(biotic(n)%file_in,                       &
             'dop_global', biotic(n)%dop_global,                &
             domain=Domain%domain2d, timelevel=1)
      endif  !}
      if (biotic(n)%do_o2_virtual_flux) then  !{
        call read_data(biotic(n)%file_in,                       &
             'o2_global_wrk', biotic(n)%o2_global_wrk,          &
             domain=Domain%domain2d, timelevel=1)
        call read_data(biotic(n)%file_in,                       &
             'o2_global', biotic(n)%o2_global,                  &
             domain=Domain%domain2d, timelevel=1)
      endif  !}

    endif  !}

  endif  !}

  if (.not. use_waterflux) then  !{
    if (biotic(n)%do_dic_virtual_flux .or.                      &
        biotic(n)%do_alk_virtual_flux .or.                      &
        biotic(n)%do_po4_virtual_flux .or.                      &
        biotic(n)%do_dop_virtual_flux .or.                      &
        biotic(n)%do_o2_virtual_flux) then  !{
      write (stdout(),'(/1x,a,es16.9,a,a)')                     &
            'Annual, global, surface mean salinity = ',         &
            biotic(n)%sal_global, ' (PSU) for instance ',       &
            trim(biotic(n)%name)
    endif  !}
    if (biotic(n)%do_dic_virtual_flux) then  !{
      write (stdout(),'(1x,a,es16.9,a,a)')                      &
            'Annual, global, surface mean DIC      = ',         &
            biotic(n)%dic_global, ' (mol/m^3) for instance ',   &
            trim(biotic(n)%name)
    endif  !}
    if (biotic(n)%do_alk_virtual_flux) then  !{
      write (stdout(),'(1x,a,es16.9,a,a)')                      &
            'Annual, global, surface mean ALK      = ',         &
            biotic(n)%alk_global, ' (eq/m^3) for instance ',    &
            trim(biotic(n)%name)
    endif  !}
    if (biotic(n)%do_po4_virtual_flux) then  !{
      write (stdout(),'(1x,a,es16.9,a,a)')                      &
            'Annual, global, surface mean PO4      = ',         &
            biotic(n)%po4_global, ' (mol/m^3) for instance ',   &
            trim(biotic(n)%name)
    endif  !}
    if (biotic(n)%do_dop_virtual_flux) then  !{
      write (stdout(),'(1x,a,es16.9,a,a)')                      &
            'Annual, global, surface mean DOP      = ',         &
            biotic(n)%dop_global, ' (mol/m^3) for instance ',   &
            trim(biotic(n)%name)
    endif  !}
    if (biotic(n)%do_o2_virtual_flux) then  !{
      write (stdout(),'(1x,a,es16.9,a,a)')                      &
            'Annual, global, surface mean O2       = ',         &
            biotic(n)%o2_global, ' (mol/m^3) for instance ',    &
            trim(biotic(n)%name)
    endif  !}

  endif  !}

enddo  !} n



!
!-----------------------------------------------------------------------
!
!       initialize some arrays which are held constant for this
!       simulation
!
!-----------------------------------------------------------------------
!

do n = 1, instances  !{
  biotic(n)%sio2(:,:,:) = biotic(n)%sio2_const
enddo  !} n

!
!       only do the following for an equilibrium run
!

if (equilibrium) then  !{

  write(stdout(),*) 
  write(stdout(),*) trim(note_header),                          &
                    'Equilibrium run'
  write(stdout(),*) 
  do n = 1, instances  !{
    biotic(n)%pco2atm(:,:) = biotic(n)%pco2atm_const
  enddo  !} n

!
!       only do the following for a non-equilibrium run
!

else  !}{

  if (control) then  !{
    write(stdout(),*) 
    write(stdout(),*) trim(note_header),                        &
                      'Control run: Equilibrium atmosphere'
    write(stdout(),*) 
    do n = 1, instances  !{
      biotic(n)%pco2atm(:,:) = biotic(n)%pco2atm_const
    enddo  !} n
  endif  !}

  write(stdout(),*) 
  write(stdout(),*) trim(note_header),                          &
                    'Scenario: ', trim(scenario), ' run'
  write(stdout(),*) 

endif  !}

!
!-----------------------------------------------------------------------
!     Set up analyses
!-----------------------------------------------------------------------
!

!
!       set the axes for analysis for a 1-level grid
!

allocate( zt_bounds(nk + 1) )

zt_bounds(1) = Grid%zt(1) - Grid%dzt(1) * 0.5
do k = 2, Grid%nk+1  !{
  zt_bounds(k)=zt_bounds(k-1) + Grid%dzt(k-1)
enddo  !} k

id_zt_1_bounds = diag_axis_init(                                &
     'zt_1_edges_' // trim(Grid%name),                          &
     zt_bounds(1:2), 'meters', 'z', 'tcell depth edges',        &
     direction=-1, set_name='ocean')

id_zt_1 = diag_axis_init(                                       &
     'zt_1_' // trim(Grid%name),                                &
     Grid%zt(1:1), 'meters', 'z', 'tcell depth',                &
     edges = id_zt_1_bounds, direction=-1, set_name='ocean')

tracer_axes_1(1:2) = grid%tracer_axes(1:2)
tracer_axes_1(3) = id_zt_1

deallocate( zt_bounds )

!
!       register the fields
!

id_sc_co2 = register_diag_field('ocean_model',                  &
     'sc_co2', grid%tracer_axes(1:2),                           &
     Time%model_time, 'Schmidt number - CO2', ' ',              &
     missing_value = -1.0e+10)

id_sc_o2 = register_diag_field('ocean_model',                   &
     'sc_o2', grid%tracer_axes(1:2),                            &
     Time%model_time, 'Schmidt number - O2', ' ',               &
     missing_value = -1.0e+10)

id_o2_sat = register_diag_field('ocean_model',                  &
     'o2_saturation', grid%tracer_axes(1:2),                    &
     Time%model_time, 'O2 saturation', ' ',                     &
     missing_value = -1.0e+10)

id_kw_co2 = register_diag_field('ocean_model',                  &
     'kw_co2', grid%tracer_axes(1:2),                           &
     Time%model_time, 'Piston velocity - CO2', ' ',             &
     missing_value = -1.0e+10)

id_kw_o2 = register_diag_field('ocean_model',                   &
     'kw_o2', grid%tracer_axes(1:2),                            &
     Time%model_time, 'Piston velocity - O2', ' ',              &
     missing_value = -1.0e+10)

do n = 1, instances  !{

  if (instances .eq. 1) then  !{
    str = ' '
  else  !}{
    str = '_' // biotic(n)%name
  endif  !}

  biotic(n)%id_csat = register_diag_field('ocean_model',        &
       'csat'//str, grid%tracer_axes(1:2),                      &
       Time%model_time, 'CO2* air', ' ',                        &
       missing_value = -1.0e+10)

  biotic(n)%id_csurf = register_diag_field('ocean_model',       &
       'csurf'//str, grid%tracer_axes(1:2),                     &
       Time%model_time, 'CO2* water', ' ',                      &
       missing_value = -1.0e+10)

  biotic(n)%id_pco2surf = register_diag_field('ocean_model',    &
       'pco2surf'//str, grid%tracer_axes(1:2),                  &
       Time%model_time, 'Oceanic pCO2', ' ',                    &
       missing_value = -1.0e+10)

  biotic(n)%id_dpco2 = register_diag_field('ocean_model',       &
       'dpco2'//str, grid%tracer_axes(1:2),                     &
       Time%model_time, 'Delta pCO2', ' ',                      &
       missing_value = -1.0e+10)

  biotic(n)%id_flux_poc = register_diag_field('ocean_model',    &
       'flux_poc'//str, grid%tracer_axes(1:2),                  &
       Time%model_time, 'POC flux', ' ',                        &
       missing_value = -1.0e+10)

  biotic(n)%id_flux_caco3 = register_diag_field('ocean_model',  &
       'flux_caco3'//str, grid%tracer_axes(1:2),                &
       Time%model_time, 'CaCO3 flux', ' ',                      &
       missing_value = -1.0e+10)

  biotic(n)%id_pco2atm = register_diag_field('ocean_model',     &
       'pco2atm'//str, grid%tracer_axes(1:2),                   &
       Time%model_time, 'Atmospheric CO_2', ' ',                &
       missing_value = -1.0e+10)

  biotic(n)%id_htotal = register_diag_field('ocean_model',      &
       'htotal'//str, tracer_axes_1(1:3),                       &
       Time%model_time, 'H+ ion concentration', ' ',            &
       missing_value = -1.0e+10)

  if (.not. use_waterflux) then  !{

    if (biotic(n)%do_dic_virtual_flux) then  !{
      biotic(n)%id_vstf_dic = register_diag_field('ocean_model',&
           'vstf_dic'//str, grid%tracer_axes(1:2),              &
           Time%model_time, 'Virtual DIC flux', ' ',            &
           missing_value = -1.0e+10)
    endif  !}
    if (biotic(n)%do_alk_virtual_flux) then  !{
      biotic(n)%id_vstf_alk = register_diag_field('ocean_model',&
           'vstf_alk'//str, grid%tracer_axes(1:2),              &
           Time%model_time, 'Virtual ALK flux', ' ',            &
           missing_value = -1.0e+10)
    endif  !}
    if (biotic(n)%do_po4_virtual_flux) then  !{
      biotic(n)%id_vstf_po4 = register_diag_field('ocean_model',&
           'vstf_po4'//str, grid%tracer_axes(1:2),              &
           Time%model_time, 'Virtual PO4 flux', ' ',            &
           missing_value = -1.0e+10)
    endif  !}
    if (biotic(n)%do_dop_virtual_flux) then  !{
      biotic(n)%id_vstf_dop = register_diag_field('ocean_model',&
           'vstf_dop'//str, grid%tracer_axes(1:2),              &
           Time%model_time, 'Virtual DOP flux', ' ',            &
           missing_value = -1.0e+10)
    endif  !}
    if (biotic(n)%do_o2_virtual_flux) then  !{
      biotic(n)%id_vstf_o2 = register_diag_field('ocean_model', &
           'vstf_o2'//str, grid%tracer_axes(1:2),               &
           Time%model_time, 'Virtual O2 flux', ' ',             &
           missing_value = -1.0e+10)
    endif  !}

  endif  !}

  biotic(n)%id_jprod = register_diag_field('ocean_model',       &
       'jprod'//str, grid%tracer_axes(1:3),                     &
       Time%model_time, 'Restoring production', ' ',            &
       missing_value = -1.0e+10)

  biotic(n)%id_jca = register_diag_field('ocean_model',         &
       'jca'//str, grid%tracer_axes(1:3),                       &
       Time%model_time, 'CaCO3 change', ' ',                    &
       missing_value = -1.0e+10)

  biotic(n)%id_jpo4 = register_diag_field('ocean_model',        &
       'jpo4'//str, grid%tracer_axes(1:3),                      &
       Time%model_time, 'PO4 source', ' ',                      &
       missing_value = -1.0e+10)

  biotic(n)%id_jdop = register_diag_field('ocean_model',        &
       'jdop'//str, grid%tracer_axes(1:3),                      &
       Time%model_time, 'DOP source', ' ',                      &
       missing_value = -1.0e+10)

  biotic(n)%id_jo2 = register_diag_field('ocean_model',         &
       'jo2'//str, grid%tracer_axes(1:3),                       &
       Time%model_time, 'O2 source', ' ',                       &
       missing_value = -1.0e+10)

  biotic(n)%id_fc = register_diag_field('ocean_model',          &
       'fc'//str, grid%tracer_axes(1:3),                        &
       Time%model_time, 'POP change', ' ',                      &
       missing_value = -1.0e+10)

  biotic(n)%id_fca = register_diag_field('ocean_model',         &
       'fca'//str, grid%tracer_axes(1:3),                       &
       Time%model_time, 'CaCO3 change', ' ',                    &
       missing_value = -1.0e+10)

  if (biotic(n)%add_phosphate) then  !{

    biotic(n)%id_jpo4_add = register_diag_field('ocean_model',  &
         'jpo4_add'//str, grid%tracer_axes(1:3),                &
         Time%model_time, 'Additional phosphate', ' ',          &
         missing_value = -1.0e+10)

  endif  !}

enddo  !} n

!
!-----------------------------------------------------------------------
!     give info
!-----------------------------------------------------------------------
!

write(stdout(),*)
write(stdout(),*) trim(note_header), 'Tracer runs initialized'
write(stdout(),*)

return

end subroutine  ocmip2_biotic_start  !}
! </SUBROUTINE> NAME="ocmip2_biotic_start"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_biotic_tracer">
!
! <DESCRIPTION>
!     Perform things that should be done in tracer, but are done here
! in order to minimize the number of hooks necessary in the MOM4 basecode
! </DESCRIPTION>
!

subroutine ocmip2_biotic_tracer  !{

!
!-----------------------------------------------------------------------
!     modules (have to come first)
!-----------------------------------------------------------------------
!

use mpp_mod, only : mpp_sum

!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'ocmip2_biotic_tracer'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header =                                &
     '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

!
!-----------------------------------------------------------------------
!     arguments
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!     local variables
!-----------------------------------------------------------------------
!

integer :: i
integer :: j
integer :: n
real    :: total_alkalinity
real    :: total_dic
real    :: total_dop
real    :: total_o2
real    :: total_phosphate
real    :: total_salinity
real    :: total_temp

!
!-----------------------------------------------------------------------
!     accumulate global annual means
!-----------------------------------------------------------------------
!

if (.not. use_waterflux) then  !{

  do n = 1, instances  !{

    if (biotic(n)%do_dic_virtual_flux .or.                      &
        biotic(n)%do_alk_virtual_flux .or.                      &
        biotic(n)%do_po4_virtual_flux .or.                      &
        biotic(n)%do_dop_virtual_flux .or.                      &
        biotic(n)%do_o2_virtual_flux) then  !{
      do j = jsc, jec  !{
        do i = isc, iec  !{
          biotic(n)%sal_global_wrk = biotic(n)%sal_global_wrk + &
               t_prog(indsal)%field(i,j,1,tau) *                &
               grid%tmask(i,j,1) * grid%dat(i,j) * dtts
        enddo  !} i
      enddo  !} j
    endif  !}

    if (biotic(n)%do_dic_virtual_flux) then  !{
      do j = jsc, jec  !{
        do i = isc, iec  !{
          biotic(n)%dic_global_wrk = biotic(n)%dic_global_wrk + &
               t_prog(biotic(n)%ind_dic)%field(i,j,1,tau) *     &
               grid%tmask(i,j,1) * grid%dat(i,j) * dtts
        enddo  !} i
      enddo  !} j
    endif  !}
    if (biotic(n)%do_alk_virtual_flux) then  !{
      do j = jsc, jec  !{
        do i = isc, iec  !{
          biotic(n)%alk_global_wrk = biotic(n)%alk_global_wrk + &
               t_prog(biotic(n)%ind_alk)%field(i,j,1,tau) *     &
               grid%tmask(i,j,1) * grid%dat(i,j) * dtts
        enddo  !} i
      enddo  !} j
    endif  !}
    if (biotic(n)%do_po4_virtual_flux) then  !{
      do j = jsc, jec  !{
        do i = isc, iec  !{
          biotic(n)%po4_global_wrk = biotic(n)%po4_global_wrk + &
               t_prog(biotic(n)%ind_po4)%field(i,j,1,tau) *     &
               grid%tmask(i,j,1) * grid%dat(i,j) * dtts
        enddo  !} i
      enddo  !} j
    endif  !}
    if (biotic(n)%do_dop_virtual_flux) then  !{
      do j = jsc, jec  !{
        do i = isc, iec  !{
          biotic(n)%dop_global_wrk = biotic(n)%dop_global_wrk + &
               t_prog(biotic(n)%ind_dop)%field(i,j,1,tau) *     &
               grid%tmask(i,j,1) * grid%dat(i,j) * dtts
        enddo  !} i
      enddo  !} j
    endif  !}
    if (biotic(n)%do_o2_virtual_flux) then  !{
      do j = jsc, jec  !{
        do i = isc, iec  !{
          biotic(n)%o2_global_wrk = biotic(n)%o2_global_wrk + &
               t_prog(biotic(n)%ind_o2)%field(i,j,1,tau) *     &
               grid%tmask(i,j,1) * grid%dat(i,j) * dtts
        enddo  !} i
      enddo  !} j
    endif  !}

  enddo  !} n

  do n = 1, instances  !{
    if (biotic(n)%do_dic_virtual_flux .or.                      &
        biotic(n)%do_alk_virtual_flux .or.                      &
        biotic(n)%do_po4_virtual_flux .or.                      &
        biotic(n)%do_dop_virtual_flux .or.                      &
        biotic(n)%do_o2_virtual_flux) then  !{
      call mpp_sum(biotic(n)%sal_global_wrk)
    endif  !}
    if (biotic(n)%do_dic_virtual_flux) then  !{
      call mpp_sum(biotic(n)%dic_global_wrk)
    endif  !}
    if (biotic(n)%do_alk_virtual_flux) then  !{
      call mpp_sum(biotic(n)%alk_global_wrk)
    endif  !}
    if (biotic(n)%do_po4_virtual_flux) then  !{
      call mpp_sum(biotic(n)%po4_global_wrk)
    endif  !}
    if (biotic(n)%do_dop_virtual_flux) then  !{
      call mpp_sum(biotic(n)%dop_global_wrk)
    endif  !}
    if (biotic(n)%do_o2_virtual_flux) then  !{
      call mpp_sum(biotic(n)%o2_global_wrk)
    endif  !}
    biotic(n)%global_wrk_duration = biotic(n)%global_wrk_duration + dtts
  enddo  !} n

!
!----------------------------------------------------------------------
!       calculate global means of at the end of the year
!----------------------------------------------------------------------
!

  if (end_of_year) then  !{

    do n = 1, instances  !{
      if (biotic(n)%do_dic_virtual_flux .or.                    &
          biotic(n)%do_alk_virtual_flux .or.                    &
          biotic(n)%do_po4_virtual_flux .or.                    &
          biotic(n)%do_dop_virtual_flux .or.                    &
          biotic(n)%do_o2_virtual_flux) then  !{
        biotic(n)%sal_global = biotic(n)%sal_global_wrk /       &
             biotic(n)%global_wrk_duration / grid%tcella(1)
      endif  !}
      if (biotic(n)%do_dic_virtual_flux) then  !{
        biotic(n)%dic_global = biotic(n)%dic_global_wrk /       &
             biotic(n)%global_wrk_duration / grid%tcella(1)
      endif  !}
      if (biotic(n)%do_alk_virtual_flux) then  !{
        biotic(n)%alk_global = biotic(n)%alk_global_wrk /       &
             biotic(n)%global_wrk_duration / grid%tcella(1)
      endif  !}
      if (biotic(n)%do_po4_virtual_flux) then  !{
        biotic(n)%po4_global = biotic(n)%po4_global_wrk /       &
             biotic(n)%global_wrk_duration / grid%tcella(1)
      endif  !}
      if (biotic(n)%do_dop_virtual_flux) then  !{
        biotic(n)%dop_global = biotic(n)%dop_global_wrk /       &
             biotic(n)%global_wrk_duration / grid%tcella(1)
      endif  !}
      if (biotic(n)%do_o2_virtual_flux) then  !{
        biotic(n)%o2_global = biotic(n)%o2_global_wrk /         &
             biotic(n)%global_wrk_duration / grid%tcella(1)
      endif  !}
    enddo  !} n

  !
  !     reset work variables to zero
  !

    do n = 1, instances  !{
      if (biotic(n)%do_dic_virtual_flux .or.                    &
          biotic(n)%do_alk_virtual_flux .or.                    &
          biotic(n)%do_po4_virtual_flux .or.                    &
          biotic(n)%do_dop_virtual_flux .or.                    &
          biotic(n)%do_o2_virtual_flux) then  !{
        biotic(n)%sal_global_wrk = 0.0
        biotic(n)%global_wrk_duration = 0.0
      endif  !}
      if (biotic(n)%do_dic_virtual_flux) then  !{
        biotic(n)%dic_global_wrk = 0.0
      endif  !}
      if (biotic(n)%do_alk_virtual_flux) then  !{
        biotic(n)%alk_global_wrk = 0.0
      endif  !}
      if (biotic(n)%do_po4_virtual_flux) then  !{
        biotic(n)%po4_global_wrk = 0.0
      endif  !}
      if (biotic(n)%do_dop_virtual_flux) then  !{
        biotic(n)%dop_global_wrk = 0.0
      endif  !}
      if (biotic(n)%do_o2_virtual_flux) then  !{
        biotic(n)%o2_global_wrk = 0.0
      endif  !}
    enddo  !} n

  endif  !}

endif  !}

return

end subroutine  ocmip2_biotic_tracer  !}
! </SUBROUTINE> NAME="ocmip2_biotic_tracer"


!#######################################################################
! <SUBROUTINE NAME="read_c14atm">
!
! <DESCRIPTION>
!     Reads temporal history of atmospheric C-14 (permil)
!
!     ARGUMENT LIST -
!
!     Note: Variable type is given in square brackets (below)
!     (r-real, i-integer, l-logical, c-character; s-scaler, a-array).
!
!     INPUT:
!
!      none
!
!     OUTPUT: 
!
!     [ia] - nc14rec    =  3-member array with the number of records
!                         of atmospheric C-14 in each of the 3
!                         latitudinal bands listed below.
!
!     [ra] - yrc14rec   =  Sequential list of times (in decimal years)
!                         for when atmospheric C-14 data is available
!
!     [ra] - atmc14rec  =  Corresponding sequential list of atmospheric 
!                         C-14 (permil), in 3 latitudinal bands:
!                         (1) 90S - 20S,
!                         (2) 20S - 20N, and
!                         (3) 20N - 90N.
!                         This record is from Enting et al. (1994).
!     
!     Reference:
!
!     Enting, I.G., T. M. L. Wigley, M. Heimann, 1994. Future Emissions 
!     and concentrations of carbon dioxide: key ocean / atmosphere / 
!     land analyses, CSIRO Aust. Div. Atmos. Res. Tech. Pap. No. 31, 
!     118 pp.
!
!     James ORR, LSCE/CEA-CNRS, Saclay, France, 17 April 1999
! </DESCRIPTION>
!
!_ ---------------------------------------------------------------------
!_ 

SUBROUTINE read_c14atm (nc14rec, yrc14rec, atmc14rec)  !{

use mpp_io_mod, only : mpp_close, mpp_open, MPP_RDONLY,         &
     MPP_ASCII, MPP_SINGLE

INTEGER maxrec, nzon
PARAMETER (maxrec=300, nzon=3)

INTEGER, intent(out)    :: nc14rec(nzon)
real, intent(out)       :: yrc14rec(maxrec,nzon)
real, intent(out)       :: atmc14rec(maxrec,nzon)

integer :: irec, il
integer, dimension(3)   :: lun

!     Open one FILE for each of 3 latitudinal bands
!     ---------------------------------------------

call mpp_open(lun(1), 'c14sth.dat', action=MPP_RDONLY,          &
              form=MPP_ASCII, threading=MPP_SINGLE)
call mpp_open(lun(2), 'c14equ.dat', action=MPP_RDONLY,          &
              form=MPP_ASCII, threading=MPP_SINGLE)
call mpp_open(lun(3), 'c14nth.dat', action=MPP_RDONLY,          &
              form=MPP_ASCII, threading=MPP_SINGLE)

!
!-----------------------------------------------------------------------
!       Get atmospheric C-14 data:
!-----------------------------------------------------------------------

do il=1,3

!
!         Skip top 4 descriptor lines

  read(lun(il),'(///)')

  do irec=1,maxrec
    read(lun(il),*,err=11,end=12) yrc14rec(irec,il), atmc14rec(irec,il)
    nc14rec(il) = irec
  end do

11 continue
  write (stdout(),*) 'Error in reading input file: unit = ', lun(il)
  write (stdout(),*) 'Sequential record number (irec) = ', irec

12 continue
  WRITE(stdout(),*) 'Atm. C-14: No. of entries in lat. band '   &
           , il, ': ', nc14rec(il)
  call mpp_close(lun(il))
end do

return
end subroutine  read_c14atm  !}
! </SUBROUTINE> NAME="read_c14atm"


!#######################################################################
! <SUBROUTINE NAME="read_co2atm">
!
! <DESCRIPTION>
!     Read temporal history of atmospheric CO2 (uatm)
!
!     Argument list -
!
!     Note: Variable TYPE is given in square brackets (below)
!     (r-REAL, i-INTEGER, l-LOGICAL, c-CHARACTER; s-scaler, a-array).
!
!     INPUT:
!
!     [cs] - futr_scen  = IPCC future scenario: either S350, S450, S550,
!                        S650, S750, DS450, or DS550 from Enting et al.
!                        (1994), or CIS9 signifying c-IS92A for
!                        IPCC (2000) run.  From 1765-1990.5, it 
!                        doesn't matter which scenario you use, i.e., 
!                        atmospheric CO2 will be the same (from a
!                        spline fit to Siple Ice core and Mauna Loa
!                        data.  Subsequently, atmospheric CO2 is 
!                        different, according to the choice given above.
!
!     OUTPUT: 
!
!     [is] - nco2rec    =  Number of records (years) for atmospheric CO2
!                         from historical (splco2.dat) plus
!                         future (stab.dat) records
!
!     [ra] - yrco2rec   =  sequential list of times (in decimal years)
!                         for WHEN atmospheric CO2 data is available
!
!     [ra] - atmco2rec  =  corresponding sequential list of atmospheric 
!                         co2 (ppm).
!
!                         This record is from Enting et al. (1994).
!     
!     Reference:
!
!     Enting, I.G., T. M. L. Wigley, M. Heimann, 1994. Future emissions 
!     and concentrations of carbon dioxide: key ocean / atmosphere / 
!     land analyses, CSIRO Aust. Div. Atmos. Res., Tech. Pap. No. 31, 
!     118 pp.
! 
!     James Orr, LSCE/CEA-CNRS, Saclay, France, 17 April 1999
! </DESCRIPTION>
!
!_ ---------------------------------------------------------------------
!_ 

SUBROUTINE read_co2atm(futr_scen, nco2rec, yrco2rec, atmco2rec)  !{

use mpp_io_mod, only : mpp_close, mpp_open, MPP_RDONLY,         &
     MPP_ASCII, MPP_SINGLE

INTEGER maxrec, nmxr
PARAMETER (maxrec=1200, nmxr=700)

CHARACTER(len=4), intent(in)    :: futr_scen
integer , intent(out)           :: nco2rec
real, intent(out)               :: yrco2rec(maxrec), atmco2rec(maxrec)

INTEGER :: lun(3), irec
INTEGER is, ifuture, ireadf
INTEGER nsipl, nstab

REAL futco2(nmxr,8)
REAL dummy

CHARACTER*4 ipcc_scen(8)

!     Note that the 1st 7 future scenarios are in FILE "stab.dat";
!     the last scenario "CIS9" is short for "CIS92A" (see "cis92a.dat").

DATA ipcc_scen/'S350', 'S450', 'S550', 'S650', 'S750',          &
               'DS45', 'DS55', 'CIS9'/


!
!     Determine index for future scenario       

ifuture = 0
DO is=1,8
  IF (futr_scen(1:4) .EQ. ipcc_scen(is))THEN
    ifuture = is
  ENDIF
END DO

!
!     Check to be sure that chosen scenario is from allowed IPCC list

IF (ifuture .eq. 0) then
    WRITE (stdout(),*)                                          &
         'Improper 1st argument for read_co2atm.f: ', futr_scen(1:4)
    WRITE (stdout(),'(a,7(1x,a))') 'You must chose from ',      &
         (ipcc_scen(is), is=1,8)
    WRITE (stdout(),*)                                          &
         ' For OCMIP-2, S650 and CIS9 are the ONLY accetable choices'
    STOP
ELSE
    WRITE (stdout(),'(2a)') 'You have chosen IPCC scenario ',   &
         ipcc_scen(ifuture)
END IF 

!
!     OPEN FILE
!     ---------

call mpp_open (lun(1), 'splco2.dat', action=MPP_RDONLY,         &
               form=MPP_ASCII, threading=MPP_SINGLE)
call mpp_open (lun(2), 'stab.dat',   action=MPP_RDONLY,         &
               form=MPP_ASCII, threading=MPP_SINGLE)
call mpp_open (lun(3), 'cis92a.dat', action=MPP_RDONLY,         &
               form=MPP_ASCII, threading=MPP_SINGLE)

!
!     skip over 1st six descriptor lines
!
!-----------------------------------------------------------------------
!       get atmospheric CO2 data
!-----------------------------------------------------------------------

WRITE(stdout(),*) '  '
WRITE(stdout(),*) '--------------------------------------------------'
WRITE(stdout(),*) 'Atm. CO2 from spline fit to Siple-Mauna Loa record'
WRITE(stdout(),*) ' And IPCC scenario ', futr_scen(1:4)
WRITE(stdout(),*) '--------------------------------------------------'

!
!     --------------------------------------------------------------
!     READ either historical or future co2 concentrations, depending
!     upon the value of ireadf, which is enabled (set to 1) during
!     the READ operation, at the END of the historical FILE (IF
!     ifuture > 0)
!     --------------------------------------------------------------

READ(lun(1),'(////)')
ireadf = 0
nsipl = 0
nstab = 0
DO irec=1,maxrec
  210 continue
  IF (ireadf .EQ. 0) THEN

!           Read from splco2.dat (historical emissions)

    READ(lun(1),*,ERR=220,END=220) yrco2rec(irec), atmco2rec(irec)
    nsipl = nsipl + 1
  ELSE IF (ireadf .EQ. 1) THEN

!           Read from stab.dat (future atm CO2 scenario)

    READ(lun(2),*,ERR=222,END=222) yrco2rec(irec)               &
                                 , (futco2(irec-nsipl,is), is=1,7)
    atmco2rec(irec) = futco2(irec-nsipl,ifuture)
    nstab = nstab + 1
  ELSE IF (ireadf .EQ. 2) THEN

!           Read from stab.dat (future atm CO2 scenario)

    READ(lun(3),*,ERR=222,END=222) yrco2rec(irec)               &
                                 , futco2(irec-nsipl,8), dummy
    atmco2rec(irec) = futco2(irec-nsipl,8)
    nstab = nstab + 1
  ENDIF

  GO TO 221

!       When end of historical co2 reached, turn on read ability
!       (ireadf>0) for future file, go back to read one or the other, 
!       THEN continue reading chosen future file until it ends

  220 continue
  IF (ifuture .gt. 0 .AND. ifuture .le. 7) THEN
    ireadf = 1

!             Read over 1st 5 lines of description + 1st line of data
!             in future file that is repeated from siple record

    READ (lun(2),'(/////)')
    GO TO 210
  ELSEIF (ifuture .EQ. 8) THEN
    ireadf = 2

!             Read over 1 Header line in future FILE cis92a.dat
!             Note that 
!               - 1st line is 1990 (not identical to historical run)
!               - records are yearly, not every 6 months as for 
!                 1 <= ifuture <=7

    READ(lun(3),'(1x)')

!           READ first line of DATA (same year as last line in splco2.dat)
!           THEN for CIS92A only, replace that with 

    WRITE (stdout(),*) 'Replace: nsipl = ', nsipl
    READ(lun(3),*)yrco2rec(nsipl), atmco2rec(nsipl), dummy
    write(stdout(),*)yrco2rec(nsipl), atmco2rec(nsipl), dummy
    GO TO 210
  ELSE
    GO TO 222
  ENDIF

! 221    nco2rec = nsipl + nstab
  221 continue
  nco2rec = irec
END DO

222  CONTINUE

!      WRITE(0,*) 'Number records in splco2.dat:', nsipl
!      WRITE(0,*) 'Number records in future (stab.dat or cis92a.dat):', nstab
!      WRITE(0,*) 'Sum                          ', nstab + nsipl

 WRITE (stdout(),*)                                             &
      'Atm. CO2:  No. of entries for 1-box atmosphere =', nco2rec

call mpp_close(lun(1))
call mpp_close(lun(2))
call mpp_close(lun(3))

RETURN
END subroutine  read_co2atm  !}
! </SUBROUTINE> NAME="read_co2atm"


!#######################################################################
! <SUBROUTINE NAME="set_array">
!
! <DESCRIPTION>
!       Set up an array covering the model domain with a user-specified
!       value, in user-specified regions. There are a given number of
!       2-d regions specified by the values slat, nlat, wlon and elon.
!       The longitudes are for a cyclic domain, and if wlon and elon
!       are on opposite sides of the cut, the correct thing will
!       be done. Elon is considered to be east of wlon, so if elon is
!       less than wlon, then the region east of elon to the cut will be
!       filled, and the region from the cut to wlon will be filled.
!
!       After setting up the array in this routine, it may prove useful
!       to allow fine-tuning the settings via an array in a namelist.
!
!       Arguments:
!         Input:
!      num_regions = number of user-specified regions which will be
!                    filled
!
!             wlon = 1-d array of western (starting) longitudes for the
!                    rectangular regions
!
!             elon = 1-d array of eastern (ending) longitudes for the
!                    rectangular regions
!
!             slat = 1-d array of southern (starting) latitudes for the
!                    rectangular regions
!
!             nlat = 1-d array of northern (ending) latitudes for the
!                    rectangular regions
!
!                       Note: if slat >= nlat, then nothing is done
!                             for that region
!
!        set_value = the value to assign to array in the user-specified
!                    regions
!
!      unset_value = the value to assign to array outside of the
!                    user-specified regions
!
!             name = character variable used in informative messages
!
!     coastal_only = true to limit changes only to coastal points
!                    (i.e., at least one bordering point is land)
!
!         Output:
!
!            array = 2-d array which will contain the set- and unset-
!                    values. The array is assumed to have a border
!                    one unit wide on all edges, ala MOM. A cyclic
!                    boundary condition will be set if requested.
! </DESCRIPTION>
!

subroutine set_array(array, isd, ied, jsd, jed,                 &
                     xt, yt, kmt,                               &
                     num_regions, wlon_in, elon_in, slat, nlat, &
                     set_value, unset_value, name,              &
                     coastal_only)  !{

!
!       local parameters
!

character(len=64), parameter    :: sub_name = 'set_array'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header =                                &
     '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

!
!       arguments
!

integer, intent(in)                             :: isd
integer, intent(in)                             :: ied
integer, intent(in)                             :: jsd
integer, intent(in)                             :: jed
integer, intent(in)                             :: num_regions
real, dimension(isd:ied,jsd:jed), intent(out)   :: array
logical, intent(in)                             :: coastal_only
real, dimension(num_regions), intent(in)        :: elon_in
integer, dimension(isd:ied,jsd:jed), intent(in) :: kmt
character*(*), intent(in)                       :: name
real, dimension(num_regions), intent(in)        :: nlat
real, intent(in)                                :: set_value
real, dimension(num_regions), intent(in)        :: slat
real, intent(in)                                :: unset_value
real, dimension(num_regions), intent(in)        :: wlon_in
real, dimension(isd:ied,jsd:jed), intent(in)    :: xt
real, dimension(isd:ied,jsd:jed), intent(in)    :: yt

!
!       local variables
!

integer :: i, j, n
real, dimension(:), allocatable :: wlon
real, dimension(:), allocatable :: elon

!
!       save the longitudes in case they need to be modified
!

allocate(wlon(num_regions))
allocate(elon(num_regions))

wlon(:) = wlon_in(:)
elon(:) = elon_in(:)

!
! loop over the regions, applying changes as necessary
!

do n = 1, num_regions  !{

  if (nlat(n) .ge. slat(n)) then  !{
    write (stdout(),*)
    write (stdout(),*) trim(note_header),                          &
                       trim(name), ' region: ', n

!
!       make sure that all longitudes are in the range [0,360]
!

    do while (wlon(n) .gt. 360.0)  !{
      wlon(n) = wlon(n) - 360.0
    enddo  !}
    do while (wlon(n) .lt. 0.0)  !{
      wlon(n) = wlon(n) + 360.0
    enddo  !}
    do while (elon(n) .gt. 360.0)  !{
      elon(n) = elon(n) - 360.0
    enddo  !}
    do while (elon(n) .lt. 0.0)  !{
      elon(n) = elon(n) + 360.0
    enddo  !}

!
!       if the southern and northern latitudes are the same, then
!       find the grid box which encompasses them ...
!

    if (slat(n) .eq. nlat(n)) then  !{

     call mpp_error(FATAL, trim(error_header) //                &
                    'Equal latitudes not supported')

    elseif (wlon(n) .eq. elon(n)) then  !}{

     call mpp_error(FATAL, trim(error_header) //                &
                    'Equal longitudes not supported')

    else  !}{

!
!       ... else find all boxes where the center lies in the
!       rectangular region
!

      do j = jsd, jed  !{
        do i = isd, ied  !{
          if (nlat(n) .ge. yt(i,j) .and.                        &
              slat(n) .le. yt(i,j) .and.                        &
              lon_between(xt(i,j), wlon(n), elon(n))) then  !{
            array(i,j) = set_value
          endif  !}
        enddo  !} i
      enddo  !} j

    endif  !}

  endif  !}

enddo  !} n

!
!       if desired only apply mask to coastal regions
!

if (coastal_only) then  !{
  do j = jsd, jed  !{
    do i = isd, ied  !{
      if (kmt(i,j) .ne. 0 .and.                         &
          array(i,j) .eq. set_value) then  !{

!
!       if all the surrounding points are ocean, then this is not
!       a coastal point, therefore reset the mask
!

        if (kmt(i-1,j) .ne. 0 .and.                     &
            kmt(i+1,j) .ne. 0 .and.                     &
            kmt(i,j-1) .ne. 0 .and.                     &
            kmt(i,j+1) .ne. 0) then  !{
          array(i,j) = unset_value
        endif  !}
      endif  !}
    enddo  !} i
  enddo  !} j
endif  !}

!
!       clean up
!

deallocate(wlon)
deallocate(elon)

return

contains

!
!       Return true if w <= x_in <= e, taking into account the
!       periodicity of longitude.
!
!       x_in    = value to test
!
!       w       = west longitude of boundary
!
!       e       = east longitude of boundary
!

function lon_between(x_in, w, e)  !{

!
!       function definition
!

logical :: lon_between

!
!       arguments
!

real, intent(in)                :: x_in
real, intent(in)                :: w
real, intent(in)                :: e

!
!       local variables
!

!real                   :: w
!real                   :: e
real                    :: x

!
!       Save input values so we may modify them safely
!

x = x_in

!
!       make sure that all longitudes are in the range [0,360]
!

do while (x .gt. 360.0)  !{
  x = x - 360.0
enddo  !}
do while (x .lt. 0.0)  !{
  x = x + 360.0
enddo  !}
 
if (w .gt. e) then  !{
  lon_between = w .le. x .or. x .le. e
else  !}{
  lon_between = w .le. x .and. x .le. e
endif  !}

return

end function  lon_between  !}

end subroutine  set_array  !}
! </SUBROUTINE> NAME="set_array"

end module  ocmip2_biotic_mod  !}
