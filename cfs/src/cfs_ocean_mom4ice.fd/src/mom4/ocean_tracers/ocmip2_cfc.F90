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
! Ocean Carbon Model Intercomparison Study II: CFC module
!</OVERVIEW>
!
!<DESCRIPTION>
!       Implementation of routines to solve the OCMIP-2 CFC
!       simulations as outlined in the CFC-HOWTO documentation,
!       revision 1.6, 1999/04/29.
!</DESCRIPTION>
!
! <INFO>
! <REFERENCE>
! http://www.ipsl.jussieu.fr/OCMIP/phase2/simulations/CFC/HOWTO-CFC.html
! </REFERENCE>
! </INFO>
!
! $Id$
!

module  ocmip2_cfc_mod  !{

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

use field_manager_mod
use ocean_tpm_util_mod
use mpp_mod, only : stdout, stdlog, mpp_error, FATAL

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

public  :: ocmip2_cfc_bbc
public  :: ocmip2_cfc_end
public  :: ocmip2_cfc_init
public  :: ocmip2_cfc_sbc
public  :: ocmip2_cfc_source
public  :: ocmip2_cfc_start
public  :: ocmip2_cfc_tracer

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
private :: read_cfc_timehist

!
!----------------------------------------------------------------------
!
!       Private parameters
!
!----------------------------------------------------------------------
!

character(len=fm_field_name_len), parameter     :: package_name = 'ocmip2_cfc'
character(len=48), parameter                    :: mod_name = 'ocmip2_cfc_mod'
character(len=fm_string_len), parameter         :: default_file_in = 'INPUT/ocmip2_cfc.res.nc'
character(len=fm_string_len), parameter         :: default_file_out = 'RESTART/ocmip2_cfc.res.nc'

integer, parameter :: max_cfc_rec = 1200

!
!----------------------------------------------------------------------
!
!       Private types
!
!----------------------------------------------------------------------
!
 
type cfc_type  !{

  real                                  :: a1_11
  real                                  :: a2_11
  real                                  :: a3_11
  real                                  :: a4_11
  real                                  :: d1_11
  real                                  :: d2_11
  real                                  :: d3_11
  real                                  :: d4_11
  real                                  :: e1_11
  real                                  :: e2_11
  real                                  :: e3_11
  real                                  :: a1_12
  real                                  :: a2_12
  real                                  :: a3_12
  real                                  :: a4_12
  real, _ALLOCATABLE, dimension(:,:)         :: csat_11 _NULL
  real, _ALLOCATABLE, dimension(:,:)         :: csat_12 _NULL
  real                                  :: d1_12
  real                                  :: d2_12
  real                                  :: d3_12
  real                                  :: d4_12
  real                                  :: e1_12
  real                                  :: e2_12
  real                                  :: e3_12
  character(len=fm_string_len)          :: file_in
  character(len=fm_string_len)          :: file_out
  integer                               :: id_csat_11 = -1
  integer                               :: id_sc_11 = -1
  integer                               :: id_kw_11 = -1
  integer                               :: id_alpha_11 = -1
  integer                               :: id_csat_12 = -1
  integer                               :: id_sc_12 = -1
  integer                               :: id_kw_12 = -1
  integer                               :: id_alpha_12 = -1
  integer                               :: ind_cfc_11
  integer                               :: ind_cfc_12
  character(len=fm_field_name_len)      :: name
  real, _ALLOCATABLE, dimension(:,:)         :: sc_11 _NULL
  real, _ALLOCATABLE, dimension(:,:)         :: sc_12 _NULL
  real, _ALLOCATABLE, dimension(:,:)         :: alpha_11 _NULL
  real, _ALLOCATABLE, dimension(:,:)         :: alpha_12 _NULL
  real, _ALLOCATABLE, dimension(:,:)         :: kw_11 _NULL
  real, _ALLOCATABLE, dimension(:,:)         :: kw_12 _NULL

end type cfc_type  !}

!
!----------------------------------------------------------------------
!
!       Public variables
!
!----------------------------------------------------------------------
!

logical, public :: do_ocmip2_cfc

!
!----------------------------------------------------------------------
!
!       Private variables
!
!----------------------------------------------------------------------
!

!
!----------------------------------------------------------------------
!
!     VARIABLES:
!        fice        : sea ice fraction (0..1) 
!        xkw         : defined as Xconv * 0.337 * (u^2 + v), where Xconv
!                      is a conversion factor, u is the instantaneous
!                      SSMI winds, and v is the variance of the winds. 
!                      units are [cm/s]
!        patm        : atmospheric pressure. units are [atm]
!
!        seaicefract_file  : bc file containing sea-ice fraction
!        pistonveloc_file  : bc file containing xkw 
!        atmpress_file     : bc file containing atm. pressure
!  
!----------------------------------------------------------------------
!

integer                                         :: atmpress_id
character(len=fm_string_len)                    :: atmpress_file    
character(len=fm_string_len)                    :: atmpress_name    
real, allocatable, dimension(:,:)               :: fice_t
real, allocatable, dimension(:,:)               :: patm_t
real, allocatable, dimension(:,:)               :: xkw_t
integer                                         :: pistonveloc_id
integer                                         :: seaicefract_id
character(len=fm_string_len)                    :: pistonveloc_file
character(len=fm_string_len)                    :: pistonveloc_name
character(len=fm_string_len)                    :: seaicefract_file
character(len=fm_string_len)                    :: seaicefract_name
type(cfc_type), allocatable, dimension(:)       :: cfc
integer                                         :: instances
integer                                         :: id_pcfc_11
integer                                         :: id_pcfc_12
real                                            :: pcfc_11_n
real                                            :: pcfc_11_s
real                                            :: pcfc_12_n
real                                            :: pcfc_12_s
real, allocatable, dimension(:,:)               :: pcfc_11
real, allocatable, dimension(:,:)               :: pcfc_12
real, allocatable, dimension(:,:)               :: interp_n
real, allocatable, dimension(:,:)               :: interp_s
integer                                         :: package_index
integer                                         :: num_cfc_rec
real, dimension(max_cfc_rec)                    :: year_cfc_rec
real, dimension(max_cfc_rec)                    :: pcfc_11_n_rec
real, dimension(max_cfc_rec)                    :: pcfc_11_s_rec
real, dimension(max_cfc_rec)                    :: pcfc_12_n_rec
real, dimension(max_cfc_rec)                    :: pcfc_12_s_rec
real                                            :: min_cfc_year
real                                            :: max_cfc_year
character(len=fm_string_len)                    :: atmcfc_file
real                                            :: pert_start_year
real                                            :: pert_start_year_model
real                                            :: pert_time

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

implicit none

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
allocate( pcfc_11(isc:iec,jsc:jec) )
allocate( pcfc_12(isc:iec,jsc:jec) )
allocate( interp_n(isc:iec,jsc:jec) )
allocate( interp_s(isc:iec,jsc:jec) )

!
!       initialize some arrays
!

xkw_t(:,:) = 0.0
patm_t(:,:) = 0.0
fice_t(:,:) = 0.0
pcfc_11(:,:) = 0.0
pcfc_12(:,:) = 0.0
interp_n(:,:) = 0.0
interp_s(:,:) = 0.0

!
!       allocate cfc array elements
!

do n = 1, instances  !{

  allocate( cfc(n)%sc_11(isc:iec,jsc:jec) )
  allocate( cfc(n)%kw_11(isc:iec,jsc:jec) )
  allocate( cfc(n)%alpha_11(isc:iec,jsc:jec) )
  allocate( cfc(n)%csat_11(isc:iec,jsc:jec) )
  allocate( cfc(n)%sc_12(isc:iec,jsc:jec) )
  allocate( cfc(n)%kw_12(isc:iec,jsc:jec) )
  allocate( cfc(n)%alpha_12(isc:iec,jsc:jec) )
  allocate( cfc(n)%csat_12(isc:iec,jsc:jec) )

enddo  !}

!
!       initialize some arrays
!

do n = 1, instances  !{

  cfc(n)%sc_11(:,:) = 0.0
  cfc(n)%kw_11(:,:) = 0.0
  cfc(n)%alpha_11(:,:) = 0.0
  cfc(n)%csat_11(:,:) = 0.0
  cfc(n)%sc_12(:,:) = 0.0
  cfc(n)%kw_12(:,:) = 0.0
  cfc(n)%alpha_12(:,:) = 0.0
  cfc(n)%csat_12(:,:) = 0.0

enddo  !} n



return
end subroutine  allocate_arrays  !}
! </SUBROUTINE> NAME="allocate_arrays"


!#######################################################################
! <SUBROUTINE NAME="bc_interp">
!
! <DESCRIPTION>
!     Interpolates atmospheric CFC-11 and CFC-12 to the timestep of the model
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
!     OUTPUT: 
!
!     [rs] - cfc11   =  Atmospheric CFC-11 (ppvt) at "year_model"
!
!     [ra] - cfc12   =  3-member array of atmospheric C-14 at
!                         year_model".  Sequentially, the 3 values
!                         correspond to forcing in 3 latitudinal bands:
!                         (1) 90S - 20S,
!                         (2) 20S - 20N, and
!                         (3) 20N - 90N.
!
!     James Orr, LSCE/CEA-CNRS, Saclay, France, 20 April 1999
! </DESCRIPTION>
!

subroutine bc_interp(year_model, pcfc_11_n, pcfc_11_s, pcfc_12_n, pcfc_12_s)  !{

implicit none 

!
!----------------------------------------------------------------
!       Arguments
!----------------------------------------------------------------
!

real, intent(in)        :: year_model
real, intent(out)       :: pcfc_11_n
real, intent(out)       :: pcfc_11_s
real, intent(out)       :: pcfc_12_n
real, intent(out)       :: pcfc_12_s

!
!----------------------------------------------------------------
!       Local variables
!----------------------------------------------------------------
!

real    :: year

integer :: iz
integer :: n
real    :: x

logical, save   :: initialized = .false.

if (.not. initialized) then
  call read_cfc_timehist
  initialized = .true.
endif

!     ------------------------------------------------------------------
!     Interpolate atmospheric CFCs to timestep (decimal years)
!     of the model
!     ------------------------------------------------------------------

if (year_model .lt. min_cfc_year) then 
  year = min_cfc_year
else if (year_model .ge. min_cfc_year .and. year_model .le. max_cfc_year) then 
  year = year_model
else if (year_model .gt. max_cfc_year) then 
  year = max_cfc_year
end if 

!     Find relative POSITION n for year_model in CFC record

call locate(year_cfc_rec, num_cfc_rec, year, n)

!     Determine linear interpolation factor "x"

x = (year - year_cfc_rec(n)) / (year_cfc_rec(n+1) - year_cfc_rec(n))

!     Perform temporal interpolation for atmospheric CO2

pcfc_11_n = pcfc_11_n_rec(n) * (1.0 - x) + pcfc_11_n_rec(n+1) * x
pcfc_11_s = pcfc_11_s_rec(n) * (1.0 - x) + pcfc_11_s_rec(n+1) * x
pcfc_12_n = pcfc_12_n_rec(n) * (1.0 - x) + pcfc_12_n_rec(n+1) * x
pcfc_12_s = pcfc_12_s_rec(n) * (1.0 - x) + pcfc_12_s_rec(n+1) * x

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

implicit none

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
! <SUBROUTINE NAME="ocmip2_cfc_bbc">
!
! <DESCRIPTION>
!     calculate the surface boundary conditions
! </DESCRIPTION>
!

subroutine ocmip2_cfc_bbc  !{

!
!-----------------------------------------------------------------------
!     modules (have to come first)
!-----------------------------------------------------------------------
!

implicit none

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

character(len=64), parameter    :: sub_name = 'ocmip2_cfc_bbc'
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

!
! =====================================================================
!     begin executable code
! =====================================================================
!

!
!   no flux bottom boundary condition is the default
!

return

end subroutine  ocmip2_cfc_bbc  !}
! </SUBROUTINE> NAME="ocmip2_cfc_bbc"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_cfc_end">
!
! <DESCRIPTION>
!     Clean up various CFC quantities for this run.
! </DESCRIPTION>
!

subroutine ocmip2_cfc_end  !{ 

!
!-----------------------------------------------------------------------
!     modules (have to come first)
!-----------------------------------------------------------------------
!

use fms_mod, only : write_data

implicit none

!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'ocmip2_cfc_end'
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

integer :: n

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
!-----------------------------------------------------------------------
!     save out the perturbation time
!-----------------------------------------------------------------------
!

write (stdout(),*)
write (stdout(),*) trim(note_header),                         &
     ' Saving additional output to a restart file...'

do n = 1, instances  !{

  call write_data(cfc(n)%file_out, 'pert_start_year', pert_start_year)
  call write_data(cfc(n)%file_out, 'pert_start_year_model', pert_start_year_model)
  call write_data(cfc(n)%file_out, 'pert_time', pert_time)

enddo  !} n

return
end subroutine  ocmip2_cfc_end  !}
! </SUBROUTINE> NAME="ocmip2_cfc_end"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_cfc_sbc">
!
! <DESCRIPTION>
!     Calculate the surface boundary conditions
! </DESCRIPTION>
!

subroutine ocmip2_cfc_sbc  !{

!
!-----------------------------------------------------------------------
!     modules (have to come first)
!-----------------------------------------------------------------------
!

use diag_manager_mod, only: send_data
use time_interp_external_mod, only: time_interp_external
use time_manager_mod, only: days_in_year, days_in_month,        &
     get_date, set_date

implicit none

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

character(len=64), parameter    :: sub_name = 'ocmip2_cfc_sbc'
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
integer :: n
integer :: m
integer :: hour
integer :: day
real    :: days_in_this_year
integer :: minute
integer :: month
integer :: num_days
integer :: second
integer :: year
real    :: ta
real    :: sal
logical :: used

!
! =====================================================================
!     begin executable code
! =====================================================================
!

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
!-----------------------------------------------------------------------
!     Update pert_time
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

!
!-----------------------------------------------------------------------
!     assign perturbation time (real years e.g. 1935.67)
!-----------------------------------------------------------------------
!	

pert_time = year + days_in_this_year/num_days -                 &
     pert_start_year_model + pert_start_year

!
!---------------------------------------------------------------------
!     calculate atmospheric pressures by calling routine bc_interp
!---------------------------------------------------------------------
!

call bc_interp(pert_time, pcfc_11_n, pcfc_11_s, pcfc_12_n, pcfc_12_s)

!
!-----------------------------------------------------------------------
!     "interpolate" in space
!-----------------------------------------------------------------------
!

do j = jsc, jec  !{
  do i = isc, iec  !{
    pcfc_11(i,j) = interp_n(i,j) * pcfc_11_n + interp_s(i,j) * pcfc_11_s
    pcfc_12(i,j) = interp_n(i,j) * pcfc_12_n + interp_s(i,j) * pcfc_12_s
  enddo  !} i
enddo  !} j

do n = 1, instances  !{

!
!---------------------------------------------------------------------
!     Calculate solubilities
!       Use Warner and Weiss (1985) DSR, vol 32, final result
!       in mol/cm3/pptv (1 part per trillion 1e-12)
!       use Bullister and Wisegavger for CCl4
!
!     the factor 1.e-09 is for the conversion from mol/(l * atm) 
!        to mol/(m3 * pptv) 
!---------------------------------------------------------------------
!

  do j = jsc, jec  !{
    do i = isc, iec  !{
      ta = (t_prog(indtemp)%field(i,j,1,taum1) + 273.15) * 0.01
      sal = t_prog(indsal)%field(i,j,1,taum1)

      cfc(n)%alpha_11(i,j) =                                                            &
           exp(cfc(n)%d1_11 + cfc(n)%d2_11 / ta + cfc(n)%d3_11 * log(ta) +              &
               cfc(n)%d4_11* ta * ta +                                                  &
               sal * ((cfc(n)%e3_11 * ta + cfc(n)%e2_11) * ta + cfc(n)%e1_11)) *        &
           1.0e-09 * grid%tmask(i,j,1)

      cfc(n)%alpha_12(i,j) =                                                            &
           exp(cfc(n)%d1_12 + cfc(n)%d2_12 / ta + cfc(n)%d3_12 * log(ta) +              &
               cfc(n)%d4_12* ta * ta +                                                  &
               sal * ((cfc(n)%e3_12 * ta + cfc(n)%e2_12) * ta + cfc(n)%e1_12)) *        &
           1.0e-09 * grid%tmask(i,j,1)
    enddo  !} i
  enddo  !} j

!
!---------------------------------------------------------------------
!     Calculate Schmidt numbers
!      use coefficients given by Zheng et al (1998), JGR vol 103, C1
!---------------------------------------------------------------------
!

  do j = jsc, jec  !{
    do i = isc, iec  !{
      cfc(n)%sc_11(i,j) = cfc(n)%a1_11 + t_prog(indtemp)%field(i,j,1,taum1) *           &
           (cfc(n)%a2_11 + t_prog(indtemp)%field(i,j,1,taum1) *                         &
            (cfc(n)%a3_11 + t_prog(indtemp)%field(i,j,1,taum1) * cfc(n)%a4_11)) *       &
           grid%tmask(i,j,1)
      cfc(n)%sc_12(i,j) = cfc(n)%a1_12 + t_prog(indtemp)%field(i,j,1,taum1) *           &
           (cfc(n)%a2_12 + t_prog(indtemp)%field(i,j,1,taum1) *                         &
            (cfc(n)%a3_12 + t_prog(indtemp)%field(i,j,1,taum1) * cfc(n)%a4_12)) *       &
           grid%tmask(i,j,1)
    enddo  !} i
  enddo  !} j 

!
!---------------------------------------------------------------------
!     calculate oceanic equilibrium concentrations
!        pcfc_11/12 are given in picoatm
!        alpha in mol/m3/ppvt
!        this gives mol/m3
!---------------------------------------------------------------------
!

  do j = jsc, jec  !{
    do i = isc, iec  !{
      cfc(n)%csat_11(i,j) = cfc(n)%alpha_11(i,j) *                      &
           pcfc_11(i,j) * patm_t(i,j) * grid%tmask(i,j,1)
      cfc(n)%csat_12(i,j) = cfc(n)%alpha_12(i,j) *                      &
           pcfc_12(i,j) * patm_t(i,j) * grid%tmask(i,j,1)
    enddo  !} i
  enddo  !} j

!
!---------------------------------------------------------------------
!     calculate piston-velocities
!      including  effect of sea-ice
!      xkw is given in m/s (converted in read_cfc_bc), therefore
!      kw is also in m/s
!---------------------------------------------------------------------
!

  do j = jsc, jec  !{
    do i = isc, iec  !{
      cfc(n)%kw_11(i,j) = (1.0 - fice_t(i,j)) * xkw_t(i,j) *             &
           sqrt(660.0 / cfc(n)%sc_11(i,j)) * grid%tmask(i,j,1)
      cfc(n)%kw_12(i,j) = (1.0 - fice_t(i,j)) * xkw_t(i,j) *             &
           sqrt(660.0 / cfc(n)%sc_12(i,j)) * grid%tmask(i,j,1)
    enddo  !} i
  enddo  !} j 

!
!---------------------------------------------------------------------
!     calculate surface fluxes for CFCs
!       kw is in m/s, csat and t are in mol/m3, therefore
!       stf is in mol/m2/s
!---------------------------------------------------------------------
!

  do j = jsc, jec  !{
    do i = isc, iec  !{
      t_prog(cfc(n)%ind_cfc_11)%stf(i,j) = cfc(n)%kw_11(i,j) *  &
            (cfc(n)%csat_11(i,j) -                              &
             t_prog(cfc(n)%ind_cfc_11)%field(i,j,1,taum1))
      t_prog(cfc(n)%ind_cfc_12)%stf(i,j) = cfc(n)%kw_12(i,j) *  &
            (cfc(n)%csat_12(i,j) -                              &
             t_prog(cfc(n)%ind_cfc_12)%field(i,j,1,taum1))
    enddo  !} i
  enddo  !} j 

enddo  !} n 

!
!-----------------------------------------------------------------------
!       Save variables for diagnostics
!-----------------------------------------------------------------------
!

if (id_pcfc_11 .gt. 0) then
  used = send_data(id_pcfc_11, pcfc_11(isc:iec,jsc:jec),              &
       time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
endif
if (id_pcfc_12 .gt. 0) then
  used = send_data(id_pcfc_12, pcfc_12(isc:iec,jsc:jec),              &
       time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
endif

do n = 1, instances  !{

  if (cfc(n)%id_alpha_11 .gt. 0) then
    used = send_data(cfc(n)%id_sc_11,                                   &
         cfc(n)%alpha_11(isc:iec,jsc:jec),                              &
         time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
  endif
  if (cfc(n)%id_sc_11 .gt. 0) then
    used = send_data(cfc(n)%id_sc_11,                                   &
         cfc(n)%sc_11(isc:iec,jsc:jec),                                 &
         time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
  endif
  if (cfc(n)%id_kw_11 .gt. 0) then
    used = send_data(cfc(n)%id_kw_11,                                   &
         cfc(n)%kw_11(isc:iec,jsc:jec),                                 &
         time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
  endif
  if (cfc(n)%id_csat_11 .gt. 0) then
    used = send_data(cfc(n)%id_csat_11,                                 &
         cfc(n)%csat_11(isc:iec,jsc:jec),                               &
         time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
  endif
  if (cfc(n)%id_alpha_12 .gt. 0) then
    used = send_data(cfc(n)%id_sc_12,                                   &
         cfc(n)%alpha_12(isc:iec,jsc:jec),                              &
         time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
  endif
  if (cfc(n)%id_sc_12 .gt. 0) then
    used = send_data(cfc(n)%id_sc_12,                                   &
         cfc(n)%sc_12(isc:iec,jsc:jec),                                 &
         time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
  endif
  if (cfc(n)%id_kw_12 .gt. 0) then
    used = send_data(cfc(n)%id_kw_12,                                   &
         cfc(n)%kw_12(isc:iec,jsc:jec),                                 &
         time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
  endif
  if (cfc(n)%id_csat_12 .gt. 0) then
    used = send_data(cfc(n)%id_csat_12,                                 &
         cfc(n)%csat_12(isc:iec,jsc:jec),                               &
         time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,1))
  endif

enddo  !} n

return

end subroutine  ocmip2_cfc_sbc  !}
! </SUBROUTINE> NAME="ocmip2_cfc_sbc"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_cfc_init">
!
! <DESCRIPTION>
!       Set up any extra fields needed by the tracer packages
!
!       Save pointers to various "types", such as Grid and Domains.
! </DESCRIPTION>

subroutine ocmip2_cfc_init  !{

implicit none

!
!-----------------------------------------------------------------------
!       local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'ocmip2_cfc_init'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header =                                &
     '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

!
!-----------------------------------------------------------------------
!     Schmidt number coefficients 
!      Use coefficients given by Zheng et al (1998), JGR vol 103, C1
!         for CFC11 and CFC12
!-----------------------------------------------------------------------
!

real, parameter :: a1_11_def = 3501.8
real, parameter :: a2_11_def = -210.31
real, parameter :: a3_11_def =    6.1851
real, parameter :: a4_11_def =   -0.07513

real, parameter :: a1_12_def = 3845.4
real, parameter :: a2_12_def = -228.95
real, parameter :: a3_12_def =    6.1908
real, parameter :: a4_12_def =   -0.067430

!
!-----------------------------------------------------------------------
!     Solubility coefficients for alpha in mol/l/atm
!      (1) for CFC11, (2) for CFC12
!     after Warner and Weiss (1985) DSR, vol 32 for CFC11 and CFC12
!-----------------------------------------------------------------------
!

real, parameter :: d1_11_def = -229.9261
real, parameter :: d2_11_def =  319.6552
real, parameter :: d3_11_def =  119.4471
real, parameter :: d4_11_def =   -1.39165
real, parameter :: e1_11_def =   -0.142382
real, parameter :: e2_11_def =    0.091459
real, parameter :: e3_11_def =   -0.0157274
 
real, parameter :: d1_12_def = -218.0971
real, parameter :: d2_12_def =  298.9702
real, parameter :: d3_12_def =  113.8049
real, parameter :: d4_12_def =   -1.39165
real, parameter :: e1_12_def =   -0.143566
real, parameter :: e2_12_def =    0.091015
real, parameter :: e3_12_def =   -0.0153924

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
character(len=256)                                      :: caller_str
character(len=fm_string_len), pointer, dimension(:)     :: good_list

!
!-----------------------------------------------------------------------
!       Check which tracer packages have been turned on
!-----------------------------------------------------------------------
!

!
!       Initialize the ocmip2 cfc package
!

package_index = otpm_set_tracer_package(package_name,           &
     caller=trim(mod_name) // '(' // trim(sub_name) // ')',     &
     file_in=default_file_in, file_out=default_file_out)

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

if (instances .eq. 0) then  !{
  write (stdout(),*) trim(note_header), ' No instances'
  do_ocmip2_cfc = .false.
else  !}{
  if (instances .eq. 1) then  !{
    write (stdout(),*) trim(note_header), ' ', instances, ' instance'
  else  !}{
    write (stdout(),*) trim(note_header), ' ', instances, ' instances'
  endif  !}
  do_ocmip2_cfc = .true.
endif  !}

!
!       Return if we don't want to use this package,
!       after changing the list back
!

if (.not. do_ocmip2_cfc) then  !{
  return
endif  !}

! after reading tracer tree
!       allocate storage for cfc array
!

allocate ( cfc(instances) )

!
!       loop over the names, saving them into the cfc array
!

do n = 1, instances  !{

  if (fm_get_value(path_to_names, name, index = n)) then  !{
    cfc(n)%name = name
  else  !}{
    write (name,*) n
    call mpp_error(FATAL, trim(error_header) //                 &
         ' Bad field name for index ' // trim(name))
  endif  !}

enddo  !}

!
!       Set up the field input
!

do n = 1, instances  !{

  name = cfc(n)%name
  if (name(1:1) .eq. '_') then  !{
    suffix = ' '
    long_suffix = ' '
  else  !}{
    suffix = '_' // name
    long_suffix = ' (' // trim(name) // ')'
  endif  !}

!
!       CFC-11
!

  cfc(n)%ind_cfc_11 = otpm_set_prog_tracer('cfc_11' // suffix, package_name,    &
       longname = 'CFC-11' // trim(long_suffix),                                &
       units = 'mol/m^3', flux_units = 'mol/m^2/s',                             &
       caller=trim(mod_name) // '(' // trim(sub_name) // ')')

!
!       CFC-12
!

  cfc(n)%ind_cfc_12 = otpm_set_prog_tracer('cfc_12' // suffix, package_name,    &
       longname = 'CFC-12' // trim(long_suffix),                                &
       units = 'mol/m^3', flux_units = 'mol/m^2/s',                             &
       caller=trim(mod_name) // '(' // trim(sub_name) // ')')

enddo  !} n

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
!       Set up the *global* CFC namelist
!-----------------------------------------------------------------------
!

caller_str=trim(mod_name) // '(' // trim(sub_name) // ')'

call otpm_start_namelist(package_name, '*global*', caller = caller_str, no_overwrite = .true., &
     check = .true.)

call otpm_set_value('atmcfc_file', 'INPUT/cfc11_cfc12_atm.dat')
call otpm_set_value('atmpress_file', 'INPUT/atmpress_ocmip2.nc')
call otpm_set_value('atmpress_name', 'atmpress')
call otpm_set_value('pistonveloc_file', 'INPUT/pistonveloc_ocmip2.nc')
call otpm_set_value('pistonveloc_name', 'pistonveloc')
call otpm_set_value('seaicefract_file', 'INPUT/f_ice_ocmip2.nc')
call otpm_set_value('seaicefract_name', 'f_ice')
call otpm_set_value('pert_first', .false.)
call otpm_set_value('pert_start_year', 1931)
call otpm_set_value('pert_start_month', 1)
call otpm_set_value('pert_start_day', 1)
call otpm_set_value('pert_start_hour', 0)
call otpm_set_value('pert_start_minute', 0)
call otpm_set_value('pert_start_second', 0)

call otpm_end_namelist(package_name, '*global*', caller = caller_str, check = .true.)

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

!
!-----------------------------------------------------------------------
!       Set up the instance CFC namelists
!-----------------------------------------------------------------------
!

do n = 1, instances  !{

  call otpm_start_namelist(package_name, cfc(n)%name, caller = caller_str, no_overwrite = .true., &
       check = .true.)

  call otpm_set_value('a1_11', a1_11_def)
  call otpm_set_value('a2_11', a2_11_def)
  call otpm_set_value('a3_11', a3_11_def)
  call otpm_set_value('a4_11', a4_11_def)

  call otpm_set_value('a1_12', a1_12_def)
  call otpm_set_value('a2_12', a2_12_def)
  call otpm_set_value('a3_12', a3_12_def)
  call otpm_set_value('a4_12', a4_12_def)

  call otpm_set_value('d1_11', d1_11_def)
  call otpm_set_value('d2_11', d2_11_def)
  call otpm_set_value('d3_11', d3_11_def)
  call otpm_set_value('d4_11', d4_11_def)

  call otpm_set_value('d1_12', d1_12_def)
  call otpm_set_value('d2_12', d2_12_def)
  call otpm_set_value('d3_12', d3_12_def)
  call otpm_set_value('d4_12', d4_12_def)

  call otpm_set_value('e1_11', e1_11_def)
  call otpm_set_value('e2_11', e2_11_def)
  call otpm_set_value('e3_11', e3_11_def)

  call otpm_set_value('e1_12', e1_12_def)
  call otpm_set_value('e2_12', e2_12_def)
  call otpm_set_value('e3_12', e3_12_def)

  call otpm_set_value('file_in', default_file_in)
  call otpm_set_value('file_out', default_file_out)

  call otpm_end_namelist(package_name, cfc(n)%name, check = .true., caller = caller_str)

enddo  !} n

return

end subroutine ocmip2_cfc_init  !}
! </SUBROUTINE> NAME="ocmip2_cfc_init"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_cfc_source">
!
! <DESCRIPTION>
!     compute the source terms for the CFCs, including boundary
!     conditions (not done in setvbc, to minimize number
!     of hooks required in MOM base code)
! </DESCRIPTION>
!

subroutine ocmip2_cfc_source  !{

!
!-----------------------------------------------------------------------
!     modules (have to come first)
!-----------------------------------------------------------------------
!

use diag_manager_mod, only: send_data
use time_manager_mod, only: get_date
use time_interp_external_mod, only: time_interp_external

implicit none

!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'ocmip2_cfc_source'
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

!integer :: day
!integer :: month
!integer :: year
!integer :: hour
!integer :: minute
!integer :: second

!
! =====================================================================
!     begin executable code
! =====================================================================
!

!
!       get the model month
!

!call get_date(time%model_time, year, month, day,                &
              !hour, minute, second)

!
!-----------------------------------------------------------------------
!     calculate the source terms for CFCs
!-----------------------------------------------------------------------
!

return

end subroutine  ocmip2_cfc_source  !}
! </SUBROUTINE> NAME="ocmip2_cfc_source"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_cfc_start">
!
! <DESCRIPTION>
! Initialize variables, read in namelists, calculate constants
! for a given run and allocate diagnostic arrays
! </DESCRIPTION>
!

subroutine ocmip2_cfc_start  !{

!
!-----------------------------------------------------------------------
!       modules (have to come first)
!-----------------------------------------------------------------------
!

use time_manager_mod, only: days_in_year, days_in_month,        &
     get_date, set_date
use time_interp_external_mod, only: init_external_field
use diag_manager_mod, only: register_diag_field, diag_axis_init
use fms_mod, only : read_data, open_namelist_file,              &
                    close_file

implicit none

!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'ocmip2_cfc_start'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header =                                &
     '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
real, parameter    :: ys = -10.0
real, parameter    :: yn =  10.0

!
!-----------------------------------------------------------------------
!       local variables
!-----------------------------------------------------------------------
!

real                                    :: check
integer                                 :: day
real                                    :: days_in_this_year
integer                                 :: hour
character(len=fm_field_name_len+3)      :: long_suffix
integer                                 :: i
integer                                 :: j
integer                                 :: m
integer                                 :: n
integer                                 :: minute
integer                                 :: month
integer                                 :: num_days
integer                                 :: second
character(len=fm_field_name_len+1)      :: suffix
integer                                 :: year
logical                                 :: pert_first
real                                    :: pert_start_day
real                                    :: pert_start_hour
real                                    :: pert_start_minute
real                                    :: pert_start_month
real                                    :: pert_start_second
character(len=256)                      :: caller_str
real                                    :: test_pert_time
integer                                 :: test_pert_start_year
integer                                 :: test_pert_start_year_model

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
                  ' Starting ', trim(package_name), ' module'

!
!-----------------------------------------------------------------------
!     dynamically allocate the global CFC arrays
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

atmcfc_file       = otpm_get_string ('atmcfc_file', scalar = .true.)
atmpress_file     = otpm_get_string ('atmpress_file', scalar = .true.)
atmpress_name     = otpm_get_string ('atmpress_name', scalar = .true.)
pistonveloc_file  = otpm_get_string ('pistonveloc_file', scalar = .true.)
pistonveloc_name  = otpm_get_string ('pistonveloc_name', scalar = .true.)
seaicefract_file  = otpm_get_string ('seaicefract_file', scalar = .true.)
seaicefract_name  = otpm_get_string ('seaicefract_name', scalar = .true.)
pert_first        = otpm_get_logical('pert_first', scalar = .true.)
pert_start_year   = otpm_get_integer('pert_start_year', scalar = .true.)
pert_start_month  = otpm_get_integer('pert_start_month', scalar = .true.)
pert_start_day    = otpm_get_integer('pert_start_day', scalar = .true.)
pert_start_hour   = otpm_get_integer('pert_start_hour', scalar = .true.)
pert_start_minute = otpm_get_integer('pert_start_minute', scalar = .true.)
pert_start_second = otpm_get_integer('pert_start_second', scalar = .true.)

call otpm_end_namelist(package_name, '*global*', caller = caller_str)

do n = 1, instances  !{

  call otpm_start_namelist(package_name, cfc(n)%name, caller = caller_str)

  cfc(n)%a1_11 = otpm_get_real('a1_11', scalar = .true.)
  cfc(n)%a2_11 = otpm_get_real('a2_11', scalar = .true.)
  cfc(n)%a3_11 = otpm_get_real('a3_11', scalar = .true.)
  cfc(n)%a4_11 = otpm_get_real('a4_11', scalar = .true.)
  cfc(n)%a1_12 = otpm_get_real('a1_12', scalar = .true.)
  cfc(n)%a2_12 = otpm_get_real('a2_12', scalar = .true.)
  cfc(n)%a3_12 = otpm_get_real('a3_12', scalar = .true.)
  cfc(n)%a4_12 = otpm_get_real('a4_12', scalar = .true.)

  cfc(n)%d1_11 = otpm_get_real('d1_11', scalar = .true.)
  cfc(n)%d2_11 = otpm_get_real('d2_11', scalar = .true.)
  cfc(n)%d3_11 = otpm_get_real('d3_11', scalar = .true.)
  cfc(n)%d4_11 = otpm_get_real('d4_11', scalar = .true.)
  cfc(n)%d1_12 = otpm_get_real('d1_12', scalar = .true.)
  cfc(n)%d2_12 = otpm_get_real('d2_12', scalar = .true.)
  cfc(n)%d3_12 = otpm_get_real('d3_12', scalar = .true.)
  cfc(n)%d4_12 = otpm_get_real('d4_12', scalar = .true.)

  cfc(n)%e1_11 = otpm_get_real('e1_11', scalar = .true.)
  cfc(n)%e2_11 = otpm_get_real('e2_11', scalar = .true.)
  cfc(n)%e3_11 = otpm_get_real('e3_11', scalar = .true.)
  cfc(n)%e1_12 = otpm_get_real('e1_12', scalar = .true.)
  cfc(n)%e2_12 = otpm_get_real('e2_12', scalar = .true.)
  cfc(n)%e3_12 = otpm_get_real('e3_12', scalar = .true.)

  cfc(n)%file_in  = otpm_get_string('file_in', scalar = .true.)
  cfc(n)%file_out = otpm_get_string('file_out', scalar = .true.)

  call otpm_end_namelist(package_name, cfc(n)%name, caller = caller_str)

enddo  !} n

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
       ' Could not open atmpress file: ' //                      &
       trim(atmpress_file))
endif  !}

pistonveloc_id = init_external_field(pistonveloc_file,          &
                                     pistonveloc_name,          &
                                     domain = Domain%domain2d)
if (pistonveloc_id .eq. 0) then  !{
  call mpp_error(FATAL, trim(error_header) //                   &
       ' Could not open pistonveloc file: ' //                   &
       trim(pistonveloc_file))
endif  !}

seaicefract_id = init_external_field(seaicefract_file,          &
                                     seaicefract_name,          &
                                     domain = Domain%domain2d)
if (seaicefract_id .eq. 0) then  !{
  call mpp_error(FATAL, trim(error_header) //                   &
       ' Could not open seaicefract file: ' //                   &
       trim(seaicefract_file))
endif  !}

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
           ' Perturbation start time is at different' //       &
           ' point in year for instance ' // trim(cfc(n)%name))
    endif  !}

    pert_start_year_model = year

    write (stdout(),*)
    write (stdout(),*) trim(note_header),                     &
                       ' Instance ', trim(cfc(n)%name)
    write (stdout(),*) 'This is the beginning of a perturbation run'
    write (stdout(),*) '   perturbation start year: ',        &
                       pert_start_year
    write (stdout(),*) '   perturbation start year (model): ',&
                       pert_start_year_model
    write (stdout(),*) '   itt : ', time%itt

!
!-----------------------------------------------------------------------
!     otherwise read info and cumulative input from a file
!       check also consistency of input 
!-----------------------------------------------------------------------
!

  else  !}{

    write(stdout(),*)
    write(stdout(),*) trim(note_header),                        &
         ' Getting pert values from restart file',               &
         ' for instance ', trim(cfc(n)%name)

!
!       save variables to check that they are the same for all files
!

    if (n .gt. 1) then  !{
      test_pert_start_year = pert_start_year
      test_pert_start_year_model = pert_start_year_model
      test_pert_time = pert_time
    endif  !}

    call read_data(cfc(n)%file_in, 'pert_start_year',                   &
         pert_start_year, timelevel=1)
    call read_data(cfc(n)%file_in, 'pert_start_year_model',             &
         pert_start_year_model, timelevel=1)
    call read_data(cfc(n)%file_in, 'pert_time',                         &
         pert_time, timelevel=1)

!
!       check that variables are the same for all files
!

    if (n .gt. 1) then  !{
      if (pert_time .ne. test_pert_time) then  !{
        call mpp_error(FATAL, trim(error_header) //                     &
             ' pert_time does not match for instances ' //               &
             trim(cfc(n-1)%name) // ' and ' // trim(cfc(n)%name))
      endif  !}
      if (pert_start_year .ne. test_pert_start_year) then  !{
        call mpp_error(FATAL, trim(error_header) //                     &
             ' pert_start_year does not match for instances ' //         &
             trim(cfc(n-1)%name) // ' and ' // trim(cfc(n)%name))
      endif  !}
      if (pert_start_year_model .ne. test_pert_start_year_model) then  !{
        call mpp_error(FATAL, trim(error_header) //                     &
             ' pert_start_year_model does not match for instances ' //   &
             trim(cfc(n-1)%name) // ' and ' // trim(cfc(n)%name))
      endif  !}
    endif  !}
        
    check = year - pert_start_year_model +                      &
            pert_start_year +                                   &
            days_in_this_year / num_days

    if (check .ne. pert_time) then  !{

      write (stdout(),*) 
      write (stdout(),*) trim(error_header)
      write (stdout(),*) 'Model time does not match ',          &
                         ' perturbation time!'
      write (stdout(),*) '     pert_start_year            = ',  &
                         pert_start_year
      write (stdout(),*) '     pert_start_year_model      = ',  &
                         pert_start_year_model
      write (stdout(),*) '     pert_time from restart    = ',   &
                         pert_time
      write (stdout(),*) '     pert_time calc from model = ',   &
                         check
      call mpp_error(FATAL, trim(error_header) // ' Time mismatch')

    endif  !}

    write (stdout(),*) 
    write (stdout(),*) trim(note_header)
    write (stdout(),*) '     pert_start_year       = ',         &
                       pert_start_year
    write (stdout(),*) '     pert_start_year_model = ',         &
                       pert_start_year_model
    write (stdout(),*) '     pert_time             = ',         &
                       pert_time
    write (stdout(),*) 

!
!-----------------------------------------------------------------------
!     end of check for pert_first
!-----------------------------------------------------------------------
!

  endif  !}

enddo  !} n

!
!-----------------------------------------------------------------------
!     determine the interpolation coefficient for the spatial
!        interpolation of the atmospheric CFC field 
!     Scheme: Two stations (41S and 45N), each representative of 
!               its own hemisphere, except between 10S and 10N WHERE 
!               values are interpolated linearly. 
!               Thus there are 3 regions:
!
!               (1) 90s to 10s, WHERE values take on the value at the 
!                   station at 41s;
!               (2) 10n to 90n, WHERE values take on the value at the 
!                   station at 45n; and
!               (3) 10s to 10n, WHERE values are interpolated
!
!-----------------------------------------------------------------------
!

do j = jsc, jec  !{
  do i = isc, iec  !{
    if (grid%yt(i,j) .ge. yn) then
      interp_n(i,j) = 1.0
      interp_s(i,j) = 0.0
    else if (grid%yt(i,j) .le. ys) then
      interp_n(i,j) = 0.0
      interp_s(i,j) = 1.0
    else 
      interp_n(i,j) = (grid%yt(i,j) - ys) / (yn - ys)
      interp_s(i,j) = 1.0 - interp_n(i,j)
    endif 
  enddo  !} i
enddo  !} j

!
!-----------------------------------------------------------------------
!     Set up analyses
!-----------------------------------------------------------------------
!

!
!       register the fields
!

id_pcfc_11 = register_diag_field('ocean_model',                 &
     'pcfc_11', grid%tracer_axes(1:2),                          &
     Time%model_time, 'Atmospheric pCFC-11', 'pptv',            &
     missing_value = -1.0e+10)

id_pcfc_12 = register_diag_field('ocean_model',                 &
     'pcfc_12', grid%tracer_axes(1:2),                          &
     Time%model_time, 'Atmospheric pCFC-12', 'pptv',            &
     missing_value = -1.0e+10)

do n = 1, instances  !{

  if (cfc(n)%name(1:1) .eq. '_') then  !{
    suffix = ' '
    long_suffix = ' '
  else  !}{
    suffix = '_' // cfc(n)%name
    long_suffix = ' (' // trim(cfc(n)%name) // ')'
  endif  !}

  cfc(n)%id_sc_11 = register_diag_field('ocean_model',                  &
       'sc_11'//trim(suffix), grid%tracer_axes(1:2),                    &
       Time%model_time,                                                 &
       'Schmidt number - CFC-11'//trim(long_suffix), ' ',               &
       missing_value = -1.0e+10)

  cfc(n)%id_kw_11 = register_diag_field('ocean_model',                  &
       'kw_11'//trim(suffix), grid%tracer_axes(1:2),                    &
       Time%model_time,                                                 &
       'Piston velocity - CFC-11'//trim(long_suffix), 'm s^-1',         &
       missing_value = -1.0e+10)

  cfc(n)%id_alpha_11 = register_diag_field('ocean_model',               &
       'alpha_11'//trim(suffix), grid%tracer_axes(1:2),                 &
       Time%model_time,                                                 &
       'Solubility CFC-11'//trim(long_suffix), 'mol m^-3 pptv^-1',      &
       missing_value = -1.0e+10)

  cfc(n)%id_csat_11 = register_diag_field('ocean_model',                &
       'csat_11'//trim(suffix), grid%tracer_axes(1:2),                  &
       Time%model_time,                                                 &
       'Saturation concentration - CFC-11'//trim(long_suffix),          &
       'mol m^-3',                                                      &
       missing_value = -1.0e+10)

  cfc(n)%id_sc_12 = register_diag_field('ocean_model',                  &
       'sc_12'//trim(suffix), grid%tracer_axes(1:2),                    &
       Time%model_time,                                                 &
       'Schmidt number - CFC-12'//trim(long_suffix), ' ',               &
       missing_value = -1.0e+10)

  cfc(n)%id_kw_12 = register_diag_field('ocean_model',                  &
       'kw_12'//trim(suffix), grid%tracer_axes(1:2),                    &
       Time%model_time,                                                 &
       'Piston velocity - CFC-12'//trim(long_suffix), 'm s^-1',         &
       missing_value = -1.0e+10)

  cfc(n)%id_alpha_12 = register_diag_field('ocean_model',               &
       'alpha_12'//trim(suffix), grid%tracer_axes(1:2),                 &
       Time%model_time,                                                 &
       'Solubility CFC-12'//trim(long_suffix), 'mol m^-3 pptv^-1',      &
       missing_value = -1.0e+10)

  cfc(n)%id_csat_12 = register_diag_field('ocean_model',                &
       'csat_12'//trim(suffix), grid%tracer_axes(1:2),                  &
       Time%model_time,                                                 &
       'Saturation concentration - CFC-12'//trim(long_suffix),          &
       'mol m^-3',                                                      &
       missing_value = -1.0e+10)

enddo  !} n

!
!-----------------------------------------------------------------------
!     give info
!-----------------------------------------------------------------------
!

write(stdout(),*)
write(stdout(),*) trim(note_header), ' Tracer runs initialized'
write(stdout(),*)

return

end subroutine  ocmip2_cfc_start  !}
! </SUBROUTINE> NAME="ocmip2_cfc_start"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_cfc_tracer">
!
! <DESCRIPTION>
!     Perform things that should be done in tracer, but are done here
! in order to minimize the number of hooks necessary in the MOM4 basecode
! </DESCRIPTION>
!

subroutine ocmip2_cfc_tracer  !{

!
!-----------------------------------------------------------------------
!     modules (have to come first)
!-----------------------------------------------------------------------
!

use mpp_mod, only : mpp_sum

implicit none

!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'ocmip2_cfc_tracer'
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

return

end subroutine  ocmip2_cfc_tracer  !}
! </SUBROUTINE> NAME="ocmip2_cfc_tracer"


!#######################################################################
! <SUBROUTINE NAME="read_cfc_timehist">
!
! <DESCRIPTION>
!
! =====================================================================
!
!     SUBROUTINE read_cfc_timehist
!
!     PURPOSE: reads in the atmospheric time histories of CFC from
!                  a specially prepared file. This file should contain
!                  five columns with the following data entries
!                  Year, CFC-11(NH),CFC-12(NH),
!                        CFC-11(SH),CFC-12(SH)
! 
! =====================================================================
!
! </DESCRIPTION>
 
subroutine read_cfc_timehist  !{

use mpp_io_mod, only : mpp_close, mpp_open, MPP_RDONLY,         &
     MPP_ASCII, MPP_SINGLE

implicit none

!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'read_cfc_timehist'
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

integer                 :: n
integer                 :: lun
character(len=1)        :: dummy
logical                 :: startdata = .false.

!
! =====================================================================
!     begin of executable code
! =====================================================================
!

if (atmcfc_file .ne. ' ') then  !{

!
!-----------------------------------------------------------------------
!     open file
!-----------------------------------------------------------------------
!

  call mpp_open(lun, atmcfc_file, action=MPP_RDONLY,         &
                form=MPP_ASCII, threading=MPP_SINGLE)

!
!-----------------------------------------------------------------------
!     discard header
!-----------------------------------------------------------------------
!

  do while (.not. startdata)              
    read(lun,'(a1)') dummy
    if (dummy .eq. 'S') then
      startdata = .true.
    endif
  enddo

!
!-----------------------------------------------------------------------
!     read in atmospheric time histories
!-----------------------------------------------------------------------
!

  n = 0
  do  !{
    n = n + 1
    if (n .gt. max_cfc_rec) then  !{
      write (stdout(),*) trim(error_header), ' Trying to read record ', n
      call mpp_error(FATAL, trim(error_header) // ' Too many records in file')
    endif  !}
    read (lun,*,end=1) year_cfc_rec(n), pcfc_11_n_rec(n), pcfc_12_n_rec(n),    &
         pcfc_11_s_rec(n), pcfc_12_s_rec(n)
  enddo  !}
1 continue
  num_cfc_rec = n - 1
  min_cfc_year = year_cfc_rec(1)
  max_cfc_year = year_cfc_rec(num_cfc_rec)

!
!-----------------------------------------------------------------------
!     release unit number
!-----------------------------------------------------------------------
!

  call mpp_close(lun)

!
!-----------------------------------------------------------------------
!     atmcfc_file has not been defined
!-----------------------------------------------------------------------
!

else  !}{
  call mpp_error(FATAL, trim(error_header) // ' Name of file ' //        &
       ' containing atmospheric time histories has not been defined')
endif  !}

return
end subroutine  read_cfc_timehist  !}
! </SUBROUTINE> NAME="read_cfc_timehist"

end module  ocmip2_cfc_mod  !}
