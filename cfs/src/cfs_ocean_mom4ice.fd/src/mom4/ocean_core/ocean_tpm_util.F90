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
module ocean_tpm_util_mod  !{
! 
!<CONTACT EMAIL="Richard.Slater@noaa.gov"> Richard D. Slater
!</CONTACT>
!
!<REVIEWER EMAIL="John.Dunne@noaa.gov"> John P. Dunne
!</REVIEWER>
!
!<OVERVIEW>
! Ocean tracer package module pointers module
!</OVERVIEW>
!
!<DESCRIPTION>
! This module allocates a suite of variables used in ocean_tpm
!</DESCRIPTION>
!
! <INFO>
! </INFO>
!
! $Id$
!

use field_manager_mod, only: fm_string_len, fm_path_name_len, fm_field_name_len, fm_type_name_len
use field_manager_mod, only: fm_get_type, fm_get_index, fm_get_length
use field_manager_mod, only: fm_get_current_list, fm_new_list, fm_change_list, fm_loop_over_list
use field_manager_mod, only: fm_new_value, fm_get_value
use field_manager_mod, only: fm_exists
use fms_mod,           only: FATAL, WARNING, NOTE, stdout
use mpp_mod,           only: mpp_error

use ocean_types_mod,   only: ocean_time_type, ocean_domain_type, ocean_grid_type
use ocean_types_mod,   only: ocean_prog_tracer_type, ocean_diag_tracer_type

implicit none

private

public  otpm_start_namelist
public  otpm_end_namelist
public  otpm_check_for_bad_fields
public  otpm_set_caller
public  otpm_reset_caller
public  otpm_set_no_overwrite
public  otpm_reset_no_overwrite
public  otpm_set_good_name_list
public  otpm_reset_good_name_list
public  otpm_set_tracer_package
public  otpm_set_prog_tracer
public  otpm_set_diag_tracer
public  otpm_get_integer
public  otpm_get_logical
public  otpm_get_real
public  otpm_get_string
public  otpm_get_integer_array
public  otpm_get_logical_array
public  otpm_get_real_array
public  otpm_get_string_array
public  otpm_set_value
public  otpm_set_value_integer_array
public  otpm_set_value_logical_array
public  otpm_set_value_real_array
public  otpm_set_value_string_array
public  otpm_set_value_integer
public  otpm_set_value_logical
public  otpm_set_value_real
public  otpm_set_value_string

!
!       Public variables (formerly in ocean_tpm_pointers)
!

type(ocean_domain_type), public, pointer, save                          :: domain => NULL()
real, public                                                            :: dtts
type(ocean_grid_type), public, pointer, save                            :: grid => NULL()
type(ocean_prog_tracer_type), public, pointer, dimension(:), save       :: t_prog => NULL()
type(ocean_diag_tracer_type), public, pointer, dimension(:), save       :: t_diag => NULL()
type(ocean_time_type), public, pointer, save                            :: time => NULL()
integer, public                                                         :: imonth
integer, public                                                         :: iyear
logical, public                                                         :: end_of_day
logical, public                                                         :: end_of_month
logical, public                                                         :: end_of_year
logical, public                                                         :: mid_month
integer, public                                                         :: isd
integer, public                                                         :: ied
integer, public                                                         :: jsd
integer, public                                                         :: jed
integer, public                                                         :: isc
integer, public                                                         :: iec
integer, public                                                         :: jsc
integer, public                                                         :: jec
integer, public                                                         :: nk
integer, public                                                         :: taum1
integer, public                                                         :: tau
integer, public                                                         :: taup1
integer, public                                                         :: indtemp
integer, public                                                         :: indsal

!
!       private parameters
!

character(len=48), parameter            :: mod_name = 'ocean_tpm_util_mod'
character(len=fm_string_len), parameter :: default_units = ' '
real, parameter                         :: default_conversion = 1.0
real, parameter                         :: default_offset = 0.0
real, parameter                         :: default_min_tracer = -1.0e+20
real, parameter                         :: default_max_tracer = +1.0e+20
logical, parameter                      :: default_use_only_advection = .false.
real, parameter                         :: default_min_range = 1.0
real, parameter                         :: default_max_range = 0.0
character(len=fm_string_len), parameter :: default_file_in = 'INPUT/ocean_tracer.res.nc'
character(len=fm_string_len), parameter :: default_name_in = ' '
real, parameter                         :: default_scale_in = 1.0
real, parameter                         :: default_additive_in = 0.0
logical, parameter                      :: default_init = .false.
logical, parameter                      :: default_const_init_tracer = .false.
real, parameter                         :: default_const_init_value = 0.0
character(len=fm_string_len), parameter :: default_file_out = 'RESTART/ocean_tracer.res.nc'
character(len=fm_string_len), parameter :: default_flux_units = ' '
real, parameter                         :: default_min_flux_range = 1.0
real, parameter                         :: default_max_flux_range = 0.0
real, parameter                         :: default_min_tracer_limit = -1.0e+20
real, parameter                         :: default_max_tracer_limit = +1.0e+20
character(len=fm_string_len), parameter :: default_vert_adv_scheme = 'quicker'
character(len=fm_string_len), parameter :: default_horiz_adv_scheme = 'quicker'

!
!       Private variables
!

character(len=128)              :: default_caller = ' '
character(len=128)              :: save_default_caller = ' '
character(len=128)              :: default_good_name_list = ' '
character(len=128)              :: save_default_good_name_list = ' '
logical                         :: default_no_overwrite = .false.
logical                         :: save_default_no_overwrite = .false.
character(len=fm_path_name_len) :: save_current_list
character(len=fm_path_name_len) :: save_path
character(len=fm_path_name_len) :: save_name

!
!        Interface definitions for overloaded routines
!

!interface  otpm_get_value  !{
  !module procedure  otpm_get_value_integer
  !module procedure  otpm_get_value_logical
  !module procedure  otpm_get_value_real
  !module procedure  otpm_get_value_string
  !module procedure  otpm_get_value_integer_array
  !module procedure  otpm_get_value_logical_array
  !module procedure  otpm_get_value_real_array
  !module procedure  otpm_get_value_string_array
!end interface  !}

interface  otpm_set_value  !{
  module procedure  otpm_set_value_integer_array
  module procedure  otpm_set_value_logical_array
  module procedure  otpm_set_value_real_array
  module procedure  otpm_set_value_string_array
  module procedure  otpm_set_value_integer
  module procedure  otpm_set_value_logical
  module procedure  otpm_set_value_real
  module procedure  otpm_set_value_string
end interface  !}

interface  set_prog_value  !{
  module procedure  set_prog_value_integer
  module procedure  set_prog_value_logical
  module procedure  set_prog_value_real
  module procedure  set_prog_value_string
end interface  !}

contains


!#######################################################################
! <SUBROUTINE NAME="otpm_set_caller">
!
! <DESCRIPTION>
! Set the default value for the optional "caller" variable used in many of these
! subroutines. If the argument is blank, then set the default to blank, otherwise
! the deault will have brackets placed around the argument.
!
! </DESCRIPTION>
!

subroutine otpm_set_caller(caller)  !{

implicit none

!
!       arguments
!

character(len=*), intent(in)          :: caller

!
!       Local parameters
!

character(len=48), parameter  :: sub_name = 'otpm_set_caller'

!
!       Local variables
!

!
!       save the default caller string
!

save_default_caller = default_caller

!
!       set the default caller string
!

if (caller .eq. ' ') then  !{
  default_caller = ' '
else  !}{
  default_caller = '[' // trim(caller) // ']'
endif  !}

return

end subroutine otpm_set_caller  !}
! </SUBROUTINE> NAME="otpm_set_caller"


!#######################################################################
! <SUBROUTINE NAME="otpm_reset_caller">
!
! <DESCRIPTION>
! Reset the default value for the optional "caller" variable used in many of these
! subroutines to blank.
!
! </DESCRIPTION>
!

subroutine otpm_reset_caller  !{

implicit none

!
!       arguments
!

!
!       Local parameters
!

character(len=48), parameter  :: sub_name = 'otpm_reset_caller'

!
!       Local variables
!

!
!       reset the default caller string
!

default_caller = save_default_caller
save_default_caller = ' '

return

end subroutine otpm_reset_caller  !}
! </SUBROUTINE> NAME="otpm_reset_caller"


!#######################################################################
! <SUBROUTINE NAME="otpm_set_good_name_list">
!
! <DESCRIPTION>
! Set the default value for the optional "good_name_list" variable used in many of these
! subroutines.
!
! </DESCRIPTION>
!

subroutine otpm_set_good_name_list(good_name_list)  !{

implicit none

!
!       arguments
!

character(len=*), intent(in)          :: good_name_list

!
!       Local parameters
!

character(len=48), parameter  :: sub_name = 'otpm_set_good_name_list'

!
!       Local variables
!

!
!       save the default good_name_list string
!

save_default_good_name_list = default_good_name_list

!
!       set the default good_name_list string
!

default_good_name_list = good_name_list

return

end subroutine otpm_set_good_name_list  !}
! </SUBROUTINE> NAME="otpm_set_good_name_list"


!#######################################################################
! <SUBROUTINE NAME="otpm_reset_good_name_list">
!
! <DESCRIPTION>
! Reset the default value for the optional "good_name_list" variable used in many of these
! subroutines to the saved value.
!
! </DESCRIPTION>
!

subroutine otpm_reset_good_name_list  !{

implicit none

!
!       arguments
!

!
!       Local parameters
!

character(len=48), parameter  :: sub_name = 'otpm_reset_good_name_list'

!
!       Local variables
!

!
!       reset the default good_name_list string
!

default_good_name_list = save_default_good_name_list
save_default_good_name_list = ' '

return

end subroutine otpm_reset_good_name_list  !}
! </SUBROUTINE> NAME="otpm_reset_good_name_list"


!#######################################################################
! <SUBROUTINE NAME="otpm_set_no_overwrite">
!
! <DESCRIPTION>
! Set the default value for the optional "no_overwrite" variable used in some of these
! subroutines.
!
! </DESCRIPTION>
!

subroutine otpm_set_no_overwrite(no_overwrite)  !{

implicit none

!
!       arguments
!

logical, intent(in)          :: no_overwrite

!
!       Local parameters
!

character(len=48), parameter  :: sub_name = 'otpm_set_no_overwrite'

!
!       Local variables
!

!
!       save the default no_overwrite string
!

save_default_no_overwrite = default_no_overwrite

!
!       set the default no_overwrite value
!

default_no_overwrite = no_overwrite

return

end subroutine otpm_set_no_overwrite  !}
! </SUBROUTINE> NAME="otpm_set_no_overwrite"


!#######################################################################
! <SUBROUTINE NAME="otpm_reset_no_overwrite">
!
! <DESCRIPTION>
! Reset the default value for the optional "no_overwrite" variable used in some of these
! subroutines to false.
!
! </DESCRIPTION>
!

subroutine otpm_reset_no_overwrite  !{

implicit none

!
!       arguments
!

!
!       Local parameters
!

character(len=48), parameter  :: sub_name = 'otpm_reset_no_overwrite'

!
!       Local variables
!

!
!       reset the default no_overwrite value
!

default_no_overwrite = save_default_no_overwrite
save_default_no_overwrite = .false.

return

end subroutine otpm_reset_no_overwrite  !}
! </SUBROUTINE> NAME="otpm_reset_no_overwrite"


!#######################################################################
! <SUBROUTINE NAME="otpm_check_for_bad_fields">
!
! <DESCRIPTION>
! Check for unrecognized fields in a list
!
! </DESCRIPTION>
!

subroutine otpm_check_for_bad_fields(list, good_fields, caller)  !{

implicit none

!
!       arguments
!

character(len=*), intent(in)                    :: list
character(len=*), intent(in), dimension(:)      :: good_fields
character(len=*), intent(in), optional          :: caller

!
!       Local parameters
!

character(len=48), parameter  :: sub_name = 'otpm_check_for_bad_fields'

!
!       Local variables
!

integer                                 :: i
integer                                 :: ind
integer                                 :: list_length
integer                                 :: good_length
character(len=fm_type_name_len)         :: typ
character(len=fm_field_name_len)        :: name
logical                                 :: found
character(len=256)                      :: error_header
character(len=256)                      :: warn_header
character(len=256)                      :: note_header
character(len=128)                      :: caller_str

!
!       set the caller string and headers
!

if (present(caller)) then  !{
  caller_str = '[' // trim(caller) // ']'
else  !}{
  caller_str = default_caller
endif  !}

error_header = '==>Error from ' // trim(mod_name) //   &
               '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a list is given (fatal if not)
!

if (list .eq. ' ') then  !{
  write (stdout(),*) trim(error_header) // ' Empty list given'
  call mpp_error(FATAL, trim(error_header) // ' Empty list given')
endif  !}

!
!       Check that we have been given a list
!

if (fm_get_type(list) .ne. 'list') then  !{
  write (stdout(),*) trim(error_header) // ' Not given a list: ' // trim(list)
  call mpp_error(FATAL, trim(error_header) // ' Not given a list: ' // trim(list))
endif  !}

!
!       Get the list length
!

list_length = fm_get_length(list)
if (list_length .lt. 0) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Problem getting length of ' // trim(list))
endif  !}

!
!       Get the number of good fields
!

good_length = size(good_fields)

if (list_length .lt. good_length) then  !{

!
!       If the list length is less than the number of good fields this is an error 
!       as the list should be fully populated and we'll check which extra fields
!       are given in good_fields
!

  write (stdout(),*) trim(error_header), 'List length < number of good fields (',       &
       list_length, ' < ', good_length, ') in list ', trim(list)

  do i = 1, good_length  !{
    if (.not. fm_exists(trim(list) // '/' // good_fields(i))) then  !{
      write (stdout(),*) trim(error_header), trim(good_fields(i)), ' is an extra good field'
    endif  !}
  enddo  !} i

  call mpp_error(FATAL, trim(error_header) //                                           &
       ' List length < number of good fields for list: ' // trim(list))

elseif (list_length .gt. good_length) then  !}{

!
!       If the list length is greater than the number of good fields this is an error
!       as the there should not be any more fields than those given in the good fields list
!       and we'll check which extra fields are given in the list
!

  write (stdout(),*) trim(warn_header), 'List length > number of good fields (',        &
       list_length, ' > ', good_length, ') in list ', trim(list)

  do while (fm_loop_over_list(list, name, typ, ind))  !{
    found = .false.
    do i = 1, good_length  !{
      found = found .or. (name .eq. good_fields(i))
    enddo  !} i
    if (.not. found) then  !{
      write (stdout(),*) trim(error_header), trim(name), ' is an extra list field'
    endif  !}
  enddo  !}

  call mpp_error(FATAL, trim(error_header) //                                           &
       ' List length > number of good fields for list: ' // trim(list))

endif  !}

!
!       If the list length equals the number of good fields then all is good
!

return

end subroutine otpm_check_for_bad_fields  !}
! </SUBROUTINE> NAME="otpm_check_for_bad_fields"


!#######################################################################
! <FUNCTION NAME="otpm_set_tracer_package">
!
! <DESCRIPTION>
! Set the values for a tracer package and return its index (0 on error)
!
! </DESCRIPTION>
!

function otpm_set_tracer_package(name, caller, units, conversion, offset,       &
     min_tracer, max_tracer, use_only_advection,                                &
     min_range, max_range, file_in, name_in,                                    &
     scale_in, additive_in, init,                                               &
     const_init_tracer, const_init_value, file_out, flux_units,                 &
     min_flux_range, max_flux_range,                                            &
     min_tracer_limit, max_tracer_limit,                                        &
     vert_adv_scheme, horiz_adv_scheme)                                         &
         result (pack_index)  !{

implicit none

!
!       Return type
!

integer :: pack_index

!
!       arguments
!

character(len=*), intent(in)            :: name
character(len=*), intent(in), optional  :: caller
character(len=*), intent(in), optional  :: units
character(len=*), intent(in), optional  :: flux_units
character(len=*), intent(in), optional  :: vert_adv_scheme
character(len=*), intent(in), optional  :: horiz_adv_scheme
real, intent(in), optional              :: conversion
character(len=*), intent(in), optional  :: file_in
character(len=*), intent(in), optional  :: file_out
logical, intent(in), optional           :: init
logical, intent(in), optional           :: use_only_advection
real, intent(in), optional              :: max_range
real, intent(in), optional              :: min_range
real, intent(in), optional              :: max_flux_range
real, intent(in), optional              :: min_flux_range
real, intent(in), optional              :: max_tracer
real, intent(in), optional              :: min_tracer
real, intent(in), optional              :: max_tracer_limit
real, intent(in), optional              :: min_tracer_limit
character(len=*), intent(in), optional  :: name_in
real, intent(in), optional              :: offset
real, intent(in), optional              :: scale_in
real, intent(in), optional              :: additive_in
logical, intent(in), optional           :: const_init_tracer 
real, intent(in), optional              :: const_init_value  

!
!       Local parameters
!

character(len=48), parameter  :: sub_name = 'otpm_set_tracer_package'

!
!       Local variables
!

character(len=fm_path_name_len)                         :: current_list
logical                                                 :: use_default
character(len=256)                                      :: error_header
character(len=256)                                      :: warn_header
character(len=256)                                      :: note_header
character(len=128)                                      :: caller_str
character(len=fm_path_name_len)                         :: tracer_package_name
character(len=fm_string_len), pointer, dimension(:)     :: good_list => NULL()

!
!       set the caller string and headers
!

if (present(caller)) then  !{
  caller_str = '[' // trim(caller) // ']'
else  !}{
  caller_str = default_caller
endif  !}

error_header = '==>Error from ' // trim(mod_name) //   &
               '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

write (stdout(),*)
write (stdout(),*) trim(note_header), ' Processing tracer package ', trim(name)

!
!       Check whether the package already exists. If so, then use that package
!

tracer_package_name = '/ocean_mod/tracer_packages/' // trim(name) // '/'
pack_index = fm_get_index(tracer_package_name)

if (pack_index .le. 0) then  !{

!
!       Set a new tracer package and get its index
!

  pack_index = fm_new_list(tracer_package_name)
  if (pack_index .le. 0) then  !{
    call mpp_error(FATAL, trim(error_header) // ' Could not set tracer package')
  endif  !}

endif  !}

!
!       Change to the new tracer package, first saving the current list
!

current_list = fm_get_current_list()
if (current_list .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not get the current list')
endif  !}

if (.not. fm_change_list(tracer_package_name)) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not change to the new list')
endif  !}

!
!       Set the array in which to save the valid names for this list,
!       used later for a consistency check. This is used in the otpm_set_value
!       routines to make the list of valid values
!

call otpm_set_good_name_list('/ocean_mod/GOOD/tracer_packages/' // trim(name) // '/good_list')

!
!       Set other defaults for the otpm_set_value routines
!

call otpm_set_no_overwrite(.true.)
call otpm_set_caller(caller_str)

!
!       Set the default number of instances (always zero)
!

call otpm_set_value('names', ' ', index = 0)

!
!       Set various values to given values, or to defaults if not given
!

if (present(units)) then  !{
  call otpm_set_value('units', units)
else  !}{
  call otpm_set_value('units', default_units, no_create = .true.)
endif  !}

if (present(conversion)) then  !{
  call otpm_set_value('conversion', conversion)
else  !}{
  call otpm_set_value('conversion', default_conversion, no_create = .true.)
endif  !}

if (present(offset)) then  !{
  call otpm_set_value('offset', offset)
else  !}{
  call otpm_set_value('offset', default_offset, no_create = .true.)
endif  !}

if (present(min_tracer)) then  !{
  call otpm_set_value('min_tracer', min_tracer)
else  !}{
  call otpm_set_value('min_tracer', default_min_tracer, no_create = .true.)
endif  !}

if (present(max_tracer)) then  !{
  call otpm_set_value('max_tracer', max_tracer)
else  !}{
  call otpm_set_value('max_tracer', default_max_tracer, no_create = .true.)
endif  !}

if (present(min_range)) then  !{
  call otpm_set_value('min_range', min_range)
else  !}{
  call otpm_set_value('min_range', default_min_range, no_create = .true.)
endif  !}

if (present(max_range)) then  !{
  call otpm_set_value('max_range', max_range)
else  !}{
  call otpm_set_value('max_range', default_max_range, no_create = .true.)
endif  !}

if (present(use_only_advection)) then  !{
  call otpm_set_value('use_only_advection', use_only_advection)
else  !}{
  call otpm_set_value('use_only_advection', default_use_only_advection, no_create = .true.)
endif  !}

if (present(file_in)) then  !{
  call otpm_set_value('file_in', file_in)
else  !}{
  call otpm_set_value('file_in', default_file_in, no_create = .true.)
endif  !}

if (present(init)) then  !{
  call otpm_set_value('init', init)
else  !}{
  call otpm_set_value('init', default_init, no_create = .true.)
endif  !}

if (present(name_in)) then  !{
  call otpm_set_value('name_in', name_in)
else  !}{
  call otpm_set_value('name_in', default_name_in, no_create = .true.)
endif  !}

if (present(scale_in)) then  !{
  call otpm_set_value('scale_in', scale_in)
else  !}{
  call otpm_set_value('scale_in', default_scale_in, no_create = .true.)
endif  !}

if (present(additive_in)) then  !{
  call otpm_set_value('additive_in', additive_in)
else  !}{
  call otpm_set_value('additive_in', default_additive_in, no_create = .true.)
endif  !}

if (present(const_init_tracer)) then  !{
  call otpm_set_value('const_init_tracer', const_init_tracer)
else  !}{
  call otpm_set_value('const_init_tracer', default_const_init_tracer, no_create = .true.)
endif  !}

if (present(const_init_value)) then  !{
  call otpm_set_value('const_init_value', const_init_value)
else  !}{
  call otpm_set_value('const_init_value', default_const_init_value, no_create = .true.)
endif  !}

if (present(file_out)) then  !{
  call otpm_set_value('file_out', file_out)
else  !}{
  call otpm_set_value('file_out', default_file_out, no_create = .true.)
endif  !}

if (present(flux_units)) then  !{
  call otpm_set_value('flux_units', flux_units)
else  !}{
  call otpm_set_value('flux_units', default_flux_units, no_create = .true.)
endif  !}

if (present(min_flux_range)) then  !{
  call otpm_set_value('min_flux_range', min_flux_range)
else  !}{
  call otpm_set_value('min_flux_range', default_min_flux_range, no_create = .true.)
endif  !}

if (present(max_flux_range)) then  !{
  call otpm_set_value('max_flux_range', max_flux_range)
else  !}{
  call otpm_set_value('max_flux_range', default_max_flux_range, no_create = .true.)
endif  !}

if (present(min_tracer_limit)) then  !{
  call otpm_set_value('min_tracer_limit', min_tracer_limit)
else  !}{
  call otpm_set_value('min_tracer_limit', default_min_tracer_limit, no_create = .true.)
endif  !}

if (present(max_tracer_limit)) then  !{
  call otpm_set_value('max_tracer_limit', max_tracer_limit)
else  !}{
  call otpm_set_value('max_tracer_limit', default_max_tracer_limit, no_create = .true.)
endif  !}

if (present(vert_adv_scheme)) then  !{
  call otpm_set_value('vertical-advection-scheme', vert_adv_scheme)
else  !}{
  call otpm_set_value('vertical-advection-scheme', default_vert_adv_scheme, no_create = .true.)
endif  !}

if (present(horiz_adv_scheme)) then  !{
  call otpm_set_value('horizontal-advection-scheme', horiz_adv_scheme)
else  !}{
  call otpm_set_value('horizontal-advection-scheme', default_horiz_adv_scheme, no_create = .true.)
endif  !}

!
!       Reset the defaults for the otpm_set_value calls
!

call otpm_reset_good_name_list
call otpm_reset_no_overwrite
call otpm_reset_caller

!
!       Change back to the saved current list
!

if (.not. fm_change_list(current_list)) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not change back to ' // trim(current_list))
endif  !}

!
!       Check for any errors in the number of fields in this list
!

if (caller_str .eq. ' ') then  !{
  caller_str = trim(mod_name) // '(' // trim(sub_name) // ')'
endif  !}
good_list => otpm_get_string_array('/ocean_mod/GOOD/tracer_packages/' // trim(name) // '/good_list', &
     caller = caller_str)
if (associated(good_list)) then  !{
  call otpm_check_for_bad_fields('/ocean_mod/tracer_packages/' // trim(name), good_list, caller = caller_str)
  deallocate(good_list)
else  !}{
  call mpp_error(FATAL, trim(error_header) // ' Empty "' // trim(name) // '" list')
endif  !}

!
!       Add the package name to the list of good packages, to be used
!       later for a consistency check
!

if (fm_new_value('/ocean_mod/GOOD/good_tracer_packages', name, append = .true.) .le. 0) then  !{
  call mpp_error(FATAL, trim(error_header) //                           &
       ' Could not add ' // trim(name) // ' to "good_tracer_packages" list')
endif  !}

return

end function otpm_set_tracer_package  !}
! </FUNCTION> NAME="otpm_set_tracer_package"


!#######################################################################
! <FUNCTION NAME="otpm_get_integer_array">
!
! <DESCRIPTION>
! Get an integer value from the Field Manager tree.
! </DESCRIPTION>
!
function otpm_get_integer_array(name, caller)            &
         result (array)  !{

implicit none

!
!       Return type
!

integer, pointer, dimension(:) :: array

!
!       arguments
!

character(len=*), intent(in)            :: name
character(len=*), intent(in), optional  :: caller

!
!       Local parameters
!

character(len=48), parameter  :: sub_name = 'otpm_get_integer_array'

!
!       Local variables
!

character(len=256)              :: error_header
character(len=256)              :: warn_header
character(len=256)              :: note_header
character(len=128)              :: caller_str
character(len=32)               :: index_str
character(len=fm_type_name_len) :: fm_type
integer                         :: i
integer                         :: length

nullify(array)

!
!       set the caller string and headers
!

if (present(caller)) then  !{
  caller_str = '[' // trim(caller) // ']'
else  !}{
  caller_str = default_caller
endif  !}

error_header = '==>Error from ' // trim(mod_name) //   &
               '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

fm_type = fm_get_type(name)
if (fm_type .eq. 'integer') then  !{
  length = fm_get_length(name)
  if (length .lt. 0) then  !{
    call mpp_error(FATAL, trim(error_header) // ' Problem getting length of ' // trim(name))
  endif  !}
  if (length .gt. 0) then  !{
    allocate(array(length))
    do i = 1, length  !{
      if (.not. fm_get_value(name, array(i), index = i)) then  !{
        write (index_str,*) '(', i, ')'
        call mpp_error(FATAL, trim(error_header) // ' Problem getting ' // trim(name) // trim(index_str))
      endif  !}
    enddo  !} i
  endif  !}
elseif (fm_type .eq. ' ') then  !}{
  call mpp_error(FATAL, trim(error_header) // ' Array does not exist: ' // trim(name))
else  !}{
 call mpp_error(FATAL, trim(error_header) // ' Wrong type for ' // trim(name) // ', found (' // trim(fm_type) // ')')
endif  !}

return

end function otpm_get_integer_array  !}
! </FUNCTION> NAME="otpm_get_integer_array"


!#######################################################################
! <FUNCTION NAME="otpm_get_logical_array">
!
! <DESCRIPTION>
! Get a logical value from the Field Manager tree.
! </DESCRIPTION>
!
function otpm_get_logical_array(name, caller)            &
         result (array)  !{

implicit none

!
!       Return type
!

logical, pointer, dimension(:) :: array

!
!       arguments
!

character(len=*), intent(in)            :: name
character(len=*), intent(in), optional  :: caller

!
!       Local parameters
!

character(len=48), parameter  :: sub_name = 'otpm_get_logical_array'

!
!       Local variables
!

character(len=256)              :: error_header
character(len=256)              :: warn_header
character(len=256)              :: note_header
character(len=128)              :: caller_str
character(len=32)               :: index_str
character(len=fm_type_name_len) :: fm_type
integer                         :: i
integer                         :: length

nullify(array)

!
!       set the caller string and headers
!

if (present(caller)) then  !{
  caller_str = '[' // trim(caller) // ']'
else  !}{
  caller_str = default_caller
endif  !}

error_header = '==>Error from ' // trim(mod_name) //   &
               '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

fm_type = fm_get_type(name)
if (fm_type .eq. 'logical') then  !{
  length = fm_get_length(name)
  if (length .lt. 0) then  !{
    call mpp_error(FATAL, trim(error_header) // ' Problem getting length of ' // trim(name))
  endif  !}
  if (length .gt. 0) then  !{
    allocate(array(length))
    do i = 1, length  !{
      if (.not. fm_get_value(name, array(i), index = i)) then  !{
        write (index_str,*) '(', i, ')'
        call mpp_error(FATAL, trim(error_header) // ' Problem getting ' // trim(name) // trim(index_str))
      endif  !}
    enddo  !} i
  endif  !}
elseif (fm_type .eq. ' ') then  !}{
  call mpp_error(FATAL, trim(error_header) // ' Array does not exist: ' // trim(name))
else  !}{
 call mpp_error(FATAL, trim(error_header) // ' Wrong type for ' // trim(name) // ', found (' // trim(fm_type) // ')')
endif  !}

return

end function otpm_get_logical_array  !}
! </FUNCTION> NAME="otpm_get_logical_array"


!#######################################################################
! <FUNCTION NAME="otpm_get_real_array">
!
! <DESCRIPTION>
! Get a real value from the Field Manager tree.
! </DESCRIPTION>
!
function otpm_get_real_array(name, caller)            &
         result (array)  !{

implicit none

!
!       Return type
!

real, pointer, dimension(:) :: array

!
!       arguments
!

character(len=*), intent(in)            :: name
character(len=*), intent(in), optional  :: caller

!
!       Local parameters
!

character(len=48), parameter  :: sub_name = 'otpm_get_real_array'

!
!       Local variables
!

character(len=256)              :: error_header
character(len=256)              :: warn_header
character(len=256)              :: note_header
character(len=128)              :: caller_str
character(len=32)               :: index_str
character(len=fm_type_name_len) :: fm_type
integer                         :: i
integer                         :: length

nullify(array)

!
!       set the caller string and headers
!

if (present(caller)) then  !{
  caller_str = '[' // trim(caller) // ']'
else  !}{
  caller_str = default_caller
endif  !}

error_header = '==>Error from ' // trim(mod_name) //   &
               '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

fm_type = fm_get_type(name)
if (fm_type .eq. 'real') then  !{
  length = fm_get_length(name)
  if (length .lt. 0) then  !{
    call mpp_error(FATAL, trim(error_header) // ' Problem getting length of ' // trim(name))
  endif  !}
  if (length .gt. 0) then  !{
    allocate(array(length))
    do i = 1, length  !{
      if (.not. fm_get_value(name, array(i), index = i)) then  !{
        write (index_str,*) '(', i, ')'
        call mpp_error(FATAL, trim(error_header) // ' Problem getting ' // trim(name) // trim(index_str))
      endif  !}
    enddo  !} i
  endif  !}
elseif (fm_type .eq. ' ') then  !}{
  call mpp_error(FATAL, trim(error_header) // ' Array does not exist: ' // trim(name))
else  !}{
 call mpp_error(FATAL, trim(error_header) // ' Wrong type for ' // trim(name) // ', found (' // trim(fm_type) // ')')
endif  !}

return

end function otpm_get_real_array  !}
! </FUNCTION> NAME="otpm_get_real_array"


!#######################################################################
! <FUNCTION NAME="otpm_get_string_array">
!
! <DESCRIPTION>
! Get a string value from the Field Manager tree.
! </DESCRIPTION>
!
function otpm_get_string_array(name, caller)            &
         result (array)  !{

implicit none

!
!       Return type
!

character(len=fm_string_len), pointer, dimension(:) :: array

!
!       arguments
!

character(len=*), intent(in)            :: name
character(len=*), intent(in), optional  :: caller

!
!       Local parameters
!

character(len=48), parameter  :: sub_name = 'otpm_get_string_array'

!
!       Local variables
!

character(len=256)              :: error_header
character(len=256)              :: warn_header
character(len=256)              :: note_header
character(len=128)              :: caller_str
character(len=32)               :: index_str
character(len=fm_type_name_len) :: fm_type
integer                         :: i
integer                         :: length

nullify(array)

!
!       set the caller string and headers
!

if (present(caller)) then  !{
  caller_str = '[' // trim(caller) // ']'
else  !}{
  caller_str = default_caller
endif  !}

error_header = '==>Error from ' // trim(mod_name) //   &
               '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

fm_type = fm_get_type(name)
if (fm_type .eq. 'string') then  !{
  length = fm_get_length(name)
  if (length .lt. 0) then  !{
    call mpp_error(FATAL, trim(error_header) // ' Problem getting length of ' // trim(name))
  endif  !}
  if (length .gt. 0) then  !{
    allocate(array(length))
    do i = 1, length  !{
      if (.not. fm_get_value(name, array(i), index = i)) then  !{
        write (index_str,*) '(', i, ')'
        call mpp_error(FATAL, trim(error_header) // ' Problem getting ' // trim(name) // trim(index_str))
      endif  !}
    enddo  !} i
  endif  !}
elseif (fm_type .eq. ' ') then  !}{
  call mpp_error(FATAL, trim(error_header) // ' Array does not exist: ' // trim(name))
else  !}{
 call mpp_error(FATAL, trim(error_header) // ' Wrong type for ' // trim(name) // ', found (' // trim(fm_type) // ')')
endif  !}

return

end function otpm_get_string_array  !}
! </FUNCTION> NAME="otpm_get_string_array"


!#######################################################################
! <FUNCTION NAME="otpm_get_integer">
!
! <DESCRIPTION>
! Get an integer value from the Field Manager tree.
! </DESCRIPTION>
!
function otpm_get_integer(name, caller, index, default_value, scalar)            &
         result (value)  !{

implicit none

!
!       Return type
!

integer :: value

!
!       arguments
!

character(len=*), intent(in)            :: name
character(len=*), intent(in), optional  :: caller
integer, intent(in), optional           :: index
integer, intent(in), optional           :: default_value
logical, intent(in), optional           :: scalar

!
!       Local parameters
!

character(len=48), parameter  :: sub_name = 'otpm_get_integer'

!
!       Local variables
!

character(len=256)              :: error_header
character(len=256)              :: warn_header
character(len=256)              :: note_header
character(len=128)              :: caller_str
integer                         :: index_t
character(len=fm_type_name_len) :: fm_type
integer                         :: field_length

!
!       set the caller string and headers
!

if (present(caller)) then  !{
  caller_str = '[' // trim(caller) // ']'
else  !}{
  caller_str = default_caller
endif  !}

error_header = '==>Error from ' // trim(mod_name) //   &
               '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

!
!       Check whether we require a scalar (length=1) and return
!       an error if we do, and it isn't
!

if (present(scalar)) then  !{
  if (scalar) then  !{
    field_length = fm_get_length(name)
    if (field_length .lt. 0) then  !{
      call mpp_error(FATAL, trim(error_header) // ' Problem getting length of ' // trim(name))
    elseif (field_length .gt. 1) then  !}{
      call mpp_error(FATAL, trim(error_header) // trim(name) // ' not scalar')
    endif  !}
  endif  !}
endif  !}

!
!       set the index
!

if (present(index)) then  !{
  index_t = index
  if (index .le. 0) then  !{
    call mpp_error(FATAL, trim(error_header) // ' Index not positive')
  endif  !}
else  !}{
  index_t = 1
endif  !}

fm_type = fm_get_type(name)
if (fm_type .eq. 'integer') then  !{
  if (.not. fm_get_value(name, value, index = index_t)) then  !{
    call mpp_error(FATAL, trim(error_header) // ' Problem getting ' // trim(name))
  endif  !}
elseif (fm_type .eq. ' ' .and. present(default_value)) then  !}{
  value = default_value
else  !}{
 call mpp_error(FATAL, trim(error_header) // ' Wrong type for ' // trim(name) // ', found (' // trim(fm_type) // ')')
endif  !}

return

end function otpm_get_integer  !}
! </FUNCTION> NAME="otpm_get_integer"


!#######################################################################
! <FUNCTION NAME="otpm_get_logical">
!
! <DESCRIPTION>
! Get a logical value from the Field Manager tree.
! </DESCRIPTION>
!
function otpm_get_logical(name, caller, index, default_value, scalar)            &
         result (value)  !{

implicit none

!
!       Return type
!

logical :: value

!
!       arguments
!

character(len=*), intent(in)            :: name
character(len=*), intent(in), optional  :: caller
integer, intent(in), optional           :: index
logical, intent(in), optional           :: default_value
logical, intent(in), optional           :: scalar

!
!       Local parameters
!

character(len=48), parameter  :: sub_name = 'otpm_get_logical'

!
!       Local variables
!

character(len=256)              :: error_header
character(len=256)              :: warn_header
character(len=256)              :: note_header
character(len=128)              :: caller_str
integer                         :: index_t
character(len=fm_type_name_len) :: fm_type
integer                         :: field_length

!
!       set the caller string and headers
!

if (present(caller)) then  !{
  caller_str = '[' // trim(caller) // ']'
else  !}{
  caller_str = default_caller
endif  !}

error_header = '==>Error from ' // trim(mod_name) //   &
               '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

!
!       Check whether we require a scalar (length=1) and return
!       an error if we do, and it isn't
!

if (present(scalar)) then  !{
  if (scalar) then  !{
    field_length = fm_get_length(name)
    if (field_length .lt. 0) then  !{
      call mpp_error(FATAL, trim(error_header) // ' Problem getting length of ' // trim(name))
    elseif (field_length .gt. 1) then  !}{
      call mpp_error(FATAL, trim(error_header) // trim(name) // ' not scalar')
    endif  !}
  endif  !}
endif  !}

!
!       set the index
!

if (present(index)) then  !{
  index_t = index
  if (index .le. 0) then  !{
    call mpp_error(FATAL, trim(error_header) // ' Index not positive')
  endif  !}
else  !}{
  index_t = 1
endif  !}

fm_type = fm_get_type(name)
if (fm_type .eq. 'logical') then  !{
  if (.not. fm_get_value(name, value, index = index_t)) then  !{
    call mpp_error(FATAL, trim(error_header) // ' Problem getting ' // trim(name))
  endif  !}
elseif (fm_type .eq. ' ' .and. present(default_value)) then  !}{
  value = default_value
else  !}{
 call mpp_error(FATAL, trim(error_header) // ' Wrong type for ' // trim(name) // ', found (' // trim(fm_type) // ')')
endif  !}

return

end function otpm_get_logical  !}
! </FUNCTION> NAME="otpm_get_logical"


!#######################################################################
! <FUNCTION NAME="otpm_get_real">
!
! <DESCRIPTION>
! Get a real value from the Field Manager tree.
! </DESCRIPTION>
!
function otpm_get_real(name, caller, index, default_value, scalar)            &
         result (value)  !{

implicit none

!
!       Return type
!

real :: value

!
!       arguments
!

character(len=*), intent(in)            :: name
character(len=*), intent(in), optional  :: caller
integer, intent(in), optional           :: index
real, intent(in), optional              :: default_value
logical, intent(in), optional           :: scalar

!
!       Local parameters
!

character(len=48), parameter  :: sub_name = 'otpm_get_real'

!
!       Local variables
!

character(len=256)              :: error_header
character(len=256)              :: warn_header
character(len=256)              :: note_header
character(len=128)              :: caller_str
integer                         :: index_t
character(len=fm_type_name_len) :: fm_type
integer                         :: field_length

!
!       set the caller string and headers
!

if (present(caller)) then  !{
  caller_str = '[' // trim(caller) // ']'
else  !}{
  caller_str = default_caller
endif  !}

error_header = '==>Error from ' // trim(mod_name) //   &
               '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

!
!       Check whether we require a scalar (length=1) and return
!       an error if we do, and it isn't
!

if (present(scalar)) then  !{
  if (scalar) then  !{
    field_length = fm_get_length(name)
    if (field_length .lt. 0) then  !{
      call mpp_error(FATAL, trim(error_header) // ' Problem getting length of ' // trim(name))
    elseif (field_length .gt. 1) then  !}{
      call mpp_error(FATAL, trim(error_header) // trim(name) // ' not scalar')
    endif  !}
  endif  !}
endif  !}

!
!       set the index
!

if (present(index)) then  !{
  index_t = index
  if (index .le. 0) then  !{
    call mpp_error(FATAL, trim(error_header) // ' Index not positive')
  endif  !}
else  !}{
  index_t = 1
endif  !}

fm_type = fm_get_type(name)
if (fm_type .eq. 'real') then  !{
  if (.not. fm_get_value(name, value, index = index_t)) then  !{
    call mpp_error(FATAL, trim(error_header) // ' Problem getting ' // trim(name))
  endif  !}
elseif (fm_type .eq. ' ' .and. present(default_value)) then  !}{
  value = default_value
else  !}{
 call mpp_error(FATAL, trim(error_header) // ' Wrong type for ' // trim(name) // ', found (' // trim(fm_type) // ')')
endif  !}

return

end function otpm_get_real  !}
! </FUNCTION> NAME="otpm_get_real"


!#######################################################################
! <FUNCTION NAME="otpm_get_string">
!
! <DESCRIPTION>
! Get a string value from the Field Manager tree.
! </DESCRIPTION>
!
function otpm_get_string(name, caller, index, default_value, scalar)            &
         result (value)  !{

implicit none

!
!       Return type
!

character(len=fm_string_len) :: value

!
!       arguments
!

character(len=*), intent(in)            :: name
character(len=*), intent(in), optional  :: caller
integer, intent(in), optional           :: index
character(len=*), intent(in), optional  :: default_value
logical, intent(in), optional           :: scalar

!
!       Local parameters
!

character(len=48), parameter  :: sub_name = 'otpm_get_string'

!
!       Local variables
!

character(len=256)              :: error_header
character(len=256)              :: warn_header
character(len=256)              :: note_header
character(len=128)              :: caller_str
integer                         :: index_t
character(len=fm_type_name_len) :: fm_type
integer                         :: field_length

!
!       set the caller string and headers
!

if (present(caller)) then  !{
  caller_str = '[' // trim(caller) // ']'
else  !}{
  caller_str = default_caller
endif  !}

error_header = '==>Error from ' // trim(mod_name) //   &
               '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

!
!       Check whether we require a scalar (length=1) and return
!       an error if we do, and it isn't
!

if (present(scalar)) then  !{
  if (scalar) then  !{
    field_length = fm_get_length(name)
    if (field_length .lt. 0) then  !{
      call mpp_error(FATAL, trim(error_header) // ' Problem getting length of ' // trim(name))
    elseif (field_length .gt. 1) then  !}{
      call mpp_error(FATAL, trim(error_header) // trim(name) // ' not scalar')
    endif  !}
  endif  !}
endif  !}

!
!       set the index
!

if (present(index)) then  !{
  index_t = index
  if (index .le. 0) then  !{
    call mpp_error(FATAL, trim(error_header) // ' Index not positive')
  endif  !}
else  !}{
  index_t = 1
endif  !}

fm_type = fm_get_type(name)
if (fm_type .eq. 'string') then  !{
  if (.not. fm_get_value(name, value, index = index_t)) then  !{
    call mpp_error(FATAL, trim(error_header) // ' Problem getting ' // trim(name))
  endif  !}
elseif (fm_type .eq. ' ' .and. present(default_value)) then  !}{
  value = default_value
else  !}{
 call mpp_error(FATAL, trim(error_header) // ' Wrong type for ' // trim(name) // ', found (' // trim(fm_type) // ')')
endif  !}

return

end function otpm_get_string  !}
! </FUNCTION> NAME="otpm_get_string"


!#######################################################################
! <SUBROUTINE NAME="otpm_set_value_integer_array">
!
! <DESCRIPTION>
! Set an integer array in the Field Manager tree.
! </DESCRIPTION>
!

subroutine otpm_set_value_integer_array(name, value, length, caller, no_overwrite, good_name_list)  !{

implicit none

!
!       arguments
!

character(len=*), intent(in)                            :: name
integer, intent(in)                                     :: length
integer, intent(in)                                     :: value(length)
character(len=*), intent(in), optional                  :: caller
logical, intent(in), optional                           :: no_overwrite
character(len=fm_path_name_len), intent(in), optional   :: good_name_list

!
!       Local parameters
!

character(len=48), parameter    :: sub_name = 'otpm_set_value_integer_array'

!
!       Local variables
!

character(len=256)              :: error_header
character(len=256)              :: warn_header
character(len=256)              :: note_header
character(len=128)              :: caller_str
character(len=32)               :: str_error
integer                         :: field_index
integer                         :: field_length
integer                         :: n
logical                         :: no_overwrite_use
character(len=fm_path_name_len) :: good_name_list_use

!
!       set the caller string and headers
!

if (present(caller)) then  !{
  caller_str = '[' // trim(caller) // ']'
else  !}{
  caller_str = default_caller
endif  !}

error_header = '==>Error from ' // trim(mod_name) //   &
               '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

!
!       check that the length is non-negative
!

if (length .lt. 0) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Negative array length')
endif  !}

!
!       check for whether to overwrite existing values
!

if (present(no_overwrite)) then  !{
  no_overwrite_use = no_overwrite
else  !}{
  no_overwrite_use = default_no_overwrite
endif  !}

!
!       check for whether to save the name in a list
!

if (present(good_name_list)) then  !{
  good_name_list_use = good_name_list
else  !}{
  good_name_list_use = default_good_name_list
endif  !}

!
!       write the data array
!

if (length .eq. 0) then  !{
  if (.not. (no_overwrite_use .and. fm_exists(name))) then  !{
    field_index = fm_new_value(name, 0, index = 0)
    if (field_index .le. 0) then  !{
      write (str_error,*) ' with length = ', length
      call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name) // trim(str_error))
    endif  !}
  endif  !}
else  !}{
  if (no_overwrite_use .and. fm_exists(name)) then  !{
    field_length = fm_get_length(name)
    if (field_length .lt. 0) then  !{
      call mpp_error(FATAL, trim(error_header) // ' Problem getting length of ' // trim(name))
    endif  !}
    do n = field_length + 1, length  !{
      field_index = fm_new_value(name, value(n), index = n)
      if (field_index .le. 0) then  !{
        write (str_error,*) ' with index = ', n
        call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name) // trim(str_error))
      endif  !}
    enddo  !} n
  else  !}{
    field_index = fm_new_value(name, value(1))
    if (field_index .le. 0) then  !{
      call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name))
    endif  !}
    do n = 2, length  !{
      field_index = fm_new_value(name, value(n), index = n)
      if (field_index .le. 0) then  !{
        write (str_error,*) ' with index = ', n
        call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name) // trim(str_error))
      endif  !}
    enddo  !} n
  endif  !}
endif  !}

!
!       Add the variable name to the list of good names, to be used
!       later for a consistency check
!

if (good_name_list_use .ne. ' ') then  !{
  if (fm_new_value(good_name_list_use, name, append = .true., create = .true.) .le. 0) then  !{
    call mpp_error(FATAL, trim(error_header) //                           &
         ' Could not add ' // trim(name) // ' to "' // trim(good_name_list_use) // '" list')
  endif  !}
endif  !}

return

end subroutine otpm_set_value_integer_array  !}
! </SUBROUTINE> NAME="otpm_set_value_integer_array"


!#######################################################################
! <SUBROUTINE NAME="otpm_set_value_logical_array">
!
! <DESCRIPTION>
! Set a logical array in the Field Manager tree.
! </DESCRIPTION>
!

subroutine otpm_set_value_logical_array(name, value, length, caller, no_overwrite, good_name_list)  !{

implicit none

!
!       arguments
!

character(len=*), intent(in)                            :: name
integer, intent(in)                                     :: length
logical, intent(in)                                     :: value(length)
character(len=*), intent(in), optional                  :: caller
logical, intent(in), optional                           :: no_overwrite
character(len=fm_path_name_len), intent(in), optional   :: good_name_list

!
!       Local parameters
!

character(len=48), parameter    :: sub_name = 'otpm_set_value_logical_array'

!
!       Local variables
!

character(len=256)              :: error_header
character(len=256)              :: warn_header
character(len=256)              :: note_header
character(len=128)              :: caller_str
character(len=32)               :: str_error
integer                         :: field_index
integer                         :: field_length
integer                         :: n
logical                         :: no_overwrite_use
character(len=fm_path_name_len) :: good_name_list_use

!
!       set the caller string and headers
!

if (present(caller)) then  !{
  caller_str = '[' // trim(caller) // ']'
else  !}{
  caller_str = default_caller
endif  !}

error_header = '==>Error from ' // trim(mod_name) //   &
               '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

!
!       check that the length is non-negative
!

if (length .lt. 0) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Negative array length')
endif  !}

!
!       check for whether to overwrite existing values
!

if (present(no_overwrite)) then  !{
  no_overwrite_use = no_overwrite
else  !}{
  no_overwrite_use = default_no_overwrite
endif  !}

!
!       check for whether to save the name in a list
!

if (present(good_name_list)) then  !{
  good_name_list_use = good_name_list
else  !}{
  good_name_list_use = default_good_name_list
endif  !}

!
!       write the data array
!

if (length .eq. 0) then  !{
  if (.not. (no_overwrite_use .and. fm_exists(name))) then  !{
    field_index = fm_new_value(name, .false., index = 0)
    if (field_index .le. 0) then  !{
      write (str_error,*) ' with length = ', length
      call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name) // trim(str_error))
    endif  !}
  endif  !}
else  !}{
  if (no_overwrite_use .and. fm_exists(name)) then  !{
    field_length = fm_get_length(name)
    if (field_length .lt. 0) then  !{
      call mpp_error(FATAL, trim(error_header) // ' Problem getting length of ' // trim(name))
    endif  !}
    do n = field_length + 1, length  !{
      field_index = fm_new_value(name, value(n), index = n)
      if (field_index .le. 0) then  !{
        write (str_error,*) ' with index = ', n
        call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name) // trim(str_error))
      endif  !}
    enddo  !} n
  else  !}{
    field_index = fm_new_value(name, value(1))
    if (field_index .le. 0) then  !{
      call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name))
    endif  !}
    do n = 2, length  !{
      field_index = fm_new_value(name, value(n), index = n)
      if (field_index .le. 0) then  !{
        write (str_error,*) ' with index = ', n
        call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name) // trim(str_error))
      endif  !}
    enddo  !} n
  endif  !}
endif  !}

!
!       Add the variable name to the list of good names, to be used
!       later for a consistency check
!

if (good_name_list_use .ne. ' ') then  !{
  if (fm_new_value(good_name_list_use, name, append = .true., create = .true.) .le. 0) then  !{
    call mpp_error(FATAL, trim(error_header) //                           &
         ' Could not add ' // trim(name) // ' to "' // trim(good_name_list_use) // '" list')
  endif  !}
endif  !}

return

end subroutine otpm_set_value_logical_array  !}
! </SUBROUTINE> NAME="otpm_set_value_logical_array"


!#######################################################################
! <SUBROUTINE NAME="otpm_set_value_real_array">
!
! <DESCRIPTION>
! Set a real array in the Field Manager tree.
! </DESCRIPTION>
!

subroutine otpm_set_value_real_array(name, value, length, caller, no_overwrite, good_name_list)  !{

implicit none

!
!       arguments
!

character(len=*), intent(in)                            :: name
integer, intent(in)                                     :: length
real, intent(in)                                        :: value(length)
character(len=*), intent(in), optional                  :: caller
logical, intent(in), optional                           :: no_overwrite
character(len=fm_path_name_len), intent(in), optional   :: good_name_list

!
!       Local parameters
!

character(len=48), parameter    :: sub_name = 'otpm_set_value_real_array'

!
!       Local variables
!

character(len=256)              :: error_header
character(len=256)              :: warn_header
character(len=256)              :: note_header
character(len=128)              :: caller_str
character(len=32)               :: str_error
integer                         :: field_index
integer                         :: field_length
integer                         :: n
logical                         :: no_overwrite_use
character(len=fm_path_name_len) :: good_name_list_use

!
!       set the caller string and headers
!

if (present(caller)) then  !{
  caller_str = '[' // trim(caller) // ']'
else  !}{
  caller_str = default_caller
endif  !}

error_header = '==>Error from ' // trim(mod_name) //   &
               '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

!
!       check that the length is non-negative
!

if (length .lt. 0) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Negative array length')
endif  !}

!
!       check for whether to overwrite existing values
!

if (present(no_overwrite)) then  !{
  no_overwrite_use = no_overwrite
else  !}{
  no_overwrite_use = default_no_overwrite
endif  !}

!
!       check for whether to save the name in a list
!

if (present(good_name_list)) then  !{
  good_name_list_use = good_name_list
else  !}{
  good_name_list_use = default_good_name_list
endif  !}

!
!       write the data array
!

if (length .eq. 0) then  !{
  if (.not. (no_overwrite_use .and. fm_exists(name))) then  !{
    field_index = fm_new_value(name, 0.0, index = 0)
    if (field_index .le. 0) then  !{
      write (str_error,*) ' with length = ', length
      call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name) // trim(str_error))
    endif  !}
  endif  !}
else  !}{
  if (no_overwrite_use .and. fm_exists(name)) then  !{
    field_length = fm_get_length(name)
    if (field_length .lt. 0) then  !{
      call mpp_error(FATAL, trim(error_header) // ' Problem getting length of ' // trim(name))
    endif  !}
    do n = field_length + 1, length  !{
      field_index = fm_new_value(name, value(n), index = n)
      if (field_index .le. 0) then  !{
        write (str_error,*) ' with index = ', n
        call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name) // trim(str_error))
      endif  !}
    enddo  !} n
  else  !}{
    field_index = fm_new_value(name, value(1))
    if (field_index .le. 0) then  !{
      call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name))
    endif  !}
    do n = 2, length  !{
      field_index = fm_new_value(name, value(n), index = n)
      if (field_index .le. 0) then  !{
        write (str_error,*) ' with index = ', n
        call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name) // trim(str_error))
      endif  !}
    enddo  !} n
  endif  !}
endif  !}

!
!       Add the variable name to the list of good names, to be used
!       later for a consistency check
!

if (good_name_list_use .ne. ' ') then  !{
  if (fm_new_value(good_name_list_use, name, append = .true., create = .true.) .le. 0) then  !{
    call mpp_error(FATAL, trim(error_header) //                           &
         ' Could not add ' // trim(name) // ' to "' // trim(good_name_list_use) // '" list')
  endif  !}
endif  !}

return

end subroutine otpm_set_value_real_array  !}
! </SUBROUTINE> NAME="otpm_set_value_real_array"


!#######################################################################
! <SUBROUTINE NAME="otpm_set_value_string_array">
!
! <DESCRIPTION>
! Set a string array in the Field Manager tree.
! </DESCRIPTION>
!

subroutine otpm_set_value_string_array(name, value, length, caller, no_overwrite, good_name_list)  !{

implicit none

!
!       arguments
!

character(len=*), intent(in)                            :: name
integer, intent(in)                                     :: length
character(len=*), intent(in)                            :: value(length)
character(len=*), intent(in), optional                  :: caller
logical, intent(in), optional                           :: no_overwrite
character(len=fm_path_name_len), intent(in), optional   :: good_name_list

!
!       Local parameters
!

character(len=48), parameter    :: sub_name = 'otpm_set_value_string_array'

!
!       Local variables
!

character(len=256)              :: error_header
character(len=256)              :: warn_header
character(len=256)              :: note_header
character(len=128)              :: caller_str
character(len=32)               :: str_error
integer                         :: field_index
integer                         :: field_length
integer                         :: n
logical                         :: no_overwrite_use
character(len=fm_path_name_len) :: good_name_list_use

!
!       set the caller string and headers
!

if (present(caller)) then  !{
  caller_str = '[' // trim(caller) // ']'
else  !}{
  caller_str = default_caller
endif  !}

error_header = '==>Error from ' // trim(mod_name) //   &
               '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

!
!       check that the length is non-negative
!

if (length .lt. 0) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Negative array length')
endif  !}

!
!       check for whether to overwrite existing values
!

if (present(no_overwrite)) then  !{
  no_overwrite_use = no_overwrite
else  !}{
  no_overwrite_use = default_no_overwrite
endif  !}

!
!       check for whether to save the name in a list
!

if (present(good_name_list)) then  !{
  good_name_list_use = good_name_list
else  !}{
  good_name_list_use = default_good_name_list
endif  !}

!
!       write the data array
!

if (length .eq. 0) then  !{
  if (.not. (no_overwrite_use .and. fm_exists(name))) then  !{
    field_index = fm_new_value(name, ' ', index = 0)
    if (field_index .le. 0) then  !{
      write (str_error,*) ' with length = ', length
      call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name) // trim(str_error))
    endif  !}
  endif  !}
else  !}{
  if (no_overwrite_use .and. fm_exists(name)) then  !{
    field_length = fm_get_length(name)
    if (field_length .lt. 0) then  !{
      call mpp_error(FATAL, trim(error_header) // ' Problem getting length of ' // trim(name))
    endif  !}
    do n = field_length + 1, length  !{
      field_index = fm_new_value(name, value(n), index = n)
      if (field_index .le. 0) then  !{
        write (str_error,*) ' with index = ', n
        call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name) // trim(str_error))
      endif  !}
    enddo  !} n
  else  !}{
    field_index = fm_new_value(name, value(1))
    if (field_index .le. 0) then  !{
      call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name))
    endif  !}
    do n = 2, length  !{
      field_index = fm_new_value(name, value(n), index = n)
      if (field_index .le. 0) then  !{
        write (str_error,*) ' with index = ', n
        call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name) // trim(str_error))
      endif  !}
    enddo  !} n
  endif  !}
endif  !}

!
!       Add the variable name to the list of good names, to be used
!       later for a consistency check
!

if (good_name_list_use .ne. ' ') then  !{
  if (fm_new_value(good_name_list_use, name, append = .true., create = .true.) .le. 0) then  !{
    call mpp_error(FATAL, trim(error_header) //                           &
         ' Could not add ' // trim(name) // ' to "' // trim(good_name_list_use) // '" list')
  endif  !}
endif  !}

return

end subroutine otpm_set_value_string_array  !}
! </SUBROUTINE> NAME="otpm_set_value_string_array"


!#######################################################################
! <SUBROUTINE NAME="otpm_set_value_integer">
!
! <DESCRIPTION>
! Set an integer value in the Field Manager tree.
! </DESCRIPTION>
!

subroutine otpm_set_value_integer(name, value, caller, index, append, no_create,        &
     no_overwrite, good_name_list)  !{

implicit none

!
!       arguments
!

character(len=*), intent(in)            :: name
integer, intent(in)                     :: value
character(len=*), intent(in), optional  :: caller
integer, intent(in), optional           :: index
logical, intent(in), optional           :: append
logical, intent(in), optional           :: no_create
logical, intent(in), optional           :: no_overwrite
character(len=*), intent(in), optional  :: good_name_list

!
!       Local parameters
!

character(len=48), parameter    :: sub_name = 'otpm_set_value_integer'

!
!       Local variables
!

character(len=256)              :: error_header
character(len=256)              :: warn_header
character(len=256)              :: note_header
character(len=128)              :: caller_str
character(len=32)               :: str_error
integer                         :: field_index
logical                         :: no_overwrite_use
integer                         :: field_length
character(len=fm_path_name_len) :: good_name_list_use
logical                         :: create
logical                         :: add_name

!
!       set the caller string and headers
!

if (present(caller)) then  !{
  caller_str = '[' // trim(caller) // ']'
else  !}{
  caller_str = default_caller
endif  !}

error_header = '==>Error from ' // trim(mod_name) //   &
               '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

!
!       check that append and index are not both given
!

if (present(index) .and. present(append)) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Append and index both given as arguments')
endif  !}

!
!       check for whether to overwrite existing values
!

if (present(no_overwrite)) then  !{
  no_overwrite_use = no_overwrite
else  !}{
  no_overwrite_use = default_no_overwrite
endif  !}

!
!       check for whether to save the name in a list
!

if (present(good_name_list)) then  !{
  good_name_list_use = good_name_list
else  !}{
  good_name_list_use = default_good_name_list
endif  !}

if (present(no_create)) then  !{
  create = .not. no_create
  if (no_create .and. (present(append) .or. present(index))) then  !{
    call mpp_error(FATAL, trim(error_header) // ' append or index are present when no_create is true for ' // trim(name))
  endif  !}
else  !}{
  create = .true.
endif  !}

add_name = .true.

if (present(index)) then  !{
  if (fm_exists(name)) then  !{
    field_length = fm_get_length(name)
    if (field_length .lt. 0) then  !{
      call mpp_error(FATAL, trim(error_header) // ' Problem getting length of ' // trim(name))
    endif  !}
    if (.not. (no_overwrite_use .and. field_length .ge. index)) then  !{
      field_index = fm_new_value(name, value, index = index)
      if (field_index .le. 0) then  !{
        write (str_error,*) ' with index = ', index
        call mpp_error(FATAL, trim(error_header) // ' Problem overwriting ' // trim(name) // trim(str_error))
      endif  !}
    endif  !}
  else  !}{
    field_index = fm_new_value(name, value, index = index)
    if (field_index .le. 0) then  !{
      write (str_error,*) ' with index = ', index
      call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name) // trim(str_error))
    endif  !}
  endif  !}
elseif (present(append)) then  !}{
  field_index = fm_new_value(name, value, append = append)
  if (field_index .le. 0) then  !{
    write (str_error,*) ' with append = ', append
    call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name) // trim(str_error))
  endif  !}
else  !}{
  if (fm_exists(name)) then  !{
    if (.not. no_overwrite_use) then  !{
      field_index = fm_new_value(name, value)
      if (field_index .le. 0) then  !{
        call mpp_error(FATAL, trim(error_header) // ' Problem overwriting ' // trim(name))
      endif  !}
    endif  !}
  elseif (create) then  !}{
    field_index = fm_new_value(name, value)
    if (field_index .le. 0) then  !{
      call mpp_error(FATAL, trim(error_header) // ' Problem creating ' // trim(name))
    endif  !}
  else  !}{
    add_name = .false.
  endif  !}
endif  !}

!
!       Add the variable name to the list of good names, to be used
!       later for a consistency check, unless the field did not exist and we did not create it
!

if (good_name_list_use .ne. ' ' .and. add_name) then  !{
  if (fm_new_value(good_name_list_use, name, append = .true., create = .true.) .le. 0) then  !{
    call mpp_error(FATAL, trim(error_header) //                           &
         ' Could not add ' // trim(name) // ' to "' // trim(good_name_list_use) // '" list')
  endif  !}
endif  !}

return

end subroutine otpm_set_value_integer  !}
! </SUBROUTINE> NAME="otpm_set_value_integer"


!#######################################################################
! <SUBROUTINE NAME="otpm_set_value_logical">
!
! <DESCRIPTION>
! Set a logical value in the Field Manager tree.
! </DESCRIPTION>
!

subroutine otpm_set_value_logical(name, value, caller, index, append, no_create,        &
     no_overwrite, good_name_list)  !{

implicit none

!
!       arguments
!

character(len=*), intent(in)            :: name
logical, intent(in)                     :: value
character(len=*), intent(in), optional  :: caller
integer, intent(in), optional           :: index
logical, intent(in), optional           :: append
logical, intent(in), optional           :: no_create
logical, intent(in), optional           :: no_overwrite
character(len=*), intent(in), optional  :: good_name_list

!
!       Local parameters
!

character(len=48), parameter    :: sub_name = 'otpm_set_value_logical'

!
!       Local variables
!

character(len=256)              :: error_header
character(len=256)              :: warn_header
character(len=256)              :: note_header
character(len=128)              :: caller_str
character(len=32)               :: str_error
integer                         :: field_index
logical                         :: no_overwrite_use
integer                         :: field_length
character(len=fm_path_name_len) :: good_name_list_use
logical                         :: create
logical                         :: add_name

!
!       set the caller string and headers
!

if (present(caller)) then  !{
  caller_str = '[' // trim(caller) // ']'
else  !}{
  caller_str = default_caller
endif  !}

error_header = '==>Error from ' // trim(mod_name) //   &
               '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

!
!       check that append and index are not both given
!

if (present(index) .and. present(append)) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Append and index both given as arguments')
endif  !}

!
!       check for whether to overwrite existing values
!

if (present(no_overwrite)) then  !{
  no_overwrite_use = no_overwrite
else  !}{
  no_overwrite_use = default_no_overwrite
endif  !}

!
!       check for whether to save the name in a list
!

if (present(good_name_list)) then  !{
  good_name_list_use = good_name_list
else  !}{
  good_name_list_use = default_good_name_list
endif  !}

if (present(no_create)) then  !{
  create = .not. no_create
  if (no_create .and. (present(append) .or. present(index))) then  !{
    call mpp_error(FATAL, trim(error_header) // ' append or index are present when no_create is true for ' // trim(name))
  endif  !}
else  !}{
  create = .true.
endif  !}

add_name = .true.

if (present(index)) then  !{
  if (fm_exists(name)) then  !{
    field_length = fm_get_length(name)
    if (field_length .lt. 0) then  !{
      call mpp_error(FATAL, trim(error_header) // ' Problem getting length of ' // trim(name))
    endif  !}
    if (.not. (no_overwrite_use .and. field_length .ge. index)) then  !{
      field_index = fm_new_value(name, value, index = index)
      if (field_index .le. 0) then  !{
        write (str_error,*) ' with index = ', index
        call mpp_error(FATAL, trim(error_header) // ' Problem overwriting ' // trim(name) // trim(str_error))
      endif  !}
    endif  !}
  else  !}{
    field_index = fm_new_value(name, value, index = index)
    if (field_index .le. 0) then  !{
      write (str_error,*) ' with index = ', index
      call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name) // trim(str_error))
    endif  !}
  endif  !}
elseif (present(append)) then  !}{
  field_index = fm_new_value(name, value, append = append)
  if (field_index .le. 0) then  !{
    write (str_error,*) ' with append = ', append
    call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name) // trim(str_error))
  endif  !}
else  !}{
  if (fm_exists(name)) then  !{
    if (.not. no_overwrite_use) then  !{
      field_index = fm_new_value(name, value)
      if (field_index .le. 0) then  !{
        call mpp_error(FATAL, trim(error_header) // ' Problem overwriting ' // trim(name))
      endif  !}
    endif  !}
  elseif (create) then  !}{
    field_index = fm_new_value(name, value)
    if (field_index .le. 0) then  !{
      call mpp_error(FATAL, trim(error_header) // ' Problem creating ' // trim(name))
    endif  !}
  else  !}{
    add_name = .false.
  endif  !}
endif  !}

!
!       Add the variable name to the list of good names, to be used
!       later for a consistency check, unless the field did not exist and we did not create it
!

if (good_name_list_use .ne. ' ' .and. add_name) then  !{
  if (fm_new_value(good_name_list_use, name, append = .true., create = .true.) .le. 0) then  !{
    call mpp_error(FATAL, trim(error_header) //                           &
         ' Could not add ' // trim(name) // ' to "' // trim(good_name_list_use) // '" list')
  endif  !}
endif  !}

return

end subroutine otpm_set_value_logical  !}
! </SUBROUTINE> NAME="otpm_set_value_logical"


!#######################################################################
! <SUBROUTINE NAME="otpm_set_value_real">
!
! <DESCRIPTION>
! Set a real value in the Field Manager tree.
! </DESCRIPTION>
!

subroutine otpm_set_value_real(name, value, caller, index, append, no_create,        &
     no_overwrite, good_name_list)  !{

implicit none

!
!       arguments
!

character(len=*), intent(in)            :: name
real, intent(in)                        :: value
character(len=*), intent(in), optional  :: caller
integer, intent(in), optional           :: index
logical, intent(in), optional           :: append
logical, intent(in), optional           :: no_create
logical, intent(in), optional           :: no_overwrite
character(len=*), intent(in), optional  :: good_name_list

!
!       Local parameters
!

character(len=48), parameter    :: sub_name = 'otpm_set_value_real'

!
!       Local variables
!

character(len=256)              :: error_header
character(len=256)              :: warn_header
character(len=256)              :: note_header
character(len=128)              :: caller_str
character(len=32)               :: str_error
integer                         :: field_index
logical                         :: no_overwrite_use
integer                         :: field_length
character(len=fm_path_name_len) :: good_name_list_use
logical                         :: create
logical                         :: add_name

!
!       set the caller string and headers
!

if (present(caller)) then  !{
  caller_str = '[' // trim(caller) // ']'
else  !}{
  caller_str = default_caller
endif  !}

error_header = '==>Error from ' // trim(mod_name) //   &
               '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

!
!       check that append and index are not both given
!

if (present(index) .and. present(append)) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Append and index both given as arguments')
endif  !}

!
!       check for whether to overwrite existing values
!

if (present(no_overwrite)) then  !{
  no_overwrite_use = no_overwrite
else  !}{
  no_overwrite_use = default_no_overwrite
endif  !}

!
!       check for whether to save the name in a list
!

if (present(good_name_list)) then  !{
  good_name_list_use = good_name_list
else  !}{
  good_name_list_use = default_good_name_list
endif  !}

if (present(no_create)) then  !{
  create = .not. no_create
  if (no_create .and. (present(append) .or. present(index))) then  !{
    call mpp_error(FATAL, trim(error_header) // ' append or index are present when no_create is true for ' // trim(name))
  endif  !}
else  !}{
  create = .true.
endif  !}

add_name = .true.

if (present(index)) then  !{
  if (fm_exists(name)) then  !{
    field_length = fm_get_length(name)
    if (field_length .lt. 0) then  !{
      call mpp_error(FATAL, trim(error_header) // ' Problem getting length of ' // trim(name))
    endif  !}
    if (.not. (no_overwrite_use .and. field_length .ge. index)) then  !{
      field_index = fm_new_value(name, value, index = index)
      if (field_index .le. 0) then  !{
        write (str_error,*) ' with index = ', index
        call mpp_error(FATAL, trim(error_header) // ' Problem overwriting ' // trim(name) // trim(str_error))
      endif  !}
    endif  !}
  else  !}{
    field_index = fm_new_value(name, value, index = index)
    if (field_index .le. 0) then  !{
      write (str_error,*) ' with index = ', index
      call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name) // trim(str_error))
    endif  !}
  endif  !}
elseif (present(append)) then  !}{
  field_index = fm_new_value(name, value, append = append)
  if (field_index .le. 0) then  !{
    write (str_error,*) ' with append = ', append
    call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name) // trim(str_error))
  endif  !}
else  !}{
  if (fm_exists(name)) then  !{
    if (.not. no_overwrite_use) then  !{
      field_index = fm_new_value(name, value)
      if (field_index .le. 0) then  !{
        call mpp_error(FATAL, trim(error_header) // ' Problem overwriting ' // trim(name))
      endif  !}
    endif  !}
  elseif (create) then  !}{
    field_index = fm_new_value(name, value)
    if (field_index .le. 0) then  !{
      call mpp_error(FATAL, trim(error_header) // ' Problem creating ' // trim(name))
    endif  !}
  else  !}{
    add_name = .false.
  endif  !}
endif  !}

!
!       Add the variable name to the list of good names, to be used
!       later for a consistency check, unless the field did not exist and we did not create it
!

if (good_name_list_use .ne. ' ' .and. add_name) then  !{
  if (fm_new_value(good_name_list_use, name, append = .true., create = .true.) .le. 0) then  !{
    call mpp_error(FATAL, trim(error_header) //                           &
         ' Could not add ' // trim(name) // ' to "' // trim(good_name_list_use) // '" list')
  endif  !}
endif  !}

return

end subroutine otpm_set_value_real  !}
! </SUBROUTINE> NAME="otpm_set_value_real"


!#######################################################################
! <SUBROUTINE NAME="otpm_set_value_string">
!
! <DESCRIPTION>
! Set a string value in the Field Manager tree.
! </DESCRIPTION>
!

subroutine otpm_set_value_string(name, value, caller, index, append, no_create,        &
     no_overwrite, good_name_list)  !{

implicit none

!
!       arguments
!

character(len=*), intent(in)            :: name
character(len=*), intent(in)            :: value
character(len=*), intent(in), optional  :: caller
integer, intent(in), optional           :: index
logical, intent(in), optional           :: append
logical, intent(in), optional           :: no_create
logical, intent(in), optional           :: no_overwrite
character(len=*), intent(in), optional  :: good_name_list

!
!       Local parameters
!

character(len=48), parameter    :: sub_name = 'otpm_set_value_string'

!
!       Local variables
!

character(len=256)              :: error_header
character(len=256)              :: warn_header
character(len=256)              :: note_header
character(len=128)              :: caller_str
character(len=32)               :: str_error
integer                         :: field_index
logical                         :: no_overwrite_use
integer                         :: field_length
character(len=fm_path_name_len) :: good_name_list_use
logical                         :: create
logical                         :: add_name

!
!       set the caller string and headers
!

if (present(caller)) then  !{
  caller_str = '[' // trim(caller) // ']'
else  !}{
  caller_str = default_caller
endif  !}

error_header = '==>Error from ' // trim(mod_name) //   &
               '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

!
!       check that append and index are not both given
!

if (present(index) .and. present(append)) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Append and index both given as arguments')
endif  !}

!
!       check for whether to overwrite existing values
!

if (present(no_overwrite)) then  !{
  no_overwrite_use = no_overwrite
else  !}{
  no_overwrite_use = default_no_overwrite
endif  !}

!
!       check for whether to save the name in a list
!

if (present(good_name_list)) then  !{
  good_name_list_use = good_name_list
else  !}{
  good_name_list_use = default_good_name_list
endif  !}

if (present(no_create)) then  !{
  create = .not. no_create
  if (no_create .and. (present(append) .or. present(index))) then  !{
    call mpp_error(FATAL, trim(error_header) // ' append or index are present when no_create is true for ' // trim(name))
  endif  !}
else  !}{
  create = .true.
endif  !}

add_name = .true.

if (present(index)) then  !{
  if (fm_exists(name)) then  !{
    field_length = fm_get_length(name)
    if (field_length .lt. 0) then  !{
      call mpp_error(FATAL, trim(error_header) // ' Problem getting length of ' // trim(name))
    endif  !}
    if (.not. (no_overwrite_use .and. field_length .ge. index)) then  !{
      field_index = fm_new_value(name, value, index = index)
      if (field_index .le. 0) then  !{
        write (str_error,*) ' with index = ', index
        call mpp_error(FATAL, trim(error_header) // ' Problem overwriting ' // trim(name) // trim(str_error))
      endif  !}
    endif  !}
  else  !}{
    field_index = fm_new_value(name, value, index = index)
    if (field_index .le. 0) then  !{
      write (str_error,*) ' with index = ', index
      call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name) // trim(str_error))
    endif  !}
  endif  !}
elseif (present(append)) then  !}{
  field_index = fm_new_value(name, value, append = append)
  if (field_index .le. 0) then  !{
    write (str_error,*) ' with append = ', append
    call mpp_error(FATAL, trim(error_header) // ' Problem setting ' // trim(name) // trim(str_error))
  endif  !}
else  !}{
  if (fm_exists(name)) then  !{
    if (.not. no_overwrite_use) then  !{
      field_index = fm_new_value(name, value)
      if (field_index .le. 0) then  !{
        call mpp_error(FATAL, trim(error_header) // ' Problem overwriting ' // trim(name))
      endif  !}
    endif  !}
  elseif (create) then  !}{
    field_index = fm_new_value(name, value)
    if (field_index .le. 0) then  !{
      call mpp_error(FATAL, trim(error_header) // ' Problem creating ' // trim(name))
    endif  !}
  else  !}{
    add_name = .false.
  endif  !}
endif  !}

!
!       Add the variable name to the list of good names, to be used
!       later for a consistency check, unless the field did not exist and we did not create it
!

if (good_name_list_use .ne. ' ' .and. add_name) then  !{
  if (fm_new_value(good_name_list_use, name, append = .true., create = .true.) .le. 0) then  !{
    call mpp_error(FATAL, trim(error_header) //                           &
         ' Could not add ' // trim(name) // ' to "' // trim(good_name_list_use) // '" list')
  endif  !}
endif  !}

return

end subroutine otpm_set_value_string  !}
! </SUBROUTINE> NAME="otpm_set_value_string"


!#######################################################################
! <SUBROUTINE NAME="set_prog_value_integer">
!
! <DESCRIPTION>
! Set an integer value for a prognostic tracer element in the Field Manager tree.
! </DESCRIPTION>
!

subroutine set_prog_value_integer(name, package_name, value, caller)  !{

implicit none

!
!       arguments
!

character(len=*), intent(in)            :: name
character(len=*), intent(in)            :: package_name
integer, intent(in)                     :: value
character(len=*), intent(in), optional  :: caller

!
!       Local parameters
!

character(len=48), parameter    :: sub_name = 'set_prog_value_integer'

!
!       Local variables
!

integer                         :: integer_value
character(len=256)              :: error_header
character(len=256)              :: warn_header
character(len=256)              :: note_header
character(len=128)              :: caller_str
integer                         :: length

!
!       set the caller string and headers
!

if (present(caller)) then  !{
  caller_str = '[' // trim(caller) // ']'
else  !}{
  caller_str = default_caller
endif  !}

error_header = '==>Error from ' // trim(mod_name) //   &
               '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

!
!       check that a package name is given (fatal if not)
!

if (package_name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty package_name given')
endif  !}

!
!       The following steps are done when setting elements in the prognostic
!       tracer. Note that the subroutine otpm_set_value should have been
!       set so as to only set the given value if there is not already a
!       value in the tracer tree (such as one specified from the field table)
!
!       The precedence of values to use is as follows:
!               1) a value set in the field table
!               2) a value present in the package defaults (either
!                       from the field table or otpm_set_tracer_package)
!               3) the value given to this subroutine
!

!       
!       First check whether there is a package default set
!

length = fm_get_length(trim(package_name) // name)
if (length .gt. 1) then  !{

  !     Error: package default is not a scalar

  call mpp_error(FATAL, trim(error_header) //                         &
       ' "' // trim(name) // '" not a scalar in: ' // trim(package_name))

elseif (length .le. 0) then  !}{

  ! Package default does not exist or is null, so use, in order:
  ! 1) value specified in field table (implicit in the otpm_set_value call)
  ! 2) value given in argument list

  call otpm_set_value(name, value)

else  !}{

  ! Package default exists, so use, in order:
  ! 1) value specified in field table (implicit in the otpm_set_value call)
  ! 2) package default

  if (fm_get_value(trim(package_name) // name, integer_value)) then  !{

    call otpm_set_value(name, integer_value)

  else  !}{

    ! This error shouldn't happen since the previous calls all
    ! show that the value exists, unless, perhaps, the type is incorrect

    call mpp_error(FATAL, trim(error_header) //                         &
         ' Could not get "' // trim(name) // '" from: ' // trim(package_name))

  endif  !}

endif  !}

return

end subroutine set_prog_value_integer  !}
! </SUBROUTINE> NAME="set_prog_value_integer"


!#######################################################################
! <SUBROUTINE NAME="set_prog_value_logical">
!
! <DESCRIPTION>
! Set a logical value for a prognostic tracer element in the Field Manager tree.
! </DESCRIPTION>
!

subroutine set_prog_value_logical(name, package_name, value, caller)  !{

implicit none

!
!       arguments
!

character(len=*), intent(in)            :: name
character(len=*), intent(in)            :: package_name
logical, intent(in)                     :: value
character(len=*), intent(in), optional  :: caller

!
!       Local parameters
!

character(len=48), parameter    :: sub_name = 'set_prog_value_logical'

!
!       Local variables
!

logical                         :: logical_value
character(len=256)              :: error_header
character(len=256)              :: warn_header
character(len=256)              :: note_header
character(len=128)              :: caller_str
integer                         :: length

!
!       set the caller string and headers
!

if (present(caller)) then  !{
  caller_str = '[' // trim(caller) // ']'
else  !}{
  caller_str = default_caller
endif  !}

error_header = '==>Error from ' // trim(mod_name) //   &
               '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

!
!       check that a package name is given (fatal if not)
!

if (package_name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty package_name given')
endif  !}

!
!       The following steps are done when setting elements in the prognostic
!       tracer. Note that the subroutine otpm_set_value should have been
!       set so as to only set the given value if there is not already a
!       value in the tracer tree (such as one specified from the field table)
!
!       The precedence of values to use is as follows:
!               1) a value set in the field table
!               2) a value present in the package defaults (either
!                       from the field table or otpm_set_tracer_package)
!               3) the value given to this subroutine
!

!       
!       First check whether there is a package default set
!

length = fm_get_length(trim(package_name) // name)
if (length .gt. 1) then  !{

  !     Error: package default is not a scalar

  call mpp_error(FATAL, trim(error_header) //                         &
       ' "' // trim(name) // '" not a scalar in: ' // trim(package_name))

elseif (length .le. 0) then  !}{

  ! Package default does not exist or is null, so use, in order:
  ! 1) value specified in field table (implicit in the otpm_set_value call)
  ! 2) value given in argument list

  call otpm_set_value(name, value)

else  !}{

  ! Package default exists, so use, in order:
  ! 1) value specified in field table (implicit in the otpm_set_value call)
  ! 2) package default

  if (fm_get_value(trim(package_name) // name, logical_value)) then  !{

    call otpm_set_value(name, logical_value)

  else  !}{

    ! This error shouldn't happen since the previous calls all
    ! show that the value exists, unless, perhaps, the type is incorrect

    call mpp_error(FATAL, trim(error_header) //                         &
         ' Could not get "' // trim(name) // '" from: ' // trim(package_name))

  endif  !}

endif  !}

return

end subroutine set_prog_value_logical  !}
! </SUBROUTINE> NAME="set_prog_value_logical"


!#######################################################################
! <SUBROUTINE NAME="set_prog_value_real">
!
! <DESCRIPTION>
! Set a real value for a prognostic tracer element in the Field Manager tree.
! </DESCRIPTION>
!

subroutine set_prog_value_real(name, package_name, value, caller)  !{

implicit none

!
!       arguments
!

character(len=*), intent(in)            :: name
character(len=*), intent(in)            :: package_name
real, intent(in)                        :: value
character(len=*), intent(in), optional  :: caller

!
!       Local parameters
!

character(len=48), parameter    :: sub_name = 'set_prog_value_real'

!
!       Local variables
!

real                            :: real_value
character(len=256)              :: error_header
character(len=256)              :: warn_header
character(len=256)              :: note_header
character(len=128)              :: caller_str
integer                         :: length

!
!       set the caller string and headers
!

if (present(caller)) then  !{
  caller_str = '[' // trim(caller) // ']'
else  !}{
  caller_str = default_caller
endif  !}

error_header = '==>Error from ' // trim(mod_name) //   &
               '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

!
!       check that a package name is given (fatal if not)
!

if (package_name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty package_name given')
endif  !}

!
!       The following steps are done when setting elements in the prognostic
!       tracer. Note that the subroutine otpm_set_value should have been
!       set so as to only set the given value if there is not already a
!       value in the tracer tree (such as one specified from the field table)
!
!       The precedence of values to use is as follows:
!               1) a value set in the field table
!               2) a value present in the package defaults (either
!                       from the field table or otpm_set_tracer_package)
!               3) the value given to this subroutine
!

!       
!       First check whether there is a package default set
!

length = fm_get_length(trim(package_name) // name)
if (length .gt. 1) then  !{

  !     Error: package default is not a scalar

  call mpp_error(FATAL, trim(error_header) //                         &
       ' "' // trim(name) // '" not a scalar in: ' // trim(package_name))

elseif (length .le. 0) then  !}{

  ! Package default does not exist or is null, so use, in order:
  ! 1) value specified in field table (implicit in the otpm_set_value call)
  ! 2) value given in argument list

  call otpm_set_value(name, value)

else  !}{

  ! Package default exists, so use, in order:
  ! 1) value specified in field table (implicit in the otpm_set_value call)
  ! 2) package default

  if (fm_get_value(trim(package_name) // name, real_value)) then  !{

    call otpm_set_value(name, real_value)

  else  !}{

    ! This error shouldn't happen since the previous calls all
    ! show that the value exists, unless, perhaps, the type is incorrect

    call mpp_error(FATAL, trim(error_header) //                         &
         ' Could not get "' // trim(name) // '" from: ' // trim(package_name))

  endif  !}

endif  !}

return

end subroutine set_prog_value_real  !}
! </SUBROUTINE> NAME="set_prog_value_real"


!#######################################################################
! <SUBROUTINE NAME="set_prog_value_string">
!
! <DESCRIPTION>
! Set a string value for a prognostic tracer element in the Field Manager tree.
! </DESCRIPTION>
!

subroutine set_prog_value_string(name, package_name, value, caller)  !{

implicit none

!
!       arguments
!

character(len=*), intent(in)            :: name
character(len=*), intent(in)            :: package_name
character(len=*), intent(in)            :: value
character(len=*), intent(in), optional  :: caller

!
!       Local parameters
!

character(len=48), parameter    :: sub_name = 'set_prog_value_string'

!
!       Local variables
!

character(len=fm_string_len)    :: string_value
character(len=256)              :: error_header
character(len=256)              :: warn_header
character(len=256)              :: note_header
character(len=128)              :: caller_str
integer                         :: length

!
!       set the caller string and headers
!

if (present(caller)) then  !{
  caller_str = '[' // trim(caller) // ']'
else  !}{
  caller_str = default_caller
endif  !}

error_header = '==>Error from ' // trim(mod_name) //   &
               '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

!
!       check that a package name is given (fatal if not)
!

if (package_name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty package_name given')
endif  !}

!
!       The following steps are done when setting elements in the prognostic
!       tracer. Note that the subroutine otpm_set_value should have been
!       set so as to only set the given value if there is not already a
!       value in the tracer tree (such as one specified from the field table)
!
!       The precedence of values to use is as follows:
!               1) a value set in the field table
!               2) a value present in the package defaults (either
!                       from the field table or otpm_set_tracer_package)
!               3) the value given to this subroutine
!

!       
!       First check whether there is a package default set
!

length = fm_get_length(trim(package_name) // name)
if (length .gt. 1) then  !{

  !     Error: package default is not a scalar

  call mpp_error(FATAL, trim(error_header) //                         &
       ' "' // trim(name) // '" not a scalar in: ' // trim(package_name))

elseif (length .le. 0) then  !}{

  ! Package default does not exist or is null, so use, in order:
  ! 1) value specified in field table (implicit in the otpm_set_value call)
  ! 2) value given in argument list

  call otpm_set_value(name, value)

else  !}{

  ! Package default exists, so use, in order:
  ! 1) value specified in field table (implicit in the otpm_set_value call)
  ! 2) package default

  if (fm_get_value(trim(package_name) // name, string_value)) then  !{

    call otpm_set_value(name, string_value)

  else  !}{

    ! This error shouldn't happen since the previous calls all
    ! show that the value exists, unless, perhaps, the type is incorrect

    call mpp_error(FATAL, trim(error_header) //                         &
         ' Could not get "' // trim(name) // '" from: ' // trim(package_name))

  endif  !}

endif  !}

return

end subroutine set_prog_value_string  !}
! </SUBROUTINE> NAME="set_prog_value_string"


!#######################################################################
! <FUNCTION NAME="otpm_set_prog_tracer">
!
! <DESCRIPTION>
! Set the values for a prog tracer and return its index (0 on error)
! </DESCRIPTION>
!
function otpm_set_prog_tracer(name, package, caller, longname, units,   &
     conversion, offset,                                                &
     min_tracer, max_tracer, use_only_advection,                        &
     min_range, max_range, file_in, name_in, scale_in,                  &
     additive_in, init,                                                 &
     const_init_tracer, const_init_value, file_out, flux_units,         &
     min_flux_range, max_flux_range,                                    &
     min_tracer_limit, max_tracer_limit,                                &
     vert_adv_scheme, horiz_adv_scheme)                                 &
         result (prog_index)  !{

implicit none

!
!       Return type
!

integer :: prog_index

!
!       arguments
!

character(len=*), intent(in)            :: name
character(len=*), intent(in)            :: package
character(len=*), intent(in), optional  :: caller
character(len=*), intent(in), optional  :: units
character(len=*), intent(in), optional  :: flux_units
character(len=*), intent(in), optional  :: longname
character(len=*), intent(in), optional  :: vert_adv_scheme
character(len=*), intent(in), optional  :: horiz_adv_scheme
real, intent(in), optional              :: conversion
character(len=*), intent(in), optional  :: file_in
character(len=*), intent(in), optional  :: file_out
logical, intent(in), optional           :: init
logical, intent(in), optional           :: use_only_advection
real, intent(in), optional              :: max_range
real, intent(in), optional              :: min_range
real, intent(in), optional              :: max_flux_range
real, intent(in), optional              :: min_flux_range
real, intent(in), optional              :: max_tracer
real, intent(in), optional              :: min_tracer
real, intent(in), optional              :: max_tracer_limit
real, intent(in), optional              :: min_tracer_limit
character(len=*), intent(in), optional  :: name_in
real, intent(in), optional              :: offset
real, intent(in), optional              :: scale_in
real, intent(in), optional              :: additive_in
logical, intent(in), optional           :: const_init_tracer 
real, intent(in), optional              :: const_init_value  

!
!       Local parameters
!

character(len=48), parameter  :: sub_name = 'otpm_set_prog_tracer'

!
!       Local variables
!

character(len=fm_path_name_len)                         :: current_list
character(len=fm_path_name_len)                         :: package_name
character(len=256)                                      :: error_header
character(len=256)                                      :: warn_header
character(len=256)                                      :: note_header
character(len=128)                                      :: caller_str
character(len=fm_path_name_len)                         :: tracer_name
character(len=fm_string_len), pointer, dimension(:)     :: good_list => NULL()

!
!       set the caller string and headers
!

if (present(caller)) then  !{
  caller_str = '[' // trim(caller) // ']'
else  !}{
  caller_str = default_caller
endif  !}

error_header = '==>Error from ' // trim(mod_name) //   &
               '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

!
!       check the package name
!

if (package .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty package given')
endif  !}
package_name = '/ocean_mod/tracer_packages/' // trim(package) // '/'
if (fm_get_type(package_name) .ne. 'list') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Package does not exist or is not a list: ' // trim(package))
endif  !}

!
!       Begin processing
!

write (stdout(),*)
write (stdout(),*) trim(note_header), ' Processing prog tracer ', trim(name)

!
!       Check whether the tracer already exists. If so, then use that tracer
!

tracer_name = '/ocean_mod/prog_tracers/' // trim(name) // '/'
prog_index = fm_get_index(tracer_name)

if (prog_index .le. 0) then  !{

!
!       Set a new prog tracer and get its index
!

  prog_index = fm_new_list(tracer_name)
  if (prog_index .le. 0) then  !{
    call mpp_error(FATAL, trim(error_header) // ' Could not set prog tracer ' // trim(name))
  endif  !}

endif  !}

!
!       Change to the new tracer, first saving the current list
!

current_list = fm_get_current_list()
if (current_list .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not get the current list')
endif  !}

if (.not. fm_change_list(tracer_name)) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not change to the new list')
endif  !}

!
!       Set the array in which to save the valid names for this list,
!       used later for a consistency check. This is used in the otpm_set_value
!       routines to make the list of valid values
!

call otpm_set_good_name_list('/ocean_mod/GOOD/prog_tracers/' // trim(name) // '/good_list')

!
!       Set other defaults for the otpm_set_value routines
!

call otpm_set_caller(caller_str)

!
!       When the following is set to true, otpm_set_value will not overwrite
!       any values already set in the tracer tree If there is no
!       value present, then a new entry will be created in the tracer tree and
!       the value supplied will be set.
!

call otpm_set_no_overwrite(.true.)

!
!       Set various values to given values, or to defaults if not given
!

!
!       The longname is distinct here in that there is no option for a package
!       default. Hence, the precedence of values is:
!               1) a value set in the field table
!               2) an optional argument given to this subroutine
!               3) the tracer name
!

if (present(longname)) then  !{
  call otpm_set_value('longname', longname)
else  !}{
  call otpm_set_value('longname', name)
endif  !}

!
!       The precedence of values to use in set_prog_value is as follows:
!               1) a value set in the field table
!               2) a value present in the package defaults (either
!                       from the field table or otpm_set_tracer_package)
!               3) the value passed to it in the argument list
!       This subroutine will preferentially supply the given optional
!       argument over the module-wide default value
!

if (present(units)) then  !{
  call set_prog_value('units', package_name, units)
else  !}{
  call set_prog_value('units', package_name, default_units)
endif  !}

if (present(conversion)) then  !}{
  call set_prog_value('conversion', package_name, conversion)
else  !}{
  call set_prog_value('conversion', package_name, default_conversion)
endif  !}

if (present(offset)) then  !}{
  call set_prog_value('offset', package_name, offset)
else  !}{
  call set_prog_value('offset', package_name, default_offset)
endif  !}

if (present(min_tracer)) then  !}{
  call set_prog_value('min_tracer', package_name, min_tracer)
else  !}{
  call set_prog_value('min_tracer', package_name, default_min_tracer)
endif  !}

if (present(max_tracer)) then  !}{
  call set_prog_value('max_tracer', package_name, max_tracer)
else  !}{
  call set_prog_value('max_tracer', package_name, default_max_tracer)
endif  !}

if (present(min_range)) then  !}{
  call set_prog_value('min_range', package_name, min_range)
else  !}{
  call set_prog_value('min_range', package_name, default_min_range)
endif  !}

if (present(max_range)) then  !}{
  call set_prog_value('max_range', package_name, max_range)
else  !}{
  call set_prog_value('max_range', package_name, default_max_range)
endif  !}

if (present(use_only_advection)) then  !}{
  call set_prog_value('use_only_advection', package_name, use_only_advection)
else  !}{
  call set_prog_value('use_only_advection', package_name, default_use_only_advection)
endif  !}

if (present(file_in)) then  !}{
  call set_prog_value('file_in', package_name, file_in)
else  !}{
  call set_prog_value('file_in', package_name, default_file_in)
endif  !}

if (present(init)) then  !}{
  call set_prog_value('init', package_name, init)
else  !}{
  call set_prog_value('init', package_name, default_init)
endif  !}

if (present(name_in)) then  !}{
  call set_prog_value('name_in', package_name, name_in)
else  !}{
  call set_prog_value('name_in', package_name, name)
endif  !}

if (present(scale_in)) then  !}{
  call set_prog_value('scale_in', package_name, scale_in)
else  !}{
  call set_prog_value('scale_in', package_name, default_scale_in)
endif  !}

if (present(additive_in)) then  !}{
  call set_prog_value('additive_in', package_name, additive_in)
else  !}{
  call set_prog_value('additive_in', package_name, default_additive_in)
endif  !}

if (present(const_init_tracer)) then  !}{
  call set_prog_value('const_init_tracer', package_name, const_init_tracer)
else  !}{
  call set_prog_value('const_init_tracer', package_name, default_const_init_tracer)
endif  !}

if (present(const_init_value)) then  !}{
  call set_prog_value('const_init_value', package_name, const_init_value)
else  !}{
  call set_prog_value('const_init_value', package_name, default_const_init_value)
endif  !}

if (present(file_out)) then  !}{
  call set_prog_value('file_out', package_name, file_out)
else  !}{
  call set_prog_value('file_out', package_name, default_file_out)
endif  !}

if (present(flux_units)) then  !}{
  call set_prog_value('flux_units', package_name, flux_units)
else  !}{
  call set_prog_value('flux_units', package_name, default_flux_units)
endif  !}

if (present(min_flux_range)) then  !}{
  call set_prog_value('min_flux_range', package_name, min_flux_range)
else  !}{
  call set_prog_value('min_flux_range', package_name, default_min_flux_range)
endif  !}

if (present(max_flux_range)) then  !}{
  call set_prog_value('max_flux_range', package_name, max_flux_range)
else  !}{
  call set_prog_value('max_flux_range', package_name, default_max_flux_range)
endif  !}

if (present(min_tracer_limit)) then  !}{
  call set_prog_value('min_tracer_limit', package_name, min_tracer_limit)
else  !}{
  call set_prog_value('min_tracer_limit', package_name, default_min_tracer_limit)
endif  !}

if (present(max_tracer_limit)) then  !}{
  call set_prog_value('max_tracer_limit', package_name, max_tracer_limit)
else  !}{
  call set_prog_value('max_tracer_limit', package_name, default_max_tracer_limit)
endif  !}

if (present(vert_adv_scheme)) then  !}{
  call set_prog_value('vertical-advection-scheme', package_name, vert_adv_scheme)
else  !}{
  call set_prog_value('vertical-advection-scheme', package_name, default_vert_adv_scheme)
endif  !}

if (present(horiz_adv_scheme)) then  !}{
  call set_prog_value('horizontal-advection-scheme', package_name, horiz_adv_scheme)
else  !}{
  call set_prog_value('horizontal-advection-scheme', package_name, default_horiz_adv_scheme)
endif  !}

!
!       Reset the defaults for the otpm_set_value calls
!

call otpm_reset_good_name_list
call otpm_reset_no_overwrite
call otpm_reset_caller

!
!       Change back to the saved current list
!

if (.not. fm_change_list(current_list)) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not change back to ' // trim(current_list))
endif  !}

!
!       Check for any errors in the number of fields in this list
!

if (caller_str .eq. ' ') then  !{
  caller_str = trim(mod_name) // '(' // trim(sub_name) // ')'
endif  !}
good_list => otpm_get_string_array('/ocean_mod/GOOD/prog_tracers/' // trim(name) // '/good_list', &
     caller = caller_str)
if (associated(good_list)) then  !{
  call otpm_check_for_bad_fields('/ocean_mod/prog_tracers/' // trim(name), good_list, caller = caller_str)
  deallocate(good_list)
else  !}{
  call mpp_error(FATAL, trim(error_header) // ' Empty "' // trim(name) // '" list')
endif  !}

!
!       Add the tracer name to the list of good tracers, to be used
!       later for a consistency check
!

if (fm_new_value('/ocean_mod/GOOD/good_prog_tracers', name, append = .true.) .le. 0) then  !{
  call mpp_error(FATAL, trim(error_header) //                           &
       ' Could not add ' // trim(name) // ' to "good_prog_tracers" list')
endif  !}

return

end function otpm_set_prog_tracer  !}
! </FUNCTION> NAME="otpm_set_prog_tracer"


!#######################################################################
! <FUNCTION NAME="otpm_set_diag_tracer">
!
! <DESCRIPTION>
! Set the values for a diag tracer and return its index (0 on error)
! </DESCRIPTION>
!
function otpm_set_diag_tracer(name, caller, longname, units,    &
     conversion, offset, min_tracer, max_tracer,                &
     min_range, max_range, file_in, name_in, scale_in,          &
     additive_in, init,                                         &
     const_init_tracer, const_init_value, file_out)             &
         result (diag_index)  !{

implicit none

!
!       Return type
!

integer :: diag_index

!
!       arguments
!

character(len=*), intent(in)            :: name
character(len=*), intent(in), optional  :: caller
character(len=*), intent(in), optional  :: units
character(len=*), intent(in), optional  :: longname
real, intent(in), optional              :: conversion
character(len=*), intent(in), optional  :: file_in
character(len=*), intent(in), optional  :: file_out
logical, intent(in), optional           :: init
real, intent(in), optional              :: max_range
real, intent(in), optional              :: min_range
real, intent(in), optional              :: max_tracer
real, intent(in), optional              :: min_tracer
character(len=*), intent(in), optional  :: name_in
real, intent(in), optional              :: offset
real, intent(in), optional              :: scale_in
real, intent(in), optional              :: additive_in
logical, intent(in), optional           :: const_init_tracer 
real, intent(in), optional              :: const_init_value  

!
!       Local parameters
!

character(len=48), parameter  :: sub_name = 'otpm_set_diag_tracer'

!
!       Local variables
!

character(len=fm_path_name_len)                         :: current_list
character(len=256)                                      :: error_header
character(len=256)                                      :: warn_header
character(len=256)                                      :: note_header
character(len=128)                                      :: caller_str
character(len=fm_string_len), pointer, dimension(:)     :: good_list => NULL()

!
!       set the caller string and headers
!

if (present(caller)) then  !{
  caller_str = '[' // trim(caller) // ']'
else  !}{
  caller_str = default_caller
endif  !}

error_header = '==>Error from ' // trim(mod_name) //   &
               '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

write (stdout(),*)
write (stdout(),*) trim(note_header), ' Processing diag tracer ', trim(name)

!
!       Check whether the tracer already exists. If so, then use that one
!

diag_index = fm_get_index('/ocean_mod/diag_tracers/' // trim(name))
if (diag_index .gt. 0) then  !{

  !write (stdout(),*) trim(note_header), ' Diag tracer already set with index ', diag_index

endif  !}

!
!       Set a new diag tracer and get its index
!

diag_index = fm_new_list('/ocean_mod/diag_tracers/' // trim(name))
if (diag_index .le. 0) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not set diag tracer')
endif  !}

!
!       Change to the new tracer, first saving the current list
!

current_list = fm_get_current_list()
if (current_list .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not get the current list')
endif  !}

if (.not. fm_change_list('/ocean_mod/diag_tracers/' // trim(name))) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not change to the new list')
endif  !}

!
!       Set the array in which to save the valid names for this list,
!       used later for a consistency check. This is used in the otpm_set_value
!       routines to make the list of valid values
!

call otpm_set_good_name_list('/ocean_mod/GOOD/diag_tracers/' // trim(name) // '/good_list')

!
!       Set other defaults for the otpm_set_value routines
!

call otpm_set_no_overwrite(.true.)
call otpm_set_caller(caller_str)

!
!       Set various values to given values, or to defaults if not given
!

if (present(longname)) then  !{
  call otpm_set_value('longname', longname)
else  !}{
  call otpm_set_value('longname', name)
endif  !}

if (present(units)) then  !{
  call otpm_set_value('units', units)
else  !}{
  call otpm_set_value('units', default_units)
endif  !}

if (present(conversion)) then  !{
  call otpm_set_value('conversion', conversion)
else  !}{
  call otpm_set_value('conversion', default_conversion)
endif  !}

if (present(offset)) then  !{
  call otpm_set_value('offset', offset)
else  !}{
  call otpm_set_value('offset', default_offset)
endif  !}

if (present(min_tracer)) then  !{
  call otpm_set_value('min_tracer', min_tracer)
else  !}{
  call otpm_set_value('min_tracer', default_min_tracer)
endif  !}

if (present(max_tracer)) then  !{
  call otpm_set_value('max_tracer', max_tracer)
else  !}{
  call otpm_set_value('max_tracer', default_max_tracer)
endif  !}

if (present(min_range)) then  !{
  call otpm_set_value('min_range', min_range)
else  !}{
  call otpm_set_value('min_range', default_min_range)
endif  !}

if (present(max_range)) then  !{
  call otpm_set_value('max_range', max_range)
else  !}{
  call otpm_set_value('max_range', default_max_range)
endif  !}

if (present(file_in)) then  !{
  call otpm_set_value('file_in', file_in)
else  !}{
  call otpm_set_value('file_in', default_file_in)
endif  !}

if (present(init)) then  !{
  call otpm_set_value('init', init)
else  !}{
  call otpm_set_value('init', default_init)
endif  !}

if (present(name_in)) then  !{
  call otpm_set_value('name_in', name_in)
else  !}{
  call otpm_set_value('name_in', name)
endif  !}

if (present(scale_in)) then  !{
  call otpm_set_value('scale_in', scale_in)
else  !}{
  call otpm_set_value('scale_in', default_scale_in)
endif  !}

if (present(additive_in)) then  !{
  call otpm_set_value('additive_in', additive_in)
else  !}{
  call otpm_set_value('additive_in', default_additive_in)
endif  !}

if (present(const_init_tracer)) then  !{
  call otpm_set_value('const_init_tracer', const_init_tracer)
else  !}{
  call otpm_set_value('const_init_tracer', default_const_init_tracer)
endif  !}

if (present(const_init_value)) then  !{
  call otpm_set_value('const_init_value', const_init_value)
else  !}{
  call otpm_set_value('const_init_value', default_const_init_value)
endif  !}

if (present(file_out)) then  !{
  call otpm_set_value('file_out', file_out)
else  !}{
  call otpm_set_value('file_out', default_file_out)
endif  !}

!
!       Reset the defaults for the otpm_set_value calls
!

call otpm_reset_good_name_list
call otpm_reset_no_overwrite
call otpm_reset_caller

!
!       Change back to the saved current list
!

if (.not. fm_change_list(current_list)) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not change back to ' // trim(current_list))
endif  !}

!
!       Check for any errors in the number of fields in this list
!

if (caller_str .eq. ' ') then  !{
  caller_str = trim(mod_name) // '(' // trim(sub_name) // ')'
endif  !}
good_list => otpm_get_string_array('/ocean_mod/GOOD/diag_tracers/' // trim(name) // '/good_list', &
     caller = caller_str)
if (associated(good_list)) then  !{
  call otpm_check_for_bad_fields('/ocean_mod/diag_tracers/' // trim(name), good_list, caller = caller_str)
  deallocate(good_list)
else  !}{
  call mpp_error(FATAL, trim(error_header) // ' Empty "' // trim(name) // '" list')
endif  !}

!
!       Add the tracer name to the list of good tracers, to be used
!       later for a consistency check
!

if (fm_new_value('/ocean_mod/GOOD/good_diag_tracers', name, append = .true.) .le. 0) then  !{
  call mpp_error(FATAL, trim(error_header) //                           &
       ' Could not add ' // trim(name) // ' to "good_diag_tracers" list')
endif  !}

return

end function otpm_set_diag_tracer  !}
! </FUNCTION> NAME="otpm_set_diag_tracer"


!#######################################################################
! <SUBROUTINE NAME="otpm_start_namelist">
!
! <DESCRIPTION>
! Set the values for a diag tracer and return its index (0 on error)
! </DESCRIPTION>
!
subroutine otpm_start_namelist(path, name, caller, no_overwrite, check)  !{

implicit none

!
!       arguments
!

character(len=*), intent(in)            :: path
character(len=*), intent(in)            :: name
character(len=*), intent(in), optional  :: caller
logical,          intent(in), optional  :: no_overwrite
logical,          intent(in), optional  :: check

!
!       Local parameters
!

character(len=48), parameter  :: sub_name = 'otpm_start_namelist'

!
!       Local variables
!

integer                         :: namelist_index
character(len=fm_path_name_len) :: path_name
character(len=256)              :: error_header
character(len=256)              :: warn_header
character(len=256)              :: note_header
character(len=128)              :: caller_str

!
!       set the caller string and headers
!

if (present(caller)) then  !{
  caller_str = '[' // trim(caller) // ']'
else  !}{
  caller_str = default_caller
endif  !}

error_header = '==>Error from ' // trim(mod_name) //   &
               '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a name is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

!
!       Concatenate the path and name
!

if (path .eq. ' ') then  !{
  path_name = name
else  !}{
  path_name = trim(path) // '/' // name
endif  !}
save_path = path
save_name = name

!
!       set the default caller string, if desired
!

if (present(caller)) then  !{
  call otpm_set_caller(caller)
else  !}{
  call otpm_reset_caller
endif  !}

!
!       set the default no_overwrite flag, if desired
!

if (present(no_overwrite)) then  !{
  call otpm_set_no_overwrite(no_overwrite)
else  !}{
  call otpm_reset_no_overwrite
endif  !}

!
!       set the default good_name_list string, if desired
!

if (present(check)) then  !{
  if (check) then  !{
    call otpm_set_good_name_list('/ocean_mod/GOOD/namelists/' // trim(path_name) // '/good_list')
  else  !}{
    call otpm_reset_good_name_list
  endif  !}
else  !}{
  call otpm_reset_good_name_list
endif  !}

!
!       Process the namelist
!

write (stdout(),*)
write (stdout(),*) trim(note_header), ' Processing namelist ', trim(path_name)

!
!       Check whether the namelist already exists. If so, then use that one
!

namelist_index = fm_get_index('/ocean_mod/namelists/' // trim(path_name))
if (namelist_index .gt. 0) then  !{

  !write (stdout(),*) trim(note_header), ' Namelist already set with index ', namelist_index

else  !}{

!
!       Set a new namelist and get its index
!

  namelist_index = fm_new_list('/ocean_mod/namelists/' // trim(path_name), create = .true.)
  if (namelist_index .le. 0) then  !{
    call mpp_error(FATAL, trim(error_header) // ' Could not set namelist ' // trim(path_name))
  endif  !}

endif  !}

!
!       Add the namelist name to the list of good namelists, to be used
!       later for a consistency check
!

if (fm_new_value('/ocean_mod/GOOD/namelists/' // trim(path) // '/good_values',    &
                 name, append = .true., create = .true.) .le. 0) then  !{
  call mpp_error(FATAL, trim(error_header) //                           &
       ' Could not add ' // trim(name) // ' to "' // trim(path) // '/good_values" list')
endif  !}

!
!       Change to the new namelist, first saving the current list
!

save_current_list = fm_get_current_list()
if (save_current_list .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not get the current list')
endif  !}

if (.not. fm_change_list('/ocean_mod/namelists/' // trim(path_name))) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not change to the namelist ' // trim(path_name))
endif  !}

return

end subroutine otpm_start_namelist  !}
! </SUBROUTINE> NAME="otpm_start_namelist"


!#######################################################################
! <SUBROUTINE NAME="otpm_end_namelist">
!
! <DESCRIPTION>
! Set the values for a diag tracer and return its index (0 on error)
! </DESCRIPTION>
!
subroutine otpm_end_namelist(path, name, caller, check)  !{

implicit none

!
!       arguments
!

character(len=*), intent(in)            :: path
character(len=*), intent(in)            :: name
character(len=*), intent(in), optional  :: caller
logical,          intent(in), optional  :: check

!
!       Local parameters
!

character(len=48), parameter  :: sub_name = 'otpm_end_namelist'

!
!       Local variables
!

character(len=fm_string_len), pointer, dimension(:)     :: good_list => NULL()
character(len=fm_path_name_len)                         :: path_name
character(len=256)                                      :: error_header
character(len=256)                                      :: warn_header
character(len=256)                                      :: note_header
character(len=128)                                      :: caller_str

!
!       set the caller string and headers
!

if (present(caller)) then  !{
  caller_str = '[' // trim(caller) // ']'
else  !}{
  caller_str = default_caller
endif  !}

error_header = '==>Error from ' // trim(mod_name) //   &
               '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
warn_header = '==>Warning from ' // trim(mod_name) //  &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
note_header = '==>Note from ' // trim(mod_name) //     &
              '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

!
!       check that a path is given (fatal if not)
!

if (name .eq. ' ') then  !{
  call mpp_error(FATAL, trim(error_header) // ' Empty name given')
endif  !}

!
!       Check that the path ane name match the preceding call to
!       otpm_start_namelist
!

if (path .ne. save_path) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Path "' // trim(path) // '" does not match saved path "' // trim(save_path) // '"')
elseif (name .ne. save_name) then  !}{
  call mpp_error(FATAL, trim(error_header) // ' Name "' // trim(name) // '" does not match saved name "' // trim(save_name) // '"')
endif  !}

!
!       Concatenate the path and name
!

if (path .eq. ' ') then  !{
  path_name = name
else  !}{
  path_name = trim(path) // '/' // name
endif  !}
save_path = ' '
save_name = ' '

!
!       Check for any errors in the number of fields in this list
!

if (present(check)) then  !{
  if (check) then  !{
    if (caller_str .eq. ' ') then  !{
      caller_str = trim(mod_name) // '(' // trim(sub_name) // ')'
    endif  !}
    good_list => otpm_get_string_array('/ocean_mod/GOOD/namelists/' // trim(path_name) // '/good_list',            &
         caller = trim(mod_name) // '(' // trim(sub_name) // ')')
    if (associated(good_list)) then  !{
      call otpm_check_for_bad_fields('/ocean_mod/namelists/' // trim(path_name), good_list, caller = caller_str)
      deallocate(good_list)
    else  !}{
      call mpp_error(FATAL, trim(error_header) // ' Empty "' // trim(path_name) // '" list')
    endif  !}
  endif  !}
endif  !}

!
!       Change back to the saved list
!

if (save_current_list .ne. ' ') then  !{
  if (.not. fm_change_list(save_current_list)) then  !{
    call mpp_error(FATAL, trim(error_header) // ' Could not change to the saved list: ' // trim(save_current_list))
  endif  !}
endif  !}
save_current_list = ' '

!
!       reset the default caller string
!

call otpm_reset_caller

!
!       reset the default no_overwrite string
!

call otpm_reset_no_overwrite

!
!       reset the default good_name_list string
!

call otpm_reset_good_name_list

return

end subroutine otpm_end_namelist  !}
! </SUBROUTINE> NAME="otpm_end_namelist"


end module ocean_tpm_util_mod  !}
