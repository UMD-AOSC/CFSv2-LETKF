!===============================================================================
module obs_tools
!-------------------------------------------------------------------------------

  use common
  use common_obs
  use common_gfs
  use common_obs_gfs
  implicit none

contains

!===============================================================================
subroutine print_obsnum (nobs, nobs_cat)
!-------------------------------------------------------------------------------

  implicit none
  integer, intent(in) :: nobs
  integer, intent(in) :: nobs_cat(obsid_num,platform_num+1)
  integer :: itype

  write (*, '(A)') '================================================================================'
  write (*, '(A,I10)') ' TOTAL NUMBER OF OBSERVATIONS:', nobs
  write (*, '(A)') '--------------------------------------------------------------------------------'
  
  write (*, '(A,8(4x,A3))') '        ', obsid_names(obsid_id2idx(obsid_atm_min): obsid_id2idx(obsid_atm_min))
  do itype = 1, platform_num
    write (*, '(A6,A2,8(I7))') platform_name(itype), ': ', nobs_cat(:,itype)
  end do
  write (*, '(A6,A2,8(I7))') 'OTHERS', ': ', nobs_cat(:,platform_num+1)
  
  write (*, '(A)') '================================================================================'

!-------------------------------------------------------------------------------
end subroutine print_obsnum
!===============================================================================

end module obs_tools
!===============================================================================
