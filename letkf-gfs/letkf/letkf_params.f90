!--------------------------------------------------------------------------------
! (GFS-LETKF) letkf_params.f90
!
! Hold the user configurable parameters for the LETKF, reads in the settings from
! a namelist
!--------------------------------------------------------------------------------
MODULE letkf_params

USE common, ONLY: r_size

IMPLICIT NONE

PUBLIC

!-------------------------------------------------------------------------------
! Details for LETKF localization and quality control
!-------------------------------------------------------------------------------
! INTEGER  :: nslots=1                 ! number of time slots for 4D-LETKF
! INTEGER  :: nbslot=1                 !STEVE: nbslot=1 for testing for GMAO example case. Normal case is nbslot=5 ! basetime slot
! REAL(r_size) :: gross_error=3.0d0    ! number of standard deviations
!                                      ! used to filter out observations
! real(r_size) :: sigma_ocn_h(2)  
! real(r_size) :: sigma_atm_h(2)

real(r_size) :: sigma_atm_h(2)
real(r_size) :: sigma_atm_v
real(r_size) :: sigma_atm_t

real(r_size) :: sigma_ocn_h(2)
real(r_size) :: sigma_ocn_v
real(r_size) :: sigma_ocn_t

!-------------------------------------------------------------------------------
! Ensemble Size
!-------------------------------------------------------------------------------
INTEGER   :: nbv=90

contains

  
  subroutine params_init()
    logical :: ex

    namelist /letkf/ nbv

    namelist /letkf_atm/ &
         sigma_atm_h, sigma_atm_v, sigma_atm_t, &
         sigma_ocn_h, sigma_ocn_v, sigma_ocn_t         


    open(99, file="letkf.nml", status="old")
    read(99, nml=letkf)

    rewind(99)   
    read(99, nml=letkf_atm)

    close(99)

    write(6,letkf)
    write(6,letkf_atm)
    
  end subroutine params_init


  
END MODULE letkf_params
