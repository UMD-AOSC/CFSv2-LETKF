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

!! Localization
!! ------------------------------------------------------------
real(r_size) :: sigma_atm_h(2)
real(r_size) :: sigma_atm_v
real(r_size) :: sigma_atm_t

real(r_size) :: sigma_ocn_h(2)
real(r_size) :: sigma_ocn_v
real(r_size) :: sigma_ocn_t

!! Inflation
!! ------------------------------------------------------------
real(r_size) :: infl_rtps
real(r_size) :: infl_mult

!! Ensemble Size
!! ------------------------------------------------------------
INTEGER   :: nbv=90


logical :: ocn_obs = .false.
logical :: atm_obs = .true.
integer :: atm_obs_plat(10)


contains

  
  subroutine params_init()
    logical :: ex

    !! parameters used by all LETKFs in the coupled system
    namelist /letkf/ nbv

    !! parameters used by just the GFS-LETKF
    namelist /letkf_atm/ &
         sigma_atm_h, sigma_atm_v, sigma_atm_t, &
         sigma_ocn_h, sigma_ocn_v, sigma_ocn_t, &
         infl_rtps, infl_mult, &
         ocn_obs, atm_obs, atm_obs_plat

    open(99, file="letkf.nml", status="old")
    read(99, nml=letkf)

    rewind(99)   
    read(99, nml=letkf_atm)

    close(99)

    write(6,letkf)
    write(6,letkf_atm)
    
  end subroutine params_init


  
END MODULE letkf_params
