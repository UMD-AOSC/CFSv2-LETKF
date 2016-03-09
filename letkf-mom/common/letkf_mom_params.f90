!--------------------------------------------------------------------------------
! (MOM-LETKF) letkf_mom_params.f90
!
! Hold the user configurable parameters for the LETKF, reads in the settings from
! a namelist
!--------------------------------------------------------------------------------
MODULE letkf_mom_params

USE common, ONLY: r_size

IMPLICIT NONE

PUBLIC

INTEGER  :: nslots=1                 ! number of time slots for 4D-LETKF
INTEGER  :: nbslot=1                 !STEVE: nbslot=1 for testing for GMAO example case. Normal case is nbslot=5 ! basetime slot
REAL(r_size) :: gross_error=3.0d0    ! number of standard deviations
                                     ! used to filter out observations

!! Localization
!! ------------------------------------------------------------
real(r_size) :: sigma_ocn_h(2)
real(r_size) :: sigma_ocn_v

real(r_size) :: sigma_atm_h(2)
real(r_size) :: sigma_atm_v

!! inflation
!! ------------------------------------------------------------
real(r_size) :: infl_rtps
real(r_size) :: infl_mult

!! Ensemble Size
!! ------------------------------------------------------------
INTEGER  :: nbv

integer  :: use_obs_plat(10) !not used right now

!-------------------------------------------------------------------------------
! From letkf_local.f90
!-------------------------------------------------------------------------------

!STEVE: Testing "Vertical Tube" localization:
!       i.e. the localization is not applied vertically
! This provides the benefit that 
! (1) the analysis only has to be computed once
! per horizontal gridpoint, thus providing a nlevX (40X) speedup
! (2) the altimetry, SST, SSH, and bottom pressure (GRACE) can now be applied
! as direct constraints on the water column.
!
! There is precedence for this as in the paper "Reconstructing the Ocean's
! Interior from Surface Data" Wang et al. (2013)
!
!STEVE: making namelist accessible:
!LOGICAL :: DO_NO_VERT_LOC=.true. !STEVE: moved to letkf_local.f90 from letkf_tools.f90
!INTEGER :: localization_method=1 !1 !(OCEAN) =0 for uniform radius (default), =1 for latitude-dependent

!-------------------------------------------------------------------------------
! From letkf_tools.f90
!-------------------------------------------------------------------------------

! > 0: globally constant covariance inflation
! < 0: 3D inflation values input from a GPV file "infl_mul.grd"

REAL(r_size) :: sp_infl_add = 0.d0 !additive inflation

!-------------------------------------------------------------------------------
! From common_mom4.f90
!-------------------------------------------------------------------------------
! For (DRIFTERS)
LOGICAL :: DO_DRIFTERS=.false.
! For (ALTIMETRY)
LOGICAL :: DO_ALTIMETRY=.false.
! For (TRIPLOAR)
LOGICAL :: DO_TRIPOLAR=.true.
! For (IRREG_GRID) for non-orthogonal/non-rectilinear grids (not fully supported)
LOGICAL :: DO_IRREG_GRID=.false.

contains

  
  subroutine params_init()
    logical :: ex

    !! parameters used by all LETKFs in the coupled system
    namelist /letkf/ nbv

    !! parameters used by just the MOM-LETKF
    namelist /letkf_ocn/ nslots, nbslot, &
         sigma_ocn_h, sigma_ocn_v, &
         sigma_atm_h, sigma_atm_v, &        
         gross_error, &
         do_drifters, do_altimetry, &
         infl_rtps, infl_mult, &
         use_obs_plat
    
    open(99, file="input.nml", status="old")

    read(99, nml=letkf)
    rewind(99)   
    read(99, nml=letkf_ocn)

    close(99)
    
    write(6,letkf)
    write(6,letkf_ocn)
    
  end subroutine params_init


  
END MODULE letkf_mom_params
