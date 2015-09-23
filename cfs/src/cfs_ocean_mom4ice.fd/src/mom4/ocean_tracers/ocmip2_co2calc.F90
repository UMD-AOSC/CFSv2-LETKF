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
!
! 
!<CONTACT EMAIL="Richard.Slater@noaa.gov"> Richard D. Slater
!</CONTACT>
!
!<REVIEWER EMAIL="John.Dunne@noaa.gov"> John P. Dunne
!</REVIEWER>
!
!<OVERVIEW>
! Surface fCO2 calculation
!</OVERVIEW>
!
!<DESCRIPTION>
! Calculate the fugacity of CO2 at the surface in thermodynamic
! equilibrium with the current alkalinity (Alk) and total dissolved
! inorganic carbon (DIC) at a particular temperature and salinity
! using an initial guess for the total hydrogen
! ion concentration (htotal)
!</DESCRIPTION>
!

!
! $Id$
!

module ocmip2_co2calc_mod  !{

!
!------------------------------------------------------------------
!
!       modules
!
!------------------------------------------------------------------
!

use time_manager_mod, only : time_type, operator( .ne. ), set_time

!
!------------------------------------------------------------------
!
!       Global definitions
!
!------------------------------------------------------------------
!

implicit none

private

public  :: ocmip2_co2calc
public  :: init_ocmip2_co2calc

logical                                 :: initialized = .false.
integer                                 :: imt
integer                                 :: km
integer                                 :: indexi
integer                                 :: indexj
integer                                 :: indexk
real                                    :: log100
real, allocatable, dimension(:,:,:)     :: pt
real, allocatable, dimension(:,:,:)     :: sit
real, allocatable, dimension(:,:,:)     :: ta
real, allocatable, dimension(:,:,:)     :: dic
real, allocatable, dimension(:,:)       :: xco2
real, allocatable, dimension(:,:,:)     :: tk
real, allocatable, dimension(:,:,:)     :: invtk
real, allocatable, dimension(:,:,:)     :: tk100
real, allocatable, dimension(:,:,:)     :: tk1002
real, allocatable, dimension(:,:,:)     :: dlogtk
real, allocatable, dimension(:,:,:)     :: is
real, allocatable, dimension(:,:,:)     :: is2
real, allocatable, dimension(:,:,:)     :: sqrtis
real, allocatable, dimension(:,:,:)     :: s2
real, allocatable, dimension(:,:,:)     :: sqrts
real, allocatable, dimension(:,:,:)     :: s15
real, allocatable, dimension(:,:,:)     :: scl
real, allocatable, dimension(:,:,:)     :: logf_of_s
real, allocatable, dimension(:,:,:)     :: ff
real, allocatable, dimension(:,:,:)     :: k0
real, allocatable, dimension(:,:,:)     :: k1
real, allocatable, dimension(:,:,:)     :: k2
real, allocatable, dimension(:,:,:)     :: kb
real, allocatable, dimension(:,:,:)     :: k1p
real, allocatable, dimension(:,:,:)     :: k2p
real, allocatable, dimension(:,:,:)     :: k3p
real, allocatable, dimension(:,:,:)     :: ksi
real, allocatable, dimension(:,:,:)     :: kw
real, allocatable, dimension(:,:,:)     :: ks
real, allocatable, dimension(:,:,:)     :: kf
real, allocatable, dimension(:,:,:)     :: bt
real, allocatable, dimension(:,:,:)     :: st
real, allocatable, dimension(:,:,:)     :: ft
real, allocatable, dimension(:,:,:)     :: htotal2
type(time_type)                         :: time_calculated

!
!-----------------------------------------------------------------------
!
!       Subroutine and function definitions
!
!-----------------------------------------------------------------------
!

contains


!#######################################################################
! <SUBROUTINE NAME="init_ocmip2_co2calc">
!
! <DESCRIPTION>
!       Initialize parameters to be used for the calculation of
! delta co2*.
!
! INPUT
!       time       = the time when this routine is called (saved so
!                       the routine may be called from multiple
!                       places, but calculations only performed once)
!
!       imt_in     = i-dimension of the arrays
!
!       km_in      = k-dimension of the arrays
!
!       t          = temperature (degrees C)
!
!       s          = salinity (PSU)
!
! OUTPUT
!
!       k2         = activity factors for carbonate species
!
!       k2              (see below)
!
!       invtk      = 1/(t+273.15)
! </DESCRIPTION>

subroutine init_ocmip2_co2calc(time, isc, iec, jsc, jec, nk,    &
     t, s, k1_out, k2_out, invtk_out)  !{

implicit none

!
!       local parameters
!

!
!       arguments
!

integer, intent(in)                                     :: isc
integer, intent(in)                                     :: iec
integer, intent(in)                                     :: jsc
integer, intent(in)                                     :: jec
integer, intent(in)                                     :: nk
real, dimension(isc:iec, jsc:jec, nk), intent(out),             &
     optional                                           :: invtk_out
real, dimension(isc:iec, jsc:jec, nk), intent(out),             &
     optional                                           :: k1_out
real, dimension(isc:iec, jsc:jec, nk), intent(out),             &
     optional                                           :: k2_out
real, dimension(isc:iec, jsc:jec, nk), intent(in)       :: s
real, dimension(isc:iec, jsc:jec, nk), intent(in)       :: t
type(time_type), intent(in)                             :: time

!
!       local variables
!

integer :: i
integer :: j
integer :: k

!
!       Initialize the module
!

if (.not. initialized) then  !{
  call alloc_ocmip2_co2calc(isc, iec, jsc, jec, nk)
endif  !}

!
!---------------------------------------------------------------------
!
!***********************************************************************
! Calculate all constants needed to convert between various measured
! carbon species. References for each equation are noted in the code.
! Once calculated, the constants are stored and passed in the common
! block "const". The original version of this code was based on
! the code by Dickson in Version 2 of "Handbook of Methods for the
! Analysis of the Various Parameters of the Carbon Dioxide System
! in Seawater", DOE, 1994 (SOP No. 3, p25-26).
!
! Derive simple terms used more than once
!

if (time .ne. time_calculated) then  !{
  do k = 1, nk  !{
    do j = jsc, jec  !{
      do i = isc, iec  !{
        tk(i,j,k)        = 273.15 + t(i,j,k)
        tk100(i,j,k)     = tk(i,j,k) / 100.0
        tk1002(i,j,k)    = tk100(i,j,k) * tk100(i,j,k)
        invtk(i,j,k)     = 1.0 / tk(i,j,k)
        dlogtk(i,j,k)    = log(tk(i,j,k))
        is(i,j,k)        = 19.924 * s(i,j,k) /                  &
                           (1000.0 - 1.005 * s(i,j,k))
        is2(i,j,k)       = is(i,j,k) * is(i,j,k)
        sqrtis(i,j,k)    = sqrt(is(i,j,k))
        s2(i,j,k)        = s(i,j,k) * s(i,j,k)
        sqrts(i,j,k)     = sqrt(s(i,j,k))
        s15(i,j,k)       = sqrts(i,j,k) ** 3
        scl(i,j,k)       = s(i,j,k) / 1.80655
        logf_of_s(i,j,k) = log(1.0 - 0.001005 * s(i,j,k))

!
! f = k0(1-pH2O)*correction term for non-ideality
!
! Weiss & Price (1980, Mar. Chem., 8, 347-359; Eq 13 with table 6
!                values)
!

        ff(i,j,k) = exp(-162.8301 + 218.2968 / tk100(i,j,k)  +  &
                        90.9241 *                               &
                        (dlogtk(i,j,k) - log100) -              &
                        1.47696 * tk1002(i,j,k) +               &
                        s(i,j,k) *                              &
                        (0.025695 - 0.025225 * tk100(i,j,k) +   &
                         0.0049867 * tk1002(i,j,k)))

!
! k0 from Weiss 1974
!

        k0(i,j,k) = exp(93.4517/tk100(i,j,k) - 60.2409 +        &
                        23.3585 * log(tk100(i,j,k)) +           &
                        s(i,j,k) *                              &
                        (0.023517 - 0.023656 * tk100(i,j,k) +   &
                         0.0047036 * tk1002(i,j,k)))

!
! k1 = [H][HCO3]/[H2CO3]
! k2 = [H][CO3]/[HCO3]
!
! Millero p.664 (1995) using Mehrbach et al. data on seawater scale 
!

        k1(i,j,k) = 10.0 ** (-(3670.7 * invtk(i,j,k) -          &
                               62.008 +                         &
                               9.7944 * dlogtk(i,j,k) -         &
                               0.0118 * s(i,j,k) +              &
                               0.000116 * s2(i,j,k)))

        k2(i,j,k) = 10.0 ** (-(1394.7 * invtk(i,j,k) + 4.777 -  &
                               0.0184 * s(i,j,k) +              &
                               0.000118 * s2(i,j,k)))

!
! kb = [H][BO2]/[HBO2]
!
! Millero p.669 (1995) using data from Dickson (1990)
!

        kb(i,j,k) = exp((-8966.90 - 2890.53 * sqrts(i,j,k) -    &
                         77.942 * s(i,j,k) +                    &
                         1.728 * s15(i,j,k) -                   &
                         0.0996 * s2(i,j,k)) * invtk(i,j,k) +   &
                        (148.0248 + 137.1942 * sqrts(i,j,k) +   &
                         1.62142 * s(i,j,k)) +                  &
                        (-24.4344 - 25.085 * sqrts(i,j,k) -     &
                         0.2474 * s(i,j,k)) * dlogtk(i,j,k) +   &
                        0.053105 * sqrts(i,j,k) * tk(i,j,k))

!
! k1p = [H][H2PO4]/[H3PO4]
!
! DOE(1994) eq 7.2.20 with footnote using data from Millero (1974)
!

        k1p(i,j,k) = exp(-4576.752 * invtk(i,j,k) + 115.525 -   &
                         18.453 * dlogtk(i,j,k) +               &
                         (-106.736 * invtk(i,j,k) + 0.69171) *  &
                         sqrts(i,j,k) +                         &
                         (-0.65643 * invtk(i,j,k) - 0.01844) * s(i,j,k))

!
! k2p = [H][HPO4]/[H2PO4]
!
! DOE(1994) eq 7.2.23 with footnote using data from Millero (1974))
!

        k2p(i,j,k) = exp(-8814.715 * invtk(i,j,k) + 172.0883 -  &
                         27.927 * dlogtk(i,j,k) +               &
                         (-160.340 * invtk(i,j,k) + 1.3566) *   &
                         sqrts(i,j,k) +                         &
                         (0.37335 * invtk(i,j,k) - 0.05778) * s(i,j,k))

!
!-----------------------------------------------------------------------
! k3p = [H][PO4]/[HPO4]
!
! DOE(1994) eq 7.2.26 with footnote using data from Millero (1974)
!

        k3p(i,j,k) = exp(-3070.75 * invtk(i,j,k) - 18.141 +     &
                         (17.27039 * invtk(i,j,k) + 2.81197) *  &
                         sqrts(i,j,k) +                         &
                         (-44.99486 * invtk(i,j,k) - 0.09984) * &
                         s(i,j,k))

!
!-----------------------------------------------------------------------
! ksi = [H][SiO(OH)3]/[Si(OH)4]
!
! Millero p.671 (1995) using data from Yao and Millero (1995)
!

        ksi(i,j,k) = exp(-8904.2 * invtk(i,j,k) + 117.385 -     &
                         19.334 * dlogtk(i,j,k) +               &
                         (-458.79 * invtk(i,j,k) + 3.5913) *    &
                         sqrtis(i,j,k) +                        &
                         (188.74 * invtk(i,j,k) - 1.5998) *     &
                         is(i,j,k) +                            &
                         (-12.1652 * invtk(i,j,k) + 0.07871) *  &
                         is2(i,j,k) + logf_of_s(i,j,k))

!
!-----------------------------------------------------------------------
! kw = [H][OH]
!
! Millero p.670 (1995) using composite data
!

        kw(i,j,k) = exp(-13847.26 * invtk(i,j,k) + 148.9652 -   &
                        23.6521 * dlogtk(i,j,k) +               &
                        (118.67 * invtk(i,j,k) - 5.977 +        &
                         1.0495 * dlogtk(i,j,k)) *              &
                        sqrts(i,j,k) - 0.01615 * s(i,j,k))

!
!-----------------------------------------------------------------------
! ks = [H][SO4]/[HSO4]
!
! Dickson (1990, J. chem. Thermodynamics 22, 113)
!

        ks(i,j,k) = exp(-4276.1 * invtk(i,j,k) + 141.328 -      &
                        23.093 * dlogtk(i,j,k) +                &
                        (-13856.0 * invtk(i,j,k) + 324.57 -     &
                         47.986 * dlogtk(i,j,k)) *              &
                        sqrtis(i,j,k) +                         &
                        (35474.0 * invtk(i,j,k) - 771.54 +      &
                         114.723 * dlogtk(i,j,k)) * is(i,j,k) - &
                        2698.0 * invtk(i,j,k) *                 &
                        sqrtis(i,j,k) ** 3 +                    &
                        1776.0 * invtk(i,j,k) * is2(i,j,k) +    &
                        logf_of_s(i,j,k))

!
!-----------------------------------------------------------------------
! kf = [H][F]/[HF]
!
! Dickson and Riley (1979) -- change pH scale to total
!

        kf(i,j,k) = exp(1590.2 * invtk(i,j,k) - 12.641 +        &
                        1.525 * sqrtis(i,j,k) +                 &
                        logf_of_s(i,j,k) +                      &
                        log(1.0 + (0.1400 / 96.062) *           &
                            scl(i,j,k) / ks(i,j,k)))

!
!-----------------------------------------------------------------------
! Calculate concentrations for borate, sulfate, and fluoride
!
! Uppstrom (1974)
!

        bt(i,j,k) = 0.000232 / 10.811 * scl(i,j,k)

!
! Morris & Riley (1966)
!

        st(i,j,k) = 0.14 / 96.062 * scl(i,j,k)

!
! Riley (1965)
!

        ft(i,j,k) = 0.000067 / 18.9984 * scl(i,j,k)

      enddo  !} i
    enddo  !} j
  enddo  !} k

!
!       Save the time
!

  time_calculated = time

endif  !}

!
!       set some stuff to pass back, if requested
!

if (present(k1_out)) then  !{
  do k = 1, nk  !{
    do j = jsc, jec  !{
      do i = isc, iec  !{
        k1_out(i,j,k) = k1(i,j,k)
      enddo  !} i
    enddo  !} j
  enddo  !} k
endif  !}

if (present(k2_out)) then  !{
  do k = 1, nk  !{
    do j = jsc, jec  !{
      do i = isc, iec  !{
        k2_out(i,j,k) = k2(i,j,k)
      enddo  !} i
    enddo  !} j
  enddo  !} k
endif  !}

if (present(invtk_out)) then  !{
  do k = 1, nk  !{
    do j = jsc, jec  !{
      do i = isc, iec  !{
        invtk_out(i,j,k)  = invtk(i,j,k)
      enddo  !} i
    enddo  !} j
  enddo  !} k
endif  !}

return

end subroutine  init_ocmip2_co2calc  !}
! </SUBROUTINE> NAME="init_ocmip2_co2calc"


!#######################################################################
! <SUBROUTINE NAME="ocmip2_co2calc">
!
! <DESCRIPTION>
!       Calculate delta co2* from total alkalinity and total CO2 at
! temperature (t), salinity (s) and "atmpres" atmosphere total pressure.
! It is assumed that init_ocmip2_co2calc has already been called with
! the T and S to calculate the various coefficients.
!
! INPUT
!
!       imt_in     = i-dimension of the arrays
!
!       km_in      = k-dimension of the arrays
!
!       mask       = land mask array (0.0 = land)
!
!       dic_in     = total inorganic carbon (mol/m^3) 
!                    where 1 T = 1 metric ton = 1000 kg
!
!       ta_in      = total alkalinity (eq/m^3) 
!
!       pt_in      = inorganic phosphate (mol/m^3) 
!
!       sit_in     = inorganic silicate (mol/m^3) 
!
!       htotallo   = lower limit of htotal range
!
!       htotalhi   = upper limit of htotal range
!
!       htotal     = H+ concentraion
!
!       xco2_in    = atmospheric mole fraction CO2 in dry air (ppmv) 
!
!       atmpres    = atmospheric pressure in atmospheres
!                    (1 atm = 1013.25 mbar)
!
! OUTPUT
!       co2star    = CO2*water (mol/m^3)
!       co2starair = CO2*air (mol/m^3)
!       dco2star   = Delta CO2* (mol/m^3)
!       pco2surf   = oceanic pCO2 (ppmv)
!
!       dpco2      = Delta pCO2, i.e, pCO2ocn - pCO2atm (ppmv)
!
! IMPORTANT: Some words about units - (JCO, 4/4/1999)
!
!     - Models carry tracers in mol/m^3 (on a per volume basis)
!
!     - Conversely, this routine, which was written by observationalists
!       (C. Sabine and R. Key), passes input arguments in umol/kg  
!       (i.e., on a per mass basis)
!
!     - I have changed things slightly so that input arguments are
!       in mol/m^3,
!
!     - Thus, all input concentrations (dic, ta, pt, and st) should be 
!       given in mol/m^3; output arguments "co2star" and "co2starair"  
!       and "dco2star" are likewise be in mol/m^3.
!
! FILES and PROGRAMS NEEDED: drtsafe, ta_iter_1
! </DESCRIPTION>

subroutine ocmip2_co2calc(isc, iec, jsc, jec, nk, mask,         &
     dic_in, ta_in, pt_in, sit_in, htotallo, htotalhi, htotal,  &
     xco2_in, atmpres, co2star, co2starair, dco2star,           &
     pCO2surf, dpco2)  !{

implicit none

!
!       local parameters
!

real, parameter :: permil = 1.0 / 1024.5
real, parameter :: permeg = 1.e-6
real, parameter :: xacc = 1.0e-10

!
!       arguments
!

integer, intent(in)                                     :: isc
integer, intent(in)                                     :: iec
integer, intent(in)                                     :: jsc
integer, intent(in)                                     :: jec
integer, intent(in)                                     :: nk
real, intent(in), dimension(isc:iec,jsc:jec,nk)         :: mask
real, intent(in), dimension(isc:iec,jsc:jec,nk)         :: dic_in
real, intent(in), dimension(isc:iec,jsc:jec,nk)         :: ta_in
real, intent(in), dimension(isc:iec,jsc:jec,nk)         :: pt_in
real, intent(in), dimension(isc:iec,jsc:jec,nk)         :: sit_in
real, intent(in), dimension(isc:iec,jsc:jec,nk)         :: htotallo
real, intent(in), dimension(isc:iec,jsc:jec,nk)         :: htotalhi
real, intent(inout), dimension(isc:iec,jsc:jec,nk)      :: htotal
real, intent(in), dimension(isc:iec,jsc:jec)            :: xco2_in
real, intent(in), dimension(isc:iec,jsc:jec)            :: atmpres
real, intent(out), dimension(isc:iec,jsc:jec,nk)        :: co2star
real, intent(out), dimension(isc:iec,jsc:jec)           :: co2starair
real, intent(out), dimension(isc:iec,jsc:jec)           :: dco2star
real, intent(out), dimension(isc:iec,jsc:jec)           :: pCO2surf
real, intent(out), dimension(isc:iec,jsc:jec)           :: dpco2

!
!       local variables
!

integer :: i, j, k
!real   :: x1, x2

!
!-----------------------------------------------------------------------
!       Change units from the input of mol/m^3 -> mol/kg:
!       (1 mol/m^3)  x (1 m^3/1024.5 kg)
!       where the ocean's mean surface density is 1024.5 kg/m^3
!       Note: mol/kg are actually what the body of this routine uses 
!       for calculations.  Units are reconverted back to mol/m^3 at the 
!       end of this routine.
!---------------------------------------------------------------------
!
!       To convert input in mol/m^3 -> mol/kg 
!

do k = 1, nk  !{
  do j = jsc, jec  !{
    do i = isc, iec  !{
      sit(i,j,k) = sit_in(i,j,k) * permil
      ta(i,j,k)  = ta_in(i,j,k)  * permil
      dic(i,j,k) = dic_in(i,j,k) * permil
      pt(i,j,k)  = pt_in(i,j,k)  * permil
    enddo  !} i
  enddo  !} j
enddo  !} k

!
!---------------------------------------------------------------------
!       Change units from uatm to atm. That is, atm is what the body of 
!       this routine uses for calculations.
!       Note: units are reconverted bac to uatm at END of this routine.
!---------------------------------------------------------------------
!
!       To convert input in uatm -> atm
!

do j = jsc, jec  !{
  do i = isc, iec  !{
    xco2(i,j) = xco2_in(i,j) * permeg
  enddo  !} i
enddo  !} j

!
!***********************************************************************
!
! Calculate [H+] total when DIC and TA are known at T, S and 1 atm.
! The solution converges to err of xacc. The solution must be within
! the range x1 to x2.
!
! If DIC and TA are known then either a root finding or iterative method
! must be used to calculate htotal. In this case we use the
! Newton-Raphson "safe" method taken from "Numerical Recipes"
! (function "rtsafe.f" with error trapping removed).
!
! As currently set, this procedure iterates about 12 times. The x1
! and x2 values set below will accomodate ANY oceanographic values.
! If an initial guess of the pH is known, then the number of
! iterations can be reduced to about 5 by narrowing the gap between
! x1 and x2. It is recommended that the first few time steps be run
! with x1 and x2 set as below. After that, set x1 and x2 to the
! previous value of the pH +/- ~0.5. The current setting of xacc will
! result in co2star accurate to 3 significant figures (xx.y). Making
! xacc bigger will result in faster convergence also, but this is not
! recommended (xacc of 10**-9 drops precision to 2 significant
! figures).
!

do k = 1, nk  !{
  do j = jsc, jec  !{
    do i = isc, iec  !{
      if (mask(i,j,k) .ne. 0.0) then  !{
        indexi = i              ! this is used in ta_iter_1 below
        indexj = j              ! this is used in ta_iter_1 below
        indexk = k              ! this is used in ta_iter_1 below
        htotal(i,j,k) = drtsafe(htotalhi(i,j,k), htotallo(i,j,k), xacc)
      endif  !}
    enddo  !} i
  enddo  !} j
enddo  !} k

!
! Calculate [CO2*] as defined in DOE Methods Handbook 1994 Ver.2, 
! ORNL/CDIAC-74, Dickson and Goyet, eds. (Ch 2 p 10, Eq A.49)
!

do k = 1, nk  !{
  do j = jsc, jec  !{
    do i = isc, iec  !{
      htotal2(i,j,k)    = htotal(i,j,k) * htotal(i,j,k)
      co2star(i,j,k)    = dic(i,j,k) * htotal2(i,j,k) /         &
                          (htotal2(i,j,k) +                     &
                           k1(i,j,k) * htotal(i,j,k) +          &
                           k1(i,j,k) * k2(i,j,k))
    enddo  !} i
  enddo  !} j
enddo  !} k
!ph         = -log10(htotal)

do j = jsc, jec  !{
  do i = isc, iec  !{
    co2starair(i,j) = xco2(i,j) * ff(i,j,1) * atmpres(i,j)
    dco2star(i,j)   = co2starair(i,j) - co2star(i,j,1)
  enddo  !} i
enddo  !} j

!
!---------------------------------------------------------------
!c      Add two output arguments for storing pCO2surf
!c      Should we be using K0 or ff for the solubility here?
!---------------------------------------------------------------
!

do j = jsc, jec  !{
  do i = isc, iec  !{
    pCO2surf(i,j) = co2star(i,j,1) / ff(i,j,1)
    dpCO2(i,j)    = pCO2surf(i,j) - xco2(i,j) * atmpres(i,j)
  enddo  !} i
enddo  !} j

!
!----------------------------------------------------------------
!
! Convert units of output arguments
!      Note: dco2star, co2star and co2starair are calculated in
!            mol/kg within this routine 
!      Thus Convert now from mol/kg -> mol/m^3
!

do k = 1, nk  !{
  do j = jsc, jec  !{
    do i = isc, iec  !{
      co2star(i,j,k) = co2star(i,j,k) / permil
    enddo  !} i
  enddo  !} j
enddo  !} k
do j = jsc, jec  !{
  do i = isc, iec  !{
    dco2star(i,j)   = dco2star(i,j) / permil
    co2starair(i,j) = co2starair(i,j) / permil
  enddo  !} i
enddo  !} j

!
!      Note: pCO2surf and dpCO2 are calculated in atm above. 
!      Thus convert now to uatm
!

do j = jsc, jec  !{
  do i = isc, iec  !{
    pCO2surf(i,j) = pCO2surf(i,j) / permeg
    dpCO2(i,j)    = dpCO2(i,j) / permeg
  enddo  !} i
enddo  !} j

return

end subroutine  ocmip2_co2calc  !}
! </SUBROUTINE> NAME="ocmip2_co2calc"


!#######################################################################
! <FUNCTION NAME="drtsafe">
!
! <DESCRIPTION>
!       File taken from Numerical Recipes. Modified  R. M. Key 4/94
! </DESCRIPTION>

function drtsafe(x1,x2,xacc)  !{

implicit none

!
!       arguments
!

real    :: drtsafe
real    :: x1, x2, xacc

!
!       local parameters
!

integer, parameter      :: maxit = 100

!
!       local variables
!

integer :: j
integer :: num
real    :: fl, df, fh, swap, xl, xh, dxold, dx, f, temp

call ta_iter_1(x1,fl,df)
call ta_iter_1(x2,fh,df)
if(fl .lt. 0.0) then
  xl=x1
  xh=x2
else
  xh=x1
  xl=x2
  swap=fl
  fl=fh
  fh=swap
end if
drtsafe=0.5*(x1+x2)
dxold=abs(x2-x1)
dx=dxold
call ta_iter_1(drtsafe,f,df)
do j=1,maxit  !{
  if (((drtsafe-xh)*df-f)*((drtsafe-xl)*df-f) .ge. 0.0 .or.     &
      abs(2.0*f) .gt. abs(dxold*df)) then
    dxold=dx
    dx=0.5*(xh-xl)
    drtsafe=xl+dx
    if (xl .eq. drtsafe) then
!     write (6,*) 'Exiting drtsafe at A on iteration  ', j, ', ph = ', -log10(drtsafe)
      return
    endif
  else
    dxold=dx
    dx=f/df
    temp=drtsafe
    drtsafe=drtsafe-dx
    if (temp .eq. drtsafe) then
!     write (6,*) 'Exiting drtsafe at B on iteration  ', j, ', ph = ', -log10(drtsafe)
      return
    endif
  end if
  if (abs(dx) .lt. xacc) then
!     write (6,*) 'Exiting drtsafe at C on iteration  ', j, ', ph = ', -log10(drtsafe)
    return
  endif
  call ta_iter_1(drtsafe,f,df)
  if(f .lt. 0.0) then
    xl=drtsafe
    fl=f
  else
    xh=drtsafe
    fh=f
  end if
enddo  !} j

return

end  function  drtsafe  !}
! </FUNCTION> NAME="drtsafe"


!#######################################################################
! <SUBROUTINE NAME="ta_iter_1">
!
! <DESCRIPTION>
! This routine expresses TA as a function of DIC, htotal and constants.
! It also calculates the derivative of this function with respect to 
! htotal. It is used in the iterative solution for htotal. In the call
! "x" is the input value for htotal, "fn" is the calculated value for TA
! and "df" is the value for dTA/dhtotal
! </DESCRIPTION>

subroutine ta_iter_1(x, fn, df)  !{

implicit none

!
!       arguments
!

real    :: x, fn, df

!
!       local variables
!

real    :: x2, x3, k12, k12p, k123p, c, a, a2, da, b, b2, db

x2 = x*x
x3 = x2*x
k12 = k1(indexi,indexj,indexk)*k2(indexi,indexj,indexk)
k12p = k1p(indexi,indexj,indexk)*k2p(indexi,indexj,indexk)
k123p = k12p*k3p(indexi,indexj,indexk)
c = 1.0 + st(indexi,indexj,indexk)/ks(indexi,indexj,indexk)
a = x3 + k1p(indexi,indexj,indexk)*x2 + k12p*x + k123p
a2 = a*a
da = 3.0*x2 + 2.0*k1p(indexi,indexj,indexk)*x + k12p
b = x2 + k1(indexi,indexj,indexk)*x + k12
b2 = b*b
db = 2.0*x + k1(indexi,indexj,indexk)

!
!     fn = hco3+co3+borate+oh+hpo4+2*po4+silicate+hfree+hso4+hf+h3po4-ta
!

fn = k1(indexi,indexj,indexk)*x*dic(indexi,indexj,indexk)/b +   &
     2.0*dic(indexi,indexj,indexk)*k12/b +                      &
     bt(indexi,indexj,indexk)/                                  &
       (1.0 + x/kb(indexi,indexj,indexk)) +                     &
     kw(indexi,indexj,indexk)/x +                               &
     pt(indexi,indexj,indexk)*k12p*x/a +                        &
     2.0*pt(indexi,indexj,indexk)*k123p/a +                     &
     sit(indexi,indexj,indexk)/                                 &
       (1.0 + x/ksi(indexi,indexj,indexk)) -                    &
     x/c -                                                      &
     st(indexi,indexj,indexk)/                                  &
       (1.0 + ks(indexi,indexj,indexk)/x/c) -                   &
     ft(indexi,indexj,indexk)/                                  &
       (1.0 + kf(indexi,indexj,indexk)/x) -                     &
     pt(indexi,indexj,indexk)*x3/a -                            &
     ta(indexi,indexj,indexk)

!
!     df = dfn/dx
!

df = ((k1(indexi,indexj,indexk)*dic(indexi,indexj,indexk)*b) -  &
      k1(indexi,indexj,indexk)*x*                               &
        dic(indexi,indexj,indexk)*db)/b2 -                      &
     2.0*dic(indexi,indexj,indexk)*k12*db/b2 -                  &
     bt(indexi,indexj,indexk)/kb(indexi,indexj,indexk)/         &
        (1.0+x/kb(indexi,indexj,indexk))**2 -                   &
     kw(indexi,indexj,indexk)/x2 +                              &
     (pt(indexi,indexj,indexk)*k12p*(a - x*da))/a2 -            &
     2.0*pt(indexi,indexj,indexk)*k123p*da/a2 -                 &
     sit(indexi,indexj,indexk)/ksi(indexi,indexj,indexk)/       &
        (1.0+x/ksi(indexi,indexj,indexk))**2 -                  &
     1.0/c +                                                    &
     st(indexi,indexj,indexk)*                                  &
       (1.0 + ks(indexi,indexj,indexk)/x/c)**(-2)*              &
       (ks(indexi,indexj,indexk)/c/x2) +                        &
     ft(indexi,indexj,indexk)*                                  &
       (1.0 + kf(indexi,indexj,indexk)/x)**(-2)*                &
       kf(indexi,indexj,indexk)/x2 -                            &
     pt(indexi,indexj,indexk)*x2*(3.0*a-x*da)/a2

return

end subroutine  ta_iter_1  !}
! </SUBROUTINE> NAME="ta_iter_1"


!#######################################################################
! <SUBROUTINE NAME="alloc_ocmip2_co2calc">
!
! <DESCRIPTION>
!       initialize arrays for use in ocmip2_co2calc subroutine
! </DESCRIPTION>

subroutine alloc_ocmip2_co2calc(isc, iec, jsc, jec, nk)  !{

implicit none

!
!       arguments
!

integer, intent(in)     :: isc
integer, intent(in)     :: iec
integer, intent(in)     :: jsc
integer, intent(in)     :: jec
integer, intent(in)     :: nk

!
!       local variables
!

if (.not. initialized) then  !{
  time_calculated = set_time(0)
  log100 = log(100.0)
  allocate( pt(isc:iec,jsc:jec,nk) )
  allocate( sit(isc:iec,jsc:jec,nk) )
  allocate( ta(isc:iec,jsc:jec,nk) )
  allocate( dic(isc:iec,jsc:jec,nk) )
  allocate( xco2(isc:iec,jsc:jec) )
  allocate( tk(isc:iec,jsc:jec,nk) )
  allocate( invtk(isc:iec,jsc:jec,nk) )
  allocate( tk100(isc:iec,jsc:jec,nk) )
  allocate( tk1002(isc:iec,jsc:jec,nk) )
  allocate( dlogtk(isc:iec,jsc:jec,nk) )
  allocate( is(isc:iec,jsc:jec,nk) )
  allocate( is2(isc:iec,jsc:jec,nk) )
  allocate( sqrtis(isc:iec,jsc:jec,nk) )
  allocate( s2(isc:iec,jsc:jec,nk) )
  allocate( sqrts(isc:iec,jsc:jec,nk) )
  allocate( s15(isc:iec,jsc:jec,nk) )
  allocate( scl(isc:iec,jsc:jec,nk) )
  allocate( logf_of_s(isc:iec,jsc:jec,nk) )
  allocate( ff(isc:iec,jsc:jec,nk) )
  allocate( k0(isc:iec,jsc:jec,nk) )
  allocate( k1(isc:iec,jsc:jec,nk) )
  allocate( k2(isc:iec,jsc:jec,nk) )
  allocate( kb(isc:iec,jsc:jec,nk) )
  allocate( k1p(isc:iec,jsc:jec,nk) )
  allocate( k2p(isc:iec,jsc:jec,nk) )
  allocate( k3p(isc:iec,jsc:jec,nk) )
  allocate( ksi(isc:iec,jsc:jec,nk) )
  allocate( kw(isc:iec,jsc:jec,nk) )
  allocate( ks(isc:iec,jsc:jec,nk) )
  allocate( kf(isc:iec,jsc:jec,nk) )
  allocate( bt(isc:iec,jsc:jec,nk) )
  allocate( st(isc:iec,jsc:jec,nk) )
  allocate( ft(isc:iec,jsc:jec,nk) )
  allocate( htotal2(isc:iec,jsc:jec,nk) )
  pt(:,:,:) = 0.0
  sit(:,:,:) = 0.0
  ta(:,:,:) = 0.0
  dic(:,:,:) = 0.0
  xco2(:,:) = 0.0
  tk(:,:,:) = 0.0
  invtk(:,:,:) = 0.0
  tk100(:,:,:) = 0.0
  tk1002(:,:,:) = 0.0
  dlogtk(:,:,:) = 0.0
  is(:,:,:) = 0.0
  is2(:,:,:) = 0.0
  sqrtis(:,:,:) = 0.0
  s2(:,:,:) = 0.0
  sqrts(:,:,:) = 0.0
  s15(:,:,:) = 0.0
  scl(:,:,:) = 0.0
  logf_of_s(:,:,:) = 0.0
  ff(:,:,:) = 0.0
  k0(:,:,:) = 0.0
  k1(:,:,:) = 0.0
  k2(:,:,:) = 0.0
  kb(:,:,:) = 0.0
  k1p(:,:,:) = 0.0
  k2p(:,:,:) = 0.0
  k3p(:,:,:) = 0.0
  ksi(:,:,:) = 0.0
  kw(:,:,:) = 0.0
  ks(:,:,:) = 0.0
  kf(:,:,:) = 0.0
  bt(:,:,:) = 0.0
  st(:,:,:) = 0.0
  ft(:,:,:) = 0.0
  htotal2(:,:,:) = 0.0
  initialized = .true.
endif  !}

return

end subroutine  alloc_ocmip2_co2calc  !}
! </SUBROUTINE> NAME="alloc_ocmip2_co2calc"

end module  ocmip2_co2calc_mod  !}
