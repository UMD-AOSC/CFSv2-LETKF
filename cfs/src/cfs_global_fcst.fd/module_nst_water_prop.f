MODULE module_nst_water_prop
  USE machine, ONLY : kind_phys
  USE module_nst_parameters, ONLY : t0K
  !
  PRIVATE
  PUBLIC :: rhocoef,density,sw_rad,sw_rad_Aw,sw_rad_sum,sw_rad_upper,sw_rad_upper_Aw,sw_rad_skin,grv,solar_time_from_Julian,compjd, &
            sw_ps_9b,sw_ps_9b_Aw
      
  !
  INTERFACE sw_ps_9b
     MODULE PROCEDURE sw_ps_9b
  END INTERFACE
  INTERFACE sw_ps_9b_Aw
     MODULE PROCEDURE sw_ps_9b_Aw
  END INTERFACE
  !
  INTERFACE sw_rad
     MODULE PROCEDURE sw_fairall_6exp_v1  ! sw_wick_v1
  END INTERFACE
  INTERFACE sw_rad_Aw
     MODULE PROCEDURE sw_fairall_6exp_v1_Aw
  END INTERFACE
  INTERFACE sw_rad_sum
     MODULE PROCEDURE sw_fairall_6exp_v1_sum
  END INTERFACE
  INTERFACE sw_rad_upper
     MODULE PROCEDURE sw_soloviev_3exp_v2
  END INTERFACE
  INTERFACE sw_rad_upper_Aw
     MODULE PROCEDURE sw_soloviev_3exp_v2_Aw
  END INTERFACE
  INTERFACE sw_rad_skin
     MODULE PROCEDURE sw_ohlmann_v1
  END INTERFACE
CONTAINS
  ! ------------------------------------------------------
  SUBROUTINE rhocoef(t, s, rhoref, alpha, beta)
    ! ------------------------------------------------------

    !  compute thermal expansion coefficient (alpha) 
    !  and saline contraction coefficient (beta) using 
    !  the international equation of state of sea water 
    !  (1980). Ref: pond and pickard, introduction to 
    !  dynamical oceanography, pp310.  
    !  note: compression effects are not included

    IMPLICIT NONE
    REAL(kind=kind_phys), INTENT(in)  :: t, s, rhoref 
    REAL(kind=kind_phys), INTENT(out) :: alpha, beta  
    REAL(kind=kind_phys) :: tc

    tc = t - t0K

    alpha =                                                        & 
         6.793952e-2                                              & 
         - 2.0 * 9.095290e-3 * tc     +  3.0 * 1.001685e-4 * tc**2  & 
         - 4.0 * 1.120083e-6 * tc**3  +  5.0 * 6.536332e-9 * tc**4  & 
         - 4.0899e-3 * s                                            & 
         + 2.0 * 7.6438e-5 * tc * s  -  3.0 * 8.2467e-7 * tc**2 * s & 
         + 4.0 * 5.3875e-9 * tc**3 * s                              & 
         + 1.0227e-4 * s**1.5 -  2.0 * 1.6546e-6 * tc * s**1.5

    ! NOTE: rhoref - specify 
    !
    alpha =  -alpha/rhoref

    beta  =                                             &
         8.24493e-1          -  4.0899e-3 * tc           &
         + 7.6438e-5 * tc**2 -  8.2467e-7 * tc**3        &
         + 5.3875e-9 * tc**4 -  1.5 * 5.72466e-3 * s**.5 &
         + 1.5 * 1.0227e-4 * tc * s**.5                  &
         -  1.5 * 1.6546e-6 * tc**2 * s**.5              &
         + 2.0 * 4.8314e-4 * s

    beta = beta / rhoref

  END SUBROUTINE rhocoef
  ! ----------------------------------------
  SUBROUTINE density(t, s, rho)
    ! ----------------------------------------
    IMPLICIT NONE

    ! input
    REAL(kind=kind_phys), INTENT(in)  :: t     !unit, K
    REAL(kind=kind_phys), INTENT(in)  :: s     !unit, 1/1000
    ! output
    REAL(kind=kind_phys), INTENT(out) :: rho   !unit, kg/m^3 
    ! local
    REAL(kind=kind_phys) :: tc

    ! compute density using the international equation 
    ! of state of sea water 1980, (pond and pickard, 
    ! introduction to dynamical oceanography, pp310). 
    ! compression effects are not included

    rho = 0.0
    tc = t - t0K

    !  effect of temperature on density (lines 1-3)
    !  effect of temperature and salinity on density (lines 4-8)
    rho = &
         999.842594                 +  6.793952e-2 * tc     &
         - 9.095290e-3 * tc**2        +  1.001685e-4 * tc**3     &
         - 1.120083e-6 * tc**4        +  6.536332e-9 * tc**5     &
         + 8.24493e-1 * s          -  4.0899e-3 * tc * s         &
         + 7.6438e-5 * tc**2 * s   -  8.2467e-7 * tc**3 * s      &
         + 5.3875e-9 * tc**4 * s   -  5.72466e-3 * s**1.5        &
         + 1.0227e-4 * tc * s**1.5 -  1.6546e-6 * tc**2 * s**1.5 &
         + 4.8314e-4 * s**2

  END SUBROUTINE density
  !
  !======================
  !
  elemental subroutine sw_ps_9b(z,fxp)
    !
    ! Fraction of the Solar radiation absorbed by the ocean at the depth z 
    ! following Paulson and Simpson, 1981
    !
    ! INPUT:
    ! z:       depth (m)
    !
    ! OUTPUT:
    ! fxp: Fraction of the solar radiation absorbed by the ocean at depth z (W/m^2)
    !
    IMPLICIT NONE
    REAL,INTENT(in):: z
    REAL,INTENT(out):: fxp
    REAL, DIMENSION(9), PARAMETER :: F=(/0.237,0.36,0.179,0.087,0.08,0.0246,0.025,0.007,0.0004/) &
                                ,gamma=(/34.8,2.27,3.15e-2,5.48e-3,8.32e-4,1.26e-4,3.13e-4,7.82e-5,1.44e-5/)
    !
    IF(z>0) THEN
      fxp=1.0-(F(1)*exp(-z/gamma(1))+F(2)*exp(-z/gamma(2))+F(3)*exp(-z/gamma(3))+ &
               F(4)*exp(-z/gamma(4))+F(5)*exp(-z/gamma(5))+F(6)*exp(-z/gamma(6))+ &
               F(7)*exp(-z/gamma(7))+F(8)*exp(-z/gamma(8))+F(9)*exp(-z/gamma(9)))
    ELSE
       fxp=0.
    ENDIF
    !
  end subroutine sw_ps_9b
  !
  !======================
  !
  !
  !======================
  !
  elemental subroutine sw_ps_9b_Aw(z,Aw)
    !
    ! d(fw)/d(z) for 9-band 
    !
    ! INPUT:
    ! z:       depth (m)
    !
    ! OUTPUT:
    ! fxp: Fraction of the solar radiation absorbed by the ocean at depth z (W/m^2)
    !
    IMPLICIT NONE
    REAL,INTENT(in):: z
    REAL,INTENT(out):: Aw
    REAL, DIMENSION(9), PARAMETER :: F=(/0.237,0.36,0.179,0.087,0.08,0.0246,0.025,0.007,0.0004/) &
                                ,gamma=(/34.8,2.27,3.15e-2,5.48e-3,8.32e-4,1.26e-4,3.13e-4,7.82e-5,1.44e-5/)
    !
    IF(z>0) THEN
      Aw=(F(1)/gamma(1))*exp(-z/gamma(1))+(F(2)/gamma(2))*exp(-z/gamma(2))+(F(3)/gamma(3))*exp(-z/gamma(3))+ &
         (F(1)/gamma(4))*exp(-z/gamma(4))+(F(2)/gamma(5))*exp(-z/gamma(5))+(F(6)/gamma(6))*exp(-z/gamma(6))+ &
         (F(1)/gamma(7))*exp(-z/gamma(7))+(F(2)/gamma(8))*exp(-z/gamma(8))+(F(9)/gamma(9))*exp(-z/gamma(9))
    ELSE
       Aw=0.
    ENDIF
    !
  end subroutine sw_ps_9b_Aw
  !
  !======================
  elemental SUBROUTINE sw_fairall_6exp_v1(z,fxp)
    !
    ! Fraction of the Solar radiation absorbed by the ocean at the depth z (Fairall et all, 1996, p. 1298)
    ! following Paulson and Simpson, 1981
    !
    ! INPUT:
    ! z:       depth (m)
    !
    ! OUTPUT:
    ! fxp: Fraction of the solar radiation absorbed by the ocean at depth z (W/m^2)
    !
    IMPLICIT NONE
    REAL(kind=kind_phys),INTENT(in):: z
    REAL(kind=kind_phys),INTENT(out):: fxp
    REAL(kind=kind_phys), DIMENSION(9), PARAMETER :: F=(/0.237,0.36,0.179,0.087,0.08,0.0246,0.025,0.007,0.0004/) &
         ,gamma=(/34.8,2.27,3.15e-2,5.48e-3,8.32e-4,1.26e-4,3.13e-4,7.82e-5,1.44e-5/)
    REAL(kind=kind_phys),DIMENSION(9) :: zgamma
    REAL(kind=kind_phys),DIMENSION(9) :: f_c
    !
    IF(z>0) THEN
       zgamma=z/gamma
       f_c=F*(1.-1./zgamma*(1-EXP(-zgamma)))
       fxp=SUM(f_c)
    ELSE
       fxp=0.
    ENDIF
    !
  END SUBROUTINE sw_fairall_6exp_v1
  !
  !======================
  !
  !
  elemental SUBROUTINE sw_fairall_6exp_v1_Aw(z,Aw)
    !
    ! Fraction of the Solar radiation absorbed by the ocean at the depth z (Fairall et all, 1996, p. 1298)
    ! following Paulson and Simpson, 1981
    !
    ! INPUT:
    ! z:       depth (m)
    !
    ! OUTPUT:
    ! Aw: d(fxp)/d(z)
    !
    ! fxp: Fraction of the solar radiation absorbed by the ocean at depth z (W/m^2)
    !
    IMPLICIT NONE
    REAL(kind=kind_phys),INTENT(in):: z
    REAL(kind=kind_phys),INTENT(out):: Aw
    REAL(kind=kind_phys) :: fxp
    REAL(kind=kind_phys), DIMENSION(9), PARAMETER :: F=(/0.237,0.36,0.179,0.087,0.08,0.0246,0.025,0.007,0.0004/) &
         ,gamma=(/34.8,2.27,3.15e-2,5.48e-3,8.32e-4,1.26e-4,3.13e-4,7.82e-5,1.44e-5/)
    REAL(kind=kind_phys),DIMENSION(9) :: zgamma
    REAL(kind=kind_phys),DIMENSION(9) :: f_Aw
    !
    IF(z>0) THEN
       zgamma=z/gamma
       f_Aw=(F/z)*((gamma/z)*(1-EXP(-zgamma))-EXP(-zgamma))
       Aw=SUM(f_Aw)

!      write(*,'(a,F6.2,F12.6,9F10.4)') 'z,Aw in sw_rad_Aw : ',z,Aw,f_Aw

    ELSE
       Aw=0.
    ENDIF
    !
  END SUBROUTINE sw_fairall_6exp_v1_Aw
  !
  elemental SUBROUTINE sw_fairall_6exp_v1_sum(z,sum)
    !
    ! Fraction of the Solar radiation absorbed by the ocean at the depth z (Fairall et all, 1996, p. 1298)
    ! following Paulson and Simpson, 1981
    !
    ! INPUT:
    ! z:       depth (m)
    !
    ! OUTPUT:
    ! sum: for convection depth calculation
    !
    !
    IMPLICIT NONE
    REAL(kind=kind_phys),INTENT(in):: z
    REAL(kind=kind_phys),INTENT(out):: sum
    REAL(kind=kind_phys), DIMENSION(9), PARAMETER :: gamma=(/34.8,2.27,3.15e-2,5.48e-3,8.32e-4,1.26e-4,3.13e-4,7.82e-5,1.44e-5/)
    REAL(kind=kind_phys),DIMENSION(9) :: zgamma
    REAL(kind=kind_phys),DIMENSION(9) :: f_sum
    !
!    zgamma=z/gamma
!    f_sum=(zgamma/z)*EXP(-zgamma)
!    sum=SUM(f_sum)

    sum=(1.0/gamma(1))*exp(-z/gamma(1))+(1.0/gamma(2))*exp(-z/gamma(2))+(1.0/gamma(3))*exp(-z/gamma(3))+ &
        (1.0/gamma(4))*exp(-z/gamma(4))+(1.0/gamma(5))*exp(-z/gamma(5))+(1.0/gamma(6))*exp(-z/gamma(6))+ &
        (1.0/gamma(7))*exp(-z/gamma(7))+(1.0/gamma(8))*exp(-z/gamma(8))+(1.0/gamma(9))*exp(-z/gamma(9))
    !
  END SUBROUTINE sw_fairall_6exp_v1_sum
  !
  !======================

  elemental SUBROUTINE sw_fairall_simple_v1(F_sol_0,z,dF_sol_z)
    !
    ! Solar radiation absorbed by the ocean at the depth z (Fairall et all, 1996, p. 1298)
    ! 
    ! INPUT: 
    ! F_sol_0: solar radiation at the ocean surface (W/m^2)
    ! z:       depth (m)
    !
    ! OUTPUT:
    ! dF_sol_z: solar radiation absorbed by the ocean at depth z (W/m^2)
    !
    IMPLICIT NONE
    REAL(kind=kind_phys),INTENT(in):: z,F_sol_0
    REAL(kind=kind_phys),INTENT(out):: dF_sol_z
    !
    IF(z>0) THEN
       dF_sol_z=F_sol_0*(0.137+11.0*z-6.6e-6/z*(1.-EXP(-z/8.e-4)))
    ELSE
       dF_sol_z=0.
    ENDIF
    !
  END SUBROUTINE sw_fairall_simple_v1
  !
  !======================
  !
  elemental SUBROUTINE sw_wick_v1(F_sol_0,z,dF_sol_z)
    !
    ! Solar radiation absorbed by the ocean at the depth z (Zeng and Beljaars, 2005, p.5)
    ! 
    ! INPUT: 
    ! F_sol_0: solar radiation at the ocean surface (W/m^2)
    ! z:       depth (m)
    !
    ! OUTPUT:
    ! dF_sol_z: solar radiation absorbed by the ocean at depth z (W/m^2)
    !
    IMPLICIT NONE
    REAL(kind=kind_phys),INTENT(in):: z,F_sol_0
    REAL(kind=kind_phys),INTENT(out):: dF_sol_z
    !
    IF(z>0) THEN
       dF_sol_z=F_sol_0*(0.065+11.0*z-6.6e-5/z*(1.-EXP(-z/8.e-4)))
    ELSE
       dF_sol_z=0.
    ENDIF
    !
  END SUBROUTINE sw_wick_v1
  !
  !======================
  !
  elemental SUBROUTINE sw_soloviev_3exp_v1(F_sol_0,z,dF_sol_z)
    !
    ! Solar radiation absorbed by the ocean at the depth z (Fairall et all, 1996, p. 1301)
    ! following Soloviev, 1982
    ! 
    ! INPUT: 
    ! F_sol_0: solar radiation at the ocean surface (W/m^2)
    ! z:       depth (m)
    !
    ! OUTPUT:
    ! dF_sol_z: solar radiation absorbed by the ocean at depth z (W/m^2)
    !
    IMPLICIT NONE
    REAL(kind=kind_phys),INTENT(in):: z,F_sol_0
    REAL(kind=kind_phys),INTENT(out):: dF_sol_z
    REAL(kind=kind_phys),DIMENSION(3) :: f_c
    REAL(kind=kind_phys), DIMENSION(3), PARAMETER :: f=(/0.45,0.27,0.28/) &
         ,gamma=(/12.8,0.357,0.014/)
    !
    IF(z>0) THEN
       f_c=f*gamma(1-EXP(-z/gamma))
       dF_sol_z=F_sol_0*(1.0-SUM(f_c)/z)
    ELSE
       dF_sol_z=0.
    ENDIF
    !
  END SUBROUTINE sw_soloviev_3exp_v1
  !
  !======================
  !
  elemental SUBROUTINE sw_soloviev_3exp_v2(F_sol_0,z,dF_sol_z)
    !
    ! Solar radiation absorbed by the ocean at the depth z (Fairall et all, 1996, p. 1301)
    ! following Soloviev, 1982
    ! 
    ! INPUT: 
    ! F_sol_0: solar radiation at the ocean surface (W/m^2)
    ! z:       depth (m)
    !
    ! OUTPUT:
    ! dF_sol_z: solar radiation absorbed by the ocean at depth z (W/m^2)
    !
    IMPLICIT NONE
    REAL(kind=kind_phys),INTENT(in):: z,F_sol_0
    REAL(kind=kind_phys),INTENT(out):: dF_sol_z
    !
    IF(z>0) THEN
       dF_sol_z=F_sol_0*(1.0 &
            -(0.28*0.014*(1.-exp(-z/0.014)) &
            +0.27*0.357*(1.-exp(-z/0.357)) &        
            +.45*12.82*(1.-exp(-z/12.82)))/z &
            )
    ELSE
       dF_sol_z=0.
    ENDIF
    !
  END SUBROUTINE sw_soloviev_3exp_v2

  elemental SUBROUTINE sw_soloviev_3exp_v2_Aw(z,Aw)
    !
    ! Aw = d(fxp)/d(z)
    ! following Soloviev, 1982
    !
    ! INPUT:
    ! z:       depth (m)
    !
    ! OUTPUT:
    ! Aw: d(fxp)/d(z)
    !
    IMPLICIT NONE
    REAL(kind=kind_phys),INTENT(in):: z
    REAL(kind=kind_phys),INTENT(out):: Aw
    REAL(kind=kind_phys):: fxp
    !
    IF(z>0) THEN
       fxp=(1.0 &
            -(0.28*0.014*(1.-exp(-z/0.014)) &
            + 0.27*0.357*(1.-exp(-z/0.357)) &
            + 0.45*12.82*(1.-exp(-z/12.82)))/z &
            )
       Aw=1.0-fxp-(0.28*exp(-z/0.014)+0.27*exp(-z/0.357)+0.45*exp(-z/12.82))
    ELSE
       Aw=0.
    ENDIF
  END SUBROUTINE sw_soloviev_3exp_v2_Aw
  !
  !
  !======================
  !
  elemental SUBROUTINE sw_ohlmann_v1(z,fxp)
    !
    ! Fraction of the Solar radiation absorbed by the ocean at the depth z
    !
    ! INPUT:
    ! z:       depth (m)
    !
    ! OUTPUT:
    ! fxp: Fraction of the solar radiation absorbed by the ocean at depth z (W/m^2)
    !
    IMPLICIT NONE
    REAL(kind=kind_phys),INTENT(in):: z
    REAL(kind=kind_phys),INTENT(out):: fxp
    !
    IF(z>0) THEN
       fxp=.065+11.*z-6.6e-5/z*(1.-EXP(-z/8.0e-4))
    ELSE
       fxp=0.
    ENDIF
    !
  END SUBROUTINE sw_ohlmann_v1
  !

function grv(lat)
  real(kind=kind_phys) :: lat
  real(kind=kind_phys) :: gamma,c1,c2,c3,c4,pi,phi,x
  gamma=9.7803267715
  c1=0.0052790414
  c2=0.0000232718
  c3=0.0000001262
  c4=0.0000000007
  pi=3.141593
                                                                                                                                                             
  phi=lat*pi/180
  x=sin(phi)
  grv=gamma*(1+(c1*x**2)+(c2*x**4)+(c3*x**6)+(c4*x**8))
  !print *,'grav=',grv,lat
end function grv

SUBROUTINE solar_time_from_Julian(jday,xlon,soltim)
  !
  ! Calculate solar time from the Julian date
  !
  IMPLICIT NONE
  REAL(kind=kind_phys), INTENT(in)  :: jday
  REAL(kind=kind_phys), INTENT(in)  :: xlon
  REAL(kind=kind_phys), INTENT(out) :: soltim
  REAL(kind=kind_phys)                            :: fjd,xhr,xmin,xsec,intime
  INTEGER                                        :: nn
  !
  fjd=jday-FLOOR(jday)
  fjd=jday
  xhr=FLOOR(fjd*24.0)-SIGN(12.0,fjd-0.5)
  xmin=NINT(fjd*1440.0)-(xhr+SIGN(12.0,fjd-0.5))*60
  xsec=0
  intime=xhr+xmin/60.0+xsec/3600.0+24.0
  soltim=mod(xlon/15.0+intime,24.0)*3600.0
END SUBROUTINE solar_time_from_Julian

!
!***********************************************************************
!
      subroutine compjd(jyr,jmnth,jday,jhr,jmn,jd,fjd)
!fpp$ noconcur r
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    compjd      computes julian day and fraction
!   prgmmr: kenneth campana  org: w/nmc23    date: 89-07-07
!
! abstract: computes julian day and fraction
!   from year, month, day and time utc.
!
! program history log:
!   77-05-06  ray orzol,gfdl
!   98-05-15  iredell   y2k compliance
!
! usage:    call compjd(jyr,jmnth,jday,jhr,jmn,jd,fjd)
!   input argument list:
!     jyr      - year (4 digits)
!     jmnth    - month
!     jday     - day
!     jhr      - hour
!     jmn      - minutes 
!   output argument list:
!     jd       - julian day.
!     fjd      - fraction of the julian day.
!
! subprograms called:
!   iw3jdn     compute julian day number
!
! attributes:
!   language: fortran.
!
!$$$
      use machine , only :kind_phys
      implicit none
!
      integer jyr,jmnth,jday,jhr,jmn,jd
      integer iw3jdn
      real (kind=kind_phys) fjd
      jd=iw3jdn(jyr,jmnth,jday)
      if(jhr.lt.12) then
        jd=jd-1
        fjd=0.5+jhr/24.+jmn/1440.
      else
        fjd=(jhr-12)/24.+jmn/1440.
      endif
      end subroutine compjd

END MODULE module_nst_water_prop
