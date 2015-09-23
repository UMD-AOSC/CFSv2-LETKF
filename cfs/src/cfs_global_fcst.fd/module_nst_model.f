module nst_module

!
! The module of diurnal thermocline layer model 
!
 USE MACHINE , ONLY : kind_phys
 use module_nst_parameters, only: z_w_max,z_w_min,z_w_ini,eps_z_w,eps_conv,    &
                                  eps_sfs,niter_z_w,niter_conv,niter_sfs,Ri_c, &
                                  Ri_g,omg_m,omg_sh, kw => tc_w,visw,t0K,cp_w, &
                                  z_c_max,z_c_ini,ustar_a_min,delz,exp_const,  &
                                  rad2deg,const_rot,tw_max,sst_max
 USE module_nst_water_prop, ONLY: sw_rad_skin,sw_ps_9b,sw_ps_9b_Aw
 implicit none

 contains

 subroutine dtm_1p(kdt,timestep,rich,tox,toy,I0,Q,sss,sep,Q_Ts,Hl_Ts,rho,      &
                   alpha,beta,alon,sinlat,soltim,grav,Le,d_conv,               &
                   xt,xs,xu,xv,xz,xzts,xtts)

   integer, intent(in) :: kdt
   real(kind=kind_phys), intent(in) :: timestep,rich,tox,toy,I0,Q,sss,sep,Q_Ts,&
                                       Hl_Ts,rho,alpha,beta,alon,sinlat,soltim,&
                                       grav,Le,d_conv
   real(kind=kind_phys), intent(inout) :: xt,xs,xu,xv,xz,xzts,xtts
! local variables
   real(kind=kind_phys) :: xt0,xs0,xu0,xv0,xz0,xzts0,xtts0
   real(kind=kind_phys) :: fac,t0,tau

!
! input variables
!
! timestep:       integration time step in seconds
! rich    :       critical Ri (flow dependent)
! tox     :       x wind stress                       (N*M^-2 or KG/M/S^2)
! toy     :       y wind stress                       (N*M^-2 or KG/M/S^2)
! I0      :       solar radiation flux at the surface (WM^-2)
! Q       :       non-solar heat flux at the surface  (WM^-2)
! sss     :       Salinity                            (ppt)
! sep     :       Sr(E-P)                             (ppt*M/S)
! Q_Ts    :       d(Q)/d(Ts) : Q = the sum of non-solar heat fluxes
! Hl_Ts   :       d(Hl)/d(Ts)
! rho     :       sea water density                   (KG*M^-3)
! alpha   :       thermal expansion coefficient       (1/K)
! beta    :       saline contraction coefficient      (1/ppt)
! sinlat  :       sine (lat)
! grav    :       gravity accelleration
! Le      :       Le=(2.501-.00237*tsea)*1e6
! d-conv  :       FCL thickness
!
! inout variables
!
! xt      :       DTL heat content            (M*K)
! xs      :       DTL salinity content        (M*ppt)
! xu      :       DTL x current content       (M*M/S)
! xv      :       DTL y current content       (M*M/S)
! xz      :       DTL thickness               (M)
! xzts    :       d(xz)/d(ts)                 (M/K )
! xtts    :       d(xt)/d(ts)                 (M)
!
! logical lprnt

! if (lprnt) print *,' first xt=',xt
  if ( xt <= 0.0 ) then                 ! dtl doesn't exist yet
    call dtm_onset(kdt,timestep,rich,tox,toy,I0,Q,sss,sep,Q_Ts,Hl_Ts,rho,alpha,&
                   beta,alon,sinlat,soltim,grav,Le,xt,xs,xu,xv,xz,xzts,xtts)
  elseif ( xt > 0.0 ) then              ! dtl already exists
!
! forward the system one time step 
!
    call Eulerm(kdt,timestep,rich,tox,toy,I0,Q,sss,sep,Q_Ts,Hl_Ts,rho,alpha,   &
                beta,alon,sinlat,soltim,grav,Le,d_conv,                        &
                xt,xs,xu,xv,xz,xzts,xtts)
  endif                         ! if ( xt == 0 ) then

 end subroutine dtm_1p

 subroutine Eulerm(kdt,timestep,rich,tox,toy,I0,Q,sss,sep,Q_Ts,Hl_Ts,rho,alpha,&
                   beta,alon,sinlat,soltim,grav,Le,d_conv,                     &
                   xt,xs,xu,xv,xz,xzts,xtts)

!
! subroutine Eulerm: integrate one time step with Modified Euler method
!
   integer, intent(in) :: kdt
   real(kind=kind_phys), intent(in) :: timestep,rich,tox,toy,I0,Q,sss,sep,Q_Ts,&
                                       Hl_Ts,rho,alpha,beta,alon,sinlat,soltim,&
                                       grav,Le,d_conv
   real(kind=kind_phys), intent(inout) :: xt,xs,xu,xv,xz,xzts,xtts
!  local variables
   real(kind=kind_phys) :: xt0,xs0,xu0,xv0,xz0,xzts0,xtts0
   real(kind=kind_phys) :: fw,Aw,Q_warm
   real(kind=kind_phys) :: xt1,xs1,xu1,xv1,xz1,xzts1,xtts1
   real(kind=kind_phys) :: xt2,xs2,xu2,xv2,xz2,xzts2,xtts2
   real(kind=kind_phys) :: dzw,drho,fc,t_fcl,ttop,dz
   real(kind=kind_phys) :: alat,t0,te,tau,speed
!  logical lprnt

!
! input variables
!
! timestep:       integration time step in seconds
! rich    :       critial Ri (flow/mass dependent)
! tox     :       x wind stress                       (N*M^-2 or KG/M/S^2)
! toy     :       y wind stress                       (N*M^-2 or KG/M/S^2)
! I0      :       solar radiation flux at the surface (WM^-2)
! Q       :       non-solar heat flux at the surface  (WM^-2)
! sss     :       Salinity                            (ppt)
! sep     :       Sr(E-P)                             (ppt*M/S)
! Q_Ts    :       d(Q)/d(Ts) : Q = the sum of non-solar heat fluxes
! Hl_Ts   :       d(Hl)/d(Ts)
! rho     :       sea water density                   (KG*M^-3)
! alpha   :       thermal expansion coefficient       (1/K)
! beta    :       saline contraction coefficient      (1/ppt)
! alon    :       longitude (DEG)
! sinlat  :       sine (lat)
! soltim  :       solar time
! grav    :       gravity accelleration
! Le      :       Le=(2.501-.00237*tsea)*1e6
! d_conv  :       FCL thickness                       (M)
!
! inout variables
!
! xt      :       DTL heat content                    (M*K)
! xs      :       DTL salinity content                (M*ppt)
! xu      :       DTL x current content               (M*M/S)
! xv      :       DTL y current content               (M*M/S)
! xz      :       DTL thickness                       (M)
! xzts    :       d(xz)/d(ts)                         (M/K )
! xtts    :       d(xt)/d(ts)                         (M)

  xt0   = xt
  xs0   = xs
  xu0   = xu
  xv0   = xv
  xz0   = xz
  xtts0 = xtts
  xzts0 = xzts
  speed = max(1.0e-8,xu0*xu0 + xv0*xv0)
  
  alat=asin(sinlat)*rad2deg

  fc   = const_rot*sinlat
  call sw_ps_9b(xz0,fw)
  q_warm = fw*I0-Q                                !total heat abs in warm layer
  call sw_ps_9b_Aw(xz0,Aw)
  drho  = -alpha*q_warm/(rho*cp_w) + omg_m*beta*sep
! dzw   = xz0*(tox*xu0+toy*xv0) / (rho*(xu0*xu0+xv0*xv0))                 &
!       + xz0*xz0*xz0*drho*grav / (4.0*rich*(xu0*xu0+xv0*xv0))
  dzw   = xz0 * ((tox*xu0+toy*xv0) / (rho*speed)                          &
              +   xz0*xz0*drho*grav / (4.0*rich*speed))

  xt1   = xt0   + timestep*q_warm/(rho*cp_w)
  xs1   = xs0   + timestep*sep
  xu1   = xu0   + timestep*(fc*xv0+tox/rho)
  xv1   = xv0   + timestep*(-fc*xu0+toy/rho)
  xz1   = xz0   + timestep*dzw

! if (lprnt) print *,' xt1=',xt1,' xz1=',xz1,' xz0=',xz0,' dzw=',dzw,     &
! 'timestep=',timestep,tox,toy,xu0,xv0,rho,drho,grav,rich

  if ( xt1 <= 0.0 .or. xz1 <= 0.0 .or. xz1 > z_w_max ) then
    call dtl_reset(xt,xs,xu,xv,xz,xzts,xtts)
    go to 100
  endif

! call dtm_1p_zwa(kdt,timestep,I0,Q,rho,d_conv,xt1,xs1,xu1,xv1,xz1,tr_mda,tr_fca,tr_tla,tr_mwa)

  xzts1 = xzts0 + timestep*((1.0/(xu0*xu0+xv0*xv0)) *                          &
         ( (alpha*Q_Ts/cp_w+omg_m*beta*sss*Hl_Ts/Le)*grav*xz0**3/(4.0*rich*rho)&
        +( (tox*xu0+toy*xv0)/rho+(3.0*drho-alpha*I0*Aw*xz0/(rho*cp_w))         &
                                         *grav*xz0*xz0/(4.0*rich) )*xzts0 ))
  xtts1 = xtts0 + timestep*(I0*Aw*xzts0-Q_Ts)/(rho*cp_w)
           
! if ( 2.0*xt1/xz1 > 0.001 ) then
! write(*,'(a,I5,2F8.3,4F8.2,F10.6,10F8.4)') 'Eulerm_01 : ',kdt,alat,alon,soltim/3600.,I0,Q,Q_warm,sep,&
!          2.0*xt1/xz1,2.0*xs1/xz1,2.0*xu1/xz1,2.0*xv1/xz1,xz1,xtts1,xzts1,d_conv,t_fcl,te
! endif

  call sw_ps_9b(xz1,fw)
  q_warm = fw*I0-Q                                !total heat abs in warm layer
  call sw_ps_9b_Aw(xz1,Aw)
  drho = -alpha*q_warm/(rho*cp_w) + omg_m*beta*sep
  dzw = xz1*(tox*xu1+toy*xv1) / (rho*(xu1*xu1+xv1*xv1))                      &
      + xz1*xz1*xz1*drho*grav / (4.0*rich*(xu1*xu1+xv1*xv1))

  xt2   = xt0   + timestep*q_warm/(rho*cp_w)
  xs2   = xs0   + timestep*sep
  xu2   = xu0   + timestep*(fc*xv1+tox/rho)
  xv2   = xv0   + timestep*(-fc*xu1+toy/rho)
  xz2   = xz0   + timestep*dzw

! if (lprnt) print *,' xt2=',xt2,' xz2=',xz2

  if ( xt2 <= 0.0 .or. xz2 <= 0.0 .or. xz2 > z_w_max ) then
    call dtl_reset(xt,xs,xu,xv,xz,xzts,xtts)
    go to 100
  endif

  xzts2 = xzts0 + timestep*((1.0/(xu1*xu1+xv1*xv1)) *                          &
         ( (alpha*Q_Ts/cp_w+omg_m*beta*sss*Hl_Ts/Le)*grav*xz1**3/(4.0*rich*rho)&
        +( (tox*xu1+toy*xv1)/rho+(3.0*drho-alpha*I0*Aw*xz1/(rho*cp_w))*        &
                                            grav*xz1*xz1/(4.0*rich) )*xzts1 ))
  xtts2 = xtts0 + timestep*(I0*Aw*xzts1-Q_Ts)/(rho*cp_w)

  xt   = 0.5*(xt1 + xt2)
  xs   = 0.5*(xs1 + xs2)
  xu   = 0.5*(xu1 + xu2)
  xv   = 0.5*(xv1 + xv2)
  xz   = 0.5*(xz1 + xz2)
  xzts = 0.5*(xzts1 + xzts2)
  xtts = 0.5*(xtts1 + xtts2)

  if ( xt <= 0.0 .or. xz < 0.0 .or. xz > z_w_max ) then
    call dtl_reset(xt,xs,xu,xv,xz,xzts,xtts)
  endif

! if (lprnt) print *,' xt=',xt,' xz=',xz
! if ( 2.0*xt/xz > 0.001 ) then
! write(*,'(a,I5,2F8.3,4F8.2,F10.6,10F8.4)') 'Eulerm_02 : ',kdt,alat,alon,soltim/3600.,I0,Q,Q_warm,sep,&
!          2.0*xt/xz,2.0*xs/xz,2.0*xu/xz,2.0*xv/xz,xz,xtts,xzts,d_conv,t_fcl,te
! endif
 100 continue

 end subroutine Eulerm

 subroutine dtm_1p_zwa(kdt,timestep,I0,Q,rho,d_conv,xt,xs,xu,xv,xz,tr_mda,tr_fca,tr_tla,tr_mwa)
!  apply xz adjustment:  Minimum Depth Adjustment (MDA)
!                        Free Convection Adjustment (FCA);
!                        Top Layer Adjustment (TLA);
!                        Maximum Warming Adjustment (MWA)
!   
   integer, intent(in) :: kdt
   real(kind=kind_phys), intent(in)    :: timestep,I0,Q,rho,d_conv
   real(kind=kind_phys), intent(inout) :: xt,xs,xu,xv,xz
   real(kind=kind_phys), intent(out)   :: tr_mda,tr_fca,tr_tla,tr_mwa
!  local variables
   real(kind=kind_phys) :: dz,t0,ttop0,ttop,fw,Q_warm
   real(kind=kind_phys) :: xz_fca,xz_tla,xz_mwa
!
   real(kind=kind_phys) xz_mda

   tr_mda = 0.0; tr_fca = 0.0; tr_tla = 0.0; tr_mwa = 0.0

!  Apply MDA
   if ( z_w_min > xz ) then
     xz_mda  = z_w_min
   endif
!  Apply FCA
   if ( d_conv > 0.0 ) then
     xz_fca = 2.0*xt/((2.0*xt/xz)*(1.0-d_conv/(2.0*xz)))
     tr_fca = 1.0 
     if ( xz_fca >= z_w_max ) then
       call dtl_reset_cv(xt,xs,xu,xv,xz)
       go to 10
     endif
   endif
!  Apply TLA
   dz = min(xz,max(d_conv,delz))
   call sw_ps_9b(dz,fw)
   Q_warm=fw*I0-Q                                !total heat abs in warm layer

   if ( Q_warm > 0.0 ) then
     call cal_ttop(kdt,timestep,Q_warm,rho,dz,xt,xz,ttop0)
!    ttop = (2.0*xt/xz)*(1.0-dz/(2.0*xz))
     ttop = ((xt+xt)/xz)*(1.0-dz/(xz+xz))
     if ( ttop > ttop0 ) then
       xz_tla = (xt+sqrt(xt*(xt-delz*ttop0)))/ttop0
       tr_tla = 1.0 
       if ( xz_tla >= z_w_max ) then
         call dtl_reset_cv(xt,xs,xu,xv,xz)
         go to 10
       endif
     endif
   endif

!  Apply MWA
   t0 = 2.0*xt/xz
   if ( t0 > tw_max ) then
     if ( xz >= z_w_max ) then
       call dtl_reset_cv(xt,xs,xu,xv,xz)
       go to 10
     endif
   endif

   xz = max(xz_mda,xz_fca,xz_tla,xz_mwa)

 10 continue
   
 end subroutine dtm_1p_zwa

 subroutine dtm_1p_fca(d_conv,xt,xtts,xz,xzts)

!  apply xz adjustment:  Free Convection Adjustment (FCA);
!   
   real(kind=kind_phys), intent(in)    :: d_conv,xt,xtts
   real(kind=kind_phys), intent(inout) :: xz,xzts
!  local variables
   real(kind=kind_phys) :: t_fcl,t0,xz0,dxz
!
  t0 = 2.0*xt/xz
  t_fcl = t0*(1.0-d_conv/(2.0*xz))
  xz   = 2.0*xt/t_fcl
! xzts = 2.0*xtts/t_fcl

 end subroutine dtm_1p_fca

 subroutine dtm_1p_tla(dz,te,xt,xtts,xz,xzts)

!  apply xz adjustment: Top Layer Adjustment (TLA);
! 
   real(kind=kind_phys), intent(in)    :: dz,te,xt,xtts
   real(kind=kind_phys), intent(inout) :: xz,xzts
!  local variables
   real(kind=kind_phys) tem
!
   tem = xt*(xt-dz*te)
   if (tem > 0.0) then
     xz = (xt+sqrt(xt*(xt-dz*te)))/te
   else
     xz = z_w_max
   endif
!  xzts = xtts*(1.0+0.5*(2.0*xt-dz*te)/sqrt(xt*(xt-dz*te)))/te
 end subroutine dtm_1p_tla

 subroutine dtm_1p_mwa(xt,xtts,xz,xzts)

!  apply xz adjustment: Maximum Warming Adjustment (MWA)
!
   real(kind=kind_phys), intent(in)    :: xt,xtts
   real(kind=kind_phys), intent(inout) :: xz,xzts
!  local variables
!
   xz   = 2.0*xt/tw_max
!  xzts = 2.0*xtts/tw_max
 end subroutine dtm_1p_mwa

 subroutine dtm_1p_mda(xt,xtts,xz,xzts)

!  apply xz adjustment: Minimum Depth Adjustment (MDA)
!
   real(kind=kind_phys), intent(in)    :: xt,xtts
   real(kind=kind_phys), intent(inout) :: xz,xzts
!  local variables
   real(kind=kind_phys) :: ta
!
   xz   = max(z_w_min,xz)
   ta   = 2.0*xt/xz
!  xzts = 2.0*xtts/ta

 end subroutine dtm_1p_mda

 subroutine dtm_1p_mta(dta,xt,xtts,xz,xzts)

!  apply xz adjustment: Maximum Temperature Adjustment (MTA)
!
   real(kind=kind_phys), intent(in)    :: dta,xt,xtts
   real(kind=kind_phys), intent(inout) :: xz,xzts
!  local variables
   real(kind=kind_phys) :: ta
!
   ta = max(0.0,2.0*xt/xz-dta)
   if ( ta > 0.0 ) then
     xz = 2.0*xt/ta
   else
     xz = z_w_max
   endif
!  xzts = 2.0*xtts/ta

 end subroutine dtm_1p_mta

subroutine convdepth(kdt,timestep,I0,Q,sss,sep,rho,alpha,beta,xt,xs,xz,d_conv)

!
! calculate depth for convective adjustment
!

   integer, intent(in) :: kdt
   real(kind=kind_phys), intent(in)  :: timestep,I0,Q,sss,sep,rho,alpha,beta
   real(kind=kind_phys), intent(in)  :: xt,xs,xz
   real(kind=kind_phys), intent(out) :: d_conv
   real(kind=kind_phys)              :: t,s,d_conv_ini,d_conv2,fxp,Aw,s1,s2,fac1,fac2
   integer :: n
!
! input variables
!
! timestep:       time step in seconds
! I0      :       solar radiation flux at the surface (WM^-2)
! Q       :       non-solar heat flux at the surface  (WM^-2)
! sss     :       Salinity                            (ppt)
! sep     :       Sr(E-P)                             (ppt*M/S)
! rho     :       sea water density                   (KG*M^-3)
! alpha   :       thermal expansion coefficient       (1/K)
! beta    :       saline contraction coefficient      (1/ppt)
! xt      :       initial heat  content               (K*M)
! xs      :       initial salinity content            (ppt*M)
! xz      :       initial DTL thickness               (M)
!
! output variables
!
! d_conv  :       Free convection depth               (m)

! t       :       initial diurnal warming T           (K)
! s       :       initial diurnal warming S           (ppt)

 n = 0
 t = 2.0*xt/xz
 s = 2.0*xs/xz

 s1 = alpha*rho*t-omg_m*beta*rho*s

 if ( s1 == 0.0 ) then
   d_conv = 0.0
 else

   fac1 = alpha*Q/cp_w+omg_m*beta*rho*sep
   if ( I0 <= 0.0 ) then
       d_conv2=(2.0*xz*timestep/s1)*fac1
     if ( d_conv2 > 0.0 ) then
       d_conv = sqrt(d_conv2)
     else
       d_conv = 0.0
     endif
   elseif ( I0 > 0.0 ) then

     d_conv_ini = 0.0

     iter_conv: do n = 1, niter_conv
       call sw_ps_9b(d_conv_ini,fxp)
       call sw_ps_9b_Aw(d_conv_ini,Aw)
       s2 = alpha*(Q-(fxp-Aw*d_conv_ini)*I0)/cp_w+omg_m*beta*rho*sep
       d_conv2=(2.0*xz*timestep/s1)*s2
       if ( d_conv2 < 0.0 ) then
         d_conv = 0.0
         exit iter_conv
       endif
       d_conv = sqrt(d_conv2)
       if ( abs(d_conv-d_conv_ini) < eps_conv .and. n <= niter_conv ) exit iter_conv
       d_conv_ini = d_conv
     enddo iter_conv
     d_conv = max(0.0,min(d_conv,z_w_max))
   endif        ! if ( I0 <= 0.0 ) then

 endif     ! if ( s1 == 0.0 ) then

!  if ( d_conv > 0.01 ) then
!    write(*,'(a,I4,I3,10F9.3,3F10.6,F10.1,F6.2)') ' d_conv : ',kdt,n,d_conv,d_conv_ini,Q,I0,rho,cp_w,timestep,xt,xs,xz,sep, &
!            s1,s2,d_conv2,Aw
!  endif

 end subroutine convdepth

 subroutine dtm_onset(kdt,timestep,rich,tox,toy,I0,Q,sss,sep,Q_Ts,Hl_Ts,rho, &
                      alpha,beta,alon,sinlat,soltim,grav,Le,xt,xs,xu,xv,xz,xzts,xtts)
!
! determine xz iteratively (starting wit fw = 0.5) and then update the other 6 variables
!
           
   integer,intent(in) :: kdt
   real(kind=kind_phys), intent(in) :: timestep,rich,tox,toy,I0,Q,sss,sep,Q_Ts,&
                                       Hl_Ts,rho,alpha,beta,alon,sinlat,soltim,grav,Le
   real(kind=kind_phys), intent(out) :: xt,xs,xu,xv,xz,xzts,xtts
   real(kind=kind_phys) :: xt0,xs0,xu0,xv0,xz0
   real(kind=kind_phys) :: xt1,xs1,xu1,xv1,xz1
   real(kind=kind_phys) :: fw,Aw,Q_warm,ft0,fs0,fu0,fv0,fz0,ft1,fs1,fu1,fv1,fz1
   real(kind=kind_phys) :: coeff1,coeff2,ftime,z_w,z_w_tmp,fc,warml,alat
   integer :: n
!
! input variables
!
! timestep:       time step in seconds
! tox     :       x wind stress                       (N*M^-2 or KG/M/S^2)
! toy     :       y wind stress                       (N*M^-2 or KG/M/S^2)
! I0      :       solar radiation flux at the surface (WM^-2)
! Q       :       non-solar heat flux at the surface  (WM^-2)
! sss     :       Salinity                            (ppt)
! sep     :       Sr(E-P)                             (ppt*M/S)
! rho     :       sea water density                   (KG*M^-3)
! alpha   :       thermal expansion coefficient       (1/K)
! beta    :       saline contraction coefficient      (1/ppt)
! alon    :       longitude
! sinlat  :       sine(latitude)
! grav    :       gravity accelleration
! Le      :       Le=(2.501-.00237*tsea)*1e6
!
! output variables
!
! xt      :       onset T content in DTL
! xs      :       onset S content in DTL 
! xu      :       onset u content in DTL 
! xv      :       onset v content in DTL 
! xz      :       onset DTL thickness               (M)
! xzts    :       onset d(xz)/d(ts)                 (M/K )
! xtts    :       onset d(xt)/d(ts)                 (M)

  fc=1.46/10000.0/2.0*sinlat
  alat = asin(sinlat)
!
! initializing DTL (just before the onset)
!
 xt0   = 0.0
 xs0   = 0.0
 xu0   = 0.0
 xv0   = 0.0

 z_w_tmp=z_w_ini

 call sw_ps_9b(z_w_tmp,fw)
! fw=0.5                         ! 
 Q_warm=fw*I0-Q                                !total heat abs in warm layer

 if ( abs(alat) > 1.0 ) then
   ftime=sqrt((2.0-2.0*cos(fc*timestep))/(fc*fc*timestep))
 else
   ftime=timestep
 endif

 coeff1=alpha*grav/cp_w
 coeff2=omg_m*beta*grav*rho
 warml = coeff1*Q_warm-coeff2*sep

 if ( warml > 0.0 .and. Q_warm > 0.0) then
   iters_z_w: do n = 1,niter_z_w
     if ( warml > 0.0 .and. Q_warm > 0.0 ) THEN
       z_w=sqrt(2.0*rich*ftime/rho)*sqrt(tox**2+toy**2)/sqrt(warml)
     else
       z_w = z_w_max
       exit iters_z_w
     endif

!    write(*,'(a,I4,I4,10F9.3,F9.6,F3.0)') ' z_w = ',kdt,n,z_w,z_w_tmp,timestep,Q_warm,Q,I0,fw,tox,toy,sep,warml,omg_m

     if (ABS(z_w - z_w_tmp) < eps_z_w .AND. z_w/=z_w_max .AND. n < niter_z_w) exit iters_z_w
     z_w_tmp=z_w
     call sw_ps_9b(z_w_tmp,fw)
     Q_warm = fw*I0-Q
     warml = coeff1*Q_warm-coeff2*sep
   end do iters_z_w
 else
   z_w=z_w_max
 endif

 xz0 = max(z_w,z_w_min)

!
! Update xt, xs, xu, xv
!
  if ( z_w < z_w_max .and. Q_warm > 0.0) then

    call sw_ps_9b(z_w,fw)
    q_warm=fw*I0-Q                                !total heat abs in warm layer

    ft0 = q_warm/(rho*cp_w)
    fs0 = sep
    fu0 = fc*xv0+tox/rho
    fv0 = -fc*xu0+toy/rho

    xt1 = xt0 + timestep*ft0
    xs1 = xs0 + timestep*fs0
    xu1 = xu0 + timestep*fu0
    xv1 = xv0 + timestep*fv0

    fz0 = xz0*((tox*xu1+toy*xv1)/rho+omg_m*beta*grav*sep*xz0*xz0/(4.0*rich) &
         -alpha*grav*Q_warm*xz0*xz0/(4.0*rich*cp_w*rho))/(xu1*xu1+xv1*xv1)
    xz1 = xz0 + timestep*fz0

    xz1 = max(xz1,z_w_min)

    if ( xt1 < 0.0 .or. xz1 > z_w_max ) then
      call dtl_reset(xt,xs,xu,xv,xz,xtts,xzts)
      go to 100
    endif

    call sw_ps_9b(xz1,fw)
    Q_warm=fw*I0-Q                                !total heat abs in warm layer

    ft1 = Q_warm/(rho*cp_w)
    fs1 = sep
    fu1 = fc*xv1+tox/rho
    fv1 = -fc*xu1+toy/rho

    fz1 = xz1*((tox*xu1+toy*xv1)/rho+omg_m*beta*grav*sep*xz1*xz1/(4.0*rich) &
         -alpha*grav*Q_warm*xz1*xz1/(4.0*rich*cp_w*rho))/(xu1*xu1+xv1*xv1)

    xt = xt0+0.5*timestep*(ft0+ft1)
    xs = xs0+0.5*timestep*(fs0+fs1)
    xu = xu0+0.5*timestep*(fu0+fu1)
    xv = xv0+0.5*timestep*(fv0+fv1)
    xz = xz0+0.5*timestep*(fz0+fz1)

    xz = max(xz,z_w_min)

    call sw_ps_9b_Aw(xz,Aw)

!   xzts = (Q_Ts+(cp_w*omg_m*beta*sss/(Le*alpha))*Hl_Ts)*xz/(I0*xz*Aw+2.0*q_warm-2.0*(rho*cp_w*omg_m*beta*sss/alpha)*(sep/sss))
    xzts = (Q_Ts+omg_m*rho*cp_w*beta*sss*Hl_Ts*xz/(Le*alpha))/(I0*xz*Aw+2.0*q_warm-2.0*omg_m*rho*cp_w*beta*sss*sep/(Le*alpha))
    xtts = timestep*(I0*Aw*xzts-Q_Ts)/(rho*cp_w)
  endif

  if ( xt < 0.0 .or. xz > z_w_max ) then
    call dtl_reset(xt,xs,xu,xv,xz,xtts,xzts)
  endif
 
 100 continue

 end subroutine dtm_onset

 subroutine cal_w(kdt,xz,xt,xzts,xtts,w_0,w_d)
!
! abstract: calculate w_0,w_d
!
! input variables
!
! kdt     :       the number of time step
! xt      :       DTL heat content  
! xz      :       DTL depth         
! xzts    :       d(zw)/d(ts)
! xtts    :       d(xt)/d(ts)
!
! output variables
!
! w_0     :       coefficint 1 to calculate d(Tw)/d(Ts) 
! w_d     :       coefficint 2 to calculate d(Tw)/d(Ts) 

  integer, intent(in) :: kdt
  real(kind=kind_phys), intent(in) :: xz,xt,xzts,xtts
  real(kind=kind_phys), intent(out) :: w_0,w_d

  w_0 = 2.0*(xtts-xt*xzts/xz)/xz
  w_d = (2.0*xt*xzts/xz**2-w_0)/xz

! if ( 2.0*xt/xz > 1.0 ) then
!   write(*,'(a,I4,2F9.3,4F10.4))') ' cal_w : ',kdt,xz,xt,w_0,w_d,xzts,xtts
! endif
 end subroutine cal_w


 subroutine cal_ttop(kdt,timestep,Q_warm,rho,dz,xt,xz,ttop)
!
! abstract: calculate
!
! input variables
!
! kdt      :       the number of record
! timestep :       the number of record
! Q_warm   :       total heat abs in layer dz
! rho      :       sea water density
! dz       :       dz = max(delz,d_conv) top layer thickness defined to adjust xz
! xt       :       heat content in DTL at previous time
! xz       :       DTL thickness at previous time
!
! output variables
!
! ttop     :       the diurnal warming amount at the top layer with thickness of delz

  integer, intent(in) :: kdt
  real(kind=kind_phys), intent(in) :: timestep,Q_warm,rho,dz,xt,xz
  real(kind=kind_phys), intent(out) :: ttop
  real(kind=kind_phys) :: dt_warm,t0

  dt_warm = (xt+xt)/xz
  t0 = dt_warm*(1.0-dz/(xz+xz))
  ttop = t0 + Q_warm*timestep/(rho*cp_w*dz)

 end subroutine cal_ttop

 subroutine app_sfs(kdt,xt,xs,xu,xv,alpha,beta,grav,d_1p,xz)
!
! abstract: adjust dtm-1p DTL thickness by applying Shear Flow Stability with assumed exponetial profile
!
! input variables
!
! kdt     :       the number of record
! xt      :       heat content in DTL
! xs      :       salinity content in DTL
! xu      :       u-current content in DTL
! xv      :       v-current content in DTL
! alpha
! beta
! grav
! d_1p    :       DTL depth before SFS applied  
!
! output variables
!
! xz      :       DTL depth                    

  integer, intent(in) :: kdt
  real(kind=kind_phys), intent(in) :: xt,xs,xu,xv,alpha,beta,grav,d_1p
  real(kind=kind_phys), intent(out) :: xz
! real(kind=kind_phys) :: ze,cc,xz0,L,d_sfs, tem
  real(kind=kind_phys) ::    cc,xz0,L,d_sfs,t_sfs, tem
  real(kind=kind_phys), parameter :: a2 = 0.2294, c2 = 0.3782
  integer :: n

  cc  = Ri_g/(grav*c2)

  tem = alpha*xt - beta*xs
  if (tem > 0.0) then
    d_sfs = sqrt(2.0*cc*(xu*xu+xv*xv)/tem)
  else
    d_sfs = 0.0
  endif

! xz0 = d_1p
! iter_sfs: do n = 1, niter_sfs
!   L = Int_epn(0.0,xz0,0.0,xz0,2)
!   d_sfs = cc*(xu*xu+xv*xv)/((alpha*xt-beta*xs)*L)
!   write(*,'(a,I6,I3,4F9.4))') ' app_sfs_iter : ',kdt,n,cc,L,xz0,d_sfs
!   if ( abs(d_sfs-xz0) < eps_sfs .and. n <= niter_sfs ) exit iter_sfs
!   xz0 = d_sfs
! enddo iter_sfs
  
! ze = a2*d_sfs             ! not used!

  L = Int_epn(0.0,d_sfs,0.0,d_sfs,2)

! t_sfs = xt/L
! xz = (xt+xt) / t_sfs

    xz = L + L

! write(*,'(a,I6,6F9.4))') ' app_sfs : ',kdt,xz0,d_sfs,d_1p,xz,2.0*xt/d_1p,t_sfs
 end subroutine app_sfs

 subroutine cal_tztr(kdt,xt,c_0,c_d,w_0,w_d,zc,zw,z,tztr)
!
! abstract: calculate d(Tz)/d(Ts)
!
! input variables
!
! kdt     :       the number of record
! xt      :       heat content in DTL
! xz      :       DTL depth                           (M)
! c_0     :       coefficint 1 to calculate d(Tc)/d(Ts) 
! c_d     :       coefficint 2 to calculate d(Tc)/d(Ts) 
! w_0     :       coefficint 1 to calculate d(Tw)/d(Ts) 
! w_d     :       coefficint 2 to calculate d(Tw)/d(Ts) 
!
! output variables
!
! tztr     :      d(Tz)/d(Tr) 

  integer, intent(in) :: kdt
  real(kind=kind_phys), intent(in) :: xt,c_0,c_d,w_0,w_d,zc,zw,z
  real(kind=kind_phys), intent(out) :: tztr
  real(kind=kind_phys)              :: tem

  if ( xt > 0.0 ) then
     if ( z <= zc ) then
!      tztr = 1.0/(1.0-w_0+c_0)+z*(w_d-c_d)/(1.0-w_0+c_0)
       tztr = (1.0+z*(w_d-c_d))/(1.0-w_0+c_0)
     elseif ( z > zc .and. z < zw ) then
!      tztr = (1.0+c_0)/(1.0-w_0+c_0)+z*w_d/(1.0-w_0+c_0)
       tztr = (1.0+c_0+z*w_d)/(1.0-w_0+c_0)
     elseif ( z >= zw ) then
       tztr = 1.0
     endif
   elseif ( xt == 0.0 ) then
     if ( z <= zc ) then
!      tztr = 1.0/(1.0+c_0)-z*c_d/(1.0+c_0)
       tztr = (1.0-z*c_d)/(1.0+c_0)
     else
       tztr = 1.0
     endif
   else
     tztr = 1.0
   endif

! write(*,'(a,I4,9F9.4))') ' cal_tztr : ',kdt,xt,c_0,c_d,w_0,w_d,zc,zw,z,tztr
 end subroutine cal_tztr

subroutine cool_skin(ustar_a,F_nsol,F_sol_0,evap,sss,alpha,beta,rho_w,rho_a,Ts,Q_Ts,Hl_Ts,grav,Le,deltaT_c,z_c,c_0,c_d)
!
! Upper ocean cool-skin parameterizaion, Fairall et al, 1996.
!
! INPUT:
! ustar_a : atmosphreic friction velocity at the air-sea interface (m/s)
! F_nsol  : the "nonsolar" part of the surface heat flux (W/m^s)
! F_sol_0 : solar radiation at the ocean surface (W/m^2)
! evap    : latent heat flux (W/M^2)
! sss     : ocean upper mixed layer salinity (ppu)
! alpha   : thermal expansion coefficient
! beta    : saline contraction coefficient
! rho_w   : oceanic density
! rho_a   : atmospheric density
! Ts      : oceanic surface temperature
! Q_Ts    : d(Q)/d(Ts) : Q = the sum of non-solar heat fluxes
! Hl_Ts   : d(Hl)/d(Ts)
! grav    : gravity 
! Le      : 
!
! OUTPUT:
! deltaT_c: cool-skin temperature correction (degrees K)
! z_c     : molecular sublayer (cool-skin) thickness (m)
! c_0     : coefficient1 to calculate d(Tz)/d(Ts)
! c_d     : coefficient2 to calculate d(Tz)/d(Ts)

!
  real(kind=kind_phys), intent(in) :: ustar_a,F_nsol,F_sol_0,evap,sss,alpha,beta,rho_w,rho_a,Ts,Q_Ts,Hl_Ts,grav,Le
  real(kind=kind_phys), intent(out):: deltaT_c,z_c,c_0,c_d
! declare local variables
  real(kind=kind_phys) :: a1,a2,a3,a4,A_c,B_c,zc_ts,bc1,bc2
  real(kind=kind_phys) :: xi,Hb,ustar1_a,bigc,deltaF,fxp
  real(kind=kind_phys) :: tcw,cc1,cc2,cc3,qcol,dtemp,corioli,A_w,dwat,dtmp,alfac,wetc

  tcw = 0.6

  z_c=z_c_ini                 ! initial quess

  ustar1_a=max(ustar_a,ustar_a_min)

  CALL sw_rad_skin(z_c,fxp)
  deltaF=F_sol_0*fxp

  Hb=alpha*(F_nsol-DeltaF)+beta*sss*cp_w*evap/Le
  bigc=16*grav*cp_w*(rho_w*visw)**3/(rho_a*rho_a*kw*kw)

  if ( Hb > 0 ) then
    xi=6./(1+(bigc*Hb/ustar1_a**4)**0.75)**0.3333333
  else
    xi=6.0
  endif
  z_c=min(z_c_max,xi*visw/(SQRT(rho_a/rho_w)*ustar1_a ))

  CALL sw_rad_skin(z_c,fxp)
  deltaF=F_sol_0*fxp
  deltaF=F_nsol - deltaF
  if ( deltaF > 0 ) then
    deltaT_c= deltaF * z_c / kw
  else
    deltaT_c=0.
    z_c=0.
  endif
!
! calculate c_0 & c_d
!
  if ( z_c > 0.0 ) then
    cc1 = 6.0*visw/(tcw*ustar1_a*(rho_a/rho_w)**0.5)
    cc2 = bigc*alpha/max(ustar_a,ustar_a_min)**4
    cc3 = beta*sss*cp_w/(alpha*Le)
    A_c = a2+a3/z_c**2-(a3/(a4*z_c)+a3/z_c**2)*exp(-z_c/a4)

    if ( Hb > 0.0 ) then
      bc1 = z_c**2*(Q_Ts+cc3*Hl_Ts)
      bc2 = z_c**2*F_sol_0*A_c-4.0*(cc1*tcw)**3*(Hb/alpha)**0.25/(cc2**0.75*z_c**2)
      zc_ts = bc1/bc2
!     B_c = z_c**2*(Q_Ts+cc3*Hl_Ts)/(z_c**2*F_sol_0*A_c-4.0*(cc1*tcw)**3*(Hb/alpha)**0.25/(cc2**0.75*z_c**2))     ! d(z_c)/d(Ts)
      B_c  = (Q_Ts+cc3*Hl_Ts)/(F_sol_0*A_c-4.0*(cc1*tcw)**3*(Hb/alpha)**0.25/(cc2**0.75*z_c**4))     ! d(z_c)/d(Ts)
      c_0 = (z_c*Q_Ts+(F_nsol-DeltaF-F_sol_0*A_c*z_c)*B_c)/tcw                    
      c_d = (F_sol_0*A_c*z_c*B_c-Q_Ts)/tcw                                   

    else
      B_c = 0.0
      zc_ts = 0.0
      c_0 = z_c*Q_Ts/tcw                                                
      c_d = -Q_Ts/tcw                                                  
    endif

!   if ( c_0 < 0.0 ) then
!     write(*,'(a,2F12.6,10F10.6)') ' c_0, c_d = ',c_0,c_d,B_c,zc_ts,Hb,bc1,bc2,z_c,cc1,cc2,cc3,z_c**2
!   endif

!   c_0 = z_c*Q_Ts/tcw                                                
!   c_d = -Q_Ts/tcw                                                  

  else
    c_0 = 0.0
    c_d = 0.0
  endif                      !  if ( z_c > 0.0 ) then

 end subroutine cool_skin
!
!======================
!
 real function Int_epn(z1,z2,zmx,ztr,n)
!
!  abstract: calculate a definitive integral of an exponetial curve (power of 2)
!
   real(kind_phys) :: z1,z2,zmx,ztr,zi
   real(kind_phys) :: fa,fb,fi,Int
   integer :: m,i,n

   m = nint((z2-z1)/delz)
   fa = exp(-exp_const*((z1-zmx)/(ztr-zmx))**n)
   fb = exp(-exp_const*((z2-zmx)/(ztr-zmx))**n)
   Int = 0.0
   do i = 1, m-1
     zi = z1 + delz*float(i)
     fi = exp(-exp_const*((zi-zmx)/(ztr-zmx))**n)
     Int = Int + fi
   enddo
     Int_epn = delz*((fa+fb)/2.0 + Int)
 end function Int_epn

 subroutine dtl_reset_cv(xt,xs,xu,xv,xz)
 real(kind=kind_phys), intent(inout) :: xt,xs,xu,xv,xz
    xt   =  0.0
    xs   =  0.0
    xu   =  0.0
    xv   =  0.0
    xz   = z_w_max
 end subroutine dtl_reset_cv

 subroutine dtl_reset(xt,xs,xu,xv,xz,xzts,xtts)
 real(kind=kind_phys), intent(inout) :: xt,xs,xu,xv,xz,xzts,xtts
    xt   =  0.0
    xs   =  0.0
    xu   =  0.0
    xv   =  0.0
    xz   = z_w_max
    xtts = 0.0
    xzts = 0.0
 end subroutine dtl_reset


end module nst_module

