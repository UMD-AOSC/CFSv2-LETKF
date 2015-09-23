 MODULE coupler_module
!
!---- error return code ---------
!   0: success
!   1: coupler init
!   3: Atm dim
!   4: SST
!   5: FLUX
!   6: coupler finilize
!----------------------
!
    implicit none
!
    include 'mpif.h'
!
!define domain type
    type domain2D
       integer npes
       integer,dimension(2) :: layout
       integer xdimg,ydimg
       integer,dimension(:),allocatable :: is,ie,js,je
    end type domain2D
    type(domain2D) :: Atm_domain
!
!define MPI communicator
    type :: communicator
      character(len=32) :: name
      integer,dimension(:),allocatable  :: list  !here we only need global rank list corresponding to 0--npes-1
      integer           :: count
      integer           :: id, group         ! MPI communicator and group id for this PE set.
      integer           :: master_rank       ! mast rank in MPI_COMM_WORLD
      integer           :: coupler_rank      ! coupler rank in MPI_COMM_WORLD
      integer           :: FlexLev           ! 
    end type communicator
    type(communicator) :: Atmosphere_comm
    type(communicator) :: Ocean_comm
    type(communicator) :: Coupler_comm
!
!parameters
   integer,parameter :: kind_REAL=8,kind_INTEGER=4,kind_alt_REAL=12-kind_REAL
   real,parameter :: dinr=180./3.1415926535897932
   real,parameter :: TZERO=273.16,TFREEZ=273.16-1.79
   real,parameter :: JCAL=4.1855,HF_CONV_A2O=1.E-4/JCAL
   real,parameter :: rhow=1.E3,cw=4200
!
!MPI vars
   integer :: Coupler_id=0, Atmos_id=1, Ocean_id=2
   integer :: ierr,rc,nprocs
   integer :: Coupler_rank,Atmos_rank,Ocean_rank
   integer :: MPI_kind_REAL,MPI_kind_alt_REAL
!
   integer :: VerbLev=2,nprint=6, mask_tolerance_level=2
   integer :: nprper=120,npr1st=0
!

   private


!*** need to change dt_c, not sent from ocean
   real(kind_REAL) :: dt_c=0., dt_o=0.
   integer :: ndto2dtc
!logical 
   logical :: WRITE_AM_SST=.true.
   logical :: FLUXDEBUG=.false.
   logical :: cice_cover=.false.
   logical :: lsout_cc_ocn=.false.
!
   logical :: restart, INVERTind=.true.
!
!vars
   character*1 restart_old
   character*8 date
   character*10 time
   real :: OM_ul_thickn,efoldingtime
!**** no equivalance between cstepmax and ostepmax, need check if it's necessary ***
   integer cstep,cstepbeg,cstepmax,cstepres
   integer ostepmax 
   character*20 cs
   character*20 sb, s, sr
   character*80 s1,s2
!
! Data on atmos grid
   integer anx,any,anxny,nptsmom4lnd
   real(kind=kind_REAL),dimension(:),allocatable:: ALON,ALAT
   real(kind=kind_REAL),dimension(:,:),allocatable:: t_bot, &
     u_bot,v_bot,q_bot,p_bot,p_surf,xmu,flux_dlw, flux_dsw, &
     ffmm,ffhh, z_bot,t_sfc,t_fix,snw,lprec,fi_sfc,hi_sfc,fi_fix,hi_fix
   real(kind=kind_REAL),dimension(:,:),allocatable:: &
     flux_lw,flux_sw, dtsfc,dusfc,dvsfc,dqsfc
   real(kind=kind_REAL),dimension(:,:),allocatable:: t_sst, fice, &
       hice,hsno,cice
   real(4),dimension(:,:),allocatable:: var4
   integer,dimension(:),allocatable:: nxmomlnd,nymomlnd
!
   NAMELIST /CPL_SETTINGS/ &
       restart,            & !<- no default value, must be specified
       cstepbeg,           & !<- no default value, must be specified
       cstepmax,           & !<- no default value, must be specified
       cstepres,           & !<- no default value, must be specified
       ostepmax,           & !<- see equivalence(...) above; either cstepmax or
                             !   this must be specified
       dt_c, cice_cover, &
       nprint, &
       VerbLev,nprper,npr1st, &
       WRITE_AM_SST 
!
!public method
   public coupler_init,coupler_run,coupler_finalize,GLOB_ABORT,cpl_announce
     
   contains
!
!*********************************************************************
!
   subroutine coupler_init(iret)      
!----------------------------------------------------------------------
! set up MPI for coupler
!----------------------------------------------------------------------
   implicit none
!
   integer,optional,intent(out)  :: iret
   integer ios
!
   iret=1
!
!define NPI REAL:
   if (kind_REAL.eq.8) then
      MPI_kind_REAL=MPI_REAL8
      MPI_kind_alt_REAL=MPI_REAL4
   else if (kind_REAL.eq.4) then
      MPI_kind_REAL=MPI_REAL4
      MPI_kind_alt_REAL=MPI_REAL8
   else
      write(s,'(i0)') kind_REAL
      call GLOB_ABORT(1,'CPL_INIT: illegal value of kind_REAL='//s,1)
   end if
!
!read namelist
   open(14,file='cpl_nml',status='old',iostat=ios)
   call GLOB_ABORT(ios,'Error opening file cpl_nml. C TERMINATED',1)
   read(14,NML=CPL_SETTINGS,iostat=ios)
   write(nprint,NML=CPL_SETTINGS)
! ** read nml error
!   call GLOB_ABORT(ios,'Error reading file cpl_nml. C TERMINATED',1)

!
   if (nprint.eq.7) then
     open(nprint,file='C_printout',form='formatted',status='unknown')
   else if (nprint.ne.6) then
     call date_and_time(date,time)
     open(nprint,file='C_printout.'// &
     date(5:6)//'.'//date(7:8)//'.'//date(3:4)//'_'// &
     time(1:2)//'.'//time(3:4)//'.'//time(5:6), &
     form='formatted',status='unknown')
   end if
!
    write(sb,'(i10)') cstepbeg
    write(s,'(i10)')  cstepmax
    write(sr,'(i10)') cstepres
    call CPL_ANNOUNCE('cstepbeg='//sb//'cstepmax='//s//' cstepres='//sr,0)
!
!RESTART --- restart read from file restart ='//restart,0)
! 02/07/05: value 'T' of restart had a special meaning in pre-MOM
! environment. In that environment, cstep init. value was 0 and
! if restart='T' the 0 time step did not send SST to or receive
! surface fluxes from AM (data in fluxes_for_OM file were used
! in lieu of the latter)

    restart_old='F'
    if (restart_old.ne.'F') CALL GLOB_ABORT(1, &
      'wrong restart_old value: restart_old='//restart_old,1)
!
!---  Initialize MPI ---
    call CPL_ANNOUNCE('MPI initialization to begin',1)
    call MPI_INIT(ierr)
    CALL GLOB_ABORT(ierr,'C: ABORTED upon CALL MPI_INIT',1)
    call CPL_ANNOUNCE('Coupler: back from MPI_INIT',1)
!
    call CPL_INIT
    call CPL_ANNOUNCE('back from CPL_INIT',1)
!
    call CPL_ANNOUNCE('before call CPL_INTRO',1)
    call CPL_INTRO(Atmos_id)
    call CPL_ANNOUNCE('after call CPL_INTRO',1)
    write(s,'(i2)') Atmos_id
    call CPL_ANNOUNCE('back from CPL_INTRO(Atmos_id), Atmos_id='//s,0)
!
    call CPL_INTRO(Ocean_id)
    write(s,'(i2)') Ocean_id
    call CPL_ANNOUNCE('back from CPL_INTRO(Ocean_id), Ocean_id='//s,0)
!
    call CPL_ANNOUNCE('MPI initialization completed',1)
!
    iret=0
   end subroutine coupler_init
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   subroutine cpl_init(iret) 
!----
    implicit none
!
    integer,optional,intent(in) :: iret
    integer color,key,natmpes,nompes
    integer,dimension(:),allocatable :: ids
    character(256) s
    integer :: i,Comm_Coupler,ibuffer_size=4,tag, world_group,ocean_group,nocn,natm
    integer,dimension(:),allocatable :: ibuffer
!
    Coupler_comm%name="Coupler"
    Coupler_comm%id=0
    Coupler_comm%count=1
    allocate(Coupler_comm%list(Coupler_comm%count))
!
    color=Coupler_comm%id
    key=1
    call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,key,COMM_Coupler,ierr)
    call GLOB_ABORT(ierr,'CPL_INIT: error in MPI_COMM_SPLIT',1)
    call CPL_ANNOUNCE('Coupler:back from MPI_Comm_split',0)
    call CPL_ANNOUNCE('Coupler:  IN cpl_init,before mpi_comm_split',0)
!
    call MPI_COMM_RANK(MPI_COMM_WORLD,Coupler_comm%master_rank,ierr)
    call GLOB_ABORT(ierr,'CPL_INIT: error in MPI_COMM_RANK',1)
    write(s,'(i2)') Coupler_comm%master_rank
    call CPL_ANNOUNCE('Coupler:back from MPI_Comm_rank,coupler_rank='//s,0)
!
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
    call GLOB_ABORT(ierr,'CPL_INIT: error in MPI_COMM_SIZE',1)
    write(s,'(i4)') nprocs
    call CPL_ANNOUNCE('Coupler:back from MPI_Comm_size,COMM_WORLD size='//s,0)
!
    Coupler_comm%list(1)=Coupler_comm%master_rank
    Coupler_comm%coupler_rank=Coupler_comm%master_rank
!
    allocate(ids(nprocs*2),ibuffer(ibuffer_size))
    ibuffer(2)=Coupler_rank
    ibuffer(1)=Coupler_id
    tag=1
    do i=0,nprocs-1
      if (i.ne.Coupler_comm%master_rank) then
        call MPI_SEND(ibuffer,2,MPI_INTEGER,i,tag, &
        MPI_COMM_WORLD,ierr)
    write(s,'(i3)') i
    call CPL_ANNOUNCE('Coupler: back from MPI_Comm_send,nproc='//s,0)
        if (ierr.ne.0) then
           write(s,'(i5)') i
           write(nprint,*) 'C: '//"CPL_INIT: stopped, error in MPI_SEND to proc. "//trim(s)
           CALL MPI_ABORT(MPI_COMM_WORLD,2,ierr)
        endif
       endif
   enddo
!
   ibuffer(1)=Coupler_id
   ibuffer(2)=Coupler_comm%master_rank  !global rank
!   ibuffer(3)=0                         !local rank
   call MPI_GATHER(ibuffer,2,MPI_INTEGER,ids,2,MPI_INTEGER,  &
       Coupler_rank,MPI_COMM_WORLD,ierr)
    call CPL_ANNOUNCE('Coupler: back from MPI_gather',0)
!
   write(nprint,*) 'from cpl pe ',Coupler_comm%master_rank,'ids=',ids
!
   Atmosphere_comm%name="Atmosphere"
   Atmosphere_comm%id=Atmos_id
   Atmosphere_comm%coupler_rank=Coupler_comm%coupler_rank
   Ocean_comm%name="Ocean"
   Ocean_comm%id=Ocean_id
   Ocean_comm%coupler_rank=Coupler_comm%coupler_rank
   Atmosphere_comm%count=0
   Ocean_comm%count=0
   do i=1,nprocs
     if(ids(2*i-1).eq.Atmos_id) Atmosphere_comm%count=Atmosphere_comm%count+1
     if(ids(2*i-1).eq.Ocean_id) Ocean_comm%count=Ocean_comm%count+1
   enddo
!
   if(Atmosphere_comm%count+Ocean_comm%count+Coupler_comm%count .ne. nprocs) &
     call CPL_ANNOUNCE('WARNING: there are some pes inactive',1)
!
   allocate(Atmosphere_comm%list(Atmosphere_comm%count),  &
              Ocean_comm%list(Ocean_comm%count)  )
   natm=0
   nocn=0
   do i=1,nprocs
     if(ids(2*i-1).eq.Atmos_id) then
         natm=natm+1
         Atmosphere_comm%list(natm)=ids(2*i)
     endif
     if(ids(2*i-1).eq.Ocean_id) then
       nocn=nocn+1
       Ocean_comm%list(nocn)=ids(2*i)
     endif
   enddo
!
!let ocean pe knows it list
   do i=1,Ocean_comm%count
        call MPI_SEND(Ocean_comm%list,Ocean_comm%count,MPI_INTEGER, &
          Ocean_comm%list(i),1, MPI_COMM_WORLD,ierr)
    write(s,'(i2)') i
    call CPL_ANNOUNCE('Coupler: back from MPI_Comm_send to ocean '//s,0)
   enddo


   return
  end subroutine cpl_init
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   subroutine cpl_intro(component_id)
!----
   implicit none
!
   integer,intent(inout) ::  component_id
!
   integer tag,i,istatus(MPI_STATUS_SIZE),ierr
   integer(kind_INTEGER),dimension(2) :: ibuffer
   logical component_inactive
   character(32) s
!
!
   if (component_id.le.0)  &
      call GLOB_ABORT(1,'C: CPL_INTRO: component_id nonpositive, terminated',1)
!
   component_inactive=.true. 
   if (Atmosphere_comm%id.eq.component_id .or. Ocean_comm%id.eq.component_id) &
      component_inactive=.false.
!

   if (component_inactive) then
      component_id=-component_id
      print*,'C: CPL_INTRO: component with id=',-component_id, &
       ' INACTIVE; id made =',component_id
       RETURN
   end if
  
    tag=component_id-1          ! to keep tag=1 for ocean
    ibuffer=0
    print *,'Coupler: in cmp_intro, before recv,tag=',tag
    call MPI_RECV(ibuffer,2,MPI_INTEGER,MPI_ANY_SOURCE,tag, &
      MPI_COMM_WORLD,istatus,ierr)
    call GLOB_ABORT(ierr,'CPL_INTRO: error in MPI_RECV',1)
    call CPL_ANNOUNCE('Coupler:back from MPI_recv from ATM/ocean to get local_root',0)
!
    if (ibuffer(1).ne.component_id) then
       write(nprint,'("C: CPL_INTRO: stopped, received component id value",i9,'// &
        ' "not equal to component_id argument",i9)') ibuffer(1),component_id 
       call GLOB_ABORT(1,'CPL_INTRO: error',2)
    end if
!
    if( component_id.eq.Atmosphere_comm%id ) then
       Atmosphere_comm%master_rank=ibuffer(2) 
    else if( component_id.eq.Ocean_comm%id ) then 
       Ocean_comm%master_rank=ibuffer(2) 
    else
       write(s,'(i5)') component_id
       call GLOB_ABORT(1,'CPL_INTRO: component id wrong'//trim(s),2)
    endif
!
     return
   end subroutine cpl_intro
!----------------------------------------------------------------------
   subroutine cpl_announce(s,DbgLev) 
!---
     implicit none
!
     character*(*),intent(in) :: s
     integer,intent(in) :: DbgLev
!
     call date_and_time(date,time)
     if (DbgLev.le.VerbLev)  write(nprint,*) 'C: '//time//' '//s

      return
   end subroutine cpl_announce
!----------------------------------------------------------------------
   subroutine GLOB_ABORT(ie,s,rc)
!---
     implicit none
!
     integer rc,ie,ierr
     character*(*) s
     if (ie.ne.0) then
        print*,'GLOB_ABORT: '//s//' ie,rc:',ie,rc
        if (rc.eq.0) RETURN
        CALL MPI_ABORT(MPI_COMM_WORLD,rc,ierr)
     end if
     RETURN
   end subroutine GLOB_ABORT
!----------------------------------------------------------------------
!
!*********************************************************************
!
   subroutine coupler_run(iret) 
!----------------------------------------------------------------------
! This subroutine contains all the transmits coupler has done including
!   receiving SST from OCEAN and sending to ATMOS, and sending flux to
!   OCEAN and receiving FLUX from ATMOS 
!----------------------------------------------------------------------
   implicit none
!
     integer,optional,intent(out) :: iret
     integer ios,i
     integer,dimension(:),allocatable :: ids
     integer npes,ierr,ibuf(2)
     integer,dimension(MPI_STATUS_SIZE) :: istatus
     real(kind_REAL) buf(2)
     real buff(2)
!
     iret=2
!
!-----------------------------------------------------------------------
!--- STEP 1: receive OM time step; send coupling OM and AM
!
! --1.1 get dt_c, dto from ocean
!
     call MPI_RECV(buf,2,MPI_kind_REAL,Ocean_comm%master_rank,Ocean_comm%id-1, &
        MPI_COMM_WORLD,istatus, ierr)
      call GLOB_ABORT(ierr, "error in transferring ATM anx/any",1)
      dt_c=buf(1)
      dt_o=buf(2)
!
      write(s1,'(e20.12,e20.12)') buf(1),buf(2)
      call CPL_ANNOUNCE('From Ocean, dt_cc, dto_cc: '//trim(s1),1)
!
! --1.2 get domain division from MOM4 atmos

     if (Ocean_comm%count .lt. 1)  &
       call GLOB_ABORT(2001,'No Ocean pes,domain size from Ocean_comm master',1)
     allocate( ids(Ocean_comm%count*6) )

     npes=Ocean_comm%count
     call recv_data_domain(npes, ids,ierr )
     if( ierr .ne.0 ) call GLOB_ABORT(2001,'Coupler cannot recv local'// &
       'domain size from Ocean_comm master',1)
!
!-----------------------------------------------------------------------
!--- STEP 2: recevie AM grid dimention and allocate arrays 
!---      ** currently, coupler only need to know ATM grid **
!

      call MPI_RECV(ibuf,2,MPI_INTEGER,Atmosphere_comm%master_rank,Atmosphere_comm%id-1, &
        MPI_COMM_WORLD,istatus, ierr)
      call GLOB_ABORT(ierr, "error in transferring ATM anx/any",1)
      anx=ibuf(1)
      any=ibuf(2)
      anxny=anx*any
!
      write(s,'(i3,2i6)') Atmos_id,anx,any
      call CPL_ANNOUNCE('Atmos_id, anx, any: '//s,1)
!
!      dt_c=600.
!      dt_o=3600.
      buf(1)=dt_c
      buf(2)=dt_o
!
      call MPI_SEND(buf,2,MPI_kind_REAL,Atmosphere_comm%master_rank,Atmosphere_comm%id-1, &
        MPI_COMM_WORLD,ierr)
      call GLOB_ABORT(ierr, "error in send Atmos the ice/ocean time step",1)
      write(s1,'(e20.12,e20.12)') buf(1),buf(2)
      call CPL_ANNOUNCE('Coupler: send dt_c/dt_o to ATM'//trim(s1),0)
!
! Allocate arrays:
!AM
      allocate(ALON(anx),ALAT(any))
      allocate(t_bot(anx,any), &
        u_bot(anx,any),v_bot(anx,any),q_bot(anx,any), &
        z_bot(anx,any), t_sfc(anx,any), t_fix(anx,any), &
        p_surf(anx,any),p_bot(anx,any),ffmm(anx,any),ffhh(anx,any), &
        snw(anx,any),lprec(anx,any),                                &
        fi_sfc(anx,any), hi_sfc(anx,any), fi_fix(anx,any), hi_fix(anx,any) )
      allocate(flux_lw(anx,any),flux_sw(anx,any), &
        xmu(anx,any),flux_dlw(anx,any),flux_dsw(anx,any))
      allocate(dusfc(anx,any),dvsfc(anx,any), &
        dtsfc(anx,any),dqsfc(anx,any) )
      allocate(t_sst(anx,any), fice(anx,any),hice(anx,any), &
               cice(anx,any),hsno(anx,any) )
      allocate(var4(anx,any))
!ocean 
!     dummy
!
!-- Implicit mpp_define domain
    call set_domain(Atm_domain,anx,any,npes,ids)
!-----------------------------------------------------------------------
!--- Step 3: receive ATM grid and ICE grid in order to do interpolation 
!---      **  currently, indeed no grid info need to be know since interpolation
!---          is done inside mom, but to keep GFS not changed, get ATM grid
!
!-->test: now use pseudo atmos 
!      IF (Atmos_id.gt.0) THEN
!        call CPL_ANNOUNCE('to receive AMG arrays (2 MPI calls)',3)
!       call CPL_RECV(ALON,anx,Atmos_id)
!       call CPL_RECV(ALAT,any,Atmos_id)
!       call CPL_ANNOUNCE('AMG arrays received (2 MPI calls)',1)
!      ELSE          !
                    ! this is for testing only, to allow interpolation
                    ! to work somehow in the case of missing AM
!        do j=1,any    !
!          ALAT(j)=90.+90./any-j*180./any
!        end do        !
!        do i=1,anx    !
!          ALON(i)=(i-1)*360./anx
!        end do        !
!      call CPL_ANNOUNCE( &
!     'Atmos_id<0; test AM gridpoint coordinates assigned',0)
!      END IF
!Convert rads to degs if required
!      if (ALON(anx).lt.7.) then
!        call CPL_ANNOUNCE('to convert AM grid from rads to degs',2)
!        ALON=ALON*dinr
!        ALAT=ALAT*dinr
!      end if
!
!      call PRI2D(anx,1,ALON,'AMG longitudes:')
!      call PRI2D(1,any,ALAT,'AMG latitudes:')
!<-- test
!
!-----------------------------------------------------------------------
!--- Step 4:receive sea/land mask from ATM 
!---     ** currently, no surface mask sent from OM 
!
!-->test: now use pseudo atmos 
!      IF (Atmos_id.gt.0) THEN
!        call CPL_ANNOUNCE('to receive AMG sea/land mask from AM',3)
!        call CPL_RECV(SLM_a,anxny,Atmos_id)
!        call CPL_ANNOUNCE('AMG sea/land mask received from AM',1)
!      ELSE           !
!                     ! this is for testing only, to allow interpolation
!        SLM_a=a_sea   ! to work somehow in the case of missing OM
!        CALL RANDOM_NUMBER(SLM_a)
!        SLM_a=nint(SLM_a)
!        i1=0           !
!        do j=1,any     !
!        do i=1,anx     !
!          i1=i1+SLM_a(i,j)+0.01
!        end do         !
!        end do         !
!        if (i1 .eq. a_sea*anxny) then
!          call CPL_ANNOUNCE('Atmos_id<0; AM SLM assigned sea values:',0)
!        else if (i1 .eq. (1-a_sea)*anxny) then
!          call CPL_ANNOUNCE('Atmos_id<0; AM SLM assigned land values:',0)
!        else           !
!          call CPL_ANNOUNCE('Atmos_id<0; AM SLM assigned mixed values:',0)
!        end if         !
!      END IF
!<-- test
!
!-----------------------------------------------------------------------
!--- Step 5: time loop
!---
!
!--- 5.0 open file in order to write out flux for OM and sst for ATM
!
      if (restart) then
        OPEN(10,file='fluxes_for_OM', &
            form='unformatted' ,status='old',iostat=ios)
        if (ios.ne.0) & 
          CALL GLOB_ABORT(1,'C: no file fluxes_for_OM : terminated',1)
        call read_restart(10)
        close (10)
      else
        call set_data_init()
      end if

      call lndseaindex
!
!--- 5.1 time steps
!cstepmax
!
     do cstep=cstepbeg,cstepmax  ! <- here and elsewhere ostep replaced by
                          ! cstep 8/23/05
       write(cs,'(i0)') cstep
!
!-- receive SST from OM, no interpolation,send it to ATM
!
       IF (cstep.gt.0 .or. restart_old.eq.'F') THEN

!send sst to atm
!
         call SEND_SST

!recev flux from Atmos
         ndto2dtc=nint(dt_o/dt_c+0.001)
         lsout_cc_ocn=(mod(cstep,ndto2dtc).eq.0)
         call RECV_FLUXES
!
!send fluxes to ocean

         call SEND_FLUXES
!
         call RECV_FROM_OCN
!
      ENDIF
!
      if (mod(cstep,cstepres) == 0 .and. cstep > 0 ) then

        OPEN(10,file='fluxes_for_OM_'//cs,form='unformatted', iostat=ios)

        call write_restart(10,lsout_cc_ocn)

        call CPL_ANNOUNCE( &
         'vars for AM/OM written to fluxes_for_OM file for cstep='//cs,0)
        close(10)
      endif
     enddo  !end cstep loop5.1
!--
     iret=0
     return
   end subroutine coupler_run

   subroutine write_restart(iunit,lsout_cc_ocn)
!
       implicit none
       integer iunit
       logical lsout_cc_ocn
       rewind iunit
!                                   write vars for AM
       write (iunit) t_sst
       write (iunit) fice
       write (iunit) hice
       write (iunit) hsno
!                                   write vars for OM
       write (iunit) t_sfc
       write (iunit) t_bot
       write (iunit) u_bot
       write (iunit) v_bot
       write (iunit) q_bot
       write (iunit) p_bot
       write (iunit) p_surf
       write (iunit) z_bot
       write (iunit) xmu
       write (iunit) flux_dsw
       write (iunit) flux_dlw
       write (iunit) ffmm
       write (iunit) ffhh
       write (iunit) snw
       write (iunit) lprec
       if (lsout_cc_ocn) then
         write (iunit) dusfc
         write (iunit) dvsfc
         write (iunit) dtsfc
         write (iunit) dqsfc
         write (iunit) flux_lw
         write (iunit) flux_sw
      endif
!-- 
      return
   end subroutine write_restart
!
!*********************************************************************
!
   subroutine coupler_finalize(iret) 
!----------------------------------------------------------------------
!  coupler finalize
!----------------------------------------------------------------------
!
   implicit none
   integer,optional,intent(out) :: iret
   iret=3
!
   deallocate(t_sfc)
   deallocate(t_fix)
   deallocate(t_bot)
   deallocate(u_bot)
   deallocate(v_bot)
   deallocate(p_bot)
   deallocate(q_bot)
   deallocate(z_bot)
   deallocate(p_surf)
   deallocate(ffmm)
   deallocate(ffhh)
   deallocate(snw)
   deallocate(lprec)
   deallocate(xmu)
   deallocate(flux_dsw)
   deallocate(flux_dlw)
   deallocate(flux_sw)
   deallocate(flux_lw)
   deallocate(dusfc)
   deallocate(dvsfc)
   deallocate(dtsfc)
   deallocate(dqsfc)
   deallocate(t_sst)
   deallocate(fice)
   deallocate(hice)
   deallocate(hsno)
   deallocate(cice)
   deallocate(var4)
   deallocate(fi_sfc)
   deallocate(hi_sfc)
   deallocate(fi_fix)
   deallocate(hi_fix)
!
! Finalize MPI
!      call MPI_Barrier(MPI_COMM_WORLD, iret)
!      CALL GLOB_ABORT(iret,'C: ABORTED upon CALL MPI_Barrier',1)
!      write(nprint,*)'C: to call MPI_FINALIZE, after call mpi_barrier,iret=',iret
      call MPI_FINALIZE(ierr)
      CALL GLOB_ABORT(ierr,'C: ABORTED upon CALL MPI_FINALIZE',1)
      write(nprint,*)'C: back from MPI_FINALIZE,ierr=',ierr
!
!        write(nprint,*)'C: done'
    iret=0
    return
   end subroutine coupler_finalize
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   subroutine RECV_FROM_OCN
!----------------------------------------------------------------------
!  receive SST from ocean model
!----------------------------------------------------------------------
!
    implicit none

    real(kind=kind_REAL),dimension(:,:),allocatable :: tclim_sst
    real(kind=kind_REAL) :: tice
    real(kind=kind_REAL) :: eps_fi=1.e-6
    integer i,j,k
    logical :: cice_avail=.true.

    allocate(tclim_sst(anx,any))
!
! t_sst
!
    call RECV_sglfield(t_sst,Atm_domain,Ocean_comm%id)
    write(s,'(i10)') cstep
    call CPL_ANNOUNCE('to receive OMG SST from OM(knowns anxny); anxny, cstep: '//s,3)
!
! c_sst
!
    call RECV_sglfield(tclim_sst,Atm_domain,Ocean_comm%id)
    write(s,'(i10)') cstep
    call CPL_ANNOUNCE('to receive OMG clim SST from OM(knowns anxny); anxny, cstep: '//s,3)
!
! fice
!
    call RECV_sglfield(fice,Atm_domain,Ocean_comm%id)
    write(s,'(i10)') cstep
    call CPL_ANNOUNCE('to receive OMG FICE from OM(knowns anxny); anxny, cstep: '//s,3)
!
! hice
!
    call RECV_sglfield(hice,Atm_domain,Ocean_comm%id)
    write(s,'(i10)') cstep
    call CPL_ANNOUNCE('to receive OMG HICE from OM(knowns anxny); anxny, cstep: '//s,3)
!
! hsno
!
    call RECV_sglfield(hsno,Atm_domain,Ocean_comm%id)
    write(s,'(i10)') cstep
    call CPL_ANNOUNCE('to receive OMG Hsno from OM(knowns anxny); anxny, cstep: '//s,3)
!
! cice
!
!write out partial data
!
    if (cice_avail) then
    call RECV_sglfield(cice,Atm_domain,Ocean_comm%id)
    write(s,'(i10)') cstep
    call CPL_ANNOUNCE('to receive OMG Cice from OM(knowns anxny); anxny, cstep: '//s,3)
!
    endif
    if(cice_cover) then
       do j=1,any
       do i=1,anx
          fice(i,j)=min(1.0,(fice(i,j)+2.0*cice(i,j))/3.0)
          hice(i,j)=(hice(i,j)+2.0*fice(i,j)*fice(i,j))/3.0
       enddo
       enddo
    endif
!
    call INVIND(size(t_fix,1),size(t_fix,2), t_fix)
    call INVIND(size(fi_fix,1),size(fi_fix,2), fi_fix)
    call INVIND(size(hi_fix,1),size(hi_fix,2), hi_fix)
    do k=1,nptsmom4lnd
       t_sst(nxmomlnd(k),nymomlnd(k)) = t_fix(nxmomlnd(k),nymomlnd(k))
       fice(nxmomlnd(k),nymomlnd(k))  = fi_fix(nxmomlnd(k),nymomlnd(k))
       hice(nxmomlnd(k),nymomlnd(k))  = hi_fix(nxmomlnd(k),nymomlnd(k))
    enddo

    deallocate(tclim_sst)
!
      return
   end subroutine RECV_FROM_OCN
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   SUBROUTINE SEND_SST
!----------------------------------------------------------------------
!  send SST to ATM
!----------------------------------------------------------------------
!
   implicit none
   integer :: ierr, t_sst_max,t_sst_min,i,j
!
   t_sst_max=t_sst(1,1)
   t_sst_min=t_sst(1,1)
   do j=1,size(t_sst,2)
    do i=1,size(t_sst,1)
     if(t_sst(i,j).gt. t_sst_max) t_sst_max=t_sst(i,j)
     if(t_sst(i,j).lt. t_sst_min) t_sst_min=t_sst(i,j)
    enddo
   enddo

!      
! change var from N->S(AM) to S->N(MOM)
    call INVIND(size(t_sst,1),size(t_sst,2), t_sst)

   call MPI_SEND(t_sst,size(t_sst,1)*size(t_sst,2),MPI_kind_REAL, &
       Atmosphere_comm%master_rank, &
       Atmosphere_comm%id-1, MPI_COMM_WORLD,ierr)  

   call CPL_ANNOUNCE('back from MPI_SEND(SST_a,anxny,Atmos_id) ',2)
!
!SEND OUT FICE HICE

    call INVIND(size(fice,1),size(fice,2), fice)

   call MPI_SEND(fice,size(fice,1)*size(fice,2),MPI_kind_REAL, &
       Atmosphere_comm%master_rank, &
       Atmosphere_comm%id-1, MPI_COMM_WORLD,ierr)
 
   call CPL_ANNOUNCE('back from MPI_SEND(FICE_a,anxny,Atmos_id) ',2)

    call INVIND(size(hice,1),size(hice,2), hice)

   call MPI_SEND(hice,size(hice,1)*size(hice,2),MPI_kind_REAL, &
       Atmosphere_comm%master_rank, &
       Atmosphere_comm%id-1, MPI_COMM_WORLD,ierr)

   call CPL_ANNOUNCE('back from MPI_SEND(HICE_a,anxny,Atmos_id) ',2)

    call INVIND(size(hsno,1),size(hsno,2), hsno)

   call MPI_SEND(hsno,size(hsno,1)*size(hsno,2),MPI_kind_REAL, &
       Atmosphere_comm%master_rank, &
       Atmosphere_comm%id-1, MPI_COMM_WORLD,ierr)

   call CPL_ANNOUNCE('back from MPI_SEND(HSNO_a,anxny,Atmos_id) ',3)

   RETURN
   END SUBROUTINE SEND_SST
!
!*********************************************************************
!
   SUBROUTINE RECV_FLUXES
!----------------------------------------------------------------------
!  receive fluxes from ATM
!----------------------------------------------------------------------
!
      implicit none

      call CPL_ANNOUNCE('to receive flux from AM for sea-ice model, '//cs,3)

      call CPL_RECV(t_sfc,anx,any,Atmos_id)
      call CPL_ANNOUNCE('recv fluxes from ATM: t_sfc',1)
      t_fix=t_sfc

      call CPL_RECV(t_bot,anx,any,Atmos_id)
      call CPL_ANNOUNCE('recv fluxes from ATM: t_bot',1)

      call CPL_RECV(u_bot,anx,any,Atmos_id)
      call CPL_ANNOUNCE('recv fluxes from ATM: u_bot',1)

      call CPL_RECV(v_bot,anx,any,Atmos_id)
      call CPL_ANNOUNCE('recv fluxes from ATM: v_bot',1)

      call CPL_RECV(q_bot,anx,any,Atmos_id)
      call CPL_ANNOUNCE('recv fluxes from ATM: q_bot',1)

      call CPL_RECV(p_bot,anx,any,Atmos_id)
      call CPL_ANNOUNCE('recv fluxes from ATM: p_bot',1)

      call CPL_RECV(p_surf,anx,any,Atmos_id)
      call CPL_ANNOUNCE('recv fluxes from ATM: p_surf',1)

      call CPL_RECV(z_bot,anx,any,Atmos_id)
      call CPL_ANNOUNCE('recv fluxes from ATM: z_bot',1)

      call CPL_RECV(xmu,anx,any,Atmos_id)
      call CPL_ANNOUNCE('recv fluxes from ATM: xmu',1)

      call CPL_RECV(flux_dlw,anx,any,Atmos_id)
      call CPL_ANNOUNCE('recv fluxes from ATM: dlw',1)

      call CPL_RECV(flux_dsw,anx,any,Atmos_id)
      call CPL_ANNOUNCE('recv fluxes from ATM: dsw',1)

      call CPL_RECV(ffmm,anx,any,Atmos_id)
      call maxmin(anx,any,ffmm,'in Coupler, ffmm')
      call CPL_ANNOUNCE('recv fluxes from ATM: ffmm',1)
      call CPL_RECV(ffhh,anx,any,Atmos_id)
      call maxmin(anx,any,ffhh,'in Coupler, ffhh')
      call CPL_ANNOUNCE('recv fluxes from ATM: ffhh',1)

      call CPL_RECV(snw,anx,any,Atmos_id)
      call CPL_ANNOUNCE('recv snw from ATM: snw',1)
      call CPL_RECV(lprec,anx,any,Atmos_id)
      call CPL_ANNOUNCE('recv lprec from ATM: lprec',1)

      call CPL_RECV(fi_sfc,anx,any,Atmos_id)
      call CPL_ANNOUNCE('recv fluxes from ATM: fi_sfc',1)
      fi_fix=fi_sfc

      call CPL_RECV(hi_sfc,anx,any,Atmos_id)
      call CPL_ANNOUNCE('recv fluxes from ATM: hi_sfc',1)
      hi_fix=hi_sfc

      call CPL_ANNOUNCE('end of recv fluxes from ATM: fast loop',1)

      if ( lsout_cc_ocn) then
        call CPL_RECV(dusfc,anx,any,Atmos_id)
        call CPL_ANNOUNCE('recv fluxes from ATM: dusfc',1)
        call CPL_RECV(dvsfc,anx,any,Atmos_id)
        call CPL_ANNOUNCE('recv fluxes from ATM: dvsfc',1)
        call CPL_RECV(dtsfc,anx,any,Atmos_id)
        call CPL_ANNOUNCE('recv fluxes from ATM: dtsfc',1)
        call CPL_RECV(dqsfc,anx,any,Atmos_id)
        call CPL_ANNOUNCE('recv fluxes from ATM: sqsfc',1)
        call CPL_RECV(flux_lw,anx,any,Atmos_id)
        call CPL_ANNOUNCE('recv fluxes from ATM: swsfc',1)
        call CPL_RECV(flux_sw,anx,any,Atmos_id)
        call CPL_ANNOUNCE('recv fluxes from ATM: swsfc',1)
      endif

 
      return
   END SUBROUTINE RECV_FLUXES
!
!*********************************************************************
!
   SUBROUTINE SEND_FLUXES
!----------------------------------------------------------------------
!  send fluxes to OCEAN. For each surface flux,if required, 
!     fill it with invalid values, to force OM to use its own data;
!     invert indexation if necessary, change sign / convert units 
!     if necessary 
!----------------------------------------------------------------------
!
   implicit none
!
!-- get partial data and send to Ocean pes
    call SEND_sglfield(t_sfc,Atm_domain,Ocean_comm%id)
    call SEND_sglfield(t_bot,Atm_domain,Ocean_comm%id)
    call SEND_sglfield(u_bot,Atm_domain,Ocean_comm%id)
    call SEND_sglfield(v_bot,Atm_domain,Ocean_comm%id)
    call SEND_sglfield(q_bot,Atm_domain,Ocean_comm%id)
    call SEND_sglfield(p_bot,Atm_domain,Ocean_comm%id)
    call SEND_sglfield(p_surf,Atm_domain,Ocean_comm%id)
    call SEND_sglfield(z_bot,Atm_domain,Ocean_comm%id)

    call SEND_sglfield(xmu,Atm_domain,Ocean_comm%id)
    call SEND_sglfield(flux_dlw,Atm_domain,Ocean_comm%id)
    call SEND_sglfield(flux_dsw,Atm_domain,Ocean_comm%id)
!
    call SEND_sglfield(ffmm,Atm_domain,Ocean_comm%id)
    call SEND_sglfield(ffhh,Atm_domain,Ocean_comm%id)
!
    call SEND_sglfield(snw,Atm_domain,Ocean_comm%id)
    call SEND_sglfield(lprec,Atm_domain,Ocean_comm%id)
!
     t_sfc=0.
     t_bot=0.
     u_bot=0.
     v_bot=0.
     q_bot=0.
     p_bot=0.
     p_surf=0.
     z_bot=0.
     xmu=0.
     flux_dlw=0.
     flux_dsw=0.
     ffmm=0.
     ffhh=0.
     snw=0.
     lprec=0.

    if (lsout_cc_ocn) then 
      call SEND_sglfield(dusfc,Atm_domain,Ocean_comm%id)
      call SEND_sglfield(dvsfc,Atm_domain,Ocean_comm%id)
      call SEND_sglfield(dtsfc,Atm_domain,Ocean_comm%id)
      call SEND_sglfield(dqsfc,Atm_domain,Ocean_comm%id)
      call SEND_sglfield(flux_lw,Atm_domain,Ocean_comm%id)
      call SEND_sglfield(flux_sw,Atm_domain,Ocean_comm%id)
       dusfc=0.
       dvsfc=0.
       dtsfc=0.
       dqsfc=0.
       flux_lw=0.
       flux_sw=0.
    endif

   return
   END SUBROUTINE SEND_FLUXES
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

   subroutine recv_data_domain(npes,ids,iret)
!----------------------------------------------------------------
!-- recv data domain
!----------------------------------------------------------------
!
     implicit none
     integer,intent(in) :: npes
     integer,dimension(npes*6),intent(inout) :: ids 
     integer,intent(inout) :: iret
     integer :: tag,istatus(MPI_STATUS_SIZE)
!
     tag=1    !! mpp send
     call MPI_recv(ids,npes*6,MPI_INTEGER,Ocean_comm%master_rank,tag, &
       MPI_COMM_WORLD,istatus,ierr)
     if(ierr .ne. 0) call GLOB_ABORT(2001,'Coupler cannot recv local'// &
       'domain size from Ocean_comm master',1)
    
     iret=0 
     return
    end subroutine  recv_data_domain
!----------------------------------------------------------------
    subroutine set_domain(domain,xdimg,ydimg,npes,ids)
!----------------------------------------------------------------
!-- set domain
!----------------------------------------------------------------
!
     implicit none
!     
     type(domain2D),intent(inout) ::  domain
     integer,intent(in)  :: xdimg,ydimg,npes
     integer,dimension(npes*6),intent(in)  :: ids
!
     integer :: ndiv,i
!
     domain%xdimg=xdimg
     domain%ydimg=ydimg
     domain%npes=npes
     allocate(domain%is(npes),domain%ie(npes),domain%js(npes),domain%je(npes))
!
     write(s,'(i2,i2,i2)') ids(1),ids(2),npes
     if ( ids(1)*ids(2) .ne. npes) call GLOB_ABORT(2001,'in set domain, layout'// &
        ' doesnot match npes'//s,1 )
!
     domain%layout= (/ ids(1), ids(2) /)
     do ndiv=1,npes
       if(ids(6*ndiv-5) .ne. domain%layout(1) .or. ids(6*ndiv-4).ne.domain%layout(2))  &
        call GLOB_ABORT(2001,'in set domain, layout doesnot same in all pes',1)
       domain%is(ndiv)=ids(6*ndiv-3)
       domain%ie(ndiv)=ids(6*ndiv-2)
       domain%js(ndiv)=ids(6*ndiv-1)
       domain%je(ndiv)=ids(6*ndiv)
     enddo
     print *,'in set_domain,set doamin is,ie,js,je'
!
     return
    end subroutine  set_domain
!----------------------------------------------------------------
    subroutine RECV_sglfield(var,domain,component_id)
!----------------------------------------------------------------
!-- giving var on global domain (xdimi,ydim), sebd partial domain data
!     to the Component_id pes
!----------------------------------------------------------------
!
    implicit none
    type(domain2D),intent(in) :: domain
    real(8),dimension(:,:),intent(inout) :: var(domain%xdimg,domain%ydimg)
    integer,intent(in) :: component_id
!
    real(8),dimension(:,:),allocatable :: varpart
    integer,dimension(:),allocatable :: pelist
    integer,dimension(MPI_STATUS_SIZE) :: istatus
    integer ndiv,npes,ndiv_xdim,ndiv_ydim,i,j,tag,xdim,ydim,is,js
!
    print *,'recv_sglfield,domain size=',domain%xdimg,domain%ydimg
!
    if(component_id.eq.Ocean_comm%id) then
      allocate(pelist(Ocean_comm%count))
      pelist=Ocean_comm%list
      npes=Ocean_comm%count
    else if (component_id.eq.Atmosphere_comm%id) then
      allocate(pelist(Atmosphere_comm%count))
      pelist=Atmosphere_comm%list
      npes=Atmosphere_comm%count
    else
      call GLOB_ABORT(2001,"in send_sglfield, invalid component is",1)
    endif
    print *,'in recv_sglfield, npes=',npes,'pelist=',pelist
!
    do ndiv=1,npes

      ndiv_xdim=domain%ie(ndiv)-domain%is(ndiv)+1
      ndiv_ydim=domain%je(ndiv)-domain%js(ndiv)+1
      allocate(varpart(ndiv_xdim,ndiv_ydim))
      print *,'in recv_sglfield,varpart size=',ndiv_xdim,ndiv_ydim
!
!recv OCEAN var part on Atm
      is=domain%is(ndiv)
      js=domain%js(ndiv)
      tag=1   !in mpp_send tag=1

      call MPI_RECV(varpart,ndiv_xdim*ndiv_ydim,MPI_kind_REAL,pelist(ndiv),tag,&
       MPI_COMM_WORLD,istatus,ierr)
!      print *,'print recv_sst,MPI_kind_REAL=',MPI_kind_REAL,'is=',is,'js=',js,&
!         'ie=',domain%ie(ndiv),'je=',domain%je(ndiv),'varpart=',varpart(1:5,:)

      call GLOB_ABORT(ierr,'Coupler read single field to Ocean pes',1)
!
      do j=domain%js(ndiv),domain%je(ndiv)
        do i=domain%is(ndiv),domain%ie(ndiv)
          var(i,j)=varpart(i-is+1,j-js+1)
        enddo
      enddo
!
      deallocate(varpart)
!
    enddo
!
    return
    end subroutine RECV_sglfield
!----------------------------------------------------------------
    subroutine SEND_sglfield(var,domain,component_id)
!----------------------------------------------------------------
!-- giving var on global domain (xdimi,ydim), sebd partial domain data 
!     to the Component_id pes
!----------------------------------------------------------------
!
    implicit none
    type(domain2D),intent(in) :: domain
    real(8),dimension(:,:),intent(inout) :: var(domain%xdimg,domain%ydimg)
    integer,intent(in) :: component_id
!
    real(8),dimension(:,:),allocatable :: varpart
    integer,dimension(:),allocatable :: pelist
    integer ndiv,npes,ndiv_xdim,ndiv_ydim,i,j,tag,xdim,ydim,is,js
!
    if(component_id.eq.Ocean_comm%id) then
      allocate(pelist(Ocean_comm%count))
      pelist=Ocean_comm%list
      npes=Ocean_comm%count
    else if (component_id.eq.Atmosphere_comm%id) then
      allocate(pelist(Atmosphere_comm%count))
      pelist=Atmosphere_comm%list
      npes=Atmosphere_comm%count
    else
      call GLOB_ABORT(2001,"in send_sglfield, invalid component is",1)
    endif
!      
! change var from N->S(AM) to S->N(MOM)
    call INVIND(domain%xdimg, domain%ydimg, var)
!
    do ndiv=1,npes
      ndiv_xdim=domain%ie(ndiv)-domain%is(ndiv)+1
      ndiv_ydim=domain%je(ndiv)-domain%js(ndiv)+1
      allocate(varpart(ndiv_xdim,ndiv_ydim))
!
      do j=domain%js(ndiv),domain%je(ndiv)
        do i=domain%is(ndiv),domain%ie(ndiv)
          varpart(i-domain%is(ndiv)+1,j-domain%js(ndiv)+1)=var(i,j)
        enddo
     enddo
!send to global OCEAN pelist
     is=domain%is(ndiv)
     js=domain%js(ndiv)
     tag=1   !in mpp_send tag=1
     call MPI_SEND(varpart,ndiv_xdim*ndiv_ydim,MPI_kind_REAL,pelist(ndiv),tag,&
       MPI_COMM_WORLD,ierr)
     call GLOB_ABORT(ierr,'Coupler send single field to Ocean pes',1)
!
     deallocate(varpart)
!
   enddo  !end doloop ndiv
!
   return
   end subroutine  SEND_sglfield
!----------------------------------------------------------------
    subroutine cpl_recv(var,xdim,ydim,component_id)
!----------------------------------------------------------------
!-- giving var on global domain (xdimi,ydim), sebd partial domain data
!     to the Component_id pes
!----------------------------------------------------------------
!
    implicit none
    integer,intent(in) :: xdim,ydim,component_id
    real(kind_REAL),dimension(xdim,ydim),intent(inout) :: var
    integer :: varlen
!
    integer :: sendpe,tag,ierr, istatus(MPI_STATUS_SIZE)
!
    if(component_id.eq.Ocean_comm%id) then
      sendpe=Ocean_comm%master_rank
    else if (component_id.eq.Atmosphere_comm%id) then
      sendpe=Atmosphere_comm%master_rank
    else
      call GLOB_ABORT(2001,"in send_sglfield, invalid component is",1)
    endif

    tag=component_id-1
!
    call MPI_RECV(var,xdim*ydim,MPI_kind_REAL,sendpe, tag, MPI_COMM_WORLD, &
       istatus,ierr)
    call GLOB_ABORT(ierr,'CPL_RECV: error in MPI_RECV',1)

   return
   end subroutine  cpl_recv
!----------------------------------------------------------------

!   SUBROUTINE CONVF(F,C,too_low_to_convert)
!----------------------------------------------------------------------
!  change unit
!----------------------------------------------------------------------
!
!   implicit none
!
!   real,intent(in) :: C,too_low_to_convert
!   real,intent(inout) :: F(onx,ony)
!   integer i,j
!!
!      do j=1,ony
!      do i=1,onx
!        if (F(i,j).gt.too_low_to_convert) F(i,j)=C*F(i,j)
!      end do
!      end do

!      RETURN
!      END SUBROUTINE CONVF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   SUBROUTINE PRI2D(IM,JM,F,s)
!----------------------------------------------------------------------
!  PRINT 2D array F(IM,JM) 
!----------------------------------------------------------------------
!
     implicit none
!
      integer,intent(in) ::  IM,JM
      real,dimension(IM,JM),intent(in) :: F
      character*(*) s
      integer n,i,j
      character*9 sn
      character*6 fr
      character*2 sc
!
      write(nprint,*)' '
      n=10

      if (s(1:3).eq.'SST') then
        fr='f12.7)'
        sc=' '
      else
        fr='e12.4)'
        sc='1p'
      end if

      IF (IM.eq.1) THEN
        write(nprint,'("C: '//s//', latitudes from")')
        do j=1,JM,n
          write(sn,'(i0)') min(j+n-1,JM)-j+1
          write(nprint,'("C:",i4," to",i4,'//sc//trim(sn)//fr) &
         j,min(j+n-1,JM),F(1,j:min(j+n-1,JM))
        end do
      ELSE IF (JM.eq.1) THEN
        write(nprint,'("C: '//s//', longitudes from")')
        do i=1,IM,n
          write(sn,'(i0)') min(i+n-1,IM)-i+1
          write(nprint,'("C:",i4," to",i4,'//sc//trim(sn)//fr) &
         i,min(i+n-1,IM),F(i:min(i+n-1,IM),1)
        end do
      ELSE
        do i=1,IM,n
          write(nprint,'("C: '//s//', longitudes from",i4," to",i4)') &
         i,min(i+n-1,IM)
          write(sn,'(i0)') min(i+n-1,IM)-i+1
          do j=1,JM
            write (nprint,'("C:",2i4,'//sc//trim(sn)//fr) &
           i,j,F(i:min(i+n-1,IM),j)
          end do
        end do
      ENDIF
!
      return
   END SUBROUTINE PRI2D
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   SUBROUTINE IPRI2D(IM,JM,F,s,L)
!----------------------------------------------------------------------
!  PRINT 2D array F(IM,JM)
!----------------------------------------------------------------------
!
     implicit none
!
     integer,intent(in) :: IM,JM,L
     real,dimension(IM,JM),intent(in) :: F
     character*(*),intent(in) ::  s
     integer n,i,j
     character*9 sn,sl
!
      write(sl,'(i0)') L
      write(nprint,*)' '
      n=L
      IF (IM.eq.1) THEN
        write(nprint,'("C: '//s//', latitudes from")')
        do j=1,JM,n
          write(sn,'(i0)') min(j+n-1,JM)-j+1
          write(nprint,  &
         '("C:",i4," to",i4,2x,'//trim(sn)//'i'//trim(sl)//')') &
         j,min(j+n-1,JM),F(1,j:min(j+n-1,JM))
        end do
      ELSE IF (JM.eq.1) THEN
        write(nprint,'("C: '//s//', longitudes from")')
        do i=1,IM,n
          write(sn,'(i0)') min(i+n-1,IM)-i+1
          write(nprint, &
        '("C:",i4," to",i4,2x,'//trim(sn)//'i'//trim(sl)//')') &
         i,min(i+n-1,IM),F(i:min(i+n-1,IM),1)
        end do
      ELSE
        do i=1,IM,n
          write(nprint,'("C: '//s//', longitudes from",i4," to",i4)') &
         i,min(i+n-1,IM)
          write(sn,'(i0)') min(i+n-1,IM)-i+1
          do j=1,JM
          write(nprint,'("C:",2i4,2x,'//trim(sn)//'i'//trim(sl)//')') &
           i,j,F(i:min(i+n-1,IM),j)
          end do
        end do
      END IF

      return
   END SUBROUTINE IPRI2D
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SUBROUTINE INVIND(IM,JM,F)
!----------------------------------------------------------------------
!  inverse 2D array F(IM,JM)  S->N or N->S
!----------------------------------------------------------------------
!
     implicit none
!
     integer,intent(in) :: IM,JM
     real(kind_real),dimension(IM,JM),intent(inout) :: F
     integer i,j,j1
     real z
!
      IF (.NOT. INVERTind) RETURN
      j1=JM/2
      do j=1,j1
        do i=1,IM
          z=F(i,j)
          F(i,j)=F(i,JM+1-j)
          F(i,JM+1-j)=z
        end do
      end do
      return
   END SUBROUTINE INVIND
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SUBROUTINE read_restart(iunit)
!----------------------------------------------------------------------
!  inverse 2D array F(IM,JM) from S->N to N->S
!----------------------------------------------------------------------
!
     implicit none
       integer iunit
!
       rewind iunit
!read vars for AM
       read (iunit) t_sst
       read (iunit) fice
       read (iunit) hice
       read (iunit) hsno
!read vars for OM
       read (iunit) t_sfc
       read (iunit) t_bot
       read (iunit) u_bot
       read (iunit) v_bot
       read (iunit) q_bot
       read (iunit) p_bot
       read (iunit) p_surf
       read (iunit) z_bot
       read (iunit) xmu
       read (iunit) flux_dsw
       read (iunit) flux_dlw
       read (iunit) ffmm
       read (iunit) ffhh
       read (iunit) snw
       read (iunit) lprec
!      if (lsout_cc_ocn) then
         read (iunit) dusfc
         read (iunit) dvsfc
         read (iunit) dtsfc
         read (iunit) dqsfc
         read (iunit) flux_lw
         read (iunit) flux_sw
!     endif
      print *,'READ fluxes_for_OM successfully'
      t_fix=t_sfc

      return
      end subroutine read_restart
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SUBROUTINE set_data_init
!----------------------------------------------------------------------
!  read 2D array F(IM,JM) from S->N to N->S for initial start
!----------------------------------------------------------------------
!
     implicit none
     print *,'chk anx,any:',anx,any
     open(91,file='fluxes_init_OM',access='direct', &
       recl=anx*any*8,form='unformatted')
     read(91,rec=1) t_bot
     read(91,rec=2) q_bot
     read(91,rec=3) u_bot
     read(91,rec=4) v_bot
     read(91,rec=5) lprec
     read(91,rec=6) flux_lw
     read(91,rec=7) flux_sw
     close (91)
     t_sfc=t_bot
     t_fix=t_sfc
     t_sst=t_sfc
     p_surf=1.e5
     p_bot=p_surf
     z_bot=10.0
     xmu=0.
     snw=0.
     fice=max(271.2-t_sst,0.0)
     fice=min(fice,0.90)
     hice=3.0*fice
     hsno=fice

      return
   END SUBROUTINE set_data_init
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       subroutine maxmin(xdim,ydim,x,s)

       implicit none

       integer xdim,ydim,i,j
       real(kind=kind_REAL) x(xdim,ydim),xmax,xmin
       character(*) s

      xmax=x(1,1)
      xmin=x(1,1)
      do j=1,ydim
      do i=1,xdim
       if ( xmax .lt. x(i,j) ) xmax=x(i,j)
       if ( xmin .gt. x(i,j) ) xmin=x(i,j)
      enddo
      enddo
      print *,s//' in maxmin,xdim=',xdim,'ydim=',ydim, &
         'xmax=',xmax,'xmin=',xmin

      return
      end subroutine maxmin
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       subroutine lndseaindex

       implicit none

       integer k
       open (21,file='MOM4LND_GFSOCN',form='formatted')
       rewind (21)
       read (21,*) nptsmom4lnd
       print *,'num of ls points: ',nptsmom4lnd

       allocate(nxmomlnd(nptsmom4lnd),nymomlnd(nptsmom4lnd))

       do k=1,nptsmom4lnd
          read (21,'(2i4)') nxmomlnd(k),nymomlnd(k)
       enddo

      return
      end subroutine lndseaindex
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   END MODULE coupler_module
!
!***********************************************************************
