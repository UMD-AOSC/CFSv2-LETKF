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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!       ROUTINES TO INITIALIZE/FINALIZE MPP MODULE: mpp_init, mpp_exit        !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Modified: Jun Wang
!           Jun.Wang@noaa.gov
  !    subroutine mpp_init( flags, in, out, err, log )
  !      integer, optional, intent(in) :: flags, in, out, err, log
  subroutine mpp_init( flags, localcomm ,myid, cpl_id, cpl_rank,ocean_comm,ocean_npes,pelist)
  integer, optional, intent(in) :: flags
  integer, optional, intent(in) :: localcomm

!--> cpl insertion
  integer, optional, intent(in) :: myid,cpl_id
  integer, optional, intent(inout) :: cpl_rank,ocean_comm,ocean_npes
  integer, dimension(256), optional, intent(inout) :: pelist
!<-- cpl insertion
  integer                       :: my_pe, num_pes, len, i, ocean_id=2,coupler_id=0,coupler_rank=0,ierr
!--> cpl insertion
  integer                       :: ibuf(2),dummy(2),istatus(MPI_STATUS_SIZE)
  integer                       :: color,key
!<-- cpl insertion
  logical                       :: opened
  character(len=5) :: this_pe

  if( module_is_initialized )return

  call MPI_INITIALIZED( opened, error ) !in case called from another MPI package
  if(opened .and. .NOT. PRESENT(localcomm)) call mpp_error( ERROR, 'MPP_INIT: communicator is required' )
  if( .NOT.opened ) then
    call MPI_INIT(error)
    mpp_comm_private = MPI_COMM_WORLD
!--> cpl insertion
    if(present(myid)) ocean_id=myid
    color=Ocean_id
    key=1
    call MPI_COMM_SPLIT( MPI_COMM_WORLD, color,key, ocean_comm, error )
!<-- cpl insertion
 
  else
    mpp_comm_private  = localcomm
  endif
  call MPI_COMM_RANK( mpp_comm_private, pe,   error )
  call MPI_COMM_SIZE( mpp_comm_private, npes, error )
  call MPI_COMM_SIZE( ocean_comm, ocean_npes, error )

!--> cpl insertion
  if( .NOT.opened ) then
   if(ocean_npes .ne. npes) then
    if(present(cpl_id)) coupler_id=cpl_id
    if(present(cpl_rank)) coupler_rank=cpl_rank
    if ( tag .ne.1) call mpp_error(FATAL,'in mpp_init, tag is not 1')
    call MPI_RECV(ibuf,2,MPI_INTEGER,MPI_ANY_SOURCE,tag, &
       MPI_COMM_WORLD,istatus,ierr)
    coupler_rank=ibuf(2)
    if (present(cpl_rank)) cpl_rank=coupler_rank
    if ( coupler_id .ne. ibuf(1)) call mpp_error(FATAL,'in mpp_init, coupler id is not right')
!  
    ibuf(1)=Ocean_id ! Ocean_id
    ibuf(2)=pe       ! global pe rank
    call MPI_GATHER(ibuf,2,MPI_INTEGER,dummy,2,MPI_INTEGER, &
       coupler_rank,MPI_COMM_WORLD,ierr)
   endif
  endif
!<-- cpl insertion

  allocate( request(0:npes-1) )
  request(:) = MPI_REQUEST_NULL

  module_is_initialized = .TRUE.

  !PEsets: make defaults illegal
  peset(:)%count = -1
  peset(:)%id = -1
  peset(:)%group = -1
  !0=single-PE, initialized so that count returns 1
!<-- cpl insertion
  if( npes.eq.1) then
!--> cpl insertion

    peset(0)%count = 1
    allocate( peset(0)%list(1) )
    peset(0)%list = pe
    current_peset_num = 0
    peset(0)%id = mpp_comm_private
  else if(npes.gt.1) then
    peset_num=0
    peset(0)%count = npes
    allocate( peset(0)%list(npes) )
    do i=0,npes-1
      peset(0)%list(i+1) = i
    enddo
    current_peset_num = 0
    peset(0)%id = mpp_comm_private
  endif
!<-- cpl insertion
  call MPI_COMM_GROUP( mpp_comm_private, peset(0)%group, error )
  world_peset_num = get_peset( (/(i,i=0,npes-1)/) )
  current_peset_num = world_peset_num !initialize current PEset to world
!<-- cpl insertion
!set peset 1
  if ( npes.ne.ocean_npes) then      !not stand alone 
    peset_num=1
    peset(1)%id=ocean_comm
    peset(1)%count=ocean_npes
!    call MPI_COMM_rank(ocean_comm,pe,error)
    allocate(peset(1)%list(peset(1)%count))
    call MPI_RECV(peset(1)%list,peset(1)%count,MPI_INTEGER,coupler_rank,1,MPI_COMM_WORLD,istatus,error)
!***
!    peset(1)%list=(/(i,i=0,ocean_npes-1)/)
!***
    call MPI_COMM_GROUP( ocean_comm, peset(1)%group, error )
    pelist(1:peset(1)%count)=peset(1)%list(1:peset(1)%count)
    current_peset_num=1
  endif
!<-- always set peset num to ocean 

  !initialize clocks
  call SYSTEM_CLOCK( count=tick0, count_rate=ticks_per_sec, count_max=max_ticks )
  tick_rate = 1./ticks_per_sec
  clock0 = mpp_clock_id( 'Total runtime', flags=MPP_CLOCK_SYNC )

  if( PRESENT(flags) )then
     debug   = flags.EQ.MPP_DEBUG
     verbose = flags.EQ.MPP_VERBOSE .OR. debug
  end if

  !logunit: log messages are written to configfile.out by default
  etc_unit=get_unit()
  !      write( etcfile,'(a,i4.4)' )trim(etcfile)//'.', pe
  if( pe.EQ.root_pe )open( unit=etc_unit, file=trim(etcfile), status='REPLACE' )
  call mpp_sync()                         
  if( pe.NE.root_pe )open( unit=etc_unit, file=trim(etcfile), status='UNKNOWN' )
!-->cpl deletion  if( pe.NE.root_pe )open( unit=etc_unit, file=trim(etcfile), status='OLD' )
  !if optional argument logunit=stdout, write messages to stdout instead.
  !if specifying non-defaults, you must specify units not yet in use.
  !      if( PRESENT(in) )then
  !          inquire( unit=in, opened=opened )
  !          if( opened )call mpp_error( FATAL, 'MPP_INIT: unable to open stdin.' )
  !          in_unit=in
  !      end if
  !      if( PRESENT(out) )then
  !          inquire( unit=out, opened=opened )
  !          if( opened )call mpp_error( FATAL, 'MPP_INIT: unable to open stdout.' )
  !          out_unit=out
  !      end if
  !      if( PRESENT(err) )then
  !          inquire( unit=err, opened=opened )
  !          if( opened )call mpp_error( FATAL, 'MPP_INIT: unable to open stderr.' )
  !          err_unit=err
  !      end if
  !      log_unit=get_unit()
  !      if( PRESENT(log) )then
  !          inquire( unit=log, opened=opened )
  !          if( opened .AND. log.NE.out_unit )call mpp_error( FATAL, 'MPP_INIT: unable to open stdlog.' )
  !          log_unit=log
  !      end if
  !!log_unit can be written to only from root_pe, all others write to stdout
  !      if( log_unit.NE.out_unit )then
  !          inquire( unit=log_unit, opened=opened )
  !          if( opened )call mpp_error( FATAL, 'MPP_INIT: specified unit for stdlog already in use.' )
  !          if( pe.EQ.root_pe )open( unit=log_unit, file=trim(configfile), status='REPLACE' )
  !          call mpp_sync()
  !          if( pe.NE.root_pe )open( unit=log_unit, file=trim(configfile), status='OLD' )
  !      end if

  call mpp_init_logfile()

  !messages
  if( verbose )call mpp_error( NOTE, 'MPP_INIT: initializing MPP module...' )
  if( pe.EQ.root_pe )then
     write( stdlog(),'(/a)' )'MPP module '//trim(version)//trim(tagname)
     write( stdlog(),'(a,i4)' )'MPP started with NPES=', npes
     write( stdlog(),'(a)' )'Using MPI library for message passing...'
     write( stdlog(), '(a,es12.4,a,i10,a)' ) &
          'Realtime clock resolution=', tick_rate, ' sec (', ticks_per_sec, ' ticks/sec)'
     write( stdlog(), '(a,es12.4,a,i20,a)' ) &
          'Clock rolls over after ', max_ticks*tick_rate, ' sec (', max_ticks, ' ticks)'
     write( stdlog(),'(/a)' )'MPP module '//trim(mpp_parameter_version)//trim(mpp_parameter_tagname)
     write( stdlog(),'(/a)' )'MPP module '//trim(mpp_data_version)//trim(mpp_data_tagname)
     write( stdlog(),'(/a)' )'MPP module '//trim(mpp_util_version)//trim(mpp_util_tagname)
  end if
  call mpp_clock_begin(clock0)

  return
end subroutine mpp_init

!#######################################################################
  !to be called at the end of a run
subroutine mpp_exit()
  integer :: i, j, k, n, nmax, istat
  real    :: t, tmin, tmax, tavg, tstd
  real    :: m, mmin, mmax, mavg, mstd, t_total
!-->cpl insertion
  integer :: worldpeset,penumm1,minpe,minpe1
  logical :: no_sync=.true.
!<--cpl insertion

  if( .NOT.module_is_initialized )return
  call mpp_set_current_pelist()
  call mpp_clock_end(clock0)
  t_total = clocks(clock0)%total_ticks*tick_rate
  if( clock_num.GT.0 )then
     if( ANY(clocks(1:clock_num)%detailed) )then
        call sum_clock_data; call dump_clock_summary
     end if
     if( pe.EQ.root_pe )then
        write( stdout(),'(/a,i4,a)' ) 'Tabulating mpp_clock statistics across ', npes, ' PEs...'
        if( ANY(clocks(1:clock_num)%detailed) ) &
             write( stdout(),'(a)' )'   ... see mpp_clock.out.#### for details on individual PEs.'
        write( stdout(),'(/32x,a)' ) '          tmin          tmax          tavg          tstd  tfrac grain pemin pemax'
     end if
     call FLUSH( stdout(), istat )
     call mpp_sync()
     do i = 1,clock_num
!-->cpl insertion
       if(i.eq.1) then
         penumm1=1
         minpe=peset(1)%list(1)    ! this pelist has been sorted
       else
         penumm1=clocks(i-1)%peset_num
         minpe=clocks(i-1)%peset_num
         minpe1=clocks(i)%peset_num
       endif
       if (clocks(i)%peset_num.eq.2 .and. penumm1.eq.1) then
         if( .NOT.ANY(peset(clocks(i)%peset_num)%list(:).EQ.pe-minpe) )cycle
       else if (  clocks(i)%peset_num.eq.1 .and. penumm1.eq.2) then
         if( .NOT.ANY(peset(clocks(i)%peset_num)%list(:).EQ.pe+minpe1) )cycle
       else
         if( .NOT.ANY(peset(clocks(i)%peset_num)%list(:).EQ.pe) )cycle
       endif
!<--cpl insertion
   
        call mpp_set_current_pelist( peset(clocks(i)%peset_num)%list )
        !times between mpp_clock ticks
        t = clocks(i)%total_ticks*tick_rate
        tmin = t; call mpp_min(tmin)
        tmax = t; call mpp_max(tmax)
        tavg = t; call mpp_sum(tavg); tavg = tavg/mpp_npes()
        tstd = (t-tavg)**2; call mpp_sum(tstd); tstd = sqrt( tstd/mpp_npes() )
        if( pe.EQ.root_pe )write( stdout(),'(a32,4f14.6,f7.3,3i6)' ) &
             clocks(i)%name, tmin, tmax, tavg, tstd, tavg/t_total, &
             clocks(i)%grain, minval(peset(clocks(i)%peset_num)%list), &
             maxval(peset(clocks(i)%peset_num)%list)
     end do
     if( ANY(clocks(1:clock_num)%detailed) .AND. pe.EQ.root_pe )write( stdout(),'(/32x,a)' ) &
          '       tmin       tmax       tavg       tstd       mmin       mmax       mavg       mstd  mavg/tavg'
    
  
     do i = 1,clock_num
        !messages: bytelengths and times
        if( .NOT.clocks(i)%detailed )cycle
        if( .NOT.ANY(peset(clocks(i)%peset_num)%list(:).EQ.pe) )cycle
        call mpp_set_current_pelist( peset(clocks(i)%peset_num)%list )
        do j = 1,MAX_EVENT_TYPES
           n = clocks(i)%events(j)%calls; nmax = n
           call mpp_max(nmax)
           if( nmax.NE.0 )then
              !don't divide by n because n might be 0
              m = 0
              if( n.GT.0 )m = sum(clocks(i)%events(j)%bytes(1:n))
              mmin = m; call mpp_min(mmin)
              mmax = m; call mpp_max(mmax)
              mavg = m; call mpp_sum(mavg); mavg = mavg/mpp_npes()
              mstd = (m-mavg)**2; call mpp_sum(mstd); mstd = sqrt( mstd/mpp_npes() )
              t = 0
              if( n.GT.0 )t = sum(clocks(i)%events(j)%ticks(1:n))*tick_rate
              tmin = t; call mpp_min(tmin)
              tmax = t; call mpp_max(tmax)
              tavg = t; call mpp_sum(tavg); tavg = tavg/mpp_npes()
              tstd = (t-tavg)**2; call mpp_sum(tstd); tstd = sqrt( tstd/mpp_npes() )
              if( pe.EQ.root_pe )write( stdout(),'(a32,4f11.3,5es11.3)' ) &
                   trim(clocks(i)%name)//' '//trim(clocks(i)%events(j)%name), &
                   tmin, tmax, tavg, tstd, mmin, mmax, mavg, mstd, mavg/tavg
           end if
        end do
     end do
  end if

  call mpp_max(mpp_stack_hwm)
  if( pe.EQ.root_pe )write( stdout(),* )'MPP_STACK high water mark=', mpp_stack_hwm
  call mpp_set_current_pelist(pesetnum=worldpeset,no_sync=no_sync)
!  call mpp_sync()
!     print *,'in mpp_exit, after  mpp_sync'
     print *,'in mpp_exit, after high water mark'
  if(mpp_comm_private == MPI_COMM_WORLD ) call MPI_FINALIZE(error)



  return
end subroutine mpp_exit

!#######################################################################
 !--- dummy routine, should never been called -------------------------
    subroutine mpp_malloc( ptr, newlen, len )
      integer, intent(in)    :: newlen
      integer, intent(inout) :: len
      real                   :: dummy
      pointer( ptr, dummy )

      call mpp_error(FATAL,'mpp_malloc: Should not been called when use_libMPI is defined')

      return
    end subroutine mpp_malloc

#ifdef use_MPI_GSM
!#######################################################################
 !--- routine to perform GSM allocations
 !
subroutine mpp_gsm_malloc( ptr, len )
  use mpp_data_mod,      only : peset, current_peset_num
  integer(KIND=MPI_ADDRESS_KIND), intent(in)    :: len
  real                   :: dummy
  integer                :: ierror
  !argument ptr is a cray pointer, points to a dummy argument in this routine 
  pointer( ptr, dummy )
  include "mpi_gsmf.h"
  if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_GSM_MALLOC: You must first call mpp_init.' )

  call gsm_alloc(len, GSM_SINGLERANK, 0, peset(current_peset_num)%id, ptr, ierror)
  if (ierror .eq. -1) call mpp_error( FATAL, 'MPP_GSM_MALLOC: gsm_alloc failed.' )
  return
end subroutine mpp_gsm_malloc

!#######################################################################
 !--- routine to free GSM allocations
 !
subroutine mpp_gsm_free( ptr )
  use mpp_data_mod,      only : peset, current_peset_num
  real                   :: dummy
  integer                :: ierror
  !argument ptr is a cray pointer, points to a dummy argument in this routine 
  pointer( ptr, dummy )
  include "mpi_gsmf.h"

  call gsm_free(dummy, peset(current_peset_num)%id, ierror)
  if (ierror .eq. -1) call mpp_error( FATAL, 'MPP_GSM_FREE: gsm_free failed.' )
  return
end subroutine mpp_gsm_free
#endif

!#######################################################################
    !set the mpp_stack variable to be at least n LONG words long
  subroutine mpp_set_stack_size(n)
    integer, intent(in) :: n
    character(len=8)    :: text

    if( n.GT.mpp_stack_size .AND. allocated(mpp_stack) )deallocate(mpp_stack)
    if( .NOT.allocated(mpp_stack) )then
       allocate( mpp_stack(n) )
       mpp_stack_size = n
    end if

    write( text,'(i8)' )n
    if( pe.EQ.root_pe )call mpp_error( NOTE, 'MPP_SET_STACK_SIZE: stack size set to '//text//'.' )

    return
  end subroutine mpp_set_stack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                BASIC MESSAGE PASSING ROUTINE: mpp_transmit                  !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define MPP_TRANSMIT_ mpp_transmit_real8
#define MPP_TRANSMIT_SCALAR_ mpp_transmit_real8_scalar
#define MPP_TRANSMIT_2D_ mpp_transmit_real8_2d
#define MPP_TRANSMIT_3D_ mpp_transmit_real8_3d
#define MPP_TRANSMIT_4D_ mpp_transmit_real8_4d
#define MPP_TRANSMIT_5D_ mpp_transmit_real8_5d
#define MPP_RECV_ mpp_recv_real8
#define MPP_RECV_SCALAR_ mpp_recv_real8_scalar
#define MPP_RECV_2D_ mpp_recv_real8_2d
#define MPP_RECV_3D_ mpp_recv_real8_3d
#define MPP_RECV_4D_ mpp_recv_real8_4d
#define MPP_RECV_5D_ mpp_recv_real8_5d
#define MPP_SEND_ mpp_send_real8
#define MPP_SEND_SCALAR_ mpp_send_real8_scalar
#define MPP_SEND_2D_ mpp_send_real8_2d
#define MPP_SEND_3D_ mpp_send_real8_3d
#define MPP_SEND_4D_ mpp_send_real8_4d
#define MPP_SEND_5D_ mpp_send_real8_5d
#define MPP_BROADCAST_ mpp_broadcast_real8
#define MPP_BROADCAST_SCALAR_ mpp_broadcast_real8_scalar
#define MPP_BROADCAST_2D_ mpp_broadcast_real8_2d
#define MPP_BROADCAST_3D_ mpp_broadcast_real8_3d
#define MPP_BROADCAST_4D_ mpp_broadcast_real8_4d
#define MPP_BROADCAST_5D_ mpp_broadcast_real8_5d
#define MPP_TYPE_ real(DOUBLE_KIND)
#define MPP_TYPE_BYTELEN_ 8
#define MPI_TYPE_ MPI_REAL8
#include <mpp_transmit_mpi.h>

#define MPP_TRANSMIT_ mpp_transmit_cmplx8
#define MPP_TRANSMIT_SCALAR_ mpp_transmit_cmplx8_scalar
#define MPP_TRANSMIT_2D_ mpp_transmit_cmplx8_2d
#define MPP_TRANSMIT_3D_ mpp_transmit_cmplx8_3d
#define MPP_TRANSMIT_4D_ mpp_transmit_cmplx8_4d
#define MPP_TRANSMIT_5D_ mpp_transmit_cmplx8_5d
#define MPP_RECV_ mpp_recv_cmplx8
#define MPP_RECV_SCALAR_ mpp_recv_cmplx8_scalar
#define MPP_RECV_2D_ mpp_recv_cmplx8_2d
#define MPP_RECV_3D_ mpp_recv_cmplx8_3d
#define MPP_RECV_4D_ mpp_recv_cmplx8_4d
#define MPP_RECV_5D_ mpp_recv_cmplx8_5d
#define MPP_SEND_ mpp_send_cmplx8
#define MPP_SEND_SCALAR_ mpp_send_cmplx8_scalar
#define MPP_SEND_2D_ mpp_send_cmplx8_2d
#define MPP_SEND_3D_ mpp_send_cmplx8_3d
#define MPP_SEND_4D_ mpp_send_cmplx8_4d
#define MPP_SEND_5D_ mpp_send_cmplx8_5d
#define MPP_BROADCAST_ mpp_broadcast_cmplx8
#define MPP_BROADCAST_SCALAR_ mpp_broadcast_cmplx8_scalar
#define MPP_BROADCAST_2D_ mpp_broadcast_cmplx8_2d
#define MPP_BROADCAST_3D_ mpp_broadcast_cmplx8_3d
#define MPP_BROADCAST_4D_ mpp_broadcast_cmplx8_4d
#define MPP_BROADCAST_5D_ mpp_broadcast_cmplx8_5d
#define MPP_TYPE_ complex(DOUBLE_KIND)
#define MPP_TYPE_BYTELEN_ 16
#define MPI_TYPE_ MPI_DOUBLE_COMPLEX
#include <mpp_transmit_mpi.h>

#ifndef no_4byte_reals
#define MPP_TRANSMIT_ mpp_transmit_real4
#define MPP_TRANSMIT_SCALAR_ mpp_transmit_real4_scalar
#define MPP_TRANSMIT_2D_ mpp_transmit_real4_2d
#define MPP_TRANSMIT_3D_ mpp_transmit_real4_3d
#define MPP_TRANSMIT_4D_ mpp_transmit_real4_4d
#define MPP_TRANSMIT_5D_ mpp_transmit_real4_5d
#define MPP_RECV_ mpp_recv_real4
#define MPP_RECV_SCALAR_ mpp_recv_real4_scalar
#define MPP_RECV_2D_ mpp_recv_real4_2d
#define MPP_RECV_3D_ mpp_recv_real4_3d
#define MPP_RECV_4D_ mpp_recv_real4_4d
#define MPP_RECV_5D_ mpp_recv_real4_5d
#define MPP_SEND_ mpp_send_real4
#define MPP_SEND_SCALAR_ mpp_send_real4_scalar
#define MPP_SEND_2D_ mpp_send_real4_2d
#define MPP_SEND_3D_ mpp_send_real4_3d
#define MPP_SEND_4D_ mpp_send_real4_4d
#define MPP_SEND_5D_ mpp_send_real4_5d
#define MPP_BROADCAST_ mpp_broadcast_real4
#define MPP_BROADCAST_SCALAR_ mpp_broadcast_real4_scalar
#define MPP_BROADCAST_2D_ mpp_broadcast_real4_2d
#define MPP_BROADCAST_3D_ mpp_broadcast_real4_3d
#define MPP_BROADCAST_4D_ mpp_broadcast_real4_4d
#define MPP_BROADCAST_5D_ mpp_broadcast_real4_5d
#define MPP_TYPE_ real(FLOAT_KIND)
#define MPP_TYPE_BYTELEN_ 4
#define MPI_TYPE_ MPI_REAL4
#include <mpp_transmit_mpi.h>

#define MPP_TRANSMIT_ mpp_transmit_cmplx4
#define MPP_TRANSMIT_SCALAR_ mpp_transmit_cmplx4_scalar
#define MPP_TRANSMIT_2D_ mpp_transmit_cmplx4_2d
#define MPP_TRANSMIT_3D_ mpp_transmit_cmplx4_3d
#define MPP_TRANSMIT_4D_ mpp_transmit_cmplx4_4d
#define MPP_TRANSMIT_5D_ mpp_transmit_cmplx4_5d
#define MPP_RECV_ mpp_recv_cmplx4
#define MPP_RECV_SCALAR_ mpp_recv_cmplx4_scalar
#define MPP_RECV_2D_ mpp_recv_cmplx4_2d
#define MPP_RECV_3D_ mpp_recv_cmplx4_3d
#define MPP_RECV_4D_ mpp_recv_cmplx4_4d
#define MPP_RECV_5D_ mpp_recv_cmplx4_5d
#define MPP_SEND_ mpp_send_cmplx4
#define MPP_SEND_SCALAR_ mpp_send_cmplx4_scalar
#define MPP_SEND_2D_ mpp_send_cmplx4_2d
#define MPP_SEND_3D_ mpp_send_cmplx4_3d
#define MPP_SEND_4D_ mpp_send_cmplx4_4d
#define MPP_SEND_5D_ mpp_send_cmplx4_5d
#define MPP_BROADCAST_ mpp_broadcast_cmplx4
#define MPP_BROADCAST_SCALAR_ mpp_broadcast_cmplx4_scalar
#define MPP_BROADCAST_2D_ mpp_broadcast_cmplx4_2d
#define MPP_BROADCAST_3D_ mpp_broadcast_cmplx4_3d
#define MPP_BROADCAST_4D_ mpp_broadcast_cmplx4_4d
#define MPP_BROADCAST_5D_ mpp_broadcast_cmplx4_5d
#define MPP_TYPE_ complex(FLOAT_KIND)
#define MPP_TYPE_BYTELEN_ 8
#define MPI_TYPE_ MPI_COMPLEX
#include <mpp_transmit_mpi.h>
#endif 

#ifndef no_8byte_integers
#define MPP_TRANSMIT_ mpp_transmit_int8
#define MPP_TRANSMIT_SCALAR_ mpp_transmit_int8_scalar
#define MPP_TRANSMIT_2D_ mpp_transmit_int8_2d
#define MPP_TRANSMIT_3D_ mpp_transmit_int8_3d
#define MPP_TRANSMIT_4D_ mpp_transmit_int8_4d
#define MPP_TRANSMIT_5D_ mpp_transmit_int8_5d
#define MPP_RECV_ mpp_recv_int8
#define MPP_RECV_SCALAR_ mpp_recv_int8_scalar
#define MPP_RECV_2D_ mpp_recv_int8_2d
#define MPP_RECV_3D_ mpp_recv_int8_3d
#define MPP_RECV_4D_ mpp_recv_int8_4d
#define MPP_RECV_5D_ mpp_recv_int8_5d
#define MPP_SEND_ mpp_send_int8
#define MPP_SEND_SCALAR_ mpp_send_int8_scalar
#define MPP_SEND_2D_ mpp_send_int8_2d
#define MPP_SEND_3D_ mpp_send_int8_3d
#define MPP_SEND_4D_ mpp_send_int8_4d
#define MPP_SEND_5D_ mpp_send_int8_5d
#define MPP_BROADCAST_ mpp_broadcast_int8
#define MPP_BROADCAST_SCALAR_ mpp_broadcast_int8_scalar
#define MPP_BROADCAST_2D_ mpp_broadcast_int8_2d
#define MPP_BROADCAST_3D_ mpp_broadcast_int8_3d
#define MPP_BROADCAST_4D_ mpp_broadcast_int8_4d
#define MPP_BROADCAST_5D_ mpp_broadcast_int8_5d
#define MPP_TYPE_ integer(LONG_KIND)
#define MPP_TYPE_BYTELEN_ 8
#define MPI_TYPE_ MPI_INTEGER8
#include <mpp_transmit_mpi.h>
#endif

#define MPP_TRANSMIT_ mpp_transmit_int4
#define MPP_TRANSMIT_SCALAR_ mpp_transmit_int4_scalar
#define MPP_TRANSMIT_2D_ mpp_transmit_int4_2d
#define MPP_TRANSMIT_3D_ mpp_transmit_int4_3d
#define MPP_TRANSMIT_4D_ mpp_transmit_int4_4d
#define MPP_TRANSMIT_5D_ mpp_transmit_int4_5d
#define MPP_RECV_ mpp_recv_int4
#define MPP_RECV_SCALAR_ mpp_recv_int4_scalar
#define MPP_RECV_2D_ mpp_recv_int4_2d
#define MPP_RECV_3D_ mpp_recv_int4_3d
#define MPP_RECV_4D_ mpp_recv_int4_4d
#define MPP_RECV_5D_ mpp_recv_int4_5d
#define MPP_SEND_ mpp_send_int4
#define MPP_SEND_SCALAR_ mpp_send_int4_scalar
#define MPP_SEND_2D_ mpp_send_int4_2d
#define MPP_SEND_3D_ mpp_send_int4_3d
#define MPP_SEND_4D_ mpp_send_int4_4d
#define MPP_SEND_5D_ mpp_send_int4_5d
#define MPP_BROADCAST_ mpp_broadcast_int4
#define MPP_BROADCAST_SCALAR_ mpp_broadcast_int4_scalar
#define MPP_BROADCAST_2D_ mpp_broadcast_int4_2d
#define MPP_BROADCAST_3D_ mpp_broadcast_int4_3d
#define MPP_BROADCAST_4D_ mpp_broadcast_int4_4d
#define MPP_BROADCAST_5D_ mpp_broadcast_int4_5d
#define MPP_TYPE_ integer(INT_KIND)
#define MPP_TYPE_BYTELEN_ 4
#define MPI_TYPE_ MPI_INTEGER4
#include <mpp_transmit_mpi.h>

#ifndef no_8byte_integers
#define MPP_TRANSMIT_ mpp_transmit_logical8
#define MPP_TRANSMIT_SCALAR_ mpp_transmit_logical8_scalar
#define MPP_TRANSMIT_2D_ mpp_transmit_logical8_2d
#define MPP_TRANSMIT_3D_ mpp_transmit_logical8_3d
#define MPP_TRANSMIT_4D_ mpp_transmit_logical8_4d
#define MPP_TRANSMIT_5D_ mpp_transmit_logical8_5d
#define MPP_RECV_ mpp_recv_logical8
#define MPP_RECV_SCALAR_ mpp_recv_logical8_scalar
#define MPP_RECV_2D_ mpp_recv_logical8_2d
#define MPP_RECV_3D_ mpp_recv_logical8_3d
#define MPP_RECV_4D_ mpp_recv_logical8_4d
#define MPP_RECV_5D_ mpp_recv_logical8_5d
#define MPP_SEND_ mpp_send_logical8
#define MPP_SEND_SCALAR_ mpp_send_logical8_scalar
#define MPP_SEND_2D_ mpp_send_logical8_2d
#define MPP_SEND_3D_ mpp_send_logical8_3d
#define MPP_SEND_4D_ mpp_send_logical8_4d
#define MPP_SEND_5D_ mpp_send_logical8_5d
#define MPP_BROADCAST_ mpp_broadcast_logical8
#define MPP_BROADCAST_SCALAR_ mpp_broadcast_logical8_scalar
#define MPP_BROADCAST_2D_ mpp_broadcast_logical8_2d
#define MPP_BROADCAST_3D_ mpp_broadcast_logical8_3d
#define MPP_BROADCAST_4D_ mpp_broadcast_logical8_4d
#define MPP_BROADCAST_5D_ mpp_broadcast_logical8_5d
#define MPP_TYPE_ logical(LONG_KIND)
#define MPP_TYPE_BYTELEN_ 8
#define MPI_TYPE_ MPI_INTEGER8
#include <mpp_transmit_mpi.h>
#endif

#define MPP_TRANSMIT_ mpp_transmit_logical4
#define MPP_TRANSMIT_SCALAR_ mpp_transmit_logical4_scalar
#define MPP_TRANSMIT_2D_ mpp_transmit_logical4_2d
#define MPP_TRANSMIT_3D_ mpp_transmit_logical4_3d
#define MPP_TRANSMIT_4D_ mpp_transmit_logical4_4d
#define MPP_TRANSMIT_5D_ mpp_transmit_logical4_5d
#define MPP_RECV_ mpp_recv_logical4
#define MPP_RECV_SCALAR_ mpp_recv_logical4_scalar
#define MPP_RECV_2D_ mpp_recv_logical4_2d
#define MPP_RECV_3D_ mpp_recv_logical4_3d
#define MPP_RECV_4D_ mpp_recv_logical4_4d
#define MPP_RECV_5D_ mpp_recv_logical4_5d
#define MPP_SEND_ mpp_send_logical4
#define MPP_SEND_SCALAR_ mpp_send_logical4_scalar
#define MPP_SEND_2D_ mpp_send_logical4_2d
#define MPP_SEND_3D_ mpp_send_logical4_3d
#define MPP_SEND_4D_ mpp_send_logical4_4d
#define MPP_SEND_5D_ mpp_send_logical4_5d
#define MPP_BROADCAST_ mpp_broadcast_logical4
#define MPP_BROADCAST_SCALAR_ mpp_broadcast_logical4_scalar
#define MPP_BROADCAST_2D_ mpp_broadcast_logical4_2d
#define MPP_BROADCAST_3D_ mpp_broadcast_logical4_3d
#define MPP_BROADCAST_4D_ mpp_broadcast_logical4_4d
#define MPP_BROADCAST_5D_ mpp_broadcast_logical4_5d
#define MPP_TYPE_ logical(INT_KIND)
#define MPP_TYPE_BYTELEN_ 4
#define MPI_TYPE_ MPI_INTEGER4
#include <mpp_transmit_mpi.h>

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!            GLOBAL REDUCTION ROUTINES: mpp_max, mpp_sum, mpp_min             !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define MPP_REDUCE_ mpp_max_real8
#define MPP_TYPE_ real(DOUBLE_KIND)
#define MPP_TYPE_BYTELEN_ 8
#define MPI_TYPE_ MPI_REAL8
#if defined(use_GSM)
#define GSM_REDUCE_ .GT.
#include <mpp_reduce_gsm.h>
#elif defined(use_CAF)
#define CAF_REDUCE_ .GT.
#define CAFPNTR_TYPE_1D_ cafptr_r8_1d_type
#include <mpp_reduce_caf.h>
#else
#define MPI_REDUCE_ MPI_MAX
#include <mpp_reduce_mpi.h>
#endif

#ifndef no_4byte_reals
#define MPP_REDUCE_ mpp_max_real4
#define MPP_TYPE_ real(FLOAT_KIND)
#define MPP_TYPE_BYTELEN_ 4
#define MPI_TYPE_ MPI_REAL4
#if defined(use_GSM)
#define GSM_REDUCE_ .GT.
#include <mpp_reduce_gsm.h>
#elif defined(use_CAF) 
#define CAF_REDUCE_ .GT.
#define CAFPNTR_TYPE_1D_ cafptr_r4_1d_type
#include <mpp_reduce_caf.h>
#else
#define MPI_REDUCE_ MPI_MAX
#include <mpp_reduce_mpi.h>
#endif
#endif

#ifndef no_8byte_integers
#define MPP_REDUCE_ mpp_max_int8
#define MPP_TYPE_ integer(LONG_KIND)
#define MPP_TYPE_BYTELEN_ 8
#define MPI_TYPE_ MPI_INTEGER8
#if defined(use_GSM)
#define GSM_REDUCE_ .GT.
#include <mpp_reduce_gsm.h>
#elif defined(use_CAF) 
#define CAF_REDUCE_ .GT.
#define CAFPNTR_TYPE_1D_ cafptr_i8_1d_type
#include <mpp_reduce_caf.h>
#else
#define MPI_REDUCE_ MPI_MAX
#include <mpp_reduce_mpi.h>
#endif
#endif

#define MPP_REDUCE_ mpp_max_int4
#define MPP_TYPE_ integer(INT_KIND)
#define MPP_TYPE_BYTELEN_ 4
#define MPI_TYPE_ MPI_INTEGER4
#if defined(use_GSM)
#define GSM_REDUCE_ .GT.
#include <mpp_reduce_gsm.h>
#elif defined(use_CAF) 
#define CAF_REDUCE_ .GT.
#define CAFPNTR_TYPE_1D_ cafptr_i4_1d_type
#include <mpp_reduce_caf.h>
#else
#define MPI_REDUCE_ MPI_MAX
#include <mpp_reduce_mpi.h>
#endif

#define MPP_REDUCE_ mpp_min_real8
#define MPP_TYPE_ real(DOUBLE_KIND)
#define MPP_TYPE_BYTELEN_ 8
#define MPI_TYPE_ MPI_REAL8
#if defined(use_GSM)
#define GSM_REDUCE_ .LT.
#include <mpp_reduce_gsm.h>
#elif defined(use_CAF) 
#define CAF_REDUCE_ .LT.
#define CAFPNTR_TYPE_1D_ cafptr_r8_1d_type
#include <mpp_reduce_caf.h>
#else
#define MPI_REDUCE_ MPI_MIN
#include <mpp_reduce_mpi.h>
#endif

#ifndef no_4byte_reals
#define MPP_REDUCE_ mpp_min_real4
#define MPP_TYPE_ real(FLOAT_KIND)
#define MPP_TYPE_BYTELEN_ 4
#define MPI_TYPE_ MPI_REAL4
#if defined(use_GSM)
#define GSM_REDUCE_ .LT.
#include <mpp_reduce_gsm.h>
#elif defined(use_CAF) 
#define CAF_REDUCE_ .LT.
#define CAFPNTR_TYPE_1D_ cafptr_r4_1d_type
#include <mpp_reduce_caf.h>
#else
#define MPI_REDUCE_ MPI_MIN
#include <mpp_reduce_mpi.h>
#endif
#endif

#ifndef no_8byte_integers
#define MPP_REDUCE_ mpp_min_int8
#define MPP_TYPE_ integer(LONG_KIND)
#define MPP_TYPE_BYTELEN_ 8
#define MPI_TYPE_ MPI_INTEGER8
#if defined(use_GSM)
#define GSM_REDUCE_ .LT.
#include <mpp_reduce_gsm.h>
#elif defined(use_CAF) 
#define CAF_REDUCE_ .LT.
#define CAFPNTR_TYPE_1D_ cafptr_i8_1d_type
#include <mpp_reduce_caf.h>
#else
#define MPI_REDUCE_ MPI_MIN
#include <mpp_reduce_mpi.h>
#endif       
#endif       
             
#define MPP_REDUCE_ mpp_min_int4
#define MPP_TYPE_ integer(INT_KIND)
#define MPP_TYPE_BYTELEN_ 4
#define MPI_TYPE_ MPI_INTEGER4
#if defined(use_GSM)
#define GSM_REDUCE_ .LT.
#include <mpp_reduce_gsm.h>
#elif defined(use_CAF) 
#define CAF_REDUCE_ .LT.
#define CAFPNTR_TYPE_1D_ cafptr_i4_1d_type
#include <mpp_reduce_caf.h>
#else
#define MPI_REDUCE_ MPI_MIN
#include <mpp_reduce_mpi.h>
#endif

#define MPP_SUM_ mpp_sum_real8
#define MPP_SUM_SCALAR_ mpp_sum_real8_scalar
#define MPP_SUM_2D_ mpp_sum_real8_2d
#define MPP_SUM_3D_ mpp_sum_real8_3d
#define MPP_SUM_4D_ mpp_sum_real8_4d
#define MPP_SUM_5D_ mpp_sum_real8_5d
#define MPP_TYPE_ real(DOUBLE_KIND)
#define MPI_TYPE_ MPI_REAL8
#define MPP_TYPE_BYTELEN_ 8
#if defined(use_GSM)
#include <mpp_sum_gsm.h>
#elif defined(use_CAF) 
#define CAFPNTR_TYPE_1D_ cafptr_r8_1d_type
#include <mpp_sum_caf.h>
#else
#include <mpp_sum_mpi.h>
#endif

#define MPP_SUM_ mpp_sum_cmplx8
#define MPP_SUM_SCALAR_ mpp_sum_cmplx8_scalar
#define MPP_SUM_2D_ mpp_sum_cmplx8_2d
#define MPP_SUM_3D_ mpp_sum_cmplx8_3d
#define MPP_SUM_4D_ mpp_sum_cmplx8_4d
#define MPP_SUM_5D_ mpp_sum_cmplx8_5d
#define MPP_TYPE_ complex(DOUBLE_KIND)
#define MPI_TYPE_ MPI_DOUBLE_COMPLEX
#define MPP_TYPE_BYTELEN_ 16
#if defined(use_GSM)
#include <mpp_sum_gsm.h>
#elif defined(use_CAF) 
#define CAFPNTR_TYPE_1D_ cafptr_c8_1d_type
#include <mpp_sum_caf.h>
#else
#include <mpp_sum_mpi.h>
#endif

#ifndef no_4byte_reals
#define MPP_SUM_ mpp_sum_real4
#define MPP_SUM_SCALAR_ mpp_sum_real4_scalar
#define MPP_SUM_2D_ mpp_sum_real4_2d
#define MPP_SUM_3D_ mpp_sum_real4_3d
#define MPP_SUM_4D_ mpp_sum_real4_4d
#define MPP_SUM_5D_ mpp_sum_real4_5d
#define MPP_TYPE_ real(FLOAT_KIND)
#define MPI_TYPE_ MPI_REAL4
#define MPP_TYPE_BYTELEN_ 4
#if defined(use_GSM)
#include <mpp_sum_gsm.h>
#elif defined(use_CAF) 
#define CAFPNTR_TYPE_1D_ cafptr_r4_1d_type
#include <mpp_sum_caf.h>
#else
#include <mpp_sum_mpi.h>
#endif

#define MPP_SUM_ mpp_sum_cmplx4
#define MPP_SUM_SCALAR_ mpp_sum_cmplx4_scalar
#define MPP_SUM_2D_ mpp_sum_cmplx4_2d
#define MPP_SUM_3D_ mpp_sum_cmplx4_3d
#define MPP_SUM_4D_ mpp_sum_cmplx4_4d
#define MPP_SUM_5D_ mpp_sum_cmplx4_5d
#define MPP_TYPE_ complex(FLOAT_KIND)
#define MPI_TYPE_ MPI_COMPLEX
#define MPP_TYPE_BYTELEN_ 8
#if defined(use_GSM)
#include <mpp_sum_gsm.h>
#elif defined(use_CAF) 
#define CAFPNTR_TYPE_1D_ cafptr_c4_1d_type
#include <mpp_sum_caf.h>
#else
#include <mpp_sum_mpi.h>
#endif
#endif

#ifndef no_8byte_integers
#define MPP_SUM_ mpp_sum_int8
#define MPP_SUM_SCALAR_ mpp_sum_int8_scalar
#define MPP_SUM_2D_ mpp_sum_int8_2d
#define MPP_SUM_3D_ mpp_sum_int8_3d
#define MPP_SUM_4D_ mpp_sum_int8_4d
#define MPP_SUM_5D_ mpp_sum_int8_5d
#define MPP_TYPE_ integer(LONG_KIND)
#define MPI_TYPE_ MPI_INTEGER8
#define MPP_TYPE_BYTELEN_ 8
#if defined(use_GSM)
#include <mpp_sum_gsm.h>
#elif defined(use_CAF) 
#define CAFPNTR_TYPE_1D_ cafptr_i8_1d_type
#include <mpp_sum_caf.h>
#else
#include <mpp_sum_mpi.h>
#endif
#endif

#define MPP_SUM_ mpp_sum_int4
#define MPP_SUM_SCALAR_ mpp_sum_int4_scalar
#define MPP_SUM_2D_ mpp_sum_int4_2d
#define MPP_SUM_3D_ mpp_sum_int4_3d
#define MPP_SUM_4D_ mpp_sum_int4_4d
#define MPP_SUM_5D_ mpp_sum_int4_5d
#define MPP_TYPE_ integer(INT_KIND)
#define MPI_TYPE_ MPI_INTEGER4
#define MPP_TYPE_BYTELEN_ 4
#if defined(use_GSM)
#include <mpp_sum_gsm.h>
#elif defined(use_CAF) 
#define CAFPNTR_TYPE_1D_ cafptr_i4_1d_type
#include <mpp_sum_caf.h>
#else
#include <mpp_sum_mpi.h>
#endif
