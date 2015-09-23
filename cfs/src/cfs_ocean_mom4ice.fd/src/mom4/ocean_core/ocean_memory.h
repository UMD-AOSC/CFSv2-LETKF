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
! This file contains some commonly used indices 
! needed for setting up the various domains in mom4.  
!
! The STATIC_MEMORY preprocessor option sets these 
! indices according to values specified at compile time, 
! whereas dynamic memory determines the indices at runtime.  
!

#ifdef STATIC_MEMORY

  ! computational domain indices (excludes halos)
  integer, parameter :: isc = 1
  integer, parameter :: iec = NI_LOCAL_
  integer, parameter :: jsc = 1
  integer, parameter :: jec = NJ_LOCAL_

  ! data domain indices (includes halo size of 1) 
  integer, parameter :: isd = 0
  integer, parameter :: ied = NI_LOCAL_ + 1 
  integer, parameter :: jsd = 0
  integer, parameter :: jed = NJ_LOCAL_ + 1

  ! global domain indices
  integer, parameter :: isg = 1
  integer, parameter :: ieg = NI_
  integer, parameter :: jsg = 1
  integer, parameter :: jeg = NJ_
  integer, parameter :: ni = NI_
  integer, parameter :: nj = NJ_
  integer, parameter :: nk  = NK_

#else

  integer :: isd
  integer :: ied
  integer :: jsd
  integer :: jed
  integer :: isc
  integer :: iec
  integer :: jsc
  integer :: jec
  integer :: isg
  integer :: ieg
  integer :: jsg
  integer :: jeg
  integer :: ni
  integer :: nj
  integer :: nk

#endif
