      module d3d_def
      use machine
      implicit none
!
      REAL(KIND=KIND_RAD) ,ALLOCATABLE :: DT3DT(:,:,:,:)
      REAL(KIND=KIND_RAD) ,ALLOCATABLE :: DQ3DT(:,:,:,:)
      REAL(KIND=KIND_RAD) ,ALLOCATABLE :: DU3DT(:,:,:,:)
      REAL(KIND=KIND_RAD) ,ALLOCATABLE :: DV3DT(:,:,:,:)
      REAL(KIND=KIND_RAD) ,ALLOCATABLE :: upd_mf(:,:,:)
     &,                                   dwn_mf(:,:,:)
     &,                                   det_mf(:,:,:)
     &,                                   dkh(:,:,:)
     &,                                   rnp(:,:,:)
      REAL(KIND=KIND_RAD) ,ALLOCATABLE :: CLDCOV(:,:,:)
!
      contains
!
      subroutine d3d_init(lonr,lats_node_r,levs,pl_coeff,
     &                    ldiag3d,lggfs3d)
      implicit none
      integer lonr,lats_node_r,levs,pl_coeff
      logical ldiag3d, lggfs3d
      integer ngpt
!
      if (ldiag3d) then
        allocate (DT3DT(LONR,LEVS,6,lats_node_r))
        allocate (DU3DT(LONR,LEVS,4,lats_node_r))
        allocate (DV3DT(LONR,LEVS,4,lats_node_r))
      else
        allocate (DT3DT(1,1,6,1))
        allocate (DU3DT(1,1,4,1))
        allocate (DV3DT(1,1,4,1))
      endif
      if (ldiag3d .or. lggfs3d) then
        allocate (DQ3DT(LONR,LEVS,5+pl_coeff,lats_node_r))
        allocate (CLDCOV(LONR,LEVS,lats_node_r))
        allocate (upd_mf(LONR,LEVS,lats_node_r))
        allocate (dwn_mf(LONR,LEVS,lats_node_r))
        allocate (det_mf(LONR,LEVS,lats_node_r))
      else
        allocate (DQ3DT(1,1,5+pl_coeff,1))
        allocate (CLDCOV(1,1,1))
        allocate (upd_mf(1,1,1))
        allocate (dwn_mf(1,1,1))
        allocate (det_mf(1,1,1))
      endif
      if (lggfs3d) then
        allocate (dkh(LONR,LEVS,lats_node_r))
        allocate (rnp(LONR,LEVS,lats_node_r))
      else
        allocate (dkh(1,1,1))
        allocate (rnp(1,1,1))
      endif
!
      end subroutine d3d_init
      subroutine d3d_zero(ldiag3d,lggfs3d)
      implicit none
      logical ldiag3d, lggfs3d
      real, parameter :: zero=0.0
!
       if (ldiag3d) then
         DT3DT  = zero
         DU3DT  = zero
         DV3DT  = zero
       endif
       if (ldiag3d .or. lggfs3d) then
         DQ3DT  = zero
         CLDCOV = zero
         upd_mf = zero
         dwn_mf = zero
         det_mf = zero
       endif
       if (lggfs3d) then
         dkh    = zero
         rnp    = zero
       endif
!
      end subroutine d3d_zero
      end module d3d_def
