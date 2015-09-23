      module coordinate_def
!     use resol_def
      use machine
      implicit none
      save
       real(kind=kind_evod) , allocatable ::
     . AK5(:),BK5(:),CK5(:),CK(:),DBK(:),bkl(:),   		! hmhj
     . AMHYB(:,:),BMHYB(:,:),SVHYB(:),tor_hyb(:),
     . D_HYB_m(:,:,:),THREF(:),dm205_hyb(:,:,:)			! hmhj
       real(kind=kind_evod) vertcoord_id,eps_si			! hmhj

!
      real(kind=kind_evod) , allocatable :: vcoord(:,:)
      integer nvcoord, idsl, idvc, idvm
!
!sela for semilag version.
       real(kind=kind_evod) , allocatable ::
     . yecm(:,:),tecm(:,:),sv_ecm(:),
     . AM_slg(:,:),BM_slg(:,:),SV_slg(:),tor_slg(:),
     . D_slg_m(:,:,:),y_ecm(:,:),t_ecm(:,:)
!
      end module coordinate_def
