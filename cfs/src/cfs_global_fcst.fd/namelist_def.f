      module namelist_def
      use machine
      implicit none
      save
      integer nszer,nsres,nslwr,nsout,nsswr,nsdfi,nscyc,igen,jo3,ngptc
     &,       lsm,ens_mem,ncw(2),num_reduce,lsea,nsout_hf,dtphys
     &,       levwgt(2)
      real(kind=kind_evod) fhswr,fhlwr,fhrot,fhseg,fhmax,fhout,fhres,
     & fhzer,fhini,fhdfi,fhcyc,filta,crtrh(3),flgmin(2),ref_temp,
     & ccwf(2),dlqf(2),ctei_rm(2),fhgoc3d,fhout_hf,fhmax_hf,cdmbgwd(2),
     & bkgd_vdif_m, bkgd_vdif_h, hdif_fac, psautco(2), prautco(2), evpco
     &,sl_epsln,phigs,phigs_d,bkgd_vdif_s,hdif_fac2,wminco(2)
     &,wgtm(2)
      logical ldiag3d,ras,zhao_mic,sashal,newsas,crick_proof,ccnorm
      logical shal_cnv
      logical mom4ice,mstrat,trans_trac,moist_adj,lggfs3d,cal_pre
      logical climate
      logical lsfwd,lssav,lscca,lsswr,lslwr
!     logical shuff_lats_a,shuff_lats_r,reshuff_lats_a,reshuff_lats_r
      logical shuff_lats_a,shuff_lats_r
      logical hybrid,gen_coord_hybrid,zflxtvd
      logical run_enthalpy, out_virttemp, ndsl
      logical adiab,explicit,pre_rad,random_clds,old_monin,cnvgwd
      logical restart, gfsio_in, gfsio_out
      logical herm_x, herm_y, herm_z
      logical use_ufo

      logical semilag,redgg_a,lingg_a,redgg_b,lingg_b,gg_tracers

!     logical nsst_active
!     logical nsst_restart
!     logical tr_analysis

      integer nst_fcst
      logical nst_spinup

      character*20 ens_nam
!
!     Radiation control parameters
!
      logical norad_precip
      integer isol, ico2, ialb, iems, iaer, iovr_sw, iovr_lw, ictm
!
!  *******************Do we need the following?  *****************************
!     integer max_kdt_switch_1
!     integer max_kdt_switch_2
!     integer , allocatable :: kdt_switch(:,:)
!  ********************************************************

      integer isubc_sw, isubc_lw

      end module namelist_def
