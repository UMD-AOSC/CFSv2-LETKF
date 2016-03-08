####   make with mpxlf_r (threading) 
EXEC	 = global_fcst
#F77	 = $(F77M)   
#F90	 = $(F90M)   
#CC	 = $(FCC)
#
FINCS    = $(FINC)

ARCH     = $(ARCHM)
PGSZ     = $(PGSZM)
FRR 	 = $(FRRM)
FXX 	 = $(FXXM)
OPTS     = $(OPTSM)
OPTS_SER = $(OPTS_SERM)
OPTS90   = $(OPTS90M)
OPTS90A  = $(OPTS90AM)
OPTSIO	 = $(OPTSIOM)
FFLAG90  = $(OPTS90) $(FINCS) $(FINCM) $(FRR)
FFLAG90A = $(OPTS90A) $(FINCS) $(FRR)
FFLAGS   = $(OPTS) $(TRAPS) $(FXX) 
FFLAGX   = $(OPTS) $(TRAPS) $(FXX)
FFLAGIO  = $(OPTSIO) $(TRAPS) $(FXX)
FFLAGY   = $(OPTS) $(FXX)
FFLAGM   = $(OPTS) $(TRAPS) $(DEBUG) $(FXX) 
FFLAGM3  = $(OPTS) $(TRAPS) $(DEBUG) $(FRR)
FFLAGSF  = $(OPTSBT) $(FRR) 
FFLAGSI  = $(OPTSBT) $(FRR)
FFLAGB   = $(OPTSBT) $(FXX)
FFLAG_SER = $(OPTS_SER) $(DOUBLEFLAG) $(TRAPS) $(DEBUG) $(FXX)
FFLAGSC   =

LDR	 = $(LDRM)
ESMFLIB  = $(ESMFLIBM)
LDFLAGS =  $(LDFLAGSM)
LIBS =  $(LIBSM)

.SUFFIXES: .o .f .F .h


.SUFFIXES: .o .f .F .h
#
# *****************************************************************
#
OBJS0	= NameList_ESMFMod.o          \
	  Sfc_Var_ESMFMod.o           \
	  Nst_Var_ESMFMod.o           \
	  GFS_ErrMsgMod.o             \
	  Lib_ESMFStateAddGetMod.o    \
	  GFS_InternalState_ESMFMod.o \
	  GFS_ESMFStateMod.o          \
	  GFS_GetCf_ESMFMod.o         \
	  GFS_Initialize_ESMFMod.o    \
	  GFS_Run_ESMFMod.o           \
	  GFS_Finalize_ESMFMod.o      \
	  GFS_AddParameterToStateMod.o \
	  GFS_GridComp_ESMFMod.o      \
	  Grid_ESMFCreate.o           \
	  StartTimeGet_ESMF.o         \
	  Ensemble_sub.o              \
	  GEFS_Cpl_InternalState_ESMFMod.o         \
	  GEFS_CplState_ESMFMod.o                  \
	  GEFS_four2fln_gg.o                       \
	  GEFS_uvoedz.o                            \
	  GEFS_uveodz.o                            \
	  GEFS_sumfln_slg_gg.o                     \
	  GEFS_dozeuv.o                            \
	  GEFS_dezouv.o                            \
	  GEFS_Sto_Per_Scheme_Step1.o              \
	  GEFS_Sto_Per_Scheme_Step2.o              \
	  GEFS_Cpl_Run_ESMFMod.o                   \
	  GEFS_Cpl_ESMFMod.o                       \
	  GEFS_GetParameterFromStateMod.o          \
	  GEFS_CplComp_ESMFMod.o                   \
	  GFS_ESMF.o                               \
	  Cal_Sto_Coef.o                           \
	  mygather.o                               \
	  mpi_GV.o                                 \
	  GEFS_bcst_global.o

#	  GFS_ESMFStateAddGetMod.o
#	  GFS_Standalone_ESMF_ENS.o
#	  GFS_InputFile2ImportState.o 


OBJ_MOD	= machine.o             \
	  num_parthds.o         \
	  module.o              \
	  iounitdef.o           \
	  physcons.o            \
	  funcphys.o            \
	  progtm_module.o       \
	  rascnvv2.o            \
	  mersenne_twister.o    \
	  resol_def.o           \
	  gg_def.o              \
	  vert_def.o            \
	  sig_io.o              \
	  date_def.o            \
	  layout1.o             \
	  layout_grid_tracers.o \
	  namelist_def.o        \
	  namelist_soilveg.o    \
	  coordinate_def.o      \
	  tracer_const_h-new.o  \
	  reduce_lons_grid_module.o \
	  matrix_sig_def.o      \
	  mpi_def.o             \
	  semi_lag_def.o           \
	  deldifs_def.o            \
	  deldifs_hyb_slg_def.o    \
	  bfilt_def.o              \
	  sfcio_module.o           \
	  sigio_module.o           \
	  sigio_r_module.o         \
	  gfsio_module.o           \
	  gfsio_def.o              \
	  d3d_def.o                \
          nstio_module.o           \
          module_nst_parameters.o  \
          module_nst_water_prop.o  \
          module_nst_model.o       \
	  calpreciptype.o          \
	  module_bfmicrophysics.o  \
	  kinds.o                  \
	  peuc.o                   \
	  pran.o                   \
	  prana.o                  \
	  layout_lag.o             \
	  slglat.o                 \
	  slgshr.o                 \
	  pmgrid.o                 \
	  set_longgs.o             \
	  get_cd_hyb_slg.o         \
	  sltini_all.o             \
	  set_halos_lag_sn.o       \
	  getcon_lag.o             \
	  gfidi_speed.o            \
	  plot85_a_h.o             \
	  setlats_lag.o            \
	  gfidi_tracers.o          \
	  slgscan_all_redgg.o      \
	  get_topo_grid_grad.o     \
	  mpi_tracers_a_to_b.o     \
	  gloopa_tracers.o         \
	  ndslfv_base.o            \
	  ndslfv_para.o            \
	  ndslfv_ppm_pack.o        \
	  sicdif_hyb_slg.o


#	wnew.o      \

OBJS	= \
	get_cd_hyb.o get_am_bm_hyb.o\
	get_cd_hyb_gc.o get_am_bm_hyb_gc.o\
	wnew.o      \
	setlats.o\
	getysk.o\
	getysk_gc.o\
	ifshuff.o\
	shuffle_grid.o\
	digifilt.o\
	gcycle.o\
	wrtout.o\
	grid_collect.o\
	info.o\
	four_to_grid.o\
	hpmdummy.o\
	coundummy.o\
	damp_speed.o\
	deldifs.o\
	deldifs_tracers.o\
	deldifs_hyb_ecm.o\
	get_tc_ecm.o\
	fix_fields.o\
	spect_fields.o\
	delnpe.o\
	delnpo.o \
	dezouv.o\
	impadj_sig.locl.o\
	impadj_slg.locl.o\
	impadj_hyb.locl.o\
	impadj_hyb.locl_gc.o\
	dozeuv.o\
	epslon.o\
	filtr1eo.o \
	filtr2eo.o\
	get_ls_node.o\
	get_lats_node_shuff_a_fix.o\
	get_lats_node_shuff_r_fix.o\
	getcon.o \
	gfidi_sig.o\
	gfidi_hyb.o\
	gfidi_hyb_gc.o\
	gfidi_hyb_gc_h-new.o\
	gfidi_hyb_gc_ndsl.o\
	ndslfv_ppm_advecth1.o\
	ndslfv_ppm_advecth.o\
	glats.o\
	gloopa_opt_gg.o\
	lon_lat_para.o\
	compns.o\
	gozrineo.o\
	pln2eo.o \
	rms_spect.o \
	uvoedz.o                     \
	mpi_quit.o                   \
	sicdif_sig.o                 \
	sicdif_hyb.o                 \
	sicdif_hyb_gc.o              \
	dotstep_tracers.o            \
	ensemble_wrt.o               \
	gather_times.o               \
	four2fln_gg.o                \
	sortrx.o                     \
	sumfln_slg_gg.o              \
	sumder2_slg_gg.o             \
	triseori.o                   \
	uveodz.o                     \
	ysminv.o                     \
	fftpack.o                    \
	ver_sig.o      
# 	byteswap.o                   \


#	noblas.o                     \
#	byteswap.o                   \
##	wrtout_hyb.not.sure.o

##	bartrieo.o\
OBJS_PORT	= \
	fftpack.o \
	four2grid.fftpack.o \
	noblas.o\
	funcphys_subsx.o\
# sli
#	gloopa_sl.o\
#	gfidiu_sl.vcnfil.o\
#	lagxin.o                      \
#	lagyin.o                      \
#	lagzin.o                      \
#	lagxinit.o                      \
#	lagyinit.o                      \
#	lagzinit.o                      \
#	step0000.o \
#	plot.o \
#	slgscan1.o \
#	extlon.o                      \
#	extlat.o                      \
#	lcbas.o                      \
#	lcdbas.o                      \
#	sulag.o                      \
#	dep3dg.o                      \
#	int2d.o                      \
#	int3dv.o                      \
#	xywgts.o                      \
#
OBJS_RAD	=              \
	radlw_param.o          \
	radlw_datatb.o         \
	radlw_main.o           \
	radsw_param.o          \
	radsw_datatb.o         \
	radsw_main.o           \
	radiation_astronomy.o  \
	radiation_aerosols.o   \
	radiation_gases.o      \
	radiation_clouds.o     \
	radiation_surface.o    \
	grrad.o                \
	gloopr.o

#
#	astronomy.o    \
##	funcphys_subs.o

OBJS_PHY= \
	ozinterp.o          \
	ozphys.o            \
	gloopb.o            \
	gbphys.o            \
	gbphys_adv_hyb_gc.o \
	gbphys_adv_hyb_gc_h-new.o \
	dcyc2.o             \
	dcyc2.pre.rad.o     \
	set_soilveg.o       \
	sfc_drv.o           \
	sfc_land.o          \
	progt2.o            \
	sfc_sice.o          \
	sfc_ocean.o         \
	sfc_nst.o           \
	sfc_diff.o          \
	sfc_diag.o          \
	sflx.o              \
	moninp.o            \
	moninp1.o           \
	moninq.o            \
	moninq1.o           \
	tridi2t3.o          \
	gwdps.o             \
	gwdc.o              \
	sascnv.o            \
	sascnvn.o           \
	cnvc90.o            \
	shalcv.o            \
	shalcv_opr.o        \
	shalcnv.o           \
	lrgsclr.o           \
	gscond.o            \
	precpd.o            \
	mstadb.o            \
	mstadbtn.o          \
	mstcnv.o            \
	get_prs.o           \
	gsmddrive.o         \
	omegtes.o           \
	omegtes_gc.o        \
	omegas.o            \
	hyb2sig.o           \
	hyb2press.o         \
	hyb2press_gc.o      \
	sig2press.o 

###funcphys_subsx.o   only srt gpxs was called in fix_fields -table not used

OBJS_IO= \
	read_fix.o            \
	gribit.o              \
	wrt3d.o               \
	wrt3d_hyb.o           \
	wrtg3d.o              \
	wrtg3d_hyb.o          \
	wrtsfc.o              \
	para_fixio_w.o        \
	para_nstio_w.o        \
	sfcsub.o              \
	treadeo.io.o          \
	treadeo.gfsio.o       \
	grid_to_spec.o        \
	spect_to_grid.o       \
	spect_tv_enthalpy_ps.o\
	setsig.o              \
	twriteeo.o            \
	bafrio.o              \
	spect_send.o          \
	spect_write.o

#	para_fixio_all_w.o \
#	para_fixio_iop_w.o \
#	wrtsfc_all.o       \
#	wrtsfc_iop.o       \
## 	getaer.o           \
## 	conrad.o           \
## 	crhtab.o           \

OBJS_CC= \
	mpi_more.o  \
	cmp.comm.o  \
	atm.comm.o  \
	tiles.o

SRC	= $(OBJS0:.o=.f) $(OBJ_MOD:.o=.f) $(OBJS:.o=.f) $(OBJS_RAD:.o=.f) $(OBJS_PHY:.o=.f) $(OBJS_IO:.o=.f) $(OBJS_CC:.o=.f)


#
INCS = f_hpm.h mpi_inc.h function2

#
# *****************************************************************
#
all: msg $(EXEC)

$(EXEC): $(OBJ_MOD) $(OBJS_CC) $(OBJS0) $(OBJS) $(OBJS_PHY) $(OBJS_RAD) $(OBJS_IO) 
	$(LDR) $(LDFLAGS) -o $(EXEC) $(OBJ_MOD) $(OBJS_CC) $(OBJS0) $(OBJS) $(OBJS_PHY) $(OBJS_RAD) $(OBJS_IO) $(LIBS)

clean: msg
	rm -f *.o *.mod $(EXEC)

tar:
	tar -cvf tar.gfs.r4r8 $(SRC) $(INCS) $(COMS) $(OBJS_PORT:.o=.f) lonsper* res* xx* Makefile* ini.* scr.* m*real_?

.F.o:
	$(F77) $(FFLAGS) -c -d $< 
	#$(F77) $(FFLAGS) -c -d -WF,-DCLR:${RASV} $< 
.f.o:
	$(F77) $(FFLAGS) -c $< 


omegas.o:	omegas.f
		$(F77) $(FFLAGM) -c omegas.f

#
# *****************************************************************
#
$(OBJ_MOD) : machine.o


slgscan_all_redgg.o:pmgrid.o
gfidi_speed.o: resol_def.o
mpi_tracers_a_to_b.o:resol_def.o
gfidi_tracers.o:resol_def.o physcons.o
getcon_lag.o:gg_def.o
set_halos_lag_sn.o:resol_def.o
get_cd_hyb_slg.o:resol_def.o
sltini_all.o: resol_def.o
set_longgs.o:resol_def.o
gg_def.o : machine.o

get_topo_grid_grad.o : sig_io.o

tracer_const_h-new.o : physcons.o

resol_def.o : machine.o

sig_io.o : resol_def.o

module.o: machine.o

cnvc90.o:	cnvc90.f
		$(F77) $(FFLAGM) -c cnvc90.f

calpreciptype.o:	calpreciptype.f funcphys.o
		$(F77) $(FFLAGM3) -c calpreciptype.f

dcyc2.o:	dcyc2.f
		$(F77) $(FFLAGM) -c dcyc2.f

dcyc2.pre.rad.o:	dcyc2.pre.rad.f
		$(F77) $(FFLAGM) -c dcyc2.pre.rad.f

digifilt.o:	digifilt.f
		$(F77) $(FFLAGX) -c digifilt.f

funcphys_subsx.o:	funcphys_subsx.f
		$(F77) $(FFLAGM) -c funcphys_subsx.f

gbphys_adv_hyb_gc.o:	gbphys_adv_hyb_gc.f
		$(F77) $(FFLAGM) -c gbphys_adv_hyb_gc.f

gbphys_adv_hyb_gc_h-new.o:	gbphys_adv_hyb_gc_h-new.f
		$(F77) $(FFLAGM) -c gbphys_adv_hyb_gc_h-new.f

#gbphys_call.o:	gbphys_call.f
#		$(F77) $(FFLAGM) -c gbphys_call.f

gbphys.o:	gbphys.f
		$(F77) $(FFLAGM) -c gbphys.f

get_prs.o:	get_prs.f
		$(F77) $(FFLAGM) -c get_prs.f

byteswap.o:	byteswap.c
	 	$(FCC) -c $(CFLAGS) byteswap.c

filtr1eo.o:	filtr1eo.f 
		$(F77) $(FFLAG_SER)   -c filtr1eo.f 

filtr2eo.o:	filtr2eo.f 
		$(F77) $(FFLAG_SER)   -c filtr2eo.f 

gscond.o:	gscond.f
		$(F77) $(FFLAGM) -c gscond.f

gsmddrive.o:	gsmddrive.f module_bfmicrophysics.o
		$(F77) $(FFLAGM) -c gsmddrive.f

reduce_lons_grid_module.o:	reduce_lons_grid_module.f
		$(F77) $(FFLAGM) -c reduce_lons_grid_module.f

module_bfmicrophysics.o:	module_bfmicrophysics.f funcphys.o
		$(F77) $(FFLAGM) -c module_bfmicrophysics.f

gwdps.o:	gwdps.f
		$(F77) $(FFLAGM) -c gwdps.f

gwdc.o:		gwdc.f
		$(F77) $(FFLAGM) -c gwdc.f

hyb2press_gc.o:	hyb2press_gc.f
		$(F77) $(FFLAGM) -c hyb2press_gc.f

hyb2press.o:	hyb2press.f
		$(F77) $(FFLAGM) -c hyb2press.f

hyb2sig.o:	hyb2sig.f
		$(F77) $(FFLAGM) -c hyb2sig.f

lrgsclr.o:	lrgsclr.f
		$(F77) $(FFLAGM) -c lrgsclr.f

moninp.o:	moninp.f
		$(F77) $(FFLAGM) -c moninp.f

moninp1.o:	moninp1.f
		$(F77) $(FFLAGM) -c moninp1.f

moninq.o:	moninq.f
		$(F77) $(FFLAGM) -c moninq.f

moninq1.o:	moninq1.f
		$(F77) $(FFLAGM) -c moninq1.f

mstadb.o:	mstadb.f
		$(F77) $(FFLAGM) -c mstadb.f

mstcnv.o:	mstcnv.f
		$(F77) $(FFLAGM) -c mstcnv.f

omegtes_gc.o:	omegtes_gc.f
		$(F77) $(FFLAGM) -c omegtes_gc.f

omegtes.o:	omegtes.f
		$(F77) $(FFLAGM) -c omegtes.f

ozinterp.o:	ozinterp.f
		$(F77) $(FFLAGM) -c ozinterp.f

ozphys.o:	ozphys.f
		$(F77) $(FFLAGM) -c ozphys.f

precpd.o:	precpd.f
		$(F77) $(FFLAGM) -c precpd.f

rascnvv2.o:	rascnvv2.f physcons.o
		$(F77) $(FFLAGM) -c rascnvv2.f

sfc_sice.o:	sfc_sice.f
		$(F77) $(FFLAGM) -c sfc_sice.f

set_soilveg.o:	set_soilveg.f
		$(F77) $(FFLAGM) -c set_soilveg.f

namelist_soilveg.o:	namelist_soilveg.f
		$(F77) $(FFLAGM) -c namelist_soilveg.f

sfc_land.o:	sfc_land.f
		$(F77) $(FFLAGM) -c sfc_land.f

progt2.o:	progt2.f
		$(F77) $(FFLAGM) -c progt2.f

sfc_drv.o:	sfc_drv.f
		$(F77) $(FFLAGM) -c sfc_drv.f

sflx.o: 	sflx.f
		$(F77) $(FFLAGM) -c sflx.f

sfc_ocean.o:	sfc_ocean.f
		$(F77) $(FFLAGM) -c sfc_ocean.f

sfc_nst.o:	sfc_nst.f module_nst_model.o
		$(F77) $(FFLAGM) -c sfc_nst.f

sfc_diff.o:	sfc_diff.f
		$(F77) $(FFLAGM) -c sfc_diff.f

sfc_diag.o:	sfc_diag.f
		$(F77) $(FFLAGM) -c sfc_diag.f

sascnv.o:	sascnv.f
		$(F77) $(FFLAGM) -c sascnv.f

sascnvn.o:	sascnvn.f
		$(F77) $(FFLAGM) -c sascnvn.f

tridi2t3.o:	tridi2t3.f
		$(F77) $(FFLAGM) -c tridi2t3.f

shalcv.o:	shalcv.f
		$(F77) $(FFLAGM) -c shalcv.f

shalcv_opr.o:	shalcv_opr.f
		$(F77) $(FFLAGM) -c shalcv_opr.f

shalcnv.o:	shalcnv.f
		$(F77) $(FFLAGM) -c shalcnv.f

sig2press.o:	sig2press.f
		$(F77) $(FFLAGM) -c sig2press.f

# *****************************************************************

radlw_param.o:	radlw_param.f
		$(F77) $(FFLAGM) -c radlw_param.f


radlw_datatb.o:	radlw_datatb.f
		$(F77) $(FFLAGM) -c radlw_datatb.f


radlw_main.o:	radlw_main.f radlw_datatb.o
		$(F77) $(FFLAGM) -c radlw_main.f

radsw_param.o:	radsw_param.f
		$(F77) $(FFLAGM) -c radsw_param.f

radsw_datatb.o:	radsw_datatb.f radsw_param.o
		$(F77) $(FFLAGM) -c radsw_datatb.f

radsw_main.o:	radsw_main.f radsw_datatb.o
		$(F77) $(FFLAGM) -c radsw_main.f

radiation_astronomy.o:	radiation_astronomy.f
		$(F77) $(FFLAGM) -c radiation_astronomy.f

radiation_aerosols.o:	radiation_aerosols.f
		$(F77) $(FFLAGM) -c radiation_aerosols.f

radiation_gases.o:	radiation_gases.f
		$(F77) $(FFLAGM) -c radiation_gases.f

radiation_clouds.o:	radiation_clouds.f module_bfmicrophysics.o
		$(F77) $(FFLAGM) -c radiation_clouds.f

radiation_surface.o:	radiation_surface.f
		$(F77) $(FFLAGM) -c radiation_surface.f

grrad.o:	grrad.f radsw_datatb.o radlw_datatb.o radsw_main.o
		$(F77) $(FFLAGM) -c grrad.f

progtm_module.o:	progtm_module.f machine.o
		$(F77) $(FFLAGM) -c progtm_module.f

machine.o:	machine.f
		$(F77) $(FFLAGM) -c machine.f

num_parthds.o:	num_parthds.f
		$(F77) $(FFLAGM) -c num_parthds.f

kinds.o:	GEFS_Cpl_Cal_Sto_Coef.fd/kinds.f90
		$(F90) $(FFLAG90) -c GEFS_Cpl_Cal_Sto_Coef.fd/kinds.f90

peuc.o:		GEFS_Cpl_Cal_Sto_Coef.fd/peuc.f90
		$(F90) $(FFLAG90) -c GEFS_Cpl_Cal_Sto_Coef.fd/peuc.f90

pran.o:		GEFS_Cpl_Cal_Sto_Coef.fd/pran.f90
		$(F90) $(FFLAG90) -c GEFS_Cpl_Cal_Sto_Coef.fd/pran.f90

prana.o:	GEFS_Cpl_Cal_Sto_Coef.fd/prana.f90 pran.o peuc.o
		$(F90) $(FFLAG90) -c GEFS_Cpl_Cal_Sto_Coef.fd/prana.f90

#
# *****************************************************************
#
physcons.o:	physcons.f machine.o
		$(F77) $(FFLAG90) -c physcons.f

iounitdef.o:	iounitdef.f
		$(F77) $(FFLAG90) -c iounitdef.f

funcphys.o:	funcphys.f physcons.o
		$(F77) $(FFLAG90) -c funcphys.f

sfcio_module.o:	sfcio_module.f
		$(F77) $(FFLAGSF) -c sfcio_module.f

sigio_module.o:	sigio_module.f
		$(F77) $(FFLAGSI) -c sigio_module.f

gfsio_def.o:	gfsio_def.f gfsio_module.o
		$(F77) $(FFLAGSI) -c gfsio_def.f

gfsio_module.o:	gfsio_module.f
		$(F77) $(FFLAGSI) -c gfsio_module.f

sigio_r_module.o:	sigio_r_module.f sigio_module.o
		$(F77) $(FFLAGSI) -c sigio_r_module.f

bafrio.o:	bafrio.f
		$(F77) $(FFLAGB) -c bafrio.f

#
#sigdas.io.o:	sigdas.io.f
##		$(F77) $(FFLAGIO) -c sigdas.io.f

read_fix.o:	read_fix.f
		$(F77) $(FFLAGX) -c read_fix.f

softcount.o:	softcount.f
		$(F77) $(FFLAGX) -c softcount.f

gloopr.o:	gloopr.f grrad.o
		$(F77) $(FFLAG_SER) -c gloopr.f

sfcsub.o:	sfcsub.f
		$(F77) $(FFLAG_SER)   -c sfcsub.f

gcycle.o:	gcycle.f
		$(F77) $(FFLAGM) -c gcycle.f

getaer.o:	getaer.f
		$(F77) $(FFLAGIO) -c getaer.f

wrt3d.o:	wrt3d.f
		$(F77) $(FFLAGIO) -c wrt3d.f

gribit.o:	gribit.f
		$(F77) $(FFLAGIO) -c gribit.f

wrtsfc.o:	wrtsfc.f
		$(F77) $(FFLAGIO) -c wrtsfc.f

para_fixio_w.o:		para_fixio_w.f
		$(F77) $(FFLAGIO) -c para_fixio_w.f

para_nstio_w.o:		para_nstio_w.f
		$(F77) $(FFLAGIO) -c para_nstio_w.f

#para_fixio_all_w.o:		para_fixio_all_w.f
#		$(F77) $(FFLAGIO) -c para_fixio_all_w.f

#para_fixio_iop_w.o:	para_fixio_iop_w.f
#		$(F77) $(FFLAGIO) -c para_fixio_iop_w.f

#conrad.o:	conrad.f
#		$(F77) $(FFLAGIO) -c conrad.f

#crhtab.o:	crhtab.f
#		$(F77) $(FFLAGIO) -c crhtab.f

treadeo.io.o:	treadeo.io.f
		$(F77) $(FFLAGIO) -c treadeo.io.f

treadeo.gfsio.o:	treadeo.gfsio.f
		$(F77) $(FFLAGIO) -c treadeo.gfsio.f

twriteeo.o:	twriteeo.f
		$(F77) $(FFLAGIO) -c twriteeo.f

spect_send.o:	spect_send.f
		$(F77) $(FFLAGIO) -c spect_send.f

spect_write.o:	spect_write.f
		$(F77) $(FFLAGIO) -c spect_write.f

spect_to_grid.o:	spect_to_grid.f
		$(F77) $(FFLAGIO) -c spect_to_grid.f

spect_tv_enthalpy_ps.o:	spect_tv_enthalpy_ps.f
		$(F77) $(FFLAGIO) -c spect_tv_enthalpy_ps.f

grid_to_spec.o:	grid_to_spec.f
		$(F77) $(FFLAGIO) -c grid_to_spec.f

wrtout.o:	wrtout.f sig_io.o
		$(F77) $(FFLAGIO) -c wrtout.f

nstio_module.o:	 nstio_module.f
		$(F77) $(FFLAG90) -c nstio_module.f

module_nst_parameters.o:	module_nst_parameters.f physcons.o
		$(F77) $(FFLAG90) -c module_nst_parameters.f

module_nst_water_prop.o:	 module_nst_water_prop.f module_nst_parameters.o
		$(F77) $(FFLAG90) -c module_nst_water_prop.f

module_nst_model.o:		module_nst_model.f module_nst_water_prop.o
		$(F77) $(FFLAG90) -c module_nst_model.f

NameList_ESMFMod.o:	 NameList_ESMFMod.f
		$(F77) $(FFLAG90) -c NameList_ESMFMod.f

Sfc_Var_ESMFMod.o:	Sfc_Var_ESMFMod.f
		$(F77) $(FFLAG90) -c Sfc_Var_ESMFMod.f

Nst_Var_ESMFMod.o:	Nst_Var_ESMFMod.f
		$(F77) $(FFLAG90) -c Nst_Var_ESMFMod.f

Lib_ESMFStateAddGetMod.o:	Lib_ESMFStateAddGetMod.f
		$(F77) $(FFLAG90) -c Lib_ESMFStateAddGetMod.f

#GFS_ESMFStateAddGetMod.o:	 GFS_ESMFStateAddGetMod.f
#		$(F77) $(FFLAG90) -c GFS_ESMFStateAddGetMod.f

GFS_InternalState_ESMFMod.o:	 GFS_InternalState_ESMFMod.f sig_io.o rascnvv2.o
		$(F77) $(FFLAG90) -c GFS_InternalState_ESMFMod.f

GFS_ESMFStateMod.o:	 GFS_ESMFStateMod.f GFS_InternalState_ESMFMod.o
		$(F77) $(FFLAG90) -c GFS_ESMFStateMod.f

GFS_ErrMsgMod.o:	 GFS_ErrMsgMod.f
		$(F77) $(FFLAG90) -c GFS_ErrMsgMod.f

GFS_GetCf_ESMFMod.o:	 GFS_GetCf_ESMFMod.f GFS_InternalState_ESMFMod.o
		$(F77) $(FFLAG90) -c GFS_GetCf_ESMFMod.f

GFS_ESMFMod.o:	 GFS_ESMFMod.f
		$(F77) $(FFLAG90) -c GFS_ESMFMod.f

#GFS_Grid_fnl_ESMFMod.o:	 GFS_Grid_fnl_ESMFMod.f
		#$(F77) $(FFLAG90A) -c GFS_Grid_fnl_ESMFMod.f

GFS_GridComp_ESMFMod.o:	 GFS_GridComp_ESMFMod.f GFS_Initialize_ESMFMod.o GFS_Run_ESMFMod.o GFS_AddParameterToStateMod.o
		$(F77) $(FFLAG90A) -c GFS_GridComp_ESMFMod.f

GFS_Initialize_ESMFMod.o:	 GFS_Initialize_ESMFMod.f GFS_GetCf_ESMFMod.o
		$(F77) $(FFLAG90) -c GFS_Initialize_ESMFMod.f

GFS_Run_ESMFMod.o:	GFS_Run_ESMFMod.f GFS_InternalState_ESMFMod.o
		$(F77) $(FFLAG90) -c GFS_Run_ESMFMod.f

GFS_Finalize_ESMFMod.o:	GFS_Finalize_ESMFMod.f GFS_InternalState_ESMFMod.o
		$(F77) $(FFLAG90) -c GFS_Finalize_ESMFMod.f

GFS_InputFile2ImportState.o:	 GFS_InputFile2ImportState.f
		$(F77) $(FFLAG90) -c GFS_InputFile2ImportState.f

#GFS_Standalone_ESMF_ENS.o:	 GFS_Standalone_ESMF_ENS.f
#		$(F77) $(FFLAG90) -c GFS_Standalone_ESMF_ENS.f

GFS_ESMF.o:	GFS_ESMF.f GEFS_CplComp_ESMFMod.o GFS_GridComp_ESMFMod.o
		$(F77) $(FFLAG90) -c GFS_ESMF.f

Grid_ESMFCreate.o:	 Grid_ESMFCreate.f GFS_InternalState_ESMFMod.o
		$(F77) $(FFLAG90) -c Grid_ESMFCreate.f

StartTimeGet_ESMF.o:	 StartTimeGet_ESMF.f sigio_r_module.o
		$(F77) $(FFLAG90) -c StartTimeGet_ESMF.f

Ensemble_sub.o:	 Ensemble_sub.f
		$(F77) $(FFLAG90) -c Ensemble_sub.f

mpi_more.o:	mpi_more.f
		$(F77) $(FFLAGS) -c mpi_more.f

cmp.comm.o:	cmp.comm.f
		$(F77) $(FFLAGS) -c cmp.comm.f

atm.comm.o:	atm.comm.f cmp.comm.o
		$(F77) $(FFLAGS) -c atm.comm.f

tiles.o:	tiles.f atm.comm.o
		$(F77) $(FFLAGS) -c tiles.f

GEFS_Cpl_InternalState_ESMFMod.o:	GEFS_Cpl_InternalState_ESMFMod.f
		$(F77) $(FFLAG90) -c GEFS_Cpl_InternalState_ESMFMod.f

GEFS_CplState_ESMFMod.o:	GEFS_CplState_ESMFMod.f GEFS_Cpl_InternalState_ESMFMod.o atm.comm.o
		$(F77) $(FFLAG90) -c GEFS_CplState_ESMFMod.f

GEFS_Sto_Per_Scheme.o:	GEFS_Sto_Per_Scheme.f GEFS_Cpl_InterfaceState_ESMFMod.o
		$(F77) $(FFLAG90) -c GEFS_Sto_Per_Scheme.f

GEFS_Cpl_Run_ESMFMod.o:	GEFS_Cpl_Run_ESMFMod.f atm.comm.o
		$(F77) $(FFLAG90) -c GEFS_Cpl_Run_ESMFMod.f

GEFS_Cpl_ESMFMod.o:	GEFS_Cpl_ESMFMod.f GEFS_Cpl_Run_ESMFMod.o
		$(F77) $(FFLAG90) -c GEFS_Cpl_ESMFMod.f

GEFS_CplComp_ESMFMod.o:	GEFS_CplComp_ESMFMod.f GEFS_Cpl_ESMFMod.o
		$(F77) $(FFLAG90) -c GEFS_CplComp_ESMFMod.f

GFS_AddParameterToStateMod.o:	GFS_AddParameterToStateMod.f GFS_InternalState_ESMFMod.o
		$(F77) $(FFLAG90) -c GFS_AddParameterToStateMod.f

GEFS_Sto_Per_Scheme_Step1.o:	GEFS_Sto_Per_Scheme_Step1.f
		$(F77) $(FFLAG90) -c GEFS_Sto_Per_Scheme_Step1.f

GEFS_Sto_Per_Scheme_Step2.o:	GEFS_Sto_Per_Scheme_Step2.f
		$(F77) $(FFLAG90) -c GEFS_Sto_Per_Scheme_Step2.f

GEFS_GetParameterFromStateMod.o:	GEFS_GetParameterFromStateMod.f
		$(F77) $(FFLAG90) -c GEFS_GetParameterFromStateMod.f

Cal_Sto_Coef.o:		GEFS_Cpl_Cal_Sto_Coef.fd/Cal_Sto_Coef.f90 prana.o
		$(F90) $(FFLAG90) -c GEFS_Cpl_Cal_Sto_Coef.fd/Cal_Sto_Coef.f90

GEFS_bcst_global.o:	GEFS_bcst_global.f
		$(F77) $(FFLAG90) -c GEFS_bcst_global.f

mpi_GV.o:		mpi_GV.f
		$(F77) $(FFLAGS) -c mpi_GV.f

mygather.o:		mygather.c
		$(CC) $(FFLAGSC) -c mygather.c
msg:
	@echo ""
	@echo "============================================================"
	@echo " CFS (GFS)"
	@echo "============================================================"
