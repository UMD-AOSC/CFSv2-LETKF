include $(CFS_LETKF_ROOT)/config/makefile.mk

SUBDIRS = gfs mom coupler chgres cycle sfchdr sighdr

.PHONY : all dirs $(SUBDIRS)

all : dirs $(SUBDIRS)

dirs :
	@mkdir -p bin



gfs : FCC       := $(CC)
gfs : LDRM      := $(F90)
gfs : CFLAGS    := -DLINUX
gfs : FRRM      := $(F_FREE)
gfs : FINC      := -I$(ESMFMODS)
gfs : LIBSM     := -L$(NCEPLIBDIR)/lib -L$(ESMFLIBS) -mkl -lbacio_4 -lsp_d -lw3lib-2.0_d -lrt -lesmf -lstdc++ -openmp
gfs : OPTSB     := -g -O3 $(F_CNVTBE) $(F_STRICT) -heap-arrays -xHost 
gfs : OPTSBT    := $(OPTSB) -traceback
gfs : OPTSM     := $(OPTSBT) $(F_R8) -openmp
gfs : OPTSIOM   := $(OPTSBT) $(F_R8)
gfs : OPTS_SERM := $(OPTSBT) $(F_R8)
gfs : OPTS90M   := $(OPTSBT) $(F_R8) $(FRRM)
gfs : OPTS90AM  := $(OPTSBT) $(F_R8) $(FRRM)
gfs : dirs
	+make -C src/cfs_global_fcst.fd
	cp src/cfs_global_fcst.fd/global_fcst bin/



mom : dirs
	+make -C src/cfs_ocean_mom4ice.fd
	cp src/cfs_ocean_mom4ice.fd/compile/exec/cfs_ocean_mom4ice bin/




coupler : dirs
	+make -C src/cfs_mlc_coupler.fd
	cp src/cfs_mlc_coupler.fd/cfs_mlc_coupler bin/


chgres : dirs
	+make -C src/global_chgres.fd
	cp src/global_chgres.fd/global_chgres bin/

cycle : dirs
	+make -C src/global_cycle.fd
	cp src/global_cycle.fd/global_cycle bin/

sfchdr : dirs
	+make -C src/global_sfchdr.fd
	cp src/global_sfchdr.fd/global_sfchdr bin/

sighdr : dirs
	+make -C src/global_sighdr.fd
	cp src/global_sighdr.fd/global_sighdr bin/

clean :
	rm -rf bin
	+make -C src/cfs_global_fcst.fd clean
	+make -C src/cfs_mlc_coupler.fd clean
	+make -C src/cfs_ocean_mom4ice.fd clean
	+make -C src/global_chgres.fd clean
	+make -C src/global_cycle.fd clean
	+make -C src/global_sfchdr.fd clean
	+make -C src/global_sighdr.fd clean
