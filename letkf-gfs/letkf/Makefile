include ../../config/makefile.mk

PROGS = letkf obsope mean efso

OBJS    = letkf_obs.o letkf_tools.o efso_tools.o efso_nml.o

COMMON_DIR = ../common
COMMON_OBJS =   $(COMMON_DIR)/common.o         \
                $(COMMON_DIR)/common_mpi.o     \
                $(COMMON_DIR)/common_mtx.o     \
                $(COMMON_DIR)/common_letkf.o   \
                $(COMMON_DIR)/common_gfs.o     \
                $(COMMON_DIR)/common_obs_gfs.o \
                $(COMMON_DIR)/common_mpi_gfs.o \
		$(COMMON_DIR)/common_precip.o \
	  	$(COMMON_DIR)/SFMT.o           \
                $(COMMON_DIR)/netlib.o         

INCS   = -I../common -I$(NCEPLIBDIR)/incmod/sigio_4

LIBS   = -L$(NCEPLIBDIR)/lib -lsigio_4 -lsp_4 $(LIBBLAS)

FOPTS = -O2 $(F_CNVTBE)

all:	common $(PROGS)

common:
	make -C ../common

letkf.o : letkf_obs.o letkf_tools.o
letkf_obs.o : 
letkf_tools.o : letkf_obs.o efso_nml.o efso_tools.o
efso_nml.o : 
efso_tools.o : letkf_obs.o efso_nml.o 
mean.o : letkf_obs.o
efso.o : efso_nml.o efso_tools.o letkf_obs.o letkf_tools.o

%.o : %.f90
	$(FC) $(FOPTS) $(FMCMODEL) $(INCS) -c $<

% : %.o
	$(FC) $(FOPTS) $(FMCMODEL) -o $@ $(OBJS) $< $(COMMON_OBJS) $(LIBS) 

clean:
	rm -f *.o *.mod $(PROGS) letkf[0-9]* mean[0-9]* efso[0-9]*
