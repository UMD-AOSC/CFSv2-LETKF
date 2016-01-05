include $(CFS_LETKF_ROOT)/config/makefile.mk

# COMMON_DIR = $(CFS_LETKF_ROOT)/letkf-gfs/common

# COMMON_OBJS =   $(COMMON_DIR)/SFMT.o           \
# 		$(COMMON_DIR)/common.o         \
#                 $(COMMON_DIR)/common_mtx.o     \
#                 $(COMMON_DIR)/netlib.o         \
#                 $(COMMON_DIR)/common_letkf.o   \
#                 $(COMMON_DIR)/common_gfs.o     \
#                 $(COMMON_DIR)/common_obs_gfs.o

common_cfs := $(CFS_LETKF_ROOT)/common/common_cfs.a
common_gfs := $(CFS_LETKF_ROOT)/letkf-gfs/common/common_gfs.a

INCS   = -I$(NCEPLIBDIR)/incmod/sigio_4 -I$(CFS_LETKF_ROOT)/common -I$(CFS_LETKF_ROOT)/letkf-gfs/common
LIBS := -L$(NCEPLIBDIR)/lib   $(common_gfs) $(common_cfs) -lsigio_4 -lsp_4 $(LIBBLAS)
FOPTS= -O2 -traceback -g

.PHONY: all clean

all: superob

superob: superob.o superob_tools.o obs_tools.o obs_gfs_ext.o
	$(FC) $(FOPTS) $(COMMON_OBJS) $^ -o $@ $(LIBS)

superob.o : superob_tools.o
superob_tools.o : obs_tools.o obs_gfs_ext.o

%.o:%.f90
	$(F90) $(FOPTS) $(INCS) -c $<

clean:
	rm -f *.o *.mod superob
