
PROGS   = ss2grd grd2ss sscycle ss2grdp grdctl sscycle2

OBJS    = ssio_tools.o physcons.o

COMMON_DIR = ./

INCS = -I$(NCEPLIBDIR)/incmod/sigio_4 -I$(NCEPLIBDIR)/incmod/sfcio_4
LIBS = -L$(NCEPLIBDIR)/lib -lsigio_4 -lbacio_4 -lsfcio_4 -lsp_4

FOPTS  = -O2 -traceback -g
FOPTS2 = -O2

.PHONY: all clean

all: msg $(PROGS)

commonlibs= common_gfs_pres.o common_gfs.o ssio_tools.o
sscycle2.o :
ss2grd.o : common_gfs.o ssio_tools.o
grd2ss.o : common_gfs.o ssio_tools.o
ss2grdp.o : common_gfs.o ssio_tools.o
grdctl.o : common_gfs.o ssio_tools.o
sscycle.o : common.o common_gfs.o ssio_tools.o

common_gfs.o : common.o
ssio_tools.o: common_gfs_pres.o physcons.o
physcons.o : common.o
common_gfs_pres.o : common.o common_gfs.o



% : %.o $(commonlibs)
	$(FC) $(FOPTS) $^ $(LIBS) -o $@

%.o : %.f90
	$(FC) $(FOPTS) -c $<  -I$(COMMON_DIR) $(INCS) -o $@

msg:
	@echo ""
	@echo "============================================================"
	@echo " ssio"
	@echo "============================================================"

clean: msg
	rm -f *.o *.mod $(PROGS)
