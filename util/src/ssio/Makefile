.PHONY: all clean
all: msg
#	+make -C ../common

	+make -f Makefile.grd2ss
	+make -f Makefile.grdctl
	+make -f Makefile.ss2grd
	+make -f Makefile.ss2grdp
	+make -f Makefile.sscycle
	+make -f Makefile.sscycle2


clean: msg
	rm -f ss2grd grd2ss sscycle ss2grdp grdctl sscycle2 *.o *.mod

msg:
	@echo ""
	@echo "============================================================"
	@echo " ssio"
	@echo "============================================================"	
