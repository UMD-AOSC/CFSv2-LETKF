.PHONY: all  clean

all: msg netlib.o
#	+make -f Makefile.common

clean: msg
	rm -f *.mod *.o

netlib.o:netlib.f
	$(F90) -O2 -c netlib.f

msg:
	@echo ""
	@echo "============================================================"
	@echo " GFS-LETKF (common)"
	@echo "============================================================"
