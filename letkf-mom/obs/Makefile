.PHONY: all clean
all: common msg gsw_oceanographic_toolbox.o
	+make -f Makefile.obsop

clean: msg
	rm -f  obsop
	make -f Makefile.obsop clean


gsw_oceanographic_toolbox.o:gsw_oceanographic_toolbox.f90
	$(F90) -c ./gsw_oceanographic_toolbox.f90

common:
	+make -C ../common

msg:
	@echo ""
	@echo "============================================================"
	@echo " MOM-LETKF (obsop)"
	@echo "============================================================"
