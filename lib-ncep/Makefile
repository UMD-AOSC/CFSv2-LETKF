include $(CFS_LETKF_ROOT)/config/makefile.mk

TARGETS := ip sp bacio sigio w3lib-2.0 gfsio sfcio landsfcutil bufrlib

.PHONY: all clean dirs $(TARGETS)

all: msg dirs $(TARGETS)

dirs :
	@mkdir -p lib
	@mkdir -p incmod

landsfcutil: w3lib-2.0
	+make -C src/landsfcutil
	cp src/landsfcutil/lib_?/*.a lib/
	mkdir -p incmod/landsfcutil_4
	mkdir -p incmod/landsfcutil_8
	mkdir -p incmod/landsfcutil_d
	cp src/landsfcutil/lib_4/*.mod incmod/landsfcutil_4/
	cp src/landsfcutil/lib_8/*.mod incmod/landsfcutil_8/
	cp src/landsfcutil/lib_d/*.mod incmod/landsfcutil_d/

bufrlib:
	+make -C src/bufr_10.2.3
	cp src/bufr_10.2.3/*.a lib/
# mkdir -p incmod/bufrlib
# cp src/bufrlib/*.mod incmod/bufrlib

sfcio:
	+make -C src/sfcio
	cp src/sfcio/*.a lib/
	mkdir -p incmod/sfcio_4
	cp src/sfcio/*.mod incmod/sfcio_4

gfsio:
	+make -C src/gfsio
	cp src/gfsio/*.a lib/
	mkdir -p incmod/gfsio_4
	cp src/gfsio/*.mod incmod/gfsio_4

sigio:
	+make -C src/sigio
	cp src/sigio/*.a lib/
	mkdir -p incmod/sigio_4
	cp src/sigio/*.mod incmod/sigio_4

ip:
	+make -C src/ip
	cp src/ip/libip_?/*.a lib/

sp:
	+make -C src/sp
	cp src/sp/libsp_?/*.a lib/

bacio:
	+make -C src/bacio
	cp src/bacio/libbacio_?/*.a lib/

w3lib-2.0: sigio
	+make -C src/w3lib-2.0
	cp src/w3lib-2.0/lib_?/*.a lib/

clean: msg
	rm -rf lib
	rm -rf incmod
	+make -C src/bacio clean
	+make -C src/bufr_10.2.3 clean
	+make -C src/gfsio clean
	+make -C src/ip clean
	+make -C src/landsfcutil clean
	+make -C src/sfcio clean
	+make -C src/sigio clean
	+make -C src/sp clean
	+make -C src/w3lib-2.0 clean

msg:
	@echo ""
	@echo "============================================================"
	@echo " NCEP libraries"
	@echo "============================================================"
