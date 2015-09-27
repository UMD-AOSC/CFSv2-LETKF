include $(CFS_LETKF_ROOT)/config/makefile.mk

SUBDIRS = ssio

.PHONY : all dirs $(SUBDIRS)

all : dirs $(SUBDIRS)

dirs :
	@mkdir -p bin

ssio : dirs
	+make -C src/ssio
	cp src/ssio/grd2ss bin/
	cp src/ssio/grdctl bin/
	cp src/ssio/ss2grd bin/
	cp src/ssio/ss2grdp bin/
	cp src/ssio/sscycle bin/

clean :
	rm -rf bin
	+make -C src/ssio clean
