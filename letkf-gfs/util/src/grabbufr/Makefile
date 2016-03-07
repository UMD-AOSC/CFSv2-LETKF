#include $(CFS_LETKF_ROOT)/config/makefile.mk

.PHONY : all


all: msg grabbufr

msg:
	@echo ""
	@echo "============================================================"
	@echo " grabbufr"
	@echo "============================================================"

grabbufr: grabbufr.o spbufr.o
	$(FC) -o $@ $?

%.o : %.f
	$(FC) -c $< -o $@

clean: msg
	rm *.o grabbufr -f
