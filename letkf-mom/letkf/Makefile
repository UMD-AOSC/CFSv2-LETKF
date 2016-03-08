.PHONY: all clean

all: common msg
	+make -f Makefile.letkf

clean: msg
	rm -f letkf 
	+make -f Makefile.letkf clean

common:
	+make -C ../common

msg:
	@echo ""
	@echo "============================================================"
	@echo " MOM-LETKF (letkf)"
	@echo "============================================================"
