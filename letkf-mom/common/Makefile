.PHONY: all clean

all: msg netlib.o

clean: msg
	rm -f *.mod *.o

netlib.o:netlib.f
	$(F90) -O2 -c netlib.f

msg:
	@echo ""
	@echo "============================================================"
	@echo " MOM-LETKF (common)"
	@echo "============================================================"
