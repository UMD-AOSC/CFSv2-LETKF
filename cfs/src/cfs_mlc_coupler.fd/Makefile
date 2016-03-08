FFLAGS  = -O3 $(F_R8) -g -traceback $(F_FREE) -assume byterecl
F77	= $(FC)

PROG=coupler

cfs_mlc_$(PROG):msg	coupler_module.o $(PROG).o
	$(F77) $(FFLAGS)  coupler_module.o $(PROG).o -o cfs_mlc_$(PROG)


coupler_module.o:	coupler_module.f
	$(F77) $(FFLAGS) -c coupler_module.f
$(PROG).o:	$(PROG).f
	$(F77) $(FFLAGS) -c $(PROG).f

clean :msg
	rm -f cfs_mlc_coupler
	rm -f *.o
	rm -f *.mod
msg:
	@echo ""
	@echo "============================================================"
	@echo " CFS (coupler)"
	@echo "============================================================"    
