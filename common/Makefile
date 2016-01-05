include $(CFS_LETKF_ROOT)/config/makefile.mk

libs=common_cfs.a
objects=common.o common_obs.o SFMT.o kdtree.o
F90OPT=-assume byterecl -O2 -traceback -g

all: $(libs)

common_cfs.a : $(objects)
	ar -rcs common_cfs.a $(objects)

%: %.o
	$(FC) $(F90OPT) -o $@ $^

%.o : %.f90
	$(FC) $(F90OPT) -c $<

%.o : %.f
	$(FC) $(F90OPT) -c $<


.PHONY: all

.PHONY: clean

clean:
	rm -f *.o *.mod *~

veryclean: clean
	rm -f $(PROGRAMS)
