include $(CFS_LETKF_ROOT)/config/makefile.mk

.PHONY : all

all: grabbufr

grabbufr: grabbufr.o spbufr.o
	$(FC) -o $@ $?

%.o : %.f
	$(FC) -c $< -o $@

clean:
	rm *.o grabbufr
