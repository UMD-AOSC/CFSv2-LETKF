.PHONY: all clean
all: netlib.o
	+make -f Makefile.letkf
	+make -f Makefile.mean
	+make -f Makefile.obsope
	+make -f Makefile.efso

clean:
	rm -f letkf mean obsope efso
	make -f Makefile.letkf clean
	make -f Makefile.mean clean
	make -f Makefile.obsope clean
	make -f Makefile.efso clean

netlib.o:netlib.f
	mpiifort -c ./netlib.f
