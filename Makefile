.PHONY: all clean
all:
	+make -C lib-ncep
	+make -C cfs
	+make -C util
	+make -C letkf-gfs
	+make -C letkf-mom

clean:
	+make -C lib-ncep clean
	+make -C cfs clean
	+make -C util clean
	+make -C letkf-gfs clean
	+make -C letkf-mom clean
