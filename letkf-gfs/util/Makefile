
.PHONY: all clean


all:
	+make -C src/grabbufr
	cp src/grabbufr/grabbufr .

	+make -C src/dec_prepbufr
	cp src/dec_prepbufr/dec_prepbufr .

#	+make -C src/superob
#	cp src/superob/superob .

clean:
	rm grabbufr dec_prepbufr -f #superob -f
	+make -C src/grabbufr clean
	+make -C src/dec_prepbufr clean
#	+make -C src/superob clean
