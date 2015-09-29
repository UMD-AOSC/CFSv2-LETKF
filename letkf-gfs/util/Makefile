
.PHONY: all clean


all: grabbufr dec_prepbufr

grabbufr:
	make -C src/grabbufr
	cp src/grabbufr/grabbufr .

dec_prepbufr:
	make -C src/dec_prepbufr
	cp src/dec_prepbufr/dec_prepbufr .


clean:
	make -C src/grabbufr clean
	make -C src/dec_prepbufr clean
	rm grabbufr dec_prepbufr
