
.PHONY: all clean


all: grabbufr dec_prepbufr superob

grabbufr:
	make -C src/grabbufr
	cp src/grabbufr/grabbufr .

dec_prepbufr:
	make -C src/dec_prepbufr
	cp src/dec_prepbufr/dec_prepbufr .
superob:
	make -C src/superob
	cp src/superob/superob .

clean:
	make -C src/grabbufr clean
	make -C src/dec_prepbufr clean
	make -C src/superob clean
	rm grabbufr dec_prepbufr superob
