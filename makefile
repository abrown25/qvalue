SOURCES_D = src/main.d src/all_p_vals.d src/parse_arg.d src/pi0_calc.d src/threshold.d
LIB_GSL = /usr/lib/libblas.a /usr/lib/libgsl.a /usr/lib/libgslcblas.a

bin/largeQvalue : $(SOURCES_D) src/bootstrap.o src/libspline.a
	ldc -ofbin/largeQvalue -release -enable-inlining -O -w -oq -d-version=Have_largeqvalue -Isrc/ src/bootstrap.o $(SOURCES_D) src/libspline.a -L=-lgsl -L=-lgslcblas -L=-lm -L=-lblas
	rm -f src/bootstrap.o largeQvalue.o
	strip bin/largeQvalue

src/libspline.a : src/spline_src/* header/*
	gcc -Iheader/ -c src/spline_src/*.c src/spline_src/*.f
	ar rcs src/libspline.a bsplvd.o bvalue.o bvalus.o dpbfa.o dpbsl.o interv.o sgram.o sinerp.o spline_fit.o sslvrg.o stxwx.o
	rm -f bsplvd.o bvalue.o bvalus.o dpbfa.o dpbsl.o interv.o sgram.o sinerp.o spline_fit.o sslvrg.o stxwx.o

.PHONY : ldc sample clean boot test static install dmd

static	: $(SOURCES_D) $(LIB_GSL) src/bootstrap.o src/libspline.a
	ldc -ofbin/largeQvalue -release -enable-inlining -O -w -oq -d-version=Have_largeqvalue -Isrc/ $(SOURCES_D) src/bootstrap.o src/libspline.a $(LIB_GSL)
	rm -f src/bootstrap.o largeQvalue.o
	strip bin/largeQvalue

gdc : $(SOURCES_D) src/bootstrap.o src/libspline.a
	gdc -frelease -finline-functions -O3 -Werror -Wall -fversion=Have_largeQvalue $(SOURCES_D) src/bootstrap.o src/libspline.a -lblas -lgsl -lgslcblas -lcurl -o bin/largeQvalue
	rm -f src/bootstrap.o largeQvalue.o
	strip bin/largeQvalue

dmd	: $(SOURCES_D) src/bootstrap.o src/libspline.a
	dmd -release -inline -O -w -version=Have_largeqvalue $(SOURCES_D) src/bootstrap.o src/libspline.a -L-lgsl -L-lgslcblas -L-lblas -ofbin/largeQvalue
	rm -f src/bootstrap.o largeQvalue.o
	strip bin/largeQvalue

install : bin/largeQvalue largeQvalue.1
	cp -v $(shell pwd)/bin/largeQvalue /usr/local/bin/
	cp -v $(shell pwd)/largeQvalue.1 /usr/local/man/man1/

clean :
	rm -f *.o bin/largeQvalue

sample :
	bin/largeQvalue --header --col 4 --out temp --param parameter_file data/vQTLresults.txt

boot :
	bin/largeQvalue --boot --header --col 4 --out temp --param parameter_file data/vQTLresults.txt

test :
	bin/largeQvalue --header --col 4 --out temp --param parameter_file data/vQTLresults.txt
	diff temp data/sample_test | head
	diff parameter_file data/sample_param | head
	bin/largeQvalue --boot --header --col 4 --out temp --param parameter_file data/vQTLresults.txt
	diff temp data/boot_test | head
	diff parameter_file data/boot_param | head
	bin/largeQvalue --fast 0.05 --out temp --col 10 data/nominal
	diff data/nominal_results temp | head
	rm temp parameter_file
