SOURCES_D = src/largeQvalue.d src/parse_arg.d
LIB_GSL = /usr/lib/libblas/libblas.a /usr/lib/libgsl.a /usr/lib/libgslcblas.a

bin/largeQvalue : $(SOURCES_D) src/bootstrap.o src/libspline.a
	gdc -frelease -finline-functions -O3 -Werror -Wall -fversion=Have_largeQvalue $(SOURCES_D) src/bootstrap.o src/libspline.a -lblas -lgsl -lgslcblas -lm -o bin/largeQvalue
	rm -f src/bootstrap.o largeQvalue.o
	strip bin/largeQvalue

src/libspline.a : src/spline_src/* header/*
	gcc -Iheader/ -c src/spline_src/*.c src/spline_src/*.f
	ar rcs src/libspline.a bsplvd.o bvalue.o bvalus.o dpbfa.o dpbsl.o interv.o sgram.o sinerp.o spline_fit.o sslvrg.o stxwx.o
	rm -f bsplvd.o bvalue.o bvalus.o dpbfa.o dpbsl.o interv.o sgram.o sinerp.o spline_fit.o sslvrg.o stxwx.o

.PHONY : sample clean boot test static install dmd

static	: $(SOURCES_D) $(LIB_GSL) src/bootstrap.o src/libspline.a
	gdc -frelease -finline-functions -O3 -Werror -Wall -fversion=Have_largeQvalue $(SOURCES_D) src/bootstrap.o src/libspline.a $(LIB_GSL) -o bin/largeQvalue
	rm -f src/bootstrap.o largeQvalue.o
	strip bin/largeQvalue

dmd	: $(SOURCES_D) $(LIB_GSL) src/bootstrap.o src/libspline.a
	dmd -release -inline -O -w -version=Have_largeqvalue $(SOURCES_D) src/bootstrap.o src/libspline.a $(LIB_GSL) -ofbin/largeQvalue
	rm -f src/bootstrap.o largeQvalue.o
	strip bin/largeQvalue

install : bin/largeQvalue largeQvalue.1
	ln -s $(shell pwd)/bin/largeQvalue /usr/local/bin/
	ln -s $(shell pwd)/largeQvalue.1 /usr/local/man/man1/

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
	bin/largeQvalue --fast 0.05 --out temp data/nominal
	diff data/nominal_results temp | head
	rm temp parameter_file

