SOURCES_D = src/largeQvalue.d src/parse_arg.d

largeQvalue : $(SOURCES_D) src/bootstrap.o src/libspline.a
	gdc -frelease -finline-functions -O3 -Werror -Wall -fversion=Have_largeQvalue $(SOURCES_D) src/bootstrap.o src/libspline.a -lblas -lgsl -lgslcblas -lm -o largeQvalue
	rm -f src/bootstrap.o largeQvalue.o
	strip largeQvalue

src/libspline.a : src/spline_src/* header/*
	gcc -Iheader/ -c src/spline_src/*.c src/spline_src/*.f
	ar rcs src/libspline.a bsplvd.o bvalue.o bvalus.o dpbfa.o dpbsl.o interv.o sgram.o sinerp.o spline_fit.o sslvrg.o stxwx.o
	rm -f bsplvd.o bvalue.o bvalus.o dpbfa.o dpbsl.o interv.o sgram.o sinerp.o spline_fit.o sslvrg.o stxwx.o

.PHONY : sample clean boot test

clean :
	rm -f *.o largeQvalue

sample :
	largeQvalue --header --col 4 --out temp --param parameter_file data/vQTLresults.txt

boot :
	largeQvalue --boot --header --col 4 --out temp --param parameter_file data/vQTLresults.txt

test :
	largeQvalue --header --col 4 --out temp --param parameter_file data/vQTLresults.txt
	diff temp data/sample_test | head
	diff parameter_file data/sample_param | head
	largeQvalue --boot --header --col 4 --out temp --param parameter_file data/vQTLresults.txt
	diff temp data/boot_test | head
	diff parameter_file data/boot_param | head
	largeQvalue --fast 0.05 --out temp data/nominal
	diff data/nominal_results temp | head
	rm temp parameter_file

