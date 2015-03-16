largeQvalue : src/largeQvalue.d src/parse_arg.d src/bootstrap.c src/spline_fit.c src/bsplvd.f src/bvalue.f src/bvalus.f src/sgram.f src/sinerp.f src/sslvrg.f src/stxwx.f
	gcc -I /usr/include/R -c src/bootstrap.c src/spline_fit.c src/bsplvd.f src/bvalue.f src/bvalus.f src/sgram.f src/sinerp.f src/sslvrg.f src/stxwx.f
	gdc src/largeQvalue.d src/parse_arg.d bootstrap.o bsplvd.o bvalue.o bvalus.o sgram.o sinerp.o spline_fit.o sslvrg.o stxwx.o -L /usr/lib/R/lib/ -lR -lgsl -lgslcblas -lm -o largeQvalue
	rm *o

dmd : src/spline.c  src/largeQvalue.d src/parse_arg.d
	gcc -c  src/spline.c -o spline.o
	dmd -release -noboundscheck -inline -O  src/largeQvalue.d  src/parse_arg.d spline.o -L-lgsl -L-lgslcblas

.PHONY : sample clean boot

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

