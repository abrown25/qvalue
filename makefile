SOURCES_C = src/bootstrap.c src/spline_fit.c
SOURCES_D = src/largeQvalue.d src/parse_arg.d
SOURCES_F = src/bsplvd.f src/bvalue.f src/bvalus.f src/sgram.f src/sinerp.f src/sslvrg.f src/stxwx.f
OBJECTS = $(SOURCES_C:.c=.o) $(SOURCES_F:.f=.o)

largeQvalue : $(SOURCES_C) $(SOURCES_D) $(SOURCES_F)
	gcc -I /usr/include/R -c $(SOURCES_C) $(SOURCES_F)
	mv *o src/
	gdc $(SOURCES_D) $(OBJECTS) -L/usr/lib/R/lib/ -lR -lgsl -lgslcblas -lm -o largeQvalue
	rm src/*o

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

