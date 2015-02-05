largeqvalue : src/spline.c  src/largeQvalue.d src/parse_arg.d
	gcc -c src/spline.c -o spline.o
	ldc2  src/largeQvalue.d  src/parse_arg.d spline.o -L-lgsl -L-lgslcblas -O3 -of="largeqvalue"

dmd : src/spline.c  src/largeQvalue.d src/parse_arg.d
	gcc -c  src/spline.c -o spline.o
	dmd -release -noboundscheck -inline -O  src/largeQvalue.d  src/parse_arg.d spline.o -L-lgsl -L-lgslcblas

.PHONY : sample clean boot

clean :
	rm -f *.o largeqvalue

sample :
	largeqvalue --header --col 4 --out temp --param parameter_file data/vQTLresults.txt

boot :
	largeqvalue --boot --header --col 4 --out temp --param parameter_file data/vQTLresults.txt

test :
	largeqvalue --header --col 4 --out temp --param parameter_file data/vQTLresults.txt
	diff temp data/sample_test | head
	diff parameter_file data/sample_param | head
	largeqvalue --boot --header --col 4 --out temp --param parameter_file data/vQTLresults.txt
	diff temp data/boot_test | head
	diff parameter_file data/boot_param | head

