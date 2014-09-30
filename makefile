GSL = /usr/lib/libgsl.a /usr/lib/libgslcblas.a

largeQvalue : spline.c largeQvalue.d parse_arg.d
	gcc -c spline.c -o spline.o
	ldc2 largeQvalue.d parse_arg.d spline.o ${GSL} -O3 -of="largeQvalue"

dmd : spline.c largeQvalue.d parse_arg.d
	gcc -c spline.c -o spline.o
	dmd -release -noboundscheck -inline -O largeQvalue.d parse_arg.d spline.o ${GSL}

.PHONY : sample clean boot

clean :
	rm -f *.o largeQvalue

sample :
	./largeQvalue --header --col 4 --out temp --param parameter_file vQTLresults.txt

boot :
	./largeQvalue --boot --header --col 4 --out temp --param parameter_file vQTLresults.txt
