GSL = /usr/lib/libgsl.a /usr/lib/libgslcblas.a

large_q_value : spline.c large_q_value.d parse_arg.d
	gcc -c spline.c -o spline.o
	dmd -release -noboundscheck -inline -O large_q_value.d parse_arg.d spline.o ${GSL}

ldc :  spline.c large_q_value.d parse_arg.d
	gcc -c spline.c -o spline.o
	ldc2 large_q_value.d parse_arg.d spline.o ${GSL} -O3 -of="large_q_value_ldc"

.PHONY : sample clean boot

clean :
	rm -f *.o large_q_value large_q_value_ldc

sample :
	./large_q_value --header --col 4 --out temp --param parameter_file vQTLresults.txt

boot :
	./large_q_value --boot --header --col 4 --out temp --param parameter_file vQTLresults.txt
