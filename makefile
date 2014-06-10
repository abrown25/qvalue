large_q_value : spline.c large_q_value.d parse_arg.d
	gcc -c spline.c -o spline.o
	dmd large_q_value.d parse_arg.d spline.o -L-lgsl -L-lgslcblas

ldc :  spline.c large_q_value.d parse_arg.d
	gcc -c spline.c -o spline.o
	ldc2 large_q_value.d parse_arg.d spline.o -L-lgsl -L-lgslcblas -O3 -of="large_q_value_ldc"

.PHONY : sample clean

clean :
	rm *.o large_q_value large_q_value_ldc

sample :
	./large_q_value --header --col 4 --out temp --param parameter_file vQTLresults.txt
