spline_fit	:	sbart.c bsplvd.f bvalue.f sgram.f sslvrg.f stxwx.f
	gcc -I /usr/include/R/ -c spline_fit.c
	gfortran -c bsplvd.f bvalue.f bvalus.f sgram.f sinerp.f sslvrg.f stxwx.f
	gcc bsplvd.o bvalue.o bvalus.o spline_fit.o sgram.o sinerp.o sslvrg.o stxwx.o -L/usr/lib/R/lib/ -lR -lblas -lgfortran -lm -o spline_fit
	./spline_fit
