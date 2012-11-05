all: compile test

compile:
	gfortran -o test test.f cholesky.f matrix.f
	gfortran -o mt404_p06 mt404_p06.f cholesky.f matrix.f

test:
	./test
