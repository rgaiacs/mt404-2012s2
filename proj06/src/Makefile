all: compile

compile:
	gfortran -g -o sample sample.f cholesky.f matrix.f
	gfortran -g -o mt404_p06 mt404_p06.f cholesky.f matrix.f

test:
	./sample > sample.out
	diff sample.exp sample.out
