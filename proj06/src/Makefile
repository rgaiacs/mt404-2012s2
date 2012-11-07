FFLAG=-g

all: compile

compile:
	gfortran $(FFLAG) -o sample sample.f cholesky.f square_matrix.f \
	vector.f
	gfortran $(FFLAG) -o mt404_p06 mt404_p06.f cholesky.f \
	square_matrix.f vector.f

test:
	./sample > sample.out
	diff sample.exp sample.out
