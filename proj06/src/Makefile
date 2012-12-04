# Copyright (C) 2012 Raniere Silva <r.gaia.cs@gmail.com>
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or (at
# your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with Octave; see the file COPYING.  If not, see
# <http://www.gnu.org/licenses/>.

# To enable debug mode, use
#
#     $ make DEBUG=-g <args>
#
# To change the flags passed to gfortran, use
#
#     $ make FFLAGS=-ansi
FFLAGS=-pedantic -Wall -Waliasing -Wconversion-extra -Werror

benchmark: mt404_p06
	./mt404_p06 > mt404_p06.csv
	- column -s \, -t mt404_p06.csv

mt404_p06: mt404_p06.f cholesky.f square_matrix.f
	gfortran $(DEBUG) $(FFLAGS) -o mt404_p06 mt404_p06.f cholesky.f \
	square_matrix.f vector.f

test: vector.so square_matrix.so cholesky.so
	python -m doctest $(VERBOSE) vector.test
	python -m doctest $(VERBOSE) square_matrix.test
	python -m doctest $(VERBOSE) cholesky.test

vector.so: vector.f
	f2py -m vector -h vector.pyf vector.f
	sed -i 's/inout/in,out/g' vector.pyf
	f2py -c vector.pyf vector.f

square_matrix.so: square_matrix.f
	f2py -m square_matrix -h square_matrix.pyf square_matrix.f
	sed -i 's/inout/in,out/g' square_matrix.pyf
	f2py -c square_matrix.pyf square_matrix.f

cholesky.so: cholesky.f
	f2py -m cholesky -h cholesky.pyf cholesky.f square_matrix.f
	sed -i 's/inout/in,out/g' cholesky.pyf
	f2py -c cholesky.pyf cholesky.f square_matrix.f

clean:
	- rm *.so
	- rm *.pyf
	- rm mt404_p06

help:
	@echo 'This is the Makefile of Proj06.'
	@echo ''
	@echo 'Usage:'
	@echo '    - benchmark: compile and run the benchmark,'
	@echo "    - test: use f2py and Python's doctest to check for"
	@echo "    bugs in the fortran files.                        "
	@echo '    - clean: remove all aux files.'
