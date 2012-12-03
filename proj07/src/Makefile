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

lblas=-lblas
llapack=-llapack

all: mt404_p07

mt404_p07: mt404_p07.f
	gfortran $(DEBUG) mt404_p07.f $(lblas) $(llapack) -o mt404_p07 

sample: sample.f
	gfortran $(DEBUG) sample.f dmatrix.f $(lblas) $(llapack) -o sample 

debug:
	$(MAKE) $(MAKEFILE) DEBUG=-g

benchmark:
	./mt404_p07

dmatrix.so: dmatrix.f
	f2py -c -m dmatrix dmatrix.f

test: dmatrix.so
	python -m doctest dmatrix.test

help:
	@echo 'To compile this files you need a BLAS library and Lapack install.'
	@echo 'To specify where the libraries are install use:'
	@echo '    lblas=-L/path/to/blas'
	@echo '    llapack=-L/path/to/lapack'
	@echo ''
	@echo 'If you need more information on how to install the libaries, see'
	@echo 'http://gcc.gnu.org/wiki/GfortranBuild'
