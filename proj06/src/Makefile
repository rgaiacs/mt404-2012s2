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

all: compile

compile:
	gfortran $(DEBUG) -o sample sample.f cholesky.f square_matrix.f \
	vector.f
	gfortran $(DEBUG) -o mt404_p06 mt404_p06.f cholesky.f \
	square_matrix.f vector.f

debug:
	$(MAKE) $(MAKEFILE) DEBUG="-g"

test:
	./sample > sample.out
	diff sample.exp sample.out

benchmark:
	./mt404_p06 > mt404_p06.csv
	- column -s \, -t mt404_p06.csv
