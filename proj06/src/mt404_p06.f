c Copyright (C) 2012 Raniere Silva <r.gaia.cs@gmail.com>
c 
c This program is free software; you can redistribute it and/or modify
c it under the terms of the GNU General Public License as published by
c the Free Software Foundation; either version 3 of the License, or (at
c your option) any later version.
c 
c This program is distributed in the hope that it will be useful, but
c WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
c General Public License for more details.
c 
c You should have received a copy of the GNU General Public License
c along with Octave; see the file COPYING.  If not, see
c <http://www.gnu.org/licenses/>.

       integer n
       real A(3, 3)
       n = 3
       A(1, 1) = 36
       A(1, 2) = 30
       A(1, 3) = 24
       A(2, 1) = 30
       A(2, 2) = 34
       A(2, 3) = 26
       A(3, 1) = 24
       A(3, 2) = 26
       A(3, 3) = 21
       call show_matrix(A, n)
       call chol(A, n)
       call show_chol(A, n)
       stop
       end
