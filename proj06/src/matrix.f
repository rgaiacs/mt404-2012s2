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

       subroutine show_matrix(A, n)
           ! parameters
           integer n
           real A(n, n)
           ! aux var
           integer i
           integer j
           i = 1
           do while(i .le. n)
               write (*, *) (A(i, j), j = 1,n)
               i = i + 1
           end do
       end

       subroutine show_chol(A, n)
           ! parameters
           integer n
           real A(n, n)
           ! aux var
           integer i
           integer j
           i = 1
           do while(i .le. n)
               write (*, *) (0, j = 1,i-1) (A(i, j), j = i,n)
               i = i + 1
           end do
       end
