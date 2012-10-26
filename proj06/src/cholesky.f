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

       subroutine chol(A, n)
           ! parameters
           integer n
           real A(n, n)
           ! aux var
           integer i
           integer j
           integer k
           ! Cholesky: Outer product version
           ! See Golub, 145.
           k = 1
           do while(k .le. n)
               A(k, k) = sqrt(A(k, k))
               i = k + 1
               do while(i .le. n)
                   A(i, k) = A(i, k) / A(k, k)
                   i = i + 1
               end do
               j = k + 1
               do while(j .le. n)
                   i = j
                   do while(i .le. n)
                       A(i, j) = A(i, j) - A(i, k) * A(j, k)
                       i = i + 1
                   end do
                   j = j + 1
               end do
               k = k + 1
           end do
       end
