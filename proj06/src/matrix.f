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

      subroutine rand_tl(A, n)
          ! parameters
          integer n
          real A(n, n)
          ! aux var
          integer i
          integer j
          i = 1
          do while(i .le. n)
              do while(j .le. n)
                  A(i, j) = rand()
                  j = j + 1
              end do
              i = i + 1
          end do
      end

      subroutine sum_identity(A, n)
          ! parameters
          integer n
          real A(n, n)
          ! aux var
          integer i
          i = 1
          do while(i .le. n)
              A(i, i) = A(i, i) + 1
              i = i + 1
          end do
      end
      
      subroutine mt2(A, G, n)
          ! This function return A = G G^t
          ! parameters
          integer n
          real A(n, n)
          real G(n, n)
          ! aux var
          integer i
          integer j
          integer k
          i = 1
          j = 1
          k = 1
          do while(i .le. n)
              do while(j .le. n)
                  A(i, j) = 0
                  do while(k .le. n)
                      A(i, j) = A(i, j) + G(i, k) * G(j, k)
                      k = k + 1
                  end do
                  j = j + 1
              end do
              i = i + 1
          end do
      end
