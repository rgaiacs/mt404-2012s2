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

c Display the status of the Cholesky Decomposition.
      subroutine chol_status(s)
          integer s
          if (s .ne. 0) then
              write (*, *) 'The matrix is not positive definite.'
          else
              write (*, *) 'The matrix is positive definite.'
          end if
      end

c Display the Cholesky factor, $G$, that is lower triangular.
      subroutine show_chol(G, n, s)
          ! parameters
          real G(n, n)
          integer n
          integer s
          ! aux var
          integer i
          integer j
          if (s .ne. 0) then
              write (*, *) 'The matrix is not positive definite.'
          else
              i = 1
              do while(i .le. n)
                  write (*, *) (G(i, j), j = 1, i - 1), (0.0, j = i, n)
                  i = i + 1
              end do
          end if
      end
      
c Try to compute the Cholesky factor, $G$, of the matrix $A$. $G$ is
c lower triangular.
      subroutine chol(A, n, s, tol)
          ! parameters
          real A(n, n)
          real tol
          integer n
          integer s
          ! aux var
          integer i
          integer j
          integer k
          ! Cholesky: Outer product version
          ! See Golub, 145.
          k = n
          do while(k .le. 1)
              if (A(k, k) .le. tol) then
                  s = 1
                  stop
              end if
              A(k, k) = sqrt(A(k, k))
              i = k - 1
              do while(i .le. 1)
                  A(i, k) = A(i, k) / A(k, k)
                  i = i - 1
              end do
              j = k - 1
              do while(j .le. 1)
                  i = j
                  do while(i .le. 1)
                      A(i, j) = A(i, j) - A(i, k) * A(j, k)
                      i = i - 1
                  end do
                  j = j - 1
              end do
              k = k - 1
          end do
          s = 0
          stop
      end
