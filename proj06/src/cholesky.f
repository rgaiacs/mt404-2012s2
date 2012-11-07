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

      subroutine chol_status(s)
          ! Display in the screen the status of the Cholesky
          ! Decomposition.

          ! parameters
          integer s

          if (s .ne. 0) then
              write (*, *) 'The matrix is not positive definite.'
          else
              write (*, *) 'The matrix is positive definite.'
          end if
      end subroutine chol_status

      subroutine show_chol(G, lda, n, s)
          ! Display in the screen the Cholesky factor, $G : n \times n$,
          ! that is lower triangular.

          ! parameters
          integer n, lda, s
          real G(lda, *)
          ! aux var
          integer i
          integer j

          if (s .ne. 0) then
              write (*, *) 'The matrix is not positive definite.'
          else
              i = 1
              do while (i .le. n)
                  write (*, *) (G(i, j), j = 1, i), (0.0, j = i + 1, n)
                  i = i + 1
              end do
          end if
      end subroutine show_chol
      
      subroutine chol(A, lda, n, s, tol)
          ! Try to compute the Cholesky factor, $G$, of the matrix $A :
          ! n \times n$. $G$ is lower triangular.

          ! parameters
          integer n, lda, s
          real A(lda, *)
          real tol
          ! aux var
          integer i
          integer j
          integer k

          ! Cholesky: Outer product version
          ! See Golub, 145.
          k = 1
          s = 0
          do while (k .le. n)
              if (A(k, k) .le. tol) then
                  s = 1
                  goto 100
              end if
              A(k, k) = sqrt(A(k, k))
              i = k + 1
              do while (i .le. n)
                  A(i, k) = A(i, k) / A(k, k)
                  i = i + 1
              end do
              j = k + 1
              do while (j .le. n)
                  i = j
                  do while (i .le. n)
                      A(i, j) = A(i, j) - A(i, k) * A(j, k)
                      i = i + 1
                  end do
                  j = j + 1
              end do
              k = k + 1
          end do
 100  end subroutine chol

      subroutine solve_with_chol(A, x, b, lda, n, s, tol)
          ! Try to solve the linear system $A x = b$, where $A : n
          ! \times n$, using the Cholesky Decomposition.

          ! parameters
          integer lda, n, s
          real A(lda, *)
          real x(*)
          real b(*)
          real tol
          ! aux var
          real y(lda)

          call chol(A, lda, n, s, tol)
          if (s .eq. 0) then
              call solve_tl(A, y, b, lda, n)
              call transpose_smatrix(A, lda, n)
              call solve_tu(A, x, y, lda, n)
          end if
      end subroutine solve_with_chol
