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

      program main
          integer n, lda, s, n_rand, i
          integer rand_dim(100)
          real A(100, 100)
          real G(100, 100)
          real tol, erro, res
          lda = 100
          tol = 1.0E-8
          
          write (*, *) 'erro, res, opt, dim'
  900     format (e10.4, ',', e10.4, ',', i4, ',', i4)
          ! Matrix A1
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
          call test_matrix(A, lda, n, tol, erro, res)
          write (*, 900) erro, res, 1, n

          ! Matrix A2
          n = 5
          A(1, 1) =  1
          A(1, 2) = -1
          A(1, 3) =  1
          A(1, 4) = -1
          A(1, 5) =  1
          A(2, 1) = -1
          A(2, 2) =  2
          A(2, 3) = -2
          A(2, 4) =  2
          A(2, 5) = -2
          A(3, 1) =  1
          A(3, 2) = -2
          A(3, 3) =  3
          A(3, 4) = -3
          A(3, 5) =  3
          A(4, 1) = -1
          A(4, 2) =  2
          A(4, 3) = -3
          A(4, 4) =  4
          A(4, 5) = -4
          A(5, 1) =  1
          A(5, 2) = -2
          A(5, 3) =  3
          A(5, 4) = -4
          A(5, 5) =  5
          call test_matrix(A, lda, n, tol, erro, res)
          write (*, 900) erro, res, 2, n

          ! Matrix A3
          n = 3
          A(1, 1) =   3
          A(1, 2) =  -6
          A(1, 3) =   9
          A(2, 1) =  -6
          A(2, 2) =  14
          A(2, 3) = -20
          A(3, 1) =   9
          A(3, 2) = -20
          A(3, 3) =  29
          call test_matrix(A, lda, n, tol, erro, res)
          write (*, 900) erro, res, 3, n

          nrand = 3
          rand_dim(1) = 5
          rand_dim(2) = 50
          rand_dim(3) = 100

          i = 1
          do while (i .le. nrand)
              call rand_tl(G, lda, n)
              call sum_identity(G, lda, n)
              call mt2(A, G, lda, n)
              call test_matrix(A, lda, n, tol, erro, res)
              write (*, 900) erro, res, 3, n
              i =  i + 1
          end do
      end program main

      subroutine test_matrix(A, lda, n, tol, erro, res)
          ! This function test the Cholesky decomposition to solve a
          ! linear system $A x = b$, where $A : n \times n$.

          ! parameters
          integer n, lda
          real A(lda, *)
          real tol, erro, res
          ! aux var
          integer state
          real A_copy(lda, lda)
          real s(lda)
          real b(lda)
          real x(lda)
          real aux(lda)
          real num, den

          call rand_v(s, n)
          call mtv(A, s, b, lda, n)
          call copy_matrix(A, A_copy, lda, n)
          call solve_with_chol(A, x, b, lda, n, state, tol)
          if (state .eq. 0) then
              call vpov(x, s, aux, n)
              call vnorm_inf(aux, n, num)
              call vnorm_inf(x, n, den)
              erro = num / den
              call mtv(A_copy, x, s, lda, n)
              call vpov(b, s, aux, n)
              call vnorm_inf(aux, n, num)
              call vnorm_inf(b, n, den)
              res = num / den
          end if
      end subroutine test_matrix
