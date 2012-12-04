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
          integer n, lda, state, i
          integer rand_dim(100)
          double precision A(100, 100)
          double precision G(100, 100)
          double precision tol, erro, res
          lda = 100
          tol = 1.0E-8
          
          write (*, *) "Teste,Tipo,Dimens\~{a}o,Erro Relativo"//
     .",Res\'{i}duo Relativo"
          ! Matrix A1
          n = 3
          A(1, 1) = 36.0
          A(1, 2) = 30.0
          A(1, 3) = 24.0
          A(2, 1) = 30.0
          A(2, 2) = 34.0
          A(2, 3) = 26.0
          A(3, 1) = 24.0
          A(3, 2) = 26.0
          A(3, 3) = 21.0
          call test_matrix(A, lda, n, tol, state, erro, res)
          call write_info(state, 1, 1, n, erro, res)

          ! Matrix A2
          n = 5
          A(1, 1) =  1.0
          A(1, 2) = -1.0
          A(1, 3) =  1.0
          A(1, 4) = -1.0
          A(1, 5) =  1.0
          A(2, 1) = -1.0
          A(2, 2) =  2.0
          A(2, 3) = -2.0
          A(2, 4) =  2.0
          A(2, 5) = -2.0
          A(3, 1) =  1.0
          A(3, 2) = -2.0
          A(3, 3) =  3.0
          A(3, 4) = -3.0
          A(3, 5) =  3.0
          A(4, 1) = -1.0
          A(4, 2) =  2.0
          A(4, 3) = -3.0
          A(4, 4) =  4.0
          A(4, 5) = -4.0
          A(5, 1) =  1.0
          A(5, 2) = -2.0
          A(5, 3) =  3.0
          A(5, 4) = -4.0
          A(5, 5) =  5.0
          call test_matrix(A, lda, n, tol, state, erro, res)
          call write_info(state, 2, 2, n, erro, res)

          ! Matrix A3
          n = 3
          A(1, 1) =   3.0
          A(1, 2) =  -6.0
          A(1, 3) =   9.0
          A(2, 1) =  -6.0
          A(2, 2) =  14.0
          A(2, 3) = -20.0
          A(3, 1) =   9.0
          A(3, 2) = -20.0
          A(3, 3) =  29.0
          call test_matrix(A, lda, n, tol, state, erro, res)
          call write_info(state, 3, 3, n, erro, res)

          nrand = 3
          rand_dim(1) = 5
          rand_dim(2) = 50
          rand_dim(3) = 100

          i = 1
          do while (i .le. nrand)
              call rand_tl(G, lda, n)
              call sum_identity(G, lda, n)
              call mt2(A, G, lda, n)
              call test_matrix(A, lda, n, tol, state, erro, res)
              call write_info(state, i + 3, 5, rand_dim(i), erro, res)
              i =  i + 1
          end do
      end program main

      subroutine test_matrix(A, lda, n, tol, state, erro, res)
          ! This function test the Cholesky decomposition to solve a
          ! linear system $A x = b$, where $A : n \times n$.

          ! parameters
          integer n, lda, state
          double precision A(lda, *)
          double precision tol, erro, res
          ! aux var
          double precision A_copy(lda, lda)
          double precision s(lda)
          double precision b(lda)
          double precision x(lda)
          double precision aux(lda)
          double precision num, den

          call rand_v(s, n)
          call mtv(A, s, b, lda, n)
          call copy_smatrix(A, A_copy, lda, n)
          ! We need copy the matrix $A$ because it will be lost when
          ! copute the Cholesky factor.
          call solve_with_chol(A, x, b, lda, n, state, tol)
          if (state .eq. 0) then
              call vpov(x, s, aux, n)
              call vnorm_inf(aux, n, num)
              call vnorm_inf(x, n, den)
              erro = num / den
              ! Once we will not need $s$ again, we can reuse it.
              call mtv(A_copy, x, s, lda, n)
              call vpov(b, s, aux, n)
              call vnorm_inf(aux, n, num)
              call vnorm_inf(b, n, den)
              res = num / den
          end if
      end subroutine test_matrix

      subroutine write_info(state, id, test, n, erro, res)
          ! Write the information about a test.

          !parameters
          integer state, id, test, n
          double precision erro, res

900       format (i4, ',', i4, ',', i4, ',',  e10.4, ',', e10.4)
901       format (i4, ',', i4, ',', i4, ',inf,inf')
          if (state .eq. 0) then
              write (*, 900) id, test, n, erro, res
          else
              write (*, 901) id, test, n
          end if
      end subroutine write_info
