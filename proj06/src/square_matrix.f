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

      subroutine copy_smatrix(A, A_copy, lda, n)
          ! This function create a copy of $A : n \times n$ in $A_copy$.

          ! arguments
          integer, intent(in) :: n, lda
          double precision, intent (in) :: A(lda, lda)
          double precision, intent (out) :: A_copy(lda, lda)
          !aux var
          integer i, j

          i = 1
          do while (i .le. n)
              j = 1
              do while (j .le. n)
                  A_copy(i, j) = A(i, j)
                  j = j + 1
              end do
              i = i + 1
          end do
      end subroutine copy_smatrix

      subroutine show_smatrix(A, lda, n)
          ! This function pretty print on the screen the matrix $A : n
          ! \times n$.

          ! arguments
          integer, intent(in) :: n, lda
          double precision, intent(in) :: A(lda, *)
          ! aux var
          integer i, j

          i = 1
          do while (i .le. n)
              write (*, *) (A(i, j), j = 1, n)
              i = i + 1
          end do
      end subroutine show_smatrix

      subroutine transpose_smatrix(A, lda, n)
          ! This function transpose the matrix $A : n \times n$ in
          ! place.

          ! arguments
          integer, intent(in) :: lda, n
          double precision, intent(inout) :: A(lda, lda)
          ! aux var
          integer i, j
          double precision temp

          i = 1
          do while (i .le. n)
              j = i + 1
              do while (j .le. n)
                  temp = A(i, j)
                  A(i, j) = A(j, i)
                  A(j, i) = temp
                  j = j + 1
              end do
              i = i + 1
          end do
      end subroutine transpose_smatrix

      subroutine is_symmetric(A, lda, n, s)
          ! This function verify if the matrix $A : n \times n$ is
          ! symmetric.

          ! arguments
          integer, intent(in) :: lda, n
          double precision, intent(in) :: A(lda, lda)
          integer, intent(out) :: s
          ! aux var
          integer i, j

          s = 1
          i = 1
          do while (i .le. n)
              j = i + 1
              do while (j .le. n)
                  if (A(i, j) .ne. A(j, i)) then
                      s = 0
                      go to 100
                  end if
                  j = j + 1
              end do
              i = i + 1
          end do
100   end subroutine is_symmetric

      subroutine rand_tl(A, lda, n)
          ! This function put random numbers in the lower triangular
          ! part of the matrix $A : n \times n$.

          ! arguments
          integer, intent(in) :: n, lda
          double precision, intent(inout) :: A(lda, lda)
          ! aux var
          integer i, j

          i = 1
          do while (i .le. n)
              j = 1
              do while (j .le. n)
                  A(i, j) = rand()
                  j = j + 1
              end do
              i = i + 1
          end do
      end subroutine rand_tl

      subroutine sum_identity(A, lda, n)
          ! This function sum 1 at each element in the main diagonal of
          ! the matrix $A : n \times n$.

          ! arguments
          integer n, lda
          double precision A(lda, *)
          ! aux var
          integer i

          i = 1
          do while (i .le. n)
              A(i, i) = A(i, i) + 1
              i = i + 1
          end do
      end subroutine sum_identity
      
      subroutine mtv(A, x, z, lda, n)
          ! This function compute $z = A x$, where $A : n \times n$ and
          ! $x : n \times 1$.

          ! arguments
          integer n, lda
          double precision A(lda, *)
          double precision x(*)
          double precision z(*)
          ! aux var
          integer i, j

          i = 1
          do while (i .le. n)
             j = 1
             do while (j .le. n)
                 z(i) = z(i) + A(i, j) * x(j)
                 j = j + 1
             end do
             i = i + 1
         end do
      end subroutine mtv

      subroutine mt2(A, G, lda, n)
          ! This function compute $A = G G^t$, where $G : n \times n$.

          ! arguments
          integer n, lda
          double precision A(lda, *)
          double precision G(lda, *)
          ! aux var
          integer i, j, k

          i = 1
          j = 1
          k = 1
          do while (i .le. n)
              j = 1
              do while (j .le. n)
                  A(i, j) = 0
                  k = 1
                  do while (k .le. n)
                      A(i, j) = A(i, j) + G(i, k) * G(j, k)
                      k = k + 1
                  end do
                  j = j + 1
              end do
              i = i + 1
          end do
      end subroutine mt2

      subroutine solve_tl(A, x, b, lda, n)
          ! This function compute $x = A^{-1} b$, where $A : n \times n$
          ! is a lower triangular matrix.

          ! argumensts
          integer, intent(in) :: lda, n
          double precision, intent(in) :: A(lda, lda)
          double precision, intent(in) :: b(lda)
          double precision, intent(out) :: x(lda)
          ! aux var
          integer i, j
          double precision b_temp(lda)

          j = 1
          do while (j .le. n)
              b_temp(j) = b(j)
              j = j + 1
          end do
          j = 1
          do while (j .le. n)
              i = j
              x(j) = b(j) / A(j, j)
              i = i + 1
              do while (i .le. j - 1)
                  b_temp(j) = b_temp(j) - A(i, j) * x(j)
                  i = i + 1
              end do
              j = j + 1
          end do
      end subroutine solve_tl

      subroutine solve_tu(A, x, b, lda, n)
          ! This function compute $x = A^{-1} b$, where $A : n \times n$
          ! is a upper triangular matrix.

          ! arguments
          integer lda, n
          double precision A(lda, *)
          double precision x(*)
          double precision b(*)
          ! aux var
          integer i, j
          double precision b_temp(lda)

          j = 1
          do while (j .le. n)
              b_temp(j) = b(j)
              j = j + 1
          end do
          j = n
          do while (j .ge. 1)
              i = n
              do while (i .ge. j + 1)
                  b_temp(j) = b_temp(j) - A(i, j) * x(j)
                  i = i - 1
              end do
              x(j) = b(j) / A(j, j)
              j = j - 1
          end do
      end subroutine solve_tu
