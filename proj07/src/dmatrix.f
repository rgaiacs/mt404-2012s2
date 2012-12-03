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

      subroutine dcopy_matrix(A, A_copy, lda, m, n)
          ! This function create a copy of $A : m \times n$ in $A_copy$.

          ! parameters
          integer, intent(in) :: m, n, lda
          double precision, intent(in) :: A(lda, lda)
          double precision, intent(out) :: A_copy(lda, lda)
          ! aux var
          integer i, j

          i = 1
          do while (i .le. m)
              j = 1
              do while (j .le. n)
                  A_copy(i, j) = A(i, j)
                  j = j + 1
              end do
              i = i + 1
          end do
      end subroutine dcopy_matrix

      subroutine dshow_matrix(A, lda, m, n)
          ! This function pretty print on the screen the matrix $A : m
          ! \times n$.

          ! parameters
          integer, intent(in) :: m, n, lda
          double precision, intent(in) :: A(lda, *)
          ! aux var
          integer i, j

          i = 1
          do while (i .le. m)
              write (*, *) (A(i, j), j = 1, n)
              i = i + 1
          end do
      end subroutine dshow_matrix
