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

      integer n, lda, s
      real A(100, 100)
      real G(100, 100)
      real tol
      lda = 100
      tol = 1.0E-16
      
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
      call chol(A, lda, n, s, tol)
      call show_chol(A, lda, n, s)

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
      call chol(A, lda, n, s, tol)
      call show_chol(A, lda, n, s)

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
      call chol(A, lda, n, s, tol)
      call show_chol(A, lda, n, s)

      ! Matrix A4
      n = 5
      call rand_tl(G, lda, n)
      call sum_identity(G, lda, n)
      call show_matrix(G, lda, n)
      write (*, *)
      call mt2(A, G, lda, n)
      call show_matrix(A, lda, n)
      call chol(A, lda, n, s, tol)
      call chol_status(s)

      ! Matrix A5
      n = 50
      call rand_tl(G, lda, n)
      call sum_identity(G, lda, n)
      call mt2(A, G, lda, n)
      call chol(A, lda, n, s, tol)
      call chol_status(s)

      ! Matrix A6
      n = 100
      call rand_tl(G, lda, n)
      call sum_identity(G, lda, n)
      call mt2(A, G, lda, n)
      call chol(A, lda, n, s, tol)
      call chol_status(s)
      stop
      end
