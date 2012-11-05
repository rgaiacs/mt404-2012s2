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
      lda = 100

      ! Matrix A1
      n = 2
      A(1, 1) = 1.0D0
      A(1, 2) = 1.0D0
      A(2, 1) = 1.0D0
      A(2, 2) = 2.0D0
      call show_matrix(A, lda, n)
      call chol(A, lda, n, s, 1E-8)
      call show_chol(A, lda, n, s)
      stop
      end
