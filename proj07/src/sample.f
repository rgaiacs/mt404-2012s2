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
          double precision A(100, 100)
          double precision b(100, 1)
          double precision tau(100, 100)
          double precision work(100)
          integer lda, ldb, lwork, m, n, nrhs, info

          lda = 100
          ldb = 100
          nrhs = 1
          lwork = 100
          m = 3
          n = 3
          A(1, 1) = 12
          A(1, 2) = -51
          A(1, 3) = 4
          A(2, 1) = 6
          A(2, 2) = 167
          A(2, 3) = -68
          A(3, 1) = -4
          A(3, 2) = 24
          A(3, 3) = -41

          call dshow_matrix(A, lda, n)
          call dgeqrf(m, n, A, lda, tau, work, lwork, info)
          call dshow_matrix(A, lda, n)
          !call dgels('N', m, n, nrhs, A, lda, b, ldb, work, lwork, info)
      end program main
