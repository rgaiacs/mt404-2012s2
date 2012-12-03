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
          ! Example from
          ! http://www.ee.ucla.edu/~vandenbe/103/lectures/ls.pdf

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
          n = 2
          A(1, 1) =  2
          A(1, 2) =  0
          A(2, 1) = -1
          A(2, 2) =  1
          A(3, 1) =  0
          A(3, 2) =  2
          b(1, 1) =  1
          b(2, 1) =  0
          b(3, 1) =  -1

          call dshow_matrix(A, lda, m, n)
          !call dgeqrf(m, n, A, lda, tau, work, lwork, info)
          call dshow_matrix(A, lda, m, n)
          call dgels('N', m, n, nrhs, A, lda, b, ldb, work, lwork, info)
          if (info .eq. 0) then
              write (*, *) 'Sucess.'
              call dshow_matrix(b, ldb, n, 1)
          else
              write (*, *) 'Fail.'
          end if
      end program main
