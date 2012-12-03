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
          integer m, n
          double precision r

          call test01(50, 10, r)
          write (*, *) r
          call test02(50, 10, r)
          write (*, *) r
          call test03(50, 10, r)
          write (*, *) r
      end program main

      subroutine test01(m, n, r)
          ! This function match the test 1.

          ! parameters
          integer, intent(in) :: m, n
          double precision, intent(out) :: r
          ! test var
          integer lda, ldb
          double precision A(1000, 1000)
          double precision b(1000, 1)
          double precision tau(1000, 1000)
          double precision work(1000)
          integer lwork, nrhs, info
          ! aux var
          integer i, j

          lda = 1000
          ldb = 1000
          nrhs = 1
          lwork = 1000

          i = 1
          do while (i .le. n)
              b(i, 1) = 1
              j = 1
              do while (j .le. n)
                  A(i, j) = - 2 / m
                  j = j + 1
              end do
              A(i, i) = A(i, i) + 1
              i = i + 1
          end do
          i =  n + 1
          do while (i .le. m)
              b(i, 1) = 1
              j = 1
              do while (j .le. n)
                  A(i, j) = - 2 / m
                  j = j + 1
              end do
              i = i + 1
          end do
          call dgels('N', m, n, nrhs, A, lda, b, ldb, work, lwork, info)

          r = 0
          i = n + 1
          do while (i .le. m)
              r = r + b(i, 1) ** 2
              i = i + 1
          end do
      end subroutine test01

      subroutine test02(m, n, r)
          ! This function match the test 1.

          ! parameters
          integer, intent(in) :: m, n
          double precision, intent(out) :: r
          ! test var
          integer lda, ldb
          double precision A(1000, 1000)
          double precision b(1000, 1)
          double precision tau(1000, 1000)
          double precision work(1000)
          integer lwork, nrhs, info
          ! aux var
          integer i, j

          lda = 1000
          ldb = 1000
          nrhs = 1
          lwork = 1000

          i = 1
          do while (i .le. m)
              b(i, 1) = 1
              j = 1
              do while (j .le. n)
                  A(i, j) = i * j
                  j = j + 1
              end do
              i = i + 1
          end do
          call dgels('N', m, n, nrhs, A, lda, b, ldb, work, lwork, info)

          r = 0
          i = n + 1
          do while (i .le. m)
              r = r + b(i, 1) ** 2
              i = i + 1
          end do
      end subroutine test02

      subroutine test03(m, n, r)
          ! This function match the test 1.

          ! parameters
          integer, intent(in) :: m, n
          double precision, intent(out) :: r
          ! test var
          integer lda, ldb
          double precision A(1000, 1000)
          double precision b(1000, 1)
          double precision tau(1000, 1000)
          double precision work(1000)
          integer lwork, nrhs, info
          ! aux var
          integer i, j

          lda = 1000
          ldb = 1000
          nrhs = 1
          lwork = 1000


          j = 1
          b(1, 1) = 1
          do while (j .le. n)
              A(1, j) = 0
              j = j + 1
          end do
          i = 2
          do while (i .lt. m)
              b(i, 1) = 1
              j = 1
              do while (j .le. n)
                  A(i, j) = (i - 1) * j
                  j = j + 1
              end do
              i = i + 1
          end do
          j = 1
          b(m, 1) = 1
          do while (j .le. n)
              A(m, j) = 0
              j = j + 1
          end do
          call dgels('N', m, n, nrhs, A, lda, b, ldb, work, lwork, info)

          r = 0
          i = n + 1
          do while (i .le. m)
              r = r + b(i, 1) ** 2
              i = i + 1
          end do
      end subroutine test03
