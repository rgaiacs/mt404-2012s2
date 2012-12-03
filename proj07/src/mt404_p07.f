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
          call bench(1)
          call bench(2)
          call bench(3)
      end program main

      subroutine test01(m, n, r)
          ! This function match the test 1.

          ! parameters
          parameter (lda=1000, ldb=1000, lwork=1000)
          ! arguments
          integer, intent(in) :: m, n
          double precision, intent(out) :: r
          ! test var
          double precision A(lda, lda)
          double precision b(ldb, 1)
          double precision work(lwork)
          integer nrhs, info
          ! aux var
          integer i, j

          nrhs = 1

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
          parameter (lda=1000, ldb=1000, lwork=1000)
          ! arguments
          integer, intent(in) :: m, n
          double precision, intent(out) :: r
          ! test var
          double precision A(lda, lda)
          double precision b(ldb, 1)
          double precision work(lwork)
          integer nrhs, info
          ! aux var
          integer i, j

          nrhs = 1

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
          parameter (lda=1000, ldb=1000, lwork=1000)
          ! arguments
          integer, intent(in) :: m, n
          double precision, intent(out) :: r
          ! test var
          double precision A(lda, lda)
          double precision b(ldb, 1)
          double precision work(lwork)
          integer nrhs, info
          ! aux var
          integer i, j

          nrhs = 1

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

      subroutine bench(test_type)
          ! arguments
          integer test_type
          ! aux var
          integer i, j, file_id
          integer n_v2test, ldv2test
          double precision r
          integer v2test(100)

          ldv2test = 100
          n_v2test = 4
          v2test(1) = 10
          v2test(2) = 50
          v2test(3) = 100
          v2test(4) = 500

          select case(test_type)
              case(1)
                  write (*, *) 'Try to open file to test01.csv'
                  open (unit = 1, file = "test01.csv")
                  write (*, *) 'Success in open file to test01.csv'
              case(2)
                  write (*, *) 'Try to open file to test02.csv'
                  open (unit = 1, file = "test02.csv")
                  write (*, *) 'Success in open file to test02.csv'
              case(3)
                  write (*, *) 'Try to open file to test03.csv'
                  open (unit = 1, file = "test03.csv")
                  write (*, *) 'Success in open file to test03.csv'
          end select
          i = 1
          write (1, '(a,$)') ','
          do while (i .lt. n_v2test)
              write (1, '(I4,a,$)') v2test(i), ','
              i = i + 1
          end do
          write (1, '(I4)') v2test(n_v2test)

          i = 1
          do while (i .le. n_v2test)
              write (1, '(I4,a,$)') v2test(i), ','
              j = 1
              do while (j .le. i)
                  write (*, *) 'Run test for m = ', i, ' n = ', j
                  select case(test_type)
                      case(1)
                          call test01(v2test(i), v2test(j), r)
                      case(2)
                          call test02(v2test(i), v2test(j), r)
                      case(3)
                          call test03(v2test(i), v2test(j), r)
                  end select
                  if (j .eq. n_v2test) then
                      write (1, '(ES10.4)') r
                  else
                      write (1, '(ES10.4,a,$)') r, ','
                  end if
                  j = j + 1
              end do
              do while (j .le. n_v2test)
                  if (j .eq. n_v2test) then
                      write (1, '(a)') '--'
                  else
                      write (1, '(a,$)') '--,'
                  end if
                  j = j + 1
              end do
              i = i + 1
          end do
          close(1)
      end subroutine bench
