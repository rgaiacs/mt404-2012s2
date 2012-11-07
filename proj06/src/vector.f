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

      subroutine rand_v(x, n)
          ! This function put random numbers in the vector $x : n \times
          ! 1$.

          ! parameters
          integer n
          real x(n)
          ! aux var
          integer i

          i = 1
          do while (i .le. n)
              x(i) = rand()
              i = i + 1
          end do
      end subroutine rand_v

      subroutine vpv(x, y, z, n)
          ! This function compute $z = x + y$, where x, y, z are vectors
          ! of size n.

          ! parameters
          integer n
          real x(*)
          real y(*)
          real z(*)

          i = 1
          do while (i .le. n)
              z(i) = x(i) + y(i)
              i = i + 1
          end do
      end subroutine vpv

      subroutine vpov(x, y, z, n)
          ! This function compute $z = x - y$, where x, y, z are vectors
          ! of size n.

          ! parameters
          integer n
          real x(*)
          real y(*)
          real z(*)

          i = 1
          do while (i .le. n)
              z(i) = x(i) - y(i)
              i = i + 1
          end do
      end subroutine vpov

      subroutine vnorm_inf(x, n, norm)
          ! This function compute $norm = \| x \|_\infty$, where $x : n
          ! \times 1$.

          ! parameters
          integer n
          real x(*)
          real norm
          !aux var
          integer i

          i = 1
          norm = x(i)
          do while (i .le. n)
              if (norm .le. x(i)) then
                  norm = x(i)
              end if
              i = i + 1
          end do
      end subroutine vnorm_inf
