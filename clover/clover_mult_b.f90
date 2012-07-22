!===============================================================================
!
! clover_mult_b.F90
!
!-------------------------------------------------------------------------------
!
! Copyright (C) 1998-2001 Hinnerk Stueben
!
! This file is part of BQCD -- Berlin Quantum ChromoDynamics program
!
! BQCD is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! BQCD is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with BQCD. If not, see <http://www.gnu.org/licenses/>.
!
!-------------------------------------------------------------------------------




!-------------------------------------------------------------------------------
subroutine clover_mult_b(b, x, volh) ! x := B x

  implicit none

  complex(8), dimension(18, 2, *) :: b
  complex(8), dimension(4, 3, *) :: x
  integer :: volh

  integer :: i
  complex(8) :: x1, x2, x3, x4, x5, x6
  complex(8) :: y1, y2, y3, y4, y5, y6

  call timing_start(32)

  !$omp parallel do private(x1, x2, x3, x4, x5, x6, y1, y2, y3, y4, y5, y6)
  do i = 1, volh

     y1 = x(1, 1, i) + x(3, 1, i)
     y2 = x(1, 2, i) + x(3, 2, i)
     y3 = x(1, 3, i) + x(3, 3, i)
     y4 = x(2, 1, i) + x(4, 1, i)
     y5 = x(2, 2, i) + x(4, 2, i)
     y6 = x(2, 3, i) + x(4, 3, i)
     y2 = y2 - b(1,1,i) * y1
     y3 = y3 - b(2,1,i) * y1 - b(3,1,i) * y2
     y4 = y4 - b(4,1,i) * y1 - b(5,1,i) * y2 - b(6,1,i) * y3
     y5 = y5 - b(7,1,i) * y1 - b(8,1,i) * y2 - b(9,1,i) * y3 - b(10,1,i) * y4
     y6 = y6 - b(11,1,i) * y1 - b(12,1,i) * y2 - b(13,1,i) * y3 - b(14,1,i) * y4 - b(15,1,i) * y5
     x6 = y6 * aimag(b(18,1,i))
     x5 = y5 * real(b(18,1,i)) - conjg(b(15,1,i)) * x6
     x4 = y4 * aimag(b(17,1,i)) - conjg(b(10,1,i)) * x5 - conjg(b(14,1,i)) * x6
     x3 = y3 * real(b(17,1,i)) - conjg(b(6,1,i)) * x4 - conjg(b(9,1,i)) * x5 &
        - conjg(b(13,1,i)) * x6
     x2 = y2 * aimag(b(16,1,i)) - conjg(b(3,1,i)) * x3 - conjg(b(5,1,i)) * x4 &
        - conjg(b(8,1,i)) * x5 - conjg(b(12,1,i)) * x6
     x1 = y1 * real(b(16,1,i)) - conjg(b(1,1,i)) * x2 - conjg(b(2,1,i)) * x3 &
        - conjg(b(4,1,i)) * x4 - conjg(b(7,1,i)) * x5 - conjg(b(11,1,i)) * x6
     y1 = x(1, 1, i) - x(3, 1, i)
     y2 = x(1, 2, i) - x(3, 2, i)
     y3 = x(1, 3, i) - x(3, 3, i)
     y4 = x(2, 1, i) - x(4, 1, i)
     y5 = x(2, 2, i) - x(4, 2, i)
     y6 = x(2, 3, i) - x(4, 3, i)
     x(1, 1, i) = x1
     x(1, 2, i) = x2
     x(1, 3, i) = x3
     x(2, 1, i) = x4
     x(2, 2, i) = x5
     x(2, 3, i) = x6
     x(3, 1, i) = x1
     x(3, 2, i) = x2
     x(3, 3, i) = x3
     x(4, 1, i) = x4
     x(4, 2, i) = x5
     x(4, 3, i) = x6
     y2 = y2 - b(1,2,i) * y1
     y3 = y3 - b(2,2,i) * y1 - b(3,2,i) * y2
     y4 = y4 - b(4,2,i) * y1 - b(5,2,i) * y2 - b(6,2,i) * y3
     y5 = y5 - b(7,2,i) * y1 - b(8,2,i) * y2 - b(9,2,i) * y3 - b(10,2,i) * y4
     y6 = y6 - b(11,2,i) * y1 - b(12,2,i) * y2 - b(13,2,i) * y3 - b(14,2,i) * y4 - b(15,2,i) * y5
     x6 = y6 * aimag(b(18,2,i))
     x5 = y5 * real(b(18,2,i)) - conjg(b(15,2,i)) * x6
     x4 = y4 * aimag(b(17,2,i)) - conjg(b(10,2,i)) * x5 - conjg(b(14,2,i)) * x6
     x3 = y3 * real(b(17,2,i)) - conjg(b(6,2,i)) * x4 - conjg(b(9,2,i)) * x5 &
        - conjg(b(13,2,i)) * x6
     x2 = y2 * aimag(b(16,2,i)) - conjg(b(3,2,i)) * x3 - conjg(b(5,2,i)) * x4 &
        - conjg(b(8,2,i)) * x5 - conjg(b(12,2,i)) * x6
     x1 = y1 * real(b(16,2,i)) - conjg(b(1,2,i)) * x2 - conjg(b(2,2,i)) * x3 &
        - conjg(b(4,2,i)) * x4 - conjg(b(7,2,i)) * x5 - conjg(b(11,2,i)) * x6
     x(1, 1, i) = x(1, 1, i) + x1
     x(1, 2, i) = x(1, 2, i) + x2
     x(1, 3, i) = x(1, 3, i) + x3
     x(2, 1, i) = x(2, 1, i) + x4
     x(2, 2, i) = x(2, 2, i) + x5
     x(2, 3, i) = x(2, 3, i) + x6
     x(3, 1, i) = x(3, 1, i) - x1
     x(3, 2, i) = x(3, 2, i) - x2
     x(3, 3, i) = x(3, 3, i) - x3
     x(4, 1, i) = x(4, 1, i) - x4
     x(4, 2, i) = x(4, 2, i) - x5
     x(4, 3, i) = x(4, 3, i) - x6
  enddo
  call timing_stop(32)
end
!===============================================================================
