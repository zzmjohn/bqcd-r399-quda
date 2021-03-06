!===============================================================================
!
! clover_mult_a.F90
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
subroutine clover_mult_a(out, a, in, volh) ! out := A in

  implicit none

  complex(8), dimension(18, 2, *) :: a
  complex(8), dimension(4, 3, *) :: out, in
  integer :: volh

  integer :: i
  complex(8) :: x1, x2, x3, x4, x5, x6
  complex(8) :: y1, y2, y3, y4, y5, y6

  call timing_start(30)

  !$omp parallel do private(x1, x2, x3, x4, x5, x6, y1, y2, y3, y4, y5, y6)
  do i = 1, volh
     x1 = in(1, 1, i) + in(3, 1, i)
     x2 = in(1, 2, i) + in(3, 2, i)
     x3 = in(1, 3, i) + in(3, 3, i)
     x4 = in(2, 1, i) + in(4, 1, i)
     x5 = in(2, 2, i) + in(4, 2, i)
     x6 = in(2, 3, i) + in(4, 3, i)
     y1 = real(a(1,1,i)) * x1 + a(2,1,i) * x2 + a(3,1,i) * x3 &
        + a(4,1,i) * x4 + a(5,1,i) * x5 + a(6,1,i) * x6
     y2 = conjg(a(2,1,i)) * x1 + aimag(a(1,1,i)) * x2 + a(7,1,i) * x3 &
        + a(8,1,i) * x4 + a(9,1,i) * x5 + a(10,1,i) * x6
     y3 = conjg(a(3,1,i)) * x1 + conjg(a(7,1,i)) * x2 + real(a(11,1,i)) * x3 &
        + a(12,1,i) * x4 + a(13,1,i) * x5 + a(14,1,i) * x6
     y4 = conjg(a(4,1,i)) * x1 + conjg(a(8,1,i)) * x2 + conjg(a(12,1,i)) * x3 &
        + aimag(a(11,1,i)) * x4 + a(15,1,i) * x5 + a(16,1,i) * x6
     y5 = conjg(a(5,1,i)) * x1 + conjg(a(9,1,i)) * x2 + conjg(a(13,1,i)) * x3 &
        + conjg(a(15,1,i)) * x4 + real(a(17,1,i)) * x5 + a(18,1,i) * x6
     y6 = conjg(a(6,1,i)) * x1 + conjg(a(10,1,i)) * x2 + conjg(a(14,1,i)) * x3 &
        + conjg(a(16,1,i)) * x4 + conjg(a(18,1,i)) * x5 + aimag(a(17,1,i)) * x6
     out(1, 1, i) = y1
     out(1, 2, i) = y2
     out(1, 3, i) = y3
     out(2, 1, i) = y4
     out(2, 2, i) = y5
     out(2, 3, i) = y6
     out(3, 1, i) = y1
     out(3, 2, i) = y2
     out(3, 3, i) = y3
     out(4, 1, i) = y4
     out(4, 2, i) = y5
     out(4, 3, i) = y6
     x1 = in(1, 1, i) - in(3, 1, i)
     x2 = in(1, 2, i) - in(3, 2, i)
     x3 = in(1, 3, i) - in(3, 3, i)
     x4 = in(2, 1, i) - in(4, 1, i)
     x5 = in(2, 2, i) - in(4, 2, i)
     x6 = in(2, 3, i) - in(4, 3, i)
     y1 = real(a(1,2,i)) * x1 + a(2,2,i) * x2 + a(3,2,i) * x3 &
        + a(4,2,i) * x4 + a(5,2,i) * x5 + a(6,2,i) * x6
     y2 = conjg(a(2,2,i)) * x1 + aimag(a(1,2,i)) * x2 + a(7,2,i) * x3 &
        + a(8,2,i) * x4 + a(9,2,i) * x5 + a(10,2,i) * x6
     y3 = conjg(a(3,2,i)) * x1 + conjg(a(7,2,i)) * x2 + real(a(11,2,i)) * x3 &
        + a(12,2,i) * x4 + a(13,2,i) * x5 + a(14,2,i) * x6
     y4 = conjg(a(4,2,i)) * x1 + conjg(a(8,2,i)) * x2 + conjg(a(12,2,i)) * x3 &
        + aimag(a(11,2,i)) * x4 + a(15,2,i) * x5 + a(16,2,i) * x6
     y5 = conjg(a(5,2,i)) * x1 + conjg(a(9,2,i)) * x2 + conjg(a(13,2,i)) * x3 &
        + conjg(a(15,2,i)) * x4 + real(a(17,2,i)) * x5 + a(18,2,i) * x6
     y6 = conjg(a(6,2,i)) * x1 + conjg(a(10,2,i)) * x2 + conjg(a(14,2,i)) * x3 &
        + conjg(a(16,2,i)) * x4 + conjg(a(18,2,i)) * x5 + aimag(a(17,2,i)) * x6
     out(1, 1, i) = out(1, 1, i) + y1
     out(1, 2, i) = out(1, 2, i) + y2
     out(1, 3, i) = out(1, 3, i) + y3
     out(2, 1, i) = out(2, 1, i) + y4
     out(2, 2, i) = out(2, 2, i) + y5
     out(2, 3, i) = out(2, 3, i) + y6
     out(3, 1, i) = out(3, 1, i) - y1
     out(3, 2, i) = out(3, 2, i) - y2
     out(3, 3, i) = out(3, 3, i) - y3
     out(4, 1, i) = out(4, 1, i) - y4
     out(4, 2, i) = out(4, 2, i) - y5
     out(4, 3, i) = out(4, 3, i) - y6
  enddo
  call timing_stop(30)
end
!===============================================================================
