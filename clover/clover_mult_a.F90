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
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with BQCD.  If not, see <http://www.gnu.org/licenses/>.
!
!-------------------------------------------------------------------------------
# define CLOVER_AS_COMPLEX_ARRAY
# include "defs.h"
# include "clover.h"

!-------------------------------------------------------------------------------
subroutine clover_mult_a(out, a, in, volh)  ! out := A in

  implicit none

  COMPLEX, dimension(18, 2, *)        :: a
  COMPLEX, dimension(NDIRAC, NCOL, *) :: out, in
  integer                             :: volh

  integer :: i
  COMPLEX :: x1, x2, x3, x4, x5, x6
  COMPLEX :: y1, y2, y3, y4, y5, y6

  TIMING_START(timing_bin_clover_mult_a)

  !$omp parallel do private(x1, x2, x3, x4, x5, x6, y1, y2, y3, y4, y5, y6)
  do i = 1, volh
     x1 = in(SC1, i) + in(SC7, i)
     x2 = in(SC2, i) + in(SC8, i)
     x3 = in(SC3, i) + in(SC9, i)
     x4 = in(SC4, i) + in(SC10, i)
     x5 = in(SC5, i) + in(SC11, i)
     x6 = in(SC6, i) + in(SC12, i)

# define J 1
# include "clover_mult_a.h90"

     out(SC1, i) = y1
     out(SC2, i) = y2
     out(SC3, i) = y3
     out(SC4, i) = y4
     out(SC5, i) = y5
     out(SC6, i) = y6
     out(SC7, i) = y1
     out(SC8, i) = y2
     out(SC9, i) = y3
     out(SC10, i) = y4
     out(SC11, i) = y5
     out(SC12, i) = y6

     x1 = in(SC1, i) - in(SC7, i)
     x2 = in(SC2, i) - in(SC8, i)
     x3 = in(SC3, i) - in(SC9, i)
     x4 = in(SC4, i) - in(SC10, i)
     x5 = in(SC5, i) - in(SC11, i)
     x6 = in(SC6, i) - in(SC12, i)

# undef J
# define J 2
# include "clover_mult_a.h90"

     out(SC1, i) = out(SC1, i) + y1
     out(SC2, i) = out(SC2, i) + y2
     out(SC3, i) = out(SC3, i) + y3
     out(SC4, i) = out(SC4, i) + y4
     out(SC5, i) = out(SC5, i) + y5
     out(SC6, i) = out(SC6, i) + y6
     out(SC7, i) = out(SC7, i) - y1
     out(SC8, i) = out(SC8, i) - y2
     out(SC9, i) = out(SC9, i) - y3
     out(SC10, i) = out(SC10, i) - y4
     out(SC11, i) = out(SC11, i) - y5
     out(SC12, i) = out(SC12, i) - y6

  enddo

  TIMING_STOP(timing_bin_clover_mult_a)
end

!===============================================================================
