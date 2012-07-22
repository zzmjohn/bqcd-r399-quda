!===============================================================================
!
! clover_mult_ao.F90 - ao: "A overwrite"
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
subroutine clover_mult_ao(a, x, volh)  ! x := A x

  implicit none

  COMPLEX, dimension(18, 2, *)        :: a
  COMPLEX, dimension(NDIRAC, NCOL, *) :: x
  integer                             :: volh

  integer :: i
  COMPLEX :: x1, x2, x3, x4, x5, x6
  COMPLEX :: y1, y2, y3, y4, y5, y6

  TIMING_START(timing_bin_clover_mult_ao)

  !$omp parallel do private(x1, x2, x3, x4, x5, x6, y1, y2, y3, y4, y5, y6)
  do i = 1, volh
     x1 = x(SC1, i) + x(SC7, i)
     x2 = x(SC2, i) + x(SC8, i)
     x3 = x(SC3, i) + x(SC9, i)
     x4 = x(SC4, i) + x(SC10, i)
     x5 = x(SC5, i) + x(SC11, i)
     x6 = x(SC6, i) + x(SC12, i)

# define J 1
# include "clover_mult_a.h90"

     x1 = x(SC1, i) - x(SC7, i)
     x2 = x(SC2, i) - x(SC8, i)
     x3 = x(SC3, i) - x(SC9, i)
     x4 = x(SC4, i) - x(SC10, i)
     x5 = x(SC5, i) - x(SC11, i)
     x6 = x(SC6, i) - x(SC12, i)

     x(SC1, i) = y1
     x(SC2, i) = y2
     x(SC3, i) = y3
     x(SC4, i) = y4
     x(SC5, i) = y5
     x(SC6, i) = y6
     x(SC7, i) = y1
     x(SC8, i) = y2
     x(SC9, i) = y3
     x(SC10, i) = y4
     x(SC11, i) = y5
     x(SC12, i) = y6

# undef J
# define J 2
# include "clover_mult_a.h90"

     x(SC1, i) = x(SC1, i) + y1
     x(SC2, i) = x(SC2, i) + y2
     x(SC3, i) = x(SC3, i) + y3
     x(SC4, i) = x(SC4, i) + y4
     x(SC5, i) = x(SC5, i) + y5
     x(SC6, i) = x(SC6, i) + y6
     x(SC7, i) = x(SC7, i) - y1
     x(SC8, i) = x(SC8, i) - y2
     x(SC9, i) = x(SC9, i) - y3
     x(SC10, i) = x(SC10, i) - y4
     x(SC11, i) = x(SC11, i) - y5
     x(SC12, i) = x(SC12, i) - y6

  enddo

  TIMING_STOP(timing_bin_clover_mult_ao)
end

!===============================================================================
