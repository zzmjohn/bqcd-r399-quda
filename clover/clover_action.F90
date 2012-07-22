!===============================================================================
!
! clover_action.F90  -  calculates: -2 Tr(log(T_oo))
!
!-------------------------------------------------------------------------------
!
! Copyright (C) 1998-2003 Hinnerk Stueben
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
# include "defs.h"
# include "clover.h"

!-------------------------------------------------------------------------------
REAL function clover_action(b)

  use typedef_clover
  use module_vol
  implicit none

  type(type_clover_b) :: b(2, volh)
  integer             :: i
  REAL                :: s, global_sum
  

  s = ZERO

  !$omp parallel do reduction(+: s)
  do i = 1, volh
     s = s + log(det(b(1, i)) * det(b(2, i)))
  enddo

  clover_action = TWO * global_sum(s)


CONTAINS

  REAL function det(b)  ! returns (1 / det)

    type(type_clover_b) :: b

    det = B11 * B22 * B33 * B44 * B55 * B66

  end function det

end

!===============================================================================
