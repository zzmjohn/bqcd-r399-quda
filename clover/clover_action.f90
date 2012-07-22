!===============================================================================
!
! clover_action.F90 - calculates: -2 Tr(log(T_oo))
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
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with BQCD. If not, see <http://www.gnu.org/licenses/>.
!
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
real(8) function clover_action(b)

  use typedef_clover
  use module_vol
  implicit none

  type(type_clover_b) :: b(2, volh)
  integer :: i
  real(8) :: s, global_sum


  s = 0.0_8

  !$omp parallel do reduction(+: s)
  do i = 1, volh
     s = s + log(det(b(1, i)) * det(b(2, i)))
  enddo

  clover_action = 2.0_8 * global_sum(s)


CONTAINS

  real(8) function det(b) ! returns (1 / det)

    type(type_clover_b) :: b

    det = b%i11 * b%i22 * b%i33 * b%i44 * b%i55 * b%i66

  end function det

end

!===============================================================================
