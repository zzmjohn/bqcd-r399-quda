!===============================================================================
!
! comm_broadcast_single_pe.F90 - stubs for broadcast routines
!
!-------------------------------------------------------------------------------
!
! Copyright (C) 2008 Hinnerk Stueben
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
subroutine comm_broadcast_i4(i, n)

  implicit none
  integer, intent(in) :: n
  integer(4), intent(inout) :: i(n)

  return
end

!-------------------------------------------------------------------------------
subroutine comm_broadcast_i8(i, n)

  implicit none
  integer, intent(in) :: n
  integer(8), intent(inout) :: i(n)

  return
end

!==============================================================================
