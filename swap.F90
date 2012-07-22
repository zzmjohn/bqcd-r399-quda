!===============================================================================
!
! swap.F90 - swap routines for various data types
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

!-------------------------------------------------------------------------------
subroutine swap_p_g_field(u, v)

  implicit none
  P_GAUGE_FIELD :: u, v, tmp

  tmp => u
  u => v
  v => tmp

end

!-------------------------------------------------------------------------------
subroutine swap_p_sc_field(a, b)

  implicit none
  P_SPINCOL_FIELD :: a, b, tmp

  tmp => a
  a => b
  b => tmp

end

!-------------------------------------------------------------------------------
subroutine swap_p_clover_field_a(x, y)

  use typedef_clover
  implicit none
  P_CLOVER_FIELD_A :: x, y, tmp

  tmp => x
  x => y
  y => tmp

end

!-------------------------------------------------------------------------------
subroutine swap_p_clover_field_b(x, y)

  use typedef_clover
  implicit none
  P_CLOVER_FIELD_B :: x, y, tmp

  tmp => x
  x => y
  y => tmp

end

!-------------------------------------------------------------------------------
subroutine swap_real(x, y)

  implicit none
  REAL :: x, y, tmp
  
  tmp = x
  x = y
  y = tmp

end

!-------------------------------------------------------------------------------
subroutine swap_integer(x, y)

  implicit none
  integer :: x, y, tmp
  
  tmp = x
  x = y
  y = tmp

end

!===============================================================================
