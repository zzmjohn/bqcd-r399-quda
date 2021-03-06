!===============================================================================
!
! cksum_dummy.F90 - dummy routines that can replace the real C routines
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
# include "defs.h"

!-------------------------------------------------------------------------------
subroutine cksum_init()
  return
end

!-------------------------------------------------------------------------------
subroutine cksum_add(i, j)
  integer i(*)
  CHECK_SUM j
  return
end

!-------------------------------------------------------------------------------
subroutine cksum_get(sum, bytes)
  CHECK_SUM sum, bytes
  sum = 0
  bytes = 0
end

!===============================================================================
