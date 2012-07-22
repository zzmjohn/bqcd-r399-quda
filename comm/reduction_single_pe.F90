!===============================================================================
!
! reduction_single_pe.F90 - reduction operations on a single processor
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
function global_sum(local_sum)
 
  implicit none
  REAL :: global_sum, local_sum

  global_sum = local_sum
end

!-------------------------------------------------------------------------------
function global_min(local_min)
 
  implicit none
  REAL :: global_min, local_min

  global_min = local_min
end

!-------------------------------------------------------------------------------
function global_max(local_max)
 
  implicit none
  REAL :: global_max, local_max

  global_max = local_max
end

!-------------------------------------------------------------------------------
subroutine global_sum_vec(n, sum)
 
  implicit none
  integer, intent(in)    :: n
  REAL,    intent(inout) :: sum(n)

  return
end
!===============================================================================
