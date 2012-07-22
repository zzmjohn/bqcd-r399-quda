!===============================================================================
!
! dotprod.F90 - dot product for parallel computers
!
!-------------------------------------------------------------------------------
!
! Copyright (C) 1998-2006 Hinnerk Stueben
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
REAL function dotprod(a, b, n)
 
  implicit none
  integer  i, n
  REAL     a(n), b(n), s, global_sum
  
  s = ZERO
  !$omp parallel do reduction(+: s)
  do i = 1, n
     s = s + a(i) * b(i)
  enddo

  dotprod = global_sum(s)
 
end

!===============================================================================
