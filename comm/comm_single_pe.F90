!===============================================================================
!
! comm_single_pe.F90 - (dummy) routines for single CPU version
!
!-------------------------------------------------------------------------------
!
! Copyright (C) 1998-2012 Hinnerk Stueben
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
subroutine comm_init()
  return
end

!-------------------------------------------------------------------------------
subroutine comm_finalize()
  return
end

!-------------------------------------------------------------------------------
subroutine comm_abort()

  call exit(i)
end

!-------------------------------------------------------------------------------
COMM_METHOD function comm_method()

#ifdef _OPENMP
  comm_method = "single_pe + OpenMP"
#else
  comm_method = "single_pe"
#endif
end

!===============================================================================
