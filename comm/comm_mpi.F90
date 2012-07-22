!===============================================================================
!
! comm_mpi.F90 - wrapper for MPI routines
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

  implicit none
  include 'mpif.h'
  integer  ierror

  call mpi_init(ierror)
end

!-------------------------------------------------------------------------------
subroutine comm_finalize()

  implicit none
  include 'mpif.h'
  integer  ierror

  call mpi_finalize(ierror)
end

!-------------------------------------------------------------------------------
subroutine comm_abort()

  implicit none
  include 'mpif.h'
  integer ierror

  call mpi_abort(MPI_COMM_WORLD, 1, ierror)
end

!-------------------------------------------------------------------------------
COMM_METHOD function comm_method()

#ifdef _OPENMP
  comm_method = "MPI + OpenMP"
#else
  comm_method = "MPI"
#endif
end

!===============================================================================
