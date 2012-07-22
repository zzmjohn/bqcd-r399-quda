!===============================================================================
!
! field_io_single_pe.F90 - I/O routine for gauge and pseudo fermion fields
!                          (single processor version)
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
# define INCLUDE_MPIF_H

# define MPI_STATUS_SIZE  2
# define MPI_REAL8        0
# define mpi_real8        0
# define MPI_COMM_WORLD   0
# define MPI_INTEGER8     0
# define mpi_integer8     0
# define MPI_ANY_SOURCE   0

# include "field_io_mpi.F90"

!-------------------------------------------------------------------------------
subroutine mpi_type_vector(a, b, c, d, e, f)
  return
end

!-------------------------------------------------------------------------------
subroutine mpi_type_commit(a, b)
   return
end

!-------------------------------------------------------------------------------
subroutine mpi_type_free(a, b)
   return
end

!-------------------------------------------------------------------------------
subroutine mpi_ssend(a, b, c, d, e, f, g)
   call die("mpi_ssend(): MPI must not be called in single PE version")
end

!-------------------------------------------------------------------------------
subroutine mpi_recv(a, b, c, d, e, f, g, h)
   call die("mpi_recv(): MPI must not be called in single PE version")
end

!===============================================================================
