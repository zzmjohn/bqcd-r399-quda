!===============================================================================
!
! reduction_mpi.F90 - reduction operations in MPI
!
!-------------------------------------------------------------------------------
!
! Copyright (C) 1998-2005 Hinnerk Stueben
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
  include 'mpif.h'
  REAL     global_sum, local_sum
  integer  ierror

  TIMING_START(timing_bin_global_sum)

  call mpi_allreduce(local_sum, global_sum, 1, &
       BQCD_REAL, MPI_SUM, MPI_COMM_WORLD, ierror)

  TIMING_STOP(timing_bin_global_sum)
end

!-------------------------------------------------------------------------------
function global_min(local_min)
 
  implicit none
  include 'mpif.h'
  REAL     global_min, local_min
  integer  ierror

  call mpi_allreduce(local_min, global_min, 1, &
       BQCD_REAL, MPI_MIN, MPI_COMM_WORLD, ierror)
end

!-------------------------------------------------------------------------------
function global_max(local_max)
 
  implicit none
  include 'mpif.h'
  REAL     global_max, local_max
  integer  ierror

  call mpi_allreduce(local_max, global_max, 1, &
       BQCD_REAL, MPI_MAX, MPI_COMM_WORLD, ierror)
end

!-------------------------------------------------------------------------------
subroutine global_sum_vec(n, sum)
 
  implicit none
  include 'mpif.h'
  integer, intent(in)    :: n
  REAL,    intent(inout) :: sum(n)
  REAL                   :: tmp(n)
  integer  ierror

  TIMING_START(timing_bin_global_sum_vec)

  tmp = sum
  call mpi_allreduce(tmp, sum, n, BQCD_REAL, MPI_SUM, MPI_COMM_WORLD, ierror)

  TIMING_STOP(timing_bin_global_sum_vec)
end

!===============================================================================
