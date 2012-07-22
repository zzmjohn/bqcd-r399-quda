!===============================================================================
!
! pes_mpi.F90 - MPI version of shmem functions
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
integer function my_pe()

  implicit none
  include 'mpif.h'
  integer  ierror

  call mpi_comm_rank(MPI_COMM_WORLD, my_pe, ierror)
end

!-------------------------------------------------------------------------------
integer function num_pes()

  implicit none
  include 'mpif.h'
  integer  ierror

  call mpi_comm_size(MPI_COMM_WORLD, num_pes, ierror)
end

!===============================================================================
