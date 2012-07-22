!===============================================================================
!
! seed_mpi.F90
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
subroutine seed_broadcast(seed)

  use module_function_decl
  implicit none
  include 'mpif.h'
  SEED     seed
  integer  ierror

  call mpi_bcast(seed, 1, BQCD_SEED, 0, MPI_COMM_WORLD, ierror)
end  

!-------------------------------------------------------------------------------
subroutine seed_compare(seed)

  use module_function_decl
  implicit none
  include 'mpif.h'  
  SEED     seed, s
  integer  pe, status(MPI_STATUS_SIZE), ierror

  if (my_pe() /= 0) then
     call mpi_ssend(seed, 1, BQCD_SEED, 0, 0, MPI_COMM_WORLD, ierror)
  else
     do pe = 1, num_pes() - 1
        call mpi_recv(s, 1, BQCD_SEED, pe, 0, MPI_COMM_WORLD, status, ierror)
        if (s /= seed) call die('rancheck(): seeds differ')
     enddo
  endif     

end

!===============================================================================
