!===============================================================================
!
! seed_shmem.F90
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
# include "shmem.h"

!-------------------------------------------------------------------------------
subroutine seed_broadcast(seed)

  use module_function_decl
  implicit none
  include "mpp/shmem.fh"

  SEED          :: seed
  SEED,    save :: s
  integer, save :: psync(SHMEM_BCAST_SYNC_SIZE)

  psync = SHMEM_SYNC_VALUE
  s = seed

  call barrier()
  call shmem_broadcast(s, s, 1, 0, 0, 0, num_pes(), psync)
  call barrier()

  seed = s
end

!-------------------------------------------------------------------------------
subroutine seed_compare(seed)

  use module_function_decl
  implicit none

  SEED       :: seed
  SEED, save :: s
  integer    :: pe

  s = seed
  call barrier()

  if (my_pe() == 0) then
     do pe = 1, num_pes() - 1
        call shmem_get(s, s, 1, pe)
        if (s /= seed) call die('rancheck(): seeds differ')
     enddo
  endif     

  call barrier()
end

!===============================================================================
