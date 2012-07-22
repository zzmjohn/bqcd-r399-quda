!===============================================================================
!
! reduction_shmem.F90 - reduction operations in shmem
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
function global_sum(local_sum)
 
  implicit none
  include 'mpp/shmem.fh'

  REAL              :: global_sum, local_sum
  REAL, save        :: source, target
  integer           :: n_pes
  integer, external :: num_pes   
  REAL, save        :: pWrk(2 + shmem_reduce_min_wrkdata_size)
  integer, save     :: pSync(shmem_reduce_sync_size)


  TIMING_START(timing_bin_global_sum)

  n_pes = num_pes()

  if (n_pes == 1) then
     global_sum = local_sum
     return
  endif

  source = local_sum

  call shmem_real8_sum_to_all(target, source, 1, 0, 0, n_pes, pWrk, pSync)

  global_sum = target

  TIMING_STOP(timing_bin_global_sum)
end

!-------------------------------------------------------------------------------
function global_min(local_min)
 
  implicit none
  include 'mpp/shmem.fh'

  REAL              :: global_min, local_min
  REAL, save        :: source, target
  integer           :: n_pes
  integer, external :: num_pes   
  REAL, save        :: pWrk(2 + shmem_reduce_min_wrkdata_size)
  integer, save     :: pSync(shmem_reduce_sync_size)

  n_pes = num_pes()

  if (n_pes == 1) then
     global_min = local_min
     return
  endif

  source = local_min

  call shmem_real8_min_to_all(target, source, 1, 0, 0, n_pes, pWrk, pSync)

  global_min = target
end

!-------------------------------------------------------------------------------
function global_max(local_max)
 
  implicit none
  include 'mpp/shmem.fh'

  REAL              :: global_max, local_max
  REAL, save        :: source, target
  integer           :: n_pes
  integer, external :: num_pes   
  REAL, save        :: pWrk(2 + shmem_reduce_min_wrkdata_size)
  integer, save     :: pSync(shmem_reduce_sync_size)

  n_pes = num_pes()

  if (n_pes == 1) then
     global_max = local_max
     return
  endif

  source = local_max

  call shmem_real8_max_to_all(target, source, 1, 0, 0, n_pes, pWrk, pSync)

  global_max = target
end

!-------------------------------------------------------------------------------
subroutine global_sum_vec(n, sum)
 
  implicit none

  integer, intent(in)    :: n
  REAL,    intent(inout) :: sum(n)

  TIMING_START(timing_bin_global_sum_vec)

  call die("global_sum_vec(): shmem version not implemented yet.")

  TIMING_STOP(timing_bin_global_sum_vec)
end

!===============================================================================
