!===============================================================================
!
! xbound_shmem.F90 - boundary exchange with shmem
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
# include "shmem.h"

!-------------------------------------------------------------------------------
subroutine init_xbound()

  return
end

!-------------------------------------------------------------------------------
subroutine xbound_g(u, eo, mu)

  use module_function_decl
  use module_vol
  implicit none

  integer :: eo, mu, x, y, z, t
  GAUGE_FIELD :: u

  if (num_pes() == 1) return

  call barrier()

  do t = -1, 1
  do z = -1, 1
  do y = -1, 1
  do x = -1, 1
     call xch_bound(NCOL * NCOL * SIZE_COMPLEX, u(1, 1, 1, eo, mu), x, y, z, t)
     call barrier()
  enddo
  enddo
  enddo
  enddo

end

!-------------------------------------------------------------------------------
subroutine xbound_g_field(u)

  use module_function_decl
  use module_vol
  implicit none

  GAUGE_FIELD :: u
  integer :: mu, eo, x, y, z, t

  if (num_pes() == 1) return

  call barrier()

  do mu = 1, DIM
  do eo = EVEN, ODD
     call xbound_g(u, eo, mu)
  enddo
  enddo

end

!-------------------------------------------------------------------------------
subroutine xbound_sc_field(array)

  use module_function_decl
  use module_vol
  implicit none

  SPINCOL_FIELD :: array
  integer :: x, y, z, t

  if (num_pes() == 1) return

  call barrier()

  do t = -1, 1
  do z = -1, 1
  do y = -1, 1
  do x = -1, 1
     if ((abs(x) + abs(y) + abs(z) + abs(t)) == 1) then
        call xch_bound(NDIRAC * NCOL * SIZE_COMPLEX, array, x, y, z, t)
        call barrier()
     endif
  enddo
  enddo
  enddo
  enddo
end

!-------------------------------------------------------------------------------
subroutine xbound_sc2_field_i(a)

  use module_function_decl
  use module_lattice
  use module_vol
  implicit none

  SC2_FIELD :: a
  integer :: x, y, z, t

  if (num_pes() == 1) return

  call barrier()

  call xch_bound(2 * NCOL * SIZE_COMPLEX, a(1, 1, 1, 1, FWD), +1,  0,  0,  0)
  call xch_bound(2 * NCOL * SIZE_COMPLEX, a(1, 1, 1, 1, BWD), -1,  0,  0,  0)
  call xch_bound(2 * NCOL * SIZE_COMPLEX, a(1, 1, 1, 2, FWD),  0, +1,  0,  0)
  call xch_bound(2 * NCOL * SIZE_COMPLEX, a(1, 1, 1, 2, BWD),  0, -1,  0,  0)
  call xch_bound(2 * NCOL * SIZE_COMPLEX, a(1, 1, 1, 3, FWD),  0,  0, +1,  0)
  call xch_bound(2 * NCOL * SIZE_COMPLEX, a(1, 1, 1, 3, BWD),  0,  0, -1,  0)
  call xch_bound(2 * NCOL * SIZE_COMPLEX, a(1, 1, 1, 4, FWD),  0,  0,  0, +1)
  call xch_bound(2 * NCOL * SIZE_COMPLEX, a(1, 1, 1, 4, BWD),  0,  0,  0, -1)

  call barrier()

end

!-------------------------------------------------------------------------------
subroutine xch_bound(mm, array, xx, yy, zz, tt)

  use module_function_decl
  use module_nnpe
  use module_offset
  use module_lattice
  use module_vol
  implicit none
  include 'mpif.h'

  integer :: mm
  REAL, dimension (mm, volh_tot) :: array
  integer, dimension (DIM) :: dir, m, i, target, source
  integer, external :: xyzt2i
  integer :: xx, yy, zz, tt, x, y, z, t, pe, size, mu

  pe = nnpe(xx, yy, zz, tt)
  if (pe == my_pe()) return

  dir = (/ xx, yy, zz, tt /)

  do mu = 1, DIM
     if (dir(mu) /= 0) then
        m(mu) = 1
     else
        m(mu) = NH(mu)
     endif
  enddo

  size = mm
  if (dir(1) == 0) then
     size = size * NXH
     m(1) = 1
     if (dir(2) == 0) then
        size = size * N(2)
        m(2) = 1
        if (dir(3) == 0) then
           size = size * N(3)
           m(3) = 1
           if (dir(4) == 0) then
              size = size * N(4)
              m(4) = 1
           endif
        endif
     endif
  endif

  do t = 0, m(4) - 1
  do z = 0, m(3) - 1
  do y = 0, m(2) - 1
  do x = 0, m(1) - 1

     i = (/ x, y, z, t /)

     do mu = 1, DIM
        if (dir(mu) == -1) then
           target(mu) = -1
           source(mu) = N(mu) - 1
        elseif (dir(mu) == +1) then
           target(mu) = N(mu)
           source(mu) = 0
        else
           target(mu) = i(mu)
           source(mu) = i(mu)
        endif
     enddo

!!!     call shmem_get(array(1, xyzt2i(target)), &
!!!                    array(1, xyzt2i(source)), size, pe)
     call shmem_put(array(1, xyzt2i(target)), &
                    array(1, xyzt2i(source)), size, nnpe(-xx,-yy,-zz,-tt))

  enddo
  enddo
  enddo
  enddo

end

!===============================================================================
