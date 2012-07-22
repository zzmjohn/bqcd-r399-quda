!===============================================================================
!
! allocate.F90 - allocation of gauge and pseudo fermion fields
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
subroutine allocate_g_field(u)

  use module_vol
  implicit none
  P_GAUGE_FIELD :: u

  if (associated(u)) then
     call die("allocate_g_field(): memory leak")
  else
     allocate(u(NCOL, NCOL, volh_tot, EVEN:ODD, DIM))
     call conf_zero(u)
  endif
end

!-------------------------------------------------------------------------------
subroutine allocate_g_field_io(u)

  use module_lattice_io
  implicit none
  P_GAUGE_FIELD_IO :: u

  if (associated(u)) then
     call die("allocate_g_field_io(): memory leak")
  else
     allocate(u(NCOL, NCOL-1, DIM, 0:NX-1, 0:NY-1, 0:NZ-1, 0:NT-1))
  endif
end

!-------------------------------------------------------------------------------
subroutine allocate_gen_field(x)

  use module_vol
  implicit none
  P_GENERATOR_FIELD :: x

  integer :: i, eo, mu

  if (associated(x)) then
     call die("allocate_gen_field(): memory leak")
  else
     allocate(x(NGEN, volh_tot, EVEN:ODD, DIM))
     do mu = 1, DIM
        do eo = EVEN, ODD
           !$omp parallel do
           do i = 1, volh
              x(1, i, eo, mu) = ZERO
              x(2, i, eo, mu) = ZERO
              x(3, i, eo, mu) = ZERO
              x(4, i, eo, mu) = ZERO
              x(5, i, eo, mu) = ZERO
              x(6, i, eo, mu) = ZERO
              x(7, i, eo, mu) = ZERO
              x(8, i, eo, mu) = ZERO
           enddo
        enddo
     enddo
  endif
end

!-------------------------------------------------------------------------------
subroutine allocate_sc_field(x)

  use module_vol
  implicit none
  P_SPINCOL_FIELD :: x

  if (associated(x)) then
     call die("allocate_sc_field(): memory leak")
  else
     allocate(x(NDIRAC, NCOL, volh_tot))
     call sc_zero(x)
  endif
end

!-------------------------------------------------------------------------------
subroutine allocate_sc_field_io(x)

  use module_lattice_io
  implicit none
  P_SPINCOL_FIELD_IO :: x

  if (associated(x)) then
     call die("allocate_sc_field_io(): memory leak")
  else
     allocate(x(NDIRAC, NCOL, 0:NXH-1, 0:NY-1, 0:NZ-1, 0:NT-1))
  endif
end

!-------------------------------------------------------------------------------
subroutine allocate_sc_overindexed(x)

  use module_vol
  implicit none
  P_SPINCOL_OVERINDEXED :: x

  if (associated(x)) then
     call die("allocate_sc_overindexed(): memory leak")
  else
     allocate(x(SIZE_COMPLEX*NDIRAC*NCOL*volh_tot))
  endif
end

!-------------------------------------------------------------------------------
subroutine allocate_sc2_field(x)

  use module_vol
  implicit none
  P_SC2_FIELD :: x

  if (associated(x)) then
     call die("allocate_sc2_field(): memory leak")
  else
     allocate(x(2, NCOL, volh_tot, DIM, FWD:BWD))
  endif
end

!===============================================================================
