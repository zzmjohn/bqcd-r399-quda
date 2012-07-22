!===============================================================================
!
! clover_dummy.F90
!
!-------------------------------------------------------------------------------
!
! Copyright (C) 2000-2003 Hinnerk Stueben
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
subroutine allocate_clover_field_a(a)

  use typedef_clover
  use module_vol
  implicit none
  P_CLOVER_FIELD_A :: a

  call die("allocate_clover_field_a(): must not be called.")
end

!-------------------------------------------------------------------------------
subroutine allocate_clover_field_b(b)

  use typedef_clover
  use module_vol
  implicit none
  P_CLOVER_FIELD_B :: b

  call die("allocate_clover_field_b(): must not be called.")
end

!-------------------------------------------------------------------------------
subroutine clover_init(a, b, u, csw_kappa)

  use typedef_clover
  use module_vol
  implicit none

  CLOVER_FIELD_A, intent(out) :: a
  CLOVER_FIELD_B, intent(out) :: b
  GAUGE_FIELD,    intent(in)  :: u
  REAL,           intent(in)  :: csw_kappa

  call die("clover_init(): must not be called.")
end

!-------------------------------------------------------------------------------
REAL function clover_action(b)

  use typedef_clover
  use module_vol
  implicit none

  type(type_clover_b) :: b(2, volh)
  integer             :: i
  REAL                :: s, global_sum
  

  call die("clover_action(): must not be called.")
  clover_action = ZERO
end

!-------------------------------------------------------------------------------
subroutine clover_mult_a(out, a, in, volh)

  implicit none

  COMPLEX, dimension(18, 2, *)        :: a
  COMPLEX, dimension(NDIRAC, NCOL, *) :: out, in
  integer                             :: volh

  call die("clover_mult_a(): must not be called.")
end

!-------------------------------------------------------------------------------
subroutine clover_mult_ao(a, x, volh)  ! x := A x

  implicit none

  COMPLEX, dimension(18, 2, *)        :: a
  COMPLEX, dimension(NDIRAC, NCOL, *) :: x
  integer                             :: volh

  call die("clover_mult_ao(): must not be called.")
end

!-------------------------------------------------------------------------------
subroutine clover_mult_b(b, x, volh)

  implicit none

  COMPLEX, dimension(18, 2, *)        :: b
  COMPLEX, dimension(NDIRAC, NCOL, *) :: x
  integer                             :: volh

  call die("clover_mult_b(): must not be called.")
end

!-------------------------------------------------------------------------------
subroutine clover_dsd(eo, p, b, s, u)

  use typedef_clover
  use module_vol
  implicit none

  integer         :: eo
  GENERATOR_FIELD :: p
  CLOVER_FIELD_B  :: b
  REAL            :: s
  GAUGE_FIELD     :: u

  call die("clover_dsd(): must not be called.")
end

!-------------------------------------------------------------------------------
subroutine clover_dsf(eo, p, b, a, s, u)

  use module_vol
  implicit none

  integer         :: eo
  GENERATOR_FIELD :: p
  SPINCOL_FIELD   :: b, a
  REAL            :: s
  GAUGE_FIELD     :: u

  call die("clover_dsf(): must not be called.")
end

!===============================================================================
