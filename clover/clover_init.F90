!===============================================================================
!
! clover_init.F90 - calculates clover matrix and its inverse
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
subroutine clover_init(a, ainv, b, u, csw_kappa)

  use typedef_clover
  use module_vol
  implicit none

  CLOVER_FIELD_A, intent(out) :: a, ainv
  CLOVER_FIELD_B, intent(out) :: b
  GAUGE_FIELD,    intent(in)  :: u
  REAL,           intent(in)  :: csw_kappa

  integer             :: i, eo
  SU3                 :: f, g
  type(type_clover_a) :: p, q
  REAL                :: factor

  TIMING_START(timing_bin_clover_init)

  factor = -csw_kappa / EIGHT

  do eo = EVEN, ODD
     !$omp parallel do private(f, g, p, q)
     do i = 1, VOLH
        call clover_f_mu_nu(f, 2, 1, i, eo, u)

        call clover_init1(p, f)

        call clover_f_mu_nu(f, 3, 2, i, eo, u)
        call clover_f_mu_nu(g, 3, 1, i, eo, u)

        call clover_init2(p, f, g)

        call clover_f_mu_nu(f, 3, 4, i, eo, u)

        call clover_init1(q, f)

        call clover_f_mu_nu(f, 1, 4, i, eo, u)
        call clover_f_mu_nu(g, 4, 2, i, eo, u)

        call clover_init2(q, f, g)

        call clover_init3(a(1, i, eo), a(2, i, eo), p, q, factor)

        call clover_inv(b(1, i, eo), ainv(1, i, eo), a(1, i, eo))
        call clover_inv(b(2, i, eo), ainv(2, i, eo), a(2, i, eo))
     enddo
  enddo

  TIMING_STOP(timing_bin_clover_init)

end

!-------------------------------------------------------------------------------
subroutine clover_init1(a, f)

  use typedef_clover
  implicit none
  type(type_clover_a) :: a
  SU3                 :: f

  a%i11 = Re(f(1, 1))
  a%i22 = Re(f(2, 2))
  a%i33 = Re(f(3, 3))
  
  a%i44 = -a%i11
  a%i55 = -a%i22
  a%i66 = -a%i33
  
  a%i12 = f(1, 2)
  a%i13 = f(1, 3)
  a%i23 = f(2, 3)
  
  a%i45 = -a%i12
  a%i46 = -a%i13
  a%i56 = -a%i23
  
end

!-------------------------------------------------------------------------------
subroutine clover_init2(a, f, g)

  use typedef_clover
  implicit none
  type(type_clover_a) :: a
  SU3                 :: f, g

  ! statement function:
  
  COMPLEX :: i_times, c
  i_times(c) = cmplx(-aimag(c), real(c), kind = RKIND)

  a%i14 = f(1, 1) + i_times(g(1, 1))
  a%i15 = f(1, 2) + i_times(g(1, 2))
  a%i16 = f(1, 3) + i_times(g(1, 3))
  
  a%i24 = f(2, 1) + i_times(g(2, 1))
  a%i25 = f(2, 2) + i_times(g(2, 2))
  a%i26 = f(2, 3) + i_times(g(2, 3))
  
  a%i34 = f(3, 1) + i_times(g(3, 1))
  a%i35 = f(3, 2) + i_times(g(3, 2))
  a%i36 = f(3, 3) + i_times(g(3, 3))

end

!-------------------------------------------------------------------------------
subroutine clover_init3(a1, a2, p, q, s)

  use typedef_clover
  implicit none
  type(type_clover_a) :: a1, a2, p, q
  REAL                :: s

# define CLOVER_INIT_3(I, J) \
a1%i ## I ## J = s * (p%i ## I ## J + q%i ## I ## J) ; \
a2%i ## I ## J = s * (p%i ## I ## J - q%i ## I ## J)

! define => 
! a1%iIJ = s * (p%iIJ + q%iIJ) ; a2%iIJ = s * (p%iIJ - q%iIJ)

  CLOVER_INIT_3(1, 1)
  CLOVER_INIT_3(1, 2)
  CLOVER_INIT_3(1, 3)
  CLOVER_INIT_3(1, 4)
  CLOVER_INIT_3(1, 5)
  CLOVER_INIT_3(1, 6)

  CLOVER_INIT_3(2, 2)
  CLOVER_INIT_3(2, 3)
  CLOVER_INIT_3(2, 4)
  CLOVER_INIT_3(2, 5)
  CLOVER_INIT_3(2, 6)

  CLOVER_INIT_3(3, 3)
  CLOVER_INIT_3(3, 4)
  CLOVER_INIT_3(3, 5)
  CLOVER_INIT_3(3, 6)

  CLOVER_INIT_3(4, 4)
  CLOVER_INIT_3(4, 5)
  CLOVER_INIT_3(4, 6)

  CLOVER_INIT_3(5, 5)
  CLOVER_INIT_3(5, 6)

  CLOVER_INIT_3(6, 6)

  a1%i11 = a1%i11 + ONE
  a1%i22 = a1%i22 + ONE
  a1%i33 = a1%i33 + ONE
  a1%i44 = a1%i44 + ONE
  a1%i55 = a1%i55 + ONE
  a1%i66 = a1%i66 + ONE

  a2%i11 = a2%i11 + ONE
  a2%i22 = a2%i22 + ONE
  a2%i33 = a2%i33 + ONE
  a2%i44 = a2%i44 + ONE
  a2%i55 = a2%i55 + ONE
  a2%i66 = a2%i66 + ONE

end

!===============================================================================
