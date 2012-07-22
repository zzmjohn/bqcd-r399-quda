!===============================================================================
!
! Dxyz.F90 - routines needed in D.F90 (x/y/z-directions)
!
!-------------------------------------------------------------------------------
!
! Copyright (C) 1998-2002 Hinnerk Stueben
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

#ifdef DIR_X
# define GAMMA_A1(C) a(1, C, j) minus i_times(a(4, C, j))
# define GAMMA_A2(C) a(2, C, j) minus i_times(a(3, C, j))
# define GAMMA_B1(C) b(3, C, i) plus  i_times(b2_ ## C)
# define GAMMA_B2(C) b(4, C, i) plus  i_times(b1_ ## C)
#endif

#ifdef DIR_Y
# define GAMMA_A1(C) a(1, C, j) minus a(4, C, j)
# define GAMMA_A2(C) a(2, C, j) plus  a(3, C, j)
# define GAMMA_B1(C) b(3, C, i) plus  b2_ ## C
# define GAMMA_B2(C) b(4, C, i) minus b1_ ## C
#endif

#ifdef DIR_Z
# define GAMMA_A1(C) a(1, C, j) minus i_times(a(3, C, j))
# define GAMMA_A2(C) a(2, C, j) plus  i_times(a(4, C, j))
# define GAMMA_B1(C) b(3, C, i) plus  i_times(b1_ ## C)
# define GAMMA_B2(C) b(4, C, i) minus i_times(b2_ ## C)
#endif

#ifdef FORWARD
# define U(A, B) u(A, B, i)
# define minus MINUS
# define plus PLUS
#else
# define U(A, B) conjg(u(B, A, j))
# define minus PLUS
# define plus MINUS
#endif

#ifdef DAGGER
# define PLUS -
# define MINUS +
#else
# define PLUS +
# define MINUS -
#endif

!-------------------------------------------------------------------------------
subroutine NAME(b, a, u, nn, volh)

  implicit none

  COMPLEX, dimension (NDIRAC, NCOL, *), intent(inout) :: b
  COMPLEX, dimension (NDIRAC, NCOL, *), intent(in) :: a
  COMPLEX, dimension (NCOL, NCOL, *), intent(in) :: u
  INTEGER, dimension (*), intent(in) :: nn
  integer :: volh

  integer :: i, j

  COMPLEX :: a1, a2
  COMPLEX :: b1_1, b2_1
  COMPLEX :: b1_2, b2_2
  COMPLEX :: b1_3, b2_3

  ! statement function:
  
  COMPLEX :: i_times, c
  i_times(c) = cmplx(-aimag(c), real(c), kind = RKIND)

  TIMING_START(STRCAT(timing_bin_, NAME))

  !$omp parallel do private(j, a1, a2, b1_1, b2_1, b1_2, b2_2, b1_3, b2_3)
  do i = 1, volh
     j = nn(i)
 
     a1 = GAMMA_A1(1)
     a2 = GAMMA_A2(1)
     
     b1_1 = a1 * U(1, 1)
     b2_1 = a2 * U(1, 1)
     b1_2 = a1 * U(2, 1)
     b2_2 = a2 * U(2, 1)
     b1_3 = a1 * U(3, 1)
     b2_3 = a2 * U(3, 1)
     
     a1 = GAMMA_A1(2)
     a2 = GAMMA_A2(2)
     
     b1_1 = b1_1 + a1 * U(1, 2)
     b2_1 = b2_1 + a2 * U(1, 2)
     b1_2 = b1_2 + a1 * U(2, 2)
     b2_2 = b2_2 + a2 * U(2, 2)
     b1_3 = b1_3 + a1 * U(3, 2)
     b2_3 = b2_3 + a2 * U(3, 2)
     
     a1 = GAMMA_A1(3)
     a2 = GAMMA_A2(3)
     
     b1_1 = b1_1 + a1 * U(1, 3)
     b2_1 = b2_1 + a2 * U(1, 3)
     
     b(1, 1, i) = b(1, 1, i) + b1_1
     b(2, 1, i) = b(2, 1, i) + b2_1
     b(3, 1, i) = GAMMA_B1(1)
     b(4, 1, i) = GAMMA_B2(1)
     
     b1_2 = b1_2 + a1 * U(2, 3)
     b2_2 = b2_2 + a2 * U(2, 3)

     b(1, 2, i) = b(1, 2, i) + b1_2
     b(2, 2, i) = b(2, 2, i) + b2_2
     b(3, 2, i) = GAMMA_B1(2)
     b(4, 2, i) = GAMMA_B2(2)
     
     b1_3 = b1_3 + a1 * U(3, 3)
     b2_3 = b2_3 + a2 * U(3, 3)

     b(1, 3, i) = b(1, 3, i) + b1_3
     b(2, 3, i) = b(2, 3, i) + b2_3
     b(3, 3, i) = GAMMA_B1(3)
     b(4, 3, i) = GAMMA_B2(3)

  enddo

  TIMING_STOP(STRCAT(timing_bin_, NAME))

end

!===============================================================================
