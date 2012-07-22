!===============================================================================
!
! Dt.F90 - routines needed in D.F90 (t-direction)
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

#ifdef DAGGER
# define GAMMA_A1(C) a(3, C, jb)
# define GAMMA_A2(C) a(4, C, jb)
# define GAMMA_A3(C) a(1, C, jf)
# define GAMMA_A4(C) a(2, C, jf)
# define GAMMA_B1(C) f1_ ## C
# define GAMMA_B2(C) f2_ ## C
# define GAMMA_B3(C) b1_ ## C
# define GAMMA_B4(C) b2_ ## C
#else
# define GAMMA_A1(C) a(1, C, jb)
# define GAMMA_A2(C) a(2, C, jb)
# define GAMMA_A3(C) a(3, C, jf)
# define GAMMA_A4(C) a(4, C, jf)
# define GAMMA_B1(C) b1_ ## C
# define GAMMA_B2(C) b2_ ## C
# define GAMMA_B3(C) f1_ ## C
# define GAMMA_B4(C) f2_ ## C
#endif

!-------------------------------------------------------------------------------
subroutine NAME(b, a, u_e, u_o, nn_fwd, nn_bwd, volh)

  implicit none

  COMPLEX, dimension (NDIRAC, NCOL, *), intent(inout) :: b
  COMPLEX, dimension (NDIRAC, NCOL, *), intent(in) :: a
  COMPLEX, dimension (NCOL, NCOL, *), intent(in) :: u_e, u_o
  INTEGER, dimension (*), intent(in) :: nn_fwd, nn_bwd
  integer :: volh

  integer :: i, jf, jb

  COMPLEX :: a1, a2, a3, a4
  COMPLEX :: f1_1, f2_1
  COMPLEX :: f1_2, f2_2
  COMPLEX :: f1_3, f2_3
  COMPLEX :: b1_1, b2_1
  COMPLEX :: b1_2, b2_2
  COMPLEX :: b1_3, b2_3

  ! statement function:
  
  COMPLEX :: i_times, c
  i_times(c) = cmplx(-aimag(c), real(c), kind = RKIND)

  TIMING_START(STRCAT(timing_bin_, NAME))

  !$omp parallel do private(jf, jb, a1, a2, a3, a4, &
  !$omp f1_1, f2_1, f1_2, f2_2, f1_3, f2_3, &
  !$omp b1_1, b2_1, b1_2, b2_2, b1_3, b2_3)
  do i = 1, volh
     jb = nn_bwd(i)

     a1 = GAMMA_A1(1)
     a2 = GAMMA_A2(1)

     b1_1 = a1 * conjg(u_o(1, 1, jb))
     b2_1 = a2 * conjg(u_o(1, 1, jb))
     b1_2 = a1 * conjg(u_o(1, 2, jb))
     b2_2 = a2 * conjg(u_o(1, 2, jb))
     b1_3 = a1 * conjg(u_o(1, 3, jb))
     b2_3 = a2 * conjg(u_o(1, 3, jb))
     
     jf = nn_fwd(i)
     
     a3 = GAMMA_A3(1)
     a4 = GAMMA_A4(1)
     
     f1_1 = a3 * u_e(1, 1, i)
     f2_1 = a4 * u_e(1, 1, i)
     f1_2 = a3 * u_e(2, 1, i)
     f2_2 = a4 * u_e(2, 1, i)
     f1_3 = a3 * u_e(3, 1, i)
     f2_3 = a4 * u_e(3, 1, i)
     
     a1 = GAMMA_A1(2)
     a2 = GAMMA_A2(2)
     
     b1_1 = b1_1 + a1 * conjg(u_o(2, 1, jb))
     b2_1 = b2_1 + a2 * conjg(u_o(2, 1, jb))
     b1_2 = b1_2 + a1 * conjg(u_o(2, 2, jb))
     b2_2 = b2_2 + a2 * conjg(u_o(2, 2, jb))
     b1_3 = b1_3 + a1 * conjg(u_o(2, 3, jb))
     b2_3 = b2_3 + a2 * conjg(u_o(2, 3, jb))
     
     a3 = GAMMA_A3(2)
     a4 = GAMMA_A4(2)
     
     f1_1 = f1_1 + a3 * u_e(1, 2, i)
     f2_1 = f2_1 + a4 * u_e(1, 2, i)
     f1_2 = f1_2 + a3 * u_e(2, 2, i)
     f2_2 = f2_2 + a4 * u_e(2, 2, i)
     f1_3 = f1_3 + a3 * u_e(3, 2, i)
     f2_3 = f2_3 + a4 * u_e(3, 2, i)
     
     a1 = GAMMA_A1(3)
     a2 = GAMMA_A2(3)
     
     b1_1 = b1_1 + a1 * conjg(u_o(3, 1, jb))
     b2_1 = b2_1 + a2 * conjg(u_o(3, 1, jb))
     b1_2 = b1_2 + a1 * conjg(u_o(3, 2, jb))
     b2_2 = b2_2 + a2 * conjg(u_o(3, 2, jb))
     b1_3 = b1_3 + a1 * conjg(u_o(3, 3, jb))
     b2_3 = b2_3 + a2 * conjg(u_o(3, 3, jb))
     
     a3 = GAMMA_A3(3)
     a4 = GAMMA_A4(3)
     
     f1_1 = f1_1 + a3 * u_e(1, 3, i)
     f2_1 = f2_1 + a4 * u_e(1, 3, i)
     f1_2 = f1_2 + a3 * u_e(2, 3, i)
     f2_2 = f2_2 + a4 * u_e(2, 3, i)
     f1_3 = f1_3 + a3 * u_e(3, 3, i)
     f2_3 = f2_3 + a4 * u_e(3, 3, i)
     
     
     b(1, 1, i) = TWO * GAMMA_B1(1)
     b(2, 1, i) = TWO * GAMMA_B2(1)
     b(3, 1, i) = TWO * GAMMA_B3(1)
     b(4, 1, i) = TWO * GAMMA_B4(1)
     
     b(1, 2, i) = TWO * GAMMA_B1(2)
     b(2, 2, i) = TWO * GAMMA_B2(2)
     b(3, 2, i) = TWO * GAMMA_B3(2)
     b(4, 2, i) = TWO * GAMMA_B4(2)
     
     b(1, 3, i) = TWO * GAMMA_B1(3)
     b(2, 3, i) = TWO * GAMMA_B2(3)
     b(3, 3, i) = TWO * GAMMA_B3(3)
     b(4, 3, i) = TWO * GAMMA_B4(3)
     
  enddo

  TIMING_STOP(STRCAT(timing_bin_, NAME))

end

!===============================================================================
