!===============================================================================
!
! D2xyzt.F90 - routines needed in D2.F90
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
# include "defs.h"

#ifdef DIR_T

#ifdef DAGGER
# define GAMMA_AB1(C) a(3, C, jb)
# define GAMMA_AB2(C) a(4, C, jb)
# define GAMMA_AF1(C) a(1, C, jf)
# define GAMMA_AF2(C) a(2, C, jf)
# define GAMMA_B1(C) bf1_ ## C
# define GAMMA_B2(C) bf2_ ## C
# define GAMMA_B3(C) bb1_ ## C
# define GAMMA_B4(C) bb2_ ## C
#else
# define GAMMA_AB1(C) a(1, C, jb)
# define GAMMA_AB2(C) a(2, C, jb)
# define GAMMA_AF1(C) a(3, C, jf)
# define GAMMA_AF2(C) a(4, C, jf)
# define GAMMA_B1(C) bb1_ ## C
# define GAMMA_B2(C) bb2_ ## C
# define GAMMA_B3(C) bf1_ ## C
# define GAMMA_B4(C) bf2_ ## C
#endif

# define UPDATE_B(S, C) b(S, C, i) = TWO * GAMMA_B ## S(C)

#else

#ifdef DAGGER
# define PLUS -
# define MINUS +
#else
# define PLUS +
# define MINUS -
#endif

#ifdef DIR_X
# define GAMMA_AB1(C) a(1, C, jb) PLUS  i_times(a(4, C, jb))
# define GAMMA_AB2(C) a(2, C, jb) PLUS  i_times(a(3, C, jb))
# define GAMMA_AF1(C) a(1, C, jf) MINUS i_times(a(4, C, jf))
# define GAMMA_AF2(C) a(2, C, jf) MINUS i_times(a(3, C, jf))
# define GAMMA_B3(C) MINUS i_times(bb2_ ## C) PLUS  i_times(bf2_ ## C)
# define GAMMA_B4(C) MINUS i_times(bb1_ ## C) PLUS  i_times(bf1_ ## C)
#endif

#ifdef DIR_Y
# define GAMMA_AB1(C) a(1, C, jb) PLUS  a(4, C, jb)
# define GAMMA_AB2(C) a(2, C, jb) MINUS a(3, C, jb)
# define GAMMA_AF1(C) a(1, C, jf) MINUS a(4, C, jf)
# define GAMMA_AF2(C) a(2, C, jf) PLUS  a(3, C, jf)
# define GAMMA_B3(C) MINUS bb2_ ## C PLUS  bf2_ ## C
# define GAMMA_B4(C) PLUS  bb1_ ## C MINUS bf1_ ## C
#endif

#ifdef DIR_Z
# define GAMMA_AB1(C) a(1, C, jb) PLUS  i_times(a(3, C, jb))
# define GAMMA_AB2(C) a(2, C, jb) MINUS i_times(a(4, C, jb))
# define GAMMA_AF1(C) a(1, C, jf) MINUS i_times(a(3, C, jf))
# define GAMMA_AF2(C) a(2, C, jf) PLUS  i_times(a(4, C, jf))
# define GAMMA_B3(C) MINUS i_times(bb1_ ## C) PLUS  i_times(bf1_ ## C)
# define GAMMA_B4(C) PLUS  i_times(bb2_ ## C) MINUS i_times(bf2_ ## C)
#endif

# define GAMMA_B1(C) + bb1_ ## C + bf1_ ## C
# define GAMMA_B2(C) + bb2_ ## C + bf2_ ## C

# define UPDATE_B(S, C) b(S, C, i) = b(S, C, i) GAMMA_B ## S(C)

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

  COMPLEX :: ab1, ab2, af1, af2
  COMPLEX :: bf1_1, bf2_1
  COMPLEX :: bf1_2, bf2_2
  COMPLEX :: bf1_3, bf2_3
  COMPLEX :: bb1_1, bb2_1
  COMPLEX :: bb1_2, bb2_2
  COMPLEX :: bb1_3, bb2_3

  ! statement function:
  
  COMPLEX :: i_times, c
  i_times(c) = cmplx(-aimag(c), real(c), kind = RKIND)

  TIMING_START(STRCAT(timing_bin_, NAME))

  !$omp parallel do private(jf, jb, ab1, ab2, af1, af2, &
  !$omp bf1_1, bf2_1, bf1_2, bf2_2, bf1_3, bf2_3, &
  !$omp bb1_1, bb2_1, bb1_2, bb2_2, bb1_3, bb2_3)
  do i = 1, volh
     jb = nn_bwd(i)

     ab1 = GAMMA_AB1(1)
     ab2 = GAMMA_AB2(1)

     bb1_1 = ab1 * conjg(u_o(1, 1, jb))
     bb2_1 = ab2 * conjg(u_o(1, 1, jb))
     bb1_2 = ab1 * conjg(u_o(1, 2, jb))
     bb2_2 = ab2 * conjg(u_o(1, 2, jb))
     bb1_3 = ab1 * conjg(u_o(1, 3, jb))
     bb2_3 = ab2 * conjg(u_o(1, 3, jb))
     
     jf = nn_fwd(i)
     
     af1 = GAMMA_AF1(1)
     af2 = GAMMA_AF2(1)
     
     bf1_1 = af1 * u_e(1, 1, i)
     bf2_1 = af2 * u_e(1, 1, i)
     bf1_2 = af1 * u_e(2, 1, i)
     bf2_2 = af2 * u_e(2, 1, i)
     bf1_3 = af1 * u_e(3, 1, i)
     bf2_3 = af2 * u_e(3, 1, i)
     
     ab1 = GAMMA_AB1(2)
     ab2 = GAMMA_AB2(2)
     
     bb1_1 = bb1_1 + ab1 * conjg(u_o(2, 1, jb))
     bb2_1 = bb2_1 + ab2 * conjg(u_o(2, 1, jb))
     bb1_2 = bb1_2 + ab1 * conjg(u_o(2, 2, jb))
     bb2_2 = bb2_2 + ab2 * conjg(u_o(2, 2, jb))
     bb1_3 = bb1_3 + ab1 * conjg(u_o(2, 3, jb))
     bb2_3 = bb2_3 + ab2 * conjg(u_o(2, 3, jb))
     
     af1 = GAMMA_AF1(2)
     af2 = GAMMA_AF2(2)
     
     bf1_1 = bf1_1 + af1 * u_e(1, 2, i)
     bf2_1 = bf2_1 + af2 * u_e(1, 2, i)
     bf1_2 = bf1_2 + af1 * u_e(2, 2, i)
     bf2_2 = bf2_2 + af2 * u_e(2, 2, i)
     bf1_3 = bf1_3 + af1 * u_e(3, 2, i)
     bf2_3 = bf2_3 + af2 * u_e(3, 2, i)
     
     ab1 = GAMMA_AB1(3)
     ab2 = GAMMA_AB2(3)
     
     bb1_1 = bb1_1 + ab1 * conjg(u_o(3, 1, jb))
     bb2_1 = bb2_1 + ab2 * conjg(u_o(3, 1, jb))
     bb1_2 = bb1_2 + ab1 * conjg(u_o(3, 2, jb))
     bb2_2 = bb2_2 + ab2 * conjg(u_o(3, 2, jb))
     bb1_3 = bb1_3 + ab1 * conjg(u_o(3, 3, jb))
     bb2_3 = bb2_3 + ab2 * conjg(u_o(3, 3, jb))
     
     af1 = GAMMA_AF1(3)
     af2 = GAMMA_AF2(3)
     
     bf1_1 = bf1_1 + af1 * u_e(1, 3, i)
     bf2_1 = bf2_1 + af2 * u_e(1, 3, i)
     bf1_2 = bf1_2 + af1 * u_e(2, 3, i)
     bf2_2 = bf2_2 + af2 * u_e(2, 3, i)
     bf1_3 = bf1_3 + af1 * u_e(3, 3, i)
     bf2_3 = bf2_3 + af2 * u_e(3, 3, i)
     
     
     UPDATE_B(1, 1)
     UPDATE_B(2, 1)
     UPDATE_B(3, 1)
     UPDATE_B(4, 1)

     UPDATE_B(1, 2)
     UPDATE_B(2, 2)
     UPDATE_B(3, 2)
     UPDATE_B(4, 2)
     
     UPDATE_B(1, 3)
     UPDATE_B(2, 3)
     UPDATE_B(3, 3)
     UPDATE_B(4, 3)
     
  enddo

  TIMING_STOP(STRCAT(timing_bin_, NAME))

end

!===============================================================================
