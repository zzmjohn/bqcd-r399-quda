!===============================================================================
!
! DSFxyzt.F90 - routines (for standard Wilson fermions) needed in dsf.F90
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

#ifdef DIR_X
# define GAMMA_A1(C) a(1, C, J) minus i_times(a(4, C, J))
# define GAMMA_A2(C) a(2, C, J) minus i_times(a(3, C, J))
# define GAMMA_A3(C) plus i_times(a2)
# define GAMMA_A4(C) plus i_times(a1)
#endif
 
#ifdef DIR_Y
# define GAMMA_A1(C) a(1, C, J) minus a(4, C, J)
# define GAMMA_A2(C) a(2, C, J) plus  a(3, C, J)
# define GAMMA_A3(C) plus a2
# define GAMMA_A4(C) minus a1
#endif
 
#ifdef DIR_Z
# define GAMMA_A1(C) a(1, C, J) minus i_times(a(3, C, J))
# define GAMMA_A2(C) a(2, C, J) plus  i_times(a(4, C, J))
# define GAMMA_A3(C) plus i_times(a1)
# define GAMMA_A4(C) minus i_times(a2)
#endif
 
#ifdef DIR_T
#ifdef FORWARD
# define GAMMA_A1(C) ZERO
# define GAMMA_A2(C) ZERO
# define GAMMA_A3(C) TWO * a(3, C, J)
# define GAMMA_A4(C) TWO * a(4, C, J)
#else
# define GAMMA_A1(C) TWO * a(1, C, J)
# define GAMMA_A2(C) TWO * a(2, C, J)
# define GAMMA_A3(C) ZERO
# define GAMMA_A4(C) ZERO
#endif
#endif

#ifdef FORWARD
# define UU(R, A, B) uu(R, A, B)
# define plus +
# define minus -
# define I i
# define J j
#else
# define UU(R, A, B) uud(R, B, A)
# define plus -
# define minus +
# define I j
# define J i
#endif

!-------------------------------------------------------------------------------
subroutine NAME(p, b, a, s, u, nn, volh)

  implicit none

  REAL, dimension(NGEN, *), intent(inout) :: p
  REAL, intent(in) :: s
  COMPLEX, dimension (NDIRAC, NCOL, *), intent(in) :: b, a
  COMPLEX, dimension (NCOL, NCOL, *), intent(in) :: u
  INTEGER, intent(in) :: nn(*)
  integer :: volh

  integer :: i, j, ca, cb
  COMPLEX :: a1, a2, a3, a4
  SU3 :: v, w

  ! statement function:
  
  COMPLEX :: i_times, c
  i_times(c) = cmplx(-aimag(c), real(c), kind = RKIND)


  !$omp parallel do private(j,ca,cb,a1,a2,a3,a4,w,v)
  do i = 1, volh
     j = nn(i)
     do ca = 1, NCOL
        a1 = GAMMA_A1(ca)
        a2 = GAMMA_A2(ca)
        a3 = GAMMA_A3(ca)
        a4 = GAMMA_A4(ca)
        do cb = 1, NCOL
           w(ca, cb) = a1 * conjg(b(1, cb, I)) &
                     + a2 * conjg(b(2, cb, I)) &
                     + a3 * conjg(b(3, cb, I)) &
                     + a4 * conjg(b(4, cb, I))
        enddo
     enddo
     call UU(v, u(1, 1, i), w)
     call im_tr_j(p(1, i), v, minus s)
  enddo

end

!===============================================================================
