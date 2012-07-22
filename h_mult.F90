!===============================================================================
!
! h_mult.F90
!
!-------------------------------------------------------------------------------
!
! Copyright (C) 2000-2002 Hinnerk Stueben
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
subroutine h_mult_a(out, h, in, volh)  ! out := out + i h gamma_5 in

  implicit none
  COMPLEX, dimension (NDIRAC, *) :: out, in
  REAL                           :: h
  integer                        :: volh

  integer :: i

  ! statement function:
  
  COMPLEX :: i_times, c
  i_times(c) = cmplx(-aimag(c), real(c), kind = RKIND)

  TIMING_START(timing_bin_h_mult_a)

  !$omp parallel do
  do i = 1, NCOL * volh
     out(1, i) = out(1, i) + h * i_times(in(3, i))
     out(2, i) = out(2, i) + h * i_times(in(4, i))
     out(3, i) = out(3, i) + h * i_times(in(1, i))
     out(4, i) = out(4, i) + h * i_times(in(2, i))
  enddo

  TIMING_STOP(timing_bin_h_mult_a)
end

!-------------------------------------------------------------------------------
subroutine h_mult_b(h, x, volh)  ! x := (1 + i h gamma_5) x

  implicit none
  REAL                           :: h
  COMPLEX, dimension (NDIRAC, *) :: x
  integer                        :: volh 

  integer :: i
  COMPLEX :: x1, x2, x3, x4

  ! statement function:
  
  COMPLEX :: i_times, c
  i_times(c) = cmplx(-aimag(c), real(c), kind = RKIND)


  TIMING_START(timing_bin_h_mult_b)

  !$omp parallel do private(x1, x2, x3, x4)
  do i = 1, NCOL * volh
     x1 = x(1, i)
     x2 = x(2, i)
     x3 = x(3, i)
     x4 = x(4, i)
     
     x(1, i) = x(1, i) + h * i_times(x3)
     x(2, i) = x(2, i) + h * i_times(x4)
     x(3, i) = x(3, i) + h * i_times(x1)
     x(4, i) = x(4, i) + h * i_times(x2)
  enddo

  TIMING_STOP(timing_bin_h_mult_b)
end

!-------------------------------------------------------------------------------
subroutine h_mult_c(out, h, in, volh)  ! out = (1 + i h gamma_5) in

  implicit none
  COMPLEX, dimension (NDIRAC, *) :: out, in
  REAL                           :: h
  integer                        :: volh

  integer :: i

  ! statement function:
  
  COMPLEX :: i_times, c
  i_times(c) = cmplx(-aimag(c), real(c), kind = RKIND)

  TIMING_START(timing_bin_h_mult_c)

  !$omp parallel do
  do i = 1, NCOL * volh
     out(1, i) = in(1, i) + h * i_times(in(3, i))
     out(2, i) = in(2, i) + h * i_times(in(4, i))
     out(3, i) = in(3, i) + h * i_times(in(1, i))
     out(4, i) = in(4, i) + h * i_times(in(2, i))
  enddo

  TIMING_STOP(timing_bin_h_mult_c)
end

!===============================================================================
