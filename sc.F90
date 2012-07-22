!===============================================================================
!
! sc.F90  -  routines for the Spin-Colour field
!
!-------------------------------------------------------------------------------
!
! Copyright (C) 2003 Hinnerk Stueben
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
module module_sc_size

  integer :: sc_n_real     ! number of real numbers of an sc-field
  integer :: sc_n_complex  ! number of complex numbers of an sc-field

end

!-------------------------------------------------------------------------------
subroutine init_module_sc_size()

  use module_sc_size
  use module_vol
  implicit none
  
  sc_n_complex = NDIRAC * NCOL * volh
  sc_n_real = NDIRAC * NCOL * volh * SIZE_COMPLEX

end

!-------------------------------------------------------------------------------
subroutine sc_zero(out)

  use module_sc_size
  implicit none
  REAL, intent(out) :: out(*)
  integer           :: i
 
  TIMING_START(timing_bin_sc_zero)

  !$omp parallel do
  do i = 1, sc_n_real
     out(i) = ZERO
  enddo

  TIMING_STOP(timing_bin_sc_zero)
end

!-------------------------------------------------------------------------------
subroutine sc_copy(out, in)

  use module_sc_size
  implicit none
  REAL, intent(out) :: out(*)
  REAL, intent(in)  :: in(*)
  integer           :: i
 
  TIMING_START(timing_bin_sc_copy)

  !$omp parallel do
  do i = 1, sc_n_real
     out(i) = in(i)
  enddo

  TIMING_STOP(timing_bin_sc_copy)
end

!-------------------------------------------------------------------------------
subroutine sc_scale(inout, factor)

  use module_sc_size
  implicit none
  REAL, intent(inout) :: inout(*)
  REAL, intent(in)    :: factor
  integer             :: i
 
  TIMING_START(timing_bin_sc_scale)

  !$omp parallel do
  do i = 1, sc_n_real
     inout(i) = inout(i) * factor
  enddo

  TIMING_STOP(timing_bin_sc_scale)
end

!-------------------------------------------------------------------------------
subroutine sc_cax2(out, in1, a1, in2, a2)  ! out = a1 * in1 + a2 * in2

  use module_sc_size
  implicit none
  COMPLEX, intent(out) :: out(*)
  COMPLEX, intent(in)  :: in1(*), in2(*)
  COMPLEX, intent(in)  :: a1, a2
  integer              :: i
 
  TIMING_START(timing_bin_sc_cax2)

  !$omp parallel do
  do i = 1, sc_n_complex
     out(i) = a1 * in1(i) + a2 * in2(i)
  enddo

  TIMING_STOP(timing_bin_sc_cax2)
end

!-------------------------------------------------------------------------------
subroutine sc_axpy(inout, in, a)  ! inout = inout + a * in

  use module_sc_size
  implicit none
  REAL, intent(inout) :: inout(*)
  REAL, intent(in)    :: in(*)
  REAL, intent(in)    :: a
  integer             :: i
 
  TIMING_START(timing_bin_sc_axpy)

  !$omp parallel do
  do i = 1, sc_n_real
     inout(i) = inout(i) + a * in(i)
  enddo

  TIMING_STOP(timing_bin_sc_axpy)
end

!-------------------------------------------------------------------------------
subroutine sc_caxpy(inout, in, a)  ! inout = inout + a * in

  use module_sc_size
  implicit none
  COMPLEX, intent(inout) :: inout(*)
  COMPLEX, intent(in)    :: in(*)
  COMPLEX, intent(in)    :: a
  integer                :: i
 
  TIMING_START(timing_bin_sc_caxpy)

  !$omp parallel do
  do i = 1, sc_n_complex
     inout(i) = inout(i) + a * in(i)
  enddo

  TIMING_STOP(timing_bin_sc_caxpy)
end

!-------------------------------------------------------------------------------
subroutine sc_caxpy2(inout, in1, a1, in2, a2)  ! inout = inout + a1*in1 + a2*in2

  use module_sc_size
  implicit none
  COMPLEX, intent(inout) :: inout(*)
  COMPLEX, intent(in)    :: in1(*), in2(*)
  COMPLEX, intent(in)    :: a1, a2
  integer                :: i
 
  TIMING_START(timing_bin_sc_caxpy2)

  !$omp parallel do
  do i = 1, sc_n_complex
     inout(i) = inout(i) + a1 * in1(i) + a2 * in2(i)
  enddo

  TIMING_STOP(timing_bin_sc_caxpy2)
end

!-------------------------------------------------------------------------------
subroutine sc_xpby(inout, in, b)  ! inout = b * inout + in

  use module_sc_size
  implicit none
  REAL, intent(inout) :: inout(*)
  REAL, intent(in)    :: in(*)
  REAL, intent(in)    :: b
  integer             :: i
 
  TIMING_START(timing_bin_sc_xpby)

  !$omp parallel do
  do i = 1, sc_n_real
     inout(i) = b * inout(i) + in(i)
  enddo

  TIMING_STOP(timing_bin_sc_xpby)
end

!-------------------------------------------------------------------------------
subroutine sc_axpby(inout, in, b, a)  ! inout = b * inout + a * in

  use module_sc_size
  implicit none
  REAL, intent(inout) :: inout(*)
  REAL, intent(in)    :: in(*)
  REAL, intent(in)    :: b, a
  integer             :: i
 
  TIMING_START(timing_bin_sc_axpby)

  !$omp parallel do
  do i = 1, sc_n_real
     inout(i) = b * inout(i) + a * in(i)
  enddo

  TIMING_STOP(timing_bin_sc_axpby)
end

!-------------------------------------------------------------------------------
REAL function sc_norm2(in)  ! Sum_i abs(in_i)**2

  use module_sc_size
  implicit none
  REAL, intent(in) :: in(*)
  REAL             :: tmp 
  integer           :: i
 
  TIMING_START(timing_bin_sc_norm2)

  tmp = ZERO
  !$omp parallel do reduction(+: tmp)
  do i = 1, sc_n_real
     tmp = tmp + in(i)**2
  enddo
  
  sc_norm2 = tmp

  TIMING_STOP(timing_bin_sc_norm2)
end

!-------------------------------------------------------------------------------
REAL function sc_dot(x, y)  ! Sum_i [Re(x_i) * Re(y_i) + Im(x_i) * Im(y_i)]

  use module_sc_size
  implicit none
  REAL, intent(in) :: x(*), y(*)
  REAL             :: tmp 
  integer           :: i
 
  TIMING_START(timing_bin_sc_dot)

  tmp = ZERO
  !$omp parallel do reduction(+: tmp)
  do i = 1, sc_n_real
     tmp = tmp + x(i) * y(i)
  enddo

  sc_dot = tmp

  TIMING_STOP(timing_bin_sc_dot)
end

!-------------------------------------------------------------------------------
COMPLEX function sc_cdotc(x, y)  ! Sum_i conjg(x_i) * y_i

  use module_sc_size
  implicit none
  COMPLEX, intent(in) :: x(*), y(*)
  COMPLEX             :: tmp 
  integer             :: i
 
  TIMING_START(timing_bin_sc_cdotc)

  tmp = ZERO
  !$omp parallel do reduction(+: tmp)
  do i = 1, sc_n_complex
     tmp = tmp + conjg(x(i)) * y(i)
  enddo

  sc_cdotc = tmp

  TIMING_STOP(timing_bin_sc_cdotc)
end

!===============================================================================
