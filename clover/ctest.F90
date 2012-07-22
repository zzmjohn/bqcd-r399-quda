!===============================================================================
!
! ctest.F90 - test of clover matrix multiplications: is (A * inv(A) = 1) ?
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
# include "clover.h"

!-------------------------------------------------------------------------------
program ctest

  use typedef_clover
  implicit none

  integer, parameter :: volh = 1
  integer, parameter :: nz = 1

  type(type_clover_a) :: a(2, volh), ainv(2, volh)
  type(type_clover_b) :: b(2, volh)

  COMPLEX, dimension(NDIRAC, NCOL, volh) :: z, r

  integer :: i, j, s, c

  do i = 1, volh
     do j = 1, 2
        call cinit(a(j, i))
        call clover_inv(b(j, i), ainv(j, i), a(j, i))
     enddo
  enddo


  do j = 1,12
     call zinit(z, j, volh)
     !!call clover_mult_b(b, z, volh)
     call clover_mult_ao(ainv, z, volh)
     call clover_mult_a(r, a, z, volh)
     call zwrite(r, j, volh)

!     call zinit(z, j+6, volh)
!     call clover_mult_b(b, z, volh)
!     call zwrite(z, j+6, volh)

!     call zinit(z, j+6, volh)
!     call clover_mult_a(r, a, z, volh)
!     call zwrite(r, j+6, volh)


!!     call clover_mult_b(b, r, volh)
     !!call zwrite(r, j, volh)
  enddo

end

!-------------------------------------------------------------------------------
subroutine cinit(a)

  use typedef_clover
  implicit none
  type(type_clover_a) :: a
  real, intrinsic :: ranf

  A11 = ranf()
  A22 = ranf()
  A33 = ranf()
  A44 = ranf()
  A55 = ranf()
  A66 = ranf()
  
  A12 = cmplx(ranf(), ranf())   
  A13 = cmplx(ranf(), ranf())   
  A14 = cmplx(ranf(), ranf())   
  A15 = cmplx(ranf(), ranf())   
  A16 = cmplx(ranf(), ranf())   

  A23 = cmplx(ranf(), ranf())   
  A24 = cmplx(ranf(), ranf())   
  A25 = cmplx(ranf(), ranf())   
  A26 = cmplx(ranf(), ranf())   

  A34 = cmplx(ranf(), ranf())   
  A35 = cmplx(ranf(), ranf())   
  A36 = cmplx(ranf(), ranf())   

  A45 = cmplx(ranf(), ranf())   
  A46 = cmplx(ranf(), ranf())   

  A56 = cmplx(ranf(), ranf())   

end

!-------------------------------------------------------------------------------
subroutine zinit(z, j, volh)

  implicit none
  integer :: volh, i, j
  COMPLEX, dimension(NDIRAC, NCOL, volh) :: z

  z = 0
  do i = 1, volh
     if (j == 1)  z(SC1, i) = 1
     if (j == 2)  z(SC2, i) = 1
     if (j == 3)  z(SC3, i) = 1
     if (j == 4)  z(SC4, i) = 1
     if (j == 5)  z(SC5, i) = 1
     if (j == 6)  z(SC6, i) = 1
     if (j == 7)  z(SC7, i) = 1
     if (j == 8)  z(SC8, i) = 1
     if (j == 9)  z(SC9, i) = 1
     if (j == 10) z(SC10, i) = 1
     if (j == 11) z(SC11, i) = 1
     if (j == 12) z(SC12, i) = 1
  enddo

end

!-------------------------------------------------------------------------------
subroutine zwrite(z, j, volh)

  implicit none
  integer :: volh, i, j
  COMPLEX, dimension(NDIRAC, NCOL, volh) :: z

  write(6,*) "-----------------------------------------------"
  do i = 1, volh
     write(6, "(4i4,2f16.8)") j, i, SC1, z(SC1, i)
     write(6, "(4i4,2f16.8)") j, i, SC2, z(SC2, i)
     write(6, "(4i4,2f16.8)") j, i, SC3, z(SC3, i)
     write(6, "(4i4,2f16.8)") j, i, SC4, z(SC4, i)
     write(6, "(4i4,2f16.8)") j, i, SC5, z(SC5, i)
     write(6, "(4i4,2f16.8)") j, i, SC6, z(SC6, i)
     write(6,*)
     write(6, "(4i4,2f16.8)") j, i, SC7, z(SC7, i)
     write(6, "(4i4,2f16.8)") j, i, SC8, z(SC8, i)
     write(6, "(4i4,2f16.8)") j, i, SC9, z(SC9, i)
     write(6, "(4i4,2f16.8)") j, i, SC10, z(SC10, i)
     write(6, "(4i4,2f16.8)") j, i, SC11, z(SC11, i)
     write(6, "(4i4,2f16.8)") j, i, SC12, z(SC12, i)
  enddo

end

!===============================================================================
