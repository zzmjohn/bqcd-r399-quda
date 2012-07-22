!===============================================================================
!
! uuu_f90.F90 - Fortran loops for (U * U * U)
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

!-------------------------------------------------------------------------------
subroutine uuu_bwd(r, a, b, c)  ! adds backward staple:
                                ! r = r + a^\dagger b^\dagger c
  implicit none
  SU3 :: r, a, b, c
  integer :: i, j, k, m

  do i = 1, NCOL
     do j = 1, NCOL
        do k = 1, NCOL
           do m = 1, NCOL
              r(i,j) = r(i,j) + conjg(a(k,i)) * conjg(b(m,k)) * c(m,j)
           enddo
        enddo
     enddo
  enddo

end

!-------------------------------------------------------------------------------
subroutine uuu_bwd_m(r, a, b, c)  ! subtracts backward staple:
                                  ! r = r - a^\dagger b^\dagger c
  implicit none
  SU3 :: r, a, b, c
  integer :: i, j, k, m

  do i = 1, NCOL
     do j = 1, NCOL
        do k = 1, NCOL
           do m = 1, NCOL
              r(i,j) = r(i,j) - conjg(a(k,i)) * conjg(b(m,k)) * c(m,j)
           enddo
        enddo
     enddo
  enddo

end

!-------------------------------------------------------------------------------
subroutine uuu_fwd(r, a, b, c)  ! adds forward staple:
                                ! r = r + a b^\dagger c^\dagger
  implicit none
  SU3 :: r, a, b, c
  integer :: i, j, k, m

  do i = 1, NCOL
     do j = 1, NCOL
        do k = 1, NCOL
           do m = 1, NCOL
              r(i,j) = r(i,j) + a(i,k) * conjg(b(m,k)) * conjg(c(j,m))
           enddo
        enddo
     enddo
  enddo

end

!===============================================================================
