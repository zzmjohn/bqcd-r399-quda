!===============================================================================
!
! clover_uuu.F90 - multiplications of three SU(3) matrices
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
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with BQCD. If not, see <http://www.gnu.org/licenses/>.
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
subroutine clover_uuu_uuu(r, u, v, w) ! r = u * v * w

  implicit none
  complex(8), dimension (3, 3) :: r, u, v, w
  integer :: i, j, k, l

  do i = 1, 3
     do l = 1, 3
        r(i, l) = 0.0_8
        do j = 1, 3
           do k = 1, 3
              r(i, l) = r(i, l) + u(i, j) * v(j, k) * w(k, l)
           enddo
        enddo
     enddo
  enddo

end

!-------------------------------------------------------------------------------
subroutine clover_uuu_duu(r, u, v, w) ! r = u+ * v * w

  implicit none
  complex(8), dimension (3, 3) :: r, u, v, w
  integer :: i, j, k, l

  do i = 1, 3
     do l = 1, 3
        r(i, l) = 0.0_8
        do j = 1, 3
           do k = 1, 3
              r(i, l) = r(i, l) + conjg(u(j, i)) * v(j, k) * w(k, l)
           enddo
        enddo
     enddo
  enddo

end

!-------------------------------------------------------------------------------
subroutine clover_uuu_udu(r, u, v, w) ! r = u * v+ * w

  implicit none
  complex(8), dimension (3, 3) :: r, u, v, w
  integer :: i, j, k, l

  do i = 1, 3
     do l = 1, 3
        r(i, l) = 0.0_8
        do j = 1, 3
           do k = 1, 3
              r(i, l) = r(i, l) + u(i, j) * conjg(v(k, j)) * w(k, l)
           enddo
        enddo
     enddo
  enddo

end

!-------------------------------------------------------------------------------
subroutine clover_uuu_uud(r, u, v, w) ! r = u * v * w+

  implicit none
  complex(8), dimension (3, 3) :: r, u, v, w
  integer :: i, j, k, l

  do i = 1, 3
     do l = 1, 3
        r(i, l) = 0.0_8
        do j = 1, 3
           do k = 1, 3
              r(i, l) = r(i, l) + u(i, j) * v(j, k) * conjg(w(l, k))
           enddo
        enddo
     enddo
  enddo

end

!-------------------------------------------------------------------------------
subroutine clover_uuu_dud(r, u, v, w) ! r = u+ * v * w+

  implicit none
  complex(8), dimension (3, 3) :: r, u, v, w
  integer :: i, j, k, l

  do i = 1, 3
     do l = 1, 3
        r(i, l) = 0.0_8
        do j = 1, 3
           do k = 1, 3
              r(i, l) = r(i, l) + conjg(u(j, i)) * v(j, k) * conjg(w(l, k))
           enddo
        enddo
     enddo
  enddo

end

!===============================================================================
