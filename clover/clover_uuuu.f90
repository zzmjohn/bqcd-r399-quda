!===============================================================================
!
! clover_uuuu.F90 - multiplications of four SU(3) matrices
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
!
! --<-- --<--
! | | | |
! v 2 ^ v 1 ^
! | | | |
! -->-- -->--
! x
! --<-- --<--
! | | | |
! v 3 ^ v 4 ^
! | | | |
! -->-- -->--
!
! uuuu1: uuuu += u1 u2 u3+ u4+ ! + = dagger
! uuuu2: uuuu += u1 u2+ u3+ u4
! uuuu3: uuuu += u1+ u2+ u3 u4
! uuuu4: uuuu += u1+ u2 u3 u4+
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
subroutine clover_uuuu1(uuuu, u1, u2, u3, u4)

  implicit none
  complex(8), dimension (3, 3) :: uuuu, u1, u2, u3, u4
  integer :: i, j, k, l, m

  do i = 1, 3
   do m = 1, 3
    do j = 1, 3
     do k = 1, 3
      do l = 1, 3
       uuuu(i,m)= uuuu(i,m)+ u1(i,j) * u2(j,k) * conjg(u3(l,k)) * conjg(u4(m,l))
      enddo
     enddo
    enddo
   enddo
  enddo
end

!-------------------------------------------------------------------------------
subroutine clover_uuuu2(uuuu, u1, u2, u3, u4)

  implicit none
  complex(8), dimension (3, 3) :: uuuu, u1, u2, u3, u4
  integer :: i, j, k, l, m

  do i = 1, 3
   do m = 1, 3
    do j = 1, 3
     do k = 1, 3
      do l = 1, 3
       uuuu(i,m)= uuuu(i,m)+ u1(i,j) * conjg(u2(k,j)) * conjg(u3(l,k)) * u4(l,m)
      enddo
     enddo
    enddo
   enddo
  enddo
end



!-------------------------------------------------------------------------------
subroutine clover_uuuu3(uuuu, u1, u2, u3, u4)

  implicit none
  complex(8), dimension (3, 3) :: uuuu, u1, u2, u3, u4
  integer :: i, j, k, l, m

  do i = 1, 3
   do m = 1, 3
    do j = 1, 3
     do k = 1, 3
      do l = 1, 3
       uuuu(i,m)= uuuu(i,m)+ conjg(u1(j,i)) * conjg(u2(k,j)) * u3(k,l) * u4(l,m)
      enddo
     enddo
    enddo
   enddo
  enddo
end

!-------------------------------------------------------------------------------
subroutine clover_uuuu4(uuuu, u1, u2, u3, u4)

  implicit none
  complex(8), dimension (3, 3) :: uuuu, u1, u2, u3, u4
  integer :: i, j, k, l, m

  do i = 1, 3
   do m = 1, 3
    do j = 1, 3
     do k = 1, 3
      do l = 1, 3
       uuuu(i,m)= uuuu(i,m)+ conjg(u1(j,i)) * u2(j,k) * u3(k,l) * conjg(u4(m,l))
      enddo
     enddo
    enddo
   enddo
  enddo
end

!===============================================================================
