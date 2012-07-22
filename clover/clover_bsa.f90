!===============================================================================
!
! clover_bsa.F90 - calculates "B sigma A"
!
!-------------------------------------------------------------------------------
!
! Copyright (C) 1998-2003 Hinnerk Stueben
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
subroutine clover_bsa(mu, nu, w, b, a) ! w = transposed(conjg(b) sigma_mu_nu a)

  use module_vol
  implicit none

  integer :: mu, nu
  complex(8), dimension (3, 3, volh_tot) :: w
  complex(8), dimension (4, 3, volh_tot) :: b, a

  if (mu == 1) then
     if (nu == 2) then ; call clover_bsa_12(w, b, a)
     elseif (nu == 3) then ; call clover_bsa_13(w, b, a)
     elseif (nu == 4) then ; call clover_bsa_14(w, b, a) ; endif
  elseif (mu == 2) then
     if (nu == 1) then ; call clover_bsa_21(w, b, a)
     elseif (nu == 3) then ; call clover_bsa_23(w, b, a)
     elseif (nu == 4) then ; call clover_bsa_24(w, b, a) ; endif
  elseif (mu == 3) then
     if (nu == 1) then ; call clover_bsa_31(w, b, a)
     elseif (nu == 2) then ; call clover_bsa_32(w, b, a)
     elseif (nu == 4) then ; call clover_bsa_34(w, b, a) ; endif
  elseif (mu == 4) then
     if (nu == 1) then ; call clover_bsa_41(w, b, a)
     elseif (nu == 2) then ; call clover_bsa_42(w, b, a)
     elseif (nu == 3) then ; call clover_bsa_43(w, b, a) ; endif
  endif
end

!-------------------------------------------------------------------------------
subroutine clover_bsa_12(w, b, a)
  use module_vol
  implicit none
  complex(8), dimension (3, 3, volh_tot) :: w
  complex(8), dimension (4, 3, volh_tot) :: b, a
  complex(8) :: a1, a2, a3, a4
  integer :: i, ca, cb
  ! statement function:
  complex(8) :: i_times, c
  i_times(c) = cmplx(-aimag(c), real(c), kind = 8)
  !$omp parallel do private(ca, cb, a1, a2, a3, a4)
  do i = 1, volh
     do ca = 1, 3
        a1 = -a(1, ca, i)
        a2 = a(2, ca, i)
        a3 = -a(3, ca, i)
        a4 = a(4, ca, i)
        do cb = 1, 3
           w(ca, cb, i) = a1 * conjg(b(1, cb, i)) &
                        + a2 * conjg(b(2, cb, i)) &
                        + a3 * conjg(b(3, cb, i)) &
                        + a4 * conjg(b(4, cb, i))
        enddo
     enddo
  enddo
end
!-------------------------------------------------------------------------------
subroutine clover_bsa_21(w, b, a)
  use module_vol
  implicit none
  complex(8), dimension (3, 3, volh_tot) :: w
  complex(8), dimension (4, 3, volh_tot) :: b, a
  complex(8) :: a1, a2, a3, a4
  integer :: i, ca, cb
  ! statement function:
  complex(8) :: i_times, c
  i_times(c) = cmplx(-aimag(c), real(c), kind = 8)
  !$omp parallel do private(ca, cb, a1, a2, a3, a4)
  do i = 1, volh
     do ca = 1, 3
        a1 = a(1, ca, i)
        a2 = -a(2, ca, i)
        a3 = a(3, ca, i)
        a4 = -a(4, ca, i)
        do cb = 1, 3
           w(ca, cb, i) = a1 * conjg(b(1, cb, i)) &
                        + a2 * conjg(b(2, cb, i)) &
                        + a3 * conjg(b(3, cb, i)) &
                        + a4 * conjg(b(4, cb, i))
        enddo
     enddo
  enddo
end
!-------------------------------------------------------------------------------
subroutine clover_bsa_13(w, b, a)
  use module_vol
  implicit none
  complex(8), dimension (3, 3, volh_tot) :: w
  complex(8), dimension (4, 3, volh_tot) :: b, a
  complex(8) :: a1, a2, a3, a4
  integer :: i, ca, cb
  ! statement function:
  complex(8) :: i_times, c
  i_times(c) = cmplx(-aimag(c), real(c), kind = 8)
  !$omp parallel do private(ca, cb, a1, a2, a3, a4)
  do i = 1, volh
     do ca = 1, 3
        a1 = -i_times(a(2, ca, i))
        a2 = i_times(a(1, ca, i))
        a3 = -i_times(a(4, ca, i))
        a4 = i_times(a(3, ca, i))
        do cb = 1, 3
           w(ca, cb, i) = a1 * conjg(b(1, cb, i)) &
                        + a2 * conjg(b(2, cb, i)) &
                        + a3 * conjg(b(3, cb, i)) &
                        + a4 * conjg(b(4, cb, i))
        enddo
     enddo
  enddo
end
!-------------------------------------------------------------------------------
subroutine clover_bsa_31(w, b, a)
  use module_vol
  implicit none
  complex(8), dimension (3, 3, volh_tot) :: w
  complex(8), dimension (4, 3, volh_tot) :: b, a
  complex(8) :: a1, a2, a3, a4
  integer :: i, ca, cb
  ! statement function:
  complex(8) :: i_times, c
  i_times(c) = cmplx(-aimag(c), real(c), kind = 8)
  !$omp parallel do private(ca, cb, a1, a2, a3, a4)
  do i = 1, volh
     do ca = 1, 3
        a1 = i_times(a(2, ca, i))
        a2 = -i_times(a(1, ca, i))
        a3 = i_times(a(4, ca, i))
        a4 = -i_times(a(3, ca, i))
        do cb = 1, 3
           w(ca, cb, i) = a1 * conjg(b(1, cb, i)) &
                        + a2 * conjg(b(2, cb, i)) &
                        + a3 * conjg(b(3, cb, i)) &
                        + a4 * conjg(b(4, cb, i))
        enddo
     enddo
  enddo
end
!-------------------------------------------------------------------------------
subroutine clover_bsa_14(w, b, a)
  use module_vol
  implicit none
  complex(8), dimension (3, 3, volh_tot) :: w
  complex(8), dimension (4, 3, volh_tot) :: b, a
  complex(8) :: a1, a2, a3, a4
  integer :: i, ca, cb
  ! statement function:
  complex(8) :: i_times, c
  i_times(c) = cmplx(-aimag(c), real(c), kind = 8)
  !$omp parallel do private(ca, cb, a1, a2, a3, a4)
  do i = 1, volh
     do ca = 1, 3
        a1 = a(4, ca, i)
        a2 = a(3, ca, i)
        a3 = a(2, ca, i)
        a4 = a(1, ca, i)
        do cb = 1, 3
           w(ca, cb, i) = a1 * conjg(b(1, cb, i)) &
                        + a2 * conjg(b(2, cb, i)) &
                        + a3 * conjg(b(3, cb, i)) &
                        + a4 * conjg(b(4, cb, i))
        enddo
     enddo
  enddo
end
!-------------------------------------------------------------------------------
subroutine clover_bsa_41(w, b, a)
  use module_vol
  implicit none
  complex(8), dimension (3, 3, volh_tot) :: w
  complex(8), dimension (4, 3, volh_tot) :: b, a
  complex(8) :: a1, a2, a3, a4
  integer :: i, ca, cb
  ! statement function:
  complex(8) :: i_times, c
  i_times(c) = cmplx(-aimag(c), real(c), kind = 8)
  !$omp parallel do private(ca, cb, a1, a2, a3, a4)
  do i = 1, volh
     do ca = 1, 3
        a1 = -a(4, ca, i)
        a2 = -a(3, ca, i)
        a3 = -a(2, ca, i)
        a4 = -a(1, ca, i)
        do cb = 1, 3
           w(ca, cb, i) = a1 * conjg(b(1, cb, i)) &
                        + a2 * conjg(b(2, cb, i)) &
                        + a3 * conjg(b(3, cb, i)) &
                        + a4 * conjg(b(4, cb, i))
        enddo
     enddo
  enddo
end
!-------------------------------------------------------------------------------
subroutine clover_bsa_23(w, b, a)
  use module_vol
  implicit none
  complex(8), dimension (3, 3, volh_tot) :: w
  complex(8), dimension (4, 3, volh_tot) :: b, a
  complex(8) :: a1, a2, a3, a4
  integer :: i, ca, cb
  ! statement function:
  complex(8) :: i_times, c
  i_times(c) = cmplx(-aimag(c), real(c), kind = 8)
  !$omp parallel do private(ca, cb, a1, a2, a3, a4)
  do i = 1, volh
     do ca = 1, 3
        a1 = -a(2, ca, i)
        a2 = -a(1, ca, i)
        a3 = -a(4, ca, i)
        a4 = -a(3, ca, i)
        do cb = 1, 3
           w(ca, cb, i) = a1 * conjg(b(1, cb, i)) &
                        + a2 * conjg(b(2, cb, i)) &
                        + a3 * conjg(b(3, cb, i)) &
                        + a4 * conjg(b(4, cb, i))
        enddo
     enddo
  enddo
end
!-------------------------------------------------------------------------------
subroutine clover_bsa_32(w, b, a)
  use module_vol
  implicit none
  complex(8), dimension (3, 3, volh_tot) :: w
  complex(8), dimension (4, 3, volh_tot) :: b, a
  complex(8) :: a1, a2, a3, a4
  integer :: i, ca, cb
  ! statement function:
  complex(8) :: i_times, c
  i_times(c) = cmplx(-aimag(c), real(c), kind = 8)
  !$omp parallel do private(ca, cb, a1, a2, a3, a4)
  do i = 1, volh
     do ca = 1, 3
        a1 = a(2, ca, i)
        a2 = a(1, ca, i)
        a3 = a(4, ca, i)
        a4 = a(3, ca, i)
        do cb = 1, 3
           w(ca, cb, i) = a1 * conjg(b(1, cb, i)) &
                        + a2 * conjg(b(2, cb, i)) &
                        + a3 * conjg(b(3, cb, i)) &
                        + a4 * conjg(b(4, cb, i))
        enddo
     enddo
  enddo
end
!-------------------------------------------------------------------------------
subroutine clover_bsa_24(w, b, a)
  use module_vol
  implicit none
  complex(8), dimension (3, 3, volh_tot) :: w
  complex(8), dimension (4, 3, volh_tot) :: b, a
  complex(8) :: a1, a2, a3, a4
  integer :: i, ca, cb
  ! statement function:
  complex(8) :: i_times, c
  i_times(c) = cmplx(-aimag(c), real(c), kind = 8)
  !$omp parallel do private(ca, cb, a1, a2, a3, a4)
  do i = 1, volh
     do ca = 1, 3
        a1 = -i_times(a(4, ca, i))
        a2 = i_times(a(3, ca, i))
        a3 = -i_times(a(2, ca, i))
        a4 = i_times(a(1, ca, i))
        do cb = 1, 3
           w(ca, cb, i) = a1 * conjg(b(1, cb, i)) &
                        + a2 * conjg(b(2, cb, i)) &
                        + a3 * conjg(b(3, cb, i)) &
                        + a4 * conjg(b(4, cb, i))
        enddo
     enddo
  enddo
end
!-------------------------------------------------------------------------------
subroutine clover_bsa_42(w, b, a)
  use module_vol
  implicit none
  complex(8), dimension (3, 3, volh_tot) :: w
  complex(8), dimension (4, 3, volh_tot) :: b, a
  complex(8) :: a1, a2, a3, a4
  integer :: i, ca, cb
  ! statement function:
  complex(8) :: i_times, c
  i_times(c) = cmplx(-aimag(c), real(c), kind = 8)
  !$omp parallel do private(ca, cb, a1, a2, a3, a4)
  do i = 1, volh
     do ca = 1, 3
        a1 = i_times(a(4, ca, i))
        a2 = -i_times(a(3, ca, i))
        a3 = i_times(a(2, ca, i))
        a4 = -i_times(a(1, ca, i))
        do cb = 1, 3
           w(ca, cb, i) = a1 * conjg(b(1, cb, i)) &
                        + a2 * conjg(b(2, cb, i)) &
                        + a3 * conjg(b(3, cb, i)) &
                        + a4 * conjg(b(4, cb, i))
        enddo
     enddo
  enddo
end
!-------------------------------------------------------------------------------
subroutine clover_bsa_34(w, b, a)
  use module_vol
  implicit none
  complex(8), dimension (3, 3, volh_tot) :: w
  complex(8), dimension (4, 3, volh_tot) :: b, a
  complex(8) :: a1, a2, a3, a4
  integer :: i, ca, cb
  ! statement function:
  complex(8) :: i_times, c
  i_times(c) = cmplx(-aimag(c), real(c), kind = 8)
  !$omp parallel do private(ca, cb, a1, a2, a3, a4)
  do i = 1, volh
     do ca = 1, 3
        a1 = a(3, ca, i)
        a2 = -a(4, ca, i)
        a3 = a(1, ca, i)
        a4 = -a(2, ca, i)
        do cb = 1, 3
           w(ca, cb, i) = a1 * conjg(b(1, cb, i)) &
                        + a2 * conjg(b(2, cb, i)) &
                        + a3 * conjg(b(3, cb, i)) &
                        + a4 * conjg(b(4, cb, i))
        enddo
     enddo
  enddo
end
!-------------------------------------------------------------------------------
subroutine clover_bsa_43(w, b, a)
  use module_vol
  implicit none
  complex(8), dimension (3, 3, volh_tot) :: w
  complex(8), dimension (4, 3, volh_tot) :: b, a
  complex(8) :: a1, a2, a3, a4
  integer :: i, ca, cb
  ! statement function:
  complex(8) :: i_times, c
  i_times(c) = cmplx(-aimag(c), real(c), kind = 8)
  !$omp parallel do private(ca, cb, a1, a2, a3, a4)
  do i = 1, volh
     do ca = 1, 3
        a1 = -a(3, ca, i)
        a2 = a(4, ca, i)
        a3 = -a(1, ca, i)
        a4 = a(2, ca, i)
        do cb = 1, 3
           w(ca, cb, i) = a1 * conjg(b(1, cb, i)) &
                        + a2 * conjg(b(2, cb, i)) &
                        + a3 * conjg(b(3, cb, i)) &
                        + a4 * conjg(b(4, cb, i))
        enddo
     enddo
  enddo
end
!===============================================================================
