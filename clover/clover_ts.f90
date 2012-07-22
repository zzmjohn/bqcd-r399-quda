!===============================================================================
!
! clover_ts.F90 - calculates T * sigma
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
subroutine clover_ts(mu, nu, w, t) ! w = t sigma_mu_nu

  use module_vol
  implicit none

  integer :: mu, nu
  complex(8), dimension (3, 3, volh_tot) :: w
  complex(8), dimension(4, 3, 4, 3, volh) :: t

  if (mu == 1) then
     if (nu == 2) then ; call clover_ts_12(w, t)
     elseif (nu == 3) then ; call clover_ts_13(w, t)
     elseif (nu == 4) then ; call clover_ts_14(w, t) ; endif
  elseif (mu == 2) then
     if (nu == 1) then ; call clover_ts_21(w, t)
     elseif (nu == 3) then ; call clover_ts_23(w, t)
     elseif (nu == 4) then ; call clover_ts_24(w, t) ; endif
  elseif (mu == 3) then
     if (nu == 1) then ; call clover_ts_31(w, t)
     elseif (nu == 2) then ; call clover_ts_32(w, t)
     elseif (nu == 4) then ; call clover_ts_34(w, t) ; endif
  elseif (mu == 4) then
     if (nu == 1) then ; call clover_ts_41(w, t)
     elseif (nu == 2) then ; call clover_ts_42(w, t)
     elseif (nu == 3) then ; call clover_ts_43(w, t) ; endif
  endif
end

!-------------------------------------------------------------------------------
subroutine clover_ts_12(w, t)
  use module_vol
  implicit none
  complex(8), dimension (3, 3, volh_tot) :: w
  complex(8), dimension(4, 3, 4, 3, volh) :: t
  integer :: i, c1, c2
  ! statement function:
  complex(8) :: i_times, c
  i_times(c) = cmplx(-aimag(c), real(c), kind = 8)
  !$omp parallel do private(c1, c2)
  do i = 1, volh
     do c2 = 1, 3
        do c1 = 1, 3
           w(c1, c2, i) = -t(1, c1, 1, c2, i) &
                         + t(2, c1, 2, c2, i) &
                         - t(3, c1, 3, c2, i) &
                         + t(4, c1, 4, c2, i)
        enddo
     enddo
  enddo
end
!-------------------------------------------------------------------------------
subroutine clover_ts_21(w, t)
  use module_vol
  implicit none
  complex(8), dimension (3, 3, volh_tot) :: w
  complex(8), dimension(4, 3, 4, 3, volh) :: t
  integer :: i, c1, c2
  ! statement function:
  complex(8) :: i_times, c
  i_times(c) = cmplx(-aimag(c), real(c), kind = 8)
  !$omp parallel do private(c1, c2)
  do i = 1, volh
     do c2 = 1, 3
        do c1 = 1, 3
           w(c1, c2, i) = t(1, c1, 1, c2, i) &
                         - t(2, c1, 2, c2, i) &
                         + t(3, c1, 3, c2, i) &
                         - t(4, c1, 4, c2, i)
        enddo
     enddo
  enddo
end
!-------------------------------------------------------------------------------
subroutine clover_ts_13(w, t)
  use module_vol
  implicit none
  complex(8), dimension (3, 3, volh_tot) :: w
  complex(8), dimension(4, 3, 4, 3, volh) :: t
  integer :: i, c1, c2
  ! statement function:
  complex(8) :: i_times, c
  i_times(c) = cmplx(-aimag(c), real(c), kind = 8)
  !$omp parallel do private(c1, c2)
  do i = 1, volh
     do c2 = 1, 3
        do c1 = 1, 3
           w(c1, c2, i) = i_times(t(1, c1, 2, c2, i)) &
                         - i_times(t(2, c1, 1, c2, i)) &
                         + i_times(t(3, c1, 4, c2, i)) &
                         - i_times(t(4, c1, 3, c2, i))
        enddo
     enddo
  enddo
end
!-------------------------------------------------------------------------------
subroutine clover_ts_31(w, t)
  use module_vol
  implicit none
  complex(8), dimension (3, 3, volh_tot) :: w
  complex(8), dimension(4, 3, 4, 3, volh) :: t
  integer :: i, c1, c2
  ! statement function:
  complex(8) :: i_times, c
  i_times(c) = cmplx(-aimag(c), real(c), kind = 8)
  !$omp parallel do private(c1, c2)
  do i = 1, volh
     do c2 = 1, 3
        do c1 = 1, 3
           w(c1, c2, i) = -i_times(t(1, c1, 2, c2, i)) &
                         + i_times(t(2, c1, 1, c2, i)) &
                         - i_times(t(3, c1, 4, c2, i)) &
                         + i_times(t(4, c1, 3, c2, i))
        enddo
     enddo
  enddo
end
!-------------------------------------------------------------------------------
subroutine clover_ts_14(w, t)
  use module_vol
  implicit none
  complex(8), dimension (3, 3, volh_tot) :: w
  complex(8), dimension(4, 3, 4, 3, volh) :: t
  integer :: i, c1, c2
  ! statement function:
  complex(8) :: i_times, c
  i_times(c) = cmplx(-aimag(c), real(c), kind = 8)
  !$omp parallel do private(c1, c2)
  do i = 1, volh
     do c2 = 1, 3
        do c1 = 1, 3
           w(c1, c2, i) = t(1, c1, 4, c2, i) &
                        + t(2, c1, 3, c2, i) &
                        + t(3, c1, 2, c2, i) &
                        + t(4, c1, 1, c2, i)
        enddo
     enddo
  enddo
end
!-------------------------------------------------------------------------------
subroutine clover_ts_41(w, t)
  use module_vol
  implicit none
  complex(8), dimension (3, 3, volh_tot) :: w
  complex(8), dimension(4, 3, 4, 3, volh) :: t
  integer :: i, c1, c2
  ! statement function:
  complex(8) :: i_times, c
  i_times(c) = cmplx(-aimag(c), real(c), kind = 8)
  !$omp parallel do private(c1, c2)
  do i = 1, volh
     do c2 = 1, 3
        do c1 = 1, 3
           w(c1, c2, i) = -t(1, c1, 4, c2, i) &
                         - t(2, c1, 3, c2, i) &
                         - t(3, c1, 2, c2, i) &
                         - t(4, c1, 1, c2, i)
        enddo
     enddo
  enddo
end
!-------------------------------------------------------------------------------
subroutine clover_ts_23(w, t)
  use module_vol
  implicit none
  complex(8), dimension (3, 3, volh_tot) :: w
  complex(8), dimension(4, 3, 4, 3, volh) :: t
  integer :: i, c1, c2
  ! statement function:
  complex(8) :: i_times, c
  i_times(c) = cmplx(-aimag(c), real(c), kind = 8)
  !$omp parallel do private(c1, c2)
  do i = 1, volh
     do c2 = 1, 3
        do c1 = 1, 3
           w(c1, c2, i) = -t(1, c1, 2, c2, i) &
                         - t(2, c1, 1, c2, i) &
                         - t(3, c1, 4, c2, i) &
                         - t(4, c1, 3, c2, i)
        enddo
     enddo
  enddo
end
!-------------------------------------------------------------------------------
subroutine clover_ts_32(w, t)
  use module_vol
  implicit none
  complex(8), dimension (3, 3, volh_tot) :: w
  complex(8), dimension(4, 3, 4, 3, volh) :: t
  integer :: i, c1, c2
  ! statement function:
  complex(8) :: i_times, c
  i_times(c) = cmplx(-aimag(c), real(c), kind = 8)
  !$omp parallel do private(c1, c2)
  do i = 1, volh
     do c2 = 1, 3
        do c1 = 1, 3
           w(c1, c2, i) = t(1, c1, 2, c2, i) &
                        + t(2, c1, 1, c2, i) &
                        + t(3, c1, 4, c2, i) &
                        + t(4, c1, 3, c2, i)
        enddo
     enddo
  enddo
end
!-------------------------------------------------------------------------------
subroutine clover_ts_24(w, t)
  use module_vol
  implicit none
  complex(8), dimension (3, 3, volh_tot) :: w
  complex(8), dimension(4, 3, 4, 3, volh) :: t
  integer :: i, c1, c2
  ! statement function:
  complex(8) :: i_times, c
  i_times(c) = cmplx(-aimag(c), real(c), kind = 8)
  !$omp parallel do private(c1, c2)
  do i = 1, volh
     do c2 = 1, 3
        do c1 = 1, 3
           w(c1, c2, i) = i_times(t(1, c1, 4, c2, i)) &
                        - i_times(t(2, c1, 3, c2, i)) &
                        + i_times(t(3, c1, 2, c2, i)) &
                        - i_times(t(4, c1, 1, c2, i))
        enddo
     enddo
  enddo
end
!-------------------------------------------------------------------------------
subroutine clover_ts_42(w, t)
  use module_vol
  implicit none
  complex(8), dimension (3, 3, volh_tot) :: w
  complex(8), dimension(4, 3, 4, 3, volh) :: t
  integer :: i, c1, c2
  ! statement function:
  complex(8) :: i_times, c
  i_times(c) = cmplx(-aimag(c), real(c), kind = 8)
  !$omp parallel do private(c1, c2)
  do i = 1, volh
     do c2 = 1, 3
        do c1 = 1, 3
           w(c1, c2, i) = -i_times(t(1, c1, 4, c2, i)) &
                         + i_times(t(2, c1, 3, c2, i)) &
                         - i_times(t(3, c1, 2, c2, i)) &
                         + i_times(t(4, c1, 1, c2, i))
        enddo
     enddo
  enddo
end
!-------------------------------------------------------------------------------
subroutine clover_ts_34(w, t)
  use module_vol
  implicit none
  complex(8), dimension (3, 3, volh_tot) :: w
  complex(8), dimension(4, 3, 4, 3, volh) :: t
  integer :: i, c1, c2
  ! statement function:
  complex(8) :: i_times, c
  i_times(c) = cmplx(-aimag(c), real(c), kind = 8)
  !$omp parallel do private(c1, c2)
  do i = 1, volh
     do c2 = 1, 3
        do c1 = 1, 3
           w(c1, c2, i) = t(1, c1, 3, c2, i) &
                        - t(2, c1, 4, c2, i) &
                        + t(3, c1, 1, c2, i) &
                        - t(4, c1, 2, c2, i)
        enddo
     enddo
  enddo
end
!-------------------------------------------------------------------------------
subroutine clover_ts_43(w, t)
  use module_vol
  implicit none
  complex(8), dimension (3, 3, volh_tot) :: w
  complex(8), dimension(4, 3, 4, 3, volh) :: t
  integer :: i, c1, c2
  ! statement function:
  complex(8) :: i_times, c
  i_times(c) = cmplx(-aimag(c), real(c), kind = 8)
  !$omp parallel do private(c1, c2)
  do i = 1, volh
     do c2 = 1, 3
        do c1 = 1, 3
           w(c1, c2, i) = -t(1, c1, 3, c2, i) &
                         + t(2, c1, 4, c2, i) &
                         - t(3, c1, 1, c2, i) &
                         + t(4, c1, 2, c2, i)
        enddo
     enddo
  enddo
end
!===============================================================================
