!===============================================================================
!
! clover_init.F90 - calculates clover matrix and its inverse
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
subroutine clover_init(a, ainv, b, u, csw_kappa)

  use typedef_clover
  use module_vol
  implicit none

  type(type_clover_a), dimension(2, volh, 0:1), intent(out) :: a, ainv
  type(type_clover_b), dimension(2, volh, 0:1), intent(out) :: b
  complex(8), dimension (3, 3, volh_tot, 0:1, 4), intent(in) :: u
  real(8), intent(in) :: csw_kappa

  integer :: i, eo
  complex(8), dimension (3, 3) :: f, g
  type(type_clover_a) :: p, q
  real(8) :: factor

  call timing_start(29)

  factor = -csw_kappa / 8.0_8

  do eo = 0, 1
     !$omp parallel do private(f, g, p, q)
     do i = 1, VOLH
        call clover_f_mu_nu(f, 2, 1, i, eo, u)

        call clover_init1(p, f)

        call clover_f_mu_nu(f, 3, 2, i, eo, u)
        call clover_f_mu_nu(g, 3, 1, i, eo, u)

        call clover_init2(p, f, g)

        call clover_f_mu_nu(f, 3, 4, i, eo, u)

        call clover_init1(q, f)

        call clover_f_mu_nu(f, 1, 4, i, eo, u)
        call clover_f_mu_nu(g, 4, 2, i, eo, u)

        call clover_init2(q, f, g)

        call clover_init3(a(1, i, eo), a(2, i, eo), p, q, factor)

        call clover_inv(b(1, i, eo), ainv(1, i, eo), a(1, i, eo))
        call clover_inv(b(2, i, eo), ainv(2, i, eo), a(2, i, eo))
     enddo
  enddo

  call timing_stop(29)

end

!-------------------------------------------------------------------------------
subroutine clover_init1(a, f)

  use typedef_clover
  implicit none
  type(type_clover_a) :: a
  complex(8), dimension (3, 3) :: f

  a%i11 = real(f(1, 1))
  a%i22 = real(f(2, 2))
  a%i33 = real(f(3, 3))

  a%i44 = -a%i11
  a%i55 = -a%i22
  a%i66 = -a%i33

  a%i12 = f(1, 2)
  a%i13 = f(1, 3)
  a%i23 = f(2, 3)

  a%i45 = -a%i12
  a%i46 = -a%i13
  a%i56 = -a%i23

end

!-------------------------------------------------------------------------------
subroutine clover_init2(a, f, g)

  use typedef_clover
  implicit none
  type(type_clover_a) :: a
  complex(8), dimension (3, 3) :: f, g

  ! statement function:

  complex(8) :: i_times, c
  i_times(c) = cmplx(-aimag(c), real(c), kind = 8)

  a%i14 = f(1, 1) + i_times(g(1, 1))
  a%i15 = f(1, 2) + i_times(g(1, 2))
  a%i16 = f(1, 3) + i_times(g(1, 3))

  a%i24 = f(2, 1) + i_times(g(2, 1))
  a%i25 = f(2, 2) + i_times(g(2, 2))
  a%i26 = f(2, 3) + i_times(g(2, 3))

  a%i34 = f(3, 1) + i_times(g(3, 1))
  a%i35 = f(3, 2) + i_times(g(3, 2))
  a%i36 = f(3, 3) + i_times(g(3, 3))

end

!-------------------------------------------------------------------------------
subroutine clover_init3(a1, a2, p, q, s)

  use typedef_clover
  implicit none
  type(type_clover_a) :: a1, a2, p, q
  real(8) :: s





! define =>
! a1%iIJ = s * (p%iIJ + q%iIJ) ; a2%iIJ = s * (p%iIJ - q%iIJ)

  a1%i11 = s * (p%i11 + q%i11) ; a2%i11 = s * (p%i11 - q%i11)
  a1%i12 = s * (p%i12 + q%i12) ; a2%i12 = s * (p%i12 - q%i12)
  a1%i13 = s * (p%i13 + q%i13) ; a2%i13 = s * (p%i13 - q%i13)
  a1%i14 = s * (p%i14 + q%i14) ; a2%i14 = s * (p%i14 - q%i14)
  a1%i15 = s * (p%i15 + q%i15) ; a2%i15 = s * (p%i15 - q%i15)
  a1%i16 = s * (p%i16 + q%i16) ; a2%i16 = s * (p%i16 - q%i16)

  a1%i22 = s * (p%i22 + q%i22) ; a2%i22 = s * (p%i22 - q%i22)
  a1%i23 = s * (p%i23 + q%i23) ; a2%i23 = s * (p%i23 - q%i23)
  a1%i24 = s * (p%i24 + q%i24) ; a2%i24 = s * (p%i24 - q%i24)
  a1%i25 = s * (p%i25 + q%i25) ; a2%i25 = s * (p%i25 - q%i25)
  a1%i26 = s * (p%i26 + q%i26) ; a2%i26 = s * (p%i26 - q%i26)

  a1%i33 = s * (p%i33 + q%i33) ; a2%i33 = s * (p%i33 - q%i33)
  a1%i34 = s * (p%i34 + q%i34) ; a2%i34 = s * (p%i34 - q%i34)
  a1%i35 = s * (p%i35 + q%i35) ; a2%i35 = s * (p%i35 - q%i35)
  a1%i36 = s * (p%i36 + q%i36) ; a2%i36 = s * (p%i36 - q%i36)

  a1%i44 = s * (p%i44 + q%i44) ; a2%i44 = s * (p%i44 - q%i44)
  a1%i45 = s * (p%i45 + q%i45) ; a2%i45 = s * (p%i45 - q%i45)
  a1%i46 = s * (p%i46 + q%i46) ; a2%i46 = s * (p%i46 - q%i46)

  a1%i55 = s * (p%i55 + q%i55) ; a2%i55 = s * (p%i55 - q%i55)
  a1%i56 = s * (p%i56 + q%i56) ; a2%i56 = s * (p%i56 - q%i56)

  a1%i66 = s * (p%i66 + q%i66) ; a2%i66 = s * (p%i66 - q%i66)

  a1%i11 = a1%i11 + 1.0_8
  a1%i22 = a1%i22 + 1.0_8
  a1%i33 = a1%i33 + 1.0_8
  a1%i44 = a1%i44 + 1.0_8
  a1%i55 = a1%i55 + 1.0_8
  a1%i66 = a1%i66 + 1.0_8

  a2%i11 = a2%i11 + 1.0_8
  a2%i22 = a2%i22 + 1.0_8
  a2%i33 = a2%i33 + 1.0_8
  a2%i44 = a2%i44 + 1.0_8
  a2%i55 = a2%i55 + 1.0_8
  a2%i66 = a2%i66 + 1.0_8

end

!===============================================================================
