!===============================================================================
!
! clover_inv.F90 - calculates inverse of clover matrix
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
subroutine clover_inv(b, ainv, a)

  use typedef_clover
  implicit none
  type(type_clover_a), intent(inout) :: a
  type(type_clover_a), intent(out) :: ainv
  type(type_clover_b), intent(out) :: b

  real(8) :: d1, d2, d3, d4, d5

  ! statement function:

  complex(8) :: z
  real(8) :: sq
  sq(z) = (real(z)**2 + aimag(z)**2)

  d1 = a%i11 ! D1
  b%i11 = 1.0_8 / d1 ! 1 / D1

  b%i21 = conjg(a%i12) * b%i11 ! L21

  d2 = a%i22 - d1 * sq(b%i21) ! D2
  b%i22 = 1.0_8 / d2 ! 1 / D2

  b%i31 = conjg(a%i13) ! L31 D1
  b%i32 = conjg(a%i23) - b%i31 * conjg(b%i21) ! L32 D2

  b%i31 = b%i31 * b%i11 ! L31
  b%i32 = b%i32 * b%i22 ! L32

  d3 = a%i33 - d1 * sq(b%i31) - d2 * sq(b%i32) ! D3
  b%i33 = 1.0_8 / d3 ! 1 / D3

  b%i41 = conjg(a%i14) ! L41 D1
  b%i42 = conjg(a%i24) - b%i41 * conjg(b%i21) ! L42 D2
  b%i43 = conjg(a%i34) - b%i41 * conjg(b%i31) - b%i42 * conjg(b%i32) ! L43 D3

  b%i41 = b%i41 * b%i11 ! L41
  b%i42 = b%i42 * b%i22 ! L42
  b%i43 = b%i43 * b%i33 ! L43

  d4 = a%i44 - d1 * sq(b%i41) - d2 * sq(b%i42) - d3 * sq(b%i43) ! D4
  b%i44 = 1.0_8 / d4 ! 1 / D4

  b%i51 = conjg(a%i15)
  b%i52 = conjg(a%i25) - b%i51 * conjg(b%i21)
  b%i53 = conjg(a%i35) - b%i51 * conjg(b%i31) - b%i52 * conjg(b%i32)
  b%i54 = conjg(a%i45) - b%i51 * conjg(b%i41) - b%i52 * conjg(b%i42) - b%i53 * conjg(b%i43)

  b%i51 = b%i51 * b%i11
  b%i52 = b%i52 * b%i22
  b%i53 = b%i53 * b%i33
  b%i54 = b%i54 * b%i44

  d5 = a%i55 - d1 * sq(b%i51) - d2 * sq(b%i52) - d3 * sq(b%i53) - d4 * sq(b%i54)
  b%i55 = 1.0_8 / d5

  b%i61 = conjg(a%i16)
  b%i62 = conjg(a%i26) - b%i61 * conjg(b%i21)
  b%i63 = conjg(a%i36) - b%i61 * conjg(b%i31) - b%i62 * conjg(b%i32)
  b%i64 = conjg(a%i46) - b%i61 * conjg(b%i41) - b%i62 * conjg(b%i42) - b%i63 * conjg(b%i43)
  b%i65 = conjg(a%i56) - b%i61 * conjg(b%i51) - b%i62 * conjg(b%i52) - b%i63 * conjg(b%i53) - b%i64 * conjg(b%i54)

  b%i61 = b%i61 * b%i11
  b%i62 = b%i62 * b%i22
  b%i63 = b%i63 * b%i33
  b%i64 = b%i64 * b%i44
  b%i65 = b%i65 * b%i55

  b%i66 = a%i66 - d1 * sq(b%i61) - d2 * sq(b%i62) - d3 * sq(b%i63) &
            - d4 * sq(b%i64) - d5 * sq(b%i65)

  b%i66 = 1.0_8 / b%i66

  call clover_inv2(ainv, b)

  b%i11 = 0.5_8 * b%i11
  b%i22 = 0.5_8 * b%i22
  b%i33 = 0.5_8 * b%i33
  b%i44 = 0.5_8 * b%i44
  b%i55 = 0.5_8 * b%i55
  b%i66 = 0.5_8 * b%i66

  a%i11 = 0.5_8 * a%i11
  a%i12 = 0.5_8 * a%i12
  a%i13 = 0.5_8 * a%i13
  a%i14 = 0.5_8 * a%i14
  a%i15 = 0.5_8 * a%i15
  a%i16 = 0.5_8 * a%i16

  a%i22 = 0.5_8 * a%i22
  a%i23 = 0.5_8 * a%i23
  a%i24 = 0.5_8 * a%i24
  a%i25 = 0.5_8 * a%i25
  a%i26 = 0.5_8 * a%i26

  a%i33 = 0.5_8 * a%i33
  a%i34 = 0.5_8 * a%i34
  a%i35 = 0.5_8 * a%i35
  a%i36 = 0.5_8 * a%i36

  a%i44 = 0.5_8 * a%i44
  a%i45 = 0.5_8 * a%i45
  a%i46 = 0.5_8 * a%i46

  a%i55 = 0.5_8 * a%i55
  a%i56 = 0.5_8 * a%i56

  a%i66 = 0.5_8 * a%i66

end

!-------------------------------------------------------------------------------
subroutine clover_inv2(a, b)

  use typedef_clover
  implicit none
  type(type_clover_a), intent(out) :: a
  type(type_clover_b), intent(in) :: b

  complex(8), dimension(6) :: u, x, y

  call inv(1)
  a%i11 = real(x(1))
  a%i12 = x(2)
  a%i13 = x(3)
  a%i14 = x(4)
  a%i15 = x(5)
  a%i16 = x(6)

  call inv(2)
  a%i22 = real(x(2))
  a%i23 = x(3)
  a%i24 = x(4)
  a%i25 = x(5)
  a%i26 = x(6)

  call inv(3)
  a%i33 = real(x(3))
  a%i34 = x(4)
  a%i35 = x(5)
  a%i36 = x(6)

  call inv(4)
  a%i44 = real(x(4))
  a%i45 = x(5)
  a%i46 = x(6)

  call inv(5)
  a%i55 = real(x(5))
  a%i56 = x(6)

  call inv(6)
  a%i66 = real(x(6))


CONTAINS

  subroutine inv(i)

    integer :: i

    u = 0.0_8
    u(i) = 1.0_8

    y(1) = u(1)
    y(2) = u(2) - b%i21 * y(1)
    y(3) = u(3) - b%i31 * y(1) - b%i32 * y(2)
    y(4) = u(4) - b%i41 * y(1) - b%i42 * y(2) - b%i43 * y(3)
    y(5) = u(5) - b%i51 * y(1) - b%i52 * y(2) - b%i53 * y(3) - b%i54 * y(4)
    y(6) = u(6) - b%i61 * y(1) - b%i62 * y(2) - b%i63 * y(3) - b%i64 * y(4) - b%i65 * y(5)

    x(6) = y(6) * b%i66
    x(5) = y(5) * b%i55 - conjg(b%i65) * x(6)
    x(4) = y(4) * b%i44 - conjg(b%i54) * x(5) - conjg(b%i64) * x(6)
    x(3) = y(3) * b%i33 - conjg(b%i43) * x(4) - conjg(b%i53) * x(5) &
         - conjg(b%i63) * x(6)
    x(2) = y(2) * b%i22 - conjg(b%i32) * x(3) - conjg(b%i42) * x(4) &
         - conjg(b%i52) * x(5) - conjg(b%i62) * x(6)
    x(1) = y(1) * b%i11 - conjg(b%i21) * x(2) - conjg(b%i31) * x(3) &
         - conjg(b%i41) * x(4) - conjg(b%i51) * x(5) - conjg(b%i61) * x(6)

    x(1) = 0.5_8 * conjg(x(1))
    x(2) = 0.5_8 * conjg(x(2))
    x(3) = 0.5_8 * conjg(x(3))
    x(4) = 0.5_8 * conjg(x(4))
    x(5) = 0.5_8 * conjg(x(5))
    x(6) = 0.5_8 * conjg(x(6))

  end subroutine inv

end

!===============================================================================
