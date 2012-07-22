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
subroutine clover_inv(b, ainv, a)

  use typedef_clover
  implicit none
  type(type_clover_a), intent(inout) :: a
  type(type_clover_a), intent(out)   :: ainv
  type(type_clover_b), intent(out)   :: b

  REAL :: d1, d2, d3, d4, d5

  ! statement function:

  COMPLEX :: z
  REAL    :: sq
  sq(z) = (Re(z)**2 + Im(z)**2)

  d1 = A11                      ! D1
  B11 = ONE / d1                ! 1 / D1
  
  B21 = A21 * B11               ! L21

  d2  = A22 - d1 * sq(B21)      ! D2
  B22 = ONE / d2                ! 1 / D2

  B31 = A31                     ! L31 D1
  B32 = A32 - B31 * B12         ! L32 D2

  B31 = B31 * B11               ! L31
  B32 = B32 * B22               ! L32

  d3  = A33 - d1 * sq(B31) - d2 * sq(B32)    ! D3
  B33 = ONE / d3                             ! 1 / D3

  B41 = A41                                  ! L41 D1
  B42 = A42 - B41 * B12                      ! L42 D2
  B43 = A43 - B41 * B13 - B42 * B23          ! L43 D3

  B41 = B41 * B11                            ! L41
  B42 = B42 * B22                            ! L42
  B43 = B43 * B33                            ! L43

  d4  = A44 - d1 * sq(B41) - d2 * sq(B42) - d3 * sq(B43)  ! D4
  B44 = ONE / d4                                          ! 1 / D4

  B51 = A51
  B52 = A52 - B51 * B12
  B53 = A53 - B51 * B13 - B52 * B23
  B54 = A54 - B51 * B14 - B52 * B24 - B53 * B34
  
  B51 = B51 * B11
  B52 = B52 * B22
  B53 = B53 * B33
  B54 = B54 * B44
  
  d5  = A55 - d1 * sq(B51) - d2 * sq(B52) - d3 * sq(B53) - d4 * sq(B54)
  B55 = ONE / d5

  B61 = A61
  B62 = A62 - B61 * B12
  B63 = A63 - B61 * B13 - B62 * B23
  B64 = A64 - B61 * B14 - B62 * B24 - B63 * B34
  B65 = A65 - B61 * B15 - B62 * B25 - B63 * B35 - B64 * B45

  B61 = B61 * B11
  B62 = B62 * B22
  B63 = B63 * B33
  B64 = B64 * B44
  B65 = B65 * B55

  B66 = A66 - d1 * sq(B61) - d2 * sq(B62) - d3 * sq(B63) &
            - d4 * sq(B64) - d5 * sq(B65)

  B66 = ONE / B66

  call clover_inv2(ainv, b)

  B11 = HALF * B11
  B22 = HALF * B22
  B33 = HALF * B33
  B44 = HALF * B44
  B55 = HALF * B55
  B66 = HALF * B66

  A11 = HALF * A11
  A12 = HALF * A12
  A13 = HALF * A13
  A14 = HALF * A14
  A15 = HALF * A15
  A16 = HALF * A16
 
  A22 = HALF * A22
  A23 = HALF * A23
  A24 = HALF * A24
  A25 = HALF * A25
  A26 = HALF * A26
 
  A33 = HALF * A33
  A34 = HALF * A34
  A35 = HALF * A35
  A36 = HALF * A36
 
  A44 = HALF * A44
  A45 = HALF * A45
  A46 = HALF * A46
 
  A55 = HALF * A55
  A56 = HALF * A56
 
  A66 = HALF * A66

end

!-------------------------------------------------------------------------------
subroutine clover_inv2(a, b)

  use typedef_clover
  implicit none
  type(type_clover_a), intent(out) :: a
  type(type_clover_b), intent(in)  :: b

  COMPLEX, dimension(6)            :: u, x, y

  call inv(1)
  A11 = Re(x(1))
  A12 = x(2)
  A13 = x(3)
  A14 = x(4)
  A15 = x(5)
  A16 = x(6)

  call inv(2)
  A22 = Re(x(2))
  A23 = x(3)
  A24 = x(4)
  A25 = x(5)
  A26 = x(6)

  call inv(3)
  A33 = Re(x(3))
  A34 = x(4)
  A35 = x(5)
  A36 = x(6)

  call inv(4)
  A44 = Re(x(4))
  A45 = x(5)
  A46 = x(6)

  call inv(5)
  A55 = Re(x(5))
  A56 = x(6)

  call inv(6)
  A66 = Re(x(6))


CONTAINS

  subroutine inv(i)

    integer :: i

    u = ZERO
    u(i) = ONE

    y(1) = u(1)
    y(2) = u(2) - B21 * y(1)
    y(3) = u(3) - B31 * y(1) - B32 * y(2)
    y(4) = u(4) - B41 * y(1) - B42 * y(2) - B43 * y(3)
    y(5) = u(5) - B51 * y(1) - B52 * y(2) - B53 * y(3) - B54 * y(4)
    y(6) = u(6) - B61 * y(1) - B62 * y(2) - B63 * y(3) - B64 * y(4) - B65 * y(5)

    x(6) = y(6) * B66
    x(5) = y(5) * B55 - B56 * x(6)
    x(4) = y(4) * B44 - B45 * x(5) - B46 * x(6)
    x(3) = y(3) * B33 - B34 * x(4) - B35 * x(5) &
         - B36 * x(6)
    x(2) = y(2) * B22 - B23 * x(3) - B24 * x(4) &
         - B25 * x(5) - B26 * x(6)
    x(1) = y(1) * B11 - B12 * x(2) - B13 * x(3) &
         - B14 * x(4) - B15 * x(5) - B16 * x(6)

    x(1) = HALF * conjg(x(1))
    x(2) = HALF * conjg(x(2))
    x(3) = HALF * conjg(x(3))
    x(4) = HALF * conjg(x(4))
    x(5) = HALF * conjg(x(5))
    x(6) = HALF * conjg(x(6))

  end subroutine inv

end

!===============================================================================
