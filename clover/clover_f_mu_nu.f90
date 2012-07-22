!===============================================================================
!
! clover_f_mu_nu.F90 - F_mu_nu = (Q_mu_nu - h.c.) / i (missing factor 1/8)
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
!
! ^ nu
! xmp x_p (xpp) |
! |
! xm_ x xp_ x --> mu
!
! xmm x_m xpm
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
subroutine clover_f_mu_nu(f, mu, nu, x, e, u)

  use module_vol
  use module_nn
  implicit none

  complex(8), dimension (3, 3), intent(out) :: f
  integer, intent(in) :: mu, nu, x, e
  complex(8), dimension (3, 3, volh_tot, 0:1, 4), intent(in) :: u

  integer :: xmp, x_p, xm_, xp_, xmm, x_m, xpm, o

  ! statement function:

  complex(8) :: i_times, c
  i_times(c) = cmplx(-aimag(c), real(c), kind = 8)

  o = 0 + 1 - e

  xp_ = nn(x, e, mu, 0)
  xm_ = nn(x, e, mu, 1)
  x_p = nn(x, e, nu, 0)
  x_m = nn(x, e, nu, 1)

  xmp = nn(xm_, o, nu, 0)
  xmm = nn(xm_, o, nu, 1)
  xpm = nn(xp_, o, nu, 1)

  if (xmp /= nn(x_p, o, mu, 1)) call die("colver_f_mu_nu(): xmp")
  if (xmm /= nn(x_m, o, mu, 1)) call die("colver_f_mu_nu(): xmm")
  if (xpm /= nn(x_m, o, mu, 0)) call die("colver_f_mu_nu(): xpm")

  f = 0.0_8

  call clover_uuuu1(f, u(1, 1, x, e, mu), &
                       u(1, 1, xp_, o, nu), &
                       u(1, 1, x_p, o, mu), &
                       u(1, 1, x, e, nu))

  call clover_uuuu2(f, u(1, 1, x, e, nu), &
                       u(1, 1, xmp, e, mu), &
                       u(1, 1, xm_, o, nu), &
                       u(1, 1, xm_, o, mu))

  call clover_uuuu3(f, u(1, 1, xm_, o, mu), &
                       u(1, 1, xmm, e, nu), &
                       u(1, 1, xmm, e, mu), &
                       u(1, 1, x_m, o, nu))

  call clover_uuuu4(f, u(1, 1, x_m, o, nu), &
                       u(1, 1, x_m, o, mu), &
                       u(1, 1, xpm, e, nu), &
                       u(1, 1, x, e, mu))

  f(1, 1) = cmplx(2.0_8 * aimag(f(1, 1)), 0.0_8, kind = 8)
  f(2, 2) = cmplx(2.0_8 * aimag(f(2, 2)), 0.0_8, kind = 8)
  f(3, 3) = cmplx(2.0_8 * aimag(f(3, 3)), 0.0_8, kind = 8)

  f(1, 2) = i_times(conjg(f(2, 1)) - f(1, 2))
  f(1, 3) = i_times(conjg(f(3, 1)) - f(1, 3))
  f(2, 3) = i_times(conjg(f(3, 2)) - f(2, 3))

  f(2, 1) = conjg(f(1, 2))
  f(3, 1) = conjg(f(1, 3))
  f(3, 2) = conjg(f(2, 3))

end

!===============================================================================
