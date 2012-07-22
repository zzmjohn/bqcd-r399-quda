!===============================================================================
!
! staple.F90 - calculates sum of staples for one link
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

!-------------------------------------------------------------------------------
subroutine staple(uuu, u, i, e, mu)

  use module_nn
  use module_vol
  implicit none

  SU3, intent(out) :: uuu
  GAUGE_FIELD, intent(in) :: u
  integer, intent(in) :: i, e, mu
  integer :: o, nu, j1, j2, j3, j4

  o = EVEN + ODD - e
  uuu = 0

  do nu = 1, DIM
     if (nu /= mu) then

        !  (j2,o) --<--   x        nu
        !    |            |        
        !    v            ^         ^
        !    |            |         |
        !  (i,e)  -->-- (j1,o)      x-->  mu
        !    |            |
        !    ^            v
        !    |            |
        !  (j3,o) --<-- (j4,e)


        j1 = nn(i, e, mu, FWD)
        j2 = nn(i, e, nu, FWD)
        j3 = nn(i, e, nu, BWD)
        j4 = nn(j3,o, mu, FWD)

        if (j4 /= nn(j1, o, nu, BWD)) call die('staple(): j4 inconsistent')
        if (nn(j1, o, nu, FWD) /= nn(j2, o, mu, FWD)) &
           call die('staple(): j12 inconsistent')

        call uuu_fwd(uuu, u(1, 1, j1, o, nu), &
                          u(1, 1, j2, o, mu), &
                          u(1, 1, i, e, nu))

        call uuu_bwd(uuu, u(1, 1, j4, e, nu), &
                          u(1, 1, j3, o, mu), &
                          u(1, 1, j3, o, nu))
        
     endif
  enddo

end

!===============================================================================
