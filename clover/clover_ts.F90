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
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with BQCD.  If not, see <http://www.gnu.org/licenses/>.
!
!-------------------------------------------------------------------------------
# include "defs.h"

!-------------------------------------------------------------------------------
subroutine clover_ts(mu, nu, w, t)  ! w = t sigma_mu_nu

  use module_vol
  implicit none

  integer        :: mu, nu
  SU3_FIELD      :: w
  CLOVER_FIELD_C :: t

  if (mu == 1) then
     if     (nu == 2) then ; call clover_ts_12(w, t)
     elseif (nu == 3) then ; call clover_ts_13(w, t)
     elseif (nu == 4) then ; call clover_ts_14(w, t) ; endif
  elseif (mu == 2) then
     if     (nu == 1) then ; call clover_ts_21(w, t)
     elseif (nu == 3) then ; call clover_ts_23(w, t)
     elseif (nu == 4) then ; call clover_ts_24(w, t) ; endif
  elseif (mu == 3) then
     if     (nu == 1) then ; call clover_ts_31(w, t)
     elseif (nu == 2) then ; call clover_ts_32(w, t)
     elseif (nu == 4) then ; call clover_ts_34(w, t) ; endif
  elseif (mu == 4) then
     if     (nu == 1) then ; call clover_ts_41(w, t)
     elseif (nu == 2) then ; call clover_ts_42(w, t)
     elseif (nu == 3) then ; call clover_ts_43(w, t) ; endif
  endif
end

!-------------------------------------------------------------------------------
subroutine clover_ts_12(w, t)

# include "clover_ts_head.h90"
           w(c1, c2, i) = -t(1, c1, 1, c2, i) &
                         + t(2, c1, 2, c2, i) &
                         - t(3, c1, 3, c2, i) &
                         + t(4, c1, 4, c2, i)
# include "clover_ts_tail.h90"

!-------------------------------------------------------------------------------
subroutine clover_ts_21(w, t)

# include "clover_ts_head.h90"
           w(c1, c2, i) =  t(1, c1, 1, c2, i) &
                         - t(2, c1, 2, c2, i) &
                         + t(3, c1, 3, c2, i) &
                         - t(4, c1, 4, c2, i)
# include "clover_ts_tail.h90"

!-------------------------------------------------------------------------------
subroutine clover_ts_13(w, t)

# include "clover_ts_head.h90"
           w(c1, c2, i) =  i_times(t(1, c1, 2, c2, i)) &
                         - i_times(t(2, c1, 1, c2, i)) &
                         + i_times(t(3, c1, 4, c2, i)) &
                         - i_times(t(4, c1, 3, c2, i))
# include "clover_ts_tail.h90"

!-------------------------------------------------------------------------------
subroutine clover_ts_31(w, t)

# include "clover_ts_head.h90"
           w(c1, c2, i) = -i_times(t(1, c1, 2, c2, i)) &
                         + i_times(t(2, c1, 1, c2, i)) &
                         - i_times(t(3, c1, 4, c2, i)) &
                         + i_times(t(4, c1, 3, c2, i))
# include "clover_ts_tail.h90"

!-------------------------------------------------------------------------------
subroutine clover_ts_14(w, t)

# include "clover_ts_head.h90"
           w(c1, c2, i) = t(1, c1, 4, c2, i) &
                        + t(2, c1, 3, c2, i) &
                        + t(3, c1, 2, c2, i) &
                        + t(4, c1, 1, c2, i)
# include "clover_ts_tail.h90"

!-------------------------------------------------------------------------------
subroutine clover_ts_41(w, t)

# include "clover_ts_head.h90"
           w(c1, c2, i) = -t(1, c1, 4, c2, i) &
                         - t(2, c1, 3, c2, i) &
                         - t(3, c1, 2, c2, i) &
                         - t(4, c1, 1, c2, i)
# include "clover_ts_tail.h90"

!-------------------------------------------------------------------------------
subroutine clover_ts_23(w, t)

# include "clover_ts_head.h90"
           w(c1, c2, i) = -t(1, c1, 2, c2, i) &
                         - t(2, c1, 1, c2, i) &
                         - t(3, c1, 4, c2, i) &
                         - t(4, c1, 3, c2, i)
# include "clover_ts_tail.h90"

!-------------------------------------------------------------------------------
subroutine clover_ts_32(w, t)

# include "clover_ts_head.h90"
           w(c1, c2, i) = t(1, c1, 2, c2, i) &
                        + t(2, c1, 1, c2, i) &
                        + t(3, c1, 4, c2, i) &
                        + t(4, c1, 3, c2, i)
# include "clover_ts_tail.h90"

!-------------------------------------------------------------------------------
subroutine clover_ts_24(w, t)

# include "clover_ts_head.h90"
           w(c1, c2, i) = i_times(t(1, c1, 4, c2, i)) &
                        - i_times(t(2, c1, 3, c2, i)) &
                        + i_times(t(3, c1, 2, c2, i)) &
                        - i_times(t(4, c1, 1, c2, i))
# include "clover_ts_tail.h90"

!-------------------------------------------------------------------------------
subroutine clover_ts_42(w, t)

# include "clover_ts_head.h90"
           w(c1, c2, i) = -i_times(t(1, c1, 4, c2, i)) &
                         + i_times(t(2, c1, 3, c2, i)) &
                         - i_times(t(3, c1, 2, c2, i)) &
                         + i_times(t(4, c1, 1, c2, i))
# include "clover_ts_tail.h90"

!-------------------------------------------------------------------------------
subroutine clover_ts_34(w, t)

# include "clover_ts_head.h90"
           w(c1, c2, i) = t(1, c1, 3, c2, i) &
                        - t(2, c1, 4, c2, i) &
                        + t(3, c1, 1, c2, i) &
                        - t(4, c1, 2, c2, i)
# include "clover_ts_tail.h90"

!-------------------------------------------------------------------------------
subroutine clover_ts_43(w, t)

# include "clover_ts_head.h90"
           w(c1, c2, i) = -t(1, c1, 3, c2, i) &
                         + t(2, c1, 4, c2, i) &
                         - t(3, c1, 1, c2, i) &
                         + t(4, c1, 2, c2, i)
# include "clover_ts_tail.h90"

!===============================================================================
