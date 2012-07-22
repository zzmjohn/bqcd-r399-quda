!===============================================================================
!
! D2.F90 - multiplication with the Wilson hopping matrix D (or D^\dagger)
!          (optimization for Hitachi SR8000)
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
subroutine NAME(e, o, out, in, u)
 
! out := NAME in 
!
! NAME = d or d_dag
!
! out is of type "e" = EVEN or ODD
! in is of type "o" = ODD or EVEN

  use module_nn
  use module_vol
  implicit none
 
  integer :: e, o
  SPINCOL_FIELD :: out, in
  GAUGE_FIELD :: u

  TIMING_START(STRCAT(timing_bin_, NAME))

  call xbound_sc_field(in)

  call STRCAT(NAME, _t )(out, in, u(1, 1, 1, e, 4), u(1, 1, 1, o, 4), &
                                       nn(1, e, 4, FWD), nn(1, e, 4, BWD), VOLH)
  call STRCAT(NAME, _zf)(out, in, u(1, 1, 1, e, 3), u(1, 1, 1, o, 3), &
                                       nn(1, e, 3, FWD), nn(1, e, 3, BWD), VOLH)
  call STRCAT(NAME, _yf)(out, in, u(1, 1, 1, e, 2), u(1, 1, 1, o, 2), &
                                       nn(1, e, 2, FWD), nn(1, e, 2, BWD), VOLH)
  call STRCAT(NAME, _xf)(out, in, u(1, 1, 1, e, 1), u(1, 1, 1, o, 1), &
                                       nn(1, e, 1, FWD), nn(1, e, 1, BWD), VOLH)

  TIMING_STOP(STRCAT(timing_bin_, NAME))

end

!===============================================================================
