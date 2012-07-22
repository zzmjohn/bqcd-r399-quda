!===============================================================================
!
! xbound_single_pe.F90 - dummy routines for boundary exchange
!
!-------------------------------------------------------------------------------
!
! Copyright (C) 1998-2006 Hinnerk Stueben
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
subroutine init_xbound()
  return
end

!-------------------------------------------------------------------------------
subroutine xbound_g(u, eo, mu)

  use module_vol
  implicit none
  integer :: eo, mu
  GAUGE_FIELD :: u

  return
end

!-------------------------------------------------------------------------------
subroutine xbound_g_field(u)

  use module_vol
  implicit none
  GAUGE_FIELD :: u

  return
end

!-------------------------------------------------------------------------------
subroutine xbound_sc_field(array)

  use module_vol
  implicit none
  SPINCOL_FIELD :: array

  return
end

!-------------------------------------------------------------------------------
subroutine xbound_sc2_field(array)

  use module_vol
  implicit none
  SC2_FIELD :: array

  return
end

!-------------------------------------------------------------------------------
subroutine xbound_sc2_field_i(array)

  use module_vol
  implicit none
  SC2_FIELD :: array

  return
end

!===============================================================================
