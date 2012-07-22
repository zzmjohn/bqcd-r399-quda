!===============================================================================
!
! module_mre.F90
!
! Important: type(type_mre) must be defined with the "save" attribute!
!
!-------------------------------------------------------------------------------
!
! Copyright (C) 2003 Hinnerk Stueben
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
module module_mre

  integer, save :: mre_n_vec = 0

  type mre_pointer_to_sc_field
      P_SPINCOL_FIELD :: sc
  end type mre_pointer_to_sc_field

  type type_mre
      integer :: rank
      type(mre_pointer_to_sc_field), dimension(:), pointer :: vec
  end type type_mre


end
!===============================================================================
