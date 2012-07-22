!===============================================================================
!
! module_p_interface.F90  ! interfaces of pointer manipulating routines
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
module module_p_interface

  interface

     subroutine allocate_g_field(u)
       P_GAUGE_FIELD :: u
     end subroutine allocate_g_field

     subroutine allocate_g_field_io(u)
       P_GAUGE_FIELD_IO :: u
     end subroutine allocate_g_field_io

     subroutine allocate_gen_field(p)
       P_GENERATOR_FIELD :: p
     end subroutine allocate_gen_field

     subroutine allocate_sc_field(a)
       P_SPINCOL_FIELD :: a
     end subroutine allocate_sc_field

     subroutine allocate_sc_field_io(a)
       P_SPINCOL_FIELD_IO :: a
     end subroutine allocate_sc_field_io

     subroutine allocate_sc_overindexed(a)
       P_SPINCOL_OVERINDEXED :: a
     end subroutine allocate_sc_overindexed

     subroutine allocate_sc2_field(a)
       P_SC2_FIELD :: a
     end subroutine allocate_sc2_field

     subroutine allocate_clover_field_a(a)
       use typedef_clover
       P_CLOVER_FIELD_A :: a
     end subroutine allocate_clover_field_a

     subroutine allocate_clover_field_b(b)
       use typedef_clover
       P_CLOVER_FIELD_B :: b
     end subroutine allocate_clover_field_b

     subroutine swap_p_g_field(u, v)
       P_GAUGE_FIELD :: u, v
     end subroutine swap_p_g_field

     subroutine swap_p_sc_field(a, b)
       P_SPINCOL_FIELD :: a, b
     end subroutine swap_p_sc_field

     subroutine swap_p_clover_field_a(x, y)
       use typedef_clover
       P_CLOVER_FIELD_A :: x, y
     end subroutine swap_p_clover_field_a

     subroutine swap_p_clover_field_b(x, y)
       use typedef_clover
       P_CLOVER_FIELD_B :: x, y
     end subroutine swap_p_clover_field_b

  end interface

end
!===============================================================================
