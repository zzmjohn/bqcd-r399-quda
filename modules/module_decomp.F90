!===============================================================================
!
! module_decomp.F90
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
module module_decomp

  type type_decomp1 
     integer, dimension(DIM)           :: L, N, NH, NPE
     integer, dimension(DIM)           :: i_pe
     integer, dimension(DIM)           :: bc_fermions
     INTEGER, dimension(:, :), pointer :: i
  end type type_decomp1

  type type_decomp2
     type(type_decomp1)      :: std  ! "standard"
     type(type_decomp1)      :: act  ! "actual" (essentially module_lattice)
     integer, dimension(DIM) :: gamma_index
     integer, dimension(DIM) :: direction
  end type type_decomp2

  type(type_decomp2), save :: decomp

end
!===============================================================================
