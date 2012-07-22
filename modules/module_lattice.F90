!===============================================================================
!
! module_lattice.F90
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
module module_lattice

  !>> See also "module_lattice_io" !!!

  ! use a common block, because without, equivalence leads to errors
  ! with Intel Fortran90 compiler 

  integer, dimension(DIM) :: L, N, NH, NPE

  common /common_lattice/ L, N, NH, NPE

  integer :: LX, LY, LZ, LT
  integer :: NX, NY, NZ, NT, NXH

  equivalence (L(1), LX)
  equivalence (L(2), LY)
  equivalence (L(3), LZ)
  equivalence (L(4), LT)
  
  equivalence (N(1), NX)
  equivalence (N(2), NY)
  equivalence (N(3), NZ)
  equivalence (N(4), NT)

  equivalence (NH(1), NXH)

  integer, dimension(DIM), save :: bc_fermions
  integer, dimension(DIM), save :: gamma_index
  integer, dimension(DIM), save :: decomp_direction

end
!===============================================================================
