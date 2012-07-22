!===============================================================================
!
! module_function_decl.F90  -  declaration of functions
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
module module_function_decl

  ! functions in misc.F90:

  character(len=10), external :: datum
  character(len=12), external :: uhrzeit
  character(len=20), external :: rechner
  logical, external :: f_exist

  ! processes (-> comm/pes*.F90):

  integer, external :: num_pes

#ifdef CRAY
  integer, intrinsic :: my_pe  ! gives the same id as MPI_Rank in MPI_COMM_WORLD
#else
  integer, external :: my_pe
#endif

  ! global reduction (-> comm/reduction_*.F90)

  REAL, external :: dotprod
  REAL, external :: global_sum
  REAL, external :: global_min
  REAL, external :: global_max

  ! sc-field (-> sc.F90)

  REAL, external    :: sc_norm2
  REAL, external    :: sc_dot
  COMPLEX, external :: sc_cdotc

  ! ranf:

#ifdef CRAY
  real(8), intrinsic :: ranf
#else
  real(8), external :: ranf
#endif

  ! identification of D (-> d/DVersion.F90)

  integer, external :: version_of_d

  ! communication method (-> comm/comm_*.F90)

  COMM_METHOD, external :: comm_method

end
!===============================================================================
