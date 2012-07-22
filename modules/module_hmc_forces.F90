!===============================================================================
!
! module_hmc_forces.F90
!
!-------------------------------------------------------------------------------
!
! Copyright (C) 2006 Hinnerk Stueben
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
module module_hmc_forces

  P_GENERATOR_FIELD, save :: p_old
  
  integer, parameter :: n_force = 4

  integer, parameter :: i_sg = 1
  integer, parameter :: i_sd = 2
  integer, parameter :: i_sf1 = 3
  integer, parameter :: i_sf2 = 4

  REAL, save :: f_count(n_force)
  REAL, save :: f_avg(n_force)
  REAL, save :: f_max(n_force)
end

!===============================================================================
