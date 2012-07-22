!===============================================================================
!
! hmc_init_p.F90 - initialization of momenta
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
subroutine hmc_init_p(p)

  use module_vol
  implicit none

  GENERATOR_FIELD, intent(out) :: p
  integer                      :: mu, eo

  TIMING_START(timing_bin_hmc_init_p)

  do mu = 1, DIM
     do eo = EVEN, ODD
        call ran_gauss_volh(NGEN/2, p(1,1,eo,mu), ONE, eo)
     enddo
  enddo

  TIMING_STOP(timing_bin_hmc_init_p)
end

!===============================================================================
