!===============================================================================
!
! hmc_u.F90  -  U := exp(i * lambda_j * P_j * step) * U
!
!-------------------------------------------------------------------------------
!
! Copyright (C) 1998-2005 Hinnerk Stueben
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
subroutine hmc_u(p, conf, step, para)

  use typedef_hmc
  use module_switches
  use module_vol
  implicit none

  GENERATOR_FIELD, intent(in)    :: p
  type(hmc_conf),  intent(inout) :: conf
  REAL,            intent(in)    :: step
  type(hmc_para),  intent(in)    :: para

  GENERATOR                      :: q
  SU3                            :: v
  integer                        :: i, mu, eo, j

  TIMING_START(timing_bin_hmc_u)

  do mu = 1, DIM
     do eo = EVEN, ODD
        !$omp parallel do private(j, q, v)
        do i = 1, VOLH
           do j = 1, NGEN
              q(j) = p(j, i, eo, mu) * step
           enddo
           call gen2u(v, q)
           call u_update(conf%u(1, 1, i, eo, mu), v)  ! u = v * u
           call u_normalize(conf%u(1, 1, i, eo, mu))
        enddo
     enddo
  enddo
  
  TIMING_STOP(timing_bin_hmc_u)

  TIMING_START(timing_bin_hmc_xbound_g)
  call xbound_g_field(conf%u)
  TIMING_STOP(timing_bin_hmc_xbound_g)

  if (switches%clover) then
     call clover_init(conf%a, conf%i, conf%b, conf%u, para%csw_kappa)
  endif
end

!===============================================================================
