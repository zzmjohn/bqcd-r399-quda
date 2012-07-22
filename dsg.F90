!===============================================================================
!
! dsg.F90  -  p(j,x,mu) := p(j,x,mu) - step * D_{x,mu,j} S_g
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
subroutine dsg(p, u, step, beta)

  use module_hmc_forces
  use module_vol
  implicit none

  GENERATOR_FIELD, intent(inout) :: p
  GAUGE_FIELD,     intent(in)    :: u
  REAL,            intent(in)    :: step, beta
  REAL                           :: s
  SU3                            :: uuu, w
  integer                        :: mu, eo, i

  if (beta == ZERO) return

  TIMING_START(timing_bin_dsg)

  call hmc_forces_old(p)

  s = -step * beta / THREE

  do mu = 1, DIM
     do eo = EVEN, ODD
        !$omp parallel do private(uuu, w)
        do i = 1, volh
           call staple(uuu, u, i, eo, mu)
           call uu(w, u(1, 1, i, eo, mu), uuu)
           call im_tr_j(p(1, i, eo, mu), w, s)
        enddo
     enddo
  enddo

  call hmc_forces_new(p, step, i_sg)

  TIMING_STOP(timing_bin_dsg)
end

!===============================================================================
