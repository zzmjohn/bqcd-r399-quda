!===============================================================================
!
! dsd.F90  ---  p(j,x,mu) := p(j,x,mu) - step * D_{x,mu,j} S_det
!
!-------------------------------------------------------------------------------
!
! Copyright (C) 2003-2006 Hinnerk Stueben
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
subroutine dsd(p, conf, step, para)

  use typedef_hmc
  use module_hmc_forces
  use module_vol
  implicit none

  type(hmc_para),  intent(in)    :: para
  type(hmc_conf),  intent(in)    :: conf
  GENERATOR_FIELD, intent(inout) :: p
  REAL,            intent(in)    :: step
  REAL                           :: s

  s = -step * TWO * (para%csw_kappa / EIGHT)
  
  if (s /= ZERO) then
     call hmc_forces_old(p)
     call clover_dsd(ODD, p, conf%b, s, conf%u)
     call hmc_forces_new(p, step, i_sd)
  endif

end

!===============================================================================
