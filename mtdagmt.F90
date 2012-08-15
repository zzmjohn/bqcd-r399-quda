!===============================================================================
!
! mtdagmt.F90 - -> see m_tilde.F90
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
subroutine mtdagmt(out, in, para, conf)

  use typedef_hmc
  use module_p_interface
  use module_vol
  implicit none

  type(hmc_para), intent(in) :: para
  type(hmc_conf), intent(in) :: conf

  SPINCOL_FIELD         :: out, in
  P_SPINCOL_FIELD, save :: tmp   

  TIMING_START(timing_bin_mtdagmt)

  ALLOCATE_SC_FIELD(tmp)

  call mtil(tmp, in, para, conf)
  call mtil_dag(out, tmp, para, conf)

  TIMING_STOP(timing_bin_mtdagmt)
end

!===============================================================================
