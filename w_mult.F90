!===============================================================================
!
! w_mult.F90  -  W := M~ + rho (Hasenbusch improvement)
!
!-------------------------------------------------------------------------------
!
! Copyright (C) 2003-2005 Hinnerk Stueben
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
subroutine w_mult(out, in, para, conf)  ! out = W in

  use typedef_hmc
  use module_vol
  implicit none

  type(hmc_para), intent(in)  :: para
  type(hmc_conf), intent(in)  :: conf

  SPINCOL_FIELD,  intent(out) :: out
  SPINCOL_FIELD,  intent(in)  :: in

  call mtil(out, in, para, conf)
  call sc_axpy(out, in, para%rho)

end

!-------------------------------------------------------------------------------
subroutine w_mult_dag(out, in, para, conf)  ! out = W+ in

  use typedef_hmc
  use module_vol
  implicit none

  type(hmc_para), intent(in)  :: para
  type(hmc_conf), intent(in)  :: conf

  SPINCOL_FIELD,  intent(out) :: out
  SPINCOL_FIELD,  intent(in)  :: in

  call mtil_dag(out, in, para, conf)
  call sc_axpy(out, in, para%rho)

end

!-------------------------------------------------------------------------------
subroutine w_dagger_w(out, in, para, conf)  ! out = (W+ W) in

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

  call w_mult(tmp, in, para, conf)
  call w_mult_dag(out, tmp, para, conf)

  TIMING_STOP(timing_bin_mtdagmt)
end

!===============================================================================
