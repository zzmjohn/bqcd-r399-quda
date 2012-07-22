!===============================================================================
!
! hmc_leap_frog.F90 - two time scale leap frog integration
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
subroutine hmc_leap_frog(p, para, conf, sf1, sf2)

  use typedef_hmc
  use module_switches
  use module_vol
  implicit none

  GENERATOR_FIELD, intent(inout) :: p
  type(hmc_para),  intent(in)    :: para
  type(hmc_conf),  intent(inout) :: conf
  REAL,            intent(out)   :: sf1, sf2

  REAL                           :: step_ir, step_uv
  integer                        :: calc_sf, itau, i_scale

  calc_sf = 0

! first half ir and uv step

  step_ir = HALF * para%tau
  step_uv = HALF * para%tau / para%m_scale

  TIMING_START(timing_bin_hmc_half_step0)
  call hmc_integrator_p_ir(p, para, conf, step_ir, calc_sf, sf1, sf2)
  call hmc_integrator_p_uv(p, para, conf, step_uv, calc_sf, sf1, sf2)
  TIMING_STOP(timing_bin_hmc_half_step0)

  step_ir = para%tau
  step_uv = para%tau / para%m_scale

  do itau = 1, para%ntau
     do i_scale = 1, para%m_scale

        call hmc_integrator_q(p, para, conf, step_uv)

        if (itau == para%ntau .and. i_scale == para%m_scale) then
           step_ir = step_ir * HALF    ! final half steps
           step_uv = step_uv * HALF 
           calc_sf = 1                 ! calculate new S_f
        endif

        TIMING_START(timing_bin_hmc_steps)
        call hmc_integrator_p_uv(p, para, conf, step_uv, calc_sf, sf1, sf2)
        TIMING_STOP(timing_bin_hmc_steps)
     enddo
   
     TIMING_START(timing_bin_hmc_steps)
     call hmc_integrator_p_ir(p, para, conf, step_ir, calc_sf, sf1, sf2)
     TIMING_STOP(timing_bin_hmc_steps)
  enddo

end

!===============================================================================
