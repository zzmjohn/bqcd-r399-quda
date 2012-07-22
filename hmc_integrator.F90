!===============================================================================
!
! hmc_integrator.F90 - integrators for the different models
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
subroutine hmc_integrator_p_ir(p, para, conf, step, calc_sf, sf1, sf2)

  use typedef_hmc
  use module_switches
  use module_vol
  implicit none

  GENERATOR_FIELD, intent(inout) :: p
  type(hmc_para),  intent(in)    :: para
  type(hmc_conf),  intent(in)    :: conf
  REAL,            intent(in)    :: step
  integer,         intent(in)    :: calc_sf
  REAL,            intent(out)   :: sf1, sf2

  select case (para%model)
     case ("A")
        call dsf1(p, conf, step, calc_sf, sf1, para)
        call dsd(p, conf, step, para)
     case ("B")
        call dsf1(p, conf, step, calc_sf, sf1, para)
        call dsf2(p, conf, step, calc_sf, sf2, para)
        call dsd(p, conf, step, para)
     case ("C")
        call dsf2(p, conf, step, calc_sf, sf2, para)
     case default
        call die("hmc_integrator_p_ir: " // para%model // ": unknown model")
  end select   

end

!-------------------------------------------------------------------------------
subroutine hmc_integrator_p_uv(p, para, conf, step, calc_sf, sf1, sf2)

  use typedef_hmc
  use module_switches
  use module_vol
  implicit none

  GENERATOR_FIELD, intent(inout) :: p
  type(hmc_para),  intent(in)    :: para
  type(hmc_conf),  intent(in)    :: conf
  REAL,            intent(in)    :: step
  integer,         intent(in)    :: calc_sf
  REAL,            intent(out)   :: sf1, sf2

  select case (para%model)
     case ("A")
        call dsg(p, conf%u, step, para%beta)
     case ("B")
        call dsg(p, conf%u, step, para%beta)
     case ("C")
        call dsg(p, conf%u, step, para%beta)
        call dsf1(p, conf, step, calc_sf, sf1, para)
        call dsd(p, conf, step, para)
     case default
        call die("hmc_integrator_p_uv: " // para%model // ": unknown model")
  end select   

end

!-------------------------------------------------------------------------------
subroutine hmc_integrator_q(p, para, conf, step)

  use typedef_hmc
  use module_switches
  use module_vol
  implicit none

  GENERATOR_FIELD, intent(in)    :: p
  type(hmc_para),  intent(in)    :: para
  type(hmc_conf),  intent(inout) :: conf
  REAL,            intent(in)    :: step

  call hmc_u(p, conf, step, para)

end

!===============================================================================
