!===============================================================================
! 
! hmc_init_phi.F90  -  initialises phi, phi2 and calculates actions
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
subroutine hmc_init_phi(conf, para, sf1, sf2)

  use typedef_hmc
  use module_function_decl
  use module_p_interface
  use module_switches
  use module_vol
  implicit none

  type(hmc_para), intent(in)    :: para
  type(hmc_conf), intent(inout) :: conf
  REAL, intent(out)             :: sf1
  REAL, intent(out)             :: sf2

  P_SPINCOL_FIELD, save :: tmp   

  integer       :: iterations
  external      :: w_dagger_w

  TIMING_START(timing_bin_hmc_init_phi)

  sf1 = ZERO
  sf2 = ZERO

  if (.not. switches%dynamical) return

  ALLOCATE_SC_FIELD(tmp)

  call flip_bc(conf%u)


  call ran_gauss_volh(NDIRAC * NCOL, tmp, HALF, EVEN)    ! tmp = noise
  sf1 = dotprod(tmp, tmp, SIZE_SC_FIELD)
  call w_mult_dag(conf%phi, tmp, para, conf)             ! phi = W+ noise

  if (switches%hasenbusch) then
     call ran_gauss_volh(NDIRAC * NCOL, tmp, HALF, EVEN) ! tmp = noise
     sf2 = dotprod(tmp, tmp, SIZE_SC_FIELD)
     call mtil_dag(conf%phi2, tmp, para, conf)       ! phi2 = M~+ noise
     call quda_solver(w_dagger_w, tmp, conf%phi2, para, conf, iterations, para%rho)
     !call cg(w_dagger_w, tmp, conf%phi2, para, conf, iterations)
                                                     ! tmp = inv(W+ W) M~+ noise
     call w_mult(conf%phi2, tmp, para, conf)         ! phi2 = inv(W+) M~+ noise
  endif

  call flip_bc(conf%u)

  TIMING_STOP(timing_bin_hmc_init_phi)
end

!===============================================================================
