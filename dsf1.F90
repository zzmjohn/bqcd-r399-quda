!===============================================================================
!
! dsf1.F90  --  p(j,x,mu) := p(j,x,mu) - step * D_{x,mu,j} S_{f1}
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
subroutine dsf1(p, conf, step, calc_sf, sf, para)

  use typedef_hmc
  use module_hmc_forces
  use module_mre
  use module_p_interface
  use module_switches
  use module_vol
  implicit none

  GENERATOR_FIELD, intent(inout) :: p
  type(hmc_conf),  intent(in)    :: conf
  type(hmc_para),  intent(in)    :: para
  integer,         intent(in)    :: calc_sf
  REAL,            intent(in)    :: step
  REAL,            intent(out)   :: sf

  type(type_mre),  save          :: solutions
  P_SPINCOL_FIELD, save          :: a, b
  REAL, external                 :: dotprod
  integer                        :: iterations
  external                       :: w_mult
  external                       :: w_dagger_w

  ! we need to store the value of kappa and restore after the solve
  REAL                           :: kappa

  sf = ZERO

  if (switches%quenched) return

  ALLOCATE_SC_FIELD(a)
  ALLOCATE_SC_FIELD(b)

  call flip_bc(conf%u)

  call mre_get(solutions, w_mult, a, conf%phi, para, conf)
#ifdef QUDA_SOLVER
  call quda_solver(w_dagger_w, a, conf%phi, para, conf, iterations, para%rho) ! A = inv(W+ W~) Phi
#else
  call cg(w_dagger_w, a, conf%phi, para, conf, iterations) ! A = inv(W+ W~) Phi
#endif
  call mre_put(solutions, a, calc_sf)                      ! calc_sf <=> reset
  call w_mult(b, a, para, conf)                            ! B = W~ A

  if (calc_sf /= 0) sf = dotprod(b, b, SIZE_SC_FIELD)

  call hmc_forces_old(p)
  call dsf(p, conf, step, para, a, b)
  call hmc_forces_new(p, step, i_sf1)  

  !! call flip_bc(conf%u) <- done in dsf()

  call iteration_count_f1(iterations)
end

!===============================================================================
