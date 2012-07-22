!===============================================================================
!
! traces.F90
!
! calculates: Tr(inv(M))          psibar psi (pbp)
!             Tr(gamma5 inv(M))   psibar gamma5 psi (p5p)
!             Tr(inv(M+ M))       pion norm  (pinorm)
!
! traces of a matrix A are calculated with a stochastic estimator:
!
!     Tr(A) = eta+ A eta  (eta: Gaussian noise)
!
!-------------------------------------------------------------------------------
!
! Copyright (C) 2000-2005 Hinnerk Stueben
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
subroutine traces(para, conf, traj, i_ensemble1, i_ensemble2)

  use typedef_hmc
  use module_function_decl
  use module_p_interface
  use module_vol
  implicit none

  integer,        intent(in)  :: traj, i_ensemble1, i_ensemble2
  type(hmc_para), intent(in)  :: para
  type(hmc_conf), intent(in)  :: conf

  P_SPINCOL_FIELD, save       :: eta_e, eta_o, zeta_e, zeta_o

  character(len=*), parameter :: key_tr = "%tr"
  integer,        save        :: count = 0
  integer                     :: n_sc_field, size_of_trace
  integer                     :: cg_ncall, cg_niter_max, cg_niter_tot
  REAL                        :: pinorm
  REAL                        :: re_pbp, im_pbp, re_p5p, im_p5p
  COMPLEX                     :: pbp, p5p
  REAL                        :: res(5)
  

  ALLOCATE_SC_FIELD(eta_e)
  ALLOCATE_SC_FIELD(eta_o)
  ALLOCATE_SC_FIELD(zeta_e)
  ALLOCATE_SC_FIELD(zeta_o)

  count = count + 1
  n_sc_field = NDIRAC * NCOL * volh
  size_of_trace = NDIRAC * NCOL * volume

  call ran_gauss_volh(NDIRAC * NCOL, eta_e, HALF, EVEN)
  call ran_gauss_volh(NDIRAC * NCOL, eta_o, HALF, ODD)

  call init_cg_stat()
  call solve(para, conf, zeta_e, zeta_o, eta_e, eta_o)  ! zeta = inv(M) eta

  pbp = sc_cdotc(eta_e, zeta_e) + sc_cdotc(eta_o, zeta_o)

  pinorm = sc_dot(zeta_e, zeta_e) + sc_dot(zeta_o, zeta_o)

  call gamma5(zeta_e, volh)  ! zeta = gamma5 inv(M) eta
  call gamma5(zeta_o, volh)

  p5p = sc_cdotc(eta_e, zeta_e) + sc_cdotc(eta_o, zeta_o)

  res(1) = Re(pbp)
  res(2) = Im(pbp)
  res(3) = Re(p5p)
  res(4) = Im(p5p)
  res(5) = pinorm

  call global_sum_vec(5, res)

  re_pbp = res(1) / size_of_trace
  im_pbp = res(2) / size_of_trace
  re_p5p = res(3) / size_of_trace
  im_p5p = res(4) / size_of_trace
  pinorm = res(5) / size_of_trace

  call get_cg_stat(cg_ncall, cg_niter_max, cg_niter_tot)

  if (my_pe() == 0) then
     if (count == 1) write(UREC, 400) &
         "T", key_tr, "traj", "e", "f", &
         "Re(pbp)", "Im(pbp)", "Re(p5p)", "-Im(p5p)", "PionNorm", "CGiter"

     write(UREC, 410) key_tr, traj, i_ensemble1, i_ensemble2, &
         re_pbp, im_pbp, re_p5p, -im_p5p, pinorm, cg_niter_max
  endif


400 format (1x, 2a, a6, 2a3, 5a20,    a10)
410 format (1x, a4, i6, 2i3, 5g20.10, i10)

end

!-------------------------------------------------------------------------------
subroutine solve(para, conf, out_e, out_o, in_e, in_o)  ! solves:  M out = in

  use typedef_hmc
  use module_vol
  implicit none

  type(hmc_para), intent(in)  :: para
  type(hmc_conf), intent(in)  :: conf
  SPINCOL_FIELD,  intent(out) :: out_e, out_o
  SPINCOL_FIELD,  intent(in)  :: in_e, in_o
  
  REAL     :: a, b
  integer  :: iterations
  external :: mtdagmt

  b = para%kappa / (ONE + para%h**2)

  call h_mult_c(out_o, -para%h, in_o, volh)

  call d(EVEN, ODD, out_e, out_o, conf%u)

  call sc_xpby(out_e, in_e, b)

  call mtil_dag(out_o, out_e, para, conf)

  call cg(mtdagmt, out_e, out_o, para, conf, iterations)

  call d(ODD, EVEN, out_o, out_e, conf%u)

  a = ONE / (ONE + para%h**2)

  call sc_axpby(out_o, in_o, b, a)

  call h_mult_b(-para%h, out_o, volh)

end

!-------------------------------------------------------------------------------
subroutine gamma5(x, volh)

  implicit none
  COMPLEX, dimension (NDIRAC, *) :: x
  integer                        :: volh 

  integer :: i
  COMPLEX :: x1, x2, x3, x4

  !$omp parallel do private(x1, x2, x3, x4)
  do i = 1, NCOL * volh
     x1 = x(1, i)
     x2 = x(2, i)
     x3 = x(3, i)
     x4 = x(4, i)
     
     x(1, i) = x3
     x(2, i) = x4
     x(3, i) = x1
     x(4, i) = x2
  enddo
end

!===============================================================================
