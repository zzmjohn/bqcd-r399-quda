!===============================================================================
!
! cg.F90
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
module module_cg

  type type_cg_para
     real    :: rest
     integer :: maxiter
     integer :: log
  end type type_cg_para

  type(type_cg_para), save :: cg_para

  type type_cg_stat
     integer :: niter
     integer :: niter_max
     integer :: niter_tot
     integer :: ncall
  end type type_cg_stat

  type(type_cg_stat), save :: cg_stat

  integer, save :: cg_iterations_total = 0   ! used in timing.F90
end

!-------------------------------------------------------------------------------
subroutine cg(matrix_mult, x, b, para, conf, iterations) 

  ! solves "matrix_mult * x = b" and returns number of iterations

  use      module_cg
  use      module_function_decl
  use      module_p_interface
  use      module_vol
  use      typedef_hmc
  implicit none

  external                          :: matrix_mult
  SPINCOL_OVERINDEXED,  intent(out) :: x
  SPINCOL_OVERINDEXED,  intent(in)  :: b
  type(hmc_para),       intent(in)  :: para
  type(hmc_conf),       intent(in)  :: conf
  integer,              intent(out) :: iterations

  P_SPINCOL_OVERINDEXED, save :: r, p, aap

  REAL          :: ak, bk, rtr, rtrold, paap
  integer       :: i, niter
  character(72) :: msg

  TIMING_START(timing_bin_cg)

  ALLOCATE_SC_OVERINDEXED(r)
  ALLOCATE_SC_OVERINDEXED(p)
  ALLOCATE_SC_OVERINDEXED(aap)

  call matrix_mult(r, x, para, conf)

  rtrold = ZERO
  !$omp parallel do reduction(+: rtrold)
  do i = 1, size_sc_field
     r(i) = b(i) - r(i)
     p(i) = r(i)
     rtrold = rtrold + r(i)**2
  enddo

  rtrold = global_sum(rtrold)

  do niter = 1, cg_para%maxiter
     call matrix_mult(aap, p, para, conf)

     paap = sc_dot(p, aap)
     paap = global_sum(paap)

     ak = rtrold / paap 

     rtr = ZERO
     !$omp parallel do reduction(+: rtr)
     do i = 1, size_sc_field
        x(i) = x(i) + ak * p(i)
        r(i) = r(i) - ak * aap(i)
        rtr = rtr + r(i)**2
     enddo 

     rtr = global_sum(rtr)

     if (rtr <= cg_para%rest) goto 9999

     bk = rtr / rtrold
     rtrold = rtr

     call sc_xpby(p, r, bk)  ! p = r + bk * p
  enddo

  niter = niter - 1

  if (cg_para%log /= 2) then
     write(msg, *) "cg(): no convergence; rtr = ", rtr 
     call die(msg)
  endif

9999 continue

  cg_stat%ncall = cg_stat%ncall + 1
  cg_stat%niter = niter
  cg_stat%niter_max = max(cg_stat%niter_max, niter)
  cg_stat%niter_tot = cg_stat%niter_tot + niter
  cg_iterations_total = cg_iterations_total + niter

  iterations = niter

  TIMING_STOP(timing_bin_cg)
end

!-------------------------------------------------------------------------------
subroutine init_cg_para(rest, maxiter, log)

  use      module_cg
  implicit none
  real     rest
  integer  maxiter, log

  cg_para%rest = rest
  cg_para%maxiter = maxiter
  cg_para%log = log

end

!-------------------------------------------------------------------------------
subroutine init_cg_stat()

  use      module_cg
  implicit none

  cg_stat%ncall = 0
  cg_stat%niter_max = 0
  cg_stat%niter_tot = 0

end

!-------------------------------------------------------------------------------
subroutine get_cg_stat(ncall, niter_max, niter_tot)

  use      module_cg
  implicit none
  integer  ncall, niter_max, niter_tot

  ncall     = cg_stat%ncall
  niter_max = cg_stat%niter_max
  niter_tot = cg_stat%niter_tot

end

!===============================================================================
