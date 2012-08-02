!===============================================================================
!
! action.F90 - calculation of actions
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
REAL function sf(para, conf)  ! returns S_f

  use typedef_hmc
  use module_p_interface
  use module_vol
  implicit none

  type(hmc_para), intent(in) :: para
  type(hmc_conf), intent(in) :: conf

  P_SPINCOL_FIELD, save :: a, b

  REAL, external  :: dotprod, clover_action
  integer         :: iterations
  external        :: mtdagmt 

  if (para%kappa == ZERO) then
     sf = ZERO
  else
     ALLOCATE_SC_FIELD(a)
     ALLOCATE_SC_FIELD(b)
     
     call flip_bc(conf%u)

     call sc_copy(a, conf%phi)                             ! A = phi

     call quda_solver(mtdagmt, a, conf%phi, para, conf, iterations, 0.0) ! A = inv(M~+ M~) Phi
     !call cg(mtdagmt, a, conf%phi, para, conf, iterations) ! A = inv(M~+ M~) Phi
     call mtil(b, a, para, conf)                           ! B = M~ A

     sf = dotprod(b, b, SIZE_SC_FIELD)

     call flip_bc(conf%u)
  endif

  if (para%csw_kappa /= ZERO) sf = sf + clover_action(conf%b(1,1,ODD))
end

!-------------------------------------------------------------------------------
REAL function sg(u)  ! returns S_g

  use module_nn
  use module_vol
  implicit none

  GAUGE_FIELD :: u
  REAL :: plaq, global_sum, p
  SU3 :: uuu
  integer :: i, e, o, mu, nu, j1, j2
  REAL, external :: Re_Tr_uu

  TIMING_START(timing_bin_plaq)

  plaq = 0

  do mu = 1, DIM
     do e = EVEN, ODD
        o = EVEN + ODD - e
        do nu = mu + 1, DIM
           p = ZERO
           !$omp parallel do reduction(+: p) private(j1, j2, uuu)
           do i = 1, VOLH

              !  (j2,o) --<--   x        nu
              !    |            |        
              !    v            ^         ^
              !    |            |         |
              !  (i,e)  -->-- (j1,o)      x-->  mu


              j1 = nn(i, e, mu, FWD)
              j2 = nn(i, e, nu, FWD)
              
              uuu = 0
              call uuu_fwd(uuu, u(1, 1, j1, o, nu), &
                                u(1, 1, j2, o, mu), &
                                u(1, 1, i, e, nu))

              p = p + Re_Tr_uu(uuu, u(1, 1, i, e, mu))
           enddo
           !$omp end parallel do
           plaq = plaq + p
        enddo
     enddo
  enddo
              
  plaq = global_sum(plaq)

  sg = (6 * volume) - plaq / THREE

  TIMING_STOP(timing_bin_plaq)

end

!-------------------------------------------------------------------------------
REAL function sp(p)  ! returns action of momenta p

  use module_vol
  implicit none
  REAL, external :: dotprod
  GENERATOR_FIELD :: p
  
  integer :: mu, eo
  
  sp = ZERO
  do mu = 1, DIM
     do eo = EVEN, ODD
        sp = sp + dotprod(p(1, 1, eo, mu), p(1, 1, eo, mu), NGEN * volh)
     enddo
  enddo
  sp = sp * HALF

end

!===============================================================================
