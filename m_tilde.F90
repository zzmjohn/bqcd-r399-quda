!===============================================================================
!
! m_tilde.F90 -  matrix multiplications involving :
!
!    M~  := \tilde{M}
!    M~+ := \tilde{M}^\dagger
!
!    M~ = 1   - kappa^2 Deo Doe           (Wilson fermions)
!    M~ = Tee - kappa^2 Deo Inv(Too) Doe  (Wilson fermions + clover)
!    M~ = H   - kappa^2 Deo Inv(H) Doe    (Wilson fermions + external h)
!
! and 
!
!    \tilde{M}^\dagger =: M~+ = 1 - kappa^2 Doe+ Deo+
!
! subroutine mtil:     out = M~ in
! subroutine mtil_dag: out = M~+ in
! subroutine mtdagmt:  out = (M~+ M~) in  --> mtdagmt.F90
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
subroutine mtil(out, in, para, conf)

  use typedef_hmc
  use module_p_interface
  use module_vol
  implicit none

  type(hmc_para), intent(in)  :: para
  type(hmc_conf), intent(in)  :: conf

  SPINCOL_FIELD,  intent(out) :: out
  SPINCOL_FIELD,  intent(in)  :: in
  P_SPINCOL_FIELD, save       :: tmp   
  REAL                        :: b

  ALLOCATE_SC_FIELD(tmp)

  if (para%kappa /= ZERO) then
     call d(ODD, EVEN, tmp, in, conf%u)
     if (para%csw_kappa /= ZERO) call clover_mult_ao(conf%i(1,1,ODD), tmp, volh) 
     if (para%h /= ZERO)         call h_mult_b(-para%h, tmp, volh)
     call d(EVEN, ODD, out, tmp, conf%u)
  endif
  
  b = -para%kappa**2 / (ONE + para%h**2)

  if (para%csw_kappa /= ZERO) then
     call clover_mult_a(tmp, conf%a(1,1,EVEN), in, volh)
     call sc_xpby(out, tmp, b)                      ! out = tmp - kappa**2 * out
  else
     call sc_xpby(out, in, b)                       ! out = in  - kappa**2 * out
     if (para%h /= ZERO) call h_mult_a(out, para%h, in, volh)
  endif
  
end subroutine mtil

!-------------------------------------------------------------------------------
subroutine mtil_dag(out, in, para, conf)

  use typedef_hmc
  use module_p_interface
  use module_vol
  implicit none

  type(hmc_para), intent(in)  :: para
  type(hmc_conf), intent(in)  :: conf

  SPINCOL_FIELD,  intent(out) :: out
  SPINCOL_FIELD,  intent(in)  :: in
  P_SPINCOL_FIELD, save       :: tmp   
  REAL                        :: b

  ALLOCATE_SC_FIELD(tmp)

  if (para%kappa /= ZERO) then
                                 call d_dag(ODD, EVEN, tmp, in, conf%u)
     if (para%csw_kappa /= ZERO) call clover_mult_ao(conf%i(1,1,ODD), tmp, volh)
     if (para%h /= ZERO)         call h_mult_b(para%h, tmp, volh)
                                 call d_dag(EVEN, ODD, out, tmp, conf%u)
  endif

  b = -para%kappa**2 / (ONE + para%h**2)

  if (para%csw_kappa /= ZERO) then
     call clover_mult_a(tmp, conf%a(1,1,EVEN), in, volh)
     call sc_xpby(out, tmp, b)                      ! out = tmp - kappa**2 * out
  else
     call sc_xpby(out, in, b)                       ! out = in  - kappa**2 * out
     if (para%h /= ZERO) call h_mult_a(out, -para%h, in, volh)
  endif

end subroutine mtil_dag

!===============================================================================
