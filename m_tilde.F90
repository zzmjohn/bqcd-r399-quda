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
#include "quda_defs.h" 

!-------------------------------------------------------------------------------
subroutine mtil(out, in, para, conf)

  use typedef_hmc
  use module_p_interface
  use module_vol
  use typedef_quda
  implicit none

  type(hmc_para), intent(in)  :: para
  type(hmc_conf), intent(in)  :: conf

  SPINCOL_FIELD,  intent(out) :: out
  SPINCOL_FIELD,  intent(in)  :: in
  P_SPINCOL_FIELD, save       :: tmp   
  REAL                        :: b

  P_SPINCOL_FIELD       :: tmp2   
  type(quda_gauge_param) :: gauge_param
  type(quda_invert_param) :: invert_param
  integer quda

  quda = 0
  if (quda.eq.1) then
     call init_quda_gauge_param(gauge_param)
     call init_quda_invert_param(invert_param, para, 0.d0)

     call load_gauge_quda(conf%u, gauge_param)
     if (invert_param%dslash_type.eq.QUDA_CLOVER_WILSON_DSLASH) then
        call load_clover_quda(conf%a, conf%i, invert_param)
     end if

     ALLOCATE_SC_FIELD(tmp2)
  end if

  ALLOCATE_SC_FIELD(tmp)

  if (para%kappa /= ZERO) then
     call d(ODD, EVEN, tmp, in, conf%u)

     if (para%csw_kappa /= ZERO) call clover_mult_ao(conf%i(1,1,ODD), tmp, volh) 

     !write(*,*) "Aoo_inv Doe check"
     !call dslash_quda(out, in, invert_param, QUDA_ODD_PARITY)      !debug
     !call compare(tmp, out)

     if (para%h /= ZERO)         call h_mult_b(-para%h, tmp, volh)
     
     call d(EVEN, ODD, out, tmp, conf%u)

     !write(*,*) "Deo Aoo_inv Doe check"
     !invert_param%dslash_type = QUDA_WILSON_DSLASH
     !call dslash_quda(tmp, tmp, invert_param, QUDA_EVEN_PARITY)      !debug
     !invert_param%dslash_type = QUDA_CLOVER_WILSON_DSLASH
     !call compare(out, tmp)

  endif
  
  b = -para%kappa**2 / (ONE + para%h**2)

  if (para%csw_kappa /= ZERO) then
     call clover_mult_a(tmp, conf%a(1,1,EVEN), in, volh)

     !write(*,*) "Aee check"
     !call clover_quda(tmp2, in, invert_param, QUDA_EVEN_PARITY, 0)      !debug
     !call compare(tmp, tmp2)     

     call sc_xpby(out, tmp, b)                      ! out = tmp - kappa**2 * out
  else
     call sc_xpby(out, in, b)                       ! out = in  - kappa**2 * out
     if (para%h /= ZERO) call h_mult_a(out, para%h, in, volh)
  endif
  
  if (quda.eq.1) then
     write(*,*) "Aee - k^2 Deo Aoo_inv Doe check"
     call mat_quda(tmp, in, invert_param)
     call compare(out, tmp, in)
    
     if (invert_param%dslash_type.eq.QUDA_CLOVER_WILSON_DSLASH) then
        call free_clover_quda()
     endif
     
     call free_gauge_quda()
  end if

end subroutine mtil

!-------------------------------------------------------------------------------
subroutine mtil_dag(out, in, para, conf)

  use typedef_quda
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

  P_SPINCOL_FIELD       :: tmp2   
  type(quda_gauge_param) :: gauge_param
  type(quda_invert_param) :: invert_param
  integer quda

  quda = 0
  if (quda.eq.1) then
     call init_quda_gauge_param(gauge_param)
     call init_quda_invert_param(invert_param, para, 0.d0)

     call load_gauge_quda(conf%u, gauge_param)
     if (invert_param%dslash_type.eq.QUDA_CLOVER_WILSON_DSLASH) then
        call load_clover_quda(conf%a, conf%i, invert_param)
     end if

     ALLOCATE_SC_FIELD(tmp2)
  end if

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

  if (quda.eq.1) then
     invert_param%dagger = QUDA_DAG_YES
     write(*,*) "dagger(Aee - k^2 Deo Aoo_inv Doe) check"
     call mat_quda(tmp2, in, invert_param)
     call compare(out, tmp2, in)
     
     if (invert_param%dslash_type.eq.QUDA_CLOVER_WILSON_DSLASH) then
        call free_clover_quda()
     endif
     
     call free_gauge_quda()
  end if

end subroutine mtil_dag

!===============================================================================
