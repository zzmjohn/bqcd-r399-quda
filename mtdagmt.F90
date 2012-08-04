!===============================================================================
!
! mtdagmt.F90 - -> see m_tilde.F90
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
# include "quda_defs.h"

!-------------------------------------------------------------------------------
subroutine mtdagmt(out, in, para, conf)

  use typedef_hmc
  use module_p_interface
  use module_vol
  use typedef_quda
  implicit none

  type(hmc_para), intent(in) :: para
  type(hmc_conf), intent(in) :: conf

  REAL norm_b, norm_q
  integer s, c, quda
  type(quda_gauge_param) :: gauge_param
  type(quda_invert_param) :: invert_param


  SPINCOL_FIELD         :: out, in
  P_SPINCOL_FIELD, save :: tmp   
  P_SPINCOL_FIELD, save :: quda_tmp
  real(8) :: rtr

  TIMING_START(timing_bin_mtdagmt)

  ALLOCATE_SC_FIELD(tmp)

  quda = 0
  if (quda.eq.1) then
     ALLOCATE_SC_FIELD(quda_tmp)
     call init_quda_gauge_param(gauge_param)
     call init_quda_invert_param(invert_param, para, 0.d0)
     call load_gauge_quda(conf%u, gauge_param)
     if (invert_param%dslash_type.eq.QUDA_CLOVER_WILSON_DSLASH) then
        call load_clover_quda(conf%a, conf%i, invert_param)
     end if
  end if

  call mtil(tmp, in, para, conf)

  if (quda.eq.1) then
     write(*,*) "M check"
     invert_param%dagger = QUDA_DAG_NO
     call mat_quda(quda_tmp, in, invert_param)
     call compare(tmp, quda_tmp, in)
  end if

  call mtil_dag(out, tmp, para, conf)

  if (quda.eq.1) then
     write(*,*) "Mdag check"
     !invert_param%dagger = QUDA_DAG_YES
     call mat_dag_mat_quda(quda_tmp, in, invert_param)
     call compare(out, quda_tmp, in)
     

     if (invert_param%dslash_type.eq.QUDA_CLOVER_WILSON_DSLASH) then
        call free_clover_quda()
     endif
     
     call free_gauge_quda()
  end if

  TIMING_STOP(timing_bin_mtdagmt)
end

!===============================================================================
