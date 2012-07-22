!===============================================================================
!
! typedef_hmc.F90
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
module typedef_hmc

  use typedef_clover

  type hmc_para
     REAL :: beta
     REAL :: kappa
     REAL :: csw
     REAL :: csw_kappa
     REAL :: h
     REAL :: traj_length
     REAL :: tau
     REAL :: rho
     integer :: ntau
     integer :: m_scale
     character :: model
  end type hmc_para

  type hmc_para_char
     character(len = 20) :: beta
     character(len = 20) :: kappa
     character(len = 20) :: csw
     character(len = 20) :: csw_kappa
     character(len = 20) :: h
     character(len = 20) :: traj_length
     character(len = 20) :: tau
     character(len = 20) :: ntau
     character(len = 20) :: rho
     character(len = 20) :: m_scale
  end type hmc_para_char

  type hmc_out
     REAL :: exp_dh ! exp(-Delta H)
     REAL :: sg     ! without factor beta
     REAL :: sf
     integer :: accepted
     integer :: cg_ncall
     integer :: cg_niter_max
     integer :: cg_niter_tot
  end type hmc_out

  type hmc_conf
     P_GAUGE_FIELD    :: u
     P_SPINCOL_FIELD  :: phi
     P_SPINCOL_FIELD  :: phi2
     P_CLOVER_FIELD_A :: a       ! A := 1 - kappa c_sw sigma F (in 6x6 blocks)
     P_CLOVER_FIELD_A :: i       ! inverse of A
     P_CLOVER_FIELD_B :: b       ! inverse of A in (L D L+) decomposed form
     integer          :: former  ! ensemble index before tempering
  end type hmc_conf

end
!===============================================================================
