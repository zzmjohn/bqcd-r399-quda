!===============================================================================
!
! module_conf_info.F90
!
!-------------------------------------------------------------------------------
!
! Copyright (C) 2002 Hinnerk Stueben
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
module module_conf_info

  type type_conf_info
     REAL, dimension(2)      :: beta, kappa, csw, csw_kappa, h
     REAL                    :: plaq
     integer, dimension(DIM) :: L, bc_fermions
     integer                 :: rkind
     integer, dimension(2)   :: ensemble
  end type type_conf_info

  character(len = *), parameter :: k_format  = "Format"
  character(len = *), parameter :: k_prog    = "Program"
  character(len = *), parameter :: k_run     = "Run"
  character(len = *), parameter :: k_traj    = "Traj"
  character(len = *), parameter :: k_host    = "Host"
  character(len = *), parameter :: k_date    = "Date"
  character(len = *), parameter :: k_L       = "L"
  character(len = *), parameter :: k_bc      = "bc_fermions"
  character(len = *), parameter :: k_rkind   = "REAL_kind"
  character(len = *), parameter :: k_plaq    = "PlaqEnergy"

  character(len = *), parameter, dimension(2) ::  &
       k_ensemble  = (/ "ensemble        ", "former_ensemble " /), &
       k_beta      = (/ "beta            ", "former_beta     " /), &
       k_kappa     = (/ "kappa           ", "former_kappa    " /), &
       k_csw       = (/ "csw             ", "former_csw      " /), &
       k_csw_kappa = (/ "csw_kappa       ", "former_csw_kappa" /), &
       k_h         = (/ "h               ", "former_h        " /)

end
!===============================================================================
