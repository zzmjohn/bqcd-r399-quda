!===============================================================================
!
! typedef_para.F90
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
module typedef_para

   use typedef_hmc

   type type_para

      integer :: run

      integer, dimension(DIM) :: L
      integer, dimension(DIM) :: NPE
      integer, dimension(DIM) :: bc_fermions
      integer, dimension(DIM) :: gamma_index

      integer :: n_temper

      type(hmc_para), dimension(MAX_TEMPER) :: hmc
      type(hmc_para_char), dimension(MAX_TEMPER) :: c_hmc

      integer :: start
      SEED    :: seed
      integer :: swap_seq

      integer :: nforce
      integer :: ntraj
      integer :: nstd
      integer :: maxtraj

      integer :: nsave

      real    :: cg_rest
      integer :: cg_maxiter
      integer :: cg_log

      character(len = 20) c_cg_rest

      FILENAME, dimension(MAX_TEMPER) :: info_file

   end type type_para

end
!===============================================================================
