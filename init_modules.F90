!===============================================================================
!
! init_modules.F90 - initialise (some) modules
!
!-------------------------------------------------------------------------------
!
! Copyright (C) 2003-2006 Hinnerk Stueben
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
subroutine init_modules()

  call init_module_decomp()
  call init_module_lattice_io()
  call init_module_sc_size()
end

!-------------------------------------------------------------------------------
subroutine init_module_decomp()

  use module_decomp
  use module_function_decl
  use module_lattice
  use module_vol
  implicit none

  integer :: i, j, me, i_pe(DIM)
  integer :: x, y, z, t, eo
  integer :: x_std(DIM), i_std
  integer :: x_act(DIM), i_act
  integer :: me_act, me_std

  integer, external :: e_o
  integer, external :: ieo


  allocate(decomp%std%i(volh, EVEN:ODD))
  allocate(decomp%act%i(volh, EVEN:ODD))

  me = my_pe()
  call unlex(me, DIM, i_pe, NPE)

  decomp%act%L           = L
  decomp%act%NPE         = NPE
  decomp%act%N           = N
  decomp%act%NH          = NH
  decomp%act%i_pe        = i_pe
  decomp%act%bc_fermions = bc_fermions
  decomp%gamma_index     = gamma_index
  decomp%direction       = decomp_direction

  do j = 1, DIM
     i = gamma_index(j)

     decomp%std%L(i)           = L(j)
     decomp%std%NPE(i)         = NPE(j)
     decomp%std%N(i)           = N(j)
     decomp%std%i_pe(i)        = i_pe(j)
     decomp%std%bc_fermions(i) = bc_fermions(j)
  enddo

  decomp%std%NH(1) = decomp%std%N(1) / 2
  decomp%std%NH(2) = decomp%std%N(2)
  decomp%std%NH(3) = decomp%std%N(3)
  decomp%std%NH(4) = decomp%std%N(4)


  decomp%std%i = 0
  decomp%act%i = 0

  do t = 0, decomp%act%N(4) - 1
  do z = 0, decomp%act%N(3) - 1
  do y = 0, decomp%act%N(2) - 1
  do x = 0, decomp%act%N(1) - 1

     x_act(1) = x
     x_act(2) = y
     x_act(3) = z
     x_act(4) = t

     x_std(gamma_index(1)) = x_act(1)
     x_std(gamma_index(2)) = x_act(2)
     x_std(gamma_index(3)) = x_act(3)
     x_std(gamma_index(4)) = x_act(4)

     i_std = ieo(DIM, x_std, decomp%std%N) + 1
     i_act = ieo(DIM, x_act, decomp%act%N) + 1
     eo = e_o(x_act)

     decomp%std%i(i_act, eo) = i_std
     decomp%act%i(i_std, eo) = i_act
  enddo
  enddo
  enddo
  enddo

end

!-------------------------------------------------------------------------------
subroutine init_module_lattice_io()

  use module_decomp
  use module_lattice_io
  implicit none

  L   = decomp%std%L
  N   = decomp%std%N
  NH  = decomp%std%NH
  NPE = decomp%std%NPE

end

!===============================================================================
