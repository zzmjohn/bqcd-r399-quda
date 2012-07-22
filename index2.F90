!===============================================================================
!
! index2.F90 - more functions for index calculations
!              these functions use modules
!
!-------------------------------------------------------------------------------
!
! Copyright (C) 2006 Hinnerk Stueben
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
subroutine local2global(home, local, global)

  use module_lattice
  implicit none

  integer, intent(in)                 :: home   ! home process of "local"
  integer, intent(in), dimension(DIM) :: local  ! local coordinates 
  integer, intent(out),dimension(DIM) :: global ! global coordinates

  integer, external                   :: i_global, i_periodic
  integer                             :: coord_home(DIM), i
  

  call unlex(home, DIM, coord_home, npe)

  do i = 1, DIM
     global(i) = i_global(local(i), N(i), coord_home(i))
     global(i) = i_periodic(global(i), L(i))
  enddo

end

!-------------------------------------------------------------------------------
integer function e_o(local)  !// returns EVEN or ODD (0 or 1)

  use module_function_decl
  implicit none

  integer, intent(in), dimension(DIM) :: local  ! local coordinates
  integer, dimension(DIM)             :: global ! global coordinates
  integer, external                   :: i_e_o

  call local2global(my_pe(), local, global)

  e_o = i_e_o(DIM, global)
end

!===============================================================================
