!===============================================================================
!
! iteration_count.F90  -  counts iterations in Hasenbusch improvement,
!                         does not work with tempering
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
module module_iteration_count

  integer, save :: it_f1 = 0  ! iterations in dsf1()
  integer, save :: it_f2 = 0  ! iterations in dsf2()
end

!-----------------------------------------------------------------------------
subroutine iteration_count_f1(iter)
    
  use module_iteration_count
  implicit none
  integer :: iter
  
  it_f1 = it_f1 + iter
end

!-----------------------------------------------------------------------------
subroutine iteration_count_f2(iter)

  use module_iteration_count
  implicit none
  integer :: iter
 
  it_f2 = it_f2 + iter
end

!-----------------------------------------------------------------------------
subroutine iteration_count_write(unit)

  use module_counter
  use module_function_decl
  use module_iteration_count
  use module_switches

  implicit none

  integer, intent(in)     :: unit
  integer, save           :: written = 0 

  character(*), parameter :: key   = "%it"
  character(*), parameter :: fmt_h = "(1x, 2a, a6, 2a16)"
  character(*), parameter :: fmt_b = "(1x, a4, i6, 2i16)"

  if (switches%hasenbusch) then

     if (written == 0 .and. my_pe() == 0) then
        write(unit, fmt_h) "T", key, "traj", "iterations(F1)", "iterations(F2)"
     endif
     
     if (my_pe() == 0) write(unit, fmt_b) key, counter%traj, it_f1, it_f2
 
  endif 
    
  written = written + 1
  it_f1 = 0
  it_f2 = 0
end
!===============================================================================
