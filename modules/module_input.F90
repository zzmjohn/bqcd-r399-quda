!===============================================================================
!
! module_input.F90
!
!-------------------------------------------------------------------------------
!
! Copyright (C) 2003-2010 Hinnerk Stueben
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
!
! Simple module to read "keyword value" formatted input.  The keywords
! are defined in "module_input.h" together with a type definition and a
! default value.  The input is stored in the structure "input".  The
! components are called "input%<keyword>".  
!
! New input/keywords can be added or modified in "module_input.h".  It
! suffices to modify "module_input.h" except for input of type
! INPUT_ARRAY_ENSEMBLES for which (re)allocation and initialisation 
! routines have to be provided in this file (in subroutines "input_allocate"
! and "input_reallocate").
!
!-------------------------------------------------------------------------------
# include "defs.h"

# define INPUT_ARRAY_DIM integer, dimension(DIM)
# define INPUT_ARRAY_ENSEMBLES character(para_len), dimension(:), pointer
# define INPUT_DEFAULT_ENSEMBLES 1

!-------------------------------------------------------------------------------
module module_input

  implicit none

  !-----------------------------------------------------------------------------
  public :: type_input  ! data structure declaration
  public :: input       ! variable containing data structure
  public :: input_init  ! subroutine input_init()
  public :: input_read  ! subroutine input_read(file)
  public :: input_dump  ! subroutine input_dump(unit)
  !-----------------------------------------------------------------------------

  private

  integer, parameter :: comment_len  = 80
  integer, parameter :: keyword_len  = 32
  integer, parameter :: para_len     = 32
  integer, parameter :: word_len     =  8

  type type_input

integer, dimension(DIM) :: gamma_index

#define INPUT_INPUT(var, type, default) type :: var
#include "module_input.h"

  end type type_input

  type(type_input), save :: input

contains

  !=============================================================================
  subroutine input_init()

    call input_allocate(INPUT_DEFAULT_ENSEMBLES)
    call input_defaults()
  end subroutine input_init

  !-----------------------------------------------------------------------------
  subroutine input_read(file)

    character(*), intent(in) :: file
    integer :: iostat
    character(keyword_len) :: keyword

    call input_init()

    open(UINPUT, file = file, action = "read", status = "old")

    do
       read(UINPUT, *, iostat = iostat) keyword
       if (iostat /= 0) exit
       if (keyword(1:1) == '#') cycle
       backspace(UINPUT)
       select case (keyword)

           case("gamma_index")
               read(UINPUT, *)
               call warn("input_read(): gamma_index: ignored")

#undef INPUT_INPUT
#define INPUT_INPUT(var, type, default) \
              case (#var); read(UINPUT, *) keyword, STRCAT3(input,FORTRAN_PERCENT,var)
#include "module_input.h"

              case default
                 call die("input_read(): " // trim(keyword)  &
                                           // ": unknown keyword")

       end select

       if (keyword == "ensembles") then
          call input_reallocate(input%ensembles)
       endif

    enddo
    close(UINPUT)

  end subroutine input_read

  !-----------------------------------------------------------------------------
  subroutine input_dump(unit)
    integer :: unit

    if (.not. associated(input%beta)) call input_init()

#undef INPUT_INPUT
#define INPUT_INPUT(var, type, default) write(unit,*) #var, " ", STRCAT3(input,FORTRAN_PERCENT,var)
#include "module_input.h"

  end subroutine input_dump

  !-----------------------------------------------------------------------------
  subroutine input_defaults()

  input%gamma_index = (/1, 2, 3, 4/)

#undef INPUT_INPUT
#define INPUT_INPUT(var, type, default) STRCAT3(input,FORTRAN_PERCENT,var) = default
#include "module_input.h"

  end subroutine input_defaults

  !-----------------------------------------------------------------------------
  subroutine input_allocate(size)

    integer :: size

    allocate(input%beta(size))
    allocate(input%kappa(size))
    allocate(input%csw(size))
    allocate(input%h(size))
    allocate(input%hmc_trajectory_length(size))
    allocate(input%hmc_steps(size))
    allocate(input%hmc_rho(size))
    allocate(input%hmc_m_scale(size))
    allocate(input%start_info_file(size))

    input%beta = "0.0"
    input%kappa = "0.0"
    input%csw = "0.0"
    input%h = "0.0"
    input%hmc_trajectory_length = "1"
    input%hmc_steps = "0"
    input%hmc_rho = "0.0"
    input%hmc_m_scale = "1"
    input%start_info_file = ""

  end subroutine input_allocate


  !-----------------------------------------------------------------------------
  subroutine input_reallocate(size)

    integer :: size

    if (size > INPUT_DEFAULT_ENSEMBLES) then
       deallocate(input%beta)
       deallocate(input%kappa)
       deallocate(input%csw)
       deallocate(input%h)
       deallocate(input%hmc_trajectory_length)
       deallocate(input%hmc_steps)
       deallocate(input%hmc_rho)
       deallocate(input%hmc_m_scale)
       deallocate(input%start_info_file)

       call input_allocate(size)
    endif
  end subroutine input_reallocate

end module module_input
!===============================================================================
