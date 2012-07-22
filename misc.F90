!===============================================================================
!
! misc.F90 - miscellaneous (service routines)
!
!-------------------------------------------------------------------------------
!
! Copyright (C) 1998-2006 Hinnerk Stueben
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
subroutine die(msg)  ! write "msg" to stderr and abort

  implicit none
  character(len = *) msg

  write(STDERR,*) msg
  call abbruch()

end

!-------------------------------------------------------------------------------
subroutine warn(msg)  ! write "msg" to stderr and unit UREC
  implicit none
  character(len = *) msg

  write(STDERR,*) "WARNING: ", msg, " !!!"
  write(UREC,*)   "WARNING: ", msg, " !!!"
end

!-------------------------------------------------------------------------------
subroutine assertion_failed(file, line, condition)  !// used by ASSERT makro
  implicit none
  character(len = *) file, condition
  integer line

  character(len = 8) aline

  write(aline, *) line

  call die(file // " (" // trim(aline) // "): assertion failed: " // condition)

end

!-------------------------------------------------------------------------------
subroutine begin(unit, str)  ! write "begin marker"

  use module_function_decl
  implicit none
  integer :: unit
  character(len = *) :: str

  if (my_pe() == 0) write(unit, "(1x,2a)") ">Begin", str

end

!-------------------------------------------------------------------------------
subroutine end(unit, str)  ! write "end marker"

  use module_function_decl
  implicit none
  integer :: unit
  character(len = *) :: str

  if (my_pe() == 0) write(unit, "(1x,2a)") ">End", str

end

!-------------------------------------------------------------------------------
function datum()  ! returns date as: YYYY-MM-DD

  implicit none
  character(len = 10) :: datum

  call date_and_time(date = datum)

  datum(9:10) = datum(7:8)
  datum(6:7) = datum(5:6)

  datum(5:5) = "-"
  datum(8:8) = "-"

end

!-------------------------------------------------------------------------------
function uhrzeit()  ! returns time as:  hh:mm:ss.sss

  implicit none
  character(len = 12) :: uhrzeit

  call date_and_time(time = uhrzeit)

  uhrzeit(7:12) = uhrzeit(5:10)
  uhrzeit(4:5) = uhrzeit(3:4)

  uhrzeit(3:3) = ":"
  uhrzeit(6:6) = ":"

end

!-------------------------------------------------------------------------------
function f_exist(file) ! check in file exists

  implicit none
  logical :: f_exist
  character(len = *) :: file
  
  inquire(file = file, exist = f_exist)

end

!-------------------------------------------------------------------------------
subroutine open_diag()  ! open debug file on each process

  use module_function_decl
  implicit none
  FILENAME :: name

  write(name, '(i4.4)') my_pe()
  name = 'diag.' // name

  open(UDIAG, file = name)

  write(UDIAG,*) 'Output from PE ', my_pe()
  write(UDIAG,*) '~~~~~~~~~~~~~~~~~~~'
  write(UDIAG,*)

end

!-------------------------------------------------------------------------------
subroutine pos_keyword(unit, keyword)  ! positions unit at keyword

  implicit none
  integer, intent(in)               :: unit
  character(len = *), intent(in)    :: keyword 
  integer                           :: iostat
  character(len = len(keyword) + 1) :: word

  iostat = 0
  do while (iostat == 0)
     read(unit, *, iostat = iostat) word
     if (word == keyword) then
        backspace(unit)
        return
     endif
  enddo

  call die("pos_keyword(): " // keyword // ": not found")

end

!-------------------------------------------------------------------------------
subroutine read_keyword_int(unit, keyword, int, dim)  ! read integer(s) at keyw.

  implicit none
  integer, intent(in)               :: unit, dim
  character(len = *), intent(in)    :: keyword 
  integer, intent(out)              :: int(dim)
  character                         :: c

  call pos_keyword(unit, keyword)
  read(unit, *) c, int

end

!-------------------------------------------------------------------------------
subroutine read_keyword_REAL(unit, keyword, x, dim)  ! read float(s) at keyword

  implicit none
  integer, intent(in)               :: unit, dim
  character(len = *), intent(in)    :: keyword 
  REAL, intent(out)                 :: x(dim)
  character                         :: c

  call pos_keyword(unit, keyword)
  read(unit, *) c, x

end

!-------------------------------------------------------------------------------
subroutine swap_endian8(n, a)

  implicit none
  integer, parameter         :: i8 = 8
  integer, intent(in)        :: n
  integer                    :: i
  integer(i8), intent(inout) :: a(n)
  integer(i8), save          :: mask1, mask2, mask3, mask4, &
                                mask5, mask6, mask7, mask8
  integer(i8)                :: tmp
  logical, external          :: is_big_endian

  data mask1 /z'00000000000000FF'/
  data mask2 /z'000000000000FF00'/
  data mask3 /z'0000000000FF0000'/
  data mask4 /z'00000000FF000000'/
  data mask5 /z'000000FF00000000'/
  data mask6 /z'0000FF0000000000'/
  data mask7 /z'00FF000000000000'/
  data mask8 /z'FF00000000000000'/


  if (is_big_endian()) return

  do i = 1, n
     tmp = 0_i8
     tmp = ior(tmp, ishft(iand(a(i), mask1), 56_i8))
     tmp = ior(tmp, ishft(iand(a(i), mask2), 40_i8))
     tmp = ior(tmp, ishft(iand(a(i), mask3), 24_i8))
     tmp = ior(tmp, ishft(iand(a(i), mask4),  8_i8))
     tmp = ior(tmp, ishft(iand(a(i), mask5), -8_i8))
     tmp = ior(tmp, ishft(iand(a(i), mask6),-24_i8))
     tmp = ior(tmp, ishft(iand(a(i), mask7),-40_i8))
     tmp = ior(tmp, ishft(iand(a(i), mask8),-56_i8))
     a(i) = tmp
  enddo

end

!===============================================================================
