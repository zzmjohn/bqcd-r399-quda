!===============================================================================
!
! service.F90 - calls to service routines (now mostly standard Fortran)
!
!-------------------------------------------------------------------------------
!
! Copyright (C) 1998-2012 Hinnerk Stueben
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
subroutine abbruch() ! exit with error status
                     ! abort parallel program

  call comm_abort()
end

!-------------------------------------------------------------------------------
function rechner()   ! returns hostname

  character(len = 20) rechner
  character(len = 32) r

  call hostnm(r)
  rechner = r
end

!-------------------------------------------------------------------------------
SECONDS function sekunden()  ! returns wall clock seconds
                             ! SECONDS is defined in "defs.h"

#ifdef USE_MPI_WTIME

! use MPI_Wtime():

  include 'mpif.h'
  sekunden = mpi_wtime()

#else

! use Fortran90 routine (resolution might not be high enough)

  implicit none
  real(8)       :: count_to_seconds
  integer(8)    :: count, count_rate
  logical, save :: initialized = .false.

  if (.not. initialized) then
     call system_clock(count, count_rate)
     count_to_seconds = 1.0_8 / real(count_rate, kind = 8)
  endif

  call system_clock(count)
  
  sekunden = count * count_to_seconds

#endif

end

!-------------------------------------------------------------------------------
subroutine pxfgetarg(iarg, arg, length, status)

  implicit none
  integer :: iarg, length, status
  character(*) :: arg

  call get_command_argument(iarg, arg, length, status)
end

!-------------------------------------------------------------------------------
integer function ipxfargc()

  ipxfargc = command_argument_count()
end

!-------------------------------------------------------------------------------
logical function is_big_endian()

  implicit none
  integer(4), save :: i = ichar('0')  &
                        + ichar('1') * 256  &
                        + ichar('2') * 256**2  &
                        + ichar('3') * 256**3
  character(4) :: c
  equivalence (i, c)

  if (c == "0123") then
     is_big_endian = .false.
  else if (c == "3210") then
     is_big_endian = .true.
  else
     call die("is_big_endian(): unable to get endianness")
  endif

end

!===============================================================================
