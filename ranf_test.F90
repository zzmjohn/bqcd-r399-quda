!===============================================================================
!
! ranf_test - test of "ranf.F90"
!
!-------------------------------------------------------------------------------
!
! Copyright (C) 1998-2001 Hinnerk Stueben
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

program ranf_test

  implicit none
  integer(8), parameter :: null = 0
  integer(8) :: i, seed
  real(8) :: x, ranf

  write(6,400) "default seed:"

  do i = 1, 10
     call ranget(seed)
     x = ranf()
     write(6,410) i, seed, x
  enddo

  write(6,400) "seed = 4711:"

  call ranset(4711_8, null)

  do i = 1, 10
     call ranget(seed)
     x = ranf()
     write(6,410) i, seed, x
  enddo

  write(6,400) "seed varies, no skip:"

  do i = -10, 20
     call ranset(i, null)
     call ranget(seed)
     x = ranf()
     write(6,410) i, seed, x
  enddo

  write(6,400) "default seed, skip varies:"

  do i = -10, 20
     call ranset(null, i)
     call ranget(seed)
     x = ranf()
     write(6,410) i, seed, x
  enddo

  write(6,400) "large seeds, no skip:"

  do i = 0, 47
     call ranset(2_8**i, null)
     call ranget(seed)
     x = ranf()
     write(6,410) i, seed, x
  enddo

  write(6,400) "default seeds, big skips:"

  do i = 0, 47
     call ranset(null, 2_8**i)
     call ranget(seed)
     x = ranf()
     write(6,410) i, seed, x
  enddo


400 format (//1x,a//)
410 format (i6,i24,f24.16)

end

!===============================================================================
