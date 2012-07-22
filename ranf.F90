!===============================================================================
!
! ranf.F90 - Fortran90 implementation of the Cray random number generator ranf()
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
module module_ranf

  integer, parameter :: ikind = 8
  integer, parameter :: jkind = 8
  integer, parameter :: rkind = 8

  integer(ikind), parameter :: mH = 2651554_ikind
  integer(ikind), parameter :: mL = 15184245_ikind
  integer(ikind), parameter :: mask24 = 16777215_ikind
  integer(ikind), parameter :: mask48 = 281474976710655_ikind
  integer(ikind), parameter :: default_seed = 48131768981101_ikind
  integer(ikind), save      :: seed = default_seed

!  integer(ikind), save :: mH
!  integer(ikind), save :: mL
!  integer(ikind), save :: mask24
!  integer(ikind), save :: mask48
!  integer(ikind), save :: default_seed
!  integer(ikind), save :: seed
!
!  data mH /o"12072642"/
!  data mL /o"71730565"/
!  data mask24 /o"77777777"/
!  data mask48 /o"7777777777777777"/
!  data default_seed /o"1274321477413155"/
!  data seed /o"1274321477413155"/

end

!-------------------------------------------------------------------------------
function ranf()

  use module_ranf
  implicit none
  real(rkind)    :: ranf
  integer(ikind) :: seedH, seedL

  ! seed = mod(m * seed, 48):

  seedH = iand(mask24, ishft(seed, -24_jkind))
  seedL = iand(mask24, seed)
  seed  = iand(mask48, seedL * mL + ishft(seedL * mH + seedH * mL, 24_jkind))

  ! normalize result:

  ranf = real(seed, kind = rkind)
  ranf = set_exponent(fraction(ranf), exponent(ranf) - 48_jkind)

end

!-------------------------------------------------------------------------------
subroutine ranget(seed_out)

  use module_ranf
  implicit none
  integer(ikind), intent(out) :: seed_out

  seed_out = seed

end

!-------------------------------------------------------------------------------
subroutine ranset(seed_in, n_skip)

  use module_ranf
  implicit none
  integer(ikind), intent(in) :: seed_in, n_skip
  integer(ikind) :: n, mm, mmH, mmL, seedH, seedL

  if (seed_in == 0_ikind) then
     seed = default_seed
  else
     seed = iand(mask48, ibset(seed_in, 0_jkind))
  endif

  ! skip "n_skip" seeds [i.e. calculate seed = mod(m**n_skip * seed, 48)]:

  n = iand(mask48, n_skip)
  mmH = mH
  mmL = mL

  do while (n > 0_ikind)

     if (btest(n, 0_ikind)) then  ! seed = seed * mm
        seedH = iand(mask24, ishft(seed, -24_jkind))
        seedL = iand(mask24, seed)
        seed  = iand(mask48, seedL*mmL + ishft(seedL*mmH + seedH*mmL, 24_jkind))
     endif

     mm  = iand(mask48, mmL * mmL + ishft(mmH * mmL, 25_jkind))  ! mm = mm * mm
     mmH = iand(mask24, ishft(mm, -24_jkind))
     mmL = iand(mask24, mm)
     
     n = ishft(n, -1_jkind)
  enddo

end

!===============================================================================
