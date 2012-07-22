!===============================================================================
!
! ran.F90 - random number related routines
!
!-------------------------------------------------------------------------------
!
! Copyright (C) 1998-2003 Hinnerk Stueben
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
subroutine ran_gauss_volh(m, ran, var, eo)

! intitializes x with Gaussian random numbers in a manner that makes
! results independent of the lattice decomposition

  use module_decomp
  use module_function_decl
  use module_vol
  implicit none

  integer, intent(in)  :: m
  COMPLEX, intent(out) :: ran(m, volh)
  REAL,    intent(in)  :: var
  integer, intent(in)  :: eo

  integer, dimension (DIM) :: i_pe, block
  integer :: i, ii, x, y, z, t, j
  integer :: nxh, nx, ny, nz, nt, lx, ly, lz, npe(DIM)
  REAL    :: twovar, phi, r

  twovar = TWO * var

  i_pe = decomp%std%i_pe
  npe  = decomp%std%NPE

  nxh = decomp%std%NH(1)

  nx = decomp%std%N(1)
  ny = decomp%std%N(2)
  nz = decomp%std%N(3)
  nt = decomp%std%N(4)

  lx  = decomp%std%L(1)
  ly  = decomp%std%L(2)
  lz  = decomp%std%L(3)

  ! block() contains the number of random numbers to be skipped.
  ! In the calculation of block() a factor 1/2 from "even/odd volume"
  ! cancels with 2 from "two random numbers per site".

  block(1) = m * nx
  block(2) = m * lx * ny
  block(3) = m * lx * ly * nz
  block(4) = m * lx * ly * lz * nt

  i = 0
  call ranskip(block(4) * i_pe(4))
  do t = 0, nt - 1
     call ranskip(block(3) * i_pe(3))
     do z = 0, nz - 1
        call ranskip(block(2) * i_pe(2))
        do y = 0, ny - 1
           call ranskip(block(1) * i_pe(1))
           do x = 0, nxh - 1
              i = i + 1
              ii = decomp%act%i(i, eo)
              do j = 1, m
                 r = sqrt(-twovar * log(ranf()))
                 phi = TWOPI * ranf()
                 ran(j, ii) = cmplx(r * cos(phi), r * sin(phi), kind = RKIND) 
              enddo
           enddo
           call ranskip(block(1) * (npe(1) - 1 - i_pe(1)))
        enddo
        call ranskip(block(2) * (npe(2) - 1 - i_pe(2)))
     enddo
     call ranskip(block(3) * (npe(3) - 1 - i_pe(3)))
  enddo
  call ranskip(block(4) * (npe(4) - 1 - i_pe(4)))


CONTAINS

  subroutine ranskip(n)
    integer :: n
    SEED    :: seed, n_skip

    n_skip = n
    call ranget(seed)
    call ranset(seed, n_skip)
  end subroutine ranskip

end

!-------------------------------------------------------------------------------
subroutine rancheck()  ! checks if seed is identical on all PEs
 
  use module_function_decl
  implicit none
  SEED :: my_seed, seed

  call ranget(my_seed)

  if (my_pe() == 0) then
     seed = my_seed
  else 
     seed = -1
  endif

  call comm_broadcast_i8(seed, 1)

  if (seed /= my_seed) call die('rancheck(): seeds differ')

end

!-------------------------------------------------------------------------------
subroutine write_ran()  ! save state of random number generator

  use module_function_decl
  implicit none
  SEED  :: seed
  FILENAME, external :: ran_file

  call rancheck()
  call ranget(seed)
  if (my_pe() == 0) then
     open(URAN, file = ran_file(), action = "write")
     write(URAN, *) seed
     close(URAN)
  endif

end

!-------------------------------------------------------------------------------
subroutine random_sequence(r, n)  ! random permutation of [1..n]

  use module_function_decl
  implicit none
  integer, intent(in)  :: n
  integer, intent(out) :: r(n)
  integer              :: ran, seq(n), i, j, len_seq
  integer              :: ceiling

  len_seq = n

  do i = 1, len_seq
     seq(i) = i
  enddo

  do i = 1, n - 1
     ran = ceiling(len_seq * ranf())
     if (ran <= 0 .or. ran > len_seq) stop "random_sequence(): ran out of range"
     r(i) = seq(ran)
     do j = ran, len_seq - 1
        seq(j) = seq(j + 1)
     enddo
     len_seq = len_seq - 1
  enddo
  r(n) = seq(1)

end

!-------------------------------------------------------------------------------
subroutine get_a_random_seed(seed)  ! (try to) generate a random seed

  use module_function_decl
  implicit none
  SEED    :: seed
  integer :: count, rate, rate_10sec
  integer :: pe, ierror

  if (my_pe() == 0) then
     call system_clock(count = count, count_rate = rate)

     if (rate <= 0) call die("get_a_random_seed(): failed")

     rate_10sec = rate * 10
     seed = mod(count, rate_10sec)
  endif

  call seed_broadcast(seed)
end

!-------------------------------------------------------------------------------
subroutine init_ran(para, flags)  ! initilize random number generator

  use typedef_para
  use typedef_flags
  implicit none

  type(type_para)    :: para
  type(type_flags)   :: flags
  FILENAME, external :: ran_file
  SEED               :: seed, null

  if (flags%continuation_job) then
     open(URAN, file = ran_file(), action = "read", status = "old")
     read(URAN, *) para%seed
     close(URAN)
  else
     seed = para%seed
     if (seed < 0) call get_a_random_seed(para%seed)
  endif
  
  null = 0
  call ranset(para%seed, null)

end

!===============================================================================
