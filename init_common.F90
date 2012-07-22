!===============================================================================
!
! init_common.F90 - initialize common blocks (now: mostly modules)
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
subroutine init_common(para)

  use typedef_para
  implicit none
  type(type_para), intent(in) :: para

  call init_common_lattice(para%L, para%NPE, para%bc_fermions, para%gamma_index)
  call init_common_vol(para%L, para%NPE)
  call init_common_nnpe(para%NPE)
  call init_common_offset()
  call init_common_nn()
  call init_common_thread()
end

!-------------------------------------------------------------------------------
subroutine init_common_lattice(ll, pe, bc_f, gamma_i)

  use module_function_decl
  use module_lattice
  implicit none

  integer, dimension(DIM), intent(in) :: ll, pe, bc_f, gamma_i

  integer :: i, count(DIM)

  L   = ll
  NPE = pe
  N   = ll / pe
  NH  = N

  NH(1) = N(1) / 2

  bc_fermions = bc_f
  gamma_index = gamma_i

  if (mod(L(1), 2) /= 0) call die("init_common_lattice(): LX must be even")
  if (mod(L(2), 2) /= 0) call die("init_common_lattice(): LY must be even")
  if (mod(L(3), 2) /= 0) call die("init_common_lattice(): LZ must be even")
  if (mod(L(4), 2) /= 0) call die("init_common_lattice(): LT must be even")

  if (mod(N(1), 2) /= 0) call die("init_common_lattice(): NX must be even")

  if (NPE(1) * NPE(2) * NPE(3) * NPE(4) /= num_pes()) then
     call die("init_common_lattice(): N_PEs wrong")
  endif

if (N(1)*NPE(1)/=L(1)) call die("init_common_lattice(): NPEX not divider of LX")
if (N(2)*NPE(2)/=L(2)) call die("init_common_lattice(): NPEY not divider of LY")
if (N(3)*NPE(3)/=L(3)) call die("init_common_lattice(): NPEZ not divider of LZ")
if (N(4)*NPE(4)/=L(4)) call die("init_common_lattice(): NPET not divider of LT")

  count = 0
  do i = 1, DIM
     if (gamma_index(i) < 1 .or. gamma_index(i) > DIM) then
        call die("init_common_lattice(): gamma_index: out of range")
     endif
     count(gamma_index(i)) = count(gamma_index(i)) + 1
  enddo

  do i = 1, DIM
     if (count(i) /= 1) then
        call die("init_common_lattice(): gamma_index: inconsistent")
     endif
  enddo

  select case (version_of_d())
     case(3,4,22)
        do i = 1, DIM
          if (gamma_index(i) /= i) then
             call die( &
  "init_common_lattice(): gamma_index: not changeable for this version of D()")
          endif
        enddo
  end select

  do i = 1, DIM
     decomp_direction(gamma_index(i)) = i
  enddo

end

!-------------------------------------------------------------------------------
subroutine init_common_nn()

  use module_function_decl
  use module_lattice
  use module_vol
  use module_nn
  implicit none

  integer :: x, y, z, t, i, eo, mu, fb, dir, tmp
  integer, dimension (DIM) :: j, start, end
  integer, external :: e_o, xyzt2i, i_periodic
  integer, parameter :: out_of_range = 2000000000

  allocate(nn(volh_tot, EVEN:ODD, DIM, FWD:BWD))

  do fb = FWD, BWD
     do mu = 1, DIM
        do eo = EVEN, ODD
           !$omp parallel do
           do i = 1, volh_tot
              nn(i, eo, mu, fb) = out_of_range
           enddo
        enddo
     enddo
  enddo

  do mu = 1, DIM
     if (NPE(mu) == 1) then
        start(mu) = 0
        end(mu) = N(mu) - 1
     else
        start(mu) = -1
        end(mu) = N(mu)
     endif
  enddo

  do t = start(4), end(4)
  do z = start(3), end(3)
  do y = start(2), end(2)
  do x = start(1), end(1)

     j = (/x,y,z,t/)

     i = xyzt2i(j)
     eo = e_o(j)

     do fb = FWD, BWD
        if (fb == FWD) then
           dir = +1
        else
           dir = -1
        endif
        do mu = 1, DIM
           j = (/x, y, z, t/)
           j(mu) = j(mu) + dir

           if (NPE(mu) == 1) j(mu) = i_periodic(j(mu), L(mu))

           if (j(mu) < -1 .or. j(mu) > N(mu)) then
              nn(i, eo, gamma_index(mu), fb) = out_of_range
           else
              nn(i, eo, gamma_index(mu), fb) = xyzt2i(j)
           endif
        enddo
     enddo

  enddo
  enddo
  enddo
  enddo

  do fb = FWD, BWD
     do mu = 1, DIM
        do eo = EVEN, ODD
           do i = 1, volh
              tmp = nn(i, eo, mu, fb)
       if (tmp < 1 .or. tmp > volh_tot) call die("init_common_nn(): error1")
       if (num_pes() == 1 .and. tmp > volh) call die("init_common_nn(): error2")
           enddo
        enddo
     enddo
  enddo

  do fb = FWD, BWD
     do mu = 1, DIM
        do eo = EVEN, ODD
           do i = 1, volh
              tmp = nn(i, eo, mu, fb)
              if (nn(tmp, EVEN + ODD - eo, mu, FWD + BWD - fb) /= i) &
                 call die("init_common_nn(): error3")
           enddo
        enddo
     enddo
  enddo

end

!-------------------------------------------------------------------------------
subroutine init_common_nnpe(NPE)

  use module_function_decl
  use module_nnpe
  implicit none

  integer, dimension (DIM) :: NPE, i, j, start, end
  integer, external :: i_periodic, ilex
  integer :: me, pe, x, y, z, t, mu

  me = my_pe()
  nnpe(:, :, :, :) = me

  call unlex(me, DIM, i, NPE)

  do mu = 1, DIM
     if (NPE(mu) == 1) then
        start(mu) = 0
        end(mu) = 0
     else
        start(mu) = -1
        end(mu) = 1
     endif
  enddo

  do t = start(4), end(4)
  do z = start(3), end(3)
  do y = start(2), end(2)
  do x = start(1), end(1)

     j(1) = i_periodic(i(1) + x, NPE(1))
     j(2) = i_periodic(i(2) + y, NPE(2))
     j(3) = i_periodic(i(3) + z, NPE(3))
     j(4) = i_periodic(i(4) + t, NPE(4))

     pe = ilex(DIM, j, NPE)

     if (pe < 0 .or. pe >= num_pes()) then
        call die ("init_common_pe(): pe out of range")
     endif

     nnpe(x, y, z, t) = pe

  enddo
  enddo
  enddo
  enddo

  ASSERT(nnpe(0,0,0,0) == my_pe())

end

!-------------------------------------------------------------------------------
subroutine init_common_offset()

  use module_lattice
  use module_vol
  use module_offset
  implicit none

  integer, external :: n_sites
  integer :: x, y, z, t, off, off2, mu
  integer :: start(DIM), end(DIM)


  !!ASSERT(n_sites(DIM, (/0,0,0,0/), N, NPE) == vol)
  !!ASSERT(n_sites(DIM, (/0,0,0,0/), NH, NPE) == volh)

  do mu = 1, DIM
     if (NPE(mu) == 1) then
        start(mu) = 0
        end(mu) = 0
     else
        start(mu) = -1
        end(mu) = 1
     endif
  enddo

  off = n_sites(DIM,(/0,0,0,0/), NH, NPE)  ! volh

  do t = start(4), end(4)
  do z = start(3), end(3)
  do y = start(2), end(2)
  do x = start(1), end(1)

     if (x == 0 .and. y == 0 .and. z == 0 .and. t == 0) then
        offset(x,y,z,t) = 0
     else
        offset(x,y,z,t) = off
        off = off + n_sites(DIM, (/x,y,z,t/), NH, NPE)
     endif

  enddo
  enddo
  enddo
  enddo

  off2 = 1
  do mu = 1, DIM
     if (NPE(mu) == 1) then
        off2 = off2 * NH(mu)
     else
        off2 = off2 * (NH(mu) + 2)
     endif
  enddo

  ASSERT(off == off2)
  ASSERT(off <= volh_tot)

end

!-------------------------------------------------------------------------------
subroutine init_common_vol(L, NPE)

  use module_vol
  implicit none

  integer, dimension(DIM), intent(in) :: L, NPE
  integer, dimension(DIM) :: N

  N = L / NPE  

  volume   = L(1) * L(2) * L(3) * L(4)
  vol      = N(1) * N(2) * N(3) * N(4)
  volh     = vol / 2
  volh_tot = (N(1)/2 + 2) * (N(2) + 2) * (N(3) + 2) * (N(4) + 2)

  size_sc_field = SIZE_COMPLEX * NDIRAC * NCOL * volh

end

!-------------------------------------------------------------------------------
subroutine init_common_thread()

  use module_function_decl
  use module_thread
  use module_vol
  implicit none

  integer :: i, size

#ifdef _OPENMP
  integer :: omp_get_max_threads
  n_thread = omp_get_max_threads()
#else
  n_thread = 1
#endif

  if (version_of_d() /= 21 .and. &
      version_of_d() /= 22 .and. &
      version_of_d() > 2) then
     if (n_thread < 2) then
        call die("init_common_thread(): " // &
                 "need at least 2 threads for compiled version of D()") 
     endif
  endif

  allocate(xyz_start(0:n_thread-1))  
  allocate(xyz_end(0:n_thread-1))  
  allocate(t_start(0:n_thread-1))  
  allocate(t_end(0:n_thread-1))  

  if (n_thread == 1) then
     xyz_start(0) = 1
     xyz_end(0) = volh
     t_start(0) = 1
     t_end(0) = volh
  else
     xyz_start(0) = 0
     xyz_end(0) = -1
     call init_common_thread_split(n_thread - 1, volh, xyz_start(1), xyz_end(1))
     call init_common_thread_split(n_thread,     volh,   t_start,      t_end)
  endif
end

!-------------------------------------------------------------------------------
subroutine init_common_thread_split(n, size, start, end)

  implicit none
  integer, intent(in) :: n, size
  integer, intent(out) :: start(n), end(n)
  integer :: chunk, rest, i

  chunk = size / n
  rest = size - (chunk * n)

  start(1) = 1
  end(1) = chunk
  if (rest >= 1) end(1) = end(1) + 1

  do i = 2, n
     start(i) = end(i - 1) + 1
     end(i)   = end(i - 1) + chunk
     if (rest >= i) end(i) = end(i) + 1
  enddo

  ASSERT(end(n) == size)
end

!===============================================================================
