!===============================================================================
!
! xyzt2i.F90 - maps local coordinates (x,y,z,t) to even/odd index
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
integer function xyzt2i(x_in)
  
  !   x_in := (x,y,z,t) 
  !  -1 <= x(mu) <= N(mu) ; mu = 1,2,3,4
  !  xyzt2i >= 1

  use module_function_decl
  use module_lattice
  use module_offset
  implicit none

  integer, dimension (DIM), intent(in) :: x_in
  integer, dimension (DIM)             :: dir, i, m, x
  integer                              :: count, mu
  integer, external                    :: ieo, n_sites, i_periodic, ilex


  count = 0
  do mu = 1, DIM

     x(mu) = x_in(mu)

     if (x(mu) < -1 .or. x(mu) > N(mu)) then
        call die('xyzt2i(): x(mu) out of range')
     endif

     if (NPE(mu) == 1) x(mu) = i_periodic(x(mu), N(mu))

     if (x(mu) == -1) then
        dir(mu) = -1
     elseif (x(mu) == N(mu)) then
        dir(mu) = 1
     else
        dir(mu) = 0
        count = count + 1
     endif

     if (dir(mu) /= 0) then
        i(mu) = 0
        m(mu) = 1
     else
        i(mu) = x(mu)
        m(mu) = N(mu)
     endif
  enddo
  
  if (count == DIM) then
     xyzt2i = offset(0,0,0,0) + ieo(DIM, x, N) + 1
  else
     ASSERT(num_pes() /= 1)

     if (dir(1) /= 0) then
         !!ASSERT(n_sites(DIM, dir, N, NPE) == n_sites(DIM, dir, NH, NPE))
         !!ASSERT(ilex(DIM, i, m) <= n_sites(DIM, dir, NH, NPE))
         xyzt2i = offset(dir(1),dir(2),dir(3),dir(4)) + ilex(DIM, i, m) + 1
     else
         xyzt2i = offset(dir(1),dir(2),dir(3),dir(4)) + ieo(DIM, i, m) + 1
     endif
  endif

end

!-------------------------------------------------------------------------------
integer function std_xyzt2i(x)
  
  use module_lattice
  implicit none

  integer, dimension (DIM), intent(in) :: x
  integer, dimension (DIM)             :: x_act
  integer, external                    :: xyzt2i

  x_act(1) = x(gamma_index(1))
  x_act(2) = x(gamma_index(2))
  x_act(3) = x(gamma_index(3))
  x_act(4) = x(gamma_index(4))

  std_xyzt2i = xyzt2i(x_act)
end

!===============================================================================
