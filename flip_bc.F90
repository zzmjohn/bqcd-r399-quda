!===============================================================================
!
! flip_bc.F90 - flip fermionic boundary conditions
!               (ie multiplication of corresponding links with -1)
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
module module_flip_bc

  INTEGER, dimension(:, :, :), pointer, save :: flip_bc_list
  integer, dimension(DIM), save              :: flip_bc_len

end

!-------------------------------------------------------------------------------
subroutine init_flip_bc()

  use      module_flip_bc
  use      module_function_decl
  use      module_lattice
  use      module_vol
  implicit none

  integer, dimension(DIM) :: i0, i1, i_pe, j
  integer :: me, mu, x, y, z, t, i, eo, count(EVEN:ODD)
  integer, external :: xyzt2i, e_o

  allocate(flip_bc_list(volh_tot, EVEN:ODD, DIM))

  me = my_pe()
  call unlex(me, DIM, i_pe, NPE)

  do mu = 1, DIM
     count = 0
     if (bc_fermions(mu) < 0) then
        if (i_pe(mu) == (NPE(mu) - 1) .or. i_pe(mu) == 0) then
           i0 = 0
           i1 = N - 1

           if (i_pe(mu) == (NPE(mu) - 1)) then
              i0(mu) = N(mu) - 1
           else
              i0(mu) = -1
           endif
        
           i1(mu) = i0(mu)

           do t = i0(4), i1(4)
           do z = i0(3), i1(3)
           do y = i0(2), i1(2)
           do x = i0(1), i1(1)
              j = (/x, y, z, t/)
              i = xyzt2i(j)
              eo = e_o(j)

              count(eo) = count(eo) + 1
              flip_bc_list(count(eo), eo, mu) = i
           enddo
           enddo
           enddo
           enddo
             
        endif
     endif

     if (count(EVEN) /= count(ODD)) then
        call die ("init_flip_bc(): count(EVEN) /= count(ODD)")
     else
        flip_bc_len(mu) = count(EVEN)
     endif
  enddo

end

!-------------------------------------------------------------------------------
subroutine flip_bc(u)

  use      module_flip_bc
  use      module_lattice
  use      module_vol
  implicit none

  GAUGE_FIELD, intent(inout) :: u
  integer :: mu, nu, count, i, eo, c1, c2

  do mu = 1, DIM
     nu = gamma_index(mu)
     do eo = EVEN,ODD
        do count = 1, flip_bc_len(mu)
           i = flip_bc_list(count, eo, mu) 
           do c2 = 1, NCOL
              do c1 = 1, NCOL
                 u(c1, c2, i, eo, nu) = -u(c1, c2, i, eo, nu)
              enddo
           enddo
        enddo
     enddo
  enddo

end

!===============================================================================
