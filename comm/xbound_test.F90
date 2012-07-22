!===============================================================================
!
! xbound_test.F90 - test of xbound_g() 
!                   all possible dimensions must be decomposed
!
!-------------------------------------------------------------------------------
!
! Copyright (C) 2006 Hinnerk Stueben
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
subroutine xbound_test()

  use module_function_decl
  use module_lattice
  use module_vol
  implicit none
  GAUGE_FIELD :: u, v
  integer     :: x, y, z, t, j(4), i, eo, global(4)
  integer     :: is_bound_x, is_bound_y, is_bound_z, is_bound_t 
  integer     :: dx, dy, dz, dt
  integer, external :: xyzt2i, e_o
  character(16) :: status

!!  call conf_zero(u)
!!
!!  do i = 1, volh
!!     u(1, 1, i, EVEN, 1) = 123
!!     u(1, 1, i, ODD, 1) = 789
!!  enddo
!!
!!  call xbound_g(u, EVEN, 1)
!!  call xbound_g(u, ODD, 1)
!!
!!  do i = 1, volh_tot
!!     ASSERT(u(1, 1, i, EVEN, 1) == 123)     
!!     ASSERT(u(1, 1, i, ODD, 1) == 789)     
!!  enddo
!!
!!!-----------------------------------------------
  call conf_zero(u)

  u = cmplx(12345.0, 67890.0)

  call open_diag()

  do t = 0, NT - 1
     do z = 0, NZ - 1
        do y = 0, NY - 1
           do x = 0, NX - 1
              j = (/x, y, z, t/)

              i = xyzt2i(j)
              eo = e_o(j)

              call local2global(my_pe(), j, global)

              u(1, 1, i, eo, 1) = global(1)
              u(2, 2, i, eo, 1) = global(2)
              u(3, 3, i, eo, 1) = global(3)
              u(1, 2, i, eo, 1) = global(4)

              !!write(UDIAG, "(4i6, 2i8)") j, i, eo

              !!write(UDIAG, "(10i6)") j, global, i, eo
           enddo
        enddo
     enddo
  enddo

  call xbound_g(u, EVEN, 1)
  call xbound_g(u, ODD, 1)

  write(UDIAG,*)
  write(UDIAG,*)

  !!ASSERT(e_o((/0,0,0,0/)) == 0)


  do i = 1, volh_tot
     do eo = EVEN, ODD
              write(UDIAG, "(a,2i6,4f8.1)") "alles: ", i, eo, &
                   real(u(1, 1, i, eo, 1)), &
                   real(u(2, 2, i, eo, 1)), &
                   real(u(3, 3, i, eo, 1)), &
                   real(u(1, 2, i, eo, 1))
     enddo
  enddo

  write(UDIAG,*)
  write(UDIAG,*)

  do t = -1, NT
     is_bound_t = 0
     if (t == -1) is_bound_t = 1
     if (t == NT) is_bound_t = 1
     do z = -1, NZ
        is_bound_z = 0
        if (z == -1) is_bound_z = 1
        if (z == NZ) is_bound_z = 1
        do y = -1, NY
           is_bound_y = 0
           if (y == -1) is_bound_y = 1
           if (y == NY) is_bound_y = 1
           do x = -1, NX
              is_bound_x = 0
              if (x == -1) is_bound_x = 1
              if (x == NX) is_bound_x = 1

              if (is_bound_x + is_bound_y + is_bound_z + is_bound_t <= 2) then


!!           do x = NX, NX

              j = (/x, y, z, t/)
              j = (/x, y, z, t/)

              i = xyzt2i(j)
              eo = e_o(j)

              call local2global(my_pe(), j, global)

              dx = -is_bound_x; if (x == NX) dx = 1
              dy = -is_bound_y; if (y == NY) dy = 1
              dz = -is_bound_z; if (z == NZ) dz = 1
              dt = -is_bound_t; if (t == NT) dt = 1

              if (u(1, 1, i, eo, 1) == global(1) .and. &
                  u(2, 2, i, eo, 1) == global(2) .and. &
                  u(3, 3, i, eo, 1) == global(3) .and. &
                  u(1, 2, i, eo, 1) == global(4)) then
                 status = "  okay"
              !!elseif (u(1, 1, i, eo, 1) == global(1) + dx .and. &
              !!        u(2, 2, i, eo, 1) == global(2) + dy .and. &
              !!        u(3, 3, i, eo, 1) == global(3) + dz .and. &
              !!        u(1, 2, i, eo, 1) == global(4) + dt) then
              !!   status = " okay2"
              else
                 status = ""

                 dx = int(u(1, 1, i, eo, 1)) - global(1)
                 dy = int(u(2, 2, i, eo, 1)) - global(2)
                 dz = int(u(3, 3, i, eo, 1)) - global(3)
                 dt = int(u(1, 2, i, eo, 1)) - global(4)

                 write(status,"(a,4i3,a)") "  (", dx, dy, dz, dt, ")"
              endif

              
              !!if (eo == 0) write(UDIAG, "(4i6, 2x, 4i6, i8, 2i3)") j, i, eo
              !write(UDIAG, "(10i6)") j, global, i, eo
              !!ASSERT(eo == mod(4+x+y+z+t, 2))

              write(UDIAG, "(10i6,4f8.1,a)") j, global, i, eo, &
                   real(u(1, 1, i, eo, 1)), &
                   real(u(2, 2, i, eo, 1)), &
                   real(u(3, 3, i, eo, 1)), &
                   real(u(1, 2, i, eo, 1)), status

              ASSERT(u(1, 1, i, eo, 1) == global(1))
              ASSERT(u(2, 2, i, eo, 1) == global(2))
              ASSERT(u(3, 3, i, eo, 1) == global(3))
              ASSERT(u(1, 2, i, eo, 1) == global(4))

           endif
           enddo
        enddo
     enddo
  enddo

end
