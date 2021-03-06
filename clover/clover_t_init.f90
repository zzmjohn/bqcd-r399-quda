!===============================================================================
!
! clover_t_init.F90 - calculates T
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
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with BQCD. If not, see <http://www.gnu.org/licenses/>.
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
subroutine clover_t_init(t, b)

  use typedef_clover
  use module_vol
  implicit none

  complex(8), dimension(4, 3, 4, 3, volh) :: t
  type(type_clover_b) :: b(2, volh)

  complex(8), dimension (4, 3, volh_tot) :: x
!dir$ cache_align x

  integer :: c1, c2, s1, s2, i

  do c2 = 1, 3
     do s2 = 1, 4

        x = 0.0_8
        !$omp parallel do
        do i = 1, volh
           x(s2, c2, i) = 1.0_8
        enddo

        call clover_mult_b(b, x, volh)

        !$omp parallel do private(c1, s1)
        do i = 1, volh
           do c1 = 1, 3
              do s1 = 1, 4
                 t(s1, c1, s2, c2, i) = x(s1, c1, i)
              enddo
           enddo
        enddo

     enddo
  enddo

end

!===============================================================================
