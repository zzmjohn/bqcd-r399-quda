!===============================================================================
!
! allocate.F90 - allocation of gauge and pseudo fermion fields
!
!-------------------------------------------------------------------------------
!
! Copyright (C) 1998-2005 Hinnerk Stueben
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
subroutine allocate_g_field(u)

  use module_vol
  implicit none
  complex(8), dimension(:, :, :, :, :), pointer :: u

  if (associated(u)) then
     call die("allocate_g_field(): memory leak")
  else
     allocate(u(3, 3, volh_tot, 0:1, 4))
     call conf_zero(u)
  endif
end

!-------------------------------------------------------------------------------
subroutine allocate_g_field_io(u)

  use module_lattice_io
  implicit none
  complex(8), dimension(:, :, :, :, :, :, :), pointer :: u

  if (associated(u)) then
     call die("allocate_g_field_io(): memory leak")
  else
     allocate(u(3, 3 -1, 4, 0:NX-1, 0:NY-1, 0:NZ-1, 0:NT-1))
  endif
end

!-------------------------------------------------------------------------------
subroutine allocate_gen_field(x)

  use module_vol
  implicit none
  real(8), dimension(:, :, :, :), pointer :: x

  integer :: i, eo, mu

  if (associated(x)) then
     call die("allocate_gen_field(): memory leak")
  else
     allocate(x(8, volh_tot, 0:1, 4))
     do mu = 1, 4
        do eo = 0, 1
           !$omp parallel do
           do i = 1, volh
              x(1, i, eo, mu) = 0.0_8
              x(2, i, eo, mu) = 0.0_8
              x(3, i, eo, mu) = 0.0_8
              x(4, i, eo, mu) = 0.0_8
              x(5, i, eo, mu) = 0.0_8
              x(6, i, eo, mu) = 0.0_8
              x(7, i, eo, mu) = 0.0_8
              x(8, i, eo, mu) = 0.0_8
           enddo
        enddo
     enddo
  endif
end

!-------------------------------------------------------------------------------
subroutine allocate_sc_field(x)

  use module_vol
  implicit none
  complex(8), dimension(:, :, :), pointer :: x

  if (associated(x)) then
     call die("allocate_sc_field(): memory leak")
  else
     allocate(x(4, 3, volh_tot))
     call sc_zero(x)
  endif
end

!-------------------------------------------------------------------------------
subroutine allocate_sc_field_io(x)

  use module_lattice_io
  implicit none
  complex(8), dimension(:, :, :, :, :, :), pointer :: x

  if (associated(x)) then
     call die("allocate_sc_field_io(): memory leak")
  else
     allocate(x(4, 3, 0:NXH-1, 0:NY-1, 0:NZ-1, 0:NT-1))
  endif
end

!-------------------------------------------------------------------------------
subroutine allocate_sc_overindexed(x)

  use module_vol
  implicit none
  real(8), dimension(:), pointer :: x

  if (associated(x)) then
     call die("allocate_sc_overindexed(): memory leak")
  else
     allocate(x(2*4*3*volh_tot))
  endif
end

!-------------------------------------------------------------------------------
subroutine allocate_sc2_field(x)

  use module_vol
  implicit none
  complex(8), dimension(:, :, :, :, :), pointer :: x

  if (associated(x)) then
     call die("allocate_sc2_field(): memory leak")
  else
     allocate(x(2, 3, volh_tot, 4, 0:1))
  endif
end

!===============================================================================
