!===============================================================================
!
! allocate_shmem.F90 - allocation of gauge and pseudo fermion fields using shmem
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
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with BQCD.  If not, see <http://www.gnu.org/licenses/>.
!
!-------------------------------------------------------------------------------
# include "defs.h"
# include "shmem.h"

!-------------------------------------------------------------------------------
subroutine allocate_g_field(u)

  use module_vol
  implicit none
  P_GAUGE_FIELD :: u

  GAUGE_FIELD :: uu
  pointer (p_uu, uu)

  integer :: ierr

  if (associated(u)) call die("allocate_g_field(): memory leak")

  call barrier()
  call shpalloc(p_uu, SIZE_COMPLEX * NCOL * NCOL * volh_tot * 2 * DIM, ierr, 1)
  call cray_pointer_to_f90_pointer(uu)

CONTAINS

  subroutine cray_pointer_to_f90_pointer(uu)
  
    implicit none
    GAUGE_FIELD, target :: uu

    u => uu
  end subroutine cray_pointer_to_f90_pointer
  
end

!-------------------------------------------------------------------------------
subroutine allocate_g_field_io(u)

  use module_lattice_io
  use module_vol
  implicit none
  P_GAUGE_FIELD_IO :: u

  GAUGE_FIELD_IO :: uu
  pointer (p_uu, uu)
  
  integer :: ierr

  if (associated(u)) call die("allocate_g_field_io(): memory leak")

  call barrier()
  call shpalloc(p_uu, SIZE_COMPLEX * NCOL * (NCOL-1) * DIM * vol, ierr, 1)
  call cray_pointer_to_f90_pointer(uu)

CONTAINS

  subroutine cray_pointer_to_f90_pointer(uu)
  
    implicit none
    GAUGE_FIELD_IO, target :: uu

    u => uu
  end subroutine cray_pointer_to_f90_pointer
  
end

!-------------------------------------------------------------------------------
subroutine allocate_gen_field(u)

  use module_vol
  implicit none
  P_GENERATOR_FIELD :: u

  GENERATOR_FIELD :: uu
  pointer (p_uu, uu)
  
  integer :: ierr

  if (associated(u)) call die("allocate_gen_field(): memory leak")

  call barrier()
  call shpalloc(p_uu, NGEN * volh_tot * 2 * DIM, ierr, 1)
  call cray_pointer_to_f90_pointer(uu)

CONTAINS

  subroutine cray_pointer_to_f90_pointer(uu)
  
    implicit none
    GENERATOR_FIELD, target :: uu

    u => uu
  end subroutine cray_pointer_to_f90_pointer
  
end

!-------------------------------------------------------------------------------
subroutine allocate_sc_field(x)

  use module_vol
  implicit none
  P_SPINCOL_FIELD :: x

  SPINCOL_FIELD :: xx
  pointer (p_xx, xx)
  
  integer :: ierr

  if (associated(x)) call die("allocate_sc_field(): memory leak")

  call barrier()
  call shpalloc(p_xx, SIZE_COMPLEX * NDIRAC * NCOL * volh_tot, ierr, 1)
  call cray_pointer_to_f90_pointer(xx)

CONTAINS

  subroutine cray_pointer_to_f90_pointer(xx)
  
    implicit none
    SPINCOL_FIELD, target :: xx

    x => xx
  end subroutine cray_pointer_to_f90_pointer
  
end

!-------------------------------------------------------------------------------
subroutine allocate_sc2_field(x)

  use module_vol
  implicit none
  P_SC2_FIELD :: x

  SC2_FIELD :: xx
  pointer (p_xx, xx)
  
  integer :: ierr

  if (associated(x)) call die("allocate_sc2_field(): memory leak")

  call barrier()
  call shpalloc(p_xx, SIZE_COMPLEX * 2 * NCOL * volh_tot * DIM * 2, ierr, 1)
  call cray_pointer_to_f90_pointer(xx)

CONTAINS

  subroutine cray_pointer_to_f90_pointer(xx)
  
    implicit none
    SC2_FIELD, target :: xx

    x => xx
  end subroutine cray_pointer_to_f90_pointer
  
end

!-------------------------------------------------------------------------------
subroutine allocate_sc_field_io(x)

  use module_lattice_io
  use module_vol
  implicit none
  P_SPINCOL_FIELD_IO :: x

  SPINCOL_FIELD_IO :: xx
  pointer (p_xx, xx)
  
  integer :: ierr

  if (associated(x)) call die("allocate_sc_field_io(): memory leak")

  call barrier()
  call shpalloc(p_xx, SIZE_COMPLEX * NDIRAC * NCOL * volh, ierr, 1)
  call cray_pointer_to_f90_pointer(xx)

CONTAINS

  subroutine cray_pointer_to_f90_pointer(xx)
  
    implicit none
    SPINCOL_FIELD_IO, target :: xx

    x => xx
  end subroutine cray_pointer_to_f90_pointer
  
end

!-------------------------------------------------------------------------------
subroutine allocate_sc_overindexed(x)

  use module_vol
  implicit none
  P_SPINCOL_OVERINDEXED :: x

  SPINCOL_OVERINDEXED :: xx
  pointer (p_xx, xx)
  
  integer :: ierr

  if (associated(x)) call die("allocate_sc_overindexed(): memory leak")

  call barrier()
  call shpalloc(p_xx, SIZE_COMPLEX * NDIRAC * NCOL * volh_tot, ierr, 1)
  call cray_pointer_to_f90_pointer(xx)

CONTAINS

  subroutine cray_pointer_to_f90_pointer(xx)
  
    implicit none
    SPINCOL_OVERINDEXED, target :: xx

    x => xx
  end subroutine cray_pointer_to_f90_pointer
  
end

!===============================================================================
