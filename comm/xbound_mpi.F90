!===============================================================================
!
! xbound_mpi.F90 - boundary exchange with MPI
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
module module_xbound
  
  implicit none
  integer, parameter :: max_bound = 3 * 3 * 3 * 3

  type type_xbound
     integer :: i_source
     integer :: i_target
     integer :: pe_source
     integer :: pe_target
     integer :: size        ! total size
     integer :: vector_type
     integer :: block_count
     integer :: block_size
     integer :: block_stride
  end type type_xbound
end

!-------------------------------------------------------------------------------
module module_xbound_g

  use module_xbound
  implicit none
  integer, save            :: n_bound = 0
  type (type_xbound), save :: b(max_bound)
end

!-------------------------------------------------------------------------------
module module_xbound_sc

  use module_xbound
  implicit none
  integer, save            :: n_bound = 0
  integer, save            :: i_bound(DIM, FWD:BWD)
  type (type_xbound), save :: b(max_bound)
end

!-------------------------------------------------------------------------------
module module_xbound_sc2

  use module_xbound
  implicit none
  integer, save            :: n_bound = 0
  integer, save            :: i_bound(DIM, FWD:BWD)
  type (type_xbound), save :: b(max_bound)
end

!-------------------------------------------------------------------------------
subroutine init_xbound()

  implicit none
  integer, external :: version_of_d

  call init_xbound_g()
  call init_xbound_sc()
  call init_xbound_sc2()
end

!-------------------------------------------------------------------------------
subroutine init_xbound_g()

  use module_xbound_g
  implicit none

  integer :: x, y, z, t

  do t = -1, 1
  do z = -1, 1
  do y = -1, 1
  do x = -1, 1
     call init_xch_bound(n_bound, b, NCOL * NCOL * SIZE_COMPLEX, x, y, z, t)
  enddo
  enddo
  enddo
  enddo
end

!-------------------------------------------------------------------------------
subroutine init_xbound_sc()

  use module_xbound_sc
  use module_lattice
  implicit none

  integer :: mu, block


  block = NDIRAC * NCOL * SIZE_COMPLEX

  call init_xch_bound(n_bound, b, block, -1,0,0,0);  i_bound(1, BWD) = n_bound
  call init_xch_bound(n_bound, b, block, +1,0,0,0);  i_bound(1, FWD) = n_bound
  call init_xch_bound(n_bound, b, block, 0,-1,0,0);  i_bound(2, BWD) = n_bound
  call init_xch_bound(n_bound, b, block, 0,+1,0,0);  i_bound(2, FWD) = n_bound
  call init_xch_bound(n_bound, b, block, 0,0,-1,0);  i_bound(3, BWD) = n_bound
  call init_xch_bound(n_bound, b, block, 0,0,+1,0);  i_bound(3, FWD) = n_bound
  call init_xch_bound(n_bound, b, block, 0,0,0,-1);  i_bound(4, BWD) = n_bound
  call init_xch_bound(n_bound, b, block, 0,0,0,+1);  i_bound(4, FWD) = n_bound

  do mu = 1, DIM
     if (npe(mu) == 1) then
        i_bound(mu, FWD) = 0
        i_bound(mu, BWD) = 0
     endif
  enddo
end

!-------------------------------------------------------------------------------
subroutine init_xbound_sc2()

  use module_xbound_sc2
  use module_lattice
  implicit none

  integer :: mu, block


  block = 2 * NCOL * SIZE_COMPLEX

  call init_xch_bound(n_bound, b, block, -1,0,0,0);  i_bound(1, BWD) = n_bound
  call init_xch_bound(n_bound, b, block, +1,0,0,0);  i_bound(1, FWD) = n_bound
  call init_xch_bound(n_bound, b, block, 0,-1,0,0);  i_bound(2, BWD) = n_bound
  call init_xch_bound(n_bound, b, block, 0,+1,0,0);  i_bound(2, FWD) = n_bound
  call init_xch_bound(n_bound, b, block, 0,0,-1,0);  i_bound(3, BWD) = n_bound
  call init_xch_bound(n_bound, b, block, 0,0,+1,0);  i_bound(3, FWD) = n_bound
  call init_xch_bound(n_bound, b, block, 0,0,0,-1);  i_bound(4, BWD) = n_bound
  call init_xch_bound(n_bound, b, block, 0,0,0,+1);  i_bound(4, FWD) = n_bound

  do mu = 1, DIM
     if (npe(mu) == 1) then
        i_bound(mu, FWD) = 0
        i_bound(mu, BWD) = 0
     endif
  enddo
end

!-------------------------------------------------------------------------------
subroutine init_xch_bound(n_bound, b, block_size, xx, yy, zz, tt)

  use      module_xbound
  use      module_function_decl
  use      module_nnpe
  use      module_offset
  use      module_lattice
  use      module_vol
  implicit none
  include 'mpif.h'

  integer, intent(inout)            :: n_bound
  type (type_xbound), intent(inout) :: b(max_bound)
  integer, intent(in)               :: block_size, xx, yy, zz, tt

  integer, dimension (DIM) :: dir, m, i, target, source
  integer, external        :: xyzt2i, n_sites, e_o
  integer                  :: x, y, z, t, size, mu, stride, block_count, ierror
  integer                  :: tmp_type1, tmp_type2, the_type
  integer(MPI_ADDRESS_KIND):: true_lb, true_extent
  integer                  :: extent

  logical                  :: special


  if (nnpe(xx, yy, zz, tt) == my_pe()) return

  if (xx /= 0 .and. yy == 0 .and. zz /= 0 .and. tt == 0) then
     special = .true.
  else
     special = .false.
  endif


  dir = (/ xx, yy, zz, tt /)

  do mu = 1, DIM
     if (dir(mu) /= 0) then
        m(mu) = 1
     else
        m(mu) = NH(mu)
     endif

     if (dir(mu) == -1) then
        target(mu) = -1
        source(mu) = N(mu) - 1
     elseif (dir(mu) == +1) then
        target(mu) = N(mu)
        source(mu) = 0
     else
        target(mu) = 0
        source(mu) = 0
     endif
  enddo


  size = block_size
  do mu = 1, DIM
     if (dir(mu) == 0) then
        size = size * NH(mu)
        m(mu) = 1
     else
        exit
     endif
  enddo

  stride = block_size
  do mu = 1, DIM
     if (m(mu) == 1) then
        stride = stride * NH(mu)
     else
        exit
     endif
  enddo

  block_count = 1
  do mu = 1, DIM
     block_count = block_count * m(mu)
  enddo

  n_bound = n_bound + 1
  ASSERT(n_bound <= max_bound)

  if (special) then  ! (y,t)-plane

     ! MPY-type for y-line:

     block_count = NY
     size = block_size
     stride = block_size * NXH

     call mpi_type_vector(block_count, size, stride, BQCD_REAL, tmp_type1, ierror)
     call mpi_type_commit(tmp_type1, ierror)


#ifdef ALTIX
     ! use MPI-1
     call mpi_type_extent(BQCD_REAL, extent, ierror)
     call mpi_type_struct(2, (/1, 1/), (/0, extent/), (/tmp_type1, MPI_UB/), &
                          tmp_type2, ierror) 
#else
     ! use MPI-2
     call mpi_type_get_true_extent(BQCD_REAL, true_lb, true_extent, ierror)
     call mpi_type_create_resized(tmp_type1, true_lb, true_extent, tmp_type2, ierror)
#endif
     call mpi_type_commit(tmp_type2, ierror)

     ! MPI-parameters for (y,t)-plane:

     block_count = NT
     size = 1
     stride = block_size * NXH * NY * NZ
     the_type = tmp_type2
     b(n_bound)%size = block_size * NY * NT

  else

     the_type = BQCD_REAL
     b(n_bound)%size = block_count * size

  endif


  b(n_bound)%i_source  = xyzt2i(source)
  b(n_bound)%i_target  = xyzt2i(target)
  b(n_bound)%pe_source = nnpe(xx, yy, zz, tt)
  b(n_bound)%pe_target = nnpe(-xx, -yy, -zz, -tt)

  b(n_bound)%block_count = block_count
  b(n_bound)%block_size  = size
  b(n_bound)%block_stride= stride

  call mpi_type_vector(block_count, size, stride, the_type, &
                       b(n_bound)%vector_type, ierror)
  call mpi_type_commit(b(n_bound)%vector_type, ierror)

  !!if ( my_pe() == 0) write(6,*) xx,yy,zz,tt, block_count, size, stride
  !!if ( my_pe() == 0) write(6,*) xx,yy,zz,tt, b(n_bound)%i_source, b(n_bound)%i_target, nnpe(xx,yy,zz,tt), my_pe()

  !!ASSERT(b(n_bound)%size == block_size * n_sites(DIM, dir, NH, NPE))

  if (special) then
     call mpi_type_free(tmp_type1, ierror)
     call mpi_type_free(tmp_type2, ierror)
  endif
end

!-------------------------------------------------------------------------------
subroutine xbound_g_field(u)

  use module_function_decl
  use module_vol
  implicit none

  GAUGE_FIELD :: u
  integer :: mu, eo, x, y, z, t

  if (num_pes() == 1) return

  do mu = 1, DIM
  do eo = EVEN, ODD
     call xbound_g(u, eo, mu)
  enddo
  enddo

end

!-------------------------------------------------------------------------------
subroutine xbound_g(u, eo, mu)

  use module_xbound_g
  use module_function_decl
  use module_vol
  implicit none
  include 'mpif.h'

  integer :: eo, mu, i, status(MPI_STATUS_SIZE), ierror
  GAUGE_FIELD :: u

  if (num_pes() == 1) return

  do i = 1, n_bound
     call mpi_sendrecv( &
         u(1,1, b(i)%i_source, eo,mu), 1, b(i)%vector_type,  b(i)%pe_target, 0,&
         u(1,1, b(i)%i_target, eo,mu), b(i)%size, BQCD_REAL, b(i)%pe_source, 0,&
         MPI_COMM_WORLD, status, ierror)
  enddo
end

!-------------------------------------------------------------------------------
subroutine xbound_sc_field(a)

  use module_function_decl
  use module_vol
  implicit none
  include 'mpif.h'

  integer :: i, status(MPI_STATUS_SIZE), ierror
  integer :: mu, fb
  SPINCOL_FIELD :: a

  if (num_pes() == 1) return

  do mu = 1, DIM
     call xbound_sc(a, mu)
  enddo
end

!-------------------------------------------------------------------------------
subroutine xbound_sc(a, direction)

  use module_xbound_sc
  use module_function_decl
  use module_vol
  implicit none
  include 'mpif.h'

  integer :: i, status(MPI_STATUS_SIZE), ierror
  integer :: direction, mu, fb
  SPINCOL_FIELD :: a

  if (num_pes() == 1) return

  mu = direction
  do fb = FWD, BWD
     i = i_bound(mu, fb)
     if (i /= 0) then
        call mpi_sendrecv( &
         a(1,1, b(i)%i_source), 1, b(i)%vector_type,  b(i)%pe_target, 0,&
         a(1,1, b(i)%i_target), b(i)%size, BQCD_REAL, b(i)%pe_source, 0,&
         MPI_COMM_WORLD, status, ierror)
     endif
  enddo
end

!-------------------------------------------------------------------------------
subroutine xbound_sc2_field(a)

  use module_function_decl
  use module_vol
  implicit none
  include 'mpif.h'

  SC2_FIELD :: a
  integer :: mu

  if (num_pes() == 1) return

  do mu = 1, DIM
     call xbound_sc2(a, mu)
  enddo
end

!-------------------------------------------------------------------------------
subroutine xbound_sc2_field_i(a)  ! "i"mmediate MPI calls

  use module_xbound_sc2
  use module_function_decl
  use module_vol
  implicit none
  include 'mpif.h'

  SC2_FIELD :: a

  integer, parameter :: max_request = 2 * 2 * DIM

  integer :: request(max_request), status(MPI_STATUS_SIZE, max_request), ierror
  integer :: mu, fb, i, n_request

  if (num_pes() == 1) return

  n_request = 0

  do mu = 1, DIM
  do fb = FWD, BWD
     i = i_bound(mu, fb)
     if (i /= 0) then
        n_request = n_request + 1
        call mpi_irecv( &
         a(1,1, b(i)%i_target, mu,fb), b(i)%size, BQCD_REAL, b(i)%pe_source, 0,&
         MPI_COMM_WORLD, request(n_request), ierror)
     endif
  enddo
  enddo

  do mu = 1, DIM
  do fb = FWD, BWD
     i = i_bound(mu, fb)
     if (i /= 0) then
        n_request = n_request + 1
        call mpi_isend( &
         a(1,1, b(i)%i_source, mu,fb), 1, b(i)%vector_type,  b(i)%pe_target, 0,&
         MPI_COMM_WORLD, request(n_request), ierror)
     endif
  enddo
  enddo

  call mpi_waitall(n_request, request, status, ierror)
end

!-------------------------------------------------------------------------------
subroutine xbound_sc2(a, direction)

  use module_xbound_sc2
  use module_function_decl
  use module_vol
  implicit none
  include 'mpif.h'

  integer :: i, status(MPI_STATUS_SIZE), ierror
  integer :: direction, mu, fb
  SC2_FIELD :: a

  if (num_pes() == 1) return

  mu = direction
  do fb = FWD, BWD
     i = i_bound(mu, fb)
     if (i /= 0) then
        call mpi_sendrecv( &
         a(1,1, b(i)%i_source, mu,fb), 1, b(i)%vector_type,  b(i)%pe_target, 0,&
         a(1,1, b(i)%i_target, mu,fb), b(i)%size, BQCD_REAL, b(i)%pe_source, 0,&
         MPI_COMM_WORLD, status, ierror)
     endif
  enddo
end

!===============================================================================
