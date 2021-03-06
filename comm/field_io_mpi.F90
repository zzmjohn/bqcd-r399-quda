!===============================================================================
!
! field_io_mpi.F90 - I/O routine for gauge and pseudo fermion fields using MPI
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

#ifndef INCLUDE_MPIF_H
#define INCLUDE_MPIF_H include 'mpif.h'
#endif

!-------------------------------------------------------------------------------
subroutine field_io(action, m, mx, field, cksum)  ! read or write g- or sc-field

  use typedef_cksum
  use module_decomp
  use module_function_decl
  use module_lattice_io
  use module_vol
  implicit none
  INCLUDE_MPIF_H

  character(len = *) :: action
  integer            :: m, mx
  type(type_cksum)   :: cksum(0:LT - 1)

  REAL :: field(SIZE_COMPLEX * m, 0:mx - 1, 0:NY - 1, 0:NZ - 1, 0:NT - 1)

  REAL :: buffer(0:(SIZE_COMPLEX * m * mx * NPE(1))-1, 0:(NY * NPE(2))-1)

  FILENAME           :: file
  integer            :: i_pe(DIM)
  integer, external  :: i_global, ilex
  integer            :: x, y, z, t, t_global, me, pe, size, rec, recl
  integer            :: size_field
  integer            :: pe_x, pe_y, pe_z, pe_t
  CHECK_SUM          :: check_sum(2), n_bytes
  integer            :: tt, pe_tt, pe_io
  integer            :: status(MPI_STATUS_SIZE), ierror
  integer            :: count, block_length, stride, buf_type
  logical            :: io_pe


  count        = NY
  block_length = SIZE_COMPLEX * m * mx
  stride       = block_length * NPE(1)

  size         = block_length * NY               ! words in send/recv
  size_field   = block_length * NY * NZ * NT     ! words in "field"
  n_bytes      = size * NPE(1) * NPE(2) * RKIND  ! bytes in "buffer"
  recl         = n_bytes                         ! cast to standard integer

  ASSERT(mod(recl, RECL_UNIT) == 0)
  recl         = recl / RECL_UNIT

  i_pe         = decomp%std%i_pe

  call mpi_type_vector(count, block_length, stride, BQCD_REAL, buf_type, ierror)
  call mpi_type_commit(buf_type, ierror)
  
  if (action == "write") call swap_endian8(size_field, field)

  if (i_pe(1) == 0 .and. i_pe(2) == 0 .and. i_pe(3) == 0) then
     io_pe = .true.
  else
     io_pe = .false.
  endif

  pe_t = i_pe(4)
  do t = 0, NT - 1
     t_global = i_global(t, NT, i_pe(4))
     
     file = cksum(t_global)%file
     
     if (io_pe) then
        open(UCONF, file = file, action = action, form = "unformatted", &
             access = "direct", recl = recl)
     endif
     
     rec = 0
     call cksum_init()
     
     do pe_z = 0, NPE(3) - 1
        do z = 0, NZ - 1
           rec = rec + 1
           
           if (io_pe .and. action == "read") then
              read(UCONF, rec = rec) buffer
              call cksum_add(buffer, n_bytes)
           endif
           
           do pe_y = 0, NPE(2) - 1
           do pe_x = 0, NPE(1) - 1

              y = count * pe_y
              x = block_length * pe_x

              call field_io_pes(pe, pe_io, (/pe_x, pe_y, pe_z, pe_t/))

              if (io_pe) then
                 if (my_pe() /= pe_io) call die("my_pe() /= pe_io")
              endif      

              if (my_pe() == pe .and. my_pe() == pe_io) then
                 call field_io_seq(action, count, block_length, stride, &
                                   field(1,0,0,z,t), buffer(x,y))
              else
                 if (action == "read") then
                    if (my_pe() == pe_io) then
                       call mpi_ssend(buffer(x,y), 1, buf_type, &
                                      pe, 0, MPI_COMM_WORLD, ierror)
                    endif
                    if (my_pe() == pe) then
                       call mpi_recv(field(1,0,0,z,t), size, BQCD_REAL, &
                                     pe_io, 0, MPI_COMM_WORLD, status, ierror)
                    endif
                 else
                    if (my_pe() == pe_io) then
                       call mpi_recv(buffer(x,y), 1, buf_type, &
                                     pe, 0, MPI_COMM_WORLD, status, ierror)
                    endif
                    if (my_pe() == pe) then
                       call mpi_ssend(field(1,0,0,z,t), size, BQCD_REAL, &
                                      pe_io, 0, MPI_COMM_WORLD, ierror)
                    endif
                 endif
              endif
              
           enddo
           enddo

           if (io_pe .and. action == "write") then
              write(UCONF, rec = rec) buffer
              call cksum_add(buffer, n_bytes)
           endif
        enddo
     enddo

     if (io_pe) then
        close(UCONF)
        call cksum_get(check_sum(1), check_sum(2))

        if (action == "read") then

              if (check_sum(1) /= cksum(t_global)%sum) then
                 call die("field_io(): check sum error in file " // file)
              endif

        else
           
           if (my_pe() == 0) then
              cksum(t_global)%sum = check_sum(1)
              cksum(t_global)%bytes = check_sum(2)
              do pe_tt = 1, NPE(4) - 1
                 tt = i_global(t, NT, pe_tt)
                 call mpi_recv(check_sum, 2, BQCD_CHECK_SUM, &
                      MPI_ANY_SOURCE, tt, MPI_COMM_WORLD, status, ierror)
                 cksum(tt)%sum = check_sum(1)
                 cksum(tt)%bytes = check_sum(2)
              enddo
           else
              call mpi_ssend(check_sum, 2, BQCD_CHECK_SUM, 0, t_global, &
                   MPI_COMM_WORLD, ierror)
           endif
        endif
     endif
  enddo

  call swap_endian8(size_field, field)
  call mpi_type_free(buf_type, ierror)
end

!-------------------------------------------------------------------------------
subroutine field_io_seq(action, count, block_length, stride, field, buffer)

  use module_lattice_io
  use module_vol
  implicit none

  character(len = *) :: action
  integer            :: count, block_length, stride
  REAL               :: field(*)
  REAL               :: buffer(*)
  integer            :: i, j, x, y

  i = 0
  j = 0
  do y = 1, count
     do x = 1, block_length
        i = i + 1
        if (action == "read") then
           field(i) = buffer(j + x)
        else
           buffer(j + x) = field(i)
        endif
     enddo
     j = j + stride
  enddo

end

!-------------------------------------------------------------------------------
subroutine field_io_pes(pe, pe_io, x_std)

  use module_lattice  ! in contrast to the calling routine !!
  implicit none
  integer, intent(out) :: pe, pe_io
  integer, intent(in)  :: x_std(DIM)
  integer              :: x_act(DIM), x_std_io(DIM), x_act_io(DIM)
  integer, external    :: ilex

  x_std_io(1) = 0
  x_std_io(2) = 0
  x_std_io(3) = 0
  x_std_io(4) = x_std(4)

  x_act(1) = x_std(gamma_index(1))
  x_act(2) = x_std(gamma_index(2))
  x_act(3) = x_std(gamma_index(3))
  x_act(4) = x_std(gamma_index(4))

  x_act_io(1) = x_std_io(gamma_index(1))
  x_act_io(2) = x_std_io(gamma_index(2))
  x_act_io(3) = x_std_io(gamma_index(3))
  x_act_io(4) = x_std_io(gamma_index(4))

  pe    = ilex(DIM, x_act,    NPE)
  pe_io = ilex(DIM, x_act_io, NPE)

end

!===============================================================================
