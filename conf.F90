!===============================================================================
!
! conf.F90 - operations on gauge field and pseudo fermion field configurations
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
subroutine init_confs(para, conf)

  use typedef_para
  use module_p_interface
  use module_switches
  implicit none

  type(type_para)                       :: para
  type(hmc_conf), dimension(MAX_TEMPER) :: conf
  integer                               :: i

  do i = 1, para%n_temper
     call allocate_g_field(conf(i)%u)
     call allocate_sc_field(conf(i)%phi)
     if (switches%hasenbusch) call allocate_sc_field(conf(i)%phi2)
     if (para%start == START_HOT .or. para%start == START_COLD) then
        call init_u(conf(i)%u, para%start)
        conf(i)%former = i
     endif
  enddo

  if (para%start == START_CONT) call conf_read(.true., para, conf)
  if (para%start == START_FILE) call conf_read(.false., para, conf)

  do i = 1, para%n_temper
     if (para%hmc(i)%csw_kappa /= ZERO) then
        call allocate_clover_field_a(conf(i)%a)
        call allocate_clover_field_a(conf(i)%i)
        call allocate_clover_field_b(conf(i)%b)
        call clover_init(conf(i)%a, conf(i)%i, conf(i)%b, &
                         conf(i)%u, para%hmc(i)%csw_kappa)
     endif
  enddo

end

!-------------------------------------------------------------------------------
subroutine init_u(u, start)  ! initialization of u-field (at trajectory 0)

  use module_vol
  implicit none

  GAUGE_FIELD :: u
  integer :: start

  select case (start)
     case (START_HOT)
        call conf_hot(u)
     case (START_COLD)
        call conf_cold(u)
     case default
        call die("init_u(): don't know how to start")
  end select

  call xbound_g_field(u)
end

!-------------------------------------------------------------------------------
subroutine conf_check(u)  ! checks if u-field is SU(3)

  use module_vol
  implicit none

  GAUGE_FIELD, intent(in) :: u
  SU3 :: v
  SU3, parameter :: su3_one = reshape( &
                           (/ ONE,ZERO,ZERO, &
                              ZERO,ONE,ZERO, &
                              ZERO,ZERO,ONE /), &
                           (/ NCOL, NCOL /))
  REAL :: dev, d
  integer :: i, j, k, eo, mu


  dev = ZERO
  do mu = 1, DIM
     do eo = EVEN, ODD
        do k = 1, VOLH
           call su3_check_det(u(1, 1, k, eo, mu))
           call uud(v, u(1, 1, k, eo, mu), u(1, 1, k, eo, mu))
           do i = 1, NCOL
              do j = 1, NCOL
                 d = abs(Re(v(i, j)) - Re(su3_one(i, j))) & 
                   + abs(Im(v(i, j)) - Im(su3_one(i, j)))
              enddo
           enddo
           dev = max(dev, d)
        enddo
     enddo
  enddo

  if (dev > 1e-13) call die('conf_check(): dev > 1e-13')
!!write(0,'(x,a,e10.2)') 'conf_check(): max deviation is ', dev

end

!-------------------------------------------------------------------------------
subroutine conf_normalize(u)  ! normalizes u-field to SU(3)

  use module_vol
  implicit none

  GAUGE_FIELD, intent(inout) :: u
  integer :: i, eo, mu

  do mu = 1, DIM
     do eo = EVEN, ODD
        do i = 1, volh
           call u_normalize(u(1, 1, i, eo, mu))
        enddo
     enddo
  enddo

end

!-------------------------------------------------------------------------------
subroutine conf_zero(u)  ! init ("OpenMP first touch")

  use module_vol
  implicit none

  GAUGE_FIELD :: u
  integer :: i, eo, mu

  do mu = 1, DIM
     do eo = EVEN, ODD
        !$omp parallel do
        do i = 1, volh
           u(1, 1, i, eo, mu) = ZERO
           u(2, 1, i, eo, mu) = ZERO
           u(3, 1, i, eo, mu) = ZERO
           u(1, 2, i, eo, mu) = ZERO
           u(2, 2, i, eo, mu) = ZERO
           u(3, 2, i, eo, mu) = ZERO
           u(1, 3, i, eo, mu) = ZERO
           u(2, 3, i, eo, mu) = ZERO
           u(3, 3, i, eo, mu) = ZERO
        enddo
     enddo
  enddo

end

!-------------------------------------------------------------------------------
subroutine conf_cold(u)  ! cold start

  use module_vol
  implicit none

  GAUGE_FIELD :: u
  integer :: i, eo, mu

  do mu = 1, DIM
     do eo = EVEN, ODD
        do i = 1, volh
           u(1, 1, i, eo, mu) = ONE
           u(2, 1, i, eo, mu) = ZERO
           u(3, 1, i, eo, mu) = ZERO
           u(1, 2, i, eo, mu) = ZERO
           u(2, 2, i, eo, mu) = ONE
           u(3, 2, i, eo, mu) = ZERO
           u(1, 3, i, eo, mu) = ZERO
           u(2, 3, i, eo, mu) = ZERO
           u(3, 3, i, eo, mu) = ONE
        enddo
     enddo
  enddo

end

!-------------------------------------------------------------------------------
subroutine conf_hot(u)  ! hot start

  use module_vol
  implicit none

  GAUGE_FIELD :: u
  integer :: i, eo, mu

  do mu = 1, DIM
     do eo = EVEN, ODD
        call ran_gauss_volh(NCOL * NCOL, u(1, 1, 1, eo, mu), ONE, eo)
        do i = 1, volh
           call u_normalize(u(1, 1, i, eo, mu))
        enddo
     enddo
  enddo

end

!-------------------------------------------------------------------------------
subroutine conf_seq(action, u, u_io)  ! arranges u-field for i/o

  use module_lattice_io
  use module_decomp
  use module_vol
  implicit none

  character(*) :: action
  GAUGE_FIELD :: u
  GAUGE_FIELD_IO :: u_io

  integer, dimension(DIM) :: j
  integer :: x, y, z, t, i, eo, mu, c1, c2
  integer, external :: std_xyzt2i, e_o

  do t = 0, NT - 1
     do z = 0, NZ - 1
        do y = 0, NY - 1
           do x = 0, NX - 1

              j = (/x, y, z, t/)

              i = std_xyzt2i(j)
              eo = e_o(j)

              do mu = 1, DIM
                 do c2 = 1, NCOL - 1
                    do c1 = 1, NCOL
                       if (action == "read") then
                          u(c1, c2, i, eo, mu) = u_io(c1, c2, mu, x, y, z, t)
                       else
                          u_io(c1, c2, mu, x, y, z, t) = u(c1, c2, i, eo, mu)
                       endif
                    enddo
                 enddo
                 if (action == "read") call u_complete(u(1, 1, i, eo ,mu))
              enddo
           enddo
        enddo
     enddo
  enddo
end

!-------------------------------------------------------------------------------
subroutine phi_seq(action, phi, phi_io)  ! arranges phi-field for i/o

  use module_lattice_io
  use module_decomp
  use module_vol
  implicit none

  character(*) :: action
  SPINCOL_FIELD :: phi
  SPINCOL_FIELD_IO :: phi_io

  integer, dimension(DIM) :: j
  integer :: x, y, z, t, d, c, i
  integer, external :: std_xyzt2i, e_o

  do t = 0, NT - 1
     do z = 0, NZ - 1
        do y = 0, NY - 1
           do x = 0, NX - 1

              j = (/x, y, z, t/)

              if (e_o(j) == EVEN) then
                 i = std_xyzt2i(j)
                 do c = 1, NCOL
                    do d = 1, NDIRAC
                       if (action == "read") then
                          phi(d, c, i) = phi_io(d, c, x, y, z, t)
                       else
                          phi_io(d, c, x, y, z, t) = phi(d, c, i)
                       endif
                    enddo
                 enddo
              endif

           enddo
        enddo
     enddo
  enddo
end

!-------------------------------------------------------------------------------
subroutine conf_read(restart, para, conf)

  use typedef_cksum
  use typedef_para
  use module_conf_info
  use module_lattice_io
  use module_p_interface
  use module_switches
  use module_vol
  implicit none

  character(len = *), parameter              :: READ = "read"

  logical                                    :: restart
  type(type_para)                            :: para
  type(hmc_conf), dimension(MAX_TEMPER)      :: conf

  type(type_conf_info)                       :: info
  type(type_cksum), dimension(0:para%L(4)-1) :: cksum
  P_GAUGE_FIELD_IO, save                     :: u_io
  P_SPINCOL_FIELD_IO, save                   :: phi_io
  FILENAME, external                         :: u_file, phi_file, info_file
  FILENAME                                   :: file
  integer                                    :: i, t
  integer                                    :: u_m, u_mx, phi_m, phi_mx
  integer                                    :: n_u_io, n_phi
  
  ALLOCATE_G_FIELD_IO(u_io)
  ALLOCATE_SC_FIELD_IO(phi_io)

  u_m    = NCOL * (NCOL - 1) * DIM
  u_mx   = NX
  phi_m  = NDIRAC * NCOL
  phi_mx = NXH
  n_u_io = u_m * vol * SIZE_COMPLEX
  n_phi  = size_sc_field
  
  do i = 1, para%n_temper
     
     if (restart) then
        file = info_file(i)
     else
        file = para%info_file(i)
     endif

     open(UINFO, file = file, action = READ, status = "old")
     call read_conf_info_header(UINFO, info)
     call check_conf_info_header(restart, info, para)

     if (restart) then
        conf(i)%former = info%ensemble(2)
     else
        conf(i)%former = i
     endif

     ! read U

     call read_cksum(restart, UINFO, cksum, para%L(4), i, u_file)

     TIMING_START(timing_bin_u_read)
     call field_io(READ, u_m, u_mx, u_io, cksum)
     TIMING_STOP(timing_bin_u_read)

     call conf_seq(READ, conf(i)%u, u_io)
     call xbound_g_field(conf(i)%u)

     if (switches%tempering .and. switches%dynamical) then
        ! read PHI
        call read_cksum(restart, UINFO, cksum, para%L(4), i, phi_file)
        call field_io(READ, phi_m, phi_mx, phi_io, cksum)
        call phi_seq(READ, conf(i)%phi, phi_io)
        call xbound_sc_field(conf(i)%phi)
     endif

     close(UINFO)
  enddo
end

!-------------------------------------------------------------------------------
subroutine conf_write(restart, para, conf)

  use typedef_cksum
  use typedef_para
  use module_function_decl
  use module_lattice_io
  use module_p_interface
  use module_switches
  use module_vol
  implicit none

  character(len = *), parameter              :: WRITE = "write"

  logical                                    :: restart
  type(type_para)                            :: para
  type(hmc_conf), dimension(MAX_TEMPER)      :: conf

  type(type_cksum), dimension(0:para%L(4)-1) :: cksum
  P_GAUGE_FIELD_IO, save                     :: u_io
  P_SPINCOL_FIELD_IO, save                   :: phi_io
  FILENAME, external                         :: u_file, phi_file, info_file
  FILENAME, external                         :: conf_file, conf_info_file
  FILENAME                                   :: f_info
  integer                                    :: i, j, t
  integer                                    :: u_m, u_mx, phi_m, phi_mx
  integer                                    :: n_u_io, n_phi
  REAL                                       :: plaq
  REAL, external                             :: sg  


  ALLOCATE_G_FIELD_IO(u_io)
  ALLOCATE_SC_FIELD_IO(phi_io)

  u_m    = NCOL * (NCOL - 1) * DIM
  u_mx   = NX
  phi_m  = NDIRAC * NCOL
  phi_mx = NXH
  n_u_io = u_m * vol * SIZE_COMPLEX
  n_phi  = size_sc_field
  
  call check_former(para%n_temper, conf)

  do i = 1, para%n_temper
     
     j = conf(i)%former

     if (restart) then
        f_info = info_file(i)
     else
        f_info = conf_info_file(i, j)
     endif

     if (my_pe() == 0) open(UINFO, file = f_info, action = WRITE)
     plaq = sg(conf(i)%u) / (SIX * volume)
     call write_conf_info_header(para, i, j, plaq)

     ! write U

     do t = 0, para%L(4) - 1
        if (restart) then
           cksum(t)%file = u_file(i, t)
        else
           cksum(t)%file = conf_file(i, j, t)
        endif
     enddo

     call conf_seq(WRITE, conf(i)%u, u_io)

     TIMING_START(timing_bin_u_write)
     call field_io(WRITE, u_m, u_mx, u_io, cksum)
     TIMING_STOP(timing_bin_u_write)

     call write_cksum(UINFO, cksum, para%L(4))

     if (switches%tempering .and. switches%dynamical .and. restart) then
        ! write PHI
        do t = 0, para%L(4) - 1
           cksum(t)%file = phi_file(i, t)
        enddo
        call phi_seq(WRITE, conf(i)%phi, phi_io)
        call field_io(WRITE, phi_m, phi_mx, phi_io, cksum)
        call write_cksum(UINFO, cksum, para%L(4))
     endif

     if (my_pe() == 0) close(UINFO)
  enddo
end
    
!===============================================================================
