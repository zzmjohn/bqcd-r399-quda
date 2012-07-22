!===============================================================================
!
! conf_info.F90 - read/write/check file containing configuration parameters
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

!-------------------------------------------------------------------------------
subroutine write_conf_info_header(para, i_ensemble1, i_ensemble2, plaq)

  use module_conf_info
  use typedef_para
  use module_bqcd
  use module_counter
  use module_decomp
  use module_function_decl
  implicit none

  type(type_para) :: para
  REAL            :: plaq
  integer         :: i_ensemble1, i_ensemble2
  integer         :: i, e(2)

  e(1) = i_ensemble1
  e(2) = i_ensemble2

  if (my_pe() == 0) then
     call begin(UINFO, "ConfInfoHeader")
     write(UINFO,   *) k_format, conf_info_version
     write(UINFO, 400) k_prog, prog_name, prog_version
     write(UINFO,   *) k_run,  para%run
     write(UINFO,   *) k_traj, counter%traj
     write(UINFO, 405) k_host, rechner()
     write(UINFO, 400) k_date, datum(), uhrzeit()
     write(UINFO, 410) k_L, decomp%std%L
     write(UINFO, 410) k_bc, decomp%std%bc_fermions
     write(UINFO,   *) k_rkind, RKIND
     write(UINFO, 420) k_plaq, plaq

     do i = 1, 2
        write(UINFO,   *) trim(k_ensemble(i)), e(i)
        write(UINFO, 405) trim(k_beta(i)),     trim(para%c_hmc(e(i))%beta)
        write(UINFO, 405) trim(k_kappa(i)),    trim(para%c_hmc(e(i))%kappa)
        write(UINFO, 405) trim(k_csw(i)),      trim(para%c_hmc(e(i))%csw)
        write(UINFO, 405) trim(k_csw_kappa(i)),trim(para%c_hmc(e(i))%csw_kappa)
        write(UINFO, 405) trim(k_h(i)),        trim(para%c_hmc(e(i))%h)
     enddo
     call end(UINFO, "ConfInfoHeader")
  endif

400  format (3(1x,a))
405  format (2(1x,a))
410  format (1x,a,4i3)
420  format (1x,a,1x,e25.14)
end

!-------------------------------------------------------------------------------
subroutine read_conf_info_header(unit, info)

  use module_bqcd
  use module_conf_info
  implicit none

  type(type_conf_info) :: info
  integer              :: unit, v, i
  integer, external    :: pos_keyword

  call read_keyword_int(unit, k_format, v, 1)

  if (v /= conf_info_version) then
     call die("read_conf_info_header(): wrong file format")
  endif

  call read_keyword_int(unit, k_L, info%L, DIM)
  call read_keyword_int(unit, k_bc, info%bc_fermions, DIM)
  call read_keyword_int(unit, k_rkind, info%rkind, 1)

  do i = 1, 2
     call read_keyword_int (unit, k_ensemble(i), info%ensemble(i), 1)
     call read_keyword_REAL(unit, k_beta(i),     info%beta(i),     1)
     call read_keyword_REAL(unit, k_kappa(i),    info%kappa(i),    1)
     call read_keyword_REAL(unit, k_csw(i),      info%csw(i),      1)
     call read_keyword_REAL(unit, k_csw_kappa(i),info%csw_kappa(i),1)
     call read_keyword_REAL(unit, k_h(i),        info%h(i),        1)
  enddo

end  

!-------------------------------------------------------------------------------
subroutine check_conf_info_header(restart, info, para)

  use module_conf_info
  use module_decomp
  use typedef_para
  implicit none

  logical              :: restart
  type(type_conf_info) :: info
  type(type_para)      :: para
  integer              :: mu, i

  if (info%rkind /= RKIND) call die("check_conf_info_header(): RKIND wrong")

  do mu = 1, DIM
     if (info%L(mu) /= decomp%std%L(mu)) then
        call die("check_conf_info_header(): L inconsistent")
     endif
  enddo

  if (restart) then
     do mu = 1, DIM
        if (info%bc_fermions(mu) /= decomp%std%bc_fermions(mu)) then
           call die("check_conf_info_header(): bc_fermions inconsistent")
        endif
     enddo
     
     do i = 1, 2
        if (info%ensemble(i) < 1 .or. info%ensemble(i) > para%n_temper) then
           call die("check_conf_info_header(): i_ensemble out of range")
        endif
        
        if (info%beta(i) /= para%hmc(info%ensemble(i))%beta) then
           call die("check_conf_info_header(): beta inconsistent")
        endif
        
        if (info%kappa(i) /= para%hmc(info%ensemble(i))%kappa) then
           call die("check_conf_info_header(): kappa inconsistent")
        endif

        if (info%csw(i) /= para%hmc(info%ensemble(i))%csw) then
           call die("check_conf_info_header(): csw inconsistent")
        endif

        if (abs(info%csw_kappa(i) - &
                para%hmc(info%ensemble(i))%csw_kappa) > 1e-13 ) then
           call die("check_conf_info_header(): csw_kappa inconsistent")
        endif

        if (info%h(i) /= para%hmc(info%ensemble(i))%h) then
           call die("check_conf_info_header(): h inconsistent")
        endif
     enddo
  endif
     
end

!-------------------------------------------------------------------------------
subroutine read_cksum(restart, unit, cksum, LT, i_temper, file_name)

  use typedef_cksum
  implicit none
  logical                             :: restart
  integer                             :: unit, LT, i_temper, t
  FILENAME, external                  :: file_name
  type(type_cksum), dimension(0:LT-1) :: cksum

  call pos_keyword(unit, ">BeginCheckSum")
  read(unit,*)

  do t = 0, LT - 1
     read(unit,*) cksum(t)%file, cksum(t)%sum

     if (restart) then
        if (cksum(t)%file /= file_name(i_temper, t)) then
           call die("read_cksum(): file names inconsistent")
        endif
     endif
  enddo

end

!-------------------------------------------------------------------------------
subroutine write_cksum(unit, cksum, LT)

  use typedef_cksum
  use module_function_decl
  implicit none

  integer                         :: unit, LT, i
  type(type_cksum), dimension(LT) :: cksum
  
  if (my_pe() == 0) then
     call begin(unit, "CheckSum")
     do i = 1, LT
        write(unit, *) trim(cksum(i)%file), cksum(i)%sum, cksum(i)%bytes
     enddo
     call end(unit, "CheckSum")
  endif

end

!===============================================================================
