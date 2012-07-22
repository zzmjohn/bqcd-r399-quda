!===============================================================================
!
! files.F90
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
!
! restart_file:      progname.run.{res|count|ran|stop}
! restart_conf_file: progname.run.s.time.{u|phi}
! info_file:         progname.run.s.info
!
! conf_info_file:    progname.run.s1.s2.traj.info
! conf_file:         progname.run.s1.s2.traj.time.u
!
!-------------------------------------------------------------------------------
# include "defs.h"

!-------------------------------------------------------------------------------
module module_files   ! formats of file name strings

  implicit none

  ! lengths of formats of file name strings

  integer, parameter :: l_name = 2
  integer, parameter :: l_ext = 2
  integer, parameter :: l_sep = 5
  integer, parameter :: l_num = 4

  integer, parameter :: l_base = l_name + l_sep + l_num
  integer, parameter :: l_conf = l_base + 3 * (l_sep + l_num)

  ! formats of file name strings

  character(len = l_name), save :: fmt_name     = "(a"
  character(len = l_name), save :: fmt_ext      = "a)"
  character(len = l_sep), save  :: fmt_sep      = ",'.',"
  character(len = l_num), save  :: fmt_run      = "i3.3"
  character(len = l_num), save  :: fmt_ensemble = "i1.1"
  character(len = l_num), save  :: fmt_traj     = "i5.5"
  character(len = l_num), save  :: fmt_time     = "i2.2"
  
  character(len = *), parameter :: ext_count = "count"
  character(len = *), parameter :: ext_ran   = "ran"
  character(len = *), parameter :: ext_res   = "res"
  character(len = *), parameter :: ext_u     = "u"
  character(len = *), parameter :: ext_phi   = "phi"
  character(len = *), parameter :: ext_info  = "info"
  character(len = *), parameter :: ext_stop  = "STOP"

CONTAINS

  character(len=l_base) function fmt_base()
    fmt_base = fmt_name // fmt_sep // fmt_run
  end function fmt_base

  character(len=l_conf) function fmt_conf()
    fmt_conf = fmt_base() // fmt_sep // fmt_ensemble &
                          // fmt_sep // fmt_ensemble &
                          // fmt_sep // fmt_traj
  end function fmt_conf

end

!-------------------------------------------------------------------------------
FILENAME function count_file()

  use module_files
  implicit none
  FILENAME :: restart_file

  count_file = restart_file(ext_count)
end

!-------------------------------------------------------------------------------
FILENAME function ran_file()

  use module_files
  implicit none
  FILENAME :: restart_file

  ran_file = restart_file(ext_ran)
end

!-------------------------------------------------------------------------------
FILENAME function res_file()

  use module_files
  implicit none
  FILENAME :: restart_file

  res_file = restart_file(ext_res)
end

!-------------------------------------------------------------------------------
FILENAME function stop_file()

  use module_files
  implicit none
  FILENAME :: restart_file

  stop_file = restart_file(ext_stop)
end

!-------------------------------------------------------------------------------
FILENAME function u_file(i_ensemble, time)

  use module_files
  implicit none
  integer, intent(in) :: i_ensemble, time
  FILENAME :: restart_conf_file

  u_file = restart_conf_file(i_ensemble, time, ext_u)
end

!-------------------------------------------------------------------------------
FILENAME function phi_file(i_ensemble, time)

  use module_files
  implicit none
  integer, intent(in) :: i_ensemble, time
  FILENAME :: restart_conf_file

  phi_file = restart_conf_file(i_ensemble, time, ext_phi)
end


!-------------------------------------------------------------------------------
FILENAME function restart_file(ext)

  use module_bqcd
  use module_counter
  use module_files
  implicit none

  character(len = *), intent(in) :: ext
  FILENAME_FORMAT                :: fmt

  fmt = fmt_base() // fmt_sep // fmt_ext
  
  write(restart_file, fmt) prog_name, counter%run, ext

end

!-------------------------------------------------------------------------------
FILENAME function restart_conf_file(i_ensemble, time, ext)

  use module_bqcd
  use module_counter
  use module_files
  implicit none

  integer, intent(in)            :: i_ensemble, time
  character(len = *), intent(in) :: ext
  FILENAME_FORMAT                :: fmt

  fmt = fmt_base() // fmt_sep // fmt_ensemble &
                   // fmt_sep // fmt_time &
                   // fmt_sep // fmt_ext
  
  write(restart_conf_file, fmt) prog_name, counter%run, i_ensemble, time, ext

end

!-------------------------------------------------------------------------------
FILENAME function info_file(i_ensemble)

  use module_bqcd
  use module_counter
  use module_files
  implicit none

  integer, intent(in)            :: i_ensemble
  FILENAME_FORMAT                :: fmt

  fmt = fmt_base() // fmt_sep // fmt_ensemble // fmt_sep // fmt_ext
  
  write(info_file, fmt) prog_name, counter%run, i_ensemble, ext_info

end

!-------------------------------------------------------------------------------
FILENAME function conf_info_file(i_ensemble1, i_ensemble2)

  use module_bqcd
  use module_counter
  use module_files
  implicit none

  integer, intent(in)            :: i_ensemble1, i_ensemble2
  FILENAME_FORMAT                :: fmt

  fmt = fmt_conf() // fmt_sep // fmt_ext
  
  write(conf_info_file, fmt) &
       prog_name, counter%run, i_ensemble1, i_ensemble2, counter%traj, ext_info

end

!-------------------------------------------------------------------------------
FILENAME function conf_file(i_ensemble1, i_ensemble2, time)

  use module_bqcd
  use module_counter
  use module_files
  implicit none

  integer, intent(in)            :: i_ensemble1, i_ensemble2, time
  FILENAME_FORMAT                :: fmt

  fmt = fmt_conf() // fmt_sep // fmt_time // fmt_sep // fmt_ext
  
  write(conf_file, fmt) &
     prog_name, counter%run, i_ensemble1, i_ensemble2, counter%traj, time, ext_u

end

!-------------------------------------------------------------------------------
subroutine check_fmt(run, max_temper, max_traj, max_time)

  use module_files
  implicit none
  integer :: run, max_temper, max_traj, max_time

  call check_len(run,        fmt_run,      "RUN")
  call check_len(max_temper, fmt_ensemble, "TEMPER")
  call check_len(max_traj,   fmt_traj,     "TRAJ")
  call check_len(max_time,   fmt_time,     "TIME")

CONTAINS

  subroutine check_len(counter, counter_fmt, counter_name)

    implicit none
    integer :: i, len, counter
    character(len = *) :: counter_fmt, counter_name

    i = index(counter_fmt, "i")

    if (i == 0) then
       call die("check_fmt(): unable to check fmt for " // counter_name)
    endif

    read(counter_fmt(i+1:i+1), *) len

    if (counter < 0 .or. counter >= 10**len) then
       call die("check_fmt(): file name format unsuitable for " // counter_name)
    endif

  end subroutine check_len

end

!-------------------------------------------------------------------------------
subroutine set_fmt_ensemble(N_temper)

  use module_files
  implicit none

  integer, intent(in) :: N_temper

  if (N_temper < 10) then
     fmt_ensemble = "i1.1"
  else if (N_temper < 100) then
     fmt_ensemble = "i2.2"
  else
     call die ("set_fmt_ensemble(): N_temper >= 100 ???")
  endif

end


!-------------------------------------------------------------------------------
function format_ensemble()

  use module_files
  implicit none
  character(len = l_num) :: format_ensemble

  format_ensemble = fmt_ensemble
end

!-------------------------------------------------------------------------------
subroutine filename_test()

  use module_function_decl
  implicit  none
  integer   i

  FILENAME, external :: count_file
  FILENAME, external :: ran_file
  FILENAME, external :: res_file
  FILENAME, external :: stop_file
  FILENAME, external :: u_file
  FILENAME, external :: phi_file
  FILENAME, external :: restart_conf_file
  FILENAME, external :: info_file
  FILENAME, external :: conf_info_file
  FILENAME, external :: conf_file

  do i = 1, 11, 10
     call set_fmt_ensemble(i)

     PUTVAL(6, count_file())
     PUTVAL(6, ran_file())
     PUTVAL(6, res_file())
     PUTVAL(6, stop_file())
     PUTVAL(6, u_file(3, 4))
     PUTVAL(6, phi_file(5, 6))
     PUTVAL(6, restart_conf_file(1, 2, 'conf'))
     PUTVAL(6, info_file(7))
     PUTVAL(6, conf_info_file(7, 8))
     PUTVAL(6, conf_file(3, 4, 5))
  enddo

end

!===============================================================================
