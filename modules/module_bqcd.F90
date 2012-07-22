!===============================================================================
!
! module_bqcd.F90
!
!-------------------------------------------------------------------------------
!
! Copyright (C) 2003-2011 Hinnerk Stueben
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
module module_bqcd

  use module_svn
  private :: svn_date, svn_revision
  integer, parameter, private :: n = len(svn_revision)


  character(len = *), parameter :: prog_name = "bqcd"
  character(len = *), parameter :: prog_version = "benchmark2"
  character(len = n), parameter :: prog_revision  = &
                             " (revision" // trim(svn_revision(11:n-1)) // ")"

  integer, parameter :: input_version = 4
  integer, parameter :: conf_info_version = 3

end
!===============================================================================
