!===============================================================================
!
! module_svn.F90 - svn revision and date information
!
! Revision and Date kept here are *global* if this file is modified before
! "svn commit".
!
! The modification is to increment "count" (see target "commit" in ../Makefile):
!
!count= 16
!
!-------------------------------------------------------------------------------
!
! Copyright (C) 2007-2008 Hinnerk Stueben 
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

module module_svn

  character(*), parameter :: svn_revision = "$Revision: 399 $"
  character(*), parameter :: svn_date = "$Date: 2012-03-02 16:54:59 +0100 (Fri, 02 Mar 2012) $"

end

!===============================================================================

