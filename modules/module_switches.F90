!===============================================================================
!
! module_switches.F90
!
!-------------------------------------------------------------------------------
!
! Copyright (C) 2003 Hinnerk Stueben
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

module module_switches

  type type_switches
     logical :: quenched
     logical :: dynamical
     logical :: tempering
     logical :: clover
     logical :: h_ext
     logical :: hasenbusch
     logical :: hmc_test
     logical :: measure_polyakov_loop
     logical :: measure_traces
  end type type_switches

  type (type_switches), save :: switches
end
!===============================================================================
