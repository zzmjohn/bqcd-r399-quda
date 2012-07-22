!===============================================================================
!
! typedef_clover.F90
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
# include "defs.h"

!-------------------------------------------------------------------------------
module typedef_clover

  type type_clover_a
    sequence

    REAL i11
    REAL i22

    COMPLEX i12
    COMPLEX i13
    COMPLEX i14
    COMPLEX i15
    COMPLEX i16

    COMPLEX i23
    COMPLEX i24
    COMPLEX i25
    COMPLEX i26

    REAL i33
    REAL i44

    COMPLEX i34
    COMPLEX i35
    COMPLEX i36

    COMPLEX i45
    COMPLEX i46

    REAL i55
    REAL i66

    COMPLEX i56
  end type type_clover_a


  type type_clover_b
    sequence

    COMPLEX i21

    COMPLEX i31
    COMPLEX i32

    COMPLEX i41
    COMPLEX i42
    COMPLEX i43

    COMPLEX i51
    COMPLEX i52
    COMPLEX i53
    COMPLEX i54

    COMPLEX i61
    COMPLEX i62
    COMPLEX i63
    COMPLEX i64
    COMPLEX i65

    REAL i11
    REAL i22
    REAL i33
    REAL i44
    REAL i55
    REAL i66
  end type type_clover_b

end
!===============================================================================
