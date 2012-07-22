!===============================================================================
!
! D21.F90 - multiplication with the Wilson hopping matrix D (or D^\dagger)
!           projection onto 2 spincol components
!
!-------------------------------------------------------------------------------
!
! Copyright (C) 2006 Hinnerk Stueben
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
subroutine NAME(e, o, out, in, u)
 
! out := NAME in 
!
! NAME = d or d_dag
!
! out is of type "e" = EVEN or ODD
! in is of type "o" = ODD or EVEN

  use module_d21
  use module_nn
  use module_vol
  use module_p_interface
  implicit none
 
  integer :: e, o
  SPINCOL_FIELD :: out, in
  GAUGE_FIELD :: u


  TIMING_START(STRCAT(timing_bin_, NAME))

  ALLOCATE_SC2_FIELD(a)

  call STRCAT(NAME, _projection)(a, in)

!!call xbound_sc2_field(a)
  call xbound_sc2_field_i(a)

  call STRCAT(NAME, _t )(out, a(1, 1, 1, 4, FWD), a(1, 1, 1, 4, BWD),  & 
                              u(1, 1, 1, e, 4), u(1, 1, 1, o, 4),      &
                             nn(1, e, 4, FWD), nn(1, e, 4, BWD), VOLH)

  call STRCAT(NAME, _zf)(out, a(1, 1, 1, 3, FWD), a(1, 1, 1, 3, BWD),  & 
                              u(1, 1, 1, e, 3), u(1, 1, 1, o, 3),      &
                             nn(1, e, 3, FWD), nn(1, e, 3, BWD), VOLH)

  call STRCAT(NAME, _yf)(out, a(1, 1, 1, 2, FWD), a(1, 1, 1, 2, BWD),  & 
                              u(1, 1, 1, e, 2), u(1, 1, 1, o, 2),      &
                             nn(1, e, 2, FWD), nn(1, e, 2, BWD), VOLH)

  call STRCAT(NAME, _xf)(out, a(1, 1, 1, 1, FWD), a(1, 1, 1, 1, BWD),  & 
                              u(1, 1, 1, e, 1), u(1, 1, 1, o, 1),      &
                             nn(1, e, 1, FWD), nn(1, e, 1, BWD), VOLH)

  TIMING_STOP(STRCAT(timing_bin_, NAME))

end

!-------------------------------------------------------------------------------
subroutine STRCAT(NAME, _projection)(out, in)

  use module_vol
  implicit none

  SC2_FIELD, intent(out)    :: out
  SPINCOL_FIELD, intent(in) :: in
  integer                   :: i

  ! statement function:
  
  COMPLEX :: i_times, c
  i_times(c) = cmplx(-aimag(c), real(c), kind = RKIND)

  TIMING_START(timing_bin_sc2_projection)

#ifdef DAGGER
# define PLUS -
# define MINUS +
# define D_T_ONE 3
# define D_T_TWO 4
# define D_T_THREE 1 
# define D_T_FOUR 2
#else
# define PLUS +
# define MINUS -
# define D_T_ONE 1
# define D_T_TWO 2
# define D_T_THREE 3 
# define D_T_FOUR 4
#endif

  
  !$omp parallel do
  do i = 1, volh

     out(1, 1, i, 1, FWD) = in(1, 1, i) MINUS i_times(in(4, 1, i))
     out(2, 1, i, 1, FWD) = in(2, 1, i) MINUS i_times(in(3, 1, i))
     out(1, 1, i, 1, BWD) = in(1, 1, i) PLUS  i_times(in(4, 1, i))
     out(2, 1, i, 1, BWD) = in(2, 1, i) PLUS  i_times(in(3, 1, i))

     out(1, 2, i, 1, FWD) = in(1, 2, i) MINUS i_times(in(4, 2, i))
     out(2, 2, i, 1, FWD) = in(2, 2, i) MINUS i_times(in(3, 2, i))
     out(1, 2, i, 1, BWD) = in(1, 2, i) PLUS  i_times(in(4, 2, i))
     out(2, 2, i, 1, BWD) = in(2, 2, i) PLUS  i_times(in(3, 2, i))

     out(1, 3, i, 1, FWD) = in(1, 3, i) MINUS i_times(in(4, 3, i))
     out(2, 3, i, 1, FWD) = in(2, 3, i) MINUS i_times(in(3, 3, i))
     out(1, 3, i, 1, BWD) = in(1, 3, i) PLUS  i_times(in(4, 3, i))
     out(2, 3, i, 1, BWD) = in(2, 3, i) PLUS  i_times(in(3, 3, i))


     out(1, 1, i, 2, FWD) = in(1, 1, i) MINUS in(4, 1, i)
     out(2, 1, i, 2, FWD) = in(2, 1, i) PLUS  in(3, 1, i)
     out(1, 1, i, 2, BWD) = in(1, 1, i) PLUS  in(4, 1, i)
     out(2, 1, i, 2, BWD) = in(2, 1, i) MINUS in(3, 1, i)

     out(1, 2, i, 2, FWD) = in(1, 2, i) MINUS in(4, 2, i)
     out(2, 2, i, 2, FWD) = in(2, 2, i) PLUS  in(3, 2, i)
     out(1, 2, i, 2, BWD) = in(1, 2, i) PLUS  in(4, 2, i)
     out(2, 2, i, 2, BWD) = in(2, 2, i) MINUS in(3, 2, i)

     out(1, 3, i, 2, FWD) = in(1, 3, i) MINUS in(4, 3, i)
     out(2, 3, i, 2, FWD) = in(2, 3, i) PLUS  in(3, 3, i)
     out(1, 3, i, 2, BWD) = in(1, 3, i) PLUS  in(4, 3, i)
     out(2, 3, i, 2, BWD) = in(2, 3, i) MINUS in(3, 3, i)


     out(1, 1, i, 3, FWD) = in(1, 1, i) MINUS i_times(in(3, 1, i))
     out(2, 1, i, 3, FWD) = in(2, 1, i) PLUS  i_times(in(4, 1, i))
     out(1, 1, i, 3, BWD) = in(1, 1, i) PLUS  i_times(in(3, 1, i))
     out(2, 1, i, 3, BWD) = in(2, 1, i) MINUS i_times(in(4, 1, i))

     out(1, 2, i, 3, FWD) = in(1, 2, i) MINUS i_times(in(3, 2, i))
     out(2, 2, i, 3, FWD) = in(2, 2, i) PLUS  i_times(in(4, 2, i))
     out(1, 2, i, 3, BWD) = in(1, 2, i) PLUS  i_times(in(3, 2, i))
     out(2, 2, i, 3, BWD) = in(2, 2, i) MINUS i_times(in(4, 2, i))

     out(1, 3, i, 3, FWD) = in(1, 3, i) MINUS i_times(in(3, 3, i))
     out(2, 3, i, 3, FWD) = in(2, 3, i) PLUS  i_times(in(4, 3, i))
     out(1, 3, i, 3, BWD) = in(1, 3, i) PLUS  i_times(in(3, 3, i))
     out(2, 3, i, 3, BWD) = in(2, 3, i) MINUS i_times(in(4, 3, i))


     out(1, 1, i, 4, FWD) = in(D_T_THREE, 1, i)
     out(2, 1, i, 4, FWD) = in(D_T_FOUR,  1, i)
     out(1, 1, i, 4, BWD) = in(D_T_ONE,   1, i)
     out(2, 1, i, 4, BWD) = in(D_T_TWO,   1, i)
 
     out(1, 2, i, 4, FWD) = in(D_T_THREE, 2, i)
     out(2, 2, i, 4, FWD) = in(D_T_FOUR,  2, i)
     out(1, 2, i, 4, BWD) = in(D_T_ONE,   2, i)
     out(2, 2, i, 4, BWD) = in(D_T_TWO,   2, i)
 
     out(1, 3, i, 4, FWD) = in(D_T_THREE, 3, i)
     out(2, 3, i, 4, FWD) = in(D_T_FOUR,  3, i)
     out(1, 3, i, 4, BWD) = in(D_T_ONE,   3, i)
     out(2, 3, i, 4, BWD) = in(D_T_TWO,   3, i)
 
 enddo

 TIMING_STOP(timing_bin_sc2_projection)

end

!===============================================================================
