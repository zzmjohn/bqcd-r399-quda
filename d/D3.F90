!===============================================================================
!
! D3.F90 - multiplication with the Wilson hopping matrix D (or D^\dagger)
!          (optimization for Hitachi SR8000: hybrid programming model, 
!           MPI + OpenMP + overlapping communication and computation)
!
!-------------------------------------------------------------------------------
!
! Copyright (C) 1998-2010 Hinnerk Stueben
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

  use module_nn
  use module_vol
  use module_thread
  implicit none
 
  integer :: e, o
  SPINCOL_FIELD :: out, in
  GAUGE_FIELD :: u

  integer :: thread, i1, i2, omp_get_thread_num

  TIMING_START(STRCAT(timing_bin_, NAME))

  !$omp parallel private(thread, i1, i2)

  thread = omp_get_thread_num()

  !$omp barrier

  i1 = xyz_start(thread)
  i2 = xyz_end(thread)

  if (thread == 0) then
     TIMING_START(timing_bin_d_xf)
     call xbound_sc(in, 2)
  else
     call STRCAT(NAME, _xf)(out, in, u(1, 1, 1, e, 1), &
                                     u(1, 1, 1, o, 1), &
                                          nn(1, e, 1, FWD), &
                                          nn(1, e, 1, BWD), i1, i2)
  endif

  !$omp barrier

  if (thread == 0) then
     TIMING_STOP(timing_bin_d_xf)
     TIMING_START(timing_bin_d_yf)
     call xbound_sc(in, 3)
  else
     call STRCAT(NAME, _yf)(out, in, u(1, 1, 1, e, 2), & 
                                     u(1, 1, 1, o, 2), &
                                          nn(1, e, 2, FWD), &
                                          nn(1, e, 2, BWD), i1, i2)
  endif

  !$omp barrier

  if (thread == 0) then
     TIMING_STOP(timing_bin_d_yf)
     TIMING_START(timing_bin_d_zf)
     call xbound_sc(in, 4)
  else
     call STRCAT(NAME, _zf)(out, in, u(1, 1, 1, e, 3), &
                                     u(1, 1, 1, o, 3), &
                                          nn(1, e, 3, FWD), &
                                          nn(1, e, 3, BWD), i1, i2)
  endif

  !$omp barrier

#ifdef TIMING
  if (thread == 0) then
     TIMING_STOP(timing_bin_d_zf)
     TIMING_START(timing_bin_d_t)
  endif
#endif

  i1 = t_start(thread)
  i2 = t_end(thread)

     call STRCAT(NAME, _t )(out, in, u(1, 1, 1, e, 4), &
                                     u(1, 1, 1, o, 4), &
                                          nn(1, e, 4, FWD), &
                                          nn(1, e, 4, BWD), i1, i2)

  !$omp end parallel

  TIMING_STOP(timing_bin_d_t)
  TIMING_STOP(STRCAT(timing_bin_, NAME))

end

!===============================================================================
