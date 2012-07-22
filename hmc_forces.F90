!===============================================================================
!
! hmc_forces.F90  -  calculation of HMC forces in Hasenbusch improvement,
!                    does not work with tempering
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

!-----------------------------------------------------------------------------
subroutine hmc_forces_old(p)

  use module_hmc_forces
  use module_p_interface
  use module_switches
  use module_vol
  implicit none

  GENERATOR_FIELD, intent(in) :: p
  integer                     :: mu, eo, i

  if (.not. switches%hasenbusch) return

  if (.not. associated(p_old)) then
     call allocate_gen_field(p_old)
     f_count = ZERO
     f_avg = ZERO
     f_max = ZERO
  endif


  do mu = 1, DIM
     do eo = EVEN, ODD
        !$omp parallel do
        do i = 1, volh
           p_old(1, i, eo, mu) = p(1, i, eo, mu)
           p_old(2, i, eo, mu) = p(2, i, eo, mu)
           p_old(3, i, eo, mu) = p(3, i, eo, mu)
           p_old(4, i, eo, mu) = p(4, i, eo, mu)
           p_old(5, i, eo, mu) = p(5, i, eo, mu)
           p_old(6, i, eo, mu) = p(6, i, eo, mu)
           p_old(7, i, eo, mu) = p(7, i, eo, mu)
           p_old(8, i, eo, mu) = p(8, i, eo, mu)
        enddo
     enddo
  enddo

end

!-----------------------------------------------------------------------------
subroutine hmc_forces_new(p, step, which)

  use module_hmc_forces
  use module_function_decl
  use module_switches
  use module_vol
  implicit none

  GENERATOR_FIELD, intent(in) :: p
  REAL, intent(in)            :: step
  integer, intent(in)         :: which

  integer                     :: mu, eo, i
  REAL                        :: force
  
  if (.not. switches%hasenbusch) return

  force = ZERO
  do mu = 1, DIM
     do eo = EVEN, ODD
        !$omp parallel do
        do i = 1, volh
           force = force &
              + (p_old(1, i, eo, mu) - p(1, i, eo, mu))**2 &
              + (p_old(2, i, eo, mu) - p(2, i, eo, mu))**2 &
              + (p_old(3, i, eo, mu) - p(3, i, eo, mu))**2 &
              + (p_old(4, i, eo, mu) - p(4, i, eo, mu))**2 &
              + (p_old(5, i, eo, mu) - p(5, i, eo, mu))**2 &
              + (p_old(6, i, eo, mu) - p(6, i, eo, mu))**2 &
              + (p_old(7, i, eo, mu) - p(7, i, eo, mu))**2 &
              + (p_old(8, i, eo, mu) - p(8, i, eo, mu))**2
        enddo
     enddo
  enddo

  force = global_sum(force) / (NGEN * volume * DIM)
  force = sqrt(force) / abs(step)

  f_count(which) = f_count(which) + ONE
  f_avg(which) = f_avg(which) + force
  f_max(which) = max(f_max(which), force)

end


!-----------------------------------------------------------------------------
subroutine hmc_forces_write(unit)

  use module_hmc_forces
  use module_counter
  use module_function_decl
  use module_switches

  implicit none

  integer, intent(in)     :: unit
  integer, save           :: written = 0 
  integer                 :: i

  character(*), parameter :: key_avg   = "%Favg"
  character(*), parameter :: key_max   = "%Fmax"
  character(*), parameter :: fmt_h = "(1x, 2a, a6, 4a)"
  character(*), parameter :: fmt_b = "(1x, a6, i6, 4g20.10)"

  character(20), dimension(n_force) :: f_name

  if (.not. switches%hasenbusch) return

  f_name(i_sg)  = "         F_gauge"
  f_name(i_sd)  = "           F_det"
  f_name(i_sf1) = "            F_F1"
  f_name(i_sf2) = "            F_F2"

  do i = 1, n_force
     if (f_count(i) /= ZERO) then
        f_avg(i) = f_avg(i) / f_count(i)
     endif
  enddo

  if (written == 0 .and. my_pe() == 0) then
     write(unit, fmt_h) "T", key_avg, "traj", (f_name(i), i = 1, n_force)
     write(unit, fmt_h) "T", key_max, "traj", (f_name(i), i = 1, n_force)
  endif
     
  if (my_pe() == 0) then
     write(unit, fmt_b) key_avg, counter%traj, (f_avg(i), i = 1, n_force)
     write(unit, fmt_b) key_max, counter%traj, (f_max(i), i = 1, n_force)
  endif

    
  written = written + 1
  f_count = ZERO
  f_avg = ZERO
  f_max = ZERO

end

!===============================================================================
