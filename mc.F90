!===============================================================================
!
! mc.F90 - Monte Carlo loop (including tempering)
!
!-------------------------------------------------------------------------------
!
! Copyright (C) 1998-2005 Hinnerk Stueben
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
module typedef_mc_temper

  type temper_out
     REAL    :: delta_h
     integer :: pair
     integer :: swapped
  end type temper_out

end module typedef_mc_temper

!-------------------------------------------------------------------------------
subroutine mc(para, conf)

  use typedef_mc_temper
  use typedef_para
  use module_counter
  use module_function_decl
  use module_switches
  use module_vol
  implicit none

  type(type_para),                         intent(in)    :: para
  type(hmc_conf),   dimension(MAX_TEMPER), intent(inout) :: conf

  type(hmc_out),    dimension(MAX_TEMPER)  :: out
  type(temper_out), dimension(MAX_TEMPER)  :: tmpr
  integer                                  :: i, j
  integer                                  :: iforce, itraj
  integer                                  :: force_accept = 0
  character(len = *), parameter            :: key_fa   = "%fa"
  character(len = *), parameter            :: key_mc   = "%mc"
  character(len = *), parameter            :: key_swap = "%sw"
  REAL                                     :: plaq
  REAL, external                           :: sf, sg

  if (switches%hmc_test) then
     if (para%n_temper /= 1) call die("mc(): HMC-test: n_ensemble has to be 1")
     call hmc_test(para%hmc(1), conf(1))
     return
  endif

  force_accept = 1
  iforce = 0
  do while (counter%traj < 0 .and. counter%j_traj < para%ntraj)

     iforce = iforce + 1
     counter%traj = counter%traj + 1
     counter%j_traj = counter%j_traj + 1

     if (iforce == 1 .and. my_pe() == 0) then
        call begin(UREC, "ForceAcceptance")
        write(UREC, 405) "T", key_fa, "i_fa", "e", &
             "PlaqEnergy", "exp(-Delta_H)", "CGcalls", "CGitTot", "CGitMax"
     endif

     do i = 1, para%n_temper
        call hmc(para%hmc(i), conf(i), out(i), force_accept, HMC_TEST_NONE)

        plaq = out(i)%sg / (SIX * volume)

        if (my_pe() == 0) then
           write(UREC, 415) key_fa, counter%traj, i, plaq, out(i)%exp_dh, &
                out(i)%cg_ncall, out(i)%cg_niter_tot, out(i)%cg_niter_max
        endif
     enddo
  enddo

  if (iforce > 0) then
     call end(UREC, "ForceAcceptance")
  endif

405 format (1x, 2a, 2a5, 2a15,    3a8)
415 format (1x, a4, 2i5, 2f15.10, 3i8)


  force_accept = 0
  itraj = 0
  do while (counter%traj >= 0 .and. counter%traj < para%maxtraj  &
            .and. counter%j_traj < para%ntraj)

     itraj = itraj + 1
     counter%traj = counter%traj + 1
     counter%j_traj = counter%j_traj + 1
     
     if (itraj == 1 .and. my_pe() == 0) then
        call begin(UREC, "MC")

        if (para%n_temper > 1) write(UREC, 450)  &
                "T", key_swap, "traj", "ie", "e1", "e2", "Delta_H", "Acc"

        write(UREC, 400) "T", key_mc, "traj", "e", "f", &
             "PlaqEnergy", "exp(-Delta_H)", "Acc", &
             "CGcalls", "CGitTot", "CGitMax"
     endif

     do i = 1, para%n_temper
        if (itraj == 1) then
           if (switches%tempering .and. switches%dynamical) then
              out(i)%sf = sf(para%hmc(i), conf(i))
           else
              out(i)%sf = ZERO
           endif
           out(i)%sg = sg(conf(i)%u)
        endif
        
        call hmc(para%hmc(i), conf(i), out(i), force_accept, HMC_TEST_NONE)
     enddo
     
     if (para%n_temper > 1 .and. counter%traj > para%nstd) then
        call temper(para%n_temper, para%swap_seq, para%hmc, conf, out, tmpr)
        do i = 1, para%n_temper - 1
           if (my_pe() == 0) write(UREC, 460) key_swap, counter%traj, i, &
                tmpr(i)%pair, tmpr(i)%pair + 1, tmpr(i)%delta_h, tmpr(i)%swapped
        enddo
     endif
     
     do i = 1, para%n_temper
        j = conf(i)%former

        if (my_pe() == 0) write(UREC,410) key_mc, counter%traj, i, j, &
                out(i)%sg / (SIX * volume), out(i)%exp_dh, out(i)%accepted, &
                out(i)%cg_ncall, out(i)%cg_niter_tot, out(i)%cg_niter_max

        if (switches%measure_traces) &
            call traces(para%hmc(i), conf(i), counter%traj, i, j)

        if (switches%measure_polyakov_loop) &
            call polyakov_loop(conf(i), counter%traj, i, j)
    
        call cooling(conf(i)%u, counter%traj, i, j)
     enddo

     if (para%nsave > 0) then
        if (mod(counter%traj, para%nsave) == 0) then
           call conf_write(.false., para, conf)
        endif
     endif

  enddo

  if (itraj > 0) then
     call end(UREC, "MC")
  endif

400 format (1x, 2a, a6, 2a3, 2a15,    a4, 3a8)
410 format (1x, a4, i6, 2i3, 2f15.10, i4, 3i8)

450 format (1x, 2a, a6, 3a3, a18,    a4)
460 format (1x, a4, i6, 3i3, f18.10, i4)

end

!-------------------------------------------------------------------------------
subroutine temper(n_temper, swap_seq, para, conf, action, tmpr)

  use typedef_hmc
  use typedef_mc_temper
  use module_function_decl
  use module_p_interface
  implicit none
   
  integer,                               intent(in)    :: n_temper, swap_seq
  type(hmc_para),   dimension(n_temper), intent(in)    :: para
  type(hmc_conf),   dimension(n_temper), intent(inout) :: conf
  type(hmc_out),    dimension(n_temper), intent(inout) :: action
  type(temper_out), dimension(n_temper), intent(out)   :: tmpr

  integer,          dimension(n_temper)                :: pair
  integer                                              :: i_pair, n_pair, i, j
  integer                                              :: swapped
  REAL,             dimension(2, 2)                    :: hf, hg 
  REAL,             external                           :: sf, sg
  REAL                                                 :: h_old, h_new, delta_h
  REAL                                                 :: random

  if (n_temper == 1) return

  do i = 1, n_temper
     conf(i)%former = i
  enddo

  n_pair = n_temper - 1

  call swap_sequence(swap_seq, pair, n_pair)

  do i_pair = 1, n_pair
     i = pair(i_pair)
     j = pair(i_pair) + 1

     hg(1, 1) = para(i)%beta * action(i)%sg
     hg(2, 2) = para(j)%beta * action(j)%sg
     
     hg(1, 2) = para(i)%beta * action(j)%sg
     hg(2, 1) = para(j)%beta * action(i)%sg
     
     hf(1, 1) = action(i)%sf
     hf(2, 2) = action(j)%sf
     hf(1, 2) = sf(para(i), conf(j))
     hf(2, 1) = sf(para(j), conf(i))
     
     if (para(i)%kappa == ZERO) then
        if (hf(1, 1) /= ZERO) call die("temper(): hf(1, 1) /= 0 ")
        if (hf(2, 1) /= ZERO) call die("temper(): hf(2, 1) /= 0 ")
        if (hf(1, 2) /= ZERO) call die("temper(): hf(1, 2) /= 0 ")
        if (hf(2, 2) /= ZERO) call die("temper(): hf(2, 2) /= 0 ")
     endif

     h_old = hg(1, 1) + hf(1, 1) + hg(2, 2) + hf(2, 2)
     h_new = hg(1, 2) + hf(1, 2) + hg(2, 1) + hf(2, 1)
     
     delta_h = h_new - h_old

     if (ranf() < exp(-delta_h)) then
        swapped = 1
     else
        swapped = 0
     endif

     if (swapped /= 0) then
        call swap_p_g_field(conf(i)%u, conf(j)%u)
        call swap_p_sc_field(conf(i)%phi, conf(j)%phi)
        call swap_integer(conf(i)%former, conf(j)%former)
        call swap_real(action(i)%sg, action(j)%sg)
        action(i)%sf = hf(1, 2)
        action(j)%sf = hf(2, 1)
     endif

     tmpr(i_pair)%pair    = i
     tmpr(i_pair)%swapped = swapped
     tmpr(i_pair)%delta_h = delta_h
  enddo

  call check_former(n_temper, conf)

end

!-------------------------------------------------------------------------------
subroutine check_former(n_temper, conf)

  use typedef_hmc
  use module_function_decl
  implicit none

  integer,        intent(in)    :: n_temper
  type(hmc_conf), intent(inout) :: conf(n_temper)
  integer                       :: count(n_temper), i

  count = 0
  do i = 1, n_temper
     count(conf(i)%former) = count(conf(i)%former) + 1
  enddo

  do i = 1, n_temper
     if (count(i) /= 1) call die("check_former(): error")
  enddo
end

!-------------------------------------------------------------------------------
subroutine swap_sequence(type, s, n)

  implicit none
  integer, intent(in) :: type, n
  integer, intent(out) :: s(n)
  integer :: i

  select case (type)
     case (SWAP_UP)
        do i = 1, n
           s(i) = i
        enddo
     case (SWAP_DOWN)
        do i = 1, n
           s(i) = n - i + 1
        enddo
     case (SWAP_RANDOM)
        call random_sequence(s, n)
     case default
        call die("swap_sequence(): don't know how to build the sequence")
  end select

end

!===============================================================================
