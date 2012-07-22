!===============================================================================
!
! hmc.F90 - one Hybrid Monte Carlo step including Hasenbusch's and Bakeyev's 
!           accelerations
!
!// IR scale: tau
!// IR steps: ntau
!// UV scale: tau / m_scale
!// UV steps: m_scale
!//
!// models / splits of actions: 
!//
!// model A:
!//          S_UV = S_g
!//          S_IR = S_det + S_f1
!// 
!// model B:
!//          S_UV = S_g
!//          S_IR = S_det + S_f1 + S_f2
!// 
!// model C:
!//          S_UV = S_g + S_det + S_f1
!//          S_IR = S_f2
!// 
!// S_f1 = phi1+ inv(W+ W) phi1
!// S_f2 = phi2+ W inv(M~+ M~) W+ phi2
!//
!// W = M~ + rho
!//
!// => ir_steps = 1 and rho = 0 corresponds exactly to the previous verions
!//    (rho = 0 is treated as S_f2 = 0), 
!//    and especially tau and ntau have the same meaning as before
!//
!// (In the whole program phi and phi2 are treated asymmetrically.
!// The reason for this is upward compatibility with the mode 
!// "standard Wilson fermions + parallel tempering".  
!// phi is needed for the tempering decisions.)
!//
!// Flags / switches:
!//
!//   force_accept: Force acceptance after, eg, a hot or cold start.
!//   test:         For testing reversibility by forward/backward integration.
!//                 If (test /= HMC_TEST_NONE) force_accept has to be .true.
!//                 If (test == HMC_TEST_BACKWARDS) para%tau has to be reversed.
!// 
!-------------------------------------------------------------------------------
!
! Copyright (C) 1998-2006 Hinnerk Stueben
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
subroutine hmc(para, conf, out, force_accept, test)

  use typedef_hmc
  use module_function_decl
  use module_p_interface
  use module_switches
  use module_vol
  implicit none

  type(hmc_para), intent(in)    :: para
  type(hmc_conf), intent(inout) :: conf
  type(hmc_out),  intent(out)   :: out
  integer,        intent(in)    :: force_accept
  integer,        intent(in)    :: test

  P_GAUGE_FIELD,    save      :: u_bck
  P_SPINCOL_FIELD,  save      :: phi_bck
  P_CLOVER_FIELD_A, save      :: a_bck, i_bck
  P_CLOVER_FIELD_B, save      :: b_bck
  P_GENERATOR_FIELD, save     :: p

  REAL                        :: sd_old, sd_new
  REAL                        :: sf1_old, sf1_new
  REAL                        :: sf2_old, sf2_new
  REAL                        :: sg_old, sg_new
  REAL                        :: sp_old, sp_new
  REAL                        :: hg_old, hg_new
  REAL                        :: h_old, h_new
  REAL                        :: sf1, sf2
  REAL, external              :: sp, sg, clover_action


  TIMING_START(timing_bin_hmc)

  TIMING_START(timing_bin_hmc_init)

  if (.not. associated(u_bck)) then
                          ALLOCATE_G_FIELD(u_bck)
                          ALLOCATE_GEN_FIELD(p)
                          ALLOCATE_SC_FIELD(phi_bck)
     if (switches%clover) call allocate_clover_field_a(a_bck)
     if (switches%clover) call allocate_clover_field_a(i_bck)
     if (switches%clover) call allocate_clover_field_b(b_bck)
  endif

  sd_old = ZERO;  sd_new = ZERO
  sf1_old = ZERO; sf1_new = ZERO
  sf2_old = ZERO; sf2_new = ZERO
  sg_old = ZERO;  sg_new = ZERO
  sp_old = ZERO;  sp_new = ZERO
  sf1 = ZERO;     sf2 = ZERO
 
  call init_cg_stat()
  TIMING_STOP(timing_bin_hmc_init)

  if (test /= HMC_TEST_BACKWARDS) then  ! ie normally do:

     ! backups:

     TIMING_START(timing_bin_hmc_backup)

     call swap_p_sc_field(phi_bck, conf%phi)

                          u_bck = conf%u
     if (switches%clover) a_bck = conf%a
     if (switches%clover) i_bck = conf%i
     if (switches%clover) b_bck = conf%b

     TIMING_STOP(timing_bin_hmc_backup)

     ! initialize momenta p, phi, phi2 and old action:

     call hmc_init_p(p)
     call hmc_init_phi(conf, para, sf1_old, sf2_old)

     TIMING_START(timing_bin_hmc_h_old)
     if (switches%clover) sd_old = clover_action(conf%b(1,1,ODD))
     sg_old = sg(conf%u)
     sp_old = sp(p)
     hg_old = sg_old * para%beta

     h_old = sd_old + sp_old + hg_old + sf1_old + sf2_old
     TIMING_STOP(timing_bin_hmc_h_old)
  endif

  if (test == HMC_TEST_FORWARDS) then
     call hmc_test_report(test, p, conf%u, &
                          sp_old, hg_old, sf1_old, sf2_old, sd_old)
  endif

! leap frog integration:

  call hmc_leap_frog(p, para, conf, sf1, sf2)

! calculate Hamiltonian:

  TIMING_START(timing_bin_hmc_h_new)
  if (switches%clover) sd_new = clover_action(conf%b(1,1,ODD))
  sf1_new = sf1
  sf2_new = sf2
  sg_new = sg(conf%u)
  sp_new = sp(p)
  hg_new = sg_new * para%beta

  h_new = sd_new + sp_new + hg_new + sf1_new + sf2_new
  TIMING_STOP(timing_bin_hmc_h_new)

! accept new U ? :

  TIMING_START(timing_bin_hmc_rest)
  out%exp_dh = exp(h_old - h_new)

  if (force_accept /= 0) then
     out%accepted = 1
  else
     if (ranf() < out%exp_dh) then
        out%accepted = 1
     else
        out%accepted = 0
     endif
  endif

  if (out%accepted == 1) then
     out%sg = sg_new
     out%sf = sf1_new
  else
     call swap_p_sc_field(conf%phi, phi_bck)
     call swap_p_g_field(conf%u, u_bck)
     call swap_p_clover_field_a(conf%a, a_bck)
     call swap_p_clover_field_a(conf%i, i_bck)
     call swap_p_clover_field_b(conf%b, b_bck)
  endif

  call get_cg_stat(out%cg_ncall, out%cg_niter_max, out%cg_niter_tot)

  call iteration_count_write(UREC)
  call hmc_forces_write(UREC)

  if (test == HMC_TEST_BACKWARDS) then
     call hmc_test_report(test, p, conf%u, &
                          sp_new, hg_new,sf1_new, sf2_new, sd_new)
  endif

  TIMING_STOP(timing_bin_hmc_rest)

  TIMING_STOP(timing_bin_hmc)

end

!===============================================================================
