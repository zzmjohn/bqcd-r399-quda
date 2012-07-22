!===============================================================================
!
! hmc_test.F90 - forward/backward leap frog integration
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
subroutine hmc_test(para, conf)

  use typedef_hmc
  use module_function_decl
  use module_vol
  implicit none

  type(hmc_para), intent(inout) :: para
  type(hmc_conf), intent(inout) :: conf
  type(hmc_out)                 :: out

  call begin(UREC, "HMCtest")

  call hmc(para, conf, out, .true., HMC_TEST_FORWARDS)

  call write_out("forward  ")

  para%tau = -para%tau

  call hmc(para, conf, out, .true., HMC_TEST_BACKWARDS)

  para%tau = -para%tau

  call write_out("backward ")

  call end(UREC, "HMCtest")


contains

  subroutine write_out(direction)

     character(*) :: direction
     REAL :: plaq

     if (my_pe() == 0) then
        plaq = out%sg / (6 * volume)
        write(UREC, *)
        write(UREC, 400) "Direction", "PlaqEnergy", &
                         "exp(-Delta_H)", "CGcalls", "CGitTot", "CGitMax"
        write(UREC, 410) direction, plaq, out%exp_dh, &
                         out%cg_ncall, out%cg_niter_tot, out%cg_niter_max
        write(UREC, *)
     endif

400 format (1x, a, 2a15,    3a8)
410 format (1x, a, 2f15.10, 3i8)

  end subroutine write_out

end

!-------------------------------------------------------------------------------
subroutine hmc_test_report(test, p, u, hp, hg, hf1, hf2, hd)

  use module_function_decl
  use module_vol
  implicit none

  integer,         intent(in) :: test
  GENERATOR_FIELD, intent(in) :: p
  GAUGE_FIELD,     intent(in) :: u
  REAL,            intent(in) :: hp, hg, hf1, hf2, hd

  P_GENERATOR_FIELD, save     :: p_start
  P_GAUGE_FIELD,     save     :: u_start
  REAL,              save     :: hp_start
  REAL,              save     :: hg_start
  REAL,              save     :: hf1_start
  REAL,              save     :: hf2_start
  REAL,              save     :: hd_start

  REAL                        :: diff_p, diff_u
  integer                     :: i, eo, mu, j, c1, c2

  if (.not. associated(p_start)) then
     allocate(p_start(NGEN, volh_tot, EVEN:ODD, DIM))
     allocate(u_start(NCOL, NCOL, volh_tot, EVEN:ODD, DIM))
  endif

  if (test == HMC_TEST_FORWARDS) then

     p_start = p
     u_start = u
     hp_start = hp
     hg_start = hg
     hf1_start = hf1
     hf2_start = hf2
     hd_start = hd

  else if (test == HMC_TEST_BACKWARDS) then

     diff_p = ZERO
     diff_u = ZERO

     do mu = 1, DIM
        do eo = EVEN, ODD
           do i = 1, volh
              do j = 1, NGEN
                 diff_p = max(diff_p, abs(p_start(j,i,eo,mu) - p(j,i,eo,mu)))
              enddo
              do c2 = 1, NCOL
                 do c1 = 1, NCOL
                    diff_u = max(diff_u, &
                              abs(relative_change(Re(u_start(c1,c2,i,eo,mu)), &
                                                        Re(u(c1,c2,i,eo,mu)))))
                    diff_u = max(diff_u, &
                              abs(relative_change(Im(u_start(c1,c2,i,eo,mu)), &
                                                        Im(u(c1,c2,i,eo,mu)))))
                 enddo
              enddo
           enddo
        enddo
     enddo

    if (my_pe() == 0) then   
     write(UREC,  *) 
     write(UREC,400) "Configuration changes (maximal abs. relative changes):"
     write(UREC,  *) 
     write(UREC,410) "Generator field:", diff_p
     write(UREC,410) "Gauge field:    ", diff_u
     write(UREC,  *) 
     write(UREC,  *) 
     write(UREC,400) "Energy changes:"
     write(UREC,  *) 
     write(UREC,420) "Energy     ", "old value", "rel.change"
     write(UREC,  *)
     write(UREC,430) "H_generator", hp_start,  relative_change(hp_start,  hp)
     write(UREC,430) "H_gauge    ", hg_start,  relative_change(hg_start,  hg)
     write(UREC,430) "H_fermion_1", hf1_start, relative_change(hf1_start, hf1)
     write(UREC,430) "H_fermion_2", hf2_start, relative_change(hf2_start, hf2)
     write(UREC,430) "H_det      ", hd_start,  relative_change(hd_start,  hd)
     write(UREC,  *)
    endif

400 format (1x, a)
410 format (1x, a, e8.1)
420 format (1x, a, a20,    a12)
430 format (1x, a, e20.10, e12.1)

  else
     call die("hmc_test_report(): illegal test flag.")   
  endif

contains

  REAL function relative_change(old, new)

     implicit none
     REAL, intent(in) :: old, new

     if (old == ZERO) then
        relative_change = ZERO
     else
        relative_change = (new - old) / old
     endif

  end function relative_change

end
!===============================================================================
