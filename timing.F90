!===============================================================================
!
! timing.F90 - measurements of execution times and performance
!
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
module module_timing_meas

  implicit none

  type type_timing
     SECONDS time
     SECONDS total_time
     integer n_call
     integer fill_the_cache_line_1
     integer fill_the_cache_line_2
     integer fill_the_cache_line_3
  end type type_timing

  integer, parameter        :: n_timing = 55
  type(type_timing), save   :: meas(n_timing)  ! measurements
!dir$ cache_align              meas

  data meas /n_timing * type_timing(0.0, 0.0, 0, 0, 0, 0)/

end

!-------------------------------------------------------------------------------
module module_timing_write

  use      module_timing_meas
  implicit none


  character(len = 16), save :: text(n_timing)  ! region name
  integer(8), save          :: n_op(n_timing)  ! # operations

  data text(timing_bin_d_xf)       /"d_xf"/
  data text(timing_bin_d_xb)       /"d_xb"/
  data text(timing_bin_d_yf)       /"d_yf"/
  data text(timing_bin_d_yb)       /"d_yb"/
  data text(timing_bin_d_zf)       /"d_zf"/
  data text(timing_bin_d_zb)       /"d_zb"/
  data text(timing_bin_d_t)        /"d_t"/
  data text(timing_bin_d)          /"D_TOTAL"/

  data text(timing_bin_global_sum)     /"global_sum"/
  data text(timing_bin_global_sum_vec) /"global_sum_vec"/
  data text(timing_bin_sc_zero)        /"sc_zero"/
  data text(timing_bin_sc_copy)        /"sc_copy"/
  data text(timing_bin_sc_scale)       /"sc_scale"/
  data text(timing_bin_sc_axpy)        /"sc_axpy"/
  data text(timing_bin_sc_caxpy)       /"sc_caxpy"/
  data text(timing_bin_sc_caxpy2)      /"sc_caxpy2"/
  data text(timing_bin_sc_cax2)        /"sc_cax2"/
  data text(timing_bin_sc_xpby)        /"sc_xpby"/
  data text(timing_bin_sc_axpby)       /"sc_axpby"/
  data text(timing_bin_sc_norm2)       /"sc_norm2"/
  data text(timing_bin_sc_dot)         /"sc_dot"/
  data text(timing_bin_sc_cdotc)       /"sc_cdotc"/

  data text(timing_bin_plaq)       /"plaquette"/
  data text(timing_bin_cooling)    /"cooling"/
  data text(timing_bin_u_read)     /"u_read"/
  data text(timing_bin_u_write)    /"u_write"/

  data text(timing_bin_total)      /"TOTAL"/
  data text(timing_bin_hmc)        /"HMC"/
  data text(timing_bin_cg)         /"CG"/
  data text(timing_bin_mtdagmt)    /"MTDAGMT"/

  data text(timing_bin_dsf)        /"dsf"/
  data text(timing_bin_dsg)        /"dsg"/
  data text(timing_bin_hmc_u)      /"hmc_u"/
  data text(timing_bin_hmc_init_p) /"hmc_init_p"/

  data text(timing_bin_clover_init)   /"clover_init"/
  data text(timing_bin_clover_mult_a) /"clover_mult_a"/
  data text(timing_bin_clover_mult_ao)/"clover_mult_ao"/
  data text(timing_bin_clover_mult_b) /"clover_mult_b"/
  data text(timing_bin_clover_dsd)    /"clover_dsd"/
  data text(timing_bin_clover_dsf)    /"clover_dsf"/

  data text(timing_bin_hmc_init)      /"hmc_init"/
  data text(timing_bin_hmc_momenta)   /"hmc_momenta"/
  data text(timing_bin_hmc_init_phi)  /"hmc_phi"/
  data text(timing_bin_hmc_h_old)     /"hmc_h_old"/
  data text(timing_bin_hmc_backup)    /"hmc_backup"/
  data text(timing_bin_hmc_half_step0)/"hmc_half_step0"/
  data text(timing_bin_hmc_half_step1)/"hmc_half_step1"/
  data text(timing_bin_hmc_xbound_g)  /"hmc_xbound_g"/
  data text(timing_bin_hmc_steps)     /"hmc_steps"/
  data text(timing_bin_hmc_h_new)     /"hmc_h_new"/
  data text(timing_bin_hmc_rest)      /"hmc_rest"/

  data text(timing_bin_h_mult_a)      /"h_mult_a"/
  data text(timing_bin_h_mult_b)      /"h_mult_b"/
  data text(timing_bin_h_mult_c)      /"h_mult_c"/

  data text(timing_bin_sc2_projection)/"sc2_projection"/

  integer(8), parameter          :: op_add   = 2  ! operations per complex add.
  integer(8), parameter          :: op_mult  = 6  ! operations per complex mult.

  integer(8), parameter, private :: op_d_xyz = 18 * op_mult + 30 * op_add
  integer(8), parameter, private :: op_d_t   = 36 * op_mult + 24 * op_add + 24

  integer(8), parameter          :: op_d     = 6 * op_d_xyz + op_d_t

  integer(8), parameter, private :: op_uuu   = 162 * op_mult
  integer(8), parameter, private :: op_re_tr = 36
  integer(8), parameter, private :: op_plaq  = 2 * 6 * (op_uuu + op_re_tr)

  integer(8), parameter          :: op_sc_r    = NDIRAC * NCOL * SIZE_COMPLEX
  integer(8), parameter          :: op_sc_c    = NDIRAC * NCOL

  integer(8), parameter          :: op_blas_r1 = op_sc_r
  integer(8), parameter          :: op_blas_r2 = op_sc_r * 2
  integer(8), parameter          :: op_blas_r3 = op_sc_r * 3
  integer(8), parameter          :: op_blas_c1 = op_sc_c * op_mult
  integer(8), parameter          :: op_blas_c2 = op_sc_c * (op_mult + op_add)
  integer(8), parameter          :: op_blas_c3 = 2 * op_blas_c2

  integer(8), parameter          :: op_clov  = 60 * op_mult + 84 * op_add + 24
  integer(8), parameter          :: op_h_mult= NDIRAC * NCOL * 3

  integer(8), parameter          :: op_sc2_proj = 36 * op_add

  data n_op(timing_bin_d_xf)    /op_d_xyz/
  data n_op(timing_bin_d_xb)    /op_d_xyz/
  data n_op(timing_bin_d_yf)    /op_d_xyz/
  data n_op(timing_bin_d_yb)    /op_d_xyz/
  data n_op(timing_bin_d_zf)    /op_d_xyz/
  data n_op(timing_bin_d_zb)    /op_d_xyz/
  data n_op(timing_bin_d_t)     /op_d_t/
  data n_op(timing_bin_d)       /op_d/

  data n_op(timing_bin_global_sum)     /0/
  data n_op(timing_bin_global_sum_vec) /0/
  data n_op(timing_bin_sc_zero)        /0/
  data n_op(timing_bin_sc_copy)        /0/
  data n_op(timing_bin_sc_scale)       /op_blas_r1/
  data n_op(timing_bin_sc_norm2)       /op_blas_r2/
  data n_op(timing_bin_sc_dot)         /op_blas_r2/
  data n_op(timing_bin_sc_axpy)        /op_blas_r2/
  data n_op(timing_bin_sc_xpby)        /op_blas_r2/
  data n_op(timing_bin_sc_axpby)       /op_blas_r3/
  data n_op(timing_bin_sc_cdotc)       /op_blas_c2/
  data n_op(timing_bin_sc_caxpy)       /op_blas_c2/
  data n_op(timing_bin_sc_caxpy2)      /op_blas_c3/
  data n_op(timing_bin_sc_cax2)        /op_blas_c3/

  data n_op(timing_bin_plaq)       /op_plaq/
  data n_op(timing_bin_cooling)    /0/
  data n_op(timing_bin_u_read)     /0/
  data n_op(timing_bin_u_write)    /0/

  data n_op(timing_bin_total)      /0/
  data n_op(timing_bin_hmc)        /0/
  data n_op(timing_bin_cg)         /0/
  data n_op(timing_bin_mtdagmt)    /0/

  data n_op(timing_bin_dsf)        /0/
  data n_op(timing_bin_dsg)        /0/
  data n_op(timing_bin_hmc_u)      /0/

  data n_op(timing_bin_clover_init)   /0/
  data n_op(timing_bin_clover_mult_a) /op_clov/
  data n_op(timing_bin_clover_mult_ao)/op_clov/
  data n_op(timing_bin_clover_mult_b) /op_clov/
  data n_op(timing_bin_clover_dsd)    /0/
  data n_op(timing_bin_clover_dsf)    /0/

  data n_op(timing_bin_hmc_init)      /0/
  data n_op(timing_bin_hmc_momenta)   /0/
  data n_op(timing_bin_hmc_init_phi)  /0/
  data n_op(timing_bin_hmc_h_old)     /0/
  data n_op(timing_bin_hmc_backup)    /0/
  data n_op(timing_bin_hmc_half_step0)/0/
  data n_op(timing_bin_hmc_half_step1)/0/
  data n_op(timing_bin_hmc_xbound_g)  /0/
  data n_op(timing_bin_hmc_steps)     /0/
  data n_op(timing_bin_hmc_h_new)     /0/
  data n_op(timing_bin_hmc_rest)      /0/

  data n_op(timing_bin_h_mult_a)      /op_h_mult/
  data n_op(timing_bin_h_mult_b)      /op_h_mult/
  data n_op(timing_bin_h_mult_c)      /op_h_mult/
  data n_op(timing_bin_sc2_projection)/op_sc2_proj/

end

!-------------------------------------------------------------------------------
subroutine timing_start(bin)

  use      module_timing_meas
  implicit none
  integer  bin
  SECONDS  sekunden

  meas(bin)%time = sekunden()
end

!-------------------------------------------------------------------------------
subroutine timing_stop(bin)

  use      module_timing_meas
  implicit none
  integer  bin
  SECONDS  sekunden

  meas(bin)%total_time = meas(bin)%total_time + sekunden() - meas(bin)%time
  meas(bin)%n_call = meas(bin)%n_call + 1
end

!-------------------------------------------------------------------------------
subroutine timing_write(unit)

  use      module_timing_write
  use      module_cg
  use      module_function_decl
  use      module_switches
  use      module_thread
  use      module_vol
  implicit none

  integer  unit, i
  integer  ierror
  integer(8)  op_mtdagmt
  integer(8)  cg_calls, cg_iter
  REAL     mflops, mflops_mean, mflops_min, mflops_max
  REAL     total_gflops, time_mean


  character(len = 8) :: a_mflops_mean, a_mflops_min, a_mflops_max, &
                        a_total_gflops, a_time_mean, a_n_call
  
  character(*), parameter :: ifmt = "(i8)", ffmt = "(f8.2)", &
                             tab_fmt ="(2(3(1x,a),2x),1x,a)"

  if (version_of_d() >= 2) then
     n_op(timing_bin_d_xf) = n_op(timing_bin_d_xf) * 2
     n_op(timing_bin_d_yf) = n_op(timing_bin_d_yf) * 2
     n_op(timing_bin_d_zf) = n_op(timing_bin_d_zf) * 2
  endif

                       op_mtdagmt = 2 * (2 * op_d + op_blas_r2)
  if (switches%clover) op_mtdagmt = op_mtdagmt + 4 * op_clov
  if (switches%h_ext)  op_mtdagmt = op_mtdagmt + 4 * op_h_mult

  if (version_of_d() == 21 .or. version_of_d() == 22) then
     n_op(timing_bin_d_xf) = n_op(timing_bin_d_xf) - 12 * op_add
     n_op(timing_bin_d_yf) = n_op(timing_bin_d_yf) - 12 * op_add
     n_op(timing_bin_d_zf) = n_op(timing_bin_d_zf) - 12 * op_add
  endif

  cg_calls = meas(timing_bin_cg)%n_call
  cg_iter  = cg_iterations_total

  n_op(timing_bin_mtdagmt) = op_mtdagmt
  n_op(timing_bin_cg)      = cg_iter * (op_mtdagmt + 5 * op_blas_r2) &
                           + cg_calls * (op_mtdagmt + op_blas_r1)
  n_op(timing_bin_cg)      = nint( &
                             real(n_op(timing_bin_cg), kind = RKIND) / &
                             real(cg_calls), kind = RKIND)


  call begin(unit, "Timing")
  if (my_pe() == 0) then
     write(unit,"(48x,a)") "Performance"
     write(unit, tab_fmt) "region          ", "  #calls", "    time", &
          "    mean", "     min", "     max", "   Total"
     write(unit, tab_fmt) "                ", "        ", "       s", &
          " Mflop/s", " Mflop/s", " Mflop/s", " Gflop/s"
     write(unit, *)
  endif

  do i = 1, n_timing

     write(a_n_call, ifmt) meas(i)%n_call
     a_time_mean    = " "
     a_mflops_mean  = " " 
     a_mflops_min   = " " 
     a_mflops_max   = " " 
     a_total_gflops = " " 

     if (meas(i)%n_call /= 0) then  ! must be true on all PEs !!

        time_mean = global_sum(meas(i)%total_time) / num_pes()
        write(a_time_mean, ffmt) time_mean
        
        if (n_op(i) /= 0) then
           mflops      = 1d-6 * n_op(i) * volh * meas(i)%n_call
           mflops_mean = mflops / time_mean
           mflops      = mflops / meas(i)%total_time

           mflops_min = global_min(mflops)
           mflops_max = global_max(mflops)
        
           total_gflops = 1e-3 * mflops_mean * num_pes()

           write(a_mflops_mean,  ffmt) mflops_mean / n_thread
           write(a_mflops_min,   ffmt) mflops_min / n_thread
           write(a_mflops_max,   ffmt) mflops_max / n_thread
           write(a_total_gflops, ffmt) total_gflops
        endif
     endif

     if (my_pe() == 0) then
        write(unit, tab_fmt) text(i), a_n_call, a_time_mean, &
                             a_mflops_mean, a_mflops_min, a_mflops_max, &
                             a_total_gflops
     endif
  enddo

  call end(unit, "Timing")
end

!===============================================================================
