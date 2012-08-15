#/*
#===============================================================================
#
# defs.h
#
#-------------------------------------------------------------------------------
#
# Copyright (C) 1998-2006 Hinnerk Stueben
#
# This file is part of BQCD -- Berlin Quantum ChromoDynamics program
#
# BQCD is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# BQCD is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with BQCD.  If not, see <http://www.gnu.org/licenses/>.
#
#-------------------------------------------------------------------------------
#*/
#ifndef BQCD_DEFS_H
#define BQCD_DEFS_H

# define QUDA_SOLVER 1

# define MAX_TEMPER 50

# define RKIND  8

# define BQCD_REAL mpi_real8

# define BQCD_CHECK_SUM mpi_integer8
# define BQCD_SEED mpi_integer8

# define CHECK_SUM integer(8)
# define SEED integer(8)
# define SECONDS real(8)
# define COMM_METHOD character(40)

# define RECL_UNIT 1

# define DIM    4
# define NCOL   3
# define NDIRAC 4
# define NGEN   8
# define EVEN   0
# define ODD    1
# define FWD    0
# define BWD    1

# define SIZE_COMPLEX 2

# define REAL real(RKIND)
# define INTEGER integer(4)
# define COMPLEX complex(RKIND)

# define SU3 COMPLEX, dimension (NCOL, NCOL)
# define GENERATOR REAL, dimension (NGEN)

# define GAUGE_FIELD_IO  COMPLEX, dimension(NCOL, NCOL-1, DIM, 0:NX-1, 0:NY-1, 0:NZ-1, 0:NT-1)
# define SPINCOL_FIELD_IO  COMPLEX, dimension(NDIRAC, NCOL, 0:NXH-1, 0:NY-1, 0:NZ-1, 0:NT-1)

# define SU3_FIELD COMPLEX, dimension (NCOL, NCOL, volh_tot)
# define GAUGE_FIELD COMPLEX, dimension (NCOL, NCOL, volh_tot, EVEN:ODD, DIM)
# define GENERATOR_FIELD REAL, dimension (NGEN, volh_tot, EVEN:ODD, DIM)
# define SPINCOL_FIELD COMPLEX, dimension (NDIRAC, NCOL, volh_tot)
# define SC2_FIELD COMPLEX, dimension(2, NCOL, volh_tot, DIM, FWD:BWD)
# define CLOVER_FIELD_A type(type_clover_a), dimension(2, volh, EVEN:ODD)
# define CLOVER_FIELD_B type(type_clover_b), dimension(2, volh, EVEN:ODD)
# define CLOVER_FIELD_C COMPLEX, dimension(NDIRAC, NCOL, NDIRAC, NCOL, volh)

# define P_GAUGE_FIELD COMPLEX, dimension(:, :, :, :, :), pointer
# define P_GAUGE_FIELD_IO COMPLEX, dimension(:, :, :, :, :, :, :), pointer
# define P_GENERATOR_FIELD REAL, dimension(:, :, :, :), pointer
# define P_SPINCOL_FIELD COMPLEX, dimension(:, :, :), pointer
# define P_SPINCOL_FIELD_IO COMPLEX, dimension(:, :, :, :, :, :), pointer
# define P_SC2_FIELD COMPLEX, dimension(:, :, :, :, :), pointer
# define P_CLOVER_FIELD_A type(type_clover_a), dimension(:, :, :), pointer
# define P_CLOVER_FIELD_B type(type_clover_b), dimension(:, :, :), pointer

# define SPINCOL_OVERINDEXED REAL, dimension(SIZE_COMPLEX*NDIRAC*NCOL*volh_tot)
# define P_SPINCOL_OVERINDEXED REAL, dimension(:), pointer

# define FILENAME character(len=80)
# define FILENAME_FORMAT character(len=80)

# define Re(z) real(z)
# define Im(z) aimag(z)

# define CAT(A, B) A ## B
# define STRCAT(A, B) CAT(A, B)
# define STRCAT3(A, B, C) STRCAT(STRCAT(A, B), C)

# define PI    STRCAT(3.1415926535897931_, RKIND)
# define TWOPI STRCAT(6.2831853071795862_, RKIND)
# define SQRT3 STRCAT(1.7320508075688772_, RKIND)

# define ZERO  STRCAT(0.0_, RKIND)
# define ONE   STRCAT(1.0_, RKIND)
# define TWO   STRCAT(2.0_, RKIND)
# define THREE STRCAT(3.0_, RKIND)
# define FOUR  STRCAT(4.0_, RKIND)
# define SIX   STRCAT(6.0_, RKIND)
# define EIGHT STRCAT(8.0_, RKIND)

# define HALF   STRCAT(0.5_, RKIND)
# define EIGHTH STRCAT(0.125_, RKIND)

# define timing_bin_d_xf            1
# define timing_bin_d_xb            2
# define timing_bin_d_yf            3
# define timing_bin_d_yb            4
# define timing_bin_d_zf            5
# define timing_bin_d_zb            6
# define timing_bin_d_t             7
# define timing_bin_d               8
# define timing_bin_mtdagmt         9
# define timing_bin_global_sum     10
# define timing_bin_global_sum_vec 11
# define timing_bin_sc_zero        12
# define timing_bin_sc_copy        13
# define timing_bin_sc_scale       14
# define timing_bin_sc_norm2       15
# define timing_bin_sc_dot         16
# define timing_bin_sc_axpy        17
# define timing_bin_sc_xpby        18
# define timing_bin_sc_axpby       19
# define timing_bin_sc_cdotc       20
# define timing_bin_sc_caxpy       21
# define timing_bin_sc_caxpy2      22
# define timing_bin_sc_cax2        23
# define timing_bin_cg             24
# define timing_bin_hmc_init_p     25
# define timing_bin_hmc_u          26
# define timing_bin_dsg            27
# define timing_bin_dsf            28
# define timing_bin_clover_init    29
# define timing_bin_clover_mult_a  30
# define timing_bin_clover_mult_ao 31
# define timing_bin_clover_mult_b  32
# define timing_bin_clover_dsd     33
# define timing_bin_clover_dsf     34
# define timing_bin_hmc            35
# define timing_bin_plaq           36
# define timing_bin_cooling        37
# define timing_bin_u_read         38
# define timing_bin_u_write        39
# define timing_bin_total          40

# define timing_bin_hmc_init       41
# define timing_bin_hmc_momenta    42
# define timing_bin_hmc_init_phi   43
# define timing_bin_hmc_h_old      44
# define timing_bin_hmc_backup     45
# define timing_bin_hmc_half_step0 46
# define timing_bin_hmc_half_step1 47
# define timing_bin_hmc_xbound_g   48
# define timing_bin_hmc_steps      49
# define timing_bin_hmc_h_new      50
# define timing_bin_hmc_rest       51

# define timing_bin_h_mult_a       52
# define timing_bin_h_mult_b       53
# define timing_bin_h_mult_c       54

# define timing_bin_sc2_projection 55

# define timing_bin_d_dag_xf timing_bin_d_xf 
# define timing_bin_d_dag_xb timing_bin_d_xb 
# define timing_bin_d_dag_yf timing_bin_d_yf 
# define timing_bin_d_dag_yb timing_bin_d_yb 
# define timing_bin_d_dag_zf timing_bin_d_zf 
# define timing_bin_d_dag_zb timing_bin_d_zb 
# define timing_bin_d_dag_t  timing_bin_d_t  
# define timing_bin_d_dag    timing_bin_d    

#ifdef TIMING

# define TIMING_START(bin) call timing_start(bin)
# define TIMING_STOP(bin) call timing_stop(bin)
# define TIMING_WRITE(unit) call timing_write(unit)

#else

# define TIMING_START(bin)
# define TIMING_STOP(bin)
# define TIMING_WRITE(unit)

#endif


# define STDERR 0
# define UINPUT 1
# define UCONF 2
# define URAN 3
# define UCOUNT 4
# define UREC 6
# define UINFO 7
# define ULIST 8
# define UDIAG 99

# define START_HOT 0
# define START_COLD 1
# define START_CONT 2
# define START_FILE 3

# define SWAP_DOWN -1
# define SWAP_RANDOM 0
# define SWAP_UP 1

# define HMC_TEST_FORWARDS 1
# define HMC_TEST_NONE 0
# define HMC_TEST_BACKWARDS -1

# define PUTSTR(unit, str) if (my_pe() == 0) write(unit,*) str
# define PUTVAL(unit, val) if (my_pe() == 0) write(unit,*) #val, ": ", val

# define DIAGSTR(str) write(UDIAG,*) str
# define DIAGVAL(val) write(UDIAG,*) #val, ": ", val

# define ALLOCATE_G_FIELD(x) if (.not. associated(x)) call allocate_g_field(x)
# define ALLOCATE_G_FIELD_IO(x) if (.not. associated(x)) call allocate_g_field_io(x)
# define ALLOCATE_GEN_FIELD(x) if (.not. associated(x)) call allocate_gen_field(x)
# define ALLOCATE_SC_FIELD(x) if (.not. associated(x)) call allocate_sc_field(x)
# define ALLOCATE_SC_FIELD_IO(x) if (.not. associated(x)) call allocate_sc_field_io(x)
# define ALLOCATE_SC_OVERINDEXED(x) if (.not. associated(x)) call allocate_sc_overindexed(x)
# define ALLOCATE_SC2_FIELD(x) if (.not. associated(x)) call allocate_sc2_field(x)

# define ASSERT(condition) if (.not. (condition)) call assertion_failed(__FILE__, __LINE__, #condition)

#endif
