!===============================================================================
!
! invert_quda.F90
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
# include "quda_defs.h"

!-------------------------------------------------------------------------------

subroutine norm2(nrm2, x)
  use module_vol
  use module_function_decl
  implicit none
  real(8), intent(out) :: nrm2
  SPINCOL_OVERINDEXED, intent(in) :: x
  integer :: i

  nrm2 = ZERO
  !$omp parallel do reduction(+: rtrold)
  do i = 1, size_sc_field
     nrm2 = nrm2 + x(i)**2
  end do
  nrm2 = global_sum(nrm2)
end subroutine norm2


!-------------------------------------------------------------------------------
subroutine quda_solver(matrix_mult, x, b, para, conf, iterations) 

  ! calls the QUDA solver

  use      module_cg
  use      module_function_decl
  use      module_p_interface
  use      module_vol
  use      typedef_hmc
  use      typedef_quda
  implicit none

  external                          :: matrix_mult
  SPINCOL_OVERINDEXED,  intent(out) :: x
  SPINCOL_OVERINDEXED,  intent(in)  :: b
  type(hmc_para),       intent(in)  :: para
  type(hmc_conf),       intent(in)  :: conf
  integer,               intent(out) :: iterations
  integer mu, p, v

  P_SPINCOL_OVERINDEXED, save :: r, tmp

  integer       :: i, niter
  integer(4) :: device
  REAL :: rtr
  character(72) :: msg
  REAL :: rtr2
  type(quda_gauge_param) :: gauge_param
  type(quda_invert_param) :: invert_param

  ALLOCATE_SC_OVERINDEXED(r)
  ALLOCATE_SC_OVERINDEXED(tmp)

  call init_quda_gauge_param(gauge_param)
  call init_quda_invert_param(invert_param, para)

  device = 0
  call init_quda(device) ! FIXME need to sort out MPI topology

  call load_gauge_quda(conf%u, gauge_param)
  if (invert_param%dslash_type.eq.QUDA_CLOVER_WILSON_DSLASH) then
     call load_clover_quda(conf%a, conf%i, invert_param)
  end if

  TIMING_START(timing_bin_cg)

  !do i=1, size_sc_field
  !   x(i) = 0.0
  !end do
  !x(1) = 1.0

  !call norm2(rtr, x)
  !write(*,*) "BQCD src = ", rtr

  ! will use this to check result on host
  !call matrix_mult(r, x, para, conf)

  !call dslash_quda(tmp, x, invert_param, QUDA_ODD_PARITY) 
  !call mat_dag_mat_quda(tmp, x, invert_param) 

  !do i=1, size_sc_field
  !   if ( (r(i).ne.0.0) .or. (tmp(i).ne.0.0) ) then
  !      write(*,*) i, " " , r(i), " ", tmp(i)
  !   end if
  !end do
  
  !call norm2(rtr, r)
  !call norm2(rtr2, tmp)
  !write(*,*) "BQCD mat vec = ", rtr, "QUDA mat vec = ", rtr2

  !call die("dying")

  call invert_quda(x, b, invert_param)

  TIMING_STOP(timing_bin_cg)

  call norm2(rtr, x)
  write(*,*) "BQCD solution = ", rtr

  ! will use this to check result on host
  call matrix_mult(r, x, para, conf)

  do i = 1, size_sc_field
     r(i) = b(i) - r(i)
  end do
  call norm2(rtr, r)
  write(*,*) "BQCD Residual = ", rtr

  if (rtr <= cg_para%rest) goto 9999

  if (cg_para%log /= 2) then
     write(msg, *) "cg(): no convergence; rtr = ", rtr 
     call die(msg)
  endif

9999 continue

  niter = invert_param%iter

  cg_stat%ncall = cg_stat%ncall + 1
  cg_stat%niter = niter
  cg_stat%niter_max = max(cg_stat%niter_max, niter)
  cg_stat%niter_tot = cg_stat%niter_tot + niter
  cg_iterations_total = cg_iterations_total + niter

  iterations = niter

  if (invert_param%dslash_type.eq.QUDA_CLOVER_WILSON_DSLASH) then
     call free_clover_quda()
  endif

  call free_gauge_quda()
  call end_quda()

end
 
!-------------------------------------------------------------------------------
! Fill the quda_gauge_param
!-------------------------------------------------------------------------------
subroutine init_quda_gauge_param(gauge_param)

  use      module_lattice
  use      typedef_quda
  implicit none
  integer i
  integer x_face, y_face, z_face, t_face
  type(quda_gauge_param), intent(out) :: gauge_param
  character(72) :: msg

  call new_quda_gauge_param(gauge_param)

  gauge_param%x(1) = NX  
  gauge_param%x(2) = NY  
  gauge_param%x(3) = NZ  
  gauge_param%x(4) = NT

  gauge_param%anisotropy = 1.0d0
  gauge_param%tadpole_coeff = ZERO !Used for staggered only

  gauge_param%link_type = QUDA_SU3_LINKS
  gauge_param%gauge_order = QUDA_BQCD_GAUGE_ORDER

  do i=1,3
     if (bc_fermions(i).ne.1) then
        write(msg, *) "init_quda_gauge_param(): bc not supported in QUDA" 
        call die(msg)
     end if
  end do
  gauge_param%t_boundary = bc_fermions(4)

  gauge_param%cpu_prec = RKIND
  gauge_param%cuda_prec = gauge_param%cpu_prec ! match CUDA precision with BQCD precision
  gauge_param%reconstruct = QUDA_RECONSTRUCT_12
  gauge_param%cuda_prec_sloppy = gauge_param%cuda_prec
  gauge_param%reconstruct_sloppy = gauge_param%reconstruct
  gauge_param%cuda_prec_precondition = gauge_param%cuda_prec_sloppy
  gauge_param%reconstruct_precondition = gauge_param%reconstruct_sloppy
  gauge_param%gauge_fix = QUDA_GAUGE_FIXED_NO

  ! compute the required pad size (this is where the backward links on the zero boundary are stored
  x_face = NX*NY*NZ/2
  y_face = NX*NZ*NT/2
  z_face = NX*NY*NT/2
  t_face = NX*NY*NZ/2
  gauge_param%ga_pad = 0 !max(max(x_face, y_face), max(z_face, t_face))

  ! these are all set to zero and are not used by the linear solver
  gauge_param%site_ga_pad = 0
  gauge_param%staple_pad = 0
  gauge_param%llfat_ga_pad = 0
  gauge_param%mom_ga_pad = 0
  gauge_param%preserve_gauge = 0
end subroutine init_quda_gauge_param

!===============================================================================

!-------------------------------------------------------------------------------
! Fill the quda_gauge_param.  We are only setting relevent parameters
! for the current BQCD benchmark.
! -------------------------------------------------------------------------------
subroutine init_quda_invert_param(invert_param, hmc_param)

  use      module_lattice
  use      module_cg
  use      typedef_hmc
  use      typedef_quda
  implicit none

  type(quda_invert_param), intent(out) :: invert_param
  type(hmc_para),       intent(in)  :: hmc_param

  call new_quda_invert_param(invert_param)

  ! The input and output spinor field both reside on the host
  invert_param%input_location = QUDA_CPU_FIELD_LOCATION
  invert_param%output_location = QUDA_CPU_FIELD_LOCATION
     
  if (hmc_param%csw_kappa.ne.0) then
     invert_param%dslash_type = QUDA_CLOVER_WILSON_DSLASH
     invert_param%matpc_type = QUDA_MATPC_EVEN_EVEN_ASYMMETRIC
  else 
     invert_param%dslash_type = QUDA_WILSON_DSLASH
     invert_param%matpc_type = QUDA_MATPC_EVEN_EVEN
  endif
  invert_param%inv_type = QUDA_CG_INVERTER
  invert_param%inv_type_precondition = QUDA_INVALID_INVERTER
     
  invert_param%kappa = hmc_param%kappa !* 2.0d0
     
  invert_param%tol = cg_para%rest ! FIXME: QUDA uses relative, BQCD uses absolute
  invert_param%maxiter = cg_para%maxiter
  invert_param%reliable_delta = 1e-20 ! no reliable updates
  
  invert_param%solution_type = QUDA_MATPCDAG_MATPC_SOLUTION
  invert_param%solve_type = QUDA_NORMEQ_PC_SOLVE
  invert_param%dagger = QUDA_DAG_NO
  invert_param%mass_normalization = QUDA_KAPPA_NORMALIZATION
     
  invert_param%preserve_source = QUDA_PRESERVE_SOURCE_YES
     
  invert_param%cpu_prec = RKIND
  invert_param%cuda_prec = invert_param%cpu_prec
  invert_param%cuda_prec_sloppy = invert_param%cuda_prec
  invert_param%cuda_prec_precondition = invert_param%cuda_prec_sloppy
     
  invert_param%dirac_order = QUDA_QDP_DIRAC_ORDER
     
  ! BQCD uses UKQCD basis
  invert_param%gamma_basis = QUDA_UKQCD_GAMMA_BASIS !FIXME
     
  invert_param%clover_cpu_prec = RKIND
  invert_param%clover_cuda_prec = invert_param%clover_cpu_prec
  invert_param%clover_cuda_prec_sloppy = invert_param%clover_cuda_prec
  invert_param%clover_cuda_prec_precondition = invert_param%clover_cuda_prec_precondition
  
  invert_param%clover_order = QUDA_PACKED_CLOVER_ORDER ! FIXME
  invert_param%use_init_guess = QUDA_USE_INIT_GUESS_YES
  
  invert_param%verbosity = QUDA_SUMMARIZE
  
  invert_param%sp_pad = 0
  invert_param%cl_pad = 0
  
  invert_param%tune = QUDA_TUNE_NO
     
end subroutine init_quda_invert_param
