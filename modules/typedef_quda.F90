!===============================================================================
!
! typedef_quda.F90
!
!-------------------------------------------------------------------------------
!
! Copyright (C) 2012 M Clark
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
# include "quda_defs.h"

!-------------------------------------------------------------------------------

#define QUDA_MAX_DIM 6
#define QUDA_MAX_MULTI_SHIFT 32

module typedef_quda

  ! This corresponds to the QudaGaugeParam struct in quda.h
  type quda_gauge_param

     integer(4), dimension(4) :: x

     real(8) :: anisotropy    !Used for Wilson and Wilson-clover
     real(8) :: tadpole_coeff !Used for staggered only

     QudaLinkType :: link_type
     QudaGaugeFieldOrder :: gauge_order
     QudaTboundary :: t_boundary
     QudaPrecision :: cpu_prec
     QudaPrecision :: cuda_prec
     QudaReconstructType :: reconstruct
     QudaPrecision :: cuda_prec_sloppy
     QudaReconstructType :: reconstruct_sloppy
     QudaPrecision :: cuda_prec_precondition
     QudaReconstructType :: reconstruct_precondition
     QudaGaugeFixed :: gauge_fix

     integer(4) :: ga_pad

     integer(4) :: site_ga_pad ! Used by link fattening and the gauge and fermion forces

     integer(4) :: staple_pad   ! Used by link fattening
     integer(4) :: llfat_ga_pad ! Used by link fattening
     integer(4) :: mom_ga_pad   ! Used by the gauge and fermion forces
     real(8) :: gauge_gib

     integer(4) :: preserve_gauge ! Used by link fattening
    
  end type quda_gauge_param

  ! This module corresponds to the QudaInvertParam struct in quda.h
  type quda_invert_param
     
     QudaFieldLocation :: input_location  ! The location of the input field
     QudaFieldLocation :: output_location ! The location of the output field 
     
     QudaDslashType :: dslash_type
     QudaInverterType :: inv_type
     
     real(8) :: mass  ! Used for staggered only 
     real(8) :: kappa ! Used for Wilson and Wilson-clover 
     
     real(8) :: m5    ! Domain wall height 
     integer(4) :: Ls       ! Extent of the 5th dimension (for domain wall) 
     
     real(8) :: mu    ! Twisted mass parameter 
     QudaTwistFlavorType :: twist_flavor  ! Twisted mass flavor 
     
     real(8) :: tol
     integer(4) :: maxiter
     real(8) :: reliable_delta ! Reliable update tolerance 
     
     integer(4) :: num_offset ! Number of offsets in the multi-shift solver 
     
     real(8), dimension(QUDA_MAX_MULTI_SHIFT) :: offset ! Offsets for multi-shift solver 
     real(8), dimension(QUDA_MAX_MULTI_SHIFT) :: tol_offset ! Solver tolerance for each offset 
     
     QudaSolutionType :: solution_type  ! Type of system to solve 
     QudaSolveType :: solve_type        ! How to solve it 
     QudaMatPCType :: matpc_type
     QudaDagType :: dagger
     QudaMassNormalization :: mass_normalization
     
     QudaPreserveSource :: preserve_source
     
     QudaPrecision :: cpu_prec
     QudaPrecision :: cuda_prec
     QudaPrecision :: cuda_prec_sloppy
     QudaPrecision :: cuda_prec_precondition
     
     QudaDiracFieldOrder :: dirac_order
     
     ! Gamma basis of the input and output host fields 
     QudaGammaBasis :: gamma_basis
     
     QudaPrecision :: clover_cpu_prec
     QudaPrecision :: clover_cuda_prec
     QudaPrecision :: clover_cuda_prec_sloppy
     QudaPrecision :: clover_cuda_prec_precondition
     
     QudaCloverFieldOrder :: clover_order
     QudaUseInitGuess :: use_init_guess
     
     QudaVerbosity :: verbosity    
     
     integer(4) :: sp_pad
     integer(4) :: cl_pad
     
     integer(4) :: iter
     real(8) :: spinor_gib
     real(8) :: clover_gib
     real(8) :: gflops
     real(8) :: secs
     
     ! Enable auto-tuning? 
     QudaTune :: tune
     
     ! Maximum size of Krylov space used by solver 
     integer(4) :: gcr_nkrylov
     
     ! The following parameters are related to the domain-decomposed preconditioner.
     
     ! The inner Krylov solver used in the preconditioner.  Set to
     ! QUDA_INVALID_INVERTER to disable the preconditioner entirely.
     QudaInverterType :: inv_type_precondition
     
     ! Verbosity of the inner Krylov solver 
     QudaVerbosity :: verbosity_precondition
     
     ! Tolerance in the inner solver 
     real(8) :: tol_precondition
     
     ! Maximum number of iterations allowed in the inner solver 
     integer(4) :: maxiter_precondition
     
     ! Precision used in the inner solver 
     QudaPrecision :: prec_precondition
     
     ! Relaxation parameter used in GCR-DD (default = 1.0) 
     real(8) :: omega
     
     ! Number of preconditioner cycles to perform per iteration 
     integer(4) :: precondition_cycle
     
     ! Whether to use additive or multiplicative Schwarz preconditioning 
     QudaSchwarzType :: schwarz_type
     
  end type quda_invert_param
   
end module typedef_quda
!===============================================================================
