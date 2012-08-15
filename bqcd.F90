!===============================================================================
!
! bqcd.F90 - main program and read/write of parameters
!
!-------------------------------------------------------------------------------
!
! Copyright (C) 1998-2011 Hinnerk Stueben
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
program bqcd

  use typedef_flags
  use typedef_para
  use module_input
  use module_function_decl
  implicit none
  
  type(type_para)                       :: para
  type(hmc_conf), dimension(MAX_TEMPER) :: conf
  type(type_flags)                      :: flags
  SECONDS                               :: time0, sekunden


  call comm_init()

  time0 = sekunden()  ! start/initialize timer

  TIMING_START(timing_bin_total)

  call get_flags(flags)
  call begin(UREC, "Job")
  call input_read(flags%input)
  call init_para(para, flags)
  call init_counter(para, flags)
  call init_ran(para, flags)
  call init_cooling(input%measure_cooling_list)

  write(*,*) para%NPE(1), para%NPE(2), para%NPE(3), para%NPE(4) 

  call set_fmt_ensemble(para%n_temper)
  call check_fmt(para%run, para%n_temper, para%maxtraj, para%L(4) - 1)

  call init_common(para)
  call init_modules()

  call write_header(para)

  call init_flip_bc()
  call init_cg_para(para%cg_rest, para%cg_maxiter, para%cg_log)
  call init_cg_stat()
  call init_xbound()
  call init_confs(para, conf)

  call check_former(para%n_temper, conf)

#ifdef QUDA_SOLVER
  call comm_set_gridsize(para%NPE) ! awaiting the official QUDA interface for this
  call init_quda(-1) ! Must be after init_para, since this is where topology is read in
#endif

  call mc(para, conf)
  !!call xbound_test()

  call conf_write(.true., para, conf)

  call write_counter(para%maxtraj)
  call write_ran()

  TIMING_STOP(timing_bin_total)

  call write_footer(time0)
  call end(UREC, "Job")

#ifdef QUDA_SOLVER
  call end_quda()
#endif

  call comm_finalize()
end

!-------------------------------------------------------------------------------
subroutine init_para(para, flags)

  ! initialises module_para, module_switches and module_mre

  use typedef_flags
  use typedef_para
  use module_bqcd
  use module_input
  use module_mre
  use module_switches
  implicit none

  type(type_para)  :: para
  type(type_flags) :: flags
  integer          :: i
  logical          :: quenched, dynamical, clover, h_ext

  quenched   = .false.
  dynamical  = .false.
  clover     = .false.
  h_ext      = .false.

  para%run         = input%run
  para%L           = input%lattice
  para%NPE         = input%processes
  para%bc_fermions = input%boundary_conditions_fermions
  para%gamma_index = input%gamma_index
  para%n_temper    = input%ensembles
  para%nstd        = input%tempering_steps_without
  para%nforce      = input%hmc_accept_first
  para%ntraj       = input%mc_steps
  para%maxtraj     = input%mc_total_steps
  para%nsave       = input%mc_save_frequency
  para%c_cg_rest   = input%solver_rest
  para%cg_maxiter  = input%solver_maxiter
  para%cg_log      = input%solver_ignore_no_convergence
  mre_n_vec        = input%solver_mre_vectors

  call check_bc_fermions(para%bc_fermions, para%gamma_index)

  read(para%c_cg_rest, *) para%cg_rest

 if (para%n_temper <= 0) call die("init_para(): n_temper <= 0")
 if (para%n_temper > MAX_TEMPER) call die("init_para(): n_temper > max_temper")
  
  do i = 1, para%n_temper
     para%c_hmc(i)%beta        = input%beta(i)
     para%c_hmc(i)%kappa       = input%kappa(i)
     para%c_hmc(i)%csw         = input%csw(i)
     para%c_hmc(i)%h           = input%h(i)
     para%c_hmc(i)%traj_length = input%hmc_trajectory_length(i)
     para%c_hmc(i)%ntau        = input%hmc_steps(i)
     para%c_hmc(i)%rho         = input%hmc_rho(i)
     para%c_hmc(i)%m_scale     = input%hmc_m_scale(i)
     para%info_file(i)         = input%start_info_file(i)

     read(para%c_hmc(i)%beta,       *) para%hmc(i)%beta
     read(para%c_hmc(i)%kappa,      *) para%hmc(i)%kappa
     read(para%c_hmc(i)%csw,        *) para%hmc(i)%csw
     read(para%c_hmc(i)%h,          *) para%hmc(i)%h
     read(para%c_hmc(i)%traj_length,*) para%hmc(i)%traj_length
     read(para%c_hmc(i)%ntau,       *) para%hmc(i)%ntau
     read(para%c_hmc(i)%rho,        *) para%hmc(i)%rho
     read(para%c_hmc(i)%m_scale,    *) para%hmc(i)%m_scale

     if (para%hmc(i)%kappa == ZERO .and. para%hmc(i)%csw /= ZERO) then
        para%hmc(i)%csw_kappa = para%hmc(i)%csw
        para%c_hmc(i)%csw = "-1 (infinity)"
        para%hmc(i)%csw = -1
     else
        para%hmc(i)%csw_kappa = para%hmc(i)%csw * para%hmc(i)%kappa
        call check_csw(para%hmc(i)%beta,  para%hmc(i)%csw)
     endif

     para%hmc(i)%tau = para%hmc(i)%traj_length / para%hmc(i)%ntau

     write(para%c_hmc(i)%csw_kappa, "(f20.15)") para%hmc(i)%csw_kappa
     write(para%c_hmc(i)%tau,       "(f20.15)") para%hmc(i)%tau

     if (para%hmc(i)%kappa == ZERO .and. para%hmc(i)%csw == ZERO) then
        quenched = .true.
     else
        dynamical = .true.
     endif

     if (para%hmc(i)%csw /= ZERO) clover     = .true.
     if (para%hmc(i)%h   /= ZERO) h_ext      = .true.

     para%hmc(i)%model = input%hmc_model

     if (para%hmc(i)%model == "A" .and. para%hmc(i)%rho /= ZERO) then
        call warn("init_para(): model == A but rho /= 0")
     endif

     if (para%hmc(i)%model /= "A" .and. para%hmc(i)%rho == ZERO) then
        call warn("init_para(): model /= A but rho == 0")
     endif
  enddo

  select case (input%start_configuration)
     case ("hot");  para%start = START_HOT
     case ("cold"); para%start = START_COLD
     case ("file"); para%start = START_FILE
     case default  
        call die("init_para(): start_configuration must be {hot|cold|file}")
  end select

  select case (input%start_random)
     case ("random");  para%seed = -1
     case ("default"); para%seed = 0
     case default;     read(input%start_random, *) para%seed
  end select

  select case (input%tempering_swap_sequence)
     case ("up");     para%swap_seq = SWAP_UP
     case ("down");   para%swap_seq = SWAP_DOWN
     case ("random"); para%swap_seq = SWAP_RANDOM
     case default
       call die("init_para(): tempering_swap_sequence must be {up|down|random}")
  end select

 if (quenched .and. dynamical) call die("init_para(): quenched/dynamical mixed")

  if (para%nforce < 0) call die("init_para(): nforce < 0")

  if (flags%continuation_job) para%start = START_CONT

  
  switches%quenched   = quenched
  switches%dynamical  = dynamical
  switches%clover     = clover
  switches%h_ext      = h_ext
  switches%hasenbusch = (input%hmc_model /= "A")

  if (quenched) switches%hasenbusch = .false.

  switches%tempering             = .false.
  switches%measure_polyakov_loop = .false.
  switches%measure_traces        = .false.

 if (input%ensembles             >  1) switches%tempering             = .true.
 if (input%measure_polyakov_loop /= 0) switches%measure_polyakov_loop = .true.
 if (input%measure_traces        /= 0) switches%measure_traces        = .true.

  if (input%hmc_test == 0) then
     switches%hmc_test = .false.
  else
     switches%hmc_test = .true.
  endif

end

!-------------------------------------------------------------------------------
subroutine init_counter(para, flags)

  use typedef_flags
  use typedef_para
  use module_counter
  use module_function_decl
  implicit none

  type(type_para)    :: para
  type(type_flags)   :: flags
  FILENAME, external :: count_file, stop_file

  if (f_exist(stop_file())) then
     call die("init_counter(): found stop file " // stop_file())
  endif

  counter%run = para%run
  counter%j_traj = 0

  if (flags%continuation_job) then
     open(UCOUNT, file = count_file(), action = "read", status = "old")
     read(UCOUNT, *) counter%run
     read(UCOUNT, *) counter%job
     read(UCOUNT, *) counter%traj
     close(UCOUNT)

     if (counter%run /= para%run) call die("init_counter(): RUN inconsistent")
     counter%job = counter%job + 1
  else
     counter%run = para%run
     counter%job = 1
     counter%traj = -para%nforce
  endif

end

!-------------------------------------------------------------------------------
subroutine write_counter(maxtraj)

  use module_counter
  use module_function_decl
  implicit none

  integer            :: maxtraj
  FILENAME, external :: count_file, stop_file

  if (my_pe() /= 0) return

  open(UCOUNT, file = count_file(), action = "write")
  write(UCOUNT, *) counter%run, " run"
  write(UCOUNT, *) counter%job, " job"
  write(UCOUNT, *) counter%traj, " traj"
  close(UCOUNT)

  if (counter%traj >= maxtraj) then
     open(UCOUNT, file = stop_file(), status = "unknown")
     close(UCOUNT)
  endif

end

!-------------------------------------------------------------------------------
subroutine write_header(para)

  use typedef_para
  use module_bqcd
  use module_counter
  use module_function_decl
  use module_input
  use module_mre
  use module_thread
  implicit none

  type(type_para)     :: para
  integer             :: i
  character(len = 50) :: fmt
  character(len = 4), external :: format_ensemble

  if (my_pe() == 0) then

     fmt = "(1x,a," // format_ensemble() // ",2a)"

     call begin(UREC, "Header")

    if (input%comment /= "") then
     write(UREC, 405) "Comment", trim(input%comment)
    endif

     write(UREC, 400) "Program", prog_name, prog_version // prog_revision
     write(UREC,   *) "Version_of_D ", version_of_d()
     write(UREC,   *) "Communication ", trim(comm_method())
     write(UREC,   *) "Run ", para%run
     write(UREC,   *) "Job ", counter%job
     write(UREC, 405) "Host", rechner()
     write(UREC, 400) "Date", datum(), uhrzeit()
     write(UREC, 410) "L          ", para%L
     write(UREC, 410) "NPE        ", para%NPE
     write(UREC, 410) "bc_fermions", para%bc_fermions
     write(UREC, 410) "gamma_index", para%gamma_index


     write(UREC,   *) "Threads ", n_thread
     write(UREC,   *) "Start   ", para%start

     if (para%start == START_FILE) then
        do i = 1, para%n_temper
           write(UREC, fmt) "StartConf_", i, " ", trim(para%info_file(i))
        enddo
     endif

     write(UREC,   *) "Seed    ", para%seed
     write(UREC,   *) "Swap_seq", para%swap_seq
     write(UREC,   *) "N_force ", para%nforce
     write(UREC,   *) "N_traj  ", para%ntraj
     write(UREC,   *) "N_save  ", para%nsave
     write(UREC,   *) "N_temper", para%n_temper

     do i = 1, para%n_temper
        write(UREC, fmt) "beta_", i, "        ", trim(para%c_hmc(i)%beta)
        write(UREC, fmt) "kappa_", i, "       ", trim(para%c_hmc(i)%kappa)
        write(UREC, fmt) "csw_", i, "         ", trim(para%c_hmc(i)%csw)
        write(UREC, fmt) "csw_kappa_", i, "   ", trim(para%c_hmc(i)%csw_kappa)
        write(UREC, fmt) "h_", i, "           ", trim(para%c_hmc(i)%h)
        write(UREC, fmt) "tau_", i, "         ", trim(para%c_hmc(i)%tau)
        write(UREC, fmt) "N_tau_", i, "       ", trim(para%c_hmc(i)%ntau)
        write(UREC, fmt) "traj_length_", i, " ", trim(para%c_hmc(i)%traj_length)
        write(UREC, fmt) "rho_", i, "         ", trim(para%c_hmc(i)%rho)
        write(UREC, fmt) "m_scale_", i, "     ", trim(para%c_hmc(i)%m_scale)
     enddo

     write(UREC,   *) "HMC_model ", para%hmc(1)%model
     write(UREC,   *) "REAL_kind ", RKIND
     write(UREC, 405) "CG_rest ", trim(para%c_cg_rest)
     write(UREC,   *) "MRE_vectors ", mre_n_vec

     call end(UREC, "Header")

400  format (3(1x,a))
405  format (2(1x,a))
410  format (1x,a,4i5)

  endif

end

!-------------------------------------------------------------------------------
subroutine write_footer(time0)

  use module_function_decl
  use module_thread
  implicit none

  SEED    :: seed
  SECONDS :: time0, sekunden

  call ranget(seed)

  call begin(UREC, "Footer")

  if (my_pe() == 0) then
     write(UREC, 400) "Date", datum(), uhrzeit()
     write(UREC,   *) "Seed", seed
     write(UREC, 410) "CPU-Time", &
                      sekunden() - time0, "s on", num_pes() * n_thread, "CPUs"
  endif

400  format (3(1x,a))
410  format (1x,a,1x,f8.1,1x,a,1x,i5,1x,a)

  TIMING_WRITE(UREC)

  call end(UREC, "Footer")

end

!-------------------------------------------------------------------------------
subroutine get_flags(flags) ! read command line arguments [and open output file]

  use typedef_cksum
  use typedef_flags
  use module_bqcd
  use module_function_decl
  use module_input
  implicit none

  type(type_flags), intent(out) :: flags

  integer                       :: iarg, length, stat, narg
  integer, external             :: ipxfargc
  character(len = 2)            :: opt

  flags%continuation_job = .false.
  flags%show_version = .false.

  narg = ipxfargc()

  iarg = 1
  do while (iarg <= narg)
     call pxfgetarg(iarg, opt, length, stat)

     if (opt(1:1) == "-") then
        if (length > 2) call usage()

        select case (opt(2:2))
        case ("c")
           flags%continuation_job = .true.
           iarg = iarg + 1
        case ("I")
           call input_dump(6)
           call comm_finalize()
           stop
        case ("V")
           flags%show_version = .true.
           iarg = iarg + 1
        case default
           call usage
        end select
     else
        exit
     endif
  enddo
           
  if (flags%show_version) then
     call input_init()  ! some input defaults are needed in version()
     call version()
     call comm_finalize()
     stop
  endif

  call take_arg(iarg, flags%input, narg)

  if (iarg <= narg) then
     call take_arg(iarg, flags%output, narg)
     if (my_pe() == 0) then
        close(UREC)

        if (flags%continuation_job) then
           open(UREC, file = flags%output, &
                   action = "write", position = "append", status = "unknown")
        else
           open(UREC, file = flags%output,  action = "write", status = "unknown")
        endif
     endif
  endif
  if (narg >= iarg) call usage

CONTAINS

  subroutine usage()
    implicit none
    call die("Usage: " // prog_name // " [-c] [-I] [-V] input [output]")
  end subroutine usage

  subroutine version()
    implicit none
    logical, external :: is_big_endian

    if (my_pe() == 0) then
       write(6,*) "This is ", prog_name, " ", prog_version, prog_revision
       write(6,*) "   input format:    ", input_version
       write(6,*) "   conf info format:", conf_info_version
       write(6,*) "   MAX_TEMPER:      ", MAX_TEMPER
       write(6,*) "   REAL kind:       ", RKIND
       write(6,*) "   Version of D:    ", version_of_d()
       write(6,*) "   Communication:   ", trim(comm_method())
       if (is_big_endian()) then
       write(6,*) "   Endianness:      big endian"
       else
       write(6,*) "   Endianness:      little endian"
       endif

    endif
  end subroutine version

  subroutine take_arg(iarg, arg, narg)
    implicit none
    integer, intent(inout)          :: iarg
    character(len = *), intent(out) :: arg
    integer, intent(in)             :: narg
    integer                         :: length, stat

    if (iarg > narg) call usage
    call pxfgetarg(iarg, arg, length, stat)
    if (length > len(arg)) then
       call die("get_flags(): " // arg // ": argument too long")
    endif
    iarg = iarg + 1
  end subroutine take_arg

end

!===============================================================================
