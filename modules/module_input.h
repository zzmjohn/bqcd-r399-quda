#/*
#===============================================================================
#
# module_input.h
#
#-------------------------------------------------------------------------------
#
# Copyright (C) 2003-2010 Hinnerk Stueben
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

INPUT_INPUT(run, integer, 0)
INPUT_INPUT(comment, character(comment_len), "")

INPUT_INPUT(lattice, INPUT_ARRAY_DIM, (/4, 4, 4, 4/))
INPUT_INPUT(processes, INPUT_ARRAY_DIM, (/1, 1, 1, 1/))
INPUT_INPUT(boundary_conditions_fermions, INPUT_ARRAY_DIM, (/1, 1, 1, -1/))

INPUT_INPUT(ensembles, integer, INPUT_DEFAULT_ENSEMBLES)

INPUT_INPUT(beta, INPUT_ARRAY_ENSEMBLES, "0.0")
INPUT_INPUT(kappa, INPUT_ARRAY_ENSEMBLES, "0.0")
INPUT_INPUT(csw, INPUT_ARRAY_ENSEMBLES, "0.0")
INPUT_INPUT(h, INPUT_ARRAY_ENSEMBLES, "0.0")

INPUT_INPUT(tempering_swap_sequence, character(word_len), "random")
INPUT_INPUT(tempering_steps_without, integer, 0)

INPUT_INPUT(hmc_model, character, "A")
INPUT_INPUT(hmc_trajectory_length, INPUT_ARRAY_ENSEMBLES, "1")
INPUT_INPUT(hmc_steps, INPUT_ARRAY_ENSEMBLES, "0")
INPUT_INPUT(hmc_rho, INPUT_ARRAY_ENSEMBLES, "0.0")
INPUT_INPUT(hmc_m_scale, INPUT_ARRAY_ENSEMBLES, "1")
INPUT_INPUT(hmc_accept_first, integer, 0)
INPUT_INPUT(hmc_test, integer, 0)

INPUT_INPUT(start_configuration, character(word_len), "cold")
INPUT_INPUT(start_info_file, INPUT_ARRAY_ENSEMBLES, "")
INPUT_INPUT(start_random, character(para_len), "default")

INPUT_INPUT(mc_total_steps, integer, 1)
INPUT_INPUT(mc_steps, integer, 1)
INPUT_INPUT(mc_save_frequency, integer, 0)

INPUT_INPUT(solver_rest, character(para_len), "1e-8")
INPUT_INPUT(solver_maxiter, integer, 100)
INPUT_INPUT(solver_ignore_no_convergence, integer, 0)
INPUT_INPUT(solver_mre_vectors, integer, 0)

INPUT_INPUT(measure_cooling_list, FILENAME, "")
INPUT_INPUT(measure_polyakov_loop, integer, 0)
INPUT_INPUT(measure_traces, integer, 0)
