#===============================================================================
#
# modules/Makefile
#
#-------------------------------------------------------------------------------
#
# Copyright (C) 2003-2006 Hinnerk Stueben
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
#===============================================================================

include ../Makefile.var

fpp = $(FPP) -I.. $(MYFLAGS)


.SUFFIXES:
.SUFFIXES: .a .o .F90

.F90.o:
	$(fpp) $< > $*.f90
	$(F90) -c $(FFLAGS) $*.f90

MODULES_DIR = .

MODULES = \
	typedef_cksum.o \
	typedef_clover.o \
	typedef_flags.o \
	typedef_hmc.o \
	typedef_para.o \
	module_svn.o \
	module_bqcd.o \
	module_counter.o \
	module_d21.o \
	module_decomp.o \
	module_hmc_forces.o \
	module_lattice.o \
	module_lattice_io.o \
	module_nn.o \
	module_nnpe.o \
	module_offset.o \
	module_thread.o \
	module_vol.o \
	module_conf_info.o \
	module_function_decl.o \
	module_input.o \
	module_switches.o \
	module_p_interface.o \
	module_mre.o

modules: $(MODULES)

fast:
	$(MAKE)

module_input.o: module_input.F90
	$(fpp) $< | sed s/FORTRAN_PERCENT/%/g > $*.f90
	$(F90) -c $(FFLAGS) $*.f90

clean:
	rm -f *.[Tiod] *.f90 *.mod work.pc work.pcl

clobber: clean
