#===============================================================================
#
# clover/Makefile
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
#===============================================================================

include ../Makefile.var

ifdef FPP2
	fpp = $(FPP2)
else
	fpp = $(FPP)
endif


.SUFFIXES:
.SUFFIXES: .a .o .F90

.F90.o:
	$(fpp) -I.. $(MYFLAGS) $< > $*.f90
	$(F90) -c $(FFLAGS) $*.f90

MODULES_DIR = ../modules

OBJS = \
	clover_action.o \
	clover_allocate.o \
	clover_bsa.o \
	clover_d.o \
	clover_f_mu_nu.o \
	clover_init.o \
	clover_inv.o \
	clover_mult_a.o \
	clover_mult_ao.o \
	clover_mult_b.o \
	clover_t_init.o \
	clover_ts.o \
	clover_uuu.o \
	clover_uuuu.o

OBJS_CTEST = \
	ctest.o \
	clover_inv.o \
	clover_mult_a.o \
	clover_mult_ao.o \
	clover_mult_b.o

$(LIBCLOVER):

libclover.a: $(OBJS)
	$(AR) $(ARFLAGS) $@ $(OBJS)
	$(RANLIB) $@

fast:
	$(FAST_MAKE)

ctest: $(OBJS_CTEST)
	f90 -o $@ $(OBJS_CTEST)

clean:
	rm -f *.[Tiod] *.f90 *.mod core work.pc work.pcl

clobber: clean
	rm -f libclover.a ctest
