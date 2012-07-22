#===============================================================================
#
# comm/Makefile
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

OBJS_MPI = \
	dotprod.o \
	comm_mpi.o \
	allocate.o \
	field_io_mpi.o \
	pes_mpi.o \
	broadcast_mpi.o \
	reduction_mpi.o \
	seed_mpi.o \
	xbound_mpi.o

OBJS_SHMEM = \
	dotprod.o \
	comm_shmem.o \
	allocate_shmem.o \
	field_io_shmem.o \
	reduction_shmem.o \
	seed_shmem.o \
	xbound_shmem.o

OBJS_SHMEMPI = \
	dotprod.o \
	comm_shmempi.o \
	allocate_shmem.o \
	field_io_mpi.o \
	pes_mpi.o \
	broadcast_mpi.o \
	reduction_mpi.o \
	seed_mpi.o \
	xbound_shmem.o

OBJS_SINGLE_PE = \
	dotprod.o \
	allocate.o \
	comm_single_pe.o \
	field_io_single_pe.o \
	pes_single_pe.o \
	broadcast_single_pe.o \
	reduction_single_pe.o \
	seed_single_pe.o \
	xbound_single_pe.o

$(LIBCOMM):

fast:
	$(FAST_MAKE)

lib_mpi.a: $(OBJS_MPI)
	$(AR) $(ARFLAGS) $@ $(OBJS_MPI)
	$(RANLIB) $@

lib_shmem.a: $(OBJS_SHMEM)
	$(AR) $(ARFLAGS) $@ $(OBJS_SHMEM)
	$(RANLIB) $@

lib_shmempi.a: $(OBJS_SHMEMPI)
	$(AR) $(ARFLAGS) $@ $(OBJS_SHMEMPI)
	$(RANLIB) $@

lib_single_pe.a: $(OBJS_SINGLE_PE)
	$(AR) $(ARFLAGS) $@ $(OBJS_SINGLE_PE)
	$(RANLIB) $@

clobber:
	rm -f *.[Tiod] *.f90 *.mod work.pc work.pcl
	rm -f lib_mpi.a lib_shmem.a lib_shmempi.a lib_single_pe.a
