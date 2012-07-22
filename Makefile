#===============================================================================
#
# Makefile
#
#-------------------------------------------------------------------------------
#
# Copyright (C) 1998-2012 Hinnerk Stueben
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

include Makefile.var

ifdef FPP2
	fpp = $(FPP2)
else
	fpp = $(FPP)
endif


.SUFFIXES:
.SUFFIXES: .o .F90 .c

.F90.f90:
	$(fpp) $(MYFLAGS) $< > $@

.F90.o:
	$(fpp) $(MYFLAGS) $< > $*.f90
	$(F90) -c $(FFLAGS) $*.f90

MODULES_DIR  = modules

MODULES = modules/*.o

OBJS = \
	action.o \
	cg.o \
	checks.o \
	$(CKSUM_O) \
	conf.o \
	conf_info.o \
	cooling.o \
	dsd.o \
	dsg.o \
	dsf.o \
	dsf1.o \
	dsf2.o \
	files.o \
	flip_bc.o \
	hmc.o \
	hmc_init_p.o \
	hmc_init_phi.o \
	hmc_integrator.o \
	hmc_forces.o \
	hmc_leap_frog.o \
	hmc_test.o \
	hmc_u.o \
	h_mult.o \
	index.o \
	index2.o \
	init_common.o \
	init_modules.o \
	invert_quda.o \
	iteration_count.o \
	mc.o \
	misc.o \
	mre.o \
	mtdagmt.o \
	m_tilde.o \
	polyakov_loop.o \
	$(RANDOM_O) \
	sc.o \
	service.o \
	staple.o \
	su3.o \
	swap.o \
	timing.o \
	traces.o \
	uuu.o \
	w_mult.o \
	xyzt2i.o

LIBS = d/$(LIBD) comm/$(LIBCOMM) clover/$(LIBCLOVER)

bqcd: bqcd.o $(MODULES) $(OBJS) $(LIBS)
	$(F90) -o $@ bqcd.o $(MODULES) $(OBJS) $(LIBS) $(SYSLIBS) $(LDFLAGS) 

fast:
	cd modules && $(MAKE)
	cd d && $(MAKE) fast 
	cd comm && $(MAKE) fast
	cd clover && $(MAKE) fast
	$(FAST_MAKE) bqcd

clean:
	rm -f bqcd.[0-9][0-9][0-9].* diag.[0-9][0-9] core app.rif
	rm -f random_test random_test.dump random_test.out
	rm -f test_echo
	rm -f a.out out out1 out2

tidy: clean
	rm -f *.[Toid] *.f90 *.mod work.pc work.pcl

clobber: tidy
	rm -f bqcd
	$(MAKE) clobber_libs
	cd modules && $(MAKE) clean

Modules:
	cd modules && $(MAKE)

libd:
	cd d && $(MAKE)

libclover: $(MODULES)
	cd clover && $(MAKE)

libs:
	cd d && $(MAKE)
	cd comm && $(MAKE)
	cd clover && $(MAKE)

clean_libs:
	cd d && $(MAKE) clean
	cd comm && $(MAKE) clean
	cd clover && $(MAKE) clean

clobber_libs:
	cd d && $(MAKE) clobber
	cd comm && $(MAKE) clobber
	cd clover && $(MAKE) clobber

the_ranf_test: ranf.o
	$(fpp) $(MYFLAGS) ranf_test.F90 ranf_test.f90
	$(F90) ranf_test.f90 ranf.o
	./a.out | diff - ranf_test.reference

test_echo: test_echo.o service.o
	$(F90) -o $@ $(LDFLAGS) test_echo.o service.o

prep:
	rm -f Makefile.var
	ln -s platform/Makefile-$(PLATFORM).var Makefile.var

prep-gnu-mpich:
	$(MAKE) prep PLATFORM=gnu-mpich

prep-intel-mpt:
	$(MAKE) prep PLATFORM=intel-mpt

prep-ibm:
	$(MAKE) prep PLATFORM=ibm

prep-ibm-bgp:
	$(MAKE) prep PLATFORM=ibm-bgp

prep-pgi:
	$(MAKE) prep PLATFORM=pgi


commit:
	cd modules; \
	count=`awk '/!count=/ { print $$2 + 1 }' module_svn.F90`; \
	ex - '+1,$$s/!count=.*/!count= '$$count'/' +wq module_svn.F90
	svn commit .

#===============================================================================
