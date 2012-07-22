#===============================================================================
#
# d/Makefile
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

include ../Makefile.var

ifdef FPP2
	fpp = $(FPP2) -I.. $(MYFLAGS)
else
	fpp = $(FPP) -I.. $(MYFLAGS)
endif


.SUFFIXES:
.SUFFIXES: .a .o .f90 .F90

.f90.o:
	$(F90) -c $(FFLAGS) $<

MODULES_DIR = ../modules

OBJS_D1 = \
	d.o  \
	d_t.o  \
	d_zb.o \
	d_zf.o \
	d_yb.o \
	d_yf.o \
	d_xb.o \
	d_xf.o \
	d_dag.o  \
	d_dag_t.o  \
	d_dag_zb.o \
	d_dag_zf.o \
	d_dag_yb.o \
	d_dag_yf.o \
	d_dag_xb.o \
	d_dag_xf.o \
	d_version.o

OBJS_D2 = \
	d2.o  \
	d2_t.o  \
	d2_zf.o \
	d2_yf.o \
	d2_xf.o \
	d2_dag.o  \
	d2_dag_t.o  \
	d2_dag_zf.o \
	d2_dag_yf.o \
	d2_dag_xf.o \
	d2_version.o

OBJS_D21 = \
	d21.o  \
	d21_t.o  \
	d21_zf.o \
	d21_yf.o \
	d21_xf.o \
	d21_dag.o  \
	d21_dag_t.o  \
	d21_dag_zf.o \
	d21_dag_yf.o \
	d21_dag_xf.o \
	d21_version.o

OBJS_D3 = \
	d3.o  \
	d3_t.o  \
	d3_zf.o \
	d3_yf.o \
	d3_xf.o \
	d3_dag.o  \
	d3_dag_t.o  \
	d3_dag_zf.o \
	d3_dag_yf.o \
	d3_dag_xf.o \
	d3_version.o

OBJS_DSF = \
	dsf_xf.o \
	dsf_xb.o \
	dsf_yf.o \
	dsf_yb.o \
	dsf_zf.o \
	dsf_zb.o \
	dsf_tf.o \
	dsf_tb.o


#-------------------------------------------------------------------------------
$(LIBD):

libd1.a: $(OBJS_D1) $(OBJS_DSF)
	$(AR) $(ARFLAGS) $@ $(OBJS_D1) $(OBJS_DSF)
	$(RANLIB) $@

libd2.a: $(OBJS_D2) $(OBJS_DSF)
	$(AR) $(ARFLAGS) $@ $(OBJS_D2) $(OBJS_DSF)
	$(RANLIB) $@

libd21.a: $(OBJS_D21) $(OBJS_DSF)
	$(AR) $(ARFLAGS) $@ $(OBJS_D21) $(OBJS_DSF)
	$(RANLIB) $@

libd3.a: $(OBJS_D3) $(OBJS_DSF)
	$(AR) $(ARFLAGS) $@ $(OBJS_D3) $(OBJS_DSF)
	$(RANLIB) $@

fast:
	$(FAST_MAKE)

#-------------------------------------------------------------------------------
d.f90: D.F90 $(DEPENDENCIES_D)
	$(fpp) -DNAME=d D.F90 > $@

d_dag.f90: D.F90 $(DEPENDENCIES_D)
	$(fpp) -DNAME=d_dag D.F90 > $@

d_xf.f90: Dxyz.F90
	$(fpp) -DNAME=d_xf -DDIR_X -DFORWARD -UDAGGER Dxyz.F90 > $@

d_yf.f90: Dxyz.F90
	$(fpp) -DNAME=d_yf -DDIR_Y -DFORWARD -UDAGGER Dxyz.F90 > $@

d_zf.f90: Dxyz.F90
	$(fpp) -DNAME=d_zf -DDIR_Z -DFORWARD -UDAGGER Dxyz.F90 > $@

d_dag_xf.f90: Dxyz.F90
	$(fpp) -DNAME=d_dag_xf -DDIR_X -DFORWARD -DDAGGER Dxyz.F90 > $@

d_dag_yf.f90: Dxyz.F90
	$(fpp) -DNAME=d_dag_yf -DDIR_Y -DFORWARD -DDAGGER Dxyz.F90 > $@

d_dag_zf.f90: Dxyz.F90
	$(fpp) -DNAME=d_dag_zf -DDIR_Z -DFORWARD -DDAGGER Dxyz.F90 > $@

d_xb.f90: Dxyz.F90
	$(fpp) -DNAME=d_xb -DDIR_X -UFORWARD -UDAGGER Dxyz.F90 > $@

d_yb.f90: Dxyz.F90
	$(fpp) -DNAME=d_yb -DDIR_Y -UFORWARD -UDAGGER Dxyz.F90 > $@

d_zb.f90: Dxyz.F90
	$(fpp) -DNAME=d_zb -DDIR_Z -UFORWARD -UDAGGER Dxyz.F90 > $@

d_dag_xb.f90: Dxyz.F90
	$(fpp) -DNAME=d_dag_xb -DDIR_X -UFORWARD -DDAGGER Dxyz.F90 > $@

d_dag_yb.f90: Dxyz.F90
	$(fpp) -DNAME=d_dag_yb -DDIR_Y -UFORWARD -DDAGGER Dxyz.F90 > $@

d_dag_zb.f90: Dxyz.F90
	$(fpp) -DNAME=d_dag_zb -DDIR_Z -UFORWARD -DDAGGER Dxyz.F90 > $@

d_t.f90: Dt.F90
	$(fpp) -DNAME=d_t -UDAGGER Dt.F90 > $@

d_dag_t.f90: Dt.F90
	$(fpp) -DNAME=d_dag_t -DDAGGER Dt.F90 > $@

d_version.f90: DVersion.F90
	$(fpp) -DVERSION=1 DVersion.F90 > $@


#-------------------------------------------------------------------------------
d2.f90:     D2.F90 $(DEPENDENCIES_D); $(fpp) -DNAME=d     D2.F90 > $@

d2_dag.f90: D2.F90 $(DEPENDENCIES_D); $(fpp) -DNAME=d_dag D2.F90 > $@

d2_xf.f90:     D2xyzt.F90; $(fpp) -DNAME=d_xf     -DDIR_X -UDAGGER D2xyzt.F90 > $@

d2_yf.f90:     D2xyzt.F90; $(fpp) -DNAME=d_yf     -DDIR_Y -UDAGGER D2xyzt.F90 > $@

d2_zf.f90:     D2xyzt.F90; $(fpp) -DNAME=d_zf     -DDIR_Z -UDAGGER D2xyzt.F90 > $@

d2_dag_xf.f90: D2xyzt.F90; $(fpp) -DNAME=d_dag_xf -DDIR_X -DDAGGER D2xyzt.F90 > $@

d2_dag_yf.f90: D2xyzt.F90; $(fpp) -DNAME=d_dag_yf -DDIR_Y -DDAGGER D2xyzt.F90 > $@

d2_dag_zf.f90: D2xyzt.F90; $(fpp) -DNAME=d_dag_zf -DDIR_Z -DDAGGER D2xyzt.F90 > $@

d2_t.f90:      D2xyzt.F90; $(fpp) -DNAME=d_t      -DDIR_T -UDAGGER D2xyzt.F90 > $@

d2_dag_t.f90:  D2xyzt.F90; $(fpp) -DNAME=d_dag_t  -DDIR_T -DDAGGER D2xyzt.F90 > $@

d2_version.f90: DVersion.F90; $(fpp) -DVERSION=2 DVersion.F90 > $@


#-------------------------------------------------------------------------------
d21.f90:     D21.F90 $(DEPENDENCIES_D); $(fpp) -DNAME=d     -UDAGGER D21.F90 > $@

d21_dag.f90: D21.F90 $(DEPENDENCIES_D); $(fpp) -DNAME=d_dag -DDAGGER D21.F90 > $@

d21_xf.f90:     D21xyzt.F90; $(fpp) -DNAME=d_xf     -DDIR_X -UDAGGER D21xyzt.F90 > $@

d21_yf.f90:     D21xyzt.F90; $(fpp) -DNAME=d_yf     -DDIR_Y -UDAGGER D21xyzt.F90 > $@

d21_zf.f90:     D21xyzt.F90; $(fpp) -DNAME=d_zf     -DDIR_Z -UDAGGER D21xyzt.F90 > $@

d21_dag_xf.f90: D21xyzt.F90; $(fpp) -DNAME=d_dag_xf -DDIR_X -DDAGGER D21xyzt.F90 > $@

d21_dag_yf.f90: D21xyzt.F90; $(fpp) -DNAME=d_dag_yf -DDIR_Y -DDAGGER D21xyzt.F90 > $@

d21_dag_zf.f90: D21xyzt.F90; $(fpp) -DNAME=d_dag_zf -DDIR_Z -DDAGGER D21xyzt.F90 > $@

d21_t.f90:      D21xyzt.F90; $(fpp) -DNAME=d_t      -DDIR_T -UDAGGER D21xyzt.F90 > $@

d21_dag_t.f90:  D21xyzt.F90; $(fpp) -DNAME=d_dag_t  -DDIR_T -DDAGGER D21xyzt.F90 > $@

d21_version.f90: DVersion.F90; $(fpp) -DVERSION=21 DVersion.F90 > $@


#-------------------------------------------------------------------------------
d3.f90:     D3.F90 $(DEPENDENCIES_D); $(fpp) -DNAME=d     D3.F90 > $@

d3_dag.f90: D3.F90 $(DEPENDENCIES_D); $(fpp) -DNAME=d_dag D3.F90 > $@

d3_xf.f90:     D3xyzt.F90; $(fpp) -DNAME=d_xf     -DDIR_X -UDAGGER D3xyzt.F90 > $@

d3_yf.f90:     D3xyzt.F90; $(fpp) -DNAME=d_yf     -DDIR_Y -UDAGGER D3xyzt.F90 > $@

d3_zf.f90:     D3xyzt.F90; $(fpp) -DNAME=d_zf     -DDIR_Z -UDAGGER D3xyzt.F90 > $@

d3_dag_xf.f90: D3xyzt.F90; $(fpp) -DNAME=d_dag_xf -DDIR_X -DDAGGER D3xyzt.F90 > $@

d3_dag_yf.f90: D3xyzt.F90; $(fpp) -DNAME=d_dag_yf -DDIR_Y -DDAGGER D3xyzt.F90 > $@

d3_dag_zf.f90: D3xyzt.F90; $(fpp) -DNAME=d_dag_zf -DDIR_Z -DDAGGER D3xyzt.F90 > $@

d3_t.f90:      D3xyzt.F90; $(fpp) -DNAME=d_t      -DDIR_T -UDAGGER D3xyzt.F90 > $@

d3_dag_t.f90:  D3xyzt.F90; $(fpp) -DNAME=d_dag_t  -DDIR_T -DDAGGER D3xyzt.F90 > $@

d3_version.f90: DVersion.F90; $(fpp) -DVERSION=3 DVersion.F90 > $@


#-------------------------------------------------------------------------------
dsf.f90:  dsf.F90 $(DEPENDENCIES_DSF);  $(fpp) dsf.F90 > $@

dsf_xf.f90:  DSFxyzt.F90;  $(fpp) -DNAME=dsf_xf -DDIR_X -DFORWARD DSFxyzt.F90 > $@

dsf_yf.f90:  DSFxyzt.F90;  $(fpp) -DNAME=dsf_yf -DDIR_Y -DFORWARD DSFxyzt.F90 > $@

dsf_zf.f90:  DSFxyzt.F90;  $(fpp) -DNAME=dsf_zf -DDIR_Z -DFORWARD DSFxyzt.F90 > $@

dsf_tf.f90:  DSFxyzt.F90;  $(fpp) -DNAME=dsf_tf -DDIR_T -DFORWARD DSFxyzt.F90 > $@

dsf_xb.f90:  DSFxyzt.F90;  $(fpp) -DNAME=dsf_xb -DDIR_X -UFORWARD DSFxyzt.F90 > $@

dsf_yb.f90:  DSFxyzt.F90;  $(fpp) -DNAME=dsf_yb -DDIR_Y -UFORWARD DSFxyzt.F90 > $@

dsf_zb.f90:  DSFxyzt.F90;  $(fpp) -DNAME=dsf_zb -DDIR_Z -UFORWARD DSFxyzt.F90 > $@

dsf_tb.f90:  DSFxyzt.F90;  $(fpp) -DNAME=dsf_tb -DDIR_T -UFORWARD DSFxyzt.F90 > $@


#-------------------------------------------------------------------------------
clean:
	rm -f *.[Tiod] *.f90 *.mod work.pc work.pcl

clobber: clean
	rm -f libd.a libd2.a libd21.a libd3.a
