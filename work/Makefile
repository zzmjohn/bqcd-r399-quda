#===============================================================================
#
# work/Makefile
#
#-------------------------------------------------------------------------------
#
# Copyright (C) 2006-2012 Hinnerk Stueben
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

nothing:

test:
	rm -f out.TMP out.tmp
	../bqcd input.test out.TMP
	grep ' %[fim][atc]' out.TMP > out.tmp
	grep ' %[fim][atc]' out.test.reference | diff - out.tmp

quick-test:
	rm -f out.TMP out.tmp
	../bqcd input.quick-test out.TMP
	grep ' %[fim][atc]' out.TMP > out.tmp
	grep ' %[fim][atc]' out.quick-test.reference | diff - out.tmp

clean:
	rm -f bqcd.*.u
	rm -f bqcd.*.info
	rm -f bqcd.*.count
	rm -f bqcd.*.ran
	rm -f out.TMP out.tmp

#===============================================================================
