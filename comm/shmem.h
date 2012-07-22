#/*
#===============================================================================
#
# shmem.h
#
#-------------------------------------------------------------------------------
#
# Copyright (C) 2006 Hinnerk Stueben
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

#ifdef ALTIX
# define barrier shmem_barrier_all
# define shmem_broadcast shmem_broadcast8
# define shmem_get shmem_get8
# define shmem_put shmem_put8
# define shpalloc(addr, length, errcode, abort) shpalloc(addr, 2 * (length), errcode, abort)
#endif
