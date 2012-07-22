#/*
#===============================================================================
#
# clover.h
#
#-------------------------------------------------------------------------------
#
# Copyright (C) 1998-2001 Hinnerk Stueben
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

#ifdef CLOVER_AS_COMPLEX_ARRAY

# define A11 Re(a(1,J,i))
# define A22 Im(a(1,J,i))
# define A33 Re(a(11,J,i))
# define A44 Im(a(11,J,i))
# define A55 Re(a(17,J,i))
# define A66 Im(a(17,J,i))

# define A12 a(2,J,i)
# define A13 a(3,J,i)
# define A14 a(4,J,i)
# define A15 a(5,J,i)
# define A16 a(6,J,i)

# define A23 a(7,J,i)
# define A24 a(8,J,i)
# define A25 a(9,J,i)
# define A26 a(10,J,i)

# define A34 a(12,J,i)
# define A35 a(13,J,i)
# define A36 a(14,J,i)

# define A45 a(15,J,i)
# define A46 a(16,J,i)

# define A56 a(18,J,i)

# define B11 Re(b(16,J,i))
# define B22 Im(b(16,J,i))
# define B33 Re(b(17,J,i))
# define B44 Im(b(17,J,i))
# define B55 Re(b(18,J,i))
# define B66 Im(b(18,J,i))

# define B21 b(1,J,i)

# define B31 b(2,J,i)
# define B32 b(3,J,i)

# define B41 b(4,J,i)
# define B42 b(5,J,i)
# define B43 b(6,J,i)

# define B51 b(7,J,i)
# define B52 b(8,J,i)
# define B53 b(9,J,i)
# define B54 b(10,J,i)

# define B61 b(11,J,i)
# define B62 b(12,J,i)
# define B63 b(13,J,i)
# define B64 b(14,J,i)
# define B65 b(15,J,i)

#else

# define A11 a%i11
# define A22 a%i22
# define A33 a%i33
# define A44 a%i44
# define A55 a%i55
# define A66 a%i66

# define A12 a%i12
# define A13 a%i13
# define A14 a%i14
# define A15 a%i15
# define A16 a%i16

# define A23 a%i23
# define A24 a%i24
# define A25 a%i25
# define A26 a%i26

# define A34 a%i34
# define A35 a%i35
# define A36 a%i36

# define A45 a%i45
# define A46 a%i46

# define A56 a%i56

# define B11 b%i11
# define B22 b%i22
# define B33 b%i33
# define B44 b%i44
# define B55 b%i55
# define B66 b%i66

# define B21 b%i21

# define B31 b%i31
# define B32 b%i32

# define B41 b%i41
# define B42 b%i42
# define B43 b%i43

# define B51 b%i51
# define B52 b%i52
# define B53 b%i53
# define B54 b%i54

# define B61 b%i61
# define B62 b%i62
# define B63 b%i63
# define B64 b%i64
# define B65 b%i65

#endif

# define A21 conjg(A12)
# define A31 conjg(A13)
# define A41 conjg(A14)
# define A51 conjg(A15)
# define A61 conjg(A16)

# define A32 conjg(A23)
# define A42 conjg(A24)
# define A52 conjg(A25)
# define A62 conjg(A26)

# define A43 conjg(A34)
# define A53 conjg(A35)
# define A63 conjg(A36)

# define A54 conjg(A45)
# define A64 conjg(A46)

# define A65 conjg(A56)

# define B12 conjg(B21)

# define B13 conjg(B31)
# define B23 conjg(B32)

# define B14 conjg(B41)
# define B24 conjg(B42)
# define B34 conjg(B43)

# define B15 conjg(B51)
# define B25 conjg(B52)
# define B35 conjg(B53)
# define B45 conjg(B54)

# define B16 conjg(B61)
# define B26 conjg(B62)
# define B36 conjg(B63)
# define B46 conjg(B64)
# define B56 conjg(B65)

# define SC1 1, 1
# define SC2 1, 2
# define SC3 1, 3
# define SC4 2, 1
# define SC5 2, 2
# define SC6 2, 3
# define SC7 3, 1
# define SC8 3, 2
# define SC9 3, 3
# define SC10 4, 1
# define SC11 4, 2
# define SC12 4, 3
