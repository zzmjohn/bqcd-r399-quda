!===============================================================================
!
! su3.F90 - SU(3) routines
!
!-------------------------------------------------------------------------------
!
! Copyright (C) 1998-2002 Hinnerk Stueben
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
subroutine gen2u(u, h)  ! u := exp(i lambda_j h_j)
 
! adapted from:
 
!  Program qcdf90, module generator_algebra, version 4.0.0
 
!  Copyright by Indranil Dasgupta, Andrea R. Levi, Vittorio Lubicz
!  and Claudio Rebbi  -  Boston University  -  January 1996
!  This program may be freely copied and used as long as this notice
!  is retained.
 
  implicit none
 
  GENERATOR, intent(in) :: h
  SU3, intent(out) :: u
 
  REAL :: p, q, a, alpha, l1, l2, l3, l12, s, aux, c, d, cs1, cs2
  REAL :: a8, a12, a45, a67
  GENERATOR ::  h2, hs, hk, hs2
  COMPLEX :: ck1, ck2, ck3, ck4
  SU3 :: ms, mk      
  integer :: i, j, k
 
  integer, parameter :: rkind = RKIND
  COMPLEX, parameter :: iu = (ZERO, ONE)
  REAL, parameter :: eps = 0.00000001_rkind
 
  REAL, parameter :: sqrt33 = SQRT3 / THREE
  REAL, parameter :: twosqrt33 = TWO * sqrt33 

!  h2=.Sq.h, inlined:
     
         a8 = h(8)*SQRT33
         a12 = h(1)**2+h(2)**2
         a45 = h(4)**2+h(5)**2
         a67 = h(6)**2+h(7)**2
         h2(1) = 2*h(1)*a8+h(4)*h(6)+h(5)*h(7)
         h2(2) = 2*h(2)*a8+h(5)*h(6)-h(4)*h(7)
         h2(3) = 2*h(3)*a8+HALF*(a45-a67)
         h2(4) = h(4)*(h(3)-a8)+h(1)*h(6)-h(2)*h(7)
         h2(5) = h(5)*(h(3)-a8)+h(1)*h(7)+h(2)*h(6)
         h2(6) = h(6)*(-h(3)-a8)+h(1)*h(4)+h(2)*h(5)
         h2(7) = h(7)*(-h(3)-a8)+h(1)*h(5)-h(2)*h(4)
         h2(8) = (h(3)**2-h(8)**2+a12-HALF*(a45+a67))*SQRT33
     
!  q = .Tr.h, inlined:
     
         q = h(1)**2
         DO i = 2,8
            q = q+h(i)**2
         END DO
         q = TWO*q
     
!  p = (h*h2)/THREE, inlined:
     
         p = h(1)*h2(1)
         DO i = 2,8
            p = p+h(i)*h2(i)
         END DO
         p = TWO*p/THREE
     
         a = SQRT(TWO*q/THREE)
         alpha = ACOS(FOUR*p/a**3)/THREE      
         IF(alpha <= PI/6) THEN
            l1 = a*COS(alpha)
            l2 = a*COS(alpha+2*PI/3)
         ELSE
            l2 = a*COS(alpha)
            l1 = a*COS(alpha+2*PI/3)
         ENDIF
         l3 = -l1-l2  
     
         l12 = l1*l2
         s = -l1-2*l2
     
         aux = (TWO*l3*l3+l12)*(l1-l2)
         c = s*(l3*l3+TWO*l12)/aux  
         d = -THREE*s*l3/aux  
     
!  hs = c*h+d*h2, and
!  hk = h-hs, inlined:
     
         DO i = 1,8 
            hs(i) = c*h(i)+d*h2(i)
            hk(i) = h(i)-hs(i)
         END DO
     
!  hs2 = .Sq.hs, inlined:
     
         a8 = hs(8)*SQRT33
         a12 = hs(1)**2+hs(2)**2
         a45 = hs(4)**2+hs(5)**2
         a67 = hs(6)**2+hs(7)**2
         hs2(1) = 2*hs(1)*a8+hs(4)*hs(6)+hs(5)*hs(7)
         hs2(2) = 2*hs(2)*a8+hs(5)*hs(6)-hs(4)*hs(7)
         hs2(3) = 2*hs(3)*a8+HALF*(a45-a67)
         hs2(4) = hs(4)*(hs(3)-a8)+hs(1)*hs(6)-hs(2)*hs(7)
         hs2(5) = hs(5)*(hs(3)-a8)+hs(1)*hs(7)+hs(2)*hs(6)
         hs2(6) = hs(6)*(-hs(3)-a8)+hs(1)*hs(4)+hs(2)*hs(5)
         hs2(7) = hs(7)*(-hs(3)-a8)+hs(1)*hs(5)-hs(2)*hs(4)
         hs2(8) = (hs(3)**2-hs(8)**2+a12-HALF*(a45+a67))*SQRT33
     
         IF(ABS(s) > eps) THEN
            cs1 = SIN(s)/s
            cs2 = (COS(s)-1)/s**2
         ELSE
            cs1 = 1
            cs2 = -HALF
         ENDIF 
     
         ck1 = EXP(IU*l3)
         ck2 = 1/ck1**2
         ck3 = (ck2+2*ck1)/3
         IF(ABS(l3) > eps) THEN
            ck4 = (ck1-ck2)/(3*l3)
         ELSE
            ck4 = 3*IU
         ENDIF
     
!  aux = .Tr.hs, inlined:
     
         aux = hs(1)**2
         DO i = 2,8
            aux = aux+hs(i)**2
         END DO
         aux = TWO*aux
     
!  ms = UNIT+IU*cs1*(.Matrix.hs.)+cs2*(.Matrix.hs)*(.Matrix.hs), inlined:
     
         ms(1,1) = ONE+cs2*aux/THREE                                    &
            +CMPLX(cs2*(hs2(3)+SQRT33*hs2(8)),cs1*(hs(3)+SQRT33*hs(8)),RKIND)
         ms(2,2) = ONE+cs2*aux/THREE                                    &
            +CMPLX(cs2*(-hs2(3)+SQRT33*hs2(8)),cs1*(-hs(3)+SQRT33*hs(8)),RKIND)
         ms(3,3) = ONE+cs2*aux/THREE                                    &
            +CMPLX(-cs2*TWOSQRT33*hs2(8),-cs1*TWOSQRT33*hs(8),RKIND)
         ms(1,2) = CMPLX(cs2*hs2(1)+cs1*hs(2),cs1*hs(1)-cs2*hs2(2),RKIND)
         ms(2,1) = CMPLX(cs2*hs2(1)-cs1*hs(2),cs1*hs(1)+cs2*hs2(2),RKIND)
         ms(1,3) = CMPLX(cs2*hs2(4)+cs1*hs(5),cs1*hs(4)-cs2*hs2(5),RKIND)
         ms(3,1) = CMPLX(cs2*hs2(4)-cs1*hs(5),cs1*hs(4)+cs2*hs2(5),RKIND)
         ms(2,3) = CMPLX(cs2*hs2(6)+cs1*hs(7),cs1*hs(6)-cs2*hs2(7),RKIND)
         ms(3,2) = CMPLX(cs2*hs2(6)-cs1*hs(7),cs1*hs(6)+cs2*hs2(7),RKIND)
     
!  mk = ck3*UNIT+ck4*(.Matrix.hk), inlined:
     
         mk(1,1) = ck3+ck4*(hk(3)+SQRT33*hk(8))
         mk(2,2) = ck3+ck4*(-hk(3)+SQRT33*hk(8))
         mk(3,3) = ck3-ck4*TWOSQRT33*hk(8)
         mk(1,2) = ck4*CMPLX(hk(1),-hk(2),RKIND)
         mk(2,1) = ck4*CMPLX(hk(1),hk(2),RKIND)
         mk(1,3) = ck4*CMPLX(hk(4),-hk(5),RKIND)
         mk(3,1) = ck4*CMPLX(hk(4),hk(5),RKIND)
         mk(2,3) = ck4*CMPLX(hk(6),-hk(7),RKIND)
         mk(3,2) = ck4*CMPLX(hk(6),hk(7),RKIND)
     
!  u = ms*mk, inlined:
     
         DO i = 1,3 
         DO j = 1,3 
            u(i,j) = ms(i,1)*mk(1,j)
            DO k = 2,3
               u(i,j) = u(i,j)+ms(i,k)*mk(k,j)
            END DO
         END DO
         END DO
         
END

!-------------------------------------------------------------------------------
subroutine im_tr_j(p, u, s)  ! p(j) := p(j) + s * Im Tr(lambda_j U)

  implicit none

  GENERATOR :: p
  SU3 :: u
  REAL :: s

  p(1) = p(1) + s * (Im(u(1, 2)) + Im(u(2, 1)))
  p(2) = p(2) + s * (Re(u(1, 2)) - Re(u(2, 1)))
  p(3) = p(3) + s * (Im(u(1, 1)) - Im(u(2, 2)))
  p(4) = p(4) + s * (Im(u(1, 3)) + Im(u(3, 1)))
  p(5) = p(5) + s * (Re(u(1, 3)) - Re(u(3, 1)))
  p(6) = p(6) + s * (Im(u(2, 3)) + Im(u(3, 2)))
  p(7) = p(7) + s * (Re(u(2, 3)) - Re(u(3, 2)))
  p(8) = p(8) + s * (Im(u(1, 1)) + Im(u(2, 2)) - TWO * Im(u(3, 3))) / SQRT3

end

!-------------------------------------------------------------------------------
subroutine re_tr_j(p, u, s)  ! p(j) := p(j) + s * Re Tr(lambda_j U)

  implicit none

  GENERATOR :: p
  SU3 :: u
  REAL :: s

  p(1) = p(1) + s * (Re(u(1, 2)) + Re(u(2, 1)))
  p(2) = p(2) + s * (Im(u(2, 1)) - Im(u(1, 2)))
  p(3) = p(3) + s * (Re(u(1, 1)) - Re(u(2, 2)))
  p(4) = p(4) + s * (Re(u(1, 3)) + Re(u(3, 1)))
  p(5) = p(5) + s * (Im(u(3, 1)) - Im(u(1, 3)))
  p(6) = p(6) + s * (Re(u(2, 3)) + Re(u(3, 2)))
  p(7) = p(7) + s * (Im(u(3, 2)) - Im(u(2, 3)))
  p(8) = p(8) + s * (Re(u(1, 1)) + Re(u(2, 2)) - TWO * Re(u(3, 3))) / SQRT3

end

!-------------------------------------------------------------------------------
subroutine su3_check(u)  ! checks if "u" is in SU(3)

  implicit none
  SU3 :: u, v
  SU3, parameter :: su3_one = reshape( &
                           (/ ONE,ZERO,ZERO, &
                              ZERO,ONE,ZERO, &
                              ZERO,ZERO,ONE /), &
                           (/ NCOL, NCOL /))
  REAL, parameter :: eps = 1e-13
  REAL :: dev
  integer :: i, j

  call uud(v, u, u)
  dev = ZERO
  do i = 1, NCOL
     do j = 1, NCOL
        dev = dev + abs(Re(v(i, j)) - Re(su3_one(i, j))) & 
                  + abs(Im(v(i, j)) - Im(su3_one(i, j)))
     enddo
  enddo

  if (dev > eps) call die('su3_check(): dev > eps')

  call su3_check_det(u)

end

!-------------------------------------------------------------------------------
subroutine su3_check_det(u)

  implicit none
  COMPLEX :: det
  SU3 :: u
  REAL, parameter :: eps = 1e-13

  det = u(1,1) * u(2,2) * u(3,3) &
      + u(1,2) * u(2,3) * u(3,1) &
      + u(1,3) * u(2,1) * u(3,2) &
      - u(1,1) * u(2,3) * u(3,2) &
      - u(1,2) * u(2,1) * u(3,3) &
      - u(1,3) * u(2,2) * u(3,1)

  if (abs(Re(det) - ONE) > eps) call die("check_su3_det(): Re(det) /= 1")
  if (abs(Im(det))       > eps) call die("check_su3_det(): Im(det) /= 0")

end

!-------------------------------------------------------------------------------
subroutine u_add(u, v)  ! u := u + v

  implicit none

  SU3 :: u, v
  integer i, j

  do j = 1, NCOL
     do i = 1, NCOL
        u(i, j) = u(i, j) + v(i, j)
     enddo
  enddo

end

!-------------------------------------------------------------------------------
subroutine u_complete(u)  ! calculate 3rd column form the first two

  implicit none
  SU3 :: u

  u(1,3) = conjg(u(2,1) * u(3,2) - u(3,1) * u(2,2))
  u(2,3) = conjg(u(3,1) * u(1,2) - u(1,1) * u(3,2))
  u(3,3) = conjg(u(1,1) * u(2,2) - u(2,1) * u(1,2))

end

!-------------------------------------------------------------------------------
subroutine u_normalize(u)

! from qcdsf_t3e program:
!
!  u_normalize() takes a complex matrix and produces a true su3 matrix from
!  the upper 6 entries (because of FORTRAN the first two rows -> DIFFERS
!  FROM APE-PROGRAM!!)
!
!             ( * * . )
!             ( * * . )
!             ( * * . )      (right 3 completely ignored)
!
!  Normalization done by Gramm-Schmitt

  implicit none
  SU3 :: u
  COMPLEX :: f
  REAL :: len
  integer :: i

  len = real(u(1,1))**2 + aimag(u(1,1))**2 + &              ! length u_1
        real(u(2,1))**2 + aimag(u(2,1))**2 + &
        real(u(3,1))**2 + aimag(u(3,1))**2
  len = sqrt(len)
 
  do i = 1, NCOL
     u(i,1) = u(i,1) / len                                  ! normalize u_1
  enddo
 
  f = u(1,2) * conjg(u(1,1)) + &
      u(2,2) * conjg(u(2,1)) + &
      u(3,2) * conjg(u(3,1))
 
  do i = 1, NCOL
     u(i,2) = u(i,2) - f * u(i,1)                           ! orthogonalize
  enddo
 
  len = real(u(1,2))**2 + aimag(u(1,2))**2 + &              ! length u_2
        real(u(2,2))**2 + aimag(u(2,2))**2 + &
        real(u(3,2))**2 + aimag(u(3,2))**2
  len = sqrt(len)
 
  do i = 1, NCOL
     u(i,2) = u(i,2) / len                                  ! normalize u_2
  enddo
 
  call u_complete(u)

!  u(1,3) = conjg(u(2,1) * u(3,2) - u(3,1) * u(2,2))         ! calculate u_3
!  u(2,3) = conjg(u(3,1) * u(1,2) - u(1,1) * u(3,2))         ! = u_1 x u_2
!  u(3,3) = conjg(u(1,1) * u(2,2) - u(2,1) * u(1,2))

end

!-------------------------------------------------------------------------------
subroutine u_trans(u)  ! u := transpose(u)

  implicit none

  SU3, intent(inout) :: u
  COMPLEX :: tmp
  integer :: i, j

  do j = 1, NCOL
     do i = j + 1, NCOL
        tmp = u(i, j)
        u(i, j) = u(j, i)
        u(j, i) = tmp
     enddo
  enddo

end

!-------------------------------------------------------------------------------
subroutine u_update(u, v)  ! u = v * u

  implicit none

  SU3, intent(in) :: v
  SU3, intent(inout) :: u
  SU3 :: w

  w = u
  call uu(u, v, w)

end

!-------------------------------------------------------------------------------
subroutine u_update2(u, v)  ! u = u * v

  implicit none

  SU3, intent(in) :: v
  SU3, intent(inout) :: u
  SU3 :: w

  w = u
  call uu(u, w, v)

end

!-------------------------------------------------------------------------------
subroutine uu(r, a, b)  ! r = a * b

  implicit none

  SU3 :: r, a, b
  integer :: i, j

  do i = 1, NCOL
     do j = 1, NCOL
        r(i, j) = a(i, 1) * b(1, j) &
                + a(i, 2) * b(2, j) &
                + a(i, 3) * b(3, j)
     enddo
  enddo

end

!-------------------------------------------------------------------------------
subroutine uud(r, a, b)  ! U U^dagger:  r = a * b+

  implicit none

  SU3 :: r, a, b
  integer :: i, j

  do i = 1, NCOL
     do j = 1, NCOL
        r(i, j) = a(i, 1) * conjg(b(j, 1)) &
                + a(i, 2) * conjg(b(j, 2)) &
                + a(i, 3) * conjg(b(j, 3))
     enddo
  enddo

end

!-------------------------------------------------------------------------------
subroutine udu(r, a, b)  ! U^dagger U:  r = a+ * b

  implicit none

  SU3 :: r, a, b
  integer :: i, j

  do i = 1, NCOL
     do j = 1, NCOL
        r(i, j) = conjg(a(1, i)) * b(1, j) &
                + conjg(a(2, i)) * b(2, j) &
                + conjg(a(3, i)) * b(3, j)
     enddo
  enddo

end

!-------------------------------------------------------------------------------
REAL function Re_Tr_uu(u, v)  ! returns Re(Tr(u * v))

  implicit none
  SU3, intent(in) :: u, v
  REAL            :: p
  integer         :: c1, c2
  
  p = 0
  do c2 = 1, NCOL
     do c1 = 1, NCOL
        p = p + Re(u(c2, c1)) * Re(v(c1, c2)) &
              - Im(u(c2, c1)) * Im(v(c1, c2))
     enddo
  enddo
  
  Re_Tr_uu = p

end

!===============================================================================
