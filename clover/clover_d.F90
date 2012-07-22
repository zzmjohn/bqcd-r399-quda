!===============================================================================
!
! clover_d.F90 - derivative of clover term
!
!-------------------------------------------------------------------------------
!
! Copyright (C) 1998-2005 Hinnerk Stueben
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
!
!    E -- 2 -- B
!    |         |                  A = x
!    3         1     ^ nu         B = x + mu^ + nu^
!    |         |     |            C = x + mu^ - nu^
!    A -- 0 -- D     x --> mu     
!    |         |                  D = x + mu^
!    4         6                  E = x + nu^
!    |         |                  F = x - nu^
!    F -- 5 -- C
!    
!-------------------------------------------------------------------------------
# include "defs.h"

!-------------------------------------------------------------------------------
subroutine clover_dsd(eo, p, b, s, u)

  use typedef_clover
  use module_p_interface
  use module_vol
  implicit none

  integer         :: eo    ! EVEN/ODD property of "b"
  GENERATOR_FIELD :: p
  CLOVER_FIELD_B  :: b
  REAL            :: s
  GAUGE_FIELD     :: u

  P_GAUGE_FIELD, save :: w     ! use existing data structure

  CLOVER_FIELD_C  :: t
!dir$ cache_align    t


  TIMING_START(timing_bin_clover_dsd)

  ALLOCATE_G_FIELD(w)

  call clover_t_init(t, b(1, 1, eo))

  call clover_ts(1, 2, w(1, 1, 1, EVEN, 1), t) ; call xbound_g(w, EVEN, 1)
  call clover_ts(1, 3, w(1, 1, 1, ODD,  1), t) ; call xbound_g(w, ODD,  1)
  call clover_ts(1, 4, w(1, 1, 1, EVEN, 2), t) ; call xbound_g(w, EVEN, 2)
  call clover_ts(2, 3, w(1, 1, 1, ODD,  2), t) ; call xbound_g(w, ODD,  2)
  call clover_ts(2, 4, w(1, 1, 1, EVEN, 3), t) ; call xbound_g(w, EVEN, 3)
  call clover_ts(3, 4, w(1, 1, 1, ODD,  3), t) ; call xbound_g(w, ODD,  3)

  call clover_d(eo, p, s, u, w)

  TIMING_STOP(timing_bin_clover_dsd)
end

!-------------------------------------------------------------------------------
subroutine clover_dsf(eo, p, b, a, s, u)

  use module_p_interface
  use module_vol
  implicit none

  integer         :: eo    ! EVEN/ODD property of "b" and "a"
  GENERATOR_FIELD :: p
  SPINCOL_FIELD   :: b, a
  REAL            :: s
  GAUGE_FIELD     :: u

  P_GAUGE_FIELD, save :: w     ! use existing data structure


  TIMING_START(timing_bin_clover_dsf)

  ALLOCATE_G_FIELD(w)

  call clover_bsa(1, 2, w(1, 1, 1, EVEN, 1), b, a) ; call xbound_g(w, EVEN, 1)
  call clover_bsa(1, 3, w(1, 1, 1, ODD,  1), b, a) ; call xbound_g(w, ODD,  1)
  call clover_bsa(1, 4, w(1, 1, 1, EVEN, 2), b, a) ; call xbound_g(w, EVEN, 2)
  call clover_bsa(2, 3, w(1, 1, 1, ODD,  2), b, a) ; call xbound_g(w, ODD,  2)
  call clover_bsa(2, 4, w(1, 1, 1, EVEN, 3), b, a) ; call xbound_g(w, EVEN, 3)
  call clover_bsa(3, 4, w(1, 1, 1, ODD,  3), b, a) ; call xbound_g(w, ODD,  3)

  call clover_d(eo, p, s, u, w)

  TIMING_STOP(timing_bin_clover_dsf)
end

!-------------------------------------------------------------------------------
subroutine clover_d(eo, p, s, u, w)

  use module_vol
  implicit none

  integer         :: eo
  GENERATOR_FIELD :: p
  REAL            :: s
  GAUGE_FIELD     :: u, w

  call clover_d_mu_nu(eo, 1, 2, p, s, u, w(1, 1, 1, EVEN, 1))
  call clover_d_mu_nu(eo, 1, 3, p, s, u, w(1, 1, 1, ODD,  1))
  call clover_d_mu_nu(eo, 1, 4, p, s, u, w(1, 1, 1, EVEN, 2))
  call clover_d_mu_nu(eo, 2, 3, p, s, u, w(1, 1, 1, ODD,  2))
  call clover_d_mu_nu(eo, 2, 4, p, s, u, w(1, 1, 1, EVEN, 3))
  call clover_d_mu_nu(eo, 3, 4, p, s, u, w(1, 1, 1, ODD,  3))

end

!-------------------------------------------------------------------------------
subroutine clover_d_mu_nu(e, mu, nu, p, s, u, w)

  use module_vol
  implicit none

  integer         :: e, o, mu, nu
  GENERATOR_FIELD :: p
  REAL            :: s
  GAUGE_FIELD     :: u
  SU3_FIELD       :: w
  external           clover_d_same_eo, clover_d_diff_eo

  o = EVEN + ODD - e

  call clover_d_loop(e, mu, nu, p,  s, u, w, clover_d_same_eo)
  call clover_d_loop(e, nu, mu, p, -s, u, w, clover_d_same_eo)
  call clover_d_loop(o, mu, nu, p,  s, u, w, clover_d_diff_eo)
  call clover_d_loop(o, nu, mu, p, -s, u, w, clover_d_diff_eo)

end

!-------------------------------------------------------------------------------
subroutine clover_d_loop(e, mu, nu, p, s, u, w, clover_dd)

  use module_vol
  use module_nn
  implicit none

  integer         :: e, mu, nu
  GENERATOR_FIELD :: p
  REAL            :: s
  GAUGE_FIELD     :: u
  SU3_FIELD       :: w
  external           clover_dd

  integer         :: o, i, ia, ib, ic, id, ie, if, j
  GENERATOR       :: q

  o = EVEN + ODD - e

  !$omp parallel do private(ia, ib, ic, id, ie, if, j, q)
  do i = 1, volh

     id = nn(i, e, mu, FWD)
     ie = nn(i, e, nu, FWD)
     if = nn(i, e, nu, BWD)
     
     ia = i
     ib = nn(id, o, nu, FWD)
     ic = nn(id, o, nu, BWD)
     
     call clover_dd(q, s, u(1, 1, ia, e, mu), &
                          u(1, 1, id, o, nu), &
                          u(1, 1, ie, o, mu), &
                          u(1, 1, ia, e, nu), &
                          u(1, 1, if, o, nu), &
                          u(1, 1, if, o, mu), &
                          u(1, 1, ic, e, nu), &
                          w(1, 1, ia), &
                          w(1, 1, ib), &
                          w(1, 1, ic), &
                          w(1, 1, id), &
                          w(1, 1, ie), &
                          w(1, 1, if))

     do j = 1, NGEN
        p(j, i, e, mu) = p(j, i, e, mu) + q(j)
     enddo
  enddo

end

!-------------------------------------------------------------------------------
subroutine clover_d_same_eo(p, s, u0,u1,u2,u3,u4,u5,u6, wa,wb,wc, wd,we,wf)

  implicit none

  GENERATOR :: p
  REAL      :: s
  SU3       :: u0, u1, u2, u3, u4, u5, u6, wa, wb, wc, wd, we, wf
  SU3       :: r, u, v

  p = ZERO

  call uu(u, u0, u1)
  call uu(v, u3, u2)
  
  call clover_uuu_udu(r, u, v, wa) ; call re_tr_j(p, r, s)
  call clover_uuu_uud(r, wa, v, u) ; call re_tr_j(p, r, s)
  call clover_uuu_uud(r, u, wb, v) ; call re_tr_j(p, r, s)
  call clover_uuu_uud(r, v, wb, u) ; call re_tr_j(p, r, s)

  call uud(u, u0, u6)
  call udu(v, u4, u5)

  call clover_uuu_uud(r, wa, v, u) ; call re_tr_j(p, r, -s)
  call clover_uuu_udu(r, u, v, wa) ; call re_tr_j(p, r, -s)
  call clover_uuu_uud(r, v, wc, u) ; call re_tr_j(p, r, -s)
  call clover_uuu_uud(r, u, wc, v) ; call re_tr_j(p, r, -s)

end

!-------------------------------------------------------------------------------
subroutine clover_d_diff_eo(p, s, u0,u1,u2,u3,u4,u5,u6, wa,wb,wc, wd,we,wf)

  implicit none

  GENERATOR :: p
  REAL      :: s
  SU3       :: u0, u1, u2, u3, u4, u5, u6, wa, wb, wc, wd, we, wf
  SU3       :: r, u, v

  p = ZERO

  u = ZERO
  v = ZERO
  call uuu_fwd(u, u1, u2, u3)
  call uuu_fwd(v, u2, u1, u0)

  call clover_uuu_uuu(r, u0, wd, u) ; call re_tr_j(p, r, s)
  call clover_uuu_dud(r, u, wd, u0) ; call re_tr_j(p, r, s)
  call clover_uuu_dud(r, v, we, u3) ; call re_tr_j(p, r, s)
  call clover_uuu_uuu(r, u3, we, v) ; call re_tr_j(p, r, s)

  u = ZERO
  v = ZERO
  call uuu_bwd(u, u6, u5, u4)
  call uuu_fwd(v, u0, u6, u5)

  call clover_uuu_dud(r, u, wd, u0) ; call re_tr_j(p, r, -s)
  call clover_uuu_uuu(r, u0, wd, u) ; call re_tr_j(p, r, -s)
  call clover_uuu_dud(r, u4, wf, v) ; call re_tr_j(p, r, -s)
  call clover_uuu_uuu(r, v, wf, u4) ; call re_tr_j(p, r, -s)

end

!===============================================================================
