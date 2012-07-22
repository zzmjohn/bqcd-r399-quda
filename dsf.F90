!===============================================================================
!
! dsf.F90  -  kernel of: p(j,x,mu) := p(j,x,mu) - step * D_{x,mu,j} S_{f 1|2}
!
!-------------------------------------------------------------------------------
!
! Copyright (C) 2003-2006 Hinnerk Stueben
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
subroutine dsf(p, conf, step, para, a, b)

  use typedef_hmc
  use module_nn
  use module_p_interface
  use module_vol
  implicit none

  type(hmc_para),  intent(in)    :: para
  type(hmc_conf),  intent(in)    :: conf
  GENERATOR_FIELD, intent(inout) :: p
  REAL,            intent(in)    :: step
  SPINCOL_FIELD,   intent(in)    :: a
  SPINCOL_FIELD,   intent(in)    :: b

  P_GAUGE_FIELD                  :: u
  P_SPINCOL_FIELD, save          :: at, bt
  REAL                           :: s, s1, s2


  !! call flip_bc(u)  <- done in calling routine

  ALLOCATE_SC_FIELD(at)
  ALLOCATE_SC_FIELD(bt)

  u => conf%u

  if (para%kappa /= ZERO) then
     call d(ODD, EVEN, at, a, u)                         ! A~ = Doe A
     call d_dag(ODD, EVEN, bt, b, u)                     ! B~ = Deo+ B

     if (para%csw_kappa /= ZERO) then
        call clover_mult_b(conf%b(1,1,ODD), at, volh)    ! A~ = inv(Too) A~
        call clover_mult_b(conf%b(1,1,ODD), bt, volh)    ! B~ = inv(Too) B~
     endif

     if (para%h /= ZERO) then
        call h_mult_b(-para%h, at, volh)                 ! A~ ~ inv(H) A~
        call h_mult_b( para%h, bt, volh)                 ! B~ ~ inv(H+) B~
     endif

     call xbound_sc_field(a)
     call xbound_sc_field(b)
     call xbound_sc_field(at)
     call xbound_sc_field(bt)
  endif

  TIMING_START(timing_bin_dsf)

  s = -step * TWO * para%kappa**2 / (ONE + para%h**2)

  if (s /= ZERO) then
   call dsf_xf(p(1,1,EVEN,1), b, at, s, u(1,1,1,EVEN,1), nn(1,EVEN,1,FWD), VOLH)
   call dsf_xf(p(1,1,ODD ,1), bt, a, s, u(1,1,1,ODD ,1), nn(1,ODD ,1,FWD), VOLH)
   call dsf_xb(p(1,1,EVEN,1), bt, a, s, u(1,1,1,EVEN,1), nn(1,EVEN,1,FWD), VOLH)
   call dsf_xb(p(1,1,ODD ,1), b, at, s, u(1,1,1,ODD ,1), nn(1,ODD ,1,FWD), VOLH)

   call dsf_yf(p(1,1,EVEN,2), b, at, s, u(1,1,1,EVEN,2), nn(1,EVEN,2,FWD), VOLH)
   call dsf_yf(p(1,1,ODD ,2), bt, a, s, u(1,1,1,ODD ,2), nn(1,ODD ,2,FWD), VOLH)
   call dsf_yb(p(1,1,EVEN,2), bt, a, s, u(1,1,1,EVEN,2), nn(1,EVEN,2,FWD), VOLH)
   call dsf_yb(p(1,1,ODD ,2), b, at, s, u(1,1,1,ODD ,2), nn(1,ODD ,2,FWD), VOLH)

   call dsf_zf(p(1,1,EVEN,3), b, at, s, u(1,1,1,EVEN,3), nn(1,EVEN,3,FWD), VOLH)
   call dsf_zf(p(1,1,ODD ,3), bt, a, s, u(1,1,1,ODD ,3), nn(1,ODD ,3,FWD), VOLH)
   call dsf_zb(p(1,1,EVEN,3), bt, a, s, u(1,1,1,EVEN,3), nn(1,EVEN,3,FWD), VOLH)
   call dsf_zb(p(1,1,ODD ,3), b, at, s, u(1,1,1,ODD ,3), nn(1,ODD ,3,FWD), VOLH)

   call dsf_tf(p(1,1,EVEN,4), b, at, s, u(1,1,1,EVEN,4), nn(1,EVEN,4,FWD), VOLH)
   call dsf_tf(p(1,1,ODD ,4), bt, a, s, u(1,1,1,ODD ,4), nn(1,ODD ,4,FWD), VOLH)
   call dsf_tb(p(1,1,EVEN,4), bt, a, s, u(1,1,1,EVEN,4), nn(1,EVEN,4,FWD), VOLH)
   call dsf_tb(p(1,1,ODD ,4), b, at, s, u(1,1,1,ODD ,4), nn(1,ODD ,4,FWD), VOLH)
  endif

  TIMING_STOP(timing_bin_dsf)

  call flip_bc(u)

  s1 = -step * TWO * (para%csw_kappa / EIGHT)
  s2 = -step * TWO * (para%csw_kappa / EIGHT) * para%kappa**2
  
  if (s1 /= ZERO) call clover_dsf(EVEN, p, b,  a,  s1, u)
  if (s2 /= ZERO) call clover_dsf(ODD,  p, bt, at, s2, u)
end

!===============================================================================
