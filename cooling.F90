!===============================================================================
!
! cooling.F90 - measurement of the topological charge using standard cooling
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
# include "defs.h"

!-------------------------------------------------------------------------------
module module_cooling

  integer, save                        :: n_cool
  logical, dimension(:), pointer, save :: measure_q

end

!-------------------------------------------------------------------------------
subroutine init_cooling(list)

  use module_cooling
  implicit none
  character(*), intent(in) :: list
  integer                  :: i, iostat

  if (list /= " ") then
     open(ULIST, file = list, action = "read", status = "old")

     iostat = 0
     n_cool = 0
     do while (iostat == 0)
        read(ULIST, *, iostat = iostat) i
        if (i < 0) then
           call die("init_cooling(): list has negative entries")
        else
           n_cool = max(n_cool, i)
        endif
     enddo
     
     allocate(measure_q(0:n_cool))
     
     do i = 0, n_cool
        measure_q(i) = .false.
     enddo
     
     rewind(ULIST)
     
     iostat = 0
     do while (iostat == 0)
        read(ULIST, *, iostat = iostat) i
        measure_q(i) = .true.
     enddo
     
     close(ULIST)

  else
     n_cool = -1
  endif

end

!-------------------------------------------------------------------------------
subroutine cooling(u_in, traj, i_ensemble1, i_ensemble2)

  use module_cooling
  use module_function_decl
  use module_p_interface
  use module_vol
  implicit none

  integer, intent(in)           :: traj, i_ensemble1, i_ensemble2
  GAUGE_FIELD, intent(in)       :: u_in
  P_GAUGE_FIELD, save           :: u
  integer                       :: i
  character(len = *), parameter :: key = "%Qc"
  REAL                          :: q, plaq

  if (n_cool < 0) return

  TIMING_START(timing_bin_cooling)

  ALLOCATE_G_FIELD(u)

  u = u_in

  call begin(UREC, "Cooling")

  if (my_pe() == 0) then
     write(UREC, 400) "T", key, "traj", "e", "f", "i_cool", &
          "Q_cool", "PlaqEnergy"
  endif

  do i = 0, n_cool
     if (measure_q(i)) then
        call conf_check(u)
        call top_charge(q, u, plaq)
        if (my_pe() == 0) then
           write(UREC, 410) key, traj, i_ensemble1, i_ensemble2, i, q, plaq
        endif
     endif
     if (i < n_cool) call conf_relax(u)
  enddo

  call end(UREC, "Cooling")

400 format (1x, 2a, a6, 2a3, a8, a15,   a15)
410 format (1x, a4, i6, 2i3, i8, f15.6, f15.10)

  TIMING_STOP(timing_bin_cooling)

end

!-------------------------------------------------------------------------------
subroutine conf_relax(u)

  use module_vol
  implicit none

  GAUGE_FIELD, intent(inout) :: u
  SU3 :: uuu, w, a
  SU3, parameter :: su3_one = reshape( &
                           (/ ONE,ZERO,ZERO, &
                              ZERO,ONE,ZERO, &
                              ZERO,ZERO,ONE /), &
                           (/ NCOL, NCOL /))
  REAL :: p0, p1, p2, p3, fac
  REAL :: a0, a1, a2, a3
  integer :: i, eo, mu, k, c1, c2

  do mu = 1, DIM
     do eo = EVEN, ODD
        !$omp parallel do private(uuu, w, a, p0, p1, p2, p3, fac, &
        !$omp a0, a1, a2, a3, k, c1, c2)
        do i = 1, VOLH
           call staple(uuu, u, i, eo, mu)
           do k = 1, NCOL
              if (k == 1) then
                 c1 = 1
                 c2 = 2
              else if (k == 2) then
                 c1 = 1
                 c2 = 3
              else if (k == 3) then
                 c1 = 2
                 c2 = 3
              endif
              
              call uu(w, u(1, 1, i, eo, mu), uuu)  ! w = u * uuu
           
              p0 =   Re(w(c1, c1)) + Re(w(c2, c2))
              p1 = -(Im(w(c1, c2)) + Im(w(c2, c1)))
              p2 = -(Re(w(c1, c2)) - Re(w(c2, c1)))
              p3 = -(Im(w(c1, c1)) - Im(w(c2, c2)))

              fac = ONE / sqrt(p0**2 + p1**2 + p2**2 + p3**2)

              a0 = fac * p0
              a1 = fac * p1
              a2 = fac * p2
              a3 = fac * p3
              
              a = su3_one
              
              a(c1, c1) = cmplx( a0,  a3, kind = RKIND)
              a(c1, c2) = cmplx( a2,  a1, kind = RKIND)
              a(c2, c1) = cmplx(-a2,  a1, kind = RKIND)
              a(c2, c2) = cmplx( a0, -a3, kind = RKIND)
              
              call u_update(u(1, 1, i, eo, mu), a)  ! u = a * u

           enddo ! k
        enddo ! i
        call xbound_g(u, eo, mu)
     enddo ! eo
  enddo ! mu
end

!-------------------------------------------------------------------------------
subroutine top_charge(qq, u, plaq_energy)

  use module_vol
  use module_nn
  implicit none

  REAL, intent(out)       :: qq, plaq_energy
  GAUGE_FIELD, intent(in) :: u

  integer :: e, o, mu, nu, i, j1, j2, j3, j4, j5, j6, j7, c1, c2
  SU3     :: uuu, left, right
  COMPLEX, dimension(NCOL, NCOL, DIM - 1, DIM) :: ut  ! U~(x,mu,nu) - h.c.
  REAL    :: q, plaq
  REAL, external :: global_sum, Re_Tr_uu

  !----------------------------------------------------------------------
  !
  !  (j3, e) --->--- (j2, o) ---<---    x
  !     |               |               |
  !     |               |               |
  !     ^               v               ^        nu
  !     |               |               |         ^
  !     |               |               |         |
  !  (j4, o) ---<---  (i,e)  --->--- (j1, o)      x-->  mu
  !     |               |               |
  !     |               |               |
  !     v               ^               v
  !     |               |               |
  !     |               |               |
  !  (j5, e) --->--- (j6, o) ---<--- (j7, e)
  !
  !----------------------------------------------------------------------

  q = 0
  plaq = 0

  do e = EVEN, ODD
     o = EVEN + ODD - e
     !$omp parallel do reduction(+: q, plaq) private(uuu, left, right, ut, &
     !$omp mu, nu, i, j1, j2, j3, j4, j5, j6, j7, c1, c2)  
     do i = 1, VOLH
        do mu = 1, DIM - 1
           do nu = mu + 1, DIM

              j1 = nn(i,  e, mu, FWD)
              j2 = nn(i,  e, nu, FWD)
              j3 = nn(j2, o, mu, BWD)
              j4 = nn(j3, e, nu, BWD)
              j5 = nn(j4, o, nu, BWD)
              j6 = nn(j5, e, mu, FWD)
              j7 = nn(j6, o, mu, FWD)

              uuu = 0

              call uuu_fwd(uuu, u(1, 1, j1, o, nu), &
                                u(1, 1, j2, o, mu), &
                                u(1, 1, i,  e, nu))

              plaq = plaq + Re_Tr_uu(uuu, u(1, 1, i, e, mu))

              call uuu_bwd_m(uuu, u(1, 1, j7, e, nu), &
                                  u(1, 1, j6, o, mu), &
                                  u(1, 1, j6, o, nu))

              call uu(right, u(1, 1, i, e, mu), uuu)

              uuu = 0

              call uuu_fwd(uuu, u(1, 1, i,  e, nu), &
                                u(1, 1, j3, e, mu), &
                                u(1, 1, j4, o, nu))

              call uuu_bwd_m(uuu, u(1, 1, j6, o, nu), &
                                  u(1, 1, j5, e, mu), &
                                  u(1, 1, j5, e, nu))

              call uu(left, uuu, u(1, 1, j4, o, mu))

              do c2 = 1, NCOL
                 do c1 = 1, NCOL
                    ut(c1, c2, mu, nu) = right(c1, c2) - conjg(right(c2, c1)) &
                                       +  left(c1, c2) - conjg(left(c2, c1))
                 enddo
              enddo
           enddo ! nu
        enddo ! mu

        q = q + Re_Tr_uu(ut(1, 1, 1, 2), ut(1, 1, 3, 4)) &
              - Re_Tr_uu(ut(1, 1, 1, 3), ut(1, 1, 2, 4)) &
              + Re_Tr_uu(ut(1, 1, 1, 4), ut(1, 1, 2, 3))

     enddo ! i
  enddo ! e/o
  
  q = global_sum(q)
  plaq = global_sum(plaq)
  
  q = -q / (256 * PI**2)
  qq = q
  plaq_energy = ONE - plaq / (THREE * SIX * volume)
end

!===============================================================================
!===============================================================================
