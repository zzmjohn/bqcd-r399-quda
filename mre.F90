!===============================================================================
!
! mre.F90 - Chronological Inverter by Minimal Residual Extrapolation
!           hep-lat/9509012
!
!-------------------------------------------------------------------------------
!
! Copyright (C) 2003 Hinnerk Stueben
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
subroutine mre_put(basis, sc_field, reset)  ! add a solution

  use module_mre
  use module_p_interface
  use module_vol
  implicit none

  type(type_mre), intent(inout) :: basis
  SPINCOL_FIELD,  intent(in)    :: sc_field
  integer,        intent(in)    :: reset

  P_SPINCOL_FIELD               :: tmp
  integer                       :: i

  if (mre_n_vec == 0) then
     return
  endif

  call mre_allocate(basis)

  if (reset /= 0) then
     basis%rank = 0
     return
  endif

  tmp => basis%vec(mre_n_vec)%sc

  do i = mre_n_vec, 2, -1
     basis%vec(i)%sc => basis%vec(i - 1)%sc
  enddo

  basis%vec(1)%sc => tmp
  call sc_copy(basis%vec(1)%sc, sc_field)

  if (basis%rank < mre_n_vec) basis%rank = basis%rank + 1

end

!-------------------------------------------------------------------------------
subroutine mre_get(basis, matrix_mult, trial, phi, para, conf)

  ! get trial solution

  use typedef_hmc
  use module_function_decl
  use module_mre
  use module_vol
  implicit none

  type(type_mre), intent(inout) :: basis
  external                      :: matrix_mult
  SPINCOL_FIELD,  intent(out)   :: trial
  SPINCOL_FIELD,  intent(in)    :: phi
  type(hmc_para), intent(in)    :: para
  type(hmc_conf), intent(in)    :: conf

  type(type_mre), save          :: mv  ! "Matrix * v"
  COMPLEX, target               :: g(mre_n_vec, mre_n_vec + 1)
  COMPLEX, pointer              :: b(:)
  integer                       :: size_g
  integer                       :: i
  integer                       :: j
  integer                       :: s
  integer                       :: c
  integer                       :: rest

  if (mre_n_vec == 0 .or. .not. associated(basis%vec) .or. basis%rank == 0) then
     call sc_copy(trial, phi)
     return
  endif

  if (basis%rank == 1) then
     call sc_copy(trial, basis%vec(1)%sc)
     return
  endif

  size_g = mre_n_vec * (mre_n_vec + 1) * SIZE_COMPLEX

  b => g(:, mre_n_vec + 1)  ! storage arrangement for global sum
                            ! => one global sum for everything

  call mre_allocate(mv)
  call mre_gram_schmidt(basis)
  
  do i = 1, basis%rank
     b(i) = sc_cdotc(basis%vec(i)%sc, phi)
     call matrix_mult(mv%vec(i)%sc, basis%vec(i)%sc, para, conf)
  enddo

  do i = 1, basis%rank
     g(i, i) = sc_norm2(mv%vec(i)%sc)
     do j = i + 1, basis%rank
        g(i, j) = sc_cdotc(mv%vec(i)%sc, mv%vec(j)%sc)
        g(j, i) = conjg(g(i, j))
     enddo
  enddo

  call global_sum_vec(size_g, g)

  call mre_gauss_jordan(g, b, basis%rank, mre_n_vec)

  ! calculation of "trial" with doubled data re-use:

  call sc_cax2(trial, basis%vec(1)%sc, b(1), basis%vec(2)%sc, b(2))

  rest = mod(basis%rank, 2)

  do j = 3, basis%rank - rest, 2
     call sc_caxpy2(trial, basis%vec(j)%sc,   b(j), &
                           basis%vec(j+1)%sc, b(j+1))
  enddo

  if (rest == 1) then
     j = basis%rank
     call sc_caxpy(trial, basis%vec(j)%sc, b(j))
  endif

end

!-------------------------------------------------------------------------------
subroutine mre_allocate(basis)

  use module_mre
  use module_p_interface
  use module_vol
  implicit none

  type(type_mre), intent(inout) :: basis
  integer                       :: i

  if (.not. associated(basis%vec)) then
     allocate(basis%vec(mre_n_vec))
     do i = 1, mre_n_vec
        nullify(basis%vec(i)%sc)  
        call allocate_sc_field(basis%vec(i)%sc)
     enddo
     basis%rank = 0
  endif

end

!-------------------------------------------------------------------------------
subroutine mre_gram_schmidt(basis)

  ! Golub and van Loon, Matrix Computations (3rd ed.), p. 232

  use module_function_decl
  use module_mre
  use module_vol
  implicit none

  type(type_mre), intent(inout) :: basis

  integer :: k, j
  REAL    :: r_kk, r_kj

  do k = 1, basis%rank
      r_kk = sc_norm2(basis%vec(k)%sc)
      r_kk = global_sum(r_kk)
      r_kk = ONE / sqrt(r_kk)
      call sc_scale(basis%vec(k)%sc, r_kk)
      do j = k + 1, basis%rank
         r_kj = sc_dot(basis%vec(k)%sc, basis%vec(j)%sc)
         r_kj = global_sum(r_kj)
         call sc_axpy(basis%vec(j)%sc, basis%vec(k)%sc, -r_kj)
      enddo
  enddo

end

!-------------------------------------------------------------------------------
subroutine mre_gauss_jordan(a, b, n, np)

  ! Numerical Recipes in Fortran (2nd ed.), p. 30
  ! inv(a) i.e. a is not unscrambled

  implicit none

  integer, intent(in)    :: n, np
  COMPLEX, intent(inout) :: a(np,np), b(np)

  integer :: ipiv(n)
  integer :: i, j, k, l, ll
  integer :: irow, icol

  REAL    :: big, tmp
  COMPLEX :: dum
  COMPLEX :: pivinv

  ipiv = 0

  do i = 1, n

     big = ZERO
     do j = 1, n
        if (ipiv(j) /= 1) then
           do k = 1, n
              if (ipiv(k) == 0) then
                 tmp = abs(a(j, k))
                 if (tmp >= big) then
                    big = tmp
                    irow = j
                    icol = k
                 endif
              else if (ipiv(k) > 1) then
                 call die("mre_gauss_jordan(): singular matrix 1")
              endif
           enddo
        endif
     enddo

     ipiv(icol) = ipiv(icol) + 1

     if (irow /= icol) then
        do l = 1, n
           dum = a(irow, l)
           a(irow, l) = a(icol, l)
           a(icol, l) = dum
        enddo
        dum = b(irow)
        b(irow) = b(icol)
        b(icol) = dum
     endif

     if (a(icol, icol) == ZERO ) then
        call die("mre_gauss_jordan(): singular matrix 2")
     endif

     pivinv = ONE / a(icol, icol)
     !!a(icol, icol) = ONE  !! only needed for inv(a)

     do  l = 1, n
        a(icol, l) = a(icol, l) * pivinv
     enddo

     b(icol) = b(icol) * pivinv

     do ll = 1, n
        if (ll /= icol) then
           dum = a(ll, icol)
           a(ll, icol) = ZERO
           do l = 1, n
              a(ll, l) = a(ll, l) - a(icol, l) * dum
           enddo
           b(ll) = b(ll) - b(icol) * dum
        endif
     enddo

  enddo

end

!===============================================================================
