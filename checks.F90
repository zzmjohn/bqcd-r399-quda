!===============================================================================
!
! checks.F90
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
subroutine check_csw(beta, csw)

  implicit none
  REAL, intent(in) :: beta, csw
  REAL             :: g, c

  if (beta == ZERO) return  
  if (csw == ZERO) return  

  g = SIX / beta
  c = ONE - 0.454 * g - 0.175 * g**2 + 0.012 * g**3 + 0.045 * g**4
  c = c / (ONE - 0.720 * g)

  if (abs(c - csw) > 0.00005) then
     call warn("check_csw(): c_sw differs more than 0.00005 from ALPHA value")
  endif

end

!-------------------------------------------------------------------------------
subroutine check_bc_fermions(bc_fermions, gamma_index)

  ! warns if the number of anti-periodic fermionic b.c. is 1 and
  ! the anti-periodic direction is not the gamma_4 direction
           
  implicit none
  integer, dimension(DIM), intent(in) :: bc_fermions, gamma_index

  integer :: i, i_anti, count

  count = 0
  do i = 1, DIM
     if (bc_fermions(i) < 0) then
        count = count + 1
        i_anti = i
     endif
  enddo

  if (count == 1 .and. gamma_index(i_anti) /= 4) then
   call warn("check_bc_fermions(): anti-periodic b.c. not in gamma_4 direction")
  endif

end
!===============================================================================
