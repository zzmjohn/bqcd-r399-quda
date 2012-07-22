!===============================================================================
!
! index.F90 - functions for index calculations
!             these functions work "stand alone"
!
!-------------------------------------------------------------------------------
!
! Copyright (C) 1998-2006 Hinnerk Stueben
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
      integer function i_e_o(dim, i)  ! returns EVEN or ODD (0 or 1)

      implicit none
      integer dim, i(dim), d
 
      i_e_o = 0
      do d = 1, dim
         i_e_o = i_e_o + i(d)
      enddo
 
      i_e_o = mod(abs(i_e_o), 2)
      end
 
!-------------------------------------------------------------------------------
      integer function i_global(i_local, N_local, i_pe)

      implicit none
      integer i_local, N_local, i_pe

      i_global = i_pe * N_local + i_local

      end
 
!-------------------------------------------------------------------------------
      integer function i_local(i_global, N_local, i_pe)

      implicit none
      integer i_global, N_local, i_pe

      i_local = i_global - i_pe * N_local

      end
 
!-------------------------------------------------------------------------------
      integer function i_periodic(i, n)
 
      implicit none
      integer i, n

      if (i .ge. 0) then
         if (i .lt. n) then
            i_periodic = i
         elseif (i .lt. 2 * n) then
            i_periodic = i - n
         else
            i_periodic = mod(i, n)
         endif
      else
         if (i .ge. -n) then
            i_periodic = i + n
         else
            i_periodic =  i + (1 - (i + 1) / n) * n
         endif
      endif
      end
 
!-------------------------------------------------------------------------------
      integer function ieo(dim, i, n)
 
      implicit none
      integer dim, i(dim), n(dim), d, ilex
 
      ieo = ilex(dim, i, n) / 2
      end
 
!-------------------------------------------------------------------------------
      integer function ilex(dim, i, n)
 
      implicit none
      integer dim, i(dim), n(dim), d
 
      ilex = i(dim)
      do d = dim - 1, 1, -1
         ilex = ilex * n(d) + i(d)
      enddo
 
      end
 
!-------------------------------------------------------------------------------
      integer function n_sites(dim, direction, n, npe)
 
!     returns number of sites of local grid and boundaries
!     n_sites(dim, (/0, 0, ..., 0/), n, npe) is the (local) grid volume
 
      implicit none
      integer dim, direction(dim), n(dim), npe(dim), d
 
      n_sites = 1
      do d = 1, dim
         if (direction(d) .eq. 0) then
            n_sites = n_sites * n(d)
         else
            if (npe(d) .eq. 1) then ! grid not partitioned in d-direction
               n_sites = 0
               return
            endif
         endif
      enddo
 
      end
 
!-------------------------------------------------------------------------------
!!      subroutine uneo(ieo, eo, dim, i, n)  ! returns i for given (ieo, eo)
!!
!!      implicit none
!!      integer ieo, eo, dim, i(dim), n(dim), e_o
!!
!!      call unlex(2 * ieo, dim, i, n)
!!      i(1) =  i(1) + ieor(e_o(dim-1, i(2)), eo)
!!
!!      end
!!
!-------------------------------------------------------------------------------
      subroutine unlex(ilex, dim, i, n)
 
      ! remember the range of ilex: 0 <= ilex < (n(1) * ... * n(dim))
 
      implicit none
      integer ilex, dim, i(dim), n(dim), j, d
 
      j = ilex
      do d = 1, dim - 1
         i(d) = mod(j, n(d))
         j = j / n(d)  ! integer division
      enddo
      i(dim) = j
 
      end

!===============================================================================
