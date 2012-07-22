!===============================================================================
!
! polyakov_loop.F90  -  in gamma_4-direction, requires (NPE(gamma_4) == 1)
!
!-------------------------------------------------------------------------------
!
! Copyright (C) 2000-2006 Hinnerk Stueben
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
subroutine polyakov_loop(conf, traj, i_ensemble1, i_ensemble2)

  use typedef_hmc
  use module_decomp
  use module_function_decl
  use module_vol
  implicit none

  integer,        intent(in)  :: traj, i_ensemble1, i_ensemble2
  type(hmc_conf), intent(in)  :: conf

  COMPLEX                     :: pl
  REAL                        :: re_pl, im_pl
  integer                     :: x, y, z, t, i, eo, j(DIM)
  integer                     :: nx, ny, nz, nt, dir4, npe4
  integer, external           :: ieo, e_o, std_xyzt2i
  SU3                         :: u

  SU3, parameter :: su3_one = reshape( &
                           (/ ONE,ZERO,ZERO, &
                              ZERO,ONE,ZERO, &
                              ZERO,ZERO,ONE /), &
                           (/ NCOL, NCOL /))

  character(len=*), parameter :: key_pl = "%pl"
  integer,        save        :: count = 0


  count = count + 1

  dir4 = decomp%direction(4)  
  npe4 = decomp%act%npe(dir4)

  if (npe4 /= 1) then
     call die("polyakov_loop(): gamma_4 direction must not be decomposed")
  endif

  nx = decomp%std%N(1)
  ny = decomp%std%N(2)
  nz = decomp%std%N(3)
  nt = decomp%std%N(4)

  pl = 0
  !$omp parallel do reduction(+: pl) private(x, y, z, t, i, j, eo, u)
  do x = 0, nx - 1
     do y = 0, ny - 1
        do z = 0, nz - 1

           u = su3_one
           do t = 0, nt - 1
              j = (/x, y, z, t/) 

              i = std_xyzt2i(j)
              eo = e_o(j)

              call u_update2(u, conf%u(1, 1, i, eo, 4))
           enddo

           pl = pl + u(1,1) + u(2,2) + u(3,3)
        enddo
     enddo
  enddo

  call global_sum_vec(SIZE_COMPLEX, pl)

  pl = pl / (THREE * decomp%std%L(1) * decomp%std%L(2) * decomp%std%L(3))

  re_pl = Re(pl)
  im_pl = Im(pl)

  if (my_pe() == 0) then
     if (count == 1) write(UREC, 400) &
         "T", key_pl, "traj", "e", "f", "Re(Polyakov_Loop)", "Im(Polyakov_Loop)"

     write(UREC, 410) key_pl, traj, i_ensemble1, i_ensemble2, re_pl, im_pl
  endif


400 format (1x, 2a, a6, 2a3, 2a20)
410 format (1x, a4, i6, 2i3, 2g20.10)

end
