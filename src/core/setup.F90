!  Copyright 2020 University of Warwick

!  Licensed under the Apache License, Version 2.0 (the "License");
!  you may not use this file except in compliance with the License.
!  You may obtain a copy of the License at

!      http://www.apache.org/licenses/LICENSE-2.0

!  Unless required by applicable law or agreed to in writing, software
!  distributed under the License is distributed on an "AS IS" BASIS,
!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!  See the License for the specific language governing permissions and
!  limitations under the License.

MODULE setup

  USE shared_data
  USE welcome
  USE diagnostics

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: before_control, after_control
  PUBLIC :: grid

CONTAINS

  !****************************************************************************
  ! Any variables that need to have default values before the user specified
  ! ones in control.f90
  !****************************************************************************

  SUBROUTINE before_control

    ! Setup basic variables which have to have default values

    nprocx = 0
    nprocy = 0
    nprocz = 0

    time = 0.0_num
    gamma = 5.0_num / 3.0_num

  END SUBROUTINE before_control



  !****************************************************************************
  ! Variables which need to be specified after the control.f90 is run,
  ! MPI has been setup
  !****************************************************************************

  SUBROUTINE after_control

    ! Setup arrays and other variables which can only be set after
    ! user input

    IF (IAND(initial, IC_RESTART) == 0) restart_snapshot = 0

    p_visc = 0.0_num
    eta = 0.0_num
    grav = 0.0_num

    rho = 0.0_num
    energy = 0.0_num
    bx = 0.0_num
    by = 0.0_num
    bz = 0.0_num
    vx = 0.0_num
    vy = 0.0_num
    vz = 0.0_num

    va_max2 = va_max**2


  END SUBROUTINE after_control



  !****************************************************************************
  ! Stretched and staggered grid
  !****************************************************************************

  SUBROUTINE grid

    REAL(num) :: dx, dy, dz
    INTEGER :: ix, iy, iz

    length_x = x_max - x_min
    length_y = y_max - y_min
    length_z = z_max - z_min

    ! Initially assume uniform grid
    dx = length_x / REAL(nx_global, num)
    dy = length_y / REAL(ny_global, num)
    dz = length_z / REAL(nz_global, num)

    ! Grid cell boundary for x coordinates
    DO ix = -2, nx_global + 2
      xb_global(ix) = x_min + REAL(ix, num) * dx
    END DO

    IF (x_stretch) CALL stretch_x ! stretch grid ?

    ! Define position of ghost cells using sizes of adjacent cells
    IF (xbc_max == BC_PERIODIC) THEN
      xb_global(nx_global+1) = xb_global(nx_global) &
          + (xb_global(1) - xb_global(0))
      xb_global(nx_global+2) = xb_global(nx_global) &
          + (xb_global(2) - xb_global(0))
      xb_global(-1) = xb_global(0) &
          - (xb_global(nx_global) - xb_global(nx_global-1))
      xb_global(-2) = xb_global(0) &
          - (xb_global(nx_global) - xb_global(nx_global-2))
    ELSE
      xb_global(nx_global+1) = 2.0_num * xb_global(nx_global) &
          - xb_global(nx_global-1)
      xb_global(nx_global+2) = 2.0_num * xb_global(nx_global) &
          - xb_global(nx_global-2)
      xb_global(-1) = 2.0_num * xb_global(0) - xb_global(1)
      xb_global(-2) = 2.0_num * xb_global(0) - xb_global(2)
    END IF

    ! Setup local grid

    xb(-2) = xb_global(-2+n_global_min(1))

    DO ix = -1, nx + 2
      ixm = ix - 1

      xb(ix) = xb_global(ix+n_global_min(1))
      ! Cell centre
      xc(ix) = 0.5_num * (xb(ixm) + xb(ix))
      ! Cell width
      dxb(ix) = xb(ix) - xb(ixm)
    END DO

    DO ix = -1, nx + 1
      ixp = ix + 1
      ! Distance between centres
      dxc(ix) = xc(ixp) - xc(ix)
    END DO

    ! Repeat for y

    DO iy = -2, ny_global + 2
      yb_global(iy) = y_min + REAL(iy, num) * dy
    END DO

    IF (y_stretch) CALL stretch_y

    IF (ybc_max == BC_PERIODIC) THEN
      yb_global(ny_global+1) = yb_global(ny_global) &
          + (yb_global(1) - yb_global(0))
      yb_global(ny_global+2) = yb_global(ny_global) &
          + (yb_global(2) - yb_global(0))
      yb_global(-1) = yb_global(0) &
          - (yb_global(ny_global) - yb_global(ny_global-1))
      yb_global(-2) = yb_global(0) &
          - (yb_global(ny_global) - yb_global(ny_global-2))
    ELSE
      yb_global(ny_global+1) = 2.0_num * yb_global(ny_global) &
          - yb_global(ny_global-1)
      yb_global(ny_global+2) = 2.0_num * yb_global(ny_global) &
          - yb_global(ny_global-2)
      yb_global(-1) = 2.0_num * yb_global(0) - yb_global(1)
      yb_global(-2) = 2.0_num * yb_global(0) - yb_global(2)
    END IF

    yb(-2) = yb_global(-2+n_global_min(2))

    DO iy = -1, ny + 2
      iym = iy - 1

      yb(iy) = yb_global(iy+n_global_min(2))
      ! Cell centre
      yc(iy) = 0.5_num * (yb(iym) + yb(iy))
      ! Cell width
      dyb(iy) = yb(iy) - yb(iym)
    END DO

    DO iy = -1, ny + 1
      iyp = iy + 1
      dyc(iy) = yc(iyp) - yc(iy)
    END DO

    ! Repeat for z

    DO iz = -2, nz_global + 2
      zb_global(iz) = z_min + REAL(iz, num) * dz
    END DO

    IF (z_stretch) CALL stretch_z

    IF (zbc_max == BC_PERIODIC) THEN
      zb_global(nz_global+1) = zb_global(nz_global) &
          + (zb_global(1) - zb_global(0))
      zb_global(nz_global+2) = zb_global(nz_global) &
          + (zb_global(2) - zb_global(0))
      zb_global(-1) = zb_global(0) &
          - (zb_global(nz_global) - zb_global(nz_global-1))
      zb_global(-2) = zb_global(0) &
          - (zb_global(nz_global) - zb_global(nz_global-2))
    ELSE
      zb_global(nz_global+1) = 2.0_num * zb_global(nz_global) &
          - zb_global(nz_global-1)
      zb_global(nz_global+2) = 2.0_num * zb_global(nz_global) &
          - zb_global(nz_global-2)
      zb_global(-1) = 2.0_num * zb_global(0) - zb_global(1)
      zb_global(-2) = 2.0_num * zb_global(0) - zb_global(2)
    END IF

    zb(-2) = zb_global(-2+n_global_min(3))

    DO iz = -1, nz + 2
      izm = iz - 1

      zb(iz) = zb_global(iz+n_global_min(3))
      ! Cell centre
      zc(iz) = 0.5_num * (zb(izm) + zb(iz))
      ! Cell width
      dzb(iz) = zb(iz) - zb(izm)
    END DO

    DO iz = -1, nz + 1
      izp = iz + 1
      dzc(iz) = zc(izp) - zc(iz)
    END DO

    ! Define the cell area
    DO iz = -1, nz + 2
      DO iy = -1, ny + 2
        DO ix = -1, nx + 2
          cv(ix,iy,iz) = dxb(ix) * dyb(iy) * dzb(iz)
        END DO
      END DO
    END DO

  END SUBROUTINE grid



  !****************************************************************************
  ! Subroutine stretches the grid in the x direction
  ! Replace with any stretching algorithm as needed
  !****************************************************************************

  SUBROUTINE stretch_x

    REAL(num) :: width, dx, L, f, lx_new
    REAL(num), DIMENSION(:), ALLOCATABLE :: dxnew

    ALLOCATE(dxnew(-2:nx_global+2))

    ! New total length
    lx_new = 200.0_num

    ! Centre of tanh stretching in unstretched coordinates
    L = length_x / 1.5_num

    ! Width of tanh stretching in unstretched coordinates
    width = length_x / 10.0_num

    f = (lx_new - length_x) / (length_x - L) / 2.0_num

    dx = length_x / REAL(nx_global, num)
    dxnew(:) = dx + f * (1.0_num + TANH((ABS(xb_global(:)) - L) / width)) * dx

    DO ix = 1, nx_global + 2
      xb_global(ix) = xb_global(ix-1) + dxnew(ix)
    END DO

    length_x = lx_new

    DEALLOCATE(dxnew)

  END SUBROUTINE stretch_x



  !****************************************************************************
  ! Subroutine stretches the domain in the y direction
  ! Stretch domain upwards only
  !****************************************************************************

  SUBROUTINE stretch_y

    REAL(num) :: width, dy, L, f, ly_new
    REAL(num), DIMENSION(:), ALLOCATABLE :: dynew

    ALLOCATE(dynew(-2:ny_global+2))

    ! New total length
    ly_new = 100.0_num

    ! Centre of tanh stretching in unstretched coordinates
    L = length_y / 1.5_num

    ! Width of tanh stretching in unstretched coordinates
    width = length_y / 10.0_num

    f = (ly_new - length_y) / (length_y - L) / 2.0_num

    dy = length_y / REAL(ny_global, num)
    dynew(:) = dy + f * (1.0_num + TANH((ABS(yb_global(:)) - L) / width)) * dy

    DO iy = 1, ny_global + 2
      yb_global(iy) = yb_global(iy-1) + dynew(iy)
    END DO

    length_y = ly_new

    DEALLOCATE(dynew)

  END SUBROUTINE stretch_y



  !****************************************************************************
  ! Subroutine stretches the domain in the z direction
  ! Stretch domain upwards only
  !****************************************************************************

  SUBROUTINE stretch_z

    REAL(num) :: width, dz, L, f
    REAL(num), DIMENSION(:), ALLOCATABLE :: dznew

    ALLOCATE(dznew(-2:nz_global+2))

    L = 0.7_num * length_z 
    width = 0.1_num * length_z 
    dz = length_z / REAL(nz_global, num)
    f = 10.0_num
    dznew(:) = dz * (1.0_num + f * ( 1.0_num + TANH((ABS(zb_global(:)) - L) / width))) 

    DO iz = 1, nz_global + 2
      zb_global(iz) = zb_global(iz-1) + dznew(iz)
    END DO

    z_min = zb_global(0)
    z_max = zb_global(nz_global)
    length_z = z_max - z_min

    DEALLOCATE(dznew)

  END SUBROUTINE stretch_z


 !****************************************************************************
  ! Function to compare two strings.
  !****************************************************************************

  FUNCTION str_cmp(str_in, str_test)

    CHARACTER(*), INTENT(IN) :: str_in, str_test
    CHARACTER(30) :: str_trim
    LOGICAL :: str_cmp
    INTEGER :: l

    str_trim = TRIM(ADJUSTL(str_in))
    l = LEN(str_test)

    IF (l > LEN(str_in)) THEN
      str_cmp = .FALSE.
      RETURN
    END IF

    IF (str_trim(l+1:l+1) /= ' ') THEN
      str_cmp = .FALSE.
      RETURN
    END IF

    str_cmp = str_trim(1:l) == str_test

  END FUNCTION str_cmp


!   SUBROUTINE check_dims(dims)
!
!     INTEGER, INTENT(IN) :: dims(:)
!     INTEGER, DIMENSION(c_ndims) :: global_dims
!     INTEGER :: ierr
!
!     global_dims = (/ nx_global, ny_global, nz_global /)
!
!     IF (ALL(dims(1:c_ndims) == global_dims(1:c_ndims))) RETURN
!
!     IF (rank == 0) THEN
!       PRINT*, '*** ERROR ***'
!       PRINT*, 'Number of gridpoints in restart dump does not match', &
!           ' the control parameters.'
!       PRINT*, 'Control grid: ', nx_global, ',', ny_global, ',', nz_global
!       PRINT*, 'Restart dump grid: ', dims(1), ',', dims(2), ',', dims(3)
!     END IF
!
!     CALL MPI_ABORT(MPI_COMM_WORLD, ierr, errcode)
!     STOP
!
!   END SUBROUTINE check_dims

END MODULE setup
