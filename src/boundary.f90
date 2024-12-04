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

!******************************************************************************
! This module contains the boundary conditions for the entire code
! Any new boundary conditions should be added here
!******************************************************************************

MODULE boundary

  USE shared_data
  USE mpiboundary

  IMPLICIT NONE

CONTAINS

  !****************************************************************************
  ! Set up any necessary variables for the chosen boundary conditions
  !****************************************************************************

  SUBROUTINE set_boundary_conditions

    any_open = .FALSE.
    IF (xbc_min == BC_OPEN .OR. xbc_max == BC_OPEN &
        .OR. ybc_min == BC_OPEN .OR. ybc_max == BC_OPEN &
        .OR. zbc_min == BC_OPEN .OR. zbc_max == BC_OPEN) any_open = .TRUE.

  END SUBROUTINE set_boundary_conditions



  !****************************************************************************
  ! Call all of the boundaries needed by the core Lagrangian solver
  !****************************************************************************

  SUBROUTINE boundary_conditions

    CALL bfield_bcs
    CALL energy_bcs
    CALL density_bcs
    CALL velocity_bcs

  END SUBROUTINE boundary_conditions



  !****************************************************************************
  ! Boundary conditions for magnetic field through plane
  !****************************************************************************

  SUBROUTINE bfield_bcs

    CALL bfield_mpi

    IF (proc_x_min == MPI_PROC_NULL .AND. xbc_min == BC_USER) THEN
       bx(-1,1:ny,:) = bx(1,1:ny,:)
       bx(-2,1:ny,:) = bx(2,1:ny,:)
       by(0,1:ny, :) = by(1,1:ny,:)
       by(-1,1:ny,:) = by(2,1:ny,:)
       bz(0,1:ny, :) = bz(1,1:ny,:)
       bz(-1,1:ny,:) = bz(2,1:ny,:)
    END IF

    IF (proc_x_max == MPI_PROC_NULL .AND. xbc_max == BC_USER) THEN
       bx(nx+1,1:ny,:) = bx(nx,1:ny,:)
       bx(nx+2,1:ny,:) = bx(nx,1:ny,:)
       by(nx+1,1:ny,:) = by(nx,1:ny,:)
       by(nx+2,1:ny,:) = by(nx-1,1:ny,:)
       bz(nx+1,1:ny,:) = bz(nx,1:ny,:)
       bz(nx+2,1:ny,:) = bz(nx-1,1:ny,:)
    END IF

    IF (proc_y_min == MPI_PROC_NULL .AND. ybc_min == BC_USER) THEN
       bx(:,0, :)  = bx(:,1,:)
       bx(:,-1,:)  = bx(:,2,:)
       by(:,-1,:)  = by(:,1,:)
       by(:,-2,:)  = by(:,2,:)
       bz(:,0, :)  = bz(:,1,:)
       bz(:,-1,:)  = bz(:,2,:)
    END IF

    IF (proc_y_max == MPI_PROC_NULL .AND. ybc_max == BC_USER) THEN
       bx(:,ny+1,:) = bx(:,ny,:)
       bx(:,ny+2,:) = bx(:,ny-1,:)
       by(:,ny+1,:) = by(:,ny-1,:)
       by(:,ny+2,:) = by(:,ny-2,:)
       bz(:,ny+1,:) = bz(:,ny,:)
       bz(:,ny+2,:) = bz(:,ny-1,:)
    END IF

    IF (proc_z_min == MPI_PROC_NULL .AND. zbc_min == BC_USER) THEN
       bx(1:nx,1:ny,-1) = bx(1:nx,1:ny,2)
       bx(1:nx,1:ny,0 ) = bx(1:nx,1:ny,1)
       by(1:nx,1:ny,-1) = by(1:nx,1:ny,2)
       by(1:nx,1:ny,0 ) = by(1:nx,1:ny,1)
       bz(1:nx,1:ny,-1) = bz(1:nx,1:ny,1)
       bz(1:nx,1:ny,-2) = bz(1:nx,1:ny,2)
    END IF

    IF (proc_z_max == MPI_PROC_NULL .AND. zbc_max == BC_USER) THEN
       bx(1:nx,1:ny,nz+1) = bx(1:nx,1:ny,nz)
       bx(1:nx,1:ny,nz+2) = bx(1:nx,1:ny,nz-1)
       by(1:nx,1:ny,nz+1) = by(1:nx,1:ny,nz)
       by(1:nx,1:ny,nz+2) = by(1:nx,1:ny,nz-1)
       bz(1:nx,1:ny,nz+1) = bz(1:nx,1:ny,nz-1)
       bz(1:nx,1:ny,nz+2) = bz(1:nx,1:ny,nz-2)
    END IF

  END SUBROUTINE bfield_bcs



  !****************************************************************************
  ! Boundary conditions for specific internal energy
  !****************************************************************************

  SUBROUTINE energy_bcs

    CALL energy_mpi

    IF (proc_x_min == MPI_PROC_NULL .AND. xbc_min == BC_USER) THEN
      energy( 0,:,:) = energy(1,:,:)
      energy(-1,:,:) = energy(2,:,:)
    END IF

    IF (proc_x_max == MPI_PROC_NULL .AND. xbc_max == BC_USER) THEN
      energy(nx+1,:,:) = energy(nx  ,:,:)
      energy(nx+2,:,:) = energy(nx-1,:,:)
    END IF

    IF (proc_y_min == MPI_PROC_NULL .AND. ybc_min == BC_USER) THEN
      energy(:, 0,:) = energy(:,1,:)
      energy(:,-1,:) = energy(:,2,:)
    END IF

    IF (proc_y_max == MPI_PROC_NULL .AND. ybc_max == BC_USER) THEN
      energy(:,ny+1,:) = energy(:,ny  ,:)
      energy(:,ny+2,:) = energy(:,ny-1,:)
    END IF

    IF (proc_z_min == MPI_PROC_NULL .AND. zbc_min == BC_USER) THEN
      energy(:,:, 0) = energy(:,:,1)
      energy(:,:,-1) = energy(:,:,2)
    END IF

    IF (proc_z_max == MPI_PROC_NULL .AND. zbc_max == BC_USER) THEN
      energy(:,:,nz+1) = energy(:,:,nz)
      energy(:,:,nz+2) = energy(:,:,nz-1)
    END IF

  END SUBROUTINE energy_bcs



  !****************************************************************************
  ! Boundary conditions for density
  !****************************************************************************

  SUBROUTINE density_bcs

    CALL density_mpi

    IF (proc_x_min == MPI_PROC_NULL .AND. xbc_min == BC_USER) THEN
      rho( 0,:,:) = rho(1,:,:)
      rho(-1,:,:) = rho(2,:,:)
    END IF

    IF (proc_x_max == MPI_PROC_NULL .AND. xbc_max == BC_USER) THEN
      rho(nx+1,:,:) = rho(nx  ,:,:)
      rho(nx+2,:,:) = rho(nx-1,:,:)
    END IF

    IF (proc_y_min == MPI_PROC_NULL .AND. ybc_min == BC_USER) THEN
      rho(:, 0,:) = rho(:,1,:)
      rho(:,-1,:) = rho(:,2,:)
    END IF

    IF (proc_y_max == MPI_PROC_NULL .AND. ybc_max == BC_USER) THEN
      rho(:,ny+1,:) = rho(:,ny  ,:)
      rho(:,ny+2,:) = rho(:,ny-1,:)
    END IF

    IF (proc_z_min == MPI_PROC_NULL .AND. zbc_min == BC_USER) THEN
      rho(:,:, 0) = rho(:,:,1)
      rho(:,:,-1) = rho(:,:,2)
    END IF

    IF (proc_z_max == MPI_PROC_NULL .AND. zbc_max == BC_USER) THEN
      rho(:,:,nz+1) = rho(:,:,nz)
      rho(:,:,nz+2) = rho(:,:,nz-1)
    END IF

  END SUBROUTINE density_bcs


  !****************************************************************************
  ! Boundary conditions for temperature
  !****************************************************************************

  SUBROUTINE temperature_bcs

    CALL density_mpi

    IF (proc_x_min == MPI_PROC_NULL .AND. xbc_min == BC_USER) THEN
      temperature( 0,:,:) = temperature(1,:,:)
      temperature(-1,:,:) = temperature(2,:,:)
    END IF

    IF (proc_x_max == MPI_PROC_NULL .AND. xbc_max == BC_USER) THEN
      temperature(nx+1,:,:) = temperature(nx  ,:,:)
      temperature(nx+2,:,:) = temperature(nx-1,:,:)
    END IF

    IF (proc_y_min == MPI_PROC_NULL .AND. ybc_min == BC_USER) THEN
      temperature(:, 0,:) = temperature(:,1,:)
      temperature(:,-1,:) = temperature(:,2,:)
    END IF

    IF (proc_y_max == MPI_PROC_NULL .AND. ybc_max == BC_USER) THEN
      temperature(:,ny+1,:) = temperature(:,ny  ,:)
      temperature(:,ny+2,:) = temperature(:,ny-1,:)
    END IF

    IF (proc_z_min == MPI_PROC_NULL .AND. zbc_min == BC_USER) THEN
      temperature(:,:, 0) = temperature(:,:,1)
      temperature(:,:,-1) = temperature(:,:,2)
    END IF

    IF (proc_z_max == MPI_PROC_NULL .AND. zbc_max == BC_USER) THEN
      temperature(:,:,nz+1) = temperature(:,:,nz  )
      temperature(:,:,nz+2) = temperature(:,:,nz-1)
    END IF

  END SUBROUTINE temperature_bcs




  !****************************************************************************
  ! Full timestep velocity boundary conditions
  !****************************************************************************

  SUBROUTINE velocity_bcs

    INTEGER:: i,j

    CALL velocity_mpi

    IF (proc_x_min == MPI_PROC_NULL .AND. xbc_min == BC_USER) THEN
      vx(-2:0,:,:) = 0.0_num
      vy(-2:0,:,:) = 0.0_num
      vz(0,:,:) = vz(1,:,:)
      vz(-1,:,:) = vz(0,:,:)
      vz(-2,:,:) = vz(-1,:,:)
    END IF

    IF (proc_x_max == MPI_PROC_NULL .AND. xbc_max == BC_USER) THEN
      vx(nx:nx+2,:,:) = 0.0_num
      vy(nx:nx+2,:,:) = 0.0_num
      vz(nx,:,:) = vz(nx-1,:,:)
      vz(nx+1,:,:) = vz(nx,:,:)
      vz(nx+2,:,:) = vz(nx+2,:,:)
    END IF

    IF (proc_y_min == MPI_PROC_NULL .AND. ybc_min == BC_USER) THEN
      vx(:,-2:0,:) = 0.0_num
      vy(:,-2:0,:) = 0.0_num
      vz(:,0,:) = vz(:,1,:)
      vz(:,-1,:) = vz(:,0,:)
      vz(:,-2,:) = vz(:,-1,:)
    END IF

    IF (proc_y_max == MPI_PROC_NULL .AND. ybc_max == BC_USER) THEN
      vx(:,ny:ny+2,:) = 0.0_num
      vy(:,ny:ny+2,:) = 0.0_num
      vz(:,ny,:) = vz(:,ny-1,:)
      vz(:,ny+1,:) = vz(:,ny,:)
      vz(:,ny+2,:) = vz(:,ny+1,:)
    END IF


    IF (proc_z_max == MPI_PROC_NULL .AND. zbc_max == BC_USER) THEN
      vx(:,:,nz:nz+2) = 0.0_num
      vy(:,:,nz:nz+2) = 0.0_num

      do i = -2, nx+2
          do j = -2, ny+2
          vz(i, j, nz) = max(0.0,1.25*vz(i,j,nz-1))
          end do
      end do
      vz(:,:,nz+1) = vz(:,:,nz)

    END IF

    IF (proc_z_min == MPI_PROC_NULL .AND. zbc_min == BC_USER) THEN
      vx(:,:,-1) = vx(:,:,0)
      vx(:,:,-2) = vx(:,:,-1)
      vy(:,:,-1) = vy(:,:,0)
      vy(:,:,-2) = vy(:,:,-1)
      vz(:,:,-2:0) = 0.0_num
    END IF

  END SUBROUTINE velocity_bcs


  !****************************************************************************
  ! Half timestep velocity boundary conditions
  !****************************************************************************

  SUBROUTINE remap_v_bcs

    INTEGER:: i,j

    CALL remap_v_mpi

    IF (proc_x_min == MPI_PROC_NULL .AND. xbc_min == BC_USER) THEN
      vx1(-2:0,:,:) = 0.0_num
      vy1(-2:0,:,:) = 0.0_num
      vz1(0,:,:) = vz1(1,:,:)
      vz1(-1,:,:) = vz1(0,:,:)
      vz1(-2,:,:) = vz1(-1,:,:)
    END IF

    IF (proc_x_max == MPI_PROC_NULL .AND. xbc_max == BC_USER) THEN
      vx1(nx:nx+2,:,:) = 0.0_num
      vy1(nx:nx+2,:,:) = 0.0_num
      vz1(nx,:,:) = vz1(nx-1,:,:)
      vz1(nx+1,:,:) = vz1(nx,:,:)
      vz1(nx+2,:,:) = vz1(nx+2,:,:)
    END IF

    IF (proc_y_min == MPI_PROC_NULL .AND. ybc_min == BC_USER) THEN
      vx1(:,-2:0,:) = 0.0_num
      vy1(:,-2:0,:) = 0.0_num
      vz1(:,0,:) = vz1(:,1,:)
      vz1(:,-1,:) = vz1(:,0,:)
      vz1(:,-2,:) = vz1(:,-1,:)
    END IF

    IF (proc_y_max == MPI_PROC_NULL .AND. ybc_max == BC_USER) THEN
      vx1(:,ny:ny+2,:) = 0.0_num
      vy1(:,ny:ny+2,:) = 0.0_num
      vz1(:,ny,:) = vz1(:,ny-1,:)
      vz1(:,ny+1,:) = vz1(:,ny,:)
      vz1(:,ny+2,:) = vz1(:,ny+1,:)
    END IF


    IF (proc_z_max == MPI_PROC_NULL .AND. zbc_max == BC_USER) THEN
      vx1(:,:,nz:nz+2) = 0.0_num
      vy1(:,:,nz:nz+2) = 0.0_num
      do i = -2, nx+2
          do j = -2, ny+2
          vz1(i, j, nz) = max(0.0,1.25*vz1(i,j,nz-1))
          end do
      end do
      vz1(:,:,nz+1) = vz1(:,:,nz)
    END IF


    IF (proc_z_min == MPI_PROC_NULL .AND. zbc_min == BC_USER) THEN
      vx1(:,:,-1) = vx1(:,:,0)
      vx1(:,:,-2) = vx1(:,:,-1)
      vy1(:,:,-1) = vy1(:,:,0)
      vy1(:,:,-2) = vy1(:,:,-1)
      vz1(:,:,-2:0) = 0.0_num
    END IF

  END SUBROUTINE remap_v_bcs



END MODULE boundary
