
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
! Main driver routine
!******************************************************************************

PROGRAM lare3d

  USE shared_data
  USE initial_conditions
  USE setup
  USE boundary
  USE diagnostics
  USE lagran
  USE remap
  USE mpi_routines
  USE welcome
  USE normalise
  USE neutral
  USE control

  IMPLICIT NONE

  step = 0

  CALL mpi_minimal_init    ! mpi_routines.f90

  CALL before_control      ! setup.F90
  CALL user_normalisation  ! control.f90
  CALL control_variables   ! control.f90
  CALL set_output_dumps    ! control.f90
  CALL mpi_initialise      ! mpi_routines.f90
  CALL after_control       ! setup.f90

  CALL welcome_message     ! welcome.f90

  CALL setup_neutral       ! neutral.f90
  CALL normalise_transport ! normalise.f90

  CALL grid                      ! setup.f90

  if (rank == 0) print*, 'Grid set up. Attempting initial conditions...'


  IF (IAND(initial, IC_NEW) /= 0) THEN
    CALL set_initial_conditions  ! initial_conditions.f90
  END IF


  CALL set_boundary_conditions   ! boundary.f90

  CALL boundary_conditions       ! boundary.f90

  CALL eta_calc                  ! lagran.f90

  IF (eos_number /= EOS_IDEAL) CALL neutral_fraction ! neutral.f90

  IF (eos_number == EOS_IDEAL .AND. neutral_gas) xi_n = 1.0_num

  IF (cowling_resistivity) CALL perpendicular_resistivity ! neutral.f90

  IF (rank == 0) PRINT*, 'Initial conditions setup OK. Running Code'

  DO
    if (time .ge. t_end*float(diag_num)/float(ndiags)) then   ! Save diagnostic data (more frequently than snapshots)
      !CALL test1
      if (rank == 0) print*, 'Diagnostic number', diag_num, 'at time', time
      CALL output_diags(diag_num)
    end if

    if (time .ge. t_end*float(snap_num)/float(nplots)) then   ! Save diagnostic data (more frequently than snapshots)
      !CALL test1
      !if (rank == 0) print*, 'Snap Number', snap_num, 'at time', time
      CALL output_snap(snap_num)
    end if


    IF ((step >= nsteps .AND. nsteps >= 0) .OR. (time >= t_end)) EXIT
    !print*, 'step', step, time, minval(rho), dt
    step = step + 1
    IF (eos_number /= EOS_IDEAL) CALL neutral_fraction ! neutral.f90
    CALL lagrangian_step             ! lagran.f90
    CALL eulerian_remap(step)        ! remap.f90
    IF (rke) CALL energy_correction  ! diagnostics.f90

    energy = energy +  correction_factor*dt*(energy_reference-energy)  !Energy correction factor for the Newton Cooling

    CALL eta_calc                    ! lagran.f90
  END DO

  IF (rank == 0) PRINT*, 'Code Terminated normally'
  CALL mpi_close                     ! mpi_routines.f90
  CALL MPI_FINALIZE(errcode)

END PROGRAM lare3d

















