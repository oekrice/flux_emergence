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
! Controls all I/O and diagnostics. Output files are 'lare3d.dat',
! 'control.dat', 'en.dat' and a series of snapshots in 'nnnn.sdf'
!******************************************************************************

MODULE diagnostics

  USE shared_data
  USE boundary
  USE conduct
  USE netcdf

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: output_diags, output_snap, energy_correction

CONTAINS

  !****************************************************************************
  ! My diagnostic functions (same as for mf2d)
  !****************************************************************************

  SUBROUTINE output_diags(diag_num)
      !Calculates some diagnostics and saves to netcdf file as for the triangle code, which was fairly neat (if i say so myself...). Should make for easy pythonning.
      IMPLICIT NONE
      INTEGER:: diag_num
      REAL(num):: diag
      character(len=100):: filename
      integer:: ncid, nd_id, id_1!, id_2, id_3, id_4, id_5, id_6, id_7, id_8

      !Allocate diagnostic arrays
      if (diag_num == 0) then
          allocate(diag_time(0:ndiags))
      end if

      !TIME
      diag = time
      diag_time(diag_num) = time

      !ADMIN
      if (run_id < 10) then
      write (filename, "(A19,I1,A3)") "./diagnostics/run00", int(run_id), ".nc"
      else if (run_id < 100) then
          write (filename, "(A18,I2,A3)") "./diagnostics/run0", int(run_id), ".nc"
      else
          write (filename, "(A17,I3,A3)") "./diagnostics/run", int(run_id), ".nc"
      end if

      !Write to diagnostics file, using netcdf
      call try(nf90_create(trim(filename), nf90_clobber, ncid))
      call try(nf90_def_dim(ncid, 'ndiags', ndiags+1, nd_id))  !Make up fake dimensions here

      call try(nf90_def_var(ncid, 'time', nf90_double, (/nd_id/), id_1))
      call try(nf90_enddef(ncid))

      call try(nf90_put_var(ncid, id_1, diag_time))


      call try(nf90_close(ncid))

      diag_num = diag_num + 1

  END SUBROUTINE output_diags


SUBROUTINE output_snap(snap_num)
    !Exports the magnetic field at this plot_num to an appropriate netcdf file
    IMPLICIT NONE

    CHARACTER(LEN =64):: output_filename
    CHARACTER(LEN = 4):: snap_id
    INTEGER:: snap_num, proc_write
    INTEGER:: ncid, vid
    INTEGER:: xs_id, ys_id, zs_id
    INTEGER:: xc_id, yc_id, zc_id
    INTEGER:: bx_id, by_id, bz_id
    INTEGER:: jx_id, jy_id, jz_id
    INTEGER:: en_id, rho_id
    INTEGER:: vx_id, vy_id, vz_id

    REAL(num), DIMENSION(:,:,:):: jx(1:nx,0:ny,0:nz), jy(0:nx,1:ny,0:nz), jz(0:nx,0:ny,1:nz)

    write (snap_id, '(I4.4)') snap_num
    output_filename = trim(trim(data_directory)//trim(snap_id)//'.nc')

    if (rank == 0) then
    call try(nf90_create(trim(output_filename), nf90_clobber, ncid))

    call try(nf90_def_dim(ncid, 'xs', nx_global+1, xs_id))
    call try(nf90_def_dim(ncid, 'ys', ny_global+1, ys_id))
    call try(nf90_def_dim(ncid, 'zs', nz_global+1, zs_id))

    call try(nf90_def_dim(ncid, 'xc', nx_global, xc_id))
    call try(nf90_def_dim(ncid, 'yc', ny_global, yc_id))
    call try(nf90_def_dim(ncid, 'zc', nz_global, zc_id))

    call try(nf90_def_var(ncid, 'bx', nf90_double, (/xs_id ,yc_id, zc_id/), bx_id))
    call try(nf90_def_var(ncid, 'by', nf90_double, (/xc_id ,ys_id, zc_id/), by_id))
    call try(nf90_def_var(ncid, 'bz', nf90_double, (/xc_id ,yc_id, zs_id/), bz_id))

    call try(nf90_def_var(ncid, 'jx', nf90_double, (/xc_id ,ys_id, zs_id/), jx_id))
    call try(nf90_def_var(ncid, 'jy', nf90_double, (/xs_id ,yc_id, zs_id/), jy_id))
    call try(nf90_def_var(ncid, 'jz', nf90_double, (/xs_id ,ys_id, zc_id/), jz_id))

    call try(nf90_def_var(ncid, 'en', nf90_double, (/xc_id ,yc_id, zc_id/), en_id))
    call try(nf90_def_var(ncid, 'rho', nf90_double, (/xc_id ,yc_id, zc_id/), rho_id))

    call try(nf90_def_var(ncid, 'vx', nf90_double, (/xs_id ,ys_id, zs_id/), vx_id))
    call try(nf90_def_var(ncid, 'vy', nf90_double, (/xs_id ,ys_id, zs_id/), vy_id))
    call try(nf90_def_var(ncid, 'vz', nf90_double, (/xs_id ,ys_id, zs_id/), vz_id))

    call try(nf90_enddef(ncid))
    call try(nf90_close(ncid))

    end if
    call MPI_BARRIER(comm,errcode)

    !Calculate current density (as this isn't done explicitly anywhere else)
    jx(1:nx,0:ny,0:nz) = (bz(1:nx,1:ny+1,0:nz) - bz(1:nx,0:ny,0:nz))/ dyc(0) - (by(1:nx,0:ny,1:nz+1) - by(1:nx,0:ny,0:nz))/ dzc(0)

    jy(0:nx,1:ny,0:nz) = (bx(0:nx,1:ny,1:nz+1) - bx(0:nx,1:ny,0:nz)) / dzc(0) - (bz(1:nx+1,1:ny,0:nz) - bz(0:nx,1:ny,0:nz)) / dxc(0)

    jz(0:nx,0:ny,1:nz) = (by(1:nx+1,0:ny,1:nz) - by(0:nx,0:ny,1:nz)) / dxc(0) - (bx(0:nx,1:ny+1,1:nz) - bx(0:nx,0:ny,1:nz)) / dyc(0)

    !Each process writes data in turn
    do proc_write = 0 ,nproc-1
        call MPI_BARRIER(comm,errcode)

        if (rank == proc_write) then
            call try(nf90_open(trim(output_filename), nf90_write, ncid))

            call try(nf90_inq_varid(ncid, 'bx', vid))
            call try(nf90_put_var(ncid, vid, bx(0:nx,1:ny,1:nz), &
            start = (/starts(1)+1,starts(2) + 1, starts(3) + 1/),count = (/nx+1,ny,nz/)))

            call try(nf90_inq_varid(ncid, 'by', vid))
            call try(nf90_put_var(ncid, vid, by(1:nx,0:ny,1:nz), &
            start = (/starts(1)+1,starts(2) + 1, starts(3) + 1/),count = (/nx,ny+1,nz/)))

            call try(nf90_inq_varid(ncid, 'bz', vid))
            call try(nf90_put_var(ncid, vid, bz(1:nx,1:ny,0:nz), &
                 start = (/starts(1)+1,starts(2) + 1, starts(3) + 1/),count = (/nx,ny,nz+1/)))

            call try(nf90_inq_varid(ncid, 'jx', vid))
            call try(nf90_put_var(ncid, vid, jx(1:nx,0:ny,0:nz), &
            start = (/starts(1)+1,starts(2) + 1, starts(3) + 1/),count = (/nx,ny+1,nz+1/)))

            call try(nf90_inq_varid(ncid, 'jy', vid))
            call try(nf90_put_var(ncid, vid, jy(0:nx,1:ny,0:nz), &
            start = (/starts(1)+1,starts(2) + 1, starts(3) + 1/),count = (/nx+1,ny,nz+1/)))

            call try(nf90_inq_varid(ncid, 'jz', vid))
            call try(nf90_put_var(ncid, vid, jz(0:nx,0:ny,1:nz), &
                 start = (/starts(1)+1,starts(2) + 1, starts(3) + 1/),count = (/nx+1,ny+1,nz/)))

            call try(nf90_inq_varid(ncid, 'en', vid))
            call try(nf90_put_var(ncid, vid, energy(1:nx,1:ny,1:nz), &
                 start = (/starts(1)+1,starts(2) + 1, starts(3) + 1/),count = (/nx,ny,nz/)))

            call try(nf90_inq_varid(ncid, 'rho', vid))
            call try(nf90_put_var(ncid, vid, rho(1:nx,1:ny,1:nz), &
            start = (/starts(1)+1,starts(2) + 1, starts(3) + 1/),count = (/nx,ny,nz/)))

            call try(nf90_inq_varid(ncid, 'vx', vid))
            call try(nf90_put_var(ncid, vid, vx(0:nx,0:ny,0:nz), &
            start = (/starts(1)+1,starts(2) + 1, starts(3) + 1/),count = (/nx+1,ny+1,nz+1/)))

            call try(nf90_inq_varid(ncid, 'vy', vid))
            call try(nf90_put_var(ncid, vid, vy(0:nx,0:ny,0:nz), &
            start = (/starts(1)+1,starts(2) + 1, starts(3) + 1/),count = (/nx+1,ny+1,nz+1/)))

            call try(nf90_inq_varid(ncid, 'vz', vid))
            call try(nf90_put_var(ncid, vid, vz(0:nx,0:ny,0:nz), &
            start = (/starts(1)+1,starts(2) + 1, starts(3) + 1/),count = (/nx+1,ny+1,nz+1/)))

            call try(nf90_close(ncid))

        end if
        call MPI_BARRIER(comm,errcode)

    end do

    call mpi_barrier(comm, errcode)
    if (rank == 0) print*, 'Saved snapshot number', snap_num, ' at time', time, 'to file ', output_filename

    snap_num = snap_num + 1
    return


END SUBROUTINE output_snap


SUBROUTINE try(status)
  ! Catch error in reading netcdf field.
  INTEGER, INTENT(IN):: status

  if (status /= NF90_noerr) THEN
      PRINT*,TRIM(ADJUSTL(NF90_STRERROR(status)))
  end if

END SUBROUTINE try



  SUBROUTINE energy_correction

    delta_ke = -delta_ke
    WHERE (delta_ke < 0.0_num) delta_ke = 0.0_num
    delta_ke(:,:,:) = delta_ke(:,:,:) / (rho(:,:,:) * cv(:,:,:))

    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx
          energy(ix,iy,iz) = energy(ix,iy,iz) + delta_ke(ix,iy,iz)
        END DO
      END DO
    END DO

    CALL energy_bcs

  END SUBROUTINE energy_correction

END MODULE diagnostics
