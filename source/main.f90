!     This is part of the netCDF package.
!     Copyright 2006 University Corporation for Atmospheric Research/Unidata.
!     See COPYRIGHT file for conditions of use.

!     This is a very simple example which writes a 2D array of sample
!     data. To handle this in netCDF we create two shared dimensions,
!     "x" and "y", and a netCDF variable, called "data". It uses
!     parallel I/O to write the file from all processors at the same
!     time.

!     This example demonstrates the netCDF Fortran 90 API. This is part
!     of the netCDF tutorial, which can be found at:
!     http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-tutorial
      
!     Full documentation of the netCDF Fortran 90 API can be found at:
!     http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90

!     $Id: simple_xy_par_wr.f90,v 1.3 2010/06/01 15:34:49 ed Exp $

program parallelIO
  use grid
  use io
  use mpi
  use dbase
  implicit none


  ! Loop indexes, and error handling.
  integer :: x, stat, k, cr, cm, c_start, c_end 
  integer :: t_start,t_end, t_max, t_min, total_start
  real :: t_collect, t_write, t_open, t_close, t_total, rate
  logical, parameter :: parallel_write = .true.

  call system_clock(count_rate=cr)
  call system_clock(count_max=cm)
  rate = real(cr)

  call system_clock(c_start)

  t_collect = 0.0
  t_write   = 0.0
  t_open    = 0.0
  t_close   = 0.0
  t_total   = 0.0

 
  call init_mpi

  call init_grid

  call init_data

  !if(parallel_write) then
  !  call system_clock(t_start)
  !  call init_file(.true.)
  !  call system_clock(t_end)
  !  t_open = (t_end-t_start)/rate
  !  call system_clock(t_start)
  !  do k=1,NLAYER
  !    call write_field(field,nx,ny,k,.true.)
  !  enddo
  !  call system_clock(t_end)
  !  t_write = (t_end-t_start)/rate
  !  call system_clock(t_start)
  !  call check( nf90_close(ncid))
  !  call system_clock(t_end)
  !  t_close = (t_end-t_start)/rate
  !elseif(my_rank==0) then
  !  call system_clock(t_start)
  !  call init_file(.false.)
  !  call system_clock(t_end)
  !  t_open = (t_end-t_start)/rate
  !  do k=1,NLAYER 
  !    call system_clock(t_start)
  !    call collect_data
  !    call system_clock(t_end)
  !    t_collect = t_collect + (t_end-t_start)/rate
  !    call system_clock(t_start)
  !    call write_field(globalField,nx_global,ny_global,k,.true.)
  !    call system_clock(t_end)
  !    t_write = t_write + (t_end-t_start)/rate
  !  enddo
  !  call system_clock(t_end)
  !  t_write = (t_end-t_start)/rate
  !  call system_clock(t_start)
  !  call check( nf90_close(ncid))
  !  call system_clock(t_end)
  !  t_close = (t_end-t_start)/rate
  !else
  !  call system_clock(t_start)
  !  do k=1,NLAYER 
  !    call collect_data
  !  enddo
  !  call system_clock(t_end)
  !  t_collect = (t_end-t_start)/rate
  !endif

  call system_clock(c_end)

  t_total=t_open+t_write+t_close+t_collect
  t_total=(c_end-c_start)/rate
  t_min=0.0
  t_max=0.0
  call MPI_REDUCE(t_open, t_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, stat)
  call MPI_REDUCE(t_open, t_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, MPI_COMM_WORLD, stat)
  if (my_rank .eq. 0) print *, "Minimum Opening Time ", t_min
  if (my_rank .eq. 0) print *, "Maximum Opening Time ", t_max
  call MPI_REDUCE(t_write, t_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, stat)
  call MPI_REDUCE(t_write, t_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, MPI_COMM_WORLD, stat)
  if (my_rank .eq. 0) print *, "Minimum Writing Time ", t_min
  if (my_rank .eq. 0) print *, "Maximum Writing Time ", t_max
  call MPI_REDUCE(t_close, t_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, stat)
  call MPI_REDUCE(t_close, t_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, MPI_COMM_WORLD, stat)
  if (my_rank .eq. 0) print *, "Minimum Closing Time ", t_min
  if (my_rank .eq. 0) print *, "Maximum Closing Time ", t_max
  call MPI_REDUCE(t_collect, t_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, stat)
  call MPI_REDUCE(t_collect, t_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, MPI_COMM_WORLD, stat)
  if (my_rank .eq. 0) print *, "Minimum Collect Time ", t_min
  if (my_rank .eq. 0) print *, "Maximum Collect Time ", t_max
  call MPI_REDUCE(t_total, t_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, stat)
  call MPI_REDUCE(t_total, t_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, MPI_COMM_WORLD, stat)
  if (my_rank .eq. 0) print *, "------------------------------"
  if (my_rank .eq. 0) print *, "Minimum Total Time ", t_min
  if (my_rank .eq. 0) print *, "Maximum Total Time ", t_max
  if (my_rank .eq. 0) print *, "Master Opening Time ", t_open
  if (my_rank .eq. 0) print *, "Master Writing Time ", t_write
  if (my_rank .eq. 0) print *, "Master Closing Time ", t_close
  if (my_rank .eq. 0) print *, "Master Collect Time ", t_collect
  if (my_rank .eq. 0) print *, "Master Total Time ", t_max

  if (my_rank .eq. 0) print *, "*** SUCCESS writing example file ", FILE_NAME, "! "

  call finalize_mpi

end program parallelIO
