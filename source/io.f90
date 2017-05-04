module io

 use netcdf
 use mpi
 use dbase
 use grid!, only : NDIMS,NX_GLOBAL,NY_GLOBAL,ib_global,jb_global
 implicit none

 ! This is the name of the data file we will create.
 character (len = *), parameter :: FILE_NAME = "simple_xy_par.nc"
 ! When we create netCDF files, variables and dimensions, we get back
 ! an ID for each one.
 integer :: ncid, varid(NLAYER), dimids(NDIMS)
 integer :: x_dimid, y_dimid, z_dimid
 

contains

 subroutine init_file(parallel)
   implicit none
   logical, intent(in) :: parallel
   integer :: k
   character(len=4) :: layer
   ! Create the netCDF file. The NF90_NETCDF4 flag causes a
   ! HDF5/netCDF-4 file to be created. The comm and info parameters
   ! cause parallel I/O to be enabled. Use either NF90_MPIIO or
   ! NF90_MPIPOSIX to select between MPI/IO and MPI/POSIX.
   if(parallel) then
     call check(nf90_create(FILE_NAME, IOR(NF90_NETCDF4, NF90_MPIIO), ncid, &
          comm = MPI_COMM_WORLD, info = MPI_INFO_NULL))
   else
     call check(nf90_create(FILE_NAME, NF90_NETCDF4, ncid))
   endif
   ! Define the dimensions. NetCDF will hand back an ID for
   ! each. Metadata operations must take place on all processors.
   call check(nf90_def_dim(ncid, "x", NX_GLOBAL, x_dimid))
   call check(nf90_def_dim(ncid, "y", NY_GLOBAL, y_dimid))
   ! The dimids array is used to pass the IDs of the dimensions of
   ! the variables. Note that in fortran arrays are stored in
   ! column-major format.
   dimids = (/ x_dimid, y_dimid /)
   ! Define the variable. The type of the variable in this case is
   ! NF90_INT (4-byte integer).
   do k=1,NLAYER
     write(layer,"(I4.4)") k
     call check(nf90_def_var(ncid, "field_"//layer, NF90_DOUBLE, dimids, varid(k)))
   enddo
   ! End define mode. This tells netCDF we are done defining
   ! metadata. This operation is collective and all processors will
   ! write their metadata to disk.
   call check(nf90_enddef(ncid))
 end subroutine init_file

 subroutine write_field(data_out,nx,ny,k,parallel)
  implicit none
  integer, intent(in) :: nx,ny,k
  real, intent(in) :: data_out(nx,ny)
  logical, intent(in) :: parallel
  ! These will tell where in the data file this processor should
  ! write.
  integer :: start(NDIMS), count(NDIMS)
  ! Write the pretend data to the file. Each processor writes one row.
  if(parallel) then
    start = (/ ib_global, jb_global/)
    count = (/ nx, ny /)
  else
    start = (/ 1,  1/)
    count = (/ nx_global, ny_global /)
  endif
  !write(*,*) 'rank: ', my_rank, 'start: ', start, 'count: ', count
  call check(nf90_put_var(ncid, varid(k), data_out, start = start, &
       count = count))
   
 end subroutine write_field

 subroutine check(status)
   integer, intent ( in) :: status
   
   if(status /= nf90_noerr) then 
     print *, trim(nf90_strerror(status))
     stop 2
   end if
 end subroutine check  

end module io
