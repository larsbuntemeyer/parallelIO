


#F90 = gfortran
F90 = mpif90
#F90 = mpiifort


vpath %.f90 ../source
vpath %.f   ../source

# compiler flags for intel ifort with netcdf 
#CFLAGS = -g -r8 -I/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.3-parallel-bullxmpi-intel14/include
#CFLAGS = -g -check all -fpe0 -warn -traceback -debug extended -r8 -I/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.3-parallel-bullxmpi-intel14/include
CFLAGS = -g -r8 -I/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.3-parallel-bullxmpi-intel14/include
#CFLAGS = -g -fbounds-check -fbacktrace -Wall -fdefault-real-8 -fdefault-integer-8 -I/usr/include
#CFLAGS = -g -fbounds-check -fbacktrace -Wall -I/usr/include

# linker flags for netcdf
LFLAGS = -L/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.3-parallel-bullxmpi-intel14/lib -lnetcdff -L/sw/rhel6-x64/netcdf/parallel_netcdf-1.6.1-bullxmpi-intel14/lib -L/sw/rhel6-x64/netcdf/netcdf_c-4.4.0-parallel-bullxmpi-intel14/lib -Wl,-rpath,/sw/rhel6-x64/netcdf/netcdf_c-4.4.0-parallel-bullxmpi-intel14/lib -lnetcdf -L/sw/rhel6-x64/hdf5/hdf5-1.8.16-parallel-bullxmpi-intel14/lib -Wl,-rpath,/sw/rhel6-x64/hdf5/hdf5-1.8.16-parallel-bullxmpi-intel14/lib -lhdf5 -lhdf5_hl -L/sw/rhel6-x64/sys/libaec-0.3.2-gcc48/lib -Wl,-rpath,/sw/rhel6-x64/sys/libaec-0.3.2-gcc48/lib -lsz -lz -lcurl -Wl,-rpath,/sw/rhel6-x64/netcdf/parallel_netcdf-1.6.1-bullxmpi-intel14/lib -lnetcdf 
#LFLAGS = #-L/usr/lib -lnetcdff -lnetcdf 

