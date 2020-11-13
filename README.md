# parallelIO

A fortran environment to write netcdf files in parallel and serial mode using NetCDF4 and MPI.

:construction: This code is under construction and I am aware that it needs much more documentation!

## test framework

This code is a test framework for different MPI communication strategies for a typical 2d regular,
cartesian subdomain decomposition. The framework offerst two strategies:

* serial IO with different MPI collection implementations
* parallel IO with NetCDF4

I have used this framework to test and understand different strategies to implement these into the
REMO IO modules.

## build instructions

## contact

lars.buntemeyer@hzg.de
