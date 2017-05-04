module mpi
 implicit none
 include 'mpif.h'

 ! MPI stuff: number of processors, rank of this processor, and error
 ! code.
 integer :: nproc, my_rank, ierr

contains

 subroutine init_mpi
  implicit none
  integer :: ierr
  ! Initialize MPI, learn local rank and total number of processors.
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
 end subroutine init_mpi

 subroutine finalize_mpi
  implicit none
  ! MPI library must be shut down.
  call MPI_Finalize(ierr)
 end subroutine finalize_mpi

end module mpi
