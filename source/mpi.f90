module mpi
 implicit none
 include 'mpif.h'

 ! MPI stuff: number of processors, rank of this processor, and error
 ! code.
 integer, parameter :: px=4
 integer, parameter :: py=4
 
 integer :: nproc, my_rank

contains

 subroutine init_mpi
  implicit none
  integer :: i,ierr
  ! Initialize MPI, learn local rank and total number of processors.
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

  if(nproc.ne.px*py) then
    if(my_rank==0) print*, 'ERROR: px*py .ne. nproc'
    call MPI_ABORT(MPI_COMM_WORLD)
  endif

 end subroutine init_mpi


 subroutine finalize_mpi
  implicit none
  integer :: ierr
  ! MPI library must be shut down.
  call MPI_Finalize(ierr)
 end subroutine finalize_mpi

end module mpi
