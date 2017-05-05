module mpi
 implicit none
 include 'mpif.h'

 ! MPI stuff: number of processors, rank of this processor, and error
 ! code.
 integer :: nproc, my_rank, ierr, dims(2), coord(2), vu, err
 logical :: period(2), reorder
 integer, parameter :: px=4
 integer, parameter :: py=4
 integer :: all_coord(px*py,1,1)

contains

 subroutine init_mpi
  implicit none
  integer :: i,ierr
  ! Initialize MPI, learn local rank and total number of processors.
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

  if(nproc.ne.px*py) then
    if(my_rank==0) print*, 'ERROR: np*py .ne. nproc'
    call MPI_ABORT(MPI_COMM_WORLD)
  endif

  dims(1)=px
  dims(2)=py

  period(1)=.false.
  period(2)=.false.
  reorder=.true.

  call MPI_CART_CREATE(MPI_COMM_WORLD,2,dims,period,reorder,vu,err)

  do i=0,nproc-1
   if(i==my_rank) then
     call MPI_CART_COORDS(vu,my_rank,2,coord,err)
     print*,'P:',my_rank,' my coordinates are',coord
   else
     call MPI_BARRIER(vu,err)
   endif
  enddo

  call MPI_Allgather(coord, 2, MPI_INTEGER, all_coord, 2, MPI_INTEGER, &
                     MPI_COMM_WORLD, err)

  if(i==my_rank) then
    do i=0,nproc-1
     call MPI_CART_COORDS(vu,my_rank,2,coord,err)
     print*,'P:',my_rank,' my coordinates are',coord
    else
     call MPI_BARRIER(vu,err)
    enddo
  endif
   
 
 end subroutine init_mpi


 subroutine finalize_mpi
  implicit none
  ! MPI library must be shut down.
  call MPI_Finalize(ierr)
 end subroutine finalize_mpi

end module mpi
