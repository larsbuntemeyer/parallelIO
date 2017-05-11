
!-------------------------------------------------------------------

module matrix

implicit none

contains 

subroutine rowcol(rank, blocks, row, col)
implicit none
integer, intent(in)  :: rank, blocks(2)
integer, intent(out) :: row, col
col = rank/blocks(1) + 1
row = mod(rank,blocks(1)) +1
end subroutine rowcol

function allocchar2darray(n, m)
implicit none
integer, intent(in) :: n,m
character, allocatable :: allocchar2darray(:,:)
allocate(allocchar2darray(n,m))
allocchar2darray = '.'
end function allocchar2darray

subroutine printarray(data, n, m)
implicit none
character, intent(in) :: data(:,:)
integer, intent(in) :: n,m
integer :: i,j
do i=1,n
  print*, data(i,:)
enddo
end subroutine printarray

logical function isLastRow(row, blocks)
implicit none
integer, intent(in) :: row
integer, intent(in) :: blocks(2)
isLastRow = (row == blocks(1))
end function isLastRow

logical function isLastCol(col, blocks)
implicit none
integer, intent(in) :: col
integer, intent(in) :: blocks(2)
isLastCol = (col == blocks(2))
end function isLastCol

integer function typeIdx(row, col, blocks)
implicit none
integer, intent(in) :: row,col,blocks(2) 
integer :: lastrow
integer :: lastcol
lastrow = 0
lastcol = 0
if (row == blocks(1)) lastrow = 1
if (col == blocks(2)) lastcol = 1
typeIdx = lastrow*2 + lastcol
end function typeIdx


end module matrix


!-------------------------------------------------------------------

module comm

use mpi
use matrix

implicit none

contains

subroutine alltoall(myrow, mycol, rank, size, blocks, blocksize,   &
                    globalsizes, localsizes, globaldata, localdata, nguard)

implicit none

integer,   intent(in)   :: myrow, mycol, rank, size, blocksize, nguard
integer,   intent(in)   :: blocks(2), globalsizes(2), localsizes(2)
character, intent(in)   :: globaldata(:,:)
character, intent(out)  :: localdata(:,:)


integer :: sendcounts(0:size-1)
integer :: senddispls(0:size-1)
integer :: sendtypes(0:size-1)
integer :: recvcounts(0:size-1)
integer :: recvdispls(0:size-1)
integer :: recvtypes(0:size-1)

integer :: i,j,proc,ierr,row,col,idx
integer :: blocktypes(size), subtypes(size), subsizes(2), halosizes(2)
integer :: starts(2), halo_starts

sendcounts = 0
senddispls = 0
sendtypes  = MPI_CHAR

recvcounts = 0
recvdispls = 0
recvtypes  = MPI_CHAR

recvcounts(0) = 1 !localsizes(1) * localsizes(2)
recvdispls(0) = 0

halosizes(1) = localsizes(1) + 2*nguard
halosizes(2) = localsizes(2) + 2*nguard

starts = (/nguard,nguard/)
call MPI_TYPE_CREATE_SUBARRAY(2,halosizes,localsizes,starts,MPI_ORDER_FORTRAN,MPI_CHAR,recvtypes(0),ierr)
call MPI_TYPE_COMMIT(recvtypes(0),ierr)

if (rank==0) then
  ! 4 types of blocks - 
  ! blocksize*blocksize, blocksize+1*blocksize, blocksize*blocksize+1, blocksize+1*blocksize+1
  starts = (/0,0/)
  do i=0,1
    subsizes(1) = blocksize+i
    do j=0,1
      subsizes(2) = blocksize+j
      call MPI_TYPE_CREATE_SUBARRAY(2,globalsizes,subsizes,starts,MPI_ORDER_FORTRAN,MPI_CHAR,blocktypes(2*i+j+1),ierr)
      call MPI_TYPE_COMMIT(blocktypes(2*i+j+1),ierr)
    enddo
  enddo 
  ! now figure out the displacement and type of each processor's data
  do proc=0,size-1
    call rowcol(proc,blocks,row,col)
    sendcounts(proc) = 1
    senddispls(proc) = ((col-1)*blocksize*globalsizes(1) + (row-1)*blocksize) * sizeof('.')
    idx = typeIdx(row,col,blocks)
    sendtypes(proc) = blocktypes(idx+1)
  enddo
endif


call MPI_Alltoallw(globaldata, sendcounts, senddispls, sendtypes,     &
                   localdata(1,1),  recvcounts, recvdispls, recvtypes, &
                   MPI_COMM_WORLD, ierr)




end subroutine alltoall

end module comm


!-------------------------------------------------------------------

program main

use mpi
use matrix
use comm

implicit none

integer   :: rank, size, ierr
integer   :: blocks(2) = (/0,0/)
integer   :: globalsizes(2), localsizes(2)
character, allocatable :: globaldata(:,:)
character, allocatable :: localdata(:,:)
integer   :: myrow, mycol
integer   ::  i,j,proc
integer, parameter :: BLOCKSIZE = 3
integer, parameter :: NGUARD = 1 
character(len=*), parameter :: method = 'alltoall'

call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr ) 
call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr ) 


call MPI_DIMS_CREATE(size, 2, blocks, ierr)


call rowcol(rank, blocks, myrow, mycol)


!create array sizes so that last block has 1 too many rows/cols */
globalsizes(1) = blocks(1)*BLOCKSIZE + 1
globalsizes(2) = blocks(2)*BLOCKSIZE + 1

if(rank==0) then
  !globaldata = allocchar2darray(globalsizes(1),globalsizes(2))
  allocate(globaldata(globalsizes(1),globalsizes(2)))
  globaldata = '.'
  do j=1,globalsizes(2)
    do i=1,globalsizes(1)
      globaldata(i,j) = achar(iachar('a') + mod((i-1)*globalsizes(2)+j - 1,26))
    enddo
  enddo
  print*, 'Global Array'
  call printarray(globaldata,globalsizes(1),globalsizes(2))
endif

!the local chunk we'll be receiving */
localsizes = BLOCKSIZE
if (isLastRow(myrow,blocks)) localsizes(1) = localsizes(1) + 1
if (isLastCol(mycol,blocks)) localsizes(2) = localsizes(2) + 1
!localdata = allocchar2darray(localsizes(1),localsizes(2))
allocate(localdata(localsizes(1)+2*nguard,localsizes(2)+2*nguard))
localdata = '.'
!
!
call alltoall(myrow, mycol, rank, size, blocks, blocksize,   &
              globalsizes, localsizes, globaldata, localdata, nguard)
!
!
do proc=0,size-1
  if(proc==rank) then
    print*, 'Rank :', proc
    call printarray(localdata, localsizes(1)+2*nguard,localsizes(2)+2*nguard)
  endif
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
enddo


call MPI_FINALIZE(ierr);

end program main

!-------------------------------------------------------------------

