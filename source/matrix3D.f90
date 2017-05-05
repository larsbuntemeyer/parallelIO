
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

function allocchar3darray(n, m, l)
implicit none
integer, intent(in) :: n,m,l
character, allocatable :: allocchar3darray(:,:,:)
allocate(allocchar3darray(n,m,l))
allocchar3darray = '.'
end function allocchar3darray

subroutine printarray(data, n, m, l)
implicit none
character, intent(in) :: data(:,:,:)
integer, intent(in) :: n,m,l
integer :: i,j,k
do k=1,l
  print*, 'Layer:', k
  do i=1,n
    print*, data(i,:,k)
  enddo
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
                    globalsizes, localsizes, globaldata, localdata)

implicit none

integer,   intent(in)   :: myrow, mycol, rank, size, blocksize
integer,   intent(in)   :: blocks(2), globalsizes(3), localsizes(3)
character, intent(in)   :: globaldata(:,:,:)
character, intent(out)  :: localdata(:,:,:)


integer :: sendcounts(0:size-1)
integer :: senddispls(0:size-1)
integer :: sendtypes(0:size-1)
integer :: recvcounts(0:size-1)
integer :: recvdispls(0:size-1)
integer :: recvtypes(0:size-1)

integer :: i,j,proc,ierr,row,col,idx
integer :: blocktypes(size), subsizes(3)
integer :: starts(3)

sendcounts = 0
senddispls = 0
sendtypes  = MPI_CHAR

recvcounts = 0
recvdispls = 0
recvtypes  = MPI_CHAR

recvcounts(0) = localsizes(1) * localsizes(2) * localsizes(3)
recvdispls(0) = 0

! The originating process needs to allocate and fill the source array,
! and then define types defining the array chunks to send, and 
! fill out senddispls, sendcounts (1) and sendtypes.

if (rank==0) then
  ! 4 types of blocks - 
  ! blocksize*blocksize, blocksize+1*blocksize, blocksize*blocksize+1, blocksize+1*blocksize+1
  starts      = (/0,0,0/)
  subsizes(3) = localsizes(3)
  do i=0,1
    subsizes(1) = blocksize+i
    do j=0,1
      subsizes(2) = blocksize+j
      call MPI_TYPE_CREATE_SUBARRAY(3,globalsizes,subsizes,starts,MPI_ORDER_FORTRAN,MPI_CHAR,blocktypes(2*i+j+1),ierr)
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
                   localdata(1,1,1),  recvcounts, recvdispls, recvtypes, &
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
integer   :: globalsizes(3), localsizes(3)
character, pointer     :: globaldata(:,:,:)
character, allocatable :: localdata(:,:,:)
integer   :: myrow, mycol
integer   ::  i,j,k,proc
integer, parameter :: BLOCKSIZE = 3
integer, parameter :: LAYERS = 2
character(len=*), parameter :: method = 'alltoall'

call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr ) 
call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr ) 


call MPI_DIMS_CREATE(size, 2, blocks, ierr)


call rowcol(rank, blocks, myrow, mycol)


!create array sizes so that last block has 1 too many rows/cols */
globalsizes(1) = blocks(1)*BLOCKSIZE + 1
globalsizes(2) = blocks(2)*BLOCKSIZE + 1
globalsizes(3) = LAYERS 
  
!allocate(globaldata(globalsizes(1),globalsizes(2),globalsizes(3)))
globaldata => NULL ()

if(rank==0) then
  !globaldata = allocchar3darray(globalsizes(1),globalsizes(2),globalsizes(3))
  allocate(globaldata(globalsizes(1),globalsizes(2),globalsizes(3)))
  do k=1,globalsizes(3)
  do j=1,globalsizes(2)
    do i=1,globalsizes(1)
      globaldata(i,j,k) = achar(iachar('a') + mod((i-1)*globalsizes(2)+(k-1)*globalsizes(1)*globalsizes(2)+j-1,26))
    enddo
  enddo
  enddo
  print*, 'Global Array'
  print*, globaldata 
  call printarray(globaldata,globalsizes(1),globalsizes(2),globalsizes(3))
else
  allocate(globaldata(1,1,1))
endif

!the local chunk we'll be receiving */
localsizes = BLOCKSIZE
if (isLastRow(myrow,blocks)) localsizes(1) = localsizes(1) + 1
if (isLastCol(mycol,blocks)) localsizes(2) = localsizes(2) + 1
localsizes(3) = globalsizes(3)
!localdata = allocchar2darray(localsizes(1),localsizes(2))
!localdata = allocchar3darray(localsizes(1),localsizes(2),localsizes(3))
allocate(localdata(localsizes(1),localsizes(2),localsizes(3)))
!
!
call alltoall(myrow, mycol, rank, size, blocks, blocksize,   &
              globalsizes, localsizes, globaldata, localdata)
!
!
do proc=0,size-1
  if(proc==rank) then
    print*, 'Rank :', proc
    call printarray(localdata, localsizes(1),localsizes(2),localsizes(3))
  endif
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
enddo


call MPI_FINALIZE(ierr);

end program main

!-------------------------------------------------------------------

