
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
  do j=1,m
    print*, data(:,j,k)
  enddo
enddo
end subroutine printarray

logical function isLastRow(row, blocks)
implicit none
integer, intent(in) :: row
integer, intent(in) :: blocks(2)
isLastRow = (row == blocks(2))
end function isLastRow

logical function isLastCol(col, blocks)
implicit none
integer, intent(in) :: col
integer, intent(in) :: blocks(2)
isLastCol = (col == blocks(1))
end function isLastCol

integer function typeIdx(col, row, blocks)
implicit none
integer, intent(in) :: row,col,blocks(2) 
integer :: lastrow
integer :: lastcol
lastrow = 0
lastcol = 0
if (col == blocks(1)) lastcol = 1
if (row == blocks(2)) lastrow = 1
typeIdx = lastrow*2 + lastcol
end function typeIdx


end module matrix


!-------------------------------------------------------------------

module comm

use mpi
use matrix

implicit none

integer :: comm_cart
integer, allocatable :: all_coord(:)
integer, parameter :: W=1, N=2, E=3, S=4
integer :: neighBor(4)
integer :: size, rank

contains

subroutine updateBounds(myrow, mycol, rank, size, blocks, blocksize,   &
                    globalsizes, localsizes, globaldata, localdata, nguard)


implicit none
integer,   intent(in)   :: myrow, mycol, rank, size, blocksize, nguard
integer,   intent(in)   :: blocks(2), globalsizes(3), localsizes(3)
character, intent(in)   :: globaldata(:,:,:)
character, intent(out)  :: localdata(:,:,:)

integer :: proc, ierr, facetype(2), sendtag, recvtag
integer :: sizes(3), subsizes(3), starts(3)
integer :: a=1, b=2
integer status(MPI_STATUS_SIZE)

do proc=0,size
  if(proc==rank) then
    print*, rank, neighbor
  else
    call MPI_BARRIER(comm_cart, ierr)
  endif
enddo

starts=(/0,0,0/)

sizes = (/localsizes(1)+2*nguard,localsizes(2)+2*nguard,localsizes(3)/)
subsizes = (/ nguard, localsizes(2)+2*nguard, localsizes(3)/) ! xface

call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_CHAR, facetype(1), ierr)
call MPI_TYPE_COMMIT(facetype(1),ierr)

sizes = (/localsizes(1)+2*nguard,localsizes(2)+2*nguard,localsizes(3)/)
subsizes = (/localsizes(1)+2*nguard, nguard, localsizes(3)/) ! yface

call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_CHAR, facetype(2), ierr)
call MPI_TYPE_COMMIT(facetype(2),ierr)

sendtag=1
recvtag=1


! Send my boundary to North and receive from South
call MPI_Sendrecv(localdata(1,nguard+1,1),                 1, facetype(2), neighbor(N), sendtag,       &
                  localdata(1,nguard+localsizes(2)+1,1),   1, facetype(2), neighbor(S), recvtag,       &
                  comm_cart, status, ierr)

! Send my boundary to South and receive from North
call MPI_Sendrecv(localdata(1,localsizes(2)+1,1),   1, facetype(2), neighbor(S), sendtag,       &
                  localdata(1,1,1),                        1, facetype(2), neighbor(N), recvtag,       &
                  comm_cart, status, ierr)

sendtag=2
recvtag=2

! Send my boundary to East and receive from West
call MPI_Sendrecv(localdata(localsizes(1)+1,1,1),        1, facetype(1), neighbor(E), sendtag,       &
                  localdata(1,1,1),                      1, facetype(1), neighbor(W), recvtag,              &
                  comm_cart, status, ierr)

! Send my boundary to West and receive from East
call MPI_Sendrecv(localdata(1+nguard,1,1),                      1, facetype(1), neighbor(W), sendtag,         &
                  localdata(localsizes(1)+nguard+1,1,1), 1, facetype(1), neighbor(E), recvtag,  &
                  comm_cart, status, ierr)

end subroutine updateBounds



subroutine alltoall(myrow, mycol, rank, size, blocks, blocksize,   &
                    globalsizes, localsizes, globaldata, localdata, nguard)

implicit none

integer,   intent(in)   :: myrow, mycol, rank, size, blocksize, nguard
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
integer :: blocktypes(size), subsizes(3), halosizes(3)
integer :: starts(3)

sendcounts = 0
senddispls = 0
sendtypes  = MPI_CHAR

recvcounts = 0
recvdispls = 0
recvtypes  = MPI_CHAR

recvcounts(0) = 1 !localsizes(1) * localsizes(2) * localsizes(3)
recvdispls(0) = 0


halosizes(1) = localsizes(1) + 2*nguard
halosizes(2) = localsizes(2) + 2*nguard
halosizes(3) = localsizes(3) 

starts = (/nguard,nguard,0/)
call MPI_TYPE_CREATE_SUBARRAY(3,halosizes,localsizes,starts,MPI_ORDER_FORTRAN,MPI_CHAR,recvtypes(0),ierr)
call MPI_TYPE_COMMIT(recvtypes(0),ierr)

! The originating process needs to allocate and fill the source array,
! and then define types defining the array chunks to send, and 
! fill out senddispls, sendcounts (1) and sendtypes.

if (rank==0) then
  ! 4 types of blocks - 
  ! blocksize*blocksize, blocksize+1*blocksize, blocksize*blocksize+1, blocksize+1*blocksize+1
  starts      = (/0,0,0/)
  subsizes(3) = localsizes(3)
  do i=0,1
    subsizes(1) = blocksize!+i
    do j=0,1
      subsizes(2) = blocksize!+j
      call MPI_TYPE_CREATE_SUBARRAY(3,globalsizes,subsizes,starts,MPI_ORDER_FORTRAN,MPI_CHAR,blocktypes(2*j+i+1),ierr)
      call MPI_TYPE_COMMIT(blocktypes(2*j+i+1),ierr)
    enddo
  enddo 
  ! now figure out the displacement and type of each processor's data
  do proc=0,size-1
    call rowcol(proc,blocks,row,col)
    col = all_coord(2*proc+1)
    row = all_coord(2*proc+2)
    sendcounts(proc) = 1
    senddispls(proc) = ((row-1)*blocksize*globalsizes(1) + (col-1)*blocksize) * sizeof('.')
    idx = typeIdx(col,row,blocks)
    sendtypes(proc) = blocktypes(idx+1)
  enddo
endif


call MPI_Alltoallw(globaldata, sendcounts, senddispls, sendtypes,     &
                   localdata(1,1,1),  recvcounts, recvdispls, recvtypes, &
                   comm_cart, ierr)




end subroutine alltoall



subroutine collect(myrow, mycol, rank, size, blocks, blocksize,   &
                    globalsizes, localsizes, globaldata, localdata, nguard)

implicit none

integer,   intent(in)   :: myrow, mycol, rank, size, blocksize, nguard
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
integer :: blocktypes(size), subsizes(3), halosizes(3)
integer :: starts(3)

sendcounts = 0
senddispls = 0
sendtypes  = MPI_CHAR

recvcounts = 0
recvdispls = 0
recvtypes  = MPI_CHAR

recvcounts(0) = 1 !localsizes(1) * localsizes(2) * localsizes(3)
recvdispls(0) = 0


halosizes(1) = localsizes(1) + 2*nguard
halosizes(2) = localsizes(2) + 2*nguard
halosizes(3) = localsizes(3) 

starts = (/nguard,nguard,0/)
call MPI_TYPE_CREATE_SUBARRAY(3,halosizes,localsizes,starts,MPI_ORDER_FORTRAN,MPI_CHAR,recvtypes(rank),ierr)
call MPI_TYPE_COMMIT(recvtypes(rank),ierr)

! The originating process needs to allocate and fill the source array,
! and then define types defining the array chunks to send, and 
! fill out senddispls, sendcounts (1) and sendtypes.

if (rank==0) then
  ! 4 types of blocks - 
  ! blocksize*blocksize, blocksize+1*blocksize, blocksize*blocksize+1, blocksize+1*blocksize+1
  starts      = (/0,0,0/)
  subsizes(3) = localsizes(3)
  do i=0,1
    subsizes(1) = blocksize!+i
    do j=0,1
      subsizes(2) = blocksize!+j
      call MPI_TYPE_CREATE_SUBARRAY(3,globalsizes,subsizes,starts,MPI_ORDER_FORTRAN,MPI_CHAR,blocktypes(2*j+i+1),ierr)
      call MPI_TYPE_COMMIT(blocktypes(2*j+i+1),ierr)
    enddo
  enddo 
  ! now figure out the displacement and type of each processor's data
  do proc=0,size-1
    call rowcol(proc,blocks,row,col)
    col = all_coord(2*proc+1)
    row = all_coord(2*proc+2)
    recvcounts(proc) = 0
    recvdispls(proc) = ((row-1)*blocksize*globalsizes(1) + (col-1)*blocksize) * sizeof('.')
    idx = typeIdx(col,row,blocks)
    recvtypes(proc) = blocktypes(idx+1)
  enddo
endif

recvcounts(0) = size

call MPI_Alltoallw(localdata, sendcounts, senddispls, sendtypes,     &
                   globaldata,  recvcounts, recvdispls, recvtypes, &
                   comm_cart, ierr)




end subroutine collect



end module comm


!-------------------------------------------------------------------

program main

use mpi
use matrix
use comm

implicit none

integer   :: ierr
integer   :: blocks(2) = (/0,0/)
integer   :: globalsizes(3), localsizes(3)
character, pointer     :: globaldata(:,:,:)
character, allocatable :: localdata(:,:,:)
integer   :: myrow, mycol
integer   ::  i,j,k,proc
integer, parameter :: LAYERS = 2
integer, parameter :: NGUARD = 2
integer, parameter :: BLOCKSIZE = 5
integer :: dims(2), coord(2)
logical :: period(2), reorder
integer, parameter :: px=10, py=10
character(len=*), parameter :: method = 'alltoall'

call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr ) 
call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr ) 

period(1)=.false.
period(2)=.false.
reorder=.true.

allocate(all_coord(2*px*py))

call MPI_DIMS_CREATE(size, 2, blocks, ierr)

blocks(1) = px
blocks(2) = py

call MPI_CART_CREATE(MPI_COMM_WORLD,2,blocks,period,reorder,comm_cart,ierr)
call MPI_COMM_RANK(comm_cart,rank, ierr)
call MPI_CART_COORDS(comm_cart,rank,2,coord,ierr)

!call rowcol(rank, blocks, myrow, mycol)

coord = coord + 1

call MPI_Allgather(coord, 2, MPI_INTEGER, all_coord, 2, MPI_INTEGER, comm_cart, ierr)

! Left/West and right/Est neigbors
call MPI_Cart_shift(comm_cart,0,1,NeighBor(W),NeighBor(E),ierr)

! Bottom/South and Upper/North neigbors
call MPI_Cart_shift(comm_cart,1,1,NeighBor(N),NeighBor(S),ierr)

mycol = coord(1)
myrow = coord(2)

!create array sizes so that last block has 1 too many rows/cols */
globalsizes(1) = blocks(1)*BLOCKSIZE + 0 
globalsizes(2) = blocks(2)*BLOCKSIZE + 0 
globalsizes(3) = LAYERS 
  
!allocate(globaldata(globalsizes(1),globalsizes(2),globalsizes(3)))
globaldata => NULL ()

if(rank==0) then
  !globaldata = allocchar3darray(globalsizes(1),globalsizes(2),globalsizes(3))
  allocate(globaldata(globalsizes(1),globalsizes(2),globalsizes(3)))
  do k=1,globalsizes(3)
  do j=1,globalsizes(2)
    do i=1,globalsizes(1)
      globaldata(i,j,k) = achar(iachar('a') + mod((j-1)*globalsizes(1)+(k-1)*globalsizes(1)*globalsizes(2)+i-1,26))
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
localsizes = BLOCKSIZE!+2*NGUARD
if (isLastRow(myrow,blocks)) localsizes(1) = localsizes(1) + 0
if (isLastCol(mycol,blocks)) localsizes(2) = localsizes(2) + 0
localsizes(3) = globalsizes(3)
!localdata = allocchar2darray(localsizes(1),localsizes(2))
!localdata = allocchar3darray(localsizes(1),localsizes(2),localsizes(3))
allocate(localdata(localsizes(1)+2*NGUARD,localsizes(2)+2*NGUARD,localsizes(3)))
localdata = '.'
!
!
call alltoall(myrow, mycol, rank, size, blocks, blocksize,   &
              globalsizes, localsizes, globaldata, localdata, nguard)
!
do proc=0,size-1
  if(proc==rank) then
    print*, '--------------'
    print*, 'Rank :', proc
    call printarray(localdata, localsizes(1)+2*NGUARD,localsizes(2)+2*NGUARD,localsizes(3))
  endif
  call MPI_BARRIER(comm_cart, ierr)
enddo
!
call updateBounds(myrow, mycol, rank, size, blocks, blocksize,   &
              globalsizes, localsizes, globaldata, localdata, nguard)
!
do proc=0,size-1
  if(proc==rank) then
    print*, '--------------'
    print*, 'Rank :', proc
    call printarray(localdata, localsizes(1)+2*NGUARD,localsizes(2)+2*NGUARD,localsizes(3))
  endif
  call MPI_BARRIER(comm_cart, ierr)
enddo

if(rank==0) then
  globaldata = '.'
endif

if(rank==0) then
  !globaldata = allocchar3darray(globalsizes(1),globalsizes(2),globalsizes(3))
  print*, 'Global Array'
  print*, globaldata 
  call printarray(globaldata,globalsizes(1),globalsizes(2),globalsizes(3))
endif

call MPI_FINALIZE(ierr);

end program main

!-------------------------------------------------------------------

