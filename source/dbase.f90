!
!
module dbase

use grid!, only: nx,ny,NX_GLOBAL,NY_GLOBAL
use mo_decomp

implicit none

real, allocatable, dimension(:,:) :: field
real, allocatable, dimension(:,:) :: globalField

contains

subroutine init_data
 implicit none
 integer :: i,j,k
 allocate(field(ie,je))
 !call allocate_field
 do j=1,je
  do i=1,ie
    field(i,j) = real(my_rank)
  enddo
 enddo 

end subroutine init_data

subroutine allocate_field
 implicit none
 allocate(field(ie,je))
end subroutine allocate_field

subroutine collect_data
 implicit none
 real, allocatable, dimension(:,:) :: sendbuf, recvbuf
 integer  :: root, sendcount, recvcount, doublesize, ierr
 integer  :: newtype, resizedtype
 integer  :: col, row, i, j
 integer, dimension(2) :: sizes, subsizes, starts
 integer, dimension(MPI_STATUS_SIZE) :: rstatus
 integer(kind=MPI_ADDRESS_KIND) :: extent, begin
 integer, dimension(px*py) :: counts, displs
 !
 root=0
 if(my_rank==root) then
   if(.not.(allocated(globalField))) then
     allocate(globalField(NX_GLOBAL,NY_GLOBAL))
   endif
   globalField = 0.0
 else
   if(.not.(allocated(globalField))) then
     allocate(globalField(0,0))
   endif
 endif
 !
 sizes    = [NX_GLOBAL,NY_GLOBAL]     ! size of global array
 subsizes = [nx,ny]     ! size of sub-region 
 starts   = [0,0]       ! let's say we're looking at region "0"
                        ! which begins at offset [0,0]
 call MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, newtype, ierr)
 call MPI_Type_size(MPI_DOUBLE_PRECISION, doublesize, ierr)
 extent = nx*doublesize 
 begin = 0
 call MPI_Type_create_resized(newtype, begin, extent, resizedtype, ierr)
 call MPI_Type_commit(resizedtype, ierr)
 !
 counts = 1
 forall( col=1:py, row=1:px )
    displs(1+(row-1)+py*(col-1)) = (row-1) + nx*px*(col-1)
 endforall
 !
 call MPI_Gatherv( field, nx*ny, MPI_DOUBLE_PRECISION, & ! I'm sending localsize**2 chars
                   globalField, counts, displs, resizedtype,&
                   root, MPI_COMM_WORLD, ierr)
 !
 !if(my_rank==root) print*, globalField
 !
end subroutine collect_data

end module dbase

