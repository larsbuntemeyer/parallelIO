module grid

use mpi 

implicit none

! local global resolution
integer, parameter :: NDIMS=2
integer, parameter :: NX_GLOBAL=100
integer, parameter :: NY_GLOBAL=100
integer, parameter :: NZ_GLOBAL=100
integer, parameter :: NLAYER=10

! local domain resolution
integer nx,ny,nz

! local domain position
integer ib_global,jb_global,kb_global

! processor decomposition
integer, parameter :: px=10
integer, parameter :: py=10
integer, parameter :: pz=1


contains

subroutine init_grid
 implicit none
 integer :: p
 nx = nx_global/px
 ny = ny_global/py
 nz = nz_global/pz
 if (my_rank==0) write(*,*) 'init_grid - subdomain size: ',nx,ny
 kb_global = 1
 jb_global = (my_rank/px) * ny + 1
 ib_global = (my_rank - (my_rank/px)*px) * nx + 1
 write(*,*) 'my_rank: ', my_rank, 'ib, jb (global): ',ib_global,jb_global
end subroutine init_grid

end module grid
