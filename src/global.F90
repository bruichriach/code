module global
 use system
 use overloading
 
 type bound
  type(use) :: z
  real(kind=db) :: slip
 end type

 type mpi_sendrecv
  type(use), allocatable :: s(:)
  type(bound), allocatable :: r(:)
  real(kind=db), allocatable :: s_mpi(:)
  real(kind=db), allocatable :: r_mpi(:)
  integer :: s_req, r_req
 end type
  
 type grd
  real (kind = db), allocatable :: x(:,:), y(:,:)
  integer, allocatable :: i(:,:), j(:,:)
  integer :: lx, ly, nx, ny
 end type
 
 type var
  type(grd), pointer :: p
  real(kind=db), pointer :: z(:,:)
  real(kind=db), pointer :: bz(:,:)
  type(mpi_sendrecv), allocatable :: mpi(:)
  integer :: tag
  logical :: synced
 end type
  
 type hvar
  type(var) :: z
  type(grd), pointer :: p
 end type
 
 type uvar
  type(var) :: z
  type(grd), pointer :: p
 end type
 
 type vvar
  type(var) :: z
  type(grd), pointer :: p
 end type
 
 type zvar
  type(var) :: z
  type(grd), pointer :: p
 end type
 
  
 type(grd), target :: hgrid, ugrid, vgrid, zgrid
  
end module


module variables
 use global
 use params

 type(hvar) :: s
 
 type(hvar) :: depth(0:nz)
 type(hvar) :: mindepth(0:nz)
 
 type(hvar) :: h(nz)
 type(uvar) :: u(nz)
 type(vvar) :: v(nz)
 
 type(hvar) :: h_tmp(nz)
 type(uvar) :: u_tmp(nz)
 type(vvar) :: v_tmp(nz)
 
 type(hvar) :: ke(nz)
 type(hvar) :: ape(nz)
 
 type(zvar) :: zeta(nz)
 
 type(uvar) :: lapu(nz)
 type(vvar) :: lapv(nz)
 
 type(hvar) :: smag(nz)
 type(hvar) :: tension(nz)
 type(zvar) :: strain(nz)
 
 type(uvar) :: smagu(nz)
 type(vvar) :: smagv(nz)
 
 type(zvar) :: q(nz)
  
 type(hvar) :: m(nz)
 type(hvar) :: minm(nz)
 
 type(uvar) :: bfricu
 type(vvar) :: bfricv
 
 type(uvar) :: utau
 type(vvar) :: vtau
 
 type(zvar) :: f



end module
