module global
 use sys
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
  integer :: lx, ly, nx, ny, ox, oy
 end type
 
 type send_dat
  real(kind=db), allocatable :: dat(:)
 end type send_dat
 
 type out_var
  integer :: tag
  integer :: req
  integer, allocatable :: req_master(:)
  real(kind=db), allocatable :: send(:)
  type(send_dat), allocatable :: recv(:)
  real(kind=db), allocatable :: z(:,:)
  integer, pointer :: grid(:,:,:)
  character(16) name
 end type out_var
 
 type var
  type(grd), pointer :: p
  real(kind=db), pointer :: z(:,:)
  real(kind=db), pointer :: bz(:,:)
  type(out_var) :: out
  type(mpi_sendrecv), allocatable :: mpi(:)
  integer :: tag
  logical :: synced
 end type
  
 type, extends(var) :: hvar
 end type
 
 type, extends(var) :: uvar
 end type
 
 type, extends(var) :: vvar
 end type
 
 type, extends(var) :: zvar
 end type
 
  
 type(grd), target :: hgrid, ugrid, vgrid, zgrid
 
 contains
 
  subroutine print_var(dat)
  
   implicit none
   
   class(var), intent(in) :: dat
   integer :: i,j
   
   print *
   do j=dat%p%ly,dat%p%ly+dat%p%ny+1
    print *,  (dat%z(i,j) , i=dat%p%lx,dat%p%lx+dat%p%nx+1)
   end do
   print *
  
  end subroutine
  
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
 
 type(hvar), pointer :: tendh0(:),tendh1(:),tendh2(:),tendh3(:)
 type(uvar), pointer :: tendu0(:),tendu1(:),tendu2(:),tendu3(:)
 type(vvar), pointer :: tendv0(:),tendv1(:),tendv2(:),tendv3(:)
 
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


module solver_variables
 use global
 use params

 implicit none

 type(hvar) :: pres, p, r
 type(uvar) :: thavx,inthavx
 type(vvar) :: thavy,inthavy
 type(hvar) :: inty
 type(hvar) :: y(nz)
 type(hvar) :: ap, z

end module
