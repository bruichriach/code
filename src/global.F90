

#include "include.h"

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
  character(32) name
 end type out_var
 
 type var
  type(grd), pointer :: p
  real(kind=db), pointer :: z(:,:)
  type(var), pointer :: tmp
  type(var), pointer :: tend_out
  type(var), pointer :: tend1
  type(var), pointer :: tend2
  type(var), pointer :: tend3
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
 
 type(hvar) :: d(0:nz)
 type(hvar) :: mind(0:nz)
 
 type(hvar) :: h(nz)
 type(uvar) :: u(nz)
 type(vvar) :: v(nz)
 

#ifdef DO_TIME_AVERAGE
 type(uvar) :: hu(nz)
 type(vvar) :: hv(nz)
#endif
 
 type(hvar) :: minh(nz)
 
 type(uvar) :: h_u(nz)
 type(vvar) :: h_v(nz)
 type(zvar) :: h_z(nz)
 
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

#ifdef DO_TIME_AVERAGE

 type(hvar) :: s_h_h(nz)
 type(uvar) :: s_h_u(nz)
 type(vvar) :: s_h_v(nz)
 type(zvar) :: s_h_z(nz)
 
 type(hvar) :: s_hu_h(nz)
 type(uvar) :: s_hu_u(nz)
 type(vvar) :: s_hu_v(nz)
 type(zvar) :: s_hu_z(nz)

 type(hvar) :: s_hv_h(nz)
 type(uvar) :: s_hv_u(nz)
 type(vvar) :: s_hv_v(nz)
 type(zvar) :: s_hv_z(nz)
 
 type(hvar) :: s_huu_h(nz)

 type(hvar) :: s_hvv_h(nz)

 type(hvar) :: s_huv_h(nz)
 type(zvar) :: s_huv_z(nz)
 type(zvar) :: s_huv_uz(nz)
 type(zvar) :: s_huv_vz(nz)
 
 type(uvar) :: s_hm_x_u(nz)
 type(vvar) :: s_hm_y_v(nz)
 
 type(uvar) :: s_m_h(nz)
 
 type(uvar) :: s_hm_h(nz)
 
 type(uvar) :: s_utau_u
 type(vvar) :: s_vtau_v
 
 type(uvar) :: s_bfricu_u
 type(vvar) :: s_bfricv_v
 
 type(uvar) :: s_hsmagu_u(nz)
 type(vvar) :: s_hsmagv_v(nz)

 type(hvar) :: s_hq_h(nz)
 type(uvar) :: s_hq_u(nz)
 type(vvar) :: s_hq_v(nz)
 type(zvar) :: s_hq_z(nz)

 type(hvar) :: s_huq_h(nz)
 type(uvar) :: s_huq_u(nz)
 type(vvar) :: s_huq_v(nz)
 type(zvar) :: s_huq_z(nz)

 type(hvar) :: s_hvq_h(nz)
 type(uvar) :: s_hvq_u(nz)
 type(vvar) :: s_hvq_v(nz)
 type(zvar) :: s_hvq_z(nz)
 
 type(zvar) :: s_hqq_z(nz)

#endif

#ifdef ALLOW_RIGID_LID

 type(hvar) :: pres, p, r
 type(uvar) :: thavx,inthavx
 type(vvar) :: thavy,inthavy
 type(hvar) :: inty
 type(hvar) :: y(nz)
 type(hvar) :: ap, z

#endif


end module
