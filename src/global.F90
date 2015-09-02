

#include "include.h"

module global
 use params
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
#ifdef DO_TIME_AVERAGE
  type(var), pointer :: t
#endif
#ifdef DO_SHORT_AVERAGE
  type(var), pointer :: s
  integer, pointer :: counter
#endif
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
 
 
#ifdef DO_SHORT_AVERAGE
 integer, target :: counter=0
#endif
 
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
 
 type(hvar), allocatable :: d(:)
 type(hvar), allocatable :: mind(:)
 type(uvar), allocatable :: d_u(:)
 type(vvar), allocatable :: d_v(:)
 
 type(hvar), allocatable :: dd(:)
 
 type(hvar), allocatable :: h(:)
 type(uvar), allocatable :: u(:)
 type(vvar), allocatable :: v(:)
 
 type(uvar), allocatable :: h_u(:)
 type(vvar), allocatable :: h_v(:)
 type(zvar), allocatable :: h_z(:)
 
 type(hvar), allocatable :: u_h(:)
 type(vvar), allocatable :: u_v(:)
 type(zvar), allocatable :: u_z(:)
 
 type(hvar), allocatable :: v_h(:)
 type(uvar), allocatable :: v_u(:)
 type(zvar), allocatable :: v_z(:)
 
 type(uvar), allocatable :: hu(:)
 type(vvar), allocatable :: hv(:)
 
 type(hvar), allocatable :: hu_h(:)
 type(vvar), allocatable :: hu_v(:)
 type(zvar), allocatable :: hu_z(:)
 
 type(hvar), allocatable :: hv_h(:)
 type(uvar), allocatable :: hv_u(:)
 type(zvar), allocatable :: hv_z(:)
 
 type(hvar), allocatable :: huu_h(:)
 type(zvar), allocatable :: huu_z(:)

 type(hvar), allocatable :: hvv_h(:)
 type(zvar), allocatable :: hvv_z(:)
 
 type(hvar), allocatable :: huv_h(:)
 type(zvar), allocatable :: huv_z(:)
 type(zvar), allocatable :: huv_uz(:)
 type(zvar), allocatable :: huv_vz(:)
 
 type(hvar), allocatable :: huuu(:)
 type(hvar), allocatable :: huuv(:)
 type(hvar), allocatable :: huvv(:)
 type(hvar), allocatable :: hvvv(:)
 
 type(hvar), allocatable :: minh(:)
 
 type(hvar), allocatable :: h_tmp(:)
 type(uvar), allocatable :: u_tmp(:)
 type(vvar), allocatable :: v_tmp(:)
 
 type(hvar), allocatable :: ke(:)
 type(hvar), allocatable :: ape(:)
 
 type(zvar), allocatable :: zeta(:)
 
 type(uvar), allocatable :: lapu(:)
 type(vvar), allocatable :: lapv(:)
 
 type(hvar), allocatable :: smag(:)
 type(hvar), allocatable :: tension(:)
 type(zvar), allocatable :: strain(:)
 
 type(uvar), allocatable :: smagu(:)
 type(vvar), allocatable :: smagv(:)

 type(uvar), allocatable :: hsmagu(:)
 type(vvar), allocatable :: hsmagv(:)

 type(uvar), allocatable :: husmagu(:)
 type(vvar), allocatable :: hvsmagv(:)
  
 type(hvar), allocatable :: m(:)
 type(hvar), allocatable :: minm(:)
 
 type(uvar), allocatable :: hum(:)
 type(vvar), allocatable :: hvm(:)
 
 type(uvar), allocatable :: hm_x_u(:)
 type(vvar), allocatable :: hm_y_v(:)
 
 type(uvar), allocatable :: up_dm_x_u(:)
 type(vvar), allocatable :: up_dm_y_v(:)
 
 type(uvar), allocatable :: dn_dm_x_u(:)
 type(vvar), allocatable :: dn_dm_y_v(:)
 
 type(uvar), allocatable :: up_du_u(:)
 type(vvar), allocatable :: up_dv_v(:)
 
 type(uvar), allocatable :: dn_du_u(:)
 type(vvar), allocatable :: dn_dv_v(:)
 
 type(hvar), allocatable :: hm(:)
 
 type(uvar) :: bfricu
 type(vvar) :: bfricv
 
 type(uvar) :: utau
 type(vvar) :: vtau
 
 type(uvar) :: ubfricu
 type(vvar) :: vbfricv
 
 type(uvar) :: uutau
 type(vvar) :: vvtau
 
 type(zvar) :: f
 
 type(zvar), allocatable :: q(:)
 
 type(hvar), allocatable :: q_h(:)
 type(uvar), allocatable :: q_u(:)
 type(vvar), allocatable :: q_v(:)
 
 type(zvar), allocatable :: hq(:)

 type(hvar), allocatable :: hq_h(:)
 type(uvar), allocatable :: hq_u(:)
 type(vvar), allocatable :: hq_v(:)
 
 type(uvar), allocatable :: huq(:)

 type(hvar), allocatable :: huq_h(:)
 type(vvar), allocatable :: huq_v(:)
 type(zvar), allocatable :: huq_z(:)
 
 type(vvar), allocatable :: hvq(:)

 type(hvar), allocatable :: hvq_h(:)
 type(uvar), allocatable :: hvq_u(:)
 type(vvar), allocatable :: hvq_v(:)
 type(zvar), allocatable :: hvq_z(:)
 
 type(zvar), allocatable :: hqq(:)
 
 type(uvar), allocatable :: uq(:)

 type(hvar), allocatable :: uq_h(:)
 type(vvar), allocatable :: uq_v(:)
 type(zvar), allocatable :: uq_z(:)
 
 type(vvar), allocatable :: vq(:)

 type(hvar), allocatable :: vq_h(:)
 type(uvar), allocatable :: vq_u(:)
 type(vvar), allocatable :: vq_v(:)
 type(zvar), allocatable :: vq_z(:)
 
 type(zvar), allocatable :: qq(:)
 
 type(hvar), allocatable :: tendh(:)

 type(uvar), allocatable :: htendu(:)

 type(vvar), allocatable :: htendv(:)

 type(hvar), allocatable :: utendh(:)
 
 type(hvar), allocatable :: vtendh(:)
 

#ifdef ALLOW_RIGID_LID

 type(hvar) :: pres, p, r
 type(uvar) :: thavx,inthavx
 type(vvar) :: thavy,inthavy
 type(hvar) :: inty
 type(hvar), allocatable :: y(:)
 type(hvar) :: ap, z

#endif


end module
