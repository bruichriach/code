

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
 
 type(hvar), allocatable :: d(:)
 type(hvar), allocatable :: mind(:)
 type(uvar), allocatable :: d_u(:)
 type(vvar), allocatable :: d_v(:)
 
 type(hvar), allocatable :: h(:)
 type(uvar), allocatable :: u(:)
 type(vvar), allocatable :: v(:)
 

#ifdef DO_TIME_AVERAGE
 type(uvar), allocatable :: hu(:)
 type(vvar), allocatable :: hv(:)
#endif
 
 type(hvar), allocatable :: minh(:)
 
 type(uvar), allocatable :: h_u(:)
 type(vvar), allocatable :: h_v(:)
 type(zvar), allocatable :: h_z(:)
 
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
 
 type(zvar), allocatable :: q(:)
  
 type(hvar), allocatable :: m(:)
 type(hvar), allocatable :: minm(:)
 
 type(uvar) :: bfricu
 type(vvar) :: bfricv
 
 type(uvar) :: utau
 type(vvar) :: vtau
 
 type(zvar) :: f

#ifdef DO_TIME_AVERAGE

 type(hvar), allocatable :: s_h_h(:)
 type(uvar), allocatable :: s_h_u(:)
 type(vvar), allocatable :: s_h_v(:)
 type(zvar), allocatable :: s_h_z(:)
 
 type(hvar), allocatable :: s_hu_h(:)
 type(uvar), allocatable :: s_hu_u(:)
 type(vvar), allocatable :: s_hu_v(:)
 type(zvar), allocatable :: s_hu_z(:)

 type(hvar), allocatable :: s_hv_h(:)
 type(uvar), allocatable :: s_hv_u(:)
 type(vvar), allocatable :: s_hv_v(:)
 type(zvar), allocatable :: s_hv_z(:)

 type(hvar), allocatable :: s_u_h(:)
 type(uvar), allocatable :: s_u_u(:)
 type(vvar), allocatable :: s_u_v(:)
 type(zvar), allocatable :: s_u_z(:)

 type(hvar), allocatable :: s_v_h(:)
 type(uvar), allocatable :: s_v_u(:)
 type(vvar), allocatable :: s_v_v(:)
 type(zvar), allocatable :: s_v_z(:)
 
 type(hvar), allocatable :: s_huu_h(:)
 type(zvar), allocatable :: s_huu_z(:)

 type(hvar), allocatable :: s_hvv_h(:)
 type(zvar), allocatable :: s_hvv_z(:)

 type(hvar), allocatable :: s_huv_h(:)
 type(zvar), allocatable :: s_huv_z(:)
 type(zvar), allocatable :: s_huv_uz(:)
 type(zvar), allocatable :: s_huv_vz(:)
 
 type(uvar), allocatable :: s_hm_x_u(:)
 type(vvar), allocatable :: s_hm_y_v(:)
 
 type(uvar), allocatable :: s_up_dm_x_u(:)
 type(vvar), allocatable :: s_up_dm_y_v(:)
 
 type(uvar), allocatable :: s_dn_dm_x_u(:)
 type(vvar), allocatable :: s_dn_dm_y_v(:)
 
 type(hvar), allocatable :: s_dd_h(:)
 
 type(hvar), allocatable :: s_d_h(:)
 
 type(hvar), allocatable :: s_m_h(:)
 
 type(hvar), allocatable :: s_hm_h(:)
 
 type(uvar) :: s_utau_u
 type(vvar) :: s_vtau_v
 
 type(uvar) :: s_bfricu_u
 type(vvar) :: s_bfricv_v
 
 type(uvar), allocatable :: s_hsmagu_u(:)
 type(vvar), allocatable :: s_hsmagv_v(:)

 type(hvar), allocatable :: s_hq_h(:)
 type(uvar), allocatable :: s_hq_u(:)
 type(vvar), allocatable :: s_hq_v(:)
 type(zvar), allocatable :: s_hq_z(:)

 type(hvar), allocatable :: s_huq_h(:)
 type(uvar), allocatable :: s_huq_u(:)
 type(vvar), allocatable :: s_huq_v(:)
 type(zvar), allocatable :: s_huq_z(:)

 type(hvar), allocatable :: s_hvq_h(:)
 type(uvar), allocatable :: s_hvq_u(:)
 type(vvar), allocatable :: s_hvq_v(:)
 type(zvar), allocatable :: s_hvq_z(:)
 
 type(zvar), allocatable :: s_hqq_z(:)

 type(hvar), allocatable :: s_q_h(:)
 type(uvar), allocatable :: s_q_u(:)
 type(vvar), allocatable :: s_q_v(:)
 type(zvar), allocatable :: s_q_z(:)

 type(hvar), allocatable :: s_uq_h(:)
 type(uvar), allocatable :: s_uq_u(:)
 type(vvar), allocatable :: s_uq_v(:)
 type(zvar), allocatable :: s_uq_z(:)

 type(hvar), allocatable :: s_vq_h(:)
 type(uvar), allocatable :: s_vq_u(:)
 type(vvar), allocatable :: s_vq_v(:)
 type(zvar), allocatable :: s_vq_z(:)
 
 type(zvar), allocatable :: s_qq_z(:)

 type(hvar), allocatable :: s_tendh_h(:)

 type(uvar), allocatable :: s_htendu_u(:)

 type(vvar), allocatable :: s_htendv_v(:)

 type(hvar), allocatable :: s_utendh_h(:)

 type(hvar), allocatable :: s_vtendh_h(:)

#endif

#ifdef ALLOW_RIGID_LID

 type(hvar) :: pres, p, r
 type(uvar) :: thavx,inthavx
 type(vvar) :: thavy,inthavy
 type(hvar) :: inty
 type(hvar), allocatable :: y(:)
 type(hvar) :: ap, z

#endif


end module
