#include "include.h"

#ifdef DO_TIME_AVERAGE

module timeav

 implicit none

 integer, parameter :: initial_timeav_count=0
 integer :: timeav_count=initial_timeav_count
 
 contains
 
  subroutine timeav_iteratation()
   use variables
   use overload
   use grid_operate

   implicit none
   
   integer :: k
   
   timeav_count=timeav_count+1
   
   do k=1,nz
    s_h_h(k)=s_h_h(k)%bz+h(k)%bz
   end do
   do k=1,nz
    s_h_u(k)=s_h_u(k)%bz+h_u(k)%bz
   end do
   do k=1,nz
    s_h_v(k)=s_h_v(k)%bz+h_v(k)%bz
   end do
   do k=1,nz
    s_h_z(k)=s_h_z(k)%bz+h_z(k)%bz
   end do
   
   do k=1,nz
    s_hu_h(k)=s_hu_h(k)%bz+h(k)%bz*Ax(u(k))
   end do
   do k=1,nz
    s_hu_u(k)=s_hu_u(k)%bz+h_u(k)%bz*u(k)%bz
   end do
   do k=1,nz
    s_hu_v(k)=s_hu_v(k)%bz+h_v(k)%bz*Ax(Ay(u(k)))
   end do
   do k=1,nz
    s_hu_z(k)=s_hu_z(k)%bz+h_z(k)%bz*Ay(u(k))
   end do
   
   do k=1,nz
    s_hv_h(k)=s_hv_h(k)%bz+h(k)%bz*Ay(v(k))
   end do
   do k=1,nz
    s_hv_u(k)=s_hv_u(k)%bz+h_u(k)%bz*Ay(Ax(v(k)))
   end do
   do k=1,nz
    s_hv_v(k)=s_hv_v(k)%bz+h_v(k)%bz*v(k)%bz
   end do
   do k=1,nz
    s_hv_z(k)=s_hv_z(k)%bz+h_z(k)%bz*Ax(v(k))
   end do
   
   do k=1,nz
    s_huu_h(k)=s_huu_h(k)%bz+h(k)%bz*Ax(u(k)%bz**2)
   end do
   
   do k=1,nz
    s_hvv_h(k)=s_hvv_h(k)%bz+h(k)%bz*Ay(v(k)%bz**2)
   end do
   
   do k=1,nz
    s_huv_h(k)=s_huv_h(k)%bz+h(k)%bz*Ax(u(k))*Ay(v(k))
   end do
   
   do k=1,nz
    s_huv_z(k)=s_huv_z(k)%bz+h_z(k)%bz*Ay(u(k))*Ax(v(k))
   end do
   
   do k=1,nz
    s_huv_uz(k)=s_huv_uz(k)%bz+Ay(hu(k))*Ax(v(k))
   end do
   
   do k=1,nz
    s_huv_vz(k)=s_huv_vz(k)%bz+Ay(hu(k))*Ax(v(k))
   end do
   
   do k=1,nz
    s_hm_x_u(k)=s_hm_x_u(k)%bz+h_u(k)%bz*(-Gx(pres)+Gx(m(k)))
   end do
   
   do k=1,nz
    s_hm_y_v(k)=s_hm_y_v(k)%bz+h_v(k)%bz*(-Gy(pres)+Gy(m(k)))
   end do
   
   do k=1,nz
    s_m_h(k)=s_m_h(k)%bz-pres%bz+m(k)%bz
   end do
   
   do k=1,nz
    s_hm_h(k)=s_hm_h(k)%bz+h(k)%bz*(-pres%bz+m(k)%bz)
   end do
 
   s_utau_u=s_utau_u%bz+utau%bz
   s_vtau_v=s_vtau_v%bz+vtau%bz
 
   s_bfricu_u=s_bfricu_u%bz+bfricu%bz
   s_bfricv_v=s_bfricv_v%bz+bfricv%bz
   
   do k=1,nz
    s_hsmagu_u(k)=s_hsmagu_u(k)%bz+h_u(k)%bz*smagu(k)%bz
   end do
   
   do k=1,nz
    s_hsmagv_v(k)=s_hsmagv_v(k)%bz+h_v(k)%bz*smagv(k)%bz
   end do
   
   do k=1,nz
    s_hq_h(k)=s_hq_h(k)%bz+Ax(Ay(f%bz+zeta(k)%bz))
   end do
   do k=1,nz
    s_hq_u(k)=s_hq_u(k)%bz+Ay(f%bz+zeta(k)%bz)
   end do
   do k=1,nz
    s_hq_v(k)=s_hq_v(k)%bz+Ax(f%bz+zeta(k)%bz)
   end do
   do k=1,nz
    s_hq_z(k)=s_hq_z(k)%bz+f%bz+zeta(k)%bz
   end do
   
   do k=1,nz
    s_huq_h(k)=s_huq_h(k)%bz+Ax(Ay(f%bz+zeta(k)%bz)*u(k)%bz)
   end do
   do k=1,nz
    s_huq_u(k)=s_huq_u(k)%bz+Ay(f%bz+zeta(k)%bz)*u(k)%bz
   end do
   do k=1,nz
    s_huq_v(k)=s_huq_v(k)%bz+Ax((f%bz+zeta(k)%bz)*Ay(u(k)))
   end do
   do k=1,nz
    s_huq_z(k)=s_huq_z(k)%bz+(f%bz+zeta(k)%bz)*Ay(u(k))
   end do
   
   do k=1,nz
    s_hvq_h(k)=s_hvq_h(k)%bz+Ay(Ax(f%bz+zeta(k)%bz)*v(k)%bz)
   end do
   do k=1,nz
    s_hvq_u(k)=s_hvq_u(k)%bz+Ay((f%bz+zeta(k)%bz)*Ax(v(k)))
   end do
   do k=1,nz
    s_hvq_v(k)=s_hvq_v(k)%bz+Ax(f%bz+zeta(k)%bz)*v(k)%bz
   end do
   do k=1,nz
    s_hvq_z(k)=s_hvq_z(k)%bz+(f%bz+zeta(k)%bz)*Ax(v(k))
   end do
 
   do k=1,nz
    s_hqq_z(k)=s_hqq_z(k)%bz+(f%bz+zeta(k)%bz)**2*merge(0.0d0,1.0d0/h_z(k)%bz,(h_z(k)%bz == 0.0d0))
   end do
   
   
  end subroutine
  
  subroutine write_timemean(dat)
   use global
   use overload
   use writeout
   
   implicit none
   
   class(var), intent(inout) :: dat
   
   
   dat=dat%bz/int2real(timeav_count)
   call fullwrite_var('timemean',dat)
   
  end subroutine
  
  subroutine write_timemeans()
   use variables
   use parallel
   
   implicit none
   
   integer :: k
   
   if (proc_name == ens_master) then
    write(format,'(a32)') '(a17,i10,a11,i4)'
    write(*,format) 'Time Mean Count: ', timeav_count, ' ens_name: ', ens_name
   end if 
   
   
   do k=1,nz
    call write_timemean(s_h_h(k))
   end do
   do k=1,nz
    call write_timemean(s_h_u(k))
   end do
   do k=1,nz
    call write_timemean(s_h_v(k))
   end do
   do k=1,nz
    call write_timemean(s_h_z(k))
   end do
   
   do k=1,nz
    call write_timemean(s_hu_h(k))
   end do
   do k=1,nz
    call write_timemean(s_hu_u(k))
   end do
   do k=1,nz
    call write_timemean(s_hu_v(k))
   end do
   do k=1,nz
    call write_timemean(s_hu_z(k))
   end do
   
   do k=1,nz
    call write_timemean(s_hv_h(k))
   end do
   do k=1,nz
    call write_timemean(s_hv_u(k))
   end do
   do k=1,nz
    call write_timemean(s_hv_v(k))
   end do
   do k=1,nz
    call write_timemean(s_hv_z(k))
   end do
   
   do k=1,nz
    call write_timemean(s_huu_h(k))
   end do
   
   do k=1,nz
    call write_timemean(s_hvv_h(k))
   end do
   
   do k=1,nz
    call write_timemean(s_huv_h(k))
   end do
   
   do k=1,nz
    call write_timemean(s_huv_z(k))
   end do
   
   do k=1,nz
    call write_timemean(s_huv_uz(k))
   end do
   
   do k=1,nz
    call write_timemean(s_huv_vz(k))
   end do
   
   do k=1,nz
    call write_timemean(s_hm_x_u(k))
   end do
   
   do k=1,nz
    call write_timemean(s_hm_y_v(k))
   end do
   
   do k=1,nz
    call write_timemean(s_m_h(k))
   end do
 
   call write_timemean(s_utau_u)
   call write_timemean(s_vtau_v)
 
   call write_timemean(s_bfricu_u)
   call write_timemean(s_bfricv_v)
   
   do k=1,nz
    call write_timemean(s_hsmagu_u(k))
   end do
   
   do k=1,nz
    call write_timemean(s_hsmagv_v(k))
   end do
   
   do k=1,nz
    call write_timemean(s_hq_h(k))
   end do
   do k=1,nz
    call write_timemean(s_hq_u(k))
   end do
   do k=1,nz
    call write_timemean(s_hq_v(k))
   end do
   do k=1,nz
    call write_timemean(s_hq_z(k))
   end do
   
   do k=1,nz
    call write_timemean(s_huq_h(k))
   end do
   do k=1,nz
    call write_timemean(s_huq_u(k))
   end do
   do k=1,nz
    call write_timemean(s_huq_v(k))
   end do
   do k=1,nz
    call write_timemean(s_huq_z(k))
   end do
   
   do k=1,nz
    call write_timemean(s_hvq_h(k))
   end do
   do k=1,nz
    call write_timemean(s_hvq_u(k))
   end do
   do k=1,nz
    call write_timemean(s_hvq_v(k))
   end do
   do k=1,nz
    call write_timemean(s_hvq_z(k))
   end do
   
  end subroutine
   
  
  

end module

#endif