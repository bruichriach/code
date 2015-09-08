#include "include.h"

#ifdef DO_TIME_AVERAGE

module timeav

 implicit none
 
 contains

  subroutine av_iteratation(dat)
   use global
   use overload
   
   class(var) :: dat
   
   dat%t=dat%t%bz+dat%bz
#ifdef DO_SHORT_AVERAGE
   dat%s=dat%s%bz+dat%bz
#endif
   
  end subroutine
 
  subroutine timeav_iteratation()
   use variables
   use overload
   use grid_operate
   use sync

   implicit none
   
   integer :: k
   
   timeav_count=timeav_count+1
   counter=counter+1
   
   do k=1,nz
    call av_iteratation(h(k))
   end do
   
   do k=1,nz
    call av_iteratation(hu(k))
   end do
   
   do k=1,nz
    call av_iteratation(hv(k))
   end do

   do k=1,nz
    call av_iteratation(u(k))
   end do

   do k=1,nz
    call av_iteratation(v(k))
   end do
   
   do k=1,nz
    huu(k)=h(k)%bz*Ax(u(k)%bz**2)
    call av_iteratation(huu(k))
   end do
   
   do k=1,nz
    hvv(k)=h(k)%bz*Ay(v(k)%bz**2)
    call av_iteratation(hvv(k))
   end do
   
   do k=1,nz
    huv_h(k)=h(k)%bz*Ax(u(k))*Ay(v(k))
    call av_iteratation(huv_h(k))
   end do
   
   do k=1,nz
    call end_sync(u(k))
    call end_sync(v(k))
    huv_z(k)=h_z(k)%bz*Ay(u(k))*Ax(v(k))
    call av_iteratation(huv_z(k))
   end do
   
   do k=1,nz
    huv_uz(k)=Ay(hu(k))*Ax(v(k))
    call av_iteratation(huv_uz(k))
   end do
   
   do k=1,nz
    huv_vz(k)=Ay(u(k))*Ax(hv(k))
    call av_iteratation(huv_vz(k))
   end do
   
   do k=1,nz
    huuu(k)=h_u(k)%bz*u(k)%bz**3
    call av_iteratation(huuu(k))
   end do

   do k=1,nz
    call end_sync(uu(k))
    huuv(k)=h_v(k)%bz*Ax(Ay(uu(k)))*v(k)%bz
    call av_iteratation(huuv(k))
   end do
   
   do k=1,nz
    call end_sync(vv(k))
    huvv(k)=h_u(k)%bz*u(k)%bz*Ay(Ax(vv(k)))
    call av_iteratation(huvv(k))
   end do

   do k=1,nz
    hvvv(k)=h_v(k)%bz*v(k)%bz**3
    call av_iteratation(hvvv(k))
   end do
   
   do k=1,nz
    call end_sync(m(k))
    hum(k)=h_u(k)%bz*u(k)%bz*Ax(m(k))
    call av_iteratation(hum(k))
   end do

   do k=1,nz
    hvm(k)=h_v(k)%bz*v(k)%bz*Ay(m(k))
    call av_iteratation(hvm(k))
   end do
      
   do k=1,nz
    hm_x(k)=h_u(k)%bz*Gx(m(k))
    call av_iteratation(hm_x(k))
   end do
   
   do k=1,nz
    hm_y(k)=h_v(k)%bz*Gy(m(k))
    call av_iteratation(hm_y(k))
   end do
   
   do k=1,nz
    up_dm_x(k)=d_u(k)%bz*Gx(m(k))
    call av_iteratation(up_dm_x(k))
   end do
   
   do k=1,nz
    up_dm_y(k)=d_v(k)%bz*Gy(m(k))
    call av_iteratation(up_dm_y(k))
   end do
   
   do k=1,nz
    up_du(k)=d_u(k)%bz*u(k)%bz
    call av_iteratation(up_du(k))
   end do
   
   do k=1,nz
    up_dv(k)=d_v(k)%bz*v(k)%bz
    call av_iteratation(up_dv(k))
   end do
   
   do k=1,nz
    dn_dm_x(k)=d_u(k-1)%bz*Gx(m(k))
    call av_iteratation(dn_dm_x(k))
   end do
   
   do k=1,nz
    dn_dm_y(k)=d_v(k-1)%bz*Gy(m(k))
    call av_iteratation(dn_dm_y(k))
   end do
   
   do k=1,nz
    dn_du(k)=d_u(k-1)%bz*u(k)%bz
    call av_iteratation(dn_du(k))
   end do
   
   do k=1,nz
    dn_dv(k)=d_v(k-1)%bz*v(k)%bz
    call av_iteratation(dn_dv(k))
   end do
   
   do k=0,nz
    dd(k)=d(k)%bz**2
    call av_iteratation(dd(k))
   end do
   
   do k=0,nz
    call av_iteratation(d(k))
   end do
      
   do k=1,nz
    call av_iteratation(m(k))
   end do
   
   do k=1,nz
    hm(k)=h(k)%bz*m(k)%bz
    call av_iteratation(hm(k))
   end do
 
   call av_iteratation(utau)
   call av_iteratation(vtau)
 
   call av_iteratation(bfricu)
   call av_iteratation(bfricv)
 
   uutau=u(nz)%bz*utau%bz
   call av_iteratation(uutau)
   vvtau=v(nz)%bz*vtau%bz
   call av_iteratation(vvtau)
 
   ubfricu=u(1)%bz*bfricu%bz
   call av_iteratation(ubfricu)
   vbfricv=v(1)%bz*bfricv%bz
   call av_iteratation(vbfricv)
   
   do k=1,nz
    hsmagu(k)=h_u(k)%bz*smagu(k)%bz
    call av_iteratation(hsmagu(k))
   end do
   
   do k=1,nz
    hsmagv(k)=h_v(k)%bz*smagv(k)%bz
    call av_iteratation(hsmagv(k))
   end do
   
   do k=1,nz
    husmagu(k)=h_u(k)%bz*u(k)%bz*smagu(k)%bz
    call av_iteratation(husmagu(k))
   end do
   
   do k=1,nz
    hvsmagv(k)=h_v(k)%bz*v(k)%bz*smagv(k)%bz
    call av_iteratation(hvsmagv(k))
   end do
   
   do k=1,nz
    hq(k)=f%bz+zeta(k)%bz
    call av_iteratation(hq(k))
   end do
   
   do k=1,nz
    huq(k)=Ay(f%bz+zeta(k)%bz)*u(k)%bz
    call av_iteratation(huq(k))
   end do
   do k=1,nz
    huq_v(k)=Ax((f%bz+zeta(k)%bz)*Ay(u(k)))
    call av_iteratation(huq_v(k))
   end do
   
   do k=1,nz
    hvq_u(k)=Ay((f%bz+zeta(k)%bz)*Ax(v(k)))
    call av_iteratation(hvq_u(k))
   end do
   do k=1,nz
    hvq(k)=Ax(f%bz+zeta(k)%bz)*v(k)%bz
    call av_iteratation(hvq(k))
   end do
 
   do k=1,nz
    hqq(k)=(f%bz+zeta(k)%bz)**2*merge(0.0d0,1.0d0/h_z(k)%bz,(h_z(k)%bz == 0.0d0))
    call av_iteratation(hqq(k))
   end do
   
   do k=1,nz
    q_h(k)=Ax(Ay(f%bz+zeta(k)%bz))*merge(0.0d0,1.0d0/h(k)%bz,(h(k)%bz == 0.0d0))
    call av_iteratation(q_h(k))
   end do
   do k=1,nz
    call av_iteratation(q(k))
   end do
   
   do k=1,nz
    uq(k)=Ay(f%bz+zeta(k)%bz)*u(k)%bz*merge(0.0d0,1.0d0/h_u(k)%bz,(h_u(k)%bz == 0.0d0))
    call av_iteratation(uq(k))
   end do
   do k=1,nz
    uq_v(k)=Ax((f%bz+zeta(k)%bz)*Ay(u(k)))*merge(0.0d0,1.0d0/h_v(k)%bz,(h_v(k)%bz == 0.0d0))
    call av_iteratation(uq_v(k))
   end do
   
   do k=1,nz
    vq_u(k)=Ay((f%bz+zeta(k)%bz)*Ax(v(k)))*merge(0.0d0,1.0d0/h_u(k)%bz,(h_u(k)%bz == 0.0d0))
    call av_iteratation(vq_u(k))
   end do
   do k=1,nz
    vq(k)=Ax(f%bz+zeta(k)%bz)*v(k)%bz*merge(0.0d0,1.0d0/h_v(k)%bz,(h_v(k)%bz == 0.0d0))
    call av_iteratation(vq(k))
   end do
 
   do k=1,nz
    qq(k)=(f%bz+zeta(k)%bz)**2*merge(0.0d0,1.0d0/h_z(k)%bz,(h_z(k)%bz == 0.0d0))**2
    call av_iteratation(qq(k))
   end do

   do k=1,nz
    tendh(k)=h(k)%tend1%bz
    call av_iteratation(tendh(k))
   end do

   do k=1,nz
    htendu(k)=h_u(k)%bz*u(k)%tend1%bz
    call av_iteratation(htendu(k))
   end do

   do k=1,nz
    htendv(k)=h_v(k)%bz*v(k)%tend1%bz
    call av_iteratation(htendv(k))
   end do



   do k=1,nz
    utendh(k)=Ax(u(k))*h(k)%tend1%bz
    call av_iteratation(utendh(k))
   end do

   do k=1,nz
    vtendh(k)=Ay(v(k))*h(k)%tend1%bz
    call av_iteratation(vtendh(k))
   end do


  end subroutine
  
  subroutine write_timemean(dat)
   use global
   use overload
   use writeout
   
   implicit none
   
   class(var), intent(inout) :: dat
   
   
   dat%t=dat%t%bz/int2real(timeav_count)
   call fullwrite_var('timemean',dat%t)
   
  end subroutine
  
  subroutine read_timemean(dat)
   use global
   use overload
   use writeout
   
   implicit none
   
   class(var), intent(inout) :: dat
   
   
   call read_var(dat%t,'timemean')
   dat%t=dat%t%bz*int2real(timeav_count)
   
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
    call write_timemean(h(k))
   end do
   
   do k=1,nz
    call write_timemean(hu(k))
   end do
   
   do k=1,nz
    call write_timemean(hv(k))
   end do

   do k=1,nz
    call write_timemean(u(k))
   end do

   do k=1,nz
    call write_timemean(v(k))
   end do

 
   do k=1,nz
    call write_timemean(huu(k))
   end do

   do k=1,nz
    call write_timemean(hvv(k))
   end do
   
   do k=1,nz
    call write_timemean(huv_h(k))
   end do
   
   do k=1,nz
    call write_timemean(huv_z(k))
   end do
   
   do k=1,nz
    call write_timemean(huv_uz(k))
   end do
   
   do k=1,nz
    call write_timemean(huv_vz(k))
   end do
   
   do k=1,nz
    call write_timemean(huuu(k))
   end do
   
   do k=1,nz
    call write_timemean(huuv(k))
   end do
   
   do k=1,nz
    call write_timemean(huvv(k))
   end do
   
   do k=1,nz
    call write_timemean(hvvv(k))
   end do
   
   do k=1,nz
    call write_timemean(hum(k))
   end do
   
   do k=1,nz
    call write_timemean(hvm(k))
   end do
   
   do k=1,nz
    call write_timemean(hm_x(k))
   end do
   
   do k=1,nz
    call write_timemean(hm_y(k))
   end do
   
   do k=1,nz
    call write_timemean(up_dm_x(k))
   end do
   
   do k=1,nz
    call write_timemean(up_dm_y(k))
   end do
   
   do k=1,nz
    call write_timemean(up_du(k))
   end do
   
   do k=1,nz
    call write_timemean(up_dv(k))
   end do
   
   do k=1,nz
    call write_timemean(dn_dm_x(k))
   end do
   
   do k=1,nz
    call write_timemean(dn_dm_y(k))
   end do
   
   do k=1,nz
    call write_timemean(dn_du(k))
   end do
   
   do k=1,nz
    call write_timemean(dn_dv(k))
   end do
   
   do k=0,nz
    call write_timemean(dd(k))
   end do
   
   do k=0,nz
    call write_timemean(d(k))
   end do
      
   do k=1,nz
    call write_timemean(m(k))
   end do

   do k=1,nz
    call write_timemean(hm(k))
   end do
          
 
   call write_timemean(utau)
   call write_timemean(vtau)
 
   call write_timemean(bfricu)
   call write_timemean(bfricv)
   
   do k=1,nz
    call write_timemean(hsmagu(k))
   end do
   
   do k=1,nz
    call write_timemean(hsmagv(k))
   end do
 
   call write_timemean(uutau)
   call write_timemean(vvtau)
 
   call write_timemean(ubfricu)
   call write_timemean(vbfricv)
   
   do k=1,nz
    call write_timemean(husmagu(k))
   end do
   
   do k=1,nz
    call write_timemean(hvsmagv(k))
   end do
   
   do k=1,nz
    call write_timemean(hq(k))
   end do
   
   do k=1,nz
    call write_timemean(huq(k))
   end do
   do k=1,nz
    call write_timemean(huq_v(k))
   end do
   
   do k=1,nz
    call write_timemean(hvq_u(k))
   end do
   do k=1,nz
    call write_timemean(hvq(k))
   end do
  
   do k=1,nz
    call write_timemean(hqq(k))
   end do
   
   do k=1,nz
    call write_timemean(q_h(k))
   end do
   do k=1,nz
    call write_timemean(q(k))
   end do
   
   do k=1,nz
    call write_timemean(uq(k))
   end do
   do k=1,nz
    call write_timemean(uq_v(k))
   end do
   
   do k=1,nz
    call write_timemean(vq_u(k))
   end do
   do k=1,nz
    call write_timemean(vq(k))
   end do
  
   do k=1,nz
    call write_timemean(qq(k))
   end do

   do k=1,nz
    call write_timemean(tendh(k))
   end do
   do k=1,nz
    call write_timemean(htendu(k))
   end do
   do k=1,nz
    call write_timemean(htendv(k))
   end do
   do k=1,nz
    call write_timemean(utendh(k))
   end do
   do k=1,nz
    call write_timemean(vtendh(k))
   end do





  end subroutine
   
  
  
  
  
  subroutine read_timemeans()
   use variables
   use parallel
   
   implicit none
   
   integer :: k
   character(32) :: filename
   logical :: dir_e
   
   
   
  write (filename, "(a11,i4.4,a2)") './timemean/', ens_name, '/.'
  inquire( file=filename, exist=dir_e )
  if (.not.(dir_e)) then
  
   timeav_count = 0
   
  else
  
   if (proc_name == ens_master) then
    write(format,'(a64)') '(a20,i10,a11,i4)'
    write(*,format) 'Starting Mean from: ', timeav_count, ' ens_name: ', ens_name
   end if 
   
   do k=1,nz
    call read_timemean(h(k))
   end do
   
   do k=1,nz
    call read_timemean(hu(k))
   end do
   
   do k=1,nz
    call read_timemean(hv(k))
   end do

   do k=1,nz
    call read_timemean(u(k))
   end do

   do k=1,nz
    call read_timemean(v(k))
   end do

   
   do k=1,nz
    call read_timemean(huu(k))
   end do
   
   do k=1,nz
    call read_timemean(hvv(k))
   end do
   
   do k=1,nz
    call read_timemean(huv_h(k))
   end do
   
   do k=1,nz
    call read_timemean(huv_z(k))
   end do
   
   do k=1,nz
    call read_timemean(huv_uz(k))
   end do
   
   do k=1,nz
    call read_timemean(huv_vz(k))
   end do
   
   do k=1,nz
    call read_timemean(huuu(k))
   end do

   do k=1,nz
    call read_timemean(huuv(k))
   end do

   do k=1,nz
    call read_timemean(huvv(k))
   end do

   do k=1,nz
    call read_timemean(hvvv(k))
   end do
   
   do k=1,nz
    call read_timemean(hum(k))
   end do

   do k=1,nz
    call read_timemean(hvm(k))
   end do
   
   do k=1,nz
    call read_timemean(hm_x(k))
   end do
   
   do k=1,nz
    call read_timemean(hm_y(k))
   end do
   
   do k=1,nz
    call read_timemean(up_dm_x(k))
   end do
   
   do k=1,nz
    call read_timemean(up_dm_y(k))
   end do
   
   do k=1,nz
    call read_timemean(up_du(k))
   end do
   
   do k=1,nz
    call read_timemean(up_dv(k))
   end do
   
   do k=1,nz
    call read_timemean(dn_dm_x(k))
   end do
   
   do k=1,nz
    call read_timemean(dn_dm_y(k))
   end do
   
   do k=1,nz
    call read_timemean(dn_du(k))
   end do
   
   do k=1,nz
    call read_timemean(dn_dv(k))
   end do
   
   do k=0,nz
    call read_timemean(dd(k))
   end do
   
   do k=0,nz
    call read_timemean(d(k))
   end do
   
   do k=1,nz
    call read_timemean(m(k))
   end do

   do k=1,nz
    call read_timemean(hm(k))
   end do
          
 
   call read_timemean(utau)
   call read_timemean(vtau)
 
   call read_timemean(bfricu)
   call read_timemean(bfricv)
   
   do k=1,nz
    call read_timemean(hsmagu(k))
   end do
   
   do k=1,nz
    call read_timemean(hsmagv(k))
   end do
 
   call read_timemean(uutau)
   call read_timemean(vvtau)
 
   call read_timemean(ubfricu)
   call read_timemean(vbfricv)
   
   do k=1,nz
    call read_timemean(husmagu(k))
   end do
   
   do k=1,nz
    call read_timemean(hvsmagv(k))
   end do
   
   do k=1,nz
    call read_timemean(hq(k))
   end do
   
   do k=1,nz
    call read_timemean(huq(k))
   end do
   do k=1,nz
    call read_timemean(huq_v(k))
   end do
   
   do k=1,nz
    call read_timemean(hvq_u(k))
   end do
   do k=1,nz
    call read_timemean(hvq(k))
   end do
  
   do k=1,nz
    call read_timemean(hqq(k))
   end do
   
   do k=1,nz
    call read_timemean(q_h(k))
   end do
   do k=1,nz
    call read_timemean(q(k))
   end do
   
   do k=1,nz
    call read_timemean(uq(k))
   end do
   do k=1,nz
    call read_timemean(uq_v(k))
   end do
   
   do k=1,nz
    call read_timemean(vq_u(k))
   end do
   do k=1,nz
    call read_timemean(vq(k))
   end do
  
   do k=1,nz
    call read_timemean(qq(k))
   end do

   do k=1,nz
    call read_timemean(tendh(k))
   end do
   do k=1,nz
    call read_timemean(htendu(k))
   end do
   do k=1,nz
    call read_timemean(htendv(k))
   end do
   do k=1,nz
    call read_timemean(utendh(k))
   end do
   do k=1,nz
    call read_timemean(vtendh(k))
   end do

   end if


  end subroutine
  

end module

#endif
