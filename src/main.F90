program sw

#include "include.h"
 use sys
 use parallel
 use grid
 use params
 use variables
 use grid_operate
 use sync
#ifdef ALLOW_RIGID_LID
 use solver
#endif
 use writeout
 use allocation
 use operations
 use llist_ops
#ifdef DO_TIME_AVERAGE
 use timeav
#endif

 
 
  implicit none
  
  integer :: k, n
  
  call mpi_init(stat)
  
  call mpi_comm_rank(mpi_comm_world, proc_name, stat)
  call mpi_comm_size(mpi_comm_world, proc_num, stat)
  
  ens_images = proc_num/ens_num
  ens_name = proc_name/ens_images
  ens_master = ens_name*ens_images
  
  call split_domain(mx, my, lx, nx, ly, ny, north, south, east, west)
  
  call set_grids
  call init_write_grids
  

 call create_field(s,'s',.true.)
 
 call create_field(d,'d',.true.)
 call create_field(mind,.false.)
 
 call create_field(h,'h',.true.)
 call create_field(u,'u',.true.)
 call create_field(v,'v',.true.)
 
 call create_field(minh,.false.)
 
 do k=1,nz
  call make_timestep(h(k))
  call make_timestep(u(k))
  call make_timestep(v(k))
 end do
 
 call create_field(h_u,.true.)
 call create_field(h_v,.true.)
 call create_field(h_z,.false.)
 
 call create_field(h_tmp,.true.)
 call create_field(u_tmp,.false.)
 call create_field(v_tmp,.false.)
 
 call create_field(ke,'ke',.true.)
 call create_field(ape,'ape',.false.)
 
 call create_field(zeta,'z',.false.)
 
 call create_field(lapu,'lapu',.true.)
 call create_field(lapv,'lapv',.true.)
 
 call create_field(smag,'smag',.true.)
 call create_field(tension,'tension',.true.)
 call create_field(strain,'strain',.false.)
 
 call create_field(smagu,'smagu',.false.)
 call create_field(smagv,'smagv',.false.)
 
 call create_field(q,'q',.false.)
  
 call create_field(m,'m',.true.)
 call create_field(minm,.false.)
 
 call create_field(bfricu,'bfricu',.false.)
 call create_field(bfricv,'bfricv',.false.)
 
 call create_field(utau,'utau',.false.)
 call create_field(vtau,'vtau',.false.)
 
 call create_field(f,'f',.false.)
 
 
 

#ifdef DO_TIME_AVERAGE

 call create_field(hu,'hu',.true.)
 call create_field(hv,'hv',.true.)

 call create_field(s_h_h,'s_h_h',.false.)
 call create_field(s_h_u,'s_h_u',.false.)
 call create_field(s_h_v,'s_h_v',.false.)
 call create_field(s_h_z,'s_h_z',.false.)

 call create_field(s_hu_h,'s_hu_h',.false.)
 call create_field(s_hu_u,'s_hu_u',.false.)
 call create_field(s_hu_v,'s_hu_v',.false.)
 call create_field(s_hu_z,'s_hu_z',.false.)

 call create_field(s_hv_h,'s_hv_h',.false.)
 call create_field(s_hv_u,'s_hv_u',.false.)
 call create_field(s_hv_v,'s_hv_v',.false.)
 call create_field(s_hv_z,'s_hv_z',.false.)

 call create_field(s_huu_h,'s_huu_h',.false.)

 call create_field(s_hvv_h,'s_hvv_h',.false.)

 call create_field(s_huv_h,'s_huv_h',.false.)
 call create_field(s_huv_z,'s_huv_z',.false.)
 call create_field(s_huv_uz,'s_huv_uz',.false.)
 call create_field(s_huv_vz,'s_huv_vz',.false.)
 
 call create_field(s_hm_x_u,'s_hm_x_u',.false.)
 call create_field(s_hm_y_v,'s_hm_y_v',.false.)
 
 call create_field(s_m_h,'s_m_h',.false.)
 
 call create_field(s_utau_u,'s_utau_u',.false.)
 call create_field(s_vtau_v,'s_vtau_v',.false.)
 
 call create_field(s_bfricu_u,'s_bfricu_u',.false.)
 call create_field(s_bfricv_v,'s_bfricv_v',.false.)
 
 call create_field(s_hsmagu_u,'s_hsmagu_u',.false.)
 call create_field(s_hsmagv_v,'s_hsmagv_v',.false.)
 

#endif
 
 
 s=-1.0d0+(0.25d0)*(max(1.0d0-(y_dist(s,y0/2.0d0)/(0.25d0*y0))**2,0.0d0) - &
            max(1.0d0-(y_dist(s,0.0d0)/(0.25d0*y0))**2,0.0d0))
 call start_sync(s)
 
 utau=0.0d0
 vtau=0.0d0
 
 f=1.0d0
 
 call end_sync(s)
 
#ifdef ALLOW_RIGID_LID
 call set_solver()
#endif
 
 mind(0)=s
 do k=1,nz
  mind(k)=(0.75d0*int2real(k-nz))/int2real(nz)
 end do
 do k=1,nz
  minh(k)=mind(k)-mind(k-1)
 end do
 
 do k=0,nz
  d(k)=mind(k)
 end do
 call mont(mind,minm)
 do k=1,nz
  h(k)=d(k)-d(k-1)
  u(k)=0.0d0
  v(k)=0.0d0
 end do
 
 
#ifdef ALLOW_RIGID_LID
  pres=0.0d0
  r=0.0d0
  p=0.0d0
#endif
 
 
 
 do k=1,nz
  call read_var(h(k),'old')
  call read_var(u(k),'old')
  call read_var(v(k),'old')
 end do
#ifdef ALLOW_RIGID_LID
 call read_var(pres,'old')
#endif
 
 do k=1,nz
  call start_sync(h(k))
  call start_sync(u(k))
  call start_sync(v(k))
 end do
#ifdef ALLOW_RIGID_LID
 call start_sync(pres)
 call start_sync(r)
 call start_sync(p)
#endif


#ifdef ALLOW_STOCHASTIC_WIND
  call input_llist(stochwind,'stochwind')
#endif
 
 call write_main('in')
 
 call depth(h,s,d)
 call mont(d,m)
 do k=1,nz
  call start_sync(m(k))
 end do
 
 do k=1,nz
  call end_sync(h(k))
  h_u(k)=Ax(h(k))
  call start_sync(h_u(k))
  h_v(k)=Ay(h(k))
  call start_sync(h_v(k))
 end do
 
 do k=1,nz
  ke(k)=0.5d0*(Ax(u(k)*u(k))+Ay(v(k)*v(k)))
  call start_sync(ke(k))
 end do
 
 
!--------------------------------------------------------------------------!
!--------------------------------------------------------------------------!
!--------------------------------------------------------------------------!
!--------------------------------------------------------------------------!
!--------------------------------------------------------------------------!
!--------------------------------------------------------------------------!
 

 
!--------------------------------------------------------------------------------------------------------------!

! Start RUN

!--------------------------------------------------------------------------------------------------------------!

  
 do n=0,nstop
 
 
  if (proc_name == proc_master) then
   if (update(2.0d0)) then
    write(format,"(a32)") "(a33,a8,a26,a1)"
    write(*,format, ADVANCE = "NO") print_time(), ' Model: ', &
      modeltime(n), CHAR(13)
   end if
  end if 
  
  do k=1,nz
   call end_sync(h_u(k))
   call end_sync(h_v(k))
   h_z(k)=0.5d0*(Ay(h_u(k))+Ax(h_v(k)))
  end do
  
  if (n == 0) then
   do k=1,nz
    h(k)%tmp=h(k)%bz
    u(k)%tmp=u(k)%bz
    v(k)%tmp=v(k)%bz
   end do
  end if
  
 
  
  call tend_h(n)
  
  do k=1,nz
   call end_sync(u(k))
   lapu(k)=sGxx(u(k))+sGyy(u(k))
   call start_sync(lapu(k))
   call end_sync(v(k))
   lapv(k)=sGxx(v(k))+sGyy(v(k))
   call start_sync(lapv(k))
  end do
  
 
#ifdef ALLOW_STOCHASTIC_WIND
  call set_wind(n)
#endif
  
  do k=1,nz
   call end_sync(u(k))
   call end_sync(v(k))
   strain(k)=Gx(v(k))+Gy(u(k))
   tension(k)=Gx(u(k))-Gy(v(k))
   smag(k)=max(dx,dy)**4*cvar*sqrt(tension(k)%bz**2+Ax(Ay(strain(k)%bz**2)))
   call start_sync(smag(k))
  end do
  
  
  
  do k=1,nz
   call end_sync(u(k))
   call end_sync(v(k))
   zeta(k)=Gx(v(k))-Gy(u(k))
  end do
  do k=1,nz
   call end_sync(lapu(k))
   call end_sync(lapv(k))
   call end_sync(smag(k))
   smagu(k)=-Ax(smag(k))*(sGxx(lapu(k))+sGyy(lapu(k)))
  end do
  call end_sync(ke(1))
  bfricu=bf*(Ax(ke(1))/(4.0d-4))*u(1)%bz
  call tend_u(n)
  do k=1,nz
   call end_sync(lapu(k))
   call end_sync(lapv(k))
   call end_sync(smag(k))
   smagv(k)=-Ay(smag(k))*(sGxx(lapv(k))+sGyy(lapv(k)))
  end do
  call end_sync(ke(1))
  bfricv=bf*(Ay(ke(1))/(4.0d-4))*v(1)%bz
  call tend_v(n)


#ifdef ALLOW_RIGID_LID

  if (n >= 3) then
   do k=1,nz
    h_tmp(k)=ab3(h(k))
   end do
  else
   if (n == 2) then
    do k=1,nz
     h_tmp(k)=ab2(h(k))
    end do
   else
    if (n == 1) then
     do k=1,nz
      h(k)=h(k)%tmp
      h_tmp(k)=rk2(h(k))
     end do
    else
     if (n ==0) then
      do k=1,nz
       h_tmp(k)=fe(h(k))
      end do
     end if
    end if
   end if
  end if
  call thickness_correct(h_tmp,s)
   
  do k=1,nz
   call start_sync(h_tmp(k))
  end do

  if (n >= 3) then
   do k=1,nz
    u_tmp(k)=ab3(u(k))
    v_tmp(k)=ab3(v(k))
   end do
  else
   if (n == 2) then
    do k=1,nz
     u_tmp(k)=ab2(u(k))
     v_tmp(k)=ab2(v(k))
    end do
   else
    if (n == 1) then
     do k=1,nz
      u(k)=u(k)%tmp
      v(k)=v(k)%tmp
      u_tmp(k)=rk2(u(k))
      v_tmp(k)=rk2(v(k))
     end do
    else
     if (n ==0) then
      do k=1,nz
       u_tmp(k)=fe(u(k))
       v_tmp(k)=fe(v(k))
      end do
     end if
    end if
   end if
  end if
   
  
  
   inty=0.0d0
   do k=1,nz
    call end_sync(h_tmp(k))
    y(k)=-Gx(u_tmp(k)%bz*Ax(h_tmp(k)))  &
           -Gy(v_tmp(k)%bz*Ay(h_tmp(k)))
 !      inty%z=inty%z+y(k)%z
    if (n >= 3) then
     inty=inty%bz+(12.0d0/23.0d0)*y(k)%bz/dt
    else
     if (n == 2) then
      inty=inty%bz+(2.0d0/3.0d0)*y(k)%bz/dt
     else
      if (n == 1) then
       inty=inty%bz+y(k)%bz/dt
      else
       if (n ==0) then
        inty=inty%bz+2.0d0*y(k)%bz/dt
       end if
      end if
     end if
    end if
   end do
  
  
  
   
   call surfpressure(stat,1000)
   call prescorrection(n)
   
#endif
 

#ifdef ALLOW_RIGID_LID
  do k=1,nz
   ape(k)=(h(k)%bz-minh(k)%bz)*(-pres%bz+m(k)%bz-minm(k)%bz)
  end do
#else
  do k=1,nz
   ape(k)=(h(k)%bz-minh(k)%bz)*(m(k)%bz-minm(k)%bz)
  end do
#endif

  call timeav_iteratation()
  
  if (mod(n,wstep) == 0) then
   if (n /= 0) then
    call do_write(n/wstep-1)
    if (proc_name == ens_master) then
     write(format,"(a32)") "(a6,a23,f10.1,a17,i2)"
     print *
     write(*,format) num_position(n/wstep-1), ' writeout done.  time: ', cputime(), ' secs. ens_name: ', ens_name
    end if
   end if
   call write()
#ifdef ALLOW_STOCHASTIC_WIND
   call writeout_llist(n/wstep,stochwind,'stochwind')
#endif
  end if



  if (n >= 3) then
   do k=1,nz
#ifdef ALLOW_RIGID_LID
    h(k)=h_tmp(k)%bz
#else
    h(k)=ab3(h(k))
#endif
   end do
  else
   if (n == 2) then
    do k=1,nz
#ifdef ALLOW_RIGID_LID
    h(k)=h_tmp(k)%bz
#else
     h(k)=ab2(h(k))
#endif
    end do
   else
    if (n == 1) then
     do k=1,nz
#ifdef ALLOW_RIGID_LID
    h(k)=h_tmp(k)%bz
#else
      h(k)=h(k)%tmp
      h(k)=rk2(h(k))
#endif
     end do
    else
     if (n ==0) then
      do k=1,nz
#ifdef ALLOW_RIGID_LID
    h(k)=h_tmp(k)%bz
#else
       h(k)=fe(h(k))
#endif
      end do
     end if
    end if
   end if
  end if
   
  do k=1,nz
   call start_sync(h(k))
  end do
  
  
    
  call depth(h,s,d)
  
  call mont(d,m)
  
  do k=1,nz
   call start_sync(m(k))
  end do
  
  
  
  
  
  
 
  do k=1,nz
   if (n >= 3) then
     u(k)=ab3(u(k))
     v(k)=ab3(v(k))
   else
    if (n == 2) then
     u(k)=ab2(u(k))
     v(k)=ab2(v(k))
    else
     if (n == 1) then
#ifndef ALLOW_RIGID_LID
      u(k)=u(k)%tmp
      v(k)=v(k)%tmp
#endif
      u(k)=rk2(u(k))
      v(k)=rk2(v(k))
     else
      if (n ==0) then
       u(k)=fe(u(k))
       v(k)=fe(v(k))
      end if
     end if
    end if
   end if
   call start_sync(u(k))
   call start_sync(v(k))
   ke(k)=0.5d0*(Ax(u(k)%bz**2)+Ay(v(k)%bz**2))
   call start_sync(ke(k))
   call end_sync(h(k))
   h_u(k)=Ax(h(k))
   h_v(k)=Ay(h(k))
   call start_sync(h_u(k))
   call start_sync(h_v(k))
#ifdef DO_TIME_AVERAGE
   hu(k)=h_u(k)%bz*u(k)%bz
   hv(k)=h_v(k)%bz*v(k)%bz
   call start_sync(hu(k))
   call start_sync(hv(k))
#endif
  end do
  
   
     
  
  if (n /= 0) then
   do k=1,nz
    call stepforward(h(k))
    call stepforward(u(k))
    call stepforward(v(k))
   end do
  end if
  
 end do
 
 
 
 call do_write(n/wstep)
 if (proc_name == ens_master) then
  write(format,"(a32)") "(a6,a23,f10.1,a17,i2)"
  print *
  write(*,format) num_position(n/wstep), ' writeout done.  time: ', cputime(), ' secs. ens_name: ', ens_name
 end if
 
 
 call write_main('out')
 
 call write_timemeans()
 
 
 
#ifdef ALLOW_STOCHASTIC_WIND
 call output_llist(stochwind,'stochwind')
#endif
 
 
 
 
 
 
 
  
 
  
 call mpi_finalize(stat)


end program
