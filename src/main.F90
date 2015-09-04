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
  
  call init_parallel()
  
  
  call init_params()
  call mpi_barrier(mpi_comm_world,stat)
  if (proc_name == proc_master) call write_params()
   
  call split_domain()
  
  call set_grids
  call init_write_grids
  

 call create_field(s,'s',.true.)
 
 call create_field(d,'d',.true.,0)
 call create_field(mind,.false.,0)
 call create_field(d_u,'d',.true.,0)
 call create_field(d_v,'d',.true.,0)
 
 call create_field(dd,'dd',.false.,0)
 
 call create_field(h,'h',.true.,1)
 call create_field(u,'u',.true.,1)
 call create_field(v,'v',.true.,1)
 
 call create_field(minh,.false.,1)
 
 do k=1,nz
  call make_timestep(h(k))
  call make_timestep(u(k))
  call make_timestep(v(k))
 end do
 
 call create_field(h_tmp,.true.,1)
 call create_field(u_tmp,.false.,1)
 call create_field(v_tmp,.false.,1)
 
 call create_field(h_u,'h_u',.true.,1)
 call create_field(h_v,'h_v',.true.,1)
 call create_field(h_z,'h_z',.false.,1)

 call create_field(hu,'hu',.true.,1)
 call create_field(hv,'hv',.true.,1)
 
 call create_field(hu_h,'hu_h',.false.,1)
 call create_field(hu_v,'hu_v',.false.,1)
 call create_field(hu_z,'hu_z',.false.,1)
 
 call create_field(hv_h,'hv_h',.false.,1)
 call create_field(hv_u,'hv_u',.false.,1)
 call create_field(hv_z,'hv_z',.false.,1)

 call create_field(u_h,'u_h',.false.,1)
 call create_field(u_v,'u_v',.false.,1)
 call create_field(u_z,'u_z',.false.,1)

 call create_field(v_h,'v_h',.false.,1)
 call create_field(v_u,'v_u',.false.,1)
 call create_field(v_z,'v_z',.false.,1)
 
 call create_field(huu_h,'huu_h',.false.,1)
 call create_field(huu_z,'huu_z',.false.,1)

 call create_field(hvv_h,'hvv_h',.false.,1)
 call create_field(hvv_z,'hvv_z',.false.,1)
 
 call create_field(huv_h,'huv_h',.false.,1)
 call create_field(huv_z,'huv_z',.false.,1)
 call create_field(huv_uz,'huv_uz',.false.,1)
 call create_field(huv_vz,'huv_vz',.false.,1)
 
 call create_field(huuu,'huuu',.false.,1)
 call create_field(huuv,'huuv',.false.,1)
 call create_field(huvv,'huvv',.false.,1)
 call create_field(hvvv,'hvvv',.false.,1)
 
 call create_field(hum,'hum',.false.,1)
 call create_field(hvm,'hvm',.false.,1)
 
 call create_field(hm_x_u,'hm_x_u',.false.,1)
 call create_field(hm_y_v,'hm_y_v',.false.,1)
 
 call create_field(up_dm_x_u,'up_dm_x_u',.false.,1)
 call create_field(up_dm_y_v,'up_dm_y_v',.false.,1)
 
 call create_field(up_du_u,'up_du_u',.false.,1)
 call create_field(up_dv_v,'up_dv_v',.false.,1)
 
 call create_field(dn_dm_x_u,'dn_dm_x_u',.false.,1)
 call create_field(dn_dm_y_v,'dn_dm_y_v',.false.,1)
 
 call create_field(dn_du_u,'dn_du_u',.false.,1)
 call create_field(dn_dv_v,'dn_dv_v',.false.,1)
 
 call create_field(ke,'ke',.true.,1)
 call create_field(ape,'ape',.false.,1)
 
 call create_field(zeta,'z',.false.,1)
 
 call create_field(lapu,'lapu',.true.,1)
 call create_field(lapv,'lapv',.true.,1)
 
 call create_field(smag,'smag',.true.,1)
 call create_field(tension,'tension',.true.,1)
 call create_field(strain,'strain',.false.,1)
 
 call create_field(smagu,'smagu',.false.,1)
 call create_field(smagv,'smagv',.false.,1)
 
 call create_field(hsmagu,'hsmagu',.false.,1)
 call create_field(hsmagv,'hsmagv',.false.,1)
 
 call create_field(husmagu,'husmagu',.false.,1)
 call create_field(hvsmagv,'hvsmagv',.false.,1)
 
 call create_field(q_h,'q_h',.false.,1)
 call create_field(q_u,'q_u',.false.,1)
 call create_field(q_v,'q_v',.false.,1)
 call create_field(q,'q',.false.,1)

 call create_field(uq_h,'uq_h',.false.,1)
 call create_field(uq,'uq',.false.,1)
 call create_field(uq_v,'uq_v',.false.,1)
 call create_field(uq_z,'uq_z',.false.,1)

 call create_field(vq_h,'vq_h',.false.,1)
 call create_field(vq_u,'vq_u',.false.,1)
 call create_field(vq,'vq',.false.,1)
 call create_field(vq_z,'vq_z',.false.,1)
 
 call create_field(qq,'qq',.false.,1)

 call create_field(hq_h,'hq_h',.false.,1)
 call create_field(hq_u,'hq_u',.false.,1)
 call create_field(hq_v,'hq_v',.false.,1)
 call create_field(hq,'hq',.false.,1)

 call create_field(huq_h,'huq_h',.false.,1)
 call create_field(huq,'huq',.false.,1)
 call create_field(huq_v,'huq_v',.false.,1)
 call create_field(huq_z,'huq_z',.false.,1)

 call create_field(hvq_h,'hvq_h',.false.,1)
 call create_field(hvq_u,'hvq_u',.false.,1)
 call create_field(hvq,'hvq',.false.,1)
 call create_field(hvq_z,'hvq_z',.false.,1)
 
 call create_field(hqq,'hqq',.false.,1)
  
 call create_field(m,'m',.true.,1)
 call create_field(minm,.false.,1)
 
 call create_field(hm,'hm',.false.,1)
 
 call create_field(bfricu,'bfricu',.false.)
 call create_field(bfricv,'bfricv',.false.)
 
 call create_field(utau,'utau',.false.)
 call create_field(vtau,'vtau',.false.)
 
 call create_field(ubfricu,'ubfricu',.false.)
 call create_field(vbfricv,'vbfricv',.false.)
 
 call create_field(uutau,'uutau',.false.)
 call create_field(vvtau,'vvtau',.false.)
 
 call create_field(f,'f',.false.)
 
 
 
 call create_field(tendh,'tendh',.false.,1) 
 call create_field(htendu,'htendu',.false.,1)
 call create_field(htendv,'htendv',.false.,1)
 call create_field(utendh,'utendh',.false.,1)
 call create_field(vtendh,'vtendh',.false.,1)


 s=-1.0d0!+(0.25d0)*(max(1.0d0-(y_dist(s,y0/2.0d0)/(0.25d0*y0))**8,0.0d0) - &
         !   max(1.0d0-(y_dist(s,0.0d0)/(0.25d0*y0))**8,0.0d0))
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
! d(1)=d(1)%bz+(0.125d0)*(max(1.0d0-(y_dist(s,y0/2.0d0)/(0.25d0*y0))**2,0.0d0) - &
!            max(1.0d0-(y_dist(s,0.0d0)/(0.25d0*y0))**2,0.0d0))
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



 
 call depth(h,s%bz,d)
 call mont(d,m)
 do k=1,nz
  call start_sync(m(k))
 end do
 call geostrophic_balance() 

 
 
 
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
 
 call depth(h,s%bz,d)
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
 
 
 call read_timemeans()
 
 
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
   call depth(h_u,Ax(s),d_u)
   call depth(h_v,Ay(s),d_v)
   h_z(k)=0.5d0*(Ay(h_u(k))+Ax(h_v(k)))
   zeta(k)%bz=Gx(v(k))-Gy(u(k))
   q(k)=(f%bz+zeta(k)%bz)*merge(0.0d0,1.0d0/h_z(k)%bz,(h_z(k)%bz == 0.0d0))
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
  bfricu=bf*sqrt(Ax(ke(1))/umax)**bfricpoly*u(1)%bz
  call tend_u(n)
  do k=1,nz
   call end_sync(lapu(k))
   call end_sync(lapv(k))
   call end_sync(smag(k))
   smagv(k)=-Ay(smag(k))*(sGxx(lapv(k))+sGyy(lapv(k)))
  end do
  call end_sync(ke(1))
  bfricv=bf*sqrt(Ay(ke(1))/umax)**bfricpoly*v(1)%bz
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
   thavx=0.0d0
   thavy=0.0d0
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
    thavx=thavx%bz+invrho(k)*Ax(h_tmp(k))
    thavy=thavy%bz+invrho(k)*Ay(h_tmp(k))
   end do
   inthavx=merge(0.0d0,1.0d0/thavx%bz,(thavx%bz == 0.0d0))
   inthavy=merge(0.0d0,1.0d0/thavy%bz,(thavy%bz == 0.0d0))

  
  
   
   call surfpressure(stat,1000)
   call prescorrection(n)
   
#endif
 

  do k=1,nz
   ape(k)=(h(k)%bz-minh(k)%bz)*(m(k)%bz-minm(k)%bz)
  end do


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
  
  
    
  call depth(h,s%bz,d)
  
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
 
 call mpi_barrier(mpi_comm_world,stat)
 call random_seed(GET=put_seed)
 seed=put_seed(1)
 if (proc_name == proc_master) call write_params()
 
 
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
