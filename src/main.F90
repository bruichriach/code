program sw
 use sys
 use parallel
 use grid
 use params
 use variables
 use grid_operate
 use sync
 use solver
 use solver_variables
 use writeout
 use allocation
 use operations
 
 
  implicit none
  
  integer :: k
  
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
 
 call create_field(h_u,.true.)
 call create_field(h_v,.true.)
 call create_field(h_z,.false.)
 
 allocate(tendh1(nz),tendh2(nz),tendh3(nz))
 allocate(tendu1(nz),tendu2(nz),tendu3(nz))
 allocate(tendv1(nz),tendv2(nz),tendv3(nz))
 
 call create_field(tendh1,'tendh',.false.)
 call create_field(tendu1,'tendu',.false.)
 call create_field(tendv1,'tendv',.false.)
 
 call create_field(tendh2,.false.)
 call create_field(tendu2,.false.)
 call create_field(tendv2,.false.)
 
 call create_field(tendh3,.false.)
 call create_field(tendu3,.false.)
 call create_field(tendv3,.false.)
 
 call create_field(h_tmp,.false.)
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
 
 
 s=-1.0d0+(0.25d0)*(max(1.0d0-(y_dist(s,y0/2.0d0)/(0.25d0*y0))**2,0.0d0) - &
            max(1.0d0-(y_dist(s,0.0d0)/(0.25d0*y0))**2,0.0d0))
 call start_sync(s)
 
 f=1.0d0
 
 call end_sync(s)
 
 call set_solver()
 
 mind(0)=s
 do k=1,nz
  mind(k)=(0.75d0*int2real(k-nz))/int2real(nz)
 end do
 
 d=mind
 call mont(mind,minm)
 do k=1,nz
  h(k)=d(k)-d(k-1)
  u(k)=0.0d0
  v(k)=0.0d0
 end do
 
 do k=1,nz
  call read_var(h(k),'old')
  call read_var(u(k),'old')
  call read_var(v(k),'old')
 end do
 call read_var(pres,'old')
 
 do k=1,nz
  call start_sync(h(k))
  call start_sync(u(k))
  call start_sync(v(k))
 end do
 call start_sync(pres)
 
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
 
 
 do k=1,nz
  call end_sync(h_u(k))
  call end_sync(h_v(k))
  h_z(k)=0.5d0*(Ay(h_u(k))+Ax(h_v(k)))
 end do
 
!--------------------------------------------------------------------------!
!--------------------------------------------------------------------------!
!--------------------------------------------------------------------------!
!--------------------------------------------------------------------------!
!--------------------------------------------------------------------------!
!--------------------------------------------------------------------------!
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 call end_sync(h(1))
  
 v(2)=2.0d0*pi*v(1)%p%y/y0
 call start_sync(v(2))
 call end_sync(v(2))
 
 v(1)=cos(v(2)%bz)
 call start_sync(v(1))
 call end_sync(v(1))
 
 
 
 inty=Gy(v(1))
 
 call surfpressure(stat)
 
 call end_sync(pres)
 
 
 if (proc_name == 0) then
  call print_var(pres)
 end if
 
 
  
 call mpi_finalize(stat)


end program
