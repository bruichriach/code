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
 
 
  implicit none
  
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
 
 call create_field(depth,'d',.true.)
 call create_field(mindepth,.false.)
 
 call create_field(h,'h',.true.)
 call create_field(u,'u',.true.)
 call create_field(v,'v',.true.)
 
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
 
 
 
 s=-1.0d0
 call start_sync(s)
 call end_sync(s)
 
 h(1)=1.0d0
 call start_sync(h(1))
 call end_sync(h(1))
 
 call set_solver()
  
 v(2)=2.0d0*pi*v(1)%p%y/y0
 call start_sync(v(2))
 call end_sync(v(2))
 
 v(1)=cos(v(2)%bz)
 call start_sync(v(1))
 call end_sync(v(1))
 
 
 call write_main('in')
 
 inty=Gy(v(1))
 
 call surfpressure(stat)
 
 call end_sync(pres)
 
 
 if (proc_name == 0) then
  call print_var(pres)
 end if
 
 
  
 call mpi_finalize(stat)


end program
