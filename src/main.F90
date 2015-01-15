program sw
 use system
 use parallel
 use grid
 use params
 use variables
 
  implicit none
  
  integer :: i,j,k
  
  call mpi_init(stat)
  
  call mpi_comm_rank(mpi_comm_world, proc_name, stat)
  call mpi_comm_size(mpi_comm_world, proc_num, stat)
  
  ens_images = proc_num/ens_num
  ens_name = proc_name/ens_images
  ens_master = ens_name*ens_images
  
  call split_domain(mx, my, lx, nx, ly, ny, north, south, east, west)
  
  call set_grids()
  

 call create_field(s,.false.)
 
 call create_field(depth,.true.)
 call create_field(mindepth,.false.)
 
 call create_field(h,.true.)
 call create_field(u,.true.)
 call create_field(v,.true.)
 
 call create_field(h_tmp,.false.)
 call create_field(u_tmp,.false.)
 call create_field(v_tmp,.false.)
 
 call create_field(ke,.true.)
 call create_field(ape,.false.)
 
 call create_field(zeta,.false.)
 
 call create_field(lapu,.true.)
 call create_field(lapv,.true.)
 
 call create_field(smag,.true.)
 call create_field(tension,.true.)
 call create_field(strain,.false.)
 
 call create_field(smagu,.false.)
 call create_field(smagv,.false.)
 
 call create_field(q,.false.)
  
 call create_field(m,.true.)
 call create_field(minm,.false.)
 
 call create_field(bfricu,.false.)
 call create_field(bfricv,.false.)
 
 call create_field(utau,.false.)
 call create_field(vtau,.false.)
 
 call create_field(f,.false.)
  
  
  
  call mpi_finalize(stat)


end program
