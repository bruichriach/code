module grid

#include "include.h"

 implicit none
 
 
  
  interface create_field
   module procedure create_hvar, create_uvar,  &
         create_vvar, create_zvar
   module procedure create_hvar_layer, create_uvar_layer,  &
         create_vvar_layer, create_zvar_layer
  end interface
 
 contains
 
  subroutine split_domain(mx, my, lx, nx, ly, ny, north, south, east, west)
   use parallel
   
   implicit none
      
   integer :: core_x, core_y, node_x, node_y
   integer, intent(in) :: mx, my
   integer :: kx, ky, centre(2)
   integer, intent(out) :: lx, ly, nx, ny
   integer, intent(out) :: north, south, east, west
   integer :: core, node
   integer :: grid_core, grid_node
   integer, allocatable :: node_names(:,:), core_names(:,:), names(:,:)
   integer :: node_grid(mx,my)
   integer :: i, j, k, i1, j1
   
   if (max_core > proc_num) then
    core = proc_num
   else
    core = max_core
   end if
   
   node=proc_num/core
   
   if (ens_num <= node) then
    grid_node = node/ens_num
    grid_core = core
   else
    grid_node = 1
    grid_core = ens_images
   end if
   
   call split_axes(mx, my, grid_node, node_x, node_y)
   
   allocate(node_names(node_x,node_y))
   
   do j=1,node_y
    do i=1,node_x
     node_names(i,j)=(j-1)*node_x+(i-1)
    end do
   end do
   
!   if (proc_name == proc_master) print *, node_names
   
   do j=1,my
    do i=1,mx
     node_grid(i,j)=node_names(   &
     ceiling(dble(node_x)*(dble(i)-0.5d0)/dble(mx)),  &
     ceiling(dble(node_y)*(dble(j)-0.5d0)/dble(my)))
    end do
   end do
   
!   if (proc_name == proc_master) print *, node_grid
   
   kx=0
   ky=0
   do k=0,grid_node-1
    kx=max(kx,count((count((node_grid == k),2) /= 0),1))
    ky=max(ky,count((count((node_grid == k),1) /= 0),1))
   end do
   call split_axes(kx,ky, grid_core, core_x, core_y)
   
   allocate(core_names(core_x,core_y))
   
   do j=1,core_y
    do i=1,core_x
     core_names(i,j)=(j-1)*core_x+(i-1)
    end do
   end do
   
!   if (proc_name == proc_master) print *, core_names
   
   allocate(names(0:node_x*core_x+1,0:node_y*core_y+1))
   
   do j1=1,node_y
    do i1=1,node_x
     do j=1,core_y
      do i=1,core_x
       names(core_x*(i1-1)+i,core_y*(j1-1)+j) = &
         grid_core*node_names(i1,j1) + core_names(i,j)
      end do
     end do
    end do
   end do
   
   names(0,1:node_y*core_y)=names(node_x*core_x,1:node_y*core_y)
   names(node_x*core_x+1,1:node_y*core_y)=names(1,1:node_y*core_y)
   names(1:node_x*core_x,0)=names(1:node_x*core_x,node_y*core_y)
   names(1:node_x*core_x,node_y*core_y+1)=names(1:node_x*core_x,1)
         
!   if (proc_name == proc_master) print *, names
   
!   if (proc_name == proc_master) print *, shape(names)

   centre = minloc(abs(names(1:node_x*core_x,1:node_y*core_y)-proc_name))
   west=names(centre(1)-1,centre(2))
   east=names(centre(1)+1,centre(2))
   north=names(centre(1),centre(2)+1)
   south=names(centre(1),centre(2)-1)
   
!   print *, proc_name, west, east, north, south
   
   allocate(proc_grid(0:mx+1,0:my+1))
   proc_grid = -1
   
   do j=1,my
    do i=1,mx
     proc_grid(i,j)=names(   &
     ceiling(dble(node_x*core_x)*(dble(i)-0.5d0)/dble(mx)),  &
     ceiling(dble(node_y*core_y)*(dble(j)-0.5d0)/dble(my)))
    end do
   end do
   
   
   nx=count((count((proc_grid(1:mx,1:my) == proc_name),2) /= 0))
   ny=count((count((proc_grid(1:mx,1:my) == proc_name),1) /= 0))
   lx = maxloc(count((proc_grid(1:mx,1:my) == proc_name),2),1)-1
   ly = maxloc(count((proc_grid(1:mx,1:my) == proc_name),1),1)-1
   
!   if (proc_name == proc_master) print *, grid
   
   
   
!   print *, proc_name, lx, ly, nx, ny
   

   
  end subroutine
 
  subroutine split_axes(mx, my, n, nx, ny)
   use sys
   use parallel
   
   implicit none
   
   integer, intent(in) :: mx, my, n
   integer, intent(out) :: nx, ny
   logical :: factor(n)
   integer, allocatable :: factors(:)
   integer, allocatable :: factor_pairs(:,:)
   integer :: i, j
   
   do i=1,n
    factor(i)=is_factor(i,n)
   end do
   
   allocate(factors(count(factor)))
   
   j=0
   do i=1,n
    if (factor(i)) then
     j=j+1
     factors(j)=i
    end if
   end do
   
   allocate(factor_pairs((count(factor)+1)/2,2))
   
   do j=1,ubound(factor_pairs,1)
    factor_pairs(j,1)=factors(j)
    factor_pairs(j,2)=factors(ubound(factors,1)+1-j)
   end do
   
   j=minloc(ceiling(dble(min(mx,my))/dble(factor_pairs(:,1)))+  &
            ceiling(dble(max(mx,my))/dble(factor_pairs(:,2))),1)
   
   if (my == max(mx,my)) then    
    nx=minval(factor_pairs(j,:))
    ny=maxval(factor_pairs(j,:))
   else
    nx=maxval(factor_pairs(j,:))
    ny=minval(factor_pairs(j,:))
   end if
   
!   if (proc_name == proc_master) print *, nx, ny, n, mx, my
   
   
  end subroutine
  
  subroutine set_grids()
   use global
   use params
   
   integer :: i,j
   
   allocate(hgrid%x(lx+1:lx+nx,ly+1:ly+ny),   &
           hgrid%y(lx+1:lx+nx,ly+1:ly+ny))
   allocate(hgrid%i(lx+1:lx+nx,ly+1:ly+ny),   &
           hgrid%j(lx+1:lx+nx,ly+1:ly+ny))
   do j=lbound(hgrid%x,2),ubound(hgrid%x,2)
    do i=lbound(hgrid%x,1),ubound(hgrid%x,1)
     hgrid%x(i,j)=dx*(int2real(i)-0.5d0)
     hgrid%y(i,j)=dy*(int2real(j)-0.5d0)
     hgrid%i(i,j)=i
     hgrid%j(i,j)=j
    end do
   end do
   hgrid%lx=lx
   hgrid%ly=ly
   hgrid%ox=0
   hgrid%oy=0
   hgrid%nx=nx
   hgrid%ny=ny
   
   
   
   allocate(ugrid%x(lx+1:lx+nx+1,ly+1:ly+ny),   &
           ugrid%y(lx+1:lx+nx+1,ly+1:ly+ny))
   allocate(ugrid%i(lx+1:lx+nx+1,ly+1:ly+ny),   &
           ugrid%j(lx+1:lx+nx+1,ly+1:ly+ny))
   do j=lbound(ugrid%x,2),ubound(ugrid%x,2)
    do i=lbound(ugrid%x,1),ubound(ugrid%x,1)
     ugrid%x(i,j)=dx*int2real(i-1)
     ugrid%y(i,j)=dy*(int2real(j)-0.5d0)
     ugrid%i(i,j)=i
     ugrid%j(i,j)=j
    end do
   end do
   ugrid%lx=lx
   ugrid%ly=ly
   ugrid%ox=1
   ugrid%oy=0
   ugrid%nx=nx+1
   ugrid%ny=ny
   
   allocate(vgrid%x(lx+1:lx+nx,ly+1:ly+ny+1),   &
           vgrid%y(lx+1:lx+nx,ly+1:ly+ny+1))
   allocate(vgrid%i(lx+1:lx+nx,ly+1:ly+ny+1),   &
           vgrid%j(lx+1:lx+nx,ly+1:ly+ny+1))
   do j=lbound(vgrid%x,2),ubound(vgrid%x,2)
    do i=lbound(vgrid%x,1),ubound(vgrid%x,1)
     vgrid%x(i,j)=dx*(int2real(i)-0.5d0)
     vgrid%y(i,j)=dy*int2real(j-1)
     vgrid%i(i,j)=i
     vgrid%j(i,j)=j
    end do
   end do
   vgrid%lx=lx
   vgrid%ly=ly
   vgrid%ox=0
   vgrid%oy=1
   vgrid%nx=nx
   vgrid%ny=ny+1
   
   allocate(zgrid%x(lx+1:lx+nx+1,ly+1:ly+ny+1),   &
           zgrid%y(lx+1:lx+nx+1,ly+1:ly+ny+1))
   allocate(zgrid%i(lx+1:lx+nx+1,ly+1:ly+ny+1),   &
           zgrid%j(lx+1:lx+nx+1,ly+1:ly+ny+1))
   do j=lbound(zgrid%x,2),ubound(zgrid%x,2)
    do i=lbound(zgrid%x,1),ubound(zgrid%x,1)
     zgrid%x(i,j)=dx*int2real(i-1)
     zgrid%y(i,j)=dy*int2real(j-1)
     zgrid%i(i,j)=i
     zgrid%j(i,j)=j
    end do
   end do
   zgrid%lx=lx
   zgrid%ly=ly
   zgrid%ox=1
   zgrid%oy=1
   zgrid%nx=nx+1
   zgrid%ny=ny+1
  
  end subroutine
  
  
  
  subroutine create_var(dat,grid,synced)
   use global
   use params
   use parallel
   
   implicit none
   
   class(var), intent(out) :: dat
   type(grd), target, intent(in) :: grid
   integer :: i,j,k,l
   logical, intent(in) :: synced
   
   dat%p => grid
      
   allocate(dat%z(dat%p%lx:dat%p%lx+dat%p%nx+1,   &
         dat%p%ly:dat%p%ly+dat%p%ny+1))
   dat%bz =>  &
         dat%z(dat%p%lx+1:dat%p%lx+dat%p%nx,   &
         dat%p%ly+1:dat%p%ly+dat%p%ny)
   dat%bz(dat%p%lx+1:,dat%p%ly+1:) => dat%bz
   
   
   if (synced) then
       
    dat%tag=max_tag
    max_tag=max_tag+1
    allocate(dat%mpi(0:proc_num-1))
   
    do k=0,proc_num-1
     j=0
     if (north == k) j=j+dat%p%nx
     if (south == k) j=j+dat%p%nx
     if (east == k) j=j+dat%p%ny
     if (west == k) j=j+dat%p%ny
     if (j /= 0) then
      allocate(dat%mpi(k)%s(j),dat%mpi(k)%r(j),    &
          dat%mpi(k)%s_mpi(j),dat%mpi(k)%r_mpi(j))
      dat%mpi(k)%r%slip=unity
      dat%mpi(k)%r_req=mpi_request_null
      dat%mpi(k)%s_req=mpi_request_null
     end if
    end do
    
    
    do j=0,proc_num-1
     l=0
     if (north == j) then
      do i=1,dat%p%nx
       dat%mpi(j)%s(l+i) =  &
         point_at(dat%z(dat%p%lx+i,dat%p%ly+dat%p%ny-dat%p%oy))
      end do
      l=l+dat%p%nx
     end if
     if (south == j) then
      do i=1,dat%p%nx
       dat%mpi(j)%s(l+i) =  &
         point_at(dat%z(dat%p%lx+i,dat%p%ly+1+dat%p%oy))
      end do
      l=l+dat%p%nx
     end if
     if (east == j) then
      do i=1,dat%p%ny
       dat%mpi(j)%s(l+i) =  &
         point_at(dat%z(dat%p%lx+dat%p%nx-dat%p%ox,dat%p%ly+i))
      end do
      l=l+dat%p%ny
     end if
     if (west == j) then
      do i=1,dat%p%ny
       dat%mpi(j)%s(l+i) =  &
         point_at(dat%z(dat%p%lx+1+dat%p%ox,dat%p%ly+i))
      end do
      l=l+dat%p%ny
     end if
    end do
   
   
    do j=0,proc_num-1
     l=0
     if (south == j) then
      do i=1,dat%p%nx
       dat%mpi(j)%r(l+i)%z =  &
         point_at(dat%z(dat%p%lx+i,dat%p%ly))
      end do
      l=l+dat%p%nx
     end if
     if (north == j) then
      do i=1,dat%p%nx
       dat%mpi(j)%r(l+i)%z =  &
         point_at(dat%z(dat%p%lx+i,dat%p%ly+dat%p%ny+1))
      end do
      l=l+dat%p%nx
     end if
     if (west == j) then
       do i=1,dat%p%ny
         dat%mpi(j)%r(l+i)%z =  &
          point_at(dat%z(dat%p%lx,dat%p%ly+i))
       end do
      l=l+dat%p%ny
     end if
     if (east == j) then
      do i=1,dat%p%ny
        dat%mpi(j)%r(l+i)%z =  &
         point_at(dat%z(dat%p%lx+dat%p%nx+1,dat%p%ly+i))
      end do
      l=l+dat%p%ny
     end if
    end do
    
   end if
    
  end subroutine
  
  
  subroutine create_hvar(dat,synced)
   use global
   
   implicit none
   
   type(hvar), intent(out) :: dat
   logical, intent(in) :: synced
   
   call create_var(dat,hgrid,synced)
   
  end subroutine
  
  
  subroutine create_uvar(dat,synced)
   use global
   
   implicit none
   
   type(uvar), intent(out) :: dat
   logical, intent(in) :: synced
   
   call create_var(dat,ugrid,synced)
   
  end subroutine
  
  
  subroutine create_vvar(dat,synced)
   use global
   
   implicit none
   
   type(vvar), intent(out) :: dat
   logical, intent(in) :: synced
   
   call create_var(dat,vgrid,synced)
   
  end subroutine
  
  
  subroutine create_zvar(dat,synced)
   use global
   
   implicit none
   
   type(zvar), intent(inout) :: dat
   logical, intent(in) :: synced
   
   call create_var(dat,zgrid,synced)
   
  end subroutine
  
  
  subroutine create_hvar_layer(dat,synced)
   use global
   
   implicit none
   
   type(hvar), intent(inout) :: dat(:)
   logical, intent(in) :: synced
   integer :: k
   
   do k=lbound(dat,1),ubound(dat,1)
    call create_hvar(dat(k),synced)
   end do
   
  end subroutine
  
  
  subroutine create_uvar_layer(dat,synced)
   use global
   
   implicit none
   
   type(uvar), intent(inout) :: dat(:)
   logical, intent(in) :: synced
   integer :: k
   
   do k=lbound(dat,1),ubound(dat,1)
    call create_uvar(dat(k),synced)
   end do
   
  end subroutine
  
  
  subroutine create_vvar_layer(dat,synced)
   use global
   
   implicit none
   
   type(vvar), intent(inout) :: dat(:)
   logical, intent(in) :: synced
   integer :: k
   
   do k=lbound(dat,1),ubound(dat,1)
    call create_vvar(dat(k),synced)
   end do
   
  end subroutine
  
  
  subroutine create_zvar_layer(dat,synced)
   use global
   
   implicit none
   
   type(zvar), intent(inout) :: dat(:)
   logical, intent(in) :: synced
   integer :: k
   
   do k=lbound(dat,1),ubound(dat,1)
    call create_zvar(dat(k),synced)
   end do
   
  end subroutine
  
  
 
end module





module sync

#include "include.h"
 
 implicit none
 
 contains

 
 subroutine start_sync(dat)
  use global
  use parallel
  
  implicit none
 
  class(var), intent(inout) :: dat
  integer :: i,k
  
  call end_sync(dat)
  if (dat%synced) then
   do k=0,proc_num-1
    if (allocated(dat%mpi(k)%r_mpi)) then
     call mpi_wait(dat%mpi(k)%r_req,mpi_status_ignore,stat)
     call mpi_irecv(dat%mpi(k)%r_mpi(1),ubound(dat%mpi(k)%r_mpi,1),mpi_precision, &
        k,dat%tag,mpi_comm_world,dat%mpi(k)%r_req,stat)
    end if
    if (allocated(dat%mpi(k)%s_mpi)) then
     call mpi_wait(dat%mpi(k)%s_req,mpi_status_ignore,stat)
     do i=1,ubound(dat%mpi(k)%s_mpi,1)
      dat%mpi(k)%s_mpi(i)=dat%mpi(k)%s(i)%z
     end do
     call mpi_isend(dat%mpi(k)%s_mpi(1),ubound(dat%mpi(k)%s_mpi,1),mpi_precision, &
        k,dat%tag,mpi_comm_world,dat%mpi(k)%s_req,stat)
    end if
   end do
   dat%synced = .false.
  end if
  
   
 end subroutine
 
 subroutine end_sync(dat)
  use global
  use parallel
  
  implicit none
 
  class(var), intent(inout) :: dat
  integer :: i,k
 
  
  if (.not.(dat%synced)) then
   do k=0,proc_num-1
    if (allocated(dat%mpi(k)%r_mpi)) then
     call mpi_wait(dat%mpi(k)%r_req,mpi_status_ignore,stat)
     do i=1,ubound(dat%mpi(k)%r_mpi,1)
      dat%mpi(k)%r(i)%z%z=sign(dat%mpi(k)%r_mpi(i),dat%mpi(k)%r_mpi(i)*dat%mpi(k)%r(i)%slip)
     end do
    end if
   end do
   dat%synced = .true.
  end if
  
   
  
 end subroutine
 
 
  
  subroutine integer_allmax(in, out)
   use parallel
   use sys
   
   implicit none
  
   integer, intent(in) :: in
   integer, intent(out) :: out
   integer :: out_tmp(0:ens_images-1), send(0:ens_images-1), recv(0:ens_images-1), i
   
   
   out = 0
   
   do i=0,ens_images-1
    call mpi_irecv(out_tmp(i),1,     &
             MPI_INTEGER,ens_master+i,40000,    &
             MPI_COMM_WORLD,recv(i),stat)
    call mpi_isend(in,1,     &
             MPI_INTEGER,ens_master+i,40000,    &
             MPI_COMM_WORLD,send(i),stat)
   end do
   
  
   do i=0,ens_images-1
    call mpi_wait(recv(i),MPI_STATUS_IGNORE,stat)
    out=max(out,out_tmp(i))
    call mpi_waitall(ens_images,send,MPI_STATUSES_IGNORE,stat)
   end do
   
  
  end subroutine


  subroutine real_allsum(in, out)
   use parallel
   use sys

   implicit none
  
   real (kind=db), intent(in) :: in
   real (kind=db), intent(out) :: out
   real (kind=db) :: out_tmp(0:ens_images-1)
   integer :: send(0:ens_images-1), recv(0:ens_images-1), i
   
   
   out = 0.0d0
  
   mobile_tag=mod(mobile_tag+1,1000)
 
   do i=0,ens_images-1
    call mpi_irecv(out_tmp(i),1,     &
             mpi_precision,ens_master+mod(proc_name+i,ens_images),40000+mobile_tag,    &
             MPI_COMM_WORLD,recv(i),stat)
    call mpi_isend(in,1,     &
             mpi_precision,ens_master+mod(proc_name+ens_images-i,ens_images),40000+mobile_tag,    &
             MPI_COMM_WORLD,send(i),stat)
   end do
   
  
   do i=0,ens_images-1
    call mpi_wait(recv(i),MPI_STATUS_IGNORE,stat)
    out=out+out_tmp(i)
   end do
   
   call mpi_waitall(ens_images,send,MPI_STATUSES_IGNORE,stat)
 
 
  
  end subroutine
  
  
  
  subroutine real_allmax(in, out)
   use parallel
   use sys
   
   implicit none
  
   real (kind=db), intent(in) :: in
   real (kind=db), intent(out) :: out
   real (kind=db) :: out_tmp(0:ens_images-1)
   integer :: send(0:ens_images-1), recv(0:ens_images-1), i
   
   
   out = 0.0d0
   
   do i=0,ens_images-1
    call mpi_irecv(out_tmp(i),1,     &
             mpi_precision,ens_master+i,40000,    &
             MPI_COMM_WORLD,recv(i),stat)
    call mpi_isend(in,1,     &
             mpi_precision,ens_master+i,40000,    &
             MPI_COMM_WORLD,send(i),stat)
   end do
   
  
   do i=0,ens_images-1
    call mpi_wait(recv(i),MPI_STATUS_IGNORE,stat)
    out=max(out,out_tmp(i))
    call mpi_waitall(ens_images,send,MPI_STATUSES_IGNORE,stat)
   end do
   
   
  
  end subroutine
 
end module


