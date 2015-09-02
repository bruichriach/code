

#include "include.h"

module writeout_grid

 use params

 implicit none

 integer, allocatable :: hgrid_out(:,:,:),ugrid_out(:,:,:),vgrid_out(:,:,:),zgrid_out(:,:,:)
 
end module

 
MODULE writeout

 implicit none
 
 integer :: max_write_tag=0
 
 CONTAINS
 
 subroutine init_write_grids
  use writeout_grid
  use parallel
  
  implicit none
  
  integer :: k
  
  allocate(hgrid_out(mx,my,0:ens_images-1), &
      ugrid_out(mx+1,my,0:ens_images-1),  &
      vgrid_out(mx,my+1,0:ens_images-1),  &
      zgrid_out(mx+1,my+1,0:ens_images-1))
      
  do k=0,ens_images-1
   hgrid_out(:,:,k)=merge(1,0,(proc_grid(1:mx,1:my) == ens_master+k))
   ugrid_out(:,:,k)=merge(1,0,((proc_grid(1:mx+1,1:my) == ens_master+k).or.&
                           (proc_grid(0:mx,1:my) == ens_master+k)))
   vgrid_out(:,:,k)=merge(1,0,((proc_grid(1:mx,1:my+1) == ens_master+k).or.&
                           (proc_grid(1:mx,0:my) == ens_master+k)))
   zgrid_out(:,:,k)=merge(1,0,(((proc_grid(1:mx+1,1:my+1) == ens_master+k).or.&
                           (proc_grid(0:mx,1:my+1) == ens_master+k)).or. &
                           ((proc_grid(1:mx+1,0:my) == ens_master+k).or. &
                           (proc_grid(0:mx,0:my) == ens_master+k))))
  end do
   
 
 end subroutine
 
 subroutine fullwrite_var(folder,dat)
  use writeout_grid
  use variables
  use parallel
  
  implicit none
 
  class(var), intent(inout) :: dat
  character(*) :: folder
  character(32) :: filename
  logical :: dir_e
  
  
   call mpi_barrier(mpi_comm_world,stat)
   call init_var_write(dat)
   if (proc_name == ens_master) then
    inquire( file=folder, exist=dir_e )
    do while (.not.(dir_e))   
     call system('mkdir -p '//folder,stat)
     inquire( file=folder, exist=dir_e )
    end do
    write (filename, "(a1,i4.4)") '/', ens_name
    inquire( file=folder//filename, exist=dir_e )
    do while (.not.(dir_e))   
     call system('mkdir -p '//folder//filename,stat)
     inquire( file=folder//filename, exist=dir_e )
    end do
   end if
   write (filename, "(a32)") folder//filename
   call mpi_barrier(mpi_comm_world,stat)
   call end_var_write(dat,adjustl(filename),0)
   call mpi_barrier(mpi_comm_world,stat)
  
 end subroutine 
 
 subroutine read_var(dat,folder)
  use writeout_grid
  use parallel
  use variables
  use grid
  
  implicit none
  
  class(var), intent(inout) :: dat
  character(*), intent(in) :: folder
  character(32) :: foldername, fullname
  integer :: i,j,k,l
  logical :: file_exist
  
   
  
  
  write (foldername, "(a28,i4.4)") trim(folder)//'/', ens_name
  write (fullname, "(a32)") trim(adjustl(foldername))//'/'//adjustl(dat%out%name)
  inquire( file=adjustl(fullname), exist=file_exist)
  
  if (file_exist) then
   
  

  if (proc_name == ens_master) then
   open(unit=10,form='unformatted',file=adjustl(fullname),status='unknown')
   write (format,"(a2,i3,a8)") '( ', mx+dat%p%ox, 'd24.16 )'
   read(10) ((dat%out%z(i,j),i=1,mx+dat%p%ox),j=1,my+dat%p%oy)
   close(10)
    
  
  
  
   do k=0,ens_images-1 
    l=0
    do j=1,ubound(dat%out%z,2)
     do i=1,ubound(dat%out%z,1)
      if (dat%out%grid(i,j,k) == 1) then
       l=l+1
       dat%out%recv(k)%dat(l)=dat%out%z(i,j)
      end if
     end do
    end do
    
    call mpi_isend(dat%out%recv(k)%dat(1),ubound(dat%out%recv(k)%dat),mpi_precision,  &
       ens_master+k,dat%out%tag,  &
       mpi_comm_world,dat%out%req_master(k),stat)
   end do
  end if
  call mpi_irecv(dat%out%send(1),ubound(dat%out%send),mpi_precision,  &
       ens_master,dat%out%tag,  &
       mpi_comm_world,dat%out%req,stat)
       
  
  k=0
  call mpi_wait(dat%out%req,mpi_status_ignore,stat)
  do j=dat%p%ly+1,dat%p%ly+dat%p%ny
   do i=dat%p%lx+1,dat%p%lx+dat%p%nx
    k=k+1
    dat%bz(i,j)=dat%out%send(k)
   end do
  end do
  
  if (proc_name == ens_master) then
   call mpi_waitall(ens_images, dat%out%req_master,mpi_statuses_ignore,stat)
   print *, 'read in complete for ', dat%out%name, ' on image ', ens_name
  end if
  
  
  end if
  

  
 end subroutine
 
 subroutine write_main(folder)
  use writeout_grid
  use parallel
  use variables
 
  implicit none
 
  character(*), intent(in) :: folder
  character(32) :: filename
  integer :: k
  logical :: dir_e
  
  if (proc_name == ens_master) then
   write (filename, "(a1,i4.4)") '/', ens_name
   inquire( file=folder//filename, exist=dir_e )
   do while (.not.(dir_e))   
    call system('mkdir -p '//folder//filename,stat)
    inquire( file=folder//filename, exist=dir_e )
   end do
  end if
  
  
  do k=1,nz
   call init_var_write(h(k))
   call init_var_write(u(k))
   call init_var_write(v(k))
  end do
#ifdef ALLOW_RIGID_LID
  call init_var_write(pres)
#endif
  write (filename, "(a1,i4.4)") '/', ens_name
  do k=1,nz
   call end_var_write(h(k),folder//filename,0)
   call end_var_write(u(k),folder//filename,0)
   call end_var_write(v(k),folder//filename,0)
  end do
#ifdef ALLOW_RIGID_LID
  call end_var_write(pres,folder//filename,0)
#endif
  
 end subroutine 
  
 subroutine init_writeouts(dat,grid)
  use writeout_grid
  use parallel
  use global
  
  implicit none
 
  class(var), intent(inout) :: dat
  integer, target, intent(in) :: grid(:,:,0:)
  integer :: k
  
  max_write_tag=max_write_tag+1
  
  allocate(dat%out%send(dat%p%nx*dat%p%ny))
  dat%out%req=mpi_request_null
  if (proc_name == ens_master) then
   dat%out%grid => grid
   allocate(dat%out%z(ubound(grid,1),ubound(grid,2)))
   dat%out%z(1:ubound(grid,1),1:ubound(grid,2))=0
   allocate(dat%out%req_master(0:ens_images-1))
   dat%out%req_master(0:ens_images-1)=mpi_request_null
   allocate(dat%out%recv(0:ens_images-1))
   do k=0,ens_images-1
    allocate(dat%out%recv(k)%dat(count((grid(:,:,k) == 1))))
   end do
  end if
  dat%out%tag=10000+max_write_tag
   
  
  
 end subroutine
 
  
 
 subroutine init_var_write(dat)
  use writeout_grid
  use variables
  use parallel

  IMPLICIT NONE

  class(var), intent(inout) :: dat
  integer :: i,j,k

  
  
  k=0
  call mpi_wait(dat%out%req,mpi_status_ignore,stat)
  do j=dat%p%ly+1,dat%p%ly+dat%p%ny
   do i=dat%p%lx+1,dat%p%lx+dat%p%nx
    k=k+1
    dat%out%send(k)=dat%bz(i,j)
   end do
  end do
  
  if (proc_name == ens_master) then
   do k=0,ens_images-1
    call mpi_irecv(dat%out%recv(k)%dat(1),ubound(dat%out%recv(k)%dat),mpi_precision,  &
       ens_master+k,dat%out%tag,  &
       mpi_comm_world,dat%out%req_master(k),stat)
   end do
  end if
  call mpi_isend(dat%out%send(1),ubound(dat%out%send),mpi_precision,  &
       ens_master,dat%out%tag,  &
       mpi_comm_world,dat%out%req,stat)
  
 end subroutine
 
 subroutine init_sum_write(dat)
  use writeout_grid
  use global
  use variables
 
  implicit none
  
  class(var) :: dat
   
  dat%s%bz=dat%s%bz/dble(dat%counter)
  call init_var_write(dat%s)
  dat%s%bz=0.0d0
  
 end subroutine
 
 subroutine end_var_write(dat,foldername,fmt)
  use writeout_grid
  use parallel
  use global
  
 
  implicit none
  
  character(*), intent(in) :: foldername
  character(32) :: filename
  class(var), intent(inout) :: dat
  integer :: i,j,k,l
  integer, intent(in) :: fmt
  
  
  
  
  write (filename, "(a32)") trim(foldername)//'/'//trim(dat%out%name)
  filename=adjustl(filename)
  
  if (proc_name == ens_master) then
   do k=0,ens_images-1 
    call mpi_wait(dat%out%req_master(k),mpi_status_ignore,stat)
    l=0
    do j=1,ubound(dat%out%z,2)
     do i=1,ubound(dat%out%z,1)
      if (dat%out%grid(i,j,k) == 1) then
       l=l+1
       dat%out%z(i,j)=dat%out%recv(k)%dat(l)
      end if
     end do
    end do
   end do
   call datawrite(fmt,filename,dat%out%z)
  end if
  
 end subroutine
 
 subroutine datawrite(type,filename, variable)
  use sys

  IMPLICIT NONE
 
  integer i, j
  integer :: type
  character(*), intent(in) :: filename
  real (kind=db), intent(in), dimension(:,:) :: variable
  
  write (format,"(a2,i3,a8)") '( ', ubound(variable,1), 'd24.16 )'

  if (type == 1) then
   open(unit=10,file=trim(filename),status='unknown')
     do j=1,ubound(variable,2)
      write(10,trim(format)) (variable(i,j),i=1,ubound(variable,1))
     end do
   close(10)
  else
   stat=-1
   do while (stat /= 0)
    open(unit=10,file=trim(filename),form='unformatted',status='unknown')
     write(10,iostat=stat) ((variable(i,j),i=1,ubound(variable,1)),j=1,ubound(variable,2))
    close(10)
   end do
  end if

  return

 end subroutine
 
 
 
 
 subroutine write()
  use writeout_grid
  use global
  use variables
  use parallel
 
  implicit none
   
  integer :: k
  
  
  
  call init_var_write(utau)
  call init_var_write(vtau)
  call init_var_write(bfricu)
  call init_var_write(bfricv)
#ifdef ALLOW_RIGID_LID
  call init_var_write(pres)
#endif
  call init_var_write(d(0))
  do k=1,nz
   call init_var_write(d(k))
   call init_var_write(h(k))
   call init_var_write(u(k))
   call init_var_write(v(k))
   call init_var_write(h(k)%tend_out)
   call init_var_write(u(k)%tend_out)
   call init_var_write(v(k)%tend_out)
   call init_var_write(smagu(k))
   call init_var_write(smagv(k))
   call init_var_write(ke(k))
   call init_var_write(ape(k))
   call init_var_write(q(k))
   call init_var_write(zeta(k))
  end do
  
#ifdef DO_TIME_AVERAGE
#ifdef DO_SHORT_AVERAGE
   do k=1,nz
    call init_sum_write(h(k))
   end do
   do k=1,nz
    call init_sum_write(h_u(k))
   end do
   do k=1,nz
    call init_sum_write(h_v(k))
   end do
   do k=1,nz
    call init_sum_write(h_z(k))
   end do
   
   do k=1,nz
    call init_sum_write(hu_h(k))
   end do
   do k=1,nz
    call init_sum_write(hu(k))
   end do
   do k=1,nz
    call init_sum_write(hu_v(k))
   end do
   do k=1,nz
    call init_sum_write(hu_z(k))
   end do
   
   do k=1,nz
    call init_sum_write(hv_h(k))
   end do
   do k=1,nz
    call init_sum_write(hv_u(k))
   end do
   do k=1,nz
    call init_sum_write(hv(k))
   end do
   do k=1,nz
    call init_sum_write(hv_z(k))
   end do

   do k=1,nz
    call init_sum_write(u_h(k))
   end do
   do k=1,nz
    call init_sum_write(u(k))
   end do
   do k=1,nz
    call init_sum_write(u_v(k))
   end do
   do k=1,nz
    call init_sum_write(u_z(k))
   end do

   do k=1,nz
    call init_sum_write(v_h(k))
   end do
   do k=1,nz
    call init_sum_write(v_u(k))
   end do
   do k=1,nz
    call init_sum_write(v(k))
   end do
   do k=1,nz
    call init_sum_write(v_z(k))
   end do

   
   do k=1,nz
    call init_sum_write(huu_h(k))
   end do

   do k=1,nz
    call init_sum_write(huu_z(k))
   end do
   
   do k=1,nz
    call init_sum_write(hvv_h(k))
   end do

   do k=1,nz
    call init_sum_write(hvv_z(k))
   end do
   
   do k=1,nz
    call init_sum_write(huv_h(k))
   end do
   
   do k=1,nz
    call init_sum_write(huv_z(k))
   end do
   
   do k=1,nz
    call init_sum_write(huv_uz(k))
   end do
   
   do k=1,nz
    call init_sum_write(huv_vz(k))
   end do
   
   do k=1,nz
    call init_sum_write(huuu(k))
   end do

   do k=1,nz
    call init_sum_write(huuv(k))
   end do
   
   do k=1,nz
    call init_sum_write(huvv(k))
   end do

   do k=1,nz
    call init_sum_write(hvvv(k))
   end do
   
   do k=1,nz
    call init_sum_write(hum(k))
   end do

   do k=1,nz
    call init_sum_write(hvm(k))
   end do
      
   do k=1,nz
    call init_sum_write(hm_x_u(k))
   end do
   
   do k=1,nz
    call init_sum_write(hm_y_v(k))
   end do
   
   do k=1,nz
    call init_sum_write(up_dm_x_u(k))
   end do
   
   do k=1,nz
    call init_sum_write(up_dm_y_v(k))
   end do
   
   do k=1,nz
    call init_sum_write(up_du_u(k))
   end do
   
   do k=1,nz
    call init_sum_write(up_dv_v(k))
   end do
   
   do k=1,nz
    call init_sum_write(dn_dm_x_u(k))
   end do
   
   do k=1,nz
    call init_sum_write(dn_dm_y_v(k))
   end do
   
   do k=1,nz
    call init_sum_write(dn_du_u(k))
   end do
   
   do k=1,nz
    call init_sum_write(dn_dv_v(k))
   end do
   
   do k=0,nz
    call init_sum_write(dd(k))
   end do
   
   do k=0,nz
    call init_sum_write(d(k))
   end do
      
   do k=1,nz
    call init_sum_write(m(k))
   end do
   
   do k=1,nz
    call init_sum_write(hm(k))
   end do
 
   call init_sum_write(utau)
   call init_sum_write(vtau)
 
   call init_sum_write(bfricu)
   call init_sum_write(bfricv)
   
   do k=1,nz
    call init_sum_write(hsmagu(k))
   end do
   
   do k=1,nz
    call init_sum_write(hsmagv(k))
   end do
 
   call init_sum_write(uutau)
   call init_sum_write(vvtau)
 
   call init_sum_write(ubfricu)
   call init_sum_write(vbfricv)
   
   do k=1,nz
    call init_sum_write(husmagu(k))
   end do
   
   do k=1,nz
    call init_sum_write(hvsmagv(k))
   end do
   
   
   do k=1,nz
    call init_sum_write(hq_h(k))
   end do
   do k=1,nz
    call init_sum_write(hq_u(k))
   end do
   do k=1,nz
    call init_sum_write(hq_v(k))
   end do
   do k=1,nz
    call init_sum_write(hq(k))
   end do
   
   do k=1,nz
    call init_sum_write(huq_h(k))
   end do
   do k=1,nz
    call init_sum_write(huq(k))
   end do
   do k=1,nz
    call init_sum_write(huq_v(k))
   end do
   do k=1,nz
    call init_sum_write(huq_z(k))
   end do
   
   do k=1,nz
    call init_sum_write(hvq_h(k))
   end do
   do k=1,nz
    call init_sum_write(hvq_u(k))
   end do
   do k=1,nz
    call init_sum_write(hvq(k))
   end do
   do k=1,nz
    call init_sum_write(hvq_z(k))
   end do
 
   do k=1,nz
    call init_sum_write(hqq(k))
   end do
   
   do k=1,nz
    call init_sum_write(q(k))
   end do
   do k=1,nz
    call init_sum_write(q_u(k))
   end do
   do k=1,nz
    call init_sum_write(q_v(k))
   end do
   do k=1,nz
    call init_sum_write(q(k))
   end do
   
   do k=1,nz
    call init_sum_write(uq_h(k))
   end do
   do k=1,nz
    call init_sum_write(uq(k))
   end do
   do k=1,nz
    call init_sum_write(uq_v(k))
   end do
   do k=1,nz
    call init_sum_write(uq_z(k))
   end do
   
   do k=1,nz
    call init_sum_write(vq_h(k))
   end do
   do k=1,nz
    call init_sum_write(vq_u(k))
   end do
   do k=1,nz
    call init_sum_write(vq(k))
   end do
   do k=1,nz
    call init_sum_write(vq_z(k))
   end do
 
   do k=1,nz
    call init_sum_write(qq(k))
   end do

   do k=1,nz
    call init_sum_write(tendh(k))
   end do

   do k=1,nz
    call init_sum_write(htendu(k))
   end do

   do k=1,nz
    call init_sum_write(htendv(k))
   end do



   do k=1,nz
    call init_sum_write(utendh(k))
   end do

   do k=1,nz
    call init_sum_write(vtendh(k))
   end do
#endif
#endif
  
  
 
 
 end subroutine
 
 subroutine do_write(num)
  use parallel
  use variables
 
  implicit none
  
  integer, intent(in) :: num
  character(32) :: filename, desc
  integer :: i, k
  logical :: dir_e
  
  
  if (proc_name == ens_master) then
   write (filename, "(a7,i4.4,a2)") './data/', ens_name, '/.'
   inquire( file=filename, exist=dir_e )
   do while (.not.(dir_e))   
    write (filename, "(a14,i4.4)") 'mkdir -p data/', ens_name
    call system(filename,stat)
    write (filename, "(a7,i4.4,a2)") './data/', ens_name, '/.'
    inquire( file=filename, exist=dir_e )
   end do
   write (filename, "(a7,i4.4,a1,i4.4,a2)") './data/', ens_name, '/', num, '/.'
   inquire( file=filename, exist=dir_e )
   do while (.not.(dir_e))   
    write (filename, "(a14,i4.4,a1,i4.4)") 'mkdir -p data/', ens_name, '/', num
    call system(filename)
    write (filename, "(a7,i4.4,a1,i4.4,a2)") './data/', ens_name, '/', num, '/.'
    inquire( file=filename, exist=dir_e )
   end do
   
   write (filename, "(a5,i4.4,a1,i4.4)") 'data/', ens_name, '/', num
   call end_var_write(utau,filename,0)
   call end_var_write(vtau,filename,0)
    call end_var_write(bfricu,filename,0)
    call end_var_write(bfricv,filename,0)
#ifdef ALLOW_RIGID_LID
   call end_var_write(pres,filename,0)
#endif
   call end_var_write(d(0),filename,0)   
   do i=1,nz   
    call end_var_write(d(i),filename,0)   
    call end_var_write(h(i),filename,0)
    call end_var_write(u(i),filename,0)
    call end_var_write(v(i),filename,0)
    call end_var_write(h(i)%tend_out,filename,0)
    call end_var_write(u(i)%tend_out,filename,0)
    call end_var_write(v(i)%tend_out,filename,0)
    call end_var_write(smagu(i),filename,0)
    call end_var_write(smagv(i),filename,0)
    call end_var_write(ke(i),filename,0)
    call end_var_write(ape(i),filename,0)
    call end_var_write(q(i),filename,0)
    call end_var_write(zeta(i),filename,0)
   
   end do
   
   
  
#ifdef DO_TIME_AVERAGE
#ifdef DO_SHORT_AVERAGE
   do k=1,nz
    call end_var_write(h(k)%s,filename,0)
   end do
   do k=1,nz
    call end_var_write(h_u(k)%s,filename,0)
   end do
   do k=1,nz
    call end_var_write(h_v(k)%s,filename,0)
   end do
   do k=1,nz
    call end_var_write(h_z(k)%s,filename,0)
   end do
   
   do k=1,nz
    call end_var_write(hu_h(k)%s,filename,0)
   end do
   do k=1,nz
    call end_var_write(hu(k)%s,filename,0)
   end do
   do k=1,nz
    call end_var_write(hu_v(k)%s,filename,0)
   end do
   do k=1,nz
    call end_var_write(hu_z(k)%s,filename,0)
   end do
   
   do k=1,nz
    call end_var_write(hv_h(k)%s,filename,0)
   end do
   do k=1,nz
    call end_var_write(hv_u(k)%s,filename,0)
   end do
   do k=1,nz
    call end_var_write(hv(k)%s,filename,0)
   end do
   do k=1,nz
    call end_var_write(hv_z(k)%s,filename,0)
   end do

   do k=1,nz
    call end_var_write(u_h(k)%s,filename,0)
   end do
   do k=1,nz
    call end_var_write(u(k)%s,filename,0)
   end do
   do k=1,nz
    call end_var_write(u_v(k)%s,filename,0)
   end do
   do k=1,nz
    call end_var_write(u_z(k)%s,filename,0)
   end do

   do k=1,nz
    call end_var_write(v_h(k)%s,filename,0)
   end do
   do k=1,nz
    call end_var_write(v_u(k)%s,filename,0)
   end do
   do k=1,nz
    call end_var_write(v(k)%s,filename,0)
   end do
   do k=1,nz
    call end_var_write(v_z(k)%s,filename,0)
   end do

   
   do k=1,nz
    call end_var_write(huu_h(k)%s,filename,0)
   end do

   do k=1,nz
    call end_var_write(huu_z(k)%s,filename,0)
   end do
   
   do k=1,nz
    call end_var_write(hvv_h(k)%s,filename,0)
   end do

   do k=1,nz
    call end_var_write(hvv_z(k)%s,filename,0)
   end do
   
   do k=1,nz
    call end_var_write(huv_h(k)%s,filename,0)
   end do
   
   do k=1,nz
    call end_var_write(huv_z(k)%s,filename,0)
   end do
   
   do k=1,nz
    call end_var_write(huv_uz(k)%s,filename,0)
   end do
   
   do k=1,nz
    call end_var_write(huv_vz(k)%s,filename,0)
   end do
   
   do k=1,nz
    call end_var_write(huuu(k)%s,filename,0)
   end do

   do k=1,nz
    call end_var_write(huuv(k)%s,filename,0)
   end do
   
   do k=1,nz
    call end_var_write(huvv(k)%s,filename,0)
   end do

   do k=1,nz
    call end_var_write(hvvv(k)%s,filename,0)
   end do
   
   do k=1,nz
    call end_var_write(hum(k)%s,filename,0)
   end do

   do k=1,nz
    call end_var_write(hvm(k)%s,filename,0)
   end do
      
   do k=1,nz
    call end_var_write(hm_x_u(k)%s,filename,0)
   end do
   
   do k=1,nz
    call end_var_write(hm_y_v(k)%s,filename,0)
   end do
   
   do k=1,nz
    call end_var_write(up_dm_x_u(k)%s,filename,0)
   end do
   
   do k=1,nz
    call end_var_write(up_dm_y_v(k)%s,filename,0)
   end do
   
   do k=1,nz
    call end_var_write(up_du_u(k)%s,filename,0)
   end do
   
   do k=1,nz
    call end_var_write(up_dv_v(k)%s,filename,0)
   end do
   
   do k=1,nz
    call end_var_write(dn_dm_x_u(k)%s,filename,0)
   end do
   
   do k=1,nz
    call end_var_write(dn_dm_y_v(k)%s,filename,0)
   end do
   
   do k=1,nz
    call end_var_write(dn_du_u(k)%s,filename,0)
   end do
   
   do k=1,nz
    call end_var_write(dn_dv_v(k)%s,filename,0)
   end do
   
   do k=0,nz
    call end_var_write(dd(k)%s,filename,0)
   end do
   
   do k=0,nz
    call end_var_write(d(k)%s,filename,0)
   end do
      
   do k=1,nz
    call end_var_write(m(k)%s,filename,0)
   end do
   
   do k=1,nz
    call end_var_write(hm(k)%s,filename,0)
   end do
 
   call end_var_write(utau%s,filename,0)
   call end_var_write(vtau%s,filename,0)
 
   call end_var_write(bfricu%s,filename,0)
   call end_var_write(bfricv%s,filename,0)
   
   do k=1,nz
    call end_var_write(hsmagu(k)%s,filename,0)
   end do
   
   do k=1,nz
    call end_var_write(hsmagv(k)%s,filename,0)
   end do
 
   call end_var_write(uutau%s,filename,0)
   call end_var_write(vvtau%s,filename,0)
 
   call end_var_write(ubfricu%s,filename,0)
   call end_var_write(vbfricv%s,filename,0)
   
   do k=1,nz
    call end_var_write(husmagu(k)%s,filename,0)
   end do
   
   do k=1,nz
    call end_var_write(hvsmagv(k)%s,filename,0)
   end do
   
   do k=1,nz
    call end_var_write(hq_h(k)%s,filename,0)
   end do
   do k=1,nz
    call end_var_write(hq_u(k)%s,filename,0)
   end do
   do k=1,nz
    call end_var_write(hq_v(k)%s,filename,0)
   end do
   do k=1,nz
    call end_var_write(hq(k)%s,filename,0)
   end do
   
   do k=1,nz
    call end_var_write(huq_h(k)%s,filename,0)
   end do
   do k=1,nz
    call end_var_write(huq(k)%s,filename,0)
   end do
   do k=1,nz
    call end_var_write(huq_v(k)%s,filename,0)
   end do
   do k=1,nz
    call end_var_write(huq_z(k)%s,filename,0)
   end do
   
   do k=1,nz
    call end_var_write(hvq_h(k)%s,filename,0)
   end do
   do k=1,nz
    call end_var_write(hvq_u(k)%s,filename,0)
   end do
   do k=1,nz
    call end_var_write(hvq(k)%s,filename,0)
   end do
   do k=1,nz
    call end_var_write(hvq_z(k)%s,filename,0)
   end do
 
   do k=1,nz
    call end_var_write(hqq(k)%s,filename,0)
   end do
   
   do k=1,nz
    call end_var_write(q(k)%s,filename,0)
   end do
   do k=1,nz
    call end_var_write(q_u(k)%s,filename,0)
   end do
   do k=1,nz
    call end_var_write(q_v(k)%s,filename,0)
   end do
   do k=1,nz
    call end_var_write(q(k)%s,filename,0)
   end do
   
   do k=1,nz
    call end_var_write(uq_h(k)%s,filename,0)
   end do
   do k=1,nz
    call end_var_write(uq(k)%s,filename,0)
   end do
   do k=1,nz
    call end_var_write(uq_v(k)%s,filename,0)
   end do
   do k=1,nz
    call end_var_write(uq_z(k)%s,filename,0)
   end do
   
   do k=1,nz
    call end_var_write(vq_h(k)%s,filename,0)
   end do
   do k=1,nz
    call end_var_write(vq_u(k)%s,filename,0)
   end do
   do k=1,nz
    call end_var_write(vq(k)%s,filename,0)
   end do
   do k=1,nz
    call end_var_write(vq_z(k)%s,filename,0)
   end do
 
   do k=1,nz
    call end_var_write(qq(k)%s,filename,0)
   end do

   do k=1,nz
    call end_var_write(tendh(k)%s,filename,0)
   end do

   do k=1,nz
    call end_var_write(htendu(k)%s,filename,0)
   end do

   do k=1,nz
    call end_var_write(htendv(k)%s,filename,0)
   end do



   do k=1,nz
    call end_var_write(utendh(k)%s,filename,0)
   end do

   do k=1,nz
    call end_var_write(vtendh(k)%s,filename,0)
   end do
#endif
#endif
   
   open(unit=10,file=trim(filename)//'/params.txt',status='unknown',form='formatted')
   
   write(desc,"(a32)") 'wind energy flux'
   format="(a32,a1,e23.16)"
   write(10,format) adjustl(desc), ':', sum(0.5d0*(utau%out%z(1:mx,1:my)*u(3)%out%z(1:mx,1:my)+ &
              utau%out%z(2:mx+1,1:my)*u(3)%out%z(2:mx+1,1:my))+  &
              0.5d0*(vtau%out%z(1:mx,1:my)*v(3)%out%z(1:mx,1:my)+  &
              vtau%out%z(1:mx,2:my+1)*v(3)%out%z(1:mx,2:my+1)))
   
   write(desc,"(a32)") 'bottom friction energy flux'
   format="(a32,a1,e23.16)"
   write(10,format) adjustl(desc), ':', -sum(0.5d0*(bfricu%out%z(1:mx,1:my)*u(1)%out%z(1:mx,1:my)+ &
              bfricu%out%z(2:mx+1,1:my)*u(1)%out%z(2:mx+1,1:my))+  &
              0.5d0*(bfricv%out%z(1:mx,1:my)*v(1)%out%z(1:mx,1:my)+  &
              bfricv%out%z(1:mx,2:my+1)*v(1)%out%z(1:mx,2:my+1)))
   
   do i=1,nz
    write(desc,"(a31,i1)") 'smagorisky energy flux, layer ', i
    format="(a32,a1,e23.16)"
    write(10,format) adjustl(desc), ':', sum(h(i)%out%z(1:mx,1:my)*(0.5d0*(smagu(i)%out%z(1:mx,1:my)*u(i)%out%z(1:mx,1:my)+ &
              smagu(i)%out%z(2:mx+1,1:my)*u(i)%out%z(2:mx+1,1:my))+  &
              0.5d0*(smagv(i)%out%z(1:mx,1:my)*v(i)%out%z(1:mx,1:my)+  &
              smagv(i)%out%z(1:mx,2:my+1)*v(i)%out%z(1:mx,2:my+1))))
   end do

   write(desc,"(a32)") 'minimum layer thickness'
   format="(a32,a1,e23.16)"
   write(10,format) adjustl(desc), ':', minval( (/ (minval(h(i)%out%z(1:mx,1:my)) , i=1,nz) /) )
              
   close(10)
  end if
  
 end subroutine

END MODULE
  
