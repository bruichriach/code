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
   hgrid_out(:,:,k)=merge(1,0,(proc_grid(1:mx,1:my) == k))
   ugrid_out(:,:,k)=merge(1,0,((proc_grid(1:mx+1,1:my) == k).or.&
                           (proc_grid(0:mx,1:my) == k)))
   vgrid_out(:,:,k)=merge(1,0,((proc_grid(1:mx,1:my+1) == k).or.&
                           (proc_grid(1:mx,0:my) == k)))
   zgrid_out(:,:,k)=merge(1,0,(((proc_grid(1:mx+1,1:my+1) == k).or.&
                           (proc_grid(0:mx,1:my+1) == k)).or. &
                           ((proc_grid(1:mx+1,0:my) == k).or. &
                           (proc_grid(0:mx,0:my) == k))))
  end do
   
 
 end subroutine
 
 subroutine write_config(dat)
  use writeout_grid
  use variables
  use parallel
  
  implicit none
 
  type(var), intent(inout) :: dat
  character(32) :: filename
  logical :: dir_e
  
  
   call mpi_barrier(mpi_comm_world,stat)
   call init_var_write(dat)
   if (proc_name == ens_master) then
    write (filename, "(a1,i4.4)") '/', ens_name
    inquire( file='global'//filename, exist=dir_e )
    do while (.not.(dir_e))   
     call system('mkdir -p '//'global'//filename,stat)
     inquire( file='global'//filename, exist=dir_e )
    end do
   end if
   write (filename, "(a14,i4.4)") 'global/', ens_name
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
  character(32) :: format, foldername, fullname
  integer :: i,j,k,l
  logical :: file_exist
  
   
  
  
  write (foldername, "(a28,i4.4)") trim(folder)//'/', ens_name
  write (fullname, "(a32)") trim(adjustl(foldername))//'/'//adjustl(dat%out%name)
  inquire( file=adjustl(fullname), exist=file_exist)
  
  if (file_exist) then
   
  

  if (proc_name == ens_master) then
   open(unit=10,form='unformatted',file=adjustl(fullname),status='unknown')
   write (format,"(a2,i3,a8)") '( ', mx+dat%p%ox, 'd24.16 )'
   do j=1,my+dat%p%oy
    read(10) (dat%out%z(i,j),i=1,mx+dat%p%ox)
   end do
   close(10)
    
  
  
   do k=0,ens_images-1 
    l=0
    do j=dat%p%ly+1,dat%p%ly+dat%p%ny
     do i=dat%p%lx+1,dat%p%lx+dat%p%nx
      l=l+1
      dat%out%recv(k)%dat(l)=dat%out%z(i,j)
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
    dat%z(i,j)=dat%out%send(k)
   end do
  end do
  
  if (proc_name == ens_master) then
   call mpi_waitall(ens_images, dat%out%req_master,mpi_statuses_ignore,stat)
   print *, 'read in complete for ', dat%out%name
  end if
  
  
  end if
  

  
 end subroutine
 
 subroutine write_main(folder)
  use writeout_grid
  use parallel
  use variables
  use solver_variables
 
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
  call init_var_write(pres)
  write (filename, "(a1,i4.4)") '/', ens_name
  do k=1,nz
   call end_var_write(h(k),folder//filename,0)
   call end_var_write(u(k),folder//filename,0)
   call end_var_write(v(k),folder//filename,0)
  end do
  call end_var_write(pres,folder//filename,0)
  
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

  IMPLICIT NONE
 
  integer i, j
  integer :: type
  character(*), intent(in) :: filename
  real (kind=8), intent(in), dimension(:,:) :: variable
  character(32) :: format
  
  write (format,"(a2,i3,a8)") '( ', ubound(variable,1), 'd24.16 )'

  if (type == 1) then
   open(unit=10,file=trim(filename),status='unknown')
     do j=1,ubound(variable,2)
      write(10,trim(format)) (variable(i,j),i=1,ubound(variable,1))
     end do
   close(10)
  else
   open(unit=10,file=trim(filename),form='unformatted',status='unknown')
     do j=1,ubound(variable,2)
      write(10) (variable(i,j),i=1,ubound(variable,1))
     end do
   close(10)
  end if

  return

 end subroutine
 
 
 
 
 subroutine write()
  use writeout_grid
  use global
  use solver_variables
  use variables
  use parallel
 
  implicit none
   
  integer :: i
  
  
  
  call init_var_write(utau)
  call init_var_write(vtau)
  call init_var_write(pres)
  call init_var_write(depth(0))
  do i=1,nz
   call init_var_write(depth(i))
   call init_var_write(h(i))
   call init_var_write(u(i))
   call init_var_write(v(i))
   call init_var_write(tendh1(i))
   call init_var_write(tendu1(i))
   call init_var_write(tendv1(i))
   call init_var_write(smagu(i))
   call init_var_write(smagv(i))
   call init_var_write(ke(i))
   call init_var_write(ape(i))
   call init_var_write(q(i))
   call init_var_write(zeta(i))
  end do
  
  
 
 
 end subroutine
 
 subroutine do_write(num)
  use parallel
  use variables
  use solver_variables
 
  implicit none
  
  integer, intent(in) :: num
  character(32) :: filename
  integer :: i
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
   call end_var_write(pres,filename,0)
   call end_var_write(depth(0),filename,0)   
   do i=1,nz   
    call end_var_write(depth(i),filename,0)   
    call end_var_write(h(i),filename,0)
    call end_var_write(u(i),filename,0)
    call end_var_write(v(i),filename,0)
    call end_var_write(tendh1(i),filename,0)
    call end_var_write(tendu1(i),filename,0)
    call end_var_write(tendv1(i),filename,0)
    call end_var_write(smagu(i),filename,0)
    call end_var_write(smagv(i),filename,0)
    call end_var_write(ke(i),filename,0)
    call end_var_write(ape(i),filename,0)
    call end_var_write(q(i),filename,0)
    call end_var_write(zeta(i),filename,0)
   
   end do
  end if
  
 end subroutine

END MODULE
  
