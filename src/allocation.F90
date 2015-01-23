module allocation
   use grid
   
   implicit none

  
  interface create_field
   module procedure create_hvar, create_uvar,  &
         create_vvar, create_zvar
   module procedure create_written_hvar, create_written_uvar,  &
         create_written_vvar, create_written_zvar
   module procedure create_hvar_layer, create_uvar_layer,  &
         create_vvar_layer, create_zvar_layer
   module procedure create_written_hvar_layer, create_written_uvar_layer,  &
         create_written_vvar_layer, create_written_zvar_layer
  end interface
  
  contains

  
  subroutine create_var(dat,grid,synced)
   use global
   use params
   use parallel
   use writeout
   
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
  
  
  
  
  
  
  subroutine create_written_hvar(dat,name,synced)
   use global
   use writeout_grid
   use writeout
   
   implicit none
   
   type(hvar), intent(out) :: dat
   logical, intent(in) :: synced
   character(*), intent(in) :: name
   
      
   call create_var(dat,hgrid,synced)
   write (dat%out%name, "(a16)") trim(name)//".dat"
   dat%out%name=adjustl(dat%out%name)
   call init_writeouts(dat,hgrid_out)
   
  end subroutine
  
  
  subroutine create_written_uvar(dat,name,synced)
   use global
   use writeout_grid
   use writeout
   
   implicit none
   
   type(uvar), intent(out) :: dat
   logical, intent(in) :: synced
   character(*), intent(in) :: name
   
   
   
   call create_var(dat,ugrid,synced)
   write (dat%out%name, "(a8)") name//".dat"
   dat%out%name=adjustl(dat%out%name)
   call init_writeouts(dat,ugrid_out)
   
  end subroutine
  
  
  subroutine create_written_vvar(dat,name,synced)
   use global
   use writeout_grid
   use writeout
   
   implicit none
   
   type(vvar), intent(out) :: dat
   logical, intent(in) :: synced
   character(*), intent(in) :: name
   
   
   
   call create_var(dat,vgrid,synced)
   write (dat%out%name, "(a8)") name//".dat"
   dat%out%name=adjustl(dat%out%name)
   call init_writeouts(dat,vgrid_out)
   
  end subroutine
  
  
  subroutine create_written_zvar(dat,name,synced)
   use global
   use writeout_grid
   use writeout
   
   implicit none
   
   type(zvar), intent(inout) :: dat
   logical, intent(in) :: synced
   character(*), intent(in) :: name
   
   
   
   call create_var(dat,zgrid,synced)
   write (dat%out%name, "(a8)") name//".dat"
   dat%out%name=adjustl(dat%out%name)
   call init_writeouts(dat,zgrid_out)
   
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
  
  
  
  
  
  
  
  
  subroutine create_written_hvar_layer(dat,name,synced)
   use global
   use writeout_grid
   use writeout
   
   implicit none
   
   type(hvar), intent(inout) :: dat(:)
   logical, intent(in) :: synced
   integer :: k
   character(*), intent(in) :: name
   
   
   
   do k=lbound(dat,1),ubound(dat,1)
    call create_hvar(dat(k),synced)
    write (dat(k)%out%name, "(a3,i1,a4)") name//"_",k,".dat"
    dat(k)%out%name=adjustl(dat(k)%out%name)
    call init_writeouts(dat(k),hgrid_out)
   end do
   
  end subroutine
  
  
  subroutine create_written_uvar_layer(dat,name,synced)
   use global
   use writeout_grid
   use writeout
   
   implicit none
   
   type(uvar), intent(inout) :: dat(:)
   logical, intent(in) :: synced
   integer :: k
   character(*), intent(in) :: name
   
   
   
   do k=lbound(dat,1),ubound(dat,1)
    call create_uvar(dat(k),synced)
    write (dat(k)%out%name, "(a3,i1,a4)") name//"_",k,".dat"
    dat(k)%out%name=adjustl(dat(k)%out%name)
    call init_writeouts(dat(k),ugrid_out)
   end do
   
  end subroutine
  
  
  subroutine create_written_vvar_layer(dat,name,synced)
   use global
   use writeout_grid
   use writeout
   
   implicit none
   
   type(vvar), intent(inout) :: dat(:)
   logical, intent(in) :: synced
   integer :: k
   character(*), intent(in) :: name
   
   
   do k=lbound(dat,1),ubound(dat,1)
    call create_vvar(dat(k),synced)
    write (dat(k)%out%name, "(a3,i1,a4)") name//"_",k,".dat"
    dat(k)%out%name=adjustl(dat(k)%out%name)
    call init_writeouts(dat(k),vgrid_out)
   end do
   
  end subroutine
  
  
  subroutine create_written_zvar_layer(dat,name,synced)
   use global
   use writeout_grid
   use writeout
   
   implicit none
   
   type(zvar), intent(inout) :: dat(:)
   logical, intent(in) :: synced
   integer :: k
   character(*), intent(in) :: name
   
   
   
   do k=lbound(dat,1),ubound(dat,1)
    call create_zvar(dat(k),synced)
    write (dat(k)%out%name, "(a3,i1,a4)") name//"_",k,".dat"
    dat(k)%out%name=adjustl(dat(k)%out%name)
    call init_writeouts(dat(k),zgrid_out)
   end do
   
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
