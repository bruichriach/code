module allocation

#include "include.h"

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
  interface make_timestep
   module procedure make_htimestep, make_utimestep, make_vtimestep
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
      
   if (synced) then
    allocate(dat%z(dat%p%lx:dat%p%lx+dat%p%nx+1,   &
         dat%p%ly:dat%p%ly+dat%p%ny+1))
    dat%bz =>  &
         dat%z(dat%p%lx+1:dat%p%lx+dat%p%nx,   &
         dat%p%ly+1:dat%p%ly+dat%p%ny)
    dat%bz(dat%p%lx+1:,dat%p%ly+1:) => dat%bz
   else
    allocate(dat%bz(dat%p%lx+1:dat%p%lx+dat%p%nx,   &
         dat%p%ly+1:dat%p%ly+dat%p%ny))
   end if
   dat%bz=0.0d0
   
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
     if (allocated(dat%mpi(j)%s)) then
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
      if (l /= ubound(dat%mpi(j)%s,1)) then
       print *, "OMG SEND ERROR!!!"
       stop
      end if
     end if
    end do
   
   
    do j=0,proc_num-1
     if (allocated(dat%mpi(j)%r)) then
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
      if (l /= ubound(dat%mpi(j)%r,1)) then
       print *, "OMG RECIEVE ERROR!!!"
       stop
      end if
     end if
    end do
    
   end if
    
  end subroutine
  
  
  
  
  subroutine make_htimestep(dat)
   use global
   use writeout
   use writeout_grid
   
   implicit none
   
   type(hvar), intent(inout) :: dat
   
   
   allocate(dat%tmp)
   allocate(dat%tend1)
   allocate(dat%tend2)
   allocate(dat%tend3)
   allocate(dat%tend_out)
   call create_var(dat%tmp,hgrid,.false.)
   call create_var(dat%tend1,hgrid,.false.)
   call create_var(dat%tend2,hgrid,.false.)
   call create_var(dat%tend3,hgrid,.false.)
   call create_var(dat%tend_out,hgrid,.false.)
   dat%tend_out%bz(dat%p%lx+1:,dat%p%ly+1:) => dat%tend1%bz
   call init_writeouts(dat%tend_out,hgrid_out)
   write (dat%tend_out%out%name, "(a16)") "tend"//adjustl(trim(dat%out%name))
   dat%tend_out%out%name=adjustl(dat%tend_out%out%name)
         
         
  end subroutine
  
  subroutine make_utimestep(dat)
   use global
   use writeout
   use writeout_grid
   
   implicit none
   
   type(uvar), intent(inout) :: dat
   
   
   allocate(dat%tmp)
   allocate(dat%tend1)
   allocate(dat%tend2)
   allocate(dat%tend3)
   allocate(dat%tend_out)
   call create_var(dat%tmp,ugrid,.false.)
   call create_var(dat%tend1,ugrid,.false.)
   call create_var(dat%tend2,ugrid,.false.)
   call create_var(dat%tend3,ugrid,.false.)
   call create_var(dat%tend_out,ugrid,.false.)
   dat%tend_out%bz(dat%p%lx+1:,dat%p%ly+1:) => dat%tend1%bz
   call init_writeouts(dat%tend_out,ugrid_out)
   write (dat%tend_out%out%name, "(a16)") "tend"//adjustl(trim(dat%out%name))
   dat%tend_out%out%name=adjustl(dat%tend_out%out%name)
         
         
  end subroutine
  
  subroutine make_vtimestep(dat)
   use global
   use writeout
   use writeout_grid
   
   implicit none
   
   type(vvar), intent(inout) :: dat
   
   
   allocate(dat%tmp)
   allocate(dat%tend1)
   allocate(dat%tend2)
   allocate(dat%tend3)
   allocate(dat%tend_out)
   call create_var(dat%tmp,vgrid,.false.)
   call create_var(dat%tend1,vgrid,.false.)
   call create_var(dat%tend2,vgrid,.false.)
   call create_var(dat%tend3,vgrid,.false.)
   call create_var(dat%tend_out,vgrid,.false.)
   dat%tend_out%bz(dat%p%lx+1:,dat%p%ly+1:) => dat%tend1%bz
   call init_writeouts(dat%tend_out,vgrid_out)
   write (dat%tend_out%out%name, "(a16)") "tend"//adjustl(trim(dat%out%name))
   dat%tend_out%out%name=adjustl(dat%tend_out%out%name)
         
         
  end subroutine
  
  
 
 
 elemental subroutine stepforward(dat)
  use global
  
  implicit none
  
  class(var), intent(inout) :: dat
  
  
   dat%tmp => dat%tend3
   dat%tend3 => dat%tend2
   dat%tend2 => dat%tend1
   dat%tend1 => dat%tmp
   dat%tend_out%bz(dat%p%lx+1:,dat%p%ly+1:) => dat%tend1%bz
  
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
   write (dat%out%name, "(a32)") adjustl(name)//".dat"
   dat%out%name=adjustl(dat%out%name)
   call init_writeouts(dat,hgrid_out)
   allocate(dat%t)
   call create_var(dat%t,hgrid,.false.)
   write (dat%t%out%name, "(a32)") "s_"//adjustl(name)//".dat"
   dat%t%out%name=adjustl(dat%t%out%name)
   call init_writeouts(dat%t,hgrid_out)
   allocate(dat%s)
   call create_var(dat%s,hgrid,.false.)
   write (dat%s%out%name, "(a32)") "s_"//adjustl(name)//".dat"
   dat%s%out%name=adjustl(dat%s%out%name)
   call init_writeouts(dat%s,hgrid_out)
   dat%counter=>counter
   
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
   write (dat%out%name, "(a32)") adjustl(name)//".dat"
   dat%out%name=adjustl(dat%out%name)
   call init_writeouts(dat,ugrid_out)
   allocate(dat%t)
   call create_var(dat%t,ugrid,.false.)
   write (dat%t%out%name, "(a32)") "s_"//adjustl(name)//".dat"
   dat%t%out%name=adjustl(dat%t%out%name)
   call init_writeouts(dat%t,ugrid_out)
   allocate(dat%s)
   call create_var(dat%s,ugrid,.false.)
   write (dat%s%out%name, "(a32)") "s_"//adjustl(name)//".dat"
   dat%s%out%name=adjustl(dat%s%out%name)
   call init_writeouts(dat%s,ugrid_out)
   dat%counter=>counter
   
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
   write (dat%out%name, "(a32)") adjustl(name)//".dat"
   dat%out%name=adjustl(dat%out%name)
   call init_writeouts(dat,vgrid_out)
   allocate(dat%t)
   call create_var(dat%t,vgrid,.false.)
   write (dat%t%out%name, "(a32)") "s_"//adjustl(name)//".dat"
   dat%t%out%name=adjustl(dat%t%out%name)
   call init_writeouts(dat%t,vgrid_out)
   allocate(dat%s)
   call create_var(dat%s,vgrid,.false.)
   write (dat%s%out%name, "(a32)") "s_"//adjustl(name)//".dat"
   dat%s%out%name=adjustl(dat%s%out%name)
   call init_writeouts(dat%s,vgrid_out)
   dat%counter=>counter
   
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
   write (dat%out%name, "(a32)") adjustl(name)//".dat"
   dat%out%name=adjustl(dat%out%name)
   call init_writeouts(dat,zgrid_out)
   allocate(dat%t)
   call create_var(dat%t,zgrid,.false.)
   write (dat%t%out%name, "(a32)") "s_"//adjustl(name)//".dat"
   dat%t%out%name=adjustl(dat%t%out%name)
   call init_writeouts(dat%t,zgrid_out)
   allocate(dat%s)
   call create_var(dat%s,zgrid,.false.)
   write (dat%s%out%name, "(a32)") "s_"//adjustl(name)//".dat"
   dat%s%out%name=adjustl(dat%s%out%name)
   call init_writeouts(dat%s,zgrid_out)
   dat%counter=>counter
   
   
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
  
  
  
  
  
  
  
  
  subroutine create_written_hvar_layer(dat,name,synced,i)
   use global
   use writeout_grid
   use writeout
   
   implicit none
   
   type(hvar), allocatable, intent(inout) :: dat(:)
   logical, intent(in) :: synced
   integer, intent(in) :: i
   integer :: k
   character(*), intent(in) :: name
   
   
   
   allocate(dat(i:nz))
   do k=lbound(dat,1),ubound(dat,1)
    call create_hvar(dat(k),synced)
    write (dat(k)%out%name, "(a27,i1,a4)") adjustl(name)//"_",k,".dat"
    dat(k)%out%name=adjustl(dat(k)%out%name)
    call init_writeouts(dat(k),hgrid_out)
    allocate(dat(k)%t)
    call create_var(dat(k)%t,hgrid,.false.)
    write (dat(k)%t%out%name, "(a27,i1,a4)") "s_"//adjustl(name)//"_",k,".dat"
    dat(k)%t%out%name=adjustl(dat(k)%t%out%name)
    call init_writeouts(dat(k)%t,hgrid_out)
    allocate(dat(k)%s)
    call create_var(dat(k)%s,hgrid,.false.)
    write (dat(k)%s%out%name, "(a27,i1,a4)") "s_"//adjustl(name)//"_",k,".dat"
    dat(k)%s%out%name=adjustl(dat(k)%s%out%name)
    call init_writeouts(dat(k)%s,hgrid_out)
    dat(k)%counter=>counter
   end do
   
  end subroutine
  
  
  subroutine create_written_uvar_layer(dat,name,synced,i)
   use global
   use writeout_grid
   use writeout
   
   implicit none
   
   type(uvar), allocatable, intent(inout) :: dat(:)
   logical, intent(in) :: synced
   integer, intent(in) :: i
   integer :: k
   character(*), intent(in) :: name
   
   
   
   allocate(dat(i:nz))
   do k=lbound(dat,1),ubound(dat,1)
    call create_uvar(dat(k),synced)
    write (dat(k)%out%name, "(a27,i1,a4)") adjustl(name)//"_",k,".dat"
    dat(k)%out%name=adjustl(dat(k)%out%name)
    call init_writeouts(dat(k),ugrid_out)
    allocate(dat(k)%t)
    call create_var(dat(k)%t,ugrid,.false.)
    write (dat(k)%t%out%name, "(a27,i1,a4)") "s_"//adjustl(name)//"_",k,".dat"
    dat(k)%t%out%name=adjustl(dat(k)%t%out%name)
    call init_writeouts(dat(k)%t,ugrid_out)
    dat(k)%counter=>counter
    allocate(dat(k)%s)
    call create_var(dat(k)%s,ugrid,.false.)
    write (dat(k)%s%out%name, "(a27,i1,a4)") "s_"//adjustl(name)//"_",k,".dat"
    dat(k)%s%out%name=adjustl(dat(k)%s%out%name)
    call init_writeouts(dat(k)%s,ugrid_out)
    dat(k)%counter=>counter
   end do
   
  end subroutine
  
  
  subroutine create_written_vvar_layer(dat,name,synced,i)
   use global
   use writeout_grid
   use writeout
   
   implicit none
   
   type(vvar), allocatable, intent(inout) :: dat(:)
   logical, intent(in) :: synced
   integer, intent(in) :: i
   integer :: k
   character(*), intent(in) :: name
   
   
   allocate(dat(i:nz))
   do k=lbound(dat,1),ubound(dat,1)
    call create_vvar(dat(k),synced)
    write (dat(k)%out%name, "(a27,i1,a4)") adjustl(name)//"_",k,".dat"
    dat(k)%out%name=adjustl(dat(k)%out%name)
    call init_writeouts(dat(k),vgrid_out)
    allocate(dat(k)%t)
    call create_var(dat(k)%t,vgrid,.false.)
    write (dat(k)%t%out%name, "(a27,i1,a4)") "s_"//adjustl(name)//"_",k,".dat"
    dat(k)%t%out%name=adjustl(dat(k)%t%out%name)
    call init_writeouts(dat(k)%t,vgrid_out)
    allocate(dat(k)%s)
    call create_var(dat(k)%s,vgrid,.false.)
    write (dat(k)%s%out%name, "(a27,i1,a4)") "s_"//adjustl(name)//"_",k,".dat"
    dat(k)%s%out%name=adjustl(dat(k)%s%out%name)
    call init_writeouts(dat(k)%s,vgrid_out)
    dat(k)%counter=>counter
   end do
   
  end subroutine
  
  
  subroutine create_written_zvar_layer(dat,name,synced,i)
   use global
   use writeout_grid
   use writeout
   
   implicit none
   
   type(zvar), allocatable, intent(inout) :: dat(:)
   logical, intent(in) :: synced
   integer, intent(in) :: i
   integer :: k
   character(*), intent(in) :: name
   
   
   
   allocate(dat(i:nz))
   do k=lbound(dat,1),ubound(dat,1)
    call create_zvar(dat(k),synced)
    write (dat(k)%out%name, "(a27,i1,a4)") adjustl(name)//"_",k,".dat"
    dat(k)%out%name=adjustl(dat(k)%out%name)
    call init_writeouts(dat(k),zgrid_out)
    allocate(dat(k)%t)
    call create_var(dat(k)%t,zgrid,.false.)
    write (dat(k)%t%out%name, "(a27,i1,a4)") "s_"//adjustl(name)//"_",k,".dat"
    dat(k)%t%out%name=adjustl(dat(k)%t%out%name)
    call init_writeouts(dat(k)%t,zgrid_out)
    allocate(dat(k)%s)
    call create_var(dat(k)%s,zgrid,.false.)
    write (dat(k)%s%out%name, "(a27,i1,a4)") "s_"//adjustl(name)//"_",k,".dat"
    dat(k)%s%out%name=adjustl(dat(k)%s%out%name)
    call init_writeouts(dat(k)%s,zgrid_out)
    dat(k)%counter=>counter
   end do
   
  end subroutine
  
  
  
  
  
  
  
  
  
  
  
  
  subroutine create_hvar_layer(dat,synced,i)
   use global
   
   implicit none
   
   type(hvar), allocatable, intent(inout) :: dat(:)
   integer, intent(in) :: i
   logical, intent(in) :: synced
   integer :: k
   
   allocate(dat(i:nz))
   do k=lbound(dat,1),ubound(dat,1)
    call create_hvar(dat(k),synced)
   end do
   
  end subroutine
  
  
  subroutine create_uvar_layer(dat,synced,i)
   use global
   
   implicit none
   
   type(uvar), allocatable, intent(inout) :: dat(:)
   integer, intent(in) :: i
   logical, intent(in) :: synced
   integer :: k
   
   allocate(dat(i:nz))
   do k=lbound(dat,1),ubound(dat,1)
    call create_uvar(dat(k),synced)
   end do
   
  end subroutine
  
  
  subroutine create_vvar_layer(dat,synced,i)
   use global
   
   implicit none
   
   type(vvar), allocatable, intent(inout) :: dat(:)
   integer, intent(in) :: i
   logical, intent(in) :: synced
   integer :: k
   
   allocate(dat(i:nz))
   do k=lbound(dat,1),ubound(dat,1)
    call create_vvar(dat(k),synced)
   end do
   
  end subroutine
  
  
  subroutine create_zvar_layer(dat,synced,i)
   use global
   
   implicit none
   
   type(zvar), allocatable, intent(inout) :: dat(:)
   integer, intent(in) :: i
   logical, intent(in) :: synced
   integer :: k
   
   allocate(dat(i:nz))
   do k=lbound(dat,1),ubound(dat,1)
    call create_zvar(dat(k),synced)
   end do
   
  end subroutine
  
  subroutine remove_field(dat)
   use global
   
   class(var), intent(inout) :: dat
   
   nullify(dat%p)
   nullify(dat%bz)
   if (associated(dat%z)) then
    deallocate(dat%z)
   end if
   
  end subroutine
  
  
end module
