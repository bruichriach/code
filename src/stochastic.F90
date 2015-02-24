module linklists
 use global

 implicit none
 
 type llist
  real(kind=db) :: rand(6)
  real(kind=db) :: timenow = 0.0d0
  real(kind=db) :: timeend
  logical :: first, last
  type(llist), pointer :: next
  type(llist), pointer :: tmp
  type(uvar) :: u
  type(vvar) :: v
 end type
 
 type(llist), pointer :: stochwind
 
end module
  
module llist_ops
 use linklists
 
 implicit none
 
 contains
 
 subroutine create_link(linkedlist)
  use linklists
 
  implicit none
  
  type(llist), pointer :: linkedlist
  
  if (associated(linkedlist)) then
   do while (.not.(linkedlist%last))
    linkedlist => linkedlist%next
   end do
   linkedlist%tmp => linkedlist%next
   nullify(linkedlist%next)
   allocate(linkedlist%next)
   linkedlist%last=.false.
   linkedlist%next%next => linkedlist%tmp
   nullify(linkedlist%tmp)
   linkedlist => linkedlist%next
  else
   allocate(linkedlist)
   linkedlist%next => linkedlist
   linkedlist%first=.true.
  end if
  call random_number(linkedlist%rand)
  call random_number(linkedlist%timeend)
  linkedlist%last=.true.
  
 end subroutine
  
 subroutine remove_link(linkedlist)
  use linklists
  use grid
  use allocation
  
  implicit none
  
  type(llist), pointer :: linkedlist
  
  if (associated(linkedlist)) then
   if (linkedlist%next%timenow >= linkedlist%next%timeend) then
    call remove_field(linkedlist%next%u)
    call remove_field(linkedlist%next%v)
    if (.not.(associated(linkedlist,linkedlist%next))) then
     if (linkedlist%next%first) linkedlist%next%next%first=.true.
     if (linkedlist%next%last) linkedlist%last=.true.
     linkedlist%tmp => linkedlist%next%next
     nullify(linkedlist%next%next)
     deallocate(linkedlist%next)
     linkedlist%next => linkedlist%tmp
     nullify(linkedlist%tmp)
    else
     deallocate(linkedlist)
    end if
   else
    linkedlist => linkedlist%next
   end if
  end if
  
 end subroutine
 
 subroutine step_link(linkedlist,timestep)
  use linklists
  
  implicit none
  
  type(llist), pointer :: linkedlist
  real(kind=db), intent(in) :: timestep
  
  
 ! call remove_link(linkedlist)
  linkedlist => linkedlist%next
  linkedlist%timenow = linkedlist%timenow + timestep
  
 end subroutine
  
 subroutine set_wind(n)
  use global
  use linklists
  use params
  use variables
  use grid_operate
  use grid
  use allocation
 
  implicit none
  
  integer, intent(in) :: n
  real(kind=db) :: r_u(utau%p%lx+1:utau%p%lx+utau%p%nx,utau%p%ly+1:utau%p%ly+utau%p%ny)
  real(kind=db) :: r_v(vtau%p%lx+1:vtau%p%lx+vtau%p%nx,vtau%p%ly+1:vtau%p%ly+vtau%p%ny)

  if (associated(stochwind)) then
   do while (.not.(stochwind%next%last))
    call remove_link(stochwind)
   end do
   call remove_link(stochwind)
  end if
  if (mod(n,ceiling(1.0d0*8.64d4/(dt/f0))) == 0.0d0) then
   call create_link(stochwind)
   call create_field(stochwind%u,.false.)
   call create_field(stochwind%v,.false.)
   stochwind%timeend=6.048d5*(4.0d0+12.0d0*stochwind%timeend)*f0
   stochwind%rand(1)=sign(tau*0.5d0,stochwind%rand(1))+  &
              tau*0.5d0*(stochwind%rand(1))
   stochwind%rand(3)=x0*(stochwind%rand(3))
   stochwind%rand(4)=y0*(stochwind%rand(4))
   stochwind%rand(5)=2.0d0*(stochwind%rand(5)-0.5d0)
   stochwind%rand(5)=0.25d0*x0*(stochwind%rand(5))
   stochwind%rand(6)=2.0d0*(stochwind%rand(6)-0.5d0)
   stochwind%rand(6)=0.25d0*y0*(stochwind%rand(6))
   stochwind%rand(2)=(1.5d0+(x0/4.0d0-1.5d0)*(stochwind%rand(2)))
   r_u=sqrt(x_dist(utau,lim(stochwind%rand(3)+stochwind%rand(5),x0))**2 + &
            y_dist(utau,lim(stochwind%rand(4)+stochwind%rand(6),y0))**2)
   stochwind%u=-stochwind%rand(1)*merge(0.0d0,  &
          y_dist(utau,lim(stochwind%rand(4)+stochwind%rand(6),y0))/(r_u*stochwind%rand(2)),r_u == 0.0d0)    &
          *(1.0d0/(2.0d0*pi*r_u))*(1.0d0-exp(-r_u**2))*exp(-r_u**2/2.0d0)*    &
          (1.0d0-exp(-tan((pi/x0)*(x0/2.0d0-  &
          abs(x_dist(utau,lim(stochwind%rand(3)+stochwind%rand(5),x0)))))**2/0.25d0))*   &
          (1.0d0-exp(-tan((pi/y0)*(y0/2.0d0-  &
          abs(y_dist(utau,lim(stochwind%rand(4)+stochwind%rand(6),y0)))))**2/0.25d0))   
          
   r_u=sqrt(x_dist(utau,lim(stochwind%rand(3)-stochwind%rand(5),x0))**2 + &
            y_dist(utau,lim(stochwind%rand(4)-stochwind%rand(6),y0))**2)
   stochwind%u= stochwind%u%bz+stochwind%rand(1)*merge(0.0d0,  &
          y_dist(utau,lim(stochwind%rand(4)-stochwind%rand(6),y0))/(r_u*stochwind%rand(2)),r_u == 0.0d0)    &
          *(1.0d0/(2.0d0*pi*r_u))*(1.0d0-exp(-r_u**2))*exp(-r_u**2/2.0d0)*    &
          (1.0d0-exp(-tan((pi/x0)*(x0/2.0d0-  &
          abs(x_dist(utau,lim(stochwind%rand(3)-stochwind%rand(5),x0)))))**2/0.25d0))*   &
          (1.0d0-exp(-tan((pi/y0)*(y0/2.0d0-  &
          abs(y_dist(utau,lim(stochwind%rand(4)-stochwind%rand(6),y0)))))**2/0.25d0))  
                      
                      
   r_v=sqrt(x_dist(vtau,lim(stochwind%rand(3)+stochwind%rand(5),x0))**2 + &
            y_dist(vtau,lim(stochwind%rand(4)+stochwind%rand(6),y0))**2)
   stochwind%v=stochwind%rand(1)*merge(0.0d0,  &
          x_dist(vtau,lim(stochwind%rand(4)+stochwind%rand(6),y0))/(r_v*stochwind%rand(2)),r_v == 0.0d0)    &
          *(1.0d0/(2.0d0*pi*r_v))*(1.0d0-exp(-r_v**2))*exp(-r_v**2/2.0d0)*    &
          (1.0d0-exp(-tan((pi/x0)*(x0/2.0d0-  &
          abs(x_dist(vtau,lim(stochwind%rand(3)+stochwind%rand(5),x0)))))**2/0.25d0))*   &
          (1.0d0-exp(-tan((pi/y0)*(y0/2.0d0-  &
          abs(y_dist(vtau,lim(stochwind%rand(4)+stochwind%rand(6),y0)))))**2/0.25d0))   
          
   r_v=sqrt(x_dist(vtau,lim(stochwind%rand(3)-stochwind%rand(5),x0))**2 + &
            y_dist(vtau,lim(stochwind%rand(4)-stochwind%rand(6),y0))**2)
   stochwind%v= stochwind%v%bz-stochwind%rand(1)*merge(0.0d0,  &
          x_dist(vtau,lim(stochwind%rand(4)-stochwind%rand(6),y0))/(r_v*stochwind%rand(2)),r_v == 0.0d0)    &
          *(1.0d0/(2.0d0*pi*r_v))*(1.0d0-exp(-r_v**2))*exp(-r_v**2/2.0d0)*    &
          (1.0d0-exp(-tan((pi/x0)*(x0/2.0d0-  &
          abs(x_dist(vtau,lim(stochwind%rand(3)-stochwind%rand(5),x0)))))**2/0.25d0))*   &
          (1.0d0-exp(-tan((pi/y0)*(y0/2.0d0-  &
          abs(y_dist(vtau,lim(stochwind%rand(4)-stochwind%rand(6),y0)))))**2/0.25d0))  
                      
                      
  end if  

  utau%bz=0.0d0
  vtau%bz=0.0d0
  call step_link(stochwind,dt)
  utau%bz=utau%bz+   &
        sin(min(2.0d0*pi,2.0d0*pi*stochwind%timenow/stochwind%timeend))*   &
        stochwind%u%bz
  vtau%bz=vtau%bz+   &
        sin(min(2.0d0*pi,2.0d0*pi*stochwind%timenow/stochwind%timeend))*   &
        stochwind%v%bz
!  print *, ceiling(stochwind%rand), ceiling(stochwind%timenow), &
!        ceiling(stochwind%timeend)
  do while (.not.(stochwind%last))
   call step_link(stochwind,dt)
   utau%bz=utau%bz+   &
        sin(min(2.0d0*pi,2.0d0*pi*stochwind%timenow/stochwind%timeend))*   &
        stochwind%u%bz
   vtau%bz=vtau%bz+   &
        sin(min(2.0d0*pi,2.0d0*pi*stochwind%timenow/stochwind%timeend))*   &
        stochwind%v%bz
!   print *, ceiling(stochwind%rand), ceiling(stochwind%timenow), &
!        ceiling(stochwind%timeend)
  end do
  
 end subroutine
  
  
 subroutine output_llist(linkedlist,varname)
  use linklists
  use parallel
  
  implicit none
  
  type(llist), pointer :: linkedlist
  character(*), intent(in) :: varname 
  character(32) :: filename, foldername, fullname
  
  if (proc_name == ens_master) then
   write (foldername, "(a5,i4.4)") 'out/', ens_name
   write (filename, "(a32)") trim(varname)//'.dat'
   write (fullname, "(a32)") trim(adjustl(foldername))//'/'//adjustl(filename)
   if (associated(linkedlist)) then
    open(unit=10,file=adjustl(fullname),status='unknown',form='unformatted')
    
    linkedlist => linkedlist%next
    write(10) linkedlist%rand, linkedlist%timenow, linkedlist%timeend
    do while (.not.(stochwind%last))
     linkedlist => linkedlist%next
     write(10) linkedlist%rand, linkedlist%timenow, linkedlist%timeend
    end do
    
    close(10)
   end if
  end if
  
 end subroutine
  
  
 subroutine writeout_llist(num,linkedlist,varname)
  use linklists
  use parallel
  
  implicit none
  
  type(llist), pointer :: linkedlist
  character(*), intent(in) :: varname 
  character(32) :: filename, foldername, fullname
  integer, intent(in) :: num
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
 
   write (foldername, "(a5,i4.4,a1,i4.4)") 'data/', ens_name,'/', num
   write (filename, "(a32)") trim(varname)//'.dat'
   write (fullname, "(a32)") trim(adjustl(foldername))//'/'//adjustl(filename)
   if (associated(linkedlist)) then
    open(unit=10,file=adjustl(fullname),status='unknown',form='unformatted')
    
    linkedlist => linkedlist%next
    write(10) linkedlist%rand, linkedlist%timenow, linkedlist%timeend
    do while (.not.(stochwind%last))
     linkedlist => linkedlist%next
     write(10) linkedlist%rand, linkedlist%timenow, linkedlist%timeend
    end do
    
    close(10)
   end if
  end if
  
 end subroutine
   
  
 subroutine input_llist(linkedlist,varname)
  use linklists
  use variables
  use params
  use grid_operate
  use grid
  use parallel
  use allocation

  
  implicit none
  
  type(llist), pointer :: linkedlist
  type(llist), pointer :: tmplist
  character(*), intent(in) :: varname 
  logical :: file_exist
  integer :: endoffile, record=0
  real(kind=db) :: r_u(utau%p%lx+1:utau%p%lx+utau%p%nx,utau%p%ly+1:utau%p%ly+utau%p%ny)
  real(kind=db) :: r_v(vtau%p%lx+1:vtau%p%lx+vtau%p%nx,vtau%p%ly+1:vtau%p%ly+vtau%p%ny)
  character(32) :: filename, foldername, fullname
  
   write (foldername, "(a5,i4.4)") 'old/', ens_name
   write (filename, "(a32)") trim(varname)//'.dat'
   write (fullname, "(a32)") trim(adjustl(foldername))//'/'//adjustl(filename)
  inquire( file=adjustl(fullname), exist=file_exist)
  if (proc_name == proc_master) print *, "reading in wind data"
    
  if (file_exist) then
 
   open(unit=10,file=adjustl(fullname),status='unknown',form='unformatted')
   allocate(tmplist)
   read(10,iostat=endoffile) tmplist%rand, tmplist%timenow, tmplist%timeend
  if (proc_name == proc_master) print *, endoffile, tmplist%rand, tmplist%timenow, tmplist%timeend
   do while (endoffile == 0)
    if (.not.(associated(linkedlist))) then
     linkedlist => tmplist
     linkedlist%next => linkedlist
     linkedlist%first = .true.
     linkedlist%last = .true.
     nullify(tmplist)
    else
     do while (.not.(linkedlist%last))
      linkedlist => linkedlist%next
     end do      
     linkedlist%tmp => linkedlist%next
     linkedlist%next => tmplist
     linkedlist%next%last = .true.
     linkedlist%last = .false.
     linkedlist%next%next => linkedlist%tmp
     nullify(linkedlist%tmp)
     nullify(tmplist)
     linkedlist => linkedlist%next
   
    end if
    if (proc_name == proc_master) then
     record=record + 1
     print *, "reading in wind data record #", record
    end if
    call create_field(linkedlist%u,.false.)
    call create_field(linkedlist%v,.false.)
   
   r_u=sqrt(x_dist(utau,lim(linkedlist%rand(3)+linkedlist%rand(5),x0))**2 + &
            y_dist(utau,lim(linkedlist%rand(4)+linkedlist%rand(6),y0))**2)
   linkedlist%u=-linkedlist%rand(1)*merge(0.0d0,  &
          y_dist(utau,lim(linkedlist%rand(4)+linkedlist%rand(6),y0))/(r_u*linkedlist%rand(2)),r_u == 0.0d0)    &
          *(1.0d0/(2.0d0*pi*r_u))*(1.0d0-exp(-r_u**2))*exp(-r_u**2/2.0d0)*    &
          (1.0d0-exp(-tan((pi/x0)*(x0/2.0d0-  &
          abs(x_dist(utau,lim(linkedlist%rand(3)+linkedlist%rand(5),x0)))))**2/0.25d0))*   &
          (1.0d0-exp(-tan((pi/y0)*(y0/2.0d0-  &
          abs(y_dist(utau,lim(linkedlist%rand(4)+linkedlist%rand(6),y0)))))**2/0.25d0))   
          
   r_u=sqrt(x_dist(utau,lim(linkedlist%rand(3)-linkedlist%rand(5),x0))**2 + &
            y_dist(utau,lim(linkedlist%rand(4)-linkedlist%rand(6),y0))**2)
   linkedlist%u= linkedlist%u%bz+linkedlist%rand(1)*merge(0.0d0,  &
          y_dist(utau,lim(linkedlist%rand(4)-linkedlist%rand(6),y0))/(r_u*linkedlist%rand(2)),r_u == 0.0d0)    &
          *(1.0d0/(2.0d0*pi*r_u))*(1.0d0-exp(-r_u**2))*exp(-r_u**2/2.0d0)*    &
          (1.0d0-exp(-tan((pi/x0)*(x0/2.0d0-  &
          abs(x_dist(utau,lim(linkedlist%rand(3)-linkedlist%rand(5),x0)))))**2/0.25d0))*   &
          (1.0d0-exp(-tan((pi/y0)*(y0/2.0d0-  &
          abs(y_dist(utau,lim(linkedlist%rand(4)-linkedlist%rand(6),y0)))))**2/0.25d0))  
                      
                      
   r_v=sqrt(x_dist(vtau,lim(linkedlist%rand(3)+linkedlist%rand(5),x0))**2 + &
            y_dist(vtau,lim(linkedlist%rand(4)+linkedlist%rand(6),y0))**2)
   linkedlist%v=linkedlist%rand(1)*merge(0.0d0,  &
          x_dist(vtau,lim(linkedlist%rand(4)+linkedlist%rand(6),y0))/(r_v*linkedlist%rand(2)),r_v == 0.0d0)    &
          *(1.0d0/(2.0d0*pi*r_v))*(1.0d0-exp(-r_v**2))*exp(-r_v**2/2.0d0)*    &
          (1.0d0-exp(-tan((pi/x0)*(x0/2.0d0-  &
          abs(x_dist(vtau,lim(linkedlist%rand(3)+linkedlist%rand(5),x0)))))**2/0.25d0))*   &
          (1.0d0-exp(-tan((pi/y0)*(y0/2.0d0-  &
          abs(y_dist(vtau,lim(linkedlist%rand(4)+linkedlist%rand(6),y0)))))**2/0.25d0))   
          
   r_v=sqrt(x_dist(vtau,lim(linkedlist%rand(3)-linkedlist%rand(5),x0))**2 + &
            y_dist(vtau,lim(linkedlist%rand(4)-linkedlist%rand(6),y0))**2)
   linkedlist%v= linkedlist%v%bz-linkedlist%rand(1)*merge(0.0d0,  &
          x_dist(vtau,lim(linkedlist%rand(4)-linkedlist%rand(6),y0))/(r_v*linkedlist%rand(2)),r_v == 0.0d0)    &
          *(1.0d0/(2.0d0*pi*r_v))*(1.0d0-exp(-r_v**2))*exp(-r_v**2/2.0d0)*    &
          (1.0d0-exp(-tan((pi/x0)*(x0/2.0d0-  &
          abs(x_dist(vtau,lim(linkedlist%rand(3)-linkedlist%rand(5),x0)))))**2/0.25d0))*   &
          (1.0d0-exp(-tan((pi/y0)*(y0/2.0d0-  &
          abs(y_dist(vtau,lim(linkedlist%rand(4)-linkedlist%rand(6),y0)))))**2/0.25d0))  

   utau%bz=utau%bz+   &
        sin(min(2.0d0*pi,2.0d0*pi*linkedlist%timenow/linkedlist%timeend))*   &
        linkedlist%u%bz
   vtau%bz=vtau%bz+   &
        sin(min(2.0d0*pi,2.0d0*pi*linkedlist%timenow/linkedlist%timeend))*   &
        linkedlist%v%bz



    allocate(tmplist)
    read(10,iostat=endoffile) tmplist%rand, tmplist%timenow, tmplist%timeend
  if (proc_name == proc_master) print *, endoffile
   end do 
   
   deallocate(tmplist)
    
  
   close(10)
   
  end if
  
 end subroutine
 
end module
  
