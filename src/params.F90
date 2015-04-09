module params

#include "include.h"

 use sys
 
 implicit none

 integer :: mx, my
 integer :: nz
 real(kind=db) :: dx, dy
 real(kind=db) :: indx, indy
#ifdef DOUBLE_PRECISION
 real(kind=db) :: x0, y0
#else
 real(kind=db) :: x0, y0
#endif
 
 
 integer :: lx, ly, nx, ny, north, south, east, west
 
 logical :: noslip
 
 
 
 
 real (kind=db), parameter :: pi=4.0d0*atan(1.0d0)
 real (kind=db) :: omega
 real (kind=db) :: g
 
  ! Degree/radian conversion parameters
 real (kind=db), parameter :: deg2rad=pi/180.0d0, rad2deg=180.0d0/pi
 
 real (kind=db), allocatable :: gp(:) 
 real (kind=db), allocatable :: ngp(:)
#ifdef DO_TIME_AVERAGE
 real (kind=db), allocatable :: ng(:)
#endif

 real (kind=db) :: f0
 real (kind=db) :: h0
 real (kind=db) :: ld


 real (kind=db) :: maxh
 real (kind=db) :: dt






 REAL (KIND=db) :: write_time
 REAL (KIND=db) :: total_time
#ifdef DO_TIME_AVERAGE
 REAL (KIND=db) :: average_time
 
 INTEGER :: nsteps
#endif
 
 INTEGER :: nstop
 INTEGER :: wstep
#ifdef DO_TIME_AVERAGE
 INTEGER :: nstep
#endif


     
 REAL (kind=db) :: bf
 real (kind=db) :: res
 REAL (kind=db) :: tau    
 real (kind=db) :: umax
 integer :: bfricpoly
 real (kind=db) :: cvar
 
 

#ifdef DO_TIME_AVERAGE
 integer :: timeav_count
#endif
 
 contains
 
 subroutine init_params()
  use parallel
 
  implicit none
  
  integer :: ii
  character(32) :: desc
  character(1) :: colon
  logical :: file_exist
  integer :: endoffile

  mx=100
  my=150
  nz=3
  dx=0.1d0
  dy=0.1d0
  indx=1.0d0/dx
  indy=1.0d0/dy
#ifdef DOUBLE_PRECISION
  x0=dx*dble(mx)
  y0=dy*dble(my)
#else
  x0=dx*real(mx)
  y0=dy*real(my)
#endif

  noslip = .true.
 
  omega=2.0d0*pi/8.64d4
  g=9.81d0
 

#ifndef ALLOW_STATIC_LAYER
#ifdef ALLOW_RIGID_LID
  allocate(gp(nz-1),ngp(nz))
  gp = (/ 3.0d-3 , 7.0d-3 /)
  ngp = (/ gp(1:nz-1) , g /)/sum(gp(1:nz-1))
#ifdef DO_TIME_AVERAGE
  allocate(ng(nz))
  ng =  (/ (  sum(ngp(ii:nz) )  , ii = 1,nz) /)
#endif
#else
  allocate(gp(nz),ngp(nz))
  gp = (/ 3.0d-3 , 7.0d-3, 9.81d0 /)
  ngp = gp(1:nz)/sum(gp(1:nz-1))
#ifdef DO_TIME_AVERAGE
  allocate(ng(nz))
  ng =  (/ (  sum(ngp(ii:nz) )  , ii = 1,nz) /)
#endif
#endif
#else
  allocate(gp(nz),ngp(nz))
  gp = (/ 3.0d-3 , 7.0d-3 /)
  ngp = gp(1:nz)/sum(gp(1:nz))
#ifdef DO_TIME_AVERAGE
  allocate(ng(nz))
  ng = (/ (  sum(ngp(ii:nz) )  , ii = 1,nz) /)
#endif
#endif

  f0=2.0d0*omega*sin(deg2rad*70.0d0)
  h0=1.0d3
#ifndef ALLOW_STATIC_LAYER
  ld=sqrt(sum(gp(1:nz-1))*h0)/f0
#else
  ld=sqrt(sum(gp(1:nz))*h0)/f0
#endif




#ifdef ALLOW_RIGID_LID
  maxh=1.25d0!h0
  dt = (0.25d0*dx)/dsqrt(sum(ngp(1:nz-1))*maxh)
#else
  maxh=1.25d0*1.01d0
  dt = (0.25d0*dx)/dsqrt(sum(ngp)*maxh)
#endif






  write_time=3.0d1*(4.0d0*pi*sin(deg2rad*70.0d0))
  total_time=2.0d4*(4.0d0*pi*sin(deg2rad*70.0d0))
#ifdef DO_TIME_AVERAGE
  average_time=1.0d0
 
  nsteps=floor(write_time/average_time)
#endif
 
  nstop=ceiling(total_time/dt)
  wstep=floor(write_time/dt)
#ifdef DO_TIME_AVERAGE
  nstep=0
#endif


     
  bf=2.0d-2
  res=0.0d0/(f0**2*ld*h0) 
  tau=2.0d-1/(f0**2*ld*h0*1026.0d0)      
  umax=2.0d-2
  bfricpoly=2
  cvar=(1.0d0/pi**2)  
  
  

#ifdef DO_TIME_AVERAGE
 timeav_count=10561868
#endif
  
  
  inquire( file='params.txt', exist=file_exist)
    
  if (file_exist) then
   open(unit=10,action='read',file='params.txt',   &
        status='unknown',form='formatted')
   read(10,"(a32,a1)",iostat=endoffile,advance='NO') desc, colon 
   do while (endoffile /= -1)
    select case (desc)
    
     case ('x-direction grid spaces')
      format="(i4)"
      read(10,format,iostat=endoffile,advance='YES') mx
      if (proc_name == proc_master) print *, 'Reading record: ', desc//':', mx
      
     case ('y-direction grid spaces')
      format="(i4)"
      read(10,format,iostat=endoffile,advance='YES') my
      if (proc_name == proc_master) print *, 'Reading record: ', desc//':', my
      
     case ('z-direction # of layers')
      format="(i4)"
      read(10,format,iostat=endoffile,advance='YES') nz
      if (proc_name == proc_master) print *, 'Reading record: ', desc//':', nz
      
     case ('x-direction grid width')
      format="(e23.16)"
      read(10,format,iostat=endoffile,advance='YES') dx
      if (proc_name == proc_master) print *, 'Reading record: ', desc//':', dx
      
     case ('y-direction grid width')
      format="(e23.16)"
      read(10,format,iostat=endoffile,advance='YES') dy
      if (proc_name == proc_master) print *, 'Reading record: ', desc//':', dy
      
     case ('no slip')
      format="(l3)"
      read(10,format,iostat=endoffile,advance='YES') noslip
      if (proc_name == proc_master) print *, 'Reading record: ', desc//':', noslip
      
     case ('planetary rotational frequency')
      format="(e23.16)"
      read(10,format,iostat=endoffile,advance='YES') omega
      if (proc_name == proc_master) print *, 'Reading record: ', desc//':', omega
      
     case ('gravitational acceleration')
      format="(e23.16)"
      read(10,format,iostat=endoffile,advance='YES') g
      if (proc_name == proc_master) print *, 'Reading record: ', desc//':', g
      
     case ('reduced gravitational accel')
      write(format,"(a1,i1,a7)") '(',ubound(gp,1)-lbound(gp,1)+1,'e23.16)'
      read(10,format,iostat=endoffile,advance='YES')  (gp(ii) , ii=lbound(gp,1),ubound(gp,1))
      if (proc_name == proc_master) print *, 'Reading record: ', desc//':', gp
      
     case ('dimensional coriolis')
      format="(e23.16)"
      read(10,format,iostat=endoffile,advance='YES') f0
      if (proc_name == proc_master) print *, 'Reading record: ', desc//':', f0
      
     case ('dimensional depth')
      format="(e23.16)"
      read(10,format,iostat=endoffile,advance='YES') h0
      if (proc_name == proc_master) print *, 'Reading record: ', desc//':', h0
      
     case ('dimensional deformation radius')
      format="(e23.16)"
      read(10,format,iostat=endoffile,advance='YES') ld
      if (proc_name == proc_master) print *, 'Reading record: ', desc//':', ld
      
     case ('timestep')
      format="(e23.16)"
      read(10,format,iostat=endoffile,advance='YES') dt
      if (proc_name == proc_master) print *, 'Reading record: ', desc//':', dt
      
     case ('runtime')
      format="(e23.16)"
      read(10,format,iostat=endoffile,advance='YES') total_time
      if (proc_name == proc_master) print *, 'Reading record: ', desc//':', total_time
      
     case ('writetime')
      format="(e23.16)"
      read(10,format,iostat=endoffile,advance='YES') write_time
      if (proc_name == proc_master) print *, 'Reading record: ', desc//':', write_time
      
     case ('bottom friction')
      format="(e23.16)"
      read(10,format,iostat=endoffile,advance='YES') bf
      if (proc_name == proc_master) print *, 'Reading record: ', desc//':', bf

     case ('bottom friction velocity')
      format="(e23.16)"
      read(10,format,iostat=endoffile,advance='YES') umax
      if (proc_name == proc_master) print *, 'Reading record: ', desc//':', umax

     case ('bottom friction polynomial')
      format="(i4)"
      read(10,format,iostat=endoffile,advance='YES') bfricpoly
      if (proc_name == proc_master) print *, 'Reading record: ', desc//':', bfricpoly
 
     case ('wind stress')
      format="(e23.16)"
      read(10,format,iostat=endoffile,advance='YES') tau
      if (proc_name == proc_master) print *, 'Reading record: ', desc//':', tau
      
     case ('viscosity')
      format="(e23.16)"
      read(10,format,iostat=endoffile,advance='YES') cvar   
      if (proc_name == proc_master) print *, 'Reading record: ', desc//':', cvar  
      
#ifdef DO_TIME_AVERAGE
     case ('initial time sum itereation')
      format="(i16)"
      read(10,format,iostat=endoffile,advance='YES') timeav_count
      if (proc_name == proc_master) print *, 'Reading record: ', desc//':', timeav_count
#endif
      
     case default
      if (proc_name == proc_master) print *, 'No param for record: ', desc
      read(10,"(a1)",iostat=endoffile,advance='YES')
      endoffile=0
      

      
    end select
    if (endoffile == -2) then
     if (proc_name == proc_master) print *, 'Read error for record: ', desc
    end if
    read(10,"(a32,a1)",iostat=endoffile,advance='NO') desc, colon
   end do
   close(10)
   
   
  indx=1.0d0/dx
  indy=1.0d0/dy
#ifdef DOUBLE_PRECISION
  x0=dx*dble(mx)
  y0=dy*dble(my)
#else
  x0=dx*real(mx)
  y0=dy*real(my)
#endif
 

#ifndef ALLOW_STATIC_LAYER
#ifdef ALLOW_RIGID_LID
  ngp = (/ gp(1:nz-1) , g /)/sum(gp(1:nz-1))
#ifdef DO_TIME_AVERAGE
  ng =  (/ (  sum(ngp(ii:nz) )  , ii = 1,nz) /)
#endif
#else
  ngp = gp(1:nz)/sum(gp(1:nz-1))
#ifdef DO_TIME_AVERAGE
  ng =  (/ (  sum(ngp(ii:nz) )  , ii = 1,nz) /)
#endif
#endif
#else
  ngp = gp(1:nz)/sum(gp(1:nz))
#ifdef DO_TIME_AVERAGE
  ng = (/ (  sum(ngp(ii:nz) )  , ii = 1,nz) /)
#endif
#endif


#ifdef DO_TIME_AVERAGE 
  nsteps=floor(write_time/average_time)
#endif
 
  nstop=ceiling(total_time/dt)
  wstep=floor(write_time/dt)
#ifdef DO_TIME_AVERAGE
  nstep=0
#endif

  end if
  
  
 
 
 end subroutine
 
 
 subroutine write_params()
  use sys
 
  implicit none
  
  integer :: ii
  character(32) :: desc

  open(unit=10,file='params.txt',status='unknown',form='formatted')
  
  desc='x-direction grid spaces'
  format="(a32,a1,i4)"
  write(10,format) desc, ':', mx
  
  desc='y-direction grid spaces'
  format="(a32,a1,i4)"
  write(10,format) desc, ':', my
  
  desc='z-direction # of layers'
  format="(a32,a1,i4)"
  write(10,format) desc, ':', nz
  
  desc='x-direction grid width'
  format="(a32,a1,e23.16)"
  write(10,format) desc, ':', dx
  
  desc='y-direction grid width'
  format="(a32,a1,e23.16)"
  write(10,format) desc, ':', dy
  
  desc='no slip'
  format="(a32,a1,l3)"
  write(10,format) desc, ':', noslip
  
  desc='planetary rotational frequency'
  format="(a32,a1,e23.16)"
  write(10,format) desc, ':', omega
  
  desc='gravitational acceleration'
  format="(a32,a1,e23.16)"
  write(10,format) desc, ':', g
  
  desc='reduced gravitational accel'
  write(format,"(a8,i1,a7)") '(a32,a1,',ubound(gp,1)-lbound(gp,1)+1,'e23.16)'
  write(10,format) desc, ':', (gp(ii) , ii=lbound(gp,1),ubound(gp,1))
  
  desc='dimensional coriolis'
  format="(a32,a1,e23.16)"
  write(10,format) desc, ':', f0
  
  desc='dimensional depth'
  format="(a32,a1,e23.16)"
  write(10,format) desc, ':', h0
  
  desc='dimensional deformation radius'
  format="(a32,a1,e23.16)"
  write(10,format) desc, ':', ld

  desc='timestep'
  format="(a32,a1,e23.16)"
  write(10,format) desc, ':', dt

  desc='runtime'
  format="(a32,a1,e23.16)"
  write(10,format) desc, ':', total_time

  desc='writetime'
  format="(a32,a1,e23.16)"
  write(10,format) desc, ':', write_time

  desc='bottom friction'
  format="(a32,a1,e23.16)"
  write(10,format) desc, ':', bf

  desc='bottom friction velocity'
  format="(a32,a1,e23.16)"
  write(10,format) desc, ':', umax

  desc='bottom friction polynomial'
  format="(a32,a1,i4)"
  write(10,format) desc, ':', bfricpoly

  desc='wind stress'
  format="(a32,a1,e23.16)"
  write(10,format) desc, ':', tau

  desc='viscosity'
  format="(a32,a1,e23.16)"
  write(10,format) desc, ':', cvar

#ifdef DO_TIME_AVERAGE
  desc='initial time sum itereation'
  format="(a32,a1,i16)"
  write(10,format) desc, ':', timeav_count
#endif

 
 end subroutine
 
  
 
 
 
 function modeltime(n)
   
  implicit none
   
  integer, intent(in) :: n
  real(kind=8) :: t
  integer :: yrs, days, hours
  character(26) :: modeltime
   
   t=(dt/f0)*dble(n)
   yrs=floor(t/3.15576d7)
   t=t-dble(yrs)*3.15576d7
   days=floor(t/8.64d4)
   t=t-dble(days)*8.64d4
   hours=floor(t/3.6d3)
   
  
  write(modeltime, "(i4.3,a5,i4.3,a5,i3.2,a5)") &
     yrs, ' yrs,', &
     days, ' dys,', &
     hours, ' hrs.'

  end function
 
end module
