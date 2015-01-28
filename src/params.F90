module params

#include "include.h"

 use sys
 
 implicit none

 integer, parameter :: mx=100, my=150
 integer, parameter :: nz=3
 real(kind=db), parameter :: dx=0.1d0, dy=0.1d0
#ifdef DOUBLE_PRECISION
 real(kind=db), parameter :: x0=dx*dble(mx), y0=dy*dble(my)
#else
 real(kind=db), parameter :: x0=dx*real(mx), y0=dy*real(my)
#endif
 
 
 integer :: lx, ly, nx, ny, north, south, east, west
 
 logical :: noslip = .true.
 
 
 
 
 real (kind=db), parameter :: pi=4.0d0*atan(1.0d0)
 real (kind=db), parameter :: omega=2.0d0*pi/8.64d4
 real (kind=db), parameter :: g=9.81d0
 
  ! Degree/radian conversion parameters
 real (kind=db), parameter :: deg2rad=pi/180.0d0, rad2deg=180.0d0/pi
 
 integer :: ii

#ifndef ALLOW_STATIC_LAYER
#ifdef ALLOW_RIGID_LID
 real (kind=db), parameter :: gp(nz-1) = (/ 3.0d-3 , 7.0d-3 /)
 real (kind=db), parameter :: ngp(nz) = (/ gp(1:nz-1) , g /)/sum(gp(1:nz-1))
#ifdef DO_TIME_AVERAGE
 real (kind=db), parameter :: ng(nz) =  &
         (/ (  sum(ngp(ii:nz) )  , ii = 1,nz) /)
#endif
#else
 real (kind=db), parameter :: gp(nz) = (/ 3.0d-3 , 7.0d-3, 9.81d0 /)
 real (kind=db), parameter :: ngp(nz) = gp(1:nz)/sum(gp(1:nz-1))
#ifdef DO_TIME_AVERAGE
 real (kind=db), parameter :: ng(nz) =  &
         (/ (  sum(ngp(ii:nz) )  , ii = 1,nz) /)
#endif
#endif
#else
 real (kind=db), parameter :: gp(nz) = (/ 3.0d-3 , 7.0d-3 /)
 real (kind=db), parameter :: ngp(nz) = gp(1:nz)/sum(gp(1:nz))
#ifdef DO_TIME_AVERAGE
 real (kind=db), parameter :: ng(nz) =  &
         (/ (  sum(ngp(ii:nz) )  , ii = 1,nz) /)
#endif
#endif

 real (kind=db), parameter :: f0=2.0d0*omega*sin(deg2rad*70.0d0)
 real (kind=db), parameter :: h0=1.0d3
#ifndef ALLOW_STATIC_LAYER
 real (kind=db), parameter :: ld=sqrt(sum(gp(1:nz-1))*h0)/f0
#else
 real (kind=db), parameter :: ld=sqrt(sum(gp(1:nz))*h0)/f0
#endif




#ifdef ALLOW_RIGID_LID
 real (kind=db), parameter :: maxh=1.25d0!h0
 real (kind=db), parameter :: dt = (0.25d0*dx)/dsqrt(sum(ngp(1:nz-1))*maxh)
#else
 real (kind=db), parameter :: maxh=1.25d0*1.01d0
 real (kind=db), parameter :: dt = (0.25d0*dx)/dsqrt(sum(ngp)*maxh)
#endif






 REAL (KIND=db), PARAMETER :: write_time=1.0d1
 REAL (KIND=db), PARAMETER :: total_time=1.0d4
#ifdef DO_TIME_AVERAGE
 REAL (KIND=db), PARAMETER :: average_time=1.0d0
 
 INTEGER, PARAMETER :: nsteps=floor(write_time/average_time)
#endif
 
 INTEGER, PARAMETER :: nstop=ceiling(total_time/dt)
 INTEGER, PARAMETER :: wstep=floor(write_time/dt)
#ifdef DO_TIME_AVERAGE
 INTEGER, PARAMETER :: nstep=0
#endif


     
 REAL (kind=db), PARAMETER :: bf=2.0d-2
 real (kind=db), parameter :: res=0.0d0/(f0**2*ld*h0) 
 REAL (kind=db), PARAMETER :: tau=2.0e-1/(f0**2*ld*h0*1026.0d0)      
 real (kind=db), parameter :: cvar=(1.0d0/pi**2)  
 
 contains
 
 
 
 function modeltime(n)
   
  implicit none
   
  integer, intent(in) :: n
  real(kind=8) :: t
  integer :: yrs, days, hours
  character(26) :: modeltime
   
   t=(dt/f0)*dble(n)
   yrs=floor(t/3.15576e7)
   t=t-dble(yrs)*3.15576e7
   days=floor(t/8.64e4)
   t=t-dble(days)*8.64e4
   hours=floor(t/3.6e3)
   
  
  write(modeltime, "(i4.3,a5,i4.3,a5,i3.2,a5)") &
     yrs, ' yrs,', &
     days, ' dys,', &
     hours, ' hrs.'

  end function
 
end module
