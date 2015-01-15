module params
 use system
 
 implicit none

 integer, parameter :: mx=2, my=4
 integer, parameter :: nz=3
 real(kind=db), parameter :: dx=0.1d0, dy=0.1d0
#ifdef DOUBLE_PRECISION
 real(kind=db), parameter :: x0=dx*dble(mx), y0=dy*dble(my)
#else
 real(kind=db), parameter :: x0=dx*real(mx), y0=dy*real(my)
#endif
 
 
 integer :: lx, ly, nx, ny, north, south, east, west
 
 logical :: noslip = .true.
 
 
 
 
 real (kind=8), parameter :: pi=4.0d0*atan(1.0d0)
 real (kind=8), parameter :: omega=2.0d0*pi/8.64d4
 real (kind=8), parameter :: g=9.81d0
 
  ! Degree/radian conversion parameters
 real (kind=8), parameter :: deg2rad=pi/180.0d0, rad2deg=180.0d0/pi
 
 integer :: ii

#ifndef ALLOW_STATIC_LAYER
#ifdef ALLOW_RIGID_LID
 real (kind=8), parameter :: gp(nz-1) = (/ 3.0d-3 , 7.0d-3 /)
 real (kind=8), parameter :: ngp(nz) = (/ gp(1:nz-1) , g /)/sum(gp(1:nz-1))
#ifdef DO_TIME_AVERAGE
 real (kind=8), parameter :: ng(nz) =  &
         (/ (  sum(ngp(ii:nz) )  , ii = 1,nz) /)
#endif
#else
 real (kind=8), parameter :: gp(nz) = (/ 3.0d-3 , 7.0d-3, 9.81d0 /)
 real (kind=8), parameter :: ngp(nz) = gp(1:nz)/sum(gp(1:nz-1))
#ifdef DO_TIME_AVERAGE
 real (kind=8), parameter :: ng(nz) =  &
         (/ (  sum(ngp(ii:nz) )  , ii = 1,nz) /)
#endif
#endif
#else
 real (kind=8), parameter :: gp(nz) = (/ 3.0d-3 , 7.0d-3 /)
 real (kind=8), parameter :: ngp(nz) = gp(1:nz)/sum(gp(1:nz))
#ifdef DO_TIME_AVERAGE
 real (kind=8), parameter :: ng(nz) =  &
         (/ (  sum(ngp(ii:nz) )  , ii = 1,nz) /)
#endif
#endif

 real (kind=8), parameter :: f0=2.0d0*omega*sin(deg2rad*70.0d0)
 real (kind=8), parameter :: h0=1.0d3
#ifndef ALLOW_STATIC_LAYER
 real (kind=8), parameter :: ld=sqrt(sum(gp(1:nz-1))*h0)/f0
#else
 real (kind=8), parameter :: ld=sqrt(sum(gp(1:nz))*h0)/f0
#endif
 
end module
