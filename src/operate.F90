
module overload

#include "include.h"


 implicit none
 
 
 
 
 interface assignment (=)
   module procedure var_selfass
   module procedure var_ass
   module procedure var_arrayass
 end interface
 
 interface operator (+)
  module procedure var_add
  module procedure realvar_add
  module procedure realarrayvar_add
 end interface
 
 interface operator (-)
  module procedure var_sub
 end interface
 
 interface operator (*)
  module procedure var_dot
  module procedure realvar_dot
  module procedure realarrayvar_dot
 end interface
 
 
 contains
 
 !-------------------------------------------------------------------------!
 
 ! OPERATIONS
 
 !-------------------------------------------------------------------------!
 
 
 
 
 
  pure function var_add(v1,v2) result (v3)
   use global
   
   implicit none
  
   class(var), intent(in) :: v1,v2
   real(kind=db)  :: v3(v1%p%lx+1:v1%p%lx+v1%p%nx,v1%p%ly+1:v1%p%ly+v1%p%ny)
!   type(var)  :: v3
   
   v3 = v1%bz + v2%bz
   
  end function
 
  pure function realvar_add(v1,v2) result (v3)
   use global
   
   implicit none
  
   real(kind=db), intent(in) :: v1
   class(var), intent(in) :: v2
   real(kind=db)  :: v3(v2%p%lx+1:v2%p%lx+v2%p%nx,v2%p%ly+1:v2%p%ly+v2%p%ny)
!    type(var)  :: v3
   
   v3 = v1 + v2%bz
   
  end function
 
  pure function realarrayvar_add(v1,v2) result (v3)
   use global
   
   implicit none
  
   class(var), intent(in) :: v2
   real(kind=db), intent(in) :: v1(v2%p%lx+1:v2%p%lx+v2%p%nx,v2%p%ly+1:v2%p%ly+v2%p%ny)
   real(kind=db)  :: v3(v2%p%lx+1:v2%p%lx+v2%p%nx,v2%p%ly+1:v2%p%ly+v2%p%ny)
!    type(var)  :: v3
   
   v3 = v1 + v2%bz
   
  end function
 
 
 
 
 
  pure function var_sub(v1,v2) result (v3)
   use global
   
   implicit none
  
   class(var), intent(in) :: v1,v2
   real(kind=db)  :: v3(v1%p%lx+1:v1%p%lx+v1%p%nx,v1%p%ly+1:v1%p%ly+v1%p%ny)
!   type(var)  :: v3
   
   v3 = v1%bz - v2%bz
   
  end function
 
  
  
 
 
 
 
  pure function var_dot(v1,v2) result (v3)
   use global
   
   implicit none
  
   class(var), intent(in) :: v1,v2
   real(kind=db)  :: v3(v1%p%lx+1:v1%p%lx+v1%p%nx,v1%p%ly+1:v1%p%ly+v1%p%ny)
!   type(var)  :: v3
   
   v3 = v1%bz * v2%bz
   
  end function
 
  pure function realvar_dot(v1,v2) result (v3)
   use global
   
   implicit none
  
   real(kind=db), intent(in) :: v1
   class(var), intent(in) :: v2
   real(kind=db)  :: v3(v2%p%lx+1:v2%p%lx+v2%p%nx,v2%p%ly+1:v2%p%ly+v2%p%ny)
!    type(var)  :: v3
   
   v3 = v1 * v2%bz
   
  end function
 
  pure function realarrayvar_dot(v1,v2) result (v3)
   use global
   
   implicit none
  
   class(var), intent(in) :: v2
   real(kind=db), intent(in) :: v1(v2%p%lx+1:v2%p%lx+v2%p%nx,v2%p%ly+1:v2%p%ly+v2%p%ny)
   real(kind=db)  :: v3(v2%p%lx+1:v2%p%lx+v2%p%nx,v2%p%ly+1:v2%p%ly+v2%p%ny)
!    type(var)  :: v3
   
   v3 = v1 * v2%bz
   
  end function
  
  
  
 
 
 
 elemental subroutine var_selfass(dat,field)
  use global
  
  implicit none
  
  class(var), intent(inout) :: dat
  class(var), intent(in) :: field
  
  dat%bz=field%bz
  
 end subroutine
 
 elemental subroutine var_ass(dat,field)
  use global
  
  implicit none
  
  class(var), intent(inout) :: dat
  real(kind=db), intent(in) :: field
  
  dat%bz=field
  
 end subroutine
 
 pure subroutine var_arrayass(dat,field)
  use global
  
  implicit none
  
  class(var), intent(inout) :: dat
  real(kind=db), intent(in) :: field(dat%p%lx+1:dat%p%lx+dat%p%nx,dat%p%ly+1:dat%p%ly+dat%p%ny)
  
  dat%bz=field
  
 end subroutine
 
 
 
 
end module




module grid_operate

#include "include.h"

 use overload

 implicit none

 interface Gx
  module procedure :: rGx, hGx, uGx, vGx, zGx
 end interface

 interface Gy
  module procedure :: rGy, hGy, uGy, vGy, zGy
 end interface

 interface Ax
  module procedure :: rAx, hAx, uAx, vAx, zAx
 end interface

 interface Ay
  module procedure :: rAy, hAy, uAy, vAy, zAy
 end interface
  

 
 contains
 
 
 
  function sGxx(a) result (b)
   use global
   use params
   use sync
   
   implicit none
   
   class(var), intent(inout) :: a
   real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx,a%p%ly+1:a%p%ly+a%p%ny)
   
   call end_sync(a)
   b=(a%z(a%p%lx+2:a%p%lx+a%p%nx+1,a%p%ly+1:a%p%ly+a%p%ny) - &
        2.0d0*a%bz + a%z(a%p%lx:a%p%lx+a%p%nx-1,a%p%ly+1:a%p%ly+a%p%ny))/(dx)**2
 
  end function
 
  function sGyy(a) result (b)
   use global
   use params
   use sync
   
   implicit none
   
   class(var), intent(inout) :: a
   real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx,a%p%ly+1:a%p%ly+a%p%ny)
   
   
   call end_sync(a)
   b=(a%z(a%p%lx+1:a%p%lx+a%p%nx,a%p%ly+2:a%p%ly+a%p%ny+1) - &
        2.0d0*a%bz + a%z(a%p%lx+1:a%p%lx+a%p%nx,a%p%ly:a%p%ly+a%p%ny-1))/(dy)**2
 
  end function
  
  
 
  pure function sGx(a) result (b)
   use global
   use params
   
   implicit none
   
   class(var), intent(in) :: a
   real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx+1,a%p%ly+1:a%p%ly+a%p%ny)
   
   b=(a%z(a%p%lx+1:a%p%lx+a%p%nx+1,a%p%ly+1:a%p%ly+a%p%ny) - &
        a%z(a%p%lx:a%p%lx+a%p%nx,a%p%ly+1:a%p%ly+a%p%ny))/dx
 
  end function
 
  pure function sGy(a) result (b)
   use global
   use params
   
   implicit none
   
   class(var), intent(in) :: a
   real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx,a%p%ly+1:a%p%ly+a%p%ny+1)
   
   
   b=(a%z(a%p%lx+1:a%p%lx+a%p%nx,a%p%ly+1:a%p%ly+a%p%ny+1) - &
        a%z(a%p%lx+1:a%p%lx+a%p%nx,a%p%ly:a%p%ly+a%p%ny))/dy
 
  end function
 
  pure function sAx(a) result (b)
   use global
   use params
   
   implicit none
   
   class(var), intent(in) :: a
   real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx+1,a%p%ly+1:a%p%ly+a%p%ny)
   
   
   b=(a%z(a%p%lx+1:a%p%lx+a%p%nx+1,a%p%ly+1:a%p%ly+a%p%ny) + &
        a%z(a%p%lx:a%p%lx+a%p%nx,a%p%ly+1:a%p%ly+a%p%ny))/2.0d0
 
  end function
 
  pure function sAy(a) result (b)
   use global
   use params
   
   implicit none
   
   class(var), intent(in) :: a
   real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx,a%p%ly+1:a%p%ly+a%p%ny+1)
   
   
   b=(a%z(a%p%lx+1:a%p%lx+a%p%nx,a%p%ly+1:a%p%ly+a%p%ny+1) + &
        a%z(a%p%lx+1:a%p%lx+a%p%nx,a%p%ly:a%p%ly+a%p%ny))/2.0d0
 
  end function
  
  
  
  
  
 
  pure function tGx(a) result (b)
   use global
   use params
   
   implicit none
   
   class(var), intent(in) :: a
   real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx-1,a%p%ly+1:a%p%ly+a%p%ny)
   
   b=rGx(a%bz)
 
  end function
 
  pure function tGy(a) result (b)
   use global
   use params
   
   implicit none
   
   class(var), intent(in) :: a
   real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx,a%p%ly+1:a%p%ly+a%p%ny-1)
   
   
   b=rGy(a%bz)
 
  end function
 
  pure function tAx(a) result (b)
   use global
   use params
   
   implicit none
   
   class(var), intent(in) :: a
   real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx-1,a%p%ly+1:a%p%ly+a%p%ny)
   
   b=rAx(a%bz)
 
  end function
 
  pure function tAy(a) result (b)
   use global
   use params
   
   class(var), intent(in) :: a
   real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx,a%p%ly+1:a%p%ly+a%p%ny-1)
   
   b=rAy(a%bz)
 
  end function
  
  
  
  
 
  pure function rGx(a) result (b)
   use global
   use params
   
   implicit none
   
   real(kind=db), intent(in) :: a(0:,:)
   real(kind=db) :: b(ubound(a,1),ubound(a,2))
   
   b=(a(1:ubound(a,1),:) - &
        a(0:ubound(a,1)-1,:))/dx
 
  end function
  
  pure function rGy(a) result (b)
   use global
   use params
   
   implicit none
   
   real(kind=db), intent(in) :: a(:,0:)
   real(kind=db) :: b(ubound(a,1),ubound(a,2))
   
   b=(a(:,1:ubound(a,2)) - &
        a(:,0:ubound(a,2)-1))/dy
 
  end function
  
  
 
  pure function rAx(a) result (b)
   use global
   use params
   
   implicit none
   
   real(kind=db), intent(in) :: a(0:,:)
   real(kind=db) :: b(ubound(a,1),ubound(a,2))
   
   b=(a(1:ubound(a,1),:) + &
        a(0:ubound(a,1)-1,:))/2.0d0
        
  end function
  
  pure function rAy(a) result (b)
   use global
   use params
   
   implicit none
   
   real(kind=db), intent(in) :: a(:,0:)
   real(kind=db) :: b(ubound(a,1),ubound(a,2))
   
   b=(a(:,1:ubound(a,2)) + &
        a(:,0:ubound(a,2)-1))/2.0d0
 
  end function
  
  
  function hGx(a) result (b)
   use global
   use params
   use sync
   
   implicit none
   
   type(hvar), intent(inout) :: a
   real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx+1,a%p%ly+1:a%p%ly+a%p%ny)
   
   call end_sync(a)
   b=sGx(a)
   
  end function
  
  
  function hGy(a) result (b)
   use global
   use params
   use sync
   
   implicit none
   
   type(hvar), intent(inout) :: a
   real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx,a%p%ly+1:a%p%ly+a%p%ny+1)
   
   call end_sync(a)
   b=sGy(a)
   
  end function
  
  function hAx(a) result (b)
   use global
   use params
   use sync
   
   implicit none
   
   type(hvar), intent(inout) :: a
   real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx+1,a%p%ly+1:a%p%ly+a%p%ny)
   
   call end_sync(a)
   b=sAx(a)
   
  end function
  
  
  function hAy(a) result (b)
   use global
   use params
   use sync
   
   implicit none
   
   type(hvar), intent(inout) :: a
   real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx,a%p%ly+1:a%p%ly+a%p%ny+1)
   
   call end_sync(a)
   b=sAy(a)
   
  end function
  
  
  
  
  
  
  
  pure function uGx(a) result (b)
   use global
   use params
   
   implicit none
   
   type(uvar), intent(in) :: a
   real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx-1,a%p%ly+1:a%p%ly+a%p%ny)
   
   b=tGx(a)
   
  end function
  
  
  function uGy(a) result (b)
   use global
   use params
   use sync
   
   implicit none
   
   type(uvar), intent(inout) :: a
   real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx,a%p%ly+1:a%p%ly+a%p%ny+1)
   
   call end_sync(a)
   b=sGy(a)
   
  end function
  
  pure function uAx(a) result (b)
   use global
   use params
   
   implicit none
   
   type(uvar), intent(in) :: a
   real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx-1,a%p%ly+1:a%p%ly+a%p%ny)
   
   b=tAx(a)
   
  end function
  
  
  function uAy(a) result (b)
   use global
   use params
   use sync
   
   implicit none
   
   type(uvar), intent(inout) :: a
   real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx,a%p%ly+1:a%p%ly+a%p%ny+1)
   
   call end_sync(a)
   b=sAy(a)
   
  end function 
  
  
  
  
  
  
  
  
  function vGx(a) result (b)
   use global
   use params
   use sync
   
   implicit none
   
   type(vvar), intent(inout) :: a
   real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx+1,a%p%ly+1:a%p%ly+a%p%ny)
   
   call end_sync(a)
   b=sGx(a)
   
  end function
  
  
  pure function vGy(a) result (b)
   use global
   use params
   
   implicit none
   
   type(vvar), intent(in) :: a
   real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx,a%p%ly+1:a%p%ly+a%p%ny-1)
   
   b=tGy(a)
   
  end function
  
  function vAx(a) result (b)
   use global
   use params
   use sync
   
   implicit none
   
   type(vvar), intent(inout) :: a
   real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx+1,a%p%ly+1:a%p%ly+a%p%ny)
   
   call end_sync(a)
   b=sAx(a)
   
  end function
  
  
  pure function vAy(a) result (b)
   use global
   use params
   
   implicit none
   
   type(vvar), intent(in) :: a
   real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx,a%p%ly+1:a%p%ly+a%p%ny-1)
   
   b=tAy(a)
   
  end function
  
  
  
  
  
  
  
  
  pure function zGx(a) result (b)
   use global
   use params
   
   implicit none
   
   type(zvar), intent(in) :: a
   real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx-1,a%p%ly+1:a%p%ly+a%p%ny)
   
   b=tGx(a)
   
  end function
  
  
  pure function zGy(a) result (b)
   use global
   use params
   
   implicit none
   
   type(zvar), intent(in) :: a
   real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx,a%p%ly+1:a%p%ly+a%p%ny-1)
   
   b=tGy(a)
   
  end function
  
  pure function zAx(a) result (b)
   use global
   use params
   
   implicit none
   
   type(zvar), intent(in) :: a
   real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx-1,a%p%ly+1:a%p%ly+a%p%ny)
   
   b=tAx(a)
   
  end function
  
  
  pure function zAy(a) result (b)
   use global
   use params
   
   implicit none
   
   type(zvar), intent(in) :: a
   real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx,a%p%ly+1:a%p%ly+a%p%ny-1)
   
   b=tAy(a)
   
  end function
  
  
  function x_dist(a,x) result (p)
   use sys
   use global
   use params
   
   implicit none
   
   class(var), intent(in) :: a
   real(kind=db) :: x
   real(kind=db) :: p(a%p%lx+1:a%p%lx+a%p%nx,a%p%ly+1:a%p%ly+a%p%ny)
   
   p=lim(a%p%x-x,x0)
   
  end function
  
  
  function y_dist(a,y) result (p)
   use sys
   use global
   use params
   
   implicit none
   
   class(var), intent(in) :: a
   real(kind=db) :: y
   real(kind=db) :: p(a%p%lx+1:a%p%lx+a%p%nx,a%p%ly+1:a%p%ly+a%p%ny)
   
   p=lim(a%p%y-y,y0)
   
  end function
  
  
 
 
!--------------------------------------------------------------------------!

! TIMESTEPPING

!--------------------------------------------------------------------------!
 
 
 pure function fe(a) result(b)
  use params
  use global
  
  implicit none
 
  class(var), intent(in) :: a
  real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx,a%p%ly+1:a%p%ly+a%p%ny)
  
  b=a%bz + 0.5d0*dt*a%tend1%bz
  
  return
  
 end function fe
 
 
 
 
 
  

 
 pure function rk2(a) result(b)
  use params
  use global
  
  implicit none
 
  class(var), intent(in) :: a
  real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx,a%p%ly+1:a%p%ly+a%p%ny)
 
  b = a%bz + dt*a%tend2%bz
  
  return
  
 end function rk2
 

 
 
 
 
 
 pure function ab2(a) result(b)
  use params
  use global
  
  implicit none
 
  class(var), intent(in) :: a
  real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx,a%p%ly+1:a%p%ly+a%p%ny)
  
 
  b=a%bz+dt*(1.5d0*a%tend1%bz-0.5d0*a%tend2%bz)
 
  return
  
 end function ab2
 

 
 
 
 
 pure function ab3(a) result(b)
  use params
  use global
  
  implicit none
 
  class(var), intent(in) :: a
  real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx,a%p%ly+1:a%p%ly+a%p%ny)
 
   b=a%bz+dt*((23.0d0/12.0d0)*a%tend1%bz-(4.0d0/3.0d0)*a%tend2%bz+                  &
       (5.0d0/12.0d0)*a%tend3%bz)
 
  return
  
 end function ab3
 
 
 
 
 
!--------------------------------------------------------------------------!

! TENDANCIES

!--------------------------------------------------------------------------!


 
 subroutine tend_h(n)
  use global
  use params
  use variables
  use sync
  
  implicit none
 
  integer :: i
  integer, intent(in) :: n
  real(kind=8), pointer :: tend(:,:)
  
  
   do i=1,nz
    if (n == 1) then
     tend => h(i)%tend2%bz
    else
     tend => h(i)%tend1%bz
    end if
    tend=-Gx(u(i)%bz*h_u(i)%bz)-Gy(v(i)%bz*h_v(i)%bz)
   end do
   
   nullify(tend)
  
  
 end subroutine
 
 subroutine tend_u(n)
  use global
  use params
  use variables
  use sync
  
  implicit none
 
  integer :: i
  integer, intent(in) :: n
  real(kind=8), pointer :: tend(:,:)
  
  
   do i=1,nz
    if (n == 1) then
     tend => u(i)%tend2%bz
    else
     tend => u(i)%tend1%bz
    end if
   call end_sync(v(i))
   call end_sync(m(i))
   call end_sync(ke(i))
    tend=0.0d0    &
     +smagu(i)%bz  &
     +Ay(f%bz+zeta(i)%bz)   &
      *Ay(Ax(v(i)))    &
     -Gx(m(i))    &
     -Gx(ke(i))     &
     +0.0d0
    if (i == nz) tend = tend + utau%bz*merge(0.0d0,1.0d0/h_u(nz)%bz,  &
        ( h_u(nz)%bz == 0.0d0)) 
    if (i == 1) tend = tend - bfricu%bz*merge(0.0d0,   &
          1.0d0/h_u(1)%bz,(h_u(1)%bz == 0.0d0))
   end do
   
   nullify(tend)
   
  
 end subroutine
 
 subroutine tend_v(n)
  use global
  use params
  use variables
  use sync
  
  implicit none
 
  integer :: i
  integer, intent(in) :: n
  real(kind=8), pointer :: tend(:,:)
  
  
   do i=1,nz
    if (n == 1) then
     tend => v(i)%tend2%bz
    else
     tend => v(i)%tend1%bz
    end if
   call end_sync(u(i))
   call end_sync(m(i))
   call end_sync(ke(i))
    tend=0.0d0    &
     +smagv(i)%bz   &
     -Ax(f%bz+zeta(i)%bz)   &
     *Ax(Ay(u(i)))    &
     -Gy(m(i))   &
     -Gy(ke(i))    &
     +0.0d0
    if (i == nz) tend = tend + vtau%bz*merge(0.0d0,1.0d0/h_v(nz)%bz,  &
         (h_v(nz)%bz == 0.0d0))
    if (i == 1) tend = tend - bfricv%bz*merge(0.0d0,   &
          1.0d0/h_v(1)%bz,(h_v(1)%bz == 0.0d0))
  end do
   
   nullify(tend)
  
  
 end subroutine


end module




module operations
 use global
 use overload
 
 implicit none
 
 contains
 
 subroutine depth(h,s,d)
  use params
  
  type(hvar), intent(in) :: h(1:nz), s
  type(hvar), intent(inout) :: d(0:nz)
  integer :: k
  
  d(0)=s
  do k=1,nz
   d(k)=d(k-1)+h(k)
  end do
  
 end subroutine
 
 subroutine mont(d,m)
  use params
  
  type(hvar), intent(in) :: d(0:nz)
  type(hvar), intent(inout) :: m(1:nz)
  integer :: k
  
#ifdef ALLOW_RIGID_LID
  m(nz)=0.0d0
#else
  m(nz)=ngp(nz)*d(nz)
#endif
  do k=nz-1,1,-1
   m(k)=ngp(k)*d(k)+m(k+1)
  end do
  
 end subroutine
  
  
  
 end module
 
 

