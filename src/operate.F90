
module overload

 implicit none
 
 
 
 
 interface assignment (=)
   module procedure var_selfass
   module procedure hvar_selfass
   module procedure uvar_selfass
   module procedure vvar_selfass
   module procedure zvar_selfass
   module procedure var_ass
   module procedure var_arrayass
   module procedure hvar_ass
   module procedure hvar_arrayass
   module procedure uvar_ass
   module procedure uvar_arrayass
   module procedure vvar_ass
   module procedure vvar_arrayass
   module procedure zvar_ass
   module procedure zvar_arrayass
 end interface
 
 interface operator (+)
  module procedure var_add
  module procedure hvar_add
  module procedure uvar_add
  module procedure vvar_add
  module procedure zvar_add
 end interface
 
 interface operator (-)
  module procedure var_sub
  module procedure hvar_sub
  module procedure uvar_sub
  module procedure vvar_sub
  module procedure zvar_sub
 end interface
 
 interface operator (*)
  module procedure var_dot
  module procedure realvar_dot
  module procedure realhvar_dot
  module procedure hvar_dot
  module procedure uvar_dot
  module procedure vvar_dot
  module procedure zvar_dot
 end interface
 
 
 contains
 
 !-------------------------------------------------------------------------!
 
 ! OPERATIONS
 
 !-------------------------------------------------------------------------!
 
 
 
 
 
  pure function var_add(v1,v2) result (v3)
   use global
   
   implicit none
  
   type(var), intent(in) :: v1,v2
   type(var)  :: v3
   
   v3%z = v1%z + v2%z
   
  end function
 
  elemental function hvar_add(v1,v2) result (v3)
   use global
   
   implicit none
  
   type(hvar), intent(in) :: v1,v2
   type(hvar) :: v3
   
   v3%z = v1%z + v2%z
   
  end function
 
  elemental function uvar_add(v1,v2) result (v3)
   use global
   
   implicit none
  
   type(uvar), intent(in) :: v1,v2
   type(uvar) :: v3
   
   v3%z = v1%z + v2%z
   
  end function
 
  elemental function vvar_add(v1,v2) result (v3)
   use global
   
   implicit none
  
   type(vvar), intent(in) :: v1,v2
   type(vvar) :: v3
   
   v3%z = v1%z + v2%z
   
  end function
 
  elemental function zvar_add(v1,v2) result (v3)
   use global
   
   implicit none
  
   type(zvar), intent(in) :: v1,v2
   type(zvar) :: v3
   
   v3%z = v1%z + v2%z
   
  end function
 
 
 
 
  pure function var_sub(v1,v2) result (v3)
   use global
   
   implicit none
  
   type(var), intent(in) :: v1,v2
   type(var)  :: v3
   
   v3%z = v1%z - v2%z
   
  end function
 
  elemental function hvar_sub(v1,v2) result (v3)
   use global
   
   implicit none
  
   type(hvar), intent(in) :: v1,v2
   type(hvar) :: v3
   
   v3%z = v1%z - v2%z
   
  end function
 
  elemental function uvar_sub(v1,v2) result (v3)
   use global
   
   implicit none
  
   type(uvar), intent(in) :: v1,v2
   type(uvar) :: v3
   
   v3%z = v1%z - v2%z
   
  end function
 
  elemental function vvar_sub(v1,v2) result (v3)
   use global
   
   implicit none
  
   type(vvar), intent(in) :: v1,v2
   type(vvar) :: v3
   
   v3%z = v1%z - v2%z
   
  end function
 
  elemental function zvar_sub(v1,v2) result (v3)
   use global
   
   implicit none
  
   type(zvar), intent(in) :: v1,v2
   type(zvar) :: v3
   
   v3%z = v1%z - v2%z
   
  end function
  
  
 
 
 
 
  pure function var_dot(v1,v2) result (v3)
   use global
   
   implicit none
  
   type(var), intent(in) :: v1,v2
   type(var)  :: v3
   
   v3%z = v1%z * v2%z
   
  end function
 
  pure function realvar_dot(v1,v2) result (v3)
   use global
   
   implicit none
  
   real(kind=db), intent(in) :: v1
   type(var), intent(in) :: v2
   type(var)  :: v3
   
   v3%z = v1 * v2%z
   
  end function
 
  elemental function hvar_dot(v1,v2) result (v3)
   use global
   
   implicit none
  
   type(hvar), intent(in) :: v1,v2
   type(hvar) :: v3
   
   v3%z = v1%z * v2%z
   
  end function
 
  elemental function realhvar_dot(v1,v2) result (v3)
   use global
   
   implicit none
  
   real(kind=db), intent(in) :: v1
   type(hvar), intent(in) :: v2
   type(hvar)  :: v3
   
   v3%z = v1 * v2%z
   
  end function
 
  elemental function uvar_dot(v1,v2) result (v3)
   use global
   
   implicit none
  
   type(uvar), intent(in) :: v1,v2
   type(uvar) :: v3
   
   v3%z = v1%z * v2%z
   
  end function
 
  elemental function vvar_dot(v1,v2) result (v3)
   use global
   
   implicit none
  
   type(vvar), intent(in) :: v1,v2
   type(vvar) :: v3
   
   v3%z = v1%z * v2%z
   
  end function
 
  elemental function zvar_dot(v1,v2) result (v3)
   use global
   
   implicit none
  
   type(zvar), intent(in) :: v1,v2
   type(zvar) :: v3
   
   v3%z = v1%z * v2%z
   
  end function
  
  
  
  
 
 
 
 pure subroutine var_selfass(dat,field)
  use global
  
  implicit none
  
  type(var), intent(inout) :: dat
  type(var), intent(in) :: field
  
  dat%bz=field%bz
  
 end subroutine
 
 elemental subroutine hvar_selfass(dat,field)
  use global
  
  implicit none
  
  type(hvar), intent(inout) :: dat
  type(hvar), intent(in) :: field
  
  dat%z=field%z
  
 end subroutine
 
 elemental subroutine uvar_selfass(dat,field)
  use global
  
  implicit none
  
  type(uvar), intent(inout) :: dat
  type(uvar), intent(in) :: field
  
  dat%z=field%z
  
 end subroutine
 
 elemental subroutine vvar_selfass(dat,field)
  use global
  
  implicit none
  
  type(vvar), intent(inout) :: dat
  type(vvar), intent(in) :: field
  
  dat%z=field%z
  
 end subroutine
 
 elemental subroutine zvar_selfass(dat,field)
  use global
  
  implicit none
  
  type(zvar), intent(inout) :: dat
  type(zvar), intent(in) :: field
  
  dat%z=field%z
  
 end subroutine
 
 pure subroutine var_ass(dat,field)
  use global
  
  implicit none
  
  type(var), intent(inout) :: dat
  real(kind=db), intent(in) :: field
  
  dat%bz=field
  
 end subroutine
 
 pure subroutine var_arrayass(dat,field)
  use global
  
  implicit none
  
  type(var), intent(inout) :: dat
  real(kind=db), intent(in) :: field(:,:)
  
  dat%bz=field
  
 end subroutine
 
 elemental subroutine hvar_ass(dat,field)
  use global
  
  implicit none
  
  type(hvar), intent(inout) :: dat
  real(kind=db), intent(in) :: field
  
  dat%z=field
  
 end subroutine
 
 pure subroutine hvar_arrayass(dat,field)
  use global
  
  implicit none
  
  type(hvar), intent(inout) :: dat
  real(kind=db), intent(in) :: field(:,:)
  
  dat%z=field
  
 end subroutine
 
 elemental subroutine uvar_ass(dat,field)
  use global
  
  implicit none
  
  type(uvar), intent(inout) :: dat
  real(kind=db), intent(in) :: field
  
  dat%z=field
  
 end subroutine
 
 pure subroutine uvar_arrayass(dat,field)
  use global
  
  implicit none
  
  type(uvar), intent(inout) :: dat
  real(kind=db), intent(in) :: field(:,:)
  
  dat%z=field
  
 end subroutine
 
 elemental subroutine vvar_ass(dat,field)
  use global
  
  implicit none
  
  type(vvar), intent(inout) :: dat
  real(kind=db), intent(in) :: field
  
  dat%z=field
  
 end subroutine
 
 pure subroutine vvar_arrayass(dat,field)
  use global
  
  implicit none
  
  type(vvar), intent(inout) :: dat
  real(kind=db), intent(in) :: field(:,:)
  
  dat%z=field
  
 end subroutine
 
 elemental subroutine zvar_ass(dat,field)
  use global
  
  implicit none
  
  type(zvar), intent(inout) :: dat
  real(kind=db), intent(in) :: field
  
  dat%z=field
  
 end subroutine
 
 pure subroutine zvar_arrayass(dat,field)
  use global
  
  implicit none
  
  type(zvar), intent(inout) :: dat
  real(kind=db), intent(in) :: field(:,:)
  
  dat%z=field
  
 end subroutine
 
 
 
end module




module grid_operate
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
 
  pure function sGx(a) result (b)
   use global
   use params
   
   implicit none
   
   type(var), intent(in) :: a
   real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx+1,a%p%ly+1:a%p%ly+a%p%ny)
   
   b=(a%z(a%p%lx+1:a%p%lx+a%p%nx+1,a%p%ly+1:a%p%ly+a%p%ny) - &
        a%z(a%p%lx:a%p%lx+a%p%nx,a%p%ly+1:a%p%ly+a%p%ny))/dx
 
  end function
 
  pure function sGy(a) result (b)
   use global
   use params
   
   implicit none
   
   type(var), intent(in) :: a
   real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx,a%p%ly+1:a%p%ly+a%p%ny+1)
   
   
   b=(a%z(a%p%lx+1:a%p%lx+a%p%nx,a%p%ly+1:a%p%ly+a%p%ny+1) - &
        a%z(a%p%lx+1:a%p%lx+a%p%nx,a%p%ly:a%p%ly+a%p%ny))/dy
 
  end function
 
  pure function sAx(a) result (b)
   use global
   use params
   
   implicit none
   
   type(var), intent(in) :: a
   real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx+1,a%p%ly+1:a%p%ly+a%p%ny)
   
   
   b=(a%z(a%p%lx+1:a%p%lx+a%p%nx+1,a%p%ly+1:a%p%ly+a%p%ny) + &
        a%z(a%p%lx:a%p%lx+a%p%nx,a%p%ly+1:a%p%ly+a%p%ny))/2.0d0
 
  end function
 
  pure function sAy(a) result (b)
   use global
   use params
   
   implicit none
   
   type(var), intent(in) :: a
   real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx,a%p%ly+1:a%p%ly+a%p%ny+1)
   
   
   b=(a%z(a%p%lx+1:a%p%lx+a%p%nx,a%p%ly+1:a%p%ly+a%p%ny+1) + &
        a%z(a%p%lx+1:a%p%lx+a%p%nx,a%p%ly:a%p%ly+a%p%ny))/2.0d0
 
  end function
  
  
  
  
  
 
  pure function tGx(a) result (b)
   use global
   use params
   
   implicit none
   
   type(var), intent(in) :: a
   real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx-1,a%p%ly+1:a%p%ly+a%p%ny)
   
   b=rGx(a%bz)
 
  end function
 
  pure function tGy(a) result (b)
   use global
   use params
   
   implicit none
   
   type(var), intent(in) :: a
   real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx,a%p%ly+1:a%p%ly+a%p%ny-1)
   
   
   b=rGy(a%bz)
 
  end function
 
  pure function tAx(a) result (b)
   use global
   use params
   
   implicit none
   
   type(var), intent(in) :: a
   real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx-1,a%p%ly+1:a%p%ly+a%p%ny)
   
   b=rAx(a%bz)
 
  end function
 
  pure function tAy(a) result (b)
   use global
   use params
   
   type(var), intent(in) :: a
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
  
  
  pure function hGx(a) result (b)
   use global
   use params
   
   implicit none
   
   type(hvar), intent(in) :: a
   real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx+1,a%p%ly+1:a%p%ly+a%p%ny)
   
   b=sGx(a%z)
   
  end function
  
  
  pure function hGy(a) result (b)
   use global
   use params
   
   implicit none
   
   type(hvar), intent(in) :: a
   real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx,a%p%ly+1:a%p%ly+a%p%ny+1)
   
   b=sGy(a%z)
   
  end function
  
  pure function hAx(a) result (b)
   use global
   use params
   
   implicit none
   
   type(hvar), intent(in) :: a
   real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx+1,a%p%ly+1:a%p%ly+a%p%ny)
   
   b=sAx(a%z)
   
  end function
  
  
  pure function hAy(a) result (b)
   use global
   use params
   
   implicit none
   
   type(hvar), intent(in) :: a
   real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx,a%p%ly+1:a%p%ly+a%p%ny+1)
   
   b=sAy(a%z)
   
  end function
  
  
  
  
  
  
  
  pure function uGx(a) result (b)
   use global
   use params
   
   implicit none
   
   type(uvar), intent(in) :: a
   real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx-1,a%p%ly+1:a%p%ly+a%p%ny)
   
   b=tGx(a%z)
   
  end function
  
  
  pure function uGy(a) result (b)
   use global
   use params
   
   implicit none
   
   type(uvar), intent(in) :: a
   real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx,a%p%ly+1:a%p%ly+a%p%ny+1)
   
   b=sGy(a%z)
   
  end function
  
  pure function uAx(a) result (b)
   use global
   use params
   
   implicit none
   
   type(uvar), intent(in) :: a
   real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx-1,a%p%ly+1:a%p%ly+a%p%ny)
   
   b=tAx(a%z)
   
  end function
  
  
  pure function uAy(a) result (b)
   use global
   use params
   
   implicit none
   
   type(uvar), intent(in) :: a
   real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx,a%p%ly+1:a%p%ly+a%p%ny+1)
   
   b=sAy(a%z)
   
  end function 
  
  
  
  
  
  
  
  
  pure function vGx(a) result (b)
   use global
   use params
   
   implicit none
   
   type(vvar), intent(in) :: a
   real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx+1,a%p%ly+1:a%p%ly+a%p%ny)
   
   b=sGx(a%z)
   
  end function
  
  
  pure function vGy(a) result (b)
   use global
   use params
   
   implicit none
   
   type(vvar), intent(in) :: a
   real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx,a%p%ly+1:a%p%ly+a%p%ny-1)
   
   b=tGy(a%z)
   
  end function
  
  pure function vAx(a) result (b)
   use global
   use params
   
   implicit none
   
   type(vvar), intent(in) :: a
   real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx+1,a%p%ly+1:a%p%ly+a%p%ny)
   
   b=sAx(a%z)
   
  end function
  
  
  pure function vAy(a) result (b)
   use global
   use params
   
   implicit none
   
   type(vvar), intent(in) :: a
   real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx,a%p%ly+1:a%p%ly+a%p%ny-1)
   
   b=tAy(a%z)
   
  end function
  
  
  
  
  
  
  
  
  pure function zGx(a) result (b)
   use global
   use params
   
   implicit none
   
   type(zvar), intent(in) :: a
   real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx-1,a%p%ly+1:a%p%ly+a%p%ny)
   
   b=tGx(a%z)
   
  end function
  
  
  pure function zGy(a) result (b)
   use global
   use params
   
   implicit none
   
   type(zvar), intent(in) :: a
   real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx,a%p%ly+1:a%p%ly+a%p%ny-1)
   
   b=tGy(a%z)
   
  end function
  
  pure function zAx(a) result (b)
   use global
   use params
   
   implicit none
   
   type(zvar), intent(in) :: a
   real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx-1,a%p%ly+1:a%p%ly+a%p%ny)
   
   b=tAx(a%z)
   
  end function
  
  
  pure function zAy(a) result (b)
   use global
   use params
   
   implicit none
   
   type(zvar), intent(in) :: a
   real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx,a%p%ly+1:a%p%ly+a%p%ny-1)
   
   b=tAy(a%z)
   
  end function
  
  


end module
