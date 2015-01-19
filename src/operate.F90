
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
  
  dat%z=field%z
  
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
  
  dat%z=field
  
 end subroutine
 
 pure subroutine var_arrayass(dat,field)
  use global
  
  implicit none
  
  type(var), intent(inout) :: dat
  real(kind=db), intent(in) :: field(dat%p%lx+1:dat%p%lx+dat%p%nx,   &
                 dat%p%ly+1:dat%p%ly+dat%p%ny)
  
  dat%z=field
  
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
 
 contains
 
  pure function sGx(a) result (b)
   use global
   use params
   
   type(var), intent(in) :: a
   type(var) :: b
   
   
   b%z=(a%bz(a%p%lx+1:a%p%lx+a%p%nx+1,a%p%ly+1:a%p%ly+a%p%ny) - &
        a%bz(a%p%lx:a%p%lx+a%p%nx,a%p%ly+1:a%p%ly+a%p%ny))/dx
 
  end function




end module
