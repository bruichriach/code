module parallel
 
 implicit none

include 'mpif.h'
#include "include.h"
 
 integer, parameter :: max_core=2
 integer, parameter :: ens_num=1
 integer, parameter :: proc_master=0
 integer :: proc_name, proc_num
 integer :: ens_name, ens_images, ens_master
 integer :: max_tag=0
 integer :: mobile_tag=0
 integer, allocatable :: proc_grid(:,:)
#ifdef DOUBLE_PRECISION
 integer, parameter :: mpi_precision=mpi_double_precision
#else
 integer, parameter ::  mpi_precision=mpi_real
#endif
 
 
end module






module sys

#include "include.h"
 
 implicit none
 
#ifdef DOUBLE_PRECISION
 integer, parameter :: db=8
#else
 integer, parameter :: db=4
#endif
 integer :: stat
 real(kind=db), parameter :: unity=1.0d0
 real(kind=db), target :: null_field=1.0d99

 
 type use
  real(kind=db), pointer :: z
 end type
 

 contains
 
 
 
  elemental function lim(x,x0)
  
         implicit none

         real(kind=db), intent(in) :: x, x0
         real(kind=db) :: lim

         lim=x
         do while ((lim < -x0/2.0d0).or.(lim > x0/2.0d0))
          if (lim < -x0/2.0d0) lim=lim+x0
          if (lim > x0/2.0d0) lim=lim-x0
         end do

  end function  
  
 
  function is_factor(m,n)
 
   implicit none
  
   integer, intent(in) :: m, n
   logical :: is_factor
   
    is_factor = ((n/m)*m == n)
    
    
  end function
  
  
  

 function print_time()
  implicit none
  
  real(kind=8) :: time
  integer :: seconds, minutes, hours, days
  character(33) :: print_time
  
  time=cputime()
  days=floor(time/8.64d4)
  time=time-dble(days)*8.64d4
  hours=floor(time/3.6d3)
  time=time-dble(hours)*3.6d3
  minutes=floor(time/6.0d1)
  time=time-dble(minutes)*6.0d1
  seconds=floor(time)
  
  write(print_time, "(i2.1,a5,i3.2,a5,i3.2,a6,i3.2,a6)") &
     days, ' dys,', &
     hours, ' hrs,', &
     minutes, ' mins,', &
     seconds, ' secs.'

 end function
 
 
 function cputime()
  implicit none
  
  real(kind=8) :: cputime
    
  call cpu_time(cputime)
  
 end function
 
 function update(seconds)
  implicit none
  
  logical :: update
  real(kind=8), save :: oldtime=-1.0d0
  real(kind=8), intent(in) :: seconds
  
  if (cputime() - oldtime >= seconds) then
   oldtime=cputime()
   update=.true.
  else
   update=.false.
  end if
   
 end function
 
  
  
  
  function int2real(n)
 
   implicit none
  
   integer, intent(in) :: n
   real(kind=db) :: int2real
 
#ifdef DOUBLE_PRECISION
    int2real=dble(n)
#else
    int2real=real(n)
#endif

  end function
  
  function point_at(a)
   real(kind=db), intent(in), target :: a
   type(use) :: point_at
   
   point_at%z => a
    
  end function
   


end module




module overloading
 implicit none
 
 
 interface assignment (=)
   module procedure usereal_ass
 end interface
 
 interface operator (+)
   module procedure use_add
   module procedure usereal_add
   module procedure realuse_add
 end interface
 
 interface operator (-)
   module procedure use_sub
   module procedure usereal_sub
   module procedure realuse_sub
 end interface
 
 interface operator (*)
   module procedure use_dot
   module procedure usereal_dot
   module procedure realuse_dot
 end interface
 
 interface operator (/)
   module procedure use_div
   module procedure usereal_div
   module procedure realuse_div
 end interface


 
 contains
 
 
 
 
 
  elemental subroutine usereal_ass(v3,v1)
   use sys
  
   implicit none
  
   type(use), intent(in) :: v1
   real(kind=db), intent(out) :: v3
   
   v3 = v1%z
   
  end subroutine
  
  
  
  
  
  
  
  
  
  
  
 
  elemental function use_add(v1,v2) result (v3)
   use sys
  
   implicit none
  
   type(use), intent(in) :: v1,v2
   real(kind=db)  :: v3
   
   v3 = v1%z + v2%z
   
  end function
 
  elemental function usereal_add(v1,v2) result (v3)
   use sys
  
   implicit none
  
   type(use), intent(in) :: v1
   real(kind=db), intent(in) :: v2
   real(kind=db)  :: v3
   
   v3 = v1%z + v2
   
  end function
  
 
  elemental function realuse_add(v1,v2) result (v3)
   use sys
  
   implicit none
  
   type(use), intent(in) :: v2
   real(kind=db), intent(in) :: v1
   real(kind=db)  :: v3
   
   v3 = v2%z + v1
   
  end function
  
  
  
  
  
  
  
  
 
  elemental function usereal_sub(v1,v2) result (v3)
   use sys
  
   implicit none
  
   type(use), intent(in) :: v1
   real(kind=db), intent(in) :: v2
   real(kind=db)  :: v3
   
   v3 = v1%z - v2
   
   
  end function
 
  elemental function use_sub(v1,v2) result (v3)
   use sys
  
   implicit none
  
   type(use), intent(in) :: v1,v2
   real(kind=db)  :: v3
   
   v3 = v1%z - v2%z
   
   
  end function
 
  elemental function realuse_sub(v1,v2) result (v3)
   use sys
  
   implicit none
  
   type(use), intent(in) :: v2
   real(kind=db), intent(in) :: v1
   real(kind=db)  :: v3
   
   v3 = v1 - v2%z
   
   
  end function
  
  
  
  
  
  
  
 
  elemental function usereal_dot(v1,v2) result (v3)
   use sys
  
   implicit none
  
   type(use), intent(in) :: v1
   real(kind=db), intent(in) :: v2
   real(kind=db)  :: v3
   
   v3 = v1%z * v2
   
  end function
 
  elemental function use_dot(v1,v2) result (v3)
   use sys
  
   implicit none
  
   type(use), intent(in) :: v1,v2
   real(kind=db)  :: v3
   
   v3 = v1%z * v2%z
   
  end function
 
  elemental function realuse_dot(v1,v2) result (v3)
   use sys
  
   implicit none
  
   type(use), intent(in) :: v2
   real(kind=db), intent(in) :: v1
   real(kind=db)  :: v3
   
   v3 = v2%z * v1
   
  end function
  
  
  
  
  
  
  
  
  
  
 
  elemental function usereal_div(v1,v2) result (v3)
   use sys
  
   implicit none
  
   type(use), intent(in) :: v1
   real(kind=db), intent(in) :: v2
   real(kind=db)  :: v3
   
   v3 = v1%z / v2
   
   
  end function
 
  elemental function use_div(v1,v2) result (v3)
   use sys
  
   implicit none
  
   type(use), intent(in) :: v1,v2
   real(kind=db)  :: v3
   
   v3 = v1%z / v2%z
   
   
  end function
 
  elemental function realuse_div(v1,v2) result (v3)
   use sys
  
   implicit none
  
   type(use), intent(in) :: v2
   real(kind=db), intent(in) :: v1
   real(kind=db)  :: v3
   
   v3 = v1 / v2%z
   
   
  end function
 
 
 
end module 
