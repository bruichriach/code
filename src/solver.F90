

#include "include.h"

#ifdef ALLOW_RIGID_LID

module solver

 
 
 implicit none
 
 
 
 contains
 
 subroutine set_solver()
  use variables
  use params
  use grid_operate
  use grid
  use sync
  use params
  use allocation
  
  implicit none
 
  call create_field(pres,'pres',.true.)
  call create_field(p,.true.)
  call create_field(r,.true.)
  call create_field(y,.false.,1)
  call create_field(inty,.false.)
  call create_field(thavx,.false.)
  call create_field(thavy,.false.)
  call create_field(inthavx,.false.)
  call create_field(inthavy,.false.)
  call create_field(ap,.false.)
  call create_field(z,.false.)
  
   
  
  
  
  call end_sync(s)
  thavx=-Ax(s)
  thavy=-Ay(s)
       
   
  inthavx=merge(0.0d0,1.0d0/thavx%bz,(thavx%bz == 0.0d0))
  inthavy=merge(0.0d0,1.0d0/thavy%bz,(thavy%bz == 0.0d0))
  
  
 end subroutine
 
  
 
 
 
 subroutine surfpressure(ierr,maxit)
  use variables
  use parallel
  use grid_operate
  use sync
 
 
  implicit none

  integer :: n
  integer, intent(out) :: ierr
  integer, intent(in) :: maxit
  real (kind=8) :: alpha, rdotold, rdotnew, rdotfirst, rdotfrac
  real (kind=8) :: pdotap, p_0
  real (kind=8) :: pdotap_tmp, rdotnew_tmp
   
  
   
  ierr = 0
  p_0=0
 ! p_0_tmp=sum(pres%z%z(1:pres%nx,1:pres%ny))
 ! call real_allsum(p_0_tmp, p_0)
 ! p_0=p_0/real(mx*my)
  
  
  
 ! pres%z%z%z%z=0.0d0!pres%z%z%z%z-p_0
 ! r%z%z%z%z=0.0d0
 ! p%z%z%z%z=0.0d0
 
 
  call end_sync(pres)
  ap=equation(pres)
            
  
 
 
  
  r =  inty
  
  call start_sync(r)
  
   
  rdotfrac = 1.0d-12
   
  
  
  call end_sync(r)
  z=precondition(r)
  
   
 
    
  rdotnew_tmp = sum(z%bz*r%bz)
  call real_allsum(rdotnew_tmp, rdotold)
  rdotfirst = rdotfrac*rdotold

  !if (proc_name == proc_master) print *, rdotold


  r =  inty%bz - ap%bz

  call start_sync(r)



  call end_sync(r)
  z=precondition(r)


  p = z%bz

  call start_sync(p)

  rdotnew_tmp = sum(z%bz*r%bz)
  call real_allsum(rdotnew_tmp, rdotold)


!  if (proc_name == proc_master) print *, rdotold, rdotfirst
 
 
  do n = 1, maxit
  
     
   if (abs(rdotold) /= 0.0) then
   
      
   call end_sync(p)
   ap=equation(p)
  
   
   
   pdotap_tmp=sum(p%bz*ap%bz)
   call real_allsum(pdotap_tmp, pdotap)
    
   
   alpha = rdotold/pdotap
   
     
   
     
   r=r%bz-alpha*ap%bz
   
   
   
   call start_sync(r)
   
   
  
   pres=pres%bz+alpha*p%bz
   
   
   call end_sync(r)
   z=precondition(r)
   
   
      
   rdotnew_tmp = sum(z%bz*r%bz)
   call real_allsum(rdotnew_tmp, rdotnew)
        
        
  
   p=z%bz+(rdotnew/rdotold)*p%bz
   
   call start_sync(p)
    
   
   
   rdotold=rdotnew
!   if (proc_name == proc_master) print *, n, rdotold, pdotap


  

   if (abs(rdotold) <=  rdotfirst) then
  
 
   !if (proc_name == ens_master) print *, ens_name, n, rdotold
    call start_sync(pres)
    
    
    return
   end if
   
   else
   
    if (proc_name == ens_master) then
     print *, 'exact solve'
     print *, ens_name
    end if
   
    call start_sync(pres)

    
    return
   
   end if
   
  end do
   
  if (proc_name == ens_master) then
   print *, 'no convergence'
   print *, rdotold, ens_name
  end if
  ierr=1
   
   print *, proc_name, alpha
  
  
   
    call start_sync(pres)
   
  return
 
 end subroutine
 
 
 function equation(a) result(b)
  use grid_operate
  use variables
   
  implicit none
   
  type(hvar), intent(inout) :: a
  real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx,a%p%ly+1:a%p%ly+a%p%ny)
 
 
  b=Gx(thavx%bz*Gx(a))+   &
            Gy(thavy%bz*Gy(a))
 
 end function
 
 
 function precondition(a) result(b)
  use variables
  use grid_operate
   
  implicit none
   
  type(hvar), intent(inout) :: a
  real(kind=db) :: b(a%p%lx+1:a%p%lx+a%p%nx,a%p%ly+1:a%p%ly+a%p%ny)
 
 
  b=Ax(inthavx%bz*Ax(a))+Ay(inthavy%bz*Ay(a))
 
 end function
 
 subroutine prescorrection(n)
  use global
  use grid_operate
  use sync
  use params
  use variables
  use overload
  
 
  implicit none
 
  integer :: i
  integer, intent(in) :: n
  real(kind=db), pointer :: tendu(:,:), tendv(:,:)
  
  
   do i=1,nz
    if (n == 1) then
     tendu => u(i)%tend2%bz
     tendv => v(i)%tend2%bz
    else
     tendu => u(i)%tend1%bz
     tendv => v(i)%tend1%bz
    end if
    call end_sync(pres)
    tendu=tendu + Gx(pres)
    tendv=tendv + Gy(pres)
   end do
  
  
  
  return
 
 end subroutine
 
 
 

 subroutine thickness_correct(h,s)
  use global
  use params
  use grid
  use overload
  use operations
 
  type(hvar), intent(inout) :: h(nz)
  type(hvar), intent(inout) :: s
  real(kind=db) :: d(s%p%lx+1:s%p%lx+s%p%nx,s%p%ly+1:s%p%ly+s%p%ny)
  integer :: k
  
  d=s%bz
  do k=1,nz-1
   d=d+h(k)%bz
  end do
  h(nz)=-d
  
 end subroutine
    
end module
 
#endif

