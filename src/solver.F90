

module solver
 
 
 implicit none
 
 
 
 contains
 
 subroutine set_solver()
  use solver_variables
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
  call create_field(y,.false.)
  call create_field(inty,.false.)
  call create_field(thavx,.false.)
  call create_field(thavy,.false.)
  call create_field(inthavx,.false.)
  call create_field(inthavy,.false.)
  call create_field(ap,.false.)
  call create_field(z,.false.)
  
  
  pres=0.0d0
  call start_sync(pres)
  r=0.0d0
  call start_sync(r)
  p=0.0d0
  call start_sync(p)
  
  
  
  
  call end_sync(s)
  thavx=-Ax(s)
  thavy=-Ay(s)
       
   
  inthavx=merge(0.0d0,1.0d0/thavx%bz,(thavx%bz == 0.0d0))
  inthavy=merge(0.0d0,1.0d0/thavy%bz,(thavy%bz == 0.0d0))
  
  
 end subroutine
 
  
 
 
 
 subroutine surfpressure(ierr)
  use solver_variables
  use parallel
  use grid_operate
  use sync
 
 
  implicit none

  integer :: n
  integer, intent(out) :: ierr
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
  
  ap=Gx(thavx%bz*Gx(pres))+   &
            Gy(thavy%bz*Gy(pres))
            
  
 
 
  
  r =  inty
  
  call start_sync(r)
  
   
  rdotfrac = 1.0d-16
   
  
  
  call end_sync(r)
  z=Ax(inthavx%bz*Ax(r))+Ay(inthavy%bz*Ay(r))
  
   
 
    
  rdotnew_tmp = sum(z%bz*r%bz)
  call real_allsum(rdotnew_tmp, rdotold)
  rdotfirst = rdotfrac*rdotold

  !if (proc_name == proc_master) print *, rdotold


  r =  inty%bz - ap%bz

  call start_sync(r)



  call end_sync(r)

  z=Ax(inthavx%bz*Ax(r))+Ay(inthavy%bz*Ay(r))


  p = z%bz

  call start_sync(p)

  rdotnew_tmp = sum(z%bz*r%bz)
  call real_allsum(rdotnew_tmp, rdotold)


 

 ! if (proc_name == proc_master) print *, rdotold
 
 
  do n = 1, 1000
  
     
   if (abs(rdotold) /= 0.0) then
   
      
   call end_sync(p)
   ap=Gx(thavx%bz*Gx(p))+Gy(thavy%bz*Gy(p))
  
   
   
   pdotap_tmp=sum(p%bz*ap%bz)
   call real_allsum(pdotap_tmp, pdotap)
    
   
   alpha = rdotold/pdotap
   
     
   
     
   r=r%bz-alpha*ap%bz
   
   
   
   call start_sync(r)
   
   
  
   pres=pres%bz+alpha*p%bz
   
   
   
   
   call end_sync(r)
   
   z=Ax(inthavx%bz*Ax(r))+Ay(inthavy%bz*Ay(r))
   
   
      
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
  
  
   
    call start_sync(pres)
   
  return
 
 end subroutine
 
 subroutine prescorrection(n)
  use solver_variables
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
    
end module
 

