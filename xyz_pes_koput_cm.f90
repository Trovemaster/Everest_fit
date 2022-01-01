module caoh_param


public :: caoh_init_pot,caoh_pot,caoh_potmin,poten_func

private
integer, parameter :: ark         = selected_real_kind(25,32)
real(ark) :: deg, pi, aucm, ang
integer ieq, parmax, i1(40), i2(40), i3(40)
double precision  par(3), par3min
double precision  force(40)
character(9) buf30

contains

subroutine caoh_potmin(r1,r2,x)
double precision r1,r2,x
r2=par(1)/ang
r1=par(2)/ang
x=par(3)
end subroutine caoh_potmin

subroutine caoh_init_pot(name,initok)

implicit none
integer initok
character*(*) name
integer:: ipar

initok=0
open(11,file=name)

  pi = 4.0d0 * datan2(1.0d0,1.0d0)

  deg = pi/180.0d0

  aucm=1d0/219474.6313702d0

  ang=0.52917721067d0

  !read number of potential parameters
  read (11,*,err=800,end=800)  parmax

  !read equilibrium geometry parameters
  do ieq =1, 3
    read (11,*,err=800,end=800) buf30, par(ieq)
  enddo

!! use radians
  par(3)=par(3)*deg

  !read potential parameters
  do ipar=1, parmax
    read (11,*,err=800,end=800) buf30, i1(ipar), i2(ipar), i3(ipar), force(ipar)
  enddo
  read(11,*) par3min
  par3min=par3min*deg
close(11)
initok=1
800 return
end subroutine caoh_init_pot

double precision function caoh_pot(x)
double precision x(3),xl(3)
xl(1)=x(2)
xl(2)=x(1)
xl(3)=x(3)
caoh_pot=aucm*poten_xyz(xl)
end function caoh_pot

!!function poten_xyz(local,parmax,par,i1,i2,i3,force) result (f)
function poten_xyz(pin) result (f)
  implicit none

!!  integer, parameter :: ark         = selected_real_kind(25,32)
!!  integer, parameter :: ik          = selected_int_kind(8)

!!  integer parmax, ipar, i1(40), i2(40), i3(40)
  integer ipar
!!  double precision local(3), force(39), par(3)
  double precision pin(3)
  double precision f

  integer :: p1, p2, p3

  real(ark) :: y1, y2, y3
  real(ark) :: r1,r2,alpha,vpot

!!!  real(ark), parameter :: tocm = 219474.63067_ark

!!!  pi = 4.0d0 * datan2(1.0d0,1.0d0)

!!!  deg = pi/180.0d0

!!  r1e    = par(1)
!!  r2e    = par(2)
!!  alphae = par(3)

!!! need angstrom
  r1  = pin(1) !*ang
  r2  = pin(2) !*ang
!!! cos to angle  alpha  = acos(local(3))
  !alpha  = max(acos(pin(3)),par3min)
  
  alpha  = pin(3)

!!  y1=(r1-r1e)/r1 ! streching coordinate 1
  y1=(r1-par(1))/r1 ! streching coordinate 1
!!  y2=(r2-r2e)/r2 ! stretching coordinate 2
  y2=(r2-par(2))/r2 ! stretching coordinate 2
!!! use radians y3=(alpha-alphae)*deg ! bending coordinate
!!  y3=(alpha-alphae) ! bending coordinate
  y3=(alpha-par(3)) ! bending coordinate

  vpot=0.0_ark

  ! construct potential function
  do ipar=1,parmax
    p1=i1(ipar)
    p2=i2(ipar)
    p3=i3(ipar)
    vpot=vpot+force(ipar)*y1**p1*y2**p2*y3**p3
  enddo

  f=vpot

end function poten_xyz



 !
 ! Defining potential energy function 
 !
 subroutine poten_func(force,f,ifield,r12,r32,theta)
   !
   implicit none
   !
   integer,parameter           :: N = 99
   !
   double precision,intent(in) ::  force(N),r12,r32,theta
   double precision,intent(out)::  f
   integer,intent(in)          :: ifield
   double precision            :: pin(3),f0
   integer                     :: n1,initok
   !
   character(len = 80)   :: name
   !
   N1 = size(force)
   !
   select case (ifield)
     !
   case (1)
     !
     name = 'field1.fit'
     !
   case(2)
     !
     name = 'field2.fit'
     !
   case(3)
     !
     name = 'field3.fit'
     !
   case(4)
     !
     name = 'field4.fit'
     !
   end select 
   
   call caoh_init_pot(name,initok)
   !
   pin(1) = par(1)
   pin(2) = par(2)
   pin(3) = par(3)
   !
   f0 = poten_xyz(pin)
   !
   pin(1) = r12
   pin(2) = r32
   pin(3) = theta
   !
   f = poten_xyz(pin)-f0
   !
  end subroutine poten_func




end module

double precision function pot(x)
use caoh_param
double precision x(3)
pot=caoh_pot(x)
end

subroutine init_pot(name,initok,icm)
use caoh_param
integer nrco,initok
character*(*) name
call caoh_init_pot(name,initok)
icm=1
end subroutine init_pot

subroutine potmin(r1,r2,x)
use caoh_param
double precision r1,r2,x
call caoh_potmin(r1,r2,x)
end subroutine potmin

subroutine close_pot
end subroutine close_pot
