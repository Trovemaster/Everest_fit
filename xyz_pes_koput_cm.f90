module caoh_param

use accuracy


public :: caoh_init_pot,caoh_pot,caoh_potmin,poten_xyz

private
real(ark) :: deg, aucm, ang
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
caoh_pot =  0 
!caoh_pot=aucm*poten_xyz(xl)
end function caoh_pot

!
function poten_xyz(Nterms,r1e,r2e,thetae,ipower1,ipower2,ipower3,force,r1,r2,alpha) result (f)
  implicit none
  integer(ik),intent(in) :: Nterms,ipower1(Nterms),ipower2(Nterms),ipower3(Nterms)
  real(rk),intent(in) :: r1e,r2e,thetae,force(Nterms),r1,r2,alpha

  real(rk) :: y1, y2, y3, f
  real(rk) :: vpot
  integer(ik) :: p1,p2,p3,ipar

  deg = pi/180.0_rk
  y1=(r1-r1e)/r1 ! streching coordinate 1
  y2=(r2-r2e)/r2 ! stretching coordinate 2
  y3=(alpha-thetae) ! bending coordinate
  !
  vpot=0
  !
  ! construct potential function
  do ipar=1,Nterms
    p1=ipower1(ipar)
    p2=ipower2(ipar)
    p3=ipower3(ipar)
    vpot=vpot+force(ipar)*y1**p1*y2**p2*y3**p3
  enddo

  f=vpot

end function poten_xyz




end module

double precision function pot(x)
use caoh_param
double precision x(3)
pot=caoh_pot(x)
end

subroutine init_pot(name,initok,icm)
use caoh_param
integer initok,icm
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
