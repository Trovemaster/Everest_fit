module caoh_param

integer, parameter :: ark         = selected_real_kind(25,32)
real(ark) :: deg, pi, aucm, ang
integer ipar, ieq, parmax, info, i1(40), i2(40), i3(40)
double precision  par(3), par3min
double precision  force(40)
character(9) buf30
integer parvar, iparvar(40)

public :: caoh_init_pot,caoh_pot,caoh_potmin,caoh_par_der

contains

subroutine caoh_potmin(r1,r2,x)
double precision r1,r2,x
r1=par(2)/ang
r2=par(1)/ang
x=par(3)
end subroutine caoh_potmin

subroutine caoh_init_pot(name,initok)
integer initok
character*(*) name
character(500) string
initok=0
open(11,file=name)

  pi = 4.0d0 * datan2(1.0d0,1.0d0)

  deg = pi/180.0d0

  aucm=1d0/219474.6313702d0

  ang=0.52917721067d0

  !read number of potential parameters
  read (11,*,err=800,end=800)  parmax

  parvar=0

  !read equilibrium geometry parameters
  do ieq =1, 3
    read (11,*,err=800,end=800) buf30, par(ieq)
  enddo

!! use radians
  par(3)=par(3)*deg

  !read potential parameters
  do ipar=1, parmax
!    read (11,*,err=800,end=800) buf30, i1(ipar), i2(ipar), i3(ipar), force(ipar)
    read (11,'(a)',err=800,end=800) string
    buf30=string(1:9)
    read (string(10:),*,err=800,end=800) i1(ipar), i2(ipar), i3(ipar), force(ipar)
    if(buf30(8:8).ne.' ') then
     parvar=parvar+1
     iparvar(parvar)=ipar
    endif
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

  real(ark) :: y1, y2, y3, r1e, r2e, alphae, deg, pi
  real(ark) :: r1,r2,alpha,vpot

!!!  real(ark), parameter :: tocm = 219474.63067_ark

!!!  pi = 4.0d0 * datan2(1.0d0,1.0d0)

!!!  deg = pi/180.0d0

!!  r1e    = par(1)
!!  r2e    = par(2)
!!  alphae = par(3)

!!! need angstrom
  r1  = pin(1)*ang
  r2  = pin(2)*ang
!!! cos to angle  alpha  = acos(local(3))
  alpha  = max(acos(pin(3)),par3min)

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

subroutine caoh_par_der(pin,der)
  implicit none

  integer ipar,ivar
  double precision pin(3),der(*)

  integer :: p1, p2, p3

  real(ark) :: y1, y2, y3, r1,r2,alpha

  r1  = pin(2)*ang
  r2  = pin(1)*ang
  alpha  = max(acos(pin(3)),par3min)

  y1=(r1-par(1))/r1 ! streching coordinate 1
  y2=(r2-par(2))/r2 ! stretching coordinate 2
  y3=(alpha-par(3)) ! bending coordinate

  do ivar=1,parvar
    ipar=iparvar(ivar)
    p1=i1(ipar)
    p2=i2(ipar)
    p3=i3(ipar)
    der(ivar)=y1**p1*y2**p2*y3**p3
  enddo

end subroutine caoh_par_der

subroutine caoh_par_der_lab(ilab)
implicit none
integer ipar,ivar
integer ilab(3,*)

do ivar=1,parvar
 ipar=iparvar(ivar)
! ilab(1,ivar)=i1(ipar)
! ilab(2,ivar)=i2(ipar)
! ilab(3,ivar)=i3(ipar)
 ilab(1,ivar)=ipar
 ilab(2,ivar)=0
 ilab(3,ivar)=0
enddo

end subroutine caoh_par_der_lab

end module

!!!!!!!!=======

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

subroutine pot_nvar(nvar)
use caoh_param
integer nvar
nvar=parvar
end subroutine pot_nvar

subroutine pot_der(x,der)
use caoh_param
double precision x(3),der(*)
!write(*,*) 'got to pot der ',x
call caoh_par_der(x,der)
end subroutine pot_der

subroutine pot_der_lab(ilab)
use caoh_param
integer ilab(3,*)
call caoh_par_der_lab(ilab)
end subroutine pot_der_lab
