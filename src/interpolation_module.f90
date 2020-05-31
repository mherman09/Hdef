module interpolation

public :: spline
public :: spline_interpolation

!==================================================================================================!
contains
!==================================================================================================!

subroutine spline(x,y,n,yp1,ypn,y2)
!----
! (From Press et al., Numerical Recipes)
!
! Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e., yi = f(xi),
! with x1 < x2 < ... < xN, and given values yp1 and ypn for the first derivative of the
! interpolating function at points 1 and n, respectively, this routine returns an
! array y2(1:n) of length n which contains the second derivatives of the interpolating
! function at the tabulated points xi. If yp1 and/or ypn are equal to 1 × 10^30 or larger,
! the routine is signaled to set the corresponding boundary condition for a natural spline,
! with zero second derivative on that boundary.
!
! Parameter: NMAX is the largest anticipated value of n.
!----

implicit none

! Arguments
integer :: n
double precision :: yp1, ypn, x(n), y(n), y2(n)

! Local variables
integer, parameter :: nmax = 500
integer :: i, k
double precision :: p, qn, sig, un, u(nmax)

! Set the lower boundary condition to be natural or specified first derivative
if (abs(yp1).ge.1.0d30) then
    y2(1) = 0.0d0
    u(1) = 0.0d0
else
    y2(1) = -0.5d0
    u(1) = (3.0d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
endif

! Tridiagonal algorithm
do i = 2,n-1
    sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
    p = sig*y2(i-1)+2.0d0
    y2(i) = (sig-1.0d0)/p
    u(i) = (6.0d0*( &
               (y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)) &
           )/(x(i+1)-x(i-1))-sig*u(i-1))/p
enddo

! Set the upper boundary condition to be natural or specified first derivative
if (abs(ypn).gt.1.0d30) then
    qn = 0.0d0
    un = 0.0d0
else
    qn = 0.5d0
    un = (3.0d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
endif

y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.0d0)

! Back-substitution of the tridiagonal algorithm
do k = n-1,1,-1
    y2(k) = y2(k)*y2(k+1)+u(k)
enddo

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine spline_interpolation(xa,ya,y2a,n,x,y,dy,d2y)
!----
! (From Press et al., Numerical Recipes)
!
! Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a function (with the xai’s in
! order), and given the array y2a(1:n), which is the output from spline above, and given a value
! of x, this routine returns a cubic-spline interpolated value y.
!----

implicit none

! Arguments
integer :: n
double precision :: x, y, dy, d2y, xa(n), y2a(n), ya(n)

! Local variables
integer :: k, khi, klo
double precision :: a, b, h

klo = 1
khi = n

! We will find the right place in the table by means of bisection. This is optimal if sequential
! calls to this routine are at random values of x. If sequential calls are in order, and closely
! spaced, one would do better to store previous values of klo and khi and test if they remain
! appropriate on the next call.

! khi and klo bracket the input value of x
do while (khi-klo.gt.1)
    k = (khi+klo)/2
    if (xa(k).gt.x) then
        khi = k
    else
        klo = k
    endif
enddo

h = xa(khi)-xa(klo)

! The xa's must be distinct
if (abs(h).lt.1.0d-10) then
    write(0,*) 'spline_interpolation: bad xa input'
endif

! Evaluate cubic spline polynomial value, derivative, and second derivative
a = (xa(khi)-x)/h
b = (x-xa(klo))/h
y = a*ya(klo) + b*ya(khi) + &
        ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0d0
dy = (ya(khi)-ya(klo))/(xa(khi)-xa(klo)) - &
     (3.0d0*a**2-1.0d0)*(xa(khi)-xa(klo))*y2a(klo)/6.0d0 + &
     (3.0d0*b**2-1.0d0)*(xa(khi)-xa(klo))*y2a(khi)/6.0d0
d2y = a*y2a(klo) + b*y2a(khi)

return
end subroutine

end module
