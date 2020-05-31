program main

use test, only: test_value
use io, only: stdout

use interpolation

implicit none

integer, parameter :: n = 5
double precision :: x(n), y(n), yp1, ypn, y2(n), xi, yi, dyi, d2yi

x(1) =  -3.5215170604995243d0
x(2) =  0.86276423386836054d0
x(3) =   3.7201642425940626d0
x(4) =   4.8004321268773040d0
x(5) =   5.0957172699573832d0
y(1) =   3.8992594179963991d0
y(2) =   4.8550456214252433d0
y(3) =   5.0722670175588256d0
y(4) =   2.8530417245485067d0
y(5) =   1.5917686109483873d0
yp1 =  0.1d0
ypn = -1.7d0
call spline(x,y,n,yp1,ypn,y2)
call test_value(y2(1), 6.4808193833249039d-2,'spline: constrained-slope: y2(1)')
call test_value(y2(2), 3.1873621935649632d-2,'spline: constrained-slope: y2(2)')
call test_value(y2(3),-0.55913359001896401d0,'spline: constrained-slope: y2(3)')
call test_value(y2(4), -7.8404689972381414d0,'spline: constrained-slope: y2(4)')
call test_value(y2(5),  30.044543365446572d0,'spline: constrained-slope: y2(5)')
yp1 =  1.1d30
ypn = -1.1d30
call spline(x,y,n,yp1,ypn,y2)
call test_value(y2(1),                0.0d0,'spline: unconstrained-slope: y2(1)')
call test_value(y2(2),0.15288513342863624d0,'spline: unconstrained-slope: y2(2)')
call test_value(y2(3),-1.0730680213363821d0,'spline: unconstrained-slope: y2(3)')
call test_value(y2(4),-4.4138856668627522d0,'spline: unconstrained-slope: y2(4)')
call test_value(y2(5),                0.0d0,'spline: unconstrained-slope: y2(5)')

xi = 3.72174d0
call spline_interpolation(x,y,y2,n,xi,yi,dyi,d2yi)
write(0,*) yi, dyi, d2yi
call test_value(yi,    5.0708896817184534d0,'spline_interpolation: yi')
call test_value(dyi, -0.87492661536198912d0,'spline_interpolation: dyi')
call test_value(d2yi, -1.0779411812473771d0,'spline_interpolation: d2yi')

write(stdout,*) 'interpolation_module unit test passed'
end
