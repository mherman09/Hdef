program main

use test, only: test_value
use io, only: stdout

use random

implicit none

integer :: i
double precision :: dp1, dp2
real :: real1, real2
complex :: c1, c2
complex (kind=8) :: dc1, dc2

! timeseed()
iseed = timeseed()
if (iseed.gt.0.and.iseed.lt.240000) then
    iseed = 1
else
    iseed = 0
endif
call test_value(iseed,1,'timeseed()')

! ran0()
iseed = -123456
dp1 = ran0(iseed)
call test_value(dp1,0.81047654151916504d0,'ran0(): value')
call test_value(iseed,123457,'ran0(): iseed')

! c4_normal_01
iseed = 98498171
c1 = c4_normal_01(iseed)
call test_value(c1,(-0.474812567,-0.152586848),'c4_normal_01()')

! c8_normal_01
iseed = 552792940
dc1 = c8_normal_01(iseed)
call test_value(dc1,(-0.98438105957821054d0,-1.0320012536812111d0),'c8_normal_01()')

! i4_normal_ab
iseed = 95061600
real1 = 10.0 ! mean
real2 = 6.0  ! standard deviation
i = i4_normal_ab(real1,real2,iseed)
call test_value(i,9,'i4_normal_ab(): value')
call test_value(iseed,1291390176,'i4_normal_ab(): iseed')

! i8_normal_ab

! r4_normal_01
iseed = 726391
real1 = r4_normal_01(iseed)
call test_value(real1,0.577163041,'r4_normal_01(): value')
call test_value(iseed,1858576450,'r4_normal_01(): iseed')

! r4_normal_ab
iseed = 9347
real1 = 4.0 ! mean
real2 = 1.0  ! standard deviation
real1 = r4_normal_ab(real1,real2,iseed)
call test_value(real1,1.72496796,'r4_normal_ab(): value')
call test_value(iseed,1038750240,'r4_normal_ab(): iseed')

! r4_uniform_01
iseed = 2384323
real1 = r4_uniform_01(iseed)
call test_value(real1,0.660592258,'r4_uniform_01(): value')
call test_value(iseed,1418611015,'r4_uniform_01(): iseed')

! r4vec_uniform_01
! r4vec_normal_ab

! r8_normal_01
iseed = 23823
dp1 = r8_normal_01(iseed)
call test_value(dp1,-1.2982174417976360d0,'r8_normal_01(): value')
call test_value(iseed,1341590876,'r8_normal_01(): iseed')

! r8_normal_ab
iseed = 92762192
dp1 = -45.0d0
dp2 = 7.0d0
dp1 = r8_normal_ab(dp1,dp2,iseed)
call test_value(dp1,-45.869090907966878d0,'r8_normal_ab(): value')
call test_value(iseed,1200425557,'r8_normal_ab(): iseed')

! r8_uniform_01
iseed = 142
dp1 = r8_uniform_01(iseed)
call test_value(dp1,1.1113444347797751d-3,'r8_uniform_01(): value')
call test_value(iseed,2386594,'r8_uniform_01(): iseed')


! r8mat_normal_01
! r8mat_normal_ab
! r8vec_normal_01
! r8vec_normal_ab
! r8vec_uniform_01

write(stdout,*) 'random_module unit test passed'

end
