program main

use test, only: test_value
use io, only: stdout

use earth

implicit none
double precision :: xyz_pole(3), geo_pole(3), x, y, z, lon, lat, h
integer :: ierr

! Calculate the pole of ITRF08 reference frame with respect to MORVEL56 no-net-rotation frame
! call itrf2008_wrt_morvel56_PA()
! call pole_xyz2geo(-0.190d0,-0.442d0,0.915d0,lon,lat,mag)
! write(0,*) lon,lat,mag*masa2degma
! call pole_xyz2geo(-0.252d0,-0.302d0,0.643d0,lon,lat,mag)
! write(0,*) lon,lat,mag*masa2degma

! pole_xyz2geo
xyz_pole(1) = -0.190d0
xyz_pole(2) = -0.442d0
xyz_pole(3) =  0.915d0
call pole_xyz2geo(xyz_pole(1),xyz_pole(2),xyz_pole(3),geo_pole(1),geo_pole(2),geo_pole(3),'sphere')
call test_value(geo_pole(1),-113.26114046643323d0,'pole_xyz2geo: geo_pole(1)')
call test_value(geo_pole(2), 62.264603563491839d0,'pole_xyz2geo: geo_pole(2)')
call test_value(geo_pole(3), 1.0337741532849427d0,'pole_xyz2geo: geo_pole(3)')
call test_value(geo_pole(3)*masa2degma,0.28715948702359523d0,'masa2degma')

! pole_geo2xyz
call pole_geo2xyz(geo_pole(1),geo_pole(2),geo_pole(3),xyz_pole(1),xyz_pole(2),xyz_pole(3),'sphere')
call test_value(xyz_pole(1),-0.190d0,'pole_geo2xyz: xyz_pole(1)')
call test_value(xyz_pole(2),-0.442d0,'pole_geo2xyz: xyz_pole(2)')
call test_value(xyz_pole(3), 0.915d0,'pole_geo2xyz: xyz_pole(3)')

! get_pole
call get_pole('NA','EU','MORVEL',geo_pole,ierr)
call test_value(geo_pole(1), 139.46111978128923d0,'get_pole: MORVEL NA/EU lon')
call test_value(geo_pole(2), 61.796375604003913d0,'get_pole: MORVEL NA/EU lat')
call test_value(geo_pole(3),0.21067757026307590d0,'get_pole: MORVEL NA/EU vel (deg/Ma)')

! sphere_geo2xyz
! sphere_xyz2geo

! ellipsoid_geo2xyz
lon = 1.4223653911d+02
lat = 4.3889293271d+01
h = 1.7767427755d+02
call ellipsoid_geo2xyz(lon,lat,h,x,y,z,semimajor_grs80,eccentricity_grs80)
call test_value(x,-3.6397836984d+06,'ellipsoid_geo2xyz: x')
call test_value(y,2.8195896747d+06,'ellipsoid_geo2xyz: y')
call test_value(z,4.3993581171d+06,'ellipsoid_geo2xyz: z')
! -3.6397836984E+06  2.8195896747E+06  4.3993581171E+06

! ellipsoid_xyz2geo
x = -3.6397836984d+06
y =  2.8195896747d+06
z =  4.3993581171d+06
call ellipsoid_xyz2geo(x,y,z,lon,lat,h,semimajor_grs80,eccentricity_grs80,10)
call test_value(lon,1.4223653911d+02,'ellipsoid_xyz2geo: lon')
call test_value(lat,4.3889293271d+01,'ellipsoid_xyz2geo: lat')
call test_value(h,1.7767427755d+02,'ellipsoid_xyz2geo: height')

write(stdout,*) 'earth_module unit test passed'
end
