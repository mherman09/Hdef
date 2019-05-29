program main

use test, only: test_value
use io, only: stdout

use geom

implicit none

double precision :: lon1, lat1, lon2, lat2, dist, az

lon1 = 5.0d0
lat1 = 52.0d0
dist = 123.0d0/6371.0d0
az = -50.0d0
call distaz2lola(lon1,lat1,dist,az,lon2,lat2)
call test_value(lon2,3.6015261937179939d0,'distaz2lola: lon2')
call test_value(lat2,52.702854017526889d0,'distaz2lola: lat2')

lon1 = -77.9d0
lat1 = 40.8d0
lon2 = -90.3d0
lat2 = 38.6d0
call lola2distaz(lon1,lat1,lon2,lat2,dist,az)
call test_value(dist,0.17072370996955819d0,'lola2distaz: dist')
call test_value(az,-98.965076442115560d0,'lola2distaz: az')

! str2vec
! strdip2normal
! strdip2updip
! normal2strike
! normal2updip
! pnpoly


write(stdout,*) 'geom_module unit test passed'
end
