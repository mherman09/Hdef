module earth

double precision :: radius_earth_m = 6371.0d3
double precision :: radius_earth_km = 6371.0d0

public :: radius_earth_m, radius_earth_km

!--------------------------------------------------------------------------------------------------!
contains
!--------------------------------------------------------------------------------------------------!

subroutine distaz(dist,az,lon1,lat1,lon2,lat2)
!----
! Given the latitude and longitude of two points, compute the great circle distance between them
! and the azimuth from point 1 to point 2 measured clockwise from north. Input coordinates are in
! degrees, output distance and azimuth are in radians.
!----

use trig, only: pi, d2r
implicit none

! Arguments
double precision :: dist, az, lon1, lat1, lon2, lat2

! Local variables
double precision :: lon1r, lat1r, lon2r, lat2r, colat1, colat2, dlon, dlat, a

! Check if points are polar opposite
if (dabs(lat1+lat2).le.1.0d-6.and.dabs(mod(lon1-lon2,1.8d2)).le.1.0d-6) then
    dist = pi
    az = 0.0d0
    return
endif

! Convert lon/lat to radians
lon1r = lon1*d2r
lat1r = lat1*d2r
lon2r = lon2*d2r
lat2r = lat2*d2r

! Haversine formula for distance (more accurate over short distances than Spherical Law of Cosines)
colat1 = pi/2.0d0-lat1r
colat2 = pi/2.0d0-lat2r
dlon = lon2r-lon1r
dlat = lat2r-lat1r
a = dsin(dlat/2.0d0)*dsin(dlat/2.0d0) + dcos(lat1r)*dcos(lat2r)*dsin(dlon/2.0d0)*dsin(dlon/2.0d0)
if (a.ge.1.0d0) then
    dist = 0.0d0
else
    dist = 2.0d0*datan2(dsqrt(a),dsqrt(1.0d0-a))
endif

! Calculate azimuth
if (dist.lt.1.0d-6) then
    az = 0.0d0
else
    ! Spherical Law of Sines
    az = datan2(dsin(dlon)*dcos(lat2r), dcos(lat1r)*dsin(lat2r)-dsin(lat1r)*dcos(lat2r)*dcos(dlon))
endif

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine lola(lon2,lat2,lon1,lat1,dist,az)
!----
! Compute the longitude and latitude a distance and azimuth clockwise from north away from an
! initial point. Input coordinates and azimuth are in degrees, input distance is in km, and output
! coordinates are in radians.
!----

use trig, only: d2r
implicit none

! Arguments
double precision :: lon1, lat1, lon2, lat2, dist, az

! Local variables
double precision :: distr, azr, lon1r, lat1r, lat2r

! Convert all input quantities to radians
distr = dist/radius_earth_km ! km -> rad
azr = az*d2r          ! deg -> rad
lon1r = lon1*d2r      ! deg -> rad
lat1r = lat1*d2r      ! deg -> rad

! Spherical law of cosines
lat2 = dasin(dsin(lat1r)*dcos(distr)+dcos(lat1r)*dsin(distr)*dcos(azr))
lat2r = lat2*d2r

! Spherical law of cosines and law of sines
lon2 = lon1r + datan2(dsin(azr)*dsin(distr)*dcos(lat1r), dcos(distr)-dsin(lat1r)*dsin(lat2r))

return
end subroutine

end module
