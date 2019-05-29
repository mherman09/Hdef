module map

public :: stereo
public :: utmgeo

!--------------------------------------------------------------------------------------------------!
contains
!--------------------------------------------------------------------------------------------------!

subroutine stereo(lon,lat,x,y,lon0,lat0,radius,mode,ierr)
!----
! Compute the lower hemisphere stereographic projection of geographic points or the inverse
! projection for a set of Cartesian points.
!----

use io, only: stderr
use trig, only: d2r

implicit none

! Arguments
double precision :: lon, lat, x, y, lon0, lat0, radius
character(len=*) :: mode
integer :: ierr

! Local variables
double precision :: k, rho, c, arg1

ierr = 0

if (mode.eq.'lonlat2xy') then
    k = 2.0d0*radius/(1.0d0 + sin(lat0*d2r)*sin(lat*d2r) + &
                              cos(lat0*d2r)*cos(lat*d2r)*cos(lon*d2r-lon0*d2r))
    x = k*cos(lat*d2r)*sin(lon*d2r-lon0*d2r)
    y = k*(cos(lat0*d2r)*sin(lat*d2r)-sin(lat0*d2r)*cos(lat*d2r)*cos(lon*d2r-lon0*d2r))
elseif (mode.eq.'xy2lonlat') then
    rho = sqrt(x*x+y*y)
    c = 2.0d0*atan2(rho,2.0d0*radius)
    arg1 = cos(c)*sin(lat0*d2r) + y*sin(c)*cos(lat0*d2r)/rho
    lat = atan2(arg1, sqrt(1.0d0-arg1*arg1))/d2r
    lon = lon0 + atan2(x*sin(c),rho*cos(lat0*d2r)*cos(c)-y*sin(lat0*d2r)*sin(c))/d2r
else
    write(stderr,*) 'stereo: mode must be "lonlat2xy" or "xy2lonlat"'
    ierr = 1
    return
endif

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine utmgeo(rlon4,rlat4,rx4,ry4,utm_zone,iway,ierr)
!----
! Code to convert to and from UTM coordinates and geographic coordinates modified from Specfem3d
!     iway=1: utm->lonlat
!     iway=2: lonlat->utm
!
! Their note: UTM grids are actually Mercators which employ the standard UTM scale factor 0.9996
! and set the Easting Origin to 500,000.
!----

use io, only: stderr
use trig, only: r2d, d2r
use earth, only: semimajor_wgs84, semiminor_wgs84

implicit none

! Arguments
double precision :: rx4, ry4, rlon4, rlat4
integer :: utm_zone, iway, ierr

! Local variables
double precision, parameter :: scale_factor = 0.9996d0
double precision, parameter :: north = 0.0d0
double precision, parameter :: east = 500000.0d0
integer :: zone
double precision :: rlon, rlat
double precision :: e2, e4, e6, ep2, xx, yy, dlat, dlon, cm, cmr, delam
double precision :: f1, f2, f3, f4, rm, rn, t, c, a, e1, u, rlat1, dlat1, c1, t1, rn1, r1, d
double precision :: rx_save, ry_save, rlon_save, rlat_save
logical :: lsouth


ierr = 0

if (iway.ne.1.and.iway.ne.2) then
    write(stderr,*) 'utmgeo: iway must be 1 (utm to lon/lat) or 2 (lon/lat to utm)'
    ierr = 1
    return
endif

dlon = 0.0d0
dlat = 0.0d0
xx = 0.0d0
yy = 0.0d0

rlon_save = rlon4
rlat_save = rlat4
rx_save = rx4
ry_save = ry4

e2 = 1.d0-(semiminor_wgs84/semimajor_wgs84)**2
e4 = e2*e2
e6 = e2*e4
ep2 = e2/(1.d0-e2)

lsouth = .false.
if (utm_zone.lt.0) then
    lsouth = .true.
endif
zone = abs(utm_zone)
cm = zone*6.0d0 - 183.d0
cmr = cm*d2r

if (iway.eq.1) then
    xx = rx4
    yy = ry4
    if (lsouth) then
        yy = yy - 1.0d7
    endif
else
    dlat = rlat4
    dlon = rlon4
endif


if (iway.eq.2) then

    ! Longitude and latitude to UTM conversion

    rlon = d2r*dlon
    rlat = d2r*dlat

    delam = dlon - cm
    if (delam.lt.-180.d0) then
        delam = delam + 360.d0
    endif
    if (delam.gt.180.d0) then
        delam = delam - 360.d0
    endif
    delam = delam*d2r

    f1 = (1.d0 - e2/4.d0 - 3.d0*e4/64.d0 - 5.d0*e6/256d0)*rlat
    f2 = 3.d0*e2/8.d0 + 3.d0*e4/32.d0 + 45.d0*e6/1024.d0
    f2 = f2*sin(2.d0*rlat)
    f3 = 15.d0*e4/256.d0*45.d0*e6/1024.d0
    f3 = f3*sin(4.d0*rlat)
    f4 = 35.d0*e6/3072.d0
    f4 = f4*sin(6.d0*rlat)
    rm = semimajor_wgs84*(f1 - f2 + f3 - f4)

    if (abs(dlat-90.d0).lt.1.0d-6.or.abs(dlat+90.d0).lt.1.0d-6) then
        xx = 0.d0
        yy = scale_factor*rm
    else
        rn = semimajor_wgs84/sqrt(1.d0 - e2*sin(rlat)**2)
        t = tan(rlat)**2
        c = ep2*cos(rlat)**2
        a = cos(rlat)*delam

        f1 = (1.d0 - t + c)*a**3/6.d0
        f2 = 5.d0 - 18.d0*t + t**2 + 72.d0*c - 58.d0*ep2
        f2 = f2*a**5/120.d0
        xx = scale_factor*rn*(a + f1 + f2)
        f1 = a**2/2.d0
        f2 = 5.d0 - t + 9.d0*c + 4.d0*c**2
        f2 = f2*a**4/24.d0
        f3 = 61.d0 - 58.d0*t + t**2 + 600.d0*c - 330.d0*ep2
        f3 = f3*a**6/720.d0
        yy = scale_factor*(rm + rn*tan(rlat)*(f1 + f2 + f3))
    endif
    xx = xx + east
    yy = yy + north

else

    ! UTM to longitude and latitude
    xx = xx - east
    yy = yy - north
    e1 = sqrt(1.d0 - e2)
    e1 = (1.d0 - e1)/(1.d0 + e1)
    rm = yy/scale_factor
    u = 1.d0 - e2/4.d0 - 3.d0*e4/64.d0 - 5.d0*e6/256.d0
    u = rm/(semimajor_wgs84*u)

    f1 = 3.d0*e1/2.d0 - 27.d0*e1**3.d0/32.d0
    f1 = f1*sin(2.d0*u)
    f2 = 21.d0*e1**2/16.d0 - 55.d0*e1**4/32.d0
    f2 = f2*sin(4.d0*u)
    f3 = 151.d0*e1**3.d0/96.d0
    f3 = f3*sin(6.d0*u)
    rlat1 = u + f1 + f2 + f3
    dlat1 = rlat1*r2d
    if (dlat1.ge.90.d0.or.dlat1.le.-90.d0) then
        dlat1 = dmin1(dlat1,90.d0)
        dlat1 = dmax1(dlat1,-90.d0)
        dlon = cm
    else
        c1 = ep2*cos(rlat1)**2
        t1 = tan(rlat1)**2
        f1 = 1.d0 - e2*sin(rlat1)**2
        rn1 = semimajor_wgs84/sqrt(f1)
        r1 = semimajor_wgs84*(1.d0 - e2)/sqrt(f1**3)
        d = xx/(rn1*scale_factor)

        f1 = rn1*tan(rlat1)/r1
        f2 = d**2/2.d0
        f3 = 5.d0*3.d0*t1 + 10.d0*c1 - 4.d0*c1**2 - 9.d0*ep2
        f3 = f3*d**2*d**2/24.d0
        f4 = 61.d0 + 90.d0*t1 + 298.d0*c1 + 45.d0*t1**2 - 252.d0*ep2 - 3.d0*c1**2
        f4 = f4*(d**2)**3.d0/720.d0
        rlat = rlat1 - f1*(f2 - f3 + f4)
        dlat = rlat*r2d

        f1 = 1.d0 + 2.d0*t1 + c1
        f1 = f1*d**2*d/6.d0
        f2 = 5.d0 - 2.d0*c1 + 28.d0*t1 - 3.d0*c1**2 + 8.d0*ep2 + 24.d0*t1**2
        f2 = f2*(d**2)**2*d/120.d0
        rlon = cmr + (d - f1 + f2)/cos(rlat1)
        dlon = rlon*r2d
        if (dlon.lt.-180.d0) then
            dlon = dlon + 360.d0
        endif
        if (dlon.gt.180.d0) then
            dlon = dlon - 360.d0
        endif
    endif
endif

if (iway.eq.1) then
    rlon4 = dlon
    rlat4 = dlat
    rx4 = rx_save
    ry4 = ry_save
else
    rx4 = xx
    if (lsouth) then
        yy = yy + 1.0d7
    endif
    ry4 = yy
    rlon4 = rlon_save
    rlat4 = rlat_save
endif

return
end subroutine

end module
