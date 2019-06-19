module geom

public :: lola2distaz
public :: distaz2lola

public :: str2vec
public :: strdip2normal
public :: strdip2updip
public :: normal2strike
public :: normal2updip

public :: perspective

! public :: nvec2sdvec
public :: pnpoly

!--------------------------------------------------------------------------------------------------!
contains
!--------------------------------------------------------------------------------------------------!

subroutine lola2distaz(lon1,lat1,lon2,lat2,dist,az)
!----
! Given the latitude and longitude of two points, compute the great circle distance between them
! and the azimuth from point 1 to point 2 measured clockwise from north.
!
! Inputs
!     lon1: starting longitude, in degrees
!     lat1: starting latitude, in degrees
!     lon2: ending longitude, in degrees
!     lat2: ending latitude, in degrees
!
! Outputs
!     dist: distance between points along sphere, in radians
!     az: azimuth pointing from point 1 to point 2, clockwise from north, in degrees
!----

use io, only: stderr
use trig, only: pi, d2r, r2d
implicit none

! Arguments
double precision :: dist, az, lon1, lat1, lon2, lat2

! Local variables
double precision :: lon1r, lat1r, lon2r, lat2r, colat1, colat2, dlon, dlat, a

dist = 0.0d0
az = 0.0d0

if (abs(lat1).gt.90.0d0) then
    write(stderr,*) 'distaz2lola: input latitude 1 is greater than 90'
    return
endif

if (abs(lat2).gt.90.0d0) then
    write(stderr,*) 'distaz2lola: input latitude 2 is greater than 90'
    return
endif

! Check if points are polar opposite
if (dabs(lat1+lat2).le.1.0d-6.and.dabs(mod((lon1-lon2)+1.8d2,3.6d2)).le.1.0d-6) then
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

az = az*r2d

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine distaz2lola(lon1,lat1,dist,az,lon2,lat2,dist_unit,az_unit,ierr)
!----
! Given a longitude and latitude and distance and azimuth clockwise from north away from this
! initial point, compute the ending longitude and latitude.
!
! Inputs
!     lon1: starting longitude, in degrees
!     lat1: starting latitude, in degrees
!     dist: distance between points along sphere, in dist_unit
!     az: azimuth pointing from point 1 to point 2, clockwise from north, in az_unit
!     dist_unit: radians, meters, kilometers
!     az_unit: degrees, radians
!     ierr: error status
!
! Outputs
!     lon2: ending latitude, in degrees
!     lat2: ending latitude, in degrees
!----

use io, only: stderr
use trig, only: d2r, r2d

implicit none

! Arguments
double precision :: lon1, lat1, lon2, lat2, dist, az
character(len=*) :: dist_unit, az_unit
integer :: ierr

! Local variables
double precision :: distr, azr, lon1r, lat1r


ierr = 0
lon2 = 0.0d0
lat2 = 0.0d0

if (abs(lat1).gt.90.0d0) then
    write(stderr,*) 'distaz2lola: input latitude is greater than 90'
    ierr = 1
    return
endif

! Convert all input quantities to radians
! Distance
if (dist_unit.eq.'radians') then
    distr = dist
else
    write(stderr,*) 'distaz2lola: no dist_unit named ',trim(dist_unit)
    ierr = 2
    return
endif

! Azimuth
if (az_unit.eq.'degrees') then
    azr = az*d2r          ! deg -> rad
elseif (az_unit.eq.'radians') then
    azr = az
else
    write(stderr,*) 'distaz2lola: no az_unit named ',trim(az_unit)
    ierr = 2
    return
endif

lon1r = lon1*d2r      ! deg -> rad
lat1r = lat1*d2r      ! deg -> rad

! Spherical law of cosines
lat2 = dasin(dsin(lat1r)*dcos(distr)+dcos(lat1r)*dsin(distr)*dcos(azr))

! Spherical law of cosines and law of sines
lon2 = lon1r + datan2(dsin(azr)*dsin(distr)*dcos(lat1r), dcos(distr)-dsin(lat1r)*dsin(lat2))

lon2 = lon2*r2d
lat2 = lat2*r2d

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine str2vec(str,strike_vec)
!----
! Given the strike (degrees CW from north) of a plane, calculate the vector parallel to strike
! (east, north, and vertical components).
!----

use trig, only: d2r

implicit none

! Arguments
double precision :: str, strike_vec(3)

strike_vec(1) = sin(str*d2r)
strike_vec(2) = cos(str*d2r)
strike_vec(3) = 0.0d0

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine strdip2normal(str,dip,normal_vec)
!----
! Given the strike (degrees CW from north) and dip (degrees from horizontal) of a plane, calculate
! the normal vector (east, north, and vertical components).
!----

use trig, only: d2r

implicit none

! Arguments
double precision :: str, dip, normal_vec(3)

normal_vec(1) = sin(dip*d2r)*sin((str+90.0d0)*d2r)
normal_vec(2) = sin(dip*d2r)*cos((str+90.0d0)*d2r)
normal_vec(3) = cos(dip*d2r)

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine strdip2updip(str,dip,updip_vec)
!----
! Given the strike (degrees CW from north) and dip (degrees from horizontal) of a plane, calculate
! the vector pointing up-dip (east, north, and vertical components).
!----

use trig, only: d2r

implicit none

! Arguments
double precision :: str, dip, updip_vec(3)

updip_vec(1) = -cos(dip*d2r)*sin((str+90.0d0)*d2r)
updip_vec(2) = -cos(dip*d2r)*cos((str+90.0d0)*d2r)
updip_vec(3) = sin(dip*d2r)

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine normal2strike(normal_vec,strike_vec)
use algebra, only: normalize
implicit none
! Arguments
double precision :: normal_vec(3), strike_vec(3)

if (abs(normal_vec(1)).lt.1.0d-6.and.abs(normal_vec(2)).lt.1.0d-6) then
    strike_vec(1) = 1.0d0
    strike_vec(2) = 0.0d0
    strike_vec(3) = 0.0d0
else
    strike_vec(1) = -normal_vec(2)
    strike_vec(2) = normal_vec(1)
    strike_vec(3) = 0.0d0
endif

call normalize(strike_vec)

return
end

!--------------------------------------------------------------------------------------------------!

subroutine normal2updip(normal_vec,updip_vec)
use trig, only: pi
implicit none
! Arguments
double precision :: normal_vec(3), updip_vec(3)

if (abs(normal_vec(1)).lt.1.0d-6.and.abs(normal_vec(2)).lt.1.0d-6) then
    updip_vec(1) = 0.0d0
    updip_vec(2) = 1.0d0
    updip_vec(3) = 0.0d0
else
    updip_vec(1) = normal_vec(3)*cos(atan2(normal_vec(2),normal_vec(1))+pi)
    updip_vec(2) = normal_vec(3)*sin(atan2(normal_vec(2),normal_vec(1))+pi)
    updip_vec(3) = sqrt(normal_vec(1)*normal_vec(1)+normal_vec(2)*normal_vec(2))
endif

return
end

!--------------------------------------------------------------------------------------------------!

!
! subroutine nvec2sdvec(normal_vector,strike_vector,dip_vector)
! !----
! ! Calculate the vectors parallel to strike and dip of a plane, given its normal vector.
! !----
!
! use io, only: stderr
!
! implicit none
!
! ! Arguments
! double precision :: normal_vector(3), strike_vector(3), dip_vector(3)
!
! ! Local variables
! double precision :: magnitude, strike_angle
!
! ! Initialize strike, dip vectors
! strike_vector = 0.0d0
! dip_vector = 0.0d0
!
! ! Normalize normal vector
! magnitude = normal_vector(1)*normal_vector(1) + &
!             normal_vector(2)*normal_vector(2) + &
!             normal_vector(3)*normal_vector(3)
!
! if (magnitude.lt.1.0d-10) then
!     write(stderr,*) 'normalvec2strdipvec: normal vector magnitude is zero'
!     return
! endif
!
! normal_vector = normal_vector/magnitude
!
! ! Strike is horizontal, parallel to plane
! strike_angle = atan2(normal_vector(1),normal_vector(2))
! strike_vector(1) = sin(strike_angle)
! strike_vector(2) = cos(strike_angle)
! strike_vector(3) = 0.0d0
!
! ! Up-dip is perpendicular to strike, normal
! dip_vector(1) = normal_vector(2)*strike_vector(3) - normal_vector(3)*strike_vector(2)
! dip_vector(2) = normal_vector(3)*strike_vector(1) - normal_vector(1)*strike_vector(3)
! dip_vector(3) = normal_vector(1)*strike_vector(2) - normal_vector(2)*strike_vector(1)
!
! return
! end subroutine

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

subroutine perspective(point_x,point_y,point_z,view_x,view_y,view_z,look_x,look_y,look_z,x,y)
!----
! Project a point in 3-D onto a 2-D plane, given viewing coordinates and the look vector of the
! viewer. The center of the plane corresponds to the vanishing point of the look direction.
!----

use algebra, only: dot_product, cross_product, normalize

implicit none

! Arguments
double precision :: point_x, point_y, point_z, view_x, view_y, view_z, look_x, look_y, look_z
double precision :: x, y

! Local variables
double precision :: vec(3), look_vec(3), dot, proj, ang_from_vanish, updip(3), strike(3), azimuth


! Calculate the angle between the look vector and the vector from the viewing location to the
! point of interest

! Vector viewer->point (point of interest vector)
vec(1) = point_x - view_x
vec(2) = point_y - view_y
vec(3) = point_z - view_z

! Look vector
look_vec(1) = look_x
look_vec(2) = look_y
look_vec(3) = look_z

! Dot product of point of interest vector and look vector
call dot_product(vec,look_vec,dot)
if (dot.lt.0.0d0) then
    ! point is behind the viewer
    x = 1e10
    y = 1e10
    return
endif

! Projection of point of interest vector on look vector direction
proj = dot/sqrt(look_x**2+look_y**2+look_z**2)

! Angle from vanishing point to point of interest
ang_from_vanish = proj/sqrt(vec(1)**2+vec(2)**2+vec(3)**2)
ang_from_vanish = acos(ang_from_vanish)

! Vector parallel to plane
vec = vec - proj*look_vec/sqrt(look_x**2+look_y**2+look_z**2)

! Up-dip and along-strike vectors
if (look_vec(3).ge.0.0d0) then
    call normal2updip(look_vec,updip)
    call normal2strike(look_vec,strike)
    strike = -strike
else
    call normal2updip(-look_vec,updip)
    call normal2strike(-look_vec,strike)
endif
call normalize(updip)
call normalize(strike)

! Project parallel vector onto strike and updip vectors
call dot_product(vec,strike,x)
call dot_product(vec,updip,y)
azimuth = atan2(y,x)

! Apparent position of projected point for viewer
x = tan(ang_from_vanish)*cos(azimuth)
y = tan(ang_from_vanish)*sin(azimuth)


return
end subroutine



!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

subroutine pnpoly(PX,PY,XX,YY,N,INOUT,eps)

!----
! SUBROUTINE pnpoly
! WRITTEN BY RANDOLPH FRANKLIN, UNIVERSITY OF OTTAWA, 7/70.
!
! PURPOSE: TO DETERMINE WHETHER A POINT IS INSIDE A POLYGON
!
! USAGE: CALL PNPOLY (PX, PY, XX, YY, N, INOUT )
!
! PARAMETERS:
!   PX      - X-COORDINATE OF POINT IN QUESTION.
!   PY      - Y-COORDINATE OF POINT IN QUESTION.
!   XX      - N LONG VECTOR CONTAINING X-COORDINATES OF
!             VERTICES OF POLYGON.
!   YY      - N LONG VECTOR CONTAINING Y-COORDINATES OF
!             VERTICES OF POLYGON.
!   N       - NUMBER OF VERTICES IN THE POLYGON.
!   INOUT   - THE SIGNAL RETURNED:
!             -1 IF THE POINT IS OUTSIDE OF THE POLYGON,
!              0 IF THE POINT IS ON AN EDGE OR AT A VERTEX,
!              1 IF THE POINT IS INSIDE OF THE POLYGON.
!   eps     - Threshold for distance from boundary
!
! REMARKS:
!   - THE VERTICES MAY BE LISTED CLOCKWISE OR ANTICLOCKWISE.
!   - THE FIRST MAY OPTIONALLY BE REPEATED, IF SO N MAY OPTIONALLY BE INCREASED BY 1.
!   - THE INPUT POLYGON MAY BE A COMPOUND POLYGON CONSISTING OF SEVERAL SEPARATE SUBPOLYGONS. IF SO,
!     THE FIRST VERTEX OF EACH SUBPOLYGON MUST BE REPEATED, AND WHEN CALCULATING N, THESE FIRST
!     VERTICES MUST BE COUNTED TWICE.
!   - INOUT IS THE ONLY PARAMETER WHOSE VALUE IS CHANGED.
!
! METHOD:
!   A VERTICAL LINE IS DRAWN THRU THE POINT IN QUESTION. IF IT CROSSES THE POLYGON AN ODD NUMBER OF
!   TIMES, THEN THE POINT IS INSIDE OF THE POLYGON.
!
! Modifications and translation to Modern Fortran by Matt Herman
!----

! use io, only: stderr
implicit none

! Arguments
double precision :: PX, PY, XX(N), YY(N), eps
integer :: N, INOUT

! Local variables
double precision :: X(N), Y(N), var
integer :: I, J
logical :: MX, MY, NX, NY


! Compute the coordinates of the points on the clipping path (xx,yy) relative to the point of
! interest (px,py)
do I = 1,N
    X(I) = XX(I)-PX
    Y(I) = YY(I)-PY

    ! Check whether point of interest is at a vertex of the polygon
    if (abs(X(I)).lt.eps.and.abs(Y(I)).lt.eps) then
        ! write(stderr,*) 'pnpoly: point of interest is at a vertex'
        INOUT = 0
        return
    endif
enddo


! Initialize point of interest as outside of polygon
INOUT = -1


! Loop over line segments in polygon to determine whether a vertical ray from above the
! point of interest crosses the line segment
do I = 1,N
    J = 1 + mod(I,N)

    ! Is the first point in the segment path right/east of of the point of interest?
    if (X(I).ge.0.0d0) then
        MX = .true.
    else
        MX = .false.
    endif

    ! Is the second point in the segment right/east of of the point of interest?
    if (X(J).ge.0.0d0) then
        NX = .true.
    else
        NX = .false.
    endif

    ! Is the first point in the segment above/north of of the point of interest?
    if (Y(I).ge.0.0d0) then
        MY = .true.
    else
        MY = .false.
    endif

    ! Is the second point in the segment above/north of of the point of interest?
    if (Y(J).ge.0.0d0) then
        NY = .true.
    else
        NY = .false.
    endif
    ! write(stderr,*) 'MX:',MX
    ! write(stderr,*) 'NX:',NX
    ! write(stderr,*) 'MY:',MY
    ! write(stderr,*) 'NY:',NY

    ! Check whether the point of interest is on a horizontal line
    if (dabs(Y(I)).lt.eps.and.dabs(Y(J)).lt.eps) then
        if (MX.and..not.NX .or. .not.MX.and.NX) then
            ! write(stderr,*) 'point of interest lies on this horizontal line segment'
            INOUT = 0
            return
        endif
    endif

    ! Check whether the point of interest is on a vertical line
    if (dabs(X(I)).lt.eps.and.dabs(X(J)).lt.eps) then
        if (MY.and..not.NY .or. .not.MY.and.NY) then
            ! write(stderr,*) 'point of interest lies on this vertical line segment'
            INOUT = 0
            return
        endif
    endif

    ! If the line segment is below, left, or right of the point of interest, go to the next segment
    if (MX.and.NX) then
        ! write(stderr,*) 'line segment is right of the point of interest'
        cycle
    elseif (.not.(MX.or.NX)) then
        ! write(stderr,*) 'line segment is left of the point of interest'
        cycle
    elseif (.not.(MY.or.NY)) then
        ! write(stderr,*) 'line segment is below the point of interest'
        cycle
    endif

    ! Check whether the line segment is above the point
    if (MY.and.NY.and..not.MX.and.NX) then
        ! write(stderr,*) 'line segment is above the point of interest'
        INOUT = -INOUT
        cycle
    else
        ! write(stderr,*) 'checking whether line segment is above the point of interest'
        var = (Y(I)*X(J)-X(I)*Y(J))/(X(J)-X(I))
    endif

    if (var.lt.-eps) then
        ! write(stderr,*) 'line segment is below the point of interest'
        cycle
    elseif (dabs(var).lt.eps) then
        ! write(stderr,*) 'point of interest is on the line segment'
        INOUT = 0
        return
    else
        ! write(stderr,*) 'line segment is above the point of interest'
        INOUT = -INOUT
    endif
enddo

return
end subroutine


end module
