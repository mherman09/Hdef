module geom

public :: coords2distaz
public :: distaz2coords
public :: pnpoly

!--------------------------------------------------------------------------------------------------!
contains
!--------------------------------------------------------------------------------------------------!

subroutine coords2distaz(lon1,lat1,lon2,lat2,dist,az)
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

subroutine distaz2coords(lon1,lat1,dist,az,lon2,lat2,dist_unit)
!----
! Compute the longitude and latitude a distance and azimuth clockwise from north away from an
! initial point. Input coordinates and azimuth are in degrees, input distance is in km, and output
! coordinates are in radians.
!----

use io, only: stderr
use trig, only: d2r
use earth, only: radius_earth_m, radius_earth_km

implicit none

! Arguments
double precision :: lon1, lat1, lon2, lat2, dist, az
character(len=*) :: dist_unit

! Local variables
double precision :: distr, azr, lon1r, lat1r, lat2r

! Convert all input quantities to radians
if (trim(dist_unit).eq.'km'.or.trim(dist_unit).eq.'kilometers') then
    distr = dist/radius_earth_km ! km -> rad
elseif (trim(dist_unit).eq.'m'.or.trim(dist_unit).eq.'meters') then
    distr = dist/radius_earth_m ! m -> rad
elseif (trim(dist_unit).eq.'rad'.or.trim(dist_unit).eq.'radians') then
    distr = dist
elseif (trim(dist_unit).eq.'deg'.or.trim(dist_unit).eq.'degrees') then
    distr = dist*d2r ! degrees -> rad
else
    write(stderr,*) 'distaz2coords: no distance units called "',trim(dist_unit),'"'
    distr = 0.0d0
endif
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

!--------------------------------------------------------------------------------------------------!

subroutine pnpoly(PX,PY,XX,YY,N,INOUT)

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
!   YY      - N LONG VECTOR CONTAING Y-COORDINATES OF
!             VERTICES OF POLYGON.
!   N       - NUMBER OF VERTICES IN THE POLYGON.
!   INOUT   - THE SIGNAL RETURNED:
!             -1 IF THE POINT IS OUTSIDE OF THE POLYGON,
!              0 IF THE POINT IS ON AN EDGE OR AT A VERTEX,
!              1 IF THE POINT IS INSIDE OF THE POLYGON.
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
double precision :: PX, PY, XX(N), YY(N)
integer :: N, INOUT

! Local variables
double precision :: X(N), Y(N), var
integer :: I, J
logical :: MX, MY, NX, NY
double precision, parameter :: eps = 1.0d-8


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
