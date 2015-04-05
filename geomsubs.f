      SUBROUTINE distaz(dist,az,lon1,lat1,lon2,lat2)
C----
C Given the latitude and longitude of two points, compute the great
C circle distance between them (in radians), and the azimuth from
C point 1 to point 2, in radians clockwise from north.
C
C (single precision)
C----
      IMPLICIT none
      REAL pi,d2r
      PARAMETER (pi=3.14159265,d2r=pi/1.8e2)
      REAL lon1,lat1,lon2,lat2,colat1,colat2,dlon,dlat,a
      REAL dist,az
C----
C Check if points are polar opposite
C----
      if (abs(lat1+lat2).le.1.0e-6.and.
     1                  abs(mod(lon1-lon2,1.8e2)).le.1.0e-6) then
          dist = pi
          az = 0.0e0
          goto 11
      endif
C----
C   12/13/12 - use the Haversine formula to get distance:
C              more accurate over short distances
C----
      colat1 = (9.0e1-lat1)*d2r
      colat2 = (9.0e1-lat2)*d2r
      dlon = (lon2-lon1)*d2r
      dlat = (lat2-lat1)*d2r
C----
C Spherical Law of Cosines to get distance
C----
C      arg = cos(colat1)*cos(colat2) +
C     1            sin(colat1)*sin(colat2)*cos(dlon)
C      if (arg.gt.1.0) arg = 1.0
C      dist = acos(arg)
C
C----
C Haversine Formula to get distance
C----
      a = sin(dlat/2.0e0)*sin(dlat/2.0e0) +
     1   cos(lat1*d2r)*cos(lat2*d2r)*sin(dlon/2.0e0)*sin(dlon/2.0e0)
      if (a.ge.1.0) then
          dist = 0.0e0
      else
          dist = 2.0e0*atan2(sqrt(a),sqrt(1.0-a))
      endif
C----
C Spherical Law of Sines to get azimuth
C----
      if (dist.lt.1.0e-6) then
          az = 0.0
          goto 11
      endif

      if (lat2.gt.lat1) then
          a = sin(colat2)*sin(dlon)/sin(dist)
          if (a.gt.1.0) a = 1.0e0
          if (a.lt.-1.0) a = -1.0e0
          az = asin(a)
      elseif (lat2.lt.lat1) then
          a = sin(pi-colat2)*sin(dlon)/sin(dist)
          if (a.gt.1.0) a = 1.0e0
          if (a.lt.-1.0) a = -1.0e0
          az = asin(a)
          az = pi - az
      elseif (lon1.lt.lon2) then
          az = pi/2.0e0
      else
          az = -pi/2.0e0
      endif
      if (az.lt.0.0) az = az+2.0e0*pi                  
      
   11 continue

      RETURN
      END

c----------------------------------------------------------------------c

      SUBROUTINE ddistaz(dist,az,lon1,lat1,lon2,lat2)
C----
C Given the latitude and longitude of two points, compute the great
C circle distance between them (in radians), and the azimuth from
C point 1 to point 2, in radians clockwise from north.
C
C (Double precision)
C----
      IMPLICIT none
      REAL*8 pi,d2r
      PARAMETER (pi=4.0d0*datan(1.0d0),d2r=pi/1.8d2)
      REAL*8 lon1,lat1,lon2,lat2,colat1,colat2,dlon,dlat,a
      REAL*8 dist,az
C----
C Check if points are polar opposite
C----
      if (dabs(lat1+lat2).le.1.0d-6.and.
     1                  dabs(mod(lon1-lon2,1.8d2)).le.1.0d-6) then
          dist = pi
          az = 0.0d0
          goto 11
      endif
C----
C   12/13/12 - use the Haversine formula to get distance:
C              more accurate over short distances
C----
      colat1 = (90.d0-lat1)*d2r
      colat2 = (90.d0-lat2)*d2r
      dlon = (lon2-lon1)*d2r
      dlat = (lat2-lat1)*d2r
C----
C Spherical Law of Cosines to get distance
C----
C      arg = cos(colat1)*cos(colat2) +
C     1            sin(colat1)*sin(colat2)*cos(dlon)
C      if (arg.gt.1.0) arg = 1.0
C      dist = acos(arg)
C----
C Haversine Formula to get distance
C----
      a = dsin(dlat/2.0d0)*dsin(dlat/2.0d0) +
     1   dcos(lat1*d2r)*dcos(lat2*d2r)*dsin(dlon/2.0d0)*dsin(dlon/2.0d0)
      if (a.ge.1.0d0) then
          dist = 0.0d0
      else
          dist = 2.0d0*datan2(dsqrt(a),dsqrt(1.0d0-a))
      endif
C----
C Spherical Law of Sines to get azimuth
C----
      if (dist.lt.1.0d-6) then
          az = 0.0d0
          goto 11
      endif

C      if (lat2.gt.lat1) then
C          a = dsin(colat2)*dsin(dlon)/dsin(dist)
C          if (a.gt.1.0) a = 1.0d0
C          if (a.lt.-1.0) a = -1.0d0
C          az = asin(a)
C      elseif (lat2.lt.lat1) then
C          a = dsin(pi-colat2)*dsin(dlon)/dsin(dist)
C          if (a.gt.1.0) a = 1.0d0
C          if (a.lt.-1.0) a = -1.0d0
C          az = dasin(a)
C          az = pi - az
C      elseif (lon1.lt.lon2) then
C          az = pi/2.0d0
C      else
C          az = -pi/2.0d0
C      endif
C      if (az.lt.0.0) az = az+2.0d0*pi   
      az = datan2(dsin(dlon)*dcos(lat2*d2r),
     1                      dcos(lat1*d2r)*dsin(lat2*d2r)
     2                       - dsin(lat1*d2r)*dcos(lat2*d2r)*dcos(dlon))
   11 continue

      RETURN
      END

c----------------------------------------------------------------------c

      SUBROUTINE dlola(lon2,lat2,lon1,lat1,dist,az)
C----
C Subroutine for computing final (lon,lat) from (lon,lat), (dist,az)
C Units: dist (km), az (deg)
C----
      IMPLICIT none
      REAL*8 pi,d2r,r2d
      PARAMETER (pi=4.0d0*atan(1.0d0),d2r=pi/1.8d2,r2d=1.8d2/pi)
      REAL*8 lon1,lat1,lon2,lat2,dist,az

      dist = dist/6371.0d0 ! km -> rad
      az = az*d2r          ! deg -> rad
      lat1 = lat1*d2r      ! deg -> rad
      lon1 = lon1*d2r      ! deg -> rad

      lat2 = dasin(dsin(lat1)*dcos(dist)+dcos(lat1)*dsin(dist)*dcos(az))

      lon2 = lon1 + datan2(dsin(az)*dsin(dist)*dcos(lat1),
     1                           dcos(dist)-dsin(lat1)*dsin(lat2))

C     Return input values in initial units
      dist = dist*6371.0d0
      az = az*r2d
      lat1 = lat1*r2d
      lon1 = lon1*r2d

      RETURN
      END

c----------------------------------------------------------------------c

      SUBROUTINE lola(lon2,lat2,lon1,lat1,dist,az)
      IMPLICIT none
      REAL lon1,lat1,lon2,lat2,dist,az

      dist = dist/6371.0
      az = az*3.14159265/180.0
      lat1 = lat1*3.14159265/180.0
      lon1 = lon1*3.14159265/180.0

      lat2 = asin(sin(lat1)*cos(dist)+cos(lat1)*sin(dist)*cos(az))

      lon2 = lon1 + atan2(sin(az)*sin(dist)*cos(lat1),
     1                           cos(dist)-sin(lat1)*sin(lat2))


      RETURN
      END

c----------------------------------------------------------------------c

      SUBROUTINE lonlat2merc(x,y,lon,lat,c)

c Given lon/lat of a point, compute its location on a mercator map
c "c" is the scale factor, in km/radian      

      IMPLICIT none
      REAL x,y,lon,lat,c,pi,r2d,d2r
      PARAMETER (pi=3.14159265,r2d=180/pi,d2r=pi/180)

      lon = lon*d2r
      lat = lat*d2r

      x = c*lon
      y = c*log(tan(pi/4+lat/2))

      lon = lon*r2d
      lat = lat*r2d

      RETURN
      END

c----------------------------------------------------------------------c

      SUBROUTINE merc2lonlat(lon,lat,x,y,c)

c Given x,y position on a mercator projection, compute the lon/lat
c "c" is the scale factor, in km/radian

      IMPLICIT none
      REAL lon,lat,x,y,c,pi,r2d
      PARAMETER (pi=3.14159265,r2d=180.0/pi)

      lon = x/c
      lat = 2*atan(exp(y/c))-pi/2

      lon = lon*r2d
      lat = lat*r2d
      
      RETURN
      END

c----------------------------------------------------------------------c

      SUBROUTINE NE2rtheta(r,theta,N,E)

c Given north and east vector components, calculate radius and azimuth,
c CW from north

      IMPLICIT none
      REAL r,theta,N,E

      r = sqrt(N*N + E*E)
      theta = atan2(E,N)

      RETURN
      END

c----------------------------------------------------------------------c
