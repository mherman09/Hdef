      SUBROUTINE distaz(dist,az,lon1,lat1,lon2,lat2)
!----
! Given the latitude and longitude of two points, compute the great
! circle distance between them (in radians), and the azimuth from
! point 1 to point 2, in radians clockwise from north.
!
! (single precision)
!----
      IMPLICIT none
      REAL pi,d2r
      PARAMETER (pi=3.14159265,d2r=pi/1.8e2)
      REAL lon1,lat1,lon2,lat2,colat1,colat2,dlon,dlat,a
      REAL dist,az
!----
! Check if points are polar opposite
!----
      if (abs(lat1+lat2).le.1.0e-6.and.
     1                  abs(mod(lon1-lon2,1.8e2)).le.1.0e-6) then
          dist = pi
          az = 0.0e0
          goto 11
      endif
!----
!   12/13/12 - use the Haversine formula to get distance:
!              more accurate over short distances
!----
      colat1 = (9.0e1-lat1)*d2r
      colat2 = (9.0e1-lat2)*d2r
      dlon = (lon2-lon1)*d2r
      dlat = (lat2-lat1)*d2r
!----
! Spherical Law of Cosines to get distance
!----
!      arg = cos(colat1)*cos(colat2) +
!     1            sin(colat1)*sin(colat2)*cos(dlon)
!      if (arg.gt.1.0) arg = 1.0
!      dist = acos(arg)
!
!----
! Haversine Formula to get distance
!----
      a = sin(dlat/2.0e0)*sin(dlat/2.0e0) +
     1   cos(lat1*d2r)*cos(lat2*d2r)*sin(dlon/2.0e0)*sin(dlon/2.0e0)
      if (a.ge.1.0) then
          dist = 0.0e0
      else
          dist = 2.0e0*atan2(sqrt(a),sqrt(1.0-a))
      endif
!----
! Spherical Law of Sines to get azimuth
!----
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

!----------------------------------------------------------------------c

      SUBROUTINE ddistaz(dist,az,lon1,lat1,lon2,lat2)
!----
! Given the latitude and longitude of two points, compute the great
! circle distance between them (in radians), and the azimuth from
! point 1 to point 2, in radians clockwise from north.
!
! (Double precision)
!----
      IMPLICIT none
      REAL*8 pi,d2r
      PARAMETER (pi=4.0d0*datan(1.0d0),d2r=pi/1.8d2)
      REAL*8 lon1,lat1,lon2,lat2,colat1,colat2,dlon,dlat,a
      REAL*8 dist,az
!----
! Check if points are polar opposite
!----
      if (dabs(lat1+lat2).le.1.0d-6.and.
     1                  dabs(mod(lon1-lon2,1.8d2)).le.1.0d-6) then
          dist = pi
          az = 0.0d0
          goto 11
      endif
!----
!   12/13/12 - use the Haversine formula to get distance:
!              more accurate over short distances
!----
      colat1 = (90.d0-lat1)*d2r
      colat2 = (90.d0-lat2)*d2r
      dlon = (lon2-lon1)*d2r
      dlat = (lat2-lat1)*d2r
!----
! Spherical Law of Cosines to get distance
!----
!      arg = cos(colat1)*cos(colat2) +
!     1            sin(colat1)*sin(colat2)*cos(dlon)
!      if (arg.gt.1.0) arg = 1.0
!      dist = acos(arg)
!----
! Haversine Formula to get distance
!----
      a = dsin(dlat/2.0d0)*dsin(dlat/2.0d0) +
     1   dcos(lat1*d2r)*dcos(lat2*d2r)*dsin(dlon/2.0d0)*dsin(dlon/2.0d0)
      if (a.ge.1.0d0) then
          dist = 0.0d0
      else
          dist = 2.0d0*datan2(dsqrt(a),dsqrt(1.0d0-a))
      endif
!----
! Spherical Law of Sines to get azimuth
!----
      if (dist.lt.1.0d-6) then
          az = 0.0d0
          goto 11
      endif

!      if (lat2.gt.lat1) then
!          a = dsin(colat2)*dsin(dlon)/dsin(dist)
!          if (a.gt.1.0) a = 1.0d0
!          if (a.lt.-1.0) a = -1.0d0
!          az = asin(a)
!      elseif (lat2.lt.lat1) then
!          a = dsin(pi-colat2)*dsin(dlon)/dsin(dist)
!          if (a.gt.1.0) a = 1.0d0
!          if (a.lt.-1.0) a = -1.0d0
!          az = dasin(a)
!          az = pi - az
!      elseif (lon1.lt.lon2) then
!          az = pi/2.0d0
!      else
!          az = -pi/2.0d0
!      endif
!      if (az.lt.0.0) az = az+2.0d0*pi   
      az = datan2(dsin(dlon)*dcos(lat2*d2r),
     1                      dcos(lat1*d2r)*dsin(lat2*d2r)
     2                       - dsin(lat1*d2r)*dcos(lat2*d2r)*dcos(dlon))
   11 continue

      RETURN
      END

!----------------------------------------------------------------------c

      SUBROUTINE dlola(lon2,lat2,lon1,lat1,dist,az)
!----
! Subroutine for computing final (lon,lat) from (lon,lat), (dist,az)
! Units: dist (km), az (deg)
!----
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

!     Return input values in initial units
      dist = dist*6371.0d0
      az = az*r2d
      lat1 = lat1*r2d
      lon1 = lon1*r2d

      RETURN
      END

!----------------------------------------------------------------------c

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

!----------------------------------------------------------------------c

      SUBROUTINE lonlat2merc(x,y,lon,lat,c)

! Given lon/lat of a point, compute its location on a mercator map
! "c" is the scale factor, in km/radian      

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

!----------------------------------------------------------------------c

      SUBROUTINE merc2lonlat(lon,lat,x,y,c)

! Given x,y position on a mercator projection, compute the lon/lat
! "c" is the scale factor, in km/radian

      IMPLICIT none
      REAL lon,lat,x,y,c,pi,r2d
      PARAMETER (pi=3.14159265,r2d=180.0/pi)

      lon = x/c
      lat = 2*atan(exp(y/c))-pi/2

      lon = lon*r2d
      lat = lat*r2d
      
      RETURN
      END

!----------------------------------------------------------------------c

      SUBROUTINE NE2rtheta(r,theta,N,E)

! Given north and east vector components, calculate radius and azimuth,
! CW from north

      IMPLICIT none
      REAL r,theta,N,E

      r = sqrt(N*N + E*E)
      theta = atan2(E,N)

      RETURN
      END

!----------------------------------------------------------------------c

      SUBROUTINE utmgeo(rlon4,rlat4,rx4,ry4,utmzon,iway)
!----
! Code modified from Specfem3d
!----
      IMPLICIT none
      REAL*8 pi,d2r,r2d
      PARAMETER (pi=3.14159265d0,d2r=pi/1.8d2,r2d=1.8d2/pi)
      REAL*8 rx4,ry4,rlon4,rlat4
      INTEGER utmzon,iway ! iway=1: utm->lonlat iway=2: lonlat->utm
! WGS84 (World Geodetic System 1984)
      REAL*8 semmaj,semmin
      PARAMETER (semmaj=6378137.0d0,semmin=6356752.314245d0)
      REAL*8 scfa
      PARAMETER (scfa=0.9996d0)
! UTM grids are actually Mercators which employ the standard UTM scale factor 0.9996 and set the Easting Origin to 500,000.
      REAL*8 north,east
      PARAMETER (north=0.0d0,east=500000.0d0)

      INTEGER zone
      REAL*8 rlon,rlat
      REAL*8 e2,e4,e6,ep2,xx,yy,dlat,dlon,cm,cmr,delam
      REAL*8 f1,f2,f3,f4,rm,rn,t,c,a,e1,u,rlat1,dlat1,c1,t1,rn1,r1,d
      REAL*8 rx_save,ry_save,rlon_save,rlat_save
      LOGICAL lsouth

      dlon = 0
      dlat = 0
      xx = 0
      yy = 0

      rlon_save = rlon4
      rlat_save = rlat4
      rx_save = rx4
      ry_save = ry4

      e2=1.d0-(semmin/semmaj)**2
      e4=e2*e2
      e6=e2*e4
      ep2=e2/(1.d0-e2)

      lsouth = .false.
      if (utmzon.lt.0) lsouth = .true.
      zone = abs(utmzon)
      cm = zone*6.0d0 - 183.d0
      cmr = cm*d2r

      if (iway.eq.1) then
          xx = rx4
          yy = ry4
          if (lsouth) yy = yy - 1.d7
      else
          dlat = rlat4
          dlon = rlon4
      endif
!
!---- Lat/Lon to UTM conversion
!
      if (iway.eq.2) then
          rlon = d2r*dlon
          rlat = d2r*dlat

          delam = dlon - cm
          if (delam.lt.-180.d0) delam = delam + 360.d0
          if (delam.gt.180.d0) delam = delam - 360.d0
          delam = delam*d2r

          f1 = (1.d0 - e2/4.d0 - 3.d0*e4/64.d0 - 5.d0*e6/256d0)*rlat
          f2 = 3.d0*e2/8.d0 + 3.d0*e4/32.d0 + 45.d0*e6/1024.d0
          f2 = f2*sin(2.d0*rlat)
          f3 = 15.d0*e4/256.d0*45.d0*e6/1024.d0
          f3 = f3*sin(4.d0*rlat)
          f4 = 35.d0*e6/3072.d0
          f4 = f4*sin(6.d0*rlat)
          rm = semmaj*(f1 - f2 + f3 - f4)
          if (dabs(dlat-90.d0).lt.1.0d-6.or.
     1                          dabs(dlat+90.d0).lt.1.0d-6) then
              xx = 0.d0
              yy = scfa*rm
          else
              rn = semmaj/sqrt(1.d0 - e2*sin(rlat)**2)
              t = tan(rlat)**2
              c = ep2*cos(rlat)**2
              a = cos(rlat)*delam

              f1 = (1.d0 - t + c)*a**3/6.d0
              f2 = 5.d0 - 18.d0*t + t**2 + 72.d0*c - 58.d0*ep2
              f2 = f2*a**5/120.d0
              xx = scfa*rn*(a + f1 + f2)
              f1 = a**2/2.d0
              f2 = 5.d0 - t + 9.d0*c + 4.d0*c**2
              f2 = f2*a**4/24.d0
              f3 = 61.d0 - 58.d0*t + t**2 + 600.d0*c - 330.d0*ep2
              f3 = f3*a**6/720.d0
              yy = scfa*(rm + rn*tan(rlat)*(f1 + f2 + f3))
          endif
          xx = xx + east
          yy = yy + north
!
!---- UTM to Lat/Lon conversion
!
      else
          xx = xx - east
          yy = yy - north
          e1 = sqrt(1.d0 - e2)
          e1 = (1.d0 - e1)/(1.d0 + e1)
          rm = yy/scfa
          u = 1.d0 - e2/4.d0 - 3.d0*e4/64.d0 - 5.d0*e6/256.d0
          u = rm/(semmaj*u)

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
              rn1 = semmaj/sqrt(f1)
              r1 = semmaj*(1.d0 - e2)/sqrt(f1**3)
              d = xx/(rn1*scfa)

              f1 = rn1*tan(rlat1)/r1
              f2 = d**2/2.d0
              f3 = 5.d0*3.d0*t1 + 10.d0*c1 - 4.d0*c1**2 - 9.d0*ep2
              f3 = f3*d**2*d**2/24.d0
              f4 = 61.d0 + 90.d0*t1 + 298.d0*c1 + 45.d0*t1**2
     1                                         - 252.d0*ep2 - 3.d0*c1**2
              f4 = f4*(d**2)**3.d0/720.d0
              rlat = rlat1 - f1*(f2 - f3 + f4)
              dlat = rlat*r2d
      
              f1 = 1.d0 + 2.d0*t1 + c1
              f1 = f1*d**2*d/6.d0
              f2 = 5.d0 - 2.d0*c1 + 28.d0*t1 - 3.d0*c1**2 + 8.d0*ep2
     1                                                    + 24.d0*t1**2
              f2 = f2*(d**2)**2*d/120.d0
              rlon = cmr + (d - f1 + f2)/cos(rlat1)
              dlon = rlon*r2d
              if (dlon.lt.-180.d0) dlon = dlon + 360.d0
              if (dlon.gt.180.d0) dlon = dlon - 360.d0
          endif
      endif

!
!----- output
!
      if (iway.eq.1) then
          rlon4 = dlon
          rlat4 = dlat
          rx4 = rx_save
          ry4 = ry_save
      else
          rx4 = xx
          if (lsouth) yy = yy + 1.d7
          ry4 = yy
          rlon4 = rlon_save
          rlat4 = rlat_save
      endif

      RETURN
      END





