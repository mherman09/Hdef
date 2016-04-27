      PROGRAM grid
C----
C Create grids primarily for use with o92util
C----
      IMPLICIT NONE
      REAL*8 pi,d2r
      PARAMETER (pi=4.0d0*datan(1.0d0),d2r=pi/180.0d0)
      CHARACTER*40 ofile
      REAL*8 x1,x2,y1,y2,dx,dy,dz,dd,ref(5),dist,az,lo,la
      REAL*8 x,y,z
      INTEGER i,j,nx,ny,p,xsec,xz

C----
C Get parameters from the command line
C----
      call gcmdln(x1,x2,nx,dx,y1,y2,ny,dy,z,ofile,p,ref,xsec,xz)
C Check that limits and increments are properly defined
      if (x2.le.x1) call usage('!! Error: X2 must be greater than X1')
      if (y2.lt.y1) call usage('!! Error: Y2 must be greater than Y1')
      if (dx.le.0.0d0.and.nx.le.0.0d0) then
          call usage('!! Error: DX/NX unspecified')
      endif
      if (dy.le.0.0d0.and.ny.le.0.0d0) then
          !call usage('!! Error: DY/NY unspecified')
          ny = 1
      endif
      if (ref(5).gt.8.5d1) then
          call usage('!! Error: fault dip greater than 85')
      endif
C If output file unspecified, print result to standard output
      if (ofile.eq.'none'.and.p.eq.0) p = 1

C----
C Calculate number of points and grid increments
C----
      if (dx.gt.0.0d0) then
          nx = int((x2-x1)/dx+1.0d-6)+1
      elseif (nx.le.1) then
          call usage('!! Error: NX must be 2 or greater')
      else
          dx = (x2-x1)/dble(nx-1)
      endif

      if (dabs(y2-y1).lt.1.0d-6) then
          ny = 1
      elseif (dy.gt.0.0d0) then
          ny = int((y2-y1)/dy+1.0d-6)+1
      elseif (ny.lt.1) then
          call usage('!! Error: NY must be 1 or greater')
      elseif (ny.eq.1) then
          write(0,*) '!! Warning: Using NY=1 with Y1 not equal to Y2'
          dy = 0.0d0
      else
          dy = (y2-y1)/dble(ny-1)
      endif

C----
C Generate grid
C----
      if (p.eq.0) open(unit=101,file=ofile,status='unknown')
      do 16 i = 0,nx-1
          x = x1 + dble(i)*dx
          if (xsec.eq.1) then
              call dlola(lo,la,ref(1),ref(2),x,ref(3))
          endif
          do 15 j = 0,ny-1
              y = y1 + dble(j)*dy
C             CALCULATE Z ON DIPPING GRID
              if (ref(5).gt.0.0d0) then
                  call ddistaz(dist,az,ref(1),ref(2),x,y)
                  dist = dist*6.371d3
                  az = (ref(4)+90.0d0)*d2r-az
                  dd = dcos(az)*dist
                  dz = dd*dtan(ref(5)*d2r)
                  z = ref(3) + dz
              endif
              if (p.eq.0.and.xsec.eq.0) then
                  write(101,9999) x,y,z
              elseif (p.eq.1.and.xsec.eq.0) then
                  write(*,9999) x,y,z
              elseif (p.eq.0.and.xsec.eq.1.and.xz.eq.0) then
                  write(101,9999) lo,la,y
              elseif (p.eq.1.and.xsec.eq.1.and.xz.eq.0) then
                  write(*,9999) lo,la,y
              elseif (p.eq.0.and.xsec.eq.1.and.xz.eq.1) then
                  write(101,9999) x,y
              elseif (p.eq.1.and.xsec.eq.1.and.xz.eq.1) then
                  write(*,9999) x,y
              endif
   15     continue
   16 continue
      stop

 9999 format (3F16.8)

      END 

C======================================================================C

      SUBROUTINE ddistaz(dist,az,lon1,lat1,lon2,lat2)
C----
C Given the latitude and longitude of two points, compute the great
C circle distance between them (in radians), and the azimuth from
C point 1 to point 2, in radians clockwise from north.
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
     1                  dabs(mod(lon1-lon2,1.8d2)).le.1.0d-6.and.
     2                  dabs(lat1).gt.1.0d-6.and.
     3                  dabs(lon1).gt.1.0d-6) then
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

      lat2 = lat2*r2d
      lon2 = lon2*r2d
C     Return input values in initial units
      dist = dist*6371.0d0
      az = az*r2d
      lat1 = lat1*r2d
      lon1 = lon1*r2d
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE gcmdln(x1,x2,nx,dx,y1,y2,ny,dy,z,ofile,p,ref,xsec,xz)
      IMPLICIT none
      CHARACTER*40 tag,ofile
      REAL*8 x1,x2,y1,y2,z,dx,dy,ref(5)
      INTEGER narg,i,j,nx,ny,p,xsec,xz
C Initialize variable values
      x1 = 0.0d0
      x2 = 0.0d0
      nx = -1
      dx = -1.0d0
      y1 = 0.0d0
      y2 = 0.0d0
      ny = -1
      dy = -1.0d0
      z = 0.0d0
      ofile = 'none'
      p = 0
      xsec = 0
      xz = 0
      do 100 i = 1,5
          ref(i) = -1.0d0
  100 continue
      narg = iargc()
      if (narg.eq.0) then
          call usage('!! Error: no command line arguments specified')
      endif
      i = 0
   11 i = i + 1
      if (i.gt.narg) goto 12
      call getarg(i,tag)
      if (tag(1:5).eq.'-xsec') then
          xsec = 1
          do 101 j = 1,3
              call getarg(i+j,tag)
              read (tag,'(BN,F12.0)') ref(j)
  101     continue
          i = i + 3
      elseif (tag(1:3).eq.'-xz') then
          xz = 1
      elseif (tag(1:2).eq.'-x') then
          i = i + 1
          call getarg(i,tag)
          read (tag,'(BN,F12.0)') x1
          i = i + 1
          call getarg(i,tag)
          read (tag,'(BN,F12.0)') x2
      elseif (tag(1:2).eq.'-y') then
          i = i + 1
          call getarg(i,tag)
          read (tag,'(BN,F12.0)') y1
          i = i + 1
          call getarg(i,tag)
          read (tag,'(BN,F12.0)') y2
      elseif (tag(1:2).eq.'-z') then
          i = i + 1
          call getarg(i,tag)
          read (tag,'(BN,F12.0)') z
      elseif (tag(1:3).eq.'-nx') then
          i = i + 1
          call getarg(i,tag)
          read (tag,'(BN,I7)') nx
      elseif (tag(1:3).eq.'-ny') then
          i = i + 1
          call getarg(i,tag)
          read (tag,'(BN,I7)') ny
      elseif (tag(1:3).eq.'-dx') then
          i = i + 1
          call getarg(i,tag)
          read (tag,'(BN,F12.0)') dx
      elseif (tag(1:3).eq.'-dy') then
          i = i + 1
          call getarg(i,tag)
          read (tag,'(BN,F12.0)') dy
      elseif (tag(1:4).eq.'-dip') then
          do 102 j = 1,5
              call getarg(i+j,tag)
              read (tag,'(BN,F12.0)') ref(j)
  102     continue
          i = i + 5
      elseif (tag(1:2).eq.'-o') then
          i = i + 1
          call getarg(i,ofile)
      elseif (tag(1:2).eq.'-p') then
          p = 1
      elseif (tag(1:2).eq.'-h') then
          call usage(' ')
      else
          call usage('!! Error: no option '//tag)
      endif
      goto 11
   12 continue
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE usage(str)
      IMPLICIT none
      INTEGER lstr
      CHARACTER str*(*)
      if (str.ne.' ') then
          lstr = len(str)
          write(*,*) str(1:lstr)
          write(*,*)
      endif
      write(*,*)
     1 'Usage: grid -x X1 X2 [-nx NX|-dx DX] ',
     2            '-y Y1 Y2 [-ny NY|-dy DY]'
      write(*,*)
     1 '            [-z Z|-dip X0 Y0 Z0 STR ',
     2              'DIP] [-xsec X0 Y0 AZ] [-xz]'
      write(*,*)
     1 '            [-o OFILE] [-p] [-h]'
      write(*,*)
      write(*,*)
     1 '-x X1 X2    First column (x) limits'
      write(*,*)
     1 '-nx NX      Number of x grid points'
      write(*,*)
     1 '-dx DX      Increment in x direction (overrides -nx)'
      write(*,*)
     1 '-y Y1 Y2    Second column (y) limits'
      write(*,*)
     1 '-ny NY      Number of y grid points'
      write(*,*)
     1 '-dy DY      Increment in y direction (overrides -ny)'
      write(*,*)
     1 '-z Z        Third column (z) value'
      write(*,*)
     1 '-dip X0 Y0 Z0 STR DIP   Put grid onto plane with STR/DIP ',
     2                'and containing (X0,Y0,Z0)'
      write(*,*)
     1 '-xsec X0 Y0 AZ          Create a vertical cross section ',
     2                 'through (X0,Y0) with strike AZ'
      write(*,*)
     1 '-xz                     Print cross section x-z instead of ',
     2                            'lon lat z'
      write(*,*)
     1 '-o OFILE    Output to file (default prints to standard ',
     2              'output)'
      write (*,*)
     1 '-p          Print results to standard output (overrides -o)'
      write (*,*)
     1 '-h/-?       Online help (this screen)'
      write (*,*)
      STOP
      END
