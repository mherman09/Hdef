      PROGRAM sphfinrot
C----
C Given a point (ilon,ilat), a pole of rotation (plon,plat), and a
C finite angle of rotation, in degrees, calculate the new point
C (olon,olat).
C
C Based on results found at:
C inside.mines.edu/~gmurray/ArbitraryAxisRotation
C----

      IMPLICIT none
      REAL*8 olon,olat,ilon,ilat,plon,plat,angle
      REAL*8 ix,iy,iz,px,py,pz,ox,oy,oz
      CHARACTER*30 ifile,ofile
      INTEGER user

      call gcmdln(ifile,ofile,user)

      if (user.eq.1) then
          print *,'Enter lon_init lat_init pole_lon pole_lat angle:'
          read *,ilon,ilat,plon,plat,angle
          call rotate(olon,olat,ilon,ilat,plon,plat,angle)
          write (*,9999), olon,olat
      else
          open (unit=11,file=ifile,status='old')
          open (unit=12,file=ofile,status='unknown')
  101     read (11,*,end=102) ilon,ilat,plon,plat,angle
              call rotate(olon,olat,ilon,ilat,plon,plat,angle)
              write (12,9999) olon,olat
              goto 101
  102     continue
      endif

 9999 format(2F14.6)

      END

C----------------------------------------------------------------------c

      SUBROUTINE rotate(olon,olat,ilon,ilat,plon,plat,angle)
      IMPLICIT NONE
      REAL*8 olon,olat,ilon,ilat,plon,plat,angle
      REAL*8 ix,iy,iz,px,py,pz,ox,oy,oz
      REAL*8 pi,d2r,r2d
      PARAMETER (pi=4.0d0*datan(1.0d0),d2r=pi/1.8d2,r2d=1.8d2/pi)

      call lonlat2xyz(ix,iy,iz,ilon,ilat)
      call lonlat2xyz(px,py,pz,plon,plat)

      ox = px*(px*ix+py*iy+pz*iz)*(1.0d0-dcos(angle*d2r))
     1         + ix*dcos(angle*d2r)
     2         + (-pz*iy+py*iz)*dsin(angle*d2r)
      oy = py*(px*ix+py*iy+pz*iz)*(1.0d0-dcos(angle*d2r))
     1         + iy*dcos(angle*d2r)
     2         + ( pz*ix-px*iz)*dsin(angle*d2r)
      oz = pz*(px*ix+py*iy+pz*iz)*(1.0d0-dcos(angle*d2r))
     1         + iz*dcos(angle*d2r)
     2         + (-py*ix+px*iy)*dsin(angle*d2r)

      olon = datan2(oy,ox)*r2d
      olat = dasin(oz)*r2d

      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE lonlat2xyz(x,y,z,lon,lat)

      IMPLICIT none
      REAL*8 pi,d2r
      PARAMETER (pi=4.0d0*datan(1.0d0),d2r=pi/1.8d2)
      REAL*8 lon,lat,x,y,z

      x = dcos(lon*d2r)*dcos(lat*d2r)
      y = dsin(lon*d2r)*dcos(lat*d2r)
      z = dsin(lat*d2r)

      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE gcmdln(ifile,ofile,user)

      IMPLICIT NONE
      CHARACTER*30 ifile,ofile,tag
      INTEGER i,narg,user

      user = 0
      ifile = 'sphfinrot.in'
      ofile = 'sphfinrot.out'

      narg = iargc()
      if (narg.eq.0) call usage()

      i = 0
 9998 i = i + 1
      if (i.gt.narg) goto 9999
          call getarg(i,tag)
          if (tag(1:2).eq.'-f') then
              i = i + 1
              call getarg(i,ifile)
          elseif (tag(1:2).eq.'-o') then
              i = i + 1
              call getarg(i,ofile)
          elseif (tag(1:2).eq.'-d') then
              write(*,*) 'Running with default file names'
          elseif (tag(1:2).eq.'-u') then
              user = 1
          elseif (tag(1:2).eq.'-h'.or.tag(1:2).eq.'-?') then
              call usage()
          endif
          goto 9998
 9999 continue

      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE usage()
      IMPLICIT none

      write(*,*)
     1 'Usage: sphfinrot -f [IFILE] -o [OFILE] -d -u -h/-?'
      write(*,*)
     1 '  -f [IFILE] (Default sphfinrot.in) name of input file'
      write(*,*)
     1 '  -o [OFILE] (Default sphfinrot.out) name of output file'
      write(*,*)
     1 '  -d         Run with default file names'
      write (*,*)
     1 '  -u         Prompt user to enter information through standard'
      write (*,*)
     1 '                 input for single calculation'
      write (*,*)
     1 '  -h/-?      Online help (this screen)'
      write (*,*) ''
      write (*,*)
     1 '  sphfinrot makes a counter-clockwise finite rotation of a',
     2     ' point around a pole'
      write (*,*) ''
      write (*,*)
     1 '    Input file'
      write (*,*)
     1 '        lon_init lat_init pole_lon pole_lat angle'
      write (*,*)
     1 '         :  :'
      write(*,*) ''
      write (*,*)
     1 '    Output file'
      write (*,*)
     1 '        lon_final lat_final'
      write (*,*)
     1 '         :  :'
      write (*,*) ''

      STOP
      END

