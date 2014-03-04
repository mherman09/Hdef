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

      call gcmdln(ifile,ofile)

      if (ifile.eq.'none') then
          print *,'To see options, use -h flag'
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
      PARAMETER (pi=4.0d0*atan(1.0d0),d2r=pi/1.8d2,r2d=1.8d2/pi)

      call lonlat2xyz(ix,iy,iz,ilon,ilat)
      call lonlat2xyz(px,py,pz,plon,plat)

      ox = px*(px*ix+py*iy+pz*iz)*(1-cos(angle*d2r))
     1         + ix*cos(angle*d2r)
     2         + (-pz*iy+py*iz)*sin(angle*d2r)
      oy = py*(px*ix+py*iy+pz*iz)*(1-cos(angle*d2r))
     1         + iy*cos(angle*d2r)
     2         + ( pz*ix-px*iz)*sin(angle*d2r)
      oz = pz*(px*ix+py*iy+pz*iz)*(1-cos(angle*d2r))
     1         + iz*cos(angle*d2r)
     2         + (-py*ix+px*iy)*sin(angle*d2r)

      olon = atan2(oy,ox)*r2d
      olat = asin(oz)*r2d

      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE lonlat2xyz(x,y,z,lon,lat)

      IMPLICIT none
      REAL*8 pi,d2r
      PARAMETER (pi=4.0d0*atan(1.0d0),d2r=pi/1.8d2)
      REAL*8 lon,lat,x,y,z

      x = cos(lon*d2r)*cos(lat*d2r)
      y = sin(lon*d2r)*cos(lat*d2r)
      z = sin(lat*d2r)

      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE gcmdln(ifile,ofile)

      IMPLICIT NONE
      CHARACTER*30 ifile,ofile,tag
      INTEGER i,narg

      ifile = 'none'
      ofile = 'sphfinrot.out'

      narg = iargc()
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
     1 'Usage: sphfinrot -f [IFILE] -o [OFILE] -h/-?'
      write(*,*)
     1 '  -f [IFILE] (Default none) name of input file'
      write(*,*)
     1 '                            if no file name given, fault ',
     2                            'parameters prompted for manual input'
      write(*,*)
     1 '  -o [OFILE] (Default sphfinrot.out) name of output file'
      write (*,*)
     1 '  -h/-?        help'
      write (*,*) ''
      write (*,*)
     1 '  sphfinrot makes a counter-clockwise finite rotation of a',
     2     ' point around a pole'
      write (*,*) ''
      write (*,*)
     1 '    input file'
      write (*,*)
     1 '        lon_init lat_init pole_lon pole_lat angle'
      write (*,*)
     1 '         :  :'
      write (*,*)
     1 '    output file'
      write (*,*)
     1 '        lon_final lat_final'
      write (*,*)
     1 '         :  :'
      write (*,*) ''

      STOP
      END

