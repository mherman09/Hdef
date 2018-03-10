      PROGRAM grid
C----
C Create grids primarily for use with o92util
C----
      IMPLICIT NONE
      REAL*8 pi,d2r
      PARAMETER (pi=4.0d0*datan(1.0d0),d2r=pi/180.0d0)
      CHARACTER*40 ofile
      CHARACTER*200 clip,line
      REAL*8 x1,x2,y1,y2,dx,dy,dz,dd,ref(5),dist,az,lo,la
      REAL*8 x,y,z,xx(1000),yy(1000),zsave
      INTEGER i,j,k,nx,ny,p,xsec,xz,nn,io,datain,expflg

C----
C Get parameters from the command line
C     x1:    First defined x value (default: 0)
C     x2:    Second defined x value (default: 0)
C     nx:    Number of x increments (default: -1)
C     dx:    x increment, overrides nx (default: -1.0)
C     y1:    First defined y value (default: 0)
C     y2:    Second defined y value (default: 0)
C     ny:    Number of y increments (default: -1)
C     dy:    y increment, overrides ny (default: -1.0)
C     z:     z value
C     ofile: Name of output file
C     p:     Standard output flag (default: 0)
C     ref:   
C     xsec:  
C     xz:    
C     clip:  
C----
      call gcmdln(x1,x2,nx,dx,y1,y2,ny,dy,z,ofile,p,ref,xsec,xz,clip,
     1            datain,expflg)
      if (datain.eq.1) goto 9001
C Check that limits and increments are properly defined
      if (x2.lt.x1) call usage('!! Error: X2 must be greater than X1')
      if (y2.lt.y1) call usage('!! Error: Y2 must be greater than Y1')
      if (dx.le.0.0d0.and.nx.le.0.0d0) then
          call usage('!! Error: DX/NX unspecified')
      endif
      if (dx.le.0.0d0.and.nx.le.0.0d0) then
          nx = 1
      endif
      if (dy.le.0.0d0.and.ny.le.0.0d0) then
          !call usage('!! Error: DY/NY unspecified')
          ny = 1
      endif
      if (ref(5).gt.8.5d1) then
          call usage('!! Error: fault dip greater than 85')
      endif
 9001 if (datain.eq.1.and.clip.eq.'none') then
          call usage('!! Error: can only check if input is inside '//
     1               'polygon along with -clip')
      endif
C If output file unspecified, print result to standard output
      if (ofile.eq.'none'.and.p.eq.0) p = 1
C If clip file is defined, read it
      if (clip.ne.'none') then
          open(unit=201,file=clip,status='old')
          i = 1
 1001     read(201,*,end=1002) xx(i),yy(i)
              i = i + 1
              goto 1001
 1002     continue
          nn = i - 1
      endif
      close(201)

C----
C Check if standard input is inside polygon
C----
      ! standard input = 5
      if (datain.eq.1) then
 1003     read(5,'(A)',end=1004) line
              read(line,*) x,y
              call pnpoly(x,y,xx,yy,nn,io)
              if (io.eq.1) then
                  write(*,*) trim(line)
              endif
              goto 1003
 1004     continue
          stop
      endif

C----
C Calculate number of points and grid increments
C----
      ! Loop is over nx and incrementing is done with dx, so always calculate both
      if (dx.gt.0.0d0) then
          if (mod(expflg,10).eq.0) then
              nx = int((x2-x1)/dx+1.0d-6)+1
          else
              if (dx.le.1.0d0) then
                  call usage('!! Error: in exponential mode, DX must '//
     1                       'be greater than 1')
              endif
              nx = int(dlog(x2/x1)/dlog(dx))+1
          endif
      elseif (nx.le.1) then
          call usage('!! Error: NX must be 2 or greater')
      else
          if (mod(expflg,10).eq.0) then
              dx = (x2-x1)/dble(nx-1)
          else
              dx = dexp(dlog(x2/x1)/dble(nx-1))
          endif
      endif

      ! Loop is over ny and incrementing is done with dy, so always calculate both
      if (dabs(y2-y1).lt.1.0d-6) then
          ny = 1
      elseif (dy.gt.0.0d0) then
          if (expflg.lt.10) then
              ny = int((y2-y1)/dy+1.0d-6)+1
          else
              if (dy.le.1.0d0) then
                  call usage('!! Error: in exponential mode, DY must '//
     1                       'be greater than 1')
              endif
              ny = int(dlog(y2/y1)/dlog(dy))+1
          endif
      elseif (ny.lt.1) then
          call usage('!! Error: NY must be 1 or greater')
      elseif (ny.eq.1) then
          write(0,*) '!! Warning: Using NY=1 with Y1 not equal to Y2'
          dy = 0.0d0
      else
          if (expflg.lt.10) then
              dy = (y2-y1)/dble(ny-1)
          else
              dy = dexp(dlog(y2/y1)/dble(ny-1))
          endif
      endif

C----
C Generate grid
C----
      zsave = z
      if (p.eq.0) open(unit=101,file=ofile,status='unknown')
      do 16 i = 0,nx-1
          ! Using equal or fractional spacing?
          if (mod(expflg,10).eq.0) then
              x = x1 + dble(i)*dx
          else
              x = x1
              do 201 k = 1,i
                  x = x*dx
  201         continue
          endif
          if (xsec.eq.1) then
              call dlola(lo,la,ref(1),ref(2),x,ref(3))
          endif
          do 15 j = 0,ny-1
              if (expflg.lt.10) then
                  y = y1 + dble(j)*dy
              else
                  y = y1
                  do 202 k = 1,j
                      y = y*dy
  202             continue
              endif
C             CALCULATE Z ON DIPPING GRID
              if (ref(5).gt.0.0d0) then
                  call ddistaz(dist,az,ref(1),ref(2),x,y)
                  dist = dist*6.371d3
                  az = (ref(4)+90.0d0)*d2r-az
                  dd = dcos(az)*dist
                  dz = dd*dtan(ref(5)*d2r)
                  z = ref(3) + dz
              endif
              if (clip.ne.'none') then
                  zsave = z
                  call pnpoly(x,y,xx,yy,nn,io)
                  if (io.eq.-1) z = -12345d0
              endif
              !!!!!!!!!!!!!!!!!!if (dabs(x).gt.1e6.or.dabs(y).gt.1e6) then
              if (p.eq.0.and.xsec.eq.0) then
                  if (dabs(x).lt.1e5.and.dabs(y).lt.1e5) then
                      write(101,9999) x,y,z
                  else
                      write(101,9998) x,y,z
                  endif
              elseif (p.eq.1.and.xsec.eq.0) then
                  if (dabs(x).lt.1e5.and.dabs(y).lt.1e5) then
                      write(*,9999) x,y,z
                  else
                      write(*,9998) x,y,z
                  endif
              elseif (p.eq.0.and.xsec.eq.1.and.xz.eq.0) then
                  write(101,9999) lo,la,y
              elseif (p.eq.1.and.xsec.eq.1.and.xz.eq.0) then
                  write(*,9999) lo,la,y
              elseif (p.eq.0.and.xsec.eq.1.and.xz.eq.1) then
                  write(101,9999) x,y
              elseif (p.eq.1.and.xsec.eq.1.and.xz.eq.1) then
                  write(*,9999) x,y
              endif
              z = zsave
   15     continue
   16 continue
      stop

 9998 format (3(1PE16.8))
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
C                                                                       
C        SUBROUTINE pnpoly                                              
C                                                                       
C        PURPOSE                                                        
C           TO DETERMINE WHETHER A POINT IS INSIDE A POLYGON            
C                                                                       
C        USAGE                                                          
C           CALL PNPOLY (PX, PY, XX, YY, N, INOUT )                     
C                                                                       
C        DESCRIPTION OF THE PARAMETERS                                  
C           PX      - X-COORDINATE OF POINT IN QUESTION.                
C           PY      - Y-COORDINATE OF POINT IN QUESTION.                
C           XX      - N LONG VECTOR CONTAINING X-COORDINATES OF         
C                     VERTICES OF POLYGON.                              
C           YY      - N LONG VECTOR CONTAING Y-COORDINATES OF           
C                     VERTICES OF POLYGON.                              
C           N       - NUMBER OF VERTICES IN THE POLYGON.                
C           INOUT   - THE SIGNAL RETURNED:                              
C                     -1 IF THE POINT IS OUTSIDE OF THE POLYGON,        
C                      0 IF THE POINT IS ON AN EDGE OR AT A VERTEX,     
C                      1 IF THE POINT IS INSIDE OF THE POLYGON.         
C                                                                       
C        REMARKS                                                        
C           THE VERTICES MAY BE LISTED CLOCKWISE OR ANTICLOCKWISE.      
C           THE FIRST MAY OPTIONALLY BE REPEATED, IF SO N MAY           
C           OPTIONALLY BE INCREASED BY 1.                               
C           THE INPUT POLYGON MAY BE A COMPOUND POLYGON CONSISTING      
C           OF SEVERAL SEPARATE SUBPOLYGONS. IF SO, THE FIRST VERTEX    
C           OF EACH SUBPOLYGON MUST BE REPEATED, AND WHEN CALCULATING   
C           N, THESE FIRST VERTICES MUST BE COUNTED TWICE.              
C           INOUT IS THE ONLY PARAMETER WHOSE VALUE IS CHANGED.         
C           THE SIZE OF THE ARRAYS MUST BE INCREASED IF N > MAXDIM      
C           WRITTEN BY RANDOLPH FRANKLIN, UNIVERSITY OF OTTAWA, 7/70.   
C                                                                       
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  
C           NONE                                                        
C                                                                       
C        METHOD                                                         
C           A VERTICAL LINE IS DRAWN THRU THE POINT IN QUESTION. IF IT  
C           CROSSES THE POLYGON AN ODD NUMBER OF TIMES, THEN THE        
C           POINT IS INSIDE OF THE POLYGON.                             
C                                                                       
C     ..................................................................
C                                                                       
      SUBROUTINE pnpoly(PX,PY,XX,YY,N,INOUT)
      IMPLICIT NONE
      INTEGER MAXDIM
      PARAMETER (MAXDIM=1000)
      REAL*8 X(MAXDIM),Y(MAXDIM),XX(N),YY(N),PX,PY,var
      INTEGER N,INOUT,I,J
      LOGICAL MX,MY,NX,NY

C Check that number of points in clipping path is less than MAXDIM
      IF (N.GT.MAXDIM) THEN
          WRITE(0,7) N
7         FORMAT('WARNING:',I5,' TOO GREAT FOR THIS VERSION OF PNPOLY')
          RETURN
      ENDIF

C Compute the position of the points on the clipping path (xx,yy)
C relative to the point of interest (px,py)
      DO 1 I=1,N
          X(I)=XX(I)-PX
          Y(I)=YY(I)-PY
    1 CONTINUE

      INOUT=-1
      DO 2 I=1,N
          J=1+MOD(I,N)
          MX=X(I).GE.0.0 ! Is the point on the clipping path east of the point of interest?
          NX=X(J).GE.0.0 ! Is the next point on the clipping path east of the point of interest?
          MY=Y(I).GE.0.0 ! Is the point on the clipping path north of the point of interest?
          NY=Y(J).GE.0.0 ! Is the next point on the clipping path north of the point of interest?
          IF(.NOT.((MY.OR.NY).AND.(MX.OR.NX)).OR.(MX.AND.NX)) THEN
              GOTO 2
          ENDIF
          IF(.NOT.(MY.AND.NY.AND.(MX.OR.NX).AND..NOT.(MX.AND.NX))) THEN
              GOTO 3
          ENDIF
          INOUT=-INOUT
          GOTO 2
3         var = (Y(I)*X(J)-X(I)*Y(J))/(X(J)-X(I))
          if (var.lt.0.0d0) then
              goto 2
          elseif (dabs(var).lt.1.0d-8) then
              goto 4
          else
              goto 5
          endif
C3         IF ((Y(I)*X(J)-X(I)*Y(J))/(X(J)-X(I))) 2,4,5
4         INOUT=0
          RETURN
5         INOUT=-INOUT
2     CONTINUE
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE gcmdln(x1,x2,nx,dx,y1,y2,ny,dy,z,ofile,p,ref,xsec,xz,
     1                  clip,datain,expflg)
      IMPLICIT none
      CHARACTER*40 tag,ofile
      CHARACTER*200 clip
      REAL*8 x1,x2,y1,y2,z,dx,dy,ref(5)
      INTEGER narg,i,j,nx,ny,p,xsec,xz,datain,expflg
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
      clip = 'none'
      datain = 0
      expflg = 0
      do 100 i = 1,5
          ref(i) = -1.0d0
  100 continue
      narg = iargc()
      if (narg.eq.0) then
          call usage('')
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
      elseif (tag(1:4).eq.'-in?') then
          datain = 1
      elseif (tag(1:5).eq.'-clip') then
          i = i + 1
          call getarg(i,clip)
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
      elseif (tag(1:4).eq.'-exp') then
          i = i + 1
          call getarg(i,tag)
          if(tag.eq.'x') then
              expflg = expflg + 1
          elseif (tag.eq.'y') then
              expflg = expflg + 10
          else
              call usage('!! Error: no -exp option'//trim(tag))
          endif
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
      elseif (tag(1:2).eq.'-R') then
          call getarg(i,tag)
          j = index(tag,'-R')
          tag(j:j+1) = ''
          j = index(tag,'/')
          tag(j:j) = ' '
          j = index(tag,'/')
          tag(j:j) = ' '
          j = index(tag,'/')
          tag(j:j) = ' '
          read(tag,*) x1,x2,y1,y2
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
     1 'Usage: grid -x X1 X2 -nx NX|-dx DX ',
     2            '[-y Y1 Y2 -ny NY|-dy DY] [-RX1/X2/Y1/Y2]'
      write(*,*)
     1 '            [-z Z|-dip X0 Y0 Z0 STR ',
     2              'DIP|-xsec X0 Y0 AZ [-xz]] [-exp x|y]'
      write(*,*)
     1 '            [-clip XYFILE] [-o OFILE] [-p] [-h]'
      write(*,*)
      write(*,*)
     1 '-x X1 X2               First column (x) limits'
      write(*,*)
     1 '-nx NX                 Number of x grid points'
      write(*,*)
     1 '-dx DX                 Increment in x direction (overrides -nx)'
      write(*,*)
     1 '-y Y1 Y2               Second column (y) limits'
      write(*,*)
     1 '-ny NY                 Number of y grid points'
      write(*,*)
     1 '-dy DY                 Increment in y direction (overrides -ny)'
      write(*,*)
     1 '-RX1/X2/Y1/Y2          Define limits with GMT limits format'
      write(*,*)
     1 '-z Z                   Third column (z) value'
      write(*,*)
     1 '-dip X0 Y0 Z0 STR DIP  Put grid onto plane with STR/DIP ',
     2                'and containing (X0,Y0,Z0)'
      write(*,*)
     1 '-xsec X0 Y0 AZ         Create a vertical cross section ',
     2                 'through (X0,Y0) with strike AZ'
      write(*,*)
     1 '-xz                    Print cross section x-z instead of ',
     2                            'lon lat z'
      write(*,*)
     1 '-clip XYFILE           Clip points to lie in interior of ',
     2                            'XYFILE (outside points are -12345)'
      write(*,*)
     1 '-in?                   Check if input points (from standard ',
     2                          'input) are in polygon (requires -clip',
     3                          ' to be set)'
      write(*,*)
     1 '-o OFILE               Output to file (default prints to ',
     2              'standard output)'
      write(*,*)
     1 '-exp x|y               Make grid increment x|y exponential'
      write(*,*)
     1 '                           nx/ny => equal fraction spacing'
      write(*,*)
     1 '                           dx/dy => fractional increment'
      write (*,*)
     1 '-p                     Print results to standard output ',
     2                               '(overrides -o)'
      write (*,*)
     1 '-h/-?                  Online help (this screen)'
      write (*,*)
      STOP
      END
