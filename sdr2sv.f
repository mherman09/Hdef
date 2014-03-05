      PROGRAM sdr2sv
C----
C Given a list of strike, dip, and rake double couples, calculate the 
C slip vectors for the two nodal planes.
C----
      IMPLICIT none
      CHARACTER*30 ifile,ofile
      REAL*8 strin,dipin,rakin
      REAL*8 svaz1,svpl1,svaz2,svpl2
      INTEGER user

      call gcmdln(ifile,ofile,user)

      if (user.eq.1) then
          print *,'Enter strike, dip, and rake:'
          read *, strin,dipin,rakin
          call sdr2svgrit(svaz1,svpl1,svaz2,svpl2,strin,dipin,rakin)
          write (*,8887)
          write (*,8888) svaz1,svpl1
          write (*,8889) svaz2,svpl2
      else
          open(unit=11,file=ifile,status='old')
          open(unit=12,file=ofile,status='unknown')
  101     read (11,*,end=102) strin,dipin,rakin
              call sdr2svgrit(svaz1,svpl1,svaz2,svpl2,strin,dipin,rakin)
              write (12,9999) svaz1,svpl1,svaz2,svpl2
              goto 101
  102     continue
      endif

 8887 format(14X,'      AZIMUTH(DEG)  INCLINATION(DEG)')
 8888 format('Slip vector 1:',2(10X,F8.3))
 8889 format('Slip vector 2:',2(10X,F8.3))
 9999 format(4F14.6)

      END

C======================================================================C

      SUBROUTINE sdr2svgrit(svaz1,svpl1,svaz2,svpl2,strin,dipin,rakin)
C----
C Calculations
C----
      IMPLICIT none
      REAL*8 pi,d2r,r2d
      PARAMETER (pi=4.0d0*datan(1.0d0),d2r=pi/1.8d2,r2d=1.8d2/pi)
      REAL*8 strin,dipin,rakin,str1,dip1,rak1,str2,dip2,rak2
      REAL*8 x1,y1,z1,TH,x2,y2,z2,N1,E1,N2,E2
      REAL*8 svaz1,svpl1,svaz2,svpl2

C Convert to radians
      str1 = strin*d2r
      dip1 = dipin*d2r
      rak1 = rakin*d2r

C----
C CALCULATE THE STRIKE, DIP, RAKE OF OTHER NODAL PLANE
C----
C Direction of motion of hanging wall on input plane.
C x is in the direction of strike, y is in the updip direction,
C z is vertical up
      X1 = dcos(rak1)
      Y1 = dsin(rak1)*dcos(dip1)
      Z1 = dsin(rak1)*dsin(dip1)

C Horizontal CCW angle from strike
      TH = datan2(Y1,X1)

C Vector (X1,Y1,Z1) is the normal to the second plane
      if (Z1.ge.0.0d0) then
          str2 = str1 - TH - 0.5*pi
          dip2 = 0.5d0*pi - dataN2(Z1,dsqrt(X1*X1+Y1*Y1))
          rak2 = dacos(dsin(str2-str1)*dcos(0.5d0*pi-dip1))
      else
          str2 = str1 - TH + 0.5d0*pi
          dip2 = 0.5d0*pi + datan2(Z1,dsqrt(X1*X1+Y1*Y1))
          rak2 = -pi+dacos(dsin(str2-str1)*dcos(0.5d0*pi-dip1))
      endif

      call principalval(str1)
      call principalval(str2)

      X2 = dcos(rak2)
      Y2 = dsin(rak2)*dcos(dip2)
      Z2 = dsin(rak2)*dsin(dip2)

C----
C CALCULATE SLIP VECTORS FROM STRIKE, DIP, RAKE
C----
      N1 = X1*dcos(str1) + Y1*dsin(str1)
      E1 = X1*dsin(str1) - Y1*dcos(str1)

      N2 = X2*dcos(str2) + Y2*dsin(str2)
      E2 = X2*dsin(str2) - Y2*dcos(str2)

      svaz1 = datan2(E1,N1)
      svpl1 = datan2(Z1,dsqrt(E1*E1+N1*N1))

      svaz2 = datan2(E2,N2)
      svpl2 = datan2(Z2,dsqrt(E2*E2+N2*N2))

      svaz1 = svaz1*r2d
      svpl1 = svpl1*r2d
      svaz2 = svaz2*r2d
      svpl2 = svpl2*r2d

      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE principalval(arg)

      IMPLICIT none
      REAL*8 arg,pi,tpi
      PARAMETER (pi=4.0d0*datan(1.0d0),tpi=2.0d0*pi)
      
   11 if (arg.lt.0.0d0) then
          arg = arg + tpi
          goto 11
      endif
   12 if (arg.gt.tpi) then
          arg = arg - tpi
          goto 12
      endif

      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE gcmdln(ifile,ofile,user)
      IMPLICIT none
      CHARACTER*30 ifile,ofile,tag
      INTEGER i,narg,user

      user = 0
      ifile = 'sdr2sv.in'
      ofile = 'sdr2sv.out'

      narg = iargc()
      if (narg.eq.0) call usage()

      i = 0
   11 i = i + 1
      if (i.gt.narg) goto 99
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
          goto 11
   99 continue

      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE usage()
      IMPLICIT none

      write(*,*)
     1 'Usage: sdr2sv -f [IFILE] -o [OFILE] -d -u -h/-?'
      write(*,*)
     1 '  -f [IFILE] (Default sdr2sv.in) name of input file'
      write(*,*)
     1 '  -o [OFILE] (Default sdr2sv.out) name of output file'
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
     1 '  sdr2sv calculates the azimuth and inclination (measured up',
     2   ' from horizontal)'
      write(*,*)
     1 '      of slip vectors from double couple source parameters'
      write (*,*) ''
      write (*,*)
     1 '    Input file'
      write (*,*)
     1 '      strike dip rake'
      write (*,*)
     1 '         :  :'
      write (*,*) ''
      write (*,*)
     1 '    Output file'
      write (*,*)
     1 '      az1 incl1 az2 incl2 (in degrees)'
      write (*,*)
     1 '         :  :'
      write (*,*) ''

      STOP
      END
