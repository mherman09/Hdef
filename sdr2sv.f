      PROGRAM sdr2sv
C----
C Given a list of strike, dip, and rake double couples, calculate the 
C slip vectors for the two nodal planes.
C----
      IMPLICIT none
      CHARACTER*30 ifile,ofile
      REAL strin,dipin,rakin
      REAL svaz1,svpl1,svaz2,svpl2

      call gcmdln(ifile,ofile)


      if (ifile.eq.'none') then
          print *,'To see options, use -h flag'
          print *,'Enter strike, dip, and rake:'
          read *, strin,dipin,rakin
          call sdr2svgrit(svaz1,svpl1,svaz2,svpl2,strin,dipin,rakin)
          print *,'                Azimuth (deg)   Inclination (deg)'
          print *,'Slip vector 1:',svaz1,svpl1
          print *,'Slip vector 2:',svaz2,svpl2
      else
          open(unit=11,file=ifile,status='old')
          open(unit=12,file=ofile,status='unknown')
  101     read (11,*,end=102) strin,dipin,rakin
              call sdr2svgrit(svaz1,svpl1,svaz2,svpl2,strin,dipin,rakin)
              write (12,*) svaz1,svpl1,svaz2,svpl2
              goto 101
  102     continue
      endif

      END

C======================================================================C

      SUBROUTINE sdr2svgrit(svaz1,svpl1,svaz2,svpl2,strin,dipin,rakin)
C----
C Calculations
C----
      IMPLICIT none
      REAL pi,d2r,r2d
      PARAMETER (pi=3.14159265,d2r=pi/180.,r2d=180./pi)
      REAL strin,dipin,rakin,str1,dip1,rak1,str2,dip2,rak2
      REAL x1,y1,z1,TH,x2,y2,z2,N1,E1,N2,E2
      REAL svaz1,svpl1,svaz2,svpl2

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
      X1 = cos(rak1)
      Y1 = sin(rak1)*cos(dip1)
      Z1 = sin(rak1)*sin(dip1)

C Horizontal CCW angle from strike
      TH = ataN2(Y1,X1)

C Vector (X1,Y1,Z1) is the normal to the second plane
      if (Z1.ge.0.0) then
          str2 = str1 - TH - 0.5*pi
          dip2 = 0.5*pi - ataN2(Z1,sqrt(X1*X1+Y1*Y1))
          rak2 = acos(sin(str2-str1)*cos(0.5*pi-dip1))
      else
          str2 = str1 - TH + 0.5*pi
          dip2 = 0.5*pi + ataN2(Z1,sqrt(X1*X1+Y1*Y1))
          rak2 = -pi+acos(sin(str2-str1)*cos(0.5*pi-dip1))
      endif

      call principalval(str1)
      call principalval(str2)

      X2 = cos(rak2)
      Y2 = sin(rak2)*cos(dip2)
      Z2 = sin(rak2)*sin(dip2)

C----
C CALCULATE SLIP VECTORS FROM STRIKE, DIP, RAKE
C----
      N1 = X1*cos(str1) + Y1*sin(str1)
      E1 = X1*sin(str1) - Y1*cos(str1)

      N2 = X2*cos(str2) + Y2*sin(str2)
      E2 = X2*sin(str2) - Y2*cos(str2)

      svaz1 = ataN2(E1,N1)
      svpl1 = ataN2(Z1,sqrt(E1*E1+N1*N1))

      svaz2 = ataN2(E2,N2)
      svpl2 = ataN2(Z2,sqrt(E2*E2+N2*N2))

      svaz1 = svaz1*r2d
      svpl1 = svpl1*r2d
      svaz2 = svaz2*r2d
      svpl2 = svpl2*r2d

      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE principalval(arg)

      IMPLICIT none
      REAL arg
      
   11 if (arg.lt.0.0) then
          arg = arg + 2.0*3.14159265
          goto 11
      endif
   12 if (arg.gt.2.0*3.14159265) then
          arg = arg - 2.0*3.14159265
          goto 12
      endif

      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE gcmdln(ifile,ofile)
      IMPLICIT none
      CHARACTER*30 ifile,ofile,tag
      INTEGER i,narg

      ifile = 'none'
      ofile = 'sdr2sv.out'

      narg = iargc()
      i = 0
   11 i = i + 1
      if (i.gt.narg) goto 99
          call getarg(i,tag)
          if (tag(1:2).eq.'-i') then
              i = i + 1
              call getarg(i,ifile)
          elseif (tag(1:2).eq.'-o') then
              i = i + 1
              call getarg(i,ofile)
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
     1 'Usage: sdr2sv -f [IFILE] -o [OFILE] -h/-?'
      write(*,*)
     1 '  -f [IFILE] (Default none) name of input file'
      write(*,*)
     1 '                            if no file name given, fault ',
     2                            'parameters prompted for manual input'
      write(*,*)
     1 '  -o [OFILE] (Default sdr2sv.out) name of output file'
      write (*,*)
     1 '  -h/-?        help'
      write (*,*) ''
      write (*,*)
     1 '  sdr2sv calculates the azimuth and inclination of slip ',
     2   'vectors from a double couple source'
      write (*,*) ''
      write (*,*)
     1 '    input file'
      write (*,*)
     1 '        strike dip rake'
      write (*,*)
     1 '         :  :'
      write (*,*)
     1 '    output file'
      write (*,*)
     1 '        az1 incl1 az2 incl2 (in degrees)'
      write (*,*)
     1 '         :  :'
      write (*,*) ''

      STOP
      END
