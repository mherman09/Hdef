      PROGRAM main
C----
C Compute relative plate velocities given a pair of plates or Euler
C pole (if plates are given, use poles from MORVEL; Demets et al., 2010)
C----
      IMPLICIT none
      REAL*8 pi,d2r,r2d
      PARAMETER (pi=4.0d0*datan(1.0d0),d2r=pi/1.8d2,r2d=1.8d2/pi)
      CHARACTER*80 ifile,ofile,polec
      INTEGER i,isave
      REAL*8 pole(3),stlo,stla,x,y,z,r(3),v(3),lo,la,ve,vn,vz

C----
C Get command line arguments
C----
      call gcmdln(ifile,ofile,polec)

C----
C Get pole
C----
      i = index(polec,'/')
      if (i.gt.0) then
          polec(i:i) = ' '
          isave = i
      endif
      i = index(polec,'/')
      if (i.gt.0) then
          polec(i:i) = ' '
          i = -i
      else
          i = isave
      endif
      if (i.gt.0) then
          call morvel(pole,polec)
      elseif (i.lt.0) then
          read(polec,*) pole(1),pole(2),pole(3)
      else
          call usage('')
      endif
      x = pole(3)*dcos(pole(2)*d2r)*dcos(pole(1)*d2r)
      y = pole(3)*dcos(pole(2)*d2r)*dsin(pole(1)*d2r)
      z = pole(3)*dsin(pole(2)*d2r)
      pole(1) = x
      pole(2) = y
      pole(3) = z

C----
C Calculate predicted velocities
C----
      if (ifile.ne.'none') open(unit=11,file=ifile,status='old')
  102 if (ifile.eq.'none') then
          read(*,*,end=101) stlo,stla
      else
          read(11,*,end=101) stlo,stla
      endif
      la = stla*d2r
      lo = stlo*d2r
      r(1) = dcos(stla*d2r)*dcos(stlo*d2r)
      r(2) = dcos(stla*d2r)*dsin(stlo*d2r)
      r(3) = dsin(stla*d2r)
      call cross(v,pole,r)
      vn = -(v(1)*dsin(la)*dcos(lo)+v(2)*dsin(la)*dsin(lo)
     1                                                   -v(3)*dcos(la))
      ve =  -v(1)         *dsin(lo)+v(2)         *dcos(lo)
      vz =   v(1)*dcos(la)*dcos(lo)+v(2)*dcos(la)*dsin(lo)
     1                                                    +v(3)*dsin(la)
      vn = vn*111.19d0
      ve = ve*111.19d0
      print *,stlo,stla,ve,vn
      goto 102
  101 continue
      END

C----------------------------------------------------------------------C

      SUBROUTINE morvel(pole,polec)
      IMPLICIT none
      CHARACTER*80 polec,plate(2)
      REAL*8 pole(3),x(2,3)
      INTEGER i
      read(polec,*) plate(1),plate(2)
      do 101 i = 1,2
          if (plate(i).eq.'AM') then
              x(i,1) = 0.04829
              x(i,2) = -0.37592
              x(i,3) = 0.84817
          elseif (plate(i).eq.'AN') then
              x(i,1) = 0.07204
              x(i,2) = -0.35489
              x(i,3) = 0.81026
          elseif (plate(i).eq.'AR') then
              x(i,1) = 0.48451
              x(i,2) = -0.31742
              x(i,3) = 1.00388
          elseif (plate(i).eq.'AU') then
              x(i,1) = 0.53463
              x(i,2) = 0.05931
              x(i,3) = 0.93478
          elseif (plate(i).eq.'CA') then
              x(i,1) = 0.11032
              x(i,2) = -0.49683
              x(i,3) = 0.74786
          elseif (plate(i).eq.'CO') then
              x(i,1) = -0.48084
              x(i,2) = -1.14521
              x(i,3) = 1.12543
          elseif (plate(i).eq.'CP') then
              x(i,1) = 0.52040
              x(i,2) = -0.09291
              x(i,3) = 1.00856
          elseif (plate(i).eq.'EU') then
              x(i,1) = 0.07941
              x(i,2) = -0.40365
              x(i,3) = 0.75055
          elseif (plate(i).eq.'IN') then
              x(i,1) = 0.46728
              x(i,2) = -0.28306
              x(i,3) = 1.00170
          elseif (plate(i).eq.'JF') then
              x(i,1) = 0.49382
              x(i,2) = 0.38353
              x(i,3) = -0.00681
          elseif (plate(i).eq.'LW') then
              x(i,1) = 0.18270
              x(i,2) = -0.42829
              x(i,3) = 0.80758
          elseif (plate(i).eq.'NA') then
              x(i,1) = 0.15483
              x(i,2) = -0.46836
              x(i,3) = 0.56521
          elseif (plate(i).eq.'NB') then
              x(i,1) = 0.19329
              x(i,2) = -0.44604
              x(i,3) = 0.79881
          elseif (plate(i).eq.'MQ') then
              x(i,1) = 0.85484
              x(i,2) = -0.11982
              x(i,3) = 1.44874
          elseif (plate(i).eq.'NZ') then
              x(i,1) = 0.02874
              x(i,2) = -0.73542
              x(i,3) = 1.08528
          elseif (plate(i).eq.'PS') then
              x(i,1) = 0.66046
              x(i,2) = -0.59193
              x(i,3) = -0.07185
          elseif (plate(i).eq.'RI') then
              x(i,1) = -1.14352
              x(i,2) = -4.32660
              x(i,3) = 2.15245
          elseif (plate(i).eq.'SA') then
              x(i,1) = 0.08198
              x(i,2) = -0.35590
              x(i,3) = 0.54095
          elseif (plate(i).eq.'SC') then
              x(i,1) = 0.08341
              x(i,2) = -0.39306
              x(i,3) = 0.63896
          elseif (plate(i).eq.'SM') then
              x(i,1) = 0.14189
              x(i,2) = -0.48048
              x(i,3) = 0.84261
          elseif (plate(i).eq.'SR') then
              x(i,1) = 0.08813
              x(i,2) = -0.34742
              x(i,3) = 0.52526
          elseif (plate(i).eq.'SU') then
              x(i,1) = 0.10210
              x(i,2) = -0.47858
              x(i,3) = 0.84111
          elseif (plate(i).eq.'SW') then
              x(i,1) = 1.06493
              x(i,2) = -0.97101
              x(i,3) = -0.09671
          elseif (plate(i).eq.'YZ') then
              x(i,1) = 0.05326
              x(i,2) = -0.39837
              x(i,3) = 0.88007
          elseif (plate(i).eq.'PA') then
              x(i,1) = 0.0
              x(i,2) = 0.0
              x(i,3) = 0.0
          else
              write(*,*) '!! Error: no plate named '//trim(plate(i))
              write(*,*) 'AM: Amur'
              write(*,*) 'AN: Antarctica'
              write(*,*) 'AR: Arabia'
              write(*,*) 'AU: Australia'
              write(*,*) 'CA: Caribbean'
              write(*,*) 'CO: Cocos'
              write(*,*) 'CP: Capricorn'
              write(*,*) 'EU: Eurasia'
              write(*,*) 'IN: India'
              write(*,*) 'JF: Juan de Fuca'
              write(*,*) 'LW: Lwandle'
              write(*,*) 'MQ: Macquarie'
              write(*,*) 'NA: North America'
              write(*,*) 'NB: Nubia'
              write(*,*) 'NZ: Nazca'
              write(*,*) 'PA: Pacific'
              write(*,*) 'PS: Philippine Sea'
              write(*,*) 'RI: Rivera'
              write(*,*) 'SA: South America'
              write(*,*) 'SC: Scotia'
              write(*,*) 'SM: Somalia'
              write(*,*) 'SR: Sur'
              write(*,*) 'SU: Sunda'
              write(*,*) 'SW: Sandwich'
              write(*,*) 'YZ: Yangtze'
              call usage('')
          endif
  101 continue
      x(1,1) = x(2,1) - x(1,1)
      x(1,2) = x(2,2) - x(1,2)
      x(1,3) = x(2,3) - x(1,3)
      pole(3) = dsqrt(x(1,1)*x(1,1)+x(1,2)*x(1,2)+x(1,3)*x(1,3))
      pole(2) = dacos(x(1,3)/pole(3))
      pole(2) = 90.0d0 - pole(2)*1.8d2/3.14159265d0
      pole(1) = datan2(x(1,2),x(1,1))*1.8d2/3.14159265d0
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE lalomag2xyz(x,y,z,lat,lon,mag)
      IMPLICIT none
      REAL*8 pi,d2r
      PARAMETER (pi=4.0d0*atan(1.0d0),d2r=pi/180.0d0)
      REAL*8 lat,lon,mag,x,y,z
      x = mag*cos(lat*d2r)*cos(lon*d2r)
      y = mag*cos(lat*d2r)*sin(lon*d2r)
      z = mag*sin(lat*d2r)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE cross(ovect,v1,v2)
      IMPLICIT none
      REAL*8 ovect(3),v1(3),v2(3)
      ovect(1) = v1(2)*v2(3) - v1(3)*v2(2)
      ovect(2) = v1(3)*v2(1) - v1(1)*v2(3)
      ovect(3) = v1(1)*v2(2) - v1(2)*v2(1)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE gcmdln(ifile,ofile,pole)
      IMPLICIT none
      CHARACTER*80 tag,ifile,ofile,pole
      INTEGER i,narg
      ifile = 'none'
      ofile = 'none'
      pole  = 'none'
      narg = iargc()
      if (narg.eq.0) call usage('')
      i = 0
  101 i = i + 1
      if (i.gt.narg) goto 102
          call getarg(i,tag)
          if (tag(1:5).eq.'-pole'.or.tag(1:6).eq.'-plate') then
              i = i + 1
              call getarg(i,pole)
          elseif (tag(1:2).eq.'-h') then
              call usage('')
          elseif (tag(1:2).eq.'-f') then
              i = i + 1
              call getarg(i,ifile)
          elseif (tag(1:2).eq.'-o') then
              i = i + 1
              call getarg(i,ofile)
          endif
      goto 101
  102 continue
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE usage(str)
      IMPLICIT none
      INTEGER d
      CHARACTER str*(*)
      d = 0
      if (str.eq.'d') d = 1
      write(*,*)
     1 'Usage: platemotion -plate P1/P2 | -pole lon/lat/vel',
     2                       ' [-f IFILE] [-o OFILE] [-h]'
      write(*,*)
     1 '  -plate P1/P2      Define plates (P2 moving w.r.t. fixed P1)'
      write(*,*)
     1 '  -pole LON/LAT/VEL Define Euler pole (velocity in deg/Ma)'
      write(*,*)
     1 '  -f IFILE      Input file (default: run with stdin)'
      write(*,*)
     1 '  -o OFILE      Output file (default: print to stdout)'
      write(*,*)
     1 '  -h            Short online help'
      write(*,*)
      STOP
      END

