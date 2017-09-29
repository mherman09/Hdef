      PROGRAM pole2vel

      IMPLICIT none
      INTEGER i
      CHARACTER*80 ifile,ofile
      INTEGER p
      REAL*8 pole(3),pxyz(3)
      INTEGER STAMAX,nsta
      PARAMETER (STAMAX=1000)
      REAL*8 sta(STAMAX,2)
      REAL*8 vel(STAMAX,3)
      INTEGER fun

C Get command line arguments
      call gcmdln(ifile,ofile,pole,p)
C      print *,'IFILE ',ifile
C      print *,'OFILE ',ofile
C      print *,'POLE ',pole(1),pole(2),pole(3)
C      print *,'PRINT ',p

C Convert Euler pole to XYZ components
      call lalomag2xyz(pxyz,pole)
C      print *,'POLE ',pxyz(1),pxyz(2),pxyz(3)

C Read station locations
      call readsta(ifile,sta,nsta)
C      print *,'NSTA ',nsta
C      do 101 i = 1,nsta
C          print *,i,sta(i,1),sta(i,2)
C  101 continue

C Compute NEZ velocity at stations
      call calcvel(sta,nsta,pxyz,vel)
C      do 102 i = 1,nsta
C          print *,i,vel(i,1),vel(i,2),vel(i,3)
C  102 continue

      if (p.eq.1.or.ofile.eq.'none') then
          fun = 6
      else
          fun = 11
          open(unit=fun,file=ofile,status='unknown')
      endif
      do 110 i = 1,nsta
          write(fun,9999) sta(i,1),sta(i,2),vel(i,1),vel(i,2),vel(i,3)
  110 continue
 9999 format(2F10.4,3F10.2)

      END

C----------------------------------------------------------------------C

      SUBROUTINE readsta(ifile,sta,nsta)
      IMPLICIT none
      INTEGER STAMAX,nsta
      PARAMETER (STAMAX=1000)
      CHARACTER*80 ifile
      REAL*8 sta(STAMAX,2)
      INTEGER funit
      if (ifile.eq.'stdin') then
          funit = 5
      else
          funit = 21
          open(unit=funit,file=ifile,status='old')
      endif
      nsta = 1
  210 read(funit,*,end=211) sta(nsta,1),sta(nsta,2)
          nsta = nsta + 1
          goto 210
  211 continue
      nsta = nsta - 1
      if(funit.eq.21)close(funit)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE lalomag2xyz(pxyz,pole)
      IMPLICIT none
      REAL*8 pi,d2r
      PARAMETER (pi=4.0d0*atan(1.0d0),d2r=pi/180.0d0)
      REAL*8 pole(3),pxyz(3)
      pxyz(1) = pole(3)*dcos(pole(2)*d2r)*dcos(pole(1)*d2r)
      pxyz(2) = pole(3)*dcos(pole(2)*d2r)*dsin(pole(1)*d2r)
      pxyz(3) = pole(3)*dsin(pole(2)*d2r)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE calcvel(sta,nsta,pxyz,vel)
      IMPLICIT none
      REAL*8 pi,d2r
      PARAMETER (pi=4.0d0*atan(1.0d0),d2r=pi/180.0d0)
      INTEGER STAMAX,nsta
      PARAMETER (STAMAX=1000)
      REAL*8 sta(STAMAX,2),lo,la,stavec(3)
      REAL*8 pxyz(3),v(3),r(3),vel(STAMAX,3)
      INTEGER i
      do 301 i = 1,nsta
          stavec(1) = sta(i,1)
          stavec(2) = sta(i,2)
          stavec(3) = 1.0d0
          call lalomag2xyz(r,stavec)
          call cross(v,pxyz,r)
          lo = sta(i,1)*d2r
          la = sta(i,2)*d2r
          vel(i,2) = -(v(1)*dsin(la)*dcos(lo)+v(2)*dsin(la)*dsin(lo)
     1                                                  -v(3)*dcos(la))
          vel(i,1) =  -v(1)         *dsin(lo)+v(2)         *dcos(lo)
          vel(i,3) =   v(1)*dcos(la)*dcos(lo)+v(2)*dcos(la)*dsin(lo)
     1                                                   +v(3)*dsin(la)
          vel(i,1) = vel(i,1)*111.19d0
          vel(i,2) = vel(i,2)*111.19d0
  301 continue
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
     1 'Usage: pole2vel -pole LON LAT WVEL -sta STAFILE -vel VFILE -p',
     2                ' -h/-?'
      write(*,*)
     1 '  -pole LON LAT WVEL   Euler pole location and angular ',
     2                                               'velocity (deg/Ma)'
      write(*,*)
     1 '  -sta STAFILE         Locations to compute velocity (',
     2                                                  'takes "stdin")'
      write(*,*)
     1 '  -vel VFILE           Velocity (ENZ) at locations in STAFILE'
      write(*,*)
     1 '  -p                   Print to standard out instead of VFILE'
      write(*,*)
     1 '  -h/-?                Online help (this screen)'
      write(*,*)
      STOP
      END

C----------------------------------------------------------------------C

      SUBROUTINE gcmdln(ifile,ofile,pole,p)
      IMPLICIT none
      CHARACTER*80 tag,ifile,ofile
      INTEGER narg,i,p
      REAL*8 pole(3)
      ifile = 'none'
      ofile = 'none'
      pole(1) = 0.0d0
      pole(2) = 0.0d0
      pole(3) = 0.0d0
      p = 0
      narg = iargc()
      if (narg.eq.0) call usage('!! Error: no command line arguments')
      i = 0
  101 i = i + 1
      if (i.gt.narg) goto 102
          call getarg(i,tag)
          if (tag(1:5).eq.'-pole') then
              i = i + 1
              call getarg(i,tag)
              read(tag,'(BN,F10.0)')pole(1)
              i = i + 1
              call getarg(i,tag)
              read(tag,'(BN,F10.0)')pole(2)
              i = i + 1
              call getarg(i,tag)
              read(tag,'(BN,F10.0)')pole(3)
          elseif (tag(1:4).eq.'-sta') then
              i = i + 1
              call getarg(i,ifile)
          elseif (tag(1:4).eq.'-vel') then
              i = i + 1
              call getarg(i,ofile)
          elseif (tag(1:2).eq.'-p') then
              p = 1
          elseif (tag(1:2).eq.'-h'.or.tag(1:2).eq.'-?') then
              call usage('')
          else
              call usage('!! Error: no option'//tag)
          endif
      goto 101
  102 continue
      RETURN
      END

