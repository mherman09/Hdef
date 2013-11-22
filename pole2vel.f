      PROGRAM pole2vel

      IMPLICIT none
      REAL*8 pi,d2r
      PARAMETER (pi=4.0d0*atan(1.0d0),d2r=pi/180.0d0)
      CHARACTER*80 ifile,ofile
      LOGICAL ex
      INTEGER pfmt
      CHARACTER*8 stnm
      REAL*8 plat,plon,pmag,p(3)
      REAL*8 stlo,stla,la,lo,uN,uE,r(3)
      REAL*8 v(3),vN,vE,vZ

C----
C Get command line arguments
C Check for existence of input file
C----
      call gcmdln(ifile,ofile,pfmt)
      inquire(file=ifile,exist=ex)
      if (.not.ex) call usage('INPUT FILE NOT PRESENT; CHECK SPELLING')

C----
C Read input file and calculate predicted velocities
C----
      open(unit=11,file=ifile,status='old')
      if (ofile.eq.'none') then
          ex = .false.
      else
          ex = .true.
          open(unit=12,file=ofile,status='unknown')
      endif

      if (pfmt.eq.1) then
          read(11,*) plat,plon,pmag
          call lalomag2xyz(p(1),p(2),p(3),plat,plon,pmag)
      else
          read(11,*) p(1),p(2),p(3)
      endif

  101 read(11,*,end=102) stnm,stla,stlo,uN,uE
          call lalomag2xyz(r(1),r(2),r(3),stla,stlo,1.0d0)
          call cross(v,p,r)
          la = stla*d2r
          lo = stlo*d2r
          vN = -(v(1)*sin(la)*cos(lo)+v(2)*sin(la)*sin(lo)-v(3)*cos(la))
          vE =  -v(1)        *sin(lo)+v(2)        *cos(lo)
          vZ =   v(1)*cos(la)*cos(lo)+v(2)*cos(la)*sin(lo)+v(3)*sin(la)
          vN = vN*111.19d0
          vE = vE*111.19d0
          if (ex) then
              write(12,9999) stnm,stla,stlo,uN,uE,vN,vE
          else
              write(6,9999) stnm,stla,stlo,uN,uE,vN,vE
          endif
          goto 101
  102 continue

 9999 format(A8,2F10.3,4F8.2)
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

      SUBROUTINE usage(text)
      IMPLICIT none
      CHARACTER text*(*)

      write(*,*)
      write(*,*) text
      write(*,*)
      write(*,*)
     1 'Usage: pole2vel -f [IFILE] -o [OFILE] -p [POLE_FORMAT]'
      write(*,*)
     1 '    -f [IFILE] (required to run) Name of input file'
      write(*,*)
     1 '    -o [OFILE] (default none) Name of output file (none =>',
     2                                         ' print to terminal)'
      write(*,*)
     1 '    -p [POLE_FORMAT] (default 1) rotation pole format'
      write(*,*)
     1 '        1: Latitude, Longitude, Angular Velocity'
      write(*,*)
     1 '        2: OmegaX, OmegaY, OmegaZ'
      write(*,*)
     1 '    -h/-? Print help'
      write(*,*)
      write(*,*)
     1 '  pole2vel takes an input file in the following format:'
      write(*,*)
     1 '    POLE_LAT POLE_LON ANGULAR_VELOCITY'
      write(*,*)
     1 '    (Alternatively: OmegaX OmegaY OmegaZ; requires -p 2)'
      write(*,*)
     1 '    STNM STLA STLO N_DISP_OBS E_DISP_OBS'
      write(*,*)
     1 '    :           :            :'
      write(*,*)
     1 '    :           :            :'
      write(*,*)
      write(*,*)
     1 '  And calculates the predicted displacements at each station',
     2   ' given the pole of rotation.'
      write(*,*)
     1 '  The following is printed to the terminal or a file:'
      write(*,*)
     1 '    STNM STLA STLO N_DISP_OBS E_DISP_OBS N_DISP_PRE E_DISP_PRE'
      write(*,*)
     1 '    :           :            :               :'
      write(*,*)
     1 '    :           :            :               :'
      write(*,*)

      STOP
      END

C----------------------------------------------------------------------C

      SUBROUTINE gcmdln(ifile,ofile,pfmt)
      IMPLICIT none
      CHARACTER*80 tag,ifile,ofile
      INTEGER narg,i,pfmt
      
      ifile = 'none'
      ofile = 'none'

C     pfmt determines input format of rotation pole
C       1: latitude, longitude, angular velocity
C       2: omegaX, omegaY, omegaZ 
      pfmt = 1

      narg = iargc()
      if (narg.eq.0) call usage('POLE2VEL REQUIRES ARGUMENTS')

      i = 0
  101 i = i + 1
      if (i.gt.narg) goto 102
          call getarg(i,tag)
          if (tag(1:2).eq.'-h'.or.tag(1:2).eq.'-?') then
              call usage('ONLINE HELP')
          elseif (tag(1:2).eq.'-f') then
              i = i + 1
              call getarg(i,ifile)
          elseif (tag(1:2).eq.'-o') then
              i = i + 1
              call getarg(i,ofile)
          elseif (tag(1:2).eq.'-p') then
              i = i + 1
              call getarg(i,tag)
              read(tag,'(BN,I1)') pfmt
          endif
      goto 101
  102 continue

      if (ifile.eq.'none') call usage('INPUT FILE UNDEFINED')
      if (pfmt.ne.1.and.pfmt.ne.2)
     1                call usage('VARIABLE PFMT MUST BE 1 OR 2')

      RETURN
      END

