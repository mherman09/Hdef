      PROGRAM pole2vel

      IMPLICIT none
      REAL*8 pi,d2r
      PARAMETER (pi=4.0d0*datan(1.0d0),d2r=pi/180.0d0)
      CHARACTER*80 ifile,ofile
      LOGICAL ex
      INTEGER pfmt,user,pr
      CHARACTER*8 stnm
      REAL*8 plat,plon,pmag,p(3)
      REAL*8 stlo,stla,la,lo,uN,uE,r(3)
      REAL*8 v(3),vN,vE,vZ

C----
C Get command line arguments
C----
      call gcmdln(ifile,ofile,pfmt,user,pr)

C----
C Calculate predicted velocities
C----
      if (user.eq.1) then
          print *,'Pole format?'
          print *,'  1: lat,lon,ang_vel(deg/Ma)'
          print *,'  2: omegax,omegay,omegaz (all deg/Ma)'
          read *,pfmt
          if (pfmt.eq.1) then
              print *,'Enter lat lon ang_vel(deg/Ma):'
              read *,plat,plon,pmag
              call lalomag2xyz(p(1),p(2),p(3),plat,plon,pmag)
          elseif (pfmt.eq.2) then
              print *,'omegax omegay omegaz:'
              read *,p(1),p(2),p(3)
          endif
          print *,'Enter lat lon:'
          read *,stla,stlo
          call lalomag2xyz(r(1),r(2),r(3),stla,stlo,1.0d0)
          call cross(v,p,r)
          la = stla*d2r
          lo = stlo*d2r
          vN = -(v(1)*dsin(la)*dcos(lo)+v(2)*dsin(la)*dsin(lo)
     1                                                  -v(3)*dcos(la))
          vE =  -v(1)         *dsin(lo)+v(2)         *dcos(lo)
          vZ =   v(1)*dcos(la)*dcos(lo)+v(2)*dcos(la)*dsin(lo)
     1                                                   +v(3)*dsin(la)
          vN = vN*111.19d0
          vE = vE*111.19d0
          write (*,8888) vN
          write (*,8889) vE
      else
          open(unit=11,file=ifile,status='old')
          if (pr.eq.0) open(unit=12,file=ofile,status='unknown')
          if (pfmt.eq.1) then
              read(11,*) plat,plon,pmag
              call lalomag2xyz(p(1),p(2),p(3),plat,plon,pmag)
          else
              read(11,*) p(1),p(2),p(3)
          endif
  101     read(11,*,end=102) stnm,stla,stlo,uN,uE
              call lalomag2xyz(r(1),r(2),r(3),stla,stlo,1.0d0)
              call cross(v,p,r)
              la = stla*d2r
              lo = stlo*d2r
              vN = -(v(1)*dsin(la)*dcos(lo)+v(2)*dsin(la)*dsin(lo)
     1                                                   -v(3)*dcos(la))
              vE =  -v(1)         *dsin(lo)+v(2)         *dcos(lo)
              vZ =   v(1)*dcos(la)*dcos(lo)+v(2)*dcos(la)*dsin(lo)
     1                                                    +v(3)*dsin(la)
              vN = vN*111.19d0
              vE = vE*111.19d0
              if (pr.eq.0) then
                  write(12,9999) stnm,stla,stlo,uN,uE,vN,vE
              else
                  write(*,9999) stnm,stla,stlo,uN,uE,vN,vE
              endif
              goto 101
  102     continue
      endif

 8888 format('N velocity =',F8.2,' km/Ma')
 8889 format('E velocity =',F8.2,' km/Ma')
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

      SUBROUTINE usage()
      IMPLICIT none

      write(*,*)
     1 'Usage: pole2vel -f [IFILE] -o [OFILE] -l [POLE_FORMAT] -d -p',
     2                ' -u -h/-?'
      write(*,*)
     1 '  -f [IFILE] (Default pole2vel.in) Name of input file'
      write(*,*)
     1 '  -o [OFILE] (default pole2vel.out) Name of output file'
      write(*,*)
     1 '  -l [POLE_FORMAT] (default 1) rotation pole format'
      write(*,*)
     1 '        1: Latitude, Longitude, Angular Velocity(deg/Ma)'
      write(*,*)
     1 '        2: OmegaX(deg/Ma), OmegaY(deg/Ma), OmegaZ(deg/Ma)'
      write(*,*)
     1 '  -d         Run with default file names'
      write(*,*)
     1 '  -p         Print to standard out instead of file'
      write (*,*)
     1 '  -u         Prompt user to enter information through standard'
      write (*,*)
     1 '                 input for single calculation'
      write(*,*)
     1 '  -h/-? Online help (this screen)'
      write(*,*)
      write(*,*)
     1 'pole2vel calculates the velocities at geographical',
     2    ' coordinates given an Euler pole.'
      write(*,*) ''
      write (*,*)
     1 '    Input file'
      write(*,*)
     1 '      pole_lat pole_lon ang_vel(deg/Ma)'
      write(*,*)
     1 '      (Alternatively: OmegaX OmegaY OmegaZ; requires -p 2)'
      write(*,*)
     1 '      stnm stla stlo N_vel_obs E_vel_obs'
      write(*,*)
     1 '                :            :'
      write(*,*) ''
      write(*,*)
     1 '    Output file'
      write(*,*)
     1 '      stnm stla stlo N_vel_obs E_vel_obs N_vel_pre E_vel_pre'
      write(*,*)
     1 '                :            :'
      write(*,*)

      STOP
      END

C----------------------------------------------------------------------C

      SUBROUTINE gcmdln(ifile,ofile,pfmt,user,pr)
      IMPLICIT none
      CHARACTER*80 tag,ifile,ofile
      INTEGER narg,i,pfmt,user,pr
      
      ifile = 'pole2vel.in'
      ofile = 'pole2vel.out'
      user = 0
      pr = 0

C     pfmt determines input format of rotation pole
C       1: latitude, longitude, angular velocity
C       2: omegaX, omegaY, omegaZ 
      pfmt = 1

      narg = iargc()
      if (narg.eq.0) call usage()

      i = 0
  101 i = i + 1
      if (i.gt.narg) goto 102
          call getarg(i,tag)
          if (tag(1:2).eq.'-h'.or.tag(1:2).eq.'-?') then
              call usage()
          elseif (tag(1:2).eq.'-f') then
              i = i + 1
              call getarg(i,ifile)
          elseif (tag(1:2).eq.'-o') then
              i = i + 1
              call getarg(i,ofile)
          elseif (tag(1:2).eq.'-l') then
              i = i + 1
              call getarg(i,tag)
              read(tag,'(BN,I1)') pfmt
              if (pfmt.ne.1.and.pfmt.ne.2) call usage
          elseif (tag(1:2).eq.'-d') then
              write(*,*) 'Running with default file names'
          elseif (tag(1:2).eq.'-u') then
              user = 1
          elseif (tag(1:2).eq.'-p') then
              pr = 1
          endif
      goto 101
  102 continue

      RETURN
      END

