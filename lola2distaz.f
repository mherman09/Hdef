      PROGRAM lola2distaz

      IMPLICIT none
      REAL*8 pi,r2d,radius
      PARAMETER (pi=4.0d0*datan(1.0d0),r2d=1.8d2/pi,radius=6371.0d0)
      REAL*8 lon1,lat1,lon2,lat2,dist,az,cmdin(4)
      CHARACTER*30 ifile,ofile
      INTEGER c,p
      LOGICAL ex
      
      call gcmdln(ifile,ofile,c,cmdin,p)
      if (c.eq.1) goto 103
      if (ifile.eq.'none') then
          write(*,*) '!! Error: Input lon1 lat1 lon2 lat2 file ',
     1               'unspecified'
          call usage('!! Use -f IFILE to specify input file')
      else
          inquire(file=ifile,EXIST=ex)
          if (.not.ex) call usage('!! Error: no input file: '//ifile)
      endif
      if (ofile.eq.'none'.and.p.eq.0) then
          write(*,*) '!! Error: No output file specified'
          call usage('!! Use -o OFILE to specify dist az output file')
      endif

  103 if (c.eq.1) then
          lon1 = cmdin(1)
          lat1 = cmdin(2)
          lon2 = cmdin(3)
          lat2 = cmdin(4)
          call ddistaz(dist,az,lon1,lat1,lon2,lat2)
          dist = dist*radius
          az = az*r2d
          write (*,8889) dist,az
      else
          open (unit=11,file=ifile,status='old')
          if (p.eq.0) open (unit=12,file=ofile,status='unknown')
  101     read (11,*,end=102) lon1,lat1,lon2,lat2
              call ddistaz(dist,az,lon1,lat1,lon2,lat2)
              dist = dist*radius
              az = az*r2d
              if (p.eq.0) then
                  write (12,9999) dist,az
              else
                  write (*,9999) dist,az
              endif
              goto 101
  102     continue
      endif

 8889 format(2F18.6)
 9999 format(2F18.6)

      END

C======================================================================C

      SUBROUTINE gcmdln(ifile,ofile,c,cmdin,p)
      IMPLICIT NONE
      CHARACTER*30 ifile,ofile,tag
      INTEGER i,narg,c,p
      REAL*8 cmdin(4)
      ifile = 'none'
      ofile = 'none'
      c = 0
      p = 0
      narg = iargc()
      if (narg.eq.0) call usage('!! Error: no command line arguments')
      i = 0
  901 i = i + 1
      if (i.gt.narg) goto 902
          call getarg(i,tag)
          if (tag(1:2).eq.'-f') then
              i = i + 1
              call getarg(i,ifile)
          elseif (tag(1:2).eq.'-o') then
              i = i + 1
              call getarg(i,ofile)
          elseif (tag(1:2).eq.'-c') then
              c = 1
              i = i + 1
              call getarg(i,tag)
              read(tag,'(BN,F12.0)') cmdin(1)
              i = i + 1
              call getarg(i,tag)
              read(tag,'(BN,F12.0)') cmdin(2)
              i = i + 1
              call getarg(i,tag)
              read(tag,'(BN,F12.0)') cmdin(3)
              i = i + 1
              call getarg(i,tag)
              read(tag,'(BN,F12.0)') cmdin(4)
          elseif (tag(1:2).eq.'-p') then
              p = 1
          elseif (tag(1:2).eq.'-h') then
              call usage(' ')
          endif
          goto 901
  902 continue
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
     1 'Usage: lola2distaz -f IFILE -o OFILE [-c LON1 LAT1 LON2 LAT2] ',
     2         '[-p] [-h]'
      write(*,*)
      write(*,*)
     1 '-f IFILE               Input file: lon1 lat1 lon2 lat2'
      write(*,*)
     1 '-o OFILE               Output file: dist(km) az(deg)'
      write(*,*)
     1 '-c LON1 LAT1 LON2 LAT2 Compute DIST AZ, ',
     2                           'print to standard output'
      write(*,*)
     1 '-p                     Print to standard output (overrides -o)'
      write(*,*)
     1 '-h                     Online help (this screen)'
      write(*,*)
      STOP
      END
