      PROGRAM lola2distaz

      IMPLICIT none
      REAL*8 pi,r2d,radius
      PARAMETER (pi=4.0d0*datan(1.0d0),r2d=1.8d2/pi,radius=6371.0d0)
      REAL*8 lon1,lat1,lon2,lat2,dist,az,cmdin(4)
      CHARACTER*30 ifile,ofile
      INTEGER c,p,input
      LOGICAL ex

C----
C Parse command line and read options
C----
      call gcmdln(ifile,ofile,c,cmdin,p)

C----
C Check inputs are valid
C----
      if (c.eq.1) goto 103
      if (ifile.eq.'none') then
          write(*,*) '!! Error: Input lon1 lat1 lon2 lat2 file ',
     1               'unspecified'
          call usage('!! Use -f IFILE to specify input file')
      elseif (ifile.ne.'stdin') then
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
          if (ifile.eq.'stdin') then
              input = 5
          else
              input = 11
              open (unit=input,file=ifile,status='old')
          endif
          if (p.eq.0) open (unit=12,file=ofile,status='unknown')
  101     read (input,*,err=991,end=102) lon1,lat1,lon2,lat2
              if (dabs(lat1).gt.90.0d0.or.dabs(lat2).gt.90.0d0) then
                  call errmes('!! Error: latitude is greater than 90')
              endif
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
          if (ifile.ne.'stdin') close(input)
      endif

      goto 1000
  991 call errmes('!! Error: problem with reading input (lon1 lat1 '//
     1                                                     'lon2 lat2)')
 8889 format(2F20.8)
 9999 format(2F20.8)
 1000 continue

      END

C======================================================================C

      SUBROUTINE errmes(text)
      IMPLICIT none
      CHARACTER text*(*)
      write(0,*) text
      STOP
      END

C----------------------------------------------------------------------C

      SUBROUTINE gcmdln(ifile,ofile,c,cmdin,p)
      IMPLICIT NONE
      CHARACTER*30 ifile,ofile,tag
      INTEGER i,narg,c,p
      REAL*8 cmdin(4)
      ifile = 'none'
      ofile = 'none'
      c = 0
      p = 1
      narg = iargc()
      if (narg.eq.0) call usage('')
      i = 0
  901 i = i + 1
      if (i.gt.narg) goto 902
          call getarg(i,tag)
          if (tag(1:2).eq.'-f') then
              i = i + 1
              call getarg(i,ifile)
          elseif (tag(1:2).eq.'-o') then
              p = 0
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
          elseif (tag(1:2).eq.'-s') then
              ifile = 'stdin'
          elseif (tag(1:2).eq.'-h') then
              call usage(' ')
          else
              call usage('!! Error: no option '//trim(tag))
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
     1 'Usage: lola2distaz -f IFILE|-s [-o OFILE|-p] [-c LON1 LAT1 ',
     2                       'LON2 LAT2] [-h]'
      write(*,*)
      write(*,*)
     1 '-f IFILE               Input file: lon1 lat1 lon2 lat2'
      write(*,*)
     1 '-s                     Read standard input (overrides -f)'
      write(*,*)
     1 '-p                     Print to standard output (default)'
      write(*,*)
     1 '-o OFILE               Output file: dist(km) az(deg)'
      write(*,*)
     1 '-c LON1 LAT1 LON2 LAT2 Compute DIST AZ, ',
     2                           'print to standard output'
      write(*,*)
     1 '-h                     Online help (this screen)'
      write(*,*)
      STOP
      END
