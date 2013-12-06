      PROGRAM distaz2lola

      IMPLICIT none
      INTEGER stdin,stdout
      REAL r2d
      PARAMETER (stdin=5,stdout=6,r2d=180.0/3.14159265)
      REAL lon1,lat1,lon2,lat2,dist,az
      CHARACTER*30 ifile,ofile
      INTEGER manual,p
      
      call gcmdln(ifile,ofile,manual,p)

      if (manual.eq.1) then
          print *,'Enter lon1 lat1 dist(km) az(deg):'
          read *,lon1,lat1,dist,az
          call lola(lon2,lat2,lon1,lat1,dist,az)
          write (6,9998) lon2*r2d,lat2*r2d
      else
          open (unit=11,file=ifile,status='old')
          open (unit=12,file=ofile,status='unknown')
  101     read (11,*,end=102) lon1,lat1,dist,az
              call lola(lon2,lat2,lon1,lat1,dist,az)
              if (p.eq.0) then
                  write (12,9999) lon2*r2d,lat2*r2d
              else
                  write (6,9999) lon2*r2d,lat2*r2d
              endif
              goto 101
  102     continue
      endif

 9998 format('(lon,lat) = (',F9.4,',',F8.4,')')
 9999 format(2F15.4)

      END

C======================================================================C

 
      SUBROUTINE gcmdln(ifile,ofile,manual,p)
      
      IMPLICIT NONE
      CHARACTER*30 ifile,ofile,tag
      INTEGER i,narg,manual,p
      
      manual = 0
      p = 0
      ifile = 'none'
      ofile = 'distaz2lola.out'
      
      narg = iargc()
      if (narg.eq.0) call usage()
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
          elseif (tag(1:2).eq.'-m') then
              manual = 1
          elseif (tag(1:2).eq.'-p') then
              p = 1
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
     1 'Usage: distaz2lola -f [IFILE] -o [OFILE] -p -m -h/-?'
      write(*,*)
     1 '  -f [IFILE]      name of input file'
      write(*,*)
     1 '  -o [OFILE]      name of output file'
      write(*,*)
     1 '                    if -o option not used, OFILE = ',
     2                   '"distaz2lola.out"'
      write(*,*)
     1 '  -p (Default no) print to standard out instead of file'
      write(*,*)
     1 '  -m (Default no) manually input origin, distance, and azimuth'
      write(*,*)
     1 '                    and print results to standard out'
      write (*,*)
     1 '  -h/-?           this online help'
      write (*,*) ''
      write (*,*)
     1 'distaz2lola calculates end coordinates given starting'
      write (*,*)
     1 '  coordinates, distance (in km) and bearing (in degrees)'
      write (*,*) ''
      write (*,*)
     1 '    input file [IFILE] format'
      write (*,*)
     1 '        lon1 lat1 dist(km) az(deg)'
      write (*,*)
     1 '         :       :'
      write (*,*)
      write (*,*)
     1 '    output file [OFILE] format'
      write (*,*)
     1 '        lon2 lat2 '
      write (*,*)
     1 '         :     :'
      write (*,*) ''

      STOP
      END
