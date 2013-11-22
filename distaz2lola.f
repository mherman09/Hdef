      PROGRAM distaz2lola

      IMPLICIT none
      REAL lon1,lat1,lon2,lat2,dist,az,r2d
      CHARACTER*30 ifile,ofile
      
      r2d = 180.0/3.14159265

      call gcmdln(ifile,ofile)

      if (ifile.eq.'none') then
          print *,'To see options, use -h flag'
          print *,'Enter lon1 lat1 dist(km) az(deg):'
          read *,lon1,lat1,dist,az
          call lola(lon2,lat2,lon1,lat1,dist,az)
          print *,'(lon2,lat2) = ',lon2*r2d,lat2*r2d
      else
          open (unit=11,file=ifile,status='old')
          open (unit=12,file=ofile,status='unknown')
  101     read (11,*,end=102) lon1,lat1,dist,az
              call lola(lon2,lat2,lon1,lat1,dist,az)
              write (12,*) lon2*r2d,lat2*r2d
              goto 101
  102     continue
      endif

      END

C======================================================================C

 
      SUBROUTINE gcmdln(ifile,ofile)
      
      IMPLICIT NONE
      CHARACTER*30 ifile,ofile,tag
      INTEGER i,narg
      
      ifile = 'none'
      ofile = 'distaz2lola.out'
      
      narg = iargc()
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
     1 'Usage: distaz2lola -f [IFILE] -o [OFILE] -h/-?'
      write(*,*)
     1 '  -f [IFILE] (Default none) name of input file'
      write(*,*)
     1 '                            if no file name given, fault ',
     2                            'parameters prompted for manual input'
      write(*,*)
     1 '  -o [OFILE] (Default lola2distaz.out) name of output file'
      write (*,*)
     1 '  -h/-?        help'
      write (*,*) ''
      write (*,*)
     1 '  distaz2lola calculates end coordinates given starting ',
     2    'coordinates, distance (in km) and bearing (in degrees)'
      write (*,*) ''
      write (*,*)
     1 '    input file'
      write (*,*)
     1 '        lon1 lat1 dist(km) az(deg)'
      write (*,*)
     1 '         :       :'
      write (*,*)
     1 '    output file'
      write (*,*)
     1 '        lon2 lat2 '
      write (*,*)
     1 '         :     :'
      write (*,*) ''

      STOP
      END
