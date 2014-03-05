      PROGRAM lola2distaz

      IMPLICIT none
      REAL*8 pi,r2d,radius
      PARAMETER (pi=4.0d0*datan(1.0d0),r2d=1.8d2/pi,radius=6371.0d0)
      REAL*8 lon1,lat1,lon2,lat2,dist,az
      CHARACTER*30 ifile,ofile
      INTEGER user,p
      
      call gcmdln(ifile,ofile,user,p)

      if (user.eq.1) then
          print *,'Enter lon1 lat1 lon2 lat2:'
          read *,lon1,lat1,lon2,lat2
          call ddistaz(dist,az,lon1,lat1,lon2,lat2)
          dist = dist*radius
          az = az*r2d
          write (*,8888) dist
          write (*,8889) az
      else
          open (unit=11,file=ifile,status='old')
          open (unit=12,file=ofile,status='unknown')
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

 8888 format('Distance =',F14.6,' km')
 8889 format('Azimuth  =',F14.6,' degrees')
 9999 format(2F14.6)

      END

C======================================================================C

      SUBROUTINE gcmdln(ifile,ofile,user,p)
      
      IMPLICIT NONE
      CHARACTER*30 ifile,ofile,tag
      INTEGER i,narg,user,p
      
      user = 0
      p = 0
      ifile = 'lola2distaz.in'
      ofile = 'lola2distaz.out'
      
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
          elseif (tag(1:2).eq.'-d') then
              write(*,*) 'Running with default file names'
          elseif (tag(1:2).eq.'-u') then
              user = 1
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
     1 'Usage: lola2distaz -f [IFILE] -o [OFILE] -d -p -u -h/-?'
      write(*,*)
     1 '  -f [IFILE] (Default lola2distaz.in) name of input file'
      write(*,*)
     1 '  -o [OFILE] (Default lola2distaz.out) name of output file'
      write(*,*)
     1 '  -d         Run with default file names'
      write(*,*)
     1 '  -p         Print to standard out instead of file'
      write (*,*)
     1 '  -u         Prompt user to enter information through standard'
      write (*,*)
     1 '                 input for single calculation'
      write (*,*)
     1 '  -h/-?      Online help (this screen)'
      write (*,*) ''
      write (*,*)
     1 '  lola2distaz calculates distance (in km) and azimuth (in ',
     2    'degrees) from coordinate pairs'
      write (*,*) ''
      write (*,*)
     1 '    Input file'
      write (*,*)
     1 '        lon1 lat1 lon2 lat2'
      write (*,*)
     1 '         :       :'
      write (*,*)
     1 '    Output file'
      write (*,*)
     1 '        dist(km) az(deg) (from coordinate 1 to 2)'
      write (*,*)
     1 '         :       :'
      write (*,*) ''

      STOP
      END
