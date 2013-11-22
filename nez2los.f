      PROGRAM nez2los
C----
C Take three-dimensional N-E-Z displacements and compute line-of-sight
C displacements given the azimuth (from north) and inclination (from
C horizontal) of the view
C----
      IMPLICIT NONE
      REAL pi,d2r
      PARAMETER (pi=4.0*atan(1.0),d2r=pi/180.0)

      REAL uLOS,stlo,stla,uN,uE,uZ
      REAL az,incl
      CHARACTER*30 ifile,ofile

C----
C Parse command line
C----
      call gcmdln(az,incl,ifile,ofile)
      az = az*d2r
      incl = incl*d2r

C----
C Control files
C----
      open (unit=10,file=ifile,status='old')
      open (unit=20,file=ofile,status='unknown')
      rewind 20

C----
C Convert NEZ to LOS
C----
   11 read (10,*,end=12) stlo,stla,uN,uE,uZ
          uLOS = uN*cos(az)*cos(incl)
     1             + uE*sin(az)*cos(incl)
     2                 - uZ*sin(incl)
          write (20,*) stlo,stla,uLOS
      goto 11

   12 continue

      END

C======================================================================C

      SUBROUTINE gcmdln(az,incl,ifile,ofile)

      IMPLICIT none
      REAL az,incl
      CHARACTER*30 tag,ifile,ofile
      INTEGER i,narg

      az = 0.0
      incl = 90.0
      ifile = 'disp.out'
      ofile = 'nez2los.out'

      narg = iargc()
      if (narg.eq.0) call usage()

      i = 0
   11 i = i + 1
      if (i.gt.narg) goto 12
          call getarg(i,tag)
          if (tag(1:2).eq.'-a') then
              i = i + 1
              call getarg(i,tag)
              read (tag,'(BN,F6.0)') az
          elseif (tag(1:2).eq.'-i') then
              i = i + 1
              call getarg(i,tag)
              read (tag,'(BN,F4.0)') incl
          elseif (tag(1:2).eq.'-f') then
              i = i + 1
              call getarg(i,ifile)
          elseif (tag(1:2).eq.'-o') then
              i = i + 1
              call getarg(i,ofile)
          elseif (tag(1:2).eq.'-h'.or.tag(1:2).eq.'-?') then
              call usage()
          endif
          goto 11
   12 continue

      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE usage()
      IMPLICIT none

      write(*,*)
     1 'Usage: nez2los -a [AZ] -i [IN] -f [IFILE] -o [OFILE] -h/-?'
      write(*,*)
     1 '  -a [AZIMUTH] look azimuth (CW from N)'
      write(*,*)
     1 '  -i [INCLINATION] look inclination (from horizontal)'
      write(*,*)
     1 '  -f [IFILE] (Default disp.out) name of input file'
      write(*,*)
     1 '  -o [OFILE] (Default nez2los.out) name of output file'
      write (*,*)
     1 '  -h/-?        help'
      write (*,*) ''
      write (*,*)
     1 '  nez2los converts NEZ to line-of-sight displacements'
      write (*,*) ''
      write (*,*)
     1 '    input file'
      write (*,*)
     1 '        stlo stla uN uE uZ'
      write (*,*)
     1 '           :         :'
      write (*,*)
     1 '    output file'
      write (*,*)
     1 '        stlo stla uLOS'
      write (*,*)
     1 '           :         :'
      write (*,*) ''

      STOP
      END

