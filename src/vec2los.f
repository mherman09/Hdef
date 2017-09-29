      PROGRAM vec2los
C----
C Take three-dimensional N-E-Z displacements and compute line-of-sight
C displacements given the azimuth (from north) and inclination (from
C horizontal) of the view
C----
      IMPLICIT NONE
      REAL*8 pi,d2r
      PARAMETER (pi=4.0d0*atan(1.0d0),d2r=pi/180.0d0)
      REAL*8 uLOS,stlo,stla,stdp,uN,uE,uZ
      REAL*8 az,incl,cosaz,cosinc,sinaz,sininc
      CHARACTER*40 ifile,ofile
      LOGICAL ex

C----
C Parse command line
C----
      call gcmdln(az,incl,ifile,ofile)
      if (az.lt.-1000.0d0.or.incl.lt.-1000.0d0) then
          write(*,*) '!! Error: azimuth or inclination are undefined'
          call usage('!! Use -a AZ or -i INCL')
      endif
      if (ifile.eq.'none') then
          write(*,*) '!! Error: Input displacement vector file ',
     1               'unspecified'
          call usage('!! Use -f IFILE to specify input file')
      else
          inquire(file=ifile,EXIST=ex)
          if (.not.ex) call usage('!! Error: no input file: '//ifile)
      endif
      if (ofile.eq.'none') then
          write(*,*) '!! Error: No output file specified'
          call usage('!! Use -o OFILE to specify LOS output file')
      endif
      cosaz = dcos(az*d2r)
      sinaz = dsin(az*d2r)
      cosinc = dcos(incl*d2r)
      sininc = dsin(incl*d2r)
      !print *,'AZ',az
      !print *,'INCL',incl
      !print *,'IFILE ',ifile 
      !print *,'OFILE ',ofile

C----
C Control files
C----
      open (unit=11,file=ifile,status='old')
      open (unit=12,file=ofile,status='unknown')

C----
C Convert NEZ to LOS
C----
   11 read (11,*,end=12) stlo,stla,stdp,uE,uN,uZ
          uLOS = uN*cosaz*cosinc
     1             + uE*sinaz*cosinc
     2                 - uZ*sininc
          write (12,9991) stlo,stla,stdp,uLOS
          goto 11
   12 continue

 9991 format(3F12.4,X,3E12.4)

      END

C======================================================================C

      SUBROUTINE gcmdln(az,incl,ifile,ofile)
      IMPLICIT none
      REAL*8 az,incl
      CHARACTER*40 tag,ifile,ofile
      INTEGER i,narg
      az = -99999.0d0
      incl = -99999.0d0
      ifile = 'none'
      ofile = 'none'
      narg = iargc()
      if (narg.eq.0) then
          call usage('!! Error: no command line arguments specified')
      endif
      i = 0
   11 i = i + 1
      if (i.gt.narg) goto 12
          call getarg(i,tag)
          if (tag(1:2).eq.'-a') then
              i = i + 1
              call getarg(i,tag)
              read (tag,'(BN,F12.0)') az
          elseif (tag(1:2).eq.'-i') then
              i = i + 1
              call getarg(i,tag)
              read (tag,'(BN,F12.0)') incl
          elseif (tag(1:2).eq.'-f') then
              i = i + 1
              call getarg(i,ifile)
          elseif (tag(1:2).eq.'-o') then
              i = i + 1
              call getarg(i,ofile)
          elseif (tag(1:2).eq.'-h') then
              call usage(' ')
          endif
          goto 11
   12 continue
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
     1 'Usage: vec2los -a AZ -i INC -f IFILE -o OFILE [-h]'
      write(*,*)
      write(*,*)
     1 '-a AZ     Look azimuth (CW from N)'
      write(*,*)
     1 '-i INC    Look inclination (from horizontal)'
      write(*,*)
     1 '-f IFILE  Input vector displacement file (stlo stla stdp uE ',
     2                                                'uN uZ)'
      write(*,*)
     1 '-o OFILE  Output LOS displacement file (stlo stla stdp uLOS)'
      write (*,*)
     1 '-h        Online help (this message)'
      write (*,*)
      STOP
      END

