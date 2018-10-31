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
      CHARACTER*256 azincfile
      LOGICAL ex
      INTEGER uin,uout

C----
C Parse command line
C----
      call gcmdln(az,incl,ifile,ofile,azincfile)
      if ((az.lt.-1000.0d0.or.incl.lt.-1000.0d0)
     1                                .and.azincfile.eq.'none') then
          write(*,*) '!! Error: azimuth or inclination are undefined'
          call usage('!! Use -a AZ, -i INCL, or -azinc FILE')
      endif
      if (ifile.eq.'none') then
          write(*,*) '!! Error: Input displacement vector file ',
     1               'unspecified'
          call usage('!! Use -f IFILE to specify input file')
      elseif (ifile.ne.'stdin') then
          inquire(file=ifile,EXIST=ex)
          if (.not.ex) call usage('!! Error: no input file: '//ifile)
      endif
      if (ofile.eq.'none') then
          write(*,*) '!! Error: No output file specified'
          call usage('!! Use -o OFILE to specify LOS output file')
      endif

C----
C Control files
C----
      if (ifile.eq.'stdin') then
          uin = 5
      else
          uin = 11
          open (unit=uin,file=ifile,status='old')
      endif
      if (ofile.eq.'stdout') then
          uout = 6
      else
          uout = 12
          open (unit=uout,file=ofile,status='unknown')
      endif

C----
C Convert NEZ to LOS
C----
      cosaz = 0.0d0
      sinaz = 0.0d0
      cosinc = 0.0d0
      sininc = 0.0d0
      if (azincfile.ne.'none') then
          open(unit=11,file=azincfile,status='old')
      else
          cosaz = dcos(az*d2r)
          sinaz = dsin(az*d2r)
          cosinc = dcos(incl*d2r)
          sininc = dsin(incl*d2r)
      endif
      !print *,'AZ',az
      !print *,'INCL',incl
      !print *,'IFILE ',ifile 
      !print *,'OFILE ',ofile
   11 read (uin,*,end=12) stlo,stla,stdp,uE,uN,uZ
          if (azincfile.ne.'none') then
              read(11,*) az,incl
              cosaz = dcos(az*d2r)
              sinaz = dsin(az*d2r)
              cosinc = dcos(incl*d2r)
              sininc = dsin(incl*d2r)
          endif
          uLOS = uN*cosaz*cosinc
     1             + uE*sinaz*cosinc
     2                 - uZ*sininc
          write (uout,9991) stlo,stla,stdp,uLOS
          goto 11
   12 continue

 9991 format(3F12.4,X,3E12.4)
      if (uin.eq.11) close(11)
      if (uout.eq.12) close(12)
      END

C======================================================================C

      SUBROUTINE gcmdln(az,incl,ifile,ofile,azincfile)
      IMPLICIT none
      REAL*8 az,incl
      CHARACTER*40 tag,ifile,ofile
      CHARACTER*256 azincfile
      INTEGER i,narg
      az = -99999.0d0
      incl = -99999.0d0
      ifile = 'stdin'
      ofile = 'stdout'
      azincfile = 'none'
      narg = iargc()
      if (narg.eq.0) then
          call usage('')
      endif
      i = 0
   11 i = i + 1
      if (i.gt.narg) goto 12
          call getarg(i,tag)
          if (trim(tag).eq.'-a') then
              i = i + 1
              call getarg(i,tag)
              read (tag,'(BN,F12.0)') az
          elseif (trim(tag).eq.'-i') then
              i = i + 1
              call getarg(i,tag)
              read (tag,'(BN,F12.0)') incl
          elseif (trim(tag).eq.'-azinc') then
              i = i + 1
              call getarg(i,azincfile)
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
     1 'Usage: vec2los -a AZ -i INC | -azinc file [-f IFILE] ',
     2                 '[-o OFILE] [-h]'
      write(*,*)
      write(*,*)
     1 '-a AZ     Look azimuth (CW from N)'
      write(*,*)
     1 '-i INC    Look inclination (from horizontal=0 to vertical=90)'
      write(*,*)
     1 '-azinc FILE Azimuth and inclination in a file'
      write(*,*)
     1 '-f IFILE  Input vector displacement file (default: stdin)',
     2            ' (stlo stla stdp uE uN uZ)'
      write(*,*)
     1 '-o OFILE  Output LOS displacement file (default: stdout)',
     2            ' (stlo stla stdp uLOS)'
      write (*,*)
     1 '-h        Online help (this message)'
      write (*,*)
      STOP
      END

