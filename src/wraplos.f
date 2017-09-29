      PROGRAM wraplos
c----
c Convert a line of sight displacement to a phase difference given
c the wavelength of the InSAR satellite.
c----
      IMPLICIT NONE
      REAL*8 pi,tpi
      PARAMETER (pi=4.0d0*datan(1.0d0),tpi=2.0d0*pi)
      REAL*8 stlo,stla,stdp,ulos,wvl,phase
      CHARACTER*40 ifile,ofile
      LOGICAL ex

      call gcmdln(wvl,ifile,ofile)
      if (wvl.le.0.0d0) then
          write(*,*) '!! Error: wavelength is undefined'
          call usage('!! Use -w WVLNTH')
      endif
      if (ifile.eq.'none') then
          write(*,*) '!! Error: Input LOS vector file unspecified'
          call usage('!! Use -f IFILE to specify input file')
      else
          inquire(file=ifile,EXIST=ex)
          if (.not.ex) call usage('!! Error: no input file: '//ifile)
      endif
      if (ofile.eq.'none') then
          write(*,*) '!! Error: No output file specified'
          call usage('!! Use -o OFILE to specify phase output file')
      endif

      open (unit=10,file=ifile,status='old')
      open (unit=20,file=ofile,status='unknown')

c----
c Convert line of sight displacement to phase difference from satellite
c----
   13 read (10,*,end=14) stlo,stla,stdp,ulos
          ulos = 100.0d0*ulos
          phase = 2.0d0*tpi*ulos/wvl
          phase = dmod(phase,tpi)
          if (phase.lt.-pi) phase = phase + tpi
          if (phase.gt.pi) phase = phase - tpi
          write (20,9991) stlo,stla,stdp,phase
          goto 13
   14 continue

 9991 format(3F12.4,X,3E12.4)

      END

C======================================================================C

      SUBROUTINE gcmdln(wvl,ifile,ofile)
      IMPLICIT none
      REAL*8 wvl
      CHARACTER*40 tag,ifile,ofile
      INTEGER i,narg
      wvl = -1.0d0
      ifile = 'none'
      ofile = 'none'
      narg = iargc()
      if (narg.eq.0) call usage(' ')
      i = 0
   11 i = i + 1
      if (i.gt.narg) goto 12
          call getarg(i,tag)
          if (tag(1:2).eq.'-w') then
              i = i + 1
              call getarg(i,tag)
              read (tag,'(BN,F12.0)') wvl
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
     1 'Usage: wraplos -w WVLNTH -f IFILE -o OFILE [-h]'
      write(*,*)
      write(*,*)
     1 '-w WVLNTH  Wavelength of radar signal (cm)'
      write(*,*)
     1 '-f IFILE   Input LOS displacement file (stlo stla stdp uLOS)'
      write(*,*)
     1 '-o OFILE   Output phase file (stlo stla stdp phase(rad))'
      write (*,*)
     1 '-h         Online help (this message)'
      write (*,*)
      STOP
      END

