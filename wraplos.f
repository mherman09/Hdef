      PROGRAM wraplos
c----
c Convert a line of sight displacement to a phase difference given
c the wavelength of the InSAR satellite.
c----
      IMPLICIT NONE
      REAL pi,d2r
      PARAMETER (pi=4.0*atan(1.0),d2r=pi/180.0)

      REAL stlo,stla,ulos,wvl,phase
      
      INTEGER narg,i
      CHARACTER*30 ifile,ofile

      call gcmdln(wvl,ifile,ofile)

      open (unit=10,file=ifile,status='old')
      open (unit=20,file=ofile,status='unknown')

c----
c Convert line of sight displacement to phase difference from satellite
c----
   13 read (10,*,end=14) stlo,stla,ulos
          ulos = 100.0*ulos
          phase = 2.0*(2.0*pi)*ulos/wvl
          phase = mod(phase,2.0*pi)
          if (phase.lt.-pi) phase = phase + 2.0*pi
          if (phase.gt.pi) phase = phase - 2.0*pi
          write (20,*) stlo,stla, phase
          goto 13
   14 continue

      END

C======================================================================C

      SUBROUTINE gcmdln(wvl,ifile,ofile)

      IMPLICIT none
      REAL wvl
      CHARACTER*30 tag,ifile,ofile
      INTEGER i,narg

      wvl = 23.6
      ifile = 'nez2los.out'
      ofile = 'wraplos.out'

      narg = iargc()
      if (narg.eq.0) call usage()

      i = 0
   11 i = i + 1
      if (i.gt.narg) goto 12
          call getarg(i,tag)
          if (tag(1:2).eq.'-w') then
              i = i + 1
              call getarg(i,tag)
              read (tag,'(BN,F6.0)') wvl
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
     1 'Usage: wraplos -w [WAVELENGTH] -f [IFILE] -o [OFILE] -h/-?'
      write(*,*)
     1 '  -w [WAVELENGTH]  wavelength of satellite ',
     2                                      'signal (in cm)'
      write(*,*)
     1 '  -f [IFILE] (Default nez2los.out) name of input file'
      write(*,*)
     1 '  -o [OFILE] (Default wraplos.out) name of output file'
      write (*,*)
     1 '  -h/-?        help'
      write (*,*) ''
      write (*,*)
     1 '  wraplos converts displacements to phase changes, given an ',
     2    'observation wavelength'
      write (*,*) ''
      write (*,*)
     1 '    input file'
      write (*,*)
     1 '        stlo stla uLOS'
      write (*,*)
     1 '         :  :'
      write (*,*)
     1 '    output file'
      write (*,*)
     1 '        stlo stla phase_change'
      write (*,*)
     1 '         :  :'
      write (*,*) ''

      STOP
      END

