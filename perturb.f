      PROGRAM main
C----
C Take a data file with uncertainties and compute random perturbations
C to the values.
C----
      IMPLICIT none
      CHARACTER*80 ifile,ofile
      INTEGER i,nvar,input,output,idum
      REAL*8 val(10),dval(10)
      EXTERNAL ran2
      REAL ran2,random

C Parse command line
      call gcmdln(ifile,ofile,nvar)
      if (ifile.eq.'stdin') then
          input = 5
      else
          input = 11
          open(unit=input,file=ifile,status='old')
      endif
      if (ofile.eq.'print'.or.ofile.eq.'stdout') then
          output = 6
      else
          output = 12
          open(unit=output,file=ofile,status='unknown')
      endif

C Initialize random number generator
      call raninit(idum)

C Perturb input data
  101 read(input,*,end=102) (val(i),i=1,nvar),(dval(i),i=1,nvar)
          do 104 i = 1,nvar
              val(i) = val(i) + dval(i)*(2.0d0*ran2(idum)-1.0d0)
  104     continue
          write(output,*) (val(i),i=1,nvar)
          goto 101
  102 continue
      END

C----------------------------------------------------------------------C

      SUBROUTINE raninit(idum)
      IMPLICIT none
      INTEGER time(3),timeseed,idum
      CHARACTER*8 timec
      call itime(time)
      write(timec,'(I0.2I0.2I0.2)') time(1),time(2),time(3)
      read(timec,*) timeseed
      if (timeseed.gt.0) timeseed = -timeseed
      idum = timeseed
      RETURN
      END

C======================================================================C

      SUBROUTINE gcmdln(ifile,ofile,nvar)
      IMPLICIT NONE
      CHARACTER*80 ifile,ofile,tag
      INTEGER i,narg,nvar
      ifile = 'none'
      ofile = 'none'
      nvar = 1
      narg = iargc()
      if (narg.eq.0) call usage('')
      i = 0
  101 i = i + 1
      if (i.gt.narg) goto 102
          call getarg(i,tag)
          if (tag(1:2).eq.'-n') then
              i = i + 1
              call getarg(i,tag)
              read(tag,*) nvar
          elseif (tag(1:2).eq.'-o') then
              i = i + 1
              call getarg(i,ofile)
          elseif (tag(1:2).eq.'-f') then
              i = i + 1
              call getarg(i,ifile)
          elseif (tag(1:2).eq.'-h') then
              call usage('')
          else
              call usage('!! Error: No option '//tag)
          endif
          goto 101
  102 continue
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE usage(str)
      IMPLICIT NONE
      CHARACTER str*(*)
      if (str.ne.'') then
          write(*,*) trim(str)
          write(*,*)
      endif
      write(*,*)
     1 'Usage: perturb -f IFILE -o OFILE -n NVAR'
      write(*,*)
      write(*,*)
     1 '-f IFILE    Input data file with uncertainties (v1..vn d1..dn)'
      write(*,*)
     1 '-o OFILE    Output perturbed data'
      write(*,*)
     1 '-n NVAR     Number of data variables'
      write(*,*)
     1 '-h          Online help'
      write(*,*)
      STOP
      END
