      PROGRAM main
C----
C Take a data file with uncertainties and compute random perturbations
C to the values.
C----
      IMPLICIT none
      CHARACTER*80 ifile,ofile
      INTEGER i,nvar,input,output,idum,seed,p,distro
      INTEGER irand
      REAL*8 val(10),dval(10),r1,r2
      EXTERNAL ran2
      REAL ran2

C Parse command line
      call gcmdln(ifile,ofile,nvar,seed,p,distro)
      if (distro.eq.2) irand = 0
      if (ifile.eq.'stdin') then
          input = 5
      else
          input = 11
          open(unit=input,file=ifile,status='old')
      endif
      if (ofile.eq.'print'.or.ofile.eq.'stdout'.or.ofile.eq.'none'
     1    .or.p.eq.1) then
          output = 6
      else
          output = 12
          open(unit=output,file=ofile,status='unknown')
      endif

C Initialize random number generator
      call raninit(idum)
      if (seed.gt.0) idum = -seed

C Perturb input data
  101 read(input,*,end=102) (val(i),i=1,nvar),(dval(i),i=1,nvar)
          do 104 i = 1,nvar
              if (distro.eq.1) then
                  val(i) = val(i) + dval(i)*(2.0d0*ran2(idum)-1.0d0)
              elseif (distro.eq.2) then
                  if (irand.eq.0) then
                      call nrand(r1,r2,idum)
                      irand = 1
                  else
                      r1 = r2
                      irand = 0
                  endif
                  val(i) = val(i) + dval(i)*r1/1.96d0
              endif
  104     continue
          write(output,*) (val(i),i=1,nvar)
          goto 101
  102 continue

      write(0,*) idum
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

      SUBROUTINE gcmdln(ifile,ofile,nvar,seed,p,distro)
      IMPLICIT NONE
      CHARACTER*80 ifile,ofile,tag
      INTEGER i,narg,nvar,seed,p,distro
      ifile  = 'none'
      ofile  = 'none'
      seed   = -1
      nvar   = 1
      p      = 0
      distro = 1
      narg = iargc()
      if (narg.eq.0) call usage('')
      i = 0
  101 i = i + 1
      if (i.gt.narg) goto 102
          call getarg(i,tag)
          if (tag(1:5).eq.'-seed') then
              i = i + 1
              call getarg(i,tag)
              read(tag,*) seed
          elseif (tag(1:7).eq.'-normal') then
              distro = 2
          elseif (tag(1:2).eq.'-n') then
              i = i + 1
              call getarg(i,tag)
              read(tag,*) nvar
          elseif (tag(1:2).eq.'-o') then
              i = i + 1
              call getarg(i,ofile)
          elseif (tag(1:2).eq.'-f') then
              i = i + 1
              call getarg(i,ifile)
          elseif (tag(1:2).eq.'-p') then
              p = 1
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
     1 'Usage: perturb -f IFILE -o OFILE -n NVAR [-seed SEED] [-p] ',
     2                '[-normal]'
      write(*,*)
      write(*,*)
     1 '-f IFILE    Input data file with uncertainties (v1..vn d1..dn)'
      write(*,*)
     1 '-o OFILE    Output perturbed data'
      write(*,*)
     1 '-n NVAR     Number of data variables'
      write(*,*)
     1 '-seed SEED  Start random numbers with integer SEED ',
     2                  '(useful for continuing chain of rands)'
      write(*,*)
     1 '-normal     Use normal distribution with 95% confidence at',
     2              'd values (default: uniform distribution)'
      write(*,*)
     1 '-p          Print results to standard output'
      write(*,*)
     1 '-h          Online help'
      write(*,*)
      STOP
      END
