      PROGRAM numint

      IMPLICIT none
      INTEGER NMAX
      PARAMETER (NMAX=100000)
      REAL*8 x(NMAX),f(NMAX),fint(NMAX),dx
      CHARACTER*30 ifile,ofile
      LOGICAL ex
      INTEGER i,ntot

C----
C Get file names from command line
C----
      write(0,*) 'WARNING: NUMINT HAS NOT BEEN TESTED AT ALL!'
      call gcmdln(ifile,ofile)
      inquire(file=ifile,EXIST=ex)
      if (.not.ex) then
          write(*,*) 'No input file:',ifile
          call usage()
      endif

C----
C Read data series
C----
      open(unit=11,file=ifile,status='old')
      i = 1
  101 read(11,*,end=102) x(i),f(i)
          i = i + 1
          goto 101
  102 continue
      ntot = i - 1

C----
C Trapezoidal integration
C----
      fint(1) = 0.0d0
      do 103 i = 2,ntot
          dx = x(i) - x(i-1)
          fint(i) = fint(i-1) + dx*0.5d0*(f(i)+f(i-1))
  103 continue

C----
C Write integrated data series
C----
      open(unit=12,file=ofile,status='unknown')
      do 104 i = 1,ntot
          write(12,*) x(i),fint(i)
  104 continue
      END

C======================================================================C

      SUBROUTINE gcmdln(ifile,ofile)
      IMPLICIT none
      CHARACTER*30 tag,ifile,ofile
      INTEGER i,narg

      ifile = 'numint.in'
      ofile = 'numint.out'

      narg = iargc()
      if (narg.eq.0) call usage()
      i = 0
  101 i = i + 1
      if (i.gt.narg) goto 102
          call getarg(i,tag)
          if (tag(1:2).eq.'-h'.or.tag(1:2).eq.'-?') then
              call usage()
          elseif (tag(1:2).eq.'-f') then
              i = i + 1
              call getarg(i,ifile)
          elseif (tag(1:2).eq.'-o') then
              i = i + 1
              call getarg(i,ofile)
          endif
          goto 101
  102 continue

      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE usage()
      IMPLICIT none
      write(*,*)
     1 'Usage: numint -f [IFILE] -o [OFILE] -h/-?'
      write(*,*)
     1 '  -f [IFILE] (Default numint.in) name of input file'
      write(*,*)
     1 '  -o [OFILE] (Default numint.out) name of output file'
      write(*,*)
     1 '  -h/-?    Online help (this screen)'
      write(*,*)
      write(*,*)
     1 'numint numerically integrates a time series using trapezoidal',
     2 ' integration'
      write(*,*)
      write(*,*)
     1 '  Input files'
      write(*,*)
     1 '    numint.in: data series'
      write(*,*)
     1 '      t1 f(t1)'
      write(*,*)
     1 '         :'
      write(*,*)
      write(*,*)
     1 '  Output files'
      write(*,*)
     1 '    numint.out: integral of data series (starting from first',
     2                   ' entry)'
      write(*,*)
     1 '      t1 F(t1)'
      write(*,*)
     1 '         :'
      STOP
      END
