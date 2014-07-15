      PROGRAM polyfit
      IMPLICIT none
      CHARACTER*30 ifile,ofile
      LOGICAL ex
      INTEGER NMAX,nobs
      PARAMETER (NMAX=40000)
      REAL*8 xobs(NMAX),yobs(NMAX)

      INTEGER ORDMAX,ord
      PARAMETER (ORDMAX=10)
      REAL*8 coeff(ORDMAX),xpre,ypre,dx,xconst,yconst
      INTEGER i,j,npre

C----
C Parse command line and check for input file
C----
      call gcmdln(ifile,ofile,ord)
      inquire(file=ifile,EXIST=ex)
      if (.not.ex) then
          print *,'No input file: ',ifile
          call usage()
      endif

C----
C Extract observed data from input time series
C----
      open(unit=11,file=ifile,status='old')
      i = 1
  101 read(11,*,end=102) xobs(i),yobs(i)
          i = i + 1
          goto 101
  102 continue
      nobs = i - 1

C----
C Compute polynomial coefficients
C----
      call polylsq(coeff,xobs,yobs,nobs,ord)

C----
C Write predicted results to ofile
C----
      open(unit=12,file=ofile,status='unknown')
      do 104 i = 1,nobs
          ypre = 0.0d0
          do 103 j = 1,ord+1
              ypre = ypre + coeff(j)*xobs(i)**(j-1)
  103     continue
          write (12,*) xobs(i),ypre
  104 continue

      END

C======================================================================C

      SUBROUTINE gcmdln(ifile,ofile,ord)
      
      IMPLICIT NONE
      CHARACTER*30 ifile,ofile,tag
      INTEGER i,narg,ord
      
      ifile = 'polyfit.in'
      ofile = 'polyfit.out'
      ord = 1
      
      narg = iargc()
      if (narg.eq.0) call usage()
      i = 0
  101 i = i + 1
      if (i.gt.narg) goto 102
          call getarg(i,tag)
          if (tag(1:2).eq.'-h'.or.tag(1:2).eq.'-?') then
              call usage()
          elseif (tag(1:4).eq.'-ord') then
              i = i + 1
              call getarg(i,tag)
              read(tag,'(BN,I3)') ord
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
     1 'Usage: polyfit -f [IFILE] -o [OFILE] -ord [ORDER]'
      write(*,*)
     1 '  -f [IFILE] (Default polyfit.in) name of input file'
      write(*,*)
     1 '  -o [OFILE] (Default polyfit.out) name of output file'
      write(*,*)
     1 '  -ord [ORDER] (Default 1) order of polynomial'
      write(*,*)
      write(*,*)
     1 'polyfit computes a least squares fit through a dataset with a'
      write(*,*)
     1 '    polynomial of order [ORDER]'
      STOP
      END




