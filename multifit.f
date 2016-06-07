      PROGRAM main
      IMPLICIT none
C Constants
      !REAL*8 pi
      !PARAMETER (pi=4.0d0*datan(1.0d0))
C Command line inputs
      CHARACTER*120 ifile
      INTEGER ipoly,isin,iexp
      REAL*8 per(3)
C Observation variables
      INTEGER nobs,OBSMAX
      PARAMETER (OBSMAX=50000)
      REAL*8 obs(OBSMAX,2)

      REAL*8 soln(12)
      !REAL*8 A(OBSMAX,4),b(OBSMAX,1),x(OBSMAX,1)
      !REAL*8 afreq,period
      !INTEGER dimx,dimy,m,n

C----
C Parse command line
C----
      call gcmdln(ifile,ipoly,isin,per,iexp)

C----
C Read in observations
C----
      call readobs(ifile,obs,nobs,OBSMAX)

C----
C Solve least squares problem
C----
      call lstsqr(soln,obs,nobs,OBSMAX,ipoly,isin,per,iexp)

C----
C Print results
C----
C Polynomial coefficients
      if (soln(1).gt.-1.0d98) write(*,1001) 'POLY 0',soln(1)
      if (soln(2).gt.-1.0d98) write(*,1001) 'POLY 1',soln(2)
      if (soln(3).gt.-1.0d98) write(*,1001) 'POLY 2',soln(3)
      if (soln(4).gt.-1.0d98) write(*,1001) 'POLY 3',soln(4)
C Sinusoid amplitude and phase
      if (soln(5).gt.-1.0d98) then
          write(*,1002) 'SINE 1AMP',soln(5)
          write(*,1003) 'SINE 1PHA',soln(6)
      endif
      if (soln(7).gt.-1.0d98) then
          write(*,1002) 'SINE 2AMP',soln(7)
          write(*,1003) 'SINE 2PHA',soln(8)
      endif
      if (soln(9).gt.-1.0d98) then
          write(*,1002) 'SINE 3AMP',soln(9)
          write(*,1003) 'SINE 3PHA',soln(10)
      endif
C Exponential amplitude and constant
      if (soln(11).gt.-1.0d98) then
          write(*,*) 'EXP AMP',soln(11)
          write(*,*) 'EXP CONST',soln(12)
      endif

 1001 format(A,1PE12.3)
 1002 format(A,1PE12.3)
 1003 format(A,F10.4)
      END

C----------------------------------------------------------------------C

      SUBROUTINE readobs(ifile,obs,nobs,OBSMAX)
      IMPLICIT none
      CHARACTER*120 ifile
      INTEGER input
      LOGICAL ex
      INTEGER i,nobs,OBSMAX
      REAL*8 obs(OBSMAX,2)
      if (ifile.eq.'stdin') then
          input = 5
      else
          inquire(file=ifile,exist=ex)
          if (.not.ex) then
              call usage('!! Error: no input file '//trim(ifile))
          endif
          input = 11
          open(unit=input,file=ifile,status='old')
      endif
      i = 1
  101 read(input,*,end=102) obs(i,1),obs(i,2)
          i = i + 1
          if (i.gt.OBSMAX) then
              call usage('!! Error: too many input points')
          endif
          goto 101
  102 continue
      nobs = i - 1
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE lstsqr(soln,obs,nobs,OBSMAX,ipoly,isin,per,iexp)
      IMPLICIT none
      REAL*8 pi
      PARAMETER (pi=4.0d0*datan(1.0d0))
      REAL*8 soln(12)
      INTEGER nobs,OBSMAX
      REAL*8 obs(OBSMAX,2)
      INTEGER ipoly,isin,iexp
      REAL*8 per(3),afreq(3)
      INTEGER i,j
      INTEGER dimx,dimy,n,m
      REAL*8 A(2*OBSMAX,12),b(2*OBSMAX,1),x(2*OBSMAX,1)
      REAL*8 C1,C2
C----
C Compute angular frequency from period
C----
      do 101 i = 1,isin
          afreq(i) = 2.0d0*pi/per(i)
  101 continue

C----
C Load model matrix, A, and observation vector, b
C----
      do 201 i = 1,nobs
          ! Polynomial
          A(i,1) = 1.0d0
          do 202 j = 1,ipoly
              A(i,j+1) = obs(i,1)**ipoly
  202     continue
          ! Sinusoid
          do 203 j = 1,isin
              A(i,ipoly+1+j) = dsin(afreq(j)*obs(i,1))
              A(i,ipoly+1+isin+j) = dcos(afreq(j)*obs(i,1))
  203     continue
          ! Exponential
          if (iexp.eq.1) then
              A(nobs+i,ipoly+1+2*isin+1) = 1.0d0
              A(nobs+i,ipoly+1+2*isin+2) = obs(i,1)
          endif
          ! Observation matrix
          b(i,1) = obs(i,2)
          if (iexp.eq.1) then
              b(nobs+i,1) = dlog(obs(i,2))
          endif
  201 continue

C----
C Solve least squares problem
C----
      dimx = 2*OBSMAX
      dimy = 12
      n    = nobs + nobs*iexp
      m    = 1+ipoly + 2*isin + 2*iexp
      call lsqsolve(A,b,x,dimx,dimy,n,m)

C----
C Save as clean solution vector
C----
      do 206 i = 1,12
          soln(i) = -1.0d99
  206 continue
      do 207 i = 1,ipoly+1
          soln(i) = x(i,1)
  207 continue
      do 208 i = 1,isin
          C1 = x(ipoly+1+i,1)
          C2 = x(ipoly+1+isin+i,1)
          soln(4+2*i-1) = dsqrt(C1*C1+C2*C2)
          soln(4+2*i)   = datan2(C2,C1)
  208 continue
      if (iexp.eq.1) then
          soln(11) = dexp(x(ipoly+1+2*isin+1,1))
          soln(12) = x(ipoly+1+2*isin+2,1)
      endif
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE lsqsolve(Ain,bin,xout,OBSMAX,FLTMAX,nobs,npar)
C----
C Solve generalized least squares problem Ax = b.
C On input, arrays Ain and bin are given maximum dimensions:
C     Ain(OBSMAX,FLTMAX)
C     bin(OBSMAX,1)
C To operate with LAPACK generalized inverse tools, convert arrays from
C maximum size to operating size:
C     Ain --> A(nobs,npar)
C     bin --> b(nobs,1)
C----
      IMPLICIT none
      INTEGER OBSMAX,FLTMAX
      REAL*8 Ain(OBSMAX,FLTMAX),bin(OBSMAX,1),xout(OBSMAX,1)
      INTEGER nobs,npar
      REAL*8 A(nobs,npar),b(nobs,1),x(nobs,1)
      INTEGER i,j
C Load existing arrays into arrays of correct dimensions
      do 432 i = 1,nobs
          b(i,1) = bin(i,1)
          do 431 j = 1,npar
              A(i,j) = Ain(i,j)
  431     continue
  432 continue
C Solve generalized least squares problem using LAPACK tools
      call solve(A,x,b,nobs,npar)
C Print results
      do 433 i = 1,npar
          xout(i,1) = x(i,1)
  433 continue
      RETURN
      END


C----------------------------------------------------------------------C

      SUBROUTINE solve(A,x,b,nobs,npar)
C----
C Determine the optimal workspace size for work array in dgels
C Solve least squares problem Ax = b for x using QR or LQ decomposition
C----
      IMPLICIT none
      INTEGER nobs,npar,nrhs,i
      REAL*8 A(nobs,npar),b(nobs,1),btmp(nobs,1),x(npar,1)
      INTEGER lwork
      INTEGER m,n,lda,ldb
      INTEGER info
      CHARACTER*1 trans
      REAL*8 work(100000)
      nrhs = 1
      trans = 'N'     ! A has form (nobs x npar), i.e. not transposed
      m = nobs        ! Number of rows in matrix A
      n = npar        ! Number of columns in matrix A
      lda = m         ! Leading dimension of A  lda >= max(1,m)
      ldb = max(m,n)  ! Leading dimension of b  ldb >= max(1,m,n)
C Compute optimal workspace
      lwork = -1
      call dgels(trans,m,n,nrhs,A,lda,b,ldb,work,lwork,info)
      lwork = int(work(1))
      !print *,'LWORK',lwork
C Copy observation vector, b to btmp (replaced in dgels)
      do 441 i = 1,nobs
          btmp(i,1) = b(i,1)
  441 continue
C Compute parameter array
      call dgels(trans,m,n,nrhs,a,lda,btmp,ldb,work,lwork,info)
      do 442 i = 1,npar
          x(i,1) = btmp(i,1)
  442 continue
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE gcmdln(ifile,ipoly,isin,per,iexp)
      IMPLICIT none
      CHARACTER*120 ifile
      INTEGER ipoly,isin,iexp
      REAL*8 per(3)
      INTEGER i,narg
      CHARACTER*80 tag
      ifile = 'none'
      ipoly = 0 ! Order of polynomial (default: compute mean)
      isin  = 0 ! Number of sinusoidal fits
      do 103 i = 1,3
          per(i) = 0.0d0
  103 continue
      iexp  = 0 ! Exponential fit flag
      narg = iargc()
      if (narg.eq.0) call usage('')
      i = 0
  101 i = i + 1
      if (i.gt.narg) goto 102
          call getarg(i,tag)
          if (tag(1:5).eq.'-poly') then
              i = i + 1
              call getarg(i,tag)
              read(tag,*) ipoly
          elseif (tag(1:4).eq.'-sin') then
              isin = isin + 1
              if (isin.gt.3) then
                  call usage('!! Error: max 3 different periods')
              endif
              i = i + 1
              call getarg(i,tag)
              read(tag,*) per(isin)
          elseif (tag(1:4).eq.'-exp') then
              iexp = 1
          elseif (tag(1:2).eq.'-f') then
              i = i + 1
              call getarg(i,ifile)
          else
              call usage('!! Error: No option '//tag)
          endif
          goto 101
  102 continue
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE usage(str)
      IMPLICIT none
      CHARACTER str*(*)
      if (str.ne.'') then
            write(0,*) trim(str)
            write(0,*)
      endif
      write(0,*)
     1 'Usage: multifit -f XYFILE [-poly ORD] [-sin PERIOD] [-exp]'
      write(0,*)
     1 '    -f XYFILE    Input x-y file'
      write(0,*)
     1 '    -poly ORD    Fit a polynomial (up to order 3)'
      write(0,*)
     1 '    -sin PER     Fit a sinusoid with period PER (can be ',
     2                           'repeated up to 3 times with ',
     2                                  'different periods)'
      write(0,*)
     1 '    -exp         Fit an exponential'
      STOP
      END

