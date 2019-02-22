      PROGRAM main
      IMPLICIT none
      INTEGER NMAX
      PARAMETER (NMAX=100000)
      REAL*8 ts(NMAX,3)
      CHARACTER*200 line,ifile,tfile
      INTEGER i,j
      INTEGER TMAX,nwind
      PARAMETER (TMAX=5000)
      REAL*8 twind(TMAX,3)
      CHARACTER*1 ttype

#ifndef USELAPACK
      write(0,*) 'polyfit_special: compile with LAPACK to use'
      stop
#endif

C      call gcmdln(ifile,wfile)
C      call readwin(wfile,twind,nwind,ttype)
C
C      open(unit=11,file=ifile,status='old')
C  103 read(11,*,end=102) np,stlo,stla
C          do 101 i = 1,np
C          read(11,*) (ts(i,j),j=1,3)
C  101     continue
C          call calcdefm(ts,np,twind,nwind,ttype,defm) ! Return values at each iwind in defm
C          call writedefm(twind,nwind,defm,stlo,stla)
C          goto 103
C  102 continue
C
C----
C Read output times, averaging window, displacement/velocity
C----
      tfile = 'times.in'
      open(unit=12,file=tfile,status='old')
      i = 1
  104 read(12,*,end=103) twind(i,1),twind(i,2),twind(i,3),ttype
          i = i + 1
          goto 104
  103 continue
      close(12)
      nwind = i - 1
      print *,'WINDOWS'
      do 105 i = 1,nwind
          print *,i,twind(i,1),twind(i,2),twind(i,3),ttype
  105 continue

C----
C Read time series
C----
      ifile = 'polyfit2.in'
      open(unit=11,file=ifile,status='old')
      ts(1,1) = -1.0d99
      i = 1
  101 read(11,'(A)',end=102) line
          j = index(line,'>')
          if (j.ne.0) then
              i = i - 1
              if (ts(1,1).gt.-1.0d98) call fitit(ts,i,twind,nwind)
              i = 0
          else
              read(line,*) ts(i,1),ts(i,2),ts(i,3)
          endif
          i = i + 1
          goto 101
  102 continue
      close(11)
      i = i - 1
      call fitit(ts,i,twind,nwind)
      END

C----------------------------------------------------------------------C

      SUBROUTINE fitit(ts,nobs,twind,nwind)
      IMPLICIT none
      INTEGER NMAX,nobs
      PARAMETER (NMAX=100000)
      REAL*8 ts(NMAX,3),tss(NMAX,3)
      INTEGER TMAX,nwind
      PARAMETER (TMAX=5000)
      REAL*8 twind(TMAX,3)
      INTEGER i,j,k
      print *,'TIME SERIES DATA'
      do 104 i = 1,nobs
          print *,i,ts(i,1),ts(i,2),ts(i,3)
  104 continue
      do 101 i = 1,nwind
          print *,'window ',i
          k = 0
          do 102 j = 1,nobs
              if (twind(i,2).le.ts(j,1).and.ts(j,1).le.twind(i,3)) then
                  k = k + 1
                  tss(k,1) = ts(j,1)
                  tss(k,2) = ts(j,2)
                  tss(k,3) = ts(j,3)
              endif
  102     continue
          do 103 j = 1,k
              print *,j,tss(j,1)
  103     continue
          if (k.gt.1) then
              call lsqsolve(tss,k)
          endif
  101 continue
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE lsqsolve(ts,nobs)
      IMPLICIT none
      INTEGER NMAX,nobs
      PARAMETER (NMAX=100000)
      REAL*8 ts(NMAX,3)
      INTEGER i,j
      INTEGER nrhs,lwmax
      PARAMETER (nrhs=1,lwmax=100000)
      CHARACTER*1 trans
      INTEGER m,n,lda,ldb,lwork,info
      REAL*8  a(nobs,2),b(nobs,nrhs),work(lwmax)

      trans = 'N' ! A matrix has form (nobs x nparam)
      m = nobs    ! Number of rows in matrix A
      n = 2       ! Number of columns in matrix A
      lda = nobs  ! Leading dimension of A  lda >= max(1,m)
      ldb = nobs  ! Leading dimension of b  ldb >= max(1,m,n)

C----
C Load model matrix (A) and observation vector (b)
C----
      do 102 i = 1,nobs
          do 101 j = 1,2
              a(i,j) = ts(i,1)**(j-1)
  101     continue
  102 continue
      do 103 i = 1,nobs
          b(i,1) = ts(i,2)
  103 continue
      lwork = -1
#ifdef USELAPACK
      call dgels(trans,m,n,nrhs,a,lda,b,ldb,work,lwork,info)
#endif
      lwork = min(lwmax,int(work(1)))
#ifdef USELAPACK
      call dgels(trans,m,n,nrhs,a,lda,b,ldb,work,lwork,info)
#endif
      do 104 i = 1,2
          print *,b(i,1)
  104 continue
C
      do 107 i = 1,nobs
          do 108 j = 1,2
              a(i,j) = ts(i,1)**(j-1)
  108     continue
  107 continue
      do 105 i = 1,nobs
          b(i,1) = ts(i,3)
  105 continue
      lwork = -1
#ifdef USELAPACK
      call dgels(trans,m,n,nrhs,a,lda,b,ldb,work,lwork,info)
#endif
      lwork = min(lwmax,int(work(1)))
#ifdef USELAPACK
      call dgels(trans,m,n,nrhs,a,lda,b,ldb,work,lwork,info)
#endif
      do 106 i = 1,2
          print *,b(i,1)
  106 continue
      RETURN
      END
