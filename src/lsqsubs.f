      SUBROUTINE polylsq(coeff,xobs,yobs,nobs,ord)
C----
C Compute the best-fitting polynomial through a set of points
C consisting of nobs ordered pairs (xobs(i),yobs(i)). The polynomial
C is of order "ord" (lin = 1, quad = 2, ...), and the (ord+1)
C coefficients are returned in the variable "coeff."
C
C   y = coeff(ord+1)*x^ord + coeff(ord)*x^(ord-1) + ... + coeff(1)
C
C Solution to system of linear equations is computed using LAPACK
C subroutine "dgels." On my computer, the LAPACK libraries are
C   /Users/mherman/Research/lapack-3.5.0/
C        liblapack.a, libmglib.a, librefblas.a
C----
      IMPLICIT none
      INTEGER OBSMAX,ORDMAX,nrhs,lwmax
      PARAMETER (OBSMAX=100000,ORDMAX=20,nrhs=1,lwmax=100000)
      REAL*8 xobs(OBSMAX),yobs(OBSMAX),coeff(ORDMAX)
      CHARACTER*1 trans
      INTEGER nobs,ord,m,n,lda,ldb,lwork,info,i,j
      REAL*8  a(nobs,ord+1),b(nobs,nrhs),work(lwmax)

      trans = 'N' ! A matrix has form (nobs x nparam)
      m = nobs    ! Number of rows in matrix A
      n = ord+1   ! Number of columns in matrix A
      lda = nobs  ! Leading dimension of A  lda >= max(1,m)
      ldb = nobs  ! Leading dimension of b  ldb >= max(1,m,n)

C----
C Load model matrix (A) and observation vector (b)
C----
      do 102 i = 1,nobs
          do 101 j = 1,ord+1
              a(i,j) = xobs(i)**(j-1)
  101     continue
  102 continue
      do 103 i = 1,nobs
          b(i,1) = yobs(i)
  103 continue

C----
C Determine the optimal workspace size for the solution
C----
      lwork = -1 
      call dgels(trans,m,n,nrhs,a,lda,b,ldb,work,lwork,info)
      lwork = min(lwmax,int(work(1)))

C----
C Compute the polynomial coefficients
C----
      call dgels(trans,m,n,nrhs,a,lda,b,ldb,work,lwork,info)
      do 104 i = 1,ord+1
          coeff(i) = b(i,1)
  104 continue
      
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE polylsqconstr(coeff,xobs,yobs,nobs,ord,xconstr,yconstr,
     1                         nconstr)
C---- 
C Compute the best-fitting polynomial through a set of points
C consisting of nobs ordered pairs (xobs(i),yobs(i)), constrained
C to pass through nconstr points (xconstr(i),yconstr(i)). The
C polynomial is of order "ord" (lin = 1, quad = 2, ...), and the
C (ord+1) coefficients are returned in the variable "coeff."
C
C Solution to system of linear equations with equality constraints is
C computed using LAPACK subroutine dgglse.
C----
      IMPLICIT none
      INTEGER OBSMAX,ORDMAX,CONSTRMAX,lwmax
      PARAMETER (OBSMAX=100000,ORDMAX=20,CONSTRMAX=21,lwmax=100000)
      REAL*8 xobs(OBSMAX),yobs(OBSMAX),coeff(ORDMAX),
     1       xconstr(CONSTRMAX),yconstr(CONSTRMAX)
      INTEGER nobs,ord,nconstr,m,n,p,lda,ldb,lwork,info,i,j
      REAL*8 a(nobs,ord+1),b(nconstr,ord+1),c(nobs),d(nconstr),
     1       work(lwmax),x(ord+1)

      if (nconstr.gt.ord+1) then
          write(*,*) 'Error: nconstr > ord+1'
          write(*,*) 'The number of equality constraints cannot be ',
     1               'larger than the order of the polynomial plus one.'
          STOP
      endif

      m = nobs      ! Number of rows in matrix A
      n = ord+1     ! Number of columns in matrices A and B
      p = nconstr   ! Number of rows in matrix B
      lda = nobs    ! Leading dimension of A  lda >= max(1,m)
      ldb = nconstr ! Leading dimension of N  ldb >= max(1,p)

C---
C Load model matrix (A), constraint matrix (B), observation vector (c),
C and the constraint RHS vector (d).
C---  
      do 101 j = 1,ord+1
          do 102 i = 1,nobs
              a(i,j) = xobs(i)**(j-1)
  102     continue
          do 103 i = 1,nconstr
              b(i,j) = xconstr(i)**(j-1)
  103     continue
  101 continue
      do 104 i = 1,nobs
          c(i) = yobs(i)
  104 continue
      do 105 i = 1,nconstr
          d(i) = yconstr(i)
  105 continue

C----
C Determine the optimal workspace size for the solution
C----
      lwork = -1
      call dgglse(m,n,p,a,lda,b,ldb,c,d,x,work,lwork,info)
      lwork = min(lwmax,int(work(1)))

C----
C Compute the polynomial coefficients
C----
      call dgglse(m,n,p,a,lda,b,ldb,c,d,x,work,lwork,info)
      do 301 i = 1,ord+1
          coeff(i) = x(i)
  301 continue

      RETURN
      END
