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
      INTEGER NMAX,ORDMAX,nobs,ord
      PARAMETER (NMAX=40000,ORDMAX=10)
      REAL*8 xobs(NMAX),yobs(NMAX),coeff(ORDMAX)
      INTEGER nrhs,lwmax
      PARAMETER (nrhs=1,lwmax=100)
      CHARACTER*1 trans
      INTEGER m,n,lda,ldb,lwork,info,i,j
      REAL*8  a(nobs,ord+1),b(nobs,nrhs),work(lwmax)

      trans = 'N' ! A matrix has form (nobs x nparam)
      m = nobs    ! Number of rows in matrix
      n = ord+1   ! Number of columns in matrix
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

