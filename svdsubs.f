C----------------------------------------------------------------------C
C Collection of classical inversion subroutines from Numerical Recipes
C----------------------------------------------------------------------C

      SUBROUTINE dsvdfit(a,y,ndata,mpars,np,mp,xhat,chisq,TOL)
C----
C Modified from Numerical Recipes "svdfit" (no weighting by stdevs)
C
C - Model matrix a(1:ndata,1:mpars)
C - Data points y(1:ndata)
C - Solution vector xhat(1:mpars)
C - Arrays u(1:ndata,1:mpars), v(1:mdata,1:mpars), w(1:mpars) provide 
C     workspace on input; on output they define the singular value 
C     decomposition
C - np, mp are the physical dimensions of the matrices a,y,xhat,u, v, w
C     - np = physical space for observations
C     - mp = physical space for parameters
C     - np >= ndata, and mp >= mpars.
C
C Singular value decomposition of (ndata-by-mpars) matrix. Program 
C returns values for mpars fit parameters in xhat, and chi-squared, 
C chisq.
C
C Uses subroutines dsvbksb, dsvdcmp
C----
      IMPLICIT none
      INTEGER ndata,mpars,np,mp
      REAL*8 a(np,mp),y(np),u(np,mp),v(mp,mp),w(mp),chisq,xhat(mp)
      REAL*8 TOL
      INTEGER j
      REAL*8 thresh,wmax

C----
C Singular value decomposition. Truncate small singular values.
C----
      call dsvdcmp(a,ndata,mpars,np,mp,u,w,v)
      wmax = 1.0d0
      do 11 j = 1,mpars
          if (w(j).gt.wmax) wmax = w(j)
   11 continue
      thresh = TOL*wmax
      do 12 j = 1,mpars
          if (w(j).lt.thresh) w(j) = 0.0d0
   12 continue
C----
C Get parameter solutions in xhat
C----
      call dsvbksb(u,w,v,ndata,mpars,np,mp,y,xhat)
C----
C Evaluate chi-squared
C----
      chisq = 0.0d0
C      do 14 i = 1,ndata
C          call funcs(x(i),afunc,ma)
C          sum = 0.0d0
C          do 13 j = 1,ma
C              sum = sum + a(j)*afunc(j)
C   13     continue
C          chisq = chisq + ((y(i)-sum)/sig(i))**2
C   14 continue
C      print *,xhat
C      do 101 i = 1,mpars
C          print *,'xhat = ',xhat(i)
C  101 continue

      RETURN
      END
 
C----------------------------------------------------------------------C

      SUBROUTINE dsvbksb(u,w,v,m,n,mp,np,b,x)
C----
C Solves
C          A*X = B
C for a vector x, where A is specified by the arrays u, w, v, as 
C returned by dsvdcmp:
C          x = V*(W^-1)*U'*B
C m and n are the logical dimensions of a, and will be equal for square
C matrices. mp and np are the physical dimensions of a. b(1:m) is the
C input right-hand side. x(1:n) is the output solution vector. No input
C quantities are destroyed, so the routine may be called with different
C b's.
C----
      IMPLICIT none
      INTEGER m,mp,n,np,NMAX
      REAL*8 b(mp),u(mp,np),v(np,np),w(np),x(np)
      PARAMETER (NMAX=1000) ! Maximum anticipated value of n
      INTEGER i,j,jj
      REAL*8 s,tmp(NMAX)

C Calculate U'*B
      do 12 j = 1,n
          s = 0.0d0
          if (w(j).ne.0.0d0) then
              do 11 i = 1,m
                  s = s + u(i,j)*b(i)
   11         continue
C Divide by wj
              s = s/w(j)
          endif
          tmp(j) = s
   12 continue
C Matrix multiply by V
      do 14 j = 1,n
          s = 0.0d0
          do 13 jj = 1,n
              s = s + v(j,jj)*tmp(jj)
   13     continue
          x(j) = s
   14 continue

      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE dsvdcmp(a,m,n,mp,np,u,w,v)
C----
C Given a matrix a(1:m,1:n), with physical dimensions mp by np, this
C routine computes its singular value decomposition:
C            A = U*W*V'
C The matrix u(1:m,1:n) replaces a on output. The diagonal matrix of
C singular values W is output as a vector w(1:n). The matrix V (NOT
C its transpose, V') is output as v(1:n,1:n).
C
C Modified to retain a, and output u separately.
C
C Uses subroutine dpythag.
C----
      IMPLICIT none
      INTEGER m,mp,n,np,NMAX
      REAL*8 a(mp,np),asave(mp,np),u(mp,np),v(np,np),w(np)
      PARAMETER (NMAX=1000) ! Maximum anticipated value of n
      INTEGER i,its,j,jj,k,l,nm
      REAL*8 anorm,c,f,g,h,s,scale,x,y,z,rv1(NMAX),dpythag

      asave = a

C Householder reduction to diagonal form
      g = 0.0d0
      scale = 0.0d0
      anorm = 0.0d0
      do 25 i = 1,n
          l = i + 1
          rv1(i) = scale*g
          g = 0.0d0
          s = 0.0d0
          scale = 0.0d0
          if (i.le.m) then
              do 11 k = i,m
                  scale = scale + dabs(a(k,i))
   11         continue
              if (scale.ne.0.0d0) then
                  do 12 k = i,m
                      a(k,i) = a(k,i)/scale
                      s = s + a(k,i)*a(k,i)
   12             continue
                  f = a(i,i)
                  g = -dsign(dsqrt(s),f)
                  h = f*g - s
                  a(i,i) = f-g
                  do 15 j = l,n
                      s = 0.0d0
                      do 13 k = i,m
                          s = s + a(k,i)*a(k,j)
   13                 continue
                      f = s/h
                      do 14 k = i,m
                          a(k,j) = a(k,j) + f*a(k,i)
   14                 continue
   15             continue
                  do 16 k = i,m
                      a(k,i) = scale*a(k,i)
   16             continue
              endif
          endif
          w(i) = scale *g
          g = 0.0d0
          s = 0.0d0
          scale = 0.0d0
          if ((i.le.m).and.(i.ne.n)) then
              do 17 k = l,n
                  scale = scale + dabs(a(i,k))
   17         continue
              if (scale.ne.0.0d0) then
                  do 18 k = l,n
                      a(i,k) = a(i,k)/scale
                      s = s + a(i,k)*a(i,k)
   18             continue
                  f = a(i,l)
                  g = -dsign(dsqrt(s),f)
                  h = f*g - s
                  a(i,l) = f - g
                  do 19 k = l,n
                      rv1(k) = a(i,k)/h
   19             continue
                  do 23 j = l,m
                      s = 0.0d0
                      do 21 k = l,n
                          s = s + a(j,k)*a(i,k)
   21                 continue
                      do 22 k = l,n
                          a(j,k) = a(j,k) + s*rv1(k)
   22                 continue
   23             continue
                  do 24 k = l,n
                      a(i,k) = scale*a(i,k)
   24             continue
              endif
          endif
          anorm = max(anorm,(dabs(w(i))+dabs(rv1(i))))
   25 continue
C Accumulation of right-hand transformations
      do 32 i = n,1,-1
          if (i.lt.n) then
              if (g.ne.0.0d0) then
C Double division to avoid underflow
                  do 26 j = l,n
                      v(j,i) = (a(i,j)/a(i,l))/g
   26             continue
                  do 29 j = l,n
                      s = 0.0d0
                      do 27 k = l,n
                          s = s + a(i,k)*v(k,j)
   27                 continue
                      do 28 k=l,n
                          v(k,j) = v(k,j) + s*v(k,i)
   28                 continue
   29             continue
              endif
              do 31 j = l,n
                  v(i,j) = 0.0d0
                  v(j,i) = 0.0d0
   31         continue
          endif
          v(i,i) = 1.0d0
          g = rv1(i)
          l = i
   32 continue
C Accumulation of left-hand transformations
      do 39 i = min(m,n),1,-1
          l = i+1
          g = w(i)
          do 33 j = l,n
              a(i,j) = 0.0d0
   33     continue
          if (g.ne.0.0d0) then
              g = 1.0d0/g
              do 36 j = l,n
                  s = 0.0
                  do 34 k = l,m
                      s = s + a(k,i)*a(k,j)
   34             continue
                  f = (s/a(i,i))*g
                  do 35 k = i,m
                      a(k,j) = a(k,j) + f*a(k,i)
   35             continue
   36         continue
              do 37 j = i,m
                  a(j,i) = a(j,i)*g
   37         continue
          else
              do 38 j = i,m
                  a(j,i) = 0.0d0
   38         continue
          endif
          a(i,i) = a(i,i) + 1.0d0
   39 continue
C Diagonalization of the bidiagonal form: loop over singular values
C and over allowed iterations
      do 49 k = n,1,-1
          do 48 its = 1,30
C Test for splitting
              do 41 l = k,1,-1
                  nm = l-1
                  if ((dabs(rv1(l))+anorm).eq.anorm)  goto 2
                  if ((dabs(w(nm))+anorm).eq.anorm)  goto 1
   41         continue
    1         c = 0.0d0
              s = 1.0d0
              do 43 i=l,k
                  f = s*rv1(i)
                  rv1(i) = c*rv1(i)
                  if ((dabs(f)+anorm).eq.anorm) goto 2
                  g = w(i)
                  h = dpythag(f,g)
                  w(i) = h
                  h = 1.0d0/h
                  c = (g*h)
                  s = -(f*h)
                  do 42 j = 1,m
                      y = a(j,nm)
                      z = a(j,i)
                      a(j,nm) = (y*c) + (z*s)
                      a(j,i) = -(y*s) + (z*c)
   42             continue
   43         continue
    2         z = w(k)
C Convergence
              if (l.eq.k) then
C Singular value made non-negative
                  if (z.lt.0.0d0) then
                      w(k) = -z
                      do 44 j = 1,n
                          v(j,k) = -v(j,k)
   44                 continue
                  endif
                  goto 3
              endif
              if (its.eq.30) stop 'no convergence in dsvdcmp'
              x = w(l)
              nm = k-1
              y = w(nm)
              g = rv1(nm)
              h = rv1(k)
              f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y)
              g = dpythag(f,1.0d0)
              f = ((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
              c = 1.0d0
              s = 1.0d0
              do 47 j = l,nm
                  i = j+1
                  g = rv1(i)
                  y = w(i)
                  h = s*g
                  g = c*g
                  z = dpythag(f,h)
                  rv1(j) = z
                  c = f/z
                  s = h/z
                  f = (x*c)+(g*s)
                  g =-(x*s)+(g*c)
                  h = y*s
                  y = y*c
                  do 45 jj = 1,n
                      x = v(jj,j)
                      z = v(jj,i)
                      v(jj,j) = (x*c) + (z*s)
                      v(jj,i) = -(x*s) + (z*c)
   45             continue
                  z = dpythag(f,h)
                  w(j) = z
                  if (z.ne.0.0d0) then
                      z = 1.0d0/z
                      c = f*z
                      s = h*z
                  endif
                  f = (c*g) + (s*y)
                  x = -(s*g) + (c*y)
                  do 46 jj = 1,m
                      y = a(jj,j)
                      z = a(jj,i)
                      a(jj,j) = (y*c) + (z*s)
                      a(jj,i) = -(y*s) + (z*c)
   46             continue
   47         continue
              rv1(l) = 0.0d0
              rv1(k) = f
              w(k) = x
   48     continue
    3     continue
   49 continue
      
      u = a
      a = asave

      RETURN
      END

C----------------------------------------------------------------------C

      FUNCTION dpythag(a,b)
      IMPLICIT none
      REAL*8 a,b,dpythag
      REAL*8 absa,absb

      absa = dabs(a)
      absb = dabs(b)
      if (absa.gt.absb) then
          dpythag = absa*dsqrt(1.0d0+(absb/absa)**2)
      else
          if (absb.eq.0.0d0) then
              dpythag = 0.0d0
          else
              dpythag = absb*sqrt(1.0d0+(absa/absb)**2)
          endif
      endif

      RETURN
      END
