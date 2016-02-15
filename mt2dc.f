      PROGRAM main
      IMPLICIT none
      INTEGER   n,lda,lwmax
      PARAMETER (n=3,lda=n,lwmax=250)
      INTEGER   info,lwork
      REAL*8    a(lda,n),w(n),work(lwmax),m(6)
      REAL*8 p(3),t(3),slip(3),normal(3),str,dip,rak
      INTEGER i
C     Moment tensor components (r=z, t=e, p=n):
C                        mrr  mtt  mpp  mrt  mrp  mtp
  101 read (*,*,END=102) m(1),m(2),m(3),m(4),m(5),m(6)
          a(3,3) =  m(1)     ! mrr = mzz
          a(2,2) =  m(2)     ! mtt = myy
          a(1,1) =  m(3)     ! mpp = mxx
          a(2,3) = -m(4)     ! mrt = myz
          a(3,2) = -m(4)
          a(1,3) =  m(5)     ! mrp = mxz
          a(3,1) =  m(5)
          a(1,2) = -m(6)     ! mtp = mxy
          a(1,2) = -m(6)
C         Query the optimal workspace
          lwork = -1
          call dsyev('Vectors','Upper',n,a,lda,w,work,lwork,info)
          lwork = min(lwmax,int(work(1)))
C         Solve eigenproblem
          call dsyev('Vectors','Upper',n,a,lda,w,work,lwork,info)
C         Compute best-fitting focal mechanism
          do 103 i = 1,3
              t(i) = a(i,1)
              p(i) = a(i,3)
  103     continue
          call pt2un(p,t,slip,normal)
          call un2sdr(slip,normal,str,dip,rak)
          write(*,1001),str,dip,rak
          goto 101
  102 continue
 1001 format(3F8.1)
      END

C----------------------------------------------------------------------C

      SUBROUTINE pt2un(p,t,slip,normal)
C----
C Compute unit normal and slip vectors from P and T axes
C----
      IMPLICIT none
      REAL*8 sqrt2
      PARAMETER (sqrt2=dsqrt(2.0d0))
      REAL*8 p(3),t(3),slip(3),normal(3)
      INTEGER i
      do 101 i = 1,3
          slip(i)   = (t(i)+p(i))/sqrt2
          normal(i) = (t(i)-p(i))/sqrt2
  101 continue
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE un2sdr(slip,normal,s,d,r)
C----
C Compute strike, dip, rake from slip and normal to plane
C----
      IMPLICIT none
      REAL*8 pi,d2r,r2d
      PARAMETER (pi=4.0d0*datan(1.0d0),d2r=pi/1.8d2,r2d=1.8d2/pi)
      REAL*8 slip(3),normal(3),s,d,r
      REAL*8 hor,st(3),ns,dum
      INTEGER i
      do 101 i = 1,3
          slip(i)   = -slip(i)
  101 continue
      if (normal(3).lt.0.0d0) then
          do 102 i = 1,3
              slip(i)   = -slip(i)
              normal(i) = -normal(i)
  102     continue
      endif
      s = datan2(normal(1),normal(2))/d2r - 90.0d0 ! strike
      if (s.lt.0.0d0) s = s + 360.0d0
      hor = dsqrt(normal(1)*normal(1)+normal(2)*normal(2))
      d = datan2(hor,normal(3))/d2r                ! dip
      st(1) = dsin(s*d2r) ! components of strike
      st(2) = dcos(s*d2r) ! components of strike
      st(3) = 0.0d0       ! components of strike
      r = dacos(st(1)*slip(1)+st(2)*slip(2))*r2d
C Check rake is pos or neg
      dum = st(1)*slip(2)-st(2)*slip(1)
      if (dum.lt.0.0d0) then
          r = -r
      endif
      RETURN
      END



