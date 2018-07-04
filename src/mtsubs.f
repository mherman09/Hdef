      SUBROUTINE sdr2sdr(line,soln)
C----
C Compute the second nodal plane of a double couple source
C algebraically.
C----
      IMPLICIT none
      CHARACTER*180 line,soln
      REAL*8 pi,d2r,r2d
      PARAMETER (pi=4.0d0*datan(1.0d0),d2r=pi/1.8d2,r2d=1.8d2/pi)
      REAL*8 str1,dip1,rak1,str,dip,rak,str2,dip2,rak2
      REAL*8 xn1,yn1,zn1,xr,yr,zr,xs2,ys2,zs2,dum1
C Get input strike, dip, rake in degrees; convert to radians
      read(line,*) str1,dip1,rak1
      str = d2r*str1
      dip = d2r*dip1
      rak = d2r*rak1
C Normal vector to input plane
      xn1 =  dsin(dip)*dcos(str)
      yn1 = -dsin(dip)*dsin(str)
      zn1 =  dcos(dip)
C Rake vector of input plane
      xr =  dcos(rak)*dsin(str) - dcos(dip)*dsin(rak)*dcos(str)
      yr =  dcos(rak)*dcos(str) + dcos(dip)*dsin(rak)*dsin(str)
      zr =  dsin(rak)*dsin(dip)
      ! If vector points downward, flip it
      if (zr.lt.0.0d0) then
          xr = -xr
          yr = -yr
          zr = -zr
          xn1 = -xn1
          yn1 = -yn1
          zn1 = -zn1
      endif
C Compute new strike
      str2 = r2d*datan2(xr,yr) - 90.0d0
      if (str2.gt.360.0d0) str2 = str2 - 360.0d0
      if (str2.lt.0.0d0) str2 = str2 + 360.0d0
C Compute new dip
      dip2 = r2d*datan(dsqrt(xr*xr+yr*yr)/zr)
C Compute new rake
      xs2 = dsin(str2*d2r)
      ys2 = dcos(str2*d2r)
      zs2 = 0.0d0
      rak2 = r2d*dacos(xn1*xs2+yn1*ys2+zn1*zs2)
C Check rake is pos or neg
      dum1 = xs2*yn1 - ys2*xn1
      if (dum1.lt.0.0d0) then
          rak2 = -rak2
      endif
C Write strike, dip, and rake in degrees
      write(soln,9999) str2,dip2,rak2
 9999 format(3F6.0)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE sdr2mij(line,soln)
C----
C Compute the (unit) moment tensor from double couple strike, dip, and
C rake angles.
C----
      IMPLICIT none
      CHARACTER*180 line,soln
      REAL*8 pi,d2r
      PARAMETER (pi=4.0d0*datan(1.0d0),d2r=pi/1.8d2)
      REAL*8 str1,dip1,rak1,str,dip,rak
      REAL*8 cs,ss,s2s,c2s,cd,sd,s2d,c2d,cr,sr
      REAL*8 mrr,mtt,mpp,mrt,mrp,mtp
C Get input strike, dip, rake in degrees; convert to radians
      read(line,*) str1,dip1,rak1
      str = d2r*str1
      dip = d2r*dip1
      rak = d2r*rak1
C Derived trigonometric variables
      cs  = dcos(str)
      ss  = dsin(str)
      s2s = dsin(2.0d0*str)
      c2s = dcos(2.0d0*str)
      cd  = dcos(dip)
      sd  = dsin(dip)
      s2d = dsin(2.0d0*dip)
      c2d = dcos(2.0d0*dip)
      cr  = dcos(rak)
      sr  = dsin(rak)
C Moment tensor components in GCMT convention (r,t,p)
      mrr =  (            s2d*sr      )
      mtt = -(sd*cr*s2s + s2d*sr*ss*ss)
      mpp =  (sd*cr*s2s - s2d*sr*cs*cs)
      mrt = -(cd*cr*cs  + c2d*sr*ss   )
      mrp =  (cd*cr*ss  - c2d*sr*cs   )
      mtp = -(sd*cr*c2s + s2d*sr*ss*cs)
C Write moment tensor components
      write(soln,9999) mrr,mtt,mpp,mrt,mrp,mtp
 9999 format(6(1PE14.6))
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE sdr2ter(line,soln)
      IMPLICIT none
      CHARACTER*180 line,soln
      REAL*8 p(3),n(3),t(3),fth,fss,fno
      call sdr2pnt(line,soln)
      read(soln,*) p(1),p(2),p(3),n(1),n(2),n(3),t(1),t(2),t(3)
      fno = datan(p(3)/dsqrt(p(1)*p(1)+p(2)*p(2)))
      fss = datan(n(3)/dsqrt(n(1)*n(1)+n(2)*n(2)))
      fth = datan(t(3)/dsqrt(t(1)*t(1)+t(2)*t(2)))
      fno = dsin(fno)*dsin(fno)
      fss = dsin(fss)*dsin(fss)
      fth = dsin(fth)*dsin(fth)
      write(soln,1001) fth,fss,fno
 1001 format(3F7.3)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE sdr2sv(line,soln)
      IMPLICIT none
      REAL*8 pi,r2d,d2r
      PARAMETER (pi=4.0d0*datan(1.0d0),d2r=pi/1.8d2,r2d=1.8d2/pi)
      CHARACTER*180 line,soln
      REAL*8 s1,d1,r1,x1,y1,z1,n1,e1
      read(line,*) s1,d1,r1
      s1 = s1*d2r
      d1 = d1*d2r
      r1 = r1*d2r
C x in strike direction, y in hor updip direction
      x1 = dcos(r1)
      y1 = dsin(r1)*dcos(d1)
      z1 = dsin(r1)*dsin(d1)
C rotate to east-north
      e1 = x1*dsin(s1) - y1*dcos(s1)
      n1 = x1*dcos(s1) + y1*dsin(s1)
      write(soln,1001) datan2(e1,n1)*r2d,
     1                 datan2(z1,dsqrt(e1*e1+n1*n1))*r2d
 1001 format(2F10.3)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE pnt2ter(line,soln)
      IMPLICIT none
      CHARACTER*180 line,soln
      REAL*8 p(3),n(3),t(3),fth,fss,fno
      read(line,*) p(1),p(2),p(3),n(1),n(2),n(3),t(1),t(2),t(3)
      fno = datan(p(3)/dsqrt(p(1)*p(1)+p(2)*p(2)))
      fss = datan(n(3)/dsqrt(n(1)*n(1)+n(2)*n(2)))
      fth = datan(t(3)/dsqrt(t(1)*t(1)+t(2)*t(2)))
      fno = dsin(fno)*dsin(fno)
      fss = dsin(fss)*dsin(fss)
      fth = dsin(fth)*dsin(fth)
      write(soln,1001) fth,fss,fno
 1001 format(3F7.3)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE pnt2mag(line,soln)
      IMPLICIT none
      CHARACTER*180 line,soln
      REAL*8 mom
      call pnt2mom(line,soln)
      read(soln,*) mom
      write(line,1001) mom
 1001 format(1P1E14.6)
      call mom2mag(line,soln)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE pnt2mom(line,soln)
      IMPLICIT none
      CHARACTER*180 line,soln
      REAL*8 dum,p,n,t,mom
      read(line,*) dum,dum,dum,dum,dum,dum,dum,dum,dum,p,n,t
      mom = 0.5d0*(dabs(p)+dabs(t))
      write(soln,1001) mom
 1001 format(1P1E14.6)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE mij2ter(line,soln)
      IMPLICIT none
      CHARACTER*180 line,soln
      REAL*8 p(3),n(3),t(3),fth,fss,fno
      call mij2pnt(line,soln)
      read(soln,*) p(1),p(2),p(3),n(1),n(2),n(3),t(1),t(2),t(3)
      fno = datan(p(3)/dsqrt(p(1)*p(1)+p(2)*p(2)))
      fss = datan(n(3)/dsqrt(n(1)*n(1)+n(2)*n(2)))
      fth = datan(t(3)/dsqrt(t(1)*t(1)+t(2)*t(2)))
      fno = dsin(fno)*dsin(fno)
      fss = dsin(fss)*dsin(fss)
      fth = dsin(fth)*dsin(fth)
      write(soln,1001) fth,fss,fno
 1001 format(3F7.3)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE mij2mom(line,soln)
C----
C Compute the seismic moment (N-m) from moment tensor elements (N-m)
C----
      IMPLICIT none
      CHARACTER*180 line,soln
      REAL*8 dum,p,n,t,mom
      call mij2pnt(line,soln)
      read(soln,*) dum,dum,dum,dum,dum,dum,dum,dum,dum,p,n,t
      mom = 0.5d0*(dabs(p)+dabs(t))
      write(soln,1001) mom
 1001 format(1P1E14.6)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE mom2mag(line,soln)
C----
C Compute the magnitude from seismic moment (N-m)
C----
      IMPLICIT none
      CHARACTER*180 line,soln
      REAL*8 mom,mag
      read(line,*) mom
      mag = 2.0d0/3.0d0*dlog10(mom*1.0d7)-10.7d0
      write(soln,1001) mag
 1001 format(1F8.2)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE mag2mom(line,soln)
      IMPLICIT none
      CHARACTER*180 line,soln
      REAL*8 mom,mag
      read(line,*) mag
      mom = 10**((mag+10.7d0)*1.5d0)*1.0d-7
      write(soln,1001) mom
 1001 format(1P1E14.6)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE mij2mag(line,soln)
      IMPLICIT none
      CHARACTER*180 line,soln
      REAL*8 mom
      call mij2mom(line,soln)
      read(soln,*) mom
      write(line,1001) mom
 1001 format(1P1E14.6)
      call mom2mag(line,soln)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE sdr2pnt(line,soln)
C----
C Compute the (unit) P, N, T axes from double couple strike, dip, and
C rake angles.
C----
      IMPLICIT none
      CHARACTER*180 line,soln
      REAL*8 mrr,mtt,mpp,mrt,mrp,mtp
C Convert to moment tensor components
      call sdr2mij(line,soln)
      read(soln,*) mrr,mtt,mpp,mrt,mrp,mtp
      write(line,9989) mrr,mtt,mpp,mrt,mrp,mtp
 9989 format(6(1P1E14.6))
C Convert to P, N, T axes
      call mij2pnt(line,soln)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE pntxyz2plaz(line,soln)
      IMPLICIT none
      REAL*8 pi,r2d
      PARAMETER (pi=4.0d0*datan(1.0d0),r2d=1.8d2/pi)
      CHARACTER*180 line,soln
      REAL*8 PNT(3,3),H
      REAL*8 pl(3),az(3),mg(3)
      INTEGER i
C Read in P, N, T vectors in Cartesian coordinates
      read(line,*) PNT(1,1),PNT(1,2),PNT(1,3),
     1             PNT(2,1),PNT(2,2),PNT(2,3),
     2             PNT(3,1),PNT(3,2),PNT(3,3)
C Compute plunge, azimuth, length
      do 103 i = 1,3
          H = dsqrt(PNT(i,1)*PNT(i,1)+PNT(i,2)*PNT(i,2))
          pl(i) = datan(PNT(i,3)/H)*r2d
          az(i) = datan2(PNT(i,1),PNT(i,2))*r2d
          if (pl(i).gt.0.0d0) then
              az(i) = az(i) + 1.8d2
              pl(i) = -pl(i)
          endif
          pl(i) = -pl(i)
          mg(i) = dsqrt(H*H+PNT(i,3)*PNT(i,3))
  103 continue
C Write to new output
      write(soln,1001) az(1),pl(1),mg(1),az(2),pl(2),mg(2),az(3),pl(3),
     1                 mg(3)
 1001 format(3(2F6.0,1P1E14.6))
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE mij2sdr(line,soln)
C----
C Compute the best fitting double couple source (strike, dip, rake) of a
C moment tensor source.
C----
      IMPLICIT none
      CHARACTER*180 line,soln
      INTEGER n,lda,lwmax
      PARAMETER (n=3,lda=n,lwmax=250)
      INTEGER info,lwork
      REAL*8 a(lda,n),w(n),work(lwmax),m(6)
      REAL*8 p(3),t(3),slip(3),normal(3),str,dip,rak,str2,dip2,rak2
      INTEGER i
C Get input moment tensor (mrr, mtt, mpp, mrt, mrp, mtp; r=z, t=s, p=e)
C Moment tensor components (r=z, t=e, p=n):
      read(line,*) m(1),m(2),m(3),m(4),m(5),m(6)
      a(3,3) =  m(1)     ! mrr =  mzz
      a(2,2) =  m(2)     ! mtt =  myy
      a(1,1) =  m(3)     ! mpp =  mxx
      a(2,3) = -m(4)     ! mrt = -myz
      a(3,2) = -m(4)
      a(1,3) =  m(5)     ! mrp =  mxz
      a(3,1) =  m(5)
      a(1,2) = -m(6)     ! mtp = -mxy
      a(1,2) = -m(6)
C Query the optimal workspace
      lwork = -1
      call dsyev('Vectors','Upper',n,a,lda,w,work,lwork,info)
      lwork = min(lwmax,int(work(1)))
C Solve eigenproblem
      call dsyev('Vectors','Upper',n,a,lda,w,work,lwork,info)
C Compute best-fitting focal mechanism
      do 103 i = 1,3
          p(i) = a(i,1)
          t(i) = a(i,3)
  103 continue
      call pt2un(p,t,slip,normal)
      call un2sdr(slip,normal,str,dip,rak)
      write(line,*) str,dip,rak
      call sdr2sdr(line,soln)
      read(soln,*) str2,dip2,rak2
      write(soln,1001) str,dip,rak,str2,dip2,rak2
 1001 format(6F6.0)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE mij2sv(line,soln)
      IMPLICIT none
      CHARACTER*180 line,soln
      REAL*8 s1,d1,r1,s2,d2,r2,az1,pl1,az2,pl2
      call mij2sdr(line,soln)
      read(soln,*) s1,d1,r1,s2,d2,r2
      write(line,*) s1,d1,r1
      call sdr2sv(line,soln)
      read(soln,*) az1,pl1
      write(line,*) s2,d2,r2
      call sdr2sv(line,soln)
      read(soln,*) az2,pl2
      write(soln,1001) az1,pl1,az2,pl2
 1001 format(2F10.3,X,2F10.3)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE pnt2sdr(line,soln)
      IMPLICIT none
      CHARACTER*180 line,soln
      REAL*8 p(3),n(3),t(3),slip(3),normal(3),s1,d1,r1,s2,d2,r2
      read(line,*) p(1),p(2),p(3),n(1),n(2),n(3),t(1),t(2),t(3)
      call pt2un(p,t,slip,normal)
      call un2sdr(slip,normal,s1,d1,r1)
      write(line,*) s1,d1,r1
      call sdr2sdr(line,soln)
      read(soln,*) s2,d2,r2
      write(soln,1001) s1,d1,r1,s2,d2,r2
 1001 format(6F6.0)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE pnt2sv(line,soln)
      IMPLICIT none
      CHARACTER*180 line,soln
      REAL*8 s1,d1,r1,s2,d2,r2,az1,pl1,az2,pl2
      call pnt2sdr(line,soln)
      read(soln,*) s1,d1,r1,s2,d2,r2
      write(line,*) s1,d1,r1
      call sdr2sv(line,soln)
      read(soln,*) az1,pl1
      write(line,*) s2,d2,r2
      call sdr2sv(line,soln)
      read(soln,*) az2,pl2
      write(soln,1001) az1,pl1,az2,pl2
 1001 format(2F10.3,X,2F10.3)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE mij2dcp(line,soln)
      IMPLICIT none
      CHARACTER*180 line,soln
      REAL*8 dum,p,n,t,dcp
      call mij2pnt(line,soln)
      read(soln,*) dum,dum,dum,dum,dum,dum,dum,dum,dum,p,n,t
      if (dabs(p).gt.dabs(t)) then
          dcp = 1.0d0 - 2.0d0*dabs(n/p)
      else
          dcp = 1.0d0 - 2.0d0*dabs(n/t)
      endif
      write(soln,1001) dcp
 1001 format(F5.2)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE mij2pnt(line,soln)
C----
C Compute P, N, and T vectors from a moment tensor input
C----
      IMPLICIT none
      CHARACTER*180 line,soln
      INTEGER n,lda,lwmax
      PARAMETER (n=3,lda=n,lwmax=250)
      INTEGER info,lwork
      REAL*8 a(lda,n),w(n),work(lwmax),m(6)
      REAL*8 p(3),o(3),t(3)
C Get input moment tensor (mrr, mtt, mpp, mrt, mrp, mtp; r=z, t=s, p=e)
C Moment tensor components (r=z, t=e, p=n):
      read(line,*) m(1),m(2),m(3),m(4),m(5),m(6)
      a(3,3) =  m(1)     ! mrr =  mzz
      a(2,2) =  m(2)     ! mtt =  myy
      a(1,1) =  m(3)     ! mpp =  mxx
      a(2,3) = -m(4)     ! mrt = -myz
      a(3,2) = -m(4)
      a(1,3) =  m(5)     ! mrp =  mxz
      a(3,1) =  m(5)
      a(1,2) = -m(6)     ! mtp = -mxy
      a(1,2) = -m(6)
C Query the optimal workspace
      lwork = -1
      call dsyev('Vectors','Upper',n,a,lda,w,work,lwork,info)
      lwork = min(lwmax,int(work(1)))
C Solve eigenproblem
      call dsyev('Vectors','Upper',n,a,lda,w,work,lwork,info)
      p(1) = a(1,1)
      p(2) = a(2,1)
      p(3) = a(3,1)
      o(1) = a(1,2)
      o(2) = a(2,2)
      o(3) = a(3,2)
      t(1) = a(1,3)
      t(2) = a(2,3)
      t(3) = a(3,3)
C Write P, N, and T vectors
      write(soln,1001) p(1),p(2),p(3),o(1),o(2),o(3),t(1),t(2),t(3),
     1                 w(1),w(2),w(3)
 1001 format(9F7.3,3(1P1E14.6))
      RETURN
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
      REAL*8 hor,st(3),dum
      INTEGER i
C      do 101 i = 1,3
C          slip(i)   = slip(i)
C  101 continue
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

C----------------------------------------------------------------------C

      SUBROUTINE pnt2dcp(line,soln)
      IMPLICIT none
      CHARACTER*180 line,soln
      REAL*8 dum,p,n,t,dcp
      read(line,*) dum,dum,dum,dum,dum,dum,dum,dum,dum,p,n,t
      if (dabs(p).gt.dabs(t)) then
          dcp = 1.0d0 - 2.0d0*dabs(n/p)
      else
          dcp = 1.0d0 - 2.0d0*dabs(n/t)
      endif
      write(soln,1001) dcp
 1001 format(F5.2)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE pnt2mij(line,soln)
      IMPLICIT none
      CHARACTER*(*) line,soln
      REAL*8 pnt(3,3),eig(3)
      REAL*8 mxx,myy,mzz,mxy,mxz,myz
      REAL*8 mrr,mtt,mpp,mrt,mrp,mtp
      read(line,*) pnt(1,1),pnt(1,2),pnt(1,3),
     1             pnt(2,1),pnt(2,2),pnt(2,3),
     2             pnt(3,1),pnt(3,2),pnt(3,3),
     3             eig(1),eig(2),eig(3)
      mxx = eig(1)*pnt(1,1)*pnt(1,1) + eig(2)*pnt(2,1)*pnt(2,1) + 
     1                                          eig(3)*pnt(3,1)*pnt(3,1)
      myy = eig(1)*pnt(1,2)*pnt(1,2) + eig(2)*pnt(2,2)*pnt(2,2) +
     1                                          eig(3)*pnt(3,2)*pnt(3,2)
      mzz = eig(1)*pnt(1,3)*pnt(1,3) + eig(2)*pnt(2,3)*pnt(2,3) +
     1                                          eig(3)*pnt(3,3)*pnt(3,3)
      mxy = eig(1)*pnt(1,1)*pnt(1,2) + eig(2)*pnt(2,1)*pnt(2,2) +
     1                                          eig(3)*pnt(3,1)*pnt(3,2)
      mxz = eig(1)*pnt(1,1)*pnt(1,3) + eig(2)*pnt(2,1)*pnt(2,3) +
     1                                          eig(3)*pnt(3,1)*pnt(3,3)
      myz = eig(1)*pnt(1,2)*pnt(1,3) + eig(2)*pnt(2,2)*pnt(2,3) +
     1                                          eig(3)*pnt(3,2)*pnt(3,3)
      mrr = mzz
      mtt = myy
      mpp = mxx
      mrt = -myz
      mrp = mxz
      mtp = -mxy
C Write moment tensor components
      write(soln,9999) mrr,mtt,mpp,mrt,mrp,mtp
 9999 format(6(1PE14.6))
      RETURN
      END
