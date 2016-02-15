      PROGRAM nodplane2
C----
C Compute the second nodal plane of a double couple source
C given strike, dip, rake
C
C x: parallel to input strike
C y: horizontal up-dip relative to input strike
C z: vertical up
C----
      IMPLICIT none
      REAL*8 pi,d2r,r2d
      PARAMETER (pi=4.0d0*datan(1.0d0),d2r=pi/1.8d2,r2d=1.8d2/pi)
      REAL*8 str1,dip1,rak1,str,dip,rak,str2,dip2,rak2
      REAL*8 xn1,yn1,zn1,xr,yr,zr,xs2,ys2,zs2,dum1,dum2
      CHARACTER*80 txt

C Read from standard input
  100 read (*,'(A)',end=101) txt
      read(txt,*) str1,dip1,rak1
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
      !print *,"NOR ",xn1,yn1,zn1
      !print *,"RAK ",xr,yr,zr

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
      !print *,"NEW RAK ",xs2,ys2,zs2
      !print *,xn1*xs2+yn1*ys2+zn1*zs2
      rak2 = r2d*dacos(xn1*xs2+yn1*ys2+zn1*zs2)

C Check rake is pos or neg
      dum1 = xs2*yn1 - ys2*xn1
      if (dum1.lt.0.0d0) then
          rak2 = -rak2
      endif

      write(*,9999) str2,dip2,rak2
      goto 100
  101 continue

9999  format(3F12.1)
      END


C      if (z.gt.0.and.x.gt.0) then
C          str2 = str1 - 90.0 - r2d*atan(y/x)
C      elseif (z.gt.0.and.x.lt.0) then
C          str2 = str1 + 90.0 + r2d*atan(abs(y/x))
C      elseif (z.lt.0.and.x.gt.0) then
C          z = -z
C          x = -x
C          y = -y
C          str2 = str1 - 180.0 - r2d*atan(x/y)
C      elseif (z.lt.0.and.x.lt.0) then
C          z = -z
C          x = -x
C          y = -y
C          str2 = str1 + 90.0 - r2d*atan(y/x)
C      elseif (x.eq.1.0) then
C          str2 = str1 - 90.0
C      endif

