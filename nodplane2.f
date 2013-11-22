      PROGRAM nodplane2

      IMPLICIT none
      REAL pi,d2r,r2d
      PARAMETER (pi=3.14159265,d2r=pi/180.0,r2d=180.0/pi)
      REAL str1,dip1,rak1,str2,dip2,rak2
      REAL x,y,z


      print *,"Strike?"
      read *,str1

      print *,"Dip?"
      read *,dip1

      print *,"Rake?"
      read *,rak1

      x = cos(rak1*d2r)
      y = sin(rak1*d2r)*cos(dip1*d2r)
      z = sin(rak1*d2r)*sin(dip1*d2r)

      if (z.gt.0.and.x.gt.0) then
          str2 = str1 - 90.0 - r2d*atan(y/x)
      elseif (z.gt.0.and.x.lt.0) then
          str2 = str1 + 90.0 + r2d*atan(abs(y/x))
      elseif (z.lt.0.and.x.gt.0) then
          z = -z
          x = -x
          y = -y
          str2 = str1 - 180.0 - r2d*atan(x/y)
      elseif (z.lt.0.and.x.lt.0) then
          z = -z
          x = -x
          y = -y
          str2 = str1 + 90.0 - r2d*atan(y/x)
      elseif (x.eq.1.0) then
          str2 = str1 - 90.0
      endif

      dip2 = 90.0 - r2d*atan(z/(sqrt(x*x+y*y)))
      rak2 = r2d*acos(sin(d2r*(str2-str1))*cos(d2r*(90.0-dip1)))

      print *,"New strike,dip,rake: ",str2,",",dip2,",",rak2


      END
