      PROGRAM makefaultstxt
      IMPLICIT none
      INTEGER f,nflt,typ,pf
      REAL lon,lat,dep,str,dip,rak,wid,len,slip,mag,mom,area,mod

      open (unit=11,file='faults.txt',status='unknown')

      print *,'How many faults?'
      read *,nflt

      f = 0
  101 f = f + 1
      if (f.gt.nflt) goto 102

      write (6,9999) f
      read *,lon,lat,dep

      write (6,9998) f
      read *,str,dip,rak

      write (6,9997) f
      read *,typ

      write (6,9996) f
      read *,mag
      call wellscoppersmith(wid,len,mag,typ)
      mom = (10**(1.5*(mag+10.7)))/1e7
      area = wid*len*1e6
      mod = 35.0e9
      slip = mom/(area*mod)
      print *,mom,area,mod

      write (6,9995) f
      read *,pf

      write (11,9994) lon,lat,dep,str,dip,rak,slip,wid,len,pf
      goto 101

  102 continue


 9999 format('Enter fault #',I2,' lon, lat, depth (km):')
 9998 format('Enter fault #',I2,' strike, dip, rake:')
 9997 format('Enter fault #',I2,' type (SS=1,RV=2,NO=3)')
 9996 format('Enter fault #',I2,' magnitude')
 9995 format('Treat fault #',I2,' as point (0) or finite (1)?')
 9994 format(F8.3,X,F8.3,X,F5.1,X,F5.1,X,F4.1,X,F6.1,X,F8.4,X,F5.1,X
     1       F6.1,X,I1)
      END

      SUBROUTINE wellscoppersmith(wid,len,mag,typ)
      IMPLICIT none
      INTEGER typ
      REAL wid,len,mag

      if (typ.eq.1) then
          len = 10.0**(-2.57+0.62*mag)
          wid = 10.0**(-0.76+0.27*mag)
      elseif (typ.eq.2) then
          len = 10.0**(-2.42+0.58*mag)
          wid = 10.0**(-1.61+0.41*mag)
      elseif (typ.eq.3) then
          len = 10.0**(-1.88+0.50*mag)
          wid = 10.0**(-1.14+0.35*mag)
      else
          len = 10.0**(-2.44+0.59*mag)
          wid = 10.0**(-1.01+0.32*mag)
      endif

      RETURN
      END
