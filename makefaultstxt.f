      PROGRAM makefaultstxt
      IMPLICIT none
      INTEGER f,nflt,typ,pf,batch
      REAL*8 lon,lat,dep,str,dip,rak,wid,len,slip,mag,mom,area,mod

      call gcmdln(batch)
      open (unit=11,file='faults.txt',status='unknown')

      if (batch.eq.0) then
          print *,'How many faults?'
          read *,nflt

          f = 0
  101     f = f + 1
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
          mom = (10.0d0**(1.5d0*(mag+10.7d0)))/1.0d7
          area = wid*len*1.0d6
          mod = 35.0d9
          slip = mom/(area*mod)
          print *,mom,area,mod
    
          write (6,9995) f
          read *,pf
    
          write (11,9994) lon,lat,dep,str,dip,rak,slip,wid,len,pf
          goto 101

  102     continue

      else
          !print *,'File must be lon,lat,dep,str,dip,rak,typ,mag'
          pf = 1
          open (unit=12,file='makefaults.in',status='old')
  103     read(12,*,end=104) lon,lat,dep,str,dip,rak,typ,mag
              call wellscoppersmith(wid,len,mag,typ)
              mom = (10.0d0**(1.5d0*(mag+10.7d0)))/1.0d7
              area = wid*len*1.0d6
              mod = 35.0d9
              slip = mom/(area*mod)
              write (11,9994) lon,lat,dep,str,dip,rak,slip,wid,len,pf
              goto 103
  104     continue
      endif

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
      REAL*8 wid,len,mag

      if (typ.eq.1) then
          len = 10.0d0**(-2.57d0+0.62d0*mag)
          wid = 10.0d0**(-0.76d0+0.27d0*mag)
      elseif (typ.eq.2) then
          len = 10.0d0**(-2.42d0+0.58d0*mag)
          wid = 10.0d0**(-1.61d0+0.41d0*mag)
      elseif (typ.eq.3) then
          len = 10.0d0**(-1.88d0+0.50d0*mag)
          wid = 10.0d0**(-1.14d0+0.35d0*mag)
      else
          len = 10.0d0**(-2.44d0+0.59d0*mag)
          wid = 10.0d0**(-1.01d0+0.32d0*mag)
      endif

      RETURN
      END

      SUBROUTINE gcmdln(batch)
      
      IMPLICIT NONE
      CHARACTER*30 tag
      INTEGER i,narg,batch
      
      batch = 0

      narg = iargc()
      i = 0
 9998 i = i + 1
      if (i.gt.narg) goto 9999
          call getarg(i,tag)
          if (tag(1:6).eq.'-batch') then
              batch = 1
          endif
          goto 9998
 9999 continue

      RETURN
      END
