      PROGRAM simpleinversion
C----
C Linear inversion Ax = y
C   x: model parameters (m x 1)
C   y: observations (n x 1)
C   A: model matrix (n x m)
C
C Assumptions about geometry:
C   - Positive x is in strike direction
C   - Positive y is in updip horizontal direction
C   - Fault segment depths are positive down
C   - Vertical displacements are positive up
C   - Trench is at y = 0, z = 0
C   - Plate boundary has constant dip
C----
      IMPLICIT none
      REAL*8 pi,d2r,r2d
      PARAMETER (pi=4.d0*atan(1.d0),d2r=pi/1.8d2,r2d=1.8d2/pi)
      INTEGER STOUT,STERR
      PARAMETER (STOUT=6,STERR=0)

      INTEGER NMAX,MMAX ! Maximum dimensions (NMAX x MMAX) of A matrix
      PARAMETER (NMAX=1000,MMAX=1000)
      REAL*8 A(NMAX,MMAX),damp,smooth

      INTEGER i,j
      INTEGER nobs,nhor,mdip,N,M

      REAL*8 xobs,yobs(NMAX),zobs,disp(NMAX)

      REAL*8 dip,rak,totaldip,ddip,totalstr
      REAL*8 xflt,yflt(MMAX),zflt(MMAX),x,y,z

      REAL*8 vp,vs,dens

      REAL*8 gfslip,dum1,dum2

      REAL*8 slip(MMAX),chisq,TOL,rms,ux,uy,uz,uynet,uznet

C----
C Load observation vector, disp ("y" from Ax = y)
C   nobs: total number of observations
C   nhor: number of horizontal observations
C
C     y_position displacement
C       (km)         (m)
C----
      open (unit=11,file='horz.dat',status='old')
      nhor = 0
      i = 0
  201 i = i + 1
      read (11,*,end=202) yobs(i),disp(i)
          yobs(i) = 1.0d3*yobs(i)
          nhor = nhor + 1
          goto 201
  202 continue
      nobs = nhor

      open (unit=12,file='vert.dat',status='old')
      i = 0
  203 i = i + 1
      read (12,*,end=204) yobs(nhor+i),disp(nhor+i)
          yobs(nhor+i) = 1.0d3*yobs(nhor+i)
          nobs = nobs + 1
          goto 203
  204 continue

      zobs = 0.0d0
      xobs = 0.0d0

C----
C Create model matrix ("A" from Ax = y)
C     totalwid: fault extent along dip
C     mdip:     # of fault segments along dip
C----
      open (unit=13,file="input.dat",status="old")
      read (13,*) dip,totaldip,mdip,totalstr
      read (13,*) damp
      read (13,*) smooth

      dip = dip*d2r
      rak = 90.d0
      totaldip = totaldip*1.0d3
      ddip = totaldip/dble(mdip)
      totalstr = totalstr*1.d3

      do 101 j = 1,mdip
          yflt(j) = -(ddip*dcos(dip)*(dble(j)-0.5d0))
          zflt(j) =   ddip*dsin(dip)*(dble(j)-0.5d0)
  101 continue

      dip = dip*r2d

C Halfspace properties
      vp = 6800.0d0
      vs = vp/dsqrt(3.0d0)
      dens = 3000.0d0

C Unit fault slip for Green Functions
      gfslip = 1.0d0

C Load A
      x = xobs
      do 304 i = 1,nobs
          if (i.le.nhor) then
              do 302 j = 1,mdip
                  y = yobs(i) - yflt(j)
                  call o92rect(dum1,A(i,j),dum2,x,y,zobs,zflt(j),
     1                          dip,rak,ddip,totalstr,gfslip,vp,vs,dens)
  302         continue
          else
              do 303 j = 1,mdip
                  y = yobs(i) - yflt(j)
                  call o92rect(dum1,dum2,A(i,j),x,y,zobs,zflt(j),
     1                          dip,rak,ddip,totalstr,gfslip,vp,vs,dens)
  303         continue
          endif
  304 continue

C----
C Add damping and/or smoothing
C----
      do 305 i = 1,mdip
          do 306 j = 1,mdip
              if (i.eq.j) then
                  A(nobs+i,j) = damp*damp*1.0d0
              else
                  A(nobs+i,j) = 0.0d0
              endif
  306     continue
  305 continue

      do 307 i = 1,mdip
          do 308 j = 1,mdip
              if (i.eq.1.and.j.eq.1) then
                  A(2*nobs+i,j) = smooth*smooth*1.0d0
              elseif (i.eq.mdip.and.j.eq.mdip) then
                  A(2*nobs+i,j) = smooth*smooth*1.0d0
              elseif (i.eq.j) then
                  A(2*nobs+i,j) = smooth*smooth*2.0d0
              elseif (i.eq.j-1) then
                  A(2*nobs+i,j) = -smooth*smooth*1.0d0
              elseif (i.eq.j+1) then
                  A(2*nobs+i,j) = -smooth*smooth*1.0d0
              else
                  A(2*nobs+i,j) = 0.0d0
              endif
  308     continue
  307 continue

C      do 309 i = 1,mdip
C          print *,(A(2*nobs+i,j),j=1,mdip)
C  309 continue

      N = nobs + mdip
      M = mdip

C----
C Perform inversion for slip, write into faults.txt format
C----
      call dsvdfit(A,disp,N,M,NMAX,MMAX,slip,chisq,TOL)
      open (unit=15,file='faults.txt',status='unknown')
      do 401 j = 1,mdip
          write (15,8999),x,yflt(j),zflt(j),0.0,dip,90.0,slip(j),
     1                    ddip,totalstr,1
  401 continue

 8999 format(3(F11.3,X),F3.1,X,2(F6.1,X),F14.3,X,2(F9.1,X),I1)

C----
C Calculate rms misfit
C----
      rms = 0.0d0
      do 501 i = 1,nobs
          if (i.le.nhor) then
              uynet = 0.0d0
              do 502 j = 1,mdip
                  y = yobs(i) - yflt(j)
                  call o92rect(dum1,uy,dum2,x,y,zobs,zflt(j),
     1                         dip,rak,ddip,totalstr,slip(j),vp,vs,dens)
                  uynet = uynet + uy
  502         continue
              rms = rms + (uynet-disp(i))*(uynet-disp(i))
          else
              uznet = 0.0d0
              do 503 j = 1,mdip
                  y = yobs(i) - yflt(j)
                  call o92rect(dum1,dum2,uz,x,y,zobs,zflt(j),
     1                         dip,rak,ddip,totalstr,slip(j),vp,vs,dens)
                  uznet = uznet + uz
  503         continue
              rms = rms + (uznet-disp(i))*(uznet-disp(i))
          endif
  501 continue
      rms = dsqrt(rms)

C----
C Write results to standard error and standard out
C----
      write (STOUT,9994) damp,smooth,rms
      write (STERR,9995) nobs
      write (STERR,9996) mdip
      write (STERR,9997) damp
      write (STERR,9998) smooth
      write (STERR,9999) rms

 9994 format(3F14.3)
 9995 format('NOBS          -',I14)
 9996 format('MDIP          -',I14)
 9997 format('DAMPING       -',F14.3)
 9998 format('SMOOTHING     -',F14.3)
 9999 format('RMS MISFIT    -',F14.3)

      END
