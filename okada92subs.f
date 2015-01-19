C======================================================================C
C SUBROUTINES TO COMPUTE DEFORMATION FROM POINT AND FINITE SOURCES IN
C AN ELASTIC HALF-SPACE.
C
C REFERENCE:
C OKADA, Y., 1992, INTERNAL DEFORMATION DUE TO SHEAR AND TENSILE
C   FAULTS IN A HALF-SPACE. BULLETIN OF THE SEISMOLOGICAL SOCIETY OF
C   AMERICA 82, PP. 1018-1040.
C
C NOTES:
C - ALL SUBROUTINE INPUTS ARE IN SI UNITS.
C - ANGLE INPUTS ARE IN DEGREES.
C - COORDINATE SYSTEM IS RELATIVE TO THE FAULT ORIENTATION:
C   - X POINTS IN THE ALONG-STRIKE DIRECTION
C   - Y POINTS IN THE HORIZONTAL UP-DIP DIRECTION
C   - DEPTHS ARE DEFINED POSITIVE DOWN
C======================================================================C

      SUBROUTINE o92pt(ux,uy,uz,x,y,stdp,evdp,dipin,rakin,area,slip,
     1                 vp,vs,dens)
C----
C COMPUTE COMPONENTS OF THE DISPLACEMENT VECTOR DUE TO A POINT SOURCE.
C----
      IMPLICIT none
      REAL*8 x,y,stdp,evdp,z,c,d
      REAL*8 dipin,rakin,area,slip,Mss,Mds
      REAL*8 vp,vs,dens
      INTEGER i,j
      REAL*8 u(6,3),f(6,3),ux,uy,uz
      REAL*8 zro
      DATA zro/0.0d0/

      REAL*8 sd,cd,s2d,c2d,cdcd,sdsd,cdsd
      REAL*8 a,CA1,CA2,CB,CC
      REAL*8 R,R2,R3,R4,R5,R7,Rd,p,q,s,t
      REAL*8 A3,A5,A7,B3,B5,B7,C3,C5,C7
      REAL*8 I1,I2,I3,I4,I5
      REAL*8 J1,J2,J3,J4
      REAL*8 K1,K2,K3
      REAL*8 U2,V2,W2
      REAL*8 U3,V3,W3
      REAL*8 xx,xy,yy,dd,xd,xc,xq,yq,dq,pq
      INTEGER thru

      COMMON /SOURCE/ sd,cd,s2d,c2d,cdcd,sdsd,cdsd
      COMMON /ELAST/  a,CA1,CA2,CB,CC
      COMMON /RVARS0/ R,R2,R3,R4,R5,R7,Rd,p,q,s,t
      COMMON /ABC0/   A3,A5,A7,B3,B5,B7,C3,C5,C7
      COMMON /IVARS0/ I1,I2,I3,I4,I5
      COMMON /JVARS0/ J1,J2,J3,J4
      COMMON /KVARS0/ K1,K2,K3
      COMMON /YVARS0/ U2,V2,W2
      COMMON /ZVARS0/ U3,V3,W3
      COMMON /PRODUCTS/ xx,xy,yy,dd,xd,xc,xq,yq,dq,pq
      COMMON /TAG/    thru

C----
C INITIALIZE COMPONENTS OF DISPLACEMENT
C   U(1-3): DISPLACEMENT FROM UNIT STRIKE-SLIP SOURCE
C   U(4-6): DISPLACEMENT FROM UNIT DIP-SLIP SOURCE
C   F(I,J): DUMMY DISPLACEMENT VARIABLES
C thru   - a tag indicating number of times through calculation;
C          used for generating mirror source
C----
      thru = 0
      do 102 i = 1,6
          do 101 j = 1,3
              u(i,j) = zro
  101     continue
  102 continue

C----
C Source and halfspace constants
C----
      call DIPVAR(dipin)
      call hafspc(vp,vs,dens)
      call moment0(Mss,Mds,slip,rakin,area)

C----
C Calculate displacement components (done in sub disp0)
C x in direction of strike
C y in updip direction (horizontal)
C z vertical
C----
      c = evdp
      z = -stdp
  103 d = c - z

      call geomvars0(x,y,c,d)

      if (thru.eq.0) then
          call disp0(f,x,y,c,d)
          do 105 i = 1,6
              do 104 j = 1,3
                  u(i,j) = u(i,j) + f(i,j)
  104         continue
  105     continue
          thru = 1
          z = stdp
          goto 103
      else
          call disp0(f,x,y,c,d)
          do 106 i = 1,6
              u(i,1) = u(i,1) - f(i,1)
  106     continue
          z = -stdp
      endif

C----
C Double double toil and trouble...put it all together...and....PRESTO!
C Combine displacement components to get displacements at the receiver
C (+x in strike direction; +y in updip horizontal direction)
C----
      ux = Mss*(u(1,1)+u(1,2)+z*u(1,3)) + Mds*(u(4,1)+u(4,2)+z*u(4,3))
      uy = Mss*(u(2,1)+u(2,2)+z*u(2,3)) + Mds*(u(5,1)+u(5,2)+z*u(5,3))
      uz = Mss*(u(3,1)+u(3,2)+z*u(3,3)) + Mds*(u(6,1)+u(6,2)+z*u(6,3))
      
      RETURN 
      END

C----------------------------------------------------------------------C

      SUBROUTINE o92ptstn(strain,x,y,stdp,evdp,dipin,rakin,area,slip,
     1                    vp,vs,dens)
C      SUBROUTINE STN0(STRAIN,X,Y,STDP,EVDP,DIPIN,RAKIN,AREA,SLIP,VP,VS,DENS)
C----
C Static elastic strain at a location (internal or surface) due to
C a point source shear dislocation in an isotropic halfspace.
C
C INPUTS (ALL IN SI UNITS):
C   - Geometric Parameters (relative to source location):
C         x (along-strike), y (horizontal up-dip), station depth,
C         event depth
C   - Source parameters:
C         fault dip (deg), fault rake (deg), fault area, slip
C   - Physical parameters of half-space:
C         vp, vs, density
C OUTPUT: strain matrix at receiver (in x-y coordinates)
C----
      IMPLICIT NONE
      REAL*8 x,y,stdp,evdp,z,c,d
      REAL*8 dipin,rakin,area,slip,Mss,Mds
      REAL*8 vp,vs,dens
      INTEGER i,j
      REAL*8 fx(6,3),fy(6,3),fz(6,3),f(6),ux(6,3),uy(6,3),uz(6,3),u(6)
      REAL*8 uxx,uxy,uxz,uyx,uyy,uyz,uzx,uzy,uzz,strain(3,3)

      REAL*8 sd,cd,s2d,c2d,cdcd,sdsd,cdsd
      REAL*8 a,CA1,CA2,CB,CC
      REAL*8 R,R2,R3,R4,R5,R7,Rd,p,q,s,t
      REAL*8 A3,A5,A7,B3,B5,B7,C3,C5,C7
      REAL*8 I1,I2,I3,I4,I5
      REAL*8 J1,J2,J3,J4
      REAL*8 K1,K2,K3
      REAL*8 U2,V2,W2
      REAL*8 U3,V3,W3
      REAL*8 xx,xy,yy,dd,xd,xc,xq,yq,dq,pq
      INTEGER thru

      COMMON /SOURCE/ sd,cd,s2d,c2d,cdcd,sdsd,cdsd
      COMMON /ELAST/ a,CA1,CA2,CB,CC
      COMMON /RVARS0/ R,R2,R3,R4,R5,R7,Rd,p,q,s,t
      COMMON /ABC0/ A3,A5,A7,B3,B5,B7,C3,C5,C7
      COMMON /IVARS0/ I1,I2,I3,I4,I5
      COMMON /JVARS0/ J1,J2,J3,J4
      COMMON /KVARS0/ K1,K2,K3
      COMMON /YVARS0/ U2,V2,W2
      COMMON /ZVARS0/ U3,V3,W3
      COMMON /PRODUCTS/ xx,xy,yy,dd,xd,xc,xq,yq,dq,pq
      COMMON /TAG/ thru
C----
      thru = 0
      do 306 i = 1,6
          u(i) = 0.0d0
          do 305 j = 1,3
              ux(i,j) = 0.0d0
              uy(i,j) = 0.0d0
              uz(i,j) = 0.0d0
  305     continue
  306 continue
C----
C Source and halfspace constants
C----
      call DIPVAR(dipin)
      call hafspc(vp,vs,dens)
      call moment0(Mss,Mds,slip,rakin,area)
C----
C x in direction of strike
C y in updip direction (horizontal)
C z vertical
C----
      c = evdp
      z = -stdp
  301 d = c - z

      call geomvars0(x,y,c,d)

      if (thru.eq.0) then
          call xderiv0(fx,x,y,c,d)
          call yderiv0(fy,x,y,c,d)
          call zderiv0(fz,x,y,c,d)
          call disp0stn(f,x,y,c,d)

           do 308 i = 1,6
               u(i) = u(i) + f(i)
               do 307 j = 1,3
                   ux(i,j) = ux(i,j) + fx(i,j)
                   uy(i,j) = uy(i,j) + fy(i,j)
                   uz(i,j) = uz(i,j) + fz(i,j)
  307          continue
  308      continue

          thru = 1
          z = stdp
          goto 301
      else
          call xderiv0(fx,x,y,c,d)
          call yderiv0(fy,x,y,c,d)
          call zderiv0(fz,x,y,c,d)

          do 309 i = 1,6
              ux(i,1) = ux(i,1) - fx(i,1)
              uy(i,1) = uy(i,1) - fy(i,1)
              uz(i,1) = uz(i,1) + fz(i,1)
  309     continue

          z = -stdp
      endif
C----
C Chug, chug, chug...partial derivatives of displacement (+x in strike
C direction, +y in updip direction)
C----
      uxx = Mss*(ux(1,1)+ux(1,2)+z*ux(1,3))
     1                                 + Mds*(ux(4,1)+ux(4,2)+z*ux(4,3))
      uyx = Mss*(ux(2,1)+ux(2,2)+z*ux(2,3))
     1                                 + Mds*(ux(5,1)+ux(5,2)+z*ux(5,3))
      uzx = Mss*(ux(3,1)+ux(3,2)+z*ux(3,3))
     1                                 + Mds*(ux(6,1)+ux(6,2)+z*ux(6,3))
      uxy = Mss*(uy(1,1)+uy(1,2)+z*uy(1,3))
     1                                 + Mds*(uy(4,1)+uy(4,2)+z*uy(4,3))
      uyy = Mss*(uy(2,1)+uy(2,2)+z*uy(2,3))
     1                                 + Mds*(uy(5,1)+uy(5,2)+z*uy(5,3))
      uzy = Mss*(uy(3,1)+uy(3,2)+z*uy(3,3))
     1                                 + Mds*(uy(6,1)+uy(6,2)+z*uy(6,3))
      uxz = Mss*(uz(1,1)+uz(1,2)+u(1)+z*uz(1,3))
     1                         + Mds*(uz(4,1)+uz(4,2)+u(4)+z*uz(4,3))
      uyz = Mss*(uz(2,1)+uz(2,2)+u(2)+z*uz(2,3))
     1                         + Mds*(uz(5,1)+uz(5,2)+u(5)+z*uz(5,3))
      uzz = Mss*(uz(3,1)+uz(3,2)+u(3)+z*uz(3,3))
     1                         + Mds*(uz(6,1)+uz(6,2)+u(6)+z*uz(6,3))
C----
C Strain tensor from partial derivatives (in coordinate system
C with +x in direction of strike, +y updip)
C----
      strain(1,1) = uxx
      strain(2,2) = uyy
      strain(3,3) = uzz
      strain(1,2) = 5.0d-1*(uxy+uyx)
      strain(1,3) = 5.0d-1*(uxz+uzx)
      strain(2,3) = 5.0d-1*(uyz+uzy)
      strain(2,1) = strain(1,2)
      strain(3,1) = strain(1,3)
      strain(3,2) = strain(2,3)

      RETURN 
      END

C----------------------------------------------------------------------C

      SUBROUTINE o92rect(ux,uy,uz,x,y,stdp,evdp,dipin,rakin,wid,len,
     1                   slip,vp,vs,dens)
C      SUBROUTINE DSPFN(UX,UY,UZ,X,Y,STDP,EVDP,DIPIN,RAKIN,WID,LEN,SLIP,VP,VS,DENS)
C----
C Static elastic displacement at a location (internal or surface) due
C to a finite rectangular source shear dislocation in a uniform, 
C isotropic halfspace.
C
C INPUTS (ALL IN SI UNITS):
C   - Geometric Parameters (relative to source location):
C         x (along-strike), y (horizontal up-dip), station depth,
C         event depth
C   - Source parameters:
C         fault dip (deg), fault rake (deg), along-dip width, along-
C         strike length, slip
C   - Physical parameters of half-space:
C         vp, vs, density
C OUTPUT: x, y, z displacements at station (in x-y coordinates, not N-E)
C----
      IMPLICIT none
      REAL*8 eps
      PARAMETER (eps=1.0d-4)
      REAL*8 x,y,stdp,evdp,z,c,d,p,q,ksi(2),eta(2)
      REAL*8 dipin,rakin,wid,len,slip,Mss,Mds
      REAL*8 vp,vs,dens
      INTEGER i,j,ii,jj,ee(2),ek(2),eq
      REAL*8 u(6,3),f(6,3),ux,uy,uz
      REAL*8 fac,zro,one
      DATA zro,one/0.0d0,1.0d0/

      REAL*8 sd,cd,s2d,c2d,cdcd,sdsd,cdsd
      REAL*8 a,CA1,CA2,CB,CC
      REAL*8 cbar,dbar,ybar
      REAL*8 R,R2,R3,R5,X11,X32,X53,Y11,Y32,Y53,Y0,Z32,Z53,Z0
      REAL*8 I1,I2,I3,I4
      REAL*8 J1,J2,J3,J4,J5,J6
      REAL*8 K1,K2,K3,K4
      REAL*8 E2,F2,G2,H2,P2,Q2
      REAL*8 E3,F3,G3,H3,P3,Q3
      REAL*8 D11,Rd,Re,Rk,logRe,logRk,TH
      INTEGER thru

      COMMON /SOURCE/ sd,cd,s2d,c2d,cdcd,sdsd,cdsd
      COMMON /ELAST/ a,CA1,CA2,CB,CC
      COMMON /BAR/ cbar,dbar,ybar
      COMMON /RVARS/ R,R2,R3,R5,X11,X32,X53,Y11,Y32,Y53,Y0,Z32,Z53,Z0
      COMMON /IVARS/ I1,I2,I3,I4
      COMMON /JVARS/ J1,J2,J3,J4,J5,J6
      COMMON /KVARS/ K1,K2,K3,K4
      COMMON /YVARS/ E2,F2,G2,H2,P2,Q2
      COMMON /ZVARS/ E3,F3,G3,H3,P3,Q3
      COMMON /MISC/ D11,Rd,Re,Rk,logRe,logRk,TH
      COMMON /TAG/ thru
C----
C Initialize displacement components
C----
      thru = 0
      do 302 i = 1,6
          do 301 j = 1,3
              u(i,j) = zro
 301      continue
 302  continue
C----
C Source and halfspace constants
C----
      call DIPVAR(dipin)
      call hafspc(vp,vs,dens)
      call moment1(Mss,Mds,slip,rakin)
C----
C Coordinates on fault plane (x,p,q) and distance from edges (ksi,eta)
C----
      ksi(1) = x+len*0.5d0
      ksi(2) = x-len*0.5d0
C      ksi(1) = x
C      ksi(2) = x-len
      if (dabs(ksi(1)).lt.eps) ksi(1) = zro
      if (dabs(ksi(2)).lt.eps) ksi(2) = zro

      c =  evdp
      z = -stdp
  303 d =  c - z
      p = y*cd + d*sd
      q = y*sd - d*cd

      eta(1) = p+wid*0.5d0
      eta(2) = p-wid*0.5d0
C      eta(1) = p
C      eta(2) = p-wid
      if (dabs(eta(1)).lt.eps) eta(1) = zro
      if (dabs(eta(2)).lt.eps) eta(2) = zro
      if (dabs(q)     .lt.eps) q      = zro
C----
C Check for singular cases (after Okada DC3D)
C----
      call rectsingular(eq,ek,ee,ksi,eta,q,eps)
      if (eq.eq.1) then
          call singulardisplacement(ux,uy,uz)
          goto 399
      endif
C----
C Integrate solutions over fault boundaries (-W/2:W/2 and -L/2:L/2) 
C----
      do 308 ii = 1,2
      do 307 jj = 1,2
          if (ii.eq.1.and.jj.eq.1) fac =  one
          if (ii.eq.1.and.jj.eq.2) fac = -one
          if (ii.eq.2.and.jj.eq.1) fac = -one
          if (ii.eq.2.and.jj.eq.2) fac =  one

          call rectvars(ksi(ii),eta(jj),q,z,ek(jj),ee(ii),eps)

          if (thru.eq.0) then
              call disp1(f,ksi(ii),eta(jj),q,z)
              do 305 i = 1,6
                  do 304 j = 1,3
                      u(i,j) = u(i,j) + fac*f(i,j)
  304             continue
  305         continue
          else
              call disp1(f,ksi(ii),eta(jj),q,z)
              do 306 i = 1,6
                  u(i,1) = u(i,1) - fac*f(i,1)
  306         continue
          endif

  307 continue
  308 continue
C----
C Mirror source
C----
      if (thru.eq.0) then
          z = stdp 
          thru = 1
          goto 303
      else
          z = -stdp
      endif
C----
C Cross your fingers, knock on wood, or pray to your God. Whatever
C works for you.....ALAKAZAM! Displacements at the receiver. Same
C as above, +x in strike direction, +y in updip horizontal direction
C----
      ux = Mss*((u(1,1)+u(1,2)+z*u(1,3)))
     1   + Mds*((u(4,1)+u(4,2)+z*u(4,3)))
      uy = Mss*((u(2,1)+u(2,2)+z*u(2,3))*cd
     1                                  - (u(3,1)+u(3,2)+z*u(3,3))*sd)
     2   + Mds*((u(5,1)+u(5,2)+z*u(5,3))*cd
     3                                  - (u(6,1)+u(6,2)+z*u(6,3))*sd)
      uz = Mss*((u(2,1)+u(2,2)-z*u(2,3))*sd
     1                                  + (u(3,1)+u(3,2)-z*u(3,3))*cd)
     2   + Mds*((u(5,1)+u(5,2)-z*u(5,3))*sd
     3                                  + (u(6,1)+u(6,2)-z*u(6,3))*cd)
  399 continue

      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE o92rectstn(strain,x,y,stdp,evdp,dipin,rakin,wid,len,
     1                      slip,vp,vs,dens)
C----
C Static elastic strain at a location (internal or surface) due to
C a finite rectangular source shear dislocation in a uniform, isotropic
C halfspace.
C
C INPUTS (ALL IN SI UNITS):
C   - Geometric Parameters (relative to source location):
C         x (along-strike), y (horizontal up-dip), station depth,
C         event depth
C   - Source parameters:
C         fault dip (deg), fault rake (deg), along-dip width, along-
C         strike length, slip
C   - Physical parameters of half-space:
C         vp, vs, density
C OUTPUT: strain at receiver
C----
      IMPLICIT NONE
      REAL*8 eps
      PARAMETER (eps=1.0d-4)
      REAL*8 x,y,stdp,evdp,z,c,d,p,q,ksi(2),eta(2)
      REAL*8 dipin,rakin,wid,len,slip,Mss,Mds
      REAL*8 vp,vs,dens
      INTEGER i,j,ii,jj,ee(2),ek(2),eq
      REAL*8 fj(6,3),fk(6,3),fl(6,3),u(6),fx(6,3),fy(6,3),fz(6,3),f(6)
      REAL*8 uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz,strain(3,3)
      REAL*8 fac,zro,one
      DATA zro,one/0.0d0,1.0d0/

      REAL*8 sd,cd,s2d,c2d,cdcd,sdsd,cdsd
      REAL*8 a,CA1,CA2,CB,CC
      REAL*8 cbar,dbar,ybar
      REAL*8 R,R2,R3,R5,X11,X32,X53,Y11,Y32,Y53,Y0,Z32,Z53,Z0
      REAL*8 I1,I2,I3,I4
      REAL*8 J1,J2,J3,J4,J5,J6
      REAL*8 K1,K2,K3,K4
      REAL*8 E2,F2,G2,H2,P2,Q2
      REAL*8 E3,F3,G3,H3,P3,Q3
      REAL*8 D11,Rd,Re,Rk,logRe,logRk,TH
      INTEGER thru

      COMMON /SOURCE/ sd,cd,s2d,c2d,cdcd,sdsd,cdsd
      COMMON /ELAST/ a,CA1,CA2,CB,CC
      COMMON /BAR/ cbar,dbar,ybar
      COMMON /RVARS/ R,R2,R3,R5,X11,X32,X53,Y11,Y32,Y53,Y0,Z32,Z53,Z0
      COMMON /IVARS/ I1,I2,I3,I4
      COMMON /JVARS/ J1,J2,J3,J4,J5,J6
      COMMON /KVARS/ K1,K2,K3,K4
      COMMON /YVARS/ E2,F2,G2,H2,P2,Q2
      COMMON /ZVARS/ E3,F3,G3,H3,P3,Q3
      COMMON /MISC/ D11,Rd,Re,Rk,logRe,logRk,TH
      COMMON /TAG/ thru
C----
C Initialize partial derivative components
C----
      thru = 0
      do 402 i = 1,6
          u(i) = zro
          do 401 j = 1,3
              fj(i,j) = zro
              fk(i,j) = zro
              fl(i,j) = zro
  401     continue
  402 continue
C----
C Source and halfspace constants
C----
      call DIPVAR(dipin)
      call hafspc(vp,vs,dens)
      call moment1(Mss,Mds,slip,rakin) 
C----
C Coordinates on fault plane (x,p,q) and distance from edges (ksi,eta)
C----
      ksi(1) = x+len*0.5d0
      ksi(2) = x-len*0.5d0
C      ksi(1) = x
C      ksi(2) = x-len
      if (dabs(ksi(1)).lt.eps) ksi(1) = zro
      if (dabs(ksi(2)).lt.eps) ksi(2) = zro

      c =  evdp
      z = -stdp
  403 d =  c - z
      p = y*cd + d*sd
      q = y*sd - d*cd

      eta(1) = p+wid*0.5d0
      eta(2) = p-wid*0.5d0
C      eta(1) = p
C      eta(2) = p-wid
      if (dabs(eta(1)).lt.eps) eta(1) = zro
      if (dabs(eta(2)).lt.eps) eta(2) = zro
      if (dabs(q)     .lt.eps) q      = zro
C----
C Check for singular cases (after Okada DC3D)
C----
      call rectsingular(eq,ek,ee,ksi,eta,q,eps)
      if (eq.eq.1) then
          call singularstrain(uxx,uxy,uxz,uyx,uyy,uyz,uzx,uzy,uzz)
          goto 499
      endif
C----
C Integrate solutions over fault boundaries (-W/2:W/2 and -L/2:L/2) 
C----
      do 408 ii = 1,2
      do 407 jj = 1,2
          if (ii.eq.1.and.jj.eq.1) fac =  one
          if (ii.eq.1.and.jj.eq.2) fac = -one
          if (ii.eq.2.and.jj.eq.1) fac = -one
          if (ii.eq.2.and.jj.eq.2) fac =  one

          call rectvars(ksi(ii),eta(jj),q,z,ek(jj),ee(ii),eps)

          if (thru.eq.0) then
              call xderiv1(fx,ksi(ii),eta(jj),q,z)
              call yderiv1(fy,ksi(ii),eta(jj),q,z)
              call zderiv1(fz,ksi(ii),eta(jj),q,z)
              call disp1stn(f ,ksi(ii),eta(jj),q,z)

              do 405 i = 1,6
                  u(i)  =  u(i) + fac*f(i)
                  do 404 j = 1,3
                      fj(i,j) = fj(i,j) + fac*fx(i,j)
                      fk(i,j) = fk(i,j) + fac*fy(i,j)
                      fl(i,j) = fl(i,j) + fac*fz(i,j)
  404             continue
  405         continue
          else
              call xderiv1(fx,ksi(ii),eta(jj),q,z)
              call yderiv1(fy,ksi(ii),eta(jj),q,z)
              call zderiv1(fz,ksi(ii),eta(jj),q,z)
              do 406 i = 1,6
                  fj(i,1) = fj(i,1) - fac*fx(i,1)
                  fk(i,1) = fk(i,1) - fac*fy(i,1)
                  fl(i,1) = fl(i,1) + fac*fz(i,1)
  406         continue
          endif

  407 continue
  408 continue
C----
C Mirror source
C----
      if (thru.eq.0) then
          z = stdp 
          thru = 1
          goto 403
      else
          z = -stdp
      endif
c----
c If this works, I'm getting a beer. 'Nuff said.
c----
      uxx = Mss*(fj(1,1)+fj(1,2)+z*fj(1,3))
     1    + Mds*(fj(4,1)+fj(4,2)+z*fj(4,3))
      uyx = Mss*((fj(2,1)+fj(2,2)+z*fj(2,3))*cd
     1                               - (fj(3,1)+fj(3,2)+z*fj(3,3))*sd)
     2    + Mds*((fj(5,1)+fj(5,2)+z*fj(5,3))*cd
     3                               - (fj(6,1)+fj(6,2)+z*fj(6,3))*sd)
      uzx = Mss*((fj(2,1)+fj(2,2)-z*fj(2,3))*sd
     1                               + (fj(3,1)+fj(3,2)-z*fj(3,3))*cd)
     2    + Mds*((fj(5,1)+fj(5,2)-z*fj(5,3))*sd
     3                               + (fj(6,1)+fj(6,2)-z*fj(6,3))*cd)
     
      uxy = Mss*(fk(1,1)+fk(1,2)+z*fk(1,3))
     1    + Mds*(fk(4,1)+fk(4,2)+z*fk(4,3))
      uyy = Mss*((fk(2,1)+fk(2,2)+z*fk(2,3))*cd
     1                               - (fk(3,1)+fk(3,2)+z*fk(3,3))*sd)
     2    + Mds*((fk(5,1)+fk(5,2)+z*fk(5,3))*cd
     3                               - (fk(6,1)+fk(6,2)+z*fk(6,3))*sd)
      uzy = Mss*((fk(2,1)+fk(2,2)-z*fk(2,3))*sd
     1                               + (fk(3,1)+fk(3,2)-z*fk(3,3))*cd)
     2    + Mds*((fk(5,1)+fk(5,2)-z*fk(5,3))*sd
     3                               + (fk(6,1)+fk(6,2)-z*fk(6,3))*cd)
     
      uxz = Mss*(fl(1,1)+fl(1,2)+u(1)+z*fl(1,3))
     1    + Mds*(fl(4,1)+fl(4,2)+u(4)+z*fl(4,3))
      uyz = Mss*((fl(2,1)+fl(2,2)+u(2)+z*fl(2,3))*cd 
     1                          - (fl(3,1)+fl(3,2)+u(3)+z*fl(3,3))*sd)
     2    + Mds*((fl(5,1)+fl(5,2)+u(5)+z*fl(5,3))*cd
     3                          - (fl(6,1)+fl(6,2)+u(6)+z*fl(6,3))*sd)
      uzz = Mss*((fl(2,1)+fl(2,2)-u(2)-z*fl(2,3))*sd
     1                          + (fl(3,1)+fl(3,2)-u(3)-z*fl(3,3))*cd)
     2    + Mds*((fl(5,1)+fl(5,2)-u(5)-z*fl(5,3))*sd
     3                          + (fl(6,1)+fl(6,2)-u(6)-z*fl(6,3))*cd)
c----
c Strain tensor from partial derivatives (still in coordinate system
c with +x in direction of strike, +y updip)
c----
  499 continue
      strain(1,1) = uxx
      strain(2,2) = uyy
      strain(3,3) = uzz
      strain(1,2) = (uxy+uyx)*0.5d0
      strain(1,3) = (uxz+uzx)*0.5d0
      strain(2,3) = (uyz+uzy)*0.5d0
      strain(2,1) = strain(1,2)
      strain(3,1) = strain(1,3)
      strain(3,2) = strain(2,3)

      RETURN
      END

C======================================================================C
C----------------------------------------------------------------------C
C-------------------------- SUBROUTINES -------------------------------C
C----------------------------------------------------------------------C
C======================================================================C

      SUBROUTINE DIPVAR(DIPIN)
C----
C PRE-CALCULATE DIP VARIABLES.
C IF (DIP.GT.89.9999D0) DIP = 89.9999D0 (TO AVOID DIVIDE BY ZERO ERROR)
C----
      IMPLICIT NONE
      REAL*8 PI,D2R
      PARAMETER (PI=4.0D0*DATAN(1.0D0),D2R=PI/1.8D2)
      REAL*8 DIP,DIPIN
      REAL*8 SD,CD,S2D,C2D,CDCD,SDSD,CDSD
      COMMON /SOURCE/ SD,CD,S2D,C2D,CDCD,SDSD,CDSD
      IF (DIPIN.GT.89.9999D0) THEN
          DIP = 89.9999D0*D2R
      ELSE
          DIP = DIPIN*D2R
      ENDIF
      SD = DSIN(DIP)
      CD = DCOS(DIP)
      S2D = DSIN(2.0D0*DIP)
      C2D = DCOS(2.0D0*DIP)
      CDCD = CD*CD
      SDSD = SD*SD
      CDSD = SD*CD
      RETURN
      END

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -C

      SUBROUTINE hafspc(vp,vs,dens)
C----
C Calculate half-space constants (Lame parameters, etc.)
C
C Note: in Poisson half-space, vs = vp/sqrt(3), lamda = mu, a = 2/3
C----
      IMPLICIT none
      REAL*8 vp,vs,dens,mu,lamda
      REAL*8 a,CA1,CA2,CB,CC
      COMMON /ELAST/ a,CA1,CA2,CB,CC

      mu = dens*vs*vs
      lamda = dens*vp*vp - 2.0d0*mu
      a   = (lamda+mu)/(lamda+2.0d0*mu)
      CA1 = (1.0d0-a)*0.5d0
      CA2 = a*0.5d0
      CB  = (1.0d0-a)/a
      CC  = 1.0d0-a

      RETURN
      END

C----------------------------------------------------------------------C
C------------------- Source Moment Subroutines ------------------------C
C----------------------------------------------------------------------C

      SUBROUTINE moment0(Mss,Mds,slip,rakin,area)
C----
C Divide slip into strike-slip and dip-slip, and calculate moment of
C each component (POINT SOURCE).
C
C Note: These values are not strictly seismic moments.
C----
      IMPLICIT none
      REAL*8 pi,d2r
      PARAMETER (pi=4.0d0*datan(1.0d0),d2r=pi/1.8d2)
      REAL*8 Mss,Mds,slip,rak,rakin,area

      rak = rakin*d2r
      Mss = slip*dcos(rak)*area/(2.0d0*pi)
      Mds = slip*dsin(rak)*area/(2.0d0*pi)

      RETURN
      END

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -C

      SUBROUTINE moment1(Mss,Mds,slip,rakin)
C----
C Divide slip into strike-slip and dip-slip, and calculate moment of
C each component (FINITE SOURCE).
C
C Note: These values are not strictly seismic moments.
C----
      IMPLICIT none
      REAL*8 pi,d2r
      PARAMETER (pi=4.0d0*datan(1.0d0),d2r=pi/1.8d2)
      REAL*8 Mss,Mds,slip,rak,rakin

      rak = rakin*d2r
      Mss = slip*dcos(rak)/(2.0d0*pi)
      Mds = slip*dsin(rak)/(2.0d0*pi)

      RETURN
      END

C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C------------------- Point Source Subroutines -------------------------C
C----------------------------------------------------------------------C
C----------------------------------------------------------------------C

      SUBROUTINE geomvars0(x,y,c,d)
C----
C Calculate geometric variables from x, y, and d (difference in depth
C between source and station)
C----
      IMPLICIT none
      REAL*8 x,y,c,d

      REAL*8 sd,cd,s2d,c2d,cdcd,sdsd,cdsd
      REAL*8 R,R2,R3,R4,R5,R7,Rd,Rd2,Rd3,Rd4,p,q,s,t
      REAL*8 A3,A5,A7,B3,B5,B7,C3,C5,C7
      REAL*8 I1,I2,I3,I4,I5
      REAL*8 J1,J2,J3,J4
      REAL*8 K1,K2,K3
      REAL*8 U2,V2,W2
      REAL*8 U3,V3,W3
      REAL*8 xx,xy,yy,dd,xd,xc,xq,yq,dq,pq

      COMMON /SOURCE/ sd,cd,s2d,c2d,cdcd,sdsd,cdsd
      COMMON /RVARS0/ R,R2,R3,R4,R5,R7,Rd,p,q,s,t
      COMMON /ABC0/   A3,A5,A7,B3,B5,B7,C3,C5,C7
      COMMON /IVARS0/ I1,I2,I3,I4,I5
      COMMON /JVARS0/ J1,J2,J3,J4
      COMMON /KVARS0/ K1,K2,K3
      COMMON /YVARS0/ U2,V2,W2
      COMMON /ZVARS0/ U3,V3,W3
      COMMON /PRODUCTS/ xx,xy,yy,dd,xd,xc,xq,yq,dq,pq
      
      xx = x*x
      xy = x*y
      yy = y*y
      dd = d*d
      xd = x*d
      xc = x*c

C R = total distance between source and station
      R  = dsqrt(xx+yy+dd)
      R2 = R*R
      R3 = R2*R
      R4 = R3*R
      R5 = R4*R
      R7 = R5*R2

      Rd = R + d
      Rd2 = Rd*Rd
      Rd3 = Rd2*Rd
      Rd4 = Rd3*Rd

C p, q, s, t = distances in rotated reference frame.
C p = fault dip-parallel distance
C q = fault dip-perpendicular distance
      p = y*cd + d*sd
      q = y*sd - d*cd
      s = p*sd + q*cd
      t = p*cd - q*sd

      xq = x*q
      yq = y*q
      dq = d*q
      pq = p*q

C----
C Other variables
C----
      A3 = 1.0d0 - 3.0d0*xx/R2
      A5 = 1.0d0 - 5.0d0*xx/R2
      A7 = 1.0d0 - 7.0d0*xx/R2
      B3 = 1.0d0 - 3.0d0*yy/R2
      B5 = 1.0d0 - 5.0d0*yy/R2 
      B7 = 1.0d0 - 7.0d0*yy/R2
      C3 = 1.0d0 - 3.0d0*dd/R2
      C5 = 1.0d0 - 5.0d0*dd/R2
      C7 = 1.0d0 - 7.0d0*dd/R2

      I1 =  y*(1.0d0/(R*Rd2) - xx*(3.0d0*R+d)/(R3*Rd3))
      I2 =  x*(1.0d0/(R*Rd2) - yy*(3.0d0*R+d)/(R3*Rd3))
      I3 =  x/R3 - I2
      I4 = -xy*(2.0d0*R+d)/(R3*Rd2)
      I5 =  1.0d0/(R*Rd) - xx*(2.0d0*R+d)/(R3*Rd2)

      J1 = -3.0d0*xy*((3.0d0*R+d)/(R3*Rd3)
     1                  - xx*(5.0d0*R2+4.0d0*R*d+dd)/(R5*Rd4))
      J2 = 1.0d0/R3 - 3.0d0/(R*Rd2)
     1         + 3.0d0*xx*yy*(5.0d0*R2+4.0d0*R*d+dd)/(R5*Rd4)
      J3 = A3/R3 - J2
      J4 = -3.0d0*xy/R5 - J1

      K1 = -y*((2.0d0*R+d)/(R3*Rd2)
     1               - xx*(8.0d0*R2+9.0d0*R*d+3.0d0*dd)/(R5*Rd3))
      K2 = -x*((2.0d0*R+d)/(R3*Rd2)
     1               - yy*(8.0d0*R2+9.0d0*R*d+3.0d0*dd)/(R5*Rd3))
      K3 = -3.0d0*x*d/R5 - K2

      U2 = sd - 5.0d0*yq/R2
      V2 = s  - 5.0d0*y*pq/R2
      W2 = sd + U2
      U3 = cd + 5.0d0*dq/R2
      V3 = t  + 5.0d0*d*pq/R2
      W3 = cd + U3

      RETURN
      END

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -C

      SUBROUTINE disp0(f,x,y,c,d)
C----
C Calculate the components of displacement (from table on p. 1025 in
C Okada, 1992) in the array f(6,3).
C----
      IMPLICIT none
      REAL*8 f(6,3),x,y,c,d

      REAL*8 sd,cd,s2d,c2d,cdcd,sdsd,cdsd
      REAL*8 a,CA1,CA2,CB,CC
      REAL*8 R,R2,R3,R4,R5,R7,Rd,p,q,s,t
      REAL*8 A3,A5,A7,B3,B5,B7,C3,C5,C7
      REAL*8 I1,I2,I3,I4,I5
      REAL*8 xx,xy,yy,dd,xd,xc,xq,yq,dq,pq
      INTEGER thru

      COMMON /SOURCE/ sd,cd,s2d,c2d,cdcd,sdsd,cdsd
      COMMON /ELAST/  a,CA1,CA2,CB,CC
      COMMON /RVARS0/ R,R2,R3,R4,R5,R7,Rd,p,q,s,t
      COMMON /ABC0/   A3,A5,A7,B3,B5,B7,C3,C5,C7
      COMMON /IVARS0/ I1,I2,I3,I4,I5
      COMMON /PRODUCTS/ xx,xy,yy,dd,xd,xc,xq,yq,dq,pq
      COMMON /TAG/    thru

C----
C When calculating terms from real source (thru = 0), include all
C components of displacement. Using mirror source (thru = 1), only
C include first column of f.
C----
      if (thru.eq.0) then
          f(1,1) =  CA1*q/R3    + CA2*3.0d0*xx*q/R5
          f(2,1) =  CA1*x*sd/R3 + CA2*3.0d0*xy*q/R5
          f(3,1) = -CA1*x*cd/R3 + CA2*3.0d0*x*d*q/R5
          f(4,1) =                CA2*3.0d0*x*pq/R5
          f(5,1) =  CA1*s/R3    + CA2*3.0d0*y*pq/R5
          f(6,1) = -CA1*t/R3    + CA2*3.0d0*d*pq/R5

          f(1,2) = -3.0d0*xx*q/R5 - CB*I1*sd
          f(2,2) = -3.0d0*xy*q/R5 - CB*I2*sd
          f(3,2) = -3.0d0*xc*q/R5 - CB*I4*sd
          f(4,2) = -3.0d0*x*pq/R5 + CB*I3*cdsd
          f(5,2) = -3.0d0*y*pq/R5 + CB*I1*cdsd
          f(6,2) = -3.0d0*c*pq/R5 + CB*I5*cdsd

          f(1,3) = -CC*A3*cd/R3 + a*3.0d0*c*q*A5/R5 
          f(2,3) =  CC*3.0d0*xy*cd/R5
     1                              + a*3.0d0*xc*(sd-5.0d0*y*q/R2)/R5
          f(3,3) = -CC*3.0d0*xy*sd/R5
     1                              + a*3.0d0*xc*(cd+5.0d0*d*q/R2)/R5
          f(4,3) =  CC*3.0d0*x*t/R5 - a*15.0d0*xc*pq/R7
          f(5,3) = -CC*(c2d-3.0d0*y*t/R2)/R3
     1                                 + a*3.0d0*c*(s-5.0d0*y*pq/R2)/R5
          f(6,3) = -CC*A3*cdsd/R3 + a*3.0d0*c*(t+5.0d0*d*pq/R2)/R5

C----
C Tensile and volume sources
C----
C          f(7,1) =  CA1*x/R3 - CA2*3.0d0*x*q*q/R5
C          f(8,1) =  CA1*t/R3 - CA2*3.0d0*y*q*q/R5
C          f(9,1) =  CA1*s/R3 - CA2*3.0d0*d*q*q/R5
C          f(10,1) = -CA1*x/R3
C          f(11,1) = -CA1*y/R3
C          f(12,1) = -CA1*d/R3
C          f(7,2) = 3.0d0 - CB*I3*sd*sd
C          f(8,2) = 3.0d0 - CB*I1*sd*sd
C          f(9,2) = 3.0d0 - CB*I5*sd*sd
C          f(10,2) =        CB*x/R3
C          f(11,2) =        CB*y/R3
C          f(12,2) =        CB*d/R3
C          f(7,3) = -CC*3.0d0*x*s + a*15.0d0*c*x*q*q/R7 - a*3.0d0*x*z/R5
C          f(8,3) =  CC*(dsin(2.0d0*dip)-3.0d0*y*s/R2)/R3
C     1             + a*3.0d0*(c*(t-y+5.0d0*y*q*q/R2) - y*z)/R5
C          f(9,3) = -CC*(1.0d0-A3*sd*sd)/R3
C     1             - a*3.0d0*(c*(s-d+5.0d0*d*q*q/R2) + d*z)/R5
C          f(10,3) = CC*3.0d0*x*d/R5
C          f(11,3) = CC*3.0d0*y*d/R5
C          f(12,3) = CC*C3/R3
      else
          f(1,1) =  CA1*q/R3    + CA2*3.0d0*xx*q/R5
          f(2,1) =  CA1*x*sd/R3 + CA2*3.0d0*xy*q/R5
          f(3,1) = -CA1*x*cd/R3 + CA2*3.0d0*x*d*q/R5
          f(4,1) =                CA2*3.0d0*x*pq/R5
          f(5,1) =  CA1*s/R3    + CA2*3.0d0*y*pq/R5
          f(6,1) = -CA1*t/R3    + CA2*3.0d0*d*pq/R5
C----
C Tensile and volume sources
C----
C          f(7,1) =  CA1*x/R3 - CA2*3.0d0*x*q*q/R5
C          f(8,1) =  CA1*t/R3 - CA2*3.0d0*y*q*q/R5
C          f(9,1) =  CA1*s/R3 - CA2*3.0d0*d*q*q/R5
C          f(10,1) = -CA1*x/R3
C          f(11,1) = -CA1*y/R3
C          f(12,1) = -CA1*d/R3
      endif


      RETURN
      END

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -C

      SUBROUTINE xderiv0(fx,x,y,c,d)
C----
C Calculate the components of displacement x-derivatives (from table on
C p. 1026 in Okada, 1992) in the array fx(6,3).
C----
      IMPLICIT none
      REAL*8 fx(6,3),x,y,c,d

      REAL*8 sd,cd,s2d,c2d,cdcd,sdsd,cdsd
      REAL*8 a,CA1,CA2,CB,CC
      REAL*8 R,R2,R3,R4,R5,R7,Rd,p,q,s,t
      REAL*8 A3,A5,A7,B3,B5,B7,C3,C5,C7
      REAL*8 J1,J2,J3,J4
      REAL*8 K1,K2,K3
      REAL*8 xx,xy,yy,dd,xd,xc,xq,yq,dq,pq
      INTEGER thru

      COMMON /SOURCE/ sd,cd,s2d,c2d,cdcd,sdsd,cdsd
      COMMON /ELAST/  a,CA1,CA2,CB,CC
      COMMON /RVARS0/ R,R2,R3,R4,R5,R7,Rd,p,q,s,t
      COMMON /ABC0/   A3,A5,A7,B3,B5,B7,C3,C5,C7
      COMMON /JVARS0/ J1,J2,J3,J4
      COMMON /KVARS0/ K1,K2,K3
      COMMON /PRODUCTS/ xx,xy,yy,dd,xd,xc,xq,yq,dq,pq
      COMMON /TAG/    thru

C----
C When calculating terms from real source (thru = 0), include all
C components of displacement. Using mirror source (thru = 1), only
C include first column of fx.
C----
      if (thru.eq.0) then
          fx(1,1) = -CA1*3.0d0*x*q/R5 + CA2*3.0d0*xq*(1.0d0+A5)/R5
          fx(2,1) =  CA1*A3*sd/R3     + CA2*3.0d0*yq*A5/R5
          fx(3,1) = -CA1*A3*cd/R3     + CA2*3.0d0*dq*A5/R5
          fx(4,1) =                     CA2*3.0d0*pq*A5/R5
          fx(5,1) = -CA1*3.0d0*x*s/R5 - CA2*15.0d0*xy*pq/R7
          fx(6,1) =  CA1*3.0d0*x*t/R5 - CA2*15.0d0*x*d*pq/R7

          fx(1,2) = -3.0d0 *xq*(1.0d0+A5)/R5 - CB*J1*sd
          fx(2,2) = -3.0d0 *yq*A5/R5         - CB*J2*sd
          fx(3,2) = -3.0d0 *c*q*A5/R5         - CB*K1*sd
          fx(4,2) = -3.0d0 *pq*A5/R5         + CB*J3*sd*cd
          fx(5,2) =  15.0d0*xy*pq/R7        + CB*J1*sd*cd
          fx(6,2) =  15.0d0*xc*pq/R7        + CB*K3*sd*cd

          fx(1,3) =  CC*3.0d0*x*(2.0d0+A5)*cd/R5
     1                                    - a*15.0d0*c*xq*(2.0d0+A7)/R7
          fx(2,3) =  CC*3.0d0*y*A5*cd/R5
     1                          + a*3.0d0*c*(A5*sd-5.0d0*yq*A7/R2)/R5
          fx(3,3) = -CC*3.0d0*y*A5*sd/R5
     1                          + a*3.0d0*c*(A5*cd+5.0d0*dq*A7/R2)/R5
          fx(4,3) =  CC*3.0d0*t*A5/R5 - a*15.0d0*c*pq*A7/R7
          fx(5,3) =  CC*3.0d0*x*(c2d-5.0d0*y*t/R2)/R5
     1                              - a*15.0d0*xc*(s-7.0d0*y*pq/R2)/R7
          fx(6,3) =  CC*3.0d0*x*(2.0d0+A5)*sd*cd/R5
     1                              - a*15.0d0*xc*(t+7.0d0*d*pq/R2)/R7

C----
C Tensile and volume sources
C----
C          fx(7,1) =  CA1*A3/R3        - CA2*3.0d0*q*q*A5/R5
C          fx(8,1) = -CA1*3.0d0*x*t/R5 + CA2*15.0d0*x*y*q*q/R7
C          fx(9,1) = -CA1*3.0d0*x*s/R5 + CA2*15.0d0*x*d*q*q/R7
C          fx(10,1) = -CA1*A3/R3
C          fx(11,1) =  CA1*3.0d0*x*y/R5
C          fx(12,1) =  CA1*3.0d0*x*d/R5
C          fx(7,2) =   3.0d0*q*q*A5/R5        - CB*J3*sd*sd
C          fx(8,2) = -15.0d0*x*y*q*q/R7       - CB*J1*sd*sd
C          fx(9,2) = -15.0d0*c*x*q*q/R7       - CB*K3*sd*sd
C          fx(10,2) =                           CB*A3/R3
C          fx(11,2) =                          -CB*3.0d0*x*y/R5
C          fx(12,2) =                          -CB*3.0d0*x*d/R5
C          fx(7,3) = -CC*3.0d0*s*A5/R5
C     1                          + a*15.0d0*c*q*q*A7/R7 - a*3.0d0*z*A5/R5
C          fx(8,3) = -CC*3.0d0*x*(s2d-5.0d0*y*s/R2)/R5
C     1        - a*15.0d0*c*x*(t-y+7.0d0*y*q*q/R2)/R7 + a*15.0d0*x*y*z/R7
C          fx(9,3) =  CC*3.0d0*x*(1.0d0-(2.0d0+A5)*sd*sd)/R5
C     1        + a*15.0d0*c*x*(s-d+7.0d0*d*q*q/R2)/R7 - a*15.0d0*x*d*z/R7
C          fx(10,3) =  CC*3.0d0*d*A5/R5
C          fx(11,3) = -CC*15.0d0*x*y*d/R7
C          fx(12,3) = -CC*3.0d0*x*C5/R5
      else
          fx(1,1) = -CA1*3.0d0*x*q/R5 + CA2*3.0d0*xq*(1.0d0+A5)/R5
          fx(2,1) =  CA1*A3*sd/R3     + CA2*3.0d0*yq*A5/R5
          fx(3,1) = -CA1*A3*cd/R3     + CA2*3.0d0*dq*A5/R5
          fx(4,1) =                     CA2*3.0d0*pq*A5/R5
          fx(5,1) = -CA1*3.0d0*x*s/R5 - CA2*15.0d0*xy*pq/R7
          fx(6,1) =  CA1*3.0d0*x*t/R5 - CA2*15.0d0*x*d*pq/R7

C----
C Tensile and volume sources
C----
C          fx(7,1) =  CA1*A3/R3        - CA2*3.0d0*q*q*A5/R5
C          fx(8,1) = -CA1*3.0d0*x*t/R5 + CA2*15.0d0*x*y*q*q/R7
C          fx(9,1) = -CA1*3.0d0*x*s/R5 + CA2*15.0d0*x*d*q*q/R7
C          fx(10,1) = -CA1*A3/R3
C          fx(11,1) =  CA1*3.0d0*x*y/R5
C          fx(12,1) =  CA1*3.0d0*x*d/R5
      endif

      RETURN
      END

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -C

      SUBROUTINE yderiv0(fy,x,y,c,d)
C----
C Calculate the components of displacement y-derivatives (from table on
C p. 1026 in Okada, 1992) in the array fx(6,3).
C----
      IMPLICIT none
      REAL*8 fy(6,3),x,y,c,d

      REAL*8 sd,cd,s2d,c2d,cdcd,sdsd,cdsd
      REAL*8 a,CA1,CA2,CB,CC
      REAL*8 R,R2,R3,R4,R5,R7,Rd,p,q,s,t
      REAL*8 A3,A5,A7,B3,B5,B7,C3,C5,C7
      REAL*8 J1,J2,J3,J4
      REAL*8 K1,K2,K3
      REAL*8 U2,V2,W2
      REAL*8 xx,xy,yy,dd,xd,xc,xq,yq,dq,pq
      INTEGER thru

      COMMON /SOURCE/ sd,cd,s2d,c2d,cdcd,sdsd,cdsd
      COMMON /ELAST/ a,CA1,CA2,CB,CC
      COMMON /RVARS0/ R,R2,R3,R4,R5,R7,Rd,p,q,s,t
      COMMON /ABC0/ A3,A5,A7,B3,B5,B7,C3,C5,C7
      COMMON /JVARS0/ J1,J2,J3,J4
      COMMON /KVARS0/ K1,K2,K3
      COMMON /YVARS0/ U2,V2,W2
      COMMON /PRODUCTS/ xx,xy,yy,dd,xd,xc,xq,yq,dq,pq
      COMMON /TAG/ thru

      if (thru.eq.0) then
          fy(1,1) =  CA1*(sd-3.0d0*y*q/R2)/R3  + CA2*3.0d0*xx*U2/R5
          fy(2,1) = -CA1*3.0d0*xy*sd/R5     + CA2*3.0d0*x*(y*U2+q)/R5
          fy(3,1) =  CA1*3.0d0*xy*cd/R5       + CA2*3.0d0*x*d*U2/R5
          fy(4,1) =                                CA2*3.0d0*x*V2/R5
          fy(5,1) =  CA1*(s2d-3.0d0*y*s/R2)/R3
     1                                         + CA2*3.0d0*(y*V2+pq)/R5
          fy(6,1) = -CA1*(c2d-3.0d0*y*t/R2)/R3 + CA2*3.0d0*d*V2/R5
C          fy(7,1) = -CA1*3.0d0*x*y/R5 - CA2*3.0d0*x*q*W2/R5
C          fy(8,1) =  CA1*(c2d-3.0d0*y*t/R2)/R3 - CA2*3.0d0*(y*q*W2+q*q)/R5
C          fy(9,1) =  CA1*(s2d-3.0d0*y*s/R2)/R3 - CA2*3.0d0*d*q*W2/R5
C          fy(10,1) =  CA1*3.0d0*x*y/R5
C          fy(11,1) = -CA1*B3/R3
C          fy(12,1) =  CA1*3.0d0*y*d/R5

          fy(1,2) = -3.0d0*xx*U2/R5                - CB*J2*sd
          fy(2,2) = -3.0d0*xy*U2/R5 - 3.0d0*x*q/R5 - CB*J4*sd
          fy(3,2) = -3.0d0*xc*U2/R5                - CB*K2*sd
          fy(4,2) = -3.0d0*x*V2/R5                  + CB*J1*sd*cd
          fy(5,2) = -3.0d0*y*V2/R5   - 3.0d0*p*q/R5 + CB*J2*sd*cd
          fy(6,2) = -3.0d0*c*V2/R5                  + CB*K1*sd*cd
C          fy(7,2) = 3.0d0*x*q*W2/R5                - CB*J1*sd*sd
C          fy(8,2) = 3.0d0*y*q*W2/R5 + 3.0d0*q*q/R5 - CB*J2*sd*sd
C          fy(9,2) = 3.0d0*c*q*W2/R5                - CB*K1*sd*sd
C          fy(10,2) =                                -CB*3.0d0*x*y/R5
C          fy(11,2) =                                 CB*B3/R3
C          fy(12,2) =                                -CB*3.0d0*y*d/R5

          fy(1,3) =  CC*3.0d0*y*A5*cd/R5
     1                          + a*3.0d0*c*(A5*sd-5.0d0*y*q*A7/R2)/R5
          fy(2,3) =  CC*3.0d0*x*B5*cd/R5
     1                             - a*15.0d0*xc*(2.0d0*y*sd+q*B7)/R7
          fy(3,3) = -CC*3.0d0*x*B5*sd/R5
     1                           + a*15.0d0*xc*(d*B7*sd-y*C7*cd)/R7
          fy(4,3) =  CC*3.0d0*x*(c2d-5.0d0*y*t/R2)/R5
     1                              - a*15.0d0*xc*(s-7.0d0*y*pq/R2)/R7
          fy(5,3) =  CC*3.0d0*(2.0d0*y*c2d+t*B5)/R5
     1              + a*3.0d0*c*(s2d-10.0d0*y*s/R2-5.0d0*pq*B7/R2)/R5
          fy(6,3) =  CC*3.0d0*y*A5*sd*cd/R5
     1               - a*3.0d0*c*((3.0d0+A5)*c2d+35.0d0*y*d*pq/R4)/R5
C          fy(7,3) = -CC*3.0d0*x*(s2d-5.0d0*y*s/R2)/R5
C     1          - a*15.0d0*c*x*(t-y+7.0d0*y*q*q/R2)/R7 + a*15.0d0*x*y*z/R7
C          fy(8,3) = -CC*3.0d0*(2.0d0*y*s2d+s*B5)/R5
C     1          - a*3.0d0*c*(2.0d0*sd*sd+10.0d0*y*(t-y)/R2
C     2             -5.0d0*q*q*B7/R2)/R5 - a*3.0d0*z*B5/R5
C          fy(9,3) =  CC*3.0d0*y*(1.0d0-A5*sd*sd)/R5
C     1          + a*3.0d0*c*((3.0d0+A5)*s2d-5.0d0*y*d*(2.0d0
C     2             -7.0d0*q*q/R2)/R2)/R5 - a*15.0d0*y*d*z/R7
C          fy(10,3) = -CC*15.0d0*x*y*d/R7
C          fy(11,3) =  CC*3.0d0*d*B5/R5
C          fy(12,3) = -CC*3.0d0*y*C5/R5
      else
          fy(1,1) =  CA1*(sd-3.0d0*y*q/R2)/R3  + CA2*3.0d0*xx*U2/R5
          fy(2,1) = -CA1*3.0d0*xy*sd/R5     + CA2*3.0d0*x*(y*U2+q)/R5
          fy(3,1) =  CA1*3.0d0*xy*cd/R5       + CA2*3.0d0*x*d*U2/R5
          fy(4,1) =                                CA2*3.0d0*x*V2/R5
          fy(5,1) =  CA1*(s2d-3.0d0*y*s/R2)/R3
     1                                         + CA2*3.0d0*(y*V2+pq)/R5
          fy(6,1) = -CA1*(c2d-3.0d0*y*t/R2)/R3 + CA2*3.0d0*d*V2/R5
C          fy(7,1) = -CA1*3.0d0*x*y/R5 - CA2*3.0d0*x*q*W2/R5
C          fy(8,1) =  CA1*(c2d-3.0d0*y*t/R2)/R3 - CA2*3.0d0*(y*q*W2+q*q)/R5
C          fy(9,1) =  CA1*(s2d-3.0d0*y*s/R2)/R3 - CA2*3.0d0*d*q*W2/R5
C          fy(10,1) =  CA1*3.0d0*x*y/R5
C          fy(11,1) = -CA1*B3/R3
C          fy(12,1) =  CA1*3.0d0*y*d/R5
      endif

      RETURN
      END

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -C

      SUBROUTINE zderiv0(fz,x,y,c,d)
C----
C Calculate the components of displacement z-derivatives (from table on
C p. 1026 in Okada, 1992) in the array fx(6,3).
C----
      IMPLICIT none
      REAL*8 fz(6,3),x,y,c,d

      REAL*8 sd,cd,s2d,c2d,cdcd,sdsd,cdsd
      REAL*8 a,CA1,CA2,CB,CC
      REAL*8 R,R2,R3,R4,R5,R7,Rd,p,q,s,t
      REAL*8 A3,A5,A7,B3,B5,B7,C3,C5,C7
      REAL*8 K1,K2,K3
      REAL*8 U3,V3,W3
      REAL*8 xx,xy,yy,dd,xd,xc,xq,yq,dq,pq
      INTEGER thru

      COMMON /SOURCE/ sd,cd,s2d,c2d,cdcd,sdsd,cdsd
      COMMON /ELAST/ a,CA1,CA2,CB,CC
      COMMON /RVARS0/ R,R2,R3,R4,R5,R7,Rd,p,q,s,t
      COMMON /ABC0/ A3,A5,A7,B3,B5,B7,C3,C5,C7
      COMMON /KVARS0/ K1,K2,K3
      COMMON /ZVARS0/ U3,V3,W3
      COMMON /PRODUCTS/ xx,xy,yy,dd,xd,xc,xq,yq,dq,pq
      COMMON /TAG/ thru

      if (thru.eq.0) then
          fz(1,1) =  CA1*(cd+3.0d0*d*q/R2)/R3  + CA2*3.0d0*xx*U3/R5
          fz(2,1) =  CA1*3.0d0*xd*sd/R5       + CA2*3.0d0*xy*U3/R5
          fz(3,1) = -CA1*3.0d0*xd*cd/R5     + CA2*3.0d0*x*(d*U3-q)/R5
          fz(4,1) =                                CA2*3.0d0*x*V3/R5
          fz(5,1) =  CA1*(c2d+3.0d0*d*s/R2)/R3 + CA2*3.0d0*y*V3/R5
          fz(6,1) =  CA1*(s2d-3.0d0*d*t/R2)/R3
     1                                      + CA2*3.0d0*(d*V3-pq)/R5
C          fz(7,1) =  CA1*3.0d0*x*d/R5 - CA2*3.0d0*x*q*W3/R5
C          fz(8,1) = -CA1*(s2d-3.0d0*d*t/R2)/R3 - CA2*3.0d0*y*q*W3/R5
C          fz(9,1) =  CA1*(c2d+3.0d0*d*s/R2)/R3
C     1                                - CA2*3.0d0*q*(d*W3-q)/R5
C          fz(10,1) = -CA1*3.0d0*x*d/R5
C          fz(11,1) = -CA1*3.0d0*y*d/R5
C          fz(12,1) =  CA1*C3/R3

          fz(1,2) = -3.0d0*xx*U3/R5 + CB*K1*sd
          fz(2,2) = -3.0d0*xy*U3/R5 + CB*K2*sd
          fz(3,2) = -3.0d0*xc*U3/R5 + CB*3.d0*xy*sd/R5
          fz(4,2) = -3.0d0*x*V3/R5   - CB*K3*sd*cd
          fz(5,2) = -3.0d0*y*V3/R5   - CB*K1*sd*cd
          fz(6,2) = -3.0d0*c*V3/R5   + CB*A3*sd*cd/R3
C          fz(7,2) = 3.0d0*x*q*W3/R5 + CB*K3*sd*sd
C          fz(8,2) = 3.0d0*y*q*W3/R5 + CB*K1*sd*sd
C          fz(9,2) = 3.0d0*c*q*W3/R5 - CB*A3*sd*sd/R3
C          fz(10,2) =  CB*3.0d0*x*d/R5
C          fz(11,2) =  CB*3.0d0*y*d/R5
C          fz(12,2) = -CB*C3/R3

          fz(1,3) = -CC*3.0d0*d*A5*cd/R5
     1                          + a*3.0d0*c*(A5*cd+5.0d0*d*q*A7/R2)/R5
          fz(2,3) =  CC*15.0d0*xy*d*cd/R7
     1                           + a*15.0d0*xc*(d*B7*sd-y*C7*cd)/R7
          fz(3,3) = -CC*15.0d0*xy*d*sd/R7
     1                             + a*15.0d0*xc*(2.0d0*d*cd-q*C7)/R7
          fz(4,3) = -CC*3.0d0*x*(s2d-5.0d0*d*t/R2)/R5
     1                              - a*15.0d0*c*x*(t+7.0d0*d*pq/R2)/R7
          fz(5,3) = -CC*3.0d0*(d*B5*c2d+y*C5*s2d)/R5 
     1               - a*3.0d0*c*((3.0d0+A5)*c2d+35.0d0*y*d*pq/R4)/R5
          fz(6,3) = -CC*3.0d0*d*A5*sd*cd/R5
     1               - a*3.0d0*c*(s2d-(10.0d0*d*t-5.0d0*pq*C7)/R2)/R5
C          fz(7,3) = -CC*3.0d0*x*(c2d+5.0d0*d*s/R2)/R5
C     1       + a*15.0d0*c*x*(s-d+7.0d0*d*q*q/R2)/R7
C     2       - a*3.0d0*x*(1.0d0+5.0d0*d*z/R2)/R5
C          fz(8,3) =  CC*3.0d0*(d*B5*s2d-y*C5*c2d)/R5
C     1       + a*3.0d0*c*((3.0d0+A5)*s2d-5.0d0*y*d*(2.0d0
C     2       -7.0d0*q*q/R2)/R2)/R5 - a*3.0d0*y*(1.0d0+5.0d0*d*z/R2)/R5
C          fz(9,3) = -CC*3.0d0*d*(1.0d0-A5*sd*sd)/R5
C     1       - a*3.0d0*c*(c2d+10.0d0*d*(s-d)/R2-5.0d0*q*q*C7/R2)/R5
C     2       - a*3.0d0*z*(1.0d0+C5)/R5
C          fz(10,3) = -CC*3.0d0*x*C5/R5
C          fz(11,3) = -CC*3.0d0*y*C5/D5
C          fz(12,3) =  CC*3.0d0*d*(2.0d0+C5)/R5
      else
          fz(1,1) =  CA1*(cd+3.0d0*d*q/R2)/R3  + CA2*3.0d0*xx*U3/R5
          fz(2,1) =  CA1*3.0d0*xd*sd/R5       + CA2*3.0d0*xy*U3/R5
          fz(3,1) = -CA1*3.0d0*xd*cd/R5     + CA2*3.0d0*x*(d*U3-q)/R5
          fz(4,1) =                                CA2*3.0d0*x*V3/R5
          fz(5,1) =  CA1*(c2d+3.0d0*d*s/R2)/R3 + CA2*3.0d0*y*V3/R5
          fz(6,1) =  CA1*(s2d-3.0d0*d*t/R2)/R3
     1                                      + CA2*3.0d0*(d*V3-pq)/R5
C          fz(7,1) =  CA1*3.0d0*x*d/R5 - CA2*3.0d0*x*q*W3/R5
C          fz(8,1) = -CA1*(s2d-3.0d0*d*t/R2)/R3 - CA2*3.0d0*y*q*W3/R5
C          fz(9,1) =  CA1*(c2d+3.0d0*d*s/R2)/R3
C     1                                - CA2*3.0d0*q*(d*W3-q)/R5
C          fz(10,1) = -CA1*3.0d0*x*d/R5
C          fz(11,1) = -CA1*3.0d0*y*d/R5
C          fz(12,1) =  CA1*C3/R3
      endif

      RETURN
      END

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -C

      SUBROUTINE disp0stn(f,x,y,c,d)
      IMPLICIT none
      REAL*8 f(6),x,y,c,d

      REAL*8 sd,cd,s2d,c2d,cdcd,sdsd,cdsd
      REAL*8 a,CA1,CA2,CB,CC
      REAL*8 R,R2,R3,R4,R5,R7,Rd,p,q,s,t
      REAL*8 A3,A5,A7,B3,B5,B7,C3,C5,C7
      REAL*8 xx,xy,yy,dd,xd,xc,xq,yq,dq,pq

      COMMON /SOURCE/ sd,cd,s2d,c2d,cdcd,sdsd,cdsd
      COMMON /ELAST/ a,CA1,CA2,CB,CC
      COMMON /RVARS0/ R,R2,R3,R4,R5,R7,Rd,p,q,s,t
      COMMON /ABC0/ A3,A5,A7,B3,B5,B7,C3,C5,C7
      COMMON /PRODUCTS/ xx,xy,yy,dd,xd,xc,xq,yq,dq,pq


           f(1) = -CC*A3*cd/R3 + a*3.0d0*c*q*A5/R5
           f(2) =  CC*3.0d0*xy*cd/R5
     1                              + a*3.0d0*xc*(sd-5.0d0*yq/R2)/R5
           f(3) = -CC*3.0d0*xy*sd/R5
     1                              + a*3.0d0*xc*(cd+5.0d0*dq/R2)/R5
           f(4) =  CC*3.0d0*x*t/R5 - a*15.0d0*xc*pq/R7
           f(5) = -CC*(c2d-3.0d0*y*t/R2)/R3
     1                                   + a*3.0d0*c*(s-5.0*y*pq/R2)/R5
           f(6) = -CC*A3*sd*cd/R3  + a*3.0d0*c*(t+5.0*d*pq/R2)/R5

      RETURN
      END


C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C------------------- Finite Source Subroutines ------------------------C
C----------------------------------------------------------------------C
C----------------------------------------------------------------------C

      SUBROUTINE rectvars(ksi,eta,q,z,ek,ee,eps)
C----
C Calculate geometric variables needed for displacements and strains
C due to finite rectangular source.
C----     
      IMPLICIT NONE
      REAL*8 sd,cd,s2d,c2d,cdcd,sdsd,cdsd
      REAL*8 a,CA1,CA2,CB,CC
      REAL*8 cbar,dbar,ybar
      REAL*8 R,R2,R3,R5,X11,X32,X53,Y11,Y32,Y53,Y0,Z32,Z53,Z0
      REAL*8 I1,I2,I3,I4
      REAL*8 J1,J2,J3,J4,J5,J6
      REAL*8 K1,K2,K3,K4
      REAL*8 E2,F2,G2,H2,P2,Q2
      REAL*8 E3,F3,G3,H3,P3,Q3
      REAL*8 D11,Rd,Re,Rk,logRe,logRk,TH
      REAL*8 h,XX
      REAL*8 ksi,eta,q,z,eps
      INTEGER ee,ek

      COMMON /SOURCE/ sd,cd,s2d,c2d,cdcd,sdsd,cdsd
      COMMON /ELAST/ a,CA1,CA2,CB,CC
      COMMON /BAR/ cbar,dbar,ybar
      COMMON /RVARS/ R,R2,R3,R5,X11,X32,X53,Y11,Y32,Y53,Y0,Z32,Z53,Z0
      COMMON /IVARS/ I1,I2,I3,I4
      COMMON /JVARS/ J1,J2,J3,J4,J5,J6
      COMMON /KVARS/ K1,K2,K3,K4
      COMMON /YVARS/ E2,F2,G2,H2,P2,Q2
      COMMON /ZVARS/ E3,F3,G3,H3,P3,Q3
      COMMON /MISC/ D11,Rd,Re,Rk,logRe,logRk,TH
C----
      ybar = eta*cd + q*sd
      dbar = eta*sd - q*cd
      cbar = dbar + z

      R  = dsqrt(ksi*ksi+eta*eta+q*q)
      R2 = R*R
      R3 = R2*R
      R5 = R3*R2
      Rd = R+dbar
      h = q*cd - z
C----
      if (dabs(q).lt.eps) then
          TH = 0.0d0
      else
          TH = datan(ksi*eta/(q*R))
      endif
C----
      if (ek.eq.1) then
          X11 = 0.0d0
          X32 = 0.0d0
          X53 = 0.0d0
          logRk = -dlog(R-ksi)
      else
          Rk = R+ksi
          X11 = 1.0d0/(R*Rk)
          X32 = (2.0d0*R+ksi)/(R3*Rk*Rk)
          X53 = (8.0d0*R2+9.0d0*R*ksi+3.0d0*ksi*ksi)/(R5*Rk*Rk*Rk)
          logRk = dlog(R+ksi)
      endif

      if (ee.eq.1) then
          Y11 = 0.0d0
          Y32 = 0.0d0
          Y53 = 0.0d0
          logRe = -dlog(R-eta)
      else
          Re = R+eta
          Y11 = 1.0d0/(R*Re)
          Y32 = (2.0d0*R+eta)/(R3*Re*Re)
          Y53 = (8.0d0*R2+9.0d0*R*eta+3.0d0*eta*eta)/(R5*Re*Re*Re)
          logRe = dlog(R+eta)
      endif
      Y0  = Y11 - ksi*ksi*Y32
      
      Z32 = sd/R3 - h*Y32
      Z53 = 3.0d0*sd/R5 - h*Y53
      Z0  = Z32 - ksi*ksi*Z53
C----
      XX = sqrt(ksi*ksi+q*q)

      I3 = ybar/(cd*Rd) - (logRe-sd*dlog(Rd))/(cdcd)
      if (dabs(ksi).lt.eps) then
          I4 = 0.0d0
      else
      I4 = sd*ksi/(Rd*cd)
     1      + 2.0d0*datan((eta*(XX+q*cd)+XX*(R+XX)*sd)
     2                               /(ksi*(R+XX)*cd))/(cdcd)
      endif
      I1 = -ksi*cd/Rd - I4*sd
      I2 = dlog(Rd) + I3*sd
C----
      D11 = 1.0d0/(R*Rd)

      K1 = ksi*(D11-Y11*sd)/cd
      K3 = (q*Y11-ybar*D11)/cd
      K2 = 1.0d0/R+K3*sd
      K4 = ksi*Y11*cd-K1*sd
      
      J2 = ksi*ybar*D11/Rd
      J5 = -(dbar+ybar*ybar/Rd)*D11
      J3 = (K1-J2*sd)/cd
      J6 = (K3-J5*sd)/cd
      J1 = J5*cd - J6*sd
      J4 = -ksi*Y11 - J2*cd + J3*sd
C----
      E2 = sd/R - ybar*q/R3
      F2 = dbar/R3 + ksi*ksi*Y32*sd
      G2 = 2.0d0*X11*sd - ybar*q*X32
      H2 = dbar*q*X32 + ksi*q*Y32*sd
      P2 = cd/R3 + q*Y32*sd
      Q2 = 3.0d0*cbar*dbar/R5 - (z*Y32+Z32+Z0)*sd

      E3 = cd/R + dbar*q/R3
      F3 = ybar/R3 + ksi*ksi*Y32*cd
      G3 = 2.0d0*X11*cd + dbar*q*X32
      H3 = ybar*q*X32 + ksi*q*Y32*cd
      P3 = sd/R3 - q*Y32*cd
      Q3 = 3.0d0*cbar*ybar/R5 + q*Y32 - (z*Y32+Z32+Z0)*cd

      RETURN
      END

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -C

      SUBROUTINE disp1(f,ksi,eta,q,z)
C----
C Components of displacement from finite source
C----
      IMPLICIT none
      REAL*8 f(6,3),ksi,eta,q,z
      REAL*8 sd,cd,s2d,c2d,cdcd,sdsd,cdsd
      REAL*8 a,CA1,CA2,CB,CC
      REAL*8 cbar,dbar,ybar
      REAL*8 R,R2,R3,R5,X11,X32,X53,Y11,Y32,Y53,Y0,Z32,Z53,Z0
      REAL*8 I1,I2,I3,I4
      REAL*8 D11,Rd,Re,Rk,logRe,logRk,TH
      INTEGER thru

      COMMON /SOURCE/ sd,cd,s2d,c2d,cdcd,sdsd,cdsd
      COMMON /ELAST/ a,CA1,CA2,CB,CC
      COMMON /BAR/ cbar,dbar,ybar
      COMMON /RVARS/ R,R2,R3,R5,X11,X32,X53,Y11,Y32,Y53,Y0,Z32,Z53,Z0
      COMMON /IVARS/ I1,I2,I3,I4
      COMMON /MISC/ D11,Rd,Re,Rk,logRe,logRk,TH
      COMMON /TAG/ thru

      if (thru.eq.0) then
      f(1,1) = TH*0.5d0    + CA2*ksi*q*Y11
      f(2,1) =               CA2*q/R
      f(3,1) = CA1*logRe   - CA2*q*q*Y11
      f(4,1) =               CA2*q/R
      f(5,1) = TH*0.5d0    + CA2*eta*q*X11
      f(6,1) = CA1*logRk   - CA2*q*q*X11

      f(1,2) = -ksi*q*Y11 - TH - CB*I1*sd
      f(2,2) = -q/R            + CB*ybar*sd/Rd
      f(3,2) =  q*q*Y11        - CB*I2*sd
      f(4,2) = -q/R            + CB*I3*cdsd
      f(5,2) = -eta*q*X11 - TH - CB*ksi*cdsd/Rd
      f(6,2) =  q*q*X11        + CB*I4*cdsd

      f(1,3) = CC*ksi*Y11*cd - a*ksi*q*Z32
      f(2,3) = CC*(cd/R+2.0d0*q*Y11*sd) - a*cbar*q/R3
      f(3,3) = CC*q*Y11*cd - a*(cbar*eta/R3-z*Y11+ksi*ksi*Z32)
      f(4,3) = CC*cd/R - q*Y11*sd - a*cbar*q/R3
      f(5,3) = CC*ybar*X11 - a*cbar*eta*q*X32
      f(6,3) = -dbar*X11 - ksi*Y11*sd - a*cbar*(X11 - q*q*X32)
      else
      f(1,1) = TH*0.5d0    + CA2*ksi*q*Y11
      f(2,1) =               CA2*q/R
      f(3,1) = CA1*logRe   - CA2*q*q*Y11
      f(4,1) =               CA2*q/R
      f(5,1) = TH*0.5d0    + CA2*eta*q*X11
      f(6,1) = CA1*logRk   - CA2*q*q*X11
      endif
      RETURN
      END

C-------------------

      SUBROUTINE xderiv1(fx,ksi,eta,q,z)
C----
C Components of x-derivatives of displacement from finite source
C----
      IMPLICIT none
      REAL*8 fx(6,3),ksi,eta,q,z
      REAL*8 sd,cd,s2d,c2d,cdcd,sdsd,cdsd
      REAL*8 a,CA1,CA2,CB,CC
      REAL*8 cbar,dbar,ybar
      REAL*8 R,R2,R3,R5,X11,X32,X53,Y11,Y32,Y53,Y0,Z32,Z53,Z0
      REAL*8 J1,J2,J3,J4,J5,J6
      INTEGER thru

      COMMON /SOURCE/ sd,cd,s2d,c2d,cdcd,sdsd,cdsd
      COMMON /ELAST/ a,CA1,CA2,CB,CC
      COMMON /BAR/ cbar,dbar,ybar
      COMMON /RVARS/ R,R2,R3,R5,X11,X32,X53,Y11,Y32,Y53,Y0,Z32,Z53,Z0
      COMMON /JVARS/ J1,J2,J3,J4,J5,J6
      COMMON /TAG/ thru

      if (thru.eq.0) then
      fx(1,1) = -CA1*q*Y11   - CA2*ksi*ksi*q*Y32
      fx(2,1) =              - CA2*ksi*q/R3
      fx(3,1) =  CA1*ksi*Y11 + CA2*ksi*q*q*Y32
      fx(4,1) =              - CA2*ksi*q/R3
      fx(5,1) = -q*Y11*0.5d0 - CA2*eta*q/R3
      fx(6,1) =  CA1/R       + CA2*q*q/R3

      fx(1,2) =  ksi*ksi*q*Y32     - CB*J1*sd
      fx(2,2) =  ksi*q/R3          - CB*J2*sd
      fx(3,2) = -ksi*q*q*Y32       - CB*J3*sd
      fx(4,2) =  ksi*q/R3          + CB*J4*cdsd
      fx(5,2) =  eta*q/R3 + q*Y11  + CB*J5*cdsd
      fx(6,2) = -q*q/R3            + CB*J6*cdsd

      fx(1,3) =  CC*Y0*cd - a*q*Z0
      fx(2,3) = -CC*ksi*(cd/R3+2.0d0*q*Y32*sd)
     1                                           + a*3.0d0*cbar*ksi*q/R5
      fx(3,3) = -CC*ksi*q*Y32*cd
     1                          + a*ksi*(3.0d0*cbar*eta/R5-z*Y32-Z32-Z0)
      fx(4,3) = -CC*ksi*cd/R3 + ksi*q*Y32*sd + a*3.0d0*cbar*ksi*q/R5
      fx(5,3) = -CC*ybar/R3 + a*3.0d0*cbar*eta*q/R5
      fx(6,3) = dbar/R3 - Y0*sd + a*cbar*(1.0d0-3.0d0*q*q/R2)/R3
      else
      fx(1,1) = -CA1*q*Y11   - CA2*ksi*ksi*q*Y32
      fx(2,1) =              - CA2*ksi*q/R3
      fx(3,1) =  CA1*ksi*Y11 + CA2*ksi*q*q*Y32
      fx(4,1) =              - CA2*ksi*q/R3
      fx(5,1) = -q*Y11*0.5d0 - CA2*eta*q/R3
      fx(6,1) =  CA1/R       + CA2*q*q/R3
      endif

      RETURN
      END

C-------------------

      SUBROUTINE yderiv1(fy,ksi,eta,q,z)
C----
C Components of y-derivatives of displacement from finite source
C----
      IMPLICIT none
      REAL*8 fy(6,3),ksi,eta,q,z
      REAL*8 sd,cd,s2d,c2d,cdcd,sdsd,cdsd
      REAL*8 a,CA1,CA2,CB,CC
      REAL*8 cbar,dbar,ybar
      REAL*8 R,R2,R3,R5,X11,X32,X53,Y11,Y32,Y53,Y0,Z32,Z53,Z0
      REAL*8 J1,J2,J3,J4,J5,J6
      REAL*8 E2,F2,G2,H2,P2,Q2
      INTEGER thru

      COMMON /SOURCE/ sd,cd,s2d,c2d,cdcd,sdsd,cdsd
      COMMON /ELAST/ a,CA1,CA2,CB,CC
      COMMON /BAR/ cbar,dbar,ybar
      COMMON /RVARS/ R,R2,R3,R5,X11,X32,X53,Y11,Y32,Y53,Y0,Z32,Z53,Z0
      COMMON /JVARS/ J1,J2,J3,J4,J5,J6
      COMMON /YVARS/ E2,F2,G2,H2,P2,Q2
      COMMON /TAG/ thru

      if (thru.eq.0) then
      fy(1,1) =  CA1*ksi*Y11*sd + dbar*X11*0.5d0 + CA2*ksi*F2
      fy(2,1) =                                    CA2*E2
      fy(3,1) =  CA1*(cd/R+q*Y11*sd)         - CA2*q*F2
      fy(4,1) =                                    CA2*E2
      fy(5,1) =  CA1*dbar*X11 + ksi*Y11*sd*0.5d0 + CA2*eta*G2
      fy(6,1) =  CA1*ybar*X11                    - CA2*q*G2

      fy(1,2) = -ksi*F2 - dbar*X11     + CB*(ksi*Y11+J4)*sd
      fy(2,2) = -E2                    + CB*(1.0d0/R+J5)*sd
      fy(3,2) =  q*F2                  - CB*(q*Y11-J6)*sd
      fy(4,2) = -E2                    + CB*J1*cdsd
      fy(5,2) = -eta*G2 - ksi*Y11*sd + CB*J2*cdsd
      fy(6,2) =  q*G2                  + CB*J3*cdsd

      fy(1,3) = -CC*ksi*P2*cd - a*ksi*Q2
      fy(2,3) =  2.0d0*CC*(dbar/R3-Y0*sd)*sd - ybar*cd/R3
     1             - a*((cbar+dbar)*sd/R3-eta/R3-3.0d0*cbar*ybar*q/R5)
      fy(3,3) = -CC*q/R3 + (ybar/R3-Y0*cd)*sd
     1                    + a*((cbar+dbar)*cd/R3+3.0d0*cbar*dbar*q/R5
     2                                             -(Y0*cd+q*Z0)*sd)
      fy(4,3) = -CC*eta/R3 + Y0*sdsd
     1                    - a*((cbar+dbar)*sd/R3-3.0d0*cbar*ybar*q/R5)
      fy(5,3) =  CC*(X11-ybar*ybar*X32)
     1                 - a*cbar*((dbar+2.0d0*q*cd)*X32-ybar*eta*q*X53)
      fy(6,3) =  ksi*P2*sd + ybar*dbar*X32
     1                   + a*cbar*((ybar+2.0d0*q*sd)*X32-ybar*q*q*X53)
      else
      fy(1,1) =  CA1*ksi*Y11*sd + dbar*X11*0.5d0 + CA2*ksi*F2
      fy(2,1) =                                    CA2*E2
      fy(3,1) =  CA1*(cd/R+q*Y11*sd)         - CA2*q*F2
      fy(4,1) =                                    CA2*E2
      fy(5,1) =  CA1*dbar*X11 + ksi*Y11*sd*0.5d0 + CA2*eta*G2
      fy(6,1) =  CA1*ybar*X11                    - CA2*q*G2
      endif
      RETURN
      END

C-------------------

      SUBROUTINE zderiv1(fz,ksi,eta,q,z)
C----
C Components of z-derivatives of displacement from finite source
C----
      IMPLICIT none
      REAL*8 fz(6,3),ksi,eta,q,z
      REAL*8 sd,cd,s2d,c2d,cdcd,sdsd,cdsd
      REAL*8 a,CA1,CA2,CB,CC
      REAL*8 cbar,dbar,ybar
      REAL*8 R,R2,R3,R5,X11,X32,X53,Y11,Y32,Y53,Y0,Z32,Z53,Z0
      REAL*8 K1,K2,K3,K4
      REAL*8 E3,F3,G3,H3,P3,Q3
      REAL*8 D11,Rd,Re,Rk,logRe,logRk,TH
      INTEGER thru

      COMMON /SOURCE/ sd,cd,s2d,c2d,cdcd,sdsd,cdsd
      COMMON /ELAST/ a,CA1,CA2,CB,CC
      COMMON /BAR/ cbar,dbar,ybar
      COMMON /RVARS/ R,R2,R3,R5,X11,X32,X53,Y11,Y32,Y53,Y0,Z32,Z53,Z0
      COMMON /KVARS/ K1,K2,K3,K4
      COMMON /ZVARS/ E3,F3,G3,H3,P3,Q3
      COMMON /MISC/ D11,Rd,Re,Rk,logRe,logRk,TH
      COMMON /TAG/ thru

      if (thru.eq.0) then
      fz(1,1) =  CA1*ksi*Y11*cd + ybar*X11*0.5d0 + CA2*ksi*F3
      fz(2,1) =                                    CA2*E3
      fz(3,1) = -CA1*(sd/R-q*Y11*cd)         - CA2*q*F3
      fz(4,1) =                                    CA2*E3
      fz(5,1) =  CA1*ybar*X11 + ksi*Y11*cd*0.5d0 + CA2*eta*G3
      fz(6,1) = -CA1*dbar*X11                    - CA2*q*G3

      fz(1,2) = -ksi*F3 - ybar*X11     + CB*K1*sd
      fz(2,2) = -E3                    + CB*ybar*D11*sd
      fz(3,2) =  q*F3                  + CB*K2*sd
      fz(4,2) = -E3                    - CB*K3*cdsd
      fz(5,2) = -eta*G3 - ksi*Y11*cd - CB*ksi*D11*cdsd
      fz(6,2) =  q*G3                  - CB*K4*cdsd

      fz(1,3) = CC*ksi*P3*cd - a*ksi*Q3
      fz(2,3) = 2.0d0*CC*(ybar/R3-Y0*cd)*sd + dbar*cd/R3
     1                    - a*((cbar+dbar)*cd/R3+3.0d0*cbar*dbar*q/R5)
      fz(3,3) = (ybar/R3-Y0*cd)*cd
     1                    - a*((cbar+dbar)*sd/R3-3.0d0*cbar*ybar*q/R5
     2                                          -Y0*sd*sd+q*Z0*cd)
      fz(4,3) = -q/R3 + Y0*sd*cd
     1                    - a*((cbar+dbar)*cd/R3+3.0d0*cbar*dbar*q/R5)
      fz(5,3) = CC*ybar*dbar*X32
     1                 - a*cbar*((ybar-2.0d0*q*sd)*X32+dbar*eta*q*X53)
      fz(6,3) = -ksi*P3*sd + X11 - dbar*dbar*X32
     1                     - a*cbar*((dbar-2.0*q*cd)*X32-dbar*q*q*X53)
      else
      fz(1,1) =  CA1*ksi*Y11*cd + ybar*X11*0.5d0 + CA2*ksi*F3
      fz(2,1) =                                    CA2*E3
      fz(3,1) = -CA1*(sd/R-q*Y11*cd)         - CA2*q*F3
      fz(4,1) =                                    CA2*E3
      fz(5,1) =  CA1*ybar*X11 + ksi*Y11*cd*0.5d0 + CA2*eta*G3
      fz(6,1) = -CA1*dbar*X11                    - CA2*q*G3
      endif

      RETURN
      END

C-------------------

      SUBROUTINE disp1stn(f,ksi,eta,q,z)
C----
C Components of displacement due to finite source for calculating
C partial derivatives
C----
      IMPLICIT none
      REAL*8 f(6),ksi,eta,q,z
      REAL*8 sd,cd,s2d,c2d,cdcd,sdsd,cdsd
      REAL*8 a,CA1,CA2,CB,CC
      REAL*8 cbar,dbar,ybar
      REAL*8 R,R2,R3,R5,X11,X32,X53,Y11,Y32,Y53,Y0,Z32,Z53,Z0

      COMMON /SOURCE/ sd,cd,s2d,c2d,cdcd,sdsd,cdsd
      COMMON /ELAST/ a,CA1,CA2,CB,CC
      COMMON /BAR/ cbar,dbar,ybar
      COMMON /RVARS/ R,R2,R3,R5,X11,X32,X53,Y11,Y32,Y53,Y0,Z32,Z53,Z0

      f(1) = CC*ksi*Y11*cd - a*ksi*q*Z32
      f(2) = CC*(cd/R+2.0d0*q*Y11*sd) - a*cbar*q/R3
      f(3) = CC*q*Y11*cd - a*(cbar*eta/R3-z*Y11+ksi*ksi*Z32)
      f(4) = CC*cd/R - q*Y11*sd - a*cbar*q/R3
      f(5) = CC*ybar*X11 - a*cbar*eta*q*X32
      f(6) = -dbar*X11 - ksi*Y11*sd - a*cbar*(X11 - q*q*X32)

      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE rectsingular(eq,ek,ee,ksi,eta,q,eps)
C----
C Check to see if observation point lies on singular locations
C----
      IMPLICIT none
      REAL*8 q,ksi(2),eta(2),eps,R12,R22,R21
      INTEGER eq,ek(2),ee(2)

      eq = 0
      ek(1) = 0
      ek(2) = 0
      ee(1) = 0
      ee(2) = 0
      
      if (q.eq.0.0d0 .and.
     1   ((ksi(1)*ksi(2).le.0.0d0.and.eta(1)*eta(2).eq.0.0d0).or.
     2    (eta(1)*eta(2).le.0.0d0.and.ksi(1)*ksi(2).eq.0.0d0))) eq = 1

      R12 = dsqrt(ksi(1)*ksi(1) + eta(2)*eta(2) + q*q)
      R22 = dsqrt(ksi(2)*ksi(2) + eta(2)*eta(2) + q*q)
      R21 = dsqrt(ksi(2)*ksi(2) + eta(1)*eta(1) + q*q)
      if (ksi(1).lt.0.0d0.and.R21+ksi(2).lt.eps) ek(1) = 1
      if (ksi(1).lt.0.0d0.and.R22+ksi(2).lt.eps) ek(2) = 1
      if (eta(1).lt.0.0d0.and.R12+eta(2).lt.eps) ee(1) = 1
      if (eta(1).lt.0.0d0.and.R22+eta(2).lt.eps) ee(2) = 1

      RETURN
      END

C----------------------------------------------------------------------C
C If observation point lies at point source or on finite fault edge, 
C return zeros.
C----
      SUBROUTINE singularstrain(uxx,uxy,uxz,uyx,uyy,uyz,uzx,uzy,uzz)
      IMPLICIT none
      REAL*8 uxx,uxy,uxz,uyx,uyy,uyz,uzx,uzy,uzz

      print *,'Observation point is singular.'
      uxx = 0.0d0
      uxy = 0.0d0
      uxz = 0.0d0
      uyx = 0.0d0
      uyy = 0.0d0
      uyz = 0.0d0
      uzx = 0.0d0
      uzy = 0.0d0
      uzz = 0.0d0

      RETURN
      END

      SUBROUTINE singulardisplacement(ux,uy,uz)
      IMPLICIT none
      REAL*8 ux,uy,uz

      print *,'Observation point is singular.'
      ux = 0.0d0
      uy = 0.0d0
      uz = 0.0d0

      RETURN
      END

