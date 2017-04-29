C======================================================================C
C Subroutines to compute deformation from point and finite sources in
C an elastic half-space.
C
C Reference:
C Okada, Y., 1992, Internal deformation due to shear and tensile
C   faults in a half-space. Bulletin of the Seismological Society of
C   America 82, pp. 1018-1040.
C
C Notes:
C - All variables are assumed to be in SI units
C - Angles input into subroutines are in degrees
C - Coordinate system is relative to the fault orientation:
C   - X points in the along-strike direction
C   - Y points in the horizontal up-dip direction
C   - Depths are defined positive DOWN
C======================================================================C

      SUBROUTINE o92ptvol(ux,uy,uz,x,y,stdp,evdp,dvol,vp,vs,dens)
C----
C Compute displacement vector at a point in anywhere in a half-space
C due to a point source
C----

      IMPLICIT none
      REAL*8 x,y,stdp,evdp,z,c,d             ! Coordinates
      REAL*8 dipin,dvol,Mvol      ! Fault parameters
      REAL*8 vp,vs,dens                      ! Half-space parameters
      INTEGER i,j
      REAL*8 u(6,3),f(6,3),ux,uy,uz
C      REAL*8 disp(3)
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
C If station is above surface, set displacements to zero and exit
C----
      if (stdp.lt.zro) then
C          disp(1) = zro
C          disp(2) = zro
C          disp(3) = zro
          ux = zro
          uy = zro
          uz = zro
          return
      endif

C----
C Initialize components of displacement, u(i,j):
C   i = 1-3: components from unit strike-slip source
C   i = 4-6: components from unit dip-slip source
C   i = 7-9: components from unit tensile source (NOT FINISHED)
C   i = 10-12: components from unit inflation source (NOT FINISHED)
C
C thru - a flag indicating number of times through calculation;
C          used for generating mirror source
C----
      thru = 0
      do 102 i = 1,6
          do 101 j = 1,3
C              ut(i,j) = zro
              u(i,j) = zro
  101     continue
  102 continue

C----
C Source and halfspace constants
C----
      dipin = 1.0d0
      call dipvar(dipin)
      call hafspc(vp,vs,dens)
      call moment0vol(Mvol,dvol,vp,vs,dens)

C----
C Calculate displacement components (done in sub disp0vol)
C   x in direction of strike
C   y in updip direction (horizontal)
C   z vertical
C----
      c = evdp
      z = -stdp
  103 d = c - z

      call geomvars0(x,y,c,d)

C f(i,j): dummy displacement components
C
      if (thru.eq.0) then
          call disp0vol(f,x,y,c,d)
          do 105 i = 1,6
              do 104 j = 1,3
C                  ut(i,j) = ut(i,j) + f(i,j)
                  u(i,j) = u(i,j) + f(i,j)
  104         continue
  105     continue
          thru = 1
          z = stdp
          goto 103
      else
          call disp0vol(f,x,y,c,d)
          do 106 i = 1,6
C              ut(i,1) = ut(i,1) - f(i,1)
              u(i,1) = u(i,1) - f(i,1)
  106     continue
          z = -stdp
      endif

C----
C Double double toil and trouble...put it all together...and....PRESTO!
C Combine displacement components to get displacements at the receiver
C (+x in strike direction; +y in updip horizontal direction)
C----
      ux = Mvol*(u(4,1)+u(4,2)+z*u(4,3))
      uy = Mvol*(u(5,1)+u(5,2)+z*u(5,3))
      uz = Mvol*(u(6,1)+u(6,2)+z*u(6,3))

      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE o92ptstnvol(strain,x,y,stdp,evdp,dvol,vp,vs,dens)
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
      REAL*8 dipin,dvol,Mvol
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

      if (stdp.lt.0.0d0) then
          do 315 i = 1,3
              do 316 j = 1,3
                  strain(i,j) = 0.0d0
 316          continue
 315      continue
          return
      endif
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
      dipin = 1.0d0
      call dipvar(dipin)
      call hafspc(vp,vs,dens)
      call moment0vol(Mvol,dvol,vp,vs,dens)
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
          call xderiv0vol(fx,x,y,c,d)
          call yderiv0vol(fy,x,y,c,d)
          call zderiv0vol(fz,x,y,c,d)
          call disp0stnvol(f,x,y,c,d)

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
          call xderiv0vol(fx,x,y,c,d)
          call yderiv0vol(fy,x,y,c,d)
          call zderiv0vol(fz,x,y,c,d)

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
      uxx = Mvol*(ux(4,1)+ux(4,2)+z*ux(4,3))
      uyx = Mvol*(ux(5,1)+ux(5,2)+z*ux(5,3))
      uzx = Mvol*(ux(6,1)+ux(6,2)+z*ux(6,3))
      uxy = Mvol*(uy(4,1)+uy(4,2)+z*uy(4,3))
      uyy = Mvol*(uy(5,1)+uy(5,2)+z*uy(5,3))
      uzy = Mvol*(uy(6,1)+uy(6,2)+z*uy(6,3))
      uxz = Mvol*(uz(4,1)+uz(4,2)+u(1)+z*uz(4,3))
      uyz = Mvol*(uz(5,1)+uz(5,2)+u(2)+z*uz(5,3))
      uzz = Mvol*(uz(6,1)+uz(6,2)+u(3)+z*uz(6,3))
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

      SUBROUTINE o92rectvol(ux,uy,uz,x,y,stdp,evdp,dipin,rakin,wid,len,
     1                   slip,vp,vs,dens)
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
      REAL*8 dipin,rakin,wid,len,slip,Mvol
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
      if (stdp.lt.0.0d0) then
          ux = 0.0d0
          uy = 0.0d0
          uz = 0.0d0
          return
      endif
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
      rakin = 0.0d0
      call dipvar(dipin)
      call hafspc(vp,vs,dens)
      call moment1vol(Mvol,slip)
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
      fac = one
      do 308 ii = 1,2
      do 307 jj = 1,2
          if (ii.eq.1.and.jj.eq.1) fac =  one
          if (ii.eq.1.and.jj.eq.2) fac = -one
          if (ii.eq.2.and.jj.eq.1) fac = -one
          if (ii.eq.2.and.jj.eq.2) fac =  one

          call rectvars(ksi(ii),eta(jj),q,z,ek(jj),ee(ii),eps)

          if (thru.eq.0) then
              call disp1vol(f,ksi(ii),eta(jj),q,z)
              do 305 i = 1,6
                  do 304 j = 1,3
                      u(i,j) = u(i,j) + fac*f(i,j)
  304             continue
  305         continue
          else
              call disp1vol(f,ksi(ii),eta(jj),q,z)
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
      ux = Mvol*((u(1,1)+u(1,2)+z*u(1,3)))
      uy = Mvol*((u(2,1)+u(2,2)+z*u(2,3))*cd
     1                                  - (u(3,1)+u(3,2)+z*u(3,3))*sd)
      uz = Mvol*((u(2,1)+u(2,2)-z*u(2,3))*sd
     1                                  + (u(3,1)+u(3,2)-z*u(3,3))*cd)
  399 continue

      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE o92rectstnvol(strain,x,y,stdp,evdp,dipin,wid,len,
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
      REAL*8 dipin,wid,len,slip,Mvol
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
      if (stdp.lt.0.0d0) then
          do 415 i = 1,3
              do 416 j = 1,3
                  strain(i,j) = 0.0d0
 416          continue
 415      continue
          return
      endif
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
      call dipvar(dipin)
      call hafspc(vp,vs,dens)
      call moment1vol(Mvol,slip)
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
      fac = one
      do 408 ii = 1,2
      do 407 jj = 1,2
          if (ii.eq.1.and.jj.eq.1) fac =  one
          if (ii.eq.1.and.jj.eq.2) fac = -one
          if (ii.eq.2.and.jj.eq.1) fac = -one
          if (ii.eq.2.and.jj.eq.2) fac =  one

          call rectvars(ksi(ii),eta(jj),q,z,ek(jj),ee(ii),eps)

          if (thru.eq.0) then
              call xderiv1vol(fx,ksi(ii),eta(jj),q,z)
              call yderiv1vol(fy,ksi(ii),eta(jj),q,z)
              call zderiv1vol(fz,ksi(ii),eta(jj),q,z)
              call disp1stnvol(f ,ksi(ii),eta(jj),q,z)

              do 405 i = 1,6
                  u(i)  =  u(i) + fac*f(i)
                  do 404 j = 1,3
                      fj(i,j) = fj(i,j) + fac*fx(i,j)
                      fk(i,j) = fk(i,j) + fac*fy(i,j)
                      fl(i,j) = fl(i,j) + fac*fz(i,j)
  404             continue
  405         continue
          else
              call xderiv1vol(fx,ksi(ii),eta(jj),q,z)
              call yderiv1vol(fy,ksi(ii),eta(jj),q,z)
              call zderiv1vol(fz,ksi(ii),eta(jj),q,z)
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
      uxx = Mvol*(fj(1,1)+fj(1,2)+z*fj(1,3))
      uyx = Mvol*((fj(2,1)+fj(2,2)+z*fj(2,3))*cd
     1                               - (fj(3,1)+fj(3,2)+z*fj(3,3))*sd)
      uzx = Mvol*((fj(2,1)+fj(2,2)-z*fj(2,3))*sd
     1                               + (fj(3,1)+fj(3,2)-z*fj(3,3))*cd)

      uxy = Mvol*(fk(1,1)+fk(1,2)+z*fk(1,3))
      uyy = Mvol*((fk(2,1)+fk(2,2)+z*fk(2,3))*cd
     1                               - (fk(3,1)+fk(3,2)+z*fk(3,3))*sd)
      uzy = Mvol*((fk(2,1)+fk(2,2)-z*fk(2,3))*sd
     1                               + (fk(3,1)+fk(3,2)-z*fk(3,3))*cd)

      uxz = Mvol*(fl(1,1)+fl(1,2)+u(1)+z*fl(1,3))
      uyz = Mvol*((fl(2,1)+fl(2,2)+u(2)+z*fl(2,3))*cd
     1                          - (fl(3,1)+fl(3,2)+u(3)+z*fl(3,3))*sd)
      uzz = Mvol*((fl(2,1)+fl(2,2)-u(2)-z*fl(2,3))*sd
     1                          + (fl(3,1)+fl(3,2)-u(3)-z*fl(3,3))*cd)
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

C----------------------------------------------------------------------C
C------------------- Source Moment Subroutines ------------------------C
C----------------------------------------------------------------------C

      SUBROUTINE moment0vol(Mvol,dvol,vp,vs,dens)
C----
C Compute point source isotropic volume change moment
C----
      IMPLICIT none
      REAL*8 Mvol,dvol,vp,vs,dens
      REAL*8 mu,lamda
      mu = dens*vs*vs
      lamda = dens*vp*vp - 2.0d0*mu
      Mvol = (lamda+2.0d0*mu)*dvol/mu
      RETURN
      END

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -C

      SUBROUTINE moment1vol(Mvol,slip)
C----
C Compute tensile moment for finite fault
C----
      IMPLICIT none
      REAL*8 pi
      PARAMETER (pi=4.0d0*datan(1.0d0))
      REAL*8 Mvol,slip
      Mvol = slip/(2.0d0*pi)
      RETURN
      END

C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C------------------- Point Source Subroutines -------------------------C
C----------------------------------------------------------------------C
C----------------------------------------------------------------------C

      SUBROUTINE disp0vol(f,x,y,c,d)
C----
C Calculate the components of displacement (from table on p. 1025 in
C Okada, 1992) in the array f(6,3).
C----
      IMPLICIT none
      REAL*8 f(6,3),x,y,c,d,z

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
      z = c - d
      if (thru.eq.0) then
          f(1,1) =  CA1*x/R3 - CA2*3.0d0*x*q*q/R5
          f(2,1) =  CA1*t/R3 - CA2*3.0d0*y*q*q/R5
          f(3,1) =  CA1*s/R3 - CA2*3.0d0*d*q*q/R5
          f(4,1) = -CA1*x/R3
          f(5,1) = -CA1*y/R3
          f(6,1) = -CA1*d/R3

          f(1,2) = 3.0d0*x*q*q/R5 - CB*I3*sdsd
          f(2,2) = 3.0d0*y*q*q/R5 - CB*I1*sdsd
          f(3,2) = 3.0d0*c*q*q/R5 - CB*I5*sdsd
          f(4,2) =                 CB*x/R3
          f(5,2) =                 CB*y/R3
          f(6,2) =                 CB*d/R3

          f(1,3) = -CC*3.0d0*x*s/R5 +a*15.0d0*c*x*q*q/R7 -a*3.0d0*x*z/R5
          f(2,3) =  CC*(s2d-3.0d0*y*s/R2)/R3
     1             + a*3.0d0*c*(t-y+5.0d0*y*q*q/R2)/R5 - a*3.0d0*y*z/R5
          f(3,3) = -CC*(1.0d0-A3*sdsd)/R3
     1             - a*3.0d0*c*(s-d+5.0d0*d*q*q/R2)/R5 + a*3.0d0*d*z/R5
          f(4,3) = CC*3.0d0*x*d/R5
          f(5,3) = CC*3.0d0*y*d/R5
          f(6,3) = CC*C3/R3
      else
          f(1,1) =  CA1*x/R3 - CA2*3.0d0*x*q*q/R5
          f(2,1) =  CA1*t/R3 - CA2*3.0d0*y*q*q/R5
          f(3,1) =  CA1*s/R3 - CA2*3.0d0*d*q*q/R5
          f(4,1) = -CA1*x/R3
          f(5,1) = -CA1*y/R3
          f(6,1) = -CA1*d/R3
      endif

      RETURN
      END

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -C

      SUBROUTINE xderiv0vol(fx,x,y,c,d)
C----
C Calculate the components of displacement x-derivatives (from table on
C p. 1026 in Okada, 1992) in the array fx(6,3).
C----
      IMPLICIT none
      REAL*8 fx(6,3),x,y,c,d,z

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
      z = c - d
      if (thru.eq.0) then
          fx(1,1) =  CA1*A3/R3        - CA2*3.0d0*q*q*A5/R5
          fx(2,1) = -CA1*3.0d0*x*t/R5 + CA2*15.0d0*xy*q*q/R7
          fx(3,1) = -CA1*3.0d0*x*s/R5 + CA2*15.0d0*xd*q*q/R7
          fx(4,1) = -CA1*A3/R3
          fx(5,1) =  CA1*3.0d0*xy/R5
          fx(6,1) =  CA1*3.0d0*xd/R5

          fx(1,2) =   3.0d0*q*q*A5/R5       - CB*J3*sdsd
          fx(2,2) = -15.0d0*xy*q*q/R7       - CB*J1*sdsd
          fx(3,2) = -15.0d0*xc*q*q/R7       - CB*K3*sdsd
          fx(4,2) =                           CB*A3/R3
          fx(5,2) =                          -CB*3.0d0*xy/R5
          fx(6,2) =                          -CB*3.0d0*xd/R5

          fx(1,3) = -CC*3.0d0*s*A5/R5
     1                          + a*15.0d0*c*q*q*A7/R7 - a*3.0d0*z*A5/R5
          fx(2,3) = -CC*3.0d0*x*(s2d-5.0d0*y*s/R2)/R5
     1        - a*15.0d0*xc*(t-y+7.0d0*y*q*q/R2)/R7 + a*15.0d0*xy*z/R7
          fx(3,3) =  CC*3.0d0*x*(1.0d0-(2.0d0+A5)*sdsd)/R5
     1        + a*15.0d0*xc*(s-d+7.0d0*d*q*q/R2)/R7 - a*15.0d0*xd*z/R7
          fx(4,3) =  CC*3.0d0*d*A5/R5
          fx(5,3) = -CC*15.0d0*xy*d/R7
          fx(6,3) = -CC*3.0d0*x*C5/R5
      else
          fx(1,1) =  CA1*A3/R3        - CA2*3.0d0*q*q*A5/R5
          fx(2,1) = -CA1*3.0d0*x*t/R5 + CA2*15.0d0*xy*q*q/R7
          fx(3,1) = -CA1*3.0d0*x*s/R5 + CA2*15.0d0*xd*q*q/R7
          fx(4,1) = -CA1*A3/R3
          fx(5,1) =  CA1*3.0d0*xy/R5
          fx(6,1) =  CA1*3.0d0*xd/R5
      endif

      RETURN
      END

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -C

      SUBROUTINE yderiv0vol(fy,x,y,c,d)
C----
C Calculate the components of displacement y-derivatives (from table on
C p. 1026 in Okada, 1992) in the array fx(6,3).
C----
      IMPLICIT none
      REAL*8 fy(6,3),x,y,c,d,z

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

      z = c - d
      if (thru.eq.0) then
          fy(1,1) = -CA1*3.0d0*xy/R5 - CA2*3.0d0*xq*W2/R5
          fy(2,1) =  CA1*(c2d-3.0d0*y*t/R2)/R3
     1                                  - CA2*3.0d0*(y*q*W2+q*q)/R5
          fy(3,1) =  CA1*(s2d-3.0d0*y*s/R2)/R3 - CA2*3.0d0*dq*W2/R5
          fy(4,1) =  CA1*3.0d0*xy/R5
          fy(5,1) = -CA1*B3/R3
          fy(6,1) =  CA1*3.0d0*y*d/R5

          fy(1,2) = 3.0d0*xq*W2/R5                - CB*J1*sdsd
          fy(2,2) = 3.0d0*yq*W2/R5 + 3.0d0*q*q/R5 - CB*J2*sdsd
          fy(3,2) = 3.0d0*c*q*W2/R5               - CB*K1*sdsd
          fy(4,2) =                                -CB*3.0d0*xy/R5
          fy(5,2) =                                 CB*B3/R3
          fy(6,2) =                                -CB*3.0d0*y*d/R5

          fy(1,3) = -CC*3.0d0*x*(s2d-5.0d0*y*s/R2)/R5
     1          - a*15.0d0*xc*(t-y+7.0d0*y*q*q/R2)/R7 + a*15.0d0*xy*z/R7
          fy(2,3) = -CC*3.0d0*(2.0d0*y*s2d+s*B5)/R5
     1          - a*3.0d0*c*(2.0d0*sdsd+10.0d0*y*(t-y)/R2
     2             -5.0d0*q*q*B7/R2)/R5 - a*3.0d0*z*B5/R5
          fy(3,3) =  CC*3.0d0*y*(1.0d0-A5*sdsd)/R5
     1          + a*3.0d0*c*((3.0d0+A5)*s2d-5.0d0*y*d*(2.0d0
     2             -7.0d0*q*q/R2)/R2)/R5 - a*15.0d0*y*d*z/R7
          fy(4,3) = -CC*15.0d0*xy*d/R7
          fy(5,3) =  CC*3.0d0*d*B5/R5
          fy(6,3) = -CC*3.0d0*y*C5/R5
      else
          fy(1,1) = -CA1*3.0d0*xy/R5 - CA2*3.0d0*xq*W2/R5
          fy(2,1) =  CA1*(c2d-3.0d0*y*t/R2)/R3
     1                                  - CA2*3.0d0*(y*q*W2+q*q)/R5
          fy(3,1) =  CA1*(s2d-3.0d0*y*s/R2)/R3 - CA2*3.0d0*dq*W2/R5
          fy(4,1) =  CA1*3.0d0*xy/R5
          fy(5,1) = -CA1*B3/R3
          fy(6,1) =  CA1*3.0d0*y*d/R5
      endif

      RETURN
      END

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -C

      SUBROUTINE zderiv0vol(fz,x,y,c,d)
C----
C Calculate the components of displacement z-derivatives (from table on
C p. 1026 in Okada, 1992) in the array fx(6,3).
C----
      IMPLICIT none
      REAL*8 fz(6,3),x,y,c,d,z

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

      z = c - d
      if (thru.eq.0) then
          fz(1,1) =  CA1*3.0d0*xd/R5 - CA2*3.0d0*xq*W3/R5
          fz(2,1) = -CA1*(s2d-3.0d0*d*t/R2)/R3 - CA2*3.0d0*yq*W3/R5
          fz(3,1) =  CA1*(c2d+3.0d0*d*s/R2)/R3
     1                                - CA2*3.0d0*q*(d*W3-q)/R5
          fz(4,1) = -CA1*3.0d0*xd/R5
          fz(5,1) = -CA1*3.0d0*y*d/R5
          fz(6,1) =  CA1*C3/R3

          fz(1,2) = 3.0d0*xq*W3/R5 + CB*K3*sdsd
          fz(2,2) = 3.0d0*yq*W3/R5 + CB*K1*sdsd
          fz(3,2) = 3.0d0*c*q*W3/R5 - CB*A3*sdsd/R3
          fz(4,2) =  CB*3.0d0*xd/R5
          fz(5,2) =  CB*3.0d0*y*d/R5
          fz(6,2) = -CB*C3/R3

          fz(1,3) = -CC*3.0d0*x*(c2d+5.0d0*d*s/R2)/R5
     1       + a*15.0d0*xc*(s-d+7.0d0*d*q*q/R2)/R7
     2       - a*3.0d0*x*(1.0d0+5.0d0*d*z/R2)/R5
          fz(2,3) =  CC*3.0d0*(d*B5*s2d-y*C5*c2d)/R5
     1       + a*3.0d0*c*((3.0d0+A5)*s2d-5.0d0*y*d*(2.0d0
     2       -7.0d0*q*q/R2)/R2)/R5 - a*3.0d0*y*(1.0d0+5.0d0*d*z/R2)/R5
          fz(3,3) = -CC*3.0d0*d*(1.0d0-A5*sdsd)/R5
     1       - a*3.0d0*c*(c2d+10.0d0*d*(s-d)/R2-5.0d0*q*q*C7/R2)/R5
     2       - a*3.0d0*z*(1.0d0+C5)/R5
          fz(4,3) = -CC*3.0d0*x*C5/R5
          fz(5,3) = -CC*3.0d0*y*C5/R5
          fz(6,3) =  CC*3.0d0*d*(2.0d0+C5)/R5
      else
          fz(1,1) =  CA1*3.0d0*xd/R5 - CA2*3.0d0*xq*W3/R5
          fz(2,1) = -CA1*(s2d-3.0d0*d*t/R2)/R3 - CA2*3.0d0*yq*W3/R5
          fz(3,1) =  CA1*(c2d+3.0d0*d*s/R2)/R3
     1                                - CA2*3.0d0*q*(d*W3-q)/R5
          fz(4,1) = -CA1*3.0d0*xd/R5
          fz(5,1) = -CA1*3.0d0*y*d/R5
          fz(6,1) =  CA1*C3/R3
      endif

      RETURN
      END

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -C

      SUBROUTINE disp0stnvol(f,x,y,c,d)
      IMPLICIT none
      REAL*8 f(6),x,y,c,d,z

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

      z = c - d
      f(1) = -CC*3.0d0*x*s/R5 +a*15.0d0*c*x*q*q/R7 -a*3.0d0*x*z/R5
      f(2) =  CC*(s2d-3.0d0*y*s/R2)/R3
     1             + a*3.0d0*c*(t-y+5.0d0*y*q*q/R2)/R5 - a*3.0d0*y*z/R5
      f(3) = -CC*(1.0d0-A3*sdsd)/R3
     1             - a*3.0d0*c*(s-d+5.0d0*d*q*q/R2)/R5 + a*3.0d0*d*z/R5
      f(4) = CC*3.0d0*x*d/R5
      f(5) = CC*3.0d0*y*d/R5
      f(6) = CC*C3/R3

      RETURN
      END

C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C------------------- Finite Source Subroutines ------------------------C
C----------------------------------------------------------------------C
C----------------------------------------------------------------------C

      SUBROUTINE disp1vol(f,ksi,eta,q,z)
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
          f(1,1) = -CA1*logRe - CA2*q*q*Y11
          f(2,1) = -CA1*logRk - CA2*q*q*X11
          f(3,1) =   TH*0.5d0 - CA2*q*(eta*X11+ksi*Y11)
          f(4,1) = 0.0d0
          f(5,1) = 0.0d0
          f(6,1) = 0.0d0

          f(1,2) =  q*q*Y11                  - CB*I3*sdsd
          f(2,2) =  q*q*X11                  + CB*ksi*sdsd/Rd
          f(3,2) =  q*(eta*X11+ksi*Y11) - TH - CB*I4*sdsd
          f(4,2) = 0.0d0
          f(5,2) = 0.0d0
          f(6,2) = 0.0d0

          f(1,3) = -CC*(sd/R+q*Y11*cd) - a*(z*Y11-q*q*Z32)
          f(2,3) =  CC*2.0d0*ksi*Y11*sd + dbar*X11
     1                                   - a*cbar*(X11-q*q*X32)
          f(3,3) =  CC*(ybar*X11+ksi*Y11*cd)
     1                                  + a*q*(cbar*eta*X32+ksi*Z32)
          f(4,3) = 0.0d0
          f(5,3) = 0.0d0
          f(6,3) = 0.0d0
      else
          f(1,1) = -CA1*logRe - CA2*q*q*Y11
          f(2,1) = -CA1*logRk - CA2*q*q*X11
          f(3,1) =   TH*0.5d0 - CA2*q*(eta*X11+ksi*Y11)
          f(4,1) = 0.0d0
          f(5,1) = 0.0d0
          f(6,1) = 0.0d0
      endif
      RETURN
      END

C-------------------

      SUBROUTINE xderiv1vol(fx,ksi,eta,q,z)
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
          fx(1,1) = -CA1*ksi*Y11 + CA2*ksi*q*q*Y32
          fx(2,1) = -CA1/R       + CA2*q*q/R3
          fx(3,1) = -CA1*q*Y11   - CA2*q*q*q*Y32
          fx(4,1) = 0.0d0
          fx(5,1) = 0.0d0
          fx(6,1) = 0.0d0

          fx(1,2) = -ksi*q*q*Y32      - CB*J4*sdsd
          fx(2,2) = -q*q/R3           - CB*J5*sdsd
          fx(3,2) = q*q*q*Y32         - CB*J6*sdsd
          fx(4,2) = 0.0d0
          fx(5,2) = 0.0d0
          fx(6,2) = 0.0d0

          fx(1,3) =  CC*ksi*sd/R3 + ksi*q*Y32*cd
     1                + a*ksi*(3.0d0*cbar*eta/R5 - 2.0d0*Z32 - Z0)
          fx(2,3) =  CC*2.0d0*Y0*sd - dbar/R3
     1                + a*cbar*(1.0d0-3.0d0*q*q/R2)/R3
          fx(3,3) = -CC*(ybar/R3-Y0*cd) - a*(3.0d0*cbar*eta*q/R5-q*Z0)
          fx(4,3) = 0.0d0
          fx(5,3) = 0.0d0
          fx(6,3) = 0.0d0
      else
          fx(1,1) = -CA1*ksi*Y11 + CA2*ksi*q*q*Y32
          fx(2,1) = -CA1/R       + CA2*q*q/R3
          fx(3,1) = -CA1*q*Y11   - CA2*q*q*q*Y32
          fx(4,1) = 0.0d0
          fx(5,1) = 0.0d0
          fx(6,1) = 0.0d0
      endif

      RETURN
      END

C-------------------

      SUBROUTINE yderiv1vol(fy,ksi,eta,q,z)
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
          fy(1,1) = -CA1*(cd/R+q*Y11*sd)             - CA2*q*F2
          fy(2,1) = -CA1*ybar*X11                    - CA2*q*G2
          fy(3,1) =  CA1*(dbar*X11+ksi*Y11*sd)       + CA2*q*H2
          fy(4,1) = 0.0d0
          fy(5,1) = 0.0d0 
          fy(6,1) = 0.0d0

          fy(1,2) =  q*F2 - CB*J1*sdsd
          fy(2,2) =  q*G2 - CB*J2*sdsd
          fy(3,2) = -q*H2 - CB*J3*sdsd
          fy(4,2) = 0.0d0
          fy(5,2) = 0.0d0
          fy(6,2) = 0.0d0

          fy(1,3) =  CC*(q/R3+Y0*cdsd)
     1                  + a*(z*cd/R3+3.0d0*cbar*dbar*q/R5-q*Z0*sd)
          fy(2,3) = -CC*2.0d0*ksi*P2*sd - ybar*dbar*X32
     1                   + a*cbar*((ybar+2.0d0*q*sd)*X32-ybar*q*q*X53)
          fy(3,3) = -CC*(ksi*P2*cd-X11+ybar*ybar*X32)
     1                 + a*cbar*((dbar+2.0d0*q*cd)*X32-ybar*eta*q*X53)
     2                 + a*ksi*Q2
          fy(4,3) = 0.0d0
          fy(5,3) = 0.0d0
          fy(6,3) = 0.0d0
      else
          fy(1,1) = -CA1*(cd/R+q*Y11*sd)             - CA2*q*F2
          fy(2,1) = -CA1*ybar*X11                    - CA2*q*G2
          fy(3,1) =  CA1*(dbar*X11+ksi*Y11*sd)       + CA2*q*H2
          fy(4,1) = 0.0d0
          fy(5,1) = 0.0d0 
          fy(6,1) = 0.0d0
      endif
      z=z

      RETURN
      END

C-------------------

      SUBROUTINE zderiv1vol(fz,ksi,eta,q,z)
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
          fz(1,1) =  CA1*(sd/R-q*Y11*cd)             - CA2*q*F3
          fz(2,1) =  CA1*dbar*X11                    - CA2*q*G3
          fz(3,1) =  CA1*(ybar*X11+ksi*Y11*cd)       - CA2*q*H3
          fz(4,1) = 0.0d0
          fz(5,1) = 0.0d0
          fz(6,1) = 0.0d0

          fz(1,2) =  q*F3 + CB*K3*sdsd
          fz(2,2) =  q*G3 + CB*ksi*D11*sdsd
          fz(3,2) = -q*H3 + CB*K4*sdsd
          fz(4,2) = 0.0d0
          fz(5,2) = 0.0d0
          fz(6,2) = 0.0d0

          fz(1,3) = -eta/R3 + Y0*cdcd
     1               - a*(z*sd/R3-3.0d0*cbar*ybar*q/R5-Y0*sdsd+q*Z0*cd)
          fz(2,3) =  CC*2.0d0*ksi*P3*sd - X11 + dbar*dbar*X32
     1                - a*cbar*((dbar-2.0d0*q*cd)*X32-dbar*q*q*X53)
          fz(3,3) =  CC*(ksi*P3*cd+ybar*dbar*X32)
     1                 + a*cbar*((ybar-2.0d0*q*sd)*X32+dbar*eta*q*X53)
     2                 + a*ksi*Q3
          fz(4,3) = 0.0d0
          fz(5,3) = 0.0d0
          fz(6,3) = 0.0d0
      else
          fz(1,1) =  CA1*(sd/R-q*Y11*cd)             - CA2*q*F3
          fz(2,1) =  CA1*dbar*X11                    - CA2*q*G3
          fz(3,1) =  CA1*(ybar*X11+ksi*Y11*cd)       - CA2*q*H3
          fz(4,1) = 0.0d0
          fz(5,1) = 0.0d0
          fz(6,1) = 0.0d0
      endif
      z=z

      RETURN
      END

C-------------------

      SUBROUTINE disp1stnvol(f,ksi,eta,q,z)
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

      f(1) = -CC*(sd/R+q*Y11*cd) - a*(z*Y11-q*q*Z32)
      f(2) =  CC*2.0d0*ksi*Y11*sd + dbar*X11
     1                                   - a*cbar*(X11-q*q*X32)
      f(3) =  CC*(ybar*X11+ksi*Y11*cd)
     1                                  + a*q*(cbar*eta*X32+ksi*Z32)
      f(4) = 0.0d0
      f(5) = 0.0d0
      f(6) = 0.0d0

      RETURN
      END

