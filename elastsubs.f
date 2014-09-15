
      SUBROUTINE strain2stress(stress,strain,vp,vs,dens)
C----
C CALCULATE STRESS FROM STRAIN ASSUMING ELASTIC, ISOTROPIC MATERIAL.
C----
      IMPLICIT NONE
      REAL*8 stress(3,3),strain(3,3),vp,vs,dens,lam,mu,diag
      
      mu  = dens*vs*vs
      lam = dens*vp*vp - 2.0d0*mu
      if (mu.lt.10.0e9) mu = 10.0e9
      if (lam.lt.10.0e9) lam = 10.0e9

      diag = strain(1,1) + strain(2,2) + strain(3,3)

      stress(1,1) = lam*diag + 2.0d0*mu*strain(1,1)
      stress(2,2) = lam*diag + 2.0d0*mu*strain(2,2)
      stress(3,3) = lam*diag + 2.0d0*mu*strain(3,3)
      stress(1,2) = 2.0d0*mu*strain(1,2)
      stress(1,3) = 2.0d0*mu*strain(1,3)
      stress(2,3) = 2.0d0*mu*strain(2,3)
      stress(2,1) = stress(1,2)
      stress(3,1) = stress(1,3)
      stress(3,2) = stress(2,3)

      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE coulomb(coul,norml,shear,stress,strin,dipin,rakin,
     1                   coeffr)
C----
C GIVEN A STRESS FIELD AND THE ORIENTATION OF A TARGET PLANE, CALCULATE
C THE NORMAL, SHEAR, AND COULOMB STRESS ON THE PLANE.
C----
      IMPLICIT none
      REAL*8 pi,d2r
      PARAMETER (pi=4.0d0*datan(1.0d0),d2r=pi/180.0d0)
      REAL*8 coul,norml,shear,stress(3,3),coeffr,mag
      REAL*8 strin,dipin,rakin,str,dip,rak,n(3),trac(3),r(3),s(3)
      INTEGER i
C----
C NORMAL VECTOR OF TARGET PLANE
C----
      str = strin*d2r
      dip = dipin*d2r
      rak = rakin*d2r
      n(1) = dsin(dip)*dsin(str+pi/2.0d0)
      n(2) = dsin(dip)*dcos(str+pi/2.0d0)
      n(3) = dcos(dip)
C----
C TRACTION ON TARGET PLANE
C----
      do 11 i = 1,3
          trac(i) = stress(i,1)*n(1)+stress(i,2)*n(2)+stress(i,3)*n(3)
   11 continue
C----
C NORMAL STRESS ON TARGET PLANE
C (POSITIVE => DILATION => COULOMB STRESS INCREASES)
C----
      norml = 0.0d0
      do 12 i = 1,3
          norml = norml + trac(i)*n(i)
   12 continue
C----
C RAKE VECTOR OF TARGET PLANE
C----
      r(1) =  dcos(str-pi/2.0d0)*dcos(rak)
     1                          + dsin(str-pi/2.0d0)*dsin(rak)*dcos(dip)
      r(2) = -dsin(str-pi/2.0d0)*dcos(rak)
     1                          + dcos(str-pi/2.0d0)*dsin(rak)*dcos(dip)
      r(3) =                                         dsin(rak)*dsin(dip)
C----
C SHEAR STRESS (VECTOR) ON TARGET PLANE
C----
      s(1) = trac(1) - norml*n(1)
      s(2) = trac(2) - norml*n(2)
      s(3) = trac(3) - norml*n(3)
C----
C COMPONENT OF SHEAR STRESS IN DIRECTION OF RAKE ON TARGET PLANE
C----
      shear = s(1)*r(1) + s(2)*r(2) + s(3)*r(3)
C----
C COULOMB STRESS CHANGE
C----
      coul = shear + coeffr*norml

      RETURN
      END




