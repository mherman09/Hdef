      PROGRAM main
      IMPLICIT none
      REAL*8 pi,d2r,pi2,const
      PARAMETER (pi=4.0d0*datan(1.0d0),d2r=pi/1.8d2,pi2=pi*0.5d0,
     1        const=1.0d0/dsqrt(2.0d0))
      REAL*8 s,d,r
      REAL*8 normal(3),slip(3),tmp(3),alpha
      REAL*8 p(3),n(3),t(3),fth,fss,fno
      INTEGER i
  101 read (*,*,END=102) s,d,r
          s = s*d2r
          d = d*d2r
          r = r*d2r
          normal(1) =  0.0d0
          normal(2) = -dsin(d)
          normal(3) =  dcos(d)
          alpha = pi2 - s
          tmp(1) = dcos(alpha)*normal(1) - dsin(alpha)*normal(2)
          tmp(2) = dsin(alpha)*normal(1) + dcos(alpha)*normal(2)
          normal(1) = tmp(1)
          normal(2) = tmp(2)
          slip(1) = dcos(r)
          slip(2) = dsin(r)*dcos(d)
          slip(3) = dsin(r)*dsin(d)
          tmp(1) = dcos(alpha)*slip(1) - dsin(alpha)*slip(2)
          tmp(2) = dsin(alpha)*slip(1) + dcos(alpha)*slip(2)
          slip(1) = tmp(1)
          slip(2) = tmp(2)
          do 103 i = 1,3
              t(i) = const*(normal(i)+slip(i))
              p(i) = const*(normal(i)-slip(i))
  103     continue
          n(1) = normal(2)*slip(3) - normal(3)*slip(2)
          n(2) = normal(3)*slip(1) - normal(1)*slip(3)
          n(3) = normal(1)*slip(2) - normal(2)*slip(1)
          !print *,'p',p(1),p(2),p(3)
          !print *,'n',n(1),n(2),n(3)
          !print *,'t',t(1),t(2),t(3)
          fno = datan(p(3)/dsqrt(p(1)*p(1)+p(2)*p(2)))
          fss = datan(n(3)/dsqrt(n(1)*n(1)+n(2)*n(2)))
          fth = datan(t(3)/dsqrt(t(1)*t(1)+t(2)*t(2)))
          fno = dsin(fno)*dsin(fno)
          fss = dsin(fss)*dsin(fss)
          fth = dsin(fth)*dsin(fth)
          write(*,'(3F10.4)') fno,fss,fth
          goto 101
  102 continue
      END
