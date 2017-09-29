      PROGRAM main

      IMPLICIT none
      REAL*8 mag,mom
  101 read (*,*,end=102) mom
          mag = 2.0d0/3.0d0*dlog10(mom) - 10.7
          write(*,1001) mag
          goto 101
  102 continue
 1001 format('Magnitude: ',F5.2)
      END
