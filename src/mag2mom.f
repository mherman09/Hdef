      PROGRAM main
      IMPLICIT none
      REAL*8 mag,mom,total
      total = 0.0d0
  101 read (*,*,end=102) mag
          mom = 10**(1.5d0*(mag+10.7d0))/1d7
          total = total + mom
          write(*,1001) mom
          goto 101
  102 continue
      write (*,1002) total
 1001 format('Moment: ',G14.6,' Nm')
 1002 format('Total Moment: ',G14.6,' Nm')
      END
