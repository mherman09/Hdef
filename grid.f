      PROGRAM grid
C----
C Create a horizontal 2-dimensional grid
C----
      IMPLICIT NONE
      REAL*8 x1,x2,y1,y2,dx,dy
      REAL*8 x,y,z
      INTEGER i,j,nx,ny,chk,p

C----
C Get parameters from the command line
C----
      call gcmdln(x1,x2,y1,y2,nx,ny,z,p)
      if (p.eq.0) then
          open (unit=101,file='grid.out',status='unknown')
          rewind 101
      endif

C----
C Calculate grid increments
C----
      if (nx.le.0) then
          print *,'Number of x grid points must be 1 or greater'
          call usage()
      elseif (nx.eq.1) then
          dx = 0.0d0
      else
          dx = (x2-x1)/dble(nx-1)
      endif

      if (ny.le.0) then
          print *,'Number of y grid points must be 1 or greater'
          call usage()
      elseif (ny.eq.1) then
          dy = 0.0d0
      else
          dy = (y2-y1)/dble(ny-1)
      endif

C----
C Generate grid
C----
      do 16 i = 0,nx-1
          x = x1 + dble(i)*dx
          do 15 j = 0,ny-1
              y = y1 + dble(j)*dy
              if (p.eq.0) then
                  write(101,9999) x,y,z
              else
                  write(*,9999) x,y,z
              endif
   15     continue
   16 continue

 9999 format (3F16.8)

      END 

C======================================================================C

      SUBROUTINE gcmdln(x1,x2,y1,y2,nx,ny,z,p)
      IMPLICIT none
      CHARACTER*20 tag
      REAL*8 x1,x2,y1,y2,z
      INTEGER narg,i,nx,ny,p

      x1 = 0.0d0
      x2 = 1.0d0
      y1 = 0.0d0
      y2 = 1.0d0
      nx = 1
      ny = 1
      z = 0.0d0
      p = 0
      
      narg = iargc()
      if (narg.eq.0) call usage()

      i = 0
   11 i = i + 1
      if (i.gt.narg) goto 12
      call getarg(i,tag)
      if (tag(1:2).eq.'-x') then
          i = i + 1
          call getarg(i,tag)
          read (tag,'(BN,F8.0)') x1
          i = i + 1
          call getarg(i,tag)
          read (tag,'(BN,F8.0)') x2
      elseif (tag(1:2).eq.'-y') then
          i = i + 1
          call getarg(i,tag)
          read (tag,'(BN,F8.0)') y1
          i = i + 1
          call getarg(i,tag)
          read (tag,'(BN,F8.0)') y2
      elseif (tag(1:2).eq.'-z') then
          i = i + 1
          call getarg(i,tag)
          read (tag,'(BN,F8.0)') z
      elseif (tag(1:3).eq.'-nx') then
          i = i + 1
          call getarg(i,tag)
          read (tag,'(BN,I6)') nx
      elseif (tag(1:3).eq.'-ny') then
          i = i + 1
          call getarg(i,tag)
          read (tag,'(BN,I6)') ny
      elseif (tag(1:2).eq.'-h'.or.tag(1:2).eq.'-?') then
          call usage()
      elseif (tag(1:2).eq.'-p') then
          p = 1
      endif
      goto 11

   12 continue

      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE usage()
      IMPLICIT none
      write(*,*)
     1 'Usage: grid -x [X1 X2] -y [Y1 Y2] -z [Z] -nx [NX] -ny [NY] ',
     2 '-p -h/-?'
      write(*,*)
     1 '  -x [X1 X2] (Default X1=0 X2=1) define x limits'
      write(*,*)
     1 '  -y [Y1 Y2] (Default Y1=0 Y2=1) define y limits'
      write(*,*)
     1 '  -z [Z]     (Default Z=0) define z'
      write(*,*)
     1 '  -nx [NX]   (Default NX=1) define number of x grid points'
      write(*,*)
     1 '  -ny [NY]   (Default NY=1) define number of y grid points'
      write(*,*)
     1 '  -p         (Default off) print output to stdout (will not ',
     2                            'write file grid.out)'
      write (*,*)
     1 '  -h/-?        help'
      write (*,*)
     1 '  grid creates an evenly spaced file grid.out'
      write (*,*)
     1 '        X1 Y1 Z'
      write (*,*)
     1 '         :  :'
      write (*,*)
     1 '        X1 Y2 Z'
      write (*,*)
     1 '         :  :'
      write (*,*)
     1 '        X2 Y2 Z'
      write (*,*) ''

      STOP
      END
