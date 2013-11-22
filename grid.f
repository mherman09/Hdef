      PROGRAM grid
C----
C Create a horizontal 2-dimensional grid
C----
      IMPLICIT NONE
      REAL x1,x2,y1,y2,dx,dy
      REAL x,y,z
      INTEGER i,j,nx,ny,chk,p

C----
C Get parameters from the command line
C----
      call gcmdln(x1,x2,y1,y2,nx,ny,z,p)
      if (p.eq.0) then
          open (unit=101,file='grid.out',status='unknown')
          rewind 101
      endif

      if (nx.eq.0) then
          dx = 0.0
      else
          dx = (x2-x1)/real(nx)
      endif
      if (ny.eq.0) then
          dy = 0.0
      else
          dy = (y2-y1)/real(ny)
      endif

      do 16 i = 0,nx
          x = x1 + real(i)*dx
          do 15 j = 0,ny
              y = y1 + real(j)*dy
              if (p.eq.0) then
                  write(101,*) x,y,z
              else
                  write(*,*) x,y,z
              endif
   15     continue
   16 continue

      END 

C======================================================================C

      SUBROUTINE gcmdln(x1,x2,y1,y2,nx,ny,z,p)
      IMPLICIT none
      CHARACTER*20 tag
      REAL x1,x2,y1,y2,z
      INTEGER narg,i,nx,ny,p

      x1 = 0.0
      x2 = 1.0
      y1 = 0.0
      y2 = 1.0
      nx = 0
      ny = 0
      z = 0.0
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
     1 '  -nx [NX]   (Default NX=1) define number of x-subdivisions'
      write(*,*)
     1 '  -ny [NY]   (Default NY=1) define number of y-subdivisions'
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
