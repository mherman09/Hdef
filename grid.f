      PROGRAM grid
C----
C Create a horizontal 2-dimensional grid
C----
      IMPLICIT NONE
      CHARACTER*40 ofile
      REAL*8 x1,x2,y1,y2,dx,dy
      REAL*8 x,y,z
      INTEGER i,j,nx,ny,p

C----
C Get parameters from the command line
C----
      call gcmdln(x1,x2,nx,dx,y1,y2,ny,dy,z,ofile,p)
      !print *,'X1',x1,'X2',x2
      !print *,'Y1',y1,'Y2',y2
      !print *,'Z',z
      !print *,'NX',nx
      !print *,'NY',ny
      !print *,'OFILE ',ofile
      !print *,'P',p
      if (ofile.eq.'none'.and.p.eq.0) p = 1

C----
C Calculate grid increments
C----
      !print *,'DX',dx
      !print *,'DY',dy
      if (dx.gt.-1.0d98) then
          nx = int((x2-x1)/dx)+1
      elseif (nx.le.0) then
          call usage('!! Error: NX must be 1 or greater')
      elseif (nx.eq.1) then
          !write(0,*) '!! Warning: NX=1; using X1 as x value'
          dx = 0.0d0
      else
          dx = (x2-x1)/dble(nx-1)
      endif
      if (dy.gt.-1.0d98) then
          ny = int((y2-y1)/dy)+1
      elseif (ny.le.0) then
          call usage('!! Error: NY must be 1 or greater')
      elseif (ny.eq.1) then
          !write(0,*) '!! Warning: NY=1; using Y1 as y value'
          dy = 0.0d0
      else
          dy = (y2-y1)/dble(ny-1)
      endif
      !print *,'DX',dx
      !print *,'DY',dy

C----
C Generate grid
C----
      if (p.eq.0) open(unit=101,file=ofile,status='unknown')
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
      stop

 9999 format (3F16.8)

      END 

C======================================================================C

      SUBROUTINE gcmdln(x1,x2,nx,dx,y1,y2,ny,dy,z,ofile,p)
      IMPLICIT none
      CHARACTER*40 tag,ofile
      REAL*8 x1,x2,y1,y2,z,dx,dy
      INTEGER narg,i,nx,ny,p
      x1 = 0.0d0
      x2 = 1.0d0
      nx = 1
      y1 = 0.0d0
      y2 = 1.0d0
      ny = 1
      dx = -1.0d99
      dy = -1.0d99
      z = 0.0d0
      ofile = 'none'
      p = 0
      narg = iargc()
      if (narg.eq.0) then
          call usage('!! Error: no command line arguments specified')
      endif
      i = 0
   11 i = i + 1
      if (i.gt.narg) goto 12
      call getarg(i,tag)
      if (tag(1:2).eq.'-x') then
          i = i + 1
          call getarg(i,tag)
          read (tag,'(BN,F12.0)') x1
          i = i + 1
          call getarg(i,tag)
          read (tag,'(BN,F12.0)') x2
      elseif (tag(1:2).eq.'-y') then
          i = i + 1
          call getarg(i,tag)
          read (tag,'(BN,F12.0)') y1
          i = i + 1
          call getarg(i,tag)
          read (tag,'(BN,F12.0)') y2
      elseif (tag(1:2).eq.'-z') then
          i = i + 1
          call getarg(i,tag)
          read (tag,'(BN,F12.0)') z
      elseif (tag(1:3).eq.'-nx') then
          i = i + 1
          call getarg(i,tag)
          read (tag,'(BN,I7)') nx
      elseif (tag(1:3).eq.'-ny') then
          i = i + 1
          call getarg(i,tag)
          read (tag,'(BN,I7)') ny
      elseif (tag(1:3).eq.'-dx') then
          i = i + 1
          call getarg(i,tag)
          read (tag,'(BN,F12.0)') dx
      elseif (tag(1:3).eq.'-dy') then
          i = i + 1
          call getarg(i,tag)
          read (tag,'(BN,F12.0)') dy
      elseif (tag(1:2).eq.'-o') then
          i = i + 1
          call getarg(i,ofile)
      elseif (tag(1:2).eq.'-p') then
          p = 1
      elseif (tag(1:2).eq.'-h') then
          call usage(' ')
      endif
      goto 11

   12 continue

      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE usage(str)
      IMPLICIT none
      INTEGER lstr
      CHARACTER str*(*)
      if (str.ne.' ') then
          lstr = len(str)
          write(*,*) str(1:lstr)
          write(*,*)
      endif
      write(*,*)
     1 'Usage: grid -x X1 X2 [-nx NX|-dx DX] ',
     2            '[-y Y1 Y2 [-ny NY|-dy DY]] [-z Z] ',
     2             '[-o OFILE] [-p] [-h]'
      write(*,*)
      write(*,*)
     1 '-x X1 X2    First column (x) limits'
      write(*,*)
     1 '-y Y1 Y2    Second column (y) limits'
      write(*,*)
     1 '-z Z        Third column (z) value'
      write(*,*)
     1 '-nx NX      Number of x grid points'
      write(*,*)
     1 '-ny NY      Number of y grid points'
      write(*,*)
     1 '-dx DX      Increment in x direction (overrides -nx)'
      write(*,*)
     1 '-dy DY      Increment in y direction (overrides -ny)'
      write (*,*)
     1 '-o OFILE    Output to file (default prints to standard ',
     2              'output)'
      write (*,*)
     1 '-p          Print results to standard output (overrides -o)'
      write (*,*)
     1 '-h/-?       Online help (this screen)'
      write (*,*)
      STOP
      END
