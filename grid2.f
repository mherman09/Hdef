      PROGRAM grid2
C----
C Create a dipping grid surface
C----
      IMPLICIT none
      REAL*8 pi,d2r,r2d
      PARAMETER (pi=4.0d0*datan(1.0d0),d2r=pi/1.8d2,r2d=1.8d2/pi)
      REAL*8 clon,clat,cdep,wid,len,str,dip
      REAL*8 dist,az,clon2,clat2,lon,lat,dep,dwid,dlen,xdist,ydist
      INTEGER nwid,nlen,iwid,ilen,p,user

      REAL*8 zro,haf,one
      DATA zro,haf,one/0.0d0,0.5d0,1.0d0/
C----
C Get parameters from the command line
C----
      call gcmdln(clon,clat,cdep,wid,len,nwid,nlen,str,dip,p,user)
      if (p.eq.0) then
          open (unit=101,file='grid2.ll',status='unknown')
          open (unit=102,file='grid2.xy',status='unknown')
          open (unit=103,file='grid2.clip',status='unknown')
          rewind 101
          rewind 102
          rewind 103
      endif
      dip = dip*d2r

C      if (user.eq.1) then
C          write(*,*) 'Enter starting and ending x-points X1 X2:'
C          read *,x1,x2
C          write(*,*) 'Enter number of grid points along x-dimension:'
C          read *,nx
C          write(*,*) 'Enter starting and ending y-points Y1 Y2:'
C          read *,y1,y2
C          write(*,*) 'Enter number of grid points along y-dimension:'
C          read *,ny
C          write(*,*) 'Enter z value for grid:'
C          read *,z
          if (p.eq.0) then
C              write(*,*) ''
C              write(*,*) 'Grid written to file grid2.out'
          endif
C      endif
      dwid = wid/(dble(nwid)-one)*dcos(dip)
      dlen = len/(dble(nlen)-one)

      do 12 iwid = 1,nwid
          dist = dwid*(dble(iwid)-haf*(dble(nwid)+one))
          ydist = dist
          az = str + 90.0d0
          dep = cdep+dist*dtan(dip)
          call dlola(clon2,clat2,clon,clat,dist,az)
          clon2 = clon2*r2d
          clat2 = clat2*r2d
          do 11 ilen = 1,nlen
              dist = dlen*(dble(ilen)-haf*(dble(nlen)+one))
              xdist = dist
              az = str
              call dlola(lon,lat,clon2,clat2,dist,az)
              lon = lon*r2d
              lat = lat*r2d
              write(101,*),lon,lat,dep
              write(102,*),xdist,ydist,dep
   11     continue
   12 continue

C----
C Write corners of grid to file grid2.clip for use with GMT
C program psclip
C----
      do 14 iwid = 1,2
          if (iwid.eq.1) then
              if (haf*wid*dsin(dip).le.cdep) then
                  dist = -haf*wid*dcos(dip)
              else
                  dist = -cdep/dtan(dip)
              endif
          endif
          if (iwid.eq.2) dist =  haf*wid*dcos(dip)
          az = str + 90.0d0
          call dlola(clon2,clat2,clon,clat,dist,az)
          clon2 = clon2*r2d
          clat2 = clat2*r2d
          az = str
          do 13 ilen = 1,2
              if (iwid.eq.1) then
                  if (ilen.eq.1) dist = -haf*len
                  if (ilen.eq.2) dist =  haf*len
              else
                  if (ilen.eq.1) dist =  haf*len
                  if (ilen.eq.2) dist = -haf*len
              endif
              call dlola(lon,lat,clon2,clat2,dist,az)
              write (103,*) lon*r2d,lat*r2d
 13       continue
 14   continue

      if (haf*wid*dsin(dip).le.cdep) then
          dist = -haf*wid*dcos(dip)
      else
          dist = -cdep/dtan(dip)
      endif
      az = str + 90.0d0
      call dlola(clon2,clat2,clon,clat,dist,az)
          clon2 = clon2*r2d
          clat2 = clat2*r2d
      dist = -haf*len
      az = str
      call dlola(lon,lat,clon2,clat2,dist,az)
      write (103,*) lon*r2d,lat*r2d

      END

C======================================================================C

      SUBROUTINE gcmdln(clon,clat,cdep,wid,len,nwid,nlen,str,dip,p,user)
      IMPLICIT none
      CHARACTER*30 tag
      REAL*8 clon,clat,cdep,wid,len,str,dip
      INTEGER nwid,nlen,narg,i,p,user

      clon = 0.0d0
      clat = 0.0d0
      cdep = 0.0d0
      wid = 1.0d0
      len = 1.0d0
      nwid = 5
      nlen = 5
      str = 0.0d0
      dip = 90.0d0
      p = 0
      user = 0
      
      narg = iargc()
      if (narg.eq.0) call usage()

      i = 0
   11 i = i + 1
      if (i.gt.narg) goto 12
      call getarg(i,tag)
      if (tag(1:2).eq.'-c') then
          i = i + 1
          call getarg(i,tag)
          read (tag,'(BN,F10.0)') clon
          i = i + 1
          call getarg(i,tag)
          read (tag,'(BN,F10.0)') clat
          i = i + 1
          call getarg(i,tag)
          read (tag,'(BN,F10.0)') cdep
      elseif (tag(1:2).eq.'-w') then
          i = i + 1
          call getarg(i,tag)
          read (tag,'(BN,F10.0)') wid
      elseif (tag(1:2).eq.'-l') then
          i = i + 1
          call getarg(i,tag)
          read (tag,'(BN,F10.0)') len
      elseif (tag(1:3).eq.'-nw') then
          i = i + 1
          call getarg(i,tag)
          read (tag,'(BN,I6)') nwid
      elseif (tag(1:3).eq.'-nl') then
          i = i + 1
          call getarg(i,tag)
          read (tag,'(BN,I6)') nlen
      elseif (tag(1:4).eq.'-str') then
          i = i + 1
          call getarg(i,tag)
          read (tag,'(BN,F10.0)') str
      elseif (tag(1:4).eq.'-dip') then
          i = i + 1
          call getarg(i,tag)
          read (tag,'(BN,F10.0)') dip
      elseif (tag(1:2).eq.'-h'.or.tag(1:2).eq.'-?') then
          call usage()
C      elseif (tag(1:2).eq.'-p') then
C          p = 1
C      elseif (tag(1:2).eq.'-u') then
C          user = 1
      endif
      goto 11

   12 continue

      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE usage()
      IMPLICIT none
      write(*,*)
     1 'Usage: grid2 -c [CLON CLAT CDEP] -w [WID] -l [LEN] -nw [NWID] ',
     2 '-nl [NLEN] -str [STR] -dip [DIP] -h/-?'
      write(*,*)
     1 '  -c [CLON CLAT CDEP] Define center point on fault'
      write(*,*)
     1 '  -w [WID]'
      write(*,*)
     1 '  -l [LEN]'
      write(*,*)
     1 '  -nw [NWID]'
      write(*,*)
     1 '  -nl [NLEN]'
      write(*,*)
     1 '  -str [STR]'
      write(*,*)
     1 '  -dip [DIP]'
C      write(*,*)
C     1 '  -p         (Default off) print output to stdout (will not ',
C     2                            'write file grid2.out)'
C      write (*,*)
C     1 '  -u         Prompt user to enter information through standard'
C      write (*,*)
C     1 '                 input for single calculation'
      write (*,*)
     1 '  -h/-?      Online help (this screen)'
      write (*,*) ''
      write (*,*)
     1 '  grid creates three files:'
      write (*,*)
     1 '    grid2.ll (lon, lat, dep)'
      write (*,*)
     1 '    grid2.xy (along strike dist, along dip dist, dep)'
      write (*,*)
     1 '    grid2.clip (clipping path for psclip)'
      write (*,*) ''

      STOP
      END
