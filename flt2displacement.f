      PROGRAM flt2displacement
C----
C Calculate the north, east and vertical displacements at receivers in
C an isotropic elastic half-space due to shear dislocations.
C----
      IMPLICIT NONE
      INTEGER STER
      PARAMETER (STER=0)
      CHARACTER*40 fltfile,mtfile,stafile,haffile,dspfile,gmtfile
      LOGICAL ex
      INTEGER auto,mt,gmt
      INTEGER nflt,FMAX
      PARAMETER (FMAX=1500)
      INTEGER flttyp
      REAL*8 evlo(FMAX),evla(FMAX),evdp(FMAX),str(FMAX),dip(FMAX),
     1       rak(FMAX),slip(FMAX),dx(FMAX),dy(FMAX),area(FMAX)
      REAL*8 stlo,stla,stdp
      REAL*8 vp,vs,dens
      INTEGER prog,prog100,progtag
      REAL*8 uNnet,uEnet,uZnet

C----
C Parse the command line, check for FFM file, and read FFM
C----
      call gcmdln(fltfile,mtfile,stafile,haffile,dspfile,flttyp,auto,mt,
     1            gmtfile,gmt)
      if (mt.eq.0) then
          inquire(file=fltfile,exist=ex)
          if (.not.ex) call usage('!! Error: no file named '//fltfile)
      else
          inquire(file=mtfile,exist=ex)
          if (.not.ex) call usage('!! Error: no file named '//mtfile)
          call mt2flt(fltfile,mtfile)
      endif
      call readfaults(fltfile,nflt,evlo,evla,evdp,str,dip,rak,slip,
     1                dy,dx,area)

C----
C If auto=0: use existing input files
C If auto=1: generate half-space and station inputs
C----
      if (auto.eq.0) then
          inquire(file=haffile,exist=ex)
          if (.not.ex) then
              write(STER,*) '!! Warning: no file named '//haffile
              write(STER,*) '!! Using default half-space parameters'
              call sethafdef(vp,vs,dens)
          else
              open (unit=23,file=haffile,status='old')
              read (23,*) vp,vs,dens
              close(23)
          endif
          inquire(file=stafile,exist=ex)
          if (.not.ex) call usage('!! Error: no file named '//stafile)
          open(unit=22,file=stafile,status='old')
      elseif (auto.eq.1) then
          call sethafdef(vp,vs,dens)
          call autogrid(stafile,nflt,evlo,evla,evdp,str,
     1                  dip,rak,dx,dy,area,slip,vp,vs,dens,flttyp)
          open (unit=22,file=stafile,status='old')
      endif

C----
C Output file
C----
      open (unit=11,file=dspfile,status='unknown')
      rewind 11

C----
C Initialize progress indicator
C----
      prog100 = 0
  112 read (22,*,end=113)
          prog100 = prog100 + 1
          goto 112
  113 continue
      rewind 22
      prog = 0
      progtag = 0
      call progbar(prog,prog100,progtag)

C----
C Calculate net displacements at each station
C----
  101 read (22,*,end=103) stlo,stla,stdp
          stdp = 1.0d3*stdp
          call calcdisp(uNnet,uEnet,uZnet,stlo,stla,stdp,nflt,
     1                  evlo,evla,evdp,str,dip,rak,dx,dy,area,slip,
     2                  vp,vs,dens,flttyp)
          write (11,9999) stlo,stla,uEnet,uNnet,uZnet
          prog = prog + 1
          call progbar(prog,prog100,progtag)
          goto 101
  103 continue
 9999 format (2F12.4,X,3F12.4)

      if (gmt.eq.1) call makegmtfile(gmtfile,nflt,evlo,evla,str,dip,dx,
     1                               dy,slip)

      END

C======================================================================C

      SUBROUTINE readfaults(fltfile,nflt,evlo,evla,evdp,str,dip,rak,
     1                      slip,dy,dx,area)

      IMPLICIT none
      CHARACTER*40 fltfile
      INTEGER f,nflt,FMAX
      PARAMETER (FMAX=1500)
      REAL*8 evlo(FMAX),evla(FMAX),evdp(FMAX),str(FMAX),
     1       dip(FMAX),rak(FMAX),slip(FMAX),dy(FMAX),dx(FMAX),
     2       area(FMAX)

      open (unit=21,file=fltfile,status='old')
      f = 0
  201 f = f + 1
      read (21,*,END=202) evlo(f),evla(f),evdp(f),str(f),
     1                    dip(f),rak(f),slip(f),dy(f),dx(f)
          evdp(f) = 1.0d3*evdp(f)
          dx(f) = 1.0d3*dx(f)
          dy(f) = 1.0d3*dy(f)
          !print *,'FAULT',f
          !print *,'EVLO',evlo(f),'EVLA',evla(f),'EVDP',evdp(f)
          !print *,'STR ',str(f), 'DIP ',dip(f), 'RAK ',rak(f)
          !print *,'DX  ',dx(f),  'DY  ',dy(f),  'SLIP',slip(f)
          area(f) = dx(f)*dy(f)
          goto 201
  202 continue
      nflt = f - 1
      close(21)

      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE calcdisp(uNnet,uEnet,uZnet,stlo,stla,stdp,nflt,
     1                    evlo,evla,evdp,str,dip,rak,dx,dy,area,slip,
     2                    vp,vs,dens,flttyp)
C----
C Compute NEZ displacements at (stlo,stla,stdp) from faults
C----
      IMPLICIT none
      REAL*8 pi,d2r
      PARAMETER (pi=4.0d0*datan(1.0d0),d2r=pi/1.8d2)
      REAL*8 uNnet,uEnet,uZnet,ux,uy,uN,uE,uz
      REAL*8 stlo,stla,stdp,dist,az,x,y
      INTEGER f,nflt,FMAX
      PARAMETER (FMAX=1500)
      REAL*8 evlo(FMAX),evla(FMAX),evdp(FMAX),str(FMAX),dip(FMAX),
     1       rak(FMAX),slip(FMAX),area(FMAX),dx(FMAX),dy(FMAX)
      INTEGER flttyp
      REAL*8 vp,vs,dens

      uNnet = 0.0d0
      uEnet = 0.0d0
      uZnet = 0.0d0
      do 102 f = 1,nflt
          call ddistaz(dist,az,evlo(f),evla(f),stlo,stla)
          dist = dist*6.371d6
          x = dist*( dcos(az-d2r*str(f)))
          y = dist*(-dsin(az-d2r*str(f)))
          if (flttyp.eq.0) then
              call o92pt(ux,uy,uz,x,y,stdp,evdp(f),dip(f),rak(f),
     1                   area(f),slip(f),vp,vs,dens)
          else
              call o92rect(ux,uy,uz,x,y,stdp,evdp(f),dip(f),rak(f),
     1                     dy(f),dx(f),slip(f),vp,vs,dens)
          endif
          call xy2NE(uN,uE,ux,uy,str(f))
          uNnet = uNnet + uN
          uEnet = uEnet + uE
          uZnet = uZnet + uz
  102 continue
      RETURN
      END

C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - C

      SUBROUTINE xy2NE(uN,uE,ux,uy,str)
      IMPLICIT none
      REAL*8 uN,uE,ux,uy,str,theta,uhor
      REAL*8 pi,d2r
      PARAMETER (pi=4.0d0*datan(1.0d0),d2r=pi/1.8d2)
      theta = datan2(uy,ux)
      uhor = dsqrt(ux*ux+uy*uy)
      theta = d2r*str - theta
      uN = uhor*dcos(theta)
      uE = uhor*dsin(theta)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE sethafdef(vp,vs,dens)
      IMPLICIT none
      REAL*8 vp,vs,dens
      vp = 6800.0d0
      vs = 3926.0d0
      dens = 2800.0d0
      RETURN
      END

C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - C

      SUBROUTINE autogrid(stafile,nflt,evlo,evla,evdp,
     1                    str,dip,rak,dx,dy,area,slip,vp,vs,dens,flttyp)
C----
C Automatically determine grid of stations to calculate displacements
C----
      IMPLICIT none
      REAL*8 pi,r2d,d2r
      PARAMETER (pi=4.0d0*atan(1.0d0),r2d=180.0d0/pi,d2r=pi/180.0d0)
      INTEGER STER
      PARAMETER (STER=0)
      CHARACTER*40 stafile
      INTEGER nflt,FMAX,i,j,NTHR,dir
      PARAMETER (FMAX=1500,NTHR=3)
      INTEGER flttyp
      REAL*8 evlo(FMAX),evla(FMAX),evdp(FMAX),str(FMAX),dip(FMAX),
     1       rak(FMAX),slip(FMAX),area(FMAX),dx(FMAX),dy(FMAX),
     2       dl(NTHR),thr(NTHR)
      REAL*8 stdp,lon,lat,avglon,avglat,az,Nmax,Smax,Emax,Wmax,
     1       vp,vs,dens,dN,dE,dZ,N,S,E,W,incr,distmax

      write(STER,'(A)') 'Determining receiver grid limits...'
      do 110 i = 1,nflt
          avglon = avglon + evlo(i)
          avglat = avglat + evla(i)
  110 continue
      avglon = avglon/dble(nflt)
      avglat = avglat/dble(nflt)
C----
C Limit maximum extent of receiver grid
C NB: Nmax and Smax are returned in RADIANS!
C----
      distmax = 750.0d0 ! km
      az = 0.0d0
      call dlola(lon,Nmax,avglon,avglat,distmax,az)
      az = 180.0d0
      call dlola(lon,Smax,avglon,avglat,distmax,az)
      Emax = avglon+distmax*3.6d2/(2.0d0*pi*6.371d3*dcos(avglat*d2r))
      Wmax = avglon-distmax*3.6d2/(2.0d0*pi*6.371d3*dcos(avglat*d2r))
C----
C Set up distance steps and corresponding displacement thresholds
C----
      dl(1) = 1.0d0 ! Amounts to move test location
      dl(2) = 0.1d0
      dl(3) = 0.05d0
      thr(1) = 0.5d0 ! Displacement thresholds
      thr(2) = 0.01d0
      thr(3) = 0.001d0
C----
C Find limits to north, south, east, and west
C Variable 'dir' controls direction: 1=N, 2=S, 3=E, 4=W
C----
      stdp = 0.0d0
      do 103 dir = 1,4
          lon = avglon
          lat = avglat
C         Take first step (by increment of dl(1))
          if (dir.eq.1) lat = dnint((lat+dl(1))*1.d2)*1.d-2
          if (dir.eq.2) lat = dnint((lat-dl(1))*1.d2)*1.d-2
          if (dir.eq.3) lon = dnint((lon+dl(1))*1.d2)*1.d-2
          if (dir.eq.4) lon = dnint((lon-dl(1))*1.d2)*1.d-2
C         Take steps of decreasing size and threshold value
          i = 1
  101     if (i.gt.NTHR) goto 102
C             Check if lon, lat exceed maximum allowed values
              if (dir.eq.1.and.lat.gt.Nmax*r2d) then
                  lat = Nmax*r2d
                  goto 102
              elseif (dir.eq.2.and.lat.lt.Smax*r2d) then
                  lat = Smax*r2d
                  goto 102
              elseif (dir.eq.3.and.lon.gt.Emax) then
                  lat = Emax
                  goto 102
              elseif (dir.eq.4.and.lon.lt.Wmax) then
                  lon = Wmax
                  goto 102
              endif
C             Calculate displacement and compare to current threshold
              call calcdisp(dN,dE,dZ,lon,lat,stdp,nflt,evlo,evla,evdp,
     1                      str,dip,rak,dx,dy,area,slip,vp,vs,dens,
     2                      flttyp)
              dN = dabs(dN)
              dE = dabs(dE)
              dZ = dabs(dZ)
              if (dN.gt.thr(i).or.dE.gt.thr(i).or.dZ.gt.thr(i)) then
                  if (dir.eq.1) lat = lat + dl(i)
                  if (dir.eq.2) lat = lat - dl(i)
                  if (dir.eq.3) lon = lon + dl(i)
                  if (dir.eq.4) lon = lon - dl(i)
                  goto 101
              else
                  i = i + 1
                  goto 101
              endif
  102     continue
          if (dir.eq.1) N = lat
          if (dir.eq.2) S = lat
          if (dir.eq.3) E = lon
          if (dir.eq.4) W = lon
  103 continue
C----
C Create gridded stafile with separation between points = 'incr'
C----
      incr = 0.05d0
      open(unit=51,file=stafile,status='unknown')
      i = 0
      lon = W
  104 if (lon.gt.E) goto 107
          i = i + 1
          lat = S
          j = 0
  105     if (lat.gt.N) goto 106
              write(51,8001) lon,lat,stdp
              lat = lat + incr
              j = j + 1
              goto 105
  106     continue
          lon = lon + incr
          goto 104
  107 continue
      open(unit=52,file='autogrid.dat',status='unknown')
      write(STER,9001) W,E,S,N
      write(STER,9002) incr
      write(STER,9003) i*j
      write(52,9001) W,E,S,N
      write(52,9002) incr
      write(52,9003) i*j
      write(52,9004) i
      write(52,9005) j
 8001 format(3F10.2)
 9001 format('Grid limits (WESN):',4F8.2)
 9002 format('Grid spacing: ',F10.2)
 9003 format('Number of points: ',I10)
 9004 format('nlon: ',I10)
 9005 format('nlat: ',I10)
      close(51)

      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE mt2flt(fltfile,mtfile)
      IMPLICIT none
      CHARACTER*40 fltfile,mtfile
      REAL*8 lon,lat,dep,str,dip,rak,mag,wid,len,vp,vs,dens,mu,slip
      INTEGER typ
      call sethafdef(vp,vs,dens)
      mu = vs*vs*dens
      open(unit=61,file=mtfile,status='old')
      open(unit=62,file=fltfile,status='unknown')
  101 read(61,*,end=102) lon,lat,dep,str,dip,rak,mag
          if (45.0d0.lt.rak.and.rak.lt.135.0d0) then
              typ = 2
          elseif (-135.0d0.lt.rak.and.rak.lt.-45.0d0) then
              typ = 3
          else
              typ = 1
          endif
          call wellscoppersmith(wid,len,mag,typ)
          slip = ((1.d1**(1.5d0*(mag+10.7d0)))*1.d-7)/(mu*wid*len*1.d6)
          write(62,*)lon,lat,dep,str,dip,rak,slip,wid,len
          goto 101
  102 continue
      close(61)
      close(62)
      RETURN
      END

C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - C

      SUBROUTINE wellscoppersmith(wid,len,mag,typ)
      IMPLICIT none
      INTEGER typ
      REAL*8 wid,len,mag

      if (typ.eq.1) then
          len = 10.0d0**(-2.57d0+0.62d0*mag)
          wid = 10.0d0**(-0.76d0+0.27d0*mag)
      elseif (typ.eq.2) then
          len = 10.0d0**(-2.42d0+0.58d0*mag)
          wid = 10.0d0**(-1.61d0+0.41d0*mag)
      elseif (typ.eq.3) then
          len = 10.0d0**(-1.88d0+0.50d0*mag)
          wid = 10.0d0**(-1.14d0+0.35d0*mag)
      else
          len = 10.0d0**(-2.44d0+0.59d0*mag)
          wid = 10.0d0**(-1.01d0+0.32d0*mag)
      endif

      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE makegmtfile(gmtfile,nflt,evlo,evla,str,dip,dx,dy,slip)
      IMPLICIT none
      REAL*8 pi,d2r
      PARAMETER (pi=4.0d0*datan(1.0d0),d2r=pi/180.0d0)
      CHARACTER*40 gmtfile
      INTEGER i,nflt,FMAX
      PARAMETER (FMAX=1500)
      REAL*8 evlo(FMAX),evla(FMAX),str(FMAX),dip(FMAX),
     1       slip(FMAX),dx(FMAX),dy(FMAX),cosd
      open(unit=71,file=gmtfile,status='unknown')
      do 301 i = 1,nflt
          cosd = dcos(dip(i)*d2r)
          write(71,1001) evlo(i),evla(i),slip(i),str(i),
     1                   dx(i)*1.0d-3,dy(i)*1.0d-3*cosd
  301 continue
      close(71)
 1001 format(2F10.3,F8.2,F8.1,2F10.3)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE progbar(prog,prog100,progtag)
      IMPLICIT none
      INTEGER prog,prog100,progtag
      CHARACTER*1 CR
      CR = char(13) ! Carriage return
      if (100*prog/prog100.ge.progtag) then
          write (*,1000,advance='no') 100*prog/prog100,CR
          progtag = progtag + 1
      endif
      if (100*prog/prog100.ge.100) write (*,1000) 100*prog/prog100
 1000 format ('flt2displacement: [',I3,'% Complete]',A)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE gcmdln(fltfile,mtfile,stafile,haffile,dspfile,flttyp,
     1                  auto,mt,gmtfile,gmt)
      IMPLICIT none
      CHARACTER*40 tag
      INTEGER narg,i
      CHARACTER*40 fltfile,mtfile,stafile,haffile,dspfile,gmtfile
      INTEGER auto,mt,flttyp,gmt

      fltfile = 'faults.txt'
      stafile = 'stations.txt'
      haffile = 'structure.txt'
      dspfile = 'disp.out'
      mtfile = 'none'
      gmtfile = 'none'
      mt = 0
      auto = 0
      gmt = 0
      flttyp = 1

      narg = iargc()
      i = 0
  101 i = i + 1
      if (i.gt.narg) goto 102
          call getarg(i,tag)
          if (tag(1:2).eq.'-h'.or.tag(1:2).eq.'-?') then
              call usage(' ')
          elseif (tag(1:4).eq.'-flt') then
              i = i + 1
              call getarg(i,fltfile)
          elseif (tag(1:4).eq.'-sta') then
              i = i + 1
              call getarg(i,stafile)
          elseif (tag(1:3).eq.'-haf') then
              i = i + 1
              call getarg(i,haffile)
          elseif (tag(1:4).eq.'-dsp') then
              i = i + 1
              call getarg(i,dspfile)
          elseif (tag(1:3).eq.'-mt') then
              mt = 1
              i = i + 1
              call getarg(i,mtfile)
          elseif (tag(1:3).eq.'-fn') then
              flttyp = 1
          elseif (tag(1:3).eq.'-pt') then
              flttyp = 0
          elseif (tag(1:4).eq.'-gmt') then
              gmt = 1
              i = i + 1
              call getarg(i,gmtfile)
          elseif (tag(1:5).eq.'-auto') then
              !i = i + 1
              !call getarg(i,tag)
              !read(tag,'(BN,I10)') auto
              auto = 1
          else
              call usage('!! Error: no option '//tag)
          endif
          goto 101
  102 continue

      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE usage(str)
      IMPLICIT none
      INTEGER STER,lstr
      PARAMETER (STER=0)
      CHARACTER str*(*)
      if (str.ne.' ') then
          lstr = len(str)
          write(STER,*)
          write(STER,*) str(1:lstr)
          write(STER,*)
      endif
      write (STER,*)
     1 'Usage: ff2displacement -flt FLTFILE -mt MTFILE -sta STAFILE ',
     2                        '-haf HAFFILE -dsp DSPFILE'
      write (STER,*)
     3 '                       -fn -pt -auto -gmt GMTFILE ',
     4                        '-h -?'
      write(STER,*)
      write(STER,*)
     1 'Compute NEZ displacements resulting from shear dislocations'
      write(STER,*)
     1 'Note: The "-mt MTFILE" option overrides "-flt FLTFILE" option'
      write(STER,*)
      write (STER,*)
     1 '-ffm FFMFILE   (faults.txt)     name of fault file'
      write (STER,*)
     1 '-mt  MTFILE                     name of fault file (psmeca ',
     2                                    '-Sa format)'
      write (STER,*)
     1 '-sta STAFILE   (stations.txt)   name of receiver location file'
      write (STER,*)
     1 '-haf HAFFILE   (structure.txt)  name of half-space parameter ',
     2                                  'file'
      write (STER,*)
     1 '-dsp DSPFILE   (disp.out)       name of output displacement ',
     2                                  'file'
      write (STER,*)
     1 '-fn            (default)        treat subfaults as ',
     2                                    'finite sources'
      write (STER,*)
     1 '-pt                             treat subfaults as ',
     2                                       'point sources'
      write (STER,*)
     1 '-auto                           automatically create ',
     2                               'receiver grid'
      write (STER,*)
     1 '-gmt GMTFILE                    create file for use with GMT ',
     2                                         '"psxy -SJ"'
      write (STER,*)
     1 '-h                              this online help'
      write (STER,*)
     1 '-?                              this online help'
      write (STER,*)
      write (STER,*)
     1 'FILE FORMATS'
      write (STER,*)
     1 '------------'
      write (STER,*)
     1 'FLTFILE: list of shear dislocations'
      write (STER,*)
     1 '    evlo evla evdp(km) str dip rak slip(m) wid(km) len(km)'
      write (STER,*)
     1 'MTFILE: list of shear dislocations'
      write (STER,*)
     1 '    evlo evla evdp(km) str dip rak mag'
      write (STER,*)
     1 'STAFILE: list of station/receiver locations and depths'
      write (STER,*)
     1 '    stlo stla stdp(km)'
      write (STER,*)
     1 'HAFFILE: half-space parameters'
      write (STER,*)
     1 '    vp(km/s) vs(km/s) dens(kg/m^3)'
      write (STER,*)
     1 'DSPFILE: list of station locations (corresponds to STAFILE) ',
     2           'and displacements'
      write (STER,*)
     1 '    stlo stla uE(m) uN(m) uZ(m)'
      write (STER,*)
     1 'GMTFILE: list of fault patches for use with "psxy -SJ"'
      write (STER,*)
     1 '    evlo evla slip(m) str len(km) hor_wid(km)'
      write (STER,*)
      write (STER,*)
     1 'If "-auto" is selected, a receiver grid is automatically'
      write (STER,*)
     1 'generated, and surface displacements are computed at each site.'
      write (STER,*)
     1 'A file "autogrid.dat" contains information about the grid.'
      write (STER,*)
      stop
      END
