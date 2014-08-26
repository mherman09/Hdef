      PROGRAM ff2displacement
C----
C Calculate the north, east and vertical displacements at receivers in
C an isotropic elastic half-space due to shear dislocations from a
C finite fault model in standard subfault format.
C----
      IMPLICIT NONE
      CHARACTER*40 ffmfile,stafile,haffile,dspfile
      LOGICAL ex
      INTEGER flttyp,auto
      INTEGER nflt,FMAX
      PARAMETER (FMAX=1500)
      REAL*8 evlo(FMAX),evla(FMAX),evdp(FMAX),str(FMAX),dip(FMAX),
     1       rak(FMAX),slip(FMAX),dx(FMAX),dy(FMAX),area(FMAX),
     2       trup(FMAX),hypolon,hypolat
      REAL*8 stlo,stla,stdp
      REAL*8 vp,vs,dens
      INTEGER prog,prog100,progtag
      REAL*8 uNnet,uEnet,uZnet

C----
C Parse the command line, check for FFM file, and read FFM
C----
      call gcmdln(ffmfile,stafile,haffile,dspfile,flttyp,auto)
      inquire(file=ffmfile,EXIST=ex)
      if (.not.ex) call usage('!! Error: no file named '//ffmfile)
      call readffm(ffmfile,nflt,hypolon,hypolat,evlo,evla,evdp,str,dip,
     1             rak,dx,dy,area,slip,trup)

C----
C If auto=0: use existing input files
C If auto=1: generate half-space and station inputs
C----
      if (auto.eq.0) then
          inquire(file=haffile,EXIST=ex)
          if (.not.ex) then
              write(0,*) '!! Warning: no file named '//haffile
              write(0,*) '!! Using default half-space parameters'
              call sethafdef(vp,vs,dens)
          else
              open (unit=23,file=haffile,status='old')
              read (23,*) vp,vs,dens
          endif
          inquire(file=stafile,EXIST=ex)
          if (.not.ex) call usage('!! Error: no file named '//stafile)
          open (unit=22,file=stafile,status='old')
      elseif (auto.eq.1) then
          call sethafdef(vp,vs,dens)
          call autogrid(stafile,nflt,hypolon,hypolat,evlo,evla,evdp,str,
     1                  dip,rak,dx,dy,area,slip,vp,vs,dens,flttyp)
          open (unit=22,file=stafile,status='old')
      endif

C----
C Output file
C----
      open (unit=10,file=dspfile,status='unknown')
      rewind 10

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
          write (10,9999) stlo,stla,uEnet,uNnet,uZnet
          prog = prog + 1
          call progbar(prog,prog100,progtag)
          goto 101
  103 continue
 9999 format (2F12.4,X,3F12.4)

      END

C======================================================================C

      SUBROUTINE readffm(ffmfile,nflt,hypolon,hypolat,evlo,evla,evdp,
     1                   str,dip,rak,dx,dy,area,slip,trup)
C----
C Read FFM in standard subfault format.
C   - number of fault segments (nseg) is in first line
C   - number (nx, ny) and dimensions (dx, dy) of subfaults and epicenter
C     location (hypolon, hypolat) are in 9-line segment headers.
C   - lat lon dep(km) slip(cm) rak str dip t_rup t_ris mo
C----
      IMPLICIT none
      CHARACTER*40 ffmfile,dum
      INTEGER nflt,FMAX
      PARAMETER (FMAX=1500)
      REAL*8 evlo(FMAX),evla(FMAX),evdp(FMAX),str(FMAX),dip(FMAX),
     1       rak(FMAX),dx(FMAX),dy(FMAX),area(FMAX),slip(FMAX),
     2       trup(FMAX),hypolon,hypolat
      INTEGER seg,nseg,ct,i,nx,ny,ptr
      CHARACTER*10 dxc,dyc
      REAL*8 dxr,dyr

      open (unit=21,file=ffmfile,status='old')
      ct = 0
      read (21,*) dum,dum,dum,dum,nseg
      do 103 seg = 1,nseg
          read (21,*) dum,dum,dum,dum,nx,dum,dxc,dum,ny,dum,dyc
          ptr = index(dxc,'km')
          dxc(ptr:ptr+1) = ''
          ptr = index(dyc,'km')
          dyc(ptr:ptr+1) = ''
          read (dxc,*) dxr
          read (dyc,*) dyr
          dxr = 1.0d3*dxr
          dyr = 1.0d3*dyr
          read (21,*) dum,dum,dum,dum,dum,dum,dum,dum,dum,dum,hypolon,
     1                dum,hypolat
          do 101 i=1,7
              read (21,*) dum
  101     continue
          do 102 i = 1,nx*ny
              read (21,*) evla(ct+i),evlo(ct+i),evdp(ct+i),slip(ct+i),
     1                    rak(ct+i),str(ct+i),dip(ct+i),trup(ct+i)
              slip(ct+i) = 1.0d-2*slip(ct+i)
              dx(ct+i) = dxr
              dy(ct+i) = dyr
              area(ct+i) = dx(ct+i)*dy(ct+i)
              evdp(ct+i) = 1.0d3*evdp(ct+i)
  102     continue
          ct = ct + nx*ny
  103 continue
      nflt = ct
      close(21)

      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE calcdisp(uNnet,uEnet,uZnet,stlo,stla,stdp,nflt,
     1                    evlo,evla,evdp,str,dip,rak,dx,dy,area,slip,
     2                    vp,vs,dens,flttyp)
C----
C Compute NEZ displacements at (stlo,stla,stdp) from FFM
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

      SUBROUTINE autogrid(stafile,nflt,hypolon,hypolat,evlo,evla,evdp,
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
      INTEGER nflt,FMAX,i,j,NTHR,flttyp,dir
      PARAMETER (FMAX=1500,NTHR=3)
      REAL*8 evlo(FMAX),evla(FMAX),evdp(FMAX),str(FMAX),dip(FMAX),
     1       rak(FMAX),slip(FMAX),area(FMAX),dx(FMAX),dy(FMAX),
     2       dl(NTHR),thr(NTHR)
      REAL*8 stdp,lon,lat,hypolon,hypolat,az,Nmax,Smax,Emax,Wmax,
     1       vp,vs,dens,dN,dE,dZ,N,S,E,W,incr,distmax

      write(STER,'(A)') 'Determining receiver grid limits...'
C----
C Limit maximum extent of receiver grid
C NB: Nmax and Smax are returned in RADIANS!
C----
      distmax = 750.0d0 ! km
      az = 0.0d0
      call dlola(lon,Nmax,hypolon,hypolat,distmax,az)
      az = 180.0d0
      call dlola(lon,Smax,hypolon,hypolat,distmax,az) 
      Emax = hypolon+distmax*3.6d2/(2.0d0*pi*6.371d3*dcos(hypolat*d2r))
      Wmax = hypolon-distmax*3.6d2/(2.0d0*pi*6.371d3*dcos(hypolat*d2r))
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
          lon = hypolon
          lat = hypolat
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
      incr = 0.10d0
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
      write(STER,9003) dnint((E-W)*(N-S)/(incr*incr))
      write(52,9001) W,E,S,N
      write(52,9002) incr
      write(52,9003) dnint((E-W)*(N-S)/(incr*incr))
      write(52,9004) i
      write(52,9005) j
 8001 format(3F10.2)
 9001 format('Grid limits (WESN):',4F8.2)
 9002 format('Grid spacing: ',F10.2)
 9003 format('Number of points: ',F10.0)
 9004 format('nlon: ',I10)
 9005 format('nlat: ',I10)
      close(51)

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
 1000 format ('ff2displacement: [',I3,'% Complete]',A)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE gcmdln(ffmfile,stafile,haffile,dspfile,flttyp,auto)
      IMPLICIT none
      CHARACTER*40 tag
      INTEGER narg,i
      CHARACTER*40 ffmfile,stafile,haffile,dspfile
      INTEGER flttyp,auto

      ffmfile = 'static_out'
      stafile = 'stations.txt'
      haffile = 'structure.txt'
      dspfile = 'disp.out'
      flttyp = 1
      auto = 0

      narg = iargc()
      i = 0
  101 i = i + 1
      if (i.gt.narg) goto 102
          call getarg(i,tag)
          if (tag(1:2).eq.'-h'.or.tag(1:2).eq.'-?') then
              call usage(' ')
          elseif (tag(1:4).eq.'-ffm') then
              i = i + 1
              call getarg(i,ffmfile)
          elseif (tag(1:4).eq.'-sta') then
              i = i + 1
              call getarg(i,stafile)
          elseif (tag(1:3).eq.'-haf') then
              i = i + 1
              call getarg(i,haffile)
          elseif (tag(1:4).eq.'-dsp') then
              i = i + 1
              call getarg(i,dspfile)
          elseif (tag(1:3).eq.'-fn') then
              flttyp = 1
          elseif (tag(1:3).eq.'-pt') then
              flttyp = 0
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
     1 'Usage: ff2displacement -ffm FFMFILE -sta STAFILE ',
     2                        '-haf HAFFILE -dsp DSPFILE ',
     3                        '-fn -pt -auto ',
     4                        '-h -?'
      write(STER,*)
      write(STER,*)
     1 'Compute NEZ displacements resulting from finite fault model'
      write(STER,*)
      write (STER,*)
     1 '-ffm FFMFILE   (static_out)     name of finite fault model file'
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
     1 '-h                              this online help'
      write (STER,*)
     1 '-?                              this online help'
      write (STER,*)
      write (STER,*)
     1 'FILE FORMATS'
      write (STER,*)
     1 '------------'
      write (STER,*)
     1 'FFMFILE: finite fault model in standard subfault format'
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
      write (STER,*)
     1 'If "-auto" is selected, a receiver grid is automatically' 
      write (STER,*)
     1 'generated, and surface displacements are computed at each site.'
      write (STER,*)
     1 'A file "autogrid.dat" contains information about the grid.'
      write (STER,*)

      STOP
      END
