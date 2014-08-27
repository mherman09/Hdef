      PROGRAM flt2displacement
C----
C Calculate the north, east and vertical displacements at receivers in
C an isotropic elastic half-space due to shear dislocations.
C----
      IMPLICIT NONE
      CHARACTER*40 fltfile,stafile,haffile,dspfile
      LOGICAL ex
      INTEGER auto
      INTEGER nflt,FMAX
      PARAMETER (FMAX=1500)
      INTEGER flttyp(FMAX)
      REAL*8 evlo(FMAX),evla(FMAX),evdp(FMAX),str(FMAX),dip(FMAX),
     1       rak(FMAX),slip(FMAX),dx(FMAX),dy(FMAX),area(FMAX)
      REAL*8 stlo,stla,stdp
      REAL*8 vp,vs,dens
      INTEGER prog,prog100,progtag
      REAL*8 uNnet,uEnet,uZnet

C----
C Parse the command line, check for FFM file, and read FFM
C----
      call gcmdln(fltfile,stafile,haffile,dspfile,auto)

C----
C Output file
C----
      open (unit=11,file='disp.out',status='unknown')
C----
C Parse fault file
C List of:
C  lon lat dep(km) str dip rak slip(m) dy(km) dx(km) flttyp
C Last column indicates fault type:
C   0 for point source
C   1 for finite source
C----
      call readfaults(nflt,evlo,evla,evdp,str,dip,rak,slip,dy,dx,area,
     1                flttyp,FMAX)
C----
C Calculate net displacements at each station
C----
  101 read (22,*,end=103) stlo,stla,stdp
          stdp = 1.0d3*stdp
          call calcdisp(uNnet,uEnet,uZnet,stlo,stla,stdp,nflt,
     1                  evlo,evla,evdp,str,dip,rak,dx,dy,area,slip,
     2                  vp,vs,dens,flttyp)
          write (11,9999) stlo,stla,uNnet,uEnet,uZnet
          goto 101
  103 continue

 9999 format (2(f10.3),3(f14.4))

      END

C======================================================================C

      SUBROUTINE readfaults(nflt,evlo,evla,evdp,str,dip,rak,slip,dy,dx,
     1                      area,flttyp,FMAX)

      IMPLICIT none
      INTEGER flt,nflt,FMAX,flttyp(FMAX)
      REAL*8 evlo(FMAX),evla(FMAX),evdp(FMAX),str(FMAX),
     1       dip(FMAX),rak(FMAX),slip(FMAX),dy(FMAX),dx(FMAX),
     2       area(FMAX)

      flt = 0
  201 flt = flt + 1
      read (21,*,END=202) evlo(flt),evla(flt),evdp(flt),str(flt),
     1                    dip(flt),rak(flt),slip(flt),dy(flt),dx(flt),
     2                    flttyp(flt)
          evdp(flt) = 1.0d3*evdp(flt)
          dx(flt) = 1.0d3*dx(flt)
          dy(flt) = 1.0d3*dy(flt)
          area(flt) = dx(flt)*dy(flt)
          goto 201
  202 continue
      nflt = flt - 1

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
      INTEGER flttyp(FMAX)
      REAL*8 vp,vs,dens

      uNnet = 0.0d0
      uEnet = 0.0d0
      uZnet = 0.0d0
      do 102 f = 1,nflt
          call ddistaz(dist,az,evlo(f),evla(f),stlo,stla)
          dist = dist*6.371d6
          x = dist*( dcos(az-d2r*str(f)))
          y = dist*(-dsin(az-d2r*str(f)))
          if (flttyp(f).eq.0) then
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

C----------------------------------------------------------------------C

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

      SUBROUTINE gcmdln(fltfile,stafile,haffile,dspfile,auto)
      IMPLICIT none
      CHARACTER*40 tag
      INTEGER narg,i
      CHARACTER*40 fltfile,stafile,haffile,dspfile
      INTEGER auto

      fltfile = 'faults.txt'
      stafile = 'stations.txt'
      haffile = 'structure.txt'
      dspfile = 'disp.out'
      auto = 0

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
      write (*,*)
     1 'Usage: flt2displacement -V -h/-?'
      write (*,*)
     1 '  -V (default false) turn on verbose operation'
      write (*,*)
     1 '  -h/-?              help'
      write (*,*) ''
      write (*,*)
     1 '  flt2displacement calculates NEZ displacements from faults at',
     2   ' user-defined locations.'
      write (*,*)
     1 '  Output file: disp.out'
      write (*,*)
     1 '  Required input files:'
      write (*,*) ''
      write (*,*)
     1 '    faults.txt: list of dislocation sources'
      write (*,*)
     1 '      evlo evla evdp str dip rak slip width length fault_type'
      write (*,*)
     1 '                (km)             (m)  (km)   (km)  (0=pt/1=fin)'
      write (*,*) ''
      write (*,*)
     1 '    stations.txt: list of station locations and depths'
      write (*,*)
     1 '      stlo stla stdp'
      write (*,*)
     1 '                (km)'
      write (*,*) ''
      write (*,*)
     1 '    structure.txt: vp vs dens (only vp/vs ratio matters)'
      write (*,*) ''
      stop
      END
