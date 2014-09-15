      PROGRAM o92util
C----
C Utility for computing displacements, strains, and stresses in an
C elastic half-space resulting from rectangular shear dislocations.
C----
      IMPLICIT none
      CHARACTER*40 ffmfile,fltfile,magfile,stafile,haffile,trgfile,
     1             dspfile,stnfile,stsfile,norfile,shrfile,coufile,
     2             gmtfile
      INTEGER flttyp,auto
      REAL*8 incr
      INTEGER nflt,FMAX
      PARAMETER (FMAX=1500)
      REAL*8 evlo(FMAX),evla(FMAX),evdp(FMAX),str(FMAX),dip(FMAX),
     1       rak(FMAX),dx(FMAX),dy(FMAX),slip(FMAX),hylo,hyla
      INTEGER nsta,ntrg
      REAL*8 vp,vs,dens
      INTEGER kffm,kflt,kmag,ksta,khaf,ktrg
      INTEGER kdsp,kstn,ksts,knor,kshr,kcou,kgmt
      COMMON /ICHECK/ kffm,kflt,kmag,ksta,khaf,ktrg
      COMMON /OCHECK/ kdsp,kstn,ksts,knor,kshr,kcou,kgmt

C----
C Parse command line
C----
      call gcmdln(ffmfile,fltfile,magfile,stafile,haffile,trgfile,
     1            dspfile,stnfile,stsfile,norfile,shrfile,coufile,
     2            gmtfile,flttyp,auto,incr)
c      call dbggcmdln(ffmfile,fltfile,magfile,stafile,haffile,
c     1                     trgfile,dspfile,stnfile,stsfile,norfile,
c     2                     shrfile,coufile,gmtfile,flttyp,auto,incr)

C----
C Check for input fault source files
C----
      call chksrc(ffmfile,fltfile,magfile)

C----
C Check if any output is specified
C----
      call chkout()

C----
C Read input fault files
C----
      call readsrc(ffmfile,fltfile,magfile,evlo,evla,evdp,str,dip,rak,
     1             dx,dy,slip,hylo,hyla,nflt)
c      call dbgreadsrc(ffmfile,fltfile,magfile,evlo,evla,evdp,str,dip,
c     1                rak,dx,dy,slip,hylo,hyla,nflt)

C----
C Check for auto flag
C     auto=0: user-defined receiver, half-space, and target parameters
C     auto=1: make receiver grid and write to STAFILE
C             use default half-space parameters
C----
      if (auto.eq.1) then
          call autodispgrid(stafile,evlo,evla,evdp,str,dip,rak,dx,dy,
     1                      slip,hylo,hyla,nflt,incr)
      endif
      call chkstahaftrg(stafile,haffile,trgfile,nsta,ntrg)
      call readhaffile(haffile,vp,vs,dens)

C----
C Compute displacements, strains, and stresses
C----
      call computedefm(stafile,trgfile,nsta,ntrg,
     1                 nflt,evlo,evla,evdp,str,dip,rak,dx,dy,slip,
     2                 vp,vs,dens,flttyp,
     3                       dspfile,stnfile,stsfile,norfile,shrfile,
     4                       coufile)

      END

C----------------------------------------------------------------------C

      SUBROUTINE chksrc(ffmfile,fltfile,magfile)
      IMPLICIT none
      CHARACTER*40 ffmfile,fltfile,magfile
      LOGICAL ex
      INTEGER kffm,kflt,kmag,ksta,khaf,ktrg
      COMMON /ICHECK/ kffm,kflt,kmag,ksta,khaf,ktrg
C----
C Check for existence of fault source files
C----
      if (kffm.eq.0.and.kflt.eq.0.and.kmag.eq.0) then
          write(*,*) '!! Error: No fault file specified'
          call usage('!! Use -ffm FFMFILE, -flt FLTFILE, or -mag '//
     1                  'MAGFILE')
      endif
      if (kffm.eq.1) then
          inquire(file=ffmfile,EXIST=ex)
          if (.not.ex) kffm = -1
      endif
      if (kflt.eq.1) then
          inquire(file=fltfile,EXIST=ex)
          if (.not.ex) kflt = -1
      endif
      if (kmag.eq.1) then
          inquire(file=magfile,EXIST=ex)
          if (.not.ex) kmag = -1
      endif
      if (kffm.ne.1.and.kflt.ne.1.and.kmag.ne.1) then
          write(*,*) '!! Error: No fault files found'
          if (kffm.eq.-1) write(*,*) '!! Looked for FFMFILE: '//ffmfile
          if (kflt.eq.-1) write(*,*) '!! Looked for FLTFILE: '//fltfile
          if (kmag.eq.-1) write(*,*) '!! Looked for MAGFILE: '//magfile
          call usage('!! Check -ffm, -flt, and -mag arguments')
      else
          if (kffm.eq.-1) write(*,*) '!! Warning: No FFMFILE: '//ffmfile
          if (kflt.eq.-1) write(*,*) '!! Warning: No FLTFILE: '//fltfile
          if (kmag.eq.-1) write(*,*) '!! Warning: No MAGFILE: '//magfile
      endif
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE chkout()
      IMPLICIT none      
      INTEGER kdsp,kstn,ksts,knor,kshr,kcou,kgmt
      COMMON /OCHECK/ kdsp,kstn,ksts,knor,kshr,kcou,kgmt
      if (kdsp.eq.0.and.kstn.eq.0) then
          write(*,*) '!! Error: No output file specified'
          call usage('!! Use -disp, -strain, -stress, '//
     1                   '-normal, -shear, or -coul')
      endif
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE chkstahaftrg(stafile,haffile,trgfile,nsta,ntrg)
      IMPLICIT none
      CHARACTER*40 stafile,haffile,trgfile
      LOGICAL ex
      INTEGER nsta,ntrg
      INTEGER kffm,kflt,kmag,ksta,khaf,ktrg
      INTEGER kdsp,kstn,ksts,knor,kshr,kcou,kgmt
      COMMON /ICHECK/ kffm,kflt,kmag,ksta,khaf,ktrg
      COMMON /OCHECK/ kdsp,kstn,ksts,knor,kshr,kcou,kgmt
C----
C Check for station and half-space files
C----
      if (ksta.eq.0) then
          write(*,*) '!! Error: No station file specified'
          call usage('!! Use -sta STAFILE or one of the -auto options')
      elseif (ksta.eq.1) then
          inquire(file=stafile,EXIST=ex)
          if (.not.ex) then
              write(*,*) '!! Error: No station file found'
              call usage('!! Looked for STAFILE: '//stafile)
          endif
      endif
      call linecount(stafile,nsta)
      if (khaf.eq.0) then
C          write(*,*) '!! Using default half-space parameters'
      elseif (khaf.eq.1) then
          inquire(file=haffile,EXIST=ex)
          if (.not.ex) then
              write(*,*) '!! Warning: No half-space file found'
              write(*,*) '!! Looked for HAFFILE: '//HAFFILE
              write(*,*) '!! Using default half-space parameters'
              khaf = 0
          endif
      endif
C----
C Check for target parameter file if resolving stress on faults
C----
      ntrg = 0
      if (knor.eq.1.or.kshr.eq.1.or.kcou.eq.1) then
          if (ktrg.eq.0) then
              write(*,*) '!! Error: -normal, -shear, and -coul '//
     1                             'require target faults'
              call usage('!! Use -trg TRGFILE')
          elseif (ktrg.eq.1) then
              inquire(file=trgfile,EXIST=ex)
              if (.not.ex) then
                  write(*,*) '!! Error: No target fault file found'
                  call usage('!! Looked for TRGFILE: '//trgfile)
              endif
          endif
          call linecount(trgfile,ntrg)
          if (ntrg.ne.nsta.and.ntrg.ne.1) then
              call usage('!! Error: NTRG must be equal to 1 or NSTA')
          endif
      endif
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE linecount(ifile,nline)
      IMPLICIT none
      CHARACTER*40 ifile
      INTEGER nline
      open(unit=41,file=ifile,status='old')
      nline = 0
  101 read(41,*,end=102)
          nline = nline + 1
          goto 101
  102 continue
      close(41)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE readhaffile(haffile,vp,vs,dens)
      IMPLICIT none
      CHARACTER*40 haffile
      REAL*8 vp,vs,dens
      INTEGER kffm,kflt,kmag,ksta,khaf,ktrg
      COMMON /ICHECK/ kffm,kflt,kmag,ksta,khaf,ktrg
      if (khaf.eq.0) then
          call autohalfspace(vp,vs,dens)
      elseif (khaf.eq.1) then
          open(unit=11,file=haffile,status='old')
          read(11,*) vp,vs,dens
      endif
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE autohalfspace(vp,vs,dens)
      IMPLICIT none
      REAL*8 vp,vs,dens
      vp = 6800.0d0
      vs = 3926.0d0
      dens = 2800.0d0
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE readsrc(ffmfile,fltfile,magfile,evlo,evla,evdp,str,dip,
     1                   rak,dx,dy,slip,hylo,hyla,nflt)
      IMPLICIT none
      CHARACTER*40 ffmfile,fltfile,magfile
      INTEGER nflt,FMAX
      PARAMETER (FMAX=1500)
      REAL*8 evlo(FMAX),evla(FMAX),evdp(FMAX),str(FMAX),dip(FMAX),
     1       rak(FMAX),dx(FMAX),dy(FMAX),slip(FMAX),hylo,hyla
      INTEGER kffm,kflt,kmag,ksta,khaf,ktrg
      INTEGER kdsp,kstn,ksts,knor,kshr,kcou,kgmt
      COMMON /ICHECK/ kffm,kflt,kmag,ksta,khaf,ktrg
      COMMON /OCHECK/ kdsp,kstn,ksts,knor,kshr,kcou,kgmt
C----
C Read fault source files
C----
      nflt = 0
      if (kffm.eq.1) then
          call readffm(ffmfile,evlo,evla,evdp,str,dip,rak,dx,dy,slip,
     1                 hylo,hyla,nflt)
      endif
      if (kflt.eq.1) then
          call readflt(fltfile,evlo,evla,evdp,str,dip,rak,dx,dy,slip,
     1                 nflt)
      endif
      if (kmag.eq.1) then
          call readmag(magfile,evlo,evla,evdp,str,dip,rak,dx,dy,slip,
     1                 nflt)
      endif
C If FFM not here, use unweighted mean coordinates as hylo and hyla
      if (kffm.ne.1) then
          call computemean(hylo,evlo,nflt)
          call computemean(hyla,evla,nflt)
      endif
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE readffm(ffmfile,evlo,evla,evdp,str,dip,rak,dx,dy,
     1                   slip,hylo,hyla,nflt)
C----
C Read FFM in standard subfault format
C----
      IMPLICIT none
      CHARACTER*40 ffmfile,dm,dxc,dyc
      INTEGER nflt,FMAX
      PARAMETER (FMAX=1500)
      REAL*8 evlo(FMAX),evla(FMAX),evdp(FMAX),str(FMAX),dip(FMAX),
     1       rak(FMAX),dx(FMAX),dy(FMAX),slip(FMAX),hylo,hyla
      INTEGER g,nseg,ct,i,nx,ny,ptr
      REAL*8 dxr,dyr
      ct = nflt
      open (unit=31,file=ffmfile,status='old')
      read (31,*) dm,dm,dm,dm,nseg
      do 103 g = 1,nseg
          read (31,*) dm,dm,dm,dm,nx,dm,dxc,dm,ny,dm,dyc
          ptr = index(dxc,'km')
          dxc(ptr:ptr+1) = ''
          ptr = index(dyc,'km')
          dyc(ptr:ptr+1) = ''
          read (dxc,*) dxr
          read (dyc,*) dyr
          dxr = 1.0d3*dxr
          dyr = 1.0d3*dyr
          read (31,*) dm,dm,dm,dm,dm,dm,dm,dm,dm,dm,hylo,dm,hyla
          do 101 i=1,7
              read (31,*) dm
  101     continue
          do 102 i = 1,nx*ny
              read (31,*) evla(ct+i),evlo(ct+i),evdp(ct+i),slip(ct+i),
     1                    rak(ct+i),str(ct+i),dip(ct+i)
              slip(ct+i) = 1.0d-2*slip(ct+i)
              dx(ct+i) = dxr
              dy(ct+i) = dyr
              evdp(ct+i) = 1.0d3*evdp(ct+i)
  102     continue
          ct = ct + nx*ny
  103 continue
      nflt = nflt + ct
      close(31)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE readflt(fltfile,evlo,evla,evdp,str,dip,rak,dx,dy,
     1                   slip,nflt)
C----
C Read fault file in (..., slip, width, length) format
C----
      IMPLICIT none
      CHARACTER*40 fltfile
      INTEGER i,f,nflt,FMAX
      PARAMETER (FMAX=1500)
      REAL*8 evlo(FMAX),evla(FMAX),evdp(FMAX),str(FMAX),dip(FMAX),
     1       rak(FMAX),dx(FMAX),dy(FMAX),slip(FMAX)
      open (unit=32,file=fltfile,status='old')
      f = 0
  201 f = f + 1
      i = nflt + f
      read (32,*,END=202) evlo(i),evla(i),evdp(i),str(i),dip(i),rak(i),
     1                    slip(i),dy(i),dx(i)
          evdp(i) = 1.0d3*evdp(i)
          dx(i) = 1.0d3*dx(i)
          dy(i) = 1.0d3*dy(i)
          goto 201
  202 continue
      nflt = i - 1
      close(32)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE readmag(magfile,evlo,evla,evdp,str,dip,rak,dx,dy,
     1                   slip,nflt)
C----
C Read fault file in (..., mag) format
C----
      IMPLICIT none
      CHARACTER*40 magfile
      INTEGER i,f,nflt,FMAX
      PARAMETER (FMAX=1500)
      REAL*8 evlo(FMAX),evla(FMAX),evdp(FMAX),str(FMAX),dip(FMAX),
     1       rak(FMAX),dx(FMAX),dy(FMAX),slip(FMAX),mag
      REAL*8 wid,len,mu
      INTEGER typ
      mu = 4.3d10 ! shear mod hard-coded, corresponds to autohalfspace params
      open (unit=33,file=magfile,status='old')
      f = 0
  201 f = f + 1
      i = nflt + f
      read (33,*,END=202) evlo(i),evla(i),evdp(i),str(i),
     1                    dip(i),rak(i),mag
          if (45.0d0.lt.rak(i).and.rak(i).lt.135.0d0) then
              typ = 2
          elseif (-135.0d0.lt.rak(i).and.rak(i).lt.-45.0d0) then
              typ = 3
          else
              typ = 1
          endif
          call wellscoppersmith(wid,len,mag,typ)
          slip(i) = ((1.0d1**(1.5d0*(mag+10.7d0)))*1.0d-7)/
     1                                            (mu*wid*len*1.0d6)
          dx(i) = len*1.0d3
          dy(i) = wid*1.0d3
          evdp(i) = 1.0d3*evdp(i)
          goto 201
  202 continue
      nflt = i - 1
      close(33)
      RETURN
      END

C----------------------------------------------------------------------C

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

      SUBROUTINE computemean(mean,array,nval)
      IMPLICIT none
      INTEGER i,nval,FMAX
      PARAMETER (FMAX=1500)
      REAL*8 mean,array(FMAX)
      mean = 0.0d0
      do 101 i = 1,nval
          mean = mean + array(i)
  101 continue
      mean = mean/dble(nval)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE computedefm(stafile,trgfile,nsta,ntrg,
     1                       nflt,evlo,evla,evdp,str,dip,rak,dx,dy,slip,
     2                       vp,vs,dens,flttyp,
     3                       dspfile,stnfile,stsfile,norfile,shrfile,
     4                       coufile)
      IMPLICIT none
      CHARACTER*40 stafile,trgfile,dspfile,stnfile,stsfile,norfile,
     4             shrfile,coufile
      INTEGER i,nsta,ntrg,progout
      REAL*8 stlo,stla,stdp,trgstr,trgdip,trgrak,frict
      INTEGER nflt,FMAX,flttyp
      PARAMETER (FMAX=1500)
      REAL*8 evlo(FMAX),evla(FMAX),evdp(FMAX),str(FMAX),dip(FMAX),
     1       rak(FMAX),dx(FMAX),dy(FMAX),slip(FMAX)
      REAL*8 vp,vs,dens
      REAL*8 uN,uE,uZ,strain(3,3),stress(3,3),norml,shear,coul
      INTEGER kffm,kflt,kmag,ksta,khaf,ktrg
      INTEGER kdsp,kstn,ksts,knor,kshr,kcou,kgmt
      COMMON /ICHECK/ kffm,kflt,kmag,ksta,khaf,ktrg
      COMMON /OCHECK/ kdsp,kstn,ksts,knor,kshr,kcou,kgmt
      if (kdsp.eq.1) open(unit=21,file=dspfile,status='unknown')
      if (kstn.eq.1) open(unit=22,file=stnfile,status='unknown')
      if (ksts.eq.1) open(unit=23,file=stsfile,status='unknown')
      if (knor.eq.1) open(unit=24,file=norfile,status='unknown')
      if (kshr.eq.1) open(unit=25,file=shrfile,status='unknown')
      if (kcou.eq.1) open(unit=26,file=coufile,status='unknown')
      i = 0
      progout = 0
      call progbar(i,nsta,progout)
      open(unit=12,file=stafile,status='old')
      if (ktrg.eq.1) open(unit=13,file=trgfile,status='old')
      if (ntrg.eq.1.and.ktrg.eq.1) read(13,*) trgstr,trgdip,trgrak,frict
      do 101 i = 1,nsta
          read(12,*) stlo,stla,stdp
          if (ntrg.eq.nsta.and.ktrg.eq.1) read(13,*) trgstr,trgdip,
     1                                               trgrak,frict
          if (kdsp.eq.1) then
              call calcdisp(uN,uE,uZ,stlo,stla,stdp,nflt,
     1                      evlo,evla,evdp,str,dip,rak,dx,dy,slip,
     2                      vp,vs,dens,flttyp)
              write(21,9990) stlo,stla,uE,uN,uZ
          endif
          if (kstn.ge.1) then
              call calcstn(strain,stlo,stla,stdp,nflt,
     1                     evlo,evla,evdp,str,dip,rak,dx,dy,slip,
     2                     vp,vs,dens,flttyp)
              if (kstn.eq.1) then
                  write(22,9991) stlo,stla,strain(1,1),strain(2,2),
     1                          strain(3,3),strain(1,2),strain(1,3),
     2                          strain(2,3)
              endif
          endif
          if (ksts.ge.1) then
              call stn2sts(stress,strain,vp,vs,dens)
              if (ksts.eq.1) then
                  write(23,9991) stlo,stla,stress(1,1),stress(2,2),
     1                          stress(3,3),stress(1,2),stress(1,3),
     2                          stress(2,3)
              endif
          endif
          if (kshr.eq.1.or.knor.eq.1.or.kcou.eq.1) then
              call coulomb(coul,norml,shear,stress,trgstr,trgdip,trgrak,
     1                     frict)
              if (kshr.eq.1) then
                  write(25,9992) stlo,stla,shear
              endif
              if (knor.eq.1) then
                  write(24,9992) stlo,stla,norml
              endif
              if (kcou.eq.1) then
                  write(26,9992) stlo,stla,coul
              endif
          endif
          call progbar(i,nsta,progout)
  101 continue
 9990 format(2F12.4,X,3F12.4)
 9991 format(2F12.4,X,6E16.6)
 9992 format(2F12.4,X,1E16.6)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE progbar(prog,prog100,progout)
      IMPLICIT none
      INTEGER prog,prog100,progout
      CHARACTER*1 CR
      CR = char(13) ! Carriage return
      if (100*prog/prog100.ge.progout) then
          write (*,1000,advance='no') 100*prog/prog100,CR
          progout = progout + 1
      endif
      if (100*prog/prog100.ge.100) write (*,1000) 100*prog/prog100
 1000 format ('Progress: [',I3,'% Complete]',A)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE calcdisp(uNnet,uEnet,uZnet,stlo,stla,stdp,nflt,
     1                    evlo,evla,evdp,str,dip,rak,dx,dy,slip,
     2                    vp,vs,dens,flttyp)
      IMPLICIT none
      REAL*8 pi,d2r
      PARAMETER (pi=4.0d0*datan(1.0d0),d2r=pi/1.8d2)
      REAL*8 uNnet,uEnet,uZnet,ux,uy,uN,uE,uz
      REAL*8 stlo,stla,stdp,dist,az,x,y
      INTEGER f,nflt,FMAX
      PARAMETER (FMAX=1500)
      REAL*8 evlo(FMAX),evla(FMAX),evdp(FMAX),str(FMAX),dip(FMAX),
     1       rak(FMAX),slip(FMAX),dx(FMAX),dy(FMAX),area
      INTEGER flttyp
      REAL*8 vp,vs,dens
C----
C Compute NEZ displacements at (stlo,stla,stdp)
C----
      uNnet = 0.0d0
      uEnet = 0.0d0
      uZnet = 0.0d0
      do 102 f = 1,nflt
          call ddistaz(dist,az,evlo(f),evla(f),stlo,stla)
          dist = dist*6.371d6
          x = dist*( dcos(az-d2r*str(f)))
          y = dist*(-dsin(az-d2r*str(f)))
          if (flttyp.eq.0) then
              area = dx(f)*dy(f)
              call o92pt(ux,uy,uz,x,y,stdp,evdp(f),dip(f),rak(f),
     1                   area,slip(f),vp,vs,dens)
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

      SUBROUTINE calcstn(strain,stlo,stla,stdp,nflt,
     1                    evlo,evla,evdp,str,dip,rak,dx,dy,slip,
     2                    vp,vs,dens,flttyp)
      IMPLICIT none
      REAL*8 pi,d2r
      PARAMETER (pi=4.0d0*datan(1.0d0),d2r=pi/1.8d2)
      REAL*8 strain(3,3),stntmp(3,3)
      REAL*8 stlo,stla,stdp,dist,az,x,y
      INTEGER f,nflt,FMAX
      PARAMETER (FMAX=1500)
      REAL*8 evlo(FMAX),evla(FMAX),evdp(FMAX),str(FMAX),dip(FMAX),
     1       rak(FMAX),slip(FMAX),dx(FMAX),dy(FMAX),area
      INTEGER flttyp
      REAL*8 vp,vs,dens
      INTEGER i,j
C----
C Compute NEZ displacements at (stlo,stla,stdp) from FFM
C----
      do 112 i = 1,3
          do 111 j = 1,3
              strain(i,j) = 0.0d0
  111     continue
  112 continue 
      do 115 f = 1,nflt
          call ddistaz(dist,az,evlo(f),evla(f),stlo,stla)
          dist = dist*6.371d6
          x = dist*( dcos(az-d2r*str(f)))
          y = dist*(-dsin(az-d2r*str(f)))
          if (flttyp.eq.0) then
              area = dx(f)*dy(f)
              call o92ptstn(stntmp,x,y,stdp,evdp(f),dip(f),rak(f),
     1                      area,slip(f),vp,vs,dens)
          else
              call o92rectstn(stntmp,x,y,stdp,evdp(f),dip(f),rak(f),
     1                        dy(f),dx(f),slip(f),vp,vs,dens)
          endif
          call rotstrain(stntmp,str(f))
          do 114 i = 1,3
              do 113 j = 1,3
                  strain(i,j) = strain(i,j) + stntmp(i,j)
  113         continue
  114     continue 
  115 continue
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE rotstrain(e,str)
      IMPLICIT none
      REAL*8 e(3,3),str,rot(3,3),rottr(3,3),tmp(3,3)
      REAL*8 pi,d2r
      PARAMETER (pi=4.0d0*datan(1.0d0),d2r=pi/1.8d2)
      rot(1,1) = dcos(d2r*str-pi/2.0d0)
      rot(2,2) = rot(1,1)
      rot(1,2) = dsin(d2r*str-pi/2.0d0)
      rot(2,1) = -rot(1,2)
      rot(1,3) = 0.0d0
      rot(2,3) = 0.0d0
      rot(3,1) = 0.0d0
      rot(3,2) = 0.0d0
      rot(3,3) = 1.0d0
      call mattr(rottr,rot)
      call matmult(tmp,rot,e)
      call matmult(e,tmp,rottr)
      RETURN
      END
      SUBROUTINE matmult(matout,mat1,mat2)
      IMPLICIT none
      REAL*8 mat1(3,3),mat2(3,3),matout(3,3)
      INTEGER i,j,k
      DO 203 i = 1,3
          DO 202 j = 1,3
              matout(i,j) = 0.0d0
              DO 201 k = 1,3
                  matout(i,j) = matout(i,j) + mat1(i,k)*mat2(k,j)
  201         CONTINUE
  202     CONTINUE
  203 CONTINUE
      RETURN
      END
      SUBROUTINE mattr(matout,matin)
      IMPLICIT none
      REAL*8 matout(3,3),matin(3,3)
      INTEGER i,j
      DO 205 i = 1,3
          DO 204 j = 1,3
              matout(i,j) = matin(j,i)
  204     CONTINUE
  205 CONTINUE
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE stn2sts(stress,strain,vp,vs,dens)
      IMPLICIT NONE
      REAL*8 stress(3,3),strain(3,3),vp,vs,dens,lam,mu,diag
C----
C CALCULATE STRESS FROM STRAIN ASSUMING ELASTIC, ISOTROPIC MATERIAL.
C----
      mu  = dens*vs*vs
      lam = dens*vp*vp - 2.0d0*mu
      if (mu.lt.10.0e9) mu = 10.0e9
      if (lam.lt.10.0e9) lam = 10.0e9
      diag = strain(1,1) + strain(2,2) + strain(3,3)
      stress(1,1) = lam*diag + 2.0d0*mu*strain(1,1)
      stress(2,2) = lam*diag + 2.0d0*mu*strain(2,2)
      stress(3,3) = lam*diag + 2.0d0*mu*strain(3,3)
      stress(1,2) = 2.0d0*mu*strain(1,2)
      stress(1,3) = 2.0d0*mu*strain(1,3)
      stress(2,3) = 2.0d0*mu*strain(2,3)
      stress(2,1) = stress(1,2)
      stress(3,1) = stress(1,3)
      stress(3,2) = stress(2,3)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE coulomb(coul,norml,shear,stress,strin,dipin,rakin,
     1                   coeffr)
C----
C GIVEN A STRESS FIELD AND THE ORIENTATION OF A TARGET PLANE, CALCULATE
C THE NORMAL, SHEAR, AND COULOMB STRESS ON THE PLANE.
C----
      IMPLICIT none
      REAL*8 pi,d2r
      PARAMETER (pi=4.0d0*datan(1.0d0),d2r=pi/180.0d0)
      REAL*8 coul,norml,shear,stress(3,3),coeffr
      REAL*8 strin,dipin,rakin,str,dip,rak,n(3),trac(3),r(3),s(3)
      INTEGER i
      str = strin*d2r
      dip = dipin*d2r
      rak = rakin*d2r
      n(1) = dsin(dip)*dsin(str+pi/2.0d0)
      n(2) = dsin(dip)*dcos(str+pi/2.0d0)
      n(3) = dcos(dip)
      do 11 i = 1,3
          trac(i) = stress(i,1)*n(1)+stress(i,2)*n(2)+stress(i,3)*n(3)
   11 continue
      norml = 0.0d0
      do 12 i = 1,3
          norml = norml + trac(i)*n(i)
   12 continue
      r(1) =  dcos(str-pi/2.0d0)*dcos(rak)
     1                          + dsin(str-pi/2.0d0)*dsin(rak)*dcos(dip)
      r(2) = -dsin(str-pi/2.0d0)*dcos(rak)
     1                          + dcos(str-pi/2.0d0)*dsin(rak)*dcos(dip)
      r(3) =                                         dsin(rak)*dsin(dip)
      s(1) = trac(1) - norml*n(1)
      s(2) = trac(2) - norml*n(2)
      s(3) = trac(3) - norml*n(3)
      shear = s(1)*r(1) + s(2)*r(2) + s(3)*r(3)
      coul = shear + coeffr*norml
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE autodispgrid(stafile,evlo,evla,evdp,str,dip,rak,dx,dy,
     1                        slip,hylo,hyla,nflt,incr)
      IMPLICIT none
      CHARACTER*40 stafile
      INTEGER nflt,FMAX
      PARAMETER (FMAX=1500)
      REAL*8 evlo(FMAX),evla(FMAX),evdp(FMAX),str(FMAX),dip(FMAX),
     1       rak(FMAX),dx(FMAX),dy(FMAX),slip(FMAX)
      REAL*8 W,E,S,N,WMAX,EMAX,SMAX,NMAX,hylo,hyla,dmax,incr
C----
C Automatically determine receiver grid, write to STAFILE
C----
      write(*,'(A)') 'Creating receiver grid in file: '//stafile
      dmax = 750.0d0 ! km
      call getmaxlims(WMAX,EMAX,SMAX,NMAX,hylo,hyla,dmax)
      call getdisplims(W,E,S,N,WMAX,EMAX,SMAX,NMAX,hylo,hyla,
     1                 evlo,evla,evdp,str,dip,rak,dx,dy,slip,nflt)
      call writestafile(stafile,W,E,S,N,incr)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE getmaxlims(WMAX,EMAX,SMAX,NMAX,hylo,hyla,dmax)
      IMPLICIT none
      REAL*8 pi,r2d,d2r
      PARAMETER (pi=4.0d0*atan(1.0d0),r2d=180.0d0/pi,d2r=pi/180.0d0)
      REAL*8 az,dmax,WMAX,EMAX,SMAX,NMAX,hylo,hyla,junk
C----
C Determine maximum geographic limits of receiver grid to W,E,S,N
C----
      WMAX = hylo-dmax*3.6d2/(2.0d0*pi*6.371d3*dcos(hyla*d2r))
      EMAX = hylo+dmax*3.6d2/(2.0d0*pi*6.371d3*dcos(hyla*d2r))
      az = 180.0d0
      call dlola(junk,SMAX,hylo,hyla,dmax,az)
      SMAX = SMAX*r2d
      az = 0.0d0
      call dlola(junk,NMAX,hylo,hyla,dmax,az)
      NMAX = NMAX*r2d
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE getdisplims(W,E,S,N,WMAX,EMAX,SMAX,NMAX,hylo,hyla,
     1                       evlo,evla,evdp,str,dip,rak,dx,dy,slip,nflt)
      IMPLICIT none
      INTEGER i,NTHR
      PARAMETER (NTHR=3)
      REAL*8 dl(NTHR),thr(NTHR)
      REAL*8 hylo,hyla,WMAX,EMAX,SMAX,NMAX,W,E,S,N,lo,la
      REAL*8 uE,uN,uZ,stdp,vp,vs,dens
      INTEGER nflt,FMAX
      PARAMETER (FMAX=1500)
      REAL*8 evlo(FMAX),evla(FMAX),evdp(FMAX),str(FMAX),dip(FMAX),
     1       rak(FMAX),dx(FMAX),dy(FMAX),slip(FMAX)
      INTEGER dir,flttyp
C----
C Get limits of receiver grid by taking steps along longitude and
C latitudes until displacements smaller than minimum threshold value.
C Variable 'dir' determines direction: 1=N, 2=S, 3=E, 4=W
C----
      dl(1) = 1.0d0 ! Amounts to move test location (deg)
      dl(2) = 0.1d0
      dl(3) = 0.05d0
      thr(1) = 0.5d0 ! Displacement thresholds (m)
      thr(2) = 0.01d0
      thr(3) = 0.001d0
      flttyp = 1
      call autohalfspace(vp,vs,dens)
      stdp = 0.0d0
      do 103 dir = 1,4
          lo = hylo
          la = hyla
C         Take first step (by increment of dl(1))
          if (dir.eq.1) la = dnint((la+dl(1))*1.0d2)*1.0d-2
          if (dir.eq.2) la = dnint((la-dl(1))*1.0d2)*1.0d-2
          if (dir.eq.3) lo = dnint((lo+dl(1))*1.0d2)*1.0d-2
          if (dir.eq.4) lo = dnint((lo-dl(1))*1.0d2)*1.0d-2
C         Take steps of decreasing size and threshold value
          i = 1
  101     if (i.gt.NTHR) goto 102
C             Check if lon, lat exceed maximum allowed values
              if (dir.eq.1.and.la.gt.NMAX) then
                  la = NMAX
                  goto 102
              elseif (dir.eq.2.and.la.lt.SMAX) then
                  la = SMAX
                  goto 102
              elseif (dir.eq.3.and.lo.gt.EMAX) then
                  lo = EMAX
                  goto 102
              elseif (dir.eq.4.and.lo.lt.WMAX) then
                  lo = WMAX
                  goto 102
              endif
C             Calculate displacement and compare to current threshold
              call calcdisp(uN,uE,uZ,lo,la,stdp,nflt,evlo,evla,evdp,
     1                      str,dip,rak,dx,dy,slip,vp,vs,dens,
     2                      flttyp)
              uN = dabs(uN)
              uE = dabs(uE)
              uZ = dabs(uZ)
              if (uN.gt.thr(i).or.uE.gt.thr(i).or.uZ.gt.thr(i)) then
                  if (dir.eq.1) la = la + dl(i)
                  if (dir.eq.2) la = la - dl(i)
                  if (dir.eq.3) lo = lo + dl(i)
                  if (dir.eq.4) lo = lo - dl(i)
                  goto 101
              else
                  i = i + 1
                  goto 101
              endif
  102     continue
          if (dir.eq.1) N = la
          if (dir.eq.2) S = la
          if (dir.eq.3) E = lo
          if (dir.eq.4) W = lo
  103 continue
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE writestafile(stafile,W,E,S,N,incr)
      IMPLICIT none
      CHARACTER*40 stafile
      REAL*8 W,E,S,N,incr,lat,lon,stdp
      INTEGER i,j
      stdp = 0.0d0
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
      write(*,9001) W,E,S,N
      write(*,9002) incr
      write(*,9003) i*j
      write(*,9004) i
      write(*,9005) j
      write(52,9001) W,E,S,N
      write(52,9002) incr
      write(52,9003) i*j
      write(52,9004) i
      write(52,9005) j
      write(*,*) 'This information is saved in autogrid.dat'
 8001 format(3F10.2)
 9001 format('Grid limits (WESN):',4F8.2)
 9002 format('Grid spacing: ',F10.2)
 9003 format('Number of points: ',I10)
 9004 format('nlon: ',I10)
 9005 format('nlat: ',I10)
      close(51)
      close(52)
      RETURN
      END
      
C----------------------------------------------------------------------C

      SUBROUTINE gcmdln(ffmfile,fltfile,magfile,stafile,haffile,trgfile,
     1                  dspfile,stnfile,stsfile,norfile,shrfile,coufile,
     2                  gmtfile,flttyp,auto,incr)
      IMPLICIT none
      CHARACTER*40 tag
      INTEGER narg,i
      CHARACTER*40 ffmfile,fltfile,magfile,stafile,haffile,trgfile,
     1             dspfile,stnfile,stsfile,norfile,shrfile,coufile,
     2             gmtfile
      INTEGER flttyp,auto
      REAL*8 incr

      INTEGER kffm,kflt,kmag,ksta,khaf,ktrg
      INTEGER kdsp,kstn,ksts,knor,kshr,kcou,kgmt
      COMMON /ICHECK/ kffm,kflt,kmag,ksta,khaf,ktrg
      COMMON /OCHECK/ kdsp,kstn,ksts,knor,kshr,kcou,kgmt

      ffmfile = 'none'
      fltfile = 'none'
      magfile = 'none'
      stafile = 'none'
      haffile = 'none'
      trgfile = 'none'
      dspfile = 'none'
      stnfile = 'none'
      stsfile = 'none'
      norfile = 'none'
      shrfile = 'none'
      coufile = 'none'
      gmtfile = 'none'
      flttyp = 1
      auto = 0
      incr = 0.1d0

C     Input file status indicators (ICHECK block variables)
C     0=undefined; 1=defined; -1=defined, but file DNE
      kffm = 0
      kflt = 0
      kmag = 0
      ksta = 0
      khaf = 0
      ktrg = 0
C     Output file status indicators (OCHECK block variables)
C     0=unused; 1=output variable; 2=use variable, do not output
      kdsp = 0
      kstn = 0
      ksts = 0
      knor = 0
      kshr = 0
      kcou = 0
      kgmt = 0

      narg = iargc()
      if (narg.eq.0) call usage('!! Error: No command line arguments '//
     1                                     'specified')
      i = 0
  101 i = i + 1
      if (i.gt.narg) goto 102
          call getarg(i,tag)
          if (tag(1:4).eq.'-ffm') then
              kffm = 1
              i = i + 1
              call getarg(i,ffmfile)
          elseif (tag(1:4).eq.'-flt') then
              kflt = 1
              i = i + 1
              call getarg(i,fltfile)
          elseif (tag(1:4).eq.'-mag') then
              kmag = 1
              i = i + 1
              call getarg(i,magfile)
          elseif (tag(1:3).eq.'-fn') then
              flttyp = 1
          elseif (tag(1:3).eq.'-pt') then
              flttyp = 0
          elseif (tag(1:4).eq.'-sta') then
              ksta = 1
              i = i + 1
              call getarg(i,stafile)
          elseif (tag(1:4).eq.'-haf') then
              khaf = 1
              i = i + 1
              call getarg(i,haffile)
          elseif (tag(1:4).eq.'-trg') then
              ktrg = 1
              i = i + 1
              call getarg(i,trgfile)
          elseif (tag(1:5).eq.'-disp') then
              kdsp = 1
              i = i + 1
              call getarg(i,dspfile)
          elseif (tag(1:7).eq.'-strain') then
              kstn = 1
              i = i + 1
              call getarg(i,stnfile)
          elseif (tag(1:7).eq.'-stress') then
              if (kstn.eq.0) kstn = 2
              ksts = 1
              i = i + 1
              call getarg(i,stsfile)
          elseif (tag(1:7).eq.'-normal') then
              if (kstn.eq.0) kstn = 2
              if (ksts.eq.0) ksts = 2
              knor = 1
              i = i + 1
              call getarg(i,norfile)
          elseif (tag(1:6).eq.'-shear') then
              if (kstn.eq.0) kstn = 2
              if (ksts.eq.0) ksts = 2
              kshr = 1
              i = i + 1
              call getarg(i,shrfile)
          elseif (tag(1:5).eq.'-coul') then
              if (kstn.eq.0) kstn = 2
              if (ksts.eq.0) ksts = 2
              kcou = 1
              i = i + 1
              call getarg(i,coufile)
          elseif (tag(1:9).eq.'-autodisp') then
              auto = 1
              kdsp = 1
              if (ksta.eq.0) stafile = 'autosta.dat'
              ksta = 1
              i = i + 1
              if (i.gt.narg) then
                  goto 102
              else
                  call getarg(i,tag)
                  if (tag(1:1).eq.'-') then
                      i = i - 1
                      goto 101
                  endif
              endif
              read(tag,'(BN,F10.0)') incr
          elseif (tag(1:4).eq.'-gmt') then
              kgmt = 1
              i = i + 1
              call getarg(i,gmtfile)
          elseif (tag(1:2).eq.'-h'.or.tag(1:5).eq.'-help') then
              call usage(' ')
          elseif (tag(1:2).eq.'-d'.or.tag(1:8).eq.'-details') then
              call usage('long')
          else
              call usage('!! Error: no option '//tag)
          endif
          goto 101
  102 continue

      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE usage(str)
      IMPLICIT NONE
      INTEGER STO,lstr
      PARAMETER (STO=6)
      CHARACTER str*(*)
      if (str.ne.' ') then
          lstr = len(str)
          write(STO,*) str(1:lstr)
          write(STO,*)
      endif
      write(STO,*)
     1 'Usage: o92util -ffm FFMFILE -flt FLTFILE -mag MAGFILE ',
     2                        '-fn -pt'
      write(STO,*)
     1 '               -sta STAFILE -haf HAFFILE -trg TRGFILE'
      write(STO,*)
     1 '               -disp DSPFILE -strain STRAFILE -stress ',
     2                             'STREFILE '
      write(STO,*)
     1 '               -normal NORMFILE -shear SHEARFILE -coul COULFILE'
      write(STO,*)
     1 '               -autodisp INCR -gmt GMTFILE -h -? -details'
      write(STO,*)
      if (str.eq.'long') write(STO,*)
      if (str.eq.'long') write(STO,*)
     1 'A utility for computing displacements, strains, and stresses ',
     2 'in an elastic'
      if (str.eq.'long') write(STO,*)
     1 'half-space resulting from rectangular shear dislocations.'
      if (str.eq.'long') write(STO,*)
      write(STO,*)
     1 '-ffm FFMFILE            Finite fault ',
     2                          'file (standard subfault format)'
      write(STO,*)
     1 '-flt FLTFILE            Fault file in "psmeca -Sa" format ',
     2                            '(with slip and dimensions)'
      if (str.eq.'long') write(STO,*)
      if (str.eq.'long') write(STO,*)
     1 '                                           along-dip along_str'
      if (str.eq.'long') write(STO,*)
     1 '    evlo evla evdp(km) str dip rak slip(m)  wid(km)   len(km)'
      if (str.eq.'long') write(STO,*)
      write(STO,*)
     1 '-mag MAGFILE            Fault file in "psmeca -Sa" format ',
     2                                '(with magnitude)'
      if (str.eq.'long') write(STO,*)
      if (str.eq.'long') write(STO,*)
     1 '    evlo evla evdp(km) str dip rak mag'
      if (str.eq.'long') write(STO,*)
     1 '    (magnitude converted to slip, width, and length using ',
     2      'Wells and'
      if (str.eq.'long') write (STO,*)
     1 '    Coppersmith (1994) empirical relations)'
      if (str.eq.'long') write(STO,*)
      write(STO,*)
     1 '-fn          (default)  Treat subfaults as finite sources'
      write(STO,*)
     1 '-pt                     Treat subfaults as point sources'
      write(STO,*)
     1 '-sta STAFILE            Receiver location file'
      if (str.eq.'long') write(STO,*)
      if (str.eq.'long') write(STO,*)
     1 '    {stlo stla stdp(km)}'
      if (str.eq.'long') write(STO,*)
      write(STO,*)
     1 '-haf HAFFILE            Elastic half-space parameter file'
      if (str.eq.'long') write(STO,*)
      if (str.eq.'long') write(STO,*)
     1 '    If -haf not used, default to vp=6800 vs=3925 dens=2800'
      if (str.eq.'long') write(STO,*)
     1 '    {vp(m/s) vs(m/s) dens(km/m^3)}'
      if (str.eq.'long') write(STO,*)
      write(STO,*)
     1 '-trg TRGFILE            Target fault ',
     2            'file (required for -normal, -shear, -coul)'
      if (str.eq.'long') write(STO,*)
      if (str.eq.'long') write(STO,*)
     1 '    If ntrg=1: all locations in STAFILE have same ',
     2                    'target parameters'
      if (str.eq.'long') write(STO,*)
     1 '    If ntrg=nsta: target fault parameters ',
     2                     'correspond to station locations'
      if (str.eq.'long') write(STO,*)
     1 '    {str dip rak frict}'
      if (str.eq.'long') write(STO,*)
      write(STO,*)
     1 '-disp DSPFILE           Displacement (ENZ)'
      if (str.eq.'long') write(STO,*)
      if (str.eq.'long') write(STO,*)
     1 '    {stlo stla u_E(m) u_N(m) u_Z(m)}'
      if (str.eq.'long') write(STO,*)
      write(STO,*)
     1 '-strain STRNFILE        Strain matrix (EE NN ZZ EN EZ NZ)'
      if (str.eq.'long') write(STO,*)
      if (str.eq.'long') write(STO,*)
     1 '    {stlo stla e_EE e_NN e_ZZ e_EN e_EZ e_NZ}'
      if (str.eq.'long') write(STO,*)
      write(STO,*)
     1 '-stress STRSFILE        Stress matrix (EE NN ZZ EN EZ NZ)'
      if (str.eq.'long') write(STO,*)
      if (str.eq.'long') write(STO,*)
     1 '    {stlo stla s_EE(Pa) s_NN(Pa) s_ZZ(Pa) s_EN(Pa) ',
     2                         's_EZ(Pa) s_NZ(Pa)}'
      if (str.eq.'long') write(STO,*)
      write(STO,*)
     1 '-normal NORMFILE        Normal stress (',
     2                                     'resolved on TRGFILE faults)'
      if (str.eq.'long') write(STO,*)
      if (str.eq.'long') write(STO,*)
     1 '    positive = dilation, negative = compression'
      if (str.eq.'long') write(STO,*)
      write(STO,*)
     1 '-shear SHEARFILE        Shear stress (',
     2                                     'resolved on TRGFILE faults)'
      write(STO,*)
     1 '-coul COULFILE          Coulomb stress (',
     2                                     'resolved on TRGFILE faults)'
      write(STO,*)
     1 '-autodisp [INCR]        Generate receiver grid (spacing ',
     2                          'INCR optional), compute displacements'
      write(STO,*)
     1 '-gmt GMTFILE            Fault file for use with "psxy -SJ ',
     2                                   '-C<cptfile>"'
      write(STO,*)
     1 '-h/-help                Short online help'
      write(STO,*)
     1 '-d/-details             Long online help'
      write(STO,*)
c      write(STO,*)
c     1 '123456789012345678901234567890123456789012345678901234567890',
c     2 '12345678901234567890'
      STOP
      END

      
