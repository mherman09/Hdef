      PROGRAM MAIN
C----
C Utility for computing displacements, strains, and stresses in an
C elastic half-space resulting from rectangular shear dislocations.
C
C Okada, Y. (1992) Internal deformation due to shear and tensile faults
C     in a half-space. BULLETIN OF THE SEISMOLOGICAL SOCIETY OF AMERICA
C     vol. 82, no. 2, pp. 1018-1040.
C----
      IMPLICIT NONE
      CHARACTER*40 ffmf,fltf,magf,haff,staf,trgf,dspf,stnf,stsf,norf,
     1             shrf,coulf,gmtf
      INTEGER flttyp,auto,xy,prog,long
      REAL*8 incr,val
      INTEGER nflt,FMAX
      PARAMETER (FMAX=1500)
      REAL*8 evlo(FMAX),evla(FMAX),evdp(FMAX),str(FMAX),dip(FMAX),
     1       rak(FMAX),dx(FMAX),dy(FMAX),slip(FMAX),hylo,hyla
      REAL*8 vp,vs,dens
      INTEGER nsta,ntrg
      INTEGER kffm,kflt,kmag,ksta,khaf,ktrg        ! kvar definitions
      INTEGER kdsp,kstn,ksts,knor,kshr,kcou,kgmt   ! in gcmdln comments
      COMMON /ICHECK/ kffm,kflt,kmag,ksta,khaf,ktrg
      COMMON /OCHECK/ kdsp,kstn,ksts,knor,kshr,kcou,kgmt

C----
C Parse command line
C----
      call gcmdln(ffmf,fltf,magf,haff,staf,trgf,dspf,stnf,stsf,norf,
     1            shrf,coulf,gmtf,flttyp,auto,incr,val,xy,prog,long)
C      call dbggcmdln(ffmf,fltf,magf,haff,staf,trgf,dspf,stnf,stsf,norf,
C     2               shrf,coulf,gmtf,flttyp,auto,incr,val,xy,prog,
C     3               long)

C----
C Check for input fault source files, and whether output is specified
C----
      call chksrc(ffmf,fltf,magf)
      call chkout()

C----
C Read input fault files
C----
      call readsrc(ffmf,fltf,magf,nflt,evlo,evla,evdp,str,dip,rak,dx,dy,
     1             slip,hylo,hyla)

C----
C Check for auto flag to define stations
C----
       if (auto.ne.0) then
           if (auto.lt.10.or.val.gt.1d90.or.incr.gt.1d90) then
               call usage('!! Error: -auto AB VAL INCR incorrect')
           endif
           call doauto(staf,nflt,evlo,evla,evdp,str,dx,dy,slip,
     1                 auto,val,incr,xy)
       endif

C----
C Check for input receiver, half-space, and target parameter files
C If resolving stresses, verify NTRG=1 or NTRG=NSTA
C----
      call chkhafstatrg(haff,staf,trgf,nsta,ntrg)
      call readhaff(haff,vp,vs,dens)

C----
C Compute displacements, strains, and stresses
C----
      call calcdefm(staf,trgf,nsta,ntrg,nflt,evlo,evla,evdp,str,dip,rak,
     1              dx,dy,slip,flttyp,vp,vs,dens,dspf,stnf,stsf,norf,
     2              shrf,coulf,xy,prog,long)

C----
C Write fault input to GMT-compatible file
C----
      if (kgmt.eq.1) then
          call writegmt(gmtf,nflt,evlo,evla,evdp,str,dip,rak,dx,dy,
     1                  slip)
      endif

      END

C----------------------------------------------------------------------C

      SUBROUTINE chksrc(ffmf,fltf,magf)
C----
C Check (a) that fault source files were defined in gcmdln and (b) that
C these fault files exist.
C----
      IMPLICIT NONE
      CHARACTER*40 ffmf,fltf,magf
      LOGICAL ex
      INTEGER kffm,kflt,kmag,ksta,khaf,ktrg
      COMMON /ICHECK/ kffm,kflt,kmag,ksta,khaf,ktrg

C Verify input fault files were defined on command line
      if (kffm.eq.0.and.kflt.eq.0.and.kmag.eq.0) then
          write(*,*) '!! Error: No fault file specified'
          call usage('!! Use -ffm FFMFILE, -flt FLTFILE, or -mag '//
     1                  'MAGFILE')
      endif

C Look for defined files in working directory
      if (kffm.eq.1) then
          inquire(file=ffmf,EXIST=ex)
          if (.not.ex) kffm = -1
      endif
      if (kflt.eq.1) then
          inquire(file=fltf,EXIST=ex)
          if (.not.ex) kflt = -1
      endif
      if (kmag.eq.1) then
          inquire(file=magf,EXIST=ex)
          if (.not.ex) kmag = -1
      endif

C Quit to usage if none of input is found. Warn (but do not quit) if
C at least one input is found but another is not.
      if (kffm.ne.1.and.kflt.ne.1.and.kmag.ne.1) then
          write(*,*) '!! Error: No fault files found'
          if (kffm.eq.-1) write(*,*) '!! Looked for FFMFILE: '//ffmf
          if (kflt.eq.-1) write(*,*) '!! Looked for FLTFILE: '//fltf
          if (kmag.eq.-1) write(*,*) '!! Looked for MAGFILE: '//magf
          call usage('!! Check -ffm, -flt, and -mag arguments')
      else
          if (kffm.eq.-1) write(*,*) '!! Warning: No FFMFILE: '//ffmf
          if (kflt.eq.-1) write(*,*) '!! Warning: No FLTFILE: '//fltf
          if (kmag.eq.-1) write(*,*) '!! Warning: No MAGFILE: '//magf
      endif
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE chkout()
C----
C Verify that an output file is specified.
C----
      IMPLICIT NONE
      INTEGER kdsp,kstn,ksts,knor,kshr,kcou,kgmt
      COMMON /OCHECK/ kdsp,kstn,ksts,knor,kshr,kcou,kgmt
      if (kdsp.eq.0.and.kstn.eq.0) then
          write(*,*) '!! Error: No output file specified'
          call usage('!! Use -disp, -strain, -stress, '//
     1                   '-normal, -shear, -coul, or an -auto option')
      endif
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE chkhafstatrg(haff,staf,trgf,nsta,ntrg)
C----
C Check (a) that station and half-space files were defined in gcmdln,
C and (b) verify that they exist. Count NSTA. If HAFFILE is not
C specified or does not exist, use default values with autohaf.
C If resolving stress on faults, check for target parameter file and
C verify that NTRG=1 or NTRG=NSTA.
C----
      IMPLICIT NONE
      CHARACTER*40 haff,staf,trgf
      LOGICAL ex
      INTEGER nsta,ntrg,ptr,i
      INTEGER kffm,kflt,kmag,ksta,khaf,ktrg
      INTEGER kdsp,kstn,ksts,knor,kshr,kcou,kgmt
      COMMON /ICHECK/ kffm,kflt,kmag,ksta,khaf,ktrg
      COMMON /OCHECK/ kdsp,kstn,ksts,knor,kshr,kcou,kgmt
      if (ksta.eq.0) then
          write(*,*) '!! Error: No station file specified'
          call usage('!! Use -sta STAFILE or an -auto option')
      elseif (ksta.eq.1) then
          inquire(file=staf,EXIST=ex)
          if (.not.ex) then
              write(*,*) '!! Error: No station file found'
              call usage('!! Looked for STAFILE: '//staf)
          endif
      endif
      call linecount(staf,nsta)
      if (khaf.eq.0) then
C          write(*,*) '!! Using default half-space parameters'
      elseif (khaf.eq.1) then
          inquire(file=haff,EXIST=ex)
          if (.not.ex) then
              write(*,*) '!! Warning: No half-space file found'
              write(*,*) '!! Looked for HAFFILE: '//haff
              write(*,*) '!! Using default half-space parameters'
              khaf = 0
          endif
      endif
C     Check for target parameter file if resolving stress on faults
      ntrg = 0
      if (knor.eq.1.or.kshr.eq.1.or.kcou.eq.1) then
          if (ktrg.eq.0) then
              write(*,*) '!! Error: -normal, -shear, and -coul '//
     1                             'require target faults'
              call usage('!! Use -trg TRGFILE or -trg S/D/R/F')
          elseif (ktrg.eq.1) then
              inquire(file=trgf,EXIST=ex)
              if (.not.ex) then
                  write(*,*) '!! Error: No target fault file found'
                  call usage('!! Looked for TRGFILE: '//trgf)
              endif
              call linecount(trgf,ntrg)
          elseif (ktrg.eq.2) then
              do 911 i = 1,3
                  ptr = index(trgf,'/')
                  trgf(ptr:ptr) = ' '
  911         continue
              ntrg = 1
          endif
          if (ntrg.ne.nsta.and.ntrg.ne.1) then
              call usage('!! Error: NTRG must be equal to 1 or NSTA')
          endif
      endif
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE linecount(ifile,nline)
      IMPLICIT NONE
      CHARACTER*40 ifile
      INTEGER nline
      open(unit=41,file=ifile,status='old')
      nline = 0
  151 read(41,*,end=152)
          nline = nline + 1
          goto 151
  152 continue
      close(41)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE readsrc(ffmf,fltf,magf,nflt,evlo,evla,evdp,str,dip,rak,
     1                   dx,dy,slip,hylo,hyla)
C----
C Read shear dislocation parameters from fault source files. Parameter
C FMAX defines maximum number of shear dislocations that can be stored.
C----
      IMPLICIT NONE
      CHARACTER*40 ffmf,fltf,magf
      INTEGER nflt,FMAX
      PARAMETER (FMAX=1500)
      REAL*8 evlo(FMAX),evla(FMAX),evdp(FMAX),str(FMAX),dip(FMAX),
     1       rak(FMAX),dx(FMAX),dy(FMAX),slip(FMAX),hylo,hyla
      INTEGER kffm,kflt,kmag,ksta,khaf,ktrg
      INTEGER kdsp,kstn,ksts,knor,kshr,kcou,kgmt
      COMMON /ICHECK/ kffm,kflt,kmag,ksta,khaf,ktrg
      COMMON /OCHECK/ kdsp,kstn,ksts,knor,kshr,kcou,kgmt
      nflt = 0
      if (kffm.eq.1) then
          call readffm(ffmf,evlo,evla,evdp,str,dip,rak,dx,dy,slip,
     1                 hylo,hyla,nflt)
      endif
      if (kflt.eq.1) then
          call readflt(fltf,evlo,evla,evdp,str,dip,rak,dx,dy,slip,
     1                 nflt)
      endif
      if (kmag.eq.1) then
          call readmag(magf,evlo,evla,evdp,str,dip,rak,dx,dy,slip,
     1                 nflt)
      endif
C     If FFM not used, use unweighted mean coordinates as hylo and hyla
      if (kffm.ne.1) then
          call calcmean(hylo,evlo,nflt)
          call calcmean(hyla,evla,nflt)
      endif
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE readffm(ffmf,evlo,evla,evdp,str,dip,rak,dx,dy,slip,
     1                   hylo,hyla,nflt)
C----
C Read shear dislocations from FFM in standard subfault format.
C----
      IMPLICIT NONE
      CHARACTER*40 ffmf,du,dxc,dyc
      INTEGER nflt,FMAX
      PARAMETER (FMAX=1500)
      REAL*8 evlo(FMAX),evla(FMAX),evdp(FMAX),str(FMAX),dip(FMAX),
     1       rak(FMAX),dx(FMAX),dy(FMAX),slip(FMAX),hylo,hyla
      INTEGER g,nseg,ct,i,nx,ny,ptr
      REAL*8 dxr,dyr
      ct = nflt
      open (unit=31,file=ffmf,status='old')
C Header indicates number of fault segments
      read (31,*) du,du,du,du,nseg
C For each segment, read header describing subfaults
      do 313 g = 1,nseg
          read (31,*) du,du,du,du,nx,du,dxc,du,ny,du,dyc
          ptr = index(dxc,'km')
          dxc(ptr:ptr+1) = ''
          ptr = index(dyc,'km')
          dyc(ptr:ptr+1) = ''
          read (dxc,*) dxr
          read (dyc,*) dyr
          dxr = 1.0d3*dxr
          dyr = 1.0d3*dyr
          read (31,*) du,du,du,du,du,du,du,du,du,du,hylo,du,hyla
          do 311 i=1,7
              read (31,*) du
  311     continue
C Read subfault data
          do 312 i = 1,nx*ny
              read (31,*) evla(ct+i),evlo(ct+i),evdp(ct+i),slip(ct+i),
     1                    rak(ct+i),str(ct+i),dip(ct+i)
              slip(ct+i) = 1.0d-2*slip(ct+i)
              dx(ct+i) = dxr
              dy(ct+i) = dyr
              evdp(ct+i) = 1.0d3*evdp(ct+i)
  312     continue
          ct = ct + nx*ny
  313 continue
      nflt = nflt + ct
      close(31)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE readflt(fltf,evlo,evla,evdp,str,dip,rak,dx,dy,slip,
     1                   nflt)
C----
C Read shear dislocations from fault file in format:
C   evlo evla evdp(km) str dip rak slip(m) dip_wid(km) str_len(km)
C----
      IMPLICIT NONE
      CHARACTER*40 fltf
      INTEGER i,f,nflt,FMAX
      PARAMETER (FMAX=1500)
      REAL*8 evlo(FMAX),evla(FMAX),evdp(FMAX),str(FMAX),dip(FMAX),
     1       rak(FMAX),dx(FMAX),dy(FMAX),slip(FMAX)
      open (unit=32,file=fltf,status='old')
      f = 0
  321 f = f + 1
      i = nflt + f
      read (32,*,END=322) evlo(i),evla(i),evdp(i),str(i),dip(i),rak(i),
     1                    slip(i),dy(i),dx(i)
          evdp(i) = 1.0d3*evdp(i)
          dx(i) = 1.0d3*dx(i)
          dy(i) = 1.0d3*dy(i)
          goto 321
  322 continue
      nflt = i - 1
      close(32)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE readmag(magf,evlo,evla,evdp,str,dip,rak,dx,dy,
     1                   slip,nflt)
C----
C Read shear dislocations from fault file in format:
C   evlo evla evdp(km) str dip rak mag
C----
      IMPLICIT NONE
      CHARACTER*40 magf
      INTEGER i,f,nflt,FMAX
      PARAMETER (FMAX=1500)
      REAL*8 evlo(FMAX),evla(FMAX),evdp(FMAX),str(FMAX),dip(FMAX),
     1       rak(FMAX),dx(FMAX),dy(FMAX),slip(FMAX),mag
      REAL*8 wid,len,mu
      INTEGER typ
      mu = 4.3d10 ! shear mod hard-coded, corresponds to autohaf params
      open (unit=33,file=magf,status='old')
      f = 0
  331 f = f + 1
      i = nflt + f
      read (33,*,END=332) evlo(i),evla(i),evdp(i),str(i),dip(i),rak(i),
     1                    mag
          if (45.0d0.lt.rak(i).and.rak(i).lt.135.0d0) then
              typ = 2 ! Reverse fault (45<rake<135)
          elseif (-135.0d0.lt.rak(i).and.rak(i).lt.-45.0d0) then
              typ = 3 ! Normal fault (-135<rake<-45)
          else
              typ = 1 ! Strike-slip faul (otherwise)
          endif
          call wellscoppersmith(wid,len,mag,typ)
          slip(i) = ((1.0d1**(1.5d0*(mag+10.7d0)))*1.0d-7)/
     1                                            (mu*wid*len*1.0d6)
          dx(i) = len*1.0d3
          dy(i) = wid*1.0d3
          evdp(i) = 1.0d3*evdp(i)
          goto 331
  332 continue
      nflt = i - 1
      close(33)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE wellscoppersmith(wid,len,mag,typ)
C----
C Use Wells and Coppersmith (1994) empirical relations to get fault
C dimensions (in kilometers) from magnitude.
C----
      IMPLICIT NONE
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

      SUBROUTINE calcmean(mean,array,nval)
      IMPLICIT NONE
      INTEGER i,nval,FMAX
      PARAMETER (FMAX=1500)
      REAL*8 mean,array(FMAX)
      mean = 0.0d0
      do 161 i = 1,nval
          mean = mean + array(i)
  161 continue
      mean = mean/dble(nval)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE readhaff(haff,vp,vs,dens)
C----
C Read vp, vs, dens from half-space file, or use default values.
C----
      IMPLICIT NONE
      CHARACTER*40 haff,ch
      REAL*8 vp,vs,dens,lamda,mu
      INTEGER kffm,kflt,kmag,ksta,khaf,ktrg
      COMMON /ICHECK/ kffm,kflt,kmag,ksta,khaf,ktrg
      if (khaf.eq.0) then
          call autohaf(vp,vs,dens)
      elseif (khaf.eq.1) then
          open(unit=11,file=haff,status='old')
          read(11,*) ch
          rewind 11
          if (ch(1:1).eq.'L'.or.ch(1:1).eq.'l') then
              read(11,*) ch,lamda,mu
              dens = 3.0d3
              vp = dsqrt((lamda+2.0d0*mu)/dens)
              vs = dsqrt(mu/dens)
          else
              read(11,*) vp,vs,dens
          endif
      endif
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE autohaf(vp,vs,dens)
      IMPLICIT NONE
      REAL*8 vp,vs,dens
      vp = 6800.0d0
      vs = 3926.0d0
      dens = 2800.0d0
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE calcdefm(staf,trgf,nsta,ntrg,nflt,evlo,evla,evdp,str,
     1                    dip,rak,dx,dy,slip,flttyp,vp,vs,dens,dspf,
     2                    stnf,stsf,norf,shrf,coulf,xy,prog,long)
C----
C Compute displacements, strains, and stresses from shear dislocations.
C----
      IMPLICIT NONE
      CHARACTER*40 staf,trgf,dspf,stnf,stsf,norf,shrf,coulf
      INTEGER i,nsta,ntrg,progout
      REAL*8 stlo,stla,stdp,trgstr,trgdip,trgrak,frict,vp,vs,dens
      INTEGER nflt,FMAX,flttyp,xy,prog,long
      PARAMETER (FMAX=1500)
      REAL*8 evlo(FMAX),evla(FMAX),evdp(FMAX),str(FMAX),dip(FMAX),
     1       rak(FMAX),dx(FMAX),dy(FMAX),slip(FMAX)
      REAL*8 uN,uE,uZ,strain(3,3),stress(3,3),norml,shear,coul
      INTEGER kffm,kflt,kmag,ksta,khaf,ktrg
      INTEGER kdsp,kstn,ksts,knor,kshr,kcou,kgmt
      COMMON /ICHECK/ kffm,kflt,kmag,ksta,khaf,ktrg
      COMMON /OCHECK/ kdsp,kstn,ksts,knor,kshr,kcou,kgmt
C     Open files specified for output
      if (kdsp.eq.1) open(unit=21,file=dspf,status='unknown')
      if (kstn.eq.1) open(unit=22,file=stnf,status='unknown')
      if (ksts.eq.1) open(unit=23,file=stsf,status='unknown')
      if (knor.eq.1) open(unit=24,file=norf,status='unknown')
      if (kshr.eq.1) open(unit=25,file=shrf,status='unknown')
      if (kcou.eq.1) open(unit=26,file=coulf,status='unknown')
C     Initialize progress indicator
      if (prog.eq.1) then
          i = 0
          progout = 0
          call progbar(i,nsta,progout)
      endif
C     If NTRG=1, read target fault parameters here
      if (ktrg.eq.2) read(trgf,*) trgstr,trgdip,trgrak,frict
      if (ktrg.eq.1) open(unit=13,file=trgf,status='old')
      if (ntrg.eq.1.and.ktrg.eq.1) read(13,*) trgstr,trgdip,trgrak,frict
C     Compute deformation for each receiver location
      open(unit=12,file=staf,status='old')
      do 201 i = 1,nsta
          read(12,*) stlo,stla,stdp
          stdp = stdp*1.0d3
C         If NTRG=1, read target fault parameters here
          if (ntrg.eq.nsta.and.ktrg.eq.1.and.nsta.ne.1) then
              read(13,*) trgstr,trgdip,trgrak,frict
          endif
C         Compute displacements
          if (kdsp.eq.1) then
              call calcdisp(uN,uE,uZ,stlo,stla,stdp,nflt,evlo,evla,evdp,
     1                      str,dip,rak,dx,dy,slip,vp,vs,dens,flttyp,xy)
              if (long.eq.0) then
                  write(21,9990) stlo,stla,uE,uN,uZ
              else
                  write(21,9993) stlo,stla,uE,uN,uZ
              endif
          endif
C         Compute (3x3) strain matrix
          if (kstn.ge.1) then
              call calcstn(strain,stlo,stla,stdp,nflt,
     1                     evlo,evla,evdp,str,dip,rak,dx,dy,slip,
     2                     vp,vs,dens,flttyp,xy)
              if (kstn.eq.1) then
                  write(22,9991) stlo,stla,strain(1,1),strain(2,2),
     1                          strain(3,3),strain(1,2),strain(1,3),
     2                          strain(2,3)
              endif
          endif
C         Compute (3x3) stress matrix
          if (ksts.ge.1) then
              call stn2sts(stress,strain,vp,vs,dens)
              if (ksts.eq.1) then
                  write(23,9991) stlo,stla,stress(1,1),stress(2,2),
     1                          stress(3,3),stress(1,2),stress(1,3),
     2                          stress(2,3)
              endif
          endif
C         Resolve stresses onto target fault planes
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
C         Update progress indicator
          if (prog.eq.1) then
              call progbar(i,nsta,progout)
          endif
  201 continue
 9990 format(2F12.4,X,3F12.4)
 9991 format(2F12.4,X,6E16.6)
 9992 format(2F12.4,X,1E16.6)
 9993 format(2F12.4,X,3F16.8)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE progbar(prog,prog100,progout)
      IMPLICIT NONE
      INTEGER prog,prog100,progout
      CHARACTER*1 CR
      CR = char(13) ! Carriage return
      if (100*prog/prog100.ge.progout) then
          write (*,1000,advance='no') 100*prog/prog100,CR
          progout = progout + 1
      endif
      if (100*prog/prog100.ge.100) write (*,1000) 100*prog/prog100
 1000 format (' o92util progress: [',I3,'% Complete]',A)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE calcdisp(uNnet,uEnet,uZnet,stlo,stla,stdp,nflt,
     1                    evlo,evla,evdp,str,dip,rak,dx,dy,slip,
     2                    vp,vs,dens,flttyp,xy)
C----
C Compute north, east, and vertical (positive up) displacements at
C (stlo,stla,stdp) from shear dislocations.
C----
      IMPLICIT NONE
      REAL*8 pi,d2r
      PARAMETER (pi=4.0d0*datan(1.0d0),d2r=pi/1.8d2)
      REAL*8 uNnet,uEnet,uZnet,ux,uy,uN,uE,uz
      REAL*8 stlo,stla,stdp,dist,az,x,y,delx,dely
      INTEGER f,nflt,FMAX
      PARAMETER (FMAX=1500)
      REAL*8 evlo(FMAX),evla(FMAX),evdp(FMAX),str(FMAX),dip(FMAX),
     1       rak(FMAX),slip(FMAX),dx(FMAX),dy(FMAX),area
      INTEGER flttyp,xy
      REAL*8 vp,vs,dens
      uNnet = 0.0d0
      uEnet = 0.0d0
      uZnet = 0.0d0
      do 211 f = 1,nflt
          if (xy.eq.0) then
              call ddistaz(dist,az,evlo(f),evla(f),stlo,stla)
              dist = dist*6.371d6
          else
              delx = (stlo-evlo(f))*1.0d3
              dely = (stla-evla(f))*1.0d3
              dist = dsqrt(delx*delx + dely*dely)
              az = datan2(delx,dely)
          endif
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
  211 continue
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE xy2NE(uN,uE,ux,uy,str)
      IMPLICIT NONE
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
     2                    vp,vs,dens,flttyp,xy)
C----
C Compute (3x3) strain matrix (x=E, y=N, positive up) at (stlo,stla,
C stdp) from shear dislocations.
C----
      IMPLICIT NONE
      REAL*8 pi,d2r
      PARAMETER (pi=4.0d0*datan(1.0d0),d2r=pi/1.8d2)
      REAL*8 strain(3,3),stntmp(3,3)
      REAL*8 stlo,stla,stdp,dist,az,x,y,delx,dely
      INTEGER f,nflt,FMAX
      PARAMETER (FMAX=1500)
      REAL*8 evlo(FMAX),evla(FMAX),evdp(FMAX),str(FMAX),dip(FMAX),
     1       rak(FMAX),slip(FMAX),dx(FMAX),dy(FMAX),area
      INTEGER flttyp,xy
      REAL*8 vp,vs,dens
      INTEGER i,j
      do 222 i = 1,3
          do 221 j = 1,3
              strain(i,j) = 0.0d0
  221     continue
  222 continue
      do 225 f = 1,nflt
          if (xy.eq.0) then
              call ddistaz(dist,az,evlo(f),evla(f),stlo,stla)
              dist = dist*6.371d6
          else
              delx = (stlo-evlo(f))*1.0d3
              dely = (stla-evla(f))*1.0d3
              dist = dsqrt(delx*delx + dely*dely)
              az = datan2(delx,dely)
          endif
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
          do 224 i = 1,3
              do 223 j = 1,3
                  strain(i,j) = strain(i,j) + stntmp(i,j)
  223         continue
  224     continue
  225 continue
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE rotstrain(e,str)
C----
C Rotate strain matrix from x=str, y=updip horizontal, z=up to
C x=E, y=N, z=up.
C----
      IMPLICIT NONE
      REAL*8 e(3,3),str,rot(3,3),rottr(3,3),tmp(3,3),pi,d2r
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
      IMPLICIT NONE
      REAL*8 mat1(3,3),mat2(3,3),matout(3,3)
      INTEGER i,j,k
      DO 223 i = 1,3
          DO 222 j = 1,3
              matout(i,j) = 0.0d0
              DO 221 k = 1,3
                  matout(i,j) = matout(i,j) + mat1(i,k)*mat2(k,j)
  221         CONTINUE
  222     CONTINUE
  223 CONTINUE
      RETURN
      END
      SUBROUTINE mattr(matout,matin)
      IMPLICIT NONE
      REAL*8 matout(3,3),matin(3,3)
      INTEGER i,j
      DO 222 i = 1,3
          DO 221 j = 1,3
              matout(i,j) = matin(j,i)
  221     CONTINUE
  222 CONTINUE
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE stn2sts(stress,strain,vp,vs,dens)
      IMPLICIT NONE
      REAL*8 stress(3,3),strain(3,3),vp,vs,dens,lam,mu,diag
C----
C Calculate (3x3) stress matrix from (3x3) strain matrix, assuming
C isotropic, elastic material.
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

      SUBROUTINE coulomb(coul,norml,shear,stress,strin,dipin,rakin,frc)
C----
C Calculate normal, shear, and Coulomb stresses resolved onto a plane
C (given strike, dip, and rake of the plane) from the (3x3) stress
C matrix.
C----
      IMPLICIT NONE
      REAL*8 pi,d2r
      PARAMETER (pi=4.0d0*datan(1.0d0),d2r=pi/180.0d0)
      REAL*8 coul,norml,shear,stress(3,3),frc
      REAL*8 strin,dipin,rakin,str,dip,rak,n(3),trac(3),r(3),s(3)
      INTEGER i
      str = strin*d2r
      dip = dipin*d2r
      rak = rakin*d2r
      n(1) = dsin(dip)*dsin(str+pi/2.0d0)
      n(2) = dsin(dip)*dcos(str+pi/2.0d0)
      n(3) = dcos(dip)
      do 241 i = 1,3
          trac(i) = stress(i,1)*n(1)+stress(i,2)*n(2)+stress(i,3)*n(3)
  241 continue
      norml = 0.0d0
      do 242 i = 1,3
          norml = norml + trac(i)*n(i)
  242 continue
      r(1) =  dcos(str-pi/2.0d0)*dcos(rak)
     1                          + dsin(str-pi/2.0d0)*dsin(rak)*dcos(dip)
      r(2) = -dsin(str-pi/2.0d0)*dcos(rak)
     1                          + dcos(str-pi/2.0d0)*dsin(rak)*dcos(dip)
      r(3) =                                         dsin(rak)*dsin(dip)
      s(1) = trac(1) - norml*n(1)
      s(2) = trac(2) - norml*n(2)
      s(3) = trac(3) - norml*n(3)
      shear = s(1)*r(1) + s(2)*r(2) + s(3)*r(3)
      coul = shear + frc*norml
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE doauto(staf,nflt,evlo,evla,evdp,str,dx,dy,slip,
     1                  auto,val,incr,xy)
C----
C Determine the location and extent of a horizontal or vertical slice
C for computing the deformation.
C----
      IMPLICIT NONE
      REAL*8 d2k,pi,d2r
      PARAMETER (pi=4.0d0*datan(1.0d0),d2r=pi/1.8d2,
     1           d2k=6.371d3*2.0d0*pi/3.6d2)
      CHARACTER*40 staf
      INTEGER nflt,FMAX
      PARAMETER (FMAX=1500)
      REAL*8 evlo(FMAX),evla(FMAX),evdp(FMAX),str(FMAX),dx(FMAX),
     1       dy(FMAX),slip(FMAX)
      INTEGER auto,xy
      REAL*8 incr,val
      REAL*8 mom,mag,str0
      REAL*8 xmin,xmax,xincr,ymin,ymax,yincr,x0,y0,d0,siz
      call centroid(x0,y0,d0,mom,nflt,evlo,evla,evdp,dx,dy,slip)
      call strike(str0,nflt,str,dx,dy,slip)
      mag = 2.0d0/3.0d0*dlog10(mom*1.0d7) - 10.7d0
      siz = 150.0d0*2.0d0**(mag-5.0d0)
      if (siz.gt.1.5d3) siz = 1.5d3
      if (auto.lt.20) then
          if (xy.eq.0) then
              xmin = x0 - siz*0.5d0/(d2k*dcos(y0*d2r))
              xmax = x0 + siz*0.5d0/(d2k*dcos(y0*d2r))
              ymin = y0 - siz*0.5d0/d2k
              ymax = y0 + siz*0.5d0/d2k
              xincr = incr/(d2k*dcos(y0*d2r))
              yincr = incr/d2k
          else
              xmin = x0 - siz*0.5d0
              xmax = x0 + siz*0.5d0
              ymin = y0 - siz*0.5d0
              ymax = y0 + siz*0.5d0
              xincr = incr
              yincr = incr
          endif
          call hgrid(staf,xmin,xmax,xincr,ymin,ymax,yincr,val)
      elseif (auto.lt.30) then
          xmin = -siz*0.5d0
          xmax = siz*0.5d0
          ymin = 0.0d0
          ymax = d0*1.0d-3 + siz*0.5d0
          if (mod(auto,10).eq.1) then
              str0 = str0 - 90.0d0
          endif
          call shift(x0,y0,str0,val,xy)
          call vgrid(staf,xmin,xmax,ymin,ymax,incr,x0,y0,str0,xy)
      endif
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE vgrid(staf,xmin,xmax,ymin,ymax,incr,x0,y0,str0,xy)
C----
C Write a vertical grid to file staf.
C----
      IMPLICIT none
      REAL*8 pi,r2d,d2r
      PARAMETER (pi=4.0d0*atan(1.0d0),r2d=180.0d0/pi,d2r=pi/180.0d0)
      CHARACTER*40 staf
      REAL*8 xmin,xmax,ymin,ymax,incr,x0,y0,str0,dx,x1,y1,z
      INTEGER xy
      open(unit=96,file=staf,status='unknown')
      z = ymin
  961 if (z.le.ymax) then
          dx = xmin
  962     if (dx.le.xmax) then
              if (xy.eq.0) then
                  call dlola(x1,y1,x0,y0,dx,str0)
                  x1 = x1*r2d
                  y1 = y1*r2d
              else
                  x1 = x0 + dx*dsin(str0*d2r)
                  y1 = y0 + dx*dcos(str0*d2r)
              endif
              write(96,*) x1,y1,z
              dx = dx + incr
              goto 962
          endif
          z = z + incr
          goto 961
      endif
      close(96)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE hgrid(staf,xmin,xmax,xincr,ymin,ymax,yincr,dep)
C----
C Write a horizontal grid to file staf.
C----
      IMPLICIT none
      CHARACTER*40 staf
      REAL*8 xmin,xmax,xincr,ymin,ymax,yincr,dep,x,y
      open(unit=95,file=staf,status='unknown')
      x = xmin
  951 if (x.le.xmax) then
          y = ymin
  952     if (y.le.ymax) then
              write(95,*) x,y,dep
              y = y + yincr
              goto 952
          endif
          x = x + xincr
          goto 951
      endif
      close(95)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE shift(x0,y0,str0,val,xy)
C----
C Shift a point (x0,y0) by val. Variable xy controls geographic or
C Cartesian coordinates
C----
      IMPLICIT none
      REAL*8 pi,r2d,d2r
      PARAMETER (pi=4.0d0*atan(1.0d0),r2d=180.0d0/pi,d2r=pi/180.0d0)
      REAL*8 x0,y0,x1,y1,str0,shftaz,val
      INTEGER xy
      shftaz = str0 + 90.0d0
      if (xy.eq.0) then
          call dlola(x1,y1,x0,y0,val,shftaz)
          x0 = x1*r2d
          y0 = y1*r2d
      else
          x0 = x0 + val*dsin(shftaz*d2r)
          y0 = y0 + val*dcos(shftaz*d2r)
      endif
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE strike(str0,nflt,str,dx,dy,slip)
C----
C Compute mean strike weighted by moment
C----
      IMPLICIT none
      INTEGER nflt,FMAX,i
      PARAMETER (FMAX=1500)
      REAL*8 str(FMAX),dx(FMAX),dy(FMAX),slip(FMAX)
      REAL*8 str0,tmp,mom,mu
      mu = 4.3d10
      str0 = 0.0d0
      mom = 0.0d0
      do 101 i = 1,nflt
          tmp = mu*dy(i)*dx(i)*slip(i)
          str0 = str0 + tmp*str(i)
          mom = mom + tmp
  101 continue
      str0 = str0/mom
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE centroid(x0,y0,d0,mom,nflt,evlo,evla,evdp,dx,dy,slip)
C----
C Compute event centroid and total moment
C----
      IMPLICIT none
      INTEGER nflt,FMAX,i
      PARAMETER (FMAX=1500)
      REAL*8 evlo(FMAX),evla(FMAX),evdp(FMAX),dx(FMAX),dy(FMAX),
     1       slip(FMAX)
      REAL*8 x0,y0,d0,mom,mu,tmp
      mu = 4.3d10
      x0 = 0.0d0
      y0 = 0.0d0
      d0 = 0.0d0
      mom = 0.0d0
      do 101 i = 1,nflt
          tmp = mu*dy(i)*dx(i)*slip(i)
          x0 = x0 + tmp*evlo(i)
          y0 = y0 + tmp*evla(i)
          d0 = d0 + tmp*evdp(i)
          mom = mom + tmp
  101 continue
      x0 = x0/mom
      y0 = y0/mom
      d0 = d0/mom
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE writegmt(gmtf,nflt,evlo,evla,evdp,str,dip,rak,dx,dy,
     1                    slip)
C----
C Write shear dislocations into format for use with psxy -SJ -C<cptfile>
C----
      IMPLICIT none
      REAL*8 pi,d2r
      PARAMETER (pi=4.0d0*atan(1.0d0),d2r=pi/180.0d0)
      CHARACTER*40 gmtf
      INTEGER i,nflt,FMAX
      PARAMETER (FMAX=1500)
      REAL*8 evlo(FMAX),evla(FMAX),evdp(FMAX),str(FMAX),dip(FMAX),
     1       rak(FMAX),dx(FMAX),dy(FMAX),slip(FMAX),wid,len
      rak(1) = rak(1)
      evdp(1) = evdp(1)
      open(unit=51,file=gmtf,status='unknown')
      do 511 i = 1,nflt
          len = dx(i)*1.0d-3
          wid = dcos(dip(i)*d2r)*dy(i)*1.0d-3
          write(51,*) evlo(i),evla(i),slip(i),str(i),len,wid
  511 continue
      close(51)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE gcmdln(ffmf,fltf,magf,haff,staf,trgf,dspf,stnf,stsf,
     1                  norf,shrf,coulf,gmtf,flttyp,auto,incr,val,xy,
     2                  prog,long)
      IMPLICIT NONE
      CHARACTER*40 tag
      INTEGER narg,i,ptr
      CHARACTER*40 ffmf,fltf,magf,haff,staf,trgf,
     1             dspf,stnf,stsf,norf,shrf,coulf,
     2             gmtf
      INTEGER flttyp,auto,xy,prog,long
      REAL*8 incr,val
      INTEGER kffm,kflt,kmag,ksta,khaf,ktrg
      INTEGER kdsp,kstn,ksts,knor,kshr,kcou,kgmt
      COMMON /ICHECK/ kffm,kflt,kmag,ksta,khaf,ktrg
      COMMON /OCHECK/ kdsp,kstn,ksts,knor,kshr,kcou,kgmt
C Initialize variables
      ffmf = 'none'
      fltf = 'none'
      magf = 'none'
      staf = 'none'
      haff = 'none'
      trgf = 'none'
      dspf = 'none'
      stnf = 'none'
      stsf = 'none'
      norf = 'none'
      shrf = 'none'
      coulf = 'none'
      gmtf = 'none'
      flttyp = 1
      auto = 0
      xy = 0
      prog = 0
      incr = 1.0d99
      val = 1.0d99
      long = 0
C Input file status indicators (ICHECK block)
C 0=unspecified; 1=specified; -1=specified, but file not found
      kffm = 0
      kflt = 0
      kmag = 0
      ksta = 0
      khaf = 0
      ktrg = 0
C Output file status indicators (OCHECK block variables)
C 0=unused; 1=output value; 2=use value, do not output
      kdsp = 0
      kstn = 0
      ksts = 0
      knor = 0
      kshr = 0
      kcou = 0
      kgmt = 0
C Parse command line
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
              call getarg(i,ffmf)
          elseif (tag(1:4).eq.'-flt') then
              kflt = 1
              i = i + 1
              call getarg(i,fltf)
          elseif (tag(1:4).eq.'-mag') then
              kmag = 1
              i = i + 1
              call getarg(i,magf)
          elseif (tag(1:3).eq.'-fn') then
              flttyp = 1
          elseif (tag(1:3).eq.'-pt') then
              flttyp = 0
          elseif (tag(1:3).eq.'-xy') then
              xy = 1
          elseif (tag(1:5).eq.'-prog') then
              prog = 1
          elseif (tag(1:4).eq.'-sta') then
              ksta = 1
              i = i + 1
              call getarg(i,staf)
          elseif (tag(1:4).eq.'-haf') then
              khaf = 1
              i = i + 1
              call getarg(i,haff)
          elseif (tag(1:4).eq.'-trg') then
              ktrg = 1
              i = i + 1
              call getarg(i,trgf)
              ptr = index(trgf,'/') ! check format of argument
              if (ptr.ne.0) then
                  ktrg = 2
              endif
          elseif (tag(1:5).eq.'-disp') then
              kdsp = 1
              i = i + 1
              call getarg(i,dspf)
          elseif (tag(1:7).eq.'-strain') then
              kstn = 1
              i = i + 1
              call getarg(i,stnf)
          elseif (tag(1:7).eq.'-stress') then
              if (kstn.eq.0) kstn = 2
              ksts = 1
              i = i + 1
              call getarg(i,stsf)
          elseif (tag(1:7).eq.'-normal') then
              if (kstn.eq.0) kstn = 2
              if (ksts.eq.0) ksts = 2
              knor = 1
              i = i + 1
              call getarg(i,norf)
          elseif (tag(1:6).eq.'-shear') then
              if (kstn.eq.0) kstn = 2
              if (ksts.eq.0) ksts = 2
              kshr = 1
              i = i + 1
              call getarg(i,shrf)
          elseif (tag(1:5).eq.'-coul') then
              if (kstn.eq.0) kstn = 2
              if (ksts.eq.0) ksts = 2
              kcou = 1
              i = i + 1
              call getarg(i,coulf)
          elseif (tag(1:5).eq.'-long') then
              long = 1
          elseif (tag(1:5).eq.'-auto') then
              auto = -1
              prog = 1
              if (ksta.eq.0) staf = 'autosta.dat'
              ksta = 1
              i = i + 1
              call getarg(i,tag)
              read(tag,'(BN,I3)') auto
              i = i + 1
              call getarg(i,tag)
              read(tag,'(BN,F10.0)') val
              i = i + 1
              call getarg(i,tag)
              read(tag,'(BN,F10.0)') incr
          elseif (tag(1:4).eq.'-gmt') then
              kgmt = 1
              i = i + 1
              call getarg(i,gmtf)
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
      INTEGER lstr
      CHARACTER str*(*)
      if (str.eq.'long') goto 100
      if (str.ne.' ') then
          lstr = len(str)
          write(*,*) str(1:lstr)
          write(*,*)
      endif
  100 write(*,*)
     1 'Usage: o92util -ffm FFMFILE -flt FLTFILE -mag MAGFILE ',
     2                        '[-fn|-pt]'
      write(*,*)
     1 '               -sta STAFILE -trg TRGFILE|S/D/R/F ',
     2                               '[-haf HAFFILE] [-xy]'
      write(*,*)
     1 '               -disp DSPFILE -strain STNFILE -stress ',
     2                             'STSFILE'
      write(*,*)
     1 '               -normal NORFILE -shear SHRFILE -coul COULFILE'
      write(*,*)
     1 '               [-auto AB VALUE INCR] [-long] [-prog] ',
     2                   '[-gmt GMTFILE]'
      write(*,*)
     1 '               [-h|-help] [-d|-details]'
      write(*,*)
      if (str.eq.'long') then
          write(*,*)
     1    'A utility for computing displacements, strains, and ',
     2    'stresses in an elastic'
          write(*,*)
     1    'half-space resulting from rectangular shear dislocations ',
     2    'using the analytical'
          write(*,*)
     1    'solutions from Okada (1992).'
          write(*,*)
          write(*,*) 'DESCRIPTION OF COMMAND LINE OPTIONS'
          write(*,9998)
          write(*,*)
      endif
      write(*,*)
     1 '-ffm FFMFILE            Finite fault ',
     2                          'file in standard subfault format'
      if (str.eq.'long') then
          write(*,*)
          write(*,*)
     1 '    Example of subfault format can be seen at:'
          write(*,*)
     1 '    http://comcat.cr.usgs.gov/product/finite-fault/usc000nzvd/',
     2      'us/1397258114263/'
          write(*,*)
     1 '        web/static2_out'
          write(*,*)
     1 '     NB: o92util assumes FFM displacements are in centimeters'
          write(*,*)
          write(*,9999)
          write(*,*)
      endif
      write(*,*)
     1 '-flt FLTFILE            Fault file with slip and dimensions ',
     2                       '(...slip wid len)'
      if (str.eq.'long') then
          write(*,*)
          write(*,*)
     1 '    Format: evlo evla evdp(km) str dip rak slip(m)  wid(km)   ',
     2                                                       'len(km)'
          write(*,*)
     1 '                                                   along-dip ',
     2                                                    'along-str'
          write(*,*)
          write(*,9999)
          write(*,*)
      endif
      write(*,*)
     1 '-mag MAGFILE            Fault file in "psmeca -Sa" format ',
     2                                '(...mag)'
      if (str.eq.'long') then
          write(*,*)
          write(*,*)
     1 '    Format: evlo evla evdp(km) str dip rak mag'
          write(*,*)
     1 '    The magnitude is converted to slip, width, and length ',
     2      'using Wells and'
          write (*,*)
     1 '    Coppersmith (1994) empirical relations, and a hard-coded ',
     2      'shear modulus'
          write (*,*)
     1 '    of 43 GPa.'
          write(*,*)
          write(*,9999)
          write(*,*)
      endif
      write(*,*)
     1 '-fn|-pt                 Treat faults as finite rectangular ',
     2                          '(default) or point'
      if (str.eq.'long') then
          write(*,*)
          write(*,9999)
          write(*,*)
      endif
      write(*,*)
     1 '-sta STAFILE            Receiver location file'
      if (str.eq.'long') then
          write(*,*)
          write(*,*)
     1 '    Format: stlo stla stdp(km)'
          write(*,*)
          write(*,9999)
          write(*,*)
      endif
      write(*,*)
     1 '-trg TRGFILE|S/D/R/F    Target fault kinematics (in ',
     2                              'file or on command line)'
      if (str.eq.'long') then
          write(*,*)
          write(*,*)
     1 '    Format: str dip rak frict (file TRGFILE)'
          write(*,*)
     1 '    If NTRG=1: all locations in STAFILE have same ',
     2                    'target parameters'
          write(*,*)
     1 '    If NTRG=NSTA: target fault parameters ',
     2                     'correspond to station locations'
          write(*,*)
          write(*,*)
     1 '    Format: str/dip/rak/frict (command line)'
          write(*,*)
     1 '      NTRG=1 automatically'
          write(*,*)
          write(*,9999)
          write(*,*)
      endif
      write(*,*)
     1 '-haf HAFFILE            Elastic half-space file'
      if (str.eq.'long') then
          write(*,*)
          write(*,*)
     1 '    Format: vp(m/s) vs(m/s) dens(kg/m^3)'
          write(*,*)
     1 '            Lame  lamda(Pa)  mu(Pa)'
          write(*,*)
     1 '    If velocity is given, Lame parameters computed:'
          write(*,*)
     1 '        mu = vs*vs*dens'
          write(*,*)
     1 '        lamda = vp*vp*dens - 2*mu'
          write(*,*)
     1 '    If "Lame" is first entry, Lame parameters are input'
          write(*,*)
     1 '    If -haf option not specified, o92util defaults to:'
          write(*,*)
     1 '        vp = 6800 m/s; vs = 3926 m/s; dens = 2800 kg/m^3'
          write(*,*)
          write(*,9999)
          write(*,*)
      endif
      write(*,*)
     1 '-xy                     Treat geographic ',
     2                               'coordinates as ',
     2                          'Cartesian x-y (km)'
      if (str.eq.'long') then
          write(*,*)
          write(*,*)
     1 '    Geographic coordinates evlo, evla, stlo, and stla are ',
     2      'converted to x-y (in km)'
          write(*,*)
     1 '    evlo -> evx; evla -> evy; stlo -> stx; stla -> sty'
          write(*,*)
     1 '    All other input formats and units are unchanged'
          write(*,*)
     1 '    NB: strike angle still defined from N=y-direction'
          write(*,*)
          write(*,9999)
          write(*,*)
      endif
      write(*,*)
     1 '-disp DSPFILE           Displacement (E N Z)'
      if (str.eq.'long') then
          write(*,*)
          write(*,*)
     1 '    Format: stlo stla u_E(m) u_N(m) u_Z(m)'
          write(*,*)
     1 '    Vertical displacement is positive up.'
          write(*,*)
          write(*,9999)
          write(*,*)
      endif
      write(*,*)
     1 '-strain STRNFILE        Strain matrix (EE NN ZZ EN EZ NZ)'
      if (str.eq.'long') then
          write(*,*)
          write(*,*)
     1 '    Format: stlo stla e_EE e_NN e_ZZ e_EN e_EZ e_NZ'
          write(*,*)
     1 '    The strain matrix is symmetric; the output contains ',
     2       'the six independent'
          write(*,*)
     1 '    components. Strains are unitless.'
          write(*,*)
          write(*,9999)
          write(*,*)
      endif
      write(*,*)
     1 '-stress STRSFILE        Stress matrix (EE NN ZZ EN EZ NZ)'
      if (str.eq.'long') then
          write(*,*)
          write(*,*)
     1 '    Format: stlo stla s_EE(Pa) s_NN(Pa) s_ZZ(Pa) s_EN(Pa) ',
     2                         's_EZ(Pa) s_NZ(Pa)'
          write(*,*)
     1 '    The stress matrix is symmetric; the output contains ',
     2       'the six independent'
          write(*,*)
     1 '    components. In an isotropic, elastic material (summation ',
     2                 'convention):'
          write(*,*)
     1 '        stress(i,j) = shear_mod*strain(i,j) + lambda*',
     2             'strain(i,i)*delta(i,j)'
          write(*,*)
          write(*,9999)
          write(*,*)
      endif
      write(*,*)
     1 '-normal NORMFILE        Normal stress (',
     2                              'requires target fault kinematics)'
      if (str.eq.'long') then
          write(*,*)
          write(*,*)
     1 '    Format: stlo stla normal(Pa)'
          write(*,*)
     1 '    The convention here is positive = dilation, negative = ',
     2                                    'compression.'
          write(*,*)
     1 '    This differs from many geological applications, where ',
     2          'positive normal'
          write(*,*)
     1 '    stress is compressive.'
          write(*,*)
          write(*,9999)
          write(*,*)
      endif
      write(*,*)
     1 '-shear SHEARFILE        Shear stress (',
     2                               'requires target fault kinematics)'
      if (str.eq.'long') then
          write(*,*)
          write(*,*)
     1 '    Format: stlo stla shear(Pa)'
          write(*,*)
     1 '    The reported shear stress is the projection of the shear ',
     2                      'traction'
          write(*,*)
     1 '    onto the rake vector of the target fault.'
          write(*,*)
          write(*,9999)
          write(*,*)
      endif
      write(*,*)
     1 '-coul COULFILE          Coulomb stress (',
     2                            'requires target fault kinematics)' 
      if (str.eq.'long') then
          write(*,*)
          write(*,*)
     1 '    Format: stlo stla coulomb(Pa)'
          write(*,*)
     1 '    coulomb = shear + frict*normal'
          write(*,*)
          write(*,9999)
          write(*,*)
      endif
      write(*,*)
     1 '-auto AB VALUE INCR     Automatically compute displacements ',
     2                             'or stresses'
      write(*,*)
     1 '                            A=1: Horizontal slice'
      write(*,*)
     1 '                            A=2: Vertical slice'
      write(*,*)
     1 '                            B=1: slice perpendicular to strike'
      write(*,*)
     1 '                            B=2: slice parallel to strike'
      write(*,*)
     1 '                            VALUE: A: depth (km) or B: ',
     2                                  'horizontal offset (km)'
      write(*,*)
     1 '                            INCR: grid increment (km)'
      if (str.eq.'long') then
          write(*,*)
          write(*,9999)
          write(*,*)
      endif
      write(*,*)
     1 '-long                   Increase displacement output precision'
      if (str.eq.'long') then
          write(*,*)
          write(*,*)
     1 '    Default format: (F12.4)    Long format: (F16.8)'
          write(*,*)
          write(*,9999)
          write(*,*)
      endif
      write(*,*)
     1 '-prog                   Display a progress indicator'
      if (str.eq.'long') then
          write(*,*)
          write(*,9999)
          write(*,*)
      endif
C      write(*,*)
C     1 '-incr INCR              Change grid spacing for -auto options ',
C     2                               '(default: 0.1)'
C      if (str.eq.'long') then
C          write(*,*)
C          write(*,9999)
C          write(*,*)
C      endif
C      write(*,*)
C     1 '-autodisp               Determine receiver grid, compute ',
C     2                          'displacements'
C      if (str.eq.'long') then
C          write(*,*)
C          write(*,*) '    Compute displacements on depth slice'
C          write(*,*)
C     1 '    Determines maximum extent of receiver grid to W, E, S, N, ',
C     2      'then writes'
C          write(*,*)
C     1 '    grid to STAFILE (if unspecified, default to autosta.dat). ',
C     2      'Maximum distance'
C          write(*,*)
C     1 '    in each direction is 750 km from hypocenter, i.e. 1500km ',
C     2             'x 1500km.'
C          write(*,*)
C          write(*,9999)
C          write(*,*)
C      endif
C      write(*,*)
C     1 '-autocoul OPT           Determine receiver grid, compute ',
C     2                         'Coulomb stresses'
C      write(*,*)
C     1 '                          OPT=1: Create a receiver grid at ',
C     2                                         'fixed depth'
C      write(*,*)
C     1 '                          OPT=2: (IN DEVELOPMENT, DOES NOT ',
C     2                            'WORK!)'
C      write(*,*)
C     1 '                            Create a receiver grid on ',
C     2                                    'finite fault (requires -ffm)'
C      if (str.eq.'long') then
C          write(*,*)
C          write(*,*) '    OPT=1: Compute Coulomb stress on depth slice'
C          write(*,*)
C     1 '    Determines maximum extent of receiver grid to W, E, S, N, ',
C     2      'then writes'
C          write(*,*)
C     1 '    grid to STAFILE (if unspecified, default to autosta.dat). ',
C     2      'Maximum distance'
C          write(*,*)
C     1 '    in each direction is 750 km from hypocenter, i.e. 1500km ',
C     2             'x 1500km.'
C          write(*,*)
C          write(*,9999)
C          write(*,*)
C      endif
C      write(*,*)
C     1 '-dep DEPTH              Change receiver depth for -auto ',
C     2                          'options (default: 0.0 km)'
C      if (str.eq.'long') then
C          write(*,*)
C          write(*,9999)
C          write(*,*)
C      endif
      write(*,*)
     1 '-gmt GMTFILE            Fault file for use with "psxy -SJ ',
     2                                   '-C<cptfile>"'
      if (str.eq.'long') then
          write(*,*)
          write(*,*)
     1 '    Format: stlo stla slip(m) str hor_len(km) hor_wid(km)'
          write(*,*)
     1 '    The horizontal length is the along-strike length ',
     2      'from the fault'
          write(*,*)
     1 '    definition. The horizontal width is the along-dip ',
     2      'width times cos(dip).'
          write(*,*)
          write(*,9999)
          write(*,*)
      endif
      write(*,*)
     1 '-h|-help                Short online help'
      write(*,*)
     1 '-d|-details             Long online help, show file formats ',
     2                                    'and units'
      write(*,*)
 9998 format('=======================================================',
     1       '=========================')
 9999 format('-------------------------------------------------------',
     1       '-------------------------')
c      write(*,*)
c     1 '123456789012345678901234567890123456789012345678901234567890',
c     2 '12345678901234567890'
      STOP
      END


      SUBROUTINE dbggcmdln(ffmf,fltf,magf,haff,staf,trgf,dspf,stnf,stsf,
     1                     norf,shrf,coulf,gmtf,flttyp,auto,incr,val,
     2                     xy,prog,long)
      IMPLICIT none
      CHARACTER*40 ffmf,fltf,magf,haff,staf,trgf,
     1             dspf,stnf,stsf,norf,shrf,coulf,
     2             gmtf
      INTEGER flttyp,auto,xy,prog,long
      REAL*8 incr,val
      INTEGER kffm,kflt,kmag,ksta,khaf,ktrg
      INTEGER kdsp,kstn,ksts,knor,kshr,kcou,kgmt
      COMMON /ICHECK/ kffm,kflt,kmag,ksta,khaf,ktrg
      COMMON /OCHECK/ kdsp,kstn,ksts,knor,kshr,kcou,kgmt
      print *,'DEBUGGING SUBROUTINE GCMDLN'
      print *,'---------------------------'
      print *,'KFFM',kffm,ffmf
      print *,'KFLT',kflt,fltf
      print *,'KMAG',kmag,magf
      print *,'KSTA',ksta,staf
      print *,'KHAF',khaf,haff
      print *,'KTRG',ktrg,trgf
      print *,'KDSP',kdsp,dspf
      print *,'KSTN',kstn,stnf
      print *,'KSTS',ksts,stsf
      print *,'KNOR',knor,norf
      print *,'KSHR',kshr,shrf
      print *,'KCOU',kcou,coulf
      print *,'KGMT',kgmt,gmtf
      print *,'FLTTYP',flttyp
      print *,'AUTO',auto
      print *,'INCR',incr
      print *,'VAL',val
      print *,'XY?',xy
      print *,'PROG',prog
      print *,'LONG',long
      RETURN
      END

C||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||C
CVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVC
C THESE ARE ALL THE OLD SUBROUTINES FOR AUTOMATIC OPTIONS
C
C      SUBROUTINE autodispgrid(staf,evlo,evla,evdp,str,dip,rak,dx,dy,
C     1                        slip,hylo,hyla,nflt,incr,val,dspf,xy)
CC----
CC Automatically determine receiver grid for displacement computation,
CC write grid to STAFILE.
CC----
C      IMPLICIT NONE
C      CHARACTER*40 staf,dspf
C      INTEGER nflt,FMAX,xy
C      PARAMETER (FMAX=1500)
C      REAL*8 evlo(FMAX),evla(FMAX),evdp(FMAX),str(FMAX),dip(FMAX),
C     1       rak(FMAX),dx(FMAX),dy(FMAX),slip(FMAX)
C      REAL*8 W,E,S,N,WMAX,EMAX,SMAX,NMAX,hylo,hyla,dmax,incr,val
C      write(*,'(A)') ' Creating receiver grid in file: '//staf
C      dmax = 750.0d0 ! km
C      call getmaxlims(WMAX,EMAX,SMAX,NMAX,hylo,hyla,dmax,xy)
C      call getautodisplims(W,E,S,N,WMAX,EMAX,SMAX,NMAX,hylo,hyla,
C     1                     evlo,evla,evdp,str,dip,rak,dx,dy,slip,
C     2                     nflt,xy,val)
C      call writestaf(staf,W,E,S,N,incr,val)
C      write(*,'(A)') ' Displacements written to file: '//dspf
C      RETURN
C      END
C
C
C      SUBROUTINE autocoulgrid(staf,evlo,evla,evdp,str,dip,rak,dx,dy,
C     1                        slip,hylo,hyla,nflt,incr,val,coulf,xy)
CC----
CC Automatically determine receiver grid for Coulomb stress computation,
CC write grid to COUFILE.
CC----
C      IMPLICIT NONE
C      CHARACTER*40 staf,coulf
C      INTEGER nflt,FMAX,xy
C      PARAMETER (FMAX=1500)
C      REAL*8 evlo(FMAX),evla(FMAX),evdp(FMAX),str(FMAX),dip(FMAX),
C     1       rak(FMAX),dx(FMAX),dy(FMAX),slip(FMAX)
C      REAL*8 W,E,S,N,WMAX,EMAX,SMAX,NMAX,hylo,hyla,dmax,incr,val
C      write(*,'(A)') ' Creating receiver grid in file: '//staf
C      dmax = 750.0d0 ! km
C      call getmaxlims(WMAX,EMAX,SMAX,NMAX,hylo,hyla,dmax,xy)
C      call getautostnlims(W,E,S,N,WMAX,EMAX,SMAX,NMAX,hylo,hyla,
C     1                    evlo,evla,evdp,str,dip,rak,dx,dy,slip,
C     2                    nflt,xy,val)
C      call writestaf(staf,W,E,S,N,incr,val)
C      write(*,'(A)') ' Coulomb stress written to file: '//coulf
C      RETURN
C      END
C
C
C      SUBROUTINE getmaxlims(WMAX,EMAX,SMAX,NMAX,hylo,hyla,dmax,xy)
C      IMPLICIT NONE
C      REAL*8 pi,r2d,d2r
C      PARAMETER (pi=4.0d0*atan(1.0d0),r2d=180.0d0/pi,d2r=pi/180.0d0)
C      REAL*8 az,dmax,WMAX,EMAX,SMAX,NMAX,hylo,hyla,junk
C      INTEGER xy
CC----
CC Determine maximum geographic limits of receiver grid to W,E,S,N.
CC----
C      if (xy.eq.0) then
C          WMAX = hylo-dmax*3.6d2/(2.0d0*pi*6.371d3*dcos(hyla*d2r))
C          EMAX = hylo+dmax*3.6d2/(2.0d0*pi*6.371d3*dcos(hyla*d2r))
C          az = 180.0d0
C          call dlola(junk,SMAX,hylo,hyla,dmax,az)
C          SMAX = SMAX*r2d
C          az = 0.0d0
C          call dlola(junk,NMAX,hylo,hyla,dmax,az)
C          NMAX = NMAX*r2d
C      else
C          WMAX = hylo - dmax
C          EMAX = hylo + dmax
C          SMAX = hyla - dmax
C          NMAX = hyla + dmax
C      endif
C      RETURN
C      END
C
C
C      SUBROUTINE getautodisplims(W,E,S,N,WMAX,EMAX,SMAX,NMAX,hylo,hyla,
C     1                           evlo,evla,evdp,str,dip,rak,dx,dy,slip,
C     2                           nflt,xy,val)
CC----
CC Get limits of receiver grid by taking steps along longitude and
CC latitudes until displacements smaller than minimum threshold value.
CC Variable 'dir' determines direction: 1=N, 2=S, 3=E, 4=W
CC----
C      IMPLICIT NONE
C      INTEGER i,NTHR
C      PARAMETER (NTHR=3)
C      REAL*8 dl(NTHR),thr(NTHR)
C      REAL*8 hylo,hyla,WMAX,EMAX,SMAX,NMAX,W,E,S,N,lo,la
C      REAL*8 uE,uN,uZ,stdp,vp,vs,dens,val
C      INTEGER nflt,FMAX,xy
C      PARAMETER (FMAX=1500)
C      REAL*8 evlo(FMAX),evla(FMAX),evdp(FMAX),str(FMAX),dip(FMAX),
C     1       rak(FMAX),dx(FMAX),dy(FMAX),slip(FMAX)
C      INTEGER dir,flttyp
C      if (xy.eq.0) then
C          dl(1) = 1.0d0 ! Amounts to move test location (deg)
C          dl(2) = 0.1d0
C          dl(3) = 0.05d0
C      else
C          dl(1) = 100.0d0 ! Amounts to move test location (km)
C          dl(2) = 10.0d0
C          dl(3) = 5.0d0
C      endif
C      thr(1) = 0.5d0 ! Displacement thresholds (m)
C      thr(2) = 0.01d0
C      thr(3) = 0.001d0
C      flttyp = 1
C      call autohaf(vp,vs,dens)
C      stdp = val
C      do 903 dir = 1,4
C          lo = hylo
C          la = hyla
C         Take first step (by increment of dl(1))
C          if (dir.eq.1) la = dnint((la+dl(1))*1.0d2)*1.0d-2
C          if (dir.eq.2) la = dnint((la-dl(1))*1.0d2)*1.0d-2
C          if (dir.eq.3) lo = dnint((lo+dl(1))*1.0d2)*1.0d-2
C          if (dir.eq.4) lo = dnint((lo-dl(1))*1.0d2)*1.0d-2
C         Take steps of decreasing size and threshold value
C          i = 1
C  901     if (i.gt.NTHR) goto 902
C             Check if lon, lat exceed maximum allowed values
C              if (dir.eq.1.and.la.gt.NMAX) then
C                  la = NMAX
C                  goto 902
C              elseif (dir.eq.2.and.la.lt.SMAX) then
C                  la = SMAX
C                  goto 902
C              elseif (dir.eq.3.and.lo.gt.EMAX) then
C                  lo = EMAX
C                  goto 902
C              elseif (dir.eq.4.and.lo.lt.WMAX) then
C                  lo = WMAX
C                  goto 902
C              endif
C             Calculate displacement and compare to current threshold
C              call calcdisp(uN,uE,uZ,lo,la,stdp,nflt,evlo,evla,evdp,
C     1                      str,dip,rak,dx,dy,slip,vp,vs,dens,
C     2                      flttyp,xy)
C              uN = dabs(uN)
C              uE = dabs(uE)
C              uZ = dabs(uZ)
C              if (uN.gt.thr(i).or.uE.gt.thr(i).or.uZ.gt.thr(i)) then
C                  if (dir.eq.1) la = la + dl(i)
C                  if (dir.eq.2) la = la - dl(i)
C                  if (dir.eq.3) lo = lo + dl(i)
C                  if (dir.eq.4) lo = lo - dl(i)
C                  goto 901
C              else
C                  i = i + 1
C                  goto 901
C              endif
C  902     continue
C          if (dir.eq.1) N = la
C          if (dir.eq.2) S = la
C          if (dir.eq.3) E = lo
C          if (dir.eq.4) W = lo
C  903 continue
C      RETURN
C      END
C
C
C      SUBROUTINE getautostnlims(W,E,S,N,WMAX,EMAX,SMAX,NMAX,hylo,hyla,
C     1                          evlo,evla,evdp,str,dip,rak,dx,dy,slip,
C     2                          nflt,xy,val)
CC----
CC Get limits of receiver grid by taking steps along longitude and
CC latitudes until strains smaller than minimum threshold value.
CC Variable 'dir' determines direction: 1=N, 2=S, 3=E, 4=W
CC----
C      IMPLICIT NONE
C      INTEGER i,NTHR
C      PARAMETER (NTHR=3)
C      REAL*8 dl(NTHR),thr(NTHR)
C      REAL*8 hylo,hyla,WMAX,EMAX,SMAX,NMAX,W,E,S,N,lo,la
C      REAL*8 stn(3,3),stdp,vp,vs,dens,stnmax,stnmin,val
C      INTEGER nflt,FMAX
C      PARAMETER (FMAX=1500)
C      REAL*8 evlo(FMAX),evla(FMAX),evdp(FMAX),str(FMAX),dip(FMAX),
C     1       rak(FMAX),dx(FMAX),dy(FMAX),slip(FMAX)
C      INTEGER dir,flttyp,xy
C      if (xy.eq.0) then
C          dl(1) = 1.0d0 ! Amounts to move test location (deg)
C          dl(2) = 0.1d0
C          dl(3) = 0.05d0
C      else
C          dl(1) = 100.0d0 ! Amounts to move test location (km)
C          dl(2) = 10.0d0
C          dl(3) = 5.0d0
C      endif
C      thr(1) = 1.0d-7 ! Strain thresholds (m)
C      thr(2) = 5.0d-8
C      thr(3) = 1.0d-8
C      flttyp = 1
C      call autohaf(vp,vs,dens)
C      stdp = val
C      do 913 dir = 1,4
C          lo = hylo
C          la = hyla
C         Take first step (by increment of dl(1))
C          if (dir.eq.1) la = dnint((la+dl(1))*1.0d2)*1.0d-2
C          if (dir.eq.2) la = dnint((la-dl(1))*1.0d2)*1.0d-2
C          if (dir.eq.3) lo = dnint((lo+dl(1))*1.0d2)*1.0d-2
C          if (dir.eq.4) lo = dnint((lo-dl(1))*1.0d2)*1.0d-2
C         Take steps of decreasing size and threshold value
C          i = 1
C  911     if (i.gt.NTHR) goto 912
C             Check if lon, lat exceed maximum allowed values
C              if (dir.eq.1.and.la.gt.NMAX) then
C                  la = NMAX
C                  goto 912
C              elseif (dir.eq.2.and.la.lt.SMAX) then
C                  la = SMAX
C                  goto 912
C              elseif (dir.eq.3.and.lo.gt.EMAX) then
C                  lo = EMAX
C                  goto 912
C              elseif (dir.eq.4.and.lo.lt.WMAX) then
C                  lo = WMAX
C                  goto 912
C              endif
C             Calculate strain and compare to current threshold
C              call calcstn(stn,lo,la,stdp,nflt,evlo,evla,evdp,
C     1                     str,dip,rak,dx,dy,slip,vp,vs,dens,
C     2                     flttyp,xy)
C              stnmax = dmax1(stn(1,1),stn(2,2),stn(3,3),stn(1,2),
C     1                       stn(1,3),stn(2,3))
C              stnmin = dmin1(stn(1,1),stn(2,2),stn(3,3),stn(1,2),
C     1                       stn(1,3),stn(2,3))
C              stnmax = dmax1(dabs(stnmax),dabs(stnmin))
C              if (stnmax.gt.thr(i)) then
C                  if (dir.eq.1) la = la + dl(i)
C                  if (dir.eq.2) la = la - dl(i)
C                  if (dir.eq.3) lo = lo + dl(i)
C                  if (dir.eq.4) lo = lo - dl(i)
C                  goto 911
C              else
C                  i = i + 1
C                  goto 911
C              endif
C  912     continue
C          if (dir.eq.1) N = la
C          if (dir.eq.2) S = la
C          if (dir.eq.3) E = lo
C          if (dir.eq.4) W = lo
C  913 continue
C      RETURN
C     END
C
C
C      SUBROUTINE writestaf(staf,W,E,S,N,incr,val)
CC----
CC Write receiver grid to STAFILE.
CC----
C      IMPLICIT NONE
C      CHARACTER*40 staf
C      REAL*8 W,E,S,N,incr,lat,lon,stdp,val
C      INTEGER i,j
C      stdp = val
C      open(unit=95,file=staf,status='unknown')
C      i = 0
C      lon = W
C  951 if (lon.gt.E) goto 954
C          i = i + 1
C          lat = S
C          j = 0
C  952     if (lat.gt.N) goto 953
C              write(95,8001) lon,lat,stdp
C              lat = lat + incr
C              j = j + 1
C              goto 952
C  953     continue
C          lon = lon + incr
C          goto 951
C  954 continue
C      open(unit=96,file='autogrid.dat',status='unknown')
C      write(*,9001) W,E,S,N
C      write(*,9002) incr
C      write(*,9003) i*j
C      write(*,9004) i
C      write(*,9005) j
C      write(96,9001) W,E,S,N
C      write(96,9002) incr
C      write(96,9003) i*j
C      write(96,9004) i
C      write(96,9005) j
C      write(*,*) 'This information is saved in autogrid.dat'
C 8001 format(3F10.2)
C 9001 format(' Grid limits (WESN):',4F8.2)
C 9002 format(' Grid spacing: ',F10.2)
C 9003 format(' Number of points: ',I10)
C 9004 format(' nlon: ',I10)
C 9005 format(' nlat: ',I10)
C      close(95)
C      close(96)
C      RETURN
C      END

