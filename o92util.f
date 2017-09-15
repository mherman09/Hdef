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
      INTEGER FMAX
      PARAMETER (FMAX=150000)

C---- Input files
      CHARACTER*40 ffmf      ! faults, static out finite fault format
      CHARACTER*40 fspf      ! faults, SRCMOD finite fault format
      CHARACTER*40 fltf      ! faults, lon-lat-dep-str-dip-rak-slip-wid-len format
      CHARACTER*40 magf      ! faults, lon-lat-dep-str-dip-rak-mag format
      CHARACTER*40 haff      ! half-space parameters
      CHARACTER*40 staf      ! station locations
      CHARACTER*40 trgf      ! target/receiver fault geometry
      CHARACTER*40 volf      ! volume sources: planar or isotropic

C---- Output files
      CHARACTER*40 dspf      ! displacement
      CHARACTER*40 stnf      ! strain tensor
      CHARACTER*40 stsf      ! stress tensor
      CHARACTER*40 norf      ! normal stress
      CHARACTER*40 shrf      ! shear stress
      CHARACTER*40 coulf     ! Coulomb stress
      CHARACTER*40 gmtf      ! GMT psxy -SJ faults
      CHARACTER*40 emprel    ! magnitude-fault area empirical relation

C---- Operation switches
      INTEGER flttyp         ! source type flag: point=0, finite=1
      INTEGER xy             ! coordinate flag: geographic=0, Cartesian=1
      INTEGER prog           ! progress indicator flag
      INTEGER long           ! higher precision output flag
      INTEGER verbos         ! verbose mode

C---- Automatic grid variables
      INTEGER auto           ! automatic grid flag: hor=1; dip-par=2, str-par=3
      REAL*8 incr            ! grid increment
      REAL*8 adist           ! grid position

C---- Input variables
      INTEGER nflt           ! number of input faults
      REAL*8 evlo(FMAX)      ! x-coordinate
      REAL*8 evla(FMAX)      ! y-coordinate
      REAL*8 evdp(FMAX)      ! z-coordinate (positive down)
      REAL*8 str(FMAX)       ! strike
      REAL*8 dip(FMAX)       ! dip
      REAL*8 rak(FMAX)       ! rake
      REAL*8 dx(FMAX)        ! along-strike length
      REAL*8 dy(FMAX)        ! along-dip width
      REAL*8 slip(FMAX)      ! slip

      REAL*8 sta(3)          ! station coordinates

      REAL*8 hylo,hyla
      REAL*8 vp,vs,dens
      INTEGER nsta,ntrg,typ(FMAX)

      REAL*8 sthr            ! threshold of slip for FFM

C---- Input file markers
      INTEGER kffm           ! 0:no 1:yes
      INTEGER kfsp           ! 0:no 1:yes
      INTEGER kflt           ! 0:no 1:yes
      INTEGER kmag           ! 0:no 1:yes
      INTEGER ksta           ! 0:no 1:yes
      INTEGER khaf           ! 0:no 1:yes
      INTEGER ktrg           ! 0:no 1:yes 2:cmdln
      INTEGER kvol           ! 0:no 1:yes

C---- Output file markers
      INTEGER kdsp           !
      INTEGER kstn           !
      INTEGER ksts           !
      INTEGER knor           !
      INTEGER kshr           !
      INTEGER kcou           !
      INTEGER kgmt           !
      COMMON /ICHECK/ kffm,kfsp,kflt,kmag,ksta,khaf,ktrg,kvol
      COMMON /OCHECK/ kdsp,kstn,ksts,knor,kshr,kcou,kgmt

C---- Program timing variables
      REAL*8 time1,time2

C----
C Parse command line
C----
      call gcmdln(ffmf,fspf,fltf,magf,haff,staf,trgf,dspf,stnf,stsf,
     1            norf,shrf,coulf,gmtf,flttyp,auto,incr,adist,xy,prog,
     2            long,emprel,volf,sthr,verbos)
      if (verbos.ge.1) then
          write(0,*) '********** RUNNING IN VERBOSE MODE **********'
          write(0,*)
          write(0,*) 'READ THE FOLLOWING FROM THE COMMAND LINE:'
          write(0,*) 'FFMF:   ',trim(ffmf)
          write(0,*) 'FSPF:   ',trim(fspf)
          write(0,*) 'FLTF:   ',trim(fltf)
          write(0,*) 'MAGF:   ',trim(magf)
          write(0,*) 'HAFF:   ',trim(haff)
          write(0,*) 'STAF:   ',trim(staf)
          write(0,*) 'TRGF:   ',trim(trgf)
          write(0,*) 'DSPF:   ',trim(dspf)
          write(0,*) 'STNF:   ',trim(stnf)
          write(0,*) 'STSF:   ',trim(stsf)
          write(0,*) 'NORF:   ',trim(norf)
          write(0,*) 'SHRF:   ',trim(shrf)
          write(0,*) 'COULF:  ',trim(coulf)
          write(0,*) 'GMTF:   ',trim(gmtf)
          write(0,*) 'FLTTYP: ',flttyp
          write(0,*) 'AUTO:   ',auto
          write(0,*) 'INCR:   ',incr
          write(0,*) 'ADIST:  ',adist
          write(0,*) 'XY:     ',xy
          write(0,*) 'PROG:   ',prog
          write(0,*) 'LONG:   ',long
          write(0,*) 'EMPREL: ',emprel
          write(0,*) 'VOLF:   ',volf
          write(0,*) 'STHR:   ',sthr
          write(0,*) 'VERBOS: ',verbos
          write(0,*)
          call cpu_time(time1)
      endif

C----
C Check for input fault source files, and whether output is specified
C----
      if (verbos.ge.1) then
          write(0,*) 'CHECKING SOURCE FAULT INPUTS'
      endif
      call chksrc(ffmf,fspf,fltf,magf,volf)
      if (verbos.ge.1) then
          write(0,*) 'SOURCE FAULTS DEFINED, FILE(S) FOUND'
          write(0,*)
      endif
      if (verbos.ge.1) then
          write(0,*) 'CHECKING DEFORMATION OUTPUTS'
      endif
      call chkout()
      if (verbos.ge.1) then
          write(0,*) 'DEFORMATION OUTPUT DEFINED'
          write(0,*)
      endif
      if (kvol.eq.1) goto 101

C----
C Read input fault files
C----
      if (verbos.ge.1) then
          write(0,*) 'READING SOURCE FAULT INPUTS'
      endif
      call readsrc(ffmf,fspf,fltf,magf,nflt,evlo,evla,evdp,str,dip,rak,
     1             dx,dy,slip,hylo,hyla,emprel,haff,sthr,verbos)
      if (verbos.ge.1) then
          write(0,*) 'SOURCE FAULT READING COMPLETE:',nflt,
     1                                         'SUBFAULTS'
          write(0,*)
      endif
      
C----
C Check for auto flag to define stations
C----
       if (auto.ne.0) then
           call doauto(staf,nflt,evlo,evla,evdp,str,dip,rak,dx,dy,slip,
     1                 auto,adist,incr,xy)
       endif

C----
C Check for input receiver, half-space, and target parameter files
C If resolving stresses, verify NTRG=1 or NTRG=NSTA
C----
      if (verbos.ge.1) then
          write(0,*) 'CHECKING STATION, HALFSPACE, AND TARGET INPUTS'
      endif
      call chkhafstatrg(haff,staf,trgf,nsta,ntrg)
      if (verbos.ge.1) then
          write(0,*) 'OTHER INPUTS DEFINED'
          write(0,*)
      endif
      call readhaff(haff,vp,vs,dens)

C----
C Read in station coordinates and see if they make sense
C----
      call chksta(sta,staf,nsta,xy)

C----
C Compute displacements, strains, and stresses
C----
      if (verbos.ge.1) then
          write(0,*) 'CALCULATING DEFORMATION'
      endif
      if (verbos.ge.2) then
          if (kdsp.ge.1) write(0,*) '    CALCULATING DISPLACEMENTS'
          if (kstn.ge.1) write(0,*) '    CALCULATING STRAINS'
          if (ksts.ge.1) write(0,*) '    CALCULATING STRESSES'
          if (knor.ge.1) write(0,*) '    CALCULATING NORMAL STRESSES'
          if (kshr.ge.1) write(0,*) '    CALCULATING SHEAR STRESSES'
          if (kcou.ge.1) write(0,*) '    CALCULATING COULOMB STRESSES'
      endif
      call calcdefm(staf,trgf,nsta,ntrg,nflt,evlo,evla,evdp,str,dip,rak,
     1              dx,dy,slip,flttyp,vp,vs,dens,dspf,stnf,stsf,norf,
     2              shrf,coulf,xy,prog,long)
      if (verbos.ge.1) then
          write(0,*) 'FINISHED WITH DEFORMATION CALCULATION'
          write(0,*)
      endif

C----
C Compute deformation for volume sources
C----
  101 if (kvol.eq.1) then
          call readvol(volf,nflt,evlo,evla,evdp,str,dip,dx,dy,slip,typ)
          call chkhafstatrg(haff,staf,trgf,nsta,ntrg)
          call readhaff(haff,vp,vs,dens)
          call voldefm(staf,trgf,nsta,ntrg,nflt,evlo,evla,evdp,str,dip,
     1                 dx,dy,slip,typ,vp,vs,dens,dspf,stnf,stsf,
     2                 norf,shrf,coulf,xy,prog,long)
      endif
 
C----
C Write fault input to GMT-compatible file
C----
      if (kgmt.eq.1) then
          call writegmt(gmtf,nflt,evlo,evla,evdp,str,dip,rak,dx,dy,
     1                  slip)
      endif

      if (verbos.ge.1) then
          call cpu_time(time2)
          write(0,9999) time2-time1,nsta*nflt
      endif
 9999 format('O92UTIL TOOK ',F10.3,' SECONDS FOR ',I12,
     1       ' FAULT-STATION PAIRS')
      END

C======================================================================C

      SUBROUTINE chksrc(ffmf,fspf,fltf,magf,volf)
C----
C Check (a) that fault source files were defined in gcmdln and (b) that
C these fault files exist.
C----
      IMPLICIT NONE
      CHARACTER*40 ffmf,fspf,fltf,magf,volf
      LOGICAL ex
      INTEGER kffm,kfsp,kflt,kmag,ksta,khaf,ktrg,kvol
      COMMON /ICHECK/ kffm,kfsp,kflt,kmag,ksta,khaf,ktrg,kvol

C Verify input fault files were defined on command line
      if ((kffm.ne.0.or.kflt.ne.0.or.kmag.ne.0.or.kfsp.ne.0)
     1                                              .and.kvol.eq.1) then
          write(*,*) '!! Error: cannot combine fault and volume sources'
          call usage('!! Run -ffm, -fsp, -flt, or -mag separately '//
     1                                                      'from -vol')
      endif
      if (kffm.eq.0.and.kflt.eq.0.and.kmag.eq.0.and.kvol.eq.0
     1                                              .and.kfsp.eq.0) then
          write(*,*) '!! Error: No source file specified'
          call usage('!! Use -ffm FFMFILE, -fsp FSPFILE, -flt '//
     1                         'FLTFILE, -mag MAGFILE, or -vol VOLFILE')
      endif

C Look for defined files in working directory
      if (kffm.eq.1) then
          inquire(file=ffmf,EXIST=ex)
          if (.not.ex) kffm = -1
      endif
      if (kfsp.eq.1) then
          inquire(file=fspf,EXIST=ex)
          if (.not.ex) kfsp = -1
      endif
      if (kflt.eq.1) then
          inquire(file=fltf,EXIST=ex)
          if (.not.ex) kflt = -1
      endif
      if (kmag.eq.1) then
          inquire(file=magf,EXIST=ex)
          if (.not.ex) kmag = -1
      endif
      if (kvol.eq.1) then
          inquire(file=volf,EXIST=ex)
          if (.not.ex) kflt = -1
      endif

C Quit to usage if none of input is found. Warn (but do not quit) if
C at least one input is found but another is not.
      if (kffm.ne.1.and.kfsp.ne.1.and.kflt.ne.1.and.kmag.ne.1
     1                                              .and.kvol.ne.1) then
          write(*,*) '!! Error: No fault files found'
          if (kffm.eq.-1) write(*,*) '!! Looked for FFMFILE: '//ffmf
          if (kfsp.eq.-1) write(*,*) '!! Looked for FSPFILE: '//fspf
          if (kflt.eq.-1) write(*,*) '!! Looked for FLTFILE: '//fltf
          if (kmag.eq.-1) write(*,*) '!! Looked for MAGFILE: '//magf
          if (kvol.eq.-1) write(*,*) '!! Looked for VOLFILE: '//volf
          call usage('!! Check -ffm, -fsp, -flt, -mag, and -vol '//
     1                                                      'arguments')
      else
          if (kffm.eq.-1) write(*,*) '!! Warning: No FFMFILE: '//ffmf
          if (kfsp.eq.-1) write(*,*) '!! Warning: No FSPFILE: '//fspf
          if (kflt.eq.-1) write(*,*) '!! Warning: No FLTFILE: '//fltf
          if (kmag.eq.-1) write(*,*) '!! Warning: No MAGFILE: '//magf
          if (kvol.eq.-1) write(*,*) '!! Warning: No VOLFILE: '//volf
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
      if (kdsp.le.0.and.kstn.le.0) then
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
      INTEGER kffm,kfsp,kflt,kmag,ksta,khaf,ktrg,kvol
      INTEGER kdsp,kstn,ksts,knor,kshr,kcou,kgmt
      COMMON /ICHECK/ kffm,kfsp,kflt,kmag,ksta,khaf,ktrg,kvol
      COMMON /OCHECK/ kdsp,kstn,ksts,knor,kshr,kcou,kgmt
C Look for station file, count number of receivers
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
C Look for half-space file, warn (not quit) if not found
      if (khaf.eq.0) then
          continue
      elseif (khaf.eq.1) then
          inquire(file=haff,EXIST=ex)
          if (.not.ex) then
              write(*,*) '!! Warning: No half-space file found'
              write(*,*) '!! Looked for HAFFILE: '//haff
              write(*,*) '!! Using default half-space parameters'
              khaf = 0
          endif
      endif
C Check for target parameter file if resolving stress on faults
      ntrg = 0
      if (knor.eq.1.or.kshr.ge.1.or.kcou.eq.1) then
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
C Verify NTRG=1 or NTRG=NSTA
          if (ntrg.ne.nsta.and.ntrg.ne.1) then
              call usage('!! Error: NTRG must be equal to 1 or NSTA')
          endif
      endif
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE linecount(ifile,nline)
C----
C Count number of lines in file
C----
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

      SUBROUTINE chksta(sta,staf,nsta,xy)
      IMPLICIT none
      REAL*8 sta(3)
      CHARACTER*40 staf 
      INTEGER i,nsta,chk,xy
      chk = 0
      open(unit=18,file=staf,status='old')
      do 181 i = 1,nsta
          read(18,*) sta(1),sta(2),sta(3)
          if ((dabs(sta(1)).gt.360.or.dabs(sta(2)).gt.90).and.
     1                                                xy.ne.1) then
              if (chk.eq.0) then
                  write(0,*)'!! Warning: longitude or latitude ',
     1                      'is very large'
                  write(0,*)'!! Did you mean to use -xy option?'
                  chk = 1
              endif
          endif
  181 continue
      close(18)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE readsrc(ffmf,fspf,fltf,magf,nflt,evlo,evla,evdp,str,
     1                   dip,rak,dx,dy,slip,hylo,hyla,emprel,haff,sthr,
     2                   verbos)
C----
C Read shear dislocation parameters from fault source files. Parameter
C FMAX defines maximum number of shear dislocations that can be stored.
C----
      IMPLICIT NONE
      CHARACTER*40 ffmf,fspf,fltf,magf,haff,emprel
      REAL*8 smax
      INTEGER i,j,nflt,FMAX
      PARAMETER (FMAX=150000)
      REAL*8 evlo(FMAX),evla(FMAX),evdp(FMAX),str(FMAX),dip(FMAX),
     1       rak(FMAX),dx(FMAX),dy(FMAX),slip(FMAX),hylo,hyla
      INTEGER verbos
      INTEGER kffm,kfsp,kflt,kmag,ksta,khaf,ktrg,kvol
      INTEGER kdsp,kstn,ksts,knor,kshr,kcou,kgmt
      COMMON /ICHECK/ kffm,kfsp,kflt,kmag,ksta,khaf,ktrg,kvol
      COMMON /OCHECK/ kdsp,kstn,ksts,knor,kshr,kcou,kgmt
      REAL*8 sthr
      nflt = 0
      if (kffm.eq.1) then
          if (verbos.ge.2) then
              write(0,*) '    READING STATIC FFM FILE'
          endif
          call readffm(ffmf,evlo,evla,evdp,str,dip,rak,dx,dy,slip,
     1                 hylo,hyla,nflt)
      endif
      if (kfsp.eq.1) then
          if (verbos.ge.2) then
              write(0,*) '    READING SRCMOD FSP FILE'
          endif
          call readfsp(fspf,evlo,evla,evdp,str,dip,rak,dx,dy,slip,
     1                 hylo,hyla,nflt)
      endif
      if (kflt.eq.1) then
          if (verbos.ge.2) then
              write(0,*) '    READING INPUT FAULTS (SLIP WID LEN)'
          endif
          call readflt(fltf,evlo,evla,evdp,str,dip,rak,dx,dy,slip,
     1                 nflt)
      endif
      if (kmag.eq.1) then
          if (verbos.ge.2) then
              write(0,*) '    READING INPUT FAULTS (MAG)'
          endif
          call readmag(magf,evlo,evla,evdp,str,dip,rak,dx,dy,slip,
     1                 nflt,emprel,haff)
      endif
      if (nflt.lt.1) then
          write(*,*) '!! Error: no input faults read'
          call usage('!! Check input files for correct format')
      endif
 
C Remove all fault sources with slip less than sthr
      smax = 0.0d0
      if (sthr.lt.-1.0d-6) then
          do 300 i = 1,nflt
              if (slip(i).gt.smax) smax = slip(i)
  300     continue
          sthr = dabs(sthr)*smax
      endif

      j = 0
      do 301 i = 1,nflt
          if (slip(i).ge.sthr) then
              j = j + 1
              evlo(j) = evlo(i)
              evla(j) = evla(i)
              evdp(j) = evdp(i)
              str(j) = str(i)
              dip(j) = dip(i)
              rak(j) = rak(i)
              dx(j) = dx(i)
              dy(j) = dy(i)
              slip(j) = slip(i)
          endif
  301 continue
      nflt = j

      if (verbos.ge.3) then
          write(0,*) '    INPUT FAULTS:'
          write(0,*) '                 EVLO'//
     1               '          EVLA'//
     2               '          EVDP'//
     3               '       STR'//
     4               '    DIP'//
     5               '    RAK'//
     6               '        SLIP'//
     7               '         WID'//
     8               '         LEN'
          do 302 i = 1,nflt
              write(0,3021) evlo(i),evla(i),evdp(i),str(i),dip(i),
     1                      rak(i),slip(i),dy(i),dx(i)
  302     continue
      endif
 3021 format(8X,2(F14.6),1F14.1,1F10.1,2F7.1,1F12.3,2F12.1)

C     If FFM not used, use unweighted mean coordinates as hylo and hyla
      if (kffm.ne.1.and.kfsp.ne.1) then
          call calcmean(hylo,evlo,nflt)
          call calcmean(hyla,evla,nflt)
      endif
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE readffm(ffmf,evlo,evla,evdp,str,dip,rak,dx,dy,slip,
     1                   hylo,hyla,nflt)
C----
C Read shear dislocations from FFM in static out subfault format.
C All units are converted to SI, angles are in degrees.
C----
      IMPLICIT NONE
      CHARACTER*40 ffmf,du,dxc,dyc
      INTEGER nflt,FMAX
      PARAMETER (FMAX=150000)
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
C Check that faults do not overflow arrays before reading
          if (ct+nx*ny.gt.FMAX) then
              write(*,*) '!! Error too many input faults (max:',FMAX,')'
              call usage('!! Reduce number of input faults or change '//
     1                             'FMAX')
          endif
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

      SUBROUTINE readfsp(fspf,evlo,evla,evdp,str,dip,rak,dx,dy,slip,
     1                   hylo,hyla,nflt)
      IMPLICIT none
      CHARACTER*40 fspf
      CHARACTER*200 line
      INTEGER i,nflt,FMAX
      PARAMETER (FMAX=150000)
      REAL*8 evlo(FMAX),evla(FMAX),evdp(FMAX),str(FMAX),dip(FMAX),
     1       rak(FMAX),dx(FMAX),dy(FMAX),slip(FMAX),hylo,hyla
      CHARACTER*1 dm
      REAL*8 dxr,dyr,strr,dipr
      i = nflt
      open(unit=35,file=fspf,status='old')
  351 read(35,'(A)',end=352) line
          if (index(line,'%').eq.1) then
              if (index(line,'Loc').gt.0) then
                  read(line,*) dm,dm,dm,dm,dm,hyla,dm,dm,hylo
              elseif (index(line,'Mech').gt.0) then
                  read(line,*) dm,dm,dm,dm,dm,strr,dm,dm,dipr
              elseif (index(line,'Invs').gt.0
     1                                  .and.index(line,'Dx').gt.0) then
                  read(line,*) dm,dm,dm,dm,dm,dxr,dm,dm,dm,dyr
              elseif (index(line,'SEGMENT').gt.0
     1                        .and.index(line,'MULTISEGMENT').eq.0) then
                  read(line,*) dm,dm,dm,dm,dm,dm,strr,dm,dm,dm,dipr
              endif
          elseif (index(line,'%').eq.0) then
              i = i + 1
              read(line,*) evla(i),evlo(i),dm,dm,evdp(i),slip(i),
     1                     rak(i)
              evdp(i) = evdp(i)*1.0d3
              dx(i) = dxr*1.0d3
              dy(i) = dyr*1.0d3
              str(i) = strr
              dip(i) = dipr
          endif
          goto 351
  352 continue
      close(35)
      nflt = nflt + i
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE readflt(fltf,evlo,evla,evdp,str,dip,rak,dx,dy,slip,
     1                   nflt)
C----
C Read shear dislocations from fault file in format:
C   evlo evla evdp(km) str dip rak slip(m) dip_wid(km) str_len(km)
C All units are converted to SI, angles are in degrees.
C----
      IMPLICIT NONE
      CHARACTER*40 fltf
      INTEGER i,f,nflt,FMAX,ln
      PARAMETER (FMAX=150000)
      REAL*8 evlo(FMAX),evla(FMAX),evdp(FMAX),str(FMAX),dip(FMAX),
     1       rak(FMAX),dx(FMAX),dy(FMAX),slip(FMAX)
      call linecount(fltf,ln)
      open (unit=32,file=fltf,status='old')
      f = 0
  321 f = f + 1
      i = nflt + f
C Check that faults do not overflow arrays before reading
      if (i.gt.FMAX) then
          write(*,*) '!! Error too many input faults (max:',FMAX,')'
          call usage('!! Reduce number of input faults or change FMAX')
      endif
      read (32,*,END=322) evlo(i),evla(i),evdp(i),str(i),dip(i),rak(i),
     1                    slip(i),dy(i),dx(i)
          evdp(i) = 1.0d3*evdp(i)
          dx(i) = 1.0d3*dx(i)
          dy(i) = 1.0d3*dy(i)
          goto 321
  322 continue
      if (f-1.ne.ln) then
          call usage('!! Error: NFLT in -flt file not equal to NLINE')
      endif
      nflt = i - 1
      close(32)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE readmag(magf,evlo,evla,evdp,str,dip,rak,dx,dy,
     1                   slip,nflt,emprel,haff)
C----
C Read shear dislocations from fault file in format:
C   evlo evla evdp(km) str dip rak mag
C All units are converted to SI, angles are in degrees.
C----
      IMPLICIT NONE
      CHARACTER*40 magf,emprel
      INTEGER i,f,nflt,FMAX,ln,p
      PARAMETER (FMAX=150000)
      REAL*8 evlo(FMAX),evla(FMAX),evdp(FMAX),str(FMAX),dip(FMAX),
     1       rak(FMAX),dx(FMAX),dy(FMAX),slip(FMAX),mag
      REAL*8 wid,len,mu
      INTEGER typ
      LOGICAL ex
      CHARACTER*40 haff
      REAL*8 vp,vs,dens
      INTEGER kffm,kfsp,kflt,kmag,ksta,khaf,ktrg,kvol
      COMMON /ICHECK/ kffm,kfsp,kflt,kmag,ksta,khaf,ktrg,kvol

      ! Determine empirical relation used
      p = index(emprel,'p')
      if (p.ne.0) then
          emprel(p:p) = ''
      endif
      if (emprel.ne.'WC'.and.emprel.ne.'MB'.and.emprel.ne.'B'
     1                                         .and.emprel.ne.'YM') then
              write(0,*) '!! No empirical relation named "',
     1                                     trim(emprel),'"'
              write(0,*) '   Using Wells and Coppersmith (1994)'
      endif

      ! Shear modulus for converting magnitude-moment to slip
      inquire(file='mu.dat',EXIST=ex) ! Not a cmdln option, but if file exists, read mu
      if (ex) then
          open(unit=61,file='mu.dat',status='old')
          read(61,*) mu
          close(61)
      else
          call readhaff(haff,vp,vs,dens)
          mu = vs*vs*dens
      endif

      call linecount(magf,ln)
      open (unit=33,file=magf,status='old')
      f = 0
  331 f = f + 1
      i = nflt + f
C Check that faults do not overflow arrays
      if (i.gt.FMAX) then
          write(*,*) '!! Error too many input faults (max:',FMAX,')'
          call usage('!! Reduce number of input faults or change FMAX')
      endif
      read (33,*,END=332) evlo(i),evla(i),evdp(i),str(i),dip(i),rak(i),
     1                    mag
          if (45.0d0.le.rak(i).and.rak(i).le.135.0d0) then
              typ = 2 ! Reverse fault (45<rake<135)
          elseif (-135.0d0.lt.rak(i).and.rak(i).lt.-45.0d0) then
              typ = 3 ! Normal fault (-135<rake<-45)
          else
              typ = 1 ! Strike-slip fault (otherwise)
          endif
          if (emprel.eq.'WC') then
              call wellscoppersmith(wid,len,mag,typ)
          elseif (emprel.eq.'MB') then
              call maiberoza(wid,len,mag,typ)
          elseif (emprel.eq.'B') then
              call blaseretal(wid,len,mag,typ)
          elseif (emprel.eq.'YM') then
              call yenma(wid,len,mag,typ)
          else
              call wellscoppersmith(wid,len,mag,typ)
          endif
          !wid = wid*fltfac
          !len = len*fltfac
          slip(i) = ((1.0d1**(1.5d0*(mag+10.7d0)))*1.0d-7)/
     1                                            (mu*wid*len*1.0d6)
          dx(i) = len*1.0d3
          dy(i) = wid*1.0d3
          evdp(i) = 1.0d3*evdp(i)
          goto 331
  332 continue
      if (f-1.ne.ln) then
          call usage('!! Error: NFLT in -mag file not equal to NLINE')
      endif
      nflt = i - 1
      close(33)
      if (p.ne.0) then
          write(0,*) '     EVLO      EVLA  EVDP(KM)    STR    DIP',
     2               '    RAK     SLIP(M)     WID(KM)',
     1               '     LEN(KM)'
          do 333 i = nflt-f+2,nflt
              write(0,9999) evlo(i),evla(i),evdp(i)*1d-3,str(i),dip(i),
     1                      rak(i),slip(i),dy(i)*1d-3,dx(i)*1d-3
  333     continue
      endif
 9999 format(2(F10.3),F10.1,3(F7.1),3(F12.3))
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

      SUBROUTINE maiberoza(wid,len,mag,typ)
C----
C Mai, P.M., Beroza, G.C. (2000). Source scaling properties from finite-
C   fault-rupture models. Bulletin of the Seismological Society of
C   America, vol. 90, no. 3, pp. 604-615.
C Note: using effective width and length scaling.
C----
      IMPLICIT none
      INTEGER typ
      REAL*8 wid,len,mag,mom,logmom
      mom = 10.0d0**(1.5d0*(mag+10.7d0))
      mom = mom*1.0d-7
      logmom = dlog10(mom)
C Strike-slip
      if (typ.eq.1) then
          len = 10.0d0**(-6.31d0+0.40d0*logmom)
          wid = 10.0d0**(-2.18d0+0.17d0*logmom)
C Dip-slip
      elseif (typ.eq.2.or.typ.eq.3) then
          len = 10.0d0**(-6.39d0+0.40d0*logmom)
          wid = 10.0d0**(-5.51d0+0.35d0*logmom)
C All
      else
          len = 10.0d0**(-6.13d0+0.39d0*logmom)
          wid = 10.0d0**(-5.05d0+0.32d0*logmom)
      endif
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE blaseretal(wid,len,mag,typ)
C----
C Blaser, L., Kruger, F., Ohrnberger, M., Scherbaum, F. (2010). Scaling
C   relations of earthquake source parameter estimates with special
C   focus on subduction environment. Bulletin of the Seismological
C   Society of America, vol. 100, no. 6, pp. 2914-2926.
C Note: using orthogonal scaling
C----
      IMPLICIT none
      INTEGER typ
      REAL*8 wid,len,mag
C Strike-slip
      if (typ.eq.1) then
          len = 10.0d0**(-2.69d0+0.64d0*mag)
          wid = 10.0d0**(-1.12d0+0.33d0*mag)
C Reverse
      elseif (typ.eq.2) then
          len = 10.0d0**(-2.37d0+0.57d0*mag)
          wid = 10.0d0**(-1.86d0+0.46d0*mag)
C Normal
      elseif (typ.eq.3) then
          len = 10.0d0**(-1.91d0+0.52d0*mag)
          wid = 10.0d0**(-1.20d0+0.36d0*mag)
C All
      else
          len = 10.0d0**(-2.31d0+0.57d0*mag)
          wid = 10.0d0**(-1.56d0+0.41d0*mag)
      endif
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE yenma(wid,len,mag,typ)
C----
C Yeb, Y.-T., Ma, K.-F. (2011). Source-scaling relationship for
C   M 4.6-8.9 earthquakes, specifically for earthquakes in the collision
C   zone of Taiwan. Bulletin of the Seismological Society of America,
C   vol. 101, no. 2, pp. 464-481.
C----
      IMPLICIT none
      INTEGER typ
      REAL*8 wid,len,mag,mom,logmom
      mom = 10.0d0**(1.5d0*(mag+10.7d0))
      mom = mom*1.0d-7
      logmom = dlog10(mom)
C Strike-slip
      if (typ.eq.1) then
          len = 10.0d0**(-8.11d0+0.50d0*logmom)
          wid = 10.0d0**(-6.67d0+0.42d0*logmom)
C Dip-slip
      elseif (typ.eq.2.or.typ.eq.3) then
          len = 10.0d0**(-6.66d0+0.42d0*logmom)
          wid = 10.0d0**(-5.76d0+0.37d0*logmom)
C All
      else
          len = 10.0d0**(-7.46d0+0.47d0*logmom)
          wid = 10.0d0**(-6.30d0+0.40d0*logmom)
      endif
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE calcmean(mean,array,nval)
C----
C Compute the arithmetic mean of array
C----
      IMPLICIT NONE
      INTEGER i,nval,FMAX
      PARAMETER (FMAX=150000)
      REAL*8 mean,array(FMAX)
      mean = 0.0d0
      do 161 i = 1,nval
          mean = mean + array(i)
  161 continue
      mean = mean/dble(nval)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE readvol(volf,nflt,evlo,evla,evdp,str,dip,dx,dy,
     1                   slip,typ)
      IMPLICIT none
      CHARACTER*40 volf,jnk
      CHARACTER*200 inline
      INTEGER nflt,FMAX
      PARAMETER (FMAX=150000)
      INTEGER typ(FMAX),i,f,ln
      REAL*8 evlo(FMAX),evla(FMAX),evdp(FMAX),str(FMAX),dip(FMAX),
     1       dx(FMAX),dy(FMAX),slip(FMAX)
      nflt = 0
      call linecount(volf,ln)
      open (unit=34,file=volf,status='old')
      f = 0
  341 f = f + 1
      i = nflt + f
C Check that faults do not overflow arrays before reading
      if (i.gt.FMAX) then
          write(*,*) '!! Error too many input faults (max:',FMAX,')'
          call usage('!! Reduce number of input faults or change FMAX')
      endif
C Read faults
      read (34,'(A200)',END=342) inline
          if(index(inline,'p').ne.0) then
              typ(i) = 0
              read(inline,*) jnk,evlo(i),evla(i),evdp(i),slip(i)
              str(i) = 0.0d0
              dip(i) = 0.0d0
              dx(i) = 0.0d0
              dy(i) = 0.0d0
          else
              typ(i) = 1
              read(inline,*) evlo(i),evla(i),evdp(i),str(i),dip(i),
     1                       slip(i),dy(i),dx(i)
              dx(i) = 1.0d3*dx(i)
              dy(i) = 1.0d3*dy(i)
          endif
          evdp(i) = 1.0d3*evdp(i)
          goto 341
  342 continue
      if (f-1.ne.ln) then
          call usage('!! Error: NFLT in -vol file not equal to NLINE')
      endif
      nflt = i - 1
!      print *,'nflt',nflt
!      do 101 i=1,nflt
!          print *,evlo(i),evla(i),evdp(i),slip(i)
!  101 continue
      close(34) 
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE readhaff(haff,vp,vs,dens)
C----
C Read vp, vs, dens from half-space file, or use default values.
C SI units.
C----
      IMPLICIT NONE
      CHARACTER*40 haff,ch
      REAL*8 vp,vs,dens,lamda,mu
      INTEGER kffm,kfsp,kflt,kmag,ksta,khaf,ktrg,kvol
      COMMON /ICHECK/ kffm,kfsp,kflt,kmag,ksta,khaf,ktrg,kvol
      if (khaf.eq.0) then
          call autohaf(vp,vs,dens)
      elseif (khaf.eq.1) then
          open(unit=11,file=haff,status='old')
          read(11,*) ch
          rewind 11
          if (ch(1:1).eq.'L'.or.ch(1:1).eq.'l') then
              ! "Lame" lambda(Pa) shr_mod(Pa)
              read(11,*) ch,lamda,mu
              dens = 3.0d3
              vp = dsqrt((lamda+2.0d0*mu)/dens)
              vs = dsqrt(mu/dens)
          elseif (ch(1:1).eq.'P'.or.ch(1:1).eq.'p') then
              ! "Poisson" shr_mod(Pa) poisson
              read(11,*) ch,mu,lamda ! lamda=poisson's ratio here!!!
              dens = 3.0d3
              vs = dsqrt(mu/dens)
              lamda = 2.0d0*mu*lamda/(1.0d0-2.0d0*lamda) ! This is actually lamda
              vp = dsqrt((lamda+2.0d0*mu)/dens)
          elseif (ch(1:1).eq.'Y'.or.ch(1:1).eq.'y') then
              ! "Young" young_mod(Pa) poisson
              read(11,*) ch,mu,lamda
              dens = 3.0d3
              lamda = mu*lamda/((1.0d0+lamda)*(1.0d0-2.0d0*lamda))
              mu = 0.25d0*(mu-3.0d0*lamda+
     1                    dsqrt(mu*mu+9.0d0*lamda*lamda+2.0d0*mu*lamda))
              vp = dsqrt((lamda+2.0d0*mu)/dens)
              vs = dsqrt(mu/dens)
          else
              ! (Default) vp(m/s) vs(m/s) dens(kg/m^3)
              read(11,*) vp,vs,dens
          endif
      endif
      rewind(11)
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
      REAL*8 pi,r2d
      PARAMETER (pi=4.0d0*atan(1.0d0),r2d=1.8d2/pi)
      CHARACTER*40 staf,trgf,dspf,stnf,stsf,norf,shrf,coulf
      INTEGER i,nsta,ntrg,progout
      REAL*8 stlo,stla,stdp,trgstr,trgdip,trgrak,frict,vp,vs,dens
      INTEGER nflt,FMAX,flttyp,xy,prog,long
      PARAMETER (FMAX=150000)
      REAL*8 evlo(FMAX),evla(FMAX),evdp(FMAX),str(FMAX),dip(FMAX),
     1       rak(FMAX),dx(FMAX),dy(FMAX),slip(FMAX)
      REAL*8 uN,uE,uZ,strain(3,3),stress(3,3),norml,shear,coul,estrs
      REAL*8 az,hdsp
      INTEGER kffm,kfsp,kflt,kmag,ksta,khaf,ktrg,kvol
      INTEGER kdsp,kstn,ksts,knor,kshr,kcou,kgmt
      COMMON /ICHECK/ kffm,kfsp,kflt,kmag,ksta,khaf,ktrg,kvol
      COMMON /OCHECK/ kdsp,kstn,ksts,knor,kshr,kcou,kgmt

C     Open files specified for output
      if (kdsp.ne.0) open(unit=21,file=dspf,status='unknown')
      if (kstn.eq.1) open(unit=22,file=stnf,status='unknown')
      if (ksts.eq.1.or.ksts.eq.3)
     1    open(unit=23,file=stsf,status='unknown')
      if (knor.eq.1) open(unit=24,file=norf,status='unknown')
      if (kshr.ge.1) open(unit=25,file=shrf,status='unknown')
      if (kcou.eq.1) open(unit=26,file=coulf,status='unknown')

C     Initialize progress indicator
      if (prog.eq.1) then
          i = 0
          progout = 0
          call progbar(i,nsta,progout)
      endif

C     If NTRG=1, either from file or command line, read target fault parameters here
      trgstr = -12345.0d0
      trgdip = -12345.0d0
      trgrak = -12345.0d0
      frict  = -12345.0d0
      if (ktrg.eq.1) then ! file
          open(unit=13,file=trgf,status='old')
          if (ntrg.eq.1) then
              read(13,*) trgstr,trgdip,trgrak,frict
          endif
      elseif (ktrg.eq.2) then ! command line
          read(trgf,*) trgstr,trgdip,trgrak,frict
      endif
      ! Check that target faults were read in correctly
      if (ntrg.eq.1) then
          if (trgrak.lt.-10000.0d0) then
              call errmes('!! Error: Not enough arguments in target '//
     1                    'fault geometry input')
          elseif (frict.lt.-10000.0d0) then
              write(0,*) 'Coefficient of friction undefined; using 0.5'
              frict = 0.5d0
          endif
      endif

C     Compute deformation for each receiver location
      open(unit=12,file=staf,status='old')
      do 201 i = 1,nsta
          read(12,*) stlo,stla,stdp
          stdp = stdp*1.0d3
C         If NTRG=NSTA, read target fault parameters here
          if (ktrg.eq.1.and.ntrg.gt.1) then ! already checked ntrg=nsta
              trgstr = -12345.0d0
              trgdip = -12345.0d0
              trgrak = -12345.0d0
              frict  = -12345.0d0
              read(13,*) trgstr,trgdip,trgrak,frict
              ! check target faults read in correctly
              if (trgrak.lt.-10000.0d0) then
                  call errmes('!! Error: Not enough arguments in '//
     1                        'target fault geometry input')
              elseif (frict.lt.-10000.0d0) then
                  write(0,*) 'Coefficient of friction not defined; '//
     1                       'using 0.5'
                  frict = 0.5d0
              endif
          endif
C         Compute displacements
          if (kdsp.ne.0) then
              call calcdisp(uN,uE,uZ,stlo,stla,stdp,nflt,evlo,evla,evdp,
     1                      str,dip,rak,dx,dy,slip,vp,vs,dens,flttyp,xy)
              if (long.eq.0.and.kdsp.eq.1) then
                  write(21,9990) stlo,stla,stdp*1d-3,uE,uN,uZ
              elseif (long.eq.1.and.kdsp.eq.1) then
                  write(21,9993) stlo,stla,stdp*1d-3,uE,uN,uZ
              elseif (kdsp.eq.2) then
                  az = atan2(uE,uN)*r2d
                  hdsp = dsqrt(uE*uE+uN*uN)
                  write(21,9990) stlo,stla,stdp*1d-3,az,hdsp,uz
              endif
          endif
C         Compute (3x3) strain matrix
          if (kstn.ge.1) then
              call calcstn(strain,stlo,stla,stdp,nflt,
     1                     evlo,evla,evdp,str,dip,rak,dx,dy,slip,
     2                     vp,vs,dens,flttyp,xy)
              if (kstn.eq.1) then
                  write(22,9991) stlo,stla,stdp*1d-3,
     1                           strain(1,1),strain(2,2),strain(3,3),
     2                           strain(1,2),strain(1,3),strain(2,3)
              endif
          endif
C         Compute (3x3) stress matrix
          if (ksts.ge.1) then
              call stn2sts(stress,strain,vp,vs,dens)
              if (ksts.eq.1) then
                  write(23,9991) stlo,stla,stdp*1d-3,
     1                           stress(1,1),stress(2,2),stress(3,3),
     2                           stress(1,2),stress(1,3),stress(2,3)
              elseif (ksts.eq.3) then
                  write(23,9991) stlo,stla,stdp*1d-3,estrs(stress)
              endif
          endif
C         Resolve stresses onto target fault planes
          if (kshr.ge.1.or.knor.eq.1.or.kcou.eq.1) then
              call coulomb(coul,norml,shear,stress,trgstr,trgdip,trgrak,
     1                     frict)
              if (kshr.ge.1) then
                  write(25,9992) stlo,stla,stdp*1d-3,shear
              endif
              if (knor.eq.1) then
                  write(24,9992) stlo,stla,stdp*1d-3,norml
              endif
              if (kcou.eq.1) then
                  write(26,9992) stlo,stla,stdp*1d-3,coul
              endif
          endif
C         Update progress indicator
          if (prog.eq.1) then
              call progbar(i,nsta,progout)
          endif
  201 continue
      close(12)
 9990 format(3F12.4,X,3F12.4)
 9991 format(3F12.4,X,6E16.6)
 9992 format(3F12.4,X,1E16.6)
 9993 format(3F12.4,X,3F16.8)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE progbar(prog,prog100,progout)
C----
C Display a progress indicator
C----
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
      PARAMETER (FMAX=150000)
      REAL*8 evlo(FMAX),evla(FMAX),evdp(FMAX),str(FMAX),dip(FMAX),
     1       rak(FMAX),slip(FMAX),dx(FMAX),dy(FMAX),area
      INTEGER flttyp,xy
      REAL*8 vp,vs,dens
C Initialize net displacements
      uNnet = 0.0d0
      uEnet = 0.0d0
      uZnet = 0.0d0
C For each fault, compute contribution to net displacement
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
C x points along strike, y points horizontal, up-dip
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
      PARAMETER (FMAX=150000)
      REAL*8 evlo(FMAX),evla(FMAX),evdp(FMAX),str(FMAX),dip(FMAX),
     1       rak(FMAX),slip(FMAX),dx(FMAX),dy(FMAX),area
      INTEGER flttyp,xy
      REAL*8 vp,vs,dens
      INTEGER i,j
C Initialize net strains
      do 222 i = 1,3
          do 221 j = 1,3
              strain(i,j) = 0.0d0
  221     continue
  222 continue
C For each fault, compute contribution to net strain
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
C x points along strike, y points horizontal, up-dip
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

      SUBROUTINE xy2NE(uN,uE,ux,uy,str)
C----
C Rotate x,y (along-strike,hor. up-dip) to E,N
C----
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

C----------------------------------------------------------------------C

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

C----------------------------------------------------------------------C

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

      FUNCTION estrs(sts)
      IMPLICIT none
      REAL*8 estrs,sts(3,3),s12,s13,s23
      s12 = sts(1,1)-sts(2,2)
      s13 = sts(1,1)-sts(3,3)
      s23 = sts(2,2)-sts(3,3)
      estrs = (1.0d0/6.0d0)*(s12*s12+s13*s13+s23*s23)
     1          + sts(1,2)*sts(1,2)+sts(1,3)*sts(1,3)+sts(2,3)*sts(2,3)
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
      INTEGER kdsp,kstn,ksts,knor,kshr,kcou,kgmt
      COMMON /OCHECK/ kdsp,kstn,ksts,knor,kshr,kcou,kgmt
      str = strin*d2r
      dip = dipin*d2r
      rak = rakin*d2r
C Unit normal vector to plane
      n(1) = dsin(dip)*dsin(str+pi/2.0d0)
      n(2) = dsin(dip)*dcos(str+pi/2.0d0)
      n(3) = dcos(dip)
C Traction is stress matrix times normal vector: t = S*n
      do 241 i = 1,3
          trac(i) = stress(i,1)*n(1)+stress(i,2)*n(2)+stress(i,3)*n(3)
  241 continue
C Normal component of traction is parallel to unit normal vector
      norml = 0.0d0
      do 242 i = 1,3
          norml = norml + trac(i)*n(i)
  242 continue
C Shear component of traction is difference between total and normal
      s(1) = trac(1) - norml*n(1)
      s(2) = trac(2) - norml*n(2)
      s(3) = trac(3) - norml*n(3)
C Compute unit slip vector
      r(1) =  dcos(str-pi/2.0d0)*dcos(rak)
     1                          + dsin(str-pi/2.0d0)*dsin(rak)*dcos(dip)
      r(2) = -dsin(str-pi/2.0d0)*dcos(rak)
     1                          + dcos(str-pi/2.0d0)*dsin(rak)*dcos(dip)
      r(3) =                                         dsin(rak)*dsin(dip)
C Project shear component on slip vector
      if (kshr.eq.2) then
          shear = dsqrt(s(1)*s(1) + s(2)*s(2) + s(3)*s(3))
      else
          shear = s(1)*r(1) + s(2)*r(2) + s(3)*r(3)
      endif
C Coulomb stress (recall sign convention: pos = dilation)
      coul = shear + frc*norml
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE doauto(staf,nflt,evlo,evla,evdp,str,dip,rak,dx,dy,slip,
     1                  auto,adist,incr,xy)
C----
C Determine the location and extent of a horizontal or vertical slice
C for computing the deformation.
C----
      IMPLICIT NONE
      REAL*8 pi,d2k
      PARAMETER (pi=4.0d0*datan(1.0d0),d2k=6.371d3*2.0d0*pi/3.6d2)
      CHARACTER*40 staf
      INTEGER nflt,FMAX
      PARAMETER (FMAX=150000)
      REAL*8 evlo(FMAX),evla(FMAX),evdp(FMAX),str(FMAX),dip(FMAX),
     1       rak(FMAX),dx(FMAX),dy(FMAX),slip(FMAX)
      INTEGER auto,xy
      REAL*8 incr,adist,mom,str0,W,E,S,N
      REAL*8 xmin,xmax,ymin,ymax,x0,y0,z0,lmin,lmax,zmin,zmax
C Get moment, moment-weighted centroid/strike of input faults
      call centroid(x0,y0,z0,mom,nflt,evlo,evla,evdp,dx,dy,slip)
      call strike(str0,nflt,str,dx,dy,slip)
C Make sure grid increments are in correct units
      if (xy.eq.0.and.incr.gt.0) then
          if (auto.eq.2.or.auto.eq.3) then
              incr = incr*d2k
          endif
      elseif (xy.eq.1.and.incr.gt.0) then
          incr = incr*d2k
      elseif (xy.eq.0.and.incr.lt.0) then
          if (auto.eq.1) then
              incr = -incr/d2k
          elseif (auto.eq.2.or.auto.eq.3) then
              incr = -incr
          endif
      elseif (xy.eq.1.and.incr.lt.0) then
          incr = -incr
      endif
C Create horizontal or vertical grid
      if (auto.eq.1) then
          call maxhlims(xmin,xmax,ymin,ymax,x0,y0,xy)
          call autohlims(W,E,S,N,xmin,xmax,ymin,ymax,x0,y0,
     2                   evlo,evla,evdp,str,dip,rak,dx,dy,slip,nflt,
     3                   xy,adist)
C Round to nearest hundredth
          W = dnint(W*1.0d2)*1.0d-2
          E = dnint(E*1.0d2)*1.0d-2
          S = dnint(S*1.0d2)*1.0d-2
          N = dnint(N*1.0d2)*1.0d-2
          write(*,2000)
          write(*,2001) x0,y0,z0*1d-3
          write(*,2002) str0
          write(*,2003)
          write(*,2004) 'HORIZONTAL'
          write(*,2005) adist
          write(*,2008) 'W: ',W
          write(*,2008) 'E: ',E
          write(*,2008) 'S: ',S
          write(*,2008) 'N: ',N
          write(*,2009) incr
          write(*,2010) 'NX: ',int((E-W)/incr)+1
          write(*,2010) 'NY: ',int((N-S)/incr)+1
          call hgrid(staf,W,E,S,N,incr,adist)
      elseif (auto.eq.2.or.auto.eq.3) then
          call autovlims(lmin,lmax,zmin,zmax,x0,y0,z0,str0,
     1                   evlo,evla,evdp,str,dip,rak,dx,dy,slip,nflt,
     2                   xy,auto)
C Round to nearest hundredth
          lmin = dnint(lmin*1.0d2)*1.0d-2
          lmax = dnint(lmax*1.0d2)*1.0d-2
          zmin = dnint(zmin*1.0d2)*1.0d-2
          zmax = dnint(zmax*1.0d2)*1.0d-2
          write(*,2000)
          write(*,2001) x0,y0,z0*1d-3
          write(*,2002) str0
          write(*,2003)
          write(*,2004) 'VERTICAL'
          if (auto.eq.2) then
              write(*,2006) str0-90.0d0
          else
              write(*,2006) str0
          endif
          write(*,2007) adist
          write(*,2008) 'LMIN: ',lmin
          write(*,2008) 'LMAX: ',lmax
          write(*,2008) 'ZMIN: ',zmin
          write(*,2008) 'ZMAX: ',zmax
          write(*,2009) incr
          write(*,2010) 'NL: ',int((lmax-lmin)/incr)+1
          write(*,2010) 'NZ: ',int((zmax-zmin)/incr)+1
          call vgrid(staf,lmin,lmax,zmin,zmax,incr,x0,y0,str0,xy,auto,
     1               adist)
      else
          call usage('!! Error: this option doesnt exist yet...sorry!')
      endif
 2000 format('FAULT PARAMETERS')
 2001 format('  CENTROID: ',3F10.4)
 2002 format('  STRIKE:   ',F10.2)
 2003 format('GRID PARAMETERS')
 2004 format('  ORIENTATION: ',A)
 2005 format('  DEPTH:    ',F10.2)
 2006 format('  GRID_AZ:  ',F10.2)
 2007 format('  HOR_SHFT: ',F10.2)
 2008 format(A5,F10.2)
 2009 format('  GRID_INCR:',F10.2)
 2010 format(A6,I8)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE maxhlims(WMAX,EMAX,SMAX,NMAX,hylo,hyla,xy)
      IMPLICIT NONE
      REAL*8 pi,r2d,d2r
      PARAMETER (pi=4.0d0*atan(1.0d0),r2d=180.0d0/pi,d2r=pi/180.0d0)
      REAL*8 az,dmax,WMAX,EMAX,SMAX,NMAX,hylo,hyla,junk
      INTEGER xy
C----
C Determine maximum geographic limits of receiver grid to W,E,S,N.
C----
      dmax = 750.0d0 ! km
      if (xy.eq.0) then
          WMAX = hylo-dmax*3.6d2/(2.0d0*pi*6.371d3*dcos(hyla*d2r))
          EMAX = hylo+dmax*3.6d2/(2.0d0*pi*6.371d3*dcos(hyla*d2r))
          az = 180.0d0
          call dlola(junk,SMAX,hylo,hyla,dmax,az) ! output in radians
          SMAX = SMAX*r2d
          az = 0.0d0
          call dlola(junk,NMAX,hylo,hyla,dmax,az)
          NMAX = NMAX*r2d
      else
          WMAX = hylo - dmax
          EMAX = hylo + dmax
          SMAX = hyla - dmax
          NMAX = hyla + dmax
      endif
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE autohlims(xmin,xmax,ymin,ymax,xmin0,xmax0,ymin0,ymax0,
     1                     x0,y0,
     2                     evlo,evla,evdp,str,dip,rak,dx,dy,slip,nflt,
     3                     xy,adist)
C----
C Get limits of receiver grid by taking steps along longitude and
C latitudes until displacements smaller than minimum threshold value.
C Variable 'dir' determines direction: 1=N, 2=S, 3=E, 4=W
C----
      IMPLICIT NONE
      INTEGER i,NTHR
      PARAMETER (NTHR=3)
      REAL*8 dl(NTHR),thr(NTHR)
      REAL*8 xmin,xmax,ymin,ymax,xmin0,xmax0,ymin0,ymax0,x,y
      REAL*8 x0,y0
      INTEGER nflt,FMAX
      PARAMETER (FMAX=150000)
      REAL*8 evlo(FMAX),evla(FMAX),evdp(FMAX),str(FMAX),dip(FMAX),
     1       rak(FMAX),dx(FMAX),dy(FMAX),slip(FMAX)
      REAL*8 uE,uN,uZ,stdp,vp,vs,dens,adist
      INTEGER xy,dir,flttyp
      if (evdp(1).gt.1.0d5) then
          xmin = x0 - 3.0d0/dcos(evla(1)*0.01745d0)
          xmax = x0 + 3.0d0/dcos(evla(1)*0.01745d0)
          ymin = y0 - 3.0d0
          ymax = y0 + 3.0d0
          return
      endif
      if (xy.eq.0) then
          dl(1) = 1.0d0 ! Amounts to move test location (deg)
          dl(2) = 0.1d0
          dl(3) = 0.05d0
      else
          dl(1) = 100.0d0 ! Amounts to move test location (km)
          dl(2) = 10.0d0
          dl(3) = 5.0d0
      endif
      thr(1) = 0.5d0 ! Displacement thresholds (m)
      thr(2) = 0.01d0
      thr(3) = 0.005d0
      flttyp = 1
      call autohaf(vp,vs,dens)
      stdp = adist*1d3
      do 903 dir = 1,4
          x = x0
          y = y0
C         Take first step (by increment of dl(1))
          call hstep(x,y,dl(1),dir)
C         Take steps of decreasing size and threshold value
          i = 1
 901      if (i.gt.NTHR) goto 902
C             Check if lon, lat exceed maximum allowed values
              if (dir.eq.1.and.y.gt.ymax0) then
                  y = ymax0
                  goto 902
              elseif (dir.eq.2.and.y.lt.ymin0) then
                  y = ymin0
                  goto 902
              elseif (dir.eq.3.and.x.gt.xmax0) then
                  x = xmax0
                  goto 902
              elseif (dir.eq.4.and.x.lt.xmin0) then
                  x = xmin0
                  goto 902
              endif
C             Calculate disp, compare each component to threshold
c             stdp is input to calcdisp in SI units
              call calcdisp(uN,uE,uZ,x,y,stdp,nflt,evlo,evla,evdp,str,
     1                      dip,rak,dx,dy,slip,vp,vs,dens,flttyp,xy)
              uN = dabs(uN)
              uE = dabs(uE)
              uZ = dabs(uZ)
              if (uN.gt.thr(i).or.uE.gt.thr(i).or.uZ.gt.thr(i)) then
                  call hstep(x,y,dl(i),dir)
              else
                  i = i + 1
              endif
              goto 901
  902     continue
          if (dir.eq.1) then
              ymax = y
          elseif (dir.eq.2) then
              ymin = y
          elseif (dir.eq.3) then
              xmax = x
          elseif (dir.eq.4) then
              xmin = x
          endif
  903 continue
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE hstep(lo,la,siz,dir)
      IMPLICIT none
      INTEGER dir
      REAL*8 lo,la,siz
      if (dir.eq.1) then
          la = la + siz
      elseif (dir.eq.2) then
          la = la - siz
      elseif (dir.eq.3) then
          lo = lo + siz
      elseif (dir.eq.4) then
          lo = lo - siz
      endif
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE autovlims(lmin,lmax,zmin,zmax,
     1                     x0,y0,z0,str0,
     1                     evlo,evla,evdp,str,dip,rak,dx,dy,slip,nflt,
     2                     xy,auto)
C----
C Get limits of receiver grid by taking steps along strike or dip and
C depth until displacements smaller than minimum threshold value.
C Variable 'dir' determines direction: 1=up, 2=down, 3=along-strike/up-dip, 4=opposite of 3
C z0 is input in meters
C----
      IMPLICIT none
      !REAL*8 pi,r2d,d2r
      !PARAMETER (pi=4.0d0*atan(1.0d0),r2d=180.0d0/pi,d2r=pi/180.0d0)
      INTEGER i,NTHR
      PARAMETER (NTHR=3)
      REAL*8 dl(NTHR),thr(NTHR)
      REAL*8 lmin,lmax,zmin,zmax,dmax,x0,y0,z0,str0,x,y,z
      REAL*8 uE,uN,uZ,vp,vs,dens,az,shftaz,l
      INTEGER nflt,FMAX,xy,dir,flttyp,auto
      PARAMETER (FMAX=150000)
      REAL*8 evlo(FMAX),evla(FMAX),evdp(FMAX),str(FMAX),dip(FMAX),
     1       rak(FMAX),dx(FMAX),dy(FMAX),slip(FMAX)
      DMAX = 750.0d0 ! max dist (km)
      dl(1) = 100.0d0 ! Amounts to move test location (km)
      dl(2) = 10.0d0
      dl(3) = 5.0d0
      thr(1) = 0.5d0 ! Displacement thresholds (m)
      thr(2) = 0.01d0
      thr(3) = 0.005d0
      if (auto.eq.2) then
          az = str0 - 90.0d0
      elseif (auto.eq.3) then
          az = str0
      endif
      shftaz = az + 90.0d0
      flttyp = 1
      call autohaf(vp,vs,dens)
      do 911 dir = 1,4
          z = z0
          x = x0
          y = y0
          l = 0.0d0
C Take first step by increment dl(1)
          call vstep(x,y,z,l,dl(1),az,dir,xy)
C Take steps of decreasing size and threshold value
          i = 1
  912     if (i.gt.NTHR) goto 913
C Check if depth, hor. extent exceed max allowed values
              if (dir.eq.1.and.z.lt.0.0d0) then
                  z = 0.0d0
                  goto 913
              elseif (dir.eq.2.and.z.gt.DMAX*1d3+z0) then
                  z = DMAX*1d3+z0
                  goto 913
              elseif (dir.eq.3.and.l.gt.DMAX) then
                  l = DMAX
                  goto 913
              elseif (dir.eq.4.and.l.gt.DMAX) then
                  l = DMAX
                  goto 913
              endif
C             Calculate disp, compare each component to threshold
c             stdp is input to calcdisp in SI units
              call calcdisp(uN,uE,uZ,x,y,z,nflt,evlo,evla,evdp,str,
     1                      dip,rak,dx,dy,slip,vp,vs,dens,flttyp,xy)
              uN = dabs(uN)
              uE = dabs(uE)
              uZ = dabs(uZ)
              !if(dir.lt.3)write(*,'(3E12.4,F12.2)'),uN,uE,uZ,z
              if (uN.gt.thr(i).or.uE.gt.thr(i).or.uZ.gt.thr(i)) then
                  call vstep(x,y,z,l,dl(i),az,dir,xy)
              else
                  i = i + 1
              endif
              goto 912
  913     continue
          if (dir.eq.1) then
              zmin = z
          elseif (dir.eq.2) then
              zmax = z*1d-3
          elseif (dir.eq.3) then
              lmax = l
          elseif (dir.eq.4) then
              lmin = -l
          endif
  911 continue
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE vstep(x,y,z,l,siz,az,dir,xy)
      IMPLICIT none
      REAL*8 pi,r2d,d2r
      PARAMETER (pi=4.0d0*atan(1.0d0),r2d=180.0d0/pi,d2r=pi/180.0d0)
      INTEGER dir,xy
      REAL*8 x,y,z,siz,az,loout,laout,l
      if (dir.eq.1) then
          z = z - siz*1e3
      elseif (dir.eq.2) then
          z = z + siz*1e3
      elseif (dir.eq.3) then
          if (xy.eq.0) then
              call dlola(loout,laout,x,y,siz,az)
              x = loout*r2d
              y = laout*r2d
          elseif (xy.eq.1) then
              x = x+siz*dsin(az*d2r)
              y = y+siz*dcos(az*d2r)
          endif
      elseif (dir.eq.4) then
          az = az + 180.0d0
          if (xy.eq.0) then
              call dlola(loout,laout,x,y,siz,az)
              x = loout*r2d
              y = laout*r2d
          elseif (xy.eq.1) then
              x = x+siz*dsin(az*d2r)
              y = y+siz*dcos(az*d2r)
          endif
      endif
      l = l + siz
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE vgrid(staf,lmin,lmax,zmin,zmax,incr,x0,y0,str0,xy,auto,
     1                 shift)
C----
C Write a vertical grid to file staf.
C----
      IMPLICIT none
      REAL*8 pi,r2d,d2r
      PARAMETER (pi=4.0d0*atan(1.0d0),r2d=180.0d0/pi,d2r=pi/180.0d0)
      CHARACTER*40 staf
      REAL*8 lmin,lmax,zmin,zmax,incr,x0,y0,str0,dx,x1,y1,z,az,shift,
     1       shftaz
      INTEGER xy,auto
      if (auto.eq.2) then
          az = str0 - 90.0d0
      elseif (auto.eq.3) then
          az = str0
      endif
C Shift center of grid by specified amount
      shftaz = az + 90.0d0
      call dlola(x1,y1,x0,y0,shift,shftaz)
      x0 = x1*r2d
      y0 = y1*r2d
      open(unit=96,file=staf,status='unknown')
      dx = lmin
  931 if (dx.le.lmax) then
          if (xy.eq.0) then
              call dlola(x1,y1,x0,y0,dx,az)
              x1 = x1*r2d
              y1 = y1*r2d
          else
              x1 = x0 + dx*dsin(az*d2r)
              y1 = y0 + dx*dcos(az*d2r)
          endif
          z = zmin
  932     if (z.le.zmax) then
              write(96,*) x1,y1,z,dx,sin(az*d2r),cos(az*d2r)
              z = z + incr
              goto 932
          endif
          dx = dx + incr
          goto 931
      endif
      close(96)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE hgrid(staf,xmin,xmax,ymin,ymax,incr,dep)
C----
C Write a horizontal grid to file staf.
C----
      IMPLICIT none
      CHARACTER*40 staf
      REAL*8 xmin,xmax,ymin,ymax,incr,dep,x,y
      open(unit=95,file=staf,status='unknown')
      x = xmin
  951 if (x.le.xmax) then
          y = ymin
  952     if (y.le.ymax) then
              write(95,*) x,y,dep
              y = y + incr
              goto 952
          endif
          x = x + incr
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
      PARAMETER (FMAX=150000)
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
      PARAMETER (FMAX=150000)
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

      SUBROUTINE voldefm(staf,trgf,nsta,ntrg,nflt,evlo,evla,evdp,str,
     1                   dip,dx,dy,slip,typ,vp,vs,dens,dspf,stnf,
     2                   stsf,norf,shrf,coulf,xy,prog,long)
      IMPLICIT none
      REAL*8 pi,r2d
      PARAMETER (pi=4.0d0*atan(1.0d0),r2d=1.8d2/pi)
      CHARACTER*40 staf,trgf,dspf,stnf,stsf,norf,shrf,coulf
      INTEGER i,nsta,ntrg,progout
      REAL*8 stlo,stla,stdp,trgstr,trgdip,trgrak,frict,vp,vs,dens
      INTEGER nflt,FMAX,xy,prog,long
      PARAMETER (FMAX=150000)
      REAL*8 evlo(FMAX),evla(FMAX),evdp(FMAX),str(FMAX),dip(FMAX),
     1       dx(FMAX),dy(FMAX),slip(FMAX)
      REAL*8 uN,uE,uZ,strain(3,3),stress(3,3),norml,shear,coul,estrs
      REAL*8 az,hdsp
      INTEGER typ(FMAX)
      INTEGER kffm,kfsp,kflt,kmag,ksta,khaf,ktrg,kvol
      INTEGER kdsp,kstn,ksts,knor,kshr,kcou,kgmt
      COMMON /ICHECK/ kffm,kfsp,kflt,kmag,ksta,khaf,ktrg,kvol
      COMMON /OCHECK/ kdsp,kstn,ksts,knor,kshr,kcou,kgmt
C     Open files specified for output
      if (kdsp.ne.0) open(unit=21,file=dspf,status='unknown')
      if (kstn.eq.1) open(unit=22,file=stnf,status='unknown')
      if (ksts.eq.1.or.ksts.eq.3)
     1    open(unit=23,file=stsf,status='unknown')
      if (knor.eq.1) open(unit=24,file=norf,status='unknown')
      if (kshr.ge.1) open(unit=25,file=shrf,status='unknown')
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
C         If NTRG=NSTA, read target fault parameters here
          if (ntrg.eq.nsta.and.ktrg.eq.1.and.nsta.ne.1) then
              read(13,*) trgstr,trgdip,trgrak,frict
          endif
C         Compute displacements
          if (kdsp.ne.0) then
              call voldisp(uN,uE,uZ,stlo,stla,stdp,nflt,evlo,evla,evdp,
     1                      str,dip,dx,dy,slip,vp,vs,dens,typ,xy)
              if (long.eq.0.and.kdsp.eq.1) then
                  write(21,9990) stlo,stla,stdp*1d-3,uE,uN,uZ
              elseif (long.eq.1.and.kdsp.eq.1) then
                  write(21,9993) stlo,stla,stdp*1d-3,uE,uN,uZ
              elseif (kdsp.eq.2) then
                  az = atan2(uE,uN)*r2d
                  hdsp = dsqrt(uE*uE+uN*uN)
                  write(21,9990) stlo,stla,stdp*1d-3,az,hdsp,uz
              endif
          endif
C         Compute (3x3) strain matrix
          if (kstn.ge.1) then
              call volstn(strain,stlo,stla,stdp,nflt,
     1                     evlo,evla,evdp,str,dip,dx,dy,slip,
     2                     vp,vs,dens,typ,xy)
              if (kstn.eq.1) then
                  write(22,9991) stlo,stla,stdp*1d-3,
     1                           strain(1,1),strain(2,2),strain(3,3),
     2                           strain(1,2),strain(1,3),strain(2,3)
              endif
          endif
C         Compute (3x3) stress matrix
          if (ksts.ge.1) then
              call stn2sts(stress,strain,vp,vs,dens)
              if (ksts.eq.1) then
                  write(23,9991) stlo,stla,stdp*1d-3,
     1                           stress(1,1),stress(2,2),stress(3,3),
     2                           stress(1,2),stress(1,3),stress(2,3)
              elseif (ksts.eq.3) then
                  write(23,9991) stlo,stla,stdp*1d-3,estrs(stress)
              endif
          endif
C         Resolve stresses onto target fault planes
          if (kshr.ge.1.or.knor.eq.1.or.kcou.eq.1) then
              call coulomb(coul,norml,shear,stress,trgstr,trgdip,trgrak,
     1                     frict)
              if (kshr.ge.1) then
                  write(25,9992) stlo,stla,stdp*1d-3,shear
              endif
              if (knor.eq.1) then
                  write(24,9992) stlo,stla,stdp*1d-3,norml
              endif
              if (kcou.eq.1) then
                  write(26,9992) stlo,stla,stdp*1d-3,coul
              endif
          endif
C         Update progress indicator
          if (prog.eq.1) then
              call progbar(i,nsta,progout)
          endif
  201 continue
      close(12)
 9990 format(3F12.4,X,3F12.4)
 9991 format(3F12.4,X,6E16.6)
 9992 format(3F12.4,X,1E16.6)
 9993 format(3F12.4,X,3F16.8)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE voldisp(uNnet,uEnet,uZnet,stlo,stla,stdp,nflt,evlo,
     1                   evla,evdp,str,dip,dx,dy,slip,vp,vs,dens,typ,xy)
C----
C Compute north, east, and vertical (positive up) displacements at
C (stlo,stla,stdp) from volume sources.
C----
      IMPLICIT NONE
      REAL*8 pi,d2r
      PARAMETER (pi=4.0d0*datan(1.0d0),d2r=pi/1.8d2)
      REAL*8 uNnet,uEnet,uZnet,ux,uy,uN,uE,uz
      REAL*8 stlo,stla,stdp,dist,az,x,y,delx,dely
      INTEGER f,nflt,FMAX
      PARAMETER (FMAX=150000)
      REAL*8 evlo(FMAX),evla(FMAX),evdp(FMAX),str(FMAX),dip(FMAX),
     1       rak(FMAX),slip(FMAX),dx(FMAX),dy(FMAX)
      INTEGER typ(FMAX),xy
      REAL*8 vp,vs,dens
C Initialize net displacements
      uNnet = 0.0d0
      uEnet = 0.0d0
      uZnet = 0.0d0
C For each fault, compute contribution to net displacement
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
C x points along strike, y points horizontal, up-dip
          x = dist*( dcos(az-d2r*str(f)))
          y = dist*(-dsin(az-d2r*str(f)))
          if (typ(f).eq.0) then
              call o92ptvol(ux,uy,uz,x,y,stdp,evdp(f),slip(f),vp,vs,
     1                      dens)
          else
              call o92rectvol(ux,uy,uz,x,y,stdp,evdp(f),dip(f),rak(f),
     1                        dy(f),dx(f),slip(f),vp,vs,dens)
          endif
          call xy2NE(uN,uE,ux,uy,str(f))
          uNnet = uNnet + uN
          uEnet = uEnet + uE
          uZnet = uZnet + uz
  211 continue
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE volstn(strain,stlo,stla,stdp,nflt,
     1                  evlo,evla,evdp,str,dip,dx,dy,slip,
     2                  vp,vs,dens,typ,xy)
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
      PARAMETER (FMAX=150000)
      REAL*8 evlo(FMAX),evla(FMAX),evdp(FMAX),str(FMAX),dip(FMAX),
     1       slip(FMAX),dx(FMAX),dy(FMAX)
      INTEGER typ(FMAX),xy
      REAL*8 vp,vs,dens
      INTEGER i,j
C Initialize net strains
      do 222 i = 1,3
          do 221 j = 1,3
              strain(i,j) = 0.0d0
  221     continue
  222 continue
C For each fault, compute contribution to net strain
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
C x points along strike, y points horizontal, up-dip
          x = dist*( dcos(az-d2r*str(f)))
          y = dist*(-dsin(az-d2r*str(f)))
          if (typ(f).eq.0) then
              call o92ptstnvol(stntmp,x,y,stdp,evdp(f),slip(f),vp,vs,
     1                         dens)
          else
              call o92rectstnvol(stntmp,x,y,stdp,evdp(f),dip(f),dy(f),
     1                           dx(f),slip(f),vp,vs,dens)
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
      PARAMETER (FMAX=150000)
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

      SUBROUTINE gcmdln(ffmf,fspf,fltf,magf,haff,staf,trgf,dspf,stnf,
     1                  stsf,norf,shrf,coulf,gmtf,flttyp,auto,incr,
     2                  adist,xy,prog,long,emprel,volf,sthr,verbos)
      IMPLICIT NONE
      CHARACTER*40 tag
      INTEGER narg,i,ptr,iargc
      CHARACTER*40 ffmf,fspf,fltf,magf,haff,staf,trgf,dspf,stnf,stsf,
     1             norf,shrf,coulf,gmtf,emprel,volf
      INTEGER flttyp,auto,xy,prog,long,verbos
      REAL*8 incr,adist,sthr
      INTEGER kffm,kfsp,kflt,kmag,ksta,khaf,ktrg,kvol
      INTEGER kdsp,kstn,ksts,knor,kshr,kcou,kgmt
      COMMON /ICHECK/ kffm,kfsp,kflt,kmag,ksta,khaf,ktrg,kvol
      COMMON /OCHECK/ kdsp,kstn,ksts,knor,kshr,kcou,kgmt
C Initialize variables
      ffmf = 'none'
      fspf = 'none'
      fltf = 'none'
      magf = 'none'
      staf = 'none'
      haff = 'none'
      trgf = 'none'
      volf = 'none'
      dspf = 'none'
      stnf = 'none'
      stsf = 'none'
      norf = 'none'
      shrf = 'none'
      coulf = 'none'
      gmtf = 'none'
      emprel = 'WC'
      flttyp = 1
      auto = 0
      xy = 0
      prog = 0
      incr = 1.0d99
      adist = 1.0d99
      long = 0
      sthr = 0.0d0
      verbos = 0
C Input file status indicators (ICHECK block)
C 0=unspecified; 1=specified; -1=specified, but file not found
      kffm = 0
      kfsp = 0
      kflt = 0
      kmag = 0
      ksta = 0
      khaf = 0
      ktrg = 0
      kvol = 0
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
      if (narg.eq.0) call usage('!! o92util: tool for computing '//
     1                    'deformation based on Okada (1992)')
      i = 0
  101 i = i + 1
      if (i.gt.narg) goto 102
          call getarg(i,tag)
          if (tag(1:4).eq.'-ffm') then
              kffm = 1
              i = i + 1
              call getarg(i,ffmf)
          elseif (tag(1:4).eq.'-fsp') then
              kfsp = 1
              i = i + 1
              call getarg(i,fspf)
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
              if (kdsp.eq.0) then
                  kdsp = 1
              elseif (kdsp.eq.-1) then
                  kdsp = 2
              endif
              i = i + 1
              call getarg(i,dspf)
          elseif (tag(1:4).eq.'-thr') then
              i = i + 1
              call getarg(i,tag)
              read(tag,'(BN,F10.0)') sthr
          elseif (tag(1:3).eq.'-az') then
              if (kdsp.eq.1) then
                  kdsp = 2
              else
                  kdsp = -1
              endif
          elseif (tag(1:7).eq.'-strain') then
              kstn = 1
              i = i + 1
              call getarg(i,stnf)
          elseif (tag(1:8).eq.'-estress') then
              if (kstn.eq.0) kstn = 2
              ksts = 3
              i = i + 1
              call getarg(i,stsf)
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
          elseif (tag(1:9).eq.'-shearmax') then
              if (kstn.eq.0) kstn = 2
              if (ksts.eq.0) ksts = 2
              kshr = 2
              i = i + 1
              call getarg(i,shrf)
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
              if (ksta.eq.0) staf = 'autosta.dat'
              ksta = 1
              ! Read auto mode
              i = i + 1
              call getarg(i,tag)
              if (tag.eq.'h') then
                  auto = 1
              elseif (tag.eq.'d') then
                  auto = 2
              elseif (tag.eq.'s') then
                  auto = 3
              else
                  call usage('!! Error: no -auto option '//tag)
              endif
              ! Read distance value
              i = i + 1
              call getarg(i,tag)
              read(tag,'(BN,F10.0)') adist
              ! Read increment
              i = i + 1
              call getarg(i,tag)
              ptr = index(tag,'d')
              if (ptr.ne.0) then
                  tag(ptr:ptr) = ''
                  read(tag,'(BN,F10.0)') incr
              endif
              ptr = index(tag,'k')
              if (ptr.ne.0) then
                  tag(ptr:ptr) = ''
                  read(tag,'(BN,F10.0)') incr
                  incr = -incr
              endif
              if (incr.gt.1.0d98) then
                  read(tag,'(BN,F10.0)') incr
                  if (auto.eq.2.or.auto.eq.3) incr = -incr
              endif
          elseif (tag(1:10).eq.'-empirical'.or.tag(1:4).eq.'-emp') then
              i = i + 1
              call getarg(i,emprel)
          elseif (tag(1:4).eq.'-gmt') then
              kgmt = 1
              i = i + 1
              call getarg(i,gmtf)
          elseif (tag(1:4).eq.'-vol') then
              kvol = 1
              i = i + 1
              call getarg(i,volf)
          elseif (tag(1:8).eq.'-verbose') then
              verbos = 1
              i = i + 1
              if (i.le.narg) then
                  call getarg(i,tag)
                  if (index(tag,'-').eq.0) then
                      read(tag,'(BN,I10)') verbos
                  else
                      i = i - 1
                  endif
              else
                  i = i - 1
              endif
          elseif (tag(1:2).eq.'-h'.or.tag(1:5).eq.'-help') then
              call usagelong(' ')
          elseif (tag(1:2).eq.'-d'.or.tag(1:8).eq.'-details') then
              call usagelong('long')
          else
              call usage('!! Error: no option '//tag)
          endif
          goto 101
  102 continue
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE errmes(text)
      IMPLICIT none
      CHARACTER text*(*)
      write(0,*) text
      STOP
      END

C----------------------------------------------------------------------C

      SUBROUTINE usage(str)
      IMPLICIT NONE
      INTEGER lstr
      CHARACTER str*(*)
      if (str.ne.' ') then
          lstr = len(str)
          write(*,*) str(1:lstr)
          write(*,*) '!! Use -h or -d to see option details'
          write(*,*)
      endif
      write(*,*)
     1 'Usage: o92util -ffm FFMFILE -fsp FSPFILE -flt FLTFILE ',
     2                        '-mag MAGFILE [-fn|-pt]'
      write(*,*)
     1 '               -sta STAFILE -trg TRGFILE|S/D/R/F ',
     2                               '[-haf HAFFILE] [-xy]'
      write(*,*)
     1 '               -disp DSPFILE [-az] -strain STNFILE -stress ',
     2                             'STSFILE -estress ESTSFILE'
      write(*,*)
     1 '               -normal NORFILE -shear[max] SHRFILE ',
     2                      '-coul COULFILE'
      write(*,*)
     1 '               [-auto h|d|s DIST INCR] [-long] [-prog] ',
     2                   '[-gmt GMTFILE]'
      write(*,*)
     1 '               [-empirical WC|MB|B|YM[p]] [-thr THR]'
      write(*,*)
     1 '               [-vol VOLFILE (!!TESTING!!)]'
      write(*,*)
     1 '               [-h|-help] [-d|-details]'
      write(*,*)
      STOP
      END

C----------------------------------------------------------------------C

      SUBROUTINE usagelong(str)
      IMPLICIT NONE
      INTEGER lstr
      CHARACTER str*(*)
      if (str.eq.'long') goto 100
      if (str.ne.' ') then
          lstr = len(str)
          write(*,*) str(1:lstr)
          write(*,*)
      endif
 100  write(*,*)
     1 'Usage: o92util -ffm FFMFILE -fsp FSPFILE -flt FLTFILE ',
     2                        '-mag MAGFILE [-fn|-pt]'
      write(*,*)
     1 '               -sta STAFILE -trg TRGFILE|S/D/R/F ',
     2                               '[-haf HAFFILE] [-xy]'
      write(*,*)
     1 '               -disp DSPFILE [-az] -strain STNFILE -stress ',
     2                             'STSFILE -estress ESTSFILE'
      write(*,*)
     1 '               -normal NORFILE -shear[max] SHRFILE -coul ',
     2                        'COULFILE'
      write(*,*)
     1 '               [-auto h|d|s DIST INCR] [-long] [-prog] ',
     2                   '[-gmt GMTFILE]'
      write(*,*)
     1 '               [-empirical WC|MB|B|YM[p]] [-thr THR]'
      write(*,*)
     1 '               [-vol VOLFILE (!!TESTING!!)]'
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
     2                          'file in static out format'
      if (str.eq.'long') then
          write(*,*)
          write(*,*)
     1 '    See Hayes (2017), The finite, kinematic rupture properties',
     2      ' of great-sized'
          write(*,*)
     1 '        earthquakes since 1990, EPSL.'
          write(*,*)
     1 '    Examples can be found at event pages on earthquake.usgs.gov'
          write(*,*)
          write(*,*)
     1 '    NB: o92util assumes static FFM displacements are in cm'
          write(*,*)
          write(*,9999)
          write(*,*)
      endif
      write(*,*)
     1 '-fsp FSPFILE            Finite fault ',
     2                          'file in SRCMOD FSP format'
      if (str.eq.'long') then
          write(*,*)
          write(*,*)
     1 '    See Hayes (2017), The finite, kinematic rupture properties',
     2      ' of great-sized'
          write(*,*)
     1 '        earthquakes since 1990, EPSL.'
          write(*,*)
     1 '    File format is described at http://equake-rc.info/SRCMOD'
          write(*,*)
     1 '    Examples can be found at event pages on earthquake.usgs.gov'
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
     1 '                                (degrees)          along-dip ',
     2                                                    'along-str'
          write(*,*)
          write(*,*)
     1 '        (or, if -xy is used: evx(km) evy(km) evdp(km)...)'
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
          write(*,*)
     1 '        (or, if -xy is used: evx(km) evy(km) evdp(km)...)'
          write(*,*)
          write(*,*)
     1 '    The magnitude is converted to slip, width, and length ',
     2      'using '
          write (*,*)
     1 '    empirical relations, and the half-space shear modulus'
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
          write(*,*)
     1 '        (or, if -xy is used: stx(km) sty(km) stdp(km))'
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
     1 '    File Format: str dip rak frict'
          write(*,*)
     1 '      If NTRG=1: all locations in STAFILE have same ',
     2                    'target parameters'
          write(*,*)
     1 '      If NTRG=NSTA: target fault parameters ',
     2                     'correspond to station locations'
          write(*,*)
          write(*,*)
     1 '    Command Line Format: str/dip/rak/frict'
          write(*,*)
     1 '      NTRG=1, and all locations have same target parameters'
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
     1 '                   OR'
          write(*,*)
     1 '            "Lame"  lamda(Pa)  shr_mod(Pa)'
          write(*,*)
     1 '                   OR'
          write(*,*)
     1 '            "Poisson"  shr_mod(Pa)  poisson'
          write(*,*)
     1 '                   OR'
          write(*,*)
     1 '            "Young" young(Pa)  poisson'
          write(*,*)
          write(*,*)
     1 '    If seismic velocities are given, elastic parameters are ',
     2     'computed:'
          write(*,*)
     1 '        mu = vs*vs*dens'
          write(*,*)
     1 '        lamda = vp*vp*dens - 2*mu'
          write(*,*)
          write(*,*)
     1 '    If -haf option not specified, o92util defaults to:'
          write(*,*)
     1 '        vp = 6800 m/s; vs = 3926 m/s; dens = 2800 kg/m^3'
          write(*,*)
     1 '        lamda = 4.3e10 Pa; mu = 4.3e10 Pa'
          write(*,*)
          write(*,9999)
          write(*,*)
      endif
      write(*,*)
     1 '-xy                     Treat geographic ',
     2                               'coordinates as ',
     2                          'Cartesian x-y (in km)'
      if (str.eq.'long') then
          write(*,*)
          write(*,*)
     1 '    evlo -> ev_x; evla -> ev_y; stlo -> st_x; stla -> st_y'
          write(*,*)
     1 '            (km)          (km)          (km)          (km)'
          write(*,*)
          write(*,*)
     1 '    All other input formats and units are unchanged'
          write(*,*)
     1 '    Strike angle still defined CW from north/pos. y'
          write(*,*)
          write(*,9999)
          write(*,*)
      endif
      write(*,*)
     1 '-disp DSPFILE           Displacement (E N Z)'
      if (str.eq.'long') then
          write(*,*)
          write(*,*)
     1 '    Format: stlo stla stdp u_E(m) u_N(m) u_Z(m)'
          write(*,*)
     1 '        (or, if -xy is used: stx(km) sty(km) stdp(km)...)'
          write(*,*)
     1 '    Vertical displacement is positive up.'
          write(*,*)
          write(*,9999)
          write(*,*)
      endif
      write(*,*)
     1 '-az                     Displacement output -> (az hdsp Z)'
      if (str.eq.'long') then
          write(*,*)
          write(*,*)
     1 '    Format: stlo stla stdp AZIM U_HOR(m) u_Z(m)'
          write(*,*)
     1 '    Vertical displacement is positive up.'
          write(*,*)
     1 '    Azimuth is clockwise from north.'
          write(*,*)
          write(*,9999)
          write(*,*)
      endif
      write(*,*)
     1 '-strain STRNFILE        Strain matrix (EE NN ZZ EN EZ NZ)'
      if (str.eq.'long') then
          write(*,*)
          write(*,*)
     1 '    Format: stlo stla stdp e_EE e_NN e_ZZ e_EN e_EZ e_NZ'
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
     1 '    Format: stlo stla stdp s_EE(Pa) s_NN(Pa) s_ZZ(Pa) s_EN(Pa)',
     2                         ' s_EZ(Pa) s_NZ(Pa)'
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
     1 '-estress STRSFILE       Effective shear stress (2nd invar. of ',
     2                                'dev. strs tensor)'
      if (str.eq.'long') then
          write(*,*)
          write(*,*)
     1 '    Format: stlo stla stdp estre(Pa)'
          write(*,*)
     1 '    The effective shear stress is a measure of the amplitude ',
     2       'of the stress tensor'
          write(*,*)
     1 '    (deviatoric part).'
          write(*,*)
     1 '        estre = 1/3*(s11*s11 + s22*s22 + s33*s33 -',
     2                      ' s11*s22 - s11*s33 - s22*s33)'
          write(*,*)
     1 '                          + (s12*s12 + s13*s13 + s23*s23)'
          write(*,*)
     1 '        estre = 1/6*((s11-s22)^2 + (s11-s33)^2',
     2                      ' + (s22-s33)^2)'
          write(*,*)
     1 '                          + (s12*s12 + s13*s13 + s23*s23)'
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
     1 '    Format: stlo stla stdp normal(Pa)'
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
     1 '-shear[max] SHEARFILE   Shear stress (',
     2                               'requires target fault kinematics)'
      if (str.eq.'long') then
          write(*,*)
          write(*,*)
     1 '    Format: stlo stla stdp shear(Pa)'
          write(*,*)
     1 '    The reported shear stress is the projection of the shear ',
     2                      'traction'
          write(*,*)
     1 '    onto the rake vector of the target fault.'
          write(*,*)
          write(*,*)
     1 '    To have output be maximum value of shear traction ',
     2       'on target geometry,'
          write(*,*)
     1 '        use -shearmax.'
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
     1 '    Format: stlo stla stdp coulomb(Pa)'
          write(*,*)
     1 '    coulomb = shear + frict*normal'
          write(*,*)
          write(*,9999)
          write(*,*)
      endif
      write(*,*)
     1 '-auto h|d|s DIST INCR   Generate automatic grid'
      if (str.eq.'long') then
          write(*,*)
          write(*,*)
     1 '    Option to automatically generate station file ',
     1       '(autosta.dat) instead of'
          write(*,*)
     1 '    creating station file and using -sta option'
          write(*,*)
          write(*,*)
     1 '    TYPE'
          write(*,*)
     1 '      h: Horizontal grid'
          write(*,*)
     1 '      d: Dip-parallel vertical grid'
          write(*,*)
     1 '      s: Strike-parallel vertical grid'
          write(*,*)
          write(*,*)
     1 '    DIST'
          write(*,*)
     1 '      h: depth (km)'
          write(*,*)
     1 '      d: horizontal, along-strike  offset (km)'
          write(*,*)
     1 '      s: horizontal, down-dip  offset (km)'
          write(*,*)
          write(*,*)
     1 '    INCR'
          write(*,*)
     1 '      Append d: degrees (default for h)'
          write(*,*)
     1 '      Append k: km (default for d, s)'
          write(*,*)
     1 '        e.g. 100, 0.1d, or 15k'
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
     1 '-empirical WC|MB|B|YM   Empirical scaling relation'
      if (str.eq.'long') then
          write(*,*)
          write(*,*)
     1 '    Choose the empirical scale relation (default Wells and ',
     2                'Coppersmith, 1994)'
          write(*,*)
     1 '        WC: Wells and Coppersmith (1994)'
          write(*,*)
     1 '        MB: Mai and Beroza (2000)'
          write(*,*)
     1 '        B:  Blaser et al. (2010)'
          write(*,*)
     1 '        YM: Yen and Ma (2011)'
          write(*,*)
     1 '    If p is appended, then the fault parameters will also be ',
     2                 'printed'
          write(*,*)
     2 '        to standard error'
          write(*,*)
          write(*,9999)
          write(*,*)
      endif
      write(*,*)
     1 '-thr THR                Minimum slip threshold'
      if (str.eq.'long') then
          write(*,*)
          write(*,*)
     1 '    Set slip less than THR to zero'
          write(*,*)
     1 '        THR>0: THR is a slip value in meters'
          write(*,*)
     1 '        -1<THR<0: |THR| is a fraction of maximum slip value'
          write(*,*)
          write(*,9999)
          write(*,*)
      endif
      write(*,*)
     1 '-vol VOLFILE            Volume source file (TESTING!!)'
      write(*,*)
     1 '    NOTE: cannot use -vol with -ffm, -flt, or -mag'
      if (str.eq.'long') then
          write(*,*)
          write(*,*)
     1 '    Format: p evlo evla evdp(km) dvol(m^3)'
          write(*,*)
     1 '                 OR'
          write(*,*)
     1 '    Format: evlo evla evdp(km) str dip disp(m)  wid(km)   ',
     2                                                       'len(km)'
          write(*,*)
     1 '                              (degrees)          along-dip ',
     2                                                    'along-str'
          write(*,*)
     1 '    To define point expansion source, line must start with "p"'
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
      STOP
      END
