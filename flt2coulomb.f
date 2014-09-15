      PROGRAM flt2coulomb
c----
c Calculate the stress changes at stations due to shear dislocations in
c an isotropic halfspace.
c
c To run, requires the following files:
c     faults.txt: list of ruptures
C     stations.txt: list of receiver lon, lat, and depth (km)
C     structure.txt: vp (m/s) vs (m/s) density (kg/m^3)
C     targetparams.txt: target fault strike, dip, rake, and
C                             coefficient of friction
c
C Produces the output files:
C     coul.out: stlo, stla, coulomb stress change
C     shear.out: stlo, stla, shear stress on target faults
C     norml.out: stlo, stla, normal stress on fault (+ compressive)
C     strain.out: stlo, stla, eEE, eNN, eZZ, eEN, eEZ, eNZ
C     stress.out: stlo, stla, sEE, sNN, sZZ, sEN, sEZ, sNZ
C
C MODIFICATIONS:
C   2013-02-27: Original file created
C   2013-07-22: Subroutine to check that input files exist
C   2013-08-08: Added verbose option
C   2013-08-08: Added help/usage option
C
C----
      IMPLICIT NONE
      REAL*8 pi,d2r
      PARAMETER (pi=4.0d0*datan(1.0d0),d2r=pi/1.8d2)

      LOGICAL verbose,trgmatchsta
      CHARACTER*20 srcfile,stafile,haffile,trgfile
      INTEGER flt,nflt,maxflt
      PARAMETER (maxflt=1000)
      INTEGER flttyp(maxflt),i,j
      REAL*8 evlo(maxflt),evla(maxflt),evdp(maxflt),str(maxflt),
     1       dip(maxflt),rak(maxflt),slip(maxflt),dx(maxflt),dy(maxflt),
     2       area(maxflt)

      INTEGER prog,prog100,progtag,stact,trgct

      REAL*8 stlo,stla,stdp,dist,az,x,y
      REAL*8 vp,vs,dens,trgstr,trgdip,trgrak,frict

      REAL*8 e(3,3),enet(3,3),stress(3,3),norml,shear,coul,n(3)

      call gcmdln(verbose)

      trgmatchsta = .false.

C----
C Input files
C----
      call checkctrlfiles(srcfile,stafile,haffile,trgfile,verbose)
      open (unit=21,file=srcfile,status='old')
      open (unit=22,file=stafile,status='old')
      open (unit=23,file=haffile,status='old')
      open (unit=24,file=trgfile,status='old')
      read (23,*) vp,vs,dens

C----
C Compare number of lines in stations.txt and targetparams.txt
C and use line count of stations.txt to initialize progress indicator
C----
      stact = 0
  112 read (22,*,end=113)
          stact = stact + 1
          goto 112
  113 continue
      trgct = 0
  114 read (24,*,end=115)
           trgct = trgct + 1
           goto 114
  115 continue
      rewind 22
      rewind 24

      prog100 = stact
      if (trgct.eq.1) then
          if (verbose) write (*,8997)
      elseif (trgct.eq.stact) then
          if (verbose) write (*,8998)
          trgmatchsta = .true.
      else
          write (*,8999)
          call usage(' ')
      endif
 8997 format ('NL in targetparams.txt = 1')
 8998 format ('NL in targetparams.txt = NL in stations.txt')
 8999 format ('Error: NL in targetparams.txt must be 1 or',
     1        ' equal to NL in stations.txt')

C----
C Output files
C----
      open (unit=11,file='coul.out',status='unknown')
      open (unit=12,file='shear.out',status='unknown')
      open (unit=13,file='norml.out',status='unknown')
      open (unit=14,file='strain.out',status='unknown')
      open (unit=15,file='stress.out',status='unknown')
C----
C Parse fault file
C List of:
C  lon lat dep(km) str dip rak slip(m) dy(km) dx(km) flttyp
C Last column indicates fault type:
C   0 for point source
C   1 for finite source
C----
      call readfaults(nflt,evlo,evla,evdp,str,dip,rak,slip,dy,dx,area,
     1                flttyp,maxflt)
C----
C Calculate net strain, stress, coulomb stress change at each station
C----
      prog = 0
      progtag = 0
      call progbar(prog,prog100,progtag)
      if (.not.trgmatchsta) read(24,*) trgstr,trgdip,trgrak,frict
  101 read (22,*,end=107) stlo,stla,stdp
          if (trgmatchsta) read(24,*) trgstr,trgdip,trgrak,frict
          stdp = 1.0d3*stdp
          do 103 i = 1,3
              do 102 j = 1,3
                  enet(i,j) = 0.0d0
  102         continue
  103     continue
   
          do 106 flt = 1,nflt
              call ddistaz(dist,az,evlo(flt),evla(flt),stlo,stla)
              dist = dist*6.371d6
              x = dist*dcos(az-d2r*str(flt))
              y = dist*(-dsin(az-d2r*str(flt)))

              if (flttyp(flt).eq.0) then
                  call o92ptstn(e,x,y,stdp,evdp(flt),dip(flt),rak(flt),
     1                          area(flt),slip(flt),vp,vs,dens)
              else
                  call o92rectstn(e,x,y,stdp,evdp(flt),dip(flt),
     1                            rak(flt),dy(flt),dx(flt),slip(flt),
     2                            vp,vs,dens)
              endif

              call rotstrain(e,str(flt))

              do 105 i = 1,3
                  do 104 j = 1,3
                      enet(i,j) = enet(i,j) + e(i,j)
  104             continue
  105         continue
  106     continue
  
          call strain2stress(stress,enet,vp,vs,dens)
          call coulomb(coul,norml,shear,stress,trgstr,trgdip,trgrak,
     1                 frict,n)

          write (11,9999) stlo,stla,coul
          write (12,9999) stlo,stla,shear
          write (13,9999) stlo,stla,norml
          write (14,*) stlo,stla,enet(1,1),enet(2,2),enet(3,3),
     1                           enet(1,2),enet(1,3),enet(2,3)
          write (15,*) stlo,stla,stress(1,1),stress(2,2),stress(3,3),
     1                           stress(1,2),stress(1,3),stress(2,3)

          prog = prog + 1
          call progbar(prog,prog100,progtag)
          goto 101
  107 continue

 9999 format (2f10.3,f20.1)
      END

C======================================================================C

      SUBROUTINE readfaults(nflt,evlo,evla,evdp,str,dip,rak,slip,dy,dx,
     1                      area,flttyp,maxflt)

      IMPLICIT none
      INTEGER flt,nflt,maxflt,flttyp(maxflt)
      REAL*8 evlo(maxflt),evla(maxflt),evdp(maxflt),str(maxflt),
     1       dip(maxflt),rak(maxflt),slip(maxflt),dy(maxflt),dx(maxflt),
     2       area(maxflt)
      
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

C----------------------------------------------------------------------C

      SUBROUTINE matmult(matout,mat1,mat2)
      IMPLICIT none
      REAL*8 mat1(3,3),mat2(3,3),matout(3,3)
      INTEGER i,j,k

      DO 1003 i = 1,3
          DO 1002 j = 1,3
              matout(i,j) = 0.0d0
              DO 1001 k = 1,3
                  matout(i,j) = matout(i,j) + mat1(i,k)*mat2(k,j)
 1001         CONTINUE
 1002     CONTINUE
 1003 CONTINUE

      RETURN
      END

      SUBROUTINE mattr(matout,matin)
      IMPLICIT none
      REAL*8 matout(3,3),matin(3,3)
      INTEGER i,j

      DO 1002 i = 1,3
          DO 1001 j = 1,3
              matout(i,j) = matin(j,i)
 1001     CONTINUE
 1002 CONTINUE

      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE progbar(prog,prog100,progtag)
      IMPLICIT none
      INTEGER prog,prog100,progtag
      CHARACTER*1 CR

      CR = char(13)
      if (100*prog/prog100.ge.progtag) then
          write (*,1000,advance='no') 100*prog/prog100,CR
          progtag = progtag + 10
      endif

      if (100*prog/prog100.ge.100) write (*,1000) 100*prog/prog100

 1000 format ('flt2coulomb: [',I3,'% Complete]',A)

      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE checkctrlfiles(srcfile,stafile,haffile,trgfile,verbose)
      IMPLICIT none
      CHARACTER*20 srcfile,stafile,haffile,trgfile
      LOGICAL ex,verbose

      srcfile = 'faults.txt'
      stafile = 'stations.txt'
      haffile = 'structure.txt'
      trgfile = 'targetparams.txt'
C----
C Finite fault model
C----
   11 if (verbose) write (*,9995,advance='no'),srcfile
      inquire(file=srcfile,EXIST=ex)
      if (.not.ex) then
          write (*,8889) srcfile
          write (*,8885)
          read *,srcfile
          if (srcfile.eq.'quit') stop
          if (srcfile.eq.'help') call usage(' ')
          goto 11
      else
           if (verbose) write(*,9999)
      endif
C----
C Stations
C----
   12 if (verbose) write (*,9996,advance='no'),stafile
      inquire (file=stafile,EXIST=ex)
      if (.not.ex) then
          write (*,8889) stafile
          write (*,8886)
          read *,stafile
          if (stafile.eq.'quit') stop
          if (stafile.eq.'help') call usage(' ')
          goto 12
      else
          if (verbose) write(*,9999)
      endif
C----
C Halfspace
C----
   13 if (verbose) write (*,9997,advance='no'),haffile
      inquire (file=haffile,EXIST=ex)
      if (.not.ex) then
          write (*,8889) haffile
          write (*,8887)
          read *,haffile
          if (haffile.eq.'quit') stop
          if (haffile.eq.'help') call usage(' ')
          goto 13
      else
           if (verbose) write(*,9999)
      endif
C----
C Target faults
C----
   14 if (verbose) write (*,9998,advance='no'),trgfile
      inquire (file=trgfile,EXIST=ex)
      if (.not.ex) then
          write (*,8889) trgfile
          write (*,8888)
          read *,trgfile
          if (trgfile.eq.'quit') stop
          if (trgfile.eq.'help') call usage(' ')
          goto 14
      else
          if (verbose) write(*,9999)
      endif

 9995 format('Looking for fault file:        ',A20)
 9996 format('Looking for station file:      ',A20)
 9997 format('Looking for half-space file:   ',A20)
 9998 format('Looking for target fault file: ',A20)
 9999 format('FOUND')
 8885 format('Enter name of fault file, "quit", or "help":')
 8886 format('Enter name of station file, "quit", or "help":')
 8887 format('Enter name of half-space file, "quit", or "help":')
 8888 format('Enter name of target fault file, "quit", or "help":')
 8889 format('No file named ',A20)

      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE gcmdln(verbose)

      IMPLICIT none
      INTEGER i,narg
      CHARACTER*25 tag
      LOGICAL verbose

      verbose = .false.

      narg = iargc()
      i = 0
  101 i = i + 1
      if (i.gt.narg) goto 102
      call getarg(i,tag)
      if (tag(1:2).eq.'-V') then
          verbose = .true.
      elseif (tag(1:2).eq.'-h'.or.tag(1:2).eq.'-?') then
          call usage(' ')
      endif
      goto 101
  102 continue
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE usage(str)
      IMPLICIT none
      INTEGER STO,lstr
      PARAMETER (STO=6)
      CHARACTER str*(*)

      if (str.ne.' ') then
          lstr = len(str)
          write(STO,*)
          write(STO,*) str(1:lstr)
          write(STO,*)
      endif
      write (STO,*)
     1 'Usage: flt2coulomb -flt FLTFILE -mt MTFILE -sta STAFILE ',
     2                    '-haf HAFFILE -trg TRGFILE '
      write (STO,*)
     3 '                   -stra STRAFILE -stre STREFILE ',
     4                    '-nor NORMLFILE -shr SHEARFILE '
      write (STO,*)
     1 '                   -coul COULFILE -fn -pt -auto -gmt GMTFILE ',
     4                        '-h -?'
      write(STO,*)
      write(STO,*)
     1 'Compute strains and stresses displacements resulting from ',
     2 'dislocations'
      write(STO,*)
     1 'Note: the "-mt MTFILE" option overrides "-flt FLTFILE"'
      write(STO,*)
      write(STO,*)
     1 'OPTION         DEFAULT_VALUE      DESCRIPTION'
      write (STO,*)
     1 '-flt FLTFILE   (faults.txt)       name of fault file'
      write (STO,*)
     1 '-mt  MTFILE                       name of fault file (psmeca ',
     2                                    '-Sa format)'
      write (STO,*)
     1 '-sta STAFILE   (stations.txt)     name of receiver location ',
     2                                            'file'
      write (STO,*)
     1 '-haf HAFFILE   (structure.txt)    name of half-space ',
     2                                            'parameter file'
      write (STO,*)
     1 '-dsp DSPFILE   (disp.out)         name of output displacement ',
     2                                  'file'
      write (STO,*)
     1 '-fn            (default)          treat subfaults as ',
     2                                    'finite sources'
      write (STO,*)
     1 '-pt                               treat subfaults as ',
     2                                       'point sources'
      write (STO,*)
     1 '-auto                             automatically create ',
     2                               'receiver grid'
      write (STO,*)
     1 '-gmt GMTFILE                      create file for use with ',
     2                                         'GMT "psxy -SJ"'
      write (STO,*)
     1 '-h                                this online help'
      write (STO,*)
     1 '-?                                this online help'
      write (STO,*)
      write (STO,*)
     1 '----- FILE FORMATS -----'
      write (STO,*)
     1 'FFMFILE: finite fault model in standard subfault format'
      write (STO,*)
     1 'STAFILE: list of station/receiver locations and depths'
      write (STO,*)
     1 '    stlo stla stdp(km)'
      write (STO,*)
     1 'HAFFILE: half-space parameters'
      write (STO,*)
     1 '    vp(km/s) vs(km/s) dens(kg/m^3)'
      write (STO,*)
     1 'DSPFILE: list of station locations (STAFILE coordinates) and ',
     2           'displacements'
      write (STO,*)
     1 '    stlo stla uE(m) uN(m) uZ(m)'
      write (STO,*)
     1 'GMTFILE: list of subfault patches for use with "psxy -SJ"'
      write (STO,*)
     1 '    evlo evla slip(m) str len(km) hor_wid(km)'
      write (STO,*)
      write (STO,*)
     1 '----- NOTES -----'
      write (STO,*)
     1 '1. If "-auto" is selected, a receiver grid is automatically ',
     2       'generated,'
      write (STO,*)
     1 '   and surface displacements are computed at each site.'
      write (STO,*)
     1 '2. A file "autogrid.dat" contains information about the grid.'
      write (STO,*)

      STOP
      END
