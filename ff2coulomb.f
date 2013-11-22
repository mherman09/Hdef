      PROGRAM ff2coulomb
C----
C Calculate the stress changes at stations due to subfault shear
C dislocations comprising a finite fault model in an isotropic
C halfspace.
C
C To run, requires the input files:
C     static_out: finite fault model in subfault format
C     stations.txt: list of receiver lon, lat, and depth (km)
C     structure.txt: vp (m/s) vs (m/s) density (kg/m^3)
C     targetparams.txt: target fault strike, dip, rake, and
C                             coefficient of friction
C
C Produces the output files:
C     coul.out: stlo, stla, coulomb stress change
C     shear.out: stlo, stla, shear stress on target faults
C     norml.out: stlo, stla, normal stress on fault (+ compressive)
C     strain.out: stlo, stla, eEE, eNN, eZZ, eEN, eEZ, eNZ
C
C MODIFICATIONS:
C   Fall 2012: Original file created
C   2013-07-23: Subroutine to check that input files exist
C   2013-08-09: Added verbose option
C   2013-08-09: Added help/usage option
C   2013-09-20: Option to define target parameters for each station
C               by making targetparams.txt same length as stations.txt
C
C----
      IMPLICIT NONE
      REAL*8 pi,d2r
      PARAMETER (pi=4.0d0*datan(1.0d0),d2r=pi/1.8d2)

      LOGICAL verbose,trgmatchsta
      INTEGER f,nflt,maxflt
      PARAMETER (maxflt=1000)
      REAL*8 evlo(maxflt),evla(maxflt),evdp(maxflt),str(maxflt),
     1       dip(maxflt),rak(maxflt),slip(maxflt),dx(maxflt),
     2       dy(maxflt),area(maxflt)

      REAL*8 stlo,stla,stdp,dist,az,x,y
      REAL*8 vp,vs,dens,trgstr,trgdip,trgrak,frict

      INTEGER prog,prog100,progtag,stact,trgct
      INTEGER flttyp,i,j
      REAL*8 e(3,3),enet(3,3),stress(3,3),norml,shear,coul
      CHARACTER*20 ffmfile,stafile,haffile,trgfile

      trgmatchsta = .false.

C----
C Read the command line to choose the fault type, and read the finite
C fault model. Default is finite source, choose -p to force point
C source.
C----
      call gcmdln(flttyp,verbose)

C----
C Input files
C----
      call checkctrlfiles(ffmfile,stafile,haffile,trgfile,verbose)
      open (unit=21,file=ffmfile,status='old')
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
          call usage()
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

C----
C Read finite fault file
C----
      call readstaticout(nflt,evlo,evla,evdp,str,dip,rak,dx,dy,area,
     1                   slip,maxflt)

C----
C Calculate strain, stress, and coulomb stress change at each station
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
C
          do 106 f = 1,nflt
              call ddistaz(dist,az,evlo(f),evla(f),stlo,stla)
              dist = dist*6.371d6
              x = dist*( dcos(az-d2r*str(f)))
              y = dist*(-dsin(az-d2r*str(f)))

              if (flttyp.eq.0) then
                  call o92ptstn(e,x,y,stdp,evdp(f),dip(f),rak(f),
     1                          area(f),slip(f),vp,vs,dens)
              else
                  call o92rectstn(e,x,y,stdp,evdp(f),dip(f),rak(f),
     1                            dy(f),dx(f),slip(f),vp,vs,dens)
              endif

              call rotstrain(e,str(f))

              do 105 i = 1,3
                  do 104 j = 1,3
                      enet(i,j) = enet(i,j) + e(i,j)
  104             continue
  105         continue
  106     continue
C
          call strain2stress(stress,enet,vp,vs,dens)
          call coulomb(coul,norml,shear,stress,trgstr,trgdip,trgrak,
     1                 frict)
C
          write (11,9999) stlo,stla,coul
          write (12,9999) stlo,stla,shear ! dot prod of shear with rake
          write (13,9999) stlo,stla,norml ! positive => dilation
          write (14,*) stlo,stla,enet(1,1),enet(2,2),enet(3,3),
     1                           enet(1,2),enet(1,3),enet(2,3)
          
          prog = prog + 1
          call progbar(prog,prog100,progtag)
          goto 101
  107 continue
C----

 9999 format (2f10.3,f15.1)

      END

C======================================================================C

      SUBROUTINE gcmdln(flttyp,verbose)
      IMPLICIT none
      CHARACTER*5 tag
      INTEGER narg,i,flttyp
      LOGICAL verbose

      flttyp = 1
      verbose = .false.

      narg = iargc()
      i = 0
   11 i = i + 1
      if (i.gt.narg) goto 12
          call getarg(i,tag)
          if (tag(1:2).eq.'-p') then
              flttyp = 0
          elseif (tag(1:2).eq.'-f') then
              flttyp = 1
          elseif (tag(1:2).eq.'-V') then
              verbose = .true.
          elseif (tag(1:2).eq.'-h'.or.tag(1:2).eq.'-?') then
              call usage() 
          endif
      goto 11
   12 continue

      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE readstaticout(nflt,evlo,evla,evdp,str,dip,rak,dx,dy,
     1                         area,slip,maxflt)
      IMPLICIT none
      CHARACTER*1 dum
      INTEGER ct,seg,nseg,nx,ny,i,maxflt,nflt,ptr
      CHARACTER*10 dxc,dyc
      REAL*8 evlo(maxflt),evla(maxflt),evdp(maxflt),str(maxflt),
     1       dip(maxflt),rak(maxflt),dx(maxflt),dy(maxflt),area(maxflt),
     2       slip(maxflt)
C----
C Read header with number of fault segments
C----
      ct = 0
      read (21,*) dum,dum,dum,dum,nseg
C----
C Each segment has a fixed line header with information about number
C and size of subfaults.
C----
      do 13 seg = 1,nseg
          read (21,*) dum,dum,dum,dum,nx,dum,dxc,dum,ny,dum,dyc
              ptr = index(dxc,'km')
              dxc(ptr:ptr+1) = ''
              ptr = index(dyc,'km')
              dyc(ptr:ptr+1) = ''
          do 11 i=1,8
              read (21,*) dum
   11     continue
C Read subfault parameters
          do 12 i = 1,nx*ny
              read (21,*) evla(ct+i),evlo(ct+i),evdp(ct+i),slip(ct+i),
     1                    rak(ct+i),str(ct+i),dip(ct+i)
C Convert to SI units
              slip(ct+i) = 1.0d-2*slip(ct+i)
              read (dxc,*) dx(ct+i)
              read (dyc,*) dy(ct+i)
              dx(ct+i) = 1.0d3*dx(ct+i)
              dy(ct+i) = 1.0d3*dy(ct+i)
              area(ct+i) = dx(ct+i)*dy(ct+i)
              evdp(ct+i) = 1.0d3*evdp(ct+i)
   12     continue
          ct = ct + nx*ny
   13 continue
      nflt = ct
 
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

 1000 format ('[',I3,'% Complete]',A)
 
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE checkctrlfiles(ffmfile,stafile,haffile,trgfile,verbose)
      IMPLICIT none
      CHARACTER*20 ffmfile,stafile,haffile,trgfile
      LOGICAL ex,verbose

      ffmfile = 'static_out'
      stafile = 'stations.txt'
      haffile = 'structure.txt'
      trgfile = 'targetparams.txt'
C----
C Finite fault model
C----
   11 if (verbose) write (*,9995,advance='no'),ffmfile
      inquire(file=ffmfile,EXIST=ex)
      if (.not.ex) then
          write (*,8889) ffmfile
          write (*,8885)
          read *,ffmfile
          if (ffmfile.eq.'quit') stop
          if (ffmfile.eq.'help') call usage()
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
          if (stafile.eq.'help') call usage()
          goto 12
      else
          if (verbose) write(*,9999)
      endif
C----
C Half-space
C----
   13 if (verbose) write (*,9997,advance='no'),haffile
      inquire (file=haffile,EXIST=ex)
      if (.not.ex) then
          write (*,8889) haffile
          write (*,8887)
          read *,haffile
          if (haffile.eq.'quit') stop
          if (haffile.eq.'help') call usage()
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
          if (trgfile.eq.'help') call usage()
          goto 14
      else
          if (verbose) write(*,9999)
      endif

 9995 format('Looking for finite fault file: ',A20)
 9996 format('Looking for station file:      ',A20)
 9997 format('Looking for half-space file:   ',A20)
 9998 format('Looking for target fault file: ',A20)
 9999 format('FOUND')
 8885 format('Enter name of finite fault file, "quit", or "help":')
 8886 format('Enter name of station file, "quit", or "help":')
 8887 format('Enter name of half-space file, "quit", or "help":')
 8888 format('Enter name of target fault file, "quit", or "help":')
 8889 format('No file named ',A20)

      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE usage()
      IMPLICIT none

      write (*,*)
     1 'Usage: ff2coulomb -p/-f -V -h/-?'
      write (*,*)
     1 '  -p/-f        (default finite) subfaults are ',
     2                            'treated as point (-p) or finite (-f)'
      write (*,*)
     1 '  -V           (default false) turn on verbose operation'
      write (*,*)
     1 '  -h/-?              help'
      write (*,*) ''
      write (*,*)
     1 '  ff2coulomb calculates strain, normal, shear, and coulomb ',
     2    'stresses from finite fault model at user-defined locations.'
      write (*,*)
     1 '  Output files: strain.out, norml.out (pos. = dilation), ',
     2    'shear.out, coul.out'
      write (*,*)
     1 '  Required input files:'
      write (*,*) ''
      write (*,*)
     1 '    static_out: finite fault model in standard subfault format'
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
      write (*,*)
     1 '    targetparams.txt: str dip rak coeff_frict'
      write (*,*)
     1 '      NTRG = 1: use str/dip/rak/frict for all locations'
      write (*,*)
     1 '      NTRG = NSTA: str/dip/rak/frict correspond to station',
     2        ' locations'
      write (*,*) ''

      STOP
      END
