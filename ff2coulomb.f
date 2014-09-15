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
C     trgparams.txt: target fault strike, dip, rake, and
C                             coefficient of friction
C
C Produces the output files:
C     coul.out: stlo, stla, coulomb stress change
C     shear.out: stlo, stla, shear stress on target faults
C     norml.out: stlo, stla, normal stress on fault (+ compressive)
C     strain.out: stlo, stla, eEE, eNN, eZZ, eEN, eEZ, eNZ
C     stress.out: stlo, stla, sEE, sNN, sZZ, sEN, sEZ, sNZ
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

      INTEGER f,nflt,maxflt,trgmatchsta
      PARAMETER (maxflt=1000)
      REAL*8 evlo(maxflt),evla(maxflt),evdp(maxflt),str(maxflt),
     1       dip(maxflt),rak(maxflt),slip(maxflt),dx(maxflt),
     2       dy(maxflt),area(maxflt)

      REAL*8 stlo,stla,stdp,dist,az,x,y
      REAL*8 vp,vs,dens,trgstr,trgdip,trgrak,frict

      INTEGER prog,prog100,progtag,stact,trgct
      INTEGER flttyp,i,j
      REAL*8 e(3,3),enet(3,3),stress(3,3),norml,shear,coul
      CHARACTER*80 ffmfile,stafile,haffile,trgfile,strafile,strefile,
     1             norfile,shrfile,coufile

C----
C Parse the command line
C----
      call gcmdln(ffmfile,stafile,haffile,trgfile,strafile,strefile,
     1            norfile,shrfile,coufile,flttyp)

C----
C Check for input files
C----
      call checkctrlfiles(ffmfile,stafile,haffile,trgfile)
      open (unit=22,file=stafile,status='old')
      open (unit=23,file=haffile,status='old')
      open (unit=24,file=trgfile,status='old')
      read (23,*) vp,vs,dens

C----
C Compare number of lines in stations.txt and targetparams.txt
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
          trgmatchsta = 0
      elseif (trgct.eq.stact) then
          trgmatchsta = 1
      else
          call usage('!! Error: NTRG must be equal to 1 or NSTA')
      endif

C----
C Output files
C----
      open (unit=11,file=coufile,status='unknown')
      open (unit=12,file=shrfile,status='unknown')
      open (unit=13,file=norfile,status='unknown')
      open (unit=14,file=strafile,status='unknown')
      open (unit=15,file=strefile,status='unknown')

C----
C Read finite fault file
C----
      call readstaticout(ffmfile,nflt,evlo,evla,evdp,str,dip,rak,dx,dy,
     1                   area,slip,maxflt)

C----
C Calculate strain, stress, and coulomb stress change at each station
C----  
      prog = 0
      progtag = 0
      call progbar(prog,prog100,progtag)
      if (trgmatchsta.eq.0) read(24,*) trgstr,trgdip,trgrak,frict
  101 read (22,*,end=107) stlo,stla,stdp
          if (trgmatchsta.eq.1) read(24,*) trgstr,trgdip,trgrak,frict
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

          call strain2stress(stress,enet,vp,vs,dens)
          call coulomb(coul,norml,shear,stress,trgstr,trgdip,trgrak,
     1                 frict)
C
          write (11,9999) stlo,stla,coul
          write (12,9999) stlo,stla,shear ! dot prod of shear with rake
          write (13,9999) stlo,stla,norml ! positive => dilation
          write (14,*) stlo,stla,enet(1,1),enet(2,2),enet(3,3),
     1                           enet(1,2),enet(1,3),enet(2,3)
          write (15,*) stlo,stla,stress(1,1),stress(2,2),stress(3,3),
     1                           stress(1,2),stress(1,3),stress(2,3)
          
          prog = prog + 1
          call progbar(prog,prog100,progtag)
          goto 101
  107 continue
C----

 9999 format (2f10.3,f15.1)

      END

C======================================================================C

      SUBROUTINE gcmdln(ffmfile,stafile,haffile,trgfile,strafile,
     1                  strefile,norfile,shrfile,coufile,flttyp)
      IMPLICIT none
      CHARACTER*80 tag
      INTEGER narg,i
      CHARACTER*80 ffmfile,stafile,haffile,trgfile,strafile,strefile,
     1             norfile,shrfile,coufile
      INTEGER flttyp

      ffmfile = 'static_out'
      stafile = 'stations.txt'
      haffile = 'structure.txt'
      trgfile = 'trgparams.txt'
      strafile = 'strain.out'
      strefile = 'stress.out'
      norfile = 'norml.out'
      shrfile = 'shear.out'
      coufile = 'coul.out'
      flttyp = 1

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
          elseif (tag(1:4).eq.'-haf') then
              i = i + 1
              call getarg(i,haffile)
          elseif (tag(1:4).eq.'-trg') then
              i = i + 1
              call getarg(i,trgfile)
          elseif (tag(1:5).eq.'-stra') then
              i = i + 1
              call getarg(i,strafile)
          elseif (tag(1:5).eq.'-stre') then
              i = i + 1
              call getarg(i,strefile)
          elseif (tag(1:4).eq.'-nor') then
              i = i + 1
              call getarg(i,norfile)
          elseif (tag(1:4).eq.'-shr') then
              i = i + 1
              call getarg(i,shrfile)
          elseif (tag(1:5).eq.'-coul') then
              i = i + 1
              call getarg(i,coufile)
          elseif (tag(1:3).eq.'-fn') then
              flttyp = 1
          elseif (tag(1:3).eq.'-pt') then
              flttyp = 0
          else
              call usage('!! Error: no option '//tag)
          endif
          goto 101
  102 continue

      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE readstaticout(ffmfile,nflt,evlo,evla,evdp,str,dip,rak,
     1                         dx,dy,area,slip,maxflt)
      IMPLICIT none
      CHARACTER*1 dum
      INTEGER ct,seg,nseg,nx,ny,i,maxflt,nflt,ptr
      CHARACTER*10 dxc,dyc
      CHARACTER*80 ffmfile
      REAL*8 evlo(maxflt),evla(maxflt),evdp(maxflt),str(maxflt),
     1       dip(maxflt),rak(maxflt),dx(maxflt),dy(maxflt),area(maxflt),
     2       slip(maxflt)
C----
C Read header with number of fault segments
C----
      open (unit=21,file=ffmfile,status='old')
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

C----------------------------------------------------------------------C

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

 1000 format ('ff2coulomb: [',I3,'% Complete]',A)
 
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE checkctrlfiles(ffmfile,stafile,haffile,trgfile)
      IMPLICIT none
      CHARACTER*20 ffmfile,stafile,haffile,trgfile
      LOGICAL ex
      inquire(file=ffmfile,EXIST=ex)
      if (.not.ex) call usage('!! Error: no file named '//ffmfile)
      inquire(file=stafile,EXIST=ex)
      if (.not.ex) call usage('!! Error: no file named '//stafile)
      inquire(file=haffile,EXIST=ex)
      if (.not.ex) call usage('!! Error: no file named '//haffile)
      inquire (file=trgfile,EXIST=ex)
      if (.not.ex) call usage('!! Error: no file named '//trgfile)
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
     1 'Usage: ff2coulomb -ffm FFMFILE ',
     2                   '-sta STAFILE -haf HAFFILE ',
     2                   '-haf HAFFILE -trg TRGFILE ',
     3                   '-stra STRAFILE -stre STREFILE ',
     4                   '-nor NORMLFILE -shr SHEARFILE ',
     5                   '-coul COULFILE ',
     6                   '-fn -pt -AUTO ',
     7                   '-h -?'
      write(STER,*)
      write(STER,*)
     1 'Compute strains and stresses resulting from finite fault model'
      write(STER,*)
      write (STER,*)
     1 '-ffm FFMFILE   (static_out)     name of finite fault model file'
      write (STER,*)
     1 '-sta STAFILE   (stations.txt)   name of receiver location file'
      write (STER,*)
     1 '-haf HAFFILE   (structure.txt)  name of half-space ',
     2                                  'parameter file'
      write (STER,*)
     1 '-trg TRGFILE   (trgparams.txt)  name of target fault file'
      write (STER,*)
     1 '-stra STRAFILE (strain.out)     name of output strain matrix ',
     2                                  'file'
      write (STER,*)
     1 '-stre STREFILE (stress.out)     name of output stress matrix ',
     2                                  'file'
      write (STER,*)
     1 '-nor NORMLFILE (norml.out)      name of output normal stress ',
     2                                  'file'
      write (STER,*)
     1 '-shr SHEARFILE (shear.out)      name of output shear stress ',
     2                                  'file'
      write (STER,*)
     1 '-coul COULFILE (coul.out)       name of output Coulomb stress ',
     2                                  'file'
      write (STER,*)
     1 '-fn            (default)        treat subfaults as ',
     2                                    'finite sources'
      write (STER,*)
     1 '-pt                             treat subfaults as ',
     2                                       'point sources'
      write (STER,*)
     1 '-AUTO'
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
     1 'TRGFILE: target fault parameters'
      write (STER,*)
     1 '    str(deg) dip(deg) rak(deg) coeff_frict'
      write (STER,*)
     1 '    NOTE: If NTRG = 1, all faults have same parameters'
      write (STER,*)
     1 '          If NTRG = NSTA, parameters correspond to location'
      write (STER,*)
     1 'STRAFILE: list of station locations (corresponds to STAFILE) ',
     2           'and strain components at receiver'
      write (STER,*)
     1 '    stlo stla eEE eNN eZZ eEN eEZ eNZ'
      write (STER,*)
     1 'STREFILE: list of station locations (corresponds to STAFILE) ',
     2           'and stress components at receiver'
      write (STER,*)
     1 '    stlo stla sEE(Pa) sNN(Pa) sZZ(Pa) sEN(Pa) sEZ(Pa) sNZ(Pa)'
      write (STER,*)
     1 'NORMLFILE: list of station locations (corresponds to STAFILE) ',
     2           'and normal stresses on receiver'
      write (STER,*)
     1 '    stlo stla norml_stress(Pa)'
      write (STER,*)
     1 '    NOTE: positive => dilation'
      write (STER,*)
     1 'SHEARFILE: list of station locations (corresponds to STAFILE) ',
     2           'and shear stresses on receiver'
      write (STER,*)
     1 '    stlo stla shear_stress(Pa)'
      write (STER,*)
     1 '    NOTE: projection of shear traction along rake direction'
      write (STER,*)
     1 'COULFILE: list of station locations (corresponds to STAFILE) ',
     2           'and Coulomb stresses on receiver'
      write (STER,*)
     1 '    stlo stla coul_stress(Pa)'
      write (STER,*)

      STOP
      END
