      PROGRAM ff2displacement
C----
C Calculate the north, east and vertical displacements at stations
C due to subfault shear dislocations comprising a finite fault model
C in an isotropic halfspace.
C
C To run, requires the input files:
C     static_out: finite fault model in subfault format
C     stations.txt: list of receiver lon, lat, and depth (km)
C     structure.txt: vp (m/s) vs (m/s) density (kg/m^3)
C     targetparameters.txt: target fault strike, dip, rake, and
C                             coefficient of friction
C
C Produces the output file:
C     disp.out: stlo, stla, uN, uE, uZ
C
C MODIFICATIONS:
C   Summer 2012: Original file created, propagation added
C   2013-07-23: Subroutine to check that input files exist
C   2013-08-09: Added verbose option
C   2013-08-09: Added help/usage option
C
C----
      IMPLICIT NONE
      REAL*8 pi,d2r
      PARAMETER (pi=4.0d0*datan(1.0d0),d2r=pi/1.8d2)

      LOGICAL verbose
      INTEGER flt,nflt,maxflt,ct
      PARAMETER (maxflt=1000)
      REAL*8 evlo(maxflt),evla(maxflt),evdp(maxflt),str(maxflt),
     1       dip(maxflt),rak(maxflt),slip(maxflt),dx(maxflt),
     2       dy(maxflt),area(maxflt),trup(maxflt)

      REAL*8 stlo,stla,stdp,dist,az,x,y
      REAL*8 vp,vs,dens

      REAL*8 time,vel,dt,lag,tt
      INTEGER flttyp

      REAL*8 ux,uy,uN,uE,uZ,uNnet,uEnet,uZnet
      CHARACTER*20 ffmfile,stafile,haffile

C----
C Parse the command line.
C time: all contributions that occur before time will be accounted for;
C       the rest will be discarded (default 10000 seconds)
C vel:  velocity of propagation for static displacement from subfault to
C       receiver (default infinite)
C dt:   to get "instantaneous" displacements, input a positive real
C       (accounts for contributions arriving between time-dt and time;
C       default 0.0); otherwise includes all contributions arriving
C       before time
C flttyp: default is finite source, choose -p to force point source
C----
      call gcmdln(time,vel,dt,flttyp,verbose)

C----
C Input files
C----
      call checkctrlfiles(ffmfile,stafile,haffile,verbose)
      open (unit=21,file=ffmfile,status='old')
      open (unit=22,file=stafile,status='old')
      open (unit=23,file=haffile,status='old')
      read (23,*) vp,vs,dens
C----
C Output files
C----
      open (unit=10,file='disp.out',status='unknown')
      rewind 10

C----
C Read finite fault file
C----
      vel = vel*1.0d3
      call readstaticout(nflt,evlo,evla,evdp,str,dip,rak,dx,dy,area,
     1                   slip,trup,maxflt)
C
C----
C Calculate net displacements at each station
C----
  101 read (22,*,end=103) stlo,stla,stdp
          stdp = 1.0d3*stdp
          uNnet = 0.0d0
          uEnet = 0.0d0
          uZnet = 0.0d0
          
          do 102 flt = 1,nflt
              call ddistaz(dist,az,evlo(flt),evla(flt),stlo,stla)
              dist = dist*6.371d6
              x = dist*( dcos(az-d2r*str(flt)))
              y = dist*(-dsin(az-d2r*str(flt)))
              
              if (vel.le.0.0) then
                  lag = 0.0d0
              else
                  lag = dist/vel
              endif
              tt = trup(flt) + lag

              if (dt.le.0.0) then
                  if (tt.lt.time.and.flttyp.eq.0) then
                      call o92pt(ux,uy,uz,x,y,stdp,evdp(flt),
     1                           dip(flt),rak(flt),area(flt),slip(flt),
     2                           vp,vs,dens)
                  elseif (tt.lt.time.and.flttyp.eq.1) then
                      ct = ct + 1
                      call o92rect(ux,uy,uz,x,y,stdp,evdp(flt),
     2                             dip(flt),rak(flt),dy(flt),dx(flt),
     3                             slip(flt),vp,vs,dens)
                  else
                      ux = 0.0d0
                      uy = 0.0d0
                      uz = 0.0d0
                  endif
              else
                  if (tt.lt.time.and.tt.gt.time-dt.and.
     1                                                 flttyp.eq.0) then
                      call o92pt(ux,uy,uz,x,y,stdp,evdp(flt),
     1                           dip(flt),rak(flt),area(flt),slip(flt),
     2                           vp,vs,dens)
                  elseif (tt.lt.time.and.tt.gt.time-dt.and.
     1                                                 flttyp.eq.1) then
                      call o92rect(ux,uy,uz,x,y,stdp,evdp(flt),
     2                             dip(flt),rak(flt),dy(flt),dx(flt),
     3                             slip(flt),vp,vs,dens)
                  else
                      ux = 0.0d0
                      uy = 0.0d0
                      uz = 0.0d0
                  endif
              endif

              call xy2NE(uN,uE,ux,uy,str(flt))
              
              uNnet = uNnet + uN
              uEnet = uEnet + uE
              uZnet = uZnet + uZ
  102     continue
          
          write (10,9999) stlo,stla,uNnet,uEnet,uZnet
          goto 101
  103 continue

 9999 format (2f10.3,3f10.3)

      END

C======================================================================C

      SUBROUTINE gcmdln(time,vel,dt,flttyp,verbose)
      IMPLICIT none
      CHARACTER*25 tag
      INTEGER narg,i,flttyp
      REAL*8 time,vel,dt,trgstr,trgdip,frict
      LOGICAL verbose

      time = 1.0d6
      vel = -1.0d0
      dt = -1.0d0
      flttyp = 1
      verbose = .false.

      narg = iargc()

      i = 0
   11 i = i + 1
      if (i.gt.narg) goto 12
          call getarg(i,tag)
          if (tag(1:2).eq.'-t'.or.tag(1:5).eq.'-time') then
              i = i + 1
              call getarg(i,tag)
              read (tag,'(BN,F10.0)') time
          elseif (tag(1:2).eq.'-v'.or.tag(1:4).eq.'-vel') then
              i = i + 1
              call getarg(i,tag)
              read (tag,'(BN,F5.0)') vel
          elseif (tag(1:2).eq.'-i'.or.tag(1:5).eq.'-inst') then
              i = i + 1
              call getarg(i,tag)
              read (tag,'(BN,F4.0)') dt
          elseif (tag(1:2).eq.'-p') then
              flttyp = 0
          elseif (tag(1:2).eq.'-f') then
              flttyp = 1
          elseif (tag(1:5).eq.'-verb') then
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
     1                         area,slip,trup,maxflt)
      
      IMPLICIT none
      CHARACTER*1 dum
      INTEGER ct,seg,nseg,nx,ny,i,maxflt,nflt,ptr
      CHARACTER*10 dxc,dyc
      REAL*8 evlo(maxflt),evla(maxflt),evdp(maxflt),str(maxflt),
     1       dip(maxflt),rak(maxflt),dx(maxflt),dy(maxflt),area(maxflt),
     2       slip(maxflt),trup(maxflt)
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
     1                    rak(ct+i),str(ct+i),dip(ct+i),trup(ct+i)
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

      SUBROUTINE checkctrlfiles(ffmfile,stafile,haffile,verbose)
      IMPLICIT none
      CHARACTER*20 ffmfile,stafile,haffile
      LOGICAL ex,verbose

      ffmfile = 'static_out'
      stafile = 'stations.txt'
      haffile = 'structure.txt'
C----
C Finite fault model
C----
   11 if (verbose) write (*,9996,advance='no'),ffmfile
      inquire(file=ffmfile,EXIST=ex)
      if (.not.ex) then
          write (*,8889) ffmfile
          write (*,8886)
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
   12 if (verbose) write (*,9997,advance='no'),stafile
      inquire (file=stafile,EXIST=ex)
      if (.not.ex) then
          write (*,8889) stafile
          write (*,8887)
          read *,stafile
          if (stafile.eq.'quit') stop
          if (stafile.eq.'help') call usage()
          goto 12
      else
          if (verbose) write(*,9999)
      endif
C----
C Halfspace
C----
   13 if (verbose) write (*,9998,advance='no'),haffile
      inquire (file=haffile,EXIST=ex)
      if (.not.ex) then
          write (*,8889) haffile
          write (*,8888)
          read *,haffile
          if (haffile.eq.'quit') stop
          if (haffile.eq.'help') call usage()
          goto 13
      else
          if (verbose) write(*,9999)
      endif

 9996 format('Looking for finite fault file: ',A20)
 9997 format('Looking for station file:      ',A20)
 9998 format('Looking for half-space file:   ',A20)
 9999 format('FOUND')
 8886 format('Enter name of finite fault file, "quit", or "help":')
 8887 format('Enter name of station file, "quit", or "help":')
 8888 format('Enter name of half-space file, "quit", or "help":')
 8889 format('No file named ',A20)

      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE usage()
      IMPLICIT none

      write (*,*)
     1 'Usage: ff2displacement -v/-vel [VEL] -t/-time [TIME] -i/-inst ',
     2 '[DT] -p/-f -verb -h/-?'
      write (*,*)
     1 '  -v/-vel [VEL]   (default -1.0) velocity (km/s) of propagation'
      write (*,*)
     1 '                              negative => infinite velocity'
      write (*,*)
     1 '  -t/-time [TIME] (default 1.0e6 sec) only displacements that ',
     2                                   'arrive before TIME contribute'
      write (*,*)
     1 '  -i/-inst [DT]   (default -1.0) include displacements ',
     2                               'arriving between TIME-DT and TIME'
      write (*,*)
     1 '                              negative => add all ',
     2                               'contributions before TIME'
      write (*,*)
     1 '  -p/-f        (default finite) subfaults are ',
     2                            'treated as point (-p) or finite (-f)'
      write (*,*)
     1 '  -verb           (default false) turn on verbose operation'
      write (*,*)
     1 '  -h/-?              help'
      write (*,*) ''
      write (*,*)
     1 '  ff2displacement calculates NEZ displacements from finite',
     2   ' fault model at user-defined locations.'
      write (*,*)
     1 '  Output file: disp.out'
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

      STOP
      END
