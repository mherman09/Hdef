      PROGRAM flt2displacement
C----
C Calculate the north, east and vertical displacements at stations
C due to shear dislocations in an isotropic halfspace.
C
C To run, requires the following files:
C     faults.txt: list of ruptures
C     stations.txt: list of receiver lon, lat, and depth (km)
C     structure.txt: vp (m/s) vs (m/s) density (kg/m^3)
C
C Produces the output file:
C     disp.out: stlo, stla, uN, uE, uZ
C
C MODIFICATIONS:
C   2013-02-21: Original file created
C   2013-07-22: Subroutine to check that input files exist
C   2013-08-08: Added verbose option
C   2013-08-08: Added help/usage option
C
C----
      IMPLICIT NONE
      REAL*8 pi,d2r
      PARAMETER (pi=4.0d0*datan(1.0d0),d2r=pi/1.8d2)

      LOGICAL verbose
      CHARACTER*20 srcfile,stafile,haffile
      INTEGER flt,nflt,maxflt
      PARAMETER (maxflt=50)
      INTEGER flttyp(maxflt),i,j
      REAL*8 evlo(maxflt),evla(maxflt),evdp(maxflt),str(maxflt),
     1       dip(maxflt),rak(maxflt),slip(maxflt),dx(maxflt),dy(maxflt),
     2       area(maxflt)

      REAL*8 stlo,stla,stdp,dist,az,x,y
      REAL*8 vp,vs,dens

      REAL*8 ux,uy,uN,uE,uZ,uNnet,uEnet,uZnet

      call gcmdln(verbose)
C----
C Input files
C----
      call checkctrlfiles(srcfile,stafile,haffile,verbose)
      open (unit=21,file=srcfile,status='old')
      open (unit=22,file=stafile,status='old')
      open (unit=23,file=haffile,status='old')
      read (23,*) vp,vs,dens
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
     1                flttyp,maxflt)
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
              x = dist*dcos(az-d2r*str(flt))
              y = dist*(-dsin(az-d2r*str(flt)))

              if (flttyp(flt).eq.0) then
                  call o92pt(ux,uy,uz,x,y,stdp,evdp(flt),dip(flt),
     1                       rak(flt),area(flt),slip(flt),vp,vs,dens)
              else
                  call o92rect(ux,uy,uz,x,y,stdp,evdp(flt),dip(flt),
     1                         rak(flt),dy(flt),dx(flt),slip(flt),
     2                         vp,vs,dens)
              endif

              call xy2NE(uN,uE,ux,uy,str(flt))

              uNnet = uNnet + uN
              uEnet = uEnet + uE
              uZnet = uZnet + uZ
  102     continue
  
          write (11,9999) stlo,stla,uNnet,uEnet,uZnet
          goto 101
  103 continue

 9999 format (2(f10.3),3(f12.3)) 

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

      SUBROUTINE checkctrlfiles(srcfile,stafile,haffile,verbose)
      IMPLICIT none
      CHARACTER*20 srcfile,stafile,haffile
      LOGICAL ex,verbose

      srcfile = 'faults.txt'
      stafile = 'stations.txt'
      haffile = 'structure.txt'
C----
C Fault file
C----
   11 if (verbose) write (*,9996,advance='no'),srcfile
      inquire(file=srcfile,EXIST=ex)
      if (.not.ex) then
          write (*,8889) srcfile
          write (*,8886)
          read *,srcfile
          if (srcfile.eq.'quit') stop
          if (srcfile.eq.'help') call usage()
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

 9996 format('Looking for fault file:      ',A20)
 9997 format('Looking for station file:    ',A20)
 9998 format('Looking for half-space file: ',A20)
 9999 format('FOUND')
 8886 format('Enter name of fault file, "quit", or "help":')
 8887 format('Enter name of station file, "quit", or "help":')
 8888 format('Enter name of half-space file, "quit" or "help":')
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
          call usage()
      endif
      goto 101
  102 continue
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE usage()
      IMPLICIT none
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
