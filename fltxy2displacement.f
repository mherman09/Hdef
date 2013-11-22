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
C----
      IMPLICIT NONE
      REAL*8 pi,d2r
      PARAMETER (pi=4.0d0*datan(1.0d0),d2r=pi/1.8d2)
      
      CHARACTER*20 srcfile,stafile,haffile
      INTEGER flt,nflt,maxflt
      PARAMETER (maxflt=50)
      INTEGER flttyp(maxflt),i,j
      REAL*8 evx(maxflt),evy(maxflt),evdp(maxflt),str(maxflt),
     1       dip(maxflt),rak(maxflt),slip(maxflt),dx(maxflt),dy(maxflt),
     2       area(maxflt)

      REAL*8 stx,sty,stdp,dist,az,x,y,delx,dely
      REAL*8 vp,vs,dens

      REAL*8 ux,uy,uN,uE,uZ,uNnet,uEnet,uZnet
C----
C Input files
C----
      call checkctrlfiles(srcfile,stafile,haffile)
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
C  x y dep(km) str dip rak slip(m) dy(km) dx(km) flttyp
C Last column indicates fault type:
C   0 for point source
C   1 for finite source
C----
      call readfaults(nflt,evx,evy,evdp,str,dip,rak,slip,dy,dx,area,
     1                flttyp,maxflt)
C----
C Calculate net displacements at each station
C----
  101 read (22,*,end=103) stx,sty,stdp
          stdp = 1.0d3*stdp
          uNnet = 0.0d0
          uEnet = 0.0d0
          uZnet = 0.0d0
          
          do 102 flt = 1,nflt
              delx = stx - evx(flt)
              dely = sty - evy(flt)
              dist = dsqrt(delx*delx+dely*dely)
              az = datan2(delx,dely)
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
  
          write (11,9999) stx,sty,uNnet,uEnet,uZnet
          goto 101
  103 continue

 9999 format (2(f10.3),3(f8.3)) 

      END

C======================================================================C

      SUBROUTINE readfaults(nflt,evx,evy,evdp,str,dip,rak,slip,dy,dx,
     1                      area,flttyp,maxflt)

      IMPLICIT none
      INTEGER flt,nflt,maxflt,flttyp(maxflt)
      REAL*8 evx(maxflt),evy(maxflt),evdp(maxflt),str(maxflt),
     1       dip(maxflt),rak(maxflt),slip(maxflt),dy(maxflt),dx(maxflt),
     2       area(maxflt)
      
      flt = 0
  201 flt = flt + 1
      read (21,*,END=202) evx(flt),evy(flt),evdp(flt),str(flt),
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

      SUBROUTINE checkctrlfiles(srcfile,stafile,haffile)
      IMPLICIT none
      CHARACTER*20 srcfile,stafile,haffile
      LOGICAL ex

      srcfile = 'faults.txt'
      stafile = 'stations.txt'
      haffile = 'structure.txt'
C----
C Fault file
C----
   11 write (*,9996,advance='no'),srcfile
      inquire(file=srcfile,EXIST=ex)
      if (.not.ex) then
          write (*,8889)
          write (*,8886)
          read *,srcfile
          if (srcfile.eq.'quit'.or.srcfile.eq.'q') stop
          goto 11
      else
           write(*,9999)
      endif
C----
C Stations
C----
   12 write (*,9997,advance='no'),stafile
      inquire (file=stafile,EXIST=ex)
      if (.not.ex) then
          write (*,8889)
          write (*,8887)
          read *,stafile
          if (stafile.eq.'quit') stop
          goto 12
      else
           write(*,9999)
      endif
C----
C Halfspace
C----
   13 write (*,9998,advance='no'),haffile
      inquire (file=haffile,EXIST=ex)
      if (.not.ex) then
          write (*,8889)
          write (*,8888)
          read *,haffile
          if (haffile.eq.'quit') stop
          goto 13
      else
           write(*,9999)
      endif

 9996 format('Looking for fault file:      ',A20)
 9997 format('Looking for station file:    ',A20)
 9998 format('Looking for half-space file: ',A20)
 9999 format('FOUND')
 8886 format('Enter name of fault file or "quit":')
 8887 format('Enter name of station file or "quit":')
 8888 format('Enter name of half-space file or "quit":')
 8889 format('NOT FOUND')

      RETURN
      END
