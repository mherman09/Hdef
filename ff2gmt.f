      PROGRAM ff2gmt
C----
C Convert the static_out finite fault model to a format that GMT
C can use for mapping using psxy with -SJ option
C----
      IMPLICIT none
      REAL*8 pi,d2r
      PARAMETER (pi=4.0d0*datan(1.0d0),d2r=pi/180.0d0)

      LOGICAL ex,verbose
      CHARACTER*40 ifile
      INTEGER flt,nflt,maxflt
      PARAMETER (maxflt=1000)
      REAL*8 evlo(maxflt),evla(maxflt),str(maxflt),dip(maxflt),
     1       rak(maxflt),U(maxflt),trup(maxflt),dx(maxflt),dy(maxflt)

      REAL*8 time,dt

C----
C Parse command line
C----
      call gcmdln(time,dt,ifile,verbose)

C----
C Use control files to manage input/output:
C   10: ff2rupture.xy:
C   11: ff2norupture.xy
C   20: static_out: subfault finite fault model
C----
      inquire(file=ifile,exist=ex)
      if (ex.and.verbose) print *,ifile,' IS HERE :-)'
      if (.not.ex) then
          print *,'WARNING: No file named ',ifile
          stop
      endif
      open (unit=20,file=ifile,status='old')
      open (unit=10,file='slip.xy',status='unknown')
      open (unit=11,file='noslip.xy',status='unknown')
      rewind 10
      rewind 11

C----
C Parse static_out file
C----
      if (verbose) print *,'Parsing finite fault model ',ifile
      call readstaticout(nflt,evlo,evla,str,dip,rak,dx,dy,U,trup,maxflt)

C----
C Write files in GMT readable format. If subfault has ruptured, include
C slip so patch can be colored, otherwise leave out so patch will be
C transparent.
C----
      do 14 flt = 1,nflt
          if (dt.le.0.0) then
              if (trup(flt).lt.time) then
                  write (10,*) evlo(flt),evla(flt),U(flt),str(flt),
     1                         dx(flt),dy(flt)*cos(d2r*dip(flt))
              else
                  write (11,*) evlo(flt),evla(flt),str(flt),dx(flt),
     1                         dy(flt)*cos(d2r*dip(flt))
              endif
          else
              if (trup(flt).lt.time.and.trup(flt).gt.time-dt) then
                  write (10,*) evlo(flt),evla(flt),U(flt),str(flt),
     1                         dx(flt),dy(flt)*cos(d2r*dip(flt))
              else
                  write (11,*) evlo(flt),evla(flt),str(flt),dx(flt),
     1                         dy(flt)*cos(d2r*dip(flt))
              endif
          endif
   14 continue

      if (verbose) print *,'Output slipped patches to slip.xy'
      if (verbose) print *,'Output unslipped patches to noslip.xy'

      END

C======================================================================c

      SUBROUTINE readstaticout(nflt,evlo,evla,str,dip,rak,dx,dy,slip,
     1                         trup,maxflt)
      IMPLICIT none
      CHARACTER*1 dum
      INTEGER ct,seg,nseg,nx,ny,i,maxflt,nflt,ptr
      CHARACTER*10 dxc,dyc
      REAL*8 evlo(maxflt),evla(maxflt),evdp(maxflt),str(maxflt),
     1       dip(maxflt),rak(maxflt),dx(maxflt),dy(maxflt),area(maxflt),
     2       slip(maxflt),ilat,ilon,idep,trup(maxflt)
C----
      ct = 0
      read (20,*) dum,dum,dum,dum,nseg
      do 13 seg = 1,nseg
          read (20,*) dum,dum,dum,dum,nx,dum,dxc,dum,ny,dum,dyc
              ptr = index(dxc,'km')
              dxc(ptr:ptr+1) = ''
              ptr = index(dyc,'km')
              dyc(ptr:ptr+1) = ''
          do 11 i=1,8
              read (20,*) dum
   11     continue
          do 12 i = 1,nx*ny
              read (20,*) evla(ct+i),evlo(ct+i),dum,slip(ct+i),
     1                    rak(ct+i),str(ct+i),dip(ct+i),trup(ct+i)
              slip(ct+i) = 1.0d-2*slip(ct+i)
              read (dxc,*) dx(ct+i)
              read (dyc,*) dy(ct+i)
   12     continue
          ct = ct + nx*ny
   13 continue
      nflt = ct
 
      RETURN
      END

C----------------------------------------------------------------------c

      SUBROUTINE gcmdln(time,dt,ifile,verbose)
      IMPLICIT none
      CHARACTER*40 tag,ifile
      INTEGER narg,i
      REAL*8 dt,time
      LOGICAL verbose

      ifile = 'static_out'
      time = -1.0
      dt = -1.0
      verbose = .false.

      narg = iargc()
      if (narg.eq.0) print *,'To see options, use -h flag'

      i = 0
   11 i = i + 1
      if (i.gt.narg) goto 12
          call getarg(i,tag)
          if (tag(1:5).eq.'-time') then
              i = i + 1
              call getarg(i,tag)
              read (tag,'(BN,F5.0)') time
          elseif (tag(1:5).eq.'-inst') then
              i = i + 1
              call getarg(i,tag)
              read (tag,'(BN,F4.0)') dt
          elseif (tag(1:2).eq.'-f') then
              i = i + 1
              call getarg(i,ifile)
          elseif (tag(1:2).eq.'-h'.or.tag(1:2).eq.'-?') then
              call usage()
          elseif (tag(1:2).eq.'-v') then
              verbose = .true.
          endif
          goto 11
   12 continue

      if (time.lt.0.0) then
          print *, 'Must input time of interest (from start of ',
     1             'rupture) with -time [TIME]'
          call usage()
      endif

      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE usage()
      IMPLICIT none

      write(*,*)
     1 'Usage: ff2gmt -time [TIME] -inst [DT] -f [IFILE] -h/-?'
      write(*,*)
     1 '  -time [TIME] subfaults rupturing before ',
     2         'TIME put in slip.xy, others in noslip.xy'
      write(*,*)
     1 '  -inst [DT] (default -1.0 s) subfaults rupturing between ',
     2         'TIME-DT and TIME put in slip.xy'
      write(*,*)
     1 '  -f [IFILE] (Default static_out) name of input file'
      write (*,*)
     1 '  -h/-?        help'
      write (*,*) ''
      write (*,*)
     1 '  ff2gmt converts standard finite fault model to two files ',
     2 'use with GMT tool psxy and -SJ option'
      write (*,*) ''
      write (*,*)
     1 '    slip.xy: lon lat slip strike length width (hor. proj.)'
      write (*,*)
     1 '                     (m)          (km)  (km)           '
      write (*,*)
     1 '    noslip.xy: lon lat strike length width (hor. proj.)'
      write (*,*)
     1 '                              (km)   (km)              '
      write (*,*) ''
      write (*,*)
     1 '  input file: standard finite fault format'
      write (*,*) ''

      STOP
      END
