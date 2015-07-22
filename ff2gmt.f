      PROGRAM MAIN
C----
C Convert a finite fault model to a format that GMT can use
C----
      IMPLICIT none
      REAL*8 pi,d2r
      PARAMETER (pi=4.0d0*datan(1.0d0),d2r=pi/180.0d0)
      CHARACTER*40 ifile,sfile,tfile,zfile,cfile,nfile
      LOGICAL ex
      INTEGER i,nflt,FMAX
      PARAMETER (FMAX=1500)
      REAL*8 evlo(FMAX),evla(FMAX),evdp(FMAX),str(FMAX),dip(FMAX),
     1       rak(FMAX),dx(FMAX),dy(FMAX),slip(FMAX),trup(FMAX),hylo,hyla
      REAL*8 time
      REAL*8 cosd

C----
C Parse command line, check for files
C----
      call gcmdln(ifile,sfile,tfile,zfile,cfile,nfile,time)
      if (ifile.eq.'none') then
          write(*,*) '!! Error: Input finite fault file unspecified'
          call usage('!! Use -f FFMFILE to specify input file')
      else
          inquire(file=ifile,EXIST=ex)
          if (.not.ex) call usage('!! Error: no input file: '//ifile)
      endif

C----
C Open active output files
C----
      if (sfile.eq.'none'.and.tfile.eq.'none'.and.zfile.eq.'none'
     1         .and.cfile.eq.'none'.and.nfile.eq.'none') then
          call usage('!! Error: no output file specified')
      endif
      if (sfile.ne.'none') open(unit=21,file=sfile,status='unknown')
      if (tfile.ne.'none') open(unit=22,file=tfile,status='unknown')
      if (zfile.ne.'none') open(unit=23,file=zfile,status='unknown')
      if (nfile.ne.'none') open(unit=25,file=nfile,status='unknown')

C----
C Read finite fault model
C----
      call readffm(ifile,evlo,evla,evdp,str,dip,rak,dx,dy,slip,trup,
     1             hylo,hyla,nflt)

C----
C Write GMT-ready files
C----
      do 102 i = 1,nflt
          cosd = dcos(d2r*dip(i))
          ! If using -trup option check here
          if (time.gt.0.0d0.and.time.lt.trup(i)) then
              if (nfile.ne.'none') then
                  write(25,*) evlo(i),evla(i),str(i),dx(i)*1d-3,
     1                                                   dy(i)*1d-3*cosd
              endif
              goto 101
          endif
          ! Slip file
          if (sfile.ne.'none') then
              write(21,*) evlo(i),evla(i),slip(i),str(i),dx(i)*1d-3,
     1                                                   dy(i)*1d-3*cosd
          endif
          ! Time file
          if (tfile.ne.'none') then
              write(22,*) evlo(i),evla(i),trup(i),str(i),dx(i)*1d-3,
     1                                                   dy(i)*1d-3*cosd
          endif
          ! Depth file
          if (zfile.ne.'none') then
              write(23,*) evlo(i),evla(i),evdp(i)*1d-3,str(i),
     1                                        dx(i)*1d-3,dy(i)*1d-3*cosd
          endif
  101     continue
  102 continue
      if (cfile.ne.'none') then
          call getclip(cfile,evlo,evla,str,dx,dy,dip,nflt)
      endif

      close(21)
      close(22)
      close(23)
      close(25)

      END

C======================================================================c

      SUBROUTINE getclip(cfile,evlo,evla,str,dx,dy,dip,nflt)
      IMPLICIT none
      REAL*8 pi,d2r
      PARAMETER (pi=4.0d0*datan(1.0d0),d2r=pi/180.0d0)
      CHARACTER*40 cfile
      INTEGER i,nflt,FMAX
      PARAMETER (FMAX=1500)
      REAL*8 evlo(FMAX),evla(FMAX),str(FMAX),dip(FMAX),
     1       dx(FMAX),dy(FMAX)
      INTEGER iL,iR,iB,iT
      REAL*8 Lmx,Rmx,Bmx,Tmx,lo1,la1,az1,az2,d1,d2
      open(unit=24,file=cfile,status='unknown')
C Print quadrilateral outline of FFM
      iL = 0
      iR = 0
      iB = 0
      iT = 0
      do 16 i = 1,nflt
          if (i.eq.1) then
              iL = i
              iR = i
              iB = i
              iT = i
              Lmx = evlo(i)
              Rmx = evlo(i)
              Bmx = evla(i)
              Tmx = evla(i)
          else
              if (evlo(i).lt.Lmx) then
                  Lmx = evlo(i)
                  iL = i
              endif
              if (evlo(i).gt.Rmx) then
                  Rmx = evlo(i)
                  iR = i
              endif
              if (evla(i).lt.Bmx) then
                  Bmx = evla(i)
                  iB = i
              endif
              if (evla(i).gt.Tmx) then
                  Tmx = evla(i)
                  iT = i
              endif
          endif
   16 continue
      if (0.0d0.le.str(1).and.str(1).lt.90.0d0) then
          ! TOP +S,+D
          az1 = str(iT)
          az2 = str(iT)-90.0d0
          d1 = 0.5d0*dx(iT)*1d-3
          d2 = 0.5d0*dy(iT)*1d-3*dcos(dip(iT)*d2r)
          call dlola(lo1,la1,evlo(iT),evla(iT),d1,az1)
          call dlola(evlo(iT),evla(iT),lo1,la1,d2,az2)
          ! RGT +S,-D
          az1 = str(iR)
          az2 = str(iR)+90.0d0
          d1 = 0.5d0*dx(iR)*1d-3
          d2 = 0.5d0*dy(iR)*1d-3*dcos(dip(iR)*d2r)
          call dlola(lo1,la1,evlo(iR),evla(iR),d1,az1)
          call dlola(evlo(iR),evla(iR),lo1,la1,d2,az2)
          ! LFT -S,+D
          az1 = 180.0d0+str(iL)
          az2 = str(iL)-90.0d0
          d1 = 0.5d0*dx(iL)*1d-3
          d2 = 0.5d0*dy(iL)*1d-3*dcos(dip(iL)*d2r)
          call dlola(lo1,la1,evlo(iL),evla(iL),d1,az1)
          call dlola(evlo(iL),evla(iL),lo1,la1,d2,az2)
          ! BOT -S,-D
          az1 = 180.0d0+str(iB)
          az2 = str(iB)+90.0d0
          d1 = 0.5d0*dx(iB)*1d-3
          d2 = 0.5d0*dy(iB)*1d-3*dcos(dip(iB)*d2r)
          call dlola(lo1,la1,evlo(iB),evla(iB),d1,az1)
          call dlola(evlo(iB),evla(iB),lo1,la1,d2,az2)
      elseif (90.0d0.le.str(1).and.str(1).lt.180.0d0) then
          ! TOP -S,+D
          az1 = 180.0d0+str(iT)
          az2 = str(iT)-90.0d0
          d1 = 0.5d0*dx(iT)*1d-3
          d2 = 0.5d0*dy(iT)*1d-3*dcos(dip(iT)*d2r)
          call dlola(lo1,la1,evlo(iT),evla(iT),d1,az1)
          call dlola(evlo(iT),evla(iT),lo1,la1,d2,az2)
          ! RGT +S,+D
          az1 = str(iR)
          az2 = str(iR)-90.0d0
          d1 = 0.5d0*dx(iR)*1d-3
          d2 = 0.5d0*dy(iR)*1d-3*dcos(dip(iR)*d2r)
          call dlola(lo1,la1,evlo(iR),evla(iR),d1,az1)
          call dlola(evlo(iR),evla(iR),lo1,la1,d2,az2)
          ! LFT -S,-D
          az1 = 180.0d0+str(iL)
          az2 = str(iL)+90.0d0
          d1 = 0.5d0*dx(iL)*1d-3
          d2 = 0.5d0*dy(iL)*1d-3*dcos(dip(iL)*d2r)
          call dlola(lo1,la1,evlo(iL),evla(iL),d1,az1)
          call dlola(evlo(iL),evla(iL),lo1,la1,d2,az2)
          ! BOT +S,-D
          az1 = str(iB)
          az2 = str(iB)+90.0d0
          d1 = 0.5d0*dx(iB)*1d-3
          d2 = 0.5d0*dy(iB)*1d-3*dcos(dip(iB)*d2r)
          call dlola(lo1,la1,evlo(iB),evla(iB),d1,az1)
          call dlola(evlo(iB),evla(iB),lo1,la1,d2,az2)
      elseif (180.0d0.le.str(1).and.str(1).lt.270.0d0
     1             .or.-180.0d0.le.str(1).and.str(1).lt.-90.0d0) then
          ! TOP -S,-D
          az1 = 180.0d0+str(iT)
          az2 = str(iT)+90.0d0
          d1 = 0.5d0*dx(iT)*1d-3
          d2 = 0.5d0*dy(iT)*1d-3*dcos(dip(iT)*d2r)
          call dlola(lo1,la1,evlo(iT),evla(iT),d1,az1)
          call dlola(evlo(iT),evla(iT),lo1,la1,d2,az2)
          ! RGT -S,+D
          az1 = 180.0d0+str(iR)
          az2 = str(iR)-90.0d0
          d1 = 0.5d0*dx(iR)*1d-3
          d2 = 0.5d0*dy(iR)*1d-3*dcos(dip(iR)*d2r)
          call dlola(lo1,la1,evlo(iR),evla(iR),d1,az1)
          call dlola(evlo(iR),evla(iR),lo1,la1,d2,az2)
          ! LFT +S,-D
          az1 = str(iL)
          az2 = str(iL)+90.0d0
          d1 = 0.5d0*dx(iL)*1d-3
          d2 = 0.5d0*dy(iL)*1d-3*dcos(dip(iL)*d2r)
          call dlola(lo1,la1,evlo(iL),evla(iL),d1,az1)
          call dlola(evlo(iL),evla(iL),lo1,la1,d2,az2)
          ! BOT +S,+D
          az1 = str(iB)
          az2 = str(iB)-90.0d0
          d1 = 0.5d0*dx(iB)*1d-3
          d2 = 0.5d0*dy(iB)*1d-3*dcos(dip(iB)*d2r)
          call dlola(lo1,la1,evlo(iB),evla(iB),d1,az1)
          call dlola(evlo(iB),evla(iB),lo1,la1,d2,az2)
      elseif (270.0d0.le.str(1).and.str(1).lt.360.0d0
     1               .or.str(1).lt.0.0d0.and.str(1).ge.-90.0d0) then
          ! TOP +S,-D
          az1 = str(iT)
          az2 = str(iT)+90.0d0
          d1 = 0.5d0*dx(iT)*1d-3
          d2 = 0.5d0*dy(iT)*1d-3*dcos(dip(iT)*d2r)
          call dlola(lo1,la1,evlo(iT),evla(iT),d1,az1)
          call dlola(evlo(iT),evla(iT),lo1,la1,d2,az2)
          ! RGT -S,-D
          az1 = 180.0d0+str(iR)
          az2 = str(iR)+90.0d0
          d1 = 0.5d0*dx(iR)*1d-3
          d2 = 0.5d0*dy(iR)*1d-3*dcos(dip(iR)*d2r)
          call dlola(lo1,la1,evlo(iR),evla(iR),d1,az1)
          call dlola(evlo(iR),evla(iR),lo1,la1,d2,az2)
          ! LFT +S,+D
          az1 = str(iL)
          az2 = str(iL)-90.0d0
          d1 = 0.5d0*dx(iL)*1d-3
          d2 = 0.5d0*dy(iL)*1d-3*dcos(dip(iL)*d2r)
          call dlola(lo1,la1,evlo(iL),evla(iL),d1,az1)
          call dlola(evlo(iL),evla(iL),lo1,la1,d2,az2)
          ! BOT -S,+D
          az1 = 180.0d0+str(iB)
          az2 = str(iB)-90.0d0
          d1 = 0.5d0*dx(iB)*1d-3
          d2 = 0.5d0*dy(iB)*1d-3*dcos(dip(iB)*d2r)
          call dlola(lo1,la1,evlo(iB),evla(iB),d1,az1)
          call dlola(evlo(iB),evla(iB),lo1,la1,d2,az2)
      endif
      write(24,*) evlo(iL),evla(iL)
      write(24,*) evlo(iT),evla(iT)
      write(24,*) evlo(iR),evla(iR)
      write(24,*) evlo(iB),evla(iB)
      write(24,*) evlo(iL),evla(iL)
      close(24)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE dlola(lon2,lat2,lon1,lat1,dist,az)
C----
C Subroutine for computing final (lon,lat) from (lon,lat), (dist,az)
C Units: dist (km), az (deg)
C----
      IMPLICIT none
      REAL*8 pi,d2r,r2d
      PARAMETER (pi=4.0d0*atan(1.0d0),d2r=pi/1.8d2,r2d=1.8d2/pi)
      REAL*8 lon1,lat1,lon2,lat2,dist,az
      dist = dist/6371.0d0 ! km -> rad
      az = az*d2r          ! deg -> rad
      lat1 = lat1*d2r      ! deg -> rad
      lon1 = lon1*d2r      ! deg -> rad
      lat2 = dasin(dsin(lat1)*dcos(dist)+dcos(lat1)*dsin(dist)*dcos(az))
      lon2 = lon1 + datan2(dsin(az)*dsin(dist)*dcos(lat1),
     1                           dcos(dist)-dsin(lat1)*dsin(lat2))
      lat2 = lat2*r2d
      lon2 = lon2*r2d
C     Return input values in initial units
      dist = dist*6371.0d0
      az = az*r2d
      lat1 = lat1*r2d
      lon1 = lon1*r2d
      RETURN
      END

C----------------------------------------------------------------------c

      SUBROUTINE readffm(ffmfile,evlo,evla,evdp,str,dip,rak,dx,dy,slip,
     1                   trup,hylo,hyla,nflt)
C----
C Read shear dislocations from FFM in standard subfault format.
C----
      IMPLICIT NONE
      CHARACTER*40 ffmfile,du,dxc,dyc
      INTEGER nflt,FMAX
      PARAMETER (FMAX=1500)
      REAL*8 evlo(FMAX),evla(FMAX),evdp(FMAX),str(FMAX),dip(FMAX),
     1       rak(FMAX),dx(FMAX),dy(FMAX),slip(FMAX),trup(FMAX),hylo,hyla
      INTEGER g,nseg,ct,i,nx,ny,ptr
      REAL*8 dxr,dyr
      ct = nflt
      open (unit=31,file=ffmfile,status='old')
      read (31,*) du,du,du,du,nseg
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
          do 312 i = 1,nx*ny
              read (31,*) evla(ct+i),evlo(ct+i),evdp(ct+i),slip(ct+i),
     1                    rak(ct+i),str(ct+i),dip(ct+i),trup(ct+i)
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

C----------------------------------------------------------------------c

      SUBROUTINE gcmdln(ifile,sfile,tfile,zfile,cfile,nfile,time)
      IMPLICIT none
      CHARACTER*40 tag,ifile,sfile,tfile,zfile,cfile,nfile
      INTEGER narg,i
      REAL*8 time
      ifile = 'none'
      sfile = 'none'
      tfile = 'none'
      zfile = 'none'
      cfile = 'none'
      nfile = 'none'
      time = -1.0d0
      narg = iargc()
      if (narg.eq.0) then
          call usage('!! Error: no command line arguments specified')
      endif
      i = 0
  901 i = i + 1
      if (i.gt.narg) goto 902
          call getarg(i,tag)
          if (tag(1:2).eq.'-f') then
              i = i + 1
              call getarg(i,ifile)
          elseif (tag(1:5).eq.'-slip') then
              i = i + 1
              call getarg(i,sfile)
          elseif (tag(1:5).eq.'-time') then
              i = i + 1
              call getarg(i,tfile)
          elseif (tag(1:4).eq.'-dep') then
              i = i + 1
              call getarg(i,zfile)
          elseif (tag(1:5).eq.'-clip') then
              i = i + 1
              call getarg(i,cfile)
          elseif (tag(1:7).eq.'-noslip') then
              i = i + 1
              call getarg(i,nfile)
          elseif (tag(1:5).eq.'-trup') then
              i = i + 1
              call getarg(i,tag)
              read (tag,'(BN,F12.0)') time
          elseif (tag(1:2).eq.'-h') then
              call usage(' ')
          else
              call usage('!! Error: no option '//tag)
          endif
          goto 901
  902 continue
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE usage(str)
      IMPLICIT none
      INTEGER lstr
      CHARACTER str*(*)
      if (str.ne.' ') then
          lstr = len(str)
          write(*,*) str(1:lstr)
          write(*,*)
      endif
      write(*,*)
     1 'Usage: ff2gmt -f FFMFILE -slip SLIPFILE|-time TIMEFILE',
     2                  '|-dep DEPFILE|-clip CLIPFILE'
      write(*,*)
     1 '              [-trup TIME] [-noslip NOSLIPFILE] [-h]'
      write(*,*)
      write(*,*)
     1 '-f FFMFILE         Finite fault model in standard subfault ',
     2                    'format'
      write(*,*)
     1 '-slip SLIPFILE     Slip value in 3rd column for use with ',
     2                    '"psxy -SJ -C<cptfile>"'
C     1 '-o SLPFILE     File with slip patches for use with "psxy ',
C     2                '-SJ -C<cptfile>"'
      write(*,*)
     1 '-time TIMEFILE     Rupture time in 3rd column for use with ',
     2                    '"psxy -SJ -C<cptfile>"'
C     1 '-tf XYTFILE    Write lon lat trup to file'
      write(*,*)
     1 '-dep DEPFILE       Depth in 3rd column for use with ',
     2                    '"psxy -SJ -C<cptfile>"'
C     1 '-z             Replace slip with depth (km) in third column'
      write(*,*)
     1 '-clip CLIPFILE     Write outline of FFM to a file'
      write(*,*)
     1 '-trup TIME         Only include subfaults that rupture before ',
     2                    'TIME'
C     1 '-t TIME        Only include patches that rupture before TIME ',
C     2                 '(default: include all subfaults)'
      write(*,*)
     1 '-noslip NSLPFILE   Write patches that slip after TIME for use ',
     2                    'with "psxy -SJ"'
C     1 '-n NOSLPFILE   Write patches that slip after TIME for use ',
C     2                'with "psxy -SJ"'
      write(*,*)
     1 '-h                 Online help (this screen)'
      STOP
      END
