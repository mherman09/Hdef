      PROGRAM MAIN
C----
C Convert a finite fault model to a format that GMT can use
C----
      IMPLICIT none
      REAL*8 pi,d2r
      PARAMETER (pi=4.0d0*datan(1.0d0),d2r=pi/180.0d0)
      CHARACTER*40 ifile,sfile,tfile,zfile,cfile,nfile,efile,dfile
      LOGICAL ex
      INTEGER i,nflt,FMAX,cf
      PARAMETER (FMAX=1500)
      REAL*8 evlo(FMAX),evla(FMAX),evdp(FMAX),str(FMAX),dip(FMAX),
     1       rak(FMAX),dx(FMAX),dy(FMAX),slip(FMAX),trup(FMAX),hylo,hyla
      REAL*8 time,cosd
      INTEGER seg(FMAX),nseg

C----
C Parse command line, check for files
C----
      call gcmdln(ifile,sfile,tfile,zfile,cfile,nfile,efile,dfile,
     1            time,cf)
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
      if (efile.ne.'none') open(unit=26,file=efile,status='unknown')
      if (dfile.ne.'none') open(unit=27,file=dfile,status='unknown')

C----
C Read finite fault model
C----
      call readffm(ifile,evlo,evla,evdp,str,dip,rak,dx,dy,slip,trup,
     1             hylo,hyla,seg,nseg,nflt)

C----
C Write output files
C----
      do 102 i = 1,nflt
          cosd = dcos(d2r*dip(i))
          ! If using -trup option check here
          if (time.gt.0.0d0.and.time.lt.trup(i)) then
              ! No-slip file
              if (nfile.ne.'none') then
                  write(25,*) evlo(i),evla(i),str(i),dx(i)*1d-3,
     1                                                   dy(i)*1d-3*cosd
              endif
          else
              ! Slip file
              if (sfile.ne.'none') then
                  write(21,*) evlo(i),evla(i),slip(i),str(i),dx(i)*1d-3,
     1                                                   dy(i)*1d-3*cosd
              endif
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
          ! Double couple file
          if (dfile.ne.'none') then
              write(27,*) evlo(i),evla(i),evdp(i),str(i),dip(i),rak(i),
     1                    slip(i)
          endif
  102 continue
      if (cfile.ne.'none') then
          call getclip(cfile,evlo,evla,str,dx,dy,dip,seg,nseg,nflt,cf)
      endif
      if (efile.ne.'none') then
          write(26,*) hylo,hyla
      endif

      close(21)
      close(22)
      close(23)
      close(25)
      close(26)

      END

C======================================================================c

      SUBROUTINE getclip(cfile,evlo,evla,str,dx,dy,dip,seg,nseg,nflt,cf)
      IMPLICIT none
      REAL*8 pi,d2r
      PARAMETER (pi=4.0d0*datan(1.0d0),d2r=pi/180.0d0)
      CHARACTER*40 cfile
      INTEGER i,nflt,FMAX,cf,nseg
      PARAMETER (FMAX=1500)
      REAL*8 evlo(FMAX),evla(FMAX),str(FMAX),dip(FMAX),
     1       dx(FMAX),dy(FMAX)
      INTEGER iL(nseg),iR(nseg),iB(nseg),iT(nseg),seg(FMAX),ichk,j
      REAL*8 Lmx,Rmx,Bmx,Tmx,lo,la,az1,az2,d1,d2
      open(unit=24,file=cfile,status='unknown')
C Print quadrilateral outline of FFM
      do 401 i = 1,nseg
          iL(i) = 0
          iR(i) = 0
          iB(i) = 0
          iT(i) = 0
  401 continue
      Lmx = 0
      Rmx = 0
      Bmx = 0
      Tmx = 0
      ichk = 1
      do 16 i = 1,nflt
          if (i.eq.1.or.seg(i).eq.ichk) then
              iL(ichk) = i
              iR(ichk) = i
              iB(ichk) = i
              iT(ichk) = i
              Lmx = evlo(i)
              Rmx = evlo(i)
              Bmx = evla(i)
              Tmx = evla(i)
              ichk = ichk + 1
          else
              if (evlo(i).lt.Lmx) then
                  Lmx = evlo(i)
                  iL(ichk-1) = i
              endif
              if (evlo(i).gt.Rmx) then
                  Rmx = evlo(i)
                  iR(ichk-1) = i
              endif
              if (evla(i).lt.Bmx) then
                  Bmx = evla(i)
                  iB(ichk-1) = i
              endif
              if (evla(i).gt.Tmx) then
                  Tmx = evla(i)
                  iT(ichk-1) = i
              endif
          endif
   16 continue
      do 402 i = 1,nseg
          if (cf.eq.2) write(24,'(A)') '>'
          if (0.0d0.le.str(iT(i)).and.str(iT(i)).lt.90.0d0) then
              ! TOP +S,+D
              az1 = str(iT(i))
              az2 = str(iT(i))-90.0d0
              d1 = 0.5d0*dx(iT(i))*1d-3
              d2 = 0.5d0*dy(iT(i))*1d-3*dcos(dip(iT(i))*d2r)
              call dlola(lo,la,evlo(iT(i)),evla(iT(i)),d1,az1)
              call dlola(evlo(iT(i)),evla(iT(i)),lo,la,d2,az2)
              ! RGT +S,-D
              az1 = str(iR(i))
              az2 = str(iR(i))+90.0d0
              d1 = 0.5d0*dx(iR(i))*1d-3
              d2 = 0.5d0*dy(iR(i))*1d-3*dcos(dip(iR(i))*d2r)
              call dlola(lo,la,evlo(iR(i)),evla(iR(i)),d1,az1)
              call dlola(evlo(iR(i)),evla(iR(i)),lo,la,d2,az2)
              ! LFT -S,+D
              az1 = 180.0d0+str(iL(i))
              az2 = str(iL(i))-90.0d0
              d1 = 0.5d0*dx(iL(i))*1d-3
              d2 = 0.5d0*dy(iL(i))*1d-3*dcos(dip(iL(i))*d2r)
              call dlola(lo,la,evlo(iL(i)),evla(iL(i)),d1,az1)
              call dlola(evlo(iL(i)),evla(iL(i)),lo,la,d2,az2)
              ! BOT -S,-D
              az1 = 180.0d0+str(iB(i))
              az2 = str(iB(i))+90.0d0
              d1 = 0.5d0*dx(iB(i))*1d-3
              d2 = 0.5d0*dy(iB(i))*1d-3*dcos(dip(iB(i))*d2r)
              call dlola(lo,la,evlo(iB(i)),evla(iB(i)),d1,az1)
              call dlola(evlo(iB(i)),evla(iB(i)),lo,la,d2,az2)
          elseif (90.0d0.le.str(1).and.str(1).lt.180.0d0) then
              ! TOP -S,+D
              az1 = 180.0d0+str(iT(i))
              az2 = str(iT(i))-90.0d0
              d1 = 0.5d0*dx(iT(i))*1d-3
              d2 = 0.5d0*dy(iT(i))*1d-3*dcos(dip(iT(i))*d2r)
              call dlola(lo,la,evlo(iT(i)),evla(iT(i)),d1,az1)
              call dlola(evlo(iT(i)),evla(iT(i)),lo,la,d2,az2)
              ! RGT +S,+D
              az1 = str(iR(i))
              az2 = str(iR(i))-90.0d0
              d1 = 0.5d0*dx(iR(i))*1d-3
              d2 = 0.5d0*dy(iR(i))*1d-3*dcos(dip(iR(i))*d2r)
              call dlola(lo,la,evlo(iR(i)),evla(iR(i)),d1,az1)
              call dlola(evlo(iR(i)),evla(iR(i)),lo,la,d2,az2)
              ! LFT -S,-D
              az1 = 180.0d0+str(iL(i))
              az2 = str(iL(i))+90.0d0
              d1 = 0.5d0*dx(iL(i))*1d-3
              d2 = 0.5d0*dy(iL(i))*1d-3*dcos(dip(iL(i))*d2r)
              call dlola(lo,la,evlo(iL(i)),evla(iL(i)),d1,az1)
              call dlola(evlo(iL(i)),evla(iL(i)),lo,la,d2,az2)
              ! BOT +S,-D
              az1 = str(iB(i))
              az2 = str(iB(i))+90.0d0
              d1 = 0.5d0*dx(iB(i))*1d-3
              d2 = 0.5d0*dy(iB(i))*1d-3*dcos(dip(iB(i))*d2r)
              call dlola(lo,la,evlo(iB(i)),evla(iB(i)),d1,az1)
              call dlola(evlo(iB(i)),evla(iB(i)),lo,la,d2,az2)
          elseif (180.0d0.le.str(1).and.str(1).lt.270.0d0
     1          .or.-180.0d0.le.str(1).and.str(1).lt.-90.0d0) then
              ! TOP -S,-D
              az1 = 180.0d0+str(iT(i))
              az2 = str(iT(i))+90.0d0
              d1 = 0.5d0*dx(iT(i))*1d-3
              d2 = 0.5d0*dy(iT(i))*1d-3*dcos(dip(iT(i))*d2r)
              call dlola(lo,la,evlo(iT(i)),evla(iT(i)),d1,az1)
              call dlola(evlo(iT(i)),evla(iT(i)),lo,la,d2,az2)
              ! RGT -S,+D
              az1 = 180.0d0+str(iR(i))
              az2 = str(iR(i))-90.0d0
              d1 = 0.5d0*dx(iR(i))*1d-3
              d2 = 0.5d0*dy(iR(i))*1d-3*dcos(dip(iR(i))*d2r)
              call dlola(lo,la,evlo(iR(i)),evla(iR(i)),d1,az1)
              call dlola(evlo(iR(i)),evla(iR(i)),lo,la,d2,az2)
              ! LFT +S,-D
              az1 = str(iL(i))
              az2 = str(iL(i))+90.0d0
              d1 = 0.5d0*dx(iL(i))*1d-3
              d2 = 0.5d0*dy(iL(i))*1d-3*dcos(dip(iL(i))*d2r)
              call dlola(lo,la,evlo(iL(i)),evla(iL(i)),d1,az1)
              call dlola(evlo(iL(i)),evla(iL(i)),lo,la,d2,az2)
              ! BOT +S,+D
              az1 = str(iB(i))
              az2 = str(iB(i))-90.0d0
              d1 = 0.5d0*dx(iB(i))*1d-3
              d2 = 0.5d0*dy(iB(i))*1d-3*dcos(dip(iB(i))*d2r)
              call dlola(lo,la,evlo(iB(i)),evla(iB(i)),d1,az1)
              call dlola(evlo(iB(i)),evla(iB(i)),lo,la,d2,az2)
          elseif (270.0d0.le.str(1).and.str(1).lt.360.0d0
     1                .or.str(1).lt.0.0d0.and.str(1).ge.-90.0d0) then
              ! TOP +S,-D
              az1 = str(iT(i))
              az2 = str(iT(i))+90.0d0
              d1 = 0.5d0*dx(iT(i))*1d-3
              d2 = 0.5d0*dy(iT(i))*1d-3*dcos(dip(iT(i))*d2r)
              call dlola(lo,la,evlo(iT(i)),evla(iT(i)),d1,az1)
              call dlola(evlo(iT(i)),evla(iT(i)),lo,la,d2,az2)
              ! RGT -S,-D
              az1 = 180.0d0+str(iR(i))
              az2 = str(iR(i))+90.0d0
              d1 = 0.5d0*dx(iR(i))*1d-3
              d2 = 0.5d0*dy(iR(i))*1d-3*dcos(dip(iR(i))*d2r)
              call dlola(lo,la,evlo(iR(i)),evla(iR(i)),d1,az1)
              call dlola(evlo(iR(i)),evla(iR(i)),lo,la,d2,az2)
              ! LFT +S,+D
              az1 = str(iL(i))
              az2 = str(iL(i))-90.0d0
              d1 = 0.5d0*dx(iL(i))*1d-3
              d2 = 0.5d0*dy(iL(i))*1d-3*dcos(dip(iL(i))*d2r)
              call dlola(lo,la,evlo(iL(i)),evla(iL(i)),d1,az1)
              call dlola(evlo(iL(i)),evla(iL(i)),lo,la,d2,az2)
              ! BOT -S,+D
              az1 = 180.0d0+str(iB(i))
              az2 = str(iB(i))-90.0d0
              d1 = 0.5d0*dx(iB(i))*1d-3
              d2 = 0.5d0*dy(iB(i))*1d-3*dcos(dip(iB(i))*d2r)
              call dlola(lo,la,evlo(iB(i)),evla(iB(i)),d1,az1)
              call dlola(evlo(iB(i)),evla(iB(i)),lo,la,d2,az2)
          endif
          if (cf.eq.2) then
              write(24,*) evlo(iL(i)),evla(iL(i))
              write(24,*) evlo(iT(i)),evla(iT(i))
              write(24,*) evlo(iR(i)),evla(iR(i))
              write(24,*) evlo(iB(i)),evla(iB(i))
              write(24,*) evlo(iL(i)),evla(iL(i))
          endif
  402 continue
      if (cf.eq.1) then
          if (nseg.eq.1) then
              write(24,*) evlo(iL(1)),evla(iL(1))
              write(24,*) evlo(iT(1)),evla(iT(1))
              write(24,*) evlo(iR(1)),evla(iR(1))
              write(24,*) evlo(iB(1)),evla(iB(1))
              write(24,*) evlo(iL(1)),evla(iL(1))
              goto 410
          endif
          lo = evlo(iL(1))
          j = iL(1)
          do 403 i = 2,nseg
              if (evlo(iL(i)).lt.lo) then
                  lo = evlo(iL(i))
                  j = iL(i)
              endif
  403     continue
          write(24,*) evlo(j),evla(j)
          la = evla(iB(1))
          j = iB(1)
          do 405 i = 2,nseg
              if (evla(iB(i)).lt.la) then
                  la = evla(iB(i))
                  j = iB(i)
              endif
  405     continue
          write(24,*) evlo(j),evla(j)
          lo = evlo(iR(1))
          j = iR(1)
          do 404 i = 2,nseg
              if (evlo(iR(i)).gt.lo) then
                  lo = evlo(iR(i))
                  j = iR(i)
              endif
  404     continue
          write(24,*) evlo(j),evla(j)
          la = evla(iT(1))
          j = iT(1)
          do 406 i = 2,nseg
              if (evla(iT(i)).gt.la) then
                  la = evla(iT(i))
                  j = iT(i)
              endif
  406     continue
          write(24,*) evlo(j),evla(j)
          lo = evlo(iL(1))
          j = iL(1)
          do 407 i = 2,nseg
              if (evlo(iL(i)).lt.lo) then
                  lo = evlo(iL(i))
                  j = iL(i)
              endif
  407     continue
          write(24,*) evlo(j),evla(j)
      endif
  410 continue
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
     1                   trup,hylo,hyla,seg,nseg,nflt)
C----
C Read shear dislocations from FFM in standard subfault format.
C----
      IMPLICIT NONE
      CHARACTER*40 ffmfile,du,dxc,dyc
      INTEGER nflt,FMAX
      PARAMETER (FMAX=1500)
      REAL*8 evlo(FMAX),evla(FMAX),evdp(FMAX),str(FMAX),dip(FMAX),
     1       rak(FMAX),dx(FMAX),dy(FMAX),slip(FMAX),trup(FMAX),hylo,hyla
      INTEGER g,nseg,ct,i,nx,ny,ptr,seg(FMAX)
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
              seg(ct+i) = g
  312     continue
          ct = ct + nx*ny
  313 continue
      nflt = nflt + ct
      close(31)
      RETURN
      END

C----------------------------------------------------------------------c

      SUBROUTINE gcmdln(ifile,sfile,tfile,zfile,cfile,nfile,efile,dfile,
     1                  time,cf)
      IMPLICIT none
      CHARACTER*40 tag,ifile,sfile,tfile,zfile,cfile,nfile,efile,dfile
      INTEGER narg,i,cf
      REAL*8 time
      ifile = 'none'
      sfile = 'none'
      tfile = 'none'
      zfile = 'none'
      cfile = 'none'
      nfile = 'none'
      efile = 'none'
      dfile = 'none'
      time = -1.0d0
      cf = 0
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
          elseif (tag(1:8).eq.'-clipseg') then
              cf = 2
              i = i + 1
              call getarg(i,cfile)
          elseif (tag(1:5).eq.'-clip') then
              cf = 1
              i = i + 1
              call getarg(i,cfile)
          elseif (tag(1:4).eq.'-epi') then
              i = i + 1
              call getarg(i,efile)
          elseif (tag(1:3).eq.'-dc') then
              i = i + 1
              call getarg(i,dfile)
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
     1 'Usage: ff2gmt -f FFMFILE -slip SLIPFILE -time TIMEFILE',
     2                  ' -dep DEPFILE'
      write(*,*)
     1 '              -clip CLIPFILE|-clipseg CLIPFILE -epi EPIFILE ',
     2                 '-dc DCFILE'
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
     1 '-clipseg CLIPFILE  Write outline of each segment in FFM to a ',
     2                    'file'
      write(*,*)
     1 '-epi EPIFILE       Extract earthquake epicenter to file'
      write(*,*)
     1 '-dc  DCFILE        Extract strike, dip, rake, slip'
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
