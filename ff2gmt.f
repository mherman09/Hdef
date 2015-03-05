      PROGRAM ff2gmt
C----
C Convert the static_out finite fault model to a format that GMT
C can use for mapping using psxy with -SJ option
C----
      IMPLICIT none
      REAL*8 pi,d2r
      PARAMETER (pi=4.0d0*datan(1.0d0),d2r=pi/180.0d0)
      CHARACTER*40 ifile,ofile,nfile,tfile,cfile
      INTEGER getdep
      LOGICAL ex
      INTEGER i,nflt,FMAX
      PARAMETER (FMAX=1500)
      REAL*8 evlo(FMAX),evla(FMAX),evdp(FMAX),str(FMAX),dip(FMAX),
     1       rak(FMAX),dx(FMAX),dy(FMAX),slip(FMAX),trup(FMAX),hylo,hyla
      REAL*8 time
      INTEGER iL,iR,iB,iT
      REAL*8 Lmx,Rmx,Bmx,Tmx

C Parse command line, check for files
      call gcmdln(ifile,ofile,time,nfile,tfile,cfile,getdep)
      if (ifile.eq.'none') then
          write(*,*) '!! Error: Input finite fault file unspecified'
          call usage('!! Use -f FFMFILE to specify input file')
      else
          inquire(file=ifile,EXIST=ex)
          if (.not.ex) call usage('!! Error: no input file: '//ifile)
      endif
      if (ofile.eq.'none') then
          write(*,*) '!! Error: No output file specified'
          call usage('!! Use -o SLPFILE to specify slip output file')
      endif

C Read finite fault model, initiate output files
      call readffm(ifile,evlo,evla,evdp,str,dip,rak,dx,dy,slip,trup,
     1             hylo,hyla,nflt)
      open(unit=10,file=ofile,status='unknown')
      if (nfile.ne.'none') open(unit=11,file=nfile,status='unknown')

C----
C Write files in GMT readable format. If subfault has ruptured, include
C slip so patch can be colored, otherwise leave out so patch will be
C transparent.
C----
      do 14 i = 1,nflt
          if (getdep.eq.1) then
              write(10,*) evlo(i),evla(i),evdp(i)*1.0d-3,str(i),
     1                    dx(i)*1.0d-3,dy(i)*dcos(d2r*dip(i))*1.0d-3
              goto 14
          endif
          if (time.le.0.0d0) then
              write(10,*) evlo(i),evla(i),slip(i),str(i),dx(i)*1.0d-3,
     1                    dy(i)*dcos(d2r*dip(i))*1.0d-3
          else
              if (trup(i).le.time) then
                  write(10,*) evlo(i),evla(i),slip(i),str(i),
     1                        dx(i)*1.0d-3,dy(i)*dcos(d2r*dip(i))*1.0d-3
              elseif (nfile.ne.'none') then
                  write(11,*) evlo(i),evla(i),str(i),dx(i)*1.0d-3,
     1                        dy(i)*dcos(d2r*dip(i))*1.0d-3
              endif
          endif
   14 continue

C Print XYT file
      if (tfile.ne.'none') then
          open(unit=13,file=tfile,status='unknown')
          do 15 i = 1,nflt
              write(13,*) evlo(i),evla(i),trup(i)
   15     continue
      endif

C Print quadrilateral outline of FFM
      if (cfile.ne.'none') then
          open(unit=14,file=cfile,status='unknown')
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
   16     continue
          write(14,*) evlo(iL),evla(iL)
          write(14,*) evlo(iT),evla(iT)
          write(14,*) evlo(iR),evla(iR)
          write(14,*) evlo(iB),evla(iB)
          write(14,*) evlo(iL),evla(iL)
      endif

      END

C======================================================================c

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

      SUBROUTINE gcmdln(ifile,ofile,time,nfile,tfile,cfile,getdep)
      IMPLICIT none
      CHARACTER*40 tag,ifile,ofile,nfile,tfile,cfile
      INTEGER narg,i,getdep
      REAL*8 time
      ifile = 'none'
      ofile = 'none'
      nfile = 'none'
      tfile = 'none'
      cfile = 'none'
      time = -1.0d0
      getdep = 0
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
          elseif (tag(1:2).eq.'-o') then
              i = i + 1
              call getarg(i,ofile)
          elseif (tag(1:5).eq.'-clip') then
              i = i + 1
              call getarg(i,cfile)
          elseif (tag(1:3).eq.'-tf') then
              i = i + 1
              call getarg(i,tfile)
          elseif (tag(1:2).eq.'-t') then
              i = i + 1
              call getarg(i,tag)
              read (tag,'(BN,F12.0)') time
          elseif (tag(1:2).eq.'-z') then
              getdep = 1
          elseif (tag(1:2).eq.'-n') then
              i = i + 1
              call getarg(i,nfile)
          elseif (tag(1:2).eq.'-h') then
              call usage(' ')
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
     1 'Usage: ff2gmt -f FFMFILE -o SLPFILE [-t TIME] [-n NOSLPFILE] ',
     2               '[-tf XYTFILE] [-clip CLIPFILE] [-h] [-z]'
      write(*,*)
      write(*,*)
     1 '-f FFMFILE     Finite fault model in standard subfault format'
      write(*,*)
     1 '-o SLPFILE     File with slip patches for use with "psxy ',
     2                '-SJ -C<cptfile>"'
      write(*,*)
     1 '-t TIME        Only include patches that rupture before TIME ',
     2                 '(default: include all subfaults)'
      write(*,*)
     1 '-n NOSLPFILE   Write patches that slip after TIME for use ',
     2                'with "psxy -SJ"'
      write(*,*)
     1 '-tf XYTFILE    Write lon lat trup to file'
      write(*,*)
     1 '-clip CLIPFILE Determine outlines of FFM and write to file'
      write(*,*)
     1 '-z             Replace slip with depth (km) in third column'
      write(*,*)
     1 '-h             Online help (this screen)'
      STOP
      END
