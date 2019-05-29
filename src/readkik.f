      PROGRAM readkik
C----
C Program to read in Japanese KiK-net files (strong motion)
C and output acceleration time series, and optional numerical
C integration to output velocity or displacement time series
C----
      IMPLICIT none
      CHARACTER*80 kikfile,accfile,velfile,dspfile
      LOGICAL ex

      INTEGER n,ntot,NMAX,nstart,nfinish
      PARAMETER (NMAX=100000)
      REAL*8 acc(NMAX),vel(NMAX),dsp(NMAX),tstart,dt,first

C----
C Parse command line
C----
      call gcmdln(kikfile,accfile,velfile,dspfile,first)
      inquire(file=kikfile,exist=ex)
      if (.not.ex) then
          write(*,9999) kikfile
          write(*,*)
          call usage()
      endif
      open(unit=11,file=kikfile,status='old')
 9999 format('Strong motion file ',A20,'does not exist')

      if (accfile.eq.'none'.and.velfile.eq.'none'.and.
     1                                        dspfile.eq.'none') then
          print *,'No output file designated!'
          call usage()
      endif

C----
C Parse the strong motion file
C----
      call readsmfile(kikfile,acc,ntot,tstart,dt)
      if (first.gt.0.0d0) then
          nstart = 1
          nfinish = int((first-tstart)/dt) + nstart
          call removemean(acc,ntot,nstart,nfinish)
      endif

C----
C Write acceleration time series
C----
      if (accfile.ne.'none') then
          open (unit=21,file=accfile,status='unknown')
          do 101 n = 1,ntot
              write(21,9998) tstart+dble(n-1)*dt,acc(n)
  101     continue
      endif
C----
C Integrate to get velocity and displacement time series
C----
      if (velfile.ne.'none'.or.dspfile.ne.'none') then
          call numint(vel,acc,dt,ntot)
      endif
      if (velfile.ne.'none') then
          open (unit=22,file=velfile,status='unknown')
          do 102 n = 1,ntot
              write(22,9998) tstart+dble(n-1)*dt,vel(n)
  102     continue
      endif
      if (dspfile.ne.'none') then
          call numint(dsp,vel,dt,ntot)
          open (unit=23,file=dspfile,status='unknown')
          do 103 n = 1,ntot
              write(23,9998) tstart+dble(n-1)*dt,dsp(n)
  103     continue
      endif
 9998 format(F10.3,E18.6)

      END

C======================================================================C

      SUBROUTINE readsmfile(kikfile,acc,ntot,tstart,dt)
      IMPLICIT none
      CHARACTER*80 kikfile
C----
C Parameters read from kikfile header
C----
      CHARACTER*20 lbl,otimec,rtimec,stnm,sampfreqc,scalefactc
      REAL*8 evla,evlo,evdp,mag,stla,stlo,stel,dur,maxacc
      INTEGER cmp
C----
C Parameters to compute
C----
      REAL*8 sampfreq,dt,scaletop,scalebot,scalefact
      REAL*8 yr0,mo0,dy0,hr0,mn0,sc0,yr1,mo1,dy1,hr1,mn1,sc1
      REAL*8 tstart
      INTEGER i

      INTEGER n,ntot,NMAX
      PARAMETER (NMAX=100000)
      REAL*8 acc(NMAX),avg

      open (unit=11,file=kikfile,status='old')
      read(11,1001) lbl,otimec
      read(11,1002) lbl,evla
      read(11,1002) lbl,evlo
      read(11,1002) lbl,evdp
      read(11,1002) lbl,mag
      read(11,1001) lbl,stnm
      read(11,1002) lbl,stla
      read(11,1002) lbl,stlo
      read(11,1002) lbl,stel
      read(11,1001) lbl,rtimec
      read(11,1001) lbl,sampfreqc
      read(11,1002) lbl,dur
      read(11,1003) lbl,cmp
      read(11,1001) lbl,scalefactc
      read(11,1002) lbl,maxacc
      read(11,*)    lbl
      read(11,*)    lbl
 1001 format(A18,A)
 1002 format(A18,F9.0)
 1003 format(A18,I3)

C----
C Start time of recording (tstart, in seconds after origin time)
C Note: Recording start time includes 15 second trigger delay
C       in data logger. The true start time is obtained by
C       subtracting 15 seconds from this "Recording start time"
C       (www.kyoshin.bosai.go.jp/kyoshin/man/knetform_en.html)
C----
      read(otimec,1111) yr0,mo0,dy0,hr0,mn0,sc0
      read(rtimec,1111) yr1,mo1,dy1,hr1,mn1,sc1
C      write(*,1112) yr0,mo0,dy0,hr0,mn0,sc0
C      write(*,1112) yr1,mo1,dy1,hr1,mn1,sc1
 1111 format(F4.0,1X,F2.0,1X,F2.0,1X,F2.0,1X,F2.0,1X,F2.0)
C 1112 format(F5.0,1X,F3.0,1X,F3.0,1X,F3.0,1X,F3.0,1X,F3.0)
      tstart = (hr1-hr0)*3.6d3 + (mn1-mn0)*6.0d1 + (sc1-sc0)
      tstart = tstart - 15.0d0

C----
C Calculate sampling interval (dt) and scale factor (scalefact)
C----
      i = index(sampfreqc,'Hz')
      sampfreqc(i:i+1) = ''
      read(sampfreqc,*) sampfreq
      dt = 1.0d0/sampfreq

      i = index(scalefactc,'(gal)/')
      scalefactc(i:i+5) = ' '
      read (scalefactc,*) scaletop,scalebot
      scalefact = scaletop/scalebot

C----
C Read in acceleration time series
C----
      avg = 0.0d0
      n = 0
  101 read(11,*,end=103) acc(n+1),acc(n+2),acc(n+3),acc(n+4),
     1                   acc(n+5),acc(n+6),acc(n+7),acc(n+8)
          do 102 i = 1,8
              acc(n+i) = acc(n+i)*scalefact
              avg = avg + acc(n+i)
  102     continue
          n = n + 8
          goto 101
  103 continue
      ntot = n
      avg = avg/dble(ntot)
C----
C Remove DC offset and convert to m/s^2 (1 gal = 0.01 m/s^2)
C----
      do 104 n = 1,ntot
          acc(n) = acc(n) - avg
          acc(n) = 1.0d-2*acc(n)
  104 continue

      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE numint(oser,iser,dt,npts)

      INTEGER n,npts,NMAX
      PARAMETER (NMAX=100000)
      REAL*8 oser(NMAX),iser(NMAX),dt

      do 101 n = 1,npts
          if (n.eq.1) then
              oser(n) = 0.0d0
          else
              oser(n) = oser(n-1) + dt*0.5d0*(iser(n-1)+iser(n))
          endif
  101 continue

      RETURN
      END
C----------------------------------------------------------------------C

      SUBROUTINE removemean(tseries,npts,nstart,nfinish)
      IMPLICIT none
      INTEGER i,npts,NMAX,nstart,nfinish
      PARAMETER (NMAX=100000)
      REAL*8 tseries(NMAX),avg

      avg = 0.0d0
      do 101 i = nstart,nfinish
          avg = avg + tseries(i)
  101 continue
      avg = avg/dble(nfinish-nstart+1)
      do 102 i = 1,npts
          tseries(i) = tseries(i) - avg
  102 continue

      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE gcmdln(kikfile,accfile,velfile,dspfile,first)
      IMPLICIT none
      CHARACTER*30 tag
      CHARACTER*80 kikfile,accfile,velfile,dspfile
      INTEGER i,narg !,termout,trap
      REAL*8 first

      kikfile = 'none'
      accfile = 'none'
      velfile = 'none'
      dspfile = 'none'
      first = -1.0d0

      narg = iargc()
      if (narg.eq.0) call usage()
      i = 0
  101 i = i + 1
      if (i.gt.narg) goto 102
          call getarg(i,tag)
          if (tag(1:2).eq.'-h'.or.tag(1:2).eq.'-?') then
              call usage()
          elseif (tag(1:4).eq.'-acc') then
              i = i + 1
              call getarg(i,accfile)
          elseif (tag(1:4).eq.'-vel') then
              i = i + 1
              call getarg(i,velfile)
          elseif (tag(1:4).eq.'-dsp') then
              i = i + 1
              call getarg(i,dspfile)
          elseif (tag(1:6).eq.'-first') then
              i = i + 1
              call getarg(i,tag)
              read(tag,'(BN,F12.0)') first
          elseif (tag(1:2).eq.'-f') then
              i = i + 1
              call getarg(i,kikfile)
          endif
          goto 101
  102 continue

      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE usage()
      IMPLICIT none
      write(*,*)
     1 'readkik -f [KIKFILE] -acc [ACCFILE] -vel [VELFILE] ',
     2                '-dsp [DSPFILE] -first [FIRST_ARR] -h/-?'
      write(*,*)
      write(*,*)
     1 '  Parse a strong motion file in KiK-net format into a two-',
     2    'column time series.'
      write(*,*)
     1 '  Options to integrate numerically to velocity and ',
     2    'displacement, and'
      write(*,*)
     1 '  remove mean from segment of time series (default ',
     2    'removes mean from entire time series).'
      write(*,*)
      write(*,*)
     1 '  -f [KIKFILE]   Strong motion data in KIK format (required)'
      write(*,*)
     1 '  -acc [ACCFILE] Write acceleration (in mgal) time series ',
     2                  'to ACCFILE'
      write(*,*)
     1 '  -vel [VELFILE] Integrate to velocity (in m/s) and write ',
     2                  'time series to VELFILE'
      write(*,*)
     1 '  -dsp [DSPFILE] Integrate to displacement (in m) and write ',
     2                  'time series to DSPFILE'
      write(*,*)
     1 ' -first [FIRST_ARR] (Default -1.0 s) Remove mean between ',
     2     'first point and FIRST_ARR time'
      write(*,*)
     1 '  -h/-?          online help (this screen)'
      write(*,*)
     1 ''
      STOP
      END
