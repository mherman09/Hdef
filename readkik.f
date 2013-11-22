      PROGRAM readkik
C----
C Program to read in Japanese KiK-net files (strong motion).
C----
      IMPLICIT none
      CHARACTER*80 ifile,ofile,dfile
      CHARACTER*18 fmtc,fmtf,fmti,lbl
      CHARACTER*20 stnm,otimec,rtimec,sfreqc,scalc
      REAL evlo,evla,evdp,stlo,stla,stdp
      REAL mag,freq,dt,dur,n,d,scal,maxacc
      INTEGER dir,i,termout,trap
      INTEGER yr0,mo0,dy0,hr0,mn0,sc0,yr1,mo1,dy1,hr1,mn1,sc1,t1

      INTEGER ct
      REAL a(32000),vel(32000),disp(32000),avg,mx
      LOGICAL ex,opt

C----
C Input the strong motion file as a command line argument
C----
      call gcmdln(ifile,ofile,dfile,termout,trap)
  100 inquire(file=ifile,exist=ex)
      if (.not.ex) then
          print *,ifile,'does not exist'
          call usage()
      endif
      open(unit=11,file=ifile,status='old')

C----
C Parse the header of the strong motion file
C----
      fmtc = '(A18,A)'
      fmtf = '(A18,F9.0)'
      fmti = '(A18,I3)'
      read(11,fmtc) lbl,otimec
      read(11,fmtf) lbl,evla
      read(11,fmtf) lbl,evlo
      read(11,fmtf) lbl,evdp
      read(11,fmtf) lbl,mag
      read(11,fmtc) lbl,stnm
      read(11,fmtf) lbl,stla
      read(11,fmtf) lbl,stlo
      read(11,fmtf) lbl,stdp
      read(11,fmtc) lbl,rtimec
      read(11,fmtc) lbl,sfreqc
      read(11,fmtf) lbl,dur
      read(11,fmti) lbl,dir
      read(11,fmtc) lbl,scalc
      read(11,fmtf) lbl,maxacc
      read(11,*)    lbl
      read(11,*)    lbl

      read(otimec,1001) yr0,mo0,dy0,hr0,mn0,sc0
      read(rtimec,1001) yr1,mo1,dy1,hr1,mn1,sc1
 1001 format(I4,1X,I2,1X,I2,1X,I2,1X,I2,1X,I2)
      t1 = (hr1-hr0)*60*60 + (mn1-mn0)*60 + (sc1-sc0)

      i = index(sfreqc,'Hz')
      sfreqc(i:i+1) = ''
      read(sfreqc,*) freq
      dt = 1.0/freq

      i = index(scalc,'(gal)/')
      scalc(i:i+5) = ' '
      read (scalc,*) n,d
      scal = n/d

C----
C Read in acceleration time series, correct to gals, remove offset
C----
      ct = 0
  101 read(11,*,end=103) a(ct+1),a(ct+2),a(ct+3),a(ct+4),a(ct+5),
     1                   a(ct+6),a(ct+7),a(ct+8)
          do 102 i = 1,8
              a(ct+i) = a(ct+i)*scal
  102     continue
          ct = ct + 8
          goto 101
  103 continue
      print *,'Raw data read and converted to gals.'

      avg = 0.0
      do 104 i = 1,ct
          avg = avg + a(i)
  104 continue
      avg = avg/real(ct)
C      mx = 0.0
      do 105 i = 1,ct
          a(i) = a(i) - avg
C          if (abs(a(i)).gt.mx) mx = abs(a(i))
  105 continue
      print *,'Offset removed.'

      do 106 i = 1,ct
          a(i) = 1.0e-2*a(i)
  106 continue

C----
C Option to print time series to a file
C----
      if (ofile.ne.'none') then
          print *,'Output file: ',ofile
      elseif (termout.eq.1) then
          do 107 i = 1,ct
               write (6,*) real(t1)+dt*real(i-1)-15.0,a(i)
  107     continue
      else
          call usage()
      endif

      open(unit=12,file=ofile,status='unknown')
      write (12,*) 0,0
      do 108 i = 1,ct
          write (12,*) real(t1)+dt*real(i-1)-15.0,a(i)
  108 continue

C----
C Numerically integrate to get velocity time series
C----
      if (trap.eq.1) then
          open (unit=13,file=dfile,status='unknown')
          do 201 i = 1,ct
              if (i.eq.1) then
                  vel(i) = 0.0
              else
                  vel(i) = vel(i-1) + dt*0.5*(a(i-1)+a(i))
              endif
  201     continue

          do 202 i = 1,ct
              if (i.eq.1) then
                  disp(i) = 0.0
              else
                  disp(i) = disp(i-1) + dt*0.5*(vel(i-1)+vel(i))
              endif
  202     continue
      endif
      write (13,*) 'V',0,0
      do 203 i = 1,ct
          write (13,*) 'V',real(t1)+dt*real(i-1)-15.0,vel(i)
  203 continue
      write (13,*) 'D',0,0
      do 204 i = 1,ct
          write (13,*) 'D',real(t1)+dt*real(i-1)-15.0,disp(i)
  204 continue

      END

C----------------------------------------------------------------------C

      SUBROUTINE gcmdln(ifile,ofile,dfile,termout,trap)
      IMPLICIT none
      CHARACTER*5 tag
      CHARACTER*80 ifile,ofile,dfile
      INTEGER i,narg,termout,trap

      ifile = 'none'
      ofile = 'none'
      dfile = 'none'
      termout = 0
      trap = 0

      narg = iargc()
      if (narg.eq.0) call usage()

      i = 0
  101 i = i + 1
      if (i.gt.narg) goto 102
      call getarg(i,tag)
      if (tag(1:2).eq.'-h'.or.tag(1:2).eq.'-?') then
          call usage()
      elseif (tag(1:2).eq.'-f') then
          i = i + 1
          call getarg(i,ifile)
      elseif (tag(1:2).eq.'-o') then
          i = i + 1
          call getarg(i,ofile)
      elseif (tag(1:2).eq.'-p') then
          termout = 1
      elseif (tag(1:2).eq.'-i') then
          trap = 1
          i = i + 1
          call getarg(i,dfile)
      endif
      goto 101
  102 continue

      RETURN
      END
      

      SUBROUTINE usage()
      IMPLICIT none
      write(*,*)
     1 'Usage: readkik -f [IN_FILE] -o [OUT_FILE] -i [DISP_FILE] -p ',
     2                '-h/-?'
      write(*,*)
     1 '  -f [IN_FILE] Strong motion data in KIK format'
      write(*,*)
     1 '  -o [OUT_FILE] For use with GMT command psxy'
      write(*,*)
     1 '      1st column: time from origin in seconds'
      write(*,*)
     1 '      2nd column: acceleration in mgal'
      write(*,*)
     1 '  -i [DISP_FILE] (default no) integrate to displacement'
      write(*,*)
     1 '  -p (default no) print to standard out instead of to file'
      write(*,*)
     1 '  -h/-? print help'
      write(*,*)
     1 ''
      write(*,*)
     1 ''
      STOP
      END
