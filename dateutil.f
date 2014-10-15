      PROGRAM dateutil
C----
C Utilities for converting to/from calendar dates and days
C----
      IMPLICIT none
      INTEGER yr1,mo1,yr2,mo2
      REAL*8 dy1,dy2,hr1,mn1,sc1,hr2,mn2,sc2,jd1,jd2,ndy,rem
      CHARACTER*40 ifile,ofile,date1,date2,arg1,arg2
      INTEGER opt,p,long
      LOGICAL ex

C----
C Parse command line and check arguments
C----
      call gcmdln(ifile,ofile,arg1,arg2,opt,long,p)
      if (ifile.eq.'none'.and.arg1.eq.'undef') then
          write(*,*) '!! Error: Input data is not specified'
          call usage('!! Use -nday[f] or -date[f]')
      elseif (arg1.eq.'undef') then
          inquire(file=ifile,EXIST=ex)
          if (.not.ex) call usage('!! Error: no input file: '//ifile)
      endif
      if (ofile.eq.'none'.and.p.eq.0) then
          write(*,*) '!! Error: no output file specified'
          call usage('!! Use -o OFILE')
      endif

C----
C If dates/days defined on command line, compute and print results
C----
      if (arg1.ne.'undef') then
          call readdatec(arg1,yr1,mo1,dy1,long)
          call date2jd(jd1,yr1,mo1,dy1)
          ! DATE1 DATE2 -> NDAY
          if (opt.eq.1) then
              call readdatec(arg2,yr2,mo2,dy2,long)
              call date2jd(jd2,yr2,mo2,dy2)
              ndy = jd2 - jd1
              write (*,9997) ndy
          ! DATE1 NDAY -> DATE2
          elseif (opt.eq.2) then
              read(arg2,'(BN,F10.0)') ndy
              jd2 = jd1 + ndy
              call jd2date(yr2,mo2,dy2,jd2)
              if (long.eq.1) then
                  call dyhrmnsc(dy2,hr2,mn2,sc2)
                  write(*,9999) yr2,mo2,dy2,hr2,mn2,sc2
              else
                  write(*,9998) yr2,mo2,dy2 
              endif
          endif
      else
          open (unit=11,file=ifile,status='old')
          if (p.eq.0) open (unit=12,file=ofile,status='unknown')
          if (opt.eq.1) then
  111         if (long.eq.0) then
                  read(11,*,end=113) yr1,mo1,dy1,yr2,mo2,dy2
              else
                  read(11,*,end=113) yr1,mo1,dy1,hr1,mn1,sc1,
     1                               yr2,mo2,dy2,hr2,mn2,sc2
                  dy1 = dy1 + hr1/2.4d1 + mn1/1.44d3 + sc1/8.64d4
                  dy2 = dy2 + hr2/2.4d1 + mn2/1.44d3 + sc2/8.64d4
              endif
              call date2jd(jd1,yr1,mo1,dy1)
              call date2jd(jd2,yr2,mo2,dy2)
              ndy = jd2 - jd1
              if (p.eq.0) then
                  write (12,9997) ndy
              else
                  write (*,9997) ndy
              endif
              goto 111
          else
  112         if (long.eq.0) then
                  read(11,*,end=113) yr1,mo1,dy1,ndy
              else
                  read(11,*,end=113) yr1,mo1,dy1,hr1,mn1,sc1,ndy
                  dy1 = dy1 + hr1/2.4d1 + mn1/1.44d3 + sc1/8.64d4
              endif
              call date2jd(jd1,yr1,mo1,dy1)
              jd2 = jd1 + ndy
              call jd2date(yr2,mo2,dy2,jd2)
              if (long.eq.1) then
                  call dyhrmnsc(dy2,hr2,mn2,sc2)
                  if (p.eq.0) then
                      write(12,9999) yr2,mo2,dy2,hr2,mn2,sc2
                  else
                      write(*,9999) yr2,mo2,dy2,hr2,mn2,sc2
                  endif
              else
                  if (p.eq.0) then
                      write(12,9998) yr2,mo2,dy2
                  else
                      write(*,9998) yr2,mo2,dy2
                  endif
              endif
              goto 112
          endif
      endif
  113 continue

 9997 format(F14.5)
 9998 format(I6,I4,1F6.1)
 9999 format(I6,I4,4F6.1)

      END

C======================================================================C

      SUBROUTINE readdatec(datec,yr,mo,dy,long)
      IMPLICIT none
      CHARACTER*40 datec
      INTEGER yr,mo,long
      REAL*8 dy,hr,mn,sc
      read(datec(1:4),'(BN,I4)') yr
      read(datec(5:6),'(BN,I2)') mo
      read(datec(7:8),'(BN,F10.0)') dy
      if (long.eq.1) then
          read(datec(9:10),'(BN,F10.0)') hr
          read(datec(11:12),'(BN,F10.0)') mn
          read(datec(13:14),'(BN,F10.0)') sc
          dy = dy + hr/2.4d1 + mn/1.44d3 + sc/8.64d4
      endif
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE dyhrmnsc(dy,hr,mn,sc)
      IMPLICIT none
      REAL*8 rem,dy,hr,mn,sc
      rem = mod(dy,1.0d0)
      dy = dy - rem
      hr = rem*2.4d1
      rem = mod(hr,1.0d0)
      hr = hr - rem
      mn = rem*6.0d1
      rem = mod(mn,1.0d0)
      mn = mn - rem
      sc = rem*6.0d1 - mod(sc,1.0d0)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE date2jd(jd,year,month,day)
C----
C Compute the Julian day from the calendar date
C Assumes all dates are after 1582.
C----
      IMPLICIT none
      INTEGER year,month,yrp,mop
      REAL*8 day,A,B,C,D,jd

      if (year.le.1582) print *,'WARNING: CALCUATION MAY NOT WORK ',
     1                          'CORRECTLY FOR YEARS BEFORE 1582!'

      if (month.eq.1.or.month.eq.2) then
          yrp = year - 1
          mop = month + 12
      else
          yrp = year
          mop = month
      endif
      if (year.lt.1582 .or.
     1    year.eq.1582.and.month.lt.10 .or.
     2    year.eq.1582.and.month.eq.10.and.day.lt.15.0d0) then
          B = 0.0d0
      else
          A = floor(dble(yrp)/100)
          B = 2.0d0 - A + floor(A*0.25d0)
      endif
      if (yrp.lt.0) then
          C = floor(365.25d0*dble(yrp) - 0.75d0)
      else
          C = floor(365.25*dble(yrp))
      endif
      D = floor(30.6001d0*(dble(mop)+1.0d0))
      jd = B+C+D+day+1720994.5

      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE jd2date(year,month,day,jd)
C----
C Compute the calendar date from the Julian day
C----
      IMPLICIT none
      REAL*8 jd,F,I,A,B,C,D,E,G,yr,mo,day
      INTEGER year,month

      jd = jd + 0.5d0
      F = jd - floor(jd)
      I = jd - F
      A = floor((I-1867216.25d0)/36524.25d0)
      if (I.gt.2299160d0) then
          B = I + 1.0d0 + A - floor(A*0.25d0)
      else
          B = I
      endif
      C = B + 1524.0d0
      D = floor((C-122.1d0)/365.25d0)
      E = floor(365.25d0*D)
      G = floor((C-E)/30.6001d0)
      day = C-E+F-floor(30.6001d0*G)
      if (G.lt.13.5d0) then
          mo = G-1.0d0
      else
          mo = G-13.0d0
      endif
      if (mo.gt.2.5d0) then
          yr = D - 4716.0d0
      else
          yr = D - 4715.0d0
      endif
      month = int(mo)
      year = int(yr)

      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE gcmdln(ifile,ofile,arg1,arg2,opt,long,p)
      IMPLICIT NONE
      CHARACTER*40 ifile,ofile,tag,arg1,arg2
      INTEGER i,narg,opt,p,long
      ifile = 'none'
      ofile = 'none'
      arg1 = 'undef'
      arg2 = 'undef'
      opt = 0
      long = 0
      p = 0
      narg = iargc()
      if (narg.eq.0) then
          call usage('!! Error: no command line arguments specified')
      endif
      i = 0
  201 i = i + 1
      if (i.gt.narg) goto 202
          call getarg(i,tag)
          if (tag(1:6).eq.'-ndayf') then
              opt = 1
              i = i + 1
              call getarg(i,ifile)
          elseif (tag(1:5).eq.'-nday') then
              opt = 1
              i = i + 1
              call getarg(i,arg1)
              i = i + 1
              call getarg(i,arg2)
              p = 1
          elseif (tag(1:6).eq.'-datef') then
              opt = 2
              i = i + 1
              call getarg(i,ifile)
          elseif (tag(1:5).eq.'-date') then
              opt = 2
              i = i + 1
              call getarg(i,arg1)
              i = i + 1
              call getarg(i,arg2)
              p = 1
          elseif (tag(1:5).eq.'-long') then
              long = 1
          elseif (tag(1:2).eq.'-o') then
              i = i + 1
              call getarg(i,ofile)
          elseif (tag(1:2).eq.'-p') then
              p = 1
          elseif (tag(1:2).eq.'-h') then
              call usage('')
          else
              call usage('!! Error: No option '//tag)
          endif
          goto 201
  202 continue

      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE usage(str)
      IMPLICIT NONE
      INTEGER lstr
      CHARACTER str*(*)
      if (str.ne.' ') then
          lstr = len(str)
          write(*,*) str(1:lstr)
          write(*,*)
      endif
      write(*,*)
     1 'Usage: dateutil -nday DATE1 DATE2 | -ndayf DATEFILE | -date ',
     2               'DATE NDAY | -datef DAYFILE'
      write(*,*)
     3 '                -o OUTFILE [-long] [-p] [-h]'
      write(*,*)
      write(*,*)
     1 '-nday DATE1 DATE2  Number of days between ',
     2                          'dates (YYYYMMDD[HHMMSS])'
      write(*,*)
     1 '-ndayf DATEFILE    Number of days between dates (YYYY MM DD ',
     2                           '[HH MM SS]) in DATEFILE'
      write(*,*)
     1 '-date DATE NDAY    Calendar date from DATE (YYYYMMDD[HHMMSS]) ',
     2                      'and # of days (NDAY)'
      write(*,*)
     1 '-dated DAYFILE     Calendar date from DATE (YYYY MM DD [HH ',
     2                         'MM SS]) and # of days (NDAY) in DAYFILE'
      write(*,*)
     1 '-long              Read/write all dates as YYYYMMDDHHMMSS or ',
     2                            'YYYY MM DD HH MM SS'
      write(*,*)
     1 '-o OUTFILE         Output of date computation'
      write (*,*)
     1 '-p                 Print dates or number of days to standard ',
     2                                       'output (overrides -o)'
      write (*,*)
     1 '-h                 Online help'
      write(*,*)
      write(*,*)
     1 '  NOTE: Dates on command line must be in YYYYMMDD[HHMMSS] ',
     2               'format, with no spaces.'
      write(*,*)
     1 '        Dates in files are delimited by spaces YYYY MM DD ',
     2                            '[HH MM SS]'
      write (*,*)
      STOP
      END
