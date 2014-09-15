      PROGRAM dateutil
C----
C Utilities for converting to/from calendar dates and days
C----
      IMPLICIT none
      INTEGER yr1,mo1,yr2,mo2
      REAL*8 dy1,dy2,jd1,jd2,ndy
      CHARACTER*30 ifile,ofile,date1,date2
      INTEGER opt,p
      LOGICAL ex

C NOTE TO MATT: PUT IN OPTION (SOMETHING LIKE -hms) TO INCLUDE
C HOUR, MINUTE, SECONDS IN COMPUTATIONS

C----
C Parse command line
C----
      call gcmdln(ifile,ofile,opt,p)
      inquire(file=ifile,EXIST=ex)
      if (.not.ex) then
          write(*,*) '!! Error: No input file '//ifile
          call usage()
      endif

C----
C Open input and output files (unless -p flag is on)
C----
      open (unit=11,file=ifile,status='old')
      if (p.eq.0) open (unit=12,file=ofile,status='unknown')

C----
C Operation differs depending on input option
C   opt=1:
C   opt=2:
C   opt=3:
C   opt=4:
C----
C     Option 1: Convert from date pairs to number of days
      if (opt.eq.1) then
  101     read (11,*,end=102) yr1,mo1,dy1,yr2,mo2,dy2
              call date2jd(jd1,yr1,mo1,dy1)
              call date2jd(jd2,yr2,mo2,dy2)
              ndy = jd2 - jd1
              if (p.eq.0) then
                  write (12,9998) ndy
              else
                  write (*,9998) ndy
              endif
              goto 101
  102     continue
C     Option 2: Convert from date and number of days to date
      elseif (opt.eq.2) then
  103     read(11,*,end=104) yr1,mo1,dy1,ndy
              call date2jd(jd1,yr1,mo1,dy1)
              jd2 = jd1 + ndy
              call jd2date(yr1,mo1,dy1,jd2)
              if (p.eq.0) then
                  write (12,9999) yr1,mo1,dy1
              else
                  write (*,9999) yr1,mo1,dy1
              endif
              goto 103
  104     continue
      elseif (opt.eq.3) then
  105     read(11,*,end=106) yr1,mo1,dy1
              call date2jd(jd1,yr1,mo1,dy1)
              if (p.eq.0) then
                  write(12,9998) jd1
              else
                  write(*,9998) jd1
              endif
              goto 105
  106     continue
      elseif (opt.eq.4) then
  107     read(11,*,end=108) jd1
              call jd2date(yr1,mo1,dy1,jd1)
              if (p.eq.0) then
                  write(12,9999) yr1,mo1,dy1
              else
                  write(*,9999) yr1,mo1,dy1
              endif
              goto 107
  108     continue
      endif

 9998 format(F14.6)
 9999 format(I7,I4,F14.6)

      END

C======================================================================C

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
 
      SUBROUTINE gcmdln(ifile,ofile,opt,p)
 
      IMPLICIT NONE
      CHARACTER*30 ifile,ofile,tag
      INTEGER i,narg,opt,p

      ifile = 'none'
      ofile = 'dateutil.out'
C     1:date->nday; 2:nday->date; 3:cal->jd; 4:jd->cal
      opt = 0
      p = 0

      narg = iargc()
      if (narg.eq.0) call usage()
      i = 0
  201 i = i + 1
      if (i.gt.narg) goto 202
          call getarg(i,tag)
          if (tag(1:2).eq.'-h'.or.tag(1:2).eq.'-?') then
              call usage()
          elseif (tag(1:2).eq.'-o') then
              i = i + 1
              call getarg(i,ofile)
          elseif (tag(1:4).eq.'-day') then
              opt = 1
              i = i + 1
              call getarg(i,ifile)
          elseif (tag(1:5).eq.'-date') then
              opt = 2
              i = i + 1
              call getarg(i,ifile)
          elseif (tag(1:4).eq.'-jd') then
              opt = 3
              i = i + 1
              call getarg(i,ifile)
          elseif (tag(1:4).eq.'-cal') then
              opt = 4
              i = i + 1
              call getarg(i,ifile)
          elseif (tag(1:2).eq.'-p') then
              p = 1
          else
              write(*,*) '!! Error: No option '//tag
              call usage()
          endif
          goto 201
  202 continue

      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE usage()
      IMPLICIT none

      write(*,*)
     1 'Usage: dateutil -day DATEFILE -date DAYFILE ',
     2                 '-jd DATEFILE -cal JDFILE'
      write(*,*)
     3 '                -o OUTFILE -p -h -?'
      write(*,*)
      write (*,*)
     1 'Multiple options for converting between calendar dates, ',
     2 'number of days, and'
      write (*,*)
     1 'Julian calendars.'
      write (*,*)
      write(*,*)
     1 'OPTION         DEFAULT_VALUE      DESCRIPTION'
      write(*,*)
     1 '-day DATEFILE                     compute number of days from ',
     2                                            'date pair'
      write(*,*)
     1 '-date DAYFILE                     compute calendar date from ',
     2                                            'date and day'
      write(*,*)
     1 '-jd CALFILE                       convert calendar date to ',
     2                                          'Julian days'
      write(*,*)
     1 '-cal JDFILE                       convert Julian days to ',
     2                                           'calendar date'
      write(*,*)
     1 '-o OUTFILE     (dateutil.out)     name of output file'
      write (*,*)
     1 '-p                                print results to standard ',
     2                                       'output'
      write (*,*)
     1 '-h                                this online help'
      write (*,*)
     1 '-?                                this online help'
      write (*,*)
      write (*,*)
     1 '----- FILE FORMATS -----'
      write (*,*)
     1 'DATEFILE: list of pairs of calendar dates'
      write (*,*)
     1 '    YYYY1 MM1 DD1 YYYY2 MM2 DD2'
      write (*,*)
     1 'DAYFILE: list of pairs of calendar dates and number of days'
      write (*,*)
     1 '    YYYY MM DD NDAY'
      write (*,*)
     1 'CALFILE: list of calendar dates'
      write (*,*)
     1 '    YYYY MM DD [HH MM SS]'
      write (*,*)
     1 'JDFILE: list of Julian days'
      write (*,*)
     1 '    JD.[  ]'
      write (*,*)
     1 'OUTFILE:'
      write (*,*)
     1 '    (-day)  NDAYS'
      write (*,*)
     1 '    (-date) YYYY MM DD'
      write (*,*)
     1 '    (-jd)   JD.[   ]'
      write (*,*)
     1 '    (-cal)  YYYY MM DD'
      write (*,*) ''

      STOP
      END
