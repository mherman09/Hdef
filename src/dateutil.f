      PROGRAM dateutil
C----
C Utilities for converting to/from calendar dates and days
C----
      IMPLICIT none
      INTEGER yr1,mo1,yr2,mo2
      REAL*8 dy1,dy2,ndy
      COMMON /DVARS/ yr1,mo1,yr2,mo2,dy1,dy2,ndy
      CHARACTER*40 ifile,ofile
      CHARACTER*240 arg,ans
      INTEGER opt,p,long,c,ifmt

C----
C Parse command line
C----
      call gcmdln(ifile,ofile,opt,c,arg,long,p)

C----
C Check that all necessary inputs/outputs defined
C----
      call chkopt(ifile,ofile,opt,p,c)

C----
C Calculate date/day
C----
      if (p.eq.0) open(unit=12,file=ofile,status='unknown')
      if (c.eq.0) then
          open(unit=11,file=ifile,status='old')
  102     read(11,'(A240)',end=101) arg
              call getfmt(arg,ifmt)
              call rdarg(arg,ifmt,opt,long)
              call getans(ans,opt,long)
              if (p.eq.0) then
                  write(12,'(A80)') ans
              else
                  write(*,'(A80)') ans
              endif
              goto 102
  101     continue
      else
          call getfmt(arg,ifmt)
          call rdarg(arg,ifmt,opt,long)
          call getans(ans,opt,long)
          write(*,*) trim(ans)
      endif
C
      END

C======================================================================C

      SUBROUTINE getfmt(arg,ifmt)
C----
C Format types:
C   ifmt=1: YYYYMMDD[HHMMSS]
C   ifmt=2: YYYY MM DD [HH MM SS]
C   ifmt=3: YYYY-MM-DDT[HH:MM:SS]
C---- 
      IMPLICIT none
      CHARACTER*240 arg
      CHARACTER*40 dum
      INTEGER ifmt,ptr
      read(arg,*) dum ! Read arg into dum with "*" to ignore leading whitespace
      ifmt = 1
      ptr = index(dum,' ')
      if (ptr.le.5) ifmt = 2
      ptr = index(dum,'T')
      if (ptr.gt.0) ifmt = 3
      ptr = index(dum,'-')
      if (ptr.gt.0) ifmt = 3
      close(41)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE rdarg(arg,ifmt,opt,long)
      IMPLICIT none
      CHARACTER*240 arg,arg1,arg2
      INTEGER opt,long,ifmt
      INTEGER yr1,mo1,yr2,mo2
      REAL*8 dy1,dy2,ndy
      REAL*8 hr1,mn1,sc1,hr2,mn2,sc2
      COMMON /DVARS/ yr1,mo1,yr2,mo2,dy1,dy2,ndy
      !print *,'IFMT',ifmt
      if (ifmt.eq.1) then
          read(arg,*) arg1,arg2
          read(arg1(1:4),*) yr1
          read(arg1(5:6),*) mo1
          read(arg1(7:8),*) dy1
          if (long.eq.1) then
              read(arg1(9:10),*) hr1
              read(arg1(11:12),*) mn1
              read(arg1(13:14),*) sc1
          endif
          if (opt.eq.1) then
              read(arg2(1:4),*) yr2
              read(arg2(5:6),*) mo2
              read(arg2(7:8),*) dy2
              if (long.eq.1) then
                  read(arg2(9:10),*) hr2
                  read(arg2(11:12),*) mn2
                  read(arg2(13:14),*) sc2
              endif
          else
              read(arg2,*) ndy
          endif
      elseif (ifmt.eq.2) then
          if (long.eq.0.and.opt.eq.1) then
              read(arg,*) yr1,mo1,dy1,yr2,mo2,dy2
          elseif (long.eq.0.and.opt.eq.2) then
              read(arg,*) yr1,mo1,dy1,ndy
          elseif (long.eq.1.and.opt.eq.1) then
              read(arg,*,end=2001) yr1,mo1,dy1,hr1,mn1,sc1,
     1                             yr2,mo2,dy2,hr2,mn2,sc2
          elseif (long.eq.1.and.opt.eq.2) then
              read(arg,*,end=2001) yr1,mo1,dy1,hr1,mn1,sc1,ndy
          endif
      elseif (ifmt.eq.3) then
          read(arg,*) arg1,arg2
          read(arg1(1:4),*) yr1
          read(arg1(6:7),*) mo1
          read(arg1(9:10),*) dy1
          if (long.eq.1) then
              read(arg1(12:13),*) hr1
              read(arg1(15:16),*) mn1
              read(arg1(18:19),*) sc1
          endif
          if (opt.eq.1) then
              read(arg2(1:4),*) yr2
              read(arg2(6:7),*) mo2
              read(arg2(9:10),*) dy2    
              if (long.eq.1) then
                  read(arg2(12:13),*) hr2
                  read(arg2(15:16),*) mn2
                  read(arg2(18:19),*) sc2
              endif
          else
              read(arg2,*) ndy
          endif
      endif
      if (opt.eq.1.and.long.eq.1) then
          dy1 = dy1 + hr1/2.4d1 + mn1/1.44d3 + sc1/8.64d4
          dy2 = dy2 + hr2/2.4d1 + mn2/1.44d3 + sc2/8.64d4
      endif
      if (opt.eq.2.and.long.eq.1) then
          dy1 = dy1 + hr1/2.4d1 + mn1/1.44d3 + sc1/8.64d4
      endif
      RETURN
 2001 write(0,*) '!! Error: read long format inputs incorrectly'
      call usage('!! Check input formatting')
      END

C----------------------------------------------------------------------C

      SUBROUTINE getans(ans,opt,long)
C----
C Compute date or number of days from input
C----
      IMPLICIT none
      CHARACTER*240 ans
      INTEGER opt,long
      INTEGER yr1,mo1,yr2,mo2
      REAL*8 dy1,dy2,ndy
      COMMON /DVARS/ yr1,mo1,yr2,mo2,dy1,dy2,ndy
      REAL*8 jd1,jd2,hr,mn,sc
      call date2jd(jd1,yr1,mo1,dy1)
      if (opt.eq.1) then
          call date2jd(jd2,yr2,mo2,dy2)
          ndy = jd2 - jd1
          write(ans,9997) ndy
      else
          jd2 = jd1 + ndy
          call jd2date(yr2,mo2,dy2,jd2)
          if (long.eq.0) then
              write(ans,9998) yr2,mo2,dy2
          elseif (long.eq.1) then
              call dyhrmnsc(dy2,hr,mn,sc)
              if (sc.gt.59.9d0) then
                  jd2 = jd1 + ndy + 1.157d-6
                  call jd2date(yr2,mo2,dy2,jd2)
                  call dyhrmnsc(dy2,hr,mn,sc)
              endif
              write(ans,9999) yr2,mo2,dy2,hr,mn,sc
          endif
      endif
 9997 format(F14.5)
 9998 format(I6,I4,1F6.1)
 9999 format(I6,I4,4F6.1)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE dyhrmnsc(dy,hr,mn,sc)
C----
C Convert decimal day to day, hour, minute, second
C----
      IMPLICIT none
      REAL*8 rem,dy,hr,mn,sc
      rem = dmod(dy,1.0d0)
      dy = dy - rem
      hr = rem*2.4d1
      rem = dmod(hr,1.0d0)
      hr = hr - rem
      mn = rem*6.0d1
      rem = dmod(mn,1.0d0)
      mn = mn - rem
      sc = rem*6.0d1 - dmod(sc,1.0d0)
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

      SUBROUTINE chkopt(ifile,ofile,opt,p,c)
      IMPLICIT none
      CHARACTER*40 ifile,ofile
      INTEGER opt,p,c
      LOGICAL ex
      if (opt.eq.0) then
          write(0,*) '!! Error: Type of date computation unspecified'
          call usage('!! Use -nday or -date')
      endif
      if (ifile.eq.'none'.and.c.eq.0) then
          write(0,*) '!! Error: Input data is unspecified'
          call usage('!! Use -f IFILE or -c ARGS')
      elseif (c.eq.0) then
          inquire(file=ifile,EXIST=ex)
          if (.not.ex) call usage('!! Error: no input file: '//ifile)
      endif
      if (ofile.eq.'none'.and.p.eq.0) then
          write(0,*) '!! Error: no output specified'
          call usage('!! Use -o OFILE or -p')
      endif
      RETURN
      END

C----------------------------------------------------------------------C
C
      SUBROUTINE gcmdln(ifile,ofile,opt,c,arg,long,p)
      IMPLICIT NONE
      CHARACTER*40 ifile,ofile,tag
      CHARACTER*240 arg
      INTEGER i,narg,j,opt,p,c,long
      ifile = 'none'
      ofile = 'none'
      opt = 0
      long = 0
      p = 0
      c = 0
      arg = 'none'
      narg = iargc()
      if (narg.eq.0) then
          call usage('!! Error: no command line arguments specified')
      endif
      i = 0
  201 i = i + 1
      if (i.gt.narg) goto 202
          call getarg(i,tag)
          if (tag(1:5).eq.'-nday') then
              opt = 1
          elseif (tag(1:5).eq.'-date') then
              opt = 2
          elseif (tag(1:5).eq.'-long') then
              long = 1
          elseif (tag(1:2).eq.'-c') then
              c = 1
              p = 1
              j = 0
  231         j = j + 1
              if (i+j.gt.narg) then
                  goto 202
              else
                  call getarg(i+j,tag)
                  if (tag(1:1).eq.'-'.and.
     1                     (ichar(tag(2:2)).lt.48.or.
     2                               ichar(tag(2:2)).gt.57)) then
                      i = i + j - 1
                      goto 201
                  else
                      if (j.eq.1) then
                          arg = tag
                      else
                          arg = arg(1:20*(j-1))//' '//tag
                      endif
                      goto 231
                  endif
              endif
          elseif (tag(1:2).eq.'-f') then
              i = i + 1
              call getarg(i,ifile)
          elseif (tag(1:2).eq.'-o') then
              i = i + 1
              call getarg(i,ofile)
          elseif (tag(1:2).eq.'-p') then
              p = 1
          elseif (tag(1:2).eq.'-h') then
              call usage(' ')
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
     1 'Usage: dateutil -f IFILE -o OFILE -nday|-date ',
     2               '[-c DATE1 DATE2|NDAY] [-long] [-p] [-h]'
      write(*,*)
      write(*,*)
     1 '-f IFILE             Input date file (format depends on -nday|',
     1                                    '-date and -long)'
      write(*,*)
     1 '-o OFILE             Output file (depends on -nday|-date and ',
     2                                            '-long)'
      write(*,*)
     1 '-nday                Compute number of days (IFILE format: ',
     2                                  'DATE1 DATE2)'
      write(*,*)
     1 '-date                Compute calendar date (IFILE format: ',
     2                                      'DATE NDAY)'
      write(*,*)
     1 '-c DATE1 DATE2|NDAY  Input dates and days on command line ',
     2                        '(overrides -f and -o)'
      write(*,*)
     1 '-long                Read/write all dates as ',
     2                      'YYYYMMDDHHMMSS instead of YYYYMMDD'
      write (*,*)
     1 '-p                   Print results to standard output ',
     2                                  '(overrides -o)'
      write (*,*)
     1 '-h                   Online help (this message)'
      write(*,*)
      write(*,*)
     1 '  dateutil recognizes three DATE formats:'
      write(*,*)
     1 '    YYYY MM HH [HH MM SS]'
      write(*,*)
     1 '    YYYYMMHH[HHMMSS]'
      write(*,*)
     1 '    YYYY-MM-HHT[HH:MM:SS]'
      write(*,*)
     1 '  NDAY format: real number (to allow fractions of a ',
     2                                        'day)'
      write(*,*)
      STOP
      END
