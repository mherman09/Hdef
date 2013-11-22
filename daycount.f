      PROGRAM daycount
 
      IMPLICIT none
      INTEGER yr1,mo1,dy1,yr2,mo2,dy2,ndy1,ndy2
      CHARACTER*30 ifile,ofile,date1,date2

      call gcmdln(ifile,ofile)
 
      if (ifile.eq.'none') then
          print *,'To see options, use -h flag'
          print *,'Enter year1 month1 day1 year2 month2 day2:'
          read *, yr1,mo1,dy1,yr2,mo2,dy2
          call ndays(ndy1,yr1,mo1,dy1)
          call ndays(ndy2,yr2,mo2,dy2)
          print *, ndy2 - ndy1
      else      
          open (unit=11,file=ifile,status='old')
          open (unit=12,file=ofile,status='unknown')
  101     read (11,*,end=102) yr1,mo1,dy1,yr2,mo2,dy2
              call ndays(ndy1,yr1,mo1,dy1)
              call ndays(ndy2,yr2,mo2,dy2)
              write (12,*) ndy2 - ndy1
              goto 101
  102     continue
      endif

      END

C======================================================================C
 
      SUBROUTINE gcmdln(ifile,ofile)
 
      IMPLICIT NONE
      CHARACTER*30 ifile,ofile,tag
      INTEGER i,narg

      ifile = 'none'
      ofile = 'daycount.out'

      narg = iargc()
      i = 0
  201 i = i + 1
      if (i.gt.narg) goto 202
          call getarg(i,tag)
          if (tag(1:2).eq.'-f') then
              i = i + 1
              call getarg(i,ifile)
          elseif (tag(1:2).eq.'-o') then
              i = i + 1
              call getarg(i,ofile)
          elseif (tag(1:2).eq.'-h'.or.tag(1:2).eq.'-?') then
              call usage()
          endif
          goto 201
  202 continue

      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE ndays(ndy,yr,mo,dy)
      
      IMPLICIT NONE
      INTEGER ndy,yr,mo,dy,lpyr
      
      ndy = dy

      if (mo.eq.1) then
          ndy = ndy
      elseif (mo.eq.2) then
          ndy = ndy + 31
      elseif (mo.eq.3) then
          ndy = ndy + 59
      elseif (mo.eq.4) then
          ndy = ndy + 90
      elseif (mo.eq.5) then
          ndy = ndy + 120
      elseif (mo.eq.6) then
          ndy = ndy + 151
      elseif (mo.eq.7) then
          ndy = ndy + 181
      elseif (mo.eq.8) then
          ndy = ndy + 212
      elseif (mo.eq.9) then
          ndy = ndy + 243
      elseif (mo.eq.10) then
          ndy = ndy + 273
      elseif (mo.eq.11) then
          ndy = ndy + 304
      elseif (mo.eq.12) then
          ndy = ndy + 334
      else
          print *,'MONTH MUST BE IN RANGE 1-12'
      endif
      
      lpyr = yr/4
      ndy = ndy + (yr-0)*365 + lpyr
      if (mod(yr,4).eq.0.and.mo.le.2) ndy = ndy - 1

      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE usage()
      IMPLICIT none

      write(*,*)
     1 'Usage: daycount -f [IFILE] -o [OFILE] -h/-?'
      write(*,*)
     1 '  -f [IFILE] (Default none) name of input file'
      write(*,*)
     1 '                            if no file name given, fault ',
     2                            'parameters prompted for manual input'
      write(*,*)
     1 '  -o [OFILE] (Default daycount.out) name of output file'
      write (*,*)
     1 '  -h/-?        help'
      write (*,*) ''
      write (*,*)
     1 '  daycount calculates number of days between calendar dates'
      write (*,*) ''
      write (*,*)
     1 '    input file'
      write (*,*)
     1 '        yr1 mo1 dy1 yr2 mo2 dy2'
      write (*,*)
     1 '         :  :'
      write (*,*)
     1 '    output file'
      write (*,*)
     1 '        no. days (starting from date1, ending at date2)'
      write (*,*)
     1 '         :  :'
      write (*,*) ''

      STOP
      END
