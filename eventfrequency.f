      PROGRAM eventfrequency
C----
C From a list of integer values, calculate the frequency of each value
C in the list
C----
      IMPLICIT none
      CHARACTER*30 ifile,ofile
      INTEGER NMAX
      PARAMETER (NMAX=1000)
      INTEGER ct(NMAX)
      INTEGER mnval,mxval,dumval,mnsrch,mxsrch,val,diff
      INTEGER i,istart,iend,itmin,itmax,j,niter,user,p,z
      
      mnsrch = 0
      mxsrch = 0

      call gcmdln(ifile,ofile,mnval,mxval,z,p,user)
      open (unit=12,file=ofile,status='unknown')
      if (mnval.eq.-99999) mnsrch = 1
      if (mxval.eq.-99999) mxsrch = 1

      if (user.eq.1) then
          write (*,*) 'Enter a list of integers. To end the list,',
     1                ' enter -99999'
          open (unit=11,file=ifile,status='unknown')
  201     read (*,*,end=202) val
              if (val.eq.-99999) goto 202
              write (11,*) val
              goto 201
  202     continue
      else
          open (unit=11,file=ifile,status='old')
      endif
      rewind 11

C----
C Find the minimum and maximum values in the list of integers
C----
      if (mnsrch.eq.1.or.mxsrch.eq.1) then
          read (11,*) dumval
          if (mnsrch.eq.1) mnval = dumval
          if (mxsrch.eq.1) mxval = dumval
  101     read (11,*,end=102) dumval
              if (dumval.lt.mnval.and.mnsrch.eq.1) mnval = dumval
              if (dumval.gt.mnval.and.mxsrch.eq.1) mxval = dumval
              goto 101
  102     continue
      endif
      diff = mxval - mnval

C----
C Print descriptive header
C----
      if (p.eq.0) then
          write (12,8889)
      else
          write (*,8889)
      endif
 8889 format('       VALUE   FREQUENCY')

C----
C Count the frequency of values in the list in intervals of 1000
C Bookkeeping is a bit annoying...
C----
      istart = 1
      iend = mxval - mnval + 1
      niter = iend/NMAX + 1
      do 107 j = 1,niter
          itmin = (j-1)*NMAX
          itmax = j*NMAX
          do 103 i = 1,NMAX
              ct(i) = 0
  103     continue
          rewind 11
  104     read (11,*,end=105) val
              if (val.ge.itmax+mnval.or.val.lt.mnval+itmin) goto 104
              val = val - itmin - mnval + 2
              ct(val) = ct(val) + 1
              goto 104
  105     continue
          do 106 i = 1,NMAX
              val = i - 2 + mnval + itmin
              if (z.eq.0) then
                  if (ct(i).gt.0.and.p.eq.0) then
                      write (12,9999) val,ct(i)
                  elseif (ct(i).gt.0.and.p.eq.1) then
                      write (*,9999) val,ct(i)
                  endif
              else
                  if (p.eq.0) then
                      write (12,9999) val,ct(i)
                  elseif (p.eq.1) then
                      write (*,9999) val,ct(i)
                  endif
              endif
  106     continue
  107 continue

 9999 format(2I12)

      END

C======================================================================C

      SUBROUTINE gcmdln(ifile,ofile,mnval,mxval,z,p,user)

      IMPLICIT none
      CHARACTER*30 tag
      CHARACTER*30 ifile,ofile
      INTEGER narg,i
      INTEGER mnval,mxval,user,p,z

      ifile = 'eventfrequency.in'
      ofile = 'eventfrequency.out'
      mnval = -99999
      mxval = -99999
      user = 0
      p = 0
      z = 0

      narg = iargc()
      if (narg.eq.0) call usage()

      i = 0
   11 i = i + 1
      if (i.gt.narg) goto 12
          call getarg(i,tag)
          if (tag(1:2).eq.'-f') then
              i = i + 1
              call getarg(i,ifile)
          elseif (tag(1:2).eq.'-o') then
              i = i + 1
              call getarg(i,ofile)
          elseif (tag(1:2).eq.'-d') then
              write(*,*) 'Running with default file names'
          elseif (tag(1:2).eq.'-u') then
              user = 1
              p = 1
          elseif (tag(1:2).eq.'-p') then
              p = 1
          elseif (tag(1:2).eq.'-z') then
              z = 1
          elseif (tag(1:4).eq.'-min') then
              i = i + 1
              call getarg(i,tag)
              read (tag,'(BN,I10)') mnval
          elseif (tag(1:4).eq.'-max') then
              i = i + 1
              call getarg(i,tag)
              read (tag,'(BN,I10)') mxval
          elseif (tag(1:2).eq.'-h'.or.tag(1:2).eq.'-?') then
              call usage()
          endif
      goto 11
   12 continue

      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE usage()
      IMPLICIT none

      write(*,*)
     1 'Usage: eventfrequency -f [IFILE] -o [OFILE] -d -u -p',
     2                        ' -min [MNVAL] -max [MXVAL] -z -h/-?'
      write(*,*)
     1 '  -f [IFILE] (Default eventfrequency.in) name of input file'
      write(*,*)
     1 '  -o [OFILE] (Default eventfrequency.out) name of output file'
      write(*,*)
     1 '  -d         Run with default file names'
      write (*,*)
     1 '  -u         Prompt user to enter information through standard'
      write (*,*)
     1 '  -p         Print to standard output instead of a file'
      write (*,*)
     1 '                 input for single calculation'
      write (*,*)
     1 '  -min [MNVAL] (Default undefined) minimum value in the list'
      write (*,*)
     1 '  -max [MXVAL] (Default undefined) maximum value in the list'
      write (*,*)
     1 '  -z         Print values with counts of zero'
      write (*,*)
     1 '  -h/-?        Online help (this screen)'
      write (*,*) ''
      
      STOP
      END
