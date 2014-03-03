      PROGRAM eventfrequency
C----
C From a list of integer values, calculate the frequency of each value
C in the list
C----
      IMPLICIT none
      CHARACTER*30 ifile
      INTEGER NMAX
      PARAMETER (NMAX=1000)
      INTEGER ct(NMAX)
      INTEGER mnval,mxval,dumval,mnsrch,mxsrch,val,diff
      INTEGER i,istart,iend,itmin,itmax,j,niter
      
      mnsrch = 0
      mxsrch = 0

      call gcmdln(ifile,mnval,mxval)
      if (mnval.eq.-99999) mnsrch = 1
      if (mxval.eq.-99999) mxsrch = 1

      open (unit=11,file=ifile,status='old')

      rewind 11
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

C      print *,mnval,mxval

      istart = 1
      iend = mxval - mnval + 1

      niter = iend/NMAX + 1

      do 107 j = 1,niter
C          print *,'J = ',j
          itmin = (j-1)*NMAX + 1
          itmax = j*NMAX
          do 103 i = 1,NMAX
              ct(i) = 0
  103     continue

          rewind 11
  104     read (11,*,end=105) val
              if (val.ge.itmax+mnval.or.val.le.itmin-mnval) goto 104
              val = val - itmin - mnval + 2
              ct(val) = ct(val) + 1
              goto 104
  105     continue

          do 106 i = 1,NMAX
              val = i - 2 + mnval + itmin
              if (ct(i).gt.0) print *,val,ct(i)
  106     continue
  107 continue

      END

C======================================================================C

      SUBROUTINE gcmdln(ifile,mnval,mxval)

      IMPLICIT none
      CHARACTER*30 tag
      CHARACTER*30 ifile
      INTEGER narg,i
      INTEGER mnval,mxval

      ifile = 'eventfrequency.in'
      mnval = -99999
      mxval = -99999

      narg = iargc()

      i = 0
   11 i = i + 1
      if (i.gt.narg) goto 12
          call getarg(i,tag)
          if (tag(1:2).eq.'-f') then
              i = i + 1
              call getarg(i,ifile)
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
     1 'Usage: eventfrequency -min [MNVAL] -max [MXVAL] -h/-?'
      write (*,*)
     1 '  -min [MNVAL] (Default undefined) minimum value in the list'
      write (*,*)
     1 '  -max [MXVAL] (Default undefined) maximum value in the list'
      write (*,*)
     1 '  -h/-?        help'
      write (*,*) ''
      
      STOP
      END
