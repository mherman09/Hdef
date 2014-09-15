      PROGRAM polyfit
      IMPLICIT none
      CHARACTER*30 ifile,ofile,cfile
      LOGICAL ex
      INTEGER NMAX,nobs
      PARAMETER (NMAX=100000)
      REAL*8 xobs(NMAX),yobs(NMAX)

      INTEGER ORDMAX,ord,CONSTRMAX,nconstr
      PARAMETER (ORDMAX=10,CONSTRMAX=11)
      REAL*8 coeff(ORDMAX),xconstr(CONSTRMAX),yconstr(CONSTRMAX),
     1       xpre,ypre
      INTEGER i,j,npre

C----
C Parse command line and check for input file
C----
      call gcmdln(ifile,ofile,cfile,ord)
      inquire(file=ifile,EXIST=ex)
      if (.not.ex) then
          print *,'No input file: ',ifile
          call usage()
      endif

C----
C Extract observed data from input time series
C----
      open(unit=11,file=ifile,status='old')
      i = 1
  101 read(11,*,end=102) xobs(i),yobs(i)
          i = i + 1
          goto 101
  102 continue
      nobs = i - 1

C----
C Check if a constraint file was defined, and compute coefficients
C----
      if (cfile.eq.'none') then
          call polylsq(coeff,xobs,yobs,nobs,ord)
      else
          open(unit=13,file=cfile,status='unknown')
          i = 1
 103      read(13,*,end=104) xconstr(i),yconstr(i)
              i = i + 1
              goto 103
 104      continue
          nconstr = i - 1
          call polylsqconstr(coeff,xobs,yobs,nobs,ord,xconstr,yconstr,
     1                         nconstr)
      endif

C----
C Write polynomial coefficients to file
C----
      open(unit=14,file='coeff.tmp',status='unknown')
      do 105 i = 1,ord+1
          write (14,*) coeff(i),'ORDER=',i-1
  105 continue

C----
C Write predicted results to ofile
C----
      open(unit=12,file=ofile,status='unknown')
      do 107 i = 1,nobs
          ypre = 0.0d0
          do 106 j = 1,ord+1
              ypre = ypre + coeff(j)*xobs(i)**(j-1)
  106     continue
          write (12,*) xobs(i),ypre
  107 continue

      END

C======================================================================C

      SUBROUTINE gcmdln(ifile,ofile,cfile,ord)
      
      IMPLICIT NONE
      CHARACTER*30 ifile,ofile,cfile,tag
      INTEGER i,narg,ord
      
      ifile = 'polyfit.in'
      ofile = 'polyfit.out'
      cfile = 'none'
      ord = 1
      
      narg = iargc()
      if (narg.eq.0) call usage()
      i = 0
  101 i = i + 1
      if (i.gt.narg) goto 102
          call getarg(i,tag)
          if (tag(1:2).eq.'-h'.or.tag(1:2).eq.'-?') then
              call usage()
          elseif (tag(1:4).eq.'-ord') then
              i = i + 1
              call getarg(i,tag)
              read(tag,'(BN,I3)') ord
          elseif (tag(1:2).eq.'-f') then
              i = i + 1
              call getarg(i,ifile)
          elseif (tag(1:2).eq.'-o') then
              i = i + 1
              call getarg(i,ofile)
          elseif (tag(1:7).eq.'-constr') then
              i = i + 1
              call getarg(i,cfile)
          endif
          goto 101
  102 continue

      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE usage()
      IMPLICIT none
      write(*,*)
     1 'Usage: polyfit -f [IFILE] -o [OFILE] -ord [ORDER] -constr ',
     2                '[CFILE] -h/-?'
      write(*,*)
     1 '  -f [IFILE] (Default polyfit.in) name of input file'
      write(*,*)
     1 '  -o [OFILE] (Default polyfit.out) name of output file'
      write(*,*)
     1 '  -ord [ORDER] (Default 1) order of polynomial'
      write(*,*)
     1 '  -constr [CFILE] (Default none) equality constraint file'
      write(*,*)
     1 '  -h/-?    Online help (this screen)'
      write(*,*)
      write(*,*)
     1 'polyfit computes a least squares polynomial fit'
      write(*,*)
      write(*,*)
     1 '  Input files'
      write(*,*)
     1 '    polyfit.in: list of ordered pairs to fit'
      write(*,*)
     1 '      xobs1 yobs1'
      write(*,*)
     1 '        :     :'
      write(*,*)
      write(*,*)
     1 '    constraints.in (optional with -constr): equality ',
     2                                                'constraints'
      write(*,*)
     1 '      xconstr1 yconstr1'
      write(*,*)
     1 '         :        :'
      write(*,*)
      write(*,*)
     1 '  Output files'
      write(*,*)
     1 '    polyfit.out: list of observed x and predicted y values'
      write(*,*)
     1 '      xobs1 ypre1'
      write(*,*)
     1 '        :     :'
      write(*,*)
      write(*,*)
     1 '    coeff.tmp: polynomial coefficients (c_0*x^0 + ',
     2                                             '...c_ord*x^ord)'
      write(*,*)
     1 '      coeff0 ORDER= order'
      write(*,*)
     1 '         :           :'
      STOP
      END




