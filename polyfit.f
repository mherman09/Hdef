      PROGRAM polyfit
      IMPLICIT none
      CHARACTER*40 ifile,ofile,cfile,pfile,rfile,MXORDC
      LOGICAL ex
      INTEGER i,nobs,MXOBS
      PARAMETER (MXOBS=100000)
      REAL*8 xobs(MXOBS),yobs(MXOBS)
      INTEGER MXORD,ord,MXCON,ncon
      PARAMETER (MXORD=20,MXCON=21)
      REAL*8 coeff(MXORD),xc(MXCON),yc(MXCON),ypre,res
      INTEGER j,p
C----
C TO-DO:
C Subroutine chkopt(ifile,ofile,cfile,ord,MXORD)
C   check existence of input files
C   check output file defined
C   check ord <= MXORD
C Subroutine readxy(ifile,xobs,yobs,nobs)
C   read in values and number of observations from ifile
C Subroutine readconstr(cfile,xc,yc,ncon)
C   read in constraints if necessary
C----

C----
C Parse command line
C----
      call gcmdln(ifile,ofile,cfile,xc,yc,ncon,pfile,rfile,ord,p)

C----
C Check for required input files and values
C----
      if (ifile.eq.'none') then
          write(*,*) '!! Error: Input x-y data is not specified'
          call usage('!! Use -f XYFILE to specify input file')
      else
          inquire(file=ifile,EXIST=ex)
          if (.not.ex) call usage('!! Error: no x-y file: '//ifile)
      endif
      if (ord.lt.0) then
          write(*,*) '!! Error: polynomial order not specified'
          call usage('!! Use -ord ORDER to specify order')
      endif
      if (ofile.eq.'none') then
          write(*,*) '!! Error: No output file specified'
          call usage('!! Use -o COEFF to specify coefficient file')
      endif
C Check for optional input files and values
      if (cfile.ne.'none') then
          inquire(file=cfile,EXIST=ex)
          if (.not.ex) call usage('!! Error: no constraint file: '
     1                                           //cfile)
      endif
      if (ord.gt.MXORD) then
          write(MXORDC,*) MXORD
          call usage('!! Error: ORD must be <='//MXORDC)
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
      close(11)
      nobs = i - 1

C----
C Check if a constraint file was defined, and compute coefficients
C----
C     Unconstrained least squares fit
      if (cfile.eq.'none'.and.xc(1).lt.-1.0d98) then
          call polylsq(coeff,xobs,yobs,nobs,ord)
C     Constrained least squares fit
      elseif (cfile.eq.'none') then
          call polylsqconstr(coeff,xobs,yobs,nobs,ord,xc,yc,ncon)
      else
          open(unit=13,file=cfile,status='old')
          i = 1
 103      read(13,*,end=104) xc(i),yc(i)
              i = i + 1
              goto 103
 104      continue
          close(13)
          ncon = i - 1
          call polylsqconstr(coeff,xobs,yobs,nobs,ord,xc,yc,ncon)
      endif

C----
C Print coefficients
C----
      if (p.eq.0) then
          open(unit=14,file=ofile,status='unknown')
          write(14,9001)
          do 105 i = 1,ord+1
              write(14,9002) i-1,coeff(i)
  105     continue
          close(14)
      else
          write(*,9001)
          do 106 i = 1,ord+1
              write(*,9002) i-1,coeff(i)
  106     continue
      endif

C----
C Print predicted results
C----
      if (pfile.ne.'none') then
          open(unit=12,file=pfile,status='unknown')
          do 108 i = 1,nobs
              ypre = 0.0d0
              do 107 j = 1,ord+1
                  ypre = ypre + coeff(j)*xobs(i)**(j-1)
  107         continue
              write (12,9003) xobs(i),ypre
  108     continue
          close(12)
      endif

      if (rfile.ne.'none') then
          open(unit=15,file=rfile,status='unknown')
          do 110 i = 1,nobs
              ypre = 0.0d0
              do 109 j = 1,ord+1
                  ypre = ypre + coeff(j)*xobs(i)**(j-1)
  109         continue
              res = yobs(i) - ypre
              write(15,9003) xobs(i),res
  110     continue
      endif

 9001 format('POLY_ORDER      COEFFICIENT')
 9002 format(I10,X,E16.6)
 9003 format(2E14.6)

      END

C======================================================================C

      SUBROUTINE gcmdln(ifile,ofile,cfile,xc,yc,ncon,pfile,rfile,ord,p)
      IMPLICIT NONE
      INTEGER MXCON
      PARAMETER (MXCON=21)
      REAL*8 xc(MXCON),yc(MXCON)
      CHARACTER*40 ifile,ofile,cfile,pfile,rfile,tag
      INTEGER i,j,narg,ord,p,ncon
      j = 0
      xc(1) = -1.0d99
      ifile = 'none'
      ofile = 'none'
      cfile = 'none'
      pfile = 'none'
      rfile = 'none'
      p = 0
      ord = -1
      ncon = 0
      narg = iargc()
      if (narg.eq.0) call usage('!! Error: no command line arguments '//
     1                                         'specified')
      i = 0
  101 i = i + 1
      if (i.gt.narg) goto 102
          call getarg(i,tag)
          if (tag(1:4).eq.'-pre') then
              i = i + 1
              call getarg(i,pfile)
          elseif (tag(1:4).eq.'-ord') then
              i = i + 1
              call getarg(i,tag)
              read(tag,'(BN,I5)') ord
          elseif (tag(1:4).eq.'-res') then
              i = i + 1
              call getarg(i,rfile)
          elseif (tag(1:3).eq.'-cf') then
              i = i + 1
              call getarg(i,cfile)
          elseif (tag(1:2).eq.'-c') then
              j = j + 1
              i = i + 1
              call getarg(i,tag)
              read(tag,'(BN,F10.0)') xc(j)
              i = i + 1
              call getarg(i,tag)
              read(tag,'(BN,F10.0)') yc(j)
          elseif (tag(1:2).eq.'-p') then
              p = 1
          elseif (tag(1:2).eq.'-o') then
              i = i + 1
              call getarg(i,ofile)
          elseif (tag(1:2).eq.'-f') then
              i = i + 1
              call getarg(i,ifile)
          elseif (tag(1:2).eq.'-h') then
              call usage(' ')
          else
              call usage('!! Error: No option '//tag)
          endif
          goto 101
  102 continue
      ncon = j
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
     1 'Usage: polyfit -f XYINPUT -o COEFFS -ord ORDER [-p] ',
     2              '[-c XC YC | -cf CONSTR]'
      write(*,*)
     1 '               [-pre PREDICTED] [-res RESFILE] [-h]'
      write(*,*)
      write(*,*)
     1 '-f XYINPUT         Input X-Y file'
      write(*,*)
     1 '-o COEFFS          File with computed polynomial coefficients'
      write(*,*)
     1 '-ord ORDER         Order of polynomial (integer)'
      write(*,*)
     1 '-p                 Print polynomial coefficients (overrides -o)'
      write(*,*)
     1 '-c XC YC           Define constraint coordinates on command ',
     2                            'line (can use multiple times)'
      write(*,*)
     1 '-cf CONSTR         File with X-Y equality constraints ',
     2                                '(overrides -c)'
      write(*,*)
     1 '-pre PREFILE       Values of polynomial at input ',
     2                                    'X-coordinates'
      write(*,*)
     1 '-res RESFILE       Residual values (observed - predicted)'
      write(*,*)
     1 '-h                 Online help'
      write(*,*)
      STOP
      END
