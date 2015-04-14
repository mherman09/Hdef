      PROGRAM MAIN
      IMPLICIT none
      CHARACTER*80 fltf,obsf,smoof
      INTEGER dimen,ios
      LOGICAL ex
      REAL*8 smooth,damp
      INTEGER OBSMAX,PARMAX
      PARAMETER (OBSMAX=10000,PARMAX=100)
      REAL*8 flt(PARMAX,7),obs(OBSMAX,6)
      REAL*8 A(OBSMAX,PARMAX),b(OBSMAX,1),x(OBSMAX,1)
      INTEGER i,j,nobs,nflt

C----
C Parse command line and check required inputs present
C----
      call gcmdln(fltf,obsf,smooth,smoof,damp,dimen)
C Fault geometry file
      if (fltf.eq.'none') then
          write(*,*) '!! Error: Input fault geometry file unspecified'
          call usage('!! Use -flt FLTFILE to specify input file')
      else
          inquire(file=fltf,EXIST=ex)
          if (.not.ex) call usage('!! Error: no input file: '//fltf)
      endif
C Displacement observation file
      if (obsf.eq.'none') then
          write(*,*) '!! Error: Observation file unspecified'
          call usage('!! Use -obs OBSFILE to specify input file')
      else
          inquire(file=obsf,EXIST=ex)
          if (.not.ex) call usage('!! Error: no input file: '//obsf)
      endif
C Problem dimension
      if (dimen.eq.0) then
          write(*,*) '!! Error: Problem dimension unspecified'
          call usage('!! Use -2 for 2D, -3 for 3D')
      endif

C----
C Count number of observations, and read in data
C----
      call linect(nobs,obsf)
      open(unit=21,file=obsf,status='old')
      do 101 i = 1,nobs
          if (dimen.eq.2) then
              ! 2D: y(m) z(m) uy(m) uz(m)
              read(21,*) (obs(i,j),j=1,4)
          elseif (dimen.eq.3) then
              ! 3D: x(m) y(m) z(m) ux(m) uy(m) uz(m)
              read(21,*) (obs(i,j),j=1,6)
          endif
  101 continue
      close(21)

C----
C Count number of fault segments, and read in data
C----
      call linect(nflt,fltf)
      open(unit=22,file=fltf,status='old')
      do 102 i = 1,nflt
          if (dimen.eq.2) then
              ! 2D: y(m) z(m) dip(deg) wid(m) len(m)
              read(22,*) (flt(i,j),j=1,5)
          elseif (dimen.eq.3) then
              ! 3D: x(m) y(m) z(m) str(deg) dip(deg) wid(m) len(m)
              read(22,*) (flt(i,j),j=1,7)
          endif
  102 continue
      close(22)

C----
C Generate Green's functions and load A matrix and b vector
C----
      if (dimen.eq.2) then
          call green2d(A,b,obs,flt,nobs,nflt)
      elseif (dimen.eq.3) then
          call green3d(A,b,obs,flt,nobs,nflt)
      endif

C----
C Add damping
C----
      call damping(A,b,nobs,nflt,damp)

C----
C Add smoothing
C----
      call smoothing(A,b,nobs,nflt,smooth,smoof)
      
C----
C Solve generalized least squares problem
C----
      call lsq(A,b,x,OBSMAX,PARMAX,nobs,nflt)
      open(unit=23,file='fltinv.out',status='unknown')
      do 103 i = 1,nflt
          write(23,*) x(i,1)
  103 continue
      END

C----------------------------------------------------------------------C

      SUBROUTINE damping(A,b,nobs,nflt,damp)
C----
C Minimize length of solution to Ax = b with weighting factor damp.
C m=nobs n=nflt d=damp
C         [ A11 A12 ... A1n ]
C     A = [ A21 A22 ... A2n ]
C         [  :           :  ]
C         [ Am1 Am2 ... Amn ]
C
C         [ d*d  0  ...  0  ]
C  Append [  0  d*d ...  0  ]
C         [                 ]
C         [  0   0  ... d*d ]
C----
      IMPLICIT none
      REAL*8 damp
      INTEGER OBSMAX,PARMAX
      PARAMETER (OBSMAX=10000,PARMAX=100)
      REAL*8 A(OBSMAX,PARMAX),b(OBSMAX,1)
      INTEGER i,j,nobs,nflt
      if (damp.le.0.0d0) return
      do 101 i = 1,nflt
          b(nobs+i,1) = 0.0d0
          do 102 j = 1,nflt
              if (i.eq.j) then
                  A(nobs+i,j) = damp*damp*1.0d0
              else
                  A(nobs+i,j) = 0.0d0
              endif
  102     continue
  101 continue
      nobs = nobs + nflt
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE smoothing(A,b,nobs,nflt,smooth,smoof)
C----
C Minimize roughness of solution to Ax = b with weighting factor smooth.
C m=nobs n=nflt d=damp
C         [ A11 A12 ... A1n ]
C     A = [ A21 A22 ... A2n ]
C         [  :           :  ]
C         [ Am1 Am2 ... Amn ]
C
C  (e.g.) [ d*d -d*d    0     0  ...   0   0  ]
C  Append [-d*d 2*d*d -d*d    0  ...   0   0  ]
C         [  0  -d*d  2*d*d -d*d ...   0   0  ]
C         [  0    0     0     0  ... -d*d d*d ]
C----
      IMPLICIT none
      CHARACTER*80 smoof
      LOGICAL ex
      REAL*8 smooth
      INTEGER OBSMAX,PARMAX
      PARAMETER (OBSMAX=10000,PARMAX=100)
      REAL*8 A(OBSMAX,PARMAX),b(OBSMAX,1)
      INTEGER nsmoo,el,nnei,nei(PARMAX)
      INTEGER i,j,nobs,nflt
      if (smooth.le.0.0d0) return 
C Check for input smoothing file
      if (smoof.eq.'none') then
          write(*,*) '!! Error: No smoothing file specified'
          call usage('!! Use -smooth CONST FILE to specify input file')
      else
          inquire(file=smoof,EXIST=ex)
          if (.not.ex) call usage('!! Error: no input file: '//smoof)
      endif
C Initialize smoothing section of A and b arrays as zeros
      do 102 i = 1,nflt
          b(nobs+i,1) = 0.0d0
          do 103 j = 1,nflt
              A(nobs+i,j) = 0.0d0
  103     continue
  102 continue
C Read smoothing element linking file and load into A
      call linect(nsmoo,smoof)
      open(unit=31,file=smoof,status='old')
      do 101 i = 1,nsmoo
          read(31,*) el,nnei,(nei(j),j=1,nnei)
          A(nobs+el,el) = dble(nnei)*smooth*smooth
          do 104 j = 1,nnei
              A(nobs+el,nei(j)) = -1.0d0*smooth*smooth
  104     continue
  101 continue
      close(31)
      nobs = nobs + nflt
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE green2d(A,b,obs,flt,nobs,nflt)
C----
C Create Green's functions for 2-D case and load model matrix
C----
      IMPLICIT none
      INTEGER OBSMAX,PARMAX
      PARAMETER (OBSMAX=10000,PARMAX=100)
      REAL*8 flt(PARMAX,7),obs(OBSMAX,6),A(OBSMAX,PARMAX),b(OBSMAX,1)
      INTEGER i,j,nobs,nflt
      REAL*8 vp,vs,dens
      REAL*8 xflt,yflt,zflt,xobs,yobs,zobs,dx,dy
      REAL*8 wid,len,dip,rak,slip
      REAL*8 gfx,gfy,gfz
C Poisson ratio = 0.25
      vp = 7.0d3
      vs = vp/dsqrt(3.0d0)
      dens = 3.0d3
C In 2-D, assume all faults, observations centered on x=0
C Load model matrix, A
      xflt = 0.0d0
      rak = 90.0d0
      xobs = 0.0d0
      slip = 1.0d0
      do 151 j = 1,nflt
          yflt = flt(j,1)
          zflt = flt(j,2)
          dip = flt(j,3)
          wid = flt(j,4)
          len = flt(j,5)
          do 152 i = 1,nobs
              yobs = obs(i,1)
              zobs = obs(i,2)
              dx = xobs - xflt
              dy = yobs - yflt
              call o92rect(gfx,gfy,gfz,dx,dy,zobs,zflt,dip,rak,wid,len,
     1                     slip,vp,vs,dens)
              A(i,j) = gfy
              A(i+nobs,j) = gfz
  152     continue
  151 continue
C Load observation vector, b
      do 153 i = 1,nobs
          b(i,1) = obs(i,3)
          b(i+nobs,1) = obs(i,4)
  153 continue
C In 2-D, each observation corresponds to two rows in A, b
      nobs = 2*nobs
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE green3d(A,b,obs,flt,nobs,nflt)
      IMPLICIT none
      REAL*8 pi,d2r
      PARAMETER (pi=4.0d0*atan(1.0d0),d2r=pi/1.8d2)
      INTEGER OBSMAX,PARMAX
      PARAMETER (OBSMAX=10000,PARMAX=100)
      REAL*8 flt(PARMAX,7),obs(OBSMAX,6),A(OBSMAX,PARMAX),b(OBSMAX,1)
      INTEGER i,j,jj,nobs,nflt
      REAL*8 vp,vs,dens
      REAL*8 xflt,yflt,zflt,xobs,yobs,zobs,dx,dy
      REAL*8 wid,len,dip,rak,slip,str
      REAL*8 gfx,gfy,gfz
      REAL*8 ddip,dstr,t,cost,sint
C Poisson ratio = 0.25
      vp = 7.0d3
      vs = vp/dsqrt(3.0d0)
      dens = 3.0d3
C Load model matrix, A
      slip = 1.0d0
      do 151 j = 1,nflt
          xflt = flt(j,1)
          yflt = flt(j,2)
          zflt = flt(j,3)
          str = flt(j,4)
          dip = flt(j,5)
          wid = flt(j,6)
          len = flt(j,7)
C         variables for rotating to fault reference frame
          t = d2r*(90.0d0-str)
          cost = dcos(t)
          sint = dsin(t)
          do 152 i = 1,nobs
              xobs = obs(i,1)
              yobs = obs(i,2)
              zobs = obs(i,3)
              dx = xobs - xflt
              dy = yobs - yflt
C             rotate to along-strike, horizontal up-dip direction
              dstr =  dx*cost + dy*sint
              ddip = -dx*sint + dy*cost
C             need Green's functions for strike- and dip-slip components
              do 153 jj = 1,2
                  if (jj.eq.1) then
                      rak = 90.0d0
                      call o92rect(gfx,gfy,gfz,dstr,ddip,zobs,zflt,dip,
     1                             rak,wid,len,slip,vp,vs,dens)
                      A(i,j) = gfx*cost - gfy*sint
                      A(i+nobs,j) = gfx*sint + gfy*cost
                      A(i+2*nobs,j) = gfz
                  else
                      rak = 0.0d0
                      call o92rect(gfx,gfy,gfz,dstr,ddip,zobs,zflt,dip,
     1                             rak,wid,len,slip,vp,vs,dens)
                      A(i,j+nflt) = gfx*cost - gfy*sint
                      A(i+nobs,j+nflt) = gfx*sint + gfy*cost
                      A(i+2*nobs,j+nflt) = gfz
                  endif
  153         continue
  152     continue
  151 continue
C Load observation vector, b
      do 154 i = 1,nobs
          b(i,1) = obs(i,4)
          b(i+nobs,1) = obs(i,5)
          b(i+2*nobs,1) = obs(i,6)
  154 continue
C In 3-D, each observation corresponds to three rows in A, b
      nobs = 3*nobs
      nflt = 2*nflt
      RETURN
      END
C----------------------------------------------------------------------C

      SUBROUTINE lsq(Ain,bin,xout,OBSMAX,PARMAX,nobs,npar)
C----
C Solve generalized least squares problem Ax = b.
C On input, arrays Ain and bin are given maximum dimensions:
C     Ain(OBSMAX,PARMAX)
C     bin(OBSMAX,1)
C To operate with LAPACK generalized inverse tools, convert arrays from
C maximum size to operating size:
C     Ain --> A(nobs,npar)
C     bin --> b(nobs,1)
C----
      IMPLICIT none
      INTEGER OBSMAX,PARMAX
      REAL*8 Ain(OBSMAX,PARMAX),bin(OBSMAX,1),xout(OBSMAX,1)
      INTEGER nobs,npar
      REAL*8 A(nobs,npar),b(nobs,1),x(nobs,1)
      INTEGER i,j
C Load existing arrays into arrays of correct dimensions
      do 101 i = 1,nobs
          b(i,1) = bin(i,1)
          do 102 j = 1,npar
              A(i,j) = Ain(i,j)
  102     continue
  101 continue
C Solve generalized least squares problem using LAPACK tools
      call solve(A,x,b,nobs,npar)
C Print results
      do 103 i = 1,npar
          xout(i,1) = x(i,1)
  103 continue
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE solve(A,x,b,nobs,npar)
C----
C Determine the optimal workspace size for work array in dgels
C Solve least squares problem Ax = b for x using QR or LQ decomposition
C----
      IMPLICIT none
      INTEGER nobs,npar,nrhs,i
      REAL*8 A(nobs,npar),b(nobs,1),btmp(nobs,1),x(npar,1)
      INTEGER lwork
      INTEGER m,n,lda,ldb
      INTEGER info
      CHARACTER*1 trans
      REAL*8 work(1000)
      nrhs = 1
      trans = 'N'     ! A has form (nobs x npar), i.e. not transposed
      m = nobs        ! Number of rows in matrix A
      n = npar        ! Number of columns in matrix A
      lda = m         ! Leading dimension of A  lda >= max(1,m)
      ldb = max(m,n)  ! Leading dimension of b  ldb >= max(1,m,n)
C Compute optimal workspace
      lwork = -1
      call dgels(trans,m,n,nrhs,A,lda,b,ldb,work,lwork,info)
      lwork = int(work(1))
C      print *,'LWORK',lwork
C Copy observation vector, b to btmp (replaced in dgels)
      do 101 i = 1,nobs
          btmp(i,1) = b(i,1)
  101 continue
C Compute parameter array
      call dgels(trans,m,n,nrhs,a,lda,btmp,ldb,work,lwork,info)
      do 102 i = 1,npar
          x(i,1) = btmp(i,1)
C          print *,x(i,1)
  102 continue
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE linect(nline,ifile)
      IMPLICIT NONE
      CHARACTER*80 ifile
      INTEGER nline
      open(unit=41,file=ifile,status='old')
      nline = 0
  101 read(41,*,end=102)
          nline = nline + 1
          goto 101
  102 continue
      close(41)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE gcmdln(fltf,obsf,smooth,smoof,damp,dimen)
      IMPLICIT NONE
      CHARACTER*80 fltf,obsf,tag,smoof
      INTEGER i,narg
      REAL*8 smooth,damp
      INTEGER dimen
      fltf = 'none'
      obsf = 'none'
      smoof = 'none'
      dimen = 0
      smooth = -1.0d0
      damp = -1.0d0
      narg = iargc()
      if (narg.eq.0) call usage('!! Error: no command line arguments')
      i = 0
  901 i = i + 1
      if (i.gt.narg) goto 902
          call getarg(i,tag)
          if (tag(1:4).eq.'-flt') then
              i = i + 1
              call getarg(i,fltf)
          elseif (tag(1:4).eq.'-obs') then
              i = i + 1
              call getarg(i,obsf)
          elseif (tag(1:7).eq.'-smooth') then
              i = i + 1
              call getarg(i,tag)
              read(tag,'(BN,F10.0)') smooth
              i = i + 1
              call getarg(i,smoof)
          elseif (tag(1:5).eq.'-damp') then
              i = i + 1
              call getarg(i,tag)
              read(tag,'(BN,F10.0)') damp
          elseif (tag(1:2).eq.'-2') then
              dimen = 2
          elseif (tag(1:2).eq.'-3') then
              dimen = 3
          elseif (tag(1:2).eq.'-h') then
              call usage(' ')
          endif
          goto 901
  902 continue
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE usage(str)
      IMPLICIT none
      INTEGER lstr
      CHARACTER str*(*)
      if (str.ne.' ') then
          lstr = len(str)
          write(*,*) str(1:lstr)
          write(*,*)
      endif
      write(*,*)
     1 'Usage: fltinv -obs OBSFILE -flt FLTFILE -2|-3 ',
     2         '[-LQR|SVD] [-smooth SMOOTH] [-smoof SMOOTHFILE] ',
     3         '[-damp DAMP] [-h]'
      write(*,*)
      write(*,*)
     1 '-obs OBSFILE       Input observations'
      write(*,*)
     1 '                       X=E, Y=N, Z pos down, uz pos up'
      write(*,*)
     1 '                       2-D: y(m) z(m) uy(m) uz(m)'
      write(*,*)
     1 '                       3-D: x(m) y(m) z(m) ux(m) uy(m) uz(m)'
      write(*,*)
     1 '-flt FLTFILE       Input fault geometries'
      write(*,*)
     1 '                       2-D: y(m) z(m) dip(deg) wid(m) len(m)'
      write(*,*)
     1 '                       3-D: x(m) y(m) z(m) str(deg) dip(deg) ',
     2                         'wid(m) len(m)'
      write(*,*)
     1 '-2|-3              Choose 2-D or 3-D problem (3-D in ',
     2                                'development)'
      write(*,*)
     1 '-smooth CONST FILE Smooth model, with weight CONST, linking ',
     2                         'file FILE (format below)'
      write(*,*)
     1 '                       EL NNEI NEI(1,NNEI)'
      write(*,*)
     1 '                        1    2   2  4     '
      write(*,*)
     1 '                        2    3   1  3  5  '
      write(*,*)
     1 '                        3    2   2  6     '
      write(*,*)
     1 '                        4    3   1  5  7  '
      write(*,*)
     1 '                        5    4   2  4  6  8'
      write(*,*)
     1 '                        6    3   3  5  9  '
      write(*,*)
     1 '                        7    2   4  8     '
      write(*,*)
     1 '                        8    3   5  7  9  '
      write(*,*)
     1 '                        9    2   6  8     '
      write(*,*)
     1 '                        *-----*-----*-----* '
      write(*,*)
     1 '                        |  1  |  2  |  3  | '
      write(*,*)
     1 '                        |     |     |     | '
      write(*,*)
     1 '                        *-----*-----*-----* '
      write(*,*)
     1 '                        |  4  |  5  |  6  | '
      write(*,*)
     1 '                        |     |     |     | '
      write(*,*)
     1 '                        *-----*-----*-----* '
      write(*,*)
     1 '                        |  7  |  8  |  9  | '
      write(*,*)
     1 '                        |     |     |     | '
      write(*,*)
     1 '                        *-----*-----*-----* '
      write(*,*)
     1 '-damp CONST         Minimize length of model, with weight CONST'
      write(*,*)
     1 '-LQR|-SVD          Solve LS problem with QR/LQ factorization ',
     2                     'or SVD (NO SVD OPTION)'
      write(*,*)
     1 '-h                 Online help (this screen)'
      write(*,*)
      STOP
      END

C dgelsd (M,N,NRHS,A,LDA,B,LDB, S, RCOND, RANK, WORK, LWORK, IWORK, INFO)
C S(min(M,N))
C RCOND = tolerance fraction (Values of S < RCOND*S(1) are set to zero)
C RANK: integer
C IWORK(3*min(M,N)*nlvl+11*min(m,n))

