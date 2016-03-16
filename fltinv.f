      PROGRAM main
C----
C Invert displacement observations for fault slip, assuming slip in an
C elastic half-space. Can run in 2-D (assumes only dip-slip) or 3-D.
C----
      IMPLICIT none
      CHARACTER*80 fltf,obsf,smoof
      REAL*8 smooth,damp
      INTEGER dimen,geo,p,dbg
      CHARACTER*80 ofile
      INTEGER OBSMAX,FLTMAX
      PARAMETER (OBSMAX=2500,FLTMAX=2000)
      REAL*8 obs(OBSMAX,6),flt(FLTMAX,7)
      REAL*8 A(OBSMAX,FLTMAX),b(OBSMAX,1),x(OBSMAX,1)
      INTEGER i,nobs,nflt
      CHARACTER*80 haff
      REAL*8 haf(3)
      INTEGER j

C----
C Parse command line and check that required values are set
C----
      call gcmdln(fltf,obsf,smoof,smooth,damp,dimen,ofile,geo,haff,p,
     1            dbg)
      if (dbg.eq.1) then
          print *,''
          print *,'****************************************'
          print *,'      FLTINV DEBUGGING OUTPUT ON'
          print *,'****************************************'
          print *,''
          print *,'COMMAND LINE INPUTS'
          print *,'-------------------'
          print *,'fltf:   ',fltf
          print *,'obsf:   ',obsf
          print *,'ofile:  ',ofile
          print *,'haff:   ',haff
          print *,'smoof:  ',smoof
          print *,'smooth: ',smooth
          print *,'damp:   ',damp
          print *,'dimen:  ',dimen
          print *,'geo:    ',geo
          print *,'p:      ',p
      endif
      call chkin(fltf,obsf,dimen)
      call chkout(ofile,p)

C----
C Read displacement and fault data
C----
      call readin(obsf,nobs,fltf,nflt,dimen,obs,flt)
      if (dbg.eq.1) then
          print *,''
          print *,'OBSERVATIONS'
          print *,'------------'
          print *,'nobs: ',nobs
          print *,'x y z ux uy uz'
          do 901 i = 1,nobs
              print *,(obs(i,j),j=1,6)
  901     continue
          print *,''
          print *,'FAULTS'
          print *,'------'
          print *,'nflt: ',nflt
          print *,'x y z s d wd ln'
          do 902 i = 1,nflt
              print *,(flt(i,j),j=1,7)
  902     continue
      endif
      call readhaf(haff,haf)

C----
C Convert geographical to Cartesian coordinates, if necessary
C----
      if (geo.eq.1.and.dimen.eq.3) then
          call geo2xy(flt,obs,nflt,nobs)
      endif

C----
C Generate Green's functions and load A matrix and b vector
C----
      if (dimen.eq.2) then
          call green2d(A,b,obs,flt,nobs,nflt,haf)
      elseif (dimen.eq.3) then
          call green3d(A,b,obs,flt,nobs,nflt,haf)
      endif
      if (dbg.eq.1) then
          print *,''
          print *,'A ='
          do 903 i = 1,nobs
              if (i.eq.1) then
                  print *,'gfx'
              elseif (i.eq.1+nobs/3) then
                  print *,'gfy'
              elseif (i.eq.1+2*nobs/3) then
                  print *,'gfz'
              endif
              print *,'[',(A(i,j),j=1,nflt/2),' |',
     1                                      (A(i,j),j=1+nflt/2,nflt),']'
  903     continue
          print *,''
          print *,'b ='
          do 904 i = 1,nobs
              if (i.eq.1) then
                  print *,'x'
              elseif (i.eq.1+nobs/3) then
                  print *,'y'
              elseif (i.eq.1+2*nobs/3) then
                  print *,'z'
              endif
              print *,'[',b(i,1),']'
  904     continue
      endif

C----
C Add damping
C----
      call damping(A,b,nobs,nflt,damp)

C----
C Add smoothing
C----
      call smoothing(A,b,nobs,nflt,smooth,smoof,dimen)

C----
C Solve generalized least squares problem
C----
      call lsqsolve(A,b,x,OBSMAX,FLTMAX,nobs,nflt)
      call chksol(x,nflt)

C----
C Write output nicely
C----
      call output(x,nflt,dimen,ofile,p)

      END

C----------------------------------------------------------------------C

      SUBROUTINE chksol(x,nflt)
      IMPLICIT none
      INTEGER OBSMAX
      PARAMETER (OBSMAX=2500)
      REAL*8 x(OBSMAX,1)
      INTEGER nflt,i
      do 601 i = 1,nflt
          if (x(i,1).gt.100.0d0) then
              write(*,*) '!! Error: found slip larger than 100 m'
              write(*,*) '!! Check inputs and units'
              write(*,*) '!! Try adding damping/smoothing'
              stop
          endif
  601 continue
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE output(x,nflt,dimen,ofile,p)
      IMPLICIT none
      INTEGER OBSMAX
      PARAMETER (OBSMAX=2500)
      REAL*8 x(OBSMAX,1)
      INTEGER nflt,dimen,p,i
      CHARACTER*80 ofile
      if (p.eq.0) then
          open(unit=23,file=ofile,status='unknown')
          if (dimen.eq.2) then
              do 103 i = 1,nflt
                  write(23,1002) x(i,1)
  103         continue
          else
              do 104 i = 1,nflt/2
                  write(23,1003) x(i,1),x(nflt/2+i,1)
  104         continue
          endif
          close(23)
      else
          if (dimen.eq.2) then
              do 105 i = 1,nflt
                  write(*,1002) x(i,1)
  105         continue
          else
              do 106 i = 1,nflt/2
                  write(*,1003) x(i,1),x(nflt/2+i,1)
  106         continue
          endif
      endif
 1002 format(F14.6)
 1003 format(2F14.6)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE chkin(fltf,obsf,dimen)
C----
C Check that input files exist and problem dimension defined
C----
      IMPLICIT none
      CHARACTER*80 fltf,obsf
      INTEGER dimen
      LOGICAL ex
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
      if (dimen.ne.2.and.dimen.ne.3) then
          write(*,*) '!! Error: Problem dimension must be 2 or 3'
          call usage('!! Use -2 for 2D, -3 for 3D')
      endif
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE chkout(ofile,p)
C----
C Check that output files defined
C----
      IMPLICIT none
      CHARACTER*80 ofile
      INTEGER p
C Output file
      if (ofile.eq.'none'.and.p.eq.0) then
          write(*,*) '!! Error: Output file not specified'
          call usage('!! Use -o OFILE to define output file or -p')
      endif
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE readin(obsf,nobs,fltf,nflt,dimen,obs,flt)
C----
C Read displacement observations from file
C----
      IMPLICIT none
      CHARACTER*80 obsf,fltf
      INTEGER nobs,nflt,dimen,i,j
      INTEGER OBSMAX,FLTMAX
      PARAMETER (OBSMAX=2500,FLTMAX=2000)
      REAL*8 obs(OBSMAX,6),flt(FLTMAX,7)
C Count number of observations
      call linect(nobs,obsf)
C Read observation locations and vector components
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
C Count number of faults in model
      call linect(nflt,fltf)
C Read fault locations and kinematics
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
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE readhaf(haff,haf)
      IMPLICIT none
      CHARACTER*80 haff
      REAL*8 haf(3),lamda,mu
      LOGICAL ex
      CHARACTER*10 ch
      if (haff.eq.'none') then
          haf(1) = 7.0d3
          haf(2) = haf(1)/dsqrt(3.0d0)
          haf(3) = 3.0d3
      else
          inquire(file=haff,EXIST=ex)
          if (.not.ex) then
              write(*,*) '!! Warning: no input file: '//haff
              write(*,*) '!! Using default half-space parameters'
              haf(1) = 7.0d3
              haf(2) = haf(1)/dsqrt(3.0d0)
              haf(3) = 3.0d3
          else
              open(unit=31,file=haff,status='old')
              read(31,*) ch
              rewind(31)
              if (ch.eq.'Lame'.or.ch.eq.'lame') then
                  read(31,*) ch,lamda,mu
                  haf(3) = 3.0d3
                  haf(1) = dsqrt(mu/haf(3))
                  haf(2) = dsqrt((lamda+2.0d0*mu)/haf(3))
              else
                  read(31,*) haf(1),haf(2),haf(3)
              endif
              close(31)
          endif
      endif
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE geo2xy(flt,obs,nflt,nobs)
C----
C Convert longitude and latitude to x and y
C----
      IMPLICIT none
      INTEGER OBSMAX,FLTMAX
      PARAMETER (OBSMAX=2500,FLTMAX=2000)
      REAL*8 obs(OBSMAX,6),flt(FLTMAX,7)
      INTEGER i,nobs,nflt
      REAL*8 dist,az,lon,lat,lon0,lat0
      ! Set origin to first fault location
      lon0 = flt(1,1)
      lat0 = flt(1,2)
      do 121 i = 1,nflt
          lon = flt(i,1)
          lat = flt(i,2)
          call ddistaz(dist,az,lon0,lat0,lon,lat)
          dist = dist*6.371d6
          flt(i,1) = dist*dsin(az)
          flt(i,2) = dist*dcos(az)
  121 continue
      do 122 i = 1,nobs
          lon = obs(i,1)
          lat = obs(i,2)
          call ddistaz(dist,az,lon0,lat0,lon,lat)
          dist = dist*6.371d6
          obs(i,1) = dist*dsin(az)
          obs(i,2) = dist*dcos(az)
  122 continue
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE green2d(A,b,obs,flt,nobs,nflt,haf)
C----
C Create Green's functions for 2-D case and load model matrix. Assume
C all faults and observations are centered on x=0.
C----
      IMPLICIT none
      INTEGER OBSMAX,FLTMAX
      PARAMETER (OBSMAX=2500,FLTMAX=2000)
      REAL*8 obs(OBSMAX,6),flt(FLTMAX,7)
      REAL*8 A(OBSMAX,FLTMAX),b(OBSMAX,1)
      INTEGER i,j,nobs,nflt
      REAL*8 haf(3),vp,vs,dens
      REAL*8 xflt,yflt,zflt,xobs,yobs,zobs,dx,dy
      REAL*8 wid,len,dip,rak,slip
      REAL*8 gfx,gfy,gfz
      CHARACTER*1 lr
C Default faults dip to left
      lr = 'l'
C Half-space properties
      vp = haf(1)
      vs = haf(2)
      dens = haf(3)
C Compute Green's functions and load model matrix, A
      xflt = 0.0d0
      xobs = 0.0d0
      rak = 90.0d0 ! only considering dip-slip
      slip = 1.0d0 ! unit slip for Green's functions
      do 151 j = 1,nflt
          yflt = flt(j,1)
          zflt = flt(j,2)
          dip = flt(j,3)
          if (dip.lt.0.0d0) then
              lr = 'r'    ! negative dip => right dipping
              dip = -dip
          endif
          wid = flt(j,4)
          len = flt(j,5)
          do 152 i = 1,nobs
              yobs = obs(i,1)
              zobs = obs(i,2)
              dx = xobs - xflt
              dy = yobs - yflt
              if (lr.eq.'r') dy = -dy
              call o92rect(gfx,gfy,gfz,dx,dy,zobs,zflt,dip,rak,wid,len,
     1                     slip,vp,vs,dens)
              A(i,j)      = gfy
              A(i+nobs,j) = gfz
  152     continue
  151 continue
C Load observation vector, b
      do 153 i = 1,nobs
          b(i,1)      = obs(i,3)
          b(i+nobs,1) = obs(i,4)
  153 continue
C In 2-D, each observation corresponds to two rows in A, b
      nobs = 2*nobs
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE green3d(A,b,obs,flt,nobs,nflt,haf)
C----
C Create Green's functions for 3-D case and load model matrix.
C----
      IMPLICIT none
      REAL*8 pi,d2r
      PARAMETER (pi=4.0d0*atan(1.0d0),d2r=pi/1.8d2)
      INTEGER OBSMAX,FLTMAX
      PARAMETER (OBSMAX=2500,FLTMAX=2000)
      REAL*8 obs(OBSMAX,6),flt(FLTMAX,7)
      REAL*8 A(OBSMAX,FLTMAX),b(OBSMAX,1)
      INTEGER i,j,jj,nobs,nflt
      REAL*8 haf(3),vp,vs,dens
      REAL*8 xflt,yflt,zflt,xobs,yobs,zobs,dx,dy
      REAL*8 wid,len,dip,rak,slip,str
      REAL*8 gfx,gfy,gfz
      REAL*8 ddip,dstr,t,cost,sint
C Half-space properties
      vp = haf(1)
      vs = haf(2)
      dens = haf(3)
      print *,haf(1),haf(2),haf(3)
C Compute Green's functions and load model matrix, A
      slip = 1.0d0 ! unit slip for Green's functions
      do 151 j = 1,nflt
          xflt = flt(j,1)
          yflt = flt(j,2)
          zflt = flt(j,3)
          str = flt(j,4)
          dip = flt(j,5)
          wid = flt(j,6)
          len = flt(j,7)
C         variables for rotating to fault reference frame
          t = d2r*(90.0d0-str) ! t for "theta"
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
                      rak = 90.0d0 ! dip-slip
                      call o92rect(gfx,gfy,gfz,dstr,ddip,zobs,zflt,dip,
     1                             rak,wid,len,slip,vp,vs,dens)
                      A(i,j)        = gfx*cost - gfy*sint
                      A(i+nobs,j)   = gfx*sint + gfy*cost
                      A(i+2*nobs,j) = gfz
                  else
                      rak = 0.0d0 ! strike-slip
                      call o92rect(gfx,gfy,gfz,dstr,ddip,zobs,zflt,dip,
     1                             rak,wid,len,slip,vp,vs,dens)
                      A(i,j+nflt)        = gfx*cost - gfy*sint
                      A(i+nobs,j+nflt)   = gfx*sint + gfy*cost
                      A(i+2*nobs,j+nflt) = gfz
                  endif
  153         continue
  152     continue
  151 continue
C Load observation vector, b
      do 154 i = 1,nobs
          b(i,1)        = obs(i,4)
          b(i+nobs,1)   = obs(i,5)
          b(i+2*nobs,1) = obs(i,6)
  154 continue
C In 3-D, each observation corresponds to three rows in A, b
      nobs = 3*nobs
      nflt = 2*nflt
      RETURN
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
      INTEGER OBSMAX,FLTMAX
      PARAMETER (OBSMAX=2500,FLTMAX=2000)
      REAL*8 A(OBSMAX,FLTMAX),b(OBSMAX,1)
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

      SUBROUTINE smoothing(A,b,nobs,nflt,smooth,smoof,dimen)
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
      INTEGER OBSMAX,FLTMAX
      PARAMETER (OBSMAX=2500,FLTMAX=2000)
      REAL*8 A(OBSMAX,FLTMAX),b(OBSMAX,1)
      INTEGER nsmoo,el,nnei,nei(FLTMAX)
      INTEGER i,j,nobs,nflt,dimen
      if (smooth.le.0.0d0) return
C Check for input smoothing file
      if (smoof.eq.'none') then
          write(*,*) '!! Error: No smoothing file specified'
          call usage('!! Use -smooth CONST FILE to specify input file')
      else
          inquire(file=smoof,EXIST=ex)
          if (.not.ex) call usage('!! Error: no input file: '//smoof)
      endif
      call linect(nsmoo,smoof)
      if (nsmoo.gt.nflt) then
          write(*,*) '!! Error: more faults to smooth than nflt'
          call usage('!! Check smoothing linking file')
      endif
C Initialize smoothing section of A and b arrays as zeros
      do 102 i = 1,nsmoo
          b(nobs+i,1) = 0.0d0
          do 103 j = 1,nflt
              A(nobs+i,j) = 0.0d0
  103     continue
  102 continue
C Read smoothing element linking file and load into A
      open(unit=31,file=smoof,status='old')
      do 101 i = 1,nsmoo
          read(31,*) el,nnei,(nei(j),j=1,nnei)
          ! Fault of interest
          A(nobs+i,el) = dble(nnei)*smooth*smooth
          if (dimen.eq.3) then
              A(nobs+i,nflt/2+el) = dble(nnei)*smooth*smooth ! nflt doubled at end of green3d
          endif
          do 104 j = 1,nnei
              A(nobs+i,nei(j)) = -1.0d0*smooth*smooth
              if (dimen.eq.3) then
                  A(nobs+i,nflt/2+nei(j)) = -1.0d0*smooth*smooth
              endif
  104     continue
  101 continue
      close(31)
      nobs = nobs + nflt
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE lsqsolve(Ain,bin,xout,OBSMAX,FLTMAX,nobs,npar)
C----
C Solve generalized least squares problem Ax = b.
C On input, arrays Ain and bin are given maximum dimensions:
C     Ain(OBSMAX,FLTMAX)
C     bin(OBSMAX,1)
C To operate with LAPACK generalized inverse tools, convert arrays from
C maximum size to operating size:
C     Ain --> A(nobs,npar)
C     bin --> b(nobs,1)
C----
      IMPLICIT none
      INTEGER OBSMAX,FLTMAX
      REAL*8 Ain(OBSMAX,FLTMAX),bin(OBSMAX,1),xout(OBSMAX,1)
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
      REAL*8 work(100000)
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
      !print *,'LWORK',lwork
C Copy observation vector, b to btmp (replaced in dgels)
      do 101 i = 1,nobs
          btmp(i,1) = b(i,1)
  101 continue
C Compute parameter array
      call dgels(trans,m,n,nrhs,a,lda,btmp,ldb,work,lwork,info)
      do 102 i = 1,npar
          x(i,1) = btmp(i,1)
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

      SUBROUTINE gcmdln(fltf,obsf,smoof,smooth,damp,dimen,ofile,geo,
     1                  haff,p,dbg)
      IMPLICIT NONE
      CHARACTER*80 fltf,obsf,tag,smoof,ofile,haff
      INTEGER i,narg
      REAL*8 smooth,damp
      INTEGER dimen,geo,p,dbg
      fltf = 'none'
      obsf = 'none'
      smoof = 'none'
      ofile = 'none'
      haff  = 'none'
      dimen = 0
      smooth = -1.0d0
      damp = -1.0d0
      geo = 0
      p = 0
      dbg = 0
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
          elseif (tag(1:6).eq.'-debug') then
              dbg = 1
          elseif (tag(1:5).eq.'-damp') then
              i = i + 1
              call getarg(i,tag)
              read(tag,'(BN,F10.0)') damp
          elseif (tag(1:2).eq.'-2') then
              dimen = 2
          elseif (tag(1:2).eq.'-3') then
              dimen = 3
          elseif (tag(1:2).eq.'-o') then
              i = i + 1
              call getarg(i,ofile)
          elseif (tag(1:4).eq.'-geo') then
              geo = 1
          elseif (tag(1:4).eq.'-haf') then
              i = i + 1
              call getarg(i,haff)
          elseif (tag(1:2).eq.'-h') then
              call usage(' ')
          elseif (tag(1:2).eq.'-d') then
              call usage('D')
          elseif (tag(1:2).eq.'-p') then
              p = 1
          else
              call usage('!! Error: No option '//tag)
          endif
          goto 901
  902 continue
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE usage(str)
      IMPLICIT none
      INTEGER lstr,d
      CHARACTER str*(*)
      if (str.eq.'D') then
          d = 1
      else
          d = 0
      endif
      if (str.ne.' '.and.d.eq.0) then
          lstr = len(str)
          write(*,*) str(1:lstr)
          write(*,*)
      endif
      write(*,*)
     1 'Usage: fltinv -obs OBSFILE -flt FLTFILE -o OFILE -2|-3 '
      write(*,*)
     1 '              [-p] [-geo] [-smooth SMOOTH FILE] ',
     2                '[-damp DAMP]'
      write(*,*)
     1 '              [-haf HAFSPC] [-h] [-d]'
      write(*,*)
      write(*,*)
     1 '-obs OBSFILE        Displacement observations'
      if (d.eq.1) then
          write(*,*)
     1 '                       2-D: y(m) z(m) uy(m) uz(m)'
          write(*,*)
     1 '                       3-D: x(m) y(m) z(m) ux(m) uy(m) uz(m)'
          write(*,*)
      endif
      write(*,*)
     1 '-flt FLTFILE        Input fault locations and geometries'
      if (d.eq.1) then
          write(*,*)
     1 '                       2-D: y(m) z(m) dip(deg) wid(m) xlen(m)'
          write(*,*)
     1 '                       3-D: x(m) y(m) z(m) str(deg) dip(deg) ',
     2                         'wid(m) len(m)'
          write(*,*)
      endif
      write(*,*)
     1 '-o OFILE            Output slip values'
      if (d.eq.1) then
          write(*,*)
     1 '                       2-D: dip_slip(m)'
          write(*,*)
     1 '                       3-D: dip_slip(m) str_slip(m)'
          write(*,*)
      endif
      write(*,*)
     1 '-2|-3               Choose 2-D or 3-D problem (3-D in ',
     2                                'development)'
      if (d.eq.1) then
          write(*,*)
      endif
      write(*,*)
     1 '-p                  Print output to standard out (overrides -o)'
      if (d.eq.1) then
          write(*,*)
      endif
      write(*,*)
     1 '-geo                Use geographic coordinates in 3-D problem'
      if (d.eq.1) then
          write(*,*)
     1 '                       X=E, Y=N, Z pos down, uz pos up'
          write(*,*)
      endif
      write(*,*)
     1 '-haf HAFSPC         Half-space parameters'
      if (d.eq.1) then
          write(*,*)
     1 '                       vp(m/s) vs(m/s) dens(kg/m^3)'
          write(*,*)
     1 '                       Lame lamda mu'
          write(*,*)
      endif
      write(*,*)
     1 '-damp CONST         Minimize length of model, with weight CONST'
      write(*,*)
     1 '-smooth CONST FILE  Smooth model, with weight CONST, linking ',
     2                         'file FILE'
      if (d.eq.1) then
          write(*,*)
     1 '                      Linking file format: flt #_neighbors ',
     2                                  'neighbor_1 neighbor_2...'
          write(*,*)
     1 '                        1    2     2  4     '
          write(*,*)
     1 '                        2    3     1  3  5  '
          write(*,*)
     1 '                        3    2     2  6     '
          write(*,*)
     1 '                        :    :      :       '
          write(*,*)
      endif
      write(*,*)
     1 '-h                  Online help'
      write(*,*)
     1 '-d                  Detailed online help'
      write(*,*)
      STOP
      END
