      PROGRAM main
C----
C Invert displacement observations for fault slip, assuming slip in an
C elastic half-space. Based on equations of Okada (1992).
C
C Inversion options include:
C     - Least squares with LAPACK tools
C     - Simulated annealing using Gibbs sampler
C----
      IMPLICIT none
C Command line inputs
      CHARACTER*80 fltf,obsf,smoof,haff,ofile,annf
      REAL*8 smooth,damp,damp0,fact
      INTEGER geo,p,invert,vrb,ostyle,misfit,getdmp
C Inversion array variables
      INTEGER nobs,OBSMAX,nflt,FLTMAX
      PARAMETER (OBSMAX=500,FLTMAX=500)
      REAL*8 obs(OBSMAX,6),flt(FLTMAX,7),gf(OBSMAX,FLTMAX,6)
      REAL*8 soln(FLTMAX,2)
      REAL*8 pre(OBSMAX,3),tot
C Other variables
      REAL*8 haf(3)
      INTEGER i,j,k
      COMMON /VERBOSE/ vrb

C----
C Parse command line and test inputs
C----
      call gcmdln(fltf,obsf,smoof,smooth,damp,ofile,p,geo,haff,invert,
     1            annf,ostyle,misfit,getdmp,fact)
      if (vrb.eq.1) then
          write(0,*)
          write(0,'("****************************************")')
          write(0,'("      FLTINV VERBOSE OUTPUT ON          ")')
          write(0,'("****************************************")')
          write(0,*)
          write(0,'("COMMAND LINE INPUTS")')
          write(0,'("-------------------")')
          write(0,'("fltf:   ",4X,A)'),fltf
          write(0,'("obsf:   ",4x,A)'),obsf
          write(0,'("ofile:  ",4X,A)'),ofile
          write(0,'("haff:   ",4X,A)'),haff
          write(0,'("smoof:  ",4X,A)'),smoof
          write(0,'("annf:   ",4X,A)'),annf
          write(0,'("smooth: ",F12.6)'),smooth
          write(0,'("damp:   ",F12.6)'),damp
          write(0,'("geo:    ",I5)'),geo
          write(0,'("p:      ",I5)'),p
          write(0,'("invert: ",I5)'),invert
          write(0,'("ostyle: ",I5)'),ostyle
          write(0,'("misfit: ",I5)'),misfit
          write(0,'("getdmp: ",I5)'),getdmp
          write(0,'("fact:   ",F12.6)'),fact
          write(0,*)
      endif
      call chkin(fltf,obsf)
      call chkout(ofile,p)

C----
C Read displacement, fault, material data
C----
      call readobs(obs,nobs,OBSMAX,obsf)
      call readflt(flt,nflt,FLTMAX,fltf)
      call readhaf(haf,haff)
      if (vrb.eq.1) then
          write(0,'("INPUT OBSERVATIONS")')
          write(0,'("------------------")')
          write(0,'("nobs: ",I5)'),nobs
          write(0,'(3A12,3A14)'),'x','y','z','ux','uy','uz'
          do 901 i = 1,nobs
              write(0,'(3F12.2,3F14.4)') (obs(i,j),j=1,6)
  901     continue
          write(0,*)
          write(0,'("INPUT FAULTS")')
          write(0,'("------------")')
          write(0,'("nflt: ",I5)'),nflt
          write(0,'(3A12,4A14)'),'x','y','z','str','dip','wid','len'
          do 902 i = 1,nflt
              write(0,'(3F12.2,4F14.4)') (flt(i,j),j=1,7)
  902     continue
          write(0,*)
          write(0,'("HALF-SPACE PARAMETERS")')
          write(0,'("---------------------")')
          write(0,'(3A14)'),'vp','vs','dens'
          write(0,'(3F14.4)') haf(1),haf(2),haf(3)
          write(0,*)
      endif

C----
C Convert geographical to Cartesian coordinates, if necessary
C----
      if (geo.eq.1) call geo2xy(obs,nobs,OBSMAX,flt,nflt,FLTMAX)

C----
C Compute Green's functions
C----
      call green(gf,obs,nobs,OBSMAX,flt,nflt,FLTMAX,haf)
      if (vrb.eq.1) then
          write(0,'("GREENS FUNCTIONS")')
          write(0,'("----------------")')
          write(0,'(2A4,6A14)') 'sta','flt','gfssx','gfssy','gfssz',
     1                                      'gfdsx','gfdsy','gfdsz'
          do 903 i = 1,nobs
              do 904 j = 1,nflt
                  write(0,'(2I4,6F14.6)') i,j,(gf(i,j,k),k=1,6)
  904         continue
  903     continue
          write(0,*)
      endif

C----
C Get best damping parameter
C----
      if (getdmp.eq.1) then
          print *,fact
          call getdamp(soln,gf,obs,nobs,OBSMAX,nflt,FLTMAX,damp0,smooth,
     1                 smoof,fact)
          damp = damp0
      endif

C----
C Run inversion
C----
      if (invert.eq.0) then
          call lstsqr(soln,gf,obs,nobs,OBSMAX,nflt,FLTMAX,damp,smooth,
     1                smoof)
      elseif (invert.eq.1) then
          call anneal(soln,gf,obs,nobs,OBSMAX,nflt,FLTMAX,damp,smooth,
     1                smoof,annf)
      else
          call usage('!! Error: inversion option does not exist')
      endif

C----
C Misfit
C----
      if (misfit.eq.1) then
          tot = 0.0d0
          do 905 i = 1,nobs
              pre(i,1) = 0.0d0
              pre(i,2) = 0.0d0
              pre(i,3) = 0.0d0
              do 906 j = 1,nflt
                  pre(i,1) = pre(i,1) + gf(i,j,1)*soln(j,1) 
     1                                             + gf(i,j,4)*soln(j,2)
                  pre(i,2) = pre(i,2) + gf(i,j,2)*soln(j,1)
     1                                             + gf(i,j,5)*soln(j,2)
                  pre(i,3) = pre(i,3) + gf(i,j,3)*soln(j,1)
     1                                             + gf(i,j,6)*soln(j,2)
  906         continue
  905     continue
          do 907 i = 1,nobs
              tot = tot + (pre(i,1)-obs(i,4))*(pre(i,1)-obs(i,4))
     1                   + (pre(i,2)-obs(i,5))*(pre(i,2)-obs(i,5))
     2                   + (pre(i,3)-obs(i,6))*(pre(i,3)-obs(i,6))
  907     continue
          pre(1,1) = 0.0d0
          do 908 i = 1,nflt
              pre(1,1) = pre(1,1) + soln(i,1)*soln(i,1)
     1                   + soln(i,2)*soln(i,2)
  908     continue
          print *,dsqrt(tot),dsqrt(pre(1,1))
      endif
C----
C Write output nicely
C----
      if (vrb.eq.1) then
          if (p.eq.0) then
              write(0,'("Solution (ss ds) in file",1X,A)') ofile
          else
              write(0,'("Solution (ss ds)")')
          endif
      endif
      call output(soln,nflt,FLTMAX,ofile,p,ostyle)

      END

C----------------------------------------------------------------------C

      SUBROUTINE chkin(fltf,obsf)
C----
C Check that input files exist
C----
      IMPLICIT none
      CHARACTER*80 fltf,obsf
      LOGICAL ex
      if (fltf.eq.'none') then
          write(*,*) '!! Error: Input fault geometry file unspecified'
          call usage('!! Use -flt FLTFILE to specify input file')
      else
          inquire(file=fltf,EXIST=ex)
          if (.not.ex) call usage('!! Error: no input file: '//
     1                                                       trim(fltf))
      endif
      if (obsf.eq.'none') then
          write(*,*) '!! Error: Observation file unspecified'
          call usage('!! Use -obs OBSFILE to specify input file')
      else
          inquire(file=obsf,EXIST=ex)
          if (.not.ex) call usage('!! Error: no input file: '//
     1                                                       trim(obsf))
      endif
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE chkout(ofile,p)
C----
C Check that output is defined
C----
      IMPLICIT none
      CHARACTER*80 ofile
      INTEGER p
      if (ofile.eq.'none'.and.p.eq.0) then
          write(*,*) '!! Error: Output file not specified'
          call usage('!! Use -o OFILE to define output file or -p')
      endif
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE readobs(obs,nobs,OBSMAX,obsf)
C----
C Read displacement observations from file
C----
      IMPLICIT none
      CHARACTER*80 obsf
      INTEGER nobs,OBSMAX,i,j
      REAL*8 obs(OBSMAX,6)
C Read in observation locations and displacement vectors
      call linect(nobs,obsf)
      if (nobs.gt.OBSMAX) call usage('!! Error: too many input data')
      open(unit=21,file=obsf,status='old')
      do 201 i = 1,nobs
          read(21,*) (obs(i,j),j=1,6) ! x y z ux uy uz (all meters)
  201 continue
      close(21)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE readflt(flt,nflt,FLTMAX,fltf)
C----
C Read model fault parameters from file
C----
      IMPLICIT none
      CHARACTER*80 fltf
      INTEGER nflt,FLTMAX,i,j
      REAL*8 flt(FLTMAX,7)
C Read model fault parameters and orientations
      call linect(nflt,fltf)
      if (nflt.gt.FLTMAX) call usage('!! Error: too many model faults')
      open(unit=22,file=fltf,status='old')
      do 202 i = 1,nflt
          read(22,*) (flt(i,j),j=1,7) ! x y z str dip wid len (m,deg)
  202 continue
      close(22)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE readhaf(haf,haff)
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
              write(*,*) '!! Warning: no input file: '//trim(haff)
              write(*,*) '!! Using default half-space parameters'
              haf(1) = 7.0d3
              haf(2) = haf(1)/dsqrt(3.0d0)
              haf(3) = 3.0d3
          else
              open(unit=23,file=haff,status='old')
              read(23,*) ch
              rewind(23)
              if (ch.eq.'Lame'.or.ch.eq.'lame') then
                  read(23,*) ch,lamda,mu
                  haf(3) = 3.0d3
                  haf(1) = dsqrt(mu/haf(3))
                  haf(2) = dsqrt((lamda+2.0d0*mu)/haf(3))
              else
                  read(23,*) haf(1),haf(2),haf(3)
              endif
              close(23)
          endif
      endif
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE geo2xy(obs,nobs,OBSMAX,flt,nflt,FLTMAX)
C----
C Convert longitude and latitude to x and y
C----
      IMPLICIT none
      INTEGER OBSMAX,FLTMAX
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

      SUBROUTINE green(gf,obs,nobs,OBSMAX,flt,nflt,FLTMAX,haf)
C----
C Compute and store static Green's functions
C----
      IMPLICIT none
      INTEGER OBSMAX,FLTMAX
      REAL*8 obs(OBSMAX,6),flt(FLTMAX,7)
      REAL*8 gf(OBSMAX,FLTMAX,6)
      INTEGER nobs,nflt
      REAL*8 haf(3)
      REAL*8 stlo,stla,stdp
      REAL*8 evlo,evla,evdp,str,dip,rak,slip,wid,len
      INTEGER i,j
      REAL*8 vp,vs,dens
C Compute vp, vs, dens from Lame parameters
      vp   = haf(1)
      vs   = haf(2)
      dens = haf(3)
C Unit slip for GF computation
      slip = 1.0d0
C GF for each flt-sta pair
      do 302 i = 1,nobs
          stlo = obs(i,1)
          stla = obs(i,2)
          stdp = obs(i,3)
          do 301 j = 1,nflt
              evlo = flt(j,1)
              evla = flt(j,2)
              evdp = flt(j,3)
              str  = flt(j,4)
              dip  = flt(j,5)
              wid  = flt(j,6)
              len  = flt(j,7)
              rak  = 0.0d0 ! strike-slip GF
              call calcdisp(gf(i,j,1),gf(i,j,2),gf(i,j,3),stlo,stla,
     1                      stdp,evlo,evla,evdp,str,dip,rak,slip,wid,
     2                      len,vp,vs,dens)
              rak  = 90.0d0 ! dip-slip GF
              call calcdisp(gf(i,j,4),gf(i,j,5),gf(i,j,6),stlo,stla,
     1                      stdp,evlo,evla,evdp,str,dip,rak,slip,wid,
     2                      len,vp,vs,dens)
  301     continue
  302 continue
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE calcdisp(ux,uy,uz,stlo,stla,stdp,evlo,evla,evdp,str,
     1                    dip,rak,slip,wid,len,vp,vs,dens)
C----
C Compute static displacement vector given a source-station pair
C All units are SI; angles inputs are in degrees
C----
      IMPLICIT none
      REAL*8 pi,d2r
      PARAMETER (pi=4.0d0*datan(1.0d0),d2r=pi/1.8d2)
      REAL*8 ux,uy,uz
      REAL*8 stlo,stla,stdp
      REAL*8 evlo,evla,evdp
      REAL*8 str,dip,rak
      REAL*8 slip,wid,len
      REAL*8 vp,vs,dens
      REAL*8 delx,dely,dist,az,x,y
      REAL*8 uxp,uyp,theta,uhor
C Distance and azimuth from source to station
      delx = stlo - evlo
      dely = stla - evla
      dist = dsqrt(delx*delx+dely*dely)
      az   = datan2(delx,dely) ! clockwise from north
C Rotate to x=along-strike, y=horizontal up-dip
      x = dist*( dcos(az-d2r*str))
      y = dist*(-dsin(az-d2r*str))
C Compute displacement
      call o92rect(uxp,uyp,uz,x,y,stdp,evdp,dip,rak,wid,len,slip,vp,vs,
     1             dens)
C Rotate back to original x-y coordinates
      theta = datan2(uyp,uxp)
      uhor = dsqrt(uxp*uxp+uyp*uyp)
      theta = d2r*str - theta
      ux = uhor*dsin(theta)
      uy = uhor*dcos(theta)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE lstsqr(soln,gf,obs,nobs,OBSMAX,nflt,FLTMAX,damp,smooth,
     1                  smoof)
C----
C Compute a least squares inversion for slip and rake
C----
      IMPLICIT none
      CHARACTER*80 smoof
      REAL*8 damp,smooth
      INTEGER OBSMAX,FLTMAX,dimx,dimy
      REAL*8 obs(OBSMAX,6)
      REAL*8 gf(OBSMAX,FLTMAX,6)
      REAL*8 A(3*OBSMAX+2*FLTMAX,2*FLTMAX)
      REAL*8 b(3*OBSMAX+2*FLTMAX,1)
      REAL*8 x(3*OBSMAX+2*FLTMAX,1)
      REAL*8 soln(FLTMAX,2)
      INTEGER i,j,nobs,nflt
      INTEGER n,m ! n=nrows, m=ncol
C----
C Load model matrix, A, and observation vector, b
C----
      do 402 i = 1,nobs
          do 401 j = 1,nflt
              A(i,j)             = gf(i,j,1) ! ssx
              A(i+nobs,j)        = gf(i,j,2) ! ssy
              A(i+2*nobs,j)      = gf(i,j,3) ! ssz
              A(i,j+nflt)        = gf(i,j,4) ! dsx
              A(i+nobs,j+nflt)   = gf(i,j,5) ! dsy
              A(i+2*nobs,j+nflt) = gf(i,j,6) ! dsz
  401     continue
  402 continue
      do 403 i = 1,nobs
          b(i,1)        = obs(i,4)
          b(i+nobs,1)   = obs(i,5)
          b(i+2*nobs,1) = obs(i,6)
  403 continue
C Change number of rows, columns
      n = 3*nobs
      m = 2*nflt
C----
C Add damping
C----
      call damping(A,b,n,m,damp)
C----
C Add smoothing
C----
      call smoothing(A,b,n,m,smooth,smoof)
C----
C Solve generalized least squares problem
C----
      dimx = 3*OBSMAX+2*FLTMAX
      dimy = 2*FLTMAX
      call lsqsolve(A,b,x,dimx,dimy,n,m)
      do 404 i = 1,nflt
          soln(i,1) = x(i,1)
          soln(i,2) = x(i+nflt,1)
  404 continue
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE damping(A,b,nobs,nflt,damp)
C----
C Minimize length of solution to Ax = b with weighting factor damp.
C----
      IMPLICIT none
      REAL*8 damp
      INTEGER OBSMAX,FLTMAX
      PARAMETER (OBSMAX=500,FLTMAX=500)
      REAL*8 A(3*OBSMAX+2*FLTMAX,2*FLTMAX)
      REAL*8 b(3*OBSMAX+2*FLTMAX,1)
      INTEGER i,j,nobs,nflt
      if (damp.le.0.0d0) return
      do 412 i = 1,nflt
          b(nobs+i,1) = 0.0d0
          do 411 j = 1,nflt
              if (i.eq.j) then
                  A(nobs+i,j) = damp*damp*1.0d0
              else
                  A(nobs+i,j) = 0.0d0
              endif
  411     continue
  412 continue
      nobs = nobs + nflt
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE smoothing(A,b,nobs,nflt,smooth,smoof)
C----
C Minimize roughness of solution to Ax = b with weighting factor smooth.
C----
      IMPLICIT none
      CHARACTER*80 smoof
      LOGICAL ex
      REAL*8 smooth
      INTEGER OBSMAX,FLTMAX
      PARAMETER (OBSMAX=500,FLTMAX=500)
      REAL*8 A(3*OBSMAX+2*FLTMAX,2*FLTMAX)
      REAL*8 b(3*OBSMAX+2*FLTMAX,1)
      INTEGER nsmoo,el,nnei,nei(FLTMAX)
      INTEGER i,j,nobs,nflt
      if (smooth.le.0.0d0) return
C Check for input smoothing file
      if (smoof.eq.'none') then
          write(*,*) '!! Error: No smoothing file specified'
          call usage('!! Use -smooth CONST FILE to specify input file')
      else
          inquire(file=smoof,EXIST=ex)
          if (.not.ex) call usage('!! Error: no input file: '//
     1                                                      trim(smoof))
      endif
      call linect(nsmoo,smoof)
      if (nsmoo.gt.nflt) then
          write(*,*) '!! Error: more faults to smooth than nflt'
          call usage('!! Check smoothing linking file')
      endif
C Initialize smoothing section of A and b arrays as zeros
      do 422 i = 1,nsmoo
          b(nobs+i,1) = 0.0d0
          do 421 j = 1,nflt
              A(nobs+i,j) = 0.0d0
  421     continue
  422 continue
C Read smoothing element linking file and load into A
      open(unit=41,file=smoof,status='old')
      do 424 i = 1,nsmoo
          read(41,*) el,nnei,(nei(j),j=1,nnei)
          ! Fault of interest
          A(nobs+i,el) = dble(nnei)*smooth*smooth
          A(nobs+i,nflt/2+el) = dble(nnei)*smooth*smooth ! nflt doubled at end of green3d
          do 423 j = 1,nnei
              A(nobs+i,nei(j)) = -1.0d0*smooth*smooth
              A(nobs+i,nflt/2+nei(j)) = -1.0d0*smooth*smooth
  423     continue
  424 continue
      close(41)
      nobs = nobs + nsmoo
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
      do 432 i = 1,nobs
          b(i,1) = bin(i,1)
          do 431 j = 1,npar
              A(i,j) = Ain(i,j)
  431     continue
  432 continue
C Solve generalized least squares problem using LAPACK tools
      call solve(A,x,b,nobs,npar)
C Print results
      do 433 i = 1,npar
          xout(i,1) = x(i,1)
  433 continue
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
      do 441 i = 1,nobs
          btmp(i,1) = b(i,1)
  441 continue
C Compute parameter array
      call dgels(trans,m,n,nrhs,a,lda,btmp,ldb,work,lwork,info)
      do 442 i = 1,npar
          x(i,1) = btmp(i,1)
  442 continue
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE getdamp(soln,gf,obs,nobs,OBSMAX,nflt,FLTMAX,damp,
     1                   smooth,smoof,fact)
      IMPLICIT none
      CHARACTER*80 smoof
      REAL*8 damp,smooth
      INTEGER OBSMAX,FLTMAX
      REAL*8 obs(OBSMAX,6)
      REAL*8 gf(OBSMAX,FLTMAX,6)
      REAL*8 soln(FLTMAX,2),misfit
      REAL*8 p(3)
      INTEGER i,j,nobs,nflt
      REAL*8 slip,slipmx,ddamp,misfit0,fact
C Compute misfit with damp = 0
      damp = 0.0d0
      call lstsqr(soln,gf,obs,nobs,OBSMAX,nflt,FLTMAX,damp,
     1            smooth,smoof)
      misfit0 = 0.0d0
      do 108 i = 1,nobs
          p(1) = 0.0d0
          p(2) = 0.0d0
          p(3) = 0.0d0
          do 107 j = 1,nflt
              p(1) = p(1) + gf(i,j,1)*soln(j,1)+gf(i,j,4)*soln(j,2)
              p(2) = p(2) + gf(i,j,2)*soln(j,1)+gf(i,j,5)*soln(j,2)
              p(3) = p(3) + gf(i,j,3)*soln(j,1)+gf(i,j,6)*soln(j,2)
  107     continue
          misfit0 = misfit0 + (p(1)-obs(i,4))*(p(1)-obs(i,4))
     1                    + (p(2)-obs(i,5))*(p(2)-obs(i,5))
     2                    + (p(3)-obs(i,6))*(p(3)-obs(i,6))
  108 continue
C Search for damping parameter
      ddamp = 1.0d0
  101 damp = damp + ddamp
          call lstsqr(soln,gf,obs,nobs,OBSMAX,nflt,FLTMAX,damp,
     1                smooth,smoof)
C Compute RMS obs-pre misfit
          misfit = 0.0d0
          do 102 i = 1,nobs
              p(1) = 0.0d0
              p(2) = 0.0d0
              p(3) = 0.0d0
              do 103 j = 1,nflt
                  p(1) = p(1) + gf(i,j,1)*soln(j,1)+gf(i,j,4)*soln(j,2)
                  p(2) = p(2) + gf(i,j,2)*soln(j,1)+gf(i,j,5)*soln(j,2)
                  p(3) = p(3) + gf(i,j,3)*soln(j,1)+gf(i,j,6)*soln(j,2)
  103         continue
              misfit = misfit + (p(1)-obs(i,4))*(p(1)-obs(i,4))
     1                        + (p(2)-obs(i,5))*(p(2)-obs(i,5))
     2                        + (p(3)-obs(i,6))*(p(3)-obs(i,6))
  102     continue
C Get maximum slip in model
          slipmx = 0.0d0
          do 105 i = 1,nflt
              slip = dsqrt(soln(i,1)*soln(i,1)+soln(i,2)*soln(i,2))
              if (slip.gt.slipmx) slipmx = slip
  105     continue
      if (misfit.le.misfit0*fact) then
          goto 101
      elseif(ddamp.gt.1.0d-4) then
          damp = damp - ddamp
          ddamp = ddamp*1.0d-1
          damp = damp - ddamp
          goto 101
      endif
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE anneal(soln,gf,obs,nobs,OBSMAX,nflt,FLTMAX,damp,smooth,
     1                  smoof,annf)
C----
C Perform a random search for the optimal solution, utilizing a
C simulated annealing algorithm.
C----
      IMPLICIT none
      REAL*8 pi,d2r
      PARAMETER (pi=4.0d0*atan(1.0d0),d2r=pi/1.8d2)
C Inversion variables
      INTEGER OBSMAX,FLTMAX
      REAL*8 obs(OBSMAX,6),gf(OBSMAX,FLTMAX,6)
      REAL*8 soln(FLTMAX,2),best(FLTMAX,2)
      INTEGER nobs,nflt
      CHARACTER*80 smoof
      REAL*8 smooth,damp
      INTEGER nsmoo,smarry(FLTMAX,8)
      REAL*8 misfit,length,rough
      REAL*8 best0(FLTMAX,2,50),obj0(50)
      INTEGER nsoln
C      CHARACTER*2 nsolnc
C Annealing-specific variables
      REAL*8 energy(500),obj,temp,objmin
      REAL*8 denom,test
      INTEGER i,j,jj
C Random number variables
      INTEGER time(3),timeseed,idum
      CHARACTER*8 timec
      EXTERNAL ran2
      REAL ran2,random
C Control file parameters
      CHARACTER*80 annf             ! Name of annealing control file 
      INTEGER      inityp           ! Initial solution type: 0=zero; 1=auto; 2=user
      REAL*8       slip0,rake0      ! Initial slip, rake values
      REAL*8       slip(FLTMAX,3)   ! Slip ranges: slipmn,slipmx,nslip,dslip
      REAL*8       rake(FLTMAX,3)   ! Rake ranges: rakemn,rakemx,nrake,drake
      REAL*8       temp0            ! Initial temperature
      REAL*8       cool             ! Cooling rate
      INTEGER      nitmax           ! Maximum iterations to run algorithm
      CHARACTER*10 algrthm          ! Type of algorithm to use
      INTEGER      reset            ! Step to reset temperature to 1
      INTEGER      tconst           ! Number of steps to keep temp at 1
C Other variables
      REAL*8 dslip(FLTMAX),drake(FLTMAX)

C Initialize random number generator
      call itime(time)
      write(timec,'(I0.2I0.2I0.2)') time(1),time(2),time(3)
      read(timec,*) timeseed
      if (timeseed.gt.0) timeseed = -timeseed
      idum = timeseed

C----
C Initialize output files
C----
      open(unit=51,file='anneal.log',status='unknown')
      open(unit=52,file='solutions.log',status='unknown')
      write(51,'(A)') 'Starting simulated annealing search'
      write(51,*)

C----
C Read control file
C----
      algrthm = 'metropolis'
      reset = 10000000
      tconst = -1
      write(51,'(A)') 'Setting up control parameters'
      call readannf(inityp,slip0,rake0,slip,rake,temp0,cool,nitmax,
     1              annf,nflt,FLTMAX,nsoln,algrthm,reset,tconst)
      write(51,'("  Initial model (0=zero,1=auto,2=user)")')
      write(51,'("      inityp:   ",I10)') inityp
      write(51,'("      slip0:    ",F10.2)') slip0
      write(51,'("      rake0:    ",F10.2)') rake0
      write(51,'("  Parameter ranges:")')
      write(51,'("      Fault slip ranges:")')
      write(51,'("         flt    slipmn    slipmx nslip")')
      do 501 i=1,nflt
          write(51,'(8X,I4,2F10.2,F6.0)') i,(slip(i,j),j=1,3)
  501 continue
      write(51,'("      Fault rake ranges:")')
      write(51,'("         flt    rakemn    rakemx nrake")')
      do 502 i=1,nflt
          write(51,'(8X,I4,2F10.2,F6.0)') i,(rake(i,j),j=1,3)
  502 continue
      write(51,'("  Cooling schedule:")')
      write(51,'("      temp0:    ",F10.4)') temp0
      write(51,'("      cool:     ",F10.6)') cool
      write(51,'("      nitmax:   ",I10)') nitmax
      write(51,'("  Number of solutions stored:")')
      write(51,'("      nsoln:    ",I10)') nsoln
      write(51,'("  Using algorithm ",A10)') algrthm
      write(51,'("  Reset period:   ",I10)') reset
      write(51,'("  Const temp steps",I10)') tconst
      write(51,*)

C----
C Set up array for computing smoothing
C----
      if (smooth.gt.0.0d0) then
          call readsmoof(smarry,nsmoo,smoof,nflt,FLTMAX)
          do 503 i = 1,nsmoo
              print *,smarry(i,1),smarry(i,2),
     1                                     (smarry(i,j),j=1,smarry(i,2))
  503     continue
      endif

C----
C Initial solution
C----
      write(51,'(A)') 'Loading initial solution'
      write(51,*)
      call model0(soln,gf,obs,nobs,OBSMAX,nflt,FLTMAX,inityp,slip0,
     1            rake0)
C Check that initial solution lies in range [slipmn,slipmx], [rakemn,rakemx]
      do 558 i = 1,nflt
          if (soln(i,1).lt.slip(i,1).or.soln(i,1).gt.slip(i,2)) then
              write(*,*) '!! Warning: initial slip:',soln(i,1)
              write(*,*) '!! outside of range [slipmn,slipmx]'
              write(*,*) '!! Choosing random slip in [',
     1                                         slip(i,1),slip(i,2),']'
              soln(i,1) = slip(i,1) + ran2(idum)*(slip(i,2)-slip(i,1))
          endif
          if (soln(i,2).lt.rake(i,1).or.soln(i,2).gt.rake(i,2)) then
              write(*,*) '!! Warning: initial rake:',rake(i,1)
              write(*,*) '!! outside of range [rakemn,rakemx]'
              write(*,*) '!! Choosing random rake in [',
     1                                         rake(i,1),rake(i,2),']'
              soln(i,2) = rake(i,1) + ran2(idum)*(rake(i,2)-rake(i,1))
          endif
  558 continue
C Initialize objective function
      do 505 j = 1,nsoln
          obj0(j) = 1d10
  505 continue

C----
C Initial solution objective function
C----
      write(51,'(A)') 'Computing initial objective function'
      call calcobj(obj,misfit,length,rough,soln,gf,obs,nobs,OBSMAX,nflt,
     1             FLTMAX,damp,smooth,smarry,nsmoo)
      objmin = obj
      obj0(1) = obj
      do 506 i = 1,nflt
          best(i,1) = soln(i,1)
          best(i,2) = soln(i,2)
          best0(i,1,1) = soln(i,1)
          best0(i,2,1) = soln(i,2)
  506 continue
      write(51,1001) 
      write(51,1002) obj,misfit,damp,length,smooth,rough
      write(51,*)
 1001 format("   objective =      misfit + (        damp*      length",
     1       ") + (      smooth*       rough)")
 1002 format(F12.6," =",F12.6," + (",F12.6,"*",F12.6,") + (",F12.6,"*",
     1       F12.6,")")

C Record solutions in solutions.log
 1003 format('Iteration: ',I4,4X,'Temperature: ',F9.6,4X,'Objective=',
     1       F12.6)
 1004 format("  flt      slip      rake")
 1005 format(I5,F10.4,F10.2)
      write(52,1003) 0,1e10,obj
      write(52,1004)
      do 504 i=1,nflt
          write(52,1005) i,soln(i,1),soln(i,2)
  504 continue
      write(52,*)

C----
C Run annealing search
C----
      write(51,'(A)') 'Starting annealing'
C Store some derived arrays
      do 507 i = 1,nflt
          dslip(i) = (slip(i,2)-slip(i,1))/slip(i,3)
          drake(i) = (rake(i,2)-rake(i,1))/rake(i,3)
  507 continue
C Iterate and search
      temp = temp0
      do 520 i = 1,nitmax
          do 514 j = 1,nflt
C------------ Metropolis-Hastings: perturb all variables by up to one increment
              if (algrthm.eq.'metropolis') then
  525             soln(j,1) = best(j,1) +
     1                                 4.0d0*dslip(j)*(ran2(idum)-0.5d0)
                  if (soln(j,1).lt.slip(j,1)
     1                              .or.soln(j,1).gt.slip(j,2)) goto 525
  526             soln(j,2) = best(j,2) +
     1                                 4.0d0*drake(j)*(ran2(idum)-0.5d0)
                  if (soln(j,2).lt.rake(j,1)
     1                              .or.soln(j,2).gt.rake(j,2)) goto 526
C------------ Gibbs sampler: create variable PDFs; directly sample
              elseif (algrthm.eq.'gibbs') then
                  denom = 0.0d0
                  do 508 jj = 1,int(slip(j,3)+0.1)
                      soln(j,1) = slip(j,1) + (dble(jj)-0.5d0)*dslip(j)
                      call calcobj(obj,misfit,length,rough,soln,gf,obs,
     1                             nobs,OBSMAX,nflt,FLTMAX,damp,smooth,
     2                             smarry,nsmoo)
                      if (obj/temp.gt.100.0d0) then
                          energy(jj) = 0.0d0
                      else
                          energy(jj) = dexp(-obj/temp)
                      endif
                      denom = denom + energy(jj)
  508             continue
                  test = 0.0d0
                  do 509 jj = 1,int(slip(j,3)+0.1)
                      test = test + energy(jj)/denom
                      if (ran2(idum).lt.test) then ! Sample PDF
                          goto 510
                      endif
  509             continue
  510             soln(j,1) = slip(j,1) + (dble(jj)-0.5d0)*dslip(j)
                  denom = 0.0d0
                  do 511 jj = 1,int(rake(j,3)+0.1)
                      soln(j,2) = rake(j,1) + (dble(jj)-0.5d0)*drake(j)
                      call calcobj(obj,misfit,length,rough,soln,gf,obs,
     1                             nobs,OBSMAX,nflt,FLTMAX,damp,smooth,
     2                             smarry,nsmoo)
                      if (obj/temp.gt.100.0d0) then
                          energy(jj) = 0.0d0
                      else
                          energy(jj) = dexp(-obj/temp)
                      endif
                      denom = denom + energy(jj)
  511             continue
                  test = 0.0d0
                  do 512 jj = 1,int(rake(j,3)+0.1)
                      test = test + energy(jj)/denom
                      if (ran2(idum).lt.test) then ! Sample PDF
                          goto 513
                      endif
  512             continue
  513             soln(j,2) = rake(j,1) + (dble(jj)-0.5d0)*drake(j)
              endif
  514     continue

C         Compute objective function for model in iteration i
          call calcobj(obj,misfit,length,rough,soln,gf,obs,nobs,OBSMAX,
     1                 nflt,FLTMAX,damp,smooth,smarry,nsmoo)
C         Record results in log files
          write(51,1003) i,temp,obj
          write(52,1003) i,temp,obj
          write(52,1004)
          do 515 j = 1,nflt
              write(52,1005) j,soln(j,1),soln(j,2)
  515     continue
          write(52,*)

C         Find the worst saved solution
          test = 0.0d0
          do 517 j = 1,nsoln
              if (obj0(j).gt.test) then
                  test = obj0(j) 
                  jj = j
              endif
  517     continue
          j = jj
C         Save solution in best models array if in top NSOLN objective values
          if (obj.lt.obj0(j)) then
              do 516 jj = 1,nflt
                  best0(jj,1,j) = soln(jj,1)
                  best0(jj,2,j) = soln(jj,2)
  516         continue
              obj0(j) = obj
          endif
C         If solution is better than previous, save it
          if (obj.lt.objmin) then
              do 519 j = 1,nflt
                  best(j,1) = soln(j,1)
                  best(j,2) = soln(j,2)
                  objmin = obj
  519         continue
C         When using Metropolis-Hastings, allow solution chance to get worse
          elseif (algrthm.eq.'metropolis') then
              test = dexp((objmin-obj)/temp)
              random = ran2(idum)
              if (random.lt.test) then
                  do 527 j = 1,nflt
                      best(j,1) = soln(j,1)
                      best(j,2) = soln(j,2)
                      objmin = obj
  527             continue
              endif
          endif
C         Cool off
          if (mod(i,reset).lt.tconst) then
              temp = temp
          elseif (temp.gt.1e-4) then
              temp = temp*cool
          else
              temp = temp*0.995d0
          endif
          if (mod(i,reset).eq.0) then
              temp = 1
          endif
  520 continue
      write(51,*)
      write(51,'("Finished annealing")')
      close(51)
      close(52)

C----
C Write nsoln best solutions to file best_nsoln.dat
C----
C      write(nsolnc,'(I0.2)') nsoln
C      open(unit=53,file='best_'//nsolnc//'.dat',status='unknown')
      open(unit=53,file='best0.dat',status='unknown')
      write(53,'(50F12.6)') (obj0(j),j=1,nsoln)
      do 521 i = 1,nflt
          write(53,'(50F12.6)') (best0(i,1,j),j=1,nsoln)
  521 continue
      do 522 i = 1,nflt
          write(53,'(50F12.6)') (best0(i,2,j),j=1,nsoln)
  522 continue
      close(53)

C----
C Return final solution (NOT NECESSARILY LOWEST OBJECTIVE VALUE!)
C----
      do 523 i = 1,nflt
          soln(i,1) = best(i,1)*dcos(best(i,2)*d2r)
          soln(i,2) = best(i,1)*dsin(best(i,2)*d2r)
  523 continue
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE readannf(inityp,slip0,rake0,slip,rake,temp0,cool,
     1                    nitmax,annf,nflt,FLTMAX,nsoln,algrthm,reset,
     2                    tconst)
C----
C Read annealing control file if it exists
C----
      IMPLICIT none
      INTEGER nflt,FLTMAX
      CHARACTER*80 annf
      INTEGER inityp,nitmax,reset,tconst
      INTEGER nsoln
      REAL*8 slip0,rake0,slip(FLTMAX,3),rake(FLTMAX,3),temp0,cool
      LOGICAL ex
      CHARACTER*120 line
      CHARACTER*1 dum
      INTEGER i,j,ns,nr
      CHARACTER*10 algrthm
C If annealing file not specified, use some default values
  551 if (annf.eq.'none') then
          inityp = 0
          slip0 = 0.0d0
          rake0 = 0.0d0
          do 552 i = 1,nflt
              slip(i,1) = 0.0d0
              slip(i,2) = 25.0d0
              slip(i,3) = 10.0d0
              rake(i,1) = -180.0d0
              rake(i,2) = 180.0d0
              rake(i,3) = 30.0d0
  552     continue
          temp0 = 1.0d0
          cool = 0.900d0
          nitmax = 250
          nsoln = 1
          return
      endif
C Check if annealing control file exists
      inquire(file=annf,exist=ex)
      if (.not.ex) then
          write(*,*) '!! Warning: no annealing file '//trim(annf)
          write(*,*) 'Using default parameters'
          annf = 'none'
          goto 551
      endif
C Initialize slip, rake arrays
      do 553 i = 1,nflt
          slip(i,1) = 0.0d0
          slip(i,2) = 25.0d0
          slip(i,3) = 25.0d0
          rake(i,1) = -180.0d0
          rake(i,2) = 180.0d0
          rake(i,3) = 36.0d0
  553 continue
C Read annealing control file
      open(unit=59,file=annf,status='old')
      read(59,*) line ! auto/user/zero
      if (line.eq.'zero') then
          inityp = 0
      elseif (line.eq.'auto') then
          inityp = 1
      elseif (line.eq.'user') then
          inityp = 2
      else
          write(*,*) '!! Error: no initialize option '//trim(line)
          call usage('!! Use auto, user, or zero')
      endif
      read(59,*) slip0,rake0
      ns = 0
      nr = 0
  554 read(59,'(A)') line
      i = index(line,'s')
      if (i.eq.1) then
          read(line,*) dum,j,slip(j,1),slip(j,2),slip(j,3)
          ns = ns + 1
          goto 554
      endif
      i = index(line,'r')
      if (i.eq.1) then
          read(line,*) dum,j,rake(j,1),rake(j,2),rake(j,3)
          nr = nr + 1
          goto 554
      endif
      if (ns.eq.0.or.nr.eq.0) then
          write(*,*) '!! Error: must define slip and rake limits'
          write(*,*) '!! s flt slpmn slipmx nslip'
          call usage('!! r flt slpmn rakemx nrake')
      endif
C If only one fault parameter limits is defined, use it for all
      if (ns.eq.1.and.nflt.ge.2) then
          do 555 i = 2,nflt
              slip(i,1) = slip(1,1)
              slip(i,2) = slip(1,2)
              slip(i,3) = slip(1,3)
  555     continue
      endif
      if (nr.eq.1.and.nflt.ge.2) then
          do 556 i = 2,nflt
              rake(i,1) = rake(1,1)
              rake(i,2) = rake(1,2)
              rake(i,3) = rake(1,3)
  556     continue
      endif
      read(line,*) temp0,cool,nitmax
      nsoln = 1
      read(59,*) nsoln
      read(59,*) algrthm
      read(59,*) reset
      read(59,*) tconst
C  557 continue
      close(59)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE readsmoof(smarry,nsmoo,smoof,nflt,FLTMAX)
      IMPLICIT none
      CHARACTER*80 smoof
      INTEGER nflt,FLTMAX,i,j
      INTEGER nsmoo,smarry(FLTMAX,8)
      call linect(nsmoo,smoof)
      if (nsmoo.gt.nflt) then
          call usage('!! Error: number of lines in smooth file > nflt')
      endif
C Read smoothing element linking file
      open(unit=63,file=smoof,status='old')
      do 561 i = 1,nsmoo
          read(63,*) smarry(i,1),smarry(i,2),
     1                                     (smarry(i,j),j=1,smarry(i,2))
  561 continue
      close(63)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE model0(soln,gf,obs,nobs,OBSMAX,nflt,FLTMAX,inityp,
     1                  slip0,rake0)
C----
C Set up the initial solution
C----
      IMPLICIT none
      REAL*8 pi,d2r
      PARAMETER (pi=4.0d0*atan(1.0d0),d2r=pi/1.8d2)
      INTEGER nobs,OBSMAX,nflt,FLTMAX,inityp,i,j
      REAL*8 obs(OBSMAX,6),gf(OBSMAX,FLTMAX,6),soln(FLTMAX,2)
      REAL*8 pre(3),op,pp,slip0,rake0,ss,ds
      if (inityp.eq.0) then
          slip0 = 0.0d0
          rake0  = 0.0d0
      elseif (inityp.eq.1) then
          ss = dcos(rake0*d2r)
          ds = dsin(rake0*d2r)
          op = 0.0d0 ! sum of dot product obs.pre
          pp = 0.0d0 ! sum of dot product pre.pre
          do 572 i = 1,nobs
              do 571 j = 1,nflt
                  pre(1) = ss*gf(i,j,1) + ds*gf(i,j,4)
                  pre(2) = ss*gf(i,j,2) + ds*gf(i,j,5)
                  pre(3) = ss*gf(i,j,3) + ds*gf(i,j,6)
                  op = op + obs(i,4)*pre(1) + obs(i,5)*pre(2) +
     1                                                   obs(i,6)*pre(3)
                  pp = pp +   pre(1)*pre(1) +   pre(2)*pre(2) +
     1                                                     pre(3)*pre(3)
  571         continue
  572     continue
          slip0 = op/pp
      elseif (inityp.eq.2) then
          ! slip0 and rake0 already defined
      else
          call usage('!! Error: no annealing option '//char(inityp))
      endif
      do 573 i = 1,nflt
          soln(i,1) = slip0
          soln(i,2) = rake0
  573 continue
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE calcobj(obj,misfit,length,rough,soln,gf,obs,nobs,
     1                   OBSMAX,nflt,FLTMAX,damp,smooth,smarry,nsmoo)
      IMPLICIT none
      REAL*8 obj
      REAL*8 misfit,length,rough,smooth,damp
      INTEGER nobs,OBSMAX,nflt,FLTMAX
      REAL*8 obs(OBSMAX,6),gf(OBSMAX,FLTMAX,6),soln(FLTMAX,2)
      INTEGER nsmoo,smarry(FLTMAX,8)
      call calcmisfit(misfit,soln,gf,obs,nobs,OBSMAX,nflt,FLTMAX)
      misfit = dsqrt(misfit/nflt) ! RMS misfit
      length = 0.0d0
      rough = 0.0d0
      if (damp.gt.0.0d0) then
          call calclength(length,soln,nflt,FLTMAX)
          length = damp*length/nflt ! Average length
      endif
      if (smooth.gt.0.0d0) then
          call calcrough(rough,soln,nflt,FLTMAX,smarry,nsmoo)
          rough = smooth*rough/nflt ! Average roughness
      endif
      obj = misfit + damp*length + smooth*rough
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE calcmisfit(misfit,soln,gf,obs,nobs,OBSMAX,nflt,FLTMAX)
      IMPLICIT none
      REAL*8 pi,d2r
      PARAMETER (pi=4.0d0*datan(1.0d0),d2r=pi/1.8d2)
      INTEGER nobs,OBSMAX,nflt,FLTMAX
      REAL*8 obs(OBSMAX,6),gf(OBSMAX,FLTMAX,6)
      REAL*8 soln(FLTMAX,2)
      REAL*8 pre(3)
      REAL*8 misfit
      INTEGER i,j
      REAL*8 cr,sr,ss(FLTMAX),ds(FLTMAX)
C Store strike-slip and dip-slip components of faults
      do 581 i = 1,nflt
          cr = dcos(soln(i,2)*d2r)
          sr = dsin(soln(i,2)*d2r)
          ss(i) = soln(i,1)*cr
          ds(i) = soln(i,1)*sr
  581 continue
C Compute misfit
      misfit = 0.0d0
      do 583 i = 1,nobs
          pre(1) = 0.0d0
          pre(2) = 0.0d0
          pre(3) = 0.0d0
          do 582 j = 1,nflt
              pre(1) = pre(1) + gf(i,j,1)*ss(j) + gf(i,j,4)*ds(j)
              pre(2) = pre(2) + gf(i,j,2)*ss(j) + gf(i,j,5)*ds(j)
              pre(3) = pre(3) + gf(i,j,3)*ss(j) + gf(i,j,6)*ds(j)
  582     continue
          misfit = misfit + (obs(i,4)-pre(1))*(obs(i,4)-pre(1))
     1                    + (obs(i,5)-pre(2))*(obs(i,5)-pre(2))
     2                    + (obs(i,6)-pre(3))*(obs(i,6)-pre(3))
  583 continue
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE calclength(length,soln,nflt,FLTMAX)
      IMPLICIT none
      INTEGER nflt,FLTMAX,i
      REAL*8 soln(FLTMAX,2),length
      length = 0.0d0
      do 584 i = 1,nflt
          length = length + soln(i,1)
  584 continue
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE calcrough(rough,soln,nflt,FLTMAX,smarry,nsmoo)
      IMPLICIT none
      REAL*8 pi,d2r
      PARAMETER (pi=4.0d0*datan(1.0d0),d2r=pi/1.8d2)
      INTEGER nflt,FLTMAX,nsmoo,smarry(FLTMAX,8),i,j,jj
      REAL*8 soln(FLTMAX,2),rough,r(2),ss,ds
      r(1) = 0.0d0
      r(2) = 0.0d0
      do 586 i = 1,nsmoo
          ss = soln(i,1)*dcos(soln(i,2)*d2r)
          ds = soln(i,1)*dsin(soln(i,2)*d2r)
          r(1) = r(1) + smarry(i,2)*ss
          r(2) = r(2) + smarry(i,2)*ds
          do 585 j = 1,smarry(i,2)
              jj = smarry(i,j+2)
              ss = soln(jj,1)*dcos(soln(jj,2)*d2r)
              ds = soln(jj,1)*dsin(soln(jj,2)*d2r)
              r(1) = r(1) - ss
              r(2) = r(2) - ds
  585     continue
  586 continue
      rough = r(1) + r(2)
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE chksol(x,nflt)
      IMPLICIT none
      INTEGER OBSMAX
      PARAMETER (OBSMAX=2500)
      REAL*8 x(OBSMAX,1)
      INTEGER nflt,i
      do 611 i = 1,nflt
          if (x(i,1).gt.100.0d0) then
              write(*,*) '!! Error: found slip larger than 100 m'
              write(*,*) '!! Check inputs and units'
              write(*,*) '!! Try adding damping/smoothing'
              stop
          endif
  611 continue
      RETURN
      END

C----------------------------------------------------------------------C

      SUBROUTINE output(soln,nflt,FLTMAX,ofile,p,ostyle)
      IMPLICIT none
      INTEGER FLTMAX
      REAL*8 soln(FLTMAX,2)
      INTEGER nflt,p,i,ostyle
      CHARACTER*80 ofile,frmt
      if (ostyle.eq.0) then
          frmt = '(2F14.6)'
      elseif (ostyle.eq.1) then
          frmt = '(2F20.10)'
      elseif (ostyle.eq.10) then
          frmt = '(2E16.6)'
      elseif (ostyle.eq.11) then
          frmt = '(2E20.10)'
      else
          frmt = '(2F14.6)'
      endif
      if (p.eq.0) then
          open(unit=23,file=ofile,status='unknown')
          do 601 i = 1,nflt
              write(23,frmt) soln(i,1),soln(i,2)
  601     continue
          close(23)
      else
          do 602 i = 1,nflt
              write(*,frmt) soln(i,1),soln(i,2)
  602     continue
      endif
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

      SUBROUTINE gcmdln(fltf,obsf,smoof,smooth,damp,ofile,p,geo,haff,
     1                  invert,annf,ostyle,misfit,getdmp,fact)
      IMPLICIT NONE
      CHARACTER*80 fltf,obsf,tag,smoof,ofile,haff,annf
      INTEGER i,narg,indx
      REAL*8 smooth,damp,fact
      INTEGER geo,p,invert,vrb,ostyle,misfit,getdmp
      COMMON /VERBOSE/ vrb
      fltf = 'none'
      obsf = 'none'
      smoof = 'none'
      ofile = 'none'
      haff  = 'none'
      annf  = 'none'
      smooth = -1.0d0
      damp = -1.0d0
      geo = 0
      p = 0
      invert = 0
      vrb = 0
      ostyle = 0
      indx = 0
      misfit = 0
      getdmp = 0
      fact = 1.5d0
      narg = iargc()
      if (narg.eq.0) call usage('!! Error: no command line arguments')
      i = 0
  101 i = i + 1
      if (i.gt.narg) goto 102
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
              indx = index(tag,'-')
              if (indx.ne.0.or.i.gt.narg) then
                  call usage('!! Error: must define smoothing file')
              endif
              read(tag,'(BN,F10.0)') smooth
              i = i + 1
              call getarg(i,smoof)
          elseif (tag(1:5).eq.'-damp') then
              i = i + 1
              call getarg(i,tag)
              read(tag,'(BN,F10.0)') damp
          elseif (tag(1:4).eq.'-lsq') then
              invert = 0
          elseif (tag(1:7).eq.'-anneal') then
              invert = 1
              i = i + 1
              if (i.gt.narg) goto 102
              call getarg(i,tag)
              indx = index(tag,'-')
              if (indx.eq.0) then
                  annf = tag
              else
                  i = i - 1
              endif
          elseif (tag(1:4).eq.'-vrb'.or.tag(1:5).eq.'-verb'.or.
     1                      tag(1:8).eq.'-verbose') then
              vrb = 1
          elseif (tag(1:2).eq.'-o') then
              i = i + 1
              call getarg(i,ofile)
          elseif (tag(1:4).eq.'-geo') then
              geo = 1
          elseif (tag(1:4).eq.'-haf') then
              i = i + 1
              call getarg(i,haff)
          elseif (tag(1:8).eq.'-declong') then
              ostyle = 1
          elseif (tag(1:4).eq.'-dec') then
              ostyle = 0
          elseif (tag(1:8).eq.'-scilong') then
              ostyle = 11
          elseif (tag(1:4).eq.'-sci') then
              ostyle = 10
          elseif (tag(1:7).eq.'-misfit') then
              misfit = 1
          elseif (tag(1:8).eq.'-getdamp') then
              getdmp = 1
              i = i + 1
              if (i.gt.narg) goto 102
              call getarg(i,tag)
              indx = index(tag,'-')
              if (indx.eq.0) then
                  read(tag,*) fact
                  if (fact.lt.1.0d0) then
                      write(0,*) '!! Warning: factor must be >= 1'
                      write(0,*) '!! Setting factor = 1'
                      fact = 1.0d0
                  endif
              else
                  i = i - 1
              endif
          elseif (tag(1:2).eq.'-h') then
              call usage(' ')
          elseif (tag(1:2).eq.'-d') then
              call usage('D')
          elseif (tag(1:2).eq.'-p') then
              p = 1
          else
              call usage('!! Error: No option '//trim(tag))
          endif
          goto 101
  102 continue
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
     1 'Usage: fltinv -obs OBSFILE -flt FLTFILE -o OFILE [-p] [-geo]'
      write(*,*)
     1 '              [-haf HAFSPC] [-damp CONST] [-smooth CONST FILE]'
      write(*,*)
     1 '              [-anneal [ANEALFILE]] [-dec[long]|-sci[long]]'
      write(*,*)
     1 '              [-misfit] [-getdamp] [-h|-d]'
      write(*,*)
      write(*,*)
     1 '-obs OBSFILE        Displacement observations'
      if (d.eq.1) then
          write(*,*)
     1 '                       x(m) y(m) z(m) ux(m) uy(m) uz(m)'
          write(*,*)
      endif
      write(*,*)
     1 '-flt FLTFILE        Input fault location and geometry'
      if (d.eq.1) then
          write(*,*)
     1 '                       x(m) y(m) z(m) str(deg) dip(deg) ',
     2                         'wid(m) len(m)'
          write(*,*)
      endif
      write(*,*)
     1 '-o OFILE            Output slip values'
      if (d.eq.1) then
          write(*,*)
     1 '                       str_slip(m) dip_slip(m)'
          write(*,*)
      endif
      write(*,*)
     1 '-p                  Print output to standard out (overrides -o)'
      if (d.eq.1) then
          write(*,*)
      endif
      write(*,*)
     1 '-geo                Use geographic coordinates'
      if (d.eq.1) then
          write(*,*)
     1 '                       X->E, Y->N'
          write(*,*)
      endif
      write(*,*)
     1 '-haf HAFSPC         Half-space parameter file'
      if (d.eq.1) then
          write(*,*)
     1 '                       vp(m/s) vs(m/s) dens(kg/m^3)'
          write(*,*)
     1 '                       Lame lamda(Pa) mu(Pa)'
          write(*,*)
      endif
      write(*,*)
     1 '-damp CONST         Minimize length of model, with weight CONST'
      if (d.eq.1) then
          write(*,*)
      endif
      write(*,*)
     1 '-smooth CONST FILE  Smooth model, with weight CONST, linking ',
     2                         'FILE'
      if (d.eq.1) then
          write(*,*)
     1 '                        Linking file format: flt nghbrs ',
     2                                  'nghbr_1-N'
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
     1 '-anneal [ANEALFILE] Implement a simulated annealing approach ',
     2                          '(NOT WORKING YET)'
      if (d.eq.1) then
          write(*,*)
     1 '                        Optional annealing control file:'
          write(*,*)
     1 '                        zero|auto|user             ',
     2                            '(initial model type)'
          write(*,*)
     1 '                        slip0 rake0                ',
     2                            '(initial model parameters)'
          write(*,*)
     1 '                        s flt slipmn slipmx nslip  ',
     2                            '(slip ranges, up to nflt)'
          write(*,*)
     1 '                        r flt rakemn rakemx nrake  ',
     2                            '(rake ranges, up to nflt)'
          write(*,*)
     1 '                        temp0 cool nitmax          ',
     2                                     '(cooling parameters)'
          write(*,*)
     1 '                        nsoln                      ',
     2                                   '(number of solutions to keep)'
          write(*,*)
      endif
      write(*,*)
     1 '-dec|-declong       Decimal or long decimal output ',
     2                                                '(F14.6|F20.10)'
      write(*,*)
     1 '-sci|-scilong       Scientific or long scientific output ',
     2                                               '(E16.6|E20.10)'
      if (d.eq.1) then
          write(*,*)
     1 '                        Format output in OFILE'
          write(*,*)
      endif
      write(*,*)
     1 '-misfit             Print least squares misfit to stdout'
      if (d.eq.1) then
          write(*,*)
      endif
      write(*,*)
     1 '-getdamp [FACT]     Determine good damping parameter'
      if (d.eq.1) then
          write(*,*)
     1 '                        Compute misfit (E0) for DAMP=0, find ',
     2                           'largest DAMP s.t. E<=E0*FACT'
          write(*,*)
     1 '                        Default: FACT=1.5'
          write(*,*)
      endif
      write(*,*)
     1 '-h                  Short online help'
      write(*,*)
     1 '-d                  Detailed online help'
      write(*,*)
      STOP
      END

