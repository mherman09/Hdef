!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!----- LINEAR LEAST-SQUARES INVERSION SUBROUTINES -------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

subroutine lsqr_invert()
!----
! Invert displacements and pre-stresses for fault slip using a linear least squares approach
!----
use command_line, only : verbosity
use arrays, only : nfaults, fault_slip
use io, only : stderr
implicit none
! Local variables
double precision, dimension(:,:), allocatable :: A, b, x
integer :: nrows, ncols, ptr_disp, ptr_stress, ptr_damp, ptr_smooth

if (verbosity.ge.1) then
    write(stderr,'("Starting lsqr_invert()")')
endif

! Determine dimensions of model matrix
call lsqr_calc_array_dimensions(nrows,ncols,ptr_disp,ptr_stress,ptr_damp,ptr_smooth)

! Allocate memory to model matrix, constraint vector, and solution vector
if (.not.allocated(A)) then
    allocate(A(nrows,ncols))
endif
if (.not.allocated(b)) then
    allocate(b(nrows,1))
endif
if (.not.allocated(x)) then
    allocate(x(ncols,1))
endif
A = 0.0d0
b = 0.0d0
x = 0.0d0

! Load model matrix, A, and constraint vector, b
call lsqr_load_arrays(A,b,nrows,ncols,ptr_disp,ptr_stress,ptr_damp,ptr_smooth)

! Solve generalized least squares problem
call lsqr_solve(x,A,b,nrows,ncols)

if (.not.allocated(fault_slip)) then
    allocate(fault_slip(nfaults))
endif
fault_slip = x(:,1)

if (verbosity.ge.1) then
    write(stderr,'("Finished lsqr_invert()")')
    write(stderr,*)
endif

return
end

!--------------------------------------------------------------------------------------------------!

subroutine lsqr_calc_array_dimensions(nrows,ncols,ptr_disp,ptr_stress,ptr_damp,ptr_smooth)
use command_line
use arrays
use io
implicit none
! I/O variables
integer :: nrows, ncols, ptr_disp, ptr_stress, ptr_damp, ptr_smooth

!  Initialize inversion array dimensions
nrows = 0
ncols = 0

! Initialize array pointers
ptr_disp = 0
ptr_stress = 0
ptr_damp = 0
ptr_smooth = 0

! Row dimensions corresponding to displacement data
if (displacement_file.ne.'none') then
    if (disp_comp.eq.'a') then
        nrows = nrows + 3*ndisplacements
    elseif (disp_comp.eq.'h') then
        nrows = nrows + 2*ndisplacements
    elseif (disp_comp.eq.'v') then
        nrows = nrows + 1*ndisplacements
    else
        call usage('!! Error: disp_comp must be "a", "h", or "v"')
    endif
    ptr_disp = 1
endif

! Row dimensions corresponding to pre-stress data
if (prestress_file.ne.'none') then
    nrows = nrows + 2*nfaults
    if (ptr_disp.eq.1) then
        ptr_stress = nrows + 1
    else
        ptr_stress = 1
    endif
endif

! Column dimensions corresponding to fault slip
if (rake_file.eq.'none') then
    ncols = 2*nfaults
else
    ncols = nfaults
endif

! Add rows for damping (one additional row per model DOF in inversion)
if (damping_constant.gt.0.0d0) then
    ptr_damp = nrows + 1
    nrows = nrows + ncols
endif

! Add rows for smoothing
if (smoothing_constant.gt.0.0d0) then
    ptr_smooth = nrows + 1
    nrows = nrows + nsmooth
endif

if (verbosity.ge.2) then
    write(stderr,'(A)') 'Model matrix dimensions'
    write(stderr,'(A,I8)') 'nrows:   ', nrows
    write(stderr,'(A,I8)') 'ncolumns:', ncols
    write(stderr,*)
endif

return
end

!--------------------------------------------------------------------------------------------------!

subroutine lsqr_load_arrays(A,b,nrows,ncols,ptr_disp,ptr_stress,ptr_damp,ptr_smooth)
use command_line
use arrays
use io
implicit none
! I/O variables
integer :: nrows, ncols, ptr_disp, ptr_stress, ptr_damp, ptr_smooth
double precision :: A(nrows,ncols)
double precision :: b(nrows,1)
! Local variables
integer :: i, j
character(len=256) :: fmt

! Load displacement data
if (displacement_file.ne.'none') then
    call lsqr_load_disp(A,b,nrows,ncols)
endif

! Load pre-stress data
if (prestress_file.ne.'none') then
    call lsqr_load_stress(A,b,nrows,ncols,ptr_stress)
endif

! Minimize model length
if (damping_constant.gt.0.0d0) then
    call lsqr_load_damp(A,b,nrows,ncols,damping_constant,ptr_damp)
endif

! Minimize model smoothness
if (smoothing_constant.gt.0.0d0) then
    call lsqr_load_smooth(A,b,nrows,ncols,smoothing_constant,ptr_smooth)
endif

if (verbosity.ge.3) then
    write(stderr,'("A =")')
    write(fmt,9001) ncols
    do i = 1,nrows
        if (i.eq.ptr_disp) write(stderr,'("displacements")')
        if (i.eq.ptr_stress) write(stderr,'("pre-stresses")')
        if (i.eq.ptr_damp) write(stderr,'("damping")')
        if (i.eq.ptr_smooth) write(stderr,'("smoothing")')
        write(0,fmt) (A(i,j),j=1,ncols)
    enddo
    write(stderr,*)
    write(stderr,'("b =")')
    do i = 1,nrows
        if (i.eq.ptr_disp) write(stderr,'("displacements")')
        if (i.eq.ptr_stress) write(stderr,'("pre-stresses")')
        if (i.eq.ptr_damp) write(stderr,'("damping")')
        if (i.eq.ptr_smooth) write(stderr,'("smoothing")')
        write(stderr,'(1PE12.4)') b(i,1)
    enddo
    write(stderr,*)
endif
9001 format('(1P',I3,'E12.4)')

return
end

!--------------------------------------------------------------------------------------------------!

subroutine lsqr_load_disp(A,b,nrows,ncols)
use command_line, only : disp_comp, rake_file, verbosity
use arrays, only : ndisplacements, nfaults, displacements, disp_gfs
use io
implicit none
! I/O variables
integer :: nrows, ncols
double precision :: A(nrows,ncols), b(nrows,1)
! Local variables
integer :: i, j, n

if (verbosity.ge.4) then
    write(stderr,'(A)') 'Starting lsqr_load_disp()'
endif

do i = 1,ndisplacements
    do j = 1,nfaults
        ! For selected components of displacements, load strike-slip source GFs into model matrix
        if (disp_comp.eq.'a') then
            A(i,j)                    = disp_gfs(i,j,1) ! ssx
            A(i+ndisplacements,j)     = disp_gfs(i,j,2) ! ssy
            A(i+2*ndisplacements,j)   = disp_gfs(i,j,3) ! ssz
        elseif (disp_comp.eq.'h') then
            A(i,j)                    = disp_gfs(i,j,1) ! ssx
            A(i+ndisplacements,j)     = disp_gfs(i,j,2) ! ssy
        elseif (disp_comp.eq.'v') then
            A(i,j)                    = disp_gfs(i,j,3) ! ssz
        endif

        ! Choose columns to load selected fault slip components
        if (rake_file.eq.'none') then
            n = nfaults
        else
            exit
        endif

        ! For selected components of displacements, load dip-slip source GFs into model matrix
        if (disp_comp.eq.'a') then
            A(i,j+n)                  = disp_gfs(i,j,4) ! dsx
            A(i+ndisplacements,j+n)   = disp_gfs(i,j,5) ! dsy
            A(i+2*ndisplacements,j+n) = disp_gfs(i,j,6) ! dsz
        elseif (disp_comp.eq.'h') then
            A(i,j+n)                  = disp_gfs(i,j,4) ! dsx
            A(i+ndisplacements,j+n)   = disp_gfs(i,j,5) ! dsy
        elseif (disp_comp.eq.'v') then
            A(i,j+n)                  = disp_gfs(i,j,6) ! dsz
        endif
    enddo

    ! Load selected components of observed displacements into constraint vector
    if (disp_comp.eq.'a') then
        b(i,1)                  = displacements(i,4) ! x
        b(i+ndisplacements,1)   = displacements(i,5) ! y
        b(i+2*ndisplacements,1) = displacements(i,6) ! z
    elseif (disp_comp.eq.'h') then
        b(i,1)                  = displacements(i,4) ! x
        b(i+ndisplacements,1)   = displacements(i,5) ! y
    elseif (disp_comp.eq.'v') then
        b(i,1)                  = displacements(i,6) ! z
    endif
enddo

if (verbosity.ge.4) then
    write(stderr,'(A)') 'Finished lsqr_load_disp()'
    write(stderr,*)
endif

return
end

!--------------------------------------------------------------------------------------------------!

subroutine lsqr_load_stress(A,b,nrows,ncols,ptr_stress)
use command_line
use arrays, only : nfaults, prestresses, stress_gfs
use io
implicit none
! I/O variables
integer :: nrows, ncols, ptr_stress
double precision :: A(nrows,ncols), b(nrows,1)
! Local variables
integer :: i, j

if (verbosity.ge.4) then
    write(stderr,'(A)') 'Finished lsqr_load_stress()'
endif

do i = 1,nfaults
    do j = 1,nfaults
        if (rake_file.eq.'none') then
            A(ptr_stress+i-1,j)                 = stress_gfs(i,j,1) ! ss->ss
            A(ptr_stress+nfaults+i-1,j)         = stress_gfs(i,j,2) ! ss->ds
            A(ptr_stress+i-1,j+nfaults)         = stress_gfs(i,j,3) ! ds->ss
            A(ptr_stress+nfaults+i-1,j+nfaults) = stress_gfs(i,j,4) ! ds->ds
        else
            A(ptr_stress+i-1,j)                 = stress_gfs(i,j,1)
            A(ptr_stress+nfaults+i-1,j)         = stress_gfs(i,j,2)
        endif
    enddo

    ! Load pre-stresses into constraint vector
    b(ptr_stress+i-1,1)         = prestresses(i,1) ! ss
    b(ptr_stress+nfaults+i-1,1) = prestresses(i,2) ! ds
enddo

if (verbosity.ge.4) then
    write(stderr,'(A)') 'Finished lsqr_load_stress()'
    write(stderr,*)
endif

return
end

!--------------------------------------------------------------------------------------------------!

subroutine lsqr_load_damp(A,b,n,m,damp,ptr)
!----
! Add equations to model matrix to minimize the length of the solution to Ax = b
! with weighting factor damp.
!----
use command_line, only : verbosity
use io
implicit none
! I/O variables
integer :: n, m, ptr
double precision :: A(n,m), b(n,1), damp
! Local variables
integer :: i, j

if (verbosity.ge.2) then
    write(stderr,'("Starting damping_linear_ls()")')
endif

do i = 0,m-1
    b(ptr+i,1) = 0.0d0
    do j = 0,m-1
        if (i.eq.j) then
            A(ptr+i,j+1) = damp*damp*1.0d0
        else
            A(ptr+i,j+1) = 0.0d0
        endif
    enddo
enddo

if (verbosity.ge.2) then
    write(stderr,'("Finished damping_linear_ls()")')
endif

return
end

!--------------------------------------------------------------------------------------------------!

subroutine lsqr_load_smooth(A,b,nrows,ncols,smoothing_constant,ptr_smooth)
!----
! Add equations to model matrix to minimize the curvature (second derivative)
! of the solution to Ax = b with weighting factor smooth.
!----
use command_line, only : rake_file
use arrays, only : nfaults, nsmooth, smoothing_pointers, smoothing_neighbors
implicit none
! I/O variables
integer :: nrows, ncols, ptr_smooth
double precision :: A(nrows,ncols), b(nrows,1), smoothing_constant
! Local variables
integer :: i, j, ifault, nneighbor, ineighbor, neighbor

! Read smoothing element linking file and load into A
do i = 1,nsmooth
    b(ptr_smooth+i-1,1) = 0.0d0
    ifault = smoothing_pointers(i,1)
    nneighbor = smoothing_pointers(i,2)
    ineighbor = smoothing_pointers(i,3)
    ! Fault to be smoothed gets weight of nneighbor
    A(ptr_smooth+i-1,ifault) = dble(nneighbor)*smoothing_constant*smoothing_constant
    if (rake_file.eq.'none') then
        A(ptr_smooth+i-1,ifault+nfaults) = dble(nneighbor)*smoothing_constant*smoothing_constant
    endif
    ! Each neighboring fault gets weight of -1
    do j = 0,nneighbor-1
        neighbor = smoothing_neighbors(ineighbor+j)
        A(ptr_smooth+i-1,neighbor) = -1.0d0*smoothing_constant*smoothing_constant
        if (rake_file.eq.'none') then
            A(ptr_smooth+i-1,neighbor+nfaults) = -1.0d0*smoothing_constant*smoothing_constant
        endif
    enddo
enddo

return
end

!----------------------------------------------------------------------C
!
!      SUBROUTINE getdamp(soln,gf,obs,nobs,OBSMAX,nflt,FLTMAX,damp,
!     1                   smooth,smoof,fact,getdmp,sgf,stscon,nsts,cmpnt)
!      IMPLICIT none
!      INTEGER getdmp
!      REAL*8 stscon
!      INTEGER nsts,cmpnt
!      CHARACTER*80 smoof
!      REAL*8 damp,smooth
!      INTEGER OBSMAX,FLTMAX
!      REAL*8 obs(OBSMAX,6),prests(FLTMAX,6)
!      REAL*8 gf(OBSMAX,FLTMAX,6),sgf(FLTMAX,FLTMAX,4)
!      REAL*8 soln(FLTMAX,2),misfit
!      REAL*8 p(3)
!      INTEGER i,j,nobs,nflt
!      REAL*8 slip,slipmx,ddamp,misfit0,fact
!      INTEGER vrb
!      COMMON /VERBOSE/ vrb
!! Compute misfit with damp = 0
!      damp = 0.0d0
!      call lstsqr(soln,gf,obs,nobs,OBSMAX,nflt,FLTMAX,damp,
!     1            smooth,smoof,sgf,stscon,nsts,prests,cmpnt)
!      misfit0 = 0.0d0
!      do 108 i = 1,nobs
!          p(1) = 0.0d0
!          p(2) = 0.0d0
!          p(3) = 0.0d0
!          do 107 j = 1,nflt
!              p(1) = p(1) + gf(i,j,1)*soln(j,1)+gf(i,j,4)*soln(j,2)
!              p(2) = p(2) + gf(i,j,2)*soln(j,1)+gf(i,j,5)*soln(j,2)
!              p(3) = p(3) + gf(i,j,3)*soln(j,1)+gf(i,j,6)*soln(j,2)
!  107     continue
!          misfit0 = misfit0 + (p(1)-obs(i,4))*(p(1)-obs(i,4))
!     1                    + (p(2)-obs(i,5))*(p(2)-obs(i,5))
!     2                    + (p(3)-obs(i,6))*(p(3)-obs(i,6))
!  108 continue
!      if (vrb.eq.1) then
!          write(0,1001) misfit0
!      endif
!! Search for damping parameter
!      ddamp = 1.0d0
!  101 damp = damp + ddamp
!          call lstsqr(soln,gf,obs,nobs,OBSMAX,nflt,FLTMAX,damp,
!     1                smooth,smoof,sgf,stscon,nsts,prests,cmpnt)
!! Compute RMS obs-pre misfit
!          misfit = 0.0d0
!          do 102 i = 1,nobs
!              p(1) = 0.0d0
!              p(2) = 0.0d0
!              p(3) = 0.0d0
!              do 103 j = 1,nflt
!                  p(1) = p(1) + gf(i,j,1)*soln(j,1)+gf(i,j,4)*soln(j,2)
!                  p(2) = p(2) + gf(i,j,2)*soln(j,1)+gf(i,j,5)*soln(j,2)
!                  p(3) = p(3) + gf(i,j,3)*soln(j,1)+gf(i,j,6)*soln(j,2)
!  103         continue
!              misfit = misfit + (p(1)-obs(i,4))*(p(1)-obs(i,4))
!     1                        + (p(2)-obs(i,5))*(p(2)-obs(i,5))
!     2                        + (p(3)-obs(i,6))*(p(3)-obs(i,6))
!  102     continue
!          if (vrb.eq.1) then
!              write(0,1002) damp,misfit
!          endif
!! Get maximum slip in model
!          slipmx = 0.0d0
!          do 105 i = 1,nflt
!              slip = dsqrt(soln(i,1)*soln(i,1)+soln(i,2)*soln(i,2))
!              if (slip.gt.slipmx) slipmx = slip
!  105     continue
!          print *,getdmp,misfit,misfit0,slipmx,fact
!      if (getdmp.eq.1.and.misfit.le.misfit0*fact) then
!          goto 101
!      elseif (getdmp.eq.2.and.slipmx.gt.fact) then
!          goto 101
!      elseif (ddamp.gt.1.0d-4) then
!          damp = damp - ddamp
!          ddamp = ddamp*1.0d-1
!          damp = damp - ddamp
!          goto 101
!      else
!          continue
!      endif
!      RETURN
! 1001 format('    Misfit with no damping is ',1PE14.6)
! 1002 format('    Misfit with damping = ',1PE10.3,' is ',1PE14.6)
!      END

!--------------------------------------------------------------------------------------------------!

subroutine lsqr_solve(x,A,b,nrows,ncols)
!----
! Solve generalized least squares problem Ax = b for x.
! Arrays A and b have dimensions: A(nrows,ncols), b(nrows,1)
!----
use command_line, only : verbosity
use io
implicit none
! I/O variables
integer :: nrows, ncols
double precision :: A(nrows,ncols), b(nrows,1), x(nrows,1)
! Local variables
integer :: i, m, n, lda, ldb
integer :: nrhs
integer :: lwork, info
character(len=1) :: trans
double precision :: work(100000)
double precision :: btmp(nrows,1)

if (verbosity.ge.2) then
    write(stderr,'(A)') 'Starting lsqr_solve()'
endif

nrhs = 1
trans = 'N'      ! A has form (nrows x ncols), i.e. not transposed
m = nrows        ! Number of rows in matrix A
n = ncols        ! Number of columns in matrix A
lda = m          ! Leading dimension of A;  lda >= max(1,m)
ldb = max(m,n)   ! Leading dimension of b;  ldb >= max(1,m,n)

! Compute optimal workspace for least squares problem
lwork = -1
call dgels(trans,m,n,nrhs,A,lda,b,ldb,work,lwork,info)
lwork = int(work(1))

! Copy displacement vector, b to btmp because it is replaced in dgels
do i = 1,nrows
    btmp(i,1) = b(i,1)
enddo

! Solve least squares problem for x
call dgels(trans,m,n,nrhs,a,lda,btmp,ldb,work,lwork,info)

! Put solution into x
do i = 1,ncols
    x(i,1) = btmp(i,1)
enddo

if (verbosity.ge.3) then
    write(stderr,'("x =")')
    do i = 1,ncols
        write(stderr,'(1PE12.4)') x(i,1)
    enddo
endif
if (verbosity.ge.2) then
    write(stderr,'(A)') 'Finished lsqr_solve()'
    write(stderr,*)
endif

return
end
