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
use arrays, only : fault_slip
use io, only : stderr
implicit none
! Local variables
double precision, dimension(:,:), allocatable :: A, b, x, Atmp
integer :: nrows, ncols, ncols2, ptr_disp, ptr_stress, ptr_damp, ptr_smooth

if (verbosity.ge.2) then
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
A = 0.0d0
b = 0.0d0

! Load model matrix, A, and constraint vector, b
ncols2 = ncols
call lsqr_load_arrays(A,b,nrows,ncols2,ptr_disp,ptr_stress,ptr_damp,ptr_smooth)

! Resize A if necessary
if (ncols2.lt.ncols) then
    if (.not.allocated(Atmp)) then
        allocate(Atmp(nrows,ncols2))
    endif
    Atmp(1:nrows,1:ncols2) = A(1:nrows,1:ncols2)
    deallocate(A)
    allocate(A(nrows,ncols2))
    A = Atmp
    deallocate(Atmp)
    ncols = ncols2
endif
!! NEED TO MATCH UP OLD FAULT LIST WITH (POTENTIALLY RE-ORDERED) NEW FAULT LIST!!!!!

! Solve generalized least squares problem
if (.not.allocated(x)) then
    allocate(x(ncols,1))
endif
call lsqr_solve(x,A,b,nrows,ncols)
! call lsqr_nnls_solve(x,A,b,nrows,ncols)

if (.not.allocated(fault_slip)) then
    allocate(fault_slip(ncols))
endif
fault_slip = x(:,1)

if (verbosity.ge.2) then
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
    nrows = nrows + len_trim(disp_comp)*ndisplacements
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

! Add rows for damping (one row per model DOF)
! Damp each slip component separately
if (damping_constant.gt.0.0d0) then
    ptr_damp = nrows + 1
    nrows = nrows + ncols
endif

! Add rows for smoothing
! Smooth each slip component separately
if (smoothing_constant.gt.0.0d0) then
    ptr_smooth = nrows + 1
    if (rake_file.eq.'none') then
        nrows = nrows + 2*nsmooth
    else
        nrows = nrows + nsmooth
    endif
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

! Apply constraints
if (slip_constraint_file.ne.'none') then
    call lsqr_load_slip_constraint(A,b,nrows,ncols,ptr_stress,ptr_smooth)
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

if (verbosity.ge.2) then
    write(stderr,'(A)') 'Starting lsqr_load_disp()'
endif

do i = 1,ndisplacements
    do j = 1,nfaults
        ! For selected components of displacements, load strike-slip (or constant rake)
        ! source GFs into model matrix
        if (disp_comp.eq.'123') then
            A(i,j)                    = disp_gfs(i,j,1) ! ssx
            A(i+ndisplacements,j)     = disp_gfs(i,j,2) ! ssy
            A(i+2*ndisplacements,j)   = disp_gfs(i,j,3) ! ssz
        elseif (disp_comp.eq.'12') then
            A(i,j)                    = disp_gfs(i,j,1) ! ssx
            A(i+ndisplacements,j)     = disp_gfs(i,j,2) ! ssy
        elseif (disp_comp.eq.'13') then
            A(i,j)                    = disp_gfs(i,j,1) ! ssx
            A(i+ndisplacements,j)     = disp_gfs(i,j,3) ! ssz
        elseif (disp_comp.eq.'23') then
            A(i,j)                    = disp_gfs(i,j,2) ! ssy
            A(i+ndisplacements,j)     = disp_gfs(i,j,3) ! ssz
        elseif (disp_comp.eq.'1') then
            A(i,j)                    = disp_gfs(i,j,1) ! ssx
        elseif (disp_comp.eq.'2') then
            A(i,j)                    = disp_gfs(i,j,2) ! ssy
        elseif (disp_comp.eq.'3') then
            A(i,j)                    = disp_gfs(i,j,3) ! ssz
        endif

        ! Choose columns to load selected fault slip components
        if (rake_file.eq.'none') then
            n = nfaults
            ! For selected components of displacements, load dip-slip source GFs into model matrix
            if (disp_comp.eq.'123') then
                A(i,j+n)                  = disp_gfs(i,j,4) ! dsx
                A(i+ndisplacements,j+n)   = disp_gfs(i,j,5) ! dsy
                A(i+2*ndisplacements,j+n) = disp_gfs(i,j,6) ! dsz
            elseif (disp_comp.eq.'12') then
                A(i,j+n)                  = disp_gfs(i,j,4) ! dsx
                A(i+ndisplacements,j+n)   = disp_gfs(i,j,5) ! dsy
            elseif (disp_comp.eq.'13') then
                A(i,j+n)                  = disp_gfs(i,j,4) ! dsx
                A(i+ndisplacements,j+n)   = disp_gfs(i,j,6) ! dsz
            elseif (disp_comp.eq.'23') then
                A(i,j+n)                  = disp_gfs(i,j,5) ! dsy
                A(i+ndisplacements,j+n)   = disp_gfs(i,j,6) ! dsz
            elseif (disp_comp.eq.'1') then
                A(i,j+n)                  = disp_gfs(i,j,4) ! dsx
            elseif (disp_comp.eq.'2') then
                A(i,j+n)                  = disp_gfs(i,j,5) ! dsy
            elseif (disp_comp.eq.'3') then
                A(i,j+n)                  = disp_gfs(i,j,6) ! dsz
            endif
        endif
    enddo

    ! Load selected components of observed displacements into constraint vector
    if (disp_comp.eq.'123') then
        b(i,1)                  = displacements(i,4) ! x
        b(i+ndisplacements,1)   = displacements(i,5) ! y
        b(i+2*ndisplacements,1) = displacements(i,6) ! z
    elseif (disp_comp.eq.'12') then
        b(i,1)                  = displacements(i,4) ! x
        b(i+ndisplacements,1)   = displacements(i,5) ! y
    elseif (disp_comp.eq.'13') then
        b(i,1)                  = displacements(i,4) ! x
        b(i+ndisplacements,1)   = displacements(i,6) ! z
    elseif (disp_comp.eq.'23') then
        b(i,1)                  = displacements(i,5) ! y
        b(i+ndisplacements,1)   = displacements(i,6) ! z
    elseif (disp_comp.eq.'1') then
        b(i,1)                  = displacements(i,4) ! x
    elseif (disp_comp.eq.'2') then
        b(i,1)                  = displacements(i,5) ! y
    elseif (disp_comp.eq.'3') then
        b(i,1)                  = displacements(i,6) ! z
    endif
enddo

if (verbosity.ge.2) then
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

if (verbosity.ge.2) then
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

if (verbosity.ge.2) then
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
    write(stderr,'(A)') "Starting lsqr_load_damp()"
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
    write(stderr,'(A)') "Finished lsqr_load_damp()"
    write(stderr,*)
endif

return
end

!--------------------------------------------------------------------------------------------------!

subroutine lsqr_load_smooth(A,b,nrows,ncols,smoothing_constant,ptr_smooth)
!----
! Add equations to model matrix to minimize the curvature (second derivative)
! of the solution to Ax = b with weighting factor smooth.
!----
use command_line, only : verbosity, rake_file
use arrays, only : nfaults, nsmooth, smoothing_pointers, smoothing_neighbors
use io, only : stderr
implicit none
! I/O variables
integer :: nrows, ncols, ptr_smooth
double precision :: A(nrows,ncols), b(nrows,1), smoothing_constant
! Local variables
integer :: i, j, ifault, nneighbor, ineighbor, neighbor

if (verbosity.ge.2) then
    write(stderr,'(A)') "Starting lsqr_load_smooth()"
endif

! Read smoothing element linking file and load into A
do i = 1,nsmooth
    b(ptr_smooth+i-1,1) = 0.0d0
    ifault = smoothing_pointers(i,1)
    nneighbor = smoothing_pointers(i,2)
    ineighbor = smoothing_pointers(i,3)
    ! Fault to be smoothed gets weight of nneighbor
    A(ptr_smooth+i-1,ifault) = dble(nneighbor)*smoothing_constant*smoothing_constant
    if (rake_file.eq.'none') then
        b(ptr_smooth+nsmooth+i-1,1) = 0.0d0
        A(ptr_smooth+nsmooth+i-1,ifault+nfaults) = &
                               dble(nneighbor)*smoothing_constant*smoothing_constant
    endif
    ! Each neighboring fault gets weight of -1
    do j = 0,nneighbor-1
        neighbor = smoothing_neighbors(ineighbor+j)
        A(ptr_smooth+i-1,neighbor) = -1.0d0*smoothing_constant*smoothing_constant
        if (rake_file.eq.'none') then
            A(ptr_smooth+nsmooth+i-1,neighbor+nfaults) = &
                                     -1.0d0*smoothing_constant*smoothing_constant
        endif
    enddo
enddo

if (verbosity.ge.2) then
    write(stderr,'(A)') "Finished lsqr_load_smooth()"
    write(stderr,*)
endif

return
end

!--------------------------------------------------------------------------------------------------!

subroutine lsqr_load_slip_constraint(A,b,nrows,ncols,ptr_stress,ptr_smooth)
use command_line, only : displacement_file, prestress_file, disp_comp, rake_file, &
                         smoothing_file, verbosity
use arrays, only : nfaults, slip_constraints, is_this_fault_constrained, &
                   ndisplacements, nsmooth
use io, only : stderr
implicit none
! I/O variables
integer :: nrows, ncols, ptr_stress, ptr_smooth
double precision :: A(nrows,ncols), b(nrows,1)
! Local variables
integer :: i, j, n, dn

if (verbosity.ge.2) then
    write(stderr,'(A)') "Starting lsqr_load_slip_constraint()"
endif

! Move contribution from constrained fault segments to observation side of equation (do not
! renumber yet.....will do that later in subroutine)
do i = 1,nfaults
    ! Is the strike-slip (or constant rake) component of this fault constrained?
    if (is_this_fault_constrained(i,1).eq.1) then
        ! Move displacements corresponding to this slip component to b vector
        if (displacement_file.ne.'none') then
            do n = 0,len_trim(disp_comp)-1
                dn = n*ndisplacements
                do j = 1,ndisplacements
                    b(j+dn,1) = b(j+dn,1) - A(j+dn,i)*slip_constraints(i,1)
                enddo
            enddo
        endif
        ! Move stresses corresponding to this slip component to b vector
        if (prestress_file.ne.'none') then
            do j = 1,nfaults
                b(ptr_stress+j-1,1) = b(ptr_stress+j-1,1) &
                                        - A(ptr_stress+j-1,i)*slip_constraints(i,1)
            enddo
        endif
        ! Move smoothing corresponding to this slip component to b vector
        if (smoothing_file.ne.'none') then
            do j = 1,nsmooth
                b(ptr_smooth+j-1,1) = b(ptr_smooth+j-1,1) &
                                        - A(ptr_smooth+j-1,i)*slip_constraints(i,1)
            enddo
        endif
        A(1,i) = -1.0d10
    endif

    ! Is the dip-slip component of this fault constrained?
    if (is_this_fault_constrained(i,2).eq.1.and.rake_file.eq.'none') then
        ! Move displacements corresponding to this slip component to b vector
        if (displacement_file.ne.'none') then
            do n = 0,len_trim(disp_comp)-1
                dn = n*ndisplacements
                do j = 1,ndisplacements
                    b(j+dn,1) = b(j+dn,1) - A(j+dn,i+nfaults)*slip_constraints(i,2)
                enddo
            enddo
        endif
        ! Move stresses corresponding to this slip component to b vector
        if (prestress_file.ne.'none') then
            do j = 1,nfaults
                b(ptr_stress+nfaults+j-1,1) = b(ptr_stress+nfaults+j-1,1) &
                                 - A(ptr_stress+nfaults+j-1,i)*slip_constraints(i,2)
            enddo
        endif
        ! Move smoothing corresponding to this slip component to b vector
        if (smoothing_file.ne.'none') then
            do j = 1,nsmooth
                b(ptr_smooth+nsmooth+j-1,1) = b(ptr_smooth+nsmooth+j-1,1) &
                        - A(ptr_smooth+nsmooth+j-1,i+nfaults)*slip_constraints(i,2)
            enddo
        endif
        A(1,i+nfaults) = -1.0d10
    endif
enddo

! Renumber rows and columns of arrays
n = 0
do i = 1,ncols
    if (A(1,i).lt.-1.0d9) then
        if (i.lt.ncols) then
            do j = i,ncols-1
                A(:,j) = A(:,j+1)
            enddo
        endif
        n = n + 1
    endif
enddo
ncols = ncols - n

if (verbosity.ge.2) then
    write(stderr,'(A)') "Finishes lsqr_load_slip_constraint()"
    write(stderr,*)
endif

return
end

!--------------------------------------------------------------------------------------------------!

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

!--------------------------------------------------------------------------------------------------!

subroutine lsqr_nnls_solve(x,A,b,nrows,ncols)
implicit none
! I/O variables
integer :: nrows, ncols
double precision :: A(nrows,ncols), b(nrows,1), x(nrows,1)
! Local variables
integer :: mode
integer, dimension(:), allocatable :: indx
double precision :: rnorm
double precision, dimension(:), allocatable :: btmp, xtmp, w

allocate(btmp(nrows))
btmp = b(:,1)
allocate(xtmp(ncols))
allocate(w(ncols))
allocate(indx(ncols))

call nnls(A, nrows, ncols, btmp, xtmp, rnorm, w, indx, mode)
x(:,1) = xtmp

return
end
