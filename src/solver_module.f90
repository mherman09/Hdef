!--------------------------------------------------------------------------------------------------!
! Module solver
!
! Routines to set up and solve linear systems of equations Ax = b
!
! Subroutines in module:
!     - load_array            Load a sub-array into a main array
!     - load_constraints      Load x constraints into A and b
!     - solve_dgels           Implement LAPACK general least-squares solver
!     - solve_nnls            Implement non-negative least squares solver
!     - solve_dgesv           Implement LAPACK general matrix solver
!     - solve_dgssv           Implement SuperLU sparse matrix solver
!     - nnls                  Non-negative least squares solver
!     - g1                    Subroutine used by nnls (orthogonal rotation matrix)
!     - h12                   Subroutine used by nnls (Householder transformation)
!
!--------------------------------------------------------------------------------------------------!

module solver

interface load_array
    module procedure load_dp_array1
    module procedure load_dp_array2
    module procedure load_int_array2
end interface

interface load_constraints
    module procedure load_dp_array_constraints
end interface

! Setup subroutines
public :: load_array
public :: load_constraints


! Solver subroutines
public :: solve_nnls

#ifdef USE_LAPACK
public :: solve_dgels
public :: solve_dgesv
public :: solve_dsysv
#endif

#ifdef USE_SUPERLU
public :: solve_dgssv
#endif


! Helper subroutines
public :: nnls
private :: g1
private :: h12

!--------------------------------------------------------------------------------------------------!
contains
!--------------------------------------------------------------------------------------------------!

subroutine load_dp_array1(array,nrows,subarray,nrows_sub,irow,subarray_lbl,ierr)
!----
! Load a double precision sub-array into the main array (1-D)
!
! Arguments:
!     array:        Main array; receives new entries
!     nrows:        Number of rows in main array
!     subarray:     Array to be loaded into main array (overwrites contents of array)
!     nrows_sub:    Number of rows in sub-array
!     irow:         Row index to start loading sub-array
!     subarray_lbl: Name of sub-array for error messages
!     ierr:         Exit status
!----

use io, only: stderr

implicit none

! Arguments
integer :: nrows
integer :: nrows_sub
integer :: irow
double precision :: array(nrows)
double precision :: subarray(nrows_sub)
character(len=*) :: subarray_lbl
integer :: ierr

! Local variables
integer :: max_row

ierr = 0

! Check that the input values do not overflow rows
max_row = irow + nrows_sub - 1
if (max_row.gt.nrows) then
    write(stderr,*) 'load_array: loading sub-array ',trim(subarray_lbl),' overflows array rows'
    write(stderr,*) 'nrows:     ',nrows
    write(stderr,*) 'nrows_sub: ',nrows_sub
    write(stderr,*) 'irow:      ',irow
    write(stderr,*) 'max_row:   ',max_row
    ierr = 1
    return
endif

! Load the input values into array
array(irow:max_row) = subarray(1:nrows_sub)

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine load_dp_array2(array,nrows,ncols,subarray,nrows_sub,ncols_sub,irow,icol,subarray_lbl,ierr)
!----
! Load a double precision sub-array into the main array (2-D array)
!
! Arguments:
!     array:        Main array; receives new entries
!     nrows:        Number of rows in main array
!     ncols:        Number of columns in main array
!     subarray:     Array to be loaded into main array (overwrites contents of array)
!     nrows_sub:    Number of rows in sub-array
!     ncols_sub:    Number of columns in sub-array
!     irow:         Row index to start loading sub-array
!     icol:         Column index to start loading sub-array
!     subarray_lbl: Name of sub-array for error messages
!     ierr:         Exit status
!----

use io, only: stderr

implicit none

! Arguments
integer :: nrows
integer :: ncols
integer :: nrows_sub
integer :: ncols_sub
integer :: irow
integer :: icol
double precision :: array(nrows,ncols)
double precision :: subarray(nrows_sub,ncols_sub)
character(len=*) :: subarray_lbl
integer :: ierr

! Local variables
integer :: max_row, max_col

ierr = 0

! Check that the input values do not overflow rows
max_row = irow + nrows_sub - 1
if (max_row.gt.nrows) then
    write(stderr,*) 'load_array: loading sub-array ',trim(subarray_lbl),' overflows array rows'
    write(stderr,*) 'nrows:     ',nrows
    write(stderr,*) 'nrows_sub: ',nrows_sub
    write(stderr,*) 'irow:      ',irow
    write(stderr,*) 'max_row:   ',max_row
    ierr = 1
    return
endif

! Check that the input values do not overflow columns
max_col = icol + ncols_sub - 1
if (max_col.gt.ncols) then
    write(stderr,*) 'load_array: loading sub-array ',trim(subarray_lbl),' overflows array cols'
    write(stderr,*) 'ncols:     ',ncols
    write(stderr,*) 'ncols_sub: ',ncols_sub
    write(stderr,*) 'icol:      ',icol
    write(stderr,*) 'max_col:   ',max_col
    ierr = 1
    return
endif

! Load the input values into array
array(irow:max_row,icol:max_col) = subarray(1:nrows_sub,1:ncols_sub)

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine load_int_array2(array,nrows,ncols,subarray,nrows_sub,ncols_sub,irow,icol,subarray_lbl,ierr)
!----
! Load an integer sub-array into the main array
!
! Arguments:
!     array:        Main array; receives new entries
!     nrows:        Number of rows in main array
!     ncols:        Number of columns in main array
!     subarray:     Array to be loaded into main array (overwrites contents of array)
!     nrows_sub:    Number of rows in sub-array
!     ncols_sub:    Number of columns in sub-array
!     irow:         Row index to start loading sub-array
!     icol:         Column index to start loading sub-array
!     subarray_lbl: Name of sub-array for error messages
!     ierr:         Exit status
!----

use io, only: stderr

implicit none

! Arguments
integer :: nrows
integer :: ncols
integer :: nrows_sub
integer :: ncols_sub
integer :: irow
integer :: icol
integer :: array(nrows,ncols)
integer :: subarray(nrows_sub,ncols_sub)
character(len=*) :: subarray_lbl
integer :: ierr

! Local variables
integer :: max_row, max_col

ierr = 0

! Check that the input values do not overflow rows
max_row = irow + nrows_sub - 1
if (max_row.gt.nrows) then
    write(stderr,*) 'load_array: loading sub-array ',trim(subarray_lbl),' overflows array rows'
    write(stderr,*) 'nrows:     ',nrows
    write(stderr,*) 'nrows_sub: ',nrows_sub
    write(stderr,*) 'irow:      ',irow
    write(stderr,*) 'max_row:   ',max_row
    ierr = 1
    return
endif

! Check that the input values do not overflow columns
max_col = icol + ncols_sub - 1
if (max_col.gt.ncols) then
    write(stderr,*) 'load_array: loading sub-array ',trim(subarray_lbl),' overflows array cols'
    write(stderr,*) 'ncols:     ',ncols
    write(stderr,*) 'ncols_sub: ',ncols_sub
    write(stderr,*) 'icol:      ',icol
    write(stderr,*) 'max_col:   ',max_col
    ierr = 1
    return
endif

! Load the input values into array
array(irow:max_row,icol:max_col) = subarray(1:nrows_sub,1:ncols_sub)

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine load_dp_array_constraints(Ain,bin,nrows,ncols,constraints,isXConstr,ierr)
!----
! Apply x-constraints to equation Ax = b by (1) moving constrained values to b vector and
! (2) renumbering A matrix column indices.
!
! Arguments
!     A:            Input: full coefficient matrix, A
!                   Output: reduced-column coefficient matrix (overwrites A)
!     b:            Input: original observation vector, b
!                   Output: observation vector modified by constraints (overwrites b)
!     nrows:        Input: number of rows in full A matrix, original b vector
!     ncols:        Input: number of columns in full A matrix, number of rows in x vector
!                   Output: number of columns in reduced-column coefficient matrix
!     constraints:  Input: values of x-constraints (ncols-length double precision array)
!     isXConstr:    Input: ncols-length logical array, true if this x value is constrained
!     ierr:         Output: exit status
!----

use io, only: stderr, stdout, verbosity

implicit none

! Arguments
integer :: nrows
integer :: ncols
double precision :: Ain(nrows,ncols)
double precision :: bin(nrows)
double precision :: constraints(ncols)
logical :: isXConstr(ncols)
integer :: ierr

! Local variables
integer :: icol, irow, nkept
logical :: deleteCol(ncols)


ierr = 0

! Move contributions from x-constraints to b vector
deleteCol = .false.
do icol = 1,ncols
    if (isXConstr(icol)) then
        ! write(0,*) 'load_constraints: moving constraint ',icol,' to b'
        do irow = 1,nrows
            bin(irow) = bin(irow) - Ain(irow,icol)*constraints(icol)
        enddo
        deleteCol(icol) = .true.
    endif
enddo

! Renumber columns of array Ain
if (verbosity.ge.2) then
    write(stdout,*) 'load_constraints: input ncols= ',ncols
endif

! Start counters for existing columns and number of shifted columns
nkept = 0
do icol = 1,ncols
    ! If this slip value has not been fixed, then add it to the modified array
    if (.not.deleteCol(icol)) then
        nkept = nkept + 1
        Ain(:,nkept) = Ain(:,icol)
    else
        ! write(0,*) 'load_constraints: deleting column ',icol,' from A'
    endif
enddo
ncols = nkept

if (verbosity.ge.2) then
    write(stdout,*) 'load_constraints: output ncols= ',nkept
endif

return
end subroutine

!--------------------------------------------------------------------------------------------------!

#ifdef USE_LAPACK

subroutine solve_dgels(Ain,bin,xout,nrows,ncols,ierr)
!----
! Solve Ax = b for x using LAPACK general least squares driver, dgels.
!
! Arguments:
!     Ain:          Input coefficient matrix
!     bin:          Input vector
!     xout:         Output solution vector
!     nrows:        Number of rows in coefficient matrix, input vector
!     ncols:        Number of columns in coefficient matrix, number of rows in solution vector
!     ierr:         Exit status
!----

use io, only : stderr

implicit none

! Arguments
integer :: nrows
integer :: ncols
double precision :: Ain(nrows,ncols)
double precision :: bin(nrows)
double precision :: xout(ncols)
integer :: ierr

! Local variables
integer :: m, n, lda, ldb, nrhs, lwork, info
character(len=1) :: trans
double precision, allocatable :: work(:)
double precision :: alocal(nrows,ncols), blocal(nrows,1)

ierr = 0

! Use local arrays for inversion
alocal = Ain
blocal(:,1) = bin

! LAPACK driver routine dgels variables
nrhs = 1
trans = 'N'      ! A has form (nrows x ncols), i.e. not transposed
m = nrows        ! Number of rows in matrix A
n = ncols        ! Number of columns in matrix A
lda = m          ! Leading dimension of A;  lda >= max(1,m)
ldb = max(m,n)   ! Leading dimension of b;  ldb >= max(1,m,n)

! Compute optimal workspace for least squares problem
if (.not.allocated(work)) then
    allocate(work(1))
endif
lwork = -1
call dgels(trans,m,n,nrhs,alocal,lda,blocal,ldb,work,lwork,info)
if (info.ne.0) then
    write(stderr,*) 'solve_dgels: info returned ',info,' indicating error in dgels'
    ierr = 1
    return
endif

! Set workspace arrays to optimal size
lwork = int(work(1))
deallocate(work)
allocate(work(lwork))

! Solve system of equations using generalized least squares algorithm
call dgels(trans,m,n,nrhs,alocal,lda,blocal,ldb,work,lwork,info)
if (info.ne.0) then
    write(stderr,*) 'solve_dgels: info returned ',info,' indicating error in dgels'
    ierr = 1
    return
endif

! Free work array memory
deallocate(work)

! Put solution into x
xout = blocal(1:ncols,1)

return
end subroutine

#endif

!--------------------------------------------------------------------------------------------------!

subroutine solve_nnls(Ain,bin,xout,nrows,ncols,ierr)
!----
! Solve Ax = b for x using non-negative least squares. (Subroutines for implementing non-
! negative least squares are included at the bottom of this module.)
!
! Arguments:
!     Ain:          Input coefficient matrix
!     bin:          Input vector
!     xout:         Output solution vector
!     nrows:        Number of rows in coefficient matrix, input vector
!     ncols:        Number of columns in coefficient matrix, number of rows in solution vector
!     ierr:         Exit status
!----

use io, only : stderr

implicit none

! Arguments
integer :: nrows
integer :: ncols
double precision :: Ain(nrows,ncols)
double precision :: bin(nrows)
double precision :: xout(ncols)
integer :: ierr

! Local variables
integer :: m, n, mode, indx(ncols)
double precision :: alocal(nrows,ncols), blocal(nrows), rnorm, w(ncols)

ierr = 0

! Use local arrays for inversion
alocal = Ain
blocal = bin

! Subroutine nnls variables
m = nrows           ! Number of rows in matrix A
n = ncols           ! Number of columns in matrix A

! Solve system of equations using non-negative least squares algorithm
call nnls(alocal,m,n,blocal,xout,rnorm,w,indx,mode)
if (mode.ne.1) then
    write(stderr,*) 'solve_nnls: mode returned ',mode,' indicating error in nnls'
    ierr = 1
    return
endif

return
end subroutine

!--------------------------------------------------------------------------------------------------!

#ifdef USE_LAPACK

subroutine solve_dgesv(Ain,bin,xout,nrows,ierr)
!----
! Solve Ax = b for x using LAPACK general matrix solver, dgesv.
!
! Arguments:
!     Ain:          Input coefficient matrix (nrows x nrows)
!     bin:          Input vector (nrows)
!     xout:         Output solution vector (nrows)
!     nrows:        Number of rows in coefficient matrix, input vector, output vector
!     ierr:         Exit status
!----

use io, only: stderr

implicit none
!
! Arguments
integer :: nrows
double precision :: Ain(nrows,nrows)
double precision :: bin(nrows)
double precision :: xout(nrows)
integer :: ierr

! Local variables
integer :: n, nrhs, lda, ipiv(nrows), ldb, info
double precision :: alocal(nrows,nrows), blocal(nrows,1)


ierr = 0

! Use local arrays for inversion
alocal = Ain
blocal(:,1) = bin

! LAPACK driver routine dgesv variables
n = nrows           ! Number of linear equations
nrhs = 1            ! Number of right hand sides
lda = n             ! Leading dimension of A;  lda >= max(1,n)
ldb = n             ! Leading dimension of b;  ldb >= max(1,n)

! Solve system of equations using dgesv algorithm
call dgesv(n,nrhs,alocal,lda,ipiv,blocal,ldb,info)
if (info.ne.0) then
    write(stderr,*) 'solve_dgesv: info returned ',info,' indicating error in dgesv'
    if (info.lt.0) then
        write(stderr,*) 'solve_dgesv: ',info,'th argument had an illegal value'
    else
        write(stderr,*) 'solve_dgesv: factor U(',info,',',info,')th is singular'
    endif
    ierr = 1
    return
endif

! Put solution into x
xout = blocal(1:nrows,1)

return
end subroutine

#endif

!--------------------------------------------------------------------------------------------------!

#ifdef USE_LAPACK

subroutine solve_dgesv_nrhs(Ain,Bin,Xout,nrows,nrhs,ierr)
!----
! Solve AX = B for X using LAPACK general matrix solver, dgesv.
!
! Arguments:
!     Ain:          Input coefficient matrix (nrows x nrows)
!     Bin:          Input matrix (nrows x nrhs)
!     Xout:         Output solution matrix (nrows x nrhs)
!     nrows:        Number of rows in coefficient matrix, observation matrix
!     nrhs:         Number of columns in observation matrix
!     ierr:         Exit status
!----

use io, only: stderr

implicit none
!
! Arguments
integer :: nrows
integer :: ncols
double precision :: Ain(nrows,nrows)
double precision :: Bin(nrows,nrhs)
double precision :: Xout(nrows,nrhs)
integer :: ierr

! Local variables
integer :: n, nrhs, lda, ipiv(nrows), ldb, info
double precision :: alocal(nrows,nrows), blocal(nrows,nrhs)


ierr = 0

! Use local arrays for inversion
alocal = Ain
blocal = Bin

! LAPACK driver routine dgesv variables
n = nrows           ! Number of linear equations
lda = n             ! Leading dimension of A;  lda >= max(1,n)
ldb = n             ! Leading dimension of b;  ldb >= max(1,n)

! Solve system of equations using dgesv algorithm
call dgesv(n,nrhs,alocal,lda,ipiv,blocal,ldb,info)
if (info.ne.0) then
    write(stderr,*) 'solve_dgesv_nrhs: info returned ',info,' indicating error in dgesv'
    if (info.lt.0) then
        write(stderr,*) 'solve_dgesv_nrhs: ',info,'th argument had an illegal value'
    else
        write(stderr,*) 'solve_dgesv_nrhs: factor U(',info,',',info,')th is singular'
    endif
    ierr = 1
    return
endif

! Put solution into x
xout = blocal(1:nrows,1:nrhs)

return
end subroutine

#endif

!--------------------------------------------------------------------------------------------------!

#ifdef USE_LAPACK

subroutine solve_dsysv(Ain,bin,xout,nrows,ierr)
!----
! Solve Ax = b for x using LAPACK symmetric matrix solver, dsysv.
!
! Arguments:
!     Ain:          Input coefficient matrix (nrows x nrows, symmetric)
!     bin:          Input vector (nrows)
!     xout:         Output solution vector (nrows)
!     nrows:        Number of rows in coefficient matrix, input vector, output vector
!     ierr:         Exit status
!----

use io, only: stderr

implicit none

! Arguments
integer :: nrows
double precision :: Ain(nrows,nrows)
double precision :: bin(nrows)
double precision :: xout(nrows)
integer :: ierr

! Local variables
integer :: n, nrhs, lda, ipiv(nrows), ldb, lwork, info
double precision :: alocal(nrows,nrows), blocal(nrows,1)
double precision, allocatable :: work(:)


ierr = 0

! Use local arrays for inversion
alocal = Ain
blocal(:,1) = bin

! LAPACK driver routine dgesv variables
n = nrows           ! Number of linear equations
nrhs = 1            ! Number of right hand sides
lda = n             ! Leading dimension of A;  lda >= max(1,n)
ldb = n             ! Leading dimension of b;  ldb >= max(1,n)

! Get optimal size of lwork array
allocate(work(1))
lwork = -1
call dsysv('Lower',n,nrhs,alocal,lda,ipiv,blocal,ldb,work,lwork,info)
if (info.ne.0) then
    write(stderr,*) 'solve_dsysv: info returned ',info,' indicating error in dgesv'
    if (info.lt.0) then
        write(stderr,*) 'solve_dsysv: ',info,'th argument had an illegal value'
    else
        write(stderr,*) 'solve_dsysv: factor U(',info,',',info,')th is singular'
    endif
    ierr = 1
    return
endif

! Resize lwork array and solve equation so that blocal = inv(cov_mat)*ddisp
lwork = int(work(1))
deallocate(work)
allocate(work(lwork))
call dsysv('Lower',n,nrhs,alocal,lda,ipiv,blocal,ldb,work,lwork,info)
if (info.ne.0) then
    write(stderr,*) 'solve_dsysv: info returned ',info,' indicating error in dgesv'
    if (info.lt.0) then
        write(stderr,*) 'solve_dsysv: ',info,'th argument had an illegal value'
    else
        write(stderr,*) 'solve_dsysv: factor U(',info,',',info,')th is singular'
    endif
    ierr = 1
    return
endif

! Put solution into x
xout = blocal(1:nrows,1)

deallocate(work)

return
end subroutine

#endif

!--------------------------------------------------------------------------------------------------!

#ifdef USE_LAPACK

subroutine solve_dsysv_nrhs(Ain,Bin,Xout,nrows,nrhs,ierr)
!----
! Solve AX = B for X using LAPACK symmetric matrix solver, dsysv.
!
! Arguments:
!     Ain:          Input coefficient matrix (nrows x nrows, symmetric)
!     Bin:          Input vector (nrows)
!     Xout:         Output solution vector (nrows)
!     nrows:        Number of rows in coefficient matrix, observation matrix, output matrix
!     nrhs:         Number of columns in observation matrix
!     ierr:         Exit status
!----

use io, only: stderr

implicit none

! Arguments
integer :: nrows
integer :: nrhs
double precision :: Ain(nrows,nrows)
double precision :: Bin(nrows,nrhs)
double precision :: Xout(nrows,nrhs)
integer :: ierr

! Local variables
integer :: n, lda, ipiv(nrows), ldb, lwork, info
double precision :: alocal(nrows,nrows), blocal(nrows,nrhs)
double precision, allocatable :: work(:)


ierr = 0

! Use local arrays for inversion
alocal = Ain
blocal = Bin

! LAPACK driver routine dgesv variables
n = nrows           ! Number of linear equations
lda = n             ! Leading dimension of A;  lda >= max(1,n)
ldb = n             ! Leading dimension of b;  ldb >= max(1,n)

! Get optimal size of lwork array
allocate(work(1))
lwork = -1
call dsysv('Lower',n,nrhs,alocal,lda,ipiv,blocal,ldb,work,lwork,info)
if (info.ne.0) then
    write(stderr,*) 'solve_dsysv_nrhs: info returned ',info,' indicating error in dgesv'
    if (info.lt.0) then
        write(stderr,*) 'solve_dsysv_nrhs: ',info,'th argument had an illegal value'
    else
        write(stderr,*) 'solve_dsysv_nrhs: factor U(',info,',',info,')th is singular'
    endif
    ierr = 1
    return
endif

! Resize lwork array and solve equation so that blocal = inv(cov_mat)*ddisp
lwork = int(work(1))
deallocate(work)
allocate(work(lwork))
call dsysv('Lower',n,nrhs,alocal,lda,ipiv,blocal,ldb,work,lwork,info)
if (info.ne.0) then
    write(stderr,*) 'solve_dsysv_nrhs: info returned ',info,' indicating error in dgesv'
    if (info.lt.0) then
        write(stderr,*) 'solve_dsysv_nrhs: ',info,'th argument had an illegal value'
    else
        write(stderr,*) 'solve_dsysv_nrhs: factor U(',info,',',info,')th is singular'
    endif
    ierr = 1
    return
endif

! Put solution into x
xout = blocal(1:nrows,1:nrhs)

deallocate(work)

return
end subroutine

#endif

!--------------------------------------------------------------------------------------------------!

#ifdef USE_SUPERLU

subroutine solve_dgssv(Ain,bin,xout,nrows,ncols,ierr)
!----
! Solve Ax = b for x using SuperLU general sparse matrix solver, dgssv.
!
! Arguments:
!     Ain:          Input coefficient matrix (must be square)
!     bin:          Input vector
!     xout:         Output solution vector
!     nrows:        Number of rows in coefficient matrix, input vector
!     ncols:        Number of columns in coefficient matrix, number of rows in solution vector
!     ierr:         Exit status
!----

use io, only : stderr

implicit none

! Arguments
integer :: nrows
integer :: ncols
double precision :: Ain(nrows,ncols)
double precision :: bin(nrows)
double precision :: xout(ncols)
integer :: ierr

! ! Local variables
integer :: i, j, info, nZero, nNonZero, iopt, nrhs, ldb
integer(8) :: factors
double precision :: blocal(nrows,1)
double precision, allocatable :: matrix_value(:)
integer, allocatable :: row_index(:), column_ptr(:)
double precision, parameter :: min_val = 1.0d-20

ierr = 0

if (nrows.ne.ncols) then
    write(stderr,*) 'solve_dgssv: nrows not equal to ncols'
    ierr = 1
    return
endif

! Save b to a local array
blocal(:,1) = bin

! ! Convert input matrix and vector into reduced column format for dgssv
if (.not.allocated(column_ptr)) then
    allocate(column_ptr(ncols+1))
endif

! Count the number of zero and non-zero elements of the array while loading the column_ptr array
nZero = 0
nNonZero = 0
do j = 1,ncols
    column_ptr(j) = nNonZero + 1 ! Counting starts at one like a normal person
    do i = 1,nrows
        if (dabs(Ain(i,j)).gt.min_val) then
            nNonZero = nNonZero + 1
        else
            nZero = nZero + 1
        endif
    enddo
enddo
if (dble(nZero)/dble(nNonZero).lt.0.1d0) then
    write(stderr,*) 'solve_dgssv: <10% of matrix entries are zero, consider using dense solver'
endif
column_ptr(ncols+1) = nNonZero + 1

! Allocate memory for values of matrix and row index
if (.not.allocated(matrix_value)) then
    allocate(matrix_value(nNonZero))
endif
if (.not.allocated(row_index)) then
    allocate(row_index(nNonZero))
endif

! Load matrix value and row_index arrays
nNonZero = 0
do j = 1,ncols
    do i = 1,nrows
        if (dabs(Ain(i,j)).gt.min_val) then
            nNonZero = nNonZero + 1 ! Counting starts at one
            matrix_value(nNonZero) = Ain(i,j)
            row_index(nNonZero) = i
        endif
    enddo
enddo

! SuperLU routine dgssv variables
nrhs = 1
ldb = nrows

! Factorize matrix
iopt = 1
call c_fortran_dgssv(iopt,nrows,nNonZero,nrhs,matrix_value,row_index,column_ptr,blocal,ldb,factors,&
                     info)

! Solve the system with factors
iopt = 2
call c_fortran_dgssv(iopt,nrows,nNonZero,nrhs,matrix_value,row_index,column_ptr,blocal,ldb,factors,&
                     info)

! Free allocated storage
iopt = 3
call c_fortran_dgssv(iopt,nrows,nNonZero,nrhs,matrix_value,row_index,column_ptr,blocal,ldb,factors,&
                     info)

if (info.ne.0) then
    write(stderr,*) 'solve_dgssv: info returned ',info,' indicating error in dgssv'
    ierr = 1
    return
endif

! Put solution into x
xout = blocal(1:ncols,1)

! Free memory
if (allocated(column_ptr)) then
    deallocate(column_ptr)
endif
if (allocated(matrix_value)) then
    deallocate(matrix_value)
endif
if (allocated(row_index)) then
    deallocate(row_index)
endif

return
end subroutine

#endif

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
! NON-NEGATIVE LEAST SQUARES SUBROUTINES - EDIT AT YOUR OWN RISK                                   !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!
!     SUBROUTINE nnls(a, m, n, b, x, rnorm, w, indx, mode)
!
!  Algorithm NNLS: NONNEGATIVE LEAST SQUARES
!
!  The original version of this code was developed by
!  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
!  1973 JUN 15, and published in the book
!  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
!  Revised FEB 1995 to accompany reprinting of the book by SIAM.
!
!  This translation into Fortran 90 by Alan Miller, February 1997
!  Latest revision - 15 April 1997

!  N.B. The following call arguments have been removed:
!       mda, zz
!
!  GIVEN AN M BY N MATRIX, A, AND AN M-VECTOR, B,  COMPUTE AN
!  N-VECTOR, X, THAT SOLVES THE LEAST SQUARES PROBLEM
!
!                   A * X = B  SUBJECT TO X  >=  0
!  ------------------------------------------------------------------
!                  Subroutine Arguments
!
!  A(), M, N   ON ENTRY, A() CONTAINS THE M BY N MATRIX, A.
!              ON EXIT, A() CONTAINS THE PRODUCT MATRIX, Q*A , WHERE Q IS AN
!              M x M ORTHOGONAL MATRIX GENERATED IMPLICITLY BY THIS SUBROUTINE.
!  B()     ON ENTRY B() CONTAINS THE M-VECTOR, B.   ON EXIT B() CONTAINS Q*B.
!  X()     ON ENTRY X() NEED NOT BE INITIALIZED.
!          ON EXIT X() WILL CONTAIN THE SOLUTION VECTOR.
!  RNORM   ON EXIT RNORM CONTAINS THE EUCLIDEAN NORM OF THE RESIDUAL VECTOR.
!  W()     AN N-ARRAY OF WORKING SPACE.  ON EXIT W() WILL CONTAIN THE DUAL
!          SOLUTION VECTOR.   W WILL SATISFY W(I) = 0. FOR ALL I IN SET P
!          AND W(I) <= 0. FOR ALL I IN SET Z
!  INDX()  AN INTEGER WORKING ARRAY OF LENGTH AT LEAST N.
!          ON EXIT THE CONTENTS OF THIS ARRAY DEFINE THE SETS P AND Z
!          AS FOLLOWS..
!              INDX(1)   THRU INDX(NSETP) = SET P.
!              INDX(IZ1) THRU INDX(IZ2)   = SET Z.
!              IZ1 = NSETP + 1 = NPP1
!              IZ2 = N
!  MODE    THIS IS A SUCCESS-FAILURE FLAG WITH THE FOLLOWING MEANINGS.
!          1   THE SOLUTION HAS BEEN COMPUTED SUCCESSFULLY.
!          2   THE DIMENSIONS OF THE PROBLEM ARE BAD.
!              EITHER M <= 0 OR N <= 0.
!          3   ITERATION COUNT EXCEEDED.  MORE THAN 3*N ITERATIONS.
!
!  ------------------------------------------------------------------
SUBROUTINE nnls (ain, m, n, bin, xout, rnorm, w, indx, mode)
!     ------------------------------------------------------------------
IMPLICIT NONE
integer, parameter :: dp = 8
INTEGER, INTENT(IN)       :: m, n
INTEGER, INTENT(OUT)      :: indx(n), mode
REAL (dp), INTENT(IN OUT) :: ain(m,n), bin(m)
REAL (dp), INTENT(OUT)    :: xout(n), rnorm, w(n)

! NB: the interface statement is not necessary when the subroutines are part of the same module

! INTERFACE
!   SUBROUTINE g1(ain, bin, cterm, sterm, sig)
!     IMPLICIT NONE
!     integer, parameter :: dp = 8
!     REAL (dp), INTENT(IN)  :: ain, bin
!     REAL (dp), INTENT(OUT) :: cterm, sterm, sig
!   END SUBROUTINE g1
!
!   SUBROUTINE h12(mode, lpivot, l1, m, u, up, c, ice, icv, ncv)
!     IMPLICIT NONE
!     integer, parameter :: dp = 8
!     INTEGER, INTENT(IN)                     :: mode, lpivot, l1, m, ice, icv, &
!                                                ncv
!     REAL (dp), DIMENSION(:), INTENT(IN OUT) :: u, c
!     REAL (dp), INTENT(IN OUT)               :: up
!   END SUBROUTINE h12
! END INTERFACE

! Local variables

INTEGER                 :: i, ii, ip, iter, itmax, iz, iz1, iz2, izmax,   &
                           j, jj, jz, l, mda, npp1, nsetp
REAL (dp), DIMENSION(m) :: zz
REAL (dp), DIMENSION(n) :: tempn
REAL (dp), DIMENSION(1) :: dummy
REAL (dp)               :: alpha, asave, cc, factor = 0.01_dp, sm, &
                           ss, t, temp, two = 2.0_dp, unorm, up, wmax,    &
                           zero = 0.0_dp, ztest
!     ------------------------------------------------------------------
mode = 1
IF (m <= 0 .OR. n <= 0) THEN
  mode = 2
  RETURN
END IF
iter = 0
itmax = 3*n

itmax=10*n

!                    INITIALIZE THE ARRAYS indx() AND X().

xout(1:n)=zero
DO i = 1,n
  indx(i) = i
END DO

iz2 = n
iz1 = 1
nsetp = 0
npp1 = 1
!                             ******  MAIN LOOP BEGINS HERE  ******
!                  QUIT IF ALL COEFFICIENTS ARE ALREADY IN THE SOLUTION.
!                        OR IF M COLS OF A HAVE BEEN TRIANGULARIZED.

30 IF (iz1 > iz2 .OR. nsetp >= m) GO TO 350

!         COMPUTE COMPONENTS OF THE DUAL (NEGATIVE GRADIENT) VECTOR W().

!$OMP parallel do private(iz,j) DEFAULT(SHARED)
DO iz = iz1,iz2
  j = indx(iz)
  w(j) = DOT_PRODUCT(ain(npp1:m,j), bin(npp1:m))
END DO
!$OMP end parallel do

!                                   FIND LARGEST POSITIVE W(J).
izmax = 0 ! <- Matt change
60 wmax = zero
DO iz = iz1,iz2
  j = indx(iz)
  IF (w(j) > wmax) THEN
    wmax = w(j)
    izmax = iz
  END IF
END DO

!             IF WMAX  <=  0. GO TO TERMINATION.
!             THIS INDICATES SATISFACTION OF THE KUHN-TUCKER CONDITIONS.

IF (wmax <= zero) GO TO 350
iz = izmax
j = indx(iz)

!     THE SIGN OF W(J) IS OK FOR J TO BE MOVED TO SET P.
!     BEGIN THE TRANSFORMATION AND CHECK NEW DIAGONAL ELEMENT TO AVOID
!     NEAR LINEAR DEPENDENCE.

asave = ain(npp1,j)
CALL h12 (1, npp1, npp1+1, m, ain(:,j), up, dummy, 1, 1, 0)
unorm = zero
IF (nsetp  /=  0) THEN
  unorm = SUM( ain(1:nsetp,j)**2 )
END IF
unorm = SQRT(unorm)
IF (unorm + ABS(ain(npp1,j))*factor - unorm  >  zero) THEN

!        COL J IS SUFFICIENTLY INDEPENDENT.  COPY B INTO ZZ, UPDATE ZZ
!        AND SOLVE FOR ZTEST ( = PROPOSED NEW VALUE FOR X(J) ).

  zz(1:m) = bin(1:m)
  CALL h12 (2, npp1, npp1+1, m, ain(:,j), up, zz, 1, 1, 1)
  ztest = zz(npp1)/ain(npp1,j)

!                                     SEE IF ZTEST IS POSITIVE

  IF (ztest > zero) GO TO 140
END IF

!     REJECT J AS A CANDIDATE TO BE MOVED FROM SET Z TO SET P.
!     RESTORE A(NPP1,J), SET W(J) = 0., AND LOOP BACK TO TEST DUAL
!     COEFFS AGAIN.

ain(npp1,j) = asave
w(j) = zero
GO TO 60

!     THE INDEX  J = indx(IZ)  HAS BEEN SELECTED TO BE MOVED FROM
!     SET Z TO SET P.    UPDATE B,  UPDATE INDICES,  APPLY HOUSEHOLDER
!     TRANSFORMATIONS TO COLS IN NEW SET Z,  ZERO SUBDIAGONAL ELTS IN
!     COL J,  SET W(J) = 0.

140 bin(1:m) = zz(1:m)

indx(iz) = indx(iz1)
indx(iz1) = j
iz1 = iz1+1
nsetp = npp1
npp1 = npp1+1

mda = SIZE(ain,1)
IF (iz1  <=  iz2) THEN
!$O M P parallel do private(jz,jj) DEFAULT(SHARED)
  DO jz = iz1,iz2
    jj = indx(jz)
    CALL h12 (2, nsetp, npp1, m, ain(:,j), up, ain(:,jj), 1, mda, 1)
  END DO
!$O M P end parallel do
END IF

IF (nsetp /= m) THEN
  ain(npp1:m,j) = zero
END IF

w(j) = zero
!                                SOLVE THE TRIANGULAR SYSTEM.
!                                STORE THE SOLUTION TEMPORARILY IN ZZ().
CALL solve_triangular(zz)

!                       ******  SECONDARY LOOP BEGINS HERE ******

!                          ITERATION COUNTER.

210 iter = iter+1
IF (iter > itmax) THEN
  mode = 3
  WRITE (*,'(/a)') ' NNLS quitting on iteration count.'
  GO TO 350
END IF

!                    SEE IF ALL NEW CONSTRAINED COEFFS ARE FEASIBLE.
!                                  IF NOT COMPUTE ALPHA.

alpha = two
DO ip = 1,nsetp
  l = indx(ip)
  IF (zz(ip)  <=  zero) THEN
    t = -xout(l)/(zz(ip)-xout(l))
    IF (alpha > t) THEN
      alpha = t
      jj = ip
    END IF
  END IF
END DO

!          IF ALL NEW CONSTRAINED COEFFS ARE FEASIBLE THEN ALPHA WILL
!          STILL = 2.    IF SO EXIT FROM SECONDARY LOOP TO MAIN LOOP.

IF (abs(alpha-two) .lt. 1.0d-8) GO TO 330

!          OTHERWISE USE ALPHA WHICH WILL BE BETWEEN 0. AND 1. TO
!          INTERPOLATE BETWEEN THE OLD X AND THE NEW ZZ.

!$OMP parallel do private(ip,l) DEFAULT(SHARED)
DO ip = 1,nsetp
  l = indx(ip)
  xout(l) = xout(l) + alpha*(zz(ip)-xout(l))
END DO
!$OMP end parallel do

!        MODIFY A AND B AND THE INDEX ARRAYS TO MOVE COEFFICIENT I
!        FROM SET P TO SET Z.

i = indx(jj)
260 xout(i) = zero

IF (jj /= nsetp) THEN
  jj = jj+1
  DO j = jj,nsetp
    ii = indx(j)
    indx(j-1) = ii
    CALL g1 (ain(j-1,ii), ain(j,ii), cc, ss, ain(j-1,ii))
    ain(j,ii) = zero
    tempn(1:n) = ain(j-1,1:n)
!                 Apply procedure G2 (CC,SS,A(J-1,L),A(J,L))
    IF(ii>1)then
      ain(j-1,1:ii-1) = cc*tempn(1:ii-1) + ss*ain(j,1:ii-1)
      ain(j,1:ii-1)   = -ss*tempn(1:ii-1) + cc*ain(j,1:ii-1)
    ENDIF
    IF(ii<n)then
      ain(j-1,ii+1:n) = cc*tempn(ii+1:n) + ss*ain(j,ii+1:n)
      ain(j,ii+1:n)   = -ss*tempn(ii+1:n) + cc*ain(j,ii+1:n)
    ENDIF
!                 Apply procedure G2 (CC,SS,B(J-1),B(J))
    temp = bin(j-1)
    bin(j-1) = cc*temp + ss*bin(j)
    bin(j)   = -ss*temp + cc*bin(j)
  END DO
END IF

npp1 = nsetp
nsetp = nsetp-1
iz1 = iz1-1
indx(iz1) = i

!        SEE IF THE REMAINING COEFFS IN SET P ARE FEASIBLE.  THEY SHOULD
!        BE BECAUSE OF THE WAY ALPHA WAS DETERMINED.
!        IF ANY ARE INFEASIBLE IT IS DUE TO ROUND-OFF ERROR.  ANY
!        THAT ARE NONPOSITIVE WILL BE SET TO ZERO
!        AND MOVED FROM SET P TO SET Z.

DO jj = 1,nsetp
  i = indx(jj)
  IF (xout(i) <= zero) GO TO 260
END DO

!         COPY B( ) INTO ZZ( ).  THEN SOLVE AGAIN AND LOOP BACK.

zz(1:m) = bin(1:m)
CALL solve_triangular(zz)
GO TO 210
!                      ******  END OF SECONDARY LOOP  ******

330 DO ip = 1,nsetp
  i = indx(ip)
  xout(i) = zz(ip)
END DO
!        ALL NEW COEFFS ARE POSITIVE.  LOOP BACK TO BEGINNING.
GO TO 30

!                        ******  END OF MAIN LOOP  ******

!                        COME TO HERE FOR TERMINATION.
!                     COMPUTE THE NORM OF THE FINAL RESIDUAL VECTOR.

350 sm = zero
IF (npp1 <= m) THEN
  sm = SUM( bin(npp1:m)**2 )
ELSE
  w(1:n) = zero
END IF
rnorm = SQRT(sm)
RETURN

CONTAINS

SUBROUTINE solve_triangular(zz)

!     THE FOLLOWING BLOCK OF CODE IS USED AS AN INTERNAL SUBROUTINE
!     TO SOLVE THE TRIANGULAR SYSTEM, PUTTING THE SOLUTION IN ZZ().

REAL (dp), INTENT(IN OUT) :: zz(:)

DO l = 1,nsetp
  ip = nsetp+1-l
  IF (l  /=  1) zz(1:ip) = zz(1:ip) - ain(1:ip,jj)*zz(ip+1)
  jj = indx(ip)
  zz(ip) = zz(ip) / ain(ip,jj)
END DO

RETURN
END SUBROUTINE solve_triangular

END SUBROUTINE nnls



SUBROUTINE g1(ain, bin, cterm, sterm, sig)

!     COMPUTE ORTHOGONAL ROTATION MATRIX..

!  The original version of this code was developed by
!  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
!  1973 JUN 12, and published in the book
!  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
!  Revised FEB 1995 to accompany reprinting of the book by SIAM.

!     COMPUTE.. MATRIX   (C, S) SO THAT (C, S)(A) = (SQRT(A**2+B**2))
!                        (-S,C)         (-S,C)(B)   (   0          )
!     COMPUTE SIG = SQRT(A**2+B**2)
!        SIG IS COMPUTED LAST TO ALLOW FOR THE POSSIBILITY THAT
!        SIG MAY BE IN THE SAME LOCATION AS A OR B .
!     ------------------------------------------------------------------
IMPLICIT NONE
integer, parameter :: dp = 8
REAL (dp), INTENT(IN)  :: ain, bin
REAL (dp), INTENT(OUT) :: cterm, sterm, sig

!     Local variables
REAL (dp) :: one = 1.0D0, xr, yr, zero = 0.0D0
!     ------------------------------------------------------------------
IF (ABS(ain) > ABS(bin)) THEN
  xr = bin / ain
  yr = SQRT(one + xr**2)
  cterm = SIGN(one/yr, ain)
  sterm = cterm * xr
  sig = ABS(ain) * yr
  RETURN
END IF

IF (abs(bin) .gt. 1.0d-8) THEN
  xr = ain / bin
  yr = SQRT(one + xr**2)
  sterm = SIGN(one/yr, bin)
  cterm = sterm * xr
  sig = ABS(bin) * yr
  RETURN
END IF

!      SIG = ZERO
cterm = zero
sterm = one
RETURN
END SUBROUTINE g1



!     SUBROUTINE h12 (mode, lpivot, l1, m, u, up, c, ice, icv, ncv)

!  CONSTRUCTION AND/OR APPLICATION OF A SINGLE
!  HOUSEHOLDER TRANSFORMATION..     Q = I + U*(U**T)/B

!  The original version of this code was developed by
!  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
!  1973 JUN 12, and published in the book
!  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
!  Revised FEB 1995 to accompany reprinting of the book by SIAM.
!     ------------------------------------------------------------------
!                     Subroutine Arguments

!     MODE   = 1 OR 2   Selects Algorithm H1 to construct and apply a
!            Householder transformation, or Algorithm H2 to apply a
!            previously constructed transformation.
!     LPIVOT IS THE INDEX OF THE PIVOT ELEMENT.
!     L1,M   IF L1  <=  M   THE TRANSFORMATION WILL BE CONSTRUCTED TO
!            ZERO ELEMENTS INDEXED FROM L1 THROUGH M.   IF L1 GT. M
!            THE SUBROUTINE DOES AN IDENTITY TRANSFORMATION.
!     U(),IUE,UP    On entry with MODE = 1, U() contains the pivot
!            vector.  IUE is the storage increment between elements.
!            On exit when MODE = 1, U() and UP contain quantities
!            defining the vector U of the Householder transformation.
!            on entry with MODE = 2, U() and UP should contain
!            quantities previously computed with MODE = 1.  These will
!            not be modified during the entry with MODE = 2.
!     C()    ON ENTRY with MODE = 1 or 2, C() CONTAINS A MATRIX WHICH
!            WILL BE REGARDED AS A SET OF VECTORS TO WHICH THE
!            HOUSEHOLDER TRANSFORMATION IS TO BE APPLIED.
!            ON EXIT C() CONTAINS THE SET OF TRANSFORMED VECTORS.
!     ICE    STORAGE INCREMENT BETWEEN ELEMENTS OF VECTORS IN C().
!     ICV    STORAGE INCREMENT BETWEEN VECTORS IN C().
!     NCV    NUMBER OF VECTORS IN C() TO BE TRANSFORMED. IF NCV  <=  0
!            NO OPERATIONS WILL BE DONE ON C().
!     ------------------------------------------------------------------
SUBROUTINE h12(mode, lpivot, l1, m, u, up, c, ice, icv, ncv)
!     ------------------------------------------------------------------

IMPLICIT NONE
integer, parameter :: dp = 8
INTEGER, INTENT(IN)                     :: mode, lpivot, l1, m, ice, icv, ncv
REAL (dp), DIMENSION(:), INTENT(IN OUT) :: u, c
REAL (dp), INTENT(IN OUT)               :: up

!  Local variables
INTEGER          :: i, i2, i3, i4, incr, j
REAL (dp)        :: b, cl, clinv, one = 1.0D0, sm
!     ------------------------------------------------------------------
IF (0 >= lpivot .OR. lpivot >= l1 .OR. l1 > m) RETURN
cl = ABS(u(lpivot))
IF (mode /= 2) THEN
!                            ****** CONSTRUCT THE TRANSFORMATION. ******
  DO j = l1, m
    cl = MAX(ABS(u(j)),cl)
  END DO
  IF (cl <= 0) RETURN
  clinv = one / cl
  sm = (u(lpivot)*clinv) ** 2 + SUM( (u(l1:m)*clinv)**2 )
  cl = cl * SQRT(sm)
  IF (u(lpivot) > 0) THEN
    cl = -cl
  END IF
  up = u(lpivot) - cl
  u(lpivot) = cl
ELSE
!            ****** APPLY THE TRANSFORMATION  I+U*(U**T)/B  TO C. ******

  IF (cl <= 0) RETURN
END IF
IF (ncv <= 0) RETURN

b = up * u(lpivot)
!                       B  MUST BE NONPOSITIVE HERE.  IF B = 0., RETURN.

IF (b < 0) THEN
  b = one / b
  i2 = 1 - icv + ice * (lpivot-1)
  incr = ice * (l1-lpivot)
  DO j = 1, ncv
    i2 = i2 + icv
    i3 = i2 + incr
    i4 = i3
    sm = c(i2) * up
    DO i = l1, m
      sm = sm + c(i3) * u(i)
      i3 = i3 + ice
    END DO
    IF (abs(sm) .gt. 1.0d-8) THEN
      sm = sm * b
      c(i2) = c(i2) + sm * up
      DO i = l1, m
        c(i4) = c(i4) + sm * u(i)
        i4 = i4 + ice
      END DO
    END IF
  END DO ! j = 1, ncv
END IF

RETURN
END SUBROUTINE h12


end module solver
