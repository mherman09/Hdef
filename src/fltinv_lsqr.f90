!--------------------------------------------------------------------------------------------------!
!--------------------------------- LINEAR LEAST-SQUARES -------------------------------------------!
!--------------------------------------------------------------------------------------------------!


subroutine invert_lsqr()
!----
! Solve for fault slip by directly solving the set of linear equations Ax=b for x.
!----

use io, only: stdout, stderr, verbosity
use solver, only: load_array, load_constraints, solve_dgels, solve_dsysv, solve_dsysv_nrhs, &
                  solve_nnls, solve_dgesv

use fltinv, only: fault, &
                  slip_constraint, &
                  rake_constraint, &
                  euler_file, &
                  npoles, &
                  gf_euler, &
                  displacement, &
                  disp_components, &
                  input_disp_unit, &
                  los, &
                  prestress, &
                  cov_matrix, &
                  gf_disp, &
                  gf_los, &
                  gf_stress, &
                  damping_constant, &
                  smoothing_constant, &
                  nsmooth, &
                  smoothing_pointers, &
                  smoothing_neighbors, &
                  lsqr_mode, &
                  fault_slip, &
                  euler_pole

implicit none

! Local variables
integer :: nrows, ncols, ptr_disp, ptr_los, ptr_stress, ptr_damp, ptr_smooth, ptr_nbr, ptr_euler
integer :: nflt, nslip, ndsp, ndsp_dof, nlos, nsts, nnbr, nfixed, nfree
integer :: i, j, ierr, icomp, iflt, inbr
double precision, allocatable :: A(:,:), b(:), x(:), atmp(:,:), btmp(:)
double precision :: damping_squared, smoothing_squared, factor
logical, allocatable :: isSlipFixed(:)
double precision, parameter :: spy = 60.0d0*60.0d0*24.0d0*365.25d0


! Initialize solution
if (fault%file.ne.'none') then
    if (.not.allocated(fault_slip)) then
        call usage('invert_lsqr: memory for fault slip output array not allocated')
    endif
    fault_slip = 0.0d0
endif
if (euler_file.ne.'none') then
    if (.not.allocated(euler_pole)) then
        call usage('invert_lsqr: memory for Euler pole output array not allocated')
    endif
    euler_pole = 0.0d0
endif


! Determine dimensions of the model matrix, A, and set the locations of components of the array

! Initialize model matrix dimensions to zero
nrows = 0
ncols = 0

! Number of faults and slip degrees of freedom in the inversion -> number of columns in A
nflt = fault%nrows
if (rake_constraint%ncols.eq.1) then
    nslip = nflt
else
    nslip = 2*nflt
endif
ncols = nslip

! Add rigid body rotations to columns
ncols = ncols + 3*npoles ! I GUARANTEE THIS IS GOING TO MAKE SOME ARRAY SIZING ISSUES DOWN THE LINE...


! Initialize sub-array pointers and dimensions
ptr_disp = 0
ptr_los = 0
ptr_stress = 0
ptr_damp = 0
ptr_smooth = 0
ndsp = 0
ndsp_dof = 0
nlos = 0
nsts = 0
ptr_euler = 0

! Row dimensions corresponding to three-component displacement data
if (displacement%file.ne.'none') then
    ptr_disp = nrows + 1
    ndsp = displacement%nrows
    ndsp_dof = len_trim(disp_components)*ndsp
    nrows = nrows + ndsp_dof
endif

! Row dimensions corresponding to LOS displacement data
if (los%file.ne.'none') then
    ptr_los = nrows + 1
    nlos = los%nrows
    nrows = nrows + nlos
endif

! Row dimensions corresponding to pre-stress data
if (prestress%file.ne.'none') then
    ptr_stress = nrows + 1
    nsts = nslip
    nrows = nrows + nsts
endif

! Add rows for damping (one row per model DOF), i.e., damp each slip component separately
if (damping_constant.gt.0.0d0) then
    ptr_damp = nrows + 1
    nrows = nrows + nslip
endif

! Add rows for smoothing, smooth each slip component separately
if (smoothing_constant.gt.0.0d0) then
    ptr_smooth = nrows + 1
    if (rake_constraint%ncols.eq.1) then
        nrows = nrows + nsmooth
    else
        nrows = nrows + 2*nsmooth
    endif
endif

! Columns of Euler pole GFs start after fault slip columns
ptr_euler = nslip + 1

! Resulting dimensions and locations of sub-arrays
if (verbosity.ge.2) then
    write(stdout,*) 'invert_lsqr: model array parameters'
    write(stdout,*) 'nrows:      ', nrows
    write(stdout,*) 'ncolumns:   ', ncols
    write(stdout,*) 'ptr_disp:   ', ptr_disp
    write(stdout,*) 'ptr_los:    ', ptr_los
    write(stdout,*) 'ptr_stress: ', ptr_stress
    write(stdout,*) 'ptr_damp:   ', ptr_damp
    write(stdout,*) 'ptr_smooth: ', ptr_smooth
    write(stdout,*) 'ptr_euler:  ', ptr_euler
endif

! Allocate memory to model matrix, A, and data vector, b
if (.not.allocated(A)) then
    allocate(A(nrows,ncols),stat=ierr)
    if (ierr.ne.0) then
        call usage('invert_lsqr: error allocating memory to A')
    endif
    if (verbosity.ge.2) then
        write(stdout,*) 'invert_lsqr: memory allocated to A(',nrows,',',ncols,')'
    endif
endif
if (.not.allocated(b)) then
    allocate(b(nrows),stat=ierr)
    if (ierr.ne.0) then
        call usage('invert_lsqr: error allocating memory to b')
    endif
    if (verbosity.ge.2) then
        write(stdout,*) 'invert_lsqr: memory allocated to b(',nrows,')'
    endif
endif
A = 0.0d0
b = 0.0d0


! Load the Green's functions and regularization arrays into the model matrix, A, and observations
! into the data vector, b

! Load three-component displacement GFs and data
if (displacement%file.ne.'none') then
    do i = 1,len_trim(disp_components)
        read(disp_components(i:i),*) icomp
        call load_array(A,nrows,ncols, &
                        gf_disp%array((icomp-1)*ndsp+1:icomp*ndsp,1:nslip),ndsp,nslip, &
                        (i-1)*ndsp+1,1,'gf_disp%array',ierr)
        if (ierr.ne.0) then
            call usage('invert_lsqr: error loading three-component displacement GFs into A')
        endif
        call load_array(b,nrows,displacement%array(:,icomp+3),ndsp,(i-1)*ndsp+1, &
                        'displacement%array',ierr)
        if (ierr.ne.0) then
            call usage('invert_lsqr: error loading three-component displacement observations into b')
        endif
    enddo
endif

! Load LOS displacement GFs and data (NEED TO INCORPORATE LOS WEIGHTS!!!!! <-THIS CAN BE DONE W COV MATRIX)
if (los%file.ne.'none') then
    call load_array(A,nrows,ncols,gf_los%array(:,1:nslip),nlos,nslip,ptr_los,1,'gf_los%array',ierr)
    if (ierr.ne.0) then
        call usage('invert_lsqr: error loading line-of-sight displacement GFs into A')
    endif
    call load_array(b,nrows,los%array(:,4),nlos,ptr_los,'los%array',ierr)
    if (ierr.ne.0) then
        call usage('invert_lsqr: error loading line-of-sight displacement observations into b')
    endif
endif

! Load rigid body rotations into model matrix
if (euler_file.ne.'none') then
    do i = 1,len_trim(disp_components)
        read(disp_components(i:i),*) icomp
        if (icomp.ne.3) then
            call load_array(A,nrows,ncols, &
                            gf_euler((icomp-1)*ndsp+1:icomp*ndsp,1:3*npoles),ndsp,3*npoles, &
                            (i-1)*ndsp+1,ptr_euler,'gf_euler%array',ierr)
            if (ierr.ne.0) then
                call usage('invert_lsqr: error loading rigid body rotation GFs into A')
            endif
        endif
    enddo
endif

! Load covariance matrix for displacements
if (displacement%file.ne.'none'.or.los%file.ne.'none') then
    allocate(atmp(ndsp_dof+nlos,nslip),stat=ierr)
    if (ierr.ne.0) then
        call usage('invert_lsqr: error allocating memory to atmp for including covariance matrix')
    endif
    allocate(btmp(ndsp_dof+nlos),stat=ierr)
    if (ierr.ne.0) then
        call usage('invert_lsqr: error allocating memory to btmp for including covariance matrix')
    endif

    ! Compute cov_matrix^-1*A and cov_matrix^1*b
    call solve_dsysv_nrhs(cov_matrix,A(1:ndsp_dof+nlos,1:nslip),atmp,ndsp_dof+nlos,nslip,ierr)
    if (ierr.ne.0) then
        call usage('invert_lsqr: error computing cov_matrix^-1*A')
    endif
    call solve_dsysv(cov_matrix,b(1:ndsp_dof+nlos),btmp,ndsp_dof+nlos,ierr)
    if (ierr.ne.0) then
        call usage('invert_lsqr: error computing cov_matrix^-1*b')
    endif
    A(1:ndsp_dof+nlos,1:nslip) = atmp
    b(1:ndsp_dof+nlos) = btmp

    deallocate(atmp)
    deallocate(btmp)
endif


! Load pre-stress data
if (prestress%file.ne.'none') then
    call load_array(A,nrows,ncols,&
                    gf_stress%array(1:nslip,1:nslip),nslip,nslip,&
                    ptr_stress,1,'gf_stress%array',ierr)
    if (ierr.ne.0) then
        call usage('error loading shear stress GFs into model matrix')
    endif
    call load_array(b,nrows,prestress%array(1:nflt,1),nflt,ptr_stress,'prestress%array(:,1)',ierr)
    if (rake_constraint%ncols.eq.2) then
        call load_array(b,nrows,prestress%array(1:nflt,2),nflt,ptr_stress+nflt, &
                        'prestress%array(:,2)',ierr)
    endif
    if (ierr.ne.0) then
        call usage('error loading pre-stresses into b vector')
    endif
endif

! Minimize model length
if (damping_constant.gt.0.0d0) then
    A(ptr_damp:ptr_damp+nslip-1,1:nslip) = 0.0d0
    damping_squared = damping_constant*damping_constant
    do i = 1,nslip
        A(ptr_damp+i-1,i) = damping_squared
    enddo
    b(ptr_damp:ptr_damp+nslip-1) = 0.0d0
endif

! Minimize model smoothness
if (smoothing_constant.gt.0.0d0) then
    A(ptr_smooth:ptr_smooth+nsmooth-1,1:nslip) = 0.0d0
    if (rake_constraint%file.eq.'none'.or.rake_constraint%ncols.eq.2) then
        A(ptr_smooth+nsmooth:ptr_smooth+2*nsmooth-1,1:nslip) = 0.0d0
    endif
    smoothing_squared = smoothing_constant*smoothing_constant
    do i = 1,nsmooth
        iflt = smoothing_pointers(i,1)
        nnbr = smoothing_pointers(i,2)
        ptr_nbr = smoothing_pointers(i,3)
        A(ptr_smooth+i-1,iflt) = dble(nnbr)*smoothing_squared
        if (rake_constraint%file.eq.'none'.or.rake_constraint%ncols.eq.2) then
            A(ptr_smooth+nsmooth+i-1,iflt+nflt) = dble(nnbr)*smoothing_squared
        endif
        do j = 1,nnbr
            inbr = smoothing_neighbors(ptr_nbr+j-1)
            A(ptr_smooth+i-1,inbr) = -1.0d0*smoothing_squared
            if (rake_constraint%file.eq.'none'.or.rake_constraint%ncols.eq.2) then
                A(ptr_smooth+nsmooth+i-1,inbr+nflt) = -1.0d0*smoothing_squared
            endif
        enddo
    enddo
    b(ptr_smooth:ptr_smooth+nsmooth-1) = 0.0d0
    if (rake_constraint%file.eq.'none'.or.rake_constraint%ncols.eq.2) then
        b(ptr_smooth+nsmooth:ptr_smooth+2*nsmooth-1) = 0.0d0
    endif
endif
if (verbosity.ge.2) then
    write(stdout,*) 'invert_lsqr: A matrix and b vector loaded'
endif
if (verbosity.ge.3) then
    write(stdout,*) 'A:'
    do i = 1,nrows
        write(stdout,*) A(i,:)
    enddo
    write(stdout,*) 'b:'
    do i = 1,nrows
        write(stdout,*) b(i)
    enddo
endif


! Fixed slip constraints
allocate(isSlipFixed(nslip),stat=ierr)
if (ierr.ne.0) then
    call usage('invert_lsqr: error allocating memory to isSlipFixed')
endif
isSlipFixed = .false.

if (slip_constraint%file.ne.'none') then
    ! Set faults that have fixed slip to T in array
    do i = 1,nflt
        if (abs(slip_constraint%array(i,1)).lt.99998.0d0) then
            isSlipFixed(i) = .true.
            ! write(0,*) 'fault ',i,' is fixed slip'
        endif
        if (rake_constraint%ncols.eq.2) then
            if (abs(slip_constraint%array(i,2)).lt.99998.0d0) then
                isSlipFixed(i+nflt) = .true.
            endif
        endif
    enddo

    ! Load these constraints into A and b vectors
    if (rake_constraint%ncols.eq.1) then
        call load_constraints(A,b,nrows,ncols, &
                              slip_constraint%array(:,1),isSlipFixed, &
                              ierr)
    elseif (rake_constraint%ncols.eq.2) then
        call load_constraints(A,b,nrows,ncols, &
                              [slip_constraint%array(:,1),slip_constraint%array(:,2)],isSlipFixed, &
                              ierr)
    else
        call usage('invert_lsqr: number of columns in rake constraints must be 1 or 2')
    endif

    ! Remove stress rows for fixed slip faults; these no longer contribute to the solution
    if (prestress%file.ne.'none') then
        nfree = 0
        nfixed = 0
        do i = 1,nslip
            if (isSlipFixed(i)) then
                ! Do not add row i to A
                nfixed = nfixed + 1
            else
                ! Add row i to A
                nfree = nfree + 1
                A(ptr_stress+nfree-1,:) = A(ptr_stress+i-1,:)
                b(ptr_stress+nfree-1) = b(ptr_stress+i-1)
            endif
        enddo

        ! Shift rows of A and b past stress block down
        if (nrows.gt.ptr_stress+nslip-1) then
            A(ptr_stress+nfree-1:nrows-nfixed,:) = A(ptr_stress+nslip:nrows,:)
            b(ptr_stress+nfree-1:nrows-nfixed) = b(ptr_stress+nslip:nrows)
        elseif (nrows.eq.ptr_stress+nslip-1) then
            ! Nothing to shift
        else
            call usage('invert_lsqr: error in indexing rows of model matrix')
        endif
        nrows = nrows - nfixed
    endif


    if (verbosity.ge.3) then
        write(stdout,*) 'After loading slip constraints:'
        write(stdout,*) 'A:'
        do i = 1,nrows
            write(stdout,*) A(i,:)
        enddo
        write(stdout,*) 'b:'
        do i = 1,nrows
            write(stdout,*) b(i)
        enddo
    endif
endif

if (nrows.eq.0) then
    call usage('invert_lsqr: no rows in model matrix A')
endif
if (ncols.eq.0) then
    call usage('invert_lsqr: no columns in model matrix A')
endif


! Solve Ax=b for x
allocate(x(ncols))
if (lsqr_mode.eq.'gels') then
    call solve_dgels(A(1:nrows,1:ncols),b(1:nrows),x,nrows,ncols,ierr)
    if (ierr.ne.0) then
        call usage('invert_lsqr: error inverting for slip with gels')
    endif

elseif (lsqr_mode.eq.'nnls') then
    call solve_nnls(A(1:nrows,1:ncols),b(1:nrows),x,nrows,ncols,ierr)
    if (ierr.ne.0) then
        call usage('invert_lsqr: error inverting for slip with nnls')
    endif

elseif (lsqr_mode.eq.'gesv') then
    if (nrows.ne.ncols) then
        write(stderr,*) 'invert_lsqr: nrows not equal to ncols'
        write(stderr,*) '    nrows=',nrows
        write(stderr,*) '    ncols=',ncols
        call usage('Inversion routine gesv requires a square A matrix')
    endif
    call solve_dgesv(A(1:nrows,1:ncols),b(1:nrows),x,nrows,ierr)
    if (ierr.ne.0) then
        call usage('invert_lsqr: error inverting for slip with gesv')
    endif

else
    write(stderr,*) 'invert_lsqr: no lsqr mode named "'//trim(lsqr_mode)//'"'
    write(stderr,*) 'Available options:'
    write(stderr,*) '    gels'
    write(stderr,*) '    nnls'
    call usage(     '    gesv')
endif

! Load solution into fault slip vector
j = 1
do iflt = 1,nflt
    if (isSlipFixed(iflt)) then
        fault_slip(iflt,1) = slip_constraint%array(iflt,1)
    else
        fault_slip(iflt,1) = x(j)
        j = j + 1
    endif
enddo
if (rake_constraint%ncols.eq.2) then
    do iflt = 1,nflt
        if (isSlipFixed(iflt+nflt)) then
            fault_slip(iflt,2) = slip_constraint%array(iflt,2)
        else
            fault_slip(iflt,2) = x(j)
            j = j + 1
        endif
    enddo
endif

! Load solution into euler pole vector, scaling by correct units
if (npoles.gt.0) then
    factor = 0.0d0
    if (input_disp_unit.eq.'m/s') then
        factor = 1.0d0
    elseif (input_disp_unit.eq.'m/yr') then
        factor = 1.0d0*spy
    elseif (input_disp_unit.eq.'mm/s') then
        factor = 1.0d3
    elseif (input_disp_unit.eq.'mm/yr') then
        factor = 1.0d3*spy
    else
        call usage('invert_lsqr: unit '//trim(input_disp_unit)//' not compatible')
    endif
    do i = 1,npoles
        euler_pole(i,1) = x(j+i-1)/factor
        euler_pole(i,2) = x(j+i  )/factor
        euler_pole(i,3) = x(j+i+1)/factor
    enddo
endif

return
end subroutine invert_lsqr
