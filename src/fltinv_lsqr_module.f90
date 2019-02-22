module lsqr_module
    double precision, allocatable :: A(:,:), Asave(:,:)
    double precision, allocatable :: b(:,:)
    double precision, allocatable :: x(:,:)
    integer :: nrows
    integer :: ncols
    integer :: ptr_disp
    integer :: ptr_los
    integer :: ptr_stress
    integer :: ptr_damp
    integer :: ptr_smooth
    character(len=8) :: lsqr_mode
    logical :: isAsaveLoaded

!--------------------------------------------------------------------------------------------------!
contains
!--------------------------------------------------------------------------------------------------!

    subroutine invert_lsqr()
    !----
    ! Invert displacements and pre-stresses for fault slip using a direct or linear least squares approach
    !----
    use io, only : stderr, verbosity
    use variable_module, only: fault, slip_constraint, rake_constraint, fault_slip, &
                               prestress, displacement, los
    use variable_module, only: inversion_mode
    implicit none
    ! Local variables
    integer :: i, j, nzeros

    if (verbosity.ge.2) then
        write(stderr,'(A)') 'invert_lsqr: starting'
    endif

    ! Determine dimensions of model matrix
    call calc_array_dimensions()

    ! Allocate memory to model matrix, constraint vector, and solution vector
    if (.not.allocated(A)) then
        allocate(A(nrows,ncols))
    endif
    if (.not.allocated(b)) then
        allocate(b(nrows,1))
    endif
    A = 0.0d0
    b = 0.0d0

    ! In annealing search for pseudo-coupling, it saves substantial effort to only load the full
    ! A matrix once, so initialize a matrix Asave
    if (inversion_mode.eq.'anneal-psc') then
        if (.not.allocated(Asave)) then
            allocate(Asave(nrows,ncols))
        endif
    else
        isAsaveLoaded = .false.
    endif

    ! Load model matrix, A (and maybe Asave), and constraint vector, b
    call load_arrays()
    open(unit=55,file='a.dat',status='unknown')
    open(unit=56,file='b.dat',status='unknown')
    do i = 1,nrows
        write(55,*) (A(i,j),j=1,ncols)
    enddo
    do i = 1,nrows
        write(56,*) b(i,1)
    enddo
    close(55)
    close(56)

    if (nrows.eq.0.and.ncols.eq.0) then
        write(0,*) 'invert_lsqr: no rows or columns in A matrix'
        return
    endif

    ! Solve generalized least squares problem
    if (.not.allocated(x)) then
        allocate(x(ncols,1))
    endif

    if (inversion_mode.eq.'anneal-psc'.or. &
           (prestress%file.ne.'none'.and.displacement%file.eq.'none'.and.los%file.eq.'none')) then
        ! In solving for slip around a locked patch (pseudo-coupling), the system is
        ! properly determined so we can use direct inversion methods rather than linear least squares.
#ifdef USESUPERLU
        call count_matrix_zeros(nzeros)
        if (dble(nzeros)/dble(nrows*ncols).gt.0.3d0) then
            call solve_dgssv() ! Use sparse solver (SuperLU)
        else
#endif
            call solve_dgesv() ! Use dense solver (LAPACK)
#ifdef USESUPERLU
        endif
#endif
    elseif (lsqr_mode.eq.'gels'.or.lsqr_mode.eq.'dgels') then
        call solve_lsqr_dgels()
    elseif (lsqr_mode.eq.'nnls') then
        call solve_lsqr_nnls()
    else
        call usage('!! Error: no lsqr_mode named '//trim(lsqr_mode))
    endif


    ! Load solution into fault_slip array
    if (.not.allocated(fault_slip)) then
        allocate(fault_slip(fault%nrecords,2))
    endif
    j = 1
    do i = 1,fault%nrecords
        if (slip_constraint%file.ne.'none') then
            if (dabs(slip_constraint%array(i,1)).lt.99998.0d0) then
                fault_slip(i,1) = slip_constraint%array(i,1)
            else
                fault_slip(i,1) = x(j,1)
                j = j + 1
            endif
        else
            fault_slip(i,1) = x(j,1)
            j = j + 1
        endif
    enddo
    if (rake_constraint%file.eq.'none'.or.rake_constraint%nfields.eq.2) then
        do i = 1,fault%nrecords
            if (slip_constraint%file.ne.'none') then
                if (dabs(slip_constraint%array(i,2)).lt.99998.0d0) then
                    fault_slip(i,2) = slip_constraint%array(i,2)
                else
                    fault_slip(i,2) = x(j,1)
                    j = j + 1
                endif
            else
                fault_slip(i,2) = x(j,1)
                j = j + 1
            endif
        enddo
    endif

    ! FREEEEEDDDOOOOOMMMMMMM
    if (inversion_mode.eq.'lsqr') then
        if (allocated(A)) then
            deallocate(A)
        endif
        if (allocated(b)) then
            deallocate(b)
        endif
        if (allocated(x)) then
            deallocate(x)
        endif
    endif

    if (verbosity.ge.2) then
        write(stderr,'(A)') 'invert_lsqr says: finished'
        write(stderr,*)
    endif

    return
    end subroutine invert_lsqr

!--------------------------------------------------------------------------------------------------!

    subroutine calc_array_dimensions()
    use io, only: stderr, verbosity
    use variable_module, only: displacement, prestress, los, fault, rake_constraint, &
                                damping_constant, smoothing_constant, smoothing, &
                                disp_components
    implicit none

    if (verbosity.ge.2) then
        write(stderr,'(A)') 'calc_array_dimensions says: starting'
    endif

    !  Initialize model matrix dimensions
    nrows = 0
    ncols = 0

    ! Initialize array pointers
    ptr_disp = 0
    ptr_los = 0
    ptr_stress = 0
    ptr_damp = 0
    ptr_smooth = 0

    ! Row dimensions corresponding to three-component displacement data
    if (displacement%file.ne.'none') then
        ptr_disp = nrows + 1
        nrows = nrows + len_trim(disp_components)*displacement%nrecords
    endif

    ! Row dimensions corresponding to LOS displacement data
    if (los%file.ne.'none') then
        ptr_los = nrows + 1
        nrows = nrows + los%nrecords
    endif

    ! Row dimensions corresponding to pre-stress data
    if (prestress%file.ne.'none') then
        ptr_stress = nrows + 1
        nrows = nrows + 2*fault%nrecords
    endif

    ! Column dimensions corresponding to fault slip
    if (rake_constraint%file.eq.'none'.or.rake_constraint%nfields.eq.2) then
        ncols = 2*fault%nrecords
    else
        ncols = fault%nrecords
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
        if (rake_constraint%file.eq.'none'.or.rake_constraint%nfields.eq.2) then
            nrows = nrows + 2*smoothing%nrecords
        else
            nrows = nrows + smoothing%nrecords
        endif
    endif

    if (verbosity.ge.2) then
        write(stderr,'(A)') 'calc_array_dimensions says: finished'
        write(stderr,'(A)') 'Model matrix parameters'
        write(stderr,'(A,I8)') 'nrows:      ', nrows
        write(stderr,'(A,I8)') 'ncolumns:   ', ncols
        write(stderr,'(A,I8)') 'ptr_disp:   ', ptr_disp
        write(stderr,'(A,I8)') 'ptr_los:    ', ptr_los
        write(stderr,'(A,I8)') 'ptr_stress: ', ptr_stress
        write(stderr,'(A,I8)') 'ptr_damp:   ', ptr_damp
        write(stderr,'(A,I8)') 'ptr_smooth: ', ptr_smooth
        write(stderr,*)
    endif

    return
    end subroutine calc_array_dimensions

!--------------------------------------------------------------------------------------------------!

    subroutine load_arrays()
    use io, only: stderr, verbosity
    use variable_module, only: displacement, los, prestress, &
                                slip_constraint, damping_constant, smoothing_constant
    implicit none
    ! Local variables
    integer :: i, j, nzeros
    character(len=256) :: fmt

    if (verbosity.ge.2) then
        write(stderr,'(A)') 'load_arrays says: starting'
    endif

    if (.not.isAsaveLoaded) then

        ! If we have not yet loaded the A matrix, then do it here

        ! Load three-component displacement data
        if (displacement%file.ne.'none') then
            call load_displacements()
        endif

        ! Load LOS displacement data
        if (los%file.ne.'none') then
            call load_los()
        endif

        ! Load pre-stress data
        if (prestress%file.ne.'none') then
            call load_prestresses()
       endif

        ! Minimize model length
        if (damping_constant.gt.0.0d0) then
            call load_damping()
        endif

        ! Minimize model smoothness
        if (smoothing_constant.gt.0.0d0) then
            call load_smoothing()
        endif

        Asave = A
        isAsaveLoaded = .true.

    else
        ! If we have already loaded the A matrix and saved it as Asave, copy it to A
        A = Asave
    endif

    ! Apply constraints
    if (slip_constraint%file.ne.'none') then
        call load_slip_constraints()
    endif

    if (nrows.eq.0) then
        call usage('!! Error: no rows in least-squares matrix')
    endif

    if (ncols.eq.0) then
        call usage('!! Error: no columns in least-squares matrix')
    endif

    if (verbosity.ge.2) then
        write(stderr,'(A)') 'load_arrays says: finished'
    endif
    if (verbosity.ge.3) then
        write(stderr,'("A =")')
        write(fmt,9001) ncols
    9001 format('(1P',I6,'E12.4)')
        do i = 1,nrows
            if (i.eq.ptr_disp) write(stderr,'("displacements")')
            if (i.eq.ptr_los) write(stderr,'("LOS displacements")')
            if (i.eq.ptr_stress) write(stderr,'("pre-stresses")')
            if (i.eq.ptr_damp) write(stderr,'("damping")')
            if (i.eq.ptr_smooth) write(stderr,'("smoothing")')
            write(0,fmt) (A(i,j),j=1,ncols)
        enddo
        write(stderr,*)
        write(stderr,'("b =")')
        do i = 1,nrows
            if (i.eq.ptr_disp) write(stderr,'("displacements")')
            if (i.eq.ptr_los) write(stderr,'("LOS displacements")')
            if (i.eq.ptr_stress) write(stderr,'("pre-stresses")')
            if (i.eq.ptr_damp) write(stderr,'("damping")')
            if (i.eq.ptr_smooth) write(stderr,'("smoothing")')
            write(stderr,'(1PE12.4)') b(i,1)
        enddo
    endif
    if (verbosity.ge.2) then
        nzeros = 0
        do i = 1,nrows
            do j = 1,ncols
                if (dabs(A(i,j)).lt.1.0d-20) then
                    nzeros = nzeros + 1
                endif
            enddo
        enddo
        write(0,*) 'load_arrays says: A matrix has ',nzeros,' zeros out of ',nrows*ncols,' total'
    endif
    if (verbosity.ge.2) then
        write(stderr,*)
    endif

    return
    end subroutine load_arrays

    !--------------------------------------------------------------------------------------------------!

    subroutine load_displacements()
    use io, only: stderr, verbosity
    use variable_module, only: displacement, fault, gf_disp, rake_constraint, disp_components
    implicit none
    integer :: i, j, ndsp, nflt
    !
    if (verbosity.ge.2) then
        write(stderr,'(A)') 'load_displacements says: starting'
    endif

    ndsp = displacement%nrecords
    nflt = fault%nrecords
    do i = 1,ndsp
        do j = 1,nflt
            ! For selected components of displacements, load strike-slip (or constant rake)
            ! source GFs into model matrix
            if (disp_components.eq.'123') then
                A(i       ,j) = gf_disp%array(i       ,j) ! strike-slip, x-displacement
                A(i+  ndsp,j) = gf_disp%array(i+  ndsp,j) ! strike-slip, y-displacement
                A(i+2*ndsp,j) = gf_disp%array(i+2*ndsp,j) ! strike-slip, z-displacement
            elseif (disp_components.eq.'12') then
                A(i       ,j) = gf_disp%array(i       ,j) ! strike-slip, x-displacement
                A(i+  ndsp,j) = gf_disp%array(i+  ndsp,j) ! strike-slip, y-displacement
            elseif (disp_components.eq.'13') then
                A(i       ,j) = gf_disp%array(i,j)        ! strike-slip, x-displacement
                A(i+  ndsp,j) = gf_disp%array(i+2*ndsp,j) ! strike-slip, z-displacement
            elseif (disp_components.eq.'23') then
                A(i       ,j) = gf_disp%array(i+  ndsp,j) ! strike-slip, y-displacement
                A(i+  ndsp,j) = gf_disp%array(i+2*ndsp,j) ! strike-slip, z-displacement
            elseif (disp_components.eq.'1') then
                A(i       ,j) = gf_disp%array(i      ,j)  ! strike-slip, x-displacement
            elseif (disp_components.eq.'2') then
                A(i       ,j) = gf_disp%array(i+  ndsp,j) ! strike-slip, y-displacement
            elseif (disp_components.eq.'3') then
                A(i       ,j) = gf_disp%array(i+2*ndsp,j) ! strike-slip, z-displacement
            endif

            ! Choose columns to load selected fault slip components
            if (rake_constraint%file.eq.'none'.or.rake_constraint%nfields.eq.2) then
                ! For selected components of displacements, load dip-slip source GFs into model matrix
                if (disp_components.eq.'123') then
                    A(i       ,j+nflt) = gf_disp%array(i       ,j+nflt) ! dip-slip, x-displacement
                    A(i+  ndsp,j+nflt) = gf_disp%array(i+  ndsp,j+nflt) ! dip-slip, y-displacement
                    A(i+2*ndsp,j+nflt) = gf_disp%array(i+2*ndsp,j+nflt) ! dip-slip, z-displacement
                elseif (disp_components.eq.'12') then
                    A(i       ,j+nflt) = gf_disp%array(i       ,j+nflt) ! dip-slip, x-displacement
                    A(i+  ndsp,j+nflt) = gf_disp%array(i+  ndsp,j+nflt) ! dip-slip, y-displacement
                elseif (disp_components.eq.'13') then
                    A(i       ,j+nflt) = gf_disp%array(i       ,j+nflt) ! dip-slip, x-displacement
                    A(i+  ndsp,j+nflt) = gf_disp%array(i+2*ndsp,j+nflt) ! dip-slip, z-displacement
                elseif (disp_components.eq.'23') then
                    A(i       ,j+nflt) = gf_disp%array(i+  ndsp,j+nflt) ! dip-slip, y-displacement
                    A(i+  ndsp,j+nflt) = gf_disp%array(i+2*ndsp,j+nflt) ! dip-slip, z-displacement
                elseif (disp_components.eq.'1') then
                    A(i       ,j+nflt) = gf_disp%array(i       ,j+nflt) ! dip-slip, x-displacement
                elseif (disp_components.eq.'2') then
                    A(i       ,j+nflt) = gf_disp%array(i+  ndsp,j+nflt) ! dip-slip, y-displacement
                elseif (disp_components.eq.'3') then
                    A(i       ,j+nflt) = gf_disp%array(i+2*ndsp,j+nflt) ! dip-slip, z-displacement
                endif
            endif
        enddo

        ! Load selected components of observed displacements into constraint vector
        if (disp_components.eq.'123') then
            b(i       ,1) = displacement%array(i,4) ! x
            b(i+  ndsp,1) = displacement%array(i,5) ! y
            b(i+2*ndsp,1) = displacement%array(i,6) ! z
        elseif (disp_components.eq.'12') then
            b(i       ,1) = displacement%array(i,4) ! x
            b(i+  ndsp,1) = displacement%array(i,5) ! y
        elseif (disp_components.eq.'13') then
            b(i       ,1) = displacement%array(i,4) ! x
            b(i+  ndsp,1) = displacement%array(i,6) ! z
        elseif (disp_components.eq.'23') then
            b(i       ,1) = displacement%array(i,5) ! y
            b(i+  ndsp,1) = displacement%array(i,6) ! z
        elseif (disp_components.eq.'1') then
            b(i       ,1) = displacement%array(i,4) ! x
        elseif (disp_components.eq.'2') then
            b(i       ,1) = displacement%array(i,5) ! y
        elseif (disp_components.eq.'3') then
            b(i       ,1) = displacement%array(i,6) ! z
        endif
    enddo

    if (verbosity.ge.2) then
        write(stderr,'(A)') 'load_displacements says: finished'
        write(stderr,*)
    endif

    return
    end subroutine load_displacements

!--------------------------------------------------------------------------------------------------!

    subroutine load_los()
    use io, only: stderr, verbosity
    use variable_module, only: los, fault, gf_los, rake_constraint, los_weight
    implicit none
    integer :: i, j, ndsp, nflt
    !
    if (verbosity.ge.2) then
        write(stderr,'(A)') 'load_los says: starting'
    endif

    ndsp = los%nrecords
    nflt = fault%nrecords
    do i = 1,ndsp
        do j = 1,nflt
            ! Load strike-slip (or constant rake) source GFs into model matrix
            A(ptr_los+i-1,j) = los_weight*gf_los%array(i,j) ! strike-slip displacement

            ! Load dip-slip source GFs into model matrix
            if (rake_constraint%file.eq.'none'.or.rake_constraint%nfields.eq.2) then
                A(ptr_los+i-1,j+nflt) = los_weight*gf_los%array(i,j+nflt) ! dip-slip displacement
            endif
        enddo

        ! Load observed displacements into constraint vector
        b(ptr_los+i-1,1) = los_weight*los%array(i,4) ! observed LOS displacement
    enddo

    if (verbosity.ge.2) then
        write(stderr,'(A)') 'load_los says: finished'
        write(stderr,*)
    endif

    return
    end subroutine load_los

!--------------------------------------------------------------------------------------------------!

    subroutine load_prestresses()
    use io, only: stderr, verbosity
    use variable_module, only : fault, prestress, rake_constraint, gf_stress, stress_weight
    implicit none
    ! Local variables
    integer :: i, j, nflt

    nflt = fault%nrecords

    if (verbosity.ge.2) then
        write(stderr,'(A)') 'load_prestresses says: starting'
    endif

    do i = 1,nflt
        do j = 1,nflt
            if (rake_constraint%file.eq.'none'.or.rake_constraint%nfields.eq.2) then
                A(ptr_stress     +i-1,j     ) = stress_weight*gf_stress%array(i     ,j     ) ! ss src -> ss traction
                A(ptr_stress+nflt+i-1,j     ) = stress_weight*gf_stress%array(i+nflt,j     ) ! ss src -> ds traction
                A(ptr_stress     +i-1,j+nflt) = stress_weight*gf_stress%array(i     ,j+nflt) ! ds src -> ss traction
                A(ptr_stress+nflt+i-1,j+nflt) = stress_weight*gf_stress%array(i+nflt,j+nflt) ! ds src -> ds traction
            else
                A(ptr_stress     +i-1,j     ) = stress_weight*gf_stress%array(i     ,j     ) ! ss src -> ss traction
                A(ptr_stress+nflt+i-1,j     ) = stress_weight*gf_stress%array(i+nflt,j     ) ! ss src -> ds traction
            endif
        enddo

        ! Load pre-stresses into constraint vector
        b(ptr_stress     +i-1,1) = stress_weight*prestress%array(i,1) ! ss traction
        b(ptr_stress+nflt+i-1,1) = stress_weight*prestress%array(i,2) ! ds traction
    enddo
    !
    if (verbosity.ge.2) then
        write(stderr,'(A)') 'load_prestresses says: finished'
        write(stderr,*)
    endif

    return
    end subroutine load_prestresses

!--------------------------------------------------------------------------------------------------!

    subroutine load_damping()
    !----
    ! Add equations to model matrix to minimize the length of the solution to Ax = b
    ! with weighting factor damping_constant*damping_constant.
    !----
    use io, only : stderr, verbosity
    use variable_module, only: damping_constant
    implicit none
    ! Local variables
    integer :: i, j
    double precision :: damping_squared

    if (verbosity.ge.2) then
        write(stderr,'(A)') "load_damping says: starting"
    endif

    damping_squared = damping_constant*damping_constant

    do i = 0,ncols-1
        b(ptr_damp+i,1) = 0.0d0
        do j = 0,ncols-1
            if (i.eq.j) then
                A(ptr_damp+i,j+1) = damping_squared
            else
                A(ptr_damp+i,j+1) = 0.0d0
            endif
        enddo
    enddo

    if (verbosity.ge.2) then
        write(stderr,'(A)') "load_damping says: finished"
        write(stderr,*)
    endif

    return
    end subroutine load_damping

!--------------------------------------------------------------------------------------------------!

    subroutine load_smoothing()
    !----
    ! Add equations to model matrix to minimize the curvature (second derivative)
    ! of the solution to Ax = b with weighting factor smooth.
    !----
    use io, only : stderr, verbosity
    use variable_module, only : smoothing, smoothing_constant, smoothing_neighbors, &
                                fault, rake_constraint
    implicit none
    ! Local variables
    integer :: i, j, ifault, nneighbor, ineighbor, neighbor, nsmooth
    double precision :: smoothing_squared

    if (verbosity.ge.2) then
        write(stderr,'(A)') "load_smoothing says: starting"
    endif

    smoothing_squared = smoothing_constant*smoothing_constant
    nsmooth = smoothing%nrecords

    ! Read smoothing element linking file and load into A
    do i = 1,smoothing%nrecords
        ! Set observation corresponding to smoothing row to zero
        b(ptr_smooth+i-1,1) = 0.0d0
        if (rake_constraint%file.eq.'none'.or.rake_constraint%nfields.eq.2) then
            b(ptr_smooth+nsmooth+i-1,1) = 0.0d0
        endif

        ifault = smoothing%intarray(i,1)
        nneighbor = smoothing%intarray(i,2)
        ineighbor = smoothing%intarray(i,3)

        ! Fault to be smoothed gets weight=nneighbor
        A(ptr_smooth+i-1,ifault) = dble(nneighbor)*smoothing_squared
        if (rake_constraint%file.eq.'none'.or.rake_constraint%nfields.eq.2) then
            A(ptr_smooth+nsmooth+i-1,ifault+fault%nrecords) = dble(nneighbor)*smoothing_squared
        endif

        ! Each neighboring fault gets weight of -1
        do j = 0,nneighbor-1
            neighbor = smoothing_neighbors(ineighbor+j)
            A(ptr_smooth+i-1,neighbor) = -1.0d0*smoothing_squared
            if (rake_constraint%file.eq.'none'.or.rake_constraint%nfields.eq.2) then
                A(ptr_smooth+nsmooth+i-1,neighbor+fault%nrecords) = -1.0d0*smoothing_squared
            endif
        enddo
    enddo

    if (verbosity.ge.2) then
        write(stderr,'(A)') "load_smoothing says: finished"
    endif

    return
    end

!--------------------------------------------------------------------------------------------------!

    subroutine load_slip_constraints()
    use io, only: stderr, verbosity
    use variable_module, only: displacement, prestress, fault, &
                                disp_components, slip_constraint, rake_constraint, &
                                smoothing
    implicit none
    ! Local variables
    integer :: i, j, n, dn, ndsp, nflt, nsmo
    integer :: cols2delete(ncols), rows2delete(nrows)
    character(len=256) :: fmt

    if (verbosity.ge.2) then
        write(stderr,'(A)') "load_slip_constraints says: starting"
    endif

    if (verbosity.ge.3) then
        write(stderr,'(A)') 'load_arrays says: arrays before modifying topology'
        write(stderr,'("A =")')
        write(fmt,9001) ncols
    9001 format('(1P',I6,'E12.4)')
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

    nflt = fault%nrecords
    ndsp = displacement%nrecords
    nsmo = smoothing%nrecords
    rows2delete = 0
    cols2delete = 0

    ! Move contribution from constrained fault segments to observation side of equation
    ! (do not renumber yet.....will do that later in subroutine)
    do i = 1,nflt
        ! Is the strike-slip (or constant rake) component of this fault constrained?
        if (dabs(slip_constraint%array(i,1)).lt.99998.0d0) then
            ! Move displacements corresponding to this slip component to b vector
            if (displacement%file.ne.'none') then
                do n = 0,len_trim(disp_components)-1
                    dn = n*ndsp
                    do j = 1,ndsp
                        b(j+dn,1) = b(j+dn,1) - A(j+dn,i)*slip_constraint%array(i,1)
                    enddo
                enddo
            endif

            ! Move stresses corresponding to this slip component to b vector
            if (prestress%file.ne.'none') then
                ! Strike-parallel or constant rake-parallel shear tractions
                do j = 1,nflt
                    b(ptr_stress+j-1,1) = b(ptr_stress+j-1,1) &
                                            - A(ptr_stress+j-1,i)*slip_constraint%array(i,1)
                enddo

                ! Dip-parallel shear tractions
                if (rake_constraint%file.eq.'none'.or.rake_constraint%nfields.eq.2) then
                    do j = 1,nflt
                        b(ptr_stress+nflt+j-1,1) = b(ptr_stress+nflt+j-1,1) &
                                            - A(ptr_stress+nflt+j-1,i)*slip_constraint%array(i,1)
                    enddo
                endif
                rows2delete(ptr_stress+i-1) = 1
            endif

            ! Move smoothing corresponding to this slip component to b vector
            if (smoothing%file.ne.'none') then
                do j = 1,nsmo
                    b(ptr_smooth+j-1,1) = b(ptr_smooth+j-1,1) &
                                            - A(ptr_smooth+j-1,i)*slip_constraint%array(i,1)
                enddo
            endif

            ! Indicate this column has been adjusted by changing cols2delete
            cols2delete(i) = 1
        endif

        ! Is the dip-slip component of this fault constrained?
        if (dabs(slip_constraint%array(i,2)).lt.99998.0d0 .and. &
                            (rake_constraint%file.eq.'none'.or.rake_constraint%nfields.eq.2)) then

            ! Move displacements corresponding to this slip component to b vector
            if (displacement%file.ne.'none') then
                do n = 0,len_trim(disp_components)-1
                    dn = n*ndsp
                    do j = 1,ndsp
                        b(j+dn,1) = b(j+dn,1) - A(j+dn,i+nflt)*slip_constraint%array(i,2)
                    enddo
                enddo
            endif

            ! Move stresses corresponding to this slip component to b vector
            if (prestress%file.ne.'none') then
                do j = 1,nflt
                    b(ptr_stress+j-1,1) = b(ptr_stress+j-1,1) &
                                            - A(ptr_stress+j-1,i+nflt)*slip_constraint%array(i,2)
                    b(ptr_stress+nflt+j-1,1) = b(ptr_stress+nflt+j-1,1) &
                                     - A(ptr_stress+nflt+j-1,i+nflt)*slip_constraint%array(i,2)
                enddo
                rows2delete(ptr_stress+i+nflt-1) = 1
            endif

            ! Move smoothing corresponding to this slip component to b vector
            if (smoothing%file.ne.'none') then
                do j = 1,nsmo
                    b(ptr_smooth+nsmo+j-1,1) = b(ptr_smooth+nsmo+j-1,1) &
                            - A(ptr_smooth+nsmo+j-1,i+nflt)*slip_constraint%array(i,2)
                enddo
            endif

            ! Indicate this column has been adjusted by making first entry ridiculous
            cols2delete(i+nflt) = 1
        endif
    enddo

    ! Renumber columns of arrays
    if (verbosity.ge.2) then
        write(stderr,'(A)') 'load_slip_constraints says: renumbering matrix columns'
        write(stderr,'(A,I10)') '    ncols=',ncols
    endif
    ! Start counters for existing columns and number of shifted columns
    j = 0
    n = 0
    do i = 1,ncols
        ! If this slip value has not been fixed, then add it to the modified array
        if (cols2delete(i).eq.0) then
            j = j + 1
            A(:,j) = A(:,i)
        else
            n = n + 1
        endif
    enddo
    ncols = ncols - n
    if (verbosity.ge.2) then
        write(stderr,'(A)') 'load_slip_constraints says: new number of columns:'
        write(stderr,'(A,I10)') '    ncols=',ncols
    endif

    ! Renumber rows of arrays
    if (verbosity.ge.2) then
        write(stderr,'(A)') 'load_slip_constraints says: renumbering matrix rows'
        write(stderr,'(A,I10)') '    nrows=',nrows
    endif
    ! Start counters for existing rows and number of shifted rows
    j = 0
    n = 0
    do i = 1,nrows
        ! If this row is retained, add it to the modified array
        if (rows2delete(i).eq.0) then
            j = j + 1
            A(j,:) = A(i,:)
            b(j,:) = b(i,:)
        else
            n = n + 1
        endif
    enddo
    nrows = nrows - n
    if (verbosity.ge.2) then
        write(stderr,'(A)') 'load_slip_constraints says: new number of rows:'
        write(stderr,'(A,I10)') '    nrows=',nrows
    endif

    if (verbosity.ge.2) then
        write(stderr,'(A)') "load_slip_constraints says: finished"
        write(stderr,*)
    endif

    return
    end

!--------------------------------------------------------------------------------------------------!
!
! !      SUBROUTINE getdamp(soln,gf,obs,nobs,OBSMAX,nflt,FLTMAX,damp,
! !     1                   smooth,smoof,fact,getdmp,sgf,stscon,nsts,cmpnt)
! !      IMPLICIT none
! !      INTEGER getdmp
! !      REAL*8 stscon
! !      INTEGER nsts,cmpnt
! !      CHARACTER*80 smoof
! !      REAL*8 damp,smooth
! !      INTEGER OBSMAX,FLTMAX
! !      REAL*8 obs(OBSMAX,6),prests(FLTMAX,6)
! !      REAL*8 gf(OBSMAX,FLTMAX,6),sgf(FLTMAX,FLTMAX,4)
! !      REAL*8 soln(FLTMAX,2),misfit
! !      REAL*8 p(3)
! !      INTEGER i,j,nobs,nflt
! !      REAL*8 slip,slipmx,ddamp,misfit0,fact
! !      INTEGER vrb
! !      COMMON /VERBOSE/ vrb
! !! Compute misfit with damp = 0
! !      damp = 0.0d0
! !      call lstsqr(soln,gf,obs,nobs,OBSMAX,nflt,FLTMAX,damp,
! !     1            smooth,smoof,sgf,stscon,nsts,prests,cmpnt)
! !      misfit0 = 0.0d0
! !      do 108 i = 1,nobs
! !          p(1) = 0.0d0
! !          p(2) = 0.0d0
! !          p(3) = 0.0d0
! !          do 107 j = 1,nflt
! !              p(1) = p(1) + gf(i,j,1)*soln(j,1)+gf(i,j,4)*soln(j,2)
! !              p(2) = p(2) + gf(i,j,2)*soln(j,1)+gf(i,j,5)*soln(j,2)
! !              p(3) = p(3) + gf(i,j,3)*soln(j,1)+gf(i,j,6)*soln(j,2)
! !  107     continue
! !          misfit0 = misfit0 + (p(1)-obs(i,4))*(p(1)-obs(i,4))
! !     1                    + (p(2)-obs(i,5))*(p(2)-obs(i,5))
! !     2                    + (p(3)-obs(i,6))*(p(3)-obs(i,6))
! !  108 continue
! !      if (vrb.eq.1) then
! !          write(0,1001) misfit0
! !      endif
! !! Search for damping parameter
! !      ddamp = 1.0d0
! !  101 damp = damp + ddamp
! !          call lstsqr(soln,gf,obs,nobs,OBSMAX,nflt,FLTMAX,damp,
! !     1                smooth,smoof,sgf,stscon,nsts,prests,cmpnt)
! !! Compute RMS obs-pre misfit
! !          misfit = 0.0d0
! !          do 102 i = 1,nobs
! !              p(1) = 0.0d0
! !              p(2) = 0.0d0
! !              p(3) = 0.0d0
! !              do 103 j = 1,nflt
! !                  p(1) = p(1) + gf(i,j,1)*soln(j,1)+gf(i,j,4)*soln(j,2)
! !                  p(2) = p(2) + gf(i,j,2)*soln(j,1)+gf(i,j,5)*soln(j,2)
! !                  p(3) = p(3) + gf(i,j,3)*soln(j,1)+gf(i,j,6)*soln(j,2)
! !  103         continue
! !              misfit = misfit + (p(1)-obs(i,4))*(p(1)-obs(i,4))
! !     1                        + (p(2)-obs(i,5))*(p(2)-obs(i,5))
! !     2                        + (p(3)-obs(i,6))*(p(3)-obs(i,6))
! !  102     continue
! !          if (vrb.eq.1) then
! !              write(0,1002) damp,misfit
! !          endif
! !! Get maximum slip in model
! !          slipmx = 0.0d0
! !          do 105 i = 1,nflt
! !              slip = dsqrt(soln(i,1)*soln(i,1)+soln(i,2)*soln(i,2))
! !              if (slip.gt.slipmx) slipmx = slip
! !  105     continue
! !          print *,getdmp,misfit,misfit0,slipmx,fact
! !      if (getdmp.eq.1.and.misfit.le.misfit0*fact) then
! !          goto 101
! !      elseif (getdmp.eq.2.and.slipmx.gt.fact) then
! !          goto 101
! !      elseif (ddamp.gt.1.0d-4) then
! !          damp = damp - ddamp
! !          ddamp = ddamp*1.0d-1
! !          damp = damp - ddamp
! !          goto 101
! !      else
! !          continue
! !      endif
! !      RETURN
! ! 1001 format('    Misfit with no damping is ',1PE14.6)
! ! 1002 format('    Misfit with damping = ',1PE10.3,' is ',1PE14.6)
! !      END
!
!--------------------------------------------------------------------------------------------------!

    subroutine solve_lsqr_dgels()
    !----
    ! Solve generalized least squares problem Ax = b for x using generalized least squares.
    !----
    use io, only : stderr, verbosity
    implicit none
    ! Local variables
    integer :: i, m, n, lda, ldb
    integer :: nrhs
    integer :: lwork, info
    character(len=1) :: trans
    double precision, allocatable :: work(:)
    double precision :: btmp(nrows,1)
    double precision :: alocal(nrows,ncols), blocal(nrows,1)

    if (verbosity.ge.2) then
        write(stderr,'(A)') 'lsqr_solve_dgels says: starting'
    endif

    ! Use local arrays of correct size for inversion
    alocal = a(1:nrows,1:ncols)
    blocal = b(1:nrows,1:1)

    ! Copy displacement vector, b to btmp because it is replaced in dgels
    do i = 1,nrows
        btmp(i,1) = blocal(i,1)
    enddo

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
#ifdef USELAPACK
    call dgels(trans,m,n,nrhs,Alocal,lda,blocal,ldb,work,lwork,info)
#endif
    lwork = int(work(1))
    if (verbosity.ge.2) then
        write(stderr,'(A,I10)') 'lsqr_solve_dgels says: lwork=',lwork
    endif
    deallocate(work)
    allocate(work(lwork))

    ! Solve least squares problem for x
#ifdef USELAPACK
    call dgels(trans,m,n,nrhs,alocal,lda,btmp,ldb,work,lwork,info)
#endif
    deallocate(work)

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
        write(stderr,'(A)') 'lsqr_solve_dgels says: finished'
        write(stderr,*)
    endif

    return
    end subroutine solve_lsqr_dgels

!--------------------------------------------------------------------------------------------------!

    subroutine solve_lsqr_nnls()
    !----
    ! Solve generalized least squares problem Ax = b for x using nonnegative least squares.
    !----
    use io, only : stderr, verbosity
    implicit none
    ! Local variables
    integer :: i, m, n
    integer :: indx(ncols), mode
    double precision :: b1col(nrows), rnorm, x1col(ncols), w(ncols)

    if (verbosity.ge.2) then
        write(stderr,'(A)') 'lsqr_solve_nnls says: starting'
    endif

    m = nrows        ! Number of rows in matrix A
    n = ncols        ! Number of columns in matrix A
    b1col = b(:,1)   ! One-column observation vector

    call nnls(a,m,n,b1col,x1col,rnorm,w,indx,mode)

    ! Put solution into x
    do i = 1,ncols
        x(i,1) = x1col(i)
    enddo

    if (verbosity.ge.3) then
        write(stderr,'("x =")')
        do i = 1,ncols
            write(stderr,'(1PE12.4)') x(i,1)
        enddo
    endif
    if (verbosity.ge.2) then
        write(stderr,'(A)') 'lsqr_solve_nnls says: finished'
        write(stderr,*)
    endif

    return
    end subroutine solve_lsqr_nnls

!--------------------------------------------------------------------------------------------------!

    subroutine solve_dgesv()
    !----
    ! Solve Ax = b for x using LU decomposition and partial pivoting (LAPACK routine GESV).
    !----
    use io, only : stderr, verbosity
    implicit none
    ! Local variables
    integer :: i, n, nrhs, lda, ipiv(nrows), ldb, info
    double precision :: alocal(nrows,ncols), blocal(nrows,1)

    if (verbosity.ge.2) then
        write(stderr,'(A)') 'solve_dgesv says: starting'
    endif
    if (nrows.ne.ncols) then
        call usage('!! Error in subroutine solve_dgesv: nrows not equal to ncols')
    endif

    ! Use local arrays of correct size for inversion
    alocal = a(1:nrows,1:ncols)
    blocal = b(1:nrows,1:1)

    n = nrows        ! Number of linear equations
    nrhs = 1         ! Number of right hand sides
    lda = n          ! Leading dimension of A;  lda >= max(1,n)
    ldb = n          ! Leading dimension of b;  ldb >= max(1,n)

    ! Solve for x (put into blocal)
#ifdef USELAPACK
    call dgesv(n, nrhs, alocal, lda, ipiv, blocal, ldb, info)
#endif
    if (info.ne.0) then
        write(0,*) 'solve_dgesv() says: info returned ',info,' indicating error in dgesv()'
        call usage('exiting program at solve_dgesv()')
    endif

    ! Put actual solution into x
    do i = 1,ncols
        x(i,1) = blocal(i,1)
    enddo

    if (verbosity.ge.3) then
        write(stderr,'("x =")')
        do i = 1,ncols
            write(stderr,'(1PE12.4)') x(i,1)
        enddo
    endif
    if (verbosity.ge.2) then
        write(stderr,'(A)') 'solve_dgesv says: finished'
        write(stderr,*)
    endif

    return
    end subroutine solve_dgesv

!--------------------------------------------------------------------------------------------------!

#ifdef USESUPERLU

    subroutine solve_dgssv()
    !----
    ! Solve Ax = b for x using LU decomposition and partial pivoting (SuperLU routine GSSV).
    !----
    use io, only : stderr, verbosity
    implicit none
    ! Local variables
    integer :: i, j, info, nZero, nNonZero, iopt, nrhs, ldb
    integer(8) :: factors
    double precision :: blocal(nrows,1)
    double precision, allocatable :: matrix_value(:)
    integer, allocatable :: row_index(:), column_ptr(:)

    if (verbosity.ge.2) then
        write(stderr,'(A)') 'solve_dgssv says: starting'
    endif
    if (nrows.ne.ncols) then
        call usage('!! Error in subroutine solve_dgssv: nrows not equal to ncols')
    endif

    blocal = b(1:nrows,1:1)

    ! Format input for dgssv in reduced column
    if (.not.allocated(column_ptr)) then
        allocate(column_ptr(ncols+1))
    endif

    nZero = 0
    nNonZero = 0
    do j = 1,ncols
        column_ptr(j) = nNonZero + 1
        do i = 1,nrows
            if (dabs(A(i,j)).gt.1.0d-10) then
                nNonZero = nNonZero + 1
            else
                nZero = nZero + 1
            endif
        enddo
    enddo
    column_ptr(ncols+1) = nNonZero + 1

    if (.not.allocated(matrix_value)) then
        allocate(matrix_value(nNonZero))
    endif
    if (.not.allocated(row_index)) then
        allocate(row_index(nNonZero))
    endif

    nNonZero = 0
    do j = 1,ncols
        do i = 1,nrows
            if (dabs(A(i,j)).gt.1.0d-10) then
                nNonZero = nNonZero + 1
                matrix_value(nNonZero) = A(i,j)
                row_index(nNonZero) = i
            endif
        enddo
    enddo

    !
    nrhs = 1
    ldb = nrows

#ifdef USELAPACK
    ! Factorize matrix
    iopt = 1
    call c_fortran_dgssv(iopt,nrows,nNonZero,nrhs,matrix_value,row_index,column_ptr,blocal,ldb,&
                       factors,info)

    ! Solve the system with factors
    iopt = 2
    call c_fortran_dgssv(iopt,nrows,nNonZero,nrhs,matrix_value,row_index,column_ptr,blocal,ldb,&
                       factors,info)

    ! Free allocated storage
    iopt = 3
    call c_fortran_dgssv(iopt,nrows,nNonZero,nrhs,matrix_value,row_index,column_ptr,blocal,ldb,&
                       factors,info)
    if (info.ne.0) then
        write(0,*) 'solve_dgssv() says: info returned ',info,' indicating error in dgssv()'
        call usage('exiting program at solve_dgssv()')
    endif
#endif

    ! Put actual solution into x
    do i = 1,ncols
        x(i,1) = blocal(i,1)
    enddo

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

    if (verbosity.ge.2) then
        write(stderr,'("x =")')
        do i = 1,ncols
            write(stderr,'(1PE12.4)') x(i,1)
        enddo
    endif
    if (verbosity.ge.2) then
        write(stderr,'(A)') 'solve_dgssv says: finished'
        write(stderr,*)
    endif

    return
    end subroutine solve_dgssv
#endif
! #endif SUPERLU

!--------------------------------------------------------------------------------------------------!

    subroutine count_matrix_zeros(n)
    implicit none
    integer :: i, j, n
    n = 0
    do i = 1,nrows
        do j = 1,ncols
            if (dabs(A(i,j)).lt.1.0d-10) then
                n = n + 1
            endif
        enddo
    enddo
    return
    end subroutine count_matrix_zeros

end module lsqr_module
