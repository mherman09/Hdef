module lsqr_module
    double precision, allocatable :: A(:,:)
    double precision, allocatable :: b(:,:)
    double precision, allocatable :: x(:,:)
    integer :: nrows
    integer :: ncols
    integer :: ptr_disp
    integer :: ptr_stress
    integer :: ptr_damp
    integer :: ptr_smooth
    character(len=8) :: lsqr_mode
    double precision :: stress_weight

!--------------------------------------------------------------------------------------------------!
contains
!--------------------------------------------------------------------------------------------------!

    subroutine invert_lsqr()
    !----
    ! Invert displacements and pre-stresses for fault slip using a linear least squares approach
    !----
    use io_module, only : stderr, verbosity
    use variable_module, only: fault, slip_constraint, rake_constraint, fault_slip
    implicit none
    ! Local variables
    integer :: i, j

    if (verbosity.ge.2) then
        write(stderr,'(A)') 'invert_lsqr says: starting'
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

    ! Load model matrix, A, and constraint vector, b
    call load_arrays()

    ! Solve generalized least squares problem
    if (.not.allocated(x)) then
        allocate(x(ncols,1))
    endif

    if (lsqr_mode.eq.'gels') then
        call solve_lsqr_dgels()
    elseif (lsqr_mode.eq.'nnls') then
        call solve_lsqr_nnls()
    else
        call print_usage('!! Error: no lsqr_mode named '//trim(lsqr_mode))
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
    if (rake_constraint%file.eq.'none') then
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

    if (verbosity.ge.2) then
        write(stderr,'(A)') 'invert_lsqr says: finished'
        write(stderr,*)
    endif

    return
    end subroutine invert_lsqr

!--------------------------------------------------------------------------------------------------!

    subroutine calc_array_dimensions()
    use io_module, only: stderr, verbosity
    use variable_module, only: displacement, prestress, fault, rake_constraint, &
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
    ptr_stress = 0
    ptr_damp = 0
    ptr_smooth = 0

    ! Row dimensions corresponding to displacement data
    if (displacement%file.ne.'none') then
        nrows = nrows + len_trim(disp_components)*displacement%nrecords
        ptr_disp = 1
    endif

    ! Row dimensions corresponding to pre-stress data
    if (prestress%file.ne.'none') then
        if (ptr_disp.eq.1) then
            ptr_stress = nrows + 1
        else
            ptr_stress = 1
        endif
        nrows = nrows + 2*fault%nrecords
    endif

    ! Column dimensions corresponding to fault slip
    if (rake_constraint%file.eq.'none') then
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
        if (rake_constraint%file.eq.'none') then
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
        write(stderr,'(A,I8)') 'ptr_stress: ', ptr_stress
        write(stderr,'(A,I8)') 'ptr_damp:   ', ptr_damp
        write(stderr,'(A,I8)') 'ptr_smooth: ', ptr_smooth
        write(stderr,*)
    endif

    return
    end subroutine calc_array_dimensions

!--------------------------------------------------------------------------------------------------!

    subroutine load_arrays()
    use io_module, only: stderr, verbosity
    use variable_module, only: displacement, prestress, &
                                slip_constraint, damping_constant, smoothing_constant
    implicit none
    ! Local variables
    integer :: i, j
    character(len=256) :: fmt

    if (verbosity.ge.2) then
        write(stderr,'(A)') 'load_arrays says: starting'
    endif

    ! Load displacement data
    if (displacement%file.ne.'none') then
        call load_displacements()
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

    ! Apply constraints
    if (slip_constraint%file.ne.'none') then
        call load_slip_constraints()
    endif

    if (verbosity.ge.2) then
        write(stderr,'(A)') 'load_arrays says: finished'
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
    endif
    9001 format('(1P',I3,'E12.4)')
    if (verbosity.ge.2) then
        write(stderr,*)
    endif

    return
    end subroutine load_arrays

    !--------------------------------------------------------------------------------------------------!

    subroutine load_displacements()
    use io_module, only: stderr, verbosity
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
            if (rake_constraint%file.eq.'none') then
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

    subroutine load_prestresses()
    use io_module, only: stderr, verbosity
    use variable_module, only : fault, prestress, rake_constraint, gf_stress
    implicit none
    ! Local variables
    integer :: i, j, nflt

    nflt = fault%nrecords

    if (verbosity.ge.2) then
        write(stderr,'(A)') 'load_prestresses says: starting'
    endif

    do i = 1,nflt
        do j = 1,nflt
            if (rake_constraint%file.eq.'none') then
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
    use io_module, only : stderr, verbosity
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
    use io_module, only : stderr, verbosity
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
        if (rake_constraint%file.eq.'none') then
            b(ptr_smooth+nsmooth+i-1,1) = 0.0d0
        endif

        ifault = smoothing%intarray(i,1)
        nneighbor = smoothing%intarray(i,2)
        ineighbor = smoothing%intarray(i,3)

        ! Fault to be smoothed gets weight=nneighbor
        A(ptr_smooth+i-1,ifault) = dble(nneighbor)*smoothing_squared
        if (rake_constraint%file.eq.'none') then
            A(ptr_smooth+nsmooth+i-1,ifault+fault%nrecords) = dble(nneighbor)*smoothing_squared
        endif

        ! Each neighboring fault gets weight of -1
        do j = 0,nneighbor-1
            neighbor = smoothing_neighbors(ineighbor+j)
            A(ptr_smooth+i-1,neighbor) = -1.0d0*smoothing_squared
            if (rake_constraint%file.eq.'none') then
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
    use io_module, only: stderr, verbosity
    use variable_module, only: displacement, prestress, fault, &
                                disp_components, slip_constraint, rake_constraint, &
                                smoothing
    implicit none
    ! Local variables
    integer :: i, j, n, dn, ndsp, nflt, nsmo

    if (verbosity.ge.2) then
        write(stderr,'(A)') "load_slip_constraints says: starting"
    endif

    nflt = fault%nrecords
    ndsp = displacement%nrecords
    nsmo = smoothing%nrecords

    ! Move contribution from constrained fault segments to observation side of equation (do not
    ! renumber yet.....will do that later in subroutine)
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
                do j = 1,nflt
                    b(ptr_stress+j-1,1) = b(ptr_stress+j-1,1) &
                                            - A(ptr_stress+j-1,i)*slip_constraint%array(i,1)
                enddo
            endif

            ! Move smoothing corresponding to this slip component to b vector
            if (smoothing%file.ne.'none') then
                do j = 1,nsmo
                    b(ptr_smooth+j-1,1) = b(ptr_smooth+j-1,1) &
                                            - A(ptr_smooth+j-1,i)*slip_constraint%array(i,1)
                enddo
            endif

            ! Indicate this column has been adjusted by making first entry ridiculous
            A(1,i) = -1.0d10
        endif

        ! Is the dip-slip component of this fault constrained?
        if (dabs(slip_constraint%array(i,2)).lt.99998.0d0 .and. &
                                            rake_constraint%file.eq.'none') then

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
                    b(ptr_stress+nflt+j-1,1) = b(ptr_stress+nflt+j-1,1) &
                                     - A(ptr_stress+nflt+j-1,i)*slip_constraint%array(i,2)
                enddo
            endif

            ! Move smoothing corresponding to this slip component to b vector
            if (smoothing%file.ne.'none') then
                do j = 1,nsmo
                    b(ptr_smooth+nsmo+j-1,1) = b(ptr_smooth+nsmo+j-1,1) &
                            - A(ptr_smooth+nsmo+j-1,i+nflt)*slip_constraint%array(i,2)
                enddo
            endif

            ! Indicate this column has been adjusted by making first entry ridiculous
            A(1,i+nflt) = -1.0d10
        endif
    enddo

    ! Renumber columns of arrays
    if (verbosity.ge.2) then
        write(stderr,'(A)') 'load_slip_constraints says: renumbering matrix columns'
    endif

    ! start a counter for number of shifted columns
    n = 0

    do i = 1,ncols
        ! For each column, check if this slip value has been fixed
        if (A(1,i).lt.-1.0d9) then
            ! If so, shift all columns to the right of selected column down by one entry
            if (i.lt.ncols) then
                do j = i,ncols-1
                    A(:,j) = A(:,j+1)
                enddo
                ! Set the last column to 0 so it will not be double-counted
                A(1,ncols-n) = 0.0d0
            endif
            n = n + 1
        endif
    enddo
    ncols = ncols - n

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
    use io_module, only : stderr, verbosity
    implicit none
    ! Local variables
    integer :: i, m, n, lda, ldb
    integer :: nrhs
    integer :: lwork, info
    character(len=1) :: trans
    double precision :: work(100000)
    double precision :: btmp(nrows,1)

    if (verbosity.ge.2) then
        write(stderr,'(A)') 'lsqr_solve_dgels says: starting'
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
    ! do i = 1,nrows
    !     write(0,*) (A(i,j),j=1,ncols)
    ! enddo

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
    use io_module, only : stderr, verbosity
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

end module lsqr_module
