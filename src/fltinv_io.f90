!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------- FLTINV INPUTS --------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

subroutine read_inputs()
!----
! Read the input files for controlling fltinv and set up control variables and arrays. Run checks
! on inputs to make sure fltinv will operate. (PROBABLY NEEDS MORE/BETTER CHECKS. ADD AS NECESSARY.)
!----

use io, only: stderr, stdout, line_count, fileExists, verbosity
use geom, only: strdip2normal
use elast, only: stress2traction, traction_components, read_halfspace_file
use tri_disloc, only: tri_geometry, tri_geo2cart

use fltinv, only: inversion_mode, &
                  fault, &
                  slip_constraint, &
                  rake_constraint, &
                  euler_file, &
                  npoles, &
                  rigid_pt_array_disp, &
                  rigid_pt_array_los, &
                  pole_array, &
                  gf_euler, &
                  displacement, &
                  disp_components, &
                  input_disp_unit, &
                  los, &
                  prestress, &
                  cov_file, &
                  cov_matrix, &
                  isCovMatrixDiagonal, &
                  gf_model, &
                  gf_disp, &
                  gf_los, &
                  gf_stress, &
                  coord_type, &
                  halfspace_file, &
                  poisson, &
                  lame, &
                  shearmod, &
                  smoothing_file, &
                  nsmooth, &
                  smoothing_pointers, &
                  smoothing_neighbors, &
                  fault_slip, &
                  euler_pole, &
                  read_fltinv_data_file

implicit none

! Local variables
integer :: i, j, ios, ierr, ii, nn, mm, ndisp_dof, nlos_dof, ndof
double precision :: dist, dp1, dp2, cov, sts(3,3), nvec(3), vec(3), pt1(3), pt2(3), pt3(3)
double precision :: geo_pt1(3), geo_pt2(3), geo_pt3(3)
double precision, allocatable :: cov_matrix_temp(:)
character(len=512) :: line
character(len=1) :: nchar, mchar
logical :: foundRigidObs, isErrorMessagePrinted


if (verbosity.ge.1) then
    write(stdout,*) 'read_inputs: starting'
endif


!----
! Check that coordinate type is specified
!----
if (coord_type.ne.'cartesian'.and.coord_type.ne.'geographic') then
    call usage('read_inputs: no coordinate type specified: use -xy or -geo option (usage:input)')
endif


!----
! Check for and read observations
!----
! Observations are required to constrain the fault slip or Euler pole
if (displacement%file.eq.'none'.and.los%file.eq.'none'.and.prestress%file.eq.'none') then
    call usage('read_inputs: no displacement, los, or pre-stress file defined (usage:input)')
endif


! Read the observation files
! The number of columns in each file must be set prior to calling read_fltinv_data_file()

! Three-component displacements
if (displacement%file.ne.'none') then

    displacement%ncols = 6 ! x y z ux uy uz

    call read_fltinv_data_file(displacement,ierr)
    if (ierr.ne.0) then
        call usage('read_inputs: problem reading displacement file (usage:input)')
    elseif (displacement%nrows.eq.0) then
        call usage('read_inputs: displacement file is empty (usage:input)')
    endif

endif

! Line-of-sight displacements
if (los%file.ne.'none') then

    los%ncols = 6 ! x y z ulos az inc

    call read_fltinv_data_file(los,ierr)
    if (ierr.ne.0) then
        call usage('read_inputs: problem reading los file (usage:input)')
    elseif (los%nrows.eq.0) then
        call usage('read_inputs: los file is empty (usage:input)')
    endif
endif

! Pre-stresses on faults
if (prestress%file.ne.'none') then

    prestress%ncols = 6 ! sxx syy szz sxy sxz syz (correspond to fault locations defined in fault file)

    call read_fltinv_data_file(prestress,ierr)
    if (ierr.ne.0) then
        call usage('read_inputs: problem reading pre-stress file (usage:input)')
    elseif (prestress%nrows.eq.0) then
        call usage('read_inputs: pre-stress file is empty (usage:input)')
    endif
endif


if (verbosity.ge.2) then
    write(stdout,*) 'read_inputs: finished reading observation files'
endif


!----
! Check for and read model parameter inputs
!----
! A fault file or Euler pole file is required
if (fault%file.eq.'none'.and.euler_file.eq.'none') then
    call usage('read_inputs: no fault file or Euler pole defined')
endif


! The fault file contains the number of sub-faults (and sub-fault locations/geometries if used)
if (fault%file.ne.'none') then

    ! The number of columns to read depends on whether fltinv is reading in Green's functions or
    ! whether they need to be calculated. Providing ANY Green's functions via input file overrides
    ! the calculation of ALL OTHER Green's functions. Explicitly define this behavior here.

    ! If Green's functions are pre-computed, only need the number of lines from the fault file
    ! so set the number of columns to 1
    if (gf_disp%file.ne.'none'.and.displacement%file.ne.'none') then
        fault%ncols = 1
        gf_model = 'precomputed'
    elseif (gf_los%file.ne.'none'.and.los%file.ne.'none') then
        fault%ncols = 1
        gf_model = 'precomputed'
    elseif (gf_stress%file.ne.'none'.and.prestress%file.ne.'none') then
        fault%ncols = 1
        gf_model = 'precomputed'

    ! Otherwise, Green's functions must be computed, so the fault file needs to have correct format
    elseif (gf_model.eq.'okada_rect') then
        fault%ncols = 7 ! evlo evla evdp str dip wid len
    elseif (gf_model.eq.'okada_pt') then
        fault%ncols = 6 ! evlo evla evdp str dip area
    elseif (gf_model.eq.'triangle') then
        fault%ncols = 9 ! v1x v1y v1z v2x v2y v2z v3x v3y v3z

    ! Need to have some Green's functions for an inversion!!!
    else
        write(stderr,*) 'read_inputs: neither GF model nor precomputed GFs are defined ',&
                        '(or did not recognize GF model "',trim(gf_model),'")'
        write(stderr,*) 'Use -gf:model to select a model:'
        write(stderr,*) '    okada_rect'
        write(stderr,*) '    okada_pt'
        write(stderr,*) '    triangle'
        call usage('Or use options such as -gf:disp_file for precomputed GFs  (usage:gf)')
    endif

    ! Read the fault data
    call read_fltinv_data_file(fault,ierr)
    if (ierr.ne.0) then
        call usage('read_inputs: problem reading fault file (usage:fault)')
    endif

    ! Allocate memory for the output fault slip array
    if (.not.allocated(fault_slip)) then
        allocate(fault_slip(fault%nrows,2),stat=ierr)
        if (ierr.ne.0) then
            call usage('read_inputs: error allocating memory to fault_slip (usage:none)')
        endif
    endif
endif

if (verbosity.ge.2) then
    write(stdout,*) 'read_inputs: finished reading fault and fault slip Greens function files'
endif


! The Euler pole file contains the set of observations points experiencing rigid body rotations
! Rigid body rotations defined by Euler poles
if (euler_file.ne.'none') then

    ! Rigid body rotational velocites correspond to units (unlike displacements, where observation
    ! units=fault slip units). These units must be defined explicitly.
    if (input_disp_unit.eq.'m') then
        call usage('read_inputs: disp units are m (displacement), incompatible with rotations '// &
                   '(usage:input)')
    elseif (input_disp_unit.eq.'mm') then
        call usage('read_inputs: disp units are mm (displacement), incompatible with rotations '// &
                   '(usage:input)')
    elseif (input_disp_unit.eq.'m/s') then
        ! Velocity: all good
    elseif (input_disp_unit.eq.'m/yr') then
        ! Velocity: all good
    elseif (input_disp_unit.eq.'mm/s') then
        ! Velocity: all good
    elseif (input_disp_unit.eq.'mm/yr') then
        ! Velocity: all good
    else
        call usage('read_inputs: no disp units named '//trim(input_disp_unit)//' (usage:input)')
    endif

    ! Check that Euler pole file exists, then open it
    if (.not.fileExists(euler_file)) then
        call usage('read_inputs: no Euler pole file named '//trim(euler_file)//' (usage:input)')
    endif
    open(unit=68,file=euler_file,status='old')
    if (verbosity.ge.2) then
        write(stderr,*) 'read_inputs: reading Euler pole file named ',trim(euler_file)
    endif

    ! Parse the Euler pole file
    ! First line is the number of Euler poles
    read(68,*,iostat=ios) npoles
    if (ios.ne.0) then
        call usage('read_inputs: error reading Euler pole file 1st line: npoles (usage:input)')
    endif
    allocate(euler_pole(npoles,3))

    ! Initialize array recording which points are affected by each rigid body rotation
    if (displacement%file.ne.'none') then
        allocate(rigid_pt_array_disp(displacement%nrows))
        rigid_pt_array_disp = 0
    endif
    if (los%file.ne.'none') then
        allocate(rigid_pt_array_los(los%nrows))
        rigid_pt_array_los = 0
    endif

    ! Initialize array for prior constraints on Euler pole location and angular velocity
    if (inversion_mode.eq.'anneal'.or.inversion_mode.eq.'anneal-psc') then
        allocate(pole_array(npoles,5))
        pole_array = 0.0d0
    endif

    ! Read prior constraints on Euler pole locations: lon lat rad min_rate max_rate
    do i = 1,npoles
        read(68,'(A)',end=4001,iostat=ios) line
        if (allocated(pole_array)) then
            read(line,*,end=4002,err=4002,iostat=ios) pole_array(i,1:5)
            if (pole_array(i,4).ge.pole_array(i,5)) then
                call usage('read_inputs: min pole velocity greater than max pole velocity '//&
                           '(usage:euler)')
            endif
        endif
    enddo
    4001 if (ios.ne.0) then
        call usage('read_inputs: reached end of Euler file before finished reading (usage:none)')
    endif
    4002 if (ios.ne.0) then
        call usage('read_inputs: error parsing Euler pole lon lat radius(km) min_rate(deg/Ma) '//&
                   'max_rate(deg/Ma) from: '//trim(line)//' (usage:none)')
    endif

    ! Read groups of points experiencing rigid rotations
    do
        read(68,'(A)',end=4003,iostat=ios) line
        read(line,*,end=4004,err=4004,iostat=ios) i,nchar,j
        if (nchar.eq.'3') then
            if (allocated(rigid_pt_array_disp)) then
                rigid_pt_array_disp(j) = i
            else
                call usage('read_inputs: specified rigid rotation for 3D motion but array '//&
                           'unallocated (usage:none)')
            endif
        elseif (nchar.eq.'l'.or.nchar.eq.'L') then
            if (allocated(rigid_pt_array_los)) then
                rigid_pt_array_los(j) = i
            else
                call usage('read_inputs: specified rigid rotation for LOS motion but array '//&
                           'unallocated (usage:none)')
            endif
        else
            call usage('read_inputs: no observation type code named '//nchar//'; use "3" or "L"'// &
                       ' (usage:none)')
        endif
    enddo
    4003 ios = 0
    4004 if (ios.ne.0) then
        call usage('read_inputs: error parsing ipole 3|L iobs from: '//trim(line)// &
                   ' (usage:none)')
    endif

    ! Check that some point is constraining the Euler pole rotation
    foundRigidObs = .false.
    if (allocated(rigid_pt_array_disp)) then
        if (maxval(rigid_pt_array_disp).gt.0) then
            foundRigidObs = .true.
        endif
    endif
    if (allocated(rigid_pt_array_los)) then
        if (maxval(rigid_pt_array_los).gt.0) then
            foundRigidObs = .true.
        endif
    endif
    if (.not.foundRigidObs) then
        call usage('read_inputs: user wants Euler pole but no observations have rigid rotations '//&
                   '(usage:euler)')
    endif

    ! Initialize Greens functions for rotations
    allocate(gf_euler(2*displacement%nrows+los%nrows,3*npoles))
    gf_euler = 0.0d0

    ! All done with rigid block file
    close(68)

    if (verbosity.ge.2) then
        write(stdout,*) 'read_inputs: finished reading Euler pole file'
    endif
else
    if (verbosity.ge.2) then
        write(stdout,*) 'read_inputs: no Euler poles being used'
    endif
    npoles = 0
endif


!----
! Check input values
!----
! Input fault depth and dimensions depend on GF model
if (fault%file.ne.'none') then

    ! Rectangular faults in an elastic half-space (okada_rect)
    if (gf_model.eq.'okada_rect') then

        ! Depth is positive down, units in meters
        if (minval(fault%array(:,3)).lt.0.0d0) then
            ! Found negative depth
            write(stderr,*) 'read_inputs: found fault depth of ',minval(fault%array(:,3))
            call usage('Depth is defined positive down for gf_model '//trim(gf_model)// &
                       '(usage:none)')
        elseif (maxval(fault%array(:,3)).lt.1000.0d0) then
            ! All depths are small, using km?
            write(stderr,*) 'read_inputs: all fault depths are less than +1000 meters (using ', &
                            'gf_model ',trim(gf_model),')'
        endif

        ! Fault dimensions are in meters
        if (maxval(fault%array(:,6)).lt.100.0d0.and.maxval(fault%array(:,7)).lt.100.0d0) then
            ! All dimensions are small, using km?
            write(stderr,*) 'read_inputs: all fault dimensions are less than 100 meters ',&
                            ' using gf_model ',trim(gf_model)
        endif

    ! Point source faults in an elastic half-space
    elseif (gf_model.eq.'okada_pt') then

        ! Depth is positive down, units in meters
        if (minval(fault%array(:,3)).lt.0.0d0) then
            ! Found negative depth
            write(stderr,*) 'read_inputs: found fault depth of ',minval(fault%array(:,3))
            call usage('Depth is defined positive down for gf_model '//trim(gf_model)// &
                       '(usage:none)')
        elseif (maxval(fault%array(:,3)).lt.1000.0d0) then
            ! All depths are small, using km?
            write(stderr,*) 'read_inputs: all fault depths are less than +1000 meters (using ', &
                            'gf_model ',trim(gf_model),')'
        endif

        ! Fault areas are in square meters
        if (maxval(fault%array(:,6)).lt.10000.0d0) then
            ! Area is small, using km^2?
            write(stderr,*) 'read_inputs: all fault areas are less than 100x100 square meters ',&
                            'using gf_model ',trim(gf_model)
        endif

    ! Triangular faults in an elastic half-space
    elseif (gf_model.eq.'triangle') then
        if (minval([fault%array(:,3),fault%array(:,6),fault%array(:,9)]).lt.0.0d0) then
            ! Found negative depth
            write(stderr,*) 'read_inputs: found triangle vertex depth of ', &
                            minval([fault%array(:,3),fault%array(:,6),fault%array(:,9)])
            call usage('Depth is defined positive down for gf_model '//trim(gf_model)// &
                       '(usage:none)')
        endif

        if (maxval([fault%array(:,3),fault%array(:,6),fault%array(:,9)]).lt.1000.0d0) then
            ! All depths are small, using km?
            write(stderr,*) 'read_inputs: all fault depths are less than +1000 meters (using ', &
                            'gf_model ',trim(gf_model),')'
        endif

    elseif (gf_model.eq.'precomputed') then
        ! Do nothing...the user is on the hook for getting this right

    else
        call usage('read_inputs: no gf_model named "'//trim(gf_model)//'"'//' (usage:gf)')
    endif
endif


! Check observation coordinates, which depend on coordinate type (Cartesian or geographic)
if (coord_type.eq.'cartesian') then

    ! Calculate distances between all faults and first displacement/los input
    ! In 'cartesian' mode, this value should be typically be 1e1-1e5 meters

    if (displacement%file.ne.'none'.and.gf_model.ne.'precomputed') then
        do i = 1,fault%nrows
            dist = (displacement%array(1,1)-fault%array(i,1))**2 + &
                                                       (displacement%array(1,2)-fault%array(i,2))**2
            if (dsqrt(dist).le.10.0d0) then
                write(stderr,*) 'read_inputs: small fault-displacement distance found'
                write(stderr,*) 'Did you mean to use the -geo flag?'
                exit
            endif
        enddo
    endif

    if (los%file.ne.'none'.and.gf_model.ne.'precomputed') then
        do i = 1,fault%nrows
            dist = (los%array(1,1)-fault%array(i,1))**2+(los%array(1,2)-fault%array(i,2))**2
            if (dsqrt(dist).le.10.0d0) then
                write(stderr,*) 'read_inputs: small fault-los distance found'
                write(stderr,*) 'Did you mean to use the -geo flag?'
                exit
            endif
        enddo
    endif

elseif (coord_type.eq.'geographic') then

    ! I am going to operate under the assumption that if the user specifically indicated to use
    ! geographic mode, then they used geographic coordinates. At some point I will make a check.
    ! ...
    ! ...Okay, FINE, I will not be lazy. Here is your stupid stupidity check.

    if (displacement%file.ne.'none') then
        if (maxval(displacement%array(:,1)).gt.360.0d0) then
            call usage('read_inputs: found displacement longitude greater than 360 (usage:none)')
        elseif (minval(displacement%array(:,1)).lt.-180.0d0) then
            call usage('read_inputs: found displacement longitude less than -180 (usage:none)')
        elseif (maxval(displacement%array(:,2)).gt.90.0d0) then
            call usage('read_inputs: found displacement latitude greater than 90 (usage:none)')
        elseif (minval(displacement%array(:,2)).lt.-90.0d0) then
            call usage('read_inputs: found displacement latitude less than -90 (usage:none)')
        endif
    endif

    if (los%file.ne.'none') then
        if (maxval(los%array(:,1)).gt.360.0d0) then
            call usage('read_inputs: found los longitude greater than 360 (usage:none)')
        elseif (minval(los%array(:,1)).lt.-180.0d0) then
            call usage('read_inputs: found los longitude less than -180 (usage:none)')
        elseif (maxval(los%array(:,2)).gt.90.0d0) then
            call usage('read_inputs: found los latitude greater than 90 (usage:none)')
        elseif (minval(los%array(:,2)).lt.-90.0d0) then
            call usage('read_inputs: found los latitude less than -90 (usage:none)')
        endif
    endif

    if (gf_model.eq.'okada_rect'.or.gf_model.eq.'okada_pt') then
        if (maxval(fault%array(:,1)).gt.360.0d0) then
            call usage('read_inputs: found fault longitude greater than 360 (usage:none)')
        elseif (minval(fault%array(:,1)).lt.-180.0d0) then
            call usage('read_inputs: found fault longitude less than -180 (usage:none)')
        elseif (maxval(fault%array(:,2)).gt.90.0d0) then
            call usage('read_inputs: found fault latitude greater than 90 (usage:none)')
        elseif (minval(fault%array(:,2)).lt.-90.0d0) then
            call usage('read_inputs: found fault latitude less than -90 (usage:none)')
        endif
    elseif (gf_model.eq.'triangle') then
        if (maxval(fault%array(:,1)).gt.360.0d0 .or. &
                maxval(fault%array(:,4)).gt.360.0d0 .or. &
                maxval(fault%array(:,7)).gt.360.0d0) then
            call usage('read_inputs: found fault longitude greater than 360 (usage:none)')
        elseif (minval(fault%array(:,1)).lt.-180.0d0 .or. &
                minval(fault%array(:,4)).lt.-180.0d0 .or. &
                minval(fault%array(:,7)).lt.-180.0d0) then
            call usage('read_inputs: found fault longitude less than -180 (usage:none)')
        elseif (maxval(fault%array(:,2)).gt.90.0d0 .or. &
                maxval(fault%array(:,5)).gt.90.0d0 .or. &
                maxval(fault%array(:,8)).gt.90.0d0) then
            call usage('read_inputs: found fault latitude greater than 90 (usage:none)')
        elseif (minval(fault%array(:,2)).lt.-90.0d0 .or. &
                minval(fault%array(:,5)).lt.-90.0d0 .or. &
                minval(fault%array(:,8)).lt.-90.0d0) then
            call usage('read_inputs: found fault latitude less than -90 (usage:none)')
        endif
    endif

else
    call usage('read_inputs: no coordinate type named "'//trim(coord_type)//'" (usage:input)')
endif


! The number of stresses must be equal to the number of faults
if (prestress%file.ne.'none') then
    if (prestress%nrows.ne.fault%nrows) then
        call usage('read_inputs: the number of pre-stresses is not equal to the number of faults'//&
                   '(usage:none)')
    endif
endif

if (verbosity.ge.2) then
    write(stdout,*) 'read_inputs: finished check on observation and fault inputs; looks okay'
endif


!----
! Convert pre-stresses from tensor components to shear tractions
!----
if (prestress%file.ne.'none') then

    if (verbosity.ge.2) then
        write(stdout,*) 'read_inputs: calculating shear tractions from pre-stresses'
    endif

    do i = 1,fault%nrows

        ! Get fault normal vector
        if (gf_model.eq.'okada_rect'.or.gf_model.eq.'okada_pt') then
            call strdip2normal(fault%array(i,4),fault%array(i,5),nvec)
        elseif (gf_model.eq.'triangle') then
            if (coord_type.eq.'cartesian') then
                call tri_geometry(nvec,vec,vec,fault%array(i,1:3),fault%array(i,4:6), &
                                  fault%array(i,7:9))
            elseif (coord_type.eq.'geographic') then
                ! Triangle points: lon lat dep(m) to x y z
                geo_pt1 = fault%array(i,1:3)
                geo_pt2 = fault%array(i,4:6)
                geo_pt3 = fault%array(i,7:9)
                call tri_geo2cart(pt1,pt2,pt3,geo_pt1,geo_pt2,geo_pt3,'m')
                call tri_geometry(nvec,vec,vec,pt1,pt2,pt3)
            endif
        else
            write(stderr,*) 'read_inputs: no gf_model named "'//trim(gf_model)//'"'
            write(stderr,*) 'Available models:'
            write(stderr,*) '    okada_rect'
            write(stderr,*) '    okada_pt'
            call usage(     '    triangle (usage:input)')
        endif

        ! Compute traction vector
        sts(1,1) = prestress%array(i,1)
        sts(2,2) = prestress%array(i,2)
        sts(3,3) = prestress%array(i,3)
        sts(1,2) = prestress%array(i,4)
        sts(2,1) = prestress%array(i,4)
        sts(1,3) = prestress%array(i,5)
        sts(3,1) = prestress%array(i,5)
        sts(2,3) = prestress%array(i,6)
        sts(3,2) = prestress%array(i,6)
        call stress2traction(sts,nvec,vec)

        ! Project traction vector onto fault strike and up-dip directions
        call traction_components(vec,nvec,dp1,prestress%array(i,1),prestress%array(i,2))
    enddo

    if (verbosity.eq.2) then
        write(stdout,*) 'read_inputs: initial shear tractions computed'
    endif
endif


!----
! Set up Green's functions arrays
!----
if (verbosity.ge.2) then
    write(stdout,*) 'read_inputs: allocating memory to Greens function arrays'
endif

! Three-component displacements
if (displacement%file.ne.'none') then
    ! Assign maximum possible dimensions for three-component displacement Green's functions
    gf_disp%ncols = 2*fault%nrows

    if (gf_disp%file.ne.'none') then
        ! Read pre-computed three-component displacement Green's functions
        call read_fltinv_data_file(gf_disp,ierr)
        if (ierr.ne.0) then
            call usage('read_inputs: problem reading gf_disp file (usage:none)')
        endif

        ! Verify that there are the correct number of rows in the array
        if (gf_disp%nrows .ne. 3*displacement%nrows) then
            call usage('read_inputs: number of lines in three-component displacement GF file '// &
                       'must be 3*ndisplacements (one line per displacement DOF) (usage:none)')
        endif
    else
        ! Displacement Green's functions need to be calculated; allocate memory to array
        gf_disp%nrows = 3*displacement%nrows
        if (allocated(gf_disp%array)) then
            deallocate(gf_disp%array)
        endif
        allocate(gf_disp%array(gf_disp%nrows,gf_disp%ncols),stat=ierr)
        if (ierr.ne.0) then
            call usage('read_inputs: error allocating memory to gf_disp%array (usage:none)')
        endif
    endif
endif

! Line-of-sight displacements
if (los%file.ne.'none') then
    ! Assign maximum possible dimensions for LOS Green's functions
    gf_los%ncols = 2*fault%nrows

    if (gf_los%file.ne.'none') then
        ! Read pre-computed LOS displacement Green's functions
        call read_fltinv_data_file(gf_los,ierr)
        if (ierr.ne.0) then
            call usage('read_inputs: problem reading gf_los file (usage:none)')
        endif

        ! Verify that there are the correct number of rows in the array
        if (gf_los%nrows .ne. los%nrows) then
            call usage('read_inputs: number of lines in LOS displacement GF file '// &
                       'must be ndisplacements (one line per displacement DOF) (usage:none)')
        endif!
    else
        ! Allocate memory to calculate Green's functions
        gf_los%nrows = los%nrows
        if (allocated(gf_los%array)) then
            deallocate(gf_los%array)
        endif
        allocate(gf_los%array(gf_los%nrows,gf_los%ncols),stat=ierr)
        if (ierr.ne.0) then
            call usage('read_inputs: error allocating memory to gf_los%array (usage:none)')
        endif
    endif
endif

! Pre-stresses
if (prestress%file.ne.'none'.or.inversion_mode.eq.'anneal-psc') then
    ! Assign maximum possible dimensions for pre-stress Green's functions
    gf_stress%ncols = 2*fault%nrows

    if (gf_stress%file.ne.'none') then
        ! Read pre-computed displacement Green's functions
        call read_fltinv_data_file(gf_stress,ierr)
        if (ierr.ne.0) then
            call usage('read_inputs: problem reading gf_stress file (usage:none)')
        endif

        ! Verify that there are the correct number of rows in the array
        if (gf_stress%nrows .ne. 2*fault%nrows) then
            call usage('read_inputs: number of lines in pre-stress GF file must be '//&
                       '2*nfaults (one line per fault slip DOF) (usage:none)')
        endif
    else
        ! Allocate memory to calculate Green's functions
        gf_stress%nrows = 2*fault%nrows
        if (allocated(gf_stress%array)) then
            deallocate(gf_stress%array)
        endif
        allocate(gf_stress%array(gf_stress%nrows,gf_stress%ncols),stat=ierr)
        if (ierr.ne.0) then
            call usage('read_inputs: error allocating memory to gf_stress%array (usage:none)')
        endif
    endif
endif

if (verbosity.ge.2) then
    write(stdout,*) 'read_inputs: Greens function arrays allocated'
endif


!----
! Set up arrays for other inversion controls
!----
! Slip magnitude constraints
if (slip_constraint%file.ne.'none') then

    ! The meaning of the slip constraint file depends on the inversion mode
    !     mode=lsqr: fix slip components set to values less than 99999
    !     mode=anneal: range of slip magnitudes to search through
    !     mode=anneal-psc: slip for locked faults (can be 99999 to set a fault to be always unlocked)

    ! Read constraints from file
    slip_constraint%ncols = 2
    call read_fltinv_data_file(slip_constraint,ierr)
    if (ierr.ne.0) then
        write(stderr,*) 'read_inputs: problem reading slip constraint file'
        write(stderr,*) 'Format:'
        call usage(     '    slip_ss slip_ds (usage:fault)')
    endif

    ! There can be 1 value for all faults or 1 value for each fault
    if (slip_constraint%nrows.ne.1.and.slip_constraint%nrows.ne.fault%nrows) then
        call usage('read_inputs: number of slip constraints must be 1 or number of faults '// &
                   '(usage:none)')
    endif

    ! If there is only 1 value, set the full constraint array to that value
    if (slip_constraint%nrows.eq.1) then
        dp1 = slip_constraint%array(1,1)
        dp2 = slip_constraint%array(1,2)
        if (allocated(slip_constraint%array)) then
            deallocate(slip_constraint%array)
        endif
        allocate(slip_constraint%array(fault%nrows,slip_constraint%ncols),stat=ierr)
        if (ierr.ne.0) then
            call usage('read_inputs: error allocating memory to slip_constraint%array (usage:none)')
        endif
        slip_constraint%nrows = fault%nrows
        slip_constraint%array(:,1) = dp1
        slip_constraint%array(:,2) = dp2
    endif
else
    if (inversion_mode.eq.'anneal') then
        call usage('read_inputs: a slip constraint file is required for annealing search '// &
                   '(usage:fault)')
    elseif (inversion_mode.eq.'anneal-psc') then
        call usage('read_inputs: a slip constraint file is required for annealing+'//&
                   'pseudo-coupling search (usage:fault)')
    endif
endif


! Rake constraints
if (rake_constraint%file.ne.'none') then

    ! The meaning of the rake constraint file depends on its format and the inversion mode
    !     mode=lsqr, ncols=1: fix the rake angle of the fault
    !     mode=lsqr, ncols=2: use these rake angles instead of 0 and 90 to compute slip GFs
    !     mode=anneal, ncols=2: range of rake angles to search through
    !     mode=anneal-psc, ncols=1: fix the rake angles of the faults and only allow pseudo-coupling
    !                               slip in this direction

    if (inversion_mode.eq.'lsqr') then
        ! Count the number of columns in the file
        open(unit=33,file=rake_constraint%file,status="old")
        read(33,'(A)') line
        read(line,*,iostat=ios) dp1, dp2
        if (ios.eq.0) then
            rake_constraint%ncols = 2
        else
            rake_constraint%ncols = 1
        endif
        close(33)
    elseif (inversion_mode.eq.'anneal') then
        rake_constraint%ncols = 2
    elseif (inversion_mode.eq.'anneal-psc') then
        rake_constraint%ncols = 1
    else
        write(stderr,*) 'read_inputs: this inversion mode does not use rake constraints'
    endif

    ! Read constraints from file
    call read_fltinv_data_file(rake_constraint,ierr)
    if (ierr.ne.0) then
        call usage('read_inputs: problem reading rake constraint file (usage:fault)')
    endif

    ! There can be 1 value for all faults or 1 value for each fault
    if (rake_constraint%nrows.ne.1.and.rake_constraint%nrows.ne.fault%nrows) then
        call usage('read_inputs: number of rake constraints must be 1 or number of faults '// &
                   '(usage:none)')
    endif

    ! If there is only 1 value, set the full constraint array to that value
    if (rake_constraint%nrows.eq.1) then
        dp1 = rake_constraint%array(1,1)
        if (rake_constraint%ncols.eq.2) then
            dp2 = rake_constraint%array(1,2)
        endif
        if (allocated(rake_constraint%array)) then
            deallocate(rake_constraint%array)
        endif
        allocate(rake_constraint%array(fault%nrows,rake_constraint%ncols),stat=ierr)
        if (ierr.ne.0) then
            call usage('read_inputs: error allocating memory to rake_constraint%array (usage:none)')
        endif
        rake_constraint%nrows = fault%nrows
        rake_constraint%array(:,1) = dp1
        if (rake_constraint%ncols.eq.2) then
            rake_constraint%array(:,2) = dp2
        endif
    endif
else
    if (inversion_mode.eq.'anneal') then
        call usage('read_inputs: a rake constraint file is required for annealing search '// &
                   '(usage:fault)')
    endif

    ! Number of rakes to compute is two, even if constraint file is not provided
    rake_constraint%ncols = 2
endif


! Smoothing
if (smoothing_file.ne.'none') then

    ! Read the smoothing file (ifault nneighbors neighbor1 neighbor2....) into two arrays:
    !     smoothing_pointers(nsmooth,3): ifault nneighbors index_in_neighbor_array
    !     smoothing_neighbors(nneighbors): single column array with neighbors only

    if (.not.fileExists(smoothing_file)) then
        call usage('read_inputs: no smoothing file found named '//trim(smoothing_file)// &
                   '(usage:lsqr)')
    endif
    nsmooth = line_count(smoothing_file)
    if (nsmooth.gt.fault%nrows) then
        call usage('read_inputs: number of faults to smooth is larger than total number of '// &
                   'faults (usage:none)')
    endif

    ! Read smoothing file and calculate pointers to neighbor array
    allocate(smoothing_pointers(nsmooth,3),stat=ierr)
    if (ierr.ne.0) then
        call usage('read_inputs: error allocating memory to smoothing_pointers (usage:none)')
    endif
    open(unit=34,file=smoothing_file,status='old')
    do i = 1,nsmooth
        read(34,*) (smoothing_pointers(i,j),j=1,2)
    enddo

    ! Calculate the index of the location of each fault's neighbors in smoothing_neighbor array
    do i = 1,nsmooth
        if (i.eq.1) then
            smoothing_pointers(i,3) = 1
        else
            smoothing_pointers(i,3) = smoothing_pointers(i-1,3) + smoothing_pointers(i-1,2)
        endif
    enddo

    ! Read smoothing file again and this time store neighbors in array
    allocate(smoothing_neighbors(smoothing_pointers(nsmooth,3)+smoothing_pointers(nsmooth,2)),&
             stat=ierr)
    if (ierr.ne.0) then
        call usage('read_inputs: error allocating memory to smoothing_neighbors (usage:none)')
    endif
    rewind(34)
    do i = 1,nsmooth
        read(34,'(A)') line
        nn = smoothing_pointers(i,2)
        ii = smoothing_pointers(i,3)
        read(line,*,iostat=ios) j,j,(smoothing_neighbors(ii+j-1),j=1,nn)
        if (ios.ne.0) then
            write(stderr,*) 'read_inputs: error parsing smoothing file line ',i
            write(stderr,*) 'Offending line: ',trim(line)
            call usage('Expected format: iflt nnbr nbr1 nbr2...nbrn (usage:none)')
        endif
    enddo
    close(34)

    if (verbosity.ge.2) then
        write(stdout,*) 'read_inputs: finished reading smoothing file "',trim(smoothing_file),'"'
    endif
    if (verbosity.ge.3) then
        write(stdout,*) 'Number of faults to smooth: ',nsmooth
        write(stdout,*) 'Smoothing pointer and neighbor arrays:'
        do i = 1,nsmooth
            write(stdout,*) (smoothing_pointers(i,j),j=1,3)
            do j = 1+smoothing_pointers(i,3),smoothing_pointers(i,2)+smoothing_pointers(i,3)
                write(stdout,*) smoothing_neighbors(j)
            enddo
        enddo
    endif
    if (verbosity.ge.2) then
        write(stdout,*)
    endif
endif


! Read half-space data or use default values (to calculate Green's functions)
if (gf_model.eq.'okada_rect'.or.gf_model.eq.'okada_pt'.or.gf_model.eq.'triangle') then
    poisson = 0.25d0
    lame = 40.0d9
    shearmod = 40.0d9
    if (halfspace_file.ne.'none') then
        call read_halfspace_file(halfspace_file,poisson,shearmod,lame)
    endif

    if (verbosity.ge.2) then
        write(stdout,*) 'Poissons ratio:  ',poisson
        write(stdout,*) 'Shear modulus:   ',shearmod
        write(stdout,*) 'Lames parameter: ',lame
    endif
endif


! Covariance of three-component and/or line-of-sight observations
if (displacement%file.ne.'none'.or.los%file.ne.'none') then

    ! Covariance matrix has nrows and ncols equal to total displacement degrees of freedom
    ndisp_dof = len_trim(disp_components)*displacement%nrows
    nlos_dof = los%nrows
    ndof = ndisp_dof + nlos_dof
    allocate(cov_matrix(ndof,ndof),stat=ierr)
    if (ierr.ne.0) then
        call usage('read_inputs: error allocating memory to cov_matrix (usage:none)')
    endif

    ! Initialize covariance matrix as the identity matrix
    cov_matrix = 0.0d0
    do i = 1,ndisp_dof+nlos_dof
        cov_matrix(i,i) = 1.0d0
    enddo

    ! Read in the covariance values
    isErrorMessagePrinted = .false.
    if (cov_file.ne.'none') then
        open(unit=72,file=cov_file,status='old')
        do
            read(72,'(A)',iostat=ios) line
            if (ios.ne.0) then
                exit
            endif

            ! Covariance file format: i_obs j_obs i_comp j_comp cov
            !     i_obs j_obs: indices corresponding to three-component or LOS displacement file
            !     i_comp j_comp: the displacement component (1, 2, or 3) or line-of-sight (L)
            !     cov: the covariance (m^2)
            read(line,*,iostat=ios) i,j,nchar,mchar,cov
            if (ios.ne.0) then
                write(stderr,*) 'read_inputs: error parsing covariance matrix file'
                write(stderr,*) 'Offending line: ',trim(line)
                call usage('Expected format: i_obs j_obs i_comp j_comp cov (usage:none)')
            endif

            if (nchar.eq.'L'.or.nchar.eq.'l') then
                nn = ndisp_dof
                mm = ndisp_dof
            else
                if (index(disp_components,nchar).le.0) then
                    if (.not.isErrorMessagePrinted) then
                        write(stderr,*) 'read_inputs: found a covariance input not used in inversion'
                        write(stderr,*) 'Line: ',trim(line)
                        write(stderr,*) 'First component (',nchar,') is not used in inversion'
                        write(stderr,*) 'It is not being loaded into the covariance matrix.'
                        isErrorMessagePrinted = .true.
                    endif
                    cycle
                endif
                if (index(disp_components,mchar).le.0) then
                    if (.not.isErrorMessagePrinted) then
                        write(stderr,*) 'read_inputs: found a covariance input not used in inversion'
                        write(stderr,*) 'Line: ',trim(line)
                        write(stderr,*) 'Second component (',mchar,') is not used in inversion'
                        write(stderr,*) 'It is not being loaded into the covariance matrix.'
                        isErrorMessagePrinted = .true.
                    endif
                    cycle
                endif
                nn = (index(disp_components,nchar)-1)*displacement%nrows
                mm = (index(disp_components,mchar)-1)*displacement%nrows
            endif
            cov_matrix(i+nn,j+mm) = cov
            cov_matrix(j+mm,i+nn) = cov
        enddo
        close(72)
        if (verbosity.ge.2) then
            write(stdout,*) 'read_inputs: finished reading covariance file "',trim(cov_file),'"'
        endif
    else
        if (verbosity.ge.2) then
            write(stdout,*) 'read_inputs: data covariance matrix is identity matrix'
        endif
    endif


    ! If covariance matrix is diagonal, then only store diagonal elements
    ! Check if matrix is diagonal
    isCovMatrixDiagonal = .true.
    do i = 1,ndof
        do j = 1,ndof
            ! If off-diagonal term is non-zero, set flag to false and exit
            if (i.ne.j) then
                if (abs(cov_matrix(i,j)).gt.1.0d-10) then
                    isCovMatrixDiagonal = .false.
                    exit
                endif
            endif
        enddo
    enddo

    ! If diagonal, resize covariance array
    if (isCovMatrixDiagonal) then
        ! Store covariance matrix in a temporary array of length ndof
        allocate(cov_matrix_temp(ndof))
        do i = 1,ndof
            cov_matrix_temp(i) = cov_matrix(i,i)
        enddo

        ! Re-size the matrix
        deallocate(cov_matrix)
        allocate(cov_matrix(ndof,1))

        ! Update the covariance matrix with stored values and get rid of the temporary array
        cov_matrix(:,1) = cov_matrix_temp
        deallocate(cov_matrix_temp)

        if (verbosity.ge.1) then
            write(stdout,*) 'read_inputs: covariance matrix is diagonal'
        endif

    else
        if (verbosity.ge.1) then
            write(stdout,*) 'read_inputs: covariance matrix is not diagonal'
        endif
    endif


    if (verbosity.ge.3) then
        write(stdout,*) 'Covariance matrix:'
        do i = 1,ndof
            write(stdout,*) cov_matrix(i,:)
        enddo
    endif
endif


if (verbosity.ge.1) then
    write(stdout,*) 'read_inputs: finished'
    write(stdout,*)
endif

return
end subroutine read_inputs






!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------- FLTINV OUTPUTS -------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

subroutine write_solution()

use trig, only: d2r
use io, only: stdout, stderr, verbosity
use earth, only: pole_xyz2geo
use misfit, only: misfit_rms

use fltinv, only: output_file, &
                  inversion_mode, &
                  fault, &
                  rake_constraint, &
                  euler_file, &
                  euler_output_file, &
                  npoles, &
                  displacement, &
                  disp_components, &
                  disp_misfit_file, &
                  gf_disp, &
                  fault_slip, &
                  euler_pole

implicit none

! Local variables
integer :: i, ierr, ounit, icmp, idsp, iflt, ndsp, nflt, nobs
double precision :: slip_mag, rms, geo_pole(3)
double precision, allocatable :: pre(:), obs(:)


if (verbosity.ge.1) then
    write(stdout,*) 'write_solution: starting'
endif


! Print RMS misfit to three-component displacements
if (disp_misfit_file.ne.'none') then

    ! Open misfit file
    open(unit=81,file=disp_misfit_file,status='unknown')

    ndsp = displacement%nrows
    nobs = len_trim(disp_components)*displacement%nrows
    nflt = fault%nrows

    ! Load predicted three-component displacements
    if (.not.allocated(pre)) then
        allocate(pre(ndsp*len_trim(disp_components)),stat=ierr)
        if (ierr.ne.0) then
            call usage('write_solution: error allocating memory to pre array (usage:none)')
        endif
    endif
    pre = 0.0d0

    do i = 1,len_trim(disp_components)
        read(disp_components(i:i),*) icmp
        do idsp = 1,ndsp
            do iflt = 1,nflt
                pre((i-1)*ndsp+idsp) = pre((i-1)*ndsp+idsp) + &
                                       gf_disp%array((icmp-1)*ndsp+idsp,     iflt)*fault_slip(iflt,1)
                pre((i-1)*ndsp+idsp) = pre((i-1)*ndsp+idsp) + &
                                       gf_disp%array((icmp-1)*ndsp+idsp,nflt+iflt)*fault_slip(iflt,2)
            enddo
        enddo
    enddo

    ! Load observed three-component displacements
    if (.not.allocated(obs)) then
        allocate(obs(nobs),stat=ierr)
        if (ierr.ne.0) then
            call usage('write_solution: error allocating memory to obs array (usage:none)')
        endif
    endif
    obs = 0.0d0

    do i = 1,len_trim(disp_components)
        read(disp_components(i:i),*) icmp
        obs((i-1)*ndsp+1:i*ndsp) = displacement%array(1:ndsp,3+icmp)
    enddo

    call misfit_rms(obs,pre,nobs,rms)
    write(81,*) rms

    close(81)

endif


! ! Line-of-sight RMS misfit
! if (los_misfit_file.ne.'none') then
!     if (verbosity.ge.1) then
!         write(stderr,'(A)') 'write_solution says: writing LOS RMS misfit to '//trim(los_misfit_file)
!     endif
!
!     write(0,*) 'write_solution: opening misfit file'
!     open(unit=81,file=los_misfit_file,status='unknown')
!
!     ! If rake is constrained, make an array of the correct size to use with misfit function
!     if (inversion_mode.eq.'lsqr'.and.rake_constraint%file.ne.'none' &
!                                                            .and.rake_constraint%nfields.eq.1) then
!         do i = 1,fault%nrecords
!             tmp_slip_array(i,1) = fault_slip(i,1) ! Green's functions already calculated for this rake
!             tmp_slip_array(i,2) = 0.0d0
!         enddo
!         write(81,*) los_misfit_l2norm(tmp_slip_array)/dsqrt(dble(fault%nrecords))
!
!     ! Otherwise, just use the misfit function directly
!     else
!         write(81,*) los_misfit_l2norm(fault_slip)/dsqrt(dble(fault%nrecords))
!     endif
!
!     write(0,*) 'write_solution: closing misfit file'
!     close(81)
! endif


! Print fault slip solution
if (fault%file.ne.'none') then
    if (verbosity.ge.1) then
        write(stdout,*) 'write_solution: writing slip solution to '//trim(output_file)
    endif
    if (output_file.eq.'stdout') then
        ounit = stdout
    else
        ounit = 99
        open (unit=ounit,file=output_file,status='unknown')
    endif

    do i = 1,fault%nrows
        slip_mag = dsqrt(fault_slip(i,1)*fault_slip(i,1)+fault_slip(i,2)*fault_slip(i,2))

        if (inversion_mode.eq.'lsqr') then
            if (rake_constraint%ncols.eq.1) then
                write(ounit,5001) fault_slip(i,1)*cos(rake_constraint%array(i,1)*d2r), &
                                  fault_slip(i,1)*sin(rake_constraint%array(i,1)*d2r)
            elseif (rake_constraint%file.ne.'none') then
                write(ounit,5001) fault_slip(i,1)*cos(rake_constraint%array(i,1)*d2r) + &
                                  fault_slip(i,2)*cos(rake_constraint%array(i,2)*d2r), &
                                  fault_slip(i,1)*sin(rake_constraint%array(i,1)*d2r) + &
                                  fault_slip(i,2)*sin(rake_constraint%array(i,2)*d2r)
            else
                write(ounit,5001) fault_slip(i,1:2)
            endif

        elseif (inversion_mode.eq.'anneal') then
            write(ounit,5001) fault_slip(i,1:2)

        elseif (inversion_mode.eq.'anneal-psc') then
            write(ounit,5001) fault_slip(i,1:2)

        else
            write(stderr,*) 'write_solution: frankly, I do not know how you got this far using an '//&
                             'inversion mode that does not seem to exist...'
            write(stderr,*) 'Available inversion modes:'
            write(stderr,*) '    lsqr'
            write(stderr,*) '    anneal'
            call usage(     '    anneal-psc (usage:general)')
        endif
    enddo
    5001 format(1P2E16.8)
    close(ounit)
endif

! Print rigid body rotation solution
if (euler_file.ne.'none') then
    if (verbosity.ge.1) then
        write(stdout,*) 'write_solution: writing Euler pole solution to '//trim(euler_output_file)
    endif
    if (euler_output_file.eq.'stdout') then
        ounit = stdout
    else
        ounit = 99
        open (unit=ounit,file=euler_output_file,status='unknown')
    endif
    do i = 1,npoles
        call pole_xyz2geo(euler_pole(i,1),euler_pole(i,2),euler_pole(i,3), &
                          geo_pole(1),geo_pole(2),geo_pole(3), 'sphere')
        write(ounit,5002) geo_pole(1:3)
    enddo
    5002 format(3F16.4)
    close(ounit)
endif


if (verbosity.ge.1) then
    write(stdout,*) 'write_solution: finished'
endif

return
end subroutine write_solution





!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!----------------------------------- PARSE COMMAND LINE -------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

subroutine gcmdln()

use io, only: stdout, verbosity, debug

use fltinv, only: output_file, &
                  inversion_mode, &
                  fault, &
                  slip_constraint, &
                  rake_constraint, &
                  euler_file, &
                  euler_output_file, &
                  displacement, &
                  disp_components, &
                  disp_misfit_file, &
                  input_disp_unit, &
                  los, &
                  prestress, &
                  cov_file, &
                  gf_model, &
                  gf_disp, &
                  gf_los, &
                  gf_stress, &
                  coord_type, &
                  halfspace_file, &
                  damping_constant, &
                  smoothing_constant, &
                  smoothing_file, &
                  lsqr_mode, &
                  anneal_init_mode, &
                  anneal_init_file, &
                  anneal_step_file, &
                  max_iteration, &
                  reset_iteration, &
                  temp_start, &
                  temp_minimum, &
                  cooling_factor, &
                  anneal_log_file, &
                  anneal_seed, &
                  modelUncertainty, &
                  min_flip, &
                  max_flip, &
                  init_fltinv_data

implicit none

! Local variables
integer :: i, narg, ios
character(len=256) :: tag
! integer :: char_index


! Initialize variables
output_file = ''
inversion_mode = ''

call init_fltinv_data(fault)
call init_fltinv_data(slip_constraint)
call init_fltinv_data(rake_constraint)

euler_file = 'none'
euler_output_file = ''

call init_fltinv_data(displacement)
disp_components = '123'
call init_fltinv_data(los)
call init_fltinv_data(prestress)
cov_file = 'none'
input_disp_unit = 'm'

gf_model = 'none'
call init_fltinv_data(gf_disp)
call init_fltinv_data(gf_los)
call init_fltinv_data(gf_stress)
coord_type = ''

halfspace_file = 'none'

damping_constant = -1.0d0
smoothing_constant = -1.0d0
smoothing_file = 'none'

lsqr_mode = 'gels'

anneal_init_mode = ''
anneal_init_file = 'none'
anneal_step_file = 'none'
max_iteration = 1000
reset_iteration = 1000000
temp_start = 2.0d0
temp_minimum = 0.00d0
cooling_factor = 0.98d0
anneal_log_file = ''
anneal_seed = 0
modelUncertainty = .false.

min_flip = 1
max_flip = 10
! pl2u = 0.01d0
! pu2l = 0.01d0

verbosity = 0
debug = .false.

disp_misfit_file = 'none'
! los_misfit_file = 'none'
!
! stress_weight = 1.0d-9
! los_weight = 1.0d0
! sts_dist = 1.0d10
! prob_lock2unlock = 0.25d0
! prob_unlock2lock = 0.10d0
! mcmc_iteration = 0
! mcmc_log_file = 'none'


narg = command_argument_count()
if (narg.eq.0) then
    call usage('')
endif

i = 1
do while (i.le.narg)

    call get_command_argument(i,tag)

    ! General options
    if (trim(tag).eq.'-o') then
        i = i + 1
        call get_command_argument(i,output_file)
    elseif (trim(tag).eq.'-mode') then
        i = i + 1
        call get_command_argument(i,inversion_mode)
!         if (trim(inversion_mode).eq.'anneal-psc') then
!             anneal_init_mode = 'unlocked'
!         endif

    ! Fault options
    elseif (trim(tag).eq.'-flt') then
        i = i + 1
        call get_command_argument(i,fault%file)
    elseif (trim(tag).eq.'-flt:slip'.or.trim(tag).eq.'-flt:slip_constraint') then
        i = i + 1
        call get_command_argument(i,slip_constraint%file)
    elseif (trim(tag).eq.'-flt:rake'.or.trim(tag).eq.'-flt:rake_constraint') then
        i = i + 1
        call get_command_argument(i,rake_constraint%file)

    ! Euler setup file
    elseif (trim(tag).eq.'-euler') then
        i = i + 1
        call get_command_argument(i,euler_file)
        i = i + 1
        call get_command_argument(i,euler_output_file)

    ! Input options
    elseif (trim(tag).eq.'-disp') then
        i = i + 1
        call get_command_argument(i,displacement%file)
    elseif (trim(tag).eq.'-disp:components') then
        i = i + 1
        call get_command_argument(i,disp_components)
    elseif (trim(tag).eq.'-los') then
        i = i + 1
        call get_command_argument(i,los%file)
    elseif (trim(tag).eq.'-prests'.or.tag.eq.'-prestress') then
        i = i + 1
        call get_command_argument(i,prestress%file)
    elseif (trim(tag).eq.'-cov') then
        i = i + 1
        call get_command_argument(i,cov_file)
    elseif (trim(tag).eq.'-disp:unit') then
        i = i + 1
        call get_command_argument(i,input_disp_unit)
!     elseif (trim(tag).eq.'-los:weight') then
!         i = i + 1
!         call get_command_argument(i,tag)
!         read(tag,*) los_weight
!     elseif (trim(tag).eq.'-los:misfit') then
!         i = i + 1
!         call get_command_argument(i,los_misfit_file)
!     elseif (trim(tag).eq.'-prests:weight') then
!         i = i + 1
!         call get_command_argument(i,tag)
!         read(tag,*) stress_weight
!     elseif (trim(tag).eq.'-prests:dist_threshold') then
!         i = i + 1
!         call get_command_argument(i,tag)
!         read(tag,*) sts_dist

    ! Output options
    elseif (trim(tag).eq.'-disp:misfit') then
        i = i + 1
        call get_command_argument(i,disp_misfit_file)

    ! Green's functions options
    elseif (trim(tag).eq.'-gf:model') then
        i = i + 1
        call get_command_argument(i,gf_model)
    elseif (trim(tag).eq.'-gf:disp') then
        i = i + 1
        call get_command_argument(i,gf_disp%file)
    elseif (trim(tag).eq.'-gf:los') then
        i = i + 1
        call get_command_argument(i,gf_los%file)
    elseif (trim(tag).eq.'-gf:stress') then
        i = i + 1
        call get_command_argument(i,gf_stress%file)
    elseif (trim(tag).eq.'-xy'.or.trim(tag).eq.'-cartesian') then
        coord_type = 'cartesian'
    elseif (trim(tag).eq.'-geo'.or.trim(tag).eq.'-geographic') then
        coord_type = 'geographic'

    ! Half-space options
    elseif (trim(tag).eq.'-haf') then
        i = i + 1
        call get_command_argument(i,halfspace_file)

    ! Regularization options
    elseif (trim(tag).eq.'-damp') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) damping_constant
    elseif (trim(tag).eq.'-smooth') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) smoothing_constant
        i = i + 1
        call get_command_argument(i,smoothing_file)

    ! Least squares options
    elseif (trim(tag).eq.'-lsqr:mode') then
        i = i + 1
        call get_command_argument(i,lsqr_mode)

    ! Simulated annealing options
    elseif (trim(tag).eq.'-anneal:init_mode') then
        i = i + 1
        call get_command_argument(i,anneal_init_mode)
        if (trim(anneal_init_mode).eq.'user') then
            i = i + 1
            call get_command_argument(i,anneal_init_file)
        endif
    elseif (trim(tag).eq.'-anneal:step') then
        i = i + 1
        call get_command_argument(i,anneal_step_file)
    elseif (trim(tag).eq.'-anneal:max_it') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) max_iteration
    elseif (trim(tag).eq.'-anneal:reset_it') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) reset_iteration
    elseif (trim(tag).eq.'-anneal:temp_start') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) temp_start
    elseif (trim(tag).eq.'-anneal:temp_min') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) temp_minimum
    elseif (trim(tag).eq.'-anneal:cool') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) cooling_factor
    elseif (trim(tag).eq.'-anneal:log_file') then
        i = i + 1
        call get_command_argument(i,anneal_log_file)
    elseif (trim(tag).eq.'-anneal:seed') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) anneal_seed
    elseif (tag.eq.'-anneal:model_uncertainty') then
        modelUncertainty = .true.
    elseif (trim(tag).eq.'-anneal-psc:min_flip') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) min_flip
    elseif (trim(tag).eq.'-anneal-psc:max_flip') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) max_flip
    ! elseif (trim(tag).eq.'-anneal-psc:pl2u') then
    !     i = i + 1
    !     call get_command_argument(i,tag)
    !     read(tag,*) pl2u
    ! elseif (trim(tag).eq.'-anneal-psc:pu2l') then
    !     i = i + 1
    !     call get_command_argument(i,tag)
    !     read(tag,*) pu2l
!     elseif (trim(tag).eq.'-anneal:mcmc') then
!         i = i + 1
!         call get_command_argument(i,tag)
!         read(tag,*) mcmc_iteration
!         i = i + 1
!         call get_command_argument(i,mcmc_log_file)
!     elseif (trim(tag).eq.'-anneal:p_unlock2lock' .or. &
!                                 trim(tag).eq.'-anneal:prob_unlock2lock') then
!         i = i + 1
!         call get_command_argument(i,tag)
!         read(tag,*) prob_unlock2lock

    ! Miscellaneous options
    elseif (trim(tag).eq.'-v') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*,iostat=ios) verbosity
        if (ios.ne.0) then
            verbosity = 1
        endif
    elseif (trim(tag).eq.'-debug') then
        debug = .true.

    else
        call usage('fltinv: No option '//trim(tag))
    endif

    i = i + 1
enddo

! Parsed command line values
if (verbosity.ge.2) then
    write(stdout,'(" Parsed command line inputs")')
    write(stdout,'(" output_file:            ",A)') trim(output_file)
    write(stdout,'(" inversion_mode:         ",A)') trim(inversion_mode)
    write(stdout,*)
    write(stdout,'(" fault%file:             ",A)') trim(fault%file)
    write(stdout,'(" slip_constraint%file:   ",A)') trim(slip_constraint%file)
    write(stdout,'(" rake_constraint%file:   ",A)') trim(rake_constraint%file)
    write(stdout,*)
    write(stdout,'(" euler_file:             ",A)') trim(euler_file)
    write(stdout,*)
    write(stdout,'(" displacement%file:      ",A)') trim(displacement%file)
    write(stdout,'(" disp_components:        ",A)') trim(disp_components)
    write(stdout,'(" disp_misfit_file:       ",A)') trim(disp_misfit_file)
    write(stdout,'(" los%file:               ",A)') trim(los%file)
    write(stdout,'(" prestress%file:         ",A)') trim(prestress%file)
    write(stdout,'(" cov_file:               ",A)') trim(cov_file)
    write(stdout,'(" input_disp_unit:        ",A)') trim(input_disp_unit)
    write(stdout,*)
    write(stdout,'(" gf_model:               ",A)') trim(gf_model)
    write(stdout,'(" gf_disp%file:           ",A)') trim(gf_disp%file)
    write(stdout,'(" gf_los%file:            ",A)') trim(gf_los%file)
    write(stdout,'(" gf_stress%file:         ",A)') trim(gf_stress%file)
    write(stdout,'(" coord_type:             ",A)') trim(coord_type)
    write(stdout,*)
    write(stdout,'(" halfspace_file:         ",A)') trim(halfspace_file)
    write(stdout,*)
!     write(stdout,'(" stress_weight:          ",1PE14.6)') stress_weight
!     write(stdout,'(" sts_dist:               ",1PE14.6)') sts_dist
!     write(stdout,'(" los_weight:             ",1PE14.6)') los_weight
!     write(stdout,'(" los_misfit_file:        ",A)') trim(los_misfit_file)
!     write(stdout,*)
    write(stdout,'(" damping_constant:       ",1PE14.6)') damping_constant
    write(stdout,'(" smoothing_constant:     ",1PE14.6)') smoothing_constant
    write(stdout,'(" smoothing_file:         ",A)') trim(smoothing_file)
    write(stdout,*)
    write(stdout,'(" lsqr_mode:              ",A)') trim(lsqr_mode)
    write(stdout,*)
    write(stdout,'(" anneal_init_mode:       ",A)') trim(anneal_init_mode)
    write(stdout,'(" anneal_init_file:       ",A)') trim(anneal_init_file)
    write(stdout,'(" anneal_step_file:       ",A)') trim(anneal_step_file)
    write(stdout,'(" max_iteration:          ",I14)') max_iteration
    write(stdout,'(" reset_iteration:        ",I14)') reset_iteration
    write(stdout,'(" temp_start:             ",1PE14.6)') temp_start
    write(stdout,'(" temp_minimum:           ",1PE14.6)') temp_minimum
    write(stdout,'(" cooling_factor:         ",1PE14.6)') cooling_factor
    write(stdout,'(" anneal_log_file:        ",A)') trim(anneal_log_file)
    write(stdout,'(" anneal_seed:            ",I14)') anneal_seed
    write(stdout,'(" modelUncertainty:       ",A)') modelUncertainty
    write(stdout,'(" min_flip:               ",I14)') min_flip
    write(stdout,'(" max_flip:               ",I14)') max_flip
    ! write(stdout,'(" pl2u:                   ",1PE14.6)') pl2u
    ! write(stdout,'(" pu2l:                   ",1PE14.6)') pu2l
!     write(stdout,'(" mcmc_iteration:          ",A)') mcmc_iteration
!     write(stdout,'(" mcmc_log_file:          ",A)') mcmc_log_file
!     write(stdout,'(" anneal_verbosity:       ",I14)') anneal_verbosity
    write(stdout,*)
endif

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine usage(string)

use io, only: stderr

implicit none

character(len=*) :: string
character(len=8) :: info
integer :: i

info = 'all'
i = len_trim(string)+1

if (string.ne.'') then
    if (index(string,'(usage:general)').ne.0) then
        i = index(string,'(usage:general)')
        info = 'general'
    elseif (index(string,'(usage:fault)').ne.0) then
        i = index(string,'(usage:fault)')
        info = 'fault'
    elseif (index(string,'(usage:euler)').ne.0) then
        i = index(string,'(usage:euler)')
        info = 'euler'
    elseif (index(string,'(usage:input)').ne.0) then
        i = index(string,'(usage:input)')
        info = 'input'
    elseif (index(string,'(usage:output)').ne.0) then
        i = index(string,'(usage:output)')
        info = 'output'
    elseif (index(string,'(usage:gf)').ne.0) then
        i = index(string,'(usage:gf)')
        info = 'gf'
    elseif (index(string,'(usage:hafspc)').ne.0) then
        i = index(string,'(usage:hafspc)')
        info = 'hafspc'
    elseif (index(string,'(usage:lsqr)').ne.0) then
        i = index(string,'(usage:lsqr)')
        info = 'lsqr'
    elseif (index(string,'(usage:anneal)').ne.0) then
        i = index(string,'(usage:anneal)')
        info = 'anneal'
    elseif (index(string,'(usage:none)').ne.0) then
        i = index(string,'(usage:none)')
        info = 'none'
    endif

    write(stderr,*) string(1:i-1)
    write(stderr,*)
    if (info.eq.'none') then
        write(stderr,*) 'See fltinv man page for details'
        write(stderr,*)
        call error_exit(1)
    endif
endif

write(stderr,*) 'Usage: fltinv ...options...'
write(stderr,*)
if (info.eq.'all'.or.info.eq.'general') then
    write(stderr,*) 'General Options'
    write(stderr,*) '-o OUTPUT_FILE               Output fault slip file'
    write(stderr,*) '-mode INVERSION_MODE         Inversion mode'
    write(stderr,*)
endif
if (info.eq.'all'.or.info.eq.'fault') then
    write(stderr,*) 'Fault Options'
    write(stderr,*) '-flt FAULT_FILE              Input fault locations and geometries'
    write(stderr,*) '-flt:slip SLIP_FILE          Slip magnitude constraints'
    write(stderr,*) '-flt:rake RAKE_FILE          Rake angle constraints'
    write(stderr,*)
endif
if (info.eq.'all'.or.info.eq.'euler') then
    write(stderr,*) 'Euler Pole Options'
    write(stderr,*) '-euler EULER_FILE POLE_FILE  Include rigid rotation in the fit [UNDER DEVELOPMENT]'
    write(stderr,*)
endif
if (info.eq.'all'.or.info.eq.'input') then
    write(stderr,*) 'Input Options'
    write(stderr,*) '-disp DISP_FILE              Input displacements'
    write(stderr,*) '-disp:components COMPNTS     Specify displacement components'
    write(stderr,*) '-los LOS_FILE                Input line-of-sight displacements'
    write(stderr,*) '-prests PRESTS_FILE          Input pre-stresses'
    write(stderr,*) '-cov COVAR_FILE              Displacement and LOS covariances'
    write(stderr,*) '-disp:unit UNITS             Units for displacement (or to specify velocity inputs)'
    ! write(stderr,*) '-prests:weight WEIGHT        Stress weighting factor'
    ! write(stderr,*) '-prests:dist_threshold DIST  Set tractions to zero at distances>DIST'
    ! write(stderr,*) '-los:misfit MISFIT_FILE      Output RMS misfit to LOS displacements'
    write(stderr,*) '-xy                          Treat input coordinates as Cartesian'
    write(stderr,*) '-geo                         Treat input coordinates as geographic'
    write(stderr,*)
endif
if (info.eq.'all'.or.info.eq.'output') then
    write(stderr,*) 'Output Options'
    write(stderr,*) '-disp:misfit MISFIT_FILE     Output RMS misfit to displacements'
    write(stderr,*)
endif
if (info.eq.'all'.or.info.eq.'gf') then
    write(stderr,*) 'Greens Functions Options'
    write(stderr,*) '-gf:model MODEL              Greens functions calculation model'
    write(stderr,*) '-gf:disp GF_DSP_FILE         Pre-computed displacement Greens functions'
    write(stderr,*) '-gf:los GF_LOS_FILE          Pre-computed LOS displacement Greens functions'
    write(stderr,*) '-gf:sts GF_STS_FILE          Pre-computed stress Greens functions'
    write(stderr,*)
endif
if (info.eq.'all'.or.info.eq.'hafspc') then
    write(stderr,*) 'Half-space Options'
    write(stderr,*) '-haf HALFSPACE_FILE          Elastic half-space parameters'
    write(stderr,*)
endif
if (info.eq.'all'.or.info.eq.'lsqr') then
    write(stderr,*) 'Regularization Options'
    write(stderr,*) '-damp DAMP                   Damping regularization'
    write(stderr,*) '-smooth SMOOTH SMOOTH_FILE   Smoothing regularization'
    write(stderr,*)
    write(stderr,*) 'Least-Squares Options'
    write(stderr,*) '-lsqr:mode MODE              Solver algorithm'
    write(stderr,*)
endif
if (info.eq.'all'.or.info.eq.'anneal') then
    write(stderr,*) 'Simulated Annealing Options'
    write(stderr,*) '-anneal:init_mode OPT [FILE] Mode to initialize solution'
    write(stderr,*) '-anneal:step STEP_FILE       Parameter step size'
    write(stderr,*) '-anneal:max_it N             Number of iterations in search'
    write(stderr,*) '-anneal:reset_it N           Reset temperature every N iterations'
    write(stderr,*) '-anneal:temp_start T0        Starting temperature'
    write(stderr,*) '-anneal:temp_min TMIN        Minimum temperature'
    write(stderr,*) '-anneal:cool COOL_FACT       Cooling factor'
    write(stderr,*) '-anneal:log_file LOG_FILE    File with annealing progress'
    write(stderr,*) '-anneal:seed SEED            Random number generator seed'
    write(stderr,*) '-anneal:model_uncertainty    Search for model uncertainty'
    write(stderr,*) '-anneal-psc:min_flip N       Minimum number of faults to flip'
    write(stderr,*) '-anneal-psc:max_flip N       Maximum number of faults to flip'
    write(stderr,*)
endif
write(stderr,*) 'Miscellaneous Options'
write(stderr,*) '-v LEVEL                     Program verbosity'
write(stderr,*)
if (info.ne.'all') then
    write(stderr,*) 'Type "fltinv" without any arguments to see all options'
    write(stderr,*)
endif
write(stderr,*) 'See man page for details'
write(stderr,*)
call error_exit(1)
end subroutine usage
