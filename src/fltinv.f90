module fltinv

! Derived data type containing file name, array size, and array information
type fltinv_data
    character(len=512) :: file
    character(len=16) :: flag
    integer :: nrows
    integer :: ncols
    double precision, allocatable :: array(:,:)
end type fltinv_data


! General program control variables
character(len=512) :: output_file                ! Name of output fault slip file
character(len=16) :: inversion_mode              ! Inversion mode: lsqr, anneal, anneal-psc

! Fault parameters
type(fltinv_data) :: fault                       ! Sub-fault locations and geometries
type(fltinv_data) :: slip_constraint             ! Fault slip magnitude constraints
type(fltinv_data) :: rake_constraint             ! Fault rake angle constraints

! Observation constraints on fault slip
type(fltinv_data) :: displacement                ! Three-component displacement observations
character(len=3) :: disp_components              ! Components of 3-D displacement to invert for
type(fltinv_data) :: los                         ! Line-of-sight displacement observations
type(fltinv_data) :: prestress                   ! Stress field at each sub-fault
character(len=512) :: cov_file                   ! File with data covariance data
double precision, allocatable :: cov_matrix(:,:) ! Data covariance matrix

! Green's function variables
character(len=16) :: gf_model                    ! Model to compute GFs: okada_rect, okada_pt, triangle, user
type(fltinv_data) :: gf_disp                     ! Three-component displacement Green's functions
type(fltinv_data) :: gf_los                      ! Line-of-sight displacement Green's functions
type(fltinv_data) :: gf_stress                   ! Shear stress Green's functions
character(len=16) :: coord_type                  ! Fault/observation coordinate type: cartesian, geographic

! Elastic properties for calculating Green's functions
character(len=512) :: halfspace_file             ! File with elastic half-space parameters
double precision :: poisson                      ! Poisson's ratio
double precision :: shearmod                     ! Shear modulus
double precision :: lame                         ! Lame's parameter

! Regularization variables
double precision :: damping_constant             ! Minimize L1 norm of solution with this weight
double precision :: smoothing_constant           ! Minimize Laplacian roughness of solution with this weight
character(len=512) :: smoothing_file             ! Smoothing network of neighboring faults
integer :: nsmooth                               ! Number of faults to smooth
integer, allocatable :: smoothing_pointers(:,:)  ! Array to store smoothing pointers
integer, allocatable :: smoothing_neighbors(:)   ! Array to store smoothing neighbors

! Least squares variables
character(len=16) :: lsqr_mode                   ! Least-squares routine: gels, nnls
double precision, allocatable :: Asave(:,:)      ! Saved model matrix

! ! Annealing control variables
character(len=16) :: anneal_init_mode            ! Starting solution: zero, mean, rand, user
character(len=512) :: anneal_init_file           ! User initialization file
character(len=512) :: anneal_step_file           ! File with fault slip and rake step sizes
double precision, allocatable :: step(:)         ! Step size array
integer :: max_iteration                         ! Number of steps in search
integer :: reset_iteration                       ! Reset temperature every N steps
double precision :: temp_start                   ! Initial temperature
double precision :: temp_minimum                 ! Temperature does not decrease below this value
double precision :: cooling_factor               ! Reduce temp by this factor every iteration
character(len=512) :: anneal_log_file            ! Annealing log file
integer :: anneal_seed                           ! Random number seed for annealing (useful for repeating runs)
integer :: min_flip                              ! Minimum number of faults to flip locked<->unlocked
integer :: max_flip                              ! Maximum number of faults to flip locked<->unlocked
! double precision :: pl2u                         ! Probability of flipping locked to unlocked
! double precision :: pu2l                         ! Probability of flipping unlocked to locked

! Output variables
double precision, allocatable :: fault_slip(:,:)

! double precision :: los_weight                   ! Weight to line-of-sight data in inversion
! double precision :: stress_weight                ! Weight to stress data in inversion

! character(len=512) :: disp_misfit_file
! character(len=512) :: los_misfit_file

!--------------------------------------------------------------------------------------------------!
contains
!--------------------------------------------------------------------------------------------------!

subroutine init_fltinv_data(val)
!----
! Reset derived data type fltinv_data variable values
!----

implicit none

! Arguments
type(fltinv_data) :: val

val%file = 'none'
val%flag = ''
val%nrows = 0
val%ncols = 0

if (allocated(val%array)) then
    deallocate(val%array)
endif

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine read_fltinv_data_file(val,ierr)
!----
! Read derived data type fltinv_data from file and set variable values
!----

use io, only: stderr, stdout, fileExists, line_count, verbosity

implicit none

! Arguments
type(fltinv_data) :: val

! Local variables
integer :: i, j, ios, ierr
character(len=32) :: fmt_str


ierr = 0

! Check that a file should be read
if (val%file.eq.'none') then
    write(stderr,*) 'read_fltinv_data_file: no file specified'
    ierr = 1
    return
endif


! Check that ncols has been specified
if (val%ncols.le.0) then
    write(stderr,*) 'read_fltinv_data_file: ncols has not been defined'
    ierr = 1
    return
endif


! Check that the file exists, then count the number of lines
if (.not.fileExists(val%file)) then
    write(stderr,*) 'read_fltinv_data_file: no file found named "',trim(val%file),'"'
    ierr = 1
    return
endif
val%nrows = line_count(val%file)


! Allocate memory for the associated array(s)
if (.not.allocated(val%array)) then
    allocate(val%array(val%nrows,val%ncols),stat=ierr)
    if (ierr.ne.0) then
        write(stderr,*) 'read_fltinv_data_file: error allocating memory to array read from file',&
                        trim(val%file)
        return
    endif
endif


! Read the file, in free format
open(unit=21,file=val%file,status='old')
do i = 1,val%nrows
    read(21,*,iostat=ios) (val%array(i,j),j=1,val%ncols)
    if (ios.ne.0) then
        write(stderr,*) 'read_fltinv_data_file: read error on file ',trim(val%file),' at line ',i
        ierr = 1
        return
    endif
enddo
close(21)


! Print finished message
if (verbosity.ge.2) then
    write(stdout,*) 'read_fltinv_data_file: finished reading file "',trim(val%file),'"'
endif
if (verbosity.ge.3) then
    write(stdout,'(" nrows: ",I5)') val%nrows
    write(fmt_str,1001) val%ncols
    1001 format('(1P',I5,'E14.6)')
    do i = 1,val%nrows
        write(stdout,fmt=fmt_str) (val%array(i,j),j=1,val%ncols)
    enddo
endif

return
end subroutine

end module





!==================================================================================================!
!==================================================================================================!
!==================================== MAIN FLTINV PROGRAM =========================================!
!==================================================================================================!
!==================================================================================================!


program main

#ifndef USE_LAPACK
    use io, only: stderr
#endif

use fltinv, only: output_file, &
                  inversion_mode

implicit none


#ifndef USE_LAPACK
    write(stderr,*) 'fltinv does not do much compiled without LAPACK libraries'
    write(stderr,*) 'In fact, all it does is print this message and exit'
    stop
#endif


call gcmdln()
if (output_file.eq.'') then
    call usage('fltinv: output file name is required')
endif
if (inversion_mode.eq.'') then
    call usage('fltinv: inversion mode is required')
endif

call read_inputs()
call calc_disp_gfs()
call calc_stress_gfs()
call run_inversion()
call write_solution()

end program


!--------------------------------------------------------------------------------------------------!


subroutine read_inputs()
!----
! Read the input files for controlling fltinv and set up control variables and arrays. Run checks
! on inputs to make sure fltinv will operate.
!----

use io, only: stderr, stdout, line_count, fileExists, verbosity
use geom, only: strdip2normal
use elast, only: stress2traction, traction_components, read_halfspace_file
use tri_disloc, only: tri_geometry, tri_geo2cart

use fltinv, only: inversion_mode, &
                  fault, &
                  slip_constraint, &
                  rake_constraint, &
                  displacement, &
                  disp_components, &
                  los, &
                  prestress, &
                  cov_file, &
                  cov_matrix, &
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
                  read_fltinv_data_file

implicit none

! Local variables
integer :: i, j, ios, ierr, ii, nn, mm, ndisp_dof, nlos_dof, ndof
double precision :: dist, dp1, dp2, cov, sts(3,3), nvec(3), vec(3), pt1(3), pt2(3), pt3(3)
character(len=512) :: line
character(len=1) :: nchar, mchar


if (verbosity.ge.1) then
    write(stdout,*) 'read_inputs: starting'
endif


!----
! Observations are required to constrain the fault slip
!----
if (displacement%file.eq.'none'.and.los%file.eq.'none'.and.prestress%file.eq.'none') then
    call usage('read_inputs: no displacement, los, or pre-stress file defined')
else

    ! The number of columns in each file must be set prior to calling read_fltinv_data_file()

    if (displacement%file.ne.'none') then
        displacement%ncols = 6 ! x y z ux uy uz
        call read_fltinv_data_file(displacement,ierr)
        if (ierr.ne.0) then
            call usage('read_inputs: problem reading displacement file')
        endif
    endif

    if (los%file.ne.'none') then
        los%ncols = 6 ! x y z ulos az inc
        call read_fltinv_data_file(los,ierr)
        if (ierr.ne.0) then
            call usage('read_inputs: problem reading los file')
        endif
    endif

    if (prestress%file.ne.'none') then
        prestress%ncols = 6 ! sxx syy szz sxy sxz syz (fault locations defined in fault file)
        call read_fltinv_data_file(prestress,ierr)
        if (ierr.ne.0) then
            call usage('read_inputs: problem reading prestress file')
        endif
    endif
endif
if (verbosity.ge.2) then
    write(stdout,*) 'read_inputs: finished reading observation files'
endif


!----
! A fault file is required for number of sub-faults and sub-fault locations/geometries (if used)
!----
if (fault%file.eq.'none') then
    call usage('read_inputs: no fault file defined')
else

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
        call usage('Or use options such as -gf:disp_file for precomputed GFs')
    endif

    ! Read the fault data
    call read_fltinv_data_file(fault,ierr)
    if (ierr.ne.0) then
        call usage('read_inputs: problem reading fault file')
    endif

    ! Allocate memory for the output fault slip array
    if (.not.allocated(fault_slip)) then
        allocate(fault_slip(fault%nrows,2),stat=ierr)
        if (ierr.ne.0) then
            call usage('read_inputs: error allocating memory to fault_slip')
        endif
    endif
endif
if (verbosity.ge.2) then
    write(stdout,*) 'read_inputs: finished reading fault and Greens function files'
endif


!----
! Check input values
!----
! Input fault depth and dimensions depend on GF model
if (gf_model.eq.'okada_rect') then
    if (minval(fault%array(:,3)).lt.0.0d0) then
        ! Depth is positive down
        write(stderr,*) 'read_inputs: found fault depth of ',minval(fault%array(:,3))
        call usage('Depth is defined positive down for gf_model '//trim(gf_model))
    elseif (maxval(fault%array(:,3)).lt.1000.0d0) then
        ! Depth units are meters
        write(stderr,*) 'read_inputs: all fault depths are less than +1000 meters (using ', &
                        'gf_model ',trim(gf_model),')'
    endif

    ! Fault dimensions are in meters
    if (maxval(fault%array(:,6)).lt.100.0d0.and.maxval(fault%array(:,7)).lt.100.0d0) then
        write(stderr,*) 'read_inputs: all fault dimensions are less than 100 meters ',&
                        ' using gf_model ',trim(gf_model)
    endif

elseif (gf_model.eq.'okada_pt') then
    if (minval(fault%array(:,3)).lt.0.0d0) then
        ! Depth is positive down
        write(stderr,*) 'read_inputs: found fault depth of ',minval(fault%array(:,3))
        call usage('Depth is defined positive down for gf_model '//trim(gf_model))
    elseif (maxval(fault%array(:,3)).lt.1000.0d0) then
        ! Depth units are meters
        write(stderr,*) 'read_inputs: all fault depths are less than +1000 meters (using ', &
                        'gf_model ',trim(gf_model),')'
    endif

    ! Fault areas are in square meters
    if (maxval(fault%array(:,6)).lt.10000.0d0) then
        write(stderr,*) 'read_inputs: all fault areas are less than 100x100 square meters ',&
                        'using gf_model ',trim(gf_model)
    endif

elseif (gf_model.eq.'triangle') then
    if (minval([fault%array(:,3),fault%array(:,6),fault%array(:,9)]).lt.0.0d0) then
        ! Depth is positive down
        write(stderr,*) 'read_inputs: found triangle vertex depth of ', &
                        minval([fault%array(:,3),fault%array(:,6),fault%array(:,9)])
        call usage('Depth is defined positive down for gf_model '//trim(gf_model))
    endif

    if (maxval([fault%array(:,3),fault%array(:,6),fault%array(:,9)]).lt.1000.0d0) then
        ! Depth units are meters
        write(stderr,*) 'read_inputs: all fault depths are less than +1000 meters (using ', &
                        'gf_model ',trim(gf_model),')'
    endif

elseif (gf_model.eq.'precomputed') then
    ! Do nothing...the user is on the hook for getting this right
else
    call usage('read_inputs: no gf_model named "'//trim(gf_model)//'"')
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
            call usage('read_inputs: found displacement longitude greater than 360')
        elseif (minval(displacement%array(:,1)).lt.-180.0d0) then
            call usage('read_inputs: found displacement longitude less than -180')
        elseif (maxval(displacement%array(:,2)).gt.90.0d0) then
            call usage('read_inputs: found displacement latitude greater than 90')
        elseif (minval(displacement%array(:,2)).lt.-90.0d0) then
            call usage('read_inputs: found displacement latitude less than -90')
        endif
    endif

    if (los%file.ne.'none') then
        if (maxval(los%array(:,1)).gt.360.0d0) then
            call usage('read_inputs: found los longitude greater than 360')
        elseif (minval(los%array(:,1)).lt.-180.0d0) then
            call usage('read_inputs: found los longitude less than -180')
        elseif (maxval(los%array(:,2)).gt.90.0d0) then
            call usage('read_inputs: found los latitude greater than 90')
        elseif (minval(los%array(:,2)).lt.-90.0d0) then
            call usage('read_inputs: found los latitude less than -90')
        endif
    endif

    if (gf_model.eq.'okada_rect'.or.gf_model.eq.'okada_pt') then
        if (maxval(fault%array(:,1)).gt.360.0d0) then
            call usage('read_inputs: found fault longitude greater than 360')
        elseif (minval(fault%array(:,1)).lt.-180.0d0) then
            call usage('read_inputs: found fault longitude less than -180')
        elseif (maxval(fault%array(:,2)).gt.90.0d0) then
            call usage('read_inputs: found fault latitude greater than 90')
        elseif (minval(fault%array(:,2)).lt.-90.0d0) then
            call usage('read_inputs: found fault latitude less than -90')
        endif
    elseif (gf_model.eq.'triangle') then
        if (maxval(fault%array(:,1)).gt.360.0d0 .or. &
                maxval(fault%array(:,4)).gt.360.0d0 .or. &
                maxval(fault%array(:,7)).gt.360.0d0) then
            call usage('read_inputs: found fault longitude greater than 360')
        elseif (minval(fault%array(:,1)).lt.-180.0d0 .or. &
                minval(fault%array(:,4)).lt.-180.0d0 .or. &
                minval(fault%array(:,7)).lt.-180.0d0) then
            call usage('read_inputs: found fault longitude less than -180')
        elseif (maxval(fault%array(:,2)).gt.90.0d0 .or. &
                maxval(fault%array(:,5)).gt.90.0d0 .or. &
                maxval(fault%array(:,8)).gt.90.0d0) then
            call usage('read_inputs: found fault latitude greater than 90')
        elseif (minval(fault%array(:,2)).lt.-90.0d0 .or. &
                minval(fault%array(:,5)).lt.-90.0d0 .or. &
                minval(fault%array(:,8)).lt.-90.0d0) then
            call usage('read_inputs: found fault latitude less than -90')
        endif
    endif

else
    call usage('read_inputs: no coordinate type named "'//trim(coord_type)//'"')
endif


! The number of stresses must be equal to the number of faults
if (prestress%file.ne.'none') then
    if (prestress%nrows.ne.fault%nrows) then
        call usage('read_inputs: the number of pre-stresses is not equal to the number of faults')
    endif
endif

if (verbosity.ge.2) then
    write(stdout,*) 'read_inputs: finished check on observation and fault inputs; looks okay'
endif


!----
! Convert pre-stresses from tensor components to shear tractions
!----
if (prestress%file.ne.'none') then

    do i = 1,fault%nrows

        ! Get fault normal vector
        if (gf_model.eq.'okada_rect'.or.gf_model.eq.'okada_pt') then
            call strdip2normal(fault%array(i,4),fault%array(i,5),nvec)
        elseif (gf_model.eq.'triangle') then
            if (coord_type.eq.'cartesian') then
                call tri_geometry(nvec,vec,vec,fault%array(i,1:3),fault%array(i,4:6),fault%array(i,7:9))
            elseif (coord_type.eq.'geographic') then
                ! Triangle points: lon lat dep(m) to x y z
                call tri_geo2cart(pt1,pt2,pt3,fault%array(i,1:3),fault%array(i,4:6),fault%array(i,7:9),'m')
                call tri_geometry(nvec,vec,vec,pt1,pt2,pt3)
            endif
        else
            write(stderr,*) 'read_inputs: no gf_model named "'//trim(gf_model)//'"'
            write(stderr,*) 'Available models:'
            write(stderr,*) '    okada_rect'
            write(stderr,*) '    okada_pt'
            call usage(     '    triangle')
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
! Three-component displacements
if (displacement%file.ne.'none') then
    ! Assign maximum possible dimensions for three-component displacement Green's functions
    gf_disp%ncols = 2*fault%nrows

    if (gf_disp%file.ne.'none') then
        ! Read pre-computed three-component displacement Green's functions
        call read_fltinv_data_file(gf_disp,ierr)
        if (ierr.ne.0) then
            call usage('read_inputs: problem reading gf_disp file')
        endif

        ! Verify that there are the correct number of rows in the array
        if (gf_disp%nrows .ne. 3*displacement%nrows) then
            call usage('read_inputs: number of lines in three-component displacement GF file '// &
                       'must be 3*ndisplacements (one line per displacement DOF)')
        endif
    else
        ! Displacement Green's functions need to be calculated; allocate memory to array
        gf_disp%nrows = 3*displacement%nrows
        if (allocated(gf_disp%array)) then
            deallocate(gf_disp%array)
        endif
        allocate(gf_disp%array(gf_disp%nrows,gf_disp%ncols),stat=ierr)
        if (ierr.ne.0) then
            call usage('read_inputs: error allocating memory to gf_disp%array')
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
            call usage('read_inputs: problem reading gf_los file')
        endif

        ! Verify that there are the correct number of rows in the array
        if (gf_los%nrows .ne. los%nrows) then
            call usage('read_inputs: number of lines in LOS displacement GF file '// &
                       'must be ndisplacements (one line per displacement DOF)')
        endif!
    else
        ! Allocate memory to calculate Green's functions
        gf_los%nrows = los%nrows
        if (allocated(gf_los%array)) then
            deallocate(gf_los%array)
        endif
        allocate(gf_los%array(gf_los%nrows,gf_los%ncols),stat=ierr)
        if (ierr.ne.0) then
            call usage('read_inputs: error allocating memory to gf_los%array')
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
            call usage('read_inputs: problem reading gf_stress file')
        endif

        ! Verify that there are the correct number of rows in the array
        if (gf_stress%nrows .ne. 2*fault%nrows) then
            call usage('read_inputs: number of lines in pre-stress GF file must be '//&
                       '2*nfaults (one line per fault slip DOF)')
        endif
    else
        ! Allocate memory to calculate Green's functions
        gf_stress%nrows = 2*fault%nrows
        if (allocated(gf_stress%array)) then
            deallocate(gf_stress%array)
        endif
        allocate(gf_stress%array(gf_stress%nrows,gf_stress%ncols),stat=ierr)
        if (ierr.ne.0) then
            call usage('read_inputs: error allocating memory to gf_stress%array')
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
        call usage(     '    slip_ss slip_ds')
    endif

    ! There can be 1 value for all faults or 1 value for each fault
    if (slip_constraint%nrows.ne.1.and.slip_constraint%nrows.ne.fault%nrows) then
        call usage('read_inputs: number of slip constraints must be 1 or number of faults')
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
            call usage('read_inputs: error allocating memory to slip_constraint%array')
        endif
        slip_constraint%nrows = fault%nrows
        slip_constraint%array(:,1) = dp1
        slip_constraint%array(:,2) = dp2
    endif
else
    if (inversion_mode.eq.'anneal') then
        call usage('read_inputs: a slip constraint file is required for annealing search')
    elseif (inversion_mode.eq.'anneal-psc') then
        call usage('read_inputs: a slip constraint file is required for annealing+'//&
                   'pseudo-coupling search')
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
        call usage('read_inputs: problem reading rake constraint file')
    endif

    ! There can be 1 value for all faults or 1 value for each fault
    if (rake_constraint%nrows.ne.1.and.rake_constraint%nrows.ne.fault%nrows) then
        call usage('read_inputs: number of rake constraints must be 1 or number of faults')
    endif

    ! If there is only 1 value, set the full constraint array to that value
    if (rake_constraint%nrows.eq.1) then
        dp1 = rake_constraint%array(1,1)
        dp2 = rake_constraint%array(1,2)
        if (allocated(rake_constraint%array)) then
            deallocate(rake_constraint%array)
        endif
        allocate(rake_constraint%array(fault%nrows,rake_constraint%ncols),stat=ierr)
        if (ierr.ne.0) then
            call usage('read_inputs: error allocating memory to rake_constraint%array')
        endif
        rake_constraint%nrows = fault%nrows
        rake_constraint%array(:,1) = dp1
        rake_constraint%array(:,2) = dp2
    endif
else
    if (inversion_mode.eq.'anneal') then
        call usage('read_inputs: a rake constraint file is required for annealing search')
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
        call usage('read_inputs: no smoothing file found named '//trim(smoothing_file))
    endif
    nsmooth = line_count(smoothing_file)
    if (nsmooth.gt.fault%nrows) then
        call usage('read_inputs: number of faults to smooth is larger than total number of faults')
    endif

    ! Read smoothing file and calculate pointers to neighbor array
    allocate(smoothing_pointers(nsmooth,3),stat=ierr)
    if (ierr.ne.0) then
        call usage('read_inputs: error allocating memory to smoothing_pointers')
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
        call usage('read_inputs: error allocating memory to smoothing_neighbors')
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
            call usage('Expected format: iflt nnbr nbr1 nbr2...nbrn')
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
        call usage('read_inputs: error allocating memory to cov_matrix')
    endif

    ! Initialize covariance matrix as the identity matrix
    cov_matrix = 0.0d0
    do i = 1,ndisp_dof+nlos_dof
        cov_matrix(i,i) = 1.0d0
    enddo

    ! Read in the covariance values
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
                call usage('Expected format: i_obs j_obs i_comp j_comp cov')
            endif

            if (nchar.eq.'L'.or.nchar.eq.'l') then
                nn = ndisp_dof
                mm = ndisp_dof
            else
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
!-------------------------------- GREEN'S FUNCTIONS SUBROUTINES -----------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

subroutine calc_disp_gfs()
!----
! Calculate the Green's functions for three-component displacements (gf_disp) or line-of-sight
! displacements (gf_los), given the specified model (gf_model).
!----

use io, only: stdout, stderr, verbosity, debug, progress_indicator
use trig, only: d2r, r2d
use algebra, only: rotate_vector_angle_axis
use earth, only: radius_earth_m
use geom, only: lola2distaz
use okada92, only: o92_pt_disp, o92_rect_disp
use tri_disloc, only: tri_center, tri_disloc_disp

use fltinv, only: inversion_mode, &
                  fault, &
                  rake_constraint, &
                  displacement, &
                  los, &
                  gf_model, &
                  gf_disp, &
                  gf_los, &
                  poisson, &
                  shearmod, &
                  lame, &
                  coord_type

implicit none

! Local variables
integer :: i, j, ierr, iTri, ndsp, nlos, nflt
double precision :: slip_mag, rak1, rak2, slip(3), mom(4), disp1(3), disp2(3)
double precision :: sta(3), sta_new(3), dist, az, dx, dy
double precision :: evlo, evla, evdp, str, dip, wid, len, area, tri(3,4), tri_new(3,4), center(3)


if (verbosity.ge.1) then
    write(stdout,*) 'calc_disp_gfs: starting'
endif

if (displacement%file.eq.'none'.and.los%file.eq.'none') then
    if (verbosity.ge.1) then
        write(stdout,*) 'calc_disp_gfs: displacements are not used in this inversion'
        write(stdout,*)
    endif
    return
endif

if (gf_model.eq.'precomputed') then
    if (verbosity.ge.1) then
        write(stdout,*) 'calc_disp_gfs: displacement Greens functions already computed'
        write(stdout,*)
    endif
    return
endif

! Set unit slip magnitude for the Green's functions calculation
slip_mag = 1.0d0

ndsp = displacement%nrows
nlos = los%nrows
nflt = fault%nrows


! Calculate displacement Green's function for each fault-station pair
do i = 1,ndsp+nlos
    if (debug) then
        write(stderr,*) 'i:',i,' n:',ndsp+nlos
    endif

    ! Station coordinates
    if (i.le.ndsp) then
        sta(1) = displacement%array(i,1)
        sta(2) = displacement%array(i,2)
        sta(3) = displacement%array(i,3)
    else
        sta(1) = los%array(i-ndsp,1)
        sta(2) = los%array(i-ndsp,2)
        sta(3) = los%array(i-ndsp,3)
    endif
    sta_new(3) = sta(3)

    do j = 1,nflt
        if (debug) then
            write(stderr,*) 'j:',j,' n:',nflt
        endif

        ! Rake angle constraints
        if (rake_constraint%file.ne.'none'.and.inversion_mode.eq.'lsqr') then
            ! In lsqr mode, rake constraints change the rake for computing Green's functions
            rak1 = rake_constraint%array(j,1)
            if (rake_constraint%ncols.eq.2) then
                rak2 = rake_constraint%array(j,2)
            else
                rak2 = 90.0d0
            endif
        else
            ! Otherwise, calculate Green's functions for strike-slip and dip-slip sources
            rak1 = 0.0d0
            rak2 = 90.0d0
        endif
        if (debug) then
            write(stderr,*) 'rak1:',rak1,' rak2:',rak2
        endif

        ! The details of calculating the Green's functions differ depending on the model, but
        ! the general steps are the same:
        !     1. Calculate the location of the station relative to the center of the fault
        !     2. Calculate the displacement for unit slip in both rake directions

        if (gf_model.eq.'okada_rect') then
            evlo = fault%array(j,1)
            evla = fault%array(j,2)
            evdp = fault%array(j,3)
            str = fault%array(j,4)
            dip = fault%array(j,5)
            wid = fault%array(j,6)
            len = fault%array(j,7)

            ! Calculate location of station relative to fault center
            if (coord_type.eq.'geographic') then
                call lola2distaz(evlo,evla,sta(1),sta(2),dist,az)
                dist = dist*radius_earth_m
            elseif (coord_type.eq.'cartesian') then
                dx = sta(1) - evlo
                dy = sta(2) - evla
                dist = dsqrt(dx*dx+dy*dy)
                az = datan2(dx,dy)*r2d
            endif

            ! Rotate coordinate axes so that x=strike, y=horizontal updip
            sta_new(1) = dist*( dcos((az-str)*d2r))
            sta_new(2) = dist*(-dsin((az-str)*d2r))

            ! Calculate displacements for both rake angles
            slip(1) = slip_mag*cos(rak1*d2r)
            slip(2) = slip_mag*sin(rak1*d2r)
            slip(3) = 0.0d0
            call o92_rect_disp(disp1,sta_new,evdp,dip,slip,wid,len,lame,shearmod)
            slip(1) = slip_mag*cos(rak2*d2r)
            slip(2) = slip_mag*sin(rak2*d2r)
            slip(3) = 0.0d0
            call o92_rect_disp(disp2,sta_new,evdp,dip,slip,wid,len,lame,shearmod)

            ! Rotate displacements back to ENZ
            call rotate_vector_angle_axis(disp1,90.0d0-str,'z',disp1,ierr)
            call rotate_vector_angle_axis(disp2,90.0d0-str,'z',disp2,ierr)

        elseif (gf_model.eq.'okada_pt') then
            evlo = fault%array(j,1)
            evla = fault%array(j,2)
            evdp = fault%array(j,3)
            str = fault%array(j,4)
            dip = fault%array(j,5)
            area = fault%array(j,6)

            ! Calculate location of station relative to fault center
            if (coord_type.eq.'geographic') then
                call lola2distaz(evlo,evla,sta(1),sta(2),dist,az)
                dist = dist*radius_earth_m
            elseif (coord_type.eq.'cartesian') then
                dx = sta(1) - evlo
                dy = sta(2) - evla
                dist = dsqrt(dx*dx+dy*dy)
                az = datan2(dx,dy)*r2d
            endif

            ! Rotate coordinate axes so that x=strike, y=horizontal updip
            sta_new(1) = dist*( dcos((az-str)*d2r))
            sta_new(2) = dist*(-dsin((az-str)*d2r))

            ! Calculate displacements for both rake angles
            mom(1) = slip_mag*cos(rak1*d2r)*area*shearmod
            mom(2) = slip_mag*sin(rak1*d2r)*area*shearmod
            mom(3:4) = 0.0d0
            call o92_pt_disp(disp1,sta_new,evdp,dip,mom,lame,shearmod)
            mom(1) = slip_mag*cos(rak2*d2r)*area*shearmod
            mom(2) = slip_mag*sin(rak2*d2r)*area*shearmod
            mom(3:4) = 0.0d0
            call o92_pt_disp(disp2,sta_new,evdp,dip,mom,lame,shearmod)

            ! Rotate displacements back to ENZ
            call rotate_vector_angle_axis(disp1,90.0d0-str,'z',disp1,ierr)
            call rotate_vector_angle_axis(disp2,90.0d0-str,'z',disp2,ierr)

        elseif (gf_model.eq.'triangle') then
            tri(1,1) = fault%array(j,1)
            tri(2,1) = fault%array(j,2)
            tri(3,1) = fault%array(j,3)
            tri(1,2) = fault%array(j,4)
            tri(2,2) = fault%array(j,5)
            tri(3,2) = fault%array(j,6)
            tri(1,3) = fault%array(j,7)
            tri(2,3) = fault%array(j,8)
            tri(3,3) = fault%array(j,9)
            tri(:,4) = 0.0d0

            ! Calculate location of station relative to fault center
            if (coord_type.eq.'geographic') then
                call tri_center(center,tri(:,1),tri(:,2),tri(:,3))
                do iTri = 1,3
                    call lola2distaz(center(1),center(2),tri(1,iTri),tri(2,iTri),dist,az)
                    dist = dist*radius_earth_m
                    tri_new(1,iTri) = dist*dsin(az*d2r)
                    tri_new(2,iTri) = dist*dcos(az*d2r)
                    tri_new(3,iTri) = tri(3,iTri)
                enddo
                call lola2distaz(center(1),center(2),sta(1),sta(2),dist,az)
                dist = dist*radius_earth_m
                sta_new(1) = dist*dsin(az*d2r)
                sta_new(2) = dist*dcos(az*d2r)
            elseif (coord_type.eq.'cartesian') then
                tri_new = tri
                sta_new = sta
            endif

            ! Calculate displacements for both rake angles
            slip(1) = -slip_mag*cos(rak1*d2r)
            slip(2) = -slip_mag*sin(rak1*d2r)
            slip(3) = 0.0d0
            call tri_disloc_disp(disp1,sta_new,tri_new,poisson,slip)
            disp1(3) = -disp1(3) ! Returns displacement with positive z down, flip it
            slip(1) = -slip_mag*cos(rak2*d2r)
            slip(2) = -slip_mag*sin(rak2*d2r)
            slip(3) = 0.0d0
            call tri_disloc_disp(disp2,sta_new,tri_new,poisson,slip)
            disp2(3) = -disp2(3) ! Returns displacement with positive z down, flip it

        else
            call usage('calc_disp_gfs: no Greens function model named "'//trim(gf_model)//'"')
        endif
        if (debug) then
            write(stderr,*) 'disp1:',disp1
            write(stderr,*) 'disp2:',disp2
        endif

        ! Load displacements into Green's function arrays
        if (displacement%file.ne.'none'.and.i.le.ndsp) then
            gf_disp%array(i       ,j     ) = disp1(1)
            gf_disp%array(i+1*ndsp,j     ) = disp1(2)
            gf_disp%array(i+2*ndsp,j     ) = disp1(3)
            gf_disp%array(i       ,j+nflt) = disp2(1)
            gf_disp%array(i+1*ndsp,j+nflt) = disp2(2)
            gf_disp%array(i+2*ndsp,j+nflt) = disp2(3)
        endif

        if (los%file.ne.'none'.and.i.gt.ndsp) then
            disp1(1) = disp1(1)*cos(los%array(i-ndsp,6)*d2r)*sin(los%array(i-ndsp,5)*d2r) + &
                       disp1(2)*cos(los%array(i-ndsp,6)*d2r)*cos(los%array(i-ndsp,5)*d2r) - &
                       disp1(3)*sin(los%array(i-ndsp,6)*d2r)
            disp2(1) = disp2(1)*cos(los%array(i-ndsp,6)*d2r)*sin(los%array(i-ndsp,5)*d2r) + &
                       disp2(2)*cos(los%array(i-ndsp,6)*d2r)*cos(los%array(i-ndsp,5)*d2r) - &
                       disp2(3)*sin(los%array(i-ndsp,6)*d2r)
            gf_los%array(i-ndsp,j     ) = disp1(1)
            gf_los%array(i-ndsp,j+nflt) = disp2(1)
        endif
    enddo

    if (verbosity.ge.1) then
        call progress_indicator(i,ndsp+nlos,'calc_disp_gfs',ierr)
        if (ierr.gt.0) then
            call usage('calc_disp_gfs: error in progress indicator subroutine')
        endif
    endif
enddo

if (verbosity.ge.1) then
    write(stdout,*) 'calc_disp_gfs: finished'
endif
if (verbosity.ge.3) then
    write(stdout,*) 'Three-component displacement GFs:'
    write(stdout,'(2A4,6A14)') 'sta','flt','dgf_ss_x','dgf_ss_y','dgf_ss_z',&
                                           'dgf_ds_x','dgf_ds_y','dgf_ds_z'
    do i = 1,ndsp
        do j = 1,nflt
            write(stdout,'(2I4,6F14.6)') i,j, &
                                         gf_disp%array(i       ,j), &
                                         gf_disp%array(i+1*ndsp,j), &
                                         gf_disp%array(i+2*ndsp,j), &
                                         gf_disp%array(i       ,j+nflt), &
                                         gf_disp%array(i+1*ndsp,j+nflt), &
                                         gf_disp%array(i+2*ndsp,j+nflt)
        enddo
    enddo
    write(stdout,*) 'Line-of-sight displacement GFs:'
    write(stdout,'(2A4,2A14)') 'sta','flt','dgf_ss_los','dgf_ds_los'
    do i = 1,nlos
        do j = 1,nflt
            write(stdout,'(2I4,2F14.6)') i,j,gf_los%array(i,j),gf_los%array(i,j+nflt)
        enddo
    enddo
endif
if (verbosity.ge.1) then
    write(stdout,*)
endif

return
end subroutine calc_disp_gfs


!--------------------------------------------------------------------------------------------------!


subroutine calc_stress_gfs()
!----
! Calculate the Green's functions for shear tractions produced by each fault on every other fault
! (gf_stress), given the specified model (gf_model).
!----

use io, only: stdout, verbosity, progress_indicator
use trig, only: d2r, r2d
use algebra, only: rotate_matrix_angle_axis
use earth, only: radius_earth_m
use elast, only: strain2stress, stress2traction, traction_components
use geom, only: lola2distaz, strdip2normal
use okada92, only: o92_pt_strain, o92_rect_strain
use tri_disloc, only: tri_center, tri_geometry, tri_geo2cart, tri_disloc_strain

use fltinv, only: inversion_mode, &
                  fault, &
                  rake_constraint, &
                  prestress, &
                  gf_model, &
                  gf_stress, &
                  poisson, &
                  shearmod, &
                  lame, &
                  coord_type

implicit none

! Local variables
integer :: i, j, ierr, iTri
double precision, parameter :: eps = 1.0d-6
double precision :: evlo, evla, evdp, str, dip, wid, len, area, tri(3,4), tri_new(3,4), center(3)
double precision :: sta(3), sta_new(3), dist, az, dx, dy, pt1(3), pt2(3), pt3(3)
double precision :: slip_mag, rak1, rak2, slip(3), mom(4)
double precision :: strain1(3,3), strain2(3,3), stress1(3,3), stress2(3,3)
double precision :: nvec(3), tvec1(3), tvec2(3), tnor, tstr1, tupd1, tstr2, tupd2


if (verbosity.ge.1) then
    write(stdout,*) 'calc_stress_gfs: starting'
endif

if (prestress%file.eq.'none'.and.inversion_mode.ne.'anneal-psc') then
    if (verbosity.ge.1) then
        write(stdout,*) 'calc_stress_gfs: shear tractions are not used in this inversion'
        write(stdout,*)
    endif
    return
endif

if (gf_model.eq.'precomputed') then
    if (verbosity.ge.1) then
        write(stdout,*) 'calc_stress_gfs: stress Greens functions already provided'
        write(stdout,*)
    endif
    return
endif


! Set unit slip magnitude for the Green's functions calculation
slip_mag = 1.0d0


! Calculate stress Green's function for each fault-fault pair
do i = 1,fault%nrows

    ! Station coordinates
    if (gf_model.eq.'okada_rect') then
        sta(1) = fault%array(i,1)
        sta(2) = fault%array(i,2)
        sta(3) = fault%array(i,3)
    elseif (gf_model.eq.'okada_pt') then
        call usage('calc_stress_gfs: cannot compute slip from shear tractions using okada_pt '//&
                   'GFs because there is a stress singularity at point source')
    elseif (gf_model.eq.'triangle') then
        call tri_center(sta,fault%array(i,1:3),fault%array(i,4:6),fault%array(i,7:9))
    endif
    sta_new(3) = sta(3)

    ! Normal to plane at target fault
    if (gf_model.eq.'okada_rect') then
        call strdip2normal(fault%array(i,4),fault%array(i,5),nvec)
    elseif (gf_model.eq.'triangle') then
        if (coord_type.eq.'geographic') then
            call tri_geo2cart(pt1,pt2,pt3,fault%array(i,1:3),fault%array(i,4:6), &
                              fault%array(i,7:9),'m')
        endif
        call tri_geometry(nvec,tvec1,tvec2,pt1,pt2,pt3)
    endif

    do j = 1,fault%nrows

        ! Rake angles for computing Green's functions
        ! Default is strike-slip and dip-slip sources
        rak1 = 0.0d0
        rak2 = 90.0d0
        if (rake_constraint%file.ne.'none') then
            if (inversion_mode.eq.'lsqr') then
                ! lsqr: constraints change the rake for computing Green's functions
                rak1 = rake_constraint%array(j,1)
                if (rake_constraint%ncols.eq.2) then
                    rak2 = rake_constraint%array(j,2)
                endif
            elseif (inversion_mode.eq.'anneal-psc') then
                ! anneal-psc: constraint changes the rake for computing Green's functions
                rak1 = rake_constraint%array(j,1)
            endif
        endif

        ! The details of calculating the Green's functions differ depending on the model, but
        ! the general steps are the same:
        !     1. Calculate the location of the station relative to the center of the fault
        !     2. Calculate the shear tractions for unit slip in both rake directions

        if (gf_model.eq.'okada_rect') then
            evlo = fault%array(j,1)
            evla = fault%array(j,2)
            evdp = fault%array(j,3)
            str = fault%array(j,4)
            dip = fault%array(j,5)
            wid = fault%array(j,6)
            len = fault%array(j,7)

            ! Calculate location of station relative to fault center
            if (coord_type.eq.'geographic') then
                call lola2distaz(evlo,evla,sta(1),sta(2),dist,az)
                dist = dist*radius_earth_m
            elseif (coord_type.eq.'cartesian') then
                dx = sta(1) - evlo
                dy = sta(2) - evla
                dist = dsqrt(dx*dx+dy*dy)
                az = datan2(dx,dy)*r2d
            endif

            ! Rotate coordinate axes so that x=strike, y=horizontal updip
            sta_new(1) = dist*( dcos((az-str)*d2r))
            sta_new(2) = dist*(-dsin((az-str)*d2r))

            ! Avoid singular solutions with measurement point lying on fault
            if (dabs(sta_new(1)).lt.eps) then
                sta_new(1) = sta_new(1) + eps
            endif
            if (dabs(sta_new(2)).lt.eps) then
                sta_new(2) = sta_new(2) + eps
            endif
            if (dabs(evdp-sta_new(3)).lt.eps) then
                sta_new(3) = evdp + eps
            endif

            ! Calculate strains for both rake angles
            slip(1) = slip_mag*cos(rak1*d2r)
            slip(2) = slip_mag*sin(rak1*d2r)
            slip(3) = 0.0d0
            call o92_rect_strain(strain1,sta_new,evdp,dip,slip,wid,len,lame,shearmod)
            slip(1) = slip_mag*cos(rak2*d2r)
            slip(2) = slip_mag*sin(rak2*d2r)
            slip(3) = 0.0d0
            call o92_rect_strain(strain2,sta_new,evdp,dip,slip,wid,len,lame,shearmod)

            ! Rotate strains back to ENZ
            call rotate_matrix_angle_axis(strain1,90.0d0-str,'z',strain1,ierr)
            call rotate_matrix_angle_axis(strain2,90.0d0-str,'z',strain2,ierr)

        elseif (gf_model.eq.'okada_pt') then
            evlo = fault%array(j,1)
            evla = fault%array(j,2)
            evdp = fault%array(j,3)
            str = fault%array(j,4)
            dip = fault%array(j,5)
            area = fault%array(j,6)

            ! Calculate location of station relative to fault center
            if (coord_type.eq.'geographic') then
                call lola2distaz(evlo,evla,sta(1),sta(2),dist,az)
                dist = dist*radius_earth_m
            elseif (coord_type.eq.'cartesian') then
                dx = sta(1) - evlo
                dy = sta(2) - evla
                dist = dsqrt(dx*dx+dy*dy)
                az = datan2(dx,dy)*r2d
            endif

            ! Rotate coordinate axes so that x=strike, y=horizontal updip
            sta_new(1) = dist*( dcos((az-str)*d2r))
            sta_new(2) = dist*(-dsin((az-str)*d2r))

            ! Avoid singular solutions with measurement point lying on fault
            if (dabs(sta_new(1)).lt.eps) then
                sta_new(1) = sta_new(1) + eps
            endif
            if (dabs(sta_new(2)).lt.eps) then
                sta_new(2) = sta_new(2) + eps
            endif
            if (dabs(evdp-sta_new(3)).lt.eps) then
                sta_new(3) = evdp + eps
            endif

            ! Calculate strains for both rake angles
            mom(1) = slip_mag*cos(rak1*d2r)*area*shearmod
            mom(2) = slip_mag*sin(rak1*d2r)*area*shearmod
            mom(3:4) = 0.0d0
            call o92_pt_strain(strain1,sta_new,evdp,dip,mom,lame,shearmod)
            mom(1) = slip_mag*cos(rak2*d2r)*area*shearmod
            mom(2) = slip_mag*sin(rak2*d2r)*area*shearmod
            mom(3:4) = 0.0d0
            call o92_pt_strain(strain2,sta_new,evdp,dip,mom,lame,shearmod)

            ! Rotate strains back to ENZ
            call rotate_matrix_angle_axis(strain1,90.0d0-str,'z',strain1,ierr)
            call rotate_matrix_angle_axis(strain2,90.0d0-str,'z',strain2,ierr)

        elseif (gf_model.eq.'triangle') then
            tri(1,1) = fault%array(j,1)
            tri(2,1) = fault%array(j,2)
            tri(3,1) = fault%array(j,3)
            tri(1,2) = fault%array(j,4)
            tri(2,2) = fault%array(j,5)
            tri(3,2) = fault%array(j,6)
            tri(1,3) = fault%array(j,7)
            tri(2,3) = fault%array(j,8)
            tri(3,3) = fault%array(j,9)
            tri(:,4) = 0.0d0

            ! Calculate location of station relative to fault center
            if (coord_type.eq.'geographic') then
                call tri_center(center,tri(:,1),tri(:,2),tri(:,3))
                do iTri = 1,3
                    call lola2distaz(center(1),center(2),tri(1,iTri),tri(2,iTri),dist,az)
                    dist = dist*radius_earth_m
                    tri_new(1,iTri) = dist*dsin(az*d2r)
                    tri_new(2,iTri) = dist*dcos(az*d2r)
                    tri_new(3,iTri) = tri(3,iTri)
                enddo
                call lola2distaz(center(1),center(2),sta(1),sta(2),dist,az)
                dist = dist*radius_earth_m
                sta_new(1) = dist*dsin(az*d2r)
                sta_new(2) = dist*dcos(az*d2r)
            elseif (coord_type.eq.'cartesian') then
                tri_new = tri
                sta_new = sta
            endif

            ! Calculate displacements for both rake angles
            slip(1) = -slip_mag*cos(rak1*d2r)
            slip(2) = -slip_mag*sin(rak1*d2r)
            slip(3) = 0.0d0
            call tri_disloc_strain(strain1,sta_new,tri_new,poisson,slip)
            strain1(1,3) = -strain1(1,3)
            strain1(3,1) = -strain1(1,3)
            strain1(2,3) = -strain1(2,3)
            strain1(3,2) = -strain1(2,3)

            slip(1) = -slip_mag*cos(rak2*d2r)
            slip(2) = -slip_mag*sin(rak2*d2r)
            slip(3) = 0.0d0
            call tri_disloc_strain(strain2,sta_new,tri_new,poisson,slip)
            strain2(1,3) = -strain2(1,3)
            strain2(3,1) = -strain2(1,3)
            strain2(2,3) = -strain2(2,3)
            strain2(3,2) = -strain2(2,3)

        else
            call usage('calc_stress_gfs: no Greens function model named "'//trim(gf_model)//'"')
        endif

        ! Calculate shear tractions on fault from strain tensor

        ! Resolve traction vector for strike-slip source
        call strain2stress(strain1,lame,shearmod,stress1)
        call stress2traction(stress1,nvec,tvec1)
        call traction_components(tvec1,nvec,tnor,tstr1,tupd1)

        ! Resolve traction vector for dip-slip source
        call strain2stress(strain2,lame,shearmod,stress2)
        call stress2traction(stress2,nvec,tvec2)
        call traction_components(tvec2,nvec,tnor,tstr2,tupd2)

        ! Load traction components into Green's function arrays
        ! We use the opposite sign because we are interested in the fault slip that reduces the
        ! traction resolved on the structure to zero, i.e.:
        !     trac_from_slip + applied_traction = 0
        gf_stress%array(i            ,j            ) = -tstr1
        gf_stress%array(i+fault%nrows,j            ) = -tupd1
        gf_stress%array(i            ,j+fault%nrows) = -tstr2
        gf_stress%array(i+fault%nrows,j+fault%nrows) = -tupd2
    enddo

    if (verbosity.ge.1) then
        call progress_indicator(i,fault%nrows,'calc_stress_gfs',ierr)
        if (ierr.gt.0) then
            call usage('calc_stress_gfs: error in progress indicator subroutine')
        endif
    endif
enddo

if (verbosity.ge.1) then
    write(stdout,*) 'calc_stress_gfs: finished'
endif
if (verbosity.ge.3) then
    write(stdout,*) 'Fault-fault shear traction GFs'
    write(stdout,'(2A4,6A14)') 'sta','flt','sgf_ss_str','sgf_ss_dip','sgf_ds_str','sgf_ds_dip'
    do i = 1,fault%nrows
        do j = 1,fault%nrows
            write(stdout,'(2I4,1P4E14.6)') i,j, &
                                           gf_stress%array(i            ,j            ), &
                                           gf_stress%array(i+fault%nrows,j            ), &
                                           gf_stress%array(i            ,j+fault%nrows), &
                                           gf_stress%array(i+fault%nrows,j+fault%nrows)
        enddo
    enddo
endif
if (verbosity.ge.1) then
    write(stdout,*)
endif

return
end subroutine calc_stress_gfs





!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!----------------------------------- INVERSION SUBROUTINES ----------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

subroutine run_inversion()

use io, only: stdout, stderr, verbosity
use solver, only: load_array, load_constraints

use fltinv, only: inversion_mode

implicit none

! Local variables


if (verbosity.ge.1) then
    write(stdout,*) 'run_inversion: starting'
endif


if (inversion_mode.eq.'lsqr') then
    call invert_lsqr()
elseif (inversion_mode.eq.'anneal') then
    call invert_anneal()
elseif (inversion_mode.eq.'anneal-psc') then
    call invert_anneal_psc()
else
    write(stderr,*) 'run_inversion: no inversion mode named '//trim(inversion_mode)
    write(stderr,*) 'Options for inversion mode:'
    write(stderr,*) '    lsqr'
    write(stderr,*) '    anneal'
    call usage(     '    anneal-psc')
endif

if (verbosity.ge.1) then
    write(stdout,*) 'run_inversion: finished'
    write(stdout,*)
endif

return
end subroutine run_inversion





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
                  displacement, &
                  disp_components, &
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
                  fault_slip

implicit none

! Local variables
integer :: nrows, ncols, ptr_disp, ptr_los, ptr_stress, ptr_damp, ptr_smooth, ptr_nbr
integer :: nflt, nslip, ndsp, ndsp_dof, nlos, nsts, nnbr, nfixed, nfree
integer :: i, j, ierr, icomp, iflt, inbr
double precision, allocatable :: A(:,:), b(:), x(:), atmp(:,:), btmp(:)
double precision :: damping_squared, smoothing_squared
logical, allocatable :: isSlipFixed(:)


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

! Load LOS displacement GFs and data (NEED TO INCORPORATE LOS WEIGHTS!!!!!)
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

! Load covariance matrix for displacements
if (displacement%file.ne.'none'.or.los%file.ne.'none') then
    allocate(atmp(ndsp_dof+nlos,ncols),stat=ierr)
    if (ierr.ne.0) then
        call usage('invert_lsqr: error allocating memory to atmp for including covariance matrix')
    endif
    allocate(btmp(ndsp_dof+nlos),stat=ierr)
    if (ierr.ne.0) then
        call usage('invert_lsqr: error allocating memory to btmp for including covariance matrix')
    endif

    ! Compute cov_matrix^-1*A and cov_matrix^1*b
    call solve_dsysv_nrhs(cov_matrix,A(1:ndsp_dof+nlos,1:ncols),atmp,ndsp_dof+nlos,ncols,ierr)
    if (ierr.ne.0) then
        call usage('invert_lsqr: error computing cov_matrix^-1*A')
    endif
    call solve_dsysv(cov_matrix,b(1:ndsp_dof+nlos),btmp,ndsp_dof+nlos,ierr)
    if (ierr.ne.0) then
        call usage('invert_lsqr: error computing cov_matrix^-1*b')
    endif
    A(1:ndsp_dof+nlos,1:ncols) = atmp
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

return
end subroutine invert_lsqr





!--------------------------------------------------------------------------------------------------!
!--------------------------------- SIMULATED ANNEALING --------------------------------------------!
!--------------------------------------------------------------------------------------------------!

subroutine invert_anneal()
!----
! Determine fault slip with a simulated annealing search, varying slip magnitude and rake angle
! randomly to find the best fitting solution.
!----

use trig, only: d2r
use annealing, only: anneal

use fltinv, only: fault, &
                  max_iteration, &
                  reset_iteration, &
                  temp_start, &
                  temp_minimum, &
                  cooling_factor, &
                  anneal_log_file, &
                  fault_slip

implicit none

! Interface to driver subroutines
interface
    subroutine anneal_init(model,n)
        integer :: n
        double precision :: model(n)
    end subroutine
    subroutine anneal_propose(model_current,model_proposed,n)
        integer :: n
        double precision :: model_current(n)
        double precision :: model_proposed(n)
    end subroutine
    function anneal_objective(model,n)
        integer :: n
        double precision :: model(n)
        double precision :: anneal_objective
    end function
end interface

! Local variables
integer :: i, ierr, nflt_dof
double precision, allocatable :: model_best(:)
logical :: saveRejected


saveRejected = .true.

nflt_dof = 2*fault%nrows
allocate(model_best(nflt_dof),stat=ierr)
if (ierr.ne.0) then
    call usage('invert_anneal: error allocating memory to model_best')
endif

! Call anneal routine with specific driver routines anneal_init, anneal_propose, and anneal_objective
call anneal(nflt_dof, &
            model_best, &
            anneal_init, &
            anneal_propose, &
            anneal_objective, &
            max_iteration, &
            reset_iteration, &
            temp_start, &
            temp_minimum, &
            cooling_factor, &
            anneal_log_file, &
            saveRejected)

! Save results for printing
do i = 1,fault%nrows
    fault_slip(i,1) = model_best(i)*cos(model_best(i+fault%nrows)*d2r)
    fault_slip(i,2) = model_best(i)*sin(model_best(i+fault%nrows)*d2r)
enddo

return
end subroutine invert_anneal

!--------------------------------------------------------------------------------------------------!

subroutine anneal_init(model,n)
!----
! Initialize the annealing variables:
!     - iseed: random number seed
!     - model array: slip magnitude and rake angle
!     - step array: standard deviations of the Gaussian step size for slip and rake
!----

use io, only: stdout, fileExists, line_count, verbosity
use random, only: iseed, timeseed, r8_uniform_01

use fltinv, only: fault, &
                  slip_constraint, &
                  rake_constraint, &
                  anneal_init_mode, &
                  anneal_init_file, &
                  anneal_step_file, &
                  step, &
                  anneal_seed

implicit none

! Arguments
integer :: n
double precision :: model(n)

! Local variables
integer :: i, ios, nflt


if (verbosity.ge.2) then
    write(stdout,*) 'anneal_init: starting'
endif

! Initialize the random number generator seed
if (anneal_seed.eq.0) then
    iseed = -timeseed()
else
    iseed = -abs(anneal_seed)
endif

! Check array dimensions; the first nflt entries in the model array are fault slip magnitude and the
! next nflt entries are rake angles
nflt = fault%nrows
if (n.ne.2*nflt) then
    call usage('anneal_init: input n not equal to 2*nflt')
endif


! Initialize the model array values
model = 0.0d0
if (anneal_init_mode.eq.'mean'.or.anneal_init_mode.eq.'half') then
    ! Set the initial solution to the magnitude and rake halfway between their limits
    do i = 1,nflt
        model(i)      = 0.5d0*(slip_constraint%array(i,1)+slip_constraint%array(i,2))
        model(i+nflt) = 0.5d0*(rake_constraint%array(i,1)+rake_constraint%array(i,2))
    enddo

elseif (anneal_init_mode.eq.'min') then
    model(1:nflt) = slip_constraint%array(:,1)
    model(nflt+1:2*nflt) = rake_constraint%array(:,1)

elseif (anneal_init_mode.eq.'max') then
    ! Set the initial solution to zero
    model(1:nflt) = slip_constraint%array(:,2)
    model(nflt+1:2*nflt) = rake_constraint%array(:,2)

elseif (anneal_init_mode.eq.'rand') then
    ! Set the initial solution to a random value within the limits of the constraint arrays
    do i = 1,nflt
        model(i)      = slip_constraint%array(i,1) + &
                        r8_uniform_01(iseed)*(slip_constraint%array(i,2)-slip_constraint%array(i,1))
        model(i+nflt) = rake_constraint%array(i,1) + &
                        r8_uniform_01(iseed)*(rake_constraint%array(i,2)-rake_constraint%array(i,1))
    enddo

elseif (anneal_init_mode.eq.'user') then
    ! Read the initial solution from a file (slip rake)
    if (anneal_init_file.eq.'none') then
        call usage('anneal_init: no initialization file defined for anneal_init_mode=user')
    elseif (.not.fileExists(anneal_init_file)) then
        call usage('anneal_init: no anneal_init_file found named "'//trim(anneal_init_file)//'"')
    endif
    if (line_count(anneal_init_file).ne.nflt) then
        call usage('anneal_init: number of lines in anneal_init_file must be equal to nflt')
    endif
    open(unit=29,file=anneal_init_file,status='old')
    do i = 1,nflt
        read(29,*,iostat=ios) model(i),model(i+nflt)
        if (ios.ne.0) then
            call usage('anneal_init: error reading anneal init file')
        endif
    enddo
    close(29)

else
    call usage('anneal_init: no initialization mode named "'//trim(anneal_init_mode)//'"')
endif


! Set the parameter step values
allocate(step(2*nflt))
if (anneal_step_file.eq.'none') then
    step(1:nflt)        = (slip_constraint%array(:,2)-slip_constraint%array(:,1))/50.0d0
    step(nflt+1:2*nflt) = (rake_constraint%array(:,2)-rake_constraint%array(:,1))/50.0d0
else
    if (.not.fileExists(anneal_step_file)) then
        call usage('anneal_init: no anneal_step_file found named "'//trim(anneal_step_file)//'"')
    endif
    if (line_count(anneal_step_file).ne.nflt) then
        call usage('anneal_init: anneal step file must have nflt lines')
    endif
    open(unit=30,file=anneal_step_file,status='old')
    do i = 1,nflt
        read(30,*,iostat=ios) step(i),step(i+nflt)
        if (ios.ne.0) then
            call usage('anneal_init: error reading anneal step file')
        endif
    enddo
    close(30)
endif

if (verbosity.ge.2) then
    write(stdout,*) 'anneal_init: finished'
endif
if (verbosity.ge.3) then
    write(stdout,*) 'Initial model:'
    do i = 1,nflt
        write(stdout,*) model(i),model(i+nflt)
    enddo
endif


return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine anneal_propose(model_in,model_out,n)
!----
! Propose a new model from a multi-dimensional normal distribution around the current model
!----

use io, only: verbosity, stdout
use random, only: iseed, r8_normal_ab

use fltinv, only: fault, &
                  slip_constraint, &
                  rake_constraint, &
                  step

implicit none

! Arguments
integer :: n
double precision :: model_in(n), model_out(n)

! Local variables
integer :: i, j, nflt


if (verbosity.ge.3) then
    write(stdout,*) 'anneal_propose: starting'
endif


nflt = fault%nrows
if (n.ne.2*nflt) then
    call usage('anneal_propose: input n not equal to 2*nflt')
endif

do i = 1,nflt
    j = i+nflt

    ! The new model is sampled from a Gaussian distribution around the current model
    model_out(i) = r8_normal_ab(model_in(i),step(i),iseed)
    model_out(j) = r8_normal_ab(model_in(j),step(j),iseed)

    ! Make sure model values do not go beyond values in constraint files

    ! Hard boundary for slip magnitude is value set in slip_constraint%array
    if (model_out(i).lt.slip_constraint%array(i,1)) then
        model_out(i) = slip_constraint%array(i,1)
    elseif (model_out(i).gt.slip_constraint%array(i,2)) then
        model_out(i) = slip_constraint%array(i,2)
    endif

    ! If range is 360, rake angle can vary continuously
    ! Otherwise, the rake constraints are hard boundaries
    if (model_out(j).lt.rake_constraint%array(i,1)) then
        if (rake_constraint%array(i,2)-rake_constraint%array(i,1).ge.360.0d0) then
            model_out(j) = model_out(j)+360.0d0
        else
            model_out(j) = rake_constraint%array(i,1)
        endif
    elseif (model_out(j).gt.rake_constraint%array(i,2)) then
        if (rake_constraint%array(i,2)-rake_constraint%array(i,1).ge.360.0d0) then
            model_out(j) = model_out(j)-360.0d0
        else
            model_out(j) = rake_constraint%array(i,2)
        endif
    endif
enddo

if (verbosity.ge.3) then
    write(stdout,*) 'anneal_propose: finished'
    write(stdout,*) 'Proposed model:'
    do i = 1,nflt
        write(stdout,*) model_out(i),model_out(i+nflt)
    enddo
endif

return
end subroutine

!--------------------------------------------------------------------------------------------------!

function anneal_objective(model,n)
!----
! The objective function is -0.5 times the chi-squared misfit
!----

use trig, only: d2r

use fltinv, only: fault, &
                  displacement, &
                  disp_components, &
                  los, &
                  cov_matrix, &
                  gf_disp, &
                  gf_los

implicit none

! Arguments
integer, intent(in) :: n
double precision :: model(n), anneal_objective

! Local variables
integer :: i, ierr, idsp, iflt, icmp, nflt, ndsp, ndsp_dof, nlos, nobs
double precision, allocatable :: obs(:), pre(:)
double precision :: slip, rake


nflt = fault%nrows
if (n.ne.2*nflt) then
    call usage('anneal_objective: input n not equal to 2*nflt')
endif

ndsp = displacement%nrows
ndsp_dof = len_trim(disp_components)*displacement%nrows
nlos = los%nrows
nobs = ndsp_dof + nlos


! Observed displacements
if (.not.allocated(obs)) then
    allocate(obs(nobs),stat=ierr)
    if (ierr.ne.0) then
        call usage('anneal_objective: error allocating memory to obs')
    endif
endif
obs = 0.0d0

! Load observed three-component displacements
if (displacement%file.ne.'none') then
    do i = 1,len_trim(disp_components)
        read(disp_components(i:i),*) icmp
        obs((i-1)*ndsp+1:i*ndsp) = displacement%array(1:ndsp,3+icmp)
    enddo
endif

! Load observed line-of-sight displacements
if (los%file.ne.'none') then
    obs(ndsp_dof+1:ndsp_dof+nlos) = los%array(1:nlos,4)
endif

! Predicted displacements
if (.not.allocated(pre)) then
    allocate(pre(nobs),stat=ierr)
    if (ierr.ne.0) then
        call usage('anneal_objective: error allocating memory to pre')
    endif
endif
pre = 0.0d0

! Load predicted three-component displacements
do i = 1,len_trim(disp_components)
    read(disp_components(i:i),*) icmp
    do idsp = 1,ndsp
        do iflt = 1,nflt
            slip = model(iflt)
            rake = model(nflt+iflt)
            pre((i-1)*ndsp+idsp) = pre((i-1)*ndsp+idsp) + &
                                   gf_disp%array((icmp-1)*ndsp+idsp,     iflt)*slip*cos(rake*d2r)
            pre((i-1)*ndsp+idsp) = pre((i-1)*ndsp+idsp) + &
                                   gf_disp%array((icmp-1)*ndsp+idsp,nflt+iflt)*slip*sin(rake*d2r)
        enddo
    enddo
enddo

! Load predicted line-of-sight displacements
do i = 1,nlos
    do iflt = 1,nflt
        slip = model(iflt)
        rake = model(nflt+iflt)
        pre(ndsp_dof+i) = pre(ndsp_dof+i) + gf_los%array(i,     iflt)*slip*cos(rake*d2r)
        pre(ndsp_dof+i) = pre(ndsp_dof+i) + gf_los%array(i,nflt+iflt)*slip*sin(rake*d2r)
    enddo
enddo


! Calculate chi-squared
call misfit_chi2(obs,pre,cov_matrix,nobs,anneal_objective)
anneal_objective = -0.5d0*anneal_objective

return
end function





!--------------------------------------------------------------------------------------------------!
!--------------------------- SIMULATED ANNEALING WITH PSEUDO-COUPLING -----------------------------!
!--------------------------------------------------------------------------------------------------!

subroutine invert_anneal_psc()

use annealing, only: anneal

use fltinv, only: fault, &
                  max_iteration, &
                  reset_iteration, &
                  temp_start, &
                  temp_minimum, &
                  cooling_factor, &
                  anneal_log_file, &
                  fault_slip

implicit none

! Interface to driver subroutines
interface
    subroutine anneal_psc_init(model,n)
        integer :: n
        integer :: model(n)
    end subroutine
    subroutine anneal_psc_propose(model_current,model_proposed,n)
        integer :: n
        integer :: model_current(n)
        integer :: model_proposed(n)
    end subroutine
    function anneal_psc_objective(model,n)
        integer :: n
        integer :: model(n)
        double precision :: anneal_psc_objective
    end function
end interface

! Local variables
integer :: ierr, nflt
integer, allocatable :: locked(:)
logical :: saveRejected


saveRejected = .true.

nflt = fault%nrows
allocate(locked(nflt),stat=ierr)
if (ierr.ne.0) then
    call usage('invert_anneal_psc: error allocating memory to locked array')
endif

! Call anneal routine with specific driver routines anneal_psc_init, anneal_psc_propose, and
! anneal_psc_objective
call anneal(nflt, &
            locked, &
            anneal_psc_init, &
            anneal_psc_propose, &
            anneal_psc_objective, &
            max_iteration, &
            reset_iteration, &
            temp_start, &
            temp_minimum, &
            cooling_factor, &
            anneal_log_file, &
            saveRejected)

! Save results for printing
fault_slip = 0.0d0

! Calculate slip in each fault
call psc_slip(locked,nflt,fault_slip)


return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine anneal_psc_init(model,n)
!----
! Initialize the annealing with pseudo-coupling variables:
!     - iseed: random number seed
!     - model array: 0=unlocked, 1=locked
!     - Maximum size Greens function matrix for pseudo-coupling inversion
!----

use io, only: stderr, stdout, verbosity, fileExists, line_count
use random, only: iseed, timeseed, r8_uniform_01

use fltinv, only: fault, &
                  slip_constraint, &
                  anneal_init_mode, &
                  anneal_init_file, &
                  anneal_seed, &
                  min_flip, &
                  max_flip

implicit none

! Arguments
integer, intent(in) :: n
integer :: model(n)

! Local variables
integer :: i, ios, nflt, n_always_unlocked
double precision :: plocked


if (verbosity.ge.2) then
    write(stdout,*) 'anneal_psc_init: starting'
endif

! Initialize the random number generator seed
if (anneal_seed.eq.0) then
    iseed = -timeseed()
else
    iseed = -abs(anneal_seed)
endif

! Check array dimensions
nflt = fault%nrows
if (nflt.ne.n) then
    call usage('anneal_psc_init: input n not equal to nflt')
endif


! Check flipping parameters
if (min_flip.lt.1) then
    min_flip = 1
endif
n_always_unlocked = 0
do i = 1,nflt
    if (abs(slip_constraint%array(i,1)).gt.99998.0d0 .or. &
        abs(slip_constraint%array(i,2)).gt.99998.0d0) then
        n_always_unlocked = n_always_unlocked + 1
    endif
enddo
if (max_flip.gt.nflt-n_always_unlocked) then
    max_flip = nflt-n_always_unlocked
endif
! if (pl2u.le.1.0d-3) then
!     pl2u = 1.0d-3
! endif
! if (pu2l.le.1.0d-3) then
!     pu2l = 1.0d-3
! endif


! Initialize the model array values
model = 0
if (anneal_init_mode.eq.'locked') then
    ! Set the initial solution to all locked
    model = 1

elseif (anneal_init_mode.eq.'unlocked') then
    ! Set the initial solution to all unlocked
    model = 0

elseif (anneal_init_mode(1:4).eq.'rand') then
    ! Set the initial solution to randomly locked, with probability set by the last characters
    ! of the anneal_init_mode variable, e.g., rand0.25 sets probability locked = 0.25
    read(anneal_init_mode(5:len_trim(anneal_init_mode)),*,iostat=ios) plocked
    if (ios.ne.0) then
        write(stderr,*) 'anneal_psc_init: could not read initial plocked, setting to 0.5'
        plocked = 0.5d0
    endif
    do i = 1,nflt
        if (r8_uniform_01(iseed).lt.plocked) then
            model(i) = 1
        endif
    enddo

elseif (anneal_init_mode.eq.'user') then
    ! Read the initial solution from a file (slip rake)
    if (anneal_init_file.eq.'none') then
        call usage('anneal_psc_init: no initialization file defined for anneal_init_mode=user')
    elseif (.not.fileExists(anneal_init_file)) then
        call usage('anneal_psc_init: no anneal_init_file found named "'//trim(anneal_init_file)//'"')
    endif
    if (line_count(anneal_init_file).ne.nflt) then
        call usage('anneal_psc_init: number of lines in anneal_init_file must be equal to nflt')
    endif
    open(unit=29,file=anneal_init_file,status='old')
    do i = 1,nflt
        read(29,*,iostat=ios) model(i)
        if (ios.ne.0) then
            call usage('anneal_psc_init: error reading anneal init file')
        endif
    enddo
    close(29)

else
    write(stderr,*) 'anneal_psc_init: no initialization mode named "'//trim(anneal_init_mode)//'"'
    write(stderr,*) 'Options for annealing with pseudo-coupling initialization:'
    write(stderr,*) '    locked'
    write(stderr,*) '    unlocked'
    write(stderr,*) '    rand'
    call usage(     '    user')
endif

if (verbosity.ge.2) then
    write(stdout,*) 'anneal_psc_init: finished'
endif

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine anneal_psc_propose(model_in,model_out,n)

use io, only: stdout, stderr, verbosity
use random, only: iseed, r8_uniform_01

use fltinv, only: fault, &
                  slip_constraint, &
                  min_flip, &
                  max_flip

implicit none

! Arguments
integer, intent(in) :: n
integer :: model_in(n), model_out(n)

! Local variables
integer :: i, j, k, nflt, iflip, nflip, rand_fault_list(n)


if (verbosity.ge.3) then
    write(stdout,*) 'anneal_psc_propose: starting'
endif

nflt = fault%nrows
if (nflt.ne.n) then
    call usage('anneal_psc_propose: input n not equal to nflt')
endif


! Randomize fault list by switching each fault with a random fault
! Initialize random fault list array in order
do i = 1,nflt
    rand_fault_list(i) = i
enddo

do i = 1,nflt
    ! Sample another fault index with uniform distribution
    j = int(r8_uniform_01(iseed)*dble(nflt))+1
    if (j.lt.1) then
        j = 1
    elseif (j.gt.nflt) then
        j = nflt
    endif

    ! Switch faults i<->j
    k = rand_fault_list(j)
    rand_fault_list(j) = rand_fault_list(i)
    rand_fault_list(i) = k
enddo


! Initialize the proposed model as the current model
model_out = model_in


! Flip between min_flip and max_flip faults from locked to unlocked or vice versa

! Uniform probability of nflip (number of faults to flip) from min_flip to max_flip
nflip = int(r8_uniform_01(iseed)*dble(max_flip-min_flip+1))+min_flip
if (nflip.gt.max_flip) then
    nflip = max_flip
elseif (nflip.lt.min_flip) then
    nflip = min_flip
endif

! Flip the first nflip faults in the randomized list (that are not always unlocked)
iflip = 0
do i = 1,nflt

    ! Check if the fault is locked/unlocked and should be flipped
    if (model_out(rand_fault_list(i)).eq.1) then
        ! Flip this fault! Locked -> unlocked
        model_out(rand_fault_list(i)) = 0

    elseif (model_out(rand_fault_list(i)).eq.0) then
        ! Flip this fault! Unlocked -> locked
        model_out(rand_fault_list(i)) = 1

    else
        write(stderr,*) 'anneal_psc_propose: invalid locking state ',model_out(rand_fault_list(i))
        call usage('Valid states are 0=unlocked or 1=locked')
    endif

    ! Only count fault as flipped if it is not always unlocked
    if (abs(slip_constraint%array(rand_fault_list(i),1)).lt.99998.0d0 .or. &
                              abs(slip_constraint%array(rand_fault_list(i),2)).lt.99998.0d0) then
        iflip = iflip + 1
        if (iflip.ge.nflip) then
            exit
        endif
    endif
enddo

! write(0,*) 'locked? flipped nflip=',nflip
! do i = 1,nflt
!     write(0,*) i,model_in(i),model_out(i)
! enddo

if (verbosity.ge.3) then
    write(stdout,*) 'anneal_psc_propose: finished'
endif

return
end subroutine

!--------------------------------------------------------------------------------------------------!

function anneal_psc_objective(model,n)

use io, only: stdout, verbosity

use fltinv, only: fault, &
                  displacement, &
                  disp_components, &
                  los, &
                  cov_matrix, &
                  gf_disp, &
                  gf_los

implicit none

! Arguments
integer, intent(in) :: n
integer :: model(n)
double precision ::  anneal_psc_objective

! Local variables
integer :: i, ierr, iflt, idsp, icmp, nflt, ndsp, ndsp_dof, nlos, nobs
double precision :: slip(n,2)
double precision, allocatable :: pre(:), obs(:)


if (verbosity.ge.3) then
    write(stdout,*) 'anneal_psc_objective: starting'
endif

nflt = fault%nrows
if (nflt.ne.n) then
    call usage('anneal_psc_objective: input n not equal to nflt')
endif


! Calculate slip in each fault
call psc_slip(model,n,slip)


! Compute the fit to observations
ndsp = displacement%nrows
ndsp_dof = len_trim(disp_components)*displacement%nrows
nlos = los%nrows
nobs = ndsp_dof + nlos


! Observed displacements
if (.not.allocated(obs)) then
    allocate(obs(nobs),stat=ierr)
    if (ierr.ne.0) then
        call usage('anneal_psc_objective: error allocating memory to obs')
    endif
endif
obs = 0.0d0

! Load observed three-component displacements
if (displacement%file.ne.'none') then
    do i = 1,len_trim(disp_components)
        read(disp_components(i:i),*) icmp
        obs((i-1)*ndsp+1:i*ndsp) = displacement%array(1:ndsp,3+icmp)
    enddo
endif

! Load observed line-of-sight displacements
if (los%file.ne.'none') then
    obs(ndsp_dof+1:ndsp_dof+nlos) = los%array(1:nlos,4)
endif

! Predicted displacements
if (.not.allocated(pre)) then
    allocate(pre(nobs),stat=ierr)
    if (ierr.ne.0) then
        call usage('anneal_psc_objective: error allocating memory to pre')
    endif
endif
pre = 0.0d0

! Load predicted three-component displacements
do i = 1,len_trim(disp_components)
    read(disp_components(i:i),*) icmp
    do idsp = 1,ndsp
        do iflt = 1,nflt
            pre((i-1)*ndsp+idsp) = pre((i-1)*ndsp+idsp) + &
                                   gf_disp%array((icmp-1)*ndsp+idsp,     iflt)*slip(iflt,1)
            pre((i-1)*ndsp+idsp) = pre((i-1)*ndsp+idsp) + &
                                   gf_disp%array((icmp-1)*ndsp+idsp,nflt+iflt)*slip(iflt,2)
        enddo
    enddo
enddo

! Load predicted line-of-sight displacements
do i = 1,nlos
    do iflt = 1,nflt
        pre(ndsp_dof+i) = pre(ndsp_dof+i) + gf_los%array(i,     iflt)*slip(iflt,1)
        pre(ndsp_dof+i) = pre(ndsp_dof+i) + gf_los%array(i,nflt+iflt)*slip(iflt,2)
    enddo
enddo

! do i = 1,nobs
!     write(0,*) i,obs(i),pre(i)
! enddo

! Calculate chi-squared
call misfit_chi2(obs,pre,cov_matrix,nobs,anneal_psc_objective)
anneal_psc_objective = -0.5d0*anneal_psc_objective

if (verbosity.ge.3) then
    write(stdout,*) 'anneal_psc_objective: objective=',anneal_psc_objective
    write(stdout,*) 'anneal_psc_objective: finished'
endif

return
end function

!--------------------------------------------------------------------------------------------------!

subroutine psc_slip(locked,n,slip)
!----
! Solve for fault slip in unlocked patches surrounding locked patches
!----

use io, only: stderr
use solver, only: load_array, load_constraints, solve_dgesv

use fltinv, only: fault, &
                  slip_constraint, &
                  rake_constraint, &
                  prestress, &
                  gf_stress, &
                  Asave

implicit none

! Arguments
integer :: n, locked(n)
double precision :: slip(n,2)

! Local variables
integer :: i, j, ierr, nflt, nrows, ncols
logical :: doInversion, isSlipFixed(2*n)
double precision :: A(2*n,2*n), b(2*n), x(2*n)


! Set dimension variables for arrays
nflt = fault%nrows
if (n.ne.nflt) then
    call usage('psc_slip: size of locked fault array is not equal to nflt')
endif

if (rake_constraint%ncols.eq.1) then
    nrows = nflt
    ncols = nflt
elseif (rake_constraint%ncols.eq.2) then
    nrows = 2*nflt
    ncols = 2*nflt
else
    call usage('psc_slip: number of rake constraints must be 1 or 2')
endif


! Set up model matrix Asave with its maximum dimensions and completely full
if (.not.allocated(Asave)) then

    ! Allocate memory for Asave model matrix
    allocate(Asave(nrows,ncols),stat=ierr)
    if (ierr.ne.0) then
        call usage('psc_slip: error allocating memory to Asave')
    endif

    ! Load shear traction Green's functions into model matrix
    call load_array(Asave,nrows,ncols,gf_stress%array(1:nrows,1:ncols),nrows,ncols,1,1,&
                    'gf_stress%array',ierr)
    if (ierr.ne.0) then
        call usage('psc_slip: error loading shear stress GFs into model matrix')
    endif
endif


! There is no need to do the inversion if all the faults are unlocked, so check first
doInversion = .false.
do i = 1,nflt
    if (locked(i).eq.1) then
        doInversion = .true.
        exit
    endif
enddo

if (.not.doInversion) then
    slip = 0.0d0
    return
endif


! Run the inversion for fault slip

! Copy saved model matrix to active model matrix
A(1:nrows,1:ncols) = Asave

! Load pre-stresses into b vector if they have been provided
if (prestress%file.ne.'none') then
    b(1:nflt) = prestress%array(:,1)
    if (rake_constraint%ncols.eq.2) then
        b(nflt+1:2*nflt) = prestress%array(:,2)
    endif
else
    b(1:nrows) = 0.0d0
endif

! Load slip constraints for locked faults
isSlipFixed = .false.

do i = 1,nflt
    if (locked(i).eq.1.and.abs(slip_constraint%array(i,1)).lt.99998.0d0) then
        isSlipFixed(i) = .true.
    endif
    if (rake_constraint%ncols.eq.2) then
        if (locked(i).eq.1.and.abs(slip_constraint%array(i,2)).lt.99998.0d0) then
            isSlipFixed(i+nflt) = .true.
        endif
    endif
enddo

if (rake_constraint%ncols.eq.1) then
    call load_constraints(A,b,nrows,ncols,slip_constraint%array(:,1),isSlipFixed,ierr)
    if (ierr.ne.0) then
        call usage('psc_slip: error loading slip constraints into A and b arrays')
    endif
elseif (rake_constraint%ncols.eq.2) then
    call load_constraints(A,b,nrows,ncols, &
                          [slip_constraint%array(:,1),slip_constraint%array(:,2)],isSlipFixed, &
                          ierr)
    if (ierr.ne.0) then
        call usage('psc_slip: error loading slip constraints into A and b arrays')
    endif
endif

! Remove stress rows for fixed slip faults; these no longer contribute to the solution
j = 0
do i = 1,nrows
    if (isSlipFixed(i)) then
        ! Do not add row i to A
    else
        ! Add row i to A
        j = j + 1
        A(j,:) = A(i,:)
        b(j) = b(i)
    endif
enddo
nrows = j

! Solve for fault slip in unlocked faults
if (nrows.ne.ncols) then
    write(stderr,*) 'psc_slip: nrows not equal to ncols'
    write(stderr,*) '    nrows=',nrows
    write(stderr,*) '    ncols=',ncols
    call usage('Inversion routine gesv requires a square A matrix')
endif
call solve_dgesv(A(1:nrows,1:ncols),b(1:nrows),x,nrows,ierr)
if (ierr.ne.0) then
    call usage('psc_slip: error inverting for slip with gesv')
endif

! Load fault slip into the output array
j = 0
do i = 1,nflt
    if (isSlipFixed(i)) then
        ! Fixed slip is read from the slip constraint array
        slip(i,1) = slip_constraint%array(i,1)
        if (rake_constraint%ncols.eq.2) then
            slip(i,2) = slip_constraint%array(i,2)
        endif
    else
        ! Unlocked slip is read from the inverted slip array
        j = j + 1
        slip(i,1) = x(j)
        if (rake_constraint%ncols.eq.2) then
            slip(i,2) = x(j+nrows/2)
        endif
    endif
    ! write(0,*) slip(i,:)
enddo


return
end subroutine psc_slip



!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------- MISFIT ROUTINES ------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

subroutine misfit_rms(obs,pre,n,rms)
!----
! Calculate the root-mean-squared misfit between an observed and predicted set of points
!
!     RMS = sqrt(sum((obs-pre)^2)/n)
!-----

implicit none

! Arguments
integer :: n
double precision :: obs(n), pre(n), rms

! Local variables
integer :: i
double precision :: dif(n)


! Initialize RMS misfit
rms = 0.0d0

! Add differences squared to RMS
do i = 1,n
    dif(i) = obs(i)-pre(i)
    dif(i) = dif(i)*dif(i)
    rms = rms + dif(i)
enddo

! Take square root
rms = sqrt(rms/n)

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine misfit_chi2(obs,pre,cov,n,chi2)

use solver, only: solve_dsysv

implicit none

! Arguments
integer :: n
double precision :: obs(n), pre(n), cov(n,n), chi2

! Local variables
integer :: i, ierr
double precision :: dif(n), vec(n)


! Initialize chi-squared value
chi2 = 0.0d0

! Calculate difference between observed and predicted
dif = obs-pre

! Calculate dif_trans*cov^(-1)*dif
call solve_dsysv(cov,dif,vec,n,ierr)
do i = 1,n
    chi2 = chi2 + dif(i)*vec(i)
enddo

return
end subroutine





!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

subroutine write_solution()

use trig, only: d2r
use io, only: stdout, stderr, verbosity

use fltinv, only: output_file, &
                  inversion_mode, &
                  fault, &
                  rake_constraint, &
                  fault_slip

implicit none

! Local variables
integer :: i, ounit
double precision :: slip_mag


if (verbosity.ge.1) then
    write(stdout,*) 'write_solution: starting'
endif
!
! ! Print RMS misfit if specified
! if (disp_misfit_file.ne.'none') then
!     if (verbosity.ge.1) then
!         write(stderr,'(A)') 'write_solution says: writing displacement RMS misfit to '// &
!                             trim(disp_misfit_file)
!     endif
!
!     write(0,*) 'write_solution: opening misfit file'
!     open(unit=81,file=disp_misfit_file,status='unknown')
!
!     ! If rake is constrained, make an array of the correct size to use with misfit function
!     if (inversion_mode.eq.'lsqr'.and.rake_constraint%file.ne.'none' &
!                                                            .and.rake_constraint%nfields.eq.1) then
!         write(0,*) 'write_solution: writing misfit for fixed rake'
!         do i = 1,fault%nrecords
!             tmp_slip_array(i,1) = fault_slip(i,1) ! Green's functions already calculated for this rake
!             tmp_slip_array(i,2) = 0.0d0
!         enddo
!         write(81,*) disp_misfit_l2norm(tmp_slip_array)/dsqrt(dble(fault%nrecords))
!
!     ! Otherwise, just use the misfit function directly
!     else
!         write(0,*) 'write_solution: writing misfit for free rake'
!         write(81,*) disp_misfit_l2norm(fault_slip)/dsqrt(dble(fault%nrecords))
!     endif
!
!     write(0,*) 'write_solution: closing misfit file'
!     close(81)
! endif
!
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
!
! Print fault slip solution
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
        call usage(     '    anneal-psc')
    endif
enddo
5001 format(1P2E16.8)

if (verbosity.ge.1) then
    write(stdout,*) 'write_solution: finished'
endif

return
end subroutine write_solution

!--------------------------------------------------------------------------------------------------!

subroutine gcmdln()

use io, only: stdout, verbosity, debug

use fltinv, only: output_file, &
                  inversion_mode, &
                  fault, &
                  slip_constraint, &
                  rake_constraint, &
                  displacement, &
                  disp_components, &
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
                  min_flip, &
                  max_flip, &
                  init_fltinv_data

implicit none

! Local variables
integer :: i, narg
character(len=256) :: tag
! integer :: char_index


! Initialize variables
output_file = ''
inversion_mode = ''

call init_fltinv_data(fault)
call init_fltinv_data(slip_constraint)
call init_fltinv_data(rake_constraint)

call init_fltinv_data(displacement)
disp_components = '123'
call init_fltinv_data(los)
call init_fltinv_data(prestress)
cov_file = 'none'

gf_model = 'none'
call init_fltinv_data(gf_disp)
call init_fltinv_data(gf_los)
call init_fltinv_data(gf_stress)
coord_type = 'cartesian'

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

min_flip = 1
max_flip = 10
! pl2u = 0.01d0
! pu2l = 0.01d0

verbosity = 0
debug = .false.

! disp_misfit_file = 'none'
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
    elseif (trim(tag).eq.'-prests') then
        i = i + 1
        call get_command_argument(i,prestress%file)
    elseif (trim(tag).eq.'-cov') then
        i = i + 1
        call get_command_argument(i,cov_file)
!     elseif (trim(tag).eq.'-los:weight') then
!         i = i + 1
!         call get_command_argument(i,tag)
!         read(tag,*) los_weight
!     elseif (trim(tag).eq.'-disp:misfit') then
!         i = i + 1
!         call get_command_argument(i,disp_misfit_file)
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
    elseif (trim(tag).eq.'-xy') then
        coord_type = 'cartesian'
    elseif (trim(tag).eq.'-geo') then
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
        read(tag,*) verbosity
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
    write(stdout,'(" displacement%file:      ",A)') trim(displacement%file)
    write(stdout,'(" disp_components:        ",A)') trim(disp_components)
    write(stdout,'(" los%file:               ",A)') trim(los%file)
    write(stdout,'(" prestress%file:         ",A)') trim(prestress%file)
    write(stdout,'(" cov_file:               ",A)') trim(cov_file)
    write(stdout,*)
    write(stdout,'(" gf_model:               ",A)') trim(gf_model)
    write(stdout,'(" gf_disp%file:           ",A)') trim(gf_disp%file)
    write(stdout,'(" gf_los%file:            ",A)') trim(gf_los%file)
    write(stdout,'(" gf_stress%file:         ",A)') trim(gf_stress%file)
    write(stdout,'(" coord_type:             ",A)') trim(coord_type)
    write(stdout,*)
    write(stdout,'(" halfspace_file:         ",A)') trim(halfspace_file)
    write(stdout,*)
!     write(stdout,'(" disp_misfit_file:       ",A)') trim(disp_misfit_file)
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
if (string.ne.'') then
    write(stderr,*) trim(string)
    write(stderr,*)
endif
write(stderr,*) 'Usage: fltinv ...options...'
write(stderr,*)
write(stderr,*) 'General Options'
write(stderr,*) '-o OUTPUT_FILE               Output fault slip file'
write(stderr,*) '-mode INVERSION_MODE         Inversion mode'
write(stderr,*)
write(stderr,*) 'Fault Options'
write(stderr,*) '-flt FAULT_FILE              Input fault locations and geometries'
write(stderr,*) '-flt:slip SLIP_FILE          Slip magnitude constraints'
write(stderr,*) '-flt:rake RAKE_FILE          Rake angle constraints'
write(stderr,*)
write(stderr,*) 'Input Options'
write(stderr,*) '-disp DISP_FILE              Input displacements'
write(stderr,*) '-disp:components COMPNTS     Specify displacement components'
write(stderr,*) '-los LOS_FILE                Input line-of-sight displacements'
write(stderr,*) '-prests PRESTS_FILE          Input pre-stresses'
write(stderr,*) '-cov COVAR_FILE              Displacement and LOS covariances'
! write(stderr,*) '-disp:misfit MISFIT_FILE     Output RMS misfit to displacements'
! write(stderr,*) '-prests:weight WEIGHT        Stress weighting factor'
! write(stderr,*) '-prests:dist_threshold DIST  Set tractions to zero at distances>DIST'
! write(stderr,*) '-los:weight LOS_WEIGHT       LOS observation weighting factor'
! write(stderr,*) '-los:misfit MISFIT_FILE      Output RMS misfit to LOS displacements'
write(stderr,*)
write(stderr,*) 'Greens Functions Options'
write(stderr,*) '-gf:model MODEL              Greens functions calculation model'
write(stderr,*) '-gf:disp GF_DSP_FILE         Pre-computed displacement Greens functions'
write(stderr,*) '-gf:los GF_LOS_FILE          Pre-computed LOS displacement Greens functions'
write(stderr,*) '-gf:sts GF_STS_FILE          Pre-computed stress Greens functions'
write(stderr,*) '-xy                          Treat input coordinates as Cartesian'
write(stderr,*) '-geo                         Treat input coordinates as geographic'
write(stderr,*)
write(stderr,*) 'Half-space Options'
write(stderr,*) '-haf HALFSPACE_FILE          Elastic half-space parameters'
write(stderr,*)
write(stderr,*) 'Regularization Options'
write(stderr,*) '-damp DAMP                   Damping regularization'
write(stderr,*) '-smooth SMOOTH SMOOTH_FILE   Smoothing regularization'
write(stderr,*)
write(stderr,*) 'Least-Squares Options'
write(stderr,*) '-lsqr:mode MODE              Solver algorithm'
write(stderr,*)
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
write(stderr,*) '-anneal-psc:min_flip N       Minimum number of faults to flip'
write(stderr,*) '-anneal-psc:max_flip N       Maximum number of faults to flip'
! write(stderr,*) '-anneal-psc:pl2u PROB        Probability of flipping locked to unlocked'
! write(stderr,*) '-anneal-psc:pu2l PROB        Probability of flipping unlocked to locked'
! write(stderr,*) '-anneal:mcmc NSTEPS LOGFILE  Run MCMC search after annealing search'
write(stderr,*)
write(stderr,*) 'Miscellaneous Options'
write(stderr,*) '-v LEVEL                     Program verbosity'
write(stderr,*)
write(stderr,*) 'See man page for details'
write(stderr,*)
stop
end subroutine usage
