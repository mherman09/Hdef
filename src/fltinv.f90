program main
implicit none

#ifndef USELAPACK
    write(0,*) 'fltinv: I do not do much compiled without the LAPACK libraries; exiting'
    stop
#endif

call initialize_fltinv_variables()
call gcmdln()

call fltinv_readin()

call calc_greens_functions()

call run_inversion()
call write_solution()

call free_memory()

end program main

!==================================================================================================!

subroutine fltinv_readin()
!----
! Read the input files for fltinv
!----

use io, only: stderr, stdout, verbosity
use variable_module, only: inversion_mode, &
                           displacement, prestress, los, fault, &
                           gf_type, gf_disp, gf_stress, gf_los, &
                           smoothing, rake_constraint, slip_constraint, &
                           halfspace, coord_type, disp_components, disp_cov_file, disp_cov_mat, &
                           read_program_data_file
use elast, only: calc_plane_unit_vectors, trac_vector, calc_traction_components
use tri_disloc, only: tri_geometry, tri_geo2cart

implicit none

! Local variables
integer :: i, j, n, m, ios
double precision :: stress(3,3), nor(3), str(3), upd(3), traction(3), traction_comp(3), dist, az, &
                    pt1(3), pt2(3), pt3(3), cov
character(len=256) :: line
character(len=1) :: nchar, mchar


if (verbosity.eq.1.or.verbosity.eq.2) then
    write(stdout,*) 'fltinv_readin: starting'
endif


!----
! Either displacement/los or pre-stress data is required to constrain the fault slip
!----
if (displacement%file.eq.'none'.and.los%file.eq.'none'.and.prestress%file.eq.'none') then

    ! No constraints on fault slip provided: exit with error message
    call usage('fltinv_readin: no displacement, los, or pre-stress file defined')

else

    ! Set number of fields for reading input data in different formats
    ! and read data from input file.

    if (displacement%file.ne.'none') then
        displacement%nfields = 6 ! x y z ux uy uz
        call read_program_data_file(displacement)
    endif

    if (prestress%file.ne.'none') then
        prestress%nfields = 6 ! sxx syy szz sxy sxz syz (fault locations already defined)
        call read_program_data_file(prestress)
    endif

    if (los%file.ne.'none') then
        los%nfields = 6 ! x y z ulos az inc
        call read_program_data_file(los)
    endif

    if (verbosity.eq.2) then
        write(stdout,*) 'fltinv_readin: finished reading displacement/pre-stress data'
    endif
endif


!----
! A fault file is required for information about sub-fault locations and geometries
!----
if (fault%file.eq.'none') then

    ! No constraints on fault slip provided: exit with error message
    call usage('fltinv_readin: no fault file defined')

else

    ! The number of fields we need to read depends on whether we are reading in Green's functions
    ! or whether we need to calculate them. Providing ANY Green's functions via input file
    ! overrides the calculation of ALL Green's functions. Explicitly define this behavior here.

    ! If Green's functions are pre-computed (and we are inverting based on the quantity), we only
    ! need the number of lines from the file
    if (gf_disp%file.ne.'none'.and.displacement%file.ne.'none') then
        fault%nfields = 1
        gf_type = 'none'
    elseif (gf_los%file.ne.'none'.and.los%file.ne.'none') then
        fault%nfields = 1
        gf_type = 'none'
    elseif (gf_stress%file.ne.'none'.and.prestress%file.ne.'none') then
        fault%nfields = 1
        gf_type = 'none'

    ! Otherwise, we need to compute the Green's functions, so the fault file needs to be in the
    ! format specified by the user.
    elseif (gf_type.eq.'okada_rect') then
        fault%nfields = 7
    elseif (gf_type.eq.'okada_pt') then
        fault%nfields = 6
    elseif (gf_type.eq.'triangle') then
        fault%nfields = 9

    ! Need to have Green's functions for an inversion!!!
    else
        call usage('fltinv_readin: neither GF computation mode (-gf:model) nor precomputed '// &
                         'GFs (e.g., -gf:disp_file) are defined')
    endif

    ! Read the fault data
    call read_program_data_file(fault)

    if (verbosity.eq.2) then
        write(stdout,*) 'fltinv_readin: finished reading fault data'
    endif

endif


!----
! A few simple sanity checks on inputs
!----
! Checks on fault depth and dimensions look a little different for each input fault type
if (gf_type.eq.'okada_rect') then

    ! Depth must be positive down
    if (minval(fault%array(:,3)).lt.0.0d0) then
        call usage('fltinv_readin: found fault depth less than zero')

    ! Depth units are meters
    elseif (maxval(fault%array(:,3)).lt.1000.0d0) then
        write(stderr,*) '!! Warning: all fault depths are less than +1000 meters'
    endif

    ! Fault dimensions are in meters
    if (maxval(fault%array(:,6)).lt.100.0d0.and.maxval(fault%array(:,7)).lt.100.0d0) then
        write(stderr,*) '!! Warning: all fault dimensions are less than 100 meters'
    endif


elseif (gf_type.eq.'okada_pt') then

    ! Depth must be positive down
    if (minval(fault%array(:,3)).lt.0.0d0) then
        call usage('fltinv_readin: found fault depth less than zero')

    ! Depth units are meters
    elseif (maxval(fault%array(:,3)).lt.1000.0d0) then
        write(stderr,*) '!! Warning: all fault depths are less than +1000 meters'
    endif

    ! Fault areas are in square meters
    if (maxval(fault%array(:,6)).lt.10000.0d0) then
        write(stderr,*) '!! Warning: all fault areas are less than 100x100 square meters'
    endif


elseif (gf_type.eq.'triangle') then

    ! Depth must be positive down
    if (minval(fault%array(:,3)).lt.0.0d0.or. &
            maxval(fault%array(:,6)).lt.0.0d0.or. &
            maxval(fault%array(:,9)).lt.0.0d0) then
        call usage('fltinv_readin: found triangle vertex depth less than zero')
    endif

    ! Depth units are meters
    if (maxval(fault%array(:,3)).lt.1000.0d0.and. &
            maxval(fault%array(:,6)).lt.1000.0d0.and. &
            maxval(fault%array(:,9)).lt.1000.0d0) then
        write(stderr,*) '!! Warning: all fault depths are less than +1000 meters'
    endif

endif

if (verbosity.eq.2) then
    write(stdout,*) 'fltinv_readin: fault depths checked'
endif

! Checks on other coordinates depend on coordinate type
if (coord_type.eq.'cartesian') then

    ! Calculate distances between all faults and first displacement/los input
    ! In 'cartesian' mode, this value should be typically be 1e1-1e5 meters

    if (displacement%file.ne.'none'.and.gf_type.ne.'none') then
        do i = 1,fault%nrecords
            dist = (displacement%array(1,1)-fault%array(i,1))**2 + &
                       (displacement%array(1,2)-fault%array(i,2))**2
            if (dsqrt(dist).le.10.0d0) then
                write(stderr,*) '!! Warning: small fault-disp distance found (want to use -geo?)'
                exit
            endif
        enddo
    endif

    if (los%file.ne.'none'.and.gf_type.ne.'none') then
        do i = 1,fault%nrecords
            dist = (los%array(1,1)-fault%array(i,1))**2+(los%array(1,2)-fault%array(i,2))**2
            if (dsqrt(dist).le.10.0d0) then
                write(stderr,*) '!! Warning: small fault-los distance found (want to use -geo?)'
                exit
            endif
        enddo
    endif


elseif (coord_type.eq.'geographic') then
    ! I am going to operate under the assumption that if the user specifically indicated to use
    ! geographic mode, then they used geographic coordinates. At some point I will make a check.
    ! Okay, FINE, I will not be lazy. Here is your stupid stupidity check.

    ! Displacement coordinates
    if (displacement%file.ne.'none') then
        if (maxval(displacement%array(:,1)).gt.360.0d0) then
            call usage('fltinv_readin: found displacement longitude greater than 360')
        elseif (maxval(displacement%array(:,1)).lt.-180.0d0) then
            call usage('fltinv_readin: found displacement longitude less than -180')
        elseif (maxval(displacement%array(:,2)).gt.90.0d0) then
            call usage('fltinv_readin: found displacement latitude greater than 90')
        elseif (maxval(displacement%array(:,2)).lt.-90.0d0) then
            call usage('fltinv_readin: found displacement latitude less than -90')
        endif
    endif

    ! LOS coordinates
    if (los%file.ne.'none') then
        if (maxval(los%array(:,1)).gt.360.0d0) then
            call usage('fltinv_readin: found los longitude greater than 360')
        elseif (maxval(los%array(:,1)).lt.-180.0d0) then
            call usage('fltinv_readin: found los longitude less than -180')
        elseif (maxval(los%array(:,2)).gt.90.0d0) then
            call usage('fltinv_readin: found los latitude greater than 90')
        elseif (maxval(los%array(:,2)).lt.-90.0d0) then
            call usage('fltinv_readin: found los latitude less than -90')
        endif
    endif

    ! Fault coordinates
    if (gf_type.eq.'okada_rect'.or.gf_type.eq.'okada_pt') then
        if (maxval(fault%array(:,1)).gt.360.0d0) then
            call usage('fltinv_readin: found fault longitude greater than 360')
        elseif (maxval(fault%array(:,1)).lt.-180.0d0) then
            call usage('fltinv_readin: found fault longitude less than -180')
        elseif (maxval(fault%array(:,2)).gt.90.0d0) then
            call usage('fltinv_readin: found fault latitude greater than 90')
        elseif (maxval(fault%array(:,2)).lt.-90.0d0) then
            call usage('fltinv_readin: found fault latitude less than -90')
        endif
    elseif (gf_type.eq.'triangle') then
        if (maxval(fault%array(:,1)).gt.360.0d0.or. &
                maxval(fault%array(:,4)).gt.360.0d0.or. &
                maxval(fault%array(:,7)).gt.360.0d0) then
            call usage('fltinv_readin: found fault longitude greater than 360')
        elseif (maxval(fault%array(:,1)).lt.-180.0d0.or. &
                maxval(fault%array(:,4)).lt.-180.0d0.or. &
                maxval(fault%array(:,7)).lt.-180.0d0) then
            call usage('fltinv_readin: found fault longitude less than -180')
        elseif (maxval(fault%array(:,2)).gt.90.0d0.or. &
                maxval(fault%array(:,5)).gt.90.0d0.or. &
                maxval(fault%array(:,8)).gt.90.0d0) then
            call usage('fltinv_readin: found fault latitude greater than 90')
        elseif (maxval(fault%array(:,2)).lt.-90.0d0.or. &
                maxval(fault%array(:,5)).lt.-90.0d0.or. &
                maxval(fault%array(:,8)).lt.-90.0d0) then
            call usage('fltinv_readin: found fault latitude less than -90')
        endif
    endif

endif

if (verbosity.eq.2) then
    write(stdout,*) 'fltinv_readin: coordinates checked'
endif


! MAKE SURE NSTRESS=NFAULT!!




! If pre-stresses are defined, calculate the shear stresses on the faults
if (prestress%file.ne.'none') then
    ! Calculate shear stresses from stress tensor
    do i = 1,fault%nrecords
        if (gf_type.eq.'okada_rect'.or.gf_type.eq.'okada_pt') then
            call calc_plane_unit_vectors(fault%array(i,4),fault%array(i,5),nor,str,upd)
        elseif (gf_type.eq.'triangle') then
            if (coord_type.eq.'cartesian') then
                call tri_geometry(nor,str,upd,fault%array(i,1:3),fault%array(i,4:6),&
                                  fault%array(i,7:9))
            elseif (coord_type.eq.'geographic') then
                ! Triangle points: lon lat dep(m) to x y z
                call tri_geo2cart(pt1,pt2,pt3,fault%array(i,1:3),fault%array(i,4:6), &
                                  fault%array(i,7:9),'m')
                ! write(0,*) 'pt1',pt1
                ! write(0,*) 'pt2',pt2
                ! write(0,*) 'pt3',pt3
                call tri_geometry(nor,str,upd,pt1,pt2,pt3)
                ! write(0,*) 'nor',nor
                ! write(0,*) 'str',str
                ! write(0,*) 'upd',upd
            endif
        endif
        stress(1,1) = prestress%array(i,1)
        stress(2,2) = prestress%array(i,2)
        stress(3,3) = prestress%array(i,3)
        stress(1,2) = prestress%array(i,4)
        stress(2,1) = prestress%array(i,4)
        stress(1,3) = prestress%array(i,5)
        stress(3,1) = prestress%array(i,5)
        stress(2,3) = prestress%array(i,6)
        stress(3,2) = prestress%array(i,6)
        traction = trac_vector(stress,nor)
        call calc_traction_components(traction,nor,str,upd,traction_comp)
        prestress%array(i,1) = traction_comp(2)
        prestress%array(i,2) = traction_comp(3)
        ! write(0,*) 'trac_ss:',traction_comp(2)
        ! write(0,*) 'trac_ds:',traction_comp(3)
    enddo

    if (verbosity.eq.2) then
        write(stdout,*) 'fltinv_readin: shear tractions computed'
    endif

endif
! stop

!----
! Displacement Green's functions?
!----
if (displacement%file.ne.'none') then
    ! Assign array for displacement Green's functions maximum possible size to start
    gf_disp%nfields = 2*fault%nrecords

    ! Read pre-computed displacement Green's functions
    if (gf_disp%file.ne.'none') then
        call read_program_data_file(gf_disp)
        ! Verify that there are the correct number of lines here
        if (gf_disp%nrecords .ne. 3*displacement%nrecords) then
            call usage('!! read_fltinv_inputs: number of lines in displacement GF file '// &
                             'must be 3*ndisplacements (one line per displacement DOF)')
        endif

        if (verbosity.eq.2) then
            write(stdout,*) 'fltinv_readin: finished reading displacement GFs'
        endif

    ! Allocate memory to calculate Green's functions
    else
        gf_disp%nrecords = 3*displacement%nrecords
        if (allocated(gf_disp%array)) then
            deallocate(gf_disp%array)
        endif
        allocate(gf_disp%array(gf_disp%nrecords,gf_disp%nfields))

        if (verbosity.eq.2) then
            write(stdout,*) 'fltinv_readin: finished allocating memory for displacement GFs'
        endif

    endif
endif

!----
! Stress Green's functions?
!----
if (prestress%file.ne.'none') then
    ! Assign array for stress Green's functions maximum possible size to start
    gf_stress%nfields = 2*fault%nrecords

    ! Read pre-computed displacement Green's functions
    if (gf_stress%file.ne.'none') then
        call read_program_data_file(gf_stress)
        ! Verify that there are the correct number of lines here
        if (gf_stress%nrecords .ne. 2*fault%nrecords) then
            call usage('!! read_fltinv_inputs: Number of lines in stress GF file must be '//&
                             '2*nfaults (one line per fault slip DOF)')
        endif

        if (verbosity.eq.2) then
            write(stdout,*) 'fltinv_readin: finished reading stress GFs'
        endif

    ! Allocate memory to calculate Green's functions
    else
        gf_stress%nrecords = 2*fault%nrecords
        if (allocated(gf_stress%array)) then
            deallocate(gf_stress%array)
        endif
        allocate(gf_stress%array(gf_stress%nrecords,gf_stress%nfields))

        if (verbosity.eq.2) then
            write(stdout,*) 'fltinv_readin: finished allocating memory for stress GFs'
        endif

    endif
endif
if (inversion_mode.eq.'anneal-psc') then
    gf_stress%nfields = 2*fault%nrecords
    gf_stress%nrecords = 2*fault%nrecords
    if (allocated(gf_stress%array)) then
        deallocate(gf_stress%array)
    endif
    allocate(gf_stress%array(gf_stress%nrecords,gf_stress%nfields))
endif

!----
! Line-of-sight displacement Green's functions?
!----
if (los%file.ne.'none') then
    ! Assign array for LOS displacement Green's functions maximum possible size to start
    gf_los%nfields = 2*fault%nrecords

    ! Read pre-computed LOS displacement Green's functions
    if (gf_los%file.ne.'none') then
        call read_program_data_file(gf_los)
        ! Verify that there are the correct number of lines here
        if (gf_los%nrecords .ne. los%nrecords) then
            call usage('!! read_fltinv_inputs: number of lines in LOS displacement GF '// &
                             'file must be ndisplacements (one line per displacement DOF)')
        endif

        if (verbosity.eq.2) then
            write(stdout,*) 'fltinv_readin: finished reading LOS displacement GFs'
        endif

    ! Allocate memory to calculate Green's functions
    else
        gf_los%nrecords = los%nrecords
        if (allocated(gf_los%array)) then
            deallocate(gf_los%array)
        endif
        allocate(gf_los%array(gf_los%nrecords,gf_los%nfields))

        if (verbosity.eq.2) then
            write(stdout,*) 'fltinv_readin: finished allocating memory for LOS displacement GFs'
        endif

    endif
endif

!----
! Smoothing
!----
if (smoothing%file.ne.'none') then
    ! Read the smoothing linking file into a pointer array (smoothing%intarray)
    smoothing%array_type = 'int'
    smoothing%nfields = 3
    call read_program_data_file(smoothing)
    if (smoothing%nrecords.gt.fault%nrecords) then
        call usage('!! read_fltinv_inputs: number of faults to smooth is larger than '//&
                         'number of faults')
    endif

    ! Read the smoothing linking file into a fault neighbors array (smoothing_neighbors)
    call read_smoothing_neighbors()
endif

!----
! Rake constraints
!----
if (rake_constraint%file.ne.'none') then
    if (inversion_mode.eq.'lsqr') then
        ! Check how many rake constraints there are (may want to rotate rakes for nnls inversion, then need two rakes)
        open(unit=33,file=rake_constraint%file,status="old")
        read(33,'(A)') line
        read(line,*,iostat=ios) dist,az
        if (ios.eq.0) then
            rake_constraint%nfields = 2
        else
            rake_constraint%nfields = 1
        endif
        close(33)
    elseif (inversion_mode.eq.'anneal') then
        rake_constraint%nfields = 2
    else
        call usage('!! read_fltinv_inputs: I do not know how many fields rake_constraint '//&
                         'should have for this inversion mode...')
    endif
    call read_program_data_file(rake_constraint)
    if (rake_constraint%nrecords.ne.1.and.rake_constraint%nrecords.ne.fault%nrecords) then
        write(0,'(A,I5,A)') '!! read_fltinv_inputs: found ',rake_constraint%nrecords,' rake '//&
                            'constraint records'
        write(0,'(A,I5,A)') '!! and ',fault%nrecords,' input faults'
        call usage('!! Number of rake constraints must be 1 or number of faults')
    endif
endif

!----
! Slip magnitude constraints
!----
if (slip_constraint%file.ne.'none') then
    slip_constraint%nfields = 2
    call read_program_data_file(slip_constraint)
    if (slip_constraint%nrecords.ne.1.and.slip_constraint%nrecords.ne.fault%nrecords) then
        call usage('!! read_fltinv_inputs: number of slip constraints must be 1 or '//&
                         'number of faults')
    endif
    if (slip_constraint%nrecords.eq.1) then
        dist = slip_constraint%array(1,1)
        az = slip_constraint%array(1,2)
        deallocate(slip_constraint%array)
        allocate(slip_constraint%array(fault%nrecords,slip_constraint%nfields))
        slip_constraint%nrecords = fault%nrecords
        slip_constraint%array(:,1) = dist
        slip_constraint%array(:,2) = az
    endif
endif

!----
! Read half-space data or use hard-coded default values
!----
if (gf_type.eq.'okada_rect'.or.gf_type.eq.'okada_pt'.or.gf_type.eq.'triangle') then
    if (halfspace%file.ne.'none') then
        if (halfspace%flag.eq.'velodens') then
            halfspace%nfields = 3
        elseif (halfspace%flag.eq.'lame') then
            halfspace%nfields = 2
        else
            call usage('fltinv_readin: no halfspace read option named '//trim(halfspace%flag))
        endif
        call read_program_data_file(halfspace)
    else
        halfspace%flag = 'velodens'
        halfspace%nrecords = 1
        halfspace%nfields = 3
        if (allocated(halfspace%array)) then
            deallocate(halfspace%array)
        endif
        allocate(halfspace%array(halfspace%nrecords,halfspace%nfields))
        halfspace%array(1,1) = 6800.0d0
        halfspace%array(1,2) = 3926.0d0
        halfspace%array(1,3) = 3000.0d0
    endif
endif


!----
! Read covariance file
!----
if (displacement%file.ne.'none') then
    ! Covariance matrix has nrows and ncols equal to number of displacements times number of components
    i = len_trim(disp_components)*displacement%nrecords
    allocate(disp_cov_mat(i,i))

    ! Initialize it as the identity matrix
    disp_cov_mat = 0.0d0
    do i = 1,displacement%nrecords
        do j = 1,len_trim(disp_components)
            disp_cov_mat(i+(j-1)*displacement%nrecords,i+(j-1)*displacement%nrecords) = 1.0d0
        enddo
    enddo

    ! Read in the covariance values
    if (disp_cov_file.ne.'none') then
        open(unit=72,file=disp_cov_file,status='old')
        do
            read(72,*,iostat=ios) i,j,nchar,mchar,cov
            n = index(disp_components,nchar)
            m = index(disp_components,mchar)
            disp_cov_mat(i+(n-1)*displacement%nrecords,j+(m-1)*displacement%nrecords) = cov
            disp_cov_mat(j+(m-1)*displacement%nrecords,i+(n-1)*displacement%nrecords) = cov
            if (ios.ne.0) then
                exit
            endif
        enddo
        close(72)
    endif
endif

if (verbosity.ge.1) then
    write(stderr,*) 'fltinv_readin: finished'
endif

return
end subroutine fltinv_readin

!--------------------------------------------------------------------------------------------------!

subroutine read_smoothing_neighbors()
!----
! Read in file defining fault neighbors for Laplacian smoothing
!----
use io, only: verbosity, stderr
use variable_module, only: smoothing, smoothing_neighbors
implicit none
! Local variables
integer :: i, j, nentries, ineighbor, nneighbor
character(len=1024) :: iline

if (verbosity.ge.2) then
    write(stderr,'(A)') "read_smoothing_neighbors says: starting"
endif

! File has already been read once and memory has been allocated to smoothing%intarray
! The first two fields are: ifault nneighbors; calculate pointer to location of neighbors array
do i = 1,smoothing%nrecords
    if (i.eq.1) then
        smoothing%intarray(i,3) = 1
    else
        smoothing%intarray(i,3) = smoothing%intarray(i-1,3) + smoothing%intarray(i-1,2)
    endif
enddo

! Allocate memory for neighbors array
if (.not.allocated(smoothing_neighbors)) then
    nentries = smoothing%intarray(smoothing%nrecords,3) + smoothing%intarray(smoothing%nrecords,2)
    allocate(smoothing_neighbors(nentries))
endif

! Re-open the smoothing file to read the neighbors
open(unit=25,file=smoothing%file,status='old')
do i = 1,smoothing%nrecords
    read(25,'(A)') iline
    nneighbor = smoothing%intarray(i,2)
    ineighbor = smoothing%intarray(i,3)
    read(iline,*) j,j,(smoothing_neighbors(ineighbor+j-1),j=1,nneighbor)
enddo
close(25)

if (verbosity.ge.2) then
    write(stderr,'(A)') "read_smoothing_neighbors says: finished"
endif
if (verbosity.ge.3) then
    write(stderr,'(A)') '       fault  nneighbors    neighbor'
    do i = 1,smoothing%nrecords
        write(stderr,'(10I12)') (smoothing%intarray(i,j),j=1,2), &
            (smoothing_neighbors(smoothing%intarray(i,3)+j-1),j=1,smoothing%intarray(i,2))
    enddo
endif
if (verbosity.ge.2) then
    write(stderr,*)
endif

return
end subroutine read_smoothing_neighbors

!--------------------------------------------------------------------------------------------------!

subroutine calc_greens_functions()
use io, only: stderr, stdout, verbosity
use variable_module, only: inversion_mode, displacement, los, prestress, &
                           gf_type, gf_disp, gf_stress, gf_los
! use gf_module, only: calc_gf_disp_okada_rect, calc_gf_stress_okada_rect, calc_gf_los_okada_rect, &
!                      calc_gf_disp_okada_pt,   calc_gf_stress_okada_pt,   calc_gf_los_okada_pt, &
!                      calc_gf_disp_tri,        calc_gf_stress_tri,        calc_gf_los_tri
implicit none

if (verbosity.eq.1.or.verbosity.eq.2) then
    write(stdout,*) 'calc_greens_functions: starting'
endif

! Displacement Green's functions
if (displacement%file.ne.'none') then

    ! Displacement Green's functions need to be calculated
    if (gf_disp%file.eq.'none') then
        if (gf_type.eq.'okada_rect') then
            call calc_gf_disp_okada_rect()
        elseif (gf_type.eq.'okada_pt') then
            call calc_gf_disp_okada_pt()
        elseif (gf_type.eq.'triangle') then
            call calc_gf_disp_tri()
        else
            call usage('!! Error: no option to calculate Greens functions called '// &
                             trim(gf_type))
        endif
    endif

endif

! Stress Green's functions
if (prestress%file.ne.'none'.or.inversion_mode.eq.'anneal-psc') then
    ! Stress Green's functions need o be calculated
    if (gf_stress%file.eq.'none') then
        if (gf_type.eq.'okada_rect') then
            call calc_gf_stress_okada_rect()
        elseif (gf_type.eq.'okada_pt') then
            call calc_gf_stress_okada_pt()
        elseif (gf_type.eq.'triangle') then
            call calc_gf_stress_tri()
        else
            call usage('!! Error: no option to calculate Greens functions called '// &
                             trim(gf_type))
        endif
    endif
endif

! LOS displacement Green's functions
if (los%file.ne.'none') then
    ! LOS displacement Green's functions need to be calculated
    if (gf_los%file.eq.'none') then
        if (gf_type.eq.'okada_rect') then
            call calc_gf_los_okada_rect()
        elseif (gf_type.eq.'okada_pt') then
            call calc_gf_los_okada_pt()
        elseif (gf_type.eq.'triangle') then
            call calc_gf_los_tri()
        else
            call usage('!! Error: no option to calculate Greens functions called '// &
                              trim(gf_type))
        endif
    endif
endif

if (verbosity.ge.1.or.verbosity.eq.2) then
    write(stderr,*) 'calc_greens_functions: finished'
endif

return
end subroutine calc_greens_functions


!--------------------------------------------------------------------------------------------------!
! RECTANGULAR DISLOCATIONS IN AN ELASTIC HALF-SPACE
!--------------------------------------------------------------------------------------------------!

    subroutine calc_gf_disp_okada_rect()
    use io, only: stderr, verbosity
    use variable_module, only: inversion_mode, coord_type, &
                               displacement, fault, gf_disp, halfspace, rake_constraint
    implicit none
    ! Local variables
    integer :: i, j
    double precision :: slip, stlo, stla, stdp, evlo, evla, evdp, str, dip, rak, wid, len
    double precision :: dx, dy, dist, az, x, y, ux, uy, uz, uxp, uyp, theta, uhor, vp, vs, dens
    double precision, parameter :: pi=datan(1.0d0)*4.0d0, d2r = pi/180.0d0

    if (verbosity.ge.2) then
        write(stderr,'(A)') 'calc_gf_disp_okada_rect says: starting'
    endif

    ! Unit slip for displacement GF computation
    slip = 1.0d0

    ! Half-space array contains vp, vs, density
    vp = halfspace%array(1,1)
    vs = halfspace%array(1,2)
    dens = halfspace%array(1,3)

    ! Displacement Green's function for each fault-station pair
    do i = 1,displacement%nrecords
        stlo = displacement%array(i,1)
        stla = displacement%array(i,2)
        stdp = displacement%array(i,3)
        if (verbosity.ge.1) then
            write(0,'(A1,A,I3,A)',advance='no') achar(13), &
                'calc_gf_disp_okada_rect progress: ',100*i/displacement%nrecords,'%'
            if (i.eq.displacement%nrecords) then
                write(0,*)
            endif
        endif

        do j = 1,fault%nrecords
            evlo = fault%array(j,1)
            evla = fault%array(j,2)
            evdp = fault%array(j,3)
            str  = fault%array(j,4)
            dip  = fault%array(j,5)
            wid  = fault%array(j,6)
            len  = fault%array(j,7)

            ! Origin is at fault coordinate, x-axis points horizontal up-dip
            if (coord_type.eq.'geographic') then
                call ddistaz(dist,az,evlo,evla,stlo,stla)
                dist = dist*6.371d6
            else
                dx = stlo - evlo
                dy = stla - evla
                dist = dsqrt(dx*dx+dy*dy)
                az = datan2(dx,dy)
            endif
            x = dist*( dcos(az-d2r*str))
            y = dist*(-dsin(az-d2r*str))

            ! Check for rake constraints; if none, calculate strike-slip GFs
            if (rake_constraint%file.ne.'none'.and.inversion_mode.eq.'lsqr') then
                if (rake_constraint%nrecords.eq.1) then
                    rak = rake_constraint%array(1,1)
                elseif (rake_constraint%nrecords.eq.fault%nrecords) then
                    rak = rake_constraint%array(j,1)
                else
                    call usage('!! Error: incorrect number of rake constraints')
                endif
            else
                rak = 0.0d0
            endif

            ! Strike-slip Green's function
            call o92rect(uxp,uyp,uz,x,y,stdp,evdp,dip,rak,wid,len,slip,vp,vs,dens)
            ! Rotate back to original x-y coordinates
            theta = datan2(uyp,uxp)
            uhor = dsqrt(uxp*uxp+uyp*uyp)
            theta = d2r*str - theta
            ux = uhor*dsin(theta)
            uy = uhor*dcos(theta)
            gf_disp%array(i,                        j) = ux
            gf_disp%array(i+1*displacement%nrecords,j) = uy
            gf_disp%array(i+2*displacement%nrecords,j) = uz

            ! Check for rake constraints; if none, calculate dip-slip GFs
            if (rake_constraint%file.ne.'none'.and.inversion_mode.eq.'lsqr' &
                    .and.rake_constraint%nfields.eq.2) then
                if (rake_constraint%nrecords.eq.1) then
                    rak = rake_constraint%array(1,2)
                elseif (rake_constraint%nrecords.eq.fault%nrecords) then
                    rak = rake_constraint%array(j,2)
                else
                    call usage('!! Error: incorrect number of rake constraints')
                endif
            else
                rak = 90.0d0
            endif

            ! Dip-slip (or second rake) Green's function
            call o92rect(uxp,uyp,uz,x,y,stdp,evdp,dip,rak,wid,len,slip,vp,vs,dens)
            ! Rotate back to original x-y coordinates
            theta = datan2(uyp,uxp)
            uhor = dsqrt(uxp*uxp+uyp*uyp)
            theta = d2r*str - theta
            ux = uhor*dsin(theta)
            uy = uhor*dcos(theta)
            gf_disp%array(i,                        j+fault%nrecords) = ux
            gf_disp%array(i+1*displacement%nrecords,j+fault%nrecords) = uy
            gf_disp%array(i+2*displacement%nrecords,j+fault%nrecords) = uz
        enddo
    enddo

    if (verbosity.ge.2) then
        write(stderr,'(A)') 'calc_gf_disp_okada_rect says: finished'
    endif
    if (verbosity.ge.3) then
        write(stderr,'(2A4,6A14)') 'sta','flt','dgf_ss_x','dgf_ss_y','dgf_ss_z',&
                                               'dgf_ds_x','dgf_ds_y','dgf_ds_z'
        do i = 1,displacement%nrecords
            do j = 1,fault%nrecords
                write(stderr,'(2I4,6F14.6)') i,j, &
                    gf_disp%array(i                        ,j), &
                    gf_disp%array(i+1*displacement%nrecords,j), &
                    gf_disp%array(i+2*displacement%nrecords,j), &
                    gf_disp%array(i                        ,j+fault%nrecords), &
                    gf_disp%array(i+1*displacement%nrecords,j+fault%nrecords), &
                    gf_disp%array(i+2*displacement%nrecords,j+fault%nrecords)
            enddo
        enddo
    endif
    if (verbosity.ge.2) then
        write(stderr,*)
    endif

    return
    end subroutine calc_gf_disp_okada_rect

!--------------------------------------------------------------------------------------------------!

    subroutine calc_gf_stress_okada_rect()
    use io, only: stderr, verbosity
    use variable_module, only: inversion_mode, coord_type, fault, rake_constraint, halfspace, &
                               gf_stress, sts_dist
    use elast
    implicit none
    ! Local variables
    double precision :: slip, stlo, stla, stdp, sta_str, sta_dip
    double precision :: evlo, evla, evdp, str, dip, rak, wid, len
    double precision :: dist, az, dx, dy, x, y, vp, vs, dens
    double precision :: unit_normal(3), unit_strike(3), unit_updip(3)
    double precision :: strain(3,3), stress(3,3), traction(3), traction_components(3)
    double precision, parameter :: pi=datan(1.0d0)*4.0d0, d2r=pi/180.0d0
    integer :: i, j, prog

    if (verbosity.ge.2) then
        write(stderr,'(A)') 'calc_gf_stress_okada_rect says: starting'
    endif
    prog = 0

    ! Unit slip for GF computation
    slip = 1.0d0

    ! Half-space array contains vp, vs, density
    vp = halfspace%array(1,1)
    vs = halfspace%array(1,2)
    dens = halfspace%array(1,3)

    ! Shear stress Green's functions for each fault-fault pair
    do i = 1,fault%nrecords
        if (verbosity.ge.1) then
            write(0,'(A1,A,I3,A)',advance='no') achar(13), &
                'calc_gf_stress_okada_rect progress: ',100*i/fault%nrecords,'%'
            if (i.eq.fault%nrecords) then
                write(0,*)
            endif
        endif
        stlo = fault%array(i,1)
        stla = fault%array(i,2)
        stdp = fault%array(i,3)
        sta_str = fault%array(i,4)
        sta_dip = fault%array(i,5)
        call calc_plane_unit_vectors(sta_str,sta_dip,unit_normal,unit_strike,unit_updip)

        do j = 1,fault%nrecords
            evlo = fault%array(j,1)
            evla = fault%array(j,2)
            evdp = fault%array(j,3)
            str  = fault%array(j,4)
            dip  = fault%array(j,5)
            wid  = fault%array(j,6)
            len  = fault%array(j,7)

            ! Rotate coordinates to x=along-strike, y=horizontal up-dip
            if (coord_type.eq.'geographic') then
                call ddistaz(dist,az,evlo,evla,stlo,stla)
                dist = dist*6.371d6
            else
                dx = stlo - evlo
                dy = stla - evla
                dist = dsqrt(dx*dx+dy*dy)
                az = datan2(dx,dy)
            endif
             x = dist*( dcos(az-d2r*str))
             y = dist*(-dsin(az-d2r*str))

            ! Avoid singular solutions with measurement point lying on fault
            if (dabs(x).lt.1.0d0) then
                x = x + 1.0d0
            endif
            if (dabs(y).lt.1.0d0) then
                y = y + 1.0d0
            endif
            if (dabs(evdp-stdp).lt.1.0d0) then
                stdp = evdp + 1.0d0
            endif

            ! Set resolved tractions to zero if distance is larger than threshold
            dist = dsqrt(x*x+ &
                         y*y+ &
                        (stdp-evdp)*(stdp-evdp))
            if (dist.gt.sts_dist) then
                gf_stress%array(i,j) = 0.0d0
                gf_stress%array(i,j+fault%nrecords) = 0.0d0
                gf_stress%array(i+fault%nrecords,j) = 0.0d0
                gf_stress%array(i+fault%nrecords,j+fault%nrecords) = 0.0d0
                cycle
            endif

            ! Check for rake constraints; if none, calculate strike-slip GFs
            if (rake_constraint%file.ne.'none'.and.inversion_mode.eq.'lsqr') then
                if (rake_constraint%nrecords.eq.1) then
                    rak = rake_constraint%array(1,1)
                elseif (rake_constraint%nrecords.eq.fault%nrecords) then
                    rak = rake_constraint%array(j,1)
                else
                    call usage('!! Error: incorrect number of rake constraints')
                endif
            else
                rak = 0.0d0
            endif

            ! Shear stress Green's function produced by strike-slip or fixed-rake source
            call o92rectstn(strain,x,y,stdp,evdp,dip,rak,wid,len,slip,vp,vs,dens)
            ! Rotate back to original coordinates
            call rotate_strain(strain,str)
            call calc_strain_to_stress(strain,vp,vs,dens,stress)
            traction = trac_vector(stress,unit_normal)
            call calc_traction_components(traction,unit_normal,unit_strike,unit_updip, &
                                          traction_components)
            ! Load stress Green's functions array
            gf_stress%array(i               ,j               ) = -traction_components(2)
            gf_stress%array(i+fault%nrecords,j               ) = -traction_components(3)

            ! Shear stress Green's function produced by dip-slip source
            if (rake_constraint%file.ne.'none'.and.inversion_mode.eq.'lsqr' &
                    .and.rake_constraint%nfields.eq.2) then
                if (rake_constraint%nrecords.eq.1) then
                    rak = rake_constraint%array(1,2)
                elseif (rake_constraint%nrecords.eq.fault%nrecords) then
                    rak = rake_constraint%array(j,2)
                else
                    call usage('!! Error: incorrect number of rake constraints')
                endif
            else
                rak = 90.0d0
            endif

            call o92rectstn(strain,x,y,stdp,evdp,dip,rak,wid,len,slip,vp,vs,dens)
            ! Rotate back to original coordinates
            call rotate_strain(strain,str)
            call calc_strain_to_stress(strain,vp,vs,dens,stress)
            traction = trac_vector(stress,unit_normal)
            call calc_traction_components(traction,unit_normal,unit_strike,unit_updip, &
                                          traction_components)
            ! Load stress Green's functions array
            gf_stress%array(i               ,j+fault%nrecords) = -traction_components(2)
            gf_stress%array(i+fault%nrecords,j+fault%nrecords) = -traction_components(3)
        enddo
    enddo

    if (verbosity.ge.2) then
        write(stderr,'(A)') "calc_gf_stress_okada_rect says: finished"
    endif
    if (verbosity.ge.3) then
        write(stderr,'(2A4,4A14)') 'sta','flt','sgf_ss_ss','sgf_ss_ds','sgf_ds_ss','sgf_ds_ds'
        do i = 1,fault%nrecords
            do j = 1,fault%nrecords
                write(stderr,'(2I4,1P4E14.6)') i,j, &
                    gf_stress%array(i               ,j), &
                    gf_stress%array(i+fault%nrecords,j), &
                    gf_stress%array(i               ,j+fault%nrecords), &
                    gf_stress%array(i+fault%nrecords,j+fault%nrecords)
            enddo
        enddo
    endif
    if (verbosity.ge.2) then
        write(stderr,*)
    endif

    return
    end subroutine calc_gf_stress_okada_rect

!--------------------------------------------------------------------------------------------------!

    subroutine calc_gf_los_okada_rect()
    use io, only: stderr, verbosity
    use variable_module, only: inversion_mode, coord_type, &
                               los, fault, gf_los, halfspace, rake_constraint
    implicit none
    ! Local variables
    integer :: i, j
    double precision :: slip, stlo, stla, stdp, evlo, evla, evdp, str, dip, rak, wid, len
    double precision :: dx, dy, dist, az, x, y, ux, uy, uz, uxp, uyp, theta, uhor, vp, vs, dens
    double precision :: lookaz, lookinc
    double precision, parameter :: pi=datan(1.0d0)*4.0d0, d2r = pi/180.0d0

    if (verbosity.ge.2) then
        write(stderr,'(A)') 'calc_gf_los_okada_rect says: starting'
    endif

    ! Unit slip for displacement GF computation
    slip = 1.0d0

    ! Half-space array contains vp, vs, density
    vp = halfspace%array(1,1)
    vs = halfspace%array(1,2)
    dens = halfspace%array(1,3)

    ! LOS displacement Green's function for each fault-station pair
    do i = 1,los%nrecords
        stlo = los%array(i,1)
        stla = los%array(i,2)
        stdp = los%array(i,3)
        lookaz   = los%array(i,5)*d2r
        lookinc  = los%array(i,6)*d2r
        do j = 1,fault%nrecords
            evlo = fault%array(j,1)
            evla = fault%array(j,2)
            evdp = fault%array(j,3)
            str  = fault%array(j,4)
            dip  = fault%array(j,5)
            wid  = fault%array(j,6)
            len  = fault%array(j,7)

            ! Origin is at fault coordinate, x-axis points horizontal up-dip
            if (coord_type.eq.'geographic') then
                call ddistaz(dist,az,evlo,evla,stlo,stla)
                dist = dist*6.371d6
            else
                dx = stlo - evlo
                dy = stla - evla
                dist = dsqrt(dx*dx+dy*dy)
                az = datan2(dx,dy)
            endif
             x = dist*( dcos(az-d2r*str))
             y = dist*(-dsin(az-d2r*str))

            ! Check for rake constraints; if none, calculate strike-slip GFs
            if (rake_constraint%file.ne.'none'.and.inversion_mode.eq.'lsqr') then
                if (rake_constraint%nrecords.eq.1) then
                    rak = rake_constraint%array(1,1)
                elseif (rake_constraint%nrecords.eq.fault%nrecords) then
                    rak = rake_constraint%array(j,1)
                else
                    call usage('!! Error: incorrect number of rake constraints')
                endif
            else
                rak = 0.0d0
            endif

            ! Strike-slip Green's function
            call o92rect(uxp,uyp,uz,x,y,stdp,evdp,dip,rak,wid,len,slip,vp,vs,dens)
            ! Rotate back to original x-y coordinates
            theta = datan2(uyp,uxp)
            uhor = dsqrt(uxp*uxp+uyp*uyp)
            theta = d2r*str - theta
            ux = uhor*dsin(theta)
            uy = uhor*dcos(theta)
            gf_los%array(i,j) = &
                         ux*sin(lookaz)*cos(lookinc) + uy*cos(lookaz)*cos(lookinc) - uz*sin(lookinc)

            ! Check for rake constraints; if none, calculate dip-slip GFs
            if (rake_constraint%file.ne.'none'.and.inversion_mode.eq.'lsqr' &
                    .and.rake_constraint%nfields.eq.2) then
                if (rake_constraint%nrecords.eq.1) then
                    rak = rake_constraint%array(1,2)
                elseif (rake_constraint%nrecords.eq.fault%nrecords) then
                    rak = rake_constraint%array(j,2)
                else
                    call usage('!! Error: incorrect number of rake constraints')
                endif
            else
                rak = 90.0d0
            endif

            call o92rect(uxp,uyp,uz,x,y,stdp,evdp,dip,rak,wid,len,slip,vp,vs,dens)
            ! Rotate back to original x-y coordinates
            theta = datan2(uyp,uxp)
            uhor = dsqrt(uxp*uxp+uyp*uyp)
            theta = d2r*str - theta
            ux = uhor*dsin(theta)
            uy = uhor*dcos(theta)
            gf_los%array(i,j+fault%nrecords) = &
                         ux*sin(lookaz)*cos(lookinc) + uy*cos(lookaz)*cos(lookinc) - uz*sin(lookinc)
        enddo
    enddo

    if (verbosity.ge.2) then
        write(stderr,'(A)') 'calc_gf_los_okada_rect says: finished'
    endif
    if (verbosity.ge.3) then
        write(stderr,'(2A4,6A14)') 'sta','flt','losgf_ss','losgf_ds'
        do i = 1,los%nrecords
            do j = 1,fault%nrecords
                write(stderr,'(2I4,6F14.6)') i,j, &
                    gf_los%array(i,j), gf_los%array(i,j+fault%nrecords)
            enddo
        enddo
    endif
    if (verbosity.ge.2) then
        write(stderr,*)
    endif

    return
    end subroutine calc_gf_los_okada_rect

!--------------------------------------------------------------------------------------------------!
! POINT DISLOCATIONS IN AN ELASTIC HALF-SPACE
!--------------------------------------------------------------------------------------------------!

    subroutine calc_gf_disp_okada_pt()
    use io, only: stderr, verbosity
    use variable_module, only: inversion_mode, coord_type, &
                               displacement, fault, gf_disp, halfspace, rake_constraint
    implicit none
    ! Local variables
    integer :: i, j
    double precision :: slip, stlo, stla, stdp, evlo, evla, evdp, str, dip, rak, area
    double precision :: dx, dy, dist, az, x, y, ux, uy, uz, uxp, uyp, theta, uhor, vp, vs, dens
    double precision, parameter :: pi=datan(1.0d0)*4.0d0, d2r = pi/180.0d0

    if (verbosity.ge.2) then
        write(stderr,'(A)') 'calc_gf_disp_okada_pt says: starting'
    endif

    ! Unit slip for displacement GF computation
    slip = 1.0d0

    ! Half-space array contains vp, vs, density
    vp = halfspace%array(1,1)
    vs = halfspace%array(1,2)
    dens = halfspace%array(1,3)

    ! Displacement Green's function for each fault-station pair
    do i = 1,displacement%nrecords
        stlo = displacement%array(i,1)
        stla = displacement%array(i,2)
        stdp = displacement%array(i,3)
        do j = 1,fault%nrecords
            evlo = fault%array(j,1)
            evla = fault%array(j,2)
            evdp = fault%array(j,3)
            str  = fault%array(j,4)
            dip  = fault%array(j,5)
            area  = fault%array(j,6)

            ! Origin is at fault coordinate, x-axis points horizontal up-dip
            if (coord_type.eq.'geographic') then
                call ddistaz(dist,az,evlo,evla,stlo,stla)
                dist = dist*6.371d6
            else
                dx = stlo - evlo
                dy = stla - evla
                dist = dsqrt(dx*dx+dy*dy)
                az = datan2(dx,dy)
            endif
             x = dist*( dcos(az-d2r*str))
             y = dist*(-dsin(az-d2r*str))

            ! Check for rake constraints; if none, calculate strike-slip GFs
            if (rake_constraint%file.ne.'none'.and.inversion_mode.eq.'lsqr') then
                if (rake_constraint%nrecords.eq.1) then
                    rak = rake_constraint%array(1,1)
                elseif (rake_constraint%nrecords.eq.fault%nrecords) then
                    rak = rake_constraint%array(j,1)
                else
                    call usage('!! calc_gf_disp_okada_pt: incorrect number of rake '// &
                                     'constraints')
                endif
            else
                rak = 0.0d0
            endif

            ! Strike-slip Green's function
            call o92pt(uxp,uyp,uz,x,y,stdp,evdp,dip,rak,area,slip,vp,vs,dens)
            ! Rotate back to original x-y coordinates
            theta = datan2(uyp,uxp)
            uhor = dsqrt(uxp*uxp+uyp*uyp)
            theta = d2r*str - theta
            ux = uhor*dsin(theta)
            uy = uhor*dcos(theta)
            gf_disp%array(i,                        j) = ux
            gf_disp%array(i+1*displacement%nrecords,j) = uy
            gf_disp%array(i+2*displacement%nrecords,j) = uz

            ! Check for rake constraints; if none, calculate dip-slip GFs
            if (rake_constraint%file.ne.'none'.and.inversion_mode.eq.'lsqr' &
                    .and.rake_constraint%nfields.eq.2) then
                if (rake_constraint%nrecords.eq.1) then
                    rak = rake_constraint%array(1,2)
                elseif (rake_constraint%nrecords.eq.fault%nrecords) then
                    rak = rake_constraint%array(j,2)
                else
                    call usage('!! calc_gf_disp_okada_pt: incorrect number of rake '// &
                                     'constraints')
                endif
            else
                rak = 90.0d0
            endif

            ! Dip-slip (or second rake) Green's function
            call o92pt(uxp,uyp,uz,x,y,stdp,evdp,dip,rak,area,slip,vp,vs,dens)
            ! Rotate back to original x-y coordinates
            theta = datan2(uyp,uxp)
            uhor = dsqrt(uxp*uxp+uyp*uyp)
            theta = d2r*str - theta
            ux = uhor*dsin(theta)
            uy = uhor*dcos(theta)
            gf_disp%array(i,                        j+fault%nrecords) = ux
            gf_disp%array(i+1*displacement%nrecords,j+fault%nrecords) = uy
            gf_disp%array(i+2*displacement%nrecords,j+fault%nrecords) = uz
        enddo
    enddo

    if (verbosity.ge.2) then
        write(stderr,'(A)') 'calc_gf_disp_okada_pt says: finished'
    endif
    if (verbosity.ge.3) then
        write(stderr,'(2A4,6A14)') 'sta','flt','dgf_ss_x','dgf_ss_y','dgf_ss_z',&
                                               'dgf_ds_x','dgf_ds_y','dgf_ds_z'
        do i = 1,displacement%nrecords
            do j = 1,fault%nrecords
                write(stderr,'(2I4,6F14.6)') i,j, &
                    gf_disp%array(i                        ,j), &
                    gf_disp%array(i+1*displacement%nrecords,j), &
                    gf_disp%array(i+2*displacement%nrecords,j), &
                    gf_disp%array(i                        ,j+fault%nrecords), &
                    gf_disp%array(i+1*displacement%nrecords,j+fault%nrecords), &
                    gf_disp%array(i+2*displacement%nrecords,j+fault%nrecords)
            enddo
        enddo
    endif
    if (verbosity.ge.2) then
        write(stderr,*)
    endif

    return
    end subroutine calc_gf_disp_okada_pt

!--------------------------------------------------------------------------------------------------!

    subroutine calc_gf_stress_okada_pt()
    use io, only: stderr, verbosity
    use variable_module, only: inversion_mode, coord_type, fault, rake_constraint, halfspace, &
                               gf_stress
    use elast
    implicit none
    ! Local variables
    double precision :: slip, stlo, stla, stdp, sta_str, sta_dip
    double precision :: evlo, evla, evdp, str, dip, rak, area
    double precision :: dist, az, dx, dy, x, y, vp, vs, dens
    double precision :: unit_normal(3), unit_strike(3), unit_updip(3)
    double precision :: strain(3,3), stress(3,3), traction(3), traction_components(3)
    double precision, parameter :: pi=datan(1.0d0)*4.0d0, d2r=pi/180.0d0
    integer :: i, j, prog

    if (verbosity.ge.2) then
        write(stderr,'(A)') 'calc_gf_stress_okada_pt says: starting'
    endif
    prog = 0

    ! Unit slip for GF computation
    slip = 1.0d0

    ! Half-space array contains vp, vs, density
    vp = halfspace%array(1,1)
    vs = halfspace%array(1,2)
    dens = halfspace%array(1,3)

    ! Shear stress Green's functions for each fault-fault pair
    do i = 1,fault%nrecords
        if (i.gt.prog) then
            if (verbosity.ge.2) then
                write(0,'(A,I3,A)') 'calc_gf_stress_okada_pt says:',100*prog/fault%nrecords,&
                               '% complete'
            endif
            prog = prog + fault%nrecords/100
        endif
        stlo = fault%array(i,1)
        stla = fault%array(i,2)
        stdp = fault%array(i,3)
        sta_str = fault%array(i,4)
        sta_dip = fault%array(i,5)
        call calc_plane_unit_vectors(sta_str,sta_dip,unit_normal,unit_strike,unit_updip)

        do j = 1,fault%nrecords
            evlo = fault%array(j,1)
            evla = fault%array(j,2)
            evdp = fault%array(j,3)
            str  = fault%array(j,4)
            dip  = fault%array(j,5)
            area  = fault%array(j,6)

            ! Rotate coordinates to x=along-strike, y=horizontal up-dip
            if (coord_type.eq.'geographic') then
                call ddistaz(dist,az,evlo,evla,stlo,stla)
                dist = dist*6.371d6
            else
                dx = stlo - evlo
                dy = stla - evla
                dist = dsqrt(dx*dx+dy*dy)
                az = datan2(dx,dy)
            endif
             x = dist*( dcos(az-d2r*str))
             y = dist*(-dsin(az-d2r*str))

            ! Avoid singular solutions with measurement point lying on fault
            if (dabs(x).lt.1.0d0) then
                x = x + 1.0d0
            endif
            if (dabs(y).lt.1.0d0) then
                y = y + 1.0d0
            endif
            if (dabs(evdp-stdp).lt.1.0d0) then
                stdp = evdp + 1.0d0
            endif

            ! Check for rake constraints; if none, calculate strike-slip GFs
            if (rake_constraint%file.ne.'none'.and.inversion_mode.eq.'lsqr') then
                if (rake_constraint%nrecords.eq.1) then
                    rak = rake_constraint%array(1,1)
                elseif (rake_constraint%nrecords.eq.fault%nrecords) then
                    rak = rake_constraint%array(j,1)
                else
                    call usage('!! calc_gf_stress_okada_pt: incorrect number of rake '// &
                                     'constraints')
                endif
            else
                rak = 0.0d0
            endif

            ! Shear stress Green's function produced by strike-slip or fixed-rake source
            call o92ptstn(strain,x,y,stdp,evdp,dip,rak,area,slip,vp,vs,dens)
            ! Rotate back to original coordinates
            call rotate_strain(strain,str)
            call calc_strain_to_stress(strain,vp,vs,dens,stress)
            traction = trac_vector(stress,unit_normal)
            call calc_traction_components(traction,unit_normal,unit_strike,unit_updip, &
                                          traction_components)
            ! Load stress Green's functions array
            gf_stress%array(i               ,j               ) = -traction_components(2)
            gf_stress%array(i+fault%nrecords,j               ) = -traction_components(3)

            ! Shear stress Green's function produced by dip-slip source
            if (rake_constraint%file.ne.'none'.and.inversion_mode.eq.'lsqr' &
                    .and.rake_constraint%nfields.eq.2) then
                if (rake_constraint%nrecords.eq.1) then
                    rak = rake_constraint%array(1,2)
                elseif (rake_constraint%nrecords.eq.fault%nrecords) then
                    rak = rake_constraint%array(j,2)
                else
                    call usage('!! calc_gf_stress_okada_pt: incorrect number of rake '// &
                                     'constraints')
                endif
            else
                rak = 90.0d0
            endif

            call o92ptstn(strain,x,y,stdp,evdp,dip,rak,area,slip,vp,vs,dens)
            ! Rotate back to original coordinates
            call rotate_strain(strain,str)
            call calc_strain_to_stress(strain,vp,vs,dens,stress)
            traction = trac_vector(stress,unit_normal)
            call calc_traction_components(traction,unit_normal,unit_strike,unit_updip, &
                                          traction_components)
            ! Load stress Green's functions array
            gf_stress%array(i               ,j+fault%nrecords) = -traction_components(2)
            gf_stress%array(i+fault%nrecords,j+fault%nrecords) = -traction_components(3)
        enddo
    enddo

    if (verbosity.ge.2) then
        write(stderr,'(A)') "calc_gf_stress_okada_pt says: finished"
    endif
    if (verbosity.ge.3) then
        write(stderr,'(2A4,4A14)') 'sta','flt','sgf_ss_ss','sgf_ss_ds','sgf_ds_ss','sgf_ds_ds'
        do i = 1,fault%nrecords
            do j = 1,fault%nrecords
                write(stderr,'(2I4,1P4E14.6)') i,j, &
                    gf_stress%array(i               ,j), &
                    gf_stress%array(i+fault%nrecords,j), &
                    gf_stress%array(i               ,j+fault%nrecords), &
                    gf_stress%array(i+fault%nrecords,j+fault%nrecords)
            enddo
        enddo
    endif
    if (verbosity.ge.2) then
        write(stderr,*)
    endif

    return
    end subroutine calc_gf_stress_okada_pt

!--------------------------------------------------------------------------------------------------!

    subroutine calc_gf_los_okada_pt()
    use io, only: stderr, verbosity
    use variable_module, only: inversion_mode, coord_type, &
                               los, fault, gf_los, halfspace, rake_constraint
    implicit none
    ! Local variables
    integer :: i, j
    double precision :: slip, stlo, stla, stdp, evlo, evla, evdp, str, dip, rak, area
    double precision :: dx, dy, dist, az, x, y, ux, uy, uz, uxp, uyp, theta, uhor, vp, vs, dens
    double precision :: lookaz, lookinc
    double precision, parameter :: pi=datan(1.0d0)*4.0d0, d2r = pi/180.0d0

    if (verbosity.ge.2) then
        write(stderr,'(A)') 'calc_gf_los_okada_pt says: starting'
    endif

    ! Unit slip for displacement GF computation
    slip = 1.0d0

    ! Half-space array contains vp, vs, density
    vp = halfspace%array(1,1)
    vs = halfspace%array(1,2)
    dens = halfspace%array(1,3)

    ! LOS displacement Green's function for each fault-station pair
    do i = 1,los%nrecords
        stlo = los%array(i,1)
        stla = los%array(i,2)
        stdp = los%array(i,3)
        lookaz   = los%array(i,5)*d2r
        lookinc  = los%array(i,6)*d2r
        do j = 1,fault%nrecords
            evlo = fault%array(j,1)
            evla = fault%array(j,2)
            evdp = fault%array(j,3)
            str  = fault%array(j,4)
            dip  = fault%array(j,5)
            area  = fault%array(j,6)

            ! Origin is at fault coordinate, x-axis points horizontal up-dip
            if (coord_type.eq.'geographic') then
                call ddistaz(dist,az,evlo,evla,stlo,stla)
                dist = dist*6.371d6
            else
                dx = stlo - evlo
                dy = stla - evla
                dist = dsqrt(dx*dx+dy*dy)
                az = datan2(dx,dy)
            endif
             x = dist*( dcos(az-d2r*str))
             y = dist*(-dsin(az-d2r*str))

            ! Check for rake constraints; if none, calculate strike-slip GFs
            if (rake_constraint%file.ne.'none'.and.inversion_mode.eq.'lsqr') then
                if (rake_constraint%nrecords.eq.1) then
                    rak = rake_constraint%array(1,1)
                elseif (rake_constraint%nrecords.eq.fault%nrecords) then
                    rak = rake_constraint%array(j,1)
                else
                    call usage('!! calc_gf_los_okada_pt: incorrect number of rake '// &
                                     'constraints')
                endif
            else
                rak = 0.0d0
            endif

            ! Strike-slip Green's function
            call o92pt(uxp,uyp,uz,x,y,stdp,evdp,dip,rak,area,slip,vp,vs,dens)
            ! Rotate back to original x-y coordinates
            theta = datan2(uyp,uxp)
            uhor = dsqrt(uxp*uxp+uyp*uyp)
            theta = d2r*str - theta
            ux = uhor*dsin(theta)
            uy = uhor*dcos(theta)
            gf_los%array(i,j) = &
                         ux*sin(lookaz)*cos(lookinc) + uy*cos(lookaz)*cos(lookinc) - uz*sin(lookinc)

            ! Check for rake constraints; if none, calculate dip-slip GFs
            if (rake_constraint%file.ne.'none'.and.inversion_mode.eq.'lsqr' &
                    .and.rake_constraint%nfields.eq.2) then
                if (rake_constraint%nrecords.eq.1) then
                    rak = rake_constraint%array(1,2)
                elseif (rake_constraint%nrecords.eq.fault%nrecords) then
                    rak = rake_constraint%array(j,2)
                else
                    call usage('!! calc_gf_los_okada_pt: incorrect number of rake '// &
                                     'constraints')
                endif
            else
                rak = 90.0d0
            endif

            call o92pt(uxp,uyp,uz,x,y,stdp,evdp,dip,rak,area,slip,vp,vs,dens)
            ! Rotate back to original x-y coordinates
            theta = datan2(uyp,uxp)
            uhor = dsqrt(uxp*uxp+uyp*uyp)
            theta = d2r*str - theta
            ux = uhor*dsin(theta)
            uy = uhor*dcos(theta)
            gf_los%array(i,j+fault%nrecords) = &
                         ux*sin(lookaz)*cos(lookinc) + uy*cos(lookaz)*cos(lookinc) - uz*sin(lookinc)
        enddo
    enddo

    if (verbosity.ge.2) then
        write(stderr,'(A)') 'calc_gf_los_okada_pt says: finished'
    endif
    if (verbosity.ge.3) then
        write(stderr,'(2A4,6A14)') 'sta','flt','losgf_ss','losgf_ds'
        do i = 1,los%nrecords
            do j = 1,fault%nrecords
                write(stderr,'(2I4,6F14.6)') i,j, &
                    gf_los%array(i,j), gf_los%array(i,j+fault%nrecords)
            enddo
        enddo
    endif
    if (verbosity.ge.2) then
        write(stderr,*)
    endif

    return
    end subroutine calc_gf_los_okada_pt

!--------------------------------------------------------------------------------------------------!
! TRIANGLE DISLOCATIONS IN AN ELASTIC HALF-SPACE
!--------------------------------------------------------------------------------------------------!

    subroutine calc_gf_disp_tri()
    use io, only: stderr, verbosity
    use variable_module, only: inversion_mode, coord_type, &
                               displacement, fault, gf_disp, halfspace, rake_constraint
    use tri_disloc, only: tri_disloc_disp, tri_center
    implicit none
    ! Local variables
    integer :: i, j, iTri
    double precision :: slip(3), slip_magnitude, rak, sta_coord(3), tri_coord(3,4)
    double precision :: sta_coord_new(3), tri_coord_new(3,4), disp(3), center(3), dist, az
    double precision :: lambda, shear_modulus, poisson
    double precision, parameter :: pi=datan(1.0d0)*4.0d0, d2r = pi/180.0d0

    if (verbosity.ge.2) then
        write(stderr,'(A)') 'calc_gf_disp_tri says: starting'
    endif

    ! Initialize rake angle
    rak = 0.0d0

    ! Unit slip for displacement GF computation
    slip_magnitude = 1.0d0

    ! Half-space array contains vp, vs, density
    shear_modulus = halfspace%array(1,2)*halfspace%array(1,2)*halfspace%array(1,3)
    lambda = halfspace%array(1,1)*halfspace%array(1,1)*halfspace%array(1,3) - 2.0d0*shear_modulus
    poisson = lambda/2.0d0/(lambda+shear_modulus)

    ! Displacement Green's function for each fault-station pair
    do i = 1,displacement%nrecords
        sta_coord(1) = displacement%array(i,1)
        sta_coord(2) = displacement%array(i,2)
        sta_coord(3) = displacement%array(i,3)
        if (verbosity.ge.1) then
            write(0,'(A1,A,I3,A)',advance='no') achar(13), &
                'calc_gf_disp_tri progress: ',100*i/displacement%nrecords,'%'
            if (i.eq.displacement%nrecords) then
                write(0,*)
            endif
        endif

        do j = 1,fault%nrecords
            tri_coord(1,1) = fault%array(j,1)
            tri_coord(2,1) = fault%array(j,2)
            tri_coord(3,1) = fault%array(j,3)
            tri_coord(1,2) = fault%array(j,4)
            tri_coord(2,2) = fault%array(j,5)
            tri_coord(3,2) = fault%array(j,6)
            tri_coord(1,3) = fault%array(j,7)
            tri_coord(2,3) = fault%array(j,8)
            tri_coord(3,3) = fault%array(j,9)
            tri_coord(:,4) = 0.0d0

            if (coord_type.eq.'geographic') then
                call tri_center(center,tri_coord(:,1),tri_coord(:,2),tri_coord(:,3))
                call ddistaz(dist,az,center(1),center(2),sta_coord(1),sta_coord(2))
                dist = dist*6.371d6
                sta_coord_new(1) = dist*dsin(az)
                sta_coord_new(2) = dist*dcos(az)
                sta_coord_new(3) = sta_coord(3)
                do iTri = 1,3
                    call ddistaz(dist,az,center(1),center(2),tri_coord(1,iTri),tri_coord(2,iTri))
                    dist = dist*6.371d6
                    tri_coord_new(1,iTri) = dist*dsin(az)
                    tri_coord_new(2,iTri) = dist*dcos(az)
                    tri_coord_new(3,iTri) = tri_coord(3,iTri)
                enddo
            else
                sta_coord_new = sta_coord
                tri_coord_new = tri_coord
            endif

            ! Check for rake constraints; if none, calculate strike-slip GFs
            if (rake_constraint%file.ne.'none'.and.inversion_mode.eq.'lsqr') then
                if (rake_constraint%nrecords.eq.1) then
                    rak = rake_constraint%array(1,1)
                elseif (rake_constraint%nrecords.eq.fault%nrecords) then
                    rak = rake_constraint%array(j,1)
                else
                    call usage('!! Error: incorrect number of rake constraints')
                endif
            else
                rak = 0.0d0
            endif
            slip(1) = -slip_magnitude*dcos(rak*d2r)   ! strike-slip, right lateral=positive
            slip(2) = -slip_magnitude*dsin(rak*d2r)   ! dip-slip, normal=positive
            slip(3) = 0.0d0                           ! tensile slip

            ! Strike-slip Green's function
            call tri_disloc_disp(disp, sta_coord_new, tri_coord_new, poisson, slip)
            ! Flip z-component
            disp(3) = -disp(3)
            gf_disp%array(i,                        j) = disp(1)
            gf_disp%array(i+1*displacement%nrecords,j) = disp(2)
            gf_disp%array(i+2*displacement%nrecords,j) = disp(3)

            ! Check for rake constraints; if none, calculate dip-slip GFs
            if (rake_constraint%file.ne.'none'.and.inversion_mode.eq.'lsqr' &
                    .and.rake_constraint%nfields.eq.2) then
                if (rake_constraint%nrecords.eq.1) then
                    rak = rake_constraint%array(1,2)
                elseif (rake_constraint%nrecords.eq.fault%nrecords) then
                    rak = rake_constraint%array(j,2)
                else
                    call usage('!! Error: incorrect number of rake constraints')
                endif
            else
                rak = 90.0d0
            endif
            slip(1) = -slip_magnitude*dcos(rak*d2r)   ! strike-slip, right lateral=positive
            slip(2) = -slip_magnitude*dsin(rak*d2r)   ! dip-slip, normal=positive
            slip(3) = 0.0d0                           ! tensile slip

            ! Dip-slip (or second rake) Green's function
            call tri_disloc_disp(disp, sta_coord_new, tri_coord_new, poisson, slip)
            ! Flip z-component
            disp(3) = -disp(3)
            gf_disp%array(i,                        j+fault%nrecords) = disp(1)
            gf_disp%array(i+1*displacement%nrecords,j+fault%nrecords) = disp(2)
            gf_disp%array(i+2*displacement%nrecords,j+fault%nrecords) = disp(3)
        enddo
    enddo

    if (verbosity.ge.2) then
        write(stderr,'(A)') 'calc_gf_disp_tri says: finished'
    endif
    if (verbosity.ge.3) then
        write(stderr,'(2A4,6A14)') 'sta','flt','dgf_ss_x','dgf_ss_y','dgf_ss_z',&
                                               'dgf_ds_x','dgf_ds_y','dgf_ds_z'
        do i = 1,displacement%nrecords
            do j = 1,fault%nrecords
                write(stderr,'(2I4,6F14.6)') i,j, &
                    gf_disp%array(i                        ,j), &
                    gf_disp%array(i+1*displacement%nrecords,j), &
                    gf_disp%array(i+2*displacement%nrecords,j), &
                    gf_disp%array(i                        ,j+fault%nrecords), &
                    gf_disp%array(i+1*displacement%nrecords,j+fault%nrecords), &
                    gf_disp%array(i+2*displacement%nrecords,j+fault%nrecords)
            enddo
        enddo
    endif
    if (verbosity.ge.2) then
        write(stderr,*)
    endif

    return
    end subroutine calc_gf_disp_tri

!--------------------------------------------------------------------------------------------------!

    subroutine calc_gf_stress_tri()
    use io, only: stderr, verbosity
    use variable_module, only: inversion_mode, coord_type, fault, rake_constraint, halfspace, &
                               gf_stress, sts_dist
    use elast
    use tri_disloc, only: tri_disloc_strain, tri_center, tri_geometry, tri_geo2cart
    implicit none
    ! Local variables
    double precision :: slip(3), slip_magnitude, rak, sta_coord(3), tri_coord(3,4), &
                        tri_coord_new(3,4), sta_coord_new(3), center(3)
    double precision :: shear_modulus, lambda, poisson
    double precision :: unit_normal(3), unit_strike(3), unit_updip(3)
    double precision :: strain(3,3), stress(3,3), traction(3), traction_components(3)
    double precision :: dist, az, pt1(3), pt2(3), pt3(3)
    double precision, parameter :: pi=datan(1.0d0)*4.0d0, d2r=pi/180.0d0
    integer :: i, j, k, prog, iTri

    if (verbosity.ge.2) then
        write(stderr,'(A)') 'calc_gf_stress_tri says: starting'
    endif
    prog = 0

    ! Initialize rake angle
    rak = 0.0d0

    ! Unit slip for GF computation
    slip_magnitude = 1.0d0

    ! Half-space array contains vp, vs, density
    shear_modulus = halfspace%array(1,2)*halfspace%array(1,2)*halfspace%array(1,3)
    lambda = halfspace%array(1,1)*halfspace%array(1,1)*halfspace%array(1,3) - 2.0d0*shear_modulus
    poisson = lambda/2.0d0/(lambda+shear_modulus)

    ! Shear stress Green's functions for each fault-fault pair
    do i = 1,fault%nrecords
        if (verbosity.ge.1) then
            write(0,'(A1,A,I3,A)',advance='no') achar(13), &
                'calc_gf_stress_tri progress: ',100*i/fault%nrecords,'%'
            if (i.eq.fault%nrecords) then
                write(0,*)
            endif
        endif
        call tri_center(sta_coord,fault%array(i,1:3),fault%array(i,4:6),fault%array(i,7:9))
        if (coord_type.eq.'cartesian') then
            call tri_geometry(unit_normal,unit_strike,unit_updip,&
                              fault%array(i,1:3),fault%array(i,4:6),fault%array(i,7:9))
        elseif (coord_type.eq.'geographic') then
            ! Triangle points: lon lat dep(m) to x y z
            call tri_geo2cart(pt1,pt2,pt3,fault%array(i,1:3),fault%array(i,4:6), &
                              fault%array(i,7:9),'m')
            call tri_geometry(unit_normal,unit_strike,unit_updip,pt1,pt2,pt3)
        endif
        ! write(0,*) 'unit_normal:',unit_normal
        ! write(0,*) 'unit_strike:',unit_strike
        ! write(0,*) 'unit_updip:',unit_updip

        do j = 1,fault%nrecords
            ! write(0,*) i,j
            tri_coord(1,1) = fault%array(j,1)
            tri_coord(2,1) = fault%array(j,2)
            tri_coord(3,1) = fault%array(j,3)
            tri_coord(1,2) = fault%array(j,4)
            tri_coord(2,2) = fault%array(j,5)
            tri_coord(3,2) = fault%array(j,6)
            tri_coord(1,3) = fault%array(j,7)
            tri_coord(2,3) = fault%array(j,8)
            tri_coord(3,3) = fault%array(j,9)
            tri_coord(:,4) = 0.0d0

            if (coord_type.eq.'geographic') then
                call tri_center(center,tri_coord(:,1),tri_coord(:,2),tri_coord(:,3))
                call ddistaz(dist,az,center(1),center(2),sta_coord(1),sta_coord(2))
                dist = dist*6.371d6
                sta_coord_new(1) = dist*dsin(az)
                sta_coord_new(2) = dist*dcos(az)
                sta_coord_new(3) = sta_coord(3)
                do iTri = 1,3
                    call ddistaz(dist,az,center(1),center(2),tri_coord(1,iTri),tri_coord(2,iTri))
                    dist = dist*6.371d6
                    tri_coord_new(1,iTri) = dist*dsin(az)
                    tri_coord_new(2,iTri) = dist*dcos(az)
                    tri_coord_new(3,iTri) = tri_coord(3,iTri)
                enddo
            else
                sta_coord_new = sta_coord
                tri_coord_new = tri_coord
            endif

            ! Set resolved tractions to zero if distance is larger than threshold
            dist = dsqrt(sta_coord_new(1)*sta_coord_new(1)+ &
                         sta_coord_new(2)*sta_coord_new(2)+ &
                        (sta_coord_new(3)-center(3))*(sta_coord(3)-center(3)))
            if (dist.gt.sts_dist) then
                gf_stress%array(i,j) = 0.0d0
                gf_stress%array(i,j+fault%nrecords) = 0.0d0
                gf_stress%array(i+fault%nrecords,j) = 0.0d0
                gf_stress%array(i+fault%nrecords,j+fault%nrecords) = 0.0d0
                cycle
            endif

            ! Avoid singular solutions with measurement point lying on fault
            do k = 1,3
                if (dabs(center(k)-sta_coord_new(k)).lt.1.0d-3) then
                    sta_coord_new(k) = center(k) + 1.0d-3
                endif
            enddo

            ! Check for rake constraints; if none, calculate strike-slip GFs
            if (rake_constraint%file.ne.'none'.and.inversion_mode.eq.'lsqr') then
                if (rake_constraint%nrecords.eq.1) then
                    rak = rake_constraint%array(1,1)
                elseif (rake_constraint%nrecords.eq.fault%nrecords) then
                    rak = rake_constraint%array(j,1)
                else
                    call usage('!! Error: incorrect number of rake constraints')
                endif
            else
                rak = 0.0d0
            endif
            slip(1) = -slip_magnitude*dcos(rak*d2r)   ! strike-slip, right lateral=positive
            slip(2) = -slip_magnitude*dsin(rak*d2r)   ! dip-slip, normal=positive
            slip(3) = 0.0d0                           ! tensile slip

            ! Shear stress Green's function produced by strike-slip or fixed-rake source
            ! write(0,*) 'sta_coord:',sta_coord_new
            ! write(0,*) 'tri_coord1:',tri_coord_new(:,1)
            ! write(0,*) 'tri_coord2:',tri_coord_new(:,2)
            ! write(0,*) 'tri_coord3:',tri_coord_new(:,3)
            ! write(0,*) 'poisson:',poisson
            ! write(0,*) 'slip:',slip
            call tri_disloc_strain(strain, sta_coord_new, tri_coord_new, poisson, slip)
            ! Flip some strain components
            strain(1,3) = -strain(1,3)
            strain(3,1) = -strain(1,3)
            strain(2,3) = -strain(2,3)
            strain(3,2) = -strain(2,3)
            ! write(0,*) strain(1,1),strain(2,2),strain(3,3)
            ! write(0,*) strain(1,2),strain(1,3),strain(2,3)
            ! Calculate tractions
            call calc_strain_to_stress(strain,halfspace%array(1,1),halfspace%array(1,2), &
                                       halfspace%array(1,3),stress)
            traction = trac_vector(stress,unit_normal)
            call calc_traction_components(traction,unit_normal,unit_strike,unit_updip, &
                                          traction_components)
            ! Load stress Green's functions array
            gf_stress%array(i               ,j               ) = -traction_components(2)
            gf_stress%array(i+fault%nrecords,j               ) = -traction_components(3)
            ! write(0,*) 'trac_ss:',-traction_components(2)
            ! write(0,*) ' trac_ds:',-traction_components(3)

            ! Shear stress Green's function produced by dip-slip source
            if (rake_constraint%file.ne.'none'.and.inversion_mode.eq.'lsqr' &
                    .and.rake_constraint%nfields.eq.2) then
                if (rake_constraint%nrecords.eq.1) then
                    rak = rake_constraint%array(1,2)
                elseif (rake_constraint%nrecords.eq.fault%nrecords) then
                    rak = rake_constraint%array(j,2)
                else
                    call usage('!! Error: incorrect number of rake constraints')
                endif
            else
                rak = 90.0d0
            endif
            slip(1) = -slip_magnitude*dcos(rak*d2r)   ! strike-slip, right lateral=positive
            slip(2) = -slip_magnitude*dsin(rak*d2r)   ! dip-slip, normal=positive
            slip(3) = 0.0d0                           ! tensile slip

            ! Shear stress Green's function produced by strike-slip or fixed-rake source
            ! write(0,*) 'sta_coord:',sta_coord_new
            ! write(0,*) 'tri_coord1:',tri_coord_new(:,1)
            ! write(0,*) 'tri_coord2:',tri_coord_new(:,2)
            ! write(0,*) 'tri_coord3:',tri_coord_new(:,3)
            ! write(0,*) 'poisson:',poisson
            ! write(0,*) 'slip:',slip
            call tri_disloc_strain(strain, sta_coord_new, tri_coord_new, poisson, slip)
            strain(1,3) = -strain(1,3)
            strain(3,1) = -strain(1,3)
            strain(2,3) = -strain(2,3)
            strain(3,2) = -strain(2,3)
            ! write(0,*) strain(1,1),strain(2,2),strain(3,3)
            ! write(0,*) strain(1,2),strain(1,3),strain(2,3)
            ! Calculate tractions
            call calc_strain_to_stress(strain,halfspace%array(1,1),halfspace%array(1,2), &
                                       halfspace%array(1,3),stress)
            traction = trac_vector(stress,unit_normal)
            call calc_traction_components(traction,unit_normal,unit_strike,unit_updip, &
                                          traction_components)
            ! Load stress Green's functions array
            gf_stress%array(i               ,j+fault%nrecords) = -traction_components(2)
            gf_stress%array(i+fault%nrecords,j+fault%nrecords) = -traction_components(3)
            ! write(0,*) 'trac_ss:',-traction_components(2)
            ! write(0,*) 'trac_ds:',-traction_components(3)
        enddo
    enddo

    if (verbosity.ge.2) then
        write(stderr,'(A)') "calc_gf_stress_tri says: finished"
    endif
    if (verbosity.ge.3) then
        write(stderr,'(2A4,4A14)') 'sta','flt','sgf_ss_ss','sgf_ss_ds','sgf_ds_ss','sgf_ds_ds'
        do i = 1,fault%nrecords
            do j = 1,fault%nrecords
                write(stderr,'(2I4,1P4E14.6)') i,j, &
                    gf_stress%array(i               ,j), &
                    gf_stress%array(i+fault%nrecords,j), &
                    gf_stress%array(i               ,j+fault%nrecords), &
                    gf_stress%array(i+fault%nrecords,j+fault%nrecords)
            enddo
        enddo
    endif
    if (verbosity.ge.2) then
        write(stderr,*)
    endif

    return
    end subroutine calc_gf_stress_tri

!--------------------------------------------------------------------------------------------------!

    subroutine calc_gf_los_tri()
    ! use io, only: stderr, verbosity
    ! use variable_module, only: inversion_mode, &
    !                            los, fault, gf_los, halfspace, rake_constraint
    implicit none
    ! ! Local variables
    ! integer :: i, j
    ! double precision :: slip, stlo, stla, stdp, evlo, evla, evdp, str, dip, rak, wid, len
    ! double precision :: dx, dy, dist, az, x, y, ux, uy, uz, uxp, uyp, theta, uhor, vp, vs, dens
    ! double precision :: lookaz, lookinc
    ! double precision, parameter :: pi=datan(1.0d0)*4.0d0, d2r = pi/180.0d0
	!
    ! if (verbosity.ge.2) then
    !     write(stderr,'(A)') 'calc_gf_los_okada_rect says: starting'
    ! endif
	!
    ! ! Unit slip for displacement GF computation
    ! slip = 1.0d0
	!
    ! ! Half-space array contains vp, vs, density
    ! vp = halfspace%array(1,1)
    ! vs = halfspace%array(1,2)
    ! dens = halfspace%array(1,3)
	!
    ! ! LOS displacement Green's function for each fault-station pair
    ! do i = 1,los%nrecords
    !     stlo = los%array(i,1)
    !     stla = los%array(i,2)
    !     stdp = los%array(i,3)
    !     lookaz   = los%array(i,5)*d2r
    !     lookinc  = los%array(i,6)*d2r
    !     do j = 1,fault%nrecords
    !         evlo = fault%array(j,1)
    !         evla = fault%array(j,2)
    !         evdp = fault%array(j,3)
    !         str  = fault%array(j,4)
    !         dip  = fault%array(j,5)
    !         wid  = fault%array(j,6)
    !         len  = fault%array(j,7)
	!
    !         ! Origin is at fault coordinate, x-axis points horizontal up-dip
    !         dx = stlo - evlo
    !         dy = stla - evla
    !         dist = dsqrt(dx*dx+dy*dy)
    !         az = datan2(dx,dy)
    !         x = dist*( dcos(az-d2r*str))
    !         y = dist*(-dsin(az-d2r*str))
	!
    !         ! Check for rake constraints; if none, calculate strike-slip GFs
    !         if (rake_constraint%file.ne.'none'.and.inversion_mode.eq.'lsqr') then
    !             if (rake_constraint%nrecords.eq.1) then
    !                 rak = rake_constraint%array(1,1)
    !             elseif (rake_constraint%nrecords.eq.fault%nrecords) then
    !                 rak = rake_constraint%array(j,1)
    !             else
    !                 call usage('!! Error: incorrect number of rake constraints')
    !             endif
    !         else
    !             rak = 0.0d0
    !         endif
	!
    !         ! Strike-slip Green's function
    !         call o92rect(uxp,uyp,uz,x,y,stdp,evdp,dip,rak,wid,len,slip,vp,vs,dens)
    !         ! Rotate back to original x-y coordinates
    !         theta = datan2(uyp,uxp)
    !         uhor = dsqrt(uxp*uxp+uyp*uyp)
    !         theta = d2r*str - theta
    !         ux = uhor*dsin(theta)
    !         uy = uhor*dcos(theta)
    !         gf_los%array(i,j) = &
    !                      ux*sin(lookaz)*cos(lookinc) + uy*cos(lookaz)*cos(lookinc) - uz*sin(lookinc)
	!
    !         ! Check for rake constraints; if none, calculate dip-slip GFs
    !         if (rake_constraint%file.ne.'none'.and.inversion_mode.eq.'lsqr' &
    !                 .and.rake_constraint%nfields.eq.2) then
    !             if (rake_constraint%nrecords.eq.1) then
    !                 rak = rake_constraint%array(1,2)
    !             elseif (rake_constraint%nrecords.eq.fault%nrecords) then
    !                 rak = rake_constraint%array(j,2)
    !             else
    !                 call usage('!! Error: incorrect number of rake constraints')
    !             endif
    !         else
    !             rak = 90.0d0
    !         endif
	!
    !         call o92rect(uxp,uyp,uz,x,y,stdp,evdp,dip,rak,wid,len,slip,vp,vs,dens)
    !         ! Rotate back to original x-y coordinates
    !         theta = datan2(uyp,uxp)
    !         uhor = dsqrt(uxp*uxp+uyp*uyp)
    !         theta = d2r*str - theta
    !         ux = uhor*dsin(theta)
    !         uy = uhor*dcos(theta)
    !         gf_los%array(i,j+fault%nrecords) = &
    !                      ux*sin(lookaz)*cos(lookinc) + uy*cos(lookaz)*cos(lookinc) - uz*sin(lookinc)
    !     enddo
    ! enddo
	!
    ! if (verbosity.ge.2) then
    !     write(stderr,'(A)') 'calc_gf_los_okada_rect says: finished'
    ! endif
    ! if (verbosity.ge.3) then
    !     write(stderr,'(2A4,6A14)') 'sta','flt','losgf_ss','losgf_ds'
    !     do i = 1,los%nrecords
    !         do j = 1,fault%nrecords
    !             write(stderr,'(2I4,6F14.6)') i,j, &
    !                 gf_los%array(i,j), gf_los%array(i,j+fault%nrecords)
    !         enddo
    !     enddo
    ! endif
    ! if (verbosity.ge.2) then
    !     write(stderr,*)
    ! endif
	!
    return
    end subroutine calc_gf_los_tri

!--------------------------------------------------------------------------------------------------!

subroutine run_inversion()
use io, only: stdout, verbosity
use variable_module, only: inversion_mode
use lsqr_module, only: invert_lsqr
use anneal_module, only: invert_anneal, invert_anneal_pseudocoupling!, invert_anneal_psc
implicit none

if (verbosity.eq.1.or.verbosity.eq.2) then
    write(stdout,*)
    write(stdout,*) 'run_inversion: starting'
endif

if (inversion_mode.eq.'lsqr') then
    call invert_lsqr()
elseif (inversion_mode.eq.'anneal') then
    call invert_anneal()
elseif (inversion_mode.eq.'anneal-psc') then
    call invert_anneal_pseudocoupling()
    ! call invert_anneal_psc()
else
    call usage('!! Error: no inversion mode named '//trim(inversion_mode))
endif

if (verbosity.eq.1.or.verbosity.eq.2) then
    write(stdout,*) 'run_inversion: finished'
endif

return
end subroutine run_inversion

!--------------------------------------------------------------------------------------------------!

subroutine free_memory()
! Deallocate arrays

use io, only: stderr, verbosity
use variable_module, only: displacement, prestress, los, fault, &
                           gf_disp, gf_stress, gf_los, &
                           smoothing, smoothing_neighbors, rake_constraint, slip_constraint, &
                           halfspace
implicit none

if (verbosity.eq.1.or.verbosity.eq.2) then
    write(stderr,*) 'free_memory: starting'
endif

if (allocated(displacement%array)) then
    deallocate(displacement%array)
endif
if (allocated(prestress%array)) then
    deallocate(prestress%array)
endif
if (allocated(los%array)) then
    deallocate(los%array)
endif
if (allocated(fault%array)) then
    deallocate(fault%array)
endif
if (allocated(gf_disp%array)) then
    deallocate(gf_disp%array)
endif
if (allocated(gf_stress%array)) then
    deallocate(gf_stress%array)
endif
if (allocated(gf_los%array)) then
    deallocate(gf_los%array)
endif
if (allocated(smoothing%intarray)) then
    deallocate(smoothing%intarray)
endif
if (allocated(smoothing_neighbors)) then
    deallocate(smoothing_neighbors)
endif
if (allocated(rake_constraint%array)) then
    deallocate(rake_constraint%array)
endif
if (allocated(slip_constraint%array)) then
    deallocate(slip_constraint%array)
endif
if (allocated(halfspace%array)) then
    deallocate(halfspace%array)
endif

if (verbosity.eq.2.or.verbosity.eq.2) then
    write(stderr,*) 'free_memory says: finished'
endif

return
end subroutine free_memory

!--------------------------------------------------------------------------------------------------!

subroutine write_solution()
use io, only: stdout, stderr, verbosity
use variable_module, only: output_file, inversion_mode, fault, fault_slip, rake_constraint, &
                           disp_misfit_file, los_misfit_file
use anneal_module, only: disp_misfit_l2norm, los_misfit_l2norm
implicit none
! Local variables
integer :: i, ounit
double precision :: slip_mag, tmp_slip_array(fault%nrecords,2)

if (verbosity.eq.1.or.verbosity.eq.2) then
    write(stderr,'(A)') 'write_solution: starting'
endif

! Print RMS misfit if specified
if (disp_misfit_file.ne.'none') then
    if (verbosity.ge.1) then
        write(stderr,'(A)') 'write_solution says: writing displacement RMS misfit to '// &
                            trim(disp_misfit_file)
    endif

    write(0,*) 'write_solution: opening misfit file'
    open(unit=81,file=disp_misfit_file,status='unknown')

    ! If rake is constrained, make an array of the correct size to use with misfit function
    if (inversion_mode.eq.'lsqr'.and.rake_constraint%file.ne.'none' &
                                                           .and.rake_constraint%nfields.eq.1) then
        write(0,*) 'write_solution: writing misfit for fixed rake'
        do i = 1,fault%nrecords
            tmp_slip_array(i,1) = fault_slip(i,1) ! Green's functions already calculated for this rake
            tmp_slip_array(i,2) = 0.0d0
        enddo
        write(81,*) disp_misfit_l2norm(tmp_slip_array)/dsqrt(dble(fault%nrecords))

    ! Otherwise, just use the misfit function directly
    else
        write(0,*) 'write_solution: writing misfit for free rake'
        write(81,*) disp_misfit_l2norm(fault_slip)/dsqrt(dble(fault%nrecords))
    endif

    write(0,*) 'write_solution: closing misfit file'
    close(81)
endif

! Line-of-sight RMS misfit
if (los_misfit_file.ne.'none') then
    if (verbosity.ge.1) then
        write(stderr,'(A)') 'write_solution says: writing LOS RMS misfit to '//trim(los_misfit_file)
    endif

    write(0,*) 'write_solution: opening misfit file'
    open(unit=81,file=los_misfit_file,status='unknown')

    ! If rake is constrained, make an array of the correct size to use with misfit function
    if (inversion_mode.eq.'lsqr'.and.rake_constraint%file.ne.'none' &
                                                           .and.rake_constraint%nfields.eq.1) then
        do i = 1,fault%nrecords
            tmp_slip_array(i,1) = fault_slip(i,1) ! Green's functions already calculated for this rake
            tmp_slip_array(i,2) = 0.0d0
        enddo
        write(81,*) los_misfit_l2norm(tmp_slip_array)/dsqrt(dble(fault%nrecords))

    ! Otherwise, just use the misfit function directly
    else
        write(81,*) los_misfit_l2norm(fault_slip)/dsqrt(dble(fault%nrecords))
    endif

    write(0,*) 'write_solution: closing misfit file'
    close(81)
endif

! Print fault slip solution
if (verbosity.ge.1) then
    write(stderr,'(A)') 'write_solution says: writing slip solution to '//trim(output_file)
endif
if (output_file.eq.'stdout') then
    ounit = stdout
else
    ounit = 99
    open (unit=ounit,file=output_file,status='unknown')
endif

do i = 1,fault%nrecords
    slip_mag = dsqrt(fault_slip(i,1)*fault_slip(i,1)+fault_slip(i,2)*fault_slip(i,2))

    if (inversion_mode.eq.'lsqr') then
        if (rake_constraint%file.eq.'none'.or.rake_constraint%nfields.eq.2) then
            if (slip_mag.lt.1.0d3) then
                write(ounit,5011) fault_slip(i,1), fault_slip(i,2)
            else
                write(ounit,5001) fault_slip(i,1), fault_slip(i,2)
            endif
        else
            slip_mag = dabs(fault_slip(i,1))
            if (slip_mag.lt.1.0d3) then
                write(ounit,5012) fault_slip(i,1)
            else
                write(ounit,5002) fault_slip(i,1)
            endif
        endif

    elseif (inversion_mode.eq.'anneal') then
        if (slip_mag.lt.1.0d3) then
            write(ounit,5011) fault_slip(i,1), fault_slip(i,2)
        else
            write(ounit,5001) fault_slip(i,1), fault_slip(i,2)
        endif

    elseif (inversion_mode.eq.'anneal-psc') then
        if (slip_mag.lt.1.0d3) then
            write(ounit,5011) fault_slip(i,1), fault_slip(i,2)
        else
            write(ounit,5001) fault_slip(i,1), fault_slip(i,2)
        endif

    else
        call usage('!! Error: frankly, I do not know how you got this far using an '//&
                         'inversion mode that does not seem to exist...')
    endif
enddo
5011 format(2F14.3)
5001 format(1P2E14.6)
5012 format(1F14.3)
5002 format(1P1E14.6)

if (allocated(fault_slip)) then
    deallocate(fault_slip)
endif

if (output_file.ne.'stdout') then
    close(ounit)
endif

if (verbosity.ge.1) then
    write(stderr,'(A)') 'write_solution says: finished'
endif

return
end subroutine write_solution

!-----------------------------------------------------------------------------__________-----------!

subroutine initialize_fltinv_variables()
use io, only: verbosity
use variable_module, only: output_file, displacement, disp_components, prestress, stress_weight, &
                           sts_dist, fault, slip_constraint, rake_constraint, &
                           gf_type, gf_disp, gf_stress, gf_los, &
                           inversion_mode, damping_constant, smoothing_constant, smoothing, &
                           coord_type, halfspace, disp_misfit_file, los_misfit_file, &
                           los, los_weight, disp_cov_file, initialize_program_data, &
                         lsqr_mode, &
                    anneal_init_mode, anneal_log_file, max_iteration, reset_iteration, &
                         temp_start, temp_minimum, cooling_factor, anneal_verbosity, &
                         anneal_init_file, prob_lock2unlock, prob_unlock2lock, &
                         mcmc_iteration, mcmc_log_file
implicit none

! Initialize program behavior variables
output_file = 'stdout'
inversion_mode = 'lsqr'
gf_type = 'none'
coord_type = 'cartesian'
disp_components = '123'
disp_misfit_file = 'none'
disp_cov_file = 'none'
los_misfit_file = 'none'

! Initialize regularization variables
damping_constant = -1.0d0
smoothing_constant = -1.0d0

! Initialize derived data variable values
call initialize_program_data(fault)
call initialize_program_data(displacement)
call initialize_program_data(prestress)
call initialize_program_data(gf_disp)
call initialize_program_data(gf_stress)
call initialize_program_data(gf_los)
call initialize_program_data(slip_constraint)
call initialize_program_data(rake_constraint)
call initialize_program_data(smoothing)
call initialize_program_data(los)

! Same type of derived data variable, but not always used
call initialize_program_data(halfspace)

! Pretty obvious, but initialize program verbosity variable
verbosity = 0

! Initialize least-squares variables
lsqr_mode = 'gels'
stress_weight = 1.0d-9
los_weight = 1.0d0

! Initialize annealing variables
anneal_init_mode = 'mean'
anneal_log_file = 'none'
max_iteration = 1000
reset_iteration = 1000000
temp_start = 2.0d0
temp_minimum = 0.00d0
cooling_factor = 0.98d0
anneal_verbosity = 0
sts_dist = 1.0d10
anneal_init_file = ''
prob_lock2unlock = 0.25d0
prob_unlock2lock = 0.10d0
mcmc_iteration = 0
mcmc_log_file = 'none'

return
end

!--------------------------------------------------------------------------------------------------!

subroutine gcmdln()
use io, only: stderr, stdout, verbosity
use variable_module, only: output_file, displacement, disp_components, prestress, stress_weight, &
                           sts_dist, fault, slip_constraint, rake_constraint, &
                           gf_type, gf_disp, gf_stress, &
                           inversion_mode, damping_constant, smoothing_constant, smoothing, &
                           coord_type, halfspace, disp_misfit_file, los_misfit_file, &
                           los, los_weight, disp_cov_file, &
                         lsqr_mode, &
                      anneal_init_mode, anneal_log_file, max_iteration, reset_iteration, &
                         temp_start, temp_minimum, cooling_factor, anneal_verbosity, &
                         anneal_control_file, anneal_init_file, prob_lock2unlock, prob_unlock2lock,&
                         mcmc_iteration,mcmc_log_file
implicit none
! Local variables
integer :: i, narg
character(len=256) :: tag
integer :: char_index

narg = command_argument_count()
if (narg.eq.0) call usage('')

i = 1
do while (i.le.narg)

    call get_command_argument(i,tag)

    ! Inversion mode
    if (trim(tag).eq.'-mode') then
        i = i + 1
        call get_command_argument(i,inversion_mode)
        if (trim(inversion_mode).eq.'anneal-psc') then
            anneal_init_mode = 'unlocked'
        endif

    ! Output options
    elseif (trim(tag).eq.'-o') then
        i = i + 1
        call get_command_argument(i,output_file)

    ! Input options
    elseif (trim(tag).eq.'-disp') then
        i = i + 1
        call get_command_argument(i,displacement%file)
    elseif (trim(tag).eq.'-disp:components') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) disp_components
    elseif (trim(tag).eq.'-disp:misfit') then
        i = i + 1
        call get_command_argument(i,disp_misfit_file)
    elseif (trim(tag).eq.'-disp:cov_file') then
        i = i + 1
        call get_command_argument(i,disp_cov_file)
    elseif (trim(tag).eq.'-los:misfit') then
        i = i + 1
        call get_command_argument(i,los_misfit_file)
    elseif (trim(tag).eq.'-prests') then
        i = i + 1
        call get_command_argument(i,prestress%file)
    elseif (trim(tag).eq.'-prests:weight') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) stress_weight
    elseif (trim(tag).eq.'-prests:dist_threshold') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) sts_dist
    elseif (trim(tag).eq.'-flt') then
        i = i + 1
        call get_command_argument(i,fault%file)
    elseif (trim(tag).eq.'-flt:rake') then
        i = i + 1
        call get_command_argument(i,rake_constraint%file)
    elseif (trim(tag).eq.'-flt:slip') then
        i = i + 1
        call get_command_argument(i,slip_constraint%file)
    elseif (trim(tag).eq.'-los') then
        i = i + 1
        call get_command_argument(i,los%file)
    elseif (trim(tag).eq.'-los:weight') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) los_weight

    ! Green's functions options
    elseif (trim(tag).eq.'-gf:model') then
        i = i + 1
        call get_command_argument(i,gf_type)
    elseif (trim(tag).eq.'-gf:disp') then
        i = i + 1
        call get_command_argument(i,gf_disp%file)
    elseif (trim(tag).eq.'-gf:stress') then
        i = i + 1
        call get_command_argument(i,gf_stress%file)

    ! Inversion options
    elseif (trim(tag).eq.'-damp') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) damping_constant
    elseif (trim(tag).eq.'-smooth') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) smoothing_constant
        i = i + 1
        call get_command_argument(i,smoothing%file)

    ! Miscellaneous options
    elseif (trim(tag).eq.'-geo') then
        coord_type = 'geographic'
    elseif (trim(tag).eq.'-haf') then
        i = i + 1
        call get_command_argument(i,halfspace%file)
        i = i + 1
        if (i.gt.narg) then
            return
        else
            call get_command_argument(i,tag)
            char_index = index(tag,'-')
            if (char_index.eq.0) then
                call get_command_argument(i,halfspace%flag)
            else
                i = i - 1
            endif
        endif
    elseif (trim(tag).eq.'-v') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) verbosity

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
    elseif (trim(tag).eq.'-anneal:max_iteration' .or. &
                               trim(tag).eq.'-anneal:it_max') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) max_iteration
    elseif (trim(tag).eq.'-anneal:reset_iteration' .or. &
                             trim(tag).eq.'-anneal:it_reset') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) reset_iteration
    elseif (trim(tag).eq.'-anneal:log_file') then
        i = i + 1
        call get_command_argument(i,anneal_log_file)
    elseif (trim(tag).eq.'-anneal:temp_start' .or. &
                            trim(tag).eq.'-anneal:temp_0') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) temp_start
    elseif (trim(tag).eq.'-anneal:temp_minimum' .or. &
                          trim(tag).eq.'-anneal:temp_min') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) temp_minimum
    elseif (trim(tag).eq.'-anneal:cooling_factor' .or. &
                                trim(tag).eq.'-anneal:cool') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) cooling_factor
    elseif (trim(tag).eq.'-anneal:p_lock2unlock' .or. &
                                trim(tag).eq.'-anneal:prob_lock2unlock') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) prob_lock2unlock
    elseif (trim(tag).eq.'-anneal:p_unlock2lock' .or. &
                                trim(tag).eq.'-anneal:prob_unlock2lock') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) prob_unlock2lock
    elseif (trim(tag).eq.'-anneal:mcmc') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) mcmc_iteration
        i = i + 1
        call get_command_argument(i,mcmc_log_file)
    elseif (trim(tag).eq.'-anneal:p_unlock2lock' .or. &
                                trim(tag).eq.'-anneal:prob_unlock2lock') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) prob_unlock2lock
    elseif (trim(tag).eq.'-anneal:verbosity') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) anneal_verbosity
    elseif (trim(tag).eq.'-anneal:control_file') then
        i = i + 1
        call get_command_argument(i,anneal_control_file)

    ! No option
    else
        call usage('!! Error: No option '//trim(tag))
    endif
    i = i + 1
enddo

! Announce verbosity mode
if (verbosity.eq.0) then
    ! fltinv running silently
elseif (verbosity.eq.1) then
    write(stdout,*) 'fltinv: verbosity level 1: major progress reports'
elseif (verbosity.eq.2) then
    write(stdout,*) 'fltinv: verbosity level 2: detailed progress reports'
elseif (verbosity.eq.3) then
    write(stdout,*) 'fltinv: verbosity level 3: parsed command line inputs'
elseif (verbosity.eq.4) then
    write(stdout,*) 'fltinv: verbosity level 4: parsed input files'
elseif (verbosity.eq.5) then
    write(stdout,*) 'fltinv: verbosity level 5: outputs'
elseif (verbosity.eq.6) then
    write(stdout,*) 'fltinv: verbosity level 6: debugging output...you asked for it...'
else
    write(stderr,*) 'gcmdln: no verbosity option ',verbosity
    call usage('')
endif

! Parsed command line values
if (verbosity.eq.3) then
    write(stdout,'("Parsed command line inputs")')
    write(stdout,'("    inversion_mode:         ",A)') trim(inversion_mode)
    write(stdout,*)
    write(stdout,'("    output_file:            ",A)') trim(output_file)
    write(stdout,*)
    write(stdout,'("    displacement%file:      ",A)') trim(displacement%file)
    write(stdout,'("    disp_components:        ",A)') trim(disp_components)
    write(stdout,'("    disp_misfit_file:       ",A)') trim(disp_misfit_file)
    write(stdout,'("    disp_cov_file:          ",A)') trim(disp_cov_file)
    write(stdout,'("    prestress%file:         ",A)') trim(prestress%file)
    write(stdout,'("    stress_weight:          ",1PE14.6)') stress_weight
    write(stdout,'("    sts_dist:               ",1PE14.6)') sts_dist
    write(stdout,'("    fault%file:             ",A)') trim(fault%file)
    write(stdout,'("    rake_constraint%file:   ",A)') trim(rake_constraint%file)
    write(stdout,'("    slip_constraint%file:   ",A)') trim(slip_constraint%file)
    write(stdout,'("    los%file:               ",A)') trim(los%file)
    write(stdout,'("    los_weight:             ",1PE14.6)') los_weight
    write(stdout,'("    los_misfit_file:        ",A)') trim(los_misfit_file)
    write(stdout,*)
    write(stdout,'("    gf_type:                ",A)') trim(gf_type)
    write(stdout,'("    gf_disp%file:           ",A)') trim(gf_disp%file)
    write(stdout,'("    gf_stress%file:         ",A)') trim(gf_stress%file)
    write(stdout,*)
    write(stdout,'("    damping_constant:       ",1PE14.6)') damping_constant
    write(stdout,'("    smoothing_constant:     ",1PE14.6)') smoothing_constant
    write(stdout,'("    smoothing_file:         ",A)') trim(smoothing%file)
    write(stdout,*)
    write(stdout,'("    coord_type:             ",A)') trim(coord_type)
    write(stdout,'("    halfspace%file:         ",A)') trim(halfspace%file)
    write(stdout,'("    halfspace%flag:         ",A)') trim(halfspace%flag)
    write(stdout,*)
    write(stdout,'("    lsqr_mode:              ",A)') trim(lsqr_mode)
    write(stdout,*)
    write(stdout,'("    anneal_init_mode:       ",A)') trim(anneal_init_mode)
    write(stdout,'("    max_iteration:          ",I14)') max_iteration
    write(stdout,'("    reset_iteration:        ",I14)') reset_iteration
    write(stdout,'("    temp_start:             ",1PE14.6)') temp_start
    write(stdout,'("    temp_minimum:           ",1PE14.6)') temp_minimum
    write(stdout,'("    anneal_log_file:        ",A)') trim(anneal_log_file)
    write(stdout,'("    cooling_factor:         ",1PE14.6)') cooling_factor
    write(stdout,'("    prob_lock2unlock:       ",1PE14.6)') prob_lock2unlock
    write(stdout,'("    prob_unlock2lock:       ",1PE14.6)') prob_lock2unlock
    write(stdout,'("    mcmc_iteration:          ",A)') mcmc_iteration
    write(stdout,'("    mcmc_log_file:          ",A)') mcmc_log_file
    write(stdout,'("    anneal_verbosity:       ",I14)') anneal_verbosity
endif

return
end subroutine gcmdln

!--------------------------------------------------------------------------------------------------!

subroutine usage(string)
!----
! Print program usage statement and exit
!----
use io, only: stderr
implicit none
character(len=*) :: string
integer :: string_length
if (string.ne.'') then
    string_length = len(string)
    write(stderr,'(A)') trim(string)
    write(stderr,'(A)')
endif
write(stderr,'(A)') 'Usage: fltinv ...options...'
write(stderr,'(A)')
write(stderr,'(A)') '-mode INVERSION_MODE         Inversion mode'
write(stderr,*)
write(stderr,'(A)') 'Output Options'
write(stderr,'(A)') '-o OUTPUT_FILE               Output file'
write(stderr,*)
write(stderr,'(A)') 'Input Options'
write(stderr,'(A)') '-disp DISP_FILE              Input displacements'
write(stderr,'(A)') '-disp:components COMPNTS     Specify displacement components'
write(stderr,'(A)') '-disp:misfit MISFIT_FILE     Output RMS misfit to displacements'
write(stderr,'(A)') '-disp:cov_file COVAR_FILE    Displacement covariances'
write(stderr,'(A)') '-prests PRESTS_FILE          Input pre-stresses'
write(stderr,'(A)') '-prests:weight WEIGHT        Stress weighting factor'
write(stderr,'(A)') '-prests:dist_threshold DIST  Set tractions to zero at distances>DIST'
write(stderr,'(A)') '-flt FAULT_FILE              Input faults'
write(stderr,'(A)') '-flt:rake RAKE_FILE          Rake angle constraints'
write(stderr,'(A)') '-flt:slip SLIP_FILE          Slip magnitude constraints'
write(stderr,'(A)') '-los LOS_FILE                Input line-of-sight displacements'
write(stderr,'(A)') '-los:weight LOS_WEIGHT       LOS observation weighting factor'
write(stderr,'(A)') '-los:misfit MISFIT_FILE      Output RMS misfit to LOS displacements'
write(stderr,*)
write(stderr,'(A)') 'Greens Functions Options'
write(stderr,'(A)') '-gf:model MODEL              Greens functions calculation model'
write(stderr,'(A)') '-gf:disp_file GF_DSP_FILE    Pre-computed displacement Greens functions'
write(stderr,'(A)') '-gf:los_file GF_LOS_FILE     Pre-computed LOS displacement Greens functions'
write(stderr,'(A)') '-gf:sts_file GF_STS_FILE     Pre-computed stress Greens functions'
write(stderr,*)
write(stderr,'(A)') 'Inversion Options'
write(stderr,'(A)') '-damp DAMP                   Damping regularization'
write(stderr,'(A)') '-smooth SMOOTH SMOOTH_FILE   Smoothing regularization'
write(stderr,*)
write(stderr,'(A)') 'Miscellaneous Options'
write(stderr,'(A)') '-geo                         Treat input coordinates as geographic'
write(stderr,'(A)') '-haf HALFSPACE_FILE [FLAG]   Elastic half-space parameters'
write(stderr,'(A)') '-v LEVEL                     Program verbosity'
write(stderr,*)
write(stderr,'(A)') 'Least-Squares Options'
write(stderr,'(A)') '-lsqr:mode MODE              Solver algorithm'
write(stderr,*)
write(stderr,'(A)') 'Simulated Annealing Options'
write(stderr,'(A)') '-anneal:init_mode MODE       Mode to initialize solution'
write(stderr,'(A)') '-anneal:max_iteration IMAX   Maximum number of iterations'
write(stderr,'(A)') '-anneal:it_reset IRESET      Reset search every IRESET iterations'
write(stderr,'(A)') '-anneal:log_file LOG_FILE    Log annealing progress'
write(stderr,'(A)') '-anneal:temp_0 START_TEMP    Starting temperature'
write(stderr,'(A)') '-anneal:temp_min MIN_TEMP    Minimum temperature'
write(stderr,'(A)') '-anneal:cool COOL_FACT       Cooling factor'
write(stderr,'(A)') '-anneal:verbosity LEVEL      Messages for annealing progress'
write(stderr,'(A)') '-anneal:p_lock2unlock P      Probability of flipping locked to unlocked'
write(stderr,'(A)') '-anneal:p_unlock2lock P      Probability of flipping unlocked to locked'
write(stderr,'(A)') '-anneal:mcmc NSTEPS LOGFILE  Run MCMC search after annealing search'
write(stderr,*)
write(stderr,'(A)') 'See man page for details'
write(stderr,*)
stop
end subroutine usage
