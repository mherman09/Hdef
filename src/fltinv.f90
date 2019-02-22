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

use io_module, only: stderr, verbosity, read_program_data_file
use variable_module, only: inversion_mode, &
                           displacement, prestress, los, fault, &
                           gf_type, gf_disp, gf_stress, gf_los, &
                           smoothing, rake_constraint, slip_constraint, &
                           halfspace, coord_type, disp_components, disp_cov_file, disp_cov_mat
use elast, only: calc_plane_unit_vectors, calc_traction, calc_traction_components
use tri_disloc_module, only: tri_geometry, tri_geo2cart

implicit none

! Local variables
integer :: i, j, n, m, ios
double precision :: stress(3,3), nor(3), str(3), upd(3), traction(3), traction_comp(3), dist, az, &
                    pt1(3), pt2(3), pt3(3), cov
character(len=256) :: line
character(len=1) :: nchar, mchar


if (verbosity.eq.2) then
    write(stderr,*) 'fltinv_readin: starting'
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
        call calc_traction(stress,nor,traction)
        call calc_traction_components(traction,nor,str,upd,traction_comp)
        prestress%array(i,1) = traction_comp(2)
        prestress%array(i,2) = traction_comp(3)
        ! write(0,*) 'trac_ss:',traction_comp(2)
        ! write(0,*) 'trac_ds:',traction_comp(3)
    enddo
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

    ! Allocate memory to calculate Green's functions
    else
        gf_disp%nrecords = 3*displacement%nrecords
        if (allocated(gf_disp%array)) then
            deallocate(gf_disp%array)
        endif
        allocate(gf_disp%array(gf_disp%nrecords,gf_disp%nfields))
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

    ! Allocate memory to calculate Green's functions
    else
        gf_stress%nrecords = 2*fault%nrecords
        if (allocated(gf_stress%array)) then
            deallocate(gf_stress%array)
        endif
        allocate(gf_stress%array(gf_stress%nrecords,gf_stress%nfields))
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

    ! Allocate memory to calculate Green's functions
    else
        gf_los%nrecords = los%nrecords
        if (allocated(gf_los%array)) then
            deallocate(gf_los%array)
        endif
        allocate(gf_los%array(gf_los%nrecords,gf_los%nfields))
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
use io_module, only: verbosity, stderr
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
use io_module, only: stderr, verbosity
use variable_module, only: inversion_mode, displacement, los, prestress, &
                           gf_type, gf_disp, gf_stress, gf_los
use gf_module, only: calc_gf_disp_okada_rect, calc_gf_stress_okada_rect, calc_gf_los_okada_rect, &
                     calc_gf_disp_okada_pt,   calc_gf_stress_okada_pt,   calc_gf_los_okada_pt, &
                     calc_gf_disp_tri,        calc_gf_stress_tri,        calc_gf_los_tri
implicit none

if (verbosity.ge.1) then
    write(stderr,'(A)') 'calc_greens_functions says: starting'
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

if (verbosity.ge.1) then
    write(stderr,'(A)') 'calc_greens_functions says: finished'
    write(stderr,*)
endif

return
end subroutine calc_greens_functions

!--------------------------------------------------------------------------------------------------!

subroutine run_inversion()
use io_module, only: stderr, verbosity
use variable_module, only: inversion_mode
use lsqr_module, only: invert_lsqr
use anneal_module, only: invert_anneal, invert_anneal_pseudocoupling
implicit none

if (verbosity.ge.1) then
    write(stderr,'(A)') 'run_inversion says: starting'
endif

if (inversion_mode.eq.'lsqr') then
    call invert_lsqr()
elseif (inversion_mode.eq.'anneal') then
    call invert_anneal()
elseif (inversion_mode.eq.'anneal-psc') then
    call invert_anneal_pseudocoupling()
else
    call usage('!! Error: no inversion mode named '//trim(inversion_mode))
endif

if (verbosity.ge.1) then
    write(stderr,'(A)') 'run_inversion says: finished'
    write(stderr,*)
endif

return
end subroutine run_inversion

!--------------------------------------------------------------------------------------------------!

subroutine free_memory()
! Deallocate arrays

use io_module, only: stderr, verbosity
use variable_module, only: displacement, prestress, los, fault, &
                           gf_disp, gf_stress, gf_los, &
                           smoothing, smoothing_neighbors, rake_constraint, slip_constraint, &
                           halfspace
implicit none

if (verbosity.eq.1) then
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

if (verbosity.eq.2) then
    write(stderr,*) 'free_memory says: finished'
endif

return
end subroutine free_memory

!--------------------------------------------------------------------------------------------------!

subroutine write_solution()
use io_module, only: stdout, stderr, verbosity
use variable_module, only: output_file, inversion_mode, fault, fault_slip, rake_constraint, &
                           disp_misfit_file, los_misfit_file
use anneal_module, only: disp_misfit_l2norm, los_misfit_l2norm
implicit none
! Local variables
integer :: i, ounit
double precision :: slip_mag, tmp_slip_array(fault%nrecords,2)

if (verbosity.ge.1) then
    write(stderr,'(A)') 'write_solution says: starting'
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
use io_module, only: verbosity, initialize_program_data
use variable_module, only: output_file, displacement, disp_components, prestress, stress_weight, &
                           sts_dist, fault, slip_constraint, rake_constraint, &
                           gf_type, gf_disp, gf_stress, gf_los, &
                           inversion_mode, damping_constant, smoothing_constant, smoothing, &
                           coord_type, halfspace, disp_misfit_file, los_misfit_file, &
                           los, los_weight, disp_cov_file
use lsqr_module, only: lsqr_mode
use anneal_module, only: anneal_init_mode, anneal_log_file, max_iteration, reset_iteration, &
                         temp_start, temp_minimum, cooling_factor, anneal_verbosity, &
                         anneal_init_file
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

return
end

!--------------------------------------------------------------------------------------------------!

subroutine gcmdln()
use io_module, only: stderr, verbosity
use variable_module, only: output_file, displacement, disp_components, prestress, stress_weight, &
                           sts_dist, fault, slip_constraint, rake_constraint, &
                           gf_type, gf_disp, gf_stress, &
                           inversion_mode, damping_constant, smoothing_constant, smoothing, &
                           coord_type, halfspace, disp_misfit_file, los_misfit_file, &
                           los, los_weight, disp_cov_file
use lsqr_module, only: lsqr_mode
use anneal_module, only: anneal_init_mode, anneal_log_file, max_iteration, reset_iteration, &
                         temp_start, temp_minimum, cooling_factor, anneal_verbosity, &
                         anneal_control_file, anneal_init_file
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

if (verbosity.ge.1) then
    write(stderr,'("fltinv verbose output turned on")')
endif
if (verbosity.ge.2) then
    write(stderr,'("Parsed command line inputs")')
    write(stderr,'("    inversion_mode:         ",A)') trim(inversion_mode)
    write(stderr,*)
    write(stderr,'("    output_file:            ",A)') trim(output_file)
    write(stderr,*)
    write(stderr,'("    displacement%file:      ",A)') trim(displacement%file)
    write(stderr,'("    disp_components:        ",A)') trim(disp_components)
    write(stderr,'("    disp_misfit_file:       ",A)') trim(disp_misfit_file)
    write(stderr,'("    disp_cov_file:          ",A)') trim(disp_cov_file)
    write(stderr,'("    prestress%file:         ",A)') trim(prestress%file)
    write(stderr,'("    stress_weight:          ",1PE14.6)') stress_weight
    write(stderr,'("    sts_dist:               ",1PE14.6)') sts_dist
    write(stderr,'("    fault%file:             ",A)') trim(fault%file)
    write(stderr,'("    rake_constraint%file:   ",A)') trim(rake_constraint%file)
    write(stderr,'("    slip_constraint%file:   ",A)') trim(slip_constraint%file)
    write(stderr,'("    los%file:               ",A)') trim(los%file)
    write(stderr,'("    los_weight:             ",1PE14.6)') los_weight
    write(stderr,'("    los_misfit_file:        ",A)') trim(los_misfit_file)
    write(stderr,*)
    write(stderr,'("    gf_type:                ",A)') trim(gf_type)
    write(stderr,'("    gf_disp%file:           ",A)') trim(gf_disp%file)
    write(stderr,'("    gf_stress%file:         ",A)') trim(gf_stress%file)
    write(stderr,*)
    write(stderr,'("    damping_constant:       ",1PE14.6)') damping_constant
    write(stderr,'("    smoothing_constant:     ",1PE14.6)') smoothing_constant
    write(stderr,'("    smoothing_file:         ",A)') trim(smoothing%file)
    write(stderr,*)
    write(stderr,'("    coord_type:             ",A)') trim(coord_type)
    write(stderr,'("    halfspace%file:         ",A)') trim(halfspace%file)
    write(stderr,'("    halfspace%flag:         ",A)') trim(halfspace%flag)
    write(stderr,*)
    write(stderr,'("    lsqr_mode:              ",A)') trim(lsqr_mode)
    write(stderr,*)
    write(stderr,'("    anneal_init_mode:       ",A)') trim(anneal_init_mode)
    write(stderr,'("    max_iteration:          ",I14)') max_iteration
    write(stderr,'("    reset_iteration:        ",I14)') reset_iteration
    write(stderr,'("    temp_start:             ",1PE14.6)') temp_start
    write(stderr,'("    temp_minimum:           ",1PE14.6)') temp_minimum
    write(stderr,'("    anneal_log_file:        ",A)') trim(anneal_log_file)
    write(stderr,'("    cooling_factor:         ",1PE14.6)') cooling_factor
    write(stderr,'("    anneal_verbosity:       ",I14)') anneal_verbosity
endif
if (verbosity.ge.1) then
    write(stderr,*)
endif

return
end subroutine gcmdln

!--------------------------------------------------------------------------------------------------!

subroutine usage(string)
!----
! Print program usage statement and exit
!----
use io_module, only: stderr
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
write(stderr,*)
write(stderr,'(A)') 'See man page for details'
write(stderr,*)
stop
end subroutine usage
