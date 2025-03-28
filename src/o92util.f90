!--------------------------------------------------------------------------------------------------!
! O92UTIL
!
! Utility for computing displacements, strains, and stresses in an elastic half-space resulting from
! point source and rectangular shear dislocations. Most of the heavy lifting is done in the Fortran
! module okada92_module.f90, which can be found in this directory.
!
! References
! Okada, Y. (1992) Internal deformation due to shear and tensile faults in a half-space. Bulletin of
! the Seismological Society of America, vol. 82, no. 2, pp. 1018-1040.
!--------------------------------------------------------------------------------------------------!

module o92util

! Fault inputs
character(len=512) :: ffm_file                    ! source faults: USGS .param format
character(len=512) :: fsp_file                    ! source faults: SRCMOD FSP format
character(len=512) :: mag_file                    ! source faults: ... mag format
character(len=512) :: flt_file                    ! source faults: ... slip wid len format
character(len=512) :: tns_file                    ! volume source: ... opening wid len
character(len=512) :: gmt_fault_file              ! for plotting with psxy -SJ
logical :: isFaultFileDefined                     ! input fault tag
logical :: isTensileFileDefined                   ! input tensile source tag
character(len=16) :: fault_type                   ! point or finite source tag
character(len=16) :: empirical_relation           ! conversion from magnitude to slip wid len
double precision :: slip_threshold                ! minimum slip to calculate (NOT USED)

! Station/target/receiver inputs
character(len=512) :: station_file                ! station locations
logical :: isStationFileDefined                   ! station location tag
double precision :: auto_depth                    ! depth of automatically generated station grid
double precision :: auto_az                       ! depth of automatically generated station grid
integer :: auto_n                                 ! number of automatically generated stations (1d)
double precision :: auto_disp_threshold           ! minimum displacement value for auto grid
character(len=1) :: auto_mode                     ! automatic station mode
integer :: n_auto_transect                        ! total number of points in initial auto transects
character(len=512) :: target_file                 ! target/receiver fault geometry
logical :: isTargetFileDefined                    ! target fault tag

! Elastic half-space inputs
character(len=512) :: halfspace_file              ! elastic half-space parameters
double precision :: poisson                       ! poisson's ratio
double precision :: lame                          ! lame parameter
double precision :: shearmod                      ! shear modulus

! Outputs
character(len=512) :: displacement_file           ! output: displacement
character(len=4) :: disp_output_mode              ! enz, amz
character(len=512) :: disp_file_save              ! temporary displacement file
character(len=512) :: strain_file                 ! output: strain tensor
character(len=512) :: rotation_file               ! output: rotation tensor
character(len=512) :: stress_file                 ! output: stress tensor
character(len=512) :: estress_file                ! output: effective (maximum) shear stress
character(len=512) :: normal_file                 ! output: normal stress
character(len=512) :: shear_file                  ! output: shear stress (resolved, maximum)
character(len=512) :: coulomb_file                ! output: coulomb stress
logical :: isOutputFileDefined                    ! output file tag
logical :: iWantDisp                              ! displacement calculation tag
logical :: iWantStrain                            ! strain calculation tag
logical :: iWantRotation                          ! rotation calculation tag
logical :: iWantStress                            ! stress calculation tag
logical :: iWantTraction                          ! traction calculation tag
character(len=16) :: coord_type                   ! cartesian-m, cartesian-km, geographic
logical :: iWantProg                              ! progress indicator tag

! Program variables
integer :: nfaults                                ! Number of fault sources
integer :: ntensile                               ! Number of tensile sources
integer :: nstations                              ! Number of stations/targets/receivers
double precision, allocatable :: faults(:,:)      ! Fault parameter array
double precision, allocatable :: tensile(:,:)     ! Tensile source parameter array
double precision, allocatable :: stations(:,:)    ! Station location array
character(len=32), allocatable :: sta_char(:)     ! Station character (for header/comment lines)
logical, allocatable :: staSkipLine(:)            ! Lines to skip
double precision, allocatable :: targets(:,:)     ! Target/receiver geometry array
double precision :: centroid(3)                   ! Moment weighted mean position of EQ
logical :: keepAllLines                           ! Default: false (ignore blank/commented lines)

! Parallel options
logical :: isParallel
integer :: nthreads

! Debugging
character(len=32) :: debug_mode                   ! Verbose output for debugging

end module

!==================================================================================================!

program main
!----
! Run program o92util:
!     - Parse command line options
!     - Check inputs/outputs
!     - Read inputs
!     - Calculate deformation
!----
use io, only: verbosity, stdout
use o92util, only: isFaultFileDefined, &
                   isTensileFileDefined, &
                   isStationFileDefined, &
                   isTargetFileDefined, &
                   iWantTraction, &
                   isOutputFileDefined, &
                   auto_mode

implicit none


! Parse the command line
call gcmdln()

! Check that input and output files are defined
if (.not.isOutputFileDefined) then
    call usage('o92util: no output defined (usage:output)')
endif
if (.not.isStationFileDefined) then
    call usage('o92util: no station file defined (usage:station)')
endif
if (iWantTraction.and..not.isTargetFileDefined) then
    call usage('o92util: no target file defined (usage:station)')
endif
if (.not.isFaultFileDefined.and..not.isTensileFileDefined) then
    call usage('o92util: no fault or tensile source file defined (usage:input)')
endif


! Here we go!
if (verbosity.ge.1) then
    write(stdout,*) 'o92util: starting'
    write(stdout,*)
endif


! Define elasticity parameters for half-space
call read_halfspace()


! Read source files - usually faults, but can do tensile sources as well
call read_faults()
call read_tensile()


! Define limits of automatic grid, if using auto mode
if (auto_mode.eq.'h'.or.auto_mode.eq.'v') then
    call auto_stations()
endif


! Read input station locations and target fault geometries
call read_stations()
call read_targets()


! Calculate displacements, strains, or stresses
call calc_deformation()


! Generate automatic station grid and recalculate displacements, strains, or stresses
if (auto_mode.eq.'h'.or.auto_mode.eq.'v') then
    call update_auto_stations()
    call calc_deformation()
    if (auto_mode.eq.'v') then
        call project_deformation()
    endif
endif


! Congratulations! You have completed the program!
if (verbosity.ge.1) then
    write(stdout,*) 'o92util: finished'
    write(stdout,*)
endif


end


!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!------------------------------------------- INPUTS -----------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

subroutine read_halfspace()
!----
! Read half-space parameters with the subroutine read_halfspace_file() in the elast module.
!----

use io, only: verbosity, stdout
use elast, only: read_halfspace_file

use o92util, only: halfspace_file, &
                   poisson, &
                   lame, &
                   shearmod, &
                   debug_mode

implicit none

integer :: ierr


if (verbosity.ge.1) then
    write(stdout,*) 'read_halfspace: starting'
endif

! Read elastic moduli
call read_halfspace_file(halfspace_file,poisson,shearmod,lame,ierr)
if (ierr.ne.0) then
    call usage('')
endif


if (debug_mode.eq.'read_halfspace'.or.debug_mode.eq.'all') then
    write(stdout,'(X,"(DEBUG)",X,A,1PE12.4)') 'read_halfspace: shear_modulus=  ',shearmod
    write(stdout,'(X,"(DEBUG)",X,A,F12.6)')   'read_halfspace: poissons_ratio= ',poisson
    write(stdout,'(X,"(DEBUG)",X,A,1PE12.4)') 'read_halfspace: lame_parameter= ',lame
endif
if (verbosity.ge.2) then
    write(stdout,'(X,A,1PE12.4)') 'read_halfspace: shear_modulus=  ',shearmod
    write(stdout,'(X,A,F12.6)')   'read_halfspace: poissons_ratio= ',poisson
    write(stdout,'(X,A,1PE12.4)') 'read_halfspace: lame_parameter= ',lame
endif
if (verbosity.ge.1) then
    write(stdout,*) 'read_halfspace: finished'
    write(stdout,*)
endif

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine read_faults()
!----
! Read the input rectangular or point fault source data in one of the following formats:
!
!     mag: lon lat dep str dip rak mag
!     flt: lon lat dep str dip rak slip wid len
!     ffm: U.S. Geological Survey .param format
!     fsp: SRCMOD FSP format
!
! using subroutines from the ffm module.
!----

use io, only: verbosity, stdout, stderr, line_count
use ffm, only: read_usgs_param, read_srcmod_fsp, read_mag, read_flt, ffm_data
use eq, only: empirical, sdr2ter, mag2mom

use o92util, only: ffm_file, &
                   fsp_file, &
                   mag_file, &
                   flt_file, &
                   gmt_fault_file, &
                   isFaultFileDefined, &
                   empirical_relation, &
                   slip_threshold, &
                   shearmod, &
                   coord_type, &
                   nfaults, &
                   faults, &
                   debug_mode

implicit none

! Local variables
integer :: ierr, i, j
logical :: areFaultsRead
type(ffm_data) :: param_data, fsp_data, mag_data, flt_data
character(len=8) :: fault_type
double precision :: fth, fss, fno, mom, area


if (verbosity.ge.1) then
    write(stdout,*) 'read_faults: starting'
endif


! Initialize variables
areFaultsRead = .false.
nfaults = 0


! Check whether faults need to be read
if(.not.isFaultFileDefined) then
    if (verbosity.ge.1) then
        write(stdout,*) 'read_faults: not using fault sources; finished'
        write(stdout,*)
    endif
    return
endif


! Read faults in one of available formats: USGS PARAM, SRCMOD FSP, MAG, or FLT

! Read faults from USGS param file
if (ffm_file.ne.'') then
    call read_usgs_param(ffm_file,param_data,ierr)
    if (ierr.eq.0) then
        areFaultsRead = .true.
    else
        call usage('read_faults: error reading fault file in USGS param format (usage:input)')
    endif
    nfaults = nfaults + param_data%nflt

    if (debug_mode.eq.'read_faults'.or.debug_mode.eq.'all') then
        write(stdout,*) '(DEBUG) read_faults: finished read_usgs_param'
        write(stdout,*) '(DEBUG)   #       lon       lat       dep     str     dip     rak      ',&
                        'slip       wid       len'
        do i = 1,param_data%nflt
            write(stdout,4316) i,(param_data%subflt(i,j),j=1,2), param_data%subflt(i,3)/1.0d3, &
                               (param_data%subflt(i,j),j=4,7),(param_data%subflt(i,j)/1.0d3,j=8,9)
        enddo
        4316 format(X,I4,2F10.3,F10.3,3F8.1,F10.4,2F10.3)
    endif
endif


! Read faults from SRCMOD FSP file
if (fsp_file.ne.'') then
    call read_srcmod_fsp(fsp_file,fsp_data,ierr)
    if (ierr.eq.0) then
        areFaultsRead = .true.
    else
        call usage('read_faults: error reading fault file in SRCMOD FSP format (usage:input)')
    endif
    nfaults = nfaults + fsp_data%nflt

    if (debug_mode.eq.'read_faults'.or.debug_mode.eq.'all') then
        write(stdout,*) '(DEBUG) read_faults: finished read_srcmod_fsp'
        write(stdout,*) '(DEBUG)   #       lon       lat       dep     str     dip     rak      ',&
                        'slip       wid       len'
        do i = 1,fsp_data%nflt
            write(stdout,4317) i,(fsp_data%subflt(i,j),j=1,2), fsp_data%subflt(i,3)/1.0d3, &
                               (fsp_data%subflt(i,j),j=4,7),(fsp_data%subflt(i,j)/1.0d3,j=8,9)
        enddo
        4317 format(X,I4,2F10.3,F10.3,3F8.1,F10.4,2F10.3)
    endif
endif


! Read faults in GMT psmeca -Sa format
if (mag_file.ne.'') then
    call read_mag(mag_file,mag_data,ierr)
    if (ierr.eq.0) then
        areFaultsRead = .true.
    else
        call usage('read_faults: error reading fault file in GMT psmeca -Sa format (usage:input)')
    endif

    ! Calculate fault slip, width, and length from magnitude
    do i = 1,mag_data%nflt

        ! Determine the type of focal mechanism for empirical relation using ternary classification
        ! of Frohlich (1992), Triangle diagrams: ternary graphs to display similarity and diversity
        ! of earthquake focal mechanisms.
        call sdr2ter(mag_data%subflt(i,4),mag_data%subflt(i,5),mag_data%subflt(i,6),fth,fss,fno)
        if (fth.gt.0.60d0) then
            fault_type = 'th'
        elseif (fss.gt.0.60d0) then
            fault_type = 'ss'
        elseif (fno.gt.0.60d0) then
            fault_type = 'no'
        else
            fault_type = 'ot'
        endif
        if (verbosity.ge.3) then
            write(stdout,*) 'read_faults: strike=',mag_data%subflt(i,4)
            write(stdout,*) 'read_faults: dip=   ',mag_data%subflt(i,5)
            write(stdout,*) 'read_faults: rake=  ',mag_data%subflt(i,6)
            write(stdout,*) 'read_faults: using fault_type="',trim(fault_type),'" with empirical ',&
                            'relation="',trim(empirical_relation),'"'
        endif

        ! Calculate fault width (down-dip) and length (along-strike) from empirical relation
        call empirical(mag_data%subflt(i,7),mag_data%subflt(i,8),mag_data%subflt(i,9),&
                       empirical_relation,fault_type)

        ! Convert fault dimensions to SI units
        mag_data%subflt(i,8) = mag_data%subflt(i,8)*1.0d3 ! km->m
        mag_data%subflt(i,9) = mag_data%subflt(i,9)*1.0d3 ! km->m

        ! Calculate slip from seismic moment, area, and shear_modulus
        call mag2mom(mag_data%subflt(i,7),mom)
        area = mag_data%subflt(i,8)*mag_data%subflt(i,9)
        mag_data%subflt(i,7) = mom/(shearmod*area)
    enddo

    nfaults = nfaults + mag_data%nflt

    ! To plot faults with psxy -SJ, write to the file here (lon lat str len proj_wid)
    if (gmt_fault_file.ne.'') then
        open(unit=6213,file=gmt_fault_file,status='unknown')
        do i = 1,mag_data%nflt
            write(6213,*) mag_data%subflt(i,1), &
                          mag_data%subflt(i,2), &
                          mag_data%subflt(i,4), &
                          mag_data%subflt(i,9)/1.0d3,&
                          mag_data%subflt(i,8)/1.0d3*cos(mag_data%subflt(i,5)*0.0174533d0)
        enddo
        close(6213)
    endif

    if (debug_mode.eq.'read_faults'.or.debug_mode.eq.'all') then
        write(stdout,*) '(DEBUG) read_faults: finished read_mag'
        write(stdout,*) '(DEBUG)   #       lon       lat       dep     str     dip     rak      ',&
                        'slip       wid       len'
        do i = 1,mag_data%nflt
            write(stdout,4318) i,(mag_data%subflt(i,j),j=1,2), mag_data%subflt(i,3)/1.0d3, &
                               (mag_data%subflt(i,j),j=4,7),(mag_data%subflt(i,j)/1.0d3,j=8,9)
        enddo
        4318 format(X,'(DEBUG)',I4,2F10.3,F10.3,3F8.1,F10.4,2F10.3)
    endif
endif


! Read faults in GMT psmeca -Sa format, except with slip wid len in place of magnitude
if (flt_file.ne.'') then
    call read_flt(flt_file,flt_data,ierr)
    if (ierr.eq.0) then
        areFaultsRead = .true.
    else
        call usage('read_faults: error reading fault file in "lo la dp str dip rak slip wid len" '//&
                   'format (usage:input)')
    endif
    nfaults = nfaults + flt_data%nflt

    if (debug_mode.eq.'read_faults'.or.debug_mode.eq.'all') then
        write(stdout,*) '(DEBUG) read_faults: finished read_flt'
        write(stdout,*) '(DEBUG)   #       lon       lat       dep     str     dip     rak      ',&
                        'slip       wid       len'
        do i = 1,flt_data%nflt
            write(stdout,4319) i,(flt_data%subflt(i,j),j=1,2), flt_data%subflt(i,3)/1.0d3, &
                               (flt_data%subflt(i,j),j=4,7), (flt_data%subflt(i,j)/1.0d3,j=8,9)
        enddo
        4319 format(X,'(DEBUG)',I4,2F10.3,F10.3,3F8.1,F10.4,2F10.3)
    endif
endif


! Sanity check: make sure something has been read before proceeding
if (.not.areFaultsRead) then
    write(stderr,*) 'read_faults: no faults were read'
    call usage('check files specified by -mag, -flt, -ffm, or -fsp  (usage:input)')
endif


! Allocate memory for master fault array
if (allocated(faults)) then
    deallocate(faults)
endif
allocate(faults(nfaults,9))


! Load the master array with the values read from files
nfaults = 0

if (fsp_file.ne.''.and.fsp_data%nflt.gt.0) then
    faults(nfaults+1:nfaults+fsp_data%nflt,:) = fsp_data%subflt(1:fsp_data%nflt,1:9)
    nfaults = nfaults + fsp_data%nflt
endif

if (ffm_file.ne.''.and.param_data%nflt.gt.0) then
    faults(nfaults+1:nfaults+param_data%nflt,:) = param_data%subflt(1:param_data%nflt,1:9)
    nfaults = nfaults + param_data%nflt
endif

if (flt_file.ne.''.and.flt_data%nflt.gt.0) then
    faults(nfaults+1:nfaults+flt_data%nflt,:) = flt_data%subflt
    nfaults = nfaults + flt_data%nflt
endif

if (mag_file.ne.''.and.mag_data%nflt.gt.0) then
    faults(nfaults+1:nfaults+mag_data%nflt,:) = mag_data%subflt
    nfaults = nfaults + mag_data%nflt
endif


! Set fault segments with low slip magnitude to zero
! The meaning of slip_threshold depends on whether it is positive or negative
if (slip_threshold.lt.0.0d0) then
    ! Negative: any fault with slip less than abs(slip_threshold) is set to zero
    slip_threshold = abs(slip_threshold)
elseif (slip_threshold.gt.0.0d0) then
    ! Positive: any fault with slip less than slip_threshold*max(slip) is set to zero
    slip_threshold = slip_threshold*maxval(faults(:,7))
endif
if (verbosity.ge.3) then
    write(stdout,'(X,A,F10.3)') 'read_faults: slip_threshold=',slip_threshold
endif
do i = 1,nfaults
    if (abs(faults(i,7)).lt.slip_threshold) then
        faults(i,7) = 0.0d0
    endif
enddo


! Check coordinates
if (coord_type.eq.'geographic'.and. &
        (maxval(abs(faults(:,1))).gt.360.0d0 .or. maxval(abs(faults(:,2))).gt.90.0d0)) then
    call usage('read_faults: found fault coordinates outside geographic range; '//&
                    'did you mean to use the -xy flag?')
endif


if (verbosity.ge.2) then
    write(stdout,*) 'read_faults: read',nfaults,' faults'
endif
if (verbosity.ge.3) then
    if (coord_type.eq.'geographic') then
        write(stdout,*) '   #         Lon         Lat   Dep(km)     Str     Dip     Rak   '//&
                        'Slip(m)   Wid(km)   Len(km)'
    elseif (coord_type.eq.'cartesian') then
        write(stdout,*) '   #       X(km)       Y(km)   Dep(km)     Str     Dip     Rak   '//&
                        'Slip(m)   Wid(km)   Len(km)'
    endif
    do i = 1,nfaults
        if (coord_type.eq.'geographic') then
            write(stdout,1101) i,faults(i,1:2),faults(i,3)/1.0d3,faults(i,4:7),faults(i,8)/1.0d3, &
                               faults(i,9)/1.0d3
        elseif (coord_type.eq.'cartesian') then
            write(stdout,1101) i,faults(i,1)/1.0d3,faults(i,2)/1.0d3,faults(i,3)/1.0d3, &
                               faults(i,4:7),faults(i,8)/1.0d3,faults(i,9)/1.0d3
        endif
        1101 format(X,I4,2F12.4,F10.3,3F8.1,F10.3,3F10.3)
    enddo
endif
if (verbosity.ge.1) then
    write(stdout,*) 'read_faults: finished'
    write(stdout,*)
endif


return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine read_tensile()
!----
! Read the input rectangular tensile source data in the format:
!     lon lat dep str dip crack wid len
!----

use io, only: verbosity, stdout, stderr, line_count
use ffm, only: read_tensile_rect, ffm_data
use eq, only: empirical, sdr2ter, mag2mom

use o92util, only: tns_file, &
                   isTensileFileDefined, &
                   coord_type, &
                   ntensile, &
                   tensile, &
                   debug_mode

implicit none

! Local variables
integer :: i, j, ierr
logical :: areFaultsRead
type(ffm_data) :: tns_data


if (verbosity.ge.1) then
    write(stdout,*) 'read_tensile: starting'
endif


! Initialize I/O status
ierr = 0


! Initialize variables
areFaultsRead = .false.
ntensile = 0
if(.not.isTensileFileDefined) then
    if (verbosity.ge.1) then
        write(stdout,*) 'read_tensile: not using tensile sources; finished'
        write(stdout,*)
    endif
    return
endif


! Read faults in "lon lat dep str dip crack wid len" format
if (tns_file.ne.'') then
    call read_tensile_rect(tns_file,tns_data,ierr)
    if (ierr.eq.0) then
        areFaultsRead = .true.
    endif
    ntensile = ntensile + tns_data%nflt

    if (debug_mode.eq.'read_tensile'.or.debug_mode.eq.'all') then
        write(stdout,*) '(DEBUG) read_tensile: finished read_tensile_rect'
        write(stdout,*) '(DEBUG)   #       lon       lat       dep     str     dip     ',&
                        'crack       wid       len'
        do i = 1,tns_data%nflt
            write(stdout,4319) i,(tns_data%subflt(i,j),j=1,2), tns_data%subflt(i,3)/1.0d3, &
                               (tns_data%subflt(i,j),j=4,6), (tns_data%subflt(i,j)/1.0d3,j=8,9)
        enddo
        4319 format(X,'(DEBUG)',I4,2F10.3,F10.3,2F8.1,F10.4,2F10.3)
    endif
endif


! Sanity check: make sure something has been read before proceeding
if (.not.areFaultsRead) then
    write(stderr,*) 'read_tensile: no sources were read'
    call usage('check files specified by -tns  (usage:input)')
endif


! Allocate memory for master fault array
if (allocated(tensile)) then
    deallocate(tensile)
endif
allocate(tensile(ntensile,8))


! Load the master array with the values read from files
ntensile = 0

if (tns_file.ne.''.and.tns_data%nflt.gt.0) then
    tensile(ntensile+1:ntensile+tns_data%nflt,:) = tns_data%subflt
    ntensile = ntensile + tns_data%nflt
endif


! Check coordinates
if (coord_type.eq.'geographic'.and. &
        (maxval(abs(tensile(:,1))).gt.360.0d0 .or. maxval(abs(tensile(:,2))).gt.90.0d0)) then
    write(stderr,*) 'read_tensile: found source coordinates outside geographic range; ',&
                    'did you mean to use the -xy flag?'
endif


if (verbosity.ge.2) then
    write(stdout,*) 'read_tensile: read',ntensile,' tensile sources'
endif
if (verbosity.ge.3) then
    if (coord_type.eq.'geographic') then
        write(stdout,*) '   #         Lon         Lat   Dep(km)     Str     Dip  '//&
                        'Crack(m)     Wid(km)     Len(km)'
    elseif (coord_type.eq.'cartesian') then
        write(stdout,*) '   #       X(km)       Y(km)   Dep(km)     Str     Dip  '//&
                        'Crack(m)     Wid(km)     Len(km)'
    endif
    do i = 1,ntensile
        if (coord_type.eq.'geographic') then
            write(stdout,1111) i,tensile(i,1:2),tensile(i,3)/1.0d3,tensile(i,4:6), &
                               tensile(i,7)/1.0d3,tensile(i,8)/1.0d3
        elseif (coord_type.eq.'cartesian') then
            write(stdout,1111) i,tensile(i,1)/1.0d3,tensile(i,2)/1.0d3,tensile(i,3)/1.0d3, &
                               tensile(i,4:6),tensile(i,7)/1.0d3,tensile(i,8)/1.0d3
        endif
        1111 format(X,I4,2F12.4,F10.3,2F8.1,F10.3,3F12.3)
    enddo
endif
if (verbosity.ge.1) then
    write(stdout,*) 'read_tensile: finished'
    write(stdout,*)
endif

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine read_stations()
!----
! Read the input station file, in format:
!
!     x1 y1 z1
!
! where (x1,y1) are geographic or cartesian coordinates (depending on coord_type), and z1 is depth
! in km, positive down.
!----

use io, only: verbosity, stdout, stderr, line_count, fileExists

use o92util, only: station_file, &
                   nstations, &
                   stations, &
                   sta_char, &
                   staSkipLine, &
                   keepAllLines, &
                   isStationFileDefined, &
                   coord_type, &
                   auto_mode, &
                   debug_mode

implicit none

! Local variables
integer :: i, j, ios
character(len=512) :: input_line
logical :: foundLineToSkip


if (verbosity.ge.1) then
    write(stdout,*) 'read_stations: starting'
endif


! Initialize I/O status
ios = 0
foundLineToSkip = .false.


! Check that station file is defined and the file exists
if (.not.isStationFileDefined) then
    call usage('read_stations: no station file defined (usage:station)')
endif
if (.not.fileExists(station_file)) then
    call usage('read_stations: no station file named '//trim(station_file)//' found (usage:station)')
endif


! Count number of stations and allocate memory to station array (x y z)
nstations = line_count(station_file)
allocate(stations(nstations,3))
allocate(sta_char(nstations))
allocate(staSkipLine(nstations))
stations = 0.0d0
sta_char = ''
staSkipLine = .false.


! Read stations
open(unit=13,file=station_file,status='old')
do i = 1,nstations
    read(13,'(A)') input_line
    input_line = adjustl(input_line)
    sta_char(i) = input_line(1:32)
    if (sta_char(i)(1:1).eq.'#'.or.sta_char(i)(1:1).eq.'>') then
        ! This line is a comment (#) or a segment header (>)
        staSkipLine(i) = .true.
        if (.not.foundLineToSkip) then
            write(*,*) 'o92util: found a line to skip in the input station file starting with '//&
                       sta_char(i)(1:1)
            if (keepAllLines) then
                write(*,*) 'o92util: keeping this line in the output file'
            else
                write(*,*) 'o92util: ignoring this line and omitting from the output file'
                write(*,*) 'To keep these lines in the output file, use "--keep-all-lines" option'
            endif
            foundLineToSkip = .true.
        endif
        cycle
    elseif (sta_char(i).eq.'') then
        ! This line is empty
        if (.not.foundLineToSkip) then
            write(*,*) 'o92util: found a blank line'
            if (keepAllLines) then
                write(*,*) 'o92util: keeping this line in the output file'
            else
                write(*,*) 'o92util: ignoring this line and omitting from the output file'
                write(*,*) 'To keep these lines in the output file, use "--keep-all-lines" option'
            endif
            foundLineToSkip = .true.
        endif
        foundLineToSkip = .true.
        cycle
    endif
    read(input_line,*,iostat=ios,err=1003,end=1004) (stations(i,j),j=1,3)
    ! print *,i,trim(input_line)
enddo
close(13)

! Error messages
1003 if (ios.ne.0) then
    write(stderr,*) 'read_stations: read error'
    call usage('offending line: '//trim(input_line)//' (usage:station)')
endif
1004 if (ios.ne.0) then
    write(stderr,*) 'read_stations: line end'
    call usage('offending line: '//trim(input_line)//' (usage:station)')
endif

if (debug_mode.eq.'read_stations'.or.debug_mode.eq.'all') then
    write(stdout,*) '(DEBUG) read_stations:'
    write(stdout,*) '(DEBUG)   #       lon       lat       dep'
    do i = 1,nstations
        write(stdout,4686) i,(stations(i,j),j=1,3)
    enddo
    4686 format(X,'(DEBUG)',I4,3F10.3)
endif


! Check coordinates
if (coord_type.eq.'geographic'.and. &
        (maxval(abs(stations(:,1))).gt.360.0d0 .or. maxval(abs(stations(:,2))).gt.90.0d0)) then
    call usage('read_stations: found station coordinates outside geographic range; '//&
                    'did you mean to use the -xy flag?')
endif


if (verbosity.ge.2) then
    write(stdout,*) 'read_stations: read',nstations,' stations'
endif
if (verbosity.ge.3) then
    if (auto_mode.eq.'') then
        if (coord_type.eq.'geographic') then
            write(stdout,*) '   #         Lon         Lat   Dep(km)'
        elseif (coord_type.eq.'cartesian') then
            write(stdout,*) '   #       X(km)       Y(km)   Dep(km)'
        endif
        do i = 1,nstations
            if (coord_type.eq.'geographic') then
                write(stdout,1211) i,stations(i,1:3)
            elseif (coord_type.eq.'cartesian') then
                write(stdout,1211) i,stations(i,1:3)
            endif
            1211 format(X,I4,2F12.4,F10.3)
        enddo
    endif
endif
if (auto_mode.ne.'') then
    write(stdout,*) 'read_stations: auto_station initial transect endpoints:'
    write(stdout,'(4X,3F10.3)') stations(1,:)
    write(stdout,'(4X,3F10.3)') stations(1001,:)
    write(stdout,'(4X,3F10.3)') stations(1002,:)
    write(stdout,'(4X,3F10.3)') stations(2002,:)
endif
if (verbosity.ge.1) then
    write(stdout,*) 'read_stations: finished'
    write(stdout,*)
endif


return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine read_targets()
!----
! Read the input target fault geometry file, in format:
!
!     str dip rak fric
!
! where str, dip, and rak are in degrees and fric is the effective coefficient of friction.
!----

use io, only: verbosity, stdout, stderr, line_count, fileExists

use o92util, only: target_file, &
                   targets, &
                   isTargetFileDefined, &
                   iWantTraction, &
                   nstations, &
                   debug_mode

implicit none


! Local variables
integer :: i, j, ios, ntargets
character(len=512) :: input_line


if (verbosity.ge.1) then
    write(stdout,*) 'read_targets: starting'
endif


! Initialize I/O status
ios = 0


! Check that tractions need to be resolved, target fault file is defined, and the file exists
if (.not.iWantTraction) then
    if (verbosity.ge.1) then
        write(stdout,*) 'read_targets: target geometry not required; finished'
        write(stdout,*)
    endif
    return
else
    if (.not.isTargetFileDefined) then
        call usage('read_targets: no target geometry file defined (usage:station)')
    endif
endif
if (.not.fileExists(target_file)) then
    call usage('read_targets: no target geometry file named '//trim(target_file)//' found '// &
               '(usage:station)')
endif


! Count number of target geometries and allocate memory to target array (str dip rak fric)
ntargets = line_count(target_file)
if (ntargets.ne.1.and.ntargets.ne.nstations) then
    call usage('read_targets: number of target geometries must be 1 or nstations (usage:station)')
endif
allocate(targets(nstations,4))


! Read target fault geometries
open(unit=14,file=target_file,status='old')
do i = 1,ntargets
    read(14,'(A)') input_line
    read(input_line,*,iostat=ios,err=1005,end=1006) (targets(i,j),j=1,4)
enddo
close(14)

! Error messages
1005 if (ios.ne.0) then
    write(stderr,*) 'read_targets: read error'
    call usage('offending line: '//trim(input_line))
endif
1006 if (ios.ne.0) then
    write(stderr,*) 'read_targets: line end'
    call usage('offending line: '//trim(input_line))
endif


! Fill target array if only one is specified
if (ntargets.eq.1) then
    do i = 2,nstations
        targets(i,:) = targets(1,:)
    enddo
endif
if (verbosity.ge.2) then
    write(stdout,*) 'read_targets: read',ntargets,' target fault geometries'
endif

if (debug_mode.eq.'read_targets'.or.debug_mode.eq.'all') then
    write(stdout,*) '(DEBUG) read_targets:'
    write(stdout,*) '(DEBUG)   #       str       dip       rak      fric'
    do i = 1,nstations
        write(stdout,4687) i,(targets(i,j),j=1,4)
    enddo
    4687 format(X,'(DEBUG)',I4,4F10.3)
endif


if (verbosity.ge.3) then
    write(stdout,*) '   #       Str       Dip       Rak      Fric'
    do i = 1,ntargets
        write(stdout,1301) i,targets(i,1:4)
        1301 format(X,I4,3F10.1,1F10.3)
    enddo
endif
if (verbosity.ge.1) then
    write(stdout,*) 'read_targets: finished'
    write(stdout,*)
endif

return
end subroutine


!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!----------------------------------------- OUTPUTS ------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

subroutine calc_deformation()
!----
! Calculate displacement, strain, stress, and resolved tractions at station locations, given fault
! and tensile sources and target geometries (if needed for tractions).
!
! The main computational subroutines are from the okada92 module.
!----

use io, only: stdout, stderr, verbosity, progress_indicator
use trig, only: d2r, r2d
use algebra, only: rotate_vector_angle_axis, rotate_matrix_angle_axis
use earth, only: radius_earth_m
use elast, only: strain2stress, stress2traction, max_shear_stress, traction_components
use eq, only: sdr2sv
use geom, only: lola2distaz, strdip2normal
use okada92, only: o92_pt_disp, &
                   o92_rect_disp, &
                   o92_pt_strain, &
                   o92_pt_rotation, &
                   o92_rect_strain, &
                   o92_rect_rotation, &
                   depthWarning

use o92util, only: iWantDisp, &
                   iWantStrain, &
                   iWantRotation, &
                   iWantStress, &
                   iWantTraction, &
                   displacement_file, &
                   disp_output_mode, &
                   strain_file, &
                   rotation_file, &
                   stress_file, &
                   estress_file, &
                   normal_file, &
                   shear_file, &
                   coulomb_file, &
                   keepAllLines, &
                   lame, &
                   shearmod, &
                   fault_type, &
                   coord_type, &
                   nfaults, &
                   ntensile, &
                   nstations, &
                   faults, &
                   tensile, &
                   stations, &
                   sta_char, &
                   staSkipLine, &
                   targets, &
                   iWantProg, &
                   debug_mode, &
                   isParallel, &
                   nthreads

implicit none

! Local variables
integer :: ierr, iSta, iFlt, file_unit
logical :: isThisUnitOpen, coordTypeWarning
double precision :: evlo, evla, evdp, str, dip, rak, slip_mag, wid, len, slip(3), mom(4)
double precision :: sta_coord(3), dist, az, warn_dist, test_dist
double precision :: disp(3), disptmp(3), stn(3,3), stntmp(3,3), sts(3,3), ests, trac(3), tshr, &
                    tshrmx, tnor, coul, rot(3,3), rottmp(3,3)
double precision :: nvec(3), svec(3), tstr, tupd


#ifdef USE_OPENMP
integer :: OMP_GET_MAX_THREADS
#endif


if (verbosity.ge.1) then
    write(stdout,*) 'calc_deformation: starting'
endif


! Initialize I/O status
ierr = 0


! Check which calculations are needed and open output files if requested

! Displacement
if (iWantDisp) then
    if (verbosity.ge.2) then
        write(stdout,*) 'calc_deformation: calculating displacements'
    endif
    if (displacement_file.ne.'') then
        open(unit=101,file=displacement_file,status='unknown')
    endif
endif

! Strain
if (iWantStrain) then
    if (verbosity.ge.2) then
        write(stdout,*) 'calc_deformation: calculating strains'
    endif
    if (strain_file.ne.'') then
        open(unit=111,file=strain_file,status='unknown')
    endif
endif

! Rotation
if (iWantRotation) then
    if (verbosity.ge.2) then
        write(stdout,*) 'calc_deformation: calculating rotations'
    endif
    if (rotation_file.ne.'') then
        open(unit=112,file=rotation_file,status='unknown')
    endif
endif

! Stress
if (iWantStress) then
    if (verbosity.ge.2) then
        write(stdout,*) 'calc_deformation: calculating stresses'
    endif
    if (stress_file.ne.'') then
        open(unit=121,file=stress_file,status='unknown')
    endif
    if (estress_file.ne.'') then
        open(unit=122,file=estress_file,status='unknown')
    endif
endif

! Resolved tractions
if (iWantTraction) then
    if (verbosity.ge.2) then
        write(stdout,*) 'calc_deformation: calculating resolved tractions'
    endif
    if (normal_file.ne.'') then
        open(unit=131,file=normal_file,status='unknown')
    endif
    if (shear_file.ne.'') then
        open(unit=132,file=shear_file,status='unknown')
    endif
    if (coulomb_file.ne.'') then
        open(unit=133,file=coulomb_file,status='unknown')
    endif
endif


! Distance to trigger coordinate type warning
warn_dist = 100.0d0
coordTypeWarning = .false.

! Negative depth warning message flag
depthWarning = .false.



! Parallelize the calculation of deformation with OpenMP.
#ifdef USE_OPENMP

if (isParallel) then
    if (nthreads.le.0) then
        nthreads = OMP_GET_MAX_THREADS()
    elseif (nthreads.gt.OMP_GET_MAX_THREADS()) then
        write(stderr,*) '***** User set nthreads to',nthreads,' on the command line'
        write(stderr,*) '***** The maximum number of threads available is',OMP_GET_MAX_THREADS()
        write(stderr,*) '***** Setting nthreads to',OMP_GET_MAX_THREADS()
        nthreads = OMP_GET_MAX_THREADS()
    endif
    if (verbosity.ge.1) then
        write(stdout,*) 'calc_deformation: running with nthreads=',nthreads
    endif
else
    nthreads = 1
    if (nstations*(nfaults+ntensile).gt.10000000) then
        write(stderr,*) '***** There are',OMP_GET_MAX_THREADS(),' threads available'
        write(stderr,*) '***** If your run takes a long time, consider using the "-parallel" option'
    endif
endif

#else

if (isParallel) then
    nthreads = 1
    write(stderr,*) 'o92util was compiled without OpenMP, so running in parallel is not possible'
    write(stderr,*)
endif

#endif



! The following block of code (enclosed between !$OMP PARALLEL and !$OMP PARALLEL END) will be
! executed by multiple threads if o92util is compiled with OpenMP enabled (-fopenmp -frecursive
! options in the GNU compiler).
!
! SHARED variables are shared across threads
! PRIVATE variables are uninitialized and isolated to each thread
! FIRSTPRIVATE variables are initialized with their existing value and isolated to each thread
!$OMP PARALLEL SHARED(nstations,stations,sta_char,staSkipLine,coord_type), &
!$OMP          SHARED(nfaults,faults,ntensile,tensile,fault_type), &
!$OMP          SHARED(targets), &
!$OMP          SHARED(lame,shearmod), &
!$OMP          SHARED(iWantDisp,iWantStrain,iWantRotation,iWantStress,iWantTraction), &
!$OMP          SHARED(debug_mode,iWantProg,depthWarning,coordTypeWarning), &
!$OMP          SHARED(displacement_file, strain_file, rotation_file, stress_file, estress_file), &
!$OMP          SHARED(normal_file, shear_file, coulomb_file, keepAllLines), &
!$OMP          SHARED(disp_output_mode), &
!$OMP&         FIRSTPRIVATE(warn_dist), &
!$OMP&         PRIVATE(iSta,sta_coord), &
!$OMP&         PRIVATE(iFlt,evlo,evla,evdp,str,dip,rak,slip_mag,wid,len,slip,mom), &
!$OMP&         PRIVATE(dist,az,test_dist), &
!$OMP&         PRIVATE(disptmp,stntmp,rottmp), &
!$OMP&         PRIVATE(svec,nvec,trac,tupd,tstr,tnor), &
!$OMP&         PRIVATE(disp,stn,rot,sts,ests,tshr,tshrmx,coul), &
!$OMP&         PRIVATE(ierr), &
!$OMP&         DEFAULT(NONE), &
!$OMP&         NUM_THREADS(nthreads)

! if (OMP_GET_THREAD_NUM().eq.0) then
!     write(stdout,*) 'OMP_NUM_THREADS=',OMP_GET_NUM_THREADS()
! endif

! The following loop is executed in parallel, with chunks of size 1, and ordering the output
! Chunks of size 1 interleaves the iterations, massively speeding up the ordered loop
!$OMP DO SCHEDULE(STATIC,1) ORDERED



! Calculate the requested quantities at each station by summing contribution from each source
do iSta = 1,nstations

    ! Initialize displacement and strain to zero
    disp = 0.0d0
    stn = 0.0d0
    rot = 0.0d0

    ! Initialize lateral station coordinates to zero
    sta_coord(1) = 0.0d0
    sta_coord(2) = 0.0d0

    ! Convert depth from km to m
    sta_coord(3) = stations(iSta,3)*1.0d3
    if (debug_mode.eq.'calc_deformation'.or.debug_mode.eq.'all') then
        write(stdout,'(X,"(DEBUG)",X,A5,I6,3F12.3)') 'ista=',iSta,stations(iSta,:)
    endif


    ! This is the old parallel construct, looping only over the faults to compute the deformation at
    ! each station. This works great for situations with multiple faults, but not so well when there
    ! are multiple stations (the typical case).
    !
    !! Parallelize the calculation of deformation at each station with OpenMP.
    !
    !! The following block of code will be executed by multiple threads if o92util is compiled with
    !! OpenMP enabled (-fopenmp -frecursive options in the GNU compiler).
    !
    !! SHARED variables are shared across threads
    !! PRIVATE variables are uninitialized and isolated to each thread
    !! FIRSTPRIVATE variables are initialized with their existing value and isolated to each thread
    !! REDUCTION defines the operation for shared (deformation) variables
    !!$OMP PARALLEL SHARED(nfaults,ntensile,faults,tensile,stations,lame,shearmod), &
    !!$OMP          SHARED(fault_type,coord_type,iWantDisp,iWantStrain,coordTypeWarning,debug_mode), &
    !!$OMP&         PRIVATE(iFlt,ierr), &
    !!$OMP&         PRIVATE(evlo,evla,evdp,str,dip,rak,slip_mag,wid,len,dist,az,test_dist), &
    !!$OMP&         PRIVATE(slip,mom,disptmp,stntmp), &
    !!$OMP&         FIRSTPRIVATE(iSta,sta_coord,warn_dist), &
    !!$OMP&         REDUCTION(+:disp), REDUCTION(+:stn), &
    !!$OMP&         DEFAULT(NONE)

    !! The following loop is executed in parallel
    !!$OMP DO


    if (.not.staSkipLine(iSta)) then

    ! Calculate deformation from each fault and tensile source at the station
    do iFlt = 1,nfaults+ntensile

        ! Set fault parameters to readable variable names
        if (iFlt.le.nfaults) then
            evlo = faults(iFlt,1)
            evla = faults(iFlt,2)
            evdp = faults(iFlt,3)
            str = faults(iFlt,4)
            dip = faults(iFlt,5)
            rak = faults(iFlt,6)
            slip_mag = faults(iFlt,7)
            wid = faults(iFlt,8)
            len = faults(iFlt,9)
        else
            evlo = tensile(iFlt,1)
            evla = tensile(iFlt,2)
            evdp = tensile(iFlt,3)
            str = tensile(iFlt,4)
            dip = tensile(iFlt,5)
            rak = 0.0d0
            slip_mag = tensile(iFlt,6)
            wid = tensile(iFlt,7)
            len = tensile(iFlt,8)
        endif
        if (debug_mode.eq.'calc_deformation'.or.debug_mode.eq.'all') then
            write(stdout,6821) 'iflt=',iFlt,evlo,evla,evdp,str,dip,rak,slip_mag,wid,len
            6821 format(X,'(DEBUG)',X,A5,I6,3F12.3,0P3F8.1,F10.4,2F12.3)
        endif


        ! Station location relative to fault at origin (ENZ coordinates)
        if (coord_type.eq.'geographic') then
            ! Calculate fault-station distance and azimuth
            call lola2distaz(evlo,evla,stations(iSta,1),stations(iSta,2),dist,az, &
                        'radians','radians',ierr)
            if (ierr.ne.0) then
                call usage('calc_deformation: error calculating fault-station position using '//&
                           'lola2distaz (usage:none)')
            endif
            ! Calculate x(East) and y(North) fault-station coordinates, in m
            sta_coord(1) = dist*radius_earth_m*sin(az)
            sta_coord(2) = dist*radius_earth_m*cos(az)

            ! If user input Cartesian coordinates but forgot the -xy flag, then the distance between
            ! coordinates (measured by pythogorean distance) will typically be larger than 100.
            if (.not.coordTypeWarning) then
                test_dist = sqrt((evlo-stations(iSta,1))**2+(evla-stations(iSta,2))**2)
                if (test_dist.gt.warn_dist) then
                    !$OMP ATOMIC WRITE
                    coordTypeWarning = .true.
                    write(stderr,*) 'calc_deformation: found large fault-station distance'
                    write(stderr,*) 'Did you mean to use the -xy flag for Cartesian coordinates?'
                endif
            endif

        elseif (coord_type.eq.'cartesian') then
            ! Calculate x and y fault-station coordinates, in m
            sta_coord(1) = (stations(iSta,1) - evlo)*1.0d3
            sta_coord(2) = (stations(iSta,2) - evla)*1.0d3

        else
            call usage('calc_deformation: no coordinate type named "'//trim(coord_type)//'" (usage:misc)')
        endif
        if (debug_mode.eq.'calc_deformation'.or.debug_mode.eq.'all') then
            write(stdout,6822) 'sta_coord (pre-rotation): ',sta_coord(:),iSta,iFlt
            6822 format(X,'(DEBUG)',5X,A26,3F14.4,' (ista=',I6,', iflt=',I6,')')
        endif


        ! Rotate coordinate axes so x is parallel to strike and y is horizontal up-dip
        call rotate_vector_angle_axis(sta_coord,str-90.0d0,'z',sta_coord,ierr)
        if (ierr.ne.0) then
            call usage('calc_deformation: error in rotating station coordinate axes to strike '//&
                       '(usage:none)')
        endif
        if (debug_mode.eq.'calc_deformation'.or.debug_mode.eq.'all') then
            write(stdout,6823) 'sta_coord (post-rotation):',sta_coord(:),iSta,iFlt
            6823 format(X,'(DEBUG)',5X,A26,3F14.4,' (ista=',I6,', iflt=',I6,')')
        endif


        ! Fault slip vector in strike-slip, dip-slip, and tensile-slip
        if (iFlt.le.nfaults) then
            slip(1) = slip_mag*cos(rak*d2r)
            slip(2) = slip_mag*sin(rak*d2r)
            slip(3) = 0.0d0
        else
            slip(1) = 0.0d0
            slip(2) = 0.0d0
            slip(3) = slip_mag
        endif
        if (fault_type.eq.'point') then
            mom(1) = slip(1)*wid*len*shearmod
            mom(2) = slip(2)*wid*len*shearmod
            mom(3:4) = 0.0d0
        elseif (fault_type.eq.'rect') then
            ! Everything is set up for rectangular sources already
        else
            call usage('calc_deformation: no fault type named "'//trim(fault_type)//'" (usage:input)')
        endif


        ! Calculate deformation for each fault-station pair and add to total at current station
        ! Calculate displacement
        if (iWantDisp) then
            if (fault_type.eq.'point') then
                call o92_pt_disp(disptmp,sta_coord,evdp,dip,mom,lame,shearmod)
            elseif (fault_type.eq.'rect') then
                call o92_rect_disp(disptmp,sta_coord,evdp,dip,slip,wid,len,lame,shearmod)
            endif
            ! Rotate displacements back to x=E, y=N, z=up
            call rotate_vector_angle_axis(disptmp,90.0d0-str,'z',disptmp,ierr)
            if (ierr.ne.0) then
                call usage('calc_deformation: error in rotating displacement vector to ENV (usage:none)')
            endif

            if (debug_mode.eq.'calc_deformation'.or.debug_mode.eq.'all') then
                write(stdout,6824) 'disptmp:  ',disptmp(:),iSta,iFlt
                6824 format(X,'(DEBUG)',5X,A8,3E14.4,18X,' (ista=',I6,', iflt=',I6,')')
            endif

            ! Add to total
            disp = disp + disptmp
        endif

        ! Calculate strain
        if (iWantStrain) then
            if (fault_type.eq.'point') then
                call o92_pt_strain(stntmp,sta_coord,evdp,dip,mom,lame,shearmod)
            elseif (fault_type.eq.'rect') then
                call o92_rect_strain(stntmp,sta_coord,evdp,dip,slip,wid,len,lame,shearmod)
            endif
            ! Rotate strain tensor back to x=E, y=N, z=up
            call rotate_matrix_angle_axis(stntmp,90.0d0-str,'z',stntmp,ierr)
            if (ierr.ne.0) then
                call usage('calc_deformation: error in rotating strain tensor to ENV (usage:none)')
            endif

            if (debug_mode.eq.'calc_deformation'.or.debug_mode.eq.'all') then
                write(stdout,'(X,"(DEBUG)",5X,A8,1P6E14.6)') 'stntmp: ', stntmp(1,1) ,stntmp(2,2), &
                                                stntmp(3,3), stntmp(1,2), stntmp(1,3), stntmp(2,3)
            endif

            ! Add to total
            stn = stn + stntmp
        endif

        ! Calculate rotation
        if (iWantRotation) then
            if (fault_type.eq.'point') then
                call o92_pt_rotation(rottmp,sta_coord,evdp,dip,mom,lame,shearmod)
            elseif (fault_type.eq.'rect') then
                call o92_rect_rotation(rottmp,sta_coord,evdp,dip,slip,wid,len,lame,shearmod)
            endif
            ! Rotate rotation tensor back to x=E, y=N, z=up
            call rotate_matrix_angle_axis(rottmp,90.0d0-str,'z',rottmp,ierr)
            if (ierr.ne.0) then
                call usage('calc_deformation: error in rotating rotation tensor to ENV (usage:none)')
            endif

            if (debug_mode.eq.'calc_deformation'.or.debug_mode.eq.'all') then
                write(stdout,'(X,"(DEBUG)",5X,A8,1P6E14.6)') 'rottmp: ', rottmp(1,1) ,rottmp(2,2), &
                                                rottmp(3,3), rottmp(1,2), rottmp(1,3), rottmp(2,3)
            endif

            ! Add to total
            rot = rot + rottmp
        endif
    enddo

    endif ! (.not.staSkipLine(iSta))


    !! End the parallel do loop
    !!$OMP END DO

    !! Stop parallelizing - the rest of the calculation for this station is serial
    !!$OMP END PARALLEL


    ! Enforce that the output is ordered the same way as it would be in a serial run.
    ! This slows things down, but is more practical for my use cases.
    !$OMP ORDERED


    if (debug_mode.eq.'calc_deformation'.or.debug_mode.eq.'all') then
        if(iWantDisp) then
            write(stdout,6826) 'disp: ',disp,iSta
            6826 format(X,'(DEBUG)',5X,A8,3E16.6,18X,' (ista=',I6,')')
        endif
        if(iWantStrain) then
            write(stdout,*) stn
        endif
    endif


    ! Displacement: ux, uy, uz  OR  az, uhor, uz
    if (displacement_file.ne.'') then
        if (staSkipLine(iSta)) then
            if (keepAllLines) then
                write(101,'(A)') sta_char(iSta)
            endif
        elseif (disp_output_mode.eq.'enz') then
            write(101,*) stations(iSta,:),disp
        elseif (disp_output_mode.eq.'amz') then
            write(101,*) stations(iSta,:), &
                         atan2(disp(1),disp(2))*r2d, &
                         sqrt(disp(1)*disp(1)+disp(2)*disp(2)), &
                         disp(3)
        else
            write(stderr,*) 'calc_deformation: no displacement output mode named "', &
                            trim(disp_output_mode),'"'
            call usage('this is probably not your fault unless you are Matt (usage:none)')
        endif
    endif

    ! Strain tensor: exx, eyy, ezz, exy, exz, eyz
    if (strain_file.ne.'') then
        if (staSkipLine(iSta)) then
            if (keepAllLines) then
                write(111,'(A)') sta_char(iSta)
            endif
        else
            write(111,*) stations(iSta,:),stn(1,1),stn(2,2),stn(3,3),stn(1,2),stn(1,3),stn(2,3)
        endif
    endif

    ! Rotation tensor: rxy, rxz, ryz
    if (rotation_file.ne.'') then
        if (staSkipLine(iSta)) then
            if (keepAllLines) then
                write(112,'(A)') sta_char(iSta)
            endif
        else
            write(112,*) stations(iSta,:),rot(1,2),rot(1,3),rot(2,3)
        endif
    endif

    ! Stress
    if (iWantStress) then
        call strain2stress(stn,lame,shearmod,sts)

        ! Stress tensor: sxx, syy, szz, sxy, sxz, syz
        if (stress_file.ne.'') then
            if (staSkipLine(iSta)) then
                if (keepAllLines) then
                    write(121,'(A)') sta_char(iSta)
                endif
            else
                write(121,*) stations(iSta,:),sts(1,1),sts(2,2),sts(3,3),sts(1,2),sts(1,3),sts(2,3)
            endif
        endif

        ! Maximum (effective) shear stress (aka, second invariant of deviatoric stress tensor): ests
        if (estress_file.ne.'') then
            call max_shear_stress(sts,ests)
            if (staSkipLine(iSta)) then
                if (keepAllLines) then
                    write(122,'(A)') sta_char(iSta)
                endif
            else
                write(122,*) stations(iSta,:),ests
            endif
        endif
    endif

    if (iWantTraction) then
        ! Calculate components of traction resolved onto plane
        call sdr2sv(targets(iSta,1),targets(iSta,2),targets(iSta,3),svec)     ! slip vector
        call strdip2normal(targets(iSta,1),targets(iSta,2),nvec)              ! normal vector
        call stress2traction(sts,nvec,trac)                                   ! traction on plane
        call traction_components(trac,nvec,tnor,tstr,tupd)                    ! traction components
        tshr = tstr*cos(targets(iSta,3)*d2r) + tupd*sin(targets(iSta,3)*d2r)  ! shear traction (projected on slip vector)

        ! Normal traction: normal (positive=dilation)
        if (normal_file.ne.'') then
            if (staSkipLine(iSta)) then
                if (keepAllLines) then
                    write(131,'(A)') sta_char(iSta)
                endif
            else
                write(131,*) stations(iSta,:),tnor
            endif
        endif

        ! Shear traction: resolved_onto_rake, max_shear_on_plane
        if (shear_file.ne.'') then
            tshrmx = sqrt(tstr*tstr+tupd*tupd)
            if (staSkipLine(iSta)) then
                if (keepAllLines) then
                    write(132,'(A)') sta_char(iSta)
                endif
            else
                write(132,*) stations(iSta,:),tshr,tshrmx
            endif
        endif

        ! Coulomb stress: coulomb
        if (coulomb_file.ne.'') then
            coul = tshr + targets(iSta,4)*tnor
            if (staSkipLine(iSta)) then
                if (keepAllLines) then
                    write(133,'(A)') sta_char(iSta)
                endif
            else
                write(133,*) stations(iSta,:),coul
            endif
        endif
    endif

    if (iWantProg) then
        call progress_indicator(iSta,nstations,'o92util calc_deformation',ierr)
        if (ierr.ne.0) then
            call usage('calc_deformation: error in progress_indicator (usage:none)')
        endif
    endif


    ! End the ordered directive, so threads can be FREE!!!!
    !$OMP END ORDERED


enddo


! End the parallel do loop
!$OMP END DO

! Stop parallelizing - the rest of the subroutine is serial
!$OMP END PARALLEL


! Close files that were opened for writing
inquire(file=displacement_file,number=file_unit,opened=isThisUnitOpen)
if (isThisUnitOpen) then
    close(file_unit)
endif
inquire(file=strain_file,number=file_unit,opened=isThisUnitOpen)
if (isThisUnitOpen) then
    close(file_unit)
endif
inquire(file=stress_file,number=file_unit,opened=isThisUnitOpen)
if (isThisUnitOpen) then
    close(file_unit)
endif
inquire(file=estress_file,number=file_unit,opened=isThisUnitOpen)
if (isThisUnitOpen) then
    close(file_unit)
endif
inquire(file=normal_file,number=file_unit,opened=isThisUnitOpen)
if (isThisUnitOpen) then
    close(file_unit)
endif
inquire(file=shear_file,number=file_unit,opened=isThisUnitOpen)
if (isThisUnitOpen) then
    close(file_unit)
endif
inquire(file=coulomb_file,number=file_unit,opened=isThisUnitOpen)
if (isThisUnitOpen) then
    close(file_unit)
endif


if (verbosity.ge.1) then
    write(stdout,*) 'calc_deformation: finished'
    write(stdout,*)
endif

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine auto_stations()
!----
! Generate transects through fault source for determining automatic grid limits.
!----

use io, only: verbosity, stdout, stderr
use earth, only: radius_earth_km
use geom, only: distaz2lola
use trig, only: r2d

use o92util, only: station_file, &
                   coord_type, &
                   nfaults, &
                   faults, &
                   centroid, &
                   auto_depth, &
                   auto_az, &
                   auto_mode, &
                   n_auto_transect, &
                   displacement_file, &
                   disp_file_save, &
                   iWantDisp, &
                   debug_mode

implicit none

! Local variables
integer :: i, ierr
double precision :: moment, lon, lat, x, y, dx, dy


if (verbosity.ge.1) then
    write(stdout,*) 'auto_stations: starting'
endif


! Initialize I/O flag
ierr = 0


! Check whether station file indicates to use auto station mode
if (station_file.ne.'o92_autosta_86_this_when_finished') then
    return
else
    open(unit=81,file=station_file,status='unknown')
endif


! Grid is centered on on the earthquake centroid; calculate the centroid (approximately)
centroid = 0.0d0
moment = 0.0d0
do i = 1,nfaults
    centroid(1) = centroid(1) + faults(i,1)*faults(i,7)*faults(i,8)*faults(i,9)
    centroid(2) = centroid(2) + faults(i,2)*faults(i,7)*faults(i,8)*faults(i,9)
    centroid(3) = centroid(3) + faults(i,3)*faults(i,7)*faults(i,8)*faults(i,9)
    moment = moment + faults(i,7)*faults(i,8)*faults(i,9)
enddo
centroid = centroid/moment
write(stdout,1001) centroid(1),centroid(2),centroid(3)/1.0d3
if (isnan(centroid(1))) then
    write(stderr,*) 'auto_stations: something went wrong calculating centroid: NaN found'
    call usage('Check input fault file to make sure it is formatted correctly (usage:none)')
endif
1001 format(' auto_stations: centroid(lon,lat,dep(km)): ',2(F10.3),F9.2)



! Automatic grid can be horizontal or vertical
if (auto_mode.eq.'h') then

    ! Automatic horizontal grid

    ! Sample along transects extending 500 km west-east, and south-north of centroid

    ! In situations where focal planes are N-S or E-W, this approach gives bad ranges,
    ! so add NW-SE and SW-NE transects
    n_auto_transect = 4004

    if (coord_type.eq.'geographic') then
        ! West-east
        if (verbosity.ge.2) then
            write(stdout,*) 'auto_stations: generating W-E transect for horizontal grid'
        endif
        dx = -500.0d0
        do while (dx.le.500.0d0)
            call distaz2lola(centroid(1),centroid(2),dx/radius_earth_km,90.0d0,lon,lat, &
                             'radians','degrees',ierr)
            if (ierr.ne.0) then
                call usage('auto_stations: error computing longitude and latitude (usage:none)')
            endif
            write(81,*) lon,lat,auto_depth
            if (debug_mode.eq.'auto_stations'.or.debug_mode.eq.'all') then
                write(stdout,'(X,"(DEBUG)",X,3F14.4)') lon,lat,auto_depth
            endif
            dx = dx + 1.0d0
        enddo

        ! South-north
        if (verbosity.ge.2) then
            write(stdout,*) 'auto_stations: generating S-N transect for horizontal grid'
        endif
        dy = -500.0d0
        do while (dy.le.500.0d0)
            call distaz2lola(centroid(1),centroid(2),dy/radius_earth_km,0.0d0,lon,lat, &
                             'radians','degrees',ierr)
            if (ierr.ne.0) then
                call usage('auto_stations: error computing longitude and latitude (usage:none)')
            endif
            write(81,*) lon,lat,auto_depth
            if (debug_mode.eq.'auto_stations'.or.debug_mode.eq.'all') then
                write(stdout,'(X,"(DEBUG)",X,3F14.4)') lon,lat,auto_depth
            endif
            dy = dy + 1.0d0
        enddo

        ! Northwest-southeast
        if (verbosity.ge.2) then
            write(stdout,*) 'auto_stations: generating NW-SE transect for horizontal grid'
        endif
        dy = -500.0d0
        do while (dy.le.500.0d0)
            call distaz2lola(centroid(1),centroid(2),dy/radius_earth_km,135.0d0,lon,lat, &
                             'radians','degrees',ierr)
            if (ierr.ne.0) then
                call usage('auto_stations: error computing longitude and latitude (usage:none)')
            endif
            write(81,*) lon,lat,auto_depth
            if (debug_mode.eq.'auto_stations'.or.debug_mode.eq.'all') then
                write(stdout,'(X,"(DEBUG)",X,3F14.4)') lon,lat,auto_depth
            endif
            dy = dy + 1.0d0
        enddo

        ! Southwest-northeast
        if (verbosity.ge.2) then
            write(stdout,*) 'auto_stations: generating SW-NE transect for horizontal grid'
        endif
        dy = -500.0d0
        do while (dy.le.500.0d0)
            call distaz2lola(centroid(1),centroid(2),dy/radius_earth_km,45.0d0,lon,lat, &
                             'radians','degrees',ierr)
            if (ierr.ne.0) then
                call usage('auto_stations: error computing longitude and latitude (usage:none)')
            endif
            write(81,*) lon,lat,auto_depth
            if (debug_mode.eq.'auto_stations'.or.debug_mode.eq.'all') then
                write(stdout,'(X,"(DEBUG)",X,3F14.4)') lon,lat,auto_depth
            endif
            dy = dy + 1.0d0
        enddo

    elseif (coord_type.eq.'cartesian') then
        ! Left-right (west-east)
        if (verbosity.ge.2) then
            write(stdout,*) 'auto_stations: generating L-R transect for horizontal grid'
        endif
        dx = -500.0d0
        do while (dx.le.500.0d0)
            x = centroid(1) + dx
            write(81,*) x,centroid(2),auto_depth
            if (debug_mode.eq.'auto_stations'.or.debug_mode.eq.'all') then
                write(stdout,'(X,"(DEBUG)",X,3F14.4)') x,centroid(2),auto_depth
            endif
            dx = dx + 1.0d0
        enddo

        ! Bottom-top (south-north)
        if (verbosity.ge.2) then
            write(stdout,*) 'auto_stations: generating B-T transect for horizontal grid'
        endif
        dy = -500.0d0
        do while (dy.le.500.0d0)
            y = centroid(2) + dy
            write(81,*) centroid(1),y,auto_depth
            if (debug_mode.eq.'auto_stations'.or.debug_mode.eq.'all') then
                write(stdout,'(X,"(DEBUG)",X,3F14.4)') centroid(1),y,auto_depth
            endif
            dy = dy + 1.0d0
        enddo

        ! Top/left-bottom/right (NW-SE)
        if (verbosity.ge.2) then
            write(stdout,*) 'auto_stations: generating TL-BR transect for horizontal grid'
        endif
        dy = -500.0d0
        do while (dy.le.500.0d0)
            x = centroid(1) + dy*0.5d0*sqrt(2.0d0)
            y = centroid(2) - dy*0.5d0*sqrt(2.0d0)
            write(81,*) x,y,auto_depth
            if (debug_mode.eq.'auto_stations'.or.debug_mode.eq.'all') then
                write(stdout,'(X,"(DEBUG)",X,3F14.4)') x,y,auto_depth
            endif
            dy = dy + 1.0d0
        enddo

        ! Bottom/left-top/right (SW-NE)
        if (verbosity.ge.2) then
            write(stdout,*) 'auto_stations: generating BL-TR transect for horizontal grid'
        endif
        dy = -500.0d0
        do while (dy.le.500.0d0)
            x = centroid(1) + dy*0.5d0*sqrt(2.0d0)
            y = centroid(2) + dy*0.5d0*sqrt(2.0d0)
            write(81,*) x,y,auto_depth
            if (debug_mode.eq.'auto_stations'.or.debug_mode.eq.'all') then
                write(stdout,'(X,"(DEBUG)",X,3F14.4)') x,y,auto_depth
            endif
            dy = dy + 1.0d0
        enddo

    endif


elseif (auto_mode.eq.'v') then

    ! Automatic vertical grid

    ! Sample along transects extending 500 km negative-positive azimuth, up-down of centroid
    n_auto_transect = 2002

    if (coord_type.eq.'geographic') then
        ! Along azimuth
        if (verbosity.ge.2) then
            write(stdout,*) 'auto_stations: generating azimuthal transect for vertical grid'
        endif
        dx = -500.0d0
        do while (dx.le.500.0d0)
            call distaz2lola(centroid(1),centroid(2),dx/radius_earth_km,auto_az,lon,lat, &
                             'radians','degrees',ierr)
            if (ierr.ne.0) then
                call usage('auto_stations: error computing longitude and latitude (usage:none)')
            endif
            write(81,*) lon,lat,centroid(3)/1.0d3
            if (debug_mode.eq.'auto_stations'.or.debug_mode.eq.'all') then
                write(stdout,'(X,"(DEBUG)",X,3F14.4)') lon,lat,centroid(3)/1.0d3
            endif
            dx = dx + 1.0d0
        enddo
    elseif (coord_type.eq.'cartesian') then
        ! Along azimuth
        if (verbosity.ge.2) then
            write(stdout,*) 'auto_stations: generating azimuthal transect for vertical grid'
        endif
        dx = -500.0d0
        do while (dx.le.500.0d0)
            x = centroid(1) + dx*sin(auto_az*r2d)
            y = centroid(2) + dx*cos(auto_az*r2d)
            write(81,*) x,y,centroid(3)/1.0d3
            if (debug_mode.eq.'auto_stations'.or.debug_mode.eq.'all') then
                write(stdout,'(X,"(DEBUG)",X,3F14.4)') x,y,centroid(3)/1.0d3
            endif
            dx = dx + 1.0d0
        enddo
    endif

    ! Vertical
    if (verbosity.ge.2) then
        write(stdout,*) 'auto_stations: generating vertical transect for vertical grid'
    endif
    dy = -500.0d0
    do while (dy.le.500.0d0)
        y = centroid(3)/1.0d3 + dy
        dy = dy + 1.0d0
        write(81,*) centroid(1),centroid(2),y
        if (debug_mode.eq.'auto_stations'.or.debug_mode.eq.'all') then
            write(stdout,'(X,"(DEBUG)",X,3F14.4)') centroid(1),centroid(2),y
        endif
    enddo

else
    call usage('auto_stations: no auto_mode named '//trim(auto_mode)//' (usage:station)')
endif


! Finished writing transects to temporary station file
close(81)


! Update the displacement file name to a temporary value, because this will be used for setting
! the automatic grid limits
disp_file_save = displacement_file
displacement_file = 'o92_autosta_disp_86_this_when_finished'
iWantDisp = .true.


if (verbosity.ge.1) then
    write(stdout,*) 'auto_stations: finished'
    write(stdout,*)
endif

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine update_auto_stations()
!----
! From the displacements that were calculated along the transects, determine the limits of the
! automatic grid, given a displacement threshold. Use these limits to define a complete grid.
!----

use io, only: verbosity, stdout
use earth, only: radius_earth_km
use geom, only: distaz2lola
use trig, only: d2r

use o92util, only: coord_type, &
                   displacement_file, &
                   disp_file_save, &
                   station_file, &
                   nstations, &
                   stations, &
                   centroid, &
                   auto_n, &
                   auto_depth, &
                   auto_az, &
                   auto_disp_threshold, &
                   auto_mode, &
                   n_auto_transect, &
                   debug_mode

implicit none

! Local variables
integer :: i, j, ierr
double precision :: disp(n_auto_transect,6), disp_mag, xmin, xmax, ymin, ymax
double precision :: xmin3(3), xmax3(3), ymin3(3), ymax3(3)
double precision :: d, lon, lat
character(len=512) :: proj_sta_file


if (verbosity.ge.1) then
    write(stdout,*) 'update_auto_stations: starting'
endif


! Station file name indicates whether using auto grid mode
if (station_file.ne.'o92_autosta_86_this_when_finished') then
    return
else
    open(unit=81,file=station_file,status='old')
    open(unit=82,file=displacement_file,status='old')
endif


! Read the displacements that were calculated from the station transects defined in auto_stations()
if (debug_mode.eq.'update_auto_stations'.or.debug_mode.eq.'all') then
    write(stdout,'(X,"(DEBUG)",X,6A14)') 'x','y','z','ux','uy','uz'
endif
do i = 1,n_auto_transect
    read(82,*) (disp(i,j),j=1,6)
    if (debug_mode.eq.'update_auto_stations'.or.debug_mode.eq.'all') then
        write(stdout,'(X,"(DEBUG)",X,1P6E14.6)') disp(i,:)
    endif
enddo


! Working backwards from ends of transects, determine when displacements exceed a threshold value
! and treat that coordinate as the new station grid boundary

! If displacements are larger than user-defined threshold, use the farthest points along the
! transect. For maintaining the elastic half-space approximation (instead of a spherical Earth),
! this is currently hard-coded to be 500 km in subroutine auto_stations().
if (auto_mode.eq.'h') then
    xmin = disp(1,1)
    xmax = disp(1001,1)
elseif (auto_mode.eq.'v') then
    xmin = -500.0d0
    xmax = 500.0d0
else
    xmin = -1d10
    xmax = 1d10
    call usage('update_auto_stations: no auto_mode named '//trim(auto_mode)//' (usage:station)')
endif
ymin = disp(1002,2)
ymax = disp(2002,2)


! Find limits for horizontal grid
if (auto_mode.eq.'h') then

    ! Initialize values for min/max values along each transect
    xmin3 = 0.0d0
    xmax3 = 0.0d0
    ymin3 = 0.0d0
    ymax3 = 0.0d0


    if (verbosity.ge.3) then
        write(stdout,*) 'update_auto_stations: looking for xmin:'
    endif
    ! Find xmin along W-E line
    do i = 1,501
        disp_mag = sqrt(disp(i,4)**2+disp(i,5)**2+disp(i,6)**2)
        if (disp_mag.ge.auto_disp_threshold) then
            if (verbosity.ge.3) then
                write(stdout,2042) 'on W-E line:  ',disp(i,1)
            endif
            xmin3(1) = disp(i,1)
            exit
        endif
    enddo

    ! Find xmin along NW-SE line
    do i = 2003,2503
        disp_mag = sqrt(disp(i,4)**2+disp(i,5)**2+disp(i,6)**2)
        if (disp_mag.ge.auto_disp_threshold) then
            if (verbosity.ge.3) then
                write(stdout,2042) 'on NW-SE line:',disp(i,1)
            endif
            xmin3(2) = disp(i,1)
            exit
        endif
    enddo

    ! Find xmin along SW-NE line
    do i = 3004,3504
        disp_mag = sqrt(disp(i,4)**2+disp(i,5)**2+disp(i,6)**2)
        if (disp_mag.ge.auto_disp_threshold) then
            if (verbosity.ge.3) then
                write(stdout,2042) 'on SW-NE line:',disp(i,1)
            endif
            xmin3(3) = disp(i,1)
            exit
        endif
    enddo


    if (verbosity.ge.3) then
        write(stdout,*) 'update_auto_stations: looking for xmax:'
    endif
    ! Find xmax along W-E line
    do i = 1001,501,-1
        disp_mag = sqrt(disp(i,4)**2+disp(i,5)**2+disp(i,6)**2)
        if (disp_mag.ge.auto_disp_threshold) then
            if (verbosity.ge.3) then
                write(stdout,2042) 'on W-E line:  ',disp(i,1)
            endif
            xmax3(1) = disp(i,1)
            exit
        endif
    enddo

    ! Find xmax along NW-SE line
    do i = 3003,2503,-1
        disp_mag = sqrt(disp(i,4)**2+disp(i,5)**2+disp(i,6)**2)
        if (disp_mag.ge.auto_disp_threshold) then
            if (verbosity.ge.3) then
                write(stdout,2042) 'on NW-SE line:',disp(i,1)
            endif
            xmax3(2) = disp(i,1)
            exit
        endif
    enddo

    ! Find xmax along SW-NE line
    do i = 4004,3504,-1
        disp_mag = sqrt(disp(i,4)**2+disp(i,5)**2+disp(i,6)**2)
        if (disp_mag.ge.auto_disp_threshold) then
            if (verbosity.ge.3) then
                write(stdout,2042) 'on SW-NE line:',disp(i,1)
            endif
            xmax3(3) = disp(i,1)
            exit
        endif
    enddo


    if (verbosity.ge.3) then
        write(stdout,*) 'update_auto_stations: looking for ymin:'
    endif
    ! Find ymin along S-N line
    do i = 1002,1502
        disp_mag = sqrt(disp(i,4)**2+disp(i,5)**2+disp(i,6)**2)
        if (disp_mag.ge.auto_disp_threshold) then
            if (verbosity.ge.3) then
                write(stdout,2042) 'on S-N line:  ',disp(i,2)
            endif
            ymin3(1) = disp(i,2)
            exit
        endif
    enddo

    ! Find ymin along NW-SE line
    do i = 3003,2503,-1
        disp_mag = sqrt(disp(i,4)**2+disp(i,5)**2+disp(i,6)**2)
        if (disp_mag.ge.auto_disp_threshold) then
            if (verbosity.ge.3) then
                write(stdout,2042) 'on NW-SE line:',disp(i,2)
            endif
            ymin3(2) = disp(i,2)
            exit
        endif
    enddo

    ! Find ymin along SW-NE line
    do i = 3004,3504
        disp_mag = sqrt(disp(i,4)**2+disp(i,5)**2+disp(i,6)**2)
        if (disp_mag.ge.auto_disp_threshold) then
            if (verbosity.ge.3) then
                write(stdout,2042) 'on SW-NE line:',disp(i,2)
            endif
            ymin3(3) = disp(i,2)
            exit
        endif
    enddo


    if (verbosity.ge.3) then
        write(stdout,*) 'update_auto_stations: looking for ymax:'
    endif
    ! Find ymax along S-N line
    do i = 2002,1502,-1
        disp_mag = sqrt(disp(i,4)**2+disp(i,5)**2+disp(i,6)**2)
        if (disp_mag.ge.auto_disp_threshold) then
            if (verbosity.ge.3) then
                write(stdout,2042) 'on S-N line:  ',disp(i,2)
            endif
            ymax3(1) = disp(i,2)
            exit
        endif
    enddo

    ! Find ymax along NW-SE line
    do i = 2003,2503
        disp_mag = sqrt(disp(i,4)**2+disp(i,5)**2+disp(i,6)**2)
        if (disp_mag.ge.auto_disp_threshold) then
            if (verbosity.ge.3) then
                write(stdout,2042) 'on NW-SE line:',disp(i,2)
            endif
            ymax3(2) = disp(i,2)
            exit
        endif
    enddo

    ! Find ymax along SW-NE line
    do i = 4004,3504,-1
        disp_mag = sqrt(disp(i,4)**2+disp(i,5)**2+disp(i,6)**2)
        if (disp_mag.ge.auto_disp_threshold) then
            if (verbosity.ge.3) then
                write(stdout,2042) 'on SW-NE line:',disp(i,2)
            endif
            ymax3(3) = disp(i,2)
            exit
        endif
    enddo

    2042 format(5X,A,F10.3)

    xmin = minval(xmin3)
    xmax = maxval(xmax3)
    ymin = minval(ymin3)
    ymax = minval(ymax3)
    write(stdout,*) 'update_auto_stations: using updated map limits for horizontal grid:'
    write(stdout,'(5X,"xmin=",X,F10.3)') xmin
    write(stdout,'(5X,"xmax=",X,F10.3)') xmax
    write(stdout,'(5X,"ymin=",X,F10.3)') ymin
    write(stdout,'(5X,"ymax=",X,F10.3)') ymax
    write(stdout,*) 'update_auto_stations: saving updated values in o92util_auto_lims.dat'
    open(unit=167,file='o92util_auto_lims.dat',status='unknown')
    write(167,*) 'xmin',xmin
    write(167,*) 'xmax',xmax
    write(167,*) 'ymin',ymin
    write(167,*) 'ymax',ymax
    close(167)



! Find limits for vertical grid
elseif (auto_mode.eq.'v') then

    do i = 1,501
        disp_mag = sqrt(disp(i,4)**2+disp(i,5)**2+disp(i,6)**2)

        ! Check whether displacement magnitude has exceed threshold
        if (disp_mag.ge.auto_disp_threshold) then
            xmin = -501.0d0+dble(i) ! distance along transect
            if (verbosity.ge.3) then
                write(stdout,*) 'update_auto_stations: transect_min=',xmin
            endif
            exit
        endif
    enddo

    ! Find xmax
    do i = 1001,501,-1
        disp_mag = sqrt(disp(i,4)**2+disp(i,5)**2+disp(i,6)**2)

        ! Check whether displacement magnitude has exceed threshold
        if (disp_mag.ge.auto_disp_threshold) then
            xmax = -501.0d0+dble(i) ! distance along transect
            if (verbosity.ge.3) then
                write(stdout,*) 'update_auto_stations: transect_max=',xmax
            endif
            exit
        endif
    enddo

    ! Find ymin
    do i = 1002,1502
        disp_mag = sqrt(disp(i,4)**2+disp(i,5)**2+disp(i,6)**2)

        ! Check whether displacement magnitude has exceed threshold
        if (disp_mag.ge.auto_disp_threshold) then
            ymin = disp(i,3)
            if (verbosity.ge.3) then
                write(stdout,*) 'update_auto_stations: zmin=',ymin
            endif
            exit
        endif
    enddo

    ! Find ymax
    do i = 2002,1502,-1
        disp_mag = sqrt(disp(i,4)**2+disp(i,5)**2+disp(i,6)**2)

        ! Check whether displacement magnitude has exceed threshold
        if (disp_mag.ge.auto_disp_threshold) then
            ymax = disp(i,3)
            if (verbosity.ge.3) then
                write(stdout,*) 'update_auto_stations: zmax=',ymax
            endif
            exit
        endif
    enddo


else

    call usage('update_auto_stations: no auto_mode named '//trim(auto_mode)//' (usage:station)')

endif ! Finished finding auto grid limits




! Finished with temporary station and displacement files - delete these
close(81,status='delete')
close(82,status='delete')



! Create new station grid from these limits
if (verbosity.ge.2) then
    write(stdout,*) 'update_auto_stations: creating new station file'
endif


! Reallocate memory to station array
deallocate(stations)
nstations = auto_n*auto_n
allocate(stations(nstations,3))


! Load the new station array with the automatic grid
if (auto_mode.eq.'h') then
    do i = 1,auto_n
        do j = 1,auto_n
            stations((i-1)*auto_n+j,1) = xmin + dble(i-1)*(xmax-xmin)/dble(auto_n-1)
            stations((i-1)*auto_n+j,2) = ymin + dble(j-1)*(ymax-ymin)/dble(auto_n-1)
            stations((i-1)*auto_n+j,3) = auto_depth
        enddo
    enddo
elseif (auto_mode.eq.'v') then
    ! Also write projected coordinates to a file
    proj_sta_file = 'station.proj'
    write(stdout,*) 'update_auto_stations: writing projected station coordinates to ', &
                    trim(proj_sta_file)
    open(unit=372,file=proj_sta_file,status='unknown')

    ! Set depth to be at surface if ymin is close to that value
    if (ymin.lt.2d3) then
        ymin = 0.0d0
    endif

    do i = 1,auto_n
        d = xmin + dble(i-1)*(xmax-xmin)/dble(auto_n-1)
        if (coord_type.eq.'geographic') then
            call distaz2lola(centroid(1),centroid(2),d/radius_earth_km,auto_az,lon,lat, &
                             'radians','degrees',ierr)
        elseif (coord_type.eq.'cartesian') then
            lon = centroid(1) + d*sin(auto_az*d2r)
            lat = centroid(2) + d*cos(auto_az*d2r)
        else
            call usage('update_auto_stations: no coord_type named '//trim(coord_type)//' (usage:misc)')
        endif
        do j = 1,auto_n
            stations((i-1)*auto_n+j,1) = lon
            stations((i-1)*auto_n+j,2) = lat
            stations((i-1)*auto_n+j,3) = ymin + dble(j-1)*(ymax-ymin)/dble(auto_n-1)
            write(372,'(3F12.4)') d,0.0d0,stations((i-1)*auto_n+j,3)
        enddo
    enddo

    close(372)
else
    call usage('update_auto_stations: no auto_mode named '//trim(auto_mode)//' (usage:station)')
endif


! Reset displacement file name
displacement_file = disp_file_save


if (verbosity.ge.1) then
    write(stdout,*) 'update_auto_stations: finished'
    write(stdout,*)
endif

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine project_deformation()
!----
! Project vectors or tensors onto vertical grid
!----

use io, only: stdout, line_count
use algebra, only: rotate_vector_angle_axis, rotate_matrix_angle_axis

use o92util, only: displacement_file, &
                   strain_file, &
                   auto_az

implicit none

! Local variables
integer :: i, j, k, ierr, ndisp, nstrain
double precision :: sta(3), disp(3), disp_rot(3), stn(3,3), stn_rot(3,3)
character(len=512) :: auto_file


! Project displacement vector
if (displacement_file.ne.'') then
    ndisp = line_count(displacement_file)

    auto_file = trim(displacement_file)//'.proj'
    write(stdout,*) 'project_deformation: writing projected displacements to file '//auto_file

    open(unit=61,file=displacement_file,status='old')
    open(unit=62,file=auto_file,status='unknown')

    do i = 1,ndisp
        read(61,*) (sta(j),j=1,3),(disp(k),k=1,3)
        call rotate_vector_angle_axis(disp,90.0d0-auto_az,'z',disp_rot,ierr)
        write(62,'(1P3E14.6)') disp_rot(1:3)
    enddo

    close(61)
    close(62)
endif


! Project strain tensor
if (strain_file.ne.'') then
    nstrain = line_count(strain_file)

    auto_file = trim(strain_file)//'.proj'
    write(stdout,*) 'project_deformation: writing projected strains to file '//auto_file

    open(unit=63,file=strain_file,status='old')
    open(unit=64,file=auto_file,status='unknown')

    do i = 1,nstrain
        read(63,*) (sta(j),j=1,3),stn(1,1),stn(2,2),stn(3,3),stn(1,2),stn(1,3),stn(2,3)
        stn(2,1) = stn(1,2)
        stn(3,1) = stn(1,3)
        stn(3,2) = stn(2,3)
        call rotate_matrix_angle_axis(stn,90.0d0-auto_az,'z',stn_rot,ierr)
        write(64,'(1P6E14.6)') stn_rot(1,1),stn_rot(2,2),stn_rot(3,3), &
                               stn_rot(1,2),stn_rot(1,3),stn_rot(2,3)
    enddo

    close(63)
    close(64)
endif

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine gcmdln()
!----
! Parse o92util command line arguments and set default values
!----

use io, only: stdout, stderr, verbosity, isNumeric
use o92util, only: ffm_file, &
                   fsp_file, &
                   mag_file, &
                   flt_file, &
                   tns_file, &
                   gmt_fault_file, &
                   isFaultFileDefined, &
                   isTensileFileDefined, &
                   fault_type, &
                   empirical_relation, &
                   slip_threshold, &
                   station_file, &
                   isStationFileDefined, &
                   auto_depth, &
                   auto_az, &
                   auto_n, &
                   auto_disp_threshold, &
                   auto_mode, &
                   target_file, &
                   isTargetFileDefined, &
                   halfspace_file, &
                   poisson, &
                   lame, &
                   shearmod, &
                   displacement_file, &
                   disp_output_mode, &
                   strain_file, &
                   rotation_file, &
                   stress_file, &
                   estress_file, &
                   normal_file, &
                   shear_file, &
                   coulomb_file, &
                   isOutputFileDefined, &
                   keepAllLines, &
                   iWantDisp, &
                   iWantStrain, &
                   iWantRotation, &
                   iWantStress, &
                   iWantTraction, &
                   isParallel, &
                   nthreads, &
                   coord_type, &
                   iWantProg, &
                   debug_mode

implicit none

! Local variables
character(len=512) tag, arg
integer :: i, ios, narg

! Initialize control variables
ios = 0

ffm_file = ''
fsp_file = ''
flt_file = ''
mag_file = ''
tns_file = ''
gmt_fault_file = ''
isFaultFileDefined = .false.
isTensileFileDefined = .false.
fault_type = 'rect'
empirical_relation = 'WC'
slip_threshold = 0.0d0

station_file = ''
isStationFileDefined = .false.
auto_mode = ''
auto_depth = 0.0d0
auto_az = 0.0d0
auto_n = 10
auto_disp_threshold = 0.001d0
target_file = ''
isTargetFileDefined = .false.

halfspace_file = ''
poisson = 0.25d0
lame = 40.0d9
shearmod = 40.0d9

displacement_file = ''
disp_output_mode = 'enz'
strain_file = ''
rotation_file = ''
stress_file = ''
estress_file = ''
normal_file = ''
shear_file = ''
coulomb_file = ''
isOutputFileDefined = .false.
keepAllLines = .false.
iWantDisp = .false.
iWantStrain = .false.
iWantStress = .false.
iWantTraction = .false.

isParallel = .false.
nthreads = 0

coord_type = 'geographic'

verbosity = 0
iWantProg = .false.
debug_mode = ''


narg = command_argument_count()
if (narg.eq.0) then
    call usage('')
endif

i = 1
do while (i.le.narg)

    call get_command_argument(i,tag)

    ! Input fault options
    if (trim(tag).eq.'-ffm') then
        i = i + 1
        call get_command_argument(i,ffm_file,status=ios)
        isFaultFileDefined = .true.

    elseif (trim(tag).eq.'-fsp') then
        i = i + 1
        call get_command_argument(i,fsp_file,status=ios)
        isFaultFileDefined = .true.

    elseif (trim(tag).eq.'-mag') then
        i = i + 1
        call get_command_argument(i,mag_file,status=ios)
        isFaultFileDefined = .true.

    elseif (trim(tag).eq.'-flt') then
        i = i + 1
        call get_command_argument(i,flt_file,status=ios)
        isFaultFileDefined = .true.

    elseif (tag.eq.'-tns') then
        i = i + 1
        call get_command_argument(i,tns_file,status=ios)
        isTensileFileDefined = .true.

    elseif (trim(tag).eq.'-fn') then
        fault_type = 'rect'
    elseif (trim(tag).eq.'-pt') then
        fault_type = 'point'

    elseif (trim(tag).eq.'-empirical'.or.trim(tag).eq.'-emp') then
          i = i + 1
          call get_command_argument(i,empirical_relation,status=ios)

    elseif (trim(tag).eq.'-thr') then
        i = i + 1
        call get_command_argument(i,arg,status=ios)
        read(arg,*,err=9001,iostat=ios) slip_threshold


    ! Input receiver options
    elseif (trim(tag).eq.'-sta') then
        i = i + 1
        call get_command_argument(i,station_file,status=ios)
        isStationFileDefined = .true.

    elseif (trim(tag).eq.'-auto') then
        station_file = 'o92_autosta_86_this_when_finished'
        isStationFileDefined = .true.
        i = i + 1
        call get_command_argument(i,auto_mode,status=ios)
        if (auto_mode.ne.'h'.and.auto_mode.ne.'v') then
            write(stderr,*) 'o92util: WARNING: you may be using deprecated "-auto" flag syntax'
            write(stderr,*) 'The auto_mode should be specified with "h" or "v"'
            write(stderr,*) 'Using "h"(orizontal) mode'
            ! call usage('o92util: auto station mode requires "-auto h" or "-auto v" (usage:station)')
            auto_mode = 'h'
        else
            i = i + 1
        endif
        call get_command_argument(i,arg,status=ios)
        if (auto_mode.eq.'h') then
            read(arg,*,err=9001,iostat=ios) auto_depth
        elseif (auto_mode.eq.'v') then
            read(arg,*,err=9001,iostat=ios) auto_az
        endif
        i = i + 1
        call get_command_argument(i,arg)
        read(arg,*) auto_n

    elseif (tag.eq.'-auto:thr') then
        i = i + 1
        call get_command_argument(i,arg,status=ios)
        read(arg,*,err=9001,iostat=ios) auto_disp_threshold


    elseif (trim(tag).eq.'-trg') then
        i = i + 1
        call get_command_argument(i,target_file,status=ios)
        isTargetFileDefined = .true.


    ! Input half-space options
    elseif (trim(tag).eq.'-haf') then
        i = i + 1
        call get_command_argument(i,halfspace_file,status=ios)


    ! Output deformation options
    elseif (trim(tag).eq.'-disp') then
        i = i + 1
        call get_command_argument(i,displacement_file,status=ios)
        isOutputFileDefined = .true.
        iWantDisp = .true.

    elseif (trim(tag).eq.'-strain') then
        i = i + 1
        call get_command_argument(i,strain_file,status=ios)
        isOutputFileDefined = .true.
        iWantStrain = .true.

    elseif (trim(tag).eq.'-rotation') then
        i = i + 1
        call get_command_argument(i,rotation_file,status=ios)
        isOutputFileDefined = .true.
        iWantRotation = .true.

    elseif (trim(tag).eq.'-stress') then
        i = i + 1
        call get_command_argument(i,stress_file,status=ios)
        isOutputFileDefined = .true.
        iWantStrain = .true.
        iWantStress = .true.

    elseif (trim(tag).eq.'-estress') then
        i = i + 1
        call get_command_argument(i,estress_file,status=ios)
        isOutputFileDefined = .true.
        iWantStrain = .true.
        iWantStress = .true.

    elseif (trim(tag).eq.'-normal') then
        i = i + 1
        call get_command_argument(i,normal_file,status=ios)
        isOutputFileDefined = .true.
        iWantStrain = .true.
        iWantStress = .true.
        iWantTraction = .true.

    elseif (trim(tag).eq.'-shear'.or.trim(tag).eq.'-shearmax') then
        i = i + 1
        call get_command_argument(i,shear_file,status=ios)
        isOutputFileDefined = .true.
        iWantStrain = .true.
        iWantStress = .true.
        iWantTraction = .true.

    elseif (trim(tag).eq.'-coul') then
        i = i + 1
        call get_command_argument(i,coulomb_file,status=ios)
        isOutputFileDefined = .true.
        iWantStrain = .true.
        iWantStress = .true.
        iWantTraction = .true.

    elseif (trim(tag).eq.'--keep-all-lines') then
        keepAllLines = .true.

    ! Miscellaneous options
    elseif (trim(tag).eq.'-parallel') then
        isParallel = .true.
        ! Check for optional argument specifying number of threads
        i = i + 1
        if (i.gt.narg) then                       ! "-parallel" was the last command line argument
            cycle
        endif
        call get_command_argument(i,arg)
        if (isNumeric(arg)) then                  ! the argument is numeric, read nthreads
            read(arg,*,iostat=ios) nthreads
            if (ios.ne.0) then                    ! trouble reading numeric argument, not an integer?
                write(stderr,*) 'o92util: could not read ',trim(arg),' as an integer NTHREADS'
                call usage('(usage:misc)')
            endif
        else                                      ! the argument is not numeric, reset i
            i = i - 1
        endif

    elseif (trim(tag).eq.'-parallel-check') then  ! Check whether compiled with OpenMP
#ifdef USE_OPENMP
        write(*,*) 'o92util: compiled with OpenMP'
#else
        write(*,*) 'o92util: compiled without OpenMP'
#endif
        stop


    elseif (trim(tag).eq.'-cartesian'.or.trim(tag).eq.'-xy') then
        coord_type = 'cartesian'
    elseif (trim(tag).eq.'-geographic'.or.trim(tag).eq.'-geo') then
        coord_type = 'geographic'

    elseif (tag.eq.'-az') then
        disp_output_mode = 'amz'

    elseif (tag.eq.'-gmt') then
        i = i + 1
        call get_command_argument(i,gmt_fault_file)

    elseif (trim(tag).eq.'-v'.or.trim(tag).eq.'-verbose'.or.trim(tag).eq.'-verbosity') then
        i = i + 1
        call get_command_argument(i,arg)
        read(arg,*,err=9001,iostat=ios) verbosity

    elseif (tag.eq.'-debug') then
        i = i + 1
        call get_command_argument(i,arg,status=ios)
        if (ios.ne.0) then
            debug_mode = 'all'
            ios = 0
        elseif (index(arg,'-').eq.1) then
            debug_mode = 'all'
            i = i - 1
        else
            read(arg,*,iostat=ios) debug_mode
        endif
        if (verbosity.eq.0) then
            verbosity = 1
        endif

    elseif (trim(tag).eq.'-prog') then
        iWantProg = .true.

    else
        call usage('o92util: no option '//trim(tag))
    endif

    9001 if (ios.ne.0) then
        call usage('o92util: error parsing "'//trim(tag)//'" flag arguments')
    endif

    i = i + 1
enddo

if (debug_mode.eq.'gcmdln'.or.debug_mode.eq.'all') then
    write(stdout,*) 'gcmdln: finished parsing the command line'
    write(stdout,*) 'ffm_file:                 ',trim(ffm_file)
    write(stdout,*) 'fsp_file:                 ',trim(fsp_file)
    write(stdout,*) 'flt_file:                 ',trim(flt_file)
    write(stdout,*) 'mag_file:                 ',trim(mag_file)
    write(stdout,*) 'tns_file:                 ',trim(tns_file)
    write(stdout,*) 'isFaultFileDefined:       ',isFaultFileDefined
    write(stdout,*) 'fault_type:               ',trim(fault_type)
    write(stdout,*) 'empirical_relation:       ',trim(empirical_relation)
    write(stdout,*) 'slip_threshold:           ',slip_threshold
    write(stdout,*) 'station_file:             ',trim(station_file)
    write(stdout,*) 'isStationFileDefined:     ',isStationFileDefined
    write(stdout,*) 'auto_depth:               ',auto_depth
    write(stdout,*) 'auto_az:                  ',auto_az
    write(stdout,*) 'auto_n:                   ',auto_n
    write(stdout,*) 'auto_disp_threshold:      ',auto_disp_threshold
    write(stdout,*) 'auto_mode:                ',trim(auto_mode)
    write(stdout,*) 'target_file:              ',trim(target_file)
    write(stdout,*) 'isTargetFileDefined:      ',isTargetFileDefined
    write(stdout,*) 'halfspace_file:           ',trim(halfspace_file)
    write(stdout,*) 'poisson:                  ',poisson
    write(stdout,*) 'lame:                     ',lame
    write(stdout,*) 'shearmod:                 ',shearmod
    write(stdout,*) 'displacement_file:        ',trim(displacement_file)
    write(stdout,*) 'disp_output_mode:         ',trim(disp_output_mode)
    write(stdout,*) 'strain_file:              ',trim(strain_file)
    write(stdout,*) 'rotation_file:            ',trim(rotation_file)
    write(stdout,*) 'stress_file:              ',trim(stress_file)
    write(stdout,*) 'estress_file:             ',trim(estress_file)
    write(stdout,*) 'normal_file:              ',trim(normal_file)
    write(stdout,*) 'shear_file:               ',trim(shear_file)
    write(stdout,*) 'coulomb_file:             ',trim(coulomb_file)
    write(stdout,*) 'keepAllLines:             ',keepAllLines
    write(stdout,*) 'isOutputFileDefined:      ',isOutputFileDefined
    write(stdout,*) 'iWantDisp:                ',iWantDisp
    write(stdout,*) 'iWantStrain:              ',iWantStrain
    write(stdout,*) 'iWantStress:              ',iWantStress
    write(stdout,*) 'iWantTraction:            ',iWantTraction
    write(stdout,*) 'isParallel:               ',isParallel
    write(stdout,*) 'nthreads:                 ',nthreads
    write(stdout,*) 'coord_type:               ',trim(coord_type)
    write(stdout,*) 'prog:                     ',iWantProg
    write(stdout,*) 'debug_mode:               ',trim(debug_mode)
    write(stdout,*)
endif

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine usage(str)

use io, only: stderr

implicit none

! Arguments
character(len=*) :: str
character(len=8) :: info
integer :: i

info = 'all'
i = len_trim(str)+1

if (str.ne.'') then
    if (index(str,'(usage:input)').ne.0) then
        i = index(str,'(usage:input)')
        info = 'input'
    elseif (index(str,'(usage:station)').ne.0) then
        i = index(str,'(usage:station)')
        info = 'station'
    elseif (index(str,'(usage:hafspc)').ne.0) then
        i = index(str,'(usage:hafspc)')
        info = 'hafspc'
    elseif (index(str,'(usage:output)').ne.0) then
        i = index(str,'(usage:output)')
        info = 'output'
    elseif (index(str,'(usage:misc)').ne.0) then
        i = index(str,'(usage:misc)')
        info = 'misc'
    elseif (index(str,'(usage:none)').ne.0) then
        i = index(str,'(usage:none)')
        info = 'none'
    endif

    write(stderr,*) str(1:i-1)
    write(stderr,*)
    if (info.eq.'none') then
        write(stderr,*) 'See o92util man page for details'
        write(stderr,*)
        call error_exit(1)
    endif
endif

write(stderr,*) 'Usage: o92util ...options...'
write(stderr,*)
if (info.eq.'all'.or.info.eq.'input') then
    write(stderr,*) 'Input fault options'
    write(stderr,*) '-ffm FFMFILE         Fault file in USGS .param format'
    write(stderr,*) '-fsp FSPFILE         Fault file in SRCMOD FSP format'
    write(stderr,*) '-mag MAGFILE         Fault file in "psmeca -Sa" format (...mag)'
    write(stderr,*) '-flt FLTFILE         Fault file with slip and dimensions (...slip wid len)'
    write(stderr,*) '-tns TNSFILE         Tensile source file (IN DEVELOPMENT)'
    write(stderr,*) '-fn|-pt              Treat faults as finite rectangular (default) or point'
    write(stderr,*) '-empirical OPT       Empirical scaling relation'
    write(stderr,*) '-thr THR             Set low slip to zero'
    write(stderr,*)
endif
if (info.eq.'all'.or.info.eq.'station') then
    write(stderr,*) 'Input target/receiver options'
    write(stderr,*) '-sta STAFILE         Station/receiver locations'
    write(stderr,*) '-auto h DEPTH N      Generate horizontal location grid'
    write(stderr,*) '-auto v AZ N         Generate vertical location grid (through centroid)'
    write(stderr,*) '-auto:thr DISP       Displacement threshold for auto grids'
    write(stderr,*) '-trg TRGFILE         Target/receiver geometry'
    write(stderr,*)
endif
if (info.eq.'all'.or.info.eq.'hafspc') then
    write(stderr,*) 'Input half-space options'
    write(stderr,*) '-haf HAFSPCFILE      Elastic half-space properties'
    write(stderr,*)
endif
if (info.eq.'all'.or.info.eq.'output') then
    write(stderr,*) 'Output options'
    write(stderr,*) '-disp DSPFILE        Displacement (E N Z)'
    write(stderr,*) '-strain STNFILE      Strain matrix (EE NN ZZ EN EZ NZ)'
    write(stderr,*) '-rotation ROTFILE    Rotation matrix (EE NN ZZ EN EZ NZ)'
    write(stderr,*) '-stress STSFILE      Stress matrix (EE NN ZZ EN EZ NZ)'
    write(stderr,*) '-estress ESTSFILE    Effective (maximum) shear stress'
    write(stderr,*) '-normal NORFILE      Normal traction on target faults (requires -trg)'
    write(stderr,*) '-shear SHRFILE       Shear traction on target faults (requires -trg)'
    write(stderr,*) '-coul COULFILE       Coulomb stress on target faults (requires -trg)'
    write(stderr,*) '--keep-all-lines     Keep comment (#), segment (>), and blank lines'
    write(stderr,*)
endif
if (info.eq.'all'.or.info.eq.'misc') then
    write(stderr,*) 'Miscellaneous options'
    write(stderr,*) '-parallel [NTHREADS] Calculate deformation in parallel'
    write(stderr,*) '-geo|-xy             Use geographic (default) or cartesian coordinates'
    write(stderr,*) '-az                  Displacement vector outputs (AZ HMAG Z)'
    write(stderr,*) '-gmt FILE            GMT psxy -SJ file (lon lat str len proj_wid)'
    write(stderr,*) '-prog                Turn on progress indicator'
    write(stderr,*) '-v LVL               Turn on verbose mode'
    write(stderr,*) '-debug [ROUTINE]     Turn on debugging'
    write(stderr,*)
endif
if (info.ne.'all') then
    write(stderr,*) 'Type "o92util" without any arguments to see all options'
endif
write(stderr,*) 'See o92util man page for details'
write(stderr,*)

call error_exit(1)
end subroutine
