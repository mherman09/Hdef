!--------------------------------------------------------------------------------------------------!
! FLTINV
!
! Program for inverting geodetic observations (and more) to determine fault slip (and more).
!
! Input datasets:
!     - 3-D displacement/velocity vectors (GPS)
!     - Line-of-sight displacement/velocity (InSAR)
!     - Data/model covariances
!     - Shear tractions on faults
!     - Model parameters/Green's functions
!
! Output quantities:
!     - Fault slip
!     - Locked faults
!     - Euler poles (rigid body rotations)
!     - Data misfit
!
! Inversion algorithms:
!     - Linear least-squares
!     - Simulated annealing
!
! Regularization options:
!     - Damping (minimizing model length)
!     - Smoothing (minimizing model roughness)
!
!--------------------------------------------------------------------------------------------------!

module fltinv
!----
! The fltinv module contains global control and data variables, as well as generic I/O subroutines
! for the program.
!----

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

! Euler pole parameters
character(len=512) :: euler_file                 !
character(len=512) :: euler_output_file          !
integer :: npoles                                !
integer, allocatable :: rigid_pt_array_disp(:)   !
integer, allocatable :: rigid_pt_array_los(:)    !
double precision, allocatable :: pole_array(:,:) !
double precision, allocatable :: gf_euler(:,:)   !

! Observation constraints on fault slip
type(fltinv_data) :: displacement                ! Three-component displacement observations
character(len=3) :: disp_components              ! Components of 3-D displacement to invert for
character(len=512) :: disp_misfit_file           ! RMS misfit to 3-D displacements
type(fltinv_data) :: los                         ! Line-of-sight displacement observations
type(fltinv_data) :: prestress                   ! Stress field at each sub-fault
character(len=512) :: cov_file                   ! File with data covariance data
double precision, allocatable :: cov_matrix(:,:) ! Data covariance matrix
character(len=16) :: input_disp_unit             !

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

! Output variables
double precision, allocatable :: fault_slip(:,:)
double precision, allocatable :: euler_pole(:,:)

! double precision :: stress_weight                ! Weight to stress data in inversion

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
character(len=512) :: input_line


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
    read(21,'(A)',iostat=ios,err=1002,end=1002) input_line
    read(input_line,*,iostat=ios) (val%array(i,j),j=1,val%ncols)
    1002 if (ios.ne.0) then
        write(stderr,*) 'read_fltinv_data_file: read error on file ',trim(val%file),' at line ',i
        write(stderr,*) 'offending line: ',trim(input_line)
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
!----
! The main fltinv program only contains a few calls to subroutines that do all the heavy lifting.
!----

#ifndef USE_LAPACK
    use io, only: stderr
#endif

use fltinv, only: output_file, &
                  euler_output_file, &
                  inversion_mode

implicit none


! Check whether LAPACK libraries are linked during compiling; fltinv will not do much without them.
#ifndef USE_LAPACK
    write(stderr,*) 'fltinv does not do much compiled without LAPACK libraries'
    write(stderr,*) 'In fact, all it does is print this message and exit'
    stop
#endif


! Parse command line arguments
call gcmdln()
if (output_file.eq.''.and.euler_output_file.eq.'') then
    call usage('fltinv: output file name is required')
endif
if (inversion_mode.eq.'') then
    call usage('fltinv: inversion mode is required')
endif


! Run the main fltinv subroutines
call read_inputs()      ! In fltinv_io.f90
call calc_disp_gfs()    ! In fltinv_gf.f90
call calc_stress_gfs()  ! In fltinv_gf.f90
call calc_euler_gfs()   ! In fltinv_gf.f90
call run_inversion()    ! In fltinv.f90 (below)
call write_solution()   ! In fltinv_io.f90

end program




!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!----------------------------------- INVERSION SUBROUTINES ----------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

subroutine run_inversion()
!----
! Run the selected inversion algorithm. Much like the main fltinv program, this subroutine is
! essentially a wrapper for starting an inversion subroutine.
!----

use io, only: stdout, stderr, verbosity
use solver, only: load_array, load_constraints

use fltinv, only: inversion_mode, &
                  euler_file

implicit none


if (verbosity.ge.1) then
    write(stdout,*) 'run_inversion: starting'
endif


! Run the selected inversion algorithm
if (inversion_mode.eq.'lsqr') then
    call invert_lsqr()                         ! In fltinv_lsqr.f90
elseif (inversion_mode.eq.'anneal') then
    call invert_anneal()                       ! In fltinv_anneal.f90
elseif (inversion_mode.eq.'anneal-psc') then
    if (euler_file.ne.'none') then
        call invert_anneal_euler_psc()         ! In fltinv_anneal_euler_psc.f90
    else
        call invert_anneal_psc()               ! In fltinv_anneal_psc.f90
    endif
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
!--------------------------------------------------------------------------------------------------!
!--------------------------------------- MISFIT ROUTINES ------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

subroutine misfit_rms(obs,pre,n,rms)
!----
! Calculate the root-mean-squared misfit between an observed and predicted set of points:
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
!----
! Calculate the chi-squared misfit between an observed and predicted set of points:
!
!     chi2 = dif_tr * cov_mat^(-1) * dif
!
! where
!
!     dif = obs-pre
!     dif_tr = transpose(dif)
!     cov_mat = covariance matrix
!-----

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
