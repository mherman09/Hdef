!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!----- PROGRAM CHECKING ROUTINES ------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

subroutine check_inputs()
!----
! Check for proper input specifications
!----
use command_line, only : fault_file, displacement_file, prestress_file, disp_comp
use io
implicit none

! Fault file is required
if (fault_file.eq.'none') then
    write(stderr,*) '!! Error: Input fault geometry file unspecified'
    call usage('!! Use -f FAULT_FILE to specify input file')
endif

! Need one of displacements or pre-stresses to invert for slip
if (displacement_file.eq.'none'.and.prestress_file.eq.'none') then
    write(stderr,*) '!! Error: No observation file specified'
    write(stderr,*) '!! Use -d DISP_FILE to specify displacement input file'
    call usage('!! Or -s PRESTS_FILE to specify pre-stress input file')
endif

! Verify displacement component variable makes sense
if (disp_comp.ne.'123'.and. &
        disp_comp.ne.'12'.and.disp_comp.ne.'13'.and.disp_comp.ne.'23'.and. &
        disp_comp.ne.'1'.and.disp_comp.ne.'2'.and.disp_comp.ne.'3') then
    write(stderr,*) '!! Error: disp_comp (',disp_comp,') must be one of ',&
                    '123, 12, 13, 23, 1, 2, or 3'
    call usage('!! Error: disp_comp set incorrectly')
endif

return
end

!--------------------------------------------------------------------------------------------------!

subroutine check_problem_parameters()
!----
! Check that problem parameters are properly defined
!----
use command_line
use io
use arrays, only : nfaults, ndisplacements
implicit none
! Local variables
integer :: n

! n: number of slip components to invert for
if (rake_file.eq.'none') then
    n = 2
else
    n = 1
endif

! Check linear least squares problem dimensions
if (trim(inversion_mode).eq.'linear') then
    if (n*nfaults.gt.len_trim(disp_comp)*ndisplacements) then
        if (damping_constant.gt.0.0d0) then
            ! okay!
        elseif (smoothing_constant.gt.0.0d0) then
            ! okay!
        elseif (prestress_file.ne.'none') then
            ! okay!
        else
            write(stderr,*) '!! Error: Number of inversion DOFs is larger than displacement input constraints'
            write(stderr,*) '!! Decrease number of sub-faults, increase number of displacements, or'
            call usage('!! add constraints: damping, smoothing, or stress minimization with pre-stresses')
        endif
    endif
else
    if (2*nfaults.gt.3*ndisplacements) then
        write(stderr,*) '!! Warning: Number of inversion DOFs is larger than displacement input constraints'
        write(stderr,*) '!! Consider decreasing number of sub-faults, increasing number of displacements, or'
        write(stderr,*) '!! adding constraints: damping, smoothing, or stress minimization with pre-stresses'
    endif
endif

return
end

!--------------------------------------------------------------------------------------------------!

subroutine check_results()
!----
! Compare inverted fault slip to observed surface displacements
!----
use command_line, only : displacement_file, disp_comp
use arrays, only : nfaults, fault_slip_nice, ndisplacements, displacements
implicit none
! Local variables
integer :: i, j, k
double precision :: fault_slip_magnitude, disp, disp_max
double precision, parameter :: factor = 100.0d0

if (displacement_file.eq.'none') then
    return
endif

! Compute maximum observed displacement
disp_max = 0.0d0
do i = 1,ndisplacements
    disp = 0.0d0
    do j = 1,len_trim(disp_comp)
        read(disp_comp(j:j),*) k
        disp = disp + displacements(i,3+k)*displacements(i,3+k)
    enddo
    if (disp.gt.disp_max) then
        disp_max = disp
    endif
enddo

! do i = 1,nfaults
!     fault_slip_magnitude = fault_slip_nice(i,1)*fault_slip_nice(i,1) &
!                              + fault_slip_nice(i,2)*fault_slip_nice(i,2)
!     if (fault_slip_magnitude.gt.mxdisp*100.0d0) then
!     endif
! !              write(0,'(A)') '!! Warning: inverted displacements are '//
! !     1                       'very large compared to the observed '//
! !     2                       'displacements.'
! !              write(0,'(A)') '!! You might want to increase the '//
! !     1                       'damping or smoothing constants'
! !              goto 912
! !          endif
! !  911 continue
! enddo

return
end

!--------------------------------------------------------------------------------------------------!

subroutine check_file_exist(file_name)
!----
! Check that a file exists
!----
implicit none
! I/O variables
character(len=*) :: file_name
! Local variables
logical :: ex

inquire(file=file_name,exist=ex)
if (.not.ex) then
    call usage('!! Error: no file found named '//trim(file_name))
endif
return
end

!--------------------------------------------------------------------------------------------------!

subroutine line_count(nlines,file_name)
!----
! Count the number of lines in a file
!----
implicit none
! I/O variables
character(len=*) :: file_name
integer :: nlines
! Local variables
integer :: ios

open(unit=41,file=file_name,status='old')
nlines= 0
do
    read(41,*,iostat=ios)
    if (ios.eq.0) then
        nlines = nlines + 1
    else
        exit
    endif
enddo
close(41)
return
end

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!----- READ DATA FROM FILES -----------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

subroutine read_faults()
!----
! Read fault locations and geometries
!----
use command_line, only : fault_file, verbosity
use arrays, only : nfaults, faults
use io
! Local variables
integer :: i, j, ios

if (verbosity.ge.2) then
    write(stderr,'("Starting read_faults()")')
endif

! Check that the fault file exists, count the number of lines
call check_file_exist(fault_file)
call line_count(nfaults,fault_file)

! Allocate memory for fault array (x, y, z, str, dip, wid, len)
if (.not.allocated(faults)) then
    allocate(faults(nfaults,7))
endif

! Read the file
open(unit=22,file=fault_file,status='old')
do i = 1,nfaults
    ! File format: x y z str dip wid len (SI units: meters, degrees)
    read(22,*,iostat=ios) (faults(i,j),j=1,7)
    if (ios.ne.0) then
        call usage('!! Read error in subroutine read_faults()')
    endif
enddo
close(22)

if (verbosity.ge.2) then
    write(stderr,'("Finished read_faults()")')
endif
if (verbosity.ge.3) then
    write(stderr,'("nfaults: ",I5)') nfaults
    write(stderr,'(7A14)') 'x','y','z','str','dip','wid','len'
    do i = 1,nfaults
        write(stderr,'(7F14.4)') (faults(i,j),j=1,7)
    enddo
endif
if (verbosity.ge.2) then
    write(stderr,*)
endif

return
end

!--------------------------------------------------------------------------------------------------!

subroutine read_displacements()
!----
! Read displacement locations and vectors
!----
use command_line, only : displacement_file, verbosity
use arrays, only : displacements, ndisplacements
use io
implicit none
! Local variables
integer :: i, j, ios

if (verbosity.ge.2) then
    write(stderr,'("Starting read_displacements()")')
endif

if (displacement_file.eq.'none') then
    if (verbosity.ge.2) then
        write(stderr,'(A)') 'User specifies no displacements to read'
        write(stderr,'(A)') 'Finished read_displacements()'
        write(stderr,*)
    endif
    return
endif

! Check that the displacement file exists, count the number of lines
call check_file_exist(displacement_file)
call line_count(ndisplacements,displacement_file)

! Allocate memory for displacement array (x, y, z, ux, uy, uz)
if (.not.allocated(displacements)) then
    allocate(displacements(ndisplacements,6))
endif

! Read the file
open(unit=21,file=displacement_file,status='old')
do i = 1,ndisplacements
    ! File format: x y z ux uy uz (SI units: meters)
    read(21,*,iostat=ios) (displacements(i,j),j=1,6)
    if (ios.ne.0) then
        call usage('!! Read error in subroutine read_displacements()')
    endif
enddo
close(21)

if (verbosity.ge.2) then
    write(stderr,'("Finished read_displacements()")')
endif
if (verbosity.ge.3) then
    write(stderr,'("ndisplacements: ",I5)') ndisplacements
    write(stderr,'(6A14)') 'x','y','z','ux','uy','uz'
    do i = 1,ndisplacements
        write(stderr,'(6F14.4)') (displacements(i,j),j=1,6)
    enddo
endif
if (verbosity.ge.2) then
    write(stderr,*)
endif

return
end

!--------------------------------------------------------------------------------------------------!

subroutine read_prestresses()
!----
! Read pre-stresses on input faults
!----
use command_line, only : prestress_file, verbosity
use arrays, only : nfaults, prestresses
use io
implicit none
! Local variables
integer :: i, j, n, ios

if (verbosity.ge.2) then
    write(stderr,'("Starting read_prestresses()")')
endif

if (prestress_file.eq.'none') then
    if (verbosity.ge.2) then
        write(stderr,'(A)') 'User specifies no pre-stresses to read'
        write(stderr,'(A)') 'Finished read_prestresses()'
        write(stderr,*)
    endif
    return
endif

! Check that the pre-stress file exists, count the number of lines, and verify
! the number of lines is the same as the number of faults
call check_file_exist(prestress_file)
call line_count(n,prestress_file)
if (n.ne.nfaults) then
    call usage('!! Error: pre-stress file must have same number of lines as fault file')
endif

! Allocate memory for prestress array (sxx, syy, szz, sxy, sxz, syzß)
if (.not.allocated(prestresses)) then
    allocate(prestresses(nfaults,6))
endif

! Read the file
open(unit=23,file=prestress_file,status='old')
do i = 1,nfaults
    ! File format: xx yy zz xy xz yz (SI units: Pa)
    read(23,*,iostat=ios) (prestresses(i,j),j=1,6)
    if (ios.ne.0) then
        call usage('!! Read error in subroutine read_prestresses()')
    endif
enddo
close(23)

if (verbosity.ge.2) then
    write(stderr,'("Finished read_prestresses()")')
endif
if (verbosity.ge.3) then
    write(stderr,'("nprestresses: ",I5)') nfaults
    write(stderr,'(6A14)') 'xx','yy','zz','xy','xz','yz'
    do i = 1,nfaults
        write(stderr,'(6E14.6)') (prestresses(i,j),j=1,6)
    enddo
endif
if (verbosity.ge.2) then
    write(stderr,*)
endif

return
end

!--------------------------------------------------------------------------------------------------!

subroutine read_halfspace()
!----
! Read half-space elastic moduli
!----
use command_line, only : halfspace_file, verbosity
use gf, only : gf_vp, gf_vs, gf_dens
use io
implicit none
! Local variables
character(len=256) :: line, moduli
double precision :: lambda, shear_modulus, poisson_ratio, young_modulus
integer :: ios

if (verbosity.ge.2) then
    write(stderr,'("Starting read_halfspace()")')
endif

! Use default values if undefined
if (halfspace_file.eq.'none') then
    if (verbosity.ge.2) then
        write(stderr,'(A)') 'User specifies no half-space information'
        write(stderr,'(A)') 'Using default half-space moduli'
    endif
    gf_vp = 6800.0d0
    gf_vs = 3926.0d0
    gf_dens = 2800.0d0
    if (verbosity.ge.3) then
        write(stderr,'(3A14)') 'vp','vs','dens'
        write(stderr,'(3F14.4)') gf_vp, gf_vs, gf_dens
    endif
    if (verbosity.ge.2) then
        write(stderr,*)
    endif
    return
endif

! Check if half-space file exists
call check_file_exist(halfspace_file)

! Parse the half-space file
open(unit=23,file=halfspace_file,status='old')
read(23,'(A)') line

! First, try reading in default format (vp vs dens)
read(line,*,iostat=ios) gf_vp, gf_vs, gf_dens
if (ios.eq.0) then
    if (verbosity.ge.2) then
        write(stderr,'("Finished read_halfspace()")')
    endif
    if (verbosity.ge.3) then
        write(stderr,'(3A14)') 'vp','vs','dens'
        write(stderr,'(3F14.4)') gf_vp, gf_vs, gf_dens
    endif
    if (verbosity.ge.2) then
        write(stderr,*)
    endif
    return
endif

! If that fails, parse the format of the line
read(line,*) moduli

! "lambda-shear" lambda shear_modulus
if (trim(moduli).eq.'lambda-shear') then
    read(line,*) moduli, lambda, shear_modulus
    gf_dens = 3.0d3
    gf_vp = dsqrt((lambda+2.0d0*shear_modulus)/gf_dens)
    gf_vs = dsqrt(shear_modulus/gf_dens)

! "shear-poisson" shear_modulus poisson_ratio
elseif (trim(moduli).eq.'shear-poisson') then
    read(line,*) moduli, shear_modulus, poisson_ratio
    lambda = 2.0d0*shear_modulus*poisson_ratio/(1.0d0-2.0d0*poisson_ratio)
    gf_dens = 3.0d3
    gf_vp = dsqrt((lambda+2.0d0*shear_modulus)/gf_dens)
    gf_vs = dsqrt(shear_modulus/gf_dens)

! "young-poisson" young_modulus poisson_ratio
elseif (trim(moduli).eq.'young-poisson') then
    read(line,*) moduli, young_modulus, poisson_ratio
    lambda = young_modulus*poisson_ratio/&
                 ((1.0d0+poisson_ratio)*(1.0d0-2.0d0*poisson_ratio))
    shear_modulus = 0.25d0*(young_modulus - 3.0d0*poisson_ratio + &
                             dsqrt(young_modulus*young_modulus + &
                                     9.0d0*poisson_ratio*poisson_ratio + &
                                     2.0d0*young_modulus*poisson_ratio))
    gf_dens = 3.0d3
    gf_vp = dsqrt((lambda+2.0d0*shear_modulus)/gf_dens)
    gf_vs = dsqrt(shear_modulus/gf_dens)

else
    call usage('!! Error: No half-space option '//trim(moduli))
endif
close(23)

if (verbosity.ge.2) then
    write(stderr,'("Finished read_halfspace()")')
endif
if (verbosity.ge.3) then
    write(stderr,'(3A14)') 'vp','vs','dens'
    write(stderr,'(3F14.4)') gf_vp, gf_vs, gf_dens
endif
if (verbosity.ge.2) then
    write(stderr,*)
endif

return
end

!--------------------------------------------------------------------------------------------------!

subroutine read_rakes()
!----
! Read in information constraining fault rakes to constant value(s)
!----
use command_line, only : verbosity, rake_file
use arrays, only : nfaults, rakes
use io
implicit none
! Local variables
integer :: i, n, ios

if (verbosity.ge.2) then
    write(stderr,'(A)') 'Starting read_rakes()'
endif

if (rake_file.eq.'none') then
    if (verbosity.ge.2) then
        write(stderr,'(A)') 'User specifies no rake constraints'
        write(stderr,'(A)') 'Finished read_rakes()'
        write(stderr,*)
    endif
    return
endif

! Allocate memory for rake constraint array (one per fault)
if (.not.allocated(rakes)) then
    allocate(rakes(nfaults))
endif

! Try to read rake as a real number first
read(rake_file,*,iostat=ios) rakes(1)

! If not a real, try to read as a file
if (ios.ne.0) then
    call check_file_exist(rake_file)
    call line_count(n,rake_file)
    if (n.ne.nfaults.and.n.ne.1) then
        call usage('!! Error: rake file must have 1 line or same number of lines as fault file')
    endif
    open(unit=24,file=rake_file,status='old')
    do i = 1,n
        read(24,*) rakes(i)
    enddo
    close(24)
else
    n = 1
endif

! If only one rake value is given, set it for all the faults
if (n.eq.1) then
    do i = 2,nfaults
        rakes(i) = rakes(1)
    enddo
endif

if (verbosity.ge.2) then
    write(stderr,'(A)') 'Finished read_rakes()'
endif
if (verbosity.ge.3) then
    write(stderr,'("nrakes: ",I10)') n
    write(stderr,'(A14)') 'rake'
    do i = 1,n
        write(stderr,'(1PE14.6)') (rakes(i))
    enddo
endif
if (verbosity.ge.2) then
    write(stderr,*)
endif

return
end

!--------------------------------------------------------------------------------------------------!

subroutine read_smoothing()
!----
! Read in file defining fault neighbors for Laplacian smoothing
!----
use command_line, only : verbosity, smoothing_file
use arrays, only : nfaults, nsmooth, smoothing_pointers, smoothing_neighbors
use io, only : stderr
implicit none
! Local variables
integer :: i, j, ineighbor, nneighbor
character(len=1024) :: iline

if (verbosity.ge.2) then
    write(stderr,'(A)') 'Starting read_smoothing()'
endif

if (smoothing_file.eq.'none') then
    if (verbosity.ge.2) then
        write(stderr,'(A)') 'User specifies no smoothing'
        write(stderr,'(A)') 'Finished read_smoothing()'
        write(stderr,*)
    endif
    return
endif

! Check that the smoothing file exists and count the number of lines
! There can be at most one set of links for each fault
call check_file_exist(smoothing_file)
call line_count(nsmooth,smoothing_file)
if (nsmooth.gt.nfaults) then
    call usage('!! Error: number of lines in smoothing file is larger than nfaults')
endif

! First, allocate memory for smoothing reference array (ifault nneighbors ptr_neighbor_array)
if (.not.allocated(smoothing_pointers)) then
    allocate(smoothing_pointers(nsmooth,3))
endif

! Read the smoothing file the first time
open(unit=25,file=smoothing_file,status='old')
do i = 1,nsmooth
    read(25,'(A)') iline
    ! First two fields: ifault nneighbors
    read(iline,*) (smoothing_pointers(i,j),j=1,2)
    ! Compute pointer to the location of the fault info in the smoothing_neighbors array
    if (i.eq.1) then
        smoothing_pointers(i,3) = 1
    else
        smoothing_pointers(i,3) = smoothing_pointers(i-1,3) + smoothing_pointers(i-1,2)
    endif
enddo
rewind(25)

! Now, allocate memory for neighbors
if (.not.allocated(smoothing_neighbors)) then
    allocate(smoothing_neighbors(smoothing_pointers(nsmooth,3)+smoothing_pointers(nsmooth,2)))
endif

! Read the smoothing file again and load the neighbors
do i = 1,nsmooth
    read(25,'(A)') iline
    nneighbor = smoothing_pointers(i,2)
    ineighbor = smoothing_pointers(i,3)
    read(iline,*) j,j,(smoothing_neighbors(ineighbor+j-1),j=1,nneighbor)
enddo
close(25)

if (verbosity.ge.2) then
    write(stderr,'("Finished read_smoothing()")')
endif
if (verbosity.ge.3) then
    write(stderr,'("nsmooth: ",I10)') nsmooth
    write(stderr,'(4A14)') 'fault','nneighbors','pointer','neighbors'
    do i = 1,nsmooth
        write(stderr,'(10I14)') (smoothing_pointers(i,j),j=1,3), &
            (smoothing_neighbors(smoothing_pointers(i,3)+j-1),j=1,smoothing_pointers(i,2))
    enddo
endif
if (verbosity.ge.2) then
    write(stderr,*)
endif

return
end

!--------------------------------------------------------------------------------------------------!

subroutine read_slip_constraints()
!----
! Read in file defining constraints on slip
!----
use command_line, only : verbosity, slip_constraint_file
use arrays, only : nfaults, nconstraints, slip_constraints, is_this_fault_constrained
use io, only : stderr
implicit none
! Local variables
integer :: i, j, ios
character(len=1024) :: iline

if (verbosity.ge.2) then
    write(stderr,'(A)') 'Starting read_slip_constraints()'
endif

if (slip_constraint_file.eq.'none') then
    if (verbosity.ge.2) then
        write(stderr,'(A)') 'User specifies no constraints on slip'
        write(stderr,'(A)') 'Finished read_slip_constraints()'
        write(stderr,*)
    endif
    return
endif

! Check that the slip constraint file exists and count the number of lines
! There can be at most one set of constraints for each fault
call check_file_exist(slip_constraint_file)
call line_count(nconstraints,slip_constraint_file)
if (nconstraints.gt.nfaults) then
    call usage('!! Error: number of constraints is larger than nfaults')
endif

! Allocate memory for slip constraint array (ss_constraint ds_constraint)
if (.not.allocated(slip_constraints)) then
    allocate(slip_constraints(nfaults,2))
endif
slip_constraints = 0.0d0

! Allocate memory for slip constraint indicator array (ss? ds?)
if (.not.allocated(is_this_fault_constrained)) then
    allocate(is_this_fault_constrained(nfaults,2))
endif
is_this_fault_constrained = 0

! Read slip constraint file
open(unit=25,file=slip_constraint_file,status='old')
do i = 1,nconstraints
    read(25,'(A)') iline
    read(iline,*) j
    if (j.gt.nfaults) then
        call usage('!! Error: read in constraint greater than nfaults')
    endif

    ! Read input (ifault [ss_]slip [ds_slip])
    read(iline,*,iostat=ios) j,slip_constraints(j,1),slip_constraints(j,2)
    is_this_fault_constrained(j,:) = 1
    if (ios.ne.0) then
        read(iline,*,iostat=ios) j,slip_constraints(j,1)
        is_this_fault_constrained(j,1) = 1
    endif
    if (ios.ne.0) then
        call usage('!! Read error in read_slip_constraints()')
    endif
enddo

if (verbosity.ge.2) then
    write(stderr,'("Finished read_slip_constraints()")')
endif
if (verbosity.ge.3) then
    write(stderr,'("nconstraints: ",I10)') nconstraints
    write(stderr,'(3A14)') 'fault','ss_constraint','ds_constraint'
    do i = 1,nfaults
        write(stderr,'(I14,2(1PE14.6))') i,slip_constraints(i,1),slip_constraints(i,2)
    enddo
endif
if (verbosity.ge.2) then
    write(stderr,*)
endif

return
end

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!----- STRESS/STRAIN MATRIX ROUTINES --------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

subroutine resolve_prestresses()
!----
! Compute strike-slip and dip-slip components of shear traction from pre-stresses
! resolved onto input faults
!----
use arrays, only : faults, prestresses, nfaults
use command_line, only : verbosity
use io
implicit none
! Local variables
integer :: i
double precision, dimension(3,3) :: stress
double precision                 :: strike, dip
double precision, dimension(3)   :: traction, traction_shear
double precision                 :: traction_normal
double precision, dimension(3)   :: normal_vector, strike_vector, updip_vector

if (verbosity.ge.1) then
    write(stderr,'("Starting resolve_prestresses()")')
endif

do i = 1,nfaults
    stress(1,1) = prestresses(i,1)
    stress(2,2) = prestresses(i,2)
    stress(3,3) = prestresses(i,3)
    stress(1,2) = prestresses(i,4)
    stress(2,1) = prestresses(i,4)
    stress(1,3) = prestresses(i,5)
    stress(3,1) = prestresses(i,5)
    stress(2,3) = prestresses(i,6)
    stress(3,2) = prestresses(i,6)
    strike = faults(i,4)
    dip = faults(i,5)

    ! Unit vectors corresponding to fault geometry
    call unit_normal(normal_vector,strike,dip)
    call unit_strike(strike_vector,strike)
    call unit_updip(updip_vector,strike,dip)

    ! Resolve traction on fault and compute shear component of traction
    call mat_vec_mul(traction,stress,normal_vector)
    call dotproduct(traction_normal,traction,normal_vector)
    call vec_add(traction_shear,traction,-traction_normal*normal_vector)

    ! Resolve shear tractions on along-strike and up-dip directions
    call dotproduct(prestresses(i,1),traction_shear,strike_vector)
    call dotproduct(prestresses(i,2),traction_shear,updip_vector)
enddo

if (verbosity.ge.1) then
    write(stderr,'("Finished resolve_prestresses()")')
endif
if (verbosity.ge.3) then
    write(stderr,'(1A4,2A14)') 'flt','prests_strike','prests_updip'
    do i = 1,nfaults
        write(stderr,'(1I4,2E14.6)') i,prestresses(i,1),prestresses(i,2)
    enddo
endif
if (verbosity.ge.1) then
    write(stderr,*)
endif
return
end

!--------------------------------------------------------------------------------------------------!

subroutine strain2stress(stress,strain,vp,vs,dens)
!----
! Calculate (3x3) stress matrix from (3x3) strain matrix for isotropic, elastic material
!----
implicit none
! I/O variables
double precision, dimension(3,3) :: stress, strain
double precision :: vp, vs, dens
! Local variables
double precision :: shear_modulus, lambda, diag

shear_modulus  = dens*vs*vs
lambda = dens*vp*vp - 2.0d0*shear_modulus
if (shear_modulus.lt.10.0e7) shear_modulus = 10.0e7
if (lambda.lt.10.0e7) lambda = 10.0e7

! Trace of strain matrix
diag = strain(1,1) + strain(2,2) + strain(3,3)

stress(1,1) = lambda*diag + 2.0d0*shear_modulus*strain(1,1)
stress(2,2) = lambda*diag + 2.0d0*shear_modulus*strain(2,2)
stress(3,3) = lambda*diag + 2.0d0*shear_modulus*strain(3,3)
stress(1,2) = 2.0d0*shear_modulus*strain(1,2)
stress(1,3) = 2.0d0*shear_modulus*strain(1,3)
stress(2,3) = 2.0d0*shear_modulus*strain(2,3)
stress(2,1) = stress(1,2)
stress(3,1) = stress(1,3)
stress(3,2) = stress(2,3)

return
end

!--------------------------------------------------------------------------------------------------!

subroutine rotstrain(strain,strike)
!----
! Rotate strain matrix from (x=str, y=updip horizontal, z=up) to (x=E, y=N, z=up)
!----
use trig
implicit none
! I/O variables
double precision, dimension(3,3) :: strain
double precision :: strike
! Local variables
double precision, dimension(3,3) :: rot, rot_transpose, tmp
rot(1,1) = dcos(d2r*strike-pi/2.0d0)
rot(2,2) = rot(1,1)
rot(1,2) = dsin(d2r*strike-pi/2.0d0)
rot(2,1) = -rot(1,2)
rot(1,3) = 0.0d0
rot(2,3) = 0.0d0
rot(3,1) = 0.0d0
rot(3,2) = 0.0d0
rot(3,3) = 1.0d0
call mattr(rot_transpose,rot)
call matmult(tmp,rot,strain)
call matmult(strain,tmp,rot_transpose)
return
end

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!----- FAULT GEOMETRY SUBROUTINES -----------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

subroutine unit_normal(normal_vector,strike_in,dip_in)
!----
! Unit vector normal to plane defined by strike and dip (x=E, y=N, z=Z)
!----
use trig
implicit none
! I/O variables
double precision, dimension(3) :: normal_vector
double precision :: strike_in, dip_in
! Local variables
double precision :: strike, dip
strike = strike_in*d2r
dip = dip_in*d2r
normal_vector(1) = dsin(dip)*dsin(strike+pi/2.0d0)
normal_vector(2) = dsin(dip)*dcos(strike+pi/2.0d0)
normal_vector(3) = dcos(dip)
return
end

!--------------------------------------------------------------------------------------------------!

subroutine unit_strike(strike_vector,strike_in)
!----
! Unit vector in along-strike direction (x=E, y=N, z=Z)
!----
use trig
implicit none
! I/O variables
double precision, dimension(3) :: strike_vector
double precision :: strike_in
! Local variables
double precision :: strike
strike = strike_in*d2r
strike_vector(1) = dsin(strike)
strike_vector(2) = dcos(strike)
strike_vector(3) = 0.0d0
return
end

!--------------------------------------------------------------------------------------------------!

subroutine unit_updip(updip_vector,strike_in,dip_in)
!----
! Unit vector in up-dip direction (x=E, y=N, z=Z)
!----
use trig
implicit none
! I/O variables
double precision, dimension(3) :: updip_vector
double precision :: strike_in, dip_in
! Local variables
double precision :: strike, dip
strike = strike_in*d2r
dip = dip_in*d2r
updip_vector(1) = dcos(dip)*dsin(strike-pi/2.0d0)
updip_vector(2) = dcos(dip)*dcos(strike-pi/2.0d0)
updip_vector(3) = dsin(dip)
return
end

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!----- GREENS FUNCTIONS ROUTINES ------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

subroutine calc_disp_gfs()
!----
! Compute displacement Green's functions
!----
use command_line, only : verbosity, rake_file, displacement_file
use arrays, only : ndisplacements, displacements, nfaults, faults, disp_gfs, rakes
use gf
use io
implicit none
! Local variables
integer :: i, j, k

if (verbosity.ge.2) then
    write(stderr,'(A)') 'Starting calc_disp_gfs()'
endif

if (displacement_file.eq.'none') then
    if (verbosity.ge.2) then
        write(stderr,'(A)') 'Not using displacements'
        write(stderr,'(A)') 'Finished calc_disp_gfs()'
        write(stderr,*)
    endif
    return
endif

if (.not.allocated(disp_gfs)) then
    allocate(disp_gfs(ndisplacements,nfaults,6))
endif

! Unit slip for GF computation
gf_slip = 1.0d0

! GF for each flt-sta pair
do i = 1,ndisplacements
    gf_stlo = displacements(i,1)
    gf_stla = displacements(i,2)
    gf_stdp = displacements(i,3)
    do j = 1,nfaults
        gf_evlo = faults(j,1)
        gf_evla = faults(j,2)
        gf_evdp = faults(j,3)
        gf_str  = faults(j,4)
        gf_dip  = faults(j,5)
        gf_wid  = faults(j,6)
        gf_len  = faults(j,7)
        if (rake_file.eq.'none') then
            ! Strike-slip Green's function
            gf_rak  = 0.0d0
            call calc_disp_okada(disp_gfs(i,j,1),disp_gfs(i,j,2),disp_gfs(i,j,3))
            ! Dip-slip Green's function
            gf_rak  = 90.0d0
            call calc_disp_okada(disp_gfs(i,j,4),disp_gfs(i,j,5),disp_gfs(i,j,6))
        else
            gf_rak = rakes(j)
            call calc_disp_okada(disp_gfs(i,j,1),disp_gfs(i,j,2),disp_gfs(i,j,3))
        endif
    enddo
enddo

if (verbosity.ge.2) then
    write(stderr,'("Finished calc_disp_gfs()")')
endif
if (verbosity.ge.3) then
    write(stderr,'(2A4,6A14)') 'sta','flt','dgf_ss_x','dgf_ss_y','dgf_ss_z',&
                                           'dgf_ds_x','dgf_ds_y','dgf_ds_z'
    do i = 1,ndisplacements
        do j = 1,nfaults
            write(stderr,'(2I4,6F14.6)') i,j,(disp_gfs(i,j,k),k=1,6)
        enddo
    enddo
endif
if (verbosity.ge.2) then
    write(stderr,*)
endif
return
end

!--------------------------------------------------------------------------------------------------!

subroutine calc_disp_okada(ux,uy,uz)
!----
! Compute static displacement vector for set of Okada (1992) variables
! All units are SI; angles inputs are in degrees
!----
use trig
use gf
! Local variables
double precision :: ux, uy, uz
double precision :: delta_x, delta_y, dist, az, x, y
double precision :: uxp, uyp, theta, uhor

! Distance and azimuth from source to station
delta_x = gf_stlo - gf_evlo
delta_y = gf_stla - gf_evla
dist = dsqrt(delta_x*delta_x+delta_y*delta_y)
az   = datan2(delta_x,delta_y) ! clockwise from north

! Rotate to x=along-strike, y=horizontal up-dip
x = dist*( dcos(az-d2r*gf_str))
y = dist*(-dsin(az-d2r*gf_str))

! Compute displacement
call o92rect(uxp,uyp,uz,x,y,gf_stdp,gf_evdp,gf_dip,gf_rak, &
             gf_wid,gf_len,gf_slip,gf_vp,gf_vs,gf_dens)

! Rotate back to original x-y coordinates
theta = datan2(uyp,uxp)
uhor = dsqrt(uxp*uxp+uyp*uyp)
theta = d2r*gf_str - theta
ux = uhor*dsin(theta)
uy = uhor*dcos(theta)
return
end

!--------------------------------------------------------------------------------------------------!

subroutine calc_stress_gfs()
!----
! Compute shear stress Green's functions
!----
use command_line, only : verbosity, rake_file, prestress_file
use arrays, only : nfaults, faults, stress_gfs, rakes
use gf
use io
implicit none
! Local variables
double precision :: ststr, stdip
integer :: i, j, k

if (verbosity.ge.2) then
    write(stderr,'(A)') 'Starting calc_stress_gfs()'
endif

if (prestress_file.eq.'none') then
    if (verbosity.ge.2) then
        write(stderr,'(A)') 'Not using pre-stresses'
        write(stderr,'(A)') 'Finished calc_stress_gfs()'
        write(stderr,*)
    endif
    return
endif

call resolve_prestresses()

if (.not.allocated(stress_gfs)) then
    allocate(stress_gfs(nfaults,nfaults,4))
endif

! Unit slip for GF computation
gf_slip = 1.0d0

! Shear stress Green's functions for each fault-fault pair
do i = 1,nfaults
    gf_stlo = faults(i,1)
    gf_stla = faults(i,2)
    gf_stdp = faults(i,3)
    ststr = faults(i,4)
    stdip = faults(i,5)
    do j = 1,nfaults
        gf_evlo = faults(j,1)
        gf_evla = faults(j,2)
        gf_evdp = faults(j,3)
        gf_str  = faults(j,4)
        gf_dip  = faults(j,5)
        gf_wid  = faults(j,6)
        gf_len  = faults(j,7)
        ! Avoid singular solutions with measurement point lying on fault
        if (dabs(gf_evlo-gf_stlo).lt.1.0d0) then
            gf_stlo = gf_evlo + 1.0d0
        endif
        if (dabs(gf_evla-gf_stla).lt.1.0d0) then
            gf_stla = gf_evla + 1.0d0
        endif
        if (dabs(gf_evdp-gf_stdp).lt.1.0d0) then
            gf_stdp = gf_evdp + 1.0d0
        endif
        if (rake_file.eq.'none') then
            ! Shear stress Green's function produced by strike-slip source
            gf_rak  = 0.0d0
            call calc_shear_okada(stress_gfs(i,j,1),stress_gfs(i,j,2),ststr,stdip)
            ! Shear stress Green's function produced by dip-slip source
            gf_rak  = 90.0d0
            call calc_shear_okada(stress_gfs(i,j,3),stress_gfs(i,j,4),ststr,stdip)
        else
            gf_rak  = rakes(j)
            call calc_shear_okada(stress_gfs(i,j,1),stress_gfs(i,j,2),ststr,stdip)
        endif
    enddo
enddo

if (verbosity.ge.1) then
    write(stderr,'("Finished calc_stress_gfs()")')
endif
if (verbosity.ge.3) then
    write(stderr,'(2A4,4A14)') 'sta','flt','sgf_ss_ss','sgf_ss_ds',&
                                           'sgf_ds_ss','sgf_ds_ds'
    do i = 1,nfaults
        do j = 1,nfaults
            write(stderr,'(2I4,1P4E14.6)') i,j,(stress_gfs(i,j,k),k=1,4)
        enddo
    enddo
endif
if (verbosity.ge.1) then
    write(stderr,*)
endif
return
end

!--------------------------------------------------------------------------------------------------!

subroutine calc_shear_okada(shear_ss,shear_ds,ststr,stdip)
!----
! Compute strike-slip and dip-slip components of shear stress
!----
use gf
implicit none
! I/O variables
double precision :: shear_ss, shear_ds, ststr, stdip
! Local variables
double precision, dimension(3,3) :: strain, stress
double precision, dimension(3)   :: normal_vector, strike_vector, updip_vector
double precision, dimension(3)   :: traction, traction_shear
double precision                 :: traction_normal

! Compute stress tensor from strain tensor
call calc_strain(strain)
call strain2stress(stress,strain,gf_vp,gf_vs,gf_dens)

! Unit vectors corresponding to fault geometry
call unit_normal(normal_vector,ststr,stdip)
call unit_strike(strike_vector,ststr)
call unit_updip(updip_vector,ststr,stdip)

! Resolve traction on fault and compute shear component of traction
call mat_vec_mul(traction,stress,normal_vector)
call dotproduct(traction_normal,traction,normal_vector)
call vec_add(traction_shear,traction,-traction_normal*normal_vector)

! Resolve shear tractions on along-strike and up-dip directions
call dotproduct(shear_ss,traction_shear,strike_vector)
call dotproduct(shear_ds,traction_shear,updip_vector)

return
end

!--------------------------------------------------------------------------------------------------!

subroutine calc_strain(strain)
use gf
use trig
implicit none
! I/O variables
double precision, dimension(*) :: strain
! Local variables
double precision :: delx, dely, dist, az, x, y

! Distance and azimuth from source to station
     delx = gf_stlo - gf_evlo
     dely = gf_stla - gf_evla
     dist = dsqrt(delx*delx+dely*dely)
     az   = datan2(delx,dely) ! clockwise from north

! Rotate to x=along-strike, y=horizontal up-dip
     x = dist*( dcos(az-d2r*gf_str))
     y = dist*(-dsin(az-d2r*gf_str))

! Compute strain in fault-centered coordinates
call o92rectstn(strain,x,y,gf_stdp,gf_evdp,gf_dip,gf_rak,gf_wid,gf_len,gf_slip,gf_vp,gf_vs,gf_dens)

! Rotate back to original x-y coordinates
call rotstrain(strain,gf_str)
return
end

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!----- LINEAR ALGEBRA ROUTINES --------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

subroutine vec_add(vector_out,vector_in_1,vector_in_2)
!----
! Add two 3x1 vectors
!----
implicit none
! I/O variables
double precision, dimension(3) :: vector_out, vector_in_1, vector_in_2
! Local variables
integer :: i
do i = 1,3
    vector_out(i) = vector_in_1(i) + vector_in_2(i)
enddo
return
end

!--------------------------------------------------------------------------------------------------!

subroutine dotproduct(scalar,vector1,vector2)
!----
! Dot product of two 3x1 vectors
!----
implicit none
! I/O variables
double precision, dimension(3) :: vector1, vector2
double precision :: scalar
! Local variables
integer :: i
scalar = 0.0d0
do i = 1,3
    scalar = scalar + vector1(i)*vector2(i)
enddo
return
end

!--------------------------------------------------------------------------------------------------!

subroutine mat_vec_mul(vector_out,matrix_in,vector_in)
!----
! Multiply a 3x3 matrix with a 3x1 vector
!----
implicit none
! I/O variables
double precision, dimension(3) :: vector_out, vector_in
double precision, dimension(3,3) :: matrix_in
! Local variables
integer :: i, j
vector_out = 0.0d0
do i = 1,3
    do j = 1,3
        vector_out(i) = vector_out(i) + matrix_in(i,j)*vector_in(j)
    enddo
enddo
return
end

!--------------------------------------------------------------------------------------------------!

subroutine matmult(matout,mat1,mat2)
!----
! Multiply two 3x3 matrices
!----
implicit none
! I/O variables
double precision, dimension(3,3) :: mat1, mat2, matout
! Local variables
integer :: i, j, k
do i = 1,3
    do j = 1,3
        matout(i,j) = 0.0d0
        do k = 1,3
            matout(i,j) = matout(i,j) + mat1(i,k)*mat2(k,j)
        enddo
    enddo
enddo
return
end

!--------------------------------------------------------------------------------------------------!

subroutine mattr(matout,matin)
!----
! Transpose a 3x3 matrix
!----
implicit none
! I/O variables
double precision, dimension(3,3) :: matout, matin
! Local variables
integer :: i, j
do i = 1,3
    do j = 1,3
        matout(i,j) = matin(j,i)
    enddo
enddo
return
end

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!----- GEOGRAPHIC SUBROUTINES ---------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

subroutine geo2xy(n,array)
!----
! Assuming first two values in array are longitude and latitude, convert them to Cartesian coords
!----
implicit none
! I/O variables
integer :: n
double precision, dimension(:,:) :: array
! Local variables
integer :: i
double precision :: lon0, lat0, lon, lat, dist, az

! Set origin to first coordinate
lon0 = array(1,1)
lat0 = array(1,2)

! Compute distance and azimuth from origin
do i = 1,n
    lon = array(i,1)
    lat = array(i,2)
    call ddistaz(dist,az,lon0,lat0,lon,lat)
    dist = dist*6.371d6
    array(i,1) = dist*dsin(az)
    array(i,2) = dist*dcos(az)
enddo

return
end

!----------------------------------------------------------------------C
!
   !   SUBROUTINE model0(soln,gf,obs,nobs,OBSMAX,nflt,FLTMAX,inityp,
!     1                  slip0,rake0)
!!----
!! Set up the initial solution
!!----
!      IMPLICIT none
!      REAL*8 pi,d2r
!      PARAMETER (pi=4.0d0*atan(1.0d0),d2r=pi/1.8d2)
!      INTEGER nobs,OBSMAX,nflt,FLTMAX,inityp,i,j
!      REAL*8 obs(OBSMAX,6),gf(OBSMAX,FLTMAX,6),soln(FLTMAX,2)
!      REAL*8 pre(3),op,pp,slip0,rake0,ss,ds
!      if (inityp.eq.0) then
!          slip0 = 0.0d0
!          rake0  = 0.0d0
!      elseif (inityp.eq.1) then
!          ss = dcos(rake0*d2r)
!          ds = dsin(rake0*d2r)
!          op = 0.0d0 ! sum of dot product obs.pre
!          pp = 0.0d0 ! sum of dot product pre.pre
!          do 572 i = 1,nobs
!              do 571 j = 1,nflt
!                  pre(1) = ss*gf(i,j,1) + ds*gf(i,j,4)
!                  pre(2) = ss*gf(i,j,2) + ds*gf(i,j,5)
!                  pre(3) = ss*gf(i,j,3) + ds*gf(i,j,6)
!                  op = op + obs(i,4)*pre(1) + obs(i,5)*pre(2) +
!     1                                                   obs(i,6)*pre(3)
!                  pp = pp +   pre(1)*pre(1) +   pre(2)*pre(2) +
!     1                                                     pre(3)*pre(3)
!  571         continue
!  572     continue
!          slip0 = op/pp
!      elseif (inityp.eq.2) then
!          ! slip0 and rake0 already defined
!      else
!          call usage('!! Error: no annealing option '//char(inityp))
!      endif
!      do 573 i = 1,nflt
!          soln(i,1) = slip0
!          soln(i,2) = rake0
!  573 continue
!      RETURN
!      END
!
!----------------------------------------------------------------------C
!
    !  SUBROUTINE calcobj(obj,misfit,length,rough,soln,gf,obs,nobs,
!     1                   OBSMAX,nflt,FLTMAX,damp,smooth,smarry,nsmoo)
!      IMPLICIT none
!      REAL*8 obj
!      REAL*8 misfit,length,rough,smooth,damp
!      INTEGER nobs,OBSMAX,nflt,FLTMAX
!      REAL*8 obs(OBSMAX,6),gf(OBSMAX,FLTMAX,6),soln(FLTMAX,2)
!      INTEGER nsmoo,smarry(FLTMAX,8)
!      call calcmisfit(misfit,soln,gf,obs,nobs,OBSMAX,nflt,FLTMAX)
!      misfit = dsqrt(misfit/nflt) ! RMS misfit
!      length = 0.0d0
!      rough = 0.0d0
!      if (damp.gt.0.0d0) then
!          call calclength(length,soln,nflt,FLTMAX)
!          length = damp*length/nflt ! Average length
!      endif
!      if (smooth.gt.0.0d0) then
!          call calcrough(rough,soln,FLTMAX,smarry,nsmoo)
!          rough = smooth*rough/nflt ! Average roughness
!      endif
!      obj = misfit + damp*length + smooth*rough
!      RETURN
!      END
!
!----------------------------------------------------------------------C
!
   !   SUBROUTINE calcmisfit(misfit,soln,gf,obs,nobs,OBSMAX,nflt,FLTMAX)
!      IMPLICIT none
!      REAL*8 pi,d2r
!      PARAMETER (pi=4.0d0*datan(1.0d0),d2r=pi/1.8d2)
!      INTEGER nobs,OBSMAX,nflt,FLTMAX
!      REAL*8 obs(OBSMAX,6),gf(OBSMAX,FLTMAX,6)
!      REAL*8 soln(FLTMAX,2)
!      REAL*8 pre(3)
!      REAL*8 misfit
!      INTEGER i,j
!      REAL*8 cr,sr,ss(FLTMAX),ds(FLTMAX)
!! Store strike-slip and dip-slip components of faults
!      do 581 i = 1,nflt
!          cr = dcos(soln(i,2)*d2r)
!          sr = dsin(soln(i,2)*d2r)
!          ss(i) = soln(i,1)*cr
!          ds(i) = soln(i,1)*sr
!  581 continue
!! Compute misfit
!      misfit = 0.0d0
!      do 583 i = 1,nobs
!          pre(1) = 0.0d0
!          pre(2) = 0.0d0
!          pre(3) = 0.0d0
!          do 582 j = 1,nflt
!              pre(1) = pre(1) + gf(i,j,1)*ss(j) + gf(i,j,4)*ds(j)
!              pre(2) = pre(2) + gf(i,j,2)*ss(j) + gf(i,j,5)*ds(j)
!              pre(3) = pre(3) + gf(i,j,3)*ss(j) + gf(i,j,6)*ds(j)
!  582     continue
!          misfit = misfit + (obs(i,4)-pre(1))*(obs(i,4)-pre(1))
!     1                    + (obs(i,5)-pre(2))*(obs(i,5)-pre(2))
!     2                    + (obs(i,6)-pre(3))*(obs(i,6)-pre(3))
!  583 continue
!      RETURN
!      END
!
!----------------------------------------------------------------------C
!
    !  SUBROUTINE calclength(length,soln,nflt,FLTMAX)
!      IMPLICIT none
!      INTEGER nflt,FLTMAX,i
!      REAL*8 soln(FLTMAX,2),length
!      length = 0.0d0
!      do 584 i = 1,nflt
!          length = length + soln(i,1)
!  584 continue
!      RETURN
!      END
!
!----------------------------------------------------------------------C
!
   !   SUBROUTINE calcrough(rough,soln,FLTMAX,smarry,nsmoo)
!      IMPLICIT none
!      REAL*8 pi,d2r
!      PARAMETER (pi=4.0d0*datan(1.0d0),d2r=pi/1.8d2)
!      INTEGER FLTMAX,nsmoo,smarry(FLTMAX,8),i,j,jj
!      REAL*8 soln(FLTMAX,2),rough,r(2),ss,ds
!      r(1) = 0.0d0
!      r(2) = 0.0d0
!      do 586 i = 1,nsmoo
!          ss = soln(i,1)*dcos(soln(i,2)*d2r)
!          ds = soln(i,1)*dsin(soln(i,2)*d2r)
!          r(1) = r(1) + smarry(i,2)*ss
!          r(2) = r(2) + smarry(i,2)*ds
!          do 585 j = 1,smarry(i,2)
!              jj = smarry(i,j+2)
!              ss = soln(jj,1)*dcos(soln(jj,2)*d2r)
!              ds = soln(jj,1)*dsin(soln(jj,2)*d2r)
!              r(1) = r(1) - ss
!              r(2) = r(2) - ds
!  585     continue
!  586 continue
!      rough = r(1) + r(2)
!      RETURN
!      END
!
!----------------------------------------------------------------------C
!
   !   SUBROUTINE chksol(x,nflt)
!      IMPLICIT none
!      INTEGER OBSMAX
!      PARAMETER (OBSMAX=2500)
!      REAL*8 x(OBSMAX,1)
!      INTEGER nflt,i
!      do 611 i = 1,nflt
!          if (x(i,1).gt.100.0d0) then
!              write(*,*) '!! Error: found slip larger than 100 m'
!              write(*,*) '!! Check inputs and units'
!              write(*,*) '!! Try adding damping/smoothing'
!              stop
!          endif
!  611 continue
!      RETURN
!      END
!
!--------------------------------------------------------------------------------------------------!

subroutine write_output()
use command_line
use io, only : stdout
use arrays, only : nfaults, fault_slip, fault_slip_nice, &
                   is_this_fault_constrained, slip_constraints
implicit none
! Local variables
integer :: i, n, output_unit

if (.not.allocated(fault_slip_nice)) then
    allocate(fault_slip_nice(nfaults,2))
endif

n = 0
do i = 1,nfaults
    if (allocated(is_this_fault_constrained)) then
        if (is_this_fault_constrained(i,1).eq.0) then
            n = n + 1
            fault_slip_nice(i,1) = fault_slip(n)
        else
            fault_slip_nice(i,1) = slip_constraints(i,1)
        endif
    else
        fault_slip_nice(i,1) = fault_slip(i)
    endif
enddo

if (rake_file.eq.'none') then
    do i = 1,nfaults
        if (allocated(is_this_fault_constrained)) then
            if (is_this_fault_constrained(i,2).eq.0) then
                n = n + 1
                fault_slip_nice(i,2) = fault_slip(n)
            else
                fault_slip_nice(i,2) = slip_constraints(i,2)
            endif
        else
            fault_slip_nice(i,2) = fault_slip(i+nfaults)
        endif
    enddo
endif

! Check that the results are reasonable
call check_results()

if (output_file.eq.'stdout') then
    output_unit = stdout
else
    output_unit = 61
    open(unit=output_unit,file=output_file,status='unknown')
endif

do i = 1,nfaults
    if (rake_file.eq.'none') then
        write(output_unit,'(1P2E14.6)') fault_slip_nice(i,1),fault_slip_nice(i,2)
    else
        write(output_unit,'(1P1E14.6)') fault_slip_nice(i,1)
    endif
enddo

if (output_file.ne.'stdout') then
    close(output_unit)
endif

return
end
