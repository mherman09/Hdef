module triutil

character(len=512) :: fault_file       ! source faults (x1 y1 z1 x2 y2 z2 x3 y3 z3 ss ds ts)
character(len=512) :: station_file     ! station locations (stlo stla stdp)
character(len=512) :: target_file      ! target/receiver fault geometry (str dip rak fric)
character(len=512) :: halfspace_file         ! half-space parameters
character(len=512) :: displacement_file      ! displacement
character(len=512) :: strain_file      ! strain tensor
character(len=512) :: stress_file      ! stress tensor
character(len=512) :: estress_file      ! effective shear stress (maximum shear stress)
character(len=512) :: normal_file      ! normal stress
character(len=512) :: shear_file      ! shear stress
character(len=512) :: coulomb_file     ! Coulomb stress
logical :: isFaultFileDefined
logical :: isStationFileDefined
logical :: isTargetFileDefined
logical :: isOutputFileDefined
logical :: iWantDisp
logical :: iWantStrain
logical :: iWantStress
logical :: iWantTraction
character(len=16) :: coord_type         ! cartesian-m, cartesian-km, geographic
integer :: nfaults
integer :: nstations
double precision, allocatable :: faults(:,:)
double precision, allocatable :: stations(:,:)
double precision, allocatable :: targets(:,:)
double precision :: poisson, lame, shearmod

end module

!==================================================================================================!

program main
!----
! Utility for computing displacements, strains, and stresses in an elastic half-space resulting
! from triangular shear dislocations. All the work is done in tri_disloc_module.f90, which is
! simply a modification of Brendan Meade's Matlab code into Fortran.
!
! Please cite: Meade, B.J. (2007). Algorithms for the calculation of exact displacements, strains,
! and stresses for Triangular Dislocation Elements in a uniform elastic half space. Computers &
! Geosciences. doi:10.1016/j.cageo.2006.12.003.
!----

use triutil, only: isOutputFileDefined
implicit none

call gcmdln()
if (.not.isOutputFileDefined) then
    call usage('triutil: no output defined')
endif

call read_faults()
call read_stations()
call read_targets()
call read_halfspace()
call calc_deformation()

end

!--------------------------------------------------------------------------------------------------!

subroutine read_faults()
!----
! Read the input fault file, in format:
!     x1 y1 z1 x2 y2 z2 x3 y3 z3 strike_slip dip_slip [tensile_slip]
!----

use io, only: stderr, line_count
use triutil, only: fault_file, nfaults, faults, isFaultFileDefined
implicit none

! Local variables
integer :: i, j, ios
character(len=512) :: input_line
logical :: doesFaultFileExist

! Initialize ios
ios = 0

! Check that fault file is defined and exists
if (.not.isFaultFileDefined) then
    call usage('read_faults: no fault file defined')
endif
inquire(exist=doesFaultFileExist,file=fault_file)
if (.not.doesFaultFileExist) then
    call usage('read_faults: no fault file named '//trim(fault_file)//' found')
endif

! Count the number of lines in the file and allocate memory for the fault array
! Fault array: x1 y1 z1 x2 y2 z2 x3 y3 z3 ss ds ts
nfaults = line_count(fault_file)
allocate(faults(nfaults,12))

! Read fault file
open(unit=12,file=fault_file,status='old')
do i = 1,nfaults
    read(12,'(A)') input_line

    ! Attempt to read 12 inputs, assuming tensile_slip is provided
    read(input_line,*,iostat=ios,err=1001) (faults(i,j),j=1,12)

    ! Attempt to read 11 inputs, assuming tensile_slip is not included
    if (ios.ne.0) then
        read(input_line,*,iostat=ios,err=1001,end=1002) (faults(i,j),j=1,11)
        faults(i,12) = 0.0d0
    endif
enddo
close(12)

! Read error messages
1001 if (ios.ne.0) then
    write(stderr,*) 'read_faults: read error'
    call usage('offending line: '//trim(input_line))
endif
1002 if (ios.ne.0) then
    write(stderr,*) 'read_faults: line end'
    call usage('offending line: '//trim(input_line))
endif

return
end subroutine read_faults

!--------------------------------------------------------------------------------------------------!

subroutine read_stations()
!----
! Read the input station file, in format:
!     x1 y1 z1
!----

use io, only: stderr, line_count
use triutil, only: station_file, nstations, stations, isStationFileDefined
implicit none

! Local variables
integer :: i, j, ios
character(len=512) :: input_line
logical :: doesStationFileExist

! Check that station file is defined and the file exists
if (.not.isStationFileDefined) then
    call usage('read_stations: no station file defined')
endif
inquire(exist=doesStationFileExist,file=station_file)
if (.not.doesStationFileExist) then
    call usage('read_stations: no station file named '//trim(station_file)//' found')
endif

! Count number of stations and allocate memory to station array (x y z)
nstations = line_count(station_file)
allocate(stations(nstations,3))

! Read stations
open(unit=13,file=station_file,status='old')
do i = 1,nstations
    read(13,'(A)') input_line
    read(input_line,*,iostat=ios,err=1003,end=1004) (stations(i,j),j=1,3)
enddo
close(13)

! Error messages
1003 if (ios.ne.0) then
    write(stderr,*) 'read_stations: read error'
    call usage('offending line: '//trim(input_line))
endif
1004 if (ios.ne.0) then
    write(stderr,*) 'read_stations: line end'
    call usage('offending line: '//trim(input_line))
endif

return
end subroutine read_stations

!--------------------------------------------------------------------------------------------------!

subroutine read_targets()
!----
! Read the input target fault geometry file, in format:
!     str dip rak fric
!----

use io, only: stderr, line_count
use triutil, only: target_file, targets, isTargetFileDefined, iWantTraction, nstations
implicit none

! Local variables
integer :: i, j, ios, ntargets
character(len=512) :: input_line
logical :: doesTargetFileExist

! Check that tractions need to be resolved, target fault file is defined, and the file exists
if (.not.iWantTraction) then
    return
else
    if (.not.isTargetFileDefined) then
        call usage('read_targets: no target geometry file defined')
    endif
endif
inquire(exist=doesTargetFileExist,file=target_file)
if (.not.doesTargetFileExist) then
    call usage('read_targets: no target geometry file named '//trim(target_file)//' found')
endif

! Count number of target geometries and allocate memory to target array (str dip rak fric)
ntargets = line_count(target_file)
if (ntargets.ne.1.and.ntargets.ne.nstations) then
    call usage('read_targets: number of target geometries must be 1 or nstations')
endif
allocate(targets(nstations,4))

! Read target fault geometries
open(unit=14,file=target_file,status='old')
do i = 1,ntargets
    read(14,'(A)') input_line
    read(input_line,*,iostat=ios,err=1005,end=1006) (targets(i,j),j=1,4)
enddo
close(14)

! Fill target array if only one is specified
if (ntargets.eq.1) then
    do i = 2,nstations
        targets(i,:) = targets(1,:)
    enddo
endif

! Error messages
1005 if (ios.ne.0) then
    write(stderr,*) 'read_targets: read error'
    call usage('offending line: '//trim(input_line))
endif
1006 if (ios.ne.0) then
    write(stderr,*) 'read_targets: line end'
    call usage('offending line: '//trim(input_line))
endif

return
end subroutine read_targets

!--------------------------------------------------------------------------------------------------!

subroutine read_halfspace()
!----
! Read half-space elastic moduli from file (if provided), in format:
!     lbl1 value1 lbl2 value2 [lbl3 value3]
!----

use io, only: stderr, line_count
use triutil, only: halfspace_file, poisson, shearmod, lame
implicit none

! Local variables
integer :: i, j
character(len=512) :: input_line
logical :: doesHalfspaceFileExist
character(len=32) :: lbl(3)
double precision :: val(3)

! If not defined, use defaults from gcmdln()
if (halfspace_file.eq.'') then
    return
endif

! Check whether half-space file exists
inquire(exist=doesHalfspaceFileExist,file=halfspace_file)
if (.not.doesHalfspaceFileExist) then
    write(stderr,*) 'read_halfspace: no half-space file named '//trim(halfspace_file)//' found'
    write(stderr,*) 'using default elastic moduli values'
    return
endif

! Initialize modulus variables
lbl = ''
val = 0.0d0

! Open the file and read the half-space elastic moduli
open(unit=15,file=halfspace_file,status='old')
read(15,'(A)') input_line
read(input_line,*,iostat=ios) lbl(1), val(1), lbl(2), val(2), lbl(3), val(3)
if (ios.ne.0) then
    read(input_line,*,iostat=ios) lbl(1), val(1), lbl(2), val(2)
endif
if (ios.ne.0) then
    write(stderr,*) 'read_halfspace: error reading modulus labels and values'
    call usage('offending line: '//trim(input_line))
endif
close(15)

do i = 1,2
    if (lbl(i).eq.'mu' .or. lbl(i).eq.'shear_modulus' .or. lbl(i).eq.'shear' .or. &
                                                                     lbl(i).eq.'shearmod') then
        lbl(i) = 'shearmod'
    elseif (lbl(i).eq.'lame') then
        lbl(i) = 'lame'
    elseif (lbl(i).eq.'poisson' .or. lbl(i).eq.'nu') then
        lbl(i) = 'poisson'
    elseif (lbl(i).eq.'young') then
        lbl(i) = 'young'
    elseif (lbl(i).eq.'vp'.or.lbl(i).eq.'vs'.or.lbl(i).eq.'dens') then
        ! Okay
    else
        write(stderr,*) 'read_halfspace: elastic modulus "',trim(lbl(i)),'" not implemented'
        call usage('specify two of: shear_modulus, lame, poisson, or young OR vp vs dens')
    endif
enddo

! Subroutines use poisson to compute displacement and strain, lame and shearmod to compute stress
do i = 1,2

    j = mod(i,2)+1

    if (lbl(i).eq.'shearmod'.and.lbl(j).eq.'lame') then
        ! Calculate Poisson's ratio from Lame's parameter and shear modulus
        shearmod = val(i)
        lame = val(j)
        poisson = lame/(2.0d0*(lame+shearmod))

    elseif (lbl(i).eq.'shearmod'.and.lbl(j).eq.'poisson') then
        ! Calculate Lame's parameter from Poisson's ratio and shear modulus
        shearmod = val(i)
        poisson = val(j)
        lame = 2.0d0*shearmod*poisson/(1.0d0-2.0d0*poisson)

    elseif (lbl(i).eq.'shearmod'.and.lbl(j).eq.'young') then
        ! Calculate Lame's parameter and Poisson's ratio from Young's modulus and shear modulus
        shearmod = val(i)
        lame = shearmod*(val(j)-2.0d0*shearmod)/(3.0d0*shearmod-val(j))
        poisson = lame/(2.0d0*(lame+shearmod))

    elseif (lbl(i).eq.'lame'.and.lbl(j).eq.'poisson') then
        ! Calculate shear modulus from Lame's parameter and Poisson's ratio
        lame = val(i)
        poisson = val(j)
        shearmod = lame*(1.0d0-2.0d0*poisson)/(2.0d0*poisson)

    elseif (lbl(i).eq.'lame'.and.lbl(j).eq.'young') then
        ! Calculate shear modulus and Poisson's ratio from Lame's parameter and Young's modulus
        lame = val(i)
        poisson = 2.0d0*lame/(val(j)+lame+sqrt(val(j)*val(j)+9.0d0*lame*lame+2.0d0*lame*val(j)))
        shearmod = lame*(1.0d0-2.0d0*poisson)/(2.0d0*poisson)

    elseif (lbl(i).eq.'poisson'.and.lbl(j).eq.'young') then
        ! Calculate shear modulus and Lame's parameter from Poisson's ratio and Young's modulus
        poisson = val(i)
        shearmod = val(j)/(2.0d0*(1.0d0+poisson))
        lame = 2.0d0*shearmod*poisson/(1.0d0-2.0d0*poisson)

    endif
enddo

if (lbl(1).eq.'vp'.and.lbl(2).eq.'vs'.and.lbl(3).eq.'dens') then
    ! Calculate shear modulus, Lame's parameter, and Poisson's ratio from seismic velocities and density
    shearmod = val(2)*val(2)*val(3)
    lame = val(3)*val(3)*val(3) + 2.0d0*shearmod
    poisson = lame/(2.0d0*(lame+shearmod))
endif

return
end subroutine read_halfspace

!--------------------------------------------------------------------------------------------------!

subroutine calc_deformation()
!----
! Calculate all requested deformation values at station locations
!----

use triutil, only: displacement_file, &
                   strain_file, &
                   stress_file, &
                   estress_file, &
                   normal_file, &
                   shear_file, &
                   coulomb_file, &
                   iWantDisp, iWantStrain, iWantStress, iWantTraction, &
                   nstations, stations, targets, &
                   lame, shearmod
use io, only: stderr, verbosity
implicit none

! Local variables
integer :: iSta, file_unit
logical :: isThisUnitOpen
double precision :: disp(3), strain(3,3), stress(3,3), estress
double precision :: trac_vector(3), shear, shearmax, normal, coulomb

! Check which calculations are needed and open output files if requested
if (iWantDisp) then
    if (verbosity.ge.1) then
        write(stderr,*) 'calc_deformation: calculating displacements'
    endif
    if (displacement_file.ne.'') then
        open(unit=101,file=displacement_file,status='unknown')
    endif
endif
if (iWantStrain) then
    if (verbosity.ge.1) then
        write(stderr,*) 'calc_deformation: calculating strains'
    endif
    if (strain_file.ne.'') then
        open(unit=111,file=strain_file,status='unknown')
    endif
endif
if (iWantStress) then
    if (verbosity.ge.1) then
        write(stderr,*) 'calc_deformation: calculating stresses'
    endif
    if (stress_file.ne.'') then
        open(unit=121,file=stress_file,status='unknown')
    endif
    if (estress_file.ne.'') then
        open(unit=122,file=estress_file,status='unknown')
    endif
endif
if (iWantTraction) then
    if (verbosity.ge.1) then
        write(stderr,*) 'calc_deformation: calculating resolved tractions'
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

! Calculate the requested quantities at each station
do iSta = 1,nstations

    if (iWantDisp) then
        call calc_displacement(disp,stations(iSta,1),stations(iSta,2),stations(iSta,3))
        ! Displacement: ux, uy, uz
        if (displacement_file.ne.'') then
            write(101,*) stations(iSta,:),disp
        endif
    endif

    if (iWantStrain) then
        call calc_strain(strain,stations(iSta,1),stations(iSta,2),stations(iSta,3))
        ! Strain tensor: exx, eyy, ezz, exy, exz, eyz
        if (strain_file.ne.'') then
            write(111,*) stations(iSta,:),strain(1,1),strain(2,2),strain(3,3), &
                         strain(1,2),strain(1,3),strain(2,3)
        endif
    endif

    if (iWantStress) then
        call strain2stress(stress,strain,lame,shearmod)
        ! Stress tensor: sxx, syy, szz, sxy, sxz, syz
        if (stress_file.ne.'') then
            write(121,*) stations(iSta,:),stress(1,1),stress(2,2),stress(3,3), &
                         stress(1,2),stress(1,3),stress(2,3)
        endif
        ! Maximum (effective) shear stress: estress
        if (estress_file.ne.'') then
            call calc_estress(estress,stress)
            write(121,*) stations(iSta,:),estress
        endif
    endif

    if (iWantTraction) then
        call calc_tractions(trac_vector,shear,shearmax,normal,coulomb,stress,targets(iSta,:))
        ! Normal traction: normal (positive=dilation)
        if (normal_file.ne.'') then
            write(131,*) stations(iSta,:),normal
        endif
        ! Shear traction: resolved_onto_rake, max_shear_on_plane
        if (shear_file.ne.'') then
            write(132,*) stations(iSta,:),shear,shearmax
        endif
        ! Coulomb stress: coulomb
        if (coulomb_file.ne.'') then
            write(133,*) stations(iSta,:),coulomb
        endif
    endif
enddo

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

return
end subroutine calc_deformation

!--------------------------------------------------------------------------------------------------!

subroutine calc_displacement(disp,x,y,z)
!----
! Calculate displacements at coordinate (x,y,z), where z is defined positive down, but the
! displacement uz is defined positive up.
!----

use triutil, only: nfaults, faults, poisson, coord_type     ! Fault and half-space variables
use tri_disloc, only: tri_disloc_disp, tri_center           ! Triangular dislocation subroutines
implicit none

! Arguments
double precision :: disp(3), x, y, z

! Local variables
double precision :: sta_coord(3), tri_coord(3,4), slip(3), tri_coord_new(3,4)
double precision :: center(3), dist, az
double precision :: disptmp(3)
integer :: iFlt, iTri

! Initialize displacement
disp = 0.0d0

! Set station location units to meters (input units are km)
if (coord_type.eq.'cartesian') then
    sta_coord(1) = x*1.0d3
    sta_coord(2) = y*1.0d3
endif
sta_coord(3) = z*1.0d3

do iFlt = 1,nfaults
    tri_coord(1,1) = faults(iFlt,1)
    tri_coord(2,1) = faults(iFlt,2)
    tri_coord(3,1) = faults(iFlt,3)
    tri_coord(1,2) = faults(iFlt,4)
    tri_coord(2,2) = faults(iFlt,5)
    tri_coord(3,2) = faults(iFlt,6)
    tri_coord(1,3) = faults(iFlt,7)
    tri_coord(2,3) = faults(iFlt,8)
    tri_coord(3,3) = faults(iFlt,9)
    slip(1) = -faults(iFlt,10)        ! triutil: ll positive;  tri_disloc_disp(): rl positive
    slip(2) = -faults(iFlt,11)        ! triutil: thr positive; tri_disloc_disp(): nor positive
    slip(3) = faults(iFlt,12)

    if (coord_type.eq.'geographic') then
        ! Convert lon lat dep(km) to x(m) y(m) z(m) from triangle center
        call tri_center(center,tri_coord(:,1),tri_coord(:,2),tri_coord(:,3))
        call ddistaz(dist,az,center(1),center(2),x,y)
        sta_coord(1) = dist*6.371d6*dsin(az)
        sta_coord(2) = dist*6.371d6*dcos(az)
        do iTri = 1,3
            call ddistaz(dist,az,center(1),center(2),tri_coord(1,iTri),tri_coord(2,iTri))
            tri_coord_new(1,iTri) = dist*6.371d6*dsin(az)
            tri_coord_new(2,iTri) = dist*6.371d6*dcos(az)
            tri_coord_new(3,iTri) = tri_coord(3,iTri)*1.0d3
        enddo
    else
        ! Convert triangle x(km) y(km) z(km) to x(m) y(m) z(m)
        tri_coord_new = tri_coord*1.0d3
    endif

    call tri_disloc_disp(disptmp, sta_coord, tri_coord_new, poisson, slip)
    disp = disp + disptmp
enddo

! Output from tri_disloc_disp is positive down, so invert vertical component
disp(3) = -disp(3)

return
end subroutine calc_displacement

!--------------------------------------------------------------------------------------------------!

subroutine calc_strain(strain,x,y,z)
!----
! Calculate strains at coordinate (x,y,z), where z is defined positive down, but the
! strain z coordinate is defined positive up.
!----

use triutil, only: nfaults, faults, poisson, coord_type   ! Fault and half-space variables
use tri_disloc, only: tri_disloc_strain, tri_center       ! Triangular dislocation subroutines
implicit none

! Arguments
double precision :: strain(3,3), x, y, z

! Local variables
double precision :: sta_coord(3), tri_coord(3,4), slip(3), tri_coord_new(3,4)
double precision :: center(3), dist, az
double precision :: straintmp(3,3)
integer :: iFlt, iTri

! Initialize displacement
strain = 0.0d0

! Set coordinates to meters (input units are geographic and/or km)
if (coord_type.eq.'cartesian') then
    sta_coord(1) = x*1.0d3
    sta_coord(2) = y*1.0d3
endif
sta_coord(3) = z*1.0d3

do iFlt = 1,nfaults
    tri_coord(1,1) = faults(iFlt,1)
    tri_coord(2,1) = faults(iFlt,2)
    tri_coord(3,1) = faults(iFlt,3)
    tri_coord(1,2) = faults(iFlt,4)
    tri_coord(2,2) = faults(iFlt,5)
    tri_coord(3,2) = faults(iFlt,6)
    tri_coord(1,3) = faults(iFlt,7)
    tri_coord(2,3) = faults(iFlt,8)
    tri_coord(3,3) = faults(iFlt,9)
    slip(1) = -faults(iFlt,10)        ! triutil: ll positive;  tri_disloc_disp(): rl positive
    slip(2) = -faults(iFlt,11)        ! triutil: thr positive; tri_disloc_disp(): nor positive
    slip(3) = faults(iFlt,12)

    if (coord_type.eq.'geographic') then
        ! Convert lon lat dep(km) to x(m) y(m) z(m) from triangle center
        call tri_center(center,tri_coord(:,1),tri_coord(:,2),tri_coord(:,3))
        call ddistaz(dist,az,center(1),center(2),x,y)
        sta_coord(1) = dist*6.371d6*dsin(az)
        sta_coord(2) = dist*6.371d6*dcos(az)
        do iTri = 1,3
            call ddistaz(dist,az,center(1),center(2),tri_coord(1,iTri),tri_coord(2,iTri))
            tri_coord_new(1,iTri) = dist*6.371d6*dsin(az)
            tri_coord_new(2,iTri) = dist*6.371d6*dcos(az)
            tri_coord_new(3,iTri) = tri_coord(3,iTri)*1.0d3
        enddo
    else
        ! Convert triangle x(km) y(km) z(km) to x(m) y(m) z(m)
        tri_coord_new = tri_coord*1.0d3
    endif

    call tri_disloc_strain(straintmp, sta_coord, tri_coord_new, poisson, slip)
    strain = strain + straintmp
enddo

! Output from tri_disloc_disp is positive down, so invert xz and yz strain components
strain(1,3) = -strain(1,3)
strain(3,1) = -strain(1,3)
strain(2,3) = -strain(2,3)
strain(3,2) = -strain(2,3)

return
end subroutine calc_strain

!--------------------------------------------------------------------------------------------------!

subroutine strain2stress(stress,strain,lame_param,shear_modulus)
!----
! Calculate stress tensor from strain tensor, assuming linear isotropic medium
!----

implicit none

! Arguments
double precision :: stress(3,3), strain(3,3), lame_param, shear_modulus

! Local variables
double precision :: diag

diag = strain(1,1) + strain(2,2) + strain(3,3)

! Diagonal terms
stress(1,1) = lame_param*diag + 2.0d0*shear_modulus*strain(1,1)
stress(2,2) = lame_param*diag + 2.0d0*shear_modulus*strain(2,2)
stress(3,3) = lame_param*diag + 2.0d0*shear_modulus*strain(3,3)

! Off-diagonal terms
stress(1,2) = 2.0d0*shear_modulus*strain(1,2)
stress(1,3) = 2.0d0*shear_modulus*strain(1,3)
stress(2,3) = 2.0d0*shear_modulus*strain(2,3)
stress(2,1) = stress(1,2)
stress(3,1) = stress(1,3)
stress(3,2) = stress(2,3)

return
end subroutine strain2stress

!--------------------------------------------------------------------------------------------------!

subroutine calc_estress(estress,stress)
!----
! Calculate maximum shear stress
!----

implicit none

! Arguments
double precision :: estress, stress(3,3)

! Local variables
double precision :: s11_s22, s11_s33, s22_s33, s12s12, s13s13, s23s23

s11_s22 = stress(1,1)-stress(2,2)
s11_s33 = stress(1,1)-stress(3,3)
s22_s33 = stress(2,2)-stress(3,3)

s12s12 = stress(1,2)*stress(1,2)
s13s13 = stress(1,3)*stress(1,3)
s23s23 = stress(2,3)*stress(2,3)

estress = (1.0d0/6.0d0)*(s11_s22*s11_s22+s11_s33*s11_s33+s22_s33*s22_s33) + s12s12+s13s13+s23s23

! Make sure it has units of Pa
estress = dsqrt(estress)

return
end subroutine calc_estress

!--------------------------------------------------------------------------------------------------!

subroutine calc_tractions(trac_vector,shear,shearmax,normal,coulomb,stress,trg)
!----
! Calculate various tractions resolved onto a fault geometry from a stress tensor
!----

use trig, only: pi, d2r
implicit none

! Arguments
double precision :: trac_vector(3), shear, shearmax, normal, coulomb, stress(3,3), trg(4)

! Local variables
integer :: i
double precision :: str, dip, rak, fric
double precision :: n(3), s(3), r(3)

str = trg(1)*d2r
dip = trg(2)*d2r
rak = trg(3)*d2r
fric = trg(4)

! Unit normal vector to target plane
n(1) = dsin(dip)*dsin(str+pi/2.0d0)
n(2) = dsin(dip)*dcos(str+pi/2.0d0)
n(3) = dcos(dip)

! Traction vector is stress matrix times normal vector: t = S*n
do i = 1,3
    trac_vector(i) = stress(i,1)*n(1) + stress(i,2)*n(2) + stress(i,3)*n(3)
enddo

! Normal component of traction is parallel to unit normal vector; take dot product
normal = 0.0d0
do i = 1,3
    normal = normal + trac_vector(i)*n(i)
enddo

! Shear component of traction is difference between total traction vector and normal traction
s(1) = trac_vector(1) - normal*n(1)
s(2) = trac_vector(2) - normal*n(2)
s(3) = trac_vector(3) - normal*n(3)
shearmax = dsqrt(s(1)*s(1) + s(2)*s(2) + s(3)*s(3))

! Compute unit slip vector (parallel to rake)
r(1) =  dcos(str-pi/2.0d0)*dcos(rak) + dsin(str-pi/2.0d0)*dsin(rak)*dcos(dip)
r(2) = -dsin(str-pi/2.0d0)*dcos(rak) + dcos(str-pi/2.0d0)*dsin(rak)*dcos(dip)
r(3) =                                                    dsin(rak)*dsin(dip)

! Project shear component on slip vector
shear = s(1)*r(1) + s(2)*r(2) + s(3)*r(3)

! Coulomb stress (recall sign convention: pos = dilation)
coulomb = shear + fric*normal

return
end subroutine calc_tractions

!--------------------------------------------------------------------------------------------------!

subroutine gcmdln()

use io, only: stderr, verbosity
use triutil, only: fault_file, &
                   station_file, &
                   target_file, &
                   halfspace_file, &
                   displacement_file, &
                   strain_file, &
                   stress_file, &
                   estress_file, &
                   normal_file, &
                   shear_file, &
                   coulomb_file, &
                   coord_type, &
                   isFaultFileDefined, &
                   isStationFileDefined, &
                   isTargetFileDefined, &
                   isOutputFileDefined, &
                   iWantDisp, &
                   iWantStrain, &
                   iWantStress, &
                   iWantTraction, &
                   poisson, &
                   lame, &
                   shearmod
implicit none

! Local variables
integer :: i, narg, ios
character(len=512) :: tag

! Initialize control variables
fault_file = ''
station_file = ''
target_file = ''
halfspace_file = ''
displacement_file = ''
strain_file = ''
stress_file = ''
estress_file = ''
normal_file = ''
shear_file = ''
coulomb_file = ''
isFaultFileDefined = .false.
isStationFileDefined = .false.
isTargetFileDefined = .false.
isOutputFileDefined = .false.
iWantDisp = .false.
iWantStrain = .false.
iWantStress = .false.
iWantTraction = .false.
poisson = 0.25d0
lame = 40.0d9
shearmod = 40.0d9
coord_type = 'geographic'
verbosity = 0

narg = command_argument_count()
if (narg.eq.0) call usage('')

! Everything is hunky-dory to start out
ios = 0

i = 1
do while (i.le.narg)

    call get_command_argument(i,tag)

    if (trim(tag).eq.'-flt') then
        i = i + 1
        call get_command_argument(i,fault_file)

        isFaultFileDefined = .true.

    elseif (trim(tag).eq.'-sta') then
        i = i + 1
        call get_command_argument(i,station_file)

        isStationFileDefined = .true.

    elseif (trim(tag).eq.'-trg') then
        i = i + 1
        call get_command_argument(i,target_file)

        isTargetFileDefined = .true.

    elseif (trim(tag).eq.'-haf') then
        i = i + 1
        call get_command_argument(i,halfspace_file)

    elseif (trim(tag).eq.'-disp') then
        i = i + 1
        call get_command_argument(i,displacement_file)

        isOutputFileDefined = .true.
        iWantDisp = .true.

    elseif (trim(tag).eq.'-strain') then
        i = i + 1
        call get_command_argument(i,strain_file)

        isOutputFileDefined = .true.
        iWantStrain = .true.

    elseif (trim(tag).eq.'-stress') then
        i = i + 1
        call get_command_argument(i,stress_file)

        isOutputFileDefined = .true.
        iWantStrain = .true.
        iWantStress = .true.

    elseif (trim(tag).eq.'-estress') then
        i = i + 1
        call get_command_argument(i,estress_file)

        isOutputFileDefined = .true.
        iWantStrain = .true.
        iWantStress = .true.

    elseif (trim(tag).eq.'-normal') then
        i = i + 1
        call get_command_argument(i,normal_file)

        isOutputFileDefined = .true.
        iWantStrain = .true.
        iWantStress = .true.
        iWantTraction = .true.

    elseif (trim(tag).eq.'-shear'.or.trim(tag).eq.'-shearmax') then
        i = i + 1
        call get_command_argument(i,shear_file)

        isOutputFileDefined = .true.
        iWantStrain = .true.
        iWantStress = .true.
        iWantTraction = .true.

    elseif (trim(tag).eq.'-coul') then
        i = i + 1
        call get_command_argument(i,coulomb_file)

        isOutputFileDefined = .true.
        iWantStrain = .true.
        iWantStress = .true.
        iWantTraction = .true.

    elseif (trim(tag).eq.'-geographic'.or.trim(tag).eq.'-geo') then
        coord_type = 'geographic'

    elseif (trim(tag).eq.'-cartesian'.or.trim(tag).eq.'-xy') then
        coord_type = 'cartesian'

    elseif (trim(tag).eq.'-v'.or.trim(tag).eq.'-verbose'.or.trim(tag).eq.'-verbosity') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*,iostat=ios,err=9001) verbosity

    else
        call usage('!! gcmdln: no option '//trim(tag))

    endif

    i = i + 1

enddo

if (verbosity.ge.1) then
    write(stderr,*) 'gcmdln read:'
    write(stderr,*) 'fault_file:          ',fault_file
    write(stderr,*) 'station_file:        ',station_file
    write(stderr,*) 'target_file:         ',target_file
    write(stderr,*) 'halfspace_file:      ',halfspace_file
    write(stderr,*) 'displacement_file:   ',displacement_file
    write(stderr,*) 'strain_file:         ',strain_file
    write(stderr,*) 'stress_file:         ',stress_file
    write(stderr,*) 'estress_file:        ',estress_file
    write(stderr,*) 'normal_file:         ',normal_file
    write(stderr,*) 'shear_file:          ',shear_file
    write(stderr,*) 'coulomb_file:        ',coulomb_file
    write(stderr,*) 'isFaultFileDefined:  ',isFaultFileDefined
    write(stderr,*) 'isStationFileDefined:',isStationFileDefined
    write(stderr,*) 'isTargetFileDefined: ',isTargetFileDefined
    write(stderr,*) 'isOutputFileDefined: ',isOutputFileDefined
    write(stderr,*) 'iWantDisp:           ',iWantDisp
    write(stderr,*) 'iWantStrain:         ',iWantStrain
    write(stderr,*) 'iWantStress:         ',iWantStress
    write(stderr,*) 'iWantTraction:       ',iWantTraction
    write(stderr,*) 'coord_type:          ',coord_type
    write(stderr,*) 'verbosity:           ',verbosity
endif

9001 if (ios.ne.0) then
    write(stderr,*) 'gcmdln: error reading string ',trim(tag)
    call usage('exiting triutil from gcmdln')
endif

return
end subroutine gcmdln

!--------------------------------------------------------------------------------------------------!

subroutine usage(str)

use io, only: stderr
implicit none

! Arguments
character(len=*) :: str

if (trim(str).ne.'') then
    write(stderr,*) trim(str)
    write(stderr,*)
endif

write(stderr,*) 'Usage: triutil ...options...'
write(stderr,*)
write(stderr,*) '-flt FLTFILE         Triangular source faults'
write(stderr,*) '-sta STAFILE         Station/receiver locations'
write(stderr,*) '-trg TRGFILE         Target/receiver kinematics'
write(stderr,*) '-haf HAFSPCFILE      Elastic half-space properties'
write(stderr,*) '-disp DSPFILE        Displacement (E N Z)'
write(stderr,*) '-strain STNFILE      Strain matrix (EE NN ZZ EN EZ NZ)'
write(stderr,*) '-stress STSFILE      Stress matrix (EE NN ZZ EN EZ NZ)'
write(stderr,*) '-estress ESTSFILE    Effective (maximum) shear stress'
write(stderr,*) '-normal NORFILE      Normal traction on target faults (requires -trg)'
write(stderr,*) '-shear[max] SHRFILE  Shear traction on target faults (requires -trg)'
write(stderr,*) '-coul COULFILE       Coulomb stress on target faults (requires -trg)'
write(stderr,*) '-geo|-xy             Use geographic (default) or cartesian coordinates'
write(stderr,*)

stop
end subroutine usage
