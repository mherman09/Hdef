module triutil_module

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
double precision :: poisson, lame, mu
integer :: verbosity

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

use triutil_module, only: isOutputFileDefined
implicit none

call gcmdln()
if (.not.isOutputFileDefined) then
    call usage('triutil: no output defined')
endif

call read_faults()
call read_stations()
call read_targets()
call calc_deformation()

end

!--------------------------------------------------------------------------------------------------!

subroutine read_faults()
use triutil_module, only: fault_file, nfaults, faults, isFaultFileDefined
implicit none
! Local variables
integer :: i, j, ios
character(len=512) :: input_line
logical :: doesFaultFileExist

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
call line_count(nfaults,fault_file)
allocate(faults(nfaults,12))

! Read fault file
open(unit=12,file=fault_file,status='old')
do i = 1,nfaults
    read(12,'(A)') input_line
    read(input_line,*,iostat=ios,err=1001,end=1002) (faults(i,j),j=1,12)
enddo
close(12)

1001 if (ios.ne.0) then
    write(0,*) 'read_faults: read error'
    call usage('offending line: '//trim(input_line))
endif
1002 if (ios.ne.0) then
    write(0,*) 'read_faults: line end'
    call usage('offending line: '//trim(input_line))
endif

return
end subroutine read_faults

!--------------------------------------------------------------------------------------------------!

subroutine read_stations()
use triutil_module, only: station_file, nstations, stations, isStationFileDefined
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
call line_count(nstations,station_file)
allocate(stations(nstations,3))

! Read stations
open(unit=13,file=station_file,status='old')
do i = 1,nstations
    read(13,'(A)') input_line
    read(input_line,*,iostat=ios,err=1003,end=1004) (stations(i,j),j=1,3)
enddo
close(13)

1003 if (ios.ne.0) then
    write(0,*) 'read_stations: read error'
    call usage('offending line: '//trim(input_line))
endif
1004 if (ios.ne.0) then
    write(0,*) 'read_stations: line end'
    call usage('offending line: '//trim(input_line))
endif

return
end subroutine read_stations

!--------------------------------------------------------------------------------------------------!

subroutine read_targets()
use triutil_module, only: target_file, targets, isTargetFileDefined, iWantTraction, nstations
implicit none
! Local variables
integer :: i, j, ios, ntargets
character(len=512) :: input_line
logical :: doesTargetFileExist

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
call line_count(ntargets,target_file)
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

1005 if (ios.ne.0) then
    write(0,*) 'read_targets: read error'
    call usage('offending line: '//trim(input_line))
endif
1006 if (ios.ne.0) then
    write(0,*) 'read_targets: line end'
    call usage('offending line: '//trim(input_line))
endif

return
end subroutine read_targets

!--------------------------------------------------------------------------------------------------!

subroutine calc_deformation()
!----
! Calculate all requested deformation values at stations
!----

use triutil_module, only: displacement_file, &
                          strain_file, &
                          stress_file, estress_file, &
                          normal_file, shear_file, &
                          coulomb_file, iWantDisp, iWantStrain, iWantStress, iWantTraction, &
                          nstations, stations, targets, verbosity, lame, mu
implicit none
integer :: iSta, file_unit
logical :: isThisUnitOpen
double precision :: disp(3), strain(3,3), stress(3,3), estress
double precision :: trac_vector(3), shear, shearmax, normal, coulomb

! Check which calculations are needed and open output files if requested
if (iWantDisp) then
    if (verbosity.ge.1) then
        write(0,*) 'calc_deformation: calculating displacements'
    endif
    if (displacement_file.ne.'') then
        open(unit=101,file=displacement_file,status='unknown')
    endif
endif
if (iWantStrain) then
    if (verbosity.ge.1) then
        write(0,*) 'calc_deformation: calculating strains'
    endif
    if (strain_file.ne.'') then
        open(unit=111,file=strain_file,status='unknown')
    endif
endif
if (iWantStress) then
    if (verbosity.ge.1) then
        write(0,*) 'calc_deformation: calculating stresses'
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
        write(0,*) 'calc_deformation: calculating resolved tractions'
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
        if (displacement_file.ne.'') then
            write(101,*) stations(iSta,:),disp
        endif
    endif

    if (iWantStrain) then
        call calc_strain(strain,stations(iSta,1),stations(iSta,2),stations(iSta,3))
        if (strain_file.ne.'') then
            write(111,*) stations(iSta,:),strain(1,1),strain(2,2),strain(3,3), &
                         strain(1,2),strain(1,3),strain(2,3)
        endif
    endif

    if (iWantStress) then
        call strain2stress(stress,strain,lame,mu)
        if (stress_file.ne.'') then
            write(121,*) stations(iSta,:),stress(1,1),stress(2,2),stress(3,3), &
                         stress(1,2),stress(1,3),stress(2,3)
        endif
        if (estress_file.ne.'') then
            call calc_estress(estress,stress)
            write(121,*) stations(iSta,:),estress
        endif
    endif

    if (iWantTraction) then
        call calc_tractions(trac_vector,shear,shearmax,normal,coulomb,stress,targets(iSta,:))
        if (normal_file.ne.'') then
            write(131,*) stations(iSta,:),normal
        endif
        if (shear_file.ne.'') then
            write(132,*) stations(iSta,:),shear,shearmax
        endif
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
use tri_disloc_module, only: tri_disloc_disp, tri_center
use triutil_module, only: nfaults, faults, poisson, coord_type
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
use tri_disloc_module, only: tri_disloc_strain, tri_center
use triutil_module, only: nfaults, faults, poisson, coord_type
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

    ! write(0,*) 'sta_coord:',sta_coord
    ! write(0,*) 'tri_coord_new(:,1):',tri_coord_new(:,1)
    ! write(0,*) 'tri_coord_new(:,2):',tri_coord_new(:,2)
    ! write(0,*) 'tri_coord_new(:,3):',tri_coord_new(:,3)
    ! write(0,*) 'poisson:',poisson
    ! write(0,*) 'slip:',slip
    call tri_disloc_strain(straintmp, sta_coord, tri_coord_new, poisson, slip)
    strain = strain + straintmp
enddo

! write(0,*) 'strain:',strain(1,1),strain(2,2),strain(3,3)
! write(0,*) 'strain:',strain(1,2),strain(1,3),strain(2,3)

! Output from tri_disloc_disp is positive down, so invert xz and yz strain components
strain(1,3) = -strain(1,3)
strain(3,1) = -strain(1,3)
strain(2,3) = -strain(2,3)
strain(3,2) = -strain(2,3)

return
end subroutine calc_strain

!--------------------------------------------------------------------------------------------------!

subroutine strain2stress(stress,strain,lame_param,shear_modulus)
implicit none
! Arguments
double precision :: stress(3,3), strain(3,3), lame_param, shear_modulus
! Local variables
double precision :: diag
diag = strain(1,1) + strain(2,2) + strain(3,3)
stress(1,1) = lame_param*diag + 2.0d0*shear_modulus*strain(1,1)
stress(2,2) = lame_param*diag + 2.0d0*shear_modulus*strain(2,2)
stress(3,3) = lame_param*diag + 2.0d0*shear_modulus*strain(3,3)
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
estress = dsqrt(estress)
return
end subroutine calc_estress

!--------------------------------------------------------------------------------------------------!

subroutine calc_tractions(trac_vector,shear,shearmax,normal,coulomb,stress,trg)
!----
! Calculate various tractions resolved onto a fault geometry from a stress tensor
!----

implicit none
! Arguments
double precision :: trac_vector(3), shear, shearmax, normal, coulomb, stress(3,3), trg(4)
! Local variables
integer :: i
double precision :: str, dip, rak, fric
double precision :: n(3), s(3), r(3)
double precision, parameter :: pi=datan(1.0d0)*4.0d0, d2r=pi/180.0d0

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

subroutine line_count(nlines,filename)
implicit none
integer nlines, ios
character(len=*) :: filename
open(unit=11,file=filename,status='old')
nlines = 0
do
    read(11,*,iostat=ios)
    if (ios.eq.0) then
        nlines = nlines + 1
    else
        exit
    endif
enddo
close(11)
return
end subroutine line_count

!--------------------------------------------------------------------------------------------------!

subroutine gcmdln()
use triutil_module, only: fault_file, station_file, target_file, halfspace_file, &
                          displacement_file, &
                          strain_file, &
                          stress_file, estress_file, &
                          normal_file, shear_file, coulomb_file, &
                          coord_type, &
             isFaultFileDefined, isStationFileDefined, isTargetFileDefined, isOutputFileDefined, &
             iWantDisp, iWantStrain, iWantStress, iWantTraction, poisson, lame, mu, verbosity
implicit none
! Local variables
integer :: i, narg, ios
character(len=512) :: tag

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
mu = 40.0d9
coord_type = 'geographic'
verbosity = 0

narg = command_argument_count()
if (narg.eq.0) call usage('')

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
    elseif (trim(tag).eq.'-poisson') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*,iostat=ios,err=9001) poisson
    elseif (trim(tag).eq.'-lame') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*,iostat=ios,err=9001) lame
    elseif (trim(tag).eq.'-mu') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*,iostat=ios,err=9001) mu
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
    write(0,*) 'gcmdln read:'
    write(0,*) 'fault_file:          ',fault_file
    write(0,*) 'station_file:        ',station_file
    write(0,*) 'target_file:         ',target_file
    write(0,*) 'halfspace_file:      ',halfspace_file
    write(0,*) 'displacement_file:   ',displacement_file
    write(0,*) 'strain_file:         ',strain_file
    write(0,*) 'stress_file:         ',stress_file
    write(0,*) 'estress_file:        ',estress_file
    write(0,*) 'normal_file:         ',normal_file
    write(0,*) 'shear_file:          ',shear_file
    write(0,*) 'coulomb_file:        ',coulomb_file
    write(0,*) 'isFaultFileDefined:  ',isFaultFileDefined
    write(0,*) 'isStationFileDefined:',isStationFileDefined
    write(0,*) 'isTargetFileDefined: ',isTargetFileDefined
    write(0,*) 'isOutputFileDefined: ',isOutputFileDefined
    write(0,*) 'iWantDisp:           ',iWantDisp
    write(0,*) 'iWantStrain:         ',iWantStrain
    write(0,*) 'iWantStress:         ',iWantStress
    write(0,*) 'iWantTraction:       ',iWantTraction
    write(0,*) 'poisson:             ',poisson
    write(0,*) 'lame:                ',lame
    write(0,*) 'mu:                  ',mu
    write(0,*) 'coord_type:          ',coord_type
    write(0,*) 'verbosity:           ',verbosity
endif

9001 if (ios.ne.0) then
    write(0,*) 'gcmdln: error reading string ',trim(tag)
    call usage('exiting triutil from gcmdln')
endif

return
end subroutine gcmdln

!--------------------------------------------------------------------------------------------------!

subroutine usage(str)
implicit none
character(len=*) :: str
if (trim(str).ne.'') then
    write(0,*) trim(str)
    write(0,*)
endif
write(0,'(A)') 'Usage: triutil ...options...'
write(0,'(A)')
write(0,'(A)') '-flt FLTFILE         Triangular source faults'
write(0,'(A)') '-sta STAFILE         Station/receiver locations'
write(0,'(A)') '-trg TRGFILE         Target/receiver kinematics (str dip rak fric)'
write(0,'(A)') '-haf HAFSPCFILE      Elastic half-space properties'
write(0,'(A)') '-poisson POISSON     Poissons ratio (default: 0.25)'
write(0,'(A)') '-lame LAME           Lames parameter (default: 40e9)'
write(0,'(A)') '-mu MU               Shear modulus (default: 40e9)'
write(0,'(A)') '-disp DSPFILE        Displacement (E N Z)'
write(0,'(A)') '-strain STNFILE      Strain matrix (EE NN ZZ EN EZ NZ)'
write(0,'(A)') '-stress STSFILE      Stress matrix (EE NN ZZ EN EZ NZ)'
write(0,'(A)') '-estress ESTSFILE    Effective (maximum) shear stress'
write(0,'(A)') '-normal NORFILE      Normal traction on target faults (requires -trg)'
write(0,'(A)') '-shear[max] SHRFILE  Shear traction on target faults (requires -trg)'
write(0,'(A)') '-coul COULFILE       Coulomb stress on target faults (requires -trg)'
write(0,'(A)') '-geographic (default)'
write(0,'(A)') '-cartesian'
write(0,*)
stop
end subroutine usage
