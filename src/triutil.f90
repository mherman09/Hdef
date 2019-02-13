module triutil_module

character(len=512) :: fault_file       ! source faults (x1 y1 z1 x2 y2 z2 x3 y3 z3 ss ds ts)
character(len=512) :: station_file     ! station locations (stlo stla stdp)
character(len=512) :: target_file      ! target/receiver fault geometry (str dip rak fric)
character(len=512) :: halfspace_file         ! half-space parameters
character(len=512) :: displacement_file      ! displacement
character(len=512) :: strain_file      ! strain tensor
character(len=512) :: stress_file      ! stress tensor
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
double precision, allocatable :: faults(:,:)
integer :: nstations
double precision, allocatable :: stations(:,:)
double precision :: poisson, lame, mu

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

if (.not.isFaultFileDefined) then
    call usage('read_faults: no fault file defined')
endif
inquire(exist=doesFaultFileExist,file=fault_file)
if (.not.doesFaultFileExist) then
    call usage('read_faults: no fault file named '//trim(fault_file)//' found')
endif

call line_count(nfaults,fault_file)
allocate(faults(nfaults,12))

open(unit=12,file=fault_file,status='old')
do i = 1,nfaults
    read(12,'(A)') input_line
    read(input_line,*,iostat=ios,err=1001,end=1002) (faults(i,j),j=1,12)
enddo
close(12)

1001 if (ios.ne.0) then
    write(0,*) 'read_faults: line end'
    call usage('offending line: '//trim(input_line))
endif
1002 if (ios.ne.0) then
    write(0,*) 'read_faults: read error'
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

! Check that stations are defined and the file exists
if (.not.isStationFileDefined) then
    call usage('read_stations: no station file defined')
endif
inquire(exist=doesStationFileExist,file=station_file)
if (.not.doesStationFileExist) then
    call usage('read_stations: no station file named '//trim(station_file)//' found')
endif

! Count number of stations and allocate memory
call line_count(nstations,station_file)
allocate(stations(nstations,3))

! Read stations and convert to Cartesian if geographic
open(unit=13,file=station_file,status='old')
do i = 1,nstations
    read(13,'(A)') input_line
    read(input_line,*,iostat=ios,err=1003,end=1004) (stations(i,j),j=1,3)
enddo
close(13)

1003 if (ios.ne.0) then
    write(0,*) 'read_stations: line end'
    call usage('offending line: '//trim(input_line))
endif
1004 if (ios.ne.0) then
    write(0,*) 'read_stations: read error'
    call usage('offending line: '//trim(input_line))
endif
return
end subroutine read_stations

!--------------------------------------------------------------------------------------------------!

subroutine calc_deformation()
use triutil_module, only: displacement_file, strain_file, stress_file, normal_file, shear_file, &
                          coulomb_file, iWantDisp, iWantStrain, iWantStress, iWantTraction, &
                          nstations, stations
implicit none
integer :: iSta
logical :: isThisUnitOpen
double precision :: disp(3)

if (iWantDisp) then
    write(0,*) 'calc_deformation: calculating displacements'
    open(unit=101,file=displacement_file,status='unknown')
endif
if (iWantStrain) then
    write(0,*) 'calc_deformation: calculating strains'
    if (strain_file.ne.'') then
        open(unit=111,file=strain_file,status='unknown')
    endif
endif
if (iWantStress) then
    write(0,*) 'calc_deformation: calculating stresses'
    if (stress_file.ne.'') then
        open(unit=121,file=stress_file,status='unknown')
    endif
endif
if (iWantTraction) then
    write(0,*) 'calc_deformation: calculating resolved tractions'
    if (normal_file.ne.'') then
        open(unit=122,file=normal_file,status='unknown')
    endif
    if (shear_file.ne.'') then
        open(unit=123,file=shear_file,status='unknown')
    endif
    if (coulomb_file.ne.'') then
        open(unit=124,file=coulomb_file,status='unknown')
    endif
endif

do iSta = 1,nstations
    if (iWantDisp) then
        call calc_displacement(disp,stations(iSta,1),stations(iSta,2),stations(iSta,3))
        write(101,*) stations(iSta,:),disp
    endif
    ! if (iWantStrain) then
    !     call calc_strain(strain,stations(iSta,:))
    ! endif
    ! if (iWantStress) then
    !     call calc_stress(stress,strain)
    ! endif
    ! if (iWantTraction) then
    !     call calc_traction(traction,stress)
    ! endif
enddo

inquire(101,opened=isThisUnitOpen)
if (isThisUnitOpen) then
    close(101)
endif
inquire(111,opened=isThisUnitOpen)
if (isThisUnitOpen) then
    close(111)
endif
inquire(121,opened=isThisUnitOpen)
if (isThisUnitOpen) then
    close(121)
endif
inquire(122,opened=isThisUnitOpen)
if (isThisUnitOpen) then
    close(122)
endif
inquire(123,opened=isThisUnitOpen)
if (isThisUnitOpen) then
    close(123)
endif
inquire(124,opened=isThisUnitOpen)
if (isThisUnitOpen) then
    close(124)
endif

return
end subroutine calc_deformation

!--------------------------------------------------------------------------------------------------!

subroutine calc_displacement(disp,x,y,z)
use tri_disloc_module, only: tri_disloc_disp, tri_center
use triutil_module, only: nfaults, faults, poisson, coord_type
implicit none
! Arguments
double precision :: disp(3), x, y, z, dist, az
! Local variables
double precision :: disptmp(3), sta_coord(3), tri_coord(3,4), tri_coord_new(3,4), slip(3), center(3)
integer :: iFlt, iTri

! Initialize displacement
disp = 0.0d0

! Set coordinates to meters
if (coord_type.eq.'geographic') then
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
    slip(1) = -faults(iFlt,10)        ! input: left-lateral positive
    slip(2) = -faults(iFlt,11)        ! input: thrust positive
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
        ! Convert x(km) y(km) z(km) to x(m) y(m) z(m)
        tri_coord_new = tri_coord*1.0d3
    endif

    write(0,*) 'iFlt',iFlt
    write(0,*) 'sta_coord',sta_coord
    write(0,*) 'tri_coord_new',tri_coord_new(:,1)
    write(0,*) '             ',tri_coord_new(:,2)
    write(0,*) '             ',tri_coord_new(:,3)
    write(0,*) 'slip',slip
    call tri_disloc_disp(disptmp, sta_coord, tri_coord_new, poisson, slip)
    disp = disp + disptmp
enddo

return
end subroutine calc_displacement

!--------------------------------------------------------------------------------------------------!

subroutine calc_strain(strain,sta_coord)
use tri_disloc_module, only: tri_disloc_strain
use triutil_module, only: nfaults, faults, poisson
implicit none
double precision :: strain(3,3), sta_coord(3)
double precision :: straintmp(3,3), tri_coord(3,4), slip(3)
integer :: iFlt

! Initialize strain
strain = 0.0d0

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
    slip(1) = faults(iFlt,10)
    slip(2) = faults(iFlt,11)
    slip(3) = faults(iFlt,12)
    call tri_disloc_strain(straintmp, sta_coord, tri_coord, poisson, slip)
enddo

return
end subroutine calc_strain

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
                          displacement_file, strain_file, stress_file, normal_file, shear_file, &
                          coulomb_file, coord_type, &
             isFaultFileDefined, isStationFileDefined, isTargetFileDefined, isOutputFileDefined, &
             iWantDisp, iWantStrain, iWantStress, iWantTraction, poisson, lame, mu
implicit none
! Local variables
integer :: i, narg
character(len=512) :: tag

fault_file = ''
station_file = ''
target_file = ''
halfspace_file = ''
displacement_file = ''
strain_file = ''
stress_file = ''
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
    elseif (trim(tag).eq.'-geographic') then
        coord_type = 'geographic'
    elseif (trim(tag).eq.'-cartesian') then
        coord_type = 'cartesian'
    elseif (trim(tag).eq.'-poisson') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) poisson
    elseif (trim(tag).eq.'-lame') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) lame
    elseif (trim(tag).eq.'-mu') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) mu
    else
        call usage('!! gcmdln: no option '//trim(tag))
    endif
    i = i + 1
enddo

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
write(0,'(A)') '-normal NORFILE      Normal traction on target faults (requires -trg)'
write(0,'(A)') '-shear[max] SHRFILE  Shear traction on target faults (requires -trg)'
write(0,'(A)') '-coul COULFILE       Coulomb stress on target faults (requires -trg)'
write(0,'(A)') '-geographic (default)'
write(0,'(A)') '-cartesian'
write(0,*)
stop
end subroutine usage
