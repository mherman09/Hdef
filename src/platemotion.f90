module platemotion

character(len=512) :: input_file
character(len=512) :: output_file
character(len=16) :: plates(2)
double precision :: pole(3)
character(len=32) :: model_name

end module

!==================================================================================================!

program main
!----
! Compute relative plate velocities given (a) plates from a hard-coded list or (b) an Euler pole
!----

use io, only: stdout, stdin
use trig, only: d2r
use algebra, only: cross_product, normalize
use earth, only: get_pole, pole_geo2xyz, deg2km, semimajor_grs80, eccentricity_grs80, ellipsoid_geo2xyz

use platemotion, only: input_file, &
                       output_file, &
                       plates, &
                       pole, &
                       model_name

implicit none

integer :: ierr, ios, luin, luout
double precision :: xyz_pole(3)
double precision :: lon, lat, lonr, latr, r(3), vel(3), vn, ve, vz
character(len=512) :: input_line


! Parse command line
call gcmdln()


! Calculate pole of rotation between plate pair
if (plates(1).ne.'') then
    call get_pole(plates(1),plates(2),model_name,pole,ierr)
    if (ierr.ne.0) then
        call usage('platemotion: problem computing pole of rotation from plate circuit')
    endif
endif

! Convert pole to Cartesian coordinates
call pole_geo2xyz(pole(1),pole(2),pole(3),xyz_pole(1),xyz_pole(2),xyz_pole(3),'sphere')


! Open input and output streams
if (input_file.eq.'') then
    luin = stdin
else
    luin = 11
    open(unit=luin,file=input_file,status='old')
endif

if (output_file.eq.'') then
    luout = stdout
else
    luout = 12
    open(unit=luout,file=output_file,status='unknown')
endif


! Calculate predicted velocities
do
    read(luin,'(A)',iostat=ios) input_line
    if (ios.ne.0) then
        exit
    endif

    ! Read longitude and latitude
    read(input_line,*) lon, lat

    ! Calculate vector from center of Earth to point
    lonr = lon*d2r
    latr = lat*d2r
    ! call ellipsoid_geo2xyz(lon,lat,0.0d0,r(1),r(2),r(3),semimajor_grs80,eccentricity_grs80)
    ! call normalize(r)
    r(1) = dcos(latr)*dcos(lonr)
    r(2) = dcos(latr)*dsin(lonr)
    r(3) = dsin(latr)

    ! Velocity at that point is cross product of angular velocity vector and position vector
    call cross_product(xyz_pole,r,vel)
    ve = -vel(1)           *dsin(lonr) + vel(2)           *dcos(lonr)
    vn = -vel(1)*dsin(latr)*dcos(lonr) - vel(2)*dsin(latr)*dsin(lonr) + vel(3)*dcos(latr)
    vz =  vel(1)*dcos(latr)*dcos(lonr) + vel(2)*dcos(latr)*dsin(lonr) + vel(3)*dsin(latr)

    ! Units are deg/Ma; convert to km/Ma = mm/yr
    ve = ve*deg2km
    vn = vn*deg2km

    ! Write location and velocity
    write(luout,*) lon, lat, ve, vn
enddo

end

!--------------------------------------------------------------------------------------------------!

subroutine gcmdln()

use earth, only: list_plate_models

use platemotion, only: input_file, &
                       output_file, &
                       plates, &
                       pole, &
                       model_name

implicit none

! Local variables
integer :: i, j, narg
character(len=512) :: tag

! Initialize control parameters
input_file = ''
output_file = ''
plates = ''
pole = 0.0d0
model_name = 'MORVEL56'

! Number of arguments
narg = command_argument_count()
if (narg.eq.0) then
  call usage('')
endif

! Read arguments
i = 1
do while (i.le.narg)

    call get_command_argument(i,tag)


    if (trim(tag).eq.'-plate'.or.trim(tag).eq.'-plates') then
        i = i + 1
        call get_command_argument(i,plates(1))
        j = index(plates(1),'/')
        plates(1)(j:j) = ' '
        read(plates(1),*) plates(2),plates(2)
        plates(1)(j:len(plates(1))) = ''

    elseif (trim(tag).eq.'-pole') then
        i = i + 1
        call get_command_argument(i,tag)
        j = index(tag,'/')
        tag(j:j) = ' '
        j = index(tag,'/')
        tag(j:j) = ' '
        read(tag,*) pole(1),pole(2),pole(3)

    elseif (trim(tag).eq.'-model') then
        i = i + 1
        call get_command_argument(i,model_name)

    elseif (trim(tag).eq.'-list') then
        i = i + 1
        if (i.gt.narg) then
            call list_plate_models('ALL')
        else
            call get_command_argument(i,tag)
            call list_plate_models(tag)
        endif
        stop

    elseif (trim(tag).eq.'-f') then
        i = i + 1
        call get_command_argument(i,input_file)

    elseif (trim(tag).eq.'-o') then
        i = i + 1
        call get_command_argument(i,output_file)

    else
        call usage('platemotion: no command line option '//trim(tag))
    endif

    i = i + 1
enddo

! write(0,*) trim(plates(1)),' ',trim(plates(2))
! write(0,*) pole
! write(0,*) model_name

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine usage(str)
use io, only: stderr
implicit none
character(len=*) :: str

if (trim(str).ne.'') then
    write(stderr,*) trim(str)
    write(stderr,*)
endif

write(stderr,*) 'Usage: platemotion -plate|-pole OPT [-model NAME] [-list [MODEL]] [-f IFILE] ',&
                '[-o OFILE]'
write(stderr,*)
write(stderr,*) '-plate P1/P2       Plates (P2 moving w.r.t. fixed P1)'
write(stderr,*) '-pole LON/LAT/VEL  Euler pole (velocity in deg/Ma)'
write(stderr,*) '-model MODEL_NAME  Plate circuit model (default: MORVEL56)'
write(stderr,*) '                       ITRF2005 (Altamimi et al., 2007)'
write(stderr,*) '                       ITRF2008 (Altamimi et al., 2012)'
write(stderr,*) '                       MORVEL (DeMets et al., 2010)'
write(stderr,*) '                       MORVEL56 (Argus et al., 2011)'
write(stderr,*) '                       NUVEL1A  (DeMets et al., 1994)'
write(stderr,*) '-list MODEL        List plates in MODEL (use ALL or leave blank to see all models)'
write(stderr,*) '-f IFILE           Input file (default: stdin)'
write(stderr,*) '-o OFILE           Output file (default: stdout)'
write(stderr,*)

stop
end subroutine
