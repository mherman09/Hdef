module sphfinrot

    character(len=512) :: input_file
    character(len=512) :: output_file

    double precision :: pole_lon
    double precision :: pole_lat
    double precision :: rotation_angle

end module

!==================================================================================================!

program main
!----
! Given a point (ilon,ilat), a pole of rotation (plon,plat), and a finite angle of rotation
! (degrees), calculate the rotated point (olon,olat).
!
! Based on equations shown by Glenn Murray:
! https://sites.google.com/site/glennmurray/Home/rotation-matrices-and-formulas/rotation-about-an-arbitrary-axis-in-3-dimensions
!----

use io, only: stdin, stderr, stdout, fileExists

use sphfinrot, only: input_file, &
                     output_file, &
                     pole_lon, &
                     pole_lat, &
                     rotation_angle

implicit none

! Local variables
integer :: luin, luout, ios
character(len=512) :: line
double precision :: input_lon, input_lat, output_lon, output_lat


call gcmdln()

! Input file
if (input_file.eq.'none'.or.input_file.eq.'stdin') then
    luin = stdin
else
    if (.not.fileExists(input_file)) then
        call usage('sphfinrot: no input file found named '//trim(input_file))
    endif
    luin = 11
    open(luin,file=input_file,status='old')
endif

! Output file
if (output_file.eq.'none'.or.output_file.eq.'stdout') then
    luout = stdout
else
    luout = 12
    open(luout,file=output_file,status='unknown')
endif

! Pole of finite rotation
if (pole_lon.lt.-1.0d9.or.pole_lat.lt.-1.0d9) then
    call usage('sphfinrot: coordinates of rotation pole not defined')
elseif (pole_lon.lt.-180.0d0) then
    call usage('sphfinrot: rotation pole longitude less than -180 degrees')
elseif (pole_lon.gt.360.0d0) then
    call usage('sphfinrot: rotation pole longitude greater than 360 degrees')
elseif (pole_lat.lt.-90.0d0) then
    call usage('sphfinrot: rotation pole latitude less than -90 degrees')
elseif (pole_lat.gt.90.0d0) then
    call usage('sphfinrot: rotation pole latitude greater than 90 degrees')
endif

! Angle
if (rotation_angle.lt.-1.0d9) then
    call usage('sphfinrot: rotation angle not defined')
endif

do
    read(luin,'(A)',iostat=ios) line
    if (ios.ne.0) then
        exit
    endif
    read(line,*,iostat=ios) input_lon,input_lat
    if (ios.ne.0) then
        write(stderr,*) 'sphfinrot: error parsing input line'
        write(stderr,*) 'Input line: ',trim(line)
        call usage('Expected format: lon lat')
    endif
    call finite_rotate(pole_lon,pole_lat,rotation_angle,input_lon,input_lat,output_lon,output_lat)
    write(luout,1001) output_lon,output_lat
enddo
1001 format(1P2E14.6)

if (input_file.ne.'none'.and.input_file.ne.'stdin') then
    close(luin)
endif
if (output_file.ne.'none'.and.output_file.ne.'stdout') then
    close(luout)
endif

end

!--------------------------------------------------------------------------------------------------!

subroutine finite_rotate(pole_lon,pole_lat,angle,lon_in,lat_in,lon_out,lat_out)

use trig, only: d2r
use earth, only: radius_earth_m, sphere_geo2xyz, sphere_xyz2geo

implicit none

! Arguments
double precision :: pole_lon, pole_lat, angle, lon_in, lat_in, lon_out, lat_out


! Local variables
double precision :: x, y, z, px, py, pz, cosr, sinr, dot, ox, oy, oz


! Convert input point to Cartesian coordinates
call sphere_geo2xyz(lon_in,lat_in,0.0d0,x,y,z,radius_earth_m)

! Convert finite rotation pole to Cartesian coordinates (unit vector)
call sphere_geo2xyz(pole_lon,pole_lat,0.0d0,px,py,pz,1.0d0)

! Compute the rotation with a matrix multiplication
cosr = dcos(angle*d2r)
sinr = dsin(angle*d2r)
dot = px*x + py*y + pz*z
ox = px*dot*(1.0d0-cosr) + x*cosr + (py*z-pz*y)*sinr
oy = py*dot*(1.0d0-cosr) + y*cosr + (pz*x-px*z)*sinr
oz = pz*dot*(1.0d0-cosr) + z*cosr + (px*y-py*x)*sinr

! Convert back to geographic coordinates
call sphere_xyz2geo(ox,oy,oz,lon_out,lat_out,dot,radius_earth_m)

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine gcmdln()

use io, only: stdout, verbosity

use sphfinrot, only: input_file, &
                     output_file, &
                     pole_lon, &
                     pole_lat, &
                     rotation_angle

implicit none

! Local variables
character(len=512) :: tag
integer :: i, narg

input_file = 'none'
output_file = 'none'
pole_lon = -1.0d10
pole_lat = -1.0d10
rotation_angle = -1.0d10

narg = command_argument_count()
if (narg.eq.0) then
    call usage('')
endif

i = 1
do while (i.le.narg)

    call get_command_argument(i,tag)

    if (trim(tag).eq.'-f') then
        i = i + 1
        call get_command_argument(i,input_file)

    elseif (trim(tag).eq.'-o') then
        i = i + 1
        call get_command_argument(i,output_file)

    elseif (trim(tag).eq.'-pole') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) pole_lon
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) pole_lat

    elseif (trim(tag).eq.'-angle') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) rotation_angle

    elseif (trim(tag).eq.'-v') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) verbosity

    else
        call usage('sphfinrot: no option '//trim(tag))
    endif

    i = i + 1
enddo

if (verbosity.ge.2) then
    write(stdout,*) 'input_file:     ',trim(input_file)
    write(stdout,*) 'output_file:    ',trim(output_file)
    write(stdout,*) 'pole_lon:       ',pole_lon
    write(stdout,*) 'pole_lon:       ',pole_lat
    write(stdout,*) 'rotation_angle: ',rotation_angle
endif


return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine usage(str)

use io, only: stderr
implicit none

character(len=*) :: str

if (str.ne.'') then
    write(stderr,*) trim(str)
    write(stderr,*)
endif
write(stderr,*) 'Usage: sphfinrot [-f IFILE] [-o OFILE] -pole PLON PLAT -angle ANGLE'
write(stderr,*)
write(stderr,*) '  -f IFILE         Input longitude and latitude (default: stdin)'
write(stderr,*) '  -o OFILE         Output longitude and latitude (default: stdout)'
write(stderr,*) '  -pole PLON PLAT  Rotation pole location'
write(stderr,*) '  -angle ANGLE     Rotation angle (counter-clockwise)'
write(stderr,*)

stop
end subroutine
