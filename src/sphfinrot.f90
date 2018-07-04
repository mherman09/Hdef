module sphfinrot_cmdln
    ! Command line variables
    character(len=1024) :: input_file
    character(len=1024) :: output_file
endmodule

module sphfinrot_vars
    ! Program variables
    double precision :: input_lon
    double precision :: input_lat
    double precision :: output_lon
    double precision :: output_lat
    double precision :: pole_lon
    double precision :: pole_lat
    double precision :: rotation_angle
endmodule

module sphfinrot_const
    ! Constants
    double precision, parameter :: pi = 4.0d0*datan(1.0d0)
    double precision, parameter :: d2r = pi/1.8d2
    double precision, parameter :: r2d = 1.8d2/pi
endmodule

!--------------------------------------------------------------------------------------------------!

program sphfinrot
!----
! Given a point (ilon,ilat), a pole of rotation (plon,plat), and a finite angle of rotation
! (degrees), calculate the rotated point (olon,olat).
!
! Based on equations shown by Glenn Murray:
! https://sites.google.com/site/glennmurray/Home/rotation-matrices-and-formulas/rotation-about-an-arbitrary-axis-in-3-dimensions
!----
use sphfinrot_cmdln
use sphfinrot_vars
implicit none
! Local variables
integer :: luin, luout, ios

call gcmdln()
write(0,*) 'input_file:     ',trim(input_file)
write(0,*) 'output_file:    ',trim(output_file)
write(0,*) 'pole_lon:       ',pole_lon
write(0,*) 'pole_lon:       ',pole_lat
write(0,*) 'rotation_angle: ',rotation_angle

call check_inputs(luin,luout)

do
    read(luin,*,iostat=ios) input_lon,input_lat
    if (ios.ne.0) then
        exit
    endif
    call rotate()
    write(luout,1001) output_lon,output_lat
enddo
1001 format(1P2E14.6)

if (input_file.ne.'stdin') then
    close(luin)
endif
if (output_file.ne.'none') then
    close(luout)
endif


      END

!--------------------------------------------------------------------------------------------------!

subroutine check_inputs(luin,luout)
use sphfinrot_cmdln
use sphfinrot_vars
implicit none
! I/O variables
integer :: luin, luout
! Local variables
logical :: ex

! Input file
if (input_file.eq.'stdin') then
    luin = 5
else
    inquire(file=input_file,exist=ex)
    if (.not.ex) then
        call usage('!! Error: no input file found named '//trim(input_file))
    endif
    luin = 11
    open(luin,file=input_file,status='old')
endif

! Output file
if (output_file.eq.'none') then
    luout = 6
else
    luout = 12
    open(luout,file=output_file,status='unknown')
endif

! Pole
if (pole_lon.lt.-1.0d9.or.pole_lat.lt.-1.0d9) then
    call usage('!! Error: location of rotation pole not defined')
endif

! Angle
if (rotation_angle.lt.-1.0d9) then
    call usage('!! Error: rotation angle not defined')
endif

return
end

!--------------------------------------------------------------------------------------------------!

subroutine rotate()
use sphfinrot_vars
use sphfinrot_const
implicit none
! Local variables
double precision :: input_x, input_y, input_z
double precision :: pole_x, pole_y, pole_z
double precision :: output_x, output_y, output_z
double precision :: cosr, sinr, dot_prod

! Convert to Cartesian coordinates
call lonlat2xyz(input_x,input_y,input_z,input_lon,input_lat)
call lonlat2xyz(pole_x,pole_y,pole_z,pole_lon,pole_lat)

! Store some variables
cosr = dcos(rotation_angle*d2r)
sinr = dsin(rotation_angle*d2r)
dot_prod = pole_x*input_x + pole_y*input_y + pole_z*input_z

! Matrix multiplication
output_x = pole_x*(dot_prod)*(1.0d0-cosr) + input_x*cosr + (-pole_z*input_y+pole_y*input_z)*sinr
output_y = pole_y*(dot_prod)*(1.0d0-cosr) + input_y*cosr + ( pole_z*input_x-pole_x*input_z)*sinr
output_z = pole_z*(dot_prod)*(1.0d0-cosr) + input_z*cosr + (-pole_y*input_x+pole_x*input_y)*sinr

! Convert back to geographic coordinates
output_lon = datan2(output_y,output_x)*r2d
output_lat = dasin(output_z)*r2d

return
end

!--------------------------------------------------------------------------------------------------!

subroutine lonlat2xyz(x,y,z,lon,lat)
use sphfinrot_const
implicit none
! I/O variables
double precision :: lon, lat ,x, y, z

x = dcos(lon*d2r)*dcos(lat*d2r)
y = dsin(lon*d2r)*dcos(lat*d2r)
z = dsin(lat*d2r)

return
end

!--------------------------------------------------------------------------------------------------!

subroutine gcmdln()
use sphfinrot_cmdln
use sphfinrot_vars, only : pole_lon, pole_lat, rotation_angle
implicit none
! Local variables
character(len=1024) :: tag
integer :: i, narg

input_file = 'none'
output_file = 'none'
pole_lon = -1.0d10
pole_lat = -1.0d10
rotation_angle = -1.0d10

narg = command_argument_count()
if (narg.eq.0) call usage('')
i = 0
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
    elseif (trim(tag).eq.'-h'.or.trim(tag).eq.'-?') then
        call usage('')
    endif
    i = i + 1
enddo
return
end

!--------------------------------------------------------------------------------------------------!

subroutine usage(string)
!----
! Print program usage statement and exit
!----
implicit none
! I/O variables
character(len=*) :: string
! Local variables
integer :: string_length

if (string.ne.'') then
    string_length = len(string)
    write(0,*) trim(string)
    write(0,*)
endif
write(0,'(A)') 'Usage: sphfinrot -f IFILE [-o OFILE] -pole PLON PLAT -angle ANGLE [-h/-?]'
write(0,'(A)')
write(0,'(A)') '  -f [IFILE]       Input file: input_lon input_lat'
write(0,'(A)') '  -o [OFILE]       Output file (default: stdout): rotated_lon rotated_lat'
write(0,'(A)') '  -pole PLON PLAT  Rotation pole location'
write(0,'(A)') '  -angle ANGLE     Rotation angle (counter-clockwise)'
write(0,'(A)') '  -h/-?            Online help (this message)'

stop
end

