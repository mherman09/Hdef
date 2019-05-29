module stereo_project

character(len=32) :: project_mode
character(len=512) :: input_file
character(len=512) :: output_file
double precision :: lon0
double precision :: lat0
double precision :: radius_sphere

end module

!==================================================================================================!

program main

use io, only: stderr, fileExists
use map, only: stereo

use stereo_project, only: project_mode, &
                          input_file, &
                          output_file, &
                          lon0, &
                          lat0, &
                          radius_sphere

implicit none

! Local variables
character(len=128) :: line
integer :: ios, ierr, luin, luout
double precision :: lon, lat, x, y


! Parse command line
call gcmdln()
!write(0,*) trim(input_file)
!write(0,*) trim(output_file)
!write(0,*) trim(mode)
!write(0,*) lon0
!write(0,*) lat0
!write(0,*) radius


if (input_file.eq.'stdin') then
    luin = 5
else
    if (.not.fileExists(input_file)) then
        call usage('stereo_project: no input file found named '//trim(input_file))
    endif
    luin = 11
    open(unit=11,file=input_file,status='old')
endif

if (trim(output_file).eq.'stdout') then
    luout = 6
else
    luout = 12
    open(unit=12,file=output_file,status='unknown')
endif

! Stereographic projection or inverse
do
    read(luin,'(A)',iostat=ios) line
    if (ios.ne.0) then
        exit
    endif

    if (project_mode.eq.'geo2stereo') then
        read(line,*,iostat=ios) lon, lat
        if (ios.ne.0) then
            write(stderr,*) 'stereo_project: error parsing lon lat from line'
            call usage(trim(line))
        endif

        call stereo(lon,lat,x,y,lon0,lat0,radius_sphere,'lonlat2xy',ierr)
        if (ierr.ne.0) then
            call usage('stereo_project: error computing stereographic projection')
        endif

        write(luout,*) x, y

    elseif (project_mode.eq.'stereo2geo') then
        read(line,*,iostat=ios) x, y
        if (ios.ne.0) then
            write(stderr,*) 'stereo_project: error parsing x y from line'
            call usage(trim(line))
        endif

        call stereo(lon,lat,x,y,lon0,lat0,radius_sphere,'xy2lonlat',ierr)
        if (ierr.ne.0) then
            call usage('stereo_project: error computing inverse stereographic projection')
        endif

        write(luout,*) lon, lat

    else
        call usage('stereo_project: project_mode must be "geo2stereo" or "stereo2geo"')
    endif
enddo

end

!--------------------------------------------------------------------------------------------------!

subroutine gcmdln()

use earth, only: radius_earth_km

use stereo_project, only: project_mode, &
                          input_file, &
                          output_file, &
                          lon0, &
                          lat0, &
                          radius_sphere

implicit none

! Local variables
integer :: i, narg
character(len=1024) :: tag


project_mode = ''
input_file = 'stdin'
output_file = 'stdout'
lon0 = 1.0d3
lat0 = 1.0d3
radius_sphere = radius_earth_km

narg = command_argument_count()
if (narg.eq.0) then
    call usage('')
endif

i = 1

do while (i.le.narg)

    call get_command_argument(i,tag)

    if (trim(tag).eq.'-geo2stereo'.or.tag.eq.'-geo2xy'.or.tag.eq.'-lonlat2xy') then
        project_mode = 'geo2stereo'
    elseif (trim(tag).eq.'-stereo2geo'.or.tag.eq.'-xy2geo'.or.tag.eq.'-xy2lonlat') then
        project_mode = 'stereo2geo'

    elseif (trim(tag).eq.'-f') then
        i = i + 1
        call get_command_argument(i,input_file)
    elseif (trim(tag).eq.'-o') then
        i = i + 1
        call get_command_argument(i,output_file)

    elseif (trim(tag).eq.'-lon0') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) lon0
    elseif (trim(tag).eq.'-lat0') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) lat0
    elseif (trim(tag).eq.'-radius') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) radius_sphere

    else
        call usage('stereo_project: no option '//trim(tag))
    endif

    i = i + 1
enddo

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine usage(string)

use io, only: stderr

implicit none

character(len=*) :: string

if (string.ne.'') then
    write(stderr,*) trim(string)
    write(stderr,*)
endif

write(stderr,*) 'Usage: stereo_project ...options...'
write(stderr,*)
write(stderr,*) '-f FILE         Input geographic or projected x-y coordinates (default: stdin)'
write(stderr,*) '-o FILE         Output coordinates (default: stdout)'
write(stderr,*) '-geo2stereo     Project geographic to x-y coordinates'
write(stderr,*) '-stereo2geo     Project x-y coordinates to geographic'
write(stderr,*) '-lon0 LON0      Center longitude of projection'
write(stderr,*) '-lat0 LAT0      Center latitude of projection'
write(stderr,*) '-radius R       Sphere radius'

stop
end subroutine
