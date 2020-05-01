module distaz2lola_prog

character(len=512) :: input_file
character(len=512) :: output_file
character(len=32) :: lon_range

end module

!==================================================================================================!

program main
!----
! Compute the longitude and latitude of points, given a starting point and distance/azimuth from
! the starting point.
!
! (Public implementation of subroutine distaz2lola() in geom_module.f90)
!----

use io, only: stderr, stdin, stdout, fileExists
use geom, only: distaz2lola
use earth, only: radius_earth_km

use distaz2lola_prog, only: input_file, &
                            output_file, &
                            lon_range

implicit none

! Local variables
integer :: luin, luout, ios, nline
double precision :: lon1, lat1, dist, az, lon2, lat2
character(len=512) :: input_line


call gcmdln()

! Check that input file exists, if defined
if (input_file.ne.'') then
    if (.not.fileExists(input_file)) then
        call usage('distaz2lola: no input file found named '//trim(input_file))
    endif
endif

! Open input and output streams
if (input_file.eq.''.or.input_file.eq.'stdin') then
    luin = stdin
else
    luin = 11
    open(unit=luin,file=input_file,status='old')
endif

if (output_file.eq.''.or.output_file.eq.'stdout') then
    luout = stdout
else
    luout = 12
    open(unit=luout,file=output_file,status='unknown')
endif


! Calculate the longitude and latitude
nline = 0
do
    ! Read input line
    read (luin,'(A)',iostat=ios) input_line
    if (ios.ne.0) then
        exit
    else
        nline = nline + 1
    endif

    ! Read lon, lat, dist, az from input line
    read (input_line,*,iostat=ios) lon1,lat1,dist,az
    if (ios.ne.0) then
        write(stderr,*) 'distaz2lola: error reading line ',nline,':'
        call usage(trim(input_line))
    endif

    ! Make sure latitudes are reasonable
    if (lat1.gt.90.0d0) then
        write(stderr,*) 'distaz2lola: latitude in line ',nline,' is greater than 90 degrees:'
        call usage(trim(input_line))
    elseif (lat1.lt.-90.0d0) then
        write(stderr,*) 'distaz2lola: latitude in line ',nline,' is less than -90 degrees:'
        call usage(trim(input_line))
    endif

    ! Compute longitude and latitude
    call distaz2lola(lon1,lat1,dist/radius_earth_km,az,lon2,lat2,'radians','degrees',ios)
    if (ios.ne.0) then
        call usage('distaz2lola: error computing longitude and latitude')
    endif

    ! Make sure lon is in specified range
    if (lon_range.eq.'NEG'.and.lon2.gt.180.0d0) then
        lon2 = lon2-360.0d0
    elseif (lon_range.eq.'POS'.and.lon2.lt.0.0d0) then
        lon2 = lon2+360.0d0
    endif

    ! Write lon/lat
    write(luout,*) lon2,lat2
enddo


! A little spring cleaning
if (input_file.ne.'') then
    close(luin)
endif
if (output_file.ne.'') then
    close(luout)
endif
if (input_file.eq.'distaz2lola_86_this_at_the_end.tmp') then
    open(unit=86,file=input_file,status='old')
    close(86,status='delete')
endif

end

!--------------------------------------------------------------------------------------------------!

subroutine gcmdln()

use distaz2lola_prog, only: input_file, &
                           output_file, &
                           lon_range

implicit none

! Local variables
integer :: i, narg
character(len=512) :: tag, tag_array(4)

! Initialize control parameters
input_file = ''
output_file = ''
lon_range = 'POS'

! Number of arguments
narg = command_argument_count()
if (narg.eq.0) then
  call usage('')
endif

! Read arguments
i = 1
do while (i.le.narg)

    call get_command_argument(i,tag)

    if (trim(tag).eq.'-f') then
        i = i + 1
        call get_command_argument(i,input_file)
        if (input_file.eq.'stdin') then
            input_file = ''
        endif
    elseif (trim(tag).eq.'-s') then
        ! Do nothing: default reads from stdin
    elseif (trim(tag).eq.'-o') then
        i = i + 1
        call get_command_argument(i,output_file)
    elseif (trim(tag).eq.'-c') then
        i = i + 1
        call get_command_argument(i,tag_array(1))
        i = i + 1
        call get_command_argument(i,tag_array(2))
        i = i + 1
        call get_command_argument(i,tag_array(3))
        i = i + 1
        call get_command_argument(i,tag_array(4))
        input_file = 'distaz2lola_86_this_at_the_end.tmp'
        open(unit=86,file=input_file,status='unknown')
        write(86,*) trim(tag_array(1)),' ', &
                    trim(tag_array(2)),' ', &
                    trim(tag_array(3)),' ', &
                    trim(tag_array(4))
        close(86)
    elseif (trim(tag).eq.'-lon') then
        i = i + 1
        call get_command_argument(i,lon_range)
    else
        call usage('distaz2lola: no option '//trim(tag))
    endif

    i = i + 1

enddo

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

write(stderr,*) 'Usage: distaz2lola -s|-f IFILE|-c LON1 LAT1 DIST AZ [-o OFILE] [-lon NEG|POS]'
write(stderr,*)
write(stderr,*) '-s                     Read standard input: lon lat dist(km) az(deg CW from N)'
write(stderr,*) '-f IFILE               Input file: lon lat dist(km) az(deg CW from N)'
write(stderr,*) '-o OFILE               Output file: lon lat'
write(stderr,*) '-c LON1 LAT1 DIST AZ   Compute LON2 LAT2 for a single DIST AZ'
write(stderr,*) '-lon NEG|POS           Output longitude range: -180,180 or 0,360(default)'
write(stderr,*)

call error_exit(1)
end subroutine
