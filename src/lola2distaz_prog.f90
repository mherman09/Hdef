module lola2distaz_prog

character(len=512) :: input_file
character(len=512) :: output_file

end module

!==================================================================================================!

program main
!----
! Compute the distance and azimuth between two geographic points.
!
! (Public implementation of subroutine lola2distaz() in geom_module.f90)
!----

use io, only: stderr, stdin, stdout, fileExists
use geom, only: lola2distaz
use earth, only: radius_earth_km

use lola2distaz_prog, only: input_file, &
                            output_file
implicit none

! Local variables
integer :: luin, luout, ios, nline
double precision :: lon1, lat1, dist, az, lon2, lat2
character(len=512) :: input_line


call gcmdln()

! Check that input file exists, if defined
if (input_file.ne.'') then
    if (.not.fileExists(input_file)) then
        call usage('lola2distaz: no input file found named '//trim(input_file))
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


! Calculate the distance and azimuth
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
    read (input_line,*,iostat=ios) lon1,lat1,lon2,lat2
    if (ios.ne.0) then
        write(stderr,*) 'lola2distaz: error reading line ',nline,':'
        call usage(trim(input_line))
    endif

    ! Make sure latitudes are reasonable
    if (lat1.gt.90.0d0) then
        write(stderr,*) 'lola2distaz: latitude 1in line ',nline,' is greater than 90 degrees:'
        call usage(trim(input_line))
    elseif (lat1.lt.-90.0d0) then
        write(stderr,*) 'lola2distaz: latitude 1 in line ',nline,' is less than -90 degrees:'
        call usage(trim(input_line))
    elseif (lat2.gt.90.0d0) then
        write(stderr,*) 'lola2distaz: latitude 2 in line ',nline,' is greater than 90 degrees:'
        call usage(trim(input_line))
    elseif (lat2.lt.-90.0d0) then
        write(stderr,*) 'lola2distaz: latitude 2 in line ',nline,' is less than -90 degrees:'
        call usage(trim(input_line))
    endif

    ! Compute longitude and latitude
    call lola2distaz(lon1,lat1,lon2,lat2,dist,az)
    dist = dist*radius_earth_km

    ! Write distance/azimuth
    write(luout,*) dist,az
enddo


! A little spring cleaning
if (input_file.ne.'') then
    close(luin)
endif
if (output_file.ne.'') then
    close(luout)
endif
if (input_file.eq.'lola2distaz_86_this_at_the_end.tmp') then
    open(unit=86,file=input_file,status='old')
    close(86,status='delete')
endif

end

!--------------------------------------------------------------------------------------------------!

subroutine gcmdln()

use lola2distaz_prog, only: input_file, &
                            output_file

implicit none

! Local variables
integer :: i, narg
character(len=512) :: tag, tag_array(4)

! Initialize control parameters
input_file = ''
output_file = ''

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
        input_file = 'lola2distaz_86_this_at_the_end.tmp'
        open(unit=86,file=input_file,status='unknown')
        write(86,*) trim(tag_array(1)),' ', &
                    trim(tag_array(2)),' ', &
                    trim(tag_array(3)),' ', &
                    trim(tag_array(4))
        close(86)
    else
        call usage('gcmdln: no option '//trim(tag))
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

write(stderr,*) 'Usage: lola2distaz -s|-f IFILE [-o OFILE] [-c LON1 LAT1 DIST AZ]'
write(stderr,*)
write(stderr,*) '-s                     Read standard input: lon1 lat1 lon2 lat2'
write(stderr,*) '-f IFILE               Input file: lon1 lat1 lon2 lat2'
write(stderr,*) '-o OFILE               Output file: dist(km) az(deg CW from N)'
write(stderr,*) '-c LON1 LAT1 LON2 LAT2 Compute DIST AZ for a single LON LAT pair'
write(stderr,*)

stop
end subroutine
