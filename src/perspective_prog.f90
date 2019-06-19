module perspective_prog

character(len=512) :: input_file
character(len=512) :: output_file

double precision :: viewer_x
double precision :: viewer_y
double precision :: viewer_z
double precision :: look_x
double precision :: look_y
double precision :: look_z

end module

!==================================================================================================!

program main

use io, only: stderr, stdin, stdout, fileExists
use geom, only: perspective

use perspective_prog, only: input_file, &
                            output_file, &
                            viewer_x, &
                            viewer_y, &
                            viewer_z, &
                            look_x, &
                            look_y, &
                            look_z

implicit none

! Local variables
integer :: luin, luout, nline, ios
character(len=512) :: input_line
double precision :: x, y, z, xproj, yproj


! Parse command line
call gcmdln()

! Check that input file exists, if defined
if (input_file.ne.'') then
    if (.not.fileExists(input_file)) then
        call usage('perspective: no input file found named '//trim(input_file))
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

    ! Read input coordinates
    read (input_line,*,iostat=ios) x,y,z
    if (ios.ne.0) then
        write(stderr,*) 'perspective: error reading line ',nline,':'
        call usage(trim(input_line))
    endif

    ! Compute projected points
    call perspective(x,y,z,viewer_x,viewer_y,viewer_z,look_x,look_y,look_z,xproj,yproj)

    ! Write projected points
    write(luout,*) xproj,yproj
enddo


! A little spring cleaning
if (input_file.ne.'') then
    close(luin)
endif
if (output_file.ne.'') then
    close(luout)
endif


end

!--------------------------------------------------------------------------------------------------!

subroutine gcmdln()

use perspective_prog, only: input_file, &
                            output_file, &
                            viewer_x, &
                            viewer_y, &
                            viewer_z, &
                            look_x, &
                            look_y, &
                            look_z

implicit none


! Local variables
integer :: i, j, narg
character(len=512) :: tag


! Initialize control parameters
input_file = ''
output_file = ''
viewer_x = 0.0d0
viewer_y = 0.0d0
viewer_z = 0.0d0
look_x = 0.0d0
look_y = 0.0d0
look_z = 0.0d0


! Number of arguments
narg = command_argument_count()
if (narg.eq.0) then
  call usage('')
endif


! Read arguments
i = 1
do while (i.le.narg)

    call get_command_argument(i,tag)

    if (tag.eq.'-f') then
        i = i + 1
        call get_command_argument(i,input_file)
    elseif (tag.eq.'-o') then
        i = i + 1
        call get_command_argument(i,output_file)
    elseif (tag.eq.'-viewer') then
        i = i + 1
        call get_command_argument(i,tag)
        j = index(tag,',')
        tag(j:j) = ' '
        j = index(tag,',')
        tag(j:j) = ' '
        read(tag,*) viewer_x,viewer_y,viewer_z
    elseif (tag.eq.'-look') then
        i = i + 1
        call get_command_argument(i,tag)
        j = index(tag,',')
        tag(j:j) = ' '
        j = index(tag,',')
        tag(j:j) = ' '
        read(tag,*) look_x,look_y,look_z
    endif

    i = i + 1

enddo

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine usage(str)

use io, only: stderr

implicit none

! Arguments
character(len=*) :: str

if (str.ne.'') then
    write(stderr,*) trim(str)
    write(stderr,*)
endif

write(stderr,*) 'Usage: perspective -viewer vx,vy,vz -look lx,ly,lz [-f IFILE] [-o OFILE]'
write(stderr,*)
write(stderr,*) '-viewer vx,vy,vz       Viewer position'
write(stderr,*) '-look lx,ly,lz         Look vector (from viewer to vanishing point)'
write(stderr,*) '-f IFILE               Input file: x y z'
write(stderr,*) '-o OFILE               Output file: projected x y'

write(stderr,*)

stop
end subroutine
