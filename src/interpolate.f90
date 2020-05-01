module interpolate_prog

character(len=512) :: anchor_point_file
character(len=512) :: interpolation_point_file
character(len=512) :: output_file
double precision :: slope1
double precision :: slope2

end module

!==================================================================================================!

program main

use interpolate_prog, only: anchor_point_file, &
                            interpolation_point_file, &
                            output_file, &
                            slope1, &
                            slope2
use interpolation, only: spline, spline_interpolation

use io, only: fileExists, line_count

implicit none

! Local variables
double precision, allocatable :: x(:), y(:), dy(:), d2y(:)
double precision :: xi, yi, dyi, d2yi
integer :: i, nanchor, ninterpolation, lu_out


! Parse command line arguments
call gcmdln()


! Check for inputs
if (.not.fileExists(anchor_point_file)) then
    call usage('interpolate: no anchor point file found named "'//trim(anchor_point_file)//'"')
endif
if (.not.fileExists(interpolation_point_file)) then
    call usage('interpolate: no interpolation point file found named "'// &
                   trim(interpolation_point_file)//'"')
endif


! Count number of inputs and allocate memory
nanchor = line_count(anchor_point_file)
ninterpolation = line_count(interpolation_point_file)

allocate(x(nanchor))
allocate(y(nanchor))
allocate(dy(nanchor))
allocate(d2y(nanchor))


! Read the anchor points
open(unit=11,file=anchor_point_file,status='old')
do i = 1,nanchor
    read(11,*) x(i),y(i)
enddo
close(11)


! Calculate the interpolation
call spline(x,y,nanchor,slope1,slope2,d2y)


! Select the points to calculate the interpolation
open(unit=12,file=interpolation_point_file,status='old')
if (output_file.ne.'') then
    lu_out = 13
else
    lu_out = 6
endif
do i = 1,ninterpolation
    read(12,*) xi
    call spline_interpolation(x,y,d2y,nanchor,xi,yi,dyi,d2yi)
    write(lu_out,*) xi,yi,dyi,d2yi
enddo
close(12)

if (lu_out.ne.6) then
    close(lu_out)
endif

end

!--------------------------------------------------------------------------------------------------!

subroutine gcmdln()

use interpolate_prog, only: anchor_point_file, &
                            interpolation_point_file, &
                            output_file, &
                            slope1, &
                            slope2

implicit none

! Local variables
integer :: i, narg
character(len=512) :: tag

! Initialize control parameters
anchor_point_file = ''
interpolation_point_file = ''
output_file = ''
slope1 = 1.01d30
slope2 = 1.01d30

! Number of arguments
narg = command_argument_count()
if (narg.eq.0) then
  call usage('')
endif

! Read arguments
i = 1
do while (i.le.narg)

    call get_command_argument(i,tag)

    if (tag.eq.'-a') then
        i = i + 1
        call get_command_argument(i,anchor_point_file)
    elseif (tag.eq.'-i') then
        i = i + 1
        call get_command_argument(i,interpolation_point_file)
    elseif (tag.eq.'-o') then
        i = i + 1
        call get_command_argument(i,output_file)
    elseif (tag.eq.'-slope1') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) slope1
    elseif (tag.eq.'-slope2') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) slope2
    else
        call usage('interpolate: no option '//trim(tag))
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

write(stderr,*) 'Usage: interpolate -a ANCHOR_POINTS -i INTERPOLATION_POINTS [-o OUTPUT_FILE]'
write(stderr,*) '       [-slope[1|2] SLOPE]'
write(stderr,*)
write(stderr,*) '-a ANCHOR_POINTS         Points in function (x,y)'
write(stderr,*) '-i INTERPOLATION_POINTS  X-coordinates to interpolate (xi)'
write(stderr,*) '-o OUTPUT_FILE           Output file: (xi,yi) (default: stdout)'
write(stderr,*) '-slope[1|2] SLOPE        Slope at first or last point'
write(stderr,*)

call error_exit(1)
end subroutine
