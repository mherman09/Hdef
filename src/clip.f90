module clip

character(len=512) :: poly_file
character(len=512) :: input_file
character(len=512) :: output_file
logical :: getPointsInsideBoundary
logical :: getPointsOutsideBoundary
logical :: getPointsOnBoundary
double precision :: epsilon
logical :: printRecord

end module clip

!==================================================================================================!

program main

use io, only: stdin, stdout, stderr, fileExists, line_count
use geom, only: pnpoly

use clip, only: poly_file, &
                input_file, &
                output_file, &
                getPointsInsideBoundary, &
                getPointsOutsideBoundary, &
                getPointsOnBoundary, &
                epsilon, &
                printRecord

implicit none

! Local variables
integer :: i, ios, luin, luout, nclip, inout, irecord
character(len=512) :: line
double precision :: x, y
double precision, allocatable :: polygon(:,:)


! Parse command line
call gcmdln()

! Make sure something needs to be printed...
if (.not.getPointsInsideBoundary .and. &
    .not.getPointsOutsideBoundary .and. &
    .not.getPointsOnBoundary) then
    call usage('clip: no output mode specified')
endif


! Read clipping polygon
if(.not.fileExists(poly_file)) then
    call usage('clip: no polygon file found named "'//trim(poly_file)//'"')
endif
nclip = line_count(poly_file)
allocate(polygon(nclip,2))
open(unit=21,file=poly_file,status='old')
do i = 1,nclip
    read(21,'(A)') line
    read(line,*,iostat=ios) polygon(i,1),polygon(i,2)
    if (ios.ne.0) then
        write(stderr,*) 'clip: error reading polygon file coordinates at line ',i
        call usage(     'Offending line: '//trim(line))
    endif
enddo
close(21)


! Determine input/output streams
if (input_file.eq.''.or.input_file.eq.'none'.or.input_file.eq.'stdin') then
    luin = stdin
else
    if (.not.fileExists(input_file)) then
        call usage('clip: no file found named "'//trim(input_file)//'"')
    endif
    luin = 11
    open(unit=luin,file=input_file,status='old')
endif

if (output_file.eq.''.or.output_file.eq.'none'.or.output_file.eq.'stdout') then
    luout = stdout
else
    luout = 12
    open(unit=luout,file=output_file,status='unknown')
endif


! Read input data and clip it
irecord = 0
do
    read(luin,'(A)',iostat=ios) line
    if (ios.ne.0) then
        exit
    endif

    irecord = irecord + 1
    read(line,*,iostat=ios) x,y
    if (ios.ne.0) then
        write(stderr,*) 'clip: error reading input file coordinates'
        call usage(     'Offending line: '//trim(line))
    endif

    call pnpoly(x,y,polygon(:,1),polygon(:,2),nclip,inout,epsilon)

    if (getPointsInsideBoundary.and.inout.eq.1) then
        if (printRecord) then
            write(luout,*) irecord,trim(line)
        else
            write(luout,*) trim(line)
        endif
    elseif (getPointsOutsideBoundary.and.inout.eq.-1) then
        if (printRecord) then
            write(luout,*) irecord,trim(line)
        else
            write(luout,*) trim(line)
        endif
    elseif (getPointsOnBoundary.and.inout.eq.0) then
        if (printRecord) then
            write(luout,*) irecord,trim(line)
        else
            write(luout,*) trim(line)
        endif
    endif
enddo



if (luin.ne.stdin) then
    close(11)
endif
if (luout.ne.stdout) then
    close(12)
endif

end

!--------------------------------------------------------------------------------------------------!

subroutine gcmdln()

use clip, only: poly_file, &
                input_file, &
                output_file, &
                getPointsInsideBoundary, &
                getPointsOutsideBoundary, &
                getPointsOnBoundary, &
                epsilon, &
                printRecord

implicit none

! Local variables
integer :: i, narg
character(len=512) :: tag


poly_file = ''
input_file = ''
output_file = ''
getPointsInsideBoundary = .false.
getPointsOutsideBoundary = .false.
getPointsOnBoundary = .false.
epsilon = 1.0d-8
printRecord = .false.

! Number of arguments
narg = command_argument_count()
if (narg.eq.0) then
  call usage('')
endif

! Read arguments
i = 1
call get_command_argument(i,poly_file)
i = i + 1

do while (i.le.narg)

    call get_command_argument(i,tag)

    if (trim(tag).eq.'-f') then
        i = i + 1
        call get_command_argument(i,input_file)

    elseif (tag.eq.'-o') then
        i = i + 1
        call get_command_argument(i,output_file)

    elseif (tag.eq.'-in') then
        getPointsInsideBoundary = .true.

    elseif (tag.eq.'-out') then
        getPointsOutsideBoundary = .true.

    elseif (tag.eq.'-on') then
        getPointsOnBoundary = .true.
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) epsilon

    elseif (tag.eq.'-NR') then
        printRecord = .true.

    else
        call usage('clip: no option '//trim(tag))
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

write(stderr,*) 'Usage: clip CLIP_FILE ...options...'
write(stderr,*)
write(stderr,*) 'CLIP_FILE         Clipping boundary'
write(stderr,*) '-f FILE           X-Y-Z data to clip (default: stdin)'
write(stderr,*) '-o FILE           Clipped data (default: stdout)'
write(stderr,*) '-in               Keep points inside clipping boundary'
write(stderr,*) '-out              Keep points outside clipping boundary'
write(stderr,*) '-on EPS           Keep points within EPS of boundary'
write(stderr,*) '-NR               Print record of clipped lines'
write(stderr,*)

stop
end subroutine
