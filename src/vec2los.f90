module vec2los

character(len=512) :: input_file
character(len=512) :: output_file
double precision :: azimuth
double precision :: inclination
character(len=512) :: look_geometry_file

end module

!==================================================================================================!

program main
!----
! Read three-dimensional E-N-Z displacements and compute line-of-sight displacements, given the
! azimuth (clockwise from north) and inclination (upward from horizontal) from the viewer to the
! location of interest.
!----

use io, only: stderr, stdin, stdout, fileExists, line_count
use trig, only: d2r

use vec2los, only: input_file, &
                   output_file, &
                   azimuth, &
                   inclination, &
                   look_geometry_file

implicit none

! Local variables
double precision :: az, inc, cosaz, cosinc, sinaz, sininc, ue, un, uz, ulos
double precision, allocatable :: vector(:,:), look(:,:)
integer :: i, luin, luout, ios, nvec, nlook
integer, parameter :: NVECMAX = 10000
character(len=256) :: input_line


! Parse command line
call gcmdln()


! Check that look geometry is defined
if ((azimuth.lt.-1.0d50 .or. inclination.lt.-1.0d50) .and. look_geometry_file.eq.'') then
    call usage('vec2los: no look geometry defined')
endif

! Check that input file exists, if defined
if (input_file.ne.'') then
    if (.not.fileExists(input_file)) then
        call usage('vec2los: no input vector displacement file found named '//trim(input_file))
    endif
endif


! Open input and output streams
if (input_file.eq.'') then
    luin = stdin
    allocate(vector(NVECMAX,6))
else
    nvec = line_count(input_file)
    luin = 11
    open(unit=luin,file=input_file,status='old')
    allocate(vector(nvec,6))
endif

if (output_file.eq.'') then
    luout = stdout
else
    luout = 12
    open(unit=luout,file=output_file,status='unknown')
endif


! Read input displacements
nvec = 0
do
    read(luin,'(A)',iostat=ios) input_line
    if (ios.ne.0) then
        exit
    else
        nvec = nvec + 1
        read(input_line,*) (vector(nvec,i),i=1,6)
    endif

    ! Maximum number of inputs is limited when reading from stdin
    if (input_file.eq.''.and.nvec.gt.NVECMAX) then
        write(stderr,*) 'vec2los: number of inputs from stdin exceeds NVECMAX'
        call usage('read input vectors from file for large numbers of observations')
    endif
enddo


if (look_geometry_file.ne.'') then

    ! If there is a look file, read it

    if (.not.fileExists(look_geometry_file)) then
        call usage('vec2los: no look geometry file found named '//trim(look_geometry_file))
    endif

    ! Number of look values must equal number of inputs
    nlook = line_count(look_geometry_file)
    if (nlook.ne.nvec) then
        call usage('vec2los: number of look geometries is not equal to number of vector inputs')
    endif

    ! Allocate memory for reading inputs
    allocate(look(nlook,2))

    ! Read look geometry from file
    open(unit=21,file=look_geometry_file,status='old')
    do i = 1,nlook
        read(21,*) look(i,1),look(i,2)
    enddo
    close(21)

else

    ! Otherwise set the look geometry to be the same for all displacements

    allocate(look(nvec,2))
    look(:,1) = azimuth
    look(:,2) = inclination

endif


! Compute the line-of-sight displacement, from the viewer to the location
! Positive implies displacement away from viewer
do i = 1,nvec
    az = look(i,1)
    inc = look(i,2)
    cosaz = dcos(look(i,1)*d2r)
    sinaz = dsin(look(i,1)*d2r)
    cosinc = dcos(look(i,2)*d2r)
    sininc = dsin(look(i,2)*d2r)

    ue = vector(i,4)
    un = vector(i,5)
    uz = vector(i,6)
    ulos = un*cosaz*cosinc + ue*sinaz*cosinc - uz*sininc
    write(luout,*) vector(i,1),vector(i,2),vector(i,3),uLOS
enddo


! FINISH HIM...IT...THE PROGRAM
if (input_file.ne.'') then
    close(luin)
endif
if (output_file.ne.'') then
    close(luout)
endif

end

!--------------------------------------------------------------------------------------------------!

subroutine gcmdln()

use vec2los, only: input_file, &
                   output_file, &
                   azimuth, &
                   inclination, &
                   look_geometry_file

implicit none

! Local variables
integer :: i, narg
character(len=512) :: tag

! Initialize control parameters
input_file = ''
output_file = ''
azimuth = -1d100
inclination = -1d100
look_geometry_file = ''

! Number of arguments
narg = command_argument_count()
if (narg.eq.0) then
  call usage('')
endif

! Read arguments
i = 1
do while (i.le.narg)

    call get_command_argument(i,tag)

    if (trim(tag).eq.'-a') then
        i = i + 1
        call get_command_argument(i,tag)
        read (tag,*) azimuth
    elseif (trim(tag).eq.'-i') then
        i = i + 1
        call get_command_argument(i,tag)
        read (tag,*) inclination
    elseif (trim(tag).eq.'-look') then
        i = i + 1
        call get_command_argument(i,look_geometry_file)
    elseif (trim(tag).eq.'-f') then
        i = i + 1
        call get_command_argument(i,input_file)
    elseif (trim(tag).eq.'-o') then
        i = i + 1
        call get_command_argument(i,output_file)
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

write(stderr,*) 'Usage: vec2los -a AZ -i INC | -look FILE [-f IFILE] [-o OFILE]'
write(stderr,*)
write(stderr,*) '-a AZ         Azimuth from viewer to ground (CW from N)'
write(stderr,*) '-i INC        Inclination from viewer to ground (horizontal=0, vertical=90)'
write(stderr,*) '-look FILE    Read azimuth and inclination from FILE'
write(stderr,*) '-f FILE       Input vector displacements (default: stdin) (lo la dp E N Z)'
write(stderr,*) '-o FILE       Output LOS displacements (default: stdout) (lo la dp LOS)'
write(stderr,*)

stop
end subroutine
