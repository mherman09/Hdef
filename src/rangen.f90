module rangen

character(len=32) :: ran_mode
double precision :: x1
double precision :: x2
double precision :: mean
double precision :: stdev
integer :: npts

end module

!==================================================================================================!

program main

use io, only: stdout
use random, only: iseed, r8_uniform_01, r8_normal_ab

use rangen, only: ran_mode, &
                  x1, &
                  x2, &
                  mean, &
                  stdev, &
                  npts

implicit none

! Local variables
integer :: i


call gcmdln()
if (ran_mode.eq.'') then
    call usage('rangen: no random number generator mode defined')
endif

do i = 1,npts
    if (ran_mode.eq.'uniform') then
        write(stdout,*) x1 + (x2-x1)*r8_uniform_01(iseed)
    elseif (ran_mode.eq.'normal') then
        write(stdout,*) r8_normal_ab(mean,stdev,iseed)
    else
        call usage('rangen: no generator mode named '//trim(ran_mode))
    endif
enddo

end

!--------------------------------------------------------------------------------------------------!

subroutine gcmdln()

use random, only: iseed, timeseed

use rangen, only: ran_mode, &
                  x1, &
                  x2, &
                  mean, &
                  stdev, &
                  npts

implicit none

! Local variables
integer :: i, narg
character(len=512) :: tag

! Initialize control parameters
ran_mode = ''
x1 = 0.0d0
x2 = 1.0d0
mean = 0.0d0
stdev = 1.0d0
npts = 100
iseed = timeseed()

! Number of arguments
narg = command_argument_count()
if (narg.eq.0) then
  call usage('')
endif

! Read arguments
i = 1
do while (i.le.narg)

    call get_command_argument(i,tag)

    if (trim(tag).eq.'-uniform') then
        ran_mode = 'uniform'
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) x1
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) x2

    elseif (trim(tag).eq.'-normal') then
        ran_mode = 'normal'
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) mean
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) stdev

    elseif (tag.eq.'-n'.or.tag.eq.'-npts') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) npts

    elseif (trim(tag).eq.'-seed') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) iseed

    else
        call usage('rangen: no option '//trim(tag))
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

write(stderr,*) 'Usage: rangen -uniform X1 X2 | -normal MEAN STD  [-n NPTS] [-seed SEED]'
write(stderr,*)
write(stderr,*) '-uniform X1 X2        Uniform distribution between X1 and X2'
write(stderr,*) '-normal MEAN STD      Gaussian distribution centered on MEAN with std dev STD'
write(stderr,*) '-n NPTS               Number of points to generate (default: 100)'
write(stderr,*) '-seed SEED            Random number seed'
write(stderr,*)

stop
end subroutine
