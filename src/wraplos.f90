module wraplos

character(len=512) :: input_file
character(len=512) :: output_file
double precision :: wavelength

end module

!==================================================================================================!

program main
!----
! Convert a line of sight displacement to a phase difference given the wavelength of the observation
!----

use trig, only: pi
use io, only: stdout, stdin, fileExists

use wraplos, only: input_file, &
                   output_file, &
                   wavelength

implicit none

! Local variables
integer :: ios, luin, luout
double precision :: phase, stlo, stla, stdp, ulos


! Parse the command line
call gcmdln()

! Check that wavelength is defined
if (wavelength.le.0.0d0) then
    call usage('wraplos: wavelength is undefined')
endif

! Check that input file exists, if defined
if (input_file.ne.'') then
    if (.not.fileExists(input_file)) then
        call usage('wraplos: no input LOS displacement file found named '//trim(input_file))
    endif
endif


! Open input and output streams
if (input_file.eq.'') then
    luin = stdin
else
    luin = 11
    open(unit=luin,file=input_file,status='old')
endif

if (output_file.eq.'') then
    luout = stdout
else
    luout = 12
    open(unit=luout,file=output_file,status='unknown')
endif


! Convert line of sight displacement to phase difference from satellite
do
    read (luin,*,iostat=ios) stlo, stla, stdp, ulos
    if (ios.ne.0) then
        exit
    endif

    phase = 4.0d0*pi*ulos/wavelength ! phase difference takes into account two-way travel
    phase = dmod(phase,2.0d0*pi)

    ! Keep phase between +/- pi
    do while (phase.le.-pi)
        phase = phase + 2.0d0*pi
    enddo
    do while (phase.ge.pi)
        phase = phase - 2.0d0*pi
    enddo
    write (luout,*) stlo,stla,stdp,phase
enddo

if (input_file.ne.'') then
    close(luin)
endif
if (output_file.ne.'') then
    close(luout)
endif


end

!--------------------------------------------------------------------------------------------------!

subroutine gcmdln()

use wraplos, only: input_file, &
                   output_file, &
                   wavelength

implicit none

! Local variables
integer :: i, narg
character(len=512) :: tag

! Initialize control parameters
input_file = ''
output_file = ''
wavelength = -1d100

! Number of arguments
narg = command_argument_count()
if (narg.eq.0) then
  call usage('')
endif

! Read arguments
i = 1
do while (i.le.narg)

    call get_command_argument(i,tag)

    if (trim(tag).eq.'-w') then
        i = i + 1
        call get_command_argument(i,tag)
        read (tag,*) wavelength
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

write(stderr,*) 'Usage: wraplos -w WAVELENGTH [-f IFILE] [-o OFILE]'
write(stderr,*)
write(stderr,*) '-w WAVELENGTH Wavelength of observation (m)'
write(stderr,*) '-f FILE       Input LOS displacements (default: stdin) (lo la dp LOS)'
write(stderr,*) '-o FILE       Output phase (default: stdout) (lo la dp phase(rad))'
write(stderr,*)

call error_exit(1)
end subroutine
