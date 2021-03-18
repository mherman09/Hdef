module mtutil

character(len=32) :: mtutil_mode
character(len=512) :: input
character(len=32) :: input_mode
character(len=512) :: output
integer :: nvals_input
integer :: nvals_output
character(len=5) :: scalar_moment_mode

end module

!==================================================================================================!

program main
!----
! Program for manipulating seismic source data, including:
!     Moment tensors
!     P, N, and T axes
!     Double couple parameters
!     Seismic moment/Moment magnitude
!     Ternary characterization (after Frohlich, 1992)
!----

use io, only: stdin, stdout, stderr, fileExists, line_count
use eq, only: mag2mom, &
              mom2mag, &
              mij2dcp, &
              mij2mag, &
              mij2mom, &
              mij2pnt, &
              mij2sdr, &
              mij2sv, &
              mij2ter, &
              pnt2dcp, &
              pnt2mag, &
              pnt2mag_nondc, &
              pnt2mij, &
              pnt2mom, &
              pnt2mom_nondc, &
              pnt2sdr, &
              pnt2sv, &
              pnt2ter, &
              sdr2mij, &
              sdr2pnt, &
              sdr2sdr, &
              sdr2sv, &
              sdr2ter, &
              sdr2kagan, &
              mij2kagan

use mtutil, only: mtutil_mode, &
                  input, &
                  input_mode, &
                  output, &
                  nvals_input, &
                  nvals_output, &
                  scalar_moment_mode
implicit none

! Local variables
integer :: i, j, ninputs, luin, luout, ios
double precision :: input_values(12), input_values2(12), output_values(12)
character(len=512) :: tmp_file, input_line


! Initialize variables
input_values = 0.0d0
output_values = 0.0d0
tmp_file = 'mtutil_86_this_file_when_program_exits.tmp'


! Parse the command line
call gcmdln()


! Count the number of input lines and set the input unit
if (input_mode.eq.'file') then
    luin = 11
    if (.not.fileExists(input)) then
        call usage('mtutil: input file "'//trim(input)//'" not found')
    endif
    ninputs = line_count(input)
elseif (input_mode.eq.'cmdln') then
    luin = 11
    ! Replace commas with spaces in command line argument
    i = index(input,',')
    do while (i.ne.0)
        input(i:i) = ' '
        i = index(input,',')
    enddo
    ! Write command line info to temporary input file
    open(unit=luin,file=tmp_file,status='unknown')
    write(luin,*) input
    close(luin)
    input = tmp_file
    ninputs = 1
elseif (input_mode.eq.'stdin') then
    luin = stdin
    ninputs = 100000000
else
    call usage('mtutil: no input_mode named "'//trim(input_mode)//'"')
endif


! Open the input file if necessary
if (input_mode.eq.'file'.or.input_mode.eq.'cmdln') then
    open(unit=luin,file=input,status='old')
endif

! Set the output unit
if (output.eq.'stdout') then
    luout = stdout
else
    luout = 12
    open(unit=luout,file=output,status='unknown')
endif


! Read the inputs and run specified calculation
do i = 1,ninputs
    read(luin,'(A)',end=1001) input_line
    if (nvals_input.ge.1) then
        read(input_line,*,iostat=ios) (input_values(j),j=1,nvals_input)
    elseif (nvals_input.le.-1) then
        read(input_line,*,iostat=ios) (input_values(j),j=1,abs(nvals_input)), &
                                      (input_values2(j),j=1,abs(nvals_input))
    else
        write(stderr,*) 'nvals_input set to illegal value'
    endif
    if (ios.ne.0) then
        write(stderr,*) 'mtutil: error parsing inputs "',trim(input_line),'"'
        call usage('')
    endif
    ! write(0,*) input_values(1:nvals_input)

    ! Convert between moment and magnitude
    if (mtutil_mode.eq.'mag2mom') then
        call mag2mom(input_values(1),output_values(1))
    elseif (mtutil_mode.eq.'mom2mag') then
        call mom2mag(input_values(1),output_values(1))

    ! Input quantity: moment tensor
    elseif (mtutil_mode.eq.'mij2dcp') then
        call mij2dcp(input_values(1),input_values(2),input_values(3),&
                     input_values(4),input_values(5),input_values(6),&
                     output_values(1))
    elseif (mtutil_mode.eq.'mij2mag') then
        call mij2mag(input_values(1),input_values(2),input_values(3),&
                     input_values(4),input_values(5),input_values(6),&
                     output_values(1),scalar_moment_mode)
        if (abs(output_values(1)).le.-99.0d0) then
            call usage('mtutil: no scalar moment mode called '//trim(scalar_moment_mode))
        endif
    elseif (mtutil_mode.eq.'mij2mom') then
        call mij2mom(input_values(1),input_values(2),input_values(3),&
                     input_values(4),input_values(5),input_values(6),&
                     output_values(1),scalar_moment_mode)
        if (abs(output_values(1)).le.1.0d-50) then
            call usage('mtutil: no scalar moment mode called '//trim(scalar_moment_mode))
        endif
    elseif (mtutil_mode.eq.'mij2pnt') then
        call mij2pnt(input_values(1),input_values(2),input_values(3),&
                     input_values(4),input_values(5),input_values(6),&
                     output_values)
    elseif (mtutil_mode.eq.'mij2sdr') then
        nvals_output = 6
        call mij2sdr(input_values(1),input_values(2),input_values(3),&
                     input_values(4),input_values(5),input_values(6),&
                     output_values(1),output_values(2),output_values(3),&
                     output_values(4),output_values(5),output_values(6))
    elseif (mtutil_mode.eq.'mij2sv') then
        nvals_output = 6
        call mij2sv(input_values(1),input_values(2),input_values(3),&
                    input_values(4),input_values(5),input_values(6),&
                    output_values(1:3),output_values(4:6))
    elseif (mtutil_mode.eq.'mij2ter') then
        call mij2ter(input_values(1),input_values(2),input_values(3),&
                     input_values(4),input_values(5),input_values(6),&
                     output_values(1),output_values(2),output_values(3))

    ! Input quantity: moment tensor eigenvectors/eigenvalues (P, N, T axes)
    elseif (mtutil_mode.eq.'pnt2dcp') then
        call pnt2dcp(input_values,output_values(1))
    elseif (mtutil_mode.eq.'pnt2mag') then
        if (scalar_moment_mode.eq.'dc') then
            call pnt2mag(input_values,output_values(1))
        elseif (scalar_moment_mode.eq.'nondc') then
            call pnt2mag_nondc(input_values,output_values(1))
        else
            call usage('mtutil: no scalar moment mode called '//trim(scalar_moment_mode))
        endif
    elseif (mtutil_mode.eq.'pnt2mij') then
        call pnt2mij(input_values,&
                     output_values(1),output_values(2),output_values(3),&
                     output_values(4),output_values(5),output_values(6))
    elseif (mtutil_mode.eq.'pnt2mom') then
        if (scalar_moment_mode.eq.'dc') then
            call pnt2mom(input_values,output_values(1))
        elseif (scalar_moment_mode.eq.'nondc') then
            call pnt2mom_nondc(input_values,output_values(1))
        else
            call usage('mtutil: no scalar moment mode called '//trim(scalar_moment_mode))
        endif
    elseif (mtutil_mode.eq.'pnt2sdr') then
        nvals_output = 6
        call pnt2sdr(input_values,&
                     output_values(1),output_values(2),output_values(3),&
                     output_values(4),output_values(5),output_values(6))
    elseif (mtutil_mode.eq.'pnt2sv') then
        nvals_output = 6
        call pnt2sv(input_values,output_values(1:3),output_values(4:6))
    elseif (mtutil_mode.eq.'pnt2ter') then
        call pnt2ter(input_values,output_values(1),output_values(2),output_values(3))

    ! Input quantity: double couple strike, dip, rake
    elseif (mtutil_mode.eq.'sdr2mij') then
        call sdr2mij(input_values(1),input_values(2),input_values(3),&
                     output_values(1),output_values(2),output_values(3),&
                     output_values(4),output_values(5),output_values(6))
    elseif (mtutil_mode.eq.'sdr2pnt') then
        call sdr2pnt(input_values(1),input_values(2),input_values(3),output_values)
    elseif (mtutil_mode.eq.'sdr2sdr') then
        nvals_output = 3
        call sdr2sdr(input_values(1),input_values(2),input_values(3),&
                     output_values(1),output_values(2),output_values(3))
    elseif (mtutil_mode.eq.'sdr2sv') then
        nvals_output = 3
        call sdr2sv(input_values(1),input_values(2),input_values(3),&
                    output_values(1:3))
    elseif (mtutil_mode.eq.'sdr2ter') then
        call sdr2ter(input_values(1),input_values(2),input_values(3),&
                     output_values(1),output_values(2),output_values(3))


    ! output -> kagan angle
    elseif (mtutil_mode.eq.'sdr2kagan') then
        call sdr2kagan(input_values(1), input_values(2), input_values(3), &
                       input_values2(1),input_values2(2),input_values2(3), &
                       output_values(1))
    elseif (mtutil_mode.eq.'mij2kagan') then
        call mij2kagan(input_values(1), input_values(2), input_values(3), &
                       input_values(4), input_values(5), input_values(6), &
                       input_values2(1),input_values2(2),input_values2(3), &
                       input_values2(4),input_values2(5),input_values2(6), &
                       output_values(1))

    else
        call usage('mtutil: no option '//trim(mtutil_mode))
    endif
    write(luout,*) output_values(1:nvals_output)
enddo
1001 continue


! Close the input file and delete it if necessary
if (input_mode.eq.'file') then
    close(luin)
elseif (input_mode.eq.'cmdln') then
    close(luin,status='delete')
endif

! Close the output file
if (output.ne.'stdout') then
    close(luout)
endif

end

!--------------------------------------------------------------------------------------------------!

subroutine gcmdln()

use io, only: isNumeric
use mtutil, only: mtutil_mode, &
                  input, &
                  input_mode, &
                  output, &
                  nvals_input, &
                  nvals_output, &
                  scalar_moment_mode
implicit none

! Local variables
integer :: i, j, narg
character(len=512) :: tag

! Initialize
mtutil_mode = ''
input = ''
input_mode = ''
output = 'stdout'
nvals_input = 0
nvals_output = 0
scalar_moment_mode = 'dc'

! Number of arguments
narg = command_argument_count()
if (narg.eq.0) then
  call usage('')
endif


! First argument defines input parameters
i = 1
call get_command_argument(i,tag)
if (trim(tag).eq.'-sdr') then
    mtutil_mode = 'sdr2'
    nvals_input = 3
elseif (trim(tag).eq.'-mij') then
    mtutil_mode = 'mij2'
    nvals_input = 6
elseif (trim(tag).eq.'-pnt') then
    mtutil_mode = 'pnt2'
    nvals_input = 12
elseif (trim(tag).eq.'-mag') then
    mtutil_mode = 'mag2'
    nvals_input = 1
elseif (trim(tag).eq.'-mom') then
    mtutil_mode = 'mom2'
    nvals_input = 1
else
    call usage('mtutil: no input option '//tag)
endif
! write(0,*) trim(mode), nvals_input


! Optional argument defines a file name or value(s) for input
i = i + 1
if (i.gt.narg) then
    call usage('mtutil: input and output are required')
endif
! Determine how to read the second argument
call get_command_argument(i,input)
j = index(input,',')
if (j.ne.0) then
    ! comma => list of numbers to be read directly
    input_mode = 'cmdln'
    i = i + 1
else
    ! no comma => single number OR file name OR blank
    if (isNumeric(input)) then
        ! argument is readable as floating point number => read directly
        input_mode = 'cmdln'
        i = i + 1
    else
        ! argument is not readable as floating point number
        j = index(input,'-')
        if (j.eq.1) then
            ! argument starts with dash => mtutil flag
            input_mode = 'stdin'
        else
            ! argument does not start with dash => file name
            input_mode = 'file'
            i = i + 1
        endif
    endif
endif
! Make sure some input mode is specified
if (input_mode.eq.'') then
    call usage('mtutil: no input mode specified')
endif
! write(0,*) trim(input_mode)


! Second argument defines output parameters
if (i.gt.narg) then
    call usage('mtutil: an output option is required')
endif
call get_command_argument(i,tag)
if (trim(tag).eq.'-sdr') then
    mtutil_mode = trim(mtutil_mode)//'sdr'
    nvals_output = 3
elseif (trim(tag).eq.'-mij') then
    mtutil_mode = trim(mtutil_mode)//'mij'
    nvals_output = 6
elseif (trim(tag).eq.'-pnt') then
    mtutil_mode = trim(mtutil_mode)//'pnt'
    nvals_output = 12
elseif (trim(tag).eq.'-mag') then
    mtutil_mode = trim(mtutil_mode)//'mag'
    nvals_output = 1
elseif (trim(tag).eq.'-mom') then
    mtutil_mode = trim(mtutil_mode)//'mom'
    nvals_output = 1
elseif (trim(tag).eq.'-dcp') then
    mtutil_mode = trim(mtutil_mode)//'dcp'
    nvals_output = 1
elseif (trim(tag).eq.'-ternary') then
    mtutil_mode = trim(mtutil_mode)//'ter'
    nvals_output = 3
elseif (trim(tag).eq.'-sv') then
    mtutil_mode = trim(mtutil_mode)//'sv'
    nvals_output = 6
elseif (trim(tag).eq.'-kagan') then
    mtutil_mode = trim(mtutil_mode)//'kagan'
    nvals_output = 1
    nvals_input = -nvals_input
else
      call usage('mtutil: no output option '//tag)
endif
! write(0,*) trim(mode)


! Optional argument defines a file name for output
i = i + 1
if (i.le.narg) then
    call get_command_argument(i,output)
    j = index(output,'-')
    if (j.eq.1) then
        ! argument starts with dash => mtutil flag, so send output to stdout
        output = 'stdout'
    else
        ! argument does not start with dash => output file name
        i = i + 1
    endif
endif


! Extra arguments
do while (i.le.narg)
    call get_command_argument(i,tag)
    if (trim(tag).eq.'-dc') then
        scalar_moment_mode = 'dc'
    elseif (trim(tag).eq.'-non-dc') then
        scalar_moment_mode = 'nondc'
    else
        call usage('mtutil: no option '//trim(tag))
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

write(stderr,*) 'Usage: mtutil -opt [INPUT] -opt [FILE] [...OPTIONS...]'
write(stderr,*)
write(stderr,*) 'Input/output:'
write(stderr,*) '-sdr      Strike, dip, and rake angles'
write(stderr,*) '-mij      Moment tensor components (mrr mtt mpp mrt mrp mtp)'
write(stderr,*) '-pnt      P, N, T axes (px,py,pz,nx,ny,nz,tx,ty,tz,pmag,nmag,tmag)'
write(stderr,*) '-mag      Moment magnitude'
write(stderr,*) '-mom      Seismic moment (Nm)'
write(stderr,*)
write(stderr,*) 'Output only:'
write(stderr,*) '-ternary  Ternary classification (fth fss fno)'
write(stderr,*) '-dcp      Double couple percentage'
write(stderr,*) '-sv       Slip vector (e n z components)'
write(stderr,*) '-kagan    Calculate angle between two moment tensors (requires double inputs)'
write(stderr,*)
write(stderr,*) 'Other options:'
write(stderr,*) '-dc       Use the double couple scalar seismic moment for -mom/-mag (default)'
write(stderr,*) '-non-dc   Use the general scalar seismic moment for -mom/-mag'
write(stderr,*)
write(stderr,*) 'Input format options:'
write(stderr,*) 'file name: read from file INPUT'
write(stderr,*) 'comma delimited list: read from command line'
write(stderr,*) '"stdin" or blank: read from standard input'
write(stderr,*)
write(stderr,*) 'Output format options:'
write(stderr,*) 'file name: write to file OUTPUT'
write(stderr,*) 'blank: print to standard output'
write(stderr,*)

call error_exit(1)
end subroutine
