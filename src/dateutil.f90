module dateutil

character(len=512) :: input_file
character(len=512) :: output_file
character(len=32) :: date_option
character(len=32) :: date_format

end module

!==================================================================================================!

program main
!----
! Utility for converting to/from calendar dates and calculating numbers of days
!----

use io, only: stdin, stdout, fileExists
use calendar, only: parse_date, date2jd, jd2date

use dateutil, only: input_file, &
                    output_file, &
                    date_option, &
                    date_format

implicit none

! Local variables
integer :: i, luin, luout, ios, date1(7), date2(7), ierr
character(len=512) :: input_line, arg1, arg2
character(len=32) :: string(14)
double precision :: nday, day1, day2


! Parse command line
call gcmdln()

! Check that all necessary inputs/outputs defined
if (date_option.eq.'') then
    call usage('dateutil: no date computation option defined')
endif

! Check that input file exists, if defined
if (input_file.ne.'') then
    if (.not.fileExists(input_file)) then
        call usage('dateutil: no input file found named '//trim(input_file))
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


! Calculate the new date or the number of days between two dates

do
    ! Read input line
    read(luin,'(A)',iostat=ios) input_line
    if (ios.ne.0) then
        exit
    endif

    ! Split it into two arguments, which depends on formatting
    ! This is an annoying block...
    if (date_format.eq.'YYYY-MM-DD') then
        read(input_line,*) arg1, arg2

    elseif (date_format.eq.'YYYY MM DD') then
        if (date_option.eq.'nday') then
            read(input_line,*) (string(i),i=1,6)
            write(arg1,*) (string(i),i=1,3)
            write(arg2,*) (string(i),i=4,6)
        elseif (date_option.eq.'date') then
            read(input_line,*) (string(i),i=1,4)
            write(arg1,*) (string(i),i=1,3)
            write(arg2,*) string(4)
        endif

    elseif (date_format.eq.'YYYYMMDD') then
        read(input_line,*) arg1, arg2

    elseif (date_format.eq.'YYYY-MM-DDTHH:MM:SS') then
        read(input_line,*) arg1, arg2

    elseif (date_format.eq.'YYYY MM DD HH MM SS') then
        if (date_option.eq.'nday') then
            read(input_line,*) (string(i),i=1,12)
            write(arg1,*) (string(i),i=1,6)
            write(arg2,*) (string(i),i=7,12)
        elseif (date_option.eq.'date') then
            read(input_line,*) (string(i),i=1,7)
            write(arg1,*) (string(i),i=1,6)
            write(arg2,*) string(7)
        endif

    elseif (date_format.eq.'YYYYMMDDHHMMSS') then
        read(input_line,*) arg1, arg2

    elseif (date_format.eq.'YYYY-MM-DDTHH:MM:SS.MMM') then
        read(input_line,*) arg1, arg2

    elseif (date_format.eq.'YYYY MM DD HH MM SS MMM') then
        if (date_option.eq.'nday') then
            read(input_line,*) (string(i),i=1,14)
            write(arg1,*) (string(i),i=1,7)
            write(arg2,*) (string(i),i=8,14)
        elseif (date_option.eq.'date') then
            read(input_line,*) (string(i),i=1,8)
            write(arg1,*) (string(i),i=1,7)
            write(arg2,*) string(8)
        endif

    elseif (date_format.eq.'YYYYMMDDHHMMSSMMM') then
        read(input_line,*) arg1, arg2

    else
        call usage('dateutil: date_format "'//trim(date_format)//'" not recognized')
    endif


    ! The first argument is always a date
    call parse_date(arg1,date1,date_format,ierr)
    call date2jd(date1,day1)

    ! The second argument depends on the mode
    if (date_option.eq.'nday') then
        call parse_date(arg2,date2,date_format,ierr)
        call date2jd(date2,day2)
        write(luout,*) day2-day1
    elseif (date_option.eq.'date') then
        read(arg2,*) nday
        day2 = day1+nday
        call jd2date(day2,date2)
        if (len(trim(date_format)).le.10) then
            write(luout,1001) (date2(i),i=1,3)
        else
            write(luout,1002) (date2(i),i=1,7)
        endif
        1001 format(I0.4,"-",I0.2,"-",I0.2)
        1002 format(I0.4,"-",I0.2,"-",I0.2,"T",I0.2,":",I0.2,":",I0.2,".",I0.3)
    endif
enddo


! Shut it down
if (input_file.ne.'') then
    close(luin)
endif
if (output_file.ne.'') then
    close(luout)
endif

end

!--------------------------------------------------------------------------------------------------!

subroutine gcmdln()

use dateutil, only: input_file, &
                    output_file, &
                    date_option, &
                    date_format

implicit none

! Local variables
integer :: i, narg
character(len=512) :: tag

! Initialize control parameters
input_file = ''
output_file = ''
date_option = ''
date_format = 'YYYY-MM-DDTHH:MM:SS'

! Number of arguments
narg = command_argument_count()
if (narg.eq.0) then
  call usage('')
endif

! Read arguments
i = 1
do while (i.le.narg)

    call get_command_argument(i,tag)

    if (trim(tag).eq.'-nday') then
        date_option = 'nday'
    elseif (tag(1:5).eq.'-date') then
        date_option = 'date'
    elseif (trim(tag).eq.'-f') then
        i = i + 1
        call get_command_argument(i,input_file)
    elseif (trim(tag).eq.'-o') then
        i = i + 1
        call get_command_argument(i,output_file)
    elseif (trim(tag).eq.'-format') then
        i = i + 1
        call get_command_argument(i,date_format)
    else
        call usage('dateutil: No option '//tag)
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

write(stderr,*) 'Usage: dateutil ...options...'
write(stderr,*)
write(stderr,*) '-f IFILE      Input file (default reads stdin)'
write(stderr,*) '-o OFILE      Output file (default prints to stdout)'
write(stderr,*) '-nday         Compute number of days from DATE1 DATE2'
write(stderr,*) '-date         Compute calendar date from DATE1 NDAY'
write(stderr,*) '-format FMT   Input date format (default: YYYY-MM-DDTHH:MM:SS)'
write(stderr,*)
write(stderr,*) 'dateutil recognizes the following date formats:'
write(stderr,*) '    YYYY-MM-DD[THH:MM:SS[.MMM]]'
write(stderr,*) '    YYYY MM DD[ HH MM SS[ MMM]]'
write(stderr,*) '    YYYYMMDD[HHMMSS[MMM]]'
write(stderr,*) 'Use quotations when using second format to read spaces correctly'
write(stderr,*)
stop
end subroutine
