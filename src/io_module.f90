module io

integer, parameter :: stdin=5
integer, parameter :: stdout=6
integer, parameter :: stderr=0

integer :: verbosity  ! 0: Silent, only error messages
                      ! 1: Major progress reports
                      ! 2: Detailed progress reports
                      ! 3: Parsed command line
                      ! 4: Parsed inputs
                      ! 5: Outputs
                      ! 6: Detailed debugging output

logical :: debug

public :: stdin
public :: stdout
public :: stderr
public :: verbosity

public :: fileExists
public :: line_count
public :: progress_indicator

!--------------------------------------------------------------------------------------------------!
contains
!--------------------------------------------------------------------------------------------------!

function fileExists(file_name)
!----
! Check that a file exists, return T or F
!----

implicit none

! I/O variables
character(len=*) :: file_name
logical :: fileExists

! Local variables
inquire(file=file_name,exist=fileExists)

return
end function fileExists

!--------------------------------------------------------------------------------------------------!

function line_count(file_name)
!----
! Count the number of lines in a file, return integer value
! Return value:
!     >0: routine executed correctly, value is number of lines in file
!      0: file not found
!     -1: file open in another unit
!----

implicit none

! I/O variables
character(len=*) :: file_name
integer :: line_count

! Local variables
integer :: ios
logical :: iopen

line_count = 0

! Check whether file is already open in another unit
inquire(file=file_name,opened=iopen)
if (iopen) then
    write(stderr,*) 'line_count: file "',trim(file_name),'" is already open in another unit'
    line_count = -1
    return
endif

if (.not.fileExists(file_name)) then
    write(stderr,*) 'line_count: file "',trim(file_name),'" not found'
    return
endif

open(unit=41,file=file_name,status='old')

do
    read(41,*,iostat=ios)
    if (ios.eq.0) then
        line_count = line_count + 1
    else
        exit
    endif
enddo

close(41)

return
end function line_count

!--------------------------------------------------------------------------------------------------!

subroutine progress_indicator(i,n,label,ierr)
!----
! Display a progress indicator showing i steps out of n total.
!----

implicit none

! Arguments
integer :: i, n, ierr
character(len=*) :: label

! Local variables
character(len=1) :: CR

ierr = 0

if (n.le.0) then
    write(stderr,*) 'progress_indicator: input n must be greater than 0'
    ierr = 1
    return
endif

CR = char(13) ! Carriage return

if (n.le.100) then
    write(stdout,1000,advance='no') label,100*i/n,CR
else
    if (mod(i,n/100).eq.0) then
        write(stdout,1000,advance='no') label,100*i/n,CR
    endif
endif

if (i.eq.n) then
    write(stdout,1000) label,100*i/n
endif

1000 format (1X,A,' progress: [',I3,'% complete]',A)

return
end subroutine

end module
