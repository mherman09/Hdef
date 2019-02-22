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

public :: stdin
public :: stdout
public :: stderr
public :: verbosity

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
!----

implicit none

! I/O variables
character(len=*) :: file_name
integer :: line_count

! Local variables
integer :: ios

line_count= 0

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

end module
