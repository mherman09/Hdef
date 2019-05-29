module test

interface test_value
    module procedure test_value_dp
    module procedure test_value_int
    module procedure test_value_real
end interface

!--------------------------------------------------------------------------------------------------!
contains
!--------------------------------------------------------------------------------------------------!

subroutine test_value_dp(actual_value,expected_value,string)
implicit none

! Arguments
double precision :: actual_value, expected_value
character(len=*) :: string

! Local variables
double precision :: diff

if (dabs(expected_value).lt.1.0d-10) then
    diff = actual_value - expected_value
else
    diff = (actual_value-expected_value)/actual_value
endif

if (dabs(diff).gt.1.0d-10) then
    write(0,*) 'test FAILED for '//trim(string)
    write(0,*) '    expected value: ',expected_value
    write(0,*) '    computed value: ',actual_value
    call error_exit(1)
else
    write(6,*) 'test passed for '//trim(string)
endif
return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine test_value_int(actual_value,expected_value,string)
implicit none

! Arguments
integer :: actual_value, expected_value
character(len=*) :: string

if (actual_value.ne.expected_value) then
    write(0,*) 'test FAILED for '//trim(string)
    write(0,*) '    expected value: ',expected_value
    write(0,*) '    computed value: ',actual_value
    call error_exit(1)
else
    write(6,*) 'test passed for '//trim(string)
endif
return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine test_value_real(actual_value,expected_value,string)
implicit none

! Arguments
real :: actual_value, expected_value
character(len=*) :: string

! Local variables
real :: diff

if (abs(expected_value).lt.1.0e-10) then
    diff = actual_value - expected_value
else
    diff = (actual_value-expected_value)/actual_value
endif

if (abs(diff).gt.1.0e-7) then
    write(0,*) 'test FAILED for '//trim(string)
    write(0,*) '    expected value: ',expected_value
    write(0,*) '    computed value: ',actual_value
    call error_exit(1)
else
    write(6,*) 'test passed for '//trim(string)
endif
return
end subroutine

end module
