module test

!--------------------------------------------------------------------------------------------------!
contains
!--------------------------------------------------------------------------------------------------!

subroutine test_value(actual_value,expected_value,string)
implicit none
double precision :: actual_value, expected_value, diff
character(len=*) :: string

diff = (actual_value-expected_value)/actual_value

if (dabs(diff).gt.1.0d-10) then
    write(0,*) 'test FAILED for '//trim(string)
    write(0,*) '    expected value: ',expected_value
    write(0,*) '    computed value: ',actual_value
    stop
else
    write(6,*) 'test passed for '//trim(string)
endif
return
end subroutine

end module
