program main

use test, only: test_value
use io

implicit none

integer :: i

call test_value(stdin,5,'stdin')
call test_value(stdout,6,'stdout')
call test_value(stderr,0,'stderr')

open(unit=11,file='io_module_unit_test',status='unknown')
write(11,*) 1
write(11,*) 2
write(11,*) 3
write(11,*) 4
write(11,*) '#'
write(11,*) '>'
write(11,*) ''
close(11)

if (fileExists('io_module_unit_test')) then
    i = 1
else
    i = 0
endif
call test_value(i,1,'fileExists')

i = line_count('io_module_unit_test')
call test_value(i,7,'line_count')

i = line_count_ignore('io_module_unit_test',2,['>','#'])
call test_value(i,4,'line_count_ignore')

! subroutine progress_indicator()

write(stdout,*) 'io_module unit test passed'
end
