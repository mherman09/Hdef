program main

use test, only: test_value
use io, only: stdout

use sort

implicit none

double precision :: array(10), array2(10)
integer :: array_int(10), array2_int(10)
real :: array_real(10), array2_real(10)

array(1)  =  14.060924487405954d0
array(2)  =  21.957872304818064d0
array(3)  =  45.959846554673334d0
array(4)  =  47.141085165458485d0
array(5)  =  0.21841768184381499d0
array(6)  =  70.945978939022865d0
array(7)  =  89.068091092152471d0
array(8)  =  67.407064819599611d0
array(9)  =  10.538482810065107d0
array(10) =  20.280598112373479d0
array2 = array
call quicksort(array2,10)
call test_value(array2(1), 0.21841768184381499d0,'quicksort_dp: array(1)')
call test_value(array2(2), 10.538482810065107d0, 'quicksort_dp: array(2)')
call test_value(array2(3), 14.060924487405954d0, 'quicksort_dp: array(3)')
call test_value(array2(4), 20.280598112373479d0, 'quicksort_dp: array(4)')
call test_value(array2(5), 21.957872304818064d0, 'quicksort_dp: array(5)')
call test_value(array2(6), 45.959846554673334d0, 'quicksort_dp: array(6)')
call test_value(array2(7), 47.141085165458485d0, 'quicksort_dp: array(7)')
call test_value(array2(8), 67.407064819599611d0, 'quicksort_dp: array(8)')
call test_value(array2(9), 70.945978939022865d0, 'quicksort_dp: array(9)')
call test_value(array2(10),89.068091092152471d0, 'quicksort_dp: array(10)')
array2 = array
call heapsort(array2,10)
call test_value(array2(1), 0.21841768184381499d0,'heapsort_dp: array(1)')
call test_value(array2(2), 10.538482810065107d0, 'heapsort_dp: array(2)')
call test_value(array2(3), 14.060924487405954d0, 'heapsort_dp: array(3)')
call test_value(array2(4), 20.280598112373479d0, 'heapsort_dp: array(4)')
call test_value(array2(5), 21.957872304818064d0, 'heapsort_dp: array(5)')
call test_value(array2(6), 45.959846554673334d0, 'heapsort_dp: array(6)')
call test_value(array2(7), 47.141085165458485d0, 'heapsort_dp: array(7)')
call test_value(array2(8), 67.407064819599611d0, 'heapsort_dp: array(8)')
call test_value(array2(9), 70.945978939022865d0, 'heapsort_dp: array(9)')
call test_value(array2(10),89.068091092152471d0, 'heapsort_dp: array(10)')

array_real(1)  =  14.0609245e0
array_real(2)  =  21.9578723e0
array_real(3)  =  45.9598466e0
array_real(4)  =  47.1410852e0
array_real(5)  =  0.218417682e0
array_real(6)  =  70.9459789e0
array_real(7)  =  89.0680911e0
array_real(8)  =  67.4070648e0
array_real(9)  =  10.5384828e0
array_real(10) =  20.2805981e0
array2_real = array_real
call quicksort(array2_real,10)
call test_value(array2_real(1), 0.218417682e0,'quicksort_real: array(1)')
call test_value(array2_real(2), 10.5384828e0, 'quicksort_real: array(2)')
call test_value(array2_real(3), 14.0609245e0, 'quicksort_real: array(3)')
call test_value(array2_real(4), 20.2805981e0, 'quicksort_real: array(4)')
call test_value(array2_real(5), 21.9578723e0, 'quicksort_real: array(5)')
call test_value(array2_real(6), 45.9598466e0, 'quicksort_real: array(6)')
call test_value(array2_real(7), 47.1410852e0, 'quicksort_real: array(7)')
call test_value(array2_real(8), 67.4070648e0, 'quicksort_real: array(8)')
call test_value(array2_real(9), 70.9459789e0, 'quicksort_real: array(9)')
call test_value(array2_real(10),89.0680911e0, 'quicksort_real: array(10)')
array2_real = array_real
call heapsort(array2_real,10)
call test_value(array2_real(1), 0.218417682e0,'heapsort_real: array(1)')
call test_value(array2_real(2), 10.5384828e0, 'heapsort_real: array(2)')
call test_value(array2_real(3), 14.0609245e0, 'heapsort_real: array(3)')
call test_value(array2_real(4), 20.2805981e0, 'heapsort_real: array(4)')
call test_value(array2_real(5), 21.9578723e0, 'heapsort_real: array(5)')
call test_value(array2_real(6), 45.9598466e0, 'heapsort_real: array(6)')
call test_value(array2_real(7), 47.1410852e0, 'heapsort_real: array(7)')
call test_value(array2_real(8), 67.4070648e0, 'heapsort_real: array(8)')
call test_value(array2_real(9), 70.9459789e0, 'heapsort_real: array(9)')
call test_value(array2_real(10),89.0680911e0, 'heapsort_real: array(10)')

array_int(1)  =  14
array_int(2)  =  22
array_int(3)  =  46
array_int(4)  =  47
array_int(5)  =   0
array_int(6)  =  71
array_int(7)  =  89
array_int(8)  =  67
array_int(9)  =  11
array_int(10) =  20
array2_int = array_int
call quicksort(array2_int,10)
call test_value(array2_int(1),  0,'quicksort_int: array(1)')
call test_value(array2_int(2), 11, 'quicksort_int: array(2)')
call test_value(array2_int(3), 14, 'quicksort_int: array(3)')
call test_value(array2_int(4), 20, 'quicksort_int: array(4)')
call test_value(array2_int(5), 22, 'quicksort_int: array(5)')
call test_value(array2_int(6), 46, 'quicksort_int: array(6)')
call test_value(array2_int(7), 47, 'quicksort_int: array(7)')
call test_value(array2_int(8), 67, 'quicksort_int: array(8)')
call test_value(array2_int(9), 71, 'quicksort_int: array(9)')
call test_value(array2_int(10),89, 'quicksort_int: array(10)')
array2_int = array_int
call heapsort(array2_int,10)
call test_value(array2_int(1),  0,'heapsort_int: array(1)')
call test_value(array2_int(2), 11, 'heapsort_int: array(2)')
call test_value(array2_int(3), 14, 'heapsort_int: array(3)')
call test_value(array2_int(4), 20, 'heapsort_int: array(4)')
call test_value(array2_int(5), 22, 'heapsort_int: array(5)')
call test_value(array2_int(6), 46, 'heapsort_int: array(6)')
call test_value(array2_int(7), 47, 'heapsort_int: array(7)')
call test_value(array2_int(8), 67, 'heapsort_int: array(8)')
call test_value(array2_int(9), 71, 'heapsort_int: array(9)')
call test_value(array2_int(10),89, 'heapsort_int: array(10)')

write(stdout,*) 'sort_module unit test passed'
end
