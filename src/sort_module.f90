module sort

interface quicksort
    module procedure quicksort_int
    module procedure quicksort_real
    module procedure quicksort_dp
end interface
interface heapsort
    module procedure heapsort_int
    module procedure heapsort_real
    module procedure heapsort_dp
end interface

!==================================================================================================!
contains
!==================================================================================================!

subroutine quicksort_int(arr, n)
!----
! Sorts an array arr(1:n) into ascending numerical order using the Quicksort algorithm. n is input;
! arr is replaced on output by its sorted rearrangement.
! Parameters: M is the size of subarrays sorted by straight insertion and NSTACK is the required
! auxiliary storage.
!
! Numerical Recipes (Press et al.): 8.2
!----

implicit none

! Arguments
integer :: n
integer :: arr(n)

! Local variables
integer, parameter :: M=7, NSTACK=50
integer :: i, ir, j, jstack, k, l, istack(NSTACK)
integer :: a, temp
logical :: skipIUpdate

! No need to sort 1-item list...
if (n.lt.2) then
    return
endif

! Initialize variables
jstack = 0
l = 1
ir = n

do
    ! Insertion sort when subarray is small enough
    if (ir-l.lt.M) then
        do j = l+1,ir
            a = arr(j)
            skipIUpdate = .false.
            do i = j-1,l,-1
                if (arr(i).le.a) then
                    skipIUpdate = .true.
                    exit
                endif
                arr(i+1) = arr(i)
            enddo
            if (.not.skipIUpdate) then
                i = l - 1
            endif
            arr(i+1) = a
        enddo

        if (jstack.eq.0) then
            return ! Everything has been sorted
        endif

        ! Pop stack and begin a new round of partitioning
        ir = istack(jstack)
        l = istack(jstack-1)
        jstack = jstack - 2

    else

        ! Choose median of left, center, and right elements as partitioning element a
        k = (l+ir)/2 ! median element -> partitioning element

        ! Also rearrange so that a(l) ≤ a(l+1) ≤ a(ir).
        ! Move partitioning element k next to element l
        temp = arr(k)
        arr(k) = arr(l+1)
        arr(l+1) = temp
        ! write(0,*) 'quicksort: rearranging elements'
        if (arr(l).gt.arr(ir)) then
            temp = arr(l)
            arr(l) = arr(ir)
            arr(ir) = temp
        endif
        if (arr(l+1).gt.arr(ir)) then
            temp = arr(l+1)
            arr(l+1) = arr(ir)
            arr(ir) = temp
        endif
        if (arr(l).gt.arr(l+1)) then
            temp = arr(l)
            arr(l) = arr(l+1)
            arr(l+1) = temp
        endif

        ! Initialize pointers for partitioning
        i = l + 1
        j = ir
        a = arr(l+1)   ! Partitioning element

        ! Put all elements < a to its left, all elements > a to its right
        do ! Beginning of innermost loop
            i = i + 1

            ! Scan up to find element > a
            if (arr(i).lt.a) then
                cycle
            endif

            ! Scan down to find element < a
            do while (arr(j).gt.a)
                j = j - 1
            enddo

            ! if (j.lt.i) then
            if (j.le.i) then
                ! Pointers crossed; exit with partitioning complete
                exit
            endif

            ! Exchange elements
            temp = arr(i)
            arr(i) = arr(j)
            arr(j) = temp
        enddo ! End of innermost loop

        ! Insert partitioning element
        arr(l+1) = arr(j)
        arr(j) = a
        jstack = jstack + 2

        ! Push pointers to larger subarray on stack, process smaller subarray immediately
        if (jstack.gt.NSTACK) then
            stop 'quicksort_real: NSTACK too small'
        endif

        if (ir-i+1.ge.j-l) then
            istack(jstack) = ir
            istack(jstack-1) = i
            ir = j - 1
        else
            istack(jstack) = j - 1
            istack(jstack-1) = l
            l = i
        endif
    endif
enddo

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine quicksort_real(arr, n)
!----
! Sorts an array arr(1:n) into ascending numerical order using the Quicksort algorithm. n is input;
! arr is replaced on output by its sorted rearrangement.
! Parameters: M is the size of subarrays sorted by straight insertion and NSTACK is the required
! auxiliary storage.
!
! Numerical Recipes (Press et al.): 8.2
!----

implicit none

! Arguments
integer :: n
real :: arr(n)

! Local variables
integer, parameter :: M=7, NSTACK=50
integer :: i, ir, j, jstack, k, l, istack(NSTACK)
real :: a, temp
logical :: skipIUpdate

! No need to sort 1-item list...
if (n.lt.2) then
    return
endif

! Initialize variables
jstack = 0
l = 1
ir = n

do
    ! Insertion sort when subarray is small enough
    if (ir-l.lt.M) then
        do j = l+1,ir
            a = arr(j)
            skipIUpdate = .false.
            do i = j-1,l,-1
                if (arr(i).le.a) then
                    skipIUpdate = .true.
                    exit
                endif
                arr(i+1) = arr(i)
            enddo
            if (.not.skipIUpdate) then
                i = l - 1
            endif
            arr(i+1) = a
        enddo

        if (jstack.eq.0) then
            return ! Everything has been sorted
        endif

        ! Pop stack and begin a new round of partitioning
        ir = istack(jstack)
        l = istack(jstack-1)
        jstack = jstack - 2

    else

        ! Choose median of left, center, and right elements as partitioning element a
        k = (l+ir)/2 ! median element -> partitioning element

        ! Also rearrange so that a(l) ≤ a(l+1) ≤ a(ir).
        ! Move partitioning element k next to element l
        temp = arr(k)
        arr(k) = arr(l+1)
        arr(l+1) = temp
        ! write(0,*) 'quicksort: rearranging elements'
        if (arr(l).gt.arr(ir)) then
            temp = arr(l)
            arr(l) = arr(ir)
            arr(ir) = temp
        endif
        if (arr(l+1).gt.arr(ir)) then
            temp = arr(l+1)
            arr(l+1) = arr(ir)
            arr(ir) = temp
        endif
        if (arr(l).gt.arr(l+1)) then
            temp = arr(l)
            arr(l) = arr(l+1)
            arr(l+1) = temp
        endif

        ! Initialize pointers for partitioning
        i = l + 1
        j = ir
        a = arr(l+1)   ! Partitioning element

        ! Put all elements < a to its left, all elements > a to its right
        do ! Beginning of innermost loop
            i = i + 1

            ! Scan up to find element > a
            if (arr(i).lt.a) then
                cycle
            endif

            ! Scan down to find element < a
            do while (arr(j).gt.a)
                j = j - 1
            enddo

            ! if (j.lt.i) then
            if (j.le.i) then
                ! Pointers crossed; exit with partitioning complete
                exit
            endif

            ! Exchange elements
            temp = arr(i)
            arr(i) = arr(j)
            arr(j) = temp
        enddo ! End of innermost loop

        ! Insert partitioning element
        arr(l+1) = arr(j)
        arr(j) = a
        jstack = jstack + 2

        ! Push pointers to larger subarray on stack, process smaller subarray immediately
        if (jstack.gt.NSTACK) then
            stop 'quicksort_real: NSTACK too small'
        endif

        if (ir-i+1.ge.j-l) then
            istack(jstack) = ir
            istack(jstack-1) = i
            ir = j - 1
        else
            istack(jstack) = j - 1
            istack(jstack-1) = l
            l = i
        endif
    endif
enddo

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine quicksort_dp(arr, n)
!----
! Sorts an array arr(1:n) into ascending numerical order using the Quicksort algorithm. n is input;
! arr is replaced on output by its sorted rearrangement.
! Parameters: M is the size of subarrays sorted by straight insertion and NSTACK is the required
! auxiliary storage.
!
! Numerical Recipes (Press et al.): 8.2
!----

implicit none

! Arguments
integer :: n
double precision :: arr(n)

! Local variables
integer, parameter :: M=7, NSTACK=50
integer :: i, ir, j, jstack, k, l, istack(NSTACK)
double precision :: a, temp
logical :: skipIUpdate

! No need to sort 1-item list...
if (n.lt.2) then
    return
endif

! Initialize variables
jstack = 0
l = 1
ir = n

do
    ! Insertion sort when subarray is small enough
    if (ir-l.lt.M) then
        do j = l+1,ir
            a = arr(j)
            skipIUpdate = .false.
            do i = j-1,l,-1
                if (arr(i).le.a) then
                    skipIUpdate = .true.
                    exit
                endif
                arr(i+1) = arr(i)
            enddo
            if (.not.skipIUpdate) then
                i = l - 1
            endif
            arr(i+1) = a
        enddo

        if (jstack.eq.0) then
            return ! Everything has been sorted
        endif

        ! Pop stack and begin a new round of partitioning
        ir = istack(jstack)
        l = istack(jstack-1)
        jstack = jstack - 2

    else

        ! Choose median of left, center, and right elements as partitioning element a
        k = (l+ir)/2 ! median element -> partitioning element

        ! Also rearrange so that a(l) ≤ a(l+1) ≤ a(ir).
        ! Move partitioning element k next to element l
        temp = arr(k)
        arr(k) = arr(l+1)
        arr(l+1) = temp
        ! write(0,*) 'quicksort: rearranging elements'
        if (arr(l).gt.arr(ir)) then
            temp = arr(l)
            arr(l) = arr(ir)
            arr(ir) = temp
        endif
        if (arr(l+1).gt.arr(ir)) then
            temp = arr(l+1)
            arr(l+1) = arr(ir)
            arr(ir) = temp
        endif
        if (arr(l).gt.arr(l+1)) then
            temp = arr(l)
            arr(l) = arr(l+1)
            arr(l+1) = temp
        endif

        ! Initialize pointers for partitioning
        i = l + 1
        j = ir
        a = arr(l+1)   ! Partitioning element

        ! Put all elements < a to its left, all elements > a to its right
        do ! Beginning of innermost loop
            i = i + 1

            ! Scan up to find element > a
            if (arr(i).lt.a) then
                cycle
            endif

            ! Scan down to find element < a
            do while (arr(j).gt.a)
                j = j - 1
            enddo

            ! if (j.lt.i) then
            if (j.le.i) then
                ! Pointers crossed; exit with partitioning complete
                exit
            endif

            ! Exchange elements
            temp = arr(i)
            arr(i) = arr(j)
            arr(j) = temp
        enddo ! End of innermost loop

        ! Insert partitioning element
        arr(l+1) = arr(j)
        arr(j) = a
        jstack = jstack + 2

        ! Push pointers to larger subarray on stack, process smaller subarray immediately
        if (jstack.gt.NSTACK) then
            stop 'quicksort_real: NSTACK too small'
        endif

        if (ir-i+1.ge.j-l) then
            istack(jstack) = ir
            istack(jstack-1) = i
            ir = j - 1
        else
            istack(jstack) = j - 1
            istack(jstack-1) = l
            l = i
        endif
    endif
enddo

return
end subroutine

!--------------------------------------------------------------------------------------------------!

! ! quicksort.f -*-f90-*-
! ! Author: t-nissie
! ! License: GPLv3
! ! Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
! !!
! recursive subroutine quicksort_int(a, first, last)
!   implicit none
!   integer ::  a(*), x, t
!   integer :: first, last
!   integer :: i, j
!
!   x = a( (first+last) / 2 )
!   i = first
!   j = last
!   do
!      do while (a(i) < x)
!         i=i+1
!      end do
!      do while (x < a(j))
!         j=j-1
!      end do
!      if (i >= j) exit
!      t = a(i);  a(i) = a(j);  a(j) = t
!      i=i+1
!      j=j-1
!   end do
!   if (first < i-1) call quicksort(a, first, i-1)
!   if (j+1 < last)  call quicksort(a, j+1, last)
! end subroutine quicksort_int

!--------------------------------------------------------------------------------------------------!

! recursive subroutine quicksort_real(a, first, last)
!   implicit none
!   real ::  a(*), x, t
!   integer :: first, last
!   integer :: i, j
!
!   x = a( (first+last) / 2 )
!   i = first
!   j = last
!   do
!      do while (a(i) < x)
!         i=i+1
!      end do
!      do while (x < a(j))
!         j=j-1
!      end do
!      if (i >= j) exit
!      t = a(i);  a(i) = a(j);  a(j) = t
!      i=i+1
!      j=j-1
!   end do
!   if (first < i-1) call quicksort(a, first, i-1)
!   if (j+1 < last)  call quicksort(a, j+1, last)
! end subroutine quicksort_real

!--------------------------------------------------------------------------------------------------!
! recursive subroutine quicksort_dp(a, first, last)
!   implicit none
!   double precision ::  a(*), x, t
!   integer :: first, last
!   integer :: i, j
!
!   x = a( (first+last) / 2 )
!   i = first
!   j = last
!   do
!      do while (a(i) < x)
!         i=i+1
!      end do
!      do while (x < a(j))
!         j=j-1
!      end do
!      if (i >= j) exit
!      t = a(i);  a(i) = a(j);  a(j) = t
!      i=i+1
!      j=j-1
!   end do
!   if (first < i-1) call quicksort(a, first, i-1)
!   if (j+1 < last)  call quicksort(a, j+1, last)
! end subroutine quicksort_dp

!--------------------------------------------------------------------------------------------------!

subroutine heapsort_int(ra, n)
!----
! Sorts an array ra(1:n) into ascending numerical order using the Heapsort algorithm. n is input;
! ra is replaced on output by its sorted rearrangement.
!
! Numerical Recipes (Press et al.): 8.3
!----

implicit none

! Arguments
integer :: n
integer :: ra(n)

! Local variables
integer :: i, ir, j, l
integer :: rra

! No need to sort 1-item list...
if (n.lt.2) then
    return
endif

! The index l will be decremented from its initial value down to 1 during the “hiring” (heap
! creation) phase. Once it reaches 1, the index ir will be decremented from its initial value down
! to 1 during the “retirement-and-promotion” (heap selection) phase.
l = n/2 + 1
ir = n

do
    if (l.gt.1) then
        ! Still in hiring phase
        l = l - 1
        rra = ra(l)
    else
        ! In retirement-and-promotion phase
        rra = ra(ir)        ! Clear a space at end of array
        ra(ir) = ra(1)      ! Retire the top of the heap into it
        ir = ir - 1         ! Decrease the size of the corporation
        if (ir.eq.1) then
            ! Done with the last promotion
            ra(1) = rra     ! The least competent worker of all
            return
       endif
   endif

    ! Whether in the hiring phase or promotion phase, we here set up to sift down element rra to
    ! its proper level
    i = l
    j = l + l
    do while (j.le.ir)
        if (j.lt.ir) then
            if (ra(j).lt.ra(j+1)) then
                j = j + 1    ! Compare to the better underling
            endif
        endif
        if (rra.lt.ra(j)) then
            ! Demote rra
            ra(i) = ra(j)
            i = j
            j = j + j
        else
            ! This is rra’s level; set j to terminate the sift-down
            j = ir + 1
        endif
    enddo
    ra(i) = rra             ! Put rra into its slot
enddo

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine heapsort_real(ra, n)
!----
! Sorts an array ra(1:n) into ascending numerical order using the Heapsort algorithm. n is input;
! ra is replaced on output by its sorted rearrangement.
!
! Numerical Recipes (Press et al.): 8.3
!----

implicit none

! Arguments
integer :: n
real :: ra(n)

! Local variables
integer :: i, ir, j, l
real :: rra

! No need to sort 1-item list...
if (n.lt.2) then
    return
endif

! The index l will be decremented from its initial value down to 1 during the “hiring” (heap
! creation) phase. Once it reaches 1, the index ir will be decremented from its initial value down
! to 1 during the “retirement-and-promotion” (heap selection) phase.
l = n/2 + 1
ir = n

do
    if (l.gt.1) then
        ! Still in hiring phase
        l = l - 1
        rra = ra(l)
    else
        ! In retirement-and-promotion phase
        rra = ra(ir)        ! Clear a space at end of array
        ra(ir) = ra(1)      ! Retire the top of the heap into it
        ir = ir - 1         ! Decrease the size of the corporation
        if (ir.eq.1) then
            ! Done with the last promotion
            ra(1) = rra     ! The least competent worker of all
            return
       endif
   endif

    ! Whether in the hiring phase or promotion phase, we here set up to sift down element rra to
    ! its proper level
    i = l
    j = l + l
    do while (j.le.ir)
        if (j.lt.ir) then
            if (ra(j).lt.ra(j+1)) then
                j = j + 1    ! Compare to the better underling
            endif
        endif
        if (rra.lt.ra(j)) then
            ! Demote rra
            ra(i) = ra(j)
            i = j
            j = j + j
        else
            ! This is rra’s level; set j to terminate the sift-down
            j = ir + 1
        endif
    enddo
    ra(i) = rra             ! Put rra into its slot
enddo

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine heapsort_dp(ra, n)
!----
! Sorts an array ra(1:n) into ascending numerical order using the Heapsort algorithm. n is input;
! ra is replaced on output by its sorted rearrangement.
!
! Numerical Recipes (Press et al.): 8.3
!----

implicit none

! Arguments
integer :: n
double precision :: ra(n)

! Local variables
integer :: i, ir, j, l
double precision :: rra

! No need to sort 1-item list...
if (n.lt.2) then
    return
endif

! The index l will be decremented from its initial value down to 1 during the “hiring” (heap
! creation) phase. Once it reaches 1, the index ir will be decremented from its initial value down
! to 1 during the “retirement-and-promotion” (heap selection) phase.
l = n/2 + 1
ir = n

do
    if (l.gt.1) then
        ! Still in hiring phase
        l = l - 1
        rra = ra(l)
    else
        ! In retirement-and-promotion phase
        rra = ra(ir)        ! Clear a space at end of array
        ra(ir) = ra(1)      ! Retire the top of the heap into it
        ir = ir - 1         ! Decrease the size of the corporation
        if (ir.eq.1) then
            ! Done with the last promotion
            ra(1) = rra     ! The least competent worker of all
            return
       endif
   endif

    ! Whether in the hiring phase or promotion phase, we here set up to sift down element rra to
    ! its proper level
    i = l
    j = l + l
    do while (j.le.ir)
        if (j.lt.ir) then
            if (ra(j).lt.ra(j+1)) then
                j = j + 1    ! Compare to the better underling
            endif
        endif
        if (rra.lt.ra(j)) then
            ! Demote rra
            ra(i) = ra(j)
            i = j
            j = j + j
        else
            ! This is rra’s level; set j to terminate the sift-down
            j = ir + 1
        endif
    enddo
    ra(i) = rra             ! Put rra into its slot
enddo

return
end subroutine


end module
