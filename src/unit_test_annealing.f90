!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!------------------ EXAMPLES OF EXTERNAL DRIVER ROUTINES FOR SIMULATED ANNEALING ------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

!--------------------------------------------------------------------------------------------------!
! EXAMPLE #1: FIND AN INTEGER POINT IN TWO DIMENSIONS
! - Objective function: obj = -0.5*((x-x0)**2/sx+(y-y0)**2/sy)
! - Initialization: model = 0
! - Proposal: model_proposed = model_current +/- 1
!--------------------------------------------------------------------------------------------------!
subroutine driver1(model_best)

use annealing                                                       ! Import the annealing module
implicit none

interface                                                           ! The driver subroutines and functions
    subroutine driver1_init(model,n)                                ! must be declared with an interface
        integer :: n                                                ! statement so they can be matched with
        integer :: model(n)                                         ! the generic routines called in anneal().
    end subroutine                                                  ! This is no longer the case if the
    subroutine driver1_propose(model_current,model_proposed,n)      ! driver routines are part of a module
        integer :: n                                                ! (see Example #3).
        integer :: model_current(n)
        integer :: model_proposed(n)
    end subroutine
    function driver1_objective(model,n)
        integer :: n
        integer :: model(n)
        double precision :: driver1_objective
    end function
    subroutine driver1_log(it,temp,obj,model_current,model_proposed,n,isAccepted,string)
        integer :: it, n
        double precision :: temp, obj
        integer :: model_current(n)
        integer :: model_proposed(n)
        logical :: isAccepted
        character(len=*) :: string
    end subroutine
end interface

integer :: n, maxit, resetit                                        ! These variable types must be exactly
integer :: model_best(2)                                            ! the same as in the corresponding
double precision :: tstart, tmin, cool                              ! version of anneal.
logical :: saveRejected

n = 2
model_best = 0
maxit = 200
resetit = -1
tstart = 1.0d0
tmin = -1.0d0
cool = 0.99d0
saveRejected = .true.

! Call anneal with specific driver routines
call anneal(n, model_best, &
            driver1_init, driver1_propose, driver1_objective, &
            maxit, resetit, tstart, tmin, cool, &
            driver1_log)

return
end subroutine

subroutine driver1_init(model,n)
!----
! Start model at origin, initialize seed to fixed value
!----
use random, only: iseed
implicit none
integer :: n                                                        ! As in the driver routine, these types
integer :: model(n)                                                 ! must be the same as in anneal
iseed = -52918
model = 0
return
end subroutine

subroutine driver1_propose(model_in,model_out,n)
!----
! Equal odds of adding or subtracting one from each model dimension
!----
use random, only: iseed, ran0
implicit none
integer :: n
integer :: model_in(n)
integer :: model_out(n)
integer :: i
do i = 1,n
    if (ran0(iseed).lt.0.50) then
        model_out(i) = model_in(i) + 1
    else
        model_out(i) = model_in(i) - 1
    endif
enddo
return
end subroutine

function driver1_objective(model,n)
!----
! Aiming for coordinate (10,6), with Gaussian probability distribution
!----
implicit none
integer :: n
integer :: model(n)
double precision :: driver1_objective
driver1_objective = 0.0d0
driver1_objective = driver1_objective + dble((model(1)-10)**2)/dble(20)
driver1_objective = driver1_objective + dble((model(2)-6)**2)/dble(8)
driver1_objective = -0.5d0*driver1_objective
return
end function

subroutine driver1_log(it,temp,obj,model_current,model_proposed,n,isAccepted,string)
integer :: it, n
double precision :: temp, obj
integer :: model_current(n)
integer :: model_proposed(n)
logical :: isAccepted
character(len=*) :: string
character(len=32) :: logfile
logfile = 'driver1.log'
if (string.eq.'init') then
    ! Open the log file
    open(unit=29,file=logfile,status='unknown')
    ! Write locked/unlocked, fault slip results to log file
    write(29,*) 'Iteration ',it,' Temperature ',temp,' Objective ',obj
    do i = 1,n
        write(29,*) model_current(i)
    enddo
elseif (string.eq.'append') then
    ! Write locked/unlocked, fault slip, old model results to log file
    write(29,*) 'Iteration ',it,' Temperature ',temp,' Objective ',obj
    do i = 1,n
        write(29,*) model_current(i),model_proposed(i)
    enddo
elseif (string.eq.'close') then
    ! Close the log file
    close(29)
endif
end subroutine

!--------------------------------------------------------------------------------------------------!
! EXAMPLE #2: FIND A DOUBLE PRECISION POINT IN TWO DIMENSIONS
! - Objective function: obj = -0.5*((x-x0)**2/sx+(y-y0)**2/sy)
! - Initialization: model = 0
! - Proposal: model_proposed = model_current + gaussian random number
!--------------------------------------------------------------------------------------------------!
subroutine driver2(model_best)

use annealing
implicit none

interface
    subroutine driver2_init(model,n)
        integer :: n
        double precision :: model(n)
    end subroutine
    subroutine driver2_propose(model_current,model_proposed,n)
        integer :: n
        double precision :: model_current(n)
        double precision :: model_proposed(n)
    end subroutine
    function driver2_objective(model,n)
        integer :: n
        double precision :: model(n)
        double precision :: driver2_objective
    end function
    subroutine driver2_log(it,temp,obj,model_current,model_proposed,n,isAccepted,string)
        integer :: it, n
        double precision :: temp, obj
        double precision :: model_current(n)
        double precision :: model_proposed(n)
        logical :: isAccepted
        character(len=*) :: string
    end subroutine
end interface

integer :: n, maxit, resetit
double precision :: model_best(2)
double precision :: tstart, tmin, cool
logical :: saveRejected

n = 2
model_best = 0.0d0
maxit = 200
resetit = -1
tstart = 1.0d0
tmin = -1.0d0
cool = 0.99d0
saveRejected = .true.

! Call anneal with specific driver routines
call anneal(n, model_best, &
            driver2_init, driver2_propose, driver2_objective, &
            maxit, resetit, tstart, tmin, cool, &
            driver2_log)

return
end subroutine

subroutine driver2_init(model,n)
!----
! Start model at origin, initialize seed to fixed value
!----
use random, only: iseed
implicit none
integer :: n
double precision :: model(n)
iseed = -829321
model = 0.0d0
return
end subroutine

subroutine driver2_propose(model_in,model_out,n)
!----
! Gaussian probability distribution with a standard deviation of 1 for each model dimension
!----
use random, only: iseed, r8_normal_ab
implicit none
integer :: n
double precision :: model_in(n)
double precision :: model_out(n)
integer :: i
do i = 1,n
    model_out(i) = model_in(i) + r8_normal_ab(0.0d0,1.0d0,iseed)
enddo
return
end subroutine

function driver2_objective(model,n)
!----
! Aiming for coordinate (10,6), with Gaussian probability distribution
!----
implicit none
integer :: n
double precision :: model(n)
double precision :: driver2_objective
driver2_objective = 0.0d0
driver2_objective = driver2_objective + (model(1)-10.0d0)**2/20.0d0
driver2_objective = driver2_objective + (model(2)-6.0d0)**2/8.0d0
driver2_objective = -0.5d0*driver2_objective
return
end function

subroutine driver2_log(it,temp,obj,model_current,model_proposed,n,isAccepted,string)
integer :: it, n
double precision :: temp, obj
double precision :: model_current(n)
double precision :: model_proposed(n)
logical :: isAccepted
character(len=*) :: string
character(len=32) :: logfile
logfile = 'driver2.log'
if (string.eq.'init') then
    ! Open the log file
    open(unit=29,file=logfile,status='unknown')
    ! Write locked/unlocked, fault slip results to log file
    write(29,*) 'Iteration ',it,' Temperature ',temp,' Objective ',obj
    do i = 1,n
        write(29,*) model_current(i)
    enddo
elseif (string.eq.'append') then
    ! Write locked/unlocked, fault slip, old model results to log file
    write(29,*) 'Iteration ',it,' Temperature ',temp,' Objective ',obj
    do i = 1,n
        write(29,*) model_current(i),model_proposed(i)
    enddo
elseif (string.eq.'close') then
    ! Close the log file
    close(29)
endif
end subroutine


!--------------------------------------------------------------------------------------------------!
! EXAMPLE #3: FIND A DOUBLE PRECISION POINT IN TWO DIMENSIONS, WITH DRIVERS IN A MODULE
! - Objective function: obj = -0.5*((x-x0)**2/x0+(y-y0)**2/y0)
! - Initialization: model = 0
! - Proposal: model_proposed = model_current + gaussian random number
!--------------------------------------------------------------------------------------------------!
module driver3_module

character(len=32) :: logfile

contains

subroutine driver3(model_best)

use annealing, only: anneal
implicit none

! interface                                                         ! When the driver routines are part
!     subroutine driver3_init(model,n)                              ! of a module, then they do not need
!         integer :: n                                              ! an explicit interface statement
!         double precision :: model(n)                              ! because they are no longer (implicitly)
!     end subroutine                                                ! external routines. (In fact, this will
!     subroutine driver3_propose(model_current,model_proposed,n)    ! not compile when these lines are
!         integer :: n                                              ! uncommented.)
!         double precision :: model_current(n)
!         double precision :: model_proposed(n)
!     end subroutine
!     function driver3_objective(model,n)
!         integer :: n
!         double precision :: model(n)
!         double precision :: driver3_objective
!     end function
! end interface

integer :: n, maxit, resetit
double precision :: model_best(2)
double precision :: tstart, tmin, cool
logical :: saveRejected

n = 2
model_best = 0.0d0
maxit = 200
resetit = -1
tstart = 1.0d0
tmin = -1.0d0
cool = 0.99d0
logfile = 'driver3.log'
saveRejected = .true.

! Call anneal with specific driver routines
call anneal(n, model_best, &
            driver3_init, driver3_propose, driver3_objective, &
            maxit, resetit, tstart, tmin, cool, &
            driver3_log)

return
end subroutine

subroutine driver3_init(model,n)
use random, only: iseed
implicit none
integer :: n
double precision :: model(n)
iseed = -829321
model = 0.0d0
return
end subroutine

subroutine driver3_propose(model_in,model_out,n)
use random, only: iseed, r8_normal_ab
implicit none
integer :: n
double precision :: model_in(n)
double precision :: model_out(n)
integer :: i
do i = 1,n
    model_out(i) = model_in(i) + r8_normal_ab(0.0d0,1.0d0,iseed)
enddo
return
end subroutine

function driver3_objective(model,n)
implicit none
integer :: n
double precision :: model(n)
double precision :: driver3_objective
driver3_objective = 0.0d0
driver3_objective = driver3_objective + (model(1)-10.0d0)**2/20.0d0
driver3_objective = driver3_objective + (model(2)-6.0d0)**2/8.0d0
driver3_objective = -0.5d0*driver3_objective
return
end function

subroutine driver3_log(it,temp,obj,model_current,model_proposed,n,isAccepted,string)
integer :: it, n
double precision :: temp, obj
double precision :: model_current(n)
double precision :: model_proposed(n)
logical :: isAccepted
character(len=*) :: string
if (string.eq.'init') then
    ! Open the log file
    open(unit=29,file=logfile,status='unknown')
    ! Write locked/unlocked, fault slip results to log file
    write(29,*) 'Iteration ',it,' Temperature ',temp,' Objective ',obj
    do i = 1,n
        write(29,*) model_current(i)
    enddo
elseif (string.eq.'append') then
    ! Write locked/unlocked, fault slip, old model results to log file
    write(29,*) 'Iteration ',it,' Temperature ',temp,' Objective ',obj
    do i = 1,n
        write(29,*) model_current(i),model_proposed(i)
    enddo
elseif (string.eq.'close') then
    ! Close the log file
    close(29)
endif
end subroutine

end module


!==================================================================================================!
!==================================================================================================!
!==================================================================================================!
!==================================================================================================!
!==================================================================================================!


program main

use test, only: test_value
use io, only: stdout
use random

use annealing
use driver3_module, only: driver3

implicit none

integer :: j, ios, it, model_best_int(2)
character(len=512) :: line
character(len=1) :: ch
double precision :: temp, obj, model_best_dp(2)


! driver1 (anneal_int_array)
call driver1(model_best_int)
open(unit=11,file='driver1.log',status='old')
do
    read(11,'(A)',iostat=ios) line
    if (ios.ne.0) then
        exit
    endif
    j = index(line,'Iteration')
    if (j.ne.0) then
        read(line,*) ch,it,ch,temp,ch,obj
    endif
enddo
call test_value(model_best_int(1),10,'driver1 (anneal_int_array): model_best_int(1)')
call test_value(model_best_int(2),6,'driver1 (anneal_int_array): model_best_int(2)')
call test_value(it,200,'driver1 (anneal_int_array): maxit')
call test_value(temp,0.13533300490703207d0,'driver1 (anneal_int_array): final temp')
call test_value(obj,-0.34999999999999998d0,'driver1 (anneal_int_array): final obj')
close(11,status='delete')


! driver2 (anneal_dp_array)
call driver2(model_best_dp)
open(unit=12,file='driver2.log',status='old')
do
    read(12,'(A)',iostat=ios) line
    if (ios.ne.0) then
        exit
    endif
    j = index(line,'Iteration')
    if (j.ne.0) then
        read(line,*) ch,it,ch,temp,ch,obj
    endif
enddo
call test_value(model_best_dp(1),9.7829045374811461d0,'driver2 (anneal_dp_array): model_best_dp(1)')
call test_value(model_best_dp(2),6.0216229846947424d0,'driver2 (anneal_dp_array): model_best_dp(2)')
call test_value(it,200,'driver2 (anneal_dp_array): maxit')
call test_value(temp,0.13533300490703207d0,'driver2 (anneal_dp_array): final temp')
call test_value(obj,-5.8331140982675264d-002,'driver2 (anneal_dp_array): final obj')
close(12,status='delete')


! driver3 (anneal_dp_array; module)
call driver3(model_best_dp)
open(unit=13,file='driver3.log',status='old')
do
    read(13,'(A)',iostat=ios) line
    if (ios.ne.0) then
        exit
    endif
    j = index(line,'Iteration')
    if (j.ne.0) then
        read(line,*) ch,it,ch,temp,ch,obj
    endif
enddo
call test_value(model_best_dp(1),9.7829045374811461d0,'driver3 (anneal_dp_array): model_best_dp(1)')
call test_value(model_best_dp(2),6.0216229846947424d0,'driver3 (anneal_dp_array): model_best_dp(2)')
call test_value(it,200,'driver3 (anneal_dp_array): maxit')
call test_value(temp,0.13533300490703207d0,'driver3 (anneal_dp_array): final temp')
call test_value(obj,-5.8331140982675264d-002,'driver3 (anneal_dp_array): final obj')
close(13,status='delete')


write(stdout,*) 'annealing_module unit test passed'
end
