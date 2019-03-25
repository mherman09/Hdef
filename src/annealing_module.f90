!--------------------------------------------------------------------------------------------------!
! Module annealing
!
! Routines to implement the simulated annealing algorithm (Kirkpatrick et al., 1983), used to
! determine a set of NPARAM model parameters that maximize some objective function.
!
! Pseudo code:
! 0) Initialize a model, M(0), and a temperature, T [init]
! 1) For i = 1, 2, ... , N
!     1) Propose a new model: M(i) [propose]
!     2) Calculate the objective function for the proposed model: obj(i) [objective]
!        (Note that the algorithm is extremely sensitive to the choice of this objective function
!        because of the exponential in the probability expression!)
!     3) Accept the proposed model with probability: min((1,exp(obj(i)-obj(i-1))/T)
!     4) Decrease T
!
! Arguments:
!     nparam:       Number of parameters in model array
!     model_best:   Model array containing best model from search
!     *init:        External subroutine to initialize model array
!     *propose:     External subroutine to propose a new model
!     *objective:   External subroutine to compute objective function for model
!     max_it:       Number of iterations to search
!     reset_it:     Number of iterations per temperature reset
!     T_start:      Starting temperature
!     T_min:        Minimum temperature
!     cool:         Cooling factor
!     log_file:     Log file for annealing results
!     saveRejected: Write rejected models to log file
!
! *The user must supply the following external routine names on the command line:
!     - subroutine init(model_init, nparam)
!     - subroutine propose(model_current, model_proposed, nparam)
!     - function objective(model, nparam)
!
! Examples of these "driver" routines are provided after the annealing module.
!
! Annealing subroutines are provided for INTEGER and DOUBLE PRECISION model arrays. These are
! interfaced to the generic subroutine call "anneal". This means that the user simply calls the
! subroutine "anneal" in their program. It also means that they must define the argument variable
! and routine types precisely the same as in either "anneal_int_array" or "anneal_dp_array". If the
! variable types do not match exactly, there will be a compiler error; unfortunately, these error
! messages are totally uninformative, e.g., if I declare reset_it as double precision in Example #1,
! gfortran gives me the error message:
!
!     "Error: There is no specific subroutine for the generic ‘anneal’ at (1)"
!
! As you can see, this is not particularly helpful. But now you know.
!
!--------------------------------------------------------------------------------------------------!

module annealing

interface anneal
    module procedure anneal_int_array
    module procedure anneal_dp_array
end interface

!--------------------------------------------------------------------------------------------------!
contains
!--------------------------------------------------------------------------------------------------!

subroutine anneal_int_array(nparam, &
                            model_best, &
                            init, &
                            propose, &
                            objective, &
                            max_it, &
                            reset_it, &
                            T_start, &
                            T_min, &
                            cool, &
                            log_file, &
                            saveRejected)
!----
! Simulated annealing algorithm with an integer model array
!----

use random, only: ran0, timeseed, iseed

implicit none

! Arguments: model array
integer :: nparam
integer :: model_best(NPARAM)

! Arguments: user-supplied external subroutines
interface
    subroutine init(model_init,n)
        integer :: n
        integer :: model_init(n)
    end subroutine
    subroutine propose(model_current,model_proposed,n)
        integer :: n
        integer :: model_current(n)
        integer :: model_proposed(n)
    end subroutine
    function objective(model,n)
        integer :: n
        integer :: model(n)
        double precision :: objective
    end function
end interface

! Arguments: annealing controls
integer :: max_it
integer :: reset_it
double precision :: T_start
double precision :: T_min
double precision :: cool
character(len=*) ::log_file
logical :: saveRejected

! Local variables
integer :: i
integer :: iparam
double precision :: T
double precision :: obj_current
double precision :: obj_proposed
double precision :: obj_best
double precision :: ptrans
integer :: model_current(NPARAM)
integer :: model_proposed(NPARAM)


! Initialize the random number generator with integer based on system time
iseed = -timeseed()

! Initialize the temperature from the argument list
T = T_start

! Initialize a model and calculate its objective function
call init(model_current,nparam)
obj_current = objective(model_current,nparam)
model_best = model_current
obj_best = obj_current

! Open log file and write first entry
if (log_file.ne.'') then
open(unit=28,file=log_file,status='unknown')
    write(28,*) 'Iteration 0  Temperature ',T,' Objective ',obj_current
    do iparam = 1,nparam
        write(28,*) model_current(iparam)
    enddo
endif

! Run the annealing search
do i = 1,max_it

    ! Propose a new model and calculate its objective function
    call propose(model_current,model_proposed,nparam)
    obj_proposed = objective(model_proposed,nparam)

    ! Save the model if it is the best overall
    if (obj_proposed.gt.obj_best) then
        model_best = model_proposed
        obj_best = obj_proposed
    endif

    ! Calculate the transition probability for current to proposed model (maximizing objective func)
    ptrans = exp((obj_proposed-obj_current)/T)
    ptrans = min(1.0d0,ptrans)

    ! Update the current model if transition is made
    if (ran0(iseed).lt.ptrans) then
        model_current = model_proposed
        obj_current = obj_proposed

        ! Save accepted models to the log file
        if (log_file.ne.'') then
            write(28,*) 'Iteration ',i,' Temperature ',T,' Objective ',obj_current
            do iparam = 1,nparam
                write(28,*) model_current(iparam)
            enddo
        endif
    else
        if (saveRejected.and.log_file.ne.'') then
            ! Save unaccepted models to the log file
            write(28,*) 'Iteration ',i,' Temperature ',T,' Objective ',obj_current,' (rejected)'
            do iparam = 1,nparam
                write(28,*) model_proposed(iparam)
            enddo
        endif
    endif

    ! Update temperature
    T = T*cool

    ! Make sure temperature stays above T_min
    if (T.lt.T_min) then
        T = T_min
    endif

    ! Reset temperature every reset_it iterations
    if (reset_it.gt.0.and.mod(i,reset_it).eq.0) then
        T = T_start
    endif
enddo

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine anneal_dp_array( nparam, &
                            model_best, &
                            init, &
                            propose, &
                            objective, &
                            max_it, &
                            reset_it, &
                            T_start, &
                            T_min, &
                            cool, &
                            log_file, &
                            saveRejected)
!----
! Simulated annealing algorithm with a double precision model array
!----

use random, only: ran0, timeseed, iseed

implicit none

! Arguments: model array
integer :: nparam
double precision :: model_best(NPARAM)

! Arguments: user-supplied external subroutines
interface
    subroutine init(model_init,n)
        integer :: n
        double precision :: model_init(n)
    end subroutine
    subroutine propose(model_current,model_proposed,n)
        integer :: n
        double precision :: model_current(n)
        double precision :: model_proposed(n)
    end subroutine
    function objective(model,n)
        integer :: n
        double precision :: model(n)
        double precision :: objective
    end function
end interface

! Arguments: annealing controls
integer :: max_it
integer :: reset_it
double precision :: T_start
double precision :: T_min
double precision :: cool
character(len=*) ::log_file
logical :: saveRejected

! Local variables
integer :: i
integer :: iparam
double precision :: T
double precision :: obj_current
double precision :: obj_proposed
double precision :: obj_best
double precision :: ptrans
double precision :: model_current(NPARAM)
double precision :: model_proposed(NPARAM)


! Initialize the random number generator with integer based on system time
iseed = -timeseed()

! Initialize the temperature from the argument list
T = T_start

! Initialize a model and calculate its objective function
call init(model_current,nparam)
obj_current = objective(model_current,nparam)
model_best = model_current
obj_best = obj_current

! Open log file and write first entry
if (log_file.ne.'') then
open(unit=28,file=log_file,status='unknown')
    write(28,*) 'Iteration 0  Temperature ',T,' Objective ',obj_current
    do iparam = 1,nparam
        write(28,*) model_current(iparam)
    enddo
endif

! Run the annealing search
do i = 1,max_it

    ! Propose a new model and calculate its objective function
    call propose(model_current,model_proposed,nparam)
    obj_proposed = objective(model_proposed,nparam)

    ! Save the model if it is the best overall
    if (obj_proposed.gt.obj_best) then
        model_best = model_proposed
        obj_best = obj_proposed
    endif

    ! Calculate the transition probability for current to proposed model (maximizing objective func)
    ptrans = exp((obj_proposed-obj_current)/T)
    ptrans = min(1.0d0,ptrans)

    ! Update the current model if transition is made
    if (ran0(iseed).lt.ptrans) then
        model_current = model_proposed
        obj_current = obj_proposed

        ! Save accepted models to the log file
        if (log_file.ne.'') then
            write(28,*) 'Iteration ',i,' Temperature ',T,' Objective ',obj_current
            do iparam = 1,nparam
                write(28,*) model_current(iparam)
            enddo
        endif
    else
        if (saveRejected.and.log_file.ne.'') then
            ! Save unaccepted models to the log file
            write(28,*) 'Iteration ',i,' Temperature ',T,' Objective ',obj_current,' (rejected)'
            do iparam = 1,nparam
                write(28,*) model_proposed(iparam)
            enddo
        endif
    endif

    ! Update temperature
    T = T*cool

    ! Make sure temperature stays above T_min
    if (T.lt.T_min) then
        T = T_min
    endif

    ! Reset temperature every reset_it iterations
    if (reset_it.gt.0.and.mod(i,reset_it).eq.0) then
        T = T_start
    endif
enddo

return
end subroutine

end module

!==================================================================================================!

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!------------------ EXAMPLES OF EXTERNAL DRIVER ROUTINES FOR SIMULATED ANNEALING ------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

!--------------------------------------------------------------------------------------------------!
! EXAMPLE #1: FIND AN INTEGER POINT IN TWO DIMENSIONS
! - Objective function: obj = -0.5*((x-x0)**2/x0+(y-y0)**2/y0)
! - Initialization: model = 0
! - Proposal: model_proposed = model_current +/- 1
!--------------------------------------------------------------------------------------------------!
subroutine driver1()

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
end interface

integer :: n, maxit, resetit                                        ! These variable types must be exactly
integer, allocatable :: model_best(:)                               ! the same as in the corresponding
double precision :: tstart, tmin, cool                              ! version of anneal.
character(len=32) :: logfile
logical :: saveRejected

n = 2
allocate(model_best(n))
maxit = 200
resetit = -1
tstart = 1.0d0
tmin = -1.0d0
cool = 0.99d0
logfile = 'driver1.log'
saveRejected = .true.

! Call anneal with specific driver routines
call anneal(n, model_best, &
            driver1_init, driver1_propose, driver1_objective, &
            maxit, resetit, tstart, tmin, cool, &
            logfile, saveRejected)

return
end subroutine

subroutine driver1_init(model,n)
implicit none
integer :: n                                                        ! As in the driver routine, these types
integer :: model(n)                                                 ! must be the same as in anneal
model = 0
return
end subroutine

subroutine driver1_propose(model_in,model_out,n)
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
implicit none
integer :: n
integer :: model(n)
double precision :: driver1_objective
driver1_objective = 0.0d0
driver1_objective = driver1_objective + dble((model(1)-10)**2)/dble(10)
driver1_objective = driver1_objective + dble((model(2)-6)**2)/dble(6)
driver1_objective = -0.5d0*driver1_objective
return
end function


!--------------------------------------------------------------------------------------------------!
! EXAMPLE #2: FIND A DOUBLE PRECISION POINT IN TWO DIMENSIONS
! - Objective function: obj = -0.5*((x-x0)**2/x0+(y-y0)**2/y0)
! - Initialization: model = 0
! - Proposal: model_proposed = model_current + gaussian random number
!--------------------------------------------------------------------------------------------------!
subroutine driver2()

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
end interface

integer :: n, maxit, resetit
double precision, allocatable :: model_best(:)
double precision :: tstart, tmin, cool
character(len=32) :: logfile
logical :: saveRejected

n = 2
allocate(model_best(n))
maxit = 200
resetit = -1
tstart = 1.0d0
tmin = -1.0d0
cool = 0.99d0
logfile = 'driver2.log'
saveRejected = .true.

! Call anneal with specific driver routines
call anneal(n, model_best, &
            driver2_init, driver2_propose, driver2_objective, &
            maxit, resetit, tstart, tmin, cool, &
            logfile, saveRejected)

return
end subroutine

subroutine driver2_init(model,n)
implicit none
integer :: n
double precision :: model(n)
model = 0.0d0
return
end subroutine

subroutine driver2_propose(model_in,model_out,n)
use random, only: iseed, ran0, r8_normal_ab
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
implicit none
integer :: n
double precision :: model(n)
double precision :: driver2_objective
driver2_objective = 0.0d0
driver2_objective = driver2_objective + (model(1)-10.0d0)**2/10.0d0
driver2_objective = driver2_objective + (model(2)-6.0d0)**2/6.0d0
driver2_objective = -0.5d0*driver2_objective
return
end function


!--------------------------------------------------------------------------------------------------!
! EXAMPLE #3: FIND A DOUBLE PRECISION POINT IN TWO DIMENSIONS, WITH DRIVERS IN A MODULE
! - Objective function: obj = -0.5*((x-x0)**2/x0+(y-y0)**2/y0)
! - Initialization: model = 0
! - Proposal: model_proposed = model_current + gaussian random number
!--------------------------------------------------------------------------------------------------!
module driver3_module

contains

subroutine driver3()

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
double precision, allocatable :: model_best(:)
double precision :: tstart, tmin, cool
character(len=32) :: logfile
logical :: saveRejected

n = 2
allocate(model_best(n))
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
            logfile, saveRejected)

return
end subroutine

subroutine driver3_init(model,n)
implicit none
integer :: n
double precision :: model(n)
model = 0.0d0
return
end subroutine

subroutine driver3_propose(model_in,model_out,n)
use random, only: iseed, ran0, r8_normal_ab
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
driver3_objective = driver3_objective + (model(1)-10.0d0)**2/10.0d0
driver3_objective = driver3_objective + (model(2)-6.0d0)**2/6.0d0
driver3_objective = -0.5d0*driver3_objective
return
end function

end module
