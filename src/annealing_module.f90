!--------------------------------------------------------------------------------------------------!
! Module: annealing
!
! Routines to implement the simulated annealing algorithm (Kirkpatrick et al., 1983), used to
! determine a set of NPARAM model parameters that maximize some objective function.
!
! Simulated annealing pseudo code:
! 0) Initialize a model, M(0), and a temperature, T [init]
! 1) For i = 1, 2, ... , N
!     1) Propose a new model: M(i) [propose]
!     2) Calculate the objective function for the proposed model: obj(i) [objective]
!        (Note that the algorithm is extremely sensitive to the choice of this objective function
!        because of the exponential in the probability expression!)
!     3) Accept the proposed model with probability: min((1,exp(obj(i)-obj(i-1))/T)
!     4) Decrease T
!
! Subroutine arguments:
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
!     *write_log:   Write  annealing results to log file
!
! *The user must supply the following external routine names in the argument list:
!     - subroutine init(model_init, nparam)
!     - subroutine propose(model_current, model_proposed, nparam)
!     - function objective(model, nparam)
!     - subroutine write_log()
!
! Examples of these "driver" routines are provided in the unit test for this module
! (unit_test_annealing.f90).
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

public :: anneal

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
                            write_log)
!----
! Simulated annealing algorithm with an integer model array
!----

use io, only: verbosity, progress_indicator, stdout
use random, only: r8_uniform_01, timeseed, iseed

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
    subroutine write_log(it,temp,obj,model_current,model_proposed,n,string)
        integer :: it, n
        double precision :: temp, obj
        integer :: model_current(n)
        integer :: model_proposed(n)
        character(len=*) :: string
    end subroutine
end interface

! Arguments: annealing controls
integer :: max_it
integer :: reset_it
double precision :: T_start
double precision :: T_min
double precision :: cool

! Local variables
integer :: i
integer :: ierr
double precision :: T
double precision :: obj_current
double precision :: obj_proposed
double precision :: obj_best
double precision :: ptrans
double precision :: ran_uniform
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
call write_log(0,T,obj_current,model_current,model_current,nparam,'init')

! Run the annealing search
do i = 1,max_it

    if (verbosity.eq.1) then
        call progress_indicator(i,max_it,'anneal',ierr)
    elseif (verbosity.ge.2) then
        write(stdout,*) 'anneal: starting iteration ',i,' of ',max_it
    endif

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
    ran_uniform = r8_uniform_01(iseed)
    if (ran_uniform.lt.ptrans) then
        model_current = model_proposed
        obj_current = obj_proposed
    endif

    ! Save results in log file
    call write_log(i,T,obj_current,model_current,model_proposed,nparam,'append')

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

! Close log file
call write_log(i,T,obj_current,model_current,model_proposed,nparam,'close')

return
end subroutine anneal_int_array

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
                            write_log)
!----
! Simulated annealing algorithm with a double precision model array
!----

use io, only: verbosity, progress_indicator, stdout
use random, only: r8_uniform_01, timeseed, iseed

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
    subroutine write_log(it,temp,obj,model_current,model_proposed,n,string)
        integer :: it, n
        double precision :: temp, obj
        double precision :: model_current(n)
        double precision :: model_proposed(n)
        character(len=*) :: string
    end subroutine
end interface

! Arguments: annealing controls
integer :: max_it
integer :: reset_it
double precision :: T_start
double precision :: T_min
double precision :: cool

! Local variables
integer :: i
integer :: ierr
double precision :: T
double precision :: obj_current
double precision :: obj_proposed
double precision :: obj_best
double precision :: ptrans
double precision :: model_current(NPARAM)
double precision :: model_proposed(NPARAM)
double precision :: ran_uniform


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
call write_log(0,T,obj_current,model_current,model_current,nparam,'init')

! Run the annealing search
do i = 1,max_it

    if (verbosity.eq.1) then
        call progress_indicator(i,max_it,'anneal',ierr)
    elseif (verbosity.ge.2) then
        write(stdout,*) 'anneal: starting iteration ',i,' of ',max_it
    endif

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
    ran_uniform = r8_uniform_01(iseed)
    if (ran_uniform.lt.ptrans) then
        model_current = model_proposed
        obj_current = obj_proposed
    endif

    ! Save results in log file
    call write_log(i,T,obj_current,model_current,model_proposed,nparam,'append')

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

! Close log file
call write_log(i,T,obj_current,model_current,model_proposed,nparam,'close')


return
end subroutine anneal_dp_array

end module

!==================================================================================================!
