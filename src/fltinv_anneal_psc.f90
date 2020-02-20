subroutine invert_anneal_psc()

use annealing, only: anneal

use fltinv, only: fault, &
                  max_iteration, &
                  reset_iteration, &
                  temp_start, &
                  temp_minimum, &
                  cooling_factor, &
                  fault_slip

implicit none

! Interface to driver subroutines
interface
    subroutine anneal_psc_init(model,n)
        integer :: n
        integer :: model(n)
    end subroutine
    subroutine anneal_psc_propose(model_current,model_proposed,n)
        integer :: n
        integer :: model_current(n)
        integer :: model_proposed(n)
    end subroutine
    function anneal_psc_objective(model,n)
        integer :: n
        integer :: model(n)
        double precision :: anneal_psc_objective
    end function
    subroutine anneal_psc_log(it,temp,obj,model_current,model_proposed,n,isAccepted,string)
        integer :: it, n
        double precision :: temp, obj
        integer :: model_current(n)
        integer :: model_proposed(n)
        logical :: isAccepted
        character(len=*) :: string
    end subroutine
end interface

! Local variables
integer :: ierr, nflt
integer, allocatable :: locked(:)
logical :: saveRejected


saveRejected = .true.

nflt = fault%nrows
allocate(locked(nflt),stat=ierr)
if (ierr.ne.0) then
    call usage('invert_anneal_psc: error allocating memory to locked array (usage:none)')
endif

! Call anneal routine with specific driver routines anneal_psc_init, anneal_psc_propose, and
! anneal_psc_objective
call anneal(nflt, &
            locked, &
            anneal_psc_init, &
            anneal_psc_propose, &
            anneal_psc_objective, &
            max_iteration, &
            reset_iteration, &
            temp_start, &
            temp_minimum, &
            cooling_factor, &
            anneal_psc_log)

! Save results for printing
fault_slip = 0.0d0

! Calculate slip in each fault
call psc_slip(locked,nflt,fault_slip)


return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine anneal_psc_init(model,n)
!----
! Initialize the annealing with pseudo-coupling variables:
!     - iseed: random number seed
!     - model array: 0=unlocked, 1=locked
!     - Maximum size Greens function matrix for pseudo-coupling inversion
!----

use io, only: stderr, stdout, verbosity, fileExists, line_count
use random, only: iseed, timeseed, r8_uniform_01

use fltinv, only: fault, &
                  slip_constraint, &
                  anneal_init_mode, &
                  anneal_init_file, &
                  anneal_seed, &
                  min_flip, &
                  max_flip

implicit none

! Arguments
integer, intent(in) :: n
integer :: model(n)

! Local variables
integer :: i, ios, nflt, n_always_unlocked
double precision :: plocked


if (verbosity.ge.2) then
    write(stdout,*) 'anneal_psc_init: starting'
endif

! Initialize the random number generator seed
if (anneal_seed.eq.0) then
    iseed = -timeseed()
else
    iseed = -abs(anneal_seed)
endif

! Check array dimensions
nflt = fault%nrows
if (nflt.ne.n) then
    call usage('anneal_psc_init: input n not equal to nflt (usage:none)')
endif


! Check flipping parameters
if (min_flip.lt.1) then
    min_flip = 1
endif
n_always_unlocked = 0
do i = 1,nflt
    if (abs(slip_constraint%array(i,1)).gt.99998.0d0 .or. &
        abs(slip_constraint%array(i,2)).gt.99998.0d0) then
        n_always_unlocked = n_always_unlocked + 1
    endif
enddo
if (max_flip.gt.nflt-n_always_unlocked) then
    max_flip = nflt-n_always_unlocked
endif
! if (pl2u.le.1.0d-3) then
!     pl2u = 1.0d-3
! endif
! if (pu2l.le.1.0d-3) then
!     pu2l = 1.0d-3
! endif


! Initialize the model array values
model = 0
if (anneal_init_mode.eq.'locked') then
    ! Set the initial solution to all locked
    model = 1

elseif (anneal_init_mode.eq.'unlocked') then
    ! Set the initial solution to all unlocked
    model = 0

elseif (anneal_init_mode(1:4).eq.'rand') then
    ! Set the initial solution to randomly locked, with probability set by the last characters
    ! of the anneal_init_mode variable, e.g., rand0.25 sets probability locked = 0.25
    read(anneal_init_mode(5:len_trim(anneal_init_mode)),*,iostat=ios) plocked
    if (ios.ne.0) then
        write(stderr,*) 'anneal_psc_init: could not read initial plocked, setting to 0.5'
        plocked = 0.5d0
    endif
    do i = 1,nflt
        if (r8_uniform_01(iseed).lt.plocked) then
            model(i) = 1
        endif
    enddo

elseif (anneal_init_mode.eq.'user') then
    ! Read the initial solution from a file (slip rake)
    if (anneal_init_file.eq.'none') then
        call usage('anneal_psc_init: no initialization file defined for anneal_init_mode=user '//&
                   '(usage:anneal)')
    elseif (.not.fileExists(anneal_init_file)) then
        call usage('anneal_psc_init: no anneal_init_file found named "'// &
                   trim(anneal_init_file)//'" (usage:anneal)')
    endif
    if (line_count(anneal_init_file).ne.nflt) then
        call usage('anneal_psc_init: number of lines in anneal_init_file must be equal to nflt '//&
                   '(usage:anneal)')
    endif
    open(unit=29,file=anneal_init_file,status='old')
    do i = 1,nflt
        read(29,*,iostat=ios) model(i)
        if (ios.ne.0) then
            call usage('anneal_psc_init: error reading anneal init file (usage:anneal)')
        endif
    enddo
    close(29)

else
    write(stderr,*) 'anneal_psc_init: no initialization mode named "'//trim(anneal_init_mode)//'"'
    write(stderr,*) 'Options for annealing with pseudo-coupling initialization:'
    write(stderr,*) '    locked'
    write(stderr,*) '    unlocked'
    write(stderr,*) '    rand'
    call usage(     '    user (usage:anneal)')
endif

if (verbosity.ge.2) then
    write(stdout,*) 'anneal_psc_init: finished'
endif

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine anneal_psc_propose(model_in,model_out,n)

use io, only: stdout, stderr, verbosity
use random, only: iseed, r8_uniform_01

use fltinv, only: fault, &
                  slip_constraint, &
                  min_flip, &
                  max_flip

implicit none

! Arguments
integer, intent(in) :: n
integer :: model_in(n), model_out(n)

! Local variables
integer :: i, j, k, nflt, iflip, nflip, rand_fault_list(n)


if (verbosity.ge.3) then
    write(stdout,*) 'anneal_psc_propose: starting'
endif

nflt = fault%nrows
if (nflt.ne.n) then
    call usage('anneal_psc_propose: input n not equal to nflt (usage:none)')
endif


! Randomize fault list by switching each fault with a random fault
! Initialize random fault list array in order
do i = 1,nflt
    rand_fault_list(i) = i
enddo

do i = 1,nflt
    ! Sample another fault index with uniform distribution
    j = int(r8_uniform_01(iseed)*dble(nflt))+1
    if (j.lt.1) then
        j = 1
    elseif (j.gt.nflt) then
        j = nflt
    endif

    ! Switch faults i<->j
    k = rand_fault_list(j)
    rand_fault_list(j) = rand_fault_list(i)
    rand_fault_list(i) = k
enddo


! Initialize the proposed model as the current model
model_out = model_in


! Flip between min_flip and max_flip faults from locked to unlocked or vice versa

! Uniform probability of nflip (number of faults to flip) from min_flip to max_flip
nflip = int(r8_uniform_01(iseed)*dble(max_flip-min_flip+1))+min_flip
if (nflip.gt.max_flip) then
    nflip = max_flip
elseif (nflip.lt.min_flip) then
    nflip = min_flip
endif

! Flip the first nflip faults in the randomized list (that are not always unlocked)
iflip = 0
do i = 1,nflt

    ! Check if the fault is locked/unlocked and should be flipped
    if (model_out(rand_fault_list(i)).eq.1) then
        ! Flip this fault! Locked -> unlocked
        model_out(rand_fault_list(i)) = 0

    elseif (model_out(rand_fault_list(i)).eq.0) then
        ! Flip this fault! Unlocked -> locked
        model_out(rand_fault_list(i)) = 1

    else
        write(stderr,*) 'anneal_psc_propose: invalid locking state ',model_out(rand_fault_list(i))
        call usage('Valid states are 0=unlocked or 1=locked (usage:none)')
    endif

    ! Only count fault as flipped if it is not always unlocked
    if (abs(slip_constraint%array(rand_fault_list(i),1)).lt.99998.0d0 .or. &
                              abs(slip_constraint%array(rand_fault_list(i),2)).lt.99998.0d0) then
        iflip = iflip + 1
        if (iflip.ge.nflip) then
            exit
        endif
    endif
enddo

! write(0,*) 'locked? flipped nflip=',nflip
! do i = 1,nflt
!     write(0,*) i,model_in(i),model_out(i)
! enddo

if (verbosity.ge.3) then
    write(stdout,*) 'anneal_psc_propose: finished'
endif

return
end subroutine

!--------------------------------------------------------------------------------------------------!

function anneal_psc_objective(model,n)

use io, only: stdout, verbosity
use misfit, only: misfit_chi2

use fltinv, only: fault, &
                  displacement, &
                  disp_components, &
                  los, &
                  cov_matrix, &
                  isCovMatrixDiagonal, &
                  gf_disp, &
                  gf_los
implicit none

! Arguments
integer, intent(in) :: n
integer :: model(n)
double precision ::  anneal_psc_objective

! Local variables
integer :: i, ierr, iflt, idsp, icmp, nflt, ndsp, ndsp_dof, nlos, nobs
double precision :: slip(n,2)
double precision, allocatable :: pre(:), obs(:)


if (verbosity.ge.3) then
    write(stdout,*) 'anneal_psc_objective: starting'
endif

nflt = fault%nrows
if (nflt.ne.n) then
    call usage('anneal_psc_objective: input n not equal to nflt (usage:none)')
endif


! Calculate slip in each fault
if (verbosity.ge.3) then
    write(stdout,*) 'anneal_psc_objective: starting psc_slip'
endif
call psc_slip(model,n,slip)
if (verbosity.ge.3) then
    write(stdout,*) 'anneal_psc_objective: finished psc_slip'
endif


! Compute the fit to observations
ndsp = displacement%nrows
ndsp_dof = len_trim(disp_components)*displacement%nrows
nlos = los%nrows
nobs = ndsp_dof + nlos


! Observed displacements
if (.not.allocated(obs)) then
    allocate(obs(nobs),stat=ierr)
    if (ierr.ne.0) then
        call usage('anneal_psc_objective: error allocating memory to obs (usage:none)')
    endif
endif
obs = 0.0d0

! Load observed three-component displacements
if (displacement%file.ne.'none') then
    do i = 1,len_trim(disp_components)
        read(disp_components(i:i),*) icmp
        obs((i-1)*ndsp+1:i*ndsp) = displacement%array(1:ndsp,3+icmp)
    enddo
endif

! Load observed line-of-sight displacements
if (los%file.ne.'none') then
    obs(ndsp_dof+1:ndsp_dof+nlos) = los%array(1:nlos,4)
endif

! Predicted displacements
if (.not.allocated(pre)) then
    allocate(pre(nobs),stat=ierr)
    if (ierr.ne.0) then
        call usage('anneal_psc_objective: error allocating memory to pre (usage:none)')
    endif
endif
pre = 0.0d0

! Load predicted three-component displacements
do i = 1,len_trim(disp_components)
    read(disp_components(i:i),*) icmp
    do idsp = 1,ndsp
        do iflt = 1,nflt
            pre((i-1)*ndsp+idsp) = pre((i-1)*ndsp+idsp) + &
                                   gf_disp%array((icmp-1)*ndsp+idsp,     iflt)*slip(iflt,1)
            pre((i-1)*ndsp+idsp) = pre((i-1)*ndsp+idsp) + &
                                   gf_disp%array((icmp-1)*ndsp+idsp,nflt+iflt)*slip(iflt,2)
        enddo
    enddo
enddo

! Load predicted line-of-sight displacements
do i = 1,nlos
    do iflt = 1,nflt
        pre(ndsp_dof+i) = pre(ndsp_dof+i) + gf_los%array(i,     iflt)*slip(iflt,1)
        pre(ndsp_dof+i) = pre(ndsp_dof+i) + gf_los%array(i,nflt+iflt)*slip(iflt,2)
    enddo
enddo

! do i = 1,nobs
!     write(0,*) i,obs(i),pre(i)
! enddo

! Calculate chi-squared
if (verbosity.ge.3) then
    write(stdout,*) 'anneal_psc_objective: starting misfit_chi2'
endif
call misfit_chi2(obs,pre,cov_matrix,isCovMatrixDiagonal,nobs,anneal_psc_objective)
anneal_psc_objective = -0.5d0*anneal_psc_objective
if (verbosity.ge.3) then
    write(stdout,*) 'anneal_psc_objective: finished misfit_chi2'
endif

if (verbosity.ge.3) then
    write(stdout,*) 'anneal_psc_objective: objective=',anneal_psc_objective
    write(stdout,*) 'anneal_psc_objective: finished'
endif

return
end function

!--------------------------------------------------------------------------------------------------!

subroutine anneal_psc_log(it,temp,obj,model_current,model_proposed,n,isModelAccepted,string)

use io, only: stdout, verbosity
use fltinv, only: fault, &
                  max_iteration, &
                  anneal_log_file

implicit none

! Arguments
integer :: it, n, model_current(n), model_proposed(n)
double precision :: temp, obj
logical :: isModelAccepted
character(len=*) :: string

! Local variables
integer :: i, nflt
double precision :: slip(n,2)
character(len=8) :: rejected_string
character(len=16) :: str, it_str, temp_str, obj_str, uncert_str


if (verbosity.ge.3) then
    write(stdout,*) 'anneal_psc_log: starting on iteration',it
endif

if (anneal_log_file.eq.'') then
    return
endif

nflt = fault%nrows

! Update annealing log file
write(it_str,'(I8)') it
write(temp_str,'(1PE12.4)') temp
write(obj_str,'(1PE12.4)') obj
write(uncert_str,'(1PE12.4)') -1.0d0

it_str = adjustl(it_str)
temp_str = adjustl(temp_str)
obj_str = adjustl(obj_str)
uncert_str = adjustl(uncert_str)


! Update annealing-with-pseudo-coupling log file
if (string.eq.'init') then
    if (verbosity.ge.3) then
        write(stdout,*) 'anneal_psc_log: initializing log file'
    endif

    ! Open the log file
    open(unit=28,file=anneal_log_file,status='unknown')

    ! Compute slip for sub-faults
    call psc_slip(model_current,n,slip)

    ! Write header
    write(str,'(I8)') max_iteration
    write(28,'(A)') '# niterations='//adjustl(str)
    write(str,'(I8)') nflt
    write(28,'(A)') '# nfaults='//adjustl(str)
    write(28,'(A)') '# format=anneal_psc'

    ! Write iteration parameters, fault slip values
    write(28,2801) it_str, temp_str, obj_str, uncert_str
    2801 format('> Iteration=',A8,X,'Temperature=',A12,X,'Objective=',A12,X,'Model_Uncertainty=',&
                A12)

    ! Write locked/unlocked, fault slip results to log file
    do i = 1,n
        write(28,2802) model_current(i),slip(i,1),slip(i,2)
    enddo
    2802 format(I4,X,1PE14.6,E14.6)

elseif (string.eq.'append') then
    ! Is this a rejected model?
    if (isModelAccepted) then
        rejected_string = 'accepted'
    else
        rejected_string = 'rejected'
    endif

    ! Compute slip for sub-faults
    call psc_slip(model_current,n,slip)

    ! Write locked/unlocked, fault slip, old model results to log file
    write(28,2803) it_str, temp_str, obj_str, uncert_str, trim(rejected_string)
    2803 format('> Iteration=',A8,X,'Temperature=',A12,X,'Objective=',A12,X,'Model_Uncertainty=',&
                A12,X,A)

    do i = 1,n
        write(28,2804) model_current(i),slip(i,1),slip(i,2),model_proposed(i)
    enddo
    2804 format(I4,X,1PE14.6,E14.6,I4)

elseif (string.eq.'close') then
    if (verbosity.ge.3) then
        write(stdout,*) 'anneal_psc_log: closing log file'
    endif

    ! Close the log file
    close(28)

else
    call usage('anneal_psc_log: no string option named '//trim(string)//' (usage:none)')
endif

if (verbosity.ge.3) then
    write(stdout,*) 'anneal_psc_log: finished'
endif


return
end subroutine

!--------------------------------------------------------------------------------------------------!


subroutine psc_slip(locked,n,slip)
!----
! Solve for fault slip in unlocked patches surrounding locked patches
!----

#ifdef USE_LAPACK

use io, only: stderr
use solver, only: load_array, load_constraints, solve_dgesv

use fltinv, only: fault, &
                  slip_constraint, &
                  rake_constraint, &
                  prestress, &
                  gf_stress, &
                  Asave

implicit none

! Arguments
integer :: n, locked(n)
double precision :: slip(n,2)

! Local variables
integer :: i, j, ierr, nflt, nrows, ncols
logical :: doInversion, isSlipFixed(2*n)
double precision :: A(2*n,2*n), b(2*n), x(2*n)


! Set dimension variables for arrays
nflt = fault%nrows
if (n.ne.nflt) then
    call usage('psc_slip: size of locked fault array is not equal to nflt (usage:none)')
endif

if (rake_constraint%ncols.eq.1) then
    nrows = nflt
    ncols = nflt
elseif (rake_constraint%ncols.eq.2) then
    nrows = 2*nflt
    ncols = 2*nflt
else
    call usage('psc_slip: number of rake constraints must be 1 or 2 (usage:fault)')
endif


! Set up model matrix Asave with its maximum dimensions and completely full
if (.not.allocated(Asave)) then

    ! Allocate memory for Asave model matrix
    allocate(Asave(nrows,ncols),stat=ierr)
    if (ierr.ne.0) then
        call usage('psc_slip: error allocating memory to Asave (usage:none)')
    endif

    ! Load shear traction Green's functions into model matrix
    call load_array(Asave,nrows,ncols,gf_stress%array(1:nrows,1:ncols),nrows,ncols,1,1,&
                    'gf_stress%array',ierr)
    if (ierr.ne.0) then
        call usage('psc_slip: error loading shear stress GFs into model matrix (usage:none)')
    endif
endif


! There is no need to do the inversion if all the faults are unlocked, so check first
doInversion = .false.
do i = 1,nflt
    if (locked(i).eq.1) then
        doInversion = .true.
        exit
    endif
enddo

if (.not.doInversion) then
    slip = 0.0d0
    return
endif

! Also check if all of the faults are locked...
doInversion = .false.
do i = 1,nflt
    if (locked(i).eq.0) then
        doInversion = .true.
        exit
    endif
enddo

if (.not.doInversion) then
    slip = slip_constraint%array
    return
endif


! Run the inversion for fault slip

! Copy saved model matrix to active model matrix
A(1:nrows,1:ncols) = Asave

! Load pre-stresses into b vector if they have been provided
if (prestress%file.ne.'none') then
    b(1:nflt) = prestress%array(:,1)
    if (rake_constraint%ncols.eq.2) then
        b(nflt+1:2*nflt) = prestress%array(:,2)
    endif
else
    b(1:nrows) = 0.0d0
endif

! Load slip constraints for locked faults
isSlipFixed = .false.

do i = 1,nflt
    if (locked(i).eq.1.and.abs(slip_constraint%array(i,1)).lt.99998.0d0) then
        isSlipFixed(i) = .true.
    endif
    if (rake_constraint%ncols.eq.2) then
        if (locked(i).eq.1.and.abs(slip_constraint%array(i,2)).lt.99998.0d0) then
            isSlipFixed(i+nflt) = .true.
        endif
    endif
enddo

if (rake_constraint%ncols.eq.1) then
    call load_constraints(A,b,nrows,ncols,slip_constraint%array(:,1),isSlipFixed,ierr)
    if (ierr.ne.0) then
        call usage('psc_slip: error loading slip constraints into A and b arrays (usage:none)')
    endif
elseif (rake_constraint%ncols.eq.2) then
    call load_constraints(A,b,nrows,ncols, &
                          [slip_constraint%array(:,1),slip_constraint%array(:,2)],isSlipFixed, &
                          ierr)
    if (ierr.ne.0) then
        call usage('psc_slip: error loading slip constraints into A and b arrays (usage:none)')
    endif
endif

! Remove stress rows for fixed slip faults; these no longer contribute to the solution
j = 0
do i = 1,nrows
    if (isSlipFixed(i)) then
        ! Do not add row i to A
    else
        ! Add row i to A
        j = j + 1
        A(j,:) = A(i,:)
        b(j) = b(i)
    endif
enddo
nrows = j

! Solve for fault slip in unlocked faults
if (nrows.ne.ncols) then
    write(stderr,*) 'psc_slip: nrows not equal to ncols'
    write(stderr,*) '    nrows=',nrows
    write(stderr,*) '    ncols=',ncols
    call usage('Inversion routine gesv requires a square A matrix (usage:none)')
endif
call solve_dgesv(A(1:nrows,1:ncols),b(1:nrows),x,nrows,ierr)
if (ierr.ne.0) then
    call usage('psc_slip: error inverting for slip with gesv (usage:none)')
endif

! Load fault slip into the output array
j = 0
do i = 1,nflt
    if (isSlipFixed(i)) then
        ! Fixed slip is read from the slip constraint array
        slip(i,1) = slip_constraint%array(i,1)
        if (rake_constraint%ncols.eq.2) then
            slip(i,2) = slip_constraint%array(i,2)
        endif
    else
        ! Unlocked slip is read from the inverted slip array
        j = j + 1
        slip(i,1) = x(j)
        if (rake_constraint%ncols.eq.2) then
            slip(i,2) = x(j+nrows/2)
        endif
    endif
    ! write(0,*) 'psc_slip',i, slip(i,:)
enddo

#else

call usage('psc_slip: unable to compute pseudo-coupling slip when compiled without LAPACK '//&
           '(usage:none)')

#endif

return
end subroutine psc_slip
