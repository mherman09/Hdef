subroutine invert_anneal()
!----
! Determine fault slip with a simulated annealing search, varying slip magnitude and rake angle
! randomly to find the best fitting solution.
!----

use trig, only: d2r
use annealing, only: anneal

use fltinv, only: fault, &
                  max_iteration, &
                  reset_iteration, &
                  temp_start, &
                  temp_minimum, &
                  cooling_factor, &
                  fault_slip, &
                  modelUncertainty

implicit none

! Interface to driver subroutines
interface
    subroutine anneal_init(model,n)
        integer :: n
        double precision :: model(n)
    end subroutine
    subroutine anneal_propose(model_current,model_proposed,n)
        integer :: n
        double precision :: model_current(n)
        double precision :: model_proposed(n)
    end subroutine
    function anneal_objective(model,n)
        integer :: n
        double precision :: model(n)
        double precision :: anneal_objective
    end function
    subroutine anneal_log(it,temp,obj,model_current,model_proposed,n,isAccepted,string)
        integer :: it, n
        double precision :: temp, obj
        double precision :: model_current(n)
        double precision :: model_proposed(n)
        logical :: isAccepted
        character(len=*) :: string
    end subroutine
end interface

! Local variables
integer :: i, ierr, nflt_dof
double precision, allocatable :: model_best(:)
logical :: saveRejected


saveRejected = .true.

! Number of parameters to include in search depends on whether searching for model uncertainty
if (modelUncertainty) then
    nflt_dof = 2*fault%nrows + 1
else
    nflt_dof = 2*fault%nrows
endif
allocate(model_best(nflt_dof),stat=ierr)
if (ierr.ne.0) then
    call usage('invert_anneal: error allocating memory to model_best (usage:none)')
endif

! Call anneal routine with specific driver routines anneal_init, anneal_propose, and anneal_objective
call anneal(nflt_dof, &
            model_best, &
            anneal_init, &
            anneal_propose, &
            anneal_objective, &
            max_iteration, &
            reset_iteration, &
            temp_start, &
            temp_minimum, &
            cooling_factor, &
            anneal_log)

! Save results for printing
do i = 1,fault%nrows
    fault_slip(i,1) = model_best(i)*cos(model_best(i+fault%nrows)*d2r)
    fault_slip(i,2) = model_best(i)*sin(model_best(i+fault%nrows)*d2r)
enddo

return
end subroutine invert_anneal

!--------------------------------------------------------------------------------------------------!

subroutine anneal_init(model,n)
!----
! Initialize the annealing variables:
!     - iseed: random number seed
!     - model array: slip magnitude and rake angle
!     - step array: standard deviations of the Gaussian step size for slip and rake
!----

use io, only: stdout, fileExists, line_count, verbosity
use random, only: iseed, timeseed, r8_uniform_01

use fltinv, only: fault, &
                  slip_constraint, &
                  rake_constraint, &
                  anneal_init_mode, &
                  anneal_init_file, &
                  anneal_step_file, &
                  step, &
                  anneal_seed, &
                  modelUncertainty

implicit none

! Arguments
integer :: n
double precision :: model(n)

! Local variables
integer :: i, ios, nflt


if (verbosity.ge.2) then
    write(stdout,*) 'anneal_init: starting'
endif

! Initialize the random number generator seed
if (anneal_seed.eq.0) then
    iseed = -timeseed()
else
    iseed = -abs(anneal_seed)
endif

! Check array dimensions:
!     The first nflt entries in the model array are fault slip magnitude
!     The next nflt entries are rake angles
!     If searching for model uncertainty, then add to the end of the array
nflt = fault%nrows
if (modelUncertainty) then
    if (n.ne.2*nflt+1) then
        call usage('anneal_init: looking for model uncertainty, input n not equal to 2*nflt+1 '//&
                   '(usage:none)')
    endif
else
    if (n.ne.2*nflt) then
        call usage('anneal_init: not looking for model uncertainty, input n not equal to 2*nflt '//&
                   '(usage:none)')
    endif
endif


! Initialize the model array values
model = 0.0d0
if (anneal_init_mode.eq.'mean'.or.anneal_init_mode.eq.'half') then
    ! Set the initial solution to the magnitude and rake halfway between their limits
    do i = 1,nflt
        model(i)      = 0.5d0*(slip_constraint%array(i,1)+slip_constraint%array(i,2))
        model(i+nflt) = 0.5d0*(rake_constraint%array(i,1)+rake_constraint%array(i,2))
    enddo

elseif (anneal_init_mode.eq.'min') then
    model(1:nflt) = slip_constraint%array(:,1)
    model(nflt+1:2*nflt) = rake_constraint%array(:,1)

elseif (anneal_init_mode.eq.'max') then
    ! Set the initial solution to zero
    model(1:nflt) = slip_constraint%array(:,2)
    model(nflt+1:2*nflt) = rake_constraint%array(:,2)

elseif (anneal_init_mode.eq.'rand') then
    ! Set the initial solution to a random value within the limits of the constraint arrays
    do i = 1,nflt
        model(i)      = slip_constraint%array(i,1) + &
                        r8_uniform_01(iseed)*(slip_constraint%array(i,2)-slip_constraint%array(i,1))
        model(i+nflt) = rake_constraint%array(i,1) + &
                        r8_uniform_01(iseed)*(rake_constraint%array(i,2)-rake_constraint%array(i,1))
    enddo

elseif (anneal_init_mode.eq.'user') then
    ! Read the initial solution from a file (slip rake)
    if (anneal_init_file.eq.'none') then
        call usage('anneal_init: no initialization file defined for anneal_init_mode=user '//&
                   '(usage:anneal)')
    elseif (.not.fileExists(anneal_init_file)) then
        call usage('anneal_init: no anneal_init_file found named "'//trim(anneal_init_file)//'" '//&
                   '(usage:anneal)')
    endif
    if (line_count(anneal_init_file).ne.nflt) then
        call usage('anneal_init: number of lines in anneal_init_file must be equal to nflt '//&
                   '(usage:anneal)')
    endif
    open(unit=29,file=anneal_init_file,status='old')
    do i = 1,nflt
        read(29,*,iostat=ios) model(i),model(i+nflt)
        if (ios.ne.0) then
            call usage('anneal_init: error reading anneal init file (usage:anneal)')
        endif
    enddo
    close(29)

else
    call usage('anneal_init: no initialization mode named "'//trim(anneal_init_mode)//'" '//&
               '(usage:anneal)')
endif


! Set the parameter step values
allocate(step(2*nflt))
if (anneal_step_file.eq.'none') then
    step(1:nflt)        = (slip_constraint%array(:,2)-slip_constraint%array(:,1))/50.0d0
    step(nflt+1:2*nflt) = (rake_constraint%array(:,2)-rake_constraint%array(:,1))/50.0d0
else
    if (.not.fileExists(anneal_step_file)) then
        call usage('anneal_init: no anneal_step_file found named "'//trim(anneal_step_file)//'" '//&
                   '(usage:anneal)')
    endif
    if (line_count(anneal_step_file).ne.nflt) then
        call usage('anneal_init: anneal step file must have nflt lines (usage:anneal)')
    endif
    open(unit=30,file=anneal_step_file,status='old')
    do i = 1,nflt
        read(30,*,iostat=ios) step(i),step(i+nflt)
        if (ios.ne.0) then
            call usage('anneal_init: error reading anneal step file (usage:anneal)')
        endif
    enddo
    close(30)
endif


! If searching for model uncertainty, set the initial value here (between 0 and 1)
if (modelUncertainty) then
    model(2*nflt+1) = r8_uniform_01(iseed)
endif


if (verbosity.ge.2) then
    write(stdout,*) 'anneal_init: finished'
endif
if (verbosity.ge.3) then
    write(stdout,*) 'Initial model:'
    do i = 1,nflt
        write(stdout,*) model(i),model(i+nflt)
    enddo
endif


return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine anneal_propose(model_in,model_out,n)
!----
! Propose a new model from a multi-dimensional normal distribution around the current model
!----

use io, only: verbosity, stdout
use random, only: iseed, r8_normal_ab

use fltinv, only: fault, &
                  slip_constraint, &
                  rake_constraint, &
                  step, &
                  modelUncertainty

implicit none

! Arguments
integer :: n
double precision :: model_in(n), model_out(n)

! Local variables
integer :: i, j, nflt


if (verbosity.ge.3) then
    write(stdout,*) 'anneal_propose: starting'
endif


nflt = fault%nrows
if (modelUncertainty) then
    if (n.ne.2*nflt+1) then
        call usage('anneal_propose: looking for model uncertainty, input n not equal to 2*nflt+1 '//&
                   '(usage:none)')
    endif
else
    if (n.ne.2*nflt) then
        call usage('anneal_propose: not looking for model uncertainty, input n not equal to '//&
                   '2*nflt (usage:none)')
    endif
endif


do i = 1,nflt
    j = i+nflt

    ! The new model is sampled from a Gaussian distribution around the current model
    model_out(i) = r8_normal_ab(model_in(i),step(i),iseed)
    model_out(j) = r8_normal_ab(model_in(j),step(j),iseed)

    ! Make sure model values do not go beyond values in constraint files

    ! Hard boundary for slip magnitude is value set in slip_constraint%array
    if (model_out(i).lt.slip_constraint%array(i,1)) then
        model_out(i) = slip_constraint%array(i,1)
    elseif (model_out(i).gt.slip_constraint%array(i,2)) then
        model_out(i) = slip_constraint%array(i,2)
    endif

    ! If range is 360, rake angle can vary continuously
    ! Otherwise, the rake constraints are hard boundaries
    if (model_out(j).lt.rake_constraint%array(i,1)) then
        if (rake_constraint%array(i,2)-rake_constraint%array(i,1).ge.360.0d0) then
            model_out(j) = model_out(j)+360.0d0
        else
            model_out(j) = rake_constraint%array(i,1)
        endif
    elseif (model_out(j).gt.rake_constraint%array(i,2)) then
        if (rake_constraint%array(i,2)-rake_constraint%array(i,1).ge.360.0d0) then
            model_out(j) = model_out(j)-360.0d0
        else
            model_out(j) = rake_constraint%array(i,2)
        endif
    endif
enddo


! If searching for model uncertainty, propose a new value here
if (modelUncertainty) then
    model_out(2*nflt+1) = r8_normal_ab(model_in(2*nflt+1),0.05d0,iseed)
    if (model_out(2*nflt+1).le.1.0d-8) then
        model_out(2*nflt+1) = 1.0d-8
    endif
endif


if (verbosity.ge.3) then
    write(stdout,*) 'anneal_propose: finished'
    write(stdout,*) 'Proposed model:'
    do i = 1,nflt
        write(stdout,*) model_out(i),model_out(i+nflt)
    enddo
endif

return
end subroutine

!--------------------------------------------------------------------------------------------------!

function anneal_objective(model,n)
!----
! The objective function is -0.5 times the chi-squared misfit.
!
! If model uncertainty (variance) is included, then the objective function needs to also
! include the pre-exponential coefficient, which includes a covariance matrix that changes with
! the model uncertainty value.
!----

use trig, only: d2r
use misfit, only: misfit_chi2

use fltinv, only: fault, &
                  displacement, &
                  disp_components, &
                  los, &
                  cov_matrix, &
                  isCovMatrixDiagonal, &
                  gf_disp, &
                  gf_los, &
                  modelUncertainty

implicit none

! Arguments
integer, intent(in) :: n
double precision :: model(n), anneal_objective

! Local variables
integer :: i, ierr, idsp, iflt, icmp, nflt, ndsp, ndsp_dof, nlos, nobs
double precision, allocatable :: obs(:), pre(:)
double precision :: slip, rake
double precision, allocatable :: model_cov_matrix(:)
double precision :: log_cov_matrix_determinant


nflt = fault%nrows
if (modelUncertainty) then
    if (n.ne.2*nflt+1) then
        call usage('anneal_objective: looking for model uncertainty, input n not equal to '//&
                   '2*nflt+1 (usage:none)')
    endif
else
    if (n.ne.2*nflt) then
        call usage('anneal_objective: not looking for model uncertainty, input n not equal to '//&
                   '2*nflt (usage:none)')
    endif
endif


ndsp = displacement%nrows
ndsp_dof = len_trim(disp_components)*displacement%nrows
nlos = los%nrows
nobs = ndsp_dof + nlos


! Observed displacements
if (.not.allocated(obs)) then
    allocate(obs(nobs),stat=ierr)
    if (ierr.ne.0) then
        call usage('anneal_objective: error allocating memory to obs (usage:none)')
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
        call usage('anneal_objective: error allocating memory to pre (usage:none)')
    endif
endif
pre = 0.0d0

! Load predicted three-component displacements
do i = 1,len_trim(disp_components)
    read(disp_components(i:i),*) icmp
    do idsp = 1,ndsp
        do iflt = 1,nflt
            slip = model(iflt)
            rake = model(nflt+iflt)
            pre((i-1)*ndsp+idsp) = pre((i-1)*ndsp+idsp) + &
                                   gf_disp%array((icmp-1)*ndsp+idsp,     iflt)*slip*cos(rake*d2r)
            pre((i-1)*ndsp+idsp) = pre((i-1)*ndsp+idsp) + &
                                   gf_disp%array((icmp-1)*ndsp+idsp,nflt+iflt)*slip*sin(rake*d2r)
        enddo
    enddo
enddo

! Load predicted line-of-sight displacements
do i = 1,nlos
    do iflt = 1,nflt
        slip = model(iflt)
        rake = model(nflt+iflt)
        pre(ndsp_dof+i) = pre(ndsp_dof+i) + gf_los%array(i,     iflt)*slip*cos(rake*d2r)
        pre(ndsp_dof+i) = pre(ndsp_dof+i) + gf_los%array(i,nflt+iflt)*slip*sin(rake*d2r)
    enddo
enddo


if (modelUncertainty) then
    if (.not.isCovMatrixDiagonal) then
        call usage('anneal_objective: covariance matrix must be diagonal to estimate model '//&
                   'uncertainty for now (usage:input)')
    endif
    if (los%file.ne.'none') then
        call usage('anneal_objective: model uncertainty cannot be estimated with LOS data yet '//&
                   '(usage:none)')
    endif

    ! Initialize model prediction error covariance matrix
    if (.not.allocated(model_cov_matrix)) then
        allocate(model_cov_matrix(nobs),stat=ierr)
        if (ierr.ne.0) then
            call usage('anneal_objective: error allocating memory to model_cov_matrix (usage:none)')
        endif
    endif
    model_cov_matrix = 0.0d0

    ! Load model error prediction covariance matrix
    do idsp = 1,ndsp
        ! Model prediction error covariance at each station is (for three-component displacements):
        !     error_cov(ista) = e^2*(obs(ista,1)^2+obs(ista,2)^2+obs(ista,3)^2)
        ! where e is the percent error
        do icmp = 1,len_trim(disp_components)
            i = (icmp-1)*ndsp + idsp
            model_cov_matrix(idsp) = model_cov_matrix(idsp) + model(2*nflt+1)**2 * obs(i)**2
        enddo
    enddo
    do icmp = 2,len_trim(disp_components)
        model_cov_matrix((icmp-1)*ndsp+1:icmp*ndsp) = model_cov_matrix(1:ndsp)
    enddo

    ! Calculate chi-squared and scale it for (Gaussian) probability
    cov_matrix(:,1) = cov_matrix(:,1) + model_cov_matrix
    call misfit_chi2(obs,pre,cov_matrix,isCovMatrixDiagonal,nobs,anneal_objective)
    cov_matrix(:,1) = cov_matrix(:,1) - model_cov_matrix

    ! Scale objective value by determinant of the covariance matrix (including both observation
    ! and model prediction error covariances) to determine its probability, e.g.:
    !
    ! Tarantola, A. (2005). Inverse Problem Theory and Methods for Model Parameter Estimation.
    !     Philadelphia, PA, USA: Society for Industrial and Applied Mathematics.
    ! Beck, J.L. (2010). Bayesian system identification based on probability logic. Structural
    !     Control and Health Monitoring 17, 825–847.
    ! Minson, S.E., Simons, M., & Beck, J.L. (2013). Bayesian inversion for finite fault earthquake
    !     source models I—theory and algorithm. Geophysical Journal International 194, 1701–1726.
    log_cov_matrix_determinant = 0.0d0
    do i = 1,nobs
        log_cov_matrix_determinant = log_cov_matrix_determinant + &
                                     log(cov_matrix(i,1) + model_cov_matrix(i))
    enddo
    anneal_objective = anneal_objective + log_cov_matrix_determinant

else
    call misfit_chi2(obs,pre,cov_matrix,isCovMatrixDiagonal,nobs,anneal_objective)
endif
anneal_objective = -0.5d0*anneal_objective

return
end function


!--------------------------------------------------------------------------------------------------!

subroutine anneal_log(it,temp,obj,model_current,model_proposed,n,isAccepted,string)

use fltinv, only: fault, &
                  max_iteration, &
                  anneal_log_file, &
                  modelUncertainty

implicit none

! Arguments
integer :: it, n
double precision :: temp, obj, model_current(n), model_proposed(n)
character(len=*) :: string
logical :: isAccepted

! Local variables
integer :: i, nflt
character(len=8) :: rejected_string
character(len=16) :: str, it_str, temp_str, obj_str, uncert_str


if (anneal_log_file.eq.'') then
    return
endif


nflt = fault%nrows

! Update annealing log file
write(it_str,'(I8)') it
write(temp_str,'(1PE12.4)') temp
write(obj_str,'(1PE12.4)') obj
if (modelUncertainty) then
    write(uncert_str,'(1PE12.4)') model_current(2*nflt+1)
else
    write(uncert_str,'(1PE12.4)') -1.0d0
endif
it_str = adjustl(it_str)
temp_str = adjustl(temp_str)
obj_str = adjustl(obj_str)
uncert_str = adjustl(uncert_str)


if (string.eq.'init') then
    ! Open the log file
    open(unit=29,file=anneal_log_file,status='unknown')

    ! Write header

    write(str,'(I8)') max_iteration
    write(29,'(A)') '# niterations='//adjustl(str)
    write(str,'(I8)') nflt
    write(29,'(A)') '# nfaults='//adjustl(str)
    write(29,'(A)') '# format=anneal'

    ! Write iteration parameters, fault slip values
    write(29,2901) it_str, temp_str, obj_str, uncert_str
    do i = 1,nflt
        write(29,2902) model_current(i)
    enddo
    do i = nflt+1,2*nflt
        write(29,2903) model_current(i)
    enddo
    2901 format('> Iteration=',A8,X,'Temperature=',A12,X,'Objective=',A12,X,'Model_Uncertainty=',&
                A12)
    2902 format(1PE14.6)
    2903 format(F12.2)

elseif (string.eq.'append') then
    ! Is this a rejected model?
    if (isAccepted) then
        rejected_string = 'accepted'
    else
        rejected_string = 'rejected'
    endif

    ! Write iteration parameters, fault slip values
    write(29,2904) it_str, temp_str, obj_str, uncert_str, trim(rejected_string)
    do i = 1,nflt
        write(29,2905) model_current(i),model_proposed(i)
    enddo
    do i = nflt+1,2*nflt
        write(29,2906) model_current(i),model_proposed(i)
    enddo
    2904 format('> Iteration=',A8,X,'Temperature=',A12,X,'Objective=',A12,X,'Model_Uncertainty=',&
                A12,X,A)
    2905 format(1PE14.6,X,E14.6)
    2906 format(F12.2,X,F12.2)

elseif (string.eq.'close') then
    ! Close the log file
    close(29)

else
    call usage('anneal_log: no string option named '//trim(string)//' (usage:none)')
endif

return
end subroutine
