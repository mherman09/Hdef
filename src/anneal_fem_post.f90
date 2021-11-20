!--------------------------------------------------------------------------------------------------!
! ANNEAL_FEM_POST
!
! Post-process a simulated annealing search output file from GTECTON. In general, the format is:
!
!     > ITERATION=I OBJECTIVE=X ACCEPTED=YES|NO
!     PARAMETER1=VALUE
!     PARAMETER2=VALUE
!     :
!
!--------------------------------------------------------------------------------------------------!

module anneal_fem_post

character(len=512) :: anneal_log_file                       ! Name of annealing search log file
integer :: lu_log_file                                      ! Logical unit of log file for interaction

integer :: nparam                                           ! Number of model parameters
integer :: noutputs                                         ! Number of statistical outputs
character(len=32), allocatable :: param_names(:)            ! Parameter name array
double precision, allocatable :: param_values(:,:)          ! Parameter value array
character(len=32), allocatable :: param_outputs(:)          ! Parameter outputs requested
character(len=512), allocatable :: param_output_files(:)    ! Parameter output file names
logical :: listParameters                                   ! Switch to print available parameters names

integer :: niterations                                      ! Number of iterations in search
integer :: burnin_iterations

character(len=512) :: best_model_file
character(len=512) :: mean_model_file
character(len=512) :: median_model_file

double precision, allocatable :: obj(:)

! integer :: nmarg1d
! integer :: marg1d_flt(1000)
! character(len=64) :: slip_marg_output_file(1000)
!
! character(len=512) :: mean_pole_output_file
! integer :: nmarg_pole
! integer :: marg_pole(20)
! character(len=64) :: pole_marg_output_file(20)
!
! integer :: nflocked
! character(len=64) :: flocked_subset_file(10)
! character(len=64) :: flocked_output_file(10)
!
! character(len=512) :: obj_file
! character(len=512) :: uncert_file
!
! logical :: doResample
! integer :: nresamp
! integer, allocatable :: samp(:)

end module


!==================================================================================================!


program main

use io, only: stderr, stdout, verbosity
! use earth, only: pole_geo2xyz, pole_xyz2geo
! use annealing, only: resample
! use random, only: iseed, timeseed
! use sort, only: heapsort

use anneal_fem_post, only: anneal_log_file, &
                           lu_log_file, &
                           nparam, &
                           noutputs, &
                           param_names, &
                           param_values, &
                           param_outputs, &
                           param_output_files, &
                           niterations, &
                           burnin_iterations, &
                           best_model_file, &
                           mean_model_file, &
                           obj
!                        plocked_output_file, &
!                        mean_slip_output_file, &
!                        nmarg1d, &
!                        marg1d_flt, &
!                        slip_marg_output_file, &
!                        mean_pole_output_file, &
!                        nmarg_pole, &
!                        marg_pole, &
!                        pole_marg_output_file, &
!                        nflocked, &
!                        flocked_subset_file, &
!                        flocked_output_file, &
!                        obj_file, &
!                        uncert_file, &
!                        doResample, &
!                        nresamp, &
!                        samp

implicit none

! Local variables
integer :: i, j
integer :: ios
integer :: ioutput
integer :: iparam
integer :: it
logical :: fileExists
logical :: foundParameterName
character(len=32) :: param_output_name
double precision :: min_obj
double precision, allocatable :: mean(:)



! Initialize status indicator
ios = 0


! Parse command line
call gcmdln()
if (verbosity.ge.1) then
    write(stdout,*) 'anneal_fem_post: starting'
endif


!----
! Open annealing log file
!----
! Check that file exists
inquire(file=anneal_log_file,exist=fileExists)
if (.not.fileExists) then
    call usage('anneal_fem_post: no anneal log file found named '//trim(anneal_log_file))
endif

! Open the log file
lu_log_file = 21
open(unit=lu_log_file,file=anneal_log_file,status='old')
if (verbosity.ge.2) then
    write(stdout,*) 'anneal_fem_post: reading log file "',trim(anneal_log_file),'"'
endif


!----
! Read annealing log file
!----
! Count the number of parameters and iterations, and set parameter names in name array
call read_model_parameters()

! Set up arrays
allocate(obj(niterations))                ! DP: Objective function values
! allocate(uncert(nit))             ! DP: Model prediction error
! allocate(isAccepted(nit))         ! Logical: Was the proposed model accepted?
! if (verbosity.ge.2) then
!     write(stdout,*) 'anneal_fem_post: finished allocating array memory'
! endif

! Read outputs
call read_models()



! ! Resample ensemble if requested
! if (doResample) then
!     if (nresamp.eq.0) then
!         nresamp = nit
!     endif
!     if (nresamp.gt.nit) then
!         call usage('anneal_fem_post: nresamp cannot be greater than nit')
!     endif
!
!     ! Read all objectives (log-posterior-probabilities)
!     if (verbosity.ge.2) then
!         write(stdout,*) 'anneal_fem_post: reading model objective functions for resampling'
!     endif
!     i = 0
!     do
!         read(lu_log_file,'(A)',iostat=ios) input_line
!         if (ios.ne.0) then
!             exit
!         endif
!         if (index(input_line,'>').eq.1) then
!             if (i.ge.1) then
!                 call parse_iteration_header(input_line,it,ddum,obj(i),ddum,ldum,ios)
!             endif
!             i = i + 1
!         endif
!     enddo
!     if (verbosity.ge.2) then
!         write(stdout,*) 'anneal_fem_post: finished reading model objective functions'
!     endif
!
!     ! Set nit to number read in previous block
!     if (i-1.ne.nit) then
!         nit = i-1
!         ! Arrays with multiple dimensions need to be reallocated, so values remain in
!         ! correct iteration index
!         deallocate(locked)
!         deallocate(slip)
!         deallocate(poles)
!         allocate(locked(nit,nflt))
!         allocate(slip(nit,nflt,2))
!         allocate(poles(nit,npoles,3))
!     endif
!
!     ! Initialize the random number seed
!     if (iseed.eq.0) then
!         iseed = timeseed()
!     endif
!
!     ! Resample the ensemble
!     if (verbosity.ge.2) then
!         write(stdout,*) 'anneal_fem_post: resampling models'
!     endif
!     allocate(samp(nit))
!     allocate(obj_resamp(nit))
!     max_obj = maxval(obj(1:nit))
!     do i = 1,nit
!         obj_resamp(i) = exp(obj(i) - max_obj)
!     enddo
!     call resample(obj_resamp,samp,nit)
!     if (verbosity.ge.2) then
!         write(stdout,*) 'anneal_fem_post: finished resampling'
!     endif
!
!     ! Save the first nresamp models
!     allocate(samp_tmp(nresamp))
!     samp_tmp = samp(1:nresamp)
!     deallocate(samp)
!     allocate(samp(nresamp))
!     samp = samp_tmp
!
!     ! Sort the sampled models for more efficient reading
!     call heapsort(samp,nresamp)
!
!     ! Reset file to read again (do not overwrite nit set above)
!     rewind(lu_log_file)
!     call read_header(idum,nflt,npoles,log_format)
!
! 	! Read search results, saving only sampled models
!     call read_resampled_models(nit,nflt,npoles,log_format,locked,obj,slip,poles,temp,uncert, &
!                                isAccepted)
!     nit = nresamp
!
! else
!
!     ! Read all of the search results
!     call read_models(nit,nflt,npoles,log_format,locked,obj,slip,poles,temp,uncert,isAccepted)
!
! endif



!----
! Print requested quantities
!----
! Marginal distributions for individual parameters
do ioutput = 1,noutputs

    ! Parse input parameter output string from command line
    i = index(param_outputs(ioutput),'-')
    param_outputs(ioutput)(i:i) = ' '
    i = index(param_outputs(ioutput),':')
    param_outputs(ioutput)(i:i) = ' '
    read(param_outputs(ioutput),*) param_output_name
    ! print *,param_outputs(ioutput)

    ! Determine which parameter index is requested
    foundParameterName = .false.
    do iparam = 1,nparam
        ! print *,'    ',iparam,param_names(iparam)
        if (param_output_name.eq.param_names(iparam)) then
            foundParameterName = .true.
            exit
        endif
    enddo
    if (foundParameterName) then
        ! print *,iparam,param_names(iparam)
    else
        write(stderr,*) 'anneal_fem_post: could not find parameter name ',trim(param_output_name)
        write(stderr,*) 'Skipping...'
        cycle
    endif

    ! Print parameter values to file
    open(unit=101,file=param_output_files(ioutput))
    do i = burnin_iterations+1,niterations
        write(101,*) param_values(iparam,i)
    enddo
    close(101)
enddo


! Best fitting model
if (best_model_file.ne.'') then
    ! write(*,*) 'anneal_fem_post: working on best-fitting output file'
    ! Find minimum chi-squared model
    min_obj = 1d10
    it = 1
    do i = 1,niterations
        if (abs(obj(i)).lt.min_obj) then
            it = i
            min_obj = abs(obj(i))
        endif
    enddo

    ! Write best model parameters to file
    write(stdout,*) 'best fit: ',it,obj(it)
    open(unit=22,file=best_model_file,status='unknown')
    do i = 1,nparam
        write(22,*) param_names(i),param_values(i,it)
    enddo
    close(22)
endif


! Mean slip (AND STANDARD DEVIATION???)
if (mean_model_file.ne.'') then
    ! Compute mean parameter values
    allocate(mean(nparam))
    mean = 0.0d0

    do i = 1,nparam
        do j = burnin_iterations+1,niterations
            mean(i) = mean(i) + param_values(i,j)
        enddo
    enddo
    mean = mean/(dble(niterations-burnin_iterations))

    ! Write mean parameter values to file
    open(unit=24,file=mean_model_file,status='unknown')
    do i = 1,nparam
        write(24,*) param_names(i),mean(i)
    enddo
    close(24)

    deallocate(mean)
endif


! ! Annealing objective versus iteration
! if (obj_file.ne.'') then
!     open(unit=29,file=obj_file,status='unknown')
!     do i = 1,nit
!         write(29,*) i,-obj(i),isAccepted(i)
!     enddo
!     close(29)
! endif
!
!
! ! Uncertainty histogram
! if (uncert_file.ne.'') then
!     open(unit=30,file=uncert_file,status='unknown')
!     do i = 1,nit
!         write(30,*) uncert(i)
!     enddo
!     close(30)
! endif
!
!
! !----
! !       CLEAN UP
! !----
! ! Deallocate arrays
! if (allocated(locked)) then
!     deallocate(locked)
! endif
! if (allocated(obj)) then
!     deallocate(obj)
! endif
! if (allocated(slip)) then
!     deallocate(slip)
! endif
!
! ! Close annealing log file
! close(lu_log_file)


if (verbosity.ge.1) then
    write(stdout,*) 'anneal_fem_post: finished'
endif

end

!--------------------------------------------------------------------------------------------------!


subroutine read_model_parameters()

use io, only: stdout, verbosity

use anneal_fem_post, only: lu_log_file, &
                           nparam, &
                           param_names, &
                           param_values, &
                           listParameters, &
                           niterations

implicit none


! Local variables
integer :: i, j
integer :: ios
logical :: paramsCounted
character(len=1) :: iteration_line
character(len=1) :: param_line
character(len=512) :: input_line


if (verbosity.ge.1) then
    write(stdout,*) 'read_model_parameters: starting'
endif


! Initialize status
ios = 0


! Initialize counts
paramsCounted = .false.
nparam = 0
niterations = 0


! Count number of iterations and parameters
do
    ! First line of each iteration is a header starting with ">"
    read(lu_log_file,'(A)',end=4001,iostat=ios) iteration_line
    niterations = niterations + 1
    if (iteration_line.ne.'>') then
        write(0,*) iteration_line
        call usage('anneal_fem_post: output file has formatting problem')
    endif

    ! Count number of parameters
    if (.not.paramsCounted) then
        do
            read(lu_log_file,'(A)',end=4002,iostat=ios) param_line
            if (param_line.eq.'>') then
                backspace(lu_log_file)
                exit
            else
                nparam = nparam + 1
            endif
        enddo
        paramsCounted = .true.
    else
        do i = 1,nparam
            read(lu_log_file,'(A)',end=4002,iostat=ios) param_line
        enddo
    endif

    if (listParameters) then
        exit
    endif
enddo
4001 ios = 0  ! Finished do loop without error
4002 if (ios.ne.0) then
    call usage('anneal_fem_post: reached end of file in middle of parameter list')
endif
write(*,*) 'nparam=     ',nparam
write(*,*) 'niterations=',niterations


! Allocate memory for parameter names and values
allocate(param_names(nparam))
allocate(param_values(nparam,niterations))


! Reset log file for reading parameter names and values
rewind(lu_log_file)


! Read parameter names
read(lu_log_file,'(A)',end=4003,iostat=ios) input_line
do i = 1,nparam
    read(lu_log_file,'(A)',end=4003,iostat=ios) input_line
    j = index(input_line,'=')
    input_line(j:j) = ' '
    read(input_line,*) param_names(i)

enddo
4003 if (ios.ne.0) then
    call usage('anneal_fem_post: reached EOF while setting parameter name array')
endif


! Print parameter names if requested
if (listParameters) then
    do i = 1,nparam
        print *,trim(param_names(i))
    enddo
    stop
endif


! Reset log file for reading
rewind(lu_log_file)


if (verbosity.ge.1) then
    write(stdout,*) 'read_model_parameters: finished'
endif


return
end subroutine


!--------------------------------------------------------------------------------------------------!


subroutine read_models()

use io, only: stdout, verbosity

use anneal_fem_post, only: lu_log_file, &
                           nparam, &
                           param_values, &
                           niterations, &
                           obj

implicit none

! Local variables
integer :: i, j, k
integer :: ios
character(len=1) :: junk
character(len=512) :: input_line
character(len=64) :: header_item



if (verbosity.ge.1) then
    write(stdout,*) 'read_models: starting'
endif


! Read all model parameter values into parameter value array
do i = 1,niterations

    ! Iterations start with a single line; read objective functions
    read(lu_log_file,'(A)',end=1001,iostat=ios) input_line
    do while (input_line.ne.'')
        read(input_line,*) header_item
        if (index(header_item,'OBJECTIVE=').ne.0) then
            j = index(header_item,'=')
            header_item(j:j) = ' '
            read(header_item,*) junk, obj(i)
            header_item(j:j) = '='
        endif
        j = index(input_line,trim(header_item))
        input_line(j:j+len_trim(header_item)) = ''
    enddo

    ! Read through parameters and save them to parameter value array
    do j = 1,nparam
        read(lu_log_file,'(A)',end=1001,iostat=ios) input_line

        k = index(input_line,'=')
        input_line(k:k) = ' '
        read(input_line,*) junk,param_values(j,i)
    enddo
enddo

1001 continue


if (verbosity.ge.1) then
    write(stdout,*) 'read_models: finished'
endif


return
end subroutine


!--------------------------------------------------------------------------------------------------!
!
!
! subroutine read_resampled_models(nit,nflt,npoles,log_format,locked,obj,slip,poles,temp, &
!                                  uncert,isAccepted)
!
! use io, only: stdout, stderr, verbosity, progress_indicator
! use trig, only: d2r
!
! use anneal_fem_post, only: lu_log_file, &
!                        nresamp, &
!                        samp
!
! implicit none
!
! ! Arguments
! integer :: nit, nflt, npoles
! character(len=*) :: log_format
! integer :: locked(nresamp,nflt)
! double precision :: obj(nresamp), slip(nresamp,nflt,2), poles(nresamp,npoles,3), temp(nresamp), &
!                     uncert(nresamp)
! logical :: isAccepted(nresamp)
!
! ! Local variables
! integer :: i, ios, it, itread, iflt, ipole, nlines_per_it, isamp
! character(len=512) :: input_line
! double precision :: rake
! integer :: locked_it(nflt)
! double precision :: obj_it, slip_it(nflt,2), poles_it(npoles,3), temp_it, uncert_it
! logical :: isAccepted_it
!
!
!
! if (verbosity.ge.1) then
!     write(stdout,*) 'read_resampled_models: starting'
! endif
!
!
! ! Read initial (0th) model
! if (log_format.eq.'anneal') then
!     nlines_per_it = 1 + nflt
! elseif (log_format.eq.'anneal_psc') then
!     nlines_per_it = 1 + nflt
! elseif (log_format.eq.'anneal_psc_euler') then
!     nlines_per_it = 1 + nflt + npoles*3
! else
!     call usage('anneal_fem_post: no log file format named '//trim(log_format))
! endif
! do i = 1,nlines_per_it
!     read(lu_log_file,'(A)') input_line
! enddo
!
!
! isamp = 1
!
! ! Read all other models
! do it = 1,nit
!
!     ! Progress indicator
!     call progress_indicator(it,nit,'anneal_fem_post: read_resampled_models',ios)
!
!     if (isamp.gt.nresamp) then
!         exit
!     endif
!
!     ! Read the header for the iteration
!     ! All iteration header lines start with a ">" in the first character and have the format:
!     !     > label=value label=value ...
!     ! Where "label" describes a parameter, and "value" is a number or string for the parameter
!     read(lu_log_file,'(A)',end=1001,iostat=ios) input_line
!     call parse_iteration_header(input_line,itread,temp_it,obj_it,uncert_it,isAccepted_it,&
!                                 ios)
!     if (ios.ne.0) then
!         write(stderr,*) 'anneal_fem_post: error parsing iteration information from line ',trim(input_line)
!         write(stderr,*) 'iteration: ',it
!         call usage('')
!     endif
!     1001 if (ios.ne.0) then
!         write(stderr,*) 'anneal_fem_post: more iterations in header (',nit,') than found in file (',it-1,')'
!         write(stderr,*) 'Continuing with iterations that have been read'
!         ios = 0
!         nit = it-1
!         exit
!     endif
!
!     ! Read fault information defined by log file format
!     if (log_format.eq.'anneal') then
!
!         ! Slip magnitudes
!         do iflt = 1,nflt
!             read(lu_log_file,'(A)',end=1002,iostat=ios) input_line
!             read(input_line,*,end=1003,err=1003,iostat=ios) slip_it(iflt,1)
!         enddo
!
!         ! Rake angles
!         do iflt = 1,nflt
!             read(lu_log_file,'(A)',end=1002,iostat=ios) input_line
!             read(input_line,*,end=1003,err=1003,iostat=ios) rake
!             ! Convert to strike- and dip-slip
!             slip_it(iflt,2) = slip_it(iflt,1)*sin(rake*d2r)
!             slip_it(iflt,1) = slip_it(iflt,1)*cos(rake*d2r)
!         enddo
!
!
!     elseif (log_format.eq.'anneal_psc') then
!         do iflt = 1,nflt
!             read(lu_log_file,'(A)',end=1002,iostat=ios) input_line
!             read(input_line,*,end=1003,err=1003,iostat=ios) locked_it(iflt), &
!                                                             slip_it(iflt,1:2)
!         enddo
!
!
!     elseif (log_format.eq.'anneal_psc_euler') then
!         do iflt = 1,nflt
!             read(lu_log_file,'(A)',end=1002,iostat=ios) input_line
!             read(input_line,*,end=1003,err=1003,iostat=ios) locked_it(iflt), &
!                                                             slip_it(iflt,1:2)
!         enddo
!         do ipole = 1,npoles
!             do i = 1,3
!                 read(lu_log_file,'(A)',end=1004,iostat=ios) input_line
!                 read(input_line,*,end=1005,err=1005,iostat=ios) poles_it(ipole,i)
!             enddo
!         enddo
!
!     else
!         call usage('anneal_fem_post: no log file format named '//trim(log_format))
!     endif
!
!
!     if (samp(isamp).ne.it) then
!         cycle
!     endif
!
!     ! Load resampled models
!     do while (samp(isamp).eq.it)
!         obj(isamp) = obj_it
!         temp(isamp) = temp_it
!         uncert(isamp) = uncert_it
!         isAccepted(isamp) = isAccepted_it
!         do iflt = 1,nflt
!             slip(isamp,iflt,1) = slip_it(iflt,1)
!             slip(isamp,iflt,2) = slip_it(iflt,2)
!             if (log_format.eq.'anneal_psc'.or.log_format.eq.'anneal_psc_euler') then
!                 locked(isamp,iflt) = locked_it(iflt)
!             endif
!         enddo
!         do ipole = 1,npoles
!             if (log_format.eq.'anneal_psc_euler') then
!                 poles(isamp,ipole,1) = poles_it(ipole,1)
!                 poles(isamp,ipole,2) = poles_it(ipole,2)
!                 poles(isamp,ipole,3) = poles_it(ipole,3)
!             endif
!         enddo
!         isamp = isamp + 1
!         if (isamp.gt.nresamp) then
!             exit
!         endif
!     enddo
!
! enddo
!
! if (verbosity.ge.1) then
!     write(stdout,*) 'read_resampled_models: finished'
! endif
!
! 1002 if (ios.ne.0) then
!     write(stderr,*) 'anneal_fem_post: reached end of anneal log file before finished reading faults'
!     write(stderr,*) 'iteration: ',it
!     write(stderr,*) 'fault: ',iflt
!     call usage('')
! endif
! 1003 if (ios.ne.0) then
!     write(stderr,*) 'anneal_fem_post: error parsing fault info from line ',trim(input_line)
!     write(stderr,*) 'iteration: ',it
!     write(stderr,*) 'fault: ',iflt
!     call usage('')
! endif
! 1004 if (ios.ne.0) then
!     write(stderr,*) 'anneal_fem_post: reached end of anneal log file before finished reading poles'
!     write(stderr,*) 'iteration: ',it
!     write(stderr,*) 'pole: ',ipole
!     call usage('')
! endif
! 1005 if (ios.ne.0) then
!     write(stderr,*) 'anneal_fem_post: error parsing pole info from line ',trim(input_line)
!     write(stderr,*) 'iteration: ',it
!     write(stderr,*) 'pole: ',ipole
!     call usage('')
! endif
!
! return
! end subroutine
!
!
!--------------------------------------------------------------------------------------------------!


subroutine gcmdln()

! use io, only: verbosity
! use random, only: iseed
!
use anneal_fem_post, only: anneal_log_file , &
                           listParameters, &
                           noutputs, &
                           param_outputs, &
                           param_output_files, &
                           burnin_iterations, &
                           best_model_file, &
                           mean_model_file, &
                           median_model_file


implicit none

! Local variables
integer :: i, narg
character(len=512) :: tag, arg
logical :: isOutputDefined


! Initialize control parameters
isOutputDefined = .false.

anneal_log_file = ''
listParameters = .false.
burnin_iterations = 0
best_model_file = ''
mean_model_file = ''
median_model_file = ''


! Number of arguments
narg = command_argument_count()
if (narg.eq.0) then
  call usage('')
endif


! Count number of quantities to extract
i = 1
noutputs = 0
do while (i.le.narg)
    call get_command_argument(i,tag)
    if (index(tag,':marg').ne.0) then
        noutputs = noutputs + 1
    elseif (index(tag,':marginal').ne.0) then
        noutputs = noutputs + 1
    endif
    i = i + 1
enddo
allocate(param_outputs(noutputs))
allocate(param_output_files(noutputs))
param_outputs = ''
param_output_files = ''


! Read arguments
i = 1
noutputs = 0
do while (i.le.narg)

    call get_command_argument(i,tag)

    if (trim(tag).eq.'-f') then
        i = i + 1
        call get_command_argument(i,anneal_log_file)

    elseif (index(tag,':marg').ne.0) then
        noutputs = noutputs + 1
        param_outputs(noutputs) = tag(1:32)
        i = i + 1
        call get_command_argument(i,param_output_files(noutputs))
        isOutputDefined = .true.
    elseif (index(tag,':marginal').ne.0) then
        noutputs = noutputs + 1
        param_outputs(noutputs) = tag(1:32)
        i = i + 1
        call get_command_argument(i,param_output_files(noutputs))
        isOutputDefined = .true.

    elseif (trim(tag).eq.'-param:list') then
        listParameters = .true.
        isOutputDefined = .true.

    elseif (trim(tag).eq.'-nburn') then
        i = i + 1
        call get_command_argument(i,arg)
        read(arg,*) burnin_iterations

    elseif (trim(tag).eq.'-best') then
        i = i + 1
        call get_command_argument(i,best_model_file)
        isOutputDefined = .true.
    elseif (trim(tag).eq.'-mean') then
        i = i + 1
        call get_command_argument(i,mean_model_file)
        isOutputDefined = .true.
    elseif (trim(tag).eq.'-median') then
        i = i + 1
        call get_command_argument(i,median_model_file)
        isOutputDefined = .true.

    else
        call usage('anneal_fem_post: no option '//trim(tag))
    endif


    i = i + 1

enddo

if (.not.isOutputDefined) then
    call usage('anneal_fem_post: no output file defined')
endif

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine usage(str)

use io, only: stderr

implicit none

character(len=*) :: str

if (str.ne.'') then
    write(stderr,*) trim(str)
    write(stderr,*)
endif

write(stderr,*) 'Usage: anneal_fem_post -f LOG_FILE [...options...]'
write(stderr,*)
write(stderr,*) '-f LOG_FILE                     Annealing log file'
write(stderr,*) '-nburn NBURN                    Number of burn-in iterations (not printed)'
write(stderr,*) '-best FILE                      Best-fitting model'
write(stderr,*) '-mean FILE                      Mean model'
write(stderr,*) '-median FILE                    Median model'
write(stderr,*) '[TODO] -obj FILE                       Iteration and objective function'
write(stderr,*) '-<param>:marg FILE              Parameter marginal distribution'
write(stderr,*) '-param:list                     List available parameters'
write(stderr,*) '[TODO] -v VERBOSITY                    Turn on verbose mode'
write(stderr,*)
write(stderr,*) 'The option "-<param>:marg" allows the user to extract marginal information from'
write(stderr,*) 'the search. It can be repeated as many times as desired for different parameters.'
write(stderr,*) '<param> must be one of the parameters in the search (use -param:list to show them)'
write(stderr,*)
! write(stderr,*) '-flocked:subset FLT_FILE FILE   Fraction of locked faults in each iteration from a ',&
!                                                  'subset of faults'
! write(stderr,*) '-uncertainty UNCERT_FILE        Uncertainty distribution'
! write(stderr,*) '-resample NRESAMPLE             Resample models proportional to probability'
! write(stderr,*) '-seed SEED                      Set random number seed'
! write(stderr,*)

call error_exit(1)
end subroutine
