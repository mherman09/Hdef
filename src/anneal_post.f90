module anneal_post

character(len=512) :: anneal_log_file
character(len=512) :: best_output_file
character(len=512) :: plocked_output_file
character(len=512) :: mean_slip_output_file
integer :: search_iterations
integer :: burnin_iterations

integer :: nmarg1d
integer :: marg1d_flt(1000)
character(len=64) :: slip_marg_output_file(1000)

integer :: nflocked
character(len=64) :: flocked_subset_file(10)
character(len=64) :: flocked_output_file(10)
end module


!==================================================================================================!


program main

use io, only: stdout, stderr

use anneal_post, only: anneal_log_file, &
                       search_iterations, &
                       burnin_iterations, &
                       best_output_file, &
                       plocked_output_file, &
                       mean_slip_output_file, &
                       nmarg1d, &
                       marg1d_flt, &
                       slip_marg_output_file, &
                       nflocked, &
                       flocked_subset_file, &
                       flocked_output_file

implicit none

! Local variables
character(len=32) :: ch32
character(len=512) :: input_line
logical :: fileExists
logical, allocatable :: isFltInSubset(:)
integer :: i, ios, iflt, nflt, it, nit, nburnin, nsubset
integer, allocatable :: locked(:,:), plocked(:)
double precision :: temp0, obj0, min_chi2, fraction_locked
double precision, allocatable :: obj(:), slip(:,:,:), mean_ss(:), mean_ds(:)


ios = 0

call gcmdln()

!----
! Open annealing log file
!----
! Check that file exists
inquire(file=anneal_log_file,exist=fileExists)
if (.not.fileExists) then
    call usage('anneal_post: no anneal log file found named '//trim(anneal_log_file))
endif

! Open the log file
open(unit=21,file=anneal_log_file,status='old')


!----
! Read annealing log file
!----
! Number of iterations, burn-in period are set on the command line
nit = search_iterations
nburnin = burnin_iterations
if (nit.lt.1) then
    call usage('anneal_post: number of iterations must be 1 or greater')
endif
if (nit.le.nburnin) then
    call usage('anneal_post: number of iterations must be greater than burn-in iterations')
endif

! Count number of faults
nflt = 0
read(21,'(A)') input_line
read(input_line,*) ch32,it,ch32,temp0,ch32,obj0
do
    read(21,'(A)') input_line
    read(input_line,*) ch32
    if (ch32.eq.'Iteration') then
        exit
    else
        nflt = nflt + 1
    endif
enddo
rewind(21)

! Set up arrays
allocate(locked(nit,nflt))        ! Int array with locked (1) and unlocked (0) faults for each iteration
allocate(obj(nit))                ! DP array with objective function values for each iteration
allocate(slip(nit,nflt,2))        ! DP array with strike- and dip-slip of each fault for each iteration

! Read search results
it = 0
do
    if (mod(it,nit/100).eq.1) then
        write(*,'(A,I5,A,A)',advance='no') 'anneal_post progress: ',100*it/nit,'%',char(13)
    endif

    read(21,'(A)',end=1001,iostat=ios) input_line
    if (it.gt.nit) then
        write(stderr,*) 'anneal_post: WARNING: you requested fewer iterations (',nit,') than are in log file'
        exit
    endif
    if (it.gt.0) then
        read(input_line,*,end=1002,err=1002,iostat=ios) ch32,it,ch32,temp0,ch32,obj(it)
    endif
    do iflt = 1,nflt
        read(21,'(A)',end=1003,iostat=ios) input_line
        if (it.gt.0) then
            read(input_line,*,end=1004,err=1004,iostat=ios) locked(it,iflt),slip(it,iflt,1:2)
        endif
    enddo
    it = it + 1
enddo
ios = 0
write(*,*) 'anneal_post progress: ',100*it/nit,'%'

1001 if (it.lt.nit) then
    write(stderr,*) 'anneal_post: you requested more iterations (',nit,') than are in log file'
    call usage('')
else
    ios = 0
endif
1002 if (ios.ne.0) then
    write(stderr,*) 'anneal_post: error parsing iteration information from line ',trim(input_line)
    write(stderr,*) 'iteration: ',it
    call usage('')
endif
1003 if (ios.ne.0) then
    write(stderr,*) 'anneal_post: reached end of anneal log file before faults finished reading'
    write(stderr,*) 'iteration: ',it
    write(stderr,*) 'fault: ',iflt
    call usage('')
endif
1004 if (ios.ne.0) then
    write(stderr,*) 'anneal_post: error parsing fault info from line ',trim(input_line)
    write(stderr,*) 'iteration: ',it
    write(stderr,*) 'fault: ',iflt
    call usage('')
endif


!----
! Print requested quantities
!----
! Best fitting model
if (best_output_file.ne.'') then
    ! Find minimum chi-squared model
    min_chi2 = 1e10
    it = 1
    do i = 1,nit
        if (abs(obj(i)).lt.min_chi2) then
            it = i
            min_chi2 = abs(obj(i))
        endif
    enddo

    ! Write locked/unlocked and slip of best model to file
    write(stdout,*) 'best fit: ',it,obj(it)
    open(unit=22,file=best_output_file,status='unknown')
    do iflt = 1,nflt
        write(22,*) locked(it,iflt),slip(it,iflt,:)
    enddo
    close(22)
endif

! Probability locked
if (plocked_output_file.ne.'') then
    ! Allocate memory to array
    allocate(plocked(nflt))

    ! Calculate number of locked instances in search
    plocked = 0
    do it = 1,nit
        if (it.gt.nburnin) then
            plocked = plocked + locked(it,:)
        endif
    enddo

    ! Write locked probability to file
    open(unit=23,file=plocked_output_file,status='unknown')
    do iflt = 1,nflt
        write(23,*) dble(plocked(iflt))/dble(nit-nburnin)
    enddo
    close(23)
    deallocate(plocked)
endif

! Mean slip (AND STANDARD DEVIATION???)
if (mean_slip_output_file.ne.'') then
    ! Compute mean slip
    allocate(mean_ss(nflt))
    allocate(mean_ds(nflt))
    mean_ss = 0.0d0
    mean_ds = 0.0d0
    do iflt = 1,nflt
        do it = 1,nit
            mean_ss(iflt) = mean_ss(iflt) + slip(it,iflt,1)
            mean_ds(iflt) = mean_ds(iflt) + slip(it,iflt,2)
        enddo
        mean_ss(iflt) = mean_ss(iflt)/dble(nit)
        mean_ds(iflt) = mean_ds(iflt)/dble(nit)
    enddo

    ! Write mean slip to file
    open(unit=24,file=mean_slip_output_file,status='unknown')
    do iflt = 1,nflt
        write(24,*) mean_ss(iflt),mean_ds(iflt)
    enddo
    close(24)

    deallocate(mean_ss)
    deallocate(mean_ds)
endif

! 1-D marginal distributions for slip components
if (nmarg1d.ge.1) then
    do i = 1,nmarg1d
        open(unit=25,file=slip_marg_output_file(i),status='unknown')
        do it = nburnin+1,nit
            write(25,*) slip(it,marg1d_flt(i),1:2)
        enddo
        close(25)
    enddo
endif

! Fraction of interface that is locked
if (nflocked.ge.1) then
    if (flocked_subset_file(1).eq.'') then
        open(unit=26,file=flocked_output_file(1),status='unknown')
        do it = nburnin+1,nit
            fraction_locked = 0.0d0
            do iflt = 1,nflt
                fraction_locked = fraction_locked + dble(locked(it,iflt))
            enddo
            write(26,*) fraction_locked/dble(nflt)
        enddo
        close(26)
    else
        allocate(isFltInSubset(nflt))

        do i = 1,nflocked
            nsubset = 0
            isFLtInSubset = .false.
            open(unit=27,file=flocked_subset_file(i),status='old')
            do
                read(27,*,iostat=ios) iflt
                if (ios.ne.0) then
                    exit
                endif
                if (iflt.lt.1) then
                    call usage('anneal_post: specified locked fault index less than 1')
                elseif (iflt.gt.nflt) then
                    call usage('anneal_post: specified locked fault index greater than nflt')
                endif
                if (.not.isFltInSubset(iflt)) then
                    nsubset = nsubset + 1
                endif
                isFltInSubset(iflt) = .true.
            enddo
            close(27)

            open(unit=26,file=flocked_output_file(i),status='unknown')
            do it = nburnin+1,nit
                fraction_locked = 0.0d0
                do iflt = 1,nflt
                    if (isFltInSubset(iflt)) then
                        fraction_locked = fraction_locked + dble(locked(it,iflt))
                    endif
                enddo
                write(26,*) fraction_locked/dble(nsubset)
            enddo
            close(26)
        enddo

        deallocate(isFltInSubset)
    endif
endif

!----
!       CLEAN UP
!----
! Deallocate arrays
if (allocated(locked)) then
    deallocate(locked)
endif
if (allocated(obj)) then
    deallocate(obj)
endif
if (allocated(slip)) then
    deallocate(slip)
endif

! Close annealing log file
close(21)

end

!--------------------------------------------------------------------------------------------------!

subroutine gcmdln()

use anneal_post, only: anneal_log_file, &
                       search_iterations, &
                       burnin_iterations, &
                       best_output_file, &
                       plocked_output_file, &
                       mean_slip_output_file, &
                       nmarg1d, &
                       marg1d_flt, &
                       slip_marg_output_file, &
                       nflocked, &
                       flocked_subset_file, &
                       flocked_output_file

implicit none

! Local variables
integer :: i, narg
character(len=512) :: tag

! Initialize control parameters
anneal_log_file = ''
search_iterations = 0
burnin_iterations = 0
best_output_file = ''
plocked_output_file = ''
mean_slip_output_file = ''
nmarg1d = 0
marg1d_flt = 0
slip_marg_output_file = ''
nflocked = 0
flocked_subset_file = ''
flocked_output_file = ''

! Number of arguments
narg = command_argument_count()
if (narg.eq.0) then
  call usage('')
endif

! Read arguments
i = 1
do while (i.le.narg)

    call get_command_argument(i,tag)

    if (trim(tag).eq.'-f') then
        i = i + 1
        call get_command_argument(i,anneal_log_file)

    elseif (trim(tag).eq.'-nit') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) search_iterations
    elseif (trim(tag).eq.'-nburn') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) burnin_iterations

    elseif (trim(tag).eq.'-best') then
        i = i + 1
        call get_command_argument(i,best_output_file)

    elseif (trim(tag).eq.'-plocked') then
        i = i + 1
        call get_command_argument(i,plocked_output_file)

    elseif (trim(tag).eq.'-slip:mean') then
        i = i + 1
        call get_command_argument(i,mean_slip_output_file)

    elseif (trim(tag).eq.'-slip:marg') then
        nmarg1d = nmarg1d + 1
        if (nmarg1d.gt.1000) then
            call usage('anneal_post: requested more than 1000 slip marginals')
        endif
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) marg1d_flt(nmarg1d)
        i = i + 1
        call get_command_argument(i,tag)
        if (len_trim(tag).gt.64) then
            call usage('anneal_post: slip marginal file must be 64 or less characters')
        endif
        slip_marg_output_file(nmarg1d) = tag(1:64)

    elseif (trim(tag).eq.'-flocked') then
        nflocked = 1
        i = i + 1
        call get_command_argument(i,flocked_output_file(1))
    elseif (trim(tag).eq.'-flocked:subset') then
        nflocked = nflocked + 1
        if (nmarg1d.gt.10) then
            call usage('anneal_post: requested more than 10 locked fraction outputs')
        endif
        i = i + 1
        call get_command_argument(i,flocked_subset_file(nflocked))
        i = i + 1
        call get_command_argument(i,flocked_output_file(nflocked))

    else
        call usage('anneal_post: no option '//trim(tag))
    endif

    i = i + 1

enddo

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine usage(str)

implicit none

integer, parameter :: stderr = 0
character(len=*) :: str

if (str.ne.'') then
    write(stderr,*) trim(str)
    write(stderr,*)
endif

write(stderr,*) 'Usage: anneal_post -f LOG_FILE -nit NITER [-nburn NBURN]'
write(stderr,*) '       [-best BEST_FILE] [-plocked PROB_FILE] [-slip:mean SLIP_FILE]'
write(stderr,*) '       [-slip:marg IFLT FILE] [-flocked[:subset FLT_FILE] FILE]'
write(stderr,*)
write(stderr,*) '-f LOG_FILE                     Annealing log file from fltinv'
write(stderr,*) '-nit NITER                      Number of iterations in search'
write(stderr,*) '-nburn NBURN                    Number of burn-in iterations in search'
write(stderr,*) '-best BEST_FILE                 Best fitting model (locked ss ds)'
write(stderr,*) '-plocked PROB_FILE              Probability locked (plocked)'
write(stderr,*) '-slip:mean SLIP_FILE            Mean slip in each fault'
write(stderr,*) '-slip:marg IFLT FILE            1-D marginal for IFLT slip (repeat for other faults)'
write(stderr,*) '-flocked FILE                   Fraction of locked faults in each iteration'
write(stderr,*) '-flocked:subset FLT_FILE FILE   Fraction of locked faults in each iteration from a ',&
                                                 'subset of faults'
write(stderr,*)

stop
end subroutine
