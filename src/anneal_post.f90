!--------------------------------------------------------------------------------------------------!
! ANNEAL_POST
!
! Post-process a simulated annealing search output file from FLTINV. The precise output format
! depends on the options used in the annealing search. In general, the format is:
!
!     # Header Line 1
!     # Header Line 2
!     :
!     > Iteration i Header
!     value(s)
!     value(s)
!     :
!
!--------------------------------------------------------------------------------------------------!

module anneal_post

character(len=512) :: anneal_log_file
integer :: lu_log_file
character(len=512) :: best_output_file
character(len=512) :: plocked_output_file
character(len=512) :: mean_slip_output_file
integer :: burnin_iterations

integer :: nmarg1d
integer :: marg1d_flt(1000)
character(len=64) :: slip_marg_output_file(1000)

character(len=512) :: mean_pole_output_file
integer :: nmarg_pole
integer :: marg_pole(20)
character(len=64) :: pole_marg_output_file(20)

integer :: nflocked
character(len=64) :: flocked_subset_file(10)
character(len=64) :: flocked_output_file(10)

character(len=512) :: obj_file
character(len=512) :: uncert_file

logical :: doResample
integer :: nresamp
integer, allocatable :: samp(:)

end module


!==================================================================================================!


program main

use io, only: stdout
use earth, only: pole_geo2xyz, pole_xyz2geo
use annealing, only: resample
use random, only: iseed, timeseed
use sort, only: heapsort

use anneal_post, only: anneal_log_file, &
                       lu_log_file, &
                       burnin_iterations, &
                       best_output_file, &
                       plocked_output_file, &
                       mean_slip_output_file, &
                       nmarg1d, &
                       marg1d_flt, &
                       slip_marg_output_file, &
                       mean_pole_output_file, &
                       nmarg_pole, &
                       marg_pole, &
                       pole_marg_output_file, &
                       nflocked, &
                       flocked_subset_file, &
                       flocked_output_file, &
                       obj_file, &
                       uncert_file, &
                       doResample, &
                       nresamp, &
                       samp

implicit none

! Local variables
character(len=512) :: input_line
character(len=32) :: log_format
logical :: fileExists
logical, allocatable :: isAccepted(:)
logical, allocatable :: isFltInSubset(:)
integer :: it, nit
integer :: i, ios, iflt, nflt, ipole, npoles, nsubset
integer, allocatable :: locked(:,:), plocked(:)
double precision, allocatable :: temp(:)
double precision, allocatable :: obj(:), obj_resamp(:)
double precision, allocatable :: uncert(:)
double precision :: min_obj, fraction_locked
double precision, allocatable :: slip(:,:,:), mean_ss(:), mean_ds(:)
double precision :: px, py, pz, plon, plat, pmag
double precision, allocatable :: poles(:,:,:), mean_px(:), mean_py(:), mean_pz(:)
double precision :: ddum
logical :: ldum

ios = 0

! Parse command line
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
lu_log_file = 21
open(unit=lu_log_file,file=anneal_log_file,status='old')


!----
! Read annealing log file
!----
! Parse log file header
call read_header(nit,nflt,npoles,log_format)


! Set up arrays
allocate(obj(nit))                ! DP: Objective function values
allocate(locked(nit,nflt))        ! Int: Locked (1) and unlocked (0) faults
allocate(slip(nit,nflt,2))        ! DP: Strike- and dip-slip of each fault
allocate(poles(nit,npoles,3))     ! DP: Pole coordinates and rotation rate
allocate(temp(nit))               ! DP: Temperature
allocate(uncert(nit))             ! DP: Model prediction error
allocate(isAccepted(nit))         ! Logical: Was the proposed model accepted?
write(*,*) 'anneal_post: finished allocating array memory'


! Resample ensemble if requested
if (doResample) then
    if (nresamp.eq.0) then
        nresamp = nit
    endif

    write(*,*) 'anneal_post: reading model objective functions for resampling'

    ! Read all objectives (log-posterior-probabilities)
    i = 0
    do
        read(lu_log_file,'(A)',iostat=ios) input_line
        if (ios.ne.0) then
            exit
        endif
        if (index(input_line,'>').eq.1) then
            if (i.ge.1) then
                call parse_iteration_header(input_line,it,ddum,obj(i),ddum,ldum,ios)
            endif
            i = i + 1
        endif
    enddo

    write(*,*) 'anneal_post: resampling models'
    if (iseed.eq.0) then
        iseed = timeseed()
    endif
    allocate(samp(nresamp))
    allocate(obj_resamp(nit))
    obj_resamp = obj - maxval(obj)
    do i = 1,nit
        obj_resamp(i) = exp(obj_resamp(i))
    enddo
    call resample(obj_resamp,samp,nit)
    call heapsort(samp,nit)

    ! Reset file to read again
    rewind(lu_log_file)
    call read_header(nit,nflt,npoles,log_format)

	! Read search results
    call read_resampled_models(nit,nflt,npoles,log_format,locked,obj,slip,poles,temp,uncert,isAccepted)
else
    ! Read search results
    call read_models(nit,nflt,npoles,log_format,locked,obj,slip,poles,temp,uncert,isAccepted)
endif



!----
! Print requested quantities
!----
! Best fitting model
if (best_output_file.ne.'') then
    write(*,*) 'anneal_post: working on best-fitting output file'
    ! Find minimum chi-squared model
    min_obj = 1e10
    it = 1
    do i = 1,nit
        if (abs(obj(i)).lt.min_obj) then
            it = i
            min_obj = abs(obj(i))
        endif
    enddo

    ! Write locked/unlocked and slip of best model to file
    write(stdout,*) 'best fit: ',it,obj(it)
    open(unit=22,file=best_output_file,status='unknown')
    do iflt = 1,nflt
        write(22,*) locked(it,iflt),slip(it,iflt,:)
    enddo
    close(22)

    ! Write locked/unlocked and slip of best model to file
    open(unit=22,file='pole_'//best_output_file,status='unknown')
    do ipole = 1,npoles
        write(22,*) poles(it,ipole,:)
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
        if (it.gt.burnin_iterations) then
            plocked = plocked + locked(it,:)
        endif
    enddo

    ! Write locked probability to file
    open(unit=23,file=plocked_output_file,status='unknown')
    do iflt = 1,nflt
        write(23,*) dble(plocked(iflt))/dble(nit-burnin_iterations)
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
            if (it.gt.burnin_iterations) then
                mean_ss(iflt) = mean_ss(iflt) + slip(it,iflt,1)
                mean_ds(iflt) = mean_ds(iflt) + slip(it,iflt,2)
            endif
        enddo
        mean_ss(iflt) = mean_ss(iflt)/dble(nit-burnin_iterations)
        mean_ds(iflt) = mean_ds(iflt)/dble(nit-burnin_iterations)
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
        if (marg1d_flt(i).gt.nflt) then
            call usage('anneal_post: requested slip in fault greater than number of faults')
        endif
        open(unit=25,file=slip_marg_output_file(i),status='unknown')
        do it = burnin_iterations+1,nit
            write(25,*) slip(it,marg1d_flt(i),1:2)
        enddo
        close(25)
    enddo
endif

! Fraction of interface that is locked
! TODO: SCALE LOCKED ZONES BY AREA
if (nflocked.ge.1) then
    if (flocked_subset_file(1).eq.'') then
        open(unit=26,file=flocked_output_file(1),status='unknown')
        do it = burnin_iterations+1,nit
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
            do it = burnin_iterations+1,nit
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


! Mean Euler pole (AND STANDARD DEVIATION???)
if (mean_pole_output_file.ne.'') then
    ! Compute mean rotation vector
    allocate(mean_px(npoles))
    allocate(mean_py(npoles))
    allocate(mean_pz(npoles))
    mean_px = 0.0d0
    mean_py = 0.0d0
    mean_pz = 0.0d0
    do ipole = 1,npoles
        do it = 1,nit
            call pole_geo2xyz(poles(it,ipole,1),poles(it,ipole,2),poles(it,ipole,3), &
                              px,py,pz,'sphere')
            mean_px(ipole) = mean_px(ipole) + px
            mean_py(ipole) = mean_py(ipole) + py
            mean_pz(ipole) = mean_pz(ipole) + pz
        enddo
        mean_px(ipole) = mean_px(ipole)/dble(nit)
        mean_py(ipole) = mean_py(ipole)/dble(nit)
        mean_pz(ipole) = mean_pz(ipole)/dble(nit)
        call pole_xyz2geo(mean_px(ipole),mean_py(ipole),mean_pz(ipole),plon,plat,pmag,'sphere')
        mean_px(ipole) = plon
        mean_py(ipole) = plat
        mean_pz(ipole) = pmag
    enddo

    ! Write mean Euler pole to file
    open(unit=27,file=mean_pole_output_file,status='unknown')
    do ipole = 1,npoles
        write(27,*) mean_px(ipole),mean_py(ipole),mean_pz(ipole)
    enddo
    close(27)

    deallocate(mean_px)
    deallocate(mean_py)
    deallocate(mean_pz)
endif

! Marginal distributions for Euler pole location and angular velocity
if (nmarg_pole.ge.1) then
    do i = 1,nmarg_pole
        if (marg_pole(i).gt.npoles) then
            call usage('anneal_post: requested information in pole greater than number of poles')
        endif
        open(unit=28,file=pole_marg_output_file(i),status='unknown')
        do it = burnin_iterations+1,nit
            write(28,*) poles(it,marg_pole(i),1:3)
        enddo
        close(28)
    enddo
endif


! Annealing objective versus iteration
if (obj_file.ne.'') then
    open(unit=29,file=obj_file,status='unknown')
    do i = 1,nit
        write(29,*) i,-obj(i),isAccepted(i)
    enddo
    close(29)
endif


! Uncertainty histogram
if (uncert_file.ne.'') then
    open(unit=30,file=uncert_file,status='unknown')
    do i = 1,nit
        write(30,*) uncert(i)
    enddo
    close(30)
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
close(lu_log_file)

end

!--------------------------------------------------------------------------------------------------!

subroutine read_models(nit,nflt,npoles,log_format,locked,obj,slip,poles,temp,uncert, &
                       isAccepted)

use io, only: stdout, stderr
use trig, only: d2r

use anneal_post, only: lu_log_file

implicit none

! Arguments
integer :: nit, nflt, npoles
character(len=*) :: log_format
integer :: locked(nit,nflt)
double precision :: obj(nit), slip(nit,nflt,2), poles(nit,npoles,3), temp(nit), uncert(nit)
logical :: isAccepted(nit)

! Local variables
integer :: i, ios, it, itread, iflt, ipole, nlines_per_it
character(len=512) :: input_line
double precision :: rake


if (log_format.eq.'anneal') then
    nlines_per_it = 1 + nflt
elseif (log_format.eq.'anneal_psc') then
    nlines_per_it = 1 + nflt
elseif (log_format.eq.'anneal_psc_euler') then
    nlines_per_it = 1 + nflt + npoles*3
else
    call usage('anneal_post: no log file format named '//trim(log_format))
endif
do i = 1,nlines_per_it
    read(lu_log_file,'(A)') input_line
enddo


! Read all other models
do it = 1,nit

    ! Progress indicator
    if (nit.ge.100) then
        if (mod(it,nit/100).eq.1.and.nit.ge.100) then
            write(*,'(A,I5,A,A)',advance='no') 'anneal_post progress: ',100*it/nit,'%',char(13)
        endif
    else
        write(*,'(A,I5,A,I5,A)',advance='no') 'anneal_post progress: ',it,'of',nit,char(13)
    endif


    ! Read the header for each iteration
    ! All iteration header lines start with a ">" in the first character and have the format:
    !     > label=value label=value ...
    ! Where "label" describes a parameter, and "value" is a number or string for the parameter
    read(lu_log_file,'(A)',end=1001,iostat=ios) input_line
    if (it.gt.0) then
        call parse_iteration_header(input_line,itread,temp(it),obj(it),uncert(it),isAccepted(it),ios)
    else
        itread = 0
    endif
    if (ios.ne.0) then
        write(stderr,*) 'anneal_post: error parsing iteration information from line ',trim(input_line)
        write(stderr,*) 'iteration: ',it
        call usage('')
    endif
    if (itread.ne.it) then
        call usage('anneal_post: iterations in log file are not in sequential order')
    endif
    1001 if (ios.ne.0) then
        write(stderr,*) 'anneal_post: more iterations indicated (',nit,') than in file (',it-1,')'
        call usage('')
    endif


    ! Read fault information defined by log file format
    if (log_format.eq.'anneal') then

        ! Slip magnitudes
        do iflt = 1,nflt
            read(lu_log_file,'(A)',end=1002,iostat=ios) input_line
            if (it.gt.0) then
                read(input_line,*,end=1003,err=1003,iostat=ios) slip(it,iflt,1)
            endif
        enddo

        ! Rake angles
        do iflt = 1,nflt
            read(lu_log_file,'(A)',end=1002,iostat=ios) input_line
            if (it.gt.0) then
                read(input_line,*,end=1003,err=1003,iostat=ios) rake
                ! Convert to strike- and dip-slip
                slip(it,iflt,2) = slip(it,iflt,1)*sin(rake*d2r)
                slip(it,iflt,1) = slip(it,iflt,1)*cos(rake*d2r)
            endif
        enddo


    elseif (log_format.eq.'anneal_psc') then
        do iflt = 1,nflt
            read(lu_log_file,'(A)',end=1002,iostat=ios) input_line
            if (it.gt.0) then
                read(input_line,*,end=1003,err=1003,iostat=ios) locked(it,iflt), slip(it,iflt,1:2)
            endif
        enddo


    elseif (log_format.eq.'anneal_psc_euler') then
        do iflt = 1,nflt
            read(lu_log_file,'(A)',end=1002,iostat=ios) input_line
            if (it.gt.0) then
                read(input_line,*,end=1003,err=1003,iostat=ios) locked(it,iflt), slip(it,iflt,1:2)
            endif
        enddo
        do ipole = 1,npoles
            do i = 1,3
                read(lu_log_file,'(A)',end=1004,iostat=ios) input_line
                if (it.gt.0) then
                    read(input_line,*,end=1005,err=1005,iostat=ios) poles(it,ipole,i)
                endif
            enddo
        enddo

    else
        call usage('anneal_post: no log file format named '//trim(log_format))
    endif

enddo
write(stdout,*) 'anneal_post: finished reading log file'
1002 if (ios.ne.0) then
    write(stderr,*) 'anneal_post: reached end of anneal log file before finished reading faults'
    write(stderr,*) 'iteration: ',it
    write(stderr,*) 'fault: ',iflt
    call usage('')
endif
1003 if (ios.ne.0) then
    write(stderr,*) 'anneal_post: error parsing fault info from line ',trim(input_line)
    write(stderr,*) 'iteration: ',it
    write(stderr,*) 'fault: ',iflt
    call usage('')
endif
1004 if (ios.ne.0) then
    write(stderr,*) 'anneal_post: reached end of anneal log file before finished reading poles'
    write(stderr,*) 'iteration: ',it
    write(stderr,*) 'pole: ',ipole
    call usage('')
endif
1005 if (ios.ne.0) then
    write(stderr,*) 'anneal_post: error parsing pole info from line ',trim(input_line)
    write(stderr,*) 'iteration: ',it
    write(stderr,*) 'pole: ',ipole
    call usage('')
endif
return
end subroutine


!--------------------------------------------------------------------------------------------------!


subroutine read_resampled_models(nit,nflt,npoles,log_format,locked,obj,slip,poles,temp, &
                                 uncert,isAccepted)

use io, only: stdout, stderr
use trig, only: d2r

use anneal_post, only: lu_log_file, &
                       nresamp, &
                       samp

implicit none

! Arguments
integer :: nit, nflt, npoles
character(len=*) :: log_format
integer :: locked(nit,nflt)
double precision :: obj(nit), slip(nit,nflt,2), poles(nit,npoles,3), temp(nit), uncert(nit)
logical :: isAccepted(nit)

! Local variables
integer :: i, ios, it, itread, iflt, ipole, nlines_per_it, isamp
character(len=512) :: input_line
double precision :: rake
integer :: locked_it(nflt)
double precision :: obj_it, slip_it(nflt,2), poles_it(npoles,3), temp_it, uncert_it
logical :: isAccepted_it



! Read initial (0th) model
if (log_format.eq.'anneal') then
    nlines_per_it = 1 + nflt
elseif (log_format.eq.'anneal_psc') then
    nlines_per_it = 1 + nflt
elseif (log_format.eq.'anneal_psc_euler') then
    nlines_per_it = 1 + nflt + npoles*3
else
    call usage('anneal_post: no log file format named '//trim(log_format))
endif
do i = 1,nlines_per_it
    read(lu_log_file,'(A)') input_line
enddo


isamp = 1

! Read all other models
do it = 1,nit

    ! Progress indicator
    if (nit.ge.100) then
        if (mod(it,nit/100).eq.1.and.nit.ge.100) then
            write(*,'(A,I5,A,A)',advance='no') 'anneal_post progress: ',100*it/nit,'%',char(13)
        endif
    else
        write(*,'(A,I5,A,I5,A)',advance='no') 'anneal_post progress: ',it,'of',nit,char(13)
    endif

    if (isamp.gt.nresamp) then
        exit
    endif

    ! Read the header for the iteration
    ! All iteration header lines start with a ">" in the first character and have the format:
    !     > label=value label=value ...
    ! Where "label" describes a parameter, and "value" is a number or string for the parameter
    read(lu_log_file,'(A)',end=1001,iostat=ios) input_line
    call parse_iteration_header(input_line,itread,temp_it,obj_it,uncert_it,isAccepted_it,&
                                ios)
    if (ios.ne.0) then
        write(stderr,*) 'anneal_post: error parsing iteration information from line ',trim(input_line)
        write(stderr,*) 'iteration: ',it
        call usage('')
    endif
    1001 if (ios.ne.0) then
        write(stderr,*) 'anneal_post: more iterations indicated (',nit,') than in file (',it-1,')'
        call usage('')
    endif

    ! Read fault information defined by log file format
    if (log_format.eq.'anneal') then

        ! Slip magnitudes
        do iflt = 1,nflt
            read(lu_log_file,'(A)',end=1002,iostat=ios) input_line
            read(input_line,*,end=1003,err=1003,iostat=ios) slip_it(iflt,1)
        enddo

        ! Rake angles
        do iflt = 1,nflt
            read(lu_log_file,'(A)',end=1002,iostat=ios) input_line
            read(input_line,*,end=1003,err=1003,iostat=ios) rake
            ! Convert to strike- and dip-slip
            slip_it(iflt,2) = slip_it(iflt,1)*sin(rake*d2r)
            slip_it(iflt,1) = slip_it(iflt,1)*cos(rake*d2r)
        enddo


    elseif (log_format.eq.'anneal_psc') then
        do iflt = 1,nflt
            read(lu_log_file,'(A)',end=1002,iostat=ios) input_line
            read(input_line,*,end=1003,err=1003,iostat=ios) locked_it(iflt), &
                                                            slip_it(iflt,1:2)
        enddo


    elseif (log_format.eq.'anneal_psc_euler') then
        do iflt = 1,nflt
            read(lu_log_file,'(A)',end=1002,iostat=ios) input_line
            read(input_line,*,end=1003,err=1003,iostat=ios) locked_it(iflt), &
                                                            slip_it(iflt,1:2)
        enddo
        do ipole = 1,npoles
            do i = 1,3
                read(lu_log_file,'(A)',end=1004,iostat=ios) input_line
                read(input_line,*,end=1005,err=1005,iostat=ios) poles_it(ipole,i)
            enddo
        enddo

    else
        call usage('anneal_post: no log file format named '//trim(log_format))
    endif


    if (samp(isamp).ne.it) then
        cycle
    endif

    ! Load resampled models
    do while (samp(isamp).eq.it)
        obj(isamp) = obj_it
        temp(isamp) = temp_it
        uncert(isamp) = uncert_it
        isAccepted(isamp) = isAccepted_it
        do iflt = 1,nflt
            slip(isamp,iflt,1) = slip_it(iflt,1)
            slip(isamp,iflt,2) = slip_it(iflt,2)
            if (log_format.eq.'anneal_psc'.or.log_format.eq.'anneal_psc_euler') then
                locked(isamp,iflt) = locked_it(iflt)
            endif
        enddo
        do ipole = 1,npoles
            if (log_format.eq.'anneal_psc_euler') then
                poles(isamp,ipole,1) = poles_it(ipole,1)
                poles(isamp,ipole,2) = poles_it(ipole,2)
                poles(isamp,ipole,3) = poles_it(ipole,3)
            endif
        enddo
        isamp = isamp + 1
        if (isamp.gt.nresamp) then
            exit
        endif
    enddo

enddo

write(stdout,*) 'anneal_post: finished reading log file'
1002 if (ios.ne.0) then
    write(stderr,*) 'anneal_post: reached end of anneal log file before finished reading faults'
    write(stderr,*) 'iteration: ',it
    write(stderr,*) 'fault: ',iflt
    call usage('')
endif
1003 if (ios.ne.0) then
    write(stderr,*) 'anneal_post: error parsing fault info from line ',trim(input_line)
    write(stderr,*) 'iteration: ',it
    write(stderr,*) 'fault: ',iflt
    call usage('')
endif
1004 if (ios.ne.0) then
    write(stderr,*) 'anneal_post: reached end of anneal log file before finished reading poles'
    write(stderr,*) 'iteration: ',it
    write(stderr,*) 'pole: ',ipole
    call usage('')
endif
1005 if (ios.ne.0) then
    write(stderr,*) 'anneal_post: error parsing pole info from line ',trim(input_line)
    write(stderr,*) 'iteration: ',it
    write(stderr,*) 'pole: ',ipole
    call usage('')
endif
return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine parse_iteration_header(input_line,it,temp,obj,uncert,isAccepted,ierr)
!----
! Parse the iteration header, which is in the format:
! > Iteration=N Temperature=T Objective=F Model_Uncertainty=U accepted|rejected
!----

use io, only: stderr

implicit none

! I/O variables
character(len=*) :: input_line
integer :: it, ierr
double precision :: temp, obj, uncert
logical :: isAccepted

! Local variables
integer :: i, ios
character(len=64) :: ch


! Initialize values
ierr = 0
it = 0
temp = 0.0d0
obj = 0.0d0
uncert = -1.0d0
isAccepted = .false.


! Make sure line is properly labeled as an iteration header
i = index(input_line,'>')
if (i.ne.1) then
    write(stderr,*) 'parse_iteration_header: iteration header line does not start with ">"'
    call usage('Offending line: '//trim(input_line))
else
    input_line(1:i) = ''
endif


! Parse the iteration header
do
    ! Read the first entry
    read(input_line,*,iostat=ios) ch
    if (ios.ne.0) then
        exit
    endif

    if (ch.eq.'accepted') then
        isAccepted = .true.
    elseif (ch.eq.'rejected') then
        isAccepted = .false.
    else
        i = index(ch,'=')
        if (ch(1:i-1).eq.'Iteration') then
            read(ch(i+1:len(trim(ch))),*) it
        elseif (ch(1:i-1).eq.'Temperature') then
            read(ch(i+1:len(trim(ch))),*) temp
        elseif (ch(1:i-1).eq.'Objective') then
            read(ch(i+1:len(trim(ch))),*) obj
        elseif (ch(1:i-1).eq.'Model_Uncertainty') then
            read(ch(i+1:len(trim(ch))),*) uncert
        else
            call usage('parse_iteration_header: could not parse this parameter: '//ch(i:i-1))
        endif
    endif

    ! Remove the entry after done parsing
    i = index(input_line,trim(ch))
    input_line(1:i+len(trim(ch))) = ''
enddo

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine read_header(nit,nflt,npoles,log_format)

use io, only: stderr, stdout

use anneal_post, only: lu_log_file, &
                       burnin_iterations

implicit none

! Arguments
integer :: nit, nflt, npoles
character(len=*) ::  log_format

! Local variables
logical :: continueReadingHeader
integer :: i, ios
character(len=512) :: input_line


! Set search variables before reading
nit = 0
nflt = 0
npoles = 0
log_format = ''


! Initialize to continue reading header
continueReadingHeader = .true.


do while (continueReadingHeader)

    read(lu_log_file,'(A)',end=4001,iostat=ios) input_line

    ! All header lines start with a "#" in the first character and have the format:
    !     # label1=value1 label2=value2 ...
    ! Where "label" describes a parameter, and "value" is a number or string for the parameter

    if (index(input_line,'#').ne.1) then
        ! Not a header line; finished parsing header
        continueReadingHeader = .false.
        backspace(lu_log_file)

    else

        ! Found header line; parse parameter and value
        if (index(input_line,'niterations').ne.0) then
            i = index(input_line,'=')
            input_line(1:i) = ''
            read(input_line,*) nit

        elseif (index(input_line,'nfaults').ne.0) then
            i = index(input_line,'=')
            input_line(1:i) = ''
            read(input_line,*) nflt

        elseif (index(input_line,'npoles').ne.0) then
            i = index(input_line,'=')
            input_line(1:i) = ''
            read(input_line,*) npoles

        elseif (index(input_line,'format').ne.0) then
            i = index(input_line,'=')
            input_line(1:i) = ''
            read(input_line,*) log_format

        else
            write(stderr,*) 'anneal_log: read header line:',trim(input_line)
            call usage('But could not parse a parameter')
        endif

    endif

enddo

4001 if (ios.ne.0) then
    call usage('anneal_post: reached end of file in header section')
endif
if (nit.le.burnin_iterations) then
    call usage('anneal_post: number of burn-in iterations is greater than total iterations')
endif

if (nit.eq.0) then
    call usage('anneal_post: read zero iterations in log file')
endif
if (nflt.eq.0) then
    call usage('anneal_post: read zero faults in log file')
endif

write(stdout,*) 'anneal_post: finished parsing header'
write(stdout,*) 'niterations=',nit
write(stdout,*) 'nfaults=',nflt
write(stdout,*) 'npoles=',npoles
write(stdout,*) 'log_format= ',log_format

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine gcmdln()

use random, only: iseed

use anneal_post, only: anneal_log_file, &
                       burnin_iterations, &
                       best_output_file, &
                       plocked_output_file, &
                       mean_slip_output_file, &
                       nmarg1d, &
                       marg1d_flt, &
                       slip_marg_output_file, &
                       mean_pole_output_file, &
                       nmarg_pole, &
                       marg_pole, &
                       pole_marg_output_file, &
                       nflocked, &
                       flocked_subset_file, &
                       flocked_output_file, &
                       obj_file, &
                       uncert_file, &
                       doResample, &
                       nresamp

implicit none

! Local variables
integer :: i, narg, ios
character(len=512) :: tag, arg
logical :: isOutputDefined


! Initialize control parameters
isOutputDefined = .false.

anneal_log_file = ''
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

mean_pole_output_file = ''
nmarg_pole = 0
marg_pole = 0
pole_marg_output_file = ''

obj_file = ''
uncert_file = ''

doResample = .false.

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

    elseif (trim(tag).eq.'-nburn') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) burnin_iterations

    elseif (trim(tag).eq.'-best') then
        i = i + 1
        call get_command_argument(i,best_output_file)
        isOutputDefined = .true.

    elseif (trim(tag).eq.'-plocked') then
        i = i + 1
        call get_command_argument(i,plocked_output_file)
        isOutputDefined = .true.

    elseif (trim(tag).eq.'-slip:mean') then
        i = i + 1
        call get_command_argument(i,mean_slip_output_file)
        isOutputDefined = .true.

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
        isOutputDefined = .true.

    elseif (trim(tag).eq.'-pole:mean') then
        i = i + 1
        call get_command_argument(i,mean_pole_output_file)
        isOutputDefined = .true.

    elseif (trim(tag).eq.'-pole:marg') then
        nmarg_pole = nmarg_pole + 1
        if (nmarg_pole.gt.20) then
            call usage('anneal_post: requested more than 1000 pole marginals')
        endif
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) marg_pole(nmarg_pole)
        i = i + 1
        call get_command_argument(i,tag)
        if (len_trim(tag).gt.64) then
            call usage('anneal_post: rotation pole marginal file must be 64 or less characters')
        endif
        pole_marg_output_file(nmarg_pole) = tag(1:64)
        isOutputDefined = .true.

    elseif (trim(tag).eq.'-flocked') then
        nflocked = 1
        i = i + 1
        call get_command_argument(i,flocked_output_file(1))
        isOutputDefined = .true.

    elseif (trim(tag).eq.'-flocked:subset') then
        nflocked = nflocked + 1
        if (nmarg1d.gt.10) then
            call usage('anneal_post: requested more than 10 locked fraction outputs')
        endif
        i = i + 1
        call get_command_argument(i,flocked_subset_file(nflocked))
        i = i + 1
        call get_command_argument(i,flocked_output_file(nflocked))
        isOutputDefined = .true.

    elseif (tag.eq.'-obj') then
        i = i + 1
        call get_command_argument(i,obj_file)
        isOutputDefined = .true.

    elseif (tag.eq.'-uncertainty') then
        i = i + 1
        call get_command_argument(i,uncert_file)
        isOutputDefined = .true.

    elseif (tag.eq.'-resample') then
        doResample = .true.
        nresamp = 0
        i = i + 1
        if (i.gt.narg) then
            i = i - 1
        else
            call get_command_argument(i,arg)
            read(arg,*,iostat=ios) iseed
            if (ios.ne.0) then
                iseed = 0
                i = i - 1
            endif
        endif

    else
        call usage('anneal_post: no option '//trim(tag))
    endif

    i = i + 1

enddo

if (.not.isOutputDefined) then
    call usage('anneal_post: no output file defined')
endif

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
write(stderr,*) '       [-slip:marg IFLT FILE] [-pole:mean POLE_FILE] [-pole:marg IPOL FILE]'
write(stderr,*) '       [-flocked[:subset FLT_FILE] FILE] [-obj OBJ_FILE]'
write(stderr,*)
write(stderr,*) '-f LOG_FILE                     Annealing log file from fltinv'
write(stderr,*) '-nit NITER                      Number of iterations in search'
write(stderr,*) '-nburn NBURN                    Number of burn-in iterations in search'
write(stderr,*) '-best BEST_FILE                 Best fitting model (locked ss ds): also pole_BEST_FILE'
write(stderr,*) '-plocked PROB_FILE              Probability locked (plocked)'
write(stderr,*) '-slip:mean SLIP_FILE            Mean slip in each fault'
write(stderr,*) '-slip:marg IFLT FILE            1-D marginal for IFLT slip (repeat for other faults)'
write(stderr,*) '-pole:mean POLE_FILE            Mean Euler poles'
write(stderr,*) '-pole:marg IPOL FILE            Marginals for IPOL rotation pole'
write(stderr,*) '-flocked FILE                   Fraction of locked faults in each iteration'
write(stderr,*) '-flocked:subset FLT_FILE FILE   Fraction of locked faults in each iteration from a ',&
                                                 'subset of faults'
write(stderr,*) '-obj OBJ_FILE                   Iteration and objective function'
write(stderr,*) '-uncertainty UNCERT_FILE        Uncertainty distribution'
write(stderr,*) '-resample [SEED]                Resample models proportional to probability (objective)'
write(stderr,*)

call error_exit(1)
end subroutine
