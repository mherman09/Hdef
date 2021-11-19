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

character(len=512) :: anneal_log_file             ! Name of annealing search log file
integer :: lu_log_file                            ! Logical unit of log file for interaction

integer :: nparam                                 ! Number of model parameters
character(len=32), allocatable :: param_names(:)  ! Parameter names
! character(len=512) :: best_output_file
! integer :: burnin_iterations
!
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

use io, only: stdout, verbosity
! use earth, only: pole_geo2xyz, pole_xyz2geo
! use annealing, only: resample
! use random, only: iseed, timeseed
! use sort, only: heapsort

use anneal_fem_post, only: anneal_log_file, &
                           lu_log_file
!                        burnin_iterations, &
!                        best_output_file, &
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
integer :: ios
logical :: fileExists
! character(len=512) :: input_line
! character(len=32) :: log_format
! logical, allocatable :: isAccepted(:)
! logical, allocatable :: isFltInSubset(:)
! integer :: it, nit
! integer :: i, iflt, nflt, ipole, npoles, nsubset
! integer, allocatable :: samp_tmp(:)
! integer, allocatable :: locked(:,:), plocked(:)
! double precision, allocatable :: temp(:)
! double precision, allocatable :: obj(:), obj_resamp(:)
! double precision :: max_obj
! double precision, allocatable :: uncert(:)
! double precision :: min_obj, fraction_locked
! double precision, allocatable :: slip(:,:,:), mean_ss(:), mean_ds(:)
! double precision :: px, py, pz, plon, plat, pmag
! double precision, allocatable :: poles(:,:,:), mean_px(:), mean_py(:), mean_pz(:)
! double precision :: ddum
! integer :: idum
! logical :: ldum

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
! Count the number of parameters and set their names in parameter name array
call set_parameters()


! ! Set up arrays
! allocate(obj(nit))                ! DP: Objective function values
! allocate(locked(nit,nflt))        ! Int: Locked (1) and unlocked (0) faults
! allocate(slip(nit,nflt,2))        ! DP: Strike- and dip-slip of each fault
! allocate(poles(nit,npoles,3))     ! DP: Pole coordinates and rotation rate
! allocate(temp(nit))               ! DP: Temperature
! allocate(uncert(nit))             ! DP: Model prediction error
! allocate(isAccepted(nit))         ! Logical: Was the proposed model accepted?
! if (verbosity.ge.2) then
!     write(stdout,*) 'anneal_fem_post: finished allocating array memory'
! endif
!
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
!
!
!
! !----
! ! Print requested quantities
! !----
! ! Best fitting model
! if (best_output_file.ne.'') then
!     ! write(*,*) 'anneal_fem_post: working on best-fitting output file'
!     ! Find minimum chi-squared model
!     min_obj = 1e10
!     it = 1
!     do i = 1,nit
!         if (abs(obj(i)).lt.min_obj) then
!             it = i
!             min_obj = abs(obj(i))
!         endif
!     enddo
!
!     ! Write locked/unlocked and slip of best model to file
!     write(stdout,*) 'best fit: ',it,obj(it)
!     open(unit=22,file=best_output_file,status='unknown')
!     do iflt = 1,nflt
!         write(22,*) locked(it,iflt),slip(it,iflt,:)
!     enddo
!     close(22)
!
!     ! Write locked/unlocked and slip of best model to file
!     open(unit=22,file='pole_'//best_output_file,status='unknown')
!     do ipole = 1,npoles
!         write(22,*) poles(it,ipole,:)
!     enddo
!     close(22)
! endif
!
!
! ! Probability locked
! if (plocked_output_file.ne.'') then
!     ! Allocate memory to array
!     allocate(plocked(nflt))
!
!     ! Calculate number of locked instances in search
!     plocked = 0
!     do it = 1,nit
!         if (it.gt.burnin_iterations) then
!             plocked = plocked + locked(it,:)
!         endif
!     enddo
!
!     ! Write locked probability to file
!     open(unit=23,file=plocked_output_file,status='unknown')
!     do iflt = 1,nflt
!         write(23,*) dble(plocked(iflt))/dble(nit-burnin_iterations)
!     enddo
!     close(23)
!     deallocate(plocked)
! endif
!
!
! ! Mean slip (AND STANDARD DEVIATION???)
! if (mean_slip_output_file.ne.'') then
!     ! Compute mean slip
!     allocate(mean_ss(nflt))
!     allocate(mean_ds(nflt))
!     mean_ss = 0.0d0
!     mean_ds = 0.0d0
!     do iflt = 1,nflt
!         do it = 1,nit
!             if (it.gt.burnin_iterations) then
!                 mean_ss(iflt) = mean_ss(iflt) + slip(it,iflt,1)
!                 mean_ds(iflt) = mean_ds(iflt) + slip(it,iflt,2)
!             endif
!         enddo
!         mean_ss(iflt) = mean_ss(iflt)/dble(nit-burnin_iterations)
!         mean_ds(iflt) = mean_ds(iflt)/dble(nit-burnin_iterations)
!     enddo
!
!     ! Write mean slip to file
!     open(unit=24,file=mean_slip_output_file,status='unknown')
!     do iflt = 1,nflt
!         write(24,*) mean_ss(iflt),mean_ds(iflt)
!     enddo
!     close(24)
!
!     deallocate(mean_ss)
!     deallocate(mean_ds)
! endif
!
!
! ! 1-D marginal distributions for slip components
! if (nmarg1d.ge.1) then
!     do i = 1,nmarg1d
!         if (marg1d_flt(i).gt.nflt) then
!             call usage('anneal_fem_post: requested slip in fault greater than number of faults')
!         endif
!         open(unit=25,file=slip_marg_output_file(i),status='unknown')
!         do it = burnin_iterations+1,nit
!             write(25,*) slip(it,marg1d_flt(i),1:2)
!         enddo
!         close(25)
!     enddo
! endif
!
! ! Fraction of interface that is locked
! ! TODO: SCALE LOCKED ZONES BY AREA
! if (nflocked.ge.1) then
!     if (flocked_subset_file(1).eq.'') then
!         open(unit=26,file=flocked_output_file(1),status='unknown')
!         do it = burnin_iterations+1,nit
!             fraction_locked = 0.0d0
!             do iflt = 1,nflt
!                 fraction_locked = fraction_locked + dble(locked(it,iflt))
!             enddo
!             write(26,*) fraction_locked/dble(nflt)
!         enddo
!         close(26)
!     else
!         allocate(isFltInSubset(nflt))
!
!         do i = 1,nflocked
!             nsubset = 0
!             isFLtInSubset = .false.
!             open(unit=27,file=flocked_subset_file(i),status='old')
!             do
!                 read(27,*,iostat=ios) iflt
!                 if (ios.ne.0) then
!                     exit
!                 endif
!                 if (iflt.lt.1) then
!                     call usage('anneal_fem_post: specified locked fault index less than 1')
!                 elseif (iflt.gt.nflt) then
!                     call usage('anneal_fem_post: specified locked fault index greater than nflt')
!                 endif
!                 if (.not.isFltInSubset(iflt)) then
!                     nsubset = nsubset + 1
!                 endif
!                 isFltInSubset(iflt) = .true.
!             enddo
!             close(27)
!
!             open(unit=26,file=flocked_output_file(i),status='unknown')
!             do it = burnin_iterations+1,nit
!                 fraction_locked = 0.0d0
!                 do iflt = 1,nflt
!                     if (isFltInSubset(iflt)) then
!                         fraction_locked = fraction_locked + dble(locked(it,iflt))
!                     endif
!                 enddo
!                 write(26,*) fraction_locked/dble(nsubset)
!             enddo
!             close(26)
!         enddo
!
!         deallocate(isFltInSubset)
!     endif
! endif
!
!
! ! Mean Euler pole (AND STANDARD DEVIATION???)
! if (mean_pole_output_file.ne.'') then
!     ! Compute mean rotation vector
!     allocate(mean_px(npoles))
!     allocate(mean_py(npoles))
!     allocate(mean_pz(npoles))
!     mean_px = 0.0d0
!     mean_py = 0.0d0
!     mean_pz = 0.0d0
!     do ipole = 1,npoles
!         do it = 1,nit
!             call pole_geo2xyz(poles(it,ipole,1),poles(it,ipole,2),poles(it,ipole,3), &
!                               px,py,pz,'sphere')
!             mean_px(ipole) = mean_px(ipole) + px
!             mean_py(ipole) = mean_py(ipole) + py
!             mean_pz(ipole) = mean_pz(ipole) + pz
!         enddo
!         mean_px(ipole) = mean_px(ipole)/dble(nit)
!         mean_py(ipole) = mean_py(ipole)/dble(nit)
!         mean_pz(ipole) = mean_pz(ipole)/dble(nit)
!         call pole_xyz2geo(mean_px(ipole),mean_py(ipole),mean_pz(ipole),plon,plat,pmag,'sphere')
!         mean_px(ipole) = plon
!         mean_py(ipole) = plat
!         mean_pz(ipole) = pmag
!     enddo
!
!     ! Write mean Euler pole to file
!     open(unit=27,file=mean_pole_output_file,status='unknown')
!     do ipole = 1,npoles
!         write(27,*) mean_px(ipole),mean_py(ipole),mean_pz(ipole)
!     enddo
!     close(27)
!
!     deallocate(mean_px)
!     deallocate(mean_py)
!     deallocate(mean_pz)
! endif
!
! ! Marginal distributions for Euler pole location and angular velocity
! if (nmarg_pole.ge.1) then
!     do i = 1,nmarg_pole
!         if (marg_pole(i).gt.npoles) then
!             call usage('anneal_fem_post: requested information in pole greater than number of poles')
!         endif
!         open(unit=28,file=pole_marg_output_file(i),status='unknown')
!         do it = burnin_iterations+1,nit
!             write(28,*) poles(it,marg_pole(i),1:3)
!         enddo
!         close(28)
!     enddo
! endif
!
!
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

! !--------------------------------------------------------------------------------------------------!
!
! subroutine read_models(nit,nflt,npoles,log_format,locked,obj,slip,poles,temp,uncert, &
!                        isAccepted)
!
! use io, only: stdout, stderr, verbosity, progress_indicator
! use trig, only: d2r
!
! use anneal_fem_post, only: lu_log_file
!
! implicit none
!
! ! Arguments
! integer :: nit, nflt, npoles
! character(len=*) :: log_format
! integer :: locked(nit,nflt)
! double precision :: obj(nit), slip(nit,nflt,2), poles(nit,npoles,3), temp(nit), uncert(nit)
! logical :: isAccepted(nit)
!
! ! Local variables
! integer :: i, ios, it, itread, iflt, ipole, nlines_per_it
! character(len=512) :: input_line
! double precision :: rake
!
!
!
! if (verbosity.ge.1) then
!     write(stdout,*) 'read_models: starting'
! endif
!
!
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
! ! Read all other models
! do it = 1,nit
!
!     ! Progress indicator
!     call progress_indicator(it,nit,'anneal_fem_post: read_models',ios)
!
!     ! Read the header for each iteration
!     ! All iteration header lines start with a ">" in the first character and have the format:
!     !     > label=value label=value ...
!     ! Where "label" describes a parameter, and "value" is a number or string for the parameter
!     read(lu_log_file,'(A)',end=1001,iostat=ios) input_line
!     if (it.gt.0) then
!         call parse_iteration_header(input_line,itread,temp(it),obj(it),uncert(it),isAccepted(it),ios)
!     else
!         itread = 0
!     endif
!     if (ios.ne.0) then
!         write(stderr,*) 'anneal_fem_post: error parsing iteration information from line ',trim(input_line)
!         write(stderr,*) 'iteration: ',it
!         call usage('')
!     endif
!     if (itread.ne.it) then
!         call usage('anneal_fem_post: iterations in log file are not in sequential order')
!     endif
!     1001 if (ios.ne.0) then
!         write(stderr,*) 'anneal_fem_post: more iterations indicated (',nit,') than in file (',it-1,')'
!         call usage('')
!     endif
!
!
!     ! Read fault information defined by log file format
!     if (log_format.eq.'anneal') then
!
!         ! Slip magnitudes
!         do iflt = 1,nflt
!             read(lu_log_file,'(A)',end=1002,iostat=ios) input_line
!             if (it.gt.0) then
!                 read(input_line,*,end=1003,err=1003,iostat=ios) slip(it,iflt,1)
!             endif
!         enddo
!
!         ! Rake angles
!         do iflt = 1,nflt
!             read(lu_log_file,'(A)',end=1002,iostat=ios) input_line
!             if (it.gt.0) then
!                 read(input_line,*,end=1003,err=1003,iostat=ios) rake
!                 ! Convert to strike- and dip-slip
!                 slip(it,iflt,2) = slip(it,iflt,1)*sin(rake*d2r)
!                 slip(it,iflt,1) = slip(it,iflt,1)*cos(rake*d2r)
!             endif
!         enddo
!
!
!     elseif (log_format.eq.'anneal_psc') then
!         do iflt = 1,nflt
!             read(lu_log_file,'(A)',end=1002,iostat=ios) input_line
!             if (it.gt.0) then
!                 read(input_line,*,end=1003,err=1003,iostat=ios) locked(it,iflt), slip(it,iflt,1:2)
!             endif
!         enddo
!
!
!     elseif (log_format.eq.'anneal_psc_euler') then
!         do iflt = 1,nflt
!             read(lu_log_file,'(A)',end=1002,iostat=ios) input_line
!             if (it.gt.0) then
!                 read(input_line,*,end=1003,err=1003,iostat=ios) locked(it,iflt), slip(it,iflt,1:2)
!             endif
!         enddo
!         do ipole = 1,npoles
!             do i = 1,3
!                 read(lu_log_file,'(A)',end=1004,iostat=ios) input_line
!                 if (it.gt.0) then
!                     read(input_line,*,end=1005,err=1005,iostat=ios) poles(it,ipole,i)
!                 endif
!             enddo
!         enddo
!
!     else
!         call usage('anneal_fem_post: no log file format named '//trim(log_format))
!     endif
!
! enddo
!
!
! if (verbosity.ge.1) then
!     write(stdout,*) 'read_models: finished'
! endif
!
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
! !--------------------------------------------------------------------------------------------------!
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
! !--------------------------------------------------------------------------------------------------!
!
!
! subroutine parse_iteration_header(input_line,it,temp,obj,uncert,isAccepted,ierr)
! !----
! ! Parse the iteration header, which is in the format:
! ! > Iteration=N Temperature=T Objective=F Model_Uncertainty=U accepted|rejected
! !----
!
! use io, only: stderr, verbosity, stdout
!
! implicit none
!
! ! I/O variables
! character(len=*) :: input_line
! integer :: it, ierr
! double precision :: temp, obj, uncert
! logical :: isAccepted
!
! ! Local variables
! integer :: i, ios
! character(len=64) :: ch
!
!
! ! Initialize values
! ierr = 0
! it = 0
! temp = 0.0d0
! obj = 0.0d0
! uncert = -1.0d0
! isAccepted = .false.
!
!
! if (verbosity.ge.3) then
!     write(stdout,*) 'parse_iteration_header: starting'
! endif
!
! ! Make sure line is properly labeled as an iteration header
! i = index(input_line,'>')
! if (i.ne.1) then
!     write(stderr,*) 'parse_iteration_header: iteration header line does not start with ">"'
!     call usage('Offending line: '//trim(input_line))
! else
!     input_line(1:i) = ''
! endif
!
!
! ! Parse the iteration header
! do
!     ! Read the first entry
!     read(input_line,*,iostat=ios) ch
!     if (ios.ne.0) then
!         exit
!     endif
!
!     if (ch.eq.'accepted') then
!         isAccepted = .true.
!     elseif (ch.eq.'rejected') then
!         isAccepted = .false.
!     else
!         i = index(ch,'=')
!         if (ch(1:i-1).eq.'Iteration') then
!             read(ch(i+1:len(trim(ch))),*) it
!         elseif (ch(1:i-1).eq.'Temperature') then
!             read(ch(i+1:len(trim(ch))),*) temp
!         elseif (ch(1:i-1).eq.'Objective') then
!             read(ch(i+1:len(trim(ch))),*) obj
!         elseif (ch(1:i-1).eq.'Model_Uncertainty') then
!             read(ch(i+1:len(trim(ch))),*) uncert
!         else
!             call usage('parse_iteration_header: could not parse this parameter: '//ch(i:i-1))
!         endif
!     endif
!
!     ! Remove the entry after done parsing
!     i = index(input_line,trim(ch))
!     input_line(1:i+len(trim(ch))) = ''
! enddo
!
! if (verbosity.ge.3) then
!     write(stdout,*) 'parse_iteration_header: finished'
! endif
!
! return
! end subroutine
!
!
!--------------------------------------------------------------------------------------------------!


subroutine set_parameters()

use io, only: stderr, stdout, verbosity

use anneal_fem_post, only: lu_log_file, &
                           nparam, &
                           param_names

implicit none

! ! Arguments
! integer :: nit, nflt, npoles
! character(len=*) ::  log_format

! Local variables
logical :: continueReading
integer :: i, j
integer :: ios
character(len=512) :: input_line


if (verbosity.ge.1) then
    write(stdout,*) 'read_header: starting'
endif


! Initialize reading parameters from first model
continueReading = .true.

! Read the first line
read(lu_log_file,'(A)',end=4001,iostat=ios) input_line

! Initialize number of parameters
nparam = 0

! Count number of parameters
do
    read(lu_log_file,'(A)',end=4001,iostat=ios) input_line
    ! Found the next model, so exit
    if (input_line(1:1).eq.'>') then
        exit
    endif
    ! Found another parameter
    nparam = nparam + 1
enddo


! Allocate parameter name array and read parameter names
allocate(param_names(nparam))
rewind(lu_log_file)
read(lu_log_file,'(A)',end=4001,iostat=ios) input_line
do i = 1,nparam
    read(lu_log_file,'(A)',end=4001,iostat=ios) input_line
    j = index(input_line,'=')
    input_line(j:j) = ' '
    read(input_line,*) param_names(i)
    print *,trim(param_names(i))
enddo


4001 if (ios.ne.0) then
    call usage('anneal_fem_post: reached end of file in header section')
endif

! if (verbosity.ge.2) then
!     write(stdout,*) 'anneal_fem_post: finished parsing header'
!     write(stdout,*) 'niterations=',nit
!     write(stdout,*) 'nfaults=',nflt
!     write(stdout,*) 'npoles=',npoles
!     write(stdout,*) 'log_format= ',log_format
! endif
! if (verbosity.ge.1) then
!     write(stdout,*) 'read_header: starting'
! endif

return
end subroutine


!--------------------------------------------------------------------------------------------------!


subroutine gcmdln()

! use io, only: verbosity
! use random, only: iseed
!
use anneal_fem_post, only: anneal_log_file !, &
!                        burnin_iterations, &
!                        best_output_file, &
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
!                        nresamp

implicit none

! Local variables
integer :: i, narg, ios
character(len=512) :: tag, arg
! logical :: isOutputDefined
!
!
! ! Initialize control parameters
! isOutputDefined = .false.
!
! anneal_log_file = ''
! burnin_iterations = 0
! best_output_file = ''
! plocked_output_file = ''
!
! mean_slip_output_file = ''
! nmarg1d = 0
! marg1d_flt = 0
! slip_marg_output_file = ''
!
! nflocked = 0
! flocked_subset_file = ''
! flocked_output_file = ''
!
! mean_pole_output_file = ''
! nmarg_pole = 0
! marg_pole = 0
! pole_marg_output_file = ''
!
! obj_file = ''
! uncert_file = ''
!
! doResample = .false.
!
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

!     elseif (trim(tag).eq.'-nburn') then
!         i = i + 1
!         call get_command_argument(i,tag)
!         read(tag,*) burnin_iterations
!
!     elseif (trim(tag).eq.'-best') then
!         i = i + 1
!         call get_command_argument(i,best_output_file)
!         isOutputDefined = .true.
!
!     elseif (trim(tag).eq.'-plocked') then
!         i = i + 1
!         call get_command_argument(i,plocked_output_file)
!         isOutputDefined = .true.
!
!     elseif (trim(tag).eq.'-slip:mean') then
!         i = i + 1
!         call get_command_argument(i,mean_slip_output_file)
!         isOutputDefined = .true.
!
!     elseif (trim(tag).eq.'-slip:marg') then
!         nmarg1d = nmarg1d + 1
!         if (nmarg1d.gt.1000) then
!             call usage('anneal_fem_post: requested more than 1000 slip marginals')
!         endif
!         i = i + 1
!         call get_command_argument(i,tag)
!         read(tag,*) marg1d_flt(nmarg1d)
!         i = i + 1
!         call get_command_argument(i,tag)
!         if (len_trim(tag).gt.64) then
!             call usage('anneal_fem_post: slip marginal file must be 64 or less characters')
!         endif
!         slip_marg_output_file(nmarg1d) = tag(1:64)
!         isOutputDefined = .true.
!
!     elseif (trim(tag).eq.'-pole:mean') then
!         i = i + 1
!         call get_command_argument(i,mean_pole_output_file)
!         isOutputDefined = .true.
!
!     elseif (trim(tag).eq.'-pole:marg') then
!         nmarg_pole = nmarg_pole + 1
!         if (nmarg_pole.gt.20) then
!             call usage('anneal_fem_post: requested more than 1000 pole marginals')
!         endif
!         i = i + 1
!         call get_command_argument(i,tag)
!         read(tag,*) marg_pole(nmarg_pole)
!         i = i + 1
!         call get_command_argument(i,tag)
!         if (len_trim(tag).gt.64) then
!             call usage('anneal_fem_post: rotation pole marginal file must be 64 or less characters')
!         endif
!         pole_marg_output_file(nmarg_pole) = tag(1:64)
!         isOutputDefined = .true.
!
!     elseif (trim(tag).eq.'-flocked') then
!         nflocked = 1
!         i = i + 1
!         call get_command_argument(i,flocked_output_file(1))
!         isOutputDefined = .true.
!
!     elseif (trim(tag).eq.'-flocked:subset') then
!         nflocked = nflocked + 1
!         if (nmarg1d.gt.10) then
!             call usage('anneal_fem_post: requested more than 10 locked fraction outputs')
!         endif
!         i = i + 1
!         call get_command_argument(i,flocked_subset_file(nflocked))
!         i = i + 1
!         call get_command_argument(i,flocked_output_file(nflocked))
!         isOutputDefined = .true.
!
!     elseif (tag.eq.'-obj') then
!         i = i + 1
!         call get_command_argument(i,obj_file)
!         isOutputDefined = .true.
!
!     elseif (tag.eq.'-uncertainty') then
!         i = i + 1
!         call get_command_argument(i,uncert_file)
!         isOutputDefined = .true.
!
!     elseif (tag.eq.'-resample') then
!         doResample = .true.
!         nresamp = 0
!         i = i + 1
!         if (i.gt.narg) then
!             i = i - 1
!         else
!             call get_command_argument(i,arg)
!             read(arg,*,iostat=ios) nresamp
!             if (ios.ne.0) then
!                 iseed = 0
!                 i = i - 1
!             endif
!         endif
!
!     elseif (tag.eq.'-seed') then
!         i = i + 1
!         call get_command_argument(i,arg)
!         read(arg,*,iostat=ios) iseed
!         if (ios.ne.0) then
!             iseed = 0
!             i = i - 1
!         endif
!
!     elseif (tag.eq.'-v'.or.tag.eq.'-verbose'.or.tag.eq.'-verbosity') then
!         i = i + 1
!         call get_command_argument(i,arg)
!         read(arg,*) verbosity
!
    else
        call usage('anneal_fem_post: no option '//trim(tag))
    endif

    i = i + 1

enddo
!
! if (.not.isOutputDefined) then
!     call usage('anneal_fem_post: no output file defined')
! endif

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
write(stderr,*) '-param:stat FILE                Parameter statistics'
write(stderr,*) '-param:list                     List available parameters'
write(stderr,*)
write(stderr,*) 'The option "-param:stat" defines statistical quantities to extract from the search'
write(stderr,*) 'It can be repeated as many times as desired'
write(stderr,*) '"param" must be one of the parameters in the search (use -param:list to show them)'
write(stderr,*) '"stat" can be one of:'
write(stderr,*) '    mean'
write(stderr,*) '    median'
write(stderr,*) '    marginal'
write(stderr,*)
! write(stderr,*) '-nit NITER                      Number of iterations in search'
! write(stderr,*) '-nburn NBURN                    Number of burn-in iterations in search'
! write(stderr,*) '-best BEST_FILE                 Best fitting model (locked ss ds): also pole_BEST_FILE'
! write(stderr,*) '-slip:mean SLIP_FILE            Mean slip in each fault'
! write(stderr,*) '-slip:marg IFLT FILE            1-D marginal for IFLT slip (repeat for other faults)'
! write(stderr,*) '-pole:mean POLE_FILE            Mean Euler poles'
! write(stderr,*) '-pole:marg IPOL FILE            Marginals for IPOL rotation pole'
! write(stderr,*) '-flocked FILE                   Fraction of locked faults in each iteration'
! write(stderr,*) '-flocked:subset FLT_FILE FILE   Fraction of locked faults in each iteration from a ',&
!                                                  'subset of faults'
! write(stderr,*) '-obj OBJ_FILE                   Iteration and objective function'
! write(stderr,*) '-uncertainty UNCERT_FILE        Uncertainty distribution'
! write(stderr,*) '-resample NRESAMPLE             Resample models proportional to probability'
! write(stderr,*) '-seed SEED                      Set random number seed'
! write(stderr,*) '-v LVL                          Turn on verbose mode'
! write(stderr,*)

call error_exit(1)
end subroutine
