!--------------------------------------------------------------------------------------------------!
!------------------ SIMULATED ANNEALING WITH PSEUDO-COUPLING & RIGID ROTATIONS --------------------!
!--------------------------------------------------------------------------------------------------!

subroutine invert_anneal_euler_psc()

use annealing, only: anneal

use fltinv, only: fault, &
                  npoles, &
                  max_iteration, &
                  reset_iteration, &
                  temp_start, &
                  temp_minimum, &
                  cooling_factor, &
                  fault_slip, &
                  euler_pole

implicit none

! Interface to driver subroutines
interface
    subroutine anneal_psc_euler_init(model,n)
        integer :: n
        integer :: model(n)
    end subroutine
    subroutine anneal_psc_euler_propose(model_current,model_proposed,n)
        integer :: n
        integer :: model_current(n)
        integer :: model_proposed(n)
    end subroutine
    function anneal_psc_euler_objective(model,n)
        integer :: n
        integer :: model(n)
        double precision :: anneal_psc_euler_objective
    end function
    subroutine anneal_psc_euler_log(it,temp,obj,model_current,model_proposed,n,isAccepted,string)
        integer :: it, n
        double precision :: temp, obj
        integer :: model_current(n)
        integer :: model_proposed(n)
        logical :: isAccepted
        character(len=*) :: string
    end subroutine
end interface

! Local variables
integer :: ierr, nflt, ntotal, i
integer, allocatable :: locked_pole(:)


nflt = fault%nrows
ntotal = nflt+3*npoles
allocate(locked_pole(ntotal),stat=ierr)
if (ierr.ne.0) then
    call usage('invert_anneal_euler_psc: error allocating memory to locked_pole array (usage:none)')
endif

! Call anneal routine with specific driver routines anneal_psc_euler_init, anneal_psc_euler_propose,
! anneal_psc_euler_objective, and anneal_psc_euler_log
call anneal(ntotal, &
            locked_pole, &
            anneal_psc_euler_init, &
            anneal_psc_euler_propose, &
            anneal_psc_euler_objective, &
            max_iteration, &
            reset_iteration, &
            temp_start, &
            temp_minimum, &
            cooling_factor, &
            anneal_psc_euler_log)

! Save results for printing
fault_slip = 0.0d0

! Calculate slip in each fault
call psc_slip(locked_pole(1:nflt),nflt,fault_slip)

! Rotation
euler_pole = 0.0d0
do i = 1,npoles
    euler_pole(i,1) = dble(locked_pole(nflt+3*i-2))/1.0d3
    euler_pole(i,2) = dble(locked_pole(nflt+3*i-1))/1.0d3
    euler_pole(i,3) = dble(locked_pole(nflt+3*i  ))/1.0d4
enddo


return
end subroutine






!--------------------------------------------------------------------------------------------------!







subroutine anneal_psc_euler_init(model,n)
!----
! Initialize the annealing with pseudo-coupling variables:
!     - iseed: random number seed
!     - model array: 0=unlocked, 1=locked and mapping Euler pole into integer values
!     - Maximum size Greens function matrix for pseudo-coupling inversion
!----

use io, only: stderr, stdout, verbosity, fileExists, line_count
use random, only: iseed, timeseed, r8_uniform_01
use geom, only: distaz2lola, lola2distaz
use earth, only: radius_earth_km

use fltinv, only: fault, &
                  slip_constraint, &
                  npoles, &
                  pole_array, &
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
double precision :: plocked, dist, az, lon, lat, rot


if (verbosity.ge.2) then
    write(stdout,*) 'anneal_psc_euler_init: starting'
endif

! Initialize the random number generator seed
if (anneal_seed.eq.0) then
    iseed = -timeseed()
else
    iseed = -abs(anneal_seed)
endif

! Check array dimensions
nflt = fault%nrows
if (nflt+3*npoles.ne.n) then
    call usage('anneal_psc_euler_init: input n not equal to nflt+3*npoles (usage:none)')
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


! Initialize the model array locking section values
model = 0
if (anneal_init_mode.eq.'locked') then
    ! Set the initial solution to all locked
    model(1:nflt) = 1

elseif (anneal_init_mode.eq.'unlocked') then
    ! Set the initial solution to all unlocked
    model(1:nflt) = 0

elseif (anneal_init_mode(1:4).eq.'rand') then
    ! Set the initial solution to randomly locked, with probability set by the last characters
    ! of the anneal_init_mode variable, e.g., rand0.25 sets probability locked = 0.25
    read(anneal_init_mode(5:len_trim(anneal_init_mode)),*,iostat=ios) plocked
    if (ios.ne.0) then
        write(stderr,*) 'anneal_psc_euler_init: could not read initial plocked, setting to 0.5'
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
        call usage('anneal_psc_euler_init: no initialization file defined for anneal_init_mode='//&
                   'user (usage:anneal)')
    elseif (.not.fileExists(anneal_init_file)) then
        call usage('anneal_psc_euler_init: no anneal_init_file found named "'//&
                   trim(anneal_init_file)//'" (usage:anneal)')
    endif
    if (line_count(anneal_init_file).ne.nflt) then
        call usage('anneal_psc_euler_init: number of lines in anneal_init_file must be equal '//&
                   'to nflt (usage:anneal')
    endif
    open(unit=29,file=anneal_init_file,status='old')
    do i = 1,nflt
        read(29,*,iostat=ios) model(i)
        if (ios.ne.0) then
            call usage('anneal_psc_euler_init: error reading anneal init file (usage:anneal)')
        endif
    enddo
    close(29)

else
    write(stderr,*) 'anneal_psc_euler_init: no initialization mode named "'//&
                    trim(anneal_init_mode)//'"'
    write(stderr,*) 'Options for annealing with pseudo-coupling & rigid rotations initialization:'
    write(stderr,*) '    locked'
    write(stderr,*) '    unlocked'
    write(stderr,*) '    rand'
    call usage(     '    user (usage:anneal)')
endif


! Initialize the model array Euler pole values randomly within prior radius
do i = 1,npoles
    dist = r8_uniform_01(iseed)*pole_array(i,3)/radius_earth_km
    az = r8_uniform_01(iseed)*360.0d0
    call distaz2lola(pole_array(i,1),pole_array(i,2),dist,az,lon,lat,'radians','degrees',ios)
    model(nflt+3*i-2) = int(lon*1.0d3)
    model(nflt+3*i-1) = int(lat*1.0d3)

    ! J/K: initialize with center point
    model(nflt+3*i-2) = int(pole_array(i,1)*1.0d3)
    model(nflt+3*i-1) = int(pole_array(i,2)*1.0d3)

    rot = pole_array(i,4) + r8_uniform_01(iseed)*(pole_array(i,5)-pole_array(i,4))
    model(nflt+3*i-0) = int(rot*1.0d4)
enddo


if (verbosity.ge.2) then
    write(stdout,*) 'anneal_psc_euler_init: finished'
endif

return
end subroutine






!--------------------------------------------------------------------------------------------------!







subroutine anneal_psc_euler_propose(model_in,model_out,n)

use io, only: stdout, stderr, verbosity
use random, only: iseed, r8_uniform_01, r8_normal_ab
use geom, only: distaz2lola, lola2distaz
use earth, only: radius_earth_km

use fltinv, only: fault, &
                  slip_constraint, &
                  npoles, &
                  pole_array, &
                  min_flip, &
                  max_flip

implicit none

! Arguments
integer, intent(in) :: n
integer :: model_in(n), model_out(n)

! Local variables
integer :: i, j, k, ios, nflt, iflip, nflip, rand_fault_list(n)
double precision :: dist, az, lon_current, lat_current, rot_current, lon, lat, rot, drot, &
                    dist_from_center


if (verbosity.ge.3) then
    write(stdout,*) 'anneal_psc_euler_propose: starting'
endif

! Check array dimensions
nflt = fault%nrows
if (nflt+3*npoles.ne.n) then
    call usage('anneal_psc_euler_propose: input n not equal to nflt+3*npoles (usage:none)')
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
        write(stderr,*) 'anneal_psc_euler_propose: invalid locking state ',model_out(rand_fault_list(i))
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


! Perturb the Euler pole location and rotation rate
do i = 1,npoles
    ! Pole coordinates
    dist_from_center = 1.0d10
    do while (dist_from_center.gt.pole_array(i,3))
        lon_current = dble(model_in(nflt+3*i-2))/1.0d3
        lat_current = dble(model_in(nflt+3*i-1))/1.0d3
        dist = r8_normal_ab(0.0d0,pole_array(i,3)/25.0d0/radius_earth_km,iseed)
        az = r8_uniform_01(iseed)*360.0d0
        call distaz2lola(lon_current,lat_current,dist,az,lon,lat,'radians','degrees',ios)
        call lola2distaz(lon,lat,pole_array(i,1),pole_array(i,2),dist_from_center,az,'radians', &
                         'degrees',ios)
        dist_from_center = dist_from_center*radius_earth_km
    enddo

    ! Rotation rate
    rot_current = dble(model_in(nflt+3*i-0))/1.0d4
    rot = 1.0d10
    do while (rot.gt.pole_array(i,5).or.rot.lt.pole_array(i,4))
        drot = (pole_array(i,5)-pole_array(i,4))/25.0d0
        rot = rot_current + drot
        rot = r8_normal_ab(rot_current,drot,iseed)
    enddo

    ! Save proposed Euler pole
    model_out(nflt+3*i-2) = int(lon*1.0d3)
    model_out(nflt+3*i-1) = int(lat*1.0d3)
    model_out(nflt+3*i-0) = int(rot*1.0d4)
enddo

if (verbosity.ge.3) then
    write(stdout,*) 'anneal_psc_euler_propose: finished'
endif

return
end subroutine






!--------------------------------------------------------------------------------------------------!







function anneal_psc_euler_objective(model,n)

use io, only: stdout, verbosity
use earth, only: pole_geo2xyz
use misfit, only: misfit_chi2

use fltinv, only: fault, &
                  npoles, &
                  rigid_pt_array_disp, &
                  displacement, &
                  disp_components, &
                  input_disp_unit, &
                  los, &
                  cov_matrix, &
                  isCovMatrixDiagonal, &
                  gf_disp, &
                  gf_los, &
                  gf_euler
implicit none

! Arguments
integer, intent(in) :: n
integer :: model(n)
double precision ::  anneal_psc_euler_objective

! Local variables
integer :: i, ierr, iflt, idsp, icmp, ipre, nflt, ndsp, ndsp_dof, nlos, nobs, ipole
double precision :: plon, plat, prot, px, py, pz, factor
double precision, allocatable :: pre(:), obs(:), slip(:,:)
double precision, parameter :: spy = 60.0d0*60.0d0*24.0d0*365.25d0


if (verbosity.ge.3) then
    write(stdout,*) 'anneal_psc_euler_objective: starting'
endif

nflt = fault%nrows
if (nflt+3*npoles.ne.n) then
    call usage('anneal_psc_euler_objective: input n not equal to nflt+3*npoles (usage:none)')
endif


! UPDATE slip_constraint%array TO REFLECT THE RELATIVE MOTION OF THE RIGID BODY

! Calculate slip in each fault
if (verbosity.ge.3) then
    write(stdout,*) 'anneal_psc_euler_objective: starting psc_slip'
endif
if (.not.allocated(slip)) then
    allocate(slip(nflt,2))
endif
call psc_slip(model(1:nflt),nflt,slip)
if (verbosity.ge.3) then
    write(stdout,*) 'anneal_psc_euler_objective: finished psc_slip'
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
        call usage('anneal_psc_euler_objective: error allocating memory to obs (usage:none)')
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
        call usage('anneal_psc_euler_objective: error allocating memory to pre (usage:none)')
    endif
endif
pre = 0.0d0

! Set unit factor adjustment for rigid rotation velocities
factor = 0.0d0
if (input_disp_unit.eq.'m/s') then
    factor = 1.0d0
elseif (input_disp_unit.eq.'m/yr') then
    factor = 1.0d0*spy
elseif (input_disp_unit.eq.'mm/s') then
    factor = 1.0d3
elseif (input_disp_unit.eq.'mm/yr') then
    factor = 1.0d3*spy
else
    call usage('anneal_psc_euler_objective: unit '//trim(input_disp_unit)//' not compatible '//&
               '(usage:input)')
endif

! Load predicted three-component displacements/velocities
do i = 1,len_trim(disp_components)
    read(disp_components(i:i),*) icmp
    do idsp = 1,ndsp

        ipre = (i-1)*ndsp+idsp

        ! Fault-generated displacements/velocities
        do iflt = 1,nflt
            pre(ipre) = pre(ipre) + gf_disp%array((icmp-1)*ndsp+idsp,     iflt)*slip(iflt,1)
            pre(ipre) = pre(ipre) + gf_disp%array((icmp-1)*ndsp+idsp,nflt+iflt)*slip(iflt,2)
        enddo

        ! Rigid body rotations
        ipole = rigid_pt_array_disp(idsp)
        if (ipole.ne.0) then
            plon = dble(model(nflt+3*ipole-2))/1.0d3
            plat = dble(model(nflt+3*ipole-1))/1.0d3
            prot = dble(model(nflt+3*ipole-0))/1.0d4
            call pole_geo2xyz(plon,plat,prot,px,py,pz,'sphere')

            ! Green's functions are in m/s, convert to input velocity
            ! Only add to horizontal velocity (vertical is zero)
            if (icmp.ne.3) then
                pre(ipre) = pre(ipre) + gf_euler((icmp-1)*ndsp+idsp,         ipole)*px*factor
                pre(ipre) = pre(ipre) + gf_euler((icmp-1)*ndsp+idsp,npoles+  ipole)*py*factor
                pre(ipre) = pre(ipre) + gf_euler((icmp-1)*ndsp+idsp,npoles+2*ipole)*pz*factor
            endif
        endif
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
    write(stdout,*) 'anneal_psc_euler_objective: starting misfit_chi2'
endif
call misfit_chi2(obs,pre,cov_matrix,isCovMatrixDiagonal,nobs,anneal_psc_euler_objective)
anneal_psc_euler_objective = -0.5d0*anneal_psc_euler_objective
if (verbosity.ge.3) then
    write(stdout,*) 'anneal_psc_euler_objective: finished misfit_chi2'
endif

if (verbosity.ge.3) then
    write(stdout,*) 'anneal_psc_euler_objective: objective=',anneal_psc_euler_objective
    write(stdout,*) 'anneal_psc_euler_objective: finished'
endif

return
end function






!--------------------------------------------------------------------------------------------------!






subroutine anneal_psc_euler_log(it,temp,obj,model_current,model_proposed,n,isModelAccepted,string)

use io, only: verbosity, stdout
use fltinv, only: anneal_log_file, &
                  max_iteration, &
                  fault

implicit none

! Arguments
integer :: it, n, model_current(n), model_proposed(n)
double precision :: temp, obj
logical :: isModelAccepted
character(len=*) :: string

! Local variables
integer :: i, nflt
double precision, allocatable :: slip(:,:)
character(len=8) :: rejected_string
character(len=16) :: str, it_str, temp_str, obj_str, uncert_str


if (anneal_log_file.eq.'') then
    return
endif

nflt = fault%nrows
allocate(slip(nflt,2))


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
        write(stdout,*) 'anneal_psc_euler_log: initializing log file'
    endif

    ! Open the log file
    open(unit=27,file=anneal_log_file,status='unknown')

    ! Compute slip for sub-faults
    call psc_slip(model_current(1:nflt),nflt,slip)
    ! do i = 1,nflt
    !     write(0,*) 'anneal_psc_euler_log',i,model_current(i),slip(i,:)
    ! enddo

    ! Write header
    write(str,'(I8)') max_iteration
    write(27,'(A)') '# niterations='//adjustl(str)
    write(str,'(I8)') nflt
    write(27,'(A)') '# nfaults='//adjustl(str)
    write(str,'(I8)') (n-nflt)/3
    write(27,'(A)') '# npoles='//adjustl(str)
    write(27,'(A)') '# format=anneal_psc_euler'

    ! Write iteration parameters, fault slip values
    write(27,2701) it_str, temp_str, obj_str, uncert_str
    2701 format('> Iteration=',A8,X,'Temperature=',A12,X,'Objective=',A12,X,'Model_Uncertainty=',&
                A12)

    ! Write locked/unlocked, fault slip results to log file
    do i = 1,nflt
        write(27,2702) model_current(i),slip(i,1),slip(i,2)
    enddo
    2702 format(I4,X,1PE14.6,E14.6)

    do i = nflt+1,n
        if (mod(i-nflt,3).eq.1.or.mod(i-nflt,3).eq.2) then
            write(27,2703) dble(model_current(i))/1.0d3
        else
            write(27,2704) dble(model_current(i))/1.0d4
        endif
    enddo
    2703 format(F12.3)
    2704 format(1PE14.6)

elseif (string.eq.'append') then
    ! Is this a rejected model?
    if (isModelAccepted) then
        rejected_string = 'accepted'
    else
        rejected_string = 'rejected'
    endif

    ! Compute slip for sub-faults
    call psc_slip(model_current(1:nflt),nflt,slip)

    ! Write locked/unlocked, fault slip, old model results to log file
    write(27,2705) it_str, temp_str, obj_str, uncert_str, trim(rejected_string)
    2705 format('> Iteration=',A8,X,'Temperature=',A12,X,'Objective=',A12,X,'Model_Uncertainty=',&
                A12,X,A)

    do i = 1,nflt
        write(27,2706) model_current(i),slip(i,1),slip(i,2),model_proposed(i)
    enddo
    2706 format(I4,X,1PE14.6,E14.6,I4)

    do i = nflt+1,n
        if (mod(i-nflt,3).eq.1.or.mod(i-nflt,3).eq.2) then
            write(27,2707) dble(model_current(i))/1.0d3,dble(model_proposed(i))/1.0d3
        else
            write(27,2708) dble(model_current(i))/1.0d4,dble(model_proposed(i))/1.0d4
        endif
    enddo
    2707 format(F12.3,F12.3)
    2708 format(1PE14.6,E14.6)

elseif (string.eq.'close') then
    ! Close the log file
    close(27)

else
    call usage('anneal_psc_euler_log: no string option named '//trim(string)//' (usage:none)')
endif

deallocate(slip)

return
end subroutine
