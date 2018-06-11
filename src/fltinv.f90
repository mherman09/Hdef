program main
!----
! Invert displacement observations for fault slip, assuming slip in an elastic half-space.
! Currently computes Green's functions from equations of Okada (1992). fltinv can also
! minimize shear stresses resolved onto faults.
!
! Inversion options:
!     - Linear least squares
!     - Simulated annealing
!     - Neighborhood algorithm
!----
use command_line
use arrays
use io
implicit none

!----
! Parse command line and check inputs
!----
call gcmdln()
call check_inputs()

!----
! Read control files
!----
if (verbosity.ge.1) then
    write(stderr,'(A)') 'Reading inversion parameters from files'
    write(stderr,*)
endif
! Read fault geometries: faults(nfaults,7)
call read_faults()
! Read displacement observations: displacements(ndisplacements,6)
call read_displacements()
! Read pre-stresses: prestresses(nfaults,6)
call read_prestresses()
! Read half-space elastic parameters
call read_halfspace()
! Read rake angle(s): rakes(nfaults)
call read_rakes()
! Read smoothing file: smoothing_pointers(nsmooth,3), smoothing_array(nneighbors_total)
call read_smoothing()
! Read slip constraint file: slip_constraints(nfaults,2), is_this_fault_constrained(nfaults,2)
call read_slip_constraints()

!----
! Post-read checks
!----
! If coordinates are specified as geographic, convert to Cartesian
if (geographic.eq.1) then
    call geo2xy(nfaults,faults)
    call geo2xy(ndisplacements,displacements)
endif

call check_problem_parameters()
if (verbosity.ge.1) then
    write(stderr,'("***** Congrats: everything is read in and appears to be in order *****")')
    write(stderr,*)
endif

!----
! Compute Green's functions
!----
! Displacement Green's functions
call calc_disp_gfs()
! Shear stress Green's functions
call calc_stress_gfs()

!----
! Run inversion
!----
if (inversion_mode.eq.'linear') then
    if (verbosity.ge.1) then
        write(stderr,'(A)') 'Performing linear least squares inversion'
    endif
    call lsqr_invert()
else
    call usage('!! Error: inversion mode '//trim(inversion_mode)//' is not available')
endif

!----
! Print misfit to standard output if requested
!----
!      if (misfit.eq.1) then
!          if (vrb.ge.1) then
!              write(0,'("COMPUTING MISFIT, PER YOUR REQUEST")')
!          endif
!          tot = 0.0d0
!          do 905 i = 1,nobs
!              pre(i,1) = 0.0d0
!              pre(i,2) = 0.0d0
!              pre(i,3) = 0.0d0
!              ! predicted displacement vector at each displacement site
!              do 906 j = 1,nflt
!                  pre(i,1) = pre(i,1) + gf(i,j,1)*soln(j,1)
!     1                                             + gf(i,j,4)*soln(j,2)
!                  pre(i,2) = pre(i,2) + gf(i,j,2)*soln(j,1)
!     1                                             + gf(i,j,5)*soln(j,2)
!                  pre(i,3) = pre(i,3) + gf(i,j,3)*soln(j,1)
!     1                                             + gf(i,j,6)*soln(j,2)
!  906         continue
!  905     continue
!          ! total misfit, squared
!          do 907 i = 1,nobs
!              tot = tot + (pre(i,1)-obs(i,4))*(pre(i,1)-obs(i,4))
!     1                   + (pre(i,2)-obs(i,5))*(pre(i,2)-obs(i,5))
!     2                   + (pre(i,3)-obs(i,6))*(pre(i,3)-obs(i,6))
!  907     continue
!          ! total model length, squared
!          pre(1,1) = 0.0d0
!          do 908 i = 1,nflt
!              pre(1,1) = pre(1,1) + soln(i,1)*soln(i,1)
!     1                   + soln(i,2)*soln(i,2)
!  908     continue
!          write(*,1010) dsqrt(tot),dsqrt(pre(1,1))
!      endif
! 1010 format(1PE14.6,1PE14.6)



!----
! Write output nicely
!----
call write_output()

end

!--------------------------------------------------------------------------------------------------!

subroutine gcmdln()
!----
! Parse the command line
!----
use command_line
use io
implicit none
! Local variables
integer :: i, narg
character(len=256) :: tag
integer :: char_index

output_file = 'stdout'
fault_file = 'none'
displacement_file = 'none'
prestress_file = 'none'
halfspace_file = 'none'
smoothing_file = 'none'
inversion_mode = 'linear'
slip_constraint_file = 'none'
damping_constant = -1.0d0
smoothing_constant = -1.0d0
geographic = 0
disp_comp = '123'
rake_file = 'none'
verbosity = 0

narg = command_argument_count()
if (narg.eq.0) call usage('')
i = 1
do while (i.le.narg)
    call get_command_argument(i,tag)
    if (trim(tag).eq.'-f'.or.trim(tag).eq.'--faults') then
        i = i + 1
        call get_command_argument(i,fault_file)
    elseif (trim(tag).eq.'-o'.or.trim(tag).eq.'--output') then
        i = i + 1
        call get_command_argument(i,output_file)
    elseif (trim(tag).eq.'-d'.or.trim(tag).eq.'--displacements') then
        i = i + 1
        call get_command_argument(i,displacement_file)
    elseif (trim(tag).eq.'-s'.or.trim(tag).eq.'--prestresses') then
        i = i + 1
        call get_command_argument(i,prestress_file)
    elseif (trim(tag).eq.'-h'.or.trim(tag).eq.'--halfspace') then
        i = i + 1
        call get_command_argument(i,halfspace_file)
    elseif (trim(tag).eq.'-m'.or.trim(tag).eq.'--mode') then
        i = i + 1
        call get_command_argument(i,inversion_mode)
    elseif (trim(tag).eq.'--damping') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) damping_constant
    elseif (trim(tag).eq.'--smoothing') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) smoothing_constant
        i = i + 1
        call get_command_argument(i,smoothing_file)
    elseif (trim(tag).eq.'--geographic') then
        geographic = 1
    elseif (trim(tag).eq.'--disp_comp') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) disp_comp
    elseif (trim(tag).eq.'--rake') then
        i = i + 1
        call get_command_argument(i,rake_file)
    elseif (trim(tag).eq.'--slip_constraint') then
        i = i + 1
        call get_command_argument(i,slip_constraint_file)
    elseif (trim(tag).eq.'-v'.or.trim(tag).eq.'--verbose') then
        verbosity = 1
        i = i + 1
        if (i.gt.narg) then
            return
        else
            call get_command_argument(i,tag)
            char_index = index(tag,'-')
            if (char_index.eq.0) then
                read(tag,'(BN,I10)') verbosity
            else
                i = i - 1
            endif
        endif
    else
        call usage('!! Error: No option '//trim(tag))
    endif
    i = i + 1
enddo

if (verbosity.ge.1) then
    write(stderr,'("fltinv verbose output turned on")')
endif
if (verbosity.ge.2) then
    write(stderr,'("Parsed command line inputs")')
    write(stderr,'("    output_file:        ",A)') trim(output_file)
    write(stderr,'("    fault_file:         ",A)') trim(fault_file)
    write(stderr,'("    displacement_file:  ",A)') trim(displacement_file)
    write(stderr,'("    prestress_file:     ",A)') trim(prestress_file)
    write(stderr,'("    halfspace_file:     ",A)') trim(halfspace_file)
    write(stderr,'("    smoothing_file:     ",A)') trim(smoothing_file)
    write(stderr,'("    inversion_mode:     ",A)') trim(inversion_mode)
    write(stderr,'("    damping constant:   ",F20.8)') damping_constant
    write(stderr,'("    smoothing constant: ",F20.8)') smoothing_constant
    write(stderr,'("    geographic flag:    ",I20)') geographic
    write(stderr,'("    disp component:     ",A)') disp_comp
    write(stderr,'("    slip constraint:    ",A)') slip_constraint_file
    write(stderr,'("    rake file:          ",A)') trim(rake_file)
    write(stderr,*)
endif
if (verbosity.ge.1) then
    write(stderr,*)
endif

return
end

!--------------------------------------------------------------------------------------------------!

subroutine usage(string)
!----
! Print program usage statement and exit
!----
use io
implicit none
character(len=*) :: string
integer :: string_length
if (string.ne.'') then
    string_length = len(string)
    write(stderr,*) trim(string)
    write(stderr,*)
endif
write(stderr,*) 'Usage: fltinv ...options...'
write(stderr,*)
write(stderr,*) '-o|--output OUTPUT_FILE'
write(stderr,*) '-m|--mode MODE'
write(stderr,*) '-f|--faults FAULT_FILE'
write(stderr,*) '-d|--displacements DISP_FILE'
write(stderr,*) '-s|--prestresses PRESTS_FILE'
write(stderr,*) '-h|--halfspace HAFSPC_FILE'
write(stderr,*) '--damping DAMP'
write(stderr,*) '--smoothing SMOOTH SMOOTH_FILE'
write(stderr,*) '--geographic'
write(stderr,*) '--disp_comp [1][2][3]'
write(stderr,*) '--rake RAKE'
write(stderr,*) '--slip_constraint CONST_FILE'
write(stderr,*) '-v|--verbose VRB_LEVEL'
stop
end
