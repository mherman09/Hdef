program main
implicit none

call initialize_fltinv_variables()
call get_command_line()

call read_fltinv_inputs()

call calc_greens_functions()

call run_inversion()
call write_solution()

end program main

!==================================================================================================!

subroutine initialize_fltinv_variables()
use io_module, only: verbosity, initialize_program_data
use variable_module, only: output_file, displacement, disp_components, prestress, stress_weight, &
                           fault, slip_constraint, rake_constraint, &
                           gf_type, gf_disp, gf_stress, gf_los, &
                           inversion_mode, damping_constant, smoothing_constant, smoothing, &
                           coord_type, halfspace, disp_misfit_file, los_misfit_file, &
                           los, los_weight
use lsqr_module, only: lsqr_mode
use anneal_module, only: anneal_init_mode, anneal_log_file, max_iteration, reset_iteration, &
                         temp_start, temp_minimum, cooling_factor
implicit none

! Initialize program behavior variables
output_file = 'stdout'
inversion_mode = 'lsqr'
gf_type = 'okada_rect'
coord_type = 'cartesian'
disp_components = '123'
disp_misfit_file = 'none'
los_misfit_file = 'none'

! Initialize regularization variables
damping_constant = -1.0d0
smoothing_constant = -1.0d0

! Initialize derived data variable values
call initialize_program_data(fault)
call initialize_program_data(displacement)
call initialize_program_data(prestress)
call initialize_program_data(gf_disp)
call initialize_program_data(gf_stress)
call initialize_program_data(gf_los)
call initialize_program_data(slip_constraint)
call initialize_program_data(rake_constraint)
call initialize_program_data(smoothing)
call initialize_program_data(los)

! Same type of derived data variable, but not always used
call initialize_program_data(halfspace)

! Pretty obvious, but initialize program verbosity variable
verbosity = 0

! Initialize least-squares variables
lsqr_mode = 'gels'
stress_weight = 1.0d-9
los_weight = 1.0d0

! Initialize annealing variables
anneal_init_mode = 'mean'
anneal_log_file = 'none'
max_iteration = 1000
reset_iteration = 1000000
temp_start = 2.0d0
temp_minimum = 0.00d0
cooling_factor = 0.98d0

return
end

!--------------------------------------------------------------------------------------------------!

subroutine get_command_line()
use io_module, only: stderr, verbosity
use variable_module, only: output_file, displacement, disp_components, prestress, stress_weight, &
                           fault, slip_constraint, rake_constraint, &
                           gf_type, gf_disp, gf_stress, &
                           inversion_mode, damping_constant, smoothing_constant, smoothing, &
                           coord_type, halfspace, disp_misfit_file, los_misfit_file, &
                           los, los_weight
use lsqr_module, only: lsqr_mode
use anneal_module, only: anneal_init_mode, anneal_log_file, max_iteration, reset_iteration, &
                         temp_start, temp_minimum, cooling_factor
implicit none
! Local variables
integer :: i, narg
character(len=256) :: tag
integer :: char_index

narg = command_argument_count()
if (narg.eq.0) call print_usage('')

i = 1
do while (i.le.narg)

    call get_command_argument(i,tag)

    ! Inversion mode
    if (trim(tag).eq.'-mode') then
        i = i + 1
        call get_command_argument(i,inversion_mode)

    ! Output options
    elseif (trim(tag).eq.'-o') then
        i = i + 1
        call get_command_argument(i,output_file)

    ! Input options
    elseif (trim(tag).eq.'-disp') then
        i = i + 1
        call get_command_argument(i,displacement%file)
    elseif (trim(tag).eq.'-disp:components') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) disp_components
    elseif (trim(tag).eq.'-disp:misfit') then
        i = i + 1
        call get_command_argument(i,disp_misfit_file)
    elseif (trim(tag).eq.'-los:misfit') then
        i = i + 1
        call get_command_argument(i,los_misfit_file)
    elseif (trim(tag).eq.'-prests') then
        i = i + 1
        call get_command_argument(i,prestress%file)
    elseif (trim(tag).eq.'-prests:weight') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) stress_weight
    elseif (trim(tag).eq.'-flt') then
        i = i + 1
        call get_command_argument(i,fault%file)
    elseif (trim(tag).eq.'-flt:rake') then
        i = i + 1
        call get_command_argument(i,rake_constraint%file)
    elseif (trim(tag).eq.'-flt:slip') then
        i = i + 1
        call get_command_argument(i,slip_constraint%file)
    elseif (trim(tag).eq.'-los') then
        i = i + 1
        call get_command_argument(i,los%file)
    elseif (trim(tag).eq.'-los:weight') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) los_weight

    ! Green's functions options
    elseif (trim(tag).eq.'-gf:model') then
        i = i + 1
        call get_command_argument(i,gf_type)
    elseif (trim(tag).eq.'-gf:disp') then
        i = i + 1
        call get_command_argument(i,gf_disp%file)
    elseif (trim(tag).eq.'-gf:stress') then
        i = i + 1
        call get_command_argument(i,gf_stress%file)

    ! Inversion options
    elseif (trim(tag).eq.'-damp') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) damping_constant
    elseif (trim(tag).eq.'-smooth') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) smoothing_constant
        i = i + 1
        call get_command_argument(i,smoothing%file)

    ! Miscellaneous options
    elseif (trim(tag).eq.'-geo') then
        coord_type = 'geographic'
    elseif (trim(tag).eq.'-haf') then
        i = i + 1
        call get_command_argument(i,halfspace%file)
        i = i + 1
        if (i.gt.narg) then
            return
        else
            call get_command_argument(i,tag)
            char_index = index(tag,'-')
            if (char_index.eq.0) then
                call get_command_argument(i,halfspace%flag)
            else
                i = i - 1
            endif
        endif
    elseif (trim(tag).eq.'-v') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) verbosity

    ! Least squares options
    elseif (trim(tag).eq.'-lsqr:mode') then
        i = i + 1
        call get_command_argument(i,lsqr_mode)

    ! Simulated annealing options
    elseif (trim(tag).eq.'-anneal:init_mode') then
        i = i + 1
        call get_command_argument(i,anneal_init_mode)
    elseif (trim(tag).eq.'-anneal:it_max') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) max_iteration
    elseif (trim(tag).eq.'-anneal:it_reset') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) reset_iteration
    elseif (trim(tag).eq.'-anneal:log_file') then
        i = i + 1
        call get_command_argument(i,anneal_log_file)
    elseif (trim(tag).eq.'-anneal:temp_0') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) temp_start
    elseif (trim(tag).eq.'-anneal:temp_min') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) temp_minimum
    elseif (trim(tag).eq.'-anneal:cool') then
        i = i + 1
        call get_command_argument(i,tag)
        read(tag,*) cooling_factor

    ! No option
    else
        call print_usage('!! Error: No option '//trim(tag))
    endif
    i = i + 1
enddo

if (verbosity.ge.1) then
    write(stderr,'("fltinv verbose output turned on")')
endif
if (verbosity.ge.2) then
    write(stderr,'("Parsed command line inputs")')
    write(stderr,'("    inversion_mode:         ",A)') trim(inversion_mode)
    write(stderr,*)
    write(stderr,'("    output_file:            ",A)') trim(output_file)
    write(stderr,*)
    write(stderr,'("    displacement%file:      ",A)') trim(displacement%file)
    write(stderr,'("    disp_components:        ",A)') trim(disp_components)
    write(stderr,'("    disp_misfit_file:       ",A)') trim(disp_misfit_file)
    write(stderr,'("    prestress%file:         ",A)') trim(prestress%file)
    write(stderr,'("    stress_weight:          ",1PE14.6)') stress_weight
    write(stderr,'("    fault%file:             ",A)') trim(fault%file)
    write(stderr,'("    rake_constraint%file:   ",A)') trim(rake_constraint%file)
    write(stderr,'("    slip_constraint%file:   ",A)') trim(slip_constraint%file)
    write(stderr,'("    los%file:               ",A)') trim(los%file)
    write(stderr,'("    los_weight:             ",1PE14.6)') los_weight
    write(stderr,'("    los_misfit_file:        ",A)') trim(los_misfit_file)
    write(stderr,*)
    write(stderr,'("    gf_type:                ",A)') trim(gf_type)
    write(stderr,'("    gf_disp%file:           ",A)') trim(gf_disp%file)
    write(stderr,'("    gf_stress%file:         ",A)') trim(gf_stress%file)
    write(stderr,*)
    write(stderr,'("    damping_constant:       ",1PE14.6)') damping_constant
    write(stderr,'("    smoothing_constant:     ",1PE14.6)') smoothing_constant
    write(stderr,'("    smoothing_file:         ",A)') trim(smoothing%file)
    write(stderr,*)
    write(stderr,'("    coord_type:             ",A)') trim(coord_type)
    write(stderr,'("    halfspace%file:         ",A)') trim(halfspace%file)
    write(stderr,'("    halfspace%flag:         ",A)') trim(halfspace%flag)
    write(stderr,*)
    write(stderr,'("    lsqr_mode:              ",A)') trim(lsqr_mode)
    write(stderr,*)
    write(stderr,'("    anneal_init_mode:       ",A)') trim(anneal_init_mode)
    write(stderr,'("    max_iteration:          ",I14)') max_iteration
    write(stderr,'("    reset_iteration:        ",I14)') reset_iteration
    write(stderr,'("    temp_start:             ",1PE14.6)') temp_start
    write(stderr,'("    temp_minimum:           ",1PE14.6)') temp_minimum
    write(stderr,'("    anneal_log_file:        ",A)') trim(anneal_log_file)
    write(stderr,'("    cooling_factor:         ",1PE14.6)') cooling_factor
endif
if (verbosity.ge.1) then
    write(stderr,*)
endif

return
end subroutine get_command_line

!--------------------------------------------------------------------------------------------------!

subroutine print_usage(string)
!----
! Print program usage statement and exit
!----
use io_module, only: stderr
implicit none
character(len=*) :: string
integer :: string_length
if (string.ne.'') then
    string_length = len(string)
    write(stderr,'(A)') trim(string)
    write(stderr,'(A)')
endif
write(stderr,'(A)') 'Usage: fltinv ...options...'
write(stderr,'(A)')
write(stderr,'(A)') '-mode INVERSION_MODE         Inversion mode'
write(stderr,*)
write(stderr,'(A)') 'Output Options'
write(stderr,'(A)') '-o OUTPUT_FILE               Output file'
write(stderr,*)
write(stderr,'(A)') 'Input Options'
write(stderr,'(A)') '-disp DISP_FILE              Input displacements'
write(stderr,'(A)') '-disp:components COMPNTS     Specify displacement components'
write(stderr,'(A)') '-disp:misfit MISFIT_FILE     Output RMS misfit to displacements'
write(stderr,'(A)') '-prests PRESTS_FILE          Input pre-stresses'
write(stderr,'(A)') '-prests:weight WEIGHT        Stress weighting factor'
write(stderr,'(A)') '-flt FAULT_FILE              Input faults'
write(stderr,'(A)') '-flt:rake RAKE_FILE          Rake angle constraints'
write(stderr,'(A)') '-flt:slip SLIP_FILE          Slip magnitude constraints'
write(stderr,'(A)') '-los LOS_FILE                Input line-of-sight displacements'
write(stderr,'(A)') '-los:weight LOS_WEIGHT       LOS observation weighting factor'
write(stderr,'(A)') '-los:misfit MISFIT_FILE      Output RMS misfit to LOS displacements'
write(stderr,*)
write(stderr,'(A)') 'Greens Functions Options'
write(stderr,'(A)') '-gf:model MODEL              Greens functions calculation model'
write(stderr,'(A)') '-gf:disp_file GF_DSP_FILE    Pre-computed displacement Greens'
write(stderr,'(A)') '-gf:sts_file GF_STS_FILE     Pre-computed stress Greens functions'
write(stderr,*)
write(stderr,'(A)') 'Inversion Options'
write(stderr,'(A)') '-damp DAMP                   Damping regularization'
write(stderr,'(A)') '-smooth SMOOTH SMOOTH_FILE   Smoothing regularization'
write(stderr,*)
write(stderr,'(A)') 'Miscellaneous Options'
write(stderr,'(A)') '-geo                         Treat input coordinates as geographic'
write(stderr,'(A)') '-haf HALFSPACE_FILE [FLAG]   Elastic half-space parameters'
write(stderr,'(A)') '-v LEVEL                     Program verbosity'
write(stderr,*)
write(stderr,'(A)') 'Least-Squares Options'
write(stderr,'(A)') '-lsqr:mode MODE              Solver algorithm'
write(stderr,*)
write(stderr,'(A)') 'Simulated Annealing Options'
write(stderr,'(A)') '-anneal:init_mode MODE       Mode to initialize solution'
write(stderr,'(A)') '-anneal:it_max IMAX          Maximum number of iterations'
write(stderr,'(A)') '-anneal:it_reset IRESET      Reset search every IRESET iterations'
write(stderr,'(A)') '-anneal:log_file LOG_FILE    Log annealing progress'
write(stderr,'(A)') '-anneal:temp_0 START_TEMP    Starting temperature'
write(stderr,'(A)') '-anneal:temp_min MIN_TEMP    Minimum temperature'
write(stderr,'(A)') '-anneal:cool COOL_FACT       Cooling factor'
write(stderr,*)
write(stderr,'(A)') 'See man page for details'
write(stderr,*)
stop
end subroutine print_usage
