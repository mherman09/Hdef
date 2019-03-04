module variable_module

use io, only: stderr, verbosity

    type program_data
        character(len=256) :: file
        character(len=8) :: flag
        character(len=8) :: array_type
        integer :: nrecords
        integer :: nfields
        integer, allocatable :: intarray(:,:)
        double precision, allocatable :: array(:,:)
    end type program_data

    ! Command line variables
    character(len=256) :: output_file             ! stdout, <file_name>
    character(len=16) :: inversion_mode           ! lsqr, anneal, neighbor
    character(len=16) :: gf_type                  ! okada_rect, okada_pt
    character(len=16) :: coord_type               ! cartesian, geographic
    character(len=3) :: disp_components           ! 123, 1, 2, 3, 12, 13, 23
    character(len=8) :: lsqr_mode

    double precision :: damping_constant
    double precision :: smoothing_constant

    ! Derived data types containing file name, array size, and array information
    type(program_data) :: fault
    type(program_data) :: displacement
    type(program_data) :: prestress
    type(program_data) :: gf_disp
    type(program_data) :: gf_stress
    type(program_data) :: gf_los
    type(program_data) :: slip_constraint
    type(program_data) :: rake_constraint
    type(program_data) :: smoothing
    type(program_data) :: los

    ! Options depending on gf_type
    type(program_data) :: halfspace

    ! Other variables
    integer, allocatable :: smoothing_neighbors(:)
    double precision :: stress_weight
    double precision :: sts_dist
    double precision, allocatable :: fault_slip(:,:)
    character(len=256) :: disp_misfit_file
    character(len=256) :: disp_cov_file
    double precision, allocatable :: disp_cov_mat(:,:)
    character(len=256) :: los_misfit_file
    double precision :: los_weight

    ! Annealing control options
    character(len=8) :: anneal_init_mode          ! zero, mean, uniform, random, user
    character(len=256) :: anneal_log_file         ! file logging annealing progress
    integer :: max_iteration                      ! maximum number of steps
    integer :: reset_iteration                    ! step to reset current solution to best, temp to temp_start
    double precision :: temp_start                ! >0: mult by initial obj value; <0: use absolute value
    double precision :: temp_minimum              ! >0: mult by initial obj value; <0: use absolute value
    double precision :: cooling_factor            ! every iteration, temp->temp*cooling_factor
    integer :: anneal_verbosity
    character(len=256) :: anneal_control_file     ! define annealing parameters in file instead of cmdln
    character(len=256) :: anneal_init_file
    double precision :: prob_lock2unlock
    double precision :: prob_unlock2lock

contains

    subroutine initialize_program_data(val)
    !----
    ! Reset derived data type program_data variable values
    !----
    implicit none
    type(program_data) :: val
    val%file = 'none'
    val%flag = ''
    val%array_type = 'dp'
    val%nrecords = 0
    val%nfields = 0
    if (allocated(val%array)) then
        deallocate(val%array)
    endif
    if (allocated(val%intarray)) then
        deallocate(val%intarray)
    endif
    return
    end subroutine initialize_program_data



    subroutine read_program_data_file(val)
    !----
    ! Read derived data type program_data from file and set variable values
    !----
    use io, only: fileExists, line_count
    implicit none
    type(program_data) :: val
    ! Local variables
    integer :: i, j, ios
    character(len=32) :: fmt_str, fmt_str_int

    ! Initialize subroutine
    if (verbosity.ge.2) then
        write(stderr,'(A)') 'read_program_data_file: starting to read file "'// &
                            trim(val%file)//'"'
    endif

    ! Check that a file should be read
    if (val%file.eq.'none') then
        if (verbosity.ge.2) then
            write(stderr,'(A)') 'read_program_data_file: user specifies no file to read'
            write(stderr,'(A)') 'read_program_data_file: finished'
            write(stderr,*)
        endif
        return
    endif

    ! Check that nfields has been specified
    if (val%nfields.le.0) then
        write(stderr,*) 'read_program_data_file: nfields has not been defined'
        stop
    endif

    ! Check that the file exists, then count the number of lines
    ! call check_file_exists(val%file)
    ! call count_lines(val%nrecords,val%file)
    if (.not.fileExists(val%file)) then
        write(0,*) 'read_program_data_file: no file found named "',trim(val%file),'"'
        stop
    endif
    val%nrecords = line_count(val%file)

    ! Allocate memory for the associated array
    if (val%array_type.eq.'dp') then
        if (.not.allocated(val%array)) then
            allocate(val%array(val%nrecords,val%nfields))
        endif
    elseif (val%array_type.eq.'int') then
        if (.not.allocated(val%intarray)) then
            allocate(val%intarray(val%nrecords,val%nfields))
        endif
    else
        write(stderr,*) 'read_program_data_file: no array type named '//trim(val%array_type)
    endif

    ! Read the file, in free format
    open(unit=21,file=val%file,status='old')
    do i = 1,val%nrecords
        if (val%array_type.eq.'dp') then
            read(21,*,iostat=ios) (val%array(i,j),j=1,val%nfields)
        elseif (val%array_type.eq.'int') then
            read(21,*,iostat=ios) (val%intarray(i,j),j=1,val%nfields)
        endif
        if (ios.ne.0) then
            write(stderr,*) 'read_program_data_file: read error on file '//trim(val%file)
        endif
    enddo
    close(21)

    ! Print finished message
    if (verbosity.ge.2) then
        write(stderr,'(A)') 'read_program_data_file: finished reading file "'// &
                            trim(val%file)//'"'
    endif
    if (verbosity.ge.3) then
        write(stderr,'("nrecords: ",I5)') val%nrecords
        write(fmt_str,1001) val%nfields
        1001 format('(1P',I5,'E14.6)')
        write(fmt_str_int,1002) val%nfields
        1002 format('(',I5,'I14)')
        do i = 1,val%nrecords
            if (val%array_type.eq.'dp') then
                write(stderr,fmt=fmt_str) (val%array(i,j),j=1,val%nfields)
            elseif (val%array_type.eq.'int') then
                write(stderr,fmt=fmt_str_int) (val%intarray(i,j),j=1,val%nfields)
            endif
        enddo
    endif
    if (verbosity.ge.2) then
        write(stderr,*)
    endif

    return
    end subroutine read_program_data_file

end module variable_module
