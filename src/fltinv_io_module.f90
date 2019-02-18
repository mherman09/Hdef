module io_module
    integer, parameter :: stdin=5
    integer, parameter :: stdout=6
    integer, parameter :: stderr=0

    integer :: verbosity  ! 1=basic progress
                          ! 2=detailed progress
                          ! 3=print command line parsing
                          ! 4=parsed inputs
                          ! 5=detailed intermediate calculations

    type program_data
        character(len=256) :: file
        character(len=8) :: flag
        character(len=8) :: array_type
        integer :: nrecords
        integer :: nfields
        integer, allocatable :: intarray(:,:)
        double precision, allocatable :: array(:,:)
    end type program_data

!--------------------------------------------------------------------------------------------------!
contains
!--------------------------------------------------------------------------------------------------!

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

!--------------------------------------------------------------------------------------------------!

    subroutine read_program_data_file(val)
    !----
    ! Read derived data type program_data from file and set variable values
    !----
    implicit none
    type(program_data) :: val
    ! Local variables
    integer :: i, j, ios
    character(len=32) :: fmt_str, fmt_str_int

    ! Initialize subroutine
    if (verbosity.ge.2) then
        write(stderr,'(A)') 'read_program_data_file says: starting to read file "'// &
                            trim(val%file)//'"'
    endif

    ! Check that a file should be read
    if (val%file.eq.'none') then
        if (verbosity.ge.2) then
            write(stderr,'(A)') 'read_program_data_file says: user specifies no file to read'
            write(stderr,'(A)') 'read_program_data_file says: finished'
            write(stderr,*)
        endif
        return
    endif

    ! Check that nfields has been specified
    if (val%nfields.le.0) then
        call print_usage('!! Error: nfields has not been defined')
    endif

    ! Check that the file exists, then count the number of lines
    call check_file_exists(val%file)
    call count_lines(val%nrecords,val%file)

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
        call print_usage('!! Error: no array type named '//trim(val%array_type))
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
            call print_usage('!! Error: read error on file '//trim(val%file))
        endif
    enddo
    close(21)

    ! Print finished message
    if (verbosity.ge.2) then
        write(stderr,'(A)') 'read_program_data_file says: finished reading file "'// &
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

!--------------------------------------------------------------------------------------------------!

    subroutine check_file_exists(file_name)
    !----
    ! Check that a file exists
    !----
    implicit none
    ! I/O variables
    character(len=*) :: file_name
    ! Local variables
    logical :: ex
    inquire(file=file_name,exist=ex)
    if (.not.ex) then
        call print_usage('!! Error: no file found named '//trim(file_name))
    endif
    return
    end subroutine check_file_exists

!--------------------------------------------------------------------------------------------------!

    subroutine count_lines(nlines,file_name)
    !----
    ! Count the number of lines in a file
    !----
    implicit none
    ! I/O variables
    character(len=*) :: file_name
    integer :: nlines
    ! Local variables
    integer :: ios
    open(unit=41,file=file_name,status='old')
    nlines= 0
    do
        read(41,*,iostat=ios)
        if (ios.eq.0) then
            nlines = nlines + 1
        else
            exit
        endif
    enddo
    close(41)
    return
    end subroutine count_lines

end module io_module
