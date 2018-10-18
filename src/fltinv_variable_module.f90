module variable_module
use io_module, only: program_data
    ! Command line variables
    character(len=256) :: output_file             ! stdout, <file_name>
    character(len=16) :: inversion_mode           ! lsqr, anneal, neighbor
    character(len=16) :: gf_type                  ! okada_rect, okada_pt
    character(len=16) :: coord_type               ! cartesian, geographic
    character(len=3) :: disp_components           ! 123, 1, 2, 3, 12, 13, 23

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
    double precision, allocatable :: fault_slip(:,:)
    character(len=256) :: disp_misfit_file
    double precision :: los_weight
end module variable_module
