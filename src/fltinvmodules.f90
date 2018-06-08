module command_line
    character(len=256) :: fault_file             ! fault locations and geometries
    character(len=256) :: displacement_file      ! observed displacements
    character(len=256) :: prestress_file         ! pre-stresses on faults
    character(len=256) :: halfspace_file         ! elastic half-space parameters
    character(len=256) :: smoothing_file         ! fault neighbor file for smoothing
    character(len=256) :: inversion_mode         ! inversion algorithm
    character(len=1) :: disp_comp                ! component(s) of displacement data to use
    character(len=256) :: rake_file              ! rake(s) of fault slip
    double precision :: damping_constant         ! reduce length of output solution
    double precision :: smoothing_constant       ! reduce curvature of solution
    integer :: geographic                        ! input coordinates are geographic
    integer :: verbosity                         ! print program progress
endmodule command_line

!--------------------------------------------------------------------------------------------------!

module arrays
    integer :: nfaults
    integer :: ndisplacements
    integer :: nignore
    integer :: nsmooth
    integer, dimension(:), allocatable :: smoothing_neighbors
    integer, dimension(:,:), allocatable :: smoothing_pointers
    double precision, dimension(:), allocatable :: rakes
    double precision, dimension(:), allocatable :: fault_slip
    double precision, dimension(:,:), allocatable :: faults
    double precision, dimension(:,:), allocatable :: displacements
    double precision, dimension(:,:), allocatable :: prestresses
    double precision, dimension(:,:,:), allocatable :: disp_gfs    ! ssx ssy ssz dsx dsy dsz
    double precision, dimension(:,:,:), allocatable :: stress_gfs  ! ss->ss ss->ds ds->ss ds->ss
endmodule arrays

!--------------------------------------------------------------------------------------------------!

module io
    integer, parameter :: stderr = 0
    integer, parameter :: stdin = 5
    integer, parameter :: stdout = 6
endmodule io

!--------------------------------------------------------------------------------------------------!

module trig
    double precision, parameter :: pi = datan(1.0d0)*4.0d0
    double precision, parameter :: r2d = 180.0d0/pi
    double precision, parameter :: d2r = pi/180.0d0
endmodule trig

!--------------------------------------------------------------------------------------------------!

module gf
    double precision :: gf_stlo
    double precision :: gf_stla
    double precision :: gf_stdp
    double precision :: gf_evlo
    double precision :: gf_evla
    double precision :: gf_evdp
    double precision :: gf_str
    double precision :: gf_dip
    double precision :: gf_rak
    double precision :: gf_slip
    double precision :: gf_wid
    double precision :: gf_len
    double precision :: gf_vp
    double precision :: gf_vs
    double precision :: gf_dens
endmodule gf
