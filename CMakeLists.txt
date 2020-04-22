cmake_minimum_required(VERSION 3.12)
# cmake -DCMAKE_INSTALL_PREFIX=`pwd` -DCMAKE_C_COMPILER=gcc-8 -DCMAKE_Fortran_COMPILER=gfortran -DLAPACK_LIB_DIR=/sw/lib/lapack ..

project(Hdef VERSION 2020.04)

# This project uses mostly Fortran and a bit of C
enable_language(C Fortran)

# Turn off OSX sysroot flags
if (NOT "${CMAKE_OSX_SYSROOT}" STREQUAL "")
    unset(CMAKE_OSX_SYSROOT CACHE)
endif()
if (NOT "${CMAKE_OSX_DEPLOYMENT_TARGET}" STREQUAL "")
    unset(CMAKE_OSX_DEPLOYMENT_TARGET CACHE)
endif()


message(STATUS "")
message(STATUS "*********************************************************************************")
message(STATUS "************************* CONFIGURING HDEF INSTALLATION *************************")
message(STATUS "*********************************************************************************")
message(STATUS "")


#####
#	INSTALLATION DIRECTORY
#####
message(STATUS "*********************************************************************************")
message(STATUS "Configuring installation directory")
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    unset(CMAKE_INSTALL_PREFIX CACHE)
    message(FATAL_ERROR "Set installation directory with: '-DCMAKE_INSTALL_PREFIX=/path/to/hdef/install'")
endif()
message(STATUS "CMAKE_INSTALL_PREFIX=" ${CMAKE_INSTALL_PREFIX})
message(STATUS "*********************************************************************************")
message(STATUS "")


#####
#	COMPILER SETUP
#####
# Hdef build must use gcc and gfortran (for now...)
message(STATUS "*********************************************************************************")
message(STATUS "Configuring compilers")

# Check for gcc
message(STATUS "CMAKE_C_COMPILER=" ${CMAKE_C_COMPILER})
if (NOT CMAKE_C_COMPILER_ID MATCHES "GNU")
    unset(CMAKE_C_COMPILER CACHE)
    unset(CMAKE_Fortran_COMPILER CACHE)
    message(STATUS "Hdef is only set up to compile with the GNU compiler set")
    message(FATAL_ERROR "Set C compiler with: '-DCMAKE_C_COMPILER=/path/to/gcc'")
endif()

# Check for gfortran
message(STATUS "CMAKE_Fortran_COMPILER=" ${CMAKE_Fortran_COMPILER})
if (NOT CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    unset(CMAKE_C_COMPILER CACHE)
    unset(CMAKE_Fortran_COMPILER CACHE)
    message(STATUS "Hdef is only set up to compile with the GNU compiler set")
    message(FATAL_ERROR "Set Fortran compiler with: '-DCMAKE_Fortran_COMPILER=/path/to/gfortran'")
endif()

# Use pre-processing commands
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp")

# Compilation flags
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -Wextra -Wunused -Wconversion -fcheck=all -fbacktrace")

message(STATUS "*********************************************************************************")
message(STATUS "")




#####
#	LAPACK
#####
message(STATUS "*********************************************************************************")
message(STATUS "Looking for LAPACK")

set(useLAPACK 0)

if ("${LAPACK_LIB_DIR}" STREQUAL "" )
    message(STATUS "Compiling Hdef without LAPACK")
    message(STATUS "Set LAPACK directory with: '-DLAPACK_LIB_DIR=/path/to/lapack/lib'")

elseif (EXISTS ${LAPACK_LIB_DIR})
    message(STATUS "Found LAPACK directory: " ${LAPACK_LIB_DIR})
    set(useLAPACK 1)
    unset(LAPACK_LIBRARIES)

    # Find lapack lib
    if (EXISTS ${LAPACK_LIB_DIR}/liblapack.a)
        set(LAPACK_LIBRARIES ${LAPACK_LIBRARIES} ${LAPACK_LIB_DIR}/liblapack.a)
    elseif (EXISTS ${LAPACK_LIB_DIR}/libreflapack.a)
        set(LAPACK_LIBRARIES ${LAPACK_LIBRARIES} ${LAPACK_LIB_DIR}/libreflapack.a)
    else()
        message(WARNING "Could not find library lapack or reflapack in " ${LAPACK_LIB_DIR} ". Compiling Hdef without LAPACK.")
        set(useLAPACK 0)
    endif()

    # Find blaslib
    if (EXISTS ${LAPACK_LIB_DIR}/libblas.a)
        set(LAPACK_LIBRARIES ${LAPACK_LIBRARIES} ${LAPACK_LIB_DIR}/libblas.a)
    elseif (EXISTS ${LAPACK_LIB_DIR}/librefblas.a)
        set(LAPACK_LIBRARIES ${LAPACK_LIBRARIES} ${LAPACK_LIB_DIR}/librefblas.a)
    else()
        message(FATAL_ERROR "Could not find library blas or refblas in " ${LAPACK_LIB_DIR} ". Compiling Hdef without LAPACK.")
        set(useLAPACK 0)
    endif()

    # Find tmglib
    if (EXISTS ${LAPACK_LIB_DIR}/libtmglib.a)
        set(LAPACK_LIBRARIES ${LAPACK_LIBRARIES} ${LAPACK_LIB_DIR}/libtmglib.a)
    elseif (EXISTS ${LAPACK_LIB_DIR}/libreftmglib.a)
        set(LAPACK_LIBRARIES ${LAPACK_LIBRARIES} ${LAPACK_LIB_DIR}/libreftmglib.a)
    else()
        message(FATAL_ERROR "Could not find library tmglib or reftmglib in " ${LAPACK_LIB_DIR} ". Compiling Hdef without LAPACK.")
        set(useLAPACK 0)
    endif()

else()
    message(WARNING "Could not find LAPACK directory: '" ${LAPACK_LIB_DIR} "'. Compiling Hdef without LAPACK")
endif()

if (useLAPACK)
    foreach (lib ${LAPACK_LIBRARIES})
        message(STATUS "Using LAPACK library: " ${lib} )
    endforeach()
    add_compile_definitions(USE_LAPACK)
endif()
message(STATUS "*********************************************************************************")
message(STATUS "")




#####
#	LIBRARIES
#####
message(STATUS "*********************************************************************************")
message(STATUS "Setting up libhdef build")

add_compile_definitions(UNIT_TEST)

# Build everything together in single hdef library
add_library(hdef
            src/algebra_module.f90
            src/annealing_module.f90
            src/calendar_module.f90
            src/earth_module.f90
            src/elast_module.f90
            src/eq_module.f90
            src/error_exit.c
            src/ffm_module.f90
            src/flexure_module.f90
            src/geom_module.f90
            src/interpolation_module.f90
            src/io_module.f90
            src/map_module.f90
            src/misfit_module.f90
            src/okada92_module.f90
            src/random_module.f90
            src/solver_module.f90
            src/test_module.f90
            src/trig_module.f90
            src/tri_disloc_module.f90
)

message(STATUS "*********************************************************************************")
message(STATUS "")


#####
#	EXECUTABLES
#####
message(STATUS "*********************************************************************************")
message(STATUS "Setting up executable builds")

# Modern Fortran files use libraries
list(APPEND f90_files
     anneal_post
     clip
     colortool
     dateutil
     distaz2lola
     ff2gmt
     fitutil
     grid
     interpolate
     lola2distaz
     mtutil
     o92util
     perspective
     platemotion
     rangen
     readGCMT
     sphfinrot
     stereo_project
     tri_tool
     triutil
     vec2los
     wraplos
)
foreach(executable ${f90_files})
    add_executable(${executable} src/${executable}.f90)
    add_dependencies(${executable} hdef)
    target_link_libraries(${executable} hdef ${LAPACK_LIBRARIES})
endforeach()

# fltinv is divided into multiple source files, so build it separately
add_executable(fltinv src/fltinv.f90
                      src/fltinv_io.f90
                      src/fltinv_gf.f90
                      src/fltinv_lsqr.f90
                      src/fltinv_anneal.f90
                      src/fltinv_anneal_psc.f90
                      src/fltinv_anneal_psc_euler.f90)
add_dependencies(fltinv hdef)
target_link_libraries(fltinv hdef ${LAPACK_LIBRARIES})

# Old-school Fortran files need only themselves
list(APPEND f77_files
     numint
     readkik
)
foreach(executable ${f77_files})
    add_executable(${executable} src/${executable}.f)
endforeach()


# Unit tests
enable_testing()
list(APPEND unit_test_files
     unit_test_annealing
     unit_test_calendar
     unit_test_earth
     unit_test_eq
     unit_test_ffm
     unit_test_geom
     unit_test_map
     unit_test_okada92
     unit_test_random
     unit_test_solver
     unit_test_tri_disloc
)
foreach(executable ${unit_test_files})
    add_executable(${executable} src/${executable}.f90)
    add_test(${executable} ${executable})
    add_dependencies(${executable} hdef)
    target_link_libraries(${executable} hdef ${LAPACK_LIBRARIES})
    target_compile_options(${executable} PUBLIC -ffree-line-length-0)
endforeach()


message(STATUS "*********************************************************************************")
message(STATUS "")


#####
#	SCRIPTS
#####
message(STATUS "*********************************************************************************")
message(STATUS "Setting up script builds")

# add_custom_target(coul_dip ALL DEPENDS coul_dip.sh)
# add_custom_command(
#     OUTPUT coul_dip.sh
#     COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/scripts/coul_dip.sh ${CMAKE_CURRENT_BINARY_DIR}/coul_dip.sh
#     MAIN_DEPENDENCY ${CMAKE_SOURCE_DIR}/scripts/coul_dip.sh
#     DEPENDS ff2gmt grid o92util
#     COMMENT "Building target coul_dip.sh")

add_custom_target(surf_disp ALL DEPENDS surf_disp.sh)
add_custom_command(
    OUTPUT surf_disp.sh
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/scripts/surf_disp.sh ${CMAKE_CURRENT_BINARY_DIR}/surf_disp.sh
    MAIN_DEPENDENCY ${CMAKE_SOURCE_DIR}/scripts/surf_disp.sh
    DEPENDS ff2gmt grid o92util
    COMMENT "Building target surf_disp.sh")

message(STATUS "*********************************************************************************")
message(STATUS "")
