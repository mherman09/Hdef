cmake_minimum_required(VERSION 3.12)
# 	cmake \
#		-DCMAKE_INSTALL_PREFIX=/path/to/hdef/installation/directory \
# 		-DCMAKE_C_COMPILER=/path/to/gcc \
#		-DCMAKE_Fortran_COMPILER=/path/to/gfortran \
#		-DLAPACK_LIB_DIR=/path/to/lapack \
#		..

project(Hdef VERSION 2021.03)

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

# Make sure that compilers are using same version
if (NOT ${CMAKE_Fortran_COMPILER_VERSION} MATCHES ${CMAKE_C_COMPILER_VERSION})
    message(STATUS "")
    message(WARNING "Compiler versions do not match; this may lead to unanticipated behavior")
    message(STATUS "CMAKE_C_COMPILER_VERSION=" ${CMAKE_C_COMPILER_VERSION})
    message(STATUS "CMAKE_Fortran_COMPILER_VERSION=" ${CMAKE_Fortran_COMPILER_VERSION})
    message(STATUS "")
else()
    message(STATUS "Compiler versions match: ${CMAKE_Fortran_COMPILER_VERSION}")
endif()

# Use pre-processing commands
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp")

# Compilation flags
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -Wextra -Wunused -Wconversion -fcheck=bounds,do,mem,pointer,recursion -fbacktrace")

message(STATUS "*********************************************************************************")
message(STATUS "")




#####
#	LAPACK
#####
message(STATUS "*********************************************************************************")
message(STATUS "Looking for LAPACK")


# Default: do not use LAPACK
set(useLAPACK 0)


# Define LAPACK libraries
if ("${LAPACK_LIB_DIR}" STREQUAL "" )
    # If LAPACK directory is not set on the command line
    message(STATUS "Compiling Hdef ***without*** LAPACK")
    message(STATUS "Set LAPACK directory with: '-DLAPACK_LIB_DIR=/path/to/lapack/lib'")

elseif (EXISTS ${LAPACK_LIB_DIR})
    # If LAPACK directory is set on the command line and exists
    message(STATUS "Found LAPACK directory: " ${LAPACK_LIB_DIR})
    set(useLAPACK 1)
    unset(LAPACK_LIBRARIES)

    # Find lapack lib
    if (EXISTS ${LAPACK_LIB_DIR}/liblapack.a)
        set(LAPACK_LIBRARIES ${LAPACK_LIBRARIES} ${LAPACK_LIB_DIR}/liblapack.a)
    elseif (EXISTS ${LAPACK_LIB_DIR}/libreflapack.a)
        set(LAPACK_LIBRARIES ${LAPACK_LIBRARIES} ${LAPACK_LIB_DIR}/libreflapack.a)
    else()
        message(WARNING "Could not find library lapack or reflapack in " ${LAPACK_LIB_DIR} ". Compiling Hdef ***without*** LAPACK.")
        set(useLAPACK 0)
    endif()

    # Find blaslib
    if (EXISTS ${LAPACK_LIB_DIR}/libblas.a)
        set(LAPACK_LIBRARIES ${LAPACK_LIBRARIES} ${LAPACK_LIB_DIR}/libblas.a)
    elseif (EXISTS ${LAPACK_LIB_DIR}/librefblas.a)
        set(LAPACK_LIBRARIES ${LAPACK_LIBRARIES} ${LAPACK_LIB_DIR}/librefblas.a)
    else()
        message(FATAL_ERROR "Could not find library blas or refblas in " ${LAPACK_LIB_DIR} ". Compiling Hdef ***without*** LAPACK.")
        set(useLAPACK 0)
    endif()

    # Find tmglib
    if (EXISTS ${LAPACK_LIB_DIR}/libtmglib.a)
        set(LAPACK_LIBRARIES ${LAPACK_LIBRARIES} ${LAPACK_LIB_DIR}/libtmglib.a)
    elseif (EXISTS ${LAPACK_LIB_DIR}/libreftmglib.a)
        set(LAPACK_LIBRARIES ${LAPACK_LIBRARIES} ${LAPACK_LIB_DIR}/libreftmglib.a)
    else()
        message(FATAL_ERROR "Could not find library tmglib or reftmglib in " ${LAPACK_LIB_DIR} ". Compiling Hdef ***without*** LAPACK.")
        set(useLAPACK 0)
    endif()

else()
    # If LAPACK directory is set on the command line, but does not exist
    message(WARNING "Could not find LAPACK directory: '" ${LAPACK_LIB_DIR} "'. Compiling Hdef ***without*** LAPACK")
endif()


# Print message if LAPACK is found
if (useLAPACK)
    foreach (lib ${LAPACK_LIBRARIES})
        message(STATUS "Using LAPACK library: " ${lib} )
    endforeach()
    add_compile_definitions(USE_LAPACK)
endif()

message(STATUS "*********************************************************************************")
message(STATUS "")




#####
#	OPENMP
#####
message(STATUS "*********************************************************************************")
message(STATUS "Looking for OpenMP libraries")

# Look for package OpenMP
find_package(OpenMP)

# Print found message during configuration
if (OpenMP_Fortran_FOUND)
    message(STATUS "Found OpenMP_Fortran_LIBRARIES: ${OpenMP_Fortran_LIBRARIES}")
    message(STATUS "Building o92util with parallel capability")
    add_compile_definitions(USE_OPENMP)
else()
    message(STATUS "Could not find package OpenMP")
    message(STATUS "Building o92util without parallel capability")
endif()

message(STATUS "*********************************************************************************")
message(STATUS "")




#####
#	LIBRARIES
#####
message(STATUS "*********************************************************************************")
message(STATUS "Setting up libhdef build")

# Compile with UNIT_TEST preprocessing directives
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
            src/geom_module.f90
            src/interpolation_module.f90
            src/io_module.f90
            src/map_module.f90
            src/misfit_module.f90
            src/okada92_module.f90
            src/random_module.f90
            src/solver_module.f90
            src/sort_module.f90
            src/test_module.f90
            src/trig_module.f90
            src/tri_disloc_module.f90
)


# Compile with OpenMP if it is found
if (OpenMP_Fortran_FOUND)
    # Fortran subroutines need to be recursive to work with OpenMP
    # The "-frecursive" option does not apply to C subroutines, so only add to Fortran flags
    # target_compile_options(hdef PRIVATE "-frecursive")
    set(APPEND CMAKE_Fortran_FLAGS "-frecursive")

    # Link the hdef library with OpenMP libraries
    target_link_libraries(hdef PRIVATE OpenMP::OpenMP_Fortran)
endif()


# Install the hdef library
install(TARGETS hdef DESTINATION lib)

message(STATUS "*********************************************************************************")
message(STATUS "")




#####
#	EXECUTABLES
#####
message(STATUS "*********************************************************************************")
message(STATUS "Setting up executable builds")

# List of modern Fortran files (.f90 extensions)
list(APPEND f90_files
     anneal_post
     clip
     colortool
     dateutil
     distaz2lola
     eqempirical
     ff2gmt
     fitutil
     grid
     interpolate
     lola2distaz
     mtutil
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

# These are built with the hdef library
foreach(executable ${f90_files})
    add_executable(${executable} src/${executable}.f90)
    add_dependencies(${executable} hdef)
    target_link_libraries(${executable} hdef ${LAPACK_LIBRARIES})
endforeach()

# If OpenMP is available, create parallel version of o92util
add_executable(o92util src/o92util.f90)
add_dependencies(o92util hdef)
target_link_libraries(o92util hdef ${LAPACK_LIBRARIES})
if (OpenMP_Fortran_FOUND)
    target_link_libraries(o92util OpenMP::OpenMP_Fortran)
    target_compile_options(o92util PUBLIC "-frecursive")
endif()
install(TARGETS o92util DESTINATION bin)

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

# Old-school Fortran files (.f extensions) need only themselves
list(APPEND f77_files
     numint
     readkik
)
foreach(executable ${f77_files})
    add_executable(${executable} src/${executable}.f)
    target_compile_options(${executable} PRIVATE "-frecursive")
endforeach()

# Install compiled executables in bin directory
list(APPEND hdef_executables ${f90_files} ${f77_files} fltinv)
foreach(executable ${hdef_executables})
    install(TARGETS ${executable} DESTINATION bin)
endforeach()

message(STATUS "*********************************************************************************")
message(STATUS "")




#####
#	SCRIPTS
#####
message(STATUS "*********************************************************************************")
message(STATUS "Setting up script builds")

# Shell scripts that use executables
list(APPEND shell_scripts
     coul_dip
     coul_hor
     coul_xsec
     eqmatch
     gmtcpt
     insar
     simplify_ffm
     surf_disp
     ternary
     trg_schem
     make_seis_movie
)

# To build and install these, define them as custom targets
foreach(s ${shell_scripts})
    add_custom_target(${s}-script ALL DEPENDS ${s}.sh)
    add_custom_command(
        OUTPUT ${s}.sh
        COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/scripts/${s}.sh ${CMAKE_CURRENT_BINARY_DIR}/${s}.sh
        MAIN_DEPENDENCY ${CMAKE_SOURCE_DIR}/scripts/${s}.sh
        DEPENDS ${hdef_executables}
        COMMENT "Building target ${s}.sh"
    )
    install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${s}.sh DESTINATION bin)
endforeach()

# Perl script for parsing ComCat CSV file
add_custom_target(parse_usgs_query-script ALL DEPENDS parse_usgs_query.pl)
add_custom_command(
    OUTPUT parse_usgs_query.pl
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/scripts/parse_usgs_query.pl ${CMAKE_CURRENT_BINARY_DIR}/parse_usgs_query.pl
    MAIN_DEPENDENCY ${CMAKE_SOURCE_DIR}/scripts/parse_usgs_query.pl
    COMMENT "Building target parse_usgs_query.pl"
)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/parse_usgs_query.pl DESTINATION bin)


message(STATUS "*********************************************************************************")
message(STATUS "")




#####
#	TESTS
#####
message(STATUS "*********************************************************************************")
message(STATUS "Setting up tests")

# Activate cmake testing in build directory
enable_testing()

# Unit tests
if (useLAPACK)
    list(APPEND unit_test_list
        unit_test_algebra
        unit_test_annealing
        unit_test_calendar
        unit_test_earth
        unit_test_elast
        unit_test_eq
        unit_test_ffm
        unit_test_geom
        unit_test_interpolation
        unit_test_io
        unit_test_map
        unit_test_misfit
        unit_test_okada92
        unit_test_random
        unit_test_solver
        unit_test_sort
        unit_test_tri_disloc
    )
else()
    list(APPEND unit_test_list
        unit_test_algebra
        unit_test_annealing
        unit_test_calendar
        unit_test_earth
        unit_test_elast
        unit_test_eq
        unit_test_ffm
        unit_test_geom
        unit_test_interpolation
        unit_test_io
        unit_test_map
        unit_test_okada92
        unit_test_random
        unit_test_sort
        unit_test_tri_disloc
    )
endif()

foreach(unit_test ${unit_test_list})
    # Make the unit test an executable and define as a test
    add_executable(${unit_test} src/${unit_test}.f90)
    add_test(${unit_test} ${unit_test})

    # All unit tests depend on hdef and some need link to LAPACK
    add_dependencies(${unit_test} hdef)
    target_link_libraries(${unit_test} hdef ${LAPACK_LIBRARIES})

    # Some of the unit tests have very long lines, so make compiler read unlimited length
    target_compile_options(${unit_test} PUBLIC -ffree-line-length-0)

    # # Install the unit tests in the test/ folder
    # install(TARGETS ${unit_test} DESTINATION test)
endforeach()


# Executable tests
if (useLAPACK)
    list(APPEND executable_test_list
        test_clip
        test_dateutil
        test_distaz2lola
        test_ff2gmt
        test_fitutil
        test_fltinv
        test_grid
        test_lola2distaz
        test_mtutil
        test_o92util
        test_platemotion
        test_rangen
        test_readGCMT
        test_sphfinrot
        test_stereo_project
        test_triutil
        test_vec2los
        test_wraplos
    )
else()
    list(APPEND executable_test_list
        test_clip
        test_dateutil
        test_distaz2lola
        test_ff2gmt
        test_grid
        test_lola2distaz
        test_mtutil
        test_o92util
        test_platemotion
        test_rangen
        test_readGCMT
        test_sphfinrot
        test_stereo_project
        test_triutil
        test_vec2los
        test_wraplos
    )
endif()
foreach(exec_test ${executable_test_list})
    # Define all executable tests in the main test target
    add_test(NAME ${exec_test} COMMAND ${CMAKE_SOURCE_DIR}/test/${exec_test}.sh)

    # # Install the executable tests in the test/ folder
    # install(PROGRAMS ${CMAKE_SOURCE_DIR}/test/${exec_test}.sh DESTINATION test)

    # Define a command to test executable directly (used in target test-exec)
    string(REGEX REPLACE test_ "" exec ${exec_test})
    add_custom_command(OUTPUT ${exec_test}
        COMMAND ${CMAKE_SOURCE_DIR}/test/${exec_test}.sh
        DEPENDS ${exec}
        COMMENT "Testing ${exec}"
    )
endforeach()
# Define a custom test target to test all hdef executables
add_custom_target(test-exec DEPENDS ${executable_test_list})

# # Also install script test_values.sh to check output values
# install(PROGRAMS ${CMAKE_SOURCE_DIR}/test/test_values.sh DESTINATION test)

message(STATUS "*********************************************************************************")
message(STATUS "")




#####
#	FINISH UP
#####
message(STATUS "*********************************************************************************")
message(STATUS "Finished configuring Hdef")
message(STATUS "")


# Add a message to print about LAPACK during build
if(useLAPACK)
    message(STATUS "Building Hdef with LAPACK")
    message(STATUS "")
    add_custom_target(hdef-lapack-check ALL DEPENDS hdef ${hdef_executables} COMMENT "Hdef built with LAPACK")
else()
    message(STATUS "Building Hdef without LAPACK")
    message(STATUS "")
    add_custom_target(hdef-lapack-check ALL DEPENDS hdef ${hdef_executables}
                      COMMENT "Hdef built without LAPACK: inversion and fitting programs will not work.
       To build with LAPACK, run cmake with '-DLAPACK_LIB_DIR=/path/to/lapack/lib'")
endif()


# Add a message to print about OpenMP during build
if(OpenMP_Fortran_FOUND)
    message(STATUS "Building Hdef with OpenMP")
    message(STATUS "")
    add_custom_target(hdef-openmp-check ALL DEPENDS hdef ${hdef_executables} COMMENT "Hdef built with OpenMP")
else()
    message(STATUS "Building Hdef without OpenMP")
    message(STATUS "")
    add_custom_target(hdef-openmp-check ALL DEPENDS hdef ${hdef_executables}
                      COMMENT "Hdef built without OpenMP: no parallel capabilities enabled")
endif()


# Print preprocessing COMPILE_DEFINITIONS
message(STATUS "Compiling Hdef with preprocessor definitions:")
get_directory_property(DIRECTORY_COMPILE_DEFINITIONS DIRECTORY ${CMAKE_SOURCE_DIR} COMPILE_DEFINITIONS)
foreach(DEF ${DIRECTORY_COMPILE_DEFINITIONS})
    message(STATUS "    ${DEF}")
endforeach()
message(STATUS "")


# Print preprocessing COMPILE_DEFINITIONS
message(STATUS "Using Fortran compiler flags: ${CMAKE_Fortran_FLAGS}")
message(STATUS "Using C compiler flags: ${CMAKE_C_FLAGS}")
message(STATUS "")


# Print simple directions
message(STATUS "Type 'make' to compile Hdef (in this directory)")
message(STATUS "Type 'make test' to run unit tests")
message(STATUS "Type 'make test-exec' to run executable tests")
message(STATUS "Type 'make install' to install Hdef in directory '${CMAKE_INSTALL_PREFIX}'")
message(STATUS "    Executables   -> ${CMAKE_INSTALL_PREFIX}/bin")
message(STATUS "    Libraries     -> ${CMAKE_INSTALL_PREFIX}/lib")
message(STATUS "    Tests         -> ${CMAKE_INSTALL_PREFIX}/test")
message(STATUS "*********************************************************************************")
message(STATUS "")
