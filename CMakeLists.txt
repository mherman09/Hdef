cmake_minimum_required(VERSION 3.12)
# cmake -DCMAKE_INSTALL_PREFIX=`pwd` -DCMAKE_C_COMPILER=gcc-8 -DCMAKE_Fortran_COMPILER=gfortran -DLAPACK_LIB_DIR=/sw/lib/lapack ..

project(Hdef VERSION 2020.03)

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
set(CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS} -cpp)

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







# message(STATUS "CMAKE_INSTALL_PREFIX=" ${CMAKE_INSTALL_PREFIX})
# install(TARGETS hdef DESTINATION lib)
# # message(STATUS "${PROJECT_SOURCE_DIR}")
# # message(STATUS "${PROJECT_BINARY_DIR}")
# install(PROGRAMS src/surf_disp.sh DESTINATION bin)
# add_executable(grid src/grid.f90)
# add_dependencies(grid algebra earth geom io trig)
# target_link_libraries(grid
#                       $<TARGET_OBJECTS:algebra>
#                       $<TARGET_OBJECTS:geom>
#                       $<TARGET_OBJECTS:io>)
#
# add_executable(fltinv src/fltinv.f90
#                       src/fltinv_io.f90
#                       src/fltinv_gf.f90
#                       src/fltinv_lsqr.f90
#                       src/fltinv_anneal.f90
#                       src/fltinv_anneal_psc.f90
#                       src/fltinv_anneal_psc_euler.f90)
# add_dependencies(fltinv algebra annealing earth elast error_exit geom io misfit okada92 random solver trig tri_disloc)
# target_link_libraries(fltinv
#                       $<TARGET_OBJECTS:algebra>
#                       $<TARGET_OBJECTS:annealing>
#                       $<TARGET_OBJECTS:earth>
#                       $<TARGET_OBJECTS:elast>
#                       $<TARGET_OBJECTS:error_exit>
#                       $<TARGET_OBJECTS:geom>
#                       $<TARGET_OBJECTS:io>
#                       $<TARGET_OBJECTS:misfit>
#                       $<TARGET_OBJECTS:okada92>
#                       $<TARGET_OBJECTS:random>
#                       $<TARGET_OBJECTS:solver>
#                       $<TARGET_OBJECTS:tri_disloc>
#                       ${LAPACK_LIBRARIES}
# )
# list(APPEND CMAKE_Fortran_FLAGS_DEBUG ${bounds})
# list(APPEND CMAKE_Fortran_FLAGS ${dialect})

# target_include_directories(Tutorial PUBLIC
#                            "${PROJECT_BINARY_DIR}"
#                            )



####################################################################################################
#
# PRINT OUT TARGET OR DIRECTORY COMPILE_DEFINITION VARIABLES
#
# get_target_property( thisProperty fltinv COMPILE_DEFINITIONS )
# message(STATUS "fltinv COMPILE_DEFINITIONS=" ${thisProperty})
#
# get_directory_property(directory_definitions COMPILE_DEFINITIONS )
# foreach (def ${directory_definitions})
#     message(STATUS "Compiling with: " ${def} )
# endforeach()
#
####################################################################################################
#
# BUILD ALL OBJECT FILES SEPARATELY
#
# add_library(algebra OBJECT src/algebra_module.f90)
# add_dependencies(algebra io trig)
#
# add_library(annealing OBJECT src/annealing_module.f90)
# add_dependencies(annealing io random)
#
# add_library(calendar OBJECT src/calendar_module.f90)
# add_dependencies(calendar io)
#
# add_library(earth OBJECT src/earth_module.f90)
# add_dependencies(earth algebra io trig)
#
# add_library(elast OBJECT src/elast_module.f90)
# add_dependencies(elast algebra geom io)
#
# add_library(eq OBJECT src/eq_module.f90)
# add_dependencies(eq io trig)
#
# add_library(error_exit OBJECT src/error_exit.c)
# #add_dependencies(error_exit)
#
# add_library(ffm OBJECT src/ffm_module.f90)
# add_dependencies(ffm geom io)
#
# add_library(flexure OBJECT src/flexure_module.f90)
# add_dependencies(flexure trig)
#
# add_library(geom OBJECT src/geom_module.f90)
# add_dependencies(geom algebra io trig)
#
# add_library(io OBJECT src/io_module.f90)
# #add_dependencies(io)
#
# add_library(map OBJECT src/map_module.f90)
# add_dependencies(map earth io trig)
#
# add_library(misfit OBJECT src/misfit_module.f90)
# add_dependencies(misfit io solver)
#
# add_library(okada92 OBJECT src/okada92_module.f90)
# add_dependencies(okada92 io trig)
#
# add_library(random OBJECT src/random_module.f90)
# #add_dependencies(random)
#
# add_library(solver OBJECT src/solver_module.f90)
# add_dependencies(solver io)
#
# add_library(test OBJECT src/test_module.f90)
# #add_dependencies(test)
#
# add_library(tri_disloc OBJECT src/tri_disloc_module.f90)
# add_dependencies(tri_disloc earth geom io trig)
#
# add_library(trig OBJECT src/trig_module.f90)
# #add_dependencies(trig)
