cmake_minimum_required(VERSION 3.5.0)

project(AmpliCI VERSION 2.2)

#set(CMAKE_C_COMPILER gcc)

# run cmake -S . -B debug -DCMAKE_BUILD_TYPE=Debug && cmake --build debug to compile debugging version
#set(CMAKE_BUILD_TYPE release)
set(CMAKE_C_STANDARD 11)

include(FindPackageHandleStandardArgs)

# loess source code is in Fortran
enable_language(Fortran)
#set(CMAKE_FORTRAN_COMPILER gfortran)

## loess library
set(F_OBJ ./libblas_lite.f ./linpack_lite.f ./loessf.f)
set(C_OBJ_LOESS ./loessc.c)

add_library(loess STATIC ${F_OBJ} ${C_OBJ_LOESS})
#target_include_directories(loess PUBLIC ${LIBR_INCLUDE_DIRS})

# prefer to use libgfortran installed with devel package first, otherwise use compiler default.
#find_library(GFORTRAN_LIBRARY gfortran)
#if(NOT GFORTRAN_LIBRARY)
#   message(STATUS "Using default ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES}")
#else()
   ## cmake with fortran will automatically link to gfortran and quadmath. Now unset it.
unset(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES)
   #target_link_libraries(loess ${GFORTRAN_LIBRARY})
   #endif()

# C compiler options
add_compile_options(-std=c99 -D_POSIX_C_SOURCE=200112L -D_XOPEN_SOURCE=500)
add_compile_options(-Wall -Wextra -pedantic)

## amplici library
set(C_OBJ_AMPLICI ./amplici.c ./ampliclust.c ./align.c ./cmdline.c ./data.c
   ./error.c ./fastq.c ./hash.c ./initialize.c ./io.c ./libamplici.c
   ./lmath.c ./model.c ./options.c ./statistics.c ./util.c ./amplici_umi.c
   ./partition.c ./ppois.c ./dpois.c ./pnorm.c ./dnorm.c ./stirlerr.c ./bd0.c
   ./pgamma.c ./fmax2.c ./lgamma.c ./gamma.c ./lgammacor.c ./chebyshev.c
   ./cospi.c ./d1mach.c)

add_library(amplici SHARED ${C_OBJ_AMPLICI})
add_library(amplici_static STATIC ${C_OBJ_AMPLICI})


target_link_libraries(amplici -lm)

## amplici executable
set(C_OBJ ./loess.c ./predict.c ./error_est.c ./run_amplici.c)
add_executable(run_AmpliCI ${C_OBJ})
add_dependencies(run_AmpliCI loess)
add_dependencies(run_AmpliCI amplici)
target_link_libraries(run_AmpliCI PRIVATE loess amplici)


### change the name of the static library 
set_target_properties(amplici_static PROPERTIES OUTPUT_NAME "amplici")
### allow ouput both static and shared library 
set_target_properties(amplici PROPERTIES CLEAN_DIRECT_OUTPUT 1)
set_target_properties(amplici_static PROPERTIES CLEAN_DIRECT_OUTPUT 1)

install(TARGETS run_AmpliCI amplici amplici_static
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib)

# Set the RPATH to the directory where the shared library is installed
set(CMAKE_INSTALL_RPATH "@executable_path/../lib")
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)

set_target_properties(run_AmpliCI PROPERTIES INSTALL_RPATH "@executable_path/../lib")
set_target_properties(amplici PROPERTIES INSTALL_RPATH "@loader_path")

if(NOT DEFINED CMAKE_INSTALL_PREFIX)
    set(CMAKE_INSTALL_PREFIX $ENV{PREFIX})
endif()
