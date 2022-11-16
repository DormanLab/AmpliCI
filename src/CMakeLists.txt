cmake_minimum_required(VERSION 3.3.0)

project(AmpliCI VERSION 2.1)

set(CMAKE_C_COMPILER gcc)

# run cmake -S . -B debug -DCMAKE_BUILD_TYPE=Debug && cmake --build debug to compile debugging version
#set(CMAKE_BUILD_TYPE release)
set(CMAKE_C_STANDARD 11)

include(FindPackageHandleStandardArgs)

# loess source code is in Fortran
enable_language(Fortran)
set(CMAKE_FORTRAN_COMPILER gfortran)

# below modified from RStudio's FindLibR.cmake to find R include and lib path
# https://github.com/rstudio/rstudio/blob/master/cmake/modules/FindLibR.cmake
#if(APPLE)  # mac
#
#   find_library(LIBR_LIBRARIES R)
#
#   if(LIBR_LIBRARIES MATCHES ".*\\.framework")
#     set(LIBR_HOME "${LIBR_LIBRARIES}/Resources" CACHE PATH "R home directory")
#     set(LIBR_INCLUDE_DIRS "${LIBR_HOME}/include" CACHE PATH "R include directory")
#     set(LIBR_LIB_DIR "${LIBR_HOME}/lib" CACHE PATH "R lib directory")
#     set(LIBR_EXECUTABLE "${LIBR_HOME}/bin/R" CACHE PATH "R executable")
#   else()
#      get_filename_component(_LIBR_LIBRARIES "${LIBR_LIBRARIES}" REALPATH)
#      get_filename_component(_LIBR_LIBRARIES_DIR "${_LIBR_LIBRARIES}" PATH)
#      set(LIBR_EXECUTABLE "${_LIBR_LIBRARIES_DIR}/../bin/R")
#      execute_process(
#         COMMAND ${LIBR_EXECUTABLE} "--slave" "--vanilla" "-e" "cat(R.home())"
#                   OUTPUT_VARIABLE LIBR_HOME
#      )
#      set(LIBR_HOME ${LIBR_HOME} CACHE PATH "R home directory")
#      set(LIBR_INCLUDE_DIRS "${LIBR_HOME}/include" CACHE PATH "R include directory")
#      set(LIBR_LIB_DIR "${LIBR_HOME}/lib" CACHE PATH "R lib directory")
#   endif()
#
#else()  # UNIX
#   # find executable
#   find_program(LIBR_EXECUTABLE R)
#   if(LIBR_EXECUTABLE-NOTFOUND)
#      message(STATUS "Unable to locate R executable")
#   endif()
#   
#   # ask R for the home path
#   if(NOT LIBR_HOME)
#      execute_process(
#         COMMAND ${LIBR_EXECUTABLE} "--slave" "--vanilla" "-e" "cat(R.home())"
#                   OUTPUT_VARIABLE LIBR_HOME
#      )
#      if(LIBR_HOME)
#        set(LIBR_HOME ${LIBR_HOME} CACHE PATH "R home directory")
#      endif()
#   endif()
#
#   # ask R for the include dir
#   if(NOT LIBR_INCLUDE_DIRS)
#      execute_process(
#         COMMAND ${LIBR_EXECUTABLE} "--slave" "--no-save" "-e" "cat(R.home('include'))"
#         OUTPUT_VARIABLE LIBR_INCLUDE_DIRS
#      )
#      if(LIBR_INCLUDE_DIRS)
#        set(LIBR_INCLUDE_DIRS ${LIBR_INCLUDE_DIRS} CACHE PATH "R include directory")
#      endif()
#   endif()
#
#   # ask R for the lib dir
#   if(NOT LIBR_LIB_DIR)
#      execute_process(
#         COMMAND ${LIBR_EXECUTABLE} "--slave" "--no-save" "-e" "cat(R.home('lib'))"
#         OUTPUT_VARIABLE LIBR_LIB_DIR
#      )
#      if(LIBR_LIB_DIR)
#         set(LIBR_LIB_DIR ${LIBR_LIB_DIR} CACHE PATH "R lib directory")
#      endif()
#   endif()
#
#
#   # sometimes lib and include are not under R.home()
#   get_filename_component(_LIBR_EXECUTABLE "${LIBR_EXECUTABLE}" PATH)
#   get_filename_component(_LIBR_EXECUTABLE_DIR "${_LIBR_EXECUTABLE}" PATH)
#   #message(STATUS "Found R upper directory: ${_LIBR_EXECUTABLE_DIR}")
#
#endif()
#
#message(STATUS "Found R: ${LIBR_EXECUTABLE}")
#message(STATUS "Found R home: ${LIBR_HOME}")
#message(STATUS "Found R include directory: ${LIBR_INCLUDE_DIRS}")
#message(STATUS "Found R lib directory: ${LIBR_LIB_DIR}")
#message(STATUS "Found R executable directory: ${_LIBR_EXECUTABLE_DIR}")
#
#include_directories(${LIBR_INCLUDE_DIRS})
#link_directories(${LIBR_LIB_DIR} ${RMATH_PATH})

## loess library
set(F_OBJ ./libblas_lite.f ./linpack_lite.f ./loessf.f)
set(C_OBJ_LOESS ./loessc.c)

add_library(loess STATIC ${F_OBJ} ${C_OBJ_LOESS})
#target_include_directories(loess PUBLIC ${LIBR_INCLUDE_DIRS})

# prefer to use libgfortran installed with devel package first, otherwise use compiler default.
find_library(GFORTRAN_LIBRARY gfortran)
if(NOT GFORTRAN_LIBRARY)
   message(STATUS "Using default ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES}")
else()
   ## cmake with fortran will automatically link to gfortran and quadmath. Now unset it.
   unset(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES)
   target_link_libraries(loess ${GFORTRAN_LIBRARY})
endif()

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

#target_include_directories(amplici PUBLIC ${LIBR_INCLUDE_DIRS})

# Rmath library
#set(RMATH_PATH ".")
#find_library(RMATH_LIBRARY NAMES Rmath rmath HINTS ${RMATH_PATH} ${LIBR_LIB_DIR} ${_LIBR_EXECUTABLE_DIR} PATH_SUFFIXES lib lib64)
#find_package_handle_standard_args(LibRmath DEFAULT_MSG RMATH_LIBRARY)
#
#if(RMATH_LIBRARY-NOTFOUND)
#   message(FATAL_ERROR "Failed to find libRmath!")
#else()
#   message(STATUS "Found Rmath library: ${RMATH_LIBRARY}")
#   target_link_libraries(amplici ${RMATH_LIBRARY})
#endif()

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
