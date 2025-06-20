cmake_minimum_required(VERSION 2.8...3.19)
project(MEX)

set(CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH}
  ${PROJECT_SOURCE_DIR}/cmake)

if(CMAKE_VERSION VERSION_LESS "3.18")
  message(WARNING "CMake version less than 3.18 detected"
    " -- using local FindMatlab module extracted from v3.18")
  # Where we keep FindMatlab module
  set(CMAKE_MODULE_PATH
    ${PROJECT_SOURCE_DIR}/cmake-compat
    ${CMAKE_MODULE_PATH})
endif()

find_package(Matlab REQUIRED MAIN_PROGRAM)
if(NOT Matlab_FOUND)
  message(FATAL_ERROR "MATLAB not found.")
else()
  message("MATLAB found.")
endif()
include_directories(
  ${Matlab_INCLUDE_DIRS}
)

find_package(GSL REQUIRED)
if(NOT GSL_FOUND)
  message(FATAL_ERROR "GSL not found.")
else()
  message("GSL found.")
endif()

find_package(MPFR)
if(MPFR_FOUND)
  message("MPFR found.")
else()
  message("MPFR not found.")
endif()

find_package(MKL)
if(MKL_FOUND)
  include_directories(
    ${MKL_INCLUDE_DIR}
    )
  link_directories(
    ${MKL_LIBRARY}
    )
endif()

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin/")
if ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  # -ffast-math is for some reason not working with gcc
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -std=c99 -Wall -O3 -march=native -fopenmp")
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
  link_directories(/NOBACKUP/davoudss/intel/lib/intel64/)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -std=c99 -Wall -fast -ipo -O3 -march=native -openmp")
endif()

# Add modules with MEX to be built
add_subdirectory("${PROJECT_SOURCE_DIR}/SE_direct/mex")

