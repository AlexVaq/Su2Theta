# basic setup for cmake
cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_INCLUDE_DIRECTORIES_PROJECT_BEFORE ON)
set(CMAKE_COLOR_MAKEFILE ON)
set(CMAKE_COLOR_MAKEFILE ON)
set(CMAKE_CXX_STANDARD_REQUIRED 17)
set(CMAKE_CXX_STANDARD 17)

set(DEFAULT_BUILD "RELEASE")
set(VALID_BUILDS DEVEL RELEASE DEBUG)
set(CMAKE_BUILD_TYPE "${DEFAULT_BUILD}" CACHE STRING  "Choose the type of build, options are: ${VALID_BUILDS}")
string(TOUPPER "${CMAKE_BUILD_TYPE}" CHECK_BUILD)
LIST(FIND VALID_BUILDS "${CHECK_BUILD}" CURRENT_BUILD)

if(CURRENT_BUILD LESS 0)
  message(SEND_ERROR "Please specify a valid CMAKE_BUILD_TYPE type! Valid build types are: " "${VALID_BUILDS}")
endif()


#
# PROJECT is Su2Theta
#
project("Su2Theta")

#######################################################################
# Su2Thta OPTIONS
# likely to be changed by users
#######################################################################

# build options
set(USE_GPU OFF CACHE BOOL "builds nVidia gpu support")
set(USE_AVX OFF CACHE BOOL "builds Avx support")
set(USE_AVX2 OFF CACHE BOOL "builds Avx2 support")
set(USE_AVX512-SKX OFF CACHE BOOL "builds Avx-512 support for Skylake")
set(USE_AVX512-KNL OFF CACHE BOOL "builds Avx-512 support for Knights Landing")
set(USE_FMA OFF CACHE BOOL "builds Fma support")
set(HDF5_PARALLEL_COMPRESSION OFF CACHE BOOL "uses experimental parallel compression, requires Hdf5 > 1.11")

# do all the build definitions
#

find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_CXX_INCLUDE_PATH})
#find_package(HDF5 REQUIRED)
#find_package (HDF5 NAMES hdf5 COMPONENTS C STATIC)
include_directories(SYSTEM ${HDF5_INCLUDE_DIR})

find_package(OpenMP REQUIRED)

find_library(HDF5_LIBRARY NAMES hdf5)
if(${USE_GPU})
  add_definitions(-DUSE_GPU)
  set(GPU_ARCH sm_52 CACHE STRING "set the GPU architecture (sm_30, sm_35, sm_50, sm_52, sm_60, sm_61, sm_62)")
  set(GPU_COMP compute_52 CACHE STRING "set the GPU compute capability (compute_30, compute_35, compute_50, compute_52, compute_60, compute_61, compute_62)")
  set(DEFARCH "520")
  find_package(CUDA REQUIRED)
  #set(CUDA_HOST_COMPILER gcc)
  set(CUDA_PROPAGATE_HOST_FLAGS ON)
  set(CUDA_NVCC_FLAGS "-arch=${GPU_COMP} -code=${GPU_ARCH} -ftz=true -prec-div=false -prec-sqrt=false -std=c++17")

  set(CUDA_NVCC_FLAGS_DEVEL "-O3 -lineinfo" CACHE STRING
      "Flags used by the CUDA compiler during regular development builds." FORCE )
  set(CUDA_NVCC_FLAGS_RELEASE "-O3 -w" CACHE STRING
    "Flags used by the C++ compiler during release builds." FORCE )
  set(CUDA_NVCC_FLAGS_HOSTDEBUG "-g -G" CACHE STRING
      "Flags used by the C++ compiler during debug builds." FORCE )
endif()

#if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
#  # using Visual Studio C++
if(${USE_FMA})
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(FMA_FLAGS "-mfma")
  elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(FMA_FLAGS "-fma")
  endif()
endif()

if(${USE_AVX512-SKX})
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(AVX_FLAGS "-mavx512f -mavx512bw")
  elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(AVX_FLAGS "-xcore-avx512")
  endif()
elseif(${USE_AVX512-KNL})
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(AVX_FLAGS "-mavx512f -mavx512er")
  elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(AVX_FLAGS "-xmic-avx512")
  endif()
elseif(${USE_AVX512F})
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(AVX_FLAGS "-mavx512f")
  elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(AVX_FLAGS "-xcommon-avx512")
  endif()
elseif(${USE_AVX2})
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(AVX_FLAGS "-mavx2")
  elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(AVX_FLAGS "-xcore-avx2")
  endif()
elseif(${USE_AVX})
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(AVX_FLAGS "-mavx")
  elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(AVX_FLAGS "-xavx")
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(SSE_FLAGS "-msse4.1")
endif()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} ${AVX_FLAGS} ${SSE_FLAGS} ${FMA_FLAGS} -mcmodel=medium -std=c++17")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} ${AVX_FLAGS} ${SSE_FLAGS} ${FMA_FLAGS} -mcmodel=medium -std=c++17")
##
##if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
##    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -inline-forceinline")
##    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -inline-forceinline")
##endif()

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wa,-q -lc++")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wa,-q -lc++")
endif()

set(CMAKE_CXX_FLAGS_DEVEL "-O3 -Wall" CACHE STRING
    "Flags used by the C++ compiler during regular development builds." FORCE )
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -w" CACHE STRING
    "Flags used by the C++ compiler during release builds." FORCE )
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -g -fno-inline" CACHE STRING
    "Flags used by the C++ compiler during debug builds." FORCE )


include_directories(include)
include_directories(SYSTEM ${CUDA_INCLUDE_DIRS})

### add tests and library
add_subdirectory(lib)
add_subdirectory(test)
