cmake_minimum_required(VERSION 3.18)
project(ugmrt_imager LANGUAGES C CXX)

# C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Options
option(USE_CUDA "Enable CUDA support" ON)
option(USE_OPENMP "Enable OpenMP support" ON)
option(USE_NOVAS "Use NOVAS instead of SLA for astrometry" OFF)

# Find required packages
find_package(Kokkos QUIET)
find_package(hpg CONFIG QUIET)

# Check for HPG
if(NOT hpg_FOUND)
    message(STATUS "HPG not found via find_package, looking for local build...")
    set(HPG_DIR "$ENV{HPG_DIR}" CACHE PATH "Path to HPG installation")
    if(EXISTS "${HPG_DIR}/lib/cmake/hpg")
        list(APPEND CMAKE_PREFIX_PATH "${HPG_DIR}")
        find_package(hpg CONFIG REQUIRED)
    else()
        message(FATAL_ERROR "HPG not found. Set HPG_DIR or install HPG.")
    endif()
endif()

# Find CFITSIO
find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
    pkg_check_modules(CFITSIO cfitsio)
endif()
if(NOT CFITSIO_FOUND)
    find_library(CFITSIO_LIBRARIES NAMES cfitsio)
    find_path(CFITSIO_INCLUDE_DIRS NAMES fitsio.h)
endif()

# Find FFTW (for CPU fallback)
if(PkgConfig_FOUND)
    pkg_check_modules(FFTW fftw3)
endif()

# OpenMP
if(USE_OPENMP)
    find_package(OpenMP)
endif()

# Compiler flags
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra -O3)
endif()

#-----------------------------------------------------------------------------
# svfits core library (C code from original svfits)
#-----------------------------------------------------------------------------
set(SVFITS_SOURCES
    svfits/svsubs.c
    svfits/utils.c
    svfits/stats.c
)

if(USE_NOVAS)
    list(APPEND SVFITS_SOURCES svfits/novas_prenut.c)
    add_definitions(-DUSE_NOVAS)
else()
    list(APPEND SVFITS_SOURCES svfits/sla_prenut.c)
endif()

add_library(svfits_core STATIC ${SVFITS_SOURCES})
target_include_directories(svfits_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/svfits
    ${CFITSIO_INCLUDE_DIRS}
)
target_link_libraries(svfits_core PUBLIC ${CFITSIO_LIBRARIES})

if(OpenMP_C_FOUND)
    target_link_libraries(svfits_core PUBLIC OpenMP::OpenMP_C)
endif()

# Suppress C warnings for legacy code
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(svfits_core PRIVATE -w)
endif()

#-----------------------------------------------------------------------------
# ugmrt_imager library
#-----------------------------------------------------------------------------
add_library(ugmrt_imager_lib STATIC
    src/svfits_reader.cpp
    src/vis_to_hpg.cpp
    src/cf_kernel.cpp
    src/fits_output.cpp
    src/ugmrt_imager.cpp
)

target_include_directories(ugmrt_imager_lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/svfits
    ${CFITSIO_INCLUDE_DIRS}
)

target_link_libraries(ugmrt_imager_lib PUBLIC
    svfits_core
    Hpg::hpg
    ${CFITSIO_LIBRARIES}
)

if(OpenMP_CXX_FOUND)
    target_link_libraries(ugmrt_imager_lib PUBLIC OpenMP::OpenMP_CXX)
endif()

#-----------------------------------------------------------------------------
# Main executable
#-----------------------------------------------------------------------------
add_executable(ugmrt_imager src/main.cpp)
target_link_libraries(ugmrt_imager PRIVATE ugmrt_imager_lib)

#-----------------------------------------------------------------------------
# Installation
#-----------------------------------------------------------------------------
install(TARGETS ugmrt_imager ugmrt_imager_lib svfits_core
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

install(DIRECTORY include/ DESTINATION include)

#-----------------------------------------------------------------------------
# Print configuration summary
#-----------------------------------------------------------------------------
message(STATUS "")
message(STATUS "========================================")
message(STATUS "uGMRT Fast Imager Configuration")
message(STATUS "========================================")
message(STATUS "  C++ Compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "  C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "  HPG: ${hpg_DIR}")
message(STATUS "  Kokkos: ${Kokkos_DIR}")
message(STATUS "  CFITSIO: ${CFITSIO_LIBRARIES}")
message(STATUS "  OpenMP: ${OpenMP_CXX_FOUND}")
message(STATUS "  CUDA: ${USE_CUDA}")
message(STATUS "  Astrometry: ${USE_NOVAS}")
message(STATUS "========================================")
message(STATUS "")
