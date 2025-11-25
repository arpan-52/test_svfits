cmake_minimum_required(VERSION 3.18)
project(ugmrt_cuda LANGUAGES C CUDA)

# CUDA standard
set(CMAKE_CUDA_STANDARD 11)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)
set(CMAKE_C_STANDARD 99)

# Options
option(USE_NOVAS "Use NOVAS library for coordinate transforms" ON)
set(NOVAS_PATH "" CACHE PATH "Path to NOVAS library (containing libnovasc.a, novas.h)")

# Check NOVAS_PATH is set when USE_NOVAS is ON
if(USE_NOVAS AND NOT NOVAS_PATH)
    message(FATAL_ERROR
        "USE_NOVAS is ON but NOVAS_PATH is not set.\n"
        "Please specify the path to NOVAS library:\n"
        "  cmake .. -DNOVAS_PATH=/path/to/novas\n"
        "Or disable NOVAS:\n"
        "  cmake .. -DUSE_NOVAS=OFF"
    )
endif()

# Find CUDA
find_package(CUDAToolkit REQUIRED)

# Find CFITSIO
find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
    pkg_check_modules(CFITSIO cfitsio)
endif()
if(NOT CFITSIO_FOUND)
    # Try hardcoded path first (like original Makefile)
    if(EXISTS "${NOVAS_PATH}/libcfitsio.a")
        set(CFITSIO_LIBRARIES "${NOVAS_PATH}/libcfitsio.a")
        set(CFITSIO_INCLUDE_DIRS "${NOVAS_PATH}")
    else()
        find_library(CFITSIO_LIBRARIES NAMES cfitsio)
        find_path(CFITSIO_INCLUDE_DIRS NAMES fitsio.h)
    endif()
endif()

# NOVAS library
if(USE_NOVAS)
    if(EXISTS "${NOVAS_PATH}/libnovasc.a")
        set(NOVAS_LIBRARY "${NOVAS_PATH}/libnovasc.a")
        set(NOVAS_INCLUDE_DIR "${NOVAS_PATH}")
        message(STATUS "Found NOVAS: ${NOVAS_LIBRARY}")
    else()
        find_library(NOVAS_LIBRARY NAMES novasc novas)
        find_path(NOVAS_INCLUDE_DIR NAMES novas.h)
    endif()

    if(NOVAS_LIBRARY)
        add_compile_definitions(USE_NOVAS)
    else()
        message(WARNING "NOVAS not found at ${NOVAS_PATH}, disabling USE_NOVAS")
        set(USE_NOVAS OFF)
    endif()
endif()

# OpenMP (optional)
find_package(OpenMP)

# Compiler flags
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -O3)
endif()

#-----------------------------------------------------------------------------
# svfits core library (C code)
#-----------------------------------------------------------------------------
set(SVFITS_SOURCES
    svfits/svsubs.c
    svfits/utils.c
    svfits/stats.c
)

# Add appropriate prenut source based on USE_NOVAS
if(USE_NOVAS)
    # Check if novas_prenut.c exists, otherwise use sla_prenut.c
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/svfits/novas_prenut.c")
        list(APPEND SVFITS_SOURCES svfits/novas_prenut.c)
    endif()
else()
    list(APPEND SVFITS_SOURCES svfits/sla_prenut.c)
endif()

add_library(svfits_core STATIC ${SVFITS_SOURCES})
target_include_directories(svfits_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/svfits
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CFITSIO_INCLUDE_DIRS}
)

if(USE_NOVAS AND NOVAS_INCLUDE_DIR)
    target_include_directories(svfits_core PUBLIC ${NOVAS_INCLUDE_DIR})
endif()

target_link_libraries(svfits_core PUBLIC ${CFITSIO_LIBRARIES} m)

if(USE_NOVAS AND NOVAS_LIBRARY)
    target_link_libraries(svfits_core PUBLIC ${NOVAS_LIBRARY})
endif()

if(OpenMP_C_FOUND)
    target_link_libraries(svfits_core PUBLIC OpenMP::OpenMP_C)
endif()

# Suppress warnings for legacy C code
target_compile_options(svfits_core PRIVATE -w)

#-----------------------------------------------------------------------------
# CUDA gridder library
#-----------------------------------------------------------------------------
add_library(cuda_gridder STATIC
    src/cuda_gridder.cu
    src/cf_generator.cu
)

target_include_directories(cuda_gridder PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(cuda_gridder PUBLIC
    CUDA::cudart
    CUDA::cufft
)

# CUDA architecture (adjust for your GPU)
# CUDA 12+ dropped compute_60/70, CUDA 13+ dropped more
# Default to modern architectures
set_target_properties(cuda_gridder PROPERTIES
    CUDA_ARCHITECTURES "80;86;89;90"
)

#-----------------------------------------------------------------------------
# svfits reader
#-----------------------------------------------------------------------------
add_library(svfits_reader STATIC
    src/svfits_reader.c
)

target_include_directories(svfits_reader PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/svfits
)

target_link_libraries(svfits_reader PUBLIC
    svfits_core
)

#-----------------------------------------------------------------------------
# FITS output
#-----------------------------------------------------------------------------
add_library(fits_output STATIC
    src/fits_output.c
)

target_include_directories(fits_output PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CFITSIO_INCLUDE_DIRS}
)

target_link_libraries(fits_output PUBLIC
    ${CFITSIO_LIBRARIES}
)

#-----------------------------------------------------------------------------
# Main executable
#-----------------------------------------------------------------------------
add_executable(ugmrt_cuda src/main.c)

target_include_directories(ugmrt_cuda PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(ugmrt_cuda PRIVATE
    cuda_gridder
    svfits_reader
    fits_output
    m
)

#-----------------------------------------------------------------------------
# CF generator utility (optional)
#-----------------------------------------------------------------------------
add_executable(generate_cf tools/generate_cf.c)
target_include_directories(generate_cf PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(generate_cf m)

#-----------------------------------------------------------------------------
# Installation
#-----------------------------------------------------------------------------
install(TARGETS ugmrt_cuda generate_cf
    RUNTIME DESTINATION bin
)

#-----------------------------------------------------------------------------
# Print configuration
#-----------------------------------------------------------------------------
message(STATUS "")
message(STATUS "========================================")
message(STATUS "uGMRT CUDA Imager Configuration")
message(STATUS "========================================")
message(STATUS "  CUDA: ${CUDAToolkit_VERSION}")
message(STATUS "  CFITSIO: ${CFITSIO_LIBRARIES}")
message(STATUS "  NOVAS: ${USE_NOVAS}")
if(USE_NOVAS)
message(STATUS "  NOVAS_PATH: ${NOVAS_PATH}")
message(STATUS "  NOVAS_LIB: ${NOVAS_LIBRARY}")
endif()
message(STATUS "  OpenMP: ${OpenMP_C_FOUND}")
message(STATUS "========================================")
message(STATUS "")
