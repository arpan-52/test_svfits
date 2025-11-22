cmake_minimum_required(VERSION 3.18)
project(ugmrt_cuda LANGUAGES C CUDA)

# CUDA standard
set(CMAKE_CUDA_STANDARD 11)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)
set(CMAKE_C_STANDARD 99)

# Find CUDA
find_package(CUDAToolkit REQUIRED)

# Find CFITSIO
find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
    pkg_check_modules(CFITSIO cfitsio)
endif()
if(NOT CFITSIO_FOUND)
    find_library(CFITSIO_LIBRARIES NAMES cfitsio)
    find_path(CFITSIO_INCLUDE_DIRS NAMES fitsio.h)
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
    svfits/sla_prenut.c
)

add_library(svfits_core STATIC ${SVFITS_SOURCES})
target_include_directories(svfits_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/svfits
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CFITSIO_INCLUDE_DIRS}
)
target_link_libraries(svfits_core PUBLIC ${CFITSIO_LIBRARIES} m)

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
)

target_include_directories(cuda_gridder PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(cuda_gridder PUBLIC
    CUDA::cudart
    CUDA::cufft
)

# CUDA architecture (adjust for your GPU)
set_target_properties(cuda_gridder PROPERTIES
    CUDA_ARCHITECTURES "60;70;75;80;86"
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
message(STATUS "  OpenMP: ${OpenMP_C_FOUND}")
message(STATUS "========================================")
message(STATUS "")
