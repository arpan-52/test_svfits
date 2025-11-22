# uGMRT Fast Imager

Direct visibility-to-image pipeline for uGMRT SPOTLIGHT burst data using GPU-accelerated gridding.

## Overview

This tool replaces the traditional multi-step pipeline:
```
Raw → svfits → FITS → CASA → MS → WSClean → Image
```

With a streamlined direct pipeline:
```
Raw → [Read + Grid in parallel] → FFT → Image
```

**Expected speedup: 5-20x** by eliminating intermediate file I/O.

## Features

- **Direct GPU gridding** using HPG (Hyperion Polyphase Gridder)
- **Full svfits processing**: burst finding, bandpass, baseline subtraction, RFI flagging
- **Zero intermediate files**: data flows directly from raw to image
- **Parallel file reading**: OpenMP-accelerated raw data reading
- **CUDA/cuFFT acceleration**: GPU-accelerated gridding and FFT
- **Standard FITS output**: Compatible with any FITS viewer

## Architecture

```
┌─────────────────────────────────────────────────────────────────────┐
│                     svfits Reader (C)                               │
│  read_slice() → half_to_float() → bandpass → baseline → UVW       │
└─────────────────────────────────────────────────────────────────────┘
                              ↓
                    VisData conversion (C++)
                              ↓
┌─────────────────────────────────────────────────────────────────────┐
│                     HPG Gridder (GPU)                               │
│  grid_visibilities() → normalize_by_weights() → apply_grid_fft()  │
└─────────────────────────────────────────────────────────────────────┘
                              ↓
                         FITS Image
```

## Requirements

### Build Dependencies
- CMake >= 3.18
- C++17 compiler (GCC >= 9 or Clang >= 10)
- Kokkos >= 4.0
- HPG (Hyperion Polyphase Gridder)
- CFITSIO
- FFTW3 (for CPU fallback)
- CUDA >= 11.0 (optional, for GPU acceleration)

### Runtime
- GPU with CUDA support (recommended)
- Or multi-core CPU with OpenMP

## Building

```bash
# Clone the repository
cd pico/ugmrt_imager

# Create build directory
mkdir build && cd build

# Configure with CMake
cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DHPG_DIR=/path/to/hpg \
    -DKokkos_DIR=/path/to/kokkos

# Build
make -j$(nproc)

# Install (optional)
make install
```

## Usage

```bash
ugmrt_imager -u svfits_par.txt -a antsamp.hdr [options]

Required:
  -u FILE    Parameter file (svfits_par.txt format)
  -a FILE    Antenna/sampler header file

Optional:
  -B FILE    Bulletin A file for DUT1 corrections
  -o FILE    Output FITS image (default: burst_image.fits)
  -n NX,NY   Grid size (default: 512,512)
  -c CELL    Cell size in arcseconds (default: 1.0)
  -s SUPPORT CF support in pixels (default: 7)
  -d DEVICE  Device: cuda, openmp, serial (default: cuda)
  -t THREADS CPU threads for reading (default: 4)
  -T THRESH  Flagging threshold in MAD (default: 5.0)
  --no-bandpass    Disable bandpass correction
  --no-baseline    Disable baseline subtraction
  --no-flag        Disable RFI flagging
```

### Example

```bash
# Process burst data with GPU gridding
ugmrt_imager \
    -u /data/bursts/svfits_par.txt \
    -a /data/config/antsamp.hdr \
    -o burst_FRB20210101.fits \
    -n 1024 \
    -c 0.5 \
    -d cuda
```

## Parameter File Format

The parameter file follows the svfits_par.txt format:

```
# Burst parameters
BURST_NAME    FRB20210101
BURST_MJD     59245.123456
BURST_DM      500.0
BURST_INTWD   0.001
BURST_FREQ    650e6

# Observing setup
FREQ_SET      550e6:750e6:4096
RA_APP        1.234
DEC_APP       0.567

# Processing options
DO_BAND       1
DO_BASE       1
DO_FLAG       1
THRESH        5.0

# File paths
NFILE         16
PATH          /data/raw/
INPUT         file00.dat,file01.dat,...
```

## Performance

| Stage | Traditional Pipeline | uGMRT Fast Imager |
|-------|---------------------|-------------------|
| Read raw | 1x | 1x |
| Process vis | 1x | 1x |
| Write FITS | 5-10 sec | **0** |
| CASA import | 30-60 sec | **0** |
| WSClean | 60-120 sec | **1-5 sec (GPU)** |
| **Total** | ~3 min | ~10-30 sec |

## Code Structure

```
ugmrt_imager/
├── CMakeLists.txt          # Build configuration
├── README.md               # This file
├── include/
│   ├── ugmrt_imager.hpp    # Main API
│   ├── svfits_reader.hpp   # svfits wrapper
│   ├── vis_to_hpg.hpp      # Visibility conversion
│   └── cf_kernel.hpp       # Convolution function
├── src/
│   ├── main.cpp            # CLI entry point
│   ├── ugmrt_imager.cpp    # Main implementation
│   ├── svfits_reader.cpp   # svfits C wrapper
│   ├── vis_to_hpg.cpp      # HPG data conversion
│   ├── cf_kernel.cpp       # PSWF kernel
│   └── fits_output.cpp     # FITS writer
└── svfits/                 # Original svfits C code
    ├── svsubs.c
    ├── utils.c
    ├── stats.c
    ├── svio.h
    └── ...
```

## License

Apache 2.0 (HPG components)
See individual files for other licenses.

## Authors

- uGMRT SPOTLIGHT Team
- HPG: NRAO (Martin Pokorny)

## References

- HPG: https://gitlab.nrao.edu/mpokorny/hpg
- Kokkos: https://github.com/kokkos/kokkos
- CFITSIO: https://heasarc.gsfc.nasa.gov/fitsio/
