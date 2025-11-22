# uGMRT CUDA Imager

Raw CUDA implementation for direct visibility-to-image pipeline. No Kokkos or HPG dependencies - just CUDA toolkit and cuFFT.

## Overview

Pipeline:
```
Raw uGMRT data → svfits processing → CUDA gridding → cuFFT → Image
```

## Features

- **Raw CUDA kernels**: Direct GPU programming, no abstraction layers
- **Full svfits processing**: Burst finding, bandpass, baseline sub, RFI flagging
- **cuFFT acceleration**: GPU-accelerated 2D FFT
- **Minimal dependencies**: CUDA toolkit + CFITSIO only
- **User-supplied CF**: Bring your own convolution kernel

## Requirements

- CUDA Toolkit >= 11.0
- CFITSIO
- GCC or compatible C compiler
- (Optional) OpenMP for parallel file reading

## Building

```bash
cd ugmrt_cuda
mkdir build && cd build
cmake ..
make -j$(nproc)
```

## Usage

### 1. Generate convolution function (or use your own)

```bash
./generate_cf -o pswf.cf -s 7 -O 128
```

This creates a prolate spheroidal wave function kernel with:
- Support = 7 pixels (full width = 15)
- Oversampling = 128x

### 2. Run the imager

```bash
./ugmrt_cuda \
    -u svfits_par.txt \
    -a antsamp.hdr \
    -k pswf.cf \
    -o burst_image.fits \
    -n 512,512 \
    -c 1.0
```

### Options

```
Required:
  -u FILE    Parameter file (svfits_par.txt format)
  -a FILE    Antenna/sampler header file
  -k FILE    Convolution kernel file

Optional:
  -o FILE    Output FITS image (default: burst_image.fits)
  -n NX,NY   Grid size (default: 512,512)
  -c CELL    Cell size in arcseconds (default: 1.0)
  -b BATCH   Visibility batch size (default: 100000)
  -T THRESH  Flagging threshold in MAD (default: 5.0)
  --no-bandpass    Disable bandpass correction
  --no-baseline    Disable baseline subtraction
  --no-flag        Disable RFI flagging
```

## Convolution Function Format

Binary file format:
```
int32:  support        (half-width in pixels, e.g., 7)
int32:  oversampling   (e.g., 128)
float[]: values        [(2*support+1)*oversampling]^2 floats
```

Total size for support=7, oversampling=128:
- Full size = (2*7+1) * 128 = 1920
- Values = 1920 * 1920 = 3,686,400 floats = 14.7 MB

### Creating your own CF

```c
// Example: create CF programmatically
int support = 7;
int oversampling = 128;
int full_size = (2 * support + 1) * oversampling;

float* cf = malloc(full_size * full_size * sizeof(float));

// Fill with your kernel values...
for (int y = 0; y < full_size; y++) {
    for (int x = 0; x < full_size; x++) {
        cf[y * full_size + x] = your_kernel_function(x, y);
    }
}

// Write to file
FILE* fp = fopen("my_kernel.cf", "wb");
fwrite(&support, sizeof(int), 1, fp);
fwrite(&oversampling, sizeof(int), 1, fp);
fwrite(cf, sizeof(float), full_size * full_size, fp);
fclose(fp);
```

## Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                    svfits_reader.c                          │
│  - Reads raw uGMRT SPOTLIGHT data                          │
│  - Applies bandpass, baseline, RFI flagging                │
│  - Computes UVW coordinates                                │
│  - Outputs CudaVisibility structs                          │
└─────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────┐
│                    cuda_gridder.cu                          │
│  CUDA Kernels:                                              │
│  - grid_visibility_kernel: Scatter visibilities to grid    │
│  - normalize_kernel: Divide by weights                     │
│  - fftshift_kernel: Center the image                       │
│  - gridding_correction_kernel: Correct for CF convolution  │
│                                                             │
│  Host Functions:                                            │
│  - cf_load(): Load convolution function                    │
│  - grid_create(): Allocate GPU grid                        │
│  - grid_visibilities(): Batch gridding                     │
│  - grid_fft(): cuFFT wrapper                               │
└─────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────┐
│                    fits_output.c                            │
│  - Writes standard FITS image with WCS headers             │
└─────────────────────────────────────────────────────────────┘
```

## Performance

Expected performance for 512x512 grid:
- Gridding: ~1-5 sec (depending on visibility count)
- FFT: < 0.1 sec (cuFFT is fast!)
- Total: ~5-10 sec for typical burst data

## CUDA Architecture

The gridding kernel uses atomic operations for thread-safe accumulation:

```cuda
__global__ void grid_visibility_kernel(...) {
    // Each thread processes one visibility
    int tid = blockIdx.x * blockDim.x + threadIdx.x;

    // Scatter to grid using convolution function
    for (int dv = -support; dv <= support; dv++) {
        for (int du = -support; du <= support; du++) {
            float cf_val = cf[...];
            atomicAdd(&grid[idx].x, vis.re * weight * cf_val);
            atomicAdd(&grid[idx].y, vis.im * weight * cf_val);
        }
    }
}
```

## Files

```
ugmrt_cuda/
├── CMakeLists.txt
├── README.md
├── include/
│   ├── cuda_types.h       # Common types
│   ├── cuda_gridder.h     # Gridding API
│   └── svfits_reader.h    # Reader API
├── src/
│   ├── cuda_gridder.cu    # CUDA kernels
│   ├── svfits_reader.c    # svfits wrapper
│   ├── fits_output.c      # FITS writer
│   └── main.c             # CLI
├── tools/
│   └── generate_cf.c      # CF generator
└── svfits/                # Original svfits C code
```

## License

Apache 2.0

## Authors

uGMRT SPOTLIGHT Team
