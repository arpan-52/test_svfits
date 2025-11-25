/**
 * @file cf_generator.cu
 * @brief W-projection CF generator (GPU implementation)
 *
 * Standard W-projection following Cornwell et al. 2008
 */

#include "cf_generator.h"
#include <cuda_runtime.h>
#include <cufft.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//-----------------------------------------------------------------------------
// Prolate Spheroidal Wave Function (PSWF) for anti-aliasing
//-----------------------------------------------------------------------------

/**
 * @brief PSWF (Prolate Spheroidal Wave Function) for anti-aliasing
 *
 * This is a standard anti-aliasing function used in CASA/libhpg.
 * Approximation using Gaussian taper for simplicity.
 * For exact PSWF, need to solve differential equation (complex).
 */
__device__ __host__ static float pswf(float r, int support) {
    // Gaussian approximation to PSWF
    // α parameter controls width (typical: α = 1.0 for support region)
    float alpha = 1.0f;
    float x = r / support;  // Normalize to [0, 1]

    if (x >= 1.0f) return 0.0f;  // Outside support

    // Gaussian taper: exp(-α * x²)
    return expf(-alpha * x * x);
}

//-----------------------------------------------------------------------------
// W-projection screen computation
//-----------------------------------------------------------------------------

/**
 * @brief CUDA kernel: Compute W-projection screen in image domain
 *
 * For each W-plane, computes:
 *   W(l,m) = exp(2πi * w * (√(1 - l² - m²) - 1)) * PSWF(r)
 *
 * Where (l,m) are direction cosines in image domain.
 */
__global__ void compute_w_screen_kernel(
    cuFloatComplex* screen,   // Output: W-screen [cf_size * cf_size * oversampling²]
    float w,                  // W value for this plane
    int support,              // CF half-width
    int oversampling,         // Oversampling factor
    int cf_size,              // 2*support + 1
    double cell_rad,          // Cell size (radians)
    int grid_size)            // Grid size (for normalization)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    int total_samples = cf_size * cf_size * oversampling * oversampling;
    if (idx >= total_samples) return;

    // Decode 4D index: [x_major, y_major, x_minor, y_minor]
    int y_minor = idx % oversampling;
    int x_minor = (idx / oversampling) % oversampling;
    int y_major = (idx / (oversampling * oversampling)) % cf_size;
    int x_major = idx / (oversampling * oversampling * cf_size);

    // Pixel coordinates (oversampled)
    float px = (x_major - support) + (x_minor - oversampling/2.0f) / (float)oversampling;
    float py = (y_major - support) + (y_minor - oversampling/2.0f) / (float)oversampling;

    // Direction cosines (l,m)
    // l = px * cell_size, m = py * cell_size
    double l = px * cell_rad;
    double m = py * cell_rad;

    // W-projection phase: exp(2πi * w * (√(1 - l² - m²) - 1))
    double lm_sq = l*l + m*m;
    double n = 0.0;

    if (lm_sq < 1.0) {
        n = sqrt(1.0 - lm_sq) - 1.0;
    } else {
        // Outside visible hemisphere - set to zero
        screen[idx].x = 0.0f;
        screen[idx].y = 0.0f;
        return;
    }

    double phase = 2.0 * M_PI * w * n;

    // Complex exponential
    cuFloatComplex w_phase;
    w_phase.x = cosf(phase);
    w_phase.y = sinf(phase);

    // Anti-aliasing taper (PSWF)
    float r = sqrtf(px*px + py*py);
    float taper = pswf(r, support);

    // Combined: W-screen * taper
    screen[idx].x = w_phase.x * taper;
    screen[idx].y = w_phase.y * taper;
}

//-----------------------------------------------------------------------------
// Copy screen to CF array (helper kernel)
//-----------------------------------------------------------------------------

/**
 * @brief Copy FFT'd screen to proper location in 6D CF array
 */
__global__ void copy_screen_to_cf(
    cuFloatComplex* cf,        // Output CF array (6D)
    cuFloatComplex* screen,    // Input screen (2D oversampled)
    int w_plane,               // W-plane index
    int cf_size,               // 2*support + 1
    int oversampling,          // Oversampling factor
    int fft_size,              // cf_size * oversampling
    int stride_x_major, int stride_y_major,
    int stride_mueller, int stride_cube,
    int stride_x_minor, int stride_y_minor,
    float norm)                // FFT normalization
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int total_samples = cf_size * cf_size * oversampling * oversampling;

    if (idx >= total_samples) return;

    // Decode 4D index
    int y_minor = idx % oversampling;
    int x_minor = (idx / oversampling) % oversampling;
    int y_major = (idx / (oversampling * oversampling)) % cf_size;
    int x_major = idx / (oversampling * oversampling * cf_size);

    // Screen index (2D oversampled grid)
    int sy = y_major * oversampling + y_minor;
    int sx = x_major * oversampling + x_minor;
    int screen_idx = sy * fft_size + sx;

    // CF 6D index
    int cf_idx = x_major * stride_x_major +
                 y_major * stride_y_major +
                 0 * stride_mueller +
                 w_plane * stride_cube +
                 x_minor * stride_x_minor +
                 y_minor * stride_y_minor;

    // Copy and normalize
    cf[cf_idx].x = screen[screen_idx].x * norm;
    cf[cf_idx].y = screen[screen_idx].y * norm;
}

//-----------------------------------------------------------------------------
// CF generation main function
//-----------------------------------------------------------------------------

HPGConvolutionFunction* cf_generate_w_projection(const CFGeneratorConfig* config) {
    printf("Generating W-projection CFs...\n");
    printf("  W-planes: %d\n", config->n_w_planes);
    printf("  Support: %d pixels\n", config->support);
    printf("  Oversampling: %dx\n", config->oversampling);
    printf("  Max W: %.2f wavelengths\n", config->max_w);

    // Allocate CF structure
    HPGConvolutionFunction* cf = (HPGConvolutionFunction*)malloc(sizeof(HPGConvolutionFunction));
    if (!cf) {
        fprintf(stderr, "Failed to allocate CF structure\n");
        return NULL;
    }

    cf->support = config->support;
    cf->oversampling = config->oversampling;
    cf->padding = config->padding;
    cf->n_mueller = config->n_mueller;
    cf->n_cube = config->n_w_planes;
    cf->cf_size = 2 * config->support + 1;

    // Calculate 6D strides (row-major, fastest varying last)
    cf->stride_y_minor = 1;
    cf->stride_x_minor = config->oversampling;
    cf->stride_cube = cf->stride_x_minor * config->oversampling;
    cf->stride_mueller = cf->stride_cube * config->n_w_planes;
    cf->stride_y_major = cf->stride_mueller * config->n_mueller;
    cf->stride_x_major = cf->stride_y_major * cf->cf_size;

    size_t n_values = cf->cf_size * cf->cf_size * config->n_mueller * config->n_w_planes *
                      config->oversampling * config->oversampling;
    cf->n_values = n_values;

    printf("  Total CF samples: %zu (%.2f MB)\n",
           n_values, n_values * sizeof(cuFloatComplex) / (1024.0 * 1024.0));

    // Allocate CF array on GPU
    CUDA_CHECK(cudaMalloc(&cf->d_values, n_values * sizeof(cuFloatComplex)));
    CUDA_CHECK(cudaMemset(cf->d_values, 0, n_values * sizeof(cuFloatComplex)));

    // Allocate host memory for one W-screen at a time
    int screen_size = cf->cf_size * cf->cf_size * config->oversampling * config->oversampling;
    cuFloatComplex* d_screen;
    CUDA_CHECK(cudaMalloc(&d_screen, screen_size * sizeof(cuFloatComplex)));

    // Generate each W-plane
    float w_min = -config->max_w;
    float w_max = config->max_w;

    for (int w_plane = 0; w_plane < config->n_w_planes; w_plane++) {
        // Calculate W value for this plane
        float w;
        if (config->n_w_planes == 1) {
            w = 0.0f;  // Single plane at W=0
        } else {
            w = w_min + (w_max - w_min) * w_plane / (config->n_w_planes - 1);
        }

        printf("  W-plane %3d/%d: w = %+10.2f wavelengths\n",
               w_plane + 1, config->n_w_planes, w);

        // Compute W-screen in image domain
        int block_size = 256;
        int n_blocks = (screen_size + block_size - 1) / block_size;

        compute_w_screen_kernel<<<n_blocks, block_size>>>(
            d_screen, w,
            config->support, config->oversampling, cf->cf_size,
            config->cell_rad, config->grid_nx
        );
        CUDA_CHECK(cudaDeviceSynchronize());

        // FFT: Image domain → UV domain
        // This gives us the convolution kernel
        cufftHandle plan;
        int fft_size = cf->cf_size * config->oversampling;
        CUFFT_CHECK(cufftPlan2d(&plan, fft_size, fft_size, CUFFT_C2C));
        CUFFT_CHECK(cufftExecC2C(plan, d_screen, d_screen, CUFFT_FORWARD));
        cufftDestroy(plan);

        // Copy and normalize screen to CF array on GPU
        // Use helper kernel below
        int copy_blocks = (screen_size + 255) / 256;
        float norm = 1.0f / (fft_size * fft_size);

        copy_screen_to_cf<<<copy_blocks, 256>>>(
            cf->d_values, d_screen,
            w_plane, cf->cf_size, config->oversampling, fft_size,
            cf->stride_x_major, cf->stride_y_major,
            cf->stride_mueller, cf->stride_cube,
            cf->stride_x_minor, cf->stride_y_minor,
            norm
        );
        CUDA_CHECK(cudaDeviceSynchronize());
    }

    cudaFree(d_screen);

    printf("CF generation complete!\n");
    return cf;
}

//-----------------------------------------------------------------------------
// W-plane selection
//-----------------------------------------------------------------------------

int cf_select_w_plane(float w, float max_w, int n_w_planes) {
    if (n_w_planes == 1) return 0;

    // Map W ∈ [-max_w, +max_w] to plane ∈ [0, n_w_planes-1]
    float w_norm = (w + max_w) / (2.0f * max_w);
    w_norm = fmaxf(0.0f, fminf(1.0f, w_norm));

    int plane = (int)(w_norm * (n_w_planes - 1) + 0.5f);
    return plane;
}

//-----------------------------------------------------------------------------
// Max W calculation
//-----------------------------------------------------------------------------

float cf_compute_max_w(const double* antenna_x, const double* antenna_y, const double* antenna_z,
                       int n_antennas, double dec_rad, double wavelength) {
    double max_baseline = 0.0;

    // Find maximum baseline length
    for (int i = 0; i < n_antennas; i++) {
        for (int j = i + 1; j < n_antennas; j++) {
            double dx = antenna_x[i] - antenna_x[j];
            double dy = antenna_y[i] - antenna_y[j];
            double dz = antenna_z[i] - antenna_z[j];
            double baseline = sqrt(dx*dx + dy*dy + dz*dz);

            if (baseline > max_baseline) {
                max_baseline = baseline;
            }
        }
    }

    // W component depends on declination
    // max_w ≈ max_baseline * sin(dec) / wavelength
    double max_w = max_baseline * fabs(sin(dec_rad)) / wavelength;

    return (float)max_w;
}
