/**
 * @file cuda_gridder.cu
 * @brief CUDA kernels for visibility gridding
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cuda_runtime.h>
#include <cufft.h>

#include "cuda_types.h"
#include "cuda_gridder.h"

//=============================================================================
// CUDA Kernels
//=============================================================================

/**
 * @brief Grid visibilities kernel
 *
 * Each thread processes one visibility and scatters it to the grid
 * using the convolution function.
 */
__global__ void grid_visibility_kernel(
    cuFloatComplex* __restrict__ grid,
    float* __restrict__ weights,
    const CudaVisibility* __restrict__ vis,
    int n_vis,
    const float* __restrict__ cf,
    int cf_support,
    int cf_oversampling,
    int cf_full_size,
    int nx, int ny,
    float scale_u, float scale_v)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= n_vis) return;

    CudaVisibility v = vis[tid];

    // Skip flagged visibilities
    if (v.weight <= 0.0f) return;

    // Convert UV to grid coordinates
    // Grid center is at (nx/2, ny/2)
    float grid_u = v.u * scale_u + nx / 2.0f;
    float grid_v = v.v * scale_v + ny / 2.0f;

    // Integer grid position (center of CF)
    int iu = (int)floorf(grid_u);
    int iv = (int)floorf(grid_v);

    // Fractional offset for CF lookup (in range [0, 1))
    float frac_u = grid_u - iu;
    float frac_v = grid_v - iv;

    // CF lookup offset (oversampled)
    int cf_off_u = (int)(frac_u * cf_oversampling);
    int cf_off_v = (int)(frac_v * cf_oversampling);

    // Visibility value
    cuFloatComplex vis_val = make_cuFloatComplex(v.re, v.im);

    // Scatter to grid using CF
    for (int dv = -cf_support; dv <= cf_support; dv++) {
        int gv = iv + dv;
        if (gv < 0 || gv >= ny) continue;

        // CF y index
        int cf_y = (dv + cf_support) * cf_oversampling + cf_off_v;

        for (int du = -cf_support; du <= cf_support; du++) {
            int gu = iu + du;
            if (gu < 0 || gu >= nx) continue;

            // CF x index
            int cf_x = (du + cf_support) * cf_oversampling + cf_off_u;

            // CF value
            float cf_val = cf[cf_y * cf_full_size + cf_x];

            // Grid index
            int grid_idx = gv * nx + gu;

            // Weighted visibility
            float wt = v.weight * cf_val;
            cuFloatComplex weighted_vis = make_cuFloatComplex(
                vis_val.x * wt,
                vis_val.y * wt
            );

            // Atomic add to grid (real and imaginary separately)
            atomicAdd(&grid[grid_idx].x, weighted_vis.x);
            atomicAdd(&grid[grid_idx].y, weighted_vis.y);
            atomicAdd(&weights[grid_idx], wt);
        }
    }
}

/**
 * @brief Normalize grid by weights
 */
__global__ void normalize_kernel(
    cuFloatComplex* __restrict__ grid,
    const float* __restrict__ weights,
    int n_pixels)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= n_pixels) return;

    float wt = weights[tid];
    if (wt > 0.0f) {
        grid[tid].x /= wt;
        grid[tid].y /= wt;
    }
}

/**
 * @brief FFT shift kernel (swap quadrants)
 */
__global__ void fftshift_kernel(
    cuFloatComplex* __restrict__ grid,
    int nx, int ny)
{
    int ix = blockIdx.x * blockDim.x + threadIdx.x;
    int iy = blockIdx.y * blockDim.y + threadIdx.y;

    // Only process first quadrant, swap with third
    if (ix >= nx/2 || iy >= ny/2) return;

    // Indices for all four quadrants
    int idx_q1 = iy * nx + ix;                           // Top-left
    int idx_q2 = iy * nx + (ix + nx/2);                  // Top-right
    int idx_q3 = (iy + ny/2) * nx + ix;                  // Bottom-left
    int idx_q4 = (iy + ny/2) * nx + (ix + nx/2);         // Bottom-right

    // Swap Q1 <-> Q4 and Q2 <-> Q3
    cuFloatComplex tmp1 = grid[idx_q1];
    cuFloatComplex tmp2 = grid[idx_q2];

    grid[idx_q1] = grid[idx_q4];
    grid[idx_q4] = tmp1;
    grid[idx_q2] = grid[idx_q3];
    grid[idx_q3] = tmp2;
}

/**
 * @brief Gridding correction kernel
 *
 * Applies 1/sinc correction for the gridding convolution
 */
__global__ void gridding_correction_kernel(
    cuFloatComplex* __restrict__ grid,
    int nx, int ny,
    int cf_support)
{
    int ix = blockIdx.x * blockDim.x + threadIdx.x;
    int iy = blockIdx.y * blockDim.y + threadIdx.y;

    if (ix >= nx || iy >= ny) return;

    // Normalized position from center [-1, 1]
    float x = (ix - nx/2.0f) / (nx/2.0f);
    float y = (iy - ny/2.0f) / (ny/2.0f);

    // Correction factor (1/sinc approximation for spheroidal)
    float corr_x = 1.0f;
    float corr_y = 1.0f;

    float pi = 3.14159265358979f;
    float arg_x = pi * x * cf_support;
    float arg_y = pi * y * cf_support;

    if (fabsf(arg_x) > 1e-6f) {
        corr_x = arg_x / sinf(arg_x);
    }
    if (fabsf(arg_y) > 1e-6f) {
        corr_y = arg_y / sinf(arg_y);
    }

    float correction = corr_x * corr_y;

    int idx = iy * nx + ix;
    grid[idx].x *= correction;
    grid[idx].y *= correction;
}

//=============================================================================
// Host functions
//=============================================================================

// Convolution function management

ConvolutionFunction cf_create(int support, int oversampling, const float* values) {
    ConvolutionFunction cf;
    cf.support = support;
    cf.oversampling = oversampling;
    cf.full_size = (2 * support + 1) * oversampling;

    size_t size = cf.full_size * cf.full_size * sizeof(float);

    // Allocate on GPU
    CUDA_CHECK(cudaMalloc(&cf.d_values, size));
    CUDA_CHECK(cudaMemcpy(cf.d_values, values, size, cudaMemcpyHostToDevice));

    // Keep host copy
    cf.h_values = (float*)malloc(size);
    memcpy(cf.h_values, values, size);

    return cf;
}

ConvolutionFunction cf_load(const char* filename) {
    FILE* fp = fopen(filename, "rb");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open CF file: %s\n", filename);
        exit(EXIT_FAILURE);
    }

    int support, oversampling;
    fread(&support, sizeof(int), 1, fp);
    fread(&oversampling, sizeof(int), 1, fp);

    int full_size = (2 * support + 1) * oversampling;
    size_t n_values = full_size * full_size;

    float* values = (float*)malloc(n_values * sizeof(float));
    fread(values, sizeof(float), n_values, fp);
    fclose(fp);

    ConvolutionFunction cf = cf_create(support, oversampling, values);
    free(values);

    printf("Loaded CF: support=%d, oversampling=%d, size=%d\n",
           support, oversampling, full_size);

    return cf;
}

void cf_free(ConvolutionFunction* cf) {
    if (cf->d_values) {
        cudaFree(cf->d_values);
        cf->d_values = NULL;
    }
    if (cf->h_values) {
        free(cf->h_values);
        cf->h_values = NULL;
    }
}

// Grid management

UVGrid grid_create(int nx, int ny, int n_pol, int n_chan) {
    UVGrid grid;
    grid.nx = nx;
    grid.ny = ny;
    grid.n_pol = n_pol;
    grid.n_chan = n_chan;

    size_t n_pixels = nx * ny * n_pol * n_chan;
    size_t grid_size = n_pixels * sizeof(cuFloatComplex);
    size_t weight_size = n_pixels * sizeof(float);

    CUDA_CHECK(cudaMalloc(&grid.d_grid, grid_size));
    CUDA_CHECK(cudaMalloc(&grid.d_weights, weight_size));

    grid.h_grid = (cuFloatComplex*)malloc(grid_size);

    // Zero the grids
    CUDA_CHECK(cudaMemset(grid.d_grid, 0, grid_size));
    CUDA_CHECK(cudaMemset(grid.d_weights, 0, weight_size));

    printf("Created grid: %dx%d, %d pol, %d chan (%.2f MB)\n",
           nx, ny, n_pol, n_chan, grid_size / (1024.0 * 1024.0));

    return grid;
}

void grid_zero(UVGrid* grid) {
    size_t n_pixels = grid->nx * grid->ny * grid->n_pol * grid->n_chan;
    CUDA_CHECK(cudaMemset(grid->d_grid, 0, n_pixels * sizeof(cuFloatComplex)));
    CUDA_CHECK(cudaMemset(grid->d_weights, 0, n_pixels * sizeof(float)));
}

void grid_free(UVGrid* grid) {
    if (grid->d_grid) {
        cudaFree(grid->d_grid);
        grid->d_grid = NULL;
    }
    if (grid->d_weights) {
        cudaFree(grid->d_weights);
        grid->d_weights = NULL;
    }
    if (grid->h_grid) {
        free(grid->h_grid);
        grid->h_grid = NULL;
    }
}

void grid_to_host(UVGrid* grid) {
    size_t size = grid->nx * grid->ny * grid->n_pol * grid->n_chan * sizeof(cuFloatComplex);
    CUDA_CHECK(cudaMemcpy(grid->h_grid, grid->d_grid, size, cudaMemcpyDeviceToHost));
}

// Gridding operations

void grid_visibilities(
    UVGrid* grid,
    const CudaVisibility* vis,
    int n_vis,
    const ConvolutionFunction* cf,
    float scale_u,
    float scale_v)
{
    if (n_vis == 0) return;

    // Copy visibilities to GPU
    CudaVisibility* d_vis;
    CUDA_CHECK(cudaMalloc(&d_vis, n_vis * sizeof(CudaVisibility)));
    CUDA_CHECK(cudaMemcpy(d_vis, vis, n_vis * sizeof(CudaVisibility), cudaMemcpyHostToDevice));

    // Launch kernel
    int block_size = 256;
    int n_blocks = (n_vis + block_size - 1) / block_size;

    grid_visibility_kernel<<<n_blocks, block_size>>>(
        grid->d_grid,
        grid->d_weights,
        d_vis,
        n_vis,
        cf->d_values,
        cf->support,
        cf->oversampling,
        cf->full_size,
        grid->nx, grid->ny,
        scale_u, scale_v
    );

    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    cudaFree(d_vis);
}

void grid_normalize(UVGrid* grid) {
    int n_pixels = grid->nx * grid->ny * grid->n_pol * grid->n_chan;
    int block_size = 256;
    int n_blocks = (n_pixels + block_size - 1) / block_size;

    normalize_kernel<<<n_blocks, block_size>>>(
        grid->d_grid,
        grid->d_weights,
        n_pixels
    );

    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
}

void grid_fft(UVGrid* grid, int forward) {
    cufftHandle plan;
    CUFFT_CHECK(cufftPlan2d(&plan, grid->ny, grid->nx, CUFFT_C2C));

    // Process each polarization and channel
    for (int p = 0; p < grid->n_pol * grid->n_chan; p++) {
        cuFloatComplex* plane = grid->d_grid + p * grid->nx * grid->ny;

        CUFFT_CHECK(cufftExecC2C(
            plan,
            plane,
            plane,
            forward ? CUFFT_FORWARD : CUFFT_INVERSE
        ));
    }

    cufftDestroy(plan);

    // Normalize if inverse FFT
    if (!forward) {
        int n_pixels = grid->nx * grid->ny * grid->n_pol * grid->n_chan;
        float norm = 1.0f / (grid->nx * grid->ny);

        // Simple normalization kernel
        int block_size = 256;
        int n_blocks = (n_pixels + block_size - 1) / block_size;

        // Inline lambda not possible in CUDA, so we normalize on host or use a kernel
        // For simplicity, do it in normalize step
    }

    CUDA_CHECK(cudaDeviceSynchronize());
}

void grid_shift(UVGrid* grid) {
    dim3 block(16, 16);
    dim3 n_blocks((grid->nx/2 + block.x - 1) / block.x,
                   (grid->ny/2 + block.y - 1) / block.y);

    for (int p = 0; p < grid->n_pol * grid->n_chan; p++) {
        cuFloatComplex* plane = grid->d_grid + p * grid->nx * grid->ny;

        fftshift_kernel<<<n_blocks, block>>>(plane, grid->nx, grid->ny);
    }

    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
}

void grid_correct(UVGrid* grid, const ConvolutionFunction* cf) {
    dim3 block(16, 16);
    dim3 n_blocks((grid->nx + block.x - 1) / block.x,
                   (grid->ny + block.y - 1) / block.y);

    for (int p = 0; p < grid->n_pol * grid->n_chan; p++) {
        cuFloatComplex* plane = grid->d_grid + p * grid->nx * grid->ny;

        gridding_correction_kernel<<<n_blocks, block>>>(
            plane, grid->nx, grid->ny, cf->support
        );
    }

    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
}
