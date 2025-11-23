/**
 * @file cuda_gridder.cu
 * @brief CUDA kernels for visibility gridding
 *
 * This implementation matches HPG (Hyperion Polyphase Gridder) exactly:
 *   - UV coordinate computation: position = grid_scale * coord + grid_size/2 (UV in wavelengths)
 *   - CF indexing: 6D [x_major, y_major, mueller, cube, x_minor, y_minor]
 *   - W-term conjugation: cf_im_factor = (pos_w ? -1 : 1) for gridding
 *   - Phase screen: cphase(phi_X + phi_Y) applied to CF
 *   - Phasor: vis * phasor * weight before gridding
 *   - Weights: sum(|CF| * vis_weight) per polarization
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cuda_runtime.h>
#include <cufft.h>

#include "cuda_types.h"
#include "cuda_gridder.h"

// Note: UV coordinates are passed in wavelengths, no C_LIGHT conversion needed

//=============================================================================
// Device helper functions (matching HPG)
//=============================================================================

/**
 * @brief Compute phase -> complex (matches HPG cphase)
 */
__device__ __forceinline__ cuFloatComplex cphase(float ph) {
    float sn, cs;
    sincosf(ph, &sn, &cs);
    return make_cuFloatComplex(cs, sn);
}

/**
 * @brief Complex magnitude
 */
__device__ __forceinline__ float cmag(cuFloatComplex c) {
    return hypotf(c.x, c.y);
}

/**
 * @brief Complex multiply
 */
__device__ __forceinline__ cuFloatComplex cmul(cuFloatComplex a, cuFloatComplex b) {
    return make_cuFloatComplex(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
}

/**
 * @brief Complex multiply-add: a * b + c
 */
__device__ __forceinline__ cuFloatComplex cfma(cuFloatComplex a, cuFloatComplex b, cuFloatComplex c) {
    return make_cuFloatComplex(a.x * b.x - a.y * b.y + c.x, a.x * b.y + a.y * b.x + c.y);
}

/**
 * @brief Compute visibility coordinate (matches HPG compute_vis_coord)
 *
 * @param g_size      Grid size
 * @param oversampling CF oversampling factor
 * @param cf_padding  CF padding (typically oversampling or 0)
 * @param cf_radius   CF half-width (support)
 * @param coord       UV coordinate in wavelengths
 * @param grid_scale  Grid scale factor (N * cell_rad)
 * @param grid_coord  Output: grid coordinate (leftmost of CF support)
 * @param cf_major    Output: CF major index
 * @param cf_minor    Output: CF minor (oversampling) index
 * @param fine_offset Output: visibility offset to nearest grid point
 */
__device__ void compute_vis_coord_hpg(
    int g_size,
    int oversampling,
    int cf_padding,
    int cf_radius,
    float coord,        // UV coordinate already in wavelengths
    float grid_scale,
    int* grid_coord,
    int* cf_major,
    int* cf_minor,
    int* fine_offset)
{
    // Position on grid: coord (wavelengths) * scale + center
    float position = grid_scale * coord + g_size / 2.0f;

    // Nearest grid point
    int g = __float2int_rn(position);  // round to nearest

    // Fine offset (visibility to nearest grid, in oversampled units)
    int fo = __float2int_rn((g - position) * oversampling);

    // Adjust grid coordinate to be leftmost of CF support
    *grid_coord = g - cf_radius;

    // Compute CF minor and major indices (matching HPG)
    if (fo >= 0) {
        *cf_minor = fo;
        *cf_major = cf_padding;
    } else {
        *cf_minor = oversampling + fo;
        *cf_major = cf_padding - 1;
    }

    *fine_offset = fo;
}

//=============================================================================
// HPG-Compatible Gridding Kernel
//=============================================================================

/**
 * @brief HPG-compatible visibility gridding kernel
 *
 * This kernel matches HPG's gridding algorithm exactly:
 *   - UV coordinate computation
 *   - 6D CF indexing
 *   - W-conjugation
 *   - Phase screen
 *   - Phasor application
 *   - Weight accumulation
 */
__global__ void grid_visibility_hpg_kernel(
    cuFloatComplex* __restrict__ grid,
    float* __restrict__ weights,
    const CudaVisibility* __restrict__ vis,
    int n_vis,
    const cuFloatComplex* __restrict__ cf,
    int cf_support,         // Half-width (radius)
    int cf_oversampling,
    int cf_padding,
    int cf_n_mueller,
    int cf_n_cube,
    int cf_stride_x_major,
    int cf_stride_y_major,
    int cf_stride_mueller,
    int cf_stride_cube,
    int cf_stride_x_minor,
    int cf_stride_y_minor,
    int nx, int ny,
    float scale_u, float scale_v)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= n_vis) return;

    CudaVisibility v = vis[tid];

    // Skip flagged visibilities
    if (v.weight <= 0.0f) return;

    // CF size (full width)
    int cf_size = 2 * cf_support + 1;

    // Compute grid and CF coordinates for U (UV already in wavelengths)
    int grid_u, cf_major_u, cf_minor_u, fine_offset_u;
    compute_vis_coord_hpg(
        nx, cf_oversampling, cf_padding, cf_support,
        v.u, scale_u,
        &grid_u, &cf_major_u, &cf_minor_u, &fine_offset_u);

    // Compute grid and CF coordinates for V
    int grid_v, cf_major_v, cf_minor_v, fine_offset_v;
    compute_vis_coord_hpg(
        ny, cf_oversampling, cf_padding, cf_support,
        v.v, scale_v,
        &grid_v, &cf_major_v, &cf_minor_v, &fine_offset_v);

    // Check if visibility is within grid bounds
    if (grid_u < 0 || grid_u + cf_size > nx) return;
    if (grid_v < 0 || grid_v + cf_size > ny) return;

    // W-term conjugation factor (matching HPG)
    // For gridding: cf_im_factor = (w > 0) ? -1 : 1
    float cf_im_factor = (v.w > 0.0f) ? -1.0f : 1.0f;

    // Compute phasor from d_phase (matching HPG)
    cuFloatComplex phasor = cphase(v.d_phase);

    // Apply phasor and weight to visibility (matching HPG: vis * phasor * weight)
    cuFloatComplex vis_val = make_cuFloatComplex(v.re, v.im);
    vis_val = cmul(vis_val, phasor);
    vis_val.x *= v.weight;
    vis_val.y *= v.weight;

    // Phase screen origin and increment (matching HPG)
    float phi0_u = -v.phase_grad_u * (cf_support * cf_oversampling - fine_offset_u);
    float phi0_v = -v.phase_grad_v * (cf_support * cf_oversampling - fine_offset_v);
    float dphi_u = v.phase_grad_u * cf_oversampling;
    float dphi_v = v.phase_grad_v * cf_oversampling;

    // CF cube index
    int cf_cube = v.cf_cube;
    if (cf_cube < 0 || cf_cube >= cf_n_cube) cf_cube = 0;

    // Mueller index (assume 0 for single polarization)
    int mueller_idx = 0;

    // Weight accumulator (sum of |CF|)
    float wt_sum = 0.0f;

    // Scatter to grid using CF (matching HPG grid_vis)
    for (int Y = 0; Y < cf_size; Y++) {
        int gv = grid_v + Y;
        float phi_Y = phi0_v + Y * dphi_v;

        // CF y major index
        int cf_y_major = cf_major_v + Y;

        for (int X = 0; X < cf_size; X++) {
            int gu = grid_u + X;
            float phi_X = phi0_u + X * dphi_u;

            // CF x major index
            int cf_x_major = cf_major_u + X;

            // Phase screen (matching HPG)
            cuFloatComplex screen = cphase(phi_X + phi_Y);

            // 6D CF index (matching HPG cf_view layout)
            int cf_idx = cf_x_major * cf_stride_x_major +
                         cf_y_major * cf_stride_y_major +
                         mueller_idx * cf_stride_mueller +
                         cf_cube * cf_stride_cube +
                         cf_minor_u * cf_stride_x_minor +
                         cf_minor_v * cf_stride_y_minor;

            // Get CF value
            cuFloatComplex cfv = cf[cf_idx];

            // Apply W-conjugation (matching HPG)
            cfv.y *= cf_im_factor;

            // Compute gridded value: cf * screen * vis (matching HPG)
            cuFloatComplex gv_contrib = cmul(cfv, screen);
            gv_contrib = cmul(gv_contrib, vis_val);

            // Grid index (2D for single pol/cube)
            int grid_idx = gv * nx + gu;

            // Atomic add to grid
            atomicAdd(&grid[grid_idx].x, gv_contrib.x);
            atomicAdd(&grid[grid_idx].y, gv_contrib.y);

            // Accumulate weight (|CF|)
            wt_sum += cmag(cfv);
        }
    }

    // Add weight contribution (matching HPG: sum(|CF|) * vis_weight)
    // HPG accumulates weights per polarization, we use single weight grid
    float total_weight = wt_sum * v.weight;
    atomicAdd(&weights[v.grid_cube], total_weight);
}

//=============================================================================
// Legacy Simple Gridding Kernel (backward compatibility)
//=============================================================================

/**
 * @brief Simple grid visibilities kernel (original implementation)
 *
 * Each thread processes one visibility and scatters it to the grid
 * using a simple real-valued convolution function.
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

    // UV coordinates are already in wavelengths - grid directly
    // grid = scale * u + N/2, where scale = N * cell_rad
    float grid_u = scale_u * v.u + nx / 2.0f;
    float grid_v = scale_v * v.v + ny / 2.0f;

    // Integer grid position (nearest)
    int iu = __float2int_rn(grid_u);
    int iv = __float2int_rn(grid_v);

    // Fractional offset (fine offset in HPG terms)
    int fine_u = __float2int_rn((iu - grid_u) * cf_oversampling);
    int fine_v = __float2int_rn((iv - grid_v) * cf_oversampling);

    // CF minor index
    int cf_minor_u = (fine_u >= 0) ? fine_u : cf_oversampling + fine_u;
    int cf_minor_v = (fine_v >= 0) ? fine_v : cf_oversampling + fine_v;

    // Apply phasor (HPG-style)
    cuFloatComplex phasor = cphase(v.d_phase);
    cuFloatComplex vis_val = make_cuFloatComplex(v.re, v.im);
    vis_val = cmul(vis_val, phasor);
    vis_val.x *= v.weight;
    vis_val.y *= v.weight;

    // W-conjugation factor
    float w_conj = (v.w > 0.0f) ? -1.0f : 1.0f;

    // Scatter to grid using CF
    for (int dv = -cf_support; dv <= cf_support; dv++) {
        int gv = iv + dv;
        if (gv < 0 || gv >= ny) continue;

        // CF y index
        int cf_y = (dv + cf_support) * cf_oversampling + cf_minor_v;

        for (int du = -cf_support; du <= cf_support; du++) {
            int gu = iu + du;
            if (gu < 0 || gu >= nx) continue;

            // CF x index
            int cf_x = (du + cf_support) * cf_oversampling + cf_minor_u;

            // CF value (real)
            float cf_val = cf[cf_y * cf_full_size + cf_x];

            // Grid index
            int grid_idx = gv * nx + gu;

            // Apply CF to visibility
            cuFloatComplex weighted_vis = make_cuFloatComplex(
                vis_val.x * cf_val,
                vis_val.y * cf_val * w_conj  // W-conjugation on imaginary
            );

            // Atomic add to grid
            atomicAdd(&grid[grid_idx].x, weighted_vis.x);
            atomicAdd(&grid[grid_idx].y, weighted_vis.y);
            atomicAdd(&weights[grid_idx], fabsf(cf_val) * v.weight);
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
 * @brief Scale grid by a constant factor (for FFT normalization)
 */
__global__ void scale_kernel(
    cuFloatComplex* __restrict__ grid,
    int n_pixels,
    float scale)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= n_pixels) return;

    grid[tid].x *= scale;
    grid[tid].y *= scale;
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

    // Normalize if inverse FFT (cuFFT doesn't include 1/N factor)
    if (!forward) {
        int n_pixels = grid->nx * grid->ny * grid->n_pol * grid->n_chan;
        float norm = 1.0f / (grid->nx * grid->ny);

        int block_size = 256;
        int n_blocks = (n_pixels + block_size - 1) / block_size;

        scale_kernel<<<n_blocks, block_size>>>(grid->d_grid, n_pixels, norm);
        CUDA_CHECK(cudaGetLastError());
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

//=============================================================================
// HPG-Compatible CF Management
//=============================================================================

HPGConvolutionFunction hpg_cf_create(
    int support,
    int oversampling,
    int padding,
    int n_mueller,
    int n_cube,
    const cuFloatComplex* values)
{
    HPGConvolutionFunction cf;
    cf.support = support;
    cf.oversampling = oversampling;
    cf.padding = padding;
    cf.n_mueller = n_mueller;
    cf.n_cube = n_cube;
    cf.cf_size = 2 * support + 1;

    // HPG layout: [x_major, y_major, mueller, cube, x_minor, y_minor]
    // Row-major, fastest varying last
    int x_major_size = cf.cf_size + 2 * padding;
    int y_major_size = cf.cf_size + 2 * padding;

    cf.stride_y_minor = 1;
    cf.stride_x_minor = oversampling;
    cf.stride_cube = oversampling * oversampling;
    cf.stride_mueller = n_cube * cf.stride_cube;
    cf.stride_y_major = n_mueller * cf.stride_mueller;
    cf.stride_x_major = y_major_size * cf.stride_y_major;

    cf.n_values = x_major_size * cf.stride_x_major;

    size_t size = cf.n_values * sizeof(cuFloatComplex);

    // Allocate on GPU
    CUDA_CHECK(cudaMalloc(&cf.d_values, size));
    CUDA_CHECK(cudaMemcpy(cf.d_values, values, size, cudaMemcpyHostToDevice));

    // Keep host copy
    cf.h_values = (cuFloatComplex*)malloc(size);
    memcpy(cf.h_values, values, size);

    printf("Created HPG CF: support=%d, oversamp=%d, mueller=%d, cube=%d (%.2f MB)\n",
           support, oversampling, n_mueller, n_cube, size / (1024.0 * 1024.0));

    return cf;
}

HPGConvolutionFunction hpg_cf_load(const char* filename) {
    FILE* fp = fopen(filename, "rb");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open HPG CF file: %s\n", filename);
        exit(EXIT_FAILURE);
    }

    // HPG CF file header
    int support, oversampling, padding, n_mueller, n_cube;
    fread(&support, sizeof(int), 1, fp);
    fread(&oversampling, sizeof(int), 1, fp);
    fread(&padding, sizeof(int), 1, fp);
    fread(&n_mueller, sizeof(int), 1, fp);
    fread(&n_cube, sizeof(int), 1, fp);

    // Calculate size
    int cf_size = 2 * support + 1;
    int x_major_size = cf_size + 2 * padding;
    int y_major_size = cf_size + 2 * padding;

    size_t n_values = (size_t)x_major_size * y_major_size * n_mueller * n_cube *
                      oversampling * oversampling;

    cuFloatComplex* values = (cuFloatComplex*)malloc(n_values * sizeof(cuFloatComplex));
    size_t read_count = fread(values, sizeof(cuFloatComplex), n_values, fp);
    fclose(fp);

    if (read_count != n_values) {
        fprintf(stderr, "Warning: Expected %zu values, read %zu\n", n_values, read_count);
    }

    HPGConvolutionFunction cf = hpg_cf_create(support, oversampling, padding, n_mueller, n_cube, values);
    free(values);

    printf("Loaded HPG CF from %s\n", filename);

    return cf;
}

void hpg_cf_free(HPGConvolutionFunction* cf) {
    if (cf->d_values) {
        cudaFree(cf->d_values);
        cf->d_values = NULL;
    }
    if (cf->h_values) {
        free(cf->h_values);
        cf->h_values = NULL;
    }
}

//=============================================================================
// HPG-Compatible Gridding Operations
//=============================================================================

void grid_visibilities_hpg(
    UVGrid* grid,
    const CudaVisibility* vis,
    int n_vis,
    const HPGConvolutionFunction* cf,
    float scale_u,
    float scale_v)
{
    if (n_vis == 0) return;

    // Copy visibilities to GPU
    CudaVisibility* d_vis;
    CUDA_CHECK(cudaMalloc(&d_vis, n_vis * sizeof(CudaVisibility)));
    CUDA_CHECK(cudaMemcpy(d_vis, vis, n_vis * sizeof(CudaVisibility), cudaMemcpyHostToDevice));

    // Launch HPG-compatible kernel
    int block_size = 256;
    int n_blocks = (n_vis + block_size - 1) / block_size;

    grid_visibility_hpg_kernel<<<n_blocks, block_size>>>(
        grid->d_grid,
        grid->d_weights,
        d_vis,
        n_vis,
        cf->d_values,
        cf->support,
        cf->oversampling,
        cf->padding,
        cf->n_mueller,
        cf->n_cube,
        cf->stride_x_major,
        cf->stride_y_major,
        cf->stride_mueller,
        cf->stride_cube,
        cf->stride_x_minor,
        cf->stride_y_minor,
        grid->nx, grid->ny,
        scale_u, scale_v
    );

    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    cudaFree(d_vis);
}

void grid_correct_hpg(UVGrid* grid, const HPGConvolutionFunction* cf) {
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
