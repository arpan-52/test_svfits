/**
 * @file cuda_gridder.h
 * @brief CUDA gridding functions
 */

#ifndef CUDA_GRIDDER_H
#define CUDA_GRIDDER_H

#include "cuda_types.h"

#ifdef __cplusplus
extern "C" {
#endif

//-----------------------------------------------------------------------------
// Convolution function management
//-----------------------------------------------------------------------------

/**
 * @brief Create convolution function from data
 * @param support Half-width in pixels
 * @param oversampling Oversampling factor
 * @param values CF values [(2*support+1)*oversamp]^2, row-major
 * @return ConvolutionFunction with data copied to GPU
 */
ConvolutionFunction cf_create(int support, int oversampling, const float* values);

/**
 * @brief Load convolution function from binary file
 *
 * File format:
 *   int32: support
 *   int32: oversampling
 *   float[]: values (full_size * full_size)
 */
ConvolutionFunction cf_load(const char* filename);

/**
 * @brief Free convolution function
 */
void cf_free(ConvolutionFunction* cf);

//-----------------------------------------------------------------------------
// Grid management
//-----------------------------------------------------------------------------

/**
 * @brief Create UV grid
 */
UVGrid grid_create(int nx, int ny, int n_pol, int n_chan);

/**
 * @brief Zero the grid
 */
void grid_zero(UVGrid* grid);

/**
 * @brief Free grid
 */
void grid_free(UVGrid* grid);

/**
 * @brief Copy grid from GPU to host
 */
void grid_to_host(UVGrid* grid);

//-----------------------------------------------------------------------------
// Gridding operations
//-----------------------------------------------------------------------------

/**
 * @brief Grid a batch of visibilities
 *
 * @param grid UV grid to accumulate into
 * @param vis Array of visibilities
 * @param n_vis Number of visibilities
 * @param cf Convolution function
 * @param grid_scale Grid scale factors [scale_u, scale_v] (wavelengths to pixels)
 */
void grid_visibilities(
    UVGrid* grid,
    const CudaVisibility* vis,
    int n_vis,
    const ConvolutionFunction* cf,
    float scale_u,
    float scale_v
);

/**
 * @brief Normalize grid by weights
 *
 * grid[i] /= weights[i] for all pixels with non-zero weight
 */
void grid_normalize(UVGrid* grid);

/**
 * @brief Apply FFT to grid (in-place)
 *
 * Uses cuFFT for 2D complex-to-complex FFT
 * @param grid UV grid
 * @param forward true for forward FFT (UV->image), false for inverse
 */
void grid_fft(UVGrid* grid, int forward);

/**
 * @brief Shift grid to center the image (fftshift)
 */
void grid_shift(UVGrid* grid);

/**
 * @brief Apply gridding correction to image
 *
 * Divides by the FT of the convolution function to correct for
 * the gridding convolution.
 *
 * @param grid Image grid (after FFT)
 * @param cf Convolution function (for correction computation)
 */
void grid_correct(UVGrid* grid, const ConvolutionFunction* cf);

#ifdef __cplusplus
}
#endif

#endif // CUDA_GRIDDER_H
