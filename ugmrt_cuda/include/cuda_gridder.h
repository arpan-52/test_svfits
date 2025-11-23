/**
 * @file cuda_gridder.h
 * @brief CUDA gridding functions
 *
 * This module provides two gridding implementations:
 *   1. Simple gridding with real-valued CF (backward compatible)
 *   2. HPG-compatible gridding with complex 6D CF (matches libhpg exactly)
 *
 * HPG compatibility includes:
 *   - 6D complex CF format [x_major, y_major, mueller, cube, x_minor, y_minor]
 *   - UV coordinate computation matching HPG compute_vis_coord
 *   - W-term conjugation (cf_im_factor = pos_w ? -1 : 1)
 *   - Phase screen support
 *   - Phasor application (vis * phasor * weight)
 *   - Weight accumulation as sum(|CF|) * vis_weight
 */

#ifndef CUDA_GRIDDER_H
#define CUDA_GRIDDER_H

#include "cuda_types.h"

#ifdef __cplusplus
extern "C" {
#endif

//-----------------------------------------------------------------------------
// Simple CF management (backward compatibility)
//-----------------------------------------------------------------------------

/**
 * @brief Create simple real-valued convolution function
 * @param support Half-width in pixels
 * @param oversampling Oversampling factor
 * @param values CF values [(2*support+1)*oversamp]^2, row-major
 * @return ConvolutionFunction with data copied to GPU
 */
ConvolutionFunction cf_create(int support, int oversampling, const float* values);

/**
 * @brief Load simple convolution function from binary file
 *
 * File format:
 *   int32: support
 *   int32: oversampling
 *   float[]: values (full_size * full_size)
 */
ConvolutionFunction cf_load(const char* filename);

/**
 * @brief Free simple convolution function
 */
void cf_free(ConvolutionFunction* cf);

//-----------------------------------------------------------------------------
// HPG-compatible CF management
//-----------------------------------------------------------------------------

/**
 * @brief Create HPG-compatible complex 6D convolution function
 *
 * HPG CF layout: [x_major, y_major, mueller, cube, x_minor, y_minor]
 *
 * @param support Half-width in pixels (cf_radius)
 * @param oversampling Oversampling factor
 * @param padding CF padding (typically 0 or oversampling)
 * @param n_mueller Number of Mueller elements (polarization products)
 * @param n_cube Number of CF cubes (e.g., W-planes)
 * @param values Complex CF values in HPG layout
 * @return HPGConvolutionFunction with data copied to GPU
 */
HPGConvolutionFunction hpg_cf_create(
    int support,
    int oversampling,
    int padding,
    int n_mueller,
    int n_cube,
    const cuFloatComplex* values
);

/**
 * @brief Load HPG-compatible CF from binary file
 *
 * File format:
 *   int32: support
 *   int32: oversampling
 *   int32: padding
 *   int32: n_mueller
 *   int32: n_cube
 *   complex[]: values in HPG 6D layout
 */
HPGConvolutionFunction hpg_cf_load(const char* filename);

/**
 * @brief Free HPG convolution function
 */
void hpg_cf_free(HPGConvolutionFunction* cf);

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
 * @brief Simple nearest-neighbor gridding (no CF, for debugging)
 *
 * Grids visibilities using nearest-neighbor without convolution function.
 * Use this to debug the pipeline without CF complications.
 */
void grid_visibilities_simple(
    UVGrid* grid,
    const CudaVisibility* vis,
    int n_vis,
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

//-----------------------------------------------------------------------------
// HPG-compatible gridding operations
//-----------------------------------------------------------------------------

/**
 * @brief Grid visibilities using HPG-compatible algorithm
 *
 * This function matches HPG (libhpg) gridding exactly:
 *   - UV coords: position = grid_scale * coord * inv_lambda + grid_size/2
 *   - CF index computation with padding and fine offset
 *   - W-term conjugation: cf.imag *= (w > 0) ? -1 : 1
 *   - Phase screen: cphase(phi_X + phi_Y)
 *   - Phasor: vis * phasor * weight before scattering
 *   - Weight: sum(|CF|) * vis_weight
 *
 * @param grid UV grid to accumulate into
 * @param vis Array of visibilities (with HPG fields populated)
 * @param n_vis Number of visibilities
 * @param cf HPG-compatible convolution function
 * @param scale_u Grid scale factor for U (radians -> pixels)
 * @param scale_v Grid scale factor for V (radians -> pixels)
 */
void grid_visibilities_hpg(
    UVGrid* grid,
    const CudaVisibility* vis,
    int n_vis,
    const HPGConvolutionFunction* cf,
    float scale_u,
    float scale_v
);

/**
 * @brief Apply gridding correction using HPG CF
 */
void grid_correct_hpg(UVGrid* grid, const HPGConvolutionFunction* cf);

//-----------------------------------------------------------------------------
// Batch gridding (for parallel processing pipeline)
//-----------------------------------------------------------------------------

// Forward declaration (actual definition in svfits_reader.h)
struct VisibilityBuffer;

/**
 * @brief Grid entire visibility buffer at once (batch mode)
 *
 * This is more efficient than incremental gridding when you have
 * all visibilities available upfront. The buffer is transferred
 * to GPU in one operation.
 *
 * @param grid UV grid to accumulate into
 * @param vis_array Visibility array (host memory)
 * @param n_vis Number of visibilities
 * @param scale_u Grid scale factor for U
 * @param scale_v Grid scale factor for V
 * @param use_cf If false, use simple nearest-neighbor gridding
 * @param cf Convolution function (ignored if use_cf is false)
 */
void grid_batch(
    UVGrid* grid,
    const CudaVisibility* vis_array,
    size_t n_vis,
    float scale_u,
    float scale_v,
    int use_cf,
    const ConvolutionFunction* cf
);

#ifdef __cplusplus
}
#endif

#endif // CUDA_GRIDDER_H
