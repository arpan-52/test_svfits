/**
 * @file cuda_types.h
 * @brief Common types for CUDA-based uGMRT imager
 *
 * This implementation is designed to match HPG (Hyperion Polyphase Gridder)
 * exactly in terms of:
 *   - CF format (6D complex array)
 *   - UV coordinate computation
 *   - W-term conjugation
 *   - Weight accumulation
 */

#ifndef CUDA_TYPES_H
#define CUDA_TYPES_H

#include <cuda_runtime.h>
#include <cuComplex.h>

#ifdef __cplusplus
extern "C" {
#endif

//-----------------------------------------------------------------------------
// Visibility data structure (matches HPG VisData)
//-----------------------------------------------------------------------------

/**
 * @brief Single visibility for gridding
 *
 * This matches HPG's VisData structure with the essential fields.
 */
typedef struct {
    float re, im;           // Visibility value (real, imaginary)
    float weight;           // Weight (negative = flagged)
    float u, v, w;          // UVW coordinates in wavelengths
    float freq;             // Frequency in Hz (for computing inv_lambda)
    float d_phase;          // Phase angle for phasor
    int channel;            // Frequency channel index
    int cf_cube;            // CF cube index (e.g., W-plane)
    int cf_grp;             // CF group index (for variable-size CF)
    int grid_cube;          // Grid cube index
    float phase_grad_u;     // Phase gradient in U
    float phase_grad_v;     // Phase gradient in V
} CudaVisibility;

//-----------------------------------------------------------------------------
// Convolution function (HPG-compatible format)
//-----------------------------------------------------------------------------

/**
 * @brief HPG-compatible convolution function
 *
 * HPG CF format is 6D complex array:
 *   [x_major, y_major, mueller, cube, x_minor, y_minor]
 *
 * Where:
 *   - x_major, y_major: Grid cell index (0 to 2*support+1 + padding)
 *   - mueller: Polarization product index
 *   - cube: CF cube (e.g., W-plane, baseline class, etc.)
 *   - x_minor, y_minor: Oversampling index (0 to oversampling-1)
 *
 * For simple use (single cube, single mueller), this reduces to:
 *   cf[x_major][y_major][0][0][x_minor][y_minor]
 *
 * Memory layout (row-major, fastest varying last):
 *   Linear index = x_major * stride_x_major + y_major * stride_y_major + ...
 */
typedef struct {
    int support;            // Half-width in pixels (cf_radius)
    int oversampling;       // Oversampling factor
    int padding;            // CF padding (typically 0 or 1)
    int n_mueller;          // Number of Mueller elements
    int n_cube;             // Number of CF cubes (e.g., W-planes)
    int cf_size;            // Full CF support size (2*support + 1)

    // Strides for 6D indexing (in complex elements)
    int stride_x_major;
    int stride_y_major;
    int stride_mueller;
    int stride_cube;
    int stride_x_minor;
    int stride_y_minor;

    cuFloatComplex* d_values;  // Complex CF values on GPU
    cuFloatComplex* h_values;  // Complex CF values on host
    size_t n_values;           // Total number of complex values
} HPGConvolutionFunction;

/**
 * @brief Simple real-valued CF (for backward compatibility)
 */
typedef struct {
    int support;            // Half-width in pixels (e.g., 7)
    int oversampling;       // Oversampling factor (e.g., 128)
    int full_size;          // (2*support + 1) * oversampling
    float* d_values;        // CF values on GPU [full_size * full_size]
    float* h_values;        // CF values on host (optional)
} ConvolutionFunction;

//-----------------------------------------------------------------------------
// Grid structure
//-----------------------------------------------------------------------------

/**
 * @brief UV grid for accumulating visibilities
 */
typedef struct {
    int nx, ny;             // Grid dimensions
    int n_chan;             // Number of frequency channels
    int n_pol;              // Number of polarizations
    cuFloatComplex* d_grid; // Complex grid on GPU [nx * ny * n_pol * n_chan]
    float* d_weights;       // Weight grid on GPU [nx * ny * n_pol * n_chan]
    cuFloatComplex* h_grid; // Grid on host (for output)
} UVGrid;

//-----------------------------------------------------------------------------
// Gridder configuration
//-----------------------------------------------------------------------------

typedef struct {
    int nx, ny;             // Grid size
    int n_chan;             // Output channels
    int n_pol;              // Output polarizations
    double cell_size_rad;   // Cell size in radians
    double ref_freq_hz;     // Reference frequency
    int batch_size;         // Visibility batch size for GPU
} GridderConfig;

//-----------------------------------------------------------------------------
// Error handling
//-----------------------------------------------------------------------------

#define CUDA_CHECK(call) \
    do { \
        cudaError_t err = call; \
        if (err != cudaSuccess) { \
            fprintf(stderr, "CUDA error at %s:%d: %s\n", \
                    __FILE__, __LINE__, cudaGetErrorString(err)); \
            exit(EXIT_FAILURE); \
        } \
    } while(0)

#define CUFFT_CHECK(call) \
    do { \
        cufftResult err = call; \
        if (err != CUFFT_SUCCESS) { \
            fprintf(stderr, "cuFFT error at %s:%d: %d\n", \
                    __FILE__, __LINE__, err); \
            exit(EXIT_FAILURE); \
        } \
    } while(0)

#ifdef __cplusplus
}
#endif

#endif // CUDA_TYPES_H
