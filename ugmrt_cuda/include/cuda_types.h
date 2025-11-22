/**
 * @file cuda_types.h
 * @brief Common types for CUDA-based uGMRT imager
 */

#ifndef CUDA_TYPES_H
#define CUDA_TYPES_H

#include <cuda_runtime.h>
#include <cuComplex.h>

#ifdef __cplusplus
extern "C" {
#endif

//-----------------------------------------------------------------------------
// Visibility data structure
//-----------------------------------------------------------------------------

/**
 * @brief Single visibility for gridding
 */
typedef struct {
    float re, im;           // Visibility value (real, imaginary)
    float weight;           // Weight (negative = flagged)
    float u, v, w;          // UVW coordinates in wavelengths
    int channel;            // Frequency channel index
} CudaVisibility;

//-----------------------------------------------------------------------------
// Convolution function (gridding kernel)
//-----------------------------------------------------------------------------

/**
 * @brief Convolution function for gridding
 *
 * The CF is stored as a 2D oversampled array.
 * For a support of S and oversampling of O:
 *   - Full size = (2*S + 1) * O per axis
 *   - Values are real (can extend to complex if needed)
 *
 * To look up CF value at fractional position (fx, fy):
 *   idx_x = (int)(fx * oversampling) + support * oversampling
 *   idx_y = (int)(fy * oversampling) + support * oversampling
 *   cf_value = values[idx_y * full_size + idx_x]
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
