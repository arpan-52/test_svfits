/**
 * @file cf_generator.h
 * @brief W-projection convolution function generator
 *
 * Implements standard W-projection algorithm following CASA/libhpg approach:
 * 1. Generate W-screen in image domain: W(l,m) = exp(2πi * w * (√(1 - l² - m²) - 1))
 * 2. FFT to UV domain to get convolution kernel
 * 3. Apply prolate spheroidal anti-aliasing
 * 4. Store in HPG 6D CF format
 */

#ifndef CF_GENERATOR_H
#define CF_GENERATOR_H

#include "cuda_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Configuration for CF generation
 */
typedef struct {
    double ra_j2000;        // Phase center RA (J2000, radians)
    double dec_j2000;       // Phase center Dec (J2000, radians)
    double freq_hz;         // Reference frequency (Hz)
    double cell_rad;        // Cell size (radians)
    int grid_nx, grid_ny;   // Grid dimensions

    int n_w_planes;         // Number of W-planes
    float max_w;            // Maximum W value (wavelengths)

    int support;            // CF half-width (pixels)
    int oversampling;       // CF oversampling factor
    int padding;            // CF padding (typically 0)

    int n_mueller;          // Number of Mueller elements (1 for single pol)
} CFGeneratorConfig;

/**
 * @brief Generate W-projection convolution functions on GPU
 *
 * Generates nW W-plane convolution functions following standard W-projection.
 * CFs are stored in HPG 6D format: [x_major, y_major, mueller, cube, x_minor, y_minor]
 *
 * @param config CF generation configuration
 * @return Allocated HPGConvolutionFunction with W-projection CFs on GPU
 */
HPGConvolutionFunction* cf_generate_w_projection(const CFGeneratorConfig* config);

/**
 * @brief Select W-plane index for a given W coordinate
 *
 * Maps W value to appropriate CF cube index.
 *
 * @param w W coordinate (wavelengths)
 * @param max_w Maximum W value (wavelengths)
 * @param n_w_planes Number of W-planes
 * @return CF cube index [0, n_w_planes-1]
 */
int cf_select_w_plane(float w, float max_w, int n_w_planes);

/**
 * @brief Compute maximum W from antenna positions
 *
 * Calculates maximum baseline W coordinate from antenna array geometry.
 *
 * @param antenna_x X positions (meters)
 * @param antenna_y Y positions (meters)
 * @param antenna_z Z positions (meters)
 * @param n_antennas Number of antennas
 * @param dec_rad Source declination (radians)
 * @param wavelength Wavelength (meters)
 * @return Maximum W coordinate (wavelengths)
 */
float cf_compute_max_w(const double* antenna_x, const double* antenna_y, const double* antenna_z,
                       int n_antennas, double dec_rad, double wavelength);

#ifdef __cplusplus
}
#endif

#endif // CF_GENERATOR_H
