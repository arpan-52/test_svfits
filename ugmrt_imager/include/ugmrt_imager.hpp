/**
 * @file ugmrt_imager.hpp
 * @brief uGMRT Fast Imager - Direct visibility to image pipeline
 *
 * This module integrates svfits visibility processing with HPG GPU gridding
 * to produce images directly from raw uGMRT SPOTLIGHT data without
 * intermediate FITS files.
 *
 * Pipeline: Raw -> Burst Find -> Bandpass -> RFI -> Baseline -> Grid -> FFT -> Image
 */

#ifndef UGMRT_IMAGER_HPP
#define UGMRT_IMAGER_HPP

#include <string>
#include <vector>
#include <array>
#include <complex>
#include <memory>

// Forward declarations for HPG types
namespace hpg {
    class Gridder;
    class GridValueArray;
    class CFArray;
    enum class Device;
}

namespace ugmrt {

//-----------------------------------------------------------------------------
// Configuration structures
//-----------------------------------------------------------------------------

/**
 * @brief Grid configuration parameters
 */
struct GridConfig {
    int nx = 512;                    // Grid size X
    int ny = 512;                    // Grid size Y
    int n_chan = 1;                  // Number of frequency channels in output
    int n_pol = 2;                   // Number of polarizations (1=Stokes I, 2=RR/LL, 4=full)
    double cell_size_asec = 1.0;     // Cell size in arcseconds
    int cf_support = 7;              // Convolution function support (pixels)
    int cf_oversampling = 128;       // Convolution function oversampling
};

/**
 * @brief Imaging parameters
 */
struct ImagingConfig {
    bool natural_weighting = true;   // Natural weighting (vs uniform)
    double robust = 0.0;             // Briggs robust parameter (-2 to 2)
    bool do_fft = true;              // Apply FFT after gridding
    bool do_shift = true;            // Shift image to center
    std::string output_fits;         // Output FITS filename
    std::string output_png;          // Optional PNG output
};

/**
 * @brief Processing flags (mirrors svfits options)
 */
struct ProcessingConfig {
    bool do_bandpass = true;         // Apply bandpass correction
    bool do_baseline = true;         // Subtract baseline (off-source mean)
    bool do_flag = true;             // Apply MAD-based RFI flagging
    float flag_threshold = 5.0f;     // Flagging threshold in MAD units
    int num_threads = 4;             // OpenMP threads for file reading
    int vis_batch_size = 10000;      // Visibilities to batch before GPU transfer
};

/**
 * @brief Complete imager configuration
 */
struct ImagerConfig {
    // Input parameters (from svfits_par.txt)
    std::string param_file;          // Parameter file path
    std::string antsamp_file;        // Antenna/sampler header file
    std::string bulletin_a;          // Optional Bulletin A for DUT1

    // Grid and imaging config
    GridConfig grid;
    ImagingConfig imaging;
    ProcessingConfig processing;

    // HPG device selection
    std::string device = "cuda";     // "cuda", "openmp", or "serial"
};

//-----------------------------------------------------------------------------
// Visibility data structure (for passing to HPG)
//-----------------------------------------------------------------------------

/**
 * @brief Processed visibility ready for gridding
 */
struct ProcessedVisibility {
    std::complex<float> vis[2];      // Visibility values per polarization
    float weight[2];                 // Weights per polarization (negative = flagged)
    double u, v, w;                  // UVW coordinates in wavelengths
    double freq_hz;                  // Frequency in Hz
    int channel;                     // Channel index
    int baseline;                    // Baseline index
    double time_mjd;                 // Timestamp (MJD)
};

//-----------------------------------------------------------------------------
// Main imager class
//-----------------------------------------------------------------------------

/**
 * @brief uGMRT Fast Imager - processes raw visibility data directly to images
 */
class UGMRTImager {
public:
    /**
     * @brief Construct imager with configuration
     */
    explicit UGMRTImager(const ImagerConfig& config);

    /**
     * @brief Destructor
     */
    ~UGMRTImager();

    /**
     * @brief Initialize the imager (loads config, sets up HPG)
     * @return true on success
     */
    bool initialize();

    /**
     * @brief Process all input files and produce image
     * @return true on success
     */
    bool run();

    /**
     * @brief Get the final image data
     * @return Pointer to image data (nx * ny * n_pol * n_chan complex values)
     */
    const std::complex<double>* get_image() const;

    /**
     * @brief Get image dimensions
     */
    std::array<int, 4> get_image_shape() const;

    /**
     * @brief Write image to FITS file
     */
    bool write_fits(const std::string& filename);

    /**
     * @brief Get processing statistics
     */
    struct Stats {
        size_t total_visibilities = 0;
        size_t flagged_visibilities = 0;
        size_t gridded_visibilities = 0;
        double read_time_sec = 0;
        double process_time_sec = 0;
        double grid_time_sec = 0;
        double fft_time_sec = 0;
    };
    Stats get_stats() const;

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

//-----------------------------------------------------------------------------
// Utility functions
//-----------------------------------------------------------------------------

/**
 * @brief Create prolate spheroidal convolution function
 */
std::unique_ptr<hpg::CFArray> create_pswf_cf(
    int support,
    int oversampling,
    int n_mueller,
    int n_cubes
);

/**
 * @brief Convert frequency to wavelength
 */
inline double freq_to_wavelength(double freq_hz) {
    constexpr double c = 299792458.0;  // Speed of light m/s
    return c / freq_hz;
}

/**
 * @brief Compute grid scale from cell size and grid dimensions
 */
inline std::array<double, 2> compute_grid_scale(
    double cell_size_rad,
    int nx, int ny
) {
    // Grid scale = 1 / (n_pixels * cell_size)
    return {1.0 / (nx * cell_size_rad), -1.0 / (ny * cell_size_rad)};
}

} // namespace ugmrt

#endif // UGMRT_IMAGER_HPP
