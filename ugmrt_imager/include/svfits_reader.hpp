/**
 * @file svfits_reader.hpp
 * @brief C++ wrapper for svfits visibility reading and processing
 *
 * This provides a clean C++ interface to the svfits C code for:
 * - Reading raw uGMRT SPOTLIGHT visibility data
 * - Burst region identification
 * - Bandpass calibration
 * - Baseline subtraction
 * - RFI flagging
 * - UVW coordinate computation
 */

#ifndef SVFITS_READER_HPP
#define SVFITS_READER_HPP

#include <string>
#include <vector>
#include <functional>
#include <complex>
#include <memory>

namespace ugmrt {

// Forward declaration
struct ProcessedVisibility;

/**
 * @brief Burst parameters from parameter file
 */
struct BurstParams {
    std::string name;
    double mjd;              // Burst MJD
    double time_ist;         // Time in IST (seconds)
    double dm;               // Dispersion measure
    double intrinsic_width;  // Intrinsic width (seconds)
    double ref_freq_hz;      // Reference frequency
    double ra_app;           // Apparent RA (radians)
    double dec_app;          // Apparent Dec (radians)
};

/**
 * @brief Frequency setup
 */
struct FreqSetup {
    double freq_start_hz;    // Start frequency
    double freq_end_hz;      // End frequency
    int n_channels;          // Number of channels
    double channel_width_hz; // Channel width

    double center_freq_hz() const {
        return (freq_start_hz + freq_end_hz) / 2.0;
    }
};

/**
 * @brief Antenna information
 */
struct AntennaInfo {
    int id;
    std::string name;
    double x, y, z;          // ECEF coordinates (meters)
};

/**
 * @brief Baseline information
 */
struct BaselineInfo {
    int ant0, ant1;
    int band0, band1;
    bool flip_sign;          // Flip imaginary part
};

/**
 * @brief Processing options
 */
struct ReaderOptions {
    bool do_bandpass = true;
    bool do_baseline = true;
    bool do_flag = true;
    float flag_threshold = 5.0f;
    int num_threads = 4;
    unsigned int antenna_mask = 0xFFFFFFFF;  // All antennas
};

/**
 * @brief Callback for processed visibilities
 *
 * Called for each visibility after all svfits processing is applied.
 * Return false to stop processing.
 */
using VisibilityCallback = std::function<bool(const ProcessedVisibility&)>;

/**
 * @brief svfits raw data reader with full processing pipeline
 */
class SvfitsReader {
public:
    /**
     * @brief Construct reader
     * @param param_file Path to svfits_par.txt
     * @param antsamp_file Path to antsamp.hdr
     * @param bulletin_a Optional path to Bulletin A file
     */
    SvfitsReader(
        const std::string& param_file,
        const std::string& antsamp_file,
        const std::string& bulletin_a = ""
    );

    ~SvfitsReader();

    /**
     * @brief Initialize reader (parse config, allocate buffers)
     */
    bool initialize();

    /**
     * @brief Get burst parameters
     */
    const BurstParams& get_burst_params() const;

    /**
     * @brief Get frequency setup
     */
    const FreqSetup& get_freq_setup() const;

    /**
     * @brief Get antenna list
     */
    const std::vector<AntennaInfo>& get_antennas() const;

    /**
     * @brief Get baseline list
     */
    const std::vector<BaselineInfo>& get_baselines() const;

    /**
     * @brief Get total number of files
     */
    int get_num_files() const;

    /**
     * @brief Get records per file that contain burst
     */
    int get_burst_records() const;

    /**
     * @brief Process all data and call callback for each visibility
     *
     * This executes the full svfits pipeline:
     * 1. Read raw data from all 16 time-multiplexed files
     * 2. Identify burst region
     * 3. Apply bandpass correction
     * 4. Apply baseline subtraction
     * 5. Apply RFI flagging
     * 6. Compute UVW coordinates
     * 7. Call callback for each visibility
     *
     * @param callback Function to call for each processed visibility
     * @param options Processing options
     * @return Number of visibilities processed
     */
    size_t process(
        VisibilityCallback callback,
        const ReaderOptions& options = ReaderOptions()
    );

    /**
     * @brief Get all processed visibilities as a vector
     * (convenience wrapper around process())
     */
    std::vector<ProcessedVisibility> get_all_visibilities(
        const ReaderOptions& options = ReaderOptions()
    );

    /**
     * @brief Get reference wavelength (meters)
     */
    double get_ref_wavelength() const;

    /**
     * @brief Get phase center (J2000 radians)
     */
    std::pair<double, double> get_phase_center() const;

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace ugmrt

#endif // SVFITS_READER_HPP
