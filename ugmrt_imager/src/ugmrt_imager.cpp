/**
 * @file ugmrt_imager.cpp
 * @brief Main uGMRT Fast Imager implementation
 *
 * Integrates svfits reading with HPG GPU gridding for fast imaging
 * of uGMRT SPOTLIGHT burst data.
 */

#include "ugmrt_imager.hpp"
#include "svfits_reader.hpp"
#include "vis_to_hpg.hpp"
#include "cf_kernel.hpp"

#include <hpg/hpg.hpp>
#include <fitsio.h>

#include <iostream>
#include <chrono>
#include <cmath>
#include <algorithm>

namespace ugmrt {

// Forward declarations
bool write_fits_image(
    const std::string& filename,
    const std::complex<double>* image_data,
    int nx, int ny, int n_pol, int n_chan,
    double ra_deg, double dec_deg,
    double cell_size_deg,
    double ref_freq_hz,
    double freq_width_hz);

void apply_gridding_correction(
    std::complex<double>* image,
    int nx, int ny,
    int support);

//-----------------------------------------------------------------------------
// Implementation class
//-----------------------------------------------------------------------------

class UGMRTImager::Impl {
public:
    ImagerConfig config_;
    std::unique_ptr<SvfitsReader> reader_;
    std::unique_ptr<hpg::Gridder> gridder_;
    std::unique_ptr<PSWFConvolutionFunction> cf_;

    // Image data (after FFT)
    std::vector<std::complex<double>> image_;

    // Statistics
    Stats stats_;

    // Timing helpers
    using Clock = std::chrono::high_resolution_clock;

    Impl(const ImagerConfig& config) : config_(config) {}

    hpg::Device get_device() const {
        if (config_.device == "cuda") {
            return hpg::Device::Cuda;
        } else if (config_.device == "openmp") {
            return hpg::Device::OpenMP;
        } else {
            return hpg::Device::Serial;
        }
    }

    bool initialize() {
        // Initialize HPG
        if (!hpg::is_initialized()) {
            if (!hpg::initialize()) {
                std::cerr << "Failed to initialize HPG" << std::endl;
                return false;
            }
        }

        // Check if requested device is available
        auto available = hpg::devices();
        hpg::Device device = get_device();
        if (available.find(device) == available.end()) {
            std::cerr << "Requested device not available, falling back to Serial" << std::endl;
            device = hpg::Device::Serial;
        }

        // Initialize svfits reader
        reader_ = std::make_unique<SvfitsReader>(
            config_.param_file,
            config_.antsamp_file,
            config_.bulletin_a
        );

        if (!reader_->initialize()) {
            std::cerr << "Failed to initialize svfits reader" << std::endl;
            return false;
        }

        // Get frequency info for grid setup
        const auto& freq = reader_->get_freq_setup();

        // Create convolution function
        cf_ = std::make_unique<PSWFConvolutionFunction>(
            config_.grid.cf_support,
            config_.grid.cf_oversampling,
            config_.grid.n_pol,  // Mueller elements
            config_.grid.n_chan  // Frequency cubes
        );

        // Create CF shape for gridder initialization
        auto cf_shape = create_cf_shape(
            config_.grid.cf_support,
            config_.grid.cf_oversampling,
            config_.grid.n_pol,
            config_.grid.n_chan
        );

        // Compute grid scale from cell size
        double cell_rad = config_.grid.cell_size_asec * M_PI / (180.0 * 3600.0);
        auto grid_scale = compute_grid_scale(cell_rad, config_.grid.nx, config_.grid.ny);

        // Mueller index configuration
        auto mueller_idx = MuellerMapping::get_mueller_indexes(config_.grid.n_pol);
        auto conj_mueller_idx = MuellerMapping::get_conj_mueller_indexes(config_.grid.n_pol);

        // Convert to HPG format
        hpg::IArrayVector mueller_indexes;
        hpg::IArrayVector conj_mueller_indexes;
        for (const auto& idx : mueller_idx) {
            mueller_indexes.push_back({idx[0], idx[1]});
        }
        for (const auto& idx : conj_mueller_idx) {
            conj_mueller_indexes.push_back({idx[0], idx[1]});
        }

        // Create gridder
        std::array<hpg::coord_t, 4> grid_size = {
            config_.grid.nx,
            config_.grid.ny,
            config_.grid.n_pol,
            config_.grid.n_chan
        };

        auto result = hpg::Gridder::create(
            device,
            4,  // max_added_tasks
            config_.processing.vis_batch_size,
            &cf_shape,
            grid_size,
            grid_scale,
            mueller_indexes,
            conj_mueller_indexes
        );

        if (hpg::is_error(result)) {
            std::cerr << "Failed to create gridder: " << hpg::get_error(result).value() << std::endl;
            return false;
        }

        gridder_ = std::make_unique<hpg::Gridder>(hpg::get_value(std::move(result)));

        // Set convolution function
        auto cf_result = gridder_->set_convolution_function(
            hpg::Device::Serial,
            std::move(*cf_)
        );

        if (hpg::is_error(cf_result)) {
            std::cerr << "Failed to set CF" << std::endl;
            return false;
        }

        std::cout << "Initialized uGMRT Imager:" << std::endl;
        std::cout << "  Device: " << config_.device << std::endl;
        std::cout << "  Grid: " << config_.grid.nx << "x" << config_.grid.ny << std::endl;
        std::cout << "  Cell size: " << config_.grid.cell_size_asec << " arcsec" << std::endl;
        std::cout << "  CF support: " << config_.grid.cf_support << " (oversampling: "
                  << config_.grid.cf_oversampling << ")" << std::endl;

        return true;
    }

    bool run() {
        auto start_time = Clock::now();

        // Create visibility batch
        VisibilityBatch batch(config_.processing.vis_batch_size);

        // Get host device for data transfers
        hpg::Device host_device = hpg::Device::Serial;

        // Get frequency setup
        const auto& freq = reader_->get_freq_setup();
        double ref_freq = freq.center_freq_hz();

        // Process all visibilities
        std::cout << "Processing visibilities..." << std::endl;

        ReaderOptions opts;
        opts.do_bandpass = config_.processing.do_bandpass;
        opts.do_baseline = config_.processing.do_baseline;
        opts.do_flag = config_.processing.do_flag;
        opts.flag_threshold = config_.processing.flag_threshold;
        opts.num_threads = config_.processing.num_threads;

        auto read_start = Clock::now();

        size_t processed = reader_->process([&](const ProcessedVisibility& vis) {
            stats_.total_visibilities++;

            // Skip flagged visibilities
            if (vis.weight[0] < 0 || vis.weight[1] < 0) {
                stats_.flagged_visibilities++;
                return true;  // Continue processing
            }

            // Map channel to output grid cube (can do channel averaging here)
            int grid_cube = 0;  // Single channel output for now

            // Add to batch
            bool is_full = batch.add(
                vis.vis[0],
                vis.vis[1],
                vis.weight[0],
                vis.weight[1],
                vis.u, vis.v, vis.w,
                vis.freq_hz,
                grid_cube,
                0,  // cf_cube
                0,  // cf_group
                0.0f,  // phase
                0.0f, 0.0f  // cf_phase_gradient
            );

            // Grid when batch is full
            if (is_full) {
                auto grid_start = Clock::now();

                auto& vis_data = batch.get_vis_data();
                auto grid_result = gridder_->grid_visibilities(
                    host_device,
                    std::move(vis_data),
                    true  // update_grid_weights
                );

                if (hpg::is_error(grid_result)) {
                    std::cerr << "Gridding error" << std::endl;
                    return false;
                }

                stats_.gridded_visibilities += batch.size();

                auto grid_end = Clock::now();
                stats_.grid_time_sec += std::chrono::duration<double>(grid_end - grid_start).count();

                batch.clear();
            }

            return true;
        }, opts);

        // Grid remaining visibilities
        if (!batch.empty()) {
            auto& vis_data = batch.get_vis_data();
            auto grid_result = gridder_->grid_visibilities(
                host_device,
                std::move(vis_data),
                true
            );
            stats_.gridded_visibilities += batch.size();
            batch.clear();
        }

        auto read_end = Clock::now();
        stats_.read_time_sec = std::chrono::duration<double>(read_end - read_start).count();
        stats_.process_time_sec = stats_.read_time_sec;  // Combined for now

        std::cout << "  Processed: " << stats_.total_visibilities << " visibilities" << std::endl;
        std::cout << "  Flagged: " << stats_.flagged_visibilities << std::endl;
        std::cout << "  Gridded: " << stats_.gridded_visibilities << std::endl;

        // Normalize by weights
        std::cout << "Normalizing..." << std::endl;
        gridder_->normalize_by_weights(1.0);

        // Apply FFT
        if (config_.imaging.do_fft) {
            std::cout << "Applying FFT..." << std::endl;
            auto fft_start = Clock::now();

            auto fft_result = gridder_->apply_grid_fft(1.0, hpg::FFTSign::FORWARD, true);
            if (hpg::is_error(fft_result)) {
                std::cerr << "FFT error" << std::endl;
                return false;
            }

            auto fft_end = Clock::now();
            stats_.fft_time_sec = std::chrono::duration<double>(fft_end - fft_start).count();
        }

        // Shift grid
        if (config_.imaging.do_shift) {
            std::cout << "Shifting grid..." << std::endl;
            gridder_->shift_grid(hpg::ShiftDirection::FORWARD);
        }

        // Extract image data
        std::cout << "Extracting image..." << std::endl;
        auto& grid = gridder_->grid();
        size_t total = config_.grid.nx * config_.grid.ny * config_.grid.n_pol * config_.grid.n_chan;
        image_.resize(total);

        // Copy grid values
        for (int cube = 0; cube < config_.grid.n_chan; cube++) {
            for (int pol = 0; pol < config_.grid.n_pol; pol++) {
                for (int y = 0; y < config_.grid.ny; y++) {
                    for (int x = 0; x < config_.grid.nx; x++) {
                        size_t idx = ((cube * config_.grid.n_pol + pol) * config_.grid.ny + y) *
                                     config_.grid.nx + x;
                        image_[idx] = grid(x, y, pol, cube);
                    }
                }
            }
        }

        // Apply gridding correction
        apply_gridding_correction(
            image_.data(),
            config_.grid.nx,
            config_.grid.ny,
            config_.grid.cf_support
        );

        // Write output
        if (!config_.imaging.output_fits.empty()) {
            auto [ra, dec] = reader_->get_phase_center();
            double ra_deg = ra * 180.0 / M_PI;
            double dec_deg = dec * 180.0 / M_PI;
            double cell_deg = config_.grid.cell_size_asec / 3600.0;

            write_fits_image(
                config_.imaging.output_fits,
                image_.data(),
                config_.grid.nx, config_.grid.ny,
                config_.grid.n_pol, config_.grid.n_chan,
                ra_deg, dec_deg, cell_deg,
                ref_freq, freq.channel_width_hz
            );
        }

        auto end_time = Clock::now();
        double total_time = std::chrono::duration<double>(end_time - start_time).count();

        std::cout << "\nTiming summary:" << std::endl;
        std::cout << "  Read + Process: " << stats_.process_time_sec << " sec" << std::endl;
        std::cout << "  Gridding: " << stats_.grid_time_sec << " sec" << std::endl;
        std::cout << "  FFT: " << stats_.fft_time_sec << " sec" << std::endl;
        std::cout << "  Total: " << total_time << " sec" << std::endl;

        return true;
    }
};

//-----------------------------------------------------------------------------
// Public interface
//-----------------------------------------------------------------------------

UGMRTImager::UGMRTImager(const ImagerConfig& config)
    : impl_(std::make_unique<Impl>(config))
{
}

UGMRTImager::~UGMRTImager() = default;

bool UGMRTImager::initialize() {
    return impl_->initialize();
}

bool UGMRTImager::run() {
    return impl_->run();
}

const std::complex<double>* UGMRTImager::get_image() const {
    return impl_->image_.data();
}

std::array<int, 4> UGMRTImager::get_image_shape() const {
    return {
        impl_->config_.grid.nx,
        impl_->config_.grid.ny,
        impl_->config_.grid.n_pol,
        impl_->config_.grid.n_chan
    };
}

bool UGMRTImager::write_fits(const std::string& filename) {
    auto [ra, dec] = impl_->reader_->get_phase_center();
    double ra_deg = ra * 180.0 / M_PI;
    double dec_deg = dec * 180.0 / M_PI;
    double cell_deg = impl_->config_.grid.cell_size_asec / 3600.0;
    const auto& freq = impl_->reader_->get_freq_setup();

    return write_fits_image(
        filename,
        impl_->image_.data(),
        impl_->config_.grid.nx,
        impl_->config_.grid.ny,
        impl_->config_.grid.n_pol,
        impl_->config_.grid.n_chan,
        ra_deg, dec_deg, cell_deg,
        freq.center_freq_hz(),
        freq.channel_width_hz
    );
}

UGMRTImager::Stats UGMRTImager::get_stats() const {
    return impl_->stats_;
}

} // namespace ugmrt
