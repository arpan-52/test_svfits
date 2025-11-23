/**
 * @file svfits_reader.cpp
 * @brief Implementation of svfits reader wrapper
 */

#include "svfits_reader.hpp"
#include "ugmrt_imager.hpp"

#include <cstring>
#include <cmath>
#include <stdexcept>
#include <iostream>

// Include svfits C headers
extern "C" {
#include "svfits/svio.h"
#include "svfits/gmrt_newcorr.h"

// Declare the C functions we need
int init_user(SvSelectionType *user, char *fname, char *anthdr,
              char *bhdrfile, char *bulletinA);
int read_slice(SvSelectionType *user, int idx, int slice, char *rbuf);
int get_file_order(SvSelectionType *user, int *order);
void init_mat(SvSelectionType *user, double tm);
float half_to_float(const unsigned short x);
double get_ha(SvSelectionType *user, double tm);
double lmst(double mjd);

#ifdef USE_NOVAS
int novas_prenut_vis(SvSelectionType *user, UvwParType *uvw, double mjd);
#else
void sla_prenut_vis(UvwParType *uvw, double mjd, double ra_app, double dec_app,
                    double epoch1);
#endif

int robust_stats(int n, float *x, float *med, float *mad);
}

namespace ugmrt {

//-----------------------------------------------------------------------------
// Implementation class
//-----------------------------------------------------------------------------

class SvfitsReader::Impl {
public:
    // Configuration
    std::string param_file_;
    std::string antsamp_file_;
    std::string bulletin_a_;

    // svfits structures
    SvSelectionType user_;
    CorrType corr_;
    ScanRecType srec_;
    ScanInfoType scan_;
    InitHdrType hdr_;

    // Derived data
    BurstParams burst_params_;
    FreqSetup freq_setup_;
    std::vector<AntennaInfo> antennas_;
    std::vector<BaselineInfo> baselines_;

    // Buffers
    std::vector<char> raw_buffer_;
    std::vector<std::vector<char>> file_buffers_;

    // Bandpass data
    std::vector<std::vector<float>> bandpass_;      // [baseline][channel]
    std::vector<std::vector<std::complex<float>>> off_source_; // [baseline][channel]

    bool initialized_ = false;

    Impl(const std::string& param_file,
         const std::string& antsamp_file,
         const std::string& bulletin_a)
        : param_file_(param_file)
        , antsamp_file_(antsamp_file)
        , bulletin_a_(bulletin_a)
    {
        // Initialize structures to zero
        memset(&user_, 0, sizeof(user_));
        memset(&corr_, 0, sizeof(corr_));
        memset(&srec_, 0, sizeof(srec_));
        memset(&scan_, 0, sizeof(scan_));
        memset(&hdr_, 0, sizeof(hdr_));

        // Link structures
        user_.corr = &corr_;
        user_.srec = &srec_;
        user_.hdr = &hdr_;
        srec_.scan = &scan_;
        srec_.corr = &corr_;
    }

    bool initialize() {
        // Call svfits init_user
        char param_file_buf[PATHLEN];
        char antsamp_buf[PATHLEN];
        char bulletin_buf[PATHLEN];

        strncpy(param_file_buf, param_file_.c_str(), PATHLEN - 1);
        strncpy(antsamp_buf, antsamp_file_.c_str(), PATHLEN - 1);

        char* bulletin_ptr = nullptr;
        if (!bulletin_a_.empty()) {
            strncpy(bulletin_buf, bulletin_a_.c_str(), PATHLEN - 1);
            bulletin_ptr = bulletin_buf;
        }

        if (init_user(&user_, param_file_buf, antsamp_buf, nullptr, bulletin_ptr) != 0) {
            std::cerr << "Failed to initialize svfits user structure" << std::endl;
            return false;
        }

        // Extract burst parameters
        burst_params_.name = user_.burst.name;
        burst_params_.mjd = user_.burst.mjd;
        burst_params_.time_ist = user_.burst.t;
        burst_params_.dm = user_.burst.DM;
        burst_params_.intrinsic_width = user_.burst.int_wd;
        burst_params_.ref_freq_hz = user_.burst.f;
        burst_params_.ra_app = user_.burst.ra_app;
        burst_params_.dec_app = user_.burst.dec_app;

        // Extract frequency setup
        freq_setup_.freq_start_hz = scan_.source.freq[0];
        freq_setup_.channel_width_hz = scan_.source.ch_width;
        freq_setup_.n_channels = corr_.daspar.channels;
        freq_setup_.freq_end_hz = freq_setup_.freq_start_hz +
            freq_setup_.n_channels * freq_setup_.channel_width_hz;

        // Extract antenna info
        for (int i = 0; i < MAX_ANTS; i++) {
            if ((user_.antmask >> i) & 1) {
                AntennaInfo ant;
                ant.id = i;
                ant.name = corr_.antenna[i].name;
                ant.x = corr_.antenna[i].bx;
                ant.y = corr_.antenna[i].by;
                ant.z = corr_.antenna[i].bz;
                antennas_.push_back(ant);
            }
        }

        // Extract baseline info from vispar
        for (int bl = 0; bl < user_.baselines; bl++) {
            BaselineInfo binfo;
            binfo.ant0 = user_.vispar.visinfo[bl].ant0;
            binfo.ant1 = user_.vispar.visinfo[bl].ant1;
            binfo.band0 = user_.vispar.visinfo[bl].band0;
            binfo.band1 = user_.vispar.visinfo[bl].band1;
            binfo.flip_sign = user_.vispar.visinfo[bl].flip != 0;
            baselines_.push_back(binfo);
        }

        // Allocate raw buffer
        size_t recl = corr_.daspar.baselines * corr_.daspar.channels * sizeof(float);
        raw_buffer_.resize(MAX_REC_PER_SLICE * recl);

        // Allocate bandpass arrays
        bandpass_.resize(corr_.daspar.baselines);
        off_source_.resize(corr_.daspar.baselines);
        for (int bl = 0; bl < corr_.daspar.baselines; bl++) {
            bandpass_[bl].resize(corr_.daspar.channels, 1.0f);
            off_source_[bl].resize(corr_.daspar.channels, {0.0f, 0.0f});
        }

        initialized_ = true;
        return true;
    }

    /**
     * @brief Compute bandpass from off-burst records
     */
    void compute_bandpass(int file_idx, int slice, char* raw_buf) {
        int channels = corr_.daspar.channels;
        int baselines = corr_.daspar.baselines;
        size_t recl = baselines * channels * sizeof(float);

        // Get burst channel range for each record
        short* start_chan = user_.bpass.start_chan;
        short* end_chan = user_.bpass.end_chan;

        // For each baseline, compute mean amplitude over off-burst records
        for (int bl = 0; bl < baselines; bl++) {
            std::vector<float> sum(channels, 0.0f);
            std::vector<int> count(channels, 0);

            for (int rec = 0; rec < MAX_REC_PER_SLICE; rec++) {
                char* rec_ptr = raw_buf + rec * recl;
                unsigned short* vis_ptr = (unsigned short*)(rec_ptr +
                    user_.vispar.visinfo[bl].off);

                for (int ch = 0; ch < channels; ch++) {
                    // Skip burst region
                    if (ch >= start_chan[rec] && ch <= end_chan[rec]) {
                        continue;
                    }

                    float re = half_to_float(vis_ptr[ch * 2]);
                    float im = half_to_float(vis_ptr[ch * 2 + 1]);
                    float amp = std::sqrt(re * re + im * im);

                    if (std::isfinite(amp) && amp > 0) {
                        sum[ch] += amp;
                        count[ch]++;
                    }
                }
            }

            // Compute normalized bandpass
            for (int ch = 0; ch < channels; ch++) {
                if (count[ch] > 0) {
                    bandpass_[bl][ch] = sum[ch] / count[ch];
                } else {
                    bandpass_[bl][ch] = 1.0f;
                }
            }
        }
    }

    /**
     * @brief Compute off-source mean for baseline subtraction
     */
    void compute_off_source(int file_idx, int slice, char* raw_buf) {
        int channels = corr_.daspar.channels;
        int baselines = corr_.daspar.baselines;
        size_t recl = baselines * channels * sizeof(float);

        short* start_chan = user_.bpass.start_chan;
        short* end_chan = user_.bpass.end_chan;

        for (int bl = 0; bl < baselines; bl++) {
            std::vector<float> sum_re(channels, 0.0f);
            std::vector<float> sum_im(channels, 0.0f);
            std::vector<int> count(channels, 0);

            for (int rec = 0; rec < MAX_REC_PER_SLICE; rec++) {
                char* rec_ptr = raw_buf + rec * recl;
                unsigned short* vis_ptr = (unsigned short*)(rec_ptr +
                    user_.vispar.visinfo[bl].off);

                for (int ch = 0; ch < channels; ch++) {
                    if (ch >= start_chan[rec] && ch <= end_chan[rec]) {
                        continue;
                    }

                    float re = half_to_float(vis_ptr[ch * 2]);
                    float im = half_to_float(vis_ptr[ch * 2 + 1]);

                    if (std::isfinite(re) && std::isfinite(im)) {
                        sum_re[ch] += re;
                        sum_im[ch] += im;
                        count[ch]++;
                    }
                }
            }

            for (int ch = 0; ch < channels; ch++) {
                if (count[ch] > 0) {
                    off_source_[bl][ch] = {
                        sum_re[ch] / count[ch],
                        sum_im[ch] / count[ch]
                    };
                } else {
                    off_source_[bl][ch] = {0.0f, 0.0f};
                }
            }
        }
    }

    /**
     * @brief Compute UVW coordinates for a baseline at given time
     */
    void compute_uvw(int bl, double mjd, double* u, double* v, double* w) {
        // Get antenna positions
        int ant0 = baselines_[bl].ant0;
        int ant1 = baselines_[bl].ant1;

        double bx = corr_.antenna[ant1].bx - corr_.antenna[ant0].bx;
        double by = corr_.antenna[ant1].by - corr_.antenna[ant0].by;
        double bz = corr_.antenna[ant1].bz - corr_.antenna[ant0].bz;

        // Get hour angle
        double ha = get_ha(&user_, mjd);

        // Phase center
        double ra = user_.burst.ra_app;
        double dec = user_.burst.dec_app;

        // Convert baseline to UVW
        double sin_ha = std::sin(ha);
        double cos_ha = std::cos(ha);
        double sin_dec = std::sin(dec);
        double cos_dec = std::cos(dec);

        // UVW rotation (equatorial to UVW)
        *u = sin_ha * bx + cos_ha * by;
        *v = -sin_dec * cos_ha * bx + sin_dec * sin_ha * by + cos_dec * bz;
        *w = cos_dec * cos_ha * bx - cos_dec * sin_ha * by + sin_dec * bz;

        // Convert to wavelengths
        double wavelength = 299792458.0 / freq_setup_.center_freq_hz();
        *u /= wavelength;
        *v /= wavelength;
        *w /= wavelength;
    }

    /**
     * @brief Process all visibilities
     */
    size_t process(VisibilityCallback callback, const ReaderOptions& options) {
        if (!initialized_) {
            throw std::runtime_error("SvfitsReader not initialized");
        }

        size_t count = 0;
        int channels = corr_.daspar.channels;
        int baselines = user_.baselines;
        int nfiles = user_.recfile.nfiles;
        size_t recl = corr_.daspar.baselines * channels * sizeof(float);

        // Get file processing order
        int file_order[MaxRecFiles];
        get_file_order(&user_, file_order);

        // Process each file
        for (int f = 0; f < nfiles; f++) {
            int file_idx = file_order[f];

            // Skip files with no burst data
            if (user_.recfile.n_rec[file_idx] <= 0) {
                continue;
            }

            // Read raw data
            if (read_slice(&user_, file_idx, 0, raw_buffer_.data()) != 0) {
                std::cerr << "Failed to read file " << file_idx << std::endl;
                continue;
            }

            // Compute bandpass and off-source if needed
            if (options.do_bandpass) {
                compute_bandpass(file_idx, 0, raw_buffer_.data());
            }
            if (options.do_baseline) {
                compute_off_source(file_idx, 0, raw_buffer_.data());
            }

            // Process records with burst signal
            int start_rec = user_.recfile.start_rec[file_idx];
            int n_rec = user_.recfile.n_rec[file_idx];

            for (int rec = start_rec; rec < start_rec + n_rec; rec++) {
                char* rec_ptr = raw_buffer_.data() + rec * recl;

                // Get timestamp for this record
                double t_rec = user_.recfile.t_start[file_idx] +
                    rec * user_.recfile.t_slice / MAX_REC_PER_SLICE;
                double mjd = user_.recfile.mjd_ref + t_rec / 86400.0;

                // Initialize coordinate matrices for this time
                init_mat(&user_, mjd);

                // Get channel range for this record
                int start_ch = user_.bpass.start_chan[rec];
                int end_ch = user_.bpass.end_chan[rec];

                // Process each baseline
                for (int bl = 0; bl < baselines; bl++) {
                    unsigned short* vis_ptr = (unsigned short*)(rec_ptr +
                        user_.vispar.visinfo[bl].off);

                    // Compute UVW
                    double u, v, w;
                    compute_uvw(bl, mjd, &u, &v, &w);

                    // Collect visibilities for flagging
                    std::vector<float> vis_amp;

                    // Process channels in burst region
                    for (int ch = start_ch; ch <= end_ch; ch++) {
                        float re = half_to_float(vis_ptr[ch * 2]);
                        float im = half_to_float(vis_ptr[ch * 2 + 1]);

                        // Skip invalid data
                        if (!std::isfinite(re) || !std::isfinite(im)) {
                            continue;
                        }

                        // Flip imaginary if needed
                        if (baselines_[bl].flip_sign) {
                            im = -im;
                        }

                        // Apply bandpass correction
                        if (options.do_bandpass && bandpass_[bl][ch] > 0) {
                            re /= bandpass_[bl][ch];
                            im /= bandpass_[bl][ch];
                        }

                        // Apply baseline subtraction
                        if (options.do_baseline) {
                            re -= off_source_[bl][ch].real();
                            im -= off_source_[bl][ch].imag();
                        }

                        // Compute weight (for flagging)
                        float weight = 1.0f;
                        if (options.do_flag) {
                            float amp = std::sqrt(re * re + im * im);
                            vis_amp.push_back(amp);
                        }

                        // Create processed visibility
                        ProcessedVisibility pvis;
                        pvis.vis[0] = {re, im};
                        pvis.vis[1] = {re, im};  // Same for both pols for now
                        pvis.weight[0] = weight;
                        pvis.weight[1] = weight;
                        pvis.u = u;
                        pvis.v = v;
                        pvis.w = w;
                        pvis.freq_hz = freq_setup_.freq_start_hz +
                            ch * freq_setup_.channel_width_hz;
                        pvis.channel = ch;
                        pvis.baseline = bl;
                        pvis.time_mjd = mjd;

                        // Apply flagging if enabled
                        if (options.do_flag && !vis_amp.empty()) {
                            float med, mad;
                            if (robust_stats(vis_amp.size(), vis_amp.data(), &med, &mad) == 0) {
                                float amp = std::sqrt(re * re + im * im);
                                if (std::abs(amp - med) > options.flag_threshold * mad) {
                                    pvis.weight[0] = -1.0f;
                                    pvis.weight[1] = -1.0f;
                                }
                            }
                        }

                        // Call callback
                        if (!callback(pvis)) {
                            return count;
                        }
                        count++;
                    }
                }
            }
        }

        return count;
    }
};

//-----------------------------------------------------------------------------
// Public interface implementation
//-----------------------------------------------------------------------------

SvfitsReader::SvfitsReader(
    const std::string& param_file,
    const std::string& antsamp_file,
    const std::string& bulletin_a)
    : impl_(std::make_unique<Impl>(param_file, antsamp_file, bulletin_a))
{
}

SvfitsReader::~SvfitsReader() = default;

bool SvfitsReader::initialize() {
    return impl_->initialize();
}

const BurstParams& SvfitsReader::get_burst_params() const {
    return impl_->burst_params_;
}

const FreqSetup& SvfitsReader::get_freq_setup() const {
    return impl_->freq_setup_;
}

const std::vector<AntennaInfo>& SvfitsReader::get_antennas() const {
    return impl_->antennas_;
}

const std::vector<BaselineInfo>& SvfitsReader::get_baselines() const {
    return impl_->baselines_;
}

int SvfitsReader::get_num_files() const {
    return impl_->user_.recfile.nfiles;
}

int SvfitsReader::get_burst_records() const {
    int total = 0;
    for (int f = 0; f < impl_->user_.recfile.nfiles; f++) {
        total += impl_->user_.recfile.n_rec[f];
    }
    return total;
}

size_t SvfitsReader::process(
    VisibilityCallback callback,
    const ReaderOptions& options)
{
    return impl_->process(callback, options);
}

std::vector<ProcessedVisibility> SvfitsReader::get_all_visibilities(
    const ReaderOptions& options)
{
    std::vector<ProcessedVisibility> result;
    process([&result](const ProcessedVisibility& vis) {
        result.push_back(vis);
        return true;
    }, options);
    return result;
}

double SvfitsReader::get_ref_wavelength() const {
    return 299792458.0 / impl_->freq_setup_.center_freq_hz();
}

std::pair<double, double> SvfitsReader::get_phase_center() const {
    // Return J2000 coordinates
    return {impl_->burst_params_.ra_app, impl_->burst_params_.dec_app};
}

} // namespace ugmrt
