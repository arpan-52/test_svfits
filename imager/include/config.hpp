#pragma once
#include <string>
#include <vector>
#include <cstdint>

enum class RefFreqMode {
    BURST_FREQ,
    HIGHEST_FREQ,
    MANUAL
};

class Config {
public:
    // I/O
    int nfiles = 16;
    std::string path;
    std::vector<std::string> input_files;
    std::string fits_output = "burst_image.fits";
    std::string antsamp_file = "antsamp.hdr";

    // Frequency
    double freq_start = 550e6;
    double freq_end = 750e6;
    int n_channels = 4096;
    double ref_freq = 650e6;
    RefFreqMode ref_freq_mode = RefFreqMode::BURST_FREQ;

    // Burst
    std::string burst_name = "FRB";
    double burst_mjd = 0.0;
    double burst_dm = 0.0;
    double burst_intwd = 0.001;
    double burst_freq = 750e6;
    int burst_beam_id = 0;
    double burst_ra_app = 0.0;
    double burst_dec_app = 0.0;
    bool update_burst = true;

    // Coordinates
    double ra_mean = 0.0;
    double dec_mean = 0.0;
    double ra_app = 0.0;
    double dec_app = 0.0;
    double epoch = 2000.0;

    // Antenna mask
    uint32_t antmask = 0x3FFFFFFF;  // All 30 antennas

    // Processing
    bool do_bandpass = true;
    bool do_baseline = true;
    bool do_flag = true;
    float flag_threshold = 5.0;
    int num_threads = 8;

    // Imaging
    int grid_nx = 512;
    int grid_ny = 512;
    double cell_size_arcsec = 1.0;
    int n_w_planes = 64;
    int cf_support = 7;
    int cf_oversampling = 128;

    // CLEAN
    float clean_gain = 0.1f;
    int clean_niter = 1000;
    float clean_threshold = 0.001f;

    // IAT-UTC
    float iatutc = 37.0f;

    // Methods
    static Config fromFile(const std::string& filename);
    void validate();
    void computeDerived();

private:
    void parseLine(const std::string& line);
    void setDefaults();
};
