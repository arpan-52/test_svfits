#pragma once
#include "config.hpp"
#include "coordinates.hpp"
#include <vector>
#include <string>
#include <cstdint>

struct Visibility {
    float re, im;
    float weight;
    float u, v, w;      // UVW in wavelengths at ref_freq
    double freq;
    double time;
    int ant0, ant1;
    int channel;
    int w_plane;        // Pre-computed W-plane index
};

struct BaselineInfo {
    int ant0, ant1;
    int band0, band1;
    bool drop;
    bool flip_conjugate;
    int index;
};

struct BurstChannel {
    int channel_start;
    int channel_end;
};

struct RawFileInfo {
    std::string filename;
    int file_index;
    double t_start;
    double mjd_ref;
    int start_rec;
    int n_rec;
};

class Reader {
public:
    Reader(const Config& config);

    void initialize();
    std::vector<Visibility> readAllVisibilities();

    const std::vector<AntennaInfo>& getAntennas() const { return antennas_; }
    double getMiddleTime() const { return middle_time_; }
    float getMaxW() const { return max_w_; }

private:
    Config config_;
    Coordinates coordinates_;

    std::vector<AntennaInfo> antennas_;
    std::vector<BaselineInfo> baselines_;
    std::vector<RawFileInfo> files_;
    std::vector<std::vector<BurstChannel>> burst_channels_;  // [file][record]

    double middle_time_;  // Middle MJD of observation
    float max_w_;

    static constexpr int REC_PER_SLICE = 50;

    // Initialization
    void loadAntsampHeader(const std::string& filename);
    void buildBaselines();
    void computeMaxW();
    void computeBurstExtent();
    void identifyBurstChannels(int file_idx);

    // Data reading
    bool readSlice(int file_idx, int slice_idx, std::vector<uint16_t>& buffer);
    double getSliceTime(int file_idx, int slice_idx);

    // Calibration
    void computeBandpass(const std::vector<uint16_t>& raw_data,
                        const std::vector<BurstChannel>& burst_ch,
                        std::vector<std::vector<float>>& bpass);
    void computeBaseline(const std::vector<uint16_t>& raw_data,
                        const std::vector<BurstChannel>& burst_ch,
                        std::vector<std::vector<std::complex<float>>>& baseline);
    void flagRFI(std::vector<Visibility>& vis);

    // Visibility extraction
    std::vector<Visibility> extractVisibilities(int file_idx, int slice_idx);

    // W-plane selection
    int selectWPlane(float w) const;
};
