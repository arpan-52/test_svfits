#include "reader.hpp"
#include "utils.hpp"
#include <fstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <sys/time.h>

Reader::Reader(const Config& config)
    : config_(config), coordinates_(config), middle_time_(0), max_w_(0)
{
    files_.resize(config.nfiles);
}

void Reader::initialize() {
    std::cout << "=== Reader Initialization ===" << std::endl;

    loadAntsampHeader(config_.antsamp_file);
    buildBaselines();
    computeBurstExtent();

    // Compute middle time
    double t_min = 1e100, t_max = -1e100;
    for (const auto& f : files_) {
        if (f.n_rec > 0) {
            t_min = std::min(t_min, f.t_start + f.start_rec * 0.00131072);
            t_max = std::max(t_max, f.t_start + (f.start_rec + f.n_rec) * 0.00131072);
        }
    }
    middle_time_ = config_.burst_mjd;  // Use burst MJD as middle time

    std::cout << "Middle time (MJD): " << middle_time_ << std::endl;

    // Initialize coordinates with middle time
    coordinates_.initialize(middle_time_);

    // Compute max W
    computeMaxW();

    std::cout << "=== Reader Initialized ===" << std::endl;
}

void Reader::loadAntsampHeader(const std::string& filename) {
    std::cout << "Loading antenna configuration from " << filename << std::endl;

    std::ifstream file(filename);
    if (!file) {
        throw std::runtime_error("Cannot open antsamp.hdr: " + filename);
    }

    std::string line;
    bool in_antenna_block = false;

    while (std::getline(file, line)) {
        if (line.find("Antenna.def") != std::string::npos) {
            in_antenna_block = true;
            continue;
        }
        if (line.find("} Antenna") != std::string::npos) {
            in_antenna_block = false;
            continue;
        }

        if (in_antenna_block && line.find("ANT") == 0) {
            auto parts = split(line, ' ');
            if (parts.size() >= 5) {
                AntennaInfo ant;
                ant.name = parts[2];  // e.g., "C00"
                ant.bx = std::stod(parts[3]);
                ant.by = std::stod(parts[4]);
                ant.bz = std::stod(parts[5]);
                ant.is_lsb = false;
                antennas_.push_back(ant);
            }
        }
    }

    std::cout << "Loaded " << antennas_.size() << " antennas" << std::endl;
}

void Reader::buildBaselines() {
    int n_ants = antennas_.size();
    int bl_idx = 0;

    for (int i = 0; i < n_ants; i++) {
        for (int j = i+1; j < n_ants; j++) {
            // Check antenna mask
            if (!((1u << i) & config_.antmask)) continue;
            if (!((1u << j) & config_.antmask)) continue;

            BaselineInfo bl;
            bl.ant0 = i;
            bl.ant1 = j;
            bl.band0 = 0;
            bl.band1 = 0;
            bl.drop = false;
            bl.flip_conjugate = false;  // Handle USB/LSB
            bl.index = bl_idx++;
            baselines_.push_back(bl);
        }
    }

    std::cout << "Built " << baselines_.size() << " baselines" << std::endl;
}

void Reader::computeMaxW() {
    double max_baseline = 0.0;

    for (size_t i = 0; i < antennas_.size(); i++) {
        for (size_t j = i+1; j < antennas_.size(); j++) {
            double dx = antennas_[i].bx - antennas_[j].bx;
            double dy = antennas_[i].by - antennas_[j].by;
            double dz = antennas_[i].bz - antennas_[j].bz;
            double baseline = std::sqrt(dx*dx + dy*dy + dz*dz);
            max_baseline = std::max(max_baseline, baseline);
        }
    }

    double wavelength = C_LIGHT / config_.ref_freq;
    max_w_ = max_baseline * std::fabs(std::sin(config_.dec_app)) / wavelength;

    std::cout << "Max baseline: " << max_baseline << " m" << std::endl;
    std::cout << "Max W: " << max_w_ << " wavelengths" << std::endl;
}

void Reader::computeBurstExtent() {
    // Simplified: assume all files have burst signal
    for (int i = 0; i < config_.nfiles; i++) {
        files_[i].filename = config_.path + "/" + config_.input_files[i];
        files_[i].file_index = i;
        files_[i].t_start = 0.0;
        files_[i].mjd_ref = config_.burst_mjd;
        files_[i].start_rec = 0;
        files_[i].n_rec = REC_PER_SLICE;  // Process one slice for now
    }

    burst_channels_.resize(config_.nfiles);
    for (int i = 0; i < config_.nfiles; i++) {
        identifyBurstChannels(i);
    }
}

void Reader::identifyBurstChannels(int file_idx) {
    burst_channels_[file_idx].resize(REC_PER_SLICE);

    // Simplified: assume burst covers all channels
    for (int rec = 0; rec < REC_PER_SLICE; rec++) {
        burst_channels_[file_idx][rec].channel_start = 0;
        burst_channels_[file_idx][rec].channel_end = config_.n_channels - 1;
    }
}

int Reader::selectWPlane(float w) const {
    if (config_.n_w_planes == 1) return 0;

    float w_norm = (w + max_w_) / (2.0f * max_w_);
    w_norm = std::clamp(w_norm, 0.0f, 1.0f);

    return static_cast<int>(std::round(w_norm * (config_.n_w_planes - 1)));
}

std::vector<Visibility> Reader::readAllVisibilities() {
    std::cout << "=== Reading Visibilities ===" << std::endl;

    std::vector<Visibility> all_vis;

    // Simplified: read first slice of first file only
    std::cout << "Reading file 0, slice 0 (simplified demo)" << std::endl;

    auto vis = extractVisibilities(0, 0);
    all_vis.insert(all_vis.end(), vis.begin(), vis.end());

    std::cout << "Total visibilities: " << all_vis.size() << std::endl;

    return all_vis;
}

std::vector<Visibility> Reader::extractVisibilities(int file_idx, int slice_idx) {
    std::vector<Visibility> vis_list;

    // Generate synthetic visibilities for demo
    double t_mid = middle_time_;

    UVWFrame uvw = coordinates_.computeUVW(t_mid, antennas_);

    double freq = config_.ref_freq;
    double wavelength = C_LIGHT / freq;

    // Create visibilities for each baseline
    for (const auto& bl : baselines_) {
        if (bl.drop) continue;

        float u_lambda = (uvw.u[bl.ant1] - uvw.u[bl.ant0]) * freq / C_LIGHT;
        float v_lambda = (uvw.v[bl.ant1] - uvw.v[bl.ant0]) * freq / C_LIGHT;
        float w_lambda = (uvw.w[bl.ant1] - uvw.w[bl.ant0]) * freq / C_LIGHT;

        Visibility vis;
        vis.re = 1.0f;  // Synthetic data
        vis.im = 0.0f;
        vis.weight = 1.0f;
        vis.u = u_lambda;
        vis.v = v_lambda;
        vis.w = w_lambda;
        vis.freq = freq;
        vis.time = t_mid;
        vis.ant0 = bl.ant0;
        vis.ant1 = bl.ant1;
        vis.channel = 0;
        vis.w_plane = selectWPlane(vis.w);

        vis_list.push_back(vis);
    }

    std::cout << "Extracted " << vis_list.size() << " visibilities" << std::endl;

    return vis_list;
}

bool Reader::readSlice(int file_idx, int slice_idx, std::vector<uint16_t>& buffer) {
    // TODO: Implement actual file reading
    return false;
}

double Reader::getSliceTime(int file_idx, int slice_idx) {
    return files_[file_idx].t_start + slice_idx * 0.65536;
}

void Reader::computeBandpass(const std::vector<uint16_t>& raw_data,
                            const std::vector<BurstChannel>& burst_ch,
                            std::vector<std::vector<float>>& bpass) {
    // TODO: Implement bandpass computation
}

void Reader::computeBaseline(const std::vector<uint16_t>& raw_data,
                            const std::vector<BurstChannel>& burst_ch,
                            std::vector<std::vector<std::complex<float>>>& baseline) {
    // TODO: Implement baseline computation
}

void Reader::flagRFI(std::vector<Visibility>& vis) {
    // TODO: Implement RFI flagging
}
