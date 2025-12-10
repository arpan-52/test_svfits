#pragma once
#include "config.hpp"
#include <vector>

struct AntennaInfo {
    std::string name;
    double bx, by, bz;  // ITRF coordinates (meters)
    int sampler_usb130;
    int sampler_usb175;
    bool is_lsb;
};

struct UVWFrame {
    std::vector<double> u, v, w;  // Per-antenna UVW (seconds)
};

class Coordinates {
public:
    explicit Coordinates(const Config& config);

    // Initialize with middle time of observation
    void initialize(double middle_mjd);

    // Compute UVW (uses pre-computed rotation matrix if J2000)
    UVWFrame computeUVW(double mjd, const std::vector<AntennaInfo>& antennas);

    // Time conversions
    double mjdToLST(double mjd) const;
    double lstToHA(double lst, double ra) const;

private:
    Config config_;

    // Pre-computed rotation matrix (set during initialize())
    double R_app_to_j2000_[3][3];
    bool use_j2000_;
    double reference_mjd_;

    // GMRT location
    static constexpr double GMRT_LON = 1.306056;  // radians (74.05° E)
    static constexpr double GMRT_LAT = 0.335567;  // radians (19.23° N)

    // Compute rotation matrices using NOVAS
    void computeRotationMatrix(double mjd);
    void applyRotation(double u_app, double v_app, double w_app,
                      double& u_j2000, double& v_j2000, double& w_j2000);
};
