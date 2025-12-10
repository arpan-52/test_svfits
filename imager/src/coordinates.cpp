#include "coordinates.hpp"
#include "utils.hpp"
#include <cmath>
#include <iostream>

Coordinates::Coordinates(const Config& config)
    : config_(config), use_j2000_(config.epoch > 0), reference_mjd_(0)
{
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            R_app_to_j2000_[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }
}

void Coordinates::initialize(double middle_mjd) {
    reference_mjd_ = middle_mjd;

    if (use_j2000_) {
        std::cout << "Computing J2000 rotation matrix for MJD "
                  << middle_mjd << std::endl;
        computeRotationMatrix(middle_mjd);
        std::cout << "Single rotation matrix computed (will be used for all visibilities)"
                  << std::endl;
    } else {
        std::cout << "Using apparent coordinates (no J2000 conversion)" << std::endl;
    }
}

void Coordinates::computeRotationMatrix(double mjd) {
    // For now, use identity matrix
    // TODO: Implement NOVAS calls for proper precession/nutation
    std::cout << "INFO: Using identity matrix - full NOVAS integration pending" << std::endl;
}

UVWFrame Coordinates::computeUVW(double mjd,
                                  const std::vector<AntennaInfo>& antennas)
{
    double lst = mjdToLST(mjd);
    double ha = lst - config_.ra_app;
    double dec = config_.dec_app;

    double sin_ha = std::sin(ha);
    double cos_ha = std::cos(ha);
    double sin_dec = std::sin(dec);
    double cos_dec = std::cos(dec);

    UVWFrame uvw;
    uvw.u.resize(antennas.size());
    uvw.v.resize(antennas.size());
    uvw.w.resize(antennas.size());

    for (size_t i = 0; i < antennas.size(); i++) {
        double bx = antennas[i].bx;
        double by = antennas[i].by;
        double bz = antennas[i].bz;

        double u_app = (bx * sin_ha + by * cos_ha) / C_LIGHT;
        double v_app = (bz * cos_dec - (bx*cos_ha - by*sin_ha) * sin_dec) / C_LIGHT;
        double w_app = ((bx*cos_ha - by*sin_ha) * cos_dec + bz * sin_dec) / C_LIGHT;

        if (use_j2000_) {
            applyRotation(u_app, v_app, w_app, uvw.u[i], uvw.v[i], uvw.w[i]);
        } else {
            uvw.u[i] = u_app;
            uvw.v[i] = v_app;
            uvw.w[i] = w_app;
        }
    }

    return uvw;
}

void Coordinates::applyRotation(double u_app, double v_app, double w_app,
                                 double& u_j2000, double& v_j2000, double& w_j2000)
{
    u_j2000 = R_app_to_j2000_[0][0] * u_app +
              R_app_to_j2000_[0][1] * v_app +
              R_app_to_j2000_[0][2] * w_app;

    v_j2000 = R_app_to_j2000_[1][0] * u_app +
              R_app_to_j2000_[1][1] * v_app +
              R_app_to_j2000_[1][2] * w_app;

    w_j2000 = R_app_to_j2000_[2][0] * u_app +
              R_app_to_j2000_[2][1] * v_app +
              R_app_to_j2000_[2][2] * w_app;
}

double Coordinates::mjdToLST(double mjd) const {
    double jd = mjd + 2400000.5;
    double t = (jd - 2451545.0) / 36525.0;

    double gmst = 280.46061837 + 360.98564736629 * (jd - 2451545.0)
                + 0.000387933 * t * t - t * t * t / 38710000.0;

    gmst = std::fmod(gmst, 360.0);
    if (gmst < 0) gmst += 360.0;

    double lst = (gmst * PI / 180.0) + GMRT_LON;
    lst = std::fmod(lst, 2.0 * PI);
    if (lst < 0) lst += 2.0 * PI;

    return lst;
}

double Coordinates::lstToHA(double lst, double ra) const {
    return lst - ra;
}
