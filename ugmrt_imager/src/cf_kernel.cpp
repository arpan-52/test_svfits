/**
 * @file cf_kernel.cpp
 * @brief Prolate Spheroidal Wave Function convolution kernel
 */

#include "cf_kernel.hpp"
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace ugmrt {

//-----------------------------------------------------------------------------
// Spheroidal function approximation (from AIPS/CASA)
//-----------------------------------------------------------------------------

/**
 * @brief Compute spheroidal function value
 *
 * This is the standard approximation used in radio astronomy gridding.
 * Based on Schwab, F.R. 1984, AJ 89, 1076
 */
double spheroidal(double nu) {
    static const double P[] = {
        0.08203343, -0.3644705, 0.627866, -0.5335581, 0.2312756,
        0.004028559, -0.03697768, 0.1021332, -0.1201436, 0.06412774
    };
    static const double Q[] = {
        1.0, 0.8212018, 0.2078043,
        1.0, 0.9599102, 0.2918724
    };

    double nuend = 1.0;
    if (std::abs(nu) > nuend) {
        return 0.0;
    }

    double nusq = nu * nu;
    int idx = (nusq <= 0.75) ? 0 : 1;

    double top, bot;
    if (idx == 0) {
        double delnusq = nusq - 0.75;
        top = P[0] + delnusq * (P[1] + delnusq * (P[2] + delnusq * (P[3] + delnusq * P[4])));
        bot = Q[0] + delnusq * (Q[1] + delnusq * Q[2]);
    } else {
        double delnusq = nusq - 0.75;
        top = P[5] + delnusq * (P[6] + delnusq * (P[7] + delnusq * (P[8] + delnusq * P[9])));
        bot = Q[3] + delnusq * (Q[4] + delnusq * Q[5]);
    }

    return (1.0 - nusq) * (top / bot);
}

/**
 * @brief Gridding correction for spheroidal function
 *
 * This is applied to the image after FFT to correct for the
 * gridding convolution.
 */
double gridding_correction(double x, int support) {
    // x is in range [0, 1] representing position in image
    // Correction is 1/FT(spheroidal)
    if (std::abs(x) < 1e-10) {
        return 1.0;
    }

    double arg = x * support;
    // Approximate correction
    double sinc = std::sin(M_PI * arg) / (M_PI * arg);
    double sph = spheroidal(x);

    if (std::abs(sph) < 1e-10) {
        return 1.0;
    }

    return 1.0 / sph;
}

//-----------------------------------------------------------------------------
// PSWF implementation
//-----------------------------------------------------------------------------

PSWF::PSWF(double m, double n)
    : m_(m), n_(n)
{
    compute_coefficients();
}

void PSWF::compute_coefficients() {
    // Polynomial coefficients for PSWF approximation
    // These are for m=6, n=0 which is standard for radio astronomy
    coeffs_ = {
        1.0,
        -0.4066667,
        0.03333333,
        -0.001269841,
        2.380952e-05,
        -2.380952e-07
    };
}

double PSWF::operator()(double x) const {
    if (std::abs(x) > 1.0) {
        return 0.0;
    }
    return spheroidal(x);
}

double PSWF::eval2d(double x, double y) const {
    return (*this)(x) * (*this)(y);
}

//-----------------------------------------------------------------------------
// PSWFConvolutionFunction implementation
//-----------------------------------------------------------------------------

PSWFConvolutionFunction::PSWFConvolutionFunction(
    int support,
    int oversampling,
    int n_mueller,
    int n_cubes)
    : support_(support)
    , oversampling_(oversampling)
    , n_mueller_(n_mueller)
    , n_cubes_(n_cubes)
{
    precompute();
}

void PSWFConvolutionFunction::precompute() {
    int full_support = 2 * support_ + 1;
    int oversampled_size = full_support * oversampling_;

    // Total size: x * y * mueller * cube
    size_t total = oversampled_size * oversampled_size * n_mueller_ * n_cubes_;
    values_.resize(total);

    // Precompute CF values
    for (int oy = 0; oy < oversampled_size; oy++) {
        for (int ox = 0; ox < oversampled_size; ox++) {
            // Convert to normalized coordinates [-1, 1]
            double x = (ox - support_ * oversampling_) / (double)(support_ * oversampling_);
            double y = (oy - support_ * oversampling_) / (double)(support_ * oversampling_);

            // Compute 2D PSWF value
            double val = pswf_.eval2d(x, y);

            // Store for all Mueller and cube indices (same value for basic gridding)
            for (int cube = 0; cube < n_cubes_; cube++) {
                for (int mueller = 0; mueller < n_mueller_; mueller++) {
                    size_t idx = index(ox, oy, mueller, cube);
                    values_[idx] = std::complex<hpg::cf_fp_t>(val, 0.0f);
                }
            }
        }
    }
}

size_t PSWFConvolutionFunction::index(int x, int y, int mueller, int cube) const {
    int full_support = 2 * support_ + 1;
    int oversampled_size = full_support * oversampling_;

    return ((cube * n_mueller_ + mueller) * oversampled_size + y) * oversampled_size + x;
}

unsigned PSWFConvolutionFunction::oversampling() const {
    return oversampling_;
}

hpg::coord_t PSWFConvolutionFunction::num_groups() const {
    return 1;  // Single CF group for basic gridding
}

std::array<hpg::coord_t, 4> PSWFConvolutionFunction::extents(hpg::coord_t grp) const {
    int full_support = 2 * support_ + 1;
    int oversampled_size = full_support * oversampling_;

    return {oversampled_size, oversampled_size, n_mueller_, n_cubes_};
}

PSWFConvolutionFunction::value_type PSWFConvolutionFunction::operator()(
    hpg::coord_t x,
    hpg::coord_t y,
    hpg::coord_t mueller,
    hpg::coord_t cube,
    hpg::coord_t grp) const
{
    int full_support = 2 * support_ + 1;
    int oversampled_size = full_support * oversampling_;

    if (x < 0 || x >= oversampled_size ||
        y < 0 || y >= oversampled_size ||
        mueller < 0 || mueller >= n_mueller_ ||
        cube < 0 || cube >= n_cubes_) {
        return value_type(0.0f, 0.0f);
    }

    return values_[index(x, y, mueller, cube)];
}

hpg::rval_t<std::string> PSWFConvolutionFunction::copy_to(
    hpg::Device device,
    hpg::Device host_device,
    hpg::coord_t grp,
    value_type* dst) const
{
    if (grp != 0) {
        return hpg::rval_t<std::string>("Invalid group index");
    }

    std::copy(values_.begin(), values_.end(), dst);
    return hpg::rval_t<std::string>();
}

std::array<hpg::coord_t, 2> PSWFConvolutionFunction::radii(hpg::coord_t grp) const {
    return {support_, support_};
}

//-----------------------------------------------------------------------------
// Helper functions
//-----------------------------------------------------------------------------

hpg::CFArrayShape create_cf_shape(
    int support,
    int oversampling,
    int n_mueller,
    int n_cubes)
{
    int full_support = 2 * support + 1;
    int oversampled_size = full_support * oversampling;

    // Create shape with 1 group
    std::vector<std::array<hpg::coord_t, 4>> group_extents;
    group_extents.push_back({oversampled_size, oversampled_size, n_mueller, n_cubes});

    return hpg::CFArrayShape(oversampling, std::move(group_extents));
}

} // namespace ugmrt
