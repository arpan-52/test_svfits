/**
 * @file vis_to_hpg.cpp
 * @brief Visibility batch and conversion utilities
 */

#include "vis_to_hpg.hpp"
#include <cstring>

namespace ugmrt {

//-----------------------------------------------------------------------------
// VisibilityBatch implementation
//-----------------------------------------------------------------------------

VisibilityBatch::VisibilityBatch(size_t max_size)
    : max_size_(max_size)
{
    vis_data_.reserve(max_size);
}

bool VisibilityBatch::add(
    std::complex<float> vis_pol0,
    std::complex<float> vis_pol1,
    float weight0,
    float weight1,
    double u, double v, double w,
    double freq_hz,
    int channel,
    int cf_cube,
    int cf_group,
    float phase,
    float cf_phase_grad_x,
    float cf_phase_grad_y)
{
    hpg::VisData<2> vis;

    // Set visibility values
    vis.m_visibilities[0] = vis_pol0;
    vis.m_visibilities[1] = vis_pol1;

    // Set weights
    vis.m_weights[0] = weight0;
    vis.m_weights[1] = weight1;

    // Set UVW coordinates
    vis.m_uvw[0] = u;
    vis.m_uvw[1] = v;
    vis.m_uvw[2] = w;

    // Set frequency
    vis.m_frequency = freq_hz;

    // Set phase
    vis.m_phase = phase;

    // Set grid cube (frequency channel in output grid)
    vis.m_grid_cube = channel;

    // Set CF indices
    vis.m_cf_index[0] = cf_cube;
    vis.m_cf_index[1] = cf_group;

    // Set CF phase gradients
    vis.m_cf_phase_gradient[0] = cf_phase_grad_x;
    vis.m_cf_phase_gradient[1] = cf_phase_grad_y;

    vis_data_.push_back(vis);

    return vis_data_.size() >= max_size_;
}

std::vector<hpg::VisData<2>>& VisibilityBatch::get_vis_data() {
    return vis_data_;
}

void VisibilityBatch::clear() {
    vis_data_.clear();
}

bool VisibilityBatch::empty() const {
    return vis_data_.empty();
}

size_t VisibilityBatch::size() const {
    return vis_data_.size();
}

bool VisibilityBatch::full() const {
    return vis_data_.size() >= max_size_;
}

//-----------------------------------------------------------------------------
// MuellerMapping implementation
//-----------------------------------------------------------------------------

std::vector<std::array<hpg::coord_t, 2>> MuellerMapping::get_mueller_indexes(int n_pol) {
    std::vector<std::array<hpg::coord_t, 2>> result;

    if (n_pol == 1) {
        // Stokes I only
        result.push_back({0, 0});
    } else if (n_pol == 2) {
        // RR, LL or XX, YY
        result.push_back({0, 0});
        result.push_back({1, 1});
    } else if (n_pol == 4) {
        // Full polarization
        result.push_back({0, 0});
        result.push_back({0, 1});
        result.push_back({1, 0});
        result.push_back({1, 1});
    }

    return result;
}

std::vector<std::array<hpg::coord_t, 2>> MuellerMapping::get_conj_mueller_indexes(int n_pol) {
    std::vector<std::array<hpg::coord_t, 2>> result;

    if (n_pol == 1) {
        result.push_back({0, 0});
    } else if (n_pol == 2) {
        result.push_back({0, 0});
        result.push_back({1, 1});
    } else if (n_pol == 4) {
        // Conjugate swaps cross-terms
        result.push_back({0, 0});
        result.push_back({1, 0});
        result.push_back({0, 1});
        result.push_back({1, 1});
    }

    return result;
}

} // namespace ugmrt
