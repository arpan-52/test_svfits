#include "gridder.hpp"
#include "utils.hpp"
#include <iostream>
#include <cmath>
#include <fftw3.h>
#include <omp.h>

Gridder::Gridder(const Config& config, const WKernel& kernel)
    : config_(config), kernel_(kernel)
{
    computeScales();
}

void Gridder::computeScales() {
    double cell_rad = config_.cell_size_arcsec * ARCSEC_TO_RAD;
    scale_u_ = config_.grid_nx * cell_rad;
    scale_v_ = config_.grid_ny * cell_rad;

    std::cout << "Grid scale: u=" << scale_u_ << ", v=" << scale_v_ << std::endl;
}

Grid Gridder::createGrid() {
    Grid grid;
    grid.nx = config_.grid_nx;
    grid.ny = config_.grid_ny;
    grid.data.resize(grid.nx * grid.ny, 0.0f);
    grid.weights.resize(grid.nx * grid.ny, 0.0f);
    return grid;
}

void Gridder::gridVisibilities(Grid& grid, const std::vector<Visibility>& vis) {
    std::cout << "=== Gridding Visibilities ===" << std::endl;
    std::cout << "Visibilities: " << vis.size() << std::endl;
    std::cout << "Grid size: " << grid.nx << " x " << grid.ny << std::endl;

    #pragma omp parallel for num_threads(config_.num_threads)
    for (size_t i = 0; i < vis.size(); i++) {
        if (vis[i].weight > 0) {
            gridVisibility(grid, vis[i]);
        }
    }

    std::cout << "=== Gridding Complete ===" << std::endl;
}

void Gridder::gridVisibility(Grid& grid, const Visibility& v) {
    float pos_u = scale_u_ * v.u + grid.nx / 2.0f;
    float pos_v = scale_v_ * v.v + grid.ny / 2.0f;

    int grid_u = static_cast<int>(std::round(pos_u));
    int grid_v = static_cast<int>(std::round(pos_v));

    int fine_u = static_cast<int>(std::round((grid_u - pos_u) * kernel_.oversampling));
    int fine_v = static_cast<int>(std::round((grid_v - pos_v) * kernel_.oversampling));

    int cf_minor_u = (fine_u >= 0) ? fine_u : kernel_.oversampling + fine_u;
    int cf_minor_v = (fine_v >= 0) ? fine_v : kernel_.oversampling + fine_v;

    float w_conj = (v.w > 0) ? -1.0f : 1.0f;

    std::complex<float> vis_val(v.re * v.weight, v.im * v.weight);

    for (int dv = -kernel_.support; dv <= kernel_.support; dv++) {
        int gv = grid_v + dv;
        if (gv < 0 || gv >= grid.ny) continue;

        int cf_y_maj = kernel_.support + dv;

        for (int du = -kernel_.support; du <= kernel_.support; du++) {
            int gu = grid_u + du;
            if (gu < 0 || gu >= grid.nx) continue;

            int cf_x_maj = kernel_.support + du;

            int cf_idx = cf_x_maj * kernel_.stride_x_major +
                        cf_y_maj * kernel_.stride_y_major +
                        v.w_plane * kernel_.stride_w_plane +
                        cf_minor_u * kernel_.stride_x_minor +
                        cf_minor_v * kernel_.stride_y_minor;

            std::complex<float> cf_val = kernel_.cf_values[cf_idx];
            cf_val.imag(cf_val.imag() * w_conj);

            std::complex<float> contrib = cf_val * vis_val;

            int grid_idx = gv * grid.nx + gu;

            #pragma omp atomic
            grid.data[grid_idx].real() += contrib.real();
            #pragma omp atomic
            grid.data[grid_idx].imag() += contrib.imag();
            #pragma omp atomic
            grid.weights[grid_idx] += std::abs(cf_val) * v.weight;
        }
    }
}

Grid Gridder::computePSF(const std::vector<Visibility>& vis) {
    std::cout << "=== Computing PSF ===" << std::endl;

    Grid psf = createGrid();

    std::vector<Visibility> unit_vis = vis;
    for (auto& v : unit_vis) {
        v.re = 1.0f;
        v.im = 0.0f;
    }

    gridVisibilities(psf, unit_vis);
    normalize(psf);
    fft(psf);
    fftshift(psf);

    std::cout << "=== PSF Complete ===" << std::endl;

    return psf;
}

void Gridder::normalize(Grid& grid) {
    std::cout << "Normalizing grid..." << std::endl;

    #pragma omp parallel for num_threads(config_.num_threads)
    for (int i = 0; i < grid.nx * grid.ny; i++) {
        if (grid.weights[i] > 0) {
            grid.data[i] /= grid.weights[i];
        }
    }
}

void Gridder::fft(Grid& grid, bool forward) {
    std::cout << "Computing FFT..." << std::endl;

    fftwf_complex* data = reinterpret_cast<fftwf_complex*>(grid.data.data());

    fftwf_plan plan = fftwf_plan_dft_2d(
        grid.ny, grid.nx,
        data, data,
        forward ? FFTW_FORWARD : FFTW_BACKWARD,
        FFTW_ESTIMATE
    );

    fftwf_execute(plan);
    fftwf_destroy_plan(plan);

    if (!forward) {
        float norm = 1.0f / (grid.nx * grid.ny);
        #pragma omp parallel for num_threads(config_.num_threads)
        for (size_t i = 0; i < grid.data.size(); i++) {
            grid.data[i] *= norm;
        }
    }
}

void Gridder::fftshift(Grid& grid) {
    int nx2 = grid.nx / 2;
    int ny2 = grid.ny / 2;

    for (int y = 0; y < ny2; y++) {
        for (int x = 0; x < nx2; x++) {
            std::swap(grid.data[y * grid.nx + x],
                     grid.data[(y + ny2) * grid.nx + (x + nx2)]);
            std::swap(grid.data[y * grid.nx + (x + nx2)],
                     grid.data[(y + ny2) * grid.nx + x]);
        }
    }
}

std::vector<float> Gridder::extractImage(const Grid& grid) {
    std::vector<float> image(grid.nx * grid.ny);

    #pragma omp parallel for num_threads(config_.num_threads)
    for (size_t i = 0; i < image.size(); i++) {
        image[i] = std::abs(grid.data[i]);
    }

    return image;
}
