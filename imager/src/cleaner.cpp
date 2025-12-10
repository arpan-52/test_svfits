#include "cleaner.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>

Cleaner::Cleaner(const Config& config)
    : config_(config)
{
}

std::vector<CleanComponent> Cleaner::clean(Grid& dirty, const Grid& psf) {
    std::cout << "=== CLEAN Deconvolution ===" << std::endl;
    std::cout << "Max iterations: " << config_.clean_niter << std::endl;
    std::cout << "Gain: " << config_.clean_gain << std::endl;
    std::cout << "Threshold: " << config_.clean_threshold << " Jy" << std::endl;

    std::vector<CleanComponent> components;
    Grid residual = dirty;

    for (int iter = 0; iter < config_.clean_niter; iter++) {
        int px, py;
        float peak_flux;
        findPeak(residual, px, py, peak_flux);

        if (std::fabs(peak_flux) < config_.clean_threshold) {
            std::cout << "CLEAN converged at iteration " << iter << std::endl;
            break;
        }

        float comp_flux = config_.clean_gain * peak_flux;
        components.push_back({px, py, comp_flux});

        subtractComponent(residual, psf, px, py, comp_flux);

        if (iter % 100 == 0) {
            std::cout << "Iteration " << iter << ": peak = " << peak_flux
                      << " Jy at (" << px << ", " << py << ")" << std::endl;
        }
    }

    dirty = residual;

    std::cout << "Total CLEAN components: " << components.size() << std::endl;
    std::cout << "=== CLEAN Complete ===" << std::endl;

    return components;
}

void Cleaner::findPeak(const Grid& image, int& x, int& y, float& val) {
    val = 0.0f;
    x = 0;
    y = 0;

    for (int iy = 0; iy < image.ny; iy++) {
        for (int ix = 0; ix < image.nx; ix++) {
            float abs_val = std::abs(image.data[iy * image.nx + ix]);
            if (abs_val > std::fabs(val)) {
                val = image.data[iy * image.nx + ix].real();
                x = ix;
                y = iy;
            }
        }
    }
}

void Cleaner::subtractComponent(Grid& image, const Grid& psf,
                                 int cx, int cy, float flux) {
    int psf_cx = psf.nx / 2;
    int psf_cy = psf.ny / 2;

    for (int y = 0; y < psf.ny; y++) {
        for (int x = 0; x < psf.nx; x++) {
            int img_x = cx + (x - psf_cx);
            int img_y = cy + (y - psf_cy);

            if (img_x >= 0 && img_x < image.nx &&
                img_y >= 0 && img_y < image.ny) {
                int psf_idx = y * psf.nx + x;
                int img_idx = img_y * image.nx + img_x;

                image.data[img_idx] -= flux * psf.data[psf_idx];
            }
        }
    }
}

Grid Cleaner::restore(const std::vector<CleanComponent>& components,
                      const Grid& residual) {
    std::cout << "Restoring CLEAN components..." << std::endl;

    Grid restored = residual;

    // Simple restoration: just add components back
    for (const auto& comp : components) {
        if (comp.x >= 0 && comp.x < restored.nx &&
            comp.y >= 0 && comp.y < restored.ny) {
            int idx = comp.y * restored.nx + comp.x;
            restored.data[idx] += comp.flux;
        }
    }

    return restored;
}

Grid Cleaner::computeCleanBeam(const Grid& psf) {
    // TODO: Fit Gaussian to PSF
    return psf;
}
