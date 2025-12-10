#include "kernel.hpp"
#include "utils.hpp"
#include <iostream>
#include <cmath>
#include <fftw3.h>

KernelGenerator::KernelGenerator(const Config& config)
    : config_(config)
{
}

WKernel KernelGenerator::generateWProjection(float max_w) {
    std::cout << "=== Generating W-Projection Kernels ===" << std::endl;
    std::cout << "W-planes: " << config_.n_w_planes << std::endl;
    std::cout << "Support: " << config_.cf_support << " pixels" << std::endl;
    std::cout << "Oversampling: " << config_.cf_oversampling << "x" << std::endl;
    std::cout << "Max W: " << max_w << " wavelengths" << std::endl;

    WKernel kernel;
    kernel.support = config_.cf_support;
    kernel.oversampling = config_.cf_oversampling;
    kernel.n_w_planes = config_.n_w_planes;
    kernel.max_w = max_w;
    kernel.cf_size = 2 * kernel.support + 1;

    // Compute strides for 6D array
    kernel.stride_y_minor = 1;
    kernel.stride_x_minor = kernel.oversampling;
    kernel.stride_w_plane = kernel.stride_x_minor * kernel.oversampling;
    kernel.stride_mueller = kernel.stride_w_plane * config_.n_w_planes;
    kernel.stride_y_major = kernel.stride_mueller * 1;  // n_mueller = 1
    kernel.stride_x_major = kernel.stride_y_major * kernel.cf_size;

    size_t n_values = kernel.cf_size * kernel.cf_size * 1 * config_.n_w_planes *
                      kernel.oversampling * kernel.oversampling;

    kernel.cf_values.resize(n_values);

    std::cout << "Total CF samples: " << n_values << std::endl;

    // Generate each W-plane
    float w_min = -max_w;
    float w_max = max_w;

    for (int w_plane = 0; w_plane < config_.n_w_planes; w_plane++) {
        float w_value = (config_.n_w_planes == 1) ? 0.0f :
            w_min + (w_max - w_min) * w_plane / (config_.n_w_planes - 1);

        std::vector<std::complex<float>> screen;
        generateWScreen(w_plane, w_value, screen);

        // FFT to UV domain
        int fft_size = kernel.cf_size * kernel.oversampling;
        fft2d(screen, fft_size, fft_size);

        // Copy to kernel array
        for (int y_maj = 0; y_maj < kernel.cf_size; y_maj++) {
            for (int x_maj = 0; x_maj < kernel.cf_size; x_maj++) {
                for (int y_min = 0; y_min < kernel.oversampling; y_min++) {
                    for (int x_min = 0; x_min < kernel.oversampling; x_min++) {
                        int screen_idx = (y_maj * kernel.oversampling + y_min) * fft_size +
                                        (x_maj * kernel.oversampling + x_min);

                        int cf_idx = x_maj * kernel.stride_x_major +
                                    y_maj * kernel.stride_y_major +
                                    w_plane * kernel.stride_w_plane +
                                    x_min * kernel.stride_x_minor +
                                    y_min * kernel.stride_y_minor;

                        float norm = 1.0f / (fft_size * fft_size);
                        kernel.cf_values[cf_idx] = screen[screen_idx] * norm;
                    }
                }
            }
        }

        if ((w_plane + 1) % 10 == 0 || w_plane == config_.n_w_planes - 1) {
            std::cout << "Generated W-plane " << (w_plane + 1) << "/"
                      << config_.n_w_planes << " (w = " << w_value << ")" << std::endl;
        }
    }

    std::cout << "=== W-Projection Kernels Generated ===" << std::endl;

    return kernel;
}

void KernelGenerator::generateWScreen(int w_plane, float w_value,
                                       std::vector<std::complex<float>>& screen) {
    int fft_size = (2 * config_.cf_support + 1) * config_.cf_oversampling;
    screen.resize(fft_size * fft_size);

    double cell_rad = config_.cell_size_arcsec * ARCSEC_TO_RAD;

    for (int y = 0; y < fft_size; y++) {
        for (int x = 0; x < fft_size; x++) {
            float px = (x - fft_size / 2.0f) / config_.cf_oversampling;
            float py = (y - fft_size / 2.0f) / config_.cf_oversampling;

            double l = px * cell_rad;
            double m = py * cell_rad;

            double lm_sq = l*l + m*m;

            if (lm_sq < 1.0) {
                double n = std::sqrt(1.0 - lm_sq) - 1.0;
                double phase = 2.0 * PI * w_value * n;

                float r = std::sqrt(px*px + py*py);
                float taper = pswf(r, config_.cf_support);

                screen[y * fft_size + x] = std::complex<float>(
                    std::cos(phase) * taper,
                    std::sin(phase) * taper
                );
            } else {
                screen[y * fft_size + x] = 0.0f;
            }
        }
    }
}

float KernelGenerator::pswf(float r, int support) {
    float alpha = 1.0f;
    float x = r / support;

    if (x >= 1.0f) return 0.0f;

    return std::exp(-alpha * x * x);
}

void KernelGenerator::fft2d(std::vector<std::complex<float>>& data, int nx, int ny) {
    fftwf_complex* in = reinterpret_cast<fftwf_complex*>(data.data());
    fftwf_plan plan = fftwf_plan_dft_2d(ny, nx, in, in, FFTW_FORWARD, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
}
