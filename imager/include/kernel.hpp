#pragma once
#include "config.hpp"
#include <vector>
#include <complex>

struct WKernel {
    int support;
    int oversampling;
    int n_w_planes;
    float max_w;

    std::vector<std::complex<float>> cf_values;

    int stride_x_major;
    int stride_y_major;
    int stride_mueller;
    int stride_w_plane;
    int stride_x_minor;
    int stride_y_minor;

    int cf_size;  // 2*support + 1
};

class KernelGenerator {
public:
    explicit KernelGenerator(const Config& config);

    WKernel generateWProjection(float max_w);

private:
    Config config_;

    void generateWScreen(int w_plane, float w_value,
                        std::vector<std::complex<float>>& screen);
    float pswf(float r, int support);
    void fft2d(std::vector<std::complex<float>>& data, int nx, int ny);
};
