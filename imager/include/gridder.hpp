#pragma once
#include "config.hpp"
#include "reader.hpp"
#include "kernel.hpp"
#include <vector>
#include <complex>

struct Grid {
    int nx, ny;
    std::vector<std::complex<float>> data;
    std::vector<float> weights;
};

class Gridder {
public:
    Gridder(const Config& config, const WKernel& kernel);

    Grid createGrid();
    void gridVisibilities(Grid& grid, const std::vector<Visibility>& vis);
    Grid computePSF(const std::vector<Visibility>& vis);

    void normalize(Grid& grid);
    void fft(Grid& grid, bool forward = true);
    void fftshift(Grid& grid);

    std::vector<float> extractImage(const Grid& grid);

private:
    Config config_;
    WKernel kernel_;
    float scale_u_, scale_v_;

    void gridVisibility(Grid& grid, const Visibility& v);
    void computeScales();
};
