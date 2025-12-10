#pragma once
#include "config.hpp"
#include "gridder.hpp"
#include <vector>

struct CleanComponent {
    int x, y;
    float flux;
};

class Cleaner {
public:
    explicit Cleaner(const Config& config);

    std::vector<CleanComponent> clean(Grid& dirty, const Grid& psf);
    Grid restore(const std::vector<CleanComponent>& components,
                 const Grid& residual);

private:
    Config config_;

    void findPeak(const Grid& image, int& x, int& y, float& val);
    void subtractComponent(Grid& image, const Grid& psf,
                          int x, int y, float flux);
    Grid computeCleanBeam(const Grid& psf);
};
