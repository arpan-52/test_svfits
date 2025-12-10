#pragma once
#include "config.hpp"
#include "gridder.hpp"
#include <string>

void writeFITS(const std::string& filename,
               const Grid& image,
               const Config& config);
