#pragma once
#include <cstdint>
#include <vector>
#include <string>
#include <complex>

// Half-float conversion
float half_to_float(uint16_t h);
uint16_t float_to_half(float f);

// Robust statistics
int robust_stats(int n, const float* x, float& median, float& mad);

// String parsing
std::vector<std::string> split(const std::string& s, char delim);
void trim(std::string& s);
bool parseKeyValue(const std::string& line, std::string& key, std::string& value);

// Constants
constexpr double C_LIGHT = 299792458.0;  // m/s
constexpr double PI = 3.141592653589793;
constexpr double ARCSEC_TO_RAD = 4.84813681109536e-6;
