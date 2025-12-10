#include "utils.hpp"
#include <algorithm>
#include <sstream>
#include <cmath>

// Half-float conversion (IEEE 754-2008 binary16)
float half_to_float(uint16_t h) {
    uint32_t sign = (h & 0x8000) << 16;
    uint32_t exp = (h & 0x7C00) >> 10;
    uint32_t mant = (h & 0x03FF);

    if (exp == 0) {
        if (mant == 0) {
            // Zero
            uint32_t f = sign;
            return *reinterpret_cast<float*>(&f);
        }
        // Denormal
        exp = 1;
        while (!(mant & 0x0400)) {
            mant <<= 1;
            exp--;
        }
        mant &= 0x03FF;
    } else if (exp == 31) {
        // Inf or NaN
        uint32_t f = sign | 0x7F800000 | (mant << 13);
        return *reinterpret_cast<float*>(&f);
    }

    uint32_t f = sign | ((exp + 112) << 23) | (mant << 13);
    return *reinterpret_cast<float*>(&f);
}

uint16_t float_to_half(float f) {
    uint32_t bits = *reinterpret_cast<uint32_t*>(&f);
    uint16_t sign = (bits >> 16) & 0x8000;
    int32_t exp = ((bits >> 23) & 0xFF) - 112;
    uint32_t mant = bits & 0x007FFFFF;

    if (exp <= 0) {
        if (exp < -10) return sign;  // Too small
        mant |= 0x00800000;
        uint16_t t = (mant >> (14 - exp));
        return sign | t;
    }
    if (exp > 30) {
        return sign | 0x7C00;  // Overflow to infinity
    }
    return sign | (exp << 10) | (mant >> 13);
}

// Robust statistics (from stats.c)
extern "C" int robust_stats(int n, float* x, float* median, float* mad);

// String utilities
std::vector<std::string> split(const std::string& s, char delim) {
    std::vector<std::string> result;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        trim(item);
        if (!item.empty()) {
            result.push_back(item);
        }
    }
    return result;
}

void trim(std::string& s) {
    // Left trim
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
        return !std::isspace(ch);
    }));
    // Right trim
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

bool parseKeyValue(const std::string& line, std::string& key, std::string& value) {
    size_t pos = line.find_first_of(" \t");
    if (pos == std::string::npos) return false;

    key = line.substr(0, pos);
    trim(key);

    value = line.substr(pos + 1);
    // Remove comment
    size_t comment_pos = value.find('!');
    if (comment_pos != std::string::npos) {
        value = value.substr(0, comment_pos);
    }
    trim(value);

    return !key.empty() && !value.empty();
}
