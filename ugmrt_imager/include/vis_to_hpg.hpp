/**
 * @file vis_to_hpg.hpp
 * @brief Convert svfits visibilities to HPG format
 */

#ifndef VIS_TO_HPG_HPP
#define VIS_TO_HPG_HPP

#include <hpg/hpg.hpp>
#include <vector>
#include <complex>

namespace ugmrt {

/**
 * @brief Visibility batch for efficient GPU transfer
 */
class VisibilityBatch {
public:
    explicit VisibilityBatch(size_t max_size = 10000);

    /**
     * @brief Add a visibility to the batch
     * @return true if batch is full after adding
     */
    bool add(
        std::complex<float> vis_pol0,
        std::complex<float> vis_pol1,
        float weight0,
        float weight1,
        double u, double v, double w,
        double freq_hz,
        int channel,
        int cf_cube,
        int cf_group,
        float phase = 0.0f,
        float cf_phase_grad_x = 0.0f,
        float cf_phase_grad_y = 0.0f
    );

    /**
     * @brief Get HPG-compatible visibility data
     */
    std::vector<hpg::VisData<2>>& get_vis_data();

    /**
     * @brief Clear the batch
     */
    void clear();

    /**
     * @brief Check if batch is empty
     */
    bool empty() const;

    /**
     * @brief Get current size
     */
    size_t size() const;

    /**
     * @brief Check if batch is full
     */
    bool full() const;

private:
    size_t max_size_;
    std::vector<hpg::VisData<2>> vis_data_;
};

/**
 * @brief Convert svfits Complex to std::complex
 */
inline std::complex<float> to_complex(float re, float im) {
    return std::complex<float>(re, im);
}

/**
 * @brief Encode baseline ID for FITS convention
 * baseline = 256 * ant1 + ant2 + 257 (for ant1 < ant2)
 */
inline int encode_baseline(int ant0, int ant1) {
    if (ant0 < ant1) {
        return 256 * (ant0 + 1) + (ant1 + 1);
    } else {
        return 256 * (ant1 + 1) + (ant0 + 1);
    }
}

/**
 * @brief Mueller index mapping for polarization products
 */
struct MuellerMapping {
    // For 2-pol (RR, LL) or (XX, YY)
    static constexpr int n_pol_2_indexes[2][2] = {
        {0, 0},  // pol 0 -> mueller row 0
        {1, 1}   // pol 1 -> mueller row 1
    };

    // For 4-pol (RR, RL, LR, LL) or (XX, XY, YX, YY)
    static constexpr int n_pol_4_indexes[4][2] = {
        {0, 0}, {0, 1}, {1, 0}, {1, 1}
    };

    /**
     * @brief Get Mueller indexes for given polarization count
     */
    static std::vector<std::array<hpg::coord_t, 2>> get_mueller_indexes(int n_pol);

    /**
     * @brief Get conjugate Mueller indexes
     */
    static std::vector<std::array<hpg::coord_t, 2>> get_conj_mueller_indexes(int n_pol);
};

} // namespace ugmrt

#endif // VIS_TO_HPG_HPP
