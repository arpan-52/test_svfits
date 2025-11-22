/**
 * @file cf_kernel.hpp
 * @brief Convolution function (gridding kernel) generation
 *
 * Implements Prolate Spheroidal Wave Function (PSWF) for optimal
 * anti-aliasing during visibility gridding.
 */

#ifndef CF_KERNEL_HPP
#define CF_KERNEL_HPP

#include <hpg/hpg.hpp>
#include <vector>
#include <complex>
#include <array>
#include <memory>

namespace ugmrt {

/**
 * @brief Prolate Spheroidal Wave Function evaluator
 *
 * The PSWF is optimal for gridding because it:
 * - Minimizes aliasing from sources outside the image
 * - Has compact support (finite extent)
 * - Is separable (2D = 1D x 1D)
 */
class PSWF {
public:
    /**
     * @brief Construct PSWF evaluator
     * @param m Shape parameter (typically 6)
     * @param n Order parameter (typically 0)
     */
    PSWF(double m = 6.0, double n = 0.0);

    /**
     * @brief Evaluate PSWF at given position
     * @param x Position in range [-1, 1]
     * @return PSWF value
     */
    double operator()(double x) const;

    /**
     * @brief Evaluate 2D separable PSWF
     * @param x X position in range [-1, 1]
     * @param y Y position in range [-1, 1]
     * @return PSWF value
     */
    double eval2d(double x, double y) const;

private:
    double m_, n_;
    std::vector<double> coeffs_;

    void compute_coefficients();
};

/**
 * @brief HPG-compatible Convolution Function Array
 *
 * Stores precomputed PSWF values in format required by HPG.
 */
class PSWFConvolutionFunction : public hpg::CFArray {
public:
    /**
     * @brief Construct convolution function
     * @param support Half-width in pixels (full width = 2*support + 1)
     * @param oversampling Subpixel oversampling factor
     * @param n_mueller Number of Mueller matrix elements
     * @param n_cubes Number of frequency cubes
     */
    PSWFConvolutionFunction(
        int support,
        int oversampling,
        int n_mueller,
        int n_cubes
    );

    /**
     * @brief Get oversampling factor
     */
    unsigned oversampling() const override;

    /**
     * @brief Get number of CF groups
     */
    hpg::coord_t num_groups() const override;

    /**
     * @brief Get extents for a group
     */
    std::array<hpg::coord_t, 4> extents(hpg::coord_t grp) const override;

    /**
     * @brief Get CF value at given indices
     */
    value_type operator()(
        hpg::coord_t x,
        hpg::coord_t y,
        hpg::coord_t mueller,
        hpg::coord_t cube,
        hpg::coord_t grp
    ) const override;

    /**
     * @brief Copy CF data to device buffer
     */
    hpg::rval_t<std::string> copy_to(
        hpg::Device device,
        hpg::Device host_device,
        hpg::coord_t grp,
        value_type* dst
    ) const override;

    /**
     * @brief Get CF radii (half-widths)
     */
    std::array<hpg::coord_t, 2> radii(hpg::coord_t grp) const override;

private:
    int support_;
    int oversampling_;
    int n_mueller_;
    int n_cubes_;
    PSWF pswf_;

    // Precomputed values: [x][y][mueller][cube]
    std::vector<std::complex<hpg::cf_fp_t>> values_;

    void precompute();
    size_t index(int x, int y, int mueller, int cube) const;
};

/**
 * @brief Create CFArray shape descriptor
 */
hpg::CFArrayShape create_cf_shape(
    int support,
    int oversampling,
    int n_mueller,
    int n_cubes
);

/**
 * @brief Simple spheroidal function approximation
 * Uses polynomial approximation from AIPS/CASA
 */
double spheroidal(double nu);

/**
 * @brief Gridding correction function
 * Returns the correction factor to apply to the image after FFT
 */
double gridding_correction(double x, int support);

} // namespace ugmrt

#endif // CF_KERNEL_HPP
