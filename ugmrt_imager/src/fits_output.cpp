/**
 * @file fits_output.cpp
 * @brief FITS image output utilities
 */

#include "ugmrt_imager.hpp"
#include <fitsio.h>
#include <cstring>
#include <cmath>
#include <stdexcept>
#include <iostream>

namespace ugmrt {

/**
 * @brief Write image to FITS file
 *
 * Creates a standard FITS image file with proper WCS headers
 * for radio astronomy imaging.
 */
bool write_fits_image(
    const std::string& filename,
    const std::complex<double>* image_data,
    int nx, int ny, int n_pol, int n_chan,
    double ra_deg, double dec_deg,
    double cell_size_deg,
    double ref_freq_hz,
    double freq_width_hz)
{
    fitsfile* fptr = nullptr;
    int status = 0;

    // Remove existing file
    std::string full_path = "!" + filename;  // ! = overwrite

    // Create FITS file
    if (fits_create_file(&fptr, full_path.c_str(), &status)) {
        std::cerr << "Error creating FITS file: " << filename << std::endl;
        return false;
    }

    // Image dimensions: [nx, ny, pol, freq]
    long naxes[4] = {nx, ny, n_pol, n_chan};
    int naxis = 4;

    // Create image HDU (use FLOAT for output)
    if (fits_create_img(fptr, FLOAT_IMG, naxis, naxes, &status)) {
        std::cerr << "Error creating image HDU" << std::endl;
        fits_close_file(fptr, &status);
        return false;
    }

    // Write WCS headers
    double crpix1 = nx / 2.0 + 1;
    double crpix2 = ny / 2.0 + 1;
    double cdelt1 = -cell_size_deg;  // RA increases left
    double cdelt2 = cell_size_deg;   // Dec increases up

    // FITS keywords
    fits_write_key(fptr, TSTRING, "OBJECT", (void*)"uGMRT_BURST", "Target name", &status);
    fits_write_key(fptr, TSTRING, "TELESCOP", (void*)"uGMRT", "Telescope name", &status);
    fits_write_key(fptr, TSTRING, "INSTRUME", (void*)"SPOTLIGHT", "Instrument", &status);
    fits_write_key(fptr, TSTRING, "ORIGIN", (void*)"ugmrt_imager", "Software", &status);

    // Axis 1: RA
    fits_write_key(fptr, TSTRING, "CTYPE1", (void*)"RA---SIN", "Coordinate type", &status);
    fits_write_key(fptr, TDOUBLE, "CRPIX1", &crpix1, "Reference pixel", &status);
    fits_write_key(fptr, TDOUBLE, "CRVAL1", &ra_deg, "Reference value (deg)", &status);
    fits_write_key(fptr, TDOUBLE, "CDELT1", &cdelt1, "Pixel size (deg)", &status);
    fits_write_key(fptr, TSTRING, "CUNIT1", (void*)"deg", "Axis unit", &status);

    // Axis 2: Dec
    fits_write_key(fptr, TSTRING, "CTYPE2", (void*)"DEC--SIN", "Coordinate type", &status);
    fits_write_key(fptr, TDOUBLE, "CRPIX2", &crpix2, "Reference pixel", &status);
    fits_write_key(fptr, TDOUBLE, "CRVAL2", &dec_deg, "Reference value (deg)", &status);
    fits_write_key(fptr, TDOUBLE, "CDELT2", &cdelt2, "Pixel size (deg)", &status);
    fits_write_key(fptr, TSTRING, "CUNIT2", (void*)"deg", "Axis unit", &status);

    // Axis 3: Stokes
    double crpix3 = 1.0;
    double crval3 = 1.0;  // Stokes I
    double cdelt3 = 1.0;
    fits_write_key(fptr, TSTRING, "CTYPE3", (void*)"STOKES", "Coordinate type", &status);
    fits_write_key(fptr, TDOUBLE, "CRPIX3", &crpix3, "Reference pixel", &status);
    fits_write_key(fptr, TDOUBLE, "CRVAL3", &crval3, "Reference value", &status);
    fits_write_key(fptr, TDOUBLE, "CDELT3", &cdelt3, "Increment", &status);

    // Axis 4: Frequency
    double crpix4 = 1.0;
    fits_write_key(fptr, TSTRING, "CTYPE4", (void*)"FREQ", "Coordinate type", &status);
    fits_write_key(fptr, TDOUBLE, "CRPIX4", &crpix4, "Reference pixel", &status);
    fits_write_key(fptr, TDOUBLE, "CRVAL4", &ref_freq_hz, "Reference value (Hz)", &status);
    fits_write_key(fptr, TDOUBLE, "CDELT4", &freq_width_hz, "Channel width (Hz)", &status);
    fits_write_key(fptr, TSTRING, "CUNIT4", (void*)"Hz", "Axis unit", &status);

    // Additional keywords
    fits_write_key(fptr, TSTRING, "BUNIT", (void*)"JY/BEAM", "Brightness unit", &status);
    double equinox = 2000.0;
    fits_write_key(fptr, TDOUBLE, "EQUINOX", &equinox, "Equinox of coordinates", &status);

    // Convert complex data to real (take real part = Stokes I for unpolarized)
    std::vector<float> real_data(nx * ny * n_pol * n_chan);
    size_t total = nx * ny * n_pol * n_chan;
    for (size_t i = 0; i < total; i++) {
        real_data[i] = static_cast<float>(image_data[i].real());
    }

    // Write image data
    long fpixel[4] = {1, 1, 1, 1};
    long nelements = nx * ny * n_pol * n_chan;
    if (fits_write_pix(fptr, TFLOAT, fpixel, nelements, real_data.data(), &status)) {
        std::cerr << "Error writing image data" << std::endl;
        fits_close_file(fptr, &status);
        return false;
    }

    // Close file
    fits_close_file(fptr, &status);

    if (status) {
        char errmsg[80];
        fits_get_errstatus(status, errmsg);
        std::cerr << "FITS error: " << errmsg << std::endl;
        return false;
    }

    std::cout << "Wrote FITS image: " << filename << " (" << nx << "x" << ny << ")" << std::endl;
    return true;
}

/**
 * @brief Apply gridding correction to image
 *
 * Corrects for the convolution with the gridding kernel by
 * dividing by the FT of the kernel.
 */
void apply_gridding_correction(
    std::complex<double>* image,
    int nx, int ny,
    int support)
{
    for (int iy = 0; iy < ny; iy++) {
        // Normalized y coordinate from center
        double y = (iy - ny / 2.0) / (ny / 2.0);

        for (int ix = 0; ix < nx; ix++) {
            // Normalized x coordinate from center
            double x = (ix - nx / 2.0) / (nx / 2.0);

            // Get correction factor
            double corr_x = gridding_correction(x, support);
            double corr_y = gridding_correction(y, support);
            double correction = corr_x * corr_y;

            // Apply correction
            size_t idx = iy * nx + ix;
            image[idx] *= correction;
        }
    }
}

} // namespace ugmrt
