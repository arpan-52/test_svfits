#include "fits_output.hpp"
#include "utils.hpp"
#include <fitsio.h>
#include <iostream>
#include <cmath>

void writeFITS(const std::string& filename,
               const Grid& image,
               const Config& config) {
    std::cout << "=== Writing FITS ===" << std::endl;
    std::cout << "Output: " << filename << std::endl;

    fitsfile* fptr;
    int status = 0;

    // Delete existing file
    remove(filename.c_str());

    // Create FITS file
    fits_create_file(&fptr, filename.c_str(), &status);
    if (status) {
        fits_report_error(stderr, status);
        return;
    }

    // Create image
    long naxes[2] = {image.nx, image.ny};
    fits_create_img(fptr, FLOAT_IMG, 2, naxes, &status);

    // Extract real image
    std::vector<float> img_data(image.nx * image.ny);
    for (size_t i = 0; i < img_data.size(); i++) {
        img_data[i] = std::abs(image.data[i]);
    }

    // Write image data
    long fpixel[2] = {1, 1};
    fits_write_pix(fptr, TFLOAT, fpixel, img_data.size(),
                   img_data.data(), &status);

    // Write WCS headers
    double cell_deg = config.cell_size_arcsec / 3600.0;

    fits_write_key(fptr, TSTRING, "CTYPE1", (void*)"RA---SIN", NULL, &status);
    fits_write_key(fptr, TSTRING, "CTYPE2", (void*)"DEC--SIN", NULL, &status);

    double crpix1 = image.nx / 2.0 + 1;
    double crpix2 = image.ny / 2.0 + 1;
    fits_write_key(fptr, TDOUBLE, "CRPIX1", &crpix1, NULL, &status);
    fits_write_key(fptr, TDOUBLE, "CRPIX2", &crpix2, NULL, &status);

    double crval1 = config.ra_app * 180.0 / PI;
    double crval2 = config.dec_app * 180.0 / PI;
    fits_write_key(fptr, TDOUBLE, "CRVAL1", &crval1, NULL, &status);
    fits_write_key(fptr, TDOUBLE, "CRVAL2", &crval2, NULL, &status);

    double cdelt1 = -cell_deg;
    double cdelt2 = cell_deg;
    fits_write_key(fptr, TDOUBLE, "CDELT1", &cdelt1, NULL, &status);
    fits_write_key(fptr, TDOUBLE, "CDELT2", &cdelt2, NULL, &status);

    fits_write_key(fptr, TSTRING, "BUNIT", (void*)"JY/BEAM", NULL, &status);

    double freq_hz = config.ref_freq;
    fits_write_key(fptr, TDOUBLE, "FREQ", &freq_hz, "Reference frequency (Hz)", &status);

    fits_close_file(fptr, &status);

    if (status) {
        fits_report_error(stderr, status);
    } else {
        std::cout << "FITS file written successfully" << std::endl;
    }

    std::cout << "=== FITS Output Complete ===" << std::endl;
}
