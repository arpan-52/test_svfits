/**
 * @file fits_output.c
 * @brief FITS image output
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fitsio.h>

#include "cuda_types.h"

/**
 * @brief Write image to FITS file
 */
int write_fits_image(
    const char* filename,
    const cuFloatComplex* image,
    int nx, int ny, int n_pol, int n_chan,
    double ra_deg, double dec_deg,
    double cell_size_deg,
    double ref_freq_hz,
    double freq_width_hz)
{
    fitsfile* fptr = NULL;
    int status = 0;

    // Remove existing file
    char full_path[1024];
    snprintf(full_path, sizeof(full_path), "!%s", filename);

    // Create FITS file
    if (fits_create_file(&fptr, full_path, &status)) {
        fprintf(stderr, "Error creating FITS file: %s\n", filename);
        return -1;
    }

    // Image dimensions
    long naxes[4] = {nx, ny, n_pol, n_chan};
    int naxis = 4;

    // Create image HDU
    if (fits_create_img(fptr, FLOAT_IMG, naxis, naxes, &status)) {
        fprintf(stderr, "Error creating image HDU\n");
        fits_close_file(fptr, &status);
        return -1;
    }

    // WCS headers
    double crpix1 = nx / 2.0 + 1;
    double crpix2 = ny / 2.0 + 1;
    double cdelt1 = -cell_size_deg;
    double cdelt2 = cell_size_deg;

    fits_write_key(fptr, TSTRING, "OBJECT", "uGMRT_BURST", "Target name", &status);
    fits_write_key(fptr, TSTRING, "TELESCOP", "uGMRT", "Telescope", &status);
    fits_write_key(fptr, TSTRING, "INSTRUME", "SPOTLIGHT", "Instrument", &status);
    fits_write_key(fptr, TSTRING, "ORIGIN", "ugmrt_cuda", "Software", &status);

    // Axis 1: RA
    fits_write_key(fptr, TSTRING, "CTYPE1", "RA---SIN", "Coordinate type", &status);
    fits_write_key(fptr, TDOUBLE, "CRPIX1", &crpix1, "Reference pixel", &status);
    fits_write_key(fptr, TDOUBLE, "CRVAL1", &ra_deg, "Reference value (deg)", &status);
    fits_write_key(fptr, TDOUBLE, "CDELT1", &cdelt1, "Pixel size (deg)", &status);
    fits_write_key(fptr, TSTRING, "CUNIT1", "deg", "Axis unit", &status);

    // Axis 2: Dec
    fits_write_key(fptr, TSTRING, "CTYPE2", "DEC--SIN", "Coordinate type", &status);
    fits_write_key(fptr, TDOUBLE, "CRPIX2", &crpix2, "Reference pixel", &status);
    fits_write_key(fptr, TDOUBLE, "CRVAL2", &dec_deg, "Reference value (deg)", &status);
    fits_write_key(fptr, TDOUBLE, "CDELT2", &cdelt2, "Pixel size (deg)", &status);
    fits_write_key(fptr, TSTRING, "CUNIT2", "deg", "Axis unit", &status);

    // Axis 3: Stokes
    double crpix3 = 1.0, crval3 = 1.0, cdelt3 = 1.0;
    fits_write_key(fptr, TSTRING, "CTYPE3", "STOKES", "Coordinate type", &status);
    fits_write_key(fptr, TDOUBLE, "CRPIX3", &crpix3, "Reference pixel", &status);
    fits_write_key(fptr, TDOUBLE, "CRVAL3", &crval3, "Reference value", &status);
    fits_write_key(fptr, TDOUBLE, "CDELT3", &cdelt3, "Increment", &status);

    // Axis 4: Frequency
    double crpix4 = 1.0;
    fits_write_key(fptr, TSTRING, "CTYPE4", "FREQ", "Coordinate type", &status);
    fits_write_key(fptr, TDOUBLE, "CRPIX4", &crpix4, "Reference pixel", &status);
    fits_write_key(fptr, TDOUBLE, "CRVAL4", &ref_freq_hz, "Reference value (Hz)", &status);
    fits_write_key(fptr, TDOUBLE, "CDELT4", &freq_width_hz, "Channel width (Hz)", &status);
    fits_write_key(fptr, TSTRING, "CUNIT4", "Hz", "Axis unit", &status);

    fits_write_key(fptr, TSTRING, "BUNIT", "JY/BEAM", "Brightness unit", &status);
    double equinox = 2000.0;
    fits_write_key(fptr, TDOUBLE, "EQUINOX", &equinox, "Equinox", &status);

    // Convert complex to real (take real part)
    long total = nx * ny * n_pol * n_chan;
    float* real_data = malloc(total * sizeof(float));
    for (long i = 0; i < total; i++) {
        real_data[i] = image[i].x;
    }

    // Write data
    long fpixel[4] = {1, 1, 1, 1};
    if (fits_write_pix(fptr, TFLOAT, fpixel, total, real_data, &status)) {
        fprintf(stderr, "Error writing image data\n");
        free(real_data);
        fits_close_file(fptr, &status);
        return -1;
    }

    free(real_data);
    fits_close_file(fptr, &status);

    if (status) {
        char errmsg[80];
        fits_get_errstatus(status, errmsg);
        fprintf(stderr, "FITS error: %s\n", errmsg);
        return -1;
    }

    printf("Wrote FITS image: %s (%dx%d)\n", filename, nx, ny);
    return 0;
}
