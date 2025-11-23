/**
 * @file main.c
 * @brief uGMRT CUDA Imager - Command line interface
 *
 * Usage:
 *   ugmrt_cuda -u param_file -a antsamp_file -k cf_file [options]
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include <math.h>

#include "cuda_types.h"
#include "cuda_gridder.h"
#include "svfits_reader.h"

//-----------------------------------------------------------------------------
// Visibility batch accumulator
//-----------------------------------------------------------------------------

typedef struct {
    CudaVisibility* vis;
    int count;
    int capacity;

    // Gridding context
    UVGrid* grid;
    ConvolutionFunction* cf;
    float scale_u, scale_v;

    // Statistics
    size_t total;
    size_t flagged;
    size_t gridded;
} BatchContext;

static int batch_callback(const CudaVisibility* vis, void* user_data) {
    BatchContext* ctx = (BatchContext*)user_data;

    ctx->total++;

    if (vis->weight < 0) {
        ctx->flagged++;
        return 0;  // Continue
    }

    // Add to batch
    ctx->vis[ctx->count++] = *vis;

    // Grid when batch is full
    if (ctx->count >= ctx->capacity) {
        grid_visibilities(ctx->grid, ctx->vis, ctx->count,
                          ctx->cf, ctx->scale_u, ctx->scale_v);
        ctx->gridded += ctx->count;
        ctx->count = 0;
    }

    return 0;  // Continue
}

//-----------------------------------------------------------------------------
// FITS output declaration
//-----------------------------------------------------------------------------

extern int write_fits_image(
    const char* filename,
    const cuFloatComplex* image,
    int nx, int ny, int n_pol, int n_chan,
    double ra_deg, double dec_deg,
    double cell_size_deg,
    double ref_freq_hz,
    double freq_width_hz);

//-----------------------------------------------------------------------------
// Main
//-----------------------------------------------------------------------------

void print_usage(const char* prog) {
    printf("uGMRT CUDA Imager v1.0\n\n");
    printf("Usage: %s -u param_file -a antsamp_file -k cf_file [options]\n\n", prog);
    printf("Required:\n");
    printf("  -u FILE    Parameter file (svfits_par.txt format)\n");
    printf("  -a FILE    Antenna/sampler header file\n");
    printf("  -k FILE    Convolution kernel file\n\n");
    printf("Optional:\n");
    printf("  -o FILE    Output FITS image (default: burst_image.fits)\n");
    printf("  -n NX,NY   Grid size (default: 512,512)\n");
    printf("  -c CELL    Cell size in arcseconds (default: 1.0)\n");
    printf("  -b BATCH   Visibility batch size (default: 100000)\n");
    printf("  -T THRESH  Flagging threshold in MAD (default: 5.0)\n");
    printf("  --no-bandpass    Disable bandpass correction\n");
    printf("  --no-baseline    Disable baseline subtraction\n");
    printf("  --no-flag        Disable RFI flagging\n");
    printf("  -h, --help       Show this help\n");
}

int main(int argc, char* argv[]) {
    // Configuration
    char param_file[1024] = "";
    char antsamp_file[1024] = "";
    char cf_file[1024] = "";
    char output_file[1024] = "burst_image.fits";
    int nx = 512, ny = 512;
    double cell_size_asec = 1.0;
    int batch_size = 100000;
    float flag_threshold = 5.0f;
    int do_bandpass = 1, do_baseline = 1, do_flag = 1;

    // Parse arguments
    static struct option long_options[] = {
        {"no-bandpass", no_argument, NULL, 1},
        {"no-baseline", no_argument, NULL, 2},
        {"no-flag", no_argument, NULL, 3},
        {"help", no_argument, NULL, 'h'},
        {NULL, 0, NULL, 0}
    };

    int opt;
    while ((opt = getopt_long(argc, argv, "u:a:k:o:n:c:b:T:h", long_options, NULL)) != -1) {
        switch (opt) {
            case 'u': strncpy(param_file, optarg, sizeof(param_file)-1); break;
            case 'a': strncpy(antsamp_file, optarg, sizeof(antsamp_file)-1); break;
            case 'k': strncpy(cf_file, optarg, sizeof(cf_file)-1); break;
            case 'o': strncpy(output_file, optarg, sizeof(output_file)-1); break;
            case 'n':
                if (sscanf(optarg, "%d,%d", &nx, &ny) != 2) {
                    nx = ny = atoi(optarg);
                }
                break;
            case 'c': cell_size_asec = atof(optarg); break;
            case 'b': batch_size = atoi(optarg); break;
            case 'T': flag_threshold = atof(optarg); break;
            case 1: do_bandpass = 0; break;
            case 2: do_baseline = 0; break;
            case 3: do_flag = 0; break;
            case 'h': print_usage(argv[0]); return 0;
            default: print_usage(argv[0]); return 1;
        }
    }

    // Check required arguments
    if (strlen(param_file) == 0 || strlen(antsamp_file) == 0 || strlen(cf_file) == 0) {
        fprintf(stderr, "Error: -u, -a, and -k options are required\n\n");
        print_usage(argv[0]);
        return 1;
    }

    // Print configuration
    printf("========================================\n");
    printf("  uGMRT CUDA Imager v1.0\n");
    printf("========================================\n");
    printf("Configuration:\n");
    printf("  Parameter file: %s\n", param_file);
    printf("  Antenna file: %s\n", antsamp_file);
    printf("  CF file: %s\n", cf_file);
    printf("  Output: %s\n", output_file);
    printf("  Grid: %dx%d\n", nx, ny);
    printf("  Cell size: %.2f arcsec\n", cell_size_asec);
    printf("  Batch size: %d\n", batch_size);
    printf("  Bandpass: %s\n", do_bandpass ? "yes" : "no");
    printf("  Baseline: %s\n", do_baseline ? "yes" : "no");
    printf("  Flagging: %s (threshold: %.1f MAD)\n", do_flag ? "yes" : "no", flag_threshold);
    printf("========================================\n\n");

    clock_t start_time = clock();

    // Load convolution function
    printf("Loading convolution function...\n");
    ConvolutionFunction cf = cf_load(cf_file);

    // Create grid
    printf("Creating UV grid...\n");
    UVGrid grid = grid_create(nx, ny, 1, 1);

    // Initialize reader
    ReaderConfig reader_config = {0};
    strncpy(reader_config.param_file, param_file, sizeof(reader_config.param_file)-1);
    strncpy(reader_config.antsamp_file, antsamp_file, sizeof(reader_config.antsamp_file)-1);
    reader_config.do_bandpass = do_bandpass;
    reader_config.do_baseline = do_baseline;
    reader_config.do_flag = do_flag;
    reader_config.flag_threshold = flag_threshold;
    reader_config.num_threads = 4;

    SvfitsReader reader = reader_create(&reader_config);
    if (reader_init(reader) != 0) {
        fprintf(stderr, "Failed to initialize reader\n");
        return 1;
    }

    // Get frequency info for grid scaling
    FreqInfo freq_info;
    reader_get_freq_info(reader, &freq_info);

    BurstInfo burst_info;
    reader_get_burst_info(reader, &burst_info);

    // Compute grid scale
    // UV in wavelengths maps to grid: grid = u * scale + N/2
    // where scale = N * cell_rad (NOT inverted!)
    double cell_rad = cell_size_asec * M_PI / (180.0 * 3600.0);
    float scale_u = nx * cell_rad;
    float scale_v = ny * cell_rad;

    // Maximum UV extent (wavelengths) that fits on grid
    double max_uv = 1.0 / (2.0 * cell_rad);
    printf("Grid scale: u=%.6f, v=%.6f (max UV: %.1f wavelengths)\n",
           scale_u, scale_v, max_uv);

    // Process visibilities
    printf("\nProcessing visibilities...\n");

    BatchContext ctx = {0};
    ctx.vis = malloc(batch_size * sizeof(CudaVisibility));
    ctx.capacity = batch_size;
    ctx.grid = &grid;
    ctx.cf = &cf;
    ctx.scale_u = scale_u;
    ctx.scale_v = scale_v;

    clock_t grid_start = clock();
    reader_process(reader, batch_callback, &ctx);

    // Grid remaining visibilities
    if (ctx.count > 0) {
        grid_visibilities(&grid, ctx.vis, ctx.count, &cf, scale_u, scale_v);
        ctx.gridded += ctx.count;
    }
    clock_t grid_end = clock();

    printf("  Total: %zu visibilities\n", ctx.total);
    printf("  Flagged: %zu (%.1f%%)\n", ctx.flagged, 100.0 * ctx.flagged / ctx.total);
    printf("  Gridded: %zu\n", ctx.gridded);
    printf("  Gridding time: %.2f sec\n", (double)(grid_end - grid_start) / CLOCKS_PER_SEC);

    // Normalize
    printf("\nNormalizing...\n");
    grid_normalize(&grid);

    // Pre-FFT shift: move DC from center to corner (FFT convention)
    printf("Pre-FFT shift (DC from center to corner)...\n");
    grid_shift(&grid);

    // FFT - INVERSE to go from UV domain to image domain
    printf("Applying FFT (inverse, UV to image)...\n");
    clock_t fft_start = clock();
    grid_fft(&grid, 0);  // Inverse FFT: UV grid -> image
    clock_t fft_end = clock();
    printf("  FFT time: %.2f sec\n", (double)(fft_end - fft_start) / CLOCKS_PER_SEC);

    // Post-FFT shift: move DC from corner back to center
    printf("Post-FFT shift (DC from corner to center)...\n");
    grid_shift(&grid);

    // Gridding correction
    printf("Applying gridding correction...\n");
    grid_correct(&grid, &cf);

    // Copy to host
    grid_to_host(&grid);

    // Write FITS
    printf("\nWriting output...\n");
    double ra_deg = burst_info.ra_rad * 180.0 / M_PI;
    double dec_deg = burst_info.dec_rad * 180.0 / M_PI;
    double cell_deg = cell_size_asec / 3600.0;
    double center_freq = (freq_info.freq_start_hz + freq_info.freq_end_hz) / 2.0;

    write_fits_image(output_file, grid.h_grid, nx, ny, 1, 1,
                     ra_deg, dec_deg, cell_deg,
                     center_freq, freq_info.channel_width_hz);

    // Cleanup
    free(ctx.vis);
    reader_free(reader);
    grid_free(&grid);
    cf_free(&cf);

    clock_t end_time = clock();
    double total_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

    printf("\n========================================\n");
    printf("Total time: %.2f sec\n", total_time);
    printf("========================================\n");
    printf("\nDone!\n");

    return 0;
}
