/**
 * @file test_reader.c
 * @brief Test for slice-based reader_process
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "svio.h"
#include "gmrt_newcorr.h"
#include "svfits_reader.h"

static int vis_count = 0;
static int callback_errors = 0;

static int vis_callback(const CudaVisibility* vis, void* user_data) {
    vis_count++;

    // Print first few visibilities
    if (vis_count <= 10) {
        printf("  vis[%d]: ch=%d u=%.3f v=%.3f re=%.6f im=%.6f\n",
               vis_count, vis->channel, vis->u, vis->v, vis->re, vis->im);
    }

    // Print progress every 10000 visibilities
    if (vis_count % 10000 == 0) {
        printf("  ... processed %d visibilities\n", vis_count);
    }

    return 0;  // Continue processing
}

int main(int argc, char *argv[]) {
    char *param_file = "svfits_par.txt";
    char *antsamp_file = "antsamp.hdr";

    if (argc > 1) param_file = argv[1];
    if (argc > 2) antsamp_file = argv[2];

    printf("Test reader_process (slice-based iteration)\n");
    printf("  param_file: %s\n", param_file);
    printf("  antsamp_file: %s\n", antsamp_file);
    printf("\n");

    // Create reader config
    ReaderConfig config;
    memset(&config, 0, sizeof(config));
    strncpy(config.param_file, param_file, sizeof(config.param_file) - 1);
    strncpy(config.antsamp_file, antsamp_file, sizeof(config.antsamp_file) - 1);

    printf("Creating reader...\n");
    SvfitsReader reader = reader_create(&config);
    if (!reader) {
        fprintf(stderr, "Failed to create reader\n");
        return 1;
    }

    printf("Initializing reader...\n");
    if (reader_init(reader) != 0) {
        fprintf(stderr, "Failed to initialize reader\n");
        reader_free(reader);
        return 1;
    }

    // Get info
    BurstInfo burst;
    FreqInfo freq;
    reader_get_burst_info(reader, &burst);
    reader_get_freq_info(reader, &freq);

    printf("\nBurst info:\n");
    printf("  Name: %s\n", burst.name);
    printf("  MJD: %.6f\n", burst.mjd);
    printf("  DM: %.2f\n", burst.dm);
    printf("  Ref freq: %.2f MHz\n", burst.ref_freq_hz / 1e6);

    printf("\nFreq info:\n");
    printf("  Start: %.2f MHz\n", freq.freq_start_hz / 1e6);
    printf("  End: %.2f MHz\n", freq.freq_end_hz / 1e6);
    printf("  Channels: %d\n", freq.n_channels);
    printf("  Channel width: %.6f MHz\n", freq.channel_width_hz / 1e6);

    printf("\nProcessing visibilities...\n");
    fflush(stdout);

    size_t count = reader_process(reader, vis_callback, NULL);

    printf("\nProcessing complete:\n");
    printf("  Total visibilities: %zu\n", count);
    printf("  Callback count: %d\n", vis_count);
    printf("  Errors: %d\n", callback_errors);

    reader_free(reader);
    printf("\nTest completed successfully!\n");

    return 0;
}
