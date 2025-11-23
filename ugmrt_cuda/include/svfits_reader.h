/**
 * @file svfits_reader.h
 * @brief C interface to svfits visibility reading
 */

#ifndef SVFITS_READER_H
#define SVFITS_READER_H

#include "cuda_types.h"

#ifdef __cplusplus
extern "C" {
#endif

//-----------------------------------------------------------------------------
// Reader configuration
//-----------------------------------------------------------------------------

typedef struct {
    char param_file[1024];      // svfits_par.txt path
    char antsamp_file[1024];    // antsamp.hdr path
    char bulletin_a[1024];      // Optional Bulletin A file

    int do_bandpass;            // Apply bandpass correction
    int do_baseline;            // Subtract baseline
    int do_flag;                // Apply RFI flagging
    float flag_threshold;       // Flagging threshold (MAD)
    int num_threads;            // OpenMP threads
} ReaderConfig;

//-----------------------------------------------------------------------------
// Burst and frequency info (output from reader)
//-----------------------------------------------------------------------------

typedef struct {
    char name[256];
    double mjd;
    double dm;
    double intrinsic_width;
    double ref_freq_hz;
    double ra_rad;              // RA in radians
    double dec_rad;             // Dec in radians
} BurstInfo;

typedef struct {
    double freq_start_hz;
    double freq_end_hz;
    int n_channels;
    double channel_width_hz;
} FreqInfo;

//-----------------------------------------------------------------------------
// Reader handle
//-----------------------------------------------------------------------------

typedef struct SvfitsReaderImpl* SvfitsReader;

/**
 * @brief Create reader
 */
SvfitsReader reader_create(const ReaderConfig* config);

/**
 * @brief Initialize reader (parse config, allocate buffers)
 * @return 0 on success, -1 on error
 */
int reader_init(SvfitsReader reader);

/**
 * @brief Get burst info
 */
void reader_get_burst_info(SvfitsReader reader, BurstInfo* info);

/**
 * @brief Get frequency info
 */
void reader_get_freq_info(SvfitsReader reader, FreqInfo* info);

/**
 * @brief Get number of baselines
 */
int reader_get_n_baselines(SvfitsReader reader);

/**
 * @brief Get reference wavelength (meters)
 */
double reader_get_wavelength(SvfitsReader reader);

/**
 * @brief Callback for processed visibilities
 *
 * @param vis Visibility data
 * @param user_data User-provided context
 * @return 0 to continue, non-zero to stop
 */
typedef int (*VisibilityCallback)(const CudaVisibility* vis, void* user_data);

/**
 * @brief Process all visibilities
 *
 * Reads raw data, applies svfits processing (bandpass, baseline, flagging),
 * and calls callback for each visibility.
 *
 * @param reader Reader handle
 * @param callback Function to call for each visibility
 * @param user_data Context passed to callback
 * @return Number of visibilities processed
 */
size_t reader_process(SvfitsReader reader, VisibilityCallback callback, void* user_data);

/**
 * @brief Free reader
 */
void reader_free(SvfitsReader reader);

//-----------------------------------------------------------------------------
// Batch/Parallel Processing API
//-----------------------------------------------------------------------------

/**
 * @brief Visibility buffer for batch processing
 */
typedef struct {
    CudaVisibility* vis;    // Visibility array
    size_t count;           // Number of visibilities
    size_t capacity;        // Allocated capacity
} VisibilityBuffer;

/**
 * @brief Create visibility buffer
 */
VisibilityBuffer* visbuf_create(size_t capacity);

/**
 * @brief Free visibility buffer
 */
void visbuf_free(VisibilityBuffer* buf);

/**
 * @brief Reset buffer (keep allocation)
 */
void visbuf_reset(VisibilityBuffer* buf);

/**
 * @brief Extract all visibilities to buffer (for batch GPU processing)
 *
 * This is more efficient than callback-based processing when you want to
 * accumulate all visibilities before sending to GPU.
 *
 * @param reader Reader handle
 * @param buf Output buffer (will be resized if needed)
 * @return Number of visibilities extracted, or -1 on error
 */
ssize_t reader_extract_all(SvfitsReader reader, VisibilityBuffer* buf);

/**
 * @brief Get estimated visibility count (for pre-allocation)
 */
size_t reader_estimate_vis_count(SvfitsReader reader);

#ifdef __cplusplus
}
#endif

#endif // SVFITS_READER_H
