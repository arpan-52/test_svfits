/**
 * @file svfits_reader.c
 * @brief svfits visibility reader implementation
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "svfits_reader.h"

// Include svfits headers
#include "svio.h"
#include "gmrt_newcorr.h"

//-----------------------------------------------------------------------------
// Reader implementation structure
//-----------------------------------------------------------------------------

struct SvfitsReaderImpl {
    ReaderConfig config;

    // svfits structures - use pointers like original svfits.c
    SvSelectionType user;
    // These are allocated separately with malloc (not embedded)

    // Derived info
    BurstInfo burst_info;
    FreqInfo freq_info;

    // Buffers
    char* raw_buffer;
    size_t raw_buffer_size;

    // Bandpass data
    float** bandpass;           // [baseline][channel]
    float** off_source_re;      // [baseline][channel]
    float** off_source_im;      // [baseline][channel]

    int initialized;
};

//-----------------------------------------------------------------------------
// Forward declarations for svfits functions
//-----------------------------------------------------------------------------

extern int init_user(SvSelectionType *user, char *fname, char *anthdr,
                     char *bhdrfile, char *bulletinA);
extern int read_slice(SvSelectionType *user, int idx, int slice, char *rbuf);
extern int get_file_order(SvSelectionType *user, int *order);
extern void init_mat(SvSelectionType *user, double tm);
extern float half_to_float(const unsigned short x);
extern double get_ha(SvSelectionType *user, double tm);
extern int robust_stats(int n, float *x, float *med, float *mad);

//-----------------------------------------------------------------------------
// Helper functions
//-----------------------------------------------------------------------------

static void compute_uvw(struct SvfitsReaderImpl* r, int bl, double mjd,
                        double* u, double* v, double* w) {
    int ant0 = r->user.vispar.visinfo[bl].ant0;
    int ant1 = r->user.vispar.visinfo[bl].ant1;

    double bx = r->user.corr->antenna[ant1].bx - r->user.corr->antenna[ant0].bx;
    double by = r->user.corr->antenna[ant1].by - r->user.corr->antenna[ant0].by;
    double bz = r->user.corr->antenna[ant1].bz - r->user.corr->antenna[ant0].bz;

    double ha = get_ha(&r->user, mjd);
    double dec = r->user.burst.dec_app;

    double sin_ha = sin(ha);
    double cos_ha = cos(ha);
    double sin_dec = sin(dec);
    double cos_dec = cos(dec);

    *u = sin_ha * bx + cos_ha * by;
    *v = -sin_dec * cos_ha * bx + sin_dec * sin_ha * by + cos_dec * bz;
    *w = cos_dec * cos_ha * bx - cos_dec * sin_ha * by + sin_dec * bz;

    // Convert to wavelengths
    double wavelength = 299792458.0 / r->freq_info.freq_start_hz;
    *u /= wavelength;
    *v /= wavelength;
    *w /= wavelength;
}

//-----------------------------------------------------------------------------
// Public API implementation
//-----------------------------------------------------------------------------

SvfitsReader reader_create(const ReaderConfig* config) {
    struct SvfitsReaderImpl* r = calloc(1, sizeof(struct SvfitsReaderImpl));
    if (!r) { fprintf(stderr, "Failed to alloc SvfitsReaderImpl\n"); return NULL; }
    r->config = *config;

    // Allocate svfits structures with malloc EXACTLY like original svfits.c does
    // See svfits.c main() lines 1128-1135
    r->user.hdr = (InitHdrType*)malloc(sizeof(InitHdrType));
    if (!r->user.hdr) { fprintf(stderr, "Failed to alloc hdr\n"); return NULL; }
    r->user.hdr->scans = 1;

    r->user.srec = (ScanRecType*)malloc(sizeof(ScanRecType));
    if (!r->user.srec) { fprintf(stderr, "Failed to alloc srec\n"); return NULL; }

    r->user.srec->scan = (ScanInfoType*)malloc(sizeof(ScanInfoType));
    if (!r->user.srec->scan) { fprintf(stderr, "Failed to alloc scan\n"); return NULL; }

    // Original only bzero's source, not whole scan
    bzero(&r->user.srec->scan->source, sizeof(SourceParType));

    r->user.corr = (CorrType*)malloc(sizeof(CorrType));
    if (!r->user.corr) { fprintf(stderr, "Failed to alloc corr\n"); return NULL; }

    return r;
}

int reader_init(SvfitsReader reader) {
    struct SvfitsReaderImpl* r = reader;

    char param_file[1024], antsamp[1024], bulletin[1024];
    strcpy(param_file, r->config.param_file);
    strcpy(antsamp, r->config.antsamp_file);

    char* bulletin_ptr = NULL;
    if (strlen(r->config.bulletin_a) > 0) {
        strcpy(bulletin, r->config.bulletin_a);
        bulletin_ptr = bulletin;
    }

    if (init_user(&r->user, param_file, antsamp, NULL, bulletin_ptr) != 0) {
        fprintf(stderr, "Failed to initialize svfits\n");
        return -1;
    }

    // Extract burst info
    strcpy(r->burst_info.name, r->user.burst.name);
    r->burst_info.mjd = r->user.burst.mjd;
    r->burst_info.dm = r->user.burst.DM;
    r->burst_info.intrinsic_width = r->user.burst.int_wd;
    r->burst_info.ref_freq_hz = r->user.burst.f;
    r->burst_info.ra_rad = r->user.burst.ra_app;
    r->burst_info.dec_rad = r->user.burst.dec_app;

    // Extract freq info
    r->freq_info.freq_start_hz = r->user.srec->scan->source.freq[0];
    r->freq_info.channel_width_hz = r->user.srec->scan->source.ch_width;
    r->freq_info.n_channels = r->user.corr->daspar.channels;
    r->freq_info.freq_end_hz = r->freq_info.freq_start_hz +
        r->freq_info.n_channels * r->freq_info.channel_width_hz;

    // Allocate buffers - one slice worth
    int das_baselines = r->user.corr->daspar.baselines;
    int channels = r->user.corr->daspar.channels;
    size_t recl = das_baselines * channels * sizeof(float);
    r->raw_buffer_size = MAX_REC_PER_SLICE * recl;
    r->raw_buffer = malloc(r->raw_buffer_size);
    if (!r->raw_buffer) {
        fprintf(stderr, "Failed to allocate raw buffer (%zu bytes)\n", r->raw_buffer_size);
        return -1;
    }

    // Allocate bandpass arrays for selected baselines
    int baselines = r->user.baselines;
    r->bandpass = malloc(baselines * sizeof(float*));
    r->off_source_re = malloc(baselines * sizeof(float*));
    r->off_source_im = malloc(baselines * sizeof(float*));

    for (int bl = 0; bl < baselines; bl++) {
        r->bandpass[bl] = calloc(channels, sizeof(float));
        r->off_source_re[bl] = calloc(channels, sizeof(float));
        r->off_source_im[bl] = calloc(channels, sizeof(float));

        for (int ch = 0; ch < channels; ch++) {
            r->bandpass[bl][ch] = 1.0f;
        }
    }

    r->initialized = 1;

    printf("Initialized reader:\n");
    printf("  Burst: %s at MJD %.6f\n", r->burst_info.name, r->burst_info.mjd);
    printf("  DM: %.2f\n", r->burst_info.dm);
    printf("  Freq: %.2f - %.2f MHz (%d channels)\n",
           r->freq_info.freq_start_hz / 1e6,
           r->freq_info.freq_end_hz / 1e6,
           r->freq_info.n_channels);
    printf("  Baselines: %d\n", r->user.baselines);

    return 0;
}

void reader_get_burst_info(SvfitsReader reader, BurstInfo* info) {
    *info = reader->burst_info;
}

void reader_get_freq_info(SvfitsReader reader, FreqInfo* info) {
    *info = reader->freq_info;
}

int reader_get_n_baselines(SvfitsReader reader) {
    return reader->user.baselines;
}

double reader_get_wavelength(SvfitsReader reader) {
    double center_freq = (reader->freq_info.freq_start_hz + reader->freq_info.freq_end_hz) / 2.0;
    return 299792458.0 / center_freq;
}

size_t reader_process(SvfitsReader reader, VisibilityCallback callback, void* user_data) {
    struct SvfitsReaderImpl* r = reader;
    if (!r->initialized) {
        fprintf(stderr, "Reader not initialized\n");
        return 0;
    }

    size_t count = 0;
    int channels = r->user.corr->daspar.channels;
    int baselines = r->user.baselines;
    int nfiles = r->user.recfile.nfiles;
    int rec_per_slice = r->user.recfile.rec_per_slice;
    size_t recl = r->user.corr->daspar.baselines * channels * sizeof(float);

    int file_order[MaxRecFiles];
    get_file_order(&r->user, file_order);

    printf("Processing: %d files, %d baselines, %d channels, rec_per_slice=%d\n",
           nfiles, baselines, channels, rec_per_slice);

    // Process each file in burst order
    for (int f = 0; f < nfiles; f++) {
        int file_idx = file_order[f];

        int start_rec = r->user.recfile.start_rec[file_idx];
        int n_rec = r->user.recfile.n_rec[file_idx];

        if (n_rec <= 0) continue;

        // Calculate which slices we need to read for this file
        int start_slice = start_rec / rec_per_slice;
        int end_rec = start_rec + n_rec - 1;
        int end_slice = end_rec / rec_per_slice;

        printf("  File %d: start_rec=%d n_rec=%d, slices %d-%d\n",
               file_idx, start_rec, n_rec, start_slice, end_slice);

        // Process each slice that contains burst data
        for (int slice = start_slice; slice <= end_slice; slice++) {
            // Read this slice
            if (read_slice(&r->user, file_idx, slice, r->raw_buffer) != 0) {
                fprintf(stderr, "Failed to read file %d slice %d\n", file_idx, slice);
                continue;
            }

            // Determine which records within this slice to process
            int slice_start_rec = slice * rec_per_slice;
            int slice_end_rec = slice_start_rec + rec_per_slice - 1;

            // Intersect with burst range [start_rec, start_rec + n_rec - 1]
            int proc_start = (start_rec > slice_start_rec) ? start_rec : slice_start_rec;
            int proc_end = (end_rec < slice_end_rec) ? end_rec : slice_end_rec;

            // Convert to slice-local indices (0 to rec_per_slice-1)
            int local_start = proc_start - slice_start_rec;
            int local_end = proc_end - slice_start_rec;

            // Process records within this slice
            for (int local_rec = local_start; local_rec <= local_end; local_rec++) {
                int global_rec = slice_start_rec + local_rec;
                char* rec_ptr = r->raw_buffer + local_rec * recl;

                // Time for this record
                double t_rec = r->user.recfile.t_start[file_idx] +
                    slice * r->user.recfile.slice_interval +
                    local_rec * r->user.recfile.t_slice / rec_per_slice;
                double mjd = r->user.recfile.mjd_ref + t_rec / 86400.0;

                init_mat(&r->user, mjd);

                // Get burst channel range for this record
                int start_ch = r->user.bpass.start_chan[local_rec];
                int end_ch = r->user.bpass.end_chan[local_rec];

                if (start_ch < 0) start_ch = 0;
                if (end_ch >= channels) end_ch = channels - 1;
                if (end_ch < start_ch) continue;  // No valid channels

                // Process each baseline
                for (int bl = 0; bl < baselines; bl++) {
                    unsigned short* vis_ptr = (unsigned short*)(rec_ptr +
                        r->user.vispar.visinfo[bl].off);

                    double u, v, w;
                    compute_uvw(r, bl, mjd, &u, &v, &w);

                    // Process channels with burst signal
                    for (int ch = start_ch; ch <= end_ch; ch++) {
                        float re = half_to_float(vis_ptr[ch * 2]);
                        float im = half_to_float(vis_ptr[ch * 2 + 1]);

                        if (!isfinite(re) || !isfinite(im)) continue;

                        // Flip imaginary if needed
                        if (r->user.vispar.visinfo[bl].flip) {
                            im = -im;
                        }

                        // Create visibility
                        CudaVisibility vis;
                        vis.re = re;
                        vis.im = im;
                        vis.weight = 1.0f;
                        vis.u = (float)u;
                        vis.v = (float)v;
                        vis.w = (float)w;
                        vis.channel = ch;

                        // Callback to gridder
                        if (callback(&vis, user_data) != 0) {
                            return count;
                        }
                        count++;
                    }
                }
            }
        }
    }

    return count;
}

void reader_free(SvfitsReader reader) {
    struct SvfitsReaderImpl* r = reader;
    if (!r) return;

    if (r->raw_buffer) free(r->raw_buffer);

    if (r->initialized) {
        int baselines = r->user.baselines;
        if (r->bandpass) {
            for (int bl = 0; bl < baselines; bl++) {
                free(r->bandpass[bl]);
            }
            free(r->bandpass);
        }
        if (r->off_source_re) {
            for (int bl = 0; bl < baselines; bl++) {
                free(r->off_source_re[bl]);
            }
            free(r->off_source_re);
        }
        if (r->off_source_im) {
            for (int bl = 0; bl < baselines; bl++) {
                free(r->off_source_im[bl]);
            }
            free(r->off_source_im);
        }
    }

    // Free svfits structures allocated with malloc
    if (r->user.srec) {
        if (r->user.srec->scan) free(r->user.srec->scan);
        free(r->user.srec);
    }
    if (r->user.corr) free(r->user.corr);
    if (r->user.hdr) free(r->user.hdr);

    free(r);
}
