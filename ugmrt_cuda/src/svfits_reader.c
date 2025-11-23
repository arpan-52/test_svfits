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

static void compute_bandpass(struct SvfitsReaderImpl* r, int file_idx, char* raw_buf) {
    int channels = r->user.corr->daspar.channels;
    int baselines = r->user.baselines;
    size_t recl = r->user.corr->daspar.baselines * channels * sizeof(float);

    short* start_chan = r->user.bpass.start_chan;
    short* end_chan = r->user.bpass.end_chan;

    for (int bl = 0; bl < baselines; bl++) {
        float* sum = calloc(channels, sizeof(float));
        int* count = calloc(channels, sizeof(int));

        for (int rec = 0; rec < MAX_REC_PER_SLICE; rec++) {
            char* rec_ptr = raw_buf + rec * recl;
            unsigned short* vis_ptr = (unsigned short*)(rec_ptr +
                r->user.vispar.visinfo[bl].off);

            for (int ch = 0; ch < channels; ch++) {
                if (ch >= start_chan[rec] && ch <= end_chan[rec]) continue;

                float re = half_to_float(vis_ptr[ch * 2]);
                float im = half_to_float(vis_ptr[ch * 2 + 1]);
                float amp = sqrtf(re * re + im * im);

                if (isfinite(amp) && amp > 0) {
                    sum[ch] += amp;
                    count[ch]++;
                }
            }
        }

        for (int ch = 0; ch < channels; ch++) {
            r->bandpass[bl][ch] = (count[ch] > 0) ? sum[ch] / count[ch] : 1.0f;
        }

        free(sum);
        free(count);
    }
}

static void compute_off_source(struct SvfitsReaderImpl* r, int file_idx, char* raw_buf) {
    int channels = r->user.corr->daspar.channels;
    int baselines = r->user.baselines;
    size_t recl = r->user.corr->daspar.baselines * channels * sizeof(float);

    short* start_chan = r->user.bpass.start_chan;
    short* end_chan = r->user.bpass.end_chan;

    for (int bl = 0; bl < baselines; bl++) {
        float* sum_re = calloc(channels, sizeof(float));
        float* sum_im = calloc(channels, sizeof(float));
        int* count = calloc(channels, sizeof(int));

        for (int rec = 0; rec < MAX_REC_PER_SLICE; rec++) {
            char* rec_ptr = raw_buf + rec * recl;
            unsigned short* vis_ptr = (unsigned short*)(rec_ptr +
                r->user.vispar.visinfo[bl].off);

            for (int ch = 0; ch < channels; ch++) {
                if (ch >= start_chan[rec] && ch <= end_chan[rec]) continue;

                float re = half_to_float(vis_ptr[ch * 2]);
                float im = half_to_float(vis_ptr[ch * 2 + 1]);

                if (isfinite(re) && isfinite(im)) {
                    sum_re[ch] += re;
                    sum_im[ch] += im;
                    count[ch]++;
                }
            }
        }

        for (int ch = 0; ch < channels; ch++) {
            if (count[ch] > 0) {
                r->off_source_re[bl][ch] = sum_re[ch] / count[ch];
                r->off_source_im[bl][ch] = sum_im[ch] / count[ch];
            } else {
                r->off_source_re[bl][ch] = 0;
                r->off_source_im[bl][ch] = 0;
            }
        }

        free(sum_re);
        free(sum_im);
        free(count);
    }
}

static void compute_uvw(struct SvfitsReaderImpl* r, int bl, double mjd,
                        double* u, double* v, double* w) {
    int ant0 = r->user.vispar.visinfo[bl].ant0;
    int ant1 = r->user.vispar.visinfo[bl].ant1;

    double bx = r->user.corr->antenna[ant1].bx - r->user.corr->antenna[ant0].bx;
    double by = r->user.corr->antenna[ant1].by - r->user.corr->antenna[ant0].by;
    double bz = r->user.corr->antenna[ant1].bz - r->user.corr->antenna[ant0].bz;

    double ha = get_ha(&r->user, mjd);
    double ra = r->user.burst.ra_app;
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

    // Don't set srec->corr here - init_user does it

    printf("reader_create: allocations done\n"); fflush(stdout);
    return r;
}

int reader_init(SvfitsReader reader) {
    struct SvfitsReaderImpl* r = reader;
    printf("reader_init: starting\n"); fflush(stdout);

    char param_file[1024], antsamp[1024], bulletin[1024];
    strcpy(param_file, r->config.param_file);
    strcpy(antsamp, r->config.antsamp_file);

    char* bulletin_ptr = NULL;
    if (strlen(r->config.bulletin_a) > 0) {
        strcpy(bulletin, r->config.bulletin_a);
        bulletin_ptr = bulletin;
    }

    printf("reader_init: calling init_user with param=%s ant=%s\n", param_file, antsamp); fflush(stdout);

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

    // Allocate buffers
    int baselines = r->user.corr->daspar.baselines;
    int channels = r->user.corr->daspar.channels;
    size_t recl = baselines * channels * sizeof(float);
    r->raw_buffer_size = MAX_REC_PER_SLICE * recl;
    r->raw_buffer = malloc(r->raw_buffer_size);

    // Allocate bandpass arrays
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
    fprintf(stderr, "DEBUG reader_process: enter\n"); fflush(stderr);
    if (!r->initialized) {
        fprintf(stderr, "Reader not initialized\n");
        return 0;
    }

    size_t count = 0;
    int channels = r->user.corr->daspar.channels;
    int baselines = r->user.baselines;
    int nfiles = r->user.recfile.nfiles;
    size_t recl = r->user.corr->daspar.baselines * channels * sizeof(float);
    fprintf(stderr, "DEBUG: channels=%d baselines=%d nfiles=%d recl=%zu\n", channels, baselines, nfiles, recl); fflush(stderr);

    int file_order[MaxRecFiles];
    fprintf(stderr, "DEBUG: calling get_file_order\n"); fflush(stderr);
    get_file_order(&r->user, file_order);
    fprintf(stderr, "DEBUG: get_file_order done\n"); fflush(stderr);

    // Process each file
    for (int f = 0; f < nfiles; f++) {
        int file_idx = file_order[f];
        fprintf(stderr, "DEBUG: file %d/%d, file_idx=%d, n_rec=%d\n", f, nfiles, file_idx, r->user.recfile.n_rec[file_idx]); fflush(stderr);

        if (r->user.recfile.n_rec[file_idx] <= 0) continue;

        // Read raw data
        fprintf(stderr, "DEBUG: calling read_slice file_idx=%d\n", file_idx); fflush(stderr);
        if (read_slice(&r->user, file_idx, 0, r->raw_buffer) != 0) {
            fprintf(stderr, "Failed to read file %d\n", file_idx);
            continue;
        }
        fprintf(stderr, "DEBUG: read_slice done\n"); fflush(stderr);

        // Compute bandpass and off-source
        fprintf(stderr, "DEBUG: compute bandpass/off-source\n"); fflush(stderr);
        if (r->config.do_bandpass) {
            compute_bandpass(r, file_idx, r->raw_buffer);
        }
        if (r->config.do_baseline) {
            compute_off_source(r, file_idx, r->raw_buffer);
        }
        fprintf(stderr, "DEBUG: bandpass/off-source done\n"); fflush(stderr);

        // Process records with burst
        int start_rec = r->user.recfile.start_rec[file_idx];
        int n_rec = r->user.recfile.n_rec[file_idx];
        fprintf(stderr, "DEBUG: start_rec=%d n_rec=%d\n", start_rec, n_rec); fflush(stderr);

        for (int rec = start_rec; rec < start_rec + n_rec; rec++) {
            char* rec_ptr = r->raw_buffer + rec * recl;

            double t_rec = r->user.recfile.t_start[file_idx] +
                rec * r->user.recfile.t_slice / MAX_REC_PER_SLICE;
            double mjd = r->user.recfile.mjd_ref + t_rec / 86400.0;

            init_mat(&r->user, mjd);

            int start_ch = r->user.bpass.start_chan[rec];
            int end_ch = r->user.bpass.end_chan[rec];
            if (rec == start_rec) {
                fprintf(stderr, "DEBUG: rec=%d start_ch=%d end_ch=%d\n", rec, start_ch, end_ch); fflush(stderr);
            }

            // Process each baseline
            for (int bl = 0; bl < baselines; bl++) {
                unsigned short* vis_ptr = (unsigned short*)(rec_ptr +
                    r->user.vispar.visinfo[bl].off);

                double u, v, w;
                compute_uvw(r, bl, mjd, &u, &v, &w);

                // Collect for flagging
                float* amps = NULL;
                int n_amps = 0;
                if (r->config.do_flag) {
                    amps = malloc((end_ch - start_ch + 1) * sizeof(float));
                }

                // Process channels
                for (int ch = start_ch; ch <= end_ch; ch++) {
                    float re = half_to_float(vis_ptr[ch * 2]);
                    float im = half_to_float(vis_ptr[ch * 2 + 1]);

                    if (!isfinite(re) || !isfinite(im)) continue;

                    // Flip imaginary if needed
                    if (r->user.vispar.visinfo[bl].flip) {
                        im = -im;
                    }

                    // Bandpass correction
                    if (r->config.do_bandpass && r->bandpass[bl][ch] > 0) {
                        re /= r->bandpass[bl][ch];
                        im /= r->bandpass[bl][ch];
                    }

                    // Baseline subtraction
                    if (r->config.do_baseline) {
                        re -= r->off_source_re[bl][ch];
                        im -= r->off_source_im[bl][ch];
                    }

                    float weight = 1.0f;

                    // Collect amplitude for flagging
                    if (r->config.do_flag) {
                        amps[n_amps++] = sqrtf(re * re + im * im);
                    }

                    // Create visibility
                    CudaVisibility vis;
                    vis.re = re;
                    vis.im = im;
                    vis.weight = weight;
                    vis.u = (float)u;
                    vis.v = (float)v;
                    vis.w = (float)w;
                    vis.channel = ch;

                    // Apply flagging
                    if (r->config.do_flag && n_amps > 0) {
                        float med, mad;
                        if (robust_stats(n_amps, amps, &med, &mad) == 0) {
                            float amp = sqrtf(re * re + im * im);
                            if (fabsf(amp - med) > r->config.flag_threshold * mad) {
                                vis.weight = -1.0f;
                            }
                        }
                    }

                    // Callback
                    if (callback(&vis, user_data) != 0) {
                        if (amps) free(amps);
                        return count;
                    }
                    count++;
                }

                if (amps) free(amps);
            }
        }
    }

    return count;
}

void reader_free(SvfitsReader reader) {
    struct SvfitsReaderImpl* r = reader;
    if (!r) return;

    if (r->raw_buffer) free(r->raw_buffer);

    int baselines = r->user.corr->daspar.baselines;
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

    // Free svfits structures allocated with malloc
    if (r->user.srec) {
        if (r->user.srec->scan) free(r->user.srec->scan);
        free(r->user.srec);
    }
    if (r->user.corr) free(r->user.corr);
    if (r->user.hdr) free(r->user.hdr);

    free(r);
}
