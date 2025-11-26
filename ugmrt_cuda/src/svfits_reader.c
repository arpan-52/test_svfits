/**
 * @file svfits_reader.c
 * @brief svfits visibility reader implementation
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "svfits_reader.h"
#include "cf_generator.h"  // For cf_select_w_plane
#include "visibility_queue.h"  // For parallel file processing

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
extern int novas_prenut_vis(SvSelectionType *user, UvwParType *uvw, double mjd);
extern float half_to_float(const unsigned short x);
extern double get_ha(SvSelectionType *user, double tm);
extern int robust_stats(int n, float *x, float *med, float *mad);
extern int make_bpass(SvSelectionType *user, BpassType *bpass, char *rbuf, int idx, int slice);
extern int get_chan_num(double trec, SvSelectionType *user, int *cs, int *ce);

//-----------------------------------------------------------------------------
// Helper functions
//-----------------------------------------------------------------------------

static void compute_uvw(struct SvfitsReaderImpl* r, int bl, double mjd,
                        double* u, double* v, double* w, double freq_hz) {
    int ant0 = r->user.vispar.visinfo[bl].ant0;
    int ant1 = r->user.vispar.visinfo[bl].ant1;

    double bx = r->user.corr->antenna[ant1].bx - r->user.corr->antenna[ant0].bx;
    double by = r->user.corr->antenna[ant1].by - r->user.corr->antenna[ant0].by;
    double bz = r->user.corr->antenna[ant1].bz - r->user.corr->antenna[ant0].bz;

    double ha = get_ha(&r->user, mjd);
    // NOTE: Use source->dec_app (pointing center), NOT burst.dec_app
    // This matches svfits svgetUvw() which uses source->dec_app for the rotation matrix
    double dec = r->user.srec->scan->source.dec_app;

    double sin_ha = sin(ha);
    double cos_ha = cos(ha);
    double sin_dec = sin(dec);
    double cos_dec = cos(dec);

    // UVW in meters
    double u_m = sin_ha * bx + cos_ha * by;
    double v_m = -sin_dec * cos_ha * bx + sin_dec * sin_ha * by + cos_dec * bz;
    double w_m = cos_dec * cos_ha * bx - cos_dec * sin_ha * by + sin_dec * bz;

    // Convert to wavelengths at this frequency
    double wavelength = 299792458.0 / freq_hz;
    double u_wl = u_m / wavelength;
    double v_wl = v_m / wavelength;
    double w_wl = w_m / wavelength;

    // Apply J2000 precession/nutation rotation
    // This matches svfits svsubs.c which calls novas_prenut_vis() when USE_NOVAS is ON
    UvwParType uvwpar;
    uvwpar.u = (float)u_wl;
    uvwpar.v = (float)v_wl;
    uvwpar.w = (float)w_wl;

    // Rotate UVW from apparent coordinates to J2000
    // Requires SvSelectionType* which has rotation matrices (initialized by init_mat)
    novas_prenut_vis(&r->user, &uvwpar, r->user.recfile.mjd_ref);

    *u = uvwpar.u;
    *v = uvwpar.v;
    *w = uvwpar.w;
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

    // Extract freq info - NOTE: ch_width is always positive, net_sign gives direction
    r->freq_info.freq_start_hz = r->user.srec->scan->source.freq[0];
    r->freq_info.channel_width_hz = r->user.srec->scan->source.ch_width;
    r->freq_info.n_channels = r->user.corr->daspar.channels;

    // Signed channel width (net_sign determines if freq increases or decreases with channel)
    int net_sign = r->user.srec->scan->source.net_sign[0];
    double signed_ch_width = r->freq_info.channel_width_hz * net_sign;
    r->freq_info.freq_end_hz = r->freq_info.freq_start_hz +
        r->freq_info.n_channels * signed_ch_width;

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

    // Allocate bpass arrays for svfits (abp and off_src)
    int baselines = r->user.baselines;
    BpassType* bpass = &r->user.bpass;
    for (int bl = 0; bl < baselines; bl++) {
        bpass->abp[bl] = calloc(channels, sizeof(float));
        bpass->off_src[bl] = calloc(channels, sizeof(Complex));
        if (!bpass->abp[bl] || !bpass->off_src[bl]) {
            fprintf(stderr, "Failed to allocate bpass arrays\n");
            return -1;
        }
        // Initialize to nominal values
        for (int ch = 0; ch < channels; ch++) {
            bpass->abp[bl][ch] = 1.0f;
            bpass->off_src[bl][ch].r = 0.0f;
            bpass->off_src[bl][ch].i = 0.0f;
        }
    }

    // Also allocate our local arrays for convenience
    r->bandpass = malloc(baselines * sizeof(float*));
    r->off_source_re = malloc(baselines * sizeof(float*));
    r->off_source_im = malloc(baselines * sizeof(float*));

    for (int bl = 0; bl < baselines; bl++) {
        r->bandpass[bl] = bpass->abp[bl];  // Point to same memory
        r->off_source_re[bl] = calloc(channels, sizeof(float));
        r->off_source_im[bl] = calloc(channels, sizeof(float));
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

//-----------------------------------------------------------------------------
// Parallel file processing support structures
//-----------------------------------------------------------------------------

typedef struct {
    struct SvfitsReaderImpl* reader;  // Shared reader config
    VisibilityQueue* queue;            // Shared queue
    int file_idx;                      // File to process
    size_t* thread_count;              // Count for this thread
    size_t* thread_flagged;            // Flagged count
} FileProcessorArgs;

// Worker thread function to process a single file
static void process_file_to_queue(FileProcessorArgs* args) {
    struct SvfitsReaderImpl* r = args->reader;
    int file_idx = args->file_idx;
    VisibilityQueue* queue = args->queue;

    int start_rec = r->user.recfile.start_rec[file_idx];
    int n_rec = r->user.recfile.n_rec[file_idx];

    if (n_rec <= 0) {
        return;
    }

    // Allocate thread-local buffer for reading
    int channels = r->user.corr->daspar.channels;
    int das_baselines = r->user.corr->daspar.baselines;
    size_t recl = das_baselines * channels * sizeof(float);
    size_t buffer_size = MAX_REC_PER_SLICE * recl;
    char* thread_buffer = malloc(buffer_size);
    if (!thread_buffer) {
        fprintf(stderr, "Thread: Failed to allocate buffer for file %d\n", file_idx);
        return;
    }

    // Thread-local bandpass (each thread needs its own)
    BpassType thread_bpass;
    memset(&thread_bpass, 0, sizeof(BpassType));
    int baselines = r->user.baselines;
    for (int bl = 0; bl < baselines; bl++) {
        thread_bpass.abp[bl] = calloc(channels, sizeof(float));
        thread_bpass.off_src[bl] = calloc(channels, sizeof(Complex));
        for (int ch = 0; ch < channels; ch++) {
            thread_bpass.abp[bl][ch] = 1.0f;
            thread_bpass.off_src[bl][ch].r = 0.0f;
            thread_bpass.off_src[bl][ch].i = 0.0f;
        }
    }

    int rec_per_slice = r->user.recfile.rec_per_slice;
    double integ = r->user.corr->daspar.lta * r->user.statime;

    // Calculate slices
    int start_slice = start_rec / rec_per_slice;
    int end_rec = start_rec + n_rec - 1;
    int end_slice = end_rec / rec_per_slice;

    size_t local_count = 0;
    size_t local_flagged = 0;

    // Process each slice
    for (int slice = start_slice; slice <= end_slice; slice++) {
        if (read_slice(&r->user, file_idx, slice, thread_buffer) != 0) {
            fprintf(stderr, "Thread: Failed to read file %d slice %d\n", file_idx, slice);
            continue;
        }

        int n_vis = make_bpass(&r->user, &thread_bpass, thread_buffer, file_idx, slice);
        if (n_vis < 0) continue;

        int slice_start_rec = slice * rec_per_slice;
        int slice_end_rec = slice_start_rec + rec_per_slice - 1;
        int proc_start = (start_rec > slice_start_rec) ? start_rec : slice_start_rec;
        int proc_end = (end_rec < slice_end_rec) ? end_rec : slice_end_rec;
        int local_start = proc_start - slice_start_rec;
        int local_end = proc_end - slice_start_rec;

        // Process records
        for (int local_rec = local_start; local_rec <= local_end; local_rec++) {
            char* rec_ptr = thread_buffer + local_rec * recl;

            double t_rec = r->user.recfile.t_start[file_idx] +
                slice * r->user.recfile.slice_interval +
                (local_rec + 0.5) * integ;
            double mjd = r->user.recfile.mjd_ref + t_rec / 86400.0;

            init_mat(&r->user, mjd);

            int start_ch = thread_bpass.start_chan[local_rec];
            int end_ch = thread_bpass.end_chan[local_rec];
            if (start_ch < 0) start_ch = 0;
            if (end_ch >= channels) end_ch = channels - 1;
            if (end_ch < start_ch) continue;

            double freq0 = r->freq_info.freq_start_hz;
            int net_sign = r->user.srec->scan->source.net_sign[0];
            double ch_width = r->freq_info.channel_width_hz * net_sign;

            // Process baselines
            for (int bl = 0; bl < baselines; bl++) {
                unsigned short* vis_ptr = (unsigned short*)(rec_ptr +
                    r->user.vispar.visinfo[bl].off);

                float* abp = thread_bpass.abp[bl];
                Complex* off_src = thread_bpass.off_src[bl];

                for (int ch = start_ch; ch <= end_ch; ch++) {
                    float re = half_to_float(vis_ptr[ch * 2]);
                    float im = half_to_float(vis_ptr[ch * 2 + 1]);

                    if (!isfinite(re) || !isfinite(im)) {
                        local_flagged++;
                        continue;
                    }

                    if (r->user.do_base) {
                        re -= off_src[ch].r;
                        im -= off_src[ch].i;
                    }

                    if (r->user.do_band && abp[ch] > 0.0f) {
                        re /= abp[ch];
                        im /= abp[ch];
                    }

                    double freq_ch = freq0 + ch * ch_width;
                    double u, v, w;
                    compute_uvw(r, bl, mjd, &u, &v, &w, freq_ch);

                    // Create visibility
                    CudaVisibility vis = {0};
                    vis.re = re;
                    vis.im = im;
                    vis.weight = 1.0f;
                    vis.u = (float)u;
                    vis.v = (float)v;
                    vis.w = (float)w;
                    vis.channel = ch;
                    vis.freq = (float)freq_ch;
                    vis.d_phase = 0.0f;

                    if (r->config.n_w_planes > 0) {
                        vis.cf_cube = cf_select_w_plane((float)w, r->config.max_w, r->config.n_w_planes);
                    } else {
                        vis.cf_cube = 0;
                    }

                    vis.cf_grp = 0;
                    vis.grid_cube = 0;
                    vis.phase_grad_u = 0.0f;
                    vis.phase_grad_v = 0.0f;

                    // Push to queue
                    vis_queue_push(queue, &vis);
                    local_count++;

                    // Hermitian conjugate
                    CudaVisibility vis_conj = {0};
                    vis_conj.re = re;
                    vis_conj.im = -im;
                    vis_conj.weight = 1.0f;
                    vis_conj.u = (float)(-u);
                    vis_conj.v = (float)(-v);
                    vis_conj.w = (float)(-w);
                    vis_conj.channel = ch;
                    vis_conj.freq = (float)freq_ch;
                    vis_conj.d_phase = 0.0f;

                    if (r->config.n_w_planes > 0) {
                        vis_conj.cf_cube = cf_select_w_plane((float)(-w), r->config.max_w, r->config.n_w_planes);
                    } else {
                        vis_conj.cf_cube = 0;
                    }

                    vis_conj.cf_grp = 0;
                    vis_conj.grid_cube = 0;
                    vis_conj.phase_grad_u = 0.0f;
                    vis_conj.phase_grad_v = 0.0f;

                    vis_queue_push(queue, &vis_conj);
                    local_count++;
                }
            }
        }
    }

    // Cleanup thread-local resources
    for (int bl = 0; bl < baselines; bl++) {
        free(thread_bpass.abp[bl]);
        free(thread_bpass.off_src[bl]);
    }
    free(thread_buffer);

    // Update counters
    *args->thread_count = local_count;
    *args->thread_flagged = local_flagged;

    printf("  Thread finished file %d: %zu visibilities, %zu flagged\n",
           file_idx, local_count, local_flagged);
}

size_t reader_process(SvfitsReader reader, VisibilityCallback callback, void* user_data) {
    struct SvfitsReaderImpl* r = reader;
    if (!r->initialized) {
        fprintf(stderr, "Reader not initialized\n");
        return 0;
    }

    int channels = r->user.corr->daspar.channels;
    int baselines = r->user.baselines;
    int nfiles = r->user.recfile.nfiles;
    int num_threads = r->config.num_threads;
    if (num_threads <= 0) num_threads = 1;

    int file_order[MaxRecFiles];
    get_file_order(&r->user, file_order);

    printf("PARALLEL Processing: %d files, %d baselines, %d channels, %d threads\n",
           nfiles, baselines, channels, num_threads);

    // Create shared visibility queue (10M visibilities = ~800MB buffer)
    VisibilityQueue queue;
    size_t queue_capacity = 10000000;  // Large buffer for maximum throughput
    if (vis_queue_init(&queue, queue_capacity) != 0) {
        fprintf(stderr, "Failed to create visibility queue\n");
        return 0;
    }

    // Set active producers
    queue.active_producers = nfiles;  // Each file is a producer

    // Thread statistics
    size_t* thread_counts = calloc(nfiles, sizeof(size_t));
    size_t* thread_flagged = calloc(nfiles, sizeof(size_t));

    // Prepare arguments for each file
    FileProcessorArgs* args = calloc(nfiles, sizeof(FileProcessorArgs));
    for (int f = 0; f < nfiles; f++) {
        args[f].reader = r;
        args[f].queue = &queue;
        args[f].file_idx = file_order[f];
        args[f].thread_count = &thread_counts[f];
        args[f].thread_flagged = &thread_flagged[f];
    }

    // Launch producer threads in parallel
    printf("Launching %d producer threads...\n", nfiles);
    #pragma omp parallel for num_threads(num_threads)
    for (int f = 0; f < nfiles; f++) {
        process_file_to_queue(&args[f]);
        vis_queue_producer_done(&queue);
    }

    // Consumer thread (main thread) - pop from queue and grid on GPU
    printf("Starting GPU consumer...\n");
    size_t count = 0;
    double u_min = 1e30, u_max = -1e30, v_min = 1e30, v_max = -1e30;

    CudaVisibility batch[10000];  // Batch size for popping
    while (!vis_queue_is_finished(&queue)) {
        size_t n = vis_queue_pop_batch(&queue, batch, 10000);

        if (n == 0) break;  // Queue finished

        // Process batch
        for (size_t i = 0; i < n; i++) {
            CudaVisibility* vis = &batch[i];

            // Track UV stats
            if (vis->u < u_min) u_min = vis->u;
            if (vis->u > u_max) u_max = vis->u;
            if (vis->v < v_min) v_min = vis->v;
            if (vis->v > v_max) v_max = vis->v;

            // Debug first few
            if (count < 5) {
                printf("  [vis %zu] u=%.1f v=%.1f w=%.1f re=%.3f im=%.3f cf_cube=%d\n",
                       count, vis->u, vis->v, vis->w, vis->re, vis->im, vis->cf_cube);
            }

            // Callback to gridder
            if (callback(vis, user_data) != 0) {
                vis_queue_destroy(&queue);
                free(thread_counts);
                free(thread_flagged);
                free(args);
                return count;
            }
            count++;
        }

        // Print progress every 1M visibilities
        if (count % 1000000 == 0) {
            printf("  Gridded %zu M visibilities...\n", count / 1000000);
        }
    }

    // Sum up thread statistics
    size_t total_flagged = 0;
    for (int f = 0; f < nfiles; f++) {
        total_flagged += thread_flagged[f];
    }

    printf("  Total: %zu visibilities\n", count);
    printf("  Flagged: %zu (%.1f%%)\n", total_flagged, 100.0 * total_flagged / (count + total_flagged));
    printf("  UV range (wavelengths): u=[%.1f, %.1f] v=[%.1f, %.1f]\n",
           u_min, u_max, v_min, v_max);

    // Cleanup
    vis_queue_destroy(&queue);
    free(thread_counts);
    free(thread_flagged);
    free(args);

    return count;
}

void reader_free(SvfitsReader reader) {
    struct SvfitsReaderImpl* r = reader;
    if (!r) return;

    if (r->raw_buffer) free(r->raw_buffer);

    if (r->initialized) {
        int baselines = r->user.baselines;
        BpassType* bpass = &r->user.bpass;

        // Free bpass arrays (abp and off_src)
        for (int bl = 0; bl < baselines; bl++) {
            if (bpass->abp[bl]) free(bpass->abp[bl]);
            if (bpass->off_src[bl]) free(bpass->off_src[bl]);
        }

        // Free our convenience pointers (bandpass points to abp, don't double-free)
        if (r->bandpass) free(r->bandpass);
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

//-----------------------------------------------------------------------------
// Batch/Parallel Processing Implementation
//-----------------------------------------------------------------------------

VisibilityBuffer* visbuf_create(size_t capacity) {
    VisibilityBuffer* buf = malloc(sizeof(VisibilityBuffer));
    if (!buf) return NULL;

    buf->vis = malloc(capacity * sizeof(CudaVisibility));
    if (!buf->vis) {
        free(buf);
        return NULL;
    }
    buf->count = 0;
    buf->capacity = capacity;
    return buf;
}

void visbuf_free(VisibilityBuffer* buf) {
    if (!buf) return;
    if (buf->vis) free(buf->vis);
    free(buf);
}

void visbuf_reset(VisibilityBuffer* buf) {
    if (buf) buf->count = 0;
}

size_t reader_estimate_vis_count(SvfitsReader reader) {
    struct SvfitsReaderImpl* r = reader;
    if (!r->initialized) return 0;

    // Estimate based on: nfiles * n_rec * baselines * channels_per_burst * 2 (hermitian)
    int nfiles = r->user.recfile.nfiles;
    int baselines = r->user.baselines;
    int channels = r->user.corr->daspar.channels;

    // Sum up n_rec across all files
    size_t total_recs = 0;
    for (int f = 0; f < nfiles; f++) {
        total_recs += r->user.recfile.n_rec[f];
    }

    // Assume ~10% of channels contain burst (conservative estimate)
    size_t est_channels = channels / 10;
    if (est_channels < 100) est_channels = 100;

    // Each visibility generates 2 entries (original + hermitian conjugate)
    size_t estimate = total_recs * baselines * est_channels * 2;

    printf("Estimated visibility count: %zu (recs=%zu, baselines=%d, est_channels=%zu)\n",
           estimate, total_recs, baselines, est_channels);

    return estimate;
}

// Internal helper to add visibility to buffer (with optional resizing)
static int visbuf_add(VisibilityBuffer* buf, const CudaVisibility* vis) {
    if (buf->count >= buf->capacity) {
        // Grow buffer by 50%
        size_t new_capacity = buf->capacity + buf->capacity / 2;
        CudaVisibility* new_vis = realloc(buf->vis, new_capacity * sizeof(CudaVisibility));
        if (!new_vis) {
            fprintf(stderr, "Failed to grow visibility buffer\n");
            return -1;
        }
        buf->vis = new_vis;
        buf->capacity = new_capacity;
    }
    buf->vis[buf->count++] = *vis;
    return 0;
}

ssize_t reader_extract_all(SvfitsReader reader, VisibilityBuffer* buf) {
    struct SvfitsReaderImpl* r = reader;
    if (!r->initialized) {
        fprintf(stderr, "Reader not initialized\n");
        return -1;
    }

    visbuf_reset(buf);

    size_t flagged = 0;
    int channels = r->user.corr->daspar.channels;
    int baselines = r->user.baselines;
    int nfiles = r->user.recfile.nfiles;
    int rec_per_slice = r->user.recfile.rec_per_slice;
    size_t recl = r->user.corr->daspar.baselines * channels * sizeof(float);
    double integ = r->user.corr->daspar.lta * r->user.statime;
    BpassType* bpass = &r->user.bpass;

    int file_order[MaxRecFiles];
    get_file_order(&r->user, file_order);

    printf("Extracting visibilities: %d files, %d baselines, %d channels\n",
           nfiles, baselines, channels);

    // Process each file in burst order
    for (int f = 0; f < nfiles; f++) {
        int file_idx = file_order[f];
        int start_rec = r->user.recfile.start_rec[file_idx];
        int n_rec = r->user.recfile.n_rec[file_idx];

        if (n_rec <= 0) continue;

        int start_slice = start_rec / rec_per_slice;
        int end_rec = start_rec + n_rec - 1;
        int end_slice = end_rec / rec_per_slice;

        // Process each slice
        for (int slice = start_slice; slice <= end_slice; slice++) {
            if (read_slice(&r->user, file_idx, slice, r->raw_buffer) != 0) {
                fprintf(stderr, "Failed to read file %d slice %d\n", file_idx, slice);
                continue;
            }

            int n_vis = make_bpass(&r->user, bpass, r->raw_buffer, file_idx, slice);
            if (n_vis < 0) continue;

            int slice_start_rec = slice * rec_per_slice;
            int slice_end_rec = slice_start_rec + rec_per_slice - 1;
            int proc_start = (start_rec > slice_start_rec) ? start_rec : slice_start_rec;
            int proc_end = (end_rec < slice_end_rec) ? end_rec : slice_end_rec;
            int local_start = proc_start - slice_start_rec;
            int local_end = proc_end - slice_start_rec;

            // Process records
            for (int local_rec = local_start; local_rec <= local_end; local_rec++) {
                char* rec_ptr = r->raw_buffer + local_rec * recl;
                double t_rec = r->user.recfile.t_start[file_idx] +
                    slice * r->user.recfile.slice_interval +
                    (local_rec + 0.5) * integ;
                double mjd = r->user.recfile.mjd_ref + t_rec / 86400.0;

                init_mat(&r->user, mjd);

                int start_ch = bpass->start_chan[local_rec];
                int end_ch = bpass->end_chan[local_rec];
                if (start_ch < 0) start_ch = 0;
                if (end_ch >= channels) end_ch = channels - 1;
                if (end_ch < start_ch) continue;

                double freq0 = r->freq_info.freq_start_hz;
                int net_sign = r->user.srec->scan->source.net_sign[0];
                double ch_width = r->freq_info.channel_width_hz * net_sign;

                // Process baselines (can be parallelized with OpenMP)
                #pragma omp parallel for num_threads(r->config.num_threads) schedule(dynamic)
                for (int bl = 0; bl < baselines; bl++) {
                    unsigned short* vis_ptr = (unsigned short*)(rec_ptr +
                        r->user.vispar.visinfo[bl].off);

                    float* abp = bpass->abp[bl];
                    Complex* off_src = bpass->off_src[bl];

                    for (int ch = start_ch; ch <= end_ch; ch++) {
                        float re = half_to_float(vis_ptr[ch * 2]);
                        float im = half_to_float(vis_ptr[ch * 2 + 1]);

                        if (!isfinite(re) || !isfinite(im)) {
                            #pragma omp atomic
                            flagged++;
                            continue;
                        }

                        if (r->user.do_base) {
                            re -= off_src[ch].r;
                            im -= off_src[ch].i;
                        }

                        if (r->user.do_band && abp[ch] > 0.0f) {
                            re /= abp[ch];
                            im /= abp[ch];
                        }

                        double freq_ch = freq0 + ch * ch_width;
                        double u, v, w;
                        compute_uvw(r, bl, mjd, &u, &v, &w, freq_ch);

                        // Create visibility
                        CudaVisibility vis = {0};
                        vis.re = re;
                        vis.im = im;
                        vis.weight = 1.0f;
                        vis.u = (float)u;
                        vis.v = (float)v;
                        vis.w = (float)w;
                        vis.channel = ch;
                        vis.freq = (float)freq_ch;

                        // Add to buffer (thread-safe with critical section)
                        #pragma omp critical
                        {
                            visbuf_add(buf, &vis);

                            // Add Hermitian conjugate
                            CudaVisibility vis_conj = {0};
                            vis_conj.re = re;
                            vis_conj.im = -im;
                            vis_conj.weight = 1.0f;
                            vis_conj.u = (float)(-u);
                            vis_conj.v = (float)(-v);
                            vis_conj.w = (float)(-w);
                            vis_conj.channel = ch;
                            vis_conj.freq = (float)freq_ch;
                            visbuf_add(buf, &vis_conj);
                        }
                    }
                }
            }
        }
    }

    printf("  Extracted: %zu visibilities\n", buf->count);
    printf("  Flagged: %zu\n", flagged);

    return (ssize_t)buf->count;
}
