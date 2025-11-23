/**
 * @file stub_prenut.c
 * @brief Stub for precession/nutation when no external library available
 *
 * This provides minimal no-op implementations for testing.
 */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "svio.h"
#include "gmrt_newcorr.h"

/**
 * Initialize precession matrices (identity for no-op)
 */
void init_mat(SvSelectionType *user, double tm) {
    /* Identity matrices */
    int i, j;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            user->i2eapp[i][j] = (i == j) ? 1.0 : 0.0;
            user->emean2i[i][j] = (i == j) ? 1.0 : 0.0;
            user->rmat[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }
}

/**
 * Apply precession to visibilities (no-op for stub)
 * Matches signature in svio.h
 */
void sla_prenut_vis(UvwParType *uvw, double mjd, double ra_app, double dec_app,
                    double epoch1) {
    /* No-op - coordinates unchanged */
    (void)uvw;
    (void)mjd;
    (void)ra_app;
    (void)dec_app;
    (void)epoch1;
}
