/**
 * @file test_init.c
 * @brief Minimal test for svfits init_user debugging
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "svio.h"
#include "gmrt_newcorr.h"

extern int init_user(SvSelectionType *user, char *fname, char *anthdr,
                     char *bhdrfile, char *bulletinA);

int main(int argc, char *argv[]) {
    char *param_file = "svfits_par.txt";
    char *antsamp_file = "antsamp.hdr";

    if (argc > 1) param_file = argv[1];
    if (argc > 2) antsamp_file = argv[2];

    printf("Test init_user\n");
    printf("  param_file: %s\n", param_file);
    printf("  antsamp_file: %s\n", antsamp_file);
    printf("\n");

    // Allocate exactly like original svfits.c does (lines 1128-1135)
    SvSelectionType user;
    memset(&user, 0, sizeof(user));

    printf("Allocating hdr...\n"); fflush(stdout);
    user.hdr = (InitHdrType*)malloc(sizeof(InitHdrType));
    if (!user.hdr) { fprintf(stderr, "Failed to alloc hdr\n"); return 1; }
    user.hdr->scans = 1;
    printf("  hdr=%p\n", (void*)user.hdr); fflush(stdout);

    printf("Allocating srec...\n"); fflush(stdout);
    user.srec = (ScanRecType*)malloc(sizeof(ScanRecType));
    if (!user.srec) { fprintf(stderr, "Failed to alloc srec\n"); return 1; }
    printf("  srec=%p\n", (void*)user.srec); fflush(stdout);

    printf("Allocating srec->scan...\n"); fflush(stdout);
    user.srec->scan = (ScanInfoType*)malloc(sizeof(ScanInfoType));
    if (!user.srec->scan) { fprintf(stderr, "Failed to alloc scan\n"); return 1; }
    printf("  srec->scan=%p\n", (void*)user.srec->scan); fflush(stdout);

    // bzero only source, like original
    bzero(&user.srec->scan->source, sizeof(SourceParType));

    printf("Allocating corr...\n"); fflush(stdout);
    user.corr = (CorrType*)malloc(sizeof(CorrType));
    if (!user.corr) { fprintf(stderr, "Failed to alloc corr\n"); return 1; }
    printf("  corr=%p\n", (void*)user.corr); fflush(stdout);

    printf("\nCalling init_user...\n"); fflush(stdout);

    int ret = init_user(&user, param_file, antsamp_file, NULL, NULL);

    printf("\ninit_user returned: %d\n", ret);

    // Cleanup
    if (user.srec && user.srec->scan) free(user.srec->scan);
    if (user.srec) free(user.srec);
    if (user.corr) free(user.corr);
    if (user.hdr) free(user.hdr);

    return ret;
}
