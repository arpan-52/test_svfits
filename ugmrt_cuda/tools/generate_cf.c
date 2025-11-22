/**
 * @file generate_cf.c
 * @brief Generate convolution function (gridding kernel) file
 *
 * Creates a binary CF file in the format expected by ugmrt_cuda:
 *   int32: support
 *   int32: oversampling
 *   float[]: values [(2*support+1)*oversamp]^2
 *
 * Uses prolate spheroidal wave function (PSWF) approximation.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/**
 * @brief Spheroidal function approximation (from AIPS/CASA)
 */
double spheroidal(double nu) {
    static const double P[] = {
        0.08203343, -0.3644705, 0.627866, -0.5335581, 0.2312756,
        0.004028559, -0.03697768, 0.1021332, -0.1201436, 0.06412774
    };
    static const double Q[] = {
        1.0, 0.8212018, 0.2078043,
        1.0, 0.9599102, 0.2918724
    };

    if (fabs(nu) > 1.0) return 0.0;

    double nusq = nu * nu;
    int idx = (nusq <= 0.75) ? 0 : 1;

    double top, bot;
    double delnusq = nusq - 0.75;

    if (idx == 0) {
        top = P[0] + delnusq * (P[1] + delnusq * (P[2] + delnusq * (P[3] + delnusq * P[4])));
        bot = Q[0] + delnusq * (Q[1] + delnusq * Q[2]);
    } else {
        top = P[5] + delnusq * (P[6] + delnusq * (P[7] + delnusq * (P[8] + delnusq * P[9])));
        bot = Q[3] + delnusq * (Q[4] + delnusq * Q[5]);
    }

    return (1.0 - nusq) * (top / bot);
}

void print_usage(const char* prog) {
    printf("Generate Convolution Function for uGMRT CUDA Imager\n\n");
    printf("Usage: %s -o output.cf [-s support] [-O oversampling]\n\n", prog);
    printf("Options:\n");
    printf("  -o FILE    Output CF file (required)\n");
    printf("  -s N       Support half-width in pixels (default: 7)\n");
    printf("  -O N       Oversampling factor (default: 128)\n");
    printf("  -h         Show this help\n");
}

int main(int argc, char* argv[]) {
    char output_file[1024] = "";
    int support = 7;
    int oversampling = 128;

    // Parse arguments
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
            strncpy(output_file, argv[++i], sizeof(output_file) - 1);
        } else if (strcmp(argv[i], "-s") == 0 && i + 1 < argc) {
            support = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-O") == 0 && i + 1 < argc) {
            oversampling = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-h") == 0) {
            print_usage(argv[0]);
            return 0;
        }
    }

    if (strlen(output_file) == 0) {
        fprintf(stderr, "Error: -o option is required\n\n");
        print_usage(argv[0]);
        return 1;
    }

    int full_size = (2 * support + 1) * oversampling;
    size_t n_values = full_size * full_size;

    printf("Generating CF:\n");
    printf("  Support: %d\n", support);
    printf("  Oversampling: %d\n", oversampling);
    printf("  Full size: %d x %d = %zu values\n", full_size, full_size, n_values);

    // Allocate
    float* values = malloc(n_values * sizeof(float));
    if (!values) {
        fprintf(stderr, "Error: Cannot allocate memory\n");
        return 1;
    }

    // Generate PSWF values
    double sum = 0.0;
    for (int iy = 0; iy < full_size; iy++) {
        // Normalized y coordinate [-1, 1]
        double y = (iy - support * oversampling) / (double)(support * oversampling);

        for (int ix = 0; ix < full_size; ix++) {
            // Normalized x coordinate [-1, 1]
            double x = (ix - support * oversampling) / (double)(support * oversampling);

            // 2D separable PSWF
            double val = spheroidal(x) * spheroidal(y);
            values[iy * full_size + ix] = (float)val;
            sum += val;
        }
    }

    // Normalize so sum = 1
    printf("  Sum before normalization: %.6f\n", sum);
    for (size_t i = 0; i < n_values; i++) {
        values[i] /= sum;
    }

    // Write to file
    FILE* fp = fopen(output_file, "wb");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open output file: %s\n", output_file);
        free(values);
        return 1;
    }

    fwrite(&support, sizeof(int), 1, fp);
    fwrite(&oversampling, sizeof(int), 1, fp);
    fwrite(values, sizeof(float), n_values, fp);
    fclose(fp);

    printf("  Output: %s (%.2f MB)\n", output_file,
           (8 + n_values * sizeof(float)) / (1024.0 * 1024.0));

    free(values);
    printf("Done!\n");

    return 0;
}
